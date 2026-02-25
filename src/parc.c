#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>
#include <float.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#define INF INT_MAX

#define BUF_SIZE (1 << 20)  // 1 MB

// ---------- Data Structure ----------
typedef struct {
	uint32_t to;
	double distance;
	double time;
} Edge;

typedef struct {
	uint32_t from;
	double distance;
	double time;
} InEdge;

typedef struct {
	Edge *out_edges;
	InEdge *in_edges;
	uint32_t out_count, in_count;
	uint32_t out_cap, in_cap;
	bool active;
	double x;
	double y;
	double dist_from_source;
	double dist_to_dest;
	double time_from_source;
	double time_to_dest;
} Node;

typedef struct {
	Node *nodes;
	uint32_t n_nodes;
	uint32_t n_edges;
} Graph;

// ---------- Heap ----------
typedef struct {
	uint32_t node;
	double cost;
} HeapNode;

typedef struct {
	HeapNode *data;
	uint32_t size;
	uint32_t capacity;
} MinHeap;

typedef struct {
	double cost;
	double resource;
	double bound_path;
	uint32_t *path;
	uint32_t path_len;
} AStarResult;

typedef struct {
	double cost;
	double resource;
	double timeFound;
	double lambda;
	int iter;
} Improvement;

typedef struct {
	Improvement *items;
	size_t count;
	size_t capacity;
} ImprovementList;

static void improvements_init(ImprovementList *list) {
	list->items = NULL;
	list->count = 0;
	list->capacity = 0;
}

static void improvements_free(ImprovementList *list) {
	free(list->items);
	list->items = NULL;
	list->count = 0;
	list->capacity = 0;
}

static void improvements_push(ImprovementList *list, double cost, double resource,
                              double timeFound, double lambda, int iter) {
	if (list->count == list->capacity) {
		size_t new_cap = list->capacity ? list->capacity * 2 : 8;
		Improvement *tmp = realloc(list->items, new_cap * sizeof(Improvement));
		if (!tmp) exit(EXIT_FAILURE);
		list->items = tmp;
		list->capacity = new_cap;
	}
	list->items[list->count++] = (Improvement){cost, resource, timeFound, lambda, iter};
}

static void improvements_write_file(const char *filename, const ImprovementList *list) {
	FILE *f = fopen(filename, "w");
	if (!f) {
		perror("Error opening improvements file");
		return;
	}
	fprintf(f, "cost resource timeFound lambda iter\n");
	for (size_t i = 0; i < list->count; i++) {
		const Improvement *imp = &list->items[i];
		fprintf(f, "%.0f %.0f %.6f %.6f %d\n",
		        imp->cost, imp->resource, imp->timeFound, imp->lambda, imp->iter);
	}
	fclose(f);
}

static double elapsed_seconds(clock_t start) {
	return ((double) (clock() - start)) / CLOCKS_PER_SEC;
}

static int ensure_output_dir(const char *dir) {
	if (!dir || dir[0] == '\0' || strcmp(dir, ".") == 0) return 0;
	struct stat st;
	if (stat(dir, &st) == 0) {
		if (S_ISDIR(st.st_mode)) return 0;
		fprintf(stderr, "Output path exists but is not a directory: %s\n", dir);
		return -1;
	}
	if (mkdir(dir, 0755) != 0) {
		perror("mkdir");
		return -1;
	}
	return 0;
}

static void build_output_path(char *dst, size_t dst_size, const char *dir, const char *filename) {
	if (!dir || dir[0] == '\0' || strcmp(dir, ".") == 0) {
		snprintf(dst, dst_size, "%s", filename);
		return;
	}
	size_t len = strlen(dir);
	if (dir[len - 1] == '/')
		snprintf(dst, dst_size, "%s%s", dir, filename);
	else
		snprintf(dst, dst_size, "%s/%s", dir, filename);
}

// ---------- Auxiliary Functions ----------
void free_astar_result(AStarResult *r) {
	if (r->path) {
		free(r->path);
		r->path = NULL;
	}
}

#define Pi 0.000000017453292519943295 // Pi / 180 / 1,000,000
#define COEF 61161607.40544 // = Earth radius M * 0.96 * 10
//NOTE: The input coordinates are in nanodegrees (to transform them in degrees they must be multiplied by 1,000,000)
// Length data are given in decimeters (0.1m)

// Evaluate the distance between two points from their coordinates
double haversine(int32_t lon1, int32_t lat1, int32_t lon2, int32_t lat2) {
	double Delta_lat_squared = pow(fabs(lat1 - lat2) * Pi, 2);
	double cos_mean_lat = cos((lat1 + lat2) * Pi / 2);
	double x = cos_mean_lat * fabs(lon1 - lon2) * Pi;
	return floor(COEF * sqrt(Delta_lat_squared + pow(x, 2)));
}

// ---------- Heap ----------
MinHeap *createHeap(uint32_t capacity) {
	MinHeap *h = malloc(sizeof(MinHeap));
	if (!h) exit(EXIT_FAILURE);
	h->data = malloc(sizeof(HeapNode) * (capacity ? capacity : 4));
	if (!h->data) exit(EXIT_FAILURE);
	h->size = 0;
	h->capacity = capacity ? capacity : 4;
	return h;
}

void swap(HeapNode *a, HeapNode *b) {
	HeapNode t = *a;
	*a = *b;
	*b = t;
}

void heapifyUp(MinHeap *h, uint32_t i) {
	while (i && h->data[i].cost < h->data[(i - 1) / 2].cost) {
		swap(&h->data[i], &h->data[(i - 1) / 2]);
		i = (i - 1) / 2;
	}
}

void heapifyDown(MinHeap *h, uint32_t i) {
	uint32_t smallest = i;
	uint32_t l = 2 * i + 1, r = 2 * i + 2;
	if (l < h->size && h->data[l].cost < h->data[smallest].cost) smallest = l;
	if (r < h->size && h->data[r].cost < h->data[smallest].cost) smallest = r;
	if (smallest != i) {
		swap(&h->data[i], &h->data[smallest]);
		heapifyDown(h, smallest);
	}
}

void push(MinHeap *h, uint32_t node, double cost) {
	if (h->size == h->capacity) {
		h->capacity = h->capacity ? h->capacity * 2 : 4;
		HeapNode *tmp = realloc(h->data, sizeof(HeapNode) * h->capacity);
		if (!tmp) exit(EXIT_FAILURE);
		h->data = tmp;
	}
	h->data[h->size++] = (HeapNode){node, cost};
	heapifyUp(h, h->size - 1);
}

HeapNode pop(MinHeap *h) {
	if (h->size == 0) return (HeapNode){UINT32_MAX, DBL_MAX};
	HeapNode root = h->data[0];
	h->data[0] = h->data[--h->size];
	heapifyDown(h, 0);
	return root;
}

bool isEmpty(MinHeap *h) { return h->size == 0; }

// ---------- Graph ----------
Graph *createGraph(uint32_t n_nodes, uint32_t n_edges) {
	Graph *g = malloc(sizeof(Graph));
	g->n_nodes = n_nodes;
	g->n_edges = n_edges;
	g->nodes = calloc(n_nodes, sizeof(Node));
	return g;
}

void addEdge(Graph *g, uint32_t from, uint32_t to, double distance, double time) {
	Node *n = &g->nodes[from];
	if (n->out_count == n->out_cap) {
		n->out_cap = n->out_cap ? n->out_cap * 2 : 4;
		n->out_edges = realloc(n->out_edges, n->out_cap * sizeof(Edge));
	}
	n->out_edges[n->out_count++] = (Edge){to, distance, time};

	Node *m = &g->nodes[to];
	if (m->in_count == m->in_cap) {
		m->in_cap = m->in_cap ? m->in_cap * 2 : 4;
		m->in_edges = realloc(m->in_edges, m->in_cap * sizeof(InEdge));
	}
	m->in_edges[m->in_count++] = (InEdge){from, distance, time};
}

AStarResult copy_res(const AStarResult *src) {
	AStarResult dst;
	dst.path_len = src->path_len;
	dst.cost = src->cost;
	dst.resource = src->resource;

	if (src->path != NULL && src->path_len > 0) {
		dst.path = malloc(src->path_len * sizeof(uint32_t));
		if (!dst.path) {
			perror("malloc failed in copy_res");
			exit(EXIT_FAILURE);
		}
		memcpy(dst.path, src->path, src->path_len * sizeof(uint32_t));
	} else {
		dst.path = NULL;
	}
	return dst;
}

void freeGraph(Graph *g) {
	if (!g) return;
	if (g->nodes) {
		for (uint32_t i = 0; i < g->n_nodes; i++) {
			free(g->nodes[i].out_edges);
			free(g->nodes[i].in_edges);
		}
		free(g->nodes);
	}
	free(g);
}

// ---------- A* ----------
AStarResult a_star_with_bound(Graph *g, uint32_t start,
                              uint32_t goal, double W, double dist_best, bool foward, uint32_t s, uint32_t d,
                              double budget, double best_feas_sol, double lambda, uint32_t destination) {
	uint32_t n = g->n_nodes;
	double *cost = malloc(sizeof(double) * n);
	uint32_t *parent = malloc(sizeof(uint32_t) * n);
	for (uint32_t i = 0; i < n; i++) {
		cost[i] = DBL_MAX;
		parent[i] = UINT32_MAX;
	}
	MinHeap *heap = createHeap(n);
	cost[start] = 0;
	push(heap, start, 0);

	while (!isEmpty(heap)) {
		HeapNode curr = pop(heap);
		uint32_t u = curr.node;
		if (!g->nodes[u].active) continue;
		if (u == goal) break; // stop as soon as you extract the goal node

		double current_cost = cost[u];
		double heuristic = 0;

		if (lambda == 0 && foward) g->nodes[u].time_from_source = cost[u];
		if (lambda == 0 && !foward) g->nodes[u].time_to_dest = cost[u];
		if (lambda == 1 && foward)
			g->nodes[u].dist_from_source = cost[u];
		if (lambda == 1 && !foward)
			g->nodes[u].dist_to_dest = cost[u];

		if ((g->nodes[u].dist_from_source + g->nodes[u].dist_to_dest >= dist_best) || (
			    g->nodes[u].time_from_source + g->nodes[u].time_to_dest > W)) {
			g->nodes[u].active = false;
			continue;
		}

		if (foward) heuristic = lambda * g->nodes[u].dist_to_dest + (1 - lambda) * g->nodes[u].time_to_dest;
		else heuristic = lambda * g->nodes[u].dist_from_source + (1 - lambda) * g->nodes[u].time_from_source;

		double estimated_total = current_cost + heuristic;

		if (lambda == 0 && estimated_total > W) {
			g->nodes[u].active = false;
			continue;
		}

		if (lambda == 1 && estimated_total >= dist_best) {
			g->nodes[u].active = false;
			continue;
		}

		if (lambda > 0 && lambda < 1 && estimated_total > (lambda * dist_best + (1 - lambda) * W))
			continue;

		for (uint32_t i = 0; i < (foward ? g->nodes[u].out_count : g->nodes[u].in_count); i++) {
			uint32_t v;
			double edge_cost;
			if (foward) {
				Edge e = g->nodes[u].out_edges[i];
				v = e.to;
				edge_cost = lambda * e.distance + (1 - lambda) * e.time;
			} else {
				InEdge e = g->nodes[u].in_edges[i];
				v = e.from;
				edge_cost = lambda * e.distance + (1 - lambda) * e.time;
			}

			if (!g->nodes[v].active) continue;

			double new_cost = cost[u] + edge_cost;
			if (new_cost < cost[v]) {
				cost[v] = new_cost;
				parent[v] = u;

				if (lambda == 0 && foward) /* && g->nodes[u].time_from_source == 0)*/
					g->nodes[v].time_from_source = cost[v];
				if (lambda == 0 && !foward)/* && g->nodes[u].time_to_dest == 0)*/ g->nodes[v].time_to_dest = cost[v];
				if (lambda == 1 && foward) /* && g->nodes[u].dist_from_source == g->nodes[u].spherical_from_source)*/
					g->nodes[v].dist_from_source = cost[v];
				if (lambda == 1 && !foward) /* && g->nodes[u].dist_to_dest == g->nodes[u].spherical_to_dest)*/
					g->nodes[v].dist_to_dest = cost[v];

				double h = 0;
				if (foward) h = lambda * g->nodes[v].dist_to_dest + (1 - lambda) * g->nodes[v].time_to_dest;
				else h = lambda * g->nodes[v].dist_from_source + (1 - lambda) * g->nodes[v].time_from_source;

				double f = new_cost + h;
				if (f > (lambda * dist_best + (1 - lambda) * W))
					continue;
				push(heap, v, f);
			}
		}
		if (u == d || u == s) {
			g->nodes[u].active = false;
			continue;
		}
	}

	g->nodes[s].active = true;
	g->nodes[d].active = true;

	if (goal == -1 && destination != -1)
		goal = destination;

	AStarResult res;
	uint32_t *path = malloc(sizeof(uint32_t) * n);
	uint32_t length = 0;
	uint32_t v = goal;
	res.cost = 0;
	res.resource = 0;
	res.bound_path = 0;
	if (goal != -1 && cost[goal] != DBL_MAX)
		while (v != start) {
			path[length++] = v;
			if (goal != -1 && v != start && v != goal)
				if (g->nodes[v].time_from_source + g->nodes[v].time_to_dest > res.bound_path)
					res.bound_path = g->nodes[v].time_from_source + g->nodes[v].time_to_dest;

			uint32_t u = parent[v];
			if (u == UINT32_MAX) break;

			double add_path_dist = DBL_MAX;
			double add_path_time = DBL_MAX;
			for (uint32_t i = 0; i < g->nodes[u].out_count; i++) {
				Edge e = g->nodes[u].out_edges[i];
				if (e.to == v) {
					if (add_path_dist > e.distance)
						add_path_dist = e.distance;
					if (add_path_time > e.time)
						add_path_time = e.time;
				}
			}
			res.cost += add_path_dist;
			res.resource += add_path_time;
			v = u;
		}

	if (goal != -1 && cost[goal] == DBL_MAX) {
		res.cost = DBL_MAX;
		res.resource = DBL_MAX;
	}
	res.path = path;
	res.path_len = length;

	free(cost);
	free(heap->data);
	free(heap);
	free(parent);
	return res;
}

// ---------- Parser ----------
void parse_args(int argc, char *argv[],
                char **input_path, uint32_t *s, uint32_t *d,
                uint64_t *W, double *time_limit, int *run_reduction_heuristic, int *max_iterations,
                int *multipleinstanceflag, char **inputinstance, char **output_dir, double *perc_red) {
	*input_path = NULL;
	*s = *d = 0;
	*W = 0;
	*time_limit = 0.0;
	*run_reduction_heuristic = 1;
	*max_iterations = 0;
	*multipleinstanceflag = 0;
	*inputinstance = NULL;
	*output_dir = (char *) ".";
	*perc_red = 0.0;
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "--input")) {
			*input_path = argv[++i];
		} else if (!strcmp(argv[i], "--s")) {
			*s = (uint32_t) atoi(argv[++i]);
		} else if (!strcmp(argv[i], "--d")) {
			*d = (uint32_t) atoi(argv[++i]);
		} else if (!strcmp(argv[i], "--W")) {
			*W = strtoull(argv[++i], NULL, 10);
		} else if (!strcmp(argv[i], "--tl")) {
			*time_limit = atof(argv[++i]);
		} else if (!strcmp(argv[i], "--redh")) {
			*run_reduction_heuristic = atoi(argv[++i]);
		} else if (!strcmp(argv[i], "--nit")) {
			*max_iterations = atoi(argv[++i]);
		} else if (!strcmp(argv[i], "--multipleinstanceflag")) {
			*multipleinstanceflag = atoi(argv[++i]);
		} else if (!strcmp(argv[i], "--inputinstance")) {
			*inputinstance = argv[++i];
		} else if (!strcmp(argv[i], "--outdir")) {
			*output_dir = argv[++i];
		} else if (!strcmp(argv[i], "--perc_red")) {
			*perc_red = atof(argv[++i]);
		} else {
			fprintf(stderr, "Unknown argument: %s\n", argv[i]);
			exit(EXIT_FAILURE);
		}
	}

	if (!*input_path) {
		fprintf(stderr,
		        "Usage: ./compute_paths "
		        "--input file --s src --d dest --W bound "
		        "[--tl time_limit_seconds] [--nit max_iterations] [--multipleinstanceflag run all instances] "
		        "[--inputinstance file run all instances] [--outdir output_directory] "
		        "[--perc_red percent_reduction]\n"
		);
		exit(EXIT_FAILURE);
	}
}

// ---------- Main ----------
int main(int argc, char *argv[]) {
	char *filename;
	char *inputinstance;
	char *output_dir;
	uint32_t s, d;
	double lambda = 0.5;
	uint64_t W_read;
	double time_limit;
	int run_reduction_heuristic = 1;
	int max_iterations;
	int multipleinstanceflag;
	double perc_red;
	double cost_spr, resource_spr, cost_spc, resource_spc;

	parse_args(argc, argv, &filename, &s, &d, &W_read, &time_limit, &run_reduction_heuristic,
	           &max_iterations, &multipleinstanceflag, &inputinstance, &output_dir, &perc_red);

	if (ensure_output_dir(output_dir) != 0) return 1;

	//-----------------------------------------START GRAPH READING
	//-----------------------------------------
	//-----------------------------------------
	//-----------------------------------------
	clock_t start_read = clock();
	clock_t end_read = clock();

	FILE *f = fopen(filename, "r");
	if (!f) {
		perror("Error opening file");
		return 1;
	}

	uint32_t n_nodes, n_edges;
	fscanf(f, "nodes %u edges %u\n", &n_nodes, &n_edges);
	fflush(f);
	fclose(f);

	Graph *g = createGraph(n_nodes, n_edges);

	FILE *fp = fopen(filename, "rb");
	if (!fp) {
		perror("fopen");
		return 1;
	}

	static unsigned char buf[BUF_SIZE];
	size_t n;
	int idx = 0;
	uint32_t id = 0;
	uint32_t node_i = 0;
	uint32_t edge_i = 0;
	uint32_t from;
	uint32_t to;
	uint32_t dist;
	uint32_t time;
	int32_t num = 0;
	int32_t vals[4];
	char header[256];
	int sign = 1;
	int reading = 0;
	int is_node = 1;

	fgets(header, sizeof(header), fp);
	sscanf(header, "nodes %u edges %u", &n_nodes, &n_edges);
	while ((n = fread(buf, 1, BUF_SIZE, fp)) > 0) {
		for (size_t i = 0; i < n; i++) {
			unsigned char c = buf[i];
			if (c == 'v' || c == 'e') {
				idx = 0;
				continue;
			}
			if (c == '-') {
				sign = -1;
				reading = 1;
				continue;
			}
			if (c >= '0' && c <= '9') {
				num = num * 10 + (c - '0');
				reading = 1;
				continue;
			}
			if (reading) {
				vals[idx++] = sign * num;
				num = 0;
				sign = 1;
				reading = 0;
				if (is_node && idx == 3) {
					id = (uint32_t) vals[0];
					g->nodes[id].x = vals[1];
					g->nodes[id].y = vals[2];
					g->nodes[id].active = true;
					g->nodes[id].time_from_source = 0;
					g->nodes[id].time_to_dest = 0;
					g->nodes[id].dist_from_source = 0;
					g->nodes[id].dist_to_dest = 0;
					node_i++;
					idx = 0;
					if (node_i == n_nodes) is_node = 0;
				} else if (!is_node && idx == 4) {
					from = (uint32_t) vals[0];
					to = (uint32_t) vals[1];
					dist = (uint32_t) vals[2];
					time = (uint32_t) vals[3];
					addEdge(g, from, to, (double) dist, (double) time);
					if (edge_i < n_nodes) {
						g->nodes[edge_i].dist_from_source = haversine(g->nodes[s].x, g->nodes[s].y, g->nodes[edge_i].x,
						                                              g->nodes[edge_i].y);
						g->nodes[edge_i].dist_to_dest = haversine(g->nodes[d].x, g->nodes[d].y, g->nodes[edge_i].x,
						                                          g->nodes[edge_i].y);
					}
					edge_i++;
					idx = 0;
				}
			}
		}
	}

	if (reading) {
		vals[idx++] = sign * num;
		if (is_node && idx == 3) {
			id = vals[0];
			g->nodes[id].x = (int32_t) vals[1];
			g->nodes[id].y = (int32_t) vals[2];
			g->nodes[id].active = true;
			g->nodes[id].time_from_source = 0;
			g->nodes[id].time_to_dest = 0;
			g->nodes[id].dist_from_source = 0;
			g->nodes[id].dist_to_dest = 0;
			node_i++;
		} else if (!is_node && idx == 4) {
			from = (uint32_t) vals[0];
			to = (uint32_t) vals[1];
			dist = (uint32_t) vals[2];
			time = (uint32_t) vals[3];
			addEdge(g, from, to, dist, time);
			if (edge_i < n_nodes) {
				g->nodes[edge_i].dist_from_source = haversine(g->nodes[s].x, g->nodes[s].y, g->nodes[edge_i].x,
				                                              g->nodes[edge_i].y);
				g->nodes[edge_i].dist_to_dest = haversine(g->nodes[d].x, g->nodes[d].y, g->nodes[edge_i].x,
				                                          g->nodes[edge_i].y);
			}
			edge_i++;
		}
	}
	fclose(fp);
	for (uint32_t i = n_edges; i < n_nodes; i++) {
		g->nodes[i].dist_from_source = haversine(g->nodes[s].x, g->nodes[s].y, g->nodes[i].x, g->nodes[i].y);
		g->nodes[i].dist_to_dest = haversine(g->nodes[d].x, g->nodes[d].y, g->nodes[i].x, g->nodes[i].y);
	}
	end_read = clock();
	printf("Graph loaded: %u nodes.\n", n_nodes);
	printf("Time read graph: %lf\n", ((double) ((double) (end_read - start_read) / CLOCKS_PER_SEC)));
	//-----------------------------------------END GRAPH READING
	//-----------------------------------------
	//-----------------------------------------
	//-----------------------------------------

	//-----------------------------------------START READ INSTANCE AND OUTPUT PREPARATION
	//-----------------------------------------
	//-----------------------------------------
	//-----------------------------------------
	FILE *input_instance;
	int num_ist = 1;
	if (multipleinstanceflag == 1) {
		input_instance = fopen(inputinstance, "r");
		if (!input_instance) {
			perror("Error opening file");
			return 1;
		}
		fscanf(input_instance, "%d\n", &num_ist);
	}

	char results_path[512];
	build_output_path(results_path, sizeof(results_path), output_dir, "results.txt");
	FILE *fout = fopen(results_path, "w");
	if (!fout) {
		perror("Error opening results file");
		return 1;
	}
	fprintf(
		fout,
		"Source Target Budget "
		"Cost_SPC Resource_SPC Cost_SPR Resource_SPR "
		"runtime_readGraph runtime_preprocess runtime_H "
		"Cost_H Resource_H "
		"BestIter BestLambda NumIter \n");
	//-----------------------------------------END READ INSTANCE AND OUTPUT PREPARATION
	//-----------------------------------------
	//-----------------------------------------
	//-----------------------------------------


	//-----------------------------------------EXECUTE INSTANCE
	//-----------------------------------------
	//-----------------------------------------
	//-----------------------------------------
	for (int ist = 0; ist < num_ist; ist++) {
		if (multipleinstanceflag == 1)
			fscanf(input_instance, "%u %u %lu\n", &s, &d, &W_read);
		printf("Instance %d source %d target %d budget %lu\n\n", ist + 1, s, d, W_read);

		char improvements_filename[512];
		char improvements_base[128];
		snprintf(improvements_base, sizeof(improvements_base),
		         "%" PRIu32 "-%" PRIu32 "-%" PRIu64 "_sol_details.txt",
		         s, d, (uint64_t) W_read);
		build_output_path(improvements_filename, sizeof(improvements_filename), output_dir, improvements_base);
		ImprovementList improvements;
		improvements_init(&improvements);

		double W = (double) W_read;
		for (uint32_t i = 0; i < n_nodes; i++) {
			g->nodes[i].active = true;
			g->nodes[i].time_from_source = 0;
			g->nodes[i].time_to_dest = 0;
			g->nodes[i].dist_from_source = haversine(g->nodes[s].x, g->nodes[s].y, g->nodes[i].x, g->nodes[i].y);
			g->nodes[i].dist_to_dest = haversine(g->nodes[d].x, g->nodes[d].y, g->nodes[i].x, g->nodes[i].y);
		}

		//-----------------------------------------START EXECUTION ALG
		//-----------------------------------------
		//-----------------------------------------
		//-----------------------------------------
		clock_t start_exec = clock();
		clock_t end_exec = clock();

		AStarResult path_inf = {0};
		AStarResult path_best = {0};
		uint32_t remaining = 0;

		//-----------------------------------------START PREPROCESS
		//-----------------------------------------
		//-----------------------------------------
		//-----------------------------------------
		clock_t start_time_preprocess = clock();
		clock_t end_time_preprocess = clock();

		AStarResult astar_spc = a_star_with_bound(g, s, d, DBL_MAX, DBL_MAX, true, s, d, W, DBL_MAX, 1, d);

		cost_spc = astar_spc.cost;
		resource_spc = astar_spc.resource;
		if (astar_spc.resource > W) {
			free_astar_result(&path_inf);
			path_inf = copy_res(&astar_spc);
			if (path_inf.resource == DBL_MAX) {
				end_time_preprocess = clock();
				end_exec = clock();
				fprintf(fout, "%u %u %.0lf 0 0 0 0 %.3f %.3f %.3f 0.0 0.0 0 0 0\n",
				        s, d, (double) W_read,
				        ((double) ((double) (end_read - start_read) / CLOCKS_PER_SEC)),
				        ((double) ((double) (end_time_preprocess - start_time_preprocess) / CLOCKS_PER_SEC)),
				        ((double) ((double) (end_time_preprocess - start_time_preprocess) / CLOCKS_PER_SEC)));
				improvements_write_file(improvements_filename, &improvements);
				improvements_free(&improvements);
				continue;
			}
		} else {
			free_astar_result(&path_best);
			path_best = copy_res(&astar_spc);
			improvements_push(&improvements, path_best.cost, path_best.resource,
			                  elapsed_seconds(start_exec), 1.0, 0);

			end_exec = clock();
			end_time_preprocess = clock();
			fprintf(fout, "%u %u %.0lf %.3f %.3f 0 0 %.3f %.3f %.3f 0.0 0.0 0 0 0\n",
			        s, d, (double) W_read,
			        cost_spc, resource_spc,
			        ((double) ((double) (end_read - start_read) / CLOCKS_PER_SEC)),
			        ((double) ((double) (end_time_preprocess - start_time_preprocess) / CLOCKS_PER_SEC)),
			        ((double) ((double) (end_time_preprocess - start_time_preprocess) / CLOCKS_PER_SEC)));

			improvements_write_file(improvements_filename, &improvements);
			improvements_free(&improvements);
			continue;
		}
		free_astar_result(&astar_spc);

		AStarResult astar_spr = a_star_with_bound(g, s, -1, W, DBL_MAX, true, s, d, W, DBL_MAX, 0, d); // forward, time

		cost_spr = astar_spr.cost;
		resource_spr = astar_spr.resource;

		if (resource_spr <= W) {
			free_astar_result(&path_best);
			path_best = copy_res(&astar_spr);
			improvements_push(&improvements, cost_spr, resource_spr,
			                  elapsed_seconds(start_exec), 0.0, 0);
		} else {
			fprintf(fout, "%u %u %.0lf %.3f %.3f 0 0 %.3f %.3f %.3f 0.0 0.0 0 0 0\n",
			        s, d, (double) W_read,
			        cost_spc, resource_spc,
			        ((double) ((double) (end_read - start_read) / CLOCKS_PER_SEC)),
			        ((double) ((double) (end_time_preprocess - start_time_preprocess) / CLOCKS_PER_SEC)),
			        ((double) ((double) (end_time_preprocess - start_time_preprocess) / CLOCKS_PER_SEC)));
			improvements_write_file(improvements_filename, &improvements);
			improvements_free(&improvements);
			continue;
		}
		free_astar_result(&astar_spr);

		// Count Remaining Nodes
		remaining = 0;
		for (uint32_t i = 0; i < g->n_nodes; i++)
			if (g->nodes[i].active)
				if (i != s && i != d) {
					if (g->nodes[i].time_to_dest + g->nodes[i].time_from_source > 0 && g->nodes[i].time_to_dest + g->
					    nodes[i].time_from_source <= W && g->nodes[i].dist_to_dest + g->nodes[i].dist_from_source <
					    path_best.cost)
						remaining++;
					else g->nodes[i].active = false;
					if (g->nodes[i].time_from_source == 0) g->nodes[i].active = false;
				}

		AStarResult astar_backward_reduction = a_star_with_bound(g, d, -1, W, path_best.cost, false, s, d, W,
		                                                         path_best.cost, 0, d); // forward, time
		free_astar_result(&astar_backward_reduction);

		// Count Remaining Nodes
		remaining = 0;
		for (uint32_t i = 0; i < g->n_nodes; i++)
			if (g->nodes[i].active)
				if (i != s && i != d) {
					if (g->nodes[i].time_to_dest + g->nodes[i].time_from_source > 0 && g->nodes[i].time_to_dest + g->
					    nodes[i].time_from_source <= W && g->nodes[i].dist_to_dest + g->nodes[i].dist_from_source <
					    path_best.cost)
						remaining++;
					else g->nodes[i].active = false;
					if (g->nodes[i].time_to_dest == 0) g->nodes[i].active = false;
				}

		AStarResult astar_forward_reduction_distance = a_star_with_bound(
			g, s, -1, W, path_best.cost, true, s, d, W, path_best.cost, 1, -1); // forward, time

		remaining = 0;
		for (uint32_t i = 0; i < g->n_nodes; i++)
			if (g->nodes[i].active)
				if (i != s && i != d) {
					if (g->nodes[i].time_to_dest + g->nodes[i].time_from_source > 0 && g->nodes[i].time_to_dest + g->
					    nodes[i].time_from_source <= W && g->nodes[i].dist_to_dest + g->nodes[i].dist_from_source <
					    path_best.cost)
						remaining++;
					else g->nodes[i].active = false;
				}
		free_astar_result(&astar_forward_reduction_distance);

		AStarResult astar_backward_reduction_distance = a_star_with_bound(
			g, d, -1, W, path_best.cost, false, s, d, W, path_best.cost, 1, -1); // forward, time
		free_astar_result(&astar_backward_reduction_distance);

		// Count Remaining Nodes
		remaining = 0;
		for (uint32_t i = 0; i < g->n_nodes; i++)
			if (g->nodes[i].active)
				if (i != s && i != d) {
					if (g->nodes[i].time_to_dest + g->nodes[i].time_from_source > 0 && g->nodes[i].time_to_dest + g->
					    nodes[i].time_from_source <= W && g->nodes[i].dist_to_dest + g->nodes[i].dist_from_source <
					    path_best.cost)
						remaining++;
					else g->nodes[i].active = false;
				}


		end_time_preprocess = clock();
		double preprocess_time = ((double) (end_time_preprocess - start_time_preprocess) / CLOCKS_PER_SEC);
		//-----------------------------------------END PREPROCESS
		//-----------------------------------------
		//-----------------------------------------
		//-----------------------------------------


		end_exec = clock();
		double time_best_solution = 0.0;
		int nIt = 0;
		int max_iterations_red = INF;
		int iter_red = 0;
		int iter_lambda = 0;
		int best_iter_lambda = 0;
		double best_lambda = 1;
		int check_connected = 1;
		while (path_inf.resource > W && elapsed_seconds(start_exec) <
		       time_limit && iter_red < max_iterations_red) {
			AStarResult astar_lambda = {0};
			if (iter_red > 0) {
				astar_lambda = a_star_with_bound(g, s, d, W, path_best.cost, true, s, d, W, path_best.cost, 0,
				                                 -1); // forward, time

				if (astar_lambda.resource == DBL_MAX) {
					free_astar_result(&astar_lambda);
					break;
				}

				if (astar_lambda.resource <= W)
					if (astar_lambda.cost < path_best.cost) {
						time_best_solution = elapsed_seconds(start_exec);
						if (time_best_solution > time_limit)
							break;

						free_astar_result(&path_best);
						path_best = copy_res(&astar_lambda);
						improvements_push(&improvements, path_best.cost, path_best.resource,
						                  time_best_solution, 0.0, nIt);
						best_iter_lambda = nIt;
						best_lambda = 0;
					}
				free_astar_result(&astar_lambda);
			}

			lambda = 0.5;
			double lambda_sup = 1, lambda_inf = 0;
			iter_lambda = 0;
			while (elapsed_seconds(start_exec) < time_limit &&
			       iter_lambda < max_iterations) {
				//Call lambda A*
				iter_lambda++;
				nIt++;
				free_astar_result(&astar_lambda);

				astar_lambda = a_star_with_bound(g, s, d, W, path_best.cost, true, s, d, W, path_best.cost,
				                                 lambda, -1); // forward, time

				if (astar_lambda.resource <= W) {
					time_best_solution = elapsed_seconds(start_exec);
					if (time_best_solution > time_limit) {
						check_connected = 0;
						break;
					}
					if (astar_lambda.cost < path_best.cost) {
						free_astar_result(&path_best);
						path_best = copy_res(&astar_lambda);
						improvements_push(&improvements, path_best.cost, path_best.resource,
						                  time_best_solution, lambda, nIt);
						best_iter_lambda = nIt;
						best_lambda = lambda;
						end_exec = clock();
					}
					lambda_inf = lambda;
					lambda = lambda_sup - (lambda_sup - lambda_inf) / 2;
				} else {
					lambda_sup = lambda;
					lambda = lambda_sup - (lambda_sup - lambda_inf) / 2;
				}
			}
			free_astar_result(&astar_lambda);

			if (check_connected == 0 || run_reduction_heuristic == 0)
				break;

			free_astar_result(&path_inf);
			path_inf = a_star_with_bound(g, s, d, W, path_best.cost, true, s, d, W, path_best.cost, 1, -1);

			if (path_inf.resource > W) {
				if (path_inf.resource == DBL_MAX) {
					break;
				}
				if (perc_red > 0.0 && perc_red <= 1.0) {
					W = path_inf.bound_path - (perc_red * path_inf.bound_path);
				} else {
					W = path_inf.bound_path - 1e-06;
				}
				remaining = 0;
				for (uint32_t i = 0; i < g->n_nodes; i++)
					if (i != s && i != d) {
						if (g->nodes[i].active && g->nodes[i].time_to_dest + g->nodes[i].time_from_source > 0 && g->
						    nodes[i].time_to_dest + g->nodes[i].time_from_source < W && g->nodes[i].dist_to_dest + g->
						    nodes[i].dist_from_source < path_best.cost)
							remaining++;
						else g->nodes[i].active = false;
					}
			} else {
				if (path_inf.cost < path_best.cost) {
					time_best_solution = elapsed_seconds(start_exec);
					if (time_best_solution > time_limit)
						break;
					free_astar_result(&path_best);
					path_best = copy_res(&path_inf);
					improvements_push(&improvements, path_best.cost, path_best.resource,
					                  time_best_solution, 1.0, nIt);
					best_iter_lambda = nIt;
					best_lambda = 1;
					end_exec = clock();
					time_best_solution = ((double) ((double) (end_exec - start_exec) / CLOCKS_PER_SEC));
				}
				break;
			}
			iter_red++;
		}
		end_exec = clock();

		fprintf(fout, "%u %u %.0lf %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %d %.3f %d\n",
		        s, d, (double) W_read,
		        cost_spc, resource_spc, cost_spr, resource_spr,
		        ((double) ((double) (end_read - start_read) / CLOCKS_PER_SEC)),
		        preprocess_time,
		        ((double) ((double) (end_exec - start_exec) / CLOCKS_PER_SEC)),
		        path_best.cost, path_best.resource, best_iter_lambda, best_lambda, nIt);

		fflush(fout);
		improvements_write_file(improvements_filename, &improvements);
		improvements_free(&improvements);

		// ---------- Cleanup ----------
		free_astar_result(&path_inf);
		free_astar_result(&path_best);
		if (multipleinstanceflag != 1)
			return 0;
	}

	fclose(fout);
	freeGraph(g);
	if (multipleinstanceflag == 1)
		fclose(input_instance);
	return 0;
}
