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

static void larac_improvements_write_file(const char *filename, const ImprovementList *list) {
	FILE *f = fopen(filename, "w");
	if (!f) {
		perror("Error opening LARAC improvements file");
		return;
	}
	fprintf(f, "cost resource timeFound iter\n");
	for (size_t i = 0; i < list->count; i++) {
		const Improvement *imp = &list->items[i];
		fprintf(f, "%.0f %.0f %.6f %d\n",
		        imp->cost, imp->resource, imp->timeFound, imp->iter);
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

// ---------- LARAC (Lagrangian Relaxation based Aggregated Cost) ----------
//
// Reference: A. Juttner, B. Szviatovszki, I. Mecs, Zs. Rajko, "Lagrange relaxation based
// method for the QoS routing problem", IEEE INFOCOM 2001 (see the PDF included in this repo).
//
// LARAC solves the same problem as PARC:
//     minimize distance(P)   subject to   time(P) <= W ,  P a path from s to d.
// It is used here purely as an independent baseline to compare against PARC: it always runs
// a plain, textbook Dijkstra on the *full* graph (it never looks at g->nodes[].active and does
// not use PARC's bidirectional bounding/graph-reduction machinery), so a difference in runtime
// or solution quality between the two reflects the algorithms themselves, not shared pruning.
//
// --- The idea in one paragraph ---
// For any lambda >= 0, define the "aggregated cost" of an edge as distance(e) + lambda*time(e).
// The path P(lambda) minimizing the sum of aggregated costs is a Lagrangian relaxation of the
// budgeted problem: for every path P, distance(P(lambda)) + lambda*time(P(lambda)) is a valid
// lower bound on the true optimum whenever time(P) <= W is enforced. LARAC searches for the
// multiplier lambda that makes this bound as tight as possible, using the classical secant
// ("line search") update: given a cost-optimal but infeasible path P_c (over budget) and a
// time-optimal feasible path P_w, it computes the lambda that makes P_c and P_w equally good
// under the aggregated cost (the slope of the segment joining them in the (time,cost) plane),
// then resolves the aggregated shortest path P_r under that lambda:
//   - if P_r cannot improve on P_c/P_w under the aggregated cost, the line search has converged:
//     P_w is returned as the LARAC solution;
//   - otherwise P_r replaces whichever of P_c (if still infeasible) or P_w (if feasible) it
//     dominates, and the process repeats.
// This is guaranteed to terminate (finitely many distinct paths), always returns a feasible
// path, and the returned path is *exactly* optimal whenever the duality gap is zero -- but on
// integer-weighted graphs a strictly positive integrality gap can remain, in which case LARAC's
// result is a heuristic (usually very close to optimal, but not certified). This is the key
// qualitative difference to keep in mind when comparing against PARC, which instead keeps
// refining/branching until it can certify optimality (subject to its own time/iteration limits).

// Plain single-pair Dijkstra minimizing alpha*distance(e) + beta*time(e) summed over the path.
// Returns the *raw* distance and time sums of the shortest path found in res.cost/res.resource
// respectively (not the aggregated alpha/beta cost), so callers can test feasibility directly.
// Uses lazy deletion (stale heap entries are skipped via the cost[] check) instead of an
// "active" flag, and never restricts itself to a subgraph: this is intentionally the simplest
// correct Dijkstra, independent of any of PARC's search-space reduction logic.
AStarResult dijkstra_linear(Graph *g, uint32_t s, uint32_t d, double alpha, double beta) {
	uint32_t n = g->n_nodes;
	double *cost = malloc(sizeof(double) * n); // aggregated cost alpha*distance + beta*time
	double *dist_acc = malloc(sizeof(double) * n); // raw distance accumulated along the path
	double *time_acc = malloc(sizeof(double) * n); // raw time accumulated along the path
	uint32_t *parent = malloc(sizeof(uint32_t) * n);
	for (uint32_t i = 0; i < n; i++) {
		cost[i] = DBL_MAX;
		dist_acc[i] = DBL_MAX;
		time_acc[i] = DBL_MAX;
		parent[i] = UINT32_MAX;
	}

	MinHeap *heap = createHeap(n);
	cost[s] = 0;
	dist_acc[s] = 0;
	time_acc[s] = 0;
	push(heap, s, 0);

	while (!isEmpty(heap)) {
		HeapNode curr = pop(heap);
		uint32_t u = curr.node;
		if (curr.cost > cost[u]) continue; // stale entry left by an earlier, worse push: skip it
		if (u == d) break; // shortest path to d is finalized as soon as it is popped

		for (uint32_t i = 0; i < g->nodes[u].out_count; i++) {
			Edge e = g->nodes[u].out_edges[i];
			uint32_t v = e.to;
			double new_cost = cost[u] + alpha * e.distance + beta * e.time;
			if (new_cost < cost[v]) {
				cost[v] = new_cost;
				dist_acc[v] = dist_acc[u] + e.distance;
				time_acc[v] = time_acc[u] + e.time;
				parent[v] = u;
				push(heap, v, new_cost);
			}
		}
	}

	AStarResult res;
	res.path = malloc(sizeof(uint32_t) * n);
	res.path_len = 0;
	res.bound_path = 0;
	if (cost[d] == DBL_MAX) {
		res.cost = DBL_MAX;
		res.resource = DBL_MAX;
	} else {
		res.cost = dist_acc[d];
		res.resource = time_acc[d];
		uint32_t v = d;
		while (v != s) {
			res.path[res.path_len++] = v;
			v = parent[v];
		}
		res.path[res.path_len++] = s;
	}

	free(cost);
	free(dist_acc);
	free(time_acc);
	free(parent);
	free(heap->data);
	free(heap);
	return res;
}

// Runs the LARAC heuristic for the WCSPP from s to d under time budget W.
// *out_iterations receives the number of aggregated-cost Dijkstra calls used after the initial
// two probes (P_c, P_w) -- this is the LARAC analogue of PARC's lambda-search iteration count.
// *out_optimal is set to 1 when the returned path is *certified* optimal (either P_c was already
// feasible, the instance is infeasible, or the duality gap closed to zero), and 0 when LARAC
// stopped only because of the iteration/time limit, meaning a (typically small) optimality gap
// may remain. res.cost == DBL_MAX signals a certified-infeasible instance (no path respects W).
AStarResult larac(Graph *g, uint32_t s, uint32_t d, double W,
                   int max_iterations, double time_limit, clock_t start_time,
                   int *out_iterations, int *out_optimal, double *out_time_best,
                   ImprovementList *out_improvements) {
	*out_iterations = 0;
	*out_optimal = 0;
	*out_time_best = -1.0;

	// P_c: the cheapest path in distance, ignoring time entirely.
	AStarResult P_c = dijkstra_linear(g, s, d, 1.0, 0.0);
	if (P_c.resource == DBL_MAX) {
		// s and d are not connected at all.
		*out_optimal = 1;
		return P_c;
	}
	if (P_c.resource <= W) {
		// The cheapest possible path already respects the budget: nothing can beat it.
		*out_time_best = elapsed_seconds(start_time);
		if (out_improvements)
			improvements_push(out_improvements, P_c.cost, P_c.resource,
			                  *out_time_best, 0.0, 0);
		*out_optimal = 1;
		return P_c;
	}

	// P_w: the fastest path in time, ignoring distance entirely.
	AStarResult P_w = dijkstra_linear(g, s, d, 0.0, 1.0);
	if (P_w.resource <= W) {
		*out_time_best = elapsed_seconds(start_time);
		if (out_improvements)
			improvements_push(out_improvements, P_w.cost, P_w.resource,
			                  *out_time_best, 0.0, 0);
	}
	if (P_w.resource > W) {
		// Even the fastest possible path violates the budget: the instance is infeasible.
		free_astar_result(&P_c);
		AStarResult infeasible = P_w;
		infeasible.cost = DBL_MAX;
		infeasible.resource = DBL_MAX;
		*out_optimal = 1;
		return infeasible;
	}

	while (*out_iterations < max_iterations && elapsed_seconds(start_time) < time_limit) {
		(*out_iterations)++;

		// Slope of the segment joining P_c and P_w in the (time, cost) plane: the Lagrange
		// multiplier making both paths equally good under the aggregated cost distance+lambda*time.
		double lambda = (P_c.cost - P_w.cost) / (P_w.resource - P_c.resource);

		AStarResult P_r = dijkstra_linear(g, s, d, 1.0, lambda);
		double aggregated_r = P_r.cost + lambda * P_r.resource;
		double aggregated_w = P_w.cost + lambda * P_w.resource; // == P_c.cost + lambda*P_c.resource

		// P_r is a shortest path under the aggregated cost, so aggregated_r can never exceed
		// aggregated_w/aggregated_c; if it cannot go strictly below it either, the line search
		// has converged and no further improvement is reachable via this lambda scan.
		if (aggregated_r >= aggregated_w - 1e-9) {
			free_astar_result(&P_c);
			free_astar_result(&P_r);
			*out_optimal = 1;
			return P_w;
		}

		if (P_r.resource <= W) {
			free_astar_result(&P_w);
			P_w = P_r; // strictly better feasible path: tighten the upper bound
			*out_time_best = elapsed_seconds(start_time);
			if (out_improvements)
				improvements_push(out_improvements, P_w.cost, P_w.resource,
				                  *out_time_best, 0.0, *out_iterations);
		} else {
			free_astar_result(&P_c);
			P_c = P_r; // strictly better infeasible path: tighten the lower bound
		}
	}

	// Iteration/time budget exhausted before the gap closed: return the best feasible path
	// found so far. out_optimal stays 0 to flag that this result is not certified optimal.
	free_astar_result(&P_c);
	return P_w;
}


// ---------- APULSE integrated on PARC graph ----------
typedef enum { APULSE_FEASIBLE, APULSE_INFEASIBLE, APULSE_TIMEOUT, APULSE_ERROR } ApulseStatus;

typedef struct {
    ApulseStatus status;
    double cost, resource, delta_t, heuristic_time, total_time;
    uint64_t expanded, generated, retained_states;
    uint64_t pruned_feasibility, pruned_incumbent, pruned_bucket;
    uint64_t settled_resource, settled_objective;
} ApulseResult;

typedef struct {
    double f, cost, resource;
    uint32_t node;
    uint64_t serial;
} ApulseLabel;

typedef struct {
    ApulseLabel *data;
    size_t size, cap;
} ApulseHeap;

static int apulse_label_less(const ApulseLabel *a, const ApulseLabel *b) {
    const double eps = 1e-9;
    if (fabs(a->f - b->f) > eps) return a->f < b->f;
    if (fabs(a->cost - b->cost) > eps) return a->cost < b->cost;
    if (fabs(a->resource - b->resource) > eps) return a->resource < b->resource;
    return a->serial < b->serial;
}

static void apulse_heap_push(ApulseHeap *h, ApulseLabel x) {
    if (h->size == h->cap) {
        h->cap = h->cap ? h->cap * 2 : 1024;
        h->data = realloc(h->data, h->cap * sizeof(*h->data));
        if (!h->data) { perror("realloc APULSE heap"); exit(EXIT_FAILURE); }
    }
    size_t i = h->size++;
    h->data[i] = x;
    while (i) {
        size_t parent = (i - 1) / 2;
        if (!apulse_label_less(&h->data[i], &h->data[parent])) break;
        ApulseLabel t = h->data[i]; h->data[i] = h->data[parent]; h->data[parent] = t;
        i = parent;
    }
}

static ApulseLabel apulse_heap_pop(ApulseHeap *h) {
    ApulseLabel root = h->data[0];
    h->data[0] = h->data[--h->size];
    size_t i = 0;
    for (;;) {
        size_t l = 2*i+1, r = l+1, m=i;
        if (l < h->size && apulse_label_less(&h->data[l], &h->data[m])) m=l;
        if (r < h->size && apulse_label_less(&h->data[r], &h->data[m])) m=r;
        if (m == i) break;
        ApulseLabel t=h->data[i]; h->data[i]=h->data[m]; h->data[m]=t; i=m;
    }
    return root;
}

typedef struct {
    uint64_t *keys;
    double *values;
    unsigned char *used;
    size_t cap, size;
} ApulseMap;

static uint64_t apulse_hash64(uint64_t x) {
    x ^= x >> 30; x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27; x *= UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

static void apulse_map_init(ApulseMap *m, size_t cap) {
    size_t c=1; while (c < cap) c <<= 1;
    m->cap=c; m->size=0;
    m->keys=malloc(c*sizeof(uint64_t)); m->values=malloc(c*sizeof(double)); m->used=calloc(c,1);
    if(!m->keys||!m->values||!m->used){perror("malloc APULSE map");exit(EXIT_FAILURE);}
}
static void apulse_map_free(ApulseMap *m){free(m->keys);free(m->values);free(m->used);memset(m,0,sizeof(*m));}
static void apulse_map_rehash(ApulseMap *m) {
    ApulseMap n; apulse_map_init(&n, m->cap*2);
    for(size_t i=0;i<m->cap;i++) if(m->used[i]){
        size_t j=(size_t)apulse_hash64(m->keys[i])&(n.cap-1);
        while(n.used[j]) j=(j+1)&(n.cap-1);
        n.used[j]=1;n.keys[j]=m->keys[i];n.values[j]=m->values[i];n.size++;
    }
    apulse_map_free(m); *m=n;
}
static int apulse_map_get(ApulseMap *m,uint64_t key,double *value){
    size_t i=(size_t)apulse_hash64(key)&(m->cap-1);
    while(m->used[i]){if(m->keys[i]==key){*value=m->values[i];return 1;}i=(i+1)&(m->cap-1);}return 0;
}
static void apulse_map_put(ApulseMap *m,uint64_t key,double value){
    if((m->size+1)*10 >= m->cap*7) apulse_map_rehash(m);
    size_t i=(size_t)apulse_hash64(key)&(m->cap-1);
    while(m->used[i]&&m->keys[i]!=key)i=(i+1)&(m->cap-1);
    if(!m->used[i]){m->used[i]=1;m->keys[i]=key;m->size++;}m->values[i]=value;
}

static int apulse_reverse_dijkstra(Graph *g, uint32_t target, int use_resource,
                                   double *dist, clock_t start, double limit, uint64_t *settled) {
    for(uint32_t i=0;i<g->n_nodes;i++) dist[i]=DBL_MAX;
    MinHeap *heap=createHeap(g->n_nodes); dist[target]=0.0; push(heap,target,0.0); *settled=0;
    while(!isEmpty(heap)){
        HeapNode cur=pop(heap); uint32_t u=cur.node;
        if(cur.cost > dist[u] + 1e-9) continue;
        (*settled)++;
        if(((*settled)&4095ULL)==0 && limit>0.0 && elapsed_seconds(start)>=limit){free(heap->data);free(heap);return 0;}
        for(uint32_t i=0;i<g->nodes[u].in_count;i++){
            InEdge e=g->nodes[u].in_edges[i]; double w=use_resource?e.time:e.distance; double nd=cur.cost+w;
            if(nd+1e-9<dist[e.from]){dist[e.from]=nd;push(heap,e.from,nd);}
        }
    }
    free(heap->data);free(heap);return 1;
}

static ApulseResult apulse(Graph *g,uint32_t source,uint32_t target,double budget,
                           uint32_t target_buckets,double min_bucket_width,double time_limit){
    ApulseResult r; memset(&r,0,sizeof(r)); r.status=APULSE_ERROR; r.cost=r.resource=DBL_MAX;
    clock_t start=clock();
    r.delta_t=fmax(min_bucket_width,budget/(double)target_buckets);
    if(r.delta_t<=0.0) r.delta_t=DBL_MIN;
    double *hres=malloc(g->n_nodes*sizeof(double)), *hobj=malloc(g->n_nodes*sizeof(double));
    if(!hres||!hobj){perror("malloc APULSE heuristics");exit(EXIT_FAILURE);}
    if(!apulse_reverse_dijkstra(g,target,1,hres,start,time_limit,&r.settled_resource) ||
       !apulse_reverse_dijkstra(g,target,0,hobj,start,time_limit,&r.settled_objective)){
        r.status=APULSE_TIMEOUT;r.total_time=elapsed_seconds(start);free(hres);free(hobj);return r;
    }
    r.heuristic_time=elapsed_seconds(start);
    if(hres[source]==DBL_MAX||hres[source]>budget+1e-9||hobj[source]==DBL_MAX){
        r.status=APULSE_INFEASIBLE;r.total_time=elapsed_seconds(start);free(hres);free(hobj);return r;
    }
    ApulseHeap q={0}; ApulseMap visited; apulse_map_init(&visited,1u<<20); uint64_t serial=0;
    apulse_heap_push(&q,(ApulseLabel){hobj[source],0.0,0.0,source,serial++});
    double incumbent=DBL_MAX,inc_res=DBL_MAX;

    /* The search is now valid and running. APULSE_ERROR is reserved only for
       an actual internal error, such as a bucket-index overflow. */
    r.status=APULSE_INFEASIBLE;

    while(q.size){
        if((r.expanded&4095ULL)==0 && time_limit>0.0 && elapsed_seconds(start)>=time_limit){r.status=APULSE_TIMEOUT;break;}
        ApulseLabel cur=apulse_heap_pop(&q);
        if(cur.resource+hres[cur.node]>budget+1e-9){r.pruned_feasibility++;continue;}
        if(cur.f>=incumbent-1e-9){r.pruned_incumbent++;continue;}
        double bd=floor(fmax(0.0,cur.resource)/r.delta_t);
        if(bd>(double)UINT32_MAX){r.status=APULSE_ERROR;break;}
        uint64_t key=((uint64_t)(uint32_t)bd<<32)|cur.node; double old;
        if(apulse_map_get(&visited,key,&old)&&cur.cost>=old-1e-9){r.pruned_bucket++;continue;}
        apulse_map_put(&visited,key,cur.cost);r.expanded++;
        if(cur.node==target){
            incumbent=cur.cost;
            inc_res=cur.resource;
            r.cost=incumbent;
            r.resource=inc_res;
            printf("  NEW_INCUMBENT source=%" PRIu32 " target=%" PRIu32
                   " budget=%.0f cost=%.6f resource=%.6f time=%.6f expanded=%" PRIu64 "\n",
                   source, target, budget, incumbent, inc_res,
                   elapsed_seconds(start), r.expanded);
            fflush(stdout);
            continue;
        }
        for(uint32_t i=0;i<g->nodes[cur.node].out_count;i++){
            Edge e=g->nodes[cur.node].out_edges[i];
            if(hres[e.to]==DBL_MAX||hobj[e.to]==DBL_MAX)continue;
            double nr=cur.resource+e.time,nc=cur.cost+e.distance;
            apulse_heap_push(&q,(ApulseLabel){nc+hobj[e.to],nc,nr,e.to,serial++});r.generated++;
        }
    }
    r.total_time=elapsed_seconds(start);r.retained_states=visited.size;
    if(r.status==APULSE_ERROR){
        /* Preserve ERROR only for a genuine internal failure. */
    } else if(incumbent<DBL_MAX){
        r.cost=incumbent;
        r.resource=inc_res;
        /* A timeout remains visible as TIMEOUT, but the incumbent is saved.
           A normally exhausted queue is reported as FEASIBLE. */
        if(r.status!=APULSE_TIMEOUT) r.status=APULSE_FEASIBLE;
    } else if(r.status!=APULSE_TIMEOUT){
        r.status=APULSE_INFEASIBLE;
    }
    free(q.data);apulse_map_free(&visited);free(hres);free(hobj);return r;
}

static const char *apulse_status_name(ApulseStatus s){
    switch(s){case APULSE_FEASIBLE:return "FEASIBLE";case APULSE_INFEASIBLE:return "INFEASIBLE";case APULSE_TIMEOUT:return "TIMEOUT";default:return "ERROR";}
}

static void write_results_header(FILE *fout, int run_parc, int run_larac, int run_apulse) {
    fprintf(fout, "Source Target Budget");
    if (run_parc) fprintf(fout, " Cost_SPC Resource_SPC Cost_SPR Resource_SPR runtime_readGraph runtime_preprocess runtime_PARC Cost_PARC Resource_PARC BestIter BestLambda NumIter");
    if (run_larac) fprintf(fout, " runtime_LARAC TimeBest_LARAC Cost_LARAC Resource_LARAC Iter_LARAC Optimal_LARAC");
    if (run_apulse) fprintf(fout, " Status_APULSE N_APULSE DeltaT_APULSE runtime_heuristics_APULSE runtime_APULSE Cost_APULSE Resource_APULSE Expanded_APULSE Generated_APULSE RetainedStates_APULSE");
    fputc('\n', fout);
}

static void write_results_row(FILE *fout,
                              int run_parc, int run_larac, int run_apulse,
                              uint32_t s, uint32_t d, uint64_t budget,
                              double cost_spc, double resource_spc, double cost_spr, double resource_spr,
                              double runtime_read_graph, double runtime_preprocess, double runtime_parc,
                              double cost_parc, double resource_parc, int best_iter, double best_lambda, int num_iter,
                              double runtime_larac, double time_best_larac, double cost_larac, double resource_larac, int iter_larac, int optimal_larac,
                              const ApulseResult *ar, uint32_t apulse_n) {
    fprintf(fout, "%" PRIu32 " %" PRIu32 " %" PRIu64, s, d, budget);
    if (run_parc) fprintf(fout, " %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %d %.3f %d", cost_spc,resource_spc,cost_spr,resource_spr,runtime_read_graph,runtime_preprocess,runtime_parc,cost_parc,resource_parc,best_iter,best_lambda,num_iter);
    if (run_larac) fprintf(fout, " %.3f %.3f %.3f %.3f %d %d",runtime_larac,time_best_larac,cost_larac,resource_larac,iter_larac,optimal_larac);
    if (run_apulse && ar) fprintf(fout, " %s %u %.9f %.6f %.6f %.3f %.3f %" PRIu64 " %" PRIu64 " %" PRIu64,
        apulse_status_name(ar->status), apulse_n, ar->delta_t, ar->heuristic_time, ar->total_time,
        ar->cost==DBL_MAX?-1.0:ar->cost, ar->resource==DBL_MAX?-1.0:ar->resource,
        ar->expanded, ar->generated, ar->retained_states);
    fputc('\n',fout);
}

// ---------- Parser ----------
void parse_args(int argc, char *argv[],
                char **input_path, uint32_t *s, uint32_t *d,
                uint64_t *W, double *time_limit, int *run_reduction_heuristic, int *max_iterations,
                int *multipleinstanceflag, char **inputinstance, char **output_dir, double *perc_red,
                int *run_larac, int *run_parc, int *run_apulse, uint32_t *apulse_n, double *apulse_min_bucket_width) {
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
	*run_larac = 1;
	*run_parc = 1;
	*run_apulse = 0;
	*apulse_n = 8192;
	*apulse_min_bucket_width = 1.0;
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
		} else if (!strcmp(argv[i], "--larac")) {
			*run_larac = atoi(argv[++i]);
		} else if (!strcmp(argv[i], "--parc")) {
			*run_parc = atoi(argv[++i]);
		} else if (!strcmp(argv[i], "--apulse")) {
			*run_apulse = atoi(argv[++i]);
		} else if (!strcmp(argv[i], "--N")) {
			*apulse_n = (uint32_t) strtoul(argv[++i], NULL, 10);
		} else if (!strcmp(argv[i], "--min-bucket-width")) {
			*apulse_min_bucket_width = atof(argv[++i]);
		} else {
			fprintf(stderr, "Unknown argument: %s\n", argv[i]);
			exit(EXIT_FAILURE);
		}
	}

	if ((*run_larac != 0 && *run_larac != 1) || (*run_parc != 0 && *run_parc != 1) || (*run_apulse != 0 && *run_apulse != 1)) {
		fprintf(stderr, "Error: --larac, --parc and --apulse must be either 0 or 1.\n");
		exit(EXIT_FAILURE);
	}

	if (*run_larac == 0 && *run_parc == 0 && *run_apulse == 0) {
		fprintf(stderr, "Error: at least one algorithm must be enabled.\n");
		exit(EXIT_FAILURE);
	}

	if (!*input_path) {
		fprintf(stderr,
		        "Usage: ./compute_paths "
		        "--input file --s src --d dest --W bound "
		        "[--tl time_limit_seconds] [--nit max_iterations] [--multipleinstanceflag run all instances] "
		        "[--inputinstance file run all instances] [--outdir output_directory] "
		        "[--perc_red percent_reduction] [--larac 0|1] [--parc 0|1] [--apulse 0|1] [--N buckets] [--min-bucket-width value]\n"
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
	int run_larac;
	int run_parc;
	int run_apulse;
	uint32_t apulse_n;
	double apulse_min_bucket_width;
	double cost_spr, resource_spr, cost_spc, resource_spc;

	parse_args(argc, argv, &filename, &s, &d, &W_read, &time_limit, &run_reduction_heuristic,
	           &max_iterations, &multipleinstanceflag, &inputinstance, &output_dir, &perc_red,
	           &run_larac, &run_parc, &run_apulse, &apulse_n, &apulse_min_bucket_width);

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
	write_results_header(fout, run_parc, run_larac, run_apulse);
	//-----------------------------------------END READ INSTANCE AND OUTPUT PREPARATION
	//-----------------------------------------
	//-----------------------------------------
	//-----------------------------------------


	//-----------------------------------------EXECUTE INSTANCE APULSE
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

		char larac_improvements_filename[512];
		char larac_improvements_base[128];
		snprintf(larac_improvements_base, sizeof(larac_improvements_base),
		         "%" PRIu32 "-%" PRIu32 "-%" PRIu64 "_larac_sol_details.txt",
		         s, d, (uint64_t) W_read);
		build_output_path(larac_improvements_filename, sizeof(larac_improvements_filename),
		                  output_dir, larac_improvements_base);
		ImprovementList larac_improvements;
		improvements_init(&larac_improvements);

		double W = (double) W_read;
		for (uint32_t i = 0; i < n_nodes; i++) {
			g->nodes[i].active = true;
			g->nodes[i].time_from_source = 0;
			g->nodes[i].time_to_dest = 0;
			g->nodes[i].dist_from_source = haversine(g->nodes[s].x, g->nodes[s].y, g->nodes[i].x, g->nodes[i].y);
			g->nodes[i].dist_to_dest = haversine(g->nodes[d].x, g->nodes[d].y, g->nodes[i].x, g->nodes[i].y);
		}

		ApulseResult apulse_res; memset(&apulse_res, 0, sizeof(apulse_res));
		apulse_res.status = APULSE_ERROR; apulse_res.cost = apulse_res.resource = DBL_MAX;
		if (run_apulse) {
			/* APULSE has an independent per-instance timer, reset inside apulse(). */
			printf("  Starting APULSE with an independent time limit of %.3f s\n", time_limit);
			apulse_res = apulse(g, s, d, W, apulse_n, apulse_min_bucket_width, time_limit);
			printf("  APULSE status=%s cost=%.3f resource=%.3f time=%.6f\n",
			       apulse_status_name(apulse_res.status),
			       apulse_res.cost == DBL_MAX ? -1.0 : apulse_res.cost,
			       apulse_res.resource == DBL_MAX ? -1.0 : apulse_res.resource, apulse_res.total_time);
		}

		//-----------------------------------------START LARAC (comparison baseline)
		//-----------------------------------------
		//-----------------------------------------
		//-----------------------------------------
		// Computed once per instance, independently of PARC: it ignores g->nodes[].active,
		// so it does not matter whether it runs before or after PARC's own search/reduction.
		// It gets its own clock so its runtime is never mixed into PARC's timings above.
		double cost_larac = -1, resource_larac = -1, runtime_larac = 0.0;
		double time_best_larac = -1.0;
		int iter_larac = 0, optimal_larac = 0;
		if (run_larac) {
			/* Reset the time limit specifically for LARAC on this instance. */
			clock_t start_larac = clock();
			printf("  Starting LARAC with an independent time limit of %.3f s\n", time_limit);
			int larac_max_iterations = max_iterations > 0 ? max_iterations : 1000;
			AStarResult larac_res = larac(g, s, d, W, larac_max_iterations, time_limit,
			                              start_larac, &iter_larac, &optimal_larac,
			                              &time_best_larac, &larac_improvements);
			runtime_larac = elapsed_seconds(start_larac);
			if (larac_res.cost != DBL_MAX) {
				cost_larac = larac_res.cost;
				resource_larac = larac_res.resource;
			}
			free_astar_result(&larac_res);
			larac_improvements_write_file(larac_improvements_filename, &larac_improvements);
		}
		//-----------------------------------------END LARAC
		//-----------------------------------------
		//-----------------------------------------
		//-----------------------------------------

		if (!run_parc) {
			write_results_row(fout, run_parc, run_larac, run_apulse,
			                  s, d, W_read,
			                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			                  runtime_larac, time_best_larac, cost_larac, resource_larac,
			                  iter_larac, optimal_larac, &apulse_res, apulse_n);
			fflush(fout);
			improvements_free(&improvements);
			improvements_free(&larac_improvements);
			continue;
		}

		//-----------------------------------------START EXECUTION ALG
		//-----------------------------------------
		//-----------------------------------------
		//-----------------------------------------
		/* Reset the time limit specifically for PARC on this instance.
		 * Time spent by APULSE and/or LARAC is not included. */
		clock_t start_exec = clock();
		clock_t end_exec = start_exec;
		printf("  Starting PARC with an independent time limit of %.3f s\n", time_limit);

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
				write_results_row(fout, run_parc, run_larac, run_apulse,
				                  s, d, W_read,
				                  0, 0, 0, 0,
				                  ((double) (end_read - start_read) / CLOCKS_PER_SEC),
				                  ((double) (end_time_preprocess - start_time_preprocess) / CLOCKS_PER_SEC),
				                  ((double) (end_exec - start_exec) / CLOCKS_PER_SEC),
				                  0, 0, 0, 0, 0,
				                  runtime_larac, time_best_larac, cost_larac, resource_larac,
				                  iter_larac, optimal_larac, &apulse_res, apulse_n);
				improvements_write_file(improvements_filename, &improvements);
				improvements_free(&improvements);
			improvements_free(&larac_improvements);
				continue;
			}
		} else {
			free_astar_result(&path_best);
			path_best = copy_res(&astar_spc);
			improvements_push(&improvements, path_best.cost, path_best.resource,
			                  elapsed_seconds(start_exec), 1.0, 0);

			end_exec = clock();
			end_time_preprocess = clock();
			write_results_row(fout, run_parc, run_larac, run_apulse,
			                  s, d, W_read,
			                  cost_spc, resource_spc, 0, 0,
			                  ((double) (end_read - start_read) / CLOCKS_PER_SEC),
			                  ((double) (end_time_preprocess - start_time_preprocess) / CLOCKS_PER_SEC),
			                  ((double) (clock() - start_exec) / CLOCKS_PER_SEC),
			                  0, 0, 0, 0, 0,
			                  runtime_larac, time_best_larac, cost_larac, resource_larac,
			                  iter_larac, optimal_larac, &apulse_res, apulse_n);

			improvements_write_file(improvements_filename, &improvements);
			improvements_free(&improvements);
			improvements_free(&larac_improvements);
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
			write_results_row(fout, run_parc, run_larac, run_apulse,
			                  s, d, W_read,
			                  cost_spc, resource_spc, 0, 0,
			                  ((double) (end_read - start_read) / CLOCKS_PER_SEC),
			                  ((double) (end_time_preprocess - start_time_preprocess) / CLOCKS_PER_SEC),
			                  ((double) (clock() - start_exec) / CLOCKS_PER_SEC),
			                  0, 0, 0, 0, 0,
			                  runtime_larac, time_best_larac, cost_larac, resource_larac,
			                  iter_larac, optimal_larac, &apulse_res, apulse_n);
			improvements_write_file(improvements_filename, &improvements);
			improvements_free(&improvements);
			improvements_free(&larac_improvements);
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

		write_results_row(fout, run_parc, run_larac, run_apulse,
		                  s, d, W_read,
		                  cost_spc, resource_spc, cost_spr, resource_spr,
		                  ((double) (end_read - start_read) / CLOCKS_PER_SEC),
		                  preprocess_time,
		                  ((double) (end_exec - start_exec) / CLOCKS_PER_SEC),
		                  path_best.cost, path_best.resource,
		                  best_iter_lambda, best_lambda, nIt,
		                  runtime_larac, time_best_larac, cost_larac, resource_larac,
		                  iter_larac, optimal_larac, &apulse_res, apulse_n);

		fflush(fout);
		improvements_write_file(improvements_filename, &improvements);
		improvements_free(&improvements);
			improvements_free(&larac_improvements);

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
