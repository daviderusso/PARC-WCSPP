# PARC-WCSPP

Pruning-Aided Resource-Constrained Search Algorithm (PARC) for the Weight Constrained Shortest Path Problem (WCSPP).

This repository contains a C implementation of the algorithm (`src/parc.c`) and a helper script (`run_parc.sh`) to run batches of instances on different maps.

## Requirements

- `gcc` with C11 support
- `bash`

## Build

Build the executable manually:

```bash
gcc -O3 src/parc.c -o parc -lm
```

The script `run_parc.sh` will compile automatically if `./parc` is missing.

## Quick Run

Run the batch script with the default maps and instance lists:

```bash
bash run_parc.sh
```

Run a single instance directly:

```bash
./parc --input map/W.xy --s 2191918 --d 1689718 --W 7501374 --tl 60 --nit 10 --redh 1 --outdir results/tmp
```

Run multiple instances from a list:

```bash
./parc --input map/W.xy --multipleinstanceflag 1 --inputinstance instances/instances_W.txt --tl 60 --nit 10 --redh 1 --outdir results/tmp
```

## Input Data

The banchmark maps comes from the three largest road networks from the 9th DIMACS implementation challenge available at http://www.dis.uniroma1.it/~challenge9

### Map files (`map/*.txt`, `map/*.xy`)

Text format with a header followed by node and edge lines.

```
nodes <N> edges <M>
v <id> <x> <y>
v <id> <x> <y>
...
e <from> <to> <distance> <time>
e <from> <to> <distance> <time>
...
```

Notes:
- Coordinates are in nanodegrees (multiply by 1,000,000 to obtain degrees).
- Distances are in decimeters (0.1 m).

### Instance lists (`instances/*.txt`)

First line: number of instances. Then one instance per line:

```
<source> <target> <budget>
```

Example:

```
160
2191918 1689718 7501374
2191918 1689718 7702103
```

## Output

Each run writes to an output directory (e.g. `results/W_rep01_nit10_perc0.0/`).

Files produced:
- `results.txt`: one line per instance with costs, resources, timings, and algorithm metadata.
- `<source>-<target>-<budget>_sol_details.txt`: improvement trace for each instance.

## Script Configuration

Edit `run_parc.sh` to change:
- `maps` and `instances` lists
- `REPS` (repetitions)
- Algorithm parameters:
  - `TL` time limit (seconds)
  - `REDH` reduction heuristic flag
  - `NIT` max iterations for lambda search
  - `PERC_RED` graph reduction percentage
  - `MULTIPLE_INSTANCE_FLAG` (1 to read instance list)

## Repository Layout

- `src/parc.c`: core implementation
- `run_parc.sh`: batch runner
- `map/`: graph files
- `instances/`: instance lists
- `results/`: outputs
