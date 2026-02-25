#!/usr/bin/env bash
set -euo pipefail

# ===================== CONFIG =====================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAP_DIR="${SCRIPT_DIR}/map"
SRC_DIR="${SCRIPT_DIR}/src"
INPUT_DIR="${SCRIPT_DIR}/instances"
RESULTS_DIR="${SCRIPT_DIR}/results"
SCRIPT_NAME="parc.c"
EXECUTABLE_NAME="parc"

# Number of repetitions for each (input, instances) pair
REPS=1

# Pair list (file names only, no paths)
maps=(
   "W.xy"
   "CTR.txt"
   "USA.txt"
)

instances=(
  "instances_W.txt"
  "instances_CTR.txt"
  "instances_USA.txt"
)

# Algorithm parameters
TL="60.0" #time limit
REDH="1"  #apply reduction
MULTIPLE_INSTANCE_FLAG="1" #1 = run multiple instance given a file with the list of instances - 0= specify source or target
NIT="10" #number of iteration for binary search
PERC_RED="0.0" #percentage of reduction applied in graph reduction in range (0.0-1.0). If = 0.0 remove only the critical node.

# ===================== CHECKS =====================

if (( ${#maps[@]} == 0 )); then
  echo "No maps configured." >&2
  exit 1
fi

if (( ${#maps[@]} != ${#instances[@]} )); then
  echo "List size mismatch: maps=${#maps[@]}, instances=${#instances[@]}" >&2
  exit 1
fi

mkdir -p "$RESULTS_DIR"

# ===================== RUNS =====================

for i in "${!maps[@]}"; do
  in_file="${maps[$i]}"
  instance_file="${instances[$i]}"

  # Build full paths
  in_path="${MAP_DIR}/${in_file}"
  instance_file_path="${INPUT_DIR}/${instance_file}"

  if [[ ! -f "$in_path" ]]; then
    echo "Missing input file: $in_path" >&2
    exit 1
  fi

  if [[ ! -f "$instance_file_path" ]]; then
    echo "Missing instance file: $instance_file_path" >&2
    exit 1
  fi

  # Base name without extension
  in_stem="${in_file%.*}"

  for ((rep=1; rep<=REPS; rep++)); do
    outdir="${RESULTS_DIR}/${in_stem}_rep$(printf "%02d" "$rep")_nit${NIT}_perc${PERC_RED}"
    mkdir -p "$outdir"

    #check if executable exist otherwise compile code
    if [[ ! -x "$EXECUTABLE_NAME" ]]; then
      echo "Executable not found. Compiling ${SCRIPT_NAME}..."
      gcc -O3 "${SRC_DIR}/${SCRIPT_NAME}" -o "$EXECUTABLE_NAME" -lm
    fi

    echo "==> Running: input=$in_path | rep=$rep/$REPS"

    ./"${EXECUTABLE_NAME}" \
      --input "$in_path" \
      --tl "$TL" \
      --redh "$REDH" \
      --nit "$NIT" \
      --multipleinstanceflag "$MULTIPLE_INSTANCE_FLAG" \
      --inputinstance "$instance_file_path" \
      --outdir "$outdir" \
      --perc_red "$PERC_RED"

  done
done
