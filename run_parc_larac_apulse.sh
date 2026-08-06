#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAP_DIR="${SCRIPT_DIR}/map"
SRC_DIR="${SCRIPT_DIR}/src"
INPUT_DIR="${SCRIPT_DIR}/instances"
RESULTS_DIR="${SCRIPT_DIR}/results"
SCRIPT_NAME="parc_larac_apulse.c"
EXECUTABLE_NAME="${SCRIPT_DIR}/parc_larac_apulse"

REPS=1
maps=(
#   "W.xy"
#   "CTR.txt"
   "USA.txt"
)

instances=(
#  "instances_W.txt"
#  "instances_CTR.txt"
  "instances_USA.txt"
)

# Common parameters
TL="60.0"
MULTIPLE_INSTANCE_FLAG="1"

# Algorithm selection: at least one must be 1
PARC="0"
LARAC="1"
APULSE="0"

# PARC parameters
REDH="1"
NIT="10"
PERC_RED="0.0"

# APULSE parameters
N="8192"
MIN_BUCKET_WIDTH="1"

if (( ${#maps[@]} == 0 )); then
  echo "No maps configured." >&2; exit 1
fi
if (( ${#maps[@]} != ${#instances[@]} )); then
  echo "List size mismatch: maps=${#maps[@]}, instances=${#instances[@]}" >&2; exit 1
fi
for flag_name in PARC LARAC APULSE; do
  flag_value="${!flag_name}"
  if [[ "$flag_value" != "0" && "$flag_value" != "1" ]]; then
    echo "$flag_name must be 0 or 1." >&2; exit 1
  fi
done
if [[ "$PARC" == "0" && "$LARAC" == "0" && "$APULSE" == "0" ]]; then
  echo "At least one algorithm must be enabled." >&2; exit 1
fi
if [[ "$N" -le 0 ]]; then echo "N must be greater than zero." >&2; exit 1; fi

mkdir -p "$RESULTS_DIR"

SOURCE_PATH="${SRC_DIR}/${SCRIPT_NAME}"
if [[ ! -x "$EXECUTABLE_NAME" || "$SOURCE_PATH" -nt "$EXECUTABLE_NAME" ]]; then
  echo "Compiling ${SCRIPT_NAME}..."
  gcc -std=c11 -O3 -DNDEBUG -march=native "$SOURCE_PATH" -o "$EXECUTABLE_NAME" -lm
fi

for i in "${!maps[@]}"; do
  in_file="${maps[$i]}"
  instance_file="${instances[$i]}"
  in_path="${MAP_DIR}/${in_file}"
  instance_path="${INPUT_DIR}/${instance_file}"
  [[ -f "$in_path" ]] || { echo "Missing input file: $in_path" >&2; exit 1; }
  [[ -f "$instance_path" ]] || { echo "Missing instance file: $instance_path" >&2; exit 1; }

  in_stem="${in_file%.*}"
  for ((rep=1; rep<=REPS; rep++)); do
    outdir="${RESULTS_DIR}/${in_stem}_rep$(printf '%02d' "$rep")_parc${PARC}_larac${LARAC}_apulse${APULSE}_nit${NIT}_N${N}"
    mkdir -p "$outdir"

    echo "==> input=$in_path | PARC=$PARC LARAC=$LARAC APULSE=$APULSE | rep=$rep/$REPS"
    "$EXECUTABLE_NAME" \
      --input "$in_path" \
      --tl "$TL" \
      --redh "$REDH" \
      --nit "$NIT" \
      --multipleinstanceflag "$MULTIPLE_INSTANCE_FLAG" \
      --inputinstance "$instance_path" \
      --outdir "$outdir" \
      --perc_red "$PERC_RED" \
      --parc "$PARC" \
      --larac "$LARAC" \
      --apulse "$APULSE" \
      --N "$N" \
      --min-bucket-width "$MIN_BUCKET_WIDTH"
  done
done
