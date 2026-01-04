#!/bin/bash
#
# vdc_random_test.sh - Random scalar field testing for vdc-del
#
# Generates random scalar fields using ijkgenscalar, picks random isovalues,
# runs vdc-del with multiple configurations using GNU Parallel, and produces
# a concise summary of manifold/angle/self-intersection stats.
#
# Usage: ./vdc_random_test.sh [options]
#   -n <N>        Number of random fields to generate (default: 5)
#   -i <N>        Number of random isovalues per field (default: 3)
#   -dim <N>      Grid dimension (default: 10)
#   -maxval <N>   Maximum scalar value (default: 100)
#   -seed <N>     Base random seed (default: random)
#   -j <N>        Number of parallel jobs (default: auto)
#   -h            Show help
#

set -o pipefail

# =============================================================================
# Configuration
# =============================================================================

VDC_EXEC="./build/vdc-del"
IJK_GEN="./tool/ijkgenscalar"
IJK_MESH="./tool/ijkmeshinfo"
TMP_DIR=$(mktemp -d)
RESULTS_DIR="$TMP_DIR/results"
NRRD_DIR="$TMP_DIR/nrrd"
mkdir -p "$RESULTS_DIR" "$NRRD_DIR"

# Defaults
NUM_FIELDS=5
NUM_ISOVALUES=3
GRID_DIM=10
MAX_VAL=100
BASE_SEED=""
NUM_JOBS=""

# Configurations: name and flags (same as vdc_test.sh)
CONFIG_NAMES=("default" "sep2_dist2" "sep2_dist3" "sep3_dist4" "sep4_dist5")
CONFIG_FLAGS=("" "-sep_split 2 -sep_dist 2" "-sep_split 2 -sep_dist 3" "-sep_split 3 -sep_dist 4" "-sep_split 4 -sep_dist 5") 

NUM_CONFIGS=${#CONFIG_NAMES[@]}

# =============================================================================
# Cleanup
# =============================================================================

cleanup() {
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

# =============================================================================
# Helper Functions
# =============================================================================

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -n <N>        Number of random fields to generate (default: $NUM_FIELDS)"
    echo "  -i <N>        Number of random isovalues per field (default: $NUM_ISOVALUES)"
    echo "  -dim <N>      Grid dimension NxNxN (default: $GRID_DIM)"
    echo "  -maxval <N>   Maximum scalar value (default: $MAX_VAL)"
    echo "  -seed <N>     Base random seed (default: random)"
    echo "  -j <N>        Number of parallel jobs (default: number of CPU cores)"
    echo "  -h            Show this help message"
    echo ""
    echo "Configurations tested:"
    for ((c=0; c<NUM_CONFIGS; c++)); do
        if [[ -z "${CONFIG_FLAGS[$c]}" ]]; then
            echo "  - ${CONFIG_NAMES[$c]}: (no flags)"
        else
            echo "  - ${CONFIG_NAMES[$c]}: ${CONFIG_FLAGS[$c]}"
        fi
    done
    echo ""
    echo "Example:"
    echo "  $0 -n 10 -i 5 -dim 15 -maxval 200 -seed 12345 -j 8"
    exit 0
}

# =============================================================================
# Parse Arguments
# =============================================================================

while [[ $# -gt 0 ]]; do
    case $1 in
        -n)     NUM_FIELDS="$2"; shift 2 ;;
        -i)     NUM_ISOVALUES="$2"; shift 2 ;;
        -dim)   GRID_DIM="$2"; shift 2 ;;
        -maxval) MAX_VAL="$2"; shift 2 ;;
        -seed)  BASE_SEED="$2"; shift 2 ;;
        -j)     NUM_JOBS="$2"; shift 2 ;;
        -h|--help) usage ;;
        *)      echo "Unknown option: $1"; usage ;;
    esac
done

# Set random seed if not specified
if [[ -z "$BASE_SEED" ]]; then
    BASE_SEED=$RANDOM
fi

# =============================================================================
# Check Prerequisites
# =============================================================================

echo "=============================================="
echo "     VDC-DEL Random Field Test Suite"
echo "=============================================="
echo ""

if [[ ! -x "$VDC_EXEC" ]]; then
    echo "ERROR: vdc-del executable not found at $VDC_EXEC"
    exit 1
fi

if [[ ! -x "$IJK_GEN" ]]; then
    echo "ERROR: ijkgenscalar not found at $IJK_GEN"
    exit 1
fi

if [[ ! -x "$IJK_MESH" ]]; then
    echo "ERROR: ijkmeshinfo not found at $IJK_MESH"
    exit 1
fi

# Check for GNU Parallel
USE_PARALLEL=true
if ! command -v parallel &> /dev/null; then
    echo "Note: GNU Parallel not found. Running in sequential mode."
    echo "      Install with: brew install parallel"
    USE_PARALLEL=false
fi

TOTAL_EXPECTED=$((NUM_FIELDS * NUM_ISOVALUES * NUM_CONFIGS))

echo "Configuration:"
echo "  Fields:       $NUM_FIELDS"
echo "  Isovalues:    $NUM_ISOVALUES per field"
echo "  Configs:      $NUM_CONFIGS per isovalue"
echo "  Grid Size:    ${GRID_DIM}x${GRID_DIM}x${GRID_DIM}"
echo "  Value Range:  0-$MAX_VAL"
echo "  Base Seed:    $BASE_SEED"
echo "  Parallel:     ${NUM_JOBS:-auto}"
echo ""

# =============================================================================
# Generate All NRRD Files First (Sequential - Fast)
# =============================================================================

echo "Generating $NUM_FIELDS random scalar fields..."
for ((f=0; f<NUM_FIELDS; f++)); do
    FIELD_SEED=$((BASE_SEED + f * 1000))
    NRRD_FILE="$NRRD_DIR/field_${f}_seed_${FIELD_SEED}.nrrd"
    $IJK_GEN -field randomint -dim 3 -asize "$GRID_DIM" -maxval "$MAX_VAL" -seed "$FIELD_SEED" "$NRRD_FILE" > /dev/null 2>&1
done
echo "Done generating fields."
echo ""

# =============================================================================
# Create Test Job List
# =============================================================================

JOBS_FILE="$TMP_DIR/jobs.txt"
> "$JOBS_FILE"

for ((f=0; f<NUM_FIELDS; f++)); do
    FIELD_SEED=$((BASE_SEED + f * 1000))
    NRRD_FILE="$NRRD_DIR/field_${f}_seed_${FIELD_SEED}.nrrd"
    
    for ((iso_i=0; iso_i<NUM_ISOVALUES; iso_i++)); do
        # Generate deterministic isovalue
        RANDOM=$((FIELD_SEED + iso_i))
        ISO_RANGE=$((MAX_VAL - 10))
        ISOVALUE="$((RANDOM % ISO_RANGE + 5)).5"
        
        for ((c=0; c<NUM_CONFIGS; c++)); do
            CONFIG_NAME="${CONFIG_NAMES[$c]}"
            CONFIG_FLAG="${CONFIG_FLAGS[$c]}"
            # Format: field_idx|iso_idx|field_seed|isovalue|config_idx|config_name|config_flags|nrrd_file
            echo "$f|$iso_i|$FIELD_SEED|$ISOVALUE|$c|$CONFIG_NAME|$CONFIG_FLAG|$NRRD_FILE"
        done
    done
done >> "$JOBS_FILE"

# =============================================================================
# Worker Function (exported for parallel)
# =============================================================================

export VDC_EXEC IJK_MESH RESULTS_DIR TMP_DIR

run_single_test() {
    IFS='|' read -r f iso_i field_seed isovalue c config_name config_flag nrrd_file <<< "$1"
    
    local mesh_file="$TMP_DIR/mesh_${f}_${iso_i}_${c}.off"
    local result_file="$RESULTS_DIR/result_${f}_${iso_i}_${c}.txt"
    
    # Run vdc-del
    local vdc_output
    if vdc_output=$($VDC_EXEC -timing_stats -o "$mesh_file" $config_flag "$isovalue" "$nrrd_file" 2>&1); then
        local verts tris
        verts=$(echo "$vdc_output" | grep "Generated.*vertices" | sed 's/.*Generated \([0-9]*\) vertices.*/\1/')
        tris=$(echo "$vdc_output" | grep "Generated.*triangles" | sed 's/.*Generated \([0-9]*\) triangles.*/\1/')
        [[ -z "$verts" ]] && verts=0
        [[ -z "$tris" ]] && tris=0
    else
        echo "FAIL|$c|$f|$field_seed|$isovalue|$config_name|$config_flag|EXEC_FAIL|0|0|0|0|OK|0|N/A|N/A|0|0|0|0|0" > "$result_file"
        rm -f "$mesh_file"
        return
    fi
    
    if [[ ! -f "$mesh_file" ]]; then
        echo "FAIL|$c|$f|$field_seed|$isovalue|$config_name|$config_flag|NO_MESH|0|0|0|0|OK|0|N/A|N/A|0|0|0|0|0" > "$result_file"
        return
    fi
    
    # Check manifold
    local nm_edges=0 nm_verts=0 orient="OK"
    local manifold_out
    manifold_out=$($IJK_MESH -oriented_manifold "$mesh_file" 2>&1)
    
    if echo "$manifold_out" | head -50 | grep -qi "num non-manifold edges"; then
        nm_edges=$(echo "$manifold_out" | head -50 | grep -i "num non-manifold edges" | sed 's/.*:[[:space:]]*\([0-9]*\).*/\1/')
        [[ -z "$nm_edges" ]] && nm_edges=0
    fi
    if echo "$manifold_out" | head -50 | grep -qi "num non-manifold vert"; then
        nm_verts=$(echo "$manifold_out" | head -50 | grep -i "num non-manifold vert" | sed 's/.*:[[:space:]]*\([0-9]*\).*/\1/')
        [[ -z "$nm_verts" ]] && nm_verts=0
    fi
    if echo "$manifold_out" | tail -100 | grep -qi "orientation mismatch"; then
        orient="FAIL"
    fi
    
    # Check self-intersections
    local selfi=0
    local selfi_out
    selfi_out=$($IJK_MESH -selfI "$mesh_file" 2>&1)
    selfi=$(echo "$selfi_out" | grep -i "pairs of intersecting triangles" | sed 's/.*:[[:space:]]*\([0-9]*\).*/\1/')
    [[ -z "$selfi" || ! "$selfi" =~ ^[0-9]+$ ]] && selfi=0
    
    # Check angles
    local min_a="N/A" max_a="N/A"
    local angle_out
    if angle_out=$($IJK_MESH -report_deep "$mesh_file" 2>&1); then
        min_a=$(echo "$angle_out" | grep "Min polygon angle" | sed 's/.*:[[:space:]]*\([0-9.]*\).*/\1/')
        max_a=$(echo "$angle_out" | grep "Max polygon angle" | sed 's/.*:[[:space:]]*\([0-9.]*\).*/\1/')
        [[ -z "$min_a" ]] && min_a="N/A"
        [[ -z "$max_a" ]] && max_a="N/A"
    fi
    
    # Get angle bucket counts
    local cnt_5=0 cnt_10=0 cnt_15=0 cnt_20=0 cnt_30=0
    for thresh in 5 10 15 20 30; do
        local cnt_out cnt
        if cnt_out=$($IJK_MESH -angle_le "$thresh" "$mesh_file" 2>&1); then
            cnt=$(echo "$cnt_out" | grep -i "polygons with angles" | sed 's/.*:[[:space:]]*\([0-9]*\).*/\1/')
            [[ -z "$cnt" ]] && cnt=0
            case $thresh in
                5)  cnt_5=$cnt ;;
                10) cnt_10=$cnt ;;
                15) cnt_15=$cnt ;;
                20) cnt_20=$cnt ;;
                30) cnt_30=$cnt ;;
            esac
        fi
    done
    
    # Determine pass/fail
    local status="PASS"
    local fail_reason=""
    if [[ "$nm_edges" -gt 0 || "$nm_verts" -gt 0 || "$orient" == "FAIL" || "$selfi" -gt 0 ]]; then
        status="FAIL"
        [[ "$nm_edges" -gt 0 ]] && fail_reason+="nm_edges=$nm_edges "
        [[ "$nm_verts" -gt 0 ]] && fail_reason+="nm_verts=$nm_verts "
        [[ "$orient" == "FAIL" ]] && fail_reason+="orient_fail "
        [[ "$selfi" -gt 0 ]] && fail_reason+="selfI=$selfi "
    fi
    
    # Output: status|config_idx|field|seed|iso|config_name|config_flag|fail_reason|verts|tris|nm_e|nm_v|orient|selfi|min_a|max_a|a5|a10|a15|a20|a30
    echo "$status|$c|$f|$field_seed|$isovalue|$config_name|$config_flag|$fail_reason|$verts|$tris|$nm_edges|$nm_verts|$orient|$selfi|$min_a|$max_a|$cnt_5|$cnt_10|$cnt_15|$cnt_20|$cnt_30" > "$result_file"
    
    rm -f "$mesh_file"
}
export -f run_single_test

# =============================================================================
# Run Tests
# =============================================================================

echo "Running $TOTAL_EXPECTED tests ($NUM_FIELDS fields × $NUM_ISOVALUES isovalues × $NUM_CONFIGS configs)..."

if [[ "$USE_PARALLEL" == "true" ]]; then
    echo "Using GNU Parallel..."
    echo ""
    
    PARALLEL_OPTS="--bar --eta"
    [[ -n "$NUM_JOBS" ]] && PARALLEL_OPTS+=" -j $NUM_JOBS"
    
    parallel $PARALLEL_OPTS run_single_test :::: "$JOBS_FILE"
else
    echo "Running sequentially..."
    echo ""
    
    TOTAL_JOBS=$(wc -l < "$JOBS_FILE")
    CURRENT_JOB=0
    
    while IFS= read -r job_line; do
        ((CURRENT_JOB++))
        # Parse job for progress display
        IFS='|' read -r f field_seed isovalue c config_name _ _ <<< "$job_line"
        printf "\r[%3d%%] f%s iso=%s cfg=%s          " "$((CURRENT_JOB * 100 / TOTAL_JOBS))" "$f" "$isovalue" "$config_name"
        run_single_test "$job_line"
    done < "$JOBS_FILE"
    
    printf "\r%-60s\r" ""
fi

echo ""
echo "All tests completed. Aggregating results..."

# =============================================================================
# Aggregate Results
# =============================================================================

# Per-configuration counters
declare -a CFG_TESTS CFG_PASS CFG_NM_EDGES CFG_NM_VERTS CFG_ORIENT_FAIL CFG_SELFI
declare -a CFG_VERTICES CFG_TRIANGLES CFG_MIN_ANGLE CFG_MAX_ANGLE
declare -a CFG_ANGLE_LE5 CFG_ANGLE_LE10 CFG_ANGLE_LE15 CFG_ANGLE_LE20 CFG_ANGLE_LE30
declare -a FAILED_TESTS

for ((c=0; c<NUM_CONFIGS; c++)); do
    CFG_TESTS[$c]=0; CFG_PASS[$c]=0; CFG_NM_EDGES[$c]=0; CFG_NM_VERTS[$c]=0
    CFG_ORIENT_FAIL[$c]=0; CFG_SELFI[$c]=0; CFG_VERTICES[$c]=0; CFG_TRIANGLES[$c]=0
    CFG_MIN_ANGLE[$c]=180.0; CFG_MAX_ANGLE[$c]=0.0
    CFG_ANGLE_LE5[$c]=0; CFG_ANGLE_LE10[$c]=0; CFG_ANGLE_LE15[$c]=0
    CFG_ANGLE_LE20[$c]=0; CFG_ANGLE_LE30[$c]=0
done

# Compare floats
float_lt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a < b) }'; }
float_gt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a > b) }'; }

# Process all result files
for result_file in "$RESULTS_DIR"/result_*.txt; do
    [[ -f "$result_file" ]] || continue
    
    IFS='|' read -r status c f field_seed isovalue config_name config_flag fail_reason \
        verts tris nm_e nm_v orient selfi min_a max_a a5 a10 a15 a20 a30 < "$result_file"
    
    CFG_TESTS[$c]=$((${CFG_TESTS[$c]} + 1))
    CFG_VERTICES[$c]=$((${CFG_VERTICES[$c]} + verts))
    CFG_TRIANGLES[$c]=$((${CFG_TRIANGLES[$c]} + tris))
    CFG_NM_EDGES[$c]=$((${CFG_NM_EDGES[$c]} + nm_e))
    CFG_NM_VERTS[$c]=$((${CFG_NM_VERTS[$c]} + nm_v))
    [[ "$orient" == "FAIL" ]] && CFG_ORIENT_FAIL[$c]=$((${CFG_ORIENT_FAIL[$c]} + 1))
    CFG_SELFI[$c]=$((${CFG_SELFI[$c]} + selfi))
    
    if [[ "$min_a" != "N/A" ]] && float_lt "$min_a" "${CFG_MIN_ANGLE[$c]}"; then
        CFG_MIN_ANGLE[$c]="$min_a"
    fi
    if [[ "$max_a" != "N/A" ]] && float_gt "$max_a" "${CFG_MAX_ANGLE[$c]}"; then
        CFG_MAX_ANGLE[$c]="$max_a"
    fi
    
    CFG_ANGLE_LE5[$c]=$((${CFG_ANGLE_LE5[$c]} + a5))
    CFG_ANGLE_LE10[$c]=$((${CFG_ANGLE_LE10[$c]} + a10))
    CFG_ANGLE_LE15[$c]=$((${CFG_ANGLE_LE15[$c]} + a15))
    CFG_ANGLE_LE20[$c]=$((${CFG_ANGLE_LE20[$c]} + a20))
    CFG_ANGLE_LE30[$c]=$((${CFG_ANGLE_LE30[$c]} + a30))
    
    if [[ "$status" == "PASS" ]]; then
        CFG_PASS[$c]=$((${CFG_PASS[$c]} + 1))
    else
        FAILED_TESTS+=("field=$f seed=$field_seed iso=$isovalue cfg=$config_name flags=\"$config_flag\" => $fail_reason")
    fi
done

# =============================================================================
# Print Summary
# =============================================================================

echo ""
echo "================================================================================"
echo "                              SUMMARY BY CONFIGURATION"
echo "================================================================================"
echo ""

# Calculate totals
GRAND_TESTS=0; GRAND_PASS=0; GRAND_NM_E=0; GRAND_NM_V=0; GRAND_ORIENT=0; GRAND_SELFI=0
GRAND_VERTS=0; GRAND_TRIS=0; GLOBAL_MIN_ANGLE=180.0; GLOBAL_MAX_ANGLE=0.0

for ((c=0; c<NUM_CONFIGS; c++)); do
    GRAND_TESTS=$((GRAND_TESTS + ${CFG_TESTS[$c]}))
    GRAND_PASS=$((GRAND_PASS + ${CFG_PASS[$c]}))
    GRAND_NM_E=$((GRAND_NM_E + ${CFG_NM_EDGES[$c]}))
    GRAND_NM_V=$((GRAND_NM_V + ${CFG_NM_VERTS[$c]}))
    GRAND_ORIENT=$((GRAND_ORIENT + ${CFG_ORIENT_FAIL[$c]}))
    GRAND_SELFI=$((GRAND_SELFI + ${CFG_SELFI[$c]}))
    GRAND_VERTS=$((GRAND_VERTS + ${CFG_VERTICES[$c]}))
    GRAND_TRIS=$((GRAND_TRIS + ${CFG_TRIANGLES[$c]}))
    if float_lt "${CFG_MIN_ANGLE[$c]}" "$GLOBAL_MIN_ANGLE"; then
        GLOBAL_MIN_ANGLE="${CFG_MIN_ANGLE[$c]}"
    fi
    if float_gt "${CFG_MAX_ANGLE[$c]}" "$GLOBAL_MAX_ANGLE"; then
        GLOBAL_MAX_ANGLE="${CFG_MAX_ANGLE[$c]}"
    fi
done

# Print per-config summary table
printf "%-14s | %6s | %6s | %5s | %5s | %5s | %6s | %7s | %7s\n" \
    "Config" "Tests" "Pass%" "NM-E" "NM-V" "OrErr" "SelfI" "MinA" "MaxA"
printf "%s\n" "---------------|--------|--------|-------|-------|-------|--------|---------|--------"

for ((c=0; c<NUM_CONFIGS; c++)); do
    pass_pct=$(awk -v p="${CFG_PASS[$c]}" -v t="${CFG_TESTS[$c]}" 'BEGIN { if (t > 0) printf "%.1f", p * 100.0 / t; else print "0.0" }')
    printf "%-14s | %6d | %5s%% | %5d | %5d | %5d | %6d | %7.2f | %7.2f\n" \
        "${CONFIG_NAMES[$c]}" "${CFG_TESTS[$c]}" "$pass_pct" \
        "${CFG_NM_EDGES[$c]}" "${CFG_NM_VERTS[$c]}" "${CFG_ORIENT_FAIL[$c]}" "${CFG_SELFI[$c]}" \
        "${CFG_MIN_ANGLE[$c]}" "${CFG_MAX_ANGLE[$c]}"
done

printf "%s\n" "---------------|--------|--------|-------|-------|-------|--------|---------|--------"
GRAND_PCT=$(awk -v p="$GRAND_PASS" -v t="$GRAND_TESTS" 'BEGIN { if (t > 0) printf "%.1f", p * 100.0 / t; else print "0.0" }')
printf "%-14s | %6d | %5s%% | %5d | %5d | %5d | %6d | %7.2f | %7.2f\n" \
    "TOTAL" "$GRAND_TESTS" "$GRAND_PCT" "$GRAND_NM_E" "$GRAND_NM_V" "$GRAND_ORIENT" "$GRAND_SELFI" "$GLOBAL_MIN_ANGLE" "$GLOBAL_MAX_ANGLE"

echo ""
echo "Mesh Totals: $GRAND_VERTS vertices, $GRAND_TRIS triangles"
echo ""

# Angle bucket summary
echo "Angle Distribution (all configs combined):"
TOTAL_LE5=0; TOTAL_LE10=0; TOTAL_LE15=0; TOTAL_LE20=0; TOTAL_LE30=0
for ((c=0; c<NUM_CONFIGS; c++)); do
    TOTAL_LE5=$((TOTAL_LE5 + ${CFG_ANGLE_LE5[$c]}))
    TOTAL_LE10=$((TOTAL_LE10 + ${CFG_ANGLE_LE10[$c]}))
    TOTAL_LE15=$((TOTAL_LE15 + ${CFG_ANGLE_LE15[$c]}))
    TOTAL_LE20=$((TOTAL_LE20 + ${CFG_ANGLE_LE20[$c]}))
    TOTAL_LE30=$((TOTAL_LE30 + ${CFG_ANGLE_LE30[$c]}))
done
printf "  <=5°: %d | <=10°: %d | <=15°: %d | <=20°: %d | <=30°: %d\n" \
    "$TOTAL_LE5" "$TOTAL_LE10" "$TOTAL_LE15" "$TOTAL_LE20" "$TOTAL_LE30"
echo ""

# Print failed tests
if [[ ${#FAILED_TESTS[@]} -gt 0 ]]; then
    echo "Failed Tests (first 10 of ${#FAILED_TESTS[@]}):"
    count=0
    for fail in "${FAILED_TESTS[@]}"; do
        ((count++))
        [[ $count -gt 10 ]] && break
        echo "  • $fail"
    done
    [[ ${#FAILED_TESTS[@]} -gt 10 ]] && echo "  ... and $((${#FAILED_TESTS[@]} - 10)) more"
    echo ""
fi

echo "================================================================================"
if [[ $GRAND_PASS -eq $GRAND_TESTS && $GRAND_TESTS -gt 0 ]]; then
    echo "  ✓ ALL $GRAND_TESTS TESTS PASSED"
else
    echo "  ✗ $((GRAND_TESTS - GRAND_PASS)) / $GRAND_TESTS TESTS FAILED"
fi
echo "================================================================================"
