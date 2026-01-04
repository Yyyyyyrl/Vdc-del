#!/bin/bash
#
# vdc_test.sh - Comprehensive testing script for vdc-del
#
# Runs vdc-del with multiple configurations on test datasets and
# generates a report with timing, manifold/orientation, self-intersection,
# and angle statistics. Supports GNU Parallel for faster execution.
#
# Usage: ./vdc_test.sh [options]
#   -j <N>    Number of parallel jobs (default: auto, requires GNU Parallel)
#   -h        Show help
#

set -o pipefail

# =============================================================================
# Configuration
# =============================================================================

VDC_EXEC="./build/vdc-del"
IJK_EXEC="tool/ijkmeshinfo"
DATA_DIR="data/volvis"
TMP_DIR=$(mktemp -d)
RESULTS_DIR="$TMP_DIR/results"
mkdir -p "$RESULTS_DIR"
ERROR_LOG="vdc_test_errors.log"

# Test datasets: name, file, isovalue 
DATASETS=("fuel" "aneurysm" "engine" "bonsai" "silicium" "neghip" "lobster" "marschnerlobb" "skull")
ISOVALUES=("70.5" "30.5" "100.5" "50.5" "120.5" "100.5" "20.5" "100.5" "45.5")

# Configurations: output filename postfix and flags
CONFIG_NAMES=("default" "sep2_dist2" "sep2_dist3" "sep3_dist4")
CONFIG_FLAGS=("" "-sep_split 2 -sep_dist 2" "-sep_split 2 -sep_dist 3" "-sep_split 3 -sep_dist 4") 

NUM_DATASETS=${#DATASETS[@]}
NUM_CONFIGS=${#CONFIG_NAMES[@]}
TOTAL_TESTS=$((NUM_DATASETS * NUM_CONFIGS))

# Parallel options
NUM_JOBS=""

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
    echo "  -j <N>    Number of parallel jobs (default: number of CPU cores)"
    echo "  -h        Show this help message"
    echo ""
    echo "Datasets tested: ${DATASETS[*]}"
    echo ""
    echo "Configurations tested:"
    for ((c=0; c<NUM_CONFIGS; c++)); do
        if [[ -z "${CONFIG_FLAGS[$c]}" ]]; then
            echo "  - ${CONFIG_NAMES[$c]}: (no flags)"
        else
            echo "  - ${CONFIG_NAMES[$c]}: ${CONFIG_FLAGS[$c]}"
        fi
    done
    exit 0
}

# =============================================================================
# Parse Arguments
# =============================================================================

while [[ $# -gt 0 ]]; do
    case $1 in
        -j)     NUM_JOBS="$2"; shift 2 ;;
        -h|--help) usage ;;
        *)      echo "Unknown option: $1"; usage ;;
    esac
done

# =============================================================================
# Worker Function (for parallel execution)
# =============================================================================

export VDC_EXEC IJK_EXEC DATA_DIR TMP_DIR RESULTS_DIR ERROR_LOG

run_single_test() {
    IFS='|' read -r d_idx c_idx dataset isovalue config_name config_flag <<< "$1"
    
    local input_file="${DATA_DIR}/${dataset}.nhdr"
    local mesh_file="$TMP_DIR/mesh_${d_idx}_${c_idx}.off"
    local result_file="$RESULTS_DIR/result_${d_idx}_${c_idx}.txt"
    
    # Run vdc-del
    local timing="N/A"
    local vdc_output
    local cmd="$VDC_EXEC -timing_stats -o $mesh_file $config_flag $isovalue $input_file"
    
    if vdc_output=$($cmd 2>&1); then
        timing=$(echo "$vdc_output" | grep "Total Processing Time" | sed 's/.*Time:[[:space:]]*\([0-9.]*\).*/\1/')
        [[ -z "$timing" ]] && timing="N/A"
    else
        echo "[$(date)] vdc-del FAILED: $dataset/$config_name" >> "$ERROR_LOG"
        echo "Command: $cmd" >> "$ERROR_LOG"
        echo "---" >> "$ERROR_LOG"
    fi
    
    if [[ ! -f "$mesh_file" ]]; then
        echo "$dataset|$config_name|$timing|N/A|N/A|N/A|N/A|N/A|N/A|N/A|N/A|N/A|N/A|N/A" > "$result_file"
        return
    fi
    
    # Check manifold
    local nm_edges=0 nm_verts=0 orient="OK"
    local tmp_out="$TMP_DIR/manifold_${d_idx}_${c_idx}.txt"
    $IJK_EXEC -oriented_manifold "$mesh_file" > "$tmp_out" 2>&1
    
    if head -100 "$tmp_out" | grep -qiE "degenerate|manifold|polygon"; then
        if head -50 "$tmp_out" | grep -qi "num non-manifold edges"; then
            nm_edges=$(head -50 "$tmp_out" | grep -i "num non-manifold edges" | sed 's/.*:[[:space:]]*\([0-9]*\).*/\1/')
            [[ -z "$nm_edges" ]] && nm_edges=0
        fi
        if head -50 "$tmp_out" | grep -qi "num non-manifold vert"; then
            nm_verts=$(head -50 "$tmp_out" | grep -i "num non-manifold vert" | sed 's/.*:[[:space:]]*\([0-9]*\).*/\1/')
            [[ -z "$nm_verts" ]] && nm_verts=0
        fi
        if tail -100 "$tmp_out" | grep -qi "orientation mismatch"; then
            orient="FAIL"
        fi
    else
        nm_edges="N/A"
        nm_verts="N/A"
        orient="N/A"
    fi
    rm -f "$tmp_out"
    
    # Check self-intersections
    local selfi=0
    local selfi_out
    selfi_out=$($IJK_EXEC -selfI "$mesh_file" 2>&1)
    if echo "$selfi_out" | grep -qi "no self intersections"; then
        selfi=0
    else
        selfi=$(echo "$selfi_out" | grep -i "pairs of intersecting triangles" | sed 's/.*:[[:space:]]*\([0-9]*\).*/\1/')
        [[ -z "$selfi" || ! "$selfi" =~ ^[0-9]+$ ]] && selfi="N/A"
    fi
    
    # Check angles
    local min_angle="N/A" max_angle="N/A"
    local angle_out
    if angle_out=$($IJK_EXEC -report_deep "$mesh_file" 2>&1); then
        min_angle=$(echo "$angle_out" | grep "Min polygon angle" | sed 's/.*:[[:space:]]*\([0-9.]*\).*/\1/')
        max_angle=$(echo "$angle_out" | grep "Max polygon angle" | sed 's/.*:[[:space:]]*\([0-9.]*\).*/\1/')
        [[ -z "$min_angle" ]] && min_angle="N/A"
        [[ -z "$max_angle" ]] && max_angle="N/A"
    fi
    
    # Get angle bucket counts
    local cnt_5=0 cnt_10=0 cnt_15=0 cnt_20=0 cnt_30=0
    for thresh in 5 10 15 20 30; do
        local cnt_out cnt
        if cnt_out=$($IJK_EXEC -angle_le "$thresh" "$mesh_file" 2>&1); then
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
    
    # Output result: dataset|config|timing|nm_e|nm_v|orient|selfi|min_a|max_a|a5|a10|a15|a20|a30
    echo "$dataset|$config_name|$timing|$nm_edges|$nm_verts|$orient|$selfi|$min_angle|$max_angle|$cnt_5|$cnt_10|$cnt_15|$cnt_20|$cnt_30" > "$result_file"
    
    rm -f "$mesh_file"
}
export -f run_single_test

# =============================================================================
# Main Execution
# =============================================================================

echo "VDC-DEL Test Suite"
echo "=================="
echo ""
echo "Executable:  $VDC_EXEC"
echo "IJK Tool:    $IJK_EXEC"
echo "Data Dir:    $DATA_DIR"
echo "Temp Dir:    $TMP_DIR"
echo "Error Log:   $ERROR_LOG"
echo ""

# Clear error log
> "$ERROR_LOG"

# Check prerequisites
if [[ ! -x "$VDC_EXEC" ]]; then
    echo "ERROR: vdc-del executable not found at $VDC_EXEC"
    exit 1
fi

if [[ ! -x "$IJK_EXEC" ]]; then
    echo "ERROR: ijkmeshinfo executable not found at $IJK_EXEC"
    exit 1
fi

if [[ ! -d "$DATA_DIR" ]]; then
    echo "ERROR: Data directory not found at $DATA_DIR"
    exit 1
fi

# Check for GNU Parallel
USE_PARALLEL=true
if ! command -v parallel &> /dev/null; then
    echo "Note: GNU Parallel not found. Running in sequential mode."
    echo "      Install with: brew install parallel"
    USE_PARALLEL=false
else
    echo "Parallel:    ${NUM_JOBS:-auto} jobs"
fi
echo ""

# =============================================================================
# Create Job List
# =============================================================================

JOBS_FILE="$TMP_DIR/jobs.txt"
> "$JOBS_FILE"

for ((d=0; d<NUM_DATASETS; d++)); do
    dataset="${DATASETS[$d]}"
    isovalue="${ISOVALUES[$d]}"
    for ((c=0; c<NUM_CONFIGS; c++)); do
        config_name="${CONFIG_NAMES[$c]}"
        config_flag="${CONFIG_FLAGS[$c]}"
        # Format: d_idx|c_idx|dataset|isovalue|config_name|config_flags
        echo "$d|$c|$dataset|$isovalue|$config_name|$config_flag"
    done
done >> "$JOBS_FILE"

# =============================================================================
# Run Tests
# =============================================================================

echo "Running $TOTAL_TESTS tests ($NUM_DATASETS datasets × $NUM_CONFIGS configs)..."

if [[ "$USE_PARALLEL" == "true" ]]; then
    echo "Using GNU Parallel..."
    echo ""
    
    PARALLEL_OPTS="--bar --eta"
    [[ -n "$NUM_JOBS" ]] && PARALLEL_OPTS+=" -j $NUM_JOBS"
    
    # GNU Parallel uses $SHELL to run jobs. On macOS this is often zsh, which
    # cannot execute exported bash functions. Force bash so `export -f` works.
    env SHELL=/bin/bash parallel $PARALLEL_OPTS run_single_test :::: "$JOBS_FILE"
else
    echo "Running sequentially..."
    echo ""
    
    TOTAL_JOBS=$(wc -l < "$JOBS_FILE")
    CURRENT_JOB=0
    
    while IFS= read -r job_line; do
        ((CURRENT_JOB++))
        IFS='|' read -r _ _ dataset _ config_name _ <<< "$job_line"
        
        # Progress bar
        local_percent=$((CURRENT_JOB * 100 / TOTAL_JOBS))
        local_filled=$((CURRENT_JOB * 30 / TOTAL_JOBS))
        bar="["
        for ((i=0; i<local_filled; i++)); do bar+="="; done
        if ((local_filled < 30)); then bar+=">"; fi
        for ((i=0; i<30-local_filled-1; i++)); do bar+=" "; done
        bar+="]"
        printf "\r%-32s %2d/%d (%3d%%) %-12s %-12s" "$bar" "$CURRENT_JOB" "$TOTAL_JOBS" "$local_percent" "$dataset" "$config_name"
        
        run_single_test "$job_line"
    done < "$JOBS_FILE"
    
    printf "\r%-80s\r" ""
fi

echo ""
echo "Tests completed! Aggregating results..."
echo ""

# =============================================================================
# Aggregate Results and Print Table
# =============================================================================

# Results storage
declare -a RESULT_DATASET RESULT_CONFIG RESULT_TIME
declare -a RESULT_NM_EDGES RESULT_NM_VERTS RESULT_ORIENT RESULT_SELFI
declare -a RESULT_MIN_ANGLE RESULT_MAX_ANGLE
declare -a RESULT_ANGLE_5 RESULT_ANGLE_10 RESULT_ANGLE_15 RESULT_ANGLE_20 RESULT_ANGLE_30

# Read results in order
idx=0
for ((d=0; d<NUM_DATASETS; d++)); do
    for ((c=0; c<NUM_CONFIGS; c++)); do
        result_file="$RESULTS_DIR/result_${d}_${c}.txt"
        if [[ -f "$result_file" ]]; then
            IFS='|' read -r dataset config timing nm_e nm_v orient selfi min_a max_a a5 a10 a15 a20 a30 < "$result_file"
            RESULT_DATASET[$idx]="$dataset"
            RESULT_CONFIG[$idx]="$config"
            RESULT_TIME[$idx]="$timing"
            RESULT_NM_EDGES[$idx]="$nm_e"
            RESULT_NM_VERTS[$idx]="$nm_v"
            RESULT_ORIENT[$idx]="$orient"
            RESULT_SELFI[$idx]="$selfi"
            RESULT_MIN_ANGLE[$idx]="$min_a"
            RESULT_MAX_ANGLE[$idx]="$max_a"
            RESULT_ANGLE_5[$idx]="$a5"
            RESULT_ANGLE_10[$idx]="$a10"
            RESULT_ANGLE_15[$idx]="$a15"
            RESULT_ANGLE_20[$idx]="$a20"
            RESULT_ANGLE_30[$idx]="$a30"
        else
            RESULT_DATASET[$idx]="${DATASETS[$d]}"
            RESULT_CONFIG[$idx]="${CONFIG_NAMES[$c]}"
            RESULT_TIME[$idx]="N/A"
            RESULT_NM_EDGES[$idx]="N/A"
            RESULT_NM_VERTS[$idx]="N/A"
            RESULT_ORIENT[$idx]="N/A"
            RESULT_SELFI[$idx]="N/A"
            RESULT_MIN_ANGLE[$idx]="N/A"
            RESULT_MAX_ANGLE[$idx]="N/A"
            RESULT_ANGLE_5[$idx]="N/A"
            RESULT_ANGLE_10[$idx]="N/A"
            RESULT_ANGLE_15[$idx]="N/A"
            RESULT_ANGLE_20[$idx]="N/A"
            RESULT_ANGLE_30[$idx]="N/A"
        fi
        ((idx++))
    done
done

# Print table
echo ""
echo "=============================================================================================================="
echo "                                       VDC-DEL TEST RESULTS"
echo "=============================================================================================================="
echo ""

# Header
printf "%-12s | %-12s | %8s | %5s | %5s | %6s | %6s | %7s | %7s | %5s | %5s | %5s | %5s | %5s\n" \
    "Dataset" "Config" "Time(s)" "NM-E" "NM-V" "Orient" "SelfI" "MinA" "MaxA" "<=5" "<=10" "<=15" "<=20" "<=30"

# Separator
printf "%s\n" "-------------|--------------|----------|-------|-------|--------|--------|---------|---------|-------|-------|-------|-------|-------"

# Data rows
for ((i=0; i<${#RESULT_DATASET[@]}; i++)); do
    printf "%-12s | %-12s | %8s | %5s | %5s | %6s | %6s | %7s | %7s | %5s | %5s | %5s | %5s | %5s\n" \
        "${RESULT_DATASET[$i]}" \
        "${RESULT_CONFIG[$i]}" \
        "${RESULT_TIME[$i]}" \
        "${RESULT_NM_EDGES[$i]}" \
        "${RESULT_NM_VERTS[$i]}" \
        "${RESULT_ORIENT[$i]}" \
        "${RESULT_SELFI[$i]}" \
        "${RESULT_MIN_ANGLE[$i]}" \
        "${RESULT_MAX_ANGLE[$i]}" \
        "${RESULT_ANGLE_5[$i]}" \
        "${RESULT_ANGLE_10[$i]}" \
        "${RESULT_ANGLE_15[$i]}" \
        "${RESULT_ANGLE_20[$i]}" \
        "${RESULT_ANGLE_30[$i]}"
done

echo ""
echo "=============================================================================================================="
echo "Legend: NM-E=Non-manifold edges, NM-V=Non-manifold vertices, SelfI=Self-intersections, MinA/MaxA=Angle range"
echo "        <=N = Number of polygons with angles <= N degrees"
echo "=============================================================================================================="

# Check for errors
if [[ -s "$ERROR_LOG" ]]; then
    echo ""
    echo "Note: Some errors occurred. See $ERROR_LOG for details."
else
    rm -f "$ERROR_LOG"
fi
