#!/bin/bash
#
# exhaustive_test.sh - Exhaustive isovalue sweep testing for vdc-del
#
# Runs vdc-del with multiple configurations across isovalue ranges and
# collects comprehensive mesh quality statistics. Outputs results to CSV files.
#
# Usage: ./exhaustive_test.sh [options]
#   -d <dataset>    Test only this dataset
#   -c <config>     Test only this configuration
#   -j <N>          Number of parallel jobs (requires GNU Parallel)
#   -o <dir>        Output directory (default: ./exhaustive_results)
#   -resume         Resume from previous run, skip completed tests
#   -test           Use minimal isovalue ranges for quick testing (3 values each)
#   -max <N>        Maximum isovalues to test per dataset-config pair
#   -h              Show help
#
# Examples:
#   ./exhaustive_test.sh                          # Run all datasets/configs
#   ./exhaustive_test.sh -d aneurysm -j 8         # Aneurysm only, 8 parallel jobs
#   ./exhaustive_test.sh -resume                  # Resume interrupted run
#

set -o pipefail

# =============================================================================
# Configuration
# =============================================================================

VDC_EXEC="./build/vdc-del"
IJK_EXEC="tool/ijkmeshinfo"
DATA_DIR="data/volvis"
OUTPUT_DIR="./exhaustive_results"
TMP_DIR=""
ERROR_LOG=""
RESUME_MODE=false

# Universal isovalue increment for all datasets
ISOVALUE_INCREMENT="0.1"

# Dataset ranges: "dataset:start_isovalue:end_isovalue"
# Modify these ranges based on your needs (all use ISOVALUE_INCREMENT)
DATASET_RANGES=(
    "aneurysm:10:200"
    "engine:50:150"
    "fuel:30:120"
    "bonsai:20:100"
    "silicium:80:160"
    "neghip:50:150"
    "lobster:10:80"
    "marschnerlobb:50:150"
    "skull:20:100"
)

# Configurations: name and flags
CONFIG_NAMES=("default" "sep2_dist2" "sep2_dist3" "sep3_dist4")
CONFIG_FLAGS=("" "-sep_split 2 -sep_dist 2" "-sep_split 2 -sep_dist 3" "-sep_split 3 -sep_dist 4")

# Filtering options
FILTER_DATASET=""
FILTER_CONFIG=""
NUM_JOBS=""
TEST_MODE=false
MAX_ISOVALUES=0

# =============================================================================
# Helper Functions
# =============================================================================

usage() {
    cat << EOF
Exhaustive VDC-DEL Testing Script

Usage: $0 [options]

Options:
  -d <dataset>    Test only this dataset (e.g., aneurysm)
  -c <config>     Test only this configuration (e.g., default, sep2_dist2)
  -j <N>          Number of parallel jobs (default: auto, requires GNU Parallel)
  -o <dir>        Output directory (default: ./exhaustive_results)
  -resume         Resume from previous run, skip completed tests
  -test           Quick test mode (only 3 isovalues per dataset)
  -max <N>        Maximum isovalues to test per dataset-config pair
  -h              Show this help message

Datasets and Isovalue Ranges (step=$ISOVALUE_INCREMENT):
EOF
    for range in "${DATASET_RANGES[@]}"; do
        IFS=':' read -r ds start end <<< "$range"
        count=$(echo "scale=0; ($end - $start) / $ISOVALUE_INCREMENT + 1" | bc)
        printf "  %-15s %s to %s = %s values\n" "$ds:" "$start" "$end" "$count"
    done
    
    echo ""
    echo "Configurations:"
    for i in "${!CONFIG_NAMES[@]}"; do
        if [[ -z "${CONFIG_FLAGS[$i]}" ]]; then
            printf "  %-15s (no flags)\n" "${CONFIG_NAMES[$i]}"
        else
            printf "  %-15s %s\n" "${CONFIG_NAMES[$i]}" "${CONFIG_FLAGS[$i]}"
        fi
    done
    
    echo ""
    echo "Output:"
    echo "  CSV files in <output_dir>/csv/<dataset>_<config>.csv"
    echo "  Summary reports in <output_dir>/summary/"
    exit 0
}

cleanup() {
    if [[ -n "$TMP_DIR" && -d "$TMP_DIR" ]]; then
        rm -rf "$TMP_DIR"
    fi
}

get_dataset_range() {
    local target="$1"
    for range in "${DATASET_RANGES[@]}"; do
        IFS=':' read -r ds start end <<< "$range"
        if [[ "$ds" == "$target" ]]; then
            echo "$start:$end:$ISOVALUE_INCREMENT"
            return 0
        fi
    done
    return 1
}

# Float comparison using bc
float_le() {
    local result
    result=$(echo "$1 <= $2" | bc -l)
    [[ "$result" -eq 1 ]]
}

# Generate list of isovalues using bc for precision
generate_isovalues() {
    local start="$1" end="$2" inc="$3"
    local val="$start"
    while float_le "$val" "$end"; do
        echo "$val"
        val=$(echo "scale=4; $val + $inc" | bc)
    done
}

# =============================================================================
# Parse Arguments
# =============================================================================

while [[ $# -gt 0 ]]; do
    case $1 in
        -d) FILTER_DATASET="$2"; shift 2 ;;
        -c) FILTER_CONFIG="$2"; shift 2 ;;
        -j) NUM_JOBS="$2"; shift 2 ;;
        -o) OUTPUT_DIR="$2"; shift 2 ;;
        -resume) RESUME_MODE=true; shift ;;
        -test) TEST_MODE=true; shift ;;
        -max) MAX_ISOVALUES="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# =============================================================================
# Setup
# =============================================================================

TMP_DIR=$(mktemp -d)
trap cleanup EXIT

mkdir -p "$OUTPUT_DIR/csv"
mkdir -p "$OUTPUT_DIR/summary"
ERROR_LOG="$OUTPUT_DIR/exhaustive_test_errors.log"

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
fi

# Check for bc
if ! command -v bc &> /dev/null; then
    echo "ERROR: 'bc' is required for floating-point calculations."
    echo "       Install with: brew install bc"
    exit 1
fi

# =============================================================================
# Export Functions and Variables for Parallel
# =============================================================================

export VDC_EXEC IJK_EXEC DATA_DIR TMP_DIR OUTPUT_DIR ERROR_LOG

# Worker function for a single test
run_single_test() {
    IFS='|' read -r dataset isovalue config_name config_flags <<< "$1"
    
    local input_file="${DATA_DIR}/${dataset}.nhdr"
    local mesh_file="$TMP_DIR/mesh_${dataset}_${config_name}_${isovalue}.off"
    
    # Run vdc-del
    local timing="N/A"
    local vertices="N/A"
    local triangles="N/A"
    local vdc_output
    local cmd="$VDC_EXEC -timing_stats -o $mesh_file $config_flags $isovalue $input_file"
    
    if vdc_output=$($cmd 2>&1); then
        timing=$(echo "$vdc_output" | grep "Total Processing Time" | sed 's/.*Time:[[:space:]]*\([0-9.]*\).*/\1/')
        [[ -z "$timing" ]] && timing="N/A"
    else
        echo "[$(date)] vdc-del FAILED: $dataset/$config_name iso=$isovalue" >> "$ERROR_LOG"
        echo "Command: $cmd" >> "$ERROR_LOG"
        echo "---" >> "$ERROR_LOG"
    fi
    
    if [[ ! -f "$mesh_file" ]]; then
        echo "$dataset,$config_name,$isovalue,$timing,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A"
        return
    fi
    
    # Get mesh size from file (OFF format: first data line has vertices faces edges)
    local mesh_header
    mesh_header=$(head -2 "$mesh_file" | tail -1)
    vertices=$(echo "$mesh_header" | awk '{print $1}')
    triangles=$(echo "$mesh_header" | awk '{print $2}')
    [[ -z "$vertices" ]] && vertices="N/A"
    [[ -z "$triangles" ]] && triangles="N/A"
    
    # Check manifold
    local nm_edges=0 nm_verts=0 orient="OK"
    local tmp_out="$TMP_DIR/manifold_${dataset}_${config_name}_${isovalue}.txt"
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
    
    # Output result line
    echo "$dataset,$config_name,$isovalue,$timing,$vertices,$triangles,$nm_edges,$nm_verts,$orient,$selfi,$min_angle,$max_angle,$cnt_5,$cnt_10,$cnt_15,$cnt_20,$cnt_30"
    
    rm -f "$mesh_file"
}
export -f run_single_test

# =============================================================================
# Main Execution
# =============================================================================

echo "======================================================================"
echo "                  VDC-DEL Exhaustive Test Suite"
echo "======================================================================"
echo ""
echo "Executable:    $VDC_EXEC"
echo "IJK Tool:      $IJK_EXEC"
echo "Data Dir:      $DATA_DIR"
echo "Output Dir:    $OUTPUT_DIR"
echo "Temp Dir:      $TMP_DIR"
echo "Error Log:     $ERROR_LOG"
echo "Resume Mode:   $RESUME_MODE"
if [[ "$USE_PARALLEL" == "true" ]]; then
    echo "Parallel:      ${NUM_JOBS:-auto} jobs"
else
    echo "Parallel:      disabled (sequential)"
fi
echo ""

# Clear error log (unless resuming)
if [[ "$RESUME_MODE" == "false" ]]; then
    > "$ERROR_LOG"
fi

# Build list of datasets to test
DATASETS_TO_TEST=()
for range in "${DATASET_RANGES[@]}"; do
    IFS=':' read -r ds _ _ <<< "$range"
    if [[ -z "$FILTER_DATASET" || "$FILTER_DATASET" == "$ds" ]]; then
        if [[ -f "${DATA_DIR}/${ds}.nhdr" ]]; then
            DATASETS_TO_TEST+=("$ds")
        else
            echo "Warning: Dataset file ${DATA_DIR}/${ds}.nhdr not found, skipping."
        fi
    fi
done

if [[ ${#DATASETS_TO_TEST[@]} -eq 0 ]]; then
    echo "ERROR: No datasets to test. Check -d option or DATASET_RANGES."
    exit 1
fi

# Build list of configs to test
CONFIGS_TO_TEST=()
CONFIG_FLAGS_TO_TEST=()
for i in "${!CONFIG_NAMES[@]}"; do
    if [[ -z "$FILTER_CONFIG" || "$FILTER_CONFIG" == "${CONFIG_NAMES[$i]}" ]]; then
        CONFIGS_TO_TEST+=("${CONFIG_NAMES[$i]}")
        CONFIG_FLAGS_TO_TEST+=("${CONFIG_FLAGS[$i]}")
    fi
done

if [[ ${#CONFIGS_TO_TEST[@]} -eq 0 ]]; then
    echo "ERROR: No configurations to test. Check -c option."
    exit 1
fi

echo "Datasets:      ${DATASETS_TO_TEST[*]}"
echo "Configs:       ${CONFIGS_TO_TEST[*]}"
echo ""

# =============================================================================
# Run Tests for Each Dataset-Config Pair
# =============================================================================

CSV_HEADER="dataset,config,isovalue,time_s,vertices,triangles,nm_edges,nm_verts,orient,self_intersections,min_angle,max_angle,angle_le_5,angle_le_10,angle_le_15,angle_le_20,angle_le_30"

for dataset in "${DATASETS_TO_TEST[@]}"; do
    range_str=$(get_dataset_range "$dataset")
    if [[ -z "$range_str" ]]; then
        echo "Warning: No range defined for $dataset, skipping."
        continue
    fi
    
    IFS=':' read -r start end inc <<< "$range_str"
    
    # In test mode, use only 3 isovalues: start, middle, end
    if [[ "$TEST_MODE" == "true" ]]; then
        mid=$(echo "scale=4; ($start + $end) / 2" | bc)
        test_isovalues=("$start" "$mid" "$end")
        iso_count=3
    else
        iso_count=$(echo "scale=0; ($end - $start) / $inc + 1" | bc)
    fi
    
    for ci in "${!CONFIGS_TO_TEST[@]}"; do
        config="${CONFIGS_TO_TEST[$ci]}"
        config_flags="${CONFIG_FLAGS_TO_TEST[$ci]}"
        
        csv_file="$OUTPUT_DIR/csv/${dataset}_${config}.csv"
        
        echo "----------------------------------------------------------------------"
        echo "Testing: $dataset + $config"
        if [[ "$TEST_MODE" == "true" ]]; then
            echo "Isovalues: TEST MODE - ${test_isovalues[*]}"
        else
            echo "Isovalues: $start to $end (step $inc) = $iso_count values"
        fi
        echo "Output: $csv_file"
        echo "----------------------------------------------------------------------"
        
        # Collect existing isovalues if resuming
        declare -A EXISTING_ISOVALS
        if [[ "$RESUME_MODE" == "true" && -f "$csv_file" ]]; then
            while IFS=',' read -r _ _ iso _; do
                EXISTING_ISOVALS["$iso"]=1
            done < <(tail -n +2 "$csv_file")
            echo "Resume: Found ${#EXISTING_ISOVALS[@]} existing results"
        fi
        
        # Write header if new file
        if [[ ! -f "$csv_file" || "$RESUME_MODE" == "false" ]]; then
            echo "$CSV_HEADER" > "$csv_file"
        fi
        
        # Generate job list
        JOBS_FILE="$TMP_DIR/jobs_${dataset}_${config}.txt"
        > "$JOBS_FILE"
        
        jobs_count=0
        if [[ "$TEST_MODE" == "true" ]]; then
            # Test mode: use pre-computed isovalues
            for isovalue in "${test_isovalues[@]}"; do
                # Skip if already computed (resume mode)
                if [[ "$RESUME_MODE" == "true" && -n "${EXISTING_ISOVALS[$isovalue]}" ]]; then
                    continue
                fi
                echo "${dataset}|${isovalue}|${config}|${config_flags}" >> "$JOBS_FILE"
                ((jobs_count++))
            done
        else
            # Normal mode: generate all isovalues
            while read -r isovalue; do
                # Skip if already computed (resume mode)
                if [[ "$RESUME_MODE" == "true" && -n "${EXISTING_ISOVALS[$isovalue]}" ]]; then
                    continue
                fi
                
                # Limit by MAX_ISOVALUES if set
                if [[ $MAX_ISOVALUES -gt 0 && $jobs_count -ge $MAX_ISOVALUES ]]; then
                    break
                fi
                
                # Format: dataset|isovalue|config_name|config_flags
                echo "${dataset}|${isovalue}|${config}|${config_flags}" >> "$JOBS_FILE"
                ((jobs_count++))
            done < <(generate_isovalues "$start" "$end" "$inc")
        fi
        
        unset EXISTING_ISOVALS
        
        if [[ $jobs_count -eq 0 ]]; then
            echo "All isovalues already computed (resume mode). Skipping."
            echo ""
            continue
        fi
        
        echo "Running $jobs_count tests..."
        
        if [[ "$USE_PARALLEL" == "true" ]]; then
            PARALLEL_OPTS="--bar --eta"
            [[ -n "$NUM_JOBS" ]] && PARALLEL_OPTS+=" -j $NUM_JOBS"
            
            # Run parallel and append to CSV
            env SHELL=/bin/bash parallel $PARALLEL_OPTS run_single_test :::: "$JOBS_FILE" >> "$csv_file"
        else
            # Sequential execution with progress
            current=0
            while IFS= read -r job_line; do
                ((current++))
                local_percent=$((current * 100 / jobs_count))
                printf "\rProgress: %d/%d (%d%%)  " "$current" "$jobs_count" "$local_percent"
                
                run_single_test "$job_line" >> "$csv_file"
            done < "$JOBS_FILE"
            printf "\r%-60s\r" ""
        fi
        
        rm -f "$JOBS_FILE"
        
        # Generate summary for this dataset-config
        summary_file="$OUTPUT_DIR/summary/${dataset}_${config}_summary.txt"
        {
            echo "======================================================================"
            echo "Summary: $dataset + $config"
            echo "======================================================================"
            echo ""
            
            # Use awk to compute statistics
            awk -F',' 'NR > 1 && $4 != "N/A" {
                count++
                
                # Timing
                if ($4 ~ /^[0-9.]+$/) {
                    time_sum += $4
                    if (time_min == "" || $4 < time_min) time_min = $4
                    if (time_max == "" || $4 > time_max) time_max = $4
                }
                
                # Self-intersections  
                if ($10 ~ /^[0-9]+$/) {
                    selfi_sum += $10
                    selfi_count++
                    if ($10 > 0) selfi_nonzero++
                    if (selfi_max == "" || $10 > selfi_max) { selfi_max = $10; selfi_max_iso = $3 }
                }
                
                # Angles
                if ($11 ~ /^[0-9.]+$/) {
                    if (min_angle_min == "" || $11 < min_angle_min) { min_angle_min = $11; min_angle_iso = $3 }
                    if (min_angle_max == "" || $11 > min_angle_max) min_angle_max = $11
                }
                if ($12 ~ /^[0-9.]+$/) {
                    if (max_angle_max == "" || $12 > max_angle_max) { max_angle_max = $12; max_angle_iso = $3 }
                }
                
                # Non-manifold
                if ($7 ~ /^[0-9]+$/) nm_edges_sum += $7
                if ($8 ~ /^[0-9]+$/) nm_verts_sum += $8
                
                # Orientation failures
                if ($9 == "FAIL") orient_fail++
                
                # Mesh size
                if ($5 ~ /^[0-9]+$/) {
                    vert_sum += $5
                    if (vert_min == "" || $5 < vert_min) vert_min = $5
                    if (vert_max == "" || $5 > vert_max) vert_max = $5
                }
                if ($6 ~ /^[0-9]+$/) {
                    tri_sum += $6
                    if (tri_min == "" || $6 < tri_min) tri_min = $6
                    if (tri_max == "" || $6 > tri_max) tri_max = $6
                }
                
                # Angle buckets
                if ($13 ~ /^[0-9]+$/) angle5_sum += $13
                if ($14 ~ /^[0-9]+$/) angle10_sum += $14
                if ($15 ~ /^[0-9]+$/) angle15_sum += $15
                if ($16 ~ /^[0-9]+$/) angle16_sum += $16
                if ($17 ~ /^[0-9]+$/) angle30_sum += $17
            }
            END {
                if (count > 0) {
                    printf "Total tests:               %d\n", count
                    printf "\n"
                    printf "Timing (seconds):\n"
                    printf "  Min:                     %.3f\n", time_min
                    printf "  Max:                     %.3f\n", time_max
                    printf "  Avg:                     %.3f\n", time_sum / count
                    printf "  Total:                   %.3f\n", time_sum
                    printf "\n"
                    printf "Mesh Size:\n"
                    printf "  Vertices (min/max):      %d / %d\n", vert_min, vert_max
                    printf "  Triangles (min/max):     %d / %d\n", tri_min, tri_max
                    printf "\n"
                    printf "Self-Intersections:\n"
                    printf "  Tests with zero:         %d / %d (%.1f%%)\n", selfi_count - selfi_nonzero, selfi_count, (selfi_count - selfi_nonzero) * 100 / selfi_count
                    printf "  Max count:               %d (at iso=%s)\n", selfi_max, selfi_max_iso
                    printf "  Total sum:               %d\n", selfi_sum
                    printf "\n"
                    printf "Angles:\n"
                    printf "  Worst min angle:         %.2f (at iso=%s)\n", min_angle_min, min_angle_iso
                    printf "  Best min angle:          %.2f\n", min_angle_max
                    printf "  Worst max angle:         %.2f (at iso=%s)\n", max_angle_max, max_angle_iso
                    printf "\n"
                    printf "Non-manifold totals:\n"
                    printf "  Edges sum:               %d\n", nm_edges_sum
                    printf "  Vertices sum:            %d\n", nm_verts_sum
                    printf "  Orientation failures:    %d\n", orient_fail
                    printf "\n"
                    printf "Angle bucket totals (polygons with angle <= N):\n"
                    printf "  <=5:   %d\n", angle5_sum
                    printf "  <=10:  %d\n", angle10_sum
                    printf "  <=15:  %d\n", angle15_sum
                    printf "  <=20:  %d\n", angle16_sum
                    printf "  <=30:  %d\n", angle30_sum
                }
            }' "$csv_file"
        } > "$summary_file"
        
        echo ""
        echo "Summary written to: $summary_file"
        echo ""
    done
done

echo "======================================================================"
echo "                        ALL TESTS COMPLETED"
echo "======================================================================"
echo ""
echo "Results:"
echo "  CSV files:    $OUTPUT_DIR/csv/"
echo "  Summaries:    $OUTPUT_DIR/summary/"
echo "  Error log:    $ERROR_LOG"
echo ""

# Report errors if any
if [[ -s "$ERROR_LOG" ]]; then
    error_count=$(wc -l < "$ERROR_LOG")
    echo "Note: $error_count error lines logged. See $ERROR_LOG"
else
    rm -f "$ERROR_LOG"
    echo "All tests completed successfully (no errors)."
fi
