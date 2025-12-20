#!/bin/bash
# Collect data from run_*.out files with outlier detection and removal

echo "=========================================="
echo "HPC Hybrid Parallel Performance Data"
echo "With Outlier Detection and Removal"
echo "=========================================="
echo ""

# Declare associative arrays to store data
declare -A all_runs       # Store all individual run times
declare -A error_values
declare -A run_counts

# Outlier detection threshold (percentage deviation from median)
OUTLIER_THRESHOLD=50  # If a value deviates more than 50% from median, mark as outlier

# Function to calculate median
calculate_median() {
    local arr=("$@")
    local count=${#arr[@]}
    if [ $count -eq 0 ]; then
        echo "0"
        return
    fi
    
    # Sort array
    IFS=$'\n' sorted=($(sort -n <<<"${arr[*]}"))
    unset IFS
    
    if [ $((count % 2)) -eq 0 ]; then
        # Even number of elements
        local mid1=$((count / 2 - 1))
        local mid2=$((count / 2))
        echo "scale=6; (${sorted[$mid1]} + ${sorted[$mid2]}) / 2" | bc
    else
        # Odd number of elements
        local mid=$((count / 2))
        echo "${sorted[$mid]}"
    fi
}

# Function to detect outliers using median-based method
detect_outliers() {
    local values=("$@")
    local count=${#values[@]}
    
    if [ $count -le 2 ]; then
        # Not enough data points to detect outliers
        echo "${values[@]}"
        return
    fi
    
    local median=$(calculate_median "${values[@]}")
    local cleaned_values=()
    local outliers=()
    
    for val in "${values[@]}"; do
        # Calculate percentage deviation from median
        local deviation=$(echo "scale=2; 100 * ($val - $median) / $median" | bc)
        # Get absolute value
        local abs_deviation=$(echo "$deviation" | sed 's/-//')
        
        # Check if deviation exceeds threshold
        local is_outlier=$(echo "$abs_deviation > $OUTLIER_THRESHOLD" | bc)
        
        if [ "$is_outlier" -eq 1 ]; then
            outliers+=("$val")
        else
            cleaned_values+=("$val")
        fi
    done
    
    # If all values are outliers (shouldn't happen), return original
    if [ ${#cleaned_values[@]} -eq 0 ]; then
        echo "${values[@]}"
    else
        echo "${cleaned_values[@]}"
    fi
}

# First pass: collect all data
echo "Phase 1: Collecting data from output files..."
for file in run_*.out; do
    if [ -f "$file" ]; then
        if [[ $file =~ run_([0-9]+)_([0-9]+)_([0-9]+)\.([0-9]+)\.out ]]; then
            grid=${BASH_REMATCH[1]}
            mpi=${BASH_REMATCH[2]}
            omp=${BASH_REMATCH[3]}
            job_id=${BASH_REMATCH[4]}
            
            key="${grid}_${mpi}_${omp}"
            
            calc_time=$(grep "Calculation time:" "$file" | awk '{print $3}')
            error=$(grep "Final max error:" "$file" | awk '{print $4}')
            
            if grep -q "Successfully completed" "$file" && [ -n "$calc_time" ]; then
                if [ -z "${all_runs[$key]}" ]; then
                    all_runs[$key]=""
                    run_counts[$key]=0
                fi
                all_runs[$key]="${all_runs[$key]} $calc_time"
                error_values[$key]=$error
            fi
        fi
    fi
done

echo "Phase 2: Detecting and removing outliers..."
echo ""

# Declare arrays for cleaned data
declare -A cleaned_runs
declare -A cleaned_counts
declare -A outlier_info

# Detect outliers for each configuration
for key in "${!all_runs[@]}"; do
    runs_array=(${all_runs[$key]})
    original_count=${#runs_array[@]}
    
    if [ $original_count -gt 0 ]; then
        # Detect and remove outliers
        cleaned_array=($(detect_outliers "${runs_array[@]}"))
        cleaned_runs[$key]="${cleaned_array[*]}"
        cleaned_counts[$key]=${#cleaned_array[@]}
        
        # Store outlier information
        outliers_removed=$((original_count - ${#cleaned_array[@]}))
        if [ $outliers_removed -gt 0 ]; then
            outlier_info[$key]="$outliers_removed outlier(s) removed"
        fi
    fi
done

echo "=== Outlier Detection Summary ==="
echo ""
for key in $(echo "${!all_runs[@]}" | tr ' ' '\n' | sort -t_ -k1,1n -k2,2n -k3,3n); do
    IFS='_' read -r grid mpi omp <<< "$key"
    runs_array=(${all_runs[$key]})
    cleaned_array=(${cleaned_runs[$key]})
    
    if [ ${#runs_array[@]} -ne ${#cleaned_array[@]} ]; then
        echo "Configuration: Grid ${grid}, ${mpi} MPI × ${omp} OMP"
        echo "  Original values: ${runs_array[*]}"
        echo "  Cleaned values:  ${cleaned_array[*]}"
        echo "  Outliers removed: $((${#runs_array[@]} - ${#cleaned_array[@]}))"
        
        # Show which values were removed
        for val in "${runs_array[@]}"; do
            if [[ ! " ${cleaned_array[*]} " =~ " ${val} " ]]; then
                echo "  ⚠️  Removed outlier: ${val}s"
            fi
        done
        echo ""
    fi
done

if [ ${#outlier_info[@]} -eq 0 ]; then
    echo "No outliers detected in any configuration."
fi

echo ""
echo "=========================================="
echo "Detailed Results (After Outlier Removal)"
echo "=========================================="
echo ""

for grid_size in 128 256 512; do
    echo "Grid Size: ${grid_size}x${grid_size}x${grid_size}"
    echo "-----------------------------------------------------------"
    
    found_data=false
    for key in $(echo "${!cleaned_runs[@]}" | tr ' ' '\n' | sort -t_ -k2,2n -k3,3n); do
        IFS='_' read -r grid mpi omp <<< "$key"
        
        if [ "$grid" = "$grid_size" ]; then
            found_data=true
            cleaned_array=(${cleaned_runs[$key]})
            count=${#cleaned_array[@]}
            
            if [ $count -gt 0 ]; then
                # Calculate average
                sum=0
                for time in "${cleaned_array[@]}"; do
                    sum=$(echo "$sum + $time" | bc)
                done
                avg_time=$(echo "scale=4; $sum / $count" | bc)
                
                # Calculate standard deviation
                std_dev="N/A"
                if [ $count -gt 1 ]; then
                    sum_sq=0
                    for time in "${cleaned_array[@]}"; do
                        diff=$(echo "$time - $avg_time" | bc)
                        sq=$(echo "$diff * $diff" | bc)
                        sum_sq=$(echo "$sum_sq + $sq" | bc)
                    done
                    variance=$(echo "scale=6; $sum_sq / $count" | bc)
                    std_dev=$(echo "scale=4; sqrt($variance)" | bc)
                fi
                
                echo "  Configuration: ${mpi} MPI × ${omp} OMP (Total: $((mpi * omp)) cores)"
                echo "    Valid runs: $count/${run_counts[$key]}"
                for i in "${!cleaned_array[@]}"; do
                    echo "    Run $((i+1)): ${cleaned_array[$i]}s"
                done
                echo "    Average: ${avg_time}s"
                [ "$std_dev" != "N/A" ] && echo "    Std Dev: ±${std_dev}s"
                [ -n "${outlier_info[$key]}" ] && echo "    ⚠️  ${outlier_info[$key]}"
                echo ""
            fi
        fi
    done
    
    if [ "$found_data" = false ]; then
        echo "  No data"
        echo ""
    fi
done

echo ""
echo "=========================================="
echo "Summary Table (Cleaned Data)"
echo "=========================================="
echo ""

printf "%-10s %-8s %-8s %-10s %-10s %-15s %-15s %-15s\n" \
    "Grid" "MPI" "OMP" "Total" "Valid/Tot" "Avg Time(s)" "Std Dev(s)" "Outliers"
echo "------------------------------------------------------------------------------------------------------"

for key in $(echo "${!cleaned_runs[@]}" | tr ' ' '\n' | sort -t_ -k1,1n -k2,2n -k3,3n); do
    IFS='_' read -r grid mpi omp <<< "$key"
    total=$((mpi * omp))
    
    cleaned_array=(${cleaned_runs[$key]})
    count=${#cleaned_array[@]}
    original_count=${run_counts[$key]:-0}
    
    if [ $count -gt 0 ]; then
        # Calculate average
        sum=0
        for time in "${cleaned_array[@]}"; do
            sum=$(echo "$sum + $time" | bc)
        done
        avg_time=$(echo "scale=4; $sum / $count" | bc)
        
        # Calculate standard deviation
        std_dev="N/A"
        if [ $count -gt 1 ]; then
            sum_sq=0
            for time in "${cleaned_array[@]}"; do
                diff=$(echo "$time - $avg_time" | bc)
                sq=$(echo "$diff * $diff" | bc)
                sum_sq=$(echo "$sum_sq + $sq" | bc)
            done
            variance=$(echo "scale=6; $sum_sq / $count" | bc)
            std_dev=$(echo "scale=4; sqrt($variance)" | bc)
        fi
        
        outliers_removed=$((original_count - count))
        outlier_str="None"
        [ $outliers_removed -gt 0 ] && outlier_str="${outliers_removed} removed"
        
        printf "%-10s %-8s %-8s %-10s %-10s %-15s %-15s %-15s\n" \
            "$grid" "$mpi" "$omp" "$total" "$count/$original_count" "$avg_time" "$std_dev" "$outlier_str"
    fi
done

echo ""
echo "=========================================="
echo "Speedup Analysis (Using Cleaned Data)"
echo "=========================================="
echo ""

for grid_size in 128 256 512; do
    echo "Grid Size: ${grid_size}x${grid_size}x${grid_size}"
    echo "-----------------------------------------------------------"
    
    baseline_time=""
    baseline_config=""
    
    # Find baseline
    for key in $(echo "${!cleaned_runs[@]}" | tr ' ' '\n' | sort -t_ -k2,2n -k3,3n); do
        IFS='_' read -r grid mpi omp <<< "$key"
        
        if [ "$grid" = "$grid_size" ]; then
            cleaned_array=(${cleaned_runs[$key]})
            count=${#cleaned_array[@]}
            
            if [ $count -gt 0 ]; then
                sum=0
                for time in "${cleaned_array[@]}"; do
                    sum=$(echo "$sum + $time" | bc)
                done
                avg_time=$(echo "scale=4; $sum / $count" | bc)
                
                if [ -z "$baseline_time" ]; then
                    baseline_time=$avg_time
                    baseline_config="${mpi} MPI × ${omp} OMP"
                fi
            fi
        fi
    done
    
    if [ -n "$baseline_time" ]; then
        echo "  Baseline: $baseline_config = ${baseline_time}s"
        echo ""
        printf "  %-25s %-15s %-15s %-15s\n" "Configuration" "Avg Time(s)" "Speedup" "Efficiency"
        echo "  -------------------------------------------------------------------------------"
        
        for key in $(echo "${!cleaned_runs[@]}" | tr ' ' '\n' | sort -t_ -k2,2n -k3,3n); do
            IFS='_' read -r grid mpi omp <<< "$key"
            
            if [ "$grid" = "$grid_size" ]; then
                cleaned_array=(${cleaned_runs[$key]})
                count=${#cleaned_array[@]}
                
                if [ $count -gt 0 ]; then
                    sum=0
                    for time in "${cleaned_array[@]}"; do
                        sum=$(echo "$sum + $time" | bc)
                    done
                    avg_time=$(echo "scale=4; $sum / $count" | bc)
                    
                    speedup=$(echo "scale=2; $baseline_time / $avg_time" | bc)
                    total_cores=$((mpi * omp))
                    efficiency=$(echo "scale=1; 100 * $speedup / $total_cores" | bc)
                    printf "  %-25s %-15s %-15s %-15s\n" \
                        "${mpi} MPI × ${omp} OMP" "$avg_time" "${speedup}x" "${efficiency}%"
                fi
            fi
        done
    else
        echo "  No data"
    fi
    echo ""
done

echo "=========================================="
echo "Analysis Complete!"
echo "=========================================="
echo ""
echo "Notes:"
echo "  - Outlier threshold: ±${OUTLIER_THRESHOLD}% deviation from median"
echo "  - Outliers are automatically detected and excluded from averages"
echo "  - At least 2 valid runs are kept for each configuration"
