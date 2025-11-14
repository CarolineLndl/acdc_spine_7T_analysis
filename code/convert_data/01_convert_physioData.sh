main_dir=/cerebro/cerebro1/dataset/spine_7T/
project_dir=$main_dir"derivatives/acdc_spine_7T_project/"
cd  $project_dir"/acdc_spine_7T_analysis/config/"
source spine_7T_env_032024.sh

ID=100

#Compress physio files ------------------------------------------------
cd $main_dir"/sourcedata/sub-"$ID"/pmu/"

# rename the files
declare -A associationMRISession=(
    ["task-rest_acq-shimBase+3mm"]="baseline_tsnr"
    ["task-motor_acq-shimBase+3mm"]="baseline_task_run1"
    ["task-motor_acq-shimBase+3mm-run02"]="baseline_task_run2"

    ["task-rest_acq-shimSlice+3mm"]="f0xyz_tsnr"
    ["task-motor_acq-shimSlice+3mm"]="f0xyz_task_run1"
    ["task-motor_acq-shimSlice+3mm-run02"]="f0xyz_task_run2"

    ["task-rest_acq-shimSlice+1mm+sms2"]="1mm_slice_SMS2_f0xyz_tsnr"
    ["task-motor_acq-shimSlice+1mm+sms2"]="1mm_slice_SMS2_f0xyz_task"
    ["task-rest_acq-shimBase+1mm+sms2"]="1mm_slice_SMS2_baseline_tsnr"
)
cd "$physio_dir" || exit

for task in "${!associationMRISession[@]}"; do
    base="${associationMRISession[$task]}"
    for ext in ext resp puls; do
        for prefix in ext pmu pulse; do
            oldfile=$(ls "${prefix}_signalep2d_bold_${base}.${ext}" 2>/dev/null)
            if [[ -f "$oldfile" ]]; then
                newfile="sub-100_${task}_bold.$ext"
                echo "Renaming $oldfile → $newfile"
                mv "$oldfile" "$newfile"
            fi
        done
    done
done



EXTENSIONS=("ext" "puls" "resp")
# Collect basenames from the known extensions only
cd "$main_dir/sourcedata/sub-$ID/pmu/" || exit 1

# Extract the ID from ext_*.ext files
basenames=$(ls *.ext 2>/dev/null | sed 's/^//' | sed 's/\.ext$//')

for base in $basenames; do
    f_ext="${base}.ext"
    f_puls="${base}.puls"
    f_resp="${base}.resp"

    # Check if all 3 files exist
    if [[ -f "$f_ext" && -f "$f_puls" && -f "$f_resp" ]]; then
        tar -czf "${base}.tar.gz" "$f_ext" "$f_puls" "$f_resp"
        echo "Created: ${base}.tar.gz"
    else
        echo "Skipping $base — missing one or more files."
    fi
done


#### rename manually the files to match BIDS convention


#Convert physio to BIDS
cd "$project_dir/acdc_spine_7T_analysis/code/convert_data/"
for archive in "$main_dir/sourcedata/sub-$ID/pmu/"*.tar.gz; do
    echo "$archive"
    python physio2bids.py -t "$archive" -s "$ID" -o "$main_dir/rawdata/" -v True
done

