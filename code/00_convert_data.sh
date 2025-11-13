main_dir=/cerebro/cerebro1/dataset/spine_7T/
project_dir=$main_dir"derivatives/acdc_spine_7T_project/"

cd  $project_dir"/acdc_spine_7T_analysis/config/"
source spine_7T_env_032024.sh

ID=101
file=acdc_spine_7T_101_20251106_091631571

# Create participant directory
cd /cerebro/cerebro1/dataset/spine_7T/sourcedata/

mkdir "sub-"$ID  # create folder
mkdir "sub-"$ID"/pmu"
mkdir "sub-"$ID"/mri"
mkdir "sub-"$ID"/behav"

# Claim and download the file
#1: login to the BIC, usually at login.bic.mni.mcgill.ca
cd /data/dicom/ 
find_mri -claim "acdc_spine_7T_"$ID

rsync -a "/data/dicom/acdc_spine_7T_"$ID"???" $main_dir"/sourcedata/sub-"$ID"/" # download the data

#rsync -a /data/dicom/acdc_spine_7T_90_20250728_163015298 $main_dir"/sourcedata/sub-"$ID"/" # download the data

# sort dicom files
# PMU and behavioral data should be copyed manually into the pmu and behav folder
cd $project_dir"/acdc_spine_7T_analysis/code/"
python sortDCM.py -d "/cerebro/cerebro1/dataset/spine_7T/sourcedata/sub-"$ID"/"$file -o "/cerebro/cerebro1/dataset/spine_7T/sourcedata/sub-"$ID"/mri/"


#Convert in BIDS
cd $project_dir"/acdc_spine_7T_analysis/code/"
dcm2bids -d $main_dir"/sourcedata/sub-$ID/mri/" -p $ID -c $project_dir"/acdc_spine_7T_analysis/config/config_bids_6Nov25.txt" -o $main_dir"/rawdata/"

#Compress physio files
cd $$main_dir"/sourcedata/$ID/pmu/"
tar -czvf $ID'_rest_run-06.tar.gz' * 

#Convert physio to BIDS
cd $project_dir"/acdc_spine_7T_analysis/code"
python physio2bids.py -t $main_dir"sourcedata/"$ID"/pmu/"$ID"_rest"$func_run".tar.gz" -s $ID -o $main_dir"/rawdata/" -v True


