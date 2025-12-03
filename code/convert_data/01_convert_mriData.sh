

cd  ${PATH_CODE}"/config/"
source spine_7T_env.sh

ID=106


# Create participant directory
cd ${PATH_DATA}"/sourcedata/"
mkdir "sub-"$ID  # create folder
mkdir "sub-"$ID"/pmu"
mkdir "sub-"$ID"/mri"
mkdir "sub-"$ID"/behav"

# Claim and download the file, this step is specific to the Neuro Brain Imaging Centre
#1: login to the BIC, usually at login.bic.mni.mcgill.ca
cd /data/dicom/ 
find_mri "acdc_spine_7T_"$ID

cd ${PATH_DATA}"/sourcedata/"
file=your_ID_data_folder #acdc_spine_7T_105_20251125_134948810
rsync -a "/data/dicom/"$file ${PATH_DATA}"/sourcedata/sub-"$ID"/" # download the data

#rsync -a /data/dicom/acdc_spine_7T_90_20250728_163015298 $main_dir"/sourcedata/sub-"$ID"/" # download the data

# sort dicom files
# PMU and behavioral data should be copyed manually into the pmu and behav folder
cd ${PATH_CODE}"/code/convert_data/"
python sortDCM.py -d ${PATH_DATA}"/sourcedata/sub-"$ID"/"$file -o ${PATH_DATA}"/sourcedata/sub-"$ID"/mri/"


#Convert in BIDS
cd ${PATH_CODE}"/code/convert_data/"
dcm2bids -d ${PATH_DATA}"/sourcedata/sub-$ID/mri/" -p $ID -c ${PATH_CODE}"/config/config_bids.txt" -o $root_dir"/spine_7t_fmri_data/"
