root_dir=/cerebro/cerebro1/dataset/spine_7t/

cd  $root_dir"/spine_7t_fmri_analysis/config/"
source spine_7T_env.sh

ID=106


# Create participant directory
cd $root_dir"/spine_7t_fmri_data/sourcedata/"
mkdir "sub-"$ID  # create folder
mkdir "sub-"$ID"/pmu"
mkdir "sub-"$ID"/mri"
mkdir "sub-"$ID"/behav"

# Claim and download the file
#1: login to the BIC, usually at login.bic.mni.mcgill.ca
cd /data/dicom/ 
find_mri "acdc_spine_7T_"$ID

cd $root_dir"/spine_7t_fmri_data/sourcedata/"
file=acdc_spine_7T_106_20251125_134832384 #acdc_spine_7T_105_20251125_134948810
rsync -a "/data/dicom/"$file $root_dir"/spine_7t_fmri_data/sourcedata/sub-"$ID"/" # download the data

#rsync -a /data/dicom/acdc_spine_7T_90_20250728_163015298 $main_dir"/sourcedata/sub-"$ID"/" # download the data

# sort dicom files
# PMU and behavioral data should be copyed manually into the pmu and behav folder
cd $root_dir"/spine_7t_fmri_analysis/code/convert_data/"
python sortDCM.py -d $root_dir"//spine_7t_fmri_data/sourcedata/sub-"$ID"/"$file -o $root_dir"/spine_7t_fmri_data/sourcedata/sub-"$ID"/mri/"


#Convert in BIDS
cd $root_dir"/spine_7t_fmri_analysis/code/convert_data/"
dcm2bids -d $root_dir"/spine_7t_fmri_data/sourcedata/sub-$ID/mri/" -p $ID -c $root_dir"/spine_7t_fmri_analysis/config/config_bids.txt" -o $root_dir"/spine_7t_fmri_data/"
