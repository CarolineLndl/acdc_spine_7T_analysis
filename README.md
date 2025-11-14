# acdc_spine_7T_func

## Overview
Processing of spinal cord functional data acquired at 7T.

---

## 1. Getting Started

### 1.1 Dependencies 🔗
Your environment should include:
- Spinal Cord Toolbox 7.1
- Conda environment: `acdc_spine_7T_analysis/config/environment.yml`
- FSL
- dcm2niix
- MATLAB (for denoising step only)

For an example on how to set up the environment, see: `acdc_spine_7T_analysis/config/spine_7T_env_032024.sh`

### 1.2 Data organization 📑
Files are organized according to the BIDS standard:
<details>
<summary>Click to expand folder tree</summary>

```
├── derivatives
│   ├── acdc_spine_7T_project
│   │   ├── acdc_spine_7T_analysis  # GitHub repository
│   │   │   ├── ...
│   │   ├── manual  # Manually corrected files
│   │   │   └── sub-100
│   │   │       ├── anat
│   │   │       │   ├── sub-100_T2s_space-orig_label-ivd_mask.nii.gz
│   │   │       │   └── sub-100_T2star_space-orig_label-ivd_mask.nii.gz
│   │   │       └── func
│   │   │           ├── task-motor_acq-shimBase+3mm
│   │   │           │   ├── sub-100_task-motor_acq-shimBase+3mm_bold_moco_mean_seg.nii.gz
│   │   │           │   ├── sub-100_task-motor_acq-shimBase+3mm_bold_tmean_centerline.csv
│   │   │           │   └── sub-100_task-motor_acq-shimBase+3mm_bold_tmean_centerline.nii.gz
│   │   │           ├── ...
│   │   └── preprocessing
│   │       ├── nov25
│   │       │   ├── QC  # QC reports
│   │       │   │   ├── ...
│   │       │   └── sub-100
│   │       │       ├── anat
│   │       │       │   ├── sct_deepseg
│   │       │       │   │   ├── sub-100_T2star_seg.json
│   │       │       │   │   └── sub-100_T2star_seg.nii.gz
│   │       │       │   ├── sct_label_vertebrae
│   │       │       │   │   ...
│   │       │       │   ├── sct_register_to_template
│   │       │       │   │   ...
│   │       │       │   └── sub-100_T2star.nii.gz
│   │       │       └── func
│   │       │           ├── task-motor_acq-shimBase+3mm
│   │       │           │   ├── sct_deepseg
│   │       │           │   │   ├── sub-100_task-motor_acq-shimBase+3mm_bold_moco_mean_seg.json
│   │       │           │   │   └── sub-100_task-motor_acq-shimBase+3mm_bold_moco_mean_seg.nii.gz
│   │       │           │   ├── sct_fmri_moco
│   │       │           │   │   ...
│   │       │           │   ├── sct_get_centerline
│   │       │           │   │   ...
│   │       │           │   ├── sct_propseg
│   │       │           │   │   ...
│   │       │           └── task-motor_acq-shimSlice+3mm
│   │       │               ...
│   │       └── ...  # Other processing versions
├── rawdata  # BIDS-compliant raw data
│   ├── dataset_description.json
│   ├── sub-100
│   │   ├── anat
│   │   │   ├── sub-100_T2star.json
│   │   │   └── sub-100_T2star.nii.gz
│   │   └── func
│   │       ├── sub-100_task-motor_acq-shimBase+3mm_bold.json
│   │       ├── sub-100_task-motor_acq-shimBase+3mm_bold.nii.gz
│   │       └── ...
├── sourcedata  # Original DICOM and behavioral data
│   ├── sub-100
│   │   ├── behav
│   │   │   ├── *.csv
│   │   │   ├── *.log
│   │   │   ├── *.psydat
│   │   │   └── ...
│   │   ├── mri
│   │   │   ├── 01-localizer_iso_ND
│   │   │   │   ├── *.dcm
│   │   │   │   └── ...
│   │   │   ├── ...
│   │   └── pmu
│   │       ├── ...
```

</details>

### 1.3 Get data into BIDS format 🗂️
Use `dcm2bids` to convert raw data:

```bash
cd $project_dir/acdc_spine_7T_analysis/code/

dcm2bids -d $main_dir/sourcedata/sub-$ID/mri/ \
          -p $ID \
          -c $project_dir/acdc_spine_7T_analysis/config/config_bids_6Nov25.txt \
          -o $main_dir/rawdata/
```

- `$ID` is the subject ID (e.g., 103)
- For full data conversion instructions, see: `/acdc_spine_7T_analysis/code/00_convert_data.sh`

---

## 2. Analysis Pipelines 📊
Files for preprocessing are in this repository.

- **code/**: Functions used by notebooks. Do not modify scripts unless necessary.
- **config/**: Configuration files for paths and parameters.
  - `config_preprocess_spine7T.json` is used by `01_spine7T_preprocessing.ipynb`
    - Modify paths line 1-5 as needed
    - Specify the participant IDs to process line 12
    - Specify the experiment tasks/acquisitions line 17-18
    - Specify file specificities for each subject if needed line 60-65 (*e.g.,* if extra run specify only the one to process)

  - `participants.tsv` contains demographical information and important info for preprocessing (*e.g.,* slice number for vertebrae labeling initiation)
- **notebooks/**: Dedicated notebooks for each analysis step. Use `verbose=True` to check outputs. Completed notebooks can be saved in HTML.
- **template images**: Used for analyses; do not modify.

### 2.1 Preprocessing 🤯
- Update `config_preprocess_spine7T.json`
- Two options to run preprocessing:
  1. **Notebook**: `notebooks/01_spine7T_preprocessing.ipynb` (recommended for QC and step-by-step checks and manual adjustments)
  2. **Script** (main path should be manually changed): `bash code/run_batch_preprocessing.sh` (runs steps automatically, less flexible)

>*You can for exemple run the script and then manually check and correct specific steps in the notebook. 
⚠️ Each step manually modified will imply that all subsequent steps need to be re-run. </span>* 
  

##### 👉 Visual check and manual corrections
<details>
<summary>Click to expand folder tree</summary>
  - **I.a Motion correction (mask)**: ✏️
  check the automatic centerline detection and the mask used for motion correction, if needed, manually correct the centerline with:
  ```
  ctrl_sc_files_, mask_sc_files_=preprocess_Sc.moco_mask(ID=ID,i_img=mean_func_f[ID][tag][run_nb],
                                                                       radius_size=25,task_name=tag,
                                                                       manual=True,
                                                                       redo_ctrl=True,
                                                                       redo_mask=True,
                                                                       verbose=verbose)
  ```

  The output files can be found in:
  ```
  $main_dir/acdc_spine_7T_project/manual/sub-<ID>/func/
      └── <task*_acq*>/
          ├── sub-<ID>_<task_acq>_bold_tmean_centerline.csv
          └── sub-<ID>_<task_acq>_bold_tmean_centerline.nii.gz

  ```
 
  - **II Segmentation** ✏️
  Check the segmentation results, if needed, manually correct the segmentation in fsleyes using the anatomical image or mean functional image as background.
 When saving the corrected segmentation, make sure to keep the same name as the original segmentation file but save it in the `manual` folder:
  ```
  $main_dir/acdc_spine_7T_project/manual/sub-<ID>/func
      └── <task*_acq*>/
          └── sub-<ID>_<task_acq>_bold_moco_mean_seg.nii.gz
  ``` 

  - **III Labeling of inter vertebral disk** ✏️
  Check the automatic labeling of the inter vertebral disks on the anatomical image, if needed (now default is manual), manually correct the labeling with:
  ```
  vert_labels_files.append(preprocess_Sc.label_vertebrae(ID=ID,
                                                               i_img=raw_anat[ID_nb],
                                                               seg_img=seg_anat_sc_files[ID_nb],
                                                               c="t2",
                                                               initz=f"{z_value},{vert}",auto=False,
                                                               redo=False,
                                                               verbose=verbose))
  ```
  The output files can be found in:
  ```
  $main_dir/acdc_spine_7T_project/manual/sub-<ID>/anat
      └── sub-<ID>_T2star_space-orig_label-ivd.nii.gz
  ``` 
</details>


##### ‼️ What we want to try to improve
> - **I. Motion correction:** try different parameters for the mask size, or different reference images (mean functional, middle volume, etc). Parameters can be easily changed in the Notebook and will be then modified as default parameters in the script.
> - **IV. Registration to template:** check if the parameters for the registration are ok. Parameters can be easily changed in the Notebook and will be then modified as default parameters in the script.

### 2.2 Denoising (TBD) 🧹

### 2.3 First-level Analysis (TBD) 📈
