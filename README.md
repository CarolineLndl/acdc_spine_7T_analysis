# Project: Spinal cord fMRI analysis at 7T

## Overview
Processing of spinal cord functional data acquired at 7T.

---

## 1. Getting Started

### 1.1 Dependencies 🔗
Your environment should include:
- Python (3.10.14 was used)
- Spinal Cord Toolbox 7.1
- Conda environment: `spine_7T_analysis/config/requirements.txt`
- FSL
- dcm2niix
- MATLAB (for denoising step only)

For an example on how to set up the environment, see: `spine_7T_analysis/config/spine_7T_env_032024.sh`


<details>
<summary>👉 How to install dependencies</summary>

#### a. Install toolboxes
**Toolboxes for preprocessing**
- Spinal Cord Toolbox 7.1: [Installation instructions](https://spinalcordtoolbox.com/en/latest/user_section/installation.html)
- FSL: see here [Installation instructions](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)

**Toolboxes for denoising:**
- Verify wich version of MATLAB is compatible with your Python version (*vis versa*): see here [Compatibility table](https://www.mathworks.com/support/requirements/python-compatibility.html)
- Install MATLAB: see here [Installation instructions](https://www.mathworks.com/help/install/)
- Install MATLAB engine for Python: see here [Installation instructions](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html)

```bash
# Example for MATLAB R2023b
LD_LIBRARY_PATH="/export01/local/matlab23b/sys/os/glnxa64:$toolbox_home/libraries"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/export01/local/matlab23b/bin/" # The LD_LIBRARY_PATH environment variable tells the system where to find shared libraries
cd /export/local/matlab23b/extern/engines #Navigate to MATLAB Folder, engines subfolder
python -m pip install matlabengine==23.2.10 #Install MATLAB engine for Python
```

#### b. Create the conda environment 

Make sure conda is installed: see here [Installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
Create the appropriate conda environment:
*If you just what to run the preprocessing you can remove matlabengine from the requirements.txt file.*

```bash
conda create --name spine_7T_env_py10 python=3.10
conda activate spine_7T_env_py10
pip install -r config/requirements.txt
```


#### c. Load your environment
This step should be done each time you open a new terminal session to work on the project. It will set up all necessary path for SCT, FSL and Matlab and activate the conda environment.
> ⚠️ *Qt version conflict: SCT and MATLAB may use incompatible Qt libraries. If you don’t need MATLAB, consider commenting out the MATLAB path in the setup script to avoid errors. If you need MATLAB for denoising, uncomment the MATLAB path, but be aware that Qt-related errors may appear when using SCT manually.*

```bash
source spine_7T_analysis/config/spine_7T_env_17112025.sh
```

</details>

### 1.2 Data organization 📑
Files are organized according to the BIDS standard:
<details>
<summary>Click to expand folder tree</summary>

```
├── spine_7T_analysis  # GitHub repository
│   ├── code
│     ├── convert_data
│     │   ├── 00_convert_mriData.sh
│     │   └── ...
│     ├── ...
│   ├── config
│     ├── config_preprocess_spine7T.json
│     ├── participants.tsv
│     └── ...
│   ├── notebooks
│     ├── 01_spine7T_preprocessing.ipynb
│     └── ...
│   ├── template_images
│     ├── ...
│   └── log
│       ├── ...
├── derivatives
│   ├── spine_7T_project
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
│   │       ├── sub-100_task-rest_acq-shimBase+3mm_physio.json
│   │       ├── sub-100_task-rest_acq-shimBase+3mm_physio.tsv.gz
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
#### Convert mri data
Use `dcm2bids` to convert raw mri data:

```bash
cd $project_dir/acdc_spine_7T_analysis/code/

dcm2bids -d $main_dir/sourcedata/sub-$ID/mri/ \
          -p $ID \
          -c $project_dir/acdc_spine_7T_analysis/config/config_bids_6Nov25.txt \
          -o $main_dir/rawdata/
```

- `$ID` is the subject ID (e.g., 103)
- For full data conversion instructions, see: `/acdc_spine_7T_analysis/code/convert_data/00_convert_mriData.sh`

#### Convert physio data
Use `/acdc_spine_7T_analysis/code/convert_data/00_convert_physioData.sh` to convert raw physio data into BIDS format.

```bash
cd $project_dir/acdc_spine_7T_analysis/code/convert_data/
bash 00_convert_physioData.sh
``` 

---

## 2. Analysis Pipelines 📊
Files for preprocessing are in this repository.

- **code/**: Functions used by notebooks. Do not modify scripts unless necessary.
  - **convert_data/**: Scripts to convert raw mri and physio data into BIDS format.
- **config/**: Configuration files for paths and parameters.
  - `config_preprocess_spine7T.json` is used by `01_spine7T_preprocessing.ipynb`
    - Modify paths line 1-6 as needed
    - Specify the participant IDs to process line 13
    - Specify the experiment tasks/acquisitions line 18-19
    - Specify file specificities for each subject if needed line 60-65 (*e.g.,* if extra run specify only the one to process)

  - `participants.tsv` contains demographical information and important info for preprocessing (*e.g.,* slice number for vertebrae labeling initiation)
- **notebooks/**: Dedicated notebooks for each analysis step. Use `verbose=True` to check outputs. Completed notebooks can be saved in HTML.
- **template images**: Used for analyses; do not modify.
- **log**: Log files generated during processing run from bash script.

### 2.1 Preprocessing 🤯
- Update `config_preprocess_spine7T.json`
- Two options to run preprocessing:
  > 1. **Script** `bash code/run_batch_preprocessing.sh` (main path should be manually changed) 
  ▸ runs steps automatically: recommanded to run all steps at once 
  ▸ less flexible than the notebook
  ▸ By default all the steps are rerun even if some outputs already exist. If manual corrections were made, these files will be used as input for subsequent steps.
  > 2. **Notebook** `notebooks/01_spine7T_preprocessing.ipynb`
  ▸ run one step at a time: recommended for QC and step-by-step checks and manual adjustments
  ▸ more flexible than the script: parameters can be easily changed and tested
  ▸ If manual corrections were made, these files will be used as input for subsequent steps.

  

*You can for exemple run the script and then manually check and correct specific steps in the notebook. 
⚠️ Each step manually modified will imply that all subsequent steps need to be re-run. </span>* 
  

##### 👉 Visual check and manual corrections
<details>
<summary>For more details, click to expand </summary>

  - **I.a Motion correction (mask)** : ✏️
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

Should be run after preprocessing.
- Update `config_preprocess_spine7T.json`
- Two options to run preprocessing:
  1. **Notebook**: `notebooks/02_spine7T_denoising.ipynb` (recommended for QC and step-by-step checks and manual adjustments)
  2. **Script** (main path should be manually changed): `bash code/run_batch_denoising.sh` (runs steps automatically, less flexible)

No manual corrections are requiered for this step.

### 2.3 First-level Analysis (TBD) 📈
