#!/usr/bin/env python
# coding: utf-8

#  Spinal Cord fMRI preprocessing 
# ____________________________________________________
# 
# ### Project: acdc_spine_7T
# ____________________________________________________
# @ author: Caroline Landelle, caroline.landelle@mcgill.ca // landelle.caroline@gmail.com
# July 2025
# 
# Description: 
# This notebook provides code for preprocessing fMRI data of spinal cord acquisition at 7T. 
# 
# Toolbox required: 
# > SpinalCordToolbox  
# > FSL (Python)  
# 
# ____________________________________________________
# 
# nb: The Philipps system includes additional "dummy scans" at the beginning of the acquisition to allow the magnetization to stabilize to a steady state. The dummy scans are not stored, so they will make the banging sound like normal scans but there will be no data associated with them.    
# 
#------------------------------------------------------------------
#------ Initialization
#------------------------------------------------------------------
# Imports
import sys, json, glob, os, re, shutil, argparse
import pandas as pd
from IPython.display import Image, display


# Get the environment variable PATH_CODE
path_code = os.environ.get("PATH_CODE")
path_data = os.environ.get("PATH_DATA")

parser = argparse.ArgumentParser()
parser.add_argument("--ids", nargs='+', default=[""])
parser.add_argument("--verbose", default="False")
parser.add_argument("--manual_centerline", default="False")
parser.add_argument("--auto_vert_labels", default="True")
parser.add_argument("--redo", default="True")
args = parser.parse_args()

IDs = args.ids
verbose = args.verbose.lower() == "true"
manual_centerline = args.manual_centerline.lower() == "true"
auto_vert_labels = args.auto_vert_labels.lower() == "true"
redo = args.redo.lower() == "true"

if isinstance(redo, bool):
    print("redo is boolean")
else:
    print("redo is not boolean")


with open(path_code + '/config/config_spine_7t_fmri.json') as config_file: # the notebook should be in 'xx/notebook/' folder #config_proprio
    config = json.load(config_file) # load config file should be open first and the path inside modified

#Import scripts
sys.path.append(path_code + "/code/") # Change this line according to your directory

from preprocess import Preprocess_main, Preprocess_Sc
import utils

#Initialize codes
Preprocess_main=Preprocess_main(config, IDs=IDs) # initialize the function
preprocess_Sc=Preprocess_Sc(config, IDs=IDs) # initialize the function
ses_name=""

participants_tsv = pd.read_csv(path_code + '/config/participants.tsv', sep='\t')

# initialize directories
preprocessing_dir=os.path.expandvars(config["preprocess_dir"]["main_dir"])
derivatives_dir=os.path.expandvars(config["derivatives_dir"])
manual_dir=os.path.expandvars(config["manual_dir"])

#------------------------------------------------------------------
#------ Preprocessing
#------------------------------------------------------------------
print("redo", redo, flush=True)
print("=== Preprocessing script Start ===", flush=True)
print("Participant(s) included : ", IDs, flush=True)
print("===================================", flush=True)
print("")

for ID_nb,ID in enumerate(IDs):
    print("", flush=True)
    print(f'=== Preprocessing start for :  {ID} ===', flush=True)

    
    #---------------Anat preprocessing ---------------------------------------------------
    raw_anat=glob.glob(preprocessing_dir.format(ID) + "/anat/" + config["preprocess_f"]["anat_raw"].format(ID,"*"))[0]
    
    #------------------------------------------------------------------
    #------ Segmentation of the anatomical image
    #------------------------------------------------------------------
    
    seg_anat_sc_file=preprocess_Sc.segmentation(ID=ID,
                                                i_img=raw_anat,
                                                img_type="anat",
                                                contrast_anat="t2",
                                                redo=redo,
                                                redo_qc=redo, # should be true if you have done manual correction
                                                verbose=verbose)
    
    print(f'=== Anat segmentation : Done {ID} ===', flush=True)

    #------------------------------------------------------------------
    #------ Vertebral labelling
    #------------------------------------------------------------------

    if auto_vert_labels:
        id_info = participants_tsv.loc[participants_tsv['participant_id'] == int(ID), ['anat_disc_ref','anat_disc_ref_z']]
        vert= id_info['anat_disc_ref'].iloc[0][-1]
        z_value = id_info['anat_disc_ref_z'].iloc[0]

    vert_labels_files=preprocess_Sc.label_vertebrae(ID=ID,
                                                    i_img=raw_anat,
                                                    seg_img=seg_anat_sc_file,
                                                    c="t2",
                                                    initz=f"{z_value},{vert}",
                                                    auto=auto_vert_labels,
                                                    redo=redo,
                                                    verbose=verbose)
        

    print(f'=== Anat vertebral labelling : Done {ID} ===', flush=True)

    #------------------------------------------------------------------
    #------ Registration in PAM50
    #------------------------------------------------------------------

    manual_seg_file=f'{manual_dir}/sub-{ID}/anat/' + os.path.basename(seg_anat_sc_file) 
    seg_anat_sc_final_file=manual_seg_file if os.path.exists(manual_seg_file) else seg_anat_sc_file
    param = "step=1,type=seg,algo=centermassrot:step=2,type=im,algo=syn,iter=5,slicewise=1,metric=CC,smooth=0"
    
    warpT2w_PAM50_files=preprocess_Sc.coreg_anat2PAM50(ID=ID,
                                                              i_img=raw_anat,
                                                              seg_img=seg_anat_sc_final_file,
                                                              labels_img=vert_labels_files,
                                                              img_type="t2",
                                                              tag='anat',
                                                              param=param,
                                                              redo=redo,
                                                              verbose=verbose)
    
    print(f'=== Registration anat to PAM50 : Done {ID} ===', flush=True)

    #---------------Func preprocessing ---------------------------------------------------
    #------ Select func data
    for task_name in config["design_exp"]["task_names"]:
        for acq_name in config["design_exp"]["acq_names"]:
            tag="task-" + task_name + "_acq-" + acq_name
            json_f=glob.glob(os.path.expandvars(config["raw_dir"]) + f'/sub-{ID}/func/sub-{ID}_{tag}_*bold.json')
            raw_func=glob.glob(os.path.expandvars(config["raw_dir"]) + f'/sub-{ID}/func/sub-{ID}_{tag}_*bold.nii.gz')  
            o_dir=preprocessing_dir.format(ID)+  "/func/" +tag + '/'

            for func_file in raw_func:
                # Check run number if multiple run exists
                match = re.search(r"_?(run-\d+)", func_file)
                if match:
                    run_name=match.group(1)
                    print(run_name)
                else:
                    run_name=""

                #------------------------------------------------------------------
                #------ Create mask around the cord for moco
                #------------------------------------------------------------------
                o_img=o_dir +  os.path.basename(func_file).split(".")[0] + "_tmean.nii.gz"
                mean_func_f=utils.tmean_img(ID=ID,i_img=func_file,o_img=o_img,verbose=False)
                ctrl_sc_file, mask_sc_file=preprocess_Sc.moco_mask(ID=ID,
                                                                       i_img=mean_func_f,
                                                                       radius_size=25,
                                                                       task_name=tag,
                                                                       manual=manual_centerline,
                                                                       redo_ctrl=redo,
                                                                       redo_mask=redo,
                                                                       verbose=verbose)

                print(mask_sc_file)
                print(f'=== Moco masks : Done  {ID} {tag} {run_name} ===', flush=True)
                      
                #------------------------------------------------------------------
                #------ Run moco
                #------------------------------------------------------------------
                params = 'poly=0,smooth=1,metric=MeanSquares,gradStep=1,sampling=0.2'
                moco_f,moco_mean_f,qc_dir=preprocess_Sc.moco(ID=ID,
                                                               i_img=func_file,
                                                               mask_img=mask_sc_file,
                                                               task_name=tag,
                                                               run_name=run_name,
                                                               params=params,
                                                               verbose=verbose,
                                                               redo=redo)

                print(f'=== Moco : Done  {ID} {tag} {run_name} ===', flush=True)

                #------------------------------------------------------------------
                #------ Run func cord and CSF segmentation
                #------------------------------------------------------------------
                # Cord segmentation
                seg_func_sc_file=preprocess_Sc.segmentation(ID=ID,
                                                             i_img=moco_mean_f,
                                                             task_name=tag,
                                                             img_type="func",
                                                             mask_qc=mask_sc_file,
                                                             redo=redo,
                                                             redo_qc=redo, # should be true if you have done manual correction
                                                             verbose=verbose)
                # csf segmentation
                preprocess_Sc.segmentation(ID=ID,
                                           i_img=moco_mean_f,
                                           task_name=tag,contrast_anat="t2s",
                                           img_type="func", 
                                           tissue="csf",
                                           redo_qc=redo, # should be true if you have done manual correction
                                           redo=redo,
                                           verbose=verbose)
                
                print(f'=== Func segmentation : Done  {ID} {tag} {run_name} ===', flush=True)

                #------------------------------------------------------------------
                #------ Registration in PAM50
                #------------------------------------------------------------------
                param="step=1,type=seg,algo=slicereg,metric=MeanSquares,smooth=2:step=2,type=im,algo=syn,metric=CC,iter=3,slicewise=1"
                func2PAM50_dir=preprocess_Sc.coreg_img2PAM50(ID=ID,
                                                             i_img=moco_mean_f,
                                                             i_seg=seg_func_sc_file,
                                                             task_name=tag,
                                                             run_name=run_name,
                                                             initwarp=warpT2w_PAM50_files[1],
                                                             initwarpinv=warpT2w_PAM50_files[0],
                                                             param=param,
                                                             redo=redo,
                                                             verbose=verbose)
                
                print(f'=== Func registration : Done  {ID} {tag} {run_name} ===')

        
    print(f'=== Preprocessing done for : {ID} ===', flush=True)
    print("=========================================", flush=True)
