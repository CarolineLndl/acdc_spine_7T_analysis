main_dir=/cerebro/cerebro1/dataset/spine_7T/derivatives/acdc_spine_7T_project/acdc_spine_7T_analysis/

cd $main_dir"/log/"
nohup python $main_dir/code/01_spine7T_preprocessing.py $main_dir"/log" .pynohup.out &
