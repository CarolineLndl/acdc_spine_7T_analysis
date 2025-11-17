root_dir=/cerebro/cerebro1/dataset/spine_7T/spine_7T_analysis/

cd $root_dir"/log/"
nohup python $root_dir/code/01_spine7T_preprocessing.py $root_dir"/log" .pynohup.out &