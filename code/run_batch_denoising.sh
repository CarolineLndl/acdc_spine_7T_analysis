mkdir -p log
cd log
timestamp=$(date +"%Y%m%d_%H%M%S")

nohup python ../code/02_spine7T_denoising.py ../log/ > "nohup_denoising_${timestamp}.out" 2>&1 &