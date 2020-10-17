#!/bin/bash
file_name=$1
echo "#!/bin/bash" > batch_"$file_name".script
echo "#SBATCH --job-name="$file_name"" >> batch_"$file_name".script
echo "#SBATCH --output=out_%j.txt" >> batch_"$file_name".script
echo "#SBATCH --partition=maloney" >> batch_"$file_name".script
echo "#SBATCH --exclude=c3079" >> batch_"$file_name".script
echo "#SBATCH --nodes=1" >> batch_"$file_name".script
echo "#SBATCH --time=10-24:00:00" >> batch_"$file_name".script
echo "#SBATCH --exclusive" >> batch_"$file_name".script
echo "#SBATCH --mem=0" >> batch_"$file_name".script
echo "module load matlab/R2019a" >> batch_"$file_name".script
echo "module load anaconda3/3.7" >> batch_"$file_name".script
echo "module unload anaconda2/2018.12" >> batch_"$file_name".script
echo "matlab -nosplash -nodesktop -r \"parpool('local',12);test_param_solvent_cdc_"$file_name";test_run_solvent_cdc_"$file_name";exit\"" >> batch_"$file_name".script
/usr/bin/sbatch batch_"$file_name".script