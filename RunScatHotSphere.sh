#$ -S /bin/bash
#$ -l h_vmem=30G
#$ -l tmem=30G
#$ -l h_rt=6:0:0
#$ -cwd
#$ -j y
#$ -N HotSphereScatter

date
hostname

source ~/.bashrc

source /SAN/inm/tools/cluster/root6/bin/thisroot.sh
source /SAN/inm/tools/cluster/Geant4/bin/geant4.sh
export PATH=/SAN/inm/tools/cluster/bin:$PATH

source /SAN/inm/tools/set_env_cluster.sh

export PATH=~/simind:$PATH

#If you need lots of local i/o create a folder on scratch 0 and stage your data to it.

/usr/bin/time -f 'elapsed time: %es\nmeory usage: %M KB\ncpu usage: %P' \
matlab -nodisplay -nosplash -nodesktop -r 'addpath(genpath('"'.'"'));CallOneForAllScript('"'/SAN/inm/FDG/amoINSERT/Week_2/20190313/PhantomWithSpheres/Correct_Spacer', 0, '/SAN/inm/FDG/amoINSERT/INSERT/PeraFiles/PeraFiles/EnergyWindows/ScatEW_PhantomWithSpheres.mat', 'Scatter_'"');exit;' \
2>&1 \
| tee ./RunHotSpheres_Scatter

qstat -j $JOB_ID

date
