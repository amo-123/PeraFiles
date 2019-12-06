#$ -S /bin/bash
#$ -l h_vmem=8G
#$ -l tmem=8G
#$ -l h_rt=6:0:0
#$ -pe smp 4
#$ -cwd
#$ -j y
#$ -N LinY

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
matlab -nodisplay -nosplash -nodesktop -r 'addpath(genpath('"'.'"'));CallOneForAllScript('"'/SAN/inm/FDG/amoINSERT/LondonData/20190909/20191127/', 1, '~', 'ldn_UEW'"');exit;' \
2>&1 \
| tee ./RunLinY_20191127_log

qstat -j $JOB_ID

date
