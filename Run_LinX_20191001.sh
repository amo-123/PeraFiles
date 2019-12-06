#$ -S /bin/bash
#$ -l h_vmem=8G
#$ -l tmem=8G
#$ -l h_rt=6:0:0
#$ -pe smp 4
#$ -cwd
#$ -j y
#$ -N LinX_20191001

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
matlab -nodisplay -nosplash -nodesktop -r 'addpath(genpath('"'.'"'));UnkilledOneForAll('"'/SAN/inm/FDG/amoINSERT/LondonData/20191001/X/', 1,'/SAN/inm/FDG/amoINSERT/LondonData/20191001/U/', '~', '20191001_'"');exit;' \
2>&1 \
| tee ./LinX_20191001_log

qstat -j $JOB_ID

date
