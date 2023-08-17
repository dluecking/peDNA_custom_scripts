#!/bin/bash -l

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=three_stations   # Job name
#SBATCH --output=out/%x.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=out/%x.%j.err # Stderr (%j expands to jobId)
#SBATCH --time=3-00:00:00   # walltime
#SBATCH --mem=100G   # memory per NODE
#SBATCH --cpus-per-task=24
#SBATCH --array=1-9

source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
conda activate bbmap

station=`cat helper_files/list_of_stations.txt | sed -n ${SLURM_ARRAY_TASK_ID}p`

for MAG in `cat helper_files/list_of_all_MAGs.txt`;
do
	bbmap.sh \
	in=reads/${station}.true.mvome.reads1.1M.fq.gz \
	in2=reads/${station}.true.mvome.reads2.1M.fq.gz \
	ref=contigs/${MAG} \
	scafstats=mappings/${station}_vs_${MAG}_scafstat.txt \
	nodisk \ 
	t=24
done
