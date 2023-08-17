#!/bin/bash -l

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=virus-activity   # Job name
#SBATCH --output=out/%x.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=out/%x.%j.err # Stderr (%j expands to jobId)
#SBATCH --time=4-00:00:00   # walltime
#SBATCH --mem=250G   # memory per NODE
#SBATCH --ntasks=24

source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
conda activate bbmap

for SAMPLE in `cat helper_files/list_of_stations.txt`;
do
	bbmap.sh \
	in=reads/${SAMPLE}.all-input.reads1.1M.fq.gz \
	in2=reads/${SAMPLE}.all-input.reads2.1M.fq.gz\
	ref=top20/virus_activity/combined_viral_seqs.fa \
    rpkm=top20/virus_activity/${SAMPLE}_virus_activity_rpkm.txt \
	t=${SLURM_CPUS_PER_TASK} \
    nodisk
done
