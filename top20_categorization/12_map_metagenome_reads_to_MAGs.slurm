#!/bin/bash -l

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=metagenome_mapping   # Job name
#SBATCH --output=out/metagenome_mapping/%x.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=out/metagenome_mapping/%x.%j.err # Stderr (%j expands to jobId)
#SBATCH --time=10:00:00   # walltime
#SBATCH --mem=50G   # memory per NODE
#SBATCH --cpus-per-task=16
#SBATCH --partition=CLUSTER
#SBATCH --array=3


source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
conda activate bbmap

STATION=`cat helper_files/list_of_stations.txt | sed -n ${SLURM_ARRAY_TASK_ID}p`

for MAG in `cat helper_files/list_of_all_MAGs.txt`;
do
    if [ ! -f metagenome_mappings/${STATION}_vs_${MAG}_scafstat.txt ]
    then
        bbmap.sh \
        in=reads_metagenome/${STATION}_microbial.reads1.1M.fq.gz \
        in2=reads_metagenome/${STATION}_microbial.reads2.1M.fq.gz \
        ref=contigs/${MAG} \
        scafstats=metagenome_mappings/${STATION}_vs_${MAG}_scafstat.txt \
        t=${SLURM_CPUS_PER_TASK} \
        nodisk \
        unpigz=t
    else
        echo "top20/metagenome_mappings/bincov/${STATION}_vs_${MAG}.bincov exists"
    fi
done
