#!/bin/bash -l

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=detailed_mapping   # Job name
#SBATCH --output=out/detailed_mapping/%x.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=out/detailed_mapping/%x.%j.err # Stderr (%j expands to jobId)
#SBATCH --time=01:00:00   # walltime
#SBATCH --mem=50G   # memory per NODE
#SBATCH --cpus-per-task=16
#SBATCH --partition=CLUSTER
#SBATCH --array=2-181%20 # needs to be 2 to 181, since first row is header

source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
conda activate bbmap

MAG=`cat top20/top20_df_mapping.tsv | cut -d ',' -f 1 | sed -n ${SLURM_ARRAY_TASK_ID}p`         # column 1 = name
STATION=`cat top20/top20_df_mapping.tsv | cut -d ',' -f 10 | sed -n ${SLURM_ARRAY_TASK_ID}p`    # column 10 = corresponding station

if [ ! -f top20/detailed_MAG_mappings/bincov/${STATION}_vs_${MAG}.bincov ]
then
    bbmap.sh \
    in=reads/${STATION}.all-input.reads1.fq.gz \
    in2=reads/${STATION}.all-input.reads2.fq.gz \
    ref=contigs/${MAG} \
    rpkm=top20/detailed_MAG_mappings/rpkm/${STATION}_vs_${MAG}.rpkm \
    basecov=top20/detailed_MAG_mappings/basecov/${STATION}_vs_${MAG}.basecov \
    bincov=top20/detailed_MAG_mappings/bincov/${STATION}_vs_${MAG}.bincov \
    t=${SLURM_CPUS_PER_TASK} \
    nodisk \
    minid=0.9 \
    unpigz=t
else
    echo "top20/detailed_MAG_mappings/bincov/${STATION}_vs_${MAG}.bincov exists"
fi