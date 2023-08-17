#!/bin/bash -l

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=true_virome_mapping_76   # Job name
#SBATCH --output=out/%x.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=out/%x.%j.err # Stderr (%j expands to jobId)
#SBATCH --time=1-00:00:00   # walltime
#SBATCH --mem=50G   # memory per NODE
#SBATCH --cpus-per-task=8
#SBATCH --array=1-9
#SBATCH --partition=CLUSTER


# source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
# conda activate bbmap

# station=`cat helper_files/list_of_stations.txt | sed -n ${SLURM_ARRAY_TASK_ID}p`

# bbmap.sh \
# in=reads/${station}.true.virome.reads1.1M.fq.gz \
# in2=reads/${station}.true.virome.reads2.1M.fq.gz \
# ref=contigs/combined_MAGs.txt \
# outm=true_virome_mappings/${station}_vs_combined_MAGs_outm.reads1.fq \ 
# outm2=true_virome_mappings/${station}_vs_combined_MAGs_outm.reads2.fq \ 
# t=8 \
# overwrite=false

# id=95
source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
conda activate bbmap

station=`cat helper_files/list_of_stations.txt | sed -n ${SLURM_ARRAY_TASK_ID}p`

bbmap.sh \
in=reads/${station}.true.virome.reads1.1M.fq.gz \
in2=reads/${station}.true.virome.reads2.1M.fq.gz \
ref=contigs/combined_MAGs.txt \
outm=true_virome_mappings/${station}_vs_combined_MAGs_outm_76.reads1.fq \
outm2=true_virome_mappings/${station}_vs_combined_MAGs_outm_76.reads2.fq \
t=8 \
overwrite=false


