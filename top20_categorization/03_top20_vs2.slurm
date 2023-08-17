#!/bin/bash -l

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=vs2_arr  # Job name
#SBATCH --output=out/vs2/%x_%A.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=out/vs2/%x_%A.%j.err # Stderr (%j expands to jobId)
#SBATCH --time=2-00:00:00   # walltime
#SBATCH --mem=40G   # memory per NODE
#SBATCH --cpus-per-task=8
#SBATCH --partition=CLUSTER
#SBATCH --array=11-116

MAG=`cat helper_files/top20_contigs_unique.txt | sed -n ${SLURM_ARRAY_TASK_ID}p`

if [ ! -f top20/vs2/${MAG}/final-viral-score.tsv ]
then
    
    source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
    conda activate vs2
   
    virsorter run --keep-original-seq -i contigs/${MAG} \
    	-w top20/vs2/${MAG} --min-score 0.5 -j 8 all --db-dir ~/dbs/vs2/db/
else
    echo "Already done!"
fi
