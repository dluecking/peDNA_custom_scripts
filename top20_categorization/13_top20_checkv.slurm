#!/bin/bash -l

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=checkv_arr  # Job name
#SBATCH --output=out/checkv/%x_%A.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=out/checkv/%x_%A.%j.err # Stderr (%j expands to jobId)
#SBATCH --time=02:00:00   # walltime
#SBATCH --mem=10G   # memory per NODE
#SBATCH --cpus-per-task=8
#SBATCH --partition=CLUSTER
#SBATCH --array=11-116%15

MAG=`cat helper_files/top20_contigs_unique.txt | sed -n ${SLURM_ARRAY_TASK_ID}p`

# only run if we have this file contains more than one row, otherwise we dont 
# have any viruses detected and we dont have to run!
LINES=$(wc -l < top20/vs2/${MAG}/final-viral-score.tsv)

if [ ${LINES} -gt 1 ]
then    
	# also we dont run if already done!
	if [ ! -f top20/checkv/${MAG}/quality_summary.tsv ]
	then
	    source /home/dlueckin/bin/miniconda3/etc/profile.d/conda.sh
	    conda activate vs2
	   
		checkv end_to_end top20/vs2/${MAG}/final-viral-combined.fa \
		top20/checkv/${MAG} -t 8 -d /home/dlueckin/dbs/checkv-db-v1.0
	else
		echo "${MAG} already done!"
	fi	    
else
    echo "No virus found in ${MAG}, skipping..."
fi
