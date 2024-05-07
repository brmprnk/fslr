#! /bin/sh
#SBATCH --partition=general --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=40000
#SBATCH --time=03:59:00
#SBATCH --job-name=fslr_brca
#SBATCH --mail-user=ibpronk
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_%j.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=slurm_%j.err # Set name of error log. %j is the Slurm jobId

ml use /opt/insy/modulefiles;

python src/fslr.py -c configs/brca.yaml