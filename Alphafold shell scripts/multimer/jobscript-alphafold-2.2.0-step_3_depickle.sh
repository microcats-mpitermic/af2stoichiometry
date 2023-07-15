#!/bin/bash -l
#SBATCH -J AF2_PAE
##SBATCH -o ./out.%j
##SBATCH -e ./err.%j
##SBATCH -D ./
#SBATCH --nodes=1             # request a full node
#SBATCH --cpus-per-task=18    
#SBATCH --mem=60000
#SBATCH --time=00:30:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=userid@mpg.de #edit to your own email adress

set -e

module purge
module load gcc/10 impi/2021.2
module load anaconda/3/2021.11

# Important:
# Set the number of OMP threads *per process* to avoid overloading of the node!
# Set number of OMP threads to fit the number of available cpus, if applicable.
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# load the user specifications
source user_parameters.inc

# grab the output folder of AlphaFold GPU job, which is the truncated path of the input .fasta file
AF_FOLDER="${OUTPUT_DIR}/$(basename ${FASTA_PATHS} .fasta)/"

echo "Running metrics extraction on ${AF_FOLDER}"
echo "Input has ${PRED} predictions per model, running alphafold_depickle.py script"

# run the depickle script that extracts the run metrics from the .pkl files
# it will print relevant file information into the slurm output file while running
srun python3 /u/${USER}/python/alphafold_depickle.py \
		--input_dir="${AF_FOLDER}" \
		--num_model=5 \
		--run_count=${PRED}