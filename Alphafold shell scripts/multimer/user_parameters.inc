## Set user parameters for alphafold, common to MSA and prediction jobs

# Fasta file of protein data, passed to alphafold via --fasta_paths
# Will also be used without .fasta as output directory name
FASTA_PATHS=/ptmp/${USER}/fasta/seq.fasta

# Output directory to write the results to, passed via --output
OUTPUT_DIR=/ptmp/${USER}/alphafold_2.2.0_output

# Preset for alphafold, passed via --db_preset
PRESET="full_dbs"

# Number of predictions per model
PRED=5