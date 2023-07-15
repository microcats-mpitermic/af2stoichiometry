## Set user parameters for alphafold, common to MSA and prediction jobs

# Comma-separated list of protein data, passed to alphafold via --fasta_paths
FASTA_PATHS=/ptmp/${USER}/fasta/seq.fasta

# Output directory to write the results to, passed via --output
OUTPUT_DIR=/ptmp/${USER}/alphafold_2.2.0_output

# Preset for alphafold, passed via --db_preset
PRESET="full_dbs"

