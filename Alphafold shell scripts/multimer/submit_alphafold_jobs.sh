#!/bin/bash

set -e

JOBID1=$(sbatch --parsable jobscript-alphafold-2.2.0-step_1-msa.sh)
JOBID2=$(sbatch --parsable --dependency=afterok:${JOBID1} --deadline=now+2weeks jobscript-alphafold-2.2.0-step_2-prediction.sh)
JOBID3=$(sbatch --parsable --dependency=afterok:${JOBID2} --deadline=now+2weeks jobscript-alphafold-2.2.0-step_3_depickle.sh)

echo "Submitted jobs"
echo "    ${JOBID1} (MSA on CPU)"
echo "    ${JOBID2} (Structure prediction on GPU)"
echo "    ${JOBID3} (AF2 metrics of run on CPU)"

