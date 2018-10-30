#!/bin/sh

#SBATCH --partition=short
#SBATCH --mem 10000

module load BEDTools
coverageBed -sorted -b $(echo $1) -a $(echo $2) > $(echo $3)
