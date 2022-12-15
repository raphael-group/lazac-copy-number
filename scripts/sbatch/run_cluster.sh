#!/bin/bash
#SBATCH -J snakemake-breaked
#SBATCH -o snakemake-breaked.out.log
#SBATCH -e snakemake-breaked.err.log
#SBATCH --nodes=1                # node count
#SBATCH --account=raphael                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time 0-03:00:00        # DAYS-HOURS:MINUTES:SECONDS

snakemake -p --profile ./snakemake/profiles/ionic
