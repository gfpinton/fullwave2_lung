#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=3g
#SBATCH --time=00:10:00
#SBATCH --job-name=pz_1


module add matlab/2019b
matlab -nodesktop < lr.m
