#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8g
#SBATCH --time=00:20:00
#SBATCH --job-name=pz_3

module add matlab/2019b
matlab -nodesktop < bf.m
