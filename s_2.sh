#!/bin/bash

for i in {1..128}
do	
	cd -- "/pine/scr/o/l/oleksii/syn_ntdX_128/txrx_"$i;
	pwd;
	sbatch -t 3:00:00 -N 1 -n 1 --mem=1g --wrap="./fullwave2_try6_nln_relaxing_pzero"
	cd ..;
done