#!/bin/bash

algos=('grover' 'DJ' 'BV' 'GHZ' 'simons')
#algos=('DJ' 'BV' 'GHZ')
reorder=('disable' 'enable')

for (( j = 0; j < 50; j++)); do
	for (( i = 5; i < 8; i++ )); do
		for algo in ${algos[@]};do
			for is_reorder in ${reorder[@]};do
				echo " run_num: $j ./big_test $algo $i $is_reorder >> output/out_${algo}_${i}_${is_reorder}.txt"
				timeout -v 15m ./big_test $algo $i $is_reorder >> output/out_${algo}_${i}_${is_reorder}.txt
				sleep 1
			done;
		done;
	done;
done;