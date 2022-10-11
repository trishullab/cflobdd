#!/bin/bash

<<<<<<< Updated upstream
#algos=('grover' 'DJ' 'BV' 'GHZ' 'simons')
algos=('DJ' 'BV' 'GHZ')
reorder=('disable')

for (( j = 0; j < 50; j++)); do
	for (( i = 4; i < 9; i++ )); do
=======
algos=('grover' 'DJ' 'BV' 'GHZ') 
#algos=('DJ' 'BV' 'GHZ')
reorder=('disable')

for (( j = 0; j < 50; j++)); do
	for (( i = 3; i < 12; i++ )); do
>>>>>>> Stashed changes
		for algo in ${algos[@]};do
			for is_reorder in ${reorder[@]};do
				val=$((j * 5));
				echo " run_num: $j ./big_test $algo $i $is_reorder >> output_nodes/out_${algo}_${i}_${is_reorder}.txt"
				timeout -v 15m ./big_test $algo $i $is_reorder $val >> output_nodes/out_${algo}_${i}_${is_reorder}.txt
				sleep 1
			done;
		done;
	done;
done;

<<<<<<< Updated upstream
for (( i = 4; i < 5; i++ )); do
	for (( j = 0; j < 50; j++)); do
		val=$((j * 5));
		echo " run_num: $j ./big_test grover $i disable $val" 
		timeout -v 90m ./big_test grover $i disable $val >> output_nodes/out_grover_${i}_disable.txt
		sleep 1
	done;
done;

for (( i = 4; i < 6; i++ )); do
	for (( j = 0; j < 50; j++)); do
		val=$((j * 5));
		echo " run_num: $j ./big_test simons $i disable $val" 
		timeout -v 90m ./big_test simons $i disable $val >> output_nodes/out_simons_${i}_disable.txt
		sleep 1
	done;
done;


=======
for (( j = 0; j < 50; j++)); do
    echo "run_num: $j ./big_test simons 4 disable >> output/out_simons_4_disable.txt"
    timeout -v 90m ./big_test simons 4 disable >> output/out_simons_4_disable.txt
    sleep 1
done;

for (( j = 0; j < 50; j++)); do
    echo "run_num: $j ./big_test simons 5 disable >> output/out_simons_5_disable.txt"
    timeout -v 90m ./big_test simons 5 disable >> output/out_simons_5_disable.txt
    sleep 1
done;

>>>>>>> Stashed changes
