#!/bin/bash

cwd=$(pwd)

for file in $cwd/output/*.txt; do
	echo "running $file"
	python3 compute_output.py $file >> results.txt
done;