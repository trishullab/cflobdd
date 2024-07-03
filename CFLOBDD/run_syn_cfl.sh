#! /bin/bash

numRuns=5

for (( c=0; c<$numRuns; c++ ))
do
	echo "Running testSyn1 15"
	timeout -v 15m ./cflobdd testSyn1_CFL 15 >> results_syn/testSyn1_CFL_15.txt
	echo "Running testSyn1 16"
	timeout -v 15m ./cflobdd testSyn1_CFL 16 >> results_syn/testSyn1_CFL_16.txt
	echo "Running testSyn1 17"
	timeout -v 15m ./cflobdd testSyn1_CFL 17 >> results_syn/testSyn1_CFL_17.txt
	echo "Running testSyn1 18"
	timeout -v 15m ./cflobdd testSyn1_CFL 18 >> results_syn/testSyn1_CFL_18.txt
	echo "Running testSyn1 19"
	timeout -v 15m ./cflobdd testSyn1_CFL 19 >> results_syn/testSyn1_CFL_19.txt
	echo "Running testSyn1 20"
	timeout -v 15m ./cflobdd testSyn1_CFL 20 >> results_syn/testSyn1_CFL_20.txt
done



for (( c=0; c<$numRuns; c++ ))
do
	echo "Running testSyn2 15"
	timeout -v 15m ./cflobdd testSyn2_CFL 15 >> results_syn/testSyn2_CFL_15.txt
	echo "Running testSyn2 16"
	timeout -v 15m ./cflobdd testSyn2_CFL 16 >> results_syn/testSyn2_CFL_16.txt
	echo "Running testSyn2 17"
	timeout -v 15m ./cflobdd testSyn2_CFL 17 >> results_syn/testSyn2_CFL_17.txt
	echo "Running testSyn2 18"
	timeout -v 15m ./cflobdd testSyn2_CFL 18 >> results_syn/testSyn2_CFL_18.txt
	echo "Running testSyn2 19"
	timeout -v 15m ./cflobdd testSyn2_CFL 19 >> results_syn/testSyn2_CFL_19.txt
	echo "Running testSyn2 20"
	timeout -v 15m ./cflobdd testSyn2_CFL 20 >> results_syn/testSyn2_CFL_20.txt
done

for (( c=0; c<$numRuns; c++ ))
do
	echo "Running testSyn3 15"
	timeout -v 15m ./cflobdd testSyn3_CFL 15 >> results_syn/testSyn3_CFL_15.txt
	echo "Running testSyn3 16"
	timeout -v 15m ./cflobdd testSyn3_CFL 16 >> results_syn/testSyn3_CFL_16.txt
	echo "Running testSyn3 17"
	timeout -v 15m ./cflobdd testSyn3_CFL 17 >> results_syn/testSyn3_CFL_17.txt
	echo "Running testSyn3 18"
	timeout -v 15m ./cflobdd testSyn3_CFL 18 >> results_syn/testSyn3_CFL_18.txt
	echo "Running testSyn3 19"
	timeout -v 15m ./cflobdd testSyn3_CFL 19 >> results_syn/testSyn3_CFL_19.txt
	echo "Running testSyn3 20"
	timeout -v 15m ./cflobdd testSyn3_CFL 20 >> results_syn/testSyn3_CFL_20.txt
done

for (( c=0; c<$numRuns; c++ ))
do
	echo "Running testSyn4 15"
	timeout -v 15m ./cflobdd testSyn4_CFL 15 >> results_syn/testSyn4_CFL_15.txt
	echo "Running testSyn4 16"
	timeout -v 15m ./cflobdd testSyn4_CFL 16 >> results_syn/testSyn4_CFL_16.txt
	echo "Running testSyn4 17"
	timeout -v 15m ./cflobdd testSyn4_CFL 17 >> results_syn/testSyn4_CFL_17.txt
	echo "Running testSyn4 18"
	timeout -v 15m ./cflobdd testSyn4_CFL 18 >> results_syn/testSyn4_CFL_18.txt
	echo "Running testSyn4 19"
	timeout -v 15m ./cflobdd testSyn4_CFL 19 >> results_syn/testSyn4_CFL_19.txt
	echo "Running testSyn4 20"
	timeout -v 15m ./cflobdd testSyn4_CFL 20 >> results_syn/testSyn4_CFL_20.txt
done

for (( c=0; c<$numRuns; c++ ))
do
	echo "Running testSyn5 15"
	timeout -v 15m ./cflobdd testSyn5_CFL 15 >> results_syn/testSyn5_CFL_15.txt
	echo "Running testSyn5 16"
	timeout -v 15m ./cflobdd testSyn5_CFL 16 >> results_syn/testSyn5_CFL_16.txt
	echo "Running testSyn5 17"
	timeout -v 15m ./cflobdd testSyn5_CFL 17 >> results_syn/testSyn5_CFL_17.txt
	echo "Running testSyn5 18"
	timeout -v 15m ./cflobdd testSyn5_CFL 18 >> results_syn/testSyn5_CFL_18.txt
	echo "Running testSyn5 19"
	timeout -v 15m ./cflobdd testSyn5_CFL 19 >> results_syn/testSyn5_CFL_19.txt
	echo "Running testSyn5 20"
	timeout -v 15m ./cflobdd testSyn5_CFL 20 >> results_syn/testSyn5_CFL_20.txt
done





