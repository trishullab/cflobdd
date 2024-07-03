#! /bin/bash

numRuns=5

for (( c=0; c<$numRuns; c++ ))
do
	echo "Running testSyn1 15"
	timeout -v 15m ./cflobdd testSyn1 15 >> results_syn/wcflobdd_testSyn1_15.txt
	echo "Running testSyn1 16"
	timeout -v 15m ./cflobdd testSyn1 16 >> results_syn/wcflobdd_testSyn1_16.txt
	echo "Running testSyn1 17"
	timeout -v 15m ./cflobdd testSyn1 17 >> results_syn/wcflobdd_testSyn1_17.txt
	echo "Running testSyn1 18"
	timeout -v 15m ./cflobdd testSyn1 18 >> results_syn/wcflobdd_testSyn1_18.txt
	echo "Running testSyn1 19"
	timeout -v 15m ./cflobdd testSyn1 19 >> results_syn/wcflobdd_testSyn1_19.txt
	echo "Running testSyn1 20"
	timeout -v 15m ./cflobdd testSyn1 20 >> results_syn/wcflobdd_testSyn1_20.txt
done

for (( c=0; c<$numRuns; c++ ))
do
	echo "Running testSyn2 15"
	timeout -v 15m ./cflobdd testSyn2 15 >> results_syn/wcflobdd_testSyn2_15.txt
	echo "Running testSyn2 16"
	timeout -v 15m ./cflobdd testSyn2 16 >> results_syn/wcflobdd_testSyn2_16.txt
	echo "Running testSyn2 17"
	timeout -v 15m ./cflobdd testSyn2 17 >> results_syn/wcflobdd_testSyn2_17.txt
	echo "Running testSyn2 18"
	timeout -v 15m ./cflobdd testSyn2 18 >> results_syn/wcflobdd_testSyn2_18.txt
	echo "Running testSyn2 19"
	timeout -v 15m ./cflobdd testSyn2 19 >> results_syn/wcflobdd_testSyn2_19.txt
	echo "Running testSyn2 20"
	timeout -v 15m ./cflobdd testSyn2 20 >> results_syn/wcflobdd_testSyn2_20.txt
done

for (( c=0; c<$numRuns; c++ ))
do
	echo "Running testSyn3 15"
	timeout -v 15m ./cflobdd testSyn3 15 >> results_syn/wcflobdd_testSyn3_15.txt
	echo "Running testSyn3 16"
	timeout -v 15m ./cflobdd testSyn3 16 >> results_syn/wcflobdd_testSyn3_16.txt
	echo "Running testSyn3 17"
	timeout -v 15m ./cflobdd testSyn3 17 >> results_syn/wcflobdd_testSyn3_17.txt
	echo "Running testSyn3 18"
	timeout -v 15m ./cflobdd testSyn3 18 >> results_syn/wcflobdd_testSyn3_18.txt
	echo "Running testSyn3 19"
	timeout -v 15m ./cflobdd testSyn3 19 >> results_syn/wcflobdd_testSyn3_19.txt
	echo "Running testSyn3 20"
	timeout -v 15m ./cflobdd testSyn3 20 >> results_syn/wcflobdd_testSyn3_20.txt
done

for (( c=0; c<$numRuns; c++ ))
do
	echo "Running testSyn4 15"
	timeout -v 15m ./cflobdd testSyn4 15 >> results_syn/wcflobdd_testSyn4_15.txt
	echo "Running testSyn4 16"
	timeout -v 15m ./cflobdd testSyn4 16 >> results_syn/wcflobdd_testSyn4_16.txt
	echo "Running testSyn4 17"
	timeout -v 15m ./cflobdd testSyn4 17 >> results_syn/wcflobdd_testSyn4_17.txt
	echo "Running testSyn4 18"
	timeout -v 15m ./cflobdd testSyn4 18 >> results_syn/wcflobdd_testSyn4_18.txt
	echo "Running testSyn4 19"
	timeout -v 15m ./cflobdd testSyn4 19 >> results_syn/wcflobdd_testSyn4_19.txt
	echo "Running testSyn4 20"
	timeout -v 15m ./cflobdd testSyn4 20 >> results_syn/wcflobdd_testSyn4_20.txt
done

for (( c=0; c<$numRuns; c++ ))
do
	echo "Running testSyn5 15"
	timeout -v 15m ./cflobdd testSyn5 15 >> results_syn/wcflobdd_testSyn5_15.txt
	echo "Running testSyn5 16"
	timeout -v 15m ./cflobdd testSyn5 16 >> results_syn/wcflobdd_testSyn5_16.txt
	echo "Running testSyn5 17"
	timeout -v 15m ./cflobdd testSyn5 17 >> results_syn/wcflobdd_testSyn5_17.txt
	echo "Running testSyn5 18"
	timeout -v 15m ./cflobdd testSyn5 18 >> results_syn/wcflobdd_testSyn5_18.txt
	echo "Running testSyn5 19"
	timeout -v 15m ./cflobdd testSyn5 19 >> results_syn/wcflobdd_testSyn5_19.txt
	echo "Running testSyn5 20"
	timeout -v 15m ./cflobdd testSyn5 20 >> results_syn/wcflobdd_testSyn5_20.txt
done




