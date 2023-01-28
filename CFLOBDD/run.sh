#!/bin/bash


for (( i = 4; i < 20; i++ )); do
    echo " running testGHZAlgo_W for $i bits"
    timeout -v 15m ./a.out testGHZAlgo_W $i >> ghz.txt
    echo " running testBVAlgo_W for $i bits"
    timeout -v 15m ./a.out testBVAlgo_W $i >> bv.txt
    echo " running testDJAlgo_W for $i bits"
    timeout -v 15m ./a.out testDJAlgo_W $i >> dj.txt
done;

