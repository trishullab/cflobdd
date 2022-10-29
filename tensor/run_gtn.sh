#!/bin/bash


#for j in {4..4};
#do
for i in {1..49};
do
  echo $i
  val=$((2 ** j))
  python3 gtn_test.py qft 16 $i >> qft_G_4.txt
done;
#done;
