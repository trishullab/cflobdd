#!/bin/bash


for j in {2..4};
do
for i in {0..49};
do
  echo $i
  val=$((2 ** j))
  python3 gtn_test.py qft $val $i >> qft_G_$j.txt
done;
done;
