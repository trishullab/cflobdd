#!/bin/bash

for j in {2..4};
do
for i in {0..49};
do
  echo $i
  #python3 gtn_test.py BV $1 >> BV_G_$1.txt
  val=$((2 ** j))
  python3 test.py qft $val  $i >> qft_$j.txt
done;
done;
