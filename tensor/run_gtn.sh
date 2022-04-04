#!/bin/bash

for i in {0..49};
do
  echo $i
  python3 gtn_test.py DJ $1 >> DJ_G_$1.txt
done;
