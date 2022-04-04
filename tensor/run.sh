#!/bin/bash

for i in {0..49};
do
  echo $i
  #python3 gtn_test.py BV $1 >> BV_G_$1.txt
  python3 test.py DJ $1 >> DJ_$1.txt
done;
