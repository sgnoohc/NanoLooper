#!/bin/bash


rm -f .jobs.txt
rm -f .nevents_*.txt

index=0
for i in $(ls /home/users/joytzphysics/VBSZZH_C2V_3/*.root);
do

    index=$((index+1))

    echo "getentriesroot $i Events > .nevents_${index}.txt 2>&1" >> .jobs.txt

done

xargs.sh .jobs.txt

cat .nevents_*.txt | grep -v "Warning" | awk '{sum += $1} END {print sum}'
