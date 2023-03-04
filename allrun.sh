#!/bin/sh

file="list_of_adsorbates.txt"
run="run.sh"
IFS=$'\n'

for line in `cat ${file}`
do
    ADSB=`echo ${line} | awk '{printf $1}'`
    ROT1=`echo ${line} | awk '{printf $2}'`
    ROT2=`echo ${line} | awk '{printf $3}'`
    pjsub $run -x "INP1=$ADSB,INP2=$ROT1,INP3=$ROT2"
    sleep 0.5
done
