#!/bin/bash
for i in job.*; do
  energy=`grep "nal e" $i | xargs | cut -f6 -d" "` 
  point=`echo $i | sed s/job.// | sed s/^0*//`
  #point=$((point+1))
  echo $point $energy
done
