#!/bin/bash
for i in job.*; do
  energy=`grep "nal e" $i | xargs | cut -f6 -d" "` 
  point=`echo $i | sed s/job.// | sed s/^0*//`
  #point=$((point+1))
  abscissa=`sed "${point}q;d" abscissa.dat | xargs | cut -f2 -d" "` 
  echo $point $abscissa $energy
done
