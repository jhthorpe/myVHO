#/bin/bash
for j in `cat dirlist` ; do
  cd $j 
  ./doit
  nl=`cat ../QNlist | wc -l`
  echo "VHO $j"
  echo "---------- ------------------"
  for i in `seq 1 $nl` ; do 
    a=`head -$i ../QNlist | tail -1`
    b="$(echo -e "${a}" | tr -d '[:space:]')"
    rm out_vho_$j_$b
    echo $a | xvho > out_vho_$j_$b 
    echo $a `grep "VPT4" out_vho_$j_$b | xargs | cut -f 4 -d " "` 
  done
  echo "Guinea $j"
  xguinea < ../runner > out_guinea_$j
  sed -n -e '/------------------/,$p' out_guinea_$j
  cd ..
done
