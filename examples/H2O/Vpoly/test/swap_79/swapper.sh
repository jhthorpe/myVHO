#/bin/bash
i=7
j=9
sed -i "s/$i /a /g"  rota
sed -i "s/$j /$i /g" rota
sed -i "s/a /$j /g"  rota
sed -i "s/$i /a /g"  quadratic
sed -i "s/$j /$i /g" quadratic
sed -i "s/a /$j /g"  quadratic
sed -i "s/$i /a /g"  xdidq
sed -i "s/$j /$i /g" xdidq
sed -i "s/a /$j /g"  xdidq
sed -i "s/$i /a /g"  ydidq
sed -i "s/$j /$i /g" ydidq
sed -i "s/a /$j /g"  ydidq
sed -i "s/$i /a /g"  xcoriolis
sed -i "s/$j /$i /g" xcoriolis
sed -i "s/a /$j /g"  xcoriolis
sed -i "s/$i /a /g"  ycoriolis
sed -i "s/$j /$i /g" ycoriolis
sed -i "s/a /$j /g"  ycoriolis
