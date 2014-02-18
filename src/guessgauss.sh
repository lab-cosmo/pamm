#!/bin/bash
## Estimate gaussians ########
## Parameters #
gmx=ndhistogram
gn=8
seed=456738
ev=10
err=0.000001
smooth=0.00000001
cpaletta=-5:12
fin="oxo.dat"
fo="gassuot.dat"
##############

i=0
ir=0
tstr=-3432943209023
while :
do
  i=$(( $i + 1 ))
  echo "Trial number $i."
  gmm -i $fin -n $gn -seed $seed -o $fo.$i -ev $ev -err $err -s $smooth
  tst = head -1 $fo.$i  | awk '{print $9}'
  echo "    loglike/nsample: $tst"
  if [ $tst -ge $tstr ]; then
    tstr=$tst
    rm $fo.$ir
    ir=$i
  else
    rm $fo.$i
  fi
done

#gmm -i $fxy -n $gn -seed $seed -o $fo -ev $ev -err $err -s $smooth
