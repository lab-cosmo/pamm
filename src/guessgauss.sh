#!/bin/bash
## Estimate gaussians ########
## Parameters #
# path to gaussian mixture executable gmm
gmx=gmm

# hard-coded parameters
err=1e-5       # relative change of log-likelihood to bail out from GM optimization
smooth=1e-5    # regularization of the GM covariance matrices
seed=456738    # base of random number seeds
#maxtry=250     # number of attempts
v=-1.1
w=2.9
dAD=2.9

# initialize parameters
read -p "Input file name ( nu mu dad format): " fin
read -p "Stride when reading input data: " ev
read -p "Number of Gaussian clusters: " gn
read -p "Number of random attempts: " maxtry
read -p "Output files prefix: " fo
##############

logbest=-3432943209023
for ((j=1; j<=$maxtry; j++)); do
  ((seed+=121))
  echo "Trial number $j."
  $gmx -i $fin -n $gn -seed $seed -o $fo.$j -ev $ev -err $err -s $smooth -rif $v,$w,$dAD 
  loglike=$( head -1 $fo.$j | awk '{print $5}' )
  echo "    loglike/nsample: $loglike"
  if [ $( echo " $loglike > $logbest " | bc ) -eq 1 ]; then
    echo "Found improved clusters"
    logbest=$loglike
    cp $fo.$j $fo.best
  fi
done
((seed+=1121))
echo "Maxmin trial."
$gmx -i $fin -n $gn -seed $seed -o $fo.maxmin -maxmin -ev $ev -err $err -s $smooth -rif $v,$w,$dAD
loglike=$( head -1 $fo.maxmin | awk '{print $5}' )
echo "    loglike/nsample: $loglike"
if [ $( echo " $loglike > $logbest " | bc ) -eq 1 ]; then
  echo "Found improved clusters"
  logbest=$loglike
  cp $fo.maxmin $fo.best
fi

echo "Redo startimg from the best trial, but using error."
$gmx -i $fin -n $gn -gf $fo.best -o $fo.init -ev $ev -err 0.0000000001 -s $smooth -rif $v,$w,$dAD
loglike=$( head -1 $fo.init | awk '{print $5}' )
echo "    loglike/nsample: $loglike"
if [ $( echo " $loglike > $logbest " | bc ) -eq 1 ]; then
  echo "Last attempt improved the LogLike!"
  logbest=$loglike
  cp $fo.init $fo.best
fi
