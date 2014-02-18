#!/bin/bash
## Estimate gaussians ########
## Parameters #
# path to gaussian mixture executable gmm
gmx=gmm

# hard-coded parameters
err=1e-5       # relative change of log-likelihood to bail out from GM optimization
smooth=1e-5    # regularization of the GM covariance matrices
seed=456738    # base of random number seeds
maxtry=100     # number of attempts

# initialize parameters
read -p "Input file name ( nu mu dad format): " fin
read -p "Stride when reading input data: " ev
read -p "Number of Gaussian clusters: " gn
read -p "Output files prefix: " fo
##############

logbest=-3432943209023
for ((j=1; j<=$maxtry; j++)); do
  ((seed+=121))
  echo "Trial number $j."
  $gmx -i $fin -n $gn -seed $seed -o $fo.$j -ev $ev -err $err -s $smooth
  loglike=$( head -1 $fo.$j | awk '{print $5}' )
  echo "    loglike/nsample: $loglike"
  if [ $( echo " $loglike > $logbest " | bc ) -eq 1 ]; then
    echo "Found improved clusters"
    logbest=$loglike
    cp $fo.$j $fo.best
  fi
done
