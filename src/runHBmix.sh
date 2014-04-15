#!/bin/bash

# makes histograms of nu, mu, dad starting from an input file containing a list in the format
# nu mu dad
# call as run.sh input prefix cutoff acceptor donor hydrogen (h)
# 
#   input    : trajectory to process
#   prefix   : outputfiles prefix
#   cutoff   : cutoff in w
#   acceptor : acceptor type
#   donor    : donor type
#   hydrogen : hydrogen type
#   h        : add h to compute the histograms

# program to call #
hban=hbanalysis
getm=getmodel
ndhisto=ndhistogram
histo=histogram

# hard-coded parameters
n=550
((n2=2*n+1))
b=0.012
#0.0142857

# input parameters
# coordinates file
input=$1
# output name prefix
prefix=$2
# cutoff
xf=$3
# acceptor,donor and hydrogen types
acc=$4
don=$5
hyd=$6

if [ $7 = 'h' ] 
 then
    # dir for the histograms
    mkdir histo
fi

# dir for the gaussians
mkdir clusters
# dir for the FE
mkdir FE

echo " First step : calculate v,w,R and the weight from the trajectory, clusterize them and get the gaussians"

if [ $7 = 'h' ] 
 then
    # calculate v,w,R and weights and save to a temporary file
    cat $input | $hban -ta $acc -th $hyd -td $don -ct $xf -w -P > tmphisto.dat
    
    cat tmphisto.dat | $getm -w -d 3 -o clusters/$prefix -v &
    
    # Get histograms
    echo " Also computing the histograms"
    # weighted
    # makes 1D histograms with normalization factor readed from the file
    awk '!/#/{print $1,$4}' tmphisto.dat | $histo -xi -$xf -xf $xf -n $n2 -t $b -whard -w > histo/$prefix.h1 &
    awk '!/#/{print $2,$4}' tmphisto.dat | $histo -xi 0 -xf $xf -n $n -t $b -whard -w > histo/$prefix.h2 &
    awk '!/#/{print $3,$4}' tmphisto.dat | $histo -xi 0 -xf $xf -n $n -t $b -whard -w > histo/$prefix.h3 &
    # makes 2D histograms with normalization factor readed from the file
    awk '!/#/{print $1,$2,$4}' tmphisto.dat | $ndhisto -d 2 -xi -$xf,0 -xf $xf,$xf -n $n2,$n -b $b,$b -g -whard -w > histo/$prefix.h12 &
    awk '!/#/{print $2,$3,$4}' tmphisto.dat | $ndhisto -d 2 -xi 0,0 -xf $xf,$xf -n $n,$n -b $b,$b -g -whard -w > histo/$prefix.h23 &
    awk '!/#/{print $1,$3,$4}' tmphisto.dat | $ndhisto -d 2 -xi -$xf,0 -xf $xf,$xf -n $n2,$n -b $b,$b -g -whard -w > histo/$prefix.h13 &

    # not weighted
    # makes 1D histograms
    awk '!/#/{print $1}' tmphisto.dat | $histo -xi -$xf -xf $xf -n $n2 -t $b -whard > histo/$prefix-nw.h1 &
    awk '!/#/{print $2}' tmphisto.dat | $histo -xi 0 -xf $xf -n $n -t $b -whard > histo/$prefix-nw.h2 &
    awk '!/#/{print $3}' tmphisto.dat | $histo -xi 0 -xf $xf -n $n -t $b -whard > histo/$prefix-nw.h3 &
    # makes 2D histograms
    awk '!/#/{print $1,$2}' tmphisto.dat | $ndhisto -d 2 -xi -$xf,0 -xf $xf,$xf -n $n2,$n -b $b,$b -g -whard > histo/$prefix-nw.h12 & 
    awk '!/#/{print $2,$3}' tmphisto.dat | $ndhisto -d 2 -xi 0,0 -xf $xf,$xf -n $n,$n -b $b,$b -g -whard > histo/$prefix-nw.h23
    awk '!/#/{print $1,$3}' tmphisto.dat | $ndhisto -d 2 -xi -$xf,0 -xf $xf,$xf -n $n2,$n -b $b,$b -g -whard > histo/$prefix-nw.h13
    
    # remove the temporary file
    rm tmphisto.dat
 else
 
    
    # decomment to choose the number of points
    # cat $input | $hban -ta OX -th HT -td OT -ct $3 -w -P | $getm -w -d 3 -nminmax 2000 -nsamples 50000 -o clusters/testooz -v

    cat $input | $hban -ta $acc -th $hyd -td $don -ct $xf -w -P | $getm -w -d 3 -o clusters/$prefix -v
fi
# First step

# Second step
echo " Second step : processing of the trajectory to find HB patterns"

# Do the hbmixture
cat $input | $hban -ta $acc -th $hyd -td $don -ct $xf -w -gf clusters/$prefix.gauss -ghb 1 -a 1 > FE/$prefix.xyz

# Get FEs

# makes 1D histograms - marginal probabilities
# Sh
awk -v class="${hyd}" '$0 ~ class {print $5}' FE/$prefix.xyz | $histo -xi 0 -xf 3.5 -n $n -t $b -whard > FE/$prefix.hh &
# Sd
awk -v class="${don}" '$0 ~ class {print $6}' FE/$prefix.xyz | $histo -xi 0 -xf 3.5 -n $n -t $b -whard > FE/$prefix.hd &
# Sa
awk -v class="${acc}" '$0 ~ class {print $7}' FE/$prefix.xyz | $histo -xi 0 -xf 3.5 -n $n -t $b -whard > FE/$prefix.ha &


# makes 2D histograms - jont probability
# Sd | Sa for don
awk -v class="${don}" '$0 ~ class {print $6, $7}' FE/$prefix.xyz | $ndhisto -d 2 -xi 0,0 -xf 3.5,3.5 -n $n,$n -b $b,$b -g -whard > FE/$prefix.hdd &
# Sd | Sa for acc 
awk -v class="${acc}" '$0 ~ class {print $6, $7}' FE/$prefix.xyz | $ndhisto -d 2 -xi 0,0 -xf 3.5,3.5 -n $n,$n -b $b,$b -g -whard > FE/$prefix.haa &
