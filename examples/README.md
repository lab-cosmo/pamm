PAMM library examples
=====================

This folder contains a few examples (currently one!) of applications
of the PAMM method. They should give an idea of the workflow for running
a PAMM-based analysis.

HBPAMM -- Hydrogen bonding in quantum BLYP water
------------------------------------------------

This example uses a data file `h2o-blyp-piglet.xyz` that contains a few
snapshots from a simulation of liquid water at 300K, using a DFT-BLYP
first-principles evaluation of the forces and the PIGLET method to model
the quantum nature of the nuclei.  Hence, both electronic and nuclear
quantum effects are considered, and contribute to define the details
of the hydrogen bond network.

The PAMM analysis consists of three steps.

1. **Evaluation of the descriptors of the hydrogen bond.**
    First, one needs to analyze the trajectory to compute the values of 
    the proton transfer coordinate nu, the symmetric stretch coordinate mu
    and the donor-acceptor distance r for each O-H...O triplet in each
    snapshot. This is done using
 
    ```bash
     ../bin/hbpamm  -td O -th H -ta O -ct 4.5 -w < h2o-blyp-piglet.xyz > h2o.nmr
    ```

    `-td` specifies the atom type(s) that are to be tested as donors, 
    `-th` specifies the label that identifies hydrogen atoms, and `-ta`
    specifies the atom type(s) that are to be tested as acceptors. 
    `-ct` discards triplets with a value of mu greater than the 
    cutoff value specified, and `-w` asks to compute a re-normalization
    weight that accounts for the trivial phase space volume factor.
    The program tries to infer the cell parameters from the header of the
    `xyz` file, expecting the format `# CELL x y z`.

2. **Optimization of the PAMM model.**
    Next, one can run the PAMM analysis on the set of HB parameters:
 
    ```bash
    ../bin/pamm -d 3 -w -o h2o -ngrid 2000 -nms 100 < h2o.nmr
    ```

    Note that most parameters in PAMM are selected automatically to 
    give reasonable defaults. Here we specify the number of grid points
    manually, since this is a very short data file and the default would
    yield a very sparse grid. 

    `pamm` generates two files, `h2o.grid` that cointains the grid points,
    that have been selected, together with the value of the KDE of the 
    probability density and with the cluster index for each point. 
    `h2o.pamm` contains the list of clusters that have been identified,
    represented by Gaussians with the given mean and covariance.

3. **Post-processing the trajectory.**
    The set of Gaussian clusters obtained by PAMM can then be used to 
    search the trajectory file for HB configurations. 

    `hbpamm` accumulates the number of such HBs that each atom is involved into, 
    and prints an `xyz` formatted file in which the first column contains sH 
    (the number of HBs in which the atom takes the role of the hydrogen), the 
    second sD (the number of HB in which the atom is the donor) and the third
    sA (number of acceptor HBs). Options are the same as for the first
    step, plus a specification of the cluster file:

    ```bash
    ../bin/hbpamm  -td O -th H -ta O -ct 4.5 -w -gf h2o.pamm < h2o-blyp-piglet.xyz > h2o.hda
    ```

    `hbpamm` automatically selects the first cluster as the one that 
    represents the HB, but one can choose another using the `-ghb` option.
  
 
The data in `h2o.hda` can be post-processed further, to obtain for instance
an histogram of the probability of having hydrogen atoms involved in 
more than one HB, or the joint probability of sA and sD. To this aim,
one can use the histogram codes that is part of the [toolbox](http://github.com/epfl-cosmo/toolbox) 
post-processing and utilities suite. 

```bash
grep H h2o.hda | awk '{print $2}' | histogram -xi 0 -xf 3 -n 300 -t 0.05 -whard > h2o.hb
grep O h2o.hda | awk '{print $3,$4}' | ndhistogram -d 2 -xi 0,0 -xf 4,4 \
           -n 200,200 -t 0.05,0.05 -whard -g -adaptive 1.0 > h2o.sasd
```

`hbpamm` also make it possible to compute *all* the HB counts between pairs of 
acceptors and donors, by specifying the `-sad` option. This generates a huge file,
that contain on each column the total number of HBs at any given time for each
A/D pair. 
This can be used to evaluate the hydrogen-bond dynamical relaxation function
(*Luzar & Chandler, Nature 1996*). Again, using the `autocorr` program from
the [toolbox](http://github.com/epfl-cosmo/toolbox) library one can compute this by

```bash
../bin/hbpamm  -td O -th H -ta O -ct 4.5 -w -gf h2o.pamm -sad < h2o-blyp-piglet.xyz > h2o.sad
for i in `seq 1 $(head -n 1 h2o.sad | wc -w)`; do awk -v i=$i '{print $i}' h2o.sad; done | \
        autocorr -maxlag 150 -runlength $(wc -l h2o.sad) > h2o.hh
```


