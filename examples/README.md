# PAMM library examples

This folder contains a few examples (currently one!) of applications
of the PAMM method. They should give an idea of the workflow for running
a PAMM-based analysis.

## HBPAMM -- Hydrogen bonding in quantum BLYP water

This example uses a data file `h2o-blyp-piglet.xyz` that contains a few
snapshots from a simulation of liquid water at 300K, using a DFT-BLYP
first-principles evaluation of the forces and the PIGLET method to model
the quantum nature of the nuclei.  Hence, both electronic and nuclear
quantum effects are considered, and contribute to define the details
of the hydrogen bond network.

The PAMM analysis consists of three steps.

1. **Evaluation of the descriptors of the hydrogen bond**
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

2. **Optimization of the PAMM model**
Next, one can run the PAMM analysis on the set of HB parameters:

```bash
    ../bin/pamm -d 3 -w -o h2o -ngrid 2000 < h2o.nmr
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

3. **Post-processing the trajectory**

