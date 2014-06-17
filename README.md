# PAMM library and tools

This repository contains a simple Fortran90 implementation of 
the Probabilistic Analysis of Molecular Motifs algorithm. 
It contains a library to compute PAMM structure counts, a 
tool to perform the analysis on an arbitrary set of structural
descriptors, and an example application to the automatic 
definition of hydrogen bond descriptors. 

Source code is stored in the `src/` directory, and the executables
will be generated inside `bin/`. `test/` contains a simple
example of the usage of these programs. 


## Compilation and installation

Compilation should be trivial, requiring only a recent version
of `gfortran` and LAPACK libraries. You can adjust the compiler
and the path of the libraries modifying `src/Makefile`. You
should then be able to compile the library with the commands

    cd src/
    make

The executables will be generated inside the `bin/` directory.
Make sure to copy them in your path, or to add the folder to
your path environment variable:

    export PATH=$PATH:$PWD/bin/

## Getting started

An example of the usage of `pamm` and `hbpamm` can be found
in the `test/` directory.
