# CHROMIBD

CHROMIBD models a set of “target chromosomes” as a mosaic of “reference chromosomes”. It can work within a genealogy and the set of reference chromosomes are the “parental chromosomes” of the target chromosome. It works also with unrelated individuals and without genealogy.

The programs are written in Fortran. A Makefile is provided to help for compilation.

Start with Make clean to remove previous file and then type Make:

Make clean
Make

The Makefile is written for the gfortran compiler but can easily be modified for other compilers.
