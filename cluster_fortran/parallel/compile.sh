#!/bin/bash
f2py --fcompiler=gfortran --f90flags="-fopenmp" -lgomp -c -m ctp cluster_toolkit_parallel.f90 only: clus_sig_p

