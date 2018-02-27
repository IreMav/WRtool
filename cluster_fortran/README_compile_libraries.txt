The WRtool uses the compiled libraries ctool.so and ctp.so; in this folder you
can find the source code. You need to compile those modules on your system
before using them, only the first time you use the tool.

Compile them with:

f2py --fcompiler=gfortran --f90flags="-fopenmp" -lgomp -c -m ctp cluster_toolkit_parallel.f90 only: clus_sig_p

f2py --fcompiler=gfortran -c -m ctool cluster_toolkit.f90 only: clus_opt adran1 gausts tsstat



