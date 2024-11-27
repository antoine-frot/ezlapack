# EzLAPACK: THE LAPACK wrapper

# What is it?

EzLAPACK is a LAPACK wrapper written in `Fortran 90`.

It's purpose is to use the speed of the LAPACK library without all this contraints

# Installation guide

The EzLAPACK wrapper can be downloaded on Github as Git repository.
```
git clone https://github.com/antoine-frot/ezlapack.git
```
The LAPACK and BLAS library should be installed.
```
apt install liblapack-dev libblas-dev
```
# User guide
The EzLAPACK wrapper should be compiled as any other module.

If you're still not using a compiler, you can install gfortran with
```
apt install gfortran
```
Install make to facilitate the compilation
```
apt install make
```
In your main program place just after program my-program-name
use `mod_ez_matmul`

# ez_matmul 
