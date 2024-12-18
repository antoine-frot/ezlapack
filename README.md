# EzLAPACK: THE LAPACK wrapper

## What is it?

EzLAPACK is a LAPACK wrapper written in `Fortran 90`.

It's purpose is to use the speed of the LAPACK library without all this contraints

## Installation guide

The EzLAPACK wrapper can be downloaded on Github as Git repository.
```
git clone https://github.com/antoine-frot/ezlapack.git
```
The LAPACK and BLAS library should be installed.
```
apt install liblapack-dev libblas-dev
```
## User guide
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

## ezmatmul 
**Purpose:**                                                                       
Generic interface for lapack matrice multiplication (*gemm)                    
i.e. performs C := alpha*op( A )*op( B ) + beta*C                              
where alpha and beta are scalars, A, B and C are matrices, and                 
op( X ) is one of op( X ) = X or op( X ) = X**T or op( X ) = X**H.             

**Subroutine:**                                                                    
Type can be real(4), real(8), complex or complex*16,                           
but all dummy arguments should have the same type.                             
                                                                               
subroutine ezmatmul	(	character*1, optional (default = 'N') :: transa,         
                      character*1, optional (default = 'N') :: transb,         
                      type,        optional (default =  1 ) :: alpha,          
                      type, dimension(:,:), contiguous      :: A,              
                      type, dimension(:,:), contiguous      :: B,              
                      type,        optional (default =  0 ) :: beta,           
                      type, dimension(:,:), contiguous      :: C,              
)                                                                              
                                                                               
transa (character*1, optional (default = 'N')): specifies the form of op( A ). 
  if transa = 'N' or 'n',  op( A ) = A.                                        
  if transa = 'T' or 't',  op( A ) = A**T.                                     
  if transa = 'C' or 'c',  op( A ) = A**H.                                     
                                                                               
transb (character*1, optional (default = 'N')): specifies the form of op( B ). 
  if transb = 'N' or 'n',  op( B ) = B.                                        
  if transb = 'T' or 't',  op( B ) = B**T.                                     
  if transb = 'C' or 'c',  op( B ) = B**H.                                     
                                                                               
alpha (type, optional (default = 1)): specifies the scalar alpha.              
                                                                               
A (type, dimension(:), contiguous): specifies the matrix A                     
                                                                               
B (type, dimension(:), contiguous): specifies the matrix B                     
                                                                               
beta (type, optional (default = 0)): specifies the scalar beta.                
                                                                               
C (type, dimension(:), contiguous): specifies the matrix C                     

**Examples:**                                                                      
                                                                               
call ezmatmul(A,B,C) computes C:= A*B with the simplicity of matmul            
while harnessing the high performance of the LAPACK library.                   
                                                                               
On the other hand, ezmatmul keep the versatility of gemm with for example:     
call ezmatmul('C','T',(5d-4,8d6),A,B,(9d2,3d-7),C)                             
