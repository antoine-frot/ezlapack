# EzLAPACK: THE LAPACK wrapper

## What is it?

EzLAPACK is a LAPACK wrapper written in `Fortran 90`.

It's purpose is to use the speed of the LAPACK library without all this contraints
, namely the number of dummy arguments and the absence of generic templates.

## Installation

### Prerequisites

You should have gfortran, make, the BLAS and LAPACK library installed.

```
sudo apt install gfortran
sudo apt install make
sudo apt install liblapack-dev libblas-dev
```

### ezLAPACK Installation

The ezLAPACK wrapper can be downloaded on Github as Git repository. 
Then you need to go inside the directory and run the installation command.

```
git clone https://github.com/antoine-frot/ezlapack.git
cd ezlapack/
make install
```
You can now use ezLAPACK in your projects !

The directory can be deleted afterward.

## User guide

At the beginning of your programm, in which you want to use ezLAPACK put the instruction `use ezlapack`

```
program program_name
    use ezlapack
    ! Your code
end program program_name
```
The ezLAPACK wrapper should be compiled as other libraries, but required the LAPACK and BLAS library just as the location of the ezlapack.mod file.

```
gfortran -J/usr/local/include -o program_name program_name.f90 -lezlapack -lblas -llapack
```

CAUTION: -lezlapack should be place before -lblas -llapack.
-J indicate where is the ezlapack.mod file placed 
by default it is /usr/local/include but this can be change in the Makefile before installation at the line `PATH_MOD = /usr/local/include`

First 'sudo make install' then make oublie pas -lezlapack -lblas et pas l'inversse et -I/usr/local/include/

## Subroutines of ezLAPACK

### ezmatmul 

**Examples:**                                                                     
                                                                               
`call ezmatmul(A, B, C)` computes `C:= A*B` with the simplicity of matmul            
while harnessing the high performance of the LAPACK library.                   
                                                                               
On the other hand, `ezmatmul` keep the versatility of `gemm` with for example:     
`call ezmatmul('C', 'T', (5d-4,8d6), A, B, (9d2,3d-7), C)`                             

**Purpose:**                                                                       
Generic interface for LAPACK matrice multiplication `gemm` (i.e. `sgemm`, `dgemm`, `cgemm`, `zgemm`)    
i.e. performs `C := alpha*op( A )*op( B ) + beta*C`                              
where `alpha` and `beta` are scalars, `A`, `B` and `C` are matrices, and                 
`op( X )` is one of `op( X ) = X` or `op( X ) = X**T` or `op( X ) = X**H`.

**Subroutine:**                                                                    
Type can be `real(4)`, `real(8)`, `complex(4)` or `complex(8)`,                           
but all dummy arguments should have the same type.                             
```                                                                              
subroutine ezmatmul ( character*1, optional (default = 'N') :: transa,         
                      character*1, optional (default = 'N') :: transb,         
                      type,        optional (default =  1 ) :: alpha,          
                      type, dimension(:,:), contiguous      :: A,              
                      type, dimension(:,:), contiguous      :: B,              
                      type,        optional (default =  0 ) :: beta,           
                      type, dimension(:,:), contiguous      :: C,              
)                                                                              
```                                                                               
```                                                                               
transa (character\*1, optional (default = 'N')): specifies the form of op( A ). 
  if transa = 'N' or 'n',  op( A ) = A.                                        
  if transa = 'T' or 't',  op( A ) = A**T.                                     
  if transa = 'C' or 'c',  op( A ) = A**H.                                     
                                                                               
transb (character\*1, optional (default = 'N')): specifies the form of op( B ). 
  if transb = 'N' or 'n',  op( B ) = B.                                        
  if transb = 'T' or 't',  op( B ) = B**T.                                     
  if transb = 'C' or 'c',  op( B ) = B**H.                                     
                                                                               
alpha (type, optional (default = 1)): specifies the scalar alpha.              
                                                                               
A (type, dimension(:), contiguous): specifies the matrix A                     
                                                                               
B (type, dimension(:), contiguous): specifies the matrix B                     
                                                                               
beta (type, optional (default = 0)): specifies the scalar beta.                
                                                                               
C (type, dimension(:), contiguous): specifies the matrix C                     
```                                                                               

## Run Tests

Go inside the ezlapack directory and run:

```
make test
```
