# EzLAPACK: The LAPACK Wrapper

## Overview

EzLAPACK is a lightweight wrapper for the LAPACK library, written in `Fortran 90`.
It provides an easy-to-use interface for leveraging LAPACK's performance while minimizing its common constraints, such as the large number of dummy arguments and the lack of generic templates.

---

## Installation

### Prerequisites

Before installing EzLAPACK, ensure that the following tools and libraries are installed on your system:

- `gfortran`
- `make`
- BLAS and LAPACK libraries

To install these on a Debian-based system (e.g., Ubuntu), run:

```bash
sudo apt install gfortran make liblapack-dev libblas-dev
```

To install these on a macOS system, run:

```bash
brew install gfortran make lapack
```

### Installing EzLAPACK

1. Clone the EzLAPACK repository from GitHub:

    ```bash
    git clone https://github.com/antoine-frot/ezlapack.git
    ```

2. Navigate to the EzLAPACK directory:

    ```bash
    cd ezlapack/
    ```

3. Install EzLAPACK using `make`:

    ```bash
    make install
    ```

4. You can now use EzLAPACK in your projects! After installation, the repository directory can be deleted if desired.

---

## User Guide

### Using EzLAPACK in Your Program

To integrate EzLAPACK, include the following statement at the beginning of your Fortran program:

```fortran
program program_name
    use ezlapack
    ! Your code here
end program program_name
```

When compiling your program, ensure that the EzLAPACK module is linked properly. For example:

```bash
gfortran -I/usr/local/include -o program_name program_name.f90 -lezlapack -lblas -llapack
```

#### Notes:
- **Library Order:** Place `-lezlapack` before `-lblas -llapack`.
- **Module Path:** The `-J` flag specifies the directory of the `ezlapack.mod` file. By default, it is `/usr/local/include`. This can be customized in the `Makefile` by modifying the `PATH_MOD` variable before installation.
```make
PATH_MOD = /your/custom/path
```

Then reinstall EzLAPACK using:

```bash
sudo make install
```

The use of `use ezlapack` and the flag `-J/usr/local/include` is essential due to the presence of the generic interface.

---

## Subroutines in EzLAPACK

### `ezmatmul`

#### Purpose:
`ezmatmul` provides a user-friendly interface for LAPACK's matrix multiplication routines (`gemm`, such as `sgemm`, `dgemm`, `cgemm`, `zgemm`).

It computes:

```text
C := alpha * op(A) * op(B) + beta * C
```

where:
- `alpha` and `beta` are scalars.
- `A`, `B`, and `C` are matrices.
- `op(X)` can be `X`, `X**T`, or `X**H`.

#### Example Usage:

1. Basic Matrix Multiplication:

    ```fortran
    call ezmatmul(A, B, C)
    ```

    This computes `C := A * B` using LAPACK's performance while maintaining the simplicity of Fortran's `matmul`.

2. Advanced Usage:

    ```fortran
    call ezmatmul('C', 'T', 8d6, A, B, -3d-7, C)
    ```

    This computes `C := 8d6 * A**H * B**T - 3d-7 * C`

#### Subroutine Definition:

```fortran
subroutine ezmatmul(
    character*1, optional (default = 'N') :: transa,
    character*1, optional (default = 'N') :: transb,
    type,        optional (default =  1 ) :: alpha,
    type, dimension(:,:), contiguous      :: A,
    type, dimension(:,:), contiguous      :: B,
    type,        optional (default =  0 ) :: beta,
    type, dimension(:,:), contiguous      :: C
)
```

#### Parameters:

- `transa`: (optional, default = `'N'`) Specifies the operation applied to matrix `A`:
    - `'N'`: No transpose, `op(A) = A`.
    - `'T'`: Transpose, `op(A) = A**T`.
    - `'C'`: Hermitian, `op(A) = A**H`.

- `transb`: (optional, default = `'N'`) Specifies the operation applied to matrix `B`:
    - `'N'`: No transpose, `op(B) = B`.
    - `'T'`: Transpose, `op(B) = B**T`.
    - `'C'`: Hermitian, `op(B) = B**H`.

- `alpha`: (optional, default = `1`) Scalar multiplier for `op(A) * op(B)`.

- `beta`: (optional, default = `0`) Scalar multiplier for `C`.

- `A`, `B`, `C`: Conitguous matrices with the same type.

#### Supported Types:
- `real(4)`
- `real(8)`
- `complex(4)`
- `complex(8)`

All arguments must have the same type.

---

## Running Tests

To verify your EzLAPACK installation, navigate to the EzLAPACK directory and run:

```bash
make test
```

This will compile and execute test cases to ensure everything works as expected.

## Coding Convention

The code follows the coding conventions outlined in the [Fortran-lang Style Guide](https://fortran-lang.org/en/learn/best_practices/style_guide/). Subroutine names are designed to closely align with the names of intrinsic Fortran subroutines, while variable names adhere to the LAPACK names.

## Contatct

For questions, suggestions, or feedback, please contact:
- Email: [Antoine Frot](mailto:antoine.frot@orange.fr)
- GitHub: [EzLAPACK Repository](https://github.com/antoine-frot/ezlapack)
