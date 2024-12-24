# EzLAPACK: The LAPACK Wrapper

## What is EzLAPACK?

EzLAPACK is a lightweight wrapper for the LAPACK library, written in `Fortran 90`. It is designed to:

- Simplify the use of the LAPACK library.
- Minimize constraints such as the large number of dummy arguments.
- Compensate for the lack of generic templates in Fortran.

With EzLAPACK, you can harness the speed and efficiency of LAPACK without the usual complexity.

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

To use EzLAPACK, include the following statement at the beginning of your Fortran program:

```fortran
program program_name
    use ezlapack
    ! Your code here
end program program_name
```

When compiling your program, ensure that the EzLAPACK module is linked properly. For example:

```bash
gfortran -J/usr/local/include -o program_name program_name.f90 -lezlapack -lblas -llapack
```

#### Notes:
- **Order of Libraries:** Place `-lezlapack` before `-lblas -llapack`.
- **Module Path:** The `-J` flag specifies the directory of the `ezlapack.mod` file. By default, it is `/usr/local/include`. This can be customized in the `Makefile` by modifying the `PATH_MOD` variable before installation.

For example, update the following line in the `Makefile`:

```make
PATH_MOD = /your/custom/path
```

Then reinstall using:

```bash
sudo make install
```

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

1. Simplified matrix multiplication:

    ```fortran
    call ezmatmul(A, B, C)
    ```

    This computes `C := A * B` using LAPACK's performance while maintaining the simplicity of Fortran's `matmul`.

2. Advanced usage:

    ```fortran
    call ezmatmul('C', 'T', (5d-4, 8d6), A, B, (9d2, 3d-7), C)
    ```

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
    - `'N'`: `op(A) = A`
    - `'T'`: `op(A) = A**T`
    - `'C'`: `op(A) = A**H`

- `transb`: (optional, default = `'N'`) Specifies the operation applied to matrix `B`:
    - `'N'`: `op(B) = B`
    - `'T'`: `op(B) = B**T`
    - `'C'`: `op(B) = B**H`

- `alpha`: (optional, default = `1`) Scalar multiplier for `op(A) * op(B)`.

- `A`: Matrix `A` (must be contiguous).

- `B`: Matrix `B` (must be contiguous).

- `beta`: (optional, default = `0`) Scalar multiplier for `C`.

- `C`: Matrix `C` (must be contiguous).

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

