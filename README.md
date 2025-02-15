# EzLAPACK: The LAPACK Wrapper

## Overview

EzLAPACK is a lightweight wrapper for the LAPACK library, written in `Fortran 90`.
It provides an easy-to-use interface for leveraging LAPACK's performance while minimizing its common constraints, such as the large number of arguments and the lack of generic templates.

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

#### Installing Locally

Install EzLAPACK locally using:

    ```bash
    make install
    ```

You can now move the directory `lib` in your projects to use it.

#### Installing Globally (Only Available for Linux)

Install EzLAPACK globally using:

    ```bash
    make install_global
    ```

By default, the static library is installed in `/usr/local/lib`, and the `.mod` files are placed in `/usr/local/include`. 
You can customize these installation paths by modifying the `PATH_LIBRARY` and `PATH_MOD` variables in the Makefile to suit your preferences.

### Uninstall EzLAPACK

Uninstall EzLAPACK using:

    ```bash
    make uninstall
    ```

After that, there will be no trace of the library remaining.

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

When compiling your program, ensure that the EzLAPACK module is linked properly by using the following command:

```bash
gfortran -o program_name program_name.f90 -I/path/to/mod/files -L/path/to/lib -lezlapack -lblas -llapack
```

#### Notes:
- **Local Installation:** `-I` and `-L` flags should point to the `lib` directory.
- **Global Installation:** Use `-I/usr/local/include` (default), `-L` flag can be omitted
- **Library Order:** Place `-lezlapack` before `-lblas -llapack`.

The use of `use ezlapack` and the flag `-J` is essential due to the presence of the generic interfaces.

---

## Subroutines in EzLAPACK

### `ezmatmul`

#### Purpose:
`ezmatmul` provides a user-friendly generic interface for LAPACK's matrix multiplication routines `*gemm` (i.e. `sgemm`, `dgemm`, `cgemm`, `zgemm`).

It computes:

```text
C := alpha*op(A)*op(B) + beta*C
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

Type can be either `real(4)`, `real(8)`, `complex(4)` or `complex(8)`, 
but all arguments should have the same type.

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

- `A`, `B`, `C`: Contiguous matrices with the same type.

---

## Additional Feature: Random Library

Fortran’s built-in `random_number` function is limited to generating random real numbers. To extend this functionality, additional subroutines have been created, particularly for testing purposes. These subroutines are included in the library and are outlined below.

### `random_complex`

```fortran
random_complex(output_complex)
```

Generate random complex numbers (single or double precision) within the unit circle. Supports scalars and arrays up to rank-2.

### `random_integer`

```fortran
random_integer(low, high, output_integer)
```

Generate a random integer between `low` and `high`, inclusive.

### `random_complex`

```fortran
random_character(character_array, output_character)
```

Randomly select and return a character from an input array of strings.

---

## Speed Test

The LAPACK library is known for its speed, but performance can vary across different computers. 
A speed test script has been created to measure LAPACK's execution time and compare it with the intrinsic `matmul` function, the `EzLAPACK` wrapper and `NumPy` from `Python`.

To run the speed test, use the following command:

```bash
make run
```

As usually observed, NumPy is as fast as EzLAPACK since it is build on LAPACK.

---

## Running Tests

To verify your EzLAPACK installation, navigate to the EzLAPACK directory and run:

```bash
make test
```

This will compile and execute test cases to ensure everything works as expected.

---

## Dependencies

This project relies on the following tools and libraries:

- Fortran Compiler: gfortran 13.3.0
- LAPACK: 3.11.0
- BLAS: 3.11.0
- Operating System: Ubuntu 22.04 (WSL)
- Bash: 5.1.16
- Make: 4.3
- Python: 3.12.3

We recommend replicating this setup for best results.

---

## Coding Convention

The code follows the coding conventions outlined in the 
[Fortran-lang Style Guide](https://fortran-lang.org/en/learn/best_practices/style_guide/).
Subroutine names are designed to closely align with the names of intrinsic Fortran subroutines, 
while variable names adhere to the LAPACK names.
The documentation adheres to the nomenclature conventions of Doxygen.

---

## Contact

For questions, suggestions, or feedback, please contact:
- Email: [Antoine Frot](mailto:antoine.frot@orange.fr)
- GitHub: [EzLAPACK Repository](https://github.com/antoine-frot/ezlapack)
