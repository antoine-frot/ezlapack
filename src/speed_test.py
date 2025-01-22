"""
Matrix Multiplication Performance Comparison: NumPy vs Ezlapack vs Lapack vs Fortran Intrinsic

This script compares the performance of matrix multiplication of:
- NumPy (Python)
- Ezlapack (Fortran wrapper for LAPACK)
- LAPACK (Library of optimized numerical algorithms)
- Fortran Intrinsic (Fortran's built-in matrix multiplication function)

The script performs matrix multiplication for varying matrix sizes (10, 100, 1000, 2000, 5000, 10000). 
The results are saved to a text file and displayed in the terminal. Additionally, a logarithmic plot of the performance comparison is generated.
"""

import numpy as np
import time
import matplotlib.pyplot as plt

# Define the matrix sizes to test
matrix_sizes = [10, 100, 1000, 2000, 5000, 10000]

def test_python_matrix_multiplication():
    """
    Tests the performance of Python's NumPy matrix multiplication 
    (np.dot) for different matrix sizes and records the time taken 
    for each multiplication.

    Returns:
        list: A list of elapsed times for each matrix size.
    """
    python_times = []
    for size in matrix_sizes:
        A = np.random.random((size, size))
        B = np.random.random((size, size))
        start_time = time.time()
        np.dot(A, B)
        elapsed_time = time.time() - start_time
        python_times.append(elapsed_time)
    return python_times

def read_fortran_results(file_path):
    """
    Reads the results of matrix multiplication performance from 
    speed_test.f90 output and returns the times for Ezlapack, Lapack, 
    and Fortran intrinsic (Matmul).

    Args:
        file_path (str): Path to speed_test.f90 output.

    Returns:
        tuple: A tuple containing three lists:
               - ezlapack_times (list of floats)
               - lapack_times (list of floats)
               - matmul_times (list of floats)
    """
    ezlapack_times = []
    lapack_times = []
    matmul_times = []
    with open(file_path, 'r') as file:
        lines = file.readlines()[2:]  # Skip headers
        for line in lines:
            _, time_ezlapack, time_lapack, time_matmul = line.split()
            ezlapack_times.append(float(time_ezlapack))
            lapack_times.append(float(time_lapack))
            matmul_times.append(float(time_matmul))
    return ezlapack_times, lapack_times, matmul_times

def write_and_compare_results(output_file, python_times, ezlapack_times, lapack_times, matmul_times):
    """
    Writes and plots the comparaison of the matrix multiplication performance 
    of NumPy, Ezlapack, Lapack, and Fortran intrinsic matrix multiplication.

    Args:
        output_file (str): The file path to save the results.
        python_times (list): List of elapsed times for Python's NumPy.
        ezlapack_times (list): List of elapsed times for Ezlapack.
        lapack_times (list): List of elapsed times for Lapack.
        matmul_times (list): List of elapsed times for Fortran intrinsic.
    """
    with open(output_file, 'w') as file:
        file.write("Matrix Multiplication Performance: NumPy vs Ezlapack vs Lapack vs Fortran Intrinsic\n")
        file.write(f"{'Matrix Size':<15}{'NumPy Time (s)':<20}{'Ezlapack Time (s)':<20}{'Lapack Time (s)':<20}{'Fortran Intrinsic Time (s)':<20}\n")
        file.write("-" * 101 + "\n")
        for size, py_time, ezlapack_time, lapack_time, matmul_time in zip(matrix_sizes, python_times, ezlapack_times, lapack_times, matmul_times):
            file.write(f"{size:<15}{py_time:<20.6f}{ezlapack_time:<20.6f}{lapack_time:<20.6f}{matmul_time:<20.6f}\n")

    # Display the content of the file
    with open(output_file, 'r') as file:
        print(file.read())

    # Plot the results
    plt.plot(matrix_sizes, python_times, label="NumPy")
    plt.plot(matrix_sizes, ezlapack_times, label="Ezlapack")
    plt.plot(matrix_sizes, lapack_times, label="Lapack")
    plt.plot(matrix_sizes, matmul_times, label="Fortran Intrinsic")
    plt.xlabel("Matrix Size")
    plt.ylabel("Time (seconds)")
    plt.title("Matrix Multiplication Performance: NumPy vs Ezlapack vs Lapack vs Fortran Intrinsic")
    plt.legend()
    plt.grid(True)
    plt.xscale("log")
    plt.yscale("log")
    plt.show()

# Main execution
if __name__ == "__main__":
    python_times = test_python_matrix_multiplication()
    ezlapack_times, lapack_times, matmul_times = read_fortran_results("src/fortran_results.txt")
    write_and_compare_results("src/speed_test.txt", python_times, ezlapack_times, lapack_times, matmul_times)

