program speed_test_ezmatmul
! Purpose: Compare wall time for matrix multiplication using
! MATMUL (intrinsic), DGEMM (lapack), and EZMATMUL (ezlapack) across different matrix sizes.

  use ezlapack, only: ezmatmul

  implicit none

  integer                              :: N, i
  real(8)                              :: t_wall_start, t_wall_end
  real(8)                              :: time_ezlapack, time_lapack, time_intrinsic
  integer                              :: start_count, end_count, count_rate
  integer                              :: size_mat(6)
  real(8), allocatable, dimension(:,:) :: A, B, C

  size_mat = (/ 10, 100, 1000, 2000, 5000, 10000/)

  call SYSTEM_CLOCK(count_rate = count_rate)

  print *, "Matrix Size ", "EZLAPACK Time (s) ", "LAPACK Time (s) ", "MatMul Time (s) "
  print *, "-------------------------------------------------------------"

  do i = 1, size(size_mat)

    N = size_mat(i)
    allocate(A(N, N), B(N, N), C(N, N))
    call random_number(A)
    call random_number(B)

    ! Matmul section
    call SYSTEM_CLOCK(start_count)
    C = matmul(A, B)
    call SYSTEM_CLOCK(end_count)
    t_wall_start = real(start_count) / real(count_rate)
    t_wall_end   = real(end_count) / real(count_rate)
    time_intrinsic = t_wall_end - t_wall_start

    ! DGEMM section
    call SYSTEM_CLOCK(start_count)
    call dgemm('N', 'N', N, N, N, 1d0, A, N, B, N, 0d0, C, N)
    call SYSTEM_CLOCK(end_count)
    t_wall_start = real(start_count) / real(count_rate)
    t_wall_end   = real(end_count) / real(count_rate)
    time_lapack = t_wall_end - t_wall_start

    ! Ezmatmul section
    call SYSTEM_CLOCK(start_count)
    call ezmatmul(A, B, C)
    call SYSTEM_CLOCK(end_count)
    t_wall_start = real(start_count) / real(count_rate)
    t_wall_end   = real(end_count) / real(count_rate)
    time_ezlapack = t_wall_end - t_wall_start

    write(*,'(I12,F17.7,F15.7,F15.7)') N, time_ezlapack, time_lapack, time_intrinsic

    deallocate(A, B, C)

  end do

end program speed_test_ezmatmul

