program speed_test_ezmatmul
! Purpose: Compare wall time and CPU time for matrix multiplication using
! MATMUL, DGEMM, and EZMATMUL across different matrix sizes.

  use ezlapack, only: ezmatmul

  implicit none

  integer                              :: N, i
  real(8)                              :: t_cpu_start_matmul, t_cpu_end_matmul
  real(8)                              :: t_cpu_start_dgemm, t_cpu_end_dgemm
  real(8)                              :: t_cpu_start_ezmatmul, t_cpu_end_ezmatmul
  real(8)                              :: t_wall_start, t_wall_end
  integer                              :: start_count, end_count, count_rate
  integer                              :: size_mat(6)
  real(8), allocatable, dimension(:,:) :: A, B, C

  size_mat = (/ 10, 100, 1000, 2000, 5000, 10000/)

  call SYSTEM_CLOCK(count_rate = count_rate)

  do i = 1, size(size_mat)

    print *, ''

    N = size_mat(i)
    write(*,'(A,I0,A,I0)') 'Size of the matrices is ', N, ' * ', N
    allocate(A(N, N), B(N, N), C(N, N))
    call random_number(A)
    call random_number(B)

    ! Matmul section
    if (N < 10000) then ! Too slow
      call cpu_time(t_cpu_start_matmul)
      call SYSTEM_CLOCK(start_count)
      C = matmul(A, B)
      call SYSTEM_CLOCK(end_count)
      call cpu_time(t_cpu_end_matmul)
      t_wall_start = real(start_count) / real(count_rate)
      t_wall_end   = real(end_count) / real(count_rate)
      write(*,'(A,E12.6,A)') 'Wall time for matmul   ', t_wall_end - t_wall_start, ' seconds'
    else
      write(*,'(A)') 'Matmul is too slow for this size of matrices'
    end if

    ! DGEMM section
    call cpu_time(t_cpu_start_dgemm)
    call SYSTEM_CLOCK(start_count)
    call dgemm('N', 'N', N, N, N, 1d0, A, N, B, N, 0d0, C, N)
    call SYSTEM_CLOCK(end_count)
    call cpu_time(t_cpu_end_dgemm)
    t_wall_start = real(start_count) / real(count_rate)
    t_wall_end   = real(end_count) / real(count_rate)
    write(*,'(A,E12.6,A)') 'Wall time for dgemm    ', t_wall_end - t_wall_start, ' seconds'

    ! Ezmatmul section
    call cpu_time(t_cpu_start_ezmatmul)
    call SYSTEM_CLOCK(start_count)
    call ezmatmul(A, B, C)
    call SYSTEM_CLOCK(end_count)
    call cpu_time(t_cpu_end_ezmatmul)
    t_wall_start = real(start_count) / real(count_rate)
    t_wall_end   = real(end_count) / real(count_rate)
    write(*,'(A,E12.6,A)') 'Wall time for ezmatmul ', t_wall_end - t_wall_start, ' seconds'

    ! Print all CPU times after wall times
    print *, ''
    if (N < 10000) then ! Too slow
      write(*,'(A,E12.6,A)') 'CPU time for matmul    ', t_cpu_end_matmul   - t_cpu_start_matmul, ' seconds'
    else
      write(*,'(A)') 'Matmul is too slow for this size of matrices'
    end if
    write(*,'(A,E12.6,A)') 'CPU time for dgemm     ', t_cpu_end_dgemm    - t_cpu_start_dgemm, ' seconds'
    write(*,'(A,E12.6,A)') 'CPU time for ezmatmul  ', t_cpu_end_ezmatmul - t_cpu_start_ezmatmul, ' seconds'

    deallocate(A, B, C)

  end do

end program speed_test_ezmatmul

