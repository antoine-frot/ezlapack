
program test_matrice_mult_d

  use mod_ezmatmul

  implicit none

  integer          :: N, i
  double precision :: t_start, t_end
  integer          :: size_mat(4)
  double precision :: alpha, beta
  double precision, allocatable, dimension(:,:) :: A, B, C

  size_mat = (/ 10, 100, 1000, 2000/)

  do i = 1, size(size_mat)
    N = size_mat(i)

    allocate(A(N, N), B(N, N), C(N, N))

    call random_number(A)
    call random_number(B)

    alpha = 1d0
    beta = 0d0

    ! Matmul section
    call cpu_time(t_start)
    C = matmul(A, B)
    call cpu_time(t_end)
    write(*,'(A30,E12.6,A10,I9)') 'CPU time for matmul ', t_end - t_start, ' s for N = ', N

    ! DGEMM section
    call cpu_time(t_start)
    call dgemm('N', 'N', N, N, N, alpha, A, N, B, N, beta, C, N)
    call cpu_time(t_end)
    write(*,'(A30,E12.6,A10,I9)') 'CPU time for dgemm ', t_end - t_start, ' s for N = ', N

    ! Wrapper section
    call cpu_time(t_start)
    call ezmatmul(A, B, C)
    call cpu_time(t_end)
    write(*,'(A30,E12.6,A10,I9)') 'CPU time for wrapper ', t_end - t_start, ' s for N = ', N

    deallocate(A, B, C)
  end do

end program test_matrice_mult_d
