program test_matrice_mult_d

  use mod_ez_matmul

  implicit none

  integer          :: N, i
  double precision :: t_start, t_end
  integer          :: size_mat(4)
  double precision :: alpha, beta
  double precision, allocatable, dimension(:,:) :: A, B, C, D

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
    call ez_matmul(A, B, C)
    call cpu_time(t_end)
    write(*,'(A30,E12.6,A10,I9)') 'CPU time for wrapper ', t_end - t_start, ' s for N = ', N

    deallocate(A, B, C)
  end do

  !Test of all default arguments
  N=10000
  print *, "Matrices de taille ", N, " * ", N
  allocate(A(N, N), B(N, N), C(N, N), D(N,N))

  call random_number(A)
  call random_number(B)

  !D = matmul(A, B)
  call dgemm('N', 'N', N, N, N, 1d0, A, N, B, N, 0d0, D, N)

  call ez_matmul(A, B, C)
  if (any(D /= C)) then
    print *, "Test failed (1)"
    print *, maxval(abs(D-C))
  end if

  call ez_matmul('N', 'N', A, B, C)
  if (any(D /= C)) then
    print *, "Test failed (2)"
    print *, maxval(abs(D-C))
  end if

  call ez_matmul(1d0, A, B, C)
  if (any(D /= C)) then
    print *, "Test failed (3)"
    print *, maxval(abs(D-C))
  end if

  call ez_matmul(A, B, 0d0, C)
  if (any(D /= C)) then
    print *, "Test failed (4)"
    print *, maxval(abs(D-C))
  end if

  call ez_matmul('N', 'N', A, B, 0d0, C)
  if (any(D /= C)) then
    print *, "Test failed (5)"
    print *, maxval(abs(D-C))
  end if

  call ez_matmul('N', 'N', 1d0, A, B, C)
  if (any(D /= C)) then
    print *, "Test failed (6)"
    print *, maxval(abs(D-C))
  end if

  call ez_matmul(1d0, A, B, 0d0, C)
  if (any(D /= C)) then
    print *, "Test failed (7)"
    print *, maxval(abs(D-C))
  end if
  
  !print *, 'C: ', C
  !print *, 'D: ', D

  !Test dgemm directly with ez_matmul

  call random_number(A)
  call random_number(B)
  call random_number(C)
  D = C

  call dgemm('T', 'T', N, N, N, 4.7d0, A, N, B, N, 6.9d0, D, N)
  call ez_matmul('T', 'T', 4.7d0, A, B, 6.9d0, C)
  

  if (any(D /= C)) then
    print *, "Test failed: dgemm doesn't match ez_matmul"
    print *, maxval(abs(D-C))
  end if

  print *, "Test finished"

end program test_matrice_mult_d

