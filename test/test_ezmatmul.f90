program test_ezmatmul
  !
  ! Test file for ezmatmul
  use module_ezmatmul,       only ezmatmul
  use module_random_complex, only random_complex

  implicit none

  logical                                       :: test_passed
  integer                                       :: N
    ! Size of the matrices
  double precision                              :: alpha, beta
  double precision, allocatable, dimension(:,:) :: A, B, C, D
  

  ! Test of all default arguments
  test_passed = .true.
  N=1000
  print *, "Matrices de taille ", N, " * ", N
  allocate(A(N, N), B(N, N), C(N, N), D(N,N))

  call random_number(A)
  call random_number(B)

  call dgemm('N', 'N', N, N, N, 1d0, A, N, B, N, 0d0, D, N)

  call ezmatmul(A, B, C)
  if (any(D /= C)) then
    print *, "Test failed (1)"
    print *, 'Maximal difference between dgemm and ezmatmum: ', maxval(abs(D-C))
    test_passed = .false.
  end if

  call ezmatmul('N', 'N', A, B, C)
  if (any(D /= C)) then
    print *, "Test failed (2)"
    print *, 'Maximal difference between dgemm and ezmatmum: ', maxval(abs(D-C))
    test_passed = .false.
  end if

  call ezmatmul(1d0, A, B, C)
  if (any(D /= C)) then
    print *, "Test failed (3)"
    print *, 'Maximal difference between dgemm and ezmatmum: ', maxval(abs(D-C))
    test_passed = .false.
  end if

  call ezmatmul(A, B, 0d0, C)
  if (any(D /= C)) then
    print *, "Test failed (4)"
    print *, 'Maximal difference between dgemm and ezmatmum: ', maxval(abs(D-C))
    test_passed = .false.
  end if

  call ezmatmul('N', 'N', A, B, 0d0, C)
  if (any(D /= C)) then
    print *, "Test failed (5)"
    print *, 'Maximal difference between dgemm and ezmatmum: ', maxval(abs(D-C))
    test_passed = .false.
  end if

  call ezmatmul('N', 'N', 1d0, A, B, C)
  if (any(D /= C)) then
    print *, "Test failed (6)"
    print *, 'Maximal difference between dgemm and ezmatmum: ', maxval(abs(D-C))
    test_passed = .false.
  end if

  call ezmatmul(1d0, A, B, 0d0, C)
  if (any(D /= C)) then
    print *, "Test failed (7)"
    print *, 'Maximal difference between dgemm and ezmatmum: ', maxval(abs(D-C))
    test_passed = .false.
  end if
  
  ! An additionnal random test with all dummy arguments
  call random_number(C)
  D = C
  call dgemm('T', 'T', N, N, N, 4.7d0, A, N, B, N, 6.9d0, D, N)
  call ezmatmul('T', 'T', 4.7d0, A, B, 6.9d0, C)
  

  if (any(D /= C)) then
    print *, "Random test with all dummy arguments failed for dgemm"
    print *, 'Maximal difference between dgemm and ezmatmum: ', maxval(abs(D-C))
    test_passed = .false.
  end if

  print *, "Test finished for dgemm"

  if (.not. test_passed) then
    print *, '!-------------!'
    print *, '! TEST FAILED !'
    print *, '!-------------!'

end program test_ezmatmul
