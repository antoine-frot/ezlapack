program test_ezmatmul

  use ezlapack

  implicit none

  logical                              :: test_passed
  character(len=1)                     :: transa, transb
  integer                              :: m, n, k
  real(8)                              :: alpha_dp, beta_dp
  real(8), allocatable, dimension(:,:) :: A_dp, B_dp, C_dp, D_dp
  

  call random_integer(1, 3000, m)
  call random_integer(1, 3000, n)
  call random_integer(1, 3000, k)
  write(*, '(A, "A: (", I5, ", ", I5, ")")') "Size of the matrices: ", m, k
  write(*, '(A, "B: (", I5, ", ", I5, ")")') "                      ", k, n
  write(*, '(A, "C: (", I5, ", ", I5, ")")') "                      ", m, n

!------------!
! DGEMM PART !
!------------!

!----------------------------!
! Handling default arguments !
!----------------------------!

  print *, "Test for dgemm"
  test_passed = .true.
  allocate(A_dp(m, k), B_dp(k, n), C_dp(m, n), D_dp(m, n))
  call random_number(A_dp)
  call random_number(B_dp)

  ! Test of all default arguments
  call dgemm('N', 'N', m, n, k, 1d0, A_dp, m, B_dp, k, 0d0, D_dp, m)

  call ezmatmul(A_dp, B_dp, C_dp)
  if (any(D_dp /= C_dp)) then
    print *, "Test failed for ezdgemm_trans_alpha_beta"
    print *, 'Maximal difference between dgemm and ezmatmul: ', maxval(abs(D_dp - C_dp))
    test_passed = .false.
  end if

  call ezmatmul('N', 'N', A_dp, B_dp, C_dp)
  if (any(D_dp /= C_dp)) then
    print *, "Test failed for ezdgemm_alpha_beta"
    print *, 'Maximal difference between dgemm and ezmatmul: ', maxval(abs(D_dp - C_dp))
    test_passed = .false.
  end if

  call ezmatmul(1d0, A_dp, B_dp, C_dp)
  if (any(D_dp /= C_dp)) then
    print *, "Test failed for ezdgemm_trans_beta"
    print *, 'Maximal difference between dgemm and ezmatmul: ', maxval(abs(D_dp - C_dp))
    test_passed = .false.
  end if

  call ezmatmul(A_dp, B_dp, 0d0, C_dp)
  if (any(D_dp /= C_dp)) then
    print *, "Test failed for ezdgemm_trans_alpha"
    print *, 'Maximal difference between dgemm and ezmatmul: ', maxval(abs(D_dp - C_dp))
    test_passed = .false.
  end if

  call ezmatmul('N', 'N', A_dp, B_dp, 0d0, C_dp)
  if (any(D_dp /= C_dp)) then
    print *, "Test failed for ezdgemm_alpha"
    print *, 'Maximal difference between dgemm and ezmatmul: ', maxval(abs(D_dp - C_dp))
    test_passed = .false.
  end if

  call ezmatmul('N', 'N', 1d0, A_dp, B_dp, C_dp)
  if (any(D_dp /= C_dp)) then
    print *, "Test failed for ezdgemm_beta"
    print *, 'Maximal difference between dgemm and ezmatmul: ', maxval(abs(D_dp - C_dp))
    test_passed = .false.
  end if

  call ezmatmul(1d0, A_dp, B_dp, 0d0, C_dp)
  if (any(D_dp /= C_dp)) then
    print *, "Test failed for ezdgemm_trans"
    print *, 'Maximal difference between dgemm and ezmatmul: ', maxval(abs(D_dp - C_dp))
    test_passed = .false.
  end if
  
  call ezmatmul('N', 'N', 1d0, A_dp, B_dp, 0d0, C_dp)
  if (any(D_dp /= C_dp)) then
    print *, "Test failed for ezdgemm"
    print *, 'Maximal difference between dgemm and ezmatmul: ', maxval(abs(D_dp - C_dp))
    test_passed = .false.
  end if
  
!--------------------------------!
! Test with all arguments random !
!--------------------------------!

  call random_character((/ 'N', 'T', 'C'/), transa)
  call random_character((/ 'N', 'T', 'C'/), transb)

  call random_number(alpha_dp)
  call random_number(beta_dp)
  alpha_dp = (alpha_dp - 0.5 ) * 10000
  beta_dp  = (beta_dp  - 0.5 ) * 10000
  call random_number(C_dp)
  D_dp = C_dp

  if (transa == 'N') then
    if (transb == 'N') then
      call dgemm(transa, transb, m, n, k, alpha_dp, A_dp, m, B_dp, k, beta_dp, D_dp, m)
    else
      deallocate(B_dp)
      allocate(B_dp(n,k))
      call dgemm(transa, transb, m, n, k, alpha_dp, A_dp, m, B_dp, n, beta_dp, D_dp, m)
    end if
  else 
    deallocate(A_dp)
    allocate(A_dp(k,m))
    if (transb == 'N') then
      call dgemm(transa, transb, m, n, k, alpha_dp, A_dp, k, B_dp, k, beta_dp, D_dp, m)
    else
      deallocate(B_dp)
      allocate(B_dp(n,k))
      call dgemm(transa, transb, m, n, k, alpha_dp, A_dp, k, B_dp, n, beta_dp, D_dp, m)
    end if
  end if

  call ezmatmul(transa, transb, alpha_dp, A_dp, B_dp, beta_dp, C_dp)

  if (any(D_dp /= C_dp)) then
    print *, "Random test with all dummy arguments failed for dgemm"
    print *, 'Maximal difference between dgemm and ezmatmul: ', maxval(abs(D_dp - C_dp))
    print *, "transa: ", trim(transa)
    print *, "transb: ", trim(transb)
    print *, "alpha_dp:", alpha_dp
    print *, "beta_dp: ", beta_dp
    test_passed = .false.
  end if

  if (.not. test_passed) then
    print *, '!-------------!'
    print *, '! TEST FAILED !'
    print *, '!-------------!'
  else
    print *, '!-------------!'
    print *, '! TEST PASSED !'
    print *, '!-------------!'
  end if
end program test_ezmatmul
