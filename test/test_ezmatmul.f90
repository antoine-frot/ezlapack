program test_ezmatmul
  ! Test file for ezmatmul
  !
  ! Tests every optional argument separately for each LAPACK matrix multicplication subroutine (sgemm, dgemm, cgemm, zgemm)
  ! For each case add a random test with random condition for every argument

  use ezlapack

  implicit none

  logical                              :: test_passed
  character(len=1)                     :: transa, transb
  integer                              :: m, n, k

  real(4)                              :: alpha_sp, beta_sp
  real(4), allocatable, dimension(:,:) :: A_sp, B_sp, C_sp, D_sp
  real(8)                              :: alpha_dp, beta_dp
  real(8), allocatable, dimension(:,:) :: A_dp, B_dp, C_dp, D_dp
  complex(4)                              :: alpha_csp, beta_csp
  complex(4), allocatable, dimension(:,:) :: A_csp, B_csp, C_csp, D_csp
  complex(8)                              :: alpha_cdp, beta_cdp
  complex(8), allocatable, dimension(:,:) :: A_cdp, B_cdp, C_cdp, D_cdp
  

  call random_integer(1, 2000, m)
  call random_integer(1, 2000, n)
  call random_integer(1, 2000, k)
  write(*, '(A, "A: (", I5, ", ", I5, ")")') "Size of the matrices: ", m, k
  write(*, '(A, "B: (", I5, ", ", I5, ")")') "                      ", k, n
  write(*, '(A, "C: (", I5, ", ", I5, ")")') "                      ", m, n

!------------!
! SGEMM PART !
!------------!

!----------------------------!
! Handling default arguments !
!----------------------------!

  print *, "Test for sgemm"
  test_passed = .true.
  allocate(A_sp(m, k), B_sp(k, n), C_sp(m, n), D_sp(m, n))
  call random_number(A_sp)
  call random_number(B_sp)

  ! Test of all default arguments
  call sgemm('N', 'N', m, n, k, 1e0, A_sp, m, B_sp, k, 0e0, D_sp, m)

  call ezmatmul(A_sp, B_sp, C_sp)
  if (any(D_sp /= C_sp)) then
    print *, "Test failed for ezsgemm_trans_alpha_beta"
    print *, 'Maximal difference between sgemm and ezmatmul: ', maxval(abs(D_sp - C_sp))
    test_passed = .false.
  end if

  call ezmatmul('N', 'N', A_sp, B_sp, C_sp)
  if (any(D_sp /= C_sp)) then
    print *, "Test failed for ezsgemm_alpha_beta"
    print *, 'Maximal difference between sgemm and ezmatmul: ', maxval(abs(D_sp - C_sp))
    test_passed = .false.
  end if

  call ezmatmul(1e0, A_sp, B_sp, C_sp)
  if (any(D_sp /= C_sp)) then
    print *, "Test failed for ezsgemm_trans_beta"
    print *, 'Maximal difference between sgemm and ezmatmul: ', maxval(abs(D_sp - C_sp))
    test_passed = .false.
  end if

  call ezmatmul(A_sp, B_sp, 0e0, C_sp)
  if (any(D_sp /= C_sp)) then
    print *, "Test failed for ezsgemm_trans_alpha"
    print *, 'Maximal difference between sgemm and ezmatmul: ', maxval(abs(D_sp - C_sp))
    test_passed = .false.
  end if

  call ezmatmul('N', 'N', A_sp, B_sp, 0e0, C_sp)
  if (any(D_sp /= C_sp)) then
    print *, "Test failed for ezsgemm_alpha"
    print *, 'Maximal difference between sgemm and ezmatmul: ', maxval(abs(D_sp - C_sp))
    test_passed = .false.
  end if

  call ezmatmul('N', 'N', 1e0, A_sp, B_sp, C_sp)
  if (any(D_sp /= C_sp)) then
    print *, "Test failed for ezsgemm_beta"
    print *, 'Maximal difference between sgemm and ezmatmul: ', maxval(abs(D_sp - C_sp))
    test_passed = .false.
  end if

  call ezmatmul(1e0, A_sp, B_sp, 0e0, C_sp)
  if (any(D_sp /= C_sp)) then
    print *, "Test failed for ezsgemm_trans"
    print *, 'Maximal difference between sgemm and ezmatmul: ', maxval(abs(D_sp - C_sp))
    test_passed = .false.
  end if
  
  call ezmatmul('N', 'N', 1e0, A_sp, B_sp, 0e0, C_sp)
  if (any(D_sp /= C_sp)) then
    print *, "Test failed for ezsgemm"
    print *, 'Maximal difference between sgemm and ezmatmul: ', maxval(abs(D_sp - C_sp))
    test_passed = .false.
  end if
  
!--------------------------------!
! Test with all arguments random !
!--------------------------------!

  call random_character((/ 'N', 'T', 'C'/), transa)
  call random_character((/ 'N', 'T', 'C'/), transb)

  call random_number(alpha_sp)
  call random_number(beta_sp)
  alpha_sp = (alpha_sp - 0.5 ) * 10000
  beta_sp  = (beta_sp  - 0.5 ) * 10000
  call random_number(C_sp)
  D_sp = C_sp

  if (transa == 'N') then
    if (transb == 'N') then
      call sgemm(transa, transb, m, n, k, alpha_sp, A_sp, m, B_sp, k, beta_sp, D_sp, m)
    else
      deallocate(B_sp)
      allocate(B_sp(n,k))
      call sgemm(transa, transb, m, n, k, alpha_sp, A_sp, m, B_sp, n, beta_sp, D_sp, m)
    end if
  else 
    deallocate(A_sp)
    allocate(A_sp(k,m))
    if (transb == 'N') then
      call sgemm(transa, transb, m, n, k, alpha_sp, A_sp, k, B_sp, k, beta_sp, D_sp, m)
    else
      deallocate(B_sp)
      allocate(B_sp(n,k))
      call sgemm(transa, transb, m, n, k, alpha_sp, A_sp, k, B_sp, n, beta_sp, D_sp, m)
    end if
  end if

  call ezmatmul(transa, transb, alpha_sp, A_sp, B_sp, beta_sp, C_sp)

  if (any(D_sp /= C_sp)) then
    print *, "Random test with all dummy arguments failed for sgemm"
    print *, 'Maximal difference between sgemm and ezmatmul: ', maxval(abs(D_sp - C_sp))
    print *, "transa: ", trim(transa)
    print *, "transb: ", trim(transb)
    print *, "alpha_sp:", alpha_sp
    print *, "beta_sp: ", beta_sp
    test_passed = .false.
  end if

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

!------------!
! CGEMM PART !
!------------!

!----------------------------!
! Handling default arguments !
!----------------------------!

  print *, "Test for cgemm"
  test_passed = .true.
  allocate(A_csp(m, k), B_csp(k, n), C_csp(m, n), D_csp(m, n))
  call random_complex(A_csp)
  call random_complex(B_csp)

  ! Test of all default arguments
  call cgemm('N', 'N', m, n, k, (1e0, 0e0), A_csp, m, B_csp, k, (0e0, 0e0), D_csp, m)

  call ezmatmul(A_csp, B_csp, C_csp)
  if (any(D_csp /= C_csp)) then
    print *, "Test failed for ezcgemm_trans_alpha_beta"
    print *, 'Maximal difference between cgemm and ezmatmul: ', maxval(abs(D_csp - C_csp))
    test_passed = .false.
  end if

  call ezmatmul('N', 'N', A_csp, B_csp, C_csp)
  if (any(D_csp /= C_csp)) then
    print *, "Test failed for ezcgemm_alpha_beta"
    print *, 'Maximal difference between cgemm and ezmatmul: ', maxval(abs(D_csp - C_csp))
    test_passed = .false.
  end if

  call ezmatmul((1e0, 0e0), A_csp, B_csp, C_csp)
  if (any(D_csp /= C_csp)) then
    print *, "Test failed for ezcgemm_trans_beta"
    print *, 'Maximal difference between cgemm and ezmatmul: ', maxval(abs(D_csp - C_csp))
    test_passed = .false.
  end if

  call ezmatmul(A_csp, B_csp, (0e0, 0e0), C_csp)
  if (any(D_csp /= C_csp)) then
    print *, "Test failed for ezcgemm_trans_alpha"
    print *, 'Maximal difference between cgemm and ezmatmul: ', maxval(abs(D_csp - C_csp))
    test_passed = .false.
  end if

  call ezmatmul('N', 'N', A_csp, B_csp, (0e0, 0e0), C_csp)
  if (any(D_csp /= C_csp)) then
    print *, "Test failed for ezcgemm_alpha"
    print *, 'Maximal difference between cgemm and ezmatmul: ', maxval(abs(D_csp - C_csp))
    test_passed = .false.
  end if

  call ezmatmul('N', 'N', (1e0, 0e0), A_csp, B_csp, C_csp)
  if (any(D_csp /= C_csp)) then
    print *, "Test failed for ezcgemm_beta"
    print *, 'Maximal difference between cgemm and ezmatmul: ', maxval(abs(D_csp - C_csp))
    test_passed = .false.
  end if

  call ezmatmul((1e0, 0e0), A_csp, B_csp, (0e0, 0e0), C_csp)
  if (any(D_csp /= C_csp)) then
    print *, "Test failed for ezcgemm_trans"
    print *, 'Maximal difference between cgemm and ezmatmul: ', maxval(abs(D_csp - C_csp))
    test_passed = .false.
  end if
  
  call ezmatmul('N', 'N', (1e0, 0e0), A_csp, B_csp, (0e0, 0e0), C_csp)
  if (any(D_csp /= C_csp)) then
    print *, "Test failed for ezcgemm"
    print *, 'Maximal difference between cgemm and ezmatmul: ', maxval(abs(D_csp - C_csp))
    test_passed = .false.
  end if
  
!--------------------------------!
! Test with all arguments random !
!--------------------------------!

  call random_character((/ 'N', 'T', 'C'/), transa)
  call random_character((/ 'N', 'T', 'C'/), transb)

  call random_complex(alpha_csp)
  call random_complex(beta_csp)
  alpha_csp = alpha_csp*10000
  beta_csp = beta_csp*10000
  call random_complex(C_csp)
  D_csp = C_csp

  if (transa == 'N') then
    if (transb == 'N') then
      call cgemm(transa, transb, m, n, k, alpha_csp, A_csp, m, B_csp, k, beta_csp, D_csp, m)
    else
      deallocate(B_csp)
      allocate(B_csp(n,k))
      call cgemm(transa, transb, m, n, k, alpha_csp, A_csp, m, B_csp, n, beta_csp, D_csp, m)
    end if
  else 
    deallocate(A_csp)
    allocate(A_csp(k,m))
    if (transb == 'N') then
      call cgemm(transa, transb, m, n, k, alpha_csp, A_csp, k, B_csp, k, beta_csp, D_csp, m)
    else
      deallocate(B_csp)
      allocate(B_csp(n,k))
      call cgemm(transa, transb, m, n, k, alpha_csp, A_csp, k, B_csp, n, beta_csp, D_csp, m)
    end if
  end if

  call ezmatmul(transa, transb, alpha_csp, A_csp, B_csp, beta_csp, C_csp)

  if (any(D_csp /= C_csp)) then
    print *, "Random test with all dummy arguments failed for cgemm"
    print *, 'Maximal difference between cgemm and ezmatmul: ', maxval(abs(D_csp - C_csp))
    print *, "transa: ", trim(transa)
    print *, "transb: ", trim(transb)
    print *, "alpha_csp:", alpha_csp
    print *, "beta_csp: ", beta_csp
    test_passed = .false.
  end if

!------------!
! ZGEMM PART !
!------------!

!----------------------------!
! Handling default arguments !
!----------------------------!

  print *, "Test for zgemm"
  test_passed = .true.
  allocate(A_cdp(m, k), B_cdp(k, n), C_cdp(m, n), D_cdp(m, n))
  call random_complex(A_cdp)
  call random_complex(B_cdp)

  ! Test of all default arguments
  call zgemm('N', 'N', m, n, k, (1d0, 0e0), A_cdp, m, B_cdp, k, (0d0, 0e0), D_cdp, m)

  call ezmatmul(A_cdp, B_cdp, C_cdp)
  if (any(D_cdp /= C_cdp)) then
    print *, "Test failed for ezzgemm_trans_alpha_beta"
    print *, 'Maximal difference between zgemm and ezmatmul: ', maxval(abs(D_cdp - C_cdp))
    test_passed = .false.
  end if

  call ezmatmul('N', 'N', A_cdp, B_cdp, C_cdp)
  if (any(D_cdp /= C_cdp)) then
    print *, "Test failed for ezzgemm_alpha_beta"
    print *, 'Maximal difference between zgemm and ezmatmul: ', maxval(abs(D_cdp - C_cdp))
    test_passed = .false.
  end if

  call ezmatmul((1d0, 0e0), A_cdp, B_cdp, C_cdp)
  if (any(D_cdp /= C_cdp)) then
    print *, "Test failed for ezzgemm_trans_beta"
    print *, 'Maximal difference between zgemm and ezmatmul: ', maxval(abs(D_cdp - C_cdp))
    test_passed = .false.
  end if

  call ezmatmul(A_cdp, B_cdp, (0d0, 0e0), C_cdp)
  if (any(D_cdp /= C_cdp)) then
    print *, "Test failed for ezzgemm_trans_alpha"
    print *, 'Maximal difference between zgemm and ezmatmul: ', maxval(abs(D_cdp - C_cdp))
    test_passed = .false.
  end if

  call ezmatmul('N', 'N', A_cdp, B_cdp, (0d0, 0e0), C_cdp)
  if (any(D_cdp /= C_cdp)) then
    print *, "Test failed for ezzgemm_alpha"
    print *, 'Maximal difference between zgemm and ezmatmul: ', maxval(abs(D_cdp - C_cdp))
    test_passed = .false.
  end if

  call ezmatmul('N', 'N', (1d0, 0e0), A_cdp, B_cdp, C_cdp)
  if (any(D_cdp /= C_cdp)) then
    print *, "Test failed for ezzgemm_beta"
    print *, 'Maximal difference between zgemm and ezmatmul: ', maxval(abs(D_cdp - C_cdp))
    test_passed = .false.
  end if

  call ezmatmul((1d0, 0e0), A_cdp, B_cdp, (0d0, 0e0), C_cdp)
  if (any(D_cdp /= C_cdp)) then
    print *, "Test failed for ezzgemm_trans"
    print *, 'Maximal difference between zgemm and ezmatmul: ', maxval(abs(D_cdp - C_cdp))
    test_passed = .false.
  end if
  
  call ezmatmul('N', 'N', (1d0, 0e0), A_cdp, B_cdp, (0d0, 0e0), C_cdp)
  if (any(D_cdp /= C_cdp)) then
    print *, "Test failed for ezzgemm"
    print *, 'Maximal difference between zgemm and ezmatmul: ', maxval(abs(D_cdp - C_cdp))
    test_passed = .false.
  end if
  
!--------------------------------!
! Test with all arguments random !
!--------------------------------!

  call random_character((/ 'N', 'T', 'C'/), transa)
  call random_character((/ 'N', 'T', 'C'/), transb)

  call random_complex(alpha_cdp)
  call random_complex(beta_cdp)
  alpha_cdp = alpha_cdp*10000
  beta_cdp  = beta_cdp*10000
  call random_complex(C_cdp)
  D_cdp = C_cdp

  if (transa == 'N') then
    if (transb == 'N') then
      call zgemm(transa, transb, m, n, k, alpha_cdp, A_cdp, m, B_cdp, k, beta_cdp, D_cdp, m)
    else
      deallocate(B_cdp)
      allocate(B_cdp(n,k))
      call zgemm(transa, transb, m, n, k, alpha_cdp, A_cdp, m, B_cdp, n, beta_cdp, D_cdp, m)
    end if
  else 
    deallocate(A_cdp)
    allocate(A_cdp(k,m))
    if (transb == 'N') then
      call zgemm(transa, transb, m, n, k, alpha_cdp, A_cdp, k, B_cdp, k, beta_cdp, D_cdp, m)
    else
      deallocate(B_cdp)
      allocate(B_cdp(n,k))
      call zgemm(transa, transb, m, n, k, alpha_cdp, A_cdp, k, B_cdp, n, beta_cdp, D_cdp, m)
    end if
  end if

  call ezmatmul(transa, transb, alpha_cdp, A_cdp, B_cdp, beta_cdp, C_cdp)

  if (any(D_cdp /= C_cdp)) then
    print *, "Random test with all dummy arguments failed for zgemm"
    print *, 'Maximal difference between zgemm and ezmatmul: ', maxval(abs(D_cdp - C_cdp))
    print *, "transa: ", trim(transa)
    print *, "transb: ", trim(transb)
    print *, "alpha_cdp:", alpha_cdp
    print *, "beta_cdp: ", beta_cdp
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
