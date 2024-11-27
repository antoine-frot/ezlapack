module mod_ez_matmul

  implicit none

  interface ez_matmul
  
    module procedure ez_sgemm,                  ez_dgemm
    module procedure ez_cgemm,                  ez_zgemm
    
    module procedure ez_sgemm_trans_alpha_beta, ez_dgemm_trans_alpha_beta
    module procedure ez_cgemm_trans_alpha_beta, ez_zgemm_trans_alpha_beta
    
    module procedure ez_sgemm_alpha_beta,       ez_dgemm_alpha_beta
    module procedure ez_cgemm_alpha_beta,       ez_zgemm_alpha_beta
    
    module procedure ez_sgemm_trans_beta,       ez_dgemm_trans_beta
    module procedure ez_cgemm_trans_beta,       ez_zgemm_trans_beta
    
    module procedure ez_sgemm_trans_alpha,      ez_dgemm_trans_alpha
    module procedure ez_cgemm_trans_alpha,      ez_zgemm_trans_alpha
    
    module procedure ez_sgemm_alpha,            ez_dgemm_alpha
    module procedure ez_cgemm_alpha,            ez_zgemm_alpha
    
    module procedure ez_sgemm_beta,             ez_dgemm_beta
    module procedure ez_cgemm_beta,             ez_zgemm_beta

    module procedure ez_sgemm_trans,            ez_dgemm_trans
    module procedure ez_cgemm_trans,            ez_zgemm_trans
  
  end interface

contains

!------------!
! SGEMM PART !
!------------!

!----------------------------!
! Handling default arguments !
!----------------------------!

  subroutine ez_sgemm_trans_alpha_beta(A, B, C)
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    call ez_sgemm('N', 'N', 1e0, A, B, 0e0, C)
  end subroutine ez_sgemm_trans_alpha_beta

  subroutine ez_sgemm_alpha_beta(transa, transb, A, B, C)
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ez_sgemm(transA, transB, 1e0, A, B, 0e0, C)
  end subroutine ez_sgemm_alpha_beta

  subroutine ez_sgemm_trans_beta(alpha, A, B, C)
    real(4), intent(in)                                :: alpha
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    call ez_sgemm('N', 'N', alpha, A, B, 0e0, C)
  end subroutine ez_sgemm_trans_beta

  subroutine ez_sgemm_trans_alpha(A, B, beta, C)
    real(4), intent(in)                                :: beta
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    call ez_sgemm('N', 'N', 1e0, A, B, beta, C)
  end subroutine ez_sgemm_trans_alpha

  subroutine ez_sgemm_alpha(transa, transb, A, B, beta, C)
    real(4), intent(in)                                :: beta
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                      :: transa, transb
    call ez_sgemm(transA, transB, 1e0, A, B, beta, C)
  end subroutine ez_sgemm_alpha

  subroutine ez_sgemm_beta(transa, transb, alpha, A, B, C)
    real(4), intent(in)                                :: alpha
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ez_sgemm(transA, transB, alpha, A, B, 0e0, C)
  end subroutine ez_sgemm_beta

  subroutine ez_sgemm_trans(alpha, A, B, beta, C)
    real(4), intent(in)                                :: alpha, beta
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    call ez_sgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ez_sgemm_trans

!------------------------------------!
! Checking dimensions and call SGEMM !
!------------------------------------!

  subroutine ez_sgemm(transa, transb, alpha, A, B, beta, C)

    implicit none

    real(4), intent(in)                                :: alpha, beta
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb

    integer                                            :: m, n, k

    ! Check dimensions
    m = size(C,1)
    k = size(C,2)

    if (transa == 'N') then
      n = size(A,2)
      if (transb == 'N') then
        if (n /= size(B,1) .or. m /= size(A,1) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,1) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      end if
    else
      n = size(A,1)
      if (transb == 'N') then
        if (n /= size(B,1) .or. m /= size(A,2) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,2) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      end if
    end if

    call sgemm(transa, transb, m, n, k, alpha, A, m, B, k, beta, C, m)

  end subroutine ez_sgemm


!------------!
! DGEMM PART !
!------------!

!----------------------------!
! Handling default arguments !
!----------------------------!

  subroutine ez_dgemm_trans_alpha_beta(A, B, C)
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    call ez_dgemm('N', 'N', 1d0, A, B, 0d0, C)
  end subroutine ez_dgemm_trans_alpha_beta

  subroutine ez_dgemm_alpha_beta(transa, transb, A, B, C)
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ez_dgemm(transA, transB, 1d0, A, B, 0d0, C)
  end subroutine ez_dgemm_alpha_beta

  subroutine ez_dgemm_trans_beta(alpha, A, B, C)
    real(8), intent(in)                                :: alpha
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    call ez_dgemm('N', 'N', alpha, A, B, 0d0, C)
  end subroutine ez_dgemm_trans_beta

  subroutine ez_dgemm_trans_alpha(A, B, beta, C)
    real(8), intent(in)                                :: beta
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    call ez_dgemm('N', 'N', 1d0, A, B, beta, C)
  end subroutine ez_dgemm_trans_alpha

  subroutine ez_dgemm_alpha(transa, transb, A, B, beta, C)
    real(8), intent(in)                                :: beta
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                      :: transa, transb
    call ez_dgemm(transA, transB, 1d0, A, B, beta, C)
  end subroutine ez_dgemm_alpha

  subroutine ez_dgemm_beta(transa, transb, alpha, A, B, C)
    real(8), intent(in)                                :: alpha
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ez_dgemm(transA, transB, alpha, A, B, 0d0, C)
  end subroutine ez_dgemm_beta

  subroutine ez_dgemm_trans(alpha, A, B, beta, C)
    real(8), intent(in)                                :: alpha, beta
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    call ez_dgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ez_dgemm_trans

!------------------------------------!
! Checking dimensions and call DGEMM !
!------------------------------------!

  subroutine ez_dgemm(transa, transb, alpha, A, B, beta, C)

    implicit none

    real(8), intent(in)                                :: alpha, beta
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb

    integer                                            :: m, n, k
    
    ! Check dimensions
    m = size(C,1)
    k = size(C,2)

    if (transa == 'N') then
      n = size(A,2)
      if (transb == 'N') then
        if (n /= size(B,1) .or. m /= size(A,1) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,1) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      end if
    else
      n = size(A,1)
      if (transb == 'N') then
        if (n /= size(B,1) .or. m /= size(A,2) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,2) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      end if
    end if

    call dgemm(transa, transb, m, n, k, alpha, A, m, B, k, beta, C, m)

  end subroutine ez_dgemm

!------------!
! CGEMM PART !
!------------!

!----------------------------!
! Handling default arguments !
!----------------------------!

  subroutine ez_cgemm_trans_alpha_beta(A, B, C)
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    call ez_cgemm('N', 'N', (1e0,0e0), A, B, (0e0,0e0), C)
  end subroutine ez_cgemm_trans_alpha_beta

  subroutine ez_cgemm_alpha_beta(transa, transb, A, B, C)
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ez_cgemm(transA, transB, (1e0,0e0), A, B, (0e0,0e0), C)
  end subroutine ez_cgemm_alpha_beta

  subroutine ez_cgemm_trans_beta(alpha, A, B, C)
    complex, intent(in)                                :: alpha
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    call ez_cgemm('N', 'N', alpha, A, B, (0e0,0e0), C)
  end subroutine ez_cgemm_trans_beta

  subroutine ez_cgemm_trans_alpha(A, B, beta, C)
    complex, intent(in)                                :: beta
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    call ez_cgemm('N', 'N', (1e0,0e0), A, B, beta, C)
  end subroutine ez_cgemm_trans_alpha

  subroutine ez_cgemm_alpha(transa, transb, A, B, beta, C)
    complex, intent(in)                                :: beta
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                      :: transa, transb
    call ez_cgemm(transA, transB, (1e0,0e0), A, B, beta, C)
  end subroutine ez_cgemm_alpha

  subroutine ez_cgemm_beta(transa, transb, alpha, A, B, C)
    complex, intent(in)                                :: alpha
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ez_cgemm(transA, transB, alpha, A, B, (0e0,0e0), C)
  end subroutine ez_cgemm_beta

  subroutine ez_cgemm_trans(alpha, A, B, beta, C)
    complex, intent(in)                                :: alpha, beta
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    call ez_cgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ez_cgemm_trans

!------------------------------------!
! Checking dimensions and call CGEMM !
!------------------------------------!

  subroutine ez_cgemm(transa, transb, alpha, A, B, beta, C)

    implicit none

    complex, intent(in)                                :: alpha, beta
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb

    integer                                            :: m, n, k
     
    ! Check dimensions
    m = size(C,1)
    k = size(C,2)

    if (transa == 'N') then
      n = size(A,2)
      if (transb == 'N') then
        if (n /= size(B,1) .or. m /= size(A,1) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,1) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      end if
    else
      n = size(A,1)
      if (transb == 'N') then
        if (n /= size(B,1) .or. m /= size(A,2) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,2) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      end if
    end if

    call cgemm(transa, transb, m, n, k, alpha, A, m, B, k, beta, C, m)

  end subroutine ez_cgemm

!------------!
! ZGEMM PART !
!------------!

!----------------------------!
! Handling default arguments !
!----------------------------!

  subroutine ez_zgemm_trans_alpha_beta(A, B, C)
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    call ez_zgemm('N', 'N', (1d0,0d0), A, B, (0d0,0d0), C)
  end subroutine ez_zgemm_trans_alpha_beta

  subroutine ez_zgemm_alpha_beta(transa, transb, A, B, C)
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ez_zgemm(transA, transB, (1d0,0d0), A, B, (0d0,0d0), C)
  end subroutine ez_zgemm_alpha_beta

  subroutine ez_zgemm_trans_beta(alpha, A, B, C)
    complex*16, intent(in)                                :: alpha
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    call ez_zgemm('N', 'N', alpha, A, B, (0d0,0d0), C)
  end subroutine ez_zgemm_trans_beta

  subroutine ez_zgemm_trans_alpha(A, B, beta, C)
    complex*16, intent(in)                                :: beta
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    call ez_zgemm('N', 'N', (1d0,0d0), A, B, beta, C)
  end subroutine ez_zgemm_trans_alpha

  subroutine ez_zgemm_alpha(transa, transb, A, B, beta, C)
    complex*16, intent(in)                                :: beta
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                      :: transa, transb
    call ez_zgemm(transA, transB, (1d0,0d0), A, B, beta, C)
  end subroutine ez_zgemm_alpha

  subroutine ez_zgemm_beta(transa, transb, alpha, A, B, C)
    complex*16, intent(in)                                :: alpha
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ez_zgemm(transA, transB, alpha, A, B, (0d0,0d0), C)
  end subroutine ez_zgemm_beta

  subroutine ez_zgemm_trans(alpha, A, B, beta, C)
    complex*16, intent(in)                                :: alpha, beta
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    call ez_zgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ez_zgemm_trans

!------------------------------------!
! Checking dimensions and call ZGEMM !
!------------------------------------!

  subroutine ez_zgemm(transa, transb, alpha, A, B, beta, C)

    implicit none

    complex*16, intent(in)                                :: alpha, beta
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb

    integer                                            :: m, n, k
     
    ! Check dimensions
    m = size(C,1)
    k = size(C,2)

    if (transa == 'N') then
      n = size(A,2)
      if (transb == 'N') then
        if (n /= size(B,1) .or. m /= size(A,1) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,1) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      end if
    else
      n = size(A,1)
      if (transb == 'N') then
        if (n /= size(B,1) .or. m /= size(A,2) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,2) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "m: ", m, "n: ", n, "k: ", k
          stop
        end if
      end if
    end if

    call zgemm(transa, transb, m, n, k, alpha, A, m, B, k, beta, C, m)

  end subroutine ez_zgemm
  
end module mod_ez_matmul
