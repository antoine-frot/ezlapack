module module_ezmatmul
!--------------------------------------------------------------------------------!
! Purpose:                                                                       ! 
! Generic interface for lapack matrice multiplication (*gemm)                    ! 
! i.e. performs C := alpha*op( A )*op( B ) + beta*C                              !     
! where alpha and beta are scalars, A, B and C are matrices, and                 ! 
! op( X ) is one of op( X ) = X or op( X ) = X**T or op( X ) = X**H.             !   
!--------------------------------------------------------------------------------!       
! Subroutine:                                                                    !         
! Type can be real(4), real(8), complex or complex*16,                           !           
! but all dummy arguments should have the same type.                             !   
!                                                                                !   
! subroutine ezmatmul	(	character*1, optional (default = 'N') :: transa,         !     
!                       character*1, optional (default = 'N') :: transb,         !     
!                       type,        optional (default =  1 ) :: alpha,          !       
!                       type, dimension(:,:), contiguous      :: A,              !       
!                       type, dimension(:,:), contiguous      :: B,              !     
!                       type,        optional (default =  0 ) :: beta,           !         
!                       type, dimension(:,:), contiguous      :: C,              ! 
! )                                                                              !  
!                                                                                ! 
! transa (character*1, optional (default = 'N')): specifies the form of op( A ). !
!   if transa = 'N' or 'n',  op( A ) = A.                                        !
!   if transa = 'T' or 't',  op( A ) = A**T.                                     !
!   if transa = 'C' or 'c',  op( A ) = A**H.                                     !
!                                                                                !
! transb (character*1, optional (default = 'N')): specifies the form of op( B ). !
!   if transb = 'N' or 'n',  op( B ) = B.                                        !
!   if transb = 'T' or 't',  op( B ) = B**T.                                     !
!   if transb = 'C' or 'c',  op( B ) = B**H.                                     !
!                                                                                !
! alpha (type, optional (default = 1)): specifies the scalar alpha.              !
!                                                                                !
! A (type, dimension(:), contiguous): specifies the matrix A                     !
!                                                                                !
! B (type, dimension(:), contiguous): specifies the matrix B                     !
!                                                                                !
! beta (type, optional (default = 0)): specifies the scalar beta.                !
!                                                                                !
! C (type, dimension(:), contiguous): specifies the matrix C                     !
!--------------------------------------------------------------------------------!
! Examples:                                                                      !
!                                                                                !
! call ezmatmul(A,B,C) computes C:= A*B with the simplicity of matmul            !
! while harnessing the high performance of the LAPACK library.                   !
!                                                                                !
! On the other hand, ezmatmul keep the versatility of gemm with for example:     !   
! call ezmatmul('C','T',(5d-4,8d6),A,B,(9d2,3d-7),C)                             !     
!--------------------------------------------------------------------------------!   
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64

  implicit none

  interface ezmatmul
    ! Generic interface
    module procedure ezsgemm,                  ezdgemm
    module procedure ezcgemm,                  ezzgemm
    
    module procedure ezsgemm_trans_alpha_beta, ezdgemm_trans_alpha_beta
    module procedure ezcgemm_trans_alpha_beta, ezzgemm_trans_alpha_beta
    
    module procedure ezsgemm_alpha_beta,       ezdgemm_alpha_beta
    module procedure ezcgemm_alpha_beta,       ezzgemm_alpha_beta
    
    module procedure ezsgemm_trans_beta,       ezdgemm_trans_beta
    module procedure ezcgemm_trans_beta,       ezzgemm_trans_beta
    
    module procedure ezsgemm_trans_alpha,      ezdgemm_trans_alpha
    module procedure ezcgemm_trans_alpha,      ezzgemm_trans_alpha
    
    module procedure ezsgemm_alpha,            ezdgemm_alpha
    module procedure ezcgemm_alpha,            ezzgemm_alpha
    
    module procedure ezsgemm_beta,             ezdgemm_beta
    module procedure ezcgemm_beta,             ezzgemm_beta

    module procedure ezsgemm_trans,            ezdgemm_trans
    module procedure ezcgemm_trans,            ezzgemm_trans
  
  end interface

contains

!------------!
! SGEMM PART !
!------------!

!----------------------------!
! Handling default arguments !
!----------------------------!

  subroutine ezsgemm_trans_alpha_beta(A, B, C)
    ! Calls the SGEMM routines with default values for 
    ! transa, transb, alpha and beta ('N', 'N', 1e0 and 0e0, respectively).
    real(sp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(sp), dimension(:,:), contiguous, intent(inout) :: C
    call ezsgemm('N', 'N', 1e0, A, B, 0e0, C)
  end subroutine ezsgemm_trans_alpha_beta

  subroutine ezsgemm_alpha_beta(transa, transb, A, B, C)
    ! Calls the SGEMM routines with default values for 
    ! alpha and beta (1e0 and 0e0, respectively).
    real(sp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(sp), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezsgemm(transA, transB, 1e0, A, B, 0e0, C)
  end subroutine ezsgemm_alpha_beta

  subroutine ezsgemm_trans_beta(alpha, A, B, C)
    ! Calls the SGEMM routines with default values for 
    ! transa, transb and beta ('N', 'N' and 0e0, respectively).
    real(sp), intent(in)                                :: alpha
    real(sp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(sp), dimension(:,:), contiguous, intent(inout) :: C
    call ezsgemm('N', 'N', alpha, A, B, 0e0, C)
  end subroutine ezsgemm_trans_beta

  subroutine ezsgemm_trans_alpha(A, B, beta, C)
    ! Calls the SGEMM routines with default values for 
    ! transa, transb and alpha  ('N', 'N' and 1e0, respectively).
    real(sp), intent(in)                                :: beta
    real(sp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(sp), dimension(:,:), contiguous, intent(inout) :: C
    call ezsgemm('N', 'N', 1e0, A, B, beta, C)
  end subroutine ezsgemm_trans_alpha

  subroutine ezsgemm_alpha(transa, transb, A, B, beta, C)
    ! Calls the SGEMM routines with default values for 
    ! alpha (1e0).
    real(sp), intent(in)                                :: beta
    real(sp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(sp), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezsgemm(transA, transB, 1e0, A, B, beta, C)
  end subroutine ezsgemm_alpha

  subroutine ezsgemm_beta(transa, transb, alpha, A, B, C)
    ! Calls the SGEMM routines with default values for 
    ! beta (0e0).
    real(sp), intent(in)                                :: alpha
    real(sp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(sp), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezsgemm(transA, transB, alpha, A, B, 0e0, C)
  end subroutine ezsgemm_beta

  subroutine ezsgemm_trans(alpha, A, B, beta, C)
    ! Calls the SGEMM routines with default values for 
    ! transa and transb ('N' and 'N', respectively).
    real(sp), intent(in)                                :: alpha, beta
    real(sp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(sp), dimension(:,:), contiguous, intent(inout) :: C
    call ezsgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ezsgemm_trans

!------------------------------------!
! Checking dimensions and call SGEMM !
!------------------------------------!

  subroutine ezsgemm(transa, transb, alpha, A, B, beta, C)

    real(sp), intent(in)                                :: alpha, beta
    real(sp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(sp), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                        :: transa, transb

    integer                                             :: m, n, k

    ! Check dimensions
    m = size(C,1)
    k = size(C,2)

    if (transa == 'N') then
      n = size(A,2)
      if (transb == 'N') then
        if (n /= size(B,1) .or. m /= size(A,1) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,1) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      end if
    else
      n = size(A,1)
      if (transb == 'N') then
        if (n /= size(B,1) .or. m /= size(A,2) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,2) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      end if
    end if

    call sgemm(transa, transb, m, n, k, alpha, A, m, B, k, beta, C, m)

  end subroutine ezsgemm


!------------!
! DGEMM PART !
!------------!

!----------------------------!
! Handling default arguments !
!----------------------------!

  subroutine ezdgemm_trans_alpha_beta(A, B, C)
    real(dp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(dp), dimension(:,:), contiguous, intent(inout) :: C
    call ezdgemm('N', 'N', 1d0, A, B, 0d0, C)
  end subroutine ezdgemm_trans_alpha_beta

  subroutine ezdgemm_alpha_beta(transa, transb, A, B, C)
    real(dp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(dp), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezdgemm(transA, transB, 1d0, A, B, 0d0, C)
  end subroutine ezdgemm_alpha_beta

  subroutine ezdgemm_trans_beta(alpha, A, B, C)
    real(dp), intent(in)                                :: alpha
    real(dp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(dp), dimension(:,:), contiguous, intent(inout) :: C
    call ezdgemm('N', 'N', alpha, A, B, 0d0, C)
  end subroutine ezdgemm_trans_beta

  subroutine ezdgemm_trans_alpha(A, B, beta, C)
    real(dp), intent(in)                                :: beta
    real(dp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(dp), dimension(:,:), contiguous, intent(inout) :: C
    call ezdgemm('N', 'N', 1d0, A, B, beta, C)
  end subroutine ezdgemm_trans_alpha

  subroutine ezdgemm_alpha(transa, transb, A, B, beta, C)
    real(dp), intent(in)                                :: beta
    real(dp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(dp), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezdgemm(transA, transB, 1d0, A, B, beta, C)
  end subroutine ezdgemm_alpha

  subroutine ezdgemm_beta(transa, transb, alpha, A, B, C)
    real(dp), intent(in)                                :: alpha
    real(dp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(dp), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezdgemm(transA, transB, alpha, A, B, 0d0, C)
  end subroutine ezdgemm_beta

  subroutine ezdgemm_trans(alpha, A, B, beta, C)
    real(dp), intent(in)                                :: alpha, beta
    real(dp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(dp), dimension(:,:), contiguous, intent(inout) :: C
    call ezdgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ezdgemm_trans

!------------------------------------!
! Checking dimensions and call DGEMM !
!------------------------------------!

  subroutine ezdgemm(transa, transb, alpha, A, B, beta, C)

    real(dp), intent(in)                                :: alpha, beta
    real(dp), dimension(:,:), contiguous, intent(in)    :: A, B
    real(dp), dimension(:,:), contiguous, intent(inout) :: C
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
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,1) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      end if
    else
      n = size(A,1)
      if (transb == 'N') then
        if (n /= size(B,1) .or. m /= size(A,2) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,2) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      end if
    end if

    call dgemm(transa, transb, m, n, k, alpha, A, m, B, k, beta, C, m)

  end subroutine ezdgemm

!------------!
! CGEMM PART !
!------------!

!----------------------------!
! Handling default arguments !
!----------------------------!

  subroutine ezcgemm_trans_alpha_beta(A, B, C)
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    call ezcgemm('N', 'N', (1e0,0e0), A, B, (0e0,0e0), C)
  end subroutine ezcgemm_trans_alpha_beta

  subroutine ezcgemm_alpha_beta(transa, transb, A, B, C)
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezcgemm(transA, transB, (1e0,0e0), A, B, (0e0,0e0), C)
  end subroutine ezcgemm_alpha_beta

  subroutine ezcgemm_trans_beta(alpha, A, B, C)
    complex, intent(in)                                :: alpha
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    call ezcgemm('N', 'N', alpha, A, B, (0e0,0e0), C)
  end subroutine ezcgemm_trans_beta

  subroutine ezcgemm_trans_alpha(A, B, beta, C)
    complex, intent(in)                                :: beta
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    call ezcgemm('N', 'N', (1e0,0e0), A, B, beta, C)
  end subroutine ezcgemm_trans_alpha

  subroutine ezcgemm_alpha(transa, transb, A, B, beta, C)
    complex, intent(in)                                :: beta
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezcgemm(transA, transB, (1e0,0e0), A, B, beta, C)
  end subroutine ezcgemm_alpha

  subroutine ezcgemm_beta(transa, transb, alpha, A, B, C)
    complex, intent(in)                                :: alpha
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezcgemm(transA, transB, alpha, A, B, (0e0,0e0), C)
  end subroutine ezcgemm_beta

  subroutine ezcgemm_trans(alpha, A, B, beta, C)
    complex, intent(in)                                :: alpha, beta
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    call ezcgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ezcgemm_trans

!------------------------------------!
! Checking dimensions and call CGEMM !
!------------------------------------!

  subroutine ezcgemm(transa, transb, alpha, A, B, beta, C)

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
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,1) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      end if
    else
      n = size(A,1)
      if (transb == 'N') then
        if (n /= size(B,1) .or. m /= size(A,2) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,2) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      end if
    end if

    call cgemm(transa, transb, m, n, k, alpha, A, m, B, k, beta, C, m)

  end subroutine ezcgemm

!------------!
! ZGEMM PART !
!------------!

!----------------------------!
! Handling default arguments !
!----------------------------!

  subroutine ezzgemm_trans_alpha_beta(A, B, C)
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    call ezzgemm('N', 'N', (1d0,0d0), A, B, (0d0,0d0), C)
  end subroutine ezzgemm_trans_alpha_beta

  subroutine ezzgemm_alpha_beta(transa, transb, A, B, C)
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                          :: transa, transb
    call ezzgemm(transA, transB, (1d0,0d0), A, B, (0d0,0d0), C)
  end subroutine ezzgemm_alpha_beta

  subroutine ezzgemm_trans_beta(alpha, A, B, C)
    complex*16, intent(in)                                :: alpha
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    call ezzgemm('N', 'N', alpha, A, B, (0d0,0d0), C)
  end subroutine ezzgemm_trans_beta

  subroutine ezzgemm_trans_alpha(A, B, beta, C)
    complex*16, intent(in)                                :: beta
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    call ezzgemm('N', 'N', (1d0,0d0), A, B, beta, C)
  end subroutine ezzgemm_trans_alpha

  subroutine ezzgemm_alpha(transa, transb, A, B, beta, C)
    complex*16, intent(in)                                :: beta
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                      :: transa, transb
    call ezzgemm(transA, transB, (1d0,0d0), A, B, beta, C)
  end subroutine ezzgemm_alpha

  subroutine ezzgemm_beta(transa, transb, alpha, A, B, C)
    complex*16, intent(in)                                :: alpha
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezzgemm(transA, transB, alpha, A, B, (0d0,0d0), C)
  end subroutine ezzgemm_beta

  subroutine ezzgemm_trans(alpha, A, B, beta, C)
    complex*16, intent(in)                                :: alpha, beta
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    call ezzgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ezzgemm_trans

!------------------------------------!
! Checking dimensions and call ZGEMM !
!------------------------------------!

  subroutine ezzgemm(transa, transb, alpha, A, B, beta, C)

    complex*16, intent(in)                                :: alpha, beta
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                          :: transa, transb

    integer                                               :: m, n, k
     
    ! Check dimensions
    m = size(C,1)
    k = size(C,2)

    if (transa == 'N') then
      n = size(A,2)
      if (transb == 'N') then
        if (n /= size(B,1) .or. m /= size(A,1) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,1) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      end if
    else
      n = size(A,1)
      if (transb == 'N') then
        if (n /= size(B,1) .or. m /= size(A,2) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,2) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          print *, "size(C,1): ", m, "size(A,2): ", n, "size(C,2): ", k
          stop
        end if
      end if
    end if

    call zgemm(transa, transb, m, n, k, alpha, A, m, B, k, beta, C, m)

  end subroutine ezzgemm
  
end module module_ezmatmul
