module module_ezmatmul
!--------------------------------------------------------------------------------
! Purpose:                                                                      
! User-friendly generic interface for lapack matrice multiplication (*gemm)     
! i.e. performs C := alpha*op( A )*op( B ) + beta*C                             
! where:                                                                        
! - alpha and beta are scalars.                                                 
! - A, B, and C are matrices.                                                   
! - op(X) can be X, X**T, or X**H.                                              
!--------------------------------------------------------------------------------       
! Subroutine:                                                                   
! Type can be either real(4), real(8), complex(4) or complex(8),                
! but all arguments should have the same type.                                  
!                                                                               
! subroutine ezmatmul	(	                                                        
!  character*1, optional (default = 'N') :: transa,                             
!  character*1, optional (default = 'N') :: transb,                             
!  type,        optional (default =  1 ) :: alpha,                              
!  type, dimension(:,:), contiguous      :: A,                                  
!  type, dimension(:,:), contiguous      :: B,                                  
!  type,        optional (default =  0 ) :: beta,                               
!  type, dimension(:,:), contiguous      :: C,                                  
! )                                                                             
!                                                                               
! - transa: (optional, default = 'N') Specifies the operation applied to matrix A:
!     - 'N': No transpose, op(A) = A.
!     - 'T': Transpose, op(A) = A**T.
!     - 'C': Hermitian, op(A) = A**H.
! 
! - transb: (optional, default = 'N') Specifies the operation applied to matrix B:
!     - 'N': No transpose, op(B) = B.
!     - 'T': Transpose, op(B) = B**T.
!     - 'C': Hermitian, op(B) = B**H.
! 
! - alpha: (optional, default = 1) Scalar multiplier for op(A) * op(B).
! 
! - beta: (optional, default = 0) Scalar multiplier for C.
! 
! - A, B, C: Contiguous matrices with the same type.
!---------------------------------------------------------------------------------
! Examples:                                                                     
!                                                                               
! call ezmatmul(A,B,C) computes C:= A*B with the simplicity of matmul           
! while harnessing the high performance of the LAPACK library.                  
!                                                                               
! On the other hand, ezmatmul keep the versatility of gemm with for example:    
! call ezmatmul('C','T',(5d-4,8d6),A,B,(9d2,3d-7),C)                            
!---------------------------------------------------------------------------------

  implicit none

  private
  public :: ezmatmul

  ! Generic interface for ezmatmul.
  ! When ezmatmul is called, the program checks which of the following subroutines has the correct arguments.
  interface ezmatmul
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
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezsgemm('N', 'N', 1e0, A, B, 0e0, C)
  end subroutine ezsgemm_trans_alpha_beta

  subroutine ezsgemm_alpha_beta(transa, transb, A, B, C)
    ! Calls the SGEMM routines with default values for 
    ! alpha and beta (1e0 and 0e0, respectively).
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezsgemm(transA, transB, 1e0, A, B, 0e0, C)
  end subroutine ezsgemm_alpha_beta

  subroutine ezsgemm_trans_beta(alpha, A, B, C)
    ! Calls the SGEMM routines with default values for 
    ! transa, transb and beta ('N', 'N' and 0e0, respectively).
    real(4), intent(in)                                :: alpha
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezsgemm('N', 'N', alpha, A, B, 0e0, C)
  end subroutine ezsgemm_trans_beta

  subroutine ezsgemm_trans_alpha(A, B, beta, C)
    ! Calls the SGEMM routines with default values for 
    ! transa, transb and alpha  ('N', 'N' and 1e0, respectively).
    real(4), intent(in)                                :: beta
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezsgemm('N', 'N', 1e0, A, B, beta, C)
  end subroutine ezsgemm_trans_alpha

  subroutine ezsgemm_alpha(transa, transb, A, B, beta, C)
    ! Calls the SGEMM routines with default values for 
    ! alpha (1e0).
    real(4), intent(in)                                :: beta
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezsgemm(transA, transB, 1e0, A, B, beta, C)
  end subroutine ezsgemm_alpha

  subroutine ezsgemm_beta(transa, transb, alpha, A, B, C)
    ! Calls the SGEMM routines with default values for 
    ! beta (0e0).
    real(4), intent(in)                                :: alpha
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezsgemm(transA, transB, alpha, A, B, 0e0, C)
  end subroutine ezsgemm_beta

  subroutine ezsgemm_trans(alpha, A, B, beta, C)
    ! Calls the SGEMM routines with default values for 
    ! transa and transb ('N' and 'N', respectively).
    real(4), intent(in)                                :: alpha, beta
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezsgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ezsgemm_trans

!------------------------------------!
! Checking dimensions and call SGEMM !
!------------------------------------!

  subroutine ezsgemm(transa, transb, alpha, A, B, beta, C)
    ! Verify that the dimensions of the matrices are correct,
    ! and call sgemm with the appropriate arguments based on transa and transb.

    character(len=1), intent(in)                       :: transa, transb
    real(4), intent(in)                                :: alpha, beta
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C

    logical                                            :: test_passed
    integer                                            :: m, n, k

    ! Check dimensions
    test_passed = .true.
    m = size(C,1)
    n = size(C,2)

    if (transa == 'N') then
      k = size(A,2)
      if (transb == 'N') then
        if (k /= size(B, 1) .or. m /= size(A, 1) .or. n /= size(B, 2)) then
          test_passed = .false.
        else
          call sgemm(transa, transb, m, n, k, alpha, A, m, B, k, beta, C, m)
        end if
      else
        if (k /= size(B, 2) .or. m /= size(A, 1) .or. n /= size(B, 1)) then
          test_passed = .false.
        else
          call sgemm(transa, transb, m, n, k, alpha, A, m, B, n, beta, C, m)
        end if
      end if
    else
      k = size(A,1)
      if (transb == 'N') then
        if (k /= size(B, 1) .or. m /= size(A, 2) .or. n /= size(B, 2)) then
          test_passed = .false.
        else
          call sgemm(transa, transb, m, n, k, alpha, A, k, B, k, beta, C, m)
        end if
      else
        if (k /= size(B,2) .or. m /= size(A, 2) .or. n /= size(B, 1)) then
          test_passed = .false.
        else
          call sgemm(transa, transb, m, n, k, alpha, A, k, B, n, beta, C, m)
        end if
      end if
    end if

    if (.not. test_passed) then
      print *, "Dimensions don't match!"
      write(*, '(A, "A: (", I5, ", ", I5, ")")') "Size of the matrices: ", m, k
      write(*, '(A, "B: (", I5, ", ", I5, ")")') "                      ", k, n
      write(*, '(A, "C: (", I5, ", ", I5, ")")') "                      ", m, n
    end if

  end subroutine ezsgemm


!------------!
! DGEMM PART !
!------------!

!----------------------------!
! Handling default arguments !
!----------------------------!

  subroutine ezdgemm_trans_alpha_beta(A, B, C)
    ! Calls the DGEMM routines with default values for 
    ! transa, transb, alpha and beta ('N', 'N', 1d0 and 0d0, respectively).
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezdgemm('N', 'N', 1d0, A, B, 0d0, C)
  end subroutine ezdgemm_trans_alpha_beta

  subroutine ezdgemm_alpha_beta(transa, transb, A, B, C)
    ! Calls the DGEMM routines with default values for 
    ! alpha and beta (1d0 and 0d0, respectively).
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezdgemm(transA, transB, 1d0, A, B, 0d0, C)
  end subroutine ezdgemm_alpha_beta

  subroutine ezdgemm_trans_beta(alpha, A, B, C)
    ! Calls the DGEMM routines with default values for 
    ! transa, transb and beta ('N', 'N' and 0d0, respectively).
    real(8), intent(in)                                :: alpha
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezdgemm('N', 'N', alpha, A, B, 0d0, C)
  end subroutine ezdgemm_trans_beta

  subroutine ezdgemm_trans_alpha(A, B, beta, C)
    ! Calls the DGEMM routines with default values for 
    ! transa, transb and alpha  ('N', 'N' and 1d0, respectively).
    real(8), intent(in)                                :: beta
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezdgemm('N', 'N', 1d0, A, B, beta, C)
  end subroutine ezdgemm_trans_alpha

  subroutine ezdgemm_alpha(transa, transb, A, B, beta, C)
    ! Calls the DGEMM routines with default values for 
    ! alpha (1d0).
    real(8), intent(in)                                :: beta
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezdgemm(transA, transB, 1d0, A, B, beta, C)
  end subroutine ezdgemm_alpha

  subroutine ezdgemm_beta(transa, transb, alpha, A, B, C)
    ! Calls the DGEMM routines with default values for 
    ! beta (0d0).
    real(8), intent(in)                                :: alpha
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezdgemm(transA, transB, alpha, A, B, 0d0, C)
  end subroutine ezdgemm_beta

  subroutine ezdgemm_trans(alpha, A, B, beta, C)
    ! Calls the DGEMM routines with default values for 
    ! transa and transb ('N' and 'N', respectively).
    real(8), intent(in)                                :: alpha, beta
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezdgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ezdgemm_trans

!------------------------------------!
! Checking dimensions and call DGEMM !
!------------------------------------!

  subroutine ezdgemm(transa, transb, alpha, A, B, beta, C)
    ! Verify that the dimensions of the matrices are correct,
    ! and call dgemm with the appropriate arguments based on transa and transb.

    character(len=1), intent(in)                       :: transa, transb
    real(8), intent(in)                                :: alpha, beta
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C

    logical                                            :: test_passed
    integer                                            :: m, n, k

    ! Check dimensions
    test_passed = .true.
    m = size(C,1)
    n = size(C,2)

    if (transa == 'N') then
      k = size(A,2)
      if (transb == 'N') then
        if (k /= size(B, 1) .or. m /= size(A, 1) .or. n /= size(B, 2)) then
          test_passed = .false.
        else
          call dgemm(transa, transb, m, n, k, alpha, A, m, B, k, beta, C, m)
        end if
      else
        if (k /= size(B, 2) .or. m /= size(A, 1) .or. n /= size(B, 1)) then
          test_passed = .false.
        else
          call dgemm(transa, transb, m, n, k, alpha, A, m, B, n, beta, C, m)
        end if
      end if
    else
      k = size(A,1)
      if (transb == 'N') then
        if (k /= size(B, 1) .or. m /= size(A, 2) .or. n /= size(B, 2)) then
          test_passed = .false.
        else
          call dgemm(transa, transb, m, n, k, alpha, A, k, B, k, beta, C, m)
        end if
      else
        if (k /= size(B,2) .or. m /= size(A, 2) .or. n /= size(B, 1)) then
          test_passed = .false.
        else
          call dgemm(transa, transb, m, n, k, alpha, A, k, B, n, beta, C, m)
        end if
      end if
    end if

    if (.not. test_passed) then
      print *, "Dimensions don't match!"
      write(*, '(A, "A: (", I5, ", ", I5, ")")') "Size of the matrices: ", m, k
      write(*, '(A, "B: (", I5, ", ", I5, ")")') "                      ", k, n
      write(*, '(A, "C: (", I5, ", ", I5, ")")') "                      ", m, n
    end if

  end subroutine ezdgemm

!------------!
! CGEMM PART !
!------------!

!----------------------------!
! Handling default arguments !
!----------------------------!

  subroutine ezcgemm_trans_alpha_beta(A, B, C)
    ! Calls the CGEMM routines with default values for 
    ! transa, transb, alpha and beta ('N', 'N', (1e0,0e0) and (0e0,0e0), respectively).
    complex(4), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezcgemm('N', 'N', (1e0,0e0), A, B, (0e0,0e0), C)
  end subroutine ezcgemm_trans_alpha_beta

  subroutine ezcgemm_alpha_beta(transa, transb, A, B, C)
    ! Calls the CGEMM routines with default values for 
    ! alpha and beta ((1e0,0e0) and (0e0,0e0), respectively).
    complex(4), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                          :: transa, transb
    call ezcgemm(transA, transB, (1e0,0e0), A, B, (0e0,0e0), C)
  end subroutine ezcgemm_alpha_beta

  subroutine ezcgemm_trans_beta(alpha, A, B, C)
    ! Calls the CGEMM routines with default values for 
    ! transa, transb and beta ('N', 'N' and (0e0,0e0), respectively).
    complex(4), intent(in)                                :: alpha
    complex(4), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezcgemm('N', 'N', alpha, A, B, (0e0,0e0), C)
  end subroutine ezcgemm_trans_beta

  subroutine ezcgemm_trans_alpha(A, B, beta, C)
    ! Calls the CGEMM routines with default values for 
    ! transa, transb and alpha ('N', 'N' and (1e0,0e0), respectively).
    complex(4), intent(in)                                :: beta
    complex(4), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezcgemm('N', 'N', (1e0,0e0), A, B, beta, C)
  end subroutine ezcgemm_trans_alpha

  subroutine ezcgemm_alpha(transa, transb, A, B, beta, C)
    ! Calls the CGEMM routines with default values for 
    ! alpha ((1e0,0e0)).
    complex(4), intent(in)                                :: beta
    complex(4), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                          :: transa, transb
    call ezcgemm(transA, transB, (1e0,0e0), A, B, beta, C)
  end subroutine ezcgemm_alpha

  subroutine ezcgemm_beta(transa, transb, alpha, A, B, C)
    ! Calls the CGEMM routines with default values for 
    ! beta ((0e0,0e0)).
    complex(4), intent(in)                                :: alpha
    complex(4), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                          :: transa, transb
    call ezcgemm(transA, transB, alpha, A, B, (0e0,0e0), C)
  end subroutine ezcgemm_beta

  subroutine ezcgemm_trans(alpha, A, B, beta, C)
    ! Calls the CGEMM routines with default values for 
    ! transa and transb ('N' and 'N', respectively).
    complex(4), intent(in)                                :: alpha, beta
    complex(4), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezcgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ezcgemm_trans

!------------------------------------!
! Checking dimensions and call CGEMM !
!------------------------------------!

  subroutine ezcgemm(transa, transb, alpha, A, B, beta, C)
    ! Verify that the dimensions of the matrices are correct,
    ! and call cgemm with the appropriate arguments based on transa and transb.

    character(len=1), intent(in)                          :: transa, transb
    complex(4), intent(in)                                :: alpha, beta
    complex(4), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(4), dimension(:,:), contiguous, intent(inout) :: C

    logical                                               :: test_passed
    integer                                               :: m, n, k

    ! Check dimensions
    test_passed = .true.
    m = size(C,1)
    n = size(C,2)

    if (transa == 'N') then
      k = size(A,2)
      if (transb == 'N') then
        if (k /= size(B, 1) .or. m /= size(A, 1) .or. n /= size(B, 2)) then
          test_passed = .false.
        else
          call cgemm(transa, transb, m, n, k, alpha, A, m, B, k, beta, C, m)
        end if
      else
        if (k /= size(B, 2) .or. m /= size(A, 1) .or. n /= size(B, 1)) then
          test_passed = .false.
        else
          call cgemm(transa, transb, m, n, k, alpha, A, m, B, n, beta, C, m)
        end if
      end if
    else
      k = size(A,1)
      if (transb == 'N') then
        if (k /= size(B, 1) .or. m /= size(A, 2) .or. n /= size(B, 2)) then
          test_passed = .false.
        else
          call cgemm(transa, transb, m, n, k, alpha, A, k, B, k, beta, C, m)
        end if
      else
        if (k /= size(B,2) .or. m /= size(A, 2) .or. n /= size(B, 1)) then
          test_passed = .false.
        else
          call cgemm(transa, transb, m, n, k, alpha, A, k, B, n, beta, C, m)
        end if
      end if
    end if

  end subroutine ezcgemm

!------------!
! ZGEMM PART !
!------------!

!----------------------------!
! Handling default arguments !
!----------------------------!

  subroutine ezzgemm_trans_alpha_beta(A, B, C)
    ! Calls the ZGEMM routines with default values for 
    ! transa, transb, alpha and beta ('N', 'N', (1d0,0d0) and (0d0,0d0), respectively).
    complex(8), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezzgemm('N', 'N', (1d0,0d0), A, B, (0d0,0d0), C)
  end subroutine ezzgemm_trans_alpha_beta

  subroutine ezzgemm_alpha_beta(transa, transb, A, B, C)
    ! Calls the ZGEMM routines with default values for 
    ! alpha and beta ((1d0,0d0) and (0d0,0d0), respectively).
    complex(8), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                          :: transa, transb
    call ezzgemm(transA, transB, (1d0,0d0), A, B, (0d0,0d0), C)
  end subroutine ezzgemm_alpha_beta

  subroutine ezzgemm_trans_beta(alpha, A, B, C)
    ! Calls the ZGEMM routines with default values for 
    ! transa, transb and beta ('N', 'N' and (0d0,0d0), respectively).
    complex(8), intent(in)                                :: alpha
    complex(8), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezzgemm('N', 'N', alpha, A, B, (0d0,0d0), C)
  end subroutine ezzgemm_trans_beta

  subroutine ezzgemm_trans_alpha(A, B, beta, C)
    ! Calls the ZGEMM routines with default values for 
    ! transa, transb and alpha ('N', 'N' and (1d0,0d0), respectively).
    complex(8), intent(in)                                :: beta
    complex(8), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezzgemm('N', 'N', (1d0,0d0), A, B, beta, C)
  end subroutine ezzgemm_trans_alpha

  subroutine ezzgemm_alpha(transa, transb, A, B, beta, C)
    ! Calls the ZGEMM routines with default values for 
    ! alpha ((1d0,0d0)).
    complex(8), intent(in)                                :: beta
    complex(8), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                      :: transa, transb
    call ezzgemm(transA, transB, (1d0,0d0), A, B, beta, C)
  end subroutine ezzgemm_alpha

  subroutine ezzgemm_beta(transa, transb, alpha, A, B, C)
    ! Calls the ZGEMM routines with default values for 
    ! beta ((0d0,0d0)).
    complex(8), intent(in)                                :: alpha
    complex(8), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                          :: transa, transb
    call ezzgemm(transA, transB, alpha, A, B, (0d0,0d0), C)
  end subroutine ezzgemm_beta

  subroutine ezzgemm_trans(alpha, A, B, beta, C)
    ! Calls the ZGEMM routines with default values for 
    ! transa and transb ('N' and 'N', respectively).
    complex(8), intent(in)                                :: alpha, beta
    complex(8), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezzgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ezzgemm_trans

!------------------------------------!
! Checking dimensions and call ZGEMM !
!------------------------------------!

  subroutine ezzgemm(transa, transb, alpha, A, B, beta, C)
    ! Verify that the dimensions of the matrices are correct,
    ! and call zgemm with the appropriate arguments based on transa and transb.

    character(len=1), intent(in)                          :: transa, transb
    complex(8), intent(in)                                :: alpha, beta
    complex(8), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(8), dimension(:,:), contiguous, intent(inout) :: C

    logical                                               :: test_passed
    integer                                               :: m, n, k

    ! Check dimensions
    test_passed = .true.
    m = size(C,1)
    n = size(C,2)

    if (transa == 'N') then
      k = size(A,2)
      if (transb == 'N') then
        if (k /= size(B, 1) .or. m /= size(A, 1) .or. n /= size(B, 2)) then
          test_passed = .false.
        else
          call zgemm(transa, transb, m, n, k, alpha, A, m, B, k, beta, C, m)
        end if
      else
        if (k /= size(B, 2) .or. m /= size(A, 1) .or. n /= size(B, 1)) then
          test_passed = .false.
        else
          call zgemm(transa, transb, m, n, k, alpha, A, m, B, n, beta, C, m)
        end if
      end if
    else
      k = size(A,1)
      if (transb == 'N') then
        if (k /= size(B, 1) .or. m /= size(A, 2) .or. n /= size(B, 2)) then
          test_passed = .false.
        else
          call zgemm(transa, transb, m, n, k, alpha, A, k, B, k, beta, C, m)
        end if
      else
        if (k /= size(B,2) .or. m /= size(A, 2) .or. n /= size(B, 1)) then
          test_passed = .false.
        else
          call zgemm(transa, transb, m, n, k, alpha, A, k, B, n, beta, C, m)
        end if
      end if
    end if

  end subroutine ezzgemm
  
end module module_ezmatmul
