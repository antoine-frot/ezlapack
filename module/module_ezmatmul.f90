!> @brief Wrapper for the BLAS ( part of LAPACK) `*gemm` routines that perform matrix-matrix multiplication.
!! These subroutines compute the operation:
!! `C = alpha * op(A) * op(B) + beta * C`
!! where `op(A)` and `op(B)` can be the transpose or the original matrices `A` and `B`.
!! It validates the dimensions of the input matrices and provides default values for
!! the scaling factors `alpha` and `beta`, as well as the transformation indicators
!! `transa` and `transb`, if they are not supplied.
!! @param[in] transa (Optional, default = 'N') Specifies the operation applied to matrix A:
!!                    - 'N': No transpose
!!                    - 'T': Transpose
!!                    - 'C': Conjugate transpose
!! @param[in] transb (Optional, default = 'N') Specifies the operation applied to matrix B:
!!                    - 'N': No transpose
!!                    - 'T': Transpose
!!                    - 'C': Conjugate transpose
!! @param[in] alpha  (Optional, default = 1) Scalar multiplier for `op(A) * op(B)`.
!! @param[in] beta   (Optional, default = 0) Scalar multiplier for `C`.
!! @param[in] A      Input matrix A with dimensions matching the specified operation.
!! @param[in] B      Input matrix B with dimensions matching the specified operation.
!! @param[in,out] C  Matrix C. On input, the initial values of C. On output, the result of
!!                   the matrix multiplication.
!! @note If matrix dimensions are incompatible, an error message is printed, and the
!!       computation is skipped.
module module_ezmatmul

  implicit none

  private
  public :: ezmatmul

  !> @brief Wrapper for the BLAS ( part of LAPACK) `*gemm` routines that perform matrix-matrix multiplication.
  !! These subroutines compute the operation:
  !! `C = alpha * op(A) * op(B) + beta * C`
  !! where `op(A)` and `op(B)` can be the transpose or the original matrices `A` and `B`.
  !! It validates the dimensions of the input matrices and provides default values for
  !! the scaling factors `alpha` and `beta`, as well as the transformation indicators
  !! `transa` and `transb`, if they are not supplied.
  !! @param[in] transa (Optional, default = 'N') Specifies the operation applied to matrix A:
  !!                    - 'N': No transpose
  !!                    - 'T': Transpose
  !!                    - 'C': Conjugate transpose
  !! @param[in] transb (Optional, default = 'N') Specifies the operation applied to matrix B:
  !!                    - 'N': No transpose
  !!                    - 'T': Transpose
  !!                    - 'C': Conjugate transpose
  !! @param[in] alpha  (Optional, default = 1) Scalar multiplier for `op(A) * op(B)`.
  !! @param[in] beta   (Optional, default = 0) Scalar multiplier for `C`.
  !! @param[in] A      Input matrix A with dimensions matching the specified operation.
  !! @param[in] B      Input matrix B with dimensions matching the specified operation.
  !! @param[in,out] C  Matrix C. On input, the initial values of C. On output, the result of
  !!                   the matrix multiplication.
  !! @note If matrix dimensions are incompatible, an error message is printed, and the
  !!       computation is skipped.
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

  !> \brief Calls the ezsgemm routines with default values for
  !> transa, transb, alpha, and beta ('N', 'N', 1e0 and 0e0, respectively).
  subroutine ezsgemm_trans_alpha_beta(A, B, C)
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezsgemm('N', 'N', 1e0, A, B, 0e0, C)
  end subroutine ezsgemm_trans_alpha_beta

  !> \brief Calls the ezsgemm routines with default values for
  !> alpha and beta (1e0 and 0e0, respectively).
  subroutine ezsgemm_alpha_beta(transa, transb, A, B, C)
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezsgemm(transA, transB, 1e0, A, B, 0e0, C)
  end subroutine ezsgemm_alpha_beta

  !> \brief Calls the ezsgemm routines with default values for
  !> transa, transb, and beta ('N', 'N' and 0e0, respectively).
  subroutine ezsgemm_trans_beta(alpha, A, B, C)
    real(4), intent(in)                                :: alpha
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezsgemm('N', 'N', alpha, A, B, 0e0, C)
  end subroutine ezsgemm_trans_beta

  !> \brief Calls the ezsgemm routines with default values for
  !> transa, transb, and alpha ('N', 'N' and 1e0, respectively).
  subroutine ezsgemm_trans_alpha(A, B, beta, C)
    real(4), intent(in)                                :: beta
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezsgemm('N', 'N', 1e0, A, B, beta, C)
  end subroutine ezsgemm_trans_alpha

  !> \brief Calls the ezsgemm routines with default values for
  !> alpha (1e0).
  subroutine ezsgemm_alpha(transa, transb, A, B, beta, C)
    real(4), intent(in)                                :: beta
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezsgemm(transA, transB, 1e0, A, B, beta, C)
  end subroutine ezsgemm_alpha

  !> \brief Calls the ezsgemm routines with default values for
  !> beta (0e0).
  subroutine ezsgemm_beta(transa, transb, alpha, A, B, C)
    real(4), intent(in)                                :: alpha
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezsgemm(transA, transB, alpha, A, B, 0e0, C)
  end subroutine ezsgemm_beta

  !> \brief Calls the ezsgemm routines with default values for
  !> transa and transb ('N' and 'N', respectively).
  subroutine ezsgemm_trans(alpha, A, B, beta, C)
    real(4), intent(in)                                :: alpha, beta
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezsgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ezsgemm_trans

!------------------------------------!
! Checking dimensions and call SGEMM !
!------------------------------------!

  !> \brief Verifies that the dimensions of the matrices are correct,
  !> and calls sgemm with the appropriate arguments based on transa and transb.
  subroutine ezsgemm(transa, transb, alpha, A, B, beta, C)

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

  !> \brief Calls the ezdgemm routines with default values for
  !> transa, transb, alpha, and beta ('N', 'N', 1d0 and 0d0, respectively).
  subroutine ezdgemm_trans_alpha_beta(A, B, C)
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezdgemm('N', 'N', 1d0, A, B, 0d0, C)
  end subroutine ezdgemm_trans_alpha_beta

  !> \brief Calls the ezdgemm routines with default values for
  !> alpha and beta (1d0 and 0d0, respectively).
  subroutine ezdgemm_alpha_beta(transa, transb, A, B, C)
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezdgemm(transA, transB, 1d0, A, B, 0d0, C)
  end subroutine ezdgemm_alpha_beta

  !> \brief Calls the ezdgemm routines with default values for
  !> transa, transb, and beta ('N', 'N' and 0d0, respectively).
  subroutine ezdgemm_trans_beta(alpha, A, B, C)
    real(8), intent(in)                                :: alpha
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezdgemm('N', 'N', alpha, A, B, 0d0, C)
  end subroutine ezdgemm_trans_beta

  !> \brief Calls the ezdgemm routines with default values for
  !> transa, transb, and alpha ('N', 'N' and 1d0, respectively).
  subroutine ezdgemm_trans_alpha(A, B, beta, C)
    real(8), intent(in)                                :: beta
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezdgemm('N', 'N', 1d0, A, B, beta, C)
  end subroutine ezdgemm_trans_alpha

  !> \brief Calls the ezdgemm routines with default values for
  !> alpha (1d0).
  subroutine ezdgemm_alpha(transa, transb, A, B, beta, C)
    real(8), intent(in)                                :: beta
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezdgemm(transA, transB, 1d0, A, B, beta, C)
  end subroutine ezdgemm_alpha

  !> \brief Calls the ezdgemm routines with default values for
  !> beta (0d0).
  subroutine ezdgemm_beta(transa, transb, alpha, A, B, C)
    real(8), intent(in)                                :: alpha
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                       :: transa, transb
    call ezdgemm(transA, transB, alpha, A, B, 0d0, C)
  end subroutine ezdgemm_beta

  !> \brief Calls the ezdgemm routines with default values for
  !> transa and transb ('N' and 'N', respectively).
  subroutine ezdgemm_trans(alpha, A, B, beta, C)
    real(8), intent(in)                                :: alpha, beta
    real(8), dimension(:,:), contiguous, intent(in)    :: A, B
    real(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezdgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ezdgemm_trans

!------------------------------------!
! Checking dimensions and call DGEMM !
!------------------------------------!

  !> \brief Verifies that the dimensions of the matrices are correct,
  !> and calls dgemm with the appropriate arguments based on transa and transb.
  subroutine ezdgemm(transa, transb, alpha, A, B, beta, C)

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

  !> \brief Calls the ezcgemm routines with default values for
  !> transa, transb, alpha, and beta ('N', 'N', (1e0,0e0) and (0e0,0e0), respectively).
  subroutine ezcgemm_trans_alpha_beta(A, B, C)
    complex(4), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezcgemm('N', 'N', (1e0,0e0), A, B, (0e0,0e0), C)
  end subroutine ezcgemm_trans_alpha_beta

  !> \brief Calls the ezcgemm routines with default values for
  !> alpha and beta ((1e0,0e0) and (0e0,0e0), respectively).
  subroutine ezcgemm_alpha_beta(transa, transb, A, B, C)
    complex(4), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                          :: transa, transb
    call ezcgemm(transA, transB, (1e0,0e0), A, B, (0e0,0e0), C)
  end subroutine ezcgemm_alpha_beta

  !> \brief Calls the ezcgemm routines with default values for
  !> transa, transb, and beta ('N', 'N' and (0e0,0e0), respectively).
  subroutine ezcgemm_trans_beta(alpha, A, B, C)
    complex(4), intent(in)                                :: alpha
    complex(4), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezcgemm('N', 'N', alpha, A, B, (0e0,0e0), C)
  end subroutine ezcgemm_trans_beta

  !> \brief Calls the ezcgemm routines with default values for
  !> transa, transb, and alpha ('N', 'N' and (1e0,0e0), respectively).
  subroutine ezcgemm_trans_alpha(A, B, beta, C)
    complex(4), intent(in)                                :: beta
    complex(4), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezcgemm('N', 'N', (1e0,0e0), A, B, beta, C)
  end subroutine ezcgemm_trans_alpha

  !> \brief Calls the ezcgemm routines with default values for
  !> alpha ((1e0,0e0)).
  subroutine ezcgemm_alpha(transa, transb, A, B, beta, C)
    complex(4), intent(in)                                :: beta
    complex(4), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                          :: transa, transb
    call ezcgemm(transA, transB, (1e0,0e0), A, B, beta, C)
  end subroutine ezcgemm_alpha

  !> \brief Calls the ezcgemm routines with default values for
  !> beta ((0e0,0e0)).
  subroutine ezcgemm_beta(transa, transb, alpha, A, B, C)
    complex(4), intent(in)                                :: alpha
    complex(4), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                          :: transa, transb
    call ezcgemm(transA, transB, alpha, A, B, (0e0,0e0), C)
  end subroutine ezcgemm_beta

  !> \brief Calls the ezcgemm routines with default values for
  !> transa and transb ('N' and 'N', respectively).
  subroutine ezcgemm_trans(alpha, A, B, beta, C)
    complex(4), intent(in)                                :: alpha, beta
    complex(4), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(4), dimension(:,:), contiguous, intent(inout) :: C
    call ezcgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ezcgemm_trans

!------------------------------------!
! Checking dimensions and call CGEMM !
!------------------------------------!

  !> \brief Verifies that the dimensions of the matrices are correct,
  !> and calls cgemm with the appropriate arguments based on transa and transb.
  subroutine ezcgemm(transa, transb, alpha, A, B, beta, C)

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

  !> \brief Calls the ezzgemm routines with default values for
  !> transa, transb, alpha, and beta ('N', 'N', (1d0,0d0) and (0d0,0d0), respectively).
  subroutine ezzgemm_trans_alpha_beta(A, B, C)
    complex(8), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezzgemm('N', 'N', (1d0,0d0), A, B, (0d0,0d0), C)
  end subroutine ezzgemm_trans_alpha_beta

  !> \brief Calls the ezzgemm routines with default values for
  !> alpha and beta ((1d0,0d0) and (0d0,0d0), respectively).
  subroutine ezzgemm_alpha_beta(transa, transb, A, B, C)
    complex(8), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                          :: transa, transb
    call ezzgemm(transA, transB, (1d0,0d0), A, B, (0d0,0d0), C)
  end subroutine ezzgemm_alpha_beta

  !> \brief Calls the ezzgemm routines with default values for
  !> transa, transb, and beta ('N', 'N' and (0d0,0d0), respectively).
  subroutine ezzgemm_trans_beta(alpha, A, B, C)
    complex(8), intent(in)                                :: alpha
    complex(8), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezzgemm('N', 'N', alpha, A, B, (0d0,0d0), C)
  end subroutine ezzgemm_trans_beta

  !> \brief Calls the ezzgemm routines with default values for
  !> transa, transb, and alpha ('N', 'N' and (1d0,0d0), respectively).
  subroutine ezzgemm_trans_alpha(A, B, beta, C)
    complex(8), intent(in)                                :: beta
    complex(8), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezzgemm('N', 'N', (1d0,0d0), A, B, beta, C)
  end subroutine ezzgemm_trans_alpha

  !> \brief Calls the ezzgemm routines with default values for
  !> alpha ((1d0,0d0)).
  subroutine ezzgemm_alpha(transa, transb, A, B, beta, C)
    complex(8), intent(in)                                :: beta
    complex(8), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                          :: transa, transb
    call ezzgemm(transA, transB, (1d0,0d0), A, B, beta, C)
  end subroutine ezzgemm_alpha

  !> \brief Calls the ezzgemm routines with default values for
  !> beta ((0d0,0d0)).
  subroutine ezzgemm_beta(transa, transb, alpha, A, B, C)
    complex(8), intent(in)                                :: alpha
    complex(8), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(8), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in)                          :: transa, transb
    call ezzgemm(transA, transB, alpha, A, B, (0d0,0d0), C)
  end subroutine ezzgemm_beta

  !> \brief Calls the ezzgemm routines with default values for
  !> transa and transb ('N' and 'N', respectively).
  subroutine ezzgemm_trans(alpha, A, B, beta, C)
    complex(8), intent(in)                                :: alpha, beta
    complex(8), dimension(:,:), contiguous, intent(in)    :: A, B
    complex(8), dimension(:,:), contiguous, intent(inout) :: C
    call ezzgemm('N', 'N', alpha, A, B, beta, C)
  end subroutine ezzgemm_trans

!------------------------------------!
! Checking dimensions and call ZGEMM !
!------------------------------------!

  !> \brief Verifies that the dimensions of the matrices are correct,
  !> and calls zgemm with the appropriate arguments based on transa and transb.
  subroutine ezzgemm(transa, transb, alpha, A, B, beta, C)
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
