module mod_ez_matmul

  implicit none

  interface ez_matmul
    module procedure ez_sgemm, ez_dgemm, ez_cgemm, ez_zgemm
    module procedure ez_dgemm_trans_alpha_beta
    module procedure ez_dgemm_alpha_beta
    module procedure ez_dgemm_trans_beta
    module procedure ez_dgemm_trans_alpha
    module procedure ez_dgemm_alpha
    module procedure ez_dgemm_beta
    module procedure ez_dgemm_trans
  end interface

contains

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

  subroutine ez_sgemm(transa, transb, alpha, A, B, beta, C)

    implicit none

    real(4), intent(in), optional                      :: alpha, beta
    real(4), dimension(:,:), contiguous, intent(in)    :: A, B
    real(4), dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in), optional             :: transa, transb

    integer                                            :: m, n, k
    real(4)                                            :: alpha_v, beta_v
    character(len=1)                                   :: transa_v, transb_v

    ! Assign default values if arguments are not provided
    if (present(transa)) then
      transa_v = transa
    else
      transa_v = 'N'
    end if

    if (present(transb)) then
      transb_v = transb
    else
      transb_v = 'N'
    end if
    
    if (present(alpha)) then
      alpha_v = alpha
    else
      alpha_v = 1.
    end if

    if (present(beta)) then
      beta_v = beta
    else
      beta_v = 0.
    end if
     
    ! Check dimensions
    m = size(C,1)
    k = size(C,2)

    if (transa_v == 'N') then
      n = size(A,2)
      if (transb_v == 'N') then
        if (n /= size(B,1) .or. m /= size(A,1) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,1) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          stop
        end if
      end if
    else
      n = size(A,1)
      if (transb_v == 'N') then
        if (n /= size(B,1) .or. m /= size(A,2) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,2) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          stop
        end if
      end if
    end if

    call sgemm(transa_v, transb_v, m, n, k, alpha, A, m, B, k, beta, C, m)

  end subroutine ez_sgemm

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
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,1) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          stop
        end if
      end if
    else
      n = size(A,1)
      if (transb == 'N') then
        if (n /= size(B,1) .or. m /= size(A,2) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,2) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          stop
        end if
      end if
    end if

    call dgemm(transa, transb, m, n, k, alpha, A, m, B, k, beta, C, m)

  end subroutine ez_dgemm

  subroutine ez_cgemm(transa, transb, alpha, A, B, beta, C)

    implicit none

    complex, intent(in), optional                      :: alpha, beta
    complex, dimension(:,:), contiguous, intent(in)    :: A, B
    complex, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in), optional             :: transa, transb

    integer                                            :: m, n, k
    complex                                            :: alpha_v, beta_v
    character(len=1)                                   :: transa_v, transb_v

    ! Assign default values if arguments are not provided
    if (present(transa)) then
      transa_v = transa
    else
      transa_v = 'N'
    end if

    if (present(transb)) then
      transb_v = transb
    else
      transb_v = 'N'
    end if
    
    if (present(alpha)) then
      alpha_v = alpha
    else
      alpha_v = (1.0, 0.0)
    end if

    if (present(beta)) then
      beta_v = beta
    else
      beta_v = (0.0, 0.0)
    end if
     
    ! Check dimensions
    m = size(C,1)
    k = size(C,2)

    if (transa_v == 'N') then
      n = size(A,2)
      if (transb_v == 'N') then
        if (n /= size(B,1) .or. m /= size(A,1) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,1) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          stop
        end if
      end if
    else
      n = size(A,1)
      if (transb_v == 'N') then
        if (n /= size(B,1) .or. m /= size(A,2) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,2) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          stop
        end if
      end if
    end if

    call cgemm(transa_v, transb_v, m, n, k, alpha, A, m, B, k, beta, C, m)

  end subroutine ez_cgemm

  subroutine ez_zgemm(transa, transb, alpha, A, B, beta, C)

    implicit none

    complex*16, intent(in), optional                      :: alpha, beta
    complex*16, dimension(:,:), contiguous, intent(in)    :: A, B
    complex*16, dimension(:,:), contiguous, intent(inout) :: C
    character(len=1), intent(in), optional             :: transa, transb

    integer                                            :: m, n, k
    complex*16                                         :: alpha_v, beta_v
    character(len=1)                                   :: transa_v, transb_v

    ! Assign default values if arguments are not provided
    if (present(transa)) then
      transa_v = transa
    else
      transa_v = 'N'
    end if

    if (present(transb)) then
      transb_v = transb
    else
      transb_v = 'N'
    end if
    
    if (present(alpha)) then
      alpha_v = alpha
    else
      alpha_v = (1d0, 0d0)
    end if

    if (present(beta)) then
      beta_v = beta
    else
      beta_v = (0d0, 0d0)
    end if
     
    ! Check dimensions
    m = size(C,1)
    k = size(C,2)

    if (transa_v == 'N') then
      n = size(A,2)
      if (transb_v == 'N') then
        if (n /= size(B,1) .or. m /= size(A,1) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,1) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          stop
        end if
      end if
    else
      n = size(A,1)
      if (transb_v == 'N') then
        if (n /= size(B,1) .or. m /= size(A,2) .or. k /= size(B,2)) then
          print *, "Dimensions don't match!"
          stop
        end if
      else
        if (n /= size(B,2) .or. m /= size(A,2) .or. k /= size(B,1)) then
          print *, "Dimensions don't match!"
          stop
        end if
      end if
    end if

    call zgemm(transa_v, transb_v, m, n, k, alpha, A, m, B, k, beta, C, m)

  end subroutine ez_zgemm

end module mod_ez_matmul
