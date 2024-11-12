subroutine matrices_mult_d(alpha, transa, A, transb, B, beta, C)
  
  implicit none

  ! Inputs
  character(len=1), intent(in)    :: transa ! 'N' or 'n' for A, else A**T
  character(len=1), intent(in)    :: transb ! 'N' or 'n' for B, else B**T

  double precision, intent(in)    :: alpha
  double precision, intent(in)    :: beta

  double precision, intent(in)    :: A(:,:)
  double precision, intent(in)    :: B(:,:)
  double precision, intent(inout) :: C(:,:)

  ! Local variables
  integer :: m, n, k
  integer :: lda, ldb, ldc

  ! Determine dimensions
  m = size(C, 1)
  n = size(C, 2)
  
  if (transa == 'N' .or. transa == 'n') then
    k = size(A, 2)
    lda = m
  else
    k = size(A, 1)
    lda = k
  end if

  if (transb == 'N' .or. transb == 'n') then
    ldb = k
  else
    ldb = n
  end if

  ldc = m

  ! Call the DGEMM routine
  call dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)

end subroutine matrices_mult_d

