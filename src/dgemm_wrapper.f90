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
    lda = size(A,1)
  else
    k = size(A, 1)
    lda = size(A,2)
  end if

  if (transb == 'N' .or. transb == 'n') then
    ldb = size(B,1)
  else
    ldb = size(B,1)
  end if

  ldc = m

  ! Debugging output to verify dimensions
  print *, 'Dimensions:'
  print *, 'm = ', m, ', n = ', n, ', k = ', k
  print *, 'lda = ', lda, ', ldb = ', ldb, ', ldc = ', ldc

  ! Ensure sizes are compatible before calling dgemm
  if (size(A,1) /= lda .or. size(B,1) /= ldb .or. size(C,1) /= ldc) then
    print *, 'Error: Mismatched dimensions'
    stop
  end if

  ! Call the DGEMM routine
  call dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)

end subroutine matrices_mult_d

