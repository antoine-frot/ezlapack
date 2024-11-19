program test_matrice_mult_d

  implicit none

  ! external         :: matrices_mult_d

  integer          :: N
  integer          :: i,p,q,r

  double precision :: t_start
  double precision :: t_end
  integer          :: size_mat(6)

  double precision :: alpha
  double precision :: beta

  double precision,allocatable :: A(:,:)
  double precision,allocatable :: B(:,:)
  double precision,allocatable :: C(:,:)

  size_mat = (/ 10, 100, 1000, 2000, 5000, 10000 /)

  do i=1,size(size_mat)
 
    N = size_mat(i)

    allocate(A(N,N),B(N,N),C(N,N))
 
    call random_number(A)
    call random_number(B)

    alpha = 1d0
    beta = 0d0

    ! matmul section

    call cpu_time(t_start)
 
    C = matmul(A,B)
 
    call cpu_time(t_end)

    write(*,'(A30,E12.6,A10,I9)') 'CPU time for matmul ',t_end - t_start,' s for N = ',N
 
    ! dgemm section

    call cpu_time(t_start)
 
    call dgemm('N','N',N,N,N,alpha,A,N,B,N,beta,C,N)
 
    call cpu_time(t_end)

    write(*,'(A30,E12.6,A10,I9)') 'CPU time for dgemm ',t_end - t_start,' s for N = ',N
 
    ! wrapper section
 
    call cpu_time(t_start)
 
    call dgemm_wrapper(alpha, beta, A, B, C)

    call cpu_time(t_end)

    write(*,'(A30,E12.6,A10,I9)') 'CPU time for wrapper ',t_end - t_start,' s for N = ',N
 
    deallocate(A,B,C)

  end do
 
  contains

  subroutine dgemm_wrapper(alpha, beta, A, B, C, transa, transb)
    implicit none
    ! Arguments
    real(8), intent(in) :: alpha, beta
    real(8), intent(in) :: A(:,:), B(:,:)    ! Input matrices
    real(8), intent(inout) :: C(:,:)         ! Result matrix
    character(len=1), intent(in), optional :: transa, transb

    ! Local variables for DGEMM
    integer :: m, n, k
    real(8) :: alpha_real, beta_real
    character(len=1) :: transa_val, transb_val

    ! Check dimensions
    m = size(A, 1)      ! Rows of A and C
    k = size(A, 2)      ! Columns of A, rows of B
    n = size(B, 2)      ! Columns of B and C
    if (size(B, 1) /= k .or. size(C, 1) /= m .or. size(C, 2) /= n) then
        print *, "Error: Matrix dimensions do not match!"
        stop
    end if
    
    ! Assign default values if arguments are not provided
    if (present(transa)) then
        transa_val = transa
    else
        transa_val = 'N'
    end if

    if (present(transb)) then
        transb_val = transb
    else
        transb_val = 'N' 
    end if

    ! Call DGEMM
    call dgemm(transa_val, transb_val, m, n, k, alpha, A, m, B, k, beta, C, m)

end subroutine dgemm_wrapper

end program
