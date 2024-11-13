program test_matrice_mult_d

  implicit none

  external         :: matrices_mult_d

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
 
    print *, A
    call random_number(A)
    call random_number(B)
    print *, A

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
 
    call matrices_mult_d(alpha, 'N', A, 'N', B, beta, C)

    call cpu_time(t_end)

    write(*,'(A30,E12.6,A10,I9)') 'CPU time for wrapper ',t_end - t_start,' s for N = ',N
 
    deallocate(A,B,C)

  end do
 
end program
