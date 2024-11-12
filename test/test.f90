program test_dgemm_wrapper
    use dgemm_module
    
    implicit none

    double precision, allocatable :: A(:,:), B(:,:), C(:,:)
    double precision :: alpha, beta
    integer :: m, n, k

    alpha = 1.0d0
    beta = 0.0d0

    m = 3      ! Number of rows of C and A (if transa = 'N')
    n = 3      ! Number of columns of C and B (if transb = 'N')
    k = 3      ! Number of columns of A (if transa = 'N') or rows of B (if transb = 'N')

    allocate(A(m, k), B(k, n), C(m, n))

    A = reshape([ &
        1.0d0, 2.0d0, 3.0d0, &
        4.0d0, 5.0d0, 6.0d0, &
        7.0d0, 8.0d0, 9.0d0], &
        shape(A))

    B = reshape([ &
        9.0d0, 8.0d0, 7.0d0, &
        6.0d0, 5.0d0, 4.0d0, &
        3.0d0, 2.0d0, 1.0d0], &
        shape(B))

    C = 0.0d0

    call matrices_mult_d(alpha, 'N', A, 'N', B, beta, C)

    print *, 'Resulting Matrix C:'
    call print_matrix(C)

    deallocate(A, B, C)

contains

    subroutine print_matrix(mat)
        double precision, intent(in) :: mat(:,:)
        integer :: i, j

        do i = 1, size(mat, 1)
            do j = 1, size(mat, 2)
                write(*, '(F6.2)', advance='no') mat(i, j)
                if (j < size(mat, 2)) then
                    write(*, '(A)', advance='no') '  '
                end if
            end do
            print *
        end do
    end subroutine print_matrix

end program test_dgemm_wrapper

