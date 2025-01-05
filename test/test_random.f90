program test_random
  ! Test file for random
  !
  ! Tests every kind of type for all dimension supported

  use ezlapack
  implicit none

  ! Variables for testing random complex numbers, integers, and characters
  complex(4), dimension(2)       :: z_scomplex_vector
  complex(4), dimension(3,2)     :: z_scomplex_matrix
  complex(8), dimension(4)       :: z_dcomplex_vector
  complex(8), dimension(2,3)     :: z_dcomplex_matrix
  integer                        :: int_output, i, j
  character(len=5)               :: char_output
  character(len=5), dimension(3) :: char_array = ['Alice', 'Betty', 'Catty']

  ! Test random single-precision complex number
  print *, '--------------------------------------'
  call random_complex(z_scomplex_vector(1))
  write(*, "(A, 2F12.6)") 'Random single-precision complex: ', z_scomplex_vector(1)
  print *, '--------------------------------------'

  ! Test random double-precision complex number
  call random_complex(z_dcomplex_vector(1))
  write(*, "(A, 2F18.10)") 'Random double-precision complex: ', z_dcomplex_vector(1)
  print *, '--------------------------------------'

  ! Test random single-precision complex rank-1 array
  call random_complex(z_scomplex_vector)
  print *, 'Random single-precision complex rank-1 array:'
  do i = 1, 2
    write(*, "(A, I1, A, A, F12.6, A, F12.6)") &
    '(', i, ') ', 'Real: ', real(z_scomplex_vector(i)), ' Imaginary: ', aimag(z_scomplex_vector(i))
  end do
  print *, '--------------------------------------'

  ! Test random double-precision complex rank-1 array
  call random_complex(z_dcomplex_vector)
  print *, 'Random double-precision complex rank-1 array:'
  do i = 1, 4
    write(*, "(A, I1, A, A, F18.15, A, F18.15)") &
    '(', i, ') ', 'Real: ', real(z_dcomplex_vector(i)), ' Imaginary: ', aimag(z_dcomplex_vector(i))
  end do
  print *, '--------------------------------------'

  ! Test random single-precision complex rank-2 array
  call random_complex(z_scomplex_matrix)
  print *, 'Random single-precision complex rank-2 array:'
  do i = 1, size(z_scomplex_matrix, dim=1)
      do j = 1, size(z_scomplex_matrix, dim=2)
          write(*, "(A, I1, A, I1, A, A, F12.6, A, F12.6)") &
               '(', i, ', ', j, ') ', 'Real: ', real(z_scomplex_matrix(i,j)), ' Imaginary: ', aimag(z_scomplex_matrix(i,j))
      end do
  end do
  print *, '--------------------------------------'

  ! Test random double-precision complex rank-2 array
  call random_complex(z_dcomplex_matrix)
  print *, 'Random double-precision complex rank-2 array:'
  do i = 1, size(z_dcomplex_matrix, dim=1)
      do j = 1, size(z_dcomplex_matrix, dim=2)
          write(*, "(A, I1, A, I1, A, A, F18.15, A, F18.15)") &
               '(', i, ', ', j, ') ', 'Real: ', real(z_dcomplex_matrix(i,j)), ' Imaginary: ', aimag(z_dcomplex_matrix(i,j))
      end do
  end do
  print *, '--------------------------------------'

  ! Test random integer generation
  call random_integer(1, 100, int_output)
  write(*, "(A, I3)") 'Random integer between 1 and 100: ', int_output
  print *, '--------------------------------------'

  ! Test random character selection
  print *, 'Testing random character selection from array: ', char_array(1), ', ', char_array(2), ', ', char_array(3)
  call random_character(char_array, char_output)
  write(*, "(A, A)") 'Random character: ', char_output
  print *, '--------------------------------------'

end program test_random
