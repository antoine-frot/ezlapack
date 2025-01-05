module module_random
  !-------------------------------------------------------------------------
  ! Purpose: Provides a generic interface for generating random complex 
  !          numbers, integers, and characters.
  !
  ! Description:
  ! This module offers utilities for generating random values with specific
  ! properties:
  ! - Random complex numbers within the unit circle.
  ! - Random integers within a user-specified range.
  ! - Random selection from an array of characters.
  !
  ! Public Subroutines:
  ! - random_complex(output_complex)
  !   Generate random complex numbers (single or double precision) within
  !   the unit circle. Supports scalars and arrays up to rank-2.
  !
  ! - random_integer(low, high, output_integer)
  !   Generate a random integer between 'low' and 'high', inclusive.
  !
  ! - random_character(character_array, output_character)
  !   Randomly select and return a character from an input array of strings.
  !-------------------------------------------------------------------------

  implicit none

  private
  public :: random_complex, random_integer, random_character
  
  ! Generic interface for random_complex
  ! When random_complex is called, the program checks which of the following subroutines has the correct arguments.
  interface random_complex
    module procedure random_scomplex, random_dcomplex
    module procedure random_scomplex_rank1, random_dcomplex_rank1
    module procedure random_scomplex_rank2, random_dcomplex_rank2
  end interface random_complex

  contains 

  subroutine random_scomplex(z)
    ! Generate a random single precision complex within the unit circle.
    complex(4), intent(out) :: z

    real(4)                 :: r, theta

    call random_number(r)
    r = sqrt(r)

    call random_number(theta)
    theta = 2.0 * 3.14159265358979323846 * theta

    z = cmplx(r * cos(theta), r * sin(theta), kind=4)
  end subroutine random_scomplex

  subroutine random_scomplex_rank1(z)
    ! Generate random single-precision complex numbers within the unit circle
    ! for a rank-1 array.
    complex(4), intent(out) :: z(:)

    integer                 :: i

    do i = 1, size(z)
      call random_scomplex(z(i))
    end do

  end subroutine random_scomplex_rank1

  subroutine random_scomplex_rank2(z)
    ! Generate random single-precision complex numbers within the unit circle
    ! for a rank-2 array.
    complex(4), intent(out) :: z(:,:)

    integer                 :: i, j

    do j = 1, size(z, dim=2)
      do i = 1, size(z, dim=1)
        call random_scomplex(z(i, j))
      end do
    end do

  end subroutine random_scomplex_rank2

  subroutine random_dcomplex(z)
    ! Generate a random double precision complex within the unit circle.
    complex(8), intent(out) :: z

    real(8)                 :: r, theta

    call random_number(r)
    r = sqrt(r)

    call random_number(theta)
    theta = 2.0 * 3.14159265358979323846 * theta

    z = cmplx(r * cos(theta), r * sin(theta), kind=8)
  end subroutine random_dcomplex

  subroutine random_dcomplex_rank1(z)
    ! Generate random double-precision complex numbers within the unit circle
    ! for a rank-1 array.
    complex(8), intent(out) :: z(:)

    integer                 :: i

    do i = 1, size(z)
      call random_dcomplex(z(i))
    end do

  end subroutine random_dcomplex_rank1

  subroutine random_dcomplex_rank2(z)
    ! Generate random double-precision complex numbers within the unit circle
    ! for a rank-2 array.
    complex(kind=8), intent(out) :: z(:,:)

    integer                      :: i, j

    do j = 1, size(z, dim=2)
      do i = 1, size(z, dim=1)
        call random_dcomplex(z(i, j))
      end do
    end do

  end subroutine random_dcomplex_rank2

  subroutine random_character(character_array, output)
    ! Randomly select and return a character from an input array of strings.
    character(len=*), intent(in)  :: character_array(:)
    character(len=*), intent(out) :: output

    integer                       :: random_int

    call random_integer(1, size(character_array), random_int)
    output = character_array(random_int)
  end subroutine random_character

  subroutine random_integer(lower_bound, upper_bound, output)
    ! Get a random integer between the lowerBound and the upperBound included.
    integer, intent(in)    :: lower_bound, upper_bound
    integer, intent(out)   :: output
    
    real                   :: rand

    call random_number(rand)
    output = lower_bound + int(rand * (upper_bound - lower_bound + 1))
  end subroutine random_integer

end module module_random
