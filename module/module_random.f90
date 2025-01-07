!> \file
!> \brief Provides a generic interface for generating random complex numbers, integers, and characters.
!> 
!> \details
!> This module offers utilities for generating random values with specific properties:
!> - Random complex numbers within the unit circle.
!> - Random integers within a user-specified range.
!> - Random selection from an array of characters.
module module_random

  implicit none

  private
  public :: random_complex, random_integer, random_character

  !> \brief Generic interface for generating random complex numbers.
  !>
  !> \details
  !> The `random_complex` interface generates random
  !> complex numbers (single or double precision) within the unit circle.
  !> Supports scalars and arrays up to rank-2.
  !>
  !> ```fortran
  !> call random_complex(z)
  !> ```
  !>
  !> \param[out] z The random complex number(s), which can be:
  !> - A single complex number.
  !> - A rank-1 array of complex numbers.
  !> - A rank-2 array of complex numbers.
  interface random_complex
    module procedure random_scomplex, random_dcomplex
    module procedure random_scomplex_rank1, random_dcomplex_rank1
    module procedure random_scomplex_rank2, random_dcomplex_rank2
  end interface random_complex

  contains

  !> \brief Generate a random single precision complex number within the unit circle.
  !> \param[out] z Single precision complex number output.
  subroutine random_scomplex(z)
    complex(4), intent(out) :: z

    real(4)                 :: r, theta

    call random_number(r)
    r = sqrt(r)

    call random_number(theta)
    theta = 2.0 * 3.14159265358979323846 * theta

    z = cmplx(r * cos(theta), r * sin(theta), kind=4)
  end subroutine random_scomplex

  !> \brief Generate random single-precision complex numbers for a rank-1 array.
  !> \param[out] z Rank-1 array of single precision complex numbers.
  subroutine random_scomplex_rank1(z)
    complex(4), intent(out) :: z(:)

    integer                 :: i

    do i = 1, size(z)
      call random_scomplex(z(i))
    end do

  end subroutine random_scomplex_rank1

  !> \brief Generate random single-precision complex numbers for a rank-2 array.
  !> \param[out] z Rank-2 array of single precision complex numbers.
  subroutine random_scomplex_rank2(z)
    complex(4), intent(out) :: z(:,:)

    integer                 :: i, j

    do j = 1, size(z, dim=2)
      do i = 1, size(z, dim=1)
        call random_scomplex(z(i, j))
      end do
    end do

  end subroutine random_scomplex_rank2

  !> \brief Generate a random double precision complex number within the unit circle.
  !> \param[out] z Double precision complex number output.
  subroutine random_dcomplex(z)
    complex(8), intent(out) :: z

    real(8)                 :: r, theta

    call random_number(r)
    r = sqrt(r)

    call random_number(theta)
    theta = 2.0 * 3.14159265358979323846 * theta

    z = cmplx(r * cos(theta), r * sin(theta), kind=8)
  end subroutine random_dcomplex

  !> \brief Generate random double-precision complex numbers for a rank-1 array.
  !> \param[out] z Rank-1 array of double precision complex numbers.
  subroutine random_dcomplex_rank1(z)
    complex(8), intent(out) :: z(:)

    integer                 :: i

    do i = 1, size(z)
      call random_dcomplex(z(i))
    end do

  end subroutine random_dcomplex_rank1

  !> \brief Generate random double-precision complex numbers for a rank-2 array.
  !> \param[out] z Rank-2 array of double precision complex numbers.
  subroutine random_dcomplex_rank2(z)
    complex(kind=8), intent(out) :: z(:,:)

    integer                      :: i, j

    do j = 1, size(z, dim=2)
      do i = 1, size(z, dim=1)
        call random_dcomplex(z(i, j))
      end do
    end do

  end subroutine random_dcomplex_rank2

  !> \brief Randomly select and return a character from an input array of strings.
  !> \param[in] character_array Input array of strings to select from.
  !> \param[out] output_character Randomly selected character from the input array.
  subroutine random_character(character_array, output_character)
    character(len=*), intent(in)  :: character_array(:)
    character(len=*), intent(out) :: output_character

    integer                       :: random_int

    call random_integer(1, size(character_array), random_int)
    output_character = character_array(random_int)
  end subroutine random_character

  !> \brief Get a random integer between the lower and upper bounds (inclusive).
  !> \param[in] lower_bound Lower bound of the random integer range.
  !> \param[in] upper_bound Upper bound of the random integer range.
  !> \param[out] output_integer Random integer within the specified range.
  subroutine random_integer(lower_bound, upper_bound, output_integer)
    integer, intent(in)    :: lower_bound, upper_bound
    integer, intent(out)   :: output_integer

    real                   :: rand

    call random_number(rand)
    output_integer = lower_bound + int(rand * (upper_bound - lower_bound + 1))
  end subroutine random_integer

end module module_random

