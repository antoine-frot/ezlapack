module random
  !------------------------------------------------!
  ! Purpose: Generic interface that generates      !
  ! a complex number of single or double precision.!
  !------------------------------------------------!

  implicit none

  private
  public :: random_complex, random_integer, random_character
  
  interface random_complex
    module procedure random_scomplex, random_dcomplex
  end interface random_complex

  subroutine random_scomplex(z)
    ! Generate a random single precision complex within the unit circle
    complex, intent(inout) :: z

    real                   :: r, theta

    call random_number(r)
    r = sqrt(r)

    call random_number(theta)
    theta = 2.0 * 3.14159265358979323846 * theta

    z = cmplx(r * cos(theta), r * sin(theta))
  end subroutine random_scomplex

  subroutine random_dcomplex(z)
    ! Generate a random double precision complex within the unit circle
    complex*16, intent(inout) :: z

    real(8)                   :: r, theta

    call random_number(r)
    r = sqrt(r)

    call random_number(theta)
    theta = 2.0 * 3.14159265358979323846 * theta

    z = cmplx(r * cos(theta), r * sin(theta))
  end subroutine random_dcomplex

  subroutine random_character(character_array, output)
    ! Get a random character of an array of characters
    character(len=*), intent(in)    :: character_array(:)
    character(len=*), intent(inout) :: output

    integer                         :: random_int

    call random_integer(1, size(character_array), random_int)
    output = character_array(random_int)
  end subroutine random_character

  subroutine random_integer(lower_bound, upper_bound, output)
    ! Get a random integer between the lowerBound and the upperBound included
    integer, intent(in)    :: lower_bound, upper_bound
    integer, intent(inout) :: output
    
    real                   :: rand

    call random_number(rand)
    output = lower_bound + int(rand * (upper_bound - lower_bound + 1))
  end subroutine random_integer

end module random
