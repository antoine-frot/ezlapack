module random
  !------------------------------------------------!
  ! Purpose: Generic interface that generates      !
  ! a complex number of single or double precision.!
  !------------------------------------------------!

  interface random_complex
    module procedure random_scomplex, randome_dcomplex
  end interface random_complex

  subroutine random_scomplex(z)
    ! Generate a random single precision complex
    implicit none

    complex, intent(inout) :: z
    real                   :: re, im

    call random_number(re) 
    call random_number(im)

    z = cmplx(re, im)
  end subroutine random_scomplex

  subroutine random_dcomplex(z)
    ! Generate a random double precision complex
    implicit none

    complex*16, intent(inout) :: z
    double precision          :: re, im

    call random_number(re) 
    call random_number(im)

    z = cmplx(re, im)
  end subroutine random_dcomplex

  subroutine random_character(character_array, output)
    ! Get a random character of an array of characters
    implicit none

    character(len=:), intent(in)    :: character_array
    character(len=*), intent(inout) :: output
    integer                         :: lengh
    integer                         :: random_int
    real                            :: random

    lengh = size(character_array)
    call random_integer(random_int)
    output = character_array
  end subroutine random_character

  subroutine random_integer(lower_bound, upper_bound, output)
    ! Get a random integer between the lowerBound and the upperBound included
    implicit none

    integer, intent(in)    :: lower_bound, upper_bound
    integer, intent(inout) :: output
    
    real                   :: rand

    call random_number(rand)
    output = lower_bound + int(rand * (upper_bound - lower_bound + 1))
  end subroutine random_integer

end module random
