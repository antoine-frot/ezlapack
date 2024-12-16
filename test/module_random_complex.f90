module random_complex
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

end module random_complex
