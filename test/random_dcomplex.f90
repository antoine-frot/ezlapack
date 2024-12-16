subroutine random_dcomplex(z)

  ! Generate a random single precision complex
  
  implicit none

  complex*16, intent(inout) :: z

  double precision          :: re, im

  call random_number(re) 
  call random_number(im)

  z = cmplx(re, im)
end subroutine random_dcomplex
