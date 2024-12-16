subroutine random_complex(z)

  ! Generate a random single precision complex
  
  implicit none

  complex, intent(inout) :: z

  real                   :: re, im

  call random_number(re) 
  call random_number(im)

  z = cmplx(re, im)
end subroutine random_complex
