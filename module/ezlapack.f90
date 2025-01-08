!> \brief Global module that specifies the other modules used in ezlapack
!> 
!> This module allows the user to import the `ezlapack` module as a whole, 
!> simplifying the process by not requiring individual imports of other modules.
!> It consolidates the necessary modules used within the `ezlapack` library.
!> 
!> \details
!> The `ezlapack` module includes:
!> - `module_ezmatmul`: A module for matrix multiplication operations.
!> - `module_random`: A module for random number generation and related utilities.
!> 
!> By using `ezlapack`, users can access all the functionality provided by 
!> the underlying modules without needing to import them separately.
module ezlapack
  use module_ezmatmul
  use module_random
end module ezlapack

