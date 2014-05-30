module D2Q9Const
  use ISO_FORTRAN_ENV
  integer, parameter :: numQ = 9, latticeDim = 2
  real(kind=real64), parameter :: weights(0:numQ - 1) = &
    (/ 4.0d0/9.0d0,  1.0d0/9.0d0,  1.0d0/9.0d0,  &
       1.0d0/9.0d0,  1.0d0/9.0d0,  1.0d0/36.0d0, &
       1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0 /)
  integer, parameter :: vectors(0:1,0:numQ - 1) = &
    reshape((/0,0, 1,0, 0,1, -1,0, 0,-1, 1,1, -1,1, -1,-1, 1,-1/), shape(vectors))
  real(kind=real64), parameter :: magnitudes(0:numQ - 1) = &
    (/0d0, 1d0, 1d0, 1d0, 1d0, dsqrt(2d0), dsqrt(2d0), dsqrt(2d0), dsqrt(2d0)/)
  integer, parameter :: reverse(0:numQ - 1) = (/0, 3, 4, 1, 2, 7, 8, 5, 6/)
end module D2Q9Const

module D2Q7Const
  use ISO_FORTRAN_ENV
  use simParam, only: pi
  integer, parameter           :: numQ = 7, latticeDim = 2
  real(kind=real64), parameter :: weights(0:numQ - 1) = &
    (/1d0/2d0, 1d0/12d0, 1d0/12d0, 1d0/12d0, 1d0/12d0, 1d0/12d0, 1d0/12d0 /)
  real(kind=real64), parameter :: vectors(0:1,0:numQ - 1) = reshape(       &
    (/ 0d0,0d0, 1d0,0d0, dcos(pi/3),dsin(pi/3), dcos(2*pi/3),dsin(2*pi/3), &
      -1d0,0d0, dcos(4*pi/3),dsin(4*pi/3), dcos(5*pi/3),dsin(5*pi/3) /),   &
      shape(vectors))
  integer, parameter :: reverse(0:numQ - 1) = (/0, 4, 5, 6, 1, 2, 3/)
end module D2Q7Const
