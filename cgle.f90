module D2Q9Const
  double precision, parameter :: pi = 4*datan(1d0)
  integer, parameter :: numQ = 9
  double precision, parameter :: weights(0:numQ - 1) = &
    (/ 4.0d0/9.0d0,  1.0d0/9.0d0,  1.0d0/9.0d0,  &
       1.0d0/9.0d0,  1.0d0/9.0d0,  1.0d0/36.0d0, &
       1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0 /)
  integer, parameter :: vectors(0:1,0:numQ - 1) = &
    reshape((/0,0,1,0,0,1,-1,0,0,-1,1,1,-1,1,-1,-1,1,-1/), shape(vectors))
  double precision, parameter :: magnitudes(0:numQ - 1) = &
    (/0d0,1d0,1d0,1d0,1d0,dsqrt(2d0),dsqrt(2d0),dsqrt(2d0),dsqrt(2d0)/)
  integer, parameter :: reverse(0:numQ - 1) = (/0,3,4,1,2,7,8,5,6/)
  integer, parameter :: latticeDim = 2
end module D2Q9Const

module D2Q7Const
  double precision, parameter :: pi = 4*datan(1d0)
  integer, parameter          :: numQ = 7
  double precision, parameter :: weights(0:numQ - 1) = &
    (/1d0/2d0,1d0/12d0,1d0/12d0,1d0/12d0,1d0/12d0,1d0/12d0,1d0/12d0 /)
  double precision, parameter :: vectors(0:1,0:numQ - 1) = reshape(        &
    (/ 0d0,0d0, 1d0,0d0, dcos(pi/3),dsin(pi/3), dcos(2*pi/3),dsin(2*pi/3), &
      -1d0,0d0, dcos(4*pi/3),dsin(4*pi/3), dcos(5*pi/3),dsin(5*pi/3) /),   &
      shape(vectors))
  double precision, parameter :: soundSpeed = 2d0
  integer, parameter :: reverse(0:numQ - 1) = (/0,4,5,6,1,2,3/)
  integer, parameter :: latticeDim = 2
  
end module D2Q7Const

module simParam
  integer, parameter          :: xDim   = 5
  integer, parameter          :: yDim   = 5
  integer, parameter          :: tMax   = 1000
  double precision, parameter :: deltaX = 0.1d0, deltaT = 0.05d0   ! dT = knudsen #
  double precision, parameter :: tau    = 0.55d0
  double precision, parameter :: lambda = 2d0/(deltaT*(2*tau - 1))
  double precision, parameter :: beta   = 2.0d-3                   ! beta = D0
end module simParam

program cgle
  use simParam
  use D2Q7Const
  implicit none

  complex(kind=kind(0d0)), allocatable :: f(:,:,:), feq(:,:,:)
  complex(kind=kind(0d0)), allocatable :: rho(:,:), u(:,:,:), uSqr(:,:)

  allocate(f(0:numQ -1, ydim, xdim), feq(0:numQ - 1, ydim, xdim))
  allocate(rho(ydim, xdim), u(0:1, ydim, xdim), uSqr(ydim, xdim))

  f = dcmplx(0d0, 0d0)
  call computeMacros(f, rho, u, uSqr)

end program cgle

subroutine computeMacros(f, rho, u, uSqr)
  use D2Q7Const, only: vectors, numQ
  use simParam,  only: xDim, yDim
  implicit none

  complex(kind=kind(0d0)), intent(in)  :: f(0:numQ - 1, ydim, xdim)
  complex(kind=kind(0d0)), intent(out) :: u(0:1, ydim, xdim), rho(ydim, xdim), usqr(ydim, xdim)
  integer :: x, y
  do x = 1, xdim
    do y = 1, ydim
      rho(y,x)  = sum(f(:, y, x))
      u(0,y,x)  = sum(f(:,y,x)*vectors(0,:))
      u(1,y,x)  = sum(f(:,y,x)*vectors(1,:))
      uSqr(y,x) = dot_product(u(:,y,x),u(:,y,x))
    end do
  end do
end subroutine computeMacros
