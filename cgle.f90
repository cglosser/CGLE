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
  integer, parameter                 :: xDim   = 5
  !for the love of god, make sure ydim is even for periodic D2Q7!
  integer, parameter                 :: yDim   = 6
  integer, parameter                 :: tMax   = 1000
  double precision, parameter        :: deltaX = 0.1d0, deltaT = 0.05d0   ! dT = knudsen #
  double precision, parameter        :: tau    = 0.55d0
  double precision, parameter        :: lambda = 2d0/(deltaT*(2*tau - 1))
  double precision, parameter        :: beta   = 2.0d-3                   ! beta = D0
  !double precision, parameter :: a = dcmplx(0.1d0, 0.0d0)
  !double precision, parameter :: d = dcmplx(0.025d0, 0.03d0)
end module simParam

program cgle
  use simParam
  use D2Q7Const
  implicit none

  double precision, allocatable :: f(:,:,:), feq(:,:,:)
  double precision, allocatable :: rho(:,:), u(:,:,:), uSqr(:,:)
  integer :: t,y

  allocate(f(0:numQ -1, ydim, xdim), feq(0:numQ - 1, ydim, xdim))
  allocate(rho(ydim, xdim), u(0:1, ydim, xdim), uSqr(ydim, xdim))

  f = 0
  f(1,1,1) = 8
  do t = 1, 10
    call stream(f)
    call computeMacros(f,rho,u,uSqr)
    do y = 1, ydim
      write(*,*) f(1,y,:)
    end do
    write(*,*)
  end do
end program cgle

subroutine computeMacros(f, rho, u, uSqr)
  use D2Q7Const, only: vectors, numQ
  use simParam,  only: xDim, yDim
  implicit none

  double precision, intent(in)  :: f(0:numQ - 1, ydim, xdim)
  double precision, intent(out) :: u(0:1, ydim, xdim), rho(ydim, xdim), usqr(ydim, xdim)
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

subroutine computeFeq(rho, feq)
  use simParam,  only: xdim, ydim, lambda, beta
  use D2Q7Const, only: numQ, latticeDim, soundSpeed
  implicit none

  double precision, intent(in)  :: rho(ydim,xdim)
  double precision, intent(out) :: feq(0:numQ - 1,ydim,xdim)
  double precision :: tmp(ydim, xdim)
  integer :: dir

  tmp(:,:) = lambda*beta*latticeDim/(numQ*soundSpeed**2)*rho(:,:)
  do dir = 1, numQ - 1
    feq(dir,:,:) = tmp(:,:)
  end do
  feq(0,:,:)  = (1 - lambda*beta*latticeDim/soundSpeed**2)*rho(:,:)
end subroutine computeFeq

subroutine stream(f)
  use simParam,  only: xdim, ydim
  use D2Q7Const, only: numQ
  implicit none

  double precision, intent(inout) :: f(0:numQ - 1,ydim,xdim)
  double precision :: horizontalBuffer(ydim), verticalBuffer(xdim)
  integer :: rowIdx

  !--stream east----------------------------------
  horizontalBuffer(:) = f(1,:,xdim)
  f(1,:,2:xdim)       = f(1,:,1:xdim - 1)
  f(1,:,1)            = horizontalBuffer(:)

  !--stream northeast----------------------------
  verticalBuffer(:) = f(2,1,:) !top row of f
  do rowIdx = 2, ydim
    if(modulo(rowIdx,2) .eq. 1) then
      horizontalBuffer(rowIdx) = f(2,rowIdx,xdim)
      f(2,rowIdx - 1,2:xdim)   = f(2,rowIdx,1:xdim-1)
      f(2,rowIdx - 1,1)        = horizontalBuffer(rowIdx)
    else
      f(2,rowIdx - 1,:) = f(2,rowIdx,:)
    end if
  end do
  f(2,ydim,2:xdim) = verticalBuffer(1:xdim - 1)
  f(2,ydim,1)      = verticalBuffer(xdim)

  !--stream northwest----------------------------
  verticalBuffer(:) = f(3,1,:) !top row of f
  do rowIdx = 2, ydim
    if(modulo(rowIdx,2) .eq. 1) then
      f(3, rowIdx - 1,:) = f(3,rowIdx,:)
    else
      horizontalBuffer(rowIdx)   = f(3,rowIdx,1)
      f(3,rowIdx - 1,1:xdim - 1) = f(3,rowIdx,2:xdim)
      f(3,rowIdx - 1,xdim)       = horizontalBuffer(rowIdx)
    end if
  end do
  f(3,ydim,:) = verticalBuffer(:)

  !--stream west---------------------------------
  horizontalBuffer(:) = f(4,:,1)
  f(4,:,1:xdim - 1)   = f(4,:,2:xdim)
  f(4,:,xdim)         = horizontalBuffer(:)

  !--stream southwest----------------------------
  verticalBuffer(:) = f(5,ydim,:) !bottom row of f
  do rowIdx = 1, ydim - 1
    if(modulo(rowIdx,2) .eq. 1) then
      f(5,rowIdx + 1,:) = f(5,rowIdx,:)
    else
      horizontalBuffer(rowIdx)   = f(5,rowIdx,1)
      f(5,rowIdx + 1,1:xdim - 1) = f(5,rowIdx,2:xdim)
      f(5,rowIdx + 1,xdim)       = horizontalBuffer(rowIdx)
    end if
  end do
  f(5,1,:) = verticalBuffer(:)
  
  !--stream southeast----------------------------
  verticalBuffer(:) = f(6,ydim,:)
  do rowIdx = 1, ydim - 1
    if(modulo(rowIdx,2) .eq. 1) then
      horizontalBuffer(rowIdx) = f(6,rowIdx,xdim)
      f(6,rowIdx + 1,2:xdim)   = f(6,rowIdx,1:xdim - 1)
      f(6,rowIdx + 1,1)        = horizontalBuffer(rowIdx)
    else
      f(6,rowIdx + 1,:) = f(6,rowIdx,:)
    end if
  end do
  f(6,1,:) = verticalBuffer(:)
end subroutine stream
