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
  !for the love of god, make sure rDim is even for periodic D2Q7!
  integer, parameter                 :: rDim   = 4
  integer, parameter                 :: cDim   = 3
  integer, parameter                 :: tMax   = 1000
  double precision, parameter        :: deltaX = 0.1d0, deltaT = 0.05d0   ! dT = knudsen #
  double precision, parameter        :: tau    = 0.55d0
  complex(kind=kind(0d0)), parameter        :: lambda = 2d0/(deltaT*(2*tau - 1))
  complex(kind=kind(0d0)), parameter        :: beta   = 2.0d-3                   ! beta = D0
  !complex(kind=kind(0d0)), parameter :: a = dcmplx(0.1d0, 0.0d0)
  !complex(kind=kind(0d0)), parameter :: d = dcmplx(0.025d0, 0.03d0)
end module simParam

program cgle
  use simParam
  use D2Q7Const
  implicit none

  complex(kind=kind(0d0)), allocatable :: f(:,:,:), feq(:,:,:)
  complex(kind=kind(0d0)), allocatable :: rho(:,:), u(:,:,:), uSqr(:,:)
  integer :: r

  allocate(f(rDim,cDim,0:numQ - 1), feq(rDim, cDim,0:numQ - 1))
  allocate(rho(rDim, cDim), u(rDim, cDim, 0:1), uSqr(rDim, cDim))

  f = 0d0
  f(:,:,6) = reshape( (/1,2,3,4,5,6,7,8,9,10,11,12/), shape(f(:,:,2)))
  do r = 1, rdim
    write(*,*) r, "|", f(r,:,6)
  end do
  write(*,*)
  call stream(f)
  do r = 1, rdim
    write(*,*) r, "|", f(r,:,6)
  end do
end program cgle

subroutine computeMacros(f, rho, u, uSqr)
  use D2Q7Const, only: vectors, numQ
  use simParam,  only: rDim, cDim
  implicit none

  complex(kind=kind(0d0)), intent(in)  :: f(rDim, cDim, 0:numQ - 1)
  complex(kind=kind(0d0)), intent(out) :: u(rDim, cDim, 0:1), rho(rDim, cDim), usqr(rDim, cDim)
  integer :: r, c
  do c = 1, cDim
    do r = 1, rDim
      rho(r,c)  = sum(f(r, c, :))
      u(r,c,0)  = sum(f(r,c,:)*vectors(0,:))
      u(r,c,1)  = sum(f(r,c,:)*vectors(1,:))
      uSqr(r,c) = dot_product(u(r,c,:),u(r,c,:))
    end do
  end do
end subroutine computeMacros

subroutine computeFeq(rho, feq)
  use simParam,  only: cDim, rDim, lambda, beta
  use D2Q7Const, only: numQ, latticeDim, soundSpeed
  implicit none

  complex(kind=kind(0d0)), intent(in)  :: rho(rDim,cDim)
  complex(kind=kind(0d0)), intent(out) :: feq(rDim,cDim,0:numQ - 1)
  complex(kind=kind(0d0)) :: tmp(rDim, cDim)
  integer :: dir

  feq(:,:,0) = (1 - lambda*beta*latticeDim/soundSpeed**2)*rho(:,:)

  tmp(:,:) = lambda*beta*latticeDim/(numQ*soundSpeed**2)*rho(:,:)
  do dir = 1, numQ - 1
    feq(:,:,dir) = tmp(:,:)
  end do
end subroutine computeFeq

subroutine stream(f) !Don't dare touch this routine. THE ORDERINGS MATTER!
  use simParam,  only: cDim, rDim
  use D2Q7Const, only: numQ
  implicit none

  complex(kind=kind(0d0)), intent(inout) :: f(rDim,cDim,0:numQ - 1)
  complex(kind=kind(0d0)) :: periodicHor(rDim), periodicVert(cDim)
  integer :: rowIdx

  !!--stream east----------------------------------
  periodicHor(:) = f(:,cDim, 1)
  f(:,2:cDim,1)  = f(:,1:cDim - 1, 1)
  f(:,1, 1)      = periodicHor

  !!--stream northeast----------------------------
  periodicVert(:) = f(1,:,2) !top row of f
  do rowIdx = 2, rDim
    if(modulo(rowIdx,2) .eq. 1) then  !odd rows move up & right
      periodicHor(rowIdx)    = f(cDim,rowIdx,2)
      f(rowIdx - 1,2:cDim,2) = f(rowIdx,1:cDim-1,2)
      f(rowIdx - 1,1,2)      = periodicHor(rowIdx)
    else                              !even rows move straight up
      f(rowIdx - 1,:,2) = f(rowIdx,:,2)
    end if
  end do
  f(rDim,2:cDim,2) = periodicVert(1:cDim - 1)
  f(rDim,1,2)      = periodicVert(cDim)

  !!--stream northwest----------------------------
  periodicVert(:) = f(1,:,3) !top row of f
  do rowIdx = 2, rDim
    if(modulo(rowIdx,2) .eq. 1) then  !odd rows move straight up
      f(rowIdx - 1,:,3) = f(rowIdx,:,3)
    else                              !even rows move up & left
      periodicHor(rowIdx)        = f(rowIdx,1,3)
      f(rowIdx - 1,1:cDim - 1,3) = f(rowIdx,2:cDim,3)
      f(rowIdx - 1,cDim,3)       = periodicHor(rowIdx)
    end if
  end do
  f(rDim,:,3) = periodicVert(:)

  !!--stream west---------------------------------
  periodicHor(:)      = f(:,1,4)
  f(:,1:cDim - 1,4)   = f(:,2:cDim,4)
  f(:,cDim,4)         = periodicHor(:)

  !!--stream southwest----------------------------
  periodicVert(:) = f(rDim,:,5) !bottom row of f
  !this index must run backwards to avoid overwriting data
  do rowIdx = rdim - 1, 1, -1
    if(modulo(rowIdx,2) .eq. 1) then  !odd rows move straight down
      f(rowIdx + 1,:,5) = f(rowIdx,:,5)
    else                              !even rows move down & left
      periodicHor(rowIdx)        = f(rowIdx,1,5)
      f(rowIdx + 1,1:cDim - 1,5) = f(rowIdx,2:cDim,5)
      f(rowIdx + 1,cDim,5)       = periodicHor(rowIdx)
    end if
  end do
  f(1,1:cdim-1,5) = periodicVert(2:cdim)
  f(1,cdim,5) = periodicVert(1)
  
  !!--stream southeast----------------------------
  periodicVert(:) = f(rDim,:,6)
  !this index must run backwards to avoid overwriting data
  do rowIdx = rdim - 1, 1, -1
    if(modulo(rowIdx,2) .eq. 1) then  !odd rows move down & right
      periodicHor(rowIdx)    = f(rowIdx,cDim,6)
      f(rowIdx + 1,2:cDim,6) = f(rowIdx,1:cDim - 1,6)
      f(rowIdx + 1,1,6)      = periodicHor(rowIdx)
    else                              !even rows move straight down
      f(rowIdx + 1,:,6) = f(rowIdx,:,6)
    end if
  end do
  f(1,:,6) = periodicVert(:)
end subroutine stream
