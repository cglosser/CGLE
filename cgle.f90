module D2Q9Const
  double precision, parameter :: pi = 4*datan(1d0)
  integer, parameter :: numQ = 9
  double precision, parameter :: weights(0:numQ - 1) = &
    (/ 4.0d0/9.0d0,  1.0d0/9.0d0,  1.0d0/9.0d0,  &
       1.0d0/9.0d0,  1.0d0/9.0d0,  1.0d0/36.0d0, &
       1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0 /)
  integer, parameter :: vectors(0:1,0:numQ - 1) = &
    reshape((/0,0, 1,0,0, 1,-1,0,0,-1, 1, 1,-1, 1,-1,-1, 1,-1/), shape(vectors))
  double precision, parameter :: magnitudes(0:numQ - 1) = &
    (/0d0, 1d0, 1d0, 1d0, 1d0, dsqrt(2d0), dsqrt(2d0), dsqrt(2d0), dsqrt(2d0)/)
  integer, parameter :: reverse(0:numQ - 1) = (/0, 3, 4, 1, 2, 7, 8, 5, 6/)
  integer, parameter :: latticeDim = 2
end module D2Q9Const

module D2Q7Const
  double precision, parameter :: pi = 4*datan(1d0)
  integer, parameter          :: numQ = 7
  double precision, parameter :: weights(0:numQ - 1) = &
    (/1d0/2d0, 1d0/12d0, 1d0/12d0, 1d0/12d0, 1d0/12d0, 1d0/12d0, 1d0/12d0 /)
  double precision, parameter :: vectors(0:1,0:numQ - 1) = reshape(        &
    (/ 0d0,0d0, 1d0,0d0, dcos(pi/3), dsin(pi/3), dcos(2*pi/3), dsin(2*pi/3), &
      -1d0,0d0, dcos(4*pi/3), dsin(4*pi/3), dcos(5*pi/3), dsin(5*pi/3) /),   &
      shape(vectors))
  double precision, parameter :: soundSpeed = 2d0
  integer, parameter :: reverse(0:numQ - 1) = (/0, 4, 5, 6, 1, 2, 3/)
  integer, parameter :: latticeDim = 2
  
end module D2Q7Const

module simParam
  !for the love of god, make sure rDim is even for periodic D2Q7!
  integer, parameter                 :: rDim   = 231
  integer, parameter                 :: cDim   = 200
  integer, parameter                 :: tMax   = 1000
  double precision, parameter        :: deltaX = 0.1d0, deltaT = 0.05d0   ! dT = knudsen #
  double precision, parameter        :: tau    = 0.55d0
  double precision, parameter        :: t0_coef= 0.3d0
  double precision, parameter        :: boxLength = 10d0
  complex(kind=kind(0d0)), parameter :: lambda = 2d0/(deltaT*(2*tau - 1))
  complex(kind=kind(0d0)), parameter :: beta   = 2.0d-3                   ! beta = D0
  complex(kind=kind(0d0)), parameter :: a = dcmplx(0.100d0, 0.24d0)
  complex(kind=kind(0d0)), parameter :: d = dcmplx(0.025d0, 0.03d0)
end module simParam

program cgle
  use simParam
  use D2Q7Const
  implicit none

  complex(kind=kind(0d0)), allocatable :: f(:, :, :), feq(:, :, :), rho(:, :), omega(:, :)
  double precision :: sumF = 0d0
  integer :: time, r, c

  allocate(f(rDim, cDim,0:numQ - 1), feq(rDim, cDim,0:numQ - 1))
  allocate(rho(rDim, cDim), omega(rDim, cDim))

  open(unit=20,file="realpart.dat")
  open(unit=30,file="imagpart.dat")

  call setInitialF(f)
  call computeMacros(f, rho, omega)
  do time = 1, tMax
    write(*,*) time, sumF
    call computeMacros(f, rho, omega)
    call computeFeq(rho, feq)
    call collide(feq, omega, f)
    call stream(f)

    sumF = sum(real(f)**2 + aimag(f)**2)
    f = f/(sqrt(sumF/(rDim*cDim*numQ)))
    if(mod(time, 500) .eq. 0) then
      do c = 1, cDim
        write(20,*) real(rho(:,c))
        write(30,*) aimag(rho(:,c))
      end do
    end if
  end do

end program cgle

subroutine computeMacros(f, rho, omega)
  use D2Q7Const, only: vectors, numQ
  use simParam,  only: rDim, cDim, deltaT, a, d
  implicit none

  complex(kind=kind(0d0)), intent(in)    :: f(rDim, cDim, 0:numQ - 1)
  complex(kind=kind(0d0)), intent(out)   :: rho(rDim, cDim), omega(rDim, cDim)
  complex(kind=kind(0d0)) :: hamiltonian
  integer :: r, c
  do c = 1, cDim
    do r = 1, rDim
      rho(r, c)   = sum(f(r, c, :))
      hamiltonian = (a - d*rho(r,c)*conjg(rho(r,c)))*rho(r,c)
      omega(r,c)  = deltaT*hamiltonian/numQ
    end do
  end do
end subroutine computeMacros

subroutine computeFeq(rho, feq)
  use simParam,  only: cDim, rDim, lambda, beta
  use D2Q7Const, only: numQ, latticeDim, soundSpeed
  implicit none

  complex(kind=kind(0d0)), intent(in)  :: rho(rDim, cDim)
  complex(kind=kind(0d0)), intent(out) :: feq(rDim, cDim,0:numQ - 1)
  complex(kind=kind(0d0)) :: tmp(rDim, cDim)
  integer :: dir

  feq(:, :,0) = (1 - lambda*beta*latticeDim/soundSpeed**2)*rho(:, :)

  tmp(:, :) = lambda*beta*latticeDim/(numQ*soundSpeed**2)*rho(:, :)
  do dir = 1, numQ - 1
    feq(:, :, dir) = tmp(:, :)
  end do
end subroutine computeFeq

subroutine stream(f) 
  ! To transform the sqare f into the hexagonal lattice it represents, 
  ! the following code is written as if even-indexed columns shift down. I.e.:
  !
  ! 1 5 9 D               1   9     Directions: 
  !                         5   D                1
  ! 2 6 A E    becomes    2   A                2   6
  !                         6   E                0
  ! 3 7 B F    ------>    3   B                3   5 
  !                         7   F                4
  ! 4 8 C G               4   C                      (0 to any equidistant)
  !                         8   G
  !
  ! Thus, the matrix must have an even number of columns, translations "up/down" 
  ! require only a circular shift, and translation in hybrid directions require
  ! special handling of every other column.
  use simParam,  only: cDim, rDim
  use D2Q7Const, only: numQ
  implicit none

  complex(kind=kind(0d0)), intent(inout) :: f(rDim, cDim,0:numQ - 1)
  integer :: colIdx

  !!--stream along 1------------------------------
  f(:,:,1) = cshift(f(:,:,1), 1, dim=1)

  !!--stream along 2------------------------------
  f(:,:,2) = cshift(f(:,:,2), 1, dim=2)
  do colIdx = 2, cDim, 2
    f(:, colIdx, 2) = cshift(f(:, colIdx, 2), 1, dim=1)
  end do

  !!--stream along 3------------------------------
  f(:,:,3) = cshift(f(:,:,3), 1, dim=2)
  do colIdx = 1, cDim - 1, 2
    f(:, colIdx, 3) = cshift(f(:, colIdx, 3), -1, dim=1)
  end do

  !!--stream along 4------------------------------
  f(:,:,4) = cshift(f(:,:,4), -1, dim=1)

  !!--stream along 5------------------------------
  f(:,:,5) = cshift(f(:,:,5), -1, dim=2)
  do colIdx = 1, cDim - 1, 2
    f(:, colIdx, 5) = cshift(f(:, colIdx, 5), -1, dim=1)
  end do

  !!--stream along 6------------------------------
  f(:,:,6) = cshift(f(:,:,6), -1, dim=2)
  do colIdx = 2, cDim, 2
    f(:, colIdx, 6) = cshift(f(:, colIdx, 6), 1, dim=1)
  end do
end subroutine stream

subroutine collide(fEq, omega, f)
  use D2Q7Const, only: numQ
  use simParam,  only: rDim, cDim, tau 
  implicit none

  complex(kind=kind(0d0)), intent(inout) :: f(rDim, cDim, 0:numQ - 1)
  complex(kind=kind(0d0)), intent(in)    :: fEq(rDim, cDim, 0:numQ - 1), &
                                            omega(rDim, cDim)
  
  integer :: dir
  do dir = 0, numQ - 1
    f(:,:,dir) = f(:,:,dir) - 1/tau*(f(:,:,dir) - fEq(:,:,dir)) + omega
  end do

  f(:,1,:)     = f(:,2,:)
  f(:, cDim,:) = f(:, cDim - 1,:)
  f(1,:,:)     = f(2,:,:)
  f(rDim,:,:)  = f(rDim - 1, :,:)
end subroutine collide

subroutine setInitialF(f)
  use D2Q7Const, only: numQ
  use simParam,  only: rDim, cDim, t0_coef, deltaX, boxLength
  implicit none

  complex(kind=kind(0d0)), intent(out) :: f(rDim, cDim, 0:numQ - 1)
  integer :: rIdx, cIdx
  double precision :: x, y
  do cIdx = 1, cDim
    x = cIdx*deltaX - boxLength
    do rIdx = 1, rDim
      y = rIdx*deltaX - boxLength
      f(rIdx, cIdx, :) = t0_coef*dcmplx(x, y)/numQ
    end do
  end do
end subroutine setInitialF
