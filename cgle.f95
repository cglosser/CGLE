program cgle
  use ISO_FORTRAN_ENV
  use simParam
  use D2Q9Const
  use CGLEvis
  implicit none

  integer :: time
  real(kind=real64) :: f_rsq = 0d0
  complex(kind=real64), allocatable :: f_density(:, :, :), feq(:, :, :)
  complex(kind=real64), allocatable :: rho(:, :), omega(:, :)

  allocate(f_density(rDim, cDim, 0:numQ - 1), feq(rDim, cDim, 0:numQ - 1))
  allocate(rho(rDim, cDim), omega(rDim, cDim))

  open(unit=20,file="realpart.dat")
  open(unit=30,file="imagpart.dat")

  call plot_init()

  !call initRandomF(f_density)
  call initSpiralF(f_density)
  do time = 1, tMax
    f_rsq = sum(real(f_density)**2 + aimag(f_density)**2)
    f_density = f_density/sqrt(normalization*f_rsq/size(f_density))

    call computeMacros(f_density, rho, omega)
    call computeFeq(rho, feq)
    call collide(feq, omega, f_density)
    call stream(f_density)

    if(mod(time, 5) .eq. 0) then
      !do c = 1, cDim
        !write(20,*) real(rho(:,c))
        !write(30,*) aimag(rho(:,c))
      !end do
      call plot_array(real(rho(:,:)))
    end if

    write(*,*) time, f_rsq
  end do

  call plot_close()

  close(20)
  close(30)

  deallocate(f_density, feq, rho, omega)
end program cgle

subroutine computeMacros(f_density, rho, omega)
  use ISO_FORTRAN_ENV
  use D2Q9Const, only: vectors, numQ
  use simParam,  only: rDim, cDim, deltaT, a, d
  implicit none

  complex(kind=real64), intent(in)    :: f_density(rDim, cDim, 0:numQ - 1)
  complex(kind=real64), intent(out)   :: rho(rDim, cDim), &
                                         omega(rDim, cDim)
  complex(kind=real64) :: hamiltonian(rDim, cDim)
  
  rho = sum(f_density, dim = 3)
  hamiltonian = (a - d*rho*conjg(rho))*rho
  omega = deltaT*hamiltonian/numQ
end subroutine computeMacros

subroutine computeFeq(rho, feq)
  use ISO_FORTRAN_ENV
  use simParam,  only: cDim, rDim, lambda, beta, soundSpeed
  use D2Q9Const, only: numQ, latticeDim
  implicit none

  complex(kind=real64), intent(in)  :: rho(rDim, cDim)
  complex(kind=real64), intent(out) :: feq(rDim, cDim,0:numQ - 1)
  complex(kind=real64) :: tmp(rDim, cDim)
  integer :: dir

  feq(:, :, 0) = (1 - 3*lambda*beta*latticeDim/(4*soundSpeed**2))*rho(:, :)

  tmp(:, :) = lambda*beta*latticeDim/((numQ - 1)*soundSpeed**2)*rho(:, :)
  do dir = 1, numQ - 1
    if (dir .lt. 5) then
      feq(:, :, dir) = tmp(:, :)
    else
      feq(:, :, dir) = tmp(:, :)/2d0
    end if
  end do
end subroutine computeFeq

subroutine stream(f_density) 
  use ISO_FORTRAN_ENV
  use simParam,  only: cDim, rDim
  use D2Q9Const, only: numQ
  implicit none

  complex(kind=real64), intent(inout) :: f_density(rDim, cDim,0:numQ - 1)

  ! Neumann boundary conditions
  f_density(2:rDim - 1, 1, :)    = f_density(2:rDim - 1, 2, :)
  f_density(2:rDim - 1, cDim, :) = f_density(2:rDim - 1, cDim - 1, :)
  f_density(1, :, :)    = f_density(2, :, :)
  f_density(rdim, :, :) = f_density(rdim - 1, :, :)

  !!--stream along 1------------------------------
  f_density(:,:,1) = cshift(f_density(:,:,1), -1, dim=1)

  !!--stream along 2------------------------------
  f_density(:,:,2) = cshift(f_density(:,:,2),  1, dim=2)

  !!--stream along 3------------------------------
  f_density(:,:,3) = cshift(f_density(:,:,3),  1, dim=1)

  !!--stream along 4------------------------------
  f_density(:,:,4) = cshift(f_density(:,:,4), -1, dim=2)

  !!--stream along 5------------------------------
  f_density(:,:,5) = cshift(f_density(:,:,5), -1, dim=1)
  f_density(:,:,5) = cshift(f_density(:,:,5),  1, dim=2)

  !!--stream along 6------------------------------
  f_density(:,:,6) = cshift(f_density(:,:,6),  1, dim=1)
  f_density(:,:,6) = cshift(f_density(:,:,6),  1, dim=2)

  !!--stream along 7------------------------------
  f_density(:,:,7) = cshift(f_density(:,:,7),  1, dim=1)
  f_density(:,:,7) = cshift(f_density(:,:,7), -1, dim=2)

  !!--stream along 7------------------------------
  f_density(:,:,8) = cshift(f_density(:,:,8), -1, dim=1)
  f_density(:,:,8) = cshift(f_density(:,:,8), -1, dim=2)
end subroutine stream

subroutine collide(fEq, omega, f_density)
  use ISO_FORTRAN_ENV
  use D2Q9Const, only: numQ
  use simParam,  only: rDim, cDim, tau 
  implicit none

  complex(kind=real64), intent(inout) :: f_density(rDim, cDim, 0:numQ - 1)
  complex(kind=real64), intent(in)    :: fEq(rDim, cDim, 0:numQ - 1), &
                                         omega(rDim, cDim)
  
  integer :: dir
  do dir = 0, numQ - 1
    f_density(:,:,dir) = f_density(:,:,dir) - 1/tau*(f_density(:,:,dir) - fEq(:,:,dir)) + omega(:, :)
  end do
end subroutine collide

subroutine initSpiralF(f_density)
  use ISO_FORTRAN_ENV
  use D2Q9Const, only: numQ
  use simParam,  only: rDim, cDim, t0_coef, deltaX, boxLength
  implicit none

  complex(kind=real64), intent(out) :: f_density(rDim, cDim, 0:numQ - 1)
  integer :: rIdx, cIdx
  real(kind=real64) :: x, y
  do cIdx = 1, cDim
    x = (cIdx - 1 - cDim/2)*deltaX
    do rIdx = 1, rDim
      y = (rIdx - 1 - rDim/2)*deltaX
      f_density(rIdx, cIdx, :) = t0_coef*dcmplx(x, y)/numQ
    end do
  end do
end subroutine initSpiralF

subroutine initRandomF(f_density)
  use ISO_FORTRAN_ENV
  use D2Q9Const, only: numQ
  use simParam,  only: rDim, cDim, t0_coef, deltaX, boxLength
  implicit none
  
  complex(kind=real64), intent(out) :: f_density(rDim, cDim, 0:numQ - 1)
  integer :: rIdx, cIdx
  real(kind=real64) :: rands(2)
  do cIdx = 1, cDim
    do rIdx = 1, rDim
      call random_number(rands)
      rands = rands - 0.5d0
      f_density(rIdx, cIdx, :) = t0_coef*dcmplx(rands(1), rands(2))/numQ
    end do
  end do
end subroutine initRandomF
