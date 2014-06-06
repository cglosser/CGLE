program cgle
  use ISO_FORTRAN_ENV
  use simParam
  use D2Q9Const
  use CGLEvis
  implicit none

  integer :: time
  real(kind=real64) :: f_rsq = 0d0, laserMask(rDim, cDim) = 0d0
  complex(kind=real64), allocatable :: f_density(:, :, :, :), feq(:, :, :, :)
  complex(kind=real64), allocatable :: psi(:,:, :), omega(:,:, :)

  allocate(f_density(rDim, cDim, 0:numQ - 1, 2), feq(rDim, cDim, 0:numQ - 1, 2))
  allocate(psi(rDim, cDim, 2), omega(rDim, cDim, 2))

  call plot_init()

  !call initRandomF(f_density)
  call initSpiralF(f_density)
  do time = 1, tMax
    call neumannBC(f_density)

    call computeMacros(f_density, laserMask, psi, omega)
    call computeFeq(psi, feq)
    call collide(feq, omega, f_density)
    call stream(f_density)

    if(mod(time, 5) .eq. 0) call plot_array(real(psi(:, :, cavity)))

    write(*,*) time, f_rsq
  end do

  call plot_close()

  deallocate(f_density, feq, psi, omega)
end program cgle

subroutine computeMacros(f_density, laserMask, psi, omega)
  use ISO_FORTRAN_ENV
  use D2Q9Const, only: vectors, numQ
  use simParam,  only: rDim, cDim, deltaT, a, d, &
                       cavity, exciton, omegaRabi
  implicit none

  complex(kind=real64), intent(in)    :: f_density(rDim, cDim, 0:numQ - 1, 2)
  real(kind=real64), intent(in)       :: laserMask(rDim, cDim)
  complex(kind=real64), intent(out)   ::       psi(rDim, cDim, 2), &
                                             omega(rDim, cDim, 2)
  complex(kind=real64) :: hamiltonian(rDim, cDim, 2)
  
  psi = sum(f_density, dim = 3)
  hamiltonian(:, :, exciton) = &
    (a(exciton) - d(exciton)*abs(psi(:, :, exciton))**2) + omegaRabi/2*psi(:, :,  cavity)
  hamiltonian(:, :, cavity)  = &
    ( a(cavity) -  d(cavity)*abs(psi(:, :, cavity))**2)  + omegaRabi/2*psi(:, :, exciton) + laserMask(:,:)
  omega = deltaT*hamiltonian/numQ
end subroutine computeMacros

subroutine computeFeq(psi, feq)
  use ISO_FORTRAN_ENV
  use simParam,  only: cDim, rDim, lambda, beta, latticeVelocity, &
                       cavity, exciton
  use D2Q9Const, only: numQ, latticeDim
  implicit none

  complex(kind=real64), intent(in)  :: psi(rDim, cDim, 2)
  complex(kind=real64), intent(out) :: feq(rDim, cDim, 0:numQ - 1, 2)
  complex(kind=real64) :: tmp(rDim, cDim, 2), lbmTerm(2)
  integer :: dir

  lbmTerm(:) = lambda*beta(:)*latticeDim/latticeVelocity**2

  feq(:, :, 0, cavity)  = (1 - 3*lbmTerm( cavity)/4)*psi(:, :, cavity)
  feq(:, :, 0, exciton) = (1 - 3*lbmTerm(exciton)/4)*psi(:, :, exciton)

  tmp(:, :, cavity)     = lbmTerm(cavity )/(numQ - 1)*psi(:, :, cavity)
  tmp(:, :, exciton)    = lbmTerm(exciton)/(numQ - 1)*psi(:, :, exciton)

  do dir = 1, 4
    feq(:, :, dir, :)     = tmp(:, :, :)
    feq(:, :, dir + 4, :) = tmp(:, :, :)/2
  end do
end subroutine computeFeq

subroutine stream(f_density) 
  use ISO_FORTRAN_ENV
  use simParam,  only: cDim, rDim
  use D2Q9Const, only: numQ
  implicit none

  complex(kind=real64), intent(inout) :: f_density(rDim, cDim,0:numQ - 1, 2)

  !!--stream along 1------------------------------
  f_density(:,:,1,:) = cshift(f_density(:,:,1,:), -1, dim=1)

  !!--stream along 2------------------------------
  f_density(:,:,2,:) = cshift(f_density(:,:,2,:),  1, dim=2)

  !!--stream along 3------------------------------
  f_density(:,:,3,:) = cshift(f_density(:,:,3,:),  1, dim=1)

  !!--stream along 4------------------------------
  f_density(:,:,4,:) = cshift(f_density(:,:,4,:), -1, dim=2)

  !!--stream along 5------------------------------
  f_density(:,:,5,:) = cshift(cshift(f_density(:,:,5,:), -1, dim=1),  1, dim=2)

  !!--stream along 6------------------------------
  f_density(:,:,6,:) = cshift(cshift(f_density(:,:,6,:),  1, dim=1),  1, dim=2)

  !!--stream along 7------------------------------
  f_density(:,:,7,:) = cshift(cshift(f_density(:,:,7,:),  1, dim=1), -1, dim=2)

  !!--stream along 7------------------------------
  f_density(:,:,8,:) = cshift(cshift(f_density(:,:,8,:), -1, dim=1), -1, dim=2)
end subroutine stream

subroutine collide(fEq, omega, f_density)
  use ISO_FORTRAN_ENV
  use D2Q9Const, only: numQ
  use simParam,  only: rDim, cDim, tau 
  implicit none

  complex(kind=real64), intent(inout) :: f_density(rDim, cDim, 0:numQ - 1, 2)
  complex(kind=real64), intent(in)    :: fEq(rDim, cDim, 0:numQ - 1, 2), &
                                         omega(rDim, cDim, 2)

  integer :: dir
  do dir = 0, numQ - 1
    f_density(:, :, dir, :) = f_density(:, :, dir, :) - &
      (f_density(:,:, dir, :) - fEq(:, :, dir, :))/tau + omega(:, :, :)
  end do
end subroutine collide

subroutine neumannBC(f_density)
  use ISO_FORTRAN_ENV
  use D2Q9Const, only: numQ
  use simParam,  only: rDim, cDim
  implicit none

  complex(kind=real64), intent(inout) :: f_density(rDim, cDim, 0:numQ - 1, 2)

  f_density(2:rDim - 1,    1, :, :) = f_density(2:rDim - 1,        2, :, :)
  f_density(2:rDim - 1, cDim, :, :) = f_density(2:rDim - 1, cDim - 1, :, :)
  f_density(   1, :, :, :)    = f_density(       2, :, :, :)
  f_density(rdim, :, :, :) = f_density(rdim - 1, :, :, :)
end subroutine neumannBC

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

function gaussian(mu, sigma, x) result(output)
  use ISO_FORTRAN_ENV
  use simparam, only: pi
  implicit none

  real(kind=real64), intent(in) :: mu, sigma, x
  real(kind=real64)             :: output

  output = (1/(sigma*dsqrt(2*pi)))*dexp(-(x - mu)**2/(2*sigma**2))
end function gaussian

