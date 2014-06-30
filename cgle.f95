program cgle
  use ISO_FORTRAN_ENV
  use simParam
  use D2Q9Const
  use CGLEvis
  implicit none

  integer :: time
  complex(kind=real64), allocatable :: f_density(:, :, :, :), feq(:, :, :, :)
  complex(kind=real64), allocatable :: psi(:, :, :), omega(:, :, :)

  call readSimParam("input.txt")
  call plot_init()

  allocate(f_density(rDim, cDim, 0:numQ - 1, numSpin))
  allocate(feq(rDim, cDim, 0:numQ - 1, numSpin))

  allocate(psi(rDim, cDim, numSpin))
  allocate(omega(rDim, cDim, numSpin))

  !call initRandomF(f_density)
  call initSpiralF(f_density)
  do time = 1, numTimesteps
    call neumannBC(f_density)

    call computeMacros(f_density, psi, omega)
    call computeFeq(psi, feq)
    call collide(feq, omega, f_density)
    call stream(f_density)

    if(mod(time, 5) .eq. 0) call plot_array(real(psi(:,:,spin_up)))

    write(*,*) time
  end do

  call plot_close()

  deallocate(f_density, feq, psi, omega)
end program cgle

subroutine computeMacros(f_density, psi, omega)
  use ISO_FORTRAN_ENV
  use D2Q9Const, only: vectors, numQ
  use simParam,  only: rDim, cDim, deltaT, a, d, numSpin
  implicit none

  complex(kind=real64), intent(in)    :: f_density(rDim, cDim, 0:numQ - 1, numSpin)
  complex(kind=real64), intent(out)   :: psi(rDim, cDim, numSpin), &
                                         omega(rDim, cDim, numSpin)
  complex(kind=real64) :: hamiltonian(rDim, cDim, numSpin)
  
  psi = sum(f_density, dim = 3)
  hamiltonian = (a - d*abs(psi)**2)*psi + 0.1*cshift(psi, 1, dim = 3) !<-- check this!
  omega = deltaT*hamiltonian/numQ
end subroutine computeMacros

subroutine computeFeq(psi, feq)
  use ISO_FORTRAN_ENV
  use simParam,  only: cDim, rDim, lambda, beta, latticeVelocity, numSpin
  use D2Q9Const, only: numQ, latticeDim
  implicit none

  complex(kind=real64), intent(in)  :: psi(rDim, cDim, numSpin)
  complex(kind=real64), intent(out) :: feq(rDim, cDim, 0:numQ - 1, numSpin)
  complex(kind=real64) :: tmp(rDim, cDim, numSpin)
  integer :: dir

  feq(:,:,0,:) =  (1 - 3*lambda*beta*latticeDim/(4*latticeVelocity**2))*psi(:,:,:)
  tmp(:,:,:)   = lambda*beta*latticeDim/((numQ - 1)*latticeVelocity**2)*psi(:,:,:)

  do dir = 1, 4
    feq(:,:,dir,:)     = tmp(:,:,:)
    feq(:,:,dir + 4,:) = tmp(:,:,:)/2
  end do
end subroutine computeFeq

subroutine stream(f_density) 
  use ISO_FORTRAN_ENV
  use simParam,  only: cDim, rDim, numSpin
  use D2Q9Const, only: numQ
  implicit none

  complex(kind=real64), intent(inout) :: f_density(rDim, cDim, 0:numQ - 1, numSpin)

  !!--stream along 1------------------------------
  f_density(:,:,1,:) = cshift(f_density(:,:,1,:), -1, dim=1)

  !!--stream along 2------------------------------
  f_density(:,:,2,:) = cshift(f_density(:,:,2,:),  1, dim=2)

  !!--stream along 3------------------------------
  f_density(:,:,3,:) = cshift(f_density(:,:,3,:),  1, dim=1)

  !!--stream along 4------------------------------
  f_density(:,:,4,:) = cshift(f_density(:,:,4,:), -1, dim=2)

  !!--stream along 5------------------------------
  f_density(:,:,5,:) = cshift(cshift(f_density(:,:,5,:), -1, dim=1), 1, dim=2)

  !!--stream along 6------------------------------
  f_density(:,:,6,:) = cshift(cshift(f_density(:,:,6,:),  1, dim=1), 1, dim=2)

  !!--stream along 7------------------------------
  f_density(:,:,7,:) = cshift(cshift(f_density(:,:,7,:),  1, dim=1), -1, dim=2)

  !!--stream along 7------------------------------
  f_density(:,:,8,:) = cshift(cshift(f_density(:,:,8,:), -1, dim=1), -1, dim=2)
end subroutine stream

subroutine collide(fEq, omega, f_density)
  use ISO_FORTRAN_ENV
  use D2Q9Const, only: numQ
  use simParam,  only: rDim, cDim, tau, numSpin
  implicit none

  complex(kind=real64), intent(inout) :: f_density(rDim, cDim, 0:numQ - 1, numSpin)
  complex(kind=real64), intent(in)    :: fEq(rDim, cDim, 0:numQ - 1, numSpin), &
                                         omega(rDim, cDim, numSpin)

  integer :: dir
  do dir = 0, numQ - 1
    f_density(:,:,dir,:) = f_density(:,:,dir,:) - &
      (f_density(:,:,dir,:) - fEq(:,:,dir,:))/tau + omega(:,:,:)
  end do
end subroutine collide

subroutine neumannBC(f_density)
  use ISO_FORTRAN_ENV
  use D2Q9Const, only: numQ
  use simParam,  only: rDim, cDim, numSpin
  implicit none

  complex(kind=real64), intent(inout) :: f_density(rDim, cDim, 0:numQ - 1, numSpin)

  f_density(2:rDim - 1, 1, :, :)    = f_density(2:rDim - 1, 2, :, :)
  f_density(2:rDim - 1, cDim, :, :) = f_density(2:rDim - 1, cDim - 1, :, :)
  f_density(1, :, :, :)    = f_density(2, :, :, :)
  f_density(rdim, :, :, :) = f_density(rdim - 1, :, :, :)
end subroutine neumannBC

subroutine initSpiralF(f_density)
  use ISO_FORTRAN_ENV
  use D2Q9Const, only: numQ
  use simParam,  only: rDim, cDim, t0_coef, deltaX, numSpin, spin_up, spin_down
  implicit none

  complex(kind=real64), intent(out) :: f_density(rDim, cDim, 0:numQ - 1, numSpin)
  integer :: rIdx, cIdx
  real(kind=real64) :: x, y
  do cIdx = 1, cDim
    x = (cIdx - 1 - cDim/2)*deltaX
    do rIdx = 1, rDim
      y = (rIdx - 1 - rDim/2)*deltaX
      f_density(rIdx, cIdx, :, spin_up)   = t0_coef*dcmplx(x, y)/numQ
      f_density(rIdx, cIdx, :, spin_down) = t0_coef*dcmplx(y, x)/numQ
    end do
  end do
end subroutine initSpiralF

subroutine initRandomF(f_density)
  use ISO_FORTRAN_ENV
  use D2Q9Const, only: numQ
  use simParam,  only: rDim, cDim, t0_coef, numSpin, spin_up, spin_down
  implicit none
  
  complex(kind=real64), intent(out) :: f_density(rDim, cDim, 0:numQ - 1, numSpin)
  integer :: rIdx, cIdx
  real(kind=real64) :: rands(2)
  do cIdx = 1, cDim
    do rIdx = 1, rDim
      call random_number(rands)
      rands = rands - 0.5d0
      f_density(rIdx, cIdx, :, spin_up) = t0_coef*dcmplx(rands(1), rands(2))/numQ

      call random_number(rands)
      rands = rands - 0.5d0
      f_density(rIdx, cIdx, :, spin_down) = t0_coef*dcmplx(rands(1), rands(2))/numQ
    end do
  end do
end subroutine initRandomF
