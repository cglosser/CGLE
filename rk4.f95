program rk4
  use ISO_FORTRAN_ENV
  use simParam
  use CGLEvis
  implicit none

  integer :: time
  complex(kind=real128), allocatable :: psi(:,:,:), hamiltonian(:,:,:)

  call readSimParam("input.txt")

  allocate(psi(rDim, cDim, 2), hamiltonian(rDim, cDim, 2))

  call plot_init()

  !call initSpiralF(psi)
  psi = 0d0
  do time = 1, numTimesteps
    call integration(psi)

    if (mod(time, 5) .eq. 0) call plot_array(real(psi(:,:,exciton), kind=real64))

    write(*,*) time, time*deltaT, sum(abs(psi)**2)
  end do

  call plot_close()

end program rk4

subroutine calcMacros(psi, hamiltonian)
  use ISO_FORTRAN_ENV
  use simParam, only: rDim, cDim, deltaX, beta, a, d, omegaRabi, exciton, cavity
  implicit none

  complex(kind=real128), intent(in)  :: psi(rDim, cDim, 2)
  complex(kind=real128), intent(out) :: hamiltonian(rDim, cDim, 2)

  hamiltonian(:,:,exciton) = beta(exciton)*laplacian(psi(:,:,exciton)) &
    + (a(exciton) - d(exciton)*abs(psi(:,:,exciton))**2)*psi(:,:,exciton) &
    + omegaRabi/2*psi(:,:,cavity)

  hamiltonian(:,:,cavity) = beta(cavity)*laplacian(psi(:,:,cavity)) &
    + (a(cavity) - d(cavity)*abs(psi(:,:,cavity))**2)*psi(:,:,cavity) &
    + omegaRabi/2*psi(:,:,exciton)

  contains
    
    function laplacian(psi) result(del_sq)
      use ISO_FORTRAN_ENV
      use simParam, only: rDim, cDim, deltaX
      implicit none

      complex(kind=real128), intent(in) :: psi(rDim, cDim)
      complex(kind=real128) :: del_sq(rDim, cDim)

      del_sq = (-2*psi + cshift(psi,-1,1) + cshift(psi,1,1) &
                -2*psi + cshift(psi,-1,2) + cshift(psi,1,2))/deltaX
    end function laplacian

end subroutine calcMacros

subroutine integration(psi)
  use ISO_FORTRAN_ENV
  use simParam, only: rDim, cDim, deltaT
  implicit none

  complex(kind=real128), intent(inout) :: psi(rDim, cDim, 2)
  complex(kind=real128) :: psiPrime(rDim, cDim, 2)
  complex(kind=real128) :: k1(rDim, cDim, 2), k2(rDim, cDim, 2)
  complex(kind=real128) :: k3(rDim, cDim, 2), k4(rDim, cDim, 2), psi_rsq

  call NeumannBCs(psi)
  call calcMacros(psi, k1)

  psiPrime = psi + deltaT/2*k1
  call NeumannBCs(psiPrime)
  call calcMacros(psiPrime, k2)

  psiPrime = psi + deltaT/2*k2
  call NeumannBCs(psiPrime)
  call calcMacros(psiPrime, k3)

  psiPrime = psi + deltaT/2*k3
  call NeumannBCs(psiPrime)
  call calcMacros(psiPrime, k4)
  psi = psi + deltaT*(k1 + 2*k2 + 2*k3 + k4)/6
end subroutine integration

subroutine NeumannBCs(psi)
  use ISO_FORTRAN_ENV
  use simParam, only: rDim, cDim
  implicit none
  complex(kind=real128), intent(inout) :: psi(rDim, cDim)

  psi(2:rDim - 1, 1)    = psi(2:rDim - 1, 2)
  psi(2:rDim - 1, cDim) = psi(2:rDim - 1, cDim - 1)
  psi(1, :)    = psi(2, :)
  psi(rdim, :) = psi(rdim - 1, :)
end subroutine NeumannBCs

subroutine initSpiralF(psi)
  use ISO_FORTRAN_ENV
  use simParam,  only: rDim, cDim, t0_coef, deltaX, boxLength
  implicit none

  complex(kind=real128), intent(out) :: psi(rDim, cDim)
  integer :: rIdx, cIdx
  real(kind=real128) :: x, y
  do cIdx = 1, cDim
    x = (cIdx - 1 - cDim/2)*deltaX
    do rIdx = 1, rDim
      y = (rIdx - 1 - rDim/2)*deltaX
      psi(rIdx, cIdx) = t0_coef*dcmplx(x, y)
    end do
  end do

  psi = psi/sqrt(sum(abs(psi)**2)/size(psi))

  !call NeumannBCs(psi)
end subroutine initSpiralF

subroutine initRandomF(psi)
  use ISO_FORTRAN_ENV
  use simParam,  only: rDim, cDim, t0_coef, deltaX, boxLength
  implicit none
  
  complex(kind=real128), intent(out) :: psi(rDim, cDim)
  integer :: rIdx, cIdx
  real(kind=real128) :: rands(2)
  do cIdx = 1, cDim
    do rIdx = 1, rDim
      call random_number(rands)
      rands = rands - 0.5d0
      psi(rIdx, cIdx) = t0_coef*dcmplx(rands(1), rands(2))
    end do
  end do
end subroutine initRandomF
