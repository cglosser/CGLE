program main
  use ISO_FORTRAN_ENV
  use simParam
  use D2Q9Const
  use cgle
#ifdef VISUALIZATION
  use CGLEvis
#endif
  implicit none

  integer :: time
  complex(kind=real64), allocatable :: f_density(:,:,:,:), feq(:,:,:,:)
  complex(kind=real64), allocatable :: psi(:,:,:), omega(:,:,:)
  complex(kind=real64), allocatable :: polarization(:,:)

  call readSimParam("input.txt")
#ifdef VISUALIZATION 
  call plot_init()
#endif

  allocate(f_density(rDim, cDim, 0:numQ - 1, numSpin))
  allocate(feq(rDim, cDim, 0:numQ - 1, numSpin))

  allocate(psi(rDim, cDim, numSpin))
  allocate(omega(rDim, cDim, numSpin))

  allocate(polarization(rDim, cDim))

  call initRandomF(f_density)
  !call initSpiralF(f_density)

  do time = 1, numTimesteps
    call neumannBC(f_density)

    call computeMacros(f_density, psi, omega)
    call computeFeq(psi, feq)
    call collide(feq, omega, f_density)
    call stream(f_density)

    polarization = spinPolarization(psi(:,:,:))
#ifdef VISUALIZATION
    if(mod(time, 5) .eq. 0) call plot_array(atan2(real(polarization),aimag(polarization)))
#endif

    write(*,*) time, minval(real(polarization)), minval(aimag(polarization))
    write(102,*) real(polarization(rDim/2, cDim/2)), aimag(polarization(rDim/2, cDim/2))
  end do

#ifdef VISUALIZATION
  call plot_close()
#endif

  deallocate(f_density, feq, psi, omega)
end program main
