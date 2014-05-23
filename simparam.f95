module simParam
use ISO_FORTRAN_ENV
! Make sure cDim is even for periodic D2Q7!
integer, parameter              :: rDim   = 100
integer, parameter              :: cDim   = 100
integer, parameter              :: tMax   = 10000
real(kind=real64), parameter    :: normalization = 1d0
real(kind=real64), parameter    :: boxLength = 10d0
real(kind=real64), parameter    :: deltaX = 0.1d0, deltaT = 0.05d0 ! dT = knudsen #
real(kind=real64), parameter    :: latticeVelocity = deltaX/deltaT   
real(kind=real64), parameter    :: tau    = 0.55d0
real(kind=real64), parameter    :: t0_coef= 0.3d0
complex(kind=real64), parameter :: lambda = 2d0/(deltaT*(2*tau - 1))
complex(kind=real64), parameter :: beta   = 2.0d-3 ! beta = diffusion coef. 
complex(kind=real64), parameter :: a = dcmplx(0.100d0, 0.00d0)
complex(kind=real64), parameter :: d = dcmplx(0.025d0, 0.03d0)
end module simParam
