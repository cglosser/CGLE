module simParam
  use ISO_FORTRAN_ENV
  implicit none

  real(kind=real64), parameter    :: pi = 4*datan(1d0)
  integer              :: rDim
  integer              :: cDim
  integer              :: numTimesteps
  real(kind=real64)    :: boxLength
  real(kind=real64)    :: deltaX, deltaT
  real(kind=real64)    :: tau
  real(kind=real64)    :: t0_coef
  real(kind=real64)    :: latticeVelocity, lambda
  real(kind=real64)    :: beta
  complex(kind=real64) :: a
  complex(kind=real64) :: d

contains

  subroutine readSimParam(fname)
    use ISO_FORTRAN_ENV
    implicit none

    character(len=*), intent(in) :: fname
    character(len=80) :: token
    integer           :: readStat
    real(kind=real64) :: realPart, imagPart
    open(unit=20, file=fname, action="read", iostat=readStat)

    do 
      read(20, *, iostat=readStat) token; backspace(20)
      if(readStat .lt. 0) exit

      select case(trim(token))
        case ("rDim")
          read(20, *, iostat=readStat) token, rDim
        case ("cDim")
          read(20, *, iostat=readStat) token, cDim
        case ("numTimesteps")
          read(20, *, iostat=readStat) token, numTimesteps
        case ("boxLength")
          read(20, *, iostat=readStat) token, boxLength
        case ("deltaX")
          read(20, *, iostat=readStat) token, deltaX
        case ("deltaT")
          read(20, *, iostat=readStat) token, deltaT
        case ("tau")
          read(20, *, iostat=readStat) token, tau
        case ("t0_coef")
          read(20, *, iostat=readStat) token, t0_coef
        case ("beta")
          read(20, *, iostat=readStat) token, beta
        case ("a")
          read(20, *, iostat=readStat) token, realPart, imagPart
          a = dcmplx(realPart, imagPart)
        case ("d")
          read(20, *, iostat=readStat) token, realPart, imagPart
          d = dcmplx(realPart, imagPart)
        case default
          read(20, *) !force I/O to advance
      end select
    end do

    latticeVelocity = deltaX/deltaT
    lambda = 2d0/(deltaT*(2*tau-1))

    close(20)
  end subroutine readSimParam

end module simParam
