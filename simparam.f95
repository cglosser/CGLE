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

  ! Named indexing constants
  integer, parameter :: numSpin = 2
  integer, parameter :: spin_up = 1, spin_down = 2

contains

  subroutine readSimParam(fname)
    use ISO_FORTRAN_ENV
    implicit none

    character(len=*), intent(in) :: fname
    character(len=80) :: token, style
    integer           :: readStat, vals(8)
    real(kind=real64) :: realPart, imagPart
    open(unit=20, file=fname, action="read", iostat=readStat)
    open(unit=30, file="simulation.log")

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

    deltaX = boxLength/max(rDim, cDim)
    latticeVelocity = deltaX/deltaT
    lambda = 2d0/(deltaT*(2*tau-1))

    call date_and_time(VALUES=vals)

    style = "(A,I4,2I2.2,I3,A,I2.2)"
    write(30,style) "Simulation started on ", vals(1:3),vals(5),":",vals(6)
    write(30,"(A,I2)") "Running in mode", 1
    write(30,*) "        delta x:", deltaX
    write(30,*) "latticeVelocity:", latticeVelocity
    write(30,*) "         lambda:", lambda

    close(20); close(30)
  end subroutine readSimParam

end module simParam
