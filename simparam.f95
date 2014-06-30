module simParam
  use ISO_FORTRAN_ENV
  implicit none

  real(kind=real64), parameter    :: pi = 4*datan(1d0)
  integer              :: rDim = -1, cDim = -1, numSpin = -1 !define these three first in the input file
  integer              :: numTimesteps
  real(kind=real64)    :: boxLength
  real(kind=real64)    :: deltaX, deltaT
  real(kind=real64)    :: tau
  real(kind=real64)    :: t0_coef
  real(kind=real64)    :: latticeVelocity, lambda
  real(kind=real64)    :: beta
  complex(kind=real64) :: potential 
  complex(kind=real64) :: nonlinear
  complex(kind=real64) :: coupling 

  ! Named indexing constants
  integer, parameter :: spin_up = 1, spin_down = 2

contains

  subroutine readSimParam(fname)
    use ISO_FORTRAN_ENV
    implicit none

    character(len=*), intent(in) :: fname
    character(len=80)  :: token, date_style
    character(len=256) :: error_str
    integer            :: readStat, vals(8)
    real(kind=real64)  :: realPart, imagPart

    open(unit=20, file=fname, action="read", iostat=readStat)
    open(unit=30, file="simulation.log")

    call date_and_time(VALUES=vals)

    do 
      read(20, *, iostat=readStat) token; backspace(20)
      if(readStat .lt. 0) exit

      select case(trim(token))
        case ("rDim")
          read(20, *, iostat=readStat) token, rDim

        case ("cDim")
          read(20, *, iostat=readStat) token, cDim

        case ("numSpin")
          read(20, *, iostat=readStat) token, numSpin

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

        case ("potential")
          read(20, *, iostat=readStat) token, realPart, imagPart
          potential = dcmplx(realPart, imagPart)

        case ("nonlinear")
          read(20, *, iostat=readStat) token, realPart, imagPart
          nonlinear = dcmplx(realPart, imagPart)

        case ("coupling")
          read(20, *, iostat=readStat) token, realPart, imagPart
          coupling = dcmplx(realPart, imagPart)
          if(numSpin .eq. -1) then
            error_str = "WARNING: multi-component coupling constant &
              &defined before the number of"//NEW_LINE('A')//"components; &
              &may produce undesired results."
            write(ERROR_UNIT, *) error_str
            write(30, *) error_str
          else if(numSpin .eq. 1 .and. coupling .ne. dcmplx(0d0, 0d0)) then
            error_str = "WARNING: multi-component coupling constant &
              &nonzero for a single-component"//NEW_LINE('A')//" simulation. &
              &I will treat it as an addition to the potential."
            write(ERROR_UNIT, *) error_str
            write(30, *) error_str
          end if

        case default
          read(20, *) !force I/O to advance
      end select
    end do

    deltaX = boxLength/max(rDim, cDim)
    latticeVelocity = deltaX/deltaT
    lambda = 2d0/(deltaT*(2*tau-1))

    date_style = "(A,I4,2I2.2,I3,A,I2.2)"
    write(30,date_style) "Simulation started on ", vals(1:3),vals(5),":",vals(6)
    write(30,"(A,I2,A)") "Running with", numSpin, " components."
    write(30,*) "        delta x:", deltaX
    write(30,*) "latticeVelocity:", latticeVelocity
    write(30,*) "         lambda:", lambda
    write(30,*) "     total time:", deltaT*numTimesteps

    close(20); close(30)
  end subroutine readSimParam

end module simParam
