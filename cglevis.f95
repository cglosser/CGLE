module CGLEvis
  use ISO_FORTRAN_ENV
  use plplot
  implicit none
  
  private

  public plot_init, plot_close, plot_array

contains

  subroutine plot_init()
    use simParam, only: rDim, cDim
    implicit none
    call plsdev("xcairo")
    call plinit()

    call plenv(0d0, rdim + 1d0, 0d0, cdim + 1d0, 0, 0)
  end subroutine plot_init

  subroutine plot_close()
    implicit none
    call plspause(.false.)
    call plend()
  end subroutine plot_close

  subroutine plot_array(f)
    use simParam, only: rDim, cDim
    implicit none
    real(kind=real64), intent(in) :: f(rDim, cDim)

    call plimage(f, 1._plflt, 1._plflt*cDim, 1._plflt, 1._plflt*rDim, &
      -10._plflt, 10._plflt, 1._plflt, 1._plflt*cDim, 1._plflt, 1._plflt*rDim)
    call plflush()
  end subroutine plot_array

end module CGLEvis
