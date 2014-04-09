module D2Q9Const
  double precision, parameter :: t(0:8) = &
    (/ 4.0d0/9.0d0,  1.0d0/9.0d0,  1.0d0/9.0d0,  &
       1.0d0/9.0d0,  1.0d0/9.0d0,  1.0d0/36.0d0, &
       1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0 /)
  integer,parameter::vectors(0:1,0:8) = &
    reshape((/0,0,1,0,0,1,-1,0,0,-1,1,1,-1,1,-1,-1,1,-1/), shape(vectors))
  integer, parameter :: reverse(0:8) = (/0,3,4,1,2,7,8,5,6/)
end module D2Q9Const

module simParam
  integer, parameter          :: xDim = 5
  integer, parameter          :: yDim = 5
  integer, parameter          :: tMax = 1000
  double precision, parameter :: omega = 1.0d0
end module simParam

program cgle
  use simParam
  use D2Q9Const
  implicit none

  double precision, allocatable, dimension(:,:,:) :: f, feq, u
  double precision, allocatable, dimension(:,:)   :: rho, uSqr

  allocate(f(0:8, ydim, xdim), fEq(0:8, ydim, xdim), u(0:1, ydim, xdim))
  allocate(rho(ydim, xdim), uSqr(ydim, xdim))

  f(:,1,1) = (/1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0, 9.0d0 /)
  call computeMacros(f, rho, uSqr, u)

  write(*,"(9E16.8)") f(:,1,1)
  write(*,"(2E16.8)") u(:,1,1)
  write(*,*) uSqr(1,1)

end program cgle

subroutine computeMacros(f, rho, uSqr, u)
  use D2Q9Const, only: vectors
  use simParam,  only: xDim, yDim
  implicit none

  double precision, intent(in)  :: f(0:8, ydim, xdim)
  double precision, intent(out) :: u(0:1, ydim, xdim), rho(ydim, xdim), usqr(ydim, xdim)
  integer :: x, y
  do x = 1, xdim
    do y = 1, ydim
      rho(y, x) = sum(f(:, y, x))
      u(0,y,x)  = sum(f(:,y,x)*vectors(0,:))
      u(1,y,x)  = sum(f(:,y,x)*vectors(1,:))
      uSqr(y,x) = dot_product(u(:,y,x),u(:,y,x))
    end do
  end do
end subroutine computeMacros
