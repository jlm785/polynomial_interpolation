program test_poly

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! variables

  integer                      ::  n                 !  order of interpolation
  integer                      ::  nd                !  order of derivative
  real(REAL64), allocatable    ::  xin(:), yin(:)    !  grid points
  real(REAL64), allocatable    ::  yex(:)            !  exact resultr
  real(REAL64)                 ::  x                 !  point of desired interpolation
  real(REAL64), allocatable    ::  y(:)              !  interpolated function and its derivatives
  real(REAL64), allocatable    ::  dy(:)             !  estimated error

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64

! counters

  integer     ::  j, k

  n = 5
  nd = n

  allocate(xin(0:n),yin(0:n))
  allocate(y(0:nd),dy(0:nd))
  allocate(yex(0:nd))

  do j = 0,n
    xin(j) = 10*ONE+j*ONE
  enddo

  do j = 0,n
    call fexact(xin(j), yex, 0)
    yin(j) = yex(0)
  enddo



  x = 12.3_REAL64
  call fexact(x, yex, nd)

  do j = 0,n
    xin(j) = xin(j) - x
  enddo

  call poly_interp(y, dy, xin, yin, n, nd)

  write(6,*)
  write(6,*) '  interpolated      exact        error      last update'
  write(6,*)
  do k = 0,nd
    write(6,'(4f14.9,f24.16,f14.8)') y(k), yex(k), y(k)-yex(k),dy(k)
  enddo

  stop

end program test_poly

subroutine fexact(x, y, nd)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, intent(in)        ::  nd                        !<  order of derivative
  real(REAL64), intent(in)   ::  x                         !<  abcissa
  real(REAL64), intent(out)  ::  y(0:nd)                   !<  f^j(x), j-th derivative
! variables

  real(REAL64)   ::  factor

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64

! counters

  integer     ::  k

  y(0) = sqrt(x)

  factor = ONE
  if(nd > 0) then
    do k = 1,nd
      factor = factor*(3*ONE-2*k)/(2*ONE)
      y(k) = factor*x**(ONE/2 - k*ONE)
    enddo
  endif

  return

end subroutine fexact

