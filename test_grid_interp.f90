program test_grid_interp

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


  integer        ::  nin
  real(REAL64)   ::  xgmin, xgmax, xdif
  integer        ::  ng, mxdng
  integer        ::  nd
  integer        ::  nord, nordp1

  real(REAL64)   ::  a, b
  real(REAL64)   ::  rmax
  integer        ::  ifail

!   real(REAL64)   ::  alfa
!   real(REAL64)   ::  fxg, dy
!   integer        ::  jlo, j0

!   integer        ::  isx, ierr
!   real(REAL64)   ::  yp(101)
!   real(REAL64)   ::  ypp(101)
!   real(REAL64)   ::  w(101,3)
!
!   real(REAL64)   ::  fgs(1000)
!   real(REAL64)   ::  xgs(1000)
!   real(REAL64)   ::  ypg(1000)
!   real(REAL64)   ::  yppg(1000)
!
!   real(REAL64)   ::  tmp
!   integer        ::  idum
!   real(REAL64)   ::  small

!   real(REAL64), external   ::  ran2
!   integer        ::  l

! allocatable arrays

  real(REAL64), allocatable   ::  xin(:), fin(:)
  real(REAL64), allocatable   ::  y(:)

  real(REAL64), allocatable   ::  xg(:), fg(:,:)
  real(REAL64), allocatable   ::  dymax(:)


  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

  integer        ::  i, j



  rmax = 20.0_REAL64

  nd = 2
  nord = 7
  nordp1 = nord+1

  nin = 101

  allocate(xin(nin), fin(nin))
  allocate(y(0:nd))

  b = 0.1
  a = rmax / (exp(b*(nin-1))-UM)
  do i = 1,nin
    j = i
!    j = nin - i + 1
    xin(i) = a*(exp(b*(j-1))-UM)
    call fexact(xin(i), y, 0)
    fin(i) = y(0)
  enddo

! swaps two values

!   tmp = xin(30)
!   xin(30) = xin(50)
!   xin(50) = tmp
!   tmp = fin(30)
!   fin(30) = fin(50)
!   fin(50) = tmp

! duplicate point

!  xin(30) = xin(29) + 1.0E-9_REAL64

!  !  adds some noise
!  idum = 12345
!  small = 0.001
!  do i = 1,nin
!    fin(i) = fin(i) + small*ran2(idum)
!  enddo

  do i = 1,nin
    write(200,'(2f16.6)') xin(i),fin(i)
  enddo

! Lagrange interpolation

  xgmin = 0.0
  xgmax = 20.0

  ng = 501
  mxdng = ng

  allocate(xg(mxdng),fg(mxdng,0:nd))
  allocate(dymax(0:nd))

  xdif = (xgmax - xgmin)/(ng-1)
  do i = 1,ng
    j = i
!    j = ng - i + 1
    xg(i) = xgmin + (j-1)*xdif
  enddo

! swaps two values

!   tmp = xg(30)
!   xg(30) = xg(50)
!   xg(50) = tmp

  call grid_interp(nordp1-1, nd, nin, xin, fin, ng, xg, fg, dymax, ifail, mxdng)

  if(ifail /=0) STOP

  write(6,'("  error estimate: ",10e10.3)') (dymax(j),j=0,nd)

  do j = 0,nd
    dymax(j) = ZERO
  enddo

  do i = 1,ng
    call fexact(xg(i), y, nd)
    do j = 0,nd
      write(100+j,'(2f16.6,f22.14)') xg(i),fg(i,j),fg(i,j) - y(j)
      dymax(j) = max(dymax(j),abs(fg(i,j) - y(j)))
    enddo
  enddo

  write(6,'("  observed error: ",10e10.3)') (dymax(j),j=0,nd)

  write(6,*)
  write(6,*) '  The results were written to tapes 200, 100, 101, 102'
  write(6,*) '  200 contains the original x_1, y_i points'
  write(6,*) '  100 contains interpolated x, y, error'
  write(6,*) '  101 and 102 the same for the first and second derivative'
  write(6,*)
  write(6,*) '  try "pl  ''fort.200'', ''fort.100''" in gnuplot'
  write(6,*) '  try "pl  ''fort.101'' u 1:3" in gnuplot'
  write(6,*)

! The following lines can be used to compare with other interpolations

! Spline interpolation

!   l = 0
!   isx = 0
!   call splift (xin,fin,yp,ypp,nin,w,ierr,isx,ZERO,ZERO,ZERO,ZERO)
!
!   if(ierr /= 1) write(6,*) '  splift  ierr = ',ierr
!
!   call splint (xin,fin,ypp,nin,xg,fgs,ypg,yppg,ng,ierr)
!
!   if(ierr /= 1) write(6,*) '  splint  ierr = ',ierr
!
!   do i = 1,ng
!     fexact = xg(i)**l * exp(-alfa*xg(i)*xg(i) / rmax*rmax)
!     write(102,'(2f16.6,f22.14)') xg(i),fgs(i),fgs(i) - fexact
!   enddo

!  do i = 1,ng
!      fexact = xg(i)**l * exp(-alfa*xg(i)*xg(i) / rmax*rmax)
!     if(l == 0) then
!       fexact = -2*(alfa*xg(i) / rmax*rmax)*exp(-alfa*xg(i)*xg(i) / rmax*rmax)
!     else
!       fexact = l * xg(i)**(l-1) * exp(-alfa*xg(i)*xg(i) / rmax*rmax) - &
!              xg(i)**l * 2*(alfa*xg(i) / rmax*rmax)*exp(-alfa*xg(i)*xg(i) / rmax*rmax)
!     endif
!!    fexact = xg(i)**l * exp(-alfa*xg(i)*xg(i) / rmax*rmax)
!     write(112,'(2f16.6,f22.14)') xg(i),ypg(i),ypg(i)-fexact
!   enddo

! numerical recipes subroutine

!   jlo = 0
!   do i = 1,ng
!     call hunt(xin,nin,xg(i),jlo)
!     j0 = min(max(jlo-(nordp1-1)/2,1),nin+1-nordp1)
!     call polint(xin(j0),fin(j0),nordp1,xg(i),fxg,dy)
!     fexact = xg(i)**l * exp(-alfa*xg(i)*xg(i) / rmax*rmax)
!     write(103,'(2f16.6,f22.14)') xg(i),fxg,fxg - fexact
!   enddo

  stop

end program test_grid_interp


!>    indexes a real array by the heapsort method
!>    adapted from http://rosettacode.org
!>    see also W. H. Preuss et al. Numerical Recipes

      subroutine sort(n,a,indx)

!     written 24 June 2013. JLM
!     Modified documentation August 2019.  JLM
!     copyright  J.L.Martins, INESC-MN.

!      version 4.94

      implicit none
      integer, parameter  :: REAL64 = selected_real_kind(12)

!     input

      integer, intent(in)        ::  n                                   !<  length of array
      real(REAL64), intent(in)   ::  a(n)                                !<  array to be indexed

!     output

      integer, intent(out)       ::  indx(n)                             !<  index of array a

!     local variables

      integer    ::  iroot,ichild,istart,ibot,indxt,ic

      if(n < 1) return

      do ichild=1,n
        indx(ichild) = ichild
      enddo

      if(n == 1) return

!     hiring phase

      ibot = n

      do istart = n/2,1,-1
        indxt = indx(istart)
        iroot = istart

!       long enough siftdown loop does not exceed ~log(n)/log(2)

        do ic = 1,n+5
          ichild = 2*iroot
          if(ichild <= ibot) then
            if(ichild < ibot) then
              if(a(indx(ichild)) < a(indx(ichild+1)))                    &
     &              ichild = ichild + 1
            endif

            if(a(indxt) < a(indx(ichild))) then
              indx(iroot) = indx(ichild)
              iroot = ichild
            else

              exit

            endif

          else

            exit

          endif

        enddo

        indx(iroot) = indxt

      enddo

!     retirement and promotion phase

      istart = 1

      do ibot = n-1,1,-1
        indxt = indx(ibot+1)
        indx(ibot+1) = indx(1)
        if(ibot == 1) then

          exit

        endif

        iroot = istart

!       long enough siftdown loop does not exceed ~log(n)/log(2)  (repeated...)

        do ic = 1,n+5
          ichild = 2*iroot
          if(ichild <= ibot) then
            if(ichild < ibot) then
              if(a(indx(ichild)) < a(indx(ichild+1)))                    &
     &              ichild = ichild + 1
            endif

            if(a(indxt) < a(indx(ichild))) then
              indx(iroot) = indx(ichild)
              iroot = ichild
            else

              exit

            endif

          else

            exit

          endif

        enddo

        indx(iroot) = indxt

      enddo

      indx(1) = indxt

      return

      end subroutine sort


subroutine fexact(x, y, nd)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, intent(in)        ::  nd                        !<  order of derivative
  real(REAL64), intent(in)   ::  x                         !<  abcissa
  real(REAL64), intent(out)  ::  y(0:nd)                   !<  f(x),  f'(x), f''(x)

  integer, parameter         ::  L = 0                     !<  hard coded positive integer
  real(REAL64), parameter    ::  ALFA = 0.4D0

  real(real64)               ::  eax, deaxdx, d2eaxdx2

! constants

  real(REAL64), parameter    ::  ZERO = 1.0_REAL64

! counters

  integer        ::  j

  eax = exp(-ALFA*x*x)
  deaxdx = -2*ALFA*x*eax
  d2eaxdx2 = -2*ALFA*eax -2*ALFA*x*deaxdx

  do j = 0,nd
    y(j) = ZERO
  enddo

  if(L == 0) then
    y(0) = eax
    if(nd > 0) then
      y(1) = deaxdx
      if(nd > 1) then
        y(2) = d2eaxdx2
      endif
    endif
  elseif(L == 1) then
    y(0) = x * eax
    if(nd > 0) then
      y(1) = eax + x*deaxdx
      if(nd > 1) then
        y(2) = 2*deaxdx + x*d2eaxdx2
      endif
    endif
  else
    y(0) = x**L * eax
    if(nd > 0) then
      y(1) = L*x**(L-1) * eax + x**L * deaxdx
      if(nd > 1) then
        y(2) = L*(L-1)*x**(L-2) * eax + 2*L*x**(L-1) * deaxdx + x**L * d2eaxdx2
      endif
    endif
  endif

  return

end subroutine fexact


