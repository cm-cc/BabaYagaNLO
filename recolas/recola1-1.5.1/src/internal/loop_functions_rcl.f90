!#####################################################################
!!
!!  File  loop_functions_rcl.f90
!!  is part of RECOLA (REcursive Computation of One Loop Amplitudes)
!!
!!  Copyright (C) 2015-2025   Stefano Actis, Ansgar Denner,
!!                            Lars Hofer, Jean-Nicolas Lang,
!!                            Andreas Scharf, Sandro Uccirati
!!
!!  RECOLA is licenced under the GNU GPL version 3,
!!         see COPYING for details.
!!
!#####################################################################

  module loop_functions_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl
  use collier

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Auxiliary subroutines and functions
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine settype (xp,x0,x1,type)

  real    (dp), intent (in)  :: xp
  complex (dp), intent (in)  :: x0,x1
  integer,      intent (out) :: type

  logical  :: xplog,x0log,x1log,xdlog,p0log,p1log,p0big,p1big
  real(dp) :: mu

  mu = max(abs(xp),abs(x0),abs(x1))

  if (mu.eq.0d0) then
    xplog = .true. ! p  = 0
    x0log = .true. ! m0 = 0
    x1log = .true. ! m1 = 0
    xdlog = .true. ! m0 = m1
    p0log = .true. ! p  = m0
    p1log = .true. ! p  = m1
  else
    xplog = abs(xp)/mu    .lt. zerocut ! p  = 0
    x0log = abs(x0)/mu    .lt. zerocut ! m0 = 0
    x1log = abs(x1)/mu    .lt. zerocut ! m1 = 0
    xdlog = abs(x0-x1)/mu .lt. zerocut ! m0 = m1
    p0log = abs(xp-x0)/mu .lt. zerocut ! p  = m0
    p1log = abs(xp-x1)/mu .lt. zerocut ! p  = m1
  endif

  p0big = .not. p0log ! p ne m0
  p1big = .not. p1log ! p ne m1

  ! define cases for two-point funtions

  if (xplog.and.x0log.and.x1log) then
        type = 1  ! p=m0=m1=0
  else
    if (xplog) then
      if     (x0log) then
        type = 2  ! p=m0=0
      elseif (x1log) then
        type = 3  ! p=m1=0
      elseif (xdlog) then
        type = 4  ! p=0,m0=m1
      else
        type = 5  ! p=0, m0 ne m1 (ne 0)
      endif
    else
      if     (x0log.and.x1log) then
        type = 6  ! p only, m0=m1=0
      elseif (x1log.and.p0log) then ! IR derivative
        type = 7  ! p=m0,m1=0
      elseif (x0log.and.p1log) then ! IR derivative
        type = 8  ! p=m1,m0=0
      elseif (x0log.and.p1big) then
        type = 9  ! m0=0, p ne m1 (ne 0)
      elseif (x1log.and.p0big) then
        type = 10 ! m1=0, p ne m0 (ne 0)
      else
        type = 11 ! general case
      endif
    endif
  endif

  end subroutine settype

!---------------------------------------------------------------------

  function cln (z,isig) result (clog)

  ! branch cut for logarithm along the real negative axis

  complex (dp), intent (in) :: z
  real    (dp), intent (in) :: isig
  complex (dp)              :: clog

  logical :: zisreal,zisnegative

  if (abs(aimag(z)).lt.zerocut) then
    zisreal = .true.
  else
    zisreal = .false.
  endif

  if (real(z).le.0d0) then
    zisnegative = .true.
  else
    zisnegative = .false.
  endif

  if (zisreal.and.zisnegative) then
    clog = log(-z) + cId0 * cmplx(sign(pi,isig),kind=dp)
  else
    clog = log( z)
  endif

  end function cln

!---------------------------------------------------------------------

  function fndd (n,x,iep) result (funct)

  ! eq (4.11) of Denner-Dittmaier, 0509141

  integer,      intent (in) :: n
  complex (dp), intent (in) :: x
  real    (dp), intent (in) :: iep
  complex (dp)              :: funct

  complex (dp) :: xm1
  integer      :: j

  xm1 = x-c1d0

  ! follow recipe by DD and choose |x|=10 for switch

  if (abs(x).lt.10d0) then
    if (abs(xm1).lt.zerocut) then
      funct = c0d0
    else
      funct = (c1d0-x**(n+1))*(cln(xm1,iep)-cln(x,iep))
    endif
    do j= 0, n
      funct = funct-x**(n-j)/cmplx((j+1),kind=dp)
    enddo
  elseif (abs(x).ge.10d0) then
    funct = cln(c1d0-c1d0/x,iep)
    do j= n+1, n+infty
      funct = funct+x**(n-j)/cmplx((j+1),kind=dp)
    enddo
  endif

  end function fndd

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! One-loop functions
!
! UV divergencies are dimensionally regulated in 4-2*ep
! dimensions; the divergency factor is defined as:
!
!   Delta_UV = 1/ep - EulerConstant + ln(4*pi)
!
! If IR divergencies are dimensionally regulated in 4-2*ep
! dimensions, divergency factor is defined as:
!
!   Delta_IR = 1/ep - EulerConstant + ln(4*pi)
!
! Correspondingly, two arbitrary squared units of mass are used,
! muUV2 (UV) and muIR2 (IR), confined to the finite part of a function
! whenever an UV/IR divergence appears.
!
! B0 function with no scales is:
!
!   B0(0,0,0) = Delta_UV - Delta_IR + ln(muUV2/muIR2)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function A0 (m2,Duv,muUV2)

  complex (dp)              :: A0
  real    (dp), intent (in) :: Duv,muUV2
  complex (dp), intent (in) :: m2

  if (abs(m2).eq.0d0) then
    A0 = c0d0
  else
    if (collier_ct) then
      call A0_cll(A0,m2)
    else
      A0 = m2 * ( Duv + 1d0 - log(m2/muUV2) )
    endif
  endif

  end function A0

!---------------------------------------------------------------------

  function B0 (p2,m02,m12,Duv,muUV2,Dir,muIR2)

  complex (dp)                        :: B0
  real    (dp),           intent (in) :: p2,Duv,muUV2
  complex (dp),           intent (in) :: m02,m12
  real    (dp), optional, intent (in) :: Dir,muIR2

  integer      :: type
  complex (dp) :: b,rt,xp,xm

  if (collier_ct) then
    call B0_cll(B0,cmplx(p2,kind=dp),m02,m12)
  else
    call settype (p2,m02,m12,type)
    select case (type)
    case (1)  ! p=m0=m1=0
      if (present(Dir).and.present(muIR2)) then
        B0 = Duv - Dir + log(muUV2/muIR2)
      else
        if (warnings(321).le.warning_limit) then
          warnings(321) = warnings(321) + 1
          call openOutput
          write(nx,*)
          write(nx,*) 'CODE ERROR 321 (loop_functions_rcl.f90):'
          write(nx,*) 'B0 with 0 scales requires values for Dir and muIR2'
          write(nx,*)
          call toomanywarnings(321)
        endif
        call istop (ifail,2)
      endif
    case (2)  ! p=m0=0
      B0 = Duv + 1d0 - log(m12/muUV2)
    case (3)  ! p=m1=0
      B0 = Duv + 1d0 - log(m02/muUV2)
    case (4)  ! p=0,m0=m1
      B0 = Duv - log(m02/muUV2)
    case (5)  ! p=0,m0,m1
      B0 = Duv + 1d0 &
               - ( m02*log(m02/muUV2) - m12*log(m12/muUV2) )/(m02-m12)
    case (6)  ! m0=m1=0
      B0 = Duv + 2d0 - cln(-p2/muUV2*c1d0,-1d0)
    case (7)  ! p=m0,m1=0
      B0 = Duv + 2d0 - log(m02/muUV2)
    case (8)  ! p=m1,m0=0
      B0 = Duv + 2d0 - log(m12/muUV2)
    case (9)  ! m0=0
      B0 = Duv + 2d0 - log(m12/muUV2) &
               + (m12/p2-1d0)*cln(1d0-p2/m12,-1d0)
    case (10) ! m1=0
      B0 = Duv + 2d0 - log(m02/muUV2) &
               + (m02/p2-1d0)*cln(1d0-p2/m02,-1d0)
    case (11) ! the general case
      b   = m12 - m02 - p2
      rt  = sqrt(b**2-4*p2*m02)
      xp = (-b+rt)/(2*p2)
      xm = (-b-rt)/(2*p2)
      B0 = Duv + log(muUV2/m02) - fndd(0,xp,1d0) - fndd(0,xm,-1d0)
    case default
      if (warnings(322).le.warning_limit) then
        warnings(322) = warnings(322) + 1
        call openOutput
        write(nx,*)
        write(nx,*) "CODE ERROR 322 (loop_functions_rcl.f90): B0 case doesn't exist"
        write(nx,*)
        call toomanywarnings(322)
      endif
      call istop (ifail,2)
    end select
  endif

  end function B0

!---------------------------------------------------------------------

  function DB0 (p2,m02,m12,Dir,muIR2)
  ! no UV singularities

  complex (dp)                        :: DB0
  real    (dp),           intent (in) :: p2
  complex (dp),           intent (in) :: m02,m12
  real    (dp), optional, intent (in) :: Dir,muIR2

  integer      :: type
  complex (dp) :: b,rt,xp,xm

  if (collier_ct) then
    call DB0_cll(DB0,cmplx(p2,kind=dp),m02,m12)
  else
    call settype (p2,m02,m12,type)
    select case (type)
    case (1)  ! p=m0=m1=0
      DB0 = c0d0
    case (2)  ! p=m0=0
      DB0 = 1/2d0/m12
    case (3)  ! p=m1=0
      DB0 = 1/2d0/m02
    case (4)  ! p=0,m0=m1
      DB0 = 1/6d0/m02
    case (5)  ! p=0,m0,m1
      DB0 = 1/2d0/(m02-m12)**3 *                         &
            ( m02**2 - m12**2 - 2*m02*m12*log(m02/m12) )
    case (6)  ! m0=m1=0
      DB0 = - 1/p2
    case (7)  ! p=m0,m1=0 is IR divergent
      if (present(Dir).and.present(muIR2)) then
        DB0 = - 1/2d0/m02*( Dir + 2 - log(m02/muIR2) )
      else
        if (warnings(323).le.warning_limit) then
          warnings(323) = warnings(323) + 1
          call openOutput
          write(nx,*)
          write(nx,*) 'CODE ERROR 323 (loop_functions_rcl.f90):'
          write(nx,*) 'IR divergent DB0 requires values for Dir and muIR2'
          write(nx,*)
          call toomanywarnings(323)
        endif
        call istop (ifail,2)
      endif
    case (8)  ! p=m1,m0=0 is IR divergent
      if (present(Dir).and.present(muIR2)) then
        DB0 = - 1/2d0/m12*( Dir + 2 - log(m12/muIR2) )
      else
        if (warnings(324).le.warning_limit) then
          warnings(324) = warnings(324) + 1
          call openOutput
          write(nx,*)
          write(nx,*) 'CODE ERROR 324 (loop_functions_rcl.f90):'
          write(nx,*) 'IR divergent DB0 requires values for Dir and muIR2'
          write(nx,*)
          call toomanywarnings(324)
        endif
        call istop (ifail,2)
      endif
    case (9)  ! m0=0
      DB0 = - 1/p2**2*( p2 + m12*cln(1d0-p2/m12,-1d0) )
    case (10) ! m1=0
      DB0 = - 1/p2**2*( p2 + m02*cln(1d0-p2/m02,-1d0) )
    case (11) ! the general case
      b  = - m12 + m02 - p2
      rt = sqrt(b**2-4*p2*m12)
      xp = (-b+rt)/(2*p2)
      xm = (-b-rt)/(2*p2)
      DB0 = 1/p2*(                                   &
            - 1d0                                    &
            + xp*(1d0-xp)/(xp-xm)*cln(1d0-1/xp,1d0)  &
            + xm*(1d0-xm)/(xm-xp)*cln(1d0-1/xm,-1d0) &
            )
    case default
      if (warnings(325).le.warning_limit) then
        warnings(325) = warnings(325) + 1
        call openOutput
        write(nx,*)
        write(nx,*) "CODE ERROR 325 (loop_functions_rcl.f90): DB0 case doesn't exist"
        write(nx,*)
        call toomanywarnings(325)
      endif
      call istop (ifail,2)
    end select
  endif

  end function DB0

!---------------------------------------------------------------------

  function B1 (p2,m02,m12,Duv,muUV2,Dir,muIR2)

  complex (dp)                        :: B1
  real    (dp),           intent (in) :: p2,Duv,muUV2
  complex (dp),           intent (in) :: m02,m12
  real    (dp), optional, intent (in) :: Dir,muIR2

  integer      :: type
  complex(dp)  :: B1_coli

  if (collier_ct) then
    call settype (p2,m02,m12,type)
    if (type.eq.1) then
      B1 = - 1/2d0 * ( DeltaUV - DeltaIR + log(muUV**2/muIR**2) + 1d0 )
    else
      B1 = B1_coli(cmplx(p2,kind=dp),m02,m12)
    endif
  else
    call settype (p2,m02,m12,type)
    select case (type)
    case (1)    ! p=m0=m1=0
      if (present(Dir).and.present(muIR2)) then
        B1 = - 1/2d0 * ( Duv - Dir + log(muUV2/muIR2) + 1d0 )
      else
        if (warnings(326).le.warning_limit) then
          warnings(326) = warnings(326) + 1
          call openOutput
          write(nx,*)
          write(nx,*) 'CODE ERROR 326 (loop_functions_rcl.f90):'
          write(nx,*) 'B1 with 0 scales requires values for Dir and muIR2'
          write(nx,*)
          call toomanywarnings(326)
        endif
        call istop (ifail,2)
      endif
    case (2:5)  ! p=0
      B1 = - 1/2d0 * B0(0d0,m02,m12,Duv,muUV2) &
           + (m12-m02)/2d0 * DB0(0d0,m02,m12)
    case (6:11) ! p ne 0
      B1 = 1/(2*p2)*(                                      &
           - ( p2 - m12 + m02 ) * B0(p2,m02,m12,Duv,muUV2) &
           + A0(m02,Duv,muUV2) - A0(m12,Duv,muUV2)         &
           )
    case default
      if (warnings(327).le.warning_limit) then
        warnings(327) = warnings(327) + 1
        call openOutput
        write(nx,*)
        write(nx,*) "CODE ERROR 327 (loop_functions_rcl.f90): B1 case doesn't exist"
        write(nx,*)
        call toomanywarnings(327)
      endif
      call istop (ifail,2)
    end select
  endif

  end function B1

!---------------------------------------------------------------------

  function DB1 (p2,m02,m12,Dir,muIR2)

  complex (dp)                        :: DB1
  real    (dp),           intent (in) :: p2
  complex (dp),           intent (in) :: m02,m12
  real    (dp), optional, intent (in) :: Dir,muIR2

  integer      :: type
  real    (dp) :: Duv,muUV2
  complex (dp) :: B0x,B1x,DB0x

  if (collier_ct) then
    call DB1_cll(DB1,cmplx(p2,kind=dp),m02,m12)
  else
    call settype (p2,m02,m12,type)
    select case (type)
    case (1)  ! p=m0=m1=0
      DB1 = c0d0
    case (2)  ! p=m0=0
      DB1 = - 1/(6*m12)
    case (3)  ! p=m1=0
      DB1 = - 1/(3*m02)
    case (4)  ! p=0,m0=m1
      DB1 = - 1/(12*m02)
    case (5)  ! p=0
      DB1 = 1/6d0/(m12-m02)**3*(                  &
            + 2*m02**2 + 5*m02*m12 - m12**2       &
            - 6*m02**2*m12/(m02-m12)*log(m02/m12) &
            )
    case (6:11) ! p neq 0, included IR case 7 p=m0,m1=0
      Duv = 0d0
      muUV2 = 1d0
      B0x = B0(p2,m02,m12,Duv,muUV2)
      B1x = B1(p2,m02,m12,Duv,muUV2)
      select case (type)
      case (7:8);   DB0x = DB0(p2,m02,m12,Dir,muIR2)
      case default; DB0x = DB0(p2,m02,m12)
      end select
      DB1 = - 1/(2*p2)*( ( p2 - m12 + m02 )*DB0x + 2*B1x + B0x )
    case default
      if (warnings(328).le.warning_limit) then
        warnings(328) = warnings(328) + 1
        call openOutput
        write(nx,*)
        write(nx,*) "CODE ERROR 328 (loop_functions_rcl.f90): DB1 case doesn't exist"
        write(nx,*)
        call toomanywarnings(328)
      endif
      call istop (ifail,2)
    end select
  endif

  end function DB1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module loop_functions_rcl

!#####################################################################

