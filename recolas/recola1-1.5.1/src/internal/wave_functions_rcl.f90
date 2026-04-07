!#####################################################################
!!
!!  File  wave_functions_rcl.f90
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

  module wave_functions_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine definewp (fl,p,pl,mass,he,wp)

  integer,      intent (in)  :: fl,he
  complex (dp), intent (in)  :: p(0:3),pl(4)
  real    (dp), intent (in)  :: mass
  complex (dp), intent (out) :: wp(0:3)

  real(dp)    :: s0,s3,p0p3
  complex(dp) :: m,pv,pT,r,php,phm,bp,bm,bhp,bhm,ap,am,mo,mo2
  logical     :: pTzero

  ! pl(1) = p(0)+p(3)
  ! pl(2) = p(0)-p(3)
  ! pl(3) = p(1)+i*p(2)
  ! pl(4) = p(1)-i*p(2)

  s0 = sign(1d0,real(p(0),kind=dp))
  s3 = sign(1d0,real(p(3),kind=dp))

  m  = cmplx(mass,kind=dp)
  pT = sqrt( pl(3)*pl(4) )     ! pT = \sqrt{p1^2+p2^2}
  pv = sqrt( pT*pT + p(3)*p(3) ) ! pv = |\vec p|

  pTzero  = abs(pT)/abs(p(0)) .lt. zerocut

  select case (parKind(fl))

  case (1)   ! scalars

    wp(0)   = c1d0
    wp(1:3) = c0d0

  case (2,3) ! vector bosons

    select case (he)

    case (0)

      mo2 = pl(1)*pl(2) - pl(3)*pl(4)
      mo = sqrt(mo2)
      wp(0) = pv/mo
      if (pTzero) then ! pT = 0 (i.e. p(1)/pT = s0, p(2)/pT = 0)
        wp(1:2) = c0d0
        wp(3)   = s3*p(0)/mo
      else             ! pT > 0 (this implies pv > 0)
        wp(1:3) = p(1:3)*p(0)/mo/pv
      endif

    case (-1,1)

      wp(0) = c0d0
      if (pTzero) then ! pT = 0 (i.e. p(1)/pT = s0, p(2)/pT = 0)
        wp(1) = - he*s3/csq2
        wp(2) = - cId0*s0/csq2
        wp(3) = c0d0
      else             ! pT > 0 (this implies pv > 0)
        wp(1) = ( - s0*he*p(1)*p(3)/pv + cId0*p(2) )/csq2/pT
        wp(2) = ( - s0*he*p(2)*p(3)/pv - cId0*p(1) )/csq2/pT
        wp(3) = s0*he*pT/csq2/pv
      endif

    case (2)

      mo2 = pl(1)*pl(2) - pl(3)*pl(4)
      mo = sqrt(mo2)
      wp(0) = p(0)/mo
      if (pTzero) then ! pT = 0 (i.e. p(1)/pT = s0, p(2)/pT = 0)
        wp(1:2) = c0d0
        wp(3)   = s3*pv/mo
      else             ! pT > 0 (this implies pv > 0)
        wp(1:3) = p(1:3)/mo
      endif

      wp(:) = wp(:)*sqrt(mo2/m**2-1d0)

    end select

  case (4,5) ! fermions - antifermions

    r  = sqrt(2*pv) ! r = \sqrt{2*|\vec p}}

    p0p3 = real(p(0),kind=dp) * real(p(3),kind=dp)

    ! php = \hat{p_+} = p_+/pT, phm = \hat{p_-} = p_-/pT
    ! bp = b_+ = \sqrt{pv+s0*p3}, bm = b_- = \sqrt{pv-s0*p3}
    ! bhp = \hat{b_+} = b_+/r, bhm = \hat{b_-} = b_-/r
    if (pTzero) then ! pT = 0 (i.e. p(1)/pT = s0, p(2)/pT = 0)
      php = s0
      phm = s0
      if (p0p3.gt.0d0) then     ! p(0)*p(3) > 0
        bp = r
        bm = c0d0
        bhp = c1d0
        bhm = c0d0
      elseif (p0p3.lt.0d0) then ! p(0)*p(3) < 0
        bp = c0d0
        bm = r
        bhp = c0d0
        bhm = c1d0
      else                      ! p(0)*p(3) = 0
        bp = c0d0
        bm = c0d0
        bhp = 1/csq2
        bhm = 1/csq2
      endif
    else             ! pT > 0 (this implies pv > 0)
      php = pl(3)/pT
      phm = pl(4)/pT
      if (p0p3.gt.0d0) then     ! p(0)*p(3) > 0
        bp = sqrt( pv + s0*p(3) ); bm = pT/bp
      elseif (p0p3.lt.0d0) then ! p(0)*p(3) < 0
        bm = sqrt( pv - s0*p(3) ); bp = pT/bm
      else                      ! p(0)*p(3) = 0
        bp = sqrt(pv); bm = bp
      endif
      bhp = bp/r
      bhm = bm/r
    endif

    if (regf(fl).le.2) then  ! massless fermions - antifermions

      select case (he)
      case ( 1)
        select case (parKind(fl))
        case (4) ! fermions
          wp(0)   = + bp
          wp(1)   = + s0 * php * bm
        case (5) ! antifermions
          wp(0)   = + s0 * php * bm
          wp(1)   = - bp
        end select
        wp(2:3) = c0d0
      case (-1)
        wp(0:1) = c0d0
        select case (parKind(fl))
        case (4) ! fermions
          wp(2)   = - s0 * phm * bm
          wp(3)   = + bp
        case (5) ! antifermions
          wp(2)   = - bp
          wp(3)   = - s0 * phm * bm
        end select
      case default
        call errorFlag(1)
      end select

    else                     ! massive fermions - antifermions

      ! ap = a_+ = \sqrt{|p0|+pv}, am = a_- = \sqrt{|p0|-pv}
      ap = sqrt( s0*p(0) + pv )
      am = m/ap

      select case (he)
      case ( 1)
        select case (parKind(fl))
        case (4) ! fermions
          wp(0) = + ap * bhp
          wp(1) = + s0 * php * ap * bhm
          wp(2) = - s0 * am * bhp
          wp(3) = - php * am * bhm
        case (5) ! antifermions
          wp(0) = + s0 * php * ap * bhm
          wp(1) = - ap * bhp
          wp(2) = + php * am * bhm
          wp(3) = - s0 * am * bhp
        end select
      case (-1)
        select case (parKind(fl))
        case (4) ! fermions
          wp(0) = + phm * am * bhm
          wp(1) = - s0 * am * bhp
          wp(2) = - s0 * phm * ap * bhm
          wp(3) = + ap * bhp
        case (5) ! antifermions
          wp(0) = - s0 * am * bhp
          wp(1) = - phm * am * bhm
          wp(2) = - ap * bhp
          wp(3) = - s0 * phm * ap * bhm
        end select
      case default
        call errorFlag(1)
      end select

    endif

  case default

    call errorFlag(2)

  end select

  contains

    subroutine errorFlag(flag)

    integer, intent (in) :: flag

    if (warnings(336).le.warning_limit) then
      warnings(336) = warnings(336) + 1
      call openOutput
      write(nx,*)
      select case (flag)
      case (1)
        write(nx,*) 'CODE ERROR 366 (wave_functions_rcl.f90):'
        write(nx,*) 'wrong value for the helicity',he
      case (2)
        write(nx,*) 'CODE ERROR 366 (wave_functions_rcl.f90):'
        write(nx,*) 'wrong value for the external particle',fl
      case default
        write(nx,*) 'CODE ERROR 366 (wave_functions_rcl.f90):'
        write(nx,*) 'undefined error flag', flag
      end select
      write(nx,*)
      call toomanywarnings(336)
    endif
    call istop (ifail,2)

    end subroutine errorFlag

  end subroutine definewp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module wave_functions_rcl

!#####################################################################


