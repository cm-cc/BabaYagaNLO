!#####################################################################
!!
!!  File  tree_vertices_rcl.f90
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

  module tree_vertices_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 2 legs
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine tree2 ( ps,p,pl,m2,den,cop,ty,wp1,wp2 )

  integer,      intent (in)  :: ty
  complex (dp), intent (in)  :: ps,m2,den,cop(4),p(0:3),pl(1:4), &
                                wp1(0:3)
  complex (dp), intent (out) :: wp2(0:3)

  complex (dp) :: cp(4),cc,spr,m,fac1,fac2,fac3,fac4

  select case (ty)
  case (1021,1041)
    if (CMscheme.eq.1) then
      m = sqrt(m2)
    else
      m = sqrt(real(m2,kind=dp))*c1d0
    endif
  end select

  ! A 2-point vertex is never the last vertex
  cp = - cop/den

  select case (ty)

  case (1001) ! s -> s
  ! currents: wp1 -> wp2
  ! coupling = i*( cop(1)*ps - cop(2) ), propagator = i/(ps-m2)

    wp2(0) = (cp(1)*ps - cp(2))*wp1(0)
    wp2(1:3) = c0d0

  case (1002) ! v -> s
  ! currents: wp1^mu -> wp2
  ! coupling = i*cop(1)*p_mu, propagator = i/(ps-m2)

    spr = p(0)*wp1(0) - sum(p(1:3)*wp1(1:3))

    wp2(0) = cp(1) * spr
    wp2(1:3) = c0d0

  case (1011) ! s -> v
  ! currents: wp1 -> wp2^mu
  ! coupling = - i*cop(1)*p_mu, propagator = - i*g_nu^mu/(ps-m2)

    cc = wp1(0) * cp(1)

    wp2 = cc * p

  case (1012) ! v -> v
  ! currents: wp1^mu -> wp2^nu
  ! coupling = - i*( ( cop(1)*p^2 - cop(2) )*g_{mu,ro} + cop(3)*p_mu*p_ro )
  ! propagator = - i*g^{ro,nu}/(ps-m2)

    spr = p(0)*wp1(0) - sum(p(1:3)*wp1(1:3))

    wp2 = (cp(1)*ps - cp(2))*wp1 + cp(3)*spr*p

  case (1021) ! f -> f, massive
  ! currents: wp1(i) -> wp2(j)
  ! coupling = i*(
  !            + cop(1)*pslash*omega_-(k,i)
  !            + cop(2)*pslash*omega_+(k,i)
  !            + cop(3)*omega_-(k,i)
  !            + cop(4)*omega_+(k,i)
  !            )
  ! propagator = i*(pslash+m)(j,k)/(ps-m2); pslash = p^mu*gamma_mu

    fac1 = + cp(1)*m + cp(3)
    fac2 = + cp(2)*m + cp(4)
    fac3 = + cp(3)*m + cp(1)*ps
    fac4 = + cp(4)*m + cp(2)*ps

    wp2(0:1) = fac4*wp1(0:1)
    wp2(2:3) = fac3*wp1(2:3)

    wp2(0) = wp2(0) + fac1*( - pl(1)*wp1(2) - pl(4)*wp1(3) )
    wp2(1) = wp2(1) + fac1*( - pl(3)*wp1(2) - pl(2)*wp1(3) )
    wp2(2) = wp2(2) + fac2*( - pl(2)*wp1(0) + pl(4)*wp1(1) )
    wp2(3) = wp2(3) + fac2*( + pl(3)*wp1(0) - pl(1)*wp1(1) )

  case (1022) ! f -> f, massless (m2=0, c3=c4=0)
  ! currents: wp1(i) -> wp2(j)
  ! coupling = i*(
  !            + cop(1)*pslash*omega_-(k,i)
  !            + cop(2)*pslash*omega_+(k,i)
  !            )
  ! propagator = i*pslash(j,k)/ps; pslash = p^mu*gamma_mu

    fac3 = + cp(1)*ps
    fac4 = + cp(2)*ps

    wp2(0:1) = fac4*wp1(0:1)
    wp2(2:3) = fac3*wp1(2:3)

  case (1023) ! f_- -> f_-, massless (m2=0, c3=c4=0)
  ! currents: wp1(i) -> wp2(j)
  ! coupling = i*cop(1)*pslash*omega_-(k,i)
  ! propagator = i*pslash(j,k)/ps; pslash = p^mu*gamma_mu

    fac3 = + cp(1)*ps

    wp2(0:1) = c0d0
    wp2(2:3) = fac3*wp1(2:3)

  case (1024) ! f_+ -> f_+, massless (m2=0, c3=c4=0)
  ! currents: wp1(i) -> wp2(j)
  ! coupling = i*cop(2)*pslash*omega_+(k,i)
  ! propagator = i*pslash(j,k)/ps; pslash = p^mu*gamma_mu

    fac4 = + cp(2)*ps

    wp2(0:1) = fac4*wp1(0:1)
    wp2(2:3) = c0d0

  case (1041) ! f~ -> f~, massive
  ! currents: wp1(i) -> wp2(j)
  ! coupling = i*(
  !            - cop(1)*pslash*omega_-(i,k)
  !            - cop(2)*pslash*omega_+(i,k)
  !            + cop(3)*omega_-(i,k)
  !            + cop(4)*omega_+(i,k)
  !            )
  ! propagator = i*(-pslash+m)(k,j)/(ps-m2); pslash = p^mu*gamma_mu

    fac1 = - cp(1)*m - cp(4)
    fac2 = - cp(2)*m - cp(3)
    fac3 = + cp(3)*m + cp(2)*ps
    fac4 = + cp(4)*m + cp(1)*ps

    wp2(0:1) = fac4*wp1(0:1)
    wp2(2:3) = fac3*wp1(2:3)

    wp2(0) = wp2(0) + fac2*( - pl(2)*wp1(2) + pl(3)*wp1(3) )
    wp2(1) = wp2(1) + fac2*( + pl(4)*wp1(2) - pl(1)*wp1(3) )
    wp2(2) = wp2(2) + fac1*( - pl(1)*wp1(0) - pl(3)*wp1(1) )
    wp2(3) = wp2(3) + fac1*( - pl(4)*wp1(0) - pl(2)*wp1(1) )

  case (1042) ! f~ -> f~, massless (m2=0, c3=c4=0)
  ! currents: wp1(i) -> wp2(j)
  ! coupling = i*(
  !            - cop(1)*pslash*omega_-(i,k)
  !            - cop(2)*pslash*omega_+(i,k)
  !            )
  ! propagator = -i*pslash(k,j)/ps; pslash = p^mu*gamma_mu

    fac3 = + cp(2)*ps
    fac4 = + cp(1)*ps

    wp2(0:1) = fac4*wp1(0:1)
    wp2(2:3) = fac3*wp1(2:3)

  case (1043) ! f~_- -> f~_-, massless (m2=0, c3=c4=0)
  ! currents: wp1(i) -> wp2(j)
  ! coupling = - i*cop(1)*pslash*omega_-(i,k)
  ! propagator = -i*pslash(k,j)/ps; pslash = p^mu*gamma_mu

    fac3 = + cp(2)*ps

    wp2(0:1) = c0d0
    wp2(2:3) = fac3*wp1(2:3)

  case (1044) ! f~_+ -> f~_+, massless (m2=0, c3=c4=0)
  ! currents: wp1(i) -> wp2(j)
  ! coupling = - i*cop(1)*pslash*omega_-(i,k)
  ! propagator = -i*pslash(k,j)/ps; pslash = p^mu*gamma_mu

    fac4 = + cp(1)*ps

    wp2(0:1) = fac4*wp1(0:1)
    wp2(2:3) = c0d0

  case default

    if (warnings(305).le.warning_limit) then
      warnings(305) = warnings(305) + 1
      call openOutput
      write(nx,*)
      write(nx,*) "CODE ERROR 305 (tree_vertices_rcl.f90): wrong 2-leg interaction"
      write(nx,*)
      call toomanywarnings(305)
    endif
    call istop (ifail,2)

 end select

  end subroutine tree2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 3 legs
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine tree3 ( last,p1,p2,pl1,pl2,m2,den,cop,ty,wp1,wp2,wp3 )

  logical,      intent (in)  :: last
  integer,      intent (in)  :: ty
  complex (dp), intent (in)  :: m2,den,cop(2),p1(0:3),p2(0:3),     &
                                pl1(1:4),pl2(1:4),wp1(0:3),wp2(0:3)
  complex (dp), intent (out) :: wp3(0:3)

  complex (dp) :: cp(2),cc1,cc2,cc,wpr,spr,spr1,spr2,              &
                  p13(0:3),p23(0:3),p12(0:3),pl(4),                &
                  prd1(0:3),prd2(0:3),                             &
                  wab02,wab13,wab12,wab03,wab20,wab31,wab21,wab30, &
                  auxa,auxb,auxc,auxd,m,wap(4)

  ! set overall factors for currents

  if (.not.last) then
    select case (ty)
    case (21,23,31,33,41,43,51,53)
      if (CMscheme.eq.1) then
        m = sqrt(m2)
      else
        m = sqrt(real(m2,kind=dp))*c1d0
      endif
    end select
  endif

  select case (ty)
  case (1:4,7:8,11:14,17:18,25:26,35:36,45:46,57:58)
 ! s + s -> s, v + v -> s, s + v -> s,  v + s -> s
 ! f_- + f~_- -> s,  f~_- + f_- -> s
 ! f_- + f~_- -> s,  f~_- + f_- -> s
 ! s + s -> v, v + v -> v, s + v -> v, v + s -> v
 ! f_- + f~_+ -> v, f~_+ + f_- -> v
 ! f_- + s -> f_+/- massless, s + f_- -> f_+/- massless
 ! f_- + v -> f_-/+ massless, v + f_- -> f_-/+ massless
 ! f~_- + s -> f~_+/- massless, s + f~_- -> f~_+/- massless
 ! f~_+ + v -> f~_+/- massless, v + f~_+ -> f~_+/- massless
    if (last) then
      cp(1) = cId0*cop(1)
    else
      cp(1) = - cop(1)/den
    endif
  case (9:10,19:20,27:28,37:38,47:48,55:56)
  ! f_+ + f~_+ -> s, f~_+ + f_+ -> s
  ! f_+ + f~_- -> v, f~_- f_+ -> v
  ! f_+ + s -> f_-/+ massless
  ! s + f_+ -> f_-/+ massless
  ! f_+ + v -> f_+/- massless
  ! v + f_+ -> f_+/- massless
  ! f~_+ + s -> f~_-/+ massless
  ! s + f~_+ -> f~_-/+ massless
  ! f~_- + v -> f~_-/+ massless
  ! v + f~_- -> f~_-/+ massless
    if (last) then
      cp(2) = cId0*cop(2)
    else
      cp(2) = - cop(2)/den
    endif
  case (5:6,15:16,21:24,31:34,41:44,51:54)
 ! f + f~ -> s,  f~ + f -> s
 ! f + f~ -> v, f~ + f -> v
 ! f + s -> f massive, f + s -> f massless
 ! s + f -> f massive, s + f -> f massless
 ! f + v -> f massive, f + v -> f massless
 ! v + f -> f massive, v + f -> f massless
 ! f~ + s -> f~ massive, f~ + s -> f~ massless
 ! s + f~ -> f~ massive, s + f~ -> f~ massless
 ! f~ + v -> f~ massive, f~ + v -> f~ massless
 ! v + f~ -> f~ massive, v + f~ -> f~ massless
    if (last) then
      cp = cId0*cop
    else
      cp = - cop/den
    endif
  end select


  ! select interaction type and build currents

  select case (ty)

  case (1) ! s + s -> s
  ! currents: wp1, wp2 -> wp3
  ! coupling = i*cop(1), propagator = i/(ps-m2)

    wp3(0) = cp(1) * wp1(0) * wp2(0)
    wp3(1:3) = c0d0

  case (2) ! v + v -> s
  ! currents: wp1^mu, wp2^nu -> wp3
  ! coupling = i*cop(1)*g_{mu,nu}, propagator = i/(ps-m2)

    wpr = wp1(0)*wp2(0) - sum(wp1(1:3)*wp2(1:3))
    wp3(0) = cp(1) * wpr
    wp3(1:3) = c0d0

  case (3) ! s + v -> s
  ! currents: wp1, wp2^mu -> wp3
  ! coupling = i*cop(1)*(p1-p3)_mu, propagator = i/(ps-m2)

    p13 = 2*p1 + p2 ! p13 = p1 - p3
    spr = p13(0)*wp2(0) - sum(p13(1:3)*wp2(1:3))
    wp3(0) = cp(1) * wp1(0) * spr
    wp3(1:3) = c0d0

  case (4) ! v + s -> s
  ! currents: wp1^mu, wp2 -> wp3
  ! coupling = i*cop(1)*(p3-p2)_mu, propagator = i/(ps-m2)

    p23 = - 2*p2 - p1 ! p23 = p3 - p2
    spr = p23(0)*wp1(0) - sum(p23(1:3)*wp1(1:3))
    wp3(0) = cp(1) * wp2(0) * spr
    wp3(1:3) = c0d0

  case (5,6) ! f + f~ -> s,  f~ + f -> s
  ! currents: wp1(i), wp2(j) -> wp3
  ! coupling = i*( cop(1)*omega_-(i,j) + cop(2)*omega_+(i,j) )
  ! propagator = i/(ps-m2)

    spr1 = wp1(2)*wp2(2) + wp1(3)*wp2(3)
    spr2 = wp1(0)*wp2(0) + wp1(1)*wp2(1)
    wp3(0) = cp(1)*spr1 + cp(2)*spr2
    wp3(1:3) = c0d0

  case (7:8) ! f_- + f~_- -> s,  f~_- + f_- -> s
  ! currents: wp1(i), wp2(j) -> wp3 (wp1(0:1)=wp2(0:1)=0)
  ! coupling = i*( cop(1)*omega_-(i,j) + cop(2)*omega_+(i,j) )
  ! propagator = i/(ps-m2)

    spr1 = wp1(2)*wp2(2) + wp1(3)*wp2(3)
    wp3(0) = cp(1) * spr1
    wp3(1:3) = c0d0

  case (9:10) ! f_+ + f~_+ -> s,  f~_+ + f_+ -> s
  ! currents: wp1(i), wp2(j) -> wp3 (wp1(2:3)=wp2(2:3)=0)
  ! coupling = i*( cop(1)*omega_-(i,j) + cop(2)*omega_+(i,j) )
  ! propagator = i/(ps-m2)

    spr2 = wp1(0)*wp2(0) + wp1(1)*wp2(1)
    wp3(0) = cp(2) * spr2
    wp3(1:3) = c0d0

  case (11) ! s + s -> v
  ! currents: wp1, wp2 -> wp3^mu
  ! coupling = i*cop(1)*(p2-p1)^nu, propagator = - i*g_nu^mu/(ps-m2)

    p12 = p1 - p2
    cc1 = wp1(0) * wp2(0) * cp(1)
    wp3 = cc1 * p12
    if (last) wp3(0) = - wp3(0)

  case (12) ! v + v -> v
  ! currents: wp1^mu, wp2^nu -> wp3^ro
  ! coupling = - i*(-cop(1))*(
  !              + g_{mu,nu}*(p2-p1)_si
  !              + g_{nu,si}*(p3-p2)_mu
  !              + g_{si,mu}*(p1-p3)_nu
  !              )
  ! propagator = - i*g^{si,ro}/(ps-m2)

    p12 = p1 - p2
    p23 = 2*p2 + p1 ! p23 = p2 - p3
    p13 = 2*p1 + p2 ! p13 = p1 - p3
    spr1 = p23(0)*wp1(0) - sum(p23(1:3)*wp1(1:3))
    spr2 = p13(0)*wp2(0) - sum(p13(1:3)*wp2(1:3))
    wpr  = wp2(0)*wp1(0) - sum(wp2(1:3)*wp1(1:3))
    cc1 = - cp(1) * spr2
    cc2 =   cp(1) * spr1
    cc  =   cp(1) * wpr
    wp3 = cc1*wp1 + cc2*wp2 + cc*p12
    if (last) wp3(0) = - wp3(0)

  case (13) ! s + v -> v
  ! currents: wp1, wp2^mu -> wp3^nu
  ! coupling = i*cop(1)*g_{mu,ro}, propagator = - i*g^{ro,nu}/(ps-m2)

    cc2 = - cp(1) * wp1(0)
    wp3 = cc2 * wp2
    if (last) wp3(0) = - wp3(0)

  case (14) ! v + s -> v
  ! currents: wp1^mu, wp2 -> wp3^nu
  ! coupling = i*cop(1)*g_{mu,ro}, propagator = - i*g^{ro,nu}/(ps-m2)

    cc1 = - cp(1)*wp2(0)
    wp3 = cc1 * wp1
    if (last) wp3(0) = - wp3(0)

  case (15) ! f + f~ -> v
  ! currents: wp1(i), wp2(j) -> wp3^mu
  ! coupling = i*gamma^ro(j,k)*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = - i*g_ro^mu/(ps-m2)

    wab02 = wp2(0)*wp1(2)
    wab13 = wp2(1)*wp1(3)
    wab12 = wp2(1)*wp1(2)
    wab03 = wp2(0)*wp1(3)
    wab20 = wp2(2)*wp1(0)
    wab31 = wp2(3)*wp1(1)
    wab21 = wp2(2)*wp1(1)
    wab30 = wp2(3)*wp1(0)
    prd1(0) =     wab02 + wab13
    prd1(1) =   - wab12 - wab03
    prd1(2) = ( - wab12 + wab03 )*cId0
    prd1(3) =   - wab02 + wab13
    prd2(0) =     wab20 + wab31
    prd2(1) =     wab30 + wab21
    prd2(2) = (   wab30 - wab21 )*cId0
    prd2(3) =     wab20 - wab31

    wp3 = cp(1)*prd1 + cp(2)*prd2
    if (last) wp3(0) = - wp3(0)

  case (16) ! f~ + f -> v
  ! currents: wp1(i), wp2(j) -> wp3^mu
  ! coupling = i*gamma^ro(i,k)*( cop(1)*omega_-(k,j) + cop(2)*omega_+(k,j) )
  ! propagator = - i*g_ro^mu/(ps-m2)

    wab02 = wp1(0)*wp2(2)
    wab13 = wp1(1)*wp2(3)
    wab12 = wp1(1)*wp2(2)
    wab03 = wp1(0)*wp2(3)
    wab20 = wp1(2)*wp2(0)
    wab31 = wp1(3)*wp2(1)
    wab21 = wp1(2)*wp2(1)
    wab30 = wp1(3)*wp2(0)
    prd1(0) =     wab02 + wab13
    prd1(1) =   - wab12 - wab03
    prd1(2) = ( - wab12 + wab03 )*cId0
    prd1(3) =   - wab02 + wab13
    prd2(0) =     wab20 + wab31
    prd2(1) =     wab30 + wab21
    prd2(2) = (   wab30 - wab21 )*cId0
    prd2(3) =     wab20 - wab31

    wp3 = cp(1)*prd1 + cp(2)*prd2
    if (last) wp3(0) = - wp3(0)

  case (17) ! f_- + f~_+ -> v
  ! currents: wp1(i), wp2(j) -> wp3^mu (wp1(0:1)=wp2(2:3)=0)
  ! coupling = i*gamma^ro(j,k)*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = - i*g_ro^mu/(ps-m2)

    wab02 = wp2(0)*wp1(2)
    wab13 = wp2(1)*wp1(3)
    wab12 = wp2(1)*wp1(2)
    wab03 = wp2(0)*wp1(3)
    prd1(0) =     wab02 + wab13
    prd1(1) =   - wab12 - wab03
    prd1(2) = ( - wab12 + wab03 )*cId0
    prd1(3) =   - wab02 + wab13

    wp3 = cp(1) * prd1
    if (last) wp3(0) = - wp3(0)

  case (18) ! f~_+ + f_- -> v
  ! currents: wp1(i), wp2(j) -> wp3^mu (wp1(2:3)=wp2(0:1)=0)
  ! coupling = i*gamma^ro(i,k)*( cop(1)*omega_-(k,j) + cop(2)*omega_+(k,j) )
  ! propagator = - i*g_ro^mu/(ps-m2)

    wab02 = wp1(0)*wp2(2)
    wab13 = wp1(1)*wp2(3)
    wab12 = wp1(1)*wp2(2)
    wab03 = wp1(0)*wp2(3)
    prd1(0) =     wab02 + wab13
    prd1(1) =   - wab12 - wab03
    prd1(2) = ( - wab12 + wab03 )*cId0
    prd1(3) =   - wab02 + wab13

    wp3 = cp(1) * prd1
    if (last) wp3(0) = - wp3(0)

  case (19) ! f_+ + f~_- -> v
  ! currents: wp1(i), wp2(j) -> wp3^mu (wp1(2:3)=wp2(0:1)=0)
  ! coupling = i*gamma^ro(j,k)*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = - i*g_ro^mu/(ps-m2)

    wab20 = wp2(2)*wp1(0)
    wab31 = wp2(3)*wp1(1)
    wab21 = wp2(2)*wp1(1)
    wab30 = wp2(3)*wp1(0)
    prd2(0) =   wab20 + wab31
    prd2(1) =   wab30 + wab21
    prd2(2) = ( wab30 - wab21 )*cId0
    prd2(3) =   wab20 - wab31

    wp3 = cp(2) * prd2
    if (last) wp3(0) = - wp3(0)

  case (20) ! f~_- f_+ -> v
  ! currents: wp1(i), wp2(j) -> wp3^mu (wp1(0:1)=wp2(2:3)=0)
  ! coupling = i*gamma^ro(i,k)*( cop(1)*omega_-(k,j) + cop(2)*omega_+(k,j) )
  ! propagator = - i*g_ro^mu/(ps-m2)

    wab20 = wp1(2)*wp2(0)
    wab31 = wp1(3)*wp2(1)
    wab21 = wp1(2)*wp2(1)
    wab30 = wp1(3)*wp2(0)
    prd2(0) =   wab20 + wab31
    prd2(1) =   wab30 + wab21
    prd2(2) = ( wab30 - wab21 )*cId0
    prd2(3) =   wab20 - wab31

    wp3 = cp(2) * prd2
    if (last) wp3(0) = - wp3(0)

  case (21) ! f + s -> f, massive
  ! currents: wp1(i), wp2 -> wp3(j)
  ! coupling = i*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*(pslash+m)(j,k)/(ps-m2); pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxa = cp(2)*wp1(0)*wp2(0)
    auxb = cp(2)*wp1(1)*wp2(0)
    auxc = cp(1)*wp1(2)*wp2(0)
    auxd = cp(1)*wp1(3)*wp2(0)
    select case (last)
    case (.false.)
      wp3(0)= - auxc*pl(1) - auxd*pl(4) + auxa*m
      wp3(1)= - auxc*pl(3) - auxd*pl(2) + auxb*m
      wp3(2)= - auxa*pl(2) + auxb*pl(4) + auxc*m
      wp3(3)=   auxa*pl(3) - auxb*pl(1) + auxd*m
    case (.true.)
      wp3(0)= auxa
      wp3(1)= auxb
      wp3(2)= auxc
      wp3(3)= auxd
    end select

  case (22) ! f + s -> f, massless
  ! currents: wp1(i), wp2 -> wp3(j)
  ! coupling = i*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*pslash(j,k)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxa = cp(2)*wp1(0)*wp2(0)
    auxb = cp(2)*wp1(1)*wp2(0)
    auxc = cp(1)*wp1(2)*wp2(0)
    auxd = cp(1)*wp1(3)*wp2(0)
    select case (last)
    case (.false.)
      wp3(0)= - auxc*pl(1) - auxd*pl(4)
      wp3(1)= - auxc*pl(3) - auxd*pl(2)
      wp3(2)= - auxa*pl(2) + auxb*pl(4)
      wp3(3)=   auxa*pl(3) - auxb*pl(1)
    case (.true.)
      wp3(0)= auxa
      wp3(1)= auxb
      wp3(2)= auxc
      wp3(3)= auxd
    end select

  case (23) ! s + f -> f, massive
  ! currents: wp1, wp2(i) -> wp3(j)
  ! coupling = i*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*(pslash+m)(j,k)/(ps-m2); pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxa = cp(2)*wp2(0)*wp1(0)
    auxb = cp(2)*wp2(1)*wp1(0)
    auxc = cp(1)*wp2(2)*wp1(0)
    auxd = cp(1)*wp2(3)*wp1(0)
    select case (last)
    case (.false.)
      wp3(0)= - auxc*pl(1) - auxd*pl(4) + auxa*m
      wp3(1)= - auxc*pl(3) - auxd*pl(2) + auxb*m
      wp3(2)= - auxa*pl(2) + auxb*pl(4) + auxc*m
      wp3(3)=   auxa*pl(3) - auxb*pl(1) + auxd*m
    case (.true.)
      wp3(0)= auxa
      wp3(1)= auxb
      wp3(2)= auxc
      wp3(3)= auxd
    end select

  case (24) ! s + f -> f, massless
  ! currents: wp1, wp2(i) -> wp3(j)
  ! coupling = i*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*pslash(j,k)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxa = cp(2)*wp2(0)*wp1(0)
    auxb = cp(2)*wp2(1)*wp1(0)
    auxc = cp(1)*wp2(2)*wp1(0)
    auxd = cp(1)*wp2(3)*wp1(0)
    select case (last)
    case (.false.)
      wp3(0)= - auxc*pl(1) - auxd*pl(4)
      wp3(1)= - auxc*pl(3) - auxd*pl(2)
      wp3(2)= - auxa*pl(2) + auxb*pl(4)
      wp3(3)=   auxa*pl(3) - auxb*pl(1)
    case (.true.)
      wp3(0)= auxa
      wp3(1)= auxb
      wp3(2)= auxc
      wp3(3)= auxd
    end select

  case (25) ! f_- + s -> f_+/-, massless
  ! currents: wp1(i), wp2 -> wp3(j) (wp1(0:1)=0)
  ! coupling = i*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*pslash(j,k)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxc = cp(1)*wp1(2)*wp2(0)
    auxd = cp(1)*wp1(3)*wp2(0)
    select case (last)
    case (.false.)
      wp3(0) = - auxc*pl(1) - auxd*pl(4)
      wp3(1) = - auxc*pl(3) - auxd*pl(2)
      wp3(2:3) = c0d0
    case (.true.)
      wp3(0:1) = c0d0
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (26) ! s + f_- -> f_+/-, massless
  ! currents: wp1, wp2(i) -> wp3(j) (wp2(0:1)=0)
  ! coupling = i*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*pslash(j,k)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxc = cp(1)*wp2(2)*wp1(0)
    auxd = cp(1)*wp2(3)*wp1(0)
    select case (last)
    case (.false.)
      wp3(0) = - auxc*pl(1) - auxd*pl(4)
      wp3(1) = - auxc*pl(3) - auxd*pl(2)
      wp3(2:3) = c0d0
    case (.true.)
      wp3(0:1) = c0d0
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (27) ! f_+ + s -> f_-/+, massless
  ! currents: wp1(i), wp2 -> wp3(j) (wp1(2:3)=0)
  ! coupling = i*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*pslash(j,k)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxa = cp(2)*wp1(0)*wp2(0)
    auxb = cp(2)*wp1(1)*wp2(0)
    select case (last)
    case (.false.)
      wp3(0:1) = c0d0
      wp3(2) = - auxa*pl(2) + auxb*pl(4)
      wp3(3) =   auxa*pl(3) - auxb*pl(1)
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2:3) = c0d0
    end select

  case (28) ! s + f_+ -> f_-/+, massless
  ! currents: wp1, wp2(i) -> wp3(j) (wp2(2:3)=0)
  ! coupling = i*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*pslash(j,k)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxa = cp(2)*wp2(0)*wp1(0)
    auxb = cp(2)*wp2(1)*wp1(0)
    select case (last)
    case (.false.)
      wp3(0:1) = c0d0
      wp3(2) = - auxa*pl(2) + auxb*pl(4)
      wp3(3) =   auxa*pl(3) - auxb*pl(1)
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2:3) = c0d0
    end select

  case (31) ! f + v -> f, massive
  ! currents: wp1(i), wp2^mu -> wp3(j)
  ! coupling = i*gamma_mu(l,k)*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*(pslash+m)(j,l)/(ps-m2); pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp2(0) + wp2(3)
    wap(2) = wp2(0) - wp2(3)
    wap(3) = wp2(1) + cId0*wp2(2)
    wap(4) = wp2(1) - cId0*wp2(2)
    auxa = cp(1)*( - wp1(2)*wap(1) - wp1(3)*wap(4) )
    auxb = cp(1)*( - wp1(2)*wap(3) - wp1(3)*wap(2) )
    auxc = cp(2)*( - wp1(0)*wap(2) + wp1(1)*wap(4) )
    auxd = cp(2)*(   wp1(0)*wap(3) - wp1(1)*wap(1) )
    select case (last)
    case (.false.)
      wp3(0) = - auxc*pl(1) - auxd*pl(4) + auxa*m
      wp3(1) = - auxc*pl(3) - auxd*pl(2) + auxb*m
      wp3(2) = - auxa*pl(2) + auxb*pl(4) + auxc*m
      wp3(3) =   auxa*pl(3) - auxb*pl(1) + auxd*m
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (32) ! f + v -> f, massless
  ! currents: wp1(i), wp2^mu -> wp3(j)
  ! coupling = i*gamma_mu(l,k)*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*pslash(j,l)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp2(0) + wp2(3)
    wap(2) = wp2(0) - wp2(3)
    wap(3) = wp2(1) + cId0*wp2(2)
    wap(4) = wp2(1) - cId0*wp2(2)
    auxa = cp(1)*( - wp1(2)*wap(1) - wp1(3)*wap(4) )
    auxb = cp(1)*( - wp1(2)*wap(3) - wp1(3)*wap(2) )
    auxc = cp(2)*( - wp1(0)*wap(2) + wp1(1)*wap(4) )
    auxd = cp(2)*(   wp1(0)*wap(3) - wp1(1)*wap(1) )
    select case (last)
    case (.false.)
      wp3(0) = - auxc*pl(1) - auxd*pl(4)
      wp3(1) = - auxc*pl(3) - auxd*pl(2)
      wp3(2) = - auxa*pl(2) + auxb*pl(4)
      wp3(3) =   auxa*pl(3) - auxb*pl(1)
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (33) ! v + f -> f, massive
  ! currents: wp1^mu, wp2(i) -> wp3(j)
  ! coupling = i*gamma_mu(l,k)*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*(pslash+m)(j,l)/(ps-m2); pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp1(0) + wp1(3)
    wap(2) = wp1(0) - wp1(3)
    wap(3) = wp1(1) + cId0*wp1(2)
    wap(4) = wp1(1) - cId0*wp1(2)
    auxa = cp(1)*( - wp2(2)*wap(1) - wp2(3)*wap(4) )
    auxb = cp(1)*( - wp2(2)*wap(3) - wp2(3)*wap(2) )
    auxc = cp(2)*( - wp2(0)*wap(2) + wp2(1)*wap(4) )
    auxd = cp(2)*(   wp2(0)*wap(3) - wp2(1)*wap(1) )
    select case (last)
    case (.false.)
      wp3(0) = - auxc*pl(1) - auxd*pl(4) + auxa*m
      wp3(1) = - auxc*pl(3) - auxd*pl(2) + auxb*m
      wp3(2) = - auxa*pl(2) + auxb*pl(4) + auxc*m
      wp3(3) =   auxa*pl(3) - auxb*pl(1) + auxd*m
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (34) ! v + f -> f, massless
  ! currents: wp1^mu, wp2(i) -> wp3(j)
  ! coupling = i*gamma_mu(l,k)*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*pslash(j,l)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp1(0) + wp1(3)
    wap(2) = wp1(0) - wp1(3)
    wap(3) = wp1(1) + cId0*wp1(2)
    wap(4) = wp1(1) - cId0*wp1(2)
    auxa = cp(1)*( - wp2(2)*wap(1) - wp2(3)*wap(4) )
    auxb = cp(1)*( - wp2(2)*wap(3) - wp2(3)*wap(2) )
    auxc = cp(2)*( - wp2(0)*wap(2) + wp2(1)*wap(4) )
    auxd = cp(2)*(   wp2(0)*wap(3) - wp2(1)*wap(1) )
    select case (last)
    case (.false.)
      wp3(0) = - auxc*pl(1) - auxd*pl(4)
      wp3(1) = - auxc*pl(3) - auxd*pl(2)
      wp3(2) = - auxa*pl(2) + auxb*pl(4)
      wp3(3) =   auxa*pl(3) - auxb*pl(1)
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (35) ! f_- + v -> f_-/+, massless
  ! currents: wp1(i), wp2^mu -> wp3(j) (wp1(0:1)=0)
  ! coupling = i*gamma_mu(l,k)*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*pslash(j,l)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp2(0) + wp2(3)
    wap(2) = wp2(0) - wp2(3)
    wap(3) = wp2(1) + cId0*wp2(2)
    wap(4) = wp2(1) - cId0*wp2(2)
    auxa = cp(1)*( - wp1(2)*wap(1) - wp1(3)*wap(4) )
    auxb = cp(1)*( - wp1(2)*wap(3) - wp1(3)*wap(2) )
    select case (last)
    case (.false.)
      wp3(0:1) = c0d0
      wp3(2) = - auxa*pl(2) + auxb*pl(4)
      wp3(3) =   auxa*pl(3) - auxb*pl(1)
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2:3) = c0d0
    end select

  case (36) ! v + f_- -> f_-/+, massless
  ! currents: wp1^mu, wp2(i) -> wp3(j) (wp2(0:1)=0)
  ! coupling = i*gamma_mu(l,k)*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*pslash(j,l)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp1(0) + wp1(3)
    wap(2) = wp1(0) - wp1(3)
    wap(3) = wp1(1) + cId0*wp1(2)
    wap(4) = wp1(1) - cId0*wp1(2)
    auxa = cp(1)*( - wp2(2)*wap(1) - wp2(3)*wap(4) )
    auxb = cp(1)*( - wp2(2)*wap(3) - wp2(3)*wap(2) )
    select case (last)
    case (.false.)
      wp3(0:1) = c0d0
      wp3(2) = - auxa*pl(2) + auxb*pl(4)
      wp3(3) =   auxa*pl(3) - auxb*pl(1)
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2:3) = c0d0
    end select

  case (37) ! f_+ + v -> f_+/-, massless
  ! currents: wp1(i), wp2^mu -> wp3(j) (wp1(2:3)=0)
  ! coupling = i*gamma_mu(l,k)*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*pslash(j,l)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp2(0) + wp2(3)
    wap(2) = wp2(0) - wp2(3)
    wap(3) = wp2(1) + cId0*wp2(2)
    wap(4) = wp2(1) - cId0*wp2(2)
    auxc = cp(2)*( - wp1(0)*wap(2) + wp1(1)*wap(4) )
    auxd = cp(2)*(   wp1(0)*wap(3) - wp1(1)*wap(1) )
    select case (last)
    case (.false.)
      wp3(0) = - auxc*pl(1) - auxd*pl(4)
      wp3(1) = - auxc*pl(3) - auxd*pl(2)
      wp3(2:3) = c0d0
    case (.true.)
      wp3(0:1) = c0d0
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (38) ! v + f_+ -> f_+/-, massless
  ! currents: wp1^mu, wp2(i) -> wp3(j) (wp2(2:3)=0)
  ! coupling = i*gamma_mu(l,k)*( cop(1)*omega_-(k,i) + cop(2)*omega_+(k,i) )
  ! propagator = i*pslash(j,l)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp1(0) + wp1(3)
    wap(2) = wp1(0) - wp1(3)
    wap(3) = wp1(1) + cId0*wp1(2)
    wap(4) = wp1(1) - cId0*wp1(2)
    auxc = cp(2)*( - wp2(0)*wap(2) + wp2(1)*wap(4) )
    auxd = cp(2)*(   wp2(0)*wap(3) - wp2(1)*wap(1) )
    select case (last)
    case (.false.)
      wp3(0) = - auxc*pl(1) - auxd*pl(4)
      wp3(1) = - auxc*pl(3) - auxd*pl(2)
      wp3(2:3) = c0d0
    case (.true.)
      wp3(0:1) = c0d0
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (41) ! f~ + s -> f~, massive
  ! currents: wp1(i), wp2 -> wp3(j)
  ! coupling = i*( cop(1)*omega_-(i,k) + cop(2)*omega_+(i,k) )
  ! propagator = i*(-pslash+m)(k,j)/(ps-m2); pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxa = cp(2)*wp1(0)*wp2(0)
    auxb = cp(2)*wp1(1)*wp2(0)
    auxc = cp(1)*wp1(2)*wp2(0)
    auxd = cp(1)*wp1(3)*wp2(0)
    select case (last)
    case (.false.)
      wp3(0) =   auxc*pl(2) - auxd*pl(3) + auxa*m
      wp3(1) = - auxc*pl(4) + auxd*pl(1) + auxb*m
      wp3(2) =   auxa*pl(1) + auxb*pl(3) + auxc*m
      wp3(3) =   auxa*pl(4) + auxb*pl(2) + auxd*m
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (42) ! f~ + s -> f~, massless
  ! currents: wp1(i), wp2 -> wp3(j)
  ! coupling = i*( cop(1)*omega_-(i,k) + cop(2)*omega_+(i,k) )
  ! propagator = -i*pslash(k,j)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxa = cp(2)*wp1(0)*wp2(0)
    auxb = cp(2)*wp1(1)*wp2(0)
    auxc = cp(1)*wp1(2)*wp2(0)
    auxd = cp(1)*wp1(3)*wp2(0)
    select case (last)
    case (.false.)
      wp3(0) =   auxc*pl(2) - auxd*pl(3)
      wp3(1) = - auxc*pl(4) + auxd*pl(1)
      wp3(2) =   auxa*pl(1) + auxb*pl(3)
      wp3(3) =   auxa*pl(4) + auxb*pl(2)
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (43) ! s + f~ -> f~, massive
  ! currents: wp1, wp2(i) -> wp3(j)
  ! coupling = i*( cop(1)*omega_-(i,k) + cop(2)*omega_+(i,k) )
  ! propagator = i*(-pslash+m)(k,j)/(ps-m2); pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxa = cp(2)*wp2(0)*wp1(0)
    auxb = cp(2)*wp2(1)*wp1(0)
    auxc = cp(1)*wp2(2)*wp1(0)
    auxd = cp(1)*wp2(3)*wp1(0)
    select case (last)
    case (.false.)
      wp3(0) =   auxc*pl(2) - auxd*pl(3) + auxa*m
      wp3(1) = - auxc*pl(4) + auxd*pl(1) + auxb*m
      wp3(2) =   auxa*pl(1) + auxb*pl(3) + auxc*m
      wp3(3) =   auxa*pl(4) + auxb*pl(2) + auxd*m
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (44) ! s + f~ -> f~, massless
  ! currents: wp1, wp2(i) -> wp3(j)
  ! coupling = i*( cop(1)*omega_-(i,k) + cop(2)*omega_+(i,k) )
  ! propagator = -i*pslash(k,j)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxa = cp(2)*wp2(0)*wp1(0)
    auxb = cp(2)*wp2(1)*wp1(0)
    auxc = cp(1)*wp2(2)*wp1(0)
    auxd = cp(1)*wp2(3)*wp1(0)
    select case (last)
    case (.false.)
      wp3(0) =   auxc*pl(2) - auxd*pl(3)
      wp3(1) = - auxc*pl(4) + auxd*pl(1)
      wp3(2) =   auxa*pl(1) + auxb*pl(3)
      wp3(3) =   auxa*pl(4) + auxb*pl(2)
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (45) ! f~_- + s -> f~_+/-, massless
  ! currents: wp1(i), wp2 -> wp3(j) (wp1(0:1)=0)
  ! coupling = i*( cop(1)*omega_-(i,k) + cop(2)*omega_+(i,k) )
  ! propagator = -i*pslash(k,j)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxc = cp(1)*wp1(2)*wp2(0)
    auxd = cp(1)*wp1(3)*wp2(0)
    select case (last)
    case (.false.)
      wp3(0) =   auxc*pl(2) - auxd*pl(3)
      wp3(1) = - auxc*pl(4) + auxd*pl(1)
      wp3(2:3) = c0d0
    case (.true.)
      wp3(0:1) = c0d0
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (46) ! s + f~_- -> f~_+/-, massless
  ! currents: wp1, wp2(i) -> wp3(j) (wp2(0:1)=0)
  ! coupling = i*( cop(1)*omega_-(i,k) + cop(2)*omega_+(i,k) )
  ! propagator = -i*pslash(k,j)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxc = cp(1)*wp2(2)*wp1(0)
    auxd = cp(1)*wp2(3)*wp1(0)
    select case (last)
    case (.false.)
      wp3(0) =   auxc*pl(2) - auxd*pl(3)
      wp3(1) = - auxc*pl(4) + auxd*pl(1)
      wp3(2:3) = c0d0
    case (.true.)
      wp3(0:1) = c0d0
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (47) ! f~_+ + s -> f~_-/+, massless
  ! currents: wp1(i), wp2 -> wp3(j) (wp1(2:3)=0)
  ! coupling = i*( cop(1)*omega_-(i,k) + cop(2)*omega_+(i,k) )
  ! propagator = -i*pslash(k,j)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxa = cp(2)*wp1(0)*wp2(0)
    auxb = cp(2)*wp1(1)*wp2(0)
    select case (last)
    case (.false.)
      wp3(0:1) = c0d0
      wp3(2) =   auxa*pl(1) + auxb*pl(3)
      wp3(3) =   auxa*pl(4) + auxb*pl(2)
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2:3) = c0d0
    end select

  case (48) ! s + f~_+ -> f~_-/+, massless
  ! currents: wp1, wp2(i) -> wp3(j) (wp2(2:3)=0)
  ! coupling = i*( cop(1)*omega_-(i,k) + cop(2)*omega_+(i,k) )
  ! propagator = -i*pslash(k,j)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    auxa = cp(2)*wp2(0)*wp1(0)
    auxb = cp(2)*wp2(1)*wp1(0)
    select case (last)
    case (.false.)
      wp3(0:1) = c0d0
      wp3(2) =   auxa*pl(1) + auxb*pl(3)
      wp3(3) =   auxa*pl(4) + auxb*pl(2)
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2:3) = c0d0
    end select

  case (51) ! f~ + v -> f~, massive
  ! currents: wp1(i), wp2^mu -> wp3(j)
  ! coupling = i*gamma_mu(i,k)*( cop(1)*omega_-(k,l) + cop(2)*omega_+(k,l) )
  ! propagator = i*(-pslash+m)(l,j)/(ps-m2); pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp2(0) + wp2(3)
    wap(2) = wp2(0) - wp2(3)
    wap(3) = wp2(1) + cId0*wp2(2)
    wap(4) = wp2(1) - cId0*wp2(2)
    auxa = cp(2)*( - wp1(2)*wap(2) + wp1(3)*wap(3) )
    auxb = cp(2)*(   wp1(2)*wap(4) - wp1(3)*wap(1) )
    auxc = cp(1)*( - wp1(0)*wap(1) - wp1(1)*wap(3) )
    auxd = cp(1)*( - wp1(0)*wap(4) - wp1(1)*wap(2) )
    select case (last)
    case (.false.)
      wp3(0) =   auxc*pl(2) - auxd*pl(3) + auxa*m
      wp3(1) = - auxc*pl(4) + auxd*pl(1) + auxb*m
      wp3(2) =   auxa*pl(1) + auxb*pl(3) + auxc*m
      wp3(3) =   auxa*pl(4) + auxb*pl(2) + auxd*m
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (52) ! f~ + v -> f~, massless
  ! currents: wp1(i), wp2^mu -> wp3(j)
  ! coupling = i*gamma_mu(i,k)*( cop(1)*omega_-(k,l) + cop(2)*omega_+(k,l) )
  ! propagator = -i*pslash(l,j)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp2(0) + wp2(3)
    wap(2) = wp2(0) - wp2(3)
    wap(3) = wp2(1) + cId0*wp2(2)
    wap(4) = wp2(1) - cId0*wp2(2)
    auxa = cp(2)*( - wp1(2)*wap(2) + wp1(3)*wap(3) )
    auxb = cp(2)*(   wp1(2)*wap(4) - wp1(3)*wap(1) )
    auxc = cp(1)*( - wp1(0)*wap(1) - wp1(1)*wap(3) )
    auxd = cp(1)*( - wp1(0)*wap(4) - wp1(1)*wap(2) )
    select case (last)
    case (.false.)
      wp3(0) =   auxc*pl(2) - auxd*pl(3)
      wp3(1) = - auxc*pl(4) + auxd*pl(1)
      wp3(2) =   auxa*pl(1) + auxb*pl(3)
      wp3(3) =   auxa*pl(4) + auxb*pl(2)
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (53) ! v + f~ -> f~, massive
  ! currents: wp1^mu, wp2(i) -> wp3(j)
  ! coupling = i*gamma_mu(i,k)*( cop(1)*omega_-(k,l) + cop(2)*omega_+(k,l) )
  ! propagator = i*(-pslash+m)(l,j)/(ps-m2); pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp1(0) + wp1(3)
    wap(2) = wp1(0) - wp1(3)
    wap(3) = wp1(1) + cId0*wp1(2)
    wap(4) = wp1(1) - cId0*wp1(2)
    auxa = cp(2)*( - wp2(2)*wap(2) + wp2(3)*wap(3) )
    auxb = cp(2)*(   wp2(2)*wap(4) - wp2(3)*wap(1) )
    auxc = cp(1)*( - wp2(0)*wap(1) - wp2(1)*wap(3) )
    auxd = cp(1)*( - wp2(0)*wap(4) - wp2(1)*wap(2) )
    select case (last)
    case (.false.)
      wp3(0) =   auxc*pl(2) - auxd*pl(3) + auxa*m
      wp3(1) = - auxc*pl(4) + auxd*pl(1) + auxb*m
      wp3(2) =   auxa*pl(1) + auxb*pl(3) + auxc*m
      wp3(3) =   auxa*pl(4) + auxb*pl(2) + auxd*m
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (54) ! v + f~ -> f~, massless
  ! currents: wp1^mu, wp2(i) -> wp3(j)
  ! coupling = i*gamma_mu(i,k)*( cop(1)*omega_-(k,l) + cop(2)*omega_+(k,l) )
  ! propagator = -i*pslash(l,j)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp1(0) + wp1(3)
    wap(2) = wp1(0) - wp1(3)
    wap(3) = wp1(1) + cId0*wp1(2)
    wap(4) = wp1(1) - cId0*wp1(2)
    auxa = cp(2)*( - wp2(2)*wap(2) + wp2(3)*wap(3) )
    auxb = cp(2)*(   wp2(2)*wap(4) - wp2(3)*wap(1) )
    auxc = cp(1)*( - wp2(0)*wap(1) - wp2(1)*wap(3) )
    auxd = cp(1)*( - wp2(0)*wap(4) - wp2(1)*wap(2) )
    select case (last)
    case (.false.)
      wp3(0) =   auxc*pl(2) - auxd*pl(3)
      wp3(1) = - auxc*pl(4) + auxd*pl(1)
      wp3(2) =   auxa*pl(1) + auxb*pl(3)
      wp3(3) =   auxa*pl(4) + auxb*pl(2)
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (55) ! f~_- + v -> f~_-/+, massless
  ! currents: wp1(i), wp2^mu -> wp3(j) (wp1(0:1)=0)
  ! coupling = i*gamma_mu(i,k)*( cop(1)*omega_-(k,l) + cop(2)*omega_+(k,l) )
  ! propagator = -i*pslash(l,j)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp2(0) + wp2(3)
    wap(2) = wp2(0) - wp2(3)
    wap(3) = wp2(1) + cId0*wp2(2)
    wap(4) = wp2(1) - cId0*wp2(2)
    auxa = cp(2)*( - wp1(2)*wap(2) + wp1(3)*wap(3) )
    auxb = cp(2)*(   wp1(2)*wap(4) - wp1(3)*wap(1) )
    select case (last)
    case (.false.)
      wp3(0:1) = c0d0
      wp3(2) =   auxa*pl(1) + auxb*pl(3)
      wp3(3) =   auxa*pl(4) + auxb*pl(2)
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2:3) = c0d0
    end select

  case (56) ! v + f~_- -> f~_-/+, massless
  ! currents: wp1^mu, wp2(i) -> wp3(j) (wp2(0:1)=0)
  ! coupling = i*gamma_mu(i,k)*( cop(1)*omega_-(k,l) + cop(2)*omega_+(k,l) )
  ! propagator = -i*pslash(l,j)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp1(0) + wp1(3)
    wap(2) = wp1(0) - wp1(3)
    wap(3) = wp1(1) + cId0*wp1(2)
    wap(4) = wp1(1) - cId0*wp1(2)
    auxa = cp(2)*( - wp2(2)*wap(2) + wp2(3)*wap(3) )
    auxb = cp(2)*(   wp2(2)*wap(4) - wp2(3)*wap(1) )
    select case (last)
    case (.false.)
      wp3(0:1) = c0d0
      wp3(2) =   auxa*pl(1) + auxb*pl(3)
      wp3(3) =   auxa*pl(4) + auxb*pl(2)
    case (.true.)
      wp3(0) = auxa
      wp3(1) = auxb
      wp3(2:3) = c0d0
    end select

  case (57) ! f~_+ + v -> f~_+/-, massless
  ! currents: wp1(i), wp2^mu -> wp3(j) (wp1(2:3)=0)
  ! coupling = i*gamma_mu(i,k)*( cop(1)*omega_-(k,l) + cop(2)*omega_+(k,l) )
  ! propagator = -i*pslash(l,j)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp2(0) + wp2(3)
    wap(2) = wp2(0) - wp2(3)
    wap(3) = wp2(1) + cId0*wp2(2)
    wap(4) = wp2(1) - cId0*wp2(2)
    auxc = cp(1)*( - wp1(0)*wap(1) - wp1(1)*wap(3) )
    auxd = cp(1)*( - wp1(0)*wap(4) - wp1(1)*wap(2) )
    select case (last)
    case (.false.)
      wp3(0) =   auxc*pl(2) - auxd*pl(3)
      wp3(1) = - auxc*pl(4) + auxd*pl(1)
      wp3(2:3) = c0d0
    case (.true.)
      wp3(0:1) = c0d0
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case (58) ! v + f~_+ -> f~_+/-, massless
  ! currents: wp1^mu, wp2(i) -> wp3(j) (wp2(2:3)=0)
  ! coupling = i*gamma_mu(i,k)*( cop(1)*omega_-(k,l) + cop(2)*omega_+(k,l) )
  ! propagator = -i*pslash(l,j)/ps; pslash = (p1+p2)^mu*gamma_mu

    pl = pl1 + pl2

    wap(1) = wp1(0) + wp1(3)
    wap(2) = wp1(0) - wp1(3)
    wap(3) = wp1(1) + cId0*wp1(2)
    wap(4) = wp1(1) - cId0*wp1(2)
    auxc = cp(1)*( - wp2(0)*wap(1) - wp2(1)*wap(3) )
    auxd = cp(1)*( - wp2(0)*wap(4) - wp2(1)*wap(2) )
    select case (last)
    case (.false.)
      wp3(0) =   auxc*pl(2) - auxd*pl(3)
      wp3(1) = - auxc*pl(4) + auxd*pl(1)
      wp3(2:3) = c0d0
    case (.true.)
      wp3(0:1) = c0d0
      wp3(2) = auxc
      wp3(3) = auxd
    end select

  case default

    if (warnings(306).le.warning_limit) then
      warnings(306) = warnings(306) + 1
      call openOutput
      write(nx,*)
      write(nx,*) "CODE ERROR 306 (tree_vertices_rcl.f90): wrong 3-leg interaction"
      write(nx,*)
      call toomanywarnings(306)
    endif
    call istop (ifail,2)

  end select


  end subroutine tree3

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 4 legs
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine tree4 ( last,den,cop,ty,wp1,wp2,wp3,wp4  )

  logical,     intent (in)  :: last
  integer,     intent (in)  :: ty
  complex(dp), intent (in)  :: den,cop(3),wp1(0:3),wp2(0:3),wp3(0:3)
  complex(dp), intent (out) :: wp4(0:3)

  complex(dp) :: cp(3),cc1,cc2,cc3,wpr,wpr1,wpr2,wpr3

  ! overall factor

  select case (ty)
  case (-4:-1,-14:-12)
    if (last) then
      cp(1) = cId0*cop(1)
    else
      cp(1) = - cop(1)/den
    endif
  case (-11)
    if (last) then
      cp = cId0*cop
    else
      cp = - cop/den
    endif
  end select

  ! select interaction type and contract currents

  select case (ty)

  case (-1) ! s + s + s -> s
  ! currents: wp1, wp2, wp3 -> wp4
  ! coupling = i*cop(1), propagator = i/(ps-m2)

    wp4(0) = cp(1) * wp1(0) * wp2(0) * wp3(0)
    wp4(1:3) = c0d0

  case (-2) ! s + v + v -> s
  ! currents: wp1, wp2^mu, wp3^nu -> wp4
  ! coupling = i*cop(1)*g_{mu,nu}, propagator = i/(ps-m2)

    wpr = wp2(0)*wp3(0) - sum(wp2(1:3)*wp3(1:3))
    wp4(0) = cp(1) * wpr * wp1(0)
    wp4(1:3) = c0d0

  case (-3) ! v + s + v -> s
  ! currents: wp1^mu, wp2, wp3^nu -> wp4
  ! coupling = i*cop(1)*g_{mu,nu}, propagator = i/(ps-m2)

    wpr = wp1(0)*wp3(0) - sum(wp1(1:3)*wp3(1:3))
    wp4(0) = cp(1) * wpr * wp2(0)
    wp4(1:3) = c0d0

  case (-4) ! v + v + s -> s
  ! currents: wp1^mu, wp2^nu, wp3 -> wp4
  ! coupling = i*cop(1)*g_{mu,nu}, propagator = i/(ps-m2)

    wpr = wp2(0)*wp1(0) - sum(wp2(1:3)*wp1(1:3))
    wp4(0) = cp(1) * wpr * wp3(0)
    wp4(1:3) = c0d0

  case (-11) ! v + v + v -> v
  ! currents: wp1^mu, wp2^nu, wp3^ro -> wp4^si
  ! coupling = i*(
  !            + cop(1)*g_{mu,al}*g_{nu,ro}
  !            + cop(2)*g_{mu,ro}*g_{nu,al}
  !            + cop(3)*g_{mu,nu}*g_{ro,al}
  !            )
  ! propagator = - i*g^{al,si}/(ps-m2)

    wpr1 = wp2(0)*wp3(0) - sum(wp2(1:3)*wp3(1:3))
    wpr2 = wp1(0)*wp3(0) - sum(wp1(1:3)*wp3(1:3))
    wpr3 = wp2(0)*wp1(0) - sum(wp2(1:3)*wp1(1:3))
    cc1 = - cp(1)*wpr1
    cc2 = - cp(2)*wpr2
    cc3 = - cp(3)*wpr3
    wp4 = cc1*wp1 + cc2*wp2 + cc3*wp3
    if (last) wp4(0) = - wp4(0)


  case (-12) ! s + s + v -> v
  ! currents: wp1, wp2, wp3^mu -> wp4^nu
  ! coupling = i*cop(1)*g_{mu,al}, propagator = i*g^{al,nu}/(ps-m2)

    cc3 = - cp(1) * wp1(0) * wp2(0)
    wp4 = cc3 * wp3
    if (last) wp4(0) = - wp4(0)

  case (-13) ! v + s + s -> v
  ! currents: wp1^mu, wp2, wp3 -> wp4^nu
  ! coupling = i*cop(1)*g_{mu,al}, propagator = i*g^{al,nu}/(ps-m2)

    cc1 = - cp(1) * wp2(0) * wp3(0)
    wp4 = cc1 * wp1
    if (last) wp4(0) = - wp4(0)

  case (-14) ! s + v + s -> v
  ! currents: wp1, wp2^mu, wp3 -> wp4^nu
  ! coupling = i*cop(1)*g_{mu,al}, propagator = i*g^{al,nu}/(ps-m2)

    cc2 = - cp(1) * wp1(0) * wp3(0)
    wp4 = cc2 * wp2
    if (last) wp4(0) = - wp4(0)

  case default

    if (warnings(307).le.warning_limit) then
      warnings(307) = warnings(307) + 1
      call openOutput
      write(nx,*)
      write(nx,*) "CODE ERROR 307 (tree_vertices_rcl.f90): wrong 4-leg interaction"
      write(nx,*)
      call toomanywarnings(307)
    endif
    call istop (ifail,2)

  end select


  end subroutine tree4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 end module tree_vertices_rcl

!#####################################################################


