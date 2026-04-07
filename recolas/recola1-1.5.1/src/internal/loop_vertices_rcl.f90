!#####################################################################
!!
!!  File  loop_vertices_rcl.f90
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

  module loop_vertices_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine loop3 ( riMaxIn,riMaxOut,p1,p2,pl1,pl2, &
                     m2,cop,ty,wpl,wpt,wp            )

  integer,     intent (in)    :: riMaxIn,riMaxOut,ty
  complex(dp), intent (in)    :: m2,cop(2),                     &
                                 p1(0:3),p2(0:3),pl1(4),pl2(4), &
                                 wpl(0:3,0:riMaxIn),wpt(0:3)
  complex(dp), intent (out)   :: wp(0:3,0:riMaxOut)

  integer                         :: j,riOut,mu,top,coefq
  complex(dp)                     :: cc1,cc2,cc3,cc4,cc5,spr,sprt, &
                                     auxa,auxb,auxc,auxd,m,        &
                                     wab02,wab13,wab12,wab03,      &
                                     wab20,wab31,wab21,wab30
  complex(dp), dimension(0:3)     :: p13,p23,p12,prd1,prd2,  &
                                     pp,vv1,vv2,vv3,vv4,vv5, &
                                     spg,wp1,wp2
  complex(dp), dimension(4)       :: pl,wap,wap1,wap2
  complex(dp), dimension(0:3,0:3) :: spgg,gg1

  select case (ty)
  case (21,23,31,33,41,43,51,53)
    if (CMscheme.eq.1) then
      m = sqrt(m2)
    else
      m = sqrt(real(m2,kind=dp))*c1d0
    endif
  end select

  ! select interaction type and build currents

  select case (ty)

  case (1) ! s + s -> s

    !  |
    !  | loop line in:
    !  V s = s0 + s1(mu1)*q(mu1) + s2(mu1,mu2)*q(mu1)*q(mu2) + ...
    !  | momentum (incoming) = q + p(el)
    !  |
    !  0---<-- tree line in:
    !  |       ss
    !  |       momentum (incoming) = p(et)
    !  |
    !  |
    !  V loop line out:
    !  | sss
    !  | momentum (incoming) = - q - p(e3) = - q - p(el) - p(et)
    !  |
    !
    !  coupling = - cop(1)

    cc1 = - cop(1) * wpt(0)
    wp(0,:) = cc1 * wpl(0,:)
    wp(1:3,:) = c0d0

  case (2) ! v + v -> s

    !  |
    !  | loop line in:
    !  V v(nu) = v0(nu) + v1(mu1,nu)*q(mu1) + v2(mu1,mu2,nu)*q(mu1)*q(mu2) + ...
    !  | momentum (incoming) = q + p(el)
    !  |
    !  0---<-- tree line in:
    !  |       vv(ro)
    !  |       momentum (incoming) = p(et)
    !  |
    !  |
    !  V loop line out:
    !  | s
    !  | momentum (incoming) = - q - p(e3) = - q - p(el) - p(et)
    !  |
    !
    !  coupling = - cop(1)

    wp2 = - cop(1) * wpt
    do j = 0,riMaxOut
      wp(0,j) = wpl(0,j)*wp2(0) - sum(wpl(1:3,j)*wp2(1:3))
    enddo
    wp(1:3,:) = c0d0

  case (3) ! s + v -> s

    !  |
    !  | loop line in:
    !  V s = s0 + s1(mu1)*q(mu1) + s2(mu1,mu2)*q(mu1)*q(mu2) + ...
    !  | momentum (incoming) = q + p(el)
    !  |
    !  0---<-- tree line in:
    !  |       v(nu)
    !  |       momentum (incoming) = p(et)
    !  |
    !  |
    !  V loop line out:
    !  | ss
    !  | momentum (incoming) = - q - p(e3) = - q - p(e1) - p(e2)
    !  |
    !
    !  coupling = - cop(1) * ( 2*q(nu) + pp(nu) )

    pp = + 2*p1 + p2

    spr = pp(0)*wpt(0) - sum(pp(1:3)*wpt(1:3))

    cc1 = -   cop(1) * spr
    cc2 = - 2*cop(1)
    wp2 = cc2 * wpt

    wp(0,0:riMaxIn) = cc1 * wpl(0,0:riMaxIn)
    do j = riMaxIn,0,-1
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
              wp(0,riOut) = wpl(0,j)*wp2(mu)
        else; wp(0,riOut) = wp(0,riOut) + wpl(0,j)*wp2(mu)
        endif
      enddo
    enddo
    wp(1:3,:) = c0d0

  case (4) ! v + s -> s

    !  |
    !  | loop line in:
    !  V v(nu) = v0(nu) + v1(mu1,nu)*q(mu1) + v2(mu1,mu2,nu)*q(mu1)*q(mu2) + ...
    !  | momentum (incoming) = q + p(el)
    !  |
    !  0---<-- tree line in:
    !  |       s
    !  |       momentum (incoming) = p(et)
    !  |
    !  |
    !  V loop line out:
    !  | ss
    !  | momentum (incoming) = - q - p(e3) = - q - p(e1) - p(e2)
    !  |
    !
    !  coupling = - cop(1) * ( - q(nu) + pp(nu) )

    pp = - 2*p2 - p1

    cc1 = - cop(1) * wpt(0)
    cc2 = - cc1

    do j = riMaxIn,0,-1
      spr = pp(0)*wpl(0,j) - sum(pp(1:3)*wpl(1:3,j))
      wp(0,j) = cc1 * spr
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
              wp(0,riOut) = cc2 * wpl(mu,j)
        else; wp(0,riOut) = wp(0,riOut) + cc2 * wpl(mu,j)
        endif
      enddo
    enddo

  case (5,6) ! f + {f} -> s, {f} + f -> s

    !  |
    !  | loop line in:
    !  V f(d1) = v0(d1) + v1(mu1,d1)*q(mu1) + v2(mu1,mu2,d1)*q(mu1)*q(mu2) + ...
    !  | momentum (incoming) = q + p(el)
    !  |
    !  0---<-- tree line in:
    !  |       {f}(d2)
    !  |       momentum (incoming) = p(et)
    !  |
    !  |
    !  V loop line out:
    !  | s
    !  | momentum (incoming) = - q - p(e3) = - q - p(e1) - p(e2)
    !  |

    wp2(0:1) = - cop(2) * wpt(0:1)
    wp2(2:3) = - cop(1) * wpt(2:3)
    do j = 0,riMaxOut
      wp(0,j) = sum(wp2(0:3)*wpl(0:3,j))
    enddo
    wp(1:3,:) = c0d0

  case (7,8) ! f_- + {f}_- -> s, {f}_- + f_- -> s

    wp2(2:3) = - cop(1) * wpt(2:3)
    wp(0,:) = wp2(2)*wpl(2,:) + wp2(3)*wpl(3,:)
    wp(1:3,:) = c0d0

  case (9,10) ! f_+ + {f}_+ -> s, {f}_+ + f_+ -> s

    wp2(0:1) = - cop(2) * wpt(0:1)
    wp(0,:) = wp2(0)*wpl(0,:) + wp2(1)*wpl(1,:)
    wp(1:3,:) = c0d0

  case (11) ! s + s -> v

    wp(:,riMaxIn+1:riMaxOut) = c0d0

    p12(:) = p1 - p2

    !  |
    !  | loop line in:
    !  V s = s0 + s1(mu1)*q(mu1) + s2(mu1,mu2)*q(mu1)*q(mu2) + ...
    !  | momentum (incoming) = q + p(el)
    !  |
    !  0---<-- tree line in:
    !  |       ss
    !  |       momentum (incoming) = p(et)
    !  |
    !  |
    !  V loop line out:
    !  | v(nu)
    !  | momentum (incoming) = - q - p(e3) = - q - p(e1) - p(e2)
    !  |
    !
    !  coupling =
    !    - cop(1) * ( q(nu) + p12(nu) ) =
    !    - cop(1) * ( q(nu) + ( p(e1,nu) - p(e2,nu) ) )
    !    - cop(1) * ( q(nu) + ( p(el,nu) - p(et,nu) ) )

    cc1 = - cop(1) * wpt(0)
    do j = riMaxIn,0,-1
      cc2 = cc1 * wpl(0,j)
      wp(:,j) = cc2 * p12(:)
      do mu= 0, 3
        riOut = incRI(mu,j)
        wp(mu,riOut) = wp(mu,riOut) + cc2 * gg(mu,mu)
      enddo
    enddo

  case (12) ! v + v -> v

    p12(:)=   p1 - p2
    p13(:)= 2*p1 + p2
    p23(:)= 2*p2 + p1

    !  |
    !  | loop line in:
    !  V v(mu) = v0(mu) + v1(mu1,mu)*q(mu1) + v2(mu1,mu2,mu)*q(mu1)*q(mu2) + ...
    !  | momentum (incoming) k1 = q + p(e1)
    !  |
    !  0---<-- tree line in:
    !  |       vv(nu)
    !  |       momentum (incoming) k2 = p(e2)
    !  |
    !  |
    !  V loop line out:
    !  | vvv(ro)
    !  | momentum (incoming) k3 = - q - p(e3)
    !  |
    !
    !  coupling =
    !    - (- cop(1)) * (
    !    + g(mu,nu)*( k1(ro) - k2(ro) )
    !    + g(nu,ro)*( k2(mu) - k3(mu) )
    !    + g(mu,ro)*( k3(nu) - k1(nu) )
    !    ) =
    !    - (- cop(1)) * (
    !    + g(mu,nu)*( +   q(ro) + p1(ro) - p2(ro) )
    !    + g(nu,ro)*( +   q(mu) + p2(mu) + p3(mu) )
    !    + g(mu,ro)*( - 2*q(nu) - p1(nu) - p3(nu) )
    !    )

    sprt = p13(0)*wpt(0) - sum(p13(1:3)*wpt(1:3))
    cc1  = cop(1) * sprt
    vv1  = - cop(1) * p23
    vv2  = - cop(1) * p12
    do mu = 0,3
      gg1(:,mu) = vv1 * wpt(mu) + wpt * vv2(mu)
    enddo

    cc3 = 2 * cop(1)
    vv3 = cc3 * wpt
    cc4 = - cop(1)
    vv4 = cc4 * wpt

    do j = riMaxIn,0,-1

      do mu= 0,3
        vv5(mu) = gg1(0,mu)*wpl(0,j) - sum(gg1(1:3,mu)*wpl(1:3,j))
      enddo

      spr = wpt(0)*wpl(0,j) - sum(wpt(1:3)*wpl(1:3,j))
      cc5 = cc4 * spr

      wp(:,j) = cc1*wpl(:,j) + vv5(:)
      do mu= 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
          wp(:,riOut) = + vv3(mu)*wpl(:,j) + vv4(:)*wpl(mu,j) &
                        + cc5*gg(:,mu)
        else
          wp(:,riOut) = wp(:,riOut)                           &
                        + vv3(mu)*wpl(:,j) + vv4(:)*wpl(mu,j) &
                        + cc5*gg(:,mu)
        endif
      enddo

    enddo

  case (13)  ! s + v -> v, v + s -> v

    !  |
    !  | loop line in:
    !  V s = s0 + s1(mu1)*q(mu1) + s2(mu1,mu2)*q(mu1)*q(mu2) + ...
    !  | momentum (incoming) = q + p(el)
    !  |
    !  0---<-- tree line in:
    !  |       v(nu)
    !  |       momentum (incoming) = p(et)
    !  |
    !  |
    !  V loop line out:
    !  | v(ro)
    !  | momentum (incoming) = - q - p(e3) = - q - p(e1) - p(e2)
    !  |
    !
    !  coupling = - cop(1) * g(nu,ro)

    vv1 = cop(1) * wpt
    do j = 0,riMaxIn
      wp(:,j) = vv1 * wpl(0,j)
    enddo

  case (14)  ! s + v -> v, v + s -> v

    !  |
    !  | loop line in:
    !  V v(nu) = v0(nu) + v1(mu1,nu)*q(mu1) + v2(mu1,mu2,nu)*q(mu1)*q(mu2) + ...
    !  | momentum (incoming) = q + p(el)
    !  |
    !  0---<-- tree line in:
    !  |       s
    !  |       momentum (incoming) = p(et)
    !  |
    !  |
    !  V loop line out:
    !  | v(ro)
    !  | momentum (incoming) = - q - p(e3) = - q - p(e1) - p(e2)
    !  |
    !
    !  coupling = - cop(1) * g(nu,ro)

    cc1 = cop(1) * wpt(0)
    wp(:,:) = cc1 * wpl(:,:)

  case (15) ! f + {f} -> v

    !  |
    !  | loop line in:
    !  V f(d1) = v0(d1) + v1(mu1,d1)*q(mu1) + v2(mu1,mu2,d1)*q(mu1)*q(mu2) + ...
    !  | momentum (incoming) = q + p(el)
    !  |
    !  0---<-- tree line in:
    !  |       {f}(d2)
    !  |       momentum (incoming) = p(et)
    !  |
    !  |
    !  V loop line out:
    !  | v(ro)
    !  | momentum (incoming) = - q - p(e3) = - q - p(e1) - p(e2)
    !  |

    wp2(0:1) = - cop(1) * wpt(0:1)
    wp2(2:3) = - cop(2) * wpt(2:3)
    do j = 0,riMaxOut
      wab02 = wp2(0)*wpl(2,j)
      wab13 = wp2(1)*wpl(3,j)
      wab12 = wp2(1)*wpl(2,j)
      wab03 = wp2(0)*wpl(3,j)
      wab20 = wp2(2)*wpl(0,j)
      wab31 = wp2(3)*wpl(1,j)
      wab21 = wp2(2)*wpl(1,j)
      wab30 = wp2(3)*wpl(0,j)
      prd1(0) =     wab02 + wab13
      prd1(1) =   - wab12 - wab03
      prd1(2) = ( - wab12 + wab03 )*cId0
      prd1(3) =   - wab02 + wab13
      prd2(0) =     wab20 + wab31
      prd2(1) =     wab30 + wab21
      prd2(2) = (   wab30 - wab21 )*cId0
      prd2(3) =     wab20 - wab31
      wp(:,j) = prd1 + prd2
    enddo

  case (16) ! {f} + f -> v

    wp2(0:1) = - cop(2) * wpt(0:1)
    wp2(2:3) = - cop(1) * wpt(2:3)
    do j = 0,riMaxOut
      wab02 = wpl(0,j)*wp2(2)
      wab13 = wpl(1,j)*wp2(3)
      wab12 = wpl(1,j)*wp2(2)
      wab03 = wpl(0,j)*wp2(3)
      wab20 = wpl(2,j)*wp2(0)
      wab31 = wpl(3,j)*wp2(1)
      wab21 = wpl(2,j)*wp2(1)
      wab30 = wpl(3,j)*wp2(0)
      prd1(0) =     wab02 + wab13
      prd1(1) =   - wab12 - wab03
      prd1(2) = ( - wab12 + wab03 )*cId0
      prd1(3) =   - wab02 + wab13
      prd2(0) =     wab20 + wab31
      prd2(1) =     wab30 + wab21
      prd2(2) = (   wab30 - wab21 )*cId0
      prd2(3) =     wab20 - wab31
      wp(:,j) = prd1 + prd2
    enddo

  case (17) ! f_- + {f}_+ -> v

    wp2(0:1) = - cop(1) * wpt(0:1)
    do j = 0,riMaxOut
      wab02 = wp2(0)*wpl(2,j)
      wab13 = wp2(1)*wpl(3,j)
      wab12 = wp2(1)*wpl(2,j)
      wab03 = wp2(0)*wpl(3,j)
      prd1(0) =     wab02 + wab13
      prd1(1) =   - wab12 - wab03
      prd1(2) = ( - wab12 + wab03 )*cId0
      prd1(3) =   - wab02 + wab13
      wp(:,j) = prd1
    enddo

  case (18) ! {f}_+ + f_- -> v

    wp2(2:3) = - cop(1) * wpt(2:3)
    do j = 0,riMaxOut
      wab02 = wpl(0,j)*wp2(2)
      wab13 = wpl(1,j)*wp2(3)
      wab12 = wpl(1,j)*wp2(2)
      wab03 = wpl(0,j)*wp2(3)
      prd1(0) =     wab02 + wab13
      prd1(1) =   - wab12 - wab03
      prd1(2) = ( - wab12 + wab03 )*cId0
      prd1(3) =   - wab02 + wab13
      wp(:,j) = prd1
    enddo

  case (19) ! f_+ + {f}_- -> v

    vv1(2:3) = - cop(2) * wpt(2:3)
    do j = 0,riMaxOut
      wab20 = vv1(2)*wpl(0,j)
      wab31 = vv1(3)*wpl(1,j)
      wab21 = vv1(2)*wpl(1,j)
      wab30 = vv1(3)*wpl(0,j)
      prd2(0) =   wab20 + wab31
      prd2(1) =   wab30 + wab21
      prd2(2) = ( wab30 - wab21 )*cId0
      prd2(3) =   wab20 - wab31
      wp(:,j) = prd2
    enddo

  case (20) ! {f}_- + f_+ -> v

    vv1(0:1) = - cop(2) * wpt(0:1)
    do j = 0,riMaxOut
      wab20 = wpl(2,j)*vv1(0)
      wab31 = wpl(3,j)*vv1(1)
      wab21 = wpl(2,j)*vv1(1)
      wab30 = wpl(3,j)*vv1(0)
      prd2(0) =   wab20 + wab31
      prd2(1) =   wab30 + wab21
      prd2(2) = ( wab30 - wab21 )*cId0
      prd2(3) =   wab20 - wab31
      wp(:,j) = prd2
    enddo

  case (21) ! f + s -> f, massive

    ! spg(i) = ( (slash{pl}+m)*psi*wp2 )_i
    ! spgg(i,mu) = ( gamma(mu)*psi*wp2 )_i
    ! psi = wp1
    ! auxa = (omega_-*psi)_0 = (psi)_0
    ! auxb = (omega_-*psi)_1 = (psi)_1
    ! auxc = (omega_+*psi)_2 = (psi)_2
    ! auxd = (omega_+*psi)_3 = (psi)_3

    pl = pl1 + pl2


    cc1 = - cop(1) * wpt(0)
    cc2 = - cop(2) * wpt(0)
    do j = riMaxIn,0,-1
      auxa = wpl(0,j)*cc2
      auxb = wpl(1,j)*cc2
      auxc = wpl(2,j)*cc1
      auxd = wpl(3,j)*cc1
      spgg(0,0) = - auxc
      spgg(1,0) = - auxd
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(0,2) = - auxd*cId0
      spgg(1,2) = + auxc*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,0) = - auxa
      spgg(3,0) = - auxb
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(2,2) = + auxb*cId0
      spgg(3,2) = - auxa*cId0
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0) = - auxc*pl(1) - auxd*pl(4) + auxa*m
      spg(1) = - auxc*pl(3) - auxd*pl(2) + auxb*m
      spg(2) = - auxa*pl(2) + auxb*pl(4) + auxc*m
      spg(3) = + auxa*pl(3) - auxb*pl(1) + auxd*m
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (22) ! f + s -> f, massless

    pl = pl1 + pl2

    cc1 = - cop(1) * wpt(0)
    cc2 = - cop(2) * wpt(0)
    do j = riMaxIn,0,-1
      auxa = wpl(0,j)*cc2
      auxb = wpl(1,j)*cc2
      auxc = wpl(2,j)*cc1
      auxd = wpl(3,j)*cc1
      spgg(0,0) = - auxc
      spgg(1,0) = - auxd
      spgg(2,0) = - auxa
      spgg(3,0) = - auxb
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(0,2) = - auxd*cId0
      spgg(1,2) = + auxc*cId0
      spgg(2,2) = + auxb*cId0
      spgg(3,2) = - auxa*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0)= - auxc*pl(1) - auxd*pl(4)
      spg(1)= - auxc*pl(3) - auxd*pl(2)
      spg(2)= - auxa*pl(2) + auxb*pl(4)
      spg(3)= + auxa*pl(3) - auxb*pl(1)
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (23) ! s + f -> f, massive

    ! spg(i) = ( (slash{pl}+m)*psi*wp1 )_i
    ! spgg(mu,i) = ( gamma(mu)*psi*wp1 )_i
    ! psi = wp2
    ! auxa = (omega_-*psi)_0 = (psi)_0
    ! auxb = (omega_-*psi)_1 = (psi)_1
    ! auxc = (omega_+*psi)_2 = (psi)_2
    ! auxd = (omega_+*psi)_3 = (psi)_3

    pl = pl1 + pl2


    wp2(0:1) = - cop(2) * wpt(0:1)
    wp2(2:3) = - cop(1) * wpt(2:3)
    do j = riMaxIn,0,-1
      auxa = wp2(0)*wpl(0,j)
      auxb = wp2(1)*wpl(0,j)
      auxc = wp2(2)*wpl(0,j)
      auxd = wp2(3)*wpl(0,j)
      spgg(0,0) = - auxc
      spgg(1,0) = - auxd
      spgg(2,0) = - auxa
      spgg(3,0) = - auxb
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(0,2) = - auxd*cId0
      spgg(1,2) = + auxc*cId0
      spgg(2,2) = + auxb*cId0
      spgg(3,2) = - auxa*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0)= - auxc*pl(1) - auxd*pl(4) + auxa*m
      spg(1)= - auxc*pl(3) - auxd*pl(2) + auxb*m
      spg(2)= - auxa*pl(2) + auxb*pl(4) + auxc*m
      spg(3)= + auxa*pl(3) - auxb*pl(1) + auxd*m
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (24) ! s + f -> f, massless

    pl = pl1 + pl2

    wp2(0:1) = - cop(2) * wpt(0:1)
    wp2(2:3) = - cop(1) * wpt(2:3)
    do j = riMaxIn,0,-1
      auxa = wp2(0)*wpl(0,j)
      auxb = wp2(1)*wpl(0,j)
      auxc = wp2(2)*wpl(0,j)
      auxd = wp2(3)*wpl(0,j)
      spgg(0,0) = - auxc
      spgg(1,0) = - auxd
      spgg(2,0) = - auxa
      spgg(3,0) = - auxb
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(0,2) = - auxd*cId0
      spgg(1,2) = + auxc*cId0
      spgg(2,2) = + auxb*cId0
      spgg(3,2) = - auxa*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0)= - auxc*pl(1) - auxd*pl(4)
      spg(1)= - auxc*pl(3) - auxd*pl(2)
      spg(2)= - auxa*pl(2) + auxb*pl(4)
      spg(3)= + auxa*pl(3) - auxb*pl(1)
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (25) ! f_- + s -> f_+, massless

    pl = pl1 + pl2

    cc1 = - cop(1) * wpt(0)
    do j = riMaxIn,0,-1
      auxc = wpl(2,j)*cc1
      auxd = wpl(3,j)*cc1
      spgg(0,0) = - auxc
      spgg(1,0) = - auxd
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(0,2) = - auxd*cId0
      spgg(1,2) = + auxc*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spg(0)= - auxc*pl(1) - auxd*pl(4)
      spg(1)= - auxc*pl(3) - auxd*pl(2)
      wp(0:1,j) = spg(0:1)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:1,riOut) = spgg(0:1,mu)
        else;  wp(0:1,riOut) = wp(0:1,riOut) + spgg(0:1,mu)
        endif
      enddo
    enddo
    wp(2:3,:) = c0d0

  case (26) ! s + f_- -> f_+, massless

    pl = pl1 + pl2

    wp2(2:3) = - cop(1) * wpt(2:3)
    do j = riMaxIn,0,-1
      auxc = wp2(2)*wpl(0,j)
      auxd = wp2(3)*wpl(0,j)
      spgg(0,0) = - auxc
      spgg(1,0) = - auxd
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(0,2) = - auxd*cId0
      spgg(1,2) = + auxc*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spg(0)= - auxc*pl(1) - auxd*pl(4)
      spg(1)= - auxc*pl(3) - auxd*pl(2)
      wp(0:1,j) = spg(0:1)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:1,riOut) = spgg(0:1,mu)
        else;  wp(0:1,riOut) = wp(0:1,riOut) + spgg(0:1,mu)
        endif
      enddo
    enddo
    wp(2:3,:) = c0d0

  case (27) ! f_+ + s -> f_-, massless

    pl = pl1 + pl2

    cc2 = - cop(2) * wpt(0)
    do j = riMaxIn,0,-1
      auxa = wpl(0,j)*cc2
      auxb = wpl(1,j)*cc2
      spgg(2,0) = - auxa
      spgg(3,0) = - auxb
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(2,2) = + auxb*cId0
      spgg(3,2) = - auxa*cId0
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(2)= - auxa*pl(2) + auxb*pl(4)
      spg(3)= + auxa*pl(3) - auxb*pl(1)
      wp(2:3,j) = spg(2:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(2:3,riOut) = spgg(2:3,mu)
        else;  wp(2:3,riOut) = wp(2:3,riOut) + spgg(2:3,mu)
        endif
      enddo
    enddo
    wp(0:1,:) = c0d0

  case (28) ! s + f_+ -> f_-, massless

    pl = pl1 + pl2

    wp2(0:1) = - cop(2) * wpt(0:1)
    do j = riMaxIn,0,-1
      auxa = wp2(0)*wpl(0,j)
      auxb = wp2(1)*wpl(0,j)
      spgg(2,0) = - auxa
      spgg(3,0) = - auxb
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(2,2) = + auxb*cId0
      spgg(3,2) = - auxa*cId0
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(2)= - auxa*pl(2) + auxb*pl(4)
      spg(3)= + auxa*pl(3) - auxb*pl(1)
      wp(2:3,j) = spg(2:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(2:3,riOut) = spgg(2:3,mu)
        else;  wp(2:3,riOut) = wp(2:3,riOut) + spgg(2:3,mu)
        endif
      enddo
    enddo
    wp(0:1,:) = c0d0

  case (31) ! f + v -> f, massive

    ! spg(i) = ( (slash{pl}+m)*slash{eps}*psi )_i
    ! spgg(i,mu) = ( gamma(mu)*slash{eps}*psi )_i
    ! eps = wp2
    ! psi = wp1
    ! auxa= (slash{eps}*omega_-*psi)_0
    ! auxb= (slash{eps}*omega_-*psi)_1
    ! auxc= (slash{eps}*omega_+*psi)_2
    ! auxd= (slash{eps}*omega_+*psi)_3

    pl = pl1 + pl2


    wap(1) = wpt(0) + wpt(3)
    wap(2) = wpt(0) - wpt(3)
    wap(3) = wpt(1) + cId0*wpt(2)
    wap(4) = wpt(1) - cId0*wpt(2)
    wap1 = - cop(1) * wap
    wap2 = - cop(2) * wap
    do j = riMaxIn,0,-1
      wp1(0:1) = wpl(0:1,j)
      wp1(2:3) = wpl(2:3,j)
      auxa  = - wp1(2)*wap1(1) - wp1(3)*wap1(4)
      auxb  = - wp1(2)*wap1(3) - wp1(3)*wap1(2)
      auxc  = - wp1(0)*wap2(2) + wp1(1)*wap2(4)
      auxd  = + wp1(0)*wap2(3) - wp1(1)*wap2(1)
      spgg(0,0) = - auxc
      spgg(1,0) = - auxd
      spgg(2,0) = - auxa
      spgg(3,0) = - auxb
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(0,2) = - auxd*cId0
      spgg(1,2) = + auxc*cId0
      spgg(2,2) = + auxb*cId0
      spgg(3,2) = - auxa*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0)= - auxc*pl(1) - auxd*pl(4) + auxa*m
      spg(1)= - auxc*pl(3) - auxd*pl(2) + auxb*m
      spg(2)= - auxa*pl(2) + auxb*pl(4) + auxc*m
      spg(3)= + auxa*pl(3) - auxb*pl(1) + auxd*m
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (32) ! f + v -> f, massless

    pl = pl1 + pl2

    wap(1) = wpt(0) + wpt(3)
    wap(2) = wpt(0) - wpt(3)
    wap(3) = wpt(1) + cId0*wpt(2)
    wap(4) = wpt(1) - cId0*wpt(2)
    wap1 = - cop(1) * wap
    wap2 = - cop(2) * wap
    do j = riMaxIn,0,-1
      wp1(0:1) = wpl(0:1,j)
      wp1(2:3) = wpl(2:3,j)
      auxa  = - wp1(2)*wap1(1) - wp1(3)*wap1(4)
      auxb  = - wp1(2)*wap1(3) - wp1(3)*wap1(2)
      auxc  = - wp1(0)*wap2(2) + wp1(1)*wap2(4)
      auxd  = + wp1(0)*wap2(3) - wp1(1)*wap2(1)
      spgg(0,0) = - auxc
      spgg(1,0) = - auxd
      spgg(2,0) = - auxa
      spgg(3,0) = - auxb
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(0,2) = - auxd*cId0
      spgg(1,2) = + auxc*cId0
      spgg(2,2) = + auxb*cId0
      spgg(3,2) = - auxa*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0)= - auxc*pl(1) - auxd*pl(4)
      spg(1)= - auxc*pl(3) - auxd*pl(2)
      spg(2)= - auxa*pl(2) + auxb*pl(4)
      spg(3)= + auxa*pl(3) - auxb*pl(1)
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (33) ! v + f -> f, massive

    ! spg(i) = ( (slash{pl}+m)*slash{eps}*psi )_i
    ! spgg(i,mu) = ( gamma(mu)*slash{eps}*psi )_i
    ! eps = wp1
    ! psi = wp2
    ! auxa= (slash{eps}*omega_-*psi)_0
    ! auxb= (slash{eps}*omega_-*psi)_1
    ! auxc= (slash{eps}*omega_+*psi)_2
    ! auxd= (slash{eps}*omega_+*psi)_3

    pl = pl1 + pl2


    wp2(0:1) = - cop(2) * wpt(0:1)
    wp2(2:3) = - cop(1) * wpt(2:3)
    do j = riMaxIn,0,-1
      wp1 = wpl(:,j)
      wap(1) = wp1(0) + wp1(3)
      wap(2) = wp1(0) - wp1(3)
      wap(3) = wp1(1) + cId0*wp1(2)
      wap(4) = wp1(1) - cId0*wp1(2)
      auxa  = - wp2(2)*wap(1) - wp2(3)*wap(4)
      auxb  = - wp2(2)*wap(3) - wp2(3)*wap(2)
      auxc  = - wp2(0)*wap(2) + wp2(1)*wap(4)
      auxd  = + wp2(0)*wap(3) - wp2(1)*wap(1)
      spgg(0,0) = - auxc
      spgg(1,0) = - auxd
      spgg(2,0) = - auxa
      spgg(3,0) = - auxb
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(0,2) = - auxd*cId0
      spgg(1,2) = + auxc*cId0
      spgg(2,2) = + auxb*cId0
      spgg(3,2) = - auxa*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0)= - auxc*pl(1) - auxd*pl(4) + auxa*m
      spg(1)= - auxc*pl(3) - auxd*pl(2) + auxb*m
      spg(2)= - auxa*pl(2) + auxb*pl(4) + auxc*m
      spg(3)= + auxa*pl(3) - auxb*pl(1) + auxd*m
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (34) ! v + f -> f, massless

    pl = pl1 + pl2

    wp2(0:1) = - cop(2) * wpt(0:1)
    wp2(2:3) = - cop(1) * wpt(2:3)
    do j = riMaxIn,0,-1
      wp1 = wpl(:,j)
      wap(1) = wp1(0) + wp1(3)
      wap(2) = wp1(0) - wp1(3)
      wap(3) = wp1(1) + cId0*wp1(2)
      wap(4) = wp1(1) - cId0*wp1(2)
      auxa  = - wp2(2)*wap(1) - wp2(3)*wap(4)
      auxb  = - wp2(2)*wap(3) - wp2(3)*wap(2)
      auxc  = - wp2(0)*wap(2) + wp2(1)*wap(4)
      auxd  = + wp2(0)*wap(3) - wp2(1)*wap(1)
      spgg(2,0) = - auxa
      spgg(3,0) = - auxb
      spgg(0,0) = - auxc
      spgg(1,0) = - auxd
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(0,2) = - auxd*cId0
      spgg(1,2) = + auxc*cId0
      spgg(2,2) = + auxb*cId0
      spgg(3,2) = - auxa*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0)= - auxc*pl(1) - auxd*pl(4)
      spg(1)= - auxc*pl(3) - auxd*pl(2)
      spg(2)= - auxa*pl(2) + auxb*pl(4)
      spg(3)= + auxa*pl(3) - auxb*pl(1)
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (35) ! f_- + v -> f_-, massless

    pl = pl1 + pl2

    wap(1) = wpt(0) + wpt(3)
    wap(2) = wpt(0) - wpt(3)
    wap(3) = wpt(1) + cId0*wpt(2)
    wap(4) = wpt(1) - cId0*wpt(2)
    wap = - cop(1) * wap
    do j = riMaxIn,0,-1
      wp1(2:3) = wpl(2:3,j)
      auxa  = - wp1(2)*wap(1) - wp1(3)*wap(4)
      auxb  = - wp1(2)*wap(3) - wp1(3)*wap(2)
      spgg(2,0) = - auxa
      spgg(3,0) = - auxb
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(2,2) = + auxb*cId0
      spgg(3,2) = - auxa*cId0
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(2)= - auxa*pl(2) + auxb*pl(4)
      spg(3)= + auxa*pl(3) - auxb*pl(1)
      wp(2:3,j) = spg(2:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(2:3,riOut) = spgg(2:3,mu)
        else;  wp(2:3,riOut) = wp(2:3,riOut) + spgg(2:3,mu)
        endif
      enddo
    enddo
    wp(0:1,:) = c0d0

  case (36) ! v + f_- -> f_-, massless

    pl = pl1 + pl2

    wp2(2:3) = - cop(1) * wpt(2:3)
    do j = riMaxIn,0,-1
      wp1      = wpl(:,j)
      wap(1) = wp1(0) + wp1(3)
      wap(2) = wp1(0) - wp1(3)
      wap(3) = wp1(1) + cId0*wp1(2)
      wap(4) = wp1(1) - cId0*wp1(2)
      auxa  = - wp2(2)*wap(1) - wp2(3)*wap(4)
      auxb  = - wp2(2)*wap(3) - wp2(3)*wap(2)
      spgg(2,0) = - auxa
      spgg(3,0) = - auxb
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(2,2) = + auxb*cId0
      spgg(3,2) = - auxa*cId0
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(2)= - auxa*pl(2) + auxb*pl(4)
      spg(3)= + auxa*pl(3) - auxb*pl(1)
      wp(2:3,j) = spg(2:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(2:3,riOut) = spgg(2:3,mu)
        else;  wp(2:3,riOut) = wp(2:3,riOut) + spgg(2:3,mu)
        endif
      enddo
    enddo
    wp(0:1,:) = c0d0

  case (37) ! f_+ + v -> f_+, massless

    pl = pl1 + pl2

    wap(1) = wpt(0) + wpt(3)
    wap(2) = wpt(0) - wpt(3)
    wap(3) = wpt(1) + cId0*wpt(2)
    wap(4) = wpt(1) - cId0*wpt(2)
    wap = - cop(2) * wap
    do j = riMaxIn,0,-1
      wp1(0:1) = wpl(0:1,j)
      auxc  = - wp1(0)*wap(2) + wp1(1)*wap(4)
      auxd  = + wp1(0)*wap(3) - wp1(1)*wap(1)
      spgg(0,0) = - auxc
      spgg(1,0) = - auxd
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(0,2) = - auxd*cId0
      spgg(1,2) = + auxc*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spg(0)= - auxc*pl(1) - auxd*pl(4)
      spg(1)= - auxc*pl(3) - auxd*pl(2)
      wp(0:1,j) = spg(0:1)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:1,riOut) = spgg(0:1,mu)
        else;  wp(0:1,riOut) = wp(0:1,riOut) + spgg(0:1,mu)
        endif
      enddo
    enddo
    wp(2:3,:) = c0d0

  case (38) ! v + f_+ -> f_+, massless

    pl = pl1 + pl2

    wp2(0:1) = - cop(2) * wpt(0:1)
    do j = riMaxIn,0,-1
      wp1      = wpl(:,j)
      wap(1) = wp1(0) + wp1(3)
      wap(2) = wp1(0) - wp1(3)
      wap(3) = wp1(1) + cId0*wp1(2)
      wap(4) = wp1(1) - cId0*wp1(2)
      auxc  = - wp2(0)*wap(2) + wp2(1)*wap(4)
      auxd  = + wp2(0)*wap(3) - wp2(1)*wap(1)
      spgg(0,0) = - auxc
      spgg(1,0) = - auxd
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(0,2) = - auxd*cId0
      spgg(1,2) = + auxc*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spg(0)= - auxc*pl(1) - auxd*pl(4)
      spg(1)= - auxc*pl(3) - auxd*pl(2)
      wp(0:1,j) = spg(0:1)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:1,riOut) = spgg(0:1,mu)
        else;  wp(0:1,riOut) = wp(0:1,riOut) + spgg(0:1,mu)
        endif
      enddo
    enddo
    wp(2:3,:) = c0d0

  case (41) ! {f} + s -> {f}, massive

    ! spg(i) = ( psibar*(-slash{pl}+m)*wp2 )_i
    ! spgg(mu,i) = ( psibar*gamma(mu)*wp2 )_i
    ! psi = wp1
    ! auxa = (psibar*omega_+)_0 = (psibar)_0
    ! auxb = (psibar*omega_+)_1 = (psibar)_1
    ! auxc = (psibar*omega_-)_2 = (psibar)_2
    ! auxd = (psibar*omega_-)_3 = (psibar)_3

    pl = pl1 + pl2


    cc1 = - cop(1) * wpt(0)
    cc2 = - cop(2) * wpt(0)
    do j = riMaxIn,0,-1
      auxa = wpl(0,j)*cc2
      auxb = wpl(1,j)*cc2
      auxc = wpl(2,j)*cc1
      auxd = wpl(3,j)*cc1
      spgg(0,0) = + auxc
      spgg(1,0) = + auxd
      spgg(2,0) = + auxa
      spgg(3,0) = + auxb
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(0,2) = + auxd*cId0
      spgg(1,2) = - auxc*cId0
      spgg(2,2) = - auxb*cId0
      spgg(3,2) = + auxa*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0)= + auxc*pl(2) - auxd*pl(3) + auxa*m
      spg(1)= - auxc*pl(4) + auxd*pl(1) + auxb*m
      spg(2)= + auxa*pl(1) + auxb*pl(3) + auxc*m
      spg(3)= + auxa*pl(4) + auxb*pl(2) + auxd*m
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (42) ! {f} + s -> {f}, massless

    pl = pl1 + pl2

    cc1 = - cop(1) * wpt(0)
    cc2 = - cop(2) * wpt(0)
    do j = riMaxIn,0,-1
      auxa = wpl(0,j)*cc2
      auxb = wpl(1,j)*cc2
      auxc = wpl(2,j)*cc1
      auxd = wpl(3,j)*cc1
      spgg(0,0) = + auxc
      spgg(1,0) = + auxd
      spgg(2,0) = + auxa
      spgg(3,0) = + auxb
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(0,2) = + auxd*cId0
      spgg(1,2) = - auxc*cId0
      spgg(2,2) = - auxb*cId0
      spgg(3,2) = + auxa*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0)= + auxc*pl(2) - auxd*pl(3)
      spg(1)= - auxc*pl(4) + auxd*pl(1)
      spg(2)= + auxa*pl(1) + auxb*pl(3)
      spg(3)= + auxa*pl(4) + auxb*pl(2)
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (43) ! s + {f} -> {f}, massive

    ! spg(i) = ( psibar*(-slash{pl}+m)*wp1 )_i
    ! spgg(mu,i) = ( psibar*gamma(mu)*wp1 )_i
    ! psi = wp2
    ! auxa = (psibar*omega_+)_0 = (psibar)_0
    ! auxb = (psibar*omega_+)_1 = (psibar)_1
    ! auxc = (psibar*omega_-)_2 = (psibar)_2
    ! auxd = (psibar*omega_-)_3 = (psibar)_3

    pl = pl1 + pl2


    wp2(0:1) = - cop(2) * wpt(0:1)
    wp2(2:3) = - cop(1) * wpt(2:3)
    do j = riMaxIn,0,-1
      auxa = wp2(0)*wpl(0,j)
      auxb = wp2(1)*wpl(0,j)
      auxc = wp2(2)*wpl(0,j)
      auxd = wp2(3)*wpl(0,j)
      spgg(0,0) = + auxc
      spgg(1,0) = + auxd
      spgg(2,0) = + auxa
      spgg(3,0) = + auxb
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(0,2) = + auxd*cId0
      spgg(1,2) = - auxc*cId0
      spgg(2,2) = - auxb*cId0
      spgg(3,2) = + auxa*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0)= + auxc*pl(2) - auxd*pl(3) + auxa*m
      spg(1)= - auxc*pl(4) + auxd*pl(1) + auxb*m
      spg(2)= + auxa*pl(1) + auxb*pl(3) + auxc*m
      spg(3)= + auxa*pl(4) + auxb*pl(2) + auxd*m
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (44) ! s + {f} -> {f}, massless

    pl = pl1 + pl2

    wp2(0:1) = - cop(2) * wpt(0:1)
    wp2(2:3) = - cop(1) * wpt(2:3)
    do j = riMaxIn,0,-1
      auxa = wp2(0)*wpl(0,j)
      auxb = wp2(1)*wpl(0,j)
      auxc = wp2(2)*wpl(0,j)
      auxd = wp2(3)*wpl(0,j)
      spgg(0,0) = + auxc
      spgg(1,0) = + auxd
      spgg(2,0) = + auxa
      spgg(3,0) = + auxb
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(0,2) = + auxd*cId0
      spgg(1,2) = - auxc*cId0
      spgg(2,2) = - auxb*cId0
      spgg(3,2) = + auxa*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0)= + auxc*pl(2) - auxd*pl(3)
      spg(1)= - auxc*pl(4) + auxd*pl(1)
      spg(2)= + auxa*pl(1) + auxb*pl(3)
      spg(3)= + auxa*pl(4) + auxb*pl(2)
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (45) ! {f}_- + s -> {f}_+, massless

    pl = pl1 + pl2

    cc1 = - cop(1) * wpt(0)
    do j = riMaxIn,0,-1
      auxc = wpl(2,j)*cc1
      auxd = wpl(3,j)*cc1
      spgg(0,0) = + auxc
      spgg(1,0) = + auxd
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(0,2) = + auxd*cId0
      spgg(1,2) = - auxc*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spg(0)= + auxc*pl(2) - auxd*pl(3)
      spg(1)= - auxc*pl(4) + auxd*pl(1)
      wp(0:1,j) = spg(0:1)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:1,riOut) = spgg(0:1,mu)
        else;  wp(0:1,riOut) = wp(0:1,riOut) + spgg(0:1,mu)
        endif
      enddo
    enddo
    wp(2:3,:) = c0d0

  case (46) ! s + {f}_- -> {f}_+, massless

    pl = pl1 + pl2

    wp2(2:3) = - cop(1) * wpt(2:3)
    do j = riMaxIn,0,-1
      auxc = wp2(2)*wpl(0,j)
      auxd = wp2(3)*wpl(0,j)
      spgg(0,0) = + auxc
      spgg(1,0) = + auxd
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(0,2) = + auxd*cId0
      spgg(1,2) = - auxc*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spg(0)= + auxc*pl(2) - auxd*pl(3)
      spg(1)= - auxc*pl(4) + auxd*pl(1)
      wp(0:1,j) = spg(0:1)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:1,riOut) = spgg(0:1,mu)
        else;  wp(0:1,riOut) = wp(0:1,riOut) + spgg(0:1,mu)
        endif
      enddo
    enddo
    wp(2:3,:) = c0d0

  case (47) ! {f}_+ + s -> {f}_-, massless

    pl = pl1 + pl2

    cc2 = - cop(2) * wpt(0)
    do j = riMaxIn,0,-1
      auxa = wpl(0,j)*cc2
      auxb = wpl(1,j)*cc2
      spgg(2,0) = + auxa
      spgg(3,0) = + auxb
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(2,2) = - auxb*cId0
      spgg(3,2) = + auxa*cId0
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(2)= + auxa*pl(1) + auxb*pl(3)
      spg(3)= + auxa*pl(4) + auxb*pl(2)
      wp(2:3,j) = spg(2:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(2:3,riOut) = spgg(2:3,mu)
        else;  wp(2:3,riOut) = wp(2:3,riOut) + spgg(2:3,mu)
        endif
      enddo
    enddo
    wp(0:1,:) = c0d0

  case (48) ! s + {f}_+ -> {f}_-, massless

    pl = pl1 + pl2

    wp2(0:1) = - cop(2) * wpt(0:1)
    do j = riMaxIn,0,-1
      auxa = wp2(0)*wpl(0,j)
      auxb = wp2(1)*wpl(0,j)
      spgg(2,0) = + auxa
      spgg(3,0) = + auxb
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(2,2) = - auxb*cId0
      spgg(3,2) = + auxa*cId0
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(2)= + auxa*pl(1) + auxb*pl(3)
      spg(3)= + auxa*pl(4) + auxb*pl(2)
      wp(2:3,j) = spg(2:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(2:3,riOut) = spgg(2:3,mu)
        else;  wp(2:3,riOut) = wp(2:3,riOut) + spgg(2:3,mu)
        endif
      enddo
    enddo
    wp(0:1,:) = c0d0

  case (51) ! {f} + v -> {f}, massive

    ! spg(i) = ( psibar*slash{eps}*(-slash{pl}+m) )_i
    ! spgg(i,mu) = ( psibar*slash{eps}*(-gamma(mu)) )_i
    ! eps = wp2
    ! bar{psi} = wp1
    ! auxa= (psibar*slash{eps}*omega_+)_0
    ! auxb= (psibar*slash{eps}*omega_+)_1
    ! auxc= (psibar*slash{eps}*omega_-)_2
    ! auxd= (psibar*slash{eps}*omega_-)_3

    pl = pl1 + pl2


    wap(1) = wpt(0) + wpt(3)
    wap(2) = wpt(0) - wpt(3)
    wap(3) = wpt(1) + cId0*wpt(2)
    wap(4) = wpt(1) - cId0*wpt(2)
    wap1 = - cop(1) * wap
    wap2 = - cop(2) * wap
    do j = riMaxIn,0,-1
      wp1 = wpl(:,j)
      auxa  = - wp1(2)*wap2(2) + wp1(3)*wap2(3)
      auxb  = + wp1(2)*wap2(4) - wp1(3)*wap2(1)
      auxc  = - wp1(0)*wap1(1) - wp1(1)*wap1(3)
      auxd  = - wp1(0)*wap1(4) - wp1(1)*wap1(2)
      spgg(0,0) = + auxc
      spgg(1,0) = + auxd
      spgg(2,0) = + auxa
      spgg(3,0) = + auxb
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(0,2) = + auxd*cId0
      spgg(1,2) = - auxc*cId0
      spgg(2,2) = - auxb*cId0
      spgg(3,2) = + auxa*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0)= + auxc*pl(2) - auxd*pl(3) + auxa*m
      spg(1)= - auxc*pl(4) + auxd*pl(1) + auxb*m
      spg(2)= + auxa*pl(1) + auxb*pl(3) + auxc*m
      spg(3)= + auxa*pl(4) + auxb*pl(2) + auxd*m
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (52) ! {f} + v -> {f}, massless

    pl = pl1 + pl2

    wap(1) = wpt(0) + wpt(3)
    wap(2) = wpt(0) - wpt(3)
    wap(3) = wpt(1) + cId0*wpt(2)
    wap(4) = wpt(1) - cId0*wpt(2)
    wap1 = - cop(1) * wap
    wap2 = - cop(2) * wap
    do j = riMaxIn,0,-1
      wp1 = wpl(:,j)
      auxa  = - wp1(2)*wap2(2) + wp1(3)*wap2(3)
      auxb  = + wp1(2)*wap2(4) - wp1(3)*wap2(1)
      auxc  = - wp1(0)*wap1(1) - wp1(1)*wap1(3)
      auxd  = - wp1(0)*wap1(4) - wp1(1)*wap1(2)
      spgg(0,0) = + auxc
      spgg(1,0) = + auxd
      spgg(2,0) = + auxa
      spgg(3,0) = + auxb
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(0,2) = + auxd*cId0
      spgg(1,2) = - auxc*cId0
      spgg(2,2) = - auxb*cId0
      spgg(3,2) = + auxa*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0)= + auxc*pl(2) - auxd*pl(3)
      spg(1)= - auxc*pl(4) + auxd*pl(1)
      spg(2)= + auxa*pl(1) + auxb*pl(3)
      spg(3)= + auxa*pl(4) + auxb*pl(2)
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (53) ! v + {f} -> {f}, massive

    ! spg(i) = ( psibar*slash{eps}*(-slash{pl}+m) )_i
    ! spgg(i,mu) = ( psibar*slash{eps}*(-gamma(mu)) )_i
    ! eps = wp1
    ! bar{psi} = wp2
    ! auxa= (psibar*slash{eps}*omega_+)_0
    ! auxb= (psibar*slash{eps}*omega_+)_1
    ! auxc= (psibar*slash{eps}*omega_-)_2
    ! auxd= (psibar*slash{eps}*omega_-)_3

    pl = pl1 + pl2


    wp2(0:1) = - cop(1) * wpt(0:1)
    wp2(2:3) = - cop(2) * wpt(2:3)
    do j = riMaxIn,0,-1
      wp1 = wpl(:,j)
      wap(1) = wp1(0) + wp1(3)
      wap(2) = wp1(0) - wp1(3)
      wap(3) = wp1(1) + cId0*wp1(2)
      wap(4) = wp1(1) - cId0*wp1(2)
      auxa  = - wp2(2)*wap(2) + wp2(3)*wap(3)
      auxb  = + wp2(2)*wap(4) - wp2(3)*wap(1)
      auxc  = - wp2(0)*wap(1) - wp2(1)*wap(3)
      auxd  = - wp2(0)*wap(4) - wp2(1)*wap(2)
      spgg(0,0) = + auxc
      spgg(1,0) = + auxd
      spgg(2,0) = + auxa
      spgg(3,0) = + auxb
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(0,2) = + auxd*cId0
      spgg(1,2) = - auxc*cId0
      spgg(2,2) = - auxb*cId0
      spgg(3,2) = + auxa*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0)= + auxc*pl(2) - auxd*pl(3) + auxa*m
      spg(1)= - auxc*pl(4) + auxd*pl(1) + auxb*m
      spg(2)= + auxa*pl(1) + auxb*pl(3) + auxc*m
      spg(3)= + auxa*pl(4) + auxb*pl(2) + auxd*m
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (54) ! v + {f} -> {f}, massless

    pl = pl1 + pl2

    wp2(0:1) = - cop(1) * wpt(0:1)
    wp2(2:3) = - cop(2) * wpt(2:3)
    do j = riMaxIn,0,-1
      wp1 = wpl(:,j)
      wap(1) = wp1(0) + wp1(3)
      wap(2) = wp1(0) - wp1(3)
      wap(3) = wp1(1) + cId0*wp1(2)
      wap(4) = wp1(1) - cId0*wp1(2)
      auxa  = - wp2(2)*wap(2) + wp2(3)*wap(3)
      auxb  = + wp2(2)*wap(4) - wp2(3)*wap(1)
      auxc  = - wp2(0)*wap(1) - wp2(1)*wap(3)
      auxd  = - wp2(0)*wap(4) - wp2(1)*wap(2)
      spgg(0,0) = + auxc
      spgg(1,0) = + auxd
      spgg(2,0) = + auxa
      spgg(3,0) = + auxb
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(0,2) = + auxd*cId0
      spgg(1,2) = - auxc*cId0
      spgg(2,2) = - auxb*cId0
      spgg(3,2) = + auxa*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(0)= + auxc*pl(2) - auxd*pl(3)
      spg(1)= - auxc*pl(4) + auxd*pl(1)
      spg(2)= + auxa*pl(1) + auxb*pl(3)
      spg(3)= + auxa*pl(4) + auxb*pl(2)
      wp(0:3,j) = spg(0:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:3,riOut) = spgg(0:3,mu)
        else;  wp(0:3,riOut) = wp(0:3,riOut) + spgg(0:3,mu)
        endif
      enddo
    enddo

  case (55) ! {f}_- + v -> {f}_-, massless

    pl = pl1 + pl2

    wap(1) = wpt(0) + wpt(3)
    wap(2) = wpt(0) - wpt(3)
    wap(3) = wpt(1) + cId0*wpt(2)
    wap(4) = wpt(1) - cId0*wpt(2)
    wap = - cop(2) * wap
    do j = riMaxIn,0,-1
      wp1(2:3) = wpl(2:3,j)
      auxa  = - wp1(2)*wap(2) + wp1(3)*wap(3)
      auxb  = + wp1(2)*wap(4) - wp1(3)*wap(1)
      spgg(2,0) = + auxa
      spgg(3,0) = + auxb
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(2,2) = - auxb*cId0
      spgg(3,2) = + auxa*cId0
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(2)= + auxa*pl(1) + auxb*pl(3)
      spg(3)= + auxa*pl(4) + auxb*pl(2)
      wp(2:3,j) = spg(2:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(2:3,riOut) = spgg(2:3,mu)
        else;  wp(2:3,riOut) = wp(2:3,riOut) + spgg(2:3,mu)
        endif
      enddo
    enddo
    wp(0:1,:) = c0d0

  case (56) ! v + {f}_- -> {f}_-, massless

    pl = pl1 + pl2

    wp2(2:3) = - cop(2) * wpt(2:3)
    do j = riMaxIn,0,-1
      wp1      = wpl(:,j)
      wap(1) = wp1(0) + wp1(3)
      wap(2) = wp1(0) - wp1(3)
      wap(3) = wp1(1) + cId0*wp1(2)
      wap(4) = wp1(1) - cId0*wp1(2)
      auxa  = - wp2(2)*wap(2) + wp2(3)*wap(3)
      auxb  = + wp2(2)*wap(4) - wp2(3)*wap(1)
      spgg(2,0) = + auxa
      spgg(3,0) = + auxb
      spgg(2,1) = - auxb
      spgg(3,1) = - auxa
      spgg(2,2) = - auxb*cId0
      spgg(3,2) = + auxa*cId0
      spgg(2,3) = - auxa
      spgg(3,3) = + auxb
      spg(2)= + auxa*pl(1) + auxb*pl(3)
      spg(3)= + auxa*pl(4) + auxb*pl(2)
      wp(2:3,j) = spg(2:3)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(2:3,riOut) = spgg(2:3,mu)
        else;  wp(2:3,riOut) = wp(2:3,riOut) + spgg(2:3,mu)
        endif
      enddo
    enddo
    wp(0:1,:) = c0d0

  case (57) ! {f}_+ + v -> {f}_+, massless

    pl = pl1 + pl2

    wap(1) = wpt(0) + wpt(3)
    wap(2) = wpt(0) - wpt(3)
    wap(3) = wpt(1) + cId0*wpt(2)
    wap(4) = wpt(1) - cId0*wpt(2)
    wap = - cop(1) * wap
    do j = riMaxIn,0,-1
      wp1(0:1) = wpl(0:1,j)
      auxc  = - wp1(0)*wap(1) - wp1(1)*wap(3)
      auxd  = - wp1(0)*wap(4) - wp1(1)*wap(2)
      spgg(0,0) = + auxc
      spgg(1,0) = + auxd
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(0,2) = + auxd*cId0
      spgg(1,2) = - auxc*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spg(0)= + auxc*pl(2) - auxd*pl(3)
      spg(1)= - auxc*pl(4) + auxd*pl(1)
      wp(0:1,j) = spg(0:1)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:1,riOut) = spgg(0:1,mu)
        else;  wp(0:1,riOut) = wp(0:1,riOut) + spgg(0:1,mu)
        endif
      enddo
    enddo
    wp(2:3,:) = c0d0

  case (58) ! v + {f}_+ -> {f}_+, massless

    pl = pl1 + pl2

    wp2(0:1) = - cop(1) * wpt(0:1)
    do j = riMaxIn,0,-1
      wp1      = wpl(:,j)
      wap(1) = wp1(0) + wp1(3)
      wap(2) = wp1(0) - wp1(3)
      wap(3) = wp1(1) + cId0*wp1(2)
      wap(4) = wp1(1) - cId0*wp1(2)
      auxc  = - wp2(0)*wap(1) - wp2(1)*wap(3)
      auxd  = - wp2(0)*wap(4) - wp2(1)*wap(2)
      spgg(0,0) = + auxc
      spgg(1,0) = + auxd
      spgg(0,1) = + auxd
      spgg(1,1) = + auxc
      spgg(0,2) = + auxd*cId0
      spgg(1,2) = - auxc*cId0
      spgg(0,3) = + auxc
      spgg(1,3) = - auxd
      spg(0)= + auxc*pl(2) - auxd*pl(3)
      spg(1)= - auxc*pl(4) + auxd*pl(1)
      wp(0:1,j) = spg(0:1)
      do mu = 0,3
        riOut = incRI(mu,j)
        if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
               wp(0:1,riOut) = spgg(0:1,mu)
        else;  wp(0:1,riOut) = wp(0:1,riOut) + spgg(0:1,mu)
        endif
      enddo
    enddo
    wp(2:3,:) = c0d0

  case (61,71) ! x + s -> x, s + x -> x, {x} + s -> {x}, s + {x} -> {x}

    cc1 = - cop(1) * wpt(0)
    wp(0,:) = cc1 * wpl(0,:)
    wp(1:3,:) = c0d0

  case (62,63,72,73) ! x + v -> x, v + x -> x, {x} + v -> {x}, v + {x} -> {x}

    select case (ty)
    case (62)
      coefq = - 1
      pp = - p1 - p2
      top = 1
    case (63)
      coefq = - 1
      pp = - p1 - p2
      top = 2
    case (72)
      coefq = + 1
      pp = + p1
      top = 1
    case (73)
      coefq =   0
      pp = + p2
      top = 2
    end select

    select case (top)
    case (1)
      !
      !  |
      !  | loop line in:
      !  V s = s0 + s1(mu1)*q(mu1) + s2(mu1,mu2)*q(mu1)*q(mu2) + ...
      !  | momentum (incoming) = q + p(el)
      !  |
      !  0---<-- tree line in:
      !  |       v(nu)
      !  |       momentum (incoming) = p(et)
      !  |
      !  |
      !  V loop line out:
      !  | ss
      !  | momentum (incoming) = - q - p(e3)
      !  |
      !
      !  coupling = - cop(1) * ( coefq*q(nu) + pp(nu) )
      !
      spr = pp(0)*wpt(0) - sum(pp(1:3)*wpt(1:3))
      cc1 = - cop(1) * spr
      cc2 = - cop(1) * coefq
      vv2 = cc2 * wpt
      do j = riMaxIn,0,-1
        wp(0,j) = cc1 * wpl(0,j)
        if ( coefq .ne. 0 ) then
          do mu = 0,3
            riOut = incRI(mu,j)
            if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
                   wp(0,riOut) = wpl(0,j)*vv2(mu)
            else;  wp(0,riOut) = wp(0,riOut) + wpl(0,j)*vv2(mu)
            endif
          enddo
        endif
      enddo
    case (2)
      !
      !  |
      !  | loop line in:
      !  V v(nu) = v0(nu) + v1(mu1,nu)*q(mu1) + v2(mu1,mu2,nu)*q(mu1)*q(mu2) + ...
      !  | momentum (incoming) = q + p(el)
      !  |
      !  0---<-- tree line in:
      !  |       s
      !  |       momentum (incoming) = p(et)
      !  |
      !  |
      !  V loop line out:
      !  | ss
      !  | momentum (incoming) = - q - p(e3)
      !  |
      !
      !  coupling = sign * - cop(1) * ( coefq*q(nu) + pp(nu) )
      !
      cc1 = - cop(1) * wpt(0)
      vv1 = cc1 * pp
      cc2 = cc1 * coefq
      do j = riMaxIn,0,-1
        wp(0,j) = vv1(0)*wpl(0,j) - sum(vv1(1:3)*wpl(1:3,j))
        if ( coefq .ne. 0 ) then
          do mu = 0,3
            riOut = incRI(mu,j)
            if ((riOut.gt.riMaxIn).and.firstRI(mu,j)) then
                   wp(0,riOut) = cc2 * wpl(mu,j)
            else;  wp(0,riOut) = wp(0,riOut) + cc2 * wpl(mu,j)
            endif
          enddo
        endif
      enddo
    end select

  case default ! -----------------------------------------------------

    if (warnings(301).le.warning_limit) then
      warnings(301) = warnings(301) + 1
      call openOutput
      write(nx,*)
      write(nx,*) "CODE ERROR 301 (loop_vertices_rcl.f90): wrong 3-leg interaction"
      write(nx,*)
      call toomanywarnings(301)
    endif
    call istop (ifail,2)

  end select  ! ------------------------------------------------------

  end subroutine loop3

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine loop4 ( riMaxOut,c,ty,wpl,wpt1,wpt2,wp )

  integer,     intent (in)    :: riMaxOut,ty
  complex(dp), intent (in)    :: c(3),wpl(0:3,0:riMaxOut), &
                                 wpt1(0:3),wpt2(0:3)
  complex(dp), intent (out)   :: wp(0:3,0:riMaxOut)

  integer     :: j
  complex(dp) :: c1,cc1,cc2,wpr,wpr1,wpr2,wpr3,wp2(0:3),wp3(0:3)

  c1 = c(1)

  ! select interaction type and contract currents

  select case (ty)

  case (-1) ! s + s + s -> s

    cc1 = - c1 * wpt1(0) * wpt2(0)
    wp(0,:) = cc1 * wpl(0,:)
    wp(1:3,:) = c0d0

  case (-2) ! s + v + v -> s

    wpr = wpt1(0)*wpt2(0) - sum(wpt1(1:3)*wpt2(1:3))
    cc1 = - c1 * wpr
    wp(0,  :) = cc1 * wpl(0,:)
    wp(1:3,:) = c0d0

  case (-3) ! v + s + v -> s

    cc1 = - c1 * wpt1(0)
    do j = 0,riMaxOut
      wpr = wpt2(0)*wpl(0,j) - sum(wpt2(1:3)*wpl(1:3,j))
      wp(0,j) = cc1 * wpr
    enddo
    wp(1:3,:) = c0d0

  case (-4) ! v + v + s -> s

    cc1 = - c1 * wpt2(0)
    do j = 0,riMaxOut
      wpr = wpt1(0)*wpl(0,j) - sum(wpt1(1:3)*wpl(1:3,j))
      wp(0,j) = cc1 * wpr
    enddo
    wp(1:3,:) = c0d0

  case (-11) ! v + v + v -> v

    wpr1 = wpt1(0)*wpt2(0) - sum(wpt1(1:3)*wpt2(1:3)) ! wp2*wp3
    cc1 = c(1) * wpr1
    wp2 = c(2) * wpt2
    wp3 = c(3) * wpt1
    do j = 0,riMaxOut
      wpr2 = wp2(0)*wpl(0,j) - sum(wp2(1:3)*wpl(1:3,j)) ! wp1*wp3*cc(2)
      wpr3 = wp3(0)*wpl(0,j) - sum(wp3(1:3)*wpl(1:3,j)) ! wp1*wp2*cc(3)
      wp(:,j) = cc1*wpl(:,j) + wpr2*wpt1 + wpr3*wpt2
    enddo

  case (-12) ! s + s + v -> v

    cc1 = c1 * wpt1(0)
    do j = 0,riMaxOut
      cc2 = cc1 * wpl(0,j)
      wp(:,j) = cc2 * wpt2(:)
    enddo

  case (-13) ! v + s + s -> v

    cc1 = c1 * wpt1(0) * wpt2(0)
    wp(:,:) = cc1 * wpl(:,:)

  case (-14) ! s + v + s -> v

    cc1 = c1 * wpt2(0)
    do j = 0,riMaxOut
      cc2 = cc1 * wpl(0,j)
      wp(:,j) = cc2 * wpt1(:)
    enddo

  case default

    if (warnings(302).le.warning_limit) then
      warnings(302) = warnings(302) + 1
      call openOutput
      write(nx,*)
      write(nx,*) "CODE ERROR 302 (loop_vertices_rcl.f90): wrong 4-leg interaction"
      write(nx,*)
      call toomanywarnings(302)
    endif
    call istop (ifail,2)

  end select

  end subroutine loop4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module loop_vertices_rcl

!#####################################################################

