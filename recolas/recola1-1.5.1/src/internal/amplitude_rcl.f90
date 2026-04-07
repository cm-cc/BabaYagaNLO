!#####################################################################
!!
!!  File  amplitude_rcl.f90
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

  module amplitude_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl
  use collier_interface_rcl
  use wave_functions_rcl
  use tables_rcl
  use model_vertices_rcl
  use tree_vertices_rcl
  use loop_vertices_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

  real(dp),    allocatable :: momenta(:,:)
  real(dp),    allocatable :: matrix2(:,:,:),matrix2h(:,:,:,:), &
                              matrix2int(:,:,:,:), &
                              matrix2cc(:,:,:,:), &
                              matrix2ccint(:,:,:,:), &
                              matrix2ccnlo(:,:,:,:), &
                              matrix2scc(:,:,:,:), &
                              matrix2sc(:,:), &
                              matrix2scm(:,:,:,:)
  complex(dp), allocatable :: matrix(:,:,:,:,:),matrixLO(:,:,:,:)
  logical                  :: momcheck

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_momenta (pr,pp)

  integer,  intent(in) :: pr
  real(dp), intent(in) :: pp(0:,:)

  integer  :: i,j,legs,lin,lout,i1,i2,i1min
  real(dp) :: p(0:size(pp,1)-1,size(pp,2)),                          &
              p0(0:size(pp,1)-1,size(pp,2)),                         &
              pIn(0:3),pOut(0:3),mOut,m,check,pv2,p0old,p0new,       &
              qmu(0:3),E,q(1:3),qs,p12(1:3),p12Xq,d3,delta,qtm,      &
              p1(1:3),p1Xq,p1q(1:3),p1t(1:3),p1ts,q1q(1:3),m1,d1,E1, &
              p2(1:3),p2Xq,p2q(1:3),p2t(1:3),p2ts,q2q(1:3),m2,d2,E2, &
              p1s,p2s

  p = pp

  lin  = legsIn(pr)
  lout = legsOut(pr)
  legs = lin + lout

  ! Check the number of external particles
  if (lin.lt.1) then
    momcheck = .false.
    if (warnings(341).le.warning_limit) then
      warnings(341) = warnings(341) + 1
      call openOutput
      write(nx,*)
      write(nx,'(1x,a,i2)') &
        'ERROR 341: Wrong number of incoming particles for process ',pr
      write(nx,*)
      call toomanywarnings(341)
    endif
    call istop (ifail,1)
  endif
  if (lout.lt.1) then
    momcheck = .false.
    if (warnings(342).le.warning_limit) then
      warnings(342) = warnings(342) + 1
      call openOutput
      write(nx,*)
      write(nx,'(1x,a,i2)') &
        'ERROR 342: Wrong number of outgoing particles for process ',pr
      write(nx,*)
      call toomanywarnings(342)
    endif
    call istop (ifail,1)
  endif
  if (legs.lt.3) then
    momcheck = .false.
    if (warnings(343).le.warning_limit) then
      warnings(343) = warnings(343) + 1
      call openOutput
      write(nx,*)
      write(nx,'(1x,a,i2)') &
        'ERROR 343: Wrong number of external particles for process ',pr
      write(nx,*)
      call toomanywarnings(343)
    endif
    call istop (ifail,1)
  endif

  ! Check that the dimension of the first index of "p" coincides
  ! with the number of external particles
  if (size(p,2).ne.legs) then
    momcheck = .false.
    if (warnings(344).le.warning_limit) then
      warnings(344) = warnings(344) + 1
      call openOutput
      write(nx,*)
      write(nx,'(1x,a,i2,a)') &
        'ERROR 344: Wrong phase-space point for process ',pr,';'
      write(nx,'(1x,a)') &
        '          the number of external momenta does not coincide'
      write(nx,'(1x,a)') &
        '          with the number of external particles'
      write(nx,*)
      call toomanywarnings(344)
    endif
    call istop (ifail,1)
  endif

  ! Check that all particles have positive energy larger than
  ! their mass
  do i = 1,legs
    m = mONS(newleg(i,pr),pr)
    if (p(0,i).le.0d0) then
      momcheck = .false.
      if (warnings(345).le.warning_limit) then
        warnings(345) = warnings(345) + 1
        call openOutput
        write(nx,*)
        write(nx,'(1x,a,i2,a)') &
          'ERROR 345: Wrong phase-space point for process ',pr,';'
        write(nx,'(1x,a,i2,a)') &
          '           particle ',i,' has not a positive energy'
        write(nx,*)
        call toomanywarnings(345)
      endif
      call istop (ifail,1)
    elseif (p(0,i).lt.m) then
      momcheck = .false.
      if (warnings(346).le.warning_limit) then
        warnings(346) = warnings(346) + 1
        call openOutput
        write(nx,*)
        write(nx,'(1x,a,i2,a)') &
          'ERROR 346: Wrong phase-space point for process ',pr,';'
        write(nx,'(1x,a,i2,a)') &
          '           energy of particle ',i,' is less than its mass'
        write(nx,*)
        write(nx,'(1x,a,i1,1x,a,1x,"(",f20.15,3(",",f21.15),")")') &
          'p',i,'=',p(0,i),p(1,i),p(2,i),p(3,i)
        write(nx,'(1x,a,i1,1x,a,2x,f20.15)') &
          'm',i,'=',m
        write(nx,*)
        write(nx,*) '           All amplitudes are set to 0'
        write(nx,*)
        call toomanywarnings(346)
      endif
      if (abs(p(0,i)-m)/(p(0,i)+m).gt.zerocheck) call istop (ifail,1)
      momenta = p
      return
    endif
  enddo

  do j = 0,3
    pIn(j) = sum(p(j,1:lin)); pOut(j) = sum(p(j,1+lin:legs))
  enddo

  mOut = 0d0
  do i = lin+1,legs
    mOut = mOut + mONS(newleg(i,pr),pr)
  enddo

  ! Check that the total incoming energy is sufficient to produce the
  ! on-shell outgoing particles
  if (pIn(0).lt.mOut) then
    momcheck = .false.
    if (warnings(347).le.warning_limit) then
      warnings(347) = warnings(347) + 1
      call openOutput
      write(nx,*)
      write(nx,'(1x,a,i2,a)') &
        'ERROR 347: Wrong phase-space point for process ',pr,';'
      write(nx,*) &
        '           not enough total incoming energy'
      write(nx,*)
      do i = 1, legs
        write(nx,'(1x,a,i1,1x,a,1x,"(",f20.15,3(",",f21.15),")")') &
          'p',i,'=',p(0,i),p(1,i),p(2,i),p(3,i)
      enddo
      write(nx,*)
      write(nx,*) '           All amplitudes are set to 0'
      write(nx,*)
      call toomanywarnings(347)
    endif
    if (abs(pIn(0)-mOut)/(pIn(0)+mOut).gt.zerocheck) &
      call istop (ifail,1)
    momenta = p
    return
  endif

  momcheck = .true.
  if (allocated(momenta)) deallocate(momenta)
  allocate (momenta(0:3,legs))
  momenta = p

  ! Check the mass shell condition of external momenta
  do i = 1, legs
    m = mONS(newleg(i,pr),pr)
    p0old = p(0,i)
    pv2 = p(1,i)*p(1,i) + p(2,i)*p(2,i) + p(3,i)*p(3,i)
    p0new = sqrt( pv2 + m*m )
    check = abs(p0old-p0new)/p0old
    if (check.gt.zerocut) then
      if (momenta_correction) then
        p(0,i) = p0new
        if (check.gt.zerocheck) then
          if (warnings(348).le.warning_limit) then
            warnings(348) = warnings(348) + 1
            call openOutput
            write(nx,*)
            write(nx,*) 'WARNING 348 !!!'
            write(nx,'(1x,a,i2,a,i2,2a)') &
              'Particle ',i,' of process ',pr,' is not on the mass shell. ', &
              'Its energy is rescaled:'
            write(nx,*) 'OLD momentum:'
            write(nx,'(1x,a,i1,1x,a,1x,"(",f20.15,3(",",f21.15),")")') &
              'p',i,'=',p0old,p(1,i),p(2,i),p(3,i)
            write(nx,*) 'NEW momentum:'
            write(nx,'(1x,a,i1,1x,a,1x,"(",f20.15,3(",",f21.15),")")') &
              'p',i,'=',p(0,i),p(1,i),p(2,i),p(3,i)
            write(nx,'(1x,a,i1,a,i1,a,i1,a,e10.3)') &
              '(p',i,'^2-m',i,'^2)/E',i,'^2 =', &
              abs(p(0,i)**2 - pv2 - m*m)/p(0,i)**2
            write(nx,*)
            call toomanywarnings(348)
          endif
        endif
      else
        if (warnings(349).le.warning_limit) then
          warnings(349) = warnings(349) + 1
          call openOutput
          write(nx,*)
          write(nx,*) 'ERROR 349: '
          write(nx,'(1x,a,i2,a,i2,2a)') &
            'Particle ',i,' of process ',pr,' is not on the mass shell. ', &
            'All amplitudes are set to 0'
          write(nx,*)
          call toomanywarnings(349)
        endif
        momcheck = .false.
        return
      endif
    endif
  enddo

  do j = 0,3
    pIn(j) = sum(p(j,1:lin)); pOut(j) = sum(p(j,1+lin:legs))
  enddo

  ! Check conservation of the four-momentum
  check = 0d0
  do j = 0,3
    check = max(check,abs(pIn(j)-pOut(j))/abs(pIn(0)))
  enddo
  if (check.gt.zerocut) then
    if (momenta_correction) then
      if (lout.eq.1) then; i1min = 1
      else;                i1min = lin+1
      endif
      i1loop: do i1 = i1min,legs
      i2loop: do i2 = i1+1,legs
        if (lout.eq.1) then
          p1 = p(1:3,1)
          p2 = p(1:3,2)
          qmu = p(:,1) + p(:,2) - pIn(:) + pOut(:)
        elseif (lout.ge.2) then
          p1 = p(1:3,i1)
          p2 = p(1:3,i2)
          qmu = p(:,i1) + p(:,i2) - pOut(:) + pIn(:)
        endif
        E = qmu(0); q = qmu(1:3); qs = dot_product(q,q)
        ! Each vector v is written as v = vq + vt, where vq=v.q/qs*q
        p1Xq = dot_product(p1,q)
        p2Xq = dot_product(p2,q)
        p12 = p1 + p2;  p12Xq = dot_product(p12,q)
        ! The new 4-vectors are
        ! q1mu = (E1,q1), q1 = q1q + q1t; q2mu = (E2,q2), q2 = q2q + q2t
        ! For momentum conservation q1t + q2t = 0 => |q1t| = |q2t| = qtm
        ! I choose q1q = sigma*p1q, q2q = sigma*p2q => sigma = qs/p12.q
        ! If statement fixed, AD 24.06.2024
!        if (abs(qs) .gt. zerocut) then; q1q = p1Xq/p12Xq * q; q2q = p2Xq/p12Xq * q
        if (abs(qs) .gt. zerocut * E**2) then
          q1q = p1Xq/p12Xq * q; q2q = p2Xq/p12Xq * q
        else;                q1q = 0d0;            q2q = 0d0
        endif
        ! d1 = m1^2 + |q1q|^2, d2 = m2^2 + |q2q|^2
        if (lout.eq.1) then
          m1 = mONS(newleg(1,pr),pr);  m2 = mONS(newleg(2,pr),pr)
        else
          m1 = mONS(newleg(i1,pr),pr); m2 = mONS(newleg(i2,pr),pr)
        endif
        d1 = m1*m1 + dot_product(q1q,q1q)
        d2 = m2*m2 + dot_product(q2q,q2q)
        d3 = E**2
        delta = d1**2 + d2**2 + d3**2 - 2*d1*d2 - 2*d1*d3 - 2*d2*d3
        if ( lout.eq.1 .or. &
             ( d3+d1-d2.ge.0d0.and. &
               d3+d2-d1.ge.0d0.and. &
               delta.ge.0d0         ) ) then
          exit i1loop
        endif
      enddo i2loop
      enddo i1loop
      if (d3+d1-d2.lt.0d0.or.d3+d2-d1.lt.0d0.or.delta.lt.0d0) then
        if (warnings(350).le.warning_limit) then
          warnings(350) = warnings(350) + 1
          call openOutput
          write(nx,*)
          write(nx,*) 'ERROR 350: '
          write(nx,'(1x,a,i2,a)') &
            'The four-momentum is not conserved in process ',pr,':'
          write(nx,'(1x,a,1x,"(",f20.15,3(",",f21.15),")")') &
            'pIn  =',pIn(0),pIn(1),pIn(2),pIn(3)
          write(nx,'(1x,a,1x,"(",f20.15,3(",",f21.15),")")') &
            'pOut =',pOut(0),pOut(1),pOut(2),pOut(3)
          write(nx,*) 'All amplitudes are set to 0'
          write(nx,*)
          call toomanywarnings(350)
          momcheck = .false.
          return
        endif
      endif
      ! q1A = E1 + qtm, q1B = E1 - qtm
      ! q2A = E2 + qtm, q2B = E2 - qtm
      ! Conditions are:
      ! q1A + q2B = E; q1B + q2A = E -> energy conservation
      ! q1A*q1B = d1;  q2A*q2B = d2  -> on-shell conditions
      ! I define x = q1A and get:
      ! E*x^2 + (d2-d1-E^2)*x + E*d1 = 0
      ! x+- = (d3 + d1 - d2 +- sq)/(2*E), sq = sqrt(delta)
      ! I choose:
      ! q1A = (d3 + d1 - d2 + sq)/(2*E) = x+
      ! q1B = d1/q1A = 2*d1*E/(d3 + d1 - d2 + sq) = x-
      ! I get
      ! q2A = (d3 + d2 - d1 + sq)/(2*E) = y+
      ! q2B = d2/q2A = 2*d2*E/(d3 + d2 - d1 + sq) = y-
      ! where y+- = (d3 + d2 - d1 +- sq)/(2*E)
      ! Then
      ! E1  = (q1A+q1B)/2 = (d3+d1-d2)/(2*E)
      ! E2  = (q2A+q2B)/2 = (d3+d2-d1)/(2*E)
      ! qtm = (q1A-q1B)/2 = (q2A-q2B)/2 = sq/(2*E)
      p0 = p
      E1 = (d3+d1-d2)/(2*E); E2 = (d3+d2-d1)/(2*E)
      qtm = sqrt(delta)/(2*E)
      p(0,i1) = E1;     p(0,i2) = E2
      ! If statement fixed, AD 24.06.2024
      ! should match if statement above!
!      if (qs.gt.0d0) then
      if (abs(qs) .gt. zerocut * E**2) then
        p1q = p1Xq/qs * q; p2q = p2Xq/qs * q
        p1t = p1 - p1q; p1ts = dot_product(p1t,p1t)
        p2t = p2 - p2q; p2ts = dot_product(p2t,p2t)
        p(1:3,i1) = q1q + qtm/sqrt(p1ts)*p1t
        p(1:3,i2) = q2q + qtm/sqrt(p2ts)*p2t
      else
        p1s = dot_product(p1,p1)
        p2s = dot_product(p2,p2)
        if (p1s.eq.0d0.and.p2s.eq.0d0) then
          p(1:3,i1) = 0d0
          p(1:3,i2) = 0d0
        elseif (p1s.ge.p2s) then
          p(1:3,i1) = + qtm/sqrt(p1s)*p1
          p(1:3,i2) = - qtm/sqrt(p1s)*p1
        else
          p(1:3,i1) = - qtm/sqrt(p2s)*p2
          p(1:3,i2) = + qtm/sqrt(p2s)*p2
        endif
      endif
      if (check.gt.zerocheck) then
        if (warnings(351).le.warning_limit) then
          warnings(351) = warnings(351) + 1
          call openOutput
          write(nx,*)
          write(nx,*) 'WARNING 351 !!!'
          write(nx,'(1x,a,i2,a)') &
            'The four-momentum is not conserved in process ',pr,':'
          write(nx,'(1x,a,1x,"(",f20.15,3(",",f21.15),")")') &
            'pIn  =',pIn(0),pIn(1),pIn(2),pIn(3)
          write(nx,'(1x,a,1x,"(",f20.15,3(",",f21.15),")")') &
            'pOut =',pOut(0),pOut(1),pOut(2),pOut(3)
          if (lout.eq.1) then
            write(nx,*) &
              'The four-momentum of the two incoming particles ', &
              'is rescaled:'
          elseif (lout.ge.2) then
            write(nx,*) &
              'The four-momentum of two outgoing particles is rescaled:'
          endif
          write(nx,*) 'OLD momenta:'
          do i = 1, legs
            write(nx,'(1x,a,i1,3x,a,1x,"(",f20.15,3(",",f21.15),")")') &
              'p',i,'=',p0(0,i),p0(1,i),p0(2,i),p0(3,i)
          enddo
          write(nx,*) 'NEW momenta:'
          do i = 1, legs
            write(nx,'(1x,a,i1,3x,a,1x,"(",f20.15,3(",",f21.15),")")') &
              'p',i,'=',p(0,i),p(1,i),p(2,i),p(3,i)
          enddo
          do i = 1, legs
            pv2 = p(1,i)*p(1,i) + p(2,i)*p(2,i) + p(3,i)*p(3,i)
            m = mONS(newleg(i,pr),pr)
            write(nx,'(1x,a,i1,a,i1,a,i1,a,e10.3)') &
              '(p',i,'^2-m',i,'^2)/E',i,'^2 =', &
              abs(p(0,i)**2 - pv2 - m*m)/p(0,i)**2
          enddo
          do j = 0,3
            pIn(j) = sum(p(j,1:lin)); pOut(j) = sum(p(j,1+lin:legs))
          enddo
          write(nx,'(1x,a,1x,"(",f20.15,3(",",f21.15),")")') &
            'pIn  =',pIn(0),pIn(1),pIn(2),pIn(3)
          write(nx,'(1x,a,1x,"(",f20.15,3(",",f21.15),")")') &
            'pOut =',pOut(0),pOut(1),pOut(2),pOut(3)
          write(nx,*)
          call toomanywarnings(351)
        endif
      endif
    else
      if (warnings(352).le.warning_limit) then
        warnings(352) = warnings(352) + 1
        call openOutput
        write(nx,*)
        write(nx,*) 'ERROR 352: '
        write(nx,'(1x,a,i2,a)') &
           'The four-momentum is not conserved in process ',pr,':'
        write(nx,'(1x,a,1x,"(",f20.15,3(",",f21.15),")")') &
           'pIn  =',pIn(0),pIn(1),pIn(2),pIn(3)
        write(nx,'(1x,a,1x,"(",f20.15,3(",",f21.15),")")') &
           'pOut =',pOut(0),pOut(1),pOut(2),pOut(3)
        write(nx,*) 'All amplitudes are set to 0'
        write(nx,*)
        call toomanywarnings(352)
        momcheck = .false.
        return
      endif
    endif
  endif

  momenta = p

  end subroutine set_momenta

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_amplitude (pr,order)

  integer,          intent(in) :: pr
  character(len=*), intent(in) :: order

  integer                  :: i,j,k,e,ii,jj,kk,jmax,kmax,s,w,legs,   &
                              j1,j2,jd,tl,tlm1,tlp1,fh,w0i,se,       &
                              semax(0:3),leg1,leg2,leg3,leg4,pa(4),  &
                              imin,imax,bimin,bimax,ih1,emin,c,      &
                              w0in(1:3),w0out,w1in,w1out,emax(0:3),  &
                              riMaxIn,riMaxOut,lp,t,lmu,b,locoef,ty, &
                              gsI,gs,cs,mo,da,ti,leg,n,ran,          &
                              cllaccuracy,cllerror,                  &
                              wid,widc,polsel,                       &
                              only_loop_prop = -1
  integer, parameter       :: nopol=-2
  integer, allocatable     :: w0l(:)
  logical                  :: computeNLO,x,last,ferlo,winit,         &
                              loop_prop(8) = .true.,                 &
                              loop_v(3:4) = .true.,                  &
                              ct_v  (2:4) = .true.,                  &
                              r2_v  (2:4) = .true.
  real(sp)                 :: timeTIin,timeTIout,timeTCin,timeTCout
  real(dp)                 :: fa,cllcritacc,check,                   &
                              mass_proj
  complex(dp)              :: iover32pi2,m2,m2p,p2p,den,co(4),fac,   &
                              wwTI,cc,spr,wp4(0:3),p4(0:3),          &
                              pol1(0:3),pol1_cv(0:3),pol1a(0:3),     &
                              pol1b(0:3),pol1a_cv(0:3),              &
                              pol1b_cv(0:3),contr_tmp(0:3),          &
                              contr_tmpa(0:3),contr_tmpb(0:3)
  complex(dp), allocatable :: p(:,:),pl(:,:),ps(:),pp(:),psp(:),     &
                              momInv(:),ms(:),momVec(:,:),           &
                              TIri(:,:),TIriUV(:),ww0(:,:),          &
                              ww1(:,:,:),ww0out(:,:),ww1out(:,:,:)


  ! Initialize matrix element
  if (.not.allocated(matrix)) then
    allocate (matrix(1:csMax,0:gsMax,1:cfMax,0:4,prTot))
    matrix = c0d0
  else
    matrix(1:pCsTot(pr),0:gsTot(0,pr),1:cfTot(pr),0,pr) = c0d0
    if (order.eq.'NLO') &
      matrix(1:pCsTot(pr),0:gsTot(1,pr),1:cfTot(pr),1:4,pr) = c0d0
  endif

  if (.not.prexists(pr)) return

  if (.not.momcheck) return

!  loop_v(4) = .false.

!  only_loop_prop = 6
!  ct_v = .false.
!  r2_v = .false.

!  only_loop_prop = 5
!  ct_v = .false.
!  r2_v = .false.

!  only_loop_prop = 4
!  ct_v(2) = .false.; ct_v(3) = .false.
!  r2_v(2) = .false.; r2_v(3) = .false.

!  only_loop_prop = 3
!  ct_v(2) = .false.; ct_v(4) = .false.
!  r2_v(2) = .false.; r2_v(4) = .false.

!  only_loop_prop = 2
!  ct_v(3) = .false.; ct_v(4) = .false.
!  r2_v(3) = .false.; r2_v(4) = .false.

!  only_loop_prop = 1
!  ct_v = .false.
!  r2_v = .false.

  do j = 1,8
    if (only_loop_prop.eq.j) then
      loop_prop    = .false.
      loop_prop(j) = .true.
    endif
  enddo

  iover32pi2 = cId0/(32d0*pi*pi)

  computeNLO = (order.eq.'NLO').and.(lpmax(pr).gt.0)

  legs = legsIn(pr) + legsOut(pr)

  tlm1 = 2**(legs-1)
  tl   = 2**legs
  tlp1 = 2**(legs+1)

  if (computeNLO) then
    allocate (p(0:3,tlp1-1))
    allocate (pl(4,tlp1-1))
    allocate (psp(tl-1))
  else
    allocate (p(0:3,tl-1))
    allocate (pl(4,tl-1))
  endif
  allocate (ps(tl-1))
  allocate (pp(tl-1))

  ! build external momenta (all "p" are incoming) and introduce:
  !
  !   pl(1) = p(0)+p(3)
  !   pl(2) = p(0)-p(3)
  !   pl(3) = p(1)+i*p(2)
  !   pl(4) = p(1)-i*p(2)
  !
  !   ps = p(0)*p(0)-p(1)*p(1)-p(2)*p(2)-p(3)*p(3)
  !      = pl(1)*pl(2)-pl(3)*pl(4)

  p    = c0d0
  pl   = c0d0
  ps   = c0d0

  do i = 1,legsIn(pr); j = 2**(newleg(i,pr)-1)
    p(:,j) = cmplx(momenta(:,i),kind=dp)
  enddo
  do i = legsIn(pr)+1,legs; j = 2**(newleg(i,pr)-1)
    p(:,j) = - cmplx(momenta(:,i),kind=dp)
  enddo

  do i= 1, legs; j = 2**(i-1)
    pl(1,j) = p(0,j) + p(3,j)
    pl(2,j) = p(0,j) - p(3,j)
    pl(3,j) = p(1,j) + cId0*p(2,j)
    pl(4,j) = conjg(pl(3,j))
    ! cmONS2 is the squared of the on-shell mass of the external
    ! particle:
    ! - it is 0 if the particle is marked as light
    ! - it is the squared of the input mass if the particle is not
    !   marked as light
    ps(j) = cmONS2(i,pr)
    jd = tl - 1 - j
    p (:,jd) = - p (:,j)
    pl(:,jd) = - pl(:,j)
    ps(  jd) =   ps(  j)
  enddo

  bin: do j = 3,tl-2
    i = levelLeg(j)
    if (i.eq.1) then
      cycle bin
    elseif (i.le.legs/2) then
      j1 = firstNumbers(j,1)
      j2 = sum(firstNumbers(j,2:i))
      p (:,j) = p(:,j1) + p(:,j2)
      pl(1,j) = p(0,j) + p(3,j)
      pl(2,j) = p(0,j) - p(3,j)
      pl(3,j) = p(1,j) + cId0*p(2,j)
      pl(4,j) = conjg(pl(3,j))
      pp(j) =           &
      + p(0,j1)*p(0,j2) &
      - p(1,j1)*p(1,j2) &
      - p(2,j1)*p(2,j2) &
      - p(3,j1)*p(3,j2)
      ps(j) = ps(j1) + ps(j2) + 2d0*pp(j)
      jd = tl - 1 - j
      p (:,jd) = - p (:,j)
      pl(:,jd) = - pl(:,j)
      ps(  jd) =   ps(  j)
    endif
  enddo bin

  ! objects with argument 2**legs-1 vanish (energy-momentum
  ! conservation

  if (computeNLO) then

    ! psp: squared momenta for tensor integrals
    psp = ps

    do j = 3,tlm1-1
      if (defresbin(j,pr)) then
        check = abs(psp(j)-pspbin(j,pr))/pspbin(j,pr)
        if (check.gt.zerocheck) then
          if (warnings(353).le.warning_limit) then
            warnings(353) = warnings(353) + 1
            call openOutput
            write(nx,*)
            write(nx,*) 'ERROR 353: set_resonant_particle_rcl called '
            write(nx,*) '           for a resonance not on the mass-shell'
            write(nx,*)
            call toomanywarnings(353)
          endif
          call istop (ifail,1)
        endif
        psp(j) = pspbin(j,pr)
        psp(tl-1-j) = pspbin(j,pr)
      endif
    enddo
    do i = 1, legs;  j = 2**(i-1)
      ! cmREG2 is the squared of the input mass of the external
      ! particle (if the particle is marked as light, mREG2 is used
      ! as IR-regulator in loop propagators only)
      psp(j) = cmREG2(i,pr)
      psp(tl-1-j) = psp(j)
    enddo

    ! additional definitions for loop legs (offsets)
    do i = 1, tl-1
      p (:,tl+i) = p (:,i)
      pl(:,tl+i) = pl(:,i)
    enddo

    call cpu_time (timeTIin)

    if (dynamic_settings.ge.1) then
      call SetDeltaUV_cll (deltaUV)
      call SetMuUV2_cll(muUV**2)
      call SetDeltaIR_cll (deltaIR,deltaIR2)
      if (reg_soft.eq.1) then; call SetMuIR2_cll(muIR**2)
      else;                    call SetMuIR2_cll(lambda**2)
      endif
    endif

    allocate (momInv(0:legs*(legs-1)/2))
    allocate (TIri(  0:ritiMax(pr),tiTot(pr)))

    do i = 1,tiTot(pr)

      ! cache
      if (mod(i-1,tiCache(pr)).eq.0) then
        n = nCacheTot(pr-1) + (i-1)/tiCache(pr) + 1
        call InitEvent_cll(n)
      endif

      leg = legsti(i,pr)
      ran = rankti(i,pr)

      allocate (ms(0:leg-1)); ms= c0d0
      do j = 1, leg
        ms(j-1) = cm2n(vmti(j,i,pr)) ! internal
      enddo

      allocate (TIriUV(0:riMax(ran)))

      select case (leg)

      case (1) ! tadpoles do not require momenta nor invariants

        call TNten_cll ( TIri(  0:riMax(ran),i), & ! output
                         TIriUV(0:riMax(ran)),   & ! output
                         ms(0:0),leg,ran         ) ! input

        if (.not.loop_prop(1)) TIri(0:riMax(ran),i) = c0d0

      case default ! arguments for more than one internal leg

        allocate (momVec(0:3,0:legs-1))

        ! momenta
        imin             = 1
        imax             = leg-1
        ii = 0
        do k = 1,legs
          if (momsti(k,i,pr).eq.imin) ii = ii + 2**(k-1)
        enddo
        momVec(0,  imin) = + p(0,ii)
        momVec(1:3,imin) = - p(1:3,ii)
        do j = imin+1, imax
          ii = 0
          do k = 1,legs
            if (momsti(k,i,pr).eq.j) ii = ii + 2**(k-1)
          enddo
          momVec(0,  j) = momVec(0,  j-1) + p(0,  ii)
          momVec(1:3,j) = momVec(1:3,j-1) - p(1:3,ii)
        enddo

        ! invariants
        bimin = 1
        bimax = leg*(leg-1)/2
        jmax = leg/2   ! leg=2*n -> jmax=n;     leg=2*n+1 -> jmax=n
        kmax = leg - 1 ! leg=2*n -> kmax=2*n-1; leg=2*n+1 -> kmax=2*n
        ii = 0
        do j = 1,jmax
          if (2*j.eq.leg) kmax = leg/2 - 1 ! leg=2*n & j=n -> kmax=n
          do k = 0,kmax
            jj = 0
            kk = 0
            do n = 1,legs
              if (momsti(n,i,pr).gt.0) then
                if (momsti(n,i,pr).le.mod(k+j,leg)) jj = jj + 2**(n-1)
                if (momsti(n,i,pr).le.k) kk = kk + 2**(n-1)
              endif
            enddo
            ii = ii + 1
            momInv(ii) = psp(abs(jj-kk))
          enddo
        enddo

        call TNten_cll ( TIri(  0:riMax(ran),i), & ! output
                         TIriUV(0:riMax(ran)),   & ! output
                         momVec(0:3,imin:imax),  & ! input
                         momInv(bimin:bimax),    & ! input
                         ms(0:imax),leg,ran      ) ! input

        do j = 1,8
          if (leg.eq.j.and.(.not.loop_prop(j))) TIri(0:riMax(ran),i) = c0d0
        enddo

        deallocate (momVec)

      end select

      deallocate (TIriUV)
      deallocate (ms)

    enddo ! tensor integrals evaluated

    deallocate (momInv)

    call GetAccFlag_cll(cllaccuracy)
    select case (cllaccuracy)
    case (-2)
      if (warnings(354).le.warning_limit) then
        warnings(354) = warnings(354) + 1
        call GetCritAcc_cll(cllcritacc)
        call openOutput
        write(nx,*)
        write(nx,*) 'WARNING 354 !!!'
        write(nx,'(1x,2a,g7.1)') &
          'Bad phase-space point: ', &
          'Accuracy of tensor integrals less than ',cllcritacc
        write(nx,*)
        if (writeMat+writeMat2.eq.0) then
          do i = 1,legs
            write(nx,'(2x,a,i1,1x,a,1x,"(",f15.10,3(",",f16.10),")")') &
              'p',i,'=',momenta(0:3,i)
          enddo
          write(nx,*)
        endif
        call toomanywarnings(354)
      endif
    end select

    call GetErrFlag_cll(cllerror)
    if (cllerror.lt.-8) then
      if (warnings(355).le.warning_limit) then
        warnings(355) = warnings(355) + 1
        call openOutput
        write(nx,*)
        write(nx,*) 'CODE ERROR 355 (amplitude_rcl):', &
                    'Collier called in a wrong way'
        write(nx,*)
        call toomanywarnings(355)
      endif
      call istop (ifail,2)

    endif

    call cpu_time (timeTIout)
    timeTI(pr) = timeTI(pr) + timeTIout - timeTIin

  endif

  call cpu_time (timeTCin)

  allocate (ww0(0:3,w0Tot(pr)))
  allocate (ww0out(0:3,0:modaTot(pr)))
  if (computeNLO) then
    allocate (ww1out(0:3,0:ritiMax(pr),0:modaTot(pr)))
  endif

  if ((dynamic_settings.ge.1).and.computeNLO) &
    call counterterms_tables

  do w = 1,w0eTot(pr)
    call definewp ( parw0e(w,pr),                                 &
                    p(:,binw0e(w,pr)),pl(:,binw0e(w,pr)),         &
                    mONS(legw0e(w,pr),pr),helw0e(w,pr),ww0(:,w) )
    if (longitudinal.eq.0) cycle
    select case (parw0e(w,pr)); case (15,16)
      if ( longitudinal.eq.111 .or.               &
           longitudinal.eq.oldleg(legw0e(w,pr),pr) ) then
        ww0(:,w) = p(:,binw0e(w,pr))/p(0,binw0e(w,pr))
      endif
    end select
  enddo

  allocate (w0l(cfTot(pr)))
  w0l(1:cfTot(pr)) = w0last(heli(legs,1:cfTot(pr),pr),pr)

  ! compute emax
  emax(0) = tlm1-1
  if (computeNLO) emax(1:3) = tl-1

  ! initialize matrix element
  if (.not.allocated(matrix)) then
    allocate (matrix(1:csMax,0:gsMax,1:cfMax,0:4,prTot))
    matrix = c0d0
  else
    matrix(1:pCsTot(pr),0:gsTot(0,pr),1:cfTot(pr),0,pr) = c0d0
    if (order.eq.'NLO') &
      matrix(1:pCsTot(pr),0:gsTot(1,pr),1:cfTot(pr),1:4,pr) = c0d0
  endif

  ! tree branches

  semax(0:1) = 0
  semax(2:3) = 1

  ! Index "c" combines the indices lp, t and fh
  cloop0: do c = 1,c0EffMax(pr)

    lp = c0TOlp(c,pr)
    if ((lp.gt.0).and.(.not.computeNLO)) cycle cloop0

    configloop0: do i = 1, cfTot(pr)

      w0i = w0l(i)

      eloop0: do e = 3,emax(lp)

        do se = 0,semax(lp)

          if (bm0min(e,i,c,pr).ne.0) then

            bm0loop: do b = bm0min(e,i,c,pr),bm0max(e,i,c,pr)

              s = sbm0(b,pr)

              mo    = mosm0(s,pr)
              leg1  = binsm0(1,s,pr)
              leg2  = binsm0(2,s,pr)
              leg3  = binsm0(3,s,pr)
              leg4  = leg1 + leg2 + leg3
              pa    = parsm0(1:4,s,pr)
              x     = xsm0(s,pr)
              gsI   = gsIncsm0(s,pr)
              cs    = cssm0(s,pr); last = (cs.ne.0)
              if (last) gs = gssm0(s,pr)

              if ((se.eq.0).and.(leg1.eq.leg4)) cycle
              if ((se.eq.1).and.(leg1.ne.leg4)) cycle

              w0in  = w0inbm0(1:3,b,pr)
              w0out = w0outbm0(b,pr)
              winit = winitbm0(b,pr)
              ty    = typebm0(b,pr)
              m2    = cm2f(pa(4))

              m2p = m2
              if (defresbin(leg4,pr)) m2p = cm2pf(pa(4))
              p2p = ps(leg4)
              if (defp2bin(leg4,pr)) p2p = p2bin(leg4,pr)
              den = p2p - m2p

              co = cosm0(:,s,pr)
              if ((dynamic_settings.eq.1).and.x.and.(lp.eq.2)) &
                co = co*cou(lp,gsI,pa)

              if (leg2.eq.0) then

                if (x.and.(lp.eq.2).and.(.not.ct_v(2))) co = c0d0
                if (x.and.(lp.eq.3).and.(.not.r2_v(2))) co = c0d0
                call tree2 (ps(leg4),p(:,leg1),pl(:,leg1),m2,den,  &
                            co(1:4),ty,ww0(:,w0in(1)),ww0out(:,mo) )
                wid  = pa(4)
                if ((polprojin(leg4,pr) .ge. -1) .and. &
                    (polprojin(leg4,pr) .le. 1) .and. &
                     parKind(wid) .eq. 3) then

                  wid = 17
                  mass_proj = sqrt(real(m2,kind=dp))
                  call definewp (wid,      &
                       p(:,leg4), &
                       pl(:,leg4),&
                       mass_proj, &
                       polprojin(leg4,pr), &
                       pol1)
                  pol1_cv(0) = pol1(0)
                  pol1_cv(1:3) = -pol1(1:3)
                  contr_tmp = sum(ww0out(:,mo)*pol1_cv(:))*conjg(pol1(:))
                  ww0out(:,mo) = -contr_tmp
                end if

                if (polprojin(leg4,pr) .eq. 3 .and. parKind(wid) .eq. 3) then
                   wid = 17
                   mass_proj = sqrt(real(m2,kind=dp))
                   call definewp (wid,       &
                        p(:,leg4), &
                        pl(:,leg4),&
                        mass_proj,-1,&
                        pol1a) !  left pol vec
                   pol1a_cv(0) = pol1a(0)
                   pol1a_cv(1:3) = -pol1a(1:3)
                   contr_tmpa = sum(ww0out(:,mo)*pol1a_cv(:))*conjg(pol1a(:))
                   call definewp (wid,       &
                        p(:,leg4), &
                        pl(:,leg4),&
                        mass_proj,+1,&
                        pol1b) ! right pol vec
                   pol1b_cv(0) = pol1b(0)
                   pol1b_cv(1:3) = -pol1b(1:3)
                   contr_tmpb = sum(ww0out(:,mo)*pol1b_cv(:))*conjg(pol1b(:))
                   ww0out(:,mo) = -contr_tmpa-contr_tmpb
                end if

                if ( (abs(polprojin(leg4,pr)).eq.1 .or. &
                     polprojin(leg4,pr) .eq. 3) .and. &
                     (parKind(wid) .eq. 4 .or. parKind(wid) .eq. 5)) then
                     write(*,*) "This part of the code needs to be generalized" // &
                                "See tree3 part where propagator is truncated."
                     stop 9
                end if

              elseif (leg3.eq.0) then

                if (x.and.(lp.eq.2).and.(.not.ct_v(3))) co = c0d0
                if (x.and.(lp.eq.3).and.(.not.r2_v(3))) co = c0d0
                call tree3 ( last,p(:,leg1),p(:,leg2),                  &
                             pl(:,leg1),pl(:,leg2),m2,den,co(1:2),ty,   &
                             ww0(:,w0in(1)),ww0(:,w0in(2)),ww0out(:,mo) )
                wid  = pa(4)

                ! Transverse polarization: modified by Giovanni to enable
                ! coherent sum of left and right pol
                if ((polprojin(leg4,pr) .ge. -1) .and. &
                    (polprojin(leg4,pr) .le. 1) .and. &
                     parKind(wid) .eq. 3) then
                  ! compute polarization for a massive vector boson (17 just a
                  ! dummy W-boson with modified mass)
                  wid = 17
                  mass_proj = sqrt(real(m2,kind=dp))
                  call definewp (wid,      &
                      p(:,leg4), &
                      pl(:,leg4),&
                      mass_proj, &
                      polprojin(leg4,pr), &
                      pol1)
                  pol1_cv(0) = pol1(0)
                  pol1_cv(1:3) = -pol1(1:3)
                  contr_tmp = sum(ww0out(:,mo)*pol1_cv(:))*conjg(pol1(:))
                  ww0out(:,mo) = -contr_tmp
                end if

                ! coherent sum of left- and right-handed polarization vector.
                if (polprojin(leg4,pr) .eq. 3 .and. parKind(wid) .eq. 3) then
                   wid = 17
                   mass_proj = sqrt(real(m2,kind=dp))
                   call definewp (wid,       &
                        p(:,leg4), &
                        pl(:,leg4),&
                        mass_proj,-1,&
                        pol1a) !  left pol vec
                   pol1a_cv(0) = pol1a(0)
                   pol1a_cv(1:3) = -pol1a(1:3)
                   contr_tmpa = sum(ww0out(:,mo)*pol1a_cv(:))*conjg(pol1a(:))
                   call definewp (wid,       &
                        p(:,leg4), &
                        pl(:,leg4),&
                        mass_proj,+1,&
                        pol1b) ! right pol vec
                   pol1b_cv(0) = pol1b(0)
                   pol1b_cv(1:3) = -pol1b(1:3)
                   contr_tmpb = sum(ww0out(:,mo)*pol1b_cv(:))*conjg(pol1b(:))
                   ww0out(:,mo) = -contr_tmpa-contr_tmpb
                end if

                ! Ghosts & Goldstone boson - transverse
                if (drop_goldstone_transverse .and. &
                    ((polprojin(leg4,pr) .eq. -1) .or. &
                     (polprojin(leg4,pr) .eq. 1)) .and. &
                     parKind(wid) .eq. 0) then
                  ww0out(:,mo) = 0
                end if
                ! Ghosts & Goldstone boson - coherent sum
                if (drop_goldstone_transverse .and. &
                    (polprojin(leg4,pr) .eq. 3) .and. &
                    parKind(wid) .eq. 0) then
                  ww0out(:,mo) = 0
                end if
                ! Ghosts & Goldstone boson - longitudinal
                if (drop_goldstone_longitudinal .and. &
                    ((polprojin(leg4,pr) .eq. 0)) .and. &
                    parKind(wid) .eq. 0) then
                  ww0out(:,mo) = 0
                end if

                ! Fermion polarization projection
                if (((polprojin(leg4,pr) .eq. -1) .or. &
                     (polprojin(leg4,pr) .eq. +1)) .and. &
                     (parKind(wid) .eq. 4 .or. parKind(wid) .eq. 5)) then
                  widc = anti(wid)
                  if (polprojin(leg4,pr) .eq. +1) then
                    if (parKind(wid) .eq. 4) then
                      polsel = -1
                    else
                      polsel = 1
                    end if
                  else
                    if (parKind(wid) .eq. 4) then
                      polsel =  1
                    else
                      polsel =  -1
                    end if
                  end if

                  mass_proj = sqrt(real(m2,kind=dp))
                  call definewp (  &
                        widc,      &
                        p(:,leg4), &
                        pl(:,leg4),&
                        mass_proj, &
                        -polsel,   &
                        pol1a)

                  call definewp (  &
                       wid,        &
                       -p(:,leg4), &
                       -pl(:,leg4),&
                       mass_proj,  &
                       polsel,     &
                       pol1b)

                  ! check
                  !call definewp ( &
                  !      widc,       &
                  !      p(:,leg4), &
                  !      pl(:,leg4),&
                  !      mass_proj, &
                  !      polsel,    &
                  !      pol1c)
                  ! call definewp ( &
                  !      wid,      &
                  !      -p(:,leg4), &
                  !      -pl(:,leg4),&
                  !      mass_proj, &
                  !      -polsel,    &
                  !      pol1d)

                  ! checks
                  !! sum(u ubar) = pslash + m
                  !do k = 1, 4
                  !  do j = 1, 4
                  !    pslm_2(k,j) = pol1a(k-1)*pol1b(j-1) + pol1c(k-1)*pol1d(j-1)
                  !  end do
                  !end do
                  !do k = 0, 3
                  !  write(*,*) "ww0out(k,mo):", ww0out(k,mo)
                  !end do

                  ! compute current and truncate full propagator
                  call tree3 ( .true.,p(:,leg1),p(:,leg2),                &
                               pl(:,leg1),pl(:,leg2),m2,den,co(1:2),ty,   &
                               ww0(:,w0in(1)),ww0(:,w0in(2)),ww0out(:,mo) )

                  ! add back propagator denominator
                  ww0out(:,mo) = cId0*ww0out(:,mo)/den

                  ! write(*,*)
                  ! do k = 1, 4
                  !   write(*,*) "sum(ww0out(0:,mo)*pslm_2(1:,k)):", &
                  !   sum(ww0out(0:,mo)*pslm_2(k, 1:))
                  ! end do

                  if (parKind(wid) .eq. 4) then
                    ww0out(:,mo) = sum(ww0out(:,mo)*pol1b(:))*pol1a(:)
                  else
                    ! this sign has been fixed by the requirement that the
                    ! full polarization sum should reproduce the original
                    ! current, i.e.:
                    ! ww0out(:,mo) = -sum(ww0out(:,mo)*pol1b(:))*pol1a(:)
                    !                -sum(ww0out(:,mo)*pol1d(:))*pol1c(:)
                    !              == ww0out_orig(:,mo)
                    ww0out(:,mo) = -sum(ww0out(:,mo)*pol1b(:))*pol1a(:)
                  end if
                end if

                if ((polprojin(leg4,pr) .eq. 3) .and. &
                    (parKind(wid) .eq. 4 .or. parKind(wid) .eq. 5)) then
                     write(*,*) "Coherent polarization sum not implemented" // &
                                " for massive fermions."
                     stop 9
                end if

                ! just for checks
                !if (polprojin(leg4,pr) .eq. 4) then
                !   wid = 17
                !   mass_proj = sqrt(real(m2,kind=dp))
                !   call definewp (wid,       &
                !        p(:,leg4), &
                !        pl(:,leg4),&
                !        mass_proj,-1,&
                !        pol1a) !  left pol vec
                !   pol1a_cv(0) = pol1a(0)
                !   pol1a_cv(1:3) = -pol1a(1:3)
                !   contr_tmpa = sum(ww0out(:,mo)*pol1a_cv(:))*conjg(pol1a(:))
                !   call definewp (wid,       &
                !        p(:,leg4), &
                !        pl(:,leg4),&
                !        mass_proj,+1,&
                !        pol1b) ! right pol vec
                !   pol1b_cv(0) = pol1b(0)
                !   pol1b_cv(1:3) = -pol1b(1:3)
                !   contr_tmpb = sum(ww0out(:,mo)*pol1b_cv(:))*conjg(pol1b(:))
                !
                !   call definewp (wid,       &
                !        p(:,leg4), &
                !        pl(:,leg4),&
                !        mass_proj,0,&
                !        pol1c) ! longitudinal pol vec
                !   pol1c_cv(0) = pol1c(0)
                !   pol1c_cv(1:3) = -pol1c(1:3)
                !   contr_tmpc = sum(ww0out(:,mo)*pol1c_cv(:))*conjg(pol1c(:))
                !
                !   call definewp (wid,       &
                !        p(:,leg4), &
                !        pl(:,leg4),&
                !        mass_proj,2,&
                !        pol1d) ! auxiliary pol vec
                !   pol1d_cv(0) = pol1d(0)
                !   pol1d_cv(1:3) = -pol1d(1:3)
                !   contr_tmpd = sum(ww0out(:,mo)*pol1d_cv(:))*conjg(pol1d(:))
                !
                !   ww0out(:,mo) = -contr_tmpa-contr_tmpb-contr_tmpc-contr_tmpd
                !end if

              else

                if (x.and.(lp.eq.2).and.(.not.ct_v(4))) co = c0d0
                if (x.and.(lp.eq.3).and.(.not.r2_v(4))) co = c0d0
                call tree4 (last,den,co(1:3),ty,           &
                            ww0(:,w0in(1)),ww0(:,w0in(2)), &
                            ww0(:,w0in(3)),ww0out(:,mo)    )
                wid  = pa(4)
                ! Transverse polarization: modified by Giovanni to enable
                ! coherent sum of left and right pol
                if ((polprojin(leg4,pr) .ne. nopol) .and. &
                    (polprojin(leg4,pr) .ne. 3) .and. &
                    (parKind(wid) .eq. 3) ) then
                   wid = 17
                   mass_proj = sqrt(real(m2,kind=dp))
                   call definewp (wid, &
                        p(:,leg4),     &
                        pl(:,leg4),    &
                        mass_proj,     &
                        polprojin(leg4,pr), &
                        pol1)
                   pol1_cv(0) = pol1(0)
                   pol1_cv(1:3) = -pol1(1:3)
                   contr_tmp = sum(ww0out(:,mo)*pol1_cv(:))*conjg(pol1(:))
                   ww0out(:,mo) = -contr_tmp
                end if
                ! coherent sum of left- and right-handed polarization vector.
                if ((polprojin(leg4,pr) .eq. 3) .and. (parKind(wid) .eq. 3)) then
                  wid = 17
                  mass_proj = sqrt(real(m2,kind=dp))
                  call definewp (wid,       &
                       p(:,leg4), &
                       pl(:,leg4),&
                       mass_proj,-1,&
                       pol1a) !  left pol vec
                  pol1a_cv(0) = pol1a(0)
                  pol1a_cv(1:3) = -pol1a(1:3)
                  contr_tmpa = sum(ww0out(:,mo)*pol1a_cv(:))*conjg(pol1a(:))

                  call definewp (wid,       &
                       p(:,leg4), &
                       pl(:,leg4),&
                       mass_proj,+1,&
                       pol1b) ! right pol vec
                  pol1b_cv(0) = pol1b(0)
                  pol1b_cv(1:3) = -pol1b(1:3)
                  contr_tmpb = sum(ww0out(:,mo)*pol1b_cv(:))*conjg(pol1b(:))
                  ww0out(:,mo) = -contr_tmpa-contr_tmpb
                end if

                ! Ghosts & Goldstone boson - transverse
                if ( drop_goldstone_transverse .and. &
                     ( (polprojin(leg4,pr) .eq. -1) .or. &
                       (polprojin(leg4,pr) .eq. 1)       ) .and. &
                     parKind(wid) .eq. 0                         ) then
                  ww0out(:,mo) = 0
                end if
                ! Ghosts & Goldstone boson - coherent sum
                if (drop_goldstone_transverse .and. &
                    (polprojin(leg4,pr) .eq. 3) .and. &
                    parKind(wid) .eq. 0) then
                  ww0out(:,mo) = 0
                end if
                ! Ghosts & Goldstone boson - longitudinal
                if (drop_goldstone_longitudinal .and. &
                    ((polprojin(leg4,pr) .eq. 0)) .and. &
                    parKind(wid) .eq. 0) then
                  ww0out(:,mo) = 0
                end if

                ! just for checks
                !if (polprojin(leg4,pr) .eq. 4) then
                !   wid = 17
                !   mass_proj = sqrt(real(m2,kind=dp))
                !   call definewp (wid,       &
                !        p(:,leg4), &
                !        pl(:,leg4),&
                !        mass_proj,-1,&
                !        pol1a) !  left pol vec
                !   pol1a_cv(0) = pol1a(0)
                !   pol1a_cv(1:3) = -pol1a(1:3)
                !   contr_tmpa = sum(ww0out(:,mo)*pol1a_cv(:))*conjg(pol1a(:))
                !
                !   call definewp (wid,       &
                !        p(:,leg4), &
                !        pl(:,leg4),&
                !        mass_proj,+1,&
                !        pol1b) ! right pol vec
                !   pol1b_cv(0) = pol1b(0)
                !   pol1b_cv(1:3) = -pol1b(1:3)
                !   contr_tmpb = sum(ww0out(:,mo)*pol1b_cv(:))*conjg(pol1b(:))
                !
                !   call definewp (wid,       &
                !        p(:,leg4), &
                !        pl(:,leg4),&
                !        mass_proj,0,&
                !        pol1c) ! longitudinal pol vec
                !   pol1c_cv(0) = pol1c(0)
                !   pol1c_cv(1:3) = -pol1c(1:3)
                !   contr_tmpc = sum(ww0out(:,mo)*pol1c_cv(:))*conjg(pol1c(:))
                !
                !   call definewp (wid,       &
                !        p(:,leg4), &
                !        pl(:,leg4),&
                !        mass_proj,2,&
                !        pol1d) ! auxiliary pol vec
                !   pol1d_cv(0) = pol1d(0)
                !   pol1d_cv(1:3) = -pol1d(1:3)
                !   contr_tmpd = sum(ww0out(:,mo)*pol1d_cv(:))*conjg(pol1d(:))
                !
                !   ww0out(:,mo) = -contr_tmpa-contr_tmpb-contr_tmpc-contr_tmpd
                !end if

              endif

              if (transverse_resonance .and. defresbin(leg4,pr)) then
                select case (pa(4))
                case (12,13,14)
                  ww0out(:,mo) = c0d0
                case (17,18,19)
                  wp4 = ww0out(:,mo)
                  p4 = p(:,leg4)
                  spr = p4(0)*wp4(0) - sum(p4(1:3)*wp4(1:3))
                  cc = spr/ps(leg4)
                  ww0out(:,mo) = wp4(:) - cc * p4(:)
                end select
              endif

              if (last) then
                wwTI = sum(ww0out(:,mo)*ww0(:,w0i))
                matrix(cs,gs,i,lp,pr) = matrix(cs,gs,i,lp,pr) + wwTI
              else
                if (winit) then
                  ww0(:,w0out) = ww0out(:,mo)
                else
                  ww0(:,w0out) = ww0(:,w0out) + ww0out(:,mo)
                endif
              endif

            enddo bm0loop

          endif

          if (bd0min(e,i,c,pr).ne.0) then

            bd0loop: do b = bd0min(e,i,c,pr),bd0max(e,i,c,pr)

              s = sbd0(b,pr)

              if ((se.eq.0).and.(sesd0(s,pr))) cycle
              if ((se.eq.1).and.(.not.sesd0(s,pr))) cycle

              da    = dasd0(s,pr)
              fa    = facsd0(s,pr)
              cs    = cssd0(s,pr); last = (cs.ne.0)
              if (last) gs = gssd0(s,pr)

              w0out = w0outbd0(b,pr)
              winit = winitbd0(b,pr)

              if (last) then
                wwTI = sum(ww0out(:,da)*ww0(:,w0i))
                matrix(cs,gs,i,lp,pr) = matrix(cs,gs,i,lp,pr) + fa * wwTI
              else
                if (winit) then
                  ww0(:,w0out) = fa * ww0out(:,da)
                else
                  ww0(:,w0out) = ww0(:,w0out) + fa * ww0out(:,da)
                endif
              endif

            enddo bd0loop

          endif

        enddo

      enddo eloop0

    enddo configloop0

  enddo cloop0

  deallocate (w0l)

  if (order.eq.'NLO') then
    if (.not.allocated(matrixLO)) &
      allocate (matrixLO(1:csMax,0:gsMax,1:cfMax,prTot))
    matrixLO(:,:,:,pr) = matrix(:,:,:,0,pr)
    do gs = 0,gsTot(0,pr)
      if (powgs(gs,0,pr).eq.0) matrix(:,gs,:,0,pr) = c0d0
    enddo
  endif


  ! loop branches

  if (computeNLO) then

    lp = 1

    ! Index "c" combine the indices t, fh and ih1
    cloop: do c = 1,cEffMax(pr)

      t   = cTOt  (c,pr)
      fh  = cTOfh (c,pr)
      ih1 = cTOih1(c,pr)

      locoef = loopCoef(t,pr)

      allocate (ww1(0:3,0:riwMax(c,pr),0:w1TotMax(c,pr)))

      if (ih1.eq.0) then
        emin = tlp1-1
      else
        emin = tl + (2*ih1-1)
      endif

      lmuloop: do lmu = minlmu(fh,t), maxlmu(fh,t)

        ww1(:,0,1) = c0d0
        ww1(lmu,0,1) = c1d0

        configloop1: do i = 1, cfTot(pr)

          eloop1: do e = emin,tlp1-1,2

            !if (bm1min(e,i,c,pr).ne.0) then
            if (bm1_b(pr)%conf(e,i,c)%bmin.ne.0) then

              ! bm1loop: do b = bm1min(e,i,c,pr),bm1max(e,i,c,pr)
              bm1loop: do b = bm1_b(pr)%conf(e,i,c)%bmin,bm1_b(pr)%conf(e,i,c)%bmax

                s = sbm1(b,pr)

                mo    = mosm1(s,pr)
                leg1  = binsm1(1,s,pr)
                leg2  = binsm1(2,s,pr)
                leg3  = binsm1(3,s,pr)
                leg4  = leg1 + leg2 + leg3
                gsI   = gsIncsm1(s,pr)
                riMaxIn  = riMax(rankInsm1(s,pr))
                riMaxOut = riMax(rankOutsm1(s,pr))
                cs    = cssm1(s,pr); last = (cs.ne.0)
                if (last) then
                  ferlo = ferloopsm1(s,pr)
                  gs = gssm1(s,pr)
                  ti = tism1(s,pr)
                endif

                w1in      = w1inbm1(b,pr)
                w0in(2:3) = w0inbm1(2:3,b,pr)
                w1out     = w1outbm1(b,pr)
                winit     = winitbm1(b,pr)
                ty        = typebm1(b,pr)

                if (leg3.eq.0) then
                  if (.not.loop_v(3)) cosm1(1:2,s,pr) = c0d0
                  m2 = cm2f(parsm1(s,pr))
                  call loop3 (riMaxIn,riMaxOut,p(:,leg1),p(:,leg2), &
                              pl(:,leg1),pl(:,leg2),m2,             &
                              cosm1(1:2,s,pr),ty,                   &
                              ww1(:,0:riMaxIn,w1in),ww0(:,w0in(2)), &
                              ww1out(:,0:riMaxOut,mo)               )
                else
                  if (.not.loop_v(4)) cosm1(1:3,s,pr) = c0d0
                  call loop4 (riMaxOut,cosm1(1:3,s,pr),ty,   &
                              ww1(:,0:riMaxIn,w1in),         &
                              ww0(:,w0in(2)),ww0(:,w0in(3)), &
                              ww1out(:,0:riMaxOut,mo)        )
                endif

                if (last) then
                  select case (legsti(ti,pr))
                  case (1,2);   fac = iover32pi2
                  case default; fac = 2*iover32pi2
                  end select
                  if (ferlo) fac = fac*locoef
                  wwTI = c0d0
                  do k = 0,riMaxOut
                    wwTI = wwTI + ww1out(lmu,k,mo) * TIri(k,ti)
                  enddo
                  matrix(cs,gs,i,lp,pr) = matrix(cs,gs,i,lp,pr) + fac * wwTI
                else
                  if (winit) then
                    ww1(:,riMaxOut+1:riwMax(c,pr),w1out) = c0d0
                    ww1(:,0:riMaxOut,w1out) = ww1out(:,0:riMaxOut,mo)
                  else
                    ww1(:,0:riMaxOut,w1out) = &
                    ww1(:,0:riMaxOut,w1out) + ww1out(:,0:riMaxOut,mo)
                  endif
                endif

              enddo bm1loop

            endif

            !if (bd1min(e,i,c,pr).ne.0) then
            if (bd1_b(pr)%conf(e,i,c)%bmin .ne. 0) then

              ! bd1loop: do b = bd1min(e,i,c,pr),bd1max(e,i,c,pr)
              bd1loop: do b = bd1_b(pr)%conf(e,i,c)%bmin,bd1_b(pr)%conf(e,i,c)%bmax

                s = sbd1(b,pr)

                da    = dasd1(s,pr)
                fa    = facsd1(s,pr)
                riMaxOut = riMax(rankOutsd1(s,pr))
                cs    = cssd1(s,pr); last = (cs.ne.0)
                if (last) then
                  ferlo = ferloopsd1(s,pr)
                  gs = gssd1(s,pr)
                  ti = tisd1(s,pr)
                endif

                w1out = w1outbd1(b,pr)
                winit = winitbd1(b,pr)

                if (last) then
                  select case (legsti(ti,pr))
                  case (1,2);   fac = iover32pi2
                  case default; fac = 2*iover32pi2
                  end select
                  if (ferlo) fac = fac*locoef
                  wwTI = c0d0
                  do k = 0,riMaxOut
                    wwTI = wwTI + ww1out(lmu,k,da) * TIri(k,ti)
                  enddo
                  matrix(cs,gs,i,lp,pr) = matrix(cs,gs,i,lp,pr) + fac * fa * wwTI
                else
                  if (winit) then
                    ww1(:,riMaxOut+1:riwMax(c,pr),w1out) = c0d0
                    ww1(:,0:riMaxOut,w1out) = fa * ww1out(:,0:riMaxOut,da)
                  else
                    ww1(:,0:riMaxOut,w1out) = &
                    ww1(:,0:riMaxOut,w1out) + fa * ww1out(:,0:riMaxOut,da)
                  endif
                endif

              enddo bd1loop

            endif

          enddo eloop1

        enddo configloop1

      enddo lmuloop

      deallocate (ww1)

    enddo cloop

  endif

  if (computeNLO) deallocate (ww1out,TIri)
  deallocate (ww0out,ww0)

  if (computeNLO) deallocate (psp)
  deallocate (pp,ps,pl,p)

  if (colour_optimization.ge.2.and.csTot(pr).gt.pCsTot(pr)) then
    do i = 1,cfTot(pr)
      do gs = 0,gsTot(0,pr)
        if (.not.comp0gs(gs,pr)) cycle
        do cs = pCsTot(pr)+1,csTot(pr)
          k = nIa(cs,pr)
          matrix(cs,gs,i,0,pr) =                                 &
          sum(matrix(pIa(1:k,cs,pr),gs,i,0,pr)*facIa(1:k,cs,pr))
        enddo
      enddo
    enddo
    if (computeNLO) then
      do i = 1,cfTot(pr)
        do gs = 0,gsTot(1,pr)
          if (.not.comp1gs(gs,pr)) cycle
          do cs = pCsTot(pr)+1,csTot(pr)
            k = nIa(cs,pr)
            do jj = 1,3
              matrix(cs,gs,i,jj,pr) =                                 &
              sum(matrix(pIa(1:k,cs,pr),gs,i,jj,pr)*facIa(1:k,cs,pr))
            enddo
          enddo
        enddo
      enddo
    endif
  endif

  do gs = 0,gsTot(0,pr)
    do i = 1,cfTot(pr)
      do cs = 1,csTot(pr)
        ! correction factor for LO charge coupling constant "e" in mixed schemes
        matrix(cs,gs,i,0,pr) = matrix(cs,gs,i,0,pr) * eLOfactor(gs,pr)
      enddo
    enddo
  enddo
  if (computeNLO) then
    do gs = 0,gsTot(1,pr)
      if (.not.comp1gs(gs,pr)) cycle
      do i = 1,cfTot(pr)
        do cs = 1,csTot(pr)
          ! correction factor for LO charge coupling constant "e" in mixed schemes
          matrix(cs,gs,i,1:3,pr) = matrix(cs,gs,i,1:3,pr) * eLOfactor(gs,pr)
          ! correction for dZe counterterm in mixed schemes
          matrix(cs,gs,i,2,pr) =   matrix(cs,gs,i,2,pr) &
                                 + matrix(cs,gs,i,0,pr) * DdZe(pr)
          ! correction factor for NLO charge coupling constant "e" in mixed schemes
          matrix(cs,gs,i,1:3,pr) = matrix(cs,gs,i,1:3,pr) * eNLOfactor(pr)
        enddo
      enddo
    enddo
  endif

  if (computeNLO) then
    do i = 1, cfTot(pr)
      do gs = 0,gsTot(1,pr)
        do cs = 1,csTot(pr)
          matrix(cs,gs,i,4,pr) = sum(matrix(cs,gs,i,1:3,pr))
        enddo
      enddo
    enddo
  endif

   als0R(0,pr) = als0
  Qren0R(0,pr) = Qren0
   Nlq0R(0,pr) = Nlq0
  if (computeNLO) then
     als0R(1,pr) = als0
    Qren0R(1,pr) = Qren0
     Nlq0R(1,pr) = Nlq0
    dZgs0R(pr) = dZgs0
  endif

  call cpu_time (timeTCout)
  timeTC(pr) = timeTC(pr) + timeTCout - timeTCin

  end subroutine compute_amplitude

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine rescale_amplitude (pr,order)

  integer,          intent(in) :: pr
  character(len=*), intent(in) :: order

  integer     :: i,gs,cs
  logical     :: computeNLO
  real(dp)    :: alsRatio
  complex(dp) :: fac,m2,DdZgs

  computeNLO = (order.eq.'NLO').and.(lpmax(pr).gt.0)

  ! running of alpha_s
  ! For each "n", matrix(:,n,:,:,pr) is proportional to als^(n/2).
  ! To get the running of alpha_s in matrix(:,n,:,:,pr), we first
  ! multiply matrix(:,n,:,:,pr) by sqrt(alsRatio)^n, where
  ! alsRatio = als2/als1 being als1 and als2 the old and new values
  ! of alpha_s respectively. The new amplitude matrix'(:,n,:,:,pr)
  ! is then:
  ! matrix'(:,n,:,:,pr) = sqrt(alsRatio)^n*matrix(:,n,:,:,pr).
  ! This is enough for all amplitudes, except for the CT amplitude:
  ! matrixNEW(:,n,:,i,pr) = matrix'(:,n,:,:,pr) for i = 0,1,3.
  ! For the CT amplitude we have to consider that also the
  ! counterterm dZgs of the strong coupling constant gs has to be
  ! rescaled. So we compute the new counterterm dZgs2 and we rescale
  ! the CT amplitude in the following way.
  ! The CT amplitude matrix(:,:,:,2,pr) is related to the Born
  ! amplitude matrix(:,:,:,0,pr) and dZgs through:
  ! matrix(:,n+2,:,2,pr) = matrix(:,n,:,0,pr)*( n*dZgs + REST )
  ! being n the number of gs-couplings in each diagram of the
  ! Born amplitude.
  ! The rescaled CT amplitude matrixNEW(:,n+2,:,2,pr) will be:
  ! matrixNEW(:,n+2,:,2,pr)
  ! =
  ! sqrt(alsRatio)^n*matrix(:,n,:,0,pr)*( n*dZgs2 + alsRatio*REST )
  ! =
  ! + sqrt(alsRatio)^n*matrix(:,n,:,0,pr)*
  !     ( n*alsRatio*dZgs1 + alsRatio*REST )
  ! + sqrt(alsRatio)^n*matrix(:,n,:,0,pr)*n*(dZgs2-alsRatio*dZgs1)
  ! =
  ! + sqrt(alsRatio)^(n+2)*matrix(:,n+2,:,2,pr)
  ! + matrix'(:,n,:,0,pr)*n*(dZgs2-alsRatio*dZgs1)
  ! =
  ! + matrix'(:,n+2,:,2,pr) + matrix'(:,n,:,0,pr) * n * DdZgs
  ! where
  ! DdZgs = dZgs2 - alsRatio*dZgs1
  ! For the Born amplitude we use matrixLO(:,n,:,pr) instead of
  ! matrix(:,n,:,0,pr)
  if ( ( als.ne. als0R(0,pr)).or. &
       (Qren.ne.Qren0R(0,pr)).or. &
       ( Nlq.ne. Nlq0R(0,pr))     ) then
    alsRatio = als/als0R(0,pr)
    ! Tree-level
    do gs = 1,gsTot(0,pr)
      fac = sqrt(alsRatio)**gs
      matrix(:,gs,:,0,pr) = matrix(:,gs,:,0,pr) * fac
      if (order.eq.'NLO') matrixLO(:,gs,:,pr) = matrixLO(:,gs,:,pr) * fac
    enddo
     als0R(0,pr) = als
    Qren0R(0,pr) = Qren
     Nlq0R(0,pr) = Nlq
    ! Loop-level
    if ( computeNLO .and.               &
         ( ( als.ne. als0R(1,pr)).or.   &
           (Qren.ne.Qren0R(1,pr)).or.   &
           ( Nlq.ne. Nlq0R(1,pr))     ) ) then
      alsRatio = als/als0R(1,pr)
      ! New counterterm dZgs
      dZgs = als/(4*pi)*(Nlq/3d0-11/2d0)* &
                        ( DeltaUV - log(Qren**2/muUV**2) )
      do i = Nlq+1,Nq
        if (CMscheme.eq.1) then
          m2 = mq2(i)*c1d0
        else
          m2 = real(mq2(i),kind=dp)*c1d0
        endif
        dZgs = dZgs + als/(12*pi)*( DeltaUV - log(m2/muUV**2) )
      enddo
      if (reguScheme.eq.1) then
        dZgs = dZgs + als/(24*pi)*Nc
      endif
      dZgs = real(dZgs,kind=dp)*c1d0
      ! Difference DdZgs = dZgs - alsRatio*dZgs0R(pr)
      DdZgs = dZgs - alsRatio*dZgs0R(pr)
      ! First rescaling
      do gs = 1,gsTot(1,pr)
        fac = sqrt(alsRatio)**gs
        matrix(:,gs,:,1:3,pr) = matrix(:,gs,:,1:3,pr) * fac
      enddo
      ! Second rescaling on the CT amplitude
      do gs = 1,gsTot(0,pr)
        matrix(:,gs+2,:,2,pr) = &
        matrix(:,gs+2,:,2,pr) + matrixLO(:,gs,:,pr) * gs * DdZgs
      enddo
       als0R(1,pr) = als
      Qren0R(1,pr) = Qren
       Nlq0R(1,pr) = Nlq
      dZgs0R(pr) = dZgs
    endif
  endif

  if (computeNLO) then
    do i = 1, cfTot(pr)
      do gs = 0,gsTot(1,pr)
        do cs = 1,csTot(pr)
          matrix(cs,gs,i,4,pr) = sum(matrix(cs,gs,i,1:3,pr))
        enddo
      enddo
    enddo
  endif

  end subroutine rescale_amplitude

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_squared_amplitude (pr,order)

  integer,          intent(in) :: pr
  character(len=*), intent(in) :: order

  integer              :: i,cs1,cs2,gs,gs1,gs2,j1,j2
  integer, allocatable :: cc(:,:)
  logical              :: computeNLO
  real(dp)             :: mat2(0:4),mat2int(1:3,1:3)


  computeNLO = (order.eq.'NLO').and.(lpmax(pr).gt.0)

  if (.not.allocated(matrix2h)) then
    allocate (matrix2h(0:gs2Max,1:cfMax,0:4,prTot))
    matrix2h = 0d0
  else
    matrix2h(0:gs2Tot(0,pr),1:cfTot(pr),0,pr) = 0d0
  end if

  if (.not.allocated(matrix2)) then
    allocate (matrix2(0:gs2Max,0:6,prTot))
    matrix2 = 0d0
  else
    matrix2(0:gs2Tot(0,pr),0,pr) = 0d0
  end if

  if (computeNLO.and.zeroLO(pr).and.writeMat2.ge.2) then
    if (.not.allocated(matrix2int)) then
      allocate (matrix2int(0:gs2Max,1:3,1:3,prTot))
      matrix2int = 0d0
    end if
  endif

  if (.not.prexists(pr)) return

  if (.not.momcheck) return

  allocate (cc(csTot(pr),csTot(pr)))

  do cs1 = 1,csTot(pr)
    do cs2 = cs1,csTot(pr)
      select case (cs1-cs2)
      case (0)
        cc(cs1,cs2) = colcoef(cs1,cs2,pr)
      case default
        cc(cs1,cs2) = 2 * colcoef(cs1,cs2,pr)
      end select
    enddo
  enddo

  do cs1 = 1,csTot(pr)
    do cs2 = cs1,csTot(pr)
      do gs1 = 0,gsTot(0,pr)
        if (.not.comp0gs(gs1,pr)) cycle
        do gs2 = 0,gsTot(0,pr)
          if (.not.comp0gs(gs2,pr)) cycle
          gs = gs1 + gs2
          do i = 1,cfTot(pr)
            mat2(0) = &
            cc(cs1,cs2) * real(conjg(matrix(cs1,gs1,i,0,pr))* &
                                     matrix(cs2,gs2,i,0,pr),kind=dp)
            matrix2h(gs,i,0,pr) = matrix2h(gs,i,0,pr) + mat2(0)
            matrix2 (gs,  0,pr) = matrix2 (gs,  0,pr) + mat2(0)
          enddo
        enddo
      enddo
    enddo
  enddo
  matrix2h(0:gs2Tot(0,pr),1:cfTot(pr),0,pr) = &
  matrix2h(0:gs2Tot(0,pr),1:cfTot(pr),0,pr) * factor(pr)
  matrix2(0:gs2Tot(0,pr),0,pr) = &
  matrix2(0:gs2Tot(0,pr),0,pr) * factor(pr)

  matrix2h(0:gs2Tot(1,pr),1:cfTot(pr),1:4,pr) = 0d0
  matrix2 (0:gs2Tot(1,pr),            1:4,pr) = 0d0
  if (computeNLO) then
    if (zeroLO(pr)) then
      if (writeMat2.ge.2) matrix2int(0:gs2Tot(1,pr),1:3,1:3,pr) = 0d0
      do cs1 = 1,csTot(pr)
        do cs2 = cs1,csTot(pr)
          do gs1 = 0,gsTot(1,pr)
            if (.not.comp1gs(gs1,pr)) cycle
            do gs2 = 0,gsTot(1,pr)
              if (.not.comp1gs(gs2,pr)) cycle
              gs = gs1 + gs2
              do i = 1,cfTot(pr)
                mat2(1:4) = &
                cc(cs1,cs2) * real(conjg(matrix(cs1,gs1,i,1:4,pr))* &
                                         matrix(cs2,gs2,i,1:4,pr),kind=dp)
                matrix2h(gs,i,1:4,pr) = matrix2h(gs,i,1:4,pr) + mat2(1:4)
                matrix2 (gs,  1:4,pr) = matrix2 (gs,  1:4,pr) + mat2(1:4)
                if (writeMat2.ge.2) then
                  do j1 = 1,3
                  do j2 = j1+1,3
                    mat2int(j1,j2) =                  &
                    cc(cs1,cs2) * real(               &
                    + conjg(matrix(cs1,gs1,i,j1,pr))* &
                            matrix(cs2,gs2,i,j2,pr)   &
                    + conjg(matrix(cs2,gs1,i,j1,pr))* &
                            matrix(cs1,gs2,i,j2,pr)   &
                    ,kind=dp)
                    mat2int(j2,j1) = mat2int(j1,j2)
                  enddo
                  enddo
                  matrix2int(gs,1:3,1:3,pr) = &
                  matrix2int(gs,1:3,1:3,pr) + mat2int(1:3,1:3)
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
      if (writeMat2.ge.2) &
        matrix2int(0:gs2Tot(1,pr),1:3,1:3,pr) = &
        matrix2int(0:gs2Tot(1,pr),1:3,1:3,pr) * factor(pr)
    else
      do cs1 = 1,csTot(pr)
        do cs2 = cs1,csTot(pr)
          do gs1 = 0,gsTot(0,pr)
            if (.not.comp0gs(gs1,pr)) cycle
            do gs2 = 0,gsTot(1,pr)
              if (.not.comp1gs(gs2,pr)) cycle
              do i = 1, cfTot(pr)
                gs = gs1 + gs2
                mat2(1:4) =                        &
                cc(cs1,cs2) * real(                &
                + conjg(matrix(cs1,gs1,i,0,  pr))* &
                        matrix(cs2,gs2,i,1:4,pr)   &
                + conjg(matrix(cs2,gs1,i,0,  pr))* &
                        matrix(cs1,gs2,i,1:4,pr)   &
                ,kind=dp)
                matrix2h(gs,i,1:4,pr) = &
                matrix2h(gs,i,1:4,pr) + mat2(1:4)
                matrix2 (gs,  1:4,pr) = &
                matrix2 (gs,  1:4,pr) + mat2(1:4)
              enddo
            enddo
          enddo
        enddo
      enddo
    endif
    matrix2h(0:gs2Tot(1,pr),1:cfTot(pr),1:4,pr) = &
    matrix2h(0:gs2Tot(1,pr),1:cfTot(pr),1:4,pr) * factor(pr)
    matrix2 (0:gs2Tot(1,pr),            1:4,pr) = &
    matrix2 (0:gs2Tot(1,pr),            1:4,pr) * factor(pr)
  endif

  deallocate (cc)

  end subroutine compute_squared_amplitude

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_squared_amplitude_cc (pr,i1,i2)

  integer,  intent(in) :: pr,i1,i2

  integer               :: legs,j1,j2,n,gs1,gs2,cs1,cs2
  real(dp), allocatable :: cc(:,:)
  complex(dp)           :: cmat


  legs = legsIn(pr) + legsOut(pr)

  j1 = newleg(i1,pr)
  j2 = newleg(i2,pr)
  if (.not.allocated(matrix2cc)) then
    allocate (matrix2cc(0:gs2Max,1:legsMax,1:legsMax,prTot))
    matrix2cc = 0d0
  else
    matrix2cc(0:gs2Tot(0,pr),j1,j2,pr) = 0d0
  endif

  if (.not.prexists(pr)) return

  if (.not.momcheck) return

  if (i1 .eq. i2) then
    matrix2cc(0:gs2Tot(0,pr),j1,j2,pr) = matrix2(0:gs2Tot(0,pr),0,pr)
  else
    allocate (cc(1:csTot(pr),1:csTot(pr)))

    cc(:,:) = colcoefc(1:csTot(pr),1:csTot(pr),j1,j2,pr)

    if (sum(abs(cc(:,:))).gt.zerocut) then

      do n = 1,cfTot(pr)
        do gs1 = 0,legs-2
          if (.not.comp0gs(gs1,pr)) cycle
          do cs1 = 1,csTot(pr)
            cmat = conjg(matrix(cs1,gs1,n,0,pr))
            do gs2 = 0,legs-2
              if (.not.comp0gs(gs2,pr)) cycle
              do cs2 = 1,csTot(pr)
                matrix2cc(gs1+gs2,j1,j2,pr) = &
                matrix2cc(gs1+gs2,j1,j2,pr) + &
                cc(cs1,cs2) * real(cmat*matrix(cs2,gs2,n,0,pr),kind=dp)
              enddo
            enddo
          enddo
        enddo
      enddo

      matrix2cc(0:gs2Tot(0,pr),j1,j2,pr) = &
      matrix2cc(0:gs2Tot(0,pr),j1,j2,pr) * factor(pr)

    endif

    deallocate (cc)

  endif

  end subroutine compute_squared_amplitude_cc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_squared_amplitude_cc_int (pr,i1,i2)

  integer,  intent(in) :: pr,i1,i2

  integer               :: legs,j1,j2,n,gs1,gs2,cs1,cs2,gs
  real(dp), allocatable :: cc(:,:)
  real(dp)              :: mat2


  legs = legsIn(pr) + legsOut(pr)

  j1 = newleg(i1,pr)
  j2 = newleg(i2,pr)
  if (.not.allocated(matrix2ccint)) then
    allocate (matrix2ccint(0:gs2Max,1:legsMax,1:legsMax,prTot))
    matrix2ccint = 0d0
  else
    matrix2ccint(0:gs2Tot(1,pr),j1,j2,pr) = 0d0
  endif

  if (.not.prexists(pr)) return

  if (.not.momcheck) return

  if (i1 .eq. i2) then
    matrix2ccint(0:gs2Tot(1,pr),j1,j2,pr) = matrix2(0:gs2Tot(1,pr),4,pr)
  else
    allocate (cc(1:csTot(pr),1:csTot(pr)))

    cc(:,:) = colcoefc(1:csTot(pr),1:csTot(pr),j1,j2,pr)

    if (sum(abs(cc(:,:))).gt.zerocut) then

      do n = 1,cfTot(pr)
        do gs1 = 0,gsTot(0,pr)
          if (.not.comp0gs(gs1,pr)) cycle
          do cs1 = 1,csTot(pr)
            do gs2 = 0,gsTot(1,pr)
              if (.not.comp1gs(gs2,pr)) cycle
              gs = gs1 + gs2
              do cs2 = 1,csTot(pr)
                mat2 =                             &
                cc(cs1,cs2) * real(                &
                + conjg(matrix(cs1,gs1,n,0,pr))* &
                        matrix(cs2,gs2,n,4,pr)   &
                + conjg(matrix(cs2,gs1,n,0,pr))* &
                        matrix(cs1,gs2,n,4,pr)   &
                  ,kind=dp)
                matrix2ccint(gs1+gs2,j1,j2,pr) =   &
                matrix2ccint(gs1+gs2,j1,j2,pr) + mat2
              enddo
            enddo
          enddo
        enddo
      enddo

      matrix2ccint(0:gs2Tot(1,pr),j1,j2,pr) = &
      matrix2ccint(0:gs2Tot(1,pr),j1,j2,pr) * factor(pr)

    endif
    deallocate (cc)
  endif

  end subroutine compute_squared_amplitude_cc_int

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_squared_amplitude_cc_nlo (pr,i1,i2)

  integer,  intent(in) :: pr,i1,i2

  integer               :: j1,j2,n,gs1,gs2,cs1,cs2
  real(dp), allocatable :: cc(:,:)
  complex(dp)           :: cmat

  j1 = newleg(i1,pr)
  j2 = newleg(i2,pr)
  if (.not.allocated(matrix2ccnlo)) then
    allocate (matrix2ccnlo(0:gs2Max,1:legsMax,1:legsMax,prTot))
    matrix2ccnlo = 0d0
  else
    matrix2ccnlo(0:gs2Tot(1,pr),j1,j2,pr) = 0d0
  endif

  if (.not.prexists(pr)) return

  if (.not.momcheck) return

  if (i1 .eq. i2) then
    matrix2ccnlo(0:gs2Tot(1,pr),j1,j2,pr) = matrix2(0:gs2Tot(1,pr),4,pr)
  else
    allocate (cc(1:csTot(pr),1:csTot(pr)))

    cc(:,:) = colcoefc(1:csTot(pr),1:csTot(pr),j1,j2,pr)

    if (sum(abs(cc(:,:))).gt.zerocut) then

      do n = 1,cfTot(pr)
        do gs1 = 0,gsTot(1,pr)
          if (.not.comp1gs(gs1,pr)) cycle
          do cs1 = 1,csTot(pr)
            cmat = conjg(matrix(cs1,gs1,n,4,pr))
            do gs2 = 0,gsTot(1,pr)
              if (.not.comp1gs(gs2,pr)) cycle
              do cs2 = 1,csTot(pr)
                matrix2ccnlo(gs1+gs2,j1,j2,pr) = &
                matrix2ccnlo(gs1+gs2,j1,j2,pr) + &
                cc(cs1,cs2) * real(cmat*matrix(cs2,gs2,n,4,pr),kind=dp)
              enddo
            enddo
          enddo
        enddo
      enddo

      matrix2ccnlo(0:gs2Tot(1,pr),j1,j2,pr) = &
      matrix2ccnlo(0:gs2Tot(1,pr),j1,j2,pr) * factor(pr)

    endif

    deallocate (cc)

  endif

  end subroutine compute_squared_amplitude_cc_nlo

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_squared_amplitude_scc (pr,i1,i2,v)

  integer,     intent(in) :: pr,i1,i2
  complex(dp), intent(in) :: v(0:)

  integer                  :: legs,j1,j2,par1,h,n,dn,gs1,gs2,cs1,cs2
  real(dp), allocatable    :: cc(:,:)
  complex(dp)              :: p1(0:3),pl1(1:4),w01(0:3,-1:1),  &
                              a(-1:1),cmatsc
  complex(dp), allocatable :: matsc(:,:,:)


  legs = legsIn(pr) + legsOut(pr)

  j1 = newleg(i1,pr)
  j2 = newleg(i2,pr)
  if (.not.allocated(matrix2scc)) then
    allocate (matrix2scc(0:gs2Max,1:legsMax,1:legsMax,prTot))
    matrix2scc = 0d0
  else
    matrix2scc(0:gs2Tot(0,pr),j1,j2,pr) = 0d0
  endif

  if (.not.prexists(pr)) return

  if (.not.momcheck) return

  allocate (cc(1:csTot(pr),1:csTot(pr)))

  cc(:,:) = colcoefc(1:csTot(pr),1:csTot(pr),j1,j2,pr)

  par1 = par(i1,pr)

  if (par1.eq.15.and.sum(abs(cc(:,:))).gt.zerocut) then

    allocate (matsc(csTot(pr),0:gsTot(0,pr),cfTot(pr)))

    p1 = cmplx(momenta(:,i1),kind=dp)
    if (i1.gt.legsIn(pr)) p1 = - p1
    pl1(1) = p1(0) + p1(3)
    pl1(2) = p1(0) - p1(3)
    pl1(3) = p1(1) + cId0*p1(2)
    pl1(4) = conjg(pl1(3))

    do h = -1,1,2
      call definewp (par1,p1,pl1,0d0,h,w01(:,h))
      a(h) = + v(0) * conjg(w01(0,h)) &
             - v(1) * conjg(w01(1,h)) &
             - v(2) * conjg(w01(2,h)) &
             - v(3) * conjg(w01(3,h))
    enddo

    do n = 1,cfTot(pr)
      if (heli(j1,n,pr).ne.-1) cycle
      dn  = dualheli(j1,n,pr)
      matsc(:,:,n) = + a(-1)*matrix(1:csTot(pr),0:gsTot(0,pr), n,0,pr) &
                     + a(+1)*matrix(1:csTot(pr),0:gsTot(0,pr),dn,0,pr)
      do gs1 = 0,gsTot(0,pr)
        if (.not.comp0gs(gs1,pr)) cycle
        do cs1 = 1,csTot(pr)
          cmatsc = conjg(matsc(cs1,gs1,n))
          do gs2 = 0,gsTot(0,pr)
            if (.not.comp0gs(gs2,pr)) cycle
            do cs2 = 1,csTot(pr)
              matrix2scc(gs1+gs2,j1,j2,pr) = &
              matrix2scc(gs1+gs2,j1,j2,pr) + &
              cc(cs1,cs2) * real(cmatsc*matsc(cs2,gs2,n),kind=dp)
            enddo
          enddo
        enddo
      enddo
    enddo

    matrix2scc(0:gs2Tot(0,pr),j1,j2,pr) = &
    matrix2scc(0:gs2Tot(0,pr),j1,j2,pr) * factor(pr)

  endif

  deallocate (cc)

  end subroutine compute_squared_amplitude_scc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_squared_amplitude_sc (pr,j,v)

  integer,     intent(in) :: pr,j
  complex(dp), intent(in) :: v(0:)

  integer                  :: legs,jj,parj,h,n,dn,gs1,gs2, &
                              cs1,cs2
  real(dp), allocatable    :: cc(:,:)
  complex(dp)              :: pj(0:3),plj(1:4),w0j(0:3,-1:1),  &
                              a(-1:1),cmats
  complex(dp), allocatable :: mats(:,:,:)


  legs = legsIn(pr) + legsOut(pr)

  if (.not.allocated(matrix2sc)) then
    allocate (matrix2sc(0:gs2Max,prTot))
    matrix2sc = 0d0
  else
    matrix2sc(0:gs2Tot(0,pr),pr) = 0d0
  endif

  if (.not.prexists(pr)) return

  if (.not.momcheck) return

  jj = newleg(j,pr)

  parj = par(j,pr)

  if (parj.eq.15.or.parj.eq.16) then

    allocate (cc(csTot(pr),csTot(pr)))

    do cs1 = 1,csTot(pr)
      do cs2 = cs1,csTot(pr)
        select case (cs1-cs2)
        case (0)
          cc(cs1,cs2) = colcoef(cs1,cs2,pr)
        case default
          cc(cs1,cs2) = 2 * colcoef(cs1,cs2,pr)
        end select
      enddo
    enddo

    allocate (mats(csTot(pr),0:gsTot(0,pr),cfTot(pr)))

    pj = cmplx(momenta(:,j),kind=dp)
    if (j.gt.legsIn(pr)) pj = - pj
    plj(1) = pj(0) + pj(3)
    plj(2) = pj(0) - pj(3)
    plj(3) = pj(1) + cId0*pj(2)
    plj(4) = conjg(plj(3))

    do h = -1,1,2
      call definewp (parj,pj,plj,0d0,h,w0j(:,h))
      a(h) = + v(0) * conjg(w0j(0,h)) &
             - v(1) * conjg(w0j(1,h)) &
             - v(2) * conjg(w0j(2,h)) &
             - v(3) * conjg(w0j(3,h))
    enddo

    do n = 1,cfTot(pr)
      if (heli(jj,n,pr).ne.-1) cycle
      dn  = dualheli(jj,n,pr)
      mats(:,:,n) =                                     &
      + a(-1)*matrix(1:csTot(pr),0:gsTot(0,pr), n,0,pr) &
      + a(+1)*matrix(1:csTot(pr),0:gsTot(0,pr),dn,0,pr)
      do gs1 = 0,gsTot(0,pr)
        if (.not.comp0gs(gs1,pr)) cycle
        do cs1 = 1,csTot(pr)
          cmats = conjg(mats(cs1,gs1,n))
          do gs2 = 0,gsTot(0,pr)
            if (.not.comp0gs(gs2,pr)) cycle
            do cs2 = cs1,csTot(pr)
              matrix2sc(gs1+gs2,pr) = &
              matrix2sc(gs1+gs2,pr) + &
              cc(cs1,cs2) * real(cmats*mats(cs2,gs2,n),kind=dp)
            enddo
          enddo
        enddo
      enddo
    enddo

    matrix2sc(0:gs2Tot(0,pr),pr) =            &
    matrix2sc(0:gs2Tot(0,pr),pr) * factor(pr)

  endif

  deallocate (cc)

  end subroutine compute_squared_amplitude_sc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_squared_amplitude_scm (pr,j)

  integer,     intent(in) :: pr,j

  integer                  :: legs,jj,parj,h,n,dn, &
                              cs1,cs2,gs1,gs2, &
                              mu,nu
  real(dp), allocatable    :: cc(:,:)
  complex(dp)              :: pj(0:3),plj(1:4),w0j(0:3,-1:1),  &
                              a(-1:1,0:3)
  legs = legsIn(pr) + legsOut(pr)

  if (.not.allocated(matrix2scm)) then
    allocate (matrix2scm(0:3,0:3,0:gs2Max,prTot))
    matrix2scm = 0d0
  else
    matrix2scm(0:3,0:3,0:gs2Tot(0,pr),pr) = 0d0
  endif

  if (.not.prexists(pr)) return

  if (.not.momcheck) return

  jj = newleg(j,pr)

  parj = par(j,pr)

  if (parj.eq.15.or.parj.eq.16) then

    allocate (cc(csTot(pr),csTot(pr)))

    do cs1 = 1,csTot(pr)
      do cs2 = cs1,csTot(pr)
        select case (cs1-cs2)
        case (0)
          cc(cs1,cs2) = colcoef(cs1,cs2,pr)
        case default
          cc(cs1,cs2) = 2 * colcoef(cs1,cs2,pr)
        end select
      enddo
    enddo

    pj = cmplx(momenta(:,j),kind=dp)
    if (j .gt. legsIn(pr)) pj = - pj
    plj(1) = pj(0) + pj(3)
    plj(2) = pj(0) - pj(3)
    plj(3) = pj(1) + cId0*pj(2)
    plj(4) = conjg(plj(3))

    do h = -1,1,2
      call definewp (parj,pj,plj,0d0,h,w0j(:,h))
      a(h,0) =  conjg(w0j(0,h))
      a(h,1) = -conjg(w0j(1,h))
      a(h,2) = -conjg(w0j(2,h))
      a(h,3) = -conjg(w0j(3,h))
    end do

    do n = 1,cfTot(pr)
      if (heli(jj,n,pr).ne.-1) cycle
      dn  = dualheli(jj,n,pr)
      do gs1 = 0,gsTot(0,pr)
        if (.not.comp0gs(gs1,pr)) cycle
        do cs1 = 1,csTot(pr)
          do gs2 = 0,gsTot(0,pr)
            if (.not.comp0gs(gs2,pr)) cycle
            do cs2 = cs1,csTot(pr)
              do mu = 1, 3
                do nu = 1, 3
                  matrix2scm(mu,nu,gs1+gs2,pr) = &
                  matrix2scm(mu,nu,gs1+gs2,pr) + &
                  cc(cs1,cs2) * real( &
                  conjg(a(-1,mu)*matrix(cs1,gs1, n,0,pr) +  &
                        a(+1,mu)*matrix(cs1,gs1,dn,0,pr)) * &
                       (a(-1,nu)*matrix(cs2,gs2, n,0,pr) +  &
                        a(+1,nu)*matrix(cs2,gs2,dn,0,pr)),  &
                  kind=dp)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    matrix2scm(:,:,0:gs2Tot(0,pr),pr) =            &
    matrix2scm(:,:,0:gs2Tot(0,pr),pr) * factor(pr)

  endif

  deallocate (cc)

  end subroutine compute_squared_amplitude_scm


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine print_process_and_momenta (pr)

  integer, intent (in) :: pr

  integer     :: i,j,legs,e1,nal
  real(dp)    :: p2,pm,rm2,rm,w
  complex(dp) :: cm2
  character   :: calpha(1:4)*15,alphaN*80,cn*2,      &
                 fmt(0:3)*6,fmtTot*80,ci*2,cpa*8,    &
                 fmtp2*80,fmtpm*80,fmtm2*80,fmtm*80

  legs = legsIn(pr) + legsOut(pr)

  call openOutput

  write(nx,*)
  write(nx,'(1x,75("x"))')
  write(nx,*)
  write(nx,*)
  write(nx,*) ' ',trim(process(pr))
  write(nx,*)

  nal = sum(Nalpha(:,pr))
  calpha(1) = 'alpha_Gf'
  calpha(2) = 'alpha(0)'
  calpha(3) = 'alpha(M_Z)'
  calpha(4) = 'alphaMSbar'
  if ( nal.eq.0 ) then
    alphaN = trim(adjustl(calpha(refscheme(pr))))//'^N'
  else
    write(cn,'(i2)') nal
    alphaN = trim(adjustl(calpha(refscheme(pr))))//'^(N-'//trim(adjustl(cn))//')'
  endif
  if (Nalpha(1,pr).gt.0) then
    write(cn,'(i2)') Nalpha(1,pr)
    alphaN = trim(alphaN)//' * alpha_Gf^'//trim(adjustl(cn))
  endif
  if (Nalpha(2,pr).gt.0) then
    write(cn,'(i2)') Nalpha(2,pr)
    alphaN = trim(alphaN)//' * alpha(0)^'//trim(adjustl(cn))
  endif
  if (Nalpha(3,pr).gt.0) then
    write(cn,'(i2)') Nalpha(3,pr)
    alphaN = trim(alphaN)//' * alpha(M_Z)^'//trim(adjustl(cn))
  endif
  if (Nalpha(4,pr).gt.0) then
    write(cn,'(i2)') Nalpha(4,pr)
    alphaN = trim(alphaN)//' * alphaMSbar^'//trim(adjustl(cn))
  endif
  write(nx,'(2x,2a)') 'alpha scheme: alpha_Born^N  = ',trim(adjustl(alphaN))
  select case (NLOscheme(pr))
  case (1); write(nx,'(16x,a)') 'alpha_NLO     = alpha_Gf'
  case (2); write(nx,'(16x,a)') 'alpha_NLO     = alpha(0)'
  case (3); write(nx,'(16x,a)') 'alpha_NLO     = alpha(M_Z)'
  case (4); write(nx,'(16x,a)') 'alpha_NLO     = alphaMSbar'
  end select
  do i = 1,4
    if ( refscheme(pr).eq.i .or. &
         NLOscheme(pr).eq.i .or. &
         Nalpha(i,pr).gt.0 ) then
      write(nx,'(4x,a10,a,g21.14)') adjustl(calpha(i)),' =',alphai(i,pr)
    endif
  enddo
  write(nx,*)

  do i = 1,legs
    write(ci,'(i2)') i
    fmt(0) = 'f14.9'
    if (abs(momenta(0,i)).ge.1d4) fmt(0) = 'e14.8'
    if (legs.lt.10) then
      fmtTot = '(2x,a,a,1x,"(",'//fmt(0)
    else
      fmtTot = '(2x,a,1x,a,1x,"(",'//fmt(0)
    endif
    do j = 1,3
      fmt(j) = 'f15.9'
      if (abs(momenta(j,i)).ge.1d4) fmt(j) = 'e15.8'
      fmtTot = trim(fmtTot)//',",",'//fmt(j)
    enddo
    fmtTot = trim(fmtTot)//',") GeV")'
    write(nx,trim(fmtTot)) 'p'//trim(adjustl(ci)),'=',momenta(0:3,i)
  enddo
  write(nx,*)

  do i = 1,resMax(pr)
    write(ci,'(i2)') i
    cpa = cpar(anti(parRes(i,pr)))
    e1 = newbin(binRes(i,pr),pr)
    if (defp2bin(e1,pr)) then
      p2 = p2bin(e1,pr)
      pm = sqrt(p2)
      fmtp2 = 'f14.9'; if (p2.ge.1d4) fmtp2 = 'e14.8'
      fmtpm = 'f14.9'; if (pm.ge.1d4) fmtpm = 'e14.8'
      write(nx,'(2x,5a)') &
        'The denominator of the propagator of resonance number ', &
        trim(adjustl(ci)),' (',trim(cpa),' particle) '
      write(nx,'(2x,a)') 'carries an off-sheel squared momentum p^2:'
      write(nx,'(4x,a,'//trim(fmtp2)//',a)') 'p^2       = ',p2,' GeV^2'
      write(nx,'(4x,a,'//trim(fmtpm)//',a)') 'sqrt(p^2) = ',pm,' GeV'
      write(nx,*)
    endif
    if (defresbin(e1,pr)) then
      cm2 = cm2pf(parRes(i,pr))
      rm2 = real(cm2,kind=dp)
      rm = sqrt(rm2)
      w = - real(aimag(cm2),kind=dp)/rm
      fmtm2 = '(4x,a,"(",e14.8,",",e15.8,")",a)'
      fmtm  = '(4x,a,g21.14,6x,a,g21.14)'
      write(nx,'(2x,5a)') &
        'The denominator of the propagator of resonance number ', &
        trim(adjustl(ci)),' (',trim(cpa),' particle) '
      write(nx,'(2x,a)') 'carries a complex mass m^2:'
      write(nx,fmtm2) 'm^2  = ',cm2,' GeV^2'
      write(nx,fmtm)  'mass =',rm,'width =',w
      write(nx,*)
    endif
  enddo

  end subroutine print_process_and_momenta

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine print_rescaling

  call openOutput

  if ((als.ne.als0).or.(Qren.ne.Qren0).or.(Nfren.ne.Nfren0)) then
    if (Nfren.eq.-1) then
      write(nx,'(2x,a)') &
      'alpha_s Renormalization Scheme: Variable flavours Scheme'
    else
      write(nx,'(2x,a,i1,a)') &
      'alpha_s Renormalization Scheme: ',Nfren,'-flavours Scheme'
    endif
    write(nx,'(2x,a,g21.14,7x,a,g21.14,a)') &
    'alpha_s(Q) =',als,'Q =',Qren,' GeV'
    write(nx,*)
  endif
  write(nx,*)

  end subroutine print_rescaling

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine print_parameters_change

  call openOutput

  if (changed_DeltaUV) then
    write(nx,'(2x,a,g21.14)') 'Delta_UV = ',DeltaUV
    write(nx,*)
  endif

  if (changed_muUV) then
    write(nx,'(2x,a,g21.14,a)') 'mu_UV = ',muUV,' GeV'
    write(nx,*)
  endif

  if (changed_DeltaIR) then
    write(nx,'(2x,a,g21.14,7x,a,g21.14)') &
      'Delta_IR^2 =',DeltaIR2,'Delta_IR = ',DeltaIR
    write(nx,*)
  endif

  if (changed_muIR) then
    write(nx,'(2x,a,g21.14,a)') 'mu_IR = ',muIR,' GeV'
    write(nx,*)
  endif

  if (changed_lambda) then
    write(nx,'(2x,a,g21.14,a)') 'Mass regulator = ',lambda,' GeV'
    write(nx,*)
  endif

  end subroutine print_parameters_change

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine print_amplitude (pr,order)

  integer,          intent(in) :: pr
  character(len=*), intent(in) :: order

  integer                :: i,j,k,l,n,legs,cs,gs,lp
  logical                :: computeNLO
  real(dp)               :: x
  character(99)          :: hconf1,hconf2,fmt
  character, allocatable :: delta(:)*6,up(:)*6,lo(:)*6

  computeNLO = (order.eq.'NLO').and.(lpmax(pr).gt.0)

  legs = legsIn(pr) + legsOut(pr)

  n = 4 - legs

  call openOutput

  h1loop: do i = 1, cfTot(pr)

    x = sum(abs(matrix(1:csTot(pr),0:gsTot(0,pr),i,0,pr)))
    if (computeNLO) &
      x = x + sum(abs(matrix(1:csTot(pr),0:gsTot(1,pr),i,4,pr)))

    if (x.eq.0d0) cycle

    select case (n)
    case (-99:-10); fmt = '(1x,a,i3,a)'
    case (-9:-1);   fmt = '(1x,a,i2,a)'
    case (0:9);     fmt = '(1x,a,i1,a)'
    case default
    end select

    write(nx,'(1x,75("-"))')
    write(nx,*)
    write(nx,*)
    write(nx,trim(fmt)) ' AMPLITUDE [GeV^',n,']'
    write(nx,*)
    write(nx,*)

    hconf1 = ''
    do j = 1,legsIn(pr)
      hconf1 = trim(hconf1)//' '//trim(cpar(par(j,pr)))
      if     (heli(newleg(j,pr),i,pr).eq.-1) then; hconf1 = trim(hconf1)//'[-]'
      elseif (heli(newleg(j,pr),i,pr).eq. 0) then; hconf1 = trim(hconf1)//'[0]'
      elseif (heli(newleg(j,pr),i,pr).eq.+1) then; hconf1 = trim(hconf1)//'[+]'
      endif
    enddo

    hconf2 = ''
    do j = legsIn(pr)+1,legs
      hconf2 = trim(hconf2)//' '//trim(cpar(anti(par(j,pr))))
      if     (heli(newleg(j,pr),i,pr).eq.-1) then; hconf2 = trim(hconf2)//'[+]'
      elseif (heli(newleg(j,pr),i,pr).eq. 0) then; hconf2 = trim(hconf2)//'[0]'
      elseif (heli(newleg(j,pr),i,pr).eq.+1) then; hconf2 = trim(hconf2)//'[-]'
      endif
    enddo

    l = len(trim(hconf1)//' ->'//trim(hconf2))

    if (l.gt.75) then
      write(nx,*) ' Helicity configuration:'
      write(nx,*) trim(hconf1),' ->'
      write(nx,*) trim(hconf2)
    elseif (l.gt.75-27) then
      write(nx,*) ' Helicity configuration:'
      write(nx,*) trim(hconf1),' ->',trim(hconf2)
    else
      write(nx,*) ' Helicity configuration:  ', &
                  trim(hconf1),' ->',trim(hconf2)
    endif
    write(nx,*)

    csloop: do cs = 1,csTot(pr)

      if (   sum(abs(matrix(cs,0:gsTot(0,pr),i,0,pr))) &
           + sum(abs(matrix(cs,0:gsTot(1,pr),i,4,pr))) &
           .eq.0d0                                     ) cycle

      allocate (delta(legs),up(legs),lo(legs))
      k = 0
      do j = 1,legs
        ! upper index is "ia", lower index is "iq"
        if (csIq(j,cs,pr).ne.0) then
          k = k + 1
          delta(k) = '  d'
          write(up(k),'(i2)') oldleg(j,pr)
          write(lo(k),'(i2)') oldleg(csIq(j,cs,pr),pr)
          up(k) = '   i'//trim(adjustl(up(k)))
          lo(k) = '   j'//trim(adjustl(lo(k)))
        endif
      enddo
      if (k.gt.0) then
        write(nx,*) '                   ',up(1:k)
        write(nx,*) ' Colour structure: ',delta(1:k)
        write(nx,*) '                   ',lo(1:k)
        write(nx,*)
      endif
      deallocate (delta,up,lo)

      if (sum(abs(matrix(cs,0:gsTot(0,pr),i,0,pr))).ne.0d0) then

        write(nx,*) '  gs |               Born Amplitude A0               '
        write(nx,*) ' ----------------------------------------------------'
        do gs = 0,gsTotEff(0,pr)
          write(nx,'(4x,i1," | (",e21.14,",",e21.14,")")') &
                   gs,matrix(cs,gs,i,0,pr)
        enddo
        write(nx,*) ' ----------------------------------------------------'
        write(nx,'( "  SUM | (",e21.14,",",e21.14,")")') &
                   sum(matrix(cs,0:gsTot(0,pr),i,0,pr))
        write(nx,*)

      endif

      if ( computeNLO .and.                                   &
           (sum(abs(matrix(cs,0:gsTot(1,pr),i,4,pr))).ne.0d0) ) then

        write(nx,*) '  gs |              1-loop Amplitude A1              '
        write(nx,*) ' ----------------------------------------------------'
        do gs = 0,gsTotEff(1,pr)
          write(nx,'(4x,i1," | (",e21.14,",",e21.14,")")') &
                   gs,matrix(cs,gs,i,4,pr)
        enddo
        write(nx,*) ' ----------------------------------------------------'
        write(nx,'( "  SUM | (",e21.14,",",e21.14,")")') &
                   sum(matrix(cs,0:gsTot(1,pr),i,4,pr))
        write(nx,*)

        if (writeMat.ge.2) then
          do lp = 1,3
            select case (lp)
            case (1)
              write(nx,*) '  gs |    4-dimensional bare-loop Amplitude A1d4     '
            case (2)
              write(nx,*) '  gs |               CT Amplitude A1ct               '
            case (3)
              write(nx,*) '  gs |               R2 Amplitude A1r2               '
            end select
            write(nx,*)   ' ----------------------------------------------------'
            do gs = 0,gsTotEff(1,pr)
              write(nx,'(4x,i1," | (",e21.14,",",e21.14,")")') &
                         gs,matrix(cs,gs,i,lp,pr)
            enddo
            write(nx,*)  ' ----------------------------------------------------'
            write(nx,'( "  SUM | (",e21.14,",",e21.14,")")') &
                         sum(matrix(cs,0:gsTot(1,pr),i,lp,pr))
            write(nx,*)
          enddo
        endif

      endif

      write(nx,*)

    enddo csloop

  enddo h1loop

  end subroutine print_amplitude

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine print_squared_amplitude (pr,order)

  integer,          intent(in) :: pr
  character(len=*), intent(in) :: order

  integer       :: i,j,l,legs,n,gs
  logical       :: computeNLO
  real(dp)      :: x,y
  character(99) :: hconf1,hconf2,fmt

  computeNLO = (order.eq.'NLO').and.(lpmax(pr).gt.0)

  legs = legsIn(pr) + legsOut(pr)

  n = 8 - 2*legs

  x = sum(abs(matrix2(0:gs2Tot(0,pr),0,pr)))
  if (computeNLO) x = x + sum(abs(matrix2(0:gs2Tot(1,pr),4,pr)))

  if (x.ne.0d0) then

    call openOutput

    if ( writeMat2.ge.3 ) then

      write(nx,'(1x,75("-"))')
      write(nx,*)
      write(nx,*)

      h2loop: do i = 1,cfTot(pr)

        y = sum(abs(matrix2h(0:gs2Tot(0,pr),i,0,pr)))
        if (computeNLO) &
          y = y + sum(abs(matrix2h(0:gs2Tot(1,pr),i,4,pr)))

        if (y.eq.0d0) cycle

        select case (n)
        case (-99:-10); fmt = '(1x,a,i3,a)'
        case (-9:-1);   fmt = '(1x,a,i2,a)'
        case (0:9);     fmt = '(1x,a,i1,a)'
        case default
        end select

        write(nx,*)
        write(nx,trim(fmt)) ' POLARIZED SQUARED AMPLITUDE [GeV^',n,']'
        write(nx,*)

        hconf1 = ''
        do j = 1,legsIn(pr)
          hconf1 = trim(hconf1)//' '//trim(cpar(par(j,pr)))
          if     (heli(newleg(j,pr),i,pr).eq.-1) then; hconf1 = trim(hconf1)//'[-]'
          elseif (heli(newleg(j,pr),i,pr).eq. 0) then; hconf1 = trim(hconf1)//'[0]'
          elseif (heli(newleg(j,pr),i,pr).eq.+1) then; hconf1 = trim(hconf1)//'[+]'
          endif
        enddo

        hconf2 = ''
        do j = legsIn(pr)+1,legs
          hconf2 = trim(hconf2)//' '//trim(cpar(anti(par(j,pr))))
          if     (heli(newleg(j,pr),i,pr).eq.-1) then; hconf2 = trim(hconf2)//'[+]'
          elseif (heli(newleg(j,pr),i,pr).eq. 0) then; hconf2 = trim(hconf2)//'[0]'
          elseif (heli(newleg(j,pr),i,pr).eq.+1) then; hconf2 = trim(hconf2)//'[-]'
          endif
        enddo

        l = len(trim(hconf1)//' ->'//trim(hconf2))

        if (l.gt.75) then
          write(nx,*) ' Helicity configuration:'
          write(nx,*) trim(hconf1),' ->'
          write(nx,*) trim(hconf2)
        elseif (l.gt.75-27) then
          write(nx,*) ' Helicity configuration:'
          write(nx,*) trim(hconf1),' ->',trim(hconf2)
        else
          write(nx,*) ' Helicity configuration:  ', &
                      trim(hconf1),' ->',trim(hconf2)
        endif
        write(nx,*)

        if ( computeNLO .and.                                 &
             sum(abs(matrix2h(0:gs2Tot(1,pr),i,4,pr))).ne.0d0 ) then

          if (zeroLO(pr)) then
            write(nx,*) '  als |       | A0h |^2                 ', &
                        '  als |       | A1h |^2       '
          else
            write(nx,*) '  als |       | A0h |^2                 ', &
                        '  als |  2*Re{ A1h * A0h^* }  '
          endif
          write(nx,*) ' -----------------------------          ', &
                      ' -----------------------------'
          do gs = 0,gs2TotEff(0,pr),2
            write(nx,'(4x,i2," | ",e21.14,12x,2x,i2," | ",e21.14)') &
                     gs/2,matrix2h(gs,i,0,pr),gs/2,matrix2h(gs,i,4,pr)
          enddo
          do gs = gs2TotEff(0,pr)+2,gs2TotEff(1,pr),2
            write(nx,'(4x,2x," | ",   21x,12x,2x,i2," | ",e21.14)') &
                     gs/2,matrix2h(gs,i,4,pr)
          enddo
          write(nx,*) ' -----------------------------          ', &
                      ' -----------------------------'
          write(nx,'( "   SUM | ",e21.14,12x,   "SUM  | ",e21.14)') &
                      sum(matrix2h(0:gs2Tot(0,pr),i,0,pr)),         &
                      sum(matrix2h(0:gs2Tot(1,pr),i,4,pr))
          write(nx,*)

        elseif (sum(abs(matrix2h(0:gs2Tot(0,pr),i,0,pr))).ne.0d0) then

          write(nx,*) '  als |       | A0h |^2       '
          write(nx,*) ' -----------------------------'
          do gs = 0,gs2TotEff(0,pr),2
            write(nx,'(4x,i2," | ",e21.14)') gs/2,matrix2h(gs,i,0,pr)
          enddo
          write(nx,*) ' -----------------------------'
          write(nx,'( "   SUM | ",e21.14)') sum(matrix2h(0:gs2Tot(0,pr),i,0,pr))
          write(nx,*)

        endif

        write(nx,*)

      enddo h2loop

    endif

    select case (n)
    case (-99:-10); fmt = '(1x,a,i3,a)'
    case (-9:-1);   fmt = '(1x,a,i2,a)'
    case (0:9);     fmt = '(1x,a,i1,a)'
    case default
    end select

    write(nx,'(1x,75("-"))')
    write(nx,*)
    write(nx,*)

    write(nx,trim(fmt)) ' UNPOLARIZED SQUARED AMPLITUDE [GeV^',n,']'
    write(nx,*)

    if (.not.computeNLO) then

      write(nx,*) '  als |        | A0 |^2       '
      write(nx,*) ' -----------------------------'
      do gs = 0,gs2TotEff(0,pr),2
        write(nx,'(4x,i2," | ",e21.14)') gs/2,matrix2(gs,0,pr)
      enddo
      write(nx,*) ' -----------------------------'
      write(nx,'( "   SUM | ",e21.14)') sum(matrix2(0:gs2Tot(0,pr),0,pr))
      write(nx,*)

    else

      if (zeroLO(pr)) then
        write(nx,*) '  als |        | A0 |^2                 ', &
                    '  als |        | A1 |^2      '
      else
        write(nx,*) '  als |        | A0 |^2                 ', &
                    '  als |   2*Re{ A1 * A0^* }   '
      endif
      write(nx,*) ' -----------------------------          ', &
                  ' -----------------------------'
      do gs = 0,gs2TotEff(0,pr),2
        write(nx,'(4x,i2," | ",e21.14,12x,2x,i2," | ",e21.14)') &
                 gs/2,matrix2(gs,0,pr),gs/2,matrix2(gs,4,pr)
      enddo
      do gs = gs2TotEff(0,pr)+2,gs2TotEff(1,pr),2
        write(nx,'(4x,2x," | ",   21x,12x,2x,i2," | ",e21.14)') &
                 gs/2,matrix2(gs,4,pr)
      enddo
      write(nx,*) ' -----------------------------          ', &
                  ' -----------------------------'
      write(nx,'( "   SUM | ",e21.14,12x,   " SUM | ",e21.14)') &
                 sum(matrix2(0:gs2Tot(0,pr),0,pr)),             &
                 sum(matrix2(0:gs2Tot(1,pr),4,pr))
      write(nx,*)

      if (writeMat2.ge.2) then
        if (zeroLO(pr)) then
          write(nx,*) '  als |     | A1d4 |^2     ', &
                            '|     | A1ct |^2     ', &
                            '|     | A1r2 |^2     '
        else
          write(nx,*) '  als |   2*Re{A1d4*A0^*}  ', &
                            '|   2*Re{A1ct*A0^*}  ', &
                            '|   2*Re{A1r2*A0^*}  '
        endif
        write(nx,*) ' --------------------------', &
                          '---------------------', &
                          '---------------------'
        do gs = 0,gs2TotEff(1,pr),2
          write(nx,'(4x,i2," | ",e18.11," | ",e18.11," | ",e18.11)') &
                   gs/2,matrix2(gs,1,pr),matrix2(gs,2,pr),           &
                        matrix2(gs,3,pr)
        enddo
        write(nx,*) ' --------------------------', &
                          '---------------------', &
                          '---------------------'
        write(nx,'( "   SUM | ",e18.11," | ",e18.11," | ",e18.11)') &
                 sum(matrix2(0:gs2Tot(1,pr),1,pr)),                 &
                 sum(matrix2(0:gs2Tot(1,pr),2,pr)),                 &
                 sum(matrix2(0:gs2Tot(1,pr),3,pr))
        write(nx,*)
        if (zeroLO(pr)) then
          write(nx,*) '  als |  2*Re{A1d4*A1ct^*} ', &
                            '|  2*Re{A1ct*A1r2^*} ', &
                            '|  2*Re{A1r2*A1d4^*} '
          write(nx,*) ' --------------------------', &
                            '---------------------', &
                            '---------------------'
          do gs = 0,gs2TotEff(1,pr),2
            write(nx,'(4x,i2," | ",e18.11," | ",e18.11," | ",e18.11)') &
                     gs/2,matrix2int(gs,1,2,pr),matrix2int(gs,2,3,pr), &
                          matrix2int(gs,3,1,pr)
          enddo
          write(nx,*) ' --------------------------', &
                            '---------------------', &
                            '---------------------'
          write(nx,'( "   SUM | ",e18.11," | ",e18.11," | ",e18.11)') &
                   sum(matrix2int(0:gs2Tot(1,pr),1,2,pr)),            &
                   sum(matrix2int(0:gs2Tot(1,pr),2,3,pr)),            &
                   sum(matrix2int(0:gs2Tot(1,pr),3,1,pr))
          write(nx,*)
        endif
      endif

    endif

    write(nx,*)
    write(nx,'(1x,75("x"))')
    write(nx,*)
    write(nx,*)

  endif

  end subroutine print_squared_amplitude

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine print_squared_amplitude_cc (pr,i1,i2)

  integer, intent (in) :: pr,i1,i2

  integer       :: legs,n,j1,j2,gs
  character(2)  :: k,l
  character(99) :: a2,fmt

  legs = legsIn(pr) + legsOut(pr)
  n = 8 - 2*legs

  j1 = newleg(i1,pr)
  j2 = newleg(i2,pr)

  if (sum(abs(matrix2cc(0:gs2Tot(0,pr),j1,j2,pr))).ne.0d0) then

    call openOutput

    select case (n)
    case (-99:-10); fmt = '(1x,a,i3,a)'
    case (-9:-1);   fmt = '(1x,a,i2,a)'
    case (0:9);     fmt = '(1x,a,i1,a)'
    case default
    end select

    write(nx,'(1x,75("-"))')
    write(nx,*)
    write(nx,*)
    write(nx,trim(fmt)) &
      ' COLOUR-CORRELATED SQUARED AMPLITUDE [GeV^',n,']'
    write(nx,*)

    write(k,'(i2)') i1
    write(l,'(i2)') i2

    a2 = '| A0c('//trim(adjustl(k))//','//trim(adjustl(l))//') |^2'

    write(nx,*) ' als |    '//trim(a2)//'    '
    write(nx,*) ' ----------------------------'
    do gs = 0,gs2TotEff(0,pr),2
      write(nx,'(4x,i1," | ",e21.14)') gs/2,matrix2cc(gs,j1,j2,pr)
    enddo
    write(nx,*) ' ----------------------------'
    write(nx,'( "  SUM | ",e21.14)') sum(matrix2cc(0:gs2Tot(0,pr),j1,j2,pr))
    write(nx,*)
    write(nx,*)

  endif

  end subroutine print_squared_amplitude_cc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine print_squared_amplitude_cc_int (pr,i1,i2)

  integer, intent (in) :: pr,i1,i2

  integer       :: legs,n,j1,j2,gs
  character(2)  :: k,l
  character(99) :: a2,fmt

  legs = legsIn(pr) + legsOut(pr)
  n = 8 - 2*legs

  j1 = newleg(i1,pr)
  j2 = newleg(i2,pr)

  if (sum(abs(matrix2ccint(0:gs2Tot(1,pr),j1,j2,pr))).ne.0d0) then

    call openOutput

    select case (n)
    case (-99:-10); fmt = '(1x,a,i3,a)'
    case (-9:-1);   fmt = '(1x,a,i2,a)'
    case (0:9);     fmt = '(1x,a,i1,a)'
    case default
    end select

    write(nx,'(1x,75("-"))')
    write(nx,*)
    write(nx,*)
    write(nx,trim(fmt)) &
      ' COLOUR-CORRELATED SQUARED AMPLITUDE [GeV^',n,']'
    write(nx,*)

    write(k,'(i2)') i1
    write(l,'(i2)') i2

    a2 = '2*Re{ A0 T^i T^j A1* }('//trim(adjustl(k))//','//trim(adjustl(l))//')'

    write(nx,*) ' als |'//trim(a2)//'    '
    write(nx,*) ' --------------------------------'
    do gs = 0,gs2TotEff(0,pr),2
      write(nx,'(4x,i1," | ",e21.14)') gs/2,matrix2ccint(gs,j1,j2,pr)
    enddo
    do gs = gs2TotEff(0,pr)+2,gs2TotEff(1,pr),2
      write(nx,'(4x,i1," | ",e21.14)') gs/2,matrix2ccint(gs,j1,j2,pr)
    enddo
    write(nx,*) ' --------------------------------'
    write(nx,'( "  SUM | ",e21.14)') sum(matrix2ccint(0:gs2Tot(1,pr),j1,j2,pr))
    write(nx,*)
    write(nx,*)

  endif

  end subroutine print_squared_amplitude_cc_int

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine print_squared_amplitude_scc (pr,i1,i2,v)

  integer, intent (in)    :: pr,i1,i2
  complex(dp), intent(in) :: v(0:)

  integer       :: legs,n,j1,j2,gs
  character(2)  :: k,l
  character(99) :: a2,fmt

  legs = legsIn(pr) + legsOut(pr)
  n = 8 - 2*legs

  j1 = newleg(i1,pr)
  j2 = newleg(i2,pr)

  if (sum(abs(matrix2scc(0:gs2Tot(0,pr),j1,j2,pr))).ne.0d0) then

    call openOutput

    select case (n)
    case (-99:-10); fmt = '(1x,a,i3,a)'
    case (-9:-1);   fmt = '(1x,a,i2,a)'
    case (0:9);     fmt = '(1x,a,i1,a)'
    case default
    end select

    write(nx,'(1x,75("-"))')
    write(nx,*)
    write(nx,*)
    write(nx,trim(fmt)) &
      ' SPIN- AND COLOUR-CORRELATED SQUARED AMPLITUDE [GeV^',n,']'
    write(nx,*)

    write(k,'(i2)') i1
    write(l,'(i2)') i2

    write(nx,*) ' Polarization vector v for particle ',trim(adjustl(k)),':'
    write(nx,*) ' v(0) = ',v(0)
    write(nx,*) ' v(1) = ',v(1)
    write(nx,*) ' v(2) = ',v(2)
    write(nx,*) ' v(3) = ',v(3)
    write(nx,*)

    a2 = '| A0sc('//trim(adjustl(k))//','//trim(adjustl(l))//') |^2'

    write(nx,*) ' als |    '//trim(a2)//'    '
    write(nx,*) ' ----------------------------'
    do gs = 0,gs2TotEff(0,pr),2
      write(nx,'(4x,i1," | ",e21.14)') gs/2,matrix2scc(gs,j1,j2,pr)
    enddo
    write(nx,*) ' ----------------------------'
    write(nx,'( "  SUM | ",e21.14)') sum(matrix2scc(0:gs2Tot(0,pr),j1,j2,pr))
    write(nx,*)
    write(nx,*)
    write(nx,'(1x,75("x"))')
    write(nx,*)
    write(nx,*)

  endif

  end subroutine print_squared_amplitude_scc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine print_squared_amplitude_sc (pr,j,v)

  integer,     intent(in) :: pr,j
  complex(dp), intent(in) :: v(0:)

  integer       :: legs,n,jj,gs
  character(2)  :: cj
  character(99) :: a2,fmt

  legs = legsIn(pr) + legsOut(pr)
  n = 8 - 2*legs

  jj = newleg(j,pr)

  call openOutput

  select case (n)
  case (-99:-10); fmt = '(1x,a,i3,a)'
  case (-9:-1);   fmt = '(1x,a,i2,a)'
  case (0:9);     fmt = '(1x,a,i1,a)'
  case default
  end select

  write(nx,'(1x,75("-"))')
  write(nx,*)
  write(nx,*)
  write(nx,trim(fmt)) ' SPIN-CORRELATED SQUARED AMPLITUDE [GeV^',n,']'
  write(nx,*)

  write(cj,'(i2)') j

  write(nx,*) ' Polarization vector v for particle ',trim(adjustl(cj)),':'
  write(nx,*) ' v(0) = ',v(0)
  write(nx,*) ' v(1) = ',v(1)
  write(nx,*) ' v(2) = ',v(2)
  write(nx,*) ' v(3) = ',v(3)
  write(nx,*)

  a2 = '| A0s('//trim(adjustl(cj))//') |^2'

  write(nx,*) ' als |    '//trim(a2)//'    '
  write(nx,*) ' ----------------------------'
  do gs = 0,gs2TotEff(0,pr),2
    write(nx,'(4x,i1," | ",e21.14)') gs/2,matrix2sc(gs,pr)
  enddo
  write(nx,*) ' ----------------------------'
  write(nx,'( "  SUM | ",e21.14)') sum(matrix2sc(0:gs2Tot(0,pr),pr))
  write(nx,*)
  write(nx,*)
  write(nx,'(1x,75("x"))')
  write(nx,*)
  write(nx,*)

  end subroutine print_squared_amplitude_sc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module amplitude_rcl

!#####################################################################


