!#####################################################################
!!
!!  File  skeleton_rcl.f90
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

  module skeleton_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl
  use model_vertices_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

  integer                   :: bLDef,wLDef,bLInc,wLInc,bMaxT,wMaxT,&
                               qflowextT
  integer,      allocatable :: parT(:),csT(:,:),binT(:),hosT(:,:), &
                               hmaT(:),xxxT(:),gsT(:),wT(:,:),     &
                               typeT(:),gsIncT(:),qflowT(:),ffT(:)
  logical,      allocatable :: noquarks(:),nogluons(:),noweaks(:), &
                               U1gT(:)
  real (dp),    allocatable :: colT(:)
  complex (dp), allocatable :: couT(:,:)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine makeSkeleton ( lp,legs,pr )

  ! lp=0: tree
  ! lp=1: loop
  ! lp=2: tree for counterterms
  ! lp=3: tree for rational terms

  integer, intent (in) :: lp,legs,pr

  integer                  :: el1,el2,sumex,b,b0,legsE,bInc,ib,      &
                              par(legs),e1,e2,e3,e4,i,check,w,xTmin, &
                              xTmax,x,gs,j,npole,                    &
                              eb,pbins(legs),pe(legs),ppar(legs),    &
                              gphOut(legs),mass(legs),hm0,step,ii,n, &
                              nprop,n1,ee(2),sumA,sumB,typ,qed0,qed1
  integer,     allocatable :: cs4(:,:),w0e(:),gph(:,:),qed(:,:)
  logical                  :: qedcheck,noqed,doubleCounting,tadpole, &
                              selfenergy,selfenergyRes,last,U1g(3),  &
                              binOverlap
  logical,     allocatable :: xpOut(:),xp(:,:),xpIR(:),              &
                              epOut(:),ep(:,:)
  real(dp),    allocatable :: cl(:)
  complex(dp), allocatable :: c(:,:,:)
  character(2)             :: cc


  if (lp.eq.1) then
    el1   = 2**legs     ! binary number of the first  loop line
    el2   = 2**(legs+1) !   "      "    "   "  second  "    "
    sumex = el1-1       ! sum of binary numbers of external lines
  endif

  if (lp.eq.1) then
    legsE = legs + 2
  else
    legsE = legs
  endif

  par(1:legs) = parT(1:legs)

  npole = resMax(pr)
  do i = 1,npole
    pbins(i) = newbin(binRes(i,pr),pr)
    pe(i)    = min(pbins(i),2**legs-1-pbins(i))
    ppar(i)  = parRes(i,pr)
  enddo
  if (npole.gt.0) then
    allocate (xpOut(npole)); xpOut = .false.
    allocate (xp(npole,wLDef)); xp = .false.
    allocate (epOut(npole)); epOut = .false.
    allocate (ep(npole,wLDef)); ep = .false.
    allocate (xpIR(npole))
    allocate (gph(legs,wLDef)); gph    = 0
    gphOut = 0
  endif

  if (.not.(loopQED.and.loopWEAK).and.lp.eq.1) allocate(qed(0:1,wLDef))

  b = 0

  w = legsE

  allocate (w0e(2**(legsE-1)-1))

  e4loop: do e4 = 3, 2**(legsE-1)-1

    if (.not.resIR) then
      eb = e4
      if (eb.ge.2**legs) eb = e4 - 2**legs
      eb = min(eb,2**legs-1-eb)
      do i = 1,npole
        n = 0
        do j = 1,legs
          n = n + vectorLeg(pe(i),j)*vectorLeg(eb,j)
        enddo
        if ( n.gt.0 .and. n.lt.levelLeg(min(eb,pe(i))) ) cycle e4loop
      enddo
    endif

    w0e(e4) = w + 1

    ! In each branch the incoming binaries are from
    ! the biggest to the smallest.

    do e2 = 1, e4-1
      e1 = e4-e2
      if ( e1.gt.e2 ) then
        binOverlap = .false.
        do i = 1, legs + 2
          if (vectorLeg(e1,i)*vectorLeg(e2,i).eq.1) then
            binOverlap = .true.
            exit
          endif
        enddo
        if (.not.binOverlap) call router3 (e1,e2)
      endif
    enddo

    do e3 = 1, e4-1
      do e2 = e3+1, e4-e3-1
        e1 = e4-e3-e2
        if ( e1.gt.e2 ) then
          binOverlap = .false.
          do i = 1, legs + 2
            if ( + vectorLeg(e1,i)*vectorLeg(e2,i)      &
                 + vectorLeg(e1,i)*vectorLeg(e3,i)      &
                 + vectorLeg(e2,i)*vectorLeg(e3,i).gt.0 ) then
              binOverlap = .true.
              exit
            endif
          enddo
          if (.not.binOverlap) call router4 (e1,e2,e3)
        endif
      enddo
    enddo

    if (lp.ge.2) then
      e1 = e4
      call router2 (e1)
    endif

  enddo e4loop

  deallocate (w0e)

  bMaxT = b
  wMaxT = w

  if (.not.(loopQED.and.loopWEAK).and.lp.eq.1) deallocate(qed)

  if (npole.gt.0) deallocate (xp,xpOut,xpIR,ep,epOut,gph)

  contains

!---------------------------------------------------------------------

  subroutine router2 (e1)

  integer, intent (in) :: e1

  integer              :: e4,i1,i1max,l1,f1,f4,x1,x4,gs1,gs4,hm4, &
                          ff4,ho4(1:size(hosT,1)),                &
                          cs1(-1:size(csT,1)-2)
  integer, allocatable :: w1(:)
  logical              :: vexist(0:2),U1g4


  b0 = b

  l1 = legsE                 ! l1=i, if e1 is a binary of the i^th
  do i = 1,legsE-1           ! external particle, otherwise l1=legsE
    if (e1.eq.2**(i-1)) then; l1 = i; exit; endif
  enddo

  if (l1.eq.legsE.and.e1.ne.2**(legsE-1)-1) then ! just for internal
                                                 ! legs

    selfenergyRes = .false.
    if (.not.resSE.and.npole.gt.0) then
      ! 2-point contribution has to be excluded for resonant particles
      do i = 1,npole
        selfenergyRes = selfenergyRes .or.         &
                        (e1.eq.pbins(i)) .or.      &
                        (e1.eq.2**legs-1-pbins(i))
      enddo
    endif

    if (.not.selfenergyRes) then

      ! here we define incoming par, colour structure,
      ! bin and histories

      i1max = 0
      ! e1 is internal
      b1: do i = 1, b
        if (e1.eq.binT(wT(4,i))) then
          do i1= 1,i1max
            if (wT(4,i).eq.w1(i1)) cycle b1
          enddo
          i1max = i1max + 1; if (i1max.eq.1) allocate (w1(b))
          w1(i1max) =  wT(4,i)
        endif
      enddo b1

      ! here we select all possible outgoing par, colour structure,
      ! bin and histories
      ! the counter is increased each time by 1 or 2

      do i1 = 1,i1max

         f1 = parT(  w1(i1))
        cs1 =  csT(:,w1(i1))
         x1 = xxxT(  w1(i1))
        gs1 =  gsT(  w1(i1))

        if (x1.eq.0) then
          x = 1 ! for insertion of counterterms or R2 terms
        else
          cycle ! no tree level 2-point vertices
        endif

        e4 = e1

        x4 = 1

        f4loop: do f4 = 1,nFs ! f4 is outgoing

          if ((Qghost(f1)+Qghost(anti(f4))).ne.0) cycle
          if ((threeQ(f1)+threeQ(anti(f4))).ne.0) cycle

          if (npole.gt.0) then
             xpOut =  xp(:,w1(i1))
             epOut =  ep(:,w1(i1))
            gphOut = gph(:,w1(i1))
            do i = 1,npole
              xpOut(i) = xpOut(i) .or. (               &
                         nmf(f4).eq.nmf(ppar(i)) .and. &
                         ( e4.eq.pbins(i) .or.         &
                           e4.eq.2**legs-1-pbins(i) )  &
                         )
            enddo
          endif

          U1g4 = U1gT(w1(i1))
          ff4 = ffT(w1(i1))

          call vert2( legsE,lp,x,f1,anti(f4),cs1, & ! in
                      vexist(0:2),typ,cs4,c,cl    ) ! out

          do gs = 0,2 ! power of gs for this vertex

            if (.not.vexist(gs)) cycle

            gs4 = gs1 + gs

            ho4 = - 1
            hm4 = 0

            b = b + 1; if (b.gt.bLDef) call reallocate_bLDef

            wT(1,b) = w1(i1)
            wT(2,b) = 0
            wT(3,b) = 0

            check = 1
            wloop: do i = w0e(e4),w
              if (  f4.ne.parT(i)) cycle wloop
              if (  e4.ne.binT(i)) cycle wloop
              do ii = 1,size(hosT,1)
                if (ho4(ii).ne.hosT(ii,i)) cycle wloop
              enddo
              if ( hm4.ne.hmaT(i)) cycle wloop
              if (  x4.ne.xxxT(i)) cycle wloop
              if ( gs4.ne. gsT(i)) cycle wloop
              if (npole.gt.0) then
                do ii = 1,npole
                  if (xpOut(ii).neqv.xp(ii,i)) cycle wloop
                  if (epOut(ii).neqv.ep(ii,i)) cycle wloop
                enddo
                do ii = 1,legs
                  if (gphOut(ii).ne.gph(ii,i)) cycle wloop
                enddo
              endif
              if (U1g4.neqv.U1gT(i)) cycle wloop
              do ii = -1,size(csT,1)-2
                if (cs4(ii,1).ne.csT(ii,i)) cycle wloop
              enddo
              if (ff4.ne.ffT(i)) cycle wloop
              ! If we are here, current "i" is the same as the
              ! present one
              check = 0
              wT(4,b) = i
              exit wloop
            enddo wloop

            if(check.ne.0) then
              w = w + 1; if (w.gt.wLDef) call reallocate_wLDef
                wT(4,b) = w
              parT(  w) = f4
               csT(:,w) = cs4(:,1)
              binT(  w) = e4
              hosT(:,w) = ho4
              hmaT(  w) = hm4
              xxxT(  w) = x4
               gsT(  w) = gs4
              U1gT(  w) = U1g4
               ffT(  w) = ff4
              if (npole.gt.0) then
                 xp(:,w) =  xpOut
                 ep(:,w) =  epOut
                gph(:,w) = gphOut
              endif
            endif

            couT(1:4,b) = c(1:4,1,gs)
            colT(b) = cl(1)
            typeT(b) = typ
            gsIncT(b) = gs

          enddo

        enddo f4loop

      enddo

      if (allocated(w1)) deallocate (w1)

    endif

  endif

  end subroutine router2

!---------------------------------------------------------------------

  subroutine router3 (e1,e2)

  integer, intent (in) :: e1,e2

  integer              :: el,et,e4,i1,i2,i1max,i2max,l1,l2,f1,f2,f4, &
                          ho(1:size(hosT,1)),ho4(1:size(hosT,1)),    &
                          hm,hm4,x1,x2,x4,gs1,gs2,gs4,ff4,k,         &
                          cs1(-1:size(csT,1)-2),                     &
                          cs2(-1:size(csT,1)-2),cf1,cf2,cll
  integer, allocatable :: w1(:),w2(:)
  logical              :: vexist(0:3),U1g4

  b0 = b

  l1 = legsE                 ! l1=i, if e1 is a binary of the i^th
  do i = 1,legsE-1           ! external particle, otherwise l1=legsE
    if (e1.eq.2**(i-1)) then; l1 = i; exit; endif
  enddo
  l2 = legsE                 ! l2=i, if e2 is a binary of the i^th
  do i = 1,legsE-1           ! external particle, otherwise l2=legsE
    if (e2.eq.2**(i-1)) then; l2 = i; exit; endif
  enddo

  if (lp.eq.1) then
    if (e1 .ge. 2**legs ) then
      el = e1
      et = e2
    else
      el = 0
      et = 0
    endif
  endif

  tadpole = .false.
  if (lp.eq.1 .and. el.ne.0) tadpole = et .eq. sumex

  if (.not.tadpole) then

    selfenergy = .false.
    if (lp.eq.1 .and. el.ne.0) then
      do i= 1, legs
        if (selfenergy) then; exit
        else; selfenergy = et .eq. (sumex-2**(i-1))
        endif
      enddo
    endif

    if (.not.selfenergy) then

      ! here we define incoming par, colour structure,
      ! bin and histories

      i1max = 0
      select case (l1-legsE)
      case (:-1) ! e1 is external
        i1max = 1; allocate (w1(1));  w1(1) = l1
      case (0) ! e1 is internal
        b1: do i = 1, b
          if (e1.eq.binT(wT(4,i))) then
            do i1= 1,i1max
              if (wT(4,i).eq.w1(i1)) cycle b1
            enddo
            i1max = i1max + 1; if (i1max.eq.1) allocate (w1(b))
            w1(i1max) =  wT(4,i)
          endif
        enddo b1
      case default
        if (warnings(361).le.warning_limit) then
          warnings(361) = warnings(361) + 1
          call openOutput
          write(nx,*)
          write(nx,*) "CODE ERROR 361 (skeleton): router3"
          write(nx,*)
          call toomanywarnings(361)
        endif
        call istop (ifail,2)
      end select

      i2max = 0
      select case (l2-legsE)
      case (:-1) ! e2 is external
        i2max = 1; allocate (w2(1));  w2(1) = l2
      case (0) ! e2 is internal
        b2: do i = 1, b
          if (e2.eq.binT(wT(4,i))) then
            do i2= 1,i2max
              if (wT(4,i).eq.w2(i2)) cycle b2
            enddo
            i2max = i2max + 1; if (i2max.eq.1) allocate (w2(b))
            w2(i2max) =  wT(4,i)
          endif
        enddo b2
      case default
        if (warnings(362).le.warning_limit) then
          warnings(362) = warnings(362) + 1
          call openOutput
          write(nx,*)
          write(nx,*) "CODE ERROR 362 (skeleton): router3"
          write(nx,*)
          call toomanywarnings(362)
        endif
        call istop (ifail,2)
      end select

      qedcheck = (.not.(loopQED.and.loopWEAK)) .and. e1.ge.2**legs

      ! here we select all possible outgoing par, colour structure,
      ! bin and histories
      ! the counter is increased each time by 1 or 2

      e4 = e1 + e2
      last = .false.; if (e4.eq.2**(legsE-1)-1) last = .true.

      i1loop: do i1 = 1,i1max

        selfenergyRes = .false.
        if (.not.resSE.and.npole.gt.0.and.e1.ge.2**legs) then
          ! 1-loop 2-point contribution has to be excluded
          ! for resonant particles
          do i = 1,npole
            selfenergyRes = selfenergyRes .or.             &
                            ( ep(i,w1(i1)) .and.           &
                              ( et.eq.pbins(i) .or.        &
                                et.eq.2**legs-1-pbins(i) ) )
          enddo
        endif

        if (.not.selfenergyRes) then

           f1 = parT(  w1(i1))
          cs1 =  csT(:,w1(i1))
           x1 = xxxT(  w1(i1))
          gs1 =  gsT(  w1(i1))

          if (lp.eq.1) then
            ho = hosT(:,w1(i1))
            call checkDoubleCounting (et,el-el1,ho,doubleCounting)
          else
            doubleCounting = .false.
          endif
          if (doubleCounting) cycle i1loop ! double counting occurs

          if (lp.eq.1) then
            if (last) then
              ho4 = ho
            else
              ho4 = ho + vectorLeg(et,1:size(ho,1)) * &
                         ( maxval(ho) + 1 )
            endif
          else
            ho4 = - 1
          endif

          if (lp.eq.1) hm = hmaT(w1(i1))

          i2loop: do i2 = 1,i2max

             f2 = parT(  w2(i2))
            cs2 =  csT(:,w2(i2))
             x2 = xxxT(  w2(i2))
            gs2 =  gsT(  w2(i2))

            if (lp.ge.2) then
              if (e1+e2.eq.2**(legsE-1)-1.and.x1+x2.eq.0) then
                xTmin = 1
              else
                xTmin = 0
              endif
              if (x1+x2.eq.0) then
                xTmax = 1 ! for insertion of counterterms or R2 terms
              elseif (x1+x2.gt.1) then
                cycle ! just one insertion allowed
              else
                xTmax = 0
              endif
            else
              xTmin = 0
              if (x1+x2.gt.0) then
                cycle ! no insertion allowed
                      ! (to be safe, it should never be the case)
              else
                xTmax = 0
              endif
            endif

            do x = xTmin,xTmax

              x4 = max(x1+x2,x)

              f4loop: do f4 = 1,nFs ! f4 is incoming

                if ((Qghost(f1)+Qghost(f2)+Qghost(f4)).ne.0) cycle
                if ((threeQ(f1)+threeQ(f2)+threeQ(f4)).ne.0) cycle

                if (e1 .ge. 2**legs ) then
                  cc = cftype2(f4)
                  if ( noquarks(pr) .and.           &
                       .not.nogluons(pr) .and.      &
                       .not.noweaks(pr) .and.       &
                       cc.ne.'q' .and. cc.ne.'q~' ) &
                    cycle f4loop
                  if ( noquarks(pr) .and. nogluons(pr)        .and.   &
                       (cc.eq.'G'.or.cc.eq.'G~'.or.cc.eq.'g')       ) &
                    cycle f4loop
                  if (noquarks(pr) .and. noweaks(pr) .and.              &
                      cc.ne.'G' .and. cc.ne.'G~'.and. cc.ne.'g' .and.   &
                      cc.ne.'q' .and. cc.ne.'q~'                      ) &
                    cycle f4loop
                endif

                ! quark-flow constraint
                if (qflowextT .ne. 0) then

                  if (e4 .lt. 2**legs) then
                    cf1 = constr_quark(e1)
                    cf2 = constr_quark(e2)
                    if (cf1 .ne. 0) then
                      if (cf2 .ne. 0) then
                        ! two incoming quark lines -> constraint must match
                        if (qflowT(cf1) .ne. cf2) then
                          cycle f4loop
                        end if
                      ! quark line e2 with no constraint combined with
                      ! constrained quark line e1
                      else if (cftype2(f2) .eq. 'q' .or. cftype2(f2) .eq. 'q~') then
                        cycle f4loop

                      ! cf1 constrained to last leg ?
                      else if (e4 .eq. 2**(legs-1)-1 .and.  2**(qflowT(cf1)-1) .ne. e4+1) then
                        cycle f4loop
                      end if
                    else if (cf2 .ne. 0) then
                      ! quark line e1 with no constraint combined with
                      ! constrained quark line e2
                      if (cftype2(f1) .eq. 'q' .or. cftype2(f1) .eq. 'q~') then
                        cycle f4loop

                      ! cf2 constrained to last leg ?
                      else if (e4 .eq. 2**(legs-1)-1 .and.  2**(qflowT(cf2)-1) .ne. e4+1) then
                        cycle f4loop
                      end if
                    end if
                  else
                    ! Check if e2 has a fermion flow constraint (, if yes cf2>0)
                    cf2 = constr_quark(e2)
                    ! We only require the ff contraint for an incoming quark
                    ! line which is exiting the loop via f2
                    if ((cftype2(f1) .eq. 'q' .or. cftype2(f1) .eq. 'q~') .and. &
                        (cftype2(f2) .eq. 'q' .or. cftype2(f2) .eq. 'q~')) then
                      ! check if loop line has open fermion flow constraint (cll > 0)
                      cll = constr_quark_loopline(size(ho),ho)
                      if (cf2 .ne. 0) then
                        ! cll and cf2 have  constraint, but not resolved -> wrong fermion flow
                        if (cll .ne. qflowT(cf2) .and. cll .ne. 0) then
                          ! write(*,*) "Preventing:", e1, e2, cf2, cll, qflowT(4)
                          ! write(*,*) "ho:", ho
                          cycle f4loop
                        ! Final quark line unresolved ??
                        else if (cll .ne. qflowT(cf2) .and. e4 .eq. 2**(legs+1)-1) then
                          write(*,*) "Unhandled case:", e1, e2, cf2, cll
                          write(*,*) "ho:", ho
                          stop
                          cycle f4loop
                        end if
                      ! cll has constraint, but not resolved -> wrong fermion flow
                      else if (cll .ne. 0) then
                        ! write(*,*) "Preventing:", e1, e2, cf2, cll, qflowT(4)
                        ! write(*,*) "ho:", ho
                        cycle f4loop
                      end if
                    end if
                  end if
                end if

                ! condition on last particle
                if ( e4.eq.2**(legsE-1)-1 .and. f4.ne.parT(legsE) ) cycle

                if (npole.gt.0) then
                  xpOut = xp(:,w1(i1)).or.xp(:,w2(i2))
                  epOut = ep(:,w1(i1)).or.ep(:,w2(i2))
                  do i = 1,npole
                    xpOut(i) = xpOut(i) .or. (               &
                               nmf(f4).eq.nmf(ppar(i)) .and. &
                               ( e4.eq.pbins(i) .or.         &
                                 e4.eq.2**legs-1-pbins(i) )  &
                               )
                    epOut(i) = epOut(i) .or.                 &
                               ( e4.eq.pbins(i) .or.         &
                                 e4.eq.2**legs-1-pbins(i) )
                  enddo
                  if ( e1.eq.2**legs .and.    &
                       (f1.eq.15.or.f1.eq.16) ) then
                    gph(1,w1(i1)) = 1
                  endif
                  gphOut = gph(:,w1(i1))
                  if ( e1.ge.2**legs.and.(.not.last).and. &
                       (f4.eq.15.or.f4.eq.16)             ) then
                    gphOut(maxval(ho)+2) = 1
                  endif
                  xpIR = .false.
                  if (e4.eq.2**(legs+1)-1) then ! last loop leg
                    hm0 = hm
                    step = nmasses + 1
                    do ii = legs,1,-1
                      mass(ii) = hm0/step**(ii-1)
                      hm0 = hm0 - mass(ii)*step**(ii-1)
                    enddo
                    nprop = maxval(ho)+1
                    do i = 1,npole
                      proploop: do n = 1,nprop
                        xpIR(i) = gph(n,w1(i1)).eq.1
                        if (.not.xpIR(i)) cycle proploop
                        if (n.eq.1) then
                          ee(1) = e2
                          ee(2) = 0
                          do ii = 1,size(ho,1)
                            if (ho(ii)+1.eq.2) ee(2) = ee(2) + 2**(ii-1)
                          enddo
                        elseif (n.eq.nprop) then
                          ee(1) = 0
                          do ii = 1,size(ho,1)
                            if (ho(ii).gt.0.and.ho(ii)+1.eq.nprop) &
                              ee(1) = ee(1) + 2**(ii-1)
                          enddo
                          ee(2) = e2
                        else
                          ee(1) = 0
                          ee(2) = 0
                          do ii = 1,size(ho,1)
                            if (ho(ii).gt.0) then
                              if (ho(ii)+1.eq.n) &
                                ee(1) = ee(1) + 2**(ii-1)
                              if (ho(ii)+1.eq.n+1) &
                                ee(2) = ee(2) + 2**(ii-1)
                            endif
                          enddo
                        endif
                        iiloop: do ii = 1,2
                          do j = 1,npole
                            xpIR(i) = (levelLeg(ee(ii)).eq.1) .or. &
                                      (ee(ii).eq.pbins(j))
                            if (xpIR(i)) exit
                          enddo
                          if (.not.xpIR(i)) cycle iiloop
                          do n1 = 1,nprop
                            select case (ii)
                            case (1)
                              sumA = 0
                              sumB = e2
                              do k = 1,size(ho,1)
                                if (ho(k).gt.0) then
                                  if (ho(k)+1.ge.n+1 .and. &
                                      ho(k)+1.le.n1        ) &
                                    sumA = sumA + 2**(k-1)
                                  if ( (ho(k)+1.ge.n+1 .and.       &
                                        ho(k)+1.le.nprop    ) .or. &
                                       (ho(k)+1.ge.1 .and.         &
                                        ho(k)+1.le.n1     )      ) &
                                    sumB = sumB + 2**(k-1)
                                endif
                              enddo
                            case (2)
                              sumA = 0
                              sumB = e2
                              do k = 1,size(ho,1)
                                if (ho(k).gt.0) then
                                  if (ho(k)+1.ge.n1+1 .and. &
                                      ho(k)+1.le.n        ) &
                                    sumA = sumA + 2**(k-1)
                                  if ( (ho(k)+1.ge.n1+1 .and.       &
                                        ho(k)+1.le.nprop     ) .or. &
                                       (ho(k)+1.ge.1 .and.          &
                                        ho(k)+1.le.n      )       ) &
                                    sumB = sumB + 2**(k-1)
                                endif
                              enddo
                            end select
                            xpIR(i) = mass(n1).eq.nmf(ppar(i)) .and. &
                                      ( sumA.eq.pbins(i) .or.        &
                                        sumB.eq.pbins(i)      )
                            if (xpIR(i)) exit
                          enddo
                          if (xpIR(i)) then
                            exit proploop
                          else
                            cycle iiloop
                          endif
                        enddo iiloop
                      enddo proploop
                    enddo
                  endif
                  do i = 1,npole
                    xpIR(i) = xpIR(i) .and. resIR
                  enddo
                  do i = 1,npole
                    if (last.and.(.not.xpOut(i)).and.(.not.xpIR(i))) &
                    cycle f4loop
                  enddo
                endif

                U1g(1) = U1gT(w1(i1))
                U1g(2) = U1gT(w2(i2))
                if ( ffT(w1(i1)).ne.0 .and. ffT(w2(i2)).ne.0 ) then
                  ff4 = 0
                else
                  ff4 = ffT(w1(i1)) + ffT(w2(i2))
                endif

                call vert3( legsE,lp,x,last,csT(-1:0,legsE), & ! in
                            f1,f2,f4,cs1,cs2,U1g(1:2),       & ! in
                            vexist(0:3),typ,cs4,U1g4,c,cl    ) ! out

                if (qedcheck) then

                  if (e1.eq.2**legs) then
                    if ( f1.eq.16         .and. &
                         parKind(f2).ge.4 .and. &
                         parKind(f4).ge.4       ) then
                      qed0 = 1
                    else
                      qed0 = 0
                    endif
                    qed1 = 0
                  else
                    qed0 = qed(0,w1(i1))
                    qed1 = qed(1,w1(i1))
                  endif

                  if     ( qed1.eq.0        .and. &
                           parKind(f1).ge.4 .and. &
                           parKind(f2).ge.4 .and. &
                           f4.eq.16               ) then
                    qed1 = 1
                  elseif ( qed1.eq.1 ) then
                    if ( parKind(f2).ge.4 .and. &
                         parKind(f4).ge.4       ) then
                      qed1 = 2
                    else
                      qed1 = 0
                    endif
                  endif

                  if (last) then
                    if     ( qed0.eq.1        .and. &
                             parKind(f1).ge.4 .and. &
                             parKind(f2).ge.4 .and. &
                             f4.eq.16               ) then
                      qed0 = 2
                    else
                      qed0 = 0
                    endif
                  endif

                  if (last) then
                    noqed = qed0.ne.2 .and. qed1.ne.2
                    if( (noqed.eqv.loopQED) .and. &
                        (noqed.eqv..not.loopWEAK) ) &
                    vexist = .false.
                  endif

                endif

                do gs = 0,3 ! power of gs for this vertex

                  if (.not.vexist(gs)) cycle

                  gs4 = gs1 + gs2 + gs

                  bInc = size(cs4,2)

                  if (lp.eq.1) then
                    if (last) then; hm4 = hm
                    else;           hm4 = sethm(nmf(f4),hm)
                    endif
                  else;             hm4 = 0
                  endif

                  ibloop: do ib = 1,bInc

                    if (colour_optimization.ge.2) then
                      do i = 1,legsE
                        if (cs4(i,ib).eq.i) cycle ibloop
                      enddo
                    endif

                    b = b + 1; if (b.gt.bLDef) call reallocate_bLDef

                    wT(1,b) = w1(i1)
                    wT(2,b) = w2(i2)
                    wT(3,b) = 0

                    check = 1
                    wloop: do i = w0e(e4),w
                      if (anti(f4).ne.parT(i)) cycle wloop
                      if (     e4 .ne.binT(i)) cycle wloop
                      do ii = 1,size(hosT,1)
                        if (ho4(ii).ne.hosT(ii,i)) cycle wloop
                      enddo
                      if (    hm4 .ne.hmaT(i)) cycle wloop
                      if (     x4 .ne.xxxT(i)) cycle wloop
                      if (    gs4 .ne. gsT(i)) cycle wloop
                      if (npole.gt.0) then
                        do ii = 1,npole
                          if (xpOut(ii).neqv.xp(ii,i)) cycle wloop
                          if (epOut(ii).neqv.ep(ii,i)) cycle wloop
                        enddo
                        do ii = 1,legs
                          if (gphOut(ii).ne.gph(ii,i)) cycle wloop
                        enddo
                      endif
                      if (U1g4.neqv.U1gT(i)) cycle wloop
                      do ii = -1,size(csT,1)-2
                        if (cs4(ii,ib).ne.csT(ii,i)) cycle wloop
                      enddo
                      if (ff4.ne.ffT(i)) cycle wloop
                      if (qedcheck) then
                        if (qed0.ne.qed(0,i)) cycle wloop
                        if (qed1.ne.qed(1,i)) cycle wloop
                      endif
                      ! If we are here, current "i" is the same as the
                      ! present one
                      check = 0
                      wT(4,b) = i
                      exit wloop
                    enddo wloop

                    if (check.ne.0) then
                      w = w + 1; if (w.gt.wLDef) call reallocate_wLDef
                        wT(4,b) = w
                      parT(  w) = anti(f4)
                       csT(:,w) = cs4(:,ib)
                      binT(  w) = e4
                      hosT(:,w) = ho4
                      hmaT(  w) = hm4
                      xxxT(  w) = x4
                       gsT(  w) = gs4
                      U1gT(  w) = U1g4
                       ffT(  w) = ff4
                      if (npole.gt.0) then
                         xp(:,w) =  xpOut
                         ep(:,w) =  epOut
                        gph(:,w) = gphOut
                      endif
                      if (qedcheck) then
                        qed(0,w) = qed0
                        qed(1,w) = qed1
                      endif
                    endif

                    couT(1:2,b) = c(1:2,ib,gs)
                    couT(3:4,b) = 0d0
                    colT(b) = cl(ib)
                    typeT(b) = typ
                    gsIncT(b) = gs

                  enddo ibloop

                enddo

              enddo f4loop

            enddo

          enddo i2loop

        endif

      enddo i1loop

      if (allocated(w2)) deallocate (w2)
      if (allocated(w1)) deallocate (w1)

    endif   ! self-energy filter

  endif     ! tadpole filter

  end subroutine router3

!---------------------------------------------------------------------

  subroutine router4 (e1,e2,e3)

  integer, intent (in) :: e1,e2,e3

  integer              :: el,eta,etb,e4,i1,i2,i3,i1max,i2max,i3max,  &
                          l1,l2,l3,f1,f2,f3,f4,ho(1:size(hosT,1)),   &
                          ho4(1:size(hosT,1)),hm,hm4,x1,x2,x3,x4,    &
                          ff4,gs1,gs2,gs3,gs4,cs1(-1:size(csT,1)-2), &
                          cs2(-1:size(csT,1)-2),cs3(-1:size(csT,1)-2)
  integer, allocatable :: w1(:),w2(:),w3(:)
  logical              :: vexist(0:4),U1g4

  l1 = legsE                 ! l1=i, if e1 is a binary of the i^th
  do i = 1,legsE-1           ! external particle, otherwise l1=legsE
    if (e1.eq.2**(i-1)) then; l1 = i; exit; endif
  enddo
  l2 = legsE                 ! l2=i, if e2 is a binary of the i^th
  do i = 1,legsE-1           ! external particle, otherwise l2=legsE
    if (e2.eq.2**(i-1)) then; l2 = i; exit; endif
  enddo
  l3 = legsE                 ! l3=i, if e3 is a binary of the i^th
  do i = 1,legsE-1           ! external particle, otherwise l3=legsE
    if (e3.eq.2**(i-1)) then; l3 = i; exit; endif
  enddo

  if (lp.eq.1) then
    if (e1 .ge. 2**legs ) then
      el  = e1
      eta = e2
      etb = e3
    else
      el  = 0
      eta = 0
      etb = 0
    endif
  endif

  selfenergy = .false.
  selfenergyRes = .false.

  if (lp.eq.1 .and. el.ne.0) then
    do i= 1, legs
      if (selfenergy) then; exit
      else; selfenergy = ((eta+etb .eq. sumex) .and.         &
                         ((eta.eq.2**(i-1)).or.(etb.eq.2**(i-1))))
      endif
    enddo
      if (.not.resSE.and.npole.gt.0.and.e1.ge.2**legs) then
      ! 1-loop 2-point contribution has to be excluded
      ! for resonant particles
      do i = 1,npole
        selfenergyRes = selfenergyRes .or. &
                        ( (eta+etb .eq. sumex) .and.         &
                          ( (eta.eq.pbins(i)).or.(etb.eq.pbins(i)) ) )
      enddo
    endif
  endif

  if (.not.selfenergy.and..not.selfenergyRes) then

    ! here we define incoming par, colour structure,
    ! bin and histories

    i1max = 0
    select case (l1-legsE)
    case (:-1) ! e1 is external
        i1max = 1; allocate (w1(1));  w1(1) = l1
    case (0) ! e1 is internal
      b1: do i = 1, b
        if (e1.eq.binT(wT(4,i))) then
            do i1= 1,i1max
              if (wT(4,i).eq.w1(i1)) cycle b1
            enddo
            i1max = i1max + 1; if (i1max.eq.1) allocate (w1(b))
            w1(i1max) = wT(4,i)
        endif
      enddo b1
      case default
        if (warnings(363).le.warning_limit) then
          warnings(363) = warnings(363) + 1
          call openOutput
          write(nx,*)
          write(nx,*) "CODE ERROR 363 (skeleton): router4"
          write(nx,*)
          call toomanywarnings(363)
        endif
        call istop (ifail,2)
    end select

    i2max = 0
    select case (l2-legsE)
    case (:-1) ! e2 is external
        i2max = 1; allocate (w2(1));  w2(1) = l2
    case (0) ! e2 is internal
      b2: do i = 1, b
        if (e2.eq.binT(wT(4,i))) then
            do i2= 1,i2max
              if (wT(4,i).eq.w2(i2)) cycle b2
            enddo
            i2max = i2max + 1; if (i2max.eq.1) allocate (w2(b))
            w2(i2max) = wT(4,i)
        endif
      enddo b2
      case default
        if (warnings(364).le.warning_limit) then
          warnings(364) = warnings(364) + 1
          call openOutput
          write(nx,*)
          write(nx,*) "CODE ERROR 364 (skeleton): router4"
          write(nx,*)
          call toomanywarnings(364)
        endif
        call istop (ifail,2)
    end select

    i3max = 0
    select case (l3-legsE)
    case (:-1) ! e3 is external
        i3max = 1; allocate (w3(1));  w3(1) = l3
    case (0) ! e3 is internal
      b3: do i = 1, b
        if (e3.eq.binT(wT(4,i))) then
            do i3= 1,i3max
              if (wT(4,i).eq.w3(i3)) cycle b3
            enddo
            i3max = i3max + 1; if (i3max.eq.1) allocate (w3(b))
            w3(i3max) = wT(4,i)
        endif
      enddo b3
      case default
        if (warnings(365).le.warning_limit) then
          warnings(365) = warnings(365) + 1
          call openOutput
          write(nx,*)
          write(nx,*) "CODE ERROR 365 (skeleton): router4"
          write(nx,*)
          call toomanywarnings(365)
        endif
        call istop (ifail,2)
    end select

    qedcheck = (.not.(loopQED.and.loopWEAK)) .and. e1.ge.2**legs

    ! Here we select all possible outgoing par, colour structure,
    ! bin and histories.
    ! The counter is increased each time

    e4 = e1 + e2 + e3
    last = .false.; if (e4.eq.2**(legsE-1)-1) last = .true.

    i1loop: do i1 = 1,i1max

       f1 = parT(  w1(i1))
      select case (cftype(f1))
      case ('x','x~','f','f~'); cycle
      end select
      cs1 =  csT(:,w1(i1))
       x1 = xxxT(  w1(i1))
      gs1 =  gsT(  w1(i1))

      if (lp.eq.1) then
        ho = hosT(:,w1(i1))
        call checkDoubleCounting (eta+etb,el-el1,ho,doubleCounting)
      else
        doubleCounting = .false.
      endif
      if (doubleCounting) cycle i1loop ! double counting occurs

      if (lp.eq.1) hm = hmaT(w1(i1))

      i2loop: do i2 = 1,i2max

         f2 = parT(  w2(i2))
        select case (cftype(f2))
        case ('x','x~','f','f~'); cycle
        end select
        cs2 =  csT(:,w2(i2))
         x2 = xxxT(  w2(i2))
        gs2 =  gsT(  w2(i2))

       i3loop: do i3 = 1,i3max

           f3 = parT(  w3(i3))
          select case (cftype(f3))
          case ('x','x~','f','f~'); cycle
          end select
          cs3 =  csT(:,w3(i3))
           x3 = xxxT(  w3(i3))
          gs3 =  gsT(  w3(i3))

          if (lp.ge.2) then
            if (e1+e2+e3.eq.2**(legsE-1)-1.and.x1+x2+x3.eq.0) then
              xTmin = 1
            else
              xTmin = 0
            endif
            if (x1+x2+x3.eq.0) then
              xTmax = 1 ! for insertion of counterterms or R2 terms
            elseif (x1+x2+x3.gt.1) then
              cycle ! just one insertion allowed
            else
              xTmax = 0
            endif
          else
            xTmin = 0
            if (x1+x2+x3.gt.0) then
              cycle ! no insertion allowed
                    ! (to be safe, it should never be the case)
            else
              xTmax = 0
            endif
          endif

          if (lp.eq.1) then
            if (last) then
              ho4 = ho
            else
              ho4 = ho + vectorLeg(eta+etb,1:size(ho,1)) * &
                         ( maxval(ho) + 1 )
            endif
          else
            ho4 = - 1
          endif

          do x = xTmin,xTmax

            x4 = max(x1+x2+x3,x)

            f4loop: do f4 = 1,nFs ! f4 is incoming

              select case (cftype(f4))
              case ('x','x~','f','f~'); cycle
              end select
              if ((Qghost(f1)+Qghost(f2)+Qghost(f3)+Qghost(f4)).ne.0) cycle
              if ((threeQ(f1)+threeQ(f2)+threeQ(f3)+threeQ(f4)).ne.0) cycle

              if (e1 .ge. 2**legs ) then
                cc = cftype2(f4)
                if ( noquarks(pr) .and.           &
                     .not.nogluons(pr) .and.      &
                     .not.noweaks(pr) .and.       &
                     cc.ne.'q' .and. cc.ne.'q~' ) &
                  cycle f4loop
                if ( noquarks(pr) .and. nogluons(pr)        .and.   &
                     (cc.eq.'G'.or.cc.eq.'G~'.or.cc.eq.'g')       ) &
                  cycle f4loop
                if (noquarks(pr) .and. noweaks(pr) .and.              &
                    cc.ne.'G' .and. cc.ne.'G~'.and. cc.ne.'g' .and.   &
                    cc.ne.'q' .and. cc.ne.'q~'                      ) &
                  cycle f4loop
              endif

              ! condition on last particle
              if ( e4.eq.2**(legsE-1)-1 .and. f4.ne.parT(legsE) ) cycle

              if (npole.gt.0) then
                xpOut = xp(:,w1(i1)).or.xp(:,w2(i2)).or.xp(:,w3(i3))
                epOut = ep(:,w1(i1)).or.ep(:,w2(i2)).or.ep(:,w3(i3))
                do i = 1,npole
                  xpOut(i) = xpOut(i) .or. (               &
                             nmf(f4).eq.nmf(ppar(i)) .and. &
                             ( e4.eq.pbins(i) .or.         &
                               e4.eq.2**legs-1-pbins(i) )  &
                             )
                  epOut(i) = epOut(i) .or.                 &
                             ( e4.eq.pbins(i) .or.         &
                               e4.eq.2**legs-1-pbins(i) )
                enddo
                if (e1.eq.2**legs          .and. &
                    (f1.eq.15.or.f1.eq.16) ) then
                  gph(1,w1(i1)) = 1
                endif
                gphOut = gph(:,w1(i1))
                if (e1.ge.2**legs.and.(.not.last).and. &
                    (f4.eq.15.or.f4.eq.16)            ) then
                  gphOut(maxval(ho)+2) = 1
                endif
                do i = 1,npole
                  xpIR(i) = xpIR(i) .and. resIR
                enddo
                do i = 1,npole
                  if (last.and.(.not.xpOut(i)).and.(.not.xpIR(i))) &
                  cycle f4loop
                enddo
              endif

              U1g(1) = U1gT(w1(i1))
              U1g(2) = U1gT(w2(i2))
              U1g(3) = U1gT(w3(i3))

              ! no fermions in four-particle vertices
              ff4 = 0

              call vert4( legsE,lp,x,last,csT(-1:0,legsE),  & ! in
                          f1,f2,f3,f4,cs1,cs2,cs3,U1g(1:3), & ! in
                          vexist(0:4),typ,cs4,U1g4,c,cl     ) ! out

              if (qedcheck) then
                if (e1.eq.2**legs) then
                  qed0 = 0
                  qed1 = 0
                else
                  qed0 = qed(0,w1(i1))
                  qed1 = qed(1,w1(i1))
                endif
                if (last) then
                  noqed = qed0.ne.2 .and. qed1.ne.2
                  if( (noqed.eqv.loopQED) .and. &
                      (noqed.eqv..not.loopWEAK) ) &
                  vexist = .false.
                endif
              endif

              do gs = 0,4        ! power of gs

                if (.not.vexist(gs)) cycle

                gs4 = gs1 + gs2 + gs3 + gs

                bInc = size(cs4,2)

                if (lp.eq.1) then
                  if (last) then; hm4 = hm
                  else;           hm4 = sethm(nmf(f4),hm)
                  endif
                else;             hm4 = 0
                endif

                ibloop: do ib = 1,bInc

                  if (colour_optimization.ge.2) then
                    do i = 1,legsE
                      if (cs4(i,ib).eq.i) cycle ibloop
                    enddo
                  endif

                  b = b + 1; if (b.gt.bLDef) call reallocate_bLDef

                  wT(1,b) = w1(i1)
                  wT(2,b) = w2(i2)
                  wT(3,b) = w3(i3)

                  check = 1
                  wloop: do i = w0e(e4),w
                    if (anti(f4).ne.parT(i)) cycle wloop
                    if (     e4 .ne.binT(i)) cycle wloop
                    do ii = 1,size(hosT,1)
                      if (ho4(ii).ne.hosT(ii,i)) cycle wloop
                    enddo
                    if (    hm4 .ne.hmaT(i)) cycle wloop
                    if (     x4 .ne.xxxT(i)) cycle wloop
                    if (    gs4 .ne. gsT(i)) cycle wloop
                    if (npole.gt.0) then
                      do ii = 1,npole
                        if (xpOut(ii).neqv.xp(ii,i)) cycle wloop
                        if (epOut(ii).neqv.ep(ii,i)) cycle wloop
                      enddo
                      do ii = 1,legs
                        if (gphOut(ii).ne.gph(ii,i)) cycle wloop
                      enddo
                    endif
                    if (U1g4.neqv.U1gT(i)) cycle wloop
                    do ii = -1,size(csT,1)-2
                      if (cs4(ii,ib).ne.csT(ii,i)) cycle wloop
                    enddo
                    if (ff4.ne.ffT(i)) cycle wloop
                    if (qedcheck) then
                      if (qed0.ne.qed(0,i)) cycle wloop
                      if (qed1.ne.qed(1,i)) cycle wloop
                    endif
                    ! If we are here, current "i" is the same as the
                    ! present one
                    check = 0
                    wT(4,b) = i
                    exit wloop
                  enddo wloop

                  if (check.ne.0) then
                    w = w + 1; if (w.gt.wLDef) call reallocate_wLDef
                      wT(4,b) = w
                    parT(  w) = anti(f4)
                     csT(:,w) = cs4(:,ib)
                    binT(  w) = e4
                    hosT(:,w) = ho4
                    hmaT(  w) = hm4
                    xxxT(  w) = x4
                     gsT(  w) = gs4
                    U1gT(  w) = U1g4
                     ffT(  w) = ff4
                    if (npole.gt.0) then
                       xp(:,w) =  xpOut
                       ep(:,w) =  epOut
                      gph(:,w) = gphOut
                    endif
                    if (qedcheck) then
                      qed(0,w) = qed0
                      qed(1,w) = qed1
                    endif
                  endif

                  couT(1:3,b) = c(1:3,ib,gs)
                  couT(4,b) = 0d0
                  colT(b) = cl(ib)
                  typeT(b) = typ
                  gsIncT(b) = gs

                enddo ibloop

              enddo

            enddo f4loop

          enddo

        enddo i3loop

      enddo i2loop

    enddo i1loop

    if (allocated(w3)) deallocate (w3)
    if (allocated(w2)) deallocate (w2)
    if (allocated(w1)) deallocate (w1)

  endif   ! self-energy filter

  end subroutine router4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function constr_quark(bin)

  integer, intent(in) :: bin
  integer             :: constr_quark,i

  constr_quark = 0
  if (iand(qflowextT, bin) .ne. 0) then
    do i = 1, legs
      if (qflowT(i) .ne. 0) then
        if (iand(2**(qflowT(i)-1),bin) .ne. 0 .and. &
            iand(2**(i-1),bin) .eq. 0) then
            constr_quark = qflowT(i)
            return
        else if (iand(2**(qflowT(i)-1),bin) .eq. 0 .and. &
                 iand(2**(i-1),bin) .ne. 0) then
            constr_quark = i
            return
        end if
      end if
    end do
  end if

  end function constr_quark

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  recursive function constr_quark_loopline(sho,ho) result(cl)

  integer, intent(in)                 :: sho
  integer, dimension(sho), intent(in) :: ho
  integer, dimension(sho)             :: hor
  integer :: mv,i,tb,cf,cl

  mv = maxval(ho)
  if (mv .eq. 0) then
    cl = 0
    return
  end if
  tb = 0
  do i = 1, sho
    if (ho(i) .eq. mv) then
      tb = tb + 2**(i-1)
    end if
  end do

  cf = constr_quark(tb)
  if (cf .eq. 0) then
    hor = ho - vectorLeg(tb,:sho)*mv
    cl = constr_quark_loopline(sho,hor)
  else
    cl = cf
  end if

  end function constr_quark_loopline


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine reallocate_bLDef

  integer                   :: bdef0,bdef,ramMb
  integer,      allocatable :: tempI(:,:)
  real (dp),    allocatable :: tempR(:,:)
  complex (dp), allocatable :: tempC(:,:)

  bdef0 = bLDef
  bdef  = bLDef + bLInc

  if (writeRAM.ge.2) then
    ram0 = ram0 + bLInc*( 4*6 + 1*0 + 8*1 + 16*4 )
    ramMb = int(real(ram0,kind=dp)/1d6) + 2
    call openOutput
    write(nx,'(2x,a)') '->  free again'
    write(nx,'(1x,a,i3,a,i8,a)',advance='no') &
      'RAM temporally used by generation of process', &
      inpr(pr),':',ramMb,' Mbytes'
  endif

  allocate (tempI(4,bdef0))
    tempI = wT; deallocate (wT)
    allocate (wT(4,bdef)); wT(:,:bdef0) = tempI(:,:)
  deallocate (tempI)

  allocate (tempC(4,bdef0))
    tempC = couT; deallocate (couT);
    allocate (couT(4,bdef)); couT(:,:bdef0) = tempC(:,:)
  deallocate (tempC)

  allocate (tempR(1,bdef0))
    tempR(1,:) = colT(:); deallocate (colT)
    allocate (colT(bdef)); colT(:bdef0) = tempR(1,:)
  deallocate (tempR)

  allocate (tempI(1,bdef0))
    tempI(1,:) = typeT(:); deallocate (typeT)
    allocate (typeT(bdef)); typeT(:bdef0) = tempI(1,:)
    tempI(1,:) = gsIncT(:); deallocate (gsIncT)
    allocate (gsIncT(bdef)); gsIncT(:bdef0) = tempI(1,:)
  deallocate (tempI)

  bLDef = bdef

  end subroutine reallocate_bLDef

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine reallocate_wLDef

  integer              :: l,wdef0,wdef,ramMb
  integer, allocatable :: tempI(:,:)
  logical, allocatable :: tempL(:,:)

  wdef0 = wLDef
  wdef  = wLDef + wLInc

  if (writeRAM.ge.2) then
    ram0 = ram0 + wLInc*( 4*(legsE+legs+7) + 1*1 )
    if (npole.gt.0) &
      ram0 = ram0 + wLInc*( 4*legs + 1*2*npole )
    if (.not.(loopQED.and.loopWEAK).and.lp.eq.1) &
      ram0 = ram0 + wLInc*( 4*2 )
    ramMb = int(real(ram0,kind=dp)/1d6) + 2
    call openOutput
    write(nx,'(2x,a)') '->  free again'
    write(nx,'(1x,a,i3,a,i8,a)',advance='no') &
      'RAM temporally used by generation of process', &
      inpr(pr),':',ramMb,' Mbytes'
  endif

  l = size(csT,1)-2; allocate (tempI(-1:l,wdef0))
    tempI = csT; deallocate (csT)
    allocate (csT(-1:l,wdef)); csT(:,:wdef0) = tempI(:,:)
  deallocate (tempI)

  l = size(hosT,1); allocate (tempI(l,wdef0))
    tempI = hosT; deallocate (hosT)
    allocate (hosT(l,wdef)); hosT(:,:wdef0) = tempI(:,:)
  deallocate (tempI)

  allocate (tempI(1,wdef0))
    tempI(1,:) = parT(:); deallocate (parT)
    allocate (parT(wdef)); parT(:wdef0) = tempI(1,:)
    tempI(1,:) = binT(:); deallocate (binT)
    allocate (binT(wdef)); binT(:wdef0) = tempI(1,:)
    tempI(1,:) = hmaT(:); deallocate (hmaT)
    allocate (hmaT(wdef)); hmaT(:wdef0) = tempI(1,:)
    tempI(1,:) = xxxT(:); deallocate (xxxT)
    allocate (xxxT(wdef)); xxxT(:wdef0) = tempI(1,:)
    tempI(1,:) = gsT(:); deallocate (gsT)
    allocate (gsT(wdef)); gsT(:wdef0) = tempI(1,:)
    tempI(1,:) = ffT(:); deallocate (ffT)
    allocate (ffT(wdef)); ffT(:wdef0) = tempI(1,:)
  deallocate (tempI)

  allocate (tempL(1,wdef0))
    tempL(1,:) = U1gT(:); deallocate (U1gT)
    allocate (U1gT(wdef)); U1gT(:wdef0) = tempL(1,:)
  deallocate (tempL)

  if (npole.gt.0) then
    l = size(gph,1); allocate (tempI(l,wdef0))
      tempI = gph; deallocate (gph)
      allocate (gph(legs,wdef)); gph(:,:wdef0) = tempI(:,:)
    deallocate (tempI)
    l = size(xp,1); allocate (tempL(l,wdef0))
      tempL = xp; deallocate (xp)
      allocate (xp(l,wdef)); xp(:,:wdef0) = tempL(:,:)
    deallocate (tempL)
    l = size(ep,1); allocate (tempL(l,wdef0))
      tempL = ep; deallocate (ep)
      allocate (ep(l,wdef)); ep(:,:wdef0) = tempL(:,:)
    deallocate (tempL)
  endif

  if (.not.(loopQED.and.loopWEAK).and.lp.eq.1) then
    allocate (tempI(0:1,wdef0))
      tempI = qed; deallocate (qed)
      allocate (qed(0:1,wdef)); qed(:,:wdef0) = tempI(:,:)
    deallocate (tempI)
  endif

  wLDef = wdef

  end subroutine reallocate_wLDef

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end subroutine makeSkeleton

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine checkDoubleCounting (et,el,ho,doubleCounting)

  ! checks whether a double counting of one-loop diagrams occurs

  integer, intent (in)   :: et,el,ho(:)
  logical, intent (out)  :: doubleCounting

  integer                :: i,o1,g3,g4,g5
  integer, dimension (2) :: g,n

  ! Each offset is a binary number, i.e. is a sum without repetition
  ! of 1,2,4,8,...,2**(N-1) (N= number of external legs). They can
  ! therefore be thought as an ordered sequence (from the smallest to
  ! the biggest) of 1,2,4,8,...,2**(N-1) with or without gaps. The
  ! first number of the sequence (the smallest) is called the
  ! identifier of the offset.

  ! Two steps:

  ! 1) First I require that just an offset containing the 1 can
  ! enter at the beginning of the loop.

  ! 2) Then I require that three of the offsets must appear in a
  ! fixed order. I choose these three as the ones with the smallest
  ! identifiers and I require that these three offsets appear in the
  ! history of the offsets in the same order as their identifiers.
  ! Ex:
  ! offSet1 = 1+8    = 9   identifier1 = 1
  ! offSet2 = 2+4+32 = 38  identifier2 = 2
  ! offSet3 = 16     = 16  identifier3 = 16
  ! offSet4 = 64     = 64  identifier4 = 64
  ! I require that o1,o2,o3 always appear in this order
  ! (first offSet1, then offSet2, then offSet3). offSet4 can
  ! appear everywhere, also in between (not at the first
  ! position, which is excluded by the first step).

  if (ho(1).eq.-1) then
    doubleCounting = .false.
  else
    select case (mod(el+et,2))
    case (1)                        ! step 1
      o1 = 0
      do i = 1,size(ho,1)
        if (ho(i).eq.1) o1 = o1 + 2**(i-1)
      enddo
      ! el can be 0 or odd. When el=0, then et contains 1 and
      ! n(1)=g(1)=1. Independently on the value of n(2), we have
      ! g3=g4, and doubleCounting is .false. as it should be.
      ! g3=g4 because el is an empty set (there is only the loop
      ! propagator, with no offset). Let us then consider
      ! always the situation when el contains 1 and et does not.
      g(:)   = firstGaps(o1,1:2)    ! first two gaps (candidates)
      n(1:2) = firstNumbers(et,1:2) ! hits on the tree
      select case (n(1)-g(1))       ! step 2
      case (0)
        select case (n(2)-g(2))
        case (0)                    ! then we need a 3rd identifier
          g3 = firstGap(o1+et)
          g4 = firstGap(el+et)
          select case (g3-g4)
          case (0)
            doubleCounting = .false. ! this means that g3 has to come
          case default               ! still...accept
            doubleCounting = .true.  ! this means that g3 has already
          end select                 ! come; reject (g3>g(1))
        case default ! the first member of the tree has filled
                     ! the first gap;
          doubleCounting = .false.
        end select
      case default
        select case (n(1)-g(2))
        case (0)
          g5 = firstGap(el)
          select case (g5-g(1))
          case (0)
            doubleCounting = .true.  ! this means that g(1) has to
          case default               ! come still...reject
            doubleCounting = .false. ! this means that g(1) has come
          end select                 ! in; accept (g(1)<g5)
        case default ! when the first gap is not filled by the tree,
                     ! no double counting
          doubleCounting = .false.
        end select
      end select
    case default
      doubleCounting = .true.
    end select
  endif

  end subroutine checkDoubleCounting

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function sethm (m,hmIn) result (hmOut)

  integer, intent (in) :: m,hmIn ! m = mass in,
                                 ! hmIn = mass history in
  integer              :: hmOut  ! hmOut = mass history out

  integer              :: n,step,lmax

  step = nmasses + 1

  lmax = maxval(legsIn+legsOut)

  if (hmIn.eq.0) then
    hmOut = 0
  else
    do n = 1,lmax
      if (hmIn.lt.step**n) then
        hmOut = hmIn + m*step**n
        exit
      endif
    enddo
  endif

  end function sethm

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module skeleton_rcl

!#####################################################################

