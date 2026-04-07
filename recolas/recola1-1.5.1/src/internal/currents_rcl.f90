!#####################################################################
!!
!!  File  currents_rcl.f90
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

  module currents_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl
  use skeleton_rcl
  use model_vertices_rcl
  use draw_current_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine generate_currents

  integer                :: i,j,k,h,n,t,ii,kk,it,pr,b,w,w0,cs,c,s,   &
                            legs,w0Def,w0Inc,lmax,lmin,lmaxE,w0eDef, &
                            confDef,cdsDef,csDef,conf,cd0s,cd1s,     &
                            tlxm1,tlx,tlxp1,tln,sav,rep,ti,tiDef,ty, &
                            many(nFs),symfac,colfac,polfac,          &
                            legsE,base,tag,htag,htt,ferhelsum,       &
                            cuts(0:3),dl,cdsMax,tlm1,tl,tlp1,lp,     &
                            h1,ih1,njwTMax,njwTStep,firstlastloopj,  &
                            e1,e2,e3,e4,p1,p2,p3,p4,w1,w2,w3,w4,     &
                            fh,fh0,fh4,m4,check,check0,lexp,i0max,   &
                            sm0,sd0,sm1,sd1,bm0,bd0,bm1,bd1,         &
                            sm0Tot,sd0Tot,sm1Tot,sd1Tot,             &
                            bm0Tot,bd0Tot,bm1Tot,bd1Tot,             &
                            floop,csign1,csign2,csign3,ferloopsign,  &
                            fersign,ghostloopsign,signcorrection,    &
                            xIn,kbin,peak,lp0,lp0step,t0,i0,j0,      &
                            j1,j2,wTj1,wTj2,moda,daj1,moj2,daj2,     &
                            krran,rankIncrease,rankIn,ranktiMax,     &
                            ng,dCs0Tot,dCsTot,iqk,iak,step,hm0,      &
                            dMax,binMax,d,dbin,i1,i2,eq,l1,l2,nMax,  &
                            startw0,ew0,startbm0,startbd0,ebm0,ebd0, &
                            ramMb
  integer, allocatable   :: park(:),ffe(:),helk(:),hmn(:),hmx(:),    &
                            ph(:),hh(:),haux(:),he0(:,:),            &
                            w0min(:,:,:,:,:),w0max(:,:,:,:,:),       &
                            sec(:),nfer(:),ferp(:),fervec(:,:),      &
                            csign(:,:),fheT(:),                      &
                            parw0(:),binw0(:),csw0(:,:),xxxw0(:),    &
                            gsw0(:),ffw0(:),wMax(:),bMax(:),         &
                            sm0j(:),sd0j(:),sm1j(:),sd1j(:),         &
                            njwT(:),jwT(:,:),jwT0(:,:),              &
                            daj(:),moj(:),jdaTOjmo(:),ej(:),csj(:),  &
                            tij(:),jfirstwT(:),jinitwT(:),           &
                            massti(:),rankwT(:),rankwMax(:),         &
                            rankInj(:),rankOutj(:),peakleg(:),       &
                            wtTOw0b(:),wtTOw1b(:,:),w1stot(:),       &
                            w1ltot(:,:),c0Eff(:,:,:),cEff(:,:,:),    &
                            bm1c(:),bd1c(:),                         &
                            bm1h1(:,:,:),bd1h1(:,:,:),               &
                            startbm1(:),startbd1(:),ebm1(:),ebd1(:), &
                            vec(:),pCs(:,:),dCs0(:,:),dCs(:,:),      &
                            iq(:),ia(:,:),diq(:),dia(:,:),           &
                            ia1(:,:),ia2(:,:),                       &
                            iabin(:,:)
  integer, parameter :: nopol=-2
  logical                :: familybound,xs,newcut=.false.,newcs,     &
                            lda,ldai,lse,lsei

  logical, allocatable   :: checkHel(:,:),compConf(:),comp0(:,:,:),  &
                            comp1(:,:,:),okwti(:,:),lightmass(:),    &
                            extern(:),U1gw0(:),winitj(:),jfilter(:), &
                            fhfilter(:), wr(:,:,:,:)
  real (dp)              :: fac,fac1,ratio,corfac(4),facr,mONS2,mREG2
  real (dp), allocatable :: facj(:),dCs0fac(:),dCsfac(:)
  complex (dp)           :: co(4)
  character(2)           :: ck1,ck2,cc


  call vertices_tables

  bLDef = 500000; bLInc = 500000
  wLDef = 400000; wLInc = 400000
  w0Def = 300000; w0Inc = 300000

  lmax = maxval(legsIn+legsOut)
  lmin = minval(legsIn+legsOut)

  if (loopMax) then; lmaxE = lmax + 2
  else;              lmaxE = lmax
  endif

  w0eDef  = 0
  confDef = 0
  cdsDef  = 0
  csDef   = 0
  do pr = 1,prTot
    w0 = 0   ! maximal number of external currents
    conf = 1 ! maximal number of helicity configurations
    cd0s = 0 ! maximal number of colour deltas in each
             ! colour structure of tree lines
    do i = 1,legsIn(pr)+legsOut(pr)
      select case (par(i,pr))
      case  (1:14); w0 = w0 + 1
      case (17:19); w0 = w0 + 3; conf = conf * 3
      case default; w0 = w0 + 2; conf = conf * 2
      end select
    enddo
    do i = 1,legsIn(pr)+legsOut(pr)
      select case (par(i,pr))
      case  (15,23:25,29:31); cd0s = cd0s + 1
      end select
    enddo
    if (loop(pr)) then
      cd1s = cd0s + 2
    else
      cd1s = cd0s
    endif
    cs = 1 ! number of final colour structures (cs = cd0s!)
    do i= 1,cd0s
      cs = cs*i
    enddo
    w0eDef  = max(w0eDef ,w0)
    confDef = max(confDef,conf)
    cdsDef  = max(cdsDef ,cd1s)
    csDef   = max(csDef  ,cs)
  enddo

  ! Look whether there are external quarks, gluons and weak
  ! particles
  allocate (noquarks(prTot))
  allocate (nogluons(prTot))
  allocate (noweaks (prTot))
  do pr = 1,prTot
    noquarks(pr) = .true.
    nogluons(pr) = .true.
    noweaks(pr)  = .true.
    do i = 1,legsIn(pr)+legsOut(pr)
      cc = cftype2(par(i,pr))
      noquarks(pr) = noquarks(pr) .and. (cc.ne.'q').and.(cc.ne.'q~')
      nogluons(pr) = nogluons(pr) .and. (cc.ne.'g')
      noweaks(pr)  = noweaks(pr)  .and.                   &
                     (cc.ne.'s')  .and. (cc.ne.'v') .and. &
                     (cc.ne.'l')  .and. (cc.ne.'l~')
    enddo
  enddo

  tlxm1 = 2**(lmax-1)
  tlx   = 2**lmax
  tlxp1 = 2**(lmax+1)
  tln   = 2**lmin

  allocate (factor(prTot))

  allocate (newleg(lmax,prTot))
  allocate (oldleg(lmax,prTot))

  allocate (newbin   (tlx-1,prTot))
  allocate (oldbin   (tlx-1,prTot))
  allocate (defp2bin (tlx-1,prTot)); defp2bin = .false.
  allocate (p2bin    (tlx-1,prTot)); p2bin = 0d0
  allocate (defresbin(tlx-1,prTot)); defresbin = .false.
  allocate (pspbin   (tlx-1,prTot)); pspbin = 0d0

  allocate (  mONS(lmax,prTot))
  allocate (cmONS2(lmax,prTot))
  allocate (cmREG2(lmax,prTot))

  allocate ( gsTot(0:1,prTot))
  allocate (gs2Tot(0:1,prTot))
  allocate ( gsTotEff(0:1,prTot))
  allocate (gs2TotEff(0:1,prTot))
  allocate (cfTot(prTot))
  allocate (csTot(prTot))
  allocate (pCsTot(prTot))
  allocate ( csIa( lmax,csDef,prTot))
  allocate ( csIq( lmax,csDef,prTot))
  allocate (  nIa(      csDef,prTot))
  allocate (  pIa(csDef,csDef,prTot))
  allocate (facIa(csDef,csDef,prTot))

  allocate (w0eTot  (             prTot))
  allocate (heli    (lmax,confDef,prTot))
  allocate (dualheli(lmax,confDef,prTot))
  allocate(polprojin(tlx-1,prTot)); polprojin=nopol

  allocate (lpmax (prTot))

  allocate (cd0sMax(prTot))

  allocate(wr(43,0:43,0:43,43))

  if (loopMax) allocate (loopCoef(nFs,prTot))

  allocate (   w0Tot(prTot))
  allocate (   w1Tot(prTot))
  allocate (bm0prTot(prTot))
  allocate (bm1prTot(prTot))

  allocate (w0last(  -1:1,prTot))
  allocate (parw0e(w0eDef,prTot))
  allocate (binw0e(w0eDef,prTot))
  allocate (legw0e(w0eDef,prTot))
  allocate (helw0e(w0eDef,prTot))

  allocate (zeroLO(prTot))

  allocate (c0EffMax(prTot))
  if (loopMax) then
    allocate (cEffMax(prTot)); cEffMax = 0
  endif

  allocate (modaTot(prTot))

  if (loopMax) then
    allocate (  tiTot(prTot))
    allocate (ritiMax(prTot))
  endif

  sm0Tot = 0
  sd0Tot = 0
  sm1Tot = 0
  sd1Tot = 0

  bm0Tot = 0
  bd0Tot = 0
  bm1Tot = 0
  bd1Tot = 0

  tiDef = 0

  if (loopMax) then
    allocate (comp0gs(0:lmax,prTot))
    comp0gs = .false.
    allocate (comp1gs(0:lmax,prTot))
    comp1gs = .false.
  else
    allocate (comp0gs(0:lmax-2,prTot))
    comp0gs = .false.
  endif

  ! The code is practically run twice.
  ! For sav=0, the dimensions of some allocatable global variables
  ! are computed and for sav=1 the variables selves are allocated and
  ! computed
  saveloop: do sav = 0,1

    if (sav.eq.0) then
      allocate (wMax(prTot)); wMax = 0
      allocate (bMax(prTot)); bMax = 0
    endif

    if (sav.eq.1) then

      if (writeRAM.ge.1) then
        ram1 = 0
        ram1 = ram1 + sm0tot*( 4*11 + 1*1 + 8*0 + 16*4 )*prTot
        ram1 = ram1 + sd0tot*( 4*3  + 1*1 + 8*1 + 16*0 )*prTot
        if (loopMax) then
          ram1 = ram1 + sm1tot*( 4*11 + 1*1 + 8*0 + 16*4 )*prTot
          ram1 = ram1 + sd1tot*( 4*5  + 1*1 + 8*1 + 16*0 )*prTot
          ram1 = ram1 + ( tiDef*(lmax+2) + tiDef + 1 )*4*prTot
        endif
        ram1 = ram1 + maxval(c0EffMax)*( 1 + 4*(tlx-3)*confDef )*4*prTot
        ram1 = ram1 + bm0tot*( 4*6 + 1*1 + 8*0 + 16*0 )*prTot
        ram1 = ram1 + bd0tot*( 4*2 + 1*1 + 8*0 + 16*0 )*prTot
        if (loopMax) then
          ram1 = ram1 + maxval(cEffMax)*                        &
                        ( 3 + 2*(tlxp1-tln-1)*confDef )*4*prTot
          ram1 = ram1 + bm1tot*( 4*6 + 1*1 + 8*0 + 16*0 )*prTot
          ram1 = ram1 + bd1tot*( 4*2 + 1*1 + 8*0 + 16*0 )*prTot
          ram1 = ram1 + maxval(cEffMax)*2*4*prTot
        endif
        ram1 = ram1 + maxval(csTot)**2*( 4 + 8*lmax**2 )*prTot
        ramMb = int(real(ram1,kind=dp)/1d6) + 2
        call openOutput
        write(nx,'(1x,a,i8,a)') &
         'RAM permanently used by process generation in this run: ', &
         ramMb,' Mbytes'
        write(nx,*)
      endif

      allocate (   mosm0(    sm0Tot,prTot))
      allocate (  binsm0(1:3,sm0Tot,prTot))
      allocate (  parsm0(1:4,sm0Tot,prTot))
      allocate (    xsm0(    sm0Tot,prTot)) ! logical
      allocate (gsIncsm0(    sm0Tot,prTot))
      allocate (   cosm0(1:4,sm0Tot,prTot)) ! complex
      allocate (   gssm0(    sm0Tot,prTot))
      allocate (   cssm0(    sm0Tot,prTot))

      allocate (   dasd0(    sd0Tot,prTot))
      allocate (   sesd0(    sd0Tot,prTot)) ! logical
      allocate (  facsd0(    sd0Tot,prTot)) ! real
      allocate (   gssd0(    sd0Tot,prTot))
      allocate (   cssd0(    sd0Tot,prTot))

      if (loopMax) then
        allocate (     mosm1(    sm1Tot,prTot))
        allocate (    binsm1(1:3,sm1Tot,prTot))
        allocate (    parsm1(    sm1Tot,prTot))
        allocate (  gsIncsm1(    sm1Tot,prTot))
        allocate (     cosm1(1:4,sm1Tot,prTot)) ! complex
        allocate ( rankInsm1(    sm1Tot,prTot))
        allocate (rankOutsm1(    sm1Tot,prTot))
        allocate (ferloopsm1(    sm1Tot,prTot)) ! logical
        allocate (     gssm1(    sm1Tot,prTot))
        allocate (     cssm1(    sm1Tot,prTot))
        allocate (     tism1(    sm1Tot,prTot))
      endif

      if (loopMax) then
        allocate (     dasd1(    sd1Tot,prTot))
        allocate (    facsd1(    sd1Tot,prTot)) ! real
        allocate (rankOutsd1(    sd1Tot,prTot))
        allocate (ferloopsd1(    sd1Tot,prTot)) ! logical
        allocate (     gssd1(    sd1Tot,prTot))
        allocate (     cssd1(    sd1Tot,prTot))
        allocate (     tisd1(    sd1Tot,prTot))
      endif

      if (loopMax) then
        allocate (legsti(     0:tiDef,prTot)); legsti = 0
        allocate (momsti(lmax,0:tiDef,prTot))
        allocate (  vmti(lmax,0:tiDef,prTot))
        allocate (rankti(     0:tiDef,prTot)); rankti = 0
      endif

      allocate (c0TOlp(maxval(c0EffMax),prTot))
      if (loopMax) then
        allocate (bm0min(3:tlx-1,confDef,maxval(c0EffMax),prTot))
        allocate (bm0max(3:tlx-1,confDef,maxval(c0EffMax),prTot))
        allocate (bd0min(3:tlx-1,confDef,maxval(c0EffMax),prTot))
        allocate (bd0max(3:tlx-1,confDef,maxval(c0EffMax),prTot))
      else
        allocate (bm0min(3:tlxm1-1,confDef,maxval(c0EffMax),prTot))
        allocate (bm0max(3:tlxm1-1,confDef,maxval(c0EffMax),prTot))
        allocate (bd0min(3:tlxm1-1,confDef,maxval(c0EffMax),prTot))
        allocate (bd0max(3:tlxm1-1,confDef,maxval(c0EffMax),prTot))
      endif

      allocate (    sbm0(    bm0Tot,prTot))
      allocate ( w0inbm0(1:3,bm0Tot,prTot))
      allocate (w0outbm0(    bm0Tot,prTot))
      allocate (winitbm0(    bm0Tot,prTot)) ! logical
      allocate ( typebm0(    bm0Tot,prTot))

      allocate (    sbd0(    bd0Tot,prTot))
      allocate (w0outbd0(    bd0Tot,prTot))
      allocate (winitbd0(    bd0Tot,prTot)) ! logical

      if (loopMax) then
        allocate (cTOt  (maxval(cEffMax),prTot))
        allocate (cTOfh (maxval(cEffMax),prTot))
        allocate (cTOih1(maxval(cEffMax),prTot))
      endif

      if (loopMax) then
        !allocate (bm1min(tln+1:tlxp1-1,confDef,maxval(cEffMax),prTot))
        !allocate (bm1max(tln+1:tlxp1-1,confDef,maxval(cEffMax),prTot))
        allocate (    sbm1(    bm1Tot,prTot))
        allocate ( w1inbm1(    bm1Tot,prTot))
        allocate ( w0inbm1(2:3,bm1Tot,prTot))
        allocate (w1outbm1(    bm1Tot,prTot))
        allocate (winitbm1(    bm1Tot,prTot)) ! logical
        allocate ( typebm1(    bm1Tot,prTot))
      endif

      if (loopMax) then
        !allocate (bd1min(tln+1:tlxp1-1,confDef,maxval(cEffMax),prTot))
        !allocate (bd1max(tln+1:tlxp1-1,confDef,maxval(cEffMax),prTot))
        allocate (    sbd1(    bd1Tot,prTot))
        allocate (w1outbd1(    bd1Tot,prTot))
        allocate (winitbd1(    bd1Tot,prTot)) ! logical

        allocate(bm1_b(prTot), bd1_b(prTot))
        do pr = 1, prTot
          allocate(bm1_b(pr)%conf(2**(legsIn(pr)+legsOut(pr))+1:   &
                                  2**(legsIn(pr)+legsOut(pr)+1)-1, &
                                  cfTot(pr), cEffMax(pr)))
          allocate(bd1_b(pr)%conf(2**(legsIn(pr)+legsOut(pr))+1:   &
                                  2**(legsIn(pr)+legsOut(pr)+1)-1, &
                                  cfTot(pr), cEffMax(pr)))
        end do
      endif

      if (loopMax) then
        allocate (w1TotMax(maxval(cEffMax),prTot))
        w1TotMax = 0
        allocate (riwMax(maxval(cEffMax),prTot))
      endif

      allocate (colcoef (maxval(csTot),maxval(csTot),          prTot))
      allocate (colcoefc(maxval(csTot),maxval(csTot),lmax,lmax,prTot))

    endif

    ! Loop over the processes
    prloop: do pr = 1,prTot

      if (sav.eq.1) then
        if (loop(pr)) then
          if ( (c0EffMax(pr).eq.0).and.(cEffMax(pr).eq.0) ) then
            if (warnings(371).le.warning_limit) then
              warnings(371) = warnings(371) + 1
              call openOutput
              write(nx,*)
              if (stop_on_nonexisting_process) then
                write(nx,'(1x,a,i0,a)') &
                  'ERROR 371: process ',pr, &
                  ' does not exist neither at LO nor at NLO'
              else
                write(nx,'(1x,a,i0,a)') &
                  'WARNING 371: process ',pr, &
                  ' does not exist neither at LO nor at NLO'
              end if
              write(nx,*)
              call toomanywarnings(371)
            endif
            if (stop_on_nonexisting_process) call istop (ifail,1)
            prexists(pr) = .false.
            cycle prloop
          endif
        else
          if ( c0EffMax(pr).eq.0 ) then
            if (warnings(372).le.warning_limit) then
              warnings(372) = warnings(372) + 1
              call openOutput
              write(nx,*)
              write(nx,'(1x,a,i3,a)') &
                'ERROR 372: process',pr,' does not exist at LO'
              write(nx,*)
              call toomanywarnings(372)
            endif
            if (stop_on_nonexisting_process) call istop (ifail,1)
            prexists(pr) = .false.
            cycle prloop
          endif
        endif
      endif

      if (sav.eq.1) then
        call texHead(pr)
      endif

      legs = legsIn(pr) + legsOut(pr)

      if (sav.eq.0) then

        ! Compute symmetry factor for identical outgoing particles
        many = 0
        symfac = 1
        do i = legsIn(pr)+1,legs
          many(par(i,pr)) = many(par(i,pr)) + 1
          symfac = symfac * many(par(i,pr))
        enddo

        ! Compute colour averaging factor for incoming particles
        colfac = 1
        do i = 1,legsIn(pr)
          select case (par(i,pr))
          case (1,6,15)
            colfac = colfac * 8
          case (23:25,29:31,35:37,41:43)
            colfac = colfac * 3
          case default
            colfac = colfac * 1
          end select
        enddo

        ! Compute spin averaging factor for incoming particles
        polfac = 1
        do i = 1,legsIn(pr)
          if (hel(i,pr).ne.111) cycle
          select case (parKind(par(i,pr)))
          case (1)     ! scalar
            polfac = polfac * 1
          case (2)     ! massless vector
            polfac = polfac * 2
          case (3)     ! massive  vector
            polfac = polfac * 3
          case (4,5)   ! fermions and anti-fermions
            polfac = polfac * 2
          end select
        enddo

        ! Combine symmetry, spin and colour factors
        factor(pr) = 1d0 / polfac / colfac / symfac

        ! Reordering of the particles
        j = 0
        do i = 1,legs
          select case (par(i,pr))
          case (11:14,16:19,20:22,26:28,32:34,38:40)
            j = j + 1
            newleg(i,pr) = j
            oldleg(j,pr) = i
          end select
        enddo
        do i = 1,legs
          select case (par(i,pr))
          case (23:25,29:31,35:37,41:43)
            j = j + 1
            newleg(i,pr) = j
            oldleg(j,pr) = i
          end select
        enddo
        do i = 1,legs
          select case (par(i,pr))
          case (15)
            j = j + 1
            newleg(i,pr) = j
            oldleg(j,pr) = i
          end select
        enddo

!        Uncomment next lines to turn off the reordering
!        do i = 1,legs
!          newleg(i,pr) = i; oldleg(i,pr) = i
!        enddo

        do e1 = 1, 2**(legs)-1
          newbin(e1,pr) = 0
          do i = 1,legs
            newbin(e1,pr) = &
            newbin(e1,pr) + 2**(newleg(i,pr)-1)*vectorLeg(e1,i)
          enddo
          oldbin(newbin(e1,pr),pr) = e1
        enddo

        do i = 1,resMax(pr)
          if (resPar(parRes(i,pr))) then
            e1 = newbin(binRes(i,pr),pr)
            e2 = 2**legs - 1 - e1
            defresbin(e1,pr) = .true.
            defresbin(e2,pr) = .true.
            pspbin(e1,pr) = real(cm2f(parRes(i,pr)),kind=dp)
            pspbin(e2,pr) = real(cm2f(parRes(i,pr)),kind=dp)
            ! the sign was fixed based on previous implementation
            polprojin(e1,pr) = polproj(i,pr)
            if (abs(polproj(i,pr)) .eq. 1) then
              polprojin(e1,pr) = -polproj(i,pr)
              polprojin(e2,pr) = polproj(i,pr)
            else
              polprojin(e2,pr) = polproj(i,pr)
            end if
          else
            if (polproj(i, pr) .ne. nopol) then
              ! warn that process should define particle as resonant
                if (warnings(376).le.warning_limit) then
                    warnings(376) = warnings(376) + 1
                    call openOutput
                    write(nx,*)
                    write(nx,'(1x,a,i3,a,i3,a)') &
                    'WARNING 373: process',pr, &
                    ' does not define particle',parRes(i,pr), &
                    ' as resonant, but has internmediate pol assigned'
                    write(nx,*)
                    call toomanywarnings(376)
                end if
            end if
          endif
        enddo
      endif

      ! translate quark-flow variables to internal (permutated) ones
      allocate(qflowT(legs))
      qflowT = 0
      qflowextT = 0
      do i = 1, legs
        if (qflow(i,pr) .ne. 0) then
          if (cftype2(par(i,pr)) .ne. 'q' .and. &
              cftype2(par(qflow(i,pr),pr)) .ne. 'q') then
              write(*,*) "cftype2(par(i,pr)):", cftype2(par(i,pr))
              write(*,*) "cftype2(par(qflow(i,pr),pr)):", &
                         cftype2(par(qflow(i,pr),pr))
            write(*,*) "Invalid fermionflow combination. ", &
                       "One of the fields is non-fermionic"
            stop
          else if (cftype2(par(i,pr)) .eq. cftype2(par(qflow(i,pr),pr))) then
            write(*,*) "Invalid fermionflow combination. Fermion number violation."
            stop
          end if
          qflowT(newleg(i,pr)) = newleg(qflow(i,pr),pr)

          if (iand(qflowextT,2**(newleg(i,pr)-1)) .eq. 0) then
            qflowextT = qflowextT + 2**(newleg(i,pr)-1)
          end if

          if (iand(qflowextT,2**(qflowT(newleg(i,pr))-1)) .eq. 0) then
            qflowextT = qflowextT + 2**(qflowT(newleg(i,pr))-1)
          end if

        else
          qflowT(newleg(i,pr)) = 0
        end if
      end do

      ! Allocate some variables
      if (loop(pr)) then; allocate (park(legs+2))
      else;               allocate (park(legs))
      endif
      allocate ( ffe(legs))
      allocate (helk(legs))
      allocate ( hmn(legs))
      allocate ( hmx(legs))
      allocate (  ph(legs))

      ! External particles and helicities
      do i = 1,legs
        park(newleg(i,pr)) = par(i,pr)
        helk(newleg(i,pr)) = hel(i,pr)
      enddo

      if (loop(pr)) then; legsE = legs + 2
      else;               legsE = legs
      endif

      do i = 1,legs
        ! mONS2 is the squared of the on-shell mass of the external
        ! particle:
        ! - it is 0 if the particle is marked as light
        ! - it is the squared of the input mass if the particle is not
        !   marked as light
        ! mREG2 is the squared of the input mass of the external
        ! particle (if the particle is marked as light, mREG2 is used
        ! as IR-regulator in loop propagators only)
        mONS2 = real(cm2f(park(i)),kind=dp)
        mREG2 = real(cm2n(nmf(park(i))),kind=dp)
        mONS(i,pr) = sqrt(mONS2)
        cmONS2(i,pr) = mONS2*c1d0
        cmREG2(i,pr) = mREG2*c1d0
      end do

      ! compute fermion flow label for external legs
      ! - bosons have ffe = 0
      ! - fermions without fermion flow label have ffe = 1
      ! - fermions with fermion flow label have ffe = 2,3,4 ...
      ffe(1:legs) = 0
      j = 1
      do i = 1,legs
        if ( cftype(park(i)).eq.'f' .or. cftype(park(i)).eq.'f~') ffe(i) = 1
        if ( ffpar(park(i)) ) then
          j = j + 1
          ffe(i) = j
        endif
      enddo

      ! Select helicity configurations
      allocate (  hh(legs))
      allocate (haux(legs))
      allocate (he0(legs,confDef))
      do i= 1, legs
        select case (parkind(park(i)))
        case (1)     ! scalar (floating point exception for ph)
          hmn(i)=  0
          hmx(i)=  0
          ph (i)=  1
        case (2,4,5) ! massless vector, fermions and anti-fermions
          hmn(i)= -1
          hmx(i)= +1
          ph (i)=  2
        case (3)     ! massive  vector
          hmn(i)= -1
          hmx(i)= +1
          ph (i)=  1
        end select
      enddo
      base = 3
      allocate (checkHel(legs,-1:+1))
      checkHel = .false.
      do i= 1, legs
        do h = hmn(i), hmx(i), ph(i)
          checkHel(i,h) = .true.
        enddo
      enddo
      tag = 0
      hloo: do htag = 0, base**legs-1
        htt = htag
        leho: do i= legs,1,-1
          haux(i) = htt/base**(i-1)
          htt     = htt - haux(i)*base**(i-1)
          hh(i)    = haux(i) - 1
          select case (checkHel(i,hh(i)))
          case (.true.)
            cycle leho
          case (.false.)
            cycle hloo
          end select
        enddo leho
        tag = tag + 1
        he0(1:legs,tag)= hh(:)
      enddo hloo
      cfTot(pr) = tag
      deallocate (haux)
      deallocate (  hh)
      deallocate (checkHel)
      allocate (compConf(confDef))
      tag = 0
      do j= 1, cfTot(pr)
        ! select which configurations has to be computed
        ferhelsum = 0
        do ii = 1,legs
          select case (park(ii))
          case(20:25,32:37)
            if ( regf(park(ii))  .le.2 .and. &
                 regf(park(ii)+6).le.2 ) then
              ferhelsum = ferhelsum + he0(ii,j)
            endif
          case(26:31,38:43)
            if ( regf(park(ii))  .le.2 .and. &
                 regf(park(ii)-6).le.2 ) then
              ferhelsum = ferhelsum + he0(ii,j)
            endif
          end select
        enddo
        if (ferhelsum.eq.0) then
          ! Configurations where the helicity sum of fermions
          ! belonging to massless families is not 0 are discarded.
          compConf(j) = .true.
          do i = 1,legs
            if (abs(helk(i)).ne.111) then
              if (helk(i).eq.he0(i,j)) then
                hmn(i) = helk(i)
                hmx(i) = helk(i)
                ph(i)  = 1
              else
                compConf(j) = .false.
              endif
            endif
          enddo
          if (compConf(j)) then
            tag = tag + 1
            heli(1:legs,tag,pr)= he0(1:legs,j)
          endif
        endif
      enddo
      cfTot(pr) = tag
      deallocate (compConf)
      deallocate (he0)

      ! Define dualheli for spin-correlation
      if (sav.eq.1) then
        do ii = 1,legs
          if (park(ii).eq.15.or.park(ii).eq.16) then
            do j= 1, cfTot(pr)
              do j1= 1, cfTot(pr)
                if ( (sum(abs( heli(1:(ii-1),j1,pr)                  &
                              -heli(1:(ii-1),j ,pr))).eq.0)    .and. &
                     (sum(abs( heli((ii+1):legs,j1,pr)               &
                              -heli((ii+1):legs,j ,pr))).eq.0) .and. &
                     (heli(ii,j1,pr).eq.-heli(ii,j,pr))          ) then
                  dualheli(ii,j,pr) = j1
                endif
              enddo
            enddo
          endif
        enddo
      endif


      ! Set some global variables
      if (loop(pr)) then
        lpmax(pr) = 3
        cuts = 1
        cuts(1) = nFs
      else
        lpmax(pr) = 0
        cuts = 1
      endif

      ! deltas for colour structures.
      ! These deltas are made just of indices of the
      ! colours of the external lines
      cd0sMax(pr) = 0
      do dl = 1,legs
        ! The number of deltas is given by the number
        ! of gluons + the number of ghosts + the number
        ! of quark-antiquark pairs. This last number can
        ! be obtained counting the quarks (each of
        ! them is for sure in a pair)
        select case (park(dl))
        case (1,6,15,23:25,29:31)
          cd0sMax(pr) = cd0sMax(pr) + 1
        end select
      enddo
      if (loop(pr)) then
        cdsMax = cd0sMax(pr) + 2
      else
        cdsMax = cd0sMax(pr)
      endif

      ! Optimization of fermion loops (loops of fermions of different
      ! families, but same masses are computed just once).
      if (loop(pr)) then
        familybound = .false.
        do i = 1,legs
          select case (park(i))
          case(13,14,18,19); familybound = .true.; exit
          end select
          do j = 1,legs; if(j.eq.i) cycle
            if (park(i).eq.20.and.park(j).eq.38) familybound = .true.
            if (park(i).eq.21.and.park(j).eq.39) familybound = .true.
            if (park(i).eq.22.and.park(j).eq.40) familybound = .true.
            if (park(i).eq.23.and.park(j).eq.41) familybound = .true.
            if (park(i).eq.24.and.park(j).eq.42) familybound = .true.
            if (park(i).eq.25.and.park(j).eq.43) familybound = .true.
            if (park(i).eq.26.and.park(j).eq.32) familybound = .true.
            if (park(i).eq.27.and.park(j).eq.33) familybound = .true.
            if (park(i).eq.28.and.park(j).eq.34) familybound = .true.
            if (park(i).eq.29.and.park(j).eq.35) familybound = .true.
            if (park(i).eq.30.and.park(j).eq.36) familybound = .true.
            if (park(i).eq.31.and.park(j).eq.37) familybound = .true.
          enddo
        enddo
        allocate (lightmass(20:43))
        allocate (extern(20:43))
        extern = .false.
        do t = 20,25
          lightmass(t) = (regf(t).le.2).and. &
                         (.not.familybound.or.regf(t+6).le.2)
          do i= 1,legs
            extern(t) = extern(t)                             .or. &
                        (park(i).eq.t  ).or.(park(i).eq.t+12) .or. &
                        (park(i).eq.t+6).or.(park(i).eq.t+18)
          enddo
        enddo
        do t = 26,31
          lightmass(t) = (regf(t).le.2).and. &
                         (.not.familybound.or.regf(t-6).le.2)
          extern   (t) = extern   (t-6)
        enddo
        do t = 32,43
          lightmass(t) = lightmass(t-12)
          extern   (t) = extern   (t-12)
        enddo
        loopCoef(1:19,pr) = 1
        do i = 20,41,3
          loopCoef(i  ,pr) = 0
          loopCoef(i+1,pr) = 0
          loopCoef(i+2,pr) = 0
          if (lightmass(i).and.(.not.extern(i))) then
            if (extern(i+1).and.lightmass(i+1)) then
              loopCoef(i+1,pr) = loopCoef(i+1,pr) + 1
            elseif (extern(i+2).and.lightmass(i+2)) then
              loopCoef(i+2,pr) = loopCoef(i+2,pr) + 1
            else
              loopCoef(i,pr) = loopCoef(i,pr) + 1
            endif
          else
            loopCoef(i,pr) = loopCoef(i,pr) + 1
          endif
          if (lightmass(i+1).and.(.not.extern(i+1))) then
            if (lightmass(i).and. loopCoef(i,pr).ne.0) then
              loopCoef(i,pr) = loopCoef(i,pr) + 1
            elseif (extern(i+2).and.lightmass(i+2)) then
              loopCoef(i+2,pr) = loopCoef(i+2,pr) + 1
            else
              loopCoef(i+1,pr) = loopCoef(i+1,pr) + 1
            endif
          else
            loopCoef(i+1,pr) = loopCoef(i+1,pr) + 1
          endif
          if (lightmass(i+2).and.(.not.extern(i+2))) then
            if (lightmass(i).and.loopCoef(i,pr).ne.0) then
              loopCoef(i,pr) = loopCoef(i,pr) + 1
            elseif (lightmass(i+1).and.loopCoef(i+1,pr).ne.0) then
              loopCoef(i+1,pr) = loopCoef(i+1,pr) + 1
            else
              loopCoef(i+2,pr) = loopCoef(i+2,pr) + 1
            endif
          else
            loopCoef(i+2,pr) = loopCoef(i+2,pr) + 1
          endif
        enddo
        deallocate (extern)
        deallocate (lightmass)
      endif

      tlm1 = 2**(legs-1)
      tl   = 2**legs
      tlp1 = 2**(legs+1)

      if (sav.eq.1) then
        wLDef = wMax (pr)
        bLDef = bMax (pr)
        w0Def = w0Tot(pr)
      endif

      if (writeRAM.ge.2) then
        ram0 = 0
        ram0 = ram0 + 4*(tl-3)*cfTot(pr)*3*cuts(1)*(lpmax(pr)+1)*2
        ram0 = ram0 + w0Def*( 4*(legsE+6) +1*1 + 8*0 + 16*0 )
        ram0 = ram0 + wLDef*( 4*(legsE+legs+7) + 1*1 + 8*0 + 16*0 )
        ram0 = ram0 + bLDef*( 4*6 + 1*0 + 8*1 + 16*4 )
        if (resMax(pr).gt.0) then
          ram0 = ram0 + wLDef*( 4*legs + 1*2*resMax(pr) + 8*0 + 16*0 )
        endif
        if (.not.(loopQED.and.loopWEAK).and.loop(pr)) then
          ram0 = ram0 + wLDef*( 4*2 + 1*0 + 8*0 + 16*0 )
        endif
        ramMb = int(real(ram0,kind=dp)/1d6) + 2
        call openOutput
        write(nx,'(1x,a,i3,a,i8,a)',advance='no') &
          'RAM temporally used by generation of process',inpr(pr),':', &
          ramMb,' Mbytes'
      endif

      ! Allocate counters for currents
      if (loop(pr)) then
        allocate (w0min(3:tl-1,cfTot(pr),-1:1,cuts(1),0:lpmax(pr)))
        allocate (w0max(3:tl-1,cfTot(pr),-1:1,cuts(1),0:lpmax(pr)))
      else
        allocate (w0min(3:tlm1-1,cfTot(pr),-1:1,cuts(1),0:lpmax(pr)))
        allocate (w0max(3:tlm1-1,cfTot(pr),-1:1,cuts(1),0:lpmax(pr)))
      endif
      w0min = 0
      w0max = 0

      ! Initialize total number of global currents
      w0Tot(pr) = 0
      if (sav.eq.1) then
        w1Tot(pr) = 0
      endif

      ! The first currents are the external particles,
      ! according to their helicities
      allocate (parw0(         w0Def))
      allocate (binw0(         w0Def))
      allocate ( csw0(-1:legsE,w0Def))
      allocate (xxxw0(         w0Def))
      allocate ( gsw0(         w0Def))
      allocate (U1gw0(         w0Def))
      allocate ( ffw0(         w0Def))

      w0last(:,pr) = 0
      do i = 1,legs
        do h = hmn(i), hmx(i), ph(i)
          w0Tot(pr) = w0Tot(pr) + 1
          if (w0Tot(pr).gt.w0Def) call reallocate_w0Def
          parw0(  w0Tot(pr)) = park(i)
          binw0(  w0Tot(pr)) = 2**(i-1)
           csw0(:,w0Tot(pr)) = 0
          select case (park(i))
          case (1,6,15);      csw0(-1:0,w0Tot(pr)) = i
          case (23:25,29:31); csw0(-1  ,w0Tot(pr)) = i
          case (35:37,41:43); csw0(   0,w0Tot(pr)) = i
          end select
          xxxw0(w0Tot(pr)) = 0
           gsw0(w0Tot(pr)) = 0
          U1gw0(w0Tot(pr)) = (park(i).eq.15)
          if (colour_optimization.ge.2) U1gw0(w0Tot(pr)) = .false.
          ffw0(w0Tot(pr)) = ffe(i)
          ! for external particles
          parw0e(w0Tot(pr),pr) = park(i)
          binw0e(w0Tot(pr),pr) = 2**(i-1)
          legw0e(w0Tot(pr),pr) = i
          helw0e(w0Tot(pr),pr) = h
          if (i.eq.legs) w0last(h,pr) = w0Tot(pr)
        enddo
      enddo

      ! w0eTot(pr) counts the number of external particles
      ! with different helicities
      w0eTot(pr) = w0Tot(pr)

      deallocate (  ph)
      deallocate ( hmx)
      deallocate ( hmn)

      ! Initialize modaTot, which counts the mother and daugther
      ! currents in the colour optimization
      modaTot(pr) = 0

      ! Initialize number of tree and loop branches for
      ! unpolarized particles
      sm0 = 0
      sd0 = 0
      if (loop(pr)) then
        sm1 = 0
        sd1 = 0
      endif

      ! Initialize total number of tree and loop branches
      bm0 = 0
      bd0 = 0
      c0EffMax(pr) = 0
      allocate (c0Eff(-1:1,nFs,0:lpmax(pr)))
      c0Eff = 0
      allocate (comp0(-1:1,nFs,0:lpmax(pr)))
      comp0 = .false.

      if (loop(pr)) then
        bm1 = 0
        bd1 = 0
        if (sav.eq.1) then
          allocate (bm1h1(0:tlm1-1,-1:1,nFs))
          allocate (bd1h1(0:tlm1-1,-1:1,nFs))
          bm1h1 = 0
          bd1h1 = 0
          allocate (bm1c(maxval(cEffMax)))
          allocate (bd1c(maxval(cEffMax)))
          bm1c = 0
          bd1c = 0
        endif
        cEffMax(pr) = 0
        allocate (comp1(0:tlxm1-1,-1:1,nFs))
        comp1(:,:,:) = .false.
      endif

      ! Some inizializations
      if (sav.eq.0) then
        if (loop(pr)) then
          ti = 0
          tiTot(pr) = 0
        endif
      elseif (sav.eq.1) then
        bm0min(:,:,:,pr) = 0
        bm0max(:,:,:,pr) = 0
        bd0min(:,:,:,pr) = 0
        bd0max(:,:,:,pr) = 0
        if (loop(pr)) then
          ! bm1min(:,:,:,pr) = 0
          ! bm1max(:,:,:,pr) = 0
          ! bd1min(:,:,:,pr) = 0
          ! bd1max(:,:,:,pr) = 0
          bm1_b(pr)%conf(:,:,:)%bmin = 0
          bm1_b(pr)%conf(:,:,:)%bmax = 0
          bd1_b(pr)%conf(:,:,:)%bmin = 0
          bd1_b(pr)%conf(:,:,:)%bmax = 0
          allocate (massti(tiDef))
        endif
      endif

      pCsTot(pr) = 0
      allocate (pCs(lmax,csDef)); pCs(legs+1:lmax,:) = 0

      if (loop(pr)) then
        ! dZaa and dZe and its derived counterterms are process-dependent
        dZaa = dZaapr(pr)
        dZe  =  dZepr(pr)
        dgpn = dgpnpr(pr)
        dgpl = dgplpr(pr)
        dgpu = dgpupr(pr)
        dgpd = dgpdpr(pr)
        dgmn = dgmnpr(pr)
        dgml = dgmlpr(pr)
        dgmu = dgmupr(pr)
        dgmd = dgmdpr(pr)
      endif

      ! lp=0 -> Born
      ! lp=1 -> Loop
      ! lp=2 -> Counterterms
      ! lp=3 -> Rational terms
      lploop: do lp = 0,lpmax(pr)

        ! Effective number of legs (for loop we have two more legs)
        if (lp.eq.1) then
          legsE = legs + 2
        else
          legsE = legs
        endif

        lexp = 2**(legsE-1)

        allocate (ferp(legsE))
        allocate (fervec(2**(legsE-1)-1,legsE))
        allocate (csign(2**(legsE-1)-1,2**(legsE-1)-1))
        allocate (peakleg(legsE))

        if (lp.eq.1) then
          if (sav.eq.1) then
            allocate (w1sTot(0:tlm1-1))
            allocate (w1lTot(0:tlm1-1,cfTot(pr)))
          endif
          allocate (cEff(0:tlm1-1,-1:1,nFs))
          cEff = 0
        endif

        ! t is the particle of the cutted loop line
        cutsloop: do t = 1,cuts(lp)

          if (lp.eq.1) then
            if (loopCoef(t,pr).eq.0) cycle cutsloop
            cc = cftype2(t)
            if ( noquarks(pr) .and.           &
                 .not.nogluons(pr) .and.      &
                 .not.noweaks(pr) .and.       &
                 cc.ne.'q' .and. cc.ne.'q~' ) &
              cycle cutsloop
            if ( noquarks(pr) .and. nogluons(pr)        .and.   &
                 (cc.eq.'G'.or.cc.eq.'G~'.or.cc.eq.'g')       ) &
              cycle cutsloop
            if (noquarks(pr) .and. noweaks(pr) .and.              &
                cc.ne.'G' .and. cc.ne.'G~'.and. cc.ne.'g' .and.   &
                cc.ne.'q' .and. cc.ne.'q~'                      ) &
              cycle cutsloop
          end if

          ! Cutted loop lines
          if (lp.eq.1) then
            park(legs+1) = t
            park(legs+2) = anti(t)
          endif

          ! Get fermion relative signs according to hep-ph/0002082.
          ! The fermion relative signs obtained in this way are not
          ! always right at loop level and must be corrected by an
          ! overall sign (losigncorrection) which is computed later.
          fervec = 0
          do e1 = 1,2**(legsE-1)-1
            do i= 1, legsE
              if ( park(i).ge.20.and.   &
                   vectorLeg(e1,i).eq.1 ) fervec(e1,i) = 1
            enddo
          enddo
          do e2= 1, lexp-1 ! no last leg, at most sum of all others
          do e1= 1, lexp-1
            csign(e1,e2) = 1
            do i = 2,legsE
            do j = 1,i-1
              csign(e1,e2) = csign(e1,e2) *                    &
                             (-1)**(fervec(e1,i)*fervec(e2,j))
            enddo
            enddo
          enddo
          enddo

          ! Create skeleton. The last branches are filtered
          ! to get the correct particle.
          ! The "T" variables are local (they are initialized
          ! for each lp, case and config)

          allocate (  parT(         wLDef))
          allocate (   csT(-1:legsE,wlDef))
          allocate (  binT(         wLDef))
          allocate (  hosT(legs,    wLDef))
          allocate (  hmaT(         wLDef))
          allocate (  xxxT(         wLDef))
          allocate (   gsT(         wLDef))
          allocate (  U1gT(         wLDef))
          allocate (   ffT(         wLDef))
          allocate (    wT(4,bLDef))
          allocate ( typeT(  bLDef))
          allocate (gsIncT(  bLDef))
          allocate (  couT(4,bLDef))
          allocate (  colT(  bLDef))

          parT(1:legsE) = park(1:legsE)
          do it = 1,legsE
            binT(it) = 2**(it-1)
            csT(:,it) = 0
            select case (park(it))
            case (1,6,15);      csT(-1:0,it) = it
            case (23:25,29:31); csT(-1  ,it) = it
            case (35:37,41:43); csT(   0,it) = it
            end select
          enddo
          hosT(:,1:legs) = - 1
          hmaT(1:legs) = 0
          if (lp.eq.1) then
            hosT(:,legs+1) = 0
            hmaT(legs+1) = nmf(park(legs+1))
            hosT(:,legs+2) = - 1
            hmaT(legs+2) = 0
          endif
          xxxT(1:legsE) = 0
           gsT(1:legsE) = 0
          do it = 1,legsE
            U1gT(it) = (parT(it).eq.15)
          enddo
          if (colour_optimization.ge.2) U1gT(1:legsE) = .false.
          ffT(1:legs) = ffe(1:legs)
          if (legsE.gt.legs) then
            ffT(legs+1:legsE) = 0
            if ( cftype(t).eq.'f' .or. cftype(t).eq.'f~') ffT(legs+1:legsE) = 1
            if ( ffpar(t) ) then
              ffT(legs+1) = maxval(ffe) + 1
              ffT(legs+2) = maxval(ffe) + 2
            endif
          endif

          ! parT is outgoing for outgoing current
          call makeSkeleton ( lp,legs,pr ) ! in

          wMax(pr) = max(wMax(pr),wMaxT)
          bMax(pr) = max(bMax(pr),bMaxT)

          ! Filter from last branch to generate currents; note that
          ! the last branch contains the correct particle only
          allocate (jfilter(bMaxT))
          if (lp.eq.0) zeroLO(pr) = .true.
          jfilter1: do j = bMaxT,1,-1 ! reverse
            w = wT(4,j)
            jfilter(j) = .false.  ! filter inizialization
            if ( binT(w).eq.lexp-1 ) then
              jfilter(j) = .true. ! starting point
              if ( loop(pr) ) then
                if ( gsT(w)+2 .gt. size(powgs,1)-1 ) then
                  if ( lp.eq.0 .and.             &
                       powgs(gsT(w),0,pr).eq.0 ) &
                     jfilter(j) = .false.
                else
                  if ( lp.eq.0 ) then
                    if( powgs(gsT(w),0,pr).eq.0 .and. &
                        powgs(gsT(w)+2,1,pr).eq.0 )   &
                      jfilter(j) = .false.
                  endif
                endif
                if ( lp.ge.1 .and.             &
                     powgs(gsT(w),1,pr).eq.0 ) &
                   jfilter(j) = .false.
              else
                if ( lp.eq.0 .and.             &
                     powgs(gsT(w),0,pr).eq.0 ) &
                   jfilter(j) = .false.
              endif
            else
              check = 1
              do k = bMaxT, j+1, -1
                if (.not.jfilter(k)) cycle
                w1 = wT(1,k)
                w2 = wT(2,k)
                w3 = wT(3,k)
                if ((w-w1) .eq. 0 .or. (w-w2) .eq. 0 .or. (w-w3) .eq. 0) then
                  check = 0
                else
                  check = 1
                end if
                if (check.eq.0) then
                  jfilter(j) = .true.
                  cycle jfilter1
                endif
              enddo
              if (check.ne.0) jfilter(j) = .false.
            endif
            if (lp.eq.0) zeroLO(pr) = zeroLO(pr) .and. .not. jfilter(j)
          enddo jfilter1


          ! Colour optimization:
          ! Branches with dab(b,pr) not zero are "daugther" of
          ! another branch, i.e. their outgoing current can be
          ! obtained from that of another branch ("mother", with
          ! mob(b,pr) not zero), multiplying it by a colour
          ! factor
          allocate (njwT(0:wMaxT))
          njwTMax = (bMaxT/wMaxT+1) * legsE
          njwTStep = njwTMax
          allocate (jwT(0:njwTMax,0:wMaxT))
          allocate ( moj(bMaxT))
          allocate ( daj(bMaxT))
          allocate (facj(bMaxT))
          allocate (jdaTOjmo(bMaxT))
          njwT = 0
          jwT  = 0
          moj  = 0
          daj  = 0
          facj = 1d0
          jdaTOjmo = 0
          moda = 0
          if (colour_optimization.ge.1) then
            do j1 = 1, bMaxT
              if (.not.jfilter(j1)) cycle
              njwT(wT(4,j1)) = njwT(wT(4,j1)) + 1
              if (njwT(wT(4,j1)).gt.njwTMax) then
                allocate (jwT0(0:njwTMax,0:wMaxT))
                jwT0 = jwT
                deallocate (jwT)
                njwTMax = njwTMax + njwTStep
                allocate (jwT(0:njwTMax,0:wMaxT))
                jwT(0:njwTMax-njwTStep,:) = jwT0
                deallocate (jwT0)
              endif
              jwT(njwT(wT(4,j1)),wT(4,j1)) = j1
              j2loop: do j2 = 1,j1-1
                if (.not.jfilter(j2)) cycle
                if (typeT(j1).ne.typeT(j2)) cycle j2loop
                if (sum(abs(couT(:,j1)-couT(:,j2))).gt.zerocut) cycle j2loop
                do ii = 1,4
                  wTj1 = wT(ii,j1)
                  wTj2 = wT(ii,j2)
                  if (ii.eq.2.and.wTj1.eq.0.and.wTj2.eq.0) cycle
                  if (ii.eq.2.and.wTj1*wTj2.eq.0.and.wTj1.ne.wTj2) cycle j2loop
                  if (ii.eq.3.and.wTj1.eq.0.and.wTj2.eq.0) cycle
                  if (ii.eq.3.and.wTj1*wTj2.eq.0.and.wTj1.ne.wTj2) cycle j2loop
                  if (parT(wTj1).ne.parT(wTj2)) cycle j2loop
                  if (binT(wTj1).ne.binT(wTj2)) cycle j2loop
                  do k = 1,size(hosT,1)
                    if (hosT(k,wTj1).ne.hosT(k,wTj2)) cycle j2loop
                  enddo
                  if (hmaT(wTj1).ne.hmaT(wTj2)) cycle j2loop
                  if (xxxT(wTj1).ne.xxxT(wTj2)) cycle j2loop
                  if ( gsT(wTj1).ne. gsT(wTj2)) cycle j2loop
                  if (U1gT(wTj1).neqv.U1gT(wTj2)) cycle j2loop
                  if ( ffT(wTj1).ne. ffT(wTj2)) cycle j2loop
                enddo
                ratio = 1d0
                do ii= 1,3
                  wTj1 = wT(ii,j1)
                  wTj2 = wT(ii,j2)
                  if (ii.eq.2.and.wTj1.eq.0.and.wTj2.eq.0) cycle
                  if (ii.eq.3.and.wTj1.eq.0.and.wTj2.eq.0) cycle
                  if (wTj1.ne.wTj2) then
                    if (njwT(wTj1).ne.njwT(wTj2)) cycle j2loop
                    do n = 1,njwT(wTj1)
                      daj1 = daj(jwT(n,wTj1))
                      moj2 = moj(jwT(n,wTj2))
                      daj2 = daj(jwT(n,wTj2))
                      if (daj1.eq.0) cycle j2loop
                      if (daj1.ne.moj2.and.daj1.ne.daj2) cycle j2loop
                      fac  = facj(jwT(n,wTj1))
                      if (n.eq.1) fac1 = fac
                      if (abs(fac-fac1)/(abs(fac)+abs(fac1)).gt.zerocut) cycle j2loop
                    enddo
                    ratio = ratio * facj(jwT(1,wT(ii,j1))) / &
                                    facj(jwT(1,wT(ii,j2)))
                  endif
                enddo
                facj(j1) = facj(j1) * ratio * colT(j1)/colT(j2)
                if (moj(j2).eq.0) then
                  moda = moda + 1
                  moj(j2) = moda
                endif
                daj(j1) = moj(j2)
                jdaTOjmo(j1) = j2
                exit j2loop
              enddo j2loop
            enddo
          endif

          modaTot(pr) = max(modaTot(pr),moda)

          ! Filter again from last branch to avoid computing
          ! currents generating just daughter currents
          jfilter2: do j = bMaxT,1,-1 ! reverse
            w = wT(4,j)
            if ( (jfilter(j)).and.(binT(w).ne.lexp-1) ) then
              check = 1
              do k = bMaxT, j+1, -1
                if (.not.jfilter(k)) cycle
                w1 = wT(1,k)
                w2 = wT(2,k)
                w3 = wT(3,k)
                if (daj(k) .gt. 0) then
                  check = 1
                else if ((w-w1) .eq. 0 .or. (w-w2) .eq. 0 .or. &
                         (w-w3) .eq. 0) then
                  check = 0
                else
                  check = 1
                end if
                if (moj(j).ne.0) check = check*(moj(j)-daj(k))
                if (check.eq.0) then
                  jfilter(j) = .true.
                  cycle jfilter2
                endif
              enddo
              if (check.ne.0) jfilter(j) = .false.
            endif
          enddo jfilter2

          ! Some initializations
          allocate (ej(1:bMaxT))
          allocate (csj(1:bMaxT)); csj = 0
          if (sav.eq.1) then
            allocate (sm0j(1:bMaxT)); sm0j = 0
            allocate (sd0j(1:bMaxT)); sd0j = 0
            if (lp.eq.1) then
              allocate (    sm1j(1:bMaxT)); sm1j = 0
              allocate (    sd1j(1:bMaxT)); sd1j = 0
              allocate ( rankInj(1:bMaxT)); rankInj = 0
              allocate (rankOutj(1:bMaxT)); rankOutj = 0
              allocate (rankwT(    wMaxT)); rankwT = 0
              allocate (     tij(1:bMaxT)); tij = 0
            endif
          endif

          ! Initializa firstlastloopj
          firstlastloopj = bMaxT + 1

          ! Loop over the branches as created by skeleton
          ! (no polarization is introduced yet).
          ! Here will be fixed the "s" (from "skeleton") global
          ! variables which do not depend on the polarization

          jloop1: do j = 1,bMaxT

            newcut = newcut.or.(j.eq.1)

            if (.not.jfilter(j)) cycle jloop1

            if (moj(j).ne.0.and.daj(j).ne.0) then
              if (warnings(373).le.warning_limit) then
                warnings(373) = warnings(373) + 1
                call openOutput
                write(nx,*)
                write(nx,*) 'CODE ERROR 373 (currents_rcl.f90):'
                write(nx,*) 'branch with mo and da:',j,moj(j),daj(j)
                write(nx,*)
                call toomanywarnings(373)
              endif
              call istop (ifail,2)
            endif

            w1 = wT(1,j)
            w2 = wT(2,j)
            w3 = wT(3,j)
            w4 = wT(4,j)

            e1 = binT(w1)
            if (w2.eq.0) then
              e2 = 0
            else
              e2 = binT(w2)
            endif
            if (w3.eq.0) then
              e3 = 0
            else
              e3 = binT(w3)
            endif
            e4 = binT(w4)

            ej(j) = e4

            if ( firstlastloopj.eq.bMaxT+1 .and. &
                 e4.eq.2**(legs+1)-1             ) firstlastloopj = j

            ! Rough counter for tensor integrals tiDef (just done for
            ! sav=0 in order to allocate global variables for tensor
            ! integrals). The exact counter for tensor integrals tiTot
            ! is computed later for sav=1
            if (sav.eq.0.and.j.ge.firstlastloopj) then
              check = 1
              do j0= firstlastloopj, j-1
                if (jfilter(j0)) then
                  check = + sum(abs(hosT(:,w4)-hosT(:,wT(4,j0)))) &
                          + abs(hmaT(w4)-hmaT(wT(4,j0)))
                  if (check.eq.0) exit
                endif
              enddo
              if (check.ne.0) then
                ti = ti + 1
                tiDef = max(tiDef,ti)
              endif
            endif

            ! Counters s are increased
            if (e4.lt.2**legs) then ! tree branch
              if (daj(j).eq.0) then ! mother branch
                s = sm0 + 1
                sm0 = s
                sm0Tot = max(sm0Tot,sm0)
              else
                s = sd0 + 1
                sd0 = s
                sd0Tot = max(sd0Tot,sd0)
              endif
            else ! loop branch
              if (daj(j).eq.0) then
                s = sm1 + 1
                sm1 = s
                sm1Tot = max(sm1Tot,sm1)
              else
                s = sd1 + 1
                sd1 = s
                sd1Tot = max(sd1Tot,sd1)
              endif
            endif

            ! colour structures, last branch only
            if (e4.eq.2**(legsE-1)-1) then
              check = 1
              do k = 1,pCsTot(pr)
                check = sum(abs(csT(1:legs,w4)-pCs(1:legs,k)))
                if (check.eq.0) then
                  csj(j) = k
                  exit
                endif
              enddo
              if (check.ne.0) then
                pCsTot(pr) = pCsTot(pr) + 1
                csj(j) = pCsTot(pr)
                pCs(1:legs,pCsTot(pr)) = csT(1:legs,w4)
              endif
            endif

            if (sav.eq.1) then

              p1 = parT(w1)
              if (w2.eq.0) then
                p2 = 0
              else
                p2 = parT(w2)
              endif
              if (w3.eq.0) then
                p3 = 0
              else
                p3 = parT(w3)
              endif
              p4 = parT(w4)

              ! xs = T -> ct or r2 coupling
              ! xs = F -> normal tree coupling
              xIn = xxxT(w1)
              if (w2.ne.0) xIn = xIn + xxxT(w2)
              if (w3.ne.0) xIn = xIn + xxxT(w3)
              if (xIn.eq.0.and.xxxT(w4).eq.1) then
                xs = .true.
              else
                xs = .false.
              endif

              ! couplings
              co(1:4) = couT(1:4,j)*colT(j)

              ! Relative fermion sign
              fersign = 1
              if (w2.ne.0) then ! no two-leg branches
                if (w3.eq.0) then ! three-leg branch
                  select case (parT(w1))
                  case (20:31)
                    select case (parT(w2))
                    case (32:43)
                      fersign = csign(e2,e1)
                    case default
                      fersign = csign(e1,e2)
                    end select
                  case default
                    fersign = csign(e1,e2)
                  end select
                else                   ! four-leg branch
                  fersign = csign(e1,e2+e3)*csign(e2,e3)
                  ! check
                  csign1 = csign(e1,e2+e3)*csign(e2,e3)
                  csign2 = csign(e2,e3+e1)*csign(e3,e1)
                  csign3 = csign(e3,e1+e2)*csign(e1,e2)
                  if ( (csign1.ne.csign2) .or. &
                       (csign1.ne.csign3) .or. &
                       (csign2.ne.csign3)      ) then
                    if (warnings(374).le.warning_limit) then
                      warnings(374) = warnings(374) + 1
                      call openOutput
                      write(nx,*)
                      write(nx,*) 'CODE ERROR 374 (currents_rcl.f90): Fermion sign'
                      write(nx,*)
                      call toomanywarnings(374)
                    endif
                    call istop (ifail,2)
                  endif
                endif
              endif

              if (e4.ge.2**legs) then ! loop branch only

                ! rankInj,rankOutj,rankwT
                if (daj(j).ne.0) then
                  rankInj(j) = rankInj(jdaTOjmo(j))
                  rankOutj(j) = rankOutj(jdaTOjmo(j))
                  rankwT(w4) = max(rankwT(w4),rankwT(wT(4,jdaTOjmo(j))))
                else
                  rankIn = rankwT(w1) ! w1 is the loop-current
                  if (w3.eq.0) then ! three-leg branch
                    select case (parT(w4))
                    case (20:43); krran = 1
                    case default; krran = 0
                    end select
                    select case (typeT(j))
                    case (1:2,5:10,13:20,21:60,61,71)
                      rankIncrease = krran
                    case (3:4,11:12,62:63,72:73)
                      rankIncrease = krran + 1
                    case default
                      if (warnings(375).le.warning_limit) then
                        warnings(375) = warnings(375) + 1
                        call openOutput
                        write(nx,*)
                        write(nx,*) 'CODE ERROR 375 (currents_rcl.f90): no typeT',typeT(j)
                        write(nx,*)
                        call toomanywarnings(375)
                      endif
                      call istop (ifail,2)
                    end select
                  else                   ! four-leg branch
                    rankIncrease = 0
                  endif
                  rankwT(w4)  = max(rankwT(w4),rankIn+rankincrease)
                  rankInj(j)  = rankIn
                  rankOutj(j) = rankIn + rankincrease
                endif

                ferloopsign    = 1
                ghostloopsign  = 1
                signcorrection = 1

                if (e4.eq.2**(legsE-1)-1) then ! last branch only

                  ! Sign for fermion loops
                  allocate ( sec(legs+1))
                  allocate (nfer(legs+1))
                  if (hosT(1,w4).ne.-1) then
                    select case (t)
                    case (20:43)
                      floop = 0
                      do ii = 1, legs
                        nfer(ii) = 0
                        sec(ii) = 0
                        do k = 1,size(hosT,1)
                          if (hosT(k,w4)+1.eq.ii) sec(ii) = sec(ii) + 2**(k-1)
                        enddo
                      enddo
                      do ii = 1, legs ! loop over propagators
                        if (sec(ii) .ne. 0) then
                          do n = 1, legs
                            nfer(ii) = nfer(ii) + fervec(sec(ii),n)
                          enddo
                        endif
                      enddo
                      do ii= 1, legs
                        floop = floop + mod (nfer(ii),2)
                      enddo
                    case default
                      floop = 1
                    end select
                  else
                    floop = 1
                  endif
                  if (floop.eq.0) then
                    ferloopsign = - 1
                  else
                    ferloopsign = + 1
                  endif
                  deallocate (nfer)
                  deallocate ( sec)

                  ! Sign for ghost loops
                  select case (t)
                  case (1:10)
                    ghostloopsign = - 1
                  case default
                    ghostloopsign = + 1
                  end select

                  ! Correction of the relative fermion signs
                  select case (ferloopsign)
                  case (1)
                    select case (park(legs))
                    case (1:19)
                      select case (t)
                      case (1:19)
                        signcorrection = +1
                      case (20:43)
                        signcorrection = -1
                      end select
                    case (20:31)
                      select case (t)
                      case (1:19)
                        signcorrection = +1
                      case (20:43)
                        signcorrection = -1
                      end select
                    case (32:43)
                      select case (t)
                      case (1:19)
                        signcorrection = -1
                      case (20:43)
                        signcorrection = +1
                      end select
                    end select
                  case default
                    select case (park(legs))
                    case (1:19)
                      signcorrection = +1
                    case (20:31)
                      signcorrection = +1
                    case (32:43)
                      signcorrection = -1
                    end select
                  end select

                  ! tensor integrals
                  check = 1
                  do k = 1, titot(pr)
                    check = + sum(abs(hosT(:legs,wT(4,j))-momsti(:legs,k,pr))) &
                            + abs(hmaT(wT(4,j))-massti(k))
                    if (check.eq.0) then
                      tij(j) = k
                      rankti(k,pr) = max(rankti(k,pr),rankOutj(j))
                      exit
                    endif
                  enddo
                  if (check.ne.0) then
                    tiTot(pr) = tiTot(pr) + 1
                    tij(j) = tiTot(pr)
                    legsti(tiTot(pr),pr) = maxval(hosT(:,w4)) + 1
                    momsti(1:size(hosT,1),tiTot(pr),pr) = hosT(:,w4)
                    massti(tiTot(pr))    = hmaT(w4)
                    rankti(tiTot(pr),pr) = rankOutj(j)
                  endif

                endif

              endif

              ! Now we fix the value of the global variables
              ! for "s" branches

              if (e4.lt.2**legs) then ! tree branch

                if (daj(j).eq.0) then ! mother branch

                  sm0j(j) = s

                  ! Mother index
                  mosm0(s,pr) = moj(j)

                  ! binaries for s
                  binsm0(1,s,pr) = e1
                  binsm0(2,s,pr) = e2
                  binsm0(3,s,pr) = e3

                  ! Fields (all fields incoming)
                  parsm0(1,s,pr) = p1
                  parsm0(2,s,pr) = p2
                  parsm0(3,s,pr) = p3
                  parsm0(4,s,pr) = anti(p4)

                  ! xs = T -> ct or r2 coupling
                  ! xs = F -> normal tree coupling
                  xsm0(s,pr) = xs

                  ! Increase of gs power
                  gsIncsm0(s,pr) = gsIncT(j)

                  ! Global coefficient
                  cosm0(:,s,pr) = co(:) * fersign

                  ! Colour structures
                  ! (just last branches have no vanishing result)
                  cssm0(s,pr) = csj(j)

                  if (e4.eq.2**(legsE-1)-1) then ! last branch

                    ! gs power for the outgoing current
                    gssm0(s,pr) = gsT(w4)

                  endif

                else

                  sd0j(j) = s

                  ! Daughter index
                  dasd0(s,pr) = daj(j)

                  ! self-energy
                  sesd0(s,pr) = e4.eq.e1

                  ! Daughter factor
                  facsd0(s,pr) = facj(j)

                  ! Colour structures
                  ! (just last branches have no vanishing result)
                  cssd0(s,pr) = csj(j)

                  if (e4.eq.2**(legsE-1)-1) then ! last branch

                    ! gs power for the outgoing current
                    gssd0(s,pr) = gsT(w4)

                  endif

                endif

              else ! loop branch

                if (daj(j).eq.0) then

                  sm1j(j) = s

                  ! Mother index
                  mosm1(s,pr) = moj(j)

                  ! Binaries
                  binsm1(1,s,pr) = e1
                  binsm1(2,s,pr) = e2
                  binsm1(3,s,pr) = e3

                  ! Field of the outgoing current
                  parsm1(s,pr) = p4

                  ! Increase of gs power
                  gsIncsm1(s,pr) = gsIncT(j)

                  ! Global coefficient
                  cosm1(:,s,pr) = co(:) * fersign * signcorrection * &
                                          ferloopsign * ghostloopsign

                  ! rankIn and rankOut
                  rankInsm1(s,pr)  = rankInj(j)
                  rankOutsm1(s,pr) = rankOutj(j)

                  ! Colour structures (just last
                  ! branches have no vanishing result)
                  cssm1(s,pr) = csj(j)

                  if (e4.eq.2**(legsE-1)-1) then ! last branch

                    ! The last branch of a fermion loop has ferloop=T
                    ferloopsm1(s,pr) = (ferloopsign.eq.-1)

                    ! gs power for the outgoing current
                    gssm1(s,pr) = gsT(w4)

                    ! Tensor integrals
                    tism1(s,pr) = tij(j)

                  endif

                else

                  sd1j(j) = s

                  ! Daughter index
                  dasd1(s,pr) = daj(j)

                  ! Daughter factor
                  facsd1(s,pr) = facj(j)

                  ! rankOut
                  rankOutsd1(s,pr) = rankOutj(j)

                  ! Colour structures (just last
                  ! branches have no vanishing result)
                  cssd1(s,pr) = csj(j)

                  if (e4.eq.2**(legsE-1)-1) then ! last branch

                    ! The last branch of a fermion loop has ferloop=T
                    ferloopsd1(s,pr) = (ferloopsign.eq.-1)

                    ! gs power for the outgoing current
                    gssd1(s,pr) = gsT(w4)

                    ! Tensor integrals
                    tisd1(s,pr) = tij(j)

                  endif

                endif

              endif

              ! draw currents
              if (daj(j).eq.0) then
                call picture(pr, lp, legs, legsE,                    &
                             e1, e2, e3, e4, p1, p2, p3, p4, t, park,&
                             cdsMax, j, wT, csT, xs, newcut, wr)
                newcut = .false.
              endif

            endif

          enddo jloop1

          deallocate (csj)

          if (sav.eq.1) then

            if (lp.eq.1) then
              deallocate (tij)
              deallocate ( rankwT)
              deallocate (rankOutj)
            endif

            ! Computation of min and max for lmu, which is the index
            ! for the spinor/polarization vector for the cutted loop
            ! line. The first argument of minlmu and maxlmu is the
            ! "fermion helicity" fh of the cutted loop line (see
            ! later).
            select case (t)
            case (1:14)
              minlmu ( :,t) = 0
              maxlmu ( :,t) = 0
            case (15:19)
              minlmu ( :,t) = 0
              maxlmu ( :,t) = 3
            case default
              minlmu ( 0,t) = 0
              maxlmu ( 0,t) = 3
              minlmu (-1,t) = 2
              maxlmu (-1,t) = 3
              minlmu (+1,t) = 0
              maxlmu (+1,t) = 1
            end select

          endif

          allocate (fheT(0:wMaxT))
          allocate (fhfilter(bMaxT))

          ! Loop over fh of the cutted loop line. fh is the
          ! "fermion helicity", which is always 0 for bosons
          ! and for massive fermions and is -1 or +1 for
          ! left-handed or right-handed massless fermions
          ! respectively.

          fhmin(t) = 0
          fhmax(t) = 0
          if (lp.eq.1) then
            select case (cftype(t))
            case ('f','f~')
              if (regf(t).le.2) then
                fhmin(t) = - 1
                fhmax(t) = + 1
              endif
            end select
          endif

          fhloop: do fh = fhmax(t),fhmin(t),-2

            if (lp.eq.1) then
              allocate (okwti(wMaxT,cfTot(pr)))
              if (sav.eq.1) then
                allocate (wtTOw1b(wMaxT,cfTot(pr)))
                allocate (rankwMax(0:tlm1-1))
                rankwMax = 0
              endif
            endif

            ! The loop over helicity configuration (configloop)
            ! is done twice for sav=1.
            ! For rep=0 we count just the currents which are common
            ! to more than one helicity; this counter (w1sTot) is
            ! not initialized at each config, but goes through.
            ! For rep=1 the remaining currents are counted: this
            ! counter (w1lTot) is initialized at each config.
            ! The final counter of the currents is a proper
            ! combination of the two

            reploop: do rep = 0, sav

              if (lp.eq.1.and.rep.eq.0) then
                okwti = .false.
                if (sav.eq.1) then
                  wtTOw1b = 0
                  wtTOw1b (legs+1,:) = 1
                  wtTOw1b (legs+2,:) = 2
                endif
              endif

              ! The first loop currents has w1sTot = 1
              ! The  last loop currents has w1sTot = 2
              if (lp.eq.1) then
                if (sav.eq.1.and.rep.eq.0) then
                  w1sTot = 2
                endif
              endif

              ! Loop over helicity configurations

              configloop: do i = 1, cfTot(pr)

                ! helicities of external particles
                ! for this configuration
                helk(:) = heli(1:legs,i,pr)

                fheT = 111 ! all fheT are first set to 111
                fhfilter = .true.
                do k = 1,legs
                  ! all external particles get fheT = 0 except
                  ! massless fermions, which get fheT = +1,-1
                  fheT(k) = 0
                  select case (cftype(park(k)))
                  case ('f','f~')
                    if (regf(park(k)).le.2) then ! massless
                      fheT(k) = helk(k)
                    endif
                  end select
                enddo
                if (lp.eq.1) then
                  fheT(legs+1) = 0
                  fheT(legs+2) = 0
                  if (regf(t).le.2) then
                    fheT(legs+1) = + fh
                    fheT(legs+2) = + fh
                  endif
                endif

                ! Here the "fermion helicity" of the particles
                ! of the each branch are computed

                do j = 1,bMaxT

                  w1 = wT(1,j)
                  w2 = wT(2,j)
                  w3 = wT(3,j)
                  w4 = wT(4,j)
                  p1 = parT(w1)
                  if (w2.ne.0) p2 = parT(w2)
                  p4 = parT(w4)
                  ck1 = cftype(p1)
                  if (w2.ne.0) ck2 = cftype(p2)
!                  m1 = regf(p1)
!                  if (w2.ne.0) m2 = regf(p2)
                  m4 = regf(p4)

                  if (w2.eq.0) then
                    fh0 = abs(fheT(w1))
                  elseif (w3.eq.0) then
                    fh0 = abs(fheT(w1)) + abs(fheT(w2))
                  else
                    fh0 = + abs(fheT(w1)) + abs(fheT(w2)) &
                          + abs(fheT(w3))
                  endif

                  ! if fermion helicity is not physically possible, then fh4 = 3
                  ! if p4 = s or v or p4 = massive f or f~, it gets fh4 = 0
                  ! if p4 = massless f or f~, it gets fh4 = +1 or -1
                  if (fh0.ge.3) then
                    ! p1 or p2 or p3 comes from a wrong branch
                    fhfilter(j) = .false.
                    fh4 = 3
                  elseif (w2.eq.0) then ! 2-leg branch
                    fh4 = fheT(w1)
                  elseif (p1.eq.18.or.p1.eq.19) then ! p1 = W
                    if (ck2.eq.'f') then ! p2 = f
                      ! p4 = f for sure
                      if (fheT(w2).eq.+1) then ! p2 is right-handed -> NO
                        fhfilter(j) = .false.
                        fh4 = 3
                      else ! p2 is left-handed or not polarized
                        if (m4.le.2) then ! p4 massless f
                          if ( lp.ne.1.and. &
                               binT(w4).eq.2**(legs-1)-1 ) then
                            ! last tree current has no propagator
                            fh4 = + 1
                          else
                            fh4 = - 1
                          endif
                        else ! p4 massive f
                          fh4 = 0
                        endif
                      endif
                    elseif (ck2.eq.'f~') then ! p2 = f~
                      ! p4 = f~ for sure
                      if (fheT(w2).eq.-1) then ! p2 is left-handed
                        fhfilter(j) = .false.
                        fh4 = 3
                      else ! p2 is right-handed or not polarized
                        if (m4.le.2) then ! p4 massless f~
                          if ( lp.ne.1 .and. &
                               binT(w4).eq.2**(legs-1)-1 ) then
                            ! last tree current has no propagator
                            fh4 = - 1
                          else
                            fh4 = + 1
                          endif
                        else ! p4 massive f~
                          fh4 = 0
                        endif
                      endif
                    else ! p2 = scalar or vector
                      ! p4 = scalar or vector
                      fh4 = 0
                    endif
                  elseif (p1.ge.15.and.p1.le.17) then ! p1 = vector (not W)
                    if (m4.le.2) then ! p4 massless
                      if (lp.ne.1.and.binT(w4).eq.2**(legs-1)-1) then
                        ! last tree current has no propagator
                        fh4 = - fheT(w2)
                      else
                        fh4 = + fheT(w2)
                      endif
                    else ! p4 massive
                      fh4 = 0
                    endif
                  elseif (p1.eq.13.or.p1.eq.14) then ! p1 = phi+/-
                    ! Here we use that in the coupling f f'~ phi+/-
                    ! the omega_+ part is proportional to m_f and
                    ! the omega_- part is proportional to m_f'~.
                    if (ck2.eq.'f') then ! p2 = f
                      ! p4 = f for sure
                      if (fheT(w2).eq.+1) then ! p2 is massless and right-handed
                        fhfilter(j) = .false.
                        fh4 = 3
                      elseif (fheT(w2).eq.-1) then ! p2 is massless and left-handed
                        if (m4.le.2) then ! p4 massless
                          fhfilter(j) = .false.
                          fh4 = 3
                        else ! p4 massive
                          fh4 = 0
                        endif
                      else ! p2 is massive and not polarized
                        if (m4.le.2) then ! p4 massless
                          if (lp.ne.1.and.binT(w4).eq.2**(legs-1)-1) then
                            ! last tree current has no propagator
                            fh4 = + 1
                          else
                            fh4 = - 1
                          endif
                        else ! p4 massive
                          fh4 = 0
                        endif
                      endif
                    elseif (ck2.eq.'f~') then ! p2 = f~
                      ! p4 = f~ for sure
                      if (fheT(w2).eq.-1) then ! p2 is massless and left-handed
                        fhfilter(j) = .false.
                        fh4 = 3
                      elseif (fheT(w2).eq.+1) then ! p2 is massless and right-handed
                        if (m4.le.2) then ! p4 massless
                          fhfilter(j) = .false.
                          fh4 = 3
                        else ! p4 massive
                          fh4 = 0
                        endif
                      else ! p2 is massive and not polarized
                        if (m4.le.2) then ! p4 massless
                          if (lp.ne.1.and.binT(w4).eq.2**(legs-1)-1) then
                            ! last tree current has no propagator
                            fh4 = - 1
                          else
                            fh4 = + 1
                          endif
                        else ! p4 massive
                          fh4 = 0
                        endif
                      endif
                    else ! p2 = s or v
                      fh4 = 0
                    endif
                  elseif (p1.le.12) then ! p1 = scalar (not phi+/-)
                    if (m4.le.2) then ! p4 massless
                      if (lp.ne.1.and.binT(w4).eq.2**(legs-1)-1) then
                        ! last tree current has no propagator
                        fh4 = + fheT(w2)
                      else
                        fh4 = - fheT(w2)
                      endif
                    else ! p4 massive
                      fh4 = 0
                    endif
                  elseif (p2.eq.18.or.p2.eq.19) then ! p2 = W
                    if (ck1.eq.'f') then ! p1 = f
                      ! p4 = f for sure
                      if (fheT(w1).eq.+1) then ! p1 is right-handed
                        fhfilter(j) = .false.
                        fh4 = 3
                      else ! p1 is left-handed or not polarized
                        if (m4.le.2) then ! p4 massless f
                          if ( lp.ne.1 .and. &
                               binT(w4).eq.2**(legs-1)-1 ) then
                            ! last tree current has no propagator
                            fh4 = + 1
                          else
                            fh4 = - 1
                          endif
                        else ! p4 massive f
                          fh4 = 0
                        endif
                      endif
                    elseif (ck1.eq.'f~') then ! p1 = f~
                      ! p4 = f~ for sure
                      if (fheT(w1).eq.-1) then ! p1 is left-handed
                        fhfilter(j) = .false.
                        fh4 = 3
                      else ! p1 is right-handed or not polarized
                        if (m4.le.2) then ! p4 massless f~
                          if ( lp.ne.1 .and. &
                               binT(w4).eq.2**(legs-1)-1 ) then
                            ! last tree current has no propagator
                            fh4 = - 1
                          else
                            fh4 = + 1
                          endif
                        else ! p4 massive f~
                          fh4 = 0
                        endif
                      endif
                    else ! p1 = scalar or vector
                      ! p4 = scalar or vector
                      fh4 = 0
                    endif
                  elseif (p2.ge.15.and.p2.le.17) then ! p2 = vector (not W)
                    if (m4.le.2) then ! p4 massless
                      if (lp.ne.1.and.binT(w4).eq.2**(legs-1)-1) then
                        ! last tree current has no propagator
                        fh4 = - fheT(w1)
                      else
                      fh4 = + fheT(w1)
                      endif
                    else ! p4 massive
                      fh4 = 0
                    endif
                  elseif (p2.eq.13.or.p2.eq.14) then ! p2 = phi+/-
                    ! Here we use that in the coupling f f'~ phi+/-
                    ! the omega_+ part is proportional to m_f and
                    ! the omega_- part is proportional to m_f'~.
                    if (ck1.eq.'f') then ! p1 = f
                      ! p4 = f for sure
                      if (fheT(w1).eq.+1) then ! p1 is massless and right-handed
                        fhfilter(j) = .false.
                        fh4 = 3
                      elseif (fheT(w1).eq.-1) then ! p1 is massless and left-handed
                        if (m4.le.2) then ! p4 massless
                          fhfilter(j) = .false.
                          fh4 = 3
                        else ! p4 massive
                          fh4 = 0
                        endif
                      else ! p1 is massive and not polarized
                        if (m4.le.2) then ! p4 massless
                          if (lp.ne.1.and.binT(w4).eq.2**(legs-1)-1) then
                            ! last tree current has no propagator
                            fh4 = + 1
                          else
                            fh4 = - 1
                          endif
                        else ! p4 massive
                          fh4 = 0
                        endif
                      endif
                    elseif (ck1.eq.'f~') then ! p1 = f~
                      ! p4 = f~ for sure
                      if (fheT(w1).eq.-1) then ! p1 is massless and left-handed
                        fhfilter(j) = .false.
                        fh4 = 3
                      elseif (fheT(w1).eq.+1) then ! p1 is massless and right-handed
                        if (m4.le.2) then ! p4 massless
                          fhfilter(j) = .false.
                          fh4 = 3
                        else ! p4 massive
                          fh4 = 0
                        endif
                      else ! p1 is massive and not polarized
                        if (m4.le.2) then ! p4 massless
                          if (lp.ne.1.and.binT(w4).eq.2**(legs-1)-1) then
                            ! last tree current has no propagator
                            fh4 = - 1
                          else
                            fh4 = + 1
                          endif
                        else ! p4 massive
                          fh4 = 0
                        endif
                      endif
                    else ! p1 = s or v
                      fh4 = 0
                    endif
                  elseif (p2.le.12) then ! p2 = scalar (not phi+/-)
                    if (m4.le.2) then ! p4 massless
                      if (lp.ne.1.and.binT(w4).eq.2**(legs-1)-1) then
                        ! last tree current has no propagator
                      fh4 = + fheT(w1)
                      else
                      fh4 = - fheT(w1)
                      endif
                    else ! p4 massive
                      fh4 = 0
                    endif
                  elseif (ck1.eq.'f'.and.ck2.eq.'f~') then ! p1 = f, p2 = f~
                    ! p4 = s or v for sure
                    if (p4.eq.18.or.p4.eq.19) then ! p4 = W
                      if ( fheT(w1).eq.+1.or. &      ! p1 right-handed
                           fheT(w2).eq.-1     ) then ! p2 left-handed
                        fhfilter(j) = .false.
                        fh4 = 3
                      else ! other cases are ok
                        fh4 = 0
                      endif
                    elseif (p4.ge.15.and.p4.le.17) then ! p4 = vector (not W)
                      if ( fheT(w1)*fheT(w2).ne.0 .and. &
                           fheT(w1)+fheT(w2).ne.0 ) then
                        fhfilter(j) = .false.
                        fh4 = 3
                      else ! other cases are ok
                        fh4 = 0
                      endif
                    elseif (p4.le.14) then ! p4 = scalar
                      ! Here we use that in the coupling f f'~ s
                      ! the omega_+ part is proportional to m_f and
                      ! the omega_- part is proportional to m_f'~.
                      if ( fheT(w1)*fheT(w2).ne.0 ) then
                        fhfilter(j) = .false.
                        fh4 = 3
                      else ! other cases are ok
                        fh4 = 0
                      endif
                    else ! p4 = ? (just to be sure, it can be omitted)
                      fh4 = 0
                    endif
                  elseif (ck2.eq.'f'.and.ck1.eq.'f~') then ! p2= f, p1= f~
                    ! p4 = s or v for sure
                    if (p4.eq.18.or.p4.eq.19) then ! p4 = W
                      if ( fheT(w2).eq.+1.or. &      ! p2 right-handed
                           fheT(w1).eq.-1     ) then ! p1 left-handed
                        fhfilter(j) = .false.
                        fh4 = 3
                      else ! other cases are ok
                        fh4 = 0
                      endif
                    elseif (p4.ge.15.and.p4.le.17) then ! p4 = vector (not W)
                      if ( fheT(w1)*fheT(w2).ne.0 .and. &
                           fheT(w1).ne.-fheT(w2) ) then
                        fhfilter(j) = .false.
                        fh4 = 3
                      else ! other cases are ok
                        fh4 = 0
                      endif
                    elseif (p4.le.14) then ! p4 = scalar
                      ! Here we use that in the coupling f f'~ s
                      ! the omega_+ part is proportional to m_f and
                      ! the omega_- part is proportional to m_f'~.
                      if ( fheT(w1)*fheT(w2).ne.0 ) then
                        fhfilter(j) = .false.
                        fh4 = 3
                      else ! other cases are ok
                        fh4 = 0
                      endif
                    else ! p4 = ? (just to be sure, it can be omitted)
                      fh4 = 0
                    endif
                  else ! all other cases (if any) are ok
                    fh4 = 0
                  endif

                  if (fh4.eq.3) then
                    if (fheT(w4).eq.111) fheT(w4) = 3
                  elseif (fh4.eq.-1.and.(fheT(w4).eq.0.or.fheT(w4).eq.+1)) then
                    ! should never be the case, just to be sure
                    fheT(w4) = 0
                  elseif (fh4.eq.+1.and.(fheT(w4).eq.0.or.fheT(w4).eq.-1)) then
                    ! should never be the case, just to be sure
                    fheT(w4) = 0
                  else
                    fheT(w4) = fh4
                  endif

                  if (binT(w4).eq.2**(legsE-1)-1) then ! last current
                    if ( fheT(w4)*fheT(legsE).ne.0 .and. &
                         fheT(w4).ne.fheT(legsE) ) then
                      fhfilter(j) = .false.
                      if (fheT(w4).eq.111) fheT(w4) = 3
                    endif
                  endif

                  if (helicity_optimization .eq. 0 ) then
                    fhfilter(j) = .true.
                    fheT(w4) = 111
                  endif

                enddo

                fhfilter1: do j = bMaxT,1,-1 ! reverse
                  if (.not.jfilter(j)) then
                    fhfilter(j) = .false.
                    cycle
                  endif
                  w = wT(4,j)
                  if ( (fhfilter(j)).and.(binT(w).ne.lexp-1) ) then
                    check = 1
                    do k = bMaxT, j+1, -1
                      if (fhfilter(k)) then
                        w1 = wT(1,k)
                        w2 = wT(2,k)
                        w3 = wT(3,k)
                        if (daj(k) .gt. 0) then
                          check = 1
                        else if ((w-w1) .eq. 0 .or. (w-w2) .eq. 0 .or. &
                                 (w-w3) .eq. 0) then
                          check = 0
                        else
                          check = 1
                        end if
                        if (moj(j).ne.0) &
                          check = check*(moj(j)-daj(k))
                        if (check.eq.0) cycle fhfilter1
                      endif
                    enddo
                    if (check.ne.0) fhfilter(j) = .false.
                  endif
                enddo fhfilter1

                if (sav.eq.rep) then

                  ! fix the map wtTOw0b for the external particles
                  ! of this lp, case and config
                  allocate (wtTOw0b(wMaxT))
                  do it = 1,legs ! "it" is the wT of the
                                 ! external particles
                    do w0 = 1, w0eTot(pr)
                      check = + abs(parw0e(w0,pr)-park(it))  &
                              + abs(binw0e(w0,pr)-2**(it-1)) &
                              + abs(helw0e(w0,pr)-helk(it))  &
                              + sum(abs(csw0(-1:legsE,w0)-csT(-1:legsE,it)))
                      if (check.eq.0) wtTOw0b (it) = w0
                    enddo
                  enddo

                  ! Some variables to help computing the counters
                  ! for currents and branches
                  startw0 = 0
                  ew0 = 0
                  if (sav.eq.1) then
                    startbm0 = 0
                    startbd0 = 0
                    ebm0 = 0
                    ebd0 = 0
                    if (lp.eq.1) then
                      allocate(startbm1(0:tlm1-1))
                      allocate(startbd1(0:tlm1-1))
                      allocate(ebm1(0:tlm1-1))
                      allocate(ebd1(0:tlm1-1))
                      startbm1 = 0
                      startbd1 = 0
                      ebm1 = 0
                      ebd1 = 0
                    endif
                  endif

                endif

                ! Initialize w1lTot
                if (lp.eq.1.and.rep.eq.1) w1lTot(:,i) = 0

                allocate (jfirstwT(wMaxT))
                jfirstwT = 0

                ! Initialization of currents.
                ! The order of the computation of the currents is:
                !  - non-self-energy mother currents
                !  - non-self-energy daughter currents
                !  - self-energy mother currents
                !  - self-energy daughter currents
                allocate (jinitwT(wMaxT))
                allocate (winitj(bMaxT))
                jinitwT = 0
                do j= 1, bMaxT
                  if (jfilter(j).and.fhfilter(j)) then
                    w4  = wT(4,j)
                    if (jinitwT(w4).eq.0) then
                      jinitwT(w4) = j
                      winitj(j) = .true.
                    else
                      ! lda  = T <=> j is daugther
                      ! lse  = T <=> j is a self-energy
                      ! ldai = T <=> jinitwT(w4) is daugther
                      ! lsei = T <=> jinitwT(w4) is a self-energy
                      lda  = daj(j).ne.0
                      lse  = binT(wT(1,j)).eq.binT(wT(4,j))
                      ldai = daj(jinitwT(w4)).ne.0
                      lsei = binT(wT(1,jinitwT(w4))).eq.binT(wT(4,jinitwT(w4)))
                      if ( &
                         ((.not.lda.and..not.lse).and.(ldai.or.lsei)) .or. &
                         ((     lda.and..not.lse).and.(lsei))         .or. &
                         ((.not.lda.and.     lse).and.(ldai.and.lsei))     &
                         ) then
                        winitj(jinitwT(w4)) = .false.
                        jinitwT(w4) = j
                        winitj(j) = .true.
                      else
                        winitj(j) = .false.
                      endif
                    endif
                  endif
                enddo

                ! Loop over the branches as created by skeleton
                ! (rerun for each helicity configuration).
                ! Here will be fixed the "b" (from "branch") global
                ! variables which depend on the polarization

                jloop2: do j = 1, bMaxT

                  if (jfilter(j).and.fhfilter(j)) then

                    w1 = wT(1,j)
                    w2 = wT(2,j)
                    w3 = wT(3,j)
                    w4 = wT(4,j)

                    e1 = binT(w1)
                    if (w2.eq.0) then
                      e2 = 0
                    else
                      e2 = binT(w2)
                    endif
                    if (w3.eq.0) then
                      e3 = 0
                    else
                      e3 = binT(w3)
                    endif
                    e4 = binT(w4)

                    ! For loop only.
                    ! h1 is the off-set of the second propagator,
                    ! which, because of the rules to avoid double
                    ! counting, is always an odd number (the off-set
                    ! of the first propagator is always 0). ih1 is
                    ! the corresponding counter. ih1 is then used to
                    ! optimize counting the currents and reduces the
                    ! memory needed for them.
                    if (e4.ge.2**legs) then
                      h1 = 0
                      do k = 1,size(hosT,1)
                        if (hosT(k,w4)+1.eq.2) h1 = h1 + 2**(k-1)
                      enddo
                      if (h1.eq.0) then
                        ih1 = 0
                      else
                        ih1 = (h1+1)/2
                      endif
                    endif

                    ! Given a local current w4, jfirstwT(w4) is
                    ! the first local branch j which has w4 as
                    ! outgoing current.
                    if (jfirstwT(w4).eq.0) jfirstwT(w4) = j

                    ! Here we compute the counting number of the
                    ! outgoing currents, for all branches except the
                    ! last ones

                    if (e4.ne.2**(legsE-1)-1) then

                      if (e4.lt.2**legs.and.sav.eq.rep) then ! tree branch

                        ! If the outgoing current of present branch
                        ! has the same properties of an already
                        ! computed current of a previous lp, t, fh or
                        ! configs,  the current number w0Tot is not
                        ! increased and wtTOw0b is fixed to the value
                        ! of w0Tot of that current. In that case, no
                        ! current needs to be computed and the branch
                        ! j is skipped. This is not done for last
                        ! branches, which are always computed.
                        ! First we define peakleg(k): if the external
                        ! leg k would contribute to the outgoing
                        ! current of this branch (this can be seen
                        ! from the binary tag of the current),
                        ! peakleg(k) is set to 1, otherwise it is 0.
                        do k = 1,legsE; kbin = 2**(k-1)
                          peakleg(k) = 0
                          do kk = 1,legsE
                            if (kbin.eq.firstNumbers(e4,kk)) &
                               peakleg(k) = 1
                          enddo
                        enddo
                        check0 = 1
                        lp0step = max(lp,1)
                        ! Just the tree level (lp0=0) and the present
                        ! level (lp0=lp) are considered
                        do lp0 = 0,lp,lp0step
                          do t0 = 1,t
                            do fh0 = fhmax(t0),fhmin(t0),-2
                              if ( lp0.eq.lp .and. &
                                    t0.eq.t  .and. &
                                   fh0.eq.fh       ) then
                                i0max = i-1
                              else
                                i0max = cfTot(pr)
                              endif
                              do i0 = 1,i0max
                                ! Here the helicities of the external
                                ! particles contributing to the old
                                ! and present currents are checked
                                peak = 0
                                do k = 1,legs
                                  peak =                                 &
                                  peak + peakleg(k) *                    &
                                         abs(heli(k,i,pr)-heli(k,i0,pr))
                                enddo
                                ! An old lp0, t0 (cut), fh0, i0
                                ! (config) can have produced the
                                ! same current as the present branch.
                                ! Now all currents for that lp0, t0,
                                ! fh0, i0 are scanned to see, whether
                                ! one of them has the same properties
                                ! as the outgoing current of the
                                ! present branch.
                                if ( (peak.eq.0) .and. &
                                     (w0min(ej(j),i0,fh0,t0,lp0).ne.0) ) then
                                  w0loop: do w0 = w0min(ej(j),i0,fh0,t0,lp0), &
                                                  w0max(ej(j),i0,fh0,t0,lp0)
                                    if (parT(w4).ne.parw0(w0)) cycle w0loop
                                    if (binT(w4).ne.binw0(w0)) cycle w0loop
                                    if (xxxT(w4).ne.xxxw0(w0)) cycle w0loop
                                    if ( gsT(w4).ne. gsw0(w0)) cycle w0loop
                                    do ii = -1,legs
                                      if (csT(ii,w4).ne.csw0(ii,w0)) cycle w0loop
                                    enddo
                                    if (U1gT(w4).neqv.U1gw0(w0)) cycle w0loop
                                    if (ffT(w4).ne.ffw0(w0)) cycle w0loop
                                    check0 = 0
                                    wtTOw0b(w4) = w0
                                    cycle jloop2
                                  enddo w0loop
                                endif
                              enddo
                            enddo
                          enddo
                        enddo
                        ! If we are here, it means that the current
                        ! was not already computed. The current
                        ! number w0Tot and wtTOw0b have to be fixed
                        ! and the current needs to be computed.
                        if (jfirstwT(w4).eq.j) then
                          ! If we are here, it means the current is
                          ! new and must be computed (and w0Tot is
                          ! increased).
                          w0Tot(pr) = w0Tot(pr) + 1
                          if (w0Tot(pr).gt.w0Def) call reallocate_w0Def
                          wtTOw0b(w4) = w0Tot(pr)
                          parw0(w0Tot(pr)) = parT(w4)
                          binw0(w0Tot(pr)) = binT(w4)
                          csw0(-1:legsE,w0Tot(pr)) = csT(-1:legsE,w4)
                          xxxw0(w0Tot(pr)) = xxxT(w4)
                           gsw0(w0Tot(pr)) =  gsT(w4)
                          U1gw0(w0Tot(pr)) = U1gT(w4)
                           ffw0(w0Tot(pr)) =  ffT(w4)
                          ! min and max currents
                          if (startw0.eq.0) then
                            w0min(ej(j),i,fh,t,lp) = w0Tot(pr)
                            startw0 = 1
                          elseif (ej(j).gt.ew0) then
                            w0max(ew0,i,fh,t,lp) = w0Tot(pr) - 1
                            w0min(ej(j),i,fh,t,lp) = w0Tot(pr)
                          endif
                          w0max(ej(j),i,fh,t,lp) = w0Tot(pr)
                          ew0 = ej(j)
                        endif

                      elseif (e4.ge.2**legs) then ! loop branch

                        check0 = 1
                        do k = 1,legsE; kbin = 2**(k-1)
                          peakleg(k) = 0
                          do kk = 1,legsE
                            if (kbin.eq.firstNumbers(e4,kk)) &
                              peakleg(k) = 1
                          enddo
                        enddo
                        ! An old i0 (config), with same lp=1, t, fh
                        ! can have produced the same current as the
                        ! present branch. Since lp, t and fh are the
                        ! same, the current index w4 is the same for
                        ! all i; therefore in order to check whether
                        ! the current has already been computed, one
                        ! must have peak=0 and check whether for the
                        ! old i0 the current w4 has been computed.
                        ! This is what okwti does.
                        do i0 = 1,i-1
                          peak = 0
                          do k = 1,legs
                            peak =                                 &
                            peak + peakleg(k) *                    &
                                   abs(heli(k,i,pr)-heli(k,i0,pr))
                          enddo
                          if ( (peak.eq.0) .and. okwti(w4,i0)) then
                            check0 = 0
                            if (sav.eq.1.and.rep.eq.0) then
                              if (wtTOw1b(w4,i0).eq.0) then
                                w1sTot(ih1) = w1sTot(ih1) + 1
                                wtTOw1b(w4, i) = w1sTot(ih1)
                                wtTOw1b(w4,i0) = w1sTot(ih1)
                              else
                                wtTOw1b(w4,i) = wtTOw1b(w4,i0)
                              endif
                            endif
                            cycle jloop2
                          endif
                        enddo
                        if (check0.ne.0.and.jfirstwT(w4).eq.j) then
                          if (rep.eq.0) then
                            okwti(w4,i) = .true.
                          else
                            if (wtTOw1b(w4,i).eq.0) then
                              w1lTot(ih1,i) = w1lTot(ih1,i) + 1
                              wtTOw1b(w4,i) = &
                              w1sTot(ih1) + w1lTot(ih1,i)
                            endif
                          endif
                        endif

                      endif

                    endif

                    ! Here we define the indices "c0Eff" and "cEff",
                    ! just for those values which produce currents,
                    ! to save memory. "c0Eff" combines for trees the
                    ! three indices lp, t and fh. "cEff" combines for
                    ! loops the three indices t, fh and ih1.
                    if (sav.eq.rep) then
                      if (e4.lt.2**legs) then ! tree branch
                        if (.not.comp0(fh,t,lp)) then
                          comp0(fh,t,lp) = .true.
                          c0EffMax(pr) = c0EffMax(pr) + 1
                          c0Eff(fh,t,lp) = c0EffMax(pr)
                        endif
                        c = c0Eff(fh,t,lp)
                      else ! loop branch
                        if (.not.comp1(ih1,fh,t)) then
                          comp1(ih1,fh,t) = .true.
                          cEffMax(pr) = cEffMax(pr) + 1
                          cEff(ih1,fh,t) = cEffMax(pr)
                        endif
                        c = cEff(ih1,fh,t)
                      endif
                    endif

                    if (sav.eq.0.and.rep.eq.0) then

                      ! Counters b are increased
                      if (e4.lt.2**legs) then ! tree branch
                        if (daj(j).eq.0) then ! mother branch
                          bm0 = bm0 + 1
                          bm0Tot = max(bm0,bm0Tot)
                        else ! daughter branch
                          bd0 = bd0 + 1
                          bd0Tot = max(bd0,bd0Tot)
                        endif
                      else ! loop branch
                        if (daj(j).eq.0) then ! mother branch
                          bm1 = bm1 + 1
                          bm1Tot = max(bm1,bm1Tot)
                        else ! daughter branch
                          bd1 = bd1 + 1
                          bd1Tot = max(bd1,bd1Tot)
                        endif
                      endif

                    elseif (sav.eq.1.and.rep.eq.0) then

                      if (e4.ge.2**legs) then ! loop branch
                        if (daj(j).eq.0) then ! mother branch
                          bm1h1(ih1,fh,t) = bm1h1(ih1,fh,t) + 1
                        else ! daughter branch
                          bd1h1(ih1,fh,t) = bd1h1(ih1,fh,t) + 1
                        endif
                      endif

                    elseif (sav.eq.1.and.rep.eq.1) then

                      if (e4.lt.2**legs) then ! tree branch
                        c0TOlp(c,pr) = lp
                      else                    ! loop branch
                        cTOt  (c,pr) = t
                        cTOfh (c,pr) = fh
                        cTOih1(c,pr) = ih1
                      endif

                      if (e4.lt.2**legs) then ! tree branch
                        if (daj(j).eq.0) then ! mother branch
                          bm0 = bm0 + 1
                          b = bm0
                          bm0Tot = max(bm0,bm0Tot)
                        else ! daughter branch
                          bd0 = bd0 + 1
                          b = bd0
                          bd0Tot = max(bd0,bd0Tot)
                        endif
                      else ! loop branch
                        if (daj(j).eq.0) then ! mother branch
                          bm1 = bm1 + 1
                          bm1Tot = max(bm1,bm1Tot)
                          bm1c(c) = bm1c(c) + 1
                          b = bm1c(c)
                          b = b + sum(bm1h1(0:ih1-1,fh,t))
                          if (fh.eq.-1) b = b + sum(bm1h1(:,1,t))
                          b = b + sum(bm1h1(:,:,1:t-1))
                        else ! daughter branch
                          bd1 = bd1 + 1
                          bd1Tot = max(bd1,bd1Tot)
                          bd1c(c) = bd1c(c) + 1
                          b = bd1c(c)
                          b = b + sum(bd1h1(0:ih1-1,fh,t))
                          if (fh.eq.-1) b = b + sum(bd1h1(:,1,t))
                          b = b + sum(bd1h1(:,:,1:t-1))
                        endif
                      endif

                      if (e4.ge.2**legs) &
                        rankwMax(ih1) = max(rankwMax(ih1),rankInj(j))

                      ! Interaction type
                      ty = typeT(j)
                      if     (typeT(j).eq.5) then  ! s = f + {f}
                        if (fheT(w1).eq.-1.and.fheT(w2).eq.-1) ty = 7
                        if (fheT(w1).eq.+1.and.fheT(w2).eq.+1) ty = 9
                      elseif (typeT(j).eq.6) then  ! s = {f} + f
                        if (fheT(w1).eq.-1.and.fheT(w2).eq.-1) ty = 8
                        if (fheT(w1).eq.+1.and.fheT(w2).eq.+1) ty = 10
                      elseif (typeT(j).eq.15) then ! v = f + {f}
                        if (fheT(w1).eq.-1.and.fheT(w2).eq.+1) ty = 17
                        if (fheT(w1).eq.+1.and.fheT(w2).eq.-1) ty = 19
                      elseif (typeT(j).eq.16) then ! v = {f} + f
                        if (fheT(w1).eq.+1.and.fheT(w2).eq.-1) ty = 18
                        if (fheT(w1).eq.-1.and.fheT(w2).eq.+1) ty = 20
                      elseif (typeT(j).eq.22) then ! f = f + s, massless
                        if (fheT(w1).eq.-1) ty = 25
                        if (fheT(w1).eq.+1) ty = 27
                      elseif (typeT(j).eq.24) then ! f = s + f, massless
                        if (fheT(w2).eq.-1) ty = 26
                        if (fheT(w2).eq.+1) ty = 28
                      elseif (typeT(j).eq.32) then ! f = f + v, massless
                        if (fheT(w1).eq.-1) ty = 35
                        if (fheT(w1).eq.+1) ty = 37
                      elseif (typeT(j).eq.34) then ! f = v + f, massless
                        if (fheT(w2).eq.-1) ty = 36
                        if (fheT(w2).eq.+1) ty = 38
                      elseif (typeT(j).eq.42) then ! {f} = {f} + s, massless
                        if (fheT(w1).eq.-1) ty = 45
                        if (fheT(w1).eq.+1) ty = 47
                      elseif (typeT(j).eq.44) then ! {f} = s + {f}, massless
                        if (fheT(w1).eq.+1) ty = 46
                        if (fheT(w1).eq.-1) ty = 48
                      elseif (typeT(j).eq.52) then ! {f} = {f} + v, massless
                        if (fheT(w1).eq.-1) ty = 55
                        if (fheT(w1).eq.+1) ty = 57
                      elseif (typeT(j).eq.54) then ! {f} = v + {f}, massless
                        if (fheT(w1).eq.+1) ty = 56
                        if (fheT(w1).eq.-1) ty = 58
                      elseif (typeT(j).eq.1022) then ! f = f, massless
                        if (fheT(w1).eq.-1) ty = 1023
                        if (fheT(w1).eq.+1) ty = 1024
                      elseif (typeT(j).eq.1042) then ! {f} = {f}, massless
                        if (fheT(w1).eq.-1) ty = 1043
                        if (fheT(w1).eq.+1) ty = 1044
                      endif

                      ! Now we fix the value of the global variables
                      ! for "b" branches

                      if (e4.lt.2**legs) then ! tree branch

                        if (daj(j).eq.0) then ! mother branch

                          ! Map from b to s branches
                          sbm0(b,pr) = sm0j(j)

                          ! Incomming and outgoing currents
                          do k = 1,3
                            if (wT(k,j).ne.0) then
                              w0inbm0(k,b,pr) = wtTOw0b(wT(k,j))
                            endif
                          enddo
                          if (e4.ne.2**(legsE-1)-1) &
                             w0outbm0(b,pr) = wtTOw0b(w4)

                          ! The current is initialized in this branch?
                          winitbm0(b,pr) = winitj(j)

                          ! Interaction type
                          typebm0(b,pr) = ty

                          ! min and max branches
                          if (startbm0.eq.0) then
                            bm0min(ej(j),i,c,pr) = b
                            startbm0 = 1
                          elseif (ej(j).gt.ebm0) then
                            bm0max(ebm0,i,c,pr) = b - 1
                            bm0min(ej(j),i,c,pr) = b
                          endif
                          bm0max(ej(j),i,c,pr) = b
                          ebm0 = ej(j)

                        else ! daughter branch

                          ! Map from b to s branches
                          sbd0(b,pr) = sd0j(j)

                          ! Outgoing current
                          if (e4.ne.2**(legsE-1)-1) &
                             w0outbd0(b,pr) = wtTOw0b(w4)

                          ! The current is initialized in this branch?
                          winitbd0(b,pr) = winitj(j)

                          ! min and max branches
                          if (startbd0.eq.0) then
                            bd0min(ej(j),i,c,pr) = b
                            startbd0 = 1
                          elseif (ej(j).gt.ebd0) then
                            bd0max(ebd0,i,c,pr) = b - 1
                            bd0min(ej(j),i,c,pr) = b
                          endif
                          bd0max(ej(j),i,c,pr) = b
                          ebd0 = ej(j)

                        endif

                      else ! loop branch

                        if (daj(j).eq.0) then ! mother branch

                          ! Map from b to s branches
                          sbm1(b,pr) = sm1j(j)

                          ! Incomming and outgoing currents
                          w1inbm1(b,pr) = wtTOw1b(wT(1,j),i)
                          do k = 2,3
                            if (wT(k,j).ne.0) then
                              w0inbm1(k,b,pr) = wtTOw0b(wT(k,j))
                            endif
                          enddo
                          if (e4.ne.2**(legsE-1)-1) &
                             w1outbm1(b,pr) = wtTOw1b(w4,i)

                          ! The current is initialized in this branch?
                          winitbm1(b,pr) = winitj(j)

                          ! Interaction type
                          typebm1(b,pr) = ty

                          ! min and max branches
                          if (startbm1(ih1).eq.0) then
                            ! bm1min(ej(j),i,c,pr) = b
                            bm1_b(pr)%conf(ej(j),i,c)%bmin = b
                            startbm1(ih1) = 1
                          elseif (ej(j).gt.ebm1(ih1)) then
                            ! bm1max(ebm1(ih1),i,c,pr) = b - 1
                            ! bm1min(ej(j),i,c,pr) = b
                            bm1_b(pr)%conf(ebm1(ih1),i,c)%bmax = b -1
                            bm1_b(pr)%conf(ej(j),i,c)%bmin = b
                          endif
                          ! bm1max(ej(j),i,c,pr) = b
                          bm1_b(pr)%conf(ej(j),i,c)%bmax = b
                          ebm1(ih1) = ej(j)

                        else ! daughter branch

                          ! Map from b to s branches
                          sbd1(b,pr) = sd1j(j)

                          ! Outgoing current
                          if (e4.ne.2**(legsE-1)-1) &
                             w1outbd1(b,pr) = wtTOw1b(w4,i)

                          ! The current is initialized in this branch?
                          winitbd1(b,pr) = winitj(j)

                          ! min and max branches
                          if (startbd1(ih1).eq.0) then
                            ! bd1min(ej(j),i,c,pr) = b
                            bd1_b(pr)%conf(ej(j),i,c)%bmin = b
                            startbd1(ih1) = 1
                          elseif (ej(j).gt.ebd1(ih1)) then
                            ! bd1max(ebd1(ih1),i,c,pr) = b - 1
                            ! bd1min(ej(j),i,c,pr) = b
                            bd1_b(pr)%conf(ebd1(ih1),i,c)%bmax = b - 1
                            bd1_b(pr)%conf(ej(j),i,c)%bmin = b
                          endif
                          ! bd1max(ej(j),i,c,pr) = b
                          bd1_b(pr)%conf(ej(j),i,c)%bmax = b
                          ebd1(ih1) = ej(j)

                        endif

                      endif

                      ! Compute comp0gs and comp1gs
                      if (e4.eq.2**(legsE-1)-1) then ! last branch
                        if (e4.lt.2**legs) then ! tree branch
                          comp0gs(gsT(w4),pr) = .true.
                        else ! loop branch
                          comp1gs(gsT(w4),pr) = .true.
                        endif
                      endif

                    endif

                  endif

                enddo jloop2

                deallocate (winitj)
                deallocate (jinitwT)

                deallocate (jfirstwT)

                if (sav.eq.1.and.rep.eq.1) then
                  if (lp.eq.1) then
                    deallocate(startbm1)
                    deallocate(startbd1)
                    deallocate(ebm1)
                    deallocate(ebd1)
                  endif
                endif

                if (sav.eq.rep) deallocate (wtTOw0b)

              enddo configloop

            enddo reploop

            if (lp.eq.1) then
              deallocate (okwti)
              if (sav.eq.1) then
                deallocate (wtTOw1b)
                do ih1 = 0,tlm1-1
                  if (cEff(ih1,fh,t).gt.0) then
                    w1TotMax(cEff(ih1,fh,t),pr) = &
                    w1sTot(ih1) + maxval(w1lTot(ih1,:))
                    riwMax(cEff(ih1,fh,t),pr) = riMax(rankwMax(ih1))
                    w1Tot(pr) = w1Tot(pr) + w1sTot(ih1) + sum(w1lTot(ih1,:))
                  endif
                enddo
              endif
            endif

            if (lp.eq.1) then
              if (sav.eq.1) then
                deallocate (rankwMax)
              endif
            endif

          enddo fhloop

          deallocate (fhfilter)
          deallocate (fheT)

          if (sav.eq.1) then
            if (lp.eq.1) then
              deallocate (rankInj)
              deallocate (sd1j)
              deallocate (sm1j)
            endif
            deallocate (sd0j)
            deallocate (sm0j)
          endif

          deallocate (ej)

          deallocate (jfilter)

          deallocate (jdaTOjmo)
          deallocate (facj)
          deallocate (daj)
          deallocate (moj)
          deallocate (jwT)
          deallocate (njwT)

          deallocate (  colT)
          deallocate (  couT)
          deallocate (gsIncT)
          deallocate ( typeT)
          deallocate (    wT)
          deallocate (ffT)
          deallocate (U1gT)
          deallocate ( gsT)
          deallocate (xxxT)
          deallocate (hmaT)
          deallocate (hosT)
          deallocate (binT)
          deallocate (csT)
          deallocate (parT)

        enddo cutsloop

        if (lp.eq.1) then
          deallocate (cEff)
          if (sav.eq.1) then
            deallocate (w1lTot)
            deallocate (w1sTot)
          endif
        endif

        deallocate (peakleg)
        deallocate (csign)
        deallocate (fervec)
        deallocate (ferp)

      enddo lploop

      if (loop(pr).and.sav.eq.1) then
        deallocate (bd1c)
        deallocate (bm1c)
        deallocate (bd1h1)
        deallocate (bm1h1)
      endif

      if (loop(pr)) deallocate (comp1)
      deallocate (comp0)
      deallocate (c0Eff)
      deallocate (w0min,w0max,parw0,binw0,csw0,xxxw0,gsw0,U1gw0,ffw0)

      if (writeRAM.ge.2) then
        call openOutput
        write(nx,'(2x,a)') '->  free again'
      endif

      !-------------------------------------------------------------
      ! csIa from pCs
      !-------------------------------------------------------------
      ! pCs's are the primary colour structures.
      ! csIa's are the final colour structures of the amplitude
      !
      ! For colour_optimization < 2, the pCs's coincide with csIa's:
      ! csIa = pCs
      !
      ! For colour_optimization = 2, in the list of pCs's are
      ! excluded those colour structures where the colour indices
      ! of gluon stay in the same delta (i.e. the pCs's do not
      ! contain d(ia_l,id_l) for any gluon l).
      ! The final colour structures csIa's are then built applying
      ! to all primary colour structures pCs's the projector
      !
      ! d(ia_l,id_l')*d(ia_l',id_l) - 1/Nc*d(ia_l,id_l)*d(ia_l',id_l')
      !
      ! for each gluon l
      !-------------------------------------------------------------

      ! ng = number of gluons in the process
      ng = 0
      do k = 1,legs
        if (cftype2(par(k,pr)).eq.'g') ng = ng + 1
      enddo

      allocate (dCs0   (lmax,2**ng))
      allocate (dCs0fac(     2**ng))
      allocate (dCs    (lmax,2**ng))
      allocate (dCsfac (     2**ng))

      csTot(pr) = pCsTot(pr)

      do i = 1,pCsTot(pr)

        ! - csIa = pCs for the structures with labels from 1 to
        !   pCsTot.
        ! - nIa(i) is the number of contributions to the final
        !   structure with label "i" coming from different pCs's.
        ! - pIa(k,i) for k=1,...,nIa(i) is the label of the
        !   primary structure from which is coming the k^th
        !   contribution to structure with label "i".
        ! - facIa(k,i) is the factor by which has to be multiplied
        !   the k^th contribution to structure with label "i".
         csIa(:,i,pr) = pCs(:,i)
          nIa(  i,pr) = 1
          pIa(1,i,pr) = i
        facIa(1,i,pr) = 1d0

        if (colour_optimization.ge.2) then

          ! the first dCs0 is the pCs, with factor 1
          dCs0Tot = 1
          dCs0   (:,1) = pCs(:,i)
          dCs0fac(  1) = 1d0

          do k = 1,legs
            if (cftype2(park(k)).eq.'g') then
              j = dCs0Tot
              do j0 = 1,dCs0Tot
                ! d(k,iqk) <=> k = dCs0(iqk,j0)
                ! d(iak,k) <=> dCs0(k,j0) = iak
                do iqk = 1,legs
                  if (dCs0(iqk,j0).eq.k) exit
                enddo
                iak = dCs0(k,j0)
                j = j + 1
                dCs0   (  :,j) = dCs0(:,j0)
                dCs0   (  k,j) = k
                dCs0   (iqk,j) = iak
                dCs0fac(    j) = - 1/Nc*dCs0fac(j0)
                if (k.eq.iqk) dCs0fac(j) = Nc*dCs0fac(j)
              enddo
              dCs0Tot = j
            endif
          enddo

          dCsTot = 0
          do j0 = 1,dCs0Tot
            check = 1
            do j = 1,dCsTot
              if (sum(abs(dCs(:,j)-dCs0(:,j0))).eq.0) then
                dCsfac(j) = dCsfac(j) + dCs0fac(j0)
                check = 0
                exit
              endif
            enddo
            if (check.ne.0) then
              dCsTot = dCsTot + 1
              dCs   (:,dCsTot) = dCs0   (:,j0)
              dCsfac(  dCsTot) = dCs0fac(  j0)
            endif
          enddo

          do j = 2,dCsTot ! We start from the 2nd dCs
            newcs = .true.
            do k = pCsTot(pr)+1,csTot(pr)
              if (sum(abs(dCs(:,j)-csIa(:,k,pr))).eq.0) then
                newcs = .false.
                nIa(k,pr) = nIa(k,pr) + 1
                n = nIa(k,pr)
                  pIa(n,k,pr) = i
                facIa(n,k,pr) = dCsfac(j)
                exit
              endif
            enddo
            if (newcs) then
              csTot(pr) = csTot(pr) + 1
              k = csTot(pr)
               csIa(:,k,pr) = dCs(:,j)
                nIa(  k,pr) = 1
                pIa(1,k,pr) = i
              facIa(1,k,pr) = dCsfac(j)
            endif
          enddo

        endif

      enddo

      deallocate (dCs0,dCs0fac,dCs,dCsfac)

      if (sav.eq.1) then

        bm0prTot(pr) = bm0Tot
        bm1prTot(pr) = bm1Tot

        !  call openOutput
        !  write(nx,*)
        !  write(nx,*) 'bm0Total =',bm0Tot
        !  write(nx,*) 'bd0Total =',bd0Tot
        !  write(nx,*)
        !if (loop(pr)) then
        !  write(nx,*) 'bm1Total =',bm1Tot
        !  write(nx,*) 'bd1Total =',bd1Tot
        !  write(nx,*)
        !endif
        !  write(nx,*) ' w0Total =',w0Tot(pr)
        !if (loop(pr)) then
        !  write(nx,*) ' w1Total =',w1Tot(pr)
        !  write(nx,*) ' tiTotal =',tiTot
        !  write(nx,*)
        !  write(nx,*) 'c0EffMax =',c0EffMax(pr)
        !  write(nx,*) 'cEffMax  =',cEffMax(pr)
        !endif
        !  write(nx,*)

        if (loop(pr)) then
          ranktiMax = maxval(rankti(0:tiTot(pr),pr))
          ritiMax(pr) = riMax(ranktiMax)
        endif

        if (loop(pr)) then
          step = nmasses + 1
          do t = 1,tiTot(pr)
            hm0 = massti(t)
            do i = lmax,1,-1
              vmti(i,t,pr) = hm0/step**(i-1)
              hm0 = hm0 - vmti(i,t,pr)*step**(i-1)
            enddo
          enddo
          deallocate (massti)
        endif

        !-------------------------------------------------------------
        ! The two representations of the colour structures
        !-------------------------------------------------------------
        ! The amplitude has the form
        !
        !   A(ia_1,...ia_N;iq_1,...,iq_N) =
        !   sum_P d(ia_P(1),iq_1) * ... * d(ia_P(N),iq_N) * A_P.
        !
        ! or equivalentely
        !
        !   A(ia_1,...ia_N;iq_1,...,iq_N) =
        !   sum_Q d(ia_1,iq_Q(N)) * ... * d(ia_L,iq_Q(N)) * A_Q,
        !
        ! N = number of gluons + number of quark-antiquark pairs
        !     (all particle are incoming).
        ! The sums run over all possible permutations P or Q of the
        ! labels.
        !
        ! Being L the  number of legs of the process, the variables
        ! csIa(1:L,P,pr) contain the informations in the first form:
        !
        ! csIa(k,P,pr) = P(k) if k is a gluon or an incoming quark
        !              = 0    otherwise
        !
        ! Here we construct the variables csIq(1:L,Q,pr) which
        ! contain the informations in the second form:
        !
        ! csIq(k,Q,pr) = Q(k) if k is a gluon or an incoming antiquark
        !              = 0    otherwise
        !-------------------------------------------------------------

        do cs = 1,csTot(pr)
          csIq(:,cs,pr) = 0
          do i = 1,legs
            if (csIa(i,cs,pr).ne.0) csIq(csIa(i,cs,pr),cs,pr) = i
          enddo
        enddo

        !-------------------------------------------------------------
        ! Contraction of colour structures for the squared amplitude
        !-------------------------------------------------------------
        ! The squared amplitude is then defined as (we use the first
        ! representation):
        !
        !   A^2 =
        !   fac*sum_(ia_1,...ia_N;iq_1,...,iq_N)
        !       |A(ia_1,...ia_N;iq_1,...,iq_N)|^2 =
        !   fac*sum_(ia_1,...ia_N;iq_1,...,iq_N) sum_P sum_P'
        !       d(ia_P(1),iq_1) * ... * d(ia_P(N),iq_N) *
        !       d(ia_P'(1),iq_1) * ... * d(ia_P'(N),iq_N) *
        !       A_P * A_P'
        !
        ! where fac = factor(pr) contains:
        ! - the symmetry factor for identical outgoing particles
        ! - the spin averaging factor for incoming particles
        ! - the colour averaging factor for incoming particles
        !
        ! N = number of gluons + number of quark-antiquark pairs
        !-------------------------------------------------------------

        ! Here we define new colour structures without 0's
        ! dMax = N
        dMax = cd0sMax(pr); binMax = 2**dMax-1
        allocate (iq(dMax),ia(dMax,csTot(pr)))
        allocate (diq(legs),dia(legs,csTot(pr))); diq = 0; dia = 0
        allocate (iabin(binMax,csTot(pr)))
        allocate (vec(dMax))
        do cs = 1,csTot(pr)
          d = 0
          do k = 1,legs
            if (csIa(k,cs,pr).ne.0) then
              d = d + 1
               iq(d) = k;  ia(d            ,cs) = csIa(k,cs,pr)
              diq(k) = d; dia(csIa(k,cs,pr),cs) = d
            endif
          enddo
          do dbin = 1,binMax
            ! For each colour structure (with index cs) of the
            ! amplitude, which is a product of dMax deltas, we
            ! consider here all possible subproducts of these deltas.
            ! We use binary notation: for example the subproduct
            ! dbin=11=1+2+8 is the product of the first
            ! (pow(1)=2^(1-1)=1), the second (pow(2)=2^(2-1)=2) and
            ! fourth (pow(4)=2^(4-1)=8) deltas.
            ! iabin(cs,dbin) is then the sum of the binary
            ! "ia"-indices of this sub product.
            ! For example, if the subproduct has "ia"-indices 2,3,4,
            ! then the binary "ia"-indices are 2,4,8 and
            ! iabin = pow(2) + pow(3) + pow(4) = 2 + 4 + 8 = 14.
            ! Analogously one should compute also the sum of the
            ! binary "iq"-indices, but this is always equal to dbin,
            ! because the "iq" indices are always in the ascending
            ! order 1,2,3,4,...
            iabin(dbin,cs) = 0
            do d = 1,dMax
              if (vectorLeg(dbin,d).eq.1) &
                iabin(dbin,cs) = iabin(dbin,cs) + 2**(ia(d,cs)-1)
            enddo
          enddo
        enddo

        ! Here we compare all possible pairs of final structures.
        ! These pairs appear in the computation of the squared
        ! amplitude: the two structures (each of them is a product
        ! of deltas) are multiplied and the deltas are then fully
        ! contracted and give a colour coefficient (colcoef). Here
        ! we make this contraction.
        ! In order to do it, we look at all subproducts of the two
        ! structures of the pair. The "iq"-indices are in both
        ! structures in the same ascending order 1,2,3,4,... per
        ! default. The point is to look at the "closed"
        ! subproducts. A subproduct is "closed" if the "ia"-indices
        ! in the first structure are the same (in a different
        ! order) as in the second structure. Example:
        ! str1= d(ia=2,iq=1)*d(ia=3,iq=2)*d(ia=4,iq=3)*d(ia=1,iq=4)
        ! str2= d(ia=3,iq=1)*d(ia=2,iq=2)*d(ia=4,iq=3)*d(ia=1,iq=4)
        ! These two substructures have 3 "closed" subproducts: the
        ! subproduct made of the first two deltas, the subproduct
        ! made of the third delta alone and the subproduct made of
        ! the fourth delta alone.
        ! Each "closed" subproduct contributes to the colour
        ! coefficient with a factor Nc=3.
        ! In order to decide if a subproduct dbin is closed, it is
        ! enough to look at the iabin(i,dbin), iabin(j,dbin) of
        ! the two structures. If they are equal, the product is
        ! closed. In order to avoid double counting, when a "closed"
        ! subproduct dbin is found, the primary binaries composing
        ! it are saved (in vec) and in the next steps just
        ! subproducts not containing these primary binaries are
        ! considered.
        do i1 = 1,csTot(pr)
        do i2 = 1,csTot(pr)
          eq = 0
          vec = 0
          dbinloop1: do dbin = 1,binMax
            do d = 1,dMax
              if (vectorLeg(dbin,d)*vec(d).eq.1) cycle dbinloop1
            enddo
            if (iabin(dbin,i1).eq.iabin(dbin,i2)) then
              eq = eq + 1
              vec(:) = vec(:) + vectorLeg(dbin,1:dMax)
            endif
          enddo dbinloop1
          colcoef(i1,i2,pr)= nCs**eq
        enddo
        enddo

        !-------------------------------------------------------------
        ! Contraction of colour structures for the squared amplitude
        ! for colour correlation
        !-------------------------------------------------------------
        ! Let us write the amplitude in the first representation
        ! using a compact notation:
        !
        !   A = sum_P str_P * A_P
        !
        ! where
        !
        ! str_P = d(ia_P(1),iq_1) * ... * d(ia_P(N),iq_N)
        !
        ! Let us pick two lines l1 and l2 and consider the two
        ! amplitudes
        !
        !   A1 = sum_P  str_P(l1)  * A_P
        !   A2 = sum_P' str_P'(l2) * A_P'
        !
        ! where str_P(l1) is the structure that one obtains from
        ! str_P when line l1 emits an extra gluon; similarly
        ! str_P'(l2) is the structure that one obtains from str_P'
        ! when line l2 emits an extra gluon.
        !
        ! Colour correlation means to compute
        !
        !   |A|_(l1,l2)^2 = c*A1*A2^* =
        !   c * sum_P sum_P' str_P(l1) * str_P'(l2) * A_P^* * A_P'
        !
        ! where c is a conventional factor depending of the type
        ! of particles l1 and l2 (which is choosen, such that colour
        ! correlation is indipendent from the normalization of the
        ! SU(Nc) generators).
        !-------------------------------------------------------------
        ! Conventions:
        !   t^a_(i,j) = 1/sqrt(2)*lambda^a_(i,j)
        ! being lambda^a_(i,j) the entry 'i,j' of the Gellmann
        ! matrix 'a'. The normalization is:
        !   Tr(t^a t^b) = delta(a,b).
        !-------------------------------------------------------------
        ! Let l be the line on which we attach the gluon (l = l1
        ! or l2).
        !
        ! - If l is an incoming quark, in str we have a delta of
        !   the type d(ia_X,iq_l). Adding a gluon with index 'a' to
        !   line l corresponds to contracting with
        !     - t^a_(iq_l,iq_l').
        !   In the colour flow representation adding a gluon with
        !   indices (ia,iq) to quark line l corresponds to
        !   contracting with
        !     str0_q(l,l',ia,iq) =
        !     - t^a_(iq_l,iq_l')*t^a_(ia,iq) =
        !     - d(ia,iq_l')*d(iq_l,iq) + 1/Nc*d(iq_l,iq_l')*d(ia,iq)
        !
        ! - If l is an incoming antiquark, in str we have a delta
        !   of the type d(ia_l,iq_X). Adding a gluon with index 'a' to
        !   line l corresponds to contracting with
        !     t^a_(ia_l',ia_l).
        !   In the colour flow representation adding a gluon with
        !   indices (ia,iq) to antiquark line l corresponds to
        !   contracting with
        !     str0_q~(l,l',ia,iq) =
        !     t^a_(ia_l',ia_l)*t^a_(ia,iq) =
        !     d(ia,ia_l)*d(ia_l',iq) - 1/Nc*d(ia_l',ia_l)*d(ia,iq)
        !
        ! - If l is a gluon, in str we have a delta of the type
        !   d(ia_l,Iq_X)*d(ia_X,iq_l). Adding a gluon with
        !   index 'a' to line l corresponds to contracting with
        !     i*f^{al' a al}.
        !   In the colour flow representation adding a gluon with
        !   indices (ia,iq) to gluon line l corresponds to
        !   contracting with
        !     str0_g(l,l',ia,iq) =
        !     i*f^{al' a al}*
        !     t^al'_(ia_l',iq_l')*t^a_(ia,iq)*t^al_(iq_l,ia_l) =
        !     + d(ia_l',iq)*d(ia,ia_l)*d(iq_l,iq_l')
        !     - d(ia_l',ia_l)*d(ia,iq_l')*d(iq_l,iq)
        !-------------------------------------------------------------
        ! Let (ia,iq) be the indices of the gluon emitted by l1
        ! and (Ia,Iq) the indices of the gluon emitted by l2.
        ! Let the two colour structures be:
        !   str1 = d(ia_P(1) ,iq_1)*...*d(ia_P(N) ,iq_N)
        !   str2 = d(Ia_P'(1),Iq_1)*...*d(Ia_P'(N),Iq_N)
        ! One "iq"-index of str1 is iq_l1 and/or one "ia"-index of
        ! str1 is ia_l1. One "iq"-index of str2 is Iq_l2 and/or one
        ! "ia"-index of str2 is Ia_l2.
        ! Without colour correlation the contraction is done in
        ! this way:
        !   colcoef = str1*str2*
        !             d(iq_1,Iq_1)*d(ia_1,Ia_1)*
        !             ...
        !             d(iq_N,Iq_N)*d(ia_N,Ia_N)
        ! With colour correlation str1 is contracted with
        ! str0_p1(l1,l1',ia,iq), thus the indices ia_l1,iq_l1 are
        ! replaced by ia_l1',iq_l1' and the indices ia,iq appear;
        ! p1 is the particle type (q,q~,g) of line l1.
        ! With colour correlation str2 is contracted with
        ! str0_p2(l2,l2',Ia,Iq), thus the indices Ia_l2,Iq_l1 are
        ! replaced by Ia_l2',Iq_l2' and the indices Ia,Iq appear;
        ! p2 is the particle type (q,q~,g) of line l2.
        ! So contraction is done in this way:
        !   colcoef = str1 * str2 *
        !             str0_p1(l1,l1',ia,iq) *
        !             str0_p2(l2,l2',Ia,Iq) *
        !             d(iq,Iq) * d(ia,Ia) *
        !             d(iq_1,Iq_1) * d(ia_1,Ia_1) *
        !             ...
        !             d(iq_N,Iq_N) * d(ia_N,Ia_N)
        !-------------------------------------------------------------

        do l1 = 1,legs; ck1 = cftype2(park(l1))
        do l2 = 1,legs; ck2 = cftype2(park(l2))

          colcoefc(:,:,l1,l2,pr) = 0d0
          if (l1.eq.l2) cycle
          if (ck1.ne.'q'.and.ck1.ne.'q~'.and.ck1.ne.'g') cycle
          if (ck2.ne.'q'.and.ck2.ne.'q~'.and.ck2.ne.'g') cycle

          do i1 = 1,csTot(pr)
          do i2 = i1,csTot(pr)

            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! q - q
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1 = d(ia_X1,iq_l1) * d(ia_X2,iq_l2) * R1
            !
            ! str1' =
            ! str1 * str0_q(l1,l1',ia,iq) * d(iq_l2',iq_l2) =
            ! R1 * d(ia_X2,iq_l2') * (
            ! - d(ia_X1,iq) * d(ia,iq_l1')
            ! + 1/Nc * d(ia,iq) * d(ia_X1,iq_l1')
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str2 = d(Ia_X1,Iq_l1) * d(Ia_X2,Iq_l2) * R2
            !
            ! str2' =
            ! str2 * str0_q(l2,l2',Ia,Iq) * d(iq_l1',iq_l1) =
            ! R2 * d(Ia_X1,Iq_l1') * (
            ! - d(Ia_X2,Iq) * d(Ia,Iq_l2')
            ! + 1/Nc * d(Ia,Iq) * d(Ia_X2,Iq_l2')
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1' * str2' *
            ! d(iq_l1',Iq_l1') * d(iq_l2',Iq_l2') *
            ! d(iq,Iq) * d(ia,Ia)
            ! =
            ! + R1 * d(ia_X1,Ia_X2) * d(ia_X2,Ia_X1) * R2
            ! - 1/Nc * R1 * d(ia_X1,Ia_X1) * d(ia_X2,Ia_X2) * R2
            ! =
            ! d(iq_l1,Iq_l1) * d(iq_l2,Iq_l2) * R1 * R2 * (
            ! + d(ia_X1,iq_l1) * d(ia_X2,iq_l2) * R1 *
            !   d(Ia_X2,Iq_l1) * d(Ia_X1,Iq_l2) * R2
            ! - 1/Nc * d(ia_X1,iq_l1) * d(ia_X2,iq_l2) * R1 *
            !          d(Ia_X1,Iq_l1) * d(Ia_X2,Iq_l2) * R2
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if     (ck1.eq.'q'.and.ck2.eq.'q') then
              ! Number of contributions
              nMax = 2
              allocate (ia1(dMax,nMax),ia2(dMax,nMax))
              do n = 1,nMax
                 ia1(:,n) = ia(:,i1)
                 ia2(:,n) = ia(:,i2)
              enddo
              ! Overall factor
              facr = 1d0/2d0/Cf
              ! Contribution 1:
              ! str1' = str1
              ! str2' = str2 [Ia_X1 <-> Ia_X2]
              ia2(diq(l1),1) = ia(diq(l2),i2)
              ia2(diq(l2),1) = ia(diq(l1),i2)
              corfac(1) = facr
              ! Contribution 2:
              ! str1' = str1
              ! str2' = str2
              corfac(2) = - 1d0/Nc * facr

            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! q~ - q~
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1 = d(ia_l1,iq_X1) * d(ia_l2,iq_X2) * R1
            !
            ! str1' =
            ! str1 * str0_q~(l1,l1',ia,iq) * d(ia_l2',ia_l2) =
            ! R1 * d(ia_l2',iq_X2) * (
            ! + d(ia_l1',iq) * d(ia,iq_X1)
            ! - 1/Nc * d(ia,iq) * d(ia_l1',iq_X1)
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str2 = d(Ia_l1,Iq_X1) * d(Ia_l2,Iq_X2) * R2
            !
            ! str2' =
            ! str2 * str0_q~(l2,l2',ia,iq) * d(Ia_l1',Ia_l1) =
            ! R2 *  d(Ia_l1',Iq_X1) * (
            ! + d(Ia_l2',Iq) *d(Ia,Iq_X2)
            ! - 1/Nc * d(Ia,Iq) * d(Ia_l2',Iq_X2)
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1' * str2' *
            ! d(ia_l1',Ia_l1') * d(ia_l2',Ia_l2') *
            ! d(iq,Iq) * d(ia,Ia)
            ! =
            ! + R1 * d(Iq_X1,iq_X2) * d(Iq_X2,iq_X1) * R2
            ! - 1/Nc * R1 * d(Iq_X1,iq_X1) * d(Iq_X2,iq_X2) * R2
            ! =
            ! d(ia_l1,Ia_l1) * d(ia_l2,Ia_l2) * (
            ! + d(ia_l1,iq_X1) * d(ia_l2,iq_X2) * R1 *
            !   d(Ia_l2,Iq_X1) * d(Ia_l1,Iq_X2) * R2
            ! - 1/Nc * d(ia_l1,iq_X1) * d(ia_l2,iq_X2) * R1 *
            !          d(Ia_l1,Iq_X1) * d(Ia_l2,Iq_X2) * R2
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            elseif (ck1.eq.'q~'.and.ck2.eq.'q~') then
              ! Number of contributions
              nMax = 2
              allocate (ia1(dMax,nMax),ia2(dMax,nMax))
              do n = 1,nMax
                 ia1(:,n) = ia(:,i1)
                 ia2(:,n) = ia(:,i2)
              enddo
              ! Overall factor
              facr = 1d0/2d0/Cf
              ! Contribution 1:
              ! str1' = str1
              ! str2' = str2 [Ia_l1 <-> Ia_l2]
              ia2(dia(l1,i2),1) = ia(dia(l2,i2),i2)
              ia2(dia(l2,i2),1) = ia(dia(l1,i2),i2)
              corfac(1) = facr
              ! Contribution 2:
              ! str1' = str1
              ! str2' = str2
              corfac(2) = - 1d0/Nc * facr

            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! q - q~
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1 = d(ia_X1,iq_l1) * d(ia_l2,iq_X2) * R1
            !
            ! str1' =
            ! str1 * str0_q(l1,l1',ia,iq) * d(ia_l2',ia_l2) =
            ! R1 * d(ia_l2',iq_X2) * (
            ! - d(ia_X1,iq) * d(ia,iq_l1')
            ! + 1/Nc * d(ia,iq) * d(ia_X1,iq_l1')
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str2 = d(Ia_X1,Iq_l1) * d(Ia_l2,Iq_X2) * R2
            !
            ! str2' =
            ! str2 * str0_q~(l2,l2',ia,iq) * d(iq_l1',iq_l1) =
            ! R2 * d(Ia_X1,Iq_l1') * (
            ! + d(Ia_l2',Iq) * d(Ia,Iq_X2)
            ! - 1/Nc * d(Ia,Iq) * d(Ia_l2',Iq_X2)
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1' * str2' *
            ! d(iq_l1',Iq_l1') * d(ia_l2',Ia_l2') *
            ! d(iq,Iq) * d(ia,Ia)
            ! =
            ! - R1 * d(ia_X1,iq_X2) * d(Ia_X1,Iq_X2) * R2
            ! + 1/Nc * R1 * d(ia_X1,Ia_X1) * d(Iq_X2,iq_X2) * R2
            ! =
            ! - d(iq_l1,Iq_l1) * d(ia_l2,Ia_l2) *
            !   d(iq,Iq) * d(ia,Ia) * (
            !     d(ia,iq) * d(ia_X1,iq_l1) * d(ia_l2,iq_X2) * R1 *
            !     d(Ia_X1,Iq) * d(Ia_l2,Iq_l1) * d(Ia,Iq_X2) * R2
            !   )
            ! + 1/Nc * d(iq_l1,Iq_l1) * d(ia_l2,Ia_l2) * (
            !     d(ia_X1,iq_l1) * d(ia_l2,iq_X2) * R1 *
            !     d(Ia_X1,Iq_l1) * d(Ia_l2,Iq_X2) * R2
            !   )
            ! =
            ! 1/Nc * d(iq_l1,Iq_l1) * d(ia_l2,Ia_l2) * (
            ! - d(ia_l2,iq_l1) * d(ia_X1,iq_X2) * R1 *
            !   d(Ia_l2,Iq_l1) * d(Ia_X1,Iq_X2) * R2
            ! + d(ia_X1,iq_l1) * d(ia_l2,iq_X2) * R1 *
            !   d(Ia_X1,Iq_l1) * d(Ia_l2,Iq_X2) * R2
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            elseif (ck1.eq.'q'.and.ck2.eq.'q~') then
              ! Number of contributions
              nMax = 2
              allocate (ia1(dMax,nMax),ia2(dMax,nMax))
              do n = 1,nMax
                 ia1(:,n) = ia(:,i1)
                 ia2(:,n) = ia(:,i2)
              enddo
              ! Overall factor
              facr = 1d0/2d0/Cf
              ! Contribution 1:
              ! str1' = str1 [ia_X1 <-> ia_l2]
              ! str2' = str2 [Ia_X1 <-> Ia_l2]
              ia1(diq(l1),1) = l2; ia1(dia(l2,i1),1) = ia(diq(l1),i1)
              ia2(diq(l1),1) = l2; ia2(dia(l2,i2),1) = ia(diq(l1),i2)
              corfac(1) = - 1d0/Nc * facr
              if (ia(diq(l1),i1).eq.l2) corfac(1) = corfac(1) * Nc
              if (ia(diq(l1),i2).eq.l2) corfac(1) = corfac(1) * Nc
              ! Contribution 2:
              ! str1' = str1
              ! str2' = str2
              corfac(2) = + 1d0/Nc * facr

            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! q~ - q
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1 = d(ia_l1,iq_X1) * d(ia_X2,iq_l2) * R1
            !
            ! str1' =
            ! R1 * d(ia_X2,iq_l2') * (
            ! str1 * str0_q~(l1,l1',ia,iq) * d(iq_l2',iq_l2) =
            ! + d(ia_l1',iq) * d(ia,iq_X1)
            ! - 1/Nc * d(ia,iq) * d(ia_l1',iq_X1)
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str2 = d(Ia_l1,Iq_X1) * d(Ia_X2,Iq_l2) * R2
            !
            ! str2' =
            ! R2 * d(Ia_l1',Iq_X1) * (
            ! str2 * str0_q(l2,l2',ia,iq) * d(ia_l1',ia_l1) =
            ! - d(Ia_X2,Iq) * d(Ia,Iq_l2')
            ! + 1/Nc * d(Ia,Iq) * d(Ia_X2,Iq_l2')
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1' * str2' *
            ! d(ia_l1',Ia_l1') * d(iq_l2',Iq_l2') *
            ! d(iq,Iq) * d(ia,Ia)
            ! =
            ! - R1 * d(ia_X2,iq_X1) * d(Ia_X2,Iq_X1) * R2
            ! + 1/Nc * R1 * d(ia_X2,Ia_X2) * d(Iq_X1,iq_X1) * R2
            ! =
            ! ...
            ! =
            ! 1/Nc * d(ia_l1,Ia_l1) * d(iq_l2,Iq_l2) * (
            ! - d(ia_X2,iq_X1) * d(ia_l1,iq_l2) * R1 *
            !   d(Ia_X2,Iq_X1) * d(Ia_l1,Iq_l2) * R2
            ! + d(ia_l1,iq_X1) * d(ia_X2,iq_l2) * R1 *
            !   d(Ia_l1,Iq_X1) * d(Ia_X2,Iq_l2) * R2
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            elseif (ck1.eq.'q~'.and.ck2.eq.'q') then
              ! Number of contributions
              nMax = 2
              allocate (ia1(dMax,nMax),ia2(dMax,nMax))
              do n = 1,nMax
                 ia1(:,n) = ia(:,i1)
                 ia2(:,n) = ia(:,i2)
              enddo
              ! Overall factor
              facr = 1d0/2d0/Cf
              ! Contribution 1:
              ! str1' = str1 [ia_l1 <-> ia_X2]
              ! str2' = str2 [Ia_l1 <-> Ia_X2]
              ia1(dia(l1,i1),1) = ia(diq(l2),i1); ia1(diq(l2),1) = l1
              ia2(dia(l1,i2),1) = ia(diq(l2),i2); ia2(diq(l2),1) = l1
              corfac(1) = - 1d0/Nc * facr
              if (ia(diq(l2),i1).eq.l1) corfac(1) = corfac(1) * Nc
              if (ia(diq(l2),i2).eq.l1) corfac(1) = corfac(1) * Nc
              ! Contribution 2:
              ! str1' = str1
              ! str2' = str2
              corfac(2) = + 1d0/Nc * facr

            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! q - g
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1 = d(ia_X1,iq_l1) *
            !        d(ia_X2,iq_l2) * d(ia_l2,iq_X3) * R1
            !
            ! str1' =
            ! str1 * str0_q(l1,l1',ia,iq) *
            !        d(iq_l2',iq_l2) * d(ia_l2',ia_l2)
            ! =
            ! R1 * d(ia_X2,iq_l2') * d(ia_l2',iq_X3) * (
            ! - d(ia,iq_l1' ) * d(iq ,ia_X1 )
            ! + 1/Nc * d(ia,iq) * d(ia_X1,iq_l1')
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str2 = d(Ia_X1,Iq_l1) *
            !        d(Ia_X2,Iq_l2) * d(Ia_l2,Iq_X3) * R2
            !
            ! str2' =
            ! str2 * str0_g(l2,l2',ia,iq) * d(Iq_l1',Iq_l1)
            ! =
            ! R2 * d(Ia_X1,Iq_l1') * (
            ! + d(Ia_l2',Iq) * d(Ia_X2,Iq_l2') * d(Ia,Iq_X3)
            ! - d(Ia_X2,Iq) * d(Ia,Iq_l2') * d(Ia_l2',Iq_X3)
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1' * str2' *
            ! d(iq_l1',Iq_l1') * d(iq_l2',Iq_l2') * d(ia_l2',Ia_l2') *
            ! d(iq,Iq) * d(ia,Ia)
            ! =
            ! - R1 * d(ia_X1,iq_X3) * d(ia_X2,Ia_X2) * d(Ia_X1,Iq_X3) * R2
            ! + R1 * d(ia_X1,Ia_X2) * d(ia_X2,Ia_X1) * d(iq_X3,Iq_X3) * R2
            ! =
            ! ...
            ! =
            ! d(iq_l1,Iq_l1) * d(iq_l2,Iq_l2) * d(ia_l2,Ia_l2) * (
            ! + d(ia_X1,iq_l1) * d(ia_X2,iq_l2) * d(ia_l2,iq_X3) * R1 *
            !   d(Ia_X2,Iq_l1) * d(Ia_X1,Iq_l2) * d(Ia_l2,Iq_X3) * R2
            ! - 1/Nc *
            !   d(ia_l2,iq_l1) * d(ia_X2,iq_l2) * d(ia_X1,iq_X3) * R1 *
            !   d(Ia_l2,Iq_l1) * d(Ia_X2,Iq_l2) * d(Ia_X1,Iq_X3) * R2
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            elseif (ck1.eq.'q'.and.ck2.eq.'g') then
              ! Number of contributions
              nMax = 2
              allocate (ia1(dMax,nMax),ia2(dMax,nMax))
              do n = 1,nMax
                 ia1(:,n) = ia(:,i1)
                 ia2(:,n) = ia(:,i2)
              enddo
              ! Overall factor
              facr = 1d0/2d0/Cf
              ! Contribution 1:
              ! str1' = str1
              ! str2' = str2 [Ia_X1 <-> Ia_X2]
              ia2(diq(l1),1) = ia(diq(l2),i2)
              ia2(diq(l2),1) = ia(diq(l1),i2)
              corfac(1) = facr
              ! Contribution 2:
              ! str1' = str1 [ia_X1 <-> ia_l2]
              ! str2' = str2 [Ia_X1 <-> Ia_l2]
              ia1(diq(l1),2) = l2; ia1(dia(l2,i1),2) = ia(diq(l1),i1)
              ia2(diq(l1),2) = l2; ia2(dia(l2,i2),2) = ia(diq(l1),i2)
              corfac(2) = - 1d0/Nc * facr
              if (ia(diq(l1),i1).eq.l2) corfac(2) = corfac(2) * Nc
              if (ia(diq(l1),i2).eq.l2) corfac(2) = corfac(2) * Nc

            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! q~ - g
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1 = d(ia_l1,iq_X1) *
            !        d(ia_X2,iq_l2) * d(ia_l2,iq_X3) * R1
            !
            ! str1' =
            ! str1 * str0_q~(l1,l1',ia,iq) *
            !        d(iq_l2',iq_l2) * d(ia_l2',ia_l2)
            ! =
            ! R1 * d(ia_X2,iq_l2') * d(ia_l2',iq_X3) * (
            ! + d(ia,iq_X1) * d(ia_l1',iq)
            ! - 1/Nc * d(ia,iq) * d(ia_l1',iq_X1)
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str2 = d(Ia_l1,Iq_X1) *
            !        d(Ia_X2,Iq_l2) * d(Ia_l2,Iq_X3) * R2
            !
            ! str2' =
            ! str2 * str0_g(l2,l2',ia,iq) * d(Ia_l1',Ia_l1)
            ! =
            ! R2 * d(Ia_l1',Iq_X1) * (
            ! + d(Ia_l2',Iq) * d(Ia_X2,Iq_l2') * d(Ia,Iq_X3)
            ! - d(Ia_X2,Iq) * d(Ia,Iq_l2') * d(Ia_l2',Iq_X3)
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1' * str2' *
            ! d(ia_l1',Ia_l1') * d(iq_l2',Iq_l2') * d(ia_l2',Ia_l2') *
            ! d(iq,Iq) * d(ia,Ia)
            ! =
            ! + R1 * d(iq_X1,Iq_X3) * d(ia_X2,Ia_X2) * d(iq_X3,Iq_X1) * R2
            ! - R1 * d(ia_X2,iq_X1) * d(Ia_X2,Iq_X1) * d(iq_X3,Iq_X3) * R2
            ! =
            ! ...
            ! =
            ! d(ia_l1,Ia_l1) * d(iq_l2,Iq_l2) * d(ia_l2,Ia_l2) * (
            ! + d(ia_l1,iq_X1) * d(ia_X2,iq_l2) * d(ia_l2,iq_X3) * R1 *
            !   d(Ia_l2,Iq_X1) * d(Ia_X2,Iq_l2) * d(Ia_l1,Iq_X3) * R2
            ! - 1/Nc *
            !   d(ia_X2,iq_X1) * d(ia_l1,iq_l2) * d(ia_l2,iq_X3) * R1 *
            !   d(Ia_X2,Iq_X1) * d(Ia_l1,Iq_l2) * d(Ia_l2,Iq_X3) * R2
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            elseif (ck1.eq.'q~'.and.ck2.eq.'g') then
              ! Number of contributions
              nMax = 2
              allocate (ia1(dMax,nMax),ia2(dMax,nMax))
              do n = 1,nMax
                 ia1(:,n) = ia(:,i1)
                 ia2(:,n) = ia(:,i2)
              enddo
              ! Overall factor
              facr = 1d0/2d0/Cf
              ! Contribution 1:
              ! str1' = str1
              ! str2' = str2 [Ia_l1 <-> Ia_l2]
              ia2(dia(l1,i2),1) = ia(dia(l2,i2),i2)
              ia2(dia(l2,i2),1) = ia(dia(l1,i2),i2)
              corfac(1) = facr
              ! Contribution 2:
              ! str1' = str1 [ia_l1 <-> ia_X2]
              ! str2' = str2 [Ia_l1 <-> Ia_X2]
              ia1(dia(l1,i1),2) = ia(diq(l2),i1); ia1(diq(l2),2) = l1
              ia2(dia(l1,i2),2) = ia(diq(l2),i2); ia2(diq(l2),2) = l1
              corfac(2) = - 1d0/Nc * facr
              if (ia(diq(l2),i1).eq.l1) corfac(2) = corfac(2) * Nc
              if (ia(diq(l2),i2).eq.l1) corfac(2) = corfac(2) * Nc

            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! g - q
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1 = d(ia_X1,iq_l1) * d(ia_l1,iq_X2) *
            !        d(ia_X3,iq_l2) * R1
            !
            ! str1' =
            ! str1 * str0_g(l1,l1',ia,iq) * d(iq_l2',iq_l2)
            ! =
            ! R1 * d(ia_X3,iq_l2') * (
            ! + d(ia_l1',iq) * d(ia_X1,iq_l1') * d(ia,iq_X2)
            ! - d(ia_X1,iq) * d(ia,iq_l1') * d(ia_l1',iq_X2)
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str2 = d(Ia_X1,Iq_l1) * d(Ia_l1,Iq_X2) *
            !        d(Ia_X3,Iq_l2) * R2
            !
            ! str2' =
            ! str2 * str0_q(l2,l2',ia,iq) *
            !        d(Iq_l1',Iq_l1) * d(Ia_l1',Ia_l1)
            ! =
            ! R2 * d(Ia_X1,Iq_l1') * d(Ia_l1',Iq_X2) * (
            ! - d(Ia,Iq_l2') * d(Ia_X3,Iq)
            ! + 1/Nc * d(Ia,Iq) * d(Ia_X3,Iq_l2')
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1' * str2' *
            ! d(iq_l1',Iq_l1') * d(ia_l1',Ia_l1') * d(iq_l2',Iq_l2') *
            ! d(iq,Iq) * d(ia,Ia)
            ! =
            ! - R1 * d(ia_X1,Ia_X1) * d(ia_X3,iq_X2) * d(Ia_X3,Iq_X2) * R2
            ! + R1 * d(ia_X1,Ia_X3) * d(iq_X2,Iq_X2) * d(ia_X3,Ia_X1) * R2
            ! =
            ! ...
            ! =
            ! d(iq_l1,Iq_l1) * d(ia_l1,Ia_l1) * d(iq_l2,Iq_l2) * (
            ! + d(ia_X1,iq_l1) * d(ia_l1,iq_X2) * d(ia_X3,iq_l2) * R1 *
            !   d(Ia_X3,Iq_l1) * d(Ia_l1,Iq_X2) * d(Ia_X1,Iq_l2) * R2
            ! - 1/Nc *
            !   d(ia_X1,iq_l1) * d(ia_X3,iq_X2) * d(ia_l1,iq_l2) * R1 *
            !   d(Ia_X1,Iq_l1) * d(Ia_X3,Iq_X2) * d(Ia_l1,Iq_l2) * R2
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            elseif (ck1.eq.'g'.and.ck2.eq.'q') then
              ! Number of contributions
              nMax = 2
              allocate (ia1(dMax,nMax),ia2(dMax,nMax))
              do n = 1,nMax
                 ia1(:,n) = ia(:,i1)
                 ia2(:,n) = ia(:,i2)
              enddo
              ! Overall factor
              facr = 1d0/2d0/Ca
              ! Contribution 1:
              ! str1' = str1
              ! str2' = str2 [Ia_X1 <-> Ia_X3]
              ia2(diq(l1),1) = ia(diq(l2),i2)
              ia2(diq(l2),1) = ia(diq(l1),i2)
              corfac(1) = facr
              ! Contribution 2:
              ! str1' = str1 [ia_X3 <-> ia_l1]
              ! str2' = str2 [Ia_X3 <-> Ia_l1]
              ia1(diq(l2),2) = l1; ia1(dia(l1,i1),2) = ia(diq(l2),i1)
              ia2(diq(l2),2) = l1; ia2(dia(l1,i2),2) = ia(diq(l2),i2)
              corfac(2) = - 1d0/Nc * facr
              if (ia(diq(l2),i1).eq.l1) corfac(2) = corfac(2) * Nc
              if (ia(diq(l2),i2).eq.l1) corfac(2) = corfac(2) * Nc

            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! g - q~
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1 = d(ia_X1,iq_l1) * d(ia_l1,iq_X2) *
            !        d(ia_l2,iq_X3) * R1
            !
            ! str1' =
            ! str1 * str0_g(l1,l1',ia,iq) * d(ia_l2',ia_l2)
            ! =
            ! R1 * d(ia_l2',iq_X3) * (
            ! + d(ia_l1',iq) * d(ia_X1,iq_l1') * d(ia,iq_X2)
            ! - d(ia_X1,iq) * d(ia,iq_l1') * d(ia_l1',iq_X2)
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str2 = d(Ia_X1,Iq_l1) * d(Ia_l1,Iq_X2) *
            !        d(Ia_l2,Iq_X3) * R2
            !
            ! str2' =
            ! str2 * str0_q~(l2,l2',ia,iq) *
            !        d(Iq_l1',Iq_l1) * d(Ia_l1',Ia_l1)
            ! =
            ! R2 * d(Ia_X1,Iq_l1') * d(Ia_l1',Iq_X2) * (
            ! + d(Ia,Iq_X3) * d(Ia_l2',Iq)
            ! - 1/Nc * d(Ia,Iq) * d(Ia_l2',Iq_X3)
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1' * str2' *
            ! d(iq_l1',Iq_l1') * d(ia_l1',Ia_l1') * d(ia_l2',Ia_l2') *
            ! d(iq,Iq) * d(ia,Ia)
            ! =
            ! + R1 * d(ia_X1,Ia_X1) * d(iq_X2,Iq_X3) * d(iq_X3,Iq_X2) * R2
            ! - R1 * d(ia_X1,iq_X3) * d(iq_X2,Iq_X2) * d(Ia_X1,Iq_X3) * R2
            ! =
            ! ...
            ! =
            ! d(iq_l1,Iq_l1) * d(ia_l1,Ia_l1) * d(ia_l2,Ia_l2) * (
            ! + d(ia_X1,iq_l1) * d(ia_l1,iq_X2) * d(ia_l2,iq_X3) * R1
            !   d(Ia_X1,Iq_l1) * d(Ia_l2,Iq_X2) * d(Ia_l1,Iq_X3) * R2
            ! - 1/Nc *
            !   d(ia_l2,iq_l1) * d(ia_l1,iq_X2) * d(ia_X1,iq_X3) * R1 *
            !   d(Ia_l2,Iq_l1) * d(Ia_l1,Iq_X2) * d(Ia_X1,Iq_X3) * R2
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            elseif (ck1.eq.'g'.and.ck2.eq.'q~') then
              ! Number of contributions
              nMax = 2
              allocate (ia1(dMax,nMax),ia2(dMax,nMax))
              do n = 1,nMax
                 ia1(:,n) = ia(:,i1)
                 ia2(:,n) = ia(:,i2)
              enddo
              ! Overall factor
              facr = 1d0/2d0/Ca
              ! Contribution 1:
              ! str1' = str1
              ! str2' = str2 [Ia_l1 <-> Ia_l2]
              ia2(dia(l1,i2),1) = ia(dia(l2,i2),i2)
              ia2(dia(l2,i2),1) = ia(dia(l1,i2),i2)
              corfac(1) = facr
              ! Contribution 2:
              ! str1' = str1 [ia_X1 <-> ia_l2]
              ! str2' = str2 [Ia_X1 <-> Ia_l2]
              ia1(diq(l1),2) = l2; ia1(dia(l2,i1),2) = ia(diq(l1),i1)
              ia2(diq(l1),2) = l2; ia2(dia(l2,i2),2) = ia(diq(l1),i2)
              corfac(2) = - 1d0/Nc * facr
              if (ia(diq(l1),i1).eq.l2) corfac(2) = corfac(2) * Nc
              if (ia(diq(l1),i2).eq.l2) corfac(2) = corfac(2) * Nc

            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! g - g
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1 = d(ia_X1,iq_l1) * d(ia_l1,iq_X2) *
            !        d(ia_X3,iq_l2) * d(ia_l2,iq_X4) * R1
            !
            ! str1' =
            ! str1 * str0_g(l1,l1',ia,iq) *
            !        d(iq_l2',iq_l2) * d(ia_l2',ia_l2)
            ! =
            ! R1 * d(ia_X3,iq_l2') * d(ia_l2',iq_X4) * (
            ! + d(ia_l1',iq) * d(ia_X1,iq_l1') * d(ia,iq_X2)
            ! - d(ia_X1,iq) * d(ia,iq_l1') * d(ia_l1',iq_X2)
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str2 = d(Ia_X1,Iq_l1) * d(Ia_l1,Iq_X2) *
            !        d(Ia_X3,Iq_l2) * d(Ia_l2,Iq_X4) * R2
            !
            ! str2' =
            ! str2 * str0_g(l2,l2',ia,iq) *
            !        d(Iq_l1',Iq_l1) * d(Ia_l1',Ia_l1)
            ! =
            ! R2 * d(Ia_X1,Iq_l1') * d(Ia_l1',Iq_X2) * (
            ! + d(Ia_l2',Iq) * d(Ia_X3,Iq_l2') * d(Ia,Iq_X4)
            ! - d(Ia_X3,Iq) * d(Ia,Iq_l2') * d(Ia_l2',Iq_X4)
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! str1' * str2' *
            ! d(iq_l1',Iq_l1') * d(ia_l1',Ia_l1') *
            ! d(iq_l2',Iq_l2') d(ia_l2',Ia_l2') *
            ! d(iq,Iq) * d(ia,Ia)
            ! =
            ! + R1 * d(ia_X1,Ia_X1) * d(iq_X2,Iq_X4) *
            !   d(ia_X3,Ia_X3) * d(iq_X4,Iq_X2) * R2
            ! - R1 * d(ia_X1,Ia_X1) * d(ia_X3,iq_X2) *
            !   d(Ia_X3,Iq_X2) * d(iq_X4,Iq_X4) * R2
            ! - R1 * d(ia_X1,iq_X4) * d(iq_X2,Iq_X2) *
            !   d(ia_X3,Ia_X3) * d(Ia_X1,Iq_X4) * R2
            ! + R1 * d(ia_X1,Ia_X3) * d(iq_X2,Iq_X2) *
            !   d(ia_X3,Ia_X1) * d(iq_X4,Iq_X4) * R2
            ! =
            ! ...
            ! =
            ! d(iq_l1,Iq_l1) * d(ia_l1,Ia_l1) *
            ! d(iq_l2,Iq_l2) * d(ia_l2,Ia_l2) * (
            ! + d(ia_X1,iq_l1) * d(ia_l1,iq_X2) *
            !   d(ia_X3,iq_l2) * d(ia_l2,iq_X4) * R1 *
            !   d(Ia_X1,Iq_l1) * d(Ia_l2,Iq_X2) *
            !   d(Ia_X3,Iq_l2) * d(Ia_l1,Iq_X4) * R2
            ! + d(ia_X1,iq_l1) * d(ia_l1,iq_X2) *
            !   d(ia_X3,iq_l2) * d(ia_l2,iq_X4) * R1 *
            !   d(Ia_X3,Iq_l1) * d(Ia_l1,Iq_X2) *
            !   d(Ia_X1,Iq_l2) * d(Ia_l2,Iq_X4) * R2
            ! - 1/Nc *
            !   d(ia_X1,iq_l1) * d(ia_X3,iq_X2) *
            !   d(ia_l1,iq_l2) * d(ia_l2,iq_X4) * R1 *
            !   d(Ia_X1,Iq_l1) * d(Ia_X3,Iq_X2) *
            !   d(Ia_l1,Iq_l2) * d(Ia_l2,Iq_X4) * R2
            ! - 1/Nc *
            !   d(ia_l2,iq_l1) * d(ia_l1,iq_X2) *
            !   d(ia_X3,iq_l2) * d(ia_X1,iq_X4) * R1 *
            !   d(Ia_l2,Iq_l1) * d(Ia_l1,Iq_X2) *
            !   d(Ia_X3,Iq_l2) * d(Ia_X1,Iq_X4) * R2
            ! )
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            elseif (ck1.eq.'g'.and.ck2.eq.'g') then
              ! Number of contributions
              nMax = 4
              allocate (ia1(dMax,nMax),ia2(dMax,nMax))
              do n = 1,nMax
                 ia1(:,n) = ia(:,i1)
                 ia2(:,n) = ia(:,i2)
              enddo
              ! Overall factor
              facr = 1d0/2d0/Ca
              ! Contribution 1:
              ! str1' = str1
              ! str2' = str2 [Ia_l1 <-> Ia_l2]
              ia2(dia(l1,i2),1) = ia(dia(l2,i2),i2)
              ia2(dia(l2,i2),1) = ia(dia(l1,i2),i2)
              corfac(1) = facr
              ! Contribution 2:
              ! str1' = str1
              ! str2' = str2 [Ia_X1 <-> Ia_X3]
              ia2(diq(l1),   2) = ia(diq(l2),i2)
              ia2(diq(l2),   2) = ia(diq(l1),i2)
              corfac(2) = facr
              ! Contribution 3:
              ! str1' = str1 [ia_l1 <-> ia_X3]
              ! str2' = str2 [Ia_l1 <-> Ia_X3]
              ia1(dia(l1,i1),3) = ia(diq(l2),i1)
              ia1(diq(l2),   3) = l1
              ia2(dia(l1,i2),3) = ia(diq(l2),i2)
              ia2(diq(l2),   3) = l1
              corfac(3) = - 1d0/Nc * facr
              if (ia(diq(l2),i1).eq.l1) corfac(3) = corfac(3) * Nc
              if (ia(diq(l2),i2).eq.l1) corfac(3) = corfac(3) * Nc
              ! Contribution 4:
              ! str1' = str1 [ia_X1 <-> ia_l2]
              ! str2' = str2 [Ia_X1 <-> Ia_l2]
              ia1(diq(l1)   ,4) = l2
              ia1(dia(l2,i1),4) = ia(diq(l1),i1)
              ia2(diq(l1)   ,4) = l2
              ia2(dia(l2,i2),4) = ia(diq(l1),i2)
              corfac(4) = - 1d0/Nc * facr
              if (ia(diq(l1),i1).eq.l2) corfac(4) = corfac(4) * Nc
              if (ia(diq(l1),i2).eq.l2) corfac(4) = corfac(4) * Nc

            endif

            do n = 1,nMax
              eq = count_cycles(n,legs)
              colcoefc(i1,i2,l1,l2,pr) = &
              colcoefc(i1,i2,l1,l2,pr) + corfac(n) * Nc**eq
            enddo

            deallocate (ia1,ia2)
            if (i1 .ne. i2) then
              colcoefc(i2,i1,l1,l2,pr) = colcoefc(i1,i2,l1,l2,pr)
            end if
          enddo
          enddo
        enddo
        enddo

        deallocate (iq,ia,diq,dia,iabin,vec)

      endif

      deallocate (pCs)

      if (sav.eq.1) then
        call texFeet(pr.eq.prTot)
      endif

      deallocate (helk,ffe,park)
      deallocate (qflowT)

    enddo prloop

      if (writeRAM.ge.2) then
        call openOutput
        write(nx,*)
      endif

    if (sav.eq.1) then

      cfMax = maxval(cfTot)
      csMax = maxval(csTot)

      do pr = 1,prTot

        legs = legsIn(pr) + legsOut(pr)
        gsTot(0,pr) = legs - 2
        gs2Tot(0,pr) = 2*gsTot(0,pr)
        if (lpmax(pr).gt.0) then
          gsTot(1,pr) = legs
          if (zeroLO(pr)) then
            gs2Tot(1,pr) = 2*gsTot(1,pr)
          else
            gs2Tot(1,pr) = gsTot(0,pr) + gsTot(1,pr)
          endif
        else
           gsTot(1,pr) = 0
          gs2Tot(1,pr) = 0
        endif

        gsTotEff(0,pr) = - 2
        do i = 1,legs
          select case (par(i,pr))
          case (15,23,24,25,29,30,31,35,36,37,41,42,43)
            gsTotEff(0,pr) = gsTotEff(0,pr) + 1
          end select
        enddo
        if ( gsTotEff(0,pr).lt. 0 ) then
          gsTotEff(:,pr) = 0
        else
          gsTotEff(1,pr) = gsTotEff(0,pr) + 2
        endif
        gs2TotEff(0,pr) = 2*gsTotEff(0,pr)
        if (zeroLO(pr)) then
          gs2TotEff(1,pr) = 2*gsTotEff(1,pr)
        else
          gs2TotEff(1,pr) = gsTotEff(0,pr) + gsTotEff(1,pr)
        endif
        if (lpmax(pr).le.0) then
           gsTotEff(1,pr) = 0
          gs2TotEff(1,pr) = 0
        endif

      enddo

       gsMax = maxval( gsTot)
      gs2Max = maxval(gs2Tot)

      if (writeRAM.ge.1) then
        ram2 = 0
        ram2 = ram2 + csMax*(gsMax+1)*cfMax*5*prTot*16 ! matrix
        ram2 = ram2 + 4*maxval(w0Tot)*16 ! ww0
        ram2 = ram2 + 4*(maxval(modaTot)+1)*16 ! ww0out
        if (loopMax) then
          ram2 = ram2 + maxval((ritiMax+1)*tiTot)*16 ! TIri
          ram2 = ram2 + 4*maxval((ritiMax+1)*(modaTot+1))*16 ! ww1out
          ram2 = ram2 + 4*maxval((riwMax+1)*(w1TotMax+1))*16 ! ww1
        endif
        ramMb = int(real(ram2,kind=dp)/1d6) + 2
        call openOutput
        write(nx,'(1x,a,i8,a)') &
         'RAM permanently used by process computation in this run:', &
         ramMb,' Mbytes'
        write(nx,*)

      endif

    endif

  enddo saveloop

  deallocate (wMax,bMax)

  deallocate (noquarks,nogluons,noweaks)

  contains

!---------------------------------------------------------------------

  subroutine reallocate_w0Def

  integer              :: l,wdef0,wdef
  integer, allocatable :: tempI(:,:)
  logical, allocatable :: tempL(:,:)

  wdef0 = w0Def
  wdef  = w0Def + w0Inc

  if (writeRAM.ge.2) then
    ram0 = ram0 + w0Inc*( 4*(legsE+6) +1*1 )
    ramMb = int(real(ram0,kind=dp)/1d6) + 2
    call openOutput
    write(nx,'(2x,a)') '->  free again'
    write(nx,'(1x,a,i3,a,i8,a)',advance='no') &
      'RAM temporally used by generation of process', &
      inpr(pr),':',ramMb,' Mbytes'
  endif

  l = size(csw0,1)-2;  allocate (tempI(-1:l,wdef0))
    tempI = csw0; deallocate (csw0)
    allocate (csw0(-1:l,wdef)); csw0(:,:wdef0) = tempI(:,:)
  deallocate (tempI)

  allocate (tempI(1,wdef0))
    tempI(1,:) = parw0(:); deallocate (parw0)
    allocate (parw0(wdef)); parw0(:wdef0) = tempI(1,:)
    tempI(1,:) = binw0(:); deallocate (binw0)
    allocate (binw0(wdef)); binw0(:wdef0) = tempI(1,:)
    tempI(1,:) = xxxw0(:); deallocate (xxxw0)
    allocate (xxxw0(wdef)); xxxw0(:wdef0) = tempI(1,:)
    tempI(1,:) = gsw0(:); deallocate (gsw0)
    allocate (gsw0(wdef)); gsw0(:wdef0) = tempI(1,:)
    tempI(1,:) = ffw0(:); deallocate (ffw0)
    allocate (ffw0(wdef)); ffw0(:wdef0) = tempI(1,:)
  deallocate (tempI)

  allocate (tempL(1,wdef0))
    tempL(1,:) = U1gw0(:); deallocate (U1gw0)
    allocate (U1gw0(wdef)); U1gw0(:wdef0) = tempL(1,:)
  deallocate (tempL)

  w0Def = wdef

  end subroutine reallocate_w0Def

  function count_cycles(n,legs) result(cycles)
    ! Computes the number of cycles in the arrays ia1 and ia2.
    ! For example
    ! ia1 = [2, 3, 1, 4]
    ! ia2 = [1, 2, 3, 4]
    ! has 2 cycles, namely 2->3->1->2 and 4->4.
    integer, intent(in) :: n,legs
    integer             :: cycles,d,dtmp,covered,ia1tmp,ia1tmp2,dia2(legs)

    ! invert ia2
    do d = 1, dMax
      dia2(ia2(d,n)) = d
    end do

    cycles = 0
    covered = 0

    ! loop over ia1
    dloop: do d = 1, dMax
      ia1tmp = ia1(d,n)
      dtmp = dia2(ia1tmp)

      ! have we counted this cycle already?
      if (iand(2**(dtmp-1),covered) .ne. 0) then
        cycle dloop
      end if

      ! follow the cycle until closed
      ia1tmp2 = ia1(dtmp,n)
      do while (ia1tmp2 .ne. ia1tmp)
        dtmp = dia2(ia1tmp2)
        ia1tmp2 = ia1(dtmp,n)
        covered = ior(covered,2**(dtmp-1))
      end do
      covered = ior(covered,2**(d-1))
      cycles = cycles + 1
    end do dloop

  end function count_cycles

!---------------------------------------------------------------------

  end subroutine generate_currents

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module currents_rcl

!#####################################################################


