      subroutine squared_matrix(model,ng,ecms,p1,p2,pbeam1,pbeam2,qph,
     >           ie,imtx,summt2,prodotto,deltah)
      implicit double precision (a-h,o-z)
      character*10 model
      character*6   ord
      character*6 eikng ! multip single
      parameter (imaxph = 40)
      integer nout,nout1
      parameter (nout   =  3)
      parameter (nout1  =  2)
      dimension iclose(40)
      dimension iclose2(2)
      double precision mass,lambda
      dimension p1(0:3),p2(0:3),pin1(0:3),pin2(0:3),qmost(0:3)
      dimension p1r(0:3),p2r(0:3),p1b(0:3),p2b(0:3),q(0:3),qf(0:3)
      dimension q1(0:3),q2(0:3),q3(0:3),q4(0:3),qq(0:3),qr(0:3),qb(0:3)
      dimension p1tmp(0:3),p2tmp(0:3),p3tmp(0:3),p4tmp(0:3)
      dimension pin1bl(0:3),pin2bl(0:3),p1bl(0:3),p2bl(0:3)
      dimension pin1b1(0:3),pin2b1(0:3),p1b1(0:3),p2b1(0:3)
      dimension pin1bt(0:3),pin2bt(0:3),p1bt(0:3),p2bt(0:3),qbt(0:3)
      dimension q1bt(0:3),newmom(0:3)
      common/maptotworad/pin1bt,pin2bt,p1bt,p2bt,qbt,q1bt
      common/counterg5/counterg4
      integer counterg4
c$$$      common/check_4g/matrix_check,eik_jac,summt2ex
c$$$      real*8 matrix_check
      double precision newmom
      dimension pin1cl(0:3),pin2cl(0:3),p1cl(0:3),p2cl(0:3),qcl(4)
      dimension pin1cl1(4),pin2cl1(4),p1cl1(4),p2cl1(4),qclu(4,2)
      dimension pphot(4,2),idphtme(2),pclout1(4,2)
      dimension pclb(4,0:3),pcla(4,0:3)
      dimension ptmp(0:3),ptmp2(0:3), q2r(0:3),q1r(0:3),pmapout(4,3)
      dimension pmapout1(4,2),qph2(imaxph,0:3)
      dimension pclin(4,2),pclout(4,3),pmapin(4,2),pchrg(4,4)
      double precision massesin(2),massesout(3),mdot,massesout1(2)
      external mdot
      external dot
      double precision dot
      dimension pbeam1(0:3),pbeam2(0:3)
      dimension qmostpt(imaxph,0:3)
      real*8 allun,qfixrad(0:3),qqtwo(0:3)
      real*8 qphr(imaxph,0:3), qtmp(imaxph,0:3)
      dimension p1rr(0:3),p2rr(0:3),pin1rr(0:3),pin2rr(0:3),prif(0:3)
      dimension pin1br(0:3),pin2br(0:3),p1br(0:3),p2br(0:3),qbr(0:3)
      dimension pin1brr(0:3),pin2brr(0:3),p1brr(0:3),p2brr(0:3)
      dimension qbrr(0:3)
      dimension pin1tm(0:3),pin2tm(0:3),p1tm(0:3),p2tm(0:3)
      dimension qph(imaxph,0:3),qphtmp(imaxph,0:3),randvec(3)
      dimension qphtest(imaxph,0:3)
      real*8 hard
      logical condrad
      integer l,phi
      character*2 fs,fsbck
      real*8 tmphardest(0:3)
      common/parameters/ame,ammu,convfac,alpha,pi
      common/qedORDER/ord
      common/finalstate/fs
      integer iffpi,iffpibck,vals,u,error
      dimension vals(2)
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi
      common/iclosest/iclose
      common/momentainitial/pin1,pin2
      common/reducedtoborn/p1b,p2b,iref
      common/momentainitialred/pin1b(0:3),pin2b(0:3),pin1r(0:3),
     >     pin2r(0:3)
      common/forborncrosssection/phsp2b,flux2b,bornme,bornmeq2,bornmez,
     >     bornmenovp
      common/bornxsecpp/bornmenovppp,bornmepp
      common/radpattern/nph(4)
      common/hparam/ampi
      common/epssoft/eps
      common/fortests/prod
      common/for_debug/ppp(0:3),denphsp,dir3(0:3),dir4(0:3),aabr,aabq,
     > aaonm
      common/singlefactors/singleeikfact(imaxph)
      common/idarkon/idarkon
      common/ialpharunning/iarun
      double complex vpolc
      external vpolc

!     imtx is to say that the matrix el. has been calculated..
! summt2 is the squared ME summed over all all spins

      common/radoption/ifisrfsr
      integer ifisrfsr

      integer iteubn
      common/teubner/iteubn
      common/scaling/scaling
*      for recola ee -> 4g
      dimension p6(0:3),p7(0:3)
      common/p6for4g/p6,p7

***   to check weird events I need some commons
*
c      real*8 deltah,deltai
      common/momentacheck1/pin1br,pin2br,p1br,p2br,qbr
      common/fh/fhmax
      double precision fhmax

** VP UNC
      dimension iswvperr(3),iswvperrt(6)
      common/iswitchvperror/iswvperr,iswvperrt
      parameter (nvps = 26)
      dimension remtx(nvps)
      common/rationsvpunc/remtx,invp
      common/rescaledhaderr/scal_err_dhad

      parameter (nvpst = 728)
      dimension remtxt(nvpst)
      integer ist(nvpst,6)
      common/rationsvpunct/remtxt,invpt,ist
***
      common/firstsqmatrix/ifirst
      data ifirst /1/
*** for mesmer matrix element
      double precision QRe,QRmu,muemuegg
      external muemuegg
      common/QReQRmu/QRe,QRmu

      if (ifirst.eq.1) then
         iii = 0
         do i1 = -1,1
         do i2 = -1,1
         do i3 = -1,1
         do i4 = -1,1
         do i5 = -1,1
         do i6 = -1,1
            if (i1.ne.0.or.i2.ne.0.or.i3.ne.0.or.
     >           i4.ne.0.or.i5.ne.0.or.i6.ne.0) then
               iii = iii + 1
               ist(iii,1) = i1
               ist(iii,2) = i2
               ist(iii,3) = i3
               ist(iii,4) = i4
               ist(iii,5) = i5
               ist(iii,6) = i6
            endif
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         ifirst = 0
      endif

      iarunbck = iarun

      imtx  = 0
      summt2   = 0.d0
      prodotto = 1.d0
      deltah   = 0.d0
      bornme   = 1.d-13
      bornmeq2 = 1.d-13
      bornmez  = 0.d0!1.d-13
      pcla     =0d0
      pclb     =0d0
      if (ie.gt.0) return

      imtx = 1
      nspinconf = 2**(4+ng)

      if(idarkon.eq.0) then
         condrad=(fs.eq.'mr'.or.fs.eq.'er'.or.fs.eq.'pr')
         if(.not.condrad) then
            if (ng.eq.0) then
               do k = 0,3
                  pin1b(k) = pin1(k)
                  pin2b(k) = pin2(k)
                  p1b(k) = p1(k)
                  p2b(k) = p2(k)
               enddo
            else
               call closest(ng,qph,pin1,pin2,p1,p2,q1,q2,q3,q4,iclose)
               do k = 0,3
                  p1tmp(k) = pin1(k) - q1(k)
                  p2tmp(k) = pin2(k) - q2(k)
                  qq(k) = p1tmp(k) + p2tmp(k)
               enddo

               if(fs.eq.'ee') then
                  call maptozero2(p1tmp,p2tmp,p1,p2,q1,q2,q3,q4,
     >                 pin1b,pin2b,p1b,p2b,iemap)
               else
                  call maptozero(p1tmp,p2tmp,p1,p2,
     >                 pin1b,pin2b,p1b,p2b,iemap)
               endif
            endif
            iq2   = 0
            iarun = 0  ! added ! ON OFF

            call born_me(pin1b,pin2b,p1b,p2b,bornmenovp,bornmeq2,
     &           vuoto,iq2,q)

            iarun = iarunbck
            iq2 = 1

            bornmenovppp=bornmenovp

            call born_me(pin1b,pin2b,p1b,p2b,bornme,bornmeq2,
     &           bornmez,iq2,q)

            bornmepp=bornme

         endif
***********************RADIATIVI
         if(condrad) then

            if(ng.eq.1)
     >           then
               do k = 0,3
                  pin1br(k) = pin1(k)
                  pin2br(k) = pin2(k)
                  p1br(k) = p1(k)
                  p2br(k) = p2(k)
                  qbr(k)  = qph(1,k)
               enddo
               call check_invariants(pin1br,pin2br,p1br,p2br,qbr)
            else
**********************MAPPING MAURO MORETTI
               call closest(ng,qph,pin1,pin2,p1,p2,q1,q2,q3,q4,iclose)

               ind=1
               call getsignalphot(pin1,pin2,p1,p2,qph,ng,ind,j)
********FROM NOW ON THE SIGNAL PHOTON IS CALLED J

               pin1cl=pin1
               pin2cl=pin2

               p1cl=p1
               p2cl=p2

               do k = 1,ng
*     * PS QUANTITIES
                  if(k.ne.j) then
!                     print*,iclose(k)
                     if (iclose(k).eq.1) then
                        pin1cl=pin1cl-qph(k,:)
                     endif
                     if (iclose(k).eq.2) then
                        pin2cl=pin2cl-qph(k,:)
                     endif
                     if (iclose(k).eq.3) then !
                        p1cl=p1cl+qph(k,:)
                     endif
                     if (iclose(k).eq.4) then
                        p2cl=p2cl+qph(k,:)
                     endif
                  endif

               enddo


               ind=1
               call maptoone_rad(pin1cl,pin2cl,p1cl,p2cl,qph(j,:),
     >              iclose(j),ind,pin1br,pin2br,p1br,p2br,qbr)

               call check_invariants(pin1br,pin2br,p1br,p2br,qbr)
               if(iscoll.eq.1.or.fscoll.eq.1) print*,'EVENTO BALORDO'
            endif
         endif
*[[[ START VP PARAMETRIC UNCERTAINTY
         IF (abs(scal_err_dhad).gt.1d-5) THEN
            if (iteubn.eq.0) then
               remtx = 1.d0
               invp = 0
               iswvperr(1) = 0
               iswvperr(2) = 0
               iswvperr(3) = 0
               do i1 = -1,1
               do i2 = -1,1
               do i3 = -1,1
                  if (i1.ne.0.or.i2.ne.0.or.i3.ne.0) then
                     invp = invp + 1
                     iswvperr(1) = i1
                     iswvperr(2) = i2
                     iswvperr(3) = i3
                     call born_me(pin1b,pin2b,p1b,p2b,xr,bq2
     .                    ,vuoto,iq2,q)
                     remtx(invp) = xr/bornme
                  endif
               enddo
               enddo
               enddo
            else  ! for hlmnt VP parametrization
               remtxt    = 1.d0
               invpt     = 0
               iswvperrt = 0
               do invpt=1,nvpst
                  iswvperrt(1) = ist(invpt,1)
                  iswvperrt(2) = ist(invpt,2)
                  iswvperrt(3) = ist(invpt,3)
                  iswvperrt(4) = ist(invpt,4)
                  iswvperrt(5) = ist(invpt,5)
                  iswvperrt(6) = ist(invpt,6)
                  call born_me(pin1b,pin2b,p1b,p2b,xr,bq2,vuoto,iq2,q)
                  remtxt(invpt) = xr/bornme
               enddo
            endif
         ENDIF
*]]]  END VP PARAMETRIC UNCERTAINTY
      endif
*******************************************************
*** DARK MATTER ON - TREE LEVEL (cmq 1 g in uscita) ***
*******************************************************
      if(idarkon.eq.1) then
         print*,'SHOULD NOT PASS HERE! Stopping...'
         stop
      endif
*
** calcola cose che servono x il matching
*
      flux2b = 8.d0*pin1b(0)*pin2b(0)
      mass = pin1b(0)+pin2b(0)
      amel = ame
      if (fs.eq.'gg') amel = 0.d0
      if (fs.eq.'mm'.or.fs.eq.'mr') amel = ammu
      if (fs.eq.'pp'.or.fs.eq.'pr') amel = ampi
      if (fs.eq.'pi') amel = ampi
      sqla = sqrt(lambda(mass**2,amel**2,amel**2))
      phsp2b = sqla/8.d0/mass**2 /(2.d0*pi) / 4.d0 ! spins...

****
*********************
**    MATCHING     **
*********************

! doing if only the k-th photon was present
!     Here i put the exact ME also for 'ps' and 'psrad' in order to correctly unweight
      if (model.eq.'matched'.or.model.eq.'ps'.or.model.eq.'psrad') then
         iarunbck = iarun
         if (ng.eq.0) then
**********************
* AGGIUNTO PER I RADIATIVI
***********************
            if(fs.eq.'mr'.or.fs.eq.'er'.or.fs.eq.'pr') then
               print*, '***********************************************'
               print*,'NON PUOI AVERE 0 FOTONI IN UN EVENTO RADIATIVO'
               print*,' ***********************************************'
               return
            else
               iref = 1
               summt2 = bornme
               imtx = 1
               iarun = iarunbck
               return
            endif
         endif

!         stop

************
c QUI NG È L'EVENTO A NG FOTONI, QUINDI NON CAMBIA NULLA SE SONO NEL RADIATIVO O MENO, È SEMPRE UN EVENTO A 1 FOTONE -> BK_ME

******* make mm and ee as if it's mr/er!
         fsbck    = fs
         extraden = 1.d0
         if (fs.eq.'mm') then
            extraden = 1.d0/ng
            fs = 'mr'
         endif
         if (fs.eq.'ee') then
            extraden = 1.d0/ng
            fs = 'er'
         endif
*******
         if (ng.eq.1) then
            do i=0,3
               q(i) = qph(1,i)
            enddo

            if(fs.ne.'mr'.and.fs.ne.'er'
     >		.and.fs.ne.'pr') then
c     * bornme already accounts for VP.
c     * I switch it off here for the matching to avoid jumps.
c     * It's anyway a higher order effect
               iarun = 0
               call bk_me(pin1,pin2,p1,p2,q,summt2,qph)
               summt2 = summt2/bornmenovp * bornme! added

               iarun = iarunbck
               call bk_me(pin1,pin2,p1,p2,q,summt2,qph)
c               summt2 = summt2/bornmenovp * bornme! added
c               call glueALPHA(ng,pin1,pin2,p1,p2,qph(1:ng,:),summt2a)
c               print*,summt2,summt2a,summt2a/summt2-1.d0
c               stop

c               fs = fsbck
c              call antenna_approx_ee(ng,pin1,pin2,p1,p2,qph(1,:),summt2)
c              call antenna_approx_gm(ng,pin1,pin2,p1,p2,qph(1,:),summt2)

               
               imtx = 1
               iarun = iarunbck
               fs = fsbck
               return
            endif

            if(fs.eq.'mr'.or.fs.eq.'er'
     >           .or.fs.eq.'pr') then
               pin1br=pin1
               pin2br=pin2
               p1br=p1
               p2br=p2
               qbr=q
               call bk_me(pin1,pin2,p1,p2,q,summt2,qph)
c               call glueALPHA(ng,pin1,pin2,p1,p2,qph(1:ng,:),summt2a)
c               print*,summt2a/summt2-1.d0,summt2a,summt2
c               stop

c               fs = fsbck
c              call antenna_approx_ee(ng,pin1,pin2,p1,p2,qph(1,:),summt2)
c              call antenna_approx(ng,pin1,pin2,p1,p2,qph(1,:),summt2)
c              call antenna_approx_gm(ng,pin1,pin2,p1,p2,qph(1,:),summt2)
               
               summt2=summt2
               imtx = 1
               fs = fsbck
              return
            endif
    	 endif
!    	 print*,iarun

         if(ng.eq.2) then

c      summt2 = 0.5d0*eeeegg(pin2,pin1,p2,p1,qph(1,:),qph(2,:))
c      return
c      call glueALPHA(ng,pin1,pin2,p1,p2,qph(1:ng,:),summt2a)
c      summt2a = 0.5d0*summt2a
c      print*,summt2,summt2a,summt2a/summt2-1.d0
c      stop
            
           do i=0,3
               q(i) = qph(1,i)
           enddo
           if(fs.eq.'er') then

             call realhard_er(pin1,pin2,p1,p2,q,summt2,qph)
             fs = fsbck
             return
           endif
           if(fs.eq.'mr') then

              if(ifisrfsr.eq.11) then
                 QRmu = 1.d0
                 QRe  = 1.d0
c                 fs = fsbck
c                 QRmu = 0.d0
c              call antenna_approx_gm(ng,pin1,pin2,p1,p2,qph(1:2,:),hard)                 
               xxx = 0.5d0*muemuegg(-p2,pin2,p1,-pin1,qph(1,:),qph(2,:)) ! is this 1/2 for photons (1/2!
               hard = xxx
c           call glueALPHA(ng,pin1,pin2,p1,p2,qph(1:ng,:),xxx)
c            xxx = xxx*0.5d0
c            call realhard_mur(pin1,pin2,p1,p2,q,yyy,qph,'full','muc')
c            print*,hard,xxx,xxx/hard-1.d0
c               stop


c            call realhard_mur(pin1,pin2,p1,p2,q,hard,qph,'full','muc')
          elseif(ifisrfsr.eq.10) then

                 QRmu = 0.d0
                 QRe  = 1.d0
            xxx = 0.5d0*muemuegg(-p2,pin2,p1,-pin1,qph(1,:),qph(2,:)) ! is this 1/2 for photons (1/2!
            hard = xxx
c             call realhard_mur(pin1,pin2,p1,p2,q,hard,qph,'isr ','muc')
          elseif(ifisrfsr.eq.01) then
                 QRmu = 1.d0
                 QRe  = 0.d0
            xxx = 0.5d0*muemuegg(-p2,pin2,p1,-pin1,qph(1,:),qph(2,:)) ! is this 1/2 for photons (1/2!
            hard = xxx
c            call realhard_mur(pin1,pin2,p1,p2,q,hard,qph,'fsr ','muc')
          endif
              summt2=hard

              fs = fsbck

              return
           endif
! AGGIUNTO PER RADIATIVE PIONS   !!!!E ANCHE PER I PIONI NORMALI (EFFETTI CHE VEDEVA FRANCESCO CON F_H)
           if(fs.eq.'pr'.or.fs.eq.'pp') then
              iarunbck=iarun
              if(fs.eq.'pp') iarun=0

              call realhard_pr(pin1,pin2,p1,p2,q,hard,qph)

              if(fs.eq.'pp') then
                 iarun=iarunbck
                 hard=hard*bornmepp/bornmenovppp
              endif
              summt2=hard

               fs = fsbck

              return
           endif

        endif

         eikng = 'multip'
         call eikonal_factor(eikng,ng,pin1,pin2,p1,p2,q,qph,eikfull)
         ei1 = pin1(0)
         ei2 = pin2(0)
         ef1 = ei1
         ef2 = ei2
         ef3 = p1(0)
         ef4 = p2(0)
         prodomega = 1.d0
         do k = 1,ng
            prodomega = prodomega*qph(k,0)
            if (iclose(k).eq.3) ef3 = ef3 + qph(k,0)
            if (iclose(k).eq.4) ef4 = ef4 + qph(k,0)
         enddo
         ei3 = ef3
         ei4 = ef4

         prodapx = 1.d0
         xjac    = 1.d0
         delta   = 0.d0
         prod    = 1.d0
         delta_eik = 0.d0
         prod_eik  = 1.d0

	if (fs.ne.'mr'.and.fs.ne.'er'.and.fs.ne.'pr'.and.fs.ne.'pp') then
           iarun=0
           do k = 1,ng

            do i=0,3
               q(i) = qph(k,i)
            enddo
**** PS QUANTITIES
            if (iclose(k).eq.1) then
               xjac = xjac/ef1/2.d0
               xap  = qph(k,0)/ef1
               apx  = (1.d0 + (1.d0 - xap)**2)/xap
               apx  = apx/(1.d0 - xap) !!
               ef1  = ef1 - qph(k,0)
            endif
            if (iclose(k).eq.2) then
               xjac = xjac/ef2/2.d0
               xap  = qph(k,0)/ef2
               apx  = (1.d0 + (1.d0 - xap)**2)/xap
               apx  = apx/(1.d0 - xap) !!
               ef2  = ef2 - qph(k,0)
            endif
            if (iclose(k).eq.3) then
               xjac = xjac/ef3/2.d0
               xap  = qph(k,0)/ef3
               apx  = (1.d0 + (1.d0 - xap)**2)/xap
               if (fs.eq.'pp') apx  = 2.d0*(1.d0 - xap)/xap
               apx  = apx/(1.d0 - xap) !!
               ef3  = ef3 - qph(k,0)
            endif
            if (iclose(k).eq.4) then
               xjac = xjac/ef4/2.d0
               xap  = qph(k,0)/ef4
               apx  = (1.d0 + (1.d0 - xap)**2)/xap
               if (fs.eq.'pp') apx  = 2.d0*(1.d0 - xap)/xap
               apx  = apx/(1.d0 - xap) !!
               ef4    = ef4 - qph(k,0)
            endif



            prodapx = prodapx * apx
! killing the event!!
            if (prodapx.lt.0.d0) then
               ie = 1
               summt2 = 0.d0
               iarun = iarunbck
               fs = fsbck
               return
            endif
***
********CARLO MAP
c            print*, prodapx

            is = 0
            if (iclose(k).eq.1.or.iclose(k).eq.2) is = 1

** doing as if only the k-th photon is present
            i1 = 0
            i2 = 0
            i3 = 0
            i4 = 0
            if (iclose(k).eq.1) i1 = 1
            if (iclose(k).eq.2) i2 = 1
            if (iclose(k).eq.3) i3 = 1
            if (iclose(k).eq.4) i4 = 1
            do ki = 0,3
               p1tmp(ki) = pin1(ki) - i1 * q(ki)
               p2tmp(ki) = pin2(ki) - i2 * q(ki)
               qq(ki) = p1tmp(ki) + p2tmp(ki)
            enddo

            call maptozero(p1tmp,p2tmp,p1,p2,pin1bl,pin2bl,
     .           p1bl,p2bl,iemap)
***
! killing the event!!
            if (iemap.eq.1) then
               ie = 1
               summt2 = 0.d0
               iarun = iarunbck
               fs = fsbck
               return
            endif
***

***
            call maptoone(is,qq,p1,p2,q,pin1r,pin2r,p1r,p2r,qr)

            iq2 = 0
            call born_me(pin1bl,pin2bl,p1bl,p2bl,bornmel,bornmelq2,
     .           bornmelz,iq2,q)


            call eikonal_factor('single',ng,pin1r,pin2r,p1r,p2r,qr,
     >           qph,eiki)
            call bk_me(pin1r,pin2r,p1r,p2r,qr,oal,qph)


            ei1r = pin1r(0)
            ei2r = pin2r(0)
            ef3r = p1r(0)
            ef4r = p2r(0)

            if (iclose(k).eq.1) then
               ef1r = ei1r - qr(0)
               xjaci = 1.d0/ei1r/2.d0
               xapi = qr(0)/ei1r
               apxi = (1.d0 + (1.d0 - xapi)**2)/xapi
               apxi = apxi/(1.d0 - xapi) !!
            elseif (iclose(k).eq.2) then
               ef2r = ei2r - qr(0)
               xjaci = 1.d0/ei2r/2.d0
               xapi = qr(0)/ei2r
               apxi = (1.d0 + (1.d0 - xapi)**2)/xapi
               apxi = apxi/(1.d0 - xapi) !!
            elseif (iclose(k).eq.3) then
               ei3r = ef3r + qr(0)
               xjaci = 1.d0/ei3r/2.d0
               xapi = qr(0)/ei3r
               apxi = (1.d0 + (1.d0 - xapi)**2)/xapi
               if (fs.eq.'pp') apxi  = 2.d0*(1.d0 - xapi)/xapi
               apxi = apxi/(1.d0 - xapi) !!
            elseif (iclose(k).eq.4) then
               ei4r = ef4r + qr(0)
               xjaci = 1.d0/ei4r/2.d0
               xapi = qr(0)/ei4r
               apxi = (1.d0 + (1.d0 - xapi)**2)/xapi
               if (fs.eq.'pp') apxi  = 2.d0*(1.d0 - xapi)/xapi
               apxi = apxi/(1.d0 - xapi) !!
            endif

            approxi =  apxi*eiki*qr(0)*xjaci
!! TO AVOID ROUND-OFFs
            if (approxi.gt.0.d0) then
               rescaletomatcheikfull=eiki/singleeikfact(k)
               deltai = (oal - approxi*bornmel)/approxi
*##
*##            deltai = deltai * rescaletomatcheikfull
*##
            else
               deltai = 0.d0
            endif

            delta  = delta + deltai

            prod = prod * (1.d0 + deltai/bornmel)
!            print*,deltai/bornmel
            deltah = deltah + deltai/bornmel
         enddo

            aaa = prodapx * prodomega
            summt2 = bornme*aaa*xjac*eikfull
            summt2 = summt2 * prod
            imtx = 1
            prodotto = prod

            iarun = iarunbck
            fs = fsbck
            return
         endif

*******RADIATIVE CHANNELS FOR PS

         if(fs.eq.'er'.or.fs.eq.'mr'.or.fs.eq.'pr'.or.fs.eq.'pp') then

            if(fs.eq.'pp') then
               iarun=0
            endif
*=========================================================
**********************MOUNT ON 2->4***********************
c$$$*=========================================================
c$$$

            call gettwophot(pin1,pin2,p1,p2,qph,ng,l1,l2)

            ei1 = pin1(0)
            ei2 = pin2(0)
            ef1 = ei1
            ef2 = ei2
            ef3 = p1(0)
            ef4 = p2(0)

            prodomega = 1.d0
            do k = 1,ng
               prodomega = prodomega*qph(k,0)
            enddo

***************ADDED TO TAKE INTO ACCOUNT THE FACT THAT THE HARDEST PHOTONS CONTRIBUTE TO THE EF_I, BUT NOT TO PRODOMEGA
            prodomega=prodomega/qph(l1,0)
            prodomega=prodomega/qph(l2,0)
**********************************

         ei3 = ef3
         ei4 = ef4
         prodapx = 1.d0
         xjac    = 1.d0
         eikfull = 1.d0
         pin1cl=pin1
         pin2cl=pin2

         p1cl=p1
         p2cl=p2

        do k = 1,ng
*     * PS QUANTITIES
           if(k.ne.l1.and.k.ne.l2) then
            if (iclose(k).eq.1) then
               xjac = xjac/ef1/2.d0
               xap  = qph(k,0)/ef1
               apx  = (1.d0 + (1.d0 - xap)**2)/xap
               apx  = apx/(1.d0 - xap) !!
               ef1  = ef1 - qph(k,0)
               pin1cl=pin1cl-qph(k,:)
            endif
            if (iclose(k).eq.2) then
               xjac = xjac/ef2/2.d0
               xap  = qph(k,0)/ef2
               apx  = (1.d0 + (1.d0 - xap)**2)/xap
               apx  = apx/(1.d0 - xap) !!
               ef2  = ef2 - qph(k,0)
               pin2cl=pin2cl-qph(k,:)
            endif

            if (iclose(k).eq.3) then
               ef3  = ef3 + qph(k,0)
               xjac = xjac/ef3/2.d0
               xap  = qph(k,0)/ef3
               apx  = (1.d0 + (1.d0 - xap)**2)/xap
               if (fs.eq.'pr'.or.fs.eq.'pp') apx = 2.d0*(1.d0 - xap)/xap
               apx  = apx/(1.d0 - xap) !!
               p1cl=p1cl+qph(k,:)
            endif
            if (iclose(k).eq.4) then
               ef4  = ef4 + qph(k,0)
               xjac = xjac/ef4/2.d0
               xap  = qph(k,0)/ef4
               apx  = (1.d0 + (1.d0 - xap)**2)/xap
               if (fs.eq.'pr'.or.fs.eq.'pp') apx = 2.d0*(1.d0 - xap)/xap
               apx  = apx/(1.d0 - xap) !!
               p2cl=p2cl+qph(k,:)
            endif
            prodapx = prodapx * apx

            call eikonal_factor_z('single',1,pin1,pin2,p1,p2,qph(k,:),
     >        xap,qph,eiksingle)

            eikfull = eikfull*eiksingle

            if (prodapx.lt.0.d0) then
               ie = 1
               summt2 = 0.d0
               iarun = iarunbck
               fs = fsbck
               return
            endif
           endif
         enddo
*
*******I DIVIDE THE EIKFULL BY (NG-1)!
        eikfull=eikfull/factorial(ng-1)

        call maptotwo_rad(pin1cl,pin2cl,p1cl,p2cl,qph(l1,:),
     >       qph(l2,:),pin1bt,pin2bt,p1bt,p2bt,qbt,q1bt)
        qph2(1,:)=qbt
        qph2(2,:)=q1bt
***************************
       call check_propagators(pin1,pin2,qph,ng,pin1bt,pin2bt,qph2,error)

       if(error.eq.1) then
          counterg4=counterg4+1
          summt2=0.d0
          fs = fsbck
          return
       endif
**********************************

        if(fs.eq.'pr'.or.fs.eq.'pp') then
           call realhard_pr(pin1bt,pin2bt,p1bt,p2bt,
     >                            qph2(1,:),bornme,qph2)
        elseif(fs.eq.'mr') then
           if(ifisrfsr.eq.11) then
              QRmu = 1.d0
              QRe  = 1.d0
              bornme = 0.5d0*
     >         muemuegg(-p2bt,pin2bt,p1bt,-pin1bt,qph2(1,:),qph2(2,:))

c              call realhard_mur(pin1bt,pin2bt,p1bt,p2bt
c     >             ,qph2(1,:),bornme,qph2,'full','muc')
           elseif(ifisrfsr.eq.10) then
              QRmu = 0.d0
              QRe  = 1.d0
              bornme = 0.5d0*
     >         muemuegg(-p2bt,pin2bt,p1bt,-pin1bt,qph2(1,:),qph2(2,:))

c              call realhard_mur(pin1bt,pin2bt,p1bt,p2bt
c     >             ,qph2(1,:),bornme,qph2,'isr ','muc')
           elseif(ifisrfsr.eq.01) then
              QRmu = 1.d0
              QRe  = 0.d0
              bornme = 0.5d0*
     >         muemuegg(-p2bt,pin2bt,p1bt,-pin1bt,qph2(1,:),qph2(2,:))

c              call realhard_mur(pin1bt,pin2bt,p1bt,p2bt
c     >             ,qph2(1,:),bornme,qph2,'fsr ','muc')
           endif
        elseif(fs.eq.'er') then
           call realhard_er(pin1bt,pin2bt,p1bt,p2bt,qph2(1,:),
     >          bornme,qph2)
        endif
        aaa = prodapx * prodomega


        summt2 = bornme * aaa * xjac * eikfull*2.d0

        summt2 = summt2 * extraden

c$$$        if (ng.gt.2.and.fs.eq.'mr') then
c$$$           fs='mm'
c$$$!     print*,ng
c$$$           call glueALPHA(ng,pin1,pin2,p1,p2,qph(1:ng,:),xxx)
c$$$           
c$$$           summt2 = summt2-xxx/factorial(ng)*ng
c$$$           
c$$$           imtx = 1
c$$$           iarun = iarunbck
c$$$           fs='mr'
c$$$!     return
c$$$        endif

        if(fs.eq.'pp') then
           iarun=iarunbck
           summt2=summt2*bornmepp/bornmenovppp
        endif
*===== ====================================================
*=========================================================
        fs = fsbck
      endif
         fs = fsbck
      endif

*******************
**     PSRAD     **
*******************
c$$$      if(model.eq.'psrad') then
c$$$
c$$$*=========================================================
c$$$**********************MOUNT ON 2->4***********************
c$$$c$$$*=========================================================
c$$$c$$$
c$$$         if(ng.eq.1) then
c$$$            do i=0,3
c$$$               q(i)=qph(1,i)
c$$$            enddo
c$$$            if(fs.eq.'mr'.or.fs.eq.'er'
c$$$     >           .or.fs.eq.'pr') then
c$$$               pin1br=pin1
c$$$               pin2br=pin2
c$$$               p1br=p1
c$$$               p2br=p2
c$$$               qbr=q
c$$$               call  bk_me(pin1,pin2,p1,p2,q,summt2,qph)
c$$$               summt2=summt2
c$$$            imtx = 1
c$$$            return
c$$$         endif
c$$$    	 endif
c$$$
c$$$
c$$$         call gettwophot(pin1,pin2,p1,p2,qph,ng,l1,l2)
c$$$
c$$$         ei1 = pin1(0)
c$$$         ei2 = pin2(0)
c$$$         ef1 = ei1
c$$$         ef2 = ei2
c$$$         ef3 = p1(0)
c$$$         ef4 = p2(0)
c$$$
c$$$
c$$$
c$$$         prodomega = 1.d0
c$$$         do k = 1,ng
c$$$            prodomega = prodomega*qph(k,0)
c$$$         enddo
c$$$
c$$$***************ADDED TO TAKE INTO ACCOUNT THE FACT THAT THE HARDEST PHOTONS CONTRIBUTE TO THE EF_I, BUT NOT TO PRODOMEGA
c$$$         prodomega=prodomega/qph(l1,0)
c$$$         prodomega=prodomega/qph(l2,0)
c$$$*****************************
c$$$
c$$$         ei3 = ef3
c$$$         ei4 = ef4
c$$$         prodapx = 1.d0
c$$$         xjac    = 1.d0
c$$$         eikfull = 1.d0
c$$$	 pin1cl=pin1
c$$$	 pin2cl=pin2
c$$$
c$$$	 p1cl=p1
c$$$	 p2cl=p2
c$$$
c$$$        do k = 1,ng
c$$$** PS QUANTITIES
c$$$            if(k.ne.l1.and.k.ne.l2) then
c$$$            if (iclose(k).eq.1) then
c$$$               xjac = xjac/ef1/2.d0
c$$$               xap  = qph(k,0)/ef1
c$$$               apx  = (1.d0 + (1.d0 - xap)**2)/xap
c$$$               apx  = apx/(1.d0 - xap) !!
c$$$               ef1  = ef1 - qph(k,0)
c$$$               pin1cl=pin1cl-qph(k,:)
c$$$            endif
c$$$            if (iclose(k).eq.2) then
c$$$               xjac = xjac/ef2/2.d0
c$$$               xap  = qph(k,0)/ef2
c$$$               apx  = (1.d0 + (1.d0 - xap)**2)/xap
c$$$               apx  = apx/(1.d0 - xap) !!
c$$$               ef2  = ef2 - qph(k,0)
c$$$               pin2cl=pin2cl-qph(k,:)
c$$$            endif
c$$$
c$$$            if (iclose(k).eq.3) then
c$$$               ef3  = ef3 + qph(k,0)
c$$$               xjac = xjac/ef3/2.d0
c$$$               xap  = qph(k,0)/ef3
c$$$               apx  = (1.d0 + (1.d0 - xap)**2)/xap
c$$$               if (fs.eq.'pr') apx = 2.d0*(1.d0 - xap)/xap
c$$$               apx  = apx/(1.d0 - xap) !!
c$$$               p1cl=p1cl+qph(k,:)
c$$$            endif
c$$$            if (iclose(k).eq.4) then
c$$$               ef4  = ef4 + qph(k,0)
c$$$               xjac = xjac/ef4/2.d0
c$$$               xap  = qph(k,0)/ef4
c$$$               apx  = (1.d0 + (1.d0 - xap)**2)/xap
c$$$               if (fs.eq.'pr') apx = 2.d0*(1.d0 - xap)/xap
c$$$               apx  = apx/(1.d0 - xap) !!
c$$$               p2cl=p2cl+qph(k,:)
c$$$            endif
c$$$            prodapx = prodapx * apx
c$$$
c$$$               call eikonal_factor_z('single',1,pin1,pin2,p1,p2,qph(k,:)
c$$$     >              ,xap,qph,eiksingle)
c$$$
c$$$     	    eikfull = eikfull*eiksingle
c$$$
c$$$            if (prodapx.lt.0.d0) then
c$$$               ie = 1
c$$$               summt2 = 0.d0
c$$$               iarun = iarunbck
c$$$               return
c$$$            endif
c$$$           endif
c$$$         enddo
c$$$*
c$$$
c$$$*******I DIVIDE THE EIKFULL BY (NG-1)!
c$$$	eikfull=eikfull/factorial(ng-1)
c$$$
c$$$         call maptotwo_rad(pin1cl,pin2cl,p1cl,p2cl,qph(l1,:),
c$$$     >        qph(l2,:),pin1bt,pin2bt,p1bt,p2bt,qbt,q1bt)
c$$$         qph2(1,:)=qbt
c$$$         qph2(2,:)=q1bt
c$$$
c$$$
c$$$***************************
c$$$         call check_propagators(pin1,pin2,qph,ng,pin1bt,pin2bt,qph2,
c$$$     >        error)
c$$$
c$$$
c$$$         if(error.eq.1) then
c$$$            counterg4=counterg4+1
c$$$            summt2=0.d0
c$$$            return
c$$$         endif
c$$$
c$$$**********************************
c$$$
c$$$         if(fs.eq.'pr') then
c$$$            call realhard_pr(pin1bt,pin2bt,p1bt,p2bt,
c$$$     >           qph2(1,:),bornme,qph2)
c$$$         elseif(fs.eq.'mr') then
c$$$            if(ifisrfsr.eq.11) then
c$$$               call realhard_mur(pin1bt,pin2bt,p1bt,p2bt
c$$$     >              ,qph2(1,:),bornme,qph2,'full','muc')
c$$$            elseif(ifisrfsr.eq.10) then
c$$$               call realhard_mur(pin1bt,pin2bt,p1bt,p2bt
c$$$     >              ,qph2(1,:),bornme,qph2,'isr ','muc')
c$$$            elseif(ifisrfsr.eq.01) then
c$$$               call realhard_mur(pin1bt,pin2bt,p1bt,p2bt
c$$$     >              ,qph2(1,:),bornme,qph2,'fsr ','muc')
c$$$            endif
c$$$         elseif(fs.eq.'er') then
c$$$            call realhard_er(pin1bt,pin2bt,p1bt,p2bt,qph2(1,:),
c$$$     >           bornme,qph2)
c$$$         endif
c$$$         aaa = prodapx * prodomega
c$$$
c$$$         summt2 = bornme * aaa * xjac * eikfull*2.d0
c$$$         if(ng.eq.2) then
c$$$            summt2=summt2/2.d0
c$$$         endif
c$$$
c$$$      endif
c$$$
c$$$*******************
c$$$c$$$**      PS       **
c$$$c$$$*******************
c$$$      if (model.eq.'ps') then
c$$$         if (ng.eq.0) then
c$$$            iref = 1
c$$$            summt2 = bornme
c$$$         imtx = 1
c$$$         iarun = iarunbck
c$$$         return
c$$$         endif
c$$$*
c$$$         ei1 = pin1(0)
c$$$         ei2 = pin2(0)
c$$$         ef1 = ei1
c$$$         ef2 = ei2
c$$$         ef3 = p1(0)
c$$$         ef4 = p2(0)
c$$$
c$$$         prodomega = 1.d0
c$$$         do k = 1,ng
c$$$            prodomega = prodomega*qph(k,0)
c$$$            if (iclose(k).eq.3) ef3 = ef3 + qph(k,0)
c$$$            if (iclose(k).eq.4) ef4 = ef4 + qph(k,0)
c$$$         enddo
c$$$
c$$$         ei3 = ef3
c$$$         ei4 = ef4
c$$$
c$$$         prodapx = 1.d0
c$$$         xjac    = 1.d0
c$$$         do k = 1,ng
c$$$** PS QUANTITIES
c$$$
c$$$            if (iclose(k).eq.1) then
c$$$               xjac = xjac/ef1/2.d0
c$$$               xap  = qph(k,0)/ef1
c$$$               apx  = (1.d0 + (1.d0 - xap)**2)/xap
c$$$               apx  = apx/(1.d0 - xap) !!
c$$$               ef1  = ef1 - qph(k,0)
c$$$            endif
c$$$            if (iclose(k).eq.2) then
c$$$               xjac = xjac/ef2/2.d0
c$$$               xap  = qph(k,0)/ef2
c$$$               apx  = (1.d0 + (1.d0 - xap)**2)/xap
c$$$               apx  = apx/(1.d0 - xap) !!
c$$$               ef2  = ef2 - qph(k,0)
c$$$            endif
c$$$            if (iclose(k).eq.3) then
c$$$               xjac = xjac/ef3/2.d0
c$$$               xap  = qph(k,0)/ef3
c$$$               apx  = (1.d0 + (1.d0 - xap)**2)/xap
c$$$               if (fs.eq.'pp') apx = 2.d0*(1.d0 - xap)/xap
c$$$               apx  = apx/(1.d0 - xap) !!
c$$$               ef3  = ef3 - qph(k,0)
c$$$            endif
c$$$            if (iclose(k).eq.4) then
c$$$               xjac = xjac/ef4/2.d0
c$$$               xap  = qph(k,0)/ef4
c$$$               apx  = (1.d0 + (1.d0 - xap)**2)/xap
c$$$               if (fs.eq.'pp') apx = 2.d0*(1.d0 - xap)/xap
c$$$               apx  = apx/(1.d0 - xap) !!
c$$$               ef4  = ef4 - qph(k,0)
c$$$            endif
c$$$            prodapx = prodapx * apx
c$$$
c$$$            if (prodapx.lt.0.d0) then
c$$$               ie = 1
c$$$               summt2 = 0.d0
c$$$               iarun = iarunbck
c$$$               return
c$$$            endif
c$$$****
c$$$!	     write(16,*)'iclose(k)',iclose(k)
c$$$!             write(16,*)'nphot',k
c$$$!             write(16,*)'jac',xjac
c$$$!             write(16,*)'xap',xap
c$$$!             write(16,*)'apx',apx
c$$$         enddo
c$$$
c$$$         if (ng.gt.1) then
c$$$            eikng = 'multip'
c$$$            call eikonal_factor(eikng,ng,pin1,pin2,p1,p2,q,qph,eikfull)
c$$$         else
c$$$            do i=0,3
c$$$               q(i) = qph(1,i)
c$$$            enddo
c$$$            eikng = 'single'
c$$$            call eikonal_factor(eikng,ng,pin1,pin2,p1,p2,q,qph,eikfull)
c$$$         endif
c$$$
c$$$         aaa = prodapx * prodomega
c$$$
c$$$c         summt2 = bornme * prodapx * eikfull * prodomega * xjac
c$$$         summt2 = bornme * aaa * eikfull * xjac
c$$$
c$$$c         print*,summt2,bornme,aaa,eikfull,xjac
c$$$!         write(16,*)'final emtx',summt2
c$$$!         write(16,*)'born',bornme
c$$$!         write(16,*)'eikfull', eikfull
c$$$         imtx = 1
c$$$         iarun = iarunbck
c$$$         return
c$$$      endif
c$$$******
c$$$      if (model.eq.'eikonal') then
c$$$         eikng = 'multip'
c$$$         call eikonal_factor(eikng,ng,pin1,pin2,p1,p2,q,qph,eikonal)
c$$$         summt2 = eikonal * bornme
c$$$         imtx = 1
c$$$         iarun = iarunbck
c$$$         return
c$$$      endif

******
!! ALPHA model not released!
c$$$      if (model.eq.'alpha') then
c$$$! per ng = 0,1 e' inutile chiamare alpha...
c$$$         if (ng.eq.0) then
c$$$            iref = 1
c$$$            summt2 = bornme
c$$$            imtx = 1
c$$$            return
c$$$         endif
c$$$         if (ng.eq.1) then
c$$$            do i=0,3
c$$$               q(i) = qph(1,i)
c$$$            enddo
c$$$c            call alpha_me('S',ng,p1,p2,qph,ialpha,emtxalpha)
c$$$
c$$$
c$$$            call bk_me(pin1,pin2,p1,p2,q,summt2)
c$$$ccc            call approxNg(ng,pin1,pin2,p1,p2,qph,summt2)
c$$$ccc            ialpha = 1
c$$$
c$$$            summt2 = summt2/factorial(ng)
c$$$            imtx = 1
c$$$            return
c$$$         endif
c$$$******
c$$$* for ALPHA, not released
c$$$c         call alpha_me('S',ng,p1,p2,qph,ialpha,emtxalpha) ! summed over spins
c$$$c         call alpha_me('R',ng,p1,p2,qph,ialpha,emtxalpha) ! random over spins
c$$$         emt1 = emtxalpha
c$$$         do k = 0,3
c$$$            ptmp(k)  = qph(1,k)
c$$$            qph(1,k) = qph(2,k)
c$$$            qph(2,k) = ptmp(k)
c$$$         enddo
c$$$         call alpha_me('S',ng,p1,p2,qph,ialpha,emtxalpha) ! random over spins
c$$$         emt2 = emtxalpha
c$$$
c$$$         print*,emt2/emt1
c$$$
c$$$
c$$$c$$$         call approxNg(ng,pin1,pin2,p1,p2,qph,emtxalpha)
c$$$c$$$         emtxalpha = emtxalpha/factorial(ng)
c$$$c$$$         ialpha = 1
c$$$
c$$$         imtx   = ialpha
c$$$         summt2 = emtxalpha
c$$$         return
c$$$      endif
******
      iarun = iarunbck
      return
      end


c     Routine to evaluate the matrix elements related to phi decays,bremstrhalung, double resonant
      subroutine squared_matrix_phi(model,ng,ecms,p1,p2,pbeam1,pbeam2,
     >     qph,ie,imtx,summt2,prodotto,deltah)
      implicit double precision (a-h,o-z)
      character*10 model
      character*6   ord
      character*6 eikng ! multip single
      parameter (imaxph = 40)
      integer nout,nout1
      parameter (nout   =  3)
      parameter (nout1  =  2)
      dimension iclose(40)
      dimension iclose2(2)
      double precision mass,lambda
      dimension p1(0:3),p2(0:3),pin1(0:3),pin2(0:3),qmost(0:3)
      dimension p1r(0:3),p2r(0:3),p1b(0:3),p2b(0:3),q(0:3),qf(0:3)
      dimension q1(0:3),q2(0:3),q3(0:3),q4(0:3),qq(0:3),qr(0:3),qb(0:3)
      dimension p1tmp(0:3),p2tmp(0:3),p3tmp(0:3),p4tmp(0:3)
      dimension pin1bl(0:3),pin2bl(0:3),p1bl(0:3),p2bl(0:3)
      dimension pin1b1(0:3),pin2b1(0:3),p1b1(0:3),p2b1(0:3)
      dimension pin1bt(0:3),pin2bt(0:3),p1bt(0:3),p2bt(0:3),qbt(0:3)
      dimension q1bt(0:3),newmom(0:3)
      common/maptotworad/pin1bt,pin2bt,p1bt,p2bt,qbt,q1bt
      common/counterg5/counterg4
      integer counterg4
c$$$      common/check_4g/matrix_check,eik_jac,summt2ex
c$$$      real*8 matrix_check
      double precision newmom
      dimension pin1cl(0:3),pin2cl(0:3),p1cl(0:3),p2cl(0:3),qcl(4)
      dimension pin1cl1(4),pin2cl1(4),p1cl1(4),p2cl1(4),qclu(4,2)
      dimension pphot(4,2),idphtme(2),pclout1(4,2)
      dimension pclb(4,0:3),pcla(4,0:3)
      dimension ptmp(0:3),ptmp2(0:3), q2r(0:3),q1r(0:3),pmapout(4,3)
      dimension pmapout1(4,2),qph2(imaxph,0:3)
      dimension pclin(4,2),pclout(4,3),pmapin(4,2),pchrg(4,4)
      double precision massesin(2),massesout(3),mdot,massesout1(2)
      external mdot
      external dot
      double precision dot
      dimension pbeam1(0:3),pbeam2(0:3)
      dimension qmostpt(imaxph,0:3)
      real*8 allun,qfixrad(0:3),qqtwo(0:3)
      real*8 qphr(imaxph,0:3), qtmp(imaxph,0:3)
      dimension p1rr(0:3),p2rr(0:3),pin1rr(0:3),pin2rr(0:3),prif(0:3)
      dimension pin1br(0:3),pin2br(0:3),p1br(0:3),p2br(0:3),qbr(0:3)
      dimension pin1brr(0:3),pin2brr(0:3),p1brr(0:3),p2brr(0:3)
      dimension qbrr(0:3)
      dimension pin1tm(0:3),pin2tm(0:3),p1tm(0:3),p2tm(0:3)
      dimension qph(imaxph,0:3),qphtmp(imaxph,0:3),randvec(3)
      dimension qphtest(imaxph,0:3)
      real*8 hard
      logical condrad
      integer l,phi
      character*2 fs,fsbck
      real*8 tmphardest(0:3)
      common/parameters/ame,ammu,convfac,alpha,pi
      common/qedORDER/ord
      common/finalstate/fs
      integer iffpi,iffpibck,vals,u,error
      dimension vals(2)
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi
      common/iclosest/iclose
      common/momentainitial/pin1,pin2
      common/reducedtoborn/p1b,p2b,iref
      common/momentainitialred/pin1b(0:3),pin2b(0:3),pin1r(0:3),
     >     pin2r(0:3)
      common/forborncrosssection/phsp2b,flux2b,bornme,bornmeq2,bornmez,
     >     bornmenovp
      common/bornxsecpp/bornmenovppp,bornmepp
      common/radpattern/nph(4)
      common/hparam/ampi
      common/epssoft/eps
      common/fortests/prod
      common/for_debug/ppp(0:3),denphsp,dir3(0:3),dir4(0:3),aabr,aabq,
     > aaonm
      common/singlefactors/singleeikfact(imaxph)
      common/idarkon/idarkon
      common/ialpharunning/iarun
      double complex vpolc
      external vpolc

!     imtx is to say that the matrix el. has been calculated..
! summt2 is the squared ME summed over all all spins

      common/radoption/ifisrfsr
      integer ifisrfsr

      integer iteubn
      common/teubner/iteubn
      common/scaling/scaling
*      for recola ee -> 4g
      dimension p6(0:3),p7(0:3)
      common/p6for4g/p6,p7

***   to check weird events I need some commons
*
c      real*8 deltah,deltai
!      common/momentacheck1/pin1br,pin2br,p1br,p2br,qbr
      common/fh/fhmax
      double precision fhmax

** VP UNC
      dimension iswvperr(3),iswvperrt(6)
      common/iswitchvperror/iswvperr,iswvperrt
      parameter (nvps = 26)
      dimension remtx(nvps)
      common/rationsvpunc/remtx,invp
      common/rescaledhaderr/scal_err_dhad

      parameter (nvpst = 728)
      dimension remtxt(nvpst)
      integer ist(nvpst,6)
      common/rationsvpunct/remtxt,invpt,ist
***
      common/firstsqmatrix/ifirst
      data ifirst /1/

      double precision fsrphi
      double complex ifiphi
      
      if(fs.ne.'pr'.and.fs.ne.'pp') then
         print*,'THIS SHOULD NOT HAPPEN'
         stop
      endif

      if(ng.eq.1)
     >     then
         do k = 0,3
            pin1br(k) = pbeam1(k)
            pin2br(k) = pbeam2(k)
            p1br(k) = p1(k)
            p2br(k) = p2(k)
            qbr(k)  = qph(1,k)
            enddo
         call FSR_pr_phi(pin1br,pin2br,p1br,p2br,qbr,fsrphi,qph)
         call IFI_pr_phi(pin1br,pin2br,p1br,p2br,qbr,ifiphi,qph)

         summt2 = real(fsrphi+2.d0*ifiphi)
!         print*,'Fxsqed',summt2
         call FSR_pr(pin1br,pin2br,p1br,p2br,qbr,fsrphi,qph)
         call IFI_pr(pin1br,pin2br,p1br,p2br,qbr,ifiphi,qph)

         summt2 = summt2 - real(fsrphi+2.d0*ifiphi)
!         print*,'full',summt2
      else
**********************MAPPING MAURO MORETTI
        
         call closest(ng,qph,pin1,pin2,p1,p2,q1,q2,q3,q4,iclose)

         ind=3
         call getsignalphot(pin1,pin2,p1,p2,qph,ng,ind,j)
********FROM NOW ON THE SIGNAL PHOTON IS CALLED J

         pin1cl=pin1
         pin2cl=pin2

         p1cl=p1
         p2cl=p2

         do k = 1,ng
*     * PS QUANTITIES
            if(k.ne.j) then
!               print*,iclose(k)
               if (iclose(k).eq.1) then
                  pin1cl=pin1cl-qph(k,:)
               endif
               if (iclose(k).eq.2) then
                  pin2cl=pin2cl-qph(k,:)
               endif
               if (iclose(k).eq.3) then !
                  p1cl=p1cl+qph(k,:)
               endif
               if (iclose(k).eq.4) then
                  p2cl=p2cl+qph(k,:)
            endif
    	 endif

         enddo

         ind=1
         call maptoone_rad(pin1cl,pin2cl,p1cl,p2cl,qph(j,:),
     >        iclose(j),ind,pin1br,pin2br,p1br,p2br,qbr)


!     Built on top of the 1gamma kinematics
         ei1 = pin1(0)
         ei2 = pin2(0)
         ef1 = ei1
         ef2 = ei2
         ef3 = p1(0)
         ef4 = p2(0)


         prodomega = 1.d0
         do k = 1,ng
            prodomega = prodomega*qph(k,0)
         enddo

***************ADDED TO TAKE INTO ACCOUNT THE FACT THAT THE HARDEST PHOTONS CONTRIBUTE TO THE EF_I, BUT NOT TO PRODOMEGA
         prodomega=prodomega/qph(j,0)
*******

         ei3 = ef3
         ei4 = ef4
         prodapx = 1.d0
         xjac    = 1.d0
         eikfull = 1.d0
         pin1cl=pin1
         pin2cl=pin2

         p1cl=p1
         p2cl=p2

         do k = 1,ng
*     * PS QUANTITIES
           if(k.ne.j) then
               if (iclose(k).eq.1) then
                  xjac = xjac/ef1/2.d0
                  xap  = qph(k,0)/ef1
                  apx  = (1.d0 + (1.d0 - xap)**2)/xap
                  apx  = apx/(1.d0 - xap) !!
                  ef1  = ef1 - qph(k,0)
                  pin1cl=pin1cl-qph(k,:)
               endif
               if (iclose(k).eq.2) then
                  xjac = xjac/ef2/2.d0
                  xap  = qph(k,0)/ef2
                  apx  = (1.d0 + (1.d0 - xap)**2)/xap
                  apx  = apx/(1.d0 - xap) !!
                  ef2  = ef2 - qph(k,0)
                  pin2cl=pin2cl-qph(k,:)
               endif

               if (iclose(k).eq.3) then
                  ef3  = ef3 + qph(k,0)
                  xjac = xjac/ef3/2.d0
                  xap  = qph(k,0)/ef3
                  apx  = (1.d0 + (1.d0 - xap)**2)/xap
               if (fs.eq.'pr'.or.fs.eq.'pp') apx = 2.d0*(1.d0 - xap)/xap
                  apx  = apx/(1.d0 - xap) !!
                  p1cl=p1cl+qph(k,:)
               endif
               if (iclose(k).eq.4) then
                  ef4  = ef4 + qph(k,0)
                  xjac = xjac/ef4/2.d0
                  xap  = qph(k,0)/ef4
                  apx  = (1.d0 + (1.d0 - xap)**2)/xap
               if (fs.eq.'pr'.or.fs.eq.'pp') apx = 2.d0*(1.d0 - xap)/xap
                  apx  = apx/(1.d0 - xap) !!
                  p2cl=p2cl+qph(k,:)
               endif
               prodapx = prodapx * apx

            call eikonal_factor_z('single',1,pin1,pin2,p1,p2,qph(k,:),
     >        xap,qph,eiksingle)

               eikfull = eikfull*eiksingle

               if (prodapx.lt.0.d0) then
                  ie = 1
                  summt2 = 0.d0
                  iarun = iarunbck
                  return
               endif
            endif
         enddo
*

*******I DIVIDE THE EIKFULL BY (NG)!
        eikfull=eikfull/factorial(ng)


        call FSR_pr_phi(pin1br,pin2br,p1br,p2br,qbr,fsrphi,qph)
        call IFI_pr_phi(pin1br,pin2br,p1br,p2br,qbr,ifiphi,qph)

        summt2 = real(fsrphi+2.d0*ifiphi)

        call FSR_pr(pin1br,pin2br,p1br,p2br,qbr,fsrphi,qph)
        call IFI_pr(pin1br,pin2br,p1br,p2br,qbr,ifiphi,qph)

        summt2 = summt2 - real(fsrphi+2.d0*ifiphi)


         aaa = prodapx * prodomega
        summt2 = summt2 * aaa * xjac * eikfull


         endif


               return
      end


****************************************************************************
* e+e- --> e+e-g
****************************************************************************
      subroutine formme(p1,p2,p3,p4,q,elmat2)
*      implicit none
      implicit double precision (Z)
      double precision elmat2AA
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision elmat2,ame,convfac,alpha,pi,e2,qph(40,0:3)
      double precision dot,elmat2cmn,e2approx,eikf
      double precision me,me2,me4,me6,ammu
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      double precision p1p22,p1p32,p1p42,p2p32,p2p42,p3p42
      double precision p1q,p2q,p3q,p4q
      double precision p1qm1,p2qm1,p3qm1,p4qm1
      double precision p1qm2,p2qm2,p3qm2,p4qm2
      double precision s12,s34,t24,t13,u14,u23
      double precision s12m1,s34m1,t24m1,t13m1
      double precision s12m2,s34m2,t24m2,t13m2
      double complex vp12,vp34,vp24,vp13,vp12d,vp34d
      double precision vp122,vp342,vp242,vp132
      double precision ptmp(0:3)
      double precision vpol
      external vpol
      double complex vpolc
      external vpolc
      integer k
      character*2 fs
      common/finalstate/fs
      common/parameters/ame,ammu,convfac,alpha,pi
      common/testsdebug/elmat2cmn,e2approx,s12,s34,t24,t13,u14,u23
      external dot

      e2 = alpha*4.d0*pi

      me  = ame
      me2 = me*me
      me4 = me2*me2
      me6 = me4*me2

      ptmp = p1+p2
      s12 = dot(ptmp,ptmp)
      s12m1 = 1.d0/s12
      s12m2 = s12m1/s12
      vp12  = vpolc(s12)
      vp122 = abs(vp12)**2

      ptmp = p3+p4
      s34 = dot(ptmp,ptmp)
      s34m1 = 1.d0/s34
      s34m2 = s34m1/s34
      vp34 = vpolc(s34)
      vp342 = abs(vp34)**2

      vp12d = dconjg(vp12)
      vp34d = dconjg(vp34)

      ptmp=p2-p4
      t24 = dot(ptmp,ptmp)
      t24m1 = 1.d0/t24
      t24m2 = t24m1/t24
      vp24 = vpolc(t24)
      vp242 = vp24*vp24

      ptmp = p1-p3
      t13= dot(ptmp,ptmp)
      t13m1 = 1.d0/t13
      t13m2 = t13m1/t13
      vp13 = vpolc(t13)
      vp132 = vp13*vp13

      p1p2  = dot(p1,p2)
      p1p22 = p1p2*p1p2

      p1p3  = dot(p1,p3)
      p1p32 = p1p3*p1p3

      p1p4  = dot(p1,p4)
      p1p42 = p1p4*p1p4

      p2p3  = dot(p2,p3)
      p2p32 = p2p3*p2p3

      p2p4  = dot(p2,p4)
      p2p42 = p2p4*p2p4

      p3p4  = dot(p3,p4)
      p3p42 = p3p4*p3p4

      p1q   = dot(p1,q)
      p2q   = dot(p2,q)
      p3q   = dot(p3,q)
      p4q   = dot(p4,q)

      p1qm1 =  1.d0/p1q
      p1qm2 = p1qm1/p1q
      p2qm1 =  1.d0/p2q
      p2qm2 = p2qm1/p2q
      p3qm1 =  1.d0/p3q
      p3qm2 = p3qm1/p3q
      p4qm1 =  1.d0/p4q
      p4qm2 = p4qm1/p4q

c      include 'form/formme_OLD.F'
c      print*,vp12,vp34
c      print*,elmat2
      include 'form/formme.F'
c      print*,elmat2

      elmat2 = elmat2 * e2**3
      return
      end
****************************************************************************
*
* e+e- --> mu+mu-g
****************************************************************************
      subroutine formmemm(p1,p2,p3,p4,q,elmat2)
c      implicit none
      implicit double precision (Z)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision elmat2,ame,convfac,alpha,pi,e2,qph(40,0:3)
      double precision dot,elmat2cmn,e2approx,eikf,elmat2bck
      double precision me,me2,me4,me6,ammu
      double precision mm,mm2,mm4,mm6
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      double precision p1p22,p1p32,p1p42,p2p32,p2p42,p3p42
      double precision p1q,p2q,p3q,p4q
      double precision p1qm1,p2qm1,p3qm1,p4qm1
      double precision p1qm2,p2qm2,p3qm2,p4qm2
      double precision s12,s34,t24,t13,u14,u23
      double precision s12m1,s34m1,t24m1,t13m1
      double precision s12m2,s34m2,t24m2,t13m2
      double precision vp24,vp13
      double precision vp122,vp342,vp242,vp132
      double precision ptmp(0:3)
      double precision vpol
      external vpol
      double complex vpolc,vp12,vp34,vp12d,vp34d,elmat2c
      external vpolc
      integer k
      character*2 fs
      common/finalstate/fs
      common/parameters/ame,ammu,convfac,alpha,pi
      common/testsdebug/elmat2cmn,e2approx,s12,s34,t24,t13,u14,u23
      external dot
      integer ifisrfsr
      common/radoption/ifisrfsr

      e2 = alpha*4.d0*pi

*******
      me  = ame
      me2 = me*me
      me4 = me2*me2
      me6 = me4*me2

      mm  = ammu
      mm2 = mm*mm
      mm4 = mm2*mm2
      mm6 = mm4*mm2

      s12   = dot(p1+p2,p1+p2)
      s12m1 = 1.d0/s12
      s12m2 = s12m1/s12
      vp12  = vpolc(s12)        !* 0.d0

      vp12d = dconjg(vp12)
c      vp122 = abs(vp12)**2 !* 0d0
      vp122 = vp12*vp12d

      s34   = dot(p3+p4,p3+p4)
      s34m1 = 1.d0/s34
      s34m2 = s34m1/s34!!  * 0.d0
      vp34  = vpolc(s34) !* 0d0

      vp34d = dconjg(vp34)
c      vp342 = abs(vp34)**2 !* 0.d0
      vp342 = vp34*vp34d

      if (ifisrfsr.eq.01) then
         vp34  = 0.d0
         vp34d = 0.d0
         vp342 = 0.d0
      elseif (ifisrfsr.eq.10) then
         vp12  = 0.d0
         vp12d = 0.d0
         vp122 = 0.d0
      endif
*
**
*
      p1p2  = dot(p1,p2)
      p1p22 = p1p2*p1p2

      p1p3  = dot(p1,p3)
      p1p32 = p1p3*p1p3

      p1p4  = dot(p1,p4)
      p1p42 = p1p4*p1p4

      p2p3  = dot(p2,p3)
      p2p32 = p2p3*p2p3

      p2p4  = dot(p2,p4)
      p2p42 = p2p4*p2p4

      p3p4  = dot(p3,p4)
      p3p42 = p3p4*p3p4

      p1q   = dot(p1,q)
      p2q   = dot(p2,q)
      p3q   = dot(p3,q)
      p4q   = dot(p4,q)
      p1qm1 =  1.d0/p1q
      p1qm2 = p1qm1/p1q
      p2qm1 =  1.d0/p2q
      p2qm2 = p2qm1/p2q
      p3qm1 =  1.d0/p3q
      p3qm2 = p3qm1/p3q
      p4qm1 =  1.d0/p4q
      p4qm2 = p4qm1/p4q
c      include 'form/formmemm_OLD.F'
c      print*,'e',vp12,vp34
c      print*,elmat2
      include 'form/formmemm.F' ! this correctly includes vp and vpd...
c      print*,elmat2

      elmat2c=vp342*Z3_ + vp122*Z7_ + vp34*vp12d*Z1_ + vp12*vp34d*Z1_

c      print*,elmat2*e2**3,elmat2c*e2**3
c      print*,vp34*vp12d*Z1_*e2**3
c      print*,vp12*vp34d*Z1_*e2**3
c      print*,vp34*vp12d*Z1_*e2**3+vp12*vp34d*Z1_*e2**3
c      print*,vp12
c      print*,vp34
c      print*,vp12d
c      print*,vp34d


      elmat2 = elmat2 * e2**3
      return
      end
ccccccccccccccccccccccccccccccccccc
****************************************************************************
****************************************************************************
* e+(pp)e-(pm) --> pi+(qp)pi-(qm)g(k)
****************************************************************************
      subroutine formmepp_c(pm,pp,qm,qp,q,elmat2)
c      implicit none
      implicit double precision (Z)
      double precision pm(0:3),pp(0:3),qm(0:3),qp(0:3),q(0:3)
      double precision elmat2,ame,convfac,alpha,pi,e2,qph(40,0:3)
      double precision dot,elmat2cmn,e2approx,eikf,elmat2bck
      double precision me,me2,me4,me6,ammu
      double precision mm,mm2,mm4,mm6
      double precision ampi,mpi,mpi2,mpi4,mpi6
      double precision pm_pp,pm_qm,pm_qp,pp_qm,pp_qp,qm_qp
      double precision pm_pp2,pm_qm2,pm_qp2,pp_qm2,pp_qp2,qm_qp2
      double precision pm_k,pp_k,qm_k,qp_k
      double precision p1q,p2q,p3q,p4q
      double precision p1qm1,p2qm1,p3qm1,p4qm1
      double precision p1qm2,p2qm2,p3qm2,p4qm2
      double precision s12,s34,t24,t13,u14,u23
      double precision s12m1,s34m1,t24m1,t13m1
      double precision s12m2,s34m2,t24m2,t13m2
      double precision vp24,vp13
      double precision vp122,vp342,vp242,vp132
      double precision ptmp(0:3)
      double precision vpol
      external vpol
      double complex vpolc,vp12,vp34
      external vpolc
      integer k
      character*2 fs
      common/finalstate/fs
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      common/testsdebug/elmat2cmn,e2approx,s12,s34,t24,t13,u14,u23
      external dot
      complex*16 im
      complex*16 cfsr,cfsrm,cfsrp,cisrm,cisrp
      complex*16 cfsrd,cfsrmd,cfsrpd,cisrmd,cisrpd
      double precision qppqmpk2,qppk2mmpi2,qmpk2mmpi2
      double precision pmmk2mme2,ppmk2mme2,qppqm2
      complex*16 interf
      complex*16 ffpi12,ffpi34
      double precision ffpi122,ffpi342
      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi
      integer idbgprnt
      common/printfordebug/idbgprnt

      e2 = alpha*4.d0*pi

      me  = ame
      me2 = me*me
      me4 = me2*me2
      me6 = me4*me2

      mm  = ammu
      mm2 = mm*mm
      mm4 = mm2*mm2
      mm6 = mm4*mm2

      mpi  = ampi
      mpi2 = mpi*mpi
      mpi4 = mpi2*mpi2
      mpi6 = mpi4*mpi2

      s12 = dot(pm+pp,pm+pp)
      s12m1 = 1.d0/s12
      s12m2 = s12m1/s12
      vp12 = vpolc(s12)
      vp122 = abs(vp12)**2
      ffpi12= (1.d0,0.d0)
      ffpi122= 1.d0
      if(iffpi.gt.0) call ffpi(s12,ffpi12,ffpi122)
*
      s34 = dot(qm+qp,qm+qp)
      s34m1 = 1.d0/s34
      s34m2 = s34m1/s34!!  * 0.d0
      vp34 = vpolc(s34)
      vp342 = abs(vp34)**2
      ffpi34= (1.d0,0.d0)
      ffpi342= 1.d0
      if(iffpi.gt.0) call ffpi(s34,ffpi34,ffpi342)
*
      qppqm2= s34
      qppqmpk2   = dot(qm+qp+q,qm+qp+q)
c      qppk2mmpi2 = dot(qp+q,qp+q) - mpi2
c      qmpk2mmpi2 = dot(qm+q,qm+q) - mpi2
c      pmmk2mme2  = dot(pm-q,pm-q) - me2
c      ppmk2mme2  = dot(pp-q,pp-q) - me2
      qppk2mmpi2 =  2.d0*dot(qp,q)
      qmpk2mmpi2 =  2.d0*dot(qm,q)
      pmmk2mme2  = -2d0*dot(pm,q)
      ppmk2mme2  = -2d0*dot(pp,q)

*
      pm_pp  = dot(pm,pp)
      pm_pp2 = pm_pp*pm_pp

      pm_qm  = dot(pm,qm)
      pm_qm2 = pm_qm*pm_qm

      pm_qp  = dot(pm,qp)
      pm_qp2 = pm_qp*pm_qp

      pp_qm  = dot(pp,qm)
      pp_qm2 = pp_qm*pp_qm

      pp_qp  = dot(pp,qp)
      pp_qp2 = pp_qp*pp_qp

      qm_qp  = dot(qm,qp)
      qm_qp2 = qm_qp*qm_qp

      pm_k  = dot(pm,q)
      pp_k  = dot(pp,q)
      qm_k  = dot(qm,q)
      qp_k  = dot(qp,q)

      im= (0.d0,1.d0)
      include 'form/formmepp_c.f'

      elmat2 = dble(interf)
      elmat2 = elmat2 * e2**3

      idbgprnt = 0
      if(elmat2.lt.0.d0) then
         print*,'p0 conservation= ',pm(0)+pp(0),qm(0)+qp(0)+q(0)
         print*,'px conservation= ',pm(1)+pp(1),qm(1)+qp(1)+q(1)
         print*,'py conservation= ',pm(2)+pp(2),qm(2)+qp(2)+q(2)
         print*,'pz conservation= ',pm(3)+pp(3),qm(3)+qp(3)+q(3)
         print*,'pm(0)^2= ',sqrt(pm(0)**2-pm(1)**2-pm(2)**2-pm(3)**2)
         print*,'pp(0)^2= ',sqrt(pp(0)**2-pp(1)**2-pp(2)**2-pp(3)**2)
         print*,'qm(0)^2= ',sqrt(qm(0)**2-qm(1)**2-qm(2)**2-qm(3)**2)
         print*,'qp(0)^2= ',sqrt(qp(0)**2-qp(1)**2-qp(2)**2-qp(3)**2)
         print*,'q(0)^2 = ',sqrt(q(0)**2-q(1)**2-q(2)**2-q(3)**2)
         print*,'real elmat2= ',elmat2
c     elmat2= 0.d0
         idbgprnt = 1
      endif

      return
      end

****************************************************************************
* e+(pp)e-(pm) --> pi+(qp)pi-(qm)g(k)
****************************************************************************
*     ATTENTION!!!!! when multiplying by vacuum polarization pay attention
*         to the imaginary parts!!!!
****************************************************************************
c$$$      subroutine formmepp(pm,pp,qm,qp,q,elmat2)
c$$$      implicit none
c$$$      double precision pm(0:3),pp(0:3),qm(0:3),qp(0:3),q(0:3)
c$$$      double precision elmat2,ame,convfac,alpha,pi,e2,qph(40,0:3)
c$$$      double precision dot,elmat2cmn,e2approx,eikf,elmat2bck
c$$$      real*16 qdot
c$$$      double precision me,me2,me4,me6,ammu
c$$$      double precision mm,mm2,mm4,mm6
c$$$      double precision ampi,mpi,mpi2,mpi4,mpi6
c$$$      double precision pm_pp,pm_qm,pm_qp,pp_qm,pp_qp,qm_qp
c$$$      double precision pm_pp2,pm_qm2,pm_qp2,pp_qm2,pp_qp2,qm_qp2
c$$$      double precision pm_k,pp_k,qm_k,qp_k
c$$$      double precision p1q,p2q,p3q,p4q
c$$$      double precision p1qm1,p2qm1,p3qm1,p4qm1
c$$$      double precision p1qm2,p2qm2,p3qm2,p4qm2
c$$$      double precision s12,s34,t24,t13,u14,u23
c$$$      double precision s12m1,s34m1,t24m1,t13m1
c$$$      double precision s12m2,s34m2,t24m2,t13m2
c$$$      double precision vp24,vp13
c$$$      double precision vp122,vp342,vp242,vp132
c$$$      double precision ptmp(0:3)
c$$$      double precision vpol
c$$$      external vpol
c$$$      double complex vpolc,vp12,vp34
c$$$      external vpolc
c$$$      integer k
c$$$      character*2 fs
c$$$      common/finalstate/fs
c$$$      common/parameters/ame,ammu,convfac,alpha,pi
c$$$      common/hparam/ampi
c$$$      common/testsdebug/elmat2cmn,e2approx,s12,s34,t24,t13,u14,u23
c$$$      external dot
c$$$      external qdot
c$$$      complex*16 im
c$$$      complex*16 cfsr,cfsrm,cfsrp,cisrm,cisrp
c$$$      complex*16 cfsrd,cfsrmd,cfsrpd,cisrmd,cisrpd
c$$$      complex*16 interf
c$$$      double precision qppqmpk2,qppk2mmpi2,qmpk2mmpi2
c$$$      double precision pmmk2mme2,ppmk2mme2,qppqm2
c$$$      double precision pm_k2,pp_k2,qm_k2,qp_k2
c$$$      real*16 qpm(0:3),qpp(0:3),qqm(0:3),qqp(0:3),qq(0:3)
c$$$      real*16 ptmps
c$$$      real*8 qelmat2
c$$$*
c$$$      e2 = alpha*4.d0*pi
c$$$*
c$$$      me  = ame
c$$$      me2 = me*me
c$$$      me4 = me2*me2
c$$$      me6 = me4*me2
c$$$
c$$$      mm  = ammu
c$$$      mm2 = mm*mm
c$$$      mm4 = mm2*mm2
c$$$      mm6 = mm4*mm2
c$$$
c$$$      mpi  = ampi
c$$$      mpi2 = mpi*mpi
c$$$      mpi4 = mpi2*mpi2
c$$$      mpi6 = mpi4*mpi2
c$$$*
c$$$      do k = 0,3
c$$$         ptmp(k) = pm(k)+pp(k)
c$$$      enddo
c$$$      s12 = dot(ptmp,ptmp)
c$$$      s12m1 = 1.d0/s12
c$$$      s12m2 = s12m1/s12
c$$$      vp12 = vpolc(s12)
c$$$      vp122 = abs(vp12)**2
c$$$*
c$$$      do k = 0,3
c$$$         ptmp(k) = qm(k)+qp(k)
c$$$      enddo
c$$$      s34 = dot(ptmp,ptmp)
c$$$      s34m1 = 1.d0/s34
c$$$      s34m2 = s34m1/s34!!  * 0.d0
c$$$      vp34 = vpolc(s34)
c$$$      vp342 = abs(vp34)**2
c$$$*
c$$$      do k = 0,3
c$$$         ptmp(k) = qm(k)+qp(k)
c$$$      enddo
c$$$      qppqm2= dot(ptmp,ptmp)
c$$$      do k = 0,3
c$$$         ptmp(k) = qm(k)+qp(k)+q(k)
c$$$      enddo
c$$$      qppqmpk2= dot(ptmp,ptmp)
c$$$      do k = 0,3
c$$$         ptmp(k) = qp(k)+q(k)
c$$$      enddo
c$$$      qppk2mmpi2= dot(ptmp,ptmp) - mpi2
c$$$      do k = 0,3
c$$$         ptmp(k) = qm(k)+q(k)
c$$$      enddo
c$$$      qmpk2mmpi2= dot(ptmp,ptmp) - mpi2
c$$$      do k = 0,3
c$$$         ptmp(k) = pm(k)-q(k)
c$$$      enddo
c$$$      pmmk2mme2= dot(ptmp,ptmp) - me2
c$$$      do k = 0,3
c$$$         ptmp(k) = pp(k)-q(k)
c$$$      enddo
c$$$      ppmk2mme2= dot(ptmp,ptmp) - me2
c$$$*
c$$$      pm_pp  = dot(pm,pp)
c$$$      pm_pp2 = pm_pp*pm_pp
c$$$
c$$$      pm_qm  = dot(pm,qm)
c$$$      pm_qm2 = pm_qm*pm_qm
c$$$
c$$$      pm_qp  = dot(pm,qp)
c$$$      pm_qp2 = pm_qp*pm_qp
c$$$
c$$$      pp_qm  = dot(pp,qm)
c$$$      pp_qm2 = pp_qm*pp_qm
c$$$
c$$$      pp_qp  = dot(pp,qp)
c$$$      pp_qp2 = pp_qp*pp_qp
c$$$
c$$$      qm_qp  = dot(qm,qp)
c$$$      qm_qp2 = qm_qp*qm_qp
c$$$
c$$$      pm_k  = dot(pm,q)
c$$$      pp_k  = dot(pp,q)
c$$$      qm_k  = dot(qm,q)
c$$$      qp_k  = dot(qp,q)
c$$$
c$$$      pm_k2  = pm_k**2
c$$$      pp_k2  = pp_k**2
c$$$      qm_k2  = qm_k**2
c$$$      qp_k2  = qp_k**2
c$$$*
c$$$      im= (0.d0,1.d0)
c$$$      include 'form/formmepp.f'
c$$$
c$$$      elmat2 = dble(interf)
c$$$
c$$$      elmat2 = elmat2 * e2**3
c$$$
c$$$      if(elmat2.lt.0.d0) then
c$$$c$$$         print*,' '
c$$$c$$$         print*,'p0 conservation= ',pm(0)+pp(0),qm(0)+qp(0)+q(0)
c$$$c$$$         print*,'px conservation= ',pm(1)+pp(1),qm(1)+qp(1)+q(1)
c$$$c$$$         print*,'py conservation= ',pm(2)+pp(2),qm(2)+qp(2)+q(2)
c$$$c$$$         print*,'pz conservation= ',pm(3)+pp(3),qm(3)+qp(3)+q(3)
c$$$c$$$         print*,'pm(0)^2= ',sqrt(pm(0)**2-pm(1)**2-pm(2)**2-pm(3)**2)
c$$$c$$$         print*,'pp(0)^2= ',sqrt(pp(0)**2-pp(1)**2-pp(2)**2-pp(3)**2)
c$$$c$$$         print*,'qm(0)^2= ',sqrt(qm(0)**2-qm(1)**2-qm(2)**2-qm(3)**2)
c$$$c$$$         print*,'qp(0)^2= ',sqrt(qp(0)**2-qp(1)**2-qp(2)**2-qp(3)**2)
c$$$c$$$         print*,'q(0)^2 = ',q(0)**2-q(1)**2-q(2)**2-q(3)**2
c$$$c$$$         print*,'q(0)   = ',q(0)
c$$$c$$$         print*,'q(1)   = ',q(1)
c$$$c$$$         print*,'q(2)   = ',q(2)
c$$$c$$$         print*,'q(3)   = ',q(3)
c$$$c$$$         print*,'real double precision elmat2= ',elmat2
c$$$         qelmat2= 0.d0
c$$$*
c$$$* calculating the matrix element with extended precision
c$$$*
c$$$         do k=0,3
c$$$            qpm(k)= pm(k)*1.q0
c$$$            qpp(k)= pp(k)*1.q0
c$$$            qqm(k)= qm(k)*1.q0
c$$$            qqp(k)= qp(k)*1.q0
c$$$            qq(k)= q(k)*1.q0
c$$$         enddo
c$$$         call formmepp_q(qpm,qpp,qqm,qqp,qq,qelmat2)
c$$$         elmat2 = qelmat2 * e2**3
c$$$c$$$         print*,'qp0 cons= ',qpm(0)+qpp(0),qqm(0)+qqp(0)+qq(0)
c$$$c$$$         print*,'qpx cons= ',qpm(1)+qpp(1),qqm(1)+qqp(1)+qq(1)
c$$$c$$$         print*,'qpy cons= ',qpm(2)+qpp(2),qqm(2)+qqp(2)+qq(2)
c$$$c$$$         print*,'qpz cons= ',qpm(3)+pp(3),qqm(3)+qqp(3)+qq(3)
c$$$c$$$         ptmps= qdot(qpm,qpm)
c$$$c$$$         print*,'qpm(0)^2= ',sqrt(ptmps)
c$$$c$$$         ptmps= qdot(qpp,qpp)
c$$$c$$$         print*,'qpp(0)^2= ',sqrt(ptmps)
c$$$c$$$         ptmps= qdot(qqm,qqm)
c$$$c$$$         print*,'qqm(0)^2= ',sqrt(ptmps)
c$$$c$$$         ptmps= qdot(qqp,qqp)
c$$$c$$$         print*,'qqp(0)^2= ',sqrt(ptmps)
c$$$c$$$         print*,'qq(0)^2 = ',qdot(qq,qq)
c$$$c$$$         print*,'qq(0)   = ',qq(0)
c$$$c$$$         print*,'qq(1)   = ',qq(1)
c$$$c$$$         print*,'qq(2)   = ',qq(2)
c$$$c$$$         print*,'qq(3)   = ',qq(3)
c$$$c$$$         print*,'real real*16 elmat2= ',elmat2
c$$$      endif
c$$$      return
c$$$      end
****************************************************************************
* e+(pp)e-(pm) --> pi+(qp)pi-(qm)g(k) with extended precision
****************************************************************************
      subroutine formmepp_q(qpm,qpp,qqm,qqp,qq,qelmat2)
      implicit none
      real*16 qpm(0:3),qpp(0:3),qqm(0:3),qqp(0:3),qq(0:3)
      real*8 qelmat2
      real*16 qdot
      external qdot
      real*16 me2,me4,mpi2,mpi4
      real*16 pm_pp,pm_qm,pm_qp,pp_qm,pp_qp,qm_qp
      real*16 pm_pp2,pm_qm2,pm_qp2,pp_qm2,pp_qp2,qm_qp2
      real*16 pm_k,pp_k,qm_k,qp_k
      real*16 s12,s34
      real*16 vp122,vp342,vp12,vp34
      real*16 ptmp(0:3)
      real*8 ds12,ds34
      real*8 vpol
      external vpol
      complex*16 vpolc,dvp12,dvp34
      external vpolc
      integer k
      complex*32 im
      complex*32 cfsr,cfsrm,cfsrp,cisrm,cisrp
      complex*32 cfsrd,cfsrmd,cfsrpd,cisrmd,cisrpd
      complex*32 interf
      real*16 qppqmpk2,qppk2mmpi2,qmpk2mmpi2
      real*16 pmmk2mme2,ppmk2mme2,qppqm2
      real*16 pm_k2,pp_k2,qm_k2,qp_k2
*
      me2= qpm(0)**2-qpm(1)**2-qpm(2)**2-qpm(3)**2
      me4 = me2*me2

      print*,'do not use this...'

      return

c$$$      mpi2 = qqm(0)**2-qqm(1)**2-qqm(2)**2-qqm(3)**2
c$$$      mpi4 = mpi2*mpi2
c$$$*
c$$$      do k = 0,3
c$$$         ptmp(k) = qpm(k)+qpp(k)
c$$$      enddo
c$$$      s12 = qdot(ptmp,ptmp)
c$$$      ds12= s12*1.d0
c$$$      dvp12 = vpolc(ds12)
c$$$      vp12= dvp12*1.q0
c$$$      vp122 = abs(dvp12)**2*1.q0
c$$$*
c$$$      do k = 0,3
c$$$         ptmp(k) = qqm(k)+qqp(k)
c$$$      enddo
c$$$      s34 = qdot(ptmp,ptmp)
c$$$      ds34 = s34*1.d0
c$$$      dvp34 = vpolc(ds34)
c$$$      vp34= dvp34*1.q0
c$$$      vp342 = abs(dvp34)**2 * 1.q0
c$$$*
c$$$      do k = 0,3
c$$$         ptmp(k) = qqm(k)+qqp(k)
c$$$      enddo
c$$$      qppqm2= qdot(ptmp,ptmp)
c$$$      do k = 0,3
c$$$         ptmp(k) = qqm(k)+qqp(k)+qq(k)
c$$$      enddo
c$$$      qppqmpk2= qdot(ptmp,ptmp)
c$$$      do k = 0,3
c$$$         ptmp(k) = qqp(k)+qq(k)
c$$$      enddo
c$$$      qppk2mmpi2= qdot(ptmp,ptmp) - mpi2
c$$$      do k = 0,3
c$$$         ptmp(k) = qqm(k)+qq(k)
c$$$      enddo
c$$$      qmpk2mmpi2= qdot(ptmp,ptmp) - mpi2
c$$$      do k = 0,3
c$$$         ptmp(k) = qpm(k)-qq(k)
c$$$      enddo
c$$$      pmmk2mme2= qdot(ptmp,ptmp) - me2
c$$$      do k = 0,3
c$$$         ptmp(k) = qpp(k)-qq(k)
c$$$      enddo
c$$$      ppmk2mme2= qdot(ptmp,ptmp) - me2
c$$$*
c$$$      pm_pp  = qdot(qpm,qpp)
c$$$      pm_pp2 = pm_pp*pm_pp
c$$$
c$$$      pm_qm  = qdot(qpm,qqm)
c$$$      pm_qm2 = pm_qm*pm_qm
c$$$
c$$$      pm_qp  = qdot(qpm,qqp)
c$$$      pm_qp2 = pm_qp*pm_qp
c$$$
c$$$      pp_qm  = qdot(qpp,qqm)
c$$$      pp_qm2 = pp_qm*pp_qm
c$$$
c$$$      pp_qp  = qdot(qpp,qqp)
c$$$      pp_qp2 = pp_qp*pp_qp
c$$$
c$$$      qm_qp  = qdot(qqm,qqp)
c$$$      qm_qp2 = qm_qp*qm_qp
c$$$
c$$$      pm_k  = qdot(qpm,qq)
c$$$      pp_k  = qdot(qpp,qq)
c$$$      qm_k  = qdot(qqm,qq)
c$$$      qp_k  = qdot(qqp,qq)
c$$$
c$$$      pm_k2  = pm_k**2
c$$$      pp_k2  = pp_k**2
c$$$      qm_k2  = qm_k**2
c$$$      qp_k2  = qp_k**2
c$$$*
c$$$      im= (0.q0,1.q0)
c$$$      interf= 0.q0
c$$$      print*,'removed the include here!!'
c$$$cc      include 'form/formmepp_q.f'
c$$$
c$$$      qelmat2 = dble(interf)
c$$$
c$$$      return
      end
****************************************************************************
* e+(pp)e-(pm) --> pi+(qp)pi-(qm)g(k)
****************************************************************************
      subroutine formmepp_isr(pm,pp,qm,qp,q,elmat2)
c     implicit none
      implicit double precision (Z)
      double precision pm(0:3),pp(0:3),qm(0:3),qp(0:3),q(0:3)
      double precision elmat2,ame,convfac,alpha,pi,e2,qph(40,0:3)
      double precision dot,elmat2cmn,e2approx,eikf,elmat2bck
      double precision me,me2,me4,me6,ammu
      double precision mm,mm2,mm4,mm6
      double precision ampi,mpi,mpi2,mpi4,mpi6
      double precision pm_pp,pm_qm,pm_qp,pp_qm,pp_qp,qm_qp
      double precision pm_pp2,pm_qm2,pm_qp2,pp_qm2,pp_qp2,qm_qp2
      double precision pm_k,pp_k,qm_k,qp_k,qm_k2,qp_k2
      double precision p1q,p2q,p3q,p4q
      double precision p1qm1,p2qm1,p3qm1,p4qm1
      double precision p1qm2,p2qm2,p3qm2,p4qm2
      double precision s12,s34,t24,t13,u14,u23
      double precision s12m1,s34m1,t24m1,t13m1
      double precision s12m2,s34m2,t24m2,t13m2
      double precision vp24,vp13
      double precision vp122,vp342,vp242,vp132
      double precision ptmp(0:3)
      double precision vpol
      external vpol
      double complex vpolc,vp12,vp34
      external vpolc
      integer k
      character*2 fs
      common/finalstate/fs
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      common/testsdebug/elmat2cmn,e2approx,s12,s34,t24,t13,u14,u23
      external dot
      complex*16 cfsr,cfsrm,cfsrp,cisrm,cisrp,cisrm2,cisrp2
      complex*16 cfsrd,cfsrmd,cfsrpd,cisrmd,cisrpd
      complex*16 interf
      double precision qppqmpk2,qppk2mmpi2,qmpk2mmpi2
      double precision pmmk2mme2,ppmk2mme2,qppqm2
      double precision s,t,s2,t2
      complex*16 ffpi34
      real*8 ffpi342
      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi

      e2 = alpha*4.d0*pi

      me  = ame
      me2 = me*me
      me4 = me2*me2
      me6 = me4*me2

      mm  = ammu
      mm2 = mm*mm
      mm4 = mm2*mm2
      mm6 = mm4*mm2

      mpi  = ampi
      mpi2 = mpi*mpi
      mpi4 = mpi2*mpi2
      mpi6 = mpi4*mpi2
*
      s = dot(pm+pp,pm+pp)
      s2= s*s
      t = dot(pm-qm,pm-qm)
      t2= t*t
      s34 = dot(qm+qp,qm+qp)
      s34m1 = 1.d0/s34
      s34m2 = s34m1/s34!!  * 0.d0
      vp34 = vpolc(s34)
      vp342 = abs(vp34)**2
*
      ffpi342= 1.d0
      if(iffpi.gt.0) call ffpi(s34,ffpi34,ffpi342)
*
      qppqm2= s34
      qppqmpk2= dot(qm+qp+q,qm+qp+q)
c      qppk2mmpi2= dot(qp+q,qp+q) - mpi2
c      qmpk2mmpi2= dot(qm+q,qm+q) - mpi2
c      pmmk2mme2= dot(pm-q,pm-q) - me2
c      ppmk2mme2= dot(pp-q,pp-q) - me2
      qppk2mmpi2 =  2.d0*dot(qp,q)
      qmpk2mmpi2 =  2.d0*dot(qm,q)
      pmmk2mme2  = -2d0*dot(pm,q)
      ppmk2mme2  = -2d0*dot(pp,q)
*
      pm_pp  = dot(pm,pp)
      pm_pp2 = pm_pp*pm_pp

      pm_qm  = dot(pm,qm)
      pm_qm2 = pm_qm*pm_qm

      pm_qp  = dot(pm,qp)
      pm_qp2 = pm_qp*pm_qp

      pp_qm  = dot(pp,qm)
      pp_qm2 = pp_qm*pp_qm

      pp_qp  = dot(pp,qp)
      pp_qp2 = pp_qp*pp_qp

      qm_qp  = dot(qm,qp)
      qm_qp2 = qm_qp*qm_qp

      pm_k  = dot(pm,q)
      pp_k  = dot(pp,q)
      qm_k  = dot(qm,q)
      qm_k2 = qm_k**2
      qp_k  = dot(qp,q)
      qp_k2 = qp_k**2

c      p1qm1 =  1.d0/p1q
c      p1qm2 = p1qm1/p1q
c      p2qm1 =  1.d0/p2q
c      p2qm2 = p2qm1/p2q
c      p3qm1 =  1.d0/p3q
c      p3qm2 = p3qm1/p3q
c      p4qm1 =  1.d0/p4q
c      p4qm2 = p4qm1/p4q

      include 'form/formmepp_isr.f'

      elmat2 = dble(interf)
      elmat2= elmat2*ffpi342

      elmat2 = elmat2 * e2**3

      if(elmat2.lt.0.d0) then
         print*,'p0 conservation= ',pm(0)+pp(0),qm(0)+qp(0)+q(0)
         print*,'px conservation= ',pm(1)+pp(1),qm(1)+qp(1)+q(1)
         print*,'py conservation= ',pm(2)+pp(2),qm(2)+qp(2)+q(2)
         print*,'pz conservation= ',pm(3)+pp(3),qm(3)+qp(3)+q(3)
         print*,'pm(0)^2= ',sqrt(pm(0)**2-pm(1)**2-pm(2)**2-pm(3)**2)
         print*,'pp(0)^2= ',sqrt(pp(0)**2-pp(1)**2-pp(2)**2-pp(3)**2)
         print*,'qm(0)^2= ',sqrt(qm(0)**2-qm(1)**2-qm(2)**2-qm(3)**2)
         print*,'qp(0)^2= ',sqrt(qp(0)**2-qp(1)**2-qp(2)**2-qp(3)**2)
         print*,'q(0)^2 = ',q(0)**2-q(1)**2-q(2)**2-q(3)**2
         print*,'q(0)   = ',q(0)
         print*,'q(1)   = ',q(1)
         print*,'q(2)   = ',q(2)
         print*,'q(3)   = ',q(3)
         print*,'real elmat2= ',elmat2
         elmat2= 0.d0
      endif

      return
      end

****************************************************************************
* e+(pp)e-(pm) --> pi+(qp)pi-(qm)g(k)
****************************************************************************
      subroutine formmepp_isr_c(pm,pp,qm,qp,q,elmat2)
      implicit none
      double precision pm(0:3),pp(0:3),qm(0:3),qp(0:3),q(0:3)
      double precision elmat2,ame,convfac,alpha,pi,e2,qph(40,0:3)
      double precision dot,elmat2cmn,e2approx,eikf,elmat2bck
      double precision me,me2,me4,me6,ammu
      double precision mm,mm2,mm4,mm6
      double precision ampi,mpi,mpi2,mpi4,mpi6
      double precision pm_pp,pm_qm,pm_qp,pp_qm,pp_qp,qm_qp
      double precision pm_pp2,pm_qm2,pm_qp2,pp_qm2,pp_qp2,qm_qp2
      double precision pm_k,pp_k,qm_k,qp_k
      double precision p1q,p2q,p3q,p4q
      double precision p1qm1,p2qm1,p3qm1,p4qm1
      double precision p1qm2,p2qm2,p3qm2,p4qm2
      double precision s12,s34,t24,t13,u14,u23
      double precision s12m1,s34m1,t24m1,t13m1
      double precision s12m2,s34m2,t24m2,t13m2
      double precision vp24,vp13
      double precision vp122,vp342,vp242,vp132
      double precision ptmp(0:3)
      double precision vpol
      external vpol
      double complex vpolc,vp12,vp34
      external vpolc
      integer k
      character*2 fs
      common/finalstate/fs
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      common/testsdebug/elmat2cmn,e2approx,s12,s34,t24,t13,u14,u23
      external dot
      complex*16 im
      complex*16 cfsr,cfsrm,cfsrp,cisrm,cisrp
      complex*16 cfsrd,cfsrmd,cfsrpd,cisrmd,cisrpd
      complex*16 interf
      double precision qppqmpk2,qppk2mmpi2,qmpk2mmpi2
      double precision pmmk2mme2,ppmk2mme2,qppqm2


      e2 = alpha*4.d0*pi

      me  = ame
      me2 = me*me
      me4 = me2*me2
      me6 = me4*me2

      mm  = ammu
      mm2 = mm*mm
      mm4 = mm2*mm2
      mm6 = mm4*mm2

      mpi  = ampi
      mpi2 = mpi*mpi
      mpi4 = mpi2*mpi2
      mpi6 = mpi4*mpi2


      do k = 0,3
         ptmp(k) = pm(k)+pp(k)
      enddo
      s12 = dot(ptmp,ptmp)
      s12m1 = 1.d0/s12
      s12m2 = s12m1/s12
      vp12 = vpolc(s12)
      vp122 = abs(vp12)**2

      s34 = dot(qm+qp,qm+qp)
      s34m1 = 1.d0/s34
      s34m2 = s34m1/s34!!  * 0.d0
      vp34 = vpolc(s34)
      vp342 = abs(vp34)**2
*
      qppqm2= s34
      qppqmpk2= dot(qm+qp+q,qm+qp+q)
c      qppk2mmpi2= dot(qp+q,qp+q) - mpi2
c      qmpk2mmpi2= dot(qm+q,qm+q) - mpi2
c      pmmk2mme2= dot(pm-q,pm-q) - me2
c      ppmk2mme2= dot(pp-q,pp-q) - me2
      qppk2mmpi2 =  2.d0*dot(qp,q)
      qmpk2mmpi2 =  2.d0*dot(qm,q)
      pmmk2mme2  = -2d0*dot(pm,q)
      ppmk2mme2  = -2d0*dot(pp,q)
*
      pm_pp  = dot(pm,pp)
      pm_pp2 = pm_pp*pm_pp

      pm_qm  = dot(pm,qm)
      pm_qm2 = pm_qm*pm_qm

      pm_qp  = dot(pm,qp)
      pm_qp2 = pm_qp*pm_qp

      pp_qm  = dot(pp,qm)
      pp_qm2 = pp_qm*pp_qm

      pp_qp  = dot(pp,qp)
      pp_qp2 = pp_qp*pp_qp

      qm_qp  = dot(qm,qp)
      qm_qp2 = qm_qp*qm_qp

      pm_k  = dot(pm,q)
      pp_k  = dot(pp,q)
      qm_k  = dot(qm,q)
      qp_k  = dot(qp,q)

c      p1qm1 =  1.d0/p1q
c      p1qm2 = p1qm1/p1q
c      p2qm1 =  1.d0/p2q
c      p2qm2 = p2qm1/p2q
c      p3qm1 =  1.d0/p3q
c      p3qm2 = p3qm1/p3q
c      p4qm1 =  1.d0/p4q
c      p4qm2 = p4qm1/p4q

      im= (0.d0,1.d0)
      include 'form/formmepp_isr_c.f'

      elmat2 = dble(interf)

      elmat2 = elmat2 * e2**3

      if(elmat2.lt.0.d0) then
         print*,'p0 conservation= ',pm(0)+pp(0),qm(0)+qp(0)+q(0)
         print*,'px conservation= ',pm(1)+pp(1),qm(1)+qp(1)+q(1)
         print*,'py conservation= ',pm(2)+pp(2),qm(2)+qp(2)+q(2)
         print*,'pz conservation= ',pm(3)+pp(3),qm(3)+qp(3)+q(3)
         print*,'pm(0)^2= ',sqrt(pm(0)**2-pm(1)**2-pm(2)**2-pm(3)**2)
         print*,'pp(0)^2= ',sqrt(pp(0)**2-pp(1)**2-pp(2)**2-pp(3)**2)
         print*,'qm(0)^2= ',sqrt(qm(0)**2-qm(1)**2-qm(2)**2-qm(3)**2)
         print*,'qp(0)^2= ',sqrt(qp(0)**2-qp(1)**2-qp(2)**2-qp(3)**2)
         print*,'q(0)^2 = ',sqrt(q(0)**2-q(1)**2-q(2)**2-q(3)**2)
         print*,'real elmat2= ',elmat2
         elmat2= 0.d0
      endif

      return
      end
****************************************************************************
      subroutine born_me(p1,p2,p3,p4,smt,smtq2,smtz,iq2,q)
      implicit double precision (a-h,l,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),ptmp(0:3)
      dimension q(0:3)
      complex*16 chis,chit,im
      parameter (im = (0.d0,1.d0))
      double complex vpolc,ec2s
      external vpolc
      character*2 fs
      common/finalstate/fs
      common/parameters/ame,ammu,convfac,alpha,pi
      common/zparameters/zm,gz,ve,ae,rv,ra,wm,s2th,gfermi,sqrt2,um
      common/bornforsvgg/borngg

      integer ifirst
      common/ifirstborn/ifirst
      data ifirst /0/

      common/hparam/ampi

      double complex ff
      character*10 what_ffpi
      integer iffpi
      common/ffpioption/iffpi,what_ffpi
      integer isQEDHO
      common/sQEDoption/isQEDHO

      smt   = 0.d0
      smtq2 = 0.d0
      smtz  = 0.d0

      if (fs.eq.'gg') then
         me2 = ame*ame
         me4 = me2*me2
         p1p2 = dot(p1,p2)
         p1p3 = dot(p1,p3)
         p1p4 = dot(p1,p4)
         p2p3 = dot(p2,p3)
         p2p4 = dot(p2,p4)
         p3p4 = dot(p3,p4)
         p2p3m1 = 1.d0/p2p3
         p2p4m1 = 1.d0/p2p4
         include 'form/borngg.F'
         e2 = 4.d0*pi*alpha
         smt = elmat2*e2*e2
         smtq2 = smt

         borngg = smt

c         S = dot(p1+p2,p1+p2)
c         beta = p1(3)/p1(0)
c         c = p4(3)/sqrt(tridot(p4,p4))
c         c2 = c*c
c         ALPHA2=ALPHA*ALPHA
c         BORN1=2.D0*PI*ALPHA2/S*BETA*BETA/
c     >        (1.D0+BETA*C)/(1.D0-BETA*C)*(1.D0+C2)
c         phsp = 3.1662869888230555d-3
c         print*,born1/smt/phsp*(8.d0 * (sqrt(s)/2.d0)**2)/(2.d0*pi/4.d0)
         return
      endif

      s = dot(p1+p2,p1+p2)
      t = dot(p1-p3,p1-p3)
      u = dot(p1-p4,p1-p4)

      em2 = ame**2
      em4 = em2*em2
      ec2 = 4.d0*pi*alpha       ! electron charge^2
      ec4 = ec2*ec2

      s2 = s*s
      t2 = t*t
      u2 = u*u

*** NEW : always running !!!
      if (fs.eq.'ee') then

         if (ifirst.eq.0) then
c     print*,' ONLY INTERFERENCE!!'

            print*,'Z exchange is OFF!!!'
            print*,'Z exchange is OFF!!!'
            print*,'Z exchange is OFF!!!'
            ifirst = 1
         endif

         ec2t  = ec2 * vpolc(t)
         ec2s  = ec2 * vpolc(s) ! * 0.d0

         sm2   = abs(ec2s)**2/s2*((t2+u2)/4.d0+2.d0*em2*s
     .        -2.d0*em4)
         tm2   = ec2t**2/t2*((s2+u2)/4.d0+2.d0*em2*t
     .        -2.d0*em4)
         st2   = ec2s*ec2t/s/t*(u2/2.d0-4.d0*em2*u+6.d0*em4 )

         smtq2 = 32.d0*(tm2 + sm2 + st2)
c         smtq2 = 32.d0*(st2)
         smt   = smtq2

      elseif (fs.eq.'mm') then
         !! WARNING::: NO Z EXCHANGE IN MU MU !!!
         mm2 = ammu*ammu
         ec2s  = ec2 * vpolc(s)

         me2   = ame*ame
         e4    = abs(ec2s)**2
         elmat2=32.d0*e4/s2*((t2+u2)/4.d0+s*(me2+mm2)-(me2+mm2)**2/2.d0)
         smtq2 = elmat2
         smt   = smtq2
         return
      elseif (fs.eq.'pp'.or.fs.eq.'pi') then
         mpi2 = ampi*ampi

         ec2s  = ec2 * vpolc(s)

         me2   = ame*ame
         e4    = abs(ec2s)**2
*
         elmat2=e4/s2*(-2.d0*u**2+4.d0*t*u-2.d0*t**2
     +                 +2.d0*s**2-8.d0*mpi2*s)    !(e4/s2 * (2 s^2 - 2 (u-t)^2 - 8 mpi^2 s)
         smtq2 = elmat2
         smt   = smtq2
*
* inclusion of the pion form factor
*
         ff  = (1.d0,0.d0)
         ff2 = 1.d0
         if(iffpi.gt.0) call ffpi(s,ff,ff2)
         smtq2 = smtq2*ff2
         smt   = smtq2

*
         return
      elseif (fs.eq.'mr'.or.fs.eq.'er') then
c         print*, 'passo qui?', smt
         return
      endif

********************
***   with Z exchange  !
      if (iq2.eq.1) then
         goto 666
****
** with Z exchange, which gives a 0.1% effect at 10 GeV for wide acceptance.
** From Caffo & Remiddi in CERN Yellow Report 'Z phyisics at LEP1' Vol 1
         chis=zm**2/16.d0/um**2*(ve**2+ae**2)*s/(s-zm**2+im*gz*zm)
         chit=zm**2/16.d0/um**2*(ve**2+ae**2)*s/2.d0/(zm**2-t-im*gz*zm)

         z = p3(3)/sqrt(tridot(p3,p3))

         z2   = z**2
         upz2 = (1.d0+z)**2

         chi2   = (abs(chis))**2
         chit2  = (abs(chit))**2
         rechi  = (chis + dconjg(chis))/2.d0
         rechit = (chit + dconjg(chit))/2.d0
         rechichit = (chis*dconjg(chit) + dconjg(chis)*chit)/2.d0

         s4  = -2.d0*rechit*upz2
         s5  =  4.d0*rechit/(1.d0-z)*(upz2 + 4.d0*(rv-ra))
         s6  =  2.d0*chit2*
     .        (upz2*(1.d0+4.d0*rv*ra)+4.d0*(1.d0-4.d0*rv*ra))
         s7  =  2.d0*rechi*((1.d0+z2)*rv+2.d0*z*ra)
         s8  = -2.d0*rechi*upz2/(1.d0-z)
         s9  = -2.d0*rechichit*upz2*(1.d0+4.d0*rv*ra)
         s10 =  chi2*((1.d0+z2)+8.d0*z*rv*ra)

         sz = (s4+s5+s6+s7+s8+s9+s10)*alpha**2/4.d0/s
         smtz = sz*32.d0*pi**2 *2.d0*s * 4.d0 !summed, not averaged...

**** ONLY PURE QED
 666     smtz = 0.d0
****
c         call recola_ampl(0,p1,p2,p3,p4,p3,bornrcl,allun)
c         bornrcl = 4.d0*bornrcl

c         print*,(smt+smtz)/bornrcl

      endif
******
** TEST 2 LOOP!
c      esoft = 1.d-4*p1(0)
c      call twoloopvirtualsoft(p1,p2,p3,p4,esoft,twoloopoverborn,pair)
c      smt = pair*smtq2 !!!*pi*pi*p1(0)**2*1024.d0
**** END TEST 2 LOOP
** TEST SOFT V+R PAIR CORRECTIONS
** virtual -- FORMULA 6b jadach
c$$$      ammu = ame
c$$$      L   = log(-t/ammu/ammu)
c$$$      L2  = L*L
c$$$      L3  = L2*L
c$$$      cf  = 383.d0/108.d0-11.d0/6.d0*pi*pi/6.d0
c$$$      z3  = 1.2020569032d0
c$$$      cfp = -z3/3.d0+3355.d0/1296.d0-19.d0/18.d0*pi**2/6.d0
c$$$
c$$$      fmu = -1.d0/36.d0*L3      ! * 0.d0
c$$$     .     +19.d0/72.d0*L2      !* 0.d0
c$$$     .     +(-1.d0/36.d0*pi*pi-265.d0/216.d0)*L+
c$$$     .     cf + 19.d0/75.d0*pi*pi
c$$$
c$$$      v = 4.d0*(alpha/pi)**2*fmu * 32.d0*tm2 ! full virtual for t channel only
c$$$      smt   = v
c$$$      smtq2 = smt
c$$$** real -- Arbuzov et al., formula 43 only t-channel
c$$$      c = p3(3)/sqrt(tridot(p3,p3))
c$$$      L = log(-t/ame/ame) ! - 1.d0 ! ?????
c$$$      L2 = L**2
c$$$      L3 = L2*L
c$$$      delta = 1d0
c$$$c      delta = 1.d0
c$$$      ld = log(delta)
c$$$cccccccccccccccccccccccc      ld = -ld
c$$$      ddlog = ddilog((1.d0+c)/2.d0)
c$$$      at = 56.d0/9.d0 -2.d0/3.d0*pi*pi+2.d0*ddlog
c$$$      r  = 1.d0/3.d0*L3
c$$$     .     +L2*(2.d0*ld-5.d0/3.d0)
c$$$     .     +L*(4.d0*ld**2-20.d0/3.d0*ld+at)
c$$$      r  = r * 1.d0/3.d0*alpha**2/pi/pi
c$$$      r  = r *32.d0*tm2
c$$$      r1 = r
c$$$* tesi Oreste
c$$$      L  = L ! + 1.d0
c$$$      z2 = pi**2/6.d0
c$$$      r  =  1.d0/18.d0*(L-2.d0*ld)**3
c$$$     .     -5.d0/18.d0*(L-2.d0*ld)**2
c$$$     .     +1.d0/9.d0*3.d0/2.d0*(56.d0/9.d0-4.d0*z2 + 0.*2*ddlog) ! Oreste
c$$$     .     *(L-2.d0*ld)
c$$$      r = r*2.d0
c$$$      r = r*alpha**2/pi/pi *32.d0*tm2
c$$$
c$$$ccc ATTENZIONE: QUELLO DI ORESTE AUMENTA AL DIMINUIRE DI DELTA!!!
c$$$ccc A ME SEMBRA POCO FISICO (INTEGRO SU UNO SPAZIO FASI PIU' AMPIO!!)
c$$$
c$$$c      print*,'TEST COPPIE!!!!'
c$$$c      print*,r-r1,L
c$$$c      r = r1
c$$$      smt = r
c$$$      smtq2 = smt
c$$$      return
************ END TEST COPPIE
      return
      end

**************************************************************************
      subroutine eikonal_factor_z(eikng,ng,p1,p2,p3,p4,q,x,qph,eikonal)
      implicit double precision (a-h,o-z)
      character*6 eikng
      parameter (imaxph = 40)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension qph(imaxph,0:3),q(0:3),peik(0:3)
      character*2 fs
      real*8 x
      common/finalstate/fs
      common/singlefactors/singleeikfact(imaxph)
      common/parameters/ame,ammu,convfac,alpha,pi

      logical cond
      common/radoption/ifisrfsr
      integer ifisrfsr

      eikonal = 1.d0
      test = 1.d0
      if (ng.eq.0) return

      nglocal = ng
      if (eikng.eq.'single') nglocal = 1

      ec2 = 4.d0*pi*alpha ! electron charge^2
      do k = 1,nglocal
         if (eikng.eq.'multip') then
            q(0) = qph(k,0)
            q(1) = qph(k,1)
            q(2) = qph(k,2)
            q(3) = qph(k,3)
         endif
         oop1dk = 1.d0/dot(p1,q)
         oop2dk = 1.d0/dot(p2,q)
         if (fs.eq.'ee'.or.fs.eq.'mm'.or.fs.eq.'pp'
     .     .or.fs.eq.'er'.or.fs.eq.'mr'.or.fs.eq.'pr') then
            oop3dk = 1.d0/dot(p3,q)
            oop4dk = 1.d0/dot(p4,q)
         elseif (fs.eq.'gg'.or.fs.eq.'pi') then
            oop3dk = 0.d0
            oop4dk = 0.d0
         endif

         cond=(fs.eq.'mr'.or.fs.eq.'pr')
         if (cond.and.ifisrfsr.eq.10) then
            oop3dk = 0.d0
            oop4dk = 0.d0
         elseif (cond.and.ifisrfsr.eq.01) then
            oop1dk = 0.d0
            oop2dk = 0.d0
         endif

         x1=x
         x2=x

*****FOR PIONS NO FURTHER CORRECTION NEEDED, EIKONAL ALREADY TAKES INTO ACCOUNT THE MASS TERM.
        if(fs.eq.'pr'.or.fs.eq.'pp') x2=0.d0
        singlephoton = -(
     .          (dot(p1,p1)*oop1dk**2
     .        + dot(p2,p2)*oop2dk**2)*(1d0-x1**2/(1d0+(1d0-x1)**2))
     .        + (dot(p3,p3)*oop3dk**2
     .        + dot(p4,p4)*oop4dk**2)*(1d0-x2**2/(1d0+(1d0-x2)**2))
     .        - 2.d0*dot(p1,p2)*oop1dk*oop2dk
     .        - 2.d0*dot(p3,p4)*oop3dk*oop4dk
     .        - 2.d0*dot(p1,p3)*oop1dk*oop3dk
     .        + 2.d0*dot(p1,p4)*oop1dk*oop4dk
     .        + 2.d0*dot(p2,p3)*oop2dk*oop3dk
     .        - 2.d0*dot(p2,p4)*oop2dk*oop4dk)*ec2

         eikonal = eikonal * singlephoton
         if (eikng.eq.'multip') singleeikfact(k) = singlephoton
      enddo

      eikonal = eikonal/factorial(nglocal)
      return
      end
**************************************************************************

**************************************************************************
      subroutine eikonal_factor(eikng,ng,p1,p2,p3,p4,q,qph,eikonal)
      implicit double precision (a-h,o-z)
      character*6 eikng
      parameter (imaxph = 40)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension qph(imaxph,0:3),q(0:3),peik(0:3)
      character*2 fs
      common/finalstate/fs
      common/singlefactors/singleeikfact(imaxph)
      common/parameters/ame,ammu,convfac,alpha,pi

      logical cond
      common/radoption/ifisrfsr
      integer ifisrfsr

      eikonal = 1.d0
      test = 1.d0
      if (ng.eq.0) return

      nglocal = ng
      if (eikng.eq.'single') nglocal = 1

      ec2 = 4.d0*pi*alpha ! electron charge^2
      do k = 1,nglocal
         if (eikng.eq.'multip') then
            q(0) = qph(k,0)
            q(1) = qph(k,1)
            q(2) = qph(k,2)
            q(3) = qph(k,3)
         endif
         oop1dk = 1.d0/dot(p1,q)
         oop2dk = 1.d0/dot(p2,q)
         if (fs.eq.'ee'.or.fs.eq.'mm'.or.fs.eq.'pp'
     .     .or.fs.eq.'er'.or.fs.eq.'mr'.or.fs.eq.'pr') then
            oop3dk = 1.d0/dot(p3,q)
            oop4dk = 1.d0/dot(p4,q)
         elseif (fs.eq.'gg'.or.fs.eq.'pi') then
            oop3dk = 0.d0
            oop4dk = 0.d0
         endif

         cond=(fs.eq.'mr'.or.fs.eq.'pr')
         if (cond.and.ifisrfsr.eq.10) then
            oop3dk = 0.d0
            oop4dk = 0.d0
         elseif (cond.and.ifisrfsr.eq.01) then
            oop1dk = 0.d0
            oop2dk = 0.d0
         endif

         singlephoton = -(
     .          dot(p1,p1)*oop1dk**2
     .        + dot(p2,p2)*oop2dk**2
     .        + dot(p3,p3)*oop3dk**2
     .        + dot(p4,p4)*oop4dk**2
     .        - 2.d0*dot(p1,p2)*oop1dk*oop2dk
     .        - 2.d0*dot(p3,p4)*oop3dk*oop4dk
     .        - 2.d0*dot(p1,p3)*oop1dk*oop3dk
     .        + 2.d0*dot(p1,p4)*oop1dk*oop4dk
     .        + 2.d0*dot(p2,p3)*oop2dk*oop3dk
     .        - 2.d0*dot(p2,p4)*oop2dk*oop4dk)*ec2

         eikonal = eikonal * singlephoton
         if (eikng.eq.'multip') singleeikfact(k) = singlephoton
      enddo

      eikonal = eikonal/factorial(nglocal)
      return
      end
**************************************************************************
*
**************************************************************************
      subroutine eikonal_factor_BCK(eikng,ng,p1,p2,p3,p4,q,qph,eikonal)
      implicit double precision (a-h,o-z)
      character*6 eikng
      parameter (imaxph = 40)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension qph(imaxph,0:3),q(0:3),peik(0:3)
      character*2 fs
      common/finalstate/fs
      common/singlefactors/singleeikfact(imaxph)
      common/parameters/ame,ammu,convfac,alpha,pi

      eikonal = 1.d0
      if (ng.eq.0) return

      nglocal = ng
      if (eikng.eq.'single') nglocal = 1

      ec2 = 4.d0*pi*alpha ! electron charge^2
      do k = 1,nglocal
         if (eikng.eq.'multip') then
            q(0) = qph(k,0)
            q(1) = qph(k,1)
            q(2) = qph(k,2)
            q(3) = qph(k,3)
         endif
         oop1dk = 1.d0/dot(p1,q)
         oop2dk = 1.d0/dot(p2,q)
         if (fs.eq.'ee'.or.fs.eq.'mm'.or.fs.eq.'pp'
     .       .or.fs.eq.'mr'.or.fs.eq.'er') then
            oop3dk = 1.d0/dot(p3,q)
            oop4dk = 1.d0/dot(p4,q)
         elseif (fs.eq.'gg'.or.fs.eq.'pi') then
            oop3dk = 0.d0
            oop4dk = 0.d0
         endif
         do i = 0,3
            peik(i)=p1(i)*oop1dk-p2(i)*oop2dk-p3(i)*oop3dk+p4(i)*oop4dk
         enddo
         singlephoton = -ec2*dot(peik,peik)
         eikonal = eikonal * singlephoton
         if (eikng.eq.'multip') singleeikfact(k) = singlephoton
      enddo
      eikonal = eikonal/factorial(nglocal)
      return
      end
**************************************************************************
*
**************************************************************************
      subroutine eikonal_noint(eikng,ng,p1,p2,p3,p4,q,qph,eikonal)
      implicit double precision (a-h,o-z)
      character*6 eikng
      parameter (imaxph = 40)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),pm(4,0:3)
      dimension qph(imaxph,0:3),q(0:3),peik(0:3),ptmp1(0:3),ptmp2(0:3)
      common/parameters/ame,ammu,convfac,alpha,pi

      eikonal = 1.d0
      if (ng.eq.0) return

      do k = 0,3
         pm(1,k) = p1(k)
         pm(2,k) = p2(k)
         pm(3,k) = p3(k)
         pm(4,k) = p4(k)
      enddo

      nglocal = ng
      if (eikng.eq.'single') nglocal = 1

      ec2 = 4.d0*pi*alpha ! electron charge^2
      do k = 1,nglocal
         singlephoton = 0.d0
         if (eikng.eq.'multip') then
            do i = 0,3
               q(i) = qph(k,i)
            enddo
         endif

         do i = 1,3
            do j = i+1,4
               do jj = 0,3
                  ptmp1(jj) = pm(i,jj)
                  ptmp2(jj) = pm(j,jj)
               enddo

               if (i.eq.1.and.j.eq.2) etaij = -1.d0
               if (i.eq.1.and.j.eq.3) etaij = -1.d0
               if (i.eq.1.and.j.eq.4) etaij =  1.d0
               if (i.eq.2.and.j.eq.3) etaij =  1.d0
               if (i.eq.2.and.j.eq.4) etaij = -1.d0
               if (i.eq.3.and.j.eq.4) etaij = -1.d0

               if (i.eq.1.and.j.eq.2) etaij = -1.d0
               if (i.eq.1.and.j.eq.3) etaij =  0.d0
               if (i.eq.1.and.j.eq.4) etaij =  0.d0
               if (i.eq.2.and.j.eq.3) etaij =  0.d0
               if (i.eq.2.and.j.eq.4) etaij =  0.d0
               if (i.eq.3.and.j.eq.4) etaij = -1.d0


               pipj = dot(ptmp1,ptmp2)
               singlephoton = singlephoton -
     .              etaij*2.d0*pipj/dot(ptmp1,q)/dot(ptmp2,q)
            enddo
         enddo
         singlephoton = singlephoton - ame**2*
     .        (1.d0/dot(p1,q)**2+1.d0/dot(p2,q)**2+
     .         1.d0/dot(p3,q)**2+1.d0/dot(p4,q)**2)
         eikonal = eikonal * singlephoton * ec2
      enddo
      eikonal = eikonal/factorial(nglocal)
      return
      end
****
      subroutine bk_me(p1,p2,p3,p4,q,bkme,qph)
      implicit double precision (a-h,o-z)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),ptmp(0:3)
      double precision qph(40,0:3),qt(0:3),p3t(0:3),p4t(0:3)
      double precision k1(0:3),k2(0:3),k3(0:3)
c      real*16 qp1(0:3),qp2(0:3),qp3(0:3),qp4(0:3),qq(0:3)
      real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
      real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
      real*8 e2,e6
      real*8 bkme1,bkme
      real*8 num,num1,num3,num2,den2,den3,den,den1
      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi

      double complex virtual1
      real*8 elmat,elmat1,maxerr
      character*2 fs,fsbck
      real*8 allun
      common/finalstate/fs

      external dot
      real*8 dot
      common/softphotonmass/phmass

      common/ialpharunning/iarun
      integer iarun

      common/radoption/ifisrfsr
      integer ifisrfsr


      e2 = 4.d0*pi*alpha
      e6 = e2*e2*e2

      fsbck = fs
      if (fs.eq.'mr') then
         fs = 'mm'
      endif
      if (fs.eq.'ee') then
         call formme(p1,p2,p3,p4,q,bkme)

! ok a 1 GeV...         f
c         call recola_ampl(1,p1,p2,p3,p4,q,realrcl,allun)
c         realrcl = realrcl*4.d0
c         print*,realrcl/bkme


      elseif (fs.eq.'mm') then
         call formmemm(p1,p2,p3,p4,q,bkme)
      elseif (fs.eq.'pp') then
         call formmepp_c(p1,p2,p3,p4,q,bkme)
      elseif (fs.eq.'pi') then
         call formmepp_isr(p1,p2,p3,p4,q,bkme)
      elseif (fs.eq.'er') then

	    call eerad_born(p1,p2,p3,p4,q,bkme,qph)

c	    print*, 'bkme', bkme
      elseif (fs.eq.'mr') then
         if(ifisrfsr.eq.11) then
            call ISR_mur(p1,p2,p3,p4,q,virtual,qph)
            bkme=virtual
            call FSR_mur(p1,p2,p3,p4,q,virtual,qph)
            bkme=bkme+virtual
            call IFI_mur(p1,p2,p3,p4,q,virtual,qph)
            bkme=bkme+2d0*virtual

 	 elseif(ifisrfsr.eq.10) then
 	    call ISR_mur(p1,p2,p3,p4,q,virtual,qph)
            bkme=virtual

         elseif(ifisrfsr.eq.01) then
            call FSR_mur(p1,p2,p3,p4,q,virtual,qph)
            bkme=virtual
         endif

      elseif (fs.eq.'pr') then

         if(ifisrfsr.eq.11) then
            call ISR_pr(p1,p2,p3,p4,q,virtual,qph)
            bkme=virtual
!            print*,'isr',virtual
            call FSR_pr(p1,p2,p3,p4,q,virtual,qph)
            bkme=bkme+virtual
!            print*,'fsr',virtual
            call IFI_pr(p1,p2,p3,p4,q,virtual1,qph)
            bkme=bkme+2d0*real(virtual1)
!            print*,'ifi',real(virtual1)

!            stop
 	 elseif(ifisrfsr.eq.10) then
            call ISR_pr(p1,p2,p3,p4,q,virtual,qph)
            bkme=virtual

         elseif(ifisrfsr.eq.01) then
            call FSR_pr(p1,p2,p3,p4,q,virtual,qph)
            bkme=virtual
         endif

      elseif (fs.eq.'gg') then

! using exact ME calculated from scratch
         call ee3gexact(p1,p2,p3,p4,q,bkme)

      endif

      fs = fsbck

      return
      end

**********************
!*****************HERE IS THE EERAD BORN CONTRIBUTIONS

      subroutine eerad_born(p1,p2,p3,p4,q,emtx,qph)
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision qph(40,0:3)
      double precision k(0:3)
c      real*16 qp1(0:3),qp2(0:3),qp3(0:3),qp4(0:3),qq(0:3)
      real*8 p1p2,p1p3,p1p4,kp1,kp2,kp3,kp4,p2p3,p2p4,p3p4
c      real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
      real*8 e2,e6
      real*8 emtx,s,s1,t,t1
      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
      double complex elmat,elmat1
      character*2 fs
      common/finalstate/fs
! using always FORM !!!!!!!!!
      common/softphotonmass/phmass
      real*8 phmass
      external dot
      real*8 dot
      integer i

      external vpolc
      double complex vpolc,alphar

      common/ialpharunning/iarun
      integer iarun

            e2 = 4.d0*pi*alpha
            e6 = e2*e2*e2
	    p1p2 = dot(p1,p2)
	    p1p3 = dot(p1,p3)
	    p1p4 = dot(p1,p4)
	    kp1  = dot(p1,q)
	    p2p3 = dot(p2,p3)
	    p2p4 = dot(p2,p4)
	    kp2  = dot(p2,q)
	    kp3  = dot(p3,q)
	    p3p4 = dot(p3,p4)
            kp4  = dot(p4,q)

            s=2d0*(ame**2+p1p2)
            s1=2d0*(ame**2+p3p4)
            t=dot(p1-p3,p1-p3)
            t1=dot(p2-p4,p2-p4)

            elmat1=(0.d0,0.d0)

            include 'eerad/born/ampTreeISRsSquared.F'
            elmat1=elmat1+elmat*vpolc(s1)*CONJG(vpolc(s1))
            include 'eerad/born/ampTreeFSRsSquared.F'
            elmat1=elmat1+elmat*vpolc(s)*CONJG(vpolc(s))
            include 'eerad/born/ampTreeINTsSquared.F'
            elmat1=elmat1+elmat*vpolc(s)*CONJG(vpolc(s1))
            include 'eerad/born/ampTreeISRtSquared.F'
            elmat1=elmat1+elmat*vpolc(t)*CONJG(vpolc(t))
            include 'eerad/born/ampTreeFSRtSquared.F'
            elmat1=elmat1+elmat*vpolc(t1)*CONJG(vpolc(t1))
            include 'eerad/born/ampTreeINTtSquared.F'
            elmat1=elmat1+elmat*vpolc(t)*CONJG(vpolc(t1))
            include 'eerad/born/ampTreeINTts1.F'
            elmat1=elmat1+elmat*vpolc(s1)*CONJG(vpolc(t))
            include 'eerad/born/ampTreeINTts2.F'
            elmat1=elmat1+elmat*vpolc(s1)*CONJG(vpolc(t1))
            include 'eerad/born/ampTreeINTts3.F'
            elmat1=elmat1+elmat*vpolc(s)*CONJG(vpolc(t))
            include 'eerad/born/ampTreeINTts4.F'
            elmat1=elmat1+elmat*vpolc(s)*CONJG(vpolc(t1))

            emtx=4d0*elmat1*e6
      end



!*****************HERE ARE THE MUMURAD BORN CONTRIBUTIONS
      subroutine ISR_mur(p1,p2,p3,p4,q,emtx,qph)
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision qph(40,0:3)
      double precision k(0:3)
c      real*16 qp1(0:3),qp2(0:3),qp3(0:3),qp4(0:3),qq(0:3)
      real*8 p1p2,p1p3,p1p4,kp1,kp2,kp3,kp4,p2p3,p2p4,p3p4
c      real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
      real*8 e2,e6
      real*8 emtx,s
      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
      real*8 elmat
      character*2 fs
      common/finalstate/fs
! using always FORM !!!!!!!!!
      common/softphotonmass/phmass
      real*8 phmass
      external dot
      real*8 dot
      integer i

      external vpolc
      double complex vpolc,alphar

      common/ialpharunning/iarun
      integer iarun

            e2 = 4.d0*pi*alpha
            e6 = e2*e2*e2
	    p1p2 = dot(p1,p2)
	    p1p3 = dot(p1,p3)
	    p1p4 = dot(p1,p4)
	    kp1  = dot(p1,q)
	    p2p3 = dot(p2,p3)
	    p2p4 = dot(p2,p4)
	    kp2  = dot(p2,q)
	    kp3  = dot(p3,q)
	    p3p4 = dot(p3,p4)
            kp4  = dot(p4,q)
            include 'mumurad/Born/BornISR.f'

            s=2.d0*(ammu**2+p3p4)

            emtx=4d0*elmat*e6*vpolc(s)*CONJG(vpolc(s))
      end

*****************

      subroutine FSR_mur(p1,p2,p3,p4,q,emtx,qph)
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision qph(40,0:3)
      double precision k(0:3)
c      real*16 qp1(0:3),qp2(0:3),qp3(0:3),qp4(0:3),qq(0:3)
      real*8 p1p2,p1p3,p1p4,kp1,kp2,kp3,kp4,p2p3,p2p4,p3p4
c      real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
      real*8 e2,e6
      real*8 emtx,s
      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
      real*8 elmat
      character*2 fs
      common/finalstate/fs
! using always FORM !!!!!!!!!
      common/softphotonmass/phmass
      real*8 phmass
      external dot
      real*8 dot
      integer i

      external vpolc
      double complex vpolc,alphar

      common/ialpharunning/iarun
      integer iarun
            e2 = 4.d0*pi*alpha
            e6 = e2*e2*e2
	    p1p2 = dot(p1,p2)
	    p1p3 = dot(p1,p3)
	    p1p4 = dot(p1,p4)
	    kp1  = dot(p1,q)
	    p2p3 = dot(p2,p3)
	    p2p4 = dot(p2,p4)
	    kp2  = dot(p2,q)
	    kp3  = dot(p3,q)
	    p3p4 = dot(p3,p4)
            kp4  = dot(p4,q)
            include 'mumurad/Born/BornFSR.F'

            s=2.d0*(ame**2+p1p2)

            emtx=4d0*elmat*e6*vpolc(s)*CONJG(vpolc(s))
      end
******************************
      subroutine IFI_mur(p1,p2,p3,p4,q,emtx,qph)
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision qph(40,0:3)
      double precision k(0:3)
c      real*16 qp1(0:3),qp2(0:3),qp3(0:3),qp4(0:3),qq(0:3)
      real*8 p1p2,p1p3,p1p4,kp1,kp2,kp3,kp4,p2p3,p2p4,p3p4
c      real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
      real*8 e2,e6
      real*8 emtx,s1,s
      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
      real*8 elmat
      character*2 fs
      common/finalstate/fs
! using always FORM !!!!!!!!!
      common/softphotonmass/phmass
      real*8 phmass
      external dot
      real*8 dot
      integer i

      external vpolc
      double complex vpolc,alphar

      common/ialpharunning/iarun
      integer iarun

            e2 = 4.d0*pi*alpha
            e6 = e2*e2*e2
	    p1p2 = dot(p1,p2)
	    p1p3 = dot(p1,p3)
	    p1p4 = dot(p1,p4)
	    kp1  = dot(p1,q)
	    p2p3 = dot(p2,p3)
	    p2p4 = dot(p2,p4)
	    kp2  = dot(p2,q)
	    kp3  = dot(p3,q)
	    p3p4 = dot(p3,p4)
            kp4  = dot(p4,q)
            include 'mumurad/Born/BornIFI.F'

            s=2.d0*(ame**2+p1p2)
            s1=2.d0*(ammu**2+p3p4)
            emtx=4d0*elmat*e6*vpolc(s)*CONJG(vpolc(s1))
      end

*****************


!*****************HERE ARE THE PION RADIATIVE BORN CONTRIBUTIONS
      subroutine ISR_pr(p1,p2,p3,p4,q,emtx,qph)
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision qph(40,0:3)
      real*8 p1p2,p1p3,p1p4,kp1,kp2,kp3,kp4,p2p3,p2p4,p3p4,s
      real*8 e2,e6
      real*8 emtx
      double precision ame,amu,ampi,convfac,alpha,pi
      common/parameters/ame,amu,convfac,alpha,pi
      common/hparam/ampi
      real*8 elmat

      external dot
      real*8 dot

      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi
      double complex ffpiis
      double precision ffpiis2

      external vpolc
      double complex vpolc,alphar
      common/ialpharunning/iarun
      integer iarun

            e2 = 4.d0*pi*alpha
            e6 = e2*e2*e2
	    p1p2 = dot(p1,p2)
	    p1p3 = dot(p1,p3)
	    p1p4 = dot(p1,p4)
	    kp1  = dot(p1,q)
	    p2p3 = dot(p2,p3)
	    p2p4 = dot(p2,p4)
	    kp2  = dot(p2,q)
	    kp3  = dot(p3,q)
	    p3p4 = dot(p3,p4)
            kp4  = dot(p4,q)


            ffpiis=(1.d0,0.d0)
            ffpiis2=1.d0
            if(iffpi.gt.0) then
               call ffpi(2.d0*(ampi**2+p3p4),ffpiis,ffpiis2)
	    endif

	    include 'pipirad/Born/ampBornISR.F'
	    s=2.d0*(ampi**2+p3p4)
	    ffpiis=ffpiis*vpolc(s)

            emtx=4d0*elmat*e6*ffpiis*CONJG(ffpiis)

      end

*****************

      subroutine FSR_pr(p1,p2,p3,p4,q,emtx,qph)
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision qph(40,0:3)
      real*8 p1p2,p1p3,p1p4,kp1,kp2,kp3,kp4,p2p3,p2p4,p3p4,s
      real*8 e2,e6
      real*8 emtx
      double precision ame,amu,ampi,convfac,alpha,pi
      common/parameters/ame,amu,convfac,alpha,pi
      common/hparam/ampi
      real*8 elmat
      character*2 fs
      common/finalstate/fs

      external dot
      real*8 dot

      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi
      double complex ffpifs
      double precision ffpifs2

      external vpolc
      double complex vpolc,alphar
      common/ialpharunning/iarun
      integer iarun

            e2 = 4.d0*pi*alpha
            e6 = e2*e2*e2
	    p1p2 = dot(p1,p2)
	    p1p3 = dot(p1,p3)
	    p1p4 = dot(p1,p4)
	    kp1  = dot(p1,q)
	    p2p3 = dot(p2,p3)
	    p2p4 = dot(p2,p4)
	    kp2  = dot(p2,q)
	    kp3  = dot(p3,q)
	    p3p4 = dot(p3,p4)
            kp4  = dot(p4,q)

            ffpifs=(1.d0,0.d0)
            ffpifs2=1.d0
            if(iffpi.gt.0) then
               call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2)
	    endif
            include 'pipirad/Born/ampBornFSR.F'

	    s=2.d0*(ame**2+p1p2)
	    ffpifs=ffpifs*vpolc(s)

            emtx=4d0*elmat*e6*ffpifs*CONJG(ffpifs)
      end

******************************
      subroutine IFI_pr(p1,p2,p3,p4,q,emtx,qph)
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision qph(40,0:3)
      real*8 p1p2,p1p3,p1p4,kp1,kp2,kp3,kp4,p2p3,p2p4,p3p4,s,s1
      real*8 e2,e6
      double complex emtx
      double precision ame,amu,ampi,convfac,alpha,pi
      common/parameters/ame,amu,convfac,alpha,pi
      common/hparam/ampi
      double complex elmat
      character*2 fs
      common/finalstate/fs

      external dot
      real*8 dot

      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi
      double complex ffpifs
      double precision ffpifs2
      double complex ffpiis
      double precision ffpiis2

      external vpolc
      double complex vpolc,alphar
      common/ialpharunning/iarun
      integer iarun

            e2 = 4.d0*pi*alpha
            e6 = e2*e2*e2
	    p1p2 = dot(p1,p2)
	    p1p3 = dot(p1,p3)
	    p1p4 = dot(p1,p4)
	    kp1  = dot(p1,q)
	    p2p3 = dot(p2,p3)
	    p2p4 = dot(p2,p4)
	    kp2  = dot(p2,q)
	    kp3  = dot(p3,q)
	    p3p4 = dot(p3,p4)
            kp4  = dot(p4,q)

            ffpifs=(1.d0,0.d0)
            ffpifs2=1.d0

            ffpiis=(1.d0,0.d0)
            ffpiis2=1.d0
            if(iffpi.gt.0) then
               call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2)
               call ffpi(2.d0*(ampi**2+p3p4),ffpiis,ffpiis2)
	    endif
            include 'pipirad/Born/ampBornIFI.F'
            ffpifs=ffpifs*vpolc(2.d0*(ame**2+p1p2))
            ffpiis=ffpiis*vpolc(2.d0*(ampi**2+p3p4))


            emtx=4d0*elmat*e6*ffpiis*CONJG(ffpifs)
!            print*,ffpiis
!!            print*,ffpifs

      end

*--------------------------------------*      
*     CORRECTIONS FROM FASTERD         *
*--------------------------------------*

      function ffphi(s,q2)
      implicit none
      double precision q2,s
      double complex ffphi,im
      double precision mk,msig,mf,gf,gsig,g1,g2
      double precision asig,aphi
      double precision ame,amu,ampi,convfac,alpha,pi
      common/parameters/ame,amu,convfac,alpha,pi
      common/hparam/ampi
      double precision mr,gr,mphi,gphi
      common/resparams/mr,gr,mphi,gphi
      external I
      double complex I
      common/indexphires/index,indexbr,indexdres
      integer index,indexbr,indexdres
      double complex ff0_II
      real*8 qq,aa,bb,amb,mm_ph,gg_ph,aa_phi
      complex*16 ff_f0,gg_f0,f1,f2
      real*8 mKp
      real*8 s1
      double complex fphikkkloe
      im=(0.d0,1.d0)
      mk=0.493677d0
      if(index.eq.1) then
         mf = 1.00291061d0
         gf = 0.108171405d0
         g1 = 59.3919713647d0
c     phi_f0(600)-pipi
         msig = 0.519709326d0
         gsig = 0.319204396d0
         g2 = 11.2156718914d0
c     phases : 1. rho and (f0 + f0(600)), 2. f0 and f0(600) 
         asig = 98.4814123d0/180.d0*pi

!!!!  From https://arxiv.org/pdf/hep-ph/0412239

         aphi=75.d0/180.d0*pi*sqrt(q2-4.d0*ampi**2)
*------KK Model      

         ffphi=exp(im*aphi)/2.d0/pi**2/mk**2*(g1/(mf**2-q2-im*mf*gf)
     >        +exp(im*asig)*g2/(msig**2-q2-im*msig*gsig))
     >        *I(mphi**2/mk**2,q2/mk**2)
         return
*------NO structure model
      elseif(index.eq.2) then
         mf =  0.98385936d0
         gf =  0.0354330317d0
         g1 =  2.83702251679d0
c     phi_f0(600)-pipi
         msig =  0.588107848d0
         gsig =  1.65341791d0
         g2 =  4.01504359106d0
c     phases : 1. rho and (f0 + f0(600)), 2. f0 and f0(600) 
         asig =  101.144768d0/180.d0*pi
         aphi = 0.d0
         ffphi = exp(im*aphi)*(g1/(mf**2-q2-im*mf*gf)
     >        +exp(im*asig)*g2/(msig**2-q2-im*msig*gsig))
!     https://arxiv.org/pdf/hep-ph/0512047  PARAMETERISATION         
      elseif(index.eq.3) then
         ffphi = -fphikkkloe(s,q2)
      endif
      return
      end

*----------MODEL OF https://arxiv.org/pdf/hep-ph/0512047 *---------      
      function  fphikkkloe(s,q2)
      implicit none
      double complex fphikkkloe
      double precision q2,s
      common/kloematparams/mpp,mp0,mkp,mk0,meta,metap,gsigpp,gsigkk,
     >     gf0pp,gf0kk,mf0,msig,b0p,b1p,b2p,lambdap,m2kk,m0kk,lambdakk
     >	 ,cf0sigma,gphikk
      double precision b0p,b1p,b2p,lambdap,m0kk,m2kk,lambdakk
      double precision mpp,mp0,mkp,mk0,meta,metap,pi,mf0,msig
      double precision gf0pp,gf0kk,gsigkk,gsigpp,gphikk
      double complex Tsigmaf0,im,I
      external Tsigmaf0,deltaB,I    
      double precision x1,x2,deltaB,fact,fk,cf0sigma
      
     
      
      im = (0.d0,1.d0)
!!! Eq. 29 b of https://arxiv.org/abs/0901.4440  
      pi = 4.d0*atan(1.d0)
      x1 = s/(mkp)**2  
      x2 = (q2)/mkp**2
!      print*,x1,x2
      fphikkkloe = gphikk/2.d0/mkp**2/pi**2*Tsigmaf0(q2)*I(x1,x2)
      
!!! 	Additional phase coming from pp
      fphikkkloe = fphikkkloe*exp(im*deltaB(q2))

!!!     Additional factor due to kk, if q2.le.4*mkp**2 the factor is real and with module as given, above is 1 and there is just a phase. THe overall phase (w.r.t. the F_pi) is therefore factorized out

      if(q2.le.4.d0*mkp**2) then
        fk = atan((m2kk**2-q2)/m0kk**2)/lambdakk       
      	fact=dsqrt(abs( 1.d0 - fk*sqrt(4.d0*mkp**2-q2) ) /
     >        ( 1.d0 + fk*sqrt(4.*mkp**2-q2) ))
      else
      	fact=1.d0
      endif
      
      fphikkkloe = fphikkkloe*fact
      
!!!! 	Here take into account the phi propagator and phig coupling
!      fphikkkloe = +gphig/(s-mphi**2+im*mphi*gphi)*fphikkkloe      
      return
      end
      
      
      function Tsigmaf0(q2)   !!!!Checked against PH.
c	Evaluating Eq. 3 (with inverse matrix) of https://arxiv.org/pdf/hep-ph/0512047
      implicit none
      double complex Tsigmaf0,PIsigma,Pif0,Pif0sigma,det,prop
      double precision q2
      double complex mat(2,2),matinv(2,2),coup1(2),coup2(2)
      integer i,j
      common/kloematparams/mpp,mp0,mkp,mk0,meta,metap,gsigpp,gsigkk,
     >     gf0pp,gf0kk,mf0,msig,b0p,b1p,b2p,lambdap,m2kk,m0kk,lambdakk
     >	 ,cf0sigma,gphikk
      double precision b0p,b1p,b2p,lambdap,m2kk,m0kk,lambdakk
      double precision mpp,mp0,mkp,mk0,meta,metap,pi,mf0,msig
      double precision gf0pp,gf0kk,gsigkk,gsigpp,cf0sigma,gphikk
      
      
      Tsigmaf0 = (0.d0,0.d0)
      coup1(1) = gf0kk
      coup1(2) = gsigkk
      
      coup2(1) = gf0pp
      coup2(2) = gsigpp
      
      mat(1,1) = prop(q2,1)
      mat(1,2) = -Pif0sigma(q2)
      mat(2,1) = -Pif0sigma(q2)
      mat(2,2) = prop(q2,2)
      
      det = mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
      
      matinv(1,1) =  mat(2,2)/det
      matinv(1,2) = -mat(1,2)/det
      matinv(2,1) = -mat(2,1)/det
      matinv(2,2) =  mat(1,1)/det
      
      do i=1,2
         do j=1,2
            Tsigmaf0 = Tsigmaf0 + coup1(i)*matinv(i,j)*coup2(j)
         enddo
      enddo
      
      return
      end
      
      

      function PIf0sigma(m2)    !!!!TESTED
!     Function Pi after Eq.5 for f0sigma from  https://arxiv.org/pdf/hep-ph/0512047
      implicit none
      
      double complex PIf0sigma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      common/kloematparams/mpp,mp0,mkp,mk0,meta,metap,gsigpp,gsigkk,
     >     gf0pp,gf0kk,mf0,msig,b0p,b1p,b2p,lambdap,m2kk,m0kk,lambdakk
     >	 ,cf0sigma,gphikk
!    
      double precision b0p,b1p,b2p,lambdap,m2kk,m0kk,lambdakk
      double precision mpp,mp0,mkp,mk0,meta,metap,pi,mf0,msig,gphikk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision gsigmaab(7),gf0ab(7),mass(7,2)
      double precision m2, gf0pp, gf0kk, gsigpp, gsigkk, Cf0sigma
      double precision arg,mp,mm,mp2,mm2,mdifp2,mdifm2
      double complex im,fterm,rho_ab
      integer i
C     
C========================================================================
C     
      PIf0sigma=(0.d0,0.d0)
      im = (0.d0,1.d0)	
      pi=4.d0*atan(1.d0)
!!!!  All to be defined in common  also c_sifma    
      mass(1,1) = mpp
      mass(1,2) = mass(1,1)
      mass(2,1) = mp0
      mass(2,2) = mass(2,1)
      mass(3,1) = mkp
      mass(3,2) = mass(3,1)
      mass(4,1) = mk0
      mass(4,2) = mass(4,1)
      mass(5,1) = Meta
      mass(5,2) = mass(5,1)
      mass(6,1) = Meta
      mass(6,2) = Metap
      mass(7,1) = Metap
      mass(7,2) = mass(7,1)
      
!!!!!!gsigkk should be given, while gf0kk should be gf0pipi/sqrt(2)

      gsigmaab(1)  =  gsigpp    ! sig pi+pi- coupling
      gsigmaab(2)  =  gsigpp/dsqrt(2.d0) ! sig pi0pi0 coupling
      gsigmaab(3)  =  gsigkk    ! sig K+K-   coupling
      gsigmaab(4)  =  gsigkk    ! sig KlKs   coupling
      gsigmaab(5)  =  gsigpp*dsqrt(2.d0)/3.d0 ! sig etaeta coupling
      gsigmaab(6)  =  gsigpp*dsqrt(2.d0)/3.d0 ! sig eta1eta coupling
      gsigmaab(7)  =  gsigpp*dsqrt(1.d0/18.d0) ! sig eta1eta1 coupling
      gf0ab(1)   =  gf0pp       ! f0 pi+pi- coupling
      gf0ab(2)   =  1.d0/dsqrt(2.d0)*gf0pp ! f0 pi0pi0 coupling
      gf0ab(3)   =  gf0kk       ! f0 K+K-   coupling
      gf0ab(4)   =  gf0kk       ! f0 KlKs   coupling
      gf0ab(5)   =  gf0kk*dsqrt(8.d0)/3.d0 ! f0 etaeta coupling
      gf0ab(6)   = -gf0kk*dsqrt(2.d0)/3.d0 ! f0 eta1eta 
      gf0ab(7)   = -gf0ab(5)    ! f0 eta1eta1
C     

C     
      do i=1,7
         mp = mass(i,1)+mass(i,2)
         mm = mass(i,1)-mass(i,2)
         
         mp2 = mp**2
         mm2 = mm**2 
         
         fterm= mp*mm/pi/m2*dlog(mass(i,2)/mass(i,1)) 
         rho_ab = dsqrt(abs((1.d0-mp2/m2)*(1.d0-mm2/m2)))
         
         mdifp2 = abs(m2-mp2)
         mdifm2 = abs(m2-mm2)
         
         if(m2.ge.mp2) then
            arg=abs((dsqrt(mdifm2)-dsqrt(mdifp2))
     >           /(dsqrt(mdifm2)+dsqrt(mdifp2)))
            PIf0sigma = PIf0sigma + gsigmaab(i)*gf0ab(i)/(16.d0*pi)
     >           *(fterm+rho_ab*(im+1.d0/pi*dlog(arg)))	
         elseif(m2.ge.mm2.and.m2.le.mp2) then
            
            arg = sqrt(mdifp2/mdifm2)
            PIf0sigma = PIf0sigma + gsigmaab(i)*gf0ab(i)/(16.d0*pi)
     >           *(fterm-rho_ab*(1.d0-2.d0/pi*atan(arg)))
         elseif(m2.le.mm2) then
            
            arg = abs((dsqrt(mdifm2)-dsqrt(mdifp2))
     >           /(dsqrt(mdifm2)+dsqrt(mdifp2)))
            PIf0sigma = PIf0sigma + gsigmaab(i)*gf0ab(i)/(16.d0*pi)
     >           *(fterm-1.d0/pi*rho_ab*dlog(arg))	
         endif
      enddo
      
      PIf0sigma = PIf0sigma + Cf0sigma
      
      return
      end   
C     

c***********************************************************************
      function PIf0(m2)         !!!!TESTED
!     Function Pi of Eq.6 for f0 from  https://arxiv.org/pdf/hep-ph/0512047
      implicit none
      double complex PIf0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      common/kloematparams/mpp,mp0,mkp,mk0,meta,metap,gsigpp,gsigkk,
     >     gf0pp,gf0kk,mf0,msig,b0p,b1p,b2p,lambdap,m2kk,m0kk,lambdakk
     >	 ,cf0sigma,gphikk
!     
      double precision b0p,b1p,b2p,lambdap,m2kk,m0kk,lambdakk
      double precision mpp,mp0,mkp,mk0,meta,metap,pi,mf0,msig,gphikk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision gf0ab(7),mass(7,2)
      double precision m2, gf0pp, gf0kk, gsigpp, gsigkk, Cf0sigma
      double precision arg,mp,mm,mp2,mm2,mdifp2,mdifm2
      double complex im,fterm,rho_ab
      integer i
C     
C========================================================================
C     
      PIf0=(0.d0,0.d0)
      im = (0.d0,1.d0)	
      pi = 4.d0*atan(1.d0)
      mass(1,1) = mpp
      mass(1,2) = mass(1,1)
      mass(2,1) = mp0
      mass(2,2) = mass(2,1)
      mass(3,1) = mkp
      mass(3,2) = mass(3,1)
      mass(4,1) = mk0
      mass(4,2) = mass(4,1)
      mass(5,1) = Meta
      mass(5,2) = mass(5,1)
      mass(6,1) = Meta
      mass(6,2) = Metap
      mass(7,1) = Metap
      mass(7,2) = mass(7,1)
      
!!!!!!gsigkk should be given, while gf0kk should be gf0pipi/sqrt(2)

      
      gf0ab(1)   =  gf0pp       ! f0 pi+pi- coupling
      gf0ab(2)   =  1.d0/dsqrt(2.d0)*gf0pp ! f0 pi0pi0 coupling
      gf0ab(3)   =  gf0kk       ! f0 K+K-   coupling
      gf0ab(4)   =  gf0kk       ! f0 KlKs   coupling
      gf0ab(5)   =  gf0kk*dsqrt(8.d0)/3.d0 ! f0 etaeta coupling
      gf0ab(6)   = -gf0kk*dsqrt(2.d0)/3.d0 ! f0 eta1eta 
      gf0ab(7)   = -gf0ab(5)    ! f0 eta1eta1
C     

C     
      do i=1,7
         mp = mass(i,1)+mass(i,2)
         mm = mass(i,1)-mass(i,2)
         
         mp2 = mp**2
         mm2 = mm**2 
         
         fterm= mp*mm/pi/m2*dlog(mass(i,2)/mass(i,1)) 
         rho_ab = dsqrt(abs((1.d0-mp2/m2)*(1.d0-mm2/m2)))
         
         mdifp2 = abs(m2-mp2)
         mdifm2 = abs(m2-mm2)
         
         if(m2.ge.mp2) then
            arg=abs((dsqrt(mdifm2)-dsqrt(mdifp2))
     >           /(dsqrt(mdifm2)+dsqrt(mdifp2)))
            PIf0 = PIf0 + gf0ab(i)**2/(16.d0*pi)
     >           *(fterm+rho_ab*(im+1.d0/pi*dlog(arg)))	
         elseif(m2.ge.mm2.and.m2.le.mp2) then
            
            arg = sqrt(mdifp2/mdifm2)
            PIf0 = PIf0 + gf0ab(i)**2/(16.d0*pi)
     >           *(fterm-rho_ab*(1.d0-2.d0/pi*atan(arg)))
         elseif(m2.le.mm2) then
            
            arg = abs((dsqrt(mdifm2)-dsqrt(mdifp2))
     >           /(dsqrt(mdifm2)+dsqrt(mdifp2)))
            PIf0 = PIf0 + gf0ab(i)**2/(16.d0*pi)
     >           *(fterm-1.d0/pi*rho_ab*dlog(arg))	
         endif
      enddo
      

      return
      end   
C     

c***********************************************************************
      function PIsigma(m2)      !!!!TESTED
!     Function Pi of Eq.6 for sigma from  https://arxiv.org/pdf/hep-ph/0512047
      implicit none

      double complex PIsigma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      common/kloematparams/mpp,mp0,mkp,mk0,meta,metap,gsigpp,gsigkk,
     >     gf0pp,gf0kk,mf0,msig,b0p,b1p,b2p,lambdap,m2kk,m0kk,lambdakk
     >	 ,cf0sigma,gphikk
!
      double precision b0p,b1p,b2p,lambdap,m2kk,m0kk,lambdakk
      double precision mpp,mp0,mkp,mk0,meta,metap,pi,mf0,msig,gphikk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision gsigmaab(7),mass(7,2)
      double precision m2, gf0pp, gf0kk, gsigpp, gsigkk, Cf0sigma
      double precision arg,mp,mm,mp2,mm2,mdifp2,mdifm2
      double complex im,fterm,rho_ab
      integer i
C     
C========================================================================
C     
      PIsigma=(0.d0,0.d0)
      im = (0.d0,1.d0)	
      pi=4.d0*atan(1.d0)
      
      mass(1,1) = mpp
      mass(1,2) = mass(1,1)
      mass(2,1) = mp0
      mass(2,2) = mass(2,1)
      mass(3,1) = mkp
      mass(3,2) = mass(3,1)
      mass(4,1) = mk0
      mass(4,2) = mass(4,1)
      mass(5,1) = Meta
      mass(5,2) = mass(5,1)
      mass(6,1) = Meta
      mass(6,2) = Metap
      mass(7,1) = Metap
      mass(7,2) = mass(7,1)
      
!!!!!!gsigkk should be given, while gf0kk should be gf0pipi/sqrt(2)

      gsigmaab(1)  =  gsigpp    ! sig pi+pi- coupling
      gsigmaab(2)  =  gsigpp/dsqrt(2.d0) ! sig pi0pi0 coupling
      gsigmaab(3)  =  gsigkk    ! sig K+K-   coupling
      gsigmaab(4)  =  gsigkk    ! sig KlKs   coupling
      gsigmaab(5)  =  gsigpp*dsqrt(2.d0)/3.d0 ! sig etaeta coupling
      gsigmaab(6)  =  gsigpp*dsqrt(2.d0)/3.d0 ! sig eta1eta coupling
      gsigmaab(7)  =  gsigpp*dsqrt(1.d0/18.d0) ! sig eta1eta1 coupling

C     
      do i=1,7
         mp = mass(i,1)+mass(i,2)
         mm = mass(i,1)-mass(i,2)
         
         mp2 = mp**2
         mm2 = mm**2 
         
         fterm= mp*mm/pi/m2*dlog(mass(i,2)/mass(i,1)) 
         rho_ab = dsqrt(abs((1.d0-mp2/m2)*(1.d0-mm2/m2)))
         
         mdifp2 = abs(m2-mp2)
         mdifm2 = abs(m2-mm2)
         
         if(m2.ge.mp2) then
            arg=abs((dsqrt(mdifm2)-dsqrt(mdifp2))
     >           /(dsqrt(mdifm2)+dsqrt(mdifp2)))
            PIsigma = PIsigma + gsigmaab(i)**2/(16.d0*pi)
     >           *(fterm+rho_ab*(im+1.d0/pi*dlog(arg)))	
         elseif(m2.ge.mm2.and.m2.le.mp2) then
            
            arg = sqrt(mdifp2/mdifm2)
            PIsigma = PIsigma + gsigmaab(i)**2/(16.d0*pi)
     >           *(fterm-rho_ab*(1.d0-2.d0/pi*atan(arg)))
         elseif(m2.le.mm2) then
            
            arg = abs((dsqrt(mdifm2)-dsqrt(mdifp2))
     >           /(dsqrt(mdifm2)+dsqrt(mdifp2)))
            PIsigma = PIsigma + gsigmaab(i)**2/(16.d0*pi)
     >           *(fterm-1.d0/pi*rho_ab*dlog(arg))	
         endif
      enddo
      

      return
      end   
C     

c***********************************************************************
      function prop(q2,ind)
!     Propagator of Eq.6 for f0 (ind=1) and sigma (ind=2) from  https://arxiv.org/pdf/hep-ph/0512047      
      implicit none
      integer ind
      double precision q2
      double complex prop,im
      common/kloematparams/mpp,mp0,mkp,mk0,meta,metap,gsigpp,gsigkk,
     >     gf0pp,gf0kk,mf0,msig,b0p,b1p,b2p,lambdap,m2kk,m0kk,lambdakk
     >	 ,cf0sigma,gphikk
!
      double precision b0p,b1p,b2p,lambdap,m2kk,m0kk,lambdakk
      double precision mpp,mp0,mkp,mk0,meta,metap,pi,mf0,msig,gphikk
      double precision gf0pp,gf0kk,gsigkk,gsigpp,cf0sigma
      external PIf0,Pisigma
      double complex PIf0,Pisigma
      
      if(ind.eq.1) then         !!!f0 res
         prop = mf0**2 - q2 +real(PIf0(mf0**2))-PIf0(q2)
      elseif(ind.eq.2) then     !!!sigma res
         prop = msig**2 - q2 +real(PIsigma(msig**2))-PIsigma(q2)
      endif
      
      return
      end

c***********************************************************************

      function deltaB(q2)
!     Eq.26 from  https://arxiv.org/pdf/hep-ph/0512047      
      implicit none
      double precision q2,deltaB
      common/kloematparams/mpp,mp0,mkp,mk0,meta,metap,gsigpp,gsigkk,
     >     gf0pp,gf0kk,mf0,msig,b0p,b1p,b2p,lambdap,m2kk,m0kk,lambdakk
     >	 ,cf0sigma,gphikk
      double precision b0p,b1p,b2p,lambdap,m2kk,m0kk,lambdakk
      double precision gf0kk,gf0pp,gsigpp,gsigkk,gphikk
      double precision mpp,mp0,mkp,mk0,meta,metap,pi,mf0,msig
      common/hparam/ampi
      double precision ampi,p_pi,p_p,cf0sigma
      external p_pi
      pi=4.d0*atan(1.d0)
      		
      	if(q2.ge.ampi**2*4.d0) then
      	  p_p=p_pi(q2)
	       deltaB= atan(-p_p/ampi*(b0p-b1p*p_p**2/ampi**2+
     >  b2p*p_p**4/ampi**4)/(1.d0+(2.d0*p_p)**2/Lambdap**2)) 	
	endif
      return
      end
      



      function I(a,b)
      implicit none
      double precision a,b
      double complex I
      external f1,f2,gg_f0
      double complex f1,f2,gg_f0,ff_f0
      
      I=1.d0/2.d0/(a-b)-2.d0/(a-b)**2*(f1(1.d0/b)-f1(1.d0/a))+
     >     a/(a-b)**2*(f2(1.d0/b)-f2(1.d0/a))

      
      return
      end

            
      function f1(a)
      implicit none
      double precision a
      double complex f1,im
      double precision ame,amu,ampi,convfac,alpha,pi
      common/parameters/ame,amu,convfac,alpha,pi
      
      im=(0.d0,1.d0)
      if(a.ge.(1.d0/4.d0)) f1=-(asin(1.d0/2.d0/sqrt(a)))**2
      if(a.le.(1.d0/4.d0)) then   
         f1=1.d0/4.d0*(log((1.d0+sqrt(1.d0-4.d0*a))**2/4.d0/a)-im*pi)**2
      endif
      
      return
      end
      
      
      function f2(a)
      implicit none
      double precision a,etap,etam
      double complex f2,im
      double precision ame,amu,ampi,convfac,alpha,pi
      common/parameters/ame,amu,convfac,alpha,pi
      
      im=(0.d0,1.d0)
      if(a.ge.(1.d0/4.d0)) f2=sqrt(4.d0*a-1.d0)*asin(1.d0/2.d0/sqrt(a))
      if(a.le.(1.d0/4.d0)) then
         f2=1.d0/2.d0*(log((1.d0+sqrt(1.d0-4.d0*a))**2/4.d0/a)-im*pi)
     >     *sqrt(1.d0-4.d0*a)
      endif
      return
      
      return
      end

      function D(ma,Q,l,k)
      implicit none
      double complex D
      double precision ma
      double precision Q(0:3),l(0:3),k(0:3)
      external dot
      double precision dot
      
      D= ma**2- (dot(Q,Q)+dot(l,l)+2.d0*dot(k,Q)+4.d0*dot(k,l))/4.d0
      return
      end

            function p_pi(q2)
      implicit none
      double precision q2,p_pi
      double precision ame,amu,ampi,convfac,alpha,pi
      common/parameters/ame,amu,convfac,alpha,pi
      common/hparam/ampi

      p_pi = 0.5d0*dsqrt(q2-4.d0*ampi**2)
      
      return
      end
      
      function den(q2)
      implicit none
      double complex den,im
      double precision q2,gammaeff
      double precision mr,gr,mphi,gphi
      common/resparams/mr,gr,mphi,gphi
      external p_pi
      double precision p_pi
      double precision ame,amu,ampi,convfac,alpha,pi
      common/parameters/ame,amu,convfac,alpha,pi
      common/hparam/ampi

!      print*,'q2',q2
      im  = (0.d0,1.d0)
!!!!  Taken from GS parameterization and Ph.
      if((q2-4.d0*ampi*ampi).ge.0.d0) then
         gammaeff = gr*(p_pi(q2)/p_pi(mr*mr))**3
     >        *dsqrt(mr**2/q2)
!         print*, gr,p_pi(q2)/p_pi(mr*mr),mr**2/q2
      else
         gammaeff=0.d0
      endif
      
      den = mr**2-q2-im*mr*gammaeff      

!      print*,'gammaeff',gammaeff
      return
      end

*-----------------------------*
*     DOUBLE RES DECAYS       *
*-----------------------------*      
      
      subroutine ff_doublephi(q2,ff_two)
      implicit none
      double precision q2
      double complex ff_two,im
      double precision ame,amu,ampi,convfac,alpha,pi
      common/parameters/ame,amu,convfac,alpha,pi
      common/hparam/ampi
      double precision mr,gr,mphi,gphi
      common/resparams/mr,gr,mphi,gphi
      common/doubleresparams/g_rpig,g_phirp,g_phiee,beta_wphi,
     >     PirhoVMD,Cmod_rpi,Cpha_rpi,delta_rho
      double precision g_rpig,g_phirp,g_phiee,beta_wphi,
     >     PirhoVMD,Cmod_rpi,Cpha_rpi,delta_rho	
      double complex c_rpi,com_rphi,f_phi
      
      im = (0.d0,1.d0)
      f_phi   = sqrt(3.d0/4.d0*g_phiee/(alpha**2*pi)/mphi)
      
!     Eq. 33 of https://arxiv.org/pdf/0901.4440
      
      ff_two  = g_phirp*f_phi*q2

      ff_two  = ff_two/(mphi**2-q2-im*mphi*gphi)
      ff_two  = ff_two*exp(im*beta_wphi)*PirhoVMD
!     Including also  γ∗ → ω′ → ρπ

      c_rpi    = cmod_rpi*exp(im*cpha_rpi)
      com_rphi = c_rpi/g_rpig/dsqrt(4.*pi*alpha)
      ff_two  = ff_two + com_rphi 
!     Additional phase
      
      ff_two = ff_two*exp(im*delta_rho)

      ff_two = ff_two/dsqrt(4.d0*pi*alpha)/4.d0*g_rpig
c$$$
c$$$      print*,'inverso del prop phi costante',
c$$$     >     1.d0/ff_two/(mphi**2-q2-im*mphi*gphi)
      

      return
      end

      
      subroutine FSR_pr_phi(p1,p2,p3,p4,q,emtx,qph)
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision qph(40,0:3),l(0:3),qb(0:3),qs(0:3)
      real*8 p1p2,p1p3,p1p4,kp1,kp2,kp3,kp4,p2p3,p2p4,p3p4,s
      real*8 e2,e6
      real*8 emtx
      double complex f1,f2,f3,cf1,cf2,cf3,im
      double complex df1,df2,df3
      double precision fv,gv,fa,fp,ma,gphig
      common/brehmstralparams/fv,gv,fa,fp,ma,gphig
      external ffphi
      double complex ffphi
      external D
      double complex D
      double precision ame,amu,ampi,convfac,alpha,pi
      common/parameters/ame,amu,convfac,alpha,pi
      common/hparam/ampi
      double precision mr,gr,mphi,gphi,geff
      common/resparams/mr,gr,mphi,gphi
      double complex ff_dphi,dp,dm
      external den
      double complex den
      double complex elmat
      character*2 fs
      common/finalstate/fs

      external dot
      real*8 dot

      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi
      double complex ffpifs
      double precision ffpifs2

      external vpolc
      double complex vpolc,alphar
      common/ialpharunning/iarun
      integer iarun
      common/indexphires/indexphires,indexbr,indexdres
      integer indexphires,indexbr,indexdres
      external p_pi
      double precision p_pi

!     index=1
      im=(0.d0,1.d0)      
      e2 = 4.d0*pi*alpha
      e6 = e2*e2*e2	   
      p1p2 = dot(p1,p2)
      p1p3 = dot(p1,p3)
      p1p4 = dot(p1,p4)
      kp1  = dot(p1,q)
      p2p3 = dot(p2,p3)
      p2p4 = dot(p2,p4) 
      kp2  = dot(p2,q)
      kp3  = dot(p3,q)
      p3p4 = dot(p3,p4)
      kp4  = dot(p4,q) 
      qb=p1+p2
      qs=p3+p4
      l=p4-p3

      ffpifs=(1.d0,0.d0)
      ffpifs2=1.d0
      if(iffpi.gt.0) then
         call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2)          
      endif
      
      
!!!   Pure FF	    
      f1=2.d0*dot(q,qb)*ffpifs/(dot(q,qb)**2-dot(q,l)**2)
      f2=-2.d0*ffpifs/(dot(q,qb)**2-dot(q,l)**2)
      f3=0.d0

     
      df1 = (0.d0,0.d0)
      df2 = (0.d0,0.d0)
      df3 = (0.d0,0.d0)
      geff=mr*gv**2/48.d0/pi/fp**4*dot(qb,qb)*
     >     (1.d0-4.d0*ampi**2/dot(qb,qb))**(3.d0/2.d0)

!     Brehmstralung corrections RPT
      df1=(fv**2-2.d0*fv*gv)/fp**2*(1.d0/mr**2+1.d0/
     >     (mr**2-dot(qb,qb)-im*mr*geff))
      df1=df1-fa**2/fp**2/ma**2*(2.d0+1.d0/D(ma,qb,l,q)
     >     /D(ma,qb,-l,q)*(dot(q,l)**2+1.d0/8.d0*((dot(qb,qb)+
     >     dot(q,qb))*(4.d0*ma**2-(dot(qb,qb)+dot(l,l)
     >     +2.d0*dot(qb,q))))))
            
      df2=-fa**2/fp**2/ma**2*(4.d0*ma**2-(dot(qb,qb)+dot(l,l)
     >     +2.d0*dot(qb,q)))
      df2=df2/D(ma,qb,l,q)/D(ma,qb,-l,q)/8.d0
      
      df3=fa**2/fp**2/ma**2*dot(q,l)/D(ma,qb,l,q)
     >     /D(ma,qb,-l,q)/2.d0

      if(indexbr.ne.0) then
         f1=f1+df1
         f2=f2+df2
         f3=f3+df3
      endif
      
!     Phi direct decays  e+ + e− → φ → (f0; f0 + σ)γ → π + π + γ
      df1=  gphig/(-dot(qb,qb)+mphi**2-im*mphi*gphi)
     >     *ffphi(dot(qb,qb),dot(qs,qs))
      df2=0.d0
      df3=0.d0
      
      if(indexphires.ne.0) then
         f1=f1+df1
         f2=f2+df2
         f3=f3+df3
      endif
!     Double resonant decay e+ + e− → (φ) → ρπ → π + π + γ 
      dm = 1.d0/den(dot(p3+q,p3+q))
      dp = 1.d0/den(dot(p4+q,p4+q))
      call ff_doublephi(dot(qb,qb),ff_dphi)
      
      df1 = + ff_dphi*((dot(q,qb)+dot(l,l))*(dp+dm)+
     >     2.d0*dot(q,l)*(dp-dm))
      df2 = - ff_dphi*(dp+dm)
      df3 = + ff_dphi*(dp-dm)

      if(indexdres.ne.0) then
         f1=f1+df1
         f2=f2+df2
         f3=f3+df3
      endif
      
      cf1=Conjg(f1)
      cf2=Conjg(f2)
      cf3=Conjg(f3)
      
      include 'pipirad/Born/ampBornFSRphi.F'
      
!      print*,'elmat',elmat
      emtx=4d0*elmat*e6
      end      
******************************   

      subroutine IFI_pr_phi(p1,p2,p3,p4,q,emtx,qph)
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision qph(40,0:3),l(0:3),qb(0:3),qs(0:3)
      real*8 p1p2,p1p3,p1p4,kp1,kp2,kp3,kp4,p2p3,p2p4,p3p4,s,s1
      real*8 e2,e6
      double complex emtx     
      double complex f1,f2,f3,cf1,cf2,cf3,im
      double complex df1,df2,df3
      external ffphi
      double complex ffphi
      common/indexphires/indexphires,indexbr,indexdres
      integer indexphires,indexbr,indexdres
      external D
      double complex D
      double precision ame,amu,ampi,convfac,alpha,pi
      common/parameters/ame,amu,convfac,alpha,pi
      common/hparam/ampi
      double precision mr,gr,mphi,gphi,geff
      common/resparams/mr,gr,mphi,gphi
      double precision fv,gv,fa,fp,ma,gphig
      common/brehmstralparams/fv,gv,fa,fp,ma,gphig
      double complex ff_dphi,dp,dm
      external den
      double complex den
      double complex elmat
      character*2 fs
      common/finalstate/fs
      
      external dot
      real*8 dot

      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi
      double complex ffpifs
      double precision ffpifs2
      double complex ffpiis
      double precision ffpiis2  

      external vpolc
      double complex vpolc,alphar
      common/ialpharunning/iarun
      integer iarun
      external p_pi
      double precision p_pi

!     index=0
      im=(0.d0,1.d0)      
      e2 = 4.d0*pi*alpha
      e6 = e2*e2*e2	   
      p1p2 = dot(p1,p2)
      p1p3 = dot(p1,p3)
      p1p4 = dot(p1,p4)
      kp1  = dot(p1,q)
      p2p3 = dot(p2,p3)
      p2p4 = dot(p2,p4) 
      kp2  = dot(p2,q)
      kp3  = dot(p3,q)
      p3p4 = dot(p3,p4)
      kp4  = dot(p4,q) 
      qb=p1+p2
      qs=p3+p4
      l=p4-p3
      
      ffpifs=(1.d0,0.d0)
      ffpifs2=1.d0
      
      ffpiis=(1.d0,0.d0)
      ffpiis2=1.d0
      if(iffpi.gt.0) then
         call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2) 
         call ffpi(2.d0*(ampi**2+p3p4),ffpiis,ffpiis2)                 
      endif
            

!!!   Pure FF	    
      f1=2.d0*dot(q,qb)*ffpifs/(dot(q,qb)**2-dot(q,l)**2)
      f2=-2.d0*ffpifs/(dot(q,qb)**2-dot(q,l)**2)
      f3=0.d0
      
      df1 = (0.d0,0.d0)
      df2 = (0.d0,0.d0)
      df3 = (0.d0,0.d0)
      
      geff=mr*gv**2/48.d0/pi/fp**4*dot(qb,qb)*
     >           (1.d0-4.d0*ampi**2/dot(qb,qb))**(3.d0/2.d0)

!!!   Brehmstralung corrections
      df1=(fv**2-2.d0*fv*gv)/fp**2*(1.d0/mr**2+1.d0/
     >     (mr**2-dot(qb,qb)-im*mr*geff))
      df1=df1-fa**2/fp**2/ma**2*(2.d0+1.d0/D(ma,qb,l,q)
     >     /D(ma,qb,-l,q)*(dot(q,l)**2+1.d0/8.d0*((dot(qb,qb)+
     >     dot(q,qb))*(4.d0*ma**2-(dot(qb,qb)+dot(l,l)
     >     +2.d0*dot(qb,q))))))
      df2=-fa**2/fp**2/ma**2*(4.d0*ma**2-(dot(qb,qb)+dot(l,l)
     >     +2.d0*dot(qb,q)))
      df2=df2/D(ma,qb,l,q)/D(ma,qb,-l,q)/8.d0
      
      df3=fa**2/fp**2/ma**2*dot(q,l)/D(ma,qb,l,q)
     >     /D(ma,qb,-l,q)/2.d0

      if(indexbr.ne.0) then
         f1=f1+df1
         f2=f2+df2
         f3=f3+df3
      endif

!     Phi direct decays  e+ + e− → φ → (f0; f0 + σ)γ → π + π + γ
      df1= gphig/(-dot(qb,qb)+mphi**2-im*mphi*gphi)
     >     *ffphi(dot(qb,qb),dot(qs,qs))
      df2= 0.d0
      df3= 0.d0

      if(indexphires.ne.0) then
         f1=f1+df1
         f2=f2+df2
         f3=f3+df3
      endif

!     Double resonant decay e+ + e− → (φ) → ρπ → π + π + γ 
      dm = 1.d0/den(dot(p3+q,p3+q))
      dp = 1.d0/den(dot(p4+q,p4+q))
      call ff_doublephi(dot(qb,qb),ff_dphi)
      
      df1 = + ff_dphi*((dot(q,qb)+dot(l,l))*(dp+dm)+
     >     2.d0*dot(q,l)*(dp-dm))
      df2 = - ff_dphi*(dp+dm)
      df3 = + ff_dphi*(dp-dm)

      if(indexdres.ne.0) then
         f1=f1+df1
         f2=f2+df2
         f3=f3+df3
      endif
!      print*,f1
      
      cf1=Conjg(f1)
      cf2=Conjg(f2)
      cf3=Conjg(f3)
      
      include 'pipirad/Born/ampBornIFIphi.F'
      
      
      emtx=4d0*elmat*e6*ffpiis
!      emtx=0.d0
      end


*-------------------------------
      subroutine rhotopig(Delta_f1rhoch,Delta_f2rhoch,
     $  Delta_f3rhoch,s1, rkpi1,rkpi2,qq)
C*****************************************************************************
c******* DOUBLE RESONANCE PART:\PHI\to\RHO^\pm\PI^mp\to\PI^+\PI^-\gamma  *****
c*****************************************************************************
	implicit none
	COMPLEX*16 Delta_f1rhoch,Delta_f2rhoch,
     $  Delta_f3rhoch
c	real*8 hv,fpi
c	common/const7/fpi
c	common/couple2/hv
        COMMON/const/s,Pi,mrho,mpi,grhoee,gammarho,alpha
        DOUBLE PRECISION s,Pi,mrho,mpi,grhoee,gammarho,alpha
	real*8 q2_Q,q2_pi,l2
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        real*8 rkpi1,rkpi2,s1,qq
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	real*8 ell,g_rho,g_phi
c        common/deltaf/q2Q,q2_pi,l2
	complex*16 cp,cm,c_sum,c_dif,ff_twores,cpom,cmom,
     $    c_omsum,c_omdif,ff_omres
	real*8 g_rhopig,g_phirpi,Prhores,beta_bro,beta_wphi,
     $            c_modrpi,c_pharpi,c_modwpi,c_phawpi
	common/vmdcom/g_rhopig,g_phirpi,Prhores,beta_bro,beta_wphi,
     $        c_modrpi,c_pharpi,c_modwpi,c_phawpi
	real*8 be,al,momega,gomega,mrhol,grhol,arga 
	COMMON/const2/be,al,momega,gomega,mrhol,grhol,arga
	double precision factfs,factcross_fs
c	real factfs,factcross_fs
	COMMON/constfs/ factfs,factcross_fs
c        REAL cut,gcut,ssin,emnsin,qmnsin,gmnsin,gmxsin
c        INTEGER nkind,nf0kind,ifif,ncrf,nr0kind,nfpkind,
c     *        nfskind,nrespikind
	REAL fppar,f0par,f0fase,vmdpar
c        COMMON/UTMP/cut(4),gcut(2),ssin,emnsin,qmnsin,
c     *             gmnsin,gmxsin,nfskind,nkind(2),nf0kind,ifif,ncrf,
c     $             nr0kind(5),nfpkind,nrespikind,fppar(13),f0par(9),
c     $             f0fase(7),vmdpar(9)

c$$$        g_rho=4.d0*hv/3.d0*dsqrt(2.d0/fpi)
        pi= 4.d0*atan(1.d0)
        alpha = 1.d0/137.03599908d0
	ell = dsqrt(4.d0*pi*alpha) ! ell=0.303
        mpi = 139.57039d-3

        print*,'ell',ell,pi,alpha
        
c$$$	g_rho=0.295d0/ell	!WE USE RELATION BETWEEN THE PARAMETERS	
c$$$	g_phi=0.811d0    ! from Gino nad Achasov works, 

        g_rhopig =0.295d0
        g_phirpi =0.811d0
	g_rho=g_rhopig/ell	!WE USE RELATION BETWEEN THE PARAMETERS	
	g_phi=g_phirpi     ! from Gino nad Achasov works, 

        print*,g_rho
	print*,g_phirpi
        
        mrho =0.77526d0 
        gammarho =0.1478d0
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        q2_Q=rkpi1+rkpi2
!        print*,'qdotqbig',q2_Q
        l2=4.d0*mpi**2-qq
!        print*,'l2',l2
        q2_pi=rkpi1-rkpi2
!        print*,'l dot qs',q2_pi
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	

	call comp(mrho,s1,q2_Q,l2,q2_pi,cp,cm)
!        print*,'cp cm',cp,cm
c	print*,'OK rhopi 1.1 cp cm',cp,cm
	c_sum = 1.d0/cp+1.d0/cm
	c_dif = 1.d0/cp-1.d0/cm

!        print*,'p e m',c_sum,c_dif
c$$$	call comp(momega,cpom,cmom)
c$$$	c_omsum = 1.d0/cpom+1.d0/cmom
c$$$	c_omdif = 1.d0/cpom-1.d0/cmom


	
	ff_twores=0.d0
c	print*,nrespikind
c	if(nrespikind.eq.1.or.nrespikind.eq.3)then

	call F_twores(s1,ff_twores)
!        print*,'ff phpokh',ff_twores*g_rho/4.d0
c	endif

	 Delta_f1rhoch = g_rho*ff_twores/4.d0*(
     $   (q2_Q+l2)*c_sum+2.d0*q2_pi*c_dif )
!         print*,'prova 1',(q2_Q+l2)*c_sum,2.d0*q2_pi*c_dif
	Delta_f2rhoch = -g_rho*ff_twores/4.d0*c_sum ! original sign

	Delta_f3rhoch = g_rho*ff_twores/4.d0*c_dif


!	print*, Delta_f1rhoch,Delta_f3rhoch,Delta_f3rhoch
c
c        Delta_f1rhoch = 0.d0
c        Delta_f2rhoch = 0.d0
c        Delta_f3rhoch = 0.d0
c
	return
	end
c------------------------------------------------------------
	subroutine comp(mm,s1,q2_Q,l2,q2_pi,cp,cm)
	implicit none
	COMPLEX*16 cp,cm
	real*8 q2_Q,q2_pi,l2,s1
	real*8 mm
c	common/const18/fa,ma
c	common/deltaf/q2_Q,q2_pi,l2	
        COMMON/const/s,Pi,mrho,mpi,grhoee,gammarho,alpha
        DOUBLE PRECISION s,Pi,mrho,mpi,grhoee,gammarho,alpha
	double precision smax,gam,smin,gammax,gammin
	real*8 be,al,momega,gomega,mrhol,grhol,arga 
	COMMON/const2/be,al,momega,gomega,mrhol,grhol,arga
	COMPLEX*16 i

	i=(0.D0,1.D0)
        
	smax = (s1+l2+2.d0*q2_Q+4.d0*q2_pi)/4.d0
	smin = (s1+l2+2.d0*q2_Q-4.d0*q2_pi)/4.d0
!        print*,'q2',smax
!        print*,'q2',smin
c$$$  
	if(mm.eq.mrho) then
	   gammax=gam(smax)
	   gammin=gam(smin)
	else if(mm.eq.momega) then	   
	gammax = gomega
	gammin=gomega
	else
	   print*, 'wrong choise VMD contribution'
	endif
c$$$          
!        print*,'gammamax',gammax
!        print*,'gammamin',gammin
        print*,'mm',mm
        cp = mm**2-smax -i*mm*gammax ! gam(smax)
 	cm = mm**2-smin -i*mm*gammin ! gam(smin)

c	print*,gammax,smax,4.d0*mpi**2,mrho 
c        print*,gammin,smin,4.d0*mpi**2,mrho
c$$$	cp = mm**2-smax -i*mm*gam(smax)
c$$$c	print*,'cp=',cp
c$$$	smin = (s+l2+2.d0*q2_Q-4.d0*q2_pi)/4.d0
c$$$c	print*,'smin=',smin
c$$$	cm = mm**2-smin -i*mm*gam(smin)
c$$$c	print*,'cm=',cm
	return 
	end

c-----------------------------------------
	DOUBLE PRECISION FUNCTION GAM(x)
	implicit none
	real*8 ps,pm,theta,x,width
        COMMON/const/s,Pi,mrho,mpi,grhoee,gammarho,alpha
        DOUBLE PRECISION s,Pi,mrho,mpi,grhoee,gammarho,alpha
        mpi=139.57039d-3
	if(THETA(x - 4.d0*mpi**2).eq.0.d0) then
	   gam = 0.d0
	else	   
c	ps = dsqrt(x/4.d0 - mpi**2)
c	pm = dsqrt(mrho**2/4.d0 - mpi**2)
	gam = width(x,mrho,gammarho)!gammarho*(ps/pm)**3 
        print*,mrho,gammarho
      endif
	return
	end


c#############################################

      DOUBLE PRECISION FUNCTION THETA(t)
c  Heviside Theta - function
      real*8 t
      IF(t.GE.0.d0)THEN
      Theta = 1.d0
      ELSE
      Theta = 0.d0
      ENDIF
      return
      end
c #############################################

*-----------------------------------------------------------       
      double precision function width(a,mass,widthfix)

      implicit none
      double precision ppi,mass,widthfix,a

      width=widthfix*(ppi(a)/ppi(mass*mass))**3
     &    *sqrt(mass**2/a)
!      print*,widthfix,ppi(a)/ppi(mass*mass),sqrt(mass**2/a)
      return
      end
*-----------------------------------------------------------       
      double precision function ppi(a)
      
      implicit none
	DOUBLE PRECISION s,Pi,mrho,mpi,grhoee,gammarho,alpha
	COMMON/const/s,Pi,mrho,mpi,grhoee,gammarho,alpha
      double precision a


      ppi=0.5d0*dsqrt(a-4.d0*mpi*mpi)

      print*,mpi
      return
      end
c*************************************************************************

	subroutine F_twores(axx,ff_twores)
	implicit none
	real*8 s,Pi,mrho,mpi,grhoee,gammarho,alpha
        COMMON/const/s,Pi,mrho,mpi,grhoee,gammarho,alpha
	real*8 ell,g_phi,fphi,axx,gam1
        real*8 cmod_rpi,cpha_rpi
	COMMON/const3/MPHI,G,GPHIG
	real*8 MPHI,G,GPHIG
	COMPLEX*16 i,g1,ff_twores,c_rpi,c_rho
        real*8 g_rhopig,g_phirpi,Prhores,beta_bro,beta_wphi,
     $            c_modrpi,c_pharpi,c_modwpi,c_phawpi
	common/vmdcom/g_rhopig,g_phirpi,Prhores,beta_bro,beta_wphi,
     $            c_modrpi,c_pharpi,c_modwpi,c_phawpi

        mphi =1.019455d0
        c_modrpi = 0.062d0
        c_pharpi = 3.140d0
        g = 4.26d-3
	i=(0.D0,1.D0)
        gphig = dsqrt(3*mphi**3*1.258404d-6/(4*alpha**2*pi))
	fphi=GPHIG/mphi**2
        
c	fphi=0.d0
        Cmod_rpi   = c_modrpi!.26040d0
        Cpha_rpi   = c_pharpi!3.1112d0
        print*,gphig
	print*,mphi
        print*,Cmod_rpi
	print*,Cpha_rpi
	c_rpi=cmod_rpi*exp(i*cpha_rpi)
	c_rho=c_rpi/g_rhopig/dsqrt(4.*pi*alpha)
c     c_rho= 0.d0
        print*,c_rpi
	print*,c_rho
        print*,mphi
	print*,g
	print*,Prhores
	print*,fphi
c	print*, fphi,g_phirpi,Prhores!,mphi,g,axx

c	                          in this code we chose 
c                              GPHIG=dsqrt(3*mphi**3*GAMMA/(4*alpha**2*pi)),
c                              i.e. ffphi=GPHIG/mphi**2
        prhores = 0.645d0
!        print*,beta_wphi
        ff_twores=fphi*axx*g_phirpi/(-axx+mphi**2-i*mphi*g)*Prhores
c     *           +c_rho
        beta_wphi = 163.d0

	ff_twores=ff_twores*zexp(i*(beta_wphi/180.d0*pi))
        ff_twores=ff_twores + c_rho
        beta_bro = 31.28d0
        print*,c_rho
	print*,beta_wphi
	print*,beta_bro
        ff_twores=ff_twores*zexp(i*(beta_bro/180.d0*pi))
	

	return
	end
*---------------------------------------------
      
*--------------------------------*
*     END FASTERD CORRECTIONS    *
*--------------------------------*


      subroutine ee3gexact(p1,p2,k1,k2,k3,ampl)
c      implicit none
      implicit double precision (Z)
      double precision A11,A12,A13,A14,A15,A16
      double precision     A22,A23,A24,A25,A26
      double precision         A33,A34,A35,A36
      double precision             A44,A45,A46
      double precision                 A55,A56
      double precision                     A66

      double precision p1(0:3),p2(0:3)
      double precision k1(0:3),k2(0:3),k3(0:3)
      double precision k1o(0:3),k2o(0:3),k3o(0:3)
      double precision ampl,ampl2
      double precision p1p2,p1k1,p1k2,p1k3
      double precision      p2k1,p2k2,p2k3
      double precision           k1k2,k1k3
      double precision                k2k3
      double precision e2,e4,e6,me2,me4,me6,dot
      external dot
      integer ifirst
      common/firstee3gexact/e2,e4,e6,me2,me4,me6,ifirst
      data ifirst /0/
      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi

      if (ifirst.eq.0) then
         e2 = 4.d0*pi*alpha
         e4 = e2*e2
         e6 = e4*e2
         me2 = ame*ame
         me4 = me2*me2
         me6 = me4*me2
c         me4 = 0.d0
c         me6 = 0.d0
         ifirst = 1
      endif

      p1p2 = dot(p1,p2)
      p1k1 = dot(p1,k1)
      p1k2 = dot(p1,k2)
      p1k3 = dot(p1,k3)

      p2k1 = dot(p2,k1)
      p2k2 = dot(p2,k2)
      p2k3 = dot(p2,k3)

      k1k2 = dot(k1,k2)
      k1k3 = dot(k1,k3)

      k2k3 = dot(k2,k3)

      include 'form/ee3gexact.F'

      ampl = -ampl

      return
      end
*****************
      SUBROUTINE EE3G(PPDK1,PMDK1,PPDK2,PMDK2,PPDK3,PMDK3,PPDPM,
     .             FUN)
      IMPLICIT REAL*8 (A-H,O-Z)
      common/funvalues/funa,funb,func,funm
      common/parameters/ame,ammu,convfac,alpha,pi
*
      AME2 = AME**2
      AME4 = AME2*ame2
      E2=ALPHA*4.D0*PI
      FUN0=0.D0
      FUNM=0.D0
*
*---- WHITHOUT MASS TERMS
c      FUN0=FUN0+(PPDK1**2+PMDK1**2)/PPDK2/PMDK2/PPDK3/PMDK3
c      FUN0=FUN0+(PPDK2**2+PMDK2**2)/PPDK1/PMDK1/PPDK3/PMDK3
c      FUN0=FUN0+(PPDK3**2+PMDK3**2)/PPDK1/PMDK1/PPDK2/PMDK2

      funa = (PPDK1**2+PMDK1**2)/PPDK2/PMDK2/PPDK3/PMDK3
      funb = (PPDK2**2+PMDK2**2)/PPDK1/PMDK1/PPDK3/PMDK3
      func = (PPDK3**2+PMDK3**2)/PPDK1/PMDK1/PPDK2/PMDK2
      fun0 = funa + funb + func

***************!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      fun0=(PPDK3**2+PMDK3**2)/PPDK1/PMDK1/PPDK2/PMDK2
***************!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUN0=2.D0*FUN0*PPDPM*(E2**3)
*----- MASS TERMS
      FUNM=FUNM+(PMDK2/PMDK3+PMDK3/PMDK2)/PPDK1/PPDK1
      FUNM=FUNM+(PMDK1/PMDK3+PMDK3/PMDK1)/PPDK2/PPDK2
      FUNM=FUNM+(PMDK1/PMDK2+PMDK2/PMDK1)/PPDK3/PPDK3
      FUNM=FUNM+(PPDK2/PPDK3+PPDK3/PPDK2)/PMDK1/PMDK1
      FUNM=FUNM+(PPDK1/PPDK3+PPDK3/PPDK1)/PMDK2/PMDK2
      FUNM=FUNM+(PPDK1/PPDK2+PPDK2/PPDK1)/PMDK3/PMDK3
      FUNM=-2.D0*AME2*(E2)**3*FUNM

*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!**!
c      FUNM=(PPDK1/PPDK3+PPDK3/PPDK1)/PMDK2/PMDK2
c      FUNM=FUNM+(PMDK1/PMDK3+PMDK3/PMDK1)/PPDK2/PPDK2
c      FUNM=-2.D0*AME2*(E2)**3*FUNM
*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!**!
*----- FULL MATRIX ELEMENT
      FUN=FUN0 + FUNM
      fun = fun * 4.d0
      RETURN
      END

*****************
      SUBROUTINE EEG0(PPDQM,PMDQM,QMDK,PPDPM,PPDK,
     #                PMDK,QPDK,S,T,U,SP,TP,UP,FLUX,SQLUMI,
     #                elm)
! essentialy from LABSPV
      implicit double precision (A-H,O-Z)
      common/parameters/ame,ammu,convfac,alpha,pi

      AME2 = AME**2
!      AME4 = AME2**2
*-----FOR ELECTRON RADIATION
      PPDK2  = PPDK**2
      PMDKM1 = 1.D0/PMDK
      PMDKM2 = PMDKM1/PMDK
      QMDKM1 = 1.D0/QMDK
      QMDKM2 = QMDKM1/QMDK
*-----FOR POSITRON RADIATION
      PMDK2  = PMDK**2
      QPDKM1 = 1.D0/QPDK
      QPDKM2 = QPDKM1/QPDK
      PPDKM1 = 1.D0/PPDK
      PPDKM2 = PPDKM1/PPDK
*
      S2  = S**2
      SP2 = SP**2
      T2  = T**2
      TP2 = TP**2
      U2  = U**2
      UP2 = UP**2
*
      W1 = S/PPDK/PMDK
      W2 = SP/QPDK/QMDK
      W3 = -T/PPDK/QPDK-TP/PMDK/QMDK+U/PPDK/QMDK+UP/PMDK/QPDK
*
      SQMBER0 = (W1+W2+W3)*(S*SP*(S*S+SP*SP)
     #                     +T*TP*(T*T+TP*TP)
     #                     +U*UP*(U*U+UP*UP))/S/SP/T/TP
      SQMBERM =
     #          PMDKM2*(SP/T  + T/SP  + 1.D0)**2
     #        + PPDKM2*(SP/TP + TP/SP + 1.D0)**2
     #        + QMDKM2*(S/T   + T/S   + 1.D0)**2
     #        + QPDKM2*(S/TP  + TP/S  + 1.D0)**2
      SQMBERM = -4.D0*AME2*SQMBERM
*
      SQMBER = (SQMBER0 + SQMBERM)*(4.D0*PI*ALPHA)**3
      sqmber = 4.d0 * sqmber ! it has to be summed over all spins! (also init.)
      elm    = sqmber
      SQLUMI = SQMBER*FLUX
      RETURN
      END
*******
      subroutine approxNg(ng,p1,p2,p3,p4,qph,me)
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),qph(40,0:3),qt(40,0:3)
      dimension tmp(0:3)
      common/parameters/ame,ammu,convfac,alpha,pi
      e2 = alpha*4.d0*pi
      nph = ng + 2
      me = 0.d0
      do i = 0,3
         tmp(i)  = p1(i) + p2(i)
         qt(1,i) = p3(i)
         qt(2,i) = p4(i)
         if (ng.gt.0) then
            do k = 1,ng
               qt(k+2,i) = qph(k,i)
            enddo
         endif
      enddo
      s = dot(tmp,tmp)

      sum  = 0.d0
      prod = 1.d0
      do k = 1,nph
         do i = 0,3
            tmp(i) = qt(k,i)
         enddo
         xi = dot(tmp,p1)
         yi = dot(tmp,p2)
         sum = sum + xi*yi * (xi**2 + yi**2)
         prod = prod * xi*yi
      enddo
      me = sum/prod * s**(nph-2) * e2**nph
      return
      end

      subroutine approxNg_SINGLE(ng,p1,p2,p3,p4,qph,me)
! only one term
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),qph(40,0:3),qt(40,0:3)
      dimension tmp(0:3)
      common/parameters/ame,ammu,convfac,alpha,pi
      e2 = alpha*4.d0*pi
      nph = ng + 2
      me = 0.d0
      do i = 0,3
         tmp(i)  = p1(i) + p2(i)
         qt(1,i) = p3(i)
         qt(2,i) = p4(i)
         if (ng.gt.0) then
            do k = 1,ng
               qt(k+2,i) = qph(k,i)
            enddo
         endif
      enddo
      s = dot(tmp,tmp)

      sum  = 0.d0
      prod = 1.d0
      do k = 2,2
         do i = 0,3
            tmp(i) = qt(k,i)
         enddo
         xi = dot(tmp,p1)
         yi = dot(tmp,p2)
         sum = sum + xi*yi * (xi**2 + yi**2)
      enddo
      do k = 1,nph
         do i = 0,3
            tmp(i) = qt(k,i)
         enddo
         xi = dot(tmp,p1)
         yi = dot(tmp,p2)
         prod = prod * xi*yi
      enddo
      me = sum/prod * s**(nph-2) * e2**nph  * 8.d0
      return
      end
******************************
