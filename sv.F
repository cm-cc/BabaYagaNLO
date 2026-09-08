!      COMPLEX VACUUMPOLARIZATION: ALPHA(0)->ALPHA(0)*VPOLC(Q2)
      FUNCTION VPOLC(Q2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 QIN4,ST24,DER4,ERRDER4,DEG4,ERRDEG4
      DIMENSION AMASSES(9)
      common/parameters/ame,ammu,convfac,alpha,pi
      common/ialpharunning/iarun
      common/teubner/iteubn
      integer inskvp
      common/nskvp/inskvp
      common/resonances/ires,ionDHAD
      common/vpolcomc/ifirst
      common/jpsiparameters/amjpsi,gtjpsi,gejpsi,effcjpsi
      common/jpsiparameters2/amjpsi2s,gtjpsi2s,gejpsi2s,effcjpsi2s
      common/jpsiparameters3/amjpsi3s,gtjpsi3s,gejpsi3s,effcjpsi3s
      common/jpsiparameters4/amjpsi4s,gtjpsi4s,gejpsi4s,effcjpsi4s
      common/jpsiparameters5/amjpsi5s,gtjpsi5s,gejpsi5s,effcjpsi5s

      data ifirst /0/

      dimension iswvperr(3),iswvperrt(6)
      common/iswitchvperror/iswvperr,iswvperrt

      common/rescaledhaderr/scal_err_dhad

      double complex vpolcnsk,vpolchlmnt,VPOLhadrnskstr2020
      external vpolcnsk,vpolchlmnt,VPOLhadrnskstr2020

      double complex vpolc,deltahad,im,cder

      common/dalphacmn/dalphaQED,der,errder,toppart
      integer ioffhvp,inlovp
      common/offhadrVP_nlovp/ioffhvp,inlovp

      parameter (nres=15)
      double precision AM(nres),GT(nres),GE(nres)
      common/resonancesparams/am,gt,ge
*     AM=HADRON MASSES, GT=TOTAL WIDTH, GE=PARTIAL E+E- OR Mu+Mu- WIDTH.
* 1=RHO, 2=OMEGA, 3=PHI, 4=J/PSI, 5=PSI(2S), 6=PSI(3770), 7=PSI(4040),
* 8=PSI(4160), 9=PSI(4415), 10=Y, 11=Y(2S), 12=Y(3S), 13=Y(4S), 14=Y(10860),
* 15=Y(11020).
      double precision EAM(nres),EGT(nres),EGE(nres)
      common/resonancesparamserrors/eam,egt,ege
** errors on the previous parameters, from PDG

****  TO BE UPDATED (e.g. from PDG)!!
      DATA AM/ ! GeV
     >     0.769D0,    ! 1.  rho        !*  0.76419892698171199
     >     0.78265D0,  ! 2.  omega      !*  0.78278988240778891
     >     1.019455D0, ! 3.  phi        !*  1.0198676384707968
     >     3.096916D0, ! 4.  J/Psi      !*  3.0969197600196230
     >     3.686109D0, ! 5.  Psi(2S)    !*  3.6861561280010915
     >     3.77315D0,  ! 6.  Psi(3770)
     >     4.039D0,    ! 7.  Psi(4040)
     >     4.191D0,    ! 8.  Psi(4160)
     >     4.421D0,    ! 9.  Psi(4415)
     >     9.4603D0,   ! 10. Y          !*   9.4603013694067908
     >     10.02326D0, ! 11. Y(2S)      !*  10.023263041544904
     >     10.35520,   ! 12. Y(3S)      !*  10.355200529202707
     >     10.5794D0,  ! 13. Y(4S)
     >     10.876D0,   ! 14. Y(10860)
     >     11.019D0/   ! 15. Y(11020)
      DATA GT/ ! GeV
     >     0.149D0,   ! 1.  rho
     >     8.49D-3,   ! 2.  omega
     >     4.266D-3,  ! 3.  phi
     >     92.9D-6,   ! 4.  J/Psi
     >     299.D-6,   ! 5.  Psi(2S)
     >     27.2D-3,   ! 6.  Psi(3770)
     >     80.D-3,    ! 7.  Psi(4040)
     >     70.D-3,    ! 8.  Psi(4160)
     >     62.D-3,    ! 9.  Psi(4415)
     >     54.02D-6,  ! 10. Y
     >     31.98D-6,  ! 11. Y(2S)
     >     20.32D-6,  ! 12. Y(3S)
     >     20.5D-3,   ! 13. Y(4S)
     >     55.D-3,    ! 14. Y(10860)
     >     79.D-3/    ! 15. Y(11020)
      DATA GE/ ! GeV
     >     7.04D-6,   ! 1.  rho
     >     0.60D-6,   ! 2.  omega
     >     1.27D-6,   ! 3.  phi
     >     5.55D-6,   ! 4.  J/Psi
     >     2.36D-6,   ! 5.  Psi(2S)
     >     0.262D-6,  ! 6.  Psi(3770)
     >     0.86D-6,   ! 7.  Psi(4040)
     >     0.48D-6,   ! 8.  Psi(4160)
     >     0.58D-6,   ! 9.  Psi(4415)
     >     1.34D-6,   ! 10. Y
     >     0.612D-6,  ! 11. Y(2S)
     >     0.443D-6,  ! 12. Y(3S)
     >     0.272D-6,  ! 13. Y(4S)
     >     0.31D-6,   ! 14. Y(10860)
     >     0.130D-6/  ! 15. Y(11020)
*** ERRORS ON THE PREVIOUS PARAMETERS, from PDG
      DATA EAM/                  ! GeV
     >     0.0009D0,   ! 1.  rho
     >     0.00012D0,  ! 2.  omega
     >     0.000019D0, ! 3.  phi
     >     0.000011D0, ! 4.  J/Psi
     >     0.000013D0, ! 5.  Psi(2S)
     >     0.00033D0,  ! 6.  Psi(3770)
     >     0.001D0,    ! 7.  Psi(4040)
     >     0.005D0,    ! 8.  Psi(4160)
     >     0.004D0,    ! 9.  Psi(4415)
     >     0.00026D0,  ! 10. Y
     >     0.00031D0,  ! 11. Y(2S)
     >     0.0005D0,   ! 12. Y(3S)
     >     0.0012D0,   ! 13. Y(4S)
     >     0.011D0,    ! 14. Y(10860)
     >     0.008D0/    ! 15. Y(11020)
      DATA EGT/ ! GeV
     >     0.001D0,  ! 1.  rho
     >     0.08D-3,  ! 2.  omega
     >     0.031D-3, ! 3.  phi
     >     2.8D-6,   ! 4.  J/Psi
     >     8.0D-6,   ! 5.  Psi(2S)
     >     1.0D-3,   ! 6.  Psi(3770)
     >    10.D-3,    ! 7.  Psi(4040)
     >    10.D-3,    ! 8.  Psi(4160)
     >    20.D-3,    ! 9.  Psi(4415)
     >     1.25D-6,  ! 10. Y
     >     2.63D-6,  ! 11. Y(2S)
     >     1.85D-6,  ! 12. Y(3S)
     >     2.5D-3,   ! 13. Y(4S)
     >    28.D-3,    ! 14. Y(10860)
     >    16.D-3/    ! 15. Y(11020)
      DATA EGE/ ! GeV
     >     0.06D-6,  ! 1.  rho
     >     0.02D-6,  ! 2.  omega
     >     0.04D-6,  ! 3.  phi
     >     0.16D-6,  ! 4.  J/Psi
     >     0.04D-6,  ! 5.  Psi(2S)
     >     0.018D-6, ! 6.  Psi(3770)
     >     0.07D-6,  ! 7.  Psi(4040)
     >     0.22D-6,  ! 8.  Psi(4160)
     >     0.07D-6,  ! 9.  Psi(4415)
     >     0.018D-6, ! 10. Y
     >     0.011D-6, ! 11. Y(2S)
     >     0.008D-6, ! 12. Y(3S)
     >     0.029D-6, ! 13. Y(4S)
     >     0.07D-6,  ! 14. Y(10860)
     >     0.03D-6/  ! 15. Y(11020)

*
* E,MU,TAU,TOP MASSES:
*
      DATA AMASSES /0.510999D-3,.1056584D0,1.77686d0,172.5D0,
     . 0.066d0,! u
     . 0.066d0,! d
     . 1.2d0,  ! c
     . 0.150d0,! s
     . 4.6d0  /! b
*
      double complex dalep,dahadr
      common/forfixedorder/dalep,dahadr

      if (ifirst.eq.0) then

c         print*,'resetting e mu masses in loops as the ones in input!'
         amasses(1) = ame
         amasses(2) = ammu

         dalep = (0.d0,0.d0)
         dhadr = (0.d0,0.d0)

! filling the common for importance sampling...
         ifirst = 1

         amjpsi   = am(4)
         gtjpsi   = gt(4)
         gejpsi   = ge(4)
         effcjpsi = 3.d0*gejpsi/amjpsi/alpha

         amjpsi2s   = am(5)
         gtjpsi2s   = gt(5)
         gejpsi2s   = ge(5)
         effcjpsi2s = 3.d0*gejpsi2s/amjpsi2s/alpha

         amjpsi3s   = am(10)
         gtjpsi3s   = gt(10)
         gejpsi3s   = ge(10)
         effcjpsi3s = 3.d0*gejpsi3s/amjpsi3s/alpha
         amjpsi4s   = am(11)
         gtjpsi4s   = gt(11)
         gejpsi4s   = ge(11)
         effcjpsi4s = 3.d0*gejpsi4s/amjpsi4s/alpha
         amjpsi5s   = am(12)
         gtjpsi5s   = gt(12)
         gejpsi5s   = ge(12)
         effcjpsi5s = 3.d0*gejpsi5s/amjpsi5s/alpha

         ires    = 0
         ionDHAD = 1
         if (ires.eq.1) ionDHAD = 0

      endif

      if (iarun.eq.0) then
         vpolc = 1.d0
         return
      endif

c$$$      if (q2.gt.15) then
c$$$cc         vpolc = 1.d0 + 1d0*(0.d0,1.d0)
c$$$         vpolc = 1.d0
c$$$c         continue
c$$$         return
c$$$      else
c$$$         vpolc = 0.d0
c$$$c         continue
c$$$         return
c$$$      endif
c      vpolc = 1.d0
c      return

**** strong2020 implementation by Fedor
      if (inskvp.eq.1) then
         SOMMA = 0.D0
         DO I=1,4
            xxx   = SUMMA(AMASSES(I),Q2,I)
            SOMMA = SOMMA+xxx
         ENDDO
         dalphaQED = ALPHA/PI*SOMMA
         toppart   = alpha/pi * xxx

         dalep = dalphaQED
         vpolc = (0.d0,0.d0)
         if (ioffhvp.eq.0) vpolc = VPOLhadrnskstr2020(q2) ! this is the hadronic part
         dahadr = vpolc

         vpolc = 1.d0/(1.d0 - vpolc - dalphaQED)
         return
      endif
**** my implementation
c      if (inskvp.eq.1) then
c         vpolc = VPOLCnsk(Q2)
c         dalep = 1.d0-1.d0/vpolc
c         return
c      endif
******************************************
      if (iteubn.eq.1) then
         vpolc = VPOLChlmnt(Q2)
         dalep = 1.d0-1.d0/vpolc
         return
      endif

      SOMMA=0.D0
      DO I=1,4
         xxx  = SUMMA(AMASSES(I),Q2,I)
         SOMMA=SOMMA+xxx
      ENDDO

      dalphaQED = ALPHA/PI*SOMMA
      toppart = alpha/pi * xxx

      AMODQ2=ABS(Q2)
      ST2 = 0.2322d0
      qin = q2/abs(q2) * sqrt(abs(q2))

      der = 0.d0
      deltahad = 0.d0
ccc      goto 333

      if (ionDHAD.eq.1.or.qin.lt.0.d0) then
         if (iteubn.eq.0) then
c            CALL DHADR5n12(QIN,ST2,DER,ERRDER,DEG,ERRDEG)
c            der12 = der
c            errder12 = errder

c            call dhadr5n17(qin,st2,der,errdersta,errdersys,deg,
c     .           errdegsta,errdegsys) ! it's the 17 version
c            errder = sqrt(errdersta**2+errdersys**2)
c            err17 = errder
c            der17 = der
c            CALL DHADR5n16(QIN,ST2,DER,ERRDER,DEG,ERRDEG) ! it's the 16 version


            der = 0.d0
            errder = 0.d0
            errderst = 0.d0
            deg = 0.d0
            errdeg = 0.d0
            errdegst = 0.d0
            if (ioffhvp.eq.0) then
               CALL DHADR5x(QIN,ST2,DER,ERRDER,errderst,DEG,ERRDEG,
     .              errdegst)   ! it's the 23 version
            endif
c            errsum = errdersta+errdersys
c            print*,' '
c            print*,qin,der17/der,errsum/errder,err17/errder
c            print*,qin,der12/der,errder12/errder
c            print*,'     ',der,der12,der17

            if (qin.lt.0.d0) then
               der = der  + iswvperr(3) * scal_err_dhad * errder
c               der = 0.d0
            else
               der = der  + iswvperr(1) * scal_err_dhad * errder
            endif
         elseif (iteubn.eq.1) then
            das = 0.d0
            ddas = 0.d0
            dat = 0.d0
            ddat = 0.d0
            if (ioffhvp.eq.0) then
               call d_alpha(dabs(qin),das,ddas,dat,ddat)
            endif
            if (qin.lt.0.d0) then
               der = alpha*(das + iswvperr(3) * scal_err_dhad * ddas)
            else
               der = alpha*(dat + iswvperr(1) * scal_err_dhad * ddat)
            endif
         endif
      endif

      if (ires.eq.1.and.qin.gt.0.d0) then ! ires = 0 by default
         im = (0.d0,1.d0)
         do j=1,nres
            AMC = AM(J)
            GEC = GE(J)
            GTC = GT(J)
            V   = 3.D0*GEC/AMC / ALPHA
            deltahad = deltahad +
     .          V * Q2 / (Q2-AMC*AMC+im*GTC*AMC)
         enddo
      endif

c      if (qin.lt.0.d0) der = 0.d0

 333  continue

      cder = deltahad + der

c     vpolc = 1.d0/(1.d0-dalphaQED-cder)
c     vpolc = 1.d0 + cder
c     vpolc = 1.d0/(1.d0-dalphaQED) + cder ! this has the same peaks as NSK routine if using the resonances!

      vpolc = 1.d0/(1.d0 - dalphaQED - cder)

      dalep  = dalphaQED
      dahadr = cder

      RETURN
      END
******
******
* COMPLEX VACUUM POLARIZATION, T. Teubner & al.
      FUNCTION VPOLChlmnt(Q2)
      implicit double precision (A-H,L,O-Z)
      parameter (np  = 7186)
c      parameter (np=6883)
      parameter (nw  = 1582)

c      parameter (npt = np) ! or
c      parameter (npt = nw) ! or
      parameter (npt = nw+np)

      dimension s(0:np-1)
      dimension st(0:npt-1)
      dimension pres(0:np-1),pre(0:np-1),pim(0:np-1)
      dimension epres(0:np-1),epre(0:np-1),epim(0:np-1),cv(0:np-1)
      double complex vpolcnsk,vpolchlmnt,deltahad,im,res
      double complex cP,ressp,restm,ctmp
      external vpolcnsk
      common/vpcnsk/s,pres,pre,pim,epres,epre,epim,cv,im,ifirstnsk

      parameter (nres=15)
      double precision AM(nres),GT(nres),GE(nres)
      common/resonancesparams/am,gt,ge
      double precision EAM(nres),EGT(nres),EGE(nres)
      common/resonancesparamserrors/eam,egt,ege

      dimension iresbh(8) ! the indexes of the resonances to be added by hand (Thomas said to include **only** 1-6)
      dimension prest(0:npt-1),pret(0:npt-1),pimt(0:npt-1)
      dimension eprest(0:npt-1),epret(0:npt-1),epimt(0:npt-1)
      common/vpcteub/st,prest,pret,pimt,eprest,epret,epimt,iresbh,
     .     ifirst,imidguess
      data ifirst/0/

      double precision wt(nw)
      common/wtabulation/wt

      double precision wtot(npt)
      integer indexes(npt)
      common/wtotcommon/wtot

      common/tabulationlimits/wmin,wmax,ic

      common/parameters/ame,ammu,convfac,alpha,pi
      common/ecms/ecms,ecmsnom

      integer iloc(1)

      dimension iswvperr(3),iswvperrteub(6)
      common/iswitchvperror/iswvperr,iswvperrteub

      common/rescaledhaderr/scal_err_dhad

*     using vp_hlmnt_v2_1.f, sent by T. Teubner, private communication

      if (ifirst.eq.0) then
***   filling the common by calling vpolcnsk once (I need the array s --> st[eubner])
         vpolchlmnt = vpolcnsk(q)
****  filling the common wtabulation
         wtot = 0.d0
         call fillwtabcommon
         wmin = wt(1)
         wmax = wt(nw)
         wtot(1:nw)     = wt**2
         wtot(nw+1:npt) = s
         call sort_vector_inverse(wtot,wtot,indexes,npt)
c         do k = 1,npt
c            print*,k,sqrt(wtot(k)),indexes(k)
c         enddo
c         stop
****
c
c I solve for P, so that is the same quantity returned by NSK
c     I also use the same evaluation points as NSK, and then linear
c     interpolation between them
c     1/(1-P) = 1/(1-danr) + res ->
c     1-P = (1-danr)/(1+(1-danr)*res)
c       P = 1 - (1-danr)/(1+(1-danr)*res)
c         = (1+(1-danr)*res - (1-danr))/(1+(1-danr)*res)
c         = ( danr + res - danr*res) / (1 + res- danr*res)
******************
* filling st, from the nsk data or Teubner tabulation points
c         st(0:np-1) = s   ! or
c         do k = 1,nw
c            st(k-1) = wt(k)*wt(k)
c         enddo ! or
         st = wtot
******************
         prest  = 0.d0
         pret   = 0.d0
         pimt   = 0.d0
         eprest = 0.d0
         epret  = 0.d0
         epimt  = 0.d0
         iresbh(1) = 4
         iresbh(2) = 5
         iresbh(3) = 10
         iresbh(4) = 11
         iresbh(5) = 12
         iresbh(6) = 13
         iresbh(7) = 14
         iresbh(8) = 15
c         print*,'Evaluating Teubner et al. VP array '
c         open(10,file='cfrvp-0-v2-3',status='unknown')
         do k = 0,npt-1
            energy = sqrt(st(k))

c            hlmnt is tabulated up to wmax and from wmin
            if (energy.le.wmax.and.energy.ge.wmin) then
               ressp = 0.d0
               restm = 0.d0
! I do it below, every time the function is called, better for error handling...
c$$$               do i=1,6
c$$$                  j = iresbh(i)
c$$$                  AMC = AM(J) !- EAM(J)
c$$$                  GEC = GE(J) !- EGE(J)
c$$$                  GTC = GT(J) !- EGT(J)
c$$$                  V   = 3.D0*GEC/AMC / ALPHA
c$$$                  restm = restm +
c$$$     .                 V * st(k) / ( st(k)-AMC*AMC+im*GTC*AMC )
c$$$                  ressp = ressp -
c$$$     .                 V * st(k) / (-st(k)-AMC*AMC+im*GTC*AMC )
c$$$               enddo
               nrflag = 1
               call vphlmntv2(energy,
     &              vprehadsp,vprehadtm,vpimhad,vprelepsp,
     &              vpreleptm,vpimlep,
     &              vpretopsp,vpretoptm,nrflag,evprehadsp,
     &              evprehadtm,evpimhad)
c               vprehadsp = vprehadsp - evprehadsp
c               vprehadtm = vprehadtm - evprehadtm
c               vpimhad   = vpimhad -  evpimhad
!     ! space like
               prest(k) = vprehadsp + vprelepsp + vpretopsp
               eprest(k) = evprehadsp
!     !     time like
               nrflag = 0
               call vphlmntv2(energy,
     &        vprehadsp,vprehadtm,vpimhad,vprelepsp,vpreleptm,vpimlep,
     &        vpretopsp,vpretoptm,nrflag,evprehadsp,evprehadtm,evpimhad)

               pret(k)  = vprehadtm + vpreleptm + vpretoptm
               pimt(k)  = vpimhad   + vpimlep

               epret(k) = evprehadtm
               epimt(k) = evpimhad

c     [[ useless, but harmless
c   this is 1/(1-P)
               ctmp = 1.d0/(1.d0 - pret(k)-im*pimt(k)) + restm ! in the time-like not resummed!
c  ctmp = 1/(1-P) -> 1-P = 1/ctmp  -> P = 1-1/ctmp
               pret(k) = DReal(1.d0-1.d0/ctmp)
               pimt(k) = DImag(1.d0-1.d0/ctmp)
c   ]]
c               write(10,*)energy,pret(k),pimt(k),pre(k),pim(k),
c     .              prest(k),pres(k)
c               write(10,*)energy,pret(k),pimt(k),prest(k)
            endif
         enddo
c         close(10)
c         print*,'done!'
c         print*,'stopping line ~420 sv.f'
c         stop
*** guessing where to start the binary search....
         k = 0
         aq2 = ecmsnom**2
         do while(aq2.gt.st(k))
            k = k + 1
         enddo
         k = k - 1
         imidguess = k ! this is correct for born s-channel and no energy spread
***************
         ic = 0
         ifirst=1
      endif
c      k   = 0
c      sk  = 0.d0
c      aq2 = abs(q2)
c      do while(aq2.gt.st(k))
c         k = k + 1
c      enddo
c      k = k - 1
******** binary search ****
      aq2     = abs(q2)
      istart  = 0
      ifinish = npt-1
c      imid = (istart + ifinish)/2
      imid = imidguess

c     do while(ifinish.ge.istart)
      do while(.not.(st(imid+1).gt.aq2.and.st(imid).le.aq2))
         if (aq2.gt.st(imid)) then
            istart  = imid + 1
         else
            ifinish = imid - 1
         endif
         imid = (istart + ifinish)/2
      enddo
      k = imid
******************
      if (q2.gt.0.d0) then
         vpre = pret(k)+(pret(k+1)-pret(k))/(st(k+1)-st(k)) * (q2-st(k))
         vpim = pimt(k)+(pimt(k+1)-pimt(k))/(st(k+1)-st(k)) * (q2-st(k))

         evpre=epret(k)+(epret(k+1)-epret(k))/(st(k+1)-st(k))*(q2-st(k))
         evpim=epimt(k)+(epimt(k+1)-epimt(k))/(st(k+1)-st(k))*(q2-st(k))

         ctmp = 0.d0
         do i=1,6
            j    = iresbh(i)
            AMC  = AM(J) + EAM(J)*iswvperrteub(4)*scal_err_dhad
            GEC  = GE(J) + EGE(J)*iswvperrteub(5)*scal_err_dhad
            GTC  = GT(J) + EGT(J)*iswvperrteub(6)*scal_err_dhad
            V    = 3.D0*GEC/AMC / ALPHA
            ctmp = ctmp +
     .           V * q2 / ( q2 - AMC*AMC+im*GTC*AMC )
         enddo

         vpre = vpre + iswvperrteub(1)*scal_err_dhad*evpre
         vpim = vpre + iswvperrteub(2)*scal_err_dhad*evpim

      else
         vpre = prest(k) +
     .        (prest(k+1)-prest(k))/(st(k+1)-st(k))*(aq2-st(k))
         evpre = eprest(k) +
     .        (eprest(k+1)-eprest(k))/(st(k+1)-st(k))*(aq2-st(k))

         vpim = 0.d0
         evpim = 0.d0

         ctmp = 0.d0

         vpre = vpre + iswvperrteub(3)*scal_err_dhad*evpre
      endif
*
! ctmp is 0 in the spacelike, not 0 in the timelike. It has not to be resummed
      VPOLChlmnt = 1.d0/(1.d0  - vpre - im*vpim) + ctmp
      return
      end
*
*     COMPLEX VACUUM POLARIZATION from NSBSK DATA: ALPHA(0)->ALPHA(0)*VPOL(Q2)
      FUNCTION VPOLCnsk(Q2)
      implicit double precision (A-H,L,O-Z)
      parameter (np=7186)
c      parameter (np=6883)
      dimension s(0:np-1)
      dimension pres(0:np-1),pre(0:np-1),pim(0:np-1)
      dimension epres(0:np-1),epre(0:np-1),epim(0:np-1),cv(0:np-1)
      double complex vpolcnsk,deltahad,im
      common/vpcnsk/s,pres,pre,pim,epres,epre,epim,cv,im,ifirst
      data ifirst/0/

      dimension iswvperr(3),iswvperrt(6)
      common/iswitchvperror/iswvperr,iswvperrt

      common/rescaledhaderr/scal_err_dhad
* from    http://cmd.inp.nsk.su/~ignatov/vpl/

      if (ifirst.eq.0) then
         im = (0.d0,1.d0)
c         open(45,file='vpol_novosibirsk_NEW.dat',status='old')
         open(45,file='vpol_novosibirsk.dat',status='old')
         do k = 0,np-1
            read(45,*)s(k),pres(k),pre(k),pim(k),epres(k),epre(k),
     .           epim(k),cv(k)
         enddo
         close(45)
         ifirst=1
      endif

      k   = 0
      sk  = 0.d0
      aq2 = abs(q2)
      do while(aq2.gt.s(k))
         k = k + 1
      enddo
      k = k - 1

      if (q2.gt.0.d0) then
         vpre = pre(k) + (pre(k+1)-pre(k))/(s(k+1)-s(k)) * (q2-s(k))
         vpim = pim(k) + (pim(k+1)-pim(k))/(s(k+1)-s(k)) * (q2-s(k))

         evpre = epre(k) + (epre(k+1)-epre(k))/(s(k+1)-s(k)) * (q2-s(k))
         evpim = epim(k) + (epim(k+1)-epim(k))/(s(k+1)-s(k)) * (q2-s(k))

         rho = cv(k) + (cv(k+1)-cv(k))/(s(k+1)-s(k)) * (q2-s(k))
         cov = rho * evpre*evpim
      else
         vpre = pres(k) + (pres(k+1)-pres(k))/(s(k+1)-s(k))*(aq2-s(k))
         vpim = 0.d0

         evpre = epres(k)+(epres(k+1)-epres(k))/(s(k+1)-s(k))*(aq2-s(k))
         evpim = 0.d0

         rho = 0.d0
         cov = 0.d0
      endif

c      print*,iswvperr

      if (abs(scal_err_dhad).gt.1d-4) then
c     vpre = vpre + scal_err_dhad * evpre
c     vpim = vpim + scal_err_dhad * evpim
*
*     USING FORMALISM FROM GLEN COWAN, "STATISTICAL DATA ANALYSIS", p.22+, OXFORD SCIENCE PUBLICATIONS
         if (q2.lt.0.d0) then

            vpre = vpre + scal_err_dhad * evpre * iswvperr(3)

         elseif (iswvperr(1).ne.0.or.iswvperr(2).ne.0) then
            x1   = vpre
            x2   = vpim
            s12  = evpre*evpre
            s22  = evpim*evpim
            rho2 = rho * rho
            sq   = sqrt((s12+s22)**2-4.d0*(1.d0-rho2)*s12*s22)
            lp   = 0.5d0*(s12 + s22 + sq)
            lm   = lp - sq
            th   = 0.5d0 * atan(2.d0*rho*evpre*evpim/(s12-s22))
            cs   = cos(th)
            ss   = sin(th)
            y1   =  cs * x1 + ss * x2
            y2   = -ss * x1 + cs * x2

            dy1 = scal_err_dhad * iswvperr(1) * sqrt(lp)
            dy2 = scal_err_dhad * iswvperr(2) * sqrt(lm)
! changing y1 and y2...
            y1   = y1 + dy1
            y2   = y2 + dy2
! ... going back to x1 and x2
            x1   = cs * y1 - ss * y2
            x2   = ss * y1 + cs * y2

!     updating the variables
            vpre = x1
            vpim = x2
         endif
      endif
*********************************************************************************************
      VPOLCnsk = 1.d0/(1.d0  - vpre - im*vpim)
      return
      end
*     COMPLEX VACUUM POLARIZATION from DATA, Fedor's version: ALPHA(0)->ALPHA(0)/(1-VPOL(Q2))
*     data from  http://cmd.inp.nsk.su/~ignatov/vpl/
*     VP fortran function by BabaYaga@NLO
*****************************************************************************************
      subroutine VPOLnskinitstr2020
      implicit none

      double precision s,pres,pre,pim,epres,epre,epim,cv
      double precision lpres,lpre,lpim
      integer ifirst,np,k
      parameter (np=7360)
      dimension s(0:np-1)
      dimension pres(0:np-1),pre(0:np-1),pim(0:np-1)
      dimension lpres(0:np-1),lpre(0:np-1),lpim(0:np-1)
      dimension epres(0:np-1),epre(0:np-1),epim(0:np-1),cv(0:np-1)
      data ifirst/0/
      common/vpoltblnsk/s,pres,pre,pim,epres,epre,epim,cv,
     &     lpres,lpre,lpim,ifirst
      double precision stmp,tmp


      if (ifirst.eq.0) then
         open(45,file='vpol_all_bare_sum_v2.9.dat',status='old')
         open(46,file='vpol_bare_lept_v2.9.dat',status='old')
         do k = 0,np-1
            read(45,*)s(k),pres(k),pre(k),pim(k),epres(k),epre(k),
     .           epim(k),cv(k)
            read(46,*)stmp,lpres(k),lpre(k),lpim(k),tmp,tmp,tmp,tmp
            if(stmp.ne.s(k)) print *,'lept. and full VP are not matched'
         enddo
         close(45)
         close(46)
         ifirst=1
      endif
      end

      FUNCTION VPOLnskstr2020(Q2)
      implicit none
      double precision s,pres,pre,pim,epres,epre,epim,cv
      double precision lpres,lpre,lpim
      integer ifirst,np
      parameter (np=7360)
      dimension s(0:np-1)
      dimension pres(0:np-1),pre(0:np-1),pim(0:np-1)
      dimension lpres(0:np-1),lpre(0:np-1),lpim(0:np-1)
      dimension epres(0:np-1),epre(0:np-1),epim(0:np-1),cv(0:np-1)
      data ifirst/0/
      common/vpoltblnsk/s,pres,pre,pim,epres,epre,epim,cv,
     &     lpres,lpre,lpim,ifirst

      double precision q2,aq2,vpre,vpim,evpre,evpim
      double complex vpolnskstr2020,im
      integer istart,ifinish,imid,k

************************************************************************
************************************************************************
      if (ifirst.eq.0) then
         call VPOLnskinitstr2020
      endif

** binary search ********************************
      aq2 = abs(q2)
      istart  = 0
      ifinish = np - 1
      imid  = np/2
      do while(.not.(s(imid+1).gt.aq2.and.s(imid).le.aq2))
         if (aq2.gt.s(imid)) then
            istart = imid + 1
         else
            ifinish = imid - 1
         endif
         imid = (istart + ifinish)/2
      enddo
      k = imid
**************************************************
      if (q2.gt.0.d0) then
         vpre = pre(k) + (pre(k+1)-pre(k))/(s(k+1)-s(k)) * (q2-s(k))
         vpim = pim(k) + (pim(k+1)-pim(k))/(s(k+1)-s(k)) * (q2-s(k))

C uncertinty of VP can be also returned
         evpre = epre(k) + (epre(k+1)-epre(k))/(s(k+1)-s(k)) * (q2-s(k))
         evpim = epim(k) + (epim(k+1)-epim(k))/(s(k+1)-s(k)) * (q2-s(k))

c         rho = cv(k) + (cv(k+1)-cv(k))/(s(k+1)-s(k)) * (q2-s(k))
c         cov = rho * evpre*evpim
      else
         vpre = pres(k) + (pres(k+1)-pres(k))/(s(k+1)-s(k))*(aq2-s(k))
         vpim = 0.d0

         evpre = epres(k)+(epres(k+1)-epres(k))/(s(k+1)-s(k))*(aq2-s(k))
         evpim = 0.d0

c         rho = 0.d0
c         cov = 0.d0
      endif
*********************************************************************************************
      im = (0.d0,1.d0)
      VPOLnskstr2020 = vpre + im*vpim
      return
      end

      FUNCTION VPOLleptnskstr2020(Q2)
      implicit none
      double precision s,pres,pre,pim,epres,epre,epim,cv
      double precision lpres,lpre,lpim
      integer ifirst,np
      parameter (np=7360)
      dimension s(0:np-1)
      dimension pres(0:np-1),pre(0:np-1),pim(0:np-1)
      dimension lpres(0:np-1),lpre(0:np-1),lpim(0:np-1)
      dimension epres(0:np-1),epre(0:np-1),epim(0:np-1),cv(0:np-1)
      data ifirst/0/
      common/vpoltblnsk/s,pres,pre,pim,epres,epre,epim,cv,
     &     lpres,lpre,lpim,ifirst

      double precision q2,aq2,vpre,vpim,evpre,evpim
      double complex vpolleptnskstr2020,im
      integer istart,ifinish,imid,k

************************************************************************
************************************************************************
      if (ifirst.eq.0) then
         call VPOLnskinitstr2020
      endif

** binary search ********************************
      aq2 = abs(q2)
      istart  = 0
      ifinish = np - 1
      imid  = np/2
      do while(.not.(s(imid+1).gt.aq2.and.s(imid).le.aq2))
         if (aq2.gt.s(imid)) then
            istart = imid + 1
         else
            ifinish = imid - 1
         endif
         imid = (istart + ifinish)/2
      enddo
      k = imid
**************************************************
      if (q2.gt.0.d0) then
        vpre = lpre(k) + (lpre(k+1)-lpre(k))/(s(k+1)-s(k)) * (q2-s(k))
        vpim = lpim(k) + (lpim(k+1)-lpim(k))/(s(k+1)-s(k)) * (q2-s(k))
      else
        vpre = lpres(k) + (lpres(k+1)-lpres(k))/(s(k+1)-s(k))*(aq2-s(k))
      endif
*********************************************************************************************
      im = (0.d0,1.d0)
      VPOLleptnskstr2020 = vpre + im*vpim
      return
      end

      FUNCTION VPOLhadrnskstr2020(Q2)
      implicit none
      double precision q2
      double complex vpolhadrnskstr2020,vpolleptnskstr2020,
     .     vpolnskstr2020

      vpolhadrnskstr2020 = vpolnskstr2020(q2)-vpolleptnskstr2020(q2)
      return
      end
***

*==============================================================================================
*     LEPTONIC AND TOP CONTRIBUTION TO VACUUM POLARIZATION
** with correct expansion for |q2| <<< m2!!
      FUNCTION SUMMA(AM,Q2,I)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 NC(9),QF2(9)
      common/ncqf2/nc,qf2
      data nc  /1.d0,1.d0,1.d0,3.d0,
     . 3.d0, !u
     . 3.d0, !d
     . 3.d0, !c
     . 3.d0, !s
     . 3.d0 /!b
      data qf2 /1.d0,1.d0,1.d0,0.44444444444444444444d0,
     . 0.44444444444444444444d0, !u
     . 0.11111111111111111111d0, !d
     . 0.44444444444444444444d0, !c
     . 0.11111111111111111111d0, !s
     . 0.11111111111111111111d0 /!b
* NC AND QF ARE COLOR FACTOR (1 FOR LEPTONS,
* 3 FOR TOP) AND CHARGE
      AM2=AM**2

      if (Q2.eq.0.d0) then
         summa = 0.d0
         return
      endif

      z6lim = 1d-8

      IF (Q2.GT.0.D0.AND.Q2.LT.(4.D0*AM2)) THEN
         x = 4.D0*AM2/Q2
         z = sqrt(1.d0/x)
         z2 = z*z
         z4 = z2*z2
         z6 = z2*z4

         if (z6.gt.z6lim) then
            SQ = DSQRT(1d0 - z2)
            s1 = -5.D0/3.D0 - 1.d0/z2
            s2 = (1.d0 + z2- 2.D0*z4)/z/z2*1.d0/SQ*DATAN(z/SQ)
            s1ps2 = s1 + s2
         else
*     with wxmaxima
* s2 = taylor((1+z^2-2*z^4)/sqrt(1-z^2)/z^3*atan(z/sqrt(1-z^2)), z , 0 , 20);
            s1ps2 = -z2*
     .        (4.d0/5.d0+12.d0/35.d0*z2+64.d0/315.d0*z4+32.d0/231.d0*z6)
         endif
         SUMMA=1.d0/3.d0*NC(I)*QF2(I)*s1ps2
      ELSE
         x = -4.D0*AM2/Q2
         if (x.gt.0.d0) then    ! space-like
            z  = sqrt(1.d0/x)
            z2 = z*z
            z4 = z2*z2
            z6 = z2*z4

            if (z6.gt.z6lim) then
               SQ     = DSQRT(1.D0 + z2)
               ARGLOG = DABS( (z - sq)/(z+sq) )
               s1 = -5.D0/3.D0 + 1.d0/z2
               s2 = 0.5d0*(1.d0 - z2 - 2.D0*z4)/SQ*DLOG(ARGLOG)/z2/z

               s1ps2 = s1 + s2
            else
*     with wxmaxima
* s2 = taylor(log( abs ((z - sqrt(1+z^2))/(z + sqrt(1+z^2)) ))/sqrt(1+z^2)*(1-z^2-2*z^4)/2/z^3,z,0,20);
               s1ps2 = z2*
     .        (4.d0/5.d0-12.d0/35.d0*z2+64.d0/315.d0*z4-32.d0/231.d0*z6)
            endif
         else                   ! for x > 0 ie q2 > 4m^2
            x = -x
            z  = sqrt(1.d0/x)
            z2 = z*z
            z4 = z2*z2

            SQ     = DSQRT(z2-1.d0)
            ARGLOG = DABS((z - sq)/(z+sq) )
            s1 = -5.D0/3.D0 - 1.d0/z2
            s2 =  0.5d0*(1.d0 + z2 - 2.D0*z4)/SQ*DLOG(ARGLOG)/z2/z
            s1ps2 = s1 + s2
         endif
         SUMMA = 1.d0/3.d0 * NC(I)*QF2(I)*s1ps2
      ENDIF
      RETURN
      END
*-----
* LEPTONIC AND TOP CONTRIBUTION TO VACUUM POLARIZATION
      FUNCTION SUMMA_OLD(AM,Q2,I)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 NC(9),QF2(9)
      common/ncqf2/nc,qf2
      data nc  /1.d0,1.d0,1.d0,3.d0,
     . 3.d0, !u
     . 3.d0, !d
     . 3.d0, !c
     . 3.d0, !s
     . 3.d0 /!b
      data qf2 /1.d0,1.d0,1.d0,0.44444444444444444444d0,
     . 0.44444444444444444444d0, !u
     . 0.11111111111111111111d0, !d
     . 0.44444444444444444444d0, !c
     . 0.11111111111111111111d0, !s
     . 0.11111111111111111111d0 /!b

* NC AND QF ARE COLOR FACTOR (1 FOR LEPTONS,
* 3 FOR TOP) AND CHARGE
*
      AM2=AM**2
      IF (Q2.GE.0.D0.AND.Q2.LT.(4.D0*AM2)) THEN
         SQ=SQRT(4.D0*AM2/Q2-1.D0)
         SUMMA_OLD=NC(I)*QF2(I)*(-5.D0/9.D0-(4.D0/3.D0)*(AM2/Q2)+
     >        (4.D0/3.D0*(AM2/Q2)**2+1.D0/3.D0*AM2/Q2-1.D0/6.D0)*
     >        4.D0/SQ*ATAN(1.D0/SQ))
        ELSE
           SQ=SQRT(1.D0-4.D0*AM2/Q2)
           ARGLOG=ABS((1.D0-SQ)/(1.D0+SQ))
*
           SUMMA_OLD=NC(I)*QF2(I)*(-5.D0/9.D0-(4.D0/3.D0)*(AM2/Q2)+
     >          (4.D0/3.D0*(AM2/Q2)**2+1.D0/3.D0*AM2/Q2-1.D0/6.D0)*
     >          2.D0/SQ*LOG(ARGLOG))
       ENDIF
*
       RETURN
       END



***********************************************************************
      subroutine svfactor(model,ng,ecms,p1,p2,eps,sv,dsv,q,qphmat)
#ifdef COLLIER
      use collier
#endif
      implicit double precision (a-h,l,o-z)
      dimension p1(0:3),p2(0:3),pin1(0:3),pin2(0:3),p1b(0:3),p2b(0:3)
      dimension p1cut1(0:3),p2cut2(0:3),p1cut(0:3),p2cut(0:3),qcut(0:3)
      dimension q(0:3)
      dimension pin1cl(0:3),pin2cl(0:3),p1cl(0:3),p2cl(0:3),qcl(4)
      dimension pin1cl1(4),pin2cl1(4),p1cl1(4),p2cl1(4),qclu(4,2)

      dimension pin1br(0:3),pin2br(0:3),p1br(0:3),p2br(0:3),qbr(0:3)
      dimension pin1brr(0:3),pin2brr(0:3),p1brr(0:3),p2brr(0:3)
      dimension qbrr(0:3),p1br1(0:3),p2br1(0:3)
      dimension pin1bbr(0:3)
      dimension p1r(0:3),p2r(0:3),q1r(0:3),q2r(0:3),qphr(40,0:3)
      common/momentacheck1/pin1br,pin2br,p1br,p2br,qbr
      common/momentacheck2/p1r,p2r,q1r,q2r,qphr

      integer errflag
      dimension qphmat(40,0:3)
      dimension iclose(40)
      common/iclosest/iclose

      character*10 model
      character*6 ord
      double precision masses(4),charges(4),kp1,kp2,kp3,kp4
      character*2 fs
      common/finalstate/fs
      common/qedORDER/ord
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      common/momentainitial/pin1,pin2
      common/reducedtoborn/p1b,p2b,iref
      common/momentainitialred/pin1b(0:3),pin2b(0:3),pin1r(0:3)
     >     ,pin2r(0:3)
      common/forborncrosssection/phsp2b,flux2b,bornme,bornmeq2,bornmez,
     >     bornmenovp
      common/ialpharunning/iarun
      common/icountsv/icsvf,ifirst
      common/eventlsv/eventlimitersv
      character*3 eventlimitersv

      data ifirst,icsvf /0,0/
      integer irecola
      common/recolamode/irecola
      common/various/beta,betafs
      common/softphotonmass/phmass
      character*100  path
      common/cmnpath/path
      double precision lambda2,mudim,mu2dim,m2min
      common/duvmudim/deltauv,mudim,mu2dim,m2min,pi216m1,ec6
      common/gammamass/lambda2
      integer ioffhvp ,inlovp
      common/offhadrVP_nlovp/ioffhvp,inlovp

      double precision vect1(0:3)
      integer efdd
      common/instabilities/instabilities
      integer instabilities
      common/checkinvariants/lim,switchcollier
      integer switchcollier
      double precision lim

      double precision egmax
      common/egmaxforrecola/egmax

      integer icll
      common/colliermode/icll

      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi

      real*8 virtplussoftpions
      external virtplussoftpions
      double complex elmat1,elmat2

      common/radoption/ifisrfsr
      integer ifisrfsr
      double precision P1oP0,P0ap,argtilde
      double precision lm1,lwmtq,z2,z3,z4

#ifdef LOOPTOOLS
#include "looptools.h"
#endif

      icsvf = icsvf + 1

      dsv = 0.d0
      sv  = 1.d0

      if (ord.eq.'born'.or.ord.eq.'struct') then
         sv = 1.d0
         return
      endif

      if (ifirst.eq.0.and.icll.eq.0) then
         ifirst = 1
         icll = 1
!	NOW DELTAU IS SET AT THE BEGINNING OF RECOLA_INT.F
!         deltauv= 0.d0

         mu2dim= 1.d0
         lambda2= phmass*phmass
         m2min= 0.d0

         pi216m1 = pi*pi*16.d0
         pi216m1 = 1.d0/pi216m1

         ec6 = 4.d0*pi*alpha
         ec6 = ec6*ec6*ec6

#ifdef LOOPTOOLS
         call ltini
         call setdelta(deltauv)
         call setmudim(mu2dim)
         call setlambda(lambda2)
         call setminmass(m2min)
#endif
c$$$#ifdef COLLIER
c$$$         call Init_cll(5,3,path(1:len_trim(path))//'cll_out')
c$$$         call SetMode_cll(1)    ! 1 is COLI, 2 DD, 3 both. 1 is faster
c$$$         call InitCacheSystem_cll(1,5)
c$$$         call SetDeltaUV_cll(deltauv)
c$$$         call SetMuUV2_cll(mu2dim)
c$$$         call SetmuIR2_cll(phmass*phmass)
c$$$c         print*, 'cll phmass',phmass, mu2dim,deltauv
c$$$#endif

      endif

      egmax=0.5d0*sqrt(dot(pin1+pin2,pin1+pin2))*eps

      iarunbck=iarun

      if(ifirst.eq.0.d0.and.icll.eq.1) then

#ifdef COLLIER
      call SwitchOffErrStop_cll()
#endif
      ifirst=1
      endif

***** a_rad
      if(ord.eq.'a_rad')then
         if (ng.ge.2) then
            sv=1d0
            return
         endif
         if (ng.eq.1) then
            if(model.ne.'psrad') then
               if(irecola.gt.0) then
                  do i=0,3
                     q(i)=qphmat(1,i)
                  enddo

*******************CHECK INVARIANTS
!                  call check_invariants(pin1,pin2,p1,p2,q)
**************************
                  if(fs.eq.'mr') then

                     if(ifisrfsr.eq.11) then

                        call virtual_MUR(pin1,pin2,p1,p2,q,
     >                       virtual,qphmat)

                     elseif(ifisrfsr.eq.10) then

                        call virtualISR_mur(pin1,pin2,p1,p2,q,
     >                       virtual,qphmat,'isr')

                     elseif(ifisrfsr.eq.01) then

                        call virtualISR_mur(pin1,pin2,p1,p2,q,
     >                       virtual,qphmat,'fsr')
                     endif



                     call bk_me(pin1,pin2,p1,p2,q,born,qphmat)


                     sv= (virtual/born+1d0
     >                 +1d0*soft_integral_muons(egmax,pin1,pin2,p1,p2))

                  elseif(fs.eq.'er') then
                     iarun=0

                     call recola_ampl(ng,pin1,pin2,p1,p2,q,brcl,
     >                    svrcl, qphmat)

                     iarun=iarunbck

                     sv= (svrcl/brcl)


!     ADDED FOR PION RADIATIVE CHANNEL
                  elseif(fs.eq.'pr') then


                         call soft_integral_pions(egmax,pin1,pin2,p1,p2,
     >                    rsofttot,rsoftISR,rsoftFSR,rsoftIFI)


                          call bk_me(pin1,pin2,p1,p2,q,born,qphmat)

                          if(ifisrfsr.eq.11) then

                             call virtual_pr(pin1,pin2,p1,p2,q,
     >                            virtual,qphmat)

                                  sv =(virtual/born+(1d0+1d0*rsofttot))
!                             sv = virtual/born
                           elseif(ifisrfsr.eq.10) then

                              call virtualISR_pr(pin1,pin2,p1,p2,q,
     >                                    virtual,qphmat)

                              call virtual_CT_pr(pin1,pin2,p1,p2,q,
     >                         virtual1,qphmat)

                              virtual=virtual+virtual1

                              sv =(virtual/born+(1d0+1d0*rsoftISR))

                           elseif(ifisrfsr.eq.01) then

                              call virtualFSR_pr(pin1,pin2,p1,p2,q,
     >                             virtual,qphmat)
                              call virtual_CT_pr(pin1,pin2,p1,p2,q,
     >                         virtual1,qphmat)

                              virtual=virtual+virtual1
                              sv =(virtual/born+(1d0+1d0*rsoftFSR))
                           endif

                        endif
c$$$                        if(eventlimitersv.eq.'unw') then
c$$$                           call new_boost(p1,p1,p1br1,1)
c$$$                           call new_boost(p1,p2,p2br1,1)
c$$$                           beta_rel= dsqrt(1.d0-dot(p2br,p2br)
c$$$     >                          /p2br1(0)**2)
c$$$                           eta_c = 2d0*pi*alpha/(2d0*beta_rel/
c$$$     >                          (1.d0+beta_rel**2))
c$$$
c$$$                           if(eta_c.ge.0.1d0) then
c$$$                              sv = sv/eta_c
c$$$                           endif
c$$$                        endif
                        return
           endif

        elseif(model.eq.'psrad') then
               aieps     = -2.d0*dlog(eps)-1.5d0+2.d0*eps-0.5d0*eps**2
               aiepsscal = -2.d0*dlog(eps)-2.d0+2.d0*eps

!     see below
               dp1p2 = dlog(2.d0*dot(pin1,pin2)/ame**2)
               dp1p3 = dlog(2.d0*dot(pin1,p1)/ame**2)
               dp1p4 = dlog(2.d0*dot(pin1,p2)/ame**2)
               dp2p3 = dlog(2.d0*dot(pin2,p1)/ame**2)
               dp2p4 = dlog(2.d0*dot(pin2,p2)/ame**2)
               dp3p4 = dlog(2.d0*dot(p1,p2)/ame**2)

               if (fs.eq.'pr'.or.fs.eq.'mr'.or.fs.eq.'er') then

                  lcollnew = dp1p2+dp1p3-dp1p4-dp2p3+dp2p4+dp3p4
                  lcollnew = 2.d0 * lcollnew - 4.d0
                  lcollnew = lcollnew/4.d0
                  arg = 2.d0*alpha/pi * lcollnew * aieps
****  NEW ***************
                  charges(1) = -1.d0
                  charges(2) = -1.d0
                  charges(3) = -1.d0
                  charges(4) = -1.d0
                  masses(1)  = ame
                  masses(2)  = ame
                  masses(3)  = ame
                  masses(4)  = ame

                  if (fs.eq.'mr') then
                     masses(3)  = ammu
                     masses(4)  = ammu

                     if(ifisrfsr.eq.10) then
		       charges(3) = 0.d0
		       charges(4) = 0.d0
		     elseif(ifisrfsr.eq.01) then
		       charges(1) = 0.d0
		       charges(2) = 0.d0
		     endif
                  endif
                  if (fs.eq.'pr') then
                     masses(3)  = ampi
                     masses(4)  = ampi

                     if(ifisrfsr.eq.10) then
		       charges(3) = 0.d0
		       charges(4) = 0.d0
		     elseif(ifisrfsr.eq.01) then
		       charges(1) = 0.d0
		       charges(2) = 0.d0
		     endif
                  endif
                  lcoll1=eikonalintegral(pin1,pin2,p1,p2,charges,masses)

                  arg2  = 2.d0*alpha/pi * lcoll1 * aieps


                  if (fs.eq.'pr') then
                     arg2 = alpha/pi*lcoll1*(aieps+aiepsscal)
                  endif
                  arg = arg2
               endif

               sv = 1.d0 - arg

c$$$*-----------------COULOMB CORRECTION AT THRESHOLD TO UNWEIGHT
c$$$               if(eventlimitersv.eq.'unw') then
c$$$                  call new_boost(p1,p1,p1br1,1)
c$$$                  call new_boost(p1,p2,p2br1,1)
c$$$                  beta_rel= dsqrt(1.d0-dot(p2br,p2br)/p2br1(0)**2)
c$$$                  eta_c = 2d0*pi*alpha/(2d0*beta_rel/(1.d0+beta_rel**2))
c$$$
c$$$                  if(eta_c.ge.0.1d0) then
c$$$                     sv = sv*eta_c
c$$$                  endif
c$$$               endif
*----------------

            endif
            return
         endif
      endif

****
**** exprad
****
      if (ord.eq.'exprad') then


****************TRY TO TAKE THE PHOTON I CALCULATE THE VIRT CORRECTIONS ON TO BE THE HARDEST (LBK THEOREM)
c$$$         if(ng.gt.3) then
c$$$            sv=0.d0
c$$$            deltasv=0.d0
c$$$            return
c$$$         endif
c$$$
c$$$         if(ng.eq.3) then
c$$$            sv=1.d0
c$$$            deltasv=0.d0
c$$$            return
c$$$         endif
                ind=1
                call getsignalphot(pin1,pin2,p1,p2,qphmat,ng,ind,j)
                pin1cl=pin1
                pin2cl=pin2

                p1cl=p1
                p2cl=p2

                do k = 1,ng
                   if(k.ne.j) then
                      if (iclose(k).eq.1) then
                         pin1cl=pin1cl-qphmat(k,:)
                      endif
                      if (iclose(k).eq.2) then
                         pin2cl=pin2cl-qphmat(k,:)
                      endif
                      if (iclose(k).eq.3) then
                         p1cl=p1cl+qphmat(k,:)
                      endif
                      if (iclose(k).eq.4) then
                         p2cl=p2cl+qphmat(k,:)
                      endif
                   endif
                enddo

                ind=1
                call maptoone_rad(pin1cl,pin2cl,p1cl,p2cl,qphmat(j,:),
     >           iclose(j),ind,pin1br,pin2br,p1br,p2br,qbr)

         aieps     = -2.d0*dlog(eps)-1.5d0+2.d0*eps-0.5d0*eps**2

         aiepsscal = -2.d0*dlog(eps)-2.d0+2.d0*eps

! with reduced to born momenta...
         Aopi = alpha/pi
         s = 4.d0*pin1br(0)*pin1br(0)

         p1mod = sqrt(tridot(p1br,p1br))
         c = p1br(3)/p1mod
         t = -s/2.d0*(1.d0 - c)
         u = -s/2.d0*(1.d0 + c)

         ame2mu = 1.d0/ame/ame
         sdif = bornme*phsp2b/flux2b

c$$$         s =dot(pin1br+pin2br,pin1br+pin2br)
c$$$         aieps   =-2.d0*dlog(eps*dsqrt(s)*0.5d0/phmass)
c$$$         aiepsscal=aieps


*     * VP STUFF
*** changed for no VP in 1-loop
         sdif = bornmenovp*phsp2b/flux2b
***
         if (fs.eq.'pr'.or.fs.eq.'mr'.or.fs.eq.'er') then

**** NEW ***************
            charges(1) = -1.d0
            charges(2) = -1.d0
            charges(3) = -1.d0
            charges(4) = -1.d0
            masses(1)  = ame
            masses(2)  = ame
            masses(3)  = ame
            masses(4)  = ame
            if (fs.eq.'mr') then
              if(ifisrfsr.eq.10) then
		    charges(3) = 0d0
		    charges(4) = 0d0
	      elseif(ifisrfsr.eq.01) then
		    charges(1) = 0d0
		    charges(2) = 0d0
	      endif
               masses(3)  = ammu
               masses(4)  = ammu
            endif
            if (fs.eq.'pr') then
              if(ifisrfsr.eq.10) then
		    charges(3) = 0d0
		    charges(4) = 0d0
	      elseif(ifisrfsr.eq.01) then
		    charges(1) = 0d0
		    charges(2) = 0d0
	      endif
               masses(3)  = ampi
               masses(4)  = ampi
            endif
            lcoll1 = eikonalintegral(pin1br,pin2br,
     .           p1br,p2br,charges,masses)
            arg2   = 2.d0*aopi * lcoll1 * aieps

******ADDING THE MISSING COLLINEAR TERM ALSO IN THE SUDAKOV
c$$$
c$$$            arg2   = 2.d0*aopi *(lcoll1*aieps +(0.5d0-eps**2/2.d0))
*********
            if (fs.eq.'pr') then
               arg2 = alpha/pi*lcoll1*(aieps+aiepsscal)
            endif
            arg = arg2

**************
            if (model.ne.'psrad') then
               virtual =0.d0
               egmaxbck=egmax
               egmax=0.5d0*dsqrt(dot(pin1br+pin2br,pin1br+pin2br))*eps
               if(fs.eq.'mr') then


   		    if(ifisrfsr.eq.11) then

                    call virtual_MUR(pin1br,pin2br,p1br,p2br,qbr,
     >                                       virtual,qphmat)

   		    elseif(ifisrfsr.eq.10) then

                    call virtualISR_mur(pin1br,pin2br,p1br,p2br,qbr,
     >                                       virtual,qphmat,'isr')

     	            elseif(ifisrfsr.eq.01) then

                    call virtualISR_mur(pin1br,pin2br,p1br,p2br,qbr,
     >                                       virtual,qphmat,'fsr')
     	 	    endif



   		    call bk_me(pin1br,pin2br,p1br,p2br,qbr,born,qphmat)

                    svpureoal= (virtual/born+
     >     1d0*soft_integral_muons(egmax,pin1br,pin2br,p1br,p2br))

!!!!  Might be interesting to add P(1)*born
!     x is the fraction of energy retained by the electron
!     only IS for the moment
c$$$                    if(ng.eq.2) then
c$$$                       x = (pin1(0) - qphmat(j,0))/pin1(0)
c$$$                       P0ap = (1.d0 + x**2)/(1-x)
c$$$
c$$$                       P1oP0 = -alpha/2.d0/pi*(2.d0*log(x)*log(1.d0-x)+
c$$$     >                      1.5d0*log(x) +((1.5d0+3.5d0*x)*log(x) +
c$$$     >                      0.5d0*(1.d0+x)*log(x)**2 +
c$$$     >                      5.d0*(1.d0-x))/P0ap)
c$$$                       svpureoal = svpureoal +P1oP0
c$$$                       
c$$$                    endif
                 elseif(fs.eq.'er') then
                    iarunbck=iarun
                    iarun=0

		    call recola_ampl(1,pin1br,pin2br,p1br,p2br,qbr,
     >                   brcl,svrcl, qphmat)

		    iarun=iarunbck

                    svpureoal= (svrcl/brcl)-1d0


!	ADDED FOR PION RADIATIVE CHANNEL
                  elseif(fs.eq.'pr') then


     		call soft_integral_pions(egmax,pin1br,pin2br,p1br,
     >     		    p2br,rsofttot,rsoftISR,rsoftFSR,rsoftIFI)

		  call bk_me(pin1br,pin2br,p1br,p2br,qbr,born,qphmat)
   		if(ifisrfsr.eq.11) then

		  call virtual_pr(pin1br,pin2br,p1br,p2br,qbr,
     >                                       virtual,qphmat)

     		  svpureoal =(virtual/born+(1d0*rsofttot))

     		elseif(ifisrfsr.eq.10) then

		  call virtualISR_pr(pin1br,pin2br,p1br,p2br,qbr,
     >                                       virtual,qphmat)

     		  call virtual_CT_pr(pin1br,pin2br,p1br,p2br,qbr,
     >		  virtual1,qphmat)

                  virtual=virtual+virtual1
     		  svpureoal =(virtual/born+(1d0*rsoftISR))
     		elseif(ifisrfsr.eq.01) then

		  call virtualFSR_pr(pin1br,pin2br,p1br,p2br,qbr,
     >                                       virtual,qphmat)
     		  call virtual_CT_pr(pin1br,pin2br,p1br,p2br,qbr,
     > 		  virtual1,qphmat)

                  virtual=virtual+virtual1



     		  svpureoal =(virtual/born+(1d0*rsoftFSR))
                 endif
                 endif
                 egmax=egmaxbck
              else

            svpureoal = -arg

         endif

!         if(j1.ne.j) print*,'arg sud',arg,svpureoal
!         eps=eps*10000.d0
         deltasv= svpureoal + arg !!! NB

*-----------------COULOMB CORRECTION AT THRESHOLD
         call new_boost(p1br,p1br,p1br1,1)
         call new_boost(p1br,p2br,p2br1,1)
         beta_rel= dsqrt(1.d0-dot(p2br,p2br)/p2br1(0)**2)
         eta_c = 2d0*pi*alpha/(2d0*beta_rel/(1.d0+beta_rel**2))

         if((eta_c.ge.0.1d0).and.ifisrfsr.ne.10) then
            deltasv=deltasv-eta_c
         endif
*----------------


***   NEW ***************
            charges(1) = -1.d0
            charges(2) = -1.d0
            charges(3) = -1.d0
            charges(4) = -1.d0
            masses(1)  = ame
            masses(2)  = ame
            masses(3)  = ame
            masses(4)  = ame
            if (fs.eq.'mr') then
              if(ifisrfsr.eq.10) then
		    charges(3) = 0d0
		    charges(4) = 0d0
	      elseif(ifisrfsr.eq.01) then
		    charges(1) = 0d0
		    charges(2) = 0d0
	      endif

               masses(3)  = ammu
               masses(4)  = ammu
            endif
            if (fs.eq.'pr') then
              if(ifisrfsr.eq.10) then
		    charges(3) = 0d0
		    charges(4) = 0d0
	      elseif(ifisrfsr.eq.01) then
		    charges(1) = 0d0
		    charges(2) = 0d0
	      endif
               masses(3)  = ampi
               masses(4)  = ampi
            endif

               lcoll1=eikonalintegral(pin1,pin2,
     .              p1,p2,charges,masses)
               arg2   = 2.d0*aopi * lcoll1 * aieps
               if (fs.eq.'pr') then
                  arg2 = alpha/pi*lcoll1*(aieps+aiepsscal)
               endif
               arg = arg2

            endif
         ffs = dexp(-arg)

*------------COULOMB CORRECTION AT THRESHOLD
         if(eta_c.ge.0.1d0.and.ifisrfsr.ne.10) then
            ffs =ffs*(1.d0+eta_c)
         endif
*------------------

         if (model.eq.'psrad') then

            deltasv = 0.d0
         endif
         sv = ffs * (1.d0 + deltasv)

!         sv = ffs * deltasv
c$$$         if(ng.eq.1) then
c$$$            sv=(1.d0-arg+arg**2/2.d0)*(1.d0+deltasv)-arg**2/2.d0*deltasv
c$$$c$$$!!!   Exact NLO
c$$$c$$$            sv = 1.d0 +deltasv -arg
c$$$c$$$!!!   Exact NNLO by https://arxiv.org/pdf/2308.09479
c$$$c$$$            sv = sv - arg**2/2.d0 !piece that goes as the square of the universal correction
c$$$c$$$            sv = sv - svpureoal*arg
c$$$c$$$            sv = sv + 2.d0*aopi*lcoll1/16.d0*(48.d0*1.202056903159594d0-
c$$$c$$$     >           24.d0*pi**2/6.d0+3.d0)*aopi
c$$$
c$$$            lm1   = log(dot(pin1+pin2,pin1+pin2)/ame**2)-1.d0
c$$$            lwmtq = -log(eps)-3.d0/4.d0
c$$$            z2 = pi**2/6.d0
c$$$            z3 = 1.202056903159594d0
c$$$            z4 = pi**4/90.d0
c$$$            sv  = aopi**2/16.d0*(-32.d0*lm1**2*lwmtq**2+
c$$$     >           16.d0*lm1*lwmtq*(7.d0-3.d0*z2)
c$$$     >           +(48.d0*z3-24.d0*z2+3.d0)*(lm1+1.d0)
c$$$     >           -129.d0*z4-12.d0*z3-96.d0*z2*log(2.d0)+
c$$$     >           128.d0*z2-215.d0/4.d0) 
c$$$
c$$$            sv = sv+svpureoal*aopi/4.d0*(-8.d0*lm1*lwmtq-6.d0*z2+14.d0)
c$$$            sv = sv +svpureoal+1.d0
c$$$         elseif(ng.eq.2) then
c$$$            sv=((1.d0-arg)*(1.d0+deltasv)+arg*deltasv)
c$$$            stilde=dot(pin1br+pin2br,pin1br+pin2br)
c$$$            s=dot(pin1+pin2,pin1+pin2)
c$$$            sv = sv + aopi*1.5d0*log(stilde/s)
c$$$         else
c$$$            sv=1.d0
c$$$         endif


         dsv = deltasv
         return
      endif


      if (ord.eq.'alpha') then

         if (ng.ge.1) then
            sv = 1.d0
            return
         endif
         if (ng.eq.0) then

            s = ecms**2
            c = p1(3)/dsqrt(tridot(p1,p1))
            t = -s/2.d0*(1.d0 - c)

            sigma0  = bornme*phsp2b/flux2b
            sigmaq2 = bornmeq2*phsp2b/flux2b
**
            if (model.ne.'ps') then
**
            if (fs.eq.'gg') then
!               if (irecola.eq.0) then

                  BORNBK=2.D0*PI*ALPHA*ALPHA/S*BETA*BETA/
     >                 (1.D0+BETA*C)/(1.D0-BETA*C)*(1.D0+C*C)
                  BORNBK = BORNBK*S*128d0*pi

                  svgg = aintsvgg(s,c,eps)/bornbk*bornme ! to restore full mass effects in the IR

                  sv = svgg/bornme/phsp2b*flux2b

c$$$               else
c$$$                  call recola_ampl(0,pin1,pin2,p1,p2,p2,brcl,svrcl
c$$$     .                           ,qphmat)
c$$$                  sv=svrcl/brcl
c$$$                  stop
c$$$               endif
            elseif(fs.eq.'pp'.or.fs.eq.'pi') then
	       if(iffpi.le.2) then

                  sv = aintsv_pions(pin1,pin2,p1,p2,eps)/sigma0

               elseif(iffpi.eq.3) then
                   sv = 1d0+virtplussoftpions(pin1,pin2,p1,p2,eps)
                endif
               smtob = 0.d0
               if (inlovp.eq.1) then
                  call get_twiceVPtimesTreeOverBorn(pin1,pin2,p1,p2,
     .                 smtob)
               endif
               sv = sv + smtob
               dsv= sv-1.d0
            else

               sv = aintsv_new(pin1,pin2,p1,p2,eps)/sigma0
c$$$               call recola_ampl(0,pin1,pin2,p1,p2,p2,brcl,svrcl
c$$$     .              ,qphmat)
c$$$               sv=svrcl/brcl
               smtob = 0.d0
               if (inlovp.eq.1) then
                  call get_twiceVPtimesTreeOverBorn(pin1,pin2,p1,p2,
     .                 smtob)
               endif
                  sv =  sv + smtob
               endif


               elseif (model.eq.'ps') then
               aieps     = -2.d0*dlog(eps)-1.5d0+2.d0*eps-0.5d0*eps**2
               aiepsscal = -2.d0*dlog(eps)-2.d0+2.d0*eps
! see below
               dp1p2 = dlog(2.d0*dot(pin1,pin2)/ame**2)
               dp1p3 = dlog(2.d0*dot(pin1,p1)/ame**2)
               dp1p4 = dlog(2.d0*dot(pin1,p2)/ame**2)
               dp2p3 = dlog(2.d0*dot(pin2,p1)/ame**2)
               dp2p4 = dlog(2.d0*dot(pin2,p2)/ame**2)
               dp3p4 = dlog(2.d0*dot(p1,p2)/ame**2)

               if (fs.eq.'ee'.or.fs.eq.'mm'
     +              .or.fs.eq.'pp'.or.fs.eq.'pi') then
                  lcollnew = dp1p2+dp1p3-dp1p4-dp2p3+dp2p4+dp3p4
                  lcollnew = 2.d0 * lcollnew - 4.d0
                  lcollnew = lcollnew/4.d0
                  arg = 2.d0*alpha/pi * lcollnew * aieps
**** NEW ***************
                  charges(1) = -1.d0
                  charges(2) = -1.d0
                  charges(3) = -1.d0
                  charges(4) = -1.d0
                  masses(1)  = ame
                  masses(2)  = ame
                  masses(3)  = ame
                  masses(4)  = ame
                  if (fs.eq.'mm') then
                     masses(3)  = ammu
                     masses(4)  = ammu
                  endif
                  if (fs.eq.'pp') then
                     masses(3)  = ampi
                     masses(4)  = ampi
                  endif
                  if (fs.eq.'pi') then
                     charges(3) = 0.d0
                     charges(4) = 0.d0
                     masses(3)  = ampi
                     masses(4)  = ampi
                  endif
                  lcoll1=eikonalintegral(pin1,pin2,p1,p2,charges,masses)
                  arg2  = 2.d0*alpha/pi * lcoll1 * aieps
                  if (fs.eq.'pp') then
                     arg2 = alpha/pi*lcoll1*(aieps+aiepsscal)
                  endif
                  arg = arg2
**************************
               elseif (fs.eq.'gg') then
                  lcollnew = dp1p2
                  lcollnew = lcollnew - 1.d0
                  arg = alpha/pi * lcollnew * aieps
               endif
               sv = 1.d0 - arg
            endif
            return
         endif
      endif

      if (ord.eq.'exp') then
         aieps     = -2.d0*dlog(eps)-1.5d0+2.d0*eps-0.5d0*eps**2
         aiepsscal = -2.d0*dlog(eps)-2.d0+2.d0*eps
! with reduced to born momenta...
         aopi = alpha/pi
         s = 4.d0*pin1b(0)*pin1b(0)

! for gg not-reduced center of mass..
ccccc         if (fs.eq.'gg') s = ecms**2

         p1mod = sqrt(tridot(p1b,p1b))
         c = p1b(3)/p1mod
         t = -s/2.d0*(1.d0 - c)
         u = -s/2.d0*(1.d0 + c)

         ame2mu = 1.d0/ame/ame
         sdif = bornme*phsp2b/flux2b

** VP STUFF
*** changed for no VP in 1-loop
         sdif = bornmenovp*phsp2b/flux2b
***
         if (fs.eq.'ee'.or.fs.eq.'mm'
     .        .or.fs.eq.'pp'.or.fs.eq.'pi') then
            scale = s*t/u
            lcoll = dlog(scale*ame2mu) - 1.d0
            arg = 2.d0 * aopi * lcoll * aieps

**** NEW ***************
            charges(1) = -1.d0
            charges(2) = -1.d0
            charges(3) = -1.d0
            charges(4) = -1.d0
            masses(1)  = ame
            masses(2)  = ame
            masses(3)  = ame
            masses(4)  = ame
            if (fs.eq.'mm') then
               masses(3)  = ammu
               masses(4)  = ammu
            endif
            if (fs.eq.'pp') then
               masses(3)  = ampi
               masses(4)  = ampi
            endif
            if (fs.eq.'pi') then
               charges(3) = 0.d0
               charges(4) = 0.d0
               masses(3)  = ampi
               masses(4)  = ampi
            endif
            lcoll1 = eikonalintegral(pin1b,pin2b,p1b,p2b,charges,masses)
            arg2   = 2.d0*aopi * lcoll1 * aieps
            if (fs.eq.'pp') then
               arg2 = alpha/pi*lcoll1*(aieps+aiepsscal)
            endif
            arg = arg2

**************************
** VP STUFF
ccc   changed for no VP in virtual
   	    iarunbck = iarun ! added
            iarun = 0           ! added ! ON OFF
            if (model.ne.'ps') then
               if (fs.eq.'pp'.or.fs.eq.'pi') then
                 if(iffpi.le.2) then
                 svpureoal=aintsv_pions(pin1b,pin2b,p1b,p2b,eps)/sdif
     .                 -1.d0
                 elseif(iffpi.eq.3) then
                    svpureoal=virtplussoftpions(pin1b,pin2b,p1b,p2b,eps)
                 endif
               else
                  svpureoal = aintsv_new(pin1b,pin2b,p1b,p2b,eps)/sdif
     .                 -1.d0
               endif
            else
            svpureoal = -arg
            endif
            iarun = iarunbck    ! added
         elseif (fs.eq.'gg') then
            scale = s
            lcoll = dlog(scale*ame2mu) - 1.d0
            arg = 1.d0 * aopi * lcoll * aieps
            svpureoal = aintsvgg(s,c,eps)/sdif - 1.d0
         endif
         deltasv = svpureoal + arg !!! NB

*** NEW SFF, formula 2.23 PhD. Thesis. in the limit dot(pi,pj) >> m**2
         dp1p2 = dlog(2.d0*dot(pin1,pin2)*ame2mu)
         dp1p3 = dlog(2.d0*dot(pin1,p1)*ame2mu)
         dp1p4 = dlog(2.d0*dot(pin1,p2)*ame2mu)
         dp2p3 = dlog(2.d0*dot(pin2,p1)*ame2mu)
         dp2p4 = dlog(2.d0*dot(pin2,p2)*ame2mu)
         dp3p4 = dlog(2.d0*dot(p1,p2)*ame2mu)

cc         dp1p2 = dlog(2.d0*dot(pin1b,pin2b)*ame2mu)
cc         dp1p3 = dlog(2.d0*dot(pin1b,p1b)*ame2mu)
cc         dp1p4 = dlog(2.d0*dot(pin1b,p2b)*ame2mu)
cc         dp2p3 = dlog(2.d0*dot(pin2b,p1b)*ame2mu)
cc         dp2p4 = dlog(2.d0*dot(pin2b,p2b)*ame2mu)
cc         dp3p4 = dlog(2.d0*dot(p1b,p2b)*ame2mu)

         argbck = arg

         if (fs.eq.'ee'.or.fs.eq.'mm'
     .        .or.fs.eq.'pp'.or.fs.eq.'pi') then
            lcollnew = dp1p2+dp1p3-dp1p4-dp2p3+dp2p4+dp3p4
            lcollnew = 2.d0 * lcollnew - 4.d0
            lcollnew = lcollnew/4.d0
            arg = 2.d0 * aopi * lcollnew * aieps

**** NEW ***************
            charges(1) = -1.d0
            charges(2) = -1.d0
            charges(3) = -1.d0
            charges(4) = -1.d0
            masses(1)  = ame
            masses(2)  = ame
            masses(3)  = ame
            masses(4)  = ame
            if (fs.eq.'mm') then
               masses(3)  = ammu
               masses(4)  = ammu
            endif
            if (fs.eq.'pp') then
               masses(3)  = ampi
               masses(4)  = ampi
            endif
            if (fs.eq.'pi') then
               charges(3) = 0.d0
               charges(4) = 0.d0
               masses(3)  = ampi
               masses(4)  = ampi
            endif
            lcoll1 = eikonalintegral(pin1,pin2,p1,p2,charges,masses)
            arg2   = 2.d0*aopi * lcoll1 * aieps
            if (fs.eq.'pp') then
               arg2 = alpha/pi*lcoll1*(aieps+aiepsscal)
            endif
            arg = arg2
**************************
         elseif (fs.eq.'gg') then
            lcollnew = dp1p2
            lcollnew = lcollnew - 1.d0
            arg = 1.d0 * aopi * lcollnew * aieps
         endif

!         deltasv = svpureoal + arg      !!! NB
         ffs = dexp(-arg)
         if (model.eq.'ps'.or.model.eq.'eikonal') then
            deltasv = 0.d0
         endif
         sv = ffs * (1.d0 + deltasv)
********TRYYYY
c$$$
c$$$         sv= 1-arg+deltasv+arg**2/2.d0-arg*deltasv
c$$$
c$$$  print*,-arg+deltasv+arg**2/2.d0-arg*deltasv-
c$$$     >        aopi*(2.d0*dlog(eps)*(dlog(s/ame**2)-1))-
c$$$     >        (aopi**2)*(2.d0*(dlog(eps)*(dlog(s/ame**2)-1))**2+
c$$$     >     2.d0*dlog(eps)*(dlog(s/ame**2)-1)*(3.d0/2.d0*dlog(s/ame**2)
c$$$     >        +2.d0*pi**2/6.d0-2.d0))
c$$$  print*,deltasv*arg
c$$$         s=dot(pin1+pin2,pin1+pin2)
c$$$         print*,2.d0*dlog(eps)*(dlog(s/ame**2)-1)*deltasv*aopi
c$$$     >    +3.d0/2.d0*(dlog(s/ame**2)-1)*deltasv*aopi
************FINISH
         dsv = deltasv
         return
      endif

      if (ord.eq.'alpha2') then
         if (fs.ne.'ee') then
            print*,'alpha2 only for bhabaha! stopping...'
            stop
         endif

         aieps = -2.d0*dlog(eps) -1.5d0+2.d0*eps-0.5d0*eps**2
         s = (pin1b(0) + pin2b(0))**2
         p1mod = sqrt(tridot(p1b,p1b))
         c = p1b(3)/p1mod

         t = -s/2.d0*(1.d0 - c)
         u = -s/2.d0*(1.d0 + c)
c         scale = s*t/u
c         lcoll = dlog(scale/ame**2) - 1.d0
c         arg = 2.d0*alpha/pi * lcoll * aieps

         sdif = bornme*phsp2b/flux2b

         svpureoal = aintsv(s,t,eps)/sdif - 1.d0
c        svpureoal = aintsv_new(pin1b,pin2b,p1b,p2b,eps)/sdif - 1.d0
cc         xxx = aintsv(s,t,eps)
c     deltasv = svpureoal + arg      !!! NB
*** NEW SFF, formula 2.23 PhD. Thesis. in the limit dot(pi,pj) >> m**2
c         dp1p2 = dlog(2.d0*dot(pin1,pin2)/ame**2)
c         dp1p3 = dlog(2.d0*dot(pin1,p1)/ame**2)
c         dp1p4 = dlog(2.d0*dot(pin1,p2)/ame**2)
c         dp2p3 = dlog(2.d0*dot(pin2,p1)/ame**2)
c         dp2p4 = dlog(2.d0*dot(pin2,p2)/ame**2)
c         dp3p4 = dlog(2.d0*dot(p1,p2)/ame**2)
c         lcollnew = dp1p2 + dp1p3 - dp1p4 - dp2p3 + dp2p4 + dp3p4
c         lcollnew = 2.d0 * lcollnew - 4.d0
c         lcollnew = lcollnew/4.d0
c         arg = 2.d0*alpha/pi * lcollnew * aieps
**** NEW ***************
         charges(1) = -1.d0
         charges(2) = -1.d0
         charges(3) = -1.d0
         charges(4) = -1.d0
         masses(1)  = ame
         masses(2)  = ame
         masses(3)  = ame
         masses(4)  = ame
         if (fs.eq.'mm'.or.fs.eq.'mr') then
            masses(3)  = ammu
            masses(4)  = ammu
         endif
         lcoll1 = eikonalintegral(pin1b,pin2b,p1b,p2b,charges,masses)
         arg2   = 2.d0*alpha/pi * lcoll1 * aieps
         arg = arg2
         deltasv = svpureoal + arg !!! NB
!     I used arg with reduced momenta to define deltasv, now I recalculate it
! with true external momenta...
c         lcoll1 = eikonalintegral(pin1,pin2,p1,p2,charges,masses)
c         arg   = 2.d0*alpha/pi * lcoll1 * aieps
**************************
c         ffs = dexp(-arg)
         if (model.eq.'ps') deltasv = 0.d0

         if (ng.le.2) then
            sv = 1.d0
         endif
         if (ng.le.1) then
            sv = sv - arg + deltasv
         endif
         if (ng.eq.0) then
            sv = sv + arg**2*0.5d0 - arg*deltasv
         endif
         return
      endif
      return
      end
******
* COMPLEX VACUUM POLARIZATION, T. Teubner & al.
      FUNCTION VPOLChlmnt_BACKUP(Q2)
      implicit double precision (A-H,L,O-Z)
      parameter (np  = 7186)
c      parameter (np=6883)
      parameter (nw  = 1582)

c      parameter (npt = np) ! or
c      parameter (npt = nw) ! or
      parameter (npt = nw+np)

      dimension s(0:np-1)
      dimension st(0:npt-1)
      dimension pres(0:np-1),pre(0:np-1),pim(0:np-1)
      dimension epres(0:np-1),epre(0:np-1),epim(0:np-1),cv(0:np-1)

      dimension iswvperr(3),iswvperrt(6)
      common/iswitchvperror/iswvperr,iswvperrt

      double complex vpolcnsk,vpolchlmnt,deltahad,im,res
      double complex cP,ressp,restm,ctmp
      external vpolcnsk
      common/vpcnsk/s,pres,pre,pim,epres,epre,epim,cv,im,ifirstnsk

      common/rescaledhaderr/scal_err_dhad

      parameter (nres=15)
      double precision AM(nres),GT(nres),GE(nres)
      common/resonancesparams/am,gt,ge
      double precision EAM(nres),EGT(nres),EGE(nres)
      common/resonancesparamserrors/eam,egt,ege

      dimension iresbh(8) ! the indexes of the resonances to be added by hand (Thomas said to include **only** 1-6)
      dimension prest(0:npt-1),pret(0:npt-1),pimt(0:npt-1)
      dimension eprest(0:npt-1),epret(0:npt-1),epimt(0:npt-1)
      common/vpcteubBCK/st,prest,pret,pimt,eprest,epret,epimt,iresbh,
     .     ifirst,imidguess
      data ifirst/0/

      double precision wt(nw)
      common/wtabulation/wt

      double precision wtot(npt)
      integer indexes(npt)
      common/wtotcommon/wtot

      common/tabulationlimits/wmin,wmax,ic

      common/parameters/ame,ammu,convfac,alpha,pi
      common/ecms/ecms,ecmsnom

      integer iloc(1)

*     using vp_hlmnt_v2_1.f, sent by T. Teubner, private communication

      if (ifirst.eq.0) then
***   filling the common by calling vpolcnsk once (I need the array s --> st[eubner])
         vpolchlmnt = vpolcnsk(q2)
****  filling the common wtabulation
         wtot = 0.d0
         call fillwtabcommon
         wmin = wt(1)
         wmax = wt(nw)
         wtot(1:nw)     = wt**2
         wtot(nw+1:npt) = s
         call sort_vector_inverse(wtot,wtot,indexes,npt)
c         do k = 1,npt
c            print*,k,sqrt(wtot(k)),indexes(k)
c         enddo
c         stop
****
c
c I solve for P, so that is the same quantity returned by NSK
c     I also use the same evaluation points as NSK, and then linear
c     interpolation between them
c     1/(1-P) = 1/(1-danr) + res ->
c     1-P = (1-danr)/(1+(1-danr)*res)
c       P = 1 - (1-danr)/(1+(1-danr)*res)
c         = (1+(1-danr)*res - (1-danr))/(1+(1-danr)*res)
c         = ( danr + res - danr*res) / (1 + res- danr*res)

******************
* filling st, from the nsk data or Teubner tabulation points
c         st(0:np-1) = s   ! or
c         do k = 1,nw
c            st(k-1) = wt(k)*wt(k)
c         enddo ! or
         st = wtot
******************
         prest  = 0.d0
         pret   = 0.d0
         pimt   = 0.d0
         eprest = 0.d0
         epret  = 0.d0
         epimt  = 0.d0
         iresbh(1) = 4
         iresbh(2) = 5
         iresbh(3) = 10
         iresbh(4) = 11
         iresbh(5) = 12
         iresbh(6) = 13
         iresbh(7) = 14
         iresbh(8) = 15
         print*,'Evaluating Teubner et al. VP array '
         open(10,file='cfrvp-0-v2-3',status='unknown')

         do k = 0,npt-1
            energy = sqrt(st(k))

c            hlmnt is tabulated up to wmax and from wmin
            if (energy.le.wmax.and.energy.ge.wmin) then
               ressp = 0.d0
               restm = 0.d0
               do i=1,6
                  j = iresbh(i)
                  AMC = AM(J) !- EAM(J)
                  GEC = GE(J) !- EGE(J)
                  GTC = GT(J) !- EGT(J)
                  V   = 3.D0*GEC/AMC / ALPHA
                  restm = restm +
     .                 V * st(k) / ( st(k)-AMC*AMC+im*GTC*AMC )
                  ressp = ressp -
     .                 V * st(k) / (-st(k)-AMC*AMC+im*GTC*AMC )
               enddo

               print*,'HERE, CHECK WHATs done below...'

               nrflag = 1
c               print*,'attenzione'
               call vphlmntv2(energy,
     &              vprehadsp,vprehadtm,vpimhad,vprelepsp,
     &              vpreleptm,vpimlep,
     &              vpretopsp,vpretoptm,nrflag,evprehadsp,
     &              evprehadtm,evpimhad)

c               vprehadsp = vprehadsp - evprehadsp
c               vprehadtm = vprehadtm - evprehadtm
c               vpimhad   = vpimhad -  evpimhad
!     ! space like
               prest(k) = vprehadsp + vprelepsp + vpretopsp
               prest(k) = prest(k)  ! + ressp ! in the space-like always resummed!
c               prest(k) = vprehadsp + ressp
!     ! time like
               nrflag = 0
               call vphlmntv2(energy,
     &        vprehadsp,vprehadtm,vpimhad,vprelepsp,vpreleptm,vpimlep,
     &        vpretopsp,vpretoptm,nrflag,evprehadsp,evprehadtm,evpimhad)

               pret(k)  = vprehadtm + vpreleptm + vpretoptm
               pimt(k)  = vpimhad   + vpimlep

c   this is 1/(1-P)
               ctmp = 1.d0/(1.d0 - pret(k)-im*pimt(k)) + restm ! in the time-like not resummed!
c   ctmp = 1/(1-P) -> 1-P = 1/ctmp  -> P = 1-1/ctmp
               pret(k) = DReal(1.d0-1.d0/ctmp)
               pimt(k) = DImag(1.d0-1.d0/ctmp)

c               write(10,*)energy,pret(k),pimt(k),pre(k),pim(k),
c     .              prest(k),pres(k)
               write(10,*)energy,pret(k),pimt(k),prest(k)

            endif
         enddo
         close(10)
         print*,'done!'
         print*,'stopping line ~420 sv.f'
         stop
*** guessing where to start the binary search....
         k = 0
         aq2 = ecmsnom**2
         do while(aq2.gt.st(k))
            k = k + 1
         enddo
         k = k - 1
         imidguess = k ! this is correct for born s-channel and no energy spread
***************
         ic = 0
         ifirst=1
      endif

c      k   = 0
c      sk  = 0.d0
c      aq2 = abs(q2)
c      do while(aq2.gt.st(k))
c         k = k + 1
c      enddo
c      k = k - 1
******** binary search ****
      aq2     = abs(q2)
      istart  = 0
      ifinish = npt-1
c      imid = (istart + ifinish)/2
      imid = imidguess

c     do while(ifinish.ge.istart)
      do while(.not.(st(imid+1).gt.aq2.and.st(imid).le.aq2))
         if (aq2.gt.st(imid)) then
            istart  = imid + 1
         else
            ifinish = imid - 1
         endif
         imid = (istart + ifinish)/2
      enddo
      k = imid
******************

      if (q2.gt.0.d0) then
         vpre = pret(k)+(pret(k+1)-pret(k))/(st(k+1)-st(k)) * (q2-st(k))
         vpim = pimt(k)+(pimt(k+1)-pimt(k))/(st(k+1)-st(k)) * (q2-st(k))

c         evpre = epret(k)+ (epret(k+1)-epret(k))/(st(k+1)-st(k))*(q2-st(k))
c         evpim = epimt(k)+ (epimt(k+1)-epimt(k))/(st(k+1)-st(k))*(q2-st(k))
      else
         vpre = prest(k) +
     .        (prest(k+1)-prest(k))/(st(k+1)-st(k))*(aq2-st(k))
         vpim = 0.d0

c         evpre=eprest(k)
c     .        +(eprest(k+1)-eprest(k))/(st(k+1)-st(k))*(aq2-st(k))
c         evpim = 0.d0
      endif
*
      VPOLChlmnt = 1.d0/(1.d0  - vpre - im*vpim)
      return
      end
