      FUNCTION AINTSVGG(S,C,eps)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/parameters/ame,ammu,convfac,alpha,pi
      common/firstaintsvgg/ifirst
      common/dalphacmn/dallep,dalhad,edalhad,toppart
      data ifirst /0/
      double complex vpolc,cresult
      external vpolc
      common/bornforsvgg/borngg

      ebeam = dsqrt(s)/2.d0
c      print*,eps
c      stop
      egmin = eps*ebeam

c      BETA=SQRT((1.D0-AME/EBEAM)*(1.D0+AME/EBEAM))
      BETA=SQRT((1.D0-AME*AME/EBEAM/EBEAM))

      E=1.D0/BETA
      AM=AME/BETA/EBEAM
*
      V=0.5D0*LOG(4.D0/AM/AM)
      U=0.5D0*LOG(2.D0*(E+C)/AM/AM)
      W=0.5D0*LOG(2.D0*(E-C)/AM/AM)
*
      PI2=PI*PI

      ARG   = 1.D0-2.D0*(E+C)/AM/AM
      ETAC  = DDILOGhere(ARG)+1.D0/6.D0*PI2
      ARG   = 1.D0-2.D0*(E-C)/AM/AM
      ETAMC = DDILOGhere(ARG)+1.D0/6.D0*PI2
*
      AK=EGMIN/BETA/EBEAM
      C2=C*C
      AM2=AM*AM
*
      DSV1=2.D0*(1.D0-2.D0*V)*(LOG(AK)+V)+1.5D0-PI2/3.D0

      DSV2=-4.D0*V*V*(3.D0-C2)-8.D0*V*C2+4.D0*U*V*(5.D0+2.D0*C+C2)
     .     +4.D0*V*W*(5.D0-2.D0*C+C2)
     .     - U *(7.D0-8.D0*C+C2)
     .     - W *(7.D0+8.D0*C+C2)
     .     + ETAC  * (5.D0+2.D0*C+C2)
     .     + ETAMC * (5.D0-2.D0*C+C2)
      DSV2=0.5D0/(1.D0+C2)*DSV2
*
      DSV3=2.D0*U/(E+C-0.5d0*AM2)+2.D0*W/(E-C-0.5d0*AM2)
     .    +(ETAC  - PI2/3.D0)*AM2*AM2/(E+C)/(E+C)/(E+C)
     .    +(ETAMC - PI2/3.D0)*AM2*AM2/(E-C)/(E-C)/(E-C)
     .    -4.D0*U*AM2/(E+C)/(E+C)-4.D0*W*AM2/(E-C)/(E-C)
     .    +2.D0*AM2/(E+C)/(E+C)+2.D0*AM2/(E-C)/(E-C)
      DSV3=0.5D0*(E-C)*(E+C)/(1+C2)*DSV3
*
      DSV=DSV1+DSV2+DSV3
      DSV=-ALPHA/PI*DSV
*
      ALPHA2=ALPHA*ALPHA
      BORN=2.D0*PI*ALPHA2/S*BETA*BETA/
     >     (1.D0+BETA*C)/(1.D0-BETA*C)*(1.D0+C2)

      AINTSVGG=BORN*(1.D0+DSV)
cc      AINTSVGG=BORN*(DSV)


! OK      print*,'here 2',(born*s*8.d0*pi*16/borngg)

      if (ifirst.eq.0) then
         ifirst = 1

c$$$         cresult = vpolc(s)
c$$$
c$$$         s0bk = 2.d0*pi*alpha**2/s*(log(s/ame/ame)-1.d0)*convfac
c$$$
c$$$         s3gbk = 2.d0*alpha**3/s * convfac
c$$$         s3gbk = s3gbk*(-(2.d0*v-1.d0)**2*(2.d0*log(ak)+1.d0)+3.d0)
c$$$
c$$$         ssvbk = 2.d0*alpha**3/s * convfac
c$$$         ssvbk = ssvbk * (2.d0*(2.d0*v-1.d0)**2*log(ak)+4.d0/3.d0*v**3
c$$$     .   +3.d0*v*v + v*(2.d0/3.d0*pi*pi-6.d0) -1.d0/12.d0*pi*pi)
c$$$
c$$$
c$$$         print*,'s0    = ',s0bk
c$$$         print*,'s_sv  = ',s0bk+ssvbk
c$$$         print*,'s_r   = ',s3gbk
c$$$         print*,'s_tot = ',s3gbk+s0bk+ssvbk
c$$$         write(*,'(A,f15.12)')'\Delta\alpha_lep(s) = ',dallep
c$$$         write(*,'(A,f15.12,A,f15.12)')'\Delta\alpha_had(s) = ',dalhad,
c$$$     .        ' +- ',edalhad
c         stop

      endif


      return
      END
***********************************************
      FUNCTION AINTSV(S,T,eps)
! essentially from LABSPV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/parameters/ame,ammu,convfac,alpha,pi
      common/ialpharunning/iarun
      COMMON/CHANNEL/ICH
      character*2 fs
      common/finalstate/fs
      double complex vpolc,vpols,vpolt
      external vpolc

      ICH = 1

      Z = 1D0 + 2*T/S
      Z2 = Z*Z

      AOPI = ALPHA/PI
      PI2 = PI*PI
      AK0 = eps
*
      TH = DACOS(Z)
      A = DSIN(TH/2.D0)
      B = DCOS(TH/2.D0)
      A2 = A**2
      B2 = B**2
      A2OB2 = A2/B2
      A4 = A2*A2
      B4 = B2*B2
      OMB4 = 1.D0 - B4
      OPB4 = 1.D0 + B4
      ALNA = DLOG(A)
      ALNA2 = ALNA**2
      ALNB = DLOG(B)
      ALNB2 = ALNB**2
      ALAOB = DLOG(A/B)
      ALAOB2 = ALAOB**2
      ame2 = ame*ame
*
      BETAE = 2*AOPI*(DLOG(S/AME2) - 1.D0)
      BETAINT = 4*AOPI*ALAOB
      BETAT= 2.D0*AOPI*(DLOG(-T/AME2) - 1.D0)
*
      CINFRA = 1.D0 + 2.D0*(BETAE+BETAINT)*LOG(AK0)

c      DDLOG1 = DDILOG(B2)
c      DDLOG2 = DDILOG(A2)
      DDLOG1 = DDILOGhere(B2)
      DDLOG2 = DDILOGhere(A2)

! using ddidlog from TOPAZ0
!      call tspence(b2,0.d0,1.d0-b2,dddlog1,ddim)
!      call tspence(a2,0.d0,1.d0-a2,dddlog2,ddim)
!!
      FAB = 2.D0*(ALNA2 - ALNB2) + DDLOG1 - DDLOG2
*
      V1GS = -Z*(ALNA2/B4+ALNB2/A4) + ALNA/B2 - ALNB/A2
      A1GS = -Z*(ALNA2/B4-ALNB2/A4) + ALNA/B2 + ALNB/A2
*
      V1GT = 8.D0*ALNA*ALNB + 0.25D0*PI2*OMB4 + OMB4/B4*ALNA2
     #    + OMB4*ALAOB2 + A2OB2*ALNA + A2*ALAOB
      A1GT = - 0.25D0*PI2*OMB4 + OMB4/B4*ALNA2
     #    - OMB4*ALAOB2 + A2OB2*ALNA - A2*ALAOB
*
      CF1 = 1.5D0*BETAE + AOPI*(
     #      2.D0*(PI2/3.D0 - 0.5D0)
     #    + 2.D0*FAB + V1GS + 2.D0*Z/(1.D0+Z*Z)*A1GS)
      CF2 = 1.5D0*BETAE + AOPI*(
     #      2.D0*(PI2/12.D0 - 0.5D0)
     #    + 2.D0*(1.5D0*ALNA - ALNA2)
     #    + 2.D0*FAB + 0.5D0*(V1GS + A1GS + V1GT + A1GT))
      CF3 = 1.5D0*BETAE + AOPI*(
     #    - 2.D0*(PI2/6.D0 + 0.5D0)
     #    + 4.D0*(1.5D0*ALNA - ALNA2)
     #    + 2.D0*FAB + V1GT + (B4-1.D0)/(B4+1.D0)*A1GT)
*
      OMZ  = 1.D0 - Z
      OMZ2 = OMZ**2
      OPZ  = 1.D0 + Z
      OPZ2 = OPZ**2
*
      vpols = vpolc(s)
      vpolt = 1.d0
      if (fs.eq.'ee') then
         vpolt = vpolc(t)
      endif
      ! s channel
      AINTSV1 = 1.D0 + Z2
      AINTSV1 = AINTSV1*(CINFRA + CF1)
      AINTSV1 = AINTSV1*ALPHA**2/4.D0/S
      AINTSV1 = AINTSV1*2*PI * abs(vpols)**2
*
      ! s-t interference
      AINTSV2 = -2.D0/OMZ*OPZ2
      AINTSV2 = AINTSV2*(CINFRA + CF2)
      AINTSV2 = AINTSV2*ALPHA**2/4.D0/S
      AINTSV2 = AINTSV2*2*PI * vpols*vpolt
*
      ! t channel
      AINTSV3 = 2.D0/OMZ2*(OPZ2 + 4.D0)
      AINTSV3 = AINTSV3*(CINFRA + CF3)
      AINTSV3 = AINTSV3*ALPHA**2/4.D0/S
      AINTSV3 = AINTSV3*2*PI * vpolt**2
*
      if (fs.eq.'ee') then
         AINTSV = AINTSV1 + AINTSV2 + AINTSV3
      else
         AINTSV = AINTSV1
      endif
      return
      END
***********************************************
      FUNCTION AINTSV_NEW(p1,p2,p3,p4,eps)
! essentially from LABSPV
#ifdef COLLIER
      use collier
#endif
      IMPLICIT DOUBLE PRECISION (A-H,L-M,O-Z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension ptmp(0:3)
      dimension p5(0:3),p6(0:3)
      dimension charges(4),masses(4)
      common/parameters/ame,ammu,convfac,alpha,pi
      common/ialpharunning/iarun
      COMMON/CHANNEL/ICH
      character*2 fs
      common/finalstate/fs
      double complex vpolc,vpols,vpolt
      external vpolc
      external lambda
#ifdef COLLIER
#include "invariants.h"
      common/collierifirst/ifirst1
      data ifirst1 /0/
      double precision phmass
      common/softphotonmass/phmass
      double complex amplborn
      double complex boxamp,vertamp,tmp2l,ctmp1,ctmp2,ctmp3
      integer hh                !,h1,h2,h3,h4
      dimension p1h(0:3),p2h(0:3),p3h(0:3),p4h(0:3)
#endif

      amfs = ame
      if (fs.eq.'mm') amfs = ammu
      amfs2 = amfs*amfs

      flux2b = 8.d0*p1(0)*p2(0)
      mass   = p1(0)+p2(0)
      sqla   = sqrt(lambda(mass**2,amfs2,amfs2))
      phsp2b = sqla/8.d0/mass**2 /(2.d0*pi) *0.25d0 ! spins...

      ame2  = ame*ame

      p3m = dsqrt(tridot(p3,p3))
      p4m = dsqrt(tridot(p4,p4))
      betafs = p3m/p3(0)

      ICH = 1

      dp1p2 = 2.d0*dot(p1,p2)
      dp1p3 = 2.d0*dot(p1,p3)
      dp1p4 = 2.d0*dot(p1,p4)
      dp2p3 = 2.d0*dot(p2,p3)
      dp2p4 = 2.d0*dot(p2,p4)
      dp3p4 = 2.d0*dot(p3,p4)

      S     =  dot(p1+p2,p1+p2)
      T     = -dp1p3

      Z = p3(3)/p3m
      Z2 = Z*Z

      AOPI = ALPHA/PI
      PI2 = PI*PI
      AK0 = eps
*
      TH = ACOS(Z)
      A = DSIN(TH/2.D0)
      B = DCOS(TH/2.D0)
      A2 = A**2
      B2 = B**2
      A2OB2 = A2/B2
      A4 = A2*A2
      B4 = B2*B2
      OMB4 = 1.D0 - B4
      OPB4 = 1.D0 + B4
      ALNA = DLOG(A)
      ALNA2 = ALNA**2
      ALNB = DLOG(B)
      ALNB2 = ALNB**2
      ALAOB = DLOG(A/B)
      ALAOB2 = ALAOB**2
*
c      BETAE = 2*AOPI*(LOG(S/AME2) - 1.D0)
c      BETAINT = 4*AOPI*ALAOB
ccc
      BETAE=2*AOPI*(0.5d0*DLOG(dp1p2/AME2)+0.5d0*dlog(dp3p4/amfs2)-1.D0)
      BETAINT=2.d0*AOPI*dlog(dp1p3/dp2p3)
ccc
      BETAT= 2.D0*AOPI*(DLOG(-T/AME2) - 1.D0)
*
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
      coeffinfra=4.d0*aopi*eikonalintegral(p1,p2,p3,p4,charges,masses)
      coeffinfra_OLD = 2.D0*(BETAE+BETAINT)

      CINFRA = 1.D0 + coeffinfra*LOG(AK0)

c      DDLOG1 = DDILOG(B2)
c      DDLOG2 = DDILOG(A2)
      DDLOG1 = DDILOGhere(B2)
      DDLOG2 = DDILOGhere(A2)

! using ddidlog from TOPAZ0
!      call tspence(b2,0.d0,1.d0-b2,dddlog1,ddim)
!      call tspence(a2,0.d0,1.d0-a2,dddlog2,ddim)
!!
      FAB = 2.D0*(ALNA2 - ALNB2) + DDLOG1 - DDLOG2
*
      V1GS = -Z*(ALNA2/B4+ALNB2/A4) + ALNA/B2 - ALNB/A2
      A1GS = -Z*(ALNA2/B4-ALNB2/A4) + ALNA/B2 + ALNB/A2
*
      V1GT = 8.D0*ALNA*ALNB + 0.25D0*PI2*OMB4 + OMB4/B4*ALNA2
     .    + OMB4*ALAOB2 + A2OB2*ALNA + A2*ALAOB
      A1GT = - 0.25D0*PI2*OMB4 + OMB4/B4*ALNA2
     .    - OMB4*ALAOB2 + A2OB2*ALNA - A2*ALAOB
*
      CF1 = 1.5D0*BETAE + AOPI*(
     .      2.D0*(PI2/3.D0 - 0.5D0)
     .    + 2.D0*FAB + V1GS + 2.D0*Z/(1.D0+Z*Z)*A1GS)
      CF2 = 1.5D0*BETAE + AOPI*(
     .      2.D0*(PI2/12.D0 - 0.5D0)
     .    + 2.D0*(1.5D0*ALNA - ALNA2)
     .    + 2.D0*FAB + 0.5D0*(V1GS + A1GS + V1GT + A1GT))
      CF3 = 1.5D0*BETAE + AOPI*(
     .    - 2.D0*(PI2/6.D0 + 0.5D0)
     .    + 4.D0*(1.5D0*ALNA - ALNA2)
     .    + 2.D0*FAB + V1GT + (B4-1.D0)/(B4+1.D0)*A1GT)
*
      OMZ  = 1.D0 - Z
      OMZ2 = OMZ**2
      OPZ  = 1.D0 + Z
      OPZ2 = OPZ**2
*
      vpols = vpolc(s)
      vpolt = 1.d0
      if (fs.eq.'ee') then
         vpolt = vpolc(t)
      endif

      ! s channel

      AINTSV1OLD = 1.D0 + Z2

      deng2 = 1.d0/S/S
      ovall = 8.d0

      p1p2 = 0.5d0*dp1p2
      p1p3 = 0.5d0*dp1p3
      p1p4 = 0.5d0*dp1p4
      p2p3 = 0.5d0*dp2p3
      p2p4 = 0.5d0*dp2p4
      p3p4 = 0.5d0*dp3p4

      AINTSV1 = ovall*deng2 * ( p1p4*p2p3 + p1p3*p2p4 + amfs2*p1p2
     &     + ame2*p3p4 + 2.D0*ame2*amfs2 )

      AINTSV1 = AINTSV1*(CINFRA + CF1)
      AINTSV1 = AINTSV1*ALPHA**2/4.D0/S
      AINTSV1 = AINTSV1*2*PI * abs(vpols)**2
*
      ! s-t interference
      AINTSV2 = -2.D0/OMZ*OPZ2
      AINTSV2 = AINTSV2*(CINFRA + CF2)
      AINTSV2 = AINTSV2*ALPHA**2/4.D0/S
      AINTSV2 = AINTSV2*2*PI * vpols*vpolt
*
      ! t channel
      AINTSV3 = 2.D0/OMZ2*(OPZ2 + 4.D0)
      AINTSV3 = AINTSV3*(CINFRA + CF3)
      AINTSV3 = AINTSV3*ALPHA**2/4.D0/S
      AINTSV3 = AINTSV3*2*PI * vpolt**2
*
      if (fs.eq.'ee') then
         AINTSV_NEW = AINTSV1 + AINTSV2 + AINTSV3
      else
         AINTSV_NEW = AINTSV1 * betafs
      endif

#ifdef COLLIER
      if (fs.eq.'mm') then
         if (ifirst1.eq.0) then
            QRe  = 1
            QRmu = 1
            Qmu  = -1.d0
            Qmu  = 1.d0         ! AAARG

            ec2 = 4.d0*pi*alpha ! electron charge^2
            ec4 = ec2*ec2
            ec  = sqrt(ec2)

            mm  = ammu
            mm2 = mm*mm
            me  = ame
            me2 = me*me

            mm4 = mm2*mm2
            mm6 = mm2*mm4
            me4 = me2*me2
            me6 = me2*me4

            pi216   = pi*pi*16.d0 ! needed in the cts...
            pi216m1 = 1.d0/pi216

* done at the begining of svfactor in sv.F
c            call Init_cll(4,2,path(1:len_trim(path))//'cll_out')
c            call SetMode_cll(1) ! 1 is COLI, 2 DD, 3 both. 1 is faster
c            call InitCacheSystem_cll(1,4)
c            phmassmm = 1d-10
c            duv    = 0.d0
c            dimmu  = 1.d0
c           call SetDeltaUV_cll(duv)
c            call SetMuUV2_cll(dimmu)
c            call SetmuIR2_cll(phmassmm*phmassmm)

            ifirst1 = 1
         endif
**************************************************************
         p1h = -p4 ! mm^2
         p2h =  p2 ! me^2
         p3h =  p3 ! mm^2
         p4h = -p1 ! me^2

         dal = 0.d0
         dah = 0.d0

         egmax = eps * sqrt(s)*0.5d0

         call loadinvariants(p1h,p2h,p3h,p4h)
         call getQEDcts(cts)
         call InitEvent_cll(1)
         call loadQEDvertexfuns
         call loadQEDboxesfuns
         sff = soft_integral_mue(egmax,p1,p2,p3,p4)

         born = 0.d0
         tm1  = 1.d0/dot(p2h-p4h,p2h-p4h)
         ampm2 = 0.d0
         do hh = 0,7            ! only half, below the factor of 2!!
c         do hh = 0, 15           ! only half, below the factor of 2!!
c         do hh = 8, 15 !15 ! only half, below the factor of 2!!
c         do hh = 15, 15 !15 ! only half, below the factor of 2!!
            h1  = mod(hh,16)/8 + 1
            h2  = mod(hh, 8)/4 + 1
            h3  = mod(hh, 4)/2 + 1
            h4  = mod(hh, 2)   + 1

            if (Qmu.lt.0.d0) then
               call getcurrentnew(Qmu,h3,h1,p3h,p1h,mm,mm,cmu,
     .              cImu,c5mu)
            else
               call getcurrentnew(Qmu,h3,h1,p1h,p3h,mm,mm,cmu,
     .              cImu,c5mu)
            endif
            call getcurrentnew(-1.d0,h4,h2,p4h,p2h,me,me,ce,
     .           cIe,c5e)

            cmuce   = dotcmom(cmu,ce)
            amplborn = cmuce*(0.d0,1.d0)*ec2*tm1

            bbar    = amplborn*conjg(amplborn)
            born    = born + bbar

            call loadspinvariants(p1h,p2h,p3h,p4h,cmu,ce)
            call getQEDvertexes(vertamp)
            call getQEDboxes(boxamp)

            oneloop = bbar * (1.d0*cts + 1.d0*sff +
     .           2.d0*(dal+dah))
     .           + 2.d0*dreal(amplborn*conjg(vertamp))
     .           + 2.d0*dreal(amplborn*conjg(boxamp))
            ampm2 = ampm2 + 1.d0*bbar + oneloop
         enddo
         born  = 2.d0*born
         ampm2 = 2.d0*ampm2
         aintsv_new = ampm2*phsp2b/flux2b*abs(vpols)**2
      endif
#endif
      return
      END
********************************************************
      FUNCTION AINTSV_PIONS(p1,p2,p3,p4,eps)
#ifdef COLLIER
      use collier
#endif
      IMPLICIT DOUBLE PRECISION (A-H,L-M,O-Z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),ptmp(0:3)
      dimension charges(4),masses(4)
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      common/ialpharunning/iarun
      COMMON/CHANNEL/ICH
      character*2 fs
      common/finalstate/fs
      double complex vpolc,vpols,vpolt
      external vpolc
      double precision egmin
      COMMON/CUTOFF/EGMIN
      double precision k1(0:3),k2(0:3),k3(0:3),k4(0:3)
      double precision messborn
      external messborn
      double precision flippedborn
      double precision kallen
      external kallen
      double precision dot
      external dot
      double precision tridot
      external tridot
*
      double complex counterterms
      double complex vertexe,soft,vertxbornd
*
      double complex B0FFsmpimpi,B0FFmpi0mpi,B1FFmpi0mpi
      double complex C22FFe,C12FFe,C11FFe,C00FFe
      double complex C1FFe,C2FFe,C0FFe
*
      double complex C22FFpi,C12FFpi,C11FFpi,C00FFpi
      double complex C1FFpi,C2FFpi,C0FFpi
      double complex c0ltd,c1ltd,c2ltd
*
      double complex D33ltd,D23ltd,D13ltd,D22ltd,D12ltd,D11ltd,D00ltd
      double complex D3ltd,D2ltd,D1ltd,D0ltd
*
      double complex D33ltc,D23ltc,D13ltc,D22ltc,D12ltc,D11ltc,D00ltc
      double complex D3ltc,D2ltc,D1ltc,D0ltc
      double complex vertex_e,vertex_pi
      double complex vertex_seagull1,vertex_seagull2,vertex_2photons
      double complex vertex_2photons_francesco
      double complex boxd,boxc
      double complex boxd_IIrouting,boxc_IIrouting
      double complex cISR,cFSR,cIFI
      double precision vISR,vFSR,vIFI

      double precision mm2,me2
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      double precision rcounterterms_e,rcounterterms_pi,rvertexmu
      double precision rvertexe,rboxd,rboxc,rsoft
      double precision rvertex,rbox
      double precision delZ12fL,delZ12fR,delZ34
      double precision rflux,virt
      double precision born
      integer i
      double precision delZ12fLcheck
      double complex B1FFmml,B0pmml,B1pmml
      double precision be1,be2,be3,be4,e1,e2,e3,e4,s,t
      double complex ffpis
      double precision ffpis2
      integer iffpi
      character*10 what_ffpi
      real*8 s_max
      common/ffpioption/iffpi,what_ffpi
      integer jfirst
      common/jfirstsv/jfirst
      data jfirst /0/
      common/cmplxffpis/ffpis
      common/abs2ffpis/ffpis2

      double precision deltauv,mu2dim,lambda2,m2min,pi216m1,ec6
      common/duvmudim/deltauv,mudim,mu2dim,m2min,pi216m1,ec6
      common/gammamass/lambda2
*
      integer j,k,nres
      double complex im

      double complex c1
      parameter (c1 = (1.d0,0.d0))

      double precision m_n(1:9),g_n(1:9)
      double complex a_n(1:9),lv2(1:9),lv2diff
      common/par_gvmd/m_n,g_n,a_n

      double complex vertex_pi_gvmd,vertex_seagull_gvmd
      double complex rcounterterms_pi_gvmd
      double complex boxt_gvmd,boxd_gvmd,boxc_gvmd

      double complex DB0p,DB0pj(1:9)
      double complex B0s,B0sj(1:9),B0sjk(1:9,1:9)
      double complex B0p,B0pj(1:9),B0lj(1:9)
      double complex C0psp,C0pspj(1:9),C0ppj(1:9)
      double complex C0jj(1:9),C1jj(1:9),C2jj(1:9)
      double complex D0jj(1:9),D2jj(1:9)
      double complex C0mms,C0mmsj(1:9),C0mmsjk(1:9,1:9)
      double complex C0pps,C0ppsj(1:9),C0ppsjk(1:9,1:9)
      double complex D0d,D0dj(1:9),D0djk(1:9,1:9)
      double complex D0u,D0uj(1:9),D0ujk(1:9,1:9)

      integer igvmdwarning
      data igvmdwarning /0/

      
      common/indexfsrdisp/iFSRdisp
      integer iFSRdisp
#ifdef COLLIER
      integer rank
      parameter (rank = 2)
      double complex cnpi(0:rank/2,0:rank,0:rank)
      double complex cnpiuv(0:rank/2,0:rank,0:rank)
      double complex cne(0:rank/2,0:rank,0:rank)
      double complex cneuv(0:rank/2,0:rank,0:rank)
      double complex dnd(0:rank/2,0:rank,0:rank,0:rank)
      double complex dnduv(0:rank/2,0:rank,0:rank,0:rank)
      double complex dnc(0:rank/2,0:rank,0:rank,0:rank)
      double complex dncuv(0:rank/2,0:rank,0:rank,0:rank)
      double complex b(0:1/2,0:1),buv(0:1/2,0:1)
      double complex b0cll,b1cll,db0cll,db1cll
      dnd=0.d0
      dnduv=0.d0

******************************************************
!      call InitEvent_cll(1)
      call GetMuUV2_cll(mu2dim)
#endif
#ifdef LOOPTOOLS
#include "looptools.h"
      call clearcache
#endif
*
      amfs = ampi
      if (fs.eq.'mm') then
         amfs = ammu
         print*,'this sould not happen - STOP'
         stop
      endif
      if (fs.eq.'pp'.or.fs.eq.'pi') amfs = ampi
*
      vertex_e= (0.d0,0.d0)
      rcounterterms_e= (0.d0,0.d0)
      vertex_pi= (0.d0,0.d0)
      vertex_seagull1= (0.d0,0.d0)
      vertex_seagull2= (0.d0,0.d0)
      vertex_2photons= (0.d0,0.d0)
      rcounterterms_pi= (0.d0,0.d0)
      boxd= (0.d0,0.d0)
      boxc= (0.d0,0.d0)
      cISR= (0.d0,0.d0)
      cFSR= (0.d0,0.d0)
      cIFI= (0.d0,0.d0)
      rsoftISR= 0.d0
      rsoftFSR= 0.d0
      rsoftIFI= 0.d0
*
      s = dot(p1+p2,p1+p2)
      t = dot(p1-p3,p1-p3)
      u = dot(p1-p4,p1-p4)
*
      p1p2 = dot(p1,p2)
      p1p3 = dot(p1,p3)
      p1p4 = dot(p1,p4)
      p2p3 = dot(p2,p3)
      p2p4 = dot(p2,p4)
      p3p4 = dot(p3,p4)
*
      p1p1 = dot(p1,p1)
      p2p2 = dot(p2,p2)
      p3p3 = dot(p3,p3)
      p4p4 = dot(p4,p4)
*
      amfs2 = amfs*amfs
      ame2  = ame*ame
      me2   = ame2
      me4   = me2*me2
      mpi2  = ampi*ampi
      mpi4  = mpi2*mpi2
*
      betafs = sqrt(tridot(p3,p3))/p3(0)
      egmin  = eps*sqrt(s)/2.d0
c      print*,'egmin= ',egmin
*
*     counterterms (according to Denner's Fortschritte)
*
** delze (QED only) according to eq. (3.32) of ArXiv:0709.1075 (Denner Fortschritte)
*
* delZe + 1/2 delaa = 0, according to their definitions
*
*     delZ_ii^fL=delZii^fR in QED  (Eq. B12) of ArXiv:0709.1075
*         expanded for small masses
*
c      delZ12fL= -alpha/4.d0/pi*( deltauv-log(me2/mu2dim)
c     +                          +4.d0 + 2.d0*log(lambda2/me2))
c      delZ12fR= delZ12fL
*
*     the derivative with 0 or lambda2 is not the same
*     in double precision the complete expression, reported below,
*     the correct numerical value is with 0.d0  in input of DB1 and DB0 and not lambda2
*     so we use the expanded version of the ct.
*
c#ifdef LOOPTOOLS
c      B1FFmml= B0i(bb1,me2,me2,0.d0)
c      B1pmml= DB1(me2,me2,0.d0)
c      B0pmml= DB0(me2,me2,0.d0)
c#endif
c#ifdef COLLIER
c      call B_cll(b,buv,c1*me2,c1*me2,c1*0.d0,1)
c      B1FFmml = b(0,1)
c      call DB0_cll(B0pmml,c1*me2,c1*me2,c1*0d0)
c      call DB1_cll(B1pmml,c1*me2,c1*me2,c1*0d0)
c#endif
c
c      delZ12fL= alpha/4.d0/pi * (2.d0 * dble(B1FFmml)+1.d0)
c     +       + alpha/4.d0/pi*me2*(4.d0*dble(B1pmml)+8.d0*dble(B0pmml))
c      delZ12fR= delZ12fL
c
      delZ12fL= -alpha/4.d0/pi*( deltauv-log(me2/mu2dim)
     +                          +4.d0 + 2.d0*log(lambda2/me2))
      delZ12fR= delZ12fL

*     counterterm along the initial state electron line
      born= (-8*t**2-8*s*t - 8*mpi4 - 8*me4 + 16*mpi2*t + 16*me2*t
     &       +8*me2*s - 16*me2*mpi2)/s/s
      born= born*(4.d0*pi*alpha)**2
*
c      rcounterterms_e= born
c     +     *2.d0*( delZ12fL +delZ12fR)*0.5d0
      rcounterterms_e= born*(delZ12fL + delZ12fR)
*
* interference born^daga * vertex_e
*
      C0FFe=  (0.d0,0.d0)
      C1FFe=  (0.d0,0.d0)
      C2FFe=  (0.d0,0.d0)
      C00FFe= (0.d0,0.d0)
      C11FFe= (0.d0,0.d0)
      C12FFe= (0.d0,0.d0)
      C22FFe= (0.d0,0.d0)

#ifdef LOOPTOOLS
      C0FFe=  C0i(cc0,me2,s,me2,0.d0,me2,me2)
      C1FFe=  C0i(cc1,me2,s,me2,0.d0,me2,me2)
      C2FFe=  C0i(cc2,me2,s,me2,0.d0,me2,me2)
      C00FFe= C0i(cc00,me2,s,me2,0.d0,me2,me2)
      C11FFe= C0i(cc11,me2,s,me2,0.d0,me2,me2)
      C12FFe= C0i(cc12,me2,s,me2,0.d0,me2,me2)
      C22FFe= C0i(cc22,me2,s,me2,0.d0,me2,me2)
#endif
#ifdef COLLIER
      call C_cll(cne,cneuv,c1*me2,c1*s,c1*me2,c1*0.d0,c1*me2,c1*me2,
     .     rank)
      C0FFe  = cne(0,0,0)
      C1FFe  = cne(0,1,0)
      C2FFe  = cne(0,0,1)
      C00FFe = cne(1,0,0)
      C11FFe = cne(0,2,0)
      C12FFe = cne(0,1,1)
      C22FFe = cne(0,0,2)
#endif

c      print*,c0ffe,c1ffe,c2ffe,c00ffe,c11ffe,c12ffe,c22ffe
c#ifdef COLLIER
c      print*,''
c#endif
c     stop
c$$$      vertxbornd = (0.d0,0.d0)
c$$$      vertxbornd = vertxbornd
c$$$     &  + C22FFe * (  - 32*me4*t - 16*me4*s + 32*mpi2*me4 + 16*me2*t**2
c$$$     &     + 16*me2*s*t + 8*me2*s**2 + 16*me2*mpi4 + 16*me2*me4 - 32*
c$$$     &    me2*mpi2*t - 32*me2*mpi2*s )
c$$$      vertxbornd = vertxbornd + C12FFe * ( 16*s*t**2 + 16*s**2*t + 16*
c$$$     &    mpi4*s - 64*me4*t - 16*me4*s - 32*mpi2*s*t + 64*mpi2*me4 + 32
c$$$     &    *me2*t**2 + 32*me2*mpi4 + 32*me2*me4 - 64*me2*mpi2*t - 32*me2
c$$$     &    *mpi2*s )
c$$$      vertxbornd = vertxbornd + C11FFe * (  - 32*me4*t - 16*me4*s + 32*
c$$$     &    mpi2*me4 + 16*me2*t**2 + 16*me2*s*t + 8*me2*s**2 + 16*me2*
c$$$     &    mpi4 + 16*me2*me4 - 32*me2*mpi2*t - 32*me2*mpi2*s )
c$$$      vertxbornd = vertxbornd + C00FFe * (  - 32*t**2 - 32*s*t - 32*
c$$$     &    mpi4 - 32*me4 + 64*mpi2*t + 64*me2*t + 32*me2*s - 64*me2*mpi2
c$$$     &     )
c$$$      vertxbornd = vertxbornd + C2FFe * ( 16*s*t**2 + 16*s**2*t + 16*
c$$$     &    mpi4*s + 64*me4*t + 48*me4*s - 32*mpi2*s*t - 64*mpi2*me4 - 32
c$$$     &    *me2*t**2 - 64*me2*s*t - 8*me2*s**2 - 32*me2*mpi4 - 32*me2*
c$$$     &    me4 + 64*me2*mpi2*t )
c$$$      vertxbornd = vertxbornd + C1FFe * ( 16*s*t**2 + 16*s**2*t + 16*
c$$$     &    mpi4*s + 64*me4*t + 48*me4*s - 32*mpi2*s*t - 64*mpi2*me4 - 32
c$$$     &    *me2*t**2 - 64*me2*s*t - 8*me2*s**2 - 32*me2*mpi4 - 32*me2*
c$$$     &    me4 + 64*me2*mpi2*t )
c$$$      vertxbornd = vertxbornd + C0FFe * ( 16*s*t**2 + 16*s**2*t + 16*
c$$$     &    mpi4*s + 64*me4*t + 48*me4*s - 32*mpi2*s*t - 64*mpi2*me4 - 32
c$$$     &    *me2*t**2 - 64*me2*s*t - 16*me2*s**2 - 32*me2*mpi4 - 32*me2*
c$$$     &     me4 + 64*me2*mpi2*t + 32*me2*mpi2*s )
c$$$      vertxbornd = vertxbornd + 16*t**2 + 16*s*t + 16*mpi4 + 16*me4 -
c$$$     &     32*mpi2*t - 32*me2*t - 16*me2*s + 32*me2*mpi2

      Z1_= - t + 2*mpi2
      Z2_= - s + 2*Z1_
      Z2_=Z2_*s
      Z1_=Z1_*t
      Z3_=Z1_ - me4 - mpi4
      Z4_=2*Z3_
      Z2_=Z2_ + Z4_
      Z2_=Z2_*me2
      Z5_=2*me4
      Z6_=t - mpi2
      Z5_=Z5_*Z6_
      Z7_=s*me4
      Z5_=Z5_ + Z7_
      Z2_=Z2_ + 2*Z5_
      Z2_=8*Z2_
      Z5_=4*me4
      Z5_=Z5_*Z6_
      Z7_=s*t
      Z1_=Z1_ - Z7_ - mpi4
      Z8_= - me4 - Z1_
      Z8_=s*Z8_
      Z9_= - s*mpi2
      Z9_=Z9_ - Z3_
      Z9_=me2*Z9_
      Z8_=2*Z9_ - Z5_ + Z8_
      Z8_=16*Z8_
      Z6_=s + 2*Z6_
      Z6_=Z6_*me2
      Z6_=Z6_ - Z7_ + Z3_
      Z7_=32*Z6_
      Z1_=3*me4 - Z1_
      Z1_=Z1_*s
      Z1_=Z1_ + Z5_
      Z5_=s + 8*t
      Z5_=Z5_*s
      Z3_= - Z5_ + 4*Z3_
      Z3_=Z3_*me2
      Z3_=Z3_ + 2*Z1_
      Z3_=8*Z3_
      Z5_=mpi2 - 2*t
      Z5_=2*Z5_ - s
      Z5_=s*Z5_
      Z4_=Z4_ + Z5_
      Z4_=me2*Z4_
      Z1_=Z4_ + Z1_
      Z1_=16*Z1_
      Z4_= - 16*Z6_
      vertxbornd= - C22FFe*Z2_ + C12FFe*Z8_ - C11FFe*Z2_ + C00FFe*Z7_
     &  + C2FFe*Z3_ + C1FFe*Z3_ + C0FFe*Z1_ + Z4_

*** STATS: original	0P 210M 80A : 290
*** STATS: optimized 0P 36M 26A : 62
*     multiplying by e^6 / 16/pi^2 / s^2 * 2   !2 is due to (A Bd + Ad B)
*     the overall phase is 1, corresponing to (-i) x (+i), the latter from Bd
      vertxbornd = vertxbornd * 2.d0
     +             * ec6 * pi216m1/s/s

      vertex_e= vertxbornd
*
      if(fs.eq.'pp') then
*  counterterm on the pion external legs
*     delZ34 is already summed over L and R
#ifdef LOOPTOOLS
         delZ34= alpha/pi*(B0(mpi2,0.d0,mpi2) + B1(mpi2,0.d0,mpi2)
     +         +mpi2*(DB0(mpi2,0.d0,mpi2) + DB1(mpi2,0.d0,mpi2)) )
#endif
#ifdef COLLIER
         call B_cll(b,buv,c1*mpi2,c1*0.d0,c1*mpi2,1)
         b0cll = b(0,0)
         b1cll = b(0,1)
         call DB0_cll(db0cll,c1*mpi2,c1*0.d0,c1*mpi2)
         call DB1_cll(db1cll,c1*mpi2,c1*0.d0,c1*mpi2)
         delZ34= alpha/pi*(b0cll + b1cll + mpi2*(db0cll + db1cll ))
#endif
ccc   delZ34 Francesco
ccc      delZ34= alpha/2.d0/pi*(B0(mpi2,0.d0,mpi2)
ccc     +                      +2.d0*mpi2*DB0(mpi2,0.d0,mpi2))
ccc   delZ34 expanded
c$$$  delZ34= -2.d0*(log(lambda2/mpi2)-log(mu2dim/mpi2)-deltauvf)
c$$$      delZ34= delZ34*alpha/4.d0/pi
         rcounterterms_pi= 2.d0 * delZ34 * born
*
* interference born^daga * vertex_pi
*
         C0FFpi=  (0.d0,0.d0)
         C1FFpi=  (0.d0,0.d0)
         C2FFpi=  (0.d0,0.d0)
         C00FFpi= (0.d0,0.d0)
         C11FFpi= (0.d0,0.d0)
         C12FFpi= (0.d0,0.d0)
         C22FFpi= (0.d0,0.d0)
         B0FFsmpimpi= (0.d0,0.d0)
#ifdef LOOPTOOLS
         C0FFpi =  C0i(cc0,mpi2,s,mpi2,0.d0,mpi2,mpi2)
         C1FFpi =  C0i(cc1,mpi2,s,mpi2,0.d0,mpi2,mpi2)
         C2FFpi =  C0i(cc2,mpi2,s,mpi2,0.d0,mpi2,mpi2)
         C00FFpi= C0i(cc00,mpi2,s,mpi2,0.d0,mpi2,mpi2)
         C11FFpi= C0i(cc11,mpi2,s,mpi2,0.d0,mpi2,mpi2)
         C12FFpi= C0i(cc12,mpi2,s,mpi2,0.d0,mpi2,mpi2)
         C22FFpi= C0i(cc22,mpi2,s,mpi2,0.d0,mpi2,mpi2)
         B0FFsmpimpi= B0i(bb0,s,mpi2,mpi2)
#endif
#ifdef COLLIER
         call C_cll(cnpi,cnpiuv,c1*mpi2,c1*s,c1*mpi2,c1*0.d0,c1*mpi2,
     .        c1*mpi2,rank)
         C0FFpi  = cnpi(0,0,0)
         C1FFpi  = cnpi(0,1,0)
         C2FFpi  = cnpi(0,0,1)
         C00FFpi = cnpi(1,0,0)
         C11FFpi = cnpi(0,2,0)
         C12FFpi = cnpi(0,1,1)
         C22FFpi = cnpi(0,0,2)
         call b_cll(b,buv,c1*s,c1*mpi2,c1*mpi2,1)
         B0FFsmpimpi= b(0,0)
#endif
*
c$$$         vertxbornd = (0.d0,0.d0)
c$$$     &  + C22FFpi * ( 8*s*t**2 + 8*s**2*t + 80*mpi4*t + 8*mpi4*s + 8*
c$$$     &    me4*s - 40*mpi2*t**2 - 56*mpi2*s*t - 40*mpi2*mpi4 - 40*mpi2*
c$$$     &    me4 - 16*me2*s*t - 8*me2*s**2 - 80*me2*mpi4 + 80*me2*mpi2*t
c$$$     &     + 56*me2*mpi2*s )
c$$$         vertxbornd = vertxbornd + C12FFpi * ( 24*s*t**2 + 24*s**2*t+160
c$$$     &    *mpi4*t + 24*mpi4*s + 24*me4*s - 80*mpi2*t**2 - 128*mpi2*s*t
c$$$     &     - 80*mpi2*mpi4 - 80*mpi2*me4 - 48*me2*s*t - 24*me2*s**2 -
c$$$     &    160*me2*mpi4 + 160*me2*mpi2*t + 128*me2*mpi2*s )
c$$$         vertxbornd = vertxbornd + C11FFpi * ( 8*s*t**2 + 8*s**2*t + 80*
c$$$     &    mpi4*t + 8*mpi4*s + 8*me4*s - 40*mpi2*t**2 - 56*mpi2*s*t - 40
c$$$     &    *mpi2*mpi4 - 40*mpi2*me4 - 16*me2*s*t - 8*me2*s**2 - 80*me2*
c$$$     &    mpi4 + 80*me2*mpi2*t + 56*me2*mpi2*s )
c$$$         vertxbornd = vertxbornd + C00FFpi * (  - 64*t**2 - 64*s*t - 64*
c$$$     &    mpi4 - 64*me4 + 128*mpi2*t + 128*me2*t + 64*me2*s - 128*me2*
c$$$     &    mpi2 )
c$$$         vertxbornd = vertxbornd + C2FFpi * ( 24*s*t**2+ 24*s**2*t+ 128*
c$$$     &    mpi4*t + 24*mpi4*s + 24*me4*s - 64*mpi2*t**2 - 112*mpi2*s*t
c$$$     &     - 64*mpi2*mpi4 - 64*mpi2*me4 - 48*me2*s*t - 24*me2*s**2 -
c$$$     &    128*me2*mpi4 + 128*me2*mpi2*t + 112*me2*mpi2*s )
c$$$         vertxbornd = vertxbornd + C1FFpi * ( 24*s*t**2+ 24*s**2*t+ 128*
c$$$     &    mpi4*t + 24*mpi4*s + 24*me4*s - 64*mpi2*t**2 - 112*mpi2*s*t
c$$$     &     - 64*mpi2*mpi4 - 64*mpi2*me4 - 48*me2*s*t - 24*me2*s**2 -
c$$$     &    128*me2*mpi4 + 128*me2*mpi2*t + 112*me2*mpi2*s )
c$$$         vertxbornd = vertxbornd + C0FFpi * ( 16*s*t**2 + 16*s**2*t+ 64*
c$$$     &    mpi4*t + 16*mpi4*s + 16*me4*s - 32*mpi2*t**2 - 64*mpi2*s*t -
c$$$     &    32*mpi2*mpi4 - 32*mpi2*me4 - 32*me2*s*t - 16*me2*s**2 - 64*
c$$$     &    me2*mpi4 + 64*me2*mpi2*t + 64*me2*mpi2*s )
c$$$         vertxbornd = vertxbornd + B0FFsmpimpi * ( 8*t**2 + 8*s*t+8*mpi4
c$$$     &     + 8*me4 - 16*mpi2*t - 16*me2*t - 8*me2*s + 16*me2*mpi2 )
c$$$         vertxbornd = vertxbornd + 4*t**2 + 4*s*t + 4*mpi4 + 4*me4 - 8*
c$$$     &        mpi2*t - 8*me2*t - 4*me2*s + 8*me2*mpi2

      Z1_= - s + 7*mpi2
      Z2_=t - me2
      Z1_=Z1_*Z2_
      Z3_= - t + 2*me2
      Z3_=Z3_*t
      Z3_=Z3_ - mpi4 - me4
      Z1_=Z3_ + Z1_
      Z1_=Z1_*s
      Z4_=2*mpi4
      Z4_=Z4_*Z2_
      Z5_=Z3_*mpi2
      Z4_=Z5_ + Z4_
      Z1_= - Z1_ + 5*Z4_
      Z1_=8*Z1_
      Z5_=3*s
      Z6_=Z5_ - 16*mpi2
      Z6_=Z2_*Z6_
      Z7_=3*Z3_
      Z6_= - Z7_ + Z6_
      Z6_=s*Z6_
      Z6_=10*Z4_ + Z6_
      Z6_=8*Z6_
      Z8_= - s + 2*mpi2
      Z8_=Z8_*Z2_
      Z8_=Z3_ + Z8_
      Z9_=64*Z8_
      Z5_= - Z5_ + 14*mpi2
      Z5_=Z5_*Z2_
      Z5_=Z7_ + Z5_
      Z5_=Z5_*s
      Z5_= - Z5_ + 8*Z4_
      Z5_=8*Z5_
      Z7_=s - 4*mpi2
      Z2_=Z2_*Z7_
      Z2_=Z2_ - Z3_
      Z2_=s*Z2_
      Z2_=2*Z4_ + Z2_
      Z2_=16*Z2_
      Z3_= - 8*Z8_
      Z4_= - 4*Z8_
      vertxbornd=C22FFpi*Z1_ + C12FFpi*Z6_ + C11FFpi*Z1_ + C00FFpi*Z9_
     &  + C2FFpi*Z5_ + C1FFpi*Z5_ + C0FFpi*Z2_ + B0FFsmpimpi*Z3_ + Z4_

***   STATS: original	0P 258M 99A : 357
***   STATS: optimized 0P 32M 19A : 51

*     multiplying by e^6 / 16/pi^2 / s^2 * 2   !2 is due to (A Bd + Ad B)
*     the overall phase is 1, corresponing to (-i) x (+i), the latter from Bd
         vertxbornd = vertxbornd * 2.d0
     +     * ec6 * pi216m1/s/s


         vertex_pi= vertxbornd
*     vertex correction with the virtual photon coming from the 4-linear vertex
#ifdef LOOPTOOLS
         B0FFmpi0mpi= B0i(bb0,mpi2,0.d0,mpi2)
         B1FFmpi0mpi= B0i(bb1,mpi2,0.d0,mpi2)
#endif
#ifdef COLLIER
         call b_cll(b,buv,c1*mpi2,c1*0.d0,c1*mpi2,1)
         B0FFmpi0mpi= b(0,0)
         B1FFmpi0mpi= b(0,1)
#endif
c         vertxbornd = (0.d0,0.d0)
c     &        + B1FFmpi0mpi * (-4*t**2-4*s*t-4*mpi4 - 4*me4 + 8*mpi2*t
c     &        + 8*me2*t + 4*me2*s - 8*me2*mpi2 )
c         vertxbornd = vertxbornd + B0FFmpi0mpi * (- 8*t**2 - 8*s*t - 8*
c     &        mpi4 - 8*me4 + 16*mpi2*t + 16*me2*t + 8*me2*s - 16*me2*mpi2 )
         Z1_= - s + 2*mpi2
         Z2_=2*t - Z1_
         Z2_=Z2_*me2
         Z1_=t - Z1_
         Z1_=Z1_*t
         Z1_= - Z2_ + Z1_ + mpi4 + me4
         Z2_= - 4*Z1_
         Z1_= - 8*Z1_
         vertxbornd=B1FFmpi0mpi*Z2_ + B0FFmpi0mpi*Z1_
*** STATS: original	0P 28M 14A : 42
*** STATS: optimized 0P 6M 6A : 12

         vertxbornd = vertxbornd * 2.d0 * (-2.d0) !one of the 2's is for interference with Born  (the phase is (+i)*(+i) the latter from Bornd)
     +        * ec6 * pi216m1/s/s

         vertex_seagull1= vertxbornd

*     vertex correction with the virtual photon coming from the 4-linear vertex
#ifdef LOOPTOOLS
         B0FFmpi0mpi= B0i(bb0,mpi2,0.d0,mpi2)
         B1FFmpi0mpi= B0i(bb1,mpi2,0.d0,mpi2)
#endif
#ifdef COLLIER
         call b_cll(b,buv,c1*mpi2,c1*0.d0,c1*mpi2,1)
         B0FFmpi0mpi= b(0,0)
         B1FFmpi0mpi= b(0,1)
#endif
c         vertxbornd = (0.d0,0.d0)
c     &        + B1FFmpi0mpi * (  - 4*t**2 - 4*s*t-4*mpi4-4*me4+8*mpi2*t
c     &        + 8*me2*t + 4*me2*s - 8*me2*mpi2 )
c         vertxbornd = vertxbornd + B0FFmpi0mpi * (- 8*t**2 - 8*s*t - 8*
c     &        mpi4 - 8*me4 + 16*mpi2*t + 16*me2*t+ 8*me2*s-16*me2*mpi2 )

         Z1_= - s + 2*mpi2
         Z2_=2*t - Z1_
         Z2_=Z2_*me2
         Z1_=t - Z1_
         Z1_=Z1_*t
         Z1_= - Z2_ + Z1_ + mpi4 + me4
         Z2_= - 4*Z1_
         Z1_= - 8*Z1_
         vertxbornd=B1FFmpi0mpi*Z2_ + B0FFmpi0mpi*Z1_
*** STATS: original	0P 28M 14A : 42
*** STATS: optimized 0P 6M 6A : 12

         vertxbornd = +vertxbornd * 2.d0 * (-2.d0) !one of the 2's is for interference with Born  (the phase is (+i)*(+i) the latter from Bornd)
     +        * ec6 * pi216m1/s/s

         vertex_seagull2= vertxbornd

*     vertex correction with two photons joining the 4-linear vertex
         vertxbornd= (0.d0,0.d0)
         C0FFe=  (0.d0,0.d0)
         C1FFe=  (0.d0,0.d0)
#ifdef LOOPTOOLS
         C0FFe=  C0i(cc0,me2,me2,s,0.d0,me2,0.d0)
         C1FFe=  C0i(cc1,me2,me2,s,0.d0,me2,0.d0)
#endif
#ifdef COLLIER
         call C_cll(cne,cneuv,c1*me2,c1*me2,c1*s,c1*0.d0,c1*me2,c1*0.d0,
     .        rank)
         C0FFe = cne(0,0,0)
         C1FFe = cne(0,1,0)
#endif
c         vertxbornd = (0.d0,0.d0)
c     &        + C1FFe * (  - 16*me4 + 16*me2*t + 8*me2*s - 16*me2*mpi2 )
c         vertxbornd = vertxbornd + C0FFe * ( 16*me4 - 16*me2*t - 8*me2*s
c     &        + 16*me2*mpi2 )

         Z1_=mpi2 - t
         Z1_=Z1_*me2
         Z1_=Z1_ + me4
         Z2_=s*me2
         Z1_= - Z2_ + 2*Z1_
         Z1_=8*Z1_
         vertxbornd= - C1FFe*Z1_ + C0FFe*Z1_
*** STATS: original	0P 14M 6A : 20
*** STATS: optimized 0P 4M 3A : 7

         vertxbornd = +vertxbornd * 2.d0 * (2.d0) !one of the 2's is for interference with Born (the phase is (-i)*(+i) the latter from Bornd)
     +        * ec6 * pi216m1/s

         vertex_2photons= vertxbornd
*
*     direct box
*
#ifdef LOOPTOOLS
         D0ltd=  D0i(dd0,mpi2,mpi2,me2,me2,s,t,0.d0,mpi2,0.d0,me2)
         D1ltd=  D0i(dd1,mpi2,mpi2,me2,me2,s,t,0.d0,mpi2,0.d0,me2)
         D2ltd=  D0i(dd2,mpi2,mpi2,me2,me2,s,t,0.d0,mpi2,0.d0,me2)
         D3ltd=  D0i(dd3,mpi2,mpi2,me2,me2,s,t,0.d0,mpi2,0.d0,me2)
         D00ltd=D0i(dd00,mpi2,mpi2,me2,me2,s,t,0.d0,mpi2,0.d0,me2)
         D11ltd=D0i(dd11,mpi2,mpi2,me2,me2,s,t,0.d0,mpi2,0.d0,me2)
         D12ltd=D0i(dd12,mpi2,mpi2,me2,me2,s,t,0.d0,mpi2,0.d0,me2)
         D13ltd=D0i(dd13,mpi2,mpi2,me2,me2,s,t,0.d0,mpi2,0.d0,me2)
         D22ltd=D0i(dd22,mpi2,mpi2,me2,me2,s,t,0.d0,mpi2,0.d0,me2)
         D23ltd=D0i(dd23,mpi2,mpi2,me2,me2,s,t,0.d0,mpi2,0.d0,me2)
         D33ltd=D0i(dd33,mpi2,mpi2,me2,me2,s,t,0.d0,mpi2,0.d0,me2)
         C0ltd= (0.d0,0.d0)
         C0ltd=  C0i(cc0,mpi2,me2,t,mpi2,0.d0,me2)
         C1ltd=  C0i(cc1,mpi2,me2,t,mpi2,0.d0,me2)
         C2ltd=  C0i(cc2,mpi2,me2,t,mpi2,0.d0,me2)
#endif
#ifdef COLLIER
         call D_cll(dnd,dnduv,c1*mpi2,c1*mpi2,c1*me2,c1*me2,c1*s,c1*t,
     .        c1*0.d0,c1*mpi2,c1*0.d0,c1*me2,rank)
         D0ltd  = dnd(0,0,0,0)
         D1ltd  = dnd(0,1,0,0)
         D2ltd  = dnd(0,0,1,0)
         D3ltd  = dnd(0,0,0,1)
         D00ltd = dnd(1,0,0,0)
         D11ltd = dnd(0,2,0,0)
         D12ltd = dnd(0,1,1,0)
         D13ltd = dnd(0,1,0,1)
         D22ltd = dnd(0,0,2,0)
         D23ltd = dnd(0,0,1,1)
         D33ltd = dnd(0,0,0,2)
         C0ltd= (0.d0,0.d0)
         call C_cll(cne,cneuv,c1*mpi2,c1*me2,c1*t,c1*mpi2,c1*0.d0,
     .        c1*me2,rank)
         C0ltd = cne(0,0,0)
         C1ltd = cne(0,1,0)
         C2ltd = cne(0,0,1)
#endif
*
c$$$         boxd = (0.d0,0.d0)
c$$$     &  + D33ltd * (  - 48*me4*t - 24*me4*s + 48*mpi2*me4 + 16*me2*t**2
c$$$     &     + 16*me2*s*t + 16*me2*mpi4 + 32*me2*me4 - 32*me2*mpi2*t )
c$$$         boxd = boxd + D23ltd * ( 16*s*t**2 + 16*s**2*t+ 16*mpi4*s + 24*
c$$$     &    me4*s - 32*mpi2*s*t - 40*me2*s*t - 20*me2*s**2 + 40*me2*mpi2*
c$$$     &    s )
c$$$         boxd = boxd + D13ltd * (  - 16*t**3 - 16*s*t**2 - 48*mpi4*t-80*
c$$$     &    me4*t - 28*me4*s + 48*mpi2*t**2 + 16*mpi2*s*t + 16*mpi2*mpi4
c$$$     &     + 80*mpi2*me4 + 64*me2*t**2 + 44*me2*s*t + 64*me2*mpi4 + 32*
c$$$     &    me2*me4 - 128*me2*mpi2*t - 20*me2*mpi2*s )
c$$$         boxd = boxd + D22ltd * ( 16*s*t**2 + 16*s**2*t + 16*mpi4*s+16*
c$$$     &    me4*s - 32*mpi2*s*t - 32*me2*s*t - 16*me2*s**2 + 32*me2*mpi2*
c$$$     &    s )
c$$$         boxd = boxd + D12ltd * ( 20*s*t**2 + 20*s**2*t + 20*mpi4*s+20*
c$$$     &    me4*s - 40*mpi2*s*t - 40*me2*s*t - 20*me2*s**2 + 40*me2*mpi2*
c$$$     &    s )
c$$$         boxd = boxd + D11ltd * (  - 4*t**3 - 4*s*t**2 - 44*mpi4*t - 12*
c$$$     &    me4*t - 4*me4*s + 28*mpi2*t**2 + 20*mpi2*s*t + 20*mpi2*mpi4
c$$$     &     + 28*mpi2*me4 + 12*me2*t**2 + 8*me2*s*t + 44*me2*mpi4 + 4*
c$$$     &    me2*me4 - 56*me2*mpi2*t - 20*me2*mpi2*s )
c$$$         boxd = boxd + D00ltd * ( 64*t**2 + 64*s*t + 64*mpi4 + 80*me4 -
c$$$     &    128*mpi2*t - 144*me2*t - 72*me2*s + 144*me2*mpi2 )
c$$$         boxd = boxd + D3ltd * (  - 16*t**3 - 16*s*t**2 - 48*mpi4*t-112*
c$$$     &    me4*t - 32*me4*s + 48*mpi2*t**2 + 16*mpi2*s*t + 16*mpi2*mpi4
c$$$     &     + 80*mpi2*me4 + 80*me2*t**2 + 48*me2*s*t + 48*me2*mpi4 + 48*
c$$$     &    me2*me4 - 128*me2*mpi2*t - 16*me2*mpi2*s )
c$$$         boxd = boxd + D2ltd * ( 16*s*t**2 + 16*s**2*t + 16*mpi4*s + 16*
c$$$     &    me4*s - 32*mpi2*s*t - 32*me2*s*t - 16*me2*s**2 + 32*me2*mpi2*
c$$$     &    s )
c$$$         boxd = boxd + D1ltd * (- 16*t**3 - 16*s*t**2 - 80*mpi4*t - 48*
c$$$     &    me4*t - 16*me4*s + 64*mpi2*t**2 + 32*mpi2*s*t + 32*mpi2*mpi4
c$$$     &     + 64*mpi2*me4 + 48*me2*t**2 + 32*me2*s*t + 80*me2*mpi4 + 16*
c$$$     &        me2*me4 - 128*me2*mpi2*t - 32*me2*mpi2*s )
c$$$         boxd = boxd + D0ltd * (  - 16*t**3 - 16*s*t**2-48*mpi4*t - 48*
c$$$     &    me4*t - 16*me4*s + 48*mpi2*t**2 + 16*mpi2*s*t + 16*mpi2*mpi4
c$$$     &     + 48*mpi2*me4 + 48*me2*t**2 + 32*me2*s*t + 48*me2*mpi4 + 16*
c$$$     &    me2*me4 - 96*me2*mpi2*t - 16*me2*mpi2*s )
c$$$         boxd = boxd + C2ltd * (  - 4*t**2 - 4*s*t - 4*mpi4 + 4*me4 + 8*
c$$$     &    mpi2*t )
c$$$         boxd = boxd + C1ltd * (  - 4*t**2 - 4*s*t - 4*mpi4 - 4*me4 + 8*
c$$$     &    mpi2*t + 8*me2*t + 4*me2*s - 8*me2*mpi2 )
c$$$         boxd = boxd + C0ltd * (  - 4*t**2 - 4*s*t - 4*mpi4 - 4*me4 + 8*
c$$$     &        mpi2*t + 8*me2*t + 4*me2*s - 8*me2*mpi2 )

         Z1_=mpi2*t
         Z2_=t**2
         Z3_= - Z2_ + 2*Z1_
         Z4_=mpi4 - Z3_
         Z5_=s*t
         Z5_=Z4_ + Z5_
         Z6_=me2*Z5_
         Z7_=mpi2 - t
         Z8_=2*Z7_
         Z9_=Z8_ - s
         Z10_=3*Z9_ + 4*me2
         Z10_=me4*Z10_
         Z6_=2*Z6_ + Z10_
         Z6_=8*Z6_
         Z10_=Z5_*s
         Z9_=Z9_*me2
         Z11_=s*Z9_
         Z12_=me4*s
         Z13_=6*Z12_ + 4*Z10_ + 5*Z11_
         Z13_=4*Z13_
         Z14_=Z1_ - Z2_
         Z14_=Z14_*s
         Z15_=3*t
         Z16_=Z15_ - mpi2
         Z17_=Z16_*mpi4
         Z18_=Z2_*mpi2
         Z19_=t**3
         Z14_= - Z19_ + Z14_ - Z17_ + 3*Z18_
         Z17_=5*mpi2
         Z20_= - Z17_ + 11*t
         Z21_=s*Z20_
         Z21_=16*Z4_ + Z21_
         Z21_=me2*Z21_
         Z22_=8*me2 + 20*Z7_ - 7*s
         Z22_=me4*Z22_
         Z21_=Z22_ + 4*Z14_ + Z21_
         Z21_=4*Z21_
         Z10_=Z10_ + Z11_ + Z12_
         Z11_=16*Z10_
         Z10_=20*Z10_
         Z12_=2*t
         Z22_=Z12_ - Z17_
         Z22_=s*Z22_
         Z23_=3*Z2_
         Z22_=Z22_ + 11*mpi4 + Z23_ - 14*Z1_
         Z22_=me2*Z22_
         Z24_= - Z2_ + 5*Z1_
         Z24_=s*Z24_
         Z25_=me2 - s
         Z15_=Z25_ - Z15_
         Z26_=7*mpi2 + Z15_
         Z26_=me4*Z26_
         Z20_= - mpi4*Z20_
         Z20_=Z26_ + Z22_ + Z24_ + Z20_ - Z19_ + 7*Z18_
         Z20_=4*Z20_
         Z22_=10*me4 + 8*Z5_ + 9*Z9_
         Z22_=8*Z22_
         Z1_=8*Z1_
         Z16_=s*Z16_
         Z2_=Z16_ + 3*mpi4 + 5*Z2_ - Z1_
         Z2_=me2*Z2_
         Z16_=3*me2 - 2*s - 7*t + Z17_
         Z16_=me4*Z16_
         Z2_=Z16_ + Z2_ + Z14_
         Z2_=16*Z2_
         Z8_= - s*Z8_
         Z1_=Z8_ + 5*mpi4 + Z23_ - Z1_
         Z1_=me2*Z1_
         Z3_=s*Z3_
         Z8_=4*mpi2 + Z15_
         Z8_=me4*Z8_
         Z15_= - 5*t + 2*mpi2
         Z15_=mpi4*Z15_
         Z1_=Z8_ + Z1_ + Z3_ + Z15_ - Z19_ + 4*Z18_
         Z1_=16*Z1_
         Z3_=Z12_ - mpi2
         Z3_=s*Z3_
         Z3_=3*Z4_ + Z3_
         Z3_=me2*Z3_
         Z4_=3*Z7_ + Z25_
         Z4_=me4*Z4_
         Z3_=Z4_ + Z3_ + Z14_
         Z3_=16*Z3_
         Z4_=me4 - Z5_
         Z4_=4*Z4_
         Z5_=me4 + Z9_ + Z5_
         Z5_=4*Z5_
         boxd=D33ltd*Z6_ + D23ltd*Z13_ + D13ltd*Z21_ + D22ltd*Z11_ +
     &        D12ltd*Z10_ + D11ltd*Z20_ + D00ltd*Z22_ + D3ltd*Z2_
     &        + D2ltd*Z11_
     &        + D1ltd*Z1_ + D0ltd*Z3_ + C2ltd*Z4_ - C1ltd*Z5_
     &        - C0ltd*Z5_
***   STATS: original	5P 337M 130A : 477
***   STATS: optimized 1P 83M 63A : 148

         boxd = boxd * 2.d0 * (-1.d0) !2 is for interference with Born^dagger
     +        * ec6* pi216m1/s
*                            (the phase is (+i)*(+i) the latter from Bornd)
*
*     crossed box
*
#ifdef LOOPTOOLS
         D0ltd=  D0i(dd0,mpi2,mpi2,me2,me2,s,u,0.d0,mpi2,0.d0,me2)
         D1ltd=  D0i(dd1,mpi2,mpi2,me2,me2,s,u,0.d0,mpi2,0.d0,me2)
         D2ltd=  D0i(dd2,mpi2,mpi2,me2,me2,s,u,0.d0,mpi2,0.d0,me2)
         D3ltd=  D0i(dd3,mpi2,mpi2,me2,me2,s,u,0.d0,mpi2,0.d0,me2)
         D00ltd=D0i(dd00,mpi2,mpi2,me2,me2,s,u,0.d0,mpi2,0.d0,me2)
         D11ltd=D0i(dd11,mpi2,mpi2,me2,me2,s,u,0.d0,mpi2,0.d0,me2)
         D12ltd=D0i(dd12,mpi2,mpi2,me2,me2,s,u,0.d0,mpi2,0.d0,me2)
         D13ltd=D0i(dd13,mpi2,mpi2,me2,me2,s,u,0.d0,mpi2,0.d0,me2)
         D22ltd=D0i(dd22,mpi2,mpi2,me2,me2,s,u,0.d0,mpi2,0.d0,me2)
         D23ltd=D0i(dd23,mpi2,mpi2,me2,me2,s,u,0.d0,mpi2,0.d0,me2)
         D33ltd=D0i(dd33,mpi2,mpi2,me2,me2,s,u,0.d0,mpi2,0.d0,me2)
         C0ltd=  C0i(cc0,mpi2,me2,u,mpi2,0.d0,me2)
         C1ltd=  C0i(cc1,mpi2,me2,u,mpi2,0.d0,me2)
         C2ltd=  C0i(cc2,mpi2,me2,u,mpi2,0.d0,me2)
#endif
#ifdef COLLIER
         call D_cll(dnd,dnduv,c1*mpi2,c1*mpi2,c1*me2,c1*me2,c1*s,c1*u,
     .        c1*0.d0,c1*mpi2,c1*0.d0,c1*me2,rank)
         D0ltd  = dnd(0,0,0,0)
         D1ltd  = dnd(0,1,0,0)
         D2ltd  = dnd(0,0,1,0)
         D3ltd  = dnd(0,0,0,1)
         D00ltd = dnd(1,0,0,0)
         D11ltd = dnd(0,2,0,0)
         D12ltd = dnd(0,1,1,0)
         D13ltd = dnd(0,1,0,1)
         D22ltd = dnd(0,0,2,0)
         D23ltd = dnd(0,0,1,1)
         D33ltd = dnd(0,0,0,2)
         call C_cll(cne,cneuv,c1*mpi2,c1*me2,c1*u,c1*mpi2,c1*0.d0,
     .        c1*me2,rank)
         C0ltd=  cne(0,0,0)
         C1ltd=  cne(0,1,0)
         C2ltd=  cne(0,0,1)
#endif
*
c$$$         boxc = (0.d0,0.d0)
c$$$     &  + D33ltd * ( 16*me4*t + 8*me4*s - 16*mpi2*me4 - 16*me2*t**2 -
c$$$     &    16*me2*s*t - 16*me2*mpi4 + 32*me2*mpi2*t )
c$$$         boxc = boxc + D23ltd * (  - 16*s*t**2 - 16*s**2*t-16*mpi4*s-8
c$$$     &    *me4*s + 32*mpi2*s*t + 24*me2*s*t + 12*me2*s**2 - 24*me2*mpi2
c$$$     &    *s )
c$$$         boxc = boxc + D13ltd * (  - 16*t**3 - 32*s*t**2 - 16*s**2*t-48*
c$$$     &    mpi4*t - 16*mpi4*s - 16*me4*t - 12*me4*s + 48*mpi2*t**2 + 48*
c$$$     &    mpi2*s*t + 16*mpi2*mpi4 + 16*mpi2*me4 + 32*me2*t**2 + 44*me2*
c$$$     &    s*t + 12*me2*s**2 + 32*me2*mpi4 - 64*me2*mpi2*t - 36*me2*mpi2
c$$$     &    *s )
c$$$         boxc = boxc + D22ltd * (  - 16*s*t**2 - 16*s**2*t - 16*mpi4*s -
c$$$     &    16*me4*s + 32*mpi2*s*t + 32*me2*s*t + 16*me2*s**2 - 32*me2*
c$$$     &    mpi2*s )
c$$$         boxc = boxc + D12ltd * (  - 20*s*t**2 - 20*s**2*t - 20*mpi4*s -
c$$$     &    20*me4*s + 40*mpi2*s*t + 40*me2*s*t + 20*me2*s**2 - 40*me2*
c$$$     &    mpi2*s )
c$$$         boxc = boxc + D11ltd * (  - 4*t**3 - 8*s*t**2 - 4*s**2*t + 20*
c$$$     &    mpi4*t - 4*mpi4*s - 12*me4*t - 8*me4*s - 4*mpi2*t**2 - 4*mpi2
c$$$     &    *s*t - 12*mpi2*mpi4 - 4*mpi2*me4 + 12*me2*t**2 + 16*me2*s*t
c$$$     &     + 4*me2*s**2 - 20*me2*mpi4 + 4*me2*me4 + 8*me2*mpi2*t + 4*
c$$$     &    me2*mpi2*s )
c$$$         boxc = boxc + D00ltd * (  - 64*t**2 - 64*s*t - 64*mpi4 - 48*me4
c$$$     &     + 128*mpi2*t + 112*me2*t + 56*me2*s - 112*me2*mpi2 )
c$$$         boxc = boxc + D3ltd * (  - 16*t**3 - 32*s*t**2-16*s**2*t - 48*
c$$$     &    mpi4*t - 16*mpi4*s + 16*me4*t + 16*me4*s + 48*mpi2*t**2 + 48*
c$$$     &    mpi2*s*t + 16*mpi2*mpi4 - 48*mpi2*me4 + 16*me2*t**2 + 16*me2*
c$$$     &    s*t - 16*me2*mpi4 - 16*me2*me4 + 16*me2*mpi2*s )
c$$$         boxc = boxc + D2ltd * (  - 16*s*t**2 - 16*s**2*t-16*mpi4*s-16
c$$$     &    *me4*s + 32*mpi2*s*t + 32*me2*s*t + 16*me2*s**2 - 32*me2*mpi2
c$$$     &    *s )
c$$$         boxc = boxc + D1ltd * (- 16*t**3 - 32*s*t**2 - 16*s**2*t - 16*
c$$$     &    mpi4*t - 16*mpi4*s - 48*me4*t - 32*me4*s + 32*mpi2*t**2 + 32*
c$$$     &    mpi2*s*t + 32*mpi2*me4 + 48*me2*t**2 + 64*me2*s*t + 16*me2*
c$$$     &    s**2 + 16*me2*mpi4 + 16*me2*me4 - 64*me2*mpi2*t - 32*me2*mpi2
c$$$     &    *s )
c$$$         boxc = boxc + D0ltd * (- 16*t**3 - 32*s*t**2 - 16*s**2*t - 48*
c$$$     &    mpi4*t - 16*mpi4*s - 48*me4*t - 32*me4*s + 48*mpi2*t**2 + 48*
c$$$     &    mpi2*s*t + 16*mpi2*mpi4 + 48*mpi2*me4 + 48*me2*t**2 + 64*me2*
c$$$     &    s*t + 16*me2*s**2 + 48*me2*mpi4 + 16*me2*me4 - 96*me2*mpi2*t
c$$$     &     - 48*me2*mpi2*s )
c$$$         boxc = boxc + C2ltd * (4*t**2+4*s*t + 4*mpi4 + 12*me4 - 8*mpi2
c$$$     &    *t - 16*me2*t - 8*me2*s + 16*me2*mpi2 )
c$$$         boxc = boxc + C1ltd * (4*t**2+4*s*t + 4*mpi4 + 4*me4 - 8*mpi2*
c$$$     &    t - 8*me2*t - 4*me2*s + 8*me2*mpi2 )
c$$$         boxc = boxc + C0ltd * (4*t**2+4*s*t + 4*mpi4 + 4*me4 - 8*mpi2*
c$$$     &        t - 8*me2*t - 4*me2*s + 8*me2*mpi2 )

         Z1_=mpi2*t
         Z2_=2*Z1_
         Z3_=s*t
         Z3_=Z3_ + mpi4
         Z4_=Z2_ - Z3_
         Z5_=t**2
         Z6_=Z4_ - Z5_
         Z7_=2*me2
         Z7_=Z6_*Z7_
         Z8_=mpi2 - t
         Z9_= - s + 2*Z8_
         Z10_= - me4*Z9_
         Z7_=Z7_ + Z10_
         Z7_=8*Z7_
         Z10_=Z6_*s
         Z9_=Z9_*me2
         Z11_=s*Z9_
         Z12_=me4*s
         Z13_= - 2*Z12_ + 4*Z10_ - 3*Z11_
         Z13_=4*Z13_
         Z14_=2*Z5_
         Z3_=Z14_ + Z3_
         Z15_=Z3_ - 3*Z1_
         Z15_=Z15_*s
         Z16_=Z5_ + 3*mpi4
         Z17_=Z16_*t
         Z18_=3*Z5_
         Z19_=Z18_ + mpi4
         Z20_=Z19_*mpi2
         Z15_=Z15_ + Z17_ - Z20_
         Z17_=Z5_ + mpi4
         Z20_=Z17_ - Z2_
         Z21_=3*s
         Z22_=Z21_ + 11*t - 9*mpi2
         Z22_=s*Z22_
         Z22_=8*Z20_ + Z22_
         Z22_=me2*Z22_
         Z21_=4*Z8_ - Z21_
         Z21_=me4*Z21_
         Z21_=Z21_ - 4*Z15_ + Z22_
         Z21_=4*Z21_
         Z10_= - Z10_ + Z11_ + Z12_
         Z11_=16*Z10_
         Z10_= - 20*Z10_
         Z12_=5*mpi4
         Z22_=s + 4*t
         Z23_=mpi2 + Z22_
         Z23_=s*Z23_
         Z2_=Z23_ + Z2_ - Z12_ + Z18_
         Z2_=me2*Z2_
         Z12_=Z12_ - Z5_
         Z12_=t*Z12_
         Z3_= - Z1_ - Z3_
         Z3_=s*Z3_
         Z16_= - mpi2*Z16_
         Z18_= - me2 + 2*s
         Z23_=Z18_ + 3*t
         Z24_= - mpi2 - Z23_
         Z24_=me4*Z24_
         Z2_=Z24_ + Z2_ + Z3_ + Z12_ + Z16_
         Z2_=4*Z2_
         Z3_= - 6*me4 + 8*Z6_ - 7*Z9_
         Z3_=8*Z3_
         Z12_=t + mpi2
         Z12_=s*Z12_
         Z12_=Z12_ - mpi4 + Z5_
         Z12_=me2*Z12_
         Z16_=3*mpi2
         Z24_= - me2 + s + t - Z16_
         Z24_=me4*Z24_
         Z12_=Z24_ + Z12_ - Z15_
         Z12_=16*Z12_
         Z24_= - mpi2 + 2*t
         Z24_=2*Z24_ + s
         Z24_=s*Z24_
         Z1_=Z24_ - 4*Z1_ + Z19_
         Z1_=me2*Z1_
         Z19_=2*mpi2
         Z5_=Z5_*Z19_
         Z4_= - Z14_ + Z4_
         Z4_=s*Z4_
         Z14_=Z19_ - Z23_
         Z14_=me4*Z14_
         Z17_= - t*Z17_
         Z1_=Z14_ + Z1_ + Z4_ + Z17_ + Z5_
         Z1_=16*Z1_
         Z4_= - Z16_ + Z22_
         Z4_=s*Z4_
         Z4_=3*Z20_ + Z4_
         Z4_=me2*Z4_
         Z5_=3*Z8_ - Z18_
         Z5_=me4*Z5_
         Z4_=Z5_ + Z4_ - Z15_
         Z4_=16*Z4_
         Z5_=3*me4 + 2*Z9_ - Z6_
         Z5_=4*Z5_
         Z6_=me4 + Z9_ - Z6_
         Z6_=4*Z6_
         boxc=D33ltd*Z7_ + D23ltd*Z13_ + D13ltd*Z21_ - D22ltd*Z11_ +
     &        D12ltd*Z10_ + D11ltd*Z2_ + D00ltd*Z3_ + D3ltd*Z12_
     &        - D2ltd*Z11_
     &        + D1ltd*Z1_ + D0ltd*Z4_ + C2ltd*Z5_ + C1ltd*Z6_
     &        + C0ltd*Z6_
*** STATS: original	5P 371M 143A : 524
*** STATS: optimized 0P 77M 67A : 144

         boxc = boxc * 2.d0 * (-1.d0) !2 is for interference with Born
     +        * ec6 * pi216m1/s
*                            (the phase is (+i)*(+i) the latter from Bornd)
*
      !!==============================================================!!
      !!                     GVMD IMPLEMENTATION                      !!
      !!==============================================================!!

      if (iffpi.eq.2) then

      !! Form factor selection

      if (what_ffpi.eq.'bwsum2') then
         call ffpi(s,ffpis,ffpis2)
         nres = 3
      elseif (what_ffpi.eq.'bwsum3') then
         call ffpi(s,ffpis,ffpis2)
         nres = 4
      elseif (what_ffpi.eq.'bwsum') then
         call ffpi(s,ffpis,ffpis2)
         nres = 6
         
      elseif (what_ffpi.eq.'cmd2') then
         call ffpi_bwsum2(s,ffpis,ffpis2)
         nres = 3
         if (igvmdwarning.eq.0) then
         print*, 'Warning: You are using two different form factors ',
     .           'for the loop amplitudes (bwsum2) and the tree-level ',
     .           'ones (cmd2), based on the same experimental data!'
         igvmdwarning = 1
         endif
      elseif (what_ffpi.eq.'cmd3') then
         call ffpi_bwsum3(s,ffpis,ffpis2)
         nres = 4
         if (igvmdwarning.eq.0) then
         print*, 'Warning: You are using two different form factors ',
     .           'for the loop amplitudes (bwsum3) and the tree-level ',
     .           'ones (cmd3), based on the same experimental data!'
         igvmdwarning = 1
         endif
      else
         call ffpi_bwsum2(s,ffpis,ffpis2)
         nres = 3
         if (igvmdwarning.eq.0) then
         print*, 'Warning: You are using two different form factors ',
     .           'for the loop amplitudes (bwsum2) and the tree-level ',
     .           'ones (', what_ffpi,'), based on data from different ',
     .           'experiments!'
         igvmdwarning = 1
         endif
      endif

      im = (0.d0,1.d0)
      do i=1,nres
      lv2(i) = m_n(i)**2 - im*m_n(i)*g_n(i)
      enddo

      !! Evaluation of loop functions (only with COLLIER)

#ifdef COLLIER

      !! 2-point functions

      call B0_cll(B0s,c1*s,c1*0.d0,c1*0.d0)
      call B0_cll(B0p,c1*mpi2,c1*0.d0,c1*mpi2)
      call DB0_cll(DB0p,c1*mpi2,c1*0.d0,c1*mpi2)

      do j=1,nres
         call B0_cll(B0lj(j),c1*0.d0,c1*lv2(j),c1*lv2(j))
         call B0_cll(B0sj(j),c1*s,c1*0.d0,c1*lv2(j))
         call B0_cll(B0pj(j),c1*mpi2,c1*lv2(j),c1*mpi2)
         call DB0_cll(DB0pj(j),c1*mpi2,c1*mpi2,c1*lv2(j))
      do k=1,nres
         call B0_cll(B0sjk(j,k),c1*s,c1*lv2(k),c1*lv2(j))
      enddo
      enddo

      !! 3-point functions

      call C0_cll(C0mms,c1*me2,c1*me2,c1*s,c1*0.d0,c1*me2,c1*0.d0)
      call C0_cll(C0pps,c1*mpi2,c1*mpi2,c1*s,c1*0.d0,c1*mpi2,c1*0.d0)
      call C0_cll(C0psp,c1*mpi2,c1*mpi2,c1*s,c1*mpi2,c1*0.d0,c1*mpi2)

      do j=1,nres
         call C0_cll(C0mmsj(j),c1*me2,c1*me2,c1*s,
     .                         c1*0.d0,c1*me2,c1*lv2(j))
         call C0_cll(C0ppsj(j),c1*mpi2,c1*mpi2,c1*s,
     .                         c1*0.d0,c1*mpi2,c1*lv2(j))
         call C0_cll(C0pspj(j),c1*mpi2,c1*mpi2,c1*s,
     .                         c1*mpi2,c1*lv2(j),c1*mpi2)
         call C0_cll(C0ppj(j),c1*0.d0,c1*mpi2,c1*mpi2,
     .                         c1*lv2(j),c1*lv2(j),c1*mpi2)
      do k=1,nres
         call C0_cll(C0mmsjk(j,k),c1*me2,c1*me2,c1*s,
     .                            c1*lv2(j),c1*me2,c1*lv2(k))
         call C0_cll(C0ppsjk(j,k),c1*mpi2,c1*mpi2,c1*s,
     .                            c1*lv2(j),c1*mpi2,c1*lv2(k))
      enddo
      enddo

      !! 4-point functions

      call D0_cll(D0d,c1*me2,c1*s,c1*mpi2,c1*t,c1*me2,c1*mpi2,
     .                c1*me2,c1*0.d0,c1*0.d0,c1*mpi2)
      call D0_cll(D0u,c1*me2,c1*s,c1*mpi2,c1*u,c1*me2,c1*mpi2,
     .                c1*me2,c1*0.d0,c1*0.d0,c1*mpi2)

      do j=1,nres
         call D0_cll(D0dj(j),c1*me2,c1*s,c1*mpi2,c1*t,c1*me2,c1*mpi2,
     .                       c1*me2,c1*lv2(j),c1*0.d0,c1*mpi2)
         call D0_cll(D0uj(j),c1*me2,c1*s,c1*mpi2,c1*u,c1*me2,c1*mpi2,
     .                       c1*me2,c1*lv2(j),c1*0.d0,c1*mpi2)
      do k=1,nres
      call D0_cll(D0djk(j,k),c1*me2,c1*s,c1*mpi2,c1*t,c1*me2,c1*mpi2,
     .                  c1*me2,c1*lv2(j),c1*lv2(k),c1*mpi2)
      call D0_cll(D0ujk(j,k),c1*me2,c1*s,c1*mpi2,c1*u,c1*me2,c1*mpi2,
     .                  c1*me2,c1*lv2(j),c1*lv2(k),c1*mpi2)
      enddo
      enddo

#endif

#ifdef LOOPTOOLS

      !! 2-point functions

      B0s = B0C(c1*s,c1*0.d0,c1*0.d0)
      B0p = B0C(c1*mpi2,c1*0.d0,c1*mpi2)
      DB0p = DB0C(c1*mpi2,c1*0.d0,c1*mpi2)

      do j=1,nres
         B0lj(j) = B0C(c1*0.d0,c1*lv2(j),c1*lv2(j))
         B0sj(j) = B0C(c1*s,c1*0.d0,c1*lv2(j))
         B0pj(j) = B0C(c1*mpi2,c1*lv2(j),c1*mpi2)
         DB0pj(j) = B0C(c1*mpi2,c1*mpi2,c1*lv2(j))
      do k=1,nres
         B0sjk(j,k) = B0C(c1*s,c1*lv2(k),c1*lv2(j))
      enddo
      enddo

      !! 3-point functions

      C0mms = C0C(c1*me2,c1*me2,c1*s,c1*0.d0,c1*me2,c1*0.d0)
      C0pps = C0C(c1*mpi2,c1*mpi2,c1*s,c1*0.d0,c1*mpi2,c1*0.d0)
      C0psp = C0C(c1*mpi2,c1*mpi2,c1*s,c1*mpi2,c1*0.d0,c1*mpi2)

      do j=1,nres
         C0mmsj(j) = C0C(c1*me2,c1*me2,c1*s,
     .                   c1*0.d0,c1*me2,c1*lv2(j))
         C0ppsj(j) = C0C(c1*mpi2,c1*mpi2,c1*s,
     .                   c1*0.d0,c1*mpi2,c1*lv2(j))
         C0pspj(j) = C0C(c1*mpi2,c1*mpi2,c1*s,
     .                   c1*mpi2,c1*lv2(j),c1*mpi2)
         C0ppj(j) = C0C(c1*0.d0,c1*mpi2+c1*1d-12,c1*mpi2,
     .                   c1*lv2(j),c1*lv2(j),c1*mpi2)
      do k=1,nres
         C0mmsjk(j,k) = C0C(c1*me2,c1*me2,c1*s,
     .                      c1*lv2(j),c1*me2,c1*lv2(k))
         C0ppsjk(j,k) = C0C(c1*mpi2,c1*mpi2,c1*s,
     .                      c1*lv2(j),c1*mpi2,c1*lv2(k))
      enddo
      enddo

      !! 4-point functions

      D0d = D0C(c1*me2,c1*s,c1*mpi2,c1*t,c1*me2,c1*mpi2,
     .          c1*me2,c1*0.d0,c1*0.d0,c1*mpi2)
      D0u = D0C(c1*me2,c1*s,c1*mpi2,c1*u,c1*me2,c1*mpi2,
     .          c1*me2,c1*0.d0,c1*0.d0,c1*mpi2)

      do j=1,nres
      D0dj(j) = D0C(c1*me2,c1*s,c1*mpi2,c1*t,c1*me2,c1*mpi2,
     .              c1*me2,c1*lv2(j),c1*0.d0,c1*mpi2)
      D0uj(j) = D0C(c1*me2,c1*s,c1*mpi2,c1*u,c1*me2,c1*mpi2,
     .              c1*me2,c1*lv2(j),c1*0.d0,c1*mpi2)
      do k=1,nres
      D0djk(j,k) = D0C(c1*me2,c1*s,c1*mpi2,c1*t,c1*me2,c1*mpi2,
     .                 c1*me2,c1*lv2(j),c1*lv2(k),c1*mpi2)
      D0ujk(j,k) = D0C(c1*me2,c1*s,c1*mpi2,c1*u,c1*me2,c1*mpi2,
     .                 c1*me2,c1*lv2(j),c1*lv2(k),c1*mpi2)
      enddo
      enddo

#endif

      rcounterterms_pi_gvmd = (0.d0,0.d0)
      vertex_pi_gvmd = (0.d0,0.d0)
      vertex_seagull_gvmd = (0.d0,0.d0)
      boxt_gvmd = (0.d0,0.d0)
      boxd_gvmd = (0.d0,0.d0)
      boxc_gvmd = (0.d0,0.d0)

      do j=1,nres
      do k=1,nres

      !!=====!!
      !! FSR !!
      !!=====!!

      !! Case 1: Same vector meson resonance in the loop

      if (j.eq.k) then

      !! Some explicit loop functions

      C0jj(j) = log(mpi2/lv2(j))/(2*mpi2) + ((lv2(j) - 2*mpi2)*
     &   log((lv2(j) + sqrt(lv2(j)*(lv2(j) - 4*mpi2)))/
     &   (2*sqrt(mpi2*lv2(j)))))/(sqrt(lv2(j)*(lv2(j) - 4*mpi2))*mpi2)

      C1jj(j) = (2*mpi2 + (lv2(j) - 2*mpi2)*log(mpi2/lv2(j)) +
     &   (2*(lv2(j)**2 - 4*lv2(j)*mpi2 + 2*mpi2**2)*
     &   log((lv2(j) + sqrt(lv2(j)*(lv2(j) - 4*mpi2)))/(2*sqrt(mpi2*
     &   lv2(j)))))/sqrt(lv2(j)*(lv2(j) - 4*mpi2)))/(4*mpi4)

      C2jj(j) = -1/mpi2 + ((-lv2(j) + mpi2)*log(mpi2/lv2(j)))/(2*mpi4) -
     &   (lv2(j)*(lv2(j) - 3*mpi2)*log((lv2(j) + sqrt(lv2(j)*
     &   (lv2(j) - 4*mpi2)))/(2*sqrt(mpi2*lv2(j))))) /
     &   (sqrt(lv2(j)*(lv2(j) - 4*mpi2))*mpi2**2)

      D0jj(j) = ((lv2(j) - 4*mpi2)*sqrt(s*(-4*mpi2 + s))*
     &   log(c1*((2*mpi2 - s + sqrt(s*(-4*mpi2 + s)))/(2*mpi2))) +
     &   2*sqrt(lv2(j)*(lv2(j) - 4*mpi2))*s*log((lv2(j) +
     &   sqrt(lv2(j)*(lv2(j) - 4*mpi2)))/(2*sqrt(mpi2*lv2(j)))))/
     &   ((lv2(j) - 4*mpi2)*s*(lv2(j) - 4*mpi2 + s)*lv2(j))

      D2jj(j) = C0pspj(j)/(-4*mpi2 + s) + log(mpi2/lv2(j))/(8*mpi2**2 -
     &   2*mpi2*s) + log(c1*((2*mpi2 - s + sqrt(s*(-4*mpi2 + s))) /
     &   (2*mpi2)))/(sqrt(s*(-4*mpi2 + s))*(lv2(j) - 4*mpi2 + s)) +
     &   (lv2(j)*(lv2(j)**2 - 8*lv2(j)*mpi2 + 8*mpi2**2 + lv2(j)*s -
     &   2*mpi2*s)*log((lv2(j) + sqrt(lv2(j)*(lv2(j) - 4*mpi2))) /
     &   (2*sqrt(mpi2*lv2(j)))))/(sqrt(lv2(j)*(lv2(j) - 4*mpi2))*mpi2*
     &   (4*mpi2 - s)*(lv2(j) - 4*mpi2 + s)*lv2(j))

      rcounterterms_pi_gvmd = rcounterterms_pi_gvmd -
     &   (16*a_n(j)**2*alpha**3*pi*(4*mpi2*s - s**2 + (t - u)**2)*
     &   (C0jj(j)*lv2(j)**3 + lv2(j)**2*(1 - 4*C0jj(j)*mpi2) +
     &   2*lv2(j)*mpi2*(-2 + B0p - B0pj(j) +
     &   (C0jj(j) + 2*DB0p - 2*DB0pj(j))*mpi2) +
     &   4*mpi2**2*(1 - B0p + B0pj(j) - 2*DB0p*mpi2 + 2*DB0pj(j)*mpi2) +
     &   lv2(j)*mpi2*log(mpi2/lv2(j))))/((lv2(j) - 2*mpi2)*mpi2*s**2)

      vertex_pi_gvmd = vertex_pi_gvmd - (16*a_n(j)**2*alpha**3*pi*
     &   (-2*lv2(j)*s + (B0pj(j)*(lv2(j)*(4*mpi2 - s)*((2*mpi2 - s)*s +
     &   (t - u)**2) - 2*mpi2*(8*mpi2 - 3*s)*((4*mpi2 - s)*s +
     &   (t - u)**2)))/(4*mpi2**2 - mpi2*s) + D0jj(j)*lv2(j)*(lv2(j) -
     &   4*mpi2 + 2*s)*(s*(2*lv2(j) - 4*mpi2 + s) - (t - u)**2) +
     &   2*C0jj(j)*lv2(j)*((-2*lv2(j) + 8*mpi2 - 3*s)*s + (t - u)**2) +
     &   2*C2jj(j)*lv2(j)*(2*mpi2*s - s**2 + (t - u)**2) -
     &   ((1 + B0lj(j))*lv2(j)*(2*mpi2*s - s**2 + (t - u)**2))/mpi2 +
     &   2*C0psp*(2*mpi2 - s)*(4*mpi2*s - s**2 + (t - u)**2) +
     &   (2*B0p*(8*mpi2 - 3*s)*(4*mpi2*s - s**2 + (t - u)**2))/
     &   (4*mpi2 - s) - (2*C0pspj(j)*(lv2(j)**2*s*(-4*mpi2 + s) +
     &   (8*mpi2**2 - 6*mpi2*s + s**2)*(4*mpi2*s - s**2 + (t - u)**2) -
     &   2*lv2(j)*(2*mpi2 - s)*(t - u)**2))/(4*mpi2 - s) -
     &   2*D2jj(j)*lv2(j)*(lv2(j) - 4*mpi2 + 2*s)*(t - u)**2))/s**2

      vertex_seagull_gvmd = vertex_seagull_gvmd +
     &   (-8*a_n(j)**2*alpha**3*pi*(2*C0jj(j)*lv2(j)**2*mpi2*s +
     &   lv2(j)*(2*mpi2*(1 + B0lj(j) - B0pj(j) - 6*C0jj(j)*mpi2 +
     &   4*C1jj(j)*mpi2)*s - (1 + B0lj(j) - B0pj(j))*(s**2 -
     &   (t - u)**2) + 2*mpi2*(s + (C0jj(j) - 2*C1jj(j))*s**2 -
     &   (C0jj(j) - 2*C1jj(j))*(t - u)**2)) + 4*(-B0p + B0pj(j))*
     &   mpi2*(4*mpi2*s - s**2 + (t - u)**2)))/(mpi2*s**2)

      endif

      !! Case 2: Different vector meson resonances in the loop

      if (j.ne.k) then

      lv2diff = lv2(j) - lv2(k)

      !! Pion field renormalisation

      rcounterterms_pi_gvmd = rcounterterms_pi_gvmd -
     &   (16*a_n(j)*a_n(k)*alpha**3*(2*B0p*lv2(j) - 2*B0pj(k)*lv2(j) -
     &   2*B0p*lv2(k) + 2*B0pj(j)*lv2(k) + DB0pj(k)*lv2(j)*lv2(k) -
     &   DB0pj(j)*lv2(k)*(lv2(j) - 4*mpi2) - 4*DB0pj(k)*lv2(j)*mpi2 +
     &   4*DB0p*lv2diff*mpi2)*pi*(4*mpi2*s - s**2 + (t - u)**2))/
     &   (lv2diff*s**2)

      !! Pion vertex

      vertex_pi_gvmd = vertex_pi_gvmd + (16*a_n(j)*a_n(k)*alpha**3*pi*
     &   (-2*B0p*lv2diff*mpi2*(8*mpi2 - 3*s) - (1 + B0lj(j))*lv2(j)*
     &   lv2(k)*(4*mpi2 - s) + (1 + B0lj(k))*lv2(j)*lv2(k)*(4*mpi2 -
     &   s) - 2*C0psp*lv2diff*mpi2*(2*mpi2 - s)*(4*mpi2 - s) -
     &   C0pspj(j)*lv2(k)*mpi2*(-lv2(j) + 4*mpi2 - 2*s)*(-2*lv2(j) +
     &   4*mpi2 - s) + C0pspj(k)*lv2(j)*mpi2*(-lv2(k) + 4*mpi2 - 2*s)*
     &   (-2*lv2(k) + 4*mpi2 - s) - B0pj(j)*lv2(k)*(16*mpi4 + lv2(j)*s -
     &   6*mpi2*(lv2(j) + s)) + B0pj(k)*lv2(j)*(16*mpi4 + lv2(k)*s -
     &   6*mpi2*(lv2(k) + s)))*(4*mpi2*s - s**2 + (t - u)**2))/
     &   (lv2diff*mpi2*(4*mpi2 - s)*s**2)

      !! Seagull vertex

      vertex_seagull_gvmd = vertex_seagull_gvmd +   !- CHANGED MARCO 03/03/26
     &   (8*a_n(j)*a_n(k)*alpha**3*pi*(B0lj(j)*lv2(j)*lv2(k) - B0lj(k)*
     &   lv2(j)*lv2(k) - B0pj(j)*lv2(j)*lv2(k) + B0pj(k)*lv2(j)*lv2(k) +
     &   4*B0p*lv2(j)*mpi2 - 4*B0pj(k)*lv2(j)*mpi2 - 4*B0p*lv2(k)*mpi2 +
     &   4*B0pj(j)*lv2(k)*mpi2)*(4*mpi2*s - s**2 + (t - u)**2))/
     &   (lv2diff*mpi2*s**2)

      endif

      !! To remove GVMD correction to FSR if iFSRdisp = 0
      if(iFSRdisp.eq.0) then
      vertex_pi_gvmd = vertex_pi
      vertex_seagull_gvmd = vertex_seagull1
      rcounterterms_pi_gvmd = rcounterterms_pi
      endif
      !!=====!!
      !! IFI !!
      !!=====!!

      !! Seagull box

      boxt_gvmd = boxt_gvmd - (128*a_n(j)*a_n(k)*alpha**3*me2*pi*
     &   (-2*B0s + 2*B0sj(j) + 2*B0sj(k) - 2*B0sjk(j,k) + C0mms*(4*me2 -
     &   2*s) - C0mmsj(j)*(lv2(j) + 4*me2 - 2*s) - C0mmsj(k)*(lv2(k) +
     &   4*me2 -2*s) + C0mmsjk(j,k)*(lv2(j) + lv2(k) + 4*me2 - 2*s))*
     &   (t - u))/((4*me2 - s)*s)

      !! Direct box

      boxd_gvmd = boxd_gvmd - (32*a_n(j)*a_n(k)*alpha**3*pi*
     &   (2*B0s*(me4 + mpi4 + me2*(2*mpi2 - s - 2*t) - 2*mpi2*t +
     &   t*(s + t)) - 2*B0sj(j)*(me4 + mpi4 + me2*(2*mpi2 - s - 2*t) -
     &   2*mpi2*t + t*(s + t)) - 2*B0sj(k)*(me4 + mpi4 + me2*(2*mpi2 -
     &   s - 2*t) - 2*mpi2*t + t*(s + t)) + 2*B0sjk(j,k)*(me4 + mpi4 +
     &   me2*(2*mpi2 - s - 2*t) - 2*mpi2*t + t*(s + t)) + 2*C0mms*
     &   (4*mpi2 - s)*(2*me4 + 2*mpi4 + me2*(4*mpi2 - s - 4*t) -
     &   4*mpi2*t + t*(s + 2*t)) - 2*C0mmsj(j)*(4*mpi2 - s)*(2*me4 +
     &   2*mpi4 + me2*(4*mpi2 - s - 4*t) - 4*mpi2*t + t*(s + 2*t)) -
     &   2*C0mmsj(k)*(4*mpi2 - s)*(2*me4 + 2*mpi4 + me2*(4*mpi2 - s -
     &   4*t) - 4*mpi2*t + t*(s + 2*t)) + 2*C0mmsjk(j,k)*(4*mpi2 - s)*
     &   (2*me4 + 2*mpi4 + me2*(4*mpi2 - s - 4*t) - 4*mpi2*t + t*(s +
     &   2*t)) + 2*D0d*(4*mpi2 - s)*(me2 + mpi2 - t)*(2*me4 + 2*mpi4 +
     &   me2*(4*mpi2 - s - 4*t) - 4*mpi2*t + t*(s + 2*t)) - 2*D0dj(j)*
     &   (4*mpi2 - s)*(me2 + mpi2 - t)*(2*me4 + 2*mpi4 - 4*mpi2*t +
     &   t*(lv2(j) + s + 2*t) - me2*(lv2(j) - 4*mpi2 + s + 4*t)) -
     &   2*D0dj(k)*(4*mpi2 - s)*(me2 + mpi2 - t)*(2*me4 + 2*mpi4 -
     &   4*mpi2*t + t*(lv2(k) + s + 2*t) - me2*(lv2(k) - 4*mpi2 + s +
     &   4*t)) + 2*D0djk(j,k)*(4*mpi2 - s)*(me2 + mpi2 - t)*
     &   (2*me4 + 2*mpi4 - 4*mpi2*t + t*(lv2(j) + lv2(k) + s + 2*t) -
     &   me2*(lv2(j) + lv2(k) - 4*mpi2 + s + 4*t)) + C0pps*(me4*
     &   (16*mpi2 - 3*s) + mpi4*(s - 16*t) - s*t*(s + 3*t) + 2*mpi2*t*
     &   (5*s + 8*t) + me2*(16*mpi4 + s*(s + 6*t) - 2*mpi2*(5*s +
     &   16*t))) - C0ppsj(j)*(-(me4*(lv2(j) - 16*mpi2 + 3*s)) +
     &   mpi4*(-lv2(j) + s - 16*t) + 2*mpi2*t*(lv2(j) + 5*s + 8*t) -
     &   t*(lv2(j)*(s + t) + s*(s + 3*t)) + me2*(16*mpi4 + lv2(j)*s +
     &   s**2 + 2*lv2(j)*t + 6*s*t - 2*mpi2*(lv2(j) + 5*s + 16*t))) -
     &   C0ppsj(k)*(-(me4*(lv2(k) - 16*mpi2 + 3*s)) + mpi4*(-lv2(k) +
     &   s - 16*t) + 2*mpi2*t*(lv2(k) + 5*s + 8*t) - t*(lv2(k)*(s + t) +
     &   s*(s + 3*t)) + me2*(16*mpi4 + lv2(k)*s + s**2 + 2*lv2(k)*t +
     &   6*s*t - 2*mpi2*(lv2(k) + 5*s + 16*t))) + C0ppsjk(j,k)*
     &   (-(me4*(lv2(j) + lv2(k) - 16*mpi2 + 3*s)) + 2*mpi2*t*(lv2(j) +
     &   lv2(k) + 5*s + 8*t) - mpi4*(lv2(j) + lv2(k) - s + 16*t) -
     &   t*(s**2 + (lv2(j) + lv2(k))*t + s*(lv2(j) + lv2(k) + 3*t)) +
     &   me2*(16*mpi4 + s**2 + 2*(lv2(j) + lv2(k))*t + s*(lv2(j) +
     &   lv2(k) + 6*t) - 2*mpi2*(lv2(j) + lv2(k) + 5*s + 16*t)))))/
     &   ((4*mpi2 - s)*s)

      !! Crossed box

      boxc_gvmd = boxc_gvmd - (32*a_n(j)*a_n(k)*alpha**3*pi*
     &   (-2*B0s*(me4 + mpi4 + me2*(2*mpi2 - s - 2*u) - 2*mpi2*u +
     &   u*(s + u)) + 2*B0sj(j)*(me4 + mpi4 + me2*(2*mpi2 - s - 2*u) -
     &   2*mpi2*u + u*(s + u)) + 2*B0sj(k)*(me4 + mpi4 + me2*(2*mpi2 -
     &   s - 2*u) - 2*mpi2*u + u*(s + u)) - 2*B0sjk(j,k)*(me4 + mpi4 +
     &   me2*(2*mpi2 - s - 2*u) - 2*mpi2*u + u*(s + u)) - 2*C0mms*
     &   (4*mpi2 - s)*(2*me4 + 2*mpi4 + me2*(4*mpi2 - s - 4*u) -
     &   4*mpi2*u + u*(s + 2*u)) + 2*C0mmsj(j)*(4*mpi2 - s)*(2*me4 +
     &   2*mpi4 + me2*(4*mpi2 - s - 4*u) - 4*mpi2*u + u*(s + 2*u)) +
     &   2*C0mmsj(k)*(4*mpi2 - s)*(2*me4 + 2*mpi4 + me2*(4*mpi2 - s -
     &   4*u) - 4*mpi2*u + u*(s + 2*u)) - 2*C0mmsjk(j,k)*(4*mpi2 - s)*
     &   (2*me4 + 2*mpi4 + me2*(4*mpi2 - s - 4*u) - 4*mpi2*u + u*(s +
     &   2*u)) - 2*D0u*(4*mpi2 - s)*(me2 + mpi2 - u)*(2*me4 + 2*mpi4 +
     &   me2*(4*mpi2 - s - 4*u) - 4*mpi2*u + u*(s + 2*u)) + 2*D0uj(j)*
     &   (4*mpi2 - s)*(me2 + mpi2 - u)*(2*me4 + 2*mpi4 - 4*mpi2*u +
     &   u*(lv2(j) + s + 2*u) - me2*(lv2(j) - 4*mpi2 + s + 4*u)) +
     &   2*D0uj(k)*(4*mpi2 - s)*(me2 + mpi2 - u)*(2*me4 + 2*mpi4 -
     &   4*mpi2*u + u*(lv2(k) + s + 2*u) - me2*(lv2(k) - 4*mpi2 + s +
     &   4*u)) - 2*D0ujk(j,k)*(4*mpi2 - s)*(me2 + mpi2 - u)*(2*me4 +
     &   2*mpi4 - 4*mpi2*u + u*(lv2(j) + lv2(k) + s + 2*u) -
     &   me2*(lv2(j) + lv2(k) - 4*mpi2 + s + 4*u)) - C0pps*(me4*
     &   (16*mpi2 - 3*s) + mpi4*(s - 16*u) - s*u*(s + 3*u) + 2*mpi2*u*
     &   (5*s + 8*u) + me2*(16*mpi4 + s*(s + 6*u) - 2*mpi2*(5*s +
     &   16*u))) + C0ppsj(j)*(-(me4*(lv2(j) - 16*mpi2 + 3*s)) + mpi4*
     &   (-lv2(j) + s - 16*u) + 2*mpi2*u*(lv2(j) + 5*s + 8*u) -
     &   u*(lv2(j)*(s + u) + s*(s + 3*u)) + me2*(16*mpi4 + lv2(j)*s +
     &   s**2 + 2*lv2(j)*u + 6*s*u - 2*mpi2*(lv2(j) + 5*s + 16*u))) +
     &   C0ppsj(k)*(-(me4*(lv2(k) - 16*mpi2 + 3*s)) + mpi4*(-lv2(k) +
     &   s - 16*u) + 2*mpi2*u*(lv2(k) + 5*s + 8*u) - u*(lv2(k)*(s + u) +
     &   s*(s + 3*u)) + me2*(16*mpi4 + lv2(k)*s + s**2 + 2*lv2(k)*u +
     &   6*s*u - 2*mpi2*(lv2(k) + 5*s + 16*u))) + C0ppsjk(j,k)*(me4*
     &   (lv2(j) + lv2(k) - 16*mpi2 + 3*s) - 2*mpi2*u*(lv2(j) + lv2(k) +
     &   5*s + 8*u) + mpi4*(lv2(j) + lv2(k) - s + 16*u) + u*(s**2 +
     &   (lv2(j) + lv2(k))*u + s*(lv2(j) + lv2(k) + 3*u)) -
     &   me2*(16*mpi4 + s**2 + 2*(lv2(j) + lv2(k))*u + s*(lv2(j) +
     &   lv2(k) + 6*u) - 2*mpi2*(lv2(j) + lv2(k) + 5*s + 16*u)))))/
     &   ((4*mpi2 - s)*s)

      enddo
      enddo

      boxt_gvmd = boxt_gvmd/ffpis
      boxd_gvmd = boxd_gvmd/ffpis
      boxc_gvmd = boxc_gvmd/ffpis

      endif
      endif

      cISR= vertex_e+rcounterterms_e
      cFSR= vertex_pi
     +      +vertex_seagull1
     +      +vertex_seagull2
     +      +rcounterterms_pi

      cIFI= boxd+boxc+vertex_2photons

      if (iffpi.eq.2) then
         cFSR = vertex_pi_gvmd
     &        + vertex_seagull_gvmd
     &        + vertex_seagull_gvmd
     &        + rcounterterms_pi_gvmd

         cIFI = boxd_gvmd + boxc_gvmd + boxt_gvmd
      endif

      vISR= dble(cISR)
      vFSR= dble(cFSR)
      vIFI= dble(cIFI)

*
c      call realsoft(p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,mpi2,me2,
c     +              rsoftISR,rsoftFSR,rsoftIFI)
*
      call soft_integral_pions(egmin,p1,p2,p3,p4,rsofttot,rsoftISR,
     .     rsoftFSR,rsoftIFI)

      rsoftISR= rsoftISR*born
      rsoftFSR= rsoftFSR*born
      rsoftIFI= rsoftIFI*born
*
      svISR= vISR+rsoftISR
      if(fs.eq.'pp') then
         svFSR= vFSR+rsoftFSR
         svIFI= vIFI+rsoftIFI
      elseif(fs.eq.'pi') then
         svFSR= 0.d0
         svIFI= 0.d0
      endif


      aintsv_pions= svISR+svFSR+svIFI
*
      vpols = vpolc(s)
*
      aintsv_pions = born + aintsv_pions !to match the meaning of
*                                        !sv in main.F (sdif = sdif * sv)
      AINTSV_PIONS = AINTSV_PIONS * abs(vpols)**2 * betafs !sqrt(1 - 4 m^2/s) from phase space
      AINTSV_PIONS = AINTSV_PIONS / pi / s * 0.0078125d0 ! .0078125 = 1/128
*
* including the form factor (completely factorized over born and sv)
*
      ffpis2 = 1.d0
      if (iffpi.gt.0) call ffpi(s,ffpis,ffpis2)
*
      aintsv_pions = aintsv_pions*ffpis2
*
      return
      end
*********************************************************
      DOUBLE PRECISION FUNCTION DDILOGhere(X)
      DOUBLE PRECISION X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
      DOUBLE PRECISION C(0:18),H,ALFA,B0,B1,B2
      DATA ZERO /0.0D0/, ONE /1.0D0/
      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/
      DATA PI3 /3.28986 81336 96453D0/, PI6 /1.64493 40668 48226D0/
      DATA C( 0) / 0.42996 69356 08137 0D0/
      DATA C( 1) / 0.40975 98753 30771 1D0/
      DATA C( 2) /-0.01858 84366 50146 0D0/
      DATA C( 3) / 0.00145 75108 40622 7D0/
      DATA C( 4) /-0.00014 30418 44423 4D0/
      DATA C( 5) / 0.00001 58841 55418 8D0/
      DATA C( 6) /-0.00000 19078 49593 9D0/
      DATA C( 7) / 0.00000 02419 51808 5D0/
      DATA C( 8) /-0.00000 00319 33412 7D0/
      DATA C( 9) / 0.00000 00043 45450 6D0/
      DATA C(10) /-0.00000 00006 05784 8D0/
      DATA C(11) / 0.00000 00000 86121 0D0/
      DATA C(12) /-0.00000 00000 12443 3D0/
      DATA C(13) / 0.00000 00000 01822 6D0/
      DATA C(14) /-0.00000 00000 00270 1D0/
      DATA C(15) / 0.00000 00000 00040 4D0/
      DATA C(16) /-0.00000 00000 00006 1D0/
      DATA C(17) / 0.00000 00000 00000 9D0/
      DATA C(18) /-0.00000 00000 00000 1D0/

      DDILOGhere = 0.d0
      IF(X .EQ. ONE) THEN
       DDILOGhere=PI6
       RETURN
      ELSE IF(X .EQ. MONE) THEN
       DDILOGhere=MALF*PI6
       RETURN
      END IF
      T=-X
      IF(T .LE. MTWO) THEN
       Y=MONE/(ONE+T)
       S=ONE
       A=-PI3+HALF*(DLOG(-T)**2-LOG(ONE+ONE/T)**2)
      ELSE IF(T .LT. MONE) THEN
       Y=MONE-T
       S=MONE
       A=DLOG(-T)
       A=-PI6+A*(A+DLOG(ONE+ONE/T))
      ELSE IF(T .LE. MALF) THEN
       Y=(MONE-T)/T
       S=ONE
       A=DLOG(-T)
       A=-PI6+A*(MALF*A+DLOG(ONE+T))
      ELSE IF(T .LT. ZERO) THEN
       Y=-T/(ONE+T)
       S=MONE
       A=HALF*DLOG(ONE+T)**2
      ELSE IF(T .LE. ONE) THEN
       Y=T
       S=ONE
       A=ZERO
      ELSE
       Y=ONE/T
       S=MONE
       A=PI6+HALF*DLOG(T)**2
      END IF

      H=Y+Y-ONE
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO I = 18,0,-1
      B0=C(I)+ALFA*B1-B2
      B2=B1
    1 B1=B0
      enddo
      DDILOGhere=-(S*(B0-H*B2)+A)
      RETURN
      END


***
      function soft_integral_gg(egmax,p1,p2,p3,p4)
! this should work also for massless momenta...
      implicit double precision (a-h,l,m,o-z)
      parameter (npart = 4)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension p(0:3),q(0:3),pmat(npart,0:3)
      double precision mass(5)
      dimension masses(npart),eta(npart)
      common/sffcommonNEW/ao2pi,masses,eta,alog2,ifirst
      data ifirst/0/
      common/parameters/ame,ammu,convfac,alpha,pi
      double precision phmass
      common/softphotonmass/phmass

* #include "invariants.h"

! ch are the charges of the field (not anti-field!!)
! Use this convention: the integer factor in front of the charge
! must be
! -1 --> for incoming particle
! -1 --> for outgoing anti-particle
! +1 --> for outgoing particle
! +1 --> for incoming anti-particle

      if (ifirst.eq.0) then

         eta(1) = -1.d0 * (1.d0)
         eta(2) = -1.d0 * (-1.d0)
         eta(3) =  0.d0
         eta(4) =  0.d0

         masses(1) = ame*ame
         masses(2) = masses(1)
         masses(3) = 0.d0
         masses(4) = 0.d0

         ao2pi = alpha/2.d0/pi
         alog2 = log(2.d0)

c         phmass = sqrt(getlambda())

         ifirst = 1
      endif

      linfra    = log(egmax/phmass)
      linfrapl2 = linfra + alog2

      pmat(1,:) = p1
      pmat(2,:) = p2
      pmat(3,:) = p3
      pmat(4,:) = p4

      npartl = npart

! off-diagonal contributions
      softint = 0.d0
      do i = 1,npartl-1
         do j = i+1,npartl
            etaij = eta(i) * eta(j)

            if (abs(etaij).gt.1.d-3) then
               call rescale_momenta_gg(npart,i,j,pmat,masses,p,q,rho)
               q2 = masses(j)
               vl = 0.5d0*(rho*rho*masses(i) - q2)
               v  = vl/(p(0) - q(0))
               arglog = 1.d0+2.d0*vl/q2

               terminfra = 0.d0
               if (arglog.gt.0.d0) then
                  terminfra = linfrapl2*log(arglog)
                  tot = terminfra + termfinite_gg(p,q,v,rho,i,j)
               else
                  tot = 0.d0
               endif

               tot = -tot * 2.d0*dot(p,q)/vl* etaij
               tot =  tot * 2.d0 ! this is the double product when
                                ! squaring the eikonal
               softint  = softint + tot
            endif
         enddo
      enddo
      soft_integral_gg = softint

! diagonal contributions
      softint = 0.d0

c     only boxes
c      goto 666

      do  i = 1,npartl
         etaii = eta(i)*eta(i)
         if (etaii.gt.1.d-3) then
            betai=1.d0-masses(i)/pmat(i,0)/pmat(i,0)
            betai=sqrt(betai)

            if (betai.gt.1.d-12) then
               bterm =  log((1.d0-betai)/(1.d0+betai))/betai
            else
               bterm = -2.d0
            endif
            term = 2.d0*linfrapl2  + bterm
            term = -2.d0 * term * etaii
            softint = softint + term
         endif
      enddo
! divide by 2 to obtain exactly what is inside {} in Dittmaier & Kraemer
! paper (hep-ph/0109062)
 666  soft_integral_gg = 0.5d0*(soft_integral_gg + softint)*ao2pi
      return
      end
***********************************************************************************************
      subroutine soft_integral_pions(egmax,p1,p2,p3,p4,rsoft,risr,rfsr,
     .     rifi)
! this should work also for massless momenta...
      implicit double precision (a-h,l,m,o-z)
      parameter (npart = 4)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension p(0:3),q(0:3),pmat(npart,0:3)
      dimension masses(npart),eta(npart)
      common/sffcommonpions/ao2pi,masses,eta,alog2,ifirst
      data ifirst/0/
      common/hparam/ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      double precision phmass
      common/softphotonmass/phmass

* #include "invariants.h"

! ch are the charges of the field (not anti-field!!)
! Use this convention: the integer factor in front of the charge
! must be
! -1 --> for incoming particle
! -1 --> for outgoing anti-particle
! +1 --> for outgoing particle
! +1 --> for incoming anti-particle

      if (ifirst.eq.0) then

         eta(1) = -1.d0 * (1.d0)
         eta(2) = -1.d0 * (-1.d0)
         eta(3) =  -eta(1)
         eta(4) =  -eta(2)

         masses(1) = ame*ame
         masses(2) = masses(1)
         masses(3) = ampi*ampi
         masses(4) = masses(3)

         ao2pi = alpha/2.d0/pi
         alog2 = log(2.d0)

         ifirst = 1
      endif

      linfra    = log(egmax/phmass)
      linfrapl2 = linfra + alog2

      pmat(1,:) = p1
      pmat(2,:) = p2
      pmat(3,:) = p3
      pmat(4,:) = p4

      npartl = npart

      rsoft = 0.d0
      risr  = 0.d0
      rfsr  = 0.d0
      rifi  = 0.d0

! off-diagonal contributions
      softint = 0.d0
      do i = 1,npartl-1
         do j = i+1,npartl
            etaij = eta(i) * eta(j)
            if (abs(etaij).gt.1.d-3) then
               call rescale_momenta_gg(npart,i,j,pmat,masses,p,q,rho)
               q2 = masses(j)
               vl = 0.5d0*(rho*rho*masses(i) - q2)
               v  = vl/(p(0) - q(0))
               arglog = 1.d0+2.d0*vl/q2

               terminfra = 0.d0
               if (arglog.gt.0.d0) then
                  terminfra = linfrapl2*log(arglog)
                  tot = terminfra + termfinite_pions(p,q,v,rho,i,j)
               else
                  tot = 0.d0
               endif

               tot = -tot * 2.d0*dot(p,q)/vl* etaij
               tot =  tot * 2.d0 ! this is the double product when
                                ! squaring the eikonal
               softint  = softint + tot

               if (i.eq.1.and.j.eq.2) then
                  risr = risr + tot
               elseif (i.eq.3.and.j.eq.4) then
                  rfsr = rfsr + tot
               else
                  rifi = rifi + tot
               endif
            endif
         enddo
      enddo
      soft_integral = softint


! diagonal contributions
      softint = 0.d0
      do  i = 1,npartl
         etaii = eta(i)*eta(i)
         if (etaii.gt.1.d-3) then
            betai=1.d0-masses(i)/pmat(i,0)/pmat(i,0)
            betai=sqrt(betai)

            if (betai.gt.1.d-12) then
               bterm =  log((1.d0-betai)/(1.d0+betai))/betai
            else
               bterm = -2.d0
            endif
            term = 2.d0*linfrapl2  + bterm
            term = -2.d0 * term * etaii
            softint = softint + term

            if (i.le.2) then
               risr = risr + term
            else
               rfsr = rfsr + term
            endif
         endif
      enddo

! divide by 2 to obtain exactly what is inside {} in Dittmaier & Kraemer
! paper (hep-ph/0109062)
 666  soft_integral = 0.5d0*(soft_integral + softint)*ao2pi

      rsoft = soft_integral
      risr = risr * 0.5d0*ao2pi
      rfsr = rfsr * 0.5d0*ao2pi
      rifi = rifi * 0.5d0*ao2pi
      return
      end
***********************************************************************************************
      function soft_integral_ee(egmax,p1,p2,p3,p4)
! this should work also for massless momenta...
      implicit double precision (a-h,l,m,o-z)
      parameter (npart = 4)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension p(0:3),q(0:3),pmat(npart,0:3)
      double precision mass(5)
      dimension masses(npart),eta(npart)
      common/sffcommonNEW/ao2pi,masses,eta,alog2,ifirst
      data ifirst/0/
      common/parameters/ame,ammu,convfac,alpha,pi
      double precision phmass
      common/softphotonmass/phmass

* #include "invariants.h"

! ch are the charges of the field (not anti-field!!)
! Use this convention: the integer factor in front of the charge
! must be
! -1 --> for incoming particle
! -1 --> for outgoing anti-particle
! +1 --> for outgoing particle
! +1 --> for incoming anti-particle

      if (ifirst.eq.0) then

         eta(1) = -1.d0 * (1.d0)
         eta(2) = -1.d0 * (-1.d0)
         eta(3) =  -eta(1)
         eta(4) =  -eta(2)

         masses(1) = ame*ame
         masses(2) = masses(1)
         masses(3) = masses(1)
         masses(4) = masses(1)

         ao2pi = alpha/2.d0/pi
         alog2 = log(2.d0)

c         phmass = sqrt(getlambda())

         ifirst = 1
      endif

      linfra    = log(egmax/phmass)
      linfrapl2 = linfra + alog2

      pmat(1,:) = p1
      pmat(2,:) = p2
      pmat(3,:) = p3
      pmat(4,:) = p4

      npartl = npart

! off-diagonal contributions
      softint = 0.d0
      do i = 1,npartl-1
         do j = i+1,npartl
            etaij = eta(i) * eta(j)

            if (abs(etaij).gt.1.d-3) then
               call rescale_momenta_gg(npart,i,j,pmat,masses,p,q,rho)
               q2 = masses(j)
               vl = 0.5d0*(rho*rho*masses(i) - q2)
               v  = vl/(p(0) - q(0))
               arglog = 1.d0+2.d0*vl/q2

               terminfra = 0.d0
               if (arglog.gt.0.d0) then
                  terminfra = linfrapl2*log(arglog)
                  tot = terminfra + termfinite_gg(p,q,v,rho,i,j)
               else
                  tot = 0.d0
               endif

               tot = -tot * 2.d0*dot(p,q)/vl* etaij
               tot =  tot * 2.d0 ! this is the double product when
                                ! squaring the eikonal
               softint  = softint + tot
            endif
         enddo
      enddo
      soft_integral_ee = softint

! diagonal contributions
      softint = 0.d0

c     only boxes
c      goto 666

      do  i = 1,npartl
         etaii = eta(i)*eta(i)
         if (etaii.gt.1.d-3) then
            betai=1.d0-masses(i)/pmat(i,0)/pmat(i,0)
            betai=sqrt(betai)

            if (betai.gt.1.d-12) then
               bterm =  log((1.d0-betai)/(1.d0+betai))/betai
            else
               bterm = -2.d0
            endif
            term = 2.d0*linfrapl2  + bterm
            term = -2.d0 * term * etaii
            softint = softint + term
         endif
      enddo
! divide by 2 to obtain exactly what is inside {} in Dittmaier & Kraemer
! paper (hep-ph/0109062)
 666  soft_integral_ee = 0.5d0*(soft_integral_ee + softint)*ao2pi
      return
      end


**********************************************************
      function soft_integral_muons(egmax,p1,p2,p3,p4)
! this should work also for massless momenta...
      implicit double precision (a-h,l,m,o-z)
      parameter (npart = 4)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension p(0:3),q(0:3),pmat(npart,0:3)
      dimension masses(npart),eta(npart)
      common/sffcommonNEW/ao2pi,masses,eta,alog2,ifirst
      data ifirst/0/
      common/parameters/ame,ammu,convfac,alpha,pi
      double precision phmass
      common/softphotonmass/phmass
      character*2 fs
      common/finalstate/fs

      common/radoption/ifisrfsr
      integer ifisrfsr

* #include "invariants.h"

! ch are the charges of the field (not anti-field!!)
! Use this convention: the integer factor in front of the charge
! must be
! -1 --> for incoming particle
! -1 --> for outgoing anti-particle
! +1 --> for outgoing particle
! +1 --> for incoming anti-particle
      if (ifirst.eq.0) then

         eta(1) = -1.d0 * (1.d0)
         eta(2) = -1.d0 * (-1.d0)
         eta(3) =  -eta(1)
         eta(4) =  -eta(2)

         masses(1) = ame*ame
         masses(2) = masses(1)
         masses(3) = ammu*ammu
         masses(4) = masses(3)

         ao2pi = alpha/2.d0/pi
         alog2 = log(2.d0)

         ifirst = 1
      endif

      linfra    = log(egmax/phmass)
      linfrapl2 = linfra + alog2

      pmat(1,:) = p1
      pmat(2,:) = p2
      pmat(3,:) = p3
      pmat(4,:) = p4

      npartl = npart

      rsoft = 0.d0
      risr  = 0.d0
      rfsr  = 0.d0
      rifi  = 0.d0

! off-diagonal contributions
      softint = 0.d0
      do i = 1,npartl-1
         do j = i+1,npartl
            etaij = eta(i) * eta(j)
            if (abs(etaij).gt.1.d-3) then
               call rescale_momenta_gg(npart,i,j,pmat,masses,p,q,rho)
               q2 = masses(j)
               vl = 0.5d0*(rho*rho*masses(i) - q2)
               v  = vl/(p(0) - q(0))
               arglog = 1.d0+2.d0*vl/q2

               terminfra = 0.d0
               if (arglog.gt.0.d0) then
                  terminfra = linfrapl2*log(arglog)
                  tot = terminfra + termfinite_muons(p,q,v,rho,i,j)
               else
                  tot = 0.d0
               endif

               tot = -tot * 2.d0*dot(p,q)/vl* etaij
               tot =  tot * 2.d0 ! this is the double product when
                                ! squaring the eikonal
               softint  = softint + tot

               if (i.eq.1.and.j.eq.2) then
                  risr = risr + tot
               elseif (i.eq.3.and.j.eq.4) then
                  rfsr = rfsr + tot
               else
                  rifi = rifi + tot
               endif
            endif
         enddo
      enddo
      soft_integral = softint


! diagonal contributions
      softint = 0.d0
      do  i = 1,npartl
         etaii = eta(i)*eta(i)
         if (etaii.gt.1.d-3) then
            betai=1.d0-masses(i)/pmat(i,0)/pmat(i,0)
            betai=sqrt(betai)

            if (betai.gt.1.d-12) then
               bterm =  log((1.d0-betai)/(1.d0+betai))/betai
            else
               bterm = -2.d0
            endif
            term = 2.d0*linfrapl2  + bterm
            term = -2.d0 * term * etaii
            softint = softint + term

            if (i.le.2) then
               risr = risr + term
            else
               rfsr = rfsr + term
            endif
         endif
      enddo


      risr = risr * 0.5d0*ao2pi
      rfsr = rfsr * 0.5d0*ao2pi
      rifi = rifi * 0.5d0*ao2pi
! divide by 2 to obtain exactly what is inside {} in Dittmaier & Kraemer
! paper (hep-ph/0109062)


      if(ifisrfsr.eq.10) then
          soft_integral_muons=risr
          return
      elseif(ifisrfsr.eq.01) then
          soft_integral_muons=rfsr
          return
      elseif(ifisrfsr.eq.11) then
         soft_integral_muons = 0.5d0*(soft_integral + softint)*ao2pi
         return
      endif

      return
      end
**********************************************************************************
      function termfinite_muons(p,q,v,rho,i,j)
      implicit double precision (a-h,m,o-z)
      dimension p(0:3),q(0:3)

      parameter (npart=4)
      dimension masses(npart),eta(npart)
      common/sffcommonNEW/ao2pi,masses,eta,alog2,ifirst
c      common/sffcommonpions/ao2pi,masses,eta,alog2,ifirst

      u0   = p(0)
      umod = sqrt(tridot(p,p))
      pp = u0 + umod
      pm = u0 - umod
      ppopm = pp*pp/rho/rho/masses(i) ! the same as pp/pm

      u0   = q(0)
      umod = sqrt(tridot(q,q))
      qp = u0 + umod
      qm = u0 - umod
      qpoqm = qp*qp/masses(j)  ! the same as qp/qm

      arg1 = (v-pm)/v
      arg2 = (v-pp)/v
      arg3 = (v-qm)/v
      arg4 = (v-qp)/v

      ddlog1 = ddiloghere(arg1)
      ddlog2 = ddiloghere(arg2)
      ddlog3 = ddiloghere(arg3) !! it works, ddlog3 was overwritten below...
      ddlog4 = ddiloghere(arg4)

c      print*,' > ',arg1,arg2,arg3,arg4
c      print*,'ddilog tp0',ddlog1,ddlog2,ddlog3,ddlog4

!     using li2 from LT
c      ddlog1 = li2(arg1)
c      ddlog2 = li2(arg2)
c      ddlog3 = li2(arg3)
c      ddlog4 = li2(arg4)

      termfinite_muons =
     .     log(ppopm)**2/4.d0 + ddlog1
     .                        + ddlog2
     .   - log(qpoqm)**2/4.d0 - ddlog3
     .                        - ddlog4

      return
      end
*************************************************************
      subroutine rescale_momenta_gg(npart,i,j,pmat,masses,p,q,rho)
      implicit double precision (a-h,m,o-z)
      dimension p(0:3),q(0:3),pmat(npart,0:3),masses(npart)
      dimension p1(0:3),p2(0:3)

      p1 = pmat(i,:)
      p2 = pmat(j,:)

      m12 = masses(i)
      m22 = masses(j)

      p1p2 = dot(p1,p2)

      rho1 = p1p2 + sqrt(p1p2**2 - m12*m22)
      rho1 = rho1 / m12
c      rho2 = p1p2 - sqrt(p1p2**2 - m12*m22)
c      rho2 = rho2 / m12
!  better numerical solution !
      rho2 = m22/m12/rho1

      if ( (rho1*p1(0)-p2(0)) .gt. 0.d0) then
         rho = rho1
      endif
      if ( (rho2*p1(0)-p2(0)) .gt. 0.d0) then
         rho = rho2
      endif

      p = rho * p1
      q = p2
      return
      end
***********************************************************************************************
      function termfinite_gg(p,q,v,rho,i,j)
      implicit double precision (a-h,m,o-z)
      dimension p(0:3),q(0:3)

      parameter (npart=4)
      dimension masses(npart),eta(npart)

      common/sffcommonNEW/ao2pi,masses,eta,alog2,ifirst

      common/firsttermfinite/ifirst2
      data ifirst2 /0/

** #include "looptools.h"

      if (ifirst2.eq.0) then
c         call ltini
         ifirst2 = 1
      endif

      u0   = p(0)
      umod = sqrt(tridot(p,p))
      pp = u0 + umod
      pm = u0 - umod
      ppopm = pp*pp/rho/rho/masses(i) ! the same as pp/pm

      u0   = q(0)
      umod = sqrt(tridot(q,q))
      qp = u0 + umod
      qm = u0 - umod
      qpoqm = qp*qp/masses(j)  ! the same as qp/qm

      arg1 = (v-pm)/v
      arg2 = (v-pp)/v
      arg3 = (v-qm)/v
      arg4 = (v-qp)/v

      ddlog1 = ddiloghere(arg1)
      ddlog2 = ddiloghere(arg2)
      ddlog3 = ddiloghere(arg3) !! it works, ddlog3 was overwritten below...
      ddlog4 = ddiloghere(arg4)

c      print*,' > ',arg1,arg2,arg3,arg4
c      print*,'ddilog tp0',ddlog1,ddlog2,ddlog3,ddlog4

!     using li2 from LT
c      ddlog1 = li2(arg1)
c      ddlog2 = li2(arg2)
c      ddlog3 = li2(arg3)
c      ddlog4 = li2(arg4)

      termfinite_gg =
     .     log(ppopm)**2/4.d0 + ddlog1
     .                        + ddlog2
     .   - log(qpoqm)**2/4.d0 - ddlog3
     .                        - ddlog4

      return
      end
***********************************************************************************************
      function termfinite_pions(p,q,v,rho,i,j)
      implicit double precision (a-h,m,o-z)
      dimension p(0:3),q(0:3)

      parameter (npart=4)
      dimension masses(npart),eta(npart)
c      common/sffcommonNEW/ao2pi,masses,eta,alog2,ifirst
      common/sffcommonpions/ao2pi,masses,eta,alog2,ifirst

      u0   = p(0)
      umod = sqrt(tridot(p,p))
      pp = u0 + umod
      pm = u0 - umod
      ppopm = pp*pp/rho/rho/masses(i) ! the same as pp/pm

      u0   = q(0)
      umod = sqrt(tridot(q,q))
      qp = u0 + umod
      qm = u0 - umod
      qpoqm = qp*qp/masses(j)  ! the same as qp/qm

      arg1 = (v-pm)/v
      arg2 = (v-pp)/v
      arg3 = (v-qm)/v
      arg4 = (v-qp)/v

      ddlog1 = ddiloghere(arg1)
      ddlog2 = ddiloghere(arg2)
      ddlog3 = ddiloghere(arg3) !! it works, ddlog3 was overwritten below...
      ddlog4 = ddiloghere(arg4)

c      print*,' > ',arg1,arg2,arg3,arg4
c      print*,'ddilog tp0',ddlog1,ddlog2,ddlog3,ddlog4

!     using li2 from LT
c      ddlog1 = li2(arg1)
c      ddlog2 = li2(arg2)
c      ddlog3 = li2(arg3)
c      ddlog4 = li2(arg4)

      termfinite_pions =
     .     log(ppopm)**2/4.d0 + ddlog1
     .                        + ddlog2
     .   - log(qpoqm)**2/4.d0 - ddlog3
     .                        - ddlog4

      return
      end
*
      subroutine realsoft(p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,mp2,me2,
     &                         rsoft,eps)

#ifdef COLLIER
      use collier
#endif

      IMPLICIT DOUBLE PRECISION (A-H,L-M,O-Z)
      real*8 p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,mp2,me2,rsoft
      complex*16 i11,i22,i33,i44,i12,i13,i14,i23,i24,i34
      real*8 S
      real*8 be1,be2,be3,be4
      real*8 e1,e2,e3,e4
      integer ISR,FSR
      common/radiation/ISR,FSR
      complex*16 soft
      real*8 deltauv,lambda2,mu2dim,m2min,egmin,ebeam
      common/loop/deltauv,mu2dim,lambda2,m2min

      real*8 bornpp
      common/pions/bornpp
      character*2 fs
      common/finalstate/fs

      real*8 esoft

      double precision ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 ecms,ecmsnom
      common/ecms/ecms,ecmsnom

#ifdef LOOPTOOLS
#include "looptools.h"
      call clearcache
      lambda2=getlambda()
#endif

#ifdef COLLIER
      call GetMuIR2_cll(lambda2)
#endif

      me2=ame*ame
      S=2.d0*me2+2.d0*p1p2

      e1 = dsqrt(s)/2.d0
      esoft=eps*dsqrt(s)/2.d0
      e2=e1
      e3=e2
      e4=e3

      be1=sqrt(1.d0-4.d0*me2/S)
      be2=be1
      be3=sqrt(1.d0-4.d0*mp2/S)
      be4=be3

      call intpipi(pi,be1,esoft,lambda2,i11)  !0709.1075  eq (7.20)
      i22=i11    !the beta factor of the electron and the positron is the same
      call intpipi(pi,be3,esoft,lambda2,i33)
      i44=i33


      call intpipjvecp(pi,p1p2,e1,be1,e2,be2,esoft,lambda2,i12) !interference initial-initial (electron)
      call intpipjvecp(pi,p3p4,e3,be3,e4,be4,esoft,lambda2,i34) !interference final-final (pion)

      call intpipj(pi,p1p3,e1,be1,ame,e3,be3,ampi,esoft,lambda2,i13)
      call intpipj(pi,p2p3,e2,be2,ame,e3,be3,ampi,esoft,lambda2,i23)
      call intpipj(pi,p1p4,e1,be1,ame,e4,be4,ampi,esoft,lambda2,i14)
      call intpipj(pi,p2p4,e1,be1,ame,e4,be4,ampi,esoft,lambda2,i24)

      if(fs.eq.'pp') then
	ISR=1
        FSR=1
      else if(fs.eq.'pi') then
	ISR=1
        FSR=0
      endif
      if(ISR.eq.1.and.FSR.eq.0) then
            soft=i11+i22-i12
      else if(ISR.eq.0.and.FSR.gt.0) then
            soft=i33+i44-i34
      else if(ISR.gt.0.and.FSR.gt.0) then
            soft= -i13-i24+i14+i23+i11+i22+i33+i44-i12-i34
      endif
      rsoft= 0.5*(soft+dconjg(soft))
      rsoft= -alpha/2.d0/pi/pi*rsoft

*
      return
      end

      subroutine intpipi(pi,bei,egmin,lambda2,iii)
      real*8 pi
      real*8 bei,egmin,lambda2,riii
      complex*16 iii
      riii= log((4.d0*egmin*egmin)/lambda2)
     +     +1.d0/bei*log((1.d0-bei)/(1.d0+bei))
      riii= pi*riii
      iii= cmplx(riii)

      return
      end
*     subroutine intpipjvecp implements Eq. 7.21 of arXiv:0709.1075,
*     for the integral of soft real radiation when vec(p_i)=-vec(p_j)
*
      subroutine intpipjvecp(pi,pipj,ei,bei,ej,bej,egmin,lambda2,iij)
      implicit none
      real*8 pi,pipj,ei,bei,ej,bej,egmin,lambda2
      complex*16 iij
      real*8 rargmyli2
      complex*16 argmyli2,li2i,li2j
      complex*16 myli2
      external myli2
*
      rargmyli2= 2.d0*bei/(1.d0+bei)
      argmyli2= cmplx(rargmyli2)
      li2i= myli2(argmyli2)
      rargmyli2= 2.d0*bej/(1.d0+bej)
      argmyli2= cmplx(rargmyli2)
      li2j= myli2(argmyli2)
*
      iij= (0.5d0,0.d0)*log((1.d0+bei)/(1.d0-bei))
     +                 *log((4.d0*egmin**2)/lambda2)
     +     -li2i
     +     -(0.25d0,0.d0)*(log((1.d0+bei)/(1.d0-bei)))**2
     +     +(0.5d0,0.d0)*log((1.d0+bej)/(1.d0-bej))
     +                  *log((4.d0*egmin**2)/lambda2)
     +     -li2j
     +     -(0.25d0,0.d0)*(log((1.d0+bej)/(1.d0-bej)))**2
      iij= iij * 2.d0*pi*pipj/(ei+ej)/bei/ei
*
      return
      end
*
*
*     subroutine intpipj implements Eqs. 7.17, 7.18, 7.19 of arXiv:0709.1075,
*     for the integral of soft real radiation, for generic p_i and p_j
*
      subroutine intpipj(pi,pipj,ei,bei,mi,ej,bej,mj,egmin,lambda2,iij)
      implicit none
      real*8 pi,pipj,ei,bei,mi,ej,bej,mj,egmin,lambda2
      complex*16 iij
      real*8 rargmyli2
      complex*16 argmyli2,li2i1,li2i2,li2j1,li2j2
      real*8 alpha,alpha1,alpha2,v,cond1,cond2
      complex*16 myli2
      external myli2
*
      if(abs(mi-mj).lt.1.d-6) then
*        alpha is searched as solution of alpha^2*m_i^2-2*alpha*pipj+m_i^2=0
*        by choosing the one that guarantees alpha > 1
         alpha1= (pipj-sqrt(pipj*pipj-mi*mi*mi*mi))/mi/mi
         alpha2= (pipj+sqrt(pipj*pipj-mi*mi*mi*mi))/mi/mi
         if(alpha1.lt.1.d0.and.alpha2.gt.1.d0) then
            alpha= alpha2
         elseif(alpha1.gt.1.d0.and.alpha2.lt.1.d0) then
            alpha= alpha1
         else
            print*,'problems with alpha in intpipj:'
            print*,'alpha1= ',alpha1
            print*,'alpha2= ',alpha2
            stop
         endif
         v= mi*mi/2.d0/ei*(alpha+1.d0)
      else
*        alpha is searched as solution of alpha^2*m_i^2-2*alpha*pipj+m_j^2=0
*        by choosing the one that guarantees (alpha*ei - ej)/ej > 0
         alpha1= (pipj-sqrt(pipj*pipj-mi*mi*mj*mj))/mi/mi
         alpha2= (pipj+sqrt(pipj*pipj-mi*mi*mj*mj))/mi/mi
         cond1= (alpha1*ei - ej)/ej
         cond2= (alpha2*ei - ej)/ej
         if(cond1.lt.0.d0.and.cond2.gt.0.d0) then
            alpha= alpha2
         elseif(cond1.gt.0.d0.and.cond2.lt.0.d0) then
            alpha= alpha1
         else
            print*,'problems with alpha in intpipj:'
            print*,'alpha1= ',alpha1
            print*,'alpha2= ',alpha2
            stop
         endif
         v= (alpha*alpha*mi*mi-mj*mj)/2.d0/(alpha*ei-ej)
      endif
*
      rargmyli2= 1.d0-alpha*ei*(1.d0+bei)/v
      argmyli2= cmplx(rargmyli2)
      li2i1= myli2(argmyli2)
*
      rargmyli2= 1.d0-alpha*ei*(1.d0-bei)/v
      argmyli2= cmplx(rargmyli2)
      li2i2= myli2(argmyli2)
*
      rargmyli2= 1.d0-ej*(1.d0+bej)/v
      argmyli2= cmplx(rargmyli2)
      li2j1= myli2(argmyli2)
*
      rargmyli2= 1.d0-ej*(1.d0-bej)/v
      argmyli2= cmplx(rargmyli2)
      li2j2= myli2(argmyli2)
*
      iij= (0.5d0,0.d0)*log((alpha*alpha*mi*mi)/(mj*mj))
     +                 *log((4.d0*egmin**2)/lambda2)
     +    +(0.25d0,0.d0)*(log((1.d0-bei)/(1.d0+bei)))**2
     +    +li2i1
     +    +li2i2
     +    -(0.25d0,0.d0)*(log((1.d0-bej)/(1.d0+bej)))**2
     +    -li2j1
     +    -li2j2

      iij= iij * 4.d0*pi*alpha*pipj/(alpha*alpha*mi*mi-mj*mj)
*
      return
      end
***************************************************************
**  DILOGARITHM STUFF
*
** if |z| > 1 calculate li2(1/z) and convert
*
      complex*16  function myli2(zz)
      implicit none
*
      complex*16 cone,zero,ii
      parameter ( cone=(1d0,0d0),
     1            zero=(0d0,0d0),
     2            ii=  (0d0,1d0) )
      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1           pi2=pi*pi)
*
      complex*16 zz
*
      complex*16 myli22
      external myli22
*
c      if(abs(zz-cone).lt.1.d-8) then
      if(abs(zz-cone).lt.1.d-12) then
         myli2=cone*pi2/6.d0
         return
      endif

      if (abs(zz).le.1.d0) then
         myli2=myli22(zz)
      else
         myli2=-log(-zz)**2/2d0-cone*pi2/6d0-myli22(cone/zz)
      endif

      end function myli2
*
** if |z|>1/2 calculate li2(1-z) and convert
*
      complex*16 function myli22(zz)
      implicit none
*
      complex*16 cone,zero,ii
      parameter ( cone=(1d0,0d0),
     1            zero=(0d0,0d0),
     2            ii=  (0d0,1d0) )
      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1           pi2=pi*pi)
*
      complex*16 zz
*
      complex*16 myli23
      external myli23
*
      if (abs(zz).lt.0.5d0) then
         myli22=myli23(zz)
      else
         myli22=cone*pi2/6d0-log(zz)*log(cone-zz)-myli23(cone-zz)
      endif
*
      end function myli22
*
**  main function for dilogarithms
*
      complex*16 function myli23(zz)
      implicit none
*
      complex*16 cone,zero,ii
      parameter ( cone=(1d0,0d0),
     1            zero=(0d0,0d0),
     2            ii=  (0d0,1d0) )
      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1           pi2=pi*pi)
*
      real*8 epsilon
      parameter (epsilon=1.d-20)
      complex*16 zz
*
      integer dim
      parameter (dim = 25)
      real*8 b(dim)
      integer i,n
      complex*16 z,add
*
* even n Bernoulli numbers, already divided by (n + 1)!
*
      data b /
     $  2.777777777777777777777777777777777777777777777777777777778d-2,
     $ -2.777777777777777777777777777777777777777777777777777777778d-4,
     $  4.724111866969009826152683295540438397581254724111866969009d-6,
     $ -9.185773074661963550852439741328630217519106407995296884185d-8,
     $  1.897886998897099907200917301927402937503947604957705967806d-9,
     $ -4.064761645144225526805909386291966674547057127439707822288d-11,
     $  8.921691020456452555217987316752748851514283613049045147810d-13,
     $ -1.993929586072107568723644347793789705630694749653880147036d-14,
     $  4.518980029619918191650476552855593228396819014466618405199d-16,
     $ -1.035651761218124701448341154221865666596091238168650515964d-17,
     $  2.395218621026186745740283743000980381678949001942974256251d-19,
     $ -5.581785874325009336283074505625419905567054667644398095136d-21,
     $  1.309150755418321285812307399186592301749849838783303836854d-22,
     $ -3.087419802426740293242279764866462431595565256132745695326d-24,
     $  7.315975652702203420357905609252148591033401063690875035693d-26,
     $ -1.740845657234000740989055147759702545340841421754271264171d-27,
     $  4.157635644613899719617899620775226673488254159511563860825d-29,
     $ -9.962148488284622103194006702455838849854860017394488768062d-31,
     $  2.394034424896165300521167987893749562934279156932915750221d-32,
     $ -5.768347355367390084291793161877654244072332317926275110062d-34,
     $  1.393179479647007977827886603911548331732411625673399565806d-35,
     $ -3.372121965485089470468473635254930958979742891656539304386d-37,
     $  8.178208777562102621764777214872834267876189462495503276198d-39,
     $ -1.987010831152385925564820669234786567541858995824743201790d-40,
     $  4.835778518040550896287059373115378207694465369420827842732d-42
     $  /
*
      z=-log(cone-zz)

      myli23=z-z**2/4d0

      if (abs(myli23).gt.epsilon) then
          do i=1,dim
             n    = 2*i
             add  = z**(n+1)*b(i)
             if (abs(add/myli23).gt.epsilon) then
                 myli23 = myli23 + add
             else
                 exit
             endif
          enddo
      endif
      end function myli23

#ifdef COLLIER
*************FROM MESMER
***************************************
      function soft_integral_mue(egmax,p1,p2,p3,p4)
! this should work also for massless momenta...
      implicit double precision (a-h,l,m,o-z)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      parameter (npart = 4)
      dimension p(0:3),q(0:3),pmat(npart,0:3)
      double precision mass(5)
      dimension masses(npart),eta(npart)
      common/sffcommonmue/ao2pi,masses,eta,alog2,ifirst
      common/pmatrix/pmat
      common/parameters/ame,ammu,convfac,alpha,pi
      data ifirst/0/

      character*2 fs
      common/finalstate/fs

      double precision phmass
      common/softphotonmass/phmass

      integer maskud(4,4),masknoud(4,4)
      common/sffmasks/maskud,masknoud
      common/sffudandnoud/sffud,sffnoud

#include "invariants.h"

! ch are the charges of the field (not anti-field!!)
! Use this convention: the integer factor in front of the charge
! must be
! -1 --> for incoming particle
! -1 --> for outgoing anti-particle
! +1 --> for outgoing particle
! +1 --> for incoming anti-particle

      if (ifirst.eq.0) then

c     eta(1) =  1.d0 * (-1.d0)
c     eta(2) = -1.d0 * (-1.d0)
c     eta(3) = -1.d0 * (-1.d0)
c     eta(4) =  1.d0 * (-1.d0)

         eta(1) = -1.d0 * (-1.d0) * QRe ! this works to cancel eps against hard
         eta(2) = +1.d0 * (-1.d0) * QRe
         eta(3) = +1.d0 * (-1.d0) * QRmu
         eta(4) = -1.d0 * (-1.d0) * QRmu
         masses(1) = ame*ame
         masses(2) = ame*ame
         masses(3) = ammu*ammu
         masses(4) = ammu*ammu

c     eta(1) = -1.d0 * Qmu     * QRmu
c     eta(2) = -1.d0 * (-1.d0) * QRe
c     eta(3) =  1.d0 * Qmu     * QRmu
c     eta(4) =  1.d0 * (-1.d0) * QRe
c     masses(1) = mm*mm
c     masses(2) = me*me
c     masses(3) = mm*mm
c     masses(4) = me*me


         ao2pi = alpha/2.d0/pi
         alog2 = log(2.d0)

***************************************
         maskud   = 0
         masknoud = 0

         maskud(1,2) = 1
         maskud(1,4) = 1
         maskud(2,3) = 1
         maskud(3,4) = 1
         maskud(2,1) = 1
         maskud(4,1) = 1
         maskud(3,2) = 1
         maskud(4,3) = 1

         masknoud(1,1) = 1
         masknoud(2,2) = 1
         masknoud(3,3) = 1
         masknoud(4,4) = 1

         masknoud(1,3) = 1
         masknoud(3,1) = 1
         masknoud(2,4) = 1
         masknoud(4,2) = 1
***************************************
         ifirst = 1
      endif

      linfra    = log(egmax/phmass)
      linfrapl2 = linfra + alog2

      pmat(1,:) = p1
      pmat(2,:) = p2
      pmat(3,:) = p3
      pmat(4,:) = p4

      npartl = npart

! off-diagonal contributions
      softint = 0.d0

      do i = 1,npartl-1
         do j = i+1,npartl
            etaij = eta(i) * eta(j)
            if (abs(etaij).gt.1.d-3) then
               call rescale_momenta_mue(npart,i,j,pmat,masses,p,q,rho)
               q2 = masses(j)
               vl = 0.5d0*(rho*rho*masses(i) - q2)
               v  = vl/(p(0) - q(0))
c               arglog = 1.d0+2.d0*vl/q2
* i.e.
               arglog=rho*rho*masses(i)/masses(j)
c               if (arglog.gt.0.d0) then
                 tot=linfrapl2*log(arglog)+termfinite_mue(p,q,v,rho,i,j)
c               else
c                  tot = 0.d0
c               endif
               tot = -tot * 2.d0*dot(p,q)/vl* etaij
               tot =  tot * 2.d0 ! this is the double product when squaring the eikonal
               softint   = softint  + tot

c               sffij(i,j) = tot
c               sffij(j,i) = tot
            endif
         enddo
      enddo
      soft_integral_mue = softint

! diagonal contributions
      softint = 0.d0
      do  i = 1,npartl
         etaii = eta(i)*eta(i)
         if (etaii.gt.1.d-3) then
c            betai=1.d0-masses(i)/pmat(i,0)/pmat(i,0)
c            betai=sqrt(betai)
            betai = sqrt(tridot(pmat(i,:),pmat(i,:)))/pmat(i,0)
            if (betai.gt.1.d-12) then
               bterm =  log((1.d0-betai)/(1.d0+betai))/betai
            else
               bterm = -2.d0
            endif
            term = 2.d0*linfrapl2  + bterm
            term = -2.d0 * term * etaii
            softint = softint + term
cc              print*, 2.d0*linfrapl2
****
c            sffij(i,i) = term
         endif
      enddo
! divide by 2 to obtain exactly what is inside {} in Dittmaier & Kraemer
!     paper (hep-ph/0109062)

      soft_integral_mue = 0.5d0*(soft_integral_mue + softint)*ao2pi

c      sffij = 0.5d0*ao2pi * sffij
c     c      sffdebugIR = 0.5d0*ao2pi*sffdebugIR
      return
      end
***********************************************************************************************
      subroutine rescale_momenta_mue(npart,i,j,pmat,masses,p,q,rho)
      implicit double precision (a-h,m,o-z)
      dimension p(0:3),q(0:3),pmat(npart,0:3),masses(npart)
      dimension p1(0:3),p2(0:3)

      p1  = pmat(i,:)
      p2  = pmat(j,:)
      m12 = masses(i)
      m22 = masses(j)

      p1p2 = dot(p1,p2)

      rho1 = p1p2 + sqrt(p1p2*p1p2 - m12*m22)
      rho1 = rho1 / m12
c      rho2 = p1p2 - sqrt(p1p2**2 - m12*m22)
c      rho2 = rho2 / m12
!  better numerical solution !
      rho2 = m22/m12/rho1
      if ( (rho1*p1(0)-p2(0)) .gt. 0.d0) then
         rho = rho1
      endif
      if ( (rho2*p1(0)-p2(0)) .gt. 0.d0) then
         rho = rho2
      endif
      p = rho * p1
      q = p2
      return
      end
***********************************************************************************************
      function termfinite_mue(p,q,v,rho,i,j)
      implicit double precision (a-h,m,o-z)
      dimension p(0:3),q(0:3)
      double complex cspen,carg
      parameter (npart=4)
      dimension masses(npart),eta(npart)
      double precision phmass
      common/photonmasslambda/phmass
      common/sffcommonmue/ao2pi,masses,eta,alog2,ifirst
      double complex myli2
      external myli2

      u0   = p(0)
      umod = sqrt(tridot(p,p))
      pp = u0 + umod
      pm = u0 - umod
      ppopm = pp*pp/rho/rho/masses(i) ! the same as pp/pm

      u0   = q(0)
      umod = sqrt(tridot(q,q))
      qp = u0 + umod
      qm = u0 - umod
      qpoqm = qp*qp/masses(j)  ! the same as qp/qm

c      arg1 = (v-pm)/v
c      arg2 = (v-pp)/v
c      arg3 = (v-qm)/v
c      arg4 = (v-qp)/v

      arg1 = 1.d0-pm/v
      arg2 = 1.d0-pp/v
      arg3 = 1.d0-qm/v
      arg4 = 1.d0-qp/v

      ddlog1 = dilogmy(arg1)
      ddlog2 = dilogmy(arg2)
      ddlog3 = dilogmy(arg3) !well it works...!
      ddlog4 = dilogmy(arg4)

c      ddlog1 = myli2(dcmplx(arg1))
c      ddlog2 = myli2(dcmplx(arg2))
c      ddlog3 = myli2(dcmplx(arg3)) !well it works...!
c      ddlog4 = myli2(dcmplx(arg4))

      termfinite_mue =
     .     0.25d0*log(ppopm)**2 + ddlog1 + ddlog2
     .   - 0.25d0*log(qpoqm)**2 - ddlog3 - ddlog4
      return
      end
*****************************************************************
*****************************************************************
      double precision FUNCTION DILOGMY(X)
      double precision X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
      double precision C(0:18),H,ALFA,B0,B1,B2
      DATA ZERO /0.0D0/, ONE /1.0D0/
      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/
      DATA PI3 /3.289868133696453D0/, PI6 /1.644934066848226D0/
      DATA C( 0) / 0.4299669356081370D0/
      DATA C( 1) / 0.4097598753307711D0/
      DATA C( 2) /-0.0185884366501460D0/
      DATA C( 3) / 0.0014575108406227D0/
      DATA C( 4) /-0.0001430418444234D0/
      DATA C( 5) / 0.0000158841554188D0/
      DATA C( 6) /-0.0000019078495939D0/
      DATA C( 7) / 0.0000002419518085D0/
      DATA C( 8) /-0.0000000319334127D0/
      DATA C( 9) / 0.0000000043454506D0/
      DATA C(10) /-0.0000000006057848D0/
      DATA C(11) / 0.0000000000861210D0/
      DATA C(12) /-0.0000000000124433D0/
      DATA C(13) / 0.0000000000018226D0/
      DATA C(14) /-0.0000000000002701D0/
      DATA C(15) / 0.0000000000000404D0/
      DATA C(16) /-0.0000000000000061D0/
      DATA C(17) / 0.0000000000000009D0/
      DATA C(18) /-0.0000000000000001D0/
      IF(X .EQ. ONE) THEN
         DILOGMY=PI6
         RETURN
      ELSE IF(X .EQ. MONE) THEN
         DILOGMY=MALF*PI6
         RETURN
      END IF
      T=-X
      IF(T .LE. MTWO) THEN
         Y=MONE/(ONE+T)
         S=ONE
         A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)
      ELSE IF(T .LT. MONE) THEN
         Y=MONE-T
         S=MONE
         A=LOG(-T)
         A=-PI6+A*(A+LOG(ONE+ONE/T))
      ELSE IF(T .LE. MALF) THEN
         Y=(MONE-T)/T
         S=ONE
         A=LOG(-T)
         A=-PI6+A*(MALF*A+LOG(ONE+T))
      ELSE IF(T .LT. ZERO) THEN
         Y=-T/(ONE+T)
         S=MONE
         A=HALF*LOG(ONE+T)**2
      ELSE IF(T .LE. ONE) THEN
         Y=T
         S=ONE
         A=ZERO
      ELSE
         Y=ONE/T
         S=MONE
         A=PI6+HALF*LOG(T)**2
      END IF
      H=Y+Y-ONE
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO  I = 18,0,-1
         B0=C(I)+ALFA*B1-B2
         B2=B1
 1       B1=B0
      enddo
      DILOGMY=-(S*(B0-H*B2)+A)
c      print*,dilogmy, li2(x)
      RETURN
      END
***********************************************************************************************
***************************************************************************
      subroutine loadinvariants(p1,p2,p3,p4)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      common/invariants/p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,tm1
      double precision dot
      external dot
      p1p2 = dot(p1,p2)
      p1p3 = dot(p1,p3)
      p1p4 = dot(p1,p4)
      p2p3 = dot(p2,p3)
      p2p4 = dot(p2,p4)
      p3p4 = dot(p3,p4)
      return
      end
***************************************************************************
      subroutine getQEDcts(cts)
      implicit none
      double precision cts,ctstmp
      double complex tmp
      double complex B1lt,dB0lt,dB1lt
!#include "invariants.h"
      double precision ec,ec2,ec4,mm,mm2,me,me2,mm4,me4,mm6,me6
      common/elchmasses/ec,ec2,ec4,mm,mm2,me,me2,mm4,me4,mm6,me6
      double precision pi216m1,pi216
      common/pi216m1/pi216m1,pi216
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu
      double precision Qmu
      common/muoncharge/Qmu

      integer ifirst
      common/qedcts/ctstmp,ifirst
      data ifirst /0/

      double precision ctstmpmu,ctstmpe
      common/qedcts_e_mu/ctstmpmu,ctstmpe

      double complex Sigrlii
      double complex dSigrlii
      double complex SigSii
      double complex dSigSii

      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi

      if (ifirst.eq.0) then
c         i     = Bget(mm2,mm2,0.d0)
c         B1lt  = Bval(bb1,i)
c        dB0lt = Bval(dbb0,i)
c         dB1lt = Bval(dbb1,i)
         call getBsforctscollier(mm2,B1lt,dB0lt,dB1lt)
         tmp   = B1lt + 0.5d0 + 2.d0*mm2*dB1lt + 4.d0*mm2*dB0lt

         cts = 2.d0*tmp *pi216m1 * ec2 * QRmu
         ctstmpmu = cts *QRmu

c         i     = Bget(me2,me2,0.d0)
c         B1lt  = Bval(bb1,i)
c         dB0lt = Bval(dbb0,i)
c         dB1lt = Bval(dbb1,i)

         call getBsforctscollier(me2,B1lt,dB0lt,dB1lt)

         tmp   = B1lt + 0.5d0 + 2.d0*me2*dB1lt + 4.d0*me2*dB0lt
         cts = cts + 2.d0*tmp *pi216m1 * ec2 * QRe
         ctstmpe = 2.d0*tmp *pi216m1 * ec2 * QRe

         ctstmp   = 2.d0 * cts
         ctstmpmu = 2.d0*ctstmpmu
         ctstmpe  = 2.d0*ctstmpe

         ifirst = 1 !!!!!!!!!!!!!!!!!!!!!!
      endif
      cts = ctstmp*1.d0
      return
      end
*****************
*****************************************************
      subroutine getBsforctscollier(m2,B1o,dB0o,dB1o)
      use collier
      implicit none
      double complex B1o,dB0o,dB1o
      double precision m2
      integer rank
      parameter (rank=1)
      double complex b(0:rank/2,0:rank),buv(0:rank/2,0:rank)

      call SwitchOffCacheSystem_cll
      call B_cll(b,buv,(1.d0,0.d0)*m2,(1.d0,0.d0)*m2,(0.d0,0.d0),1)
      B1o = b(0,1)
      call DB0_cll(dB0o,(1.d0,0.d0)*m2,(1.d0,0.d0)*m2,(0.d0,0.d0))
      call DB1_cll(dB1o,(1.d0,0.d0)*m2,(1.d0,0.d0)*m2,(0.d0,0.d0))
      call SwitchOnCacheSystem_cll
      return
      end
**************************************************************************
      subroutine Dfunscollier
      use collier
      double complex D0ltd,D1ltd,D2ltd,D3ltd
      double complex D0ltc,D1ltc,D2ltc,D3ltc
      double complex D00ltd,D11ltd,D12ltd,D13ltd
      double complex D22ltd,D23ltd,D33ltd
      double complex D00ltc,D11ltc,D12ltc,D13ltc
      double complex D22ltc,D23ltc,D33ltc
      common/boxesfuns/
     .     D0ltd,D1ltd,D2ltd,D3ltd,
     .     D0ltc,D1ltc,D2ltc,D3ltc,
     .     D00ltd,D11ltd,D12ltd,D13ltd,
     .     D22ltd,D23ltd,D33ltd,
     .     D00ltc,D11ltc,D12ltc,D13ltc,
     .     D22ltc,D23ltc,D33ltc

      integer rank ! max rank
      parameter (rank = 2)
      double complex ptmp(0:3)
      double complex dnd(0:rank/2,0:rank,0:rank,0:rank)
      double complex dnduv(0:rank/2,0:rank,0:rank,0:rank)
      double complex dnc(0:rank/2,0:rank,0:rank,0:rank)
      double complex dncuv(0:rank/2,0:rank,0:rank,0:rank)

      double complex Dtend(0:rank,0:rank,0:rank,0:rank)
      double complex Dtenc(0:rank,0:rank,0:rank,0:rank)
      double complex Dtenda(0:rank,0:rank,0:rank,0:rank)
      double complex Dtenca(0:rank,0:rank,0:rank,0:rank)
      double complex Dtenduv(0:rank,0:rank,0:rank,0:rank)
      double complex Dtencuv(0:rank,0:rank,0:rank,0:rank)
      common/boxestenfuns/Dtend,Dtenc
      data Dtend /81*(0.d0,0.d0)/
      data Dtenc /81*(0.d0,0.d0)/

      double precision k1(0:3),k2(0:3),k3(0:3)
      double complex dotten
      external dotten
      double complex dottencmom
      external dottencmom

#include "invariants.h"

      double complex d0cll,cpard(10),cparc(10),cu
      common/dcpar/cpard,cparc,cu,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         cpard(1)  = cuno*mm2
         cpard(2)  = cuno*mm2
         cpard(3)  = cuno*me2
         cpard(4)  = cuno*me2
         cpard(7)  = cuno*0.d0
         cpard(8)  = cuno*mm2
         cpard(9)  = cuno*0.d0
         cpard(10) = cuno*me2
         cparc = cpard
         cu = (1.d0,0.d0)
         ifirst   = 1
      endif
      cpard(5) = cuno*(2.d0*mm2-2.d0*p1p3)
      cpard(6) = cuno*(mm2+me2+2.d0*p3p4)
      cparc(5) = cuno*(2.d0*mm2-2.d0*p1p3)
      cparc(6) = cuno*(mm2+me2-2.d0*p2p3)

! direct box
      call D_cll(dnd,dnduv,
     .     cpard(1),
     .     cpard(2),
     .     cpard(3),
     .     cpard(4),
     .     cpard(5),
     .     cpard(6),
     .     cpard(7),
     .     cpard(8),
     .     cpard(9),
     .     cpard(10),rank)

! crossed box
      call D_cll(dnc,dncuv,
     .     cparc(1),
     .     cparc(2),
     .     cparc(3),
     .     cparc(4),
     .     cparc(5),
     .     cparc(6),
     .     cparc(7),
     .     cparc(8),
     .     cparc(9),
     .     cparc(10),rank)

      D0ltd  = dnd(0,0,0,0)
      D1ltd  = dnd(0,1,0,0)
      D2ltd  = dnd(0,0,1,0)
      D3ltd  = dnd(0,0,0,1)
      D00ltd = dnd(1,0,0,0)
      D11ltd = dnd(0,2,0,0)
      D12ltd = dnd(0,1,1,0)
      D13ltd = dnd(0,1,0,1)
      D22ltd = dnd(0,0,2,0)
      D23ltd = dnd(0,0,1,1)
      D33ltd = dnd(0,0,0,2)

      D0ltc  = dnc(0,0,0,0)
      D1ltc  = dnc(0,1,0,0)
      D2ltc  = dnc(0,0,1,0)
      D3ltc  = dnc(0,0,0,1)
      D00ltc = dnc(1,0,0,0)
      D11ltc = dnc(0,2,0,0)
      D12ltc = dnc(0,1,1,0)
      D13ltc = dnc(0,1,0,1)
      D22ltc = dnc(0,0,2,0)
      D23ltc = dnc(0,0,1,1)
      D33ltc = dnc(0,0,0,2)

      return
      end
**************************************************************************
      subroutine Cfunscollier
      use collier
      double complex C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu
      double complex C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte
      double complex C00ltmu,C00lte
#include "invariants.h"

      common/vertexfuns/
     .     C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu,
     .     C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte,
     .     C00ltmu,C00lte
      data
     .     C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu,
     .     C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte,
     .     C00ltmu,C00lte
     .     /14*(0.d0,0.d0)/

      integer rank ! max rank
      parameter (rank = 2)
      double complex cnmu(0:rank/2,0:rank,0:rank)
      double complex cnmuuv(0:rank/2,0:rank,0:rank)
      double complex cne(0:rank/2,0:rank,0:rank)
      double complex cneuv(0:rank/2,0:rank,0:rank)

      double complex cparmu(6),cpare(6)
      common/ccpar/cparmu,cpare,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         cparmu(1)  = cuno*mm2
         cparmu(3)  = cuno*mm2
         cparmu(4)  = cuno*0.d0
         cparmu(5)  = cuno*mm2
         cparmu(6)  = cuno*mm2

         cpare(1)  = cuno*me2
         cpare(3)  = cuno*me2
         cpare(4)  = cuno*0.d0
         cpare(5)  = cuno*me2
         cpare(6)  = cuno*me2
         ifirst    = 1
      endif
      cparmu(2) = cuno*(2.d0*mm2-2.d0*p1p3)
      cpare(2)  = cuno*(2.d0*me2-2.d0*p2p4)

      if (QRmu.gt.0) then
! muon
         call C_cll(cnmu,cnmuuv,
     .        cparmu(1),
     .        cparmu(2),
     .        cparmu(3),
     .        cparmu(4),
     .        cparmu(5),
     .        cparmu(6),
     .        rank)
         C0ltmu  = cnmu(0,0,0)
         C1ltmu  = cnmu(0,1,0)
         C00ltmu = cnmu(1,0,0)
         C11ltmu = cnmu(0,2,0)
         C12ltmu = cnmu(0,1,1)
      endif
! electron
      if (QRe.gt.0) then
         call C_cll(cne,cneuv,
     .        cpare(1),
     .        cpare(2),
     .        cpare(3),
     .        cpare(4),
     .        cpare(5),
     .        cpare(6),
     .        rank)
         C0lte  = cne(0,0,0)
         C1lte  = cne(0,1,0)
         C00lte = cne(1,0,0)
         C11lte = cne(0,2,0)
         C12lte = cne(0,1,1)
      endif
      return
      end

*****************
      subroutine loadQEDvertexfuns
      implicit none
      double complex C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu
      double complex C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte
      double complex C00ltmu,C00lte
!#include "invariants.h"
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu
      common/vertexfuns/
     .     C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu,
     .     C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte,
     .     C00ltmu,C00lte

      data
     .     C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu,
     .     C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte,
     .     C00ltmu,C00lte
     .     /14*(0.d0,0.d0)/


      call Cfunscollier

      return
      end
*********************************************************
********************************************************
      subroutine loadQEDboxesfuns
      implicit none
      double complex D0ltd,D1ltd,D2ltd,D3ltd
      double complex D0ltc,D1ltc,D2ltc,D3ltc
      double complex D00ltd,D11ltd,D12ltd,D13ltd
      double complex D22ltd,D23ltd,D33ltd
      double complex D00ltc,D11ltc,D12ltc,D13ltc
      double complex D22ltc,D23ltc,D33ltc
      double precision pard(10),parc(10)
!#include "invariants.h"
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu
      common/boxesfuns/
     .     D0ltd,D1ltd,D2ltd,D3ltd,
     .     D0ltc,D1ltc,D2ltc,D3ltc,
     .     D00ltd,D11ltd,D12ltd,D13ltd,
     .     D22ltd,D23ltd,D33ltd,
     .     D00ltc,D11ltc,D12ltc,D13ltc,
     .     D22ltc,D23ltc,D33ltc
      integer ifirst
      double complex cu
      common/loadboxes/pard,parc,cu,ifirst
      data ifirst /0/

      data
     .     D0ltd,D1ltd,D2ltd,D3ltd,
     .     D0ltc,D1ltc,D2ltc,D3ltc,
     .     D00ltd,D11ltd,D12ltd,D13ltd,
     .     D22ltd,D23ltd,D33ltd,
     .     D00ltc,D11ltc,D12ltc,D13ltc,
     .     D22ltc,D23ltc,D33ltc
     .     /22*(0.d0,0.d0)/

      if (QRmu.eq.0.or.QRe.eq.0) return

      call Dfunscollier
      return
      end
*********************************************
      double complex function dotcmom(p,q)
      implicit none
      double complex p(0:3),q(0:3)
      dotcmom=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
      return
      end
********************************************
***************************************************
      subroutine getcurrentnew(ch,h1,h2,p1,p2,m1,m2,c,cI,c5)
      implicit none
      integer h1,h2
      double precision p1(0:3),p2(0:3),den12,den22,m1,m2,ch
      double complex c0,c1,c2,c3
      double complex c50,c51,c52,c53
      double complex nm1
      double complex c(0:3),im,cI,c5(0:3)
      parameter (im = (0.d0,1.d0))

      den12 = p1(0) + m1
      den22 = p2(0) + m2

      if (ch.lt.0.d0) then

      if (h1.eq.1.and.h2.eq.1) then
      c0 =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      c0 = c0 + den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c1 =
     &  + im * (  - p1(2)*den22 + p2(2)*den12 )
      c1 = c1 + p1(1)*den22 + p2(1)*den12

      c2 =
     &  + im * ( p1(1)*den22 - p2(1)*den12 )
      c2 = c2 + p1(2)*den22 + p2(2)*den12

      c3 =
     &  + p1(3)*den22 + p2(3)*den12

      c50 =
     &  - p1(3)*den22 - p2(3)*den12

      c51 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c51 = c51 - p1(1)*p2(3) - p1(3)*p2(1)

      c52 =
     &  + im * (  - p1(1)*p2(3) + p1(3)*p2(1) )
      c52 = c52 - p1(2)*p2(3) - p1(3)*p2(2)

      c53 =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      c53 = c53 - den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) - p1(3)*p2(3)

      cI =
     &  + im * (  - p1(1)*p2(2) + p1(2)*p2(1) )
      cI = cI + den12*den22 - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)
      elseif (h1.eq.1.and.h2.eq.2) then
      c0 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c0 = c0 - p1(1)*p2(3) + p1(3)*p2(1)

      c1 =
     &  + p1(3)*den22 - p2(3)*den12

      c2 =
     &  + im * (  - p1(3)*den22 + p2(3)*den12 )

      c3 =
     &  + im * ( p1(2)*den22 - p2(2)*den12 )
      c3 = c3 - p1(1)*den22 + p2(1)*den12

      c50 =
     &  + im * ( p1(2)*den22 + p2(2)*den12 )
      c50 = c50 - p1(1)*den22 - p2(1)*den12

      c51 =
     &  + im * ( p1(1)*p2(2) + p1(2)*p2(1) )
      c51 = c51 - den12*den22 - p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c52 =
     &  + im * ( den12*den22 - p1(1)*p2(1) + p1(2)*p2(2) - p1(3)*p2(3)
     &     )
      c52 = c52 - p1(1)*p2(2) - p1(2)*p2(1)

      c53 =
     &  + im * ( p1(2)*p2(3) + p1(3)*p2(2) )
      c53 = c53 - p1(1)*p2(3) - p1(3)*p2(1)

      cI =
     &  + im * (  - p1(2)*p2(3) + p1(3)*p2(2) )
      cI = cI + p1(1)*p2(3) - p1(3)*p2(1)
      elseif (h1.eq.2.and.h2.eq.1) then
      c0 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c0 = c0 + p1(1)*p2(3) - p1(3)*p2(1)

      c1 =
     &  - p1(3)*den22 + p2(3)*den12

      c2 =
     &  + im * (  - p1(3)*den22 + p2(3)*den12 )

      c3 =
     &  + im * ( p1(2)*den22 - p2(2)*den12 )
      c3 = c3 + p1(1)*den22 - p2(1)*den12

      c50 =
     &  + im * (  - p1(2)*den22 - p2(2)*den12 )
      c50 = c50 - p1(1)*den22 - p2(1)*den12

      c51 =
     &  + im * (  - p1(1)*p2(2) - p1(2)*p2(1) )
      c51 = c51 - den12*den22 - p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c52 =
     &  + im * (  - den12*den22 + p1(1)*p2(1) - p1(2)*p2(2) + p1(3)*
     &    p2(3) )
      c52 = c52 - p1(1)*p2(2) - p1(2)*p2(1)

      c53 =
     &  + im * (  - p1(2)*p2(3) - p1(3)*p2(2) )
      c53 = c53 - p1(1)*p2(3) - p1(3)*p2(1)

      cI =
     &  + im * (  - p1(2)*p2(3) + p1(3)*p2(2) )
      cI = cI - p1(1)*p2(3) + p1(3)*p2(1)
      elseif (h1.eq.2.and.h2.eq.2) then
      c0 =
     &  + im * (  - p1(1)*p2(2) + p1(2)*p2(1) )
      c0 = c0 + den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c1 =
     &  + im * ( p1(2)*den22 - p2(2)*den12 )
      c1 = c1 + p1(1)*den22 + p2(1)*den12

      c2 =
     &  + im * (  - p1(1)*den22 + p2(1)*den12 )
      c2 = c2 + p1(2)*den22 + p2(2)*den12

      c3 =
     &  + p1(3)*den22 + p2(3)*den12

      c50 =
     &  + p1(3)*den22 + p2(3)*den12

      c51 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c51 = c51 + p1(1)*p2(3) + p1(3)*p2(1)

      c52 =
     &  + im * (  - p1(1)*p2(3) + p1(3)*p2(1) )
      c52 = c52 + p1(2)*p2(3) + p1(3)*p2(2)

      c53 =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      c53 = c53 + den12*den22 - p1(1)*p2(1) - p1(2)*p2(2) + p1(3)*p2(3)

      cI =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      cI = cI + den12*den22 - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)

      endif

      else ! if ch > 0

      if (h1.eq.1.and.h2.eq.1) then
      c0 =
     &  + im * (  - p1(1)*p2(2) + p1(2)*p2(1) )
      c0 = c0 + den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c1 =
     &  + im * ( p1(2)*den22 - p2(2)*den12 )
      c1 = c1 + p1(1)*den22 + p2(1)*den12

      c2 =
     &  + im * (  - p1(1)*den22 + p2(1)*den12 )
      c2 = c2 + p1(2)*den22 + p2(2)*den12

      c3 =
     &  + p1(3)*den22 + p2(3)*den12

      c50 =
     &  + p1(3)*den22 + p2(3)*den12

      c51 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c51 = c51 + p1(1)*p2(3) + p1(3)*p2(1)

      c52 =
     &  + im * (  - p1(1)*p2(3) + p1(3)*p2(1) )
      c52 = c52 + p1(2)*p2(3) + p1(3)*p2(2)

      c53 =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      c53 = c53 + den12*den22 - p1(1)*p2(1) - p1(2)*p2(2) + p1(3)*p2(3)

      cI =
     &  + im * (  - p1(1)*p2(2) + p1(2)*p2(1) )
      cI = cI - den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      elseif (h1.eq.1.and.h2.eq.2) then
      c0 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c0 = c0 + p1(1)*p2(3) - p1(3)*p2(1)

      c1 =
     &  - p1(3)*den22 + p2(3)*den12

      c2 =
     &  + im * (  - p1(3)*den22 + p2(3)*den12 )

      c3 =
     &  + im * ( p1(2)*den22 - p2(2)*den12 )
      c3 = c3 + p1(1)*den22 - p2(1)*den12

      c50 =
     &  + im * (  - p1(2)*den22 - p2(2)*den12 )
      c50 = c50 - p1(1)*den22 - p2(1)*den12

      c51 =
     &  + im * (  - p1(1)*p2(2) - p1(2)*p2(1) )
      c51 = c51 - den12*den22 - p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c52 =
     &  + im * (  - den12*den22 + p1(1)*p2(1) - p1(2)*p2(2) + p1(3)*
     &    p2(3) )
      c52 = c52 - p1(1)*p2(2) - p1(2)*p2(1)

      c53 =
     &  + im * (  - p1(2)*p2(3) - p1(3)*p2(2) )
      c53 = c53 - p1(1)*p2(3) - p1(3)*p2(1)

      cI =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      cI = cI + p1(1)*p2(3) - p1(3)*p2(1)

      elseif (h1.eq.2.and.h2.eq.1) then
      c0 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c0 = c0 - p1(1)*p2(3) + p1(3)*p2(1)

      c1 =
     &  + p1(3)*den22 - p2(3)*den12

      c2 =
     &  + im * (  - p1(3)*den22 + p2(3)*den12 )

      c3 =
     &  + im * ( p1(2)*den22 - p2(2)*den12 )
      c3 = c3 - p1(1)*den22 + p2(1)*den12

      c50 =
     &  + im * ( p1(2)*den22 + p2(2)*den12 )
      c50 = c50 - p1(1)*den22 - p2(1)*den12

      c51 =
     &  + im * ( p1(1)*p2(2) + p1(2)*p2(1) )
      c51 = c51 - den12*den22 - p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c52 =
     &  + im * ( den12*den22 - p1(1)*p2(1) + p1(2)*p2(2) - p1(3)*p2(3)
     &     )
      c52 = c52 - p1(1)*p2(2) - p1(2)*p2(1)

      c53 =
     &  + im * ( p1(2)*p2(3) + p1(3)*p2(2) )
      c53 = c53 - p1(1)*p2(3) - p1(3)*p2(1)

      cI =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      cI = cI - p1(1)*p2(3) + p1(3)*p2(1)

      elseif (h1.eq.2.and.h2.eq.2) then
      c0 =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      c0 = c0 + den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c1 =
     &  + im * (  - p1(2)*den22 + p2(2)*den12 )
      c1 = c1 + p1(1)*den22 + p2(1)*den12

      c2 =
     &  + im * ( p1(1)*den22 - p2(1)*den12 )
      c2 = c2 + p1(2)*den22 + p2(2)*den12

      c3 =
     &  + p1(3)*den22 + p2(3)*den12

      c50 =
     &  - p1(3)*den22 - p2(3)*den12

      c51 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c51 = c51 - p1(1)*p2(3) - p1(3)*p2(1)

      c52 =
     &  + im * (  - p1(1)*p2(3) + p1(3)*p2(1) )
      c52 = c52 - p1(2)*p2(3) - p1(3)*p2(2)

      c53 =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      c53 = c53 - den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) - p1(3)*p2(3)

      cI =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      cI = cI - den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)
      endif

      endif

      c(0) = c0
      c(1) = c1
      c(2) = c2
      c(3) = c3
      c5(0) = c50
      c5(1) = c51
      c5(2) = c52
      c5(3) = c53

      nm1 = 1.d0/sqrt((1.d0,0.d0)*den12*den22)

      c  = c*nm1
      c5 = c5*nm1
      cI = cI*nm1

      return
      end
********************************************
      double complex function dot1cmom(p,q)
      implicit none
      double complex q(0:3)
      double precision p(0:3)
      dot1cmom=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
      return
      end
***************************************************
*** Levi-Civita tensor contracted with 2 real and 2 COMPLEX momenta
      double complex function epsfun2cmom(q1,q2,q3,q4)
      implicit none
      double complex q3(0:3),q4(0:3)
      double precision q1(0:3),q2(0:3)
      double complex e,ev(0:3),dot1cmom
      external dot1cmom

      e = q1(1) * (q2(2)*(q3(3)*q4(0) - q3(0)*q4(3))  +
     .             q2(3)*(q3(0)*q4(2) - q3(2)*q4(0))  +
     .             q2(0)*(q3(2)*q4(3) - q3(3)*q4(2))) +
     .    q1(2) * (q2(1)*(q3(0)*q4(3) - q3(3)*q4(0))  +
     .             q2(3)*(q3(1)*q4(0) - q3(0)*q4(1))  +
     .             q2(0)*(q3(3)*q4(1) - q3(1)*q4(3))) +
     .    q1(3) * (q2(1)*(q3(2)*q4(0) - q3(0)*q4(2))  +
     .             q2(2)*(q3(0)*q4(1) - q3(1)*q4(0))  +
     .             q2(0)*(q3(1)*q4(2) - q3(2)*q4(1))) +
     .    q1(0) * (q2(1)*(q3(3)*q4(2) - q3(2)*q4(3))  +
     .             q2(2)*(q3(1)*q4(3) - q3(3)*q4(1))  +
     .             q2(3)*(q3(2)*q4(1) - q3(1)*q4(2)))

      epsfun2cmom = -(0.d0,1.d0) * e!! e chittemmuorto!!!
      return
      end

**************************************************************************
      subroutine loadspinvariants(p1,p2,p3,p4,cmu,ce)
#include "invariants.h"
!     double precision f
!     c5X are in common!
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)

      p1cmu = dot1cmom(p1,cmu)
      p2cmu = dot1cmom(p2,cmu)
      p3cmu = dot1cmom(p3,cmu)
      p4cmu = dot1cmom(p4,cmu)
      p1ce  = dot1cmom(p1,ce)
      p2ce  = dot1cmom(p2,ce)
      p3ce  = dot1cmom(p3,ce)
      p4ce  = dot1cmom(p4,ce)

      e13cmuc5e = epsfun2cmom(p1,p3,cmu,c5e)
      e14cmuc5e = epsfun2cmom(p1,p4,cmu,c5e)
      e14cec5mu = epsfun2cmom(p1,p4,ce,c5mu)
      e34cec5mu = epsfun2cmom(p3,p4,ce,c5mu)
      e12cmuc5e = epsfun2cmom(p1,p2,cmu,c5e)
      e12cec5mu = epsfun2cmom(p1,p2,ce,c5mu)
      e23cec5mu = epsfun2cmom(p2,p3,ce,c5mu)

      p1c5mu = dot1cmom(p1,c5mu)
      p2c5mu = dot1cmom(p2,c5mu)
      p3c5mu = dot1cmom(p3,c5mu)
      p4c5mu = dot1cmom(p4,c5mu)
      p1c5e  = dot1cmom(p1,c5e)
      p2c5e  = dot1cmom(p2,c5e)
      p3c5e  = dot1cmom(p3,c5e)
      p4c5e  = dot1cmom(p4,c5e)

      c5muc5e = dotcmom(c5mu,c5e)

! if needed
      cmucmn = cmu
      cecmn  = ce
      return
      end
***

*********************************************************
      subroutine getQEDvertexes(vertexes)
      implicit none
      double complex C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu
      double complex C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte
      double complex C00ltmu,C00lte,vertexes,cover,amplvertexes
      double complex amplvertexesmu,amplvertexese,vertexmu,vertexe
      integer ifirst
      data ifirst /0/
      common/ifirstvert/cover,ifirst
#include "invariants.h"
      common/vertexfuns/
     .     C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu,
     .     C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte,
     .     C00ltmu,C00lte

      double complex amplve,amplvmu
      common/vertexeslegs/amplve,amplvmu


      if (ifirst.eq.0) then
         cover = ec4*im*pi216m1*4.d0
         ifirst = 1
      endif

      if (Qmu.lt.0.d0) then
c         include 'oneloop/vertexesminus.f'
      amplvmu =
     .  + QRmu * (  - 1.D0/2.D0*cmuce + cmuce*C00ltmu + cmuce*p1p3*
     .    C12ltmu + 2*cmuce*p1p3*C1ltmu + cmuce*p1p3*C0ltmu + mm2*cmuce
     .    *C11ltmu + 2*mm2*cmuce*C1ltmu - cImu*mm*p3ce*C12ltmu - cImu*
     .    mm*p3ce*C11ltmu - cImu*mm*p3ce*C1ltmu - cImu*mm*p1ce*C12ltmu
     .     - cImu*mm*p1ce*C11ltmu - cImu*mm*p1ce*C1ltmu )
      amplve =
     . + QRe * (  - 1.D0/2.D0*cmuce + cmuce*
     .    C00lte + cmuce*p2p4*C12lte + 2*cmuce*p2p4*C1lte + cmuce*p2p4*
     .    C0lte + me2*cmuce*C11lte + 2*me2*cmuce*C1lte - cIe*me*p4cmu*
     .    C12lte - cIe*me*p4cmu*C11lte - cIe*me*p4cmu*C1lte - cIe*me*
     .    p2cmu*C12lte - cIe*me*p2cmu*C11lte - cIe*me*p2cmu*C1lte )
      else
c         include 'oneloop/vertexesplus.f'
         amplvmu =
     .  + QRmu * (  - 1.D0/2.D0*cmuce + cmuce*C00ltmu + cmuce*p1p3*
     .    C12ltmu + 2*cmuce*p1p3*C1ltmu + cmuce*p1p3*C0ltmu + mm2*cmuce
     .    *C11ltmu + 2*mm2*cmuce*C1ltmu + cImu*mm*p3ce*C12ltmu + cImu*
     .    mm*p3ce*C11ltmu + cImu*mm*p3ce*C1ltmu + cImu*mm*p1ce*C12ltmu
     .     + cImu*mm*p1ce*C11ltmu + cImu*mm*p1ce*C1ltmu )
         amplve =
     .    + QRe * (  - 1.D0/2.D0*cmuce + cmuce*
     .    C00lte + cmuce*p2p4*C12lte + 2*cmuce*p2p4*C1lte + cmuce*p2p4*
     .    C0lte + me2*cmuce*C11lte + 2*me2*cmuce*C1lte - cIe*me*p4cmu*
     .    C12lte - cIe*me*p4cmu*C11lte - cIe*me*p4cmu*C1lte - cIe*me*
     .    p2cmu*C12lte - cIe*me*p2cmu*C11lte - cIe*me*p2cmu*C1lte )
      endif
      amplve  = amplve * tm1 * cover
      amplvmu = amplvmu * tm1 * cover

      vertexes = amplve+amplvmu

c      vertexes = amplvertexes * tm1 * cover
!     !vertexes = vertexes * vpt ! done outside, it's better
      return
      end
*********************************************************
**************************************************************************
      subroutine getQEDboxes(boxes)
c      implicit none
c      double complex Z1_,Z2_,Z3_,Z4_,Z5_,Z6_,Z7_
c      double complex Z8_,Z9_,Z10_,Z11_,Z12_,Z13_
c      double complex Z14_,Z15_,Z16_,Z17_,Z18_,Z19_,Z20_
c      double complex Z21_,Z22_,Z23_,Z24_,Z25_,Z26_
c      double complex Z27_
      implicit double complex (Z)
      double complex D0ltd,D1ltd,D2ltd,D3ltd
      double complex D0ltc,D1ltc,D2ltc,D3ltc
      double complex D00ltd,D11ltd,D12ltd,D13ltd
      double complex D22ltd,D23ltd,D33ltd
      double complex D00ltc,D11ltc,D12ltc,D13ltc
      double complex D22ltc,D23ltc,D33ltc
      double complex boxes,cover,box1
      integer ifirst
      double complex cu
      double precision zerop(0:3),k1(0:3),k2(0:3),k3(0:3)
      common/ifirstbox/cover,zerop,cu,ifirst
      data ifirst /0/
#include "invariants.h"
      common/boxesfuns/
     .     D0ltd,D1ltd,D2ltd,D3ltd,
     .     D0ltc,D1ltc,D2ltc,D3ltc,
     .     D00ltd,D11ltd,D12ltd,D13ltd,
     .     D22ltd,D23ltd,D33ltd,
     .     D00ltc,D11ltc,D12ltc,D13ltc,
     .     D22ltc,D23ltc,D33ltc

      integer rank
      parameter (rank = 2)
      double complex Dtend(0:rank,0:rank,0:rank,0:rank)
      double complex Dtenc(0:rank,0:rank,0:rank,0:rank)
      double complex Dtenduv(0:rank,0:rank,0:rank,0:rank)
      double complex Dtencuv(0:rank,0:rank,0:rank,0:rank)
      common/boxestenfuns/Dtend,Dtenc
      double complex dottencmom
      external dottencmom
      if (QRmu.eq.0.or.QRe.eq.0) then
         boxes = (0.d0,0.d0)
         return
      endif

      if (ifirst.eq.0) then
         cover = ec4*im*pi216m1*2.d0
         cu    = (1.d0,0.d0)
         zerop = 0.d0
         ifirst = 1
      endif

      if (Qmu.lt.0.d0) then
         include 'oneloop/boxesminus.f'
      else
         include 'oneloop/boxesplus.f'
      endif
      boxes = boxes * cover
      return
      end
***************************************************************************
************* END FROM MESMER
#endif
      subroutine exitlooplibraries
#ifdef LOOPTOOLS
#include "looptools.h"
      call ltexi
#endif
      return
      end
**************************************************************************
      subroutine get_twiceVPtimesTreeOverBorn(p1,p2,p3,p4,smtob)
      implicit double precision (a-h,l,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),ptmp(0:3)
      complex*16 chis,chit,im
      parameter (im = (0.d0,1.d0))
      double complex vpolc,ec2s
      external vpolc
      character*2 fs
      common/finalstate/fs
      common/parameters/ame,ammu,convfac,alpha,pi
      common/zparameters/zm,gz,ve,ae,rv,ra,wm,s2th,gfermi,sqrt2,um

      integer ifirst
      common/ifirstgetdalphaFO/ifirst
      data ifirst /0/

      common/ialpharunning/iarun
      double complex dalep,dahadr
      common/forfixedorder/dalep,dahadr
      double complex cx,daleps,dahadrs,dalept,dahadrt
      double complex Das,Dat

      double complex cvpt_new,cvpt_hig,cvpt_low,alphac,alphah,
     &     alphal,alphact,alphaht,alphalt,calept1,cahadr,caltop


      smtob = 0d0
      if (fs.eq.'gg') return

      smt   = 0.d0

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

      if (fs.eq.'ee') then

         ibck = iarun
         iarun = 1
         cx = vpolc(s)
         daleps  = dalep
         dahadrs = dahadr
         cx = vpolc(t)
         dalept  = dalep
         dahadrt = dahadr
         iarun = ibck

         Dat = dalept + dahadrt
         Das = daleps + dahadrs

         ec2t  = ec2
         ec2s  = ec2

         sm2   = abs(ec2s)**2/s2*((t2+u2)/4.d0+2.d0*em2*s
     .        -2.d0*em4)
         tm2   = ec2t**2/t2*((s2+u2)/4.d0+2.d0*em2*t
     .        -2.d0*em4)
         st2   = ec2s*ec2t/s/t*(u2/2.d0-4.d0*em2*u+6.d0*em4)

         smtq2 = 2.d0*(tm2*Dat + sm2*Das + 0.5d0*st2*(Das+Dat))  ! 1/2 perche' st2 e' 2*s*t!!

         born = 32.d0*(tm2 + sm2 + st2)

         smtob   = smtq2*32.d0/born

c         call alphaQEDcx(s,cvpt new,cvpt hig,cvpt low,alphac,alphah,
c     &        alphal,alphact,alphaht,alphalt,calept,cahadr,caltop)
c         print*,alphac,alpha/(1.d0-Das-Dat)
c         stop

      elseif (fs.eq.'mm'.or.fs.eq.'pp'.or.fs.eq.'pi') then
         ibck = iarun
         iarun = 1
         cx = vpolc(s)
         daleps  = dalep
         dahadrs = dahadr
         iarun = ibck

         Das = daleps + dahadrs

c         me2   = ame*ame
c         e4    = abs(ec2s)**2
c         elmat2=32.d0*e4/s2*((t2+u2)/4.d0+s*(me2+mm2)-(me2+mm2)**2/2.d0)
         smtob   = 2.d0*Das
         return
      endif
********************
      return
      end



*****************************************
      subroutine qedvertex_isr(p1,p2,p3,p4,virtual)
#ifdef COLLIER
        use collier
#endif
        implicit none
        double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
        double precision qph(40,0:3)
        double precision virtual,virtual1
        real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
        real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45,xsum
        real*8 e2,e6,me2,me4,me6,mm2,mm4,mm6
        double precision ame,ammu,convfac,alpha,pi,ame2,ammu2,amu
        common/parameters/ame,ammu,convfac,alpha,pi
        double complex cnpi(0:1,0:2,0:2)
        double complex elmat
        double complex cnpiuv(0:1,0:2,0:2)
        external dot
        double precision dot,born
        external  soft_integral_muons
        double precision soft_integral_muons
        Complex*16 A0ae
        Complex*16 B0ax12ee
        Complex*16 B0ae0e

        Complex*16 C00aeex12e0e
	Complex*16 C0aeex12e0e
        Complex*16 C2aeex12e0e
        Complex*16 C12aeex12e0e
        Complex*16 C22aeex12e0e
        Complex*16 C12aex12e0ee
        Complex*16 DB0e0e
        complex*16 dpsi
        double precision egmax
        common/egmaxforrecola/egmax
        double complex c1
        parameter (c1 = (1.d0,0.d0))
        double complex b1
        parameter (b1 = (1.d0,0.d0))
        double complex a1
	parameter (a1 = (1.d0,0.d0))

        e2 = 4.d0*pi*alpha
        e6 = e2*e2*e2
                ame2=ame*ame
                me2=ame2
                amu=ammu
                ammu2=ammu*ammu
                mm2=ammu2
                p1p2 = dot(p1,p2)
                p1p3 = dot(p1,p3)
                p1p4 = dot(p1,p4)
                p2p3 = dot(p2,p3)
                p2p4 = dot(p2,p4)
                p3p4 = dot(p3,p4)

        x12  = 2d0*ame2+2.d0*p1p2
#ifdef COLLIER

        call DB0_cll(DB0e0e,b1*ame2,b1*0d0,b1*ame2)
        call A0_cll(A0ae,a1*ame2)
	call B0_cll(B0ae0e,b1*ame2,b1*0d0,b1*ame2)

        call B0_cll(B0ax12ee,b1*x12,b1*ame2,b1*ame2)


       call C_cll(cnpi,cnpiuv,c1*ame2,c1*ame2,c1*x12,c1*ame2,
     >     c1*0d0,c1*ame2,2)
       C00aeex12e0e=cnpi(1,0,0)
       C0aeex12e0e=cnpi(0,0,0)
       C2aeex12e0e=cnpi(0,0,1)
       C12aeex12e0e=cnpi(0,1,1)
       C22aeex12e0e=cnpi(0,0,2)

       call C_cll(cnpi,cnpiuv,c1*ame2,c1*x12,c1*ame2,
     >     c1*0d0,c1*ame2,c1*ame2,2)
       C12aex12e0ee=cnpi(0,1,1)
#endif

       dpsi =(4d0*ame2**2*DB0e0e-A0ae+ame2)/ame2


              elmat = -(((p1p4*p2p3 + p1p3*p2p4 +
     &        ame**2*(p3p4 + amu**2*2D0) +
     &        B0ax12ee*2D0*
     &         (p1p4*p2p3 + p1p3*p2p4 +
     &           ame**2*(p3p4 + amu**2*2D0)) +
     &        amu**2*C0aeex12e0e*p1p2**2*4D0 -
     &        B0ae0e*(p1p4*p2p3 + p1p3*p2p4 +
     &           ame**2*(p3p4 + amu**2*2D0))*4D0 +
     &        (ame**2*(C12aex12e0ee + C22aeex12e0e)*
     &            (ame**2*(amu**2 + p3p4) -
     &              1D0*(p1p3 - p2p3*1D0)*(p1p4 - p2p4*1D0))
     &            + ame**2*C2aeex12e0e*
     &            (ame**2*(amu**2 + p3p4) -
     &              1D0*(p1p3 - p2p3*1D0)*(p1p4 - p2p4*1D0))
     &            + C00aeex12e0e*
     &            (p1p4*p2p3 + p1p3*p2p4 +
     &              ame**2*(p3p4 + amu**2*2D0)))*4D0 +
     &        p1p2*(amu**2 + amu**2*B0ax12ee*2D0 -
     &           amu**2*B0ae0e*4D0 +
     &           amu**2*C00aeex12e0e*4D0 -
     &           ame**2*amu**2*C12aex12e0ee*4D0 -
     &           ame**2*amu**2*C2aeex12e0e*4D0 -
     &           ame**2*C12aex12e0ee*p3p4*4D0 -
     &           ame**2*C2aeex12e0e*p3p4*4D0 -
     &           ame**2*C22aeex12e0e*(amu**2 + p3p4)*4D0 +
     &           C0aeex12e0e*
     &            (p1p4*p2p3 + p1p3*p2p4 +
     &              ame**2*(p3p4 + amu**2*2D0))*4D0))*8D0)/
     &      (amu**2*2D0 + p3p4*2D0)**2)

              born=(2D0*(amu**2*p1p2 + p1p4*p2p3 + p1p3*p2p4 +
     &             ame**2*(p3p4 + amu**2*2D0)))/(amu**2 + p3p4)**2
              virtual=4d0*(elmat+dpsi*born)*e2*e2*alpha/(2d0*pi)/
     &         (4d0*born*e2*e2) +
     &             soft_integral_muons(egmax,p1,p2,p3,p4)

       end
****
***********************************

      subroutine virtualISR_mur(p1,p2,p3,p4,q,virtual,qph,ifc)

c IN QUESTA ROUTINE CALCOLO IL CONTRIBUTO DI INTERFERENZA TRA IL VIRTUALE ISR E IL BORN, NELLA PRIMA PARTE CALCOLO L'INTERFERENZA CON ISR DEL BORN, MENTRE NELLA SECONDA CALCOLO QUELLA
C CON IL BORN FSR, POI PER SIMMETRIA PORTO TUTTO QUESTO AL SIMMETRICO, OSSIA VIRTUALE FSR-BORN FSR E VIRTUALE FSR BORN ISR
#ifdef COLLIER
        use collier
#endif
	implicit double complex (z)
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision qph(40,0:3)
	double precision virtual,virtual1
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45,xsum
	real*8 e2,e6,me2,me4,me6,mm2,mm4,mm6
	double precision ame,ammu,convfac,alpha,pi,ame2,ammu2
	common/parameters/ame,ammu,convfac,alpha,pi
	Complex*16 elmat,elmat1,elmat2,elmat3,elmat4,elmat5,elmat6
        Complex*16 elmat7,elmat8,elmat9,elmat10,elmat11
	character*2 fs
	character*3 ifc
	common/finalstate/fs
	! using always FORM !!!!!!!!!
	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min
	real*16 born,oneloop
	double complex cnpi(0:1,0:2,0:2)
	double complex cnpiuv(0:1,0:2,0:2)
	double complex dnpi(0:1,0:3,0:3,0:3)
	double complex dnpiuv(0:1,0:3,0:3,0:3)
	common/scaling/scaling
	double complex scaling
	common/softphotonmass/phmass
	real*8 phmass

        common/radoption/ifisrfsr
        integer ifisrfsr

        external vpolc
        double complex vpolc,alphar

        common/ialpharunning/iarun
        integer iarun

	external dot
	real*8 dot

        Complex*16 dpsi
        Complex*16 dmu
        COmplex*16 dme

	double complex d1
	parameter (d1 = (1.d0,0.d0))
	double complex c1
	parameter (c1 = (1.d0,0.d0))
	double complex b1
	parameter (b1 = (1.d0,0.d0))
	double complex a1
 	parameter (a1 = (1.d0,0.d0))

 	Complex*16 A0ae
	Complex*16 B0ax150e
	Complex*16 B0ax250e
	Complex*16 B0ax34ee
	Complex*16 B0ax12ee
	Complex*16 B0a0ee
	Complex*16 B0ae0e
	Complex*16 DB0e0e

	Complex*16 C00aex150e0e
        Complex*16 C1aex150e0e
        Complex*16 C11aex150e0e
        Complex*16 C0aex150e0e
        Complex*16 C0ax150e0ee
        Complex*16 C2aex150e0e
        Complex*16 C0a0ex15ee0

        Complex*16 C0aex250e0e
	Complex*16 C11aex250e0e
	Complex*16 C1aex250e0e
	Complex*16 C2aex250e0e
	Complex*16 C00aex250e0e
	Complex*16 C0a0ex25ee0
       Complex*16 C0ax250e0ee

        Complex*16 C1ax15x34e0ee

        Complex*16 C00aex25x34e0e
	Complex*16 C12aex25x34e0e
        Complex*16 C11aex25x34e0e
        Complex*16 C1aex25x34e0e
        Complex*16 C2aex25x34e0e
	Complex*16 C0aex25x34e0e
	Complex*16 C1ax25x34e0ee
	Complex*16 C22aex25x34e0e
	Complex*16 C0ax25x34e0ee
       Complex*16 C0aex34x250ee

       Complex*16 C0aex15x34e0e
       Complex*16 C11aex15x34e0e
       Complex*16 C12aex15x34e0e
       Complex*16 C1aex15x34e0e
       Complex*16 C2aex15x34e0e
       Complex*16 C00aex15x34e0e
       Complex*16 C0ax15x34e0ee
       Complex*16 C22aex15x34e0e

	Complex*16 C1ax12x340eee
        Complex*16 C0a0x12x34eee

	Complex*16 C0aex12e0ee
        Complex*16 C1aex12e0ee
        Complex*16 C0aeex12e0e
        Complex*16 C00aeex12e0e
        Complex*16 C12aex12e0ee
        Complex*16 C22aeex12e0e
        Complex*16 C2aeex12e0e

	Complex*16 D001aeex340x12x15e0ee
	Complex*16 D002aeex340x12x15e0ee
	Complex*16 D003aeex340x12x15e0ee
	Complex*16 D00aeex340x12x15e0ee
	Complex*16 D3aeex340x12x15e0ee
	Complex*16 D0aeex340x12x15e0ee
	Complex*16 D33aeex340x12x15e0ee
	Complex*16 D333aeex340x12x15e0ee
	Complex*16 D12aeex340x12x15e0ee
	Complex*16 D112aeex340x12x15e0ee
	Complex*16 D22aeex340x12x15e0ee
	Complex*16 D11aeex340x12x15e0ee
	Complex*16 D2aeex340x12x15e0ee
	Complex*16 D111aeex340x12x15e0ee
	Complex*16 D122aeex340x12x15e0ee
	Complex*16 D1aeex340x12x15e0ee
	Complex*16 D123aeex340x12x15e0ee
	Complex*16 D23aeex340x12x15e0ee
	Complex*16 D13aeex340x12x15e0ee
	Complex*16 D113aeex340x12x15e0ee
	Complex*16 D223aeex340x12x15e0ee
	Complex*16 D133aeex340x12x15e0ee
	Complex*16 D233aeex340x12x15e0ee


	Complex*16 D113aeex340x12x25e0ee
	Complex*16 D112aeex340x12x25e0ee
	Complex*16 D003aeex340x12x25e0ee
	Complex*16 D002aeex340x12x25e0ee
	Complex*16 D00aeex340x12x25e0ee
	Complex*16 D11aeex340x12x25e0ee
	Complex*16 D1aeex340x12x25e0ee
	Complex*16 D133aeex340x12x25e0ee
	Complex*16 D0aeex340x12x25e0ee
	Complex*16 D3aeex340x12x25e0ee
	Complex*16 D333aeex340x12x25e0ee
	Complex*16 D33aeex340x12x25e0ee
	Complex*16 D22aeex340x12x25e0ee
	Complex*16 D12aeex340x12x25e0ee
	Complex*16 D2aeex340x12x25e0ee
	Complex*16 D122aeex340x12x25e0ee
	Complex*16 D123aeex340x12x25e0ee
	Complex*16 D13aeex340x12x25e0ee
	Complex*16 D23aeex340x12x25e0ee
	Complex*16 D223aeex340x12x25e0ee
	Complex*16 D233aeex340x12x25e0ee
	Complex*16 D111aeex340x12x25e0ee
	Complex*16 D001aeex340x12x25e0ee



	e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2
	if(ifc.eq.'isr') then
                ame2=ame*ame
                me2=ame2
		ammu2=ammu*ammu
                mm2=ammu2
		p1p2 = dot(p1,p2)
		p1p3 = dot(p1,p3)
		p1p4 = dot(p1,p4)
		kp1  = dot(p1,q)
		p2p3 = dot(p2,p3)
		p2p4 = dot(p2,p4)
		kp2  = dot(p2,q)
		kp3 = dot(p3,q)
		p3p4 = dot(p3,p4)
		kp4  = dot(p4,q)
	elseif(ifc.eq.'fsr') then
                ame2=ammu*ammu
                me2=ame2
		ammu2=ame*ame
                mm2=ammu2
		p1p2 = dot(p3,p4)
		p1p3 = dot(p1,p3)
		p1p4 = dot(p3,p2)
		kp1  = -dot(p3,q)
		p2p3 = dot(p1,p4)
		p2p4 = dot(p2,p4)
		kp2  = -dot(p4,q)
		kp3 = -dot(p1,q)
		p3p4 = dot(p1,p2)
		kp4  = -dot(p2,q)
	endif
!
	x12  = 2d0*ame2+2.d0*p1p2
	x13  = ame2+ammu2-2.d0*p1p3
	x14  = ame2+ammu2-2.d0*p1p4
	x15  = ame2-2.d0*kp1
	x23  = ame2+ammu2-2d0*p2p3
	x24  = ame2+ammu2-2d0*p2p4
	x25 = ame2-2d0*kp2
	x34 = 2d0*ammu2+2d0*p3p4
	x35 = ammu2+2d0*kp3
	x45 = ammu2+2d0*kp4


        mm4=mm2**2
        mm6=mm2**3
        me4=me2**2
        me6=me2**3
c$$$
c$$$        x34=    0.8934752924784175D0
c$$$        x12= 1.040400000000000D0
c$$$        x25=0.3947831043530346d-15
c$$$
#ifdef COLLIER
!     call SetDeltaUV_cll(1000000d0)

        call A0_cll(A0ae,a1*ame2)

	call DB0_cll(DB0e0e,b1*ame2,b1*0d0,b1*ame2)

        call B0_cll(B0ax250e,b1*x25,b1*0d0,b1*ame2)

	call B0_cll(B0ax150e,b1*x15,b1*0d0,b1*ame2)

	call B0_cll(B0ax34ee,b1*x34,b1*ame2,b1*ame2)

	call B0_cll(B0ae0e,b1*ame2,b1*0d0,b1*ame2)

	call B0_cll(B0a0ee,b1*0d0,b1*ame2,b1*ame2)

	call B0_cll(B0ax12ee,b1*x12,b1*ame2,b1*ame2)


       call C_cll(cnpi,cnpiuv,c1*ame2,c1*x15,c1*0d0,c1*ame2,
     >     c1*0d0,c1*ame2,2)
     		C1aex150e0e=cnpi(0,1,0)
       	        C00aex150e0e=cnpi(1,0,0)
      	 	C11aex150e0e=cnpi(0,2,0)
      	 	C0aex150e0e=cnpi(0,0,0)
      	 	C2aex150e0e=cnpi(0,0,1)
      	 	C0ax150e0ee=cnpi(0,0,0)
      	 	C0a0ex15ee0=cnpi(0,0,0)

       call C_cll(cnpi,cnpiuv,c1*ame2,c1*x25,c1*0d0,c1*ame2,
     >     c1*0d0,c1*ame2,2)
		C0aex250e0e=cnpi(0,0,0)
		C11aex250e0e=cnpi(0,2,0)
		C1aex250e0e=cnpi(0,1,0)
		C2aex250e0e=cnpi(0,0,1)
		C00aex250e0e=cnpi(1,0,0)
		C0ax250e0ee=cnpi(0,0,0)
                C0a0ex25ee0=cnpi(0,0,0)

       call C_cll(cnpi,cnpiuv,c1*ame2,c1*x25,c1*x34,c1*ame2,
     >       c1*0d0,c1*ame2,2)
      	C12aex25x34e0e	=cnpi(0,1,1)
       	C11aex25x34e0e =cnpi(0,2,0)
       	C1aex25x34e0e =cnpi(0,1,0)
       	C2aex25x34e0e=cnpi(0,0,1)
       	C00aex25x34e0e =cnpi(1,0,0)
       	C0aex25x34e0e =cnpi(0,0,0)
       	C0ax25x34e0ee=cnpi(0,0,0)

		C1ax25x34e0ee=cnpi(0,0,1)

		C22aex25x34e0e =cnpi(0,0,2)
		C0aex34x250ee =cnpi(0,0,0)

       call C_cll(cnpi,cnpiuv,c1*ame2,c1*x15,c1*x34,c1*ame2,
     >       c1*0d0,c1*ame2,2)
       	C0aex15x34e0e=cnpi(0,0,0)
       	C11aex15x34e0e=cnpi(0,2,0)
       	C12aex15x34e0e =cnpi(0,1,1)
       	C1aex15x34e0e =cnpi(0,1,0)
       	C2aex15x34e0e =cnpi(0,0,1)
       	C00aex15x34e0e=cnpi(1,0,0)

		C1ax15x34e0ee=cnpi(0,0,1)
		C0ax15x34e0ee =cnpi(0,0,0)
		C22aex15x34e0e =cnpi(0,0,2)

       call C_cll(cnpi,cnpiuv,c1*x12,c1*x34,c1*0d0,c1*ame2,
     >       c1*ame2,c1*ame2,2)
       	C1ax12x340eee=cnpi(0,1,0)
       	C0a0x12x34eee=cnpi(0,0,0)

       call C_cll(cnpi,cnpiuv,c1*ame2,c1*x12,c1*ame2,c1*0d0,
     >       c1*ame2,c1*ame2,2)
       	C1aex12e0ee=  cnpi(0,1,0)
       	C0aex12e0ee=  cnpi(0,0,0)
       	C12aex12e0ee=cnpi(0,1,1)

	call D_cll(dnpi,dnpiuv,d1*ame2,d1*ame2,d1*x34,d1*0d0,d1*x12,
     >    d1*x25,d1*ame2,d1*0d0,d1*ame2,d1*ame2,3)
	D113aeex340x12x25e0ee = dnpi(0,2,0,1)
	D112aeex340x12x25e0ee = dnpi(0,2,1,0)
	D003aeex340x12x25e0ee = dnpi(1,0,0,1)
	D002aeex340x12x25e0ee = dnpi(1,0,1,0)
	D00aeex340x12x25e0ee = dnpi(1,0,0,0)
	D11aeex340x12x25e0ee = dnpi(0,2,0,0)
	D1aeex340x12x25e0ee  = dnpi(0,1,0,0)
	D133aeex340x12x25e0ee = dnpi(0,1,0,2)
	D0aeex340x12x25e0ee = dnpi(0,0,0,0)
	D3aeex340x12x25e0ee = dnpi(0,0,0,1)
	D333aeex340x12x25e0ee = dnpi(0,0,0,3)
	D33aeex340x12x25e0ee = dnpi(0,0,0,2)
	D22aeex340x12x25e0ee = dnpi(0,0,2,0)
	D12aeex340x12x25e0ee = dnpi(0,1,1,0)
	D2aeex340x12x25e0ee = dnpi(0,0,1,0)
	D122aeex340x12x25e0ee = dnpi(0,1,2,0)
	D123aeex340x12x25e0ee = dnpi(0,1,1,1)
	D13aeex340x12x25e0ee = dnpi(0,1,0,1)
	D23aeex340x12x25e0ee = dnpi(0,0,1,1)
	D223aeex340x12x25e0ee = dnpi(0,0,2,1)
	D233aeex340x12x25e0ee = dnpi(0,0,1,2)
        D111aeex340x12x25e0ee = dnpi(0,3,0,0)
	D001aeex340x12x25e0ee = dnpi(1,1,0,0)


	call D_cll(dnpi,dnpiuv,d1*ame2,d1*ame2,d1*x34,d1*0d0,d1*x12,
     >    d1*x15,d1*ame2,d1*0d0,d1*ame2,d1*ame2,3)
	D001aeex340x12x15e0ee =dnpi(1,1,0,0)
	D002aeex340x12x15e0ee =dnpi(1,0,1,0)
	D003aeex340x12x15e0ee = dnpi(1,0,0,1)
	D00aeex340x12x15e0ee = dnpi(1,0,0,0)
	D3aeex340x12x15e0ee = dnpi(0,0,0,1)
	D0aeex340x12x15e0ee = dnpi(0,0,0,0)
	D33aeex340x12x15e0ee = dnpi(0,0,0,2)
	D333aeex340x12x15e0ee = dnpi(0,0,0,3)
	D12aeex340x12x15e0ee = dnpi(0,1,1,0)
	D112aeex340x12x15e0ee = dnpi(0,2,1,0)
	D22aeex340x12x15e0ee =dnpi(0,0,2,0)
	D11aeex340x12x15e0ee =dnpi(0,2,0,0)
	D2aeex340x12x15e0ee =dnpi(0,0,1,0)
	D111aeex340x12x15e0ee =dnpi(0,3,0,0)
	D122aeex340x12x15e0ee=dnpi(0,1,2,0)
	D1aeex340x12x15e0ee=dnpi(0,1,0,0)
	D123aeex340x12x15e0ee = dnpi(0,1,1,1)
	D23aeex340x12x15e0ee = dnpi(0,0,1,1)
	D13aeex340x12x15e0ee = dnpi(0,1,0,1)
	D113aeex340x12x15e0ee = dnpi(0,2,0,1)
	D223aeex340x12x15e0ee = dnpi(0,0,2,1)
	D133aeex340x12x15e0ee = dnpi(0,1,0,2)
	D233aeex340x12x15e0ee = dnpi(0,0,1,2)

#endif

****************(VIRTUALE ISR REALE ISR) CROSSATO CON BORN ISR

************NEW WITH FORM

        include 'mumurad/Virtual_new/isrisr.f'


!******************CT dpsi (fattorizza sul born con dpsi)

	dpsi =(4d0*ame2**2*DB0e0e-A0ae+ame2)/ame2

	if(ifc.eq.'isr') then
		call ISR_mur(p1,p2,p3,p4,q,elmat11,qph)
	elseif(ifc.eq.'fsr') then
		call FSR_mur(p1,p2,p3,p4,q,elmat11,qph)
	endif

	elmat11=elmat11*dpsi/(4d0*e6)

***********************

        virtual=real(elmat*vpolc(x34)*CONJG(vpolc(x34))+elmat11)






******QUI INIZIO CON (VIRTUALE ISR REALE ISR)  CROSSATO CON FSR BORN

        dme=-(ame2*(2d0*B0ae0e-1d0)+A0ae)/ame2

******+NEW WITH FORM

        include 'mumurad/Virtual_new/isrfsr.f'


******IL CT FATTORIZZA SUL BORN, CHE PERÒ IN QUESTO CASO È IFI
	dpsi =(4d0*ame2**2*DB0e0e-A0ae+ame2)/ame2

	call IFI_mur(p1,p2,p3,p4,q,elmat11,qph)
	elmat11=elmat11*dpsi/(4d0*e6)


*****************************


        if(ifisrfsr.eq.11) then

           virtual=virtual+real(elmat*vpolc(x34)*CONJG(vpolc(x12))+
     >      elmat11)

	else
	    virtual=virtual

	endif

	virtual=4d0*virtual*e6*alpha/(2d0*pi)


      end
*****************************************************


      subroutine virtualMixed_mur(p1,p2,p3,p4,q,virtual,qph,ifc)

#ifdef COLLIER
        use collier
#endif
	implicit double complex (z)
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision qph(40,0:3)
	double precision virtual,virtual1
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45,xsum
	real*8 e2,e6,me2,me4,me6,mm2,mm4,mm6
	double precision ame,ammu,convfac,alpha,pi,ame2,ammu2
	common/parameters/ame,ammu,convfac,alpha,pi
	double precision amme,mmu
	Complex*16 elmat,elmat1,elmat2,elmat3,elmat4,elmat5,elmat6
        Complex*16 elmat7,elmat8,elmat9,elmat10,elmat11
	character*3 ifc
	common/finalstate/fs
	character*2 fs

	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min

        external vpolc
        double complex vpolc,alphar

        common/ialpharunning/iarun
        integer iarun

	double complex cnpi(0:1,0:2,0:2)
	double complex cnpiuv(0:1,0:2,0:2)

	common/softphotonmass/phmass
	real*8 phmass

	external dot
	real*8 dot

	double complex b1
	parameter (b1 = (1.d0,0.d0))
	double complex c1
	parameter (c1 = (1.d0,0.d0))


	Complex*16  C0aeex12e0e
	Complex*16  C00aeex12e0e
	Complex*16  C12aeex12e0e
	Complex*16  C22aeex12e0e
	Complex*16  C2aeex12e0e
	Complex*16 C12aex12e0ee

	Complex*16  B0ae0e
	Complex*16  B0ax12ee
	Complex*16  DB0e0e

	Complex*16 A0ae

	double complex dpsi
ccc IN THIS SUBROUTINE WE CALCULATE isr_virtual real fsr crossed with born and vicevera

	e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2

	if(ifc.eq.'isr') then
		amme=ame*1d0
	 	ame2=ame*ame
                me2=ame2
		mmu=ammu*1d0
		ammu2=ammu*ammu
                mm2=ammu2
		p1p2 = dot(p1,p2)
		p1p3 = dot(p1,p3)
		p1p4 = dot(p1,p4)
		kp1  = dot(p1,q)
		p2p3 = dot(p2,p3)
		p2p4 = dot(p2,p4)
		kp2  = dot(p2,q)
		kp3 = dot(p3,q)
		p3p4 = dot(p3,p4)
		kp4  = dot(p4,q)
	elseif(ifc.eq.'fsr') then
		amme=ammu*1d0
		mmu=ame
		ame2=ammu*ammu
                me2=ame2
		ammu2=ame*ame
                mm2=ammu2
		p1p2 = dot(p3,p4)
		p1p3 = dot(p1,p3)
		p1p4 = dot(p3,p2)
		kp1  = -dot(p3,q)
		p2p3 = dot(p1,p4)
		p2p4 = dot(p2,p4)
		kp2  = -dot(p4,q)
		kp3 = -dot(p1,q)
		p3p4 = dot(p1,p2)
		kp4  = -dot(p2,q)
	endif

	x12  = 2d0*ame2+2.d0*p1p2
	x13  = ame2+ammu2-2.d0*p1p3
	x14  = ame2+ammu2-2.d0*p1p4
	x15  = ame2-2.d0*kp1
	x23  = ame2+ammu2-2d0*p2p3
	x24  = ame2+ammu2-2d0*p2p4
	x25 = ame2-2d0*kp2
	x34 = 2d0*ammu2+2d0*p3p4
	x35 = ammu2+2d0*kp3
	x45 = ammu2+2d0*kp4

        mm4=mm2**2
        mm6=mm2**3
        me4=me2**2
        me6=me2**3

#ifdef COLLIER
!     call SetDeltaUV_cll(1000000d0)

	call A0_cll(A0ae,b1*ame2)

	call DB0_cll(DB0e0e, b1*ame2,b1*0d0,b1*ame2)

        call B0_cll(B0ae0e,b1*ame2,b1*0d0,b1*ame2)

        call B0_cll(B0ax12ee,b1*x12,b1*ame2,b1*ame2)

        call C_cll(cnpi,cnpiuv,c1*ame2,c1*ame2,c1*x12,
     >       c1*ame2,c1*0d0,c1*ame2,2)

	C0aeex12e0e =  cnpi(0,0,0)
	C00aeex12e0e = cnpi(1,0,0)
	C12aeex12e0e =cnpi(0,1,1)
	C22aeex12e0e =cnpi(0,0,2)
	C2aeex12e0e =cnpi(0,0,1)

        call C_cll(cnpi,cnpiuv,c1*ame2,c1*x12,c1*ame2,
     >       c1*0d0,c1*ame2,c1*ame2,2)

	C12aex12e0ee=cnpi(0,1,1)

#endif

cc        qui ho virtuale isr reale fsr crossato con born isr e fsr
***********************

********NEW WITH FORM

	include 'mumurad/Virtual_new/mixedisr_fsr.f'

******

	virtual=real(elmat1*vpolc(x12)*CONJG(vpolc(x34))
     >		      +elmat2*vpolc(x12)*CONJG(vpolc(x12)))

*******************CT
	dpsi =(4d0*ame2**2*DB0e0e-A0ae+ame2)/ame2

	if(ifc.eq.'isr') then
 		call FSR_mur(p1,p2,p3,p4,q,elmat1,qph)
 		call IFI_mur(p1,p2,p3,p4,q,elmat2,qph)
	elseif(ifc.eq.'fsr') then
 		call ISR_mur(p1,p2,p3,p4,q,elmat1,qph)
 		call IFI_mur(p1,p2,p3,p4,q,elmat2,qph)
	endif

	virtual=virtual+dpsi*real(elmat1+elmat2)/(4d0*e6)

	virtual=4d0*virtual*e6*alpha/(2d0*pi)

	end



*****************************************
      subroutine virtualBox_mur(p1,p2,p3,p4,q,virtual,qph,ifc)

#ifdef COLLIER
        use collier
#endif
	implicit double complex (z)
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision qph(40,0:3)
	double precision virtual,virtual1
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45,xsum
	real*8 e2,e6,me2,me4,me6,mm2,mm4,mm6
	double precision ame,ammu,convfac,alpha,pi,ame2,ammu2
	common/parameters/ame,ammu,convfac,alpha,pi
	double precision amme,mmu
	Complex*16 elmat,elmat1,elmat2,elmat3,elmat4,elmat5,elmat6
        Complex*16 elmat7,elmat8,elmat9,elmat10,elmat11
	character*2 fs
	character*3 ifc
	common/finalstate/fs
	! using always FORM !!!!!!!!!
	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min
	real*16 born,oneloop
	double complex dnpi(0:1,0:2,0:2,0:2)
	double complex dnpiuv(0:1,0:2,0:2,0:2)
	common/scaling/scaling
	double complex scaling
	common/softphotonmass/phmass
	real*8 phmass

	external dot
	real*8 dot

        external vpolc
        double complex vpolc,alphar

        common/ialpharunning/iarun
        integer iarun

        Complex*16 dpsi
        Complex*16 dmu
        Complex*16 dme

        Complex*16 D0apx15epx23x34p0e0
	Complex*16 D3apx15epx23x34p0e0
	Complex*16 D00apx15epx23x34p0e0
	Complex*16 D11apx15epx23x34p0e0
	Complex*16 D13apx15epx23x34p0e0
	Complex*16 D33apx15epx23x34p0e0
	Complex*16 D22apx15epx23x34p0e0
	Complex*16 D23apx15epx23x34p0e0
	Complex*16 D2apx15epx23x34p0e0
	Complex*16 D1apx15epx23x34p0e0
	Complex*16 D12apx15epx23x34p0e0

	Complex*16 D0apx15epx24x34p0e0
	Complex*16 D3apx15epx24x34p0e0
	Complex*16 D00apx15epx24x34p0e0
	Complex*16 D11apx15epx24x34p0e0
	Complex*16 D13apx15epx24x34p0e0
	Complex*16 D33apx15epx24x34p0e0
	Complex*16 D2apx15epx24x34p0e0
	Complex*16 D23apx15epx24x34p0e0
	Complex*16 D22apx15epx24x34p0e0
	Complex*16 D12apx15epx24x34p0e0
	Complex*16 D1apx15epx24x34p0e0

	Complex*16 D0apx25epx13x34p0e0
	Complex*16 D3apx25epx13x34p0e0
	Complex*16 D00apx25epx13x34p0e0
	Complex*16 D11apx25epx13x34p0e0
	Complex*16 D13apx25epx13x34p0e0
	Complex*16 D33apx25epx13x34p0e0
	Complex*16 D2apx25epx13x34p0e0
	Complex*16 D23apx25epx13x34p0e0
	Complex*16 D22apx25epx13x34p0e0
	Complex*16 D12apx25epx13x34p0e0
	Complex*16 D1apx25epx13x34p0e0

	Complex*16 D0apx25epx14x34p0e0
	Complex*16 D3apx25epx14x34p0e0
	Complex*16 D00apx25epx14x34p0e0
	Complex*16 D11apx25epx14x34p0e0
	Complex*16 D13apx25epx14x34p0e0
	Complex*16 D33apx25epx14x34p0e0
	Complex*16 D2apx25epx14x34p0e0
	Complex*16 D22apx25epx14x34p0e0
	Complex*16 D23apx25epx14x34p0e0
	Complex*16 D1apx25epx14x34p0e0
	Complex*16 D12apx25epx14x34p0e0

	double complex d1
	parameter (d1 = (1.d0,0.d0))

	e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2


	if(ifc.eq.'isr') then
		amme=ame*1d0
		ame2=ame*ame
		mmu=ammu*1d0
		ammu2=ammu*ammu
                me2=ame2
                mm2=ammu2
		p1p2 = dot(p1,p2)
		p1p3 = dot(p1,p3)
		p1p4 = dot(p1,p4)
		kp1  = dot(p1,q)
		p2p3 = dot(p2,p3)
		p2p4 = dot(p2,p4)
		kp2  = dot(p2,q)
		kp3 = dot(p3,q)
		p3p4 = dot(p3,p4)
		kp4  = dot(p4,q)
	elseif(ifc.eq.'fsr') then
		amme=ammu*1d0
		mmu=ame
		ame2=ammu*ammu
		ammu2=ame*ame
                me2=ame2
                mm2=ammu2
		p1p2 = dot(p3,p4)
		p1p3 = dot(p1,p3)
		p1p4 = dot(p3,p2)
		kp1  = -dot(p3,q)
		p2p3 = dot(p1,p4)
		p2p4 = dot(p2,p4)
		kp2  = -dot(p4,q)
		kp3 = -dot(p1,q)
		p3p4 = dot(p1,p2)
		kp4  = -dot(p2,q)
	elseif(ifc.eq.'try') then
		amme=ame*1d0
		ame2=ame*ame
		mmu=ammu*1d0
		ammu2=ammu*ammu
                me2=ame2
                mm2=ammu2
		p1p2 = dot(p1,p2)
		p1p3 = dot(p1,p4)
		p1p4 = dot(p1,p3)
		kp1  = dot(p1,q)
		p2p3 = dot(p2,p4)
		p2p4 = dot(p2,p3)
		kp2  = dot(p2,q)
		kp3 = dot(p4,q)
		p3p4 = dot(p3,p4)
		kp4  = dot(p3,q)
	endif

	x12  = 2d0*ame2+2.d0*p1p2
	x13  = ame2+ammu2-2.d0*p1p3
	x14  = ame2+ammu2-2.d0*p1p4
	x15  = ame2-2.d0*kp1
	x23  = ame2+ammu2-2d0*p2p3
	x24  = ame2+ammu2-2d0*p2p4
	x25 = ame2-2d0*kp2
	x34 = 2d0*ammu2+2d0*p3p4
	x35 = ammu2+2d0*kp3
	x45 = ammu2+2d0*kp4


        mm4=mm2**2
        mm6=mm2**3
        me4=me2**2
        me6=me2**3
#ifdef COLLIER


	call D_cll(dnpi,dnpiuv,d1*ammu2,d1*x15,d1*ame2,d1*ammu2,d1*x23,
     >    d1*x34,d1*ammu2,d1*0d0,d1*ame2,d1*0d0,2)

        D0apx15epx23x34p0e0= dnpi(0,0,0,0)
	D3apx15epx23x34p0e0= dnpi(0,0,0,1)
	D00apx15epx23x34p0e0= dnpi(1,0,0,0)
	D11apx15epx23x34p0e0= dnpi(0,2,0,0)
	D13apx15epx23x34p0e0= dnpi(0,1,0,1)
	D33apx15epx23x34p0e0= dnpi(0,0,0,2)
	D22apx15epx23x34p0e0= dnpi(0,0,2,0)
	D23apx15epx23x34p0e0= dnpi(0,0,1,1)
	D2apx15epx23x34p0e0= dnpi(0,0,1,0)
	D1apx15epx23x34p0e0= dnpi(0,1,0,0)
	D12apx15epx23x34p0e0 = dnpi(0,1,1,0)

	call D_cll(dnpi,dnpiuv,d1*ammu2,d1*x15,d1*ame2,d1*ammu2,d1*x24,
     >    d1*x34,d1*ammu2,d1*0d0,d1*ame2,d1*0d0,2)

     	D0apx15epx24x34p0e0= dnpi(0,0,0,0)
	D3apx15epx24x34p0e0= dnpi(0,0,0,1)
	D00apx15epx24x34p0e0= dnpi(1,0,0,0)
	D11apx15epx24x34p0e0= dnpi(0,2,0,0)
	D13apx15epx24x34p0e0= dnpi(0,1,0,1)
	D33apx15epx24x34p0e0= dnpi(0,0,0,2)
	D2apx15epx24x34p0e0= dnpi(0,0,1,0)
	D23apx15epx24x34p0e0= dnpi(0,0,1,1)
	D22apx15epx24x34p0e0= dnpi(0,0,2,0)
	D12apx15epx24x34p0e0= dnpi(0,1,1,0)
	D1apx15epx24x34p0e0  = dnpi(0,1,0,0)

	call D_cll(dnpi,dnpiuv,d1*ammu2,d1*x25,d1*ame2,d1*ammu2,d1*x13,
     >    d1*x34,d1*ammu2,d1*0d0,d1*ame2,d1*0d0,2)

	D0apx25epx13x34p0e0= dnpi(0,0,0,0)
	D3apx25epx13x34p0e0= dnpi(0,0,0,1)
	D00apx25epx13x34p0e0= dnpi(1,0,0,0)
	D11apx25epx13x34p0e0= dnpi(0,2,0,0)
	D13apx25epx13x34p0e0= dnpi(0,1,0,1)
	D33apx25epx13x34p0e0= dnpi(0,0,0,2)
	D2apx25epx13x34p0e0= dnpi(0,0,1,0)
	D23apx25epx13x34p0e0= dnpi(0,0,1,1)
	D22apx25epx13x34p0e0= dnpi(0,0,2,0)
	D12apx25epx13x34p0e0= dnpi(0,1,1,0)
	D1apx25epx13x34p0e0= dnpi(0,1,0,0)

	call D_cll(dnpi,dnpiuv,d1*ammu2,d1*x25,d1*ame2,d1*ammu2,d1*x14,
     >    d1*x34,d1*ammu2,d1*0d0,d1*ame2,d1*0d0,2)

	D0apx25epx14x34p0e0= dnpi(0,0,0,0)
	D3apx25epx14x34p0e0= dnpi(0,0,0,1)
	D00apx25epx14x34p0e0= dnpi(1,0,0,0)
	D11apx25epx14x34p0e0= dnpi(0,2,0,0)
	D13apx25epx14x34p0e0= dnpi(0,1,0,1)
	D33apx25epx14x34p0e0= dnpi(0,0,0,2)
	D2apx25epx14x34p0e0= dnpi(0,0,1,0)
	D22apx25epx14x34p0e0= dnpi(0,0,2,0)
	D23apx25epx14x34p0e0= dnpi(0,0,1,1)
	D1apx25epx14x34p0e0= dnpi(0,1,0,0)
	D12apx25epx14x34p0e0= dnpi(0,1,1,0)

#endif
***************BOX CON RADIAZIONE REALE DA ISR CROSSATI CON BORN ISR

	include 'mumurad/Virtual_new/boxisrisr.f'

	virtual=real(elmat*vpolc(x34)*CONJG(vpolc(x34)))


***************BOX CON RADIAZIONE REALE DA ISR CROSSATI CON BORN FSR


	include 'mumurad/Virtual_new/boxisrfsr.f'

	virtual=virtual+real(elmat*vpolc(x34)*CONJG(vpolc(x12)))

	virtual=4d0*virtual*e6*alpha/(2d0*pi)
c	print*,virtual








	end
***************************************

      subroutine virtualPent_mur(p1,p2,p3,p4,q,virtual,qph,ifc)

#ifdef COLLIER
        use collier
#endif
	implicit none
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision pt1(0:3),pt2(0:3),pt3(0:3),pt4(0:3),qt(0:3)
	double precision qph(40,0:3)
	double precision virtual,virtual1
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*8 e2,e6
	double precision ame,ammu,convfac,alpha,pi,ame2,ammu2
	common/parameters/ame,ammu,convfac,alpha,pi

	double precision amme,mmu
	real*8 elmat,elmat1,elmat2,elmat3,elmat4,elmat5,elmat6
        Complex*16 elmat7,elmat8,elmat9,elmat10,elmat11
	character*2 fs
	character*3 ifc
	common/finalstate/fs

	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min

	double complex masses(0:4)
	double complex mominv(1:10)
	double complex momvec(0:3,4)


c	double precision routing(5)
	double complex routing(5)
	double complex average


        external vpolc
        double complex vpolc,alphar

        common/ialpharunning/iarun
        integer iarun

	integer i

	common/softphotonmass/phmass
	real*8 phmass

	external dot
	double precision dot


	double complex e1
	parameter (e1 = (1.d0,0.d0))

        common/instabilities/instabilities
        integer instabilities


	e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2

	if(ifc.eq.'isr') then
		amme=ame*1d0
		ame2=ame*ame
		mmu=ammu*1d0
		ammu2=ammu*ammu
		pt1=p1
		pt2=p2
		pt3=p3
		pt4=p4
		qt=q
	elseif(ifc.eq.'is1') then
		amme=ame*1d0
		ame2=ame*ame
		mmu=ammu*1d0
		ammu2=ammu*ammu
		pt1=p1
		pt2=p2
		pt3=p4
		pt4=p3
		qt=q
	elseif(ifc.eq.'fsr') then
		amme=ammu*1d0
		mmu=ame
		ame2=ammu*ammu
		ammu2=ame*ame
		pt1=-p3
		pt2=-p4
		pt3=-p1
		pt4=-p2
		qt=q

	elseif(ifc.eq.'fs1') then
		amme=ammu*1d0
		mmu=ame
		ame2=ammu*ammu
		ammu2=ame*ame
		pt1=-p3
		pt2=-p4
		pt3=-p2
		pt4=-p1
		qt=q
	 endif
	p1p2 = dot(pt1,pt2)
	p1p3 = dot(pt1,pt3)
	p1p4 = dot(pt1,pt4)
	kp1  = dot(pt1,qt)
	p2p3 = dot(pt2,pt3)
	p2p4 = dot(pt2,pt4)
	kp2  = dot(pt2,qt)
	kp3 = dot(pt3,qt)
	p3p4 = dot(pt3,pt4)
	kp4  = dot(pt4,qt)


!        instabilities=1



********************QUI PER LE TENSORIALI CON LA CONVENZIONE DI MATHEMATICA SULL'ESPRESSIONE INIZIALE

	momvec(:,1)=-e1*pt3
	momvec(:,2)=-e1*(1d0*pt3+1d0*pt4)
	momvec(:,3)=-e1*(1d0*pt1-1d0*qt)
	momvec(:,4)=-e1*pt1

	mominv(1)=e1*ammu2
	mominv(2)=e1*ammu2
	mominv(3)=e1*(ame2)
	mominv(4)=e1*0d0
	mominv(5)=e1*ame2
	mominv(6)=e1*(2d0*ammu2+2d0*p3p4)
	mominv(7)=e1*(ame2+ammu2-2d0*p2p4)
	mominv(8)=e1*(ame2-2d0*kp2)
	mominv(9)=e1*(ame2-2d0*kp1)
	mominv(10)=e1*(ame2+ammu2-2.d0*p1p3)

	masses(0)=e1*0d0
	masses(1)=e1*ammu2
	masses(2)=e1*0d0
	masses(3)=e1*ame2
	masses(4)=e1*ame2
	if(instabilities.eq.1) then
	call routing_mumu(amme,mmu,pt1,pt2,pt3,pt4,qt,
     >				momvec,mominv,masses,1,routing(1),ifc)
        endif
!

********************ROUTING 2 PER LE FULL TENSOR

	momvec(:,1)=e1*pt3
	momvec(:,2)=-e1*pt4
	momvec(:,3)=-e1*(1d0*pt1-1d0*pt3)
	momvec(:,4)=-e1*(1d0*pt4-1d0*pt2)

	mominv(1)=e1*ammu2
	mominv(2)=e1*(2d0*ammu2+2d0*p3p4)
	mominv(3)=e1*(ame2-2d0*kp2)
	mominv(4)=e1*0d0
	mominv(5)=e1*(ame2+ammu2-2d0*p2p4)
	mominv(6)=e1*ammu2
	mominv(7)=e1*ame2
	mominv(8)=e1*ame2
	mominv(9)=e1*(ame2+ammu2-2.d0*p1p3)
	mominv(10)=e1*(ame2-2d0*kp1)

	masses(0)=e1*ammu2
	masses(1)=e1*0d0
	masses(2)=e1*0d0
	masses(3)=e1*ame2
	masses(4)=e1*ame2
	if(instabilities.eq.1) then
	call routing_mumu(amme,mmu,pt1,pt2,pt3,pt4,qt,
     >				momvec,mominv,masses,2,routing(2),ifc)
	endif
********************ROUTING 3 PER LE FULL TENSOR

	momvec(:,1)=e1*(pt3+pt4)
	momvec(:,2)=e1*pt2
	momvec(:,3)=e1*pt4
	momvec(:,4)=e1*(pt2-qt)

	mominv(1)=e1*(2d0*ammu2+2d0*p3p4)
	mominv(2)=e1*(ame2-2d0*kp1)
	mominv(3)=e1*(ame2+ammu2-2d0*p2p4)
	mominv(4)=e1*(ame2+ammu2-2.d0*p1p3)
	mominv(5)=e1*(ame2-2d0*kp2)
	mominv(6)=e1*ame2
	mominv(7)=e1*ammu2
	mominv(8)=e1*0d0
	mominv(9)=e1*ammu2
	mominv(10)=e1*ame2

	masses(0)=e1*0d0
	masses(1)=e1*0d0
	masses(2)=e1*ame2
	masses(3)=e1*ammu2
	masses(4)=e1*ame2

	if(instabilities.eq.1) then
	call routing_mumu(amme,mmu,pt1,pt2,pt3,pt4,qt,
     >				momvec,mominv,masses,3,routing(3),ifc)
     	endif
********************ROUTING 4 PER LE FULL TENSOR

	momvec(:,1)=-e1*pt2
	momvec(:,2)=e1*(pt1-qt)
	momvec(:,3)=-e1*qt
	momvec(:,4)=e1*(pt4-pt2)

	mominv(1)=e1*ame2
	mominv(2)=e1*(2d0*ammu2+2d0*p3p4)
	mominv(3)=e1*ame2
	mominv(4)=e1*(ame2+ammu2-2.d0*p1p3)
	mominv(5)=e1*(ame2+ammu2-2d0*p2p4)
	mominv(6)=e1*(ame2-2d0*kp1)
	mominv(7)=e1*(ame2-2d0*kp2)
	mominv(8)=e1*ammu2
	mominv(9)=e1*0d0
	mominv(10)=e1*ammu2

	masses(0)=e1*ame2
	masses(1)=e1*0d0
	masses(2)=e1*0d0
	masses(3)=e1*ame2
	masses(4)=e1*ammu2

	if(instabilities.eq.1) then
	call routing_mumu(amme,mmu,pt1,pt2,pt3,pt4,qt,
     >				momvec,mominv,masses,4,routing(4),ifc)
     	endif
********************ROUTING 5 PER LE FULL TENSOR

	momvec(:,1)=e1*pt1
	momvec(:,2)=-e1*(pt2-qt)
	momvec(:,3)=e1*qt
	momvec(:,4)=e1*(pt1-pt3)

	mominv(1)=e1*ame2
	mominv(2)=e1*(2d0*ammu2+2d0*p3p4)
	mominv(3)=e1*ame2
	mominv(4)=e1*(ame2+ammu2-2d0*p2p4)
	mominv(5)=e1*(ame2+ammu2-2.d0*p1p3)
	mominv(6)=e1*(ame2-2d0*kp2)
	mominv(7)=e1*(ame2-2d0*kp1)
	mominv(8)=e1*ammu2
	mominv(9)=e1*0d0
	mominv(10)=e1*ammu2

	masses(0)=e1*ame2
	masses(1)=e1*0d0
	masses(2)=e1*0d0
	masses(3)=e1*ame2
	masses(4)=e1*ammu2

	call routing_mumu(amme,mmu,pt1,pt2,pt3,pt4,qt,
     >				momvec,mominv,masses,5,routing(5),ifc)

	average=routing(5)

	if(instabilities.eq.1) then
		call find_closest_two(routing,average)

     	endif
                virtual=4d0*average*e6*alpha/(2d0*pi)

!	write(16,*) routing(1)
!	write(17,*) routing(2)
!	write(18,*) routing(3)
!	write(19,*) routing(4)
!	write(20,*) routing(5)






	end
!----------------------------------------
!QUESTO È IL VIRTUALE COMPLETO
!----------------------------------------

      subroutine Virtual_MUR(pin1,pin2,p1,p2,q,virtual1,qphmat)

#ifdef COLLIER
      use collier
#endif
      implicit none
      double precision endtime,start
      double precision virtual1,virtual
      double precision pin1(0:3),pin2(0:3),p1(0:3),p2(0:3),q(0:3)
      double precision qphmat(40,0:3)
      real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
      real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45,xsum
      real*8 e2,e6,me2,me4,me6,mm2,mm4,mm6,amme,mmu
      double precision ame,ammu,convfac,alpha,pi,ame2,ammu2
      common/parameters/ame,ammu,convfac,alpha,pi

      external dot
      double precision dot
      common/lt/ltint,quadevents,ilt
      integer ltint,quadevents,ilt


      common/instabilities/instabilities
      integer instabilities


      virtual1=0d0
                   call virtualISR_mur(pin1,pin2,p1,p2,q,
     >                                       virtual,qphmat,'fsr')
		   	 virtual1=virtual

                         call virtualISR_mur(pin1,pin2,p1,p2,q,
     >                                       virtual,qphmat,'isr')

		   	virtual1=virtual1+virtual


		    call virtualMixed_mur(pin1,pin2,p1,p2,q,
     >                                       virtual,qphmat,'isr')

     			virtual1=virtual1+virtual

     		    call virtualMixed_mur(pin1,pin2,p1,p2,q,
     >                                       virtual,qphmat,'fsr')

			virtual1=virtual1+virtual


	    	    call virtualBox_mur(pin1,pin2,p1,p2,q,
     >                                       virtual,qphmat,'fsr')


     			virtual1=virtual1+virtual


	    	    call virtualBox_mur(pin1,pin2,p1,p2,q,
     >                                       virtual,qphmat,'isr')

                    virtual1=virtual1+virtual


                    if(ltint.eq.0) then

	    	    call virtualPent_mur(pin1,pin2,p1,p2,q,
     >                                       virtual,qphmat,'isr')
!


     			virtual1=virtual1+virtual

 	    	    call virtualPent_mur(pin1,pin2,p1,p2,q,
     >                                      virtual,qphmat,'is1')
                    virtual1=virtual1+virtual

                    endif

                    if(ltint.eq.1) then
                    call virtualPent_mur_quad(pin1,pin2,p1,p2,q,
     >                                       virtual,qphmat,'isr')
                          virtual1=virtual1+virtual
                    call virtualPent_mur_quad(pin1,pin2,p1,p2,q,
     >                                       virtual,qphmat,'is1')

                          virtual1=virtual1+virtual
                          ltint=0
                     endif

                          call virtualPent_mur(pin1,pin2,p1,p2,q,
     >                                       virtual,qphmat,'fsr')

      			  virtual1=virtual1+virtual


 	    	    call virtualPent_mur(pin1,pin2,p1,p2,q,
     >                                       virtual,qphmat,'fs1')

      			  virtual1=virtual1+virtual


      end
***************************************

      !======================!
      !  RADIATIVE PIONS-NLO !
      !======================!
      subroutine virtualISR_pr(p1,p2,p3,p4,q,virtual,qph)

c IN QUESTA ROUTINE CALCOLO IL CONTRIBUTO DI INTERFERENZA TRA IL VIRTUALE ISR E IL BORN, NELLA PRIMA PARTE CALCOLO L'INTERFERENZA CON ISR DEL BORN, MENTRE NELLA SECONDA CALCOLO QUELLA
C CON IL BORN FSR, POI PER SIMMETRIA PORTO TUTTO QUESTO AL SIMMETRICO, OSSIA VIRTUALE FSR-BORN FSR E VIRTUALE FSR BORN ISR
#ifdef COLLIER
        use collier
#endif
	implicit double complex (z)
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision qph(40,0:3)
	double precision virtual,virtual1
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
	real*8 e2,e6,me2,me4,me6,mp2,mp4,mp6
	double precision ame,ammu,convfac,alpha,pi,ame2,ampi2,ampi
	common/parameters/ame,ammu,convfac,alpha,pi
        common/hparam/ampi
	Complex*16 elmat,elmat1,elmat2,elmat3,elmat4,elmat5,elmat6
        Complex*16 elmat7,elmat8,elmat9,elmat10,elmat11
	character*2 fs
	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min
	common/softphotonmass/phmass
	real*8 phmass

	external dot
	real*8 dot

	common/radoption/ifisrfsr
        integer ifisrfsr


        double complex vpolc,alphar

        common/ialpharunning/iarun
        integer iarun

**********FF

	integer iffpi
	character*10 what_ffpi
	common/ffpioption/iffpi,what_ffpi
	double complex ffpiis
	double precision ffpiis2
	double complex ffpifs
	double precision ffpifs2


*********CLL STUFF
	double complex cnpi(0:1,0:2,0:2)
	double complex cnpiuv(0:1,0:2,0:2)
	double complex cnpi2(0:0,0:1,0:1)
	double complex cnpiuv2(0:0,0:1,0:1)
	double complex dnpi(0:1,0:3,0:3,0:3)
	double complex dnpiuv(0:1,0:3,0:3,0:3)
	double complex d1
	parameter (d1 = (1.d0,0.d0))
	double complex c1
	parameter (c1 = (1.d0,0.d0))
	double complex b1
	parameter (b1 = (1.d0,0.d0))
	double complex a1
	parameter (a1 = (1.d0,0.d0))

****BUBBLES AND TADPOLE

	Complex*16 A0ae
	Complex*16 B0ax150e
	Complex*16 B0ax250e
	Complex*16 B0a0ee
	Complex*16 B0ae0e

***TRIANGLES

	Complex*16 C00aex150e0e
	Complex*16 C0aex150e0e
	Complex*16 C11aex150e0e
	Complex*16 C1aex150e0e
	Complex*16 C2aex150e0e
	Complex*16 C0ax150e0ee

	Complex*16 C00aex250e0e
	Complex*16 C0aex250e0e
	Complex*16 C11aex250e0e
	Complex*16 C1aex250e0e
	Complex*16 C2aex250e0e
	Complex*16 C0ax250e0ee

	Complex*16 C00aex15x34e0e
	Complex*16 C0aex15x34e0e
	Complex*16 C11aex15x34e0e
	Complex*16 C1aex15x34e0e
	Complex*16 C12aex15x34e0e
	Complex*16 C22aex15x34e0e
	Complex*16 C2aex15x34e0e
	Complex*16 C0ax15x34e0ee
	Complex*16 C1ax15x34e0ee


	Complex*16 C00aex25x34e0e
	Complex*16 C0aex25x34e0e
	Complex*16 C11aex25x34e0e
	Complex*16 C1aex25x34e0e
	Complex*16 C12aex25x34e0e
	Complex*16 C22aex25x34e0e
	Complex*16 C2aex25x34e0e
	Complex*16 C0ax25x34e0ee
	Complex*16 C1ax25x34e0ee


	Complex*16 C0aex12e0ee
	Complex*16 C1aex12e0ee
	Complex*16 C0aeex12e0e
	Complex*16 C2aeex12e0e

	Complex*16 C1ax12x340eee
	Complex*16 C0a0x12x34eee

************BOXES
	Complex*16 D00aeex340x12x15e0ee
	Complex*16 D0aeex340x12x15e0ee
	Complex*16 D001aeex340x12x15e0ee
	Complex*16 D111aeex340x12x15e0ee
	Complex*16 D11aeex340x12x15e0ee
	Complex*16 D1aeex340x12x15e0ee
	Complex*16 D002aeex340x12x15e0ee
	Complex*16 D112aeex340x12x15e0ee
	Complex*16 D12aeex340x12x15e0ee
	Complex*16 D122aeex340x12x15e0ee
	Complex*16 D22aeex340x12x15e0ee
	Complex*16 D2aeex340x12x15e0ee
	Complex*16 D003aeex340x12x15e0ee
	Complex*16 D113aeex340x12x15e0ee
	Complex*16 D13aeex340x12x15e0ee
	Complex*16 D123aeex340x12x15e0ee
	Complex*16 D223aeex340x12x15e0ee
	Complex*16 D23aeex340x12x15e0ee
	Complex*16 D133aeex340x12x15e0ee
	Complex*16 D233aeex340x12x15e0ee
	Complex*16 D3aeex340x12x15e0ee

	Complex*16 D00aeex340x12x25e0ee
	Complex*16 D0aeex340x12x25e0ee
	Complex*16 D001aeex340x12x25e0ee
	Complex*16 D111aeex340x12x25e0ee
	Complex*16 D11aeex340x12x25e0ee
	Complex*16 D1aeex340x12x25e0ee
	Complex*16 D002aeex340x12x25e0ee
	Complex*16 D112aeex340x12x25e0ee
	Complex*16 D12aeex340x12x25e0ee
	Complex*16 D122aeex340x12x25e0ee
	Complex*16 D22aeex340x12x25e0ee
	Complex*16 D2aeex340x12x25e0ee
	Complex*16 D003aeex340x12x25e0ee
	Complex*16 D113aeex340x12x25e0ee
	Complex*16 D13aeex340x12x25e0ee
	Complex*16 D123aeex340x12x25e0ee
	Complex*16 D223aeex340x12x25e0ee
	Complex*16 D23aeex340x12x25e0ee
	Complex*16 D133aeex340x12x25e0ee
	Complex*16 D233aeex340x12x25e0ee
	Complex*16 D3aeex340x12x25e0ee



**********

	e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2
	ame2=ame*ame
	ampi2=ampi*ampi

        me2=ame2
        mp2=ampi2
        me4=me2**2
        me6=me2**3
        mp4=mp2**2
        mp6=mp2**3

	p1p2 = dot(p1,p2)
	p1p3 = dot(p1,p3)
	p1p4 = dot(p1,p4)
	kp1  = dot(p1,q)
	p2p3 = dot(p2,p3)
	p2p4 = dot(p2,p4)
	kp2  = dot(p2,q)
	kp3 = dot(p3,q)
	p3p4 = dot(p3,p4)
	kp4  = dot(p4,q)



	x12  = 2d0*ame2+2.d0*p1p2
	x13  = ame2+ampi2-2.d0*p1p3
	x14  = ame2+ampi2-2.d0*p1p4
	x15  = ame2-2.d0*kp1
	x23  = ame2+ampi2-2d0*p2p3
	x24  = ame2+ampi2-2d0*p2p4
	x25 = ame2-2d0*kp2
	x34 = 2d0*ampi2+2d0*p3p4
	x35 = ampi2+2d0*kp3
	x45 = ampi2+2d0*kp4

	ffpifs=(1.d0,0.d0)
	ffpifs2=1.d0

	ffpiis=(1.d0,0.d0)
	ffpiis2=1.d0
	if(iffpi.gt.0) then
		call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2)
		call ffpi(2.d0*(ampi**2+p3p4),ffpiis,ffpiis2)

	endif

******ALPHA RUNNING
	ffpiis=ffpiis*vpolc(2.d0*(ampi**2+p3p4))
	ffpiis2=ffpiis*CONJG(ffpiis)

	ffpifs=ffpifs*vpolc(2.d0*(ame**2+p1p2))
	ffpifs2=ffpifs*CONJG(ffpifs)
****************

#ifdef COLLIER
!     call SetDeltaUV_cll(1000000d0)

	call A0_cll(A0ae,a1*ame2)

        call B0_cll(B0ax250e,b1*x25,b1*0d0,b1*ame2)

	call B0_cll(B0ax150e,b1*x15,b1*0d0,b1*ame2)

	call B0_cll(B0ae0e,b1*ame2,b1*0d0,b1*ame2)

	call B0_cll(B0a0ee,b1*0d0,b1*ame2,b1*ame2)




       call C_cll(cnpi,cnpiuv,c1*ame2,c1*x15,c1*0d0,c1*ame2,
     >     c1*0d0,c1*ame2,2)
     		C00aex150e0e=cnpi(1,0,0)
		C0aex150e0e=cnpi(0,0,0)
		C11aex150e0e=cnpi(0,2,0)
		C1aex150e0e=cnpi(0,1,0)
		C2aex150e0e=cnpi(0,0,1)
		C0ax150e0ee=cnpi(0,0,0)

       call C_cll(cnpi,cnpiuv,c1*ame2,c1*x25,c1*0d0,c1*ame2,
     >     c1*0d0,c1*ame2,2)
     		C00aex250e0e=cnpi(1,0,0)
		C0aex250e0e=cnpi(0,0,0)
		C11aex250e0e=cnpi(0,2,0)
		C1aex250e0e=cnpi(0,1,0)
		C2aex250e0e=cnpi(0,0,1)
		C0ax250e0ee=cnpi(0,0,0)

       call C_cll(cnpi,cnpiuv,c1*ame2,c1*x25,c1*x34,c1*ame2,
     >       c1*0d0,c1*ame2,2)
	     	C00aex25x34e0e=cnpi(1,0,0)
		C0aex25x34e0e=cnpi(0,0,0)
		C11aex25x34e0e=cnpi(0,2,0)
		C1aex25x34e0e=cnpi(0,1,0)
		C12aex25x34e0e=cnpi(0,1,1)
		C22aex25x34e0e=cnpi(0,0,2)
		C2aex25x34e0e=cnpi(0,0,1)
		C0ax25x34e0ee=cnpi(0,0,0)
		C1ax25x34e0ee=cnpi(0,0,1)

       call C_cll(cnpi,cnpiuv,c1*ame2,c1*x15,c1*x34,c1*ame2,
     >       c1*0d0,c1*ame2,2)
	     	C00aex15x34e0e=cnpi(1,0,0)
		C0aex15x34e0e=cnpi(0,0,0)
		C11aex15x34e0e=cnpi(0,2,0)
		C1aex15x34e0e=cnpi(0,1,0)
		C12aex15x34e0e=cnpi(0,1,1)
		C22aex15x34e0e=cnpi(0,0,2)
		C2aex15x34e0e=cnpi(0,0,1)
		C0ax15x34e0ee=cnpi(0,0,0)
		C1ax15x34e0ee=cnpi(0,0,1)


       call C_cll(cnpi2,cnpiuv2,c1*x12,c1*x34,c1*0d0,c1*ame2,
     >       c1*ame2,c1*ame2,1)

		C1ax12x340eee=cnpi2(0,1,0)
		C0a0x12x34eee=cnpi2(0,0,0)

       call C_cll(cnpi2,cnpiuv2,c1*ame2,c1*x12,c1*ame2,c1*0d0,
     >       c1*ame2,c1*ame2,1)

		C0aex12e0ee=cnpi2(0,0,0)
		C1aex12e0ee=cnpi2(0,1,0)
		C0aeex12e0e=cnpi2(0,0,0)
		C2aeex12e0e=cnpi2(0,1,0)


	call D_cll(dnpi,dnpiuv,d1*ame2,d1*ame2,d1*x34,d1*0d0,d1*x12,
     >    d1*x25,d1*ame2,d1*0d0,d1*ame2,d1*ame2,3)
		D00aeex340x12x25e0ee = dnpi(1,0,0,0)
		D0aeex340x12x25e0ee  = dnpi(0,0,0,0)
		D001aeex340x12x25e0ee= dnpi(1,1,0,0)
		D111aeex340x12x25e0ee= dnpi(0,3,0,0)
		D11aeex340x12x25e0ee = dnpi(0,2,0,0)
		D1aeex340x12x25e0ee  = dnpi(0,1,0,0)
		D002aeex340x12x25e0ee= dnpi(1,0,1,0)
		D112aeex340x12x25e0ee= dnpi(0,2,1,0)
		D12aeex340x12x25e0ee = dnpi(0,1,1,0)
		D122aeex340x12x25e0ee= dnpi(0,1,2,0)
		D22aeex340x12x25e0ee = dnpi(0,0,2,0)
		D2aeex340x12x25e0ee  = dnpi(0,0,1,0)
		D003aeex340x12x25e0ee= dnpi(1,0,0,1)
		D113aeex340x12x25e0ee= dnpi(0,2,0,1)
		D13aeex340x12x25e0ee = dnpi(0,1,0,1)
		D123aeex340x12x25e0ee= dnpi(0,1,1,1)
		D223aeex340x12x25e0ee= dnpi(0,0,2,1)
		D23aeex340x12x25e0ee = dnpi(0,0,1,1)
		D133aeex340x12x25e0ee= dnpi(0,1,0,2)
		D233aeex340x12x25e0ee= dnpi(0,0,1,2)
		D3aeex340x12x25e0ee  = dnpi(0,0,0,1)


	call D_cll(dnpi,dnpiuv,d1*ame2,d1*ame2,d1*x34,d1*0d0,d1*x12,
     >    d1*x15,d1*ame2,d1*0d0,d1*ame2,d1*ame2,3)
		D00aeex340x12x15e0ee = dnpi(1,0,0,0)
		D0aeex340x12x15e0ee  = dnpi(0,0,0,0)
		D001aeex340x12x15e0ee= dnpi(1,1,0,0)
		D111aeex340x12x15e0ee= dnpi(0,3,0,0)
		D11aeex340x12x15e0ee = dnpi(0,2,0,0)
		D1aeex340x12x15e0ee  = dnpi(0,1,0,0)
		D002aeex340x12x15e0ee= dnpi(1,0,1,0)
		D112aeex340x12x15e0ee= dnpi(0,2,1,0)
		D12aeex340x12x15e0ee = dnpi(0,1,1,0)
		D122aeex340x12x15e0ee= dnpi(0,1,2,0)
		D22aeex340x12x15e0ee = dnpi(0,0,2,0)
		D2aeex340x12x15e0ee  = dnpi(0,0,1,0)
		D003aeex340x12x15e0ee= dnpi(1,0,0,1)
		D113aeex340x12x15e0ee= dnpi(0,2,0,1)
		D13aeex340x12x15e0ee = dnpi(0,1,0,1)
		D123aeex340x12x15e0ee= dnpi(0,1,1,1)
		D223aeex340x12x15e0ee= dnpi(0,0,2,1)
		D23aeex340x12x15e0ee = dnpi(0,0,1,1)
		D133aeex340x12x15e0ee= dnpi(0,1,0,2)
		D233aeex340x12x15e0ee= dnpi(0,0,1,2)
		D3aeex340x12x15e0ee  = dnpi(0,0,0,1)

#endif

****************(VIRTUALE ISR REALE ISR) CROSSATO CON BORN ISR

                include 'pipirad/Virtual_new/isrisr.f'

                virtual=real(elmat)*ffpiis2

****************(VIRTUALE ISR REALE ISR) CROSSATO CON BORN FSR

                include 'pipirad/Virtual_new/isrfsr.f'
**************************
	if(ifisrfsr.eq.11) then

           virtual=virtual+real((elmat)*ffpiis*CONJG(ffpifs))

	elseif(ifisrfsr.eq.10) then

	  virtual=virtual

	endif


	virtual=4d0*virtual*e6*alpha/(2d0*pi)



      end
*****************************************************

      subroutine virtualMixedISR_pr(p1,p2,p3,p4,q,virtual,qph)
#ifdef COLLIER
        use collier
#endif
	implicit double complex (z)
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision qph(40,0:3)
	double precision virtual,virtual1
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
	real*8 e2,e6,me2,me4,me6,mp2,mp4,mp6
	double precision ame,ammu,convfac,alpha,pi,ame2,ampi2,ampi
	common/parameters/ame,ammu,convfac,alpha,pi
        common/hparam/ampi
	Complex*16 elmat1,elmat2,elmat3,elmat
	character*2 fs
	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min
	common/softphotonmass/phmass
	real*8 phmass

	external dot
	real*8 dot

        double complex vpolc,alphar

        common/ialpharunning/iarun
        integer iarun

        common/indexfsrdisp/iFSRdisp
        integer iFSRdisp
**********FF

	integer iffpi
	character*10 what_ffpi
	common/ffpioption/iffpi,what_ffpi
	double complex ffpiis
	double precision ffpiis2
	double complex ffpifs
	double precision ffpifs2


*********CLL STUFF
	double complex cnpi(0:1,0:2,0:2)
	double complex cnpiuv(0:1,0:2,0:2)

	double complex c1
	parameter (c1 = (1.d0,0.d0))
	double complex b1
	parameter (b1 = (1.d0,0.d0))



****BUBBLES AND TRIANGLES

	Complex*16 B0ae0e
	Complex*16 B0ax12ee

	Complex*16 C00aex12e0ee
	Complex*16 C0aex12e0ee
	Complex*16 C11aex12e0ee
	Complex*16 C1aex12e0ee
	Complex*16 C12aex12e0ee

****************

	e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2
	ame2=ame*ame
	ampi2=ampi*ampi
        me2=ame2
        mp2=ampi2
        me4=me2**2
        me6=me2**3
        mp4=mp2**2
        mp6=mp2**3

	p1p2 = dot(p1,p2)
	p1p3 = dot(p1,p3)
	p1p4 = dot(p1,p4)
	kp1  = dot(p1,q)
	p2p3 = dot(p2,p3)
	p2p4 = dot(p2,p4)
	kp2  = dot(p2,q)
	kp3 = dot(p3,q)
	p3p4 = dot(p3,p4)
	kp4  = dot(p4,q)

	x12  = 2d0*ame2+2.d0*p1p2
	x13  = ame2+ampi2-2.d0*p1p3
	x14  = ame2+ampi2-2.d0*p1p4
	x15  = ame2-2.d0*kp1
	x23  = ame2+ampi2-2d0*p2p3
	x24  = ame2+ampi2-2d0*p2p4
	x25 = ame2-2d0*kp2
	x34 = 2d0*ampi2+2d0*p3p4
	x35 = ampi2+2d0*kp3
	x45 = ampi2+2d0*kp4

**************
	ffpifs=(1.d0,0.d0)
	ffpifs2=1.d0

	ffpiis=(1.d0,0.d0)
	ffpiis2=1.d0
	if(iffpi.gt.0) then
		call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2)
		call ffpi(2.d0*(ampi**2+p3p4),ffpiis,ffpiis2)

	endif

******ALPHA RUNNING
	ffpiis=ffpiis*vpolc(2.d0*(ampi**2+p3p4))
	ffpiis2=ffpiis*CONJG(ffpiis)

	ffpifs=ffpifs*vpolc(2.d0*(ame**2+p1p2))
	ffpifs2=ffpifs*CONJG(ffpifs)
****************

#ifdef COLLIER

!	call SetDeltaUV_cll(1000000d0)
        call B0_cll(B0ae0e,b1*ame2,b1*0d0,b1*ame2)

        call B0_cll(B0ax12ee,b1*x12,b1*ame2,b1*ame2)

        call C_cll(cnpi,cnpiuv,c1*ame2,c1*x12,c1*ame2,
     >       c1*0d0,c1*ame2,c1*ame2,2)

	C00aex12e0ee=cnpi(1,0,0)
	C0aex12e0ee=cnpi(0,0,0)
	C11aex12e0ee=cnpi(0,2,0)
	C1aex12e0ee=cnpi(0,1,0)
	C12aex12e0ee =cnpi(0,1,1)

#endif

****************(VIRTUALE ISR REALE FSR) CROSSATO CON BORN ISR

	include 'pipirad/Virtual_new/mixedisrfsr_isr.f'

	virtual=real(elmat*(CONJG(ffpiis)*ffpifs))
****************(VIRTUALE ISR REALE FSR) CROSSATO CON BORN FSR

	include 'pipirad/Virtual_new/mixedisrfsr_fsr.f'
**************************

	virtual=virtual+real(elmat)*ffpifs2


	virtual=4d0*virtual*e6*alpha/(2d0*pi)
	end

      subroutine virtualFSR_pr(p1,p2,p3,p4,q,virtual,qph)

#ifdef COLLIER
        use collier
#endif
	implicit double complex (z)
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision qph(40,0:3)
	double precision virtual,virtual1
        Complex*32 gvmdq,gvmdq1,gvmdq2,mq,derampl
        complex*32 gvmdq0,deltafsr,dfsr2
        complex*16 gvmd
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
	real*8 e2,e6,me2,me4,me6,mp2,mp4,mp6
	double precision ame,ammu,convfac,alpha,pi,ame2,ampi2,ampi
	common/parameters/ame,ammu,convfac,alpha,pi
        common/hparam/ampi
	Complex*16 elmat,elmat2,elmat3,elmat4,elmat5,elmat6
        Complex*16 elmat7,elmat8,elmat9,elmat10,elmat11,elmat12
        Complex*16 elmat13,elmat14,elmat15,elmat16,elmat17,elmat18
        Complex*16 elmat19,elmat20,elmat21,elmat22,elmat23
	character*2 fs
	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min
	common/softphotonmass/phmass
	real*8 phmass

	external dot
	real*8 dot

	common/radoption/ifisrfsr
        integer ifisrfsr

        double complex vpolc,alphar

        common/ialpharunning/iarun
        integer iarun
        common/indexfsrdisp/iFSRdisp
        integer iFSRdisp
**************CHECK FULVIO
        complex*32 c1Bd,c2Bd,c3Bd,c4Bd,c5Bd
        common/borncouplings/c1Bd,c2Bd,c3Bd,c4Bd,c5Bd
        complex*32 e3
**********FF

	integer iffpi
	character*10 what_ffpi
	common/ffpioption/iffpi,what_ffpi
	double complex ffpiis
	double precision ffpiis2
	double complex ffpifs
	double precision ffpifs2

**********GVMD
        integer i,j,k
        double complex ct,im,cij
        double precision m_n(1:9),g_n(1:9)
        double complex m1,m2,m3
        double complex a_n(1:9),lv2(1:9),lv2diff
        common/par_gvmd/m_n,g_n,a_n

        common/resonancesgvmd/nres
        integer nres
        real*16 delta,delta1,math,form
        
        common/gvmdchecks/iscoll,fscoll
        integer iscoll,fscoll
********FsQED
        double precision spmin,spmax
        double complex sp,ss
        double precision wsp,wss
        double complex dispsp,dispss
        double complex ffpisp,ffpiss
        complex*32 dispq,dispqct
        complex*32 dispq1,dispq1ct
        double complex dispspq,dispssq
        double precision rdiff
        common/compotontensor/compttens
        integer compttens

*********CLL STUFF
	double complex cnpi(0:1,0:2,0:2)
	double complex cnpiuv(0:1,0:2,0:2)
	double complex dnpi(0:1,0:2,0:2,0:2)
	double complex dnpiuv(0:1,0:2,0:2,0:2)

	double complex d1
	parameter (d1 = (1.d0,0.d0))
	double complex c1
	parameter (c1 = (1.d0,0.d0))
	double complex b1
	parameter (b1 = (1.d0,0.d0))
	double complex a1
	parameter (a1 = (1.d0,0.d0))

****BUBBLES AND TADPOLE
	Complex*16 A0ap
	Complex*16 B0ap0p
	Complex*16 B0ax350p
	Complex*16 B0ax450p
	Complex*16 B0ax34pp

********TRIANGLES

	Complex*16 C00apx350p0p
	Complex*16 C11apx350p0p
	Complex*16 C1apx350p0p
	Complex*16 C12apx350p0p

	Complex*16 C00apx450p0p
	Complex*16 C11apx450p0p
	Complex*16 C1apx450p0p
	Complex*16 C12apx450p0p

	Complex*16 C00apx35x12p0p
	Complex*16 C0apx35x12p0p
	Complex*16 C11apx35x12p0p
	Complex*16 C1apx35x12p0p
	Complex*16 C12apx35x12p0p
	Complex*16 C22apx35x12p0p
	Complex*16 C2apx35x12p0p

	Complex*16 C00apx45x12p0p
	Complex*16 C0apx45x12p0p
	Complex*16 C11apx45x12p0p
	Complex*16 C1apx45x12p0p
	Complex*16 C12apx45x12p0p
	Complex*16 C22apx45x12p0p
	Complex*16 C2apx45x12p0p

	Complex*16 C0appx34p0p

	Complex*16 C11ax34x120ppp
	Complex*16 C1ax34x120ppp
	Complex*16 C00a0x34x12ppp
	Complex*16 C12a0x34x12ppp

**********BOXES

	Complex*16 D00appx120x34x35p0pp
	Complex*16 D11appx120x34x35p0pp
	Complex*16 D1appx120x34x35p0pp
	Complex*16 D12appx120x34x35p0pp
	Complex*16 D22appx120x34x35p0pp
	Complex*16 D2appx120x34x35p0pp
	Complex*16 D13appx120x34x35p0pp
	Complex*16 D23appx120x34x35p0pp

	Complex*16 D00appx120x34x45p0pp
	Complex*16 D11appx120x34x45p0pp
	Complex*16 D1appx120x34x45p0pp
	Complex*16 D12appx120x34x45p0pp
	Complex*16 D22appx120x34x45p0pp
	Complex*16 D2appx120x34x45p0pp
	Complex*16 D13appx120x34x45p0pp
	Complex*16 D23appx120x34x45p0pp
**********************

	e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2
        
	ame2=ame*ame
	ampi2=ampi*ampi
        me2=ame2
        mp2=ampi2
        me4=me2**2
        me6=me2**3
        mp4=mp2**2
        mp6=mp2**3

	p1p2 = dot(p1,p2)
	p1p3 = dot(p1,p3)
	p1p4 = dot(p1,p4)
	kp1  = dot(p1,q)
	p2p3 = dot(p2,p3)
	p2p4 = dot(p2,p4)
	kp2  = dot(p2,q)
	kp3 = dot(p3,q)
	p3p4 = dot(p3,p4)
	kp4  = dot(p4,q)


	x12  = 2d0*ame2+2.d0*p1p2
	x13  = ame2+ampi2-2.d0*p1p3
	x14  = ame2+ampi2-2.d0*p1p4
	x15  = ame2-2.d0*kp1
	x23  = ame2+ampi2-2d0*p2p3
	x24  = ame2+ampi2-2d0*p2p4
	x25 = ame2-2d0*kp2
	x34 = 2d0*ampi2+2d0*p3p4
	x35 = ampi2+2d0*kp3
	x45 = ampi2+2d0*kp4


**************
	ffpifs=(1.d0,0.d0)
	ffpifs2=1.d0

	ffpiis=(1.d0,0.d0)
	ffpiis2=1.d0
	if(iffpi.gt.0) then
		call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2)
		call ffpi(2.d0*(ampi**2+p3p4),ffpiis,ffpiis2)

	endif

******ALPHA RUNNING
	ffpiis=ffpiis*vpolc(2.d0*(ampi**2+p3p4))
	ffpiis2=ffpiis*CONJG(ffpiis)

	ffpifs=ffpifs*vpolc(2.d0*(ame**2+p1p2))
	ffpifs2=ffpifs*CONJG(ffpifs)
****************

#ifdef COLLIER
	call GetMuUV2_cll(mu2dim)
        
        call GetDeltaUV_cll(deltauv)
        
!       call SetDeltaUV_cll(1000000d0)
	call A0_cll(A0ap,a1*ampi2)

        call B0_cll(B0ax350p,b1*x35,b1*0d0,b1*ampi2)

	call B0_cll(B0ax450p,b1*x45,b1*0d0,b1*ampi2)

	call B0_cll(B0ap0p,b1*ampi2,b1*0d0,b1*ampi2)

	call B0_cll(B0ax34pp,b1*x34,b1*ampi2,b1*ampi2)



       call C_cll(cnpi,cnpiuv,c1*ampi2,c1*x35,c1*0d0,c1*ampi2,
     >     c1*0d0,c1*ampi2,2)

		C00apx350p0p=cnpi(1,0,0)
		C11apx350p0p=cnpi(0,2,0)
		C1apx350p0p=cnpi(0,1,0)
		C12apx350p0p=cnpi(0,1,1)

       call C_cll(cnpi,cnpiuv,c1*ampi2,c1*x45,c1*0d0,c1*ampi2,
     >     c1*0d0,c1*ampi2,2)

		C00apx450p0p=cnpi(1,0,0)
		C11apx450p0p=cnpi(0,2,0)
		C1apx450p0p=cnpi(0,1,0)
		C12apx450p0p=cnpi(0,1,1)



       call C_cll(cnpi,cnpiuv,c1*ampi2,c1*x35,c1*x12,c1*ampi2,
     >       c1*0d0,c1*ampi2,2)
		C00apx35x12p0p=cnpi(1,0,0)
		C0apx35x12p0p=cnpi(0,0,0)
		C11apx35x12p0p=cnpi(0,2,0)
		C1apx35x12p0p=cnpi(0,1,0)
		C12apx35x12p0p=cnpi(0,1,1)
		C22apx35x12p0p=cnpi(0,0,2)
		C2apx35x12p0p=cnpi(0,0,1)

       call C_cll(cnpi,cnpiuv,c1*ampi2,c1*x45,c1*x12,c1*ampi2,
     >       c1*0d0,c1*ampi2,2)
		C00apx45x12p0p=cnpi(1,0,0)
		C0apx45x12p0p=cnpi(0,0,0)
		C11apx45x12p0p=cnpi(0,2,0)
		C1apx45x12p0p=cnpi(0,1,0)
		C12apx45x12p0p=cnpi(0,1,1)
		C22apx45x12p0p=cnpi(0,0,2)
		C2apx45x12p0p=cnpi(0,0,1)

       call C0_cll(C0appx34p0p,c1*ampi2,c1*ampi2,c1*x34,c1*ampi2,
     >		   c1*0d0,c1*ampi2)

       call C_cll(cnpi,cnpiuv,c1*x34,c1*x12,c1*0d0,c1*ampi2,
     >       c1*ampi2,c1*ampi2,2)
		C11ax34x120ppp=cnpi(0,2,0)
		C1ax34x120ppp =cnpi(0,1,0)
       call C_cll(cnpi,cnpiuv,c1*0d0,c1*x34,c1*x12,c1*ampi2,
     >       c1*ampi2,c1*ampi2,2)
		C00a0x34x12ppp=cnpi(1,0,0)
		C12a0x34x12ppp=cnpi(0,1,1)


	call D_cll(dnpi,dnpiuv,d1*ampi2,d1*ampi2,d1*x12,d1*0d0,d1*x34,
     >    d1*x35,d1*ampi2,d1*0d0,d1*ampi2,d1*ampi2,2)

		D00appx120x34x35p0pp = dnpi(1,0,0,0)
		D11appx120x34x35p0pp = dnpi(0,2,0,0)
		D1appx120x34x35p0pp  = dnpi(0,1,0,0)
		D12appx120x34x35p0pp = dnpi(0,1,1,0)
		D22appx120x34x35p0pp = dnpi(0,0,2,0)
		D2appx120x34x35p0pp  = dnpi(0,0,1,0)
		D13appx120x34x35p0pp = dnpi(0,1,0,1)
		D23appx120x34x35p0pp = dnpi(0,0,1,1)



	call D_cll(dnpi,dnpiuv,d1*ampi2,d1*ampi2,d1*x12,d1*0d0,d1*x34,
     >    d1*x45,d1*ampi2,d1*0d0,d1*ampi2,d1*ampi2,2)

		D00appx120x34x45p0pp = dnpi(1,0,0,0)
		D11appx120x34x45p0pp = dnpi(0,2,0,0)
		D1appx120x34x45p0pp  = dnpi(0,1,0,0)
		D12appx120x34x45p0pp = dnpi(0,1,1,0)
		D22appx120x34x45p0pp = dnpi(0,0,2,0)
		D2appx120x34x45p0pp  = dnpi(0,0,1,0)
		D13appx120x34x45p0pp = dnpi(0,1,0,1)
		D23appx120x34x45p0pp = dnpi(0,0,1,1)

#endif

****************(VIRTUALE FSR REALE FSR) CROSSATO CON BORN ISR

	include 'pipirad/Virtual_new/purefsrisr.f'
**************************
       if(ifisrfsr.eq.11) then

        virtual=real(elmat*CONJG(ffpiis)*ffpifs)

      elseif(ifisrfsr.eq.01) then
        virtual=0d0

      endif


****************(VIRTUALE FSR REALE FSR) CROSSATO CON BORN FSR

	include 'pipirad/Virtual_new/purefsrfsr.f'
**************************

!
 	virtual=virtual+real(elmat)*ffpifs2


*---------------------------------------------
*------------------GVMD-----------------------
*---------------------------------------------
        if(iFSRdisp.eq.1) then
           if(iffpi.eq.2) then
!!!!! GIÀ DEFINITI COSÌ NEL FF,OSSIA A_N/CT
              im=(0.d0,1.d0)           
              do i=1,nres

                 lv2(i) = m_n(i)**2 - im*m_n(i)*g_n(i)
              enddo
!!!!!!!!!NOW CALCULATION AS I IMPLEMENTED IT IN THE LATEX
              do i=1,nres
                 cij=-a_n(i)
                 do k=1,nres

                    if(k.ne.i) then
                       cij = cij + 2.d0*a_n(k)*lv2(k)/(lv2(i)-lv2(k))
                    endif
                 enddo
                 delta=1.q-10
                 
                 if(fscoll.eq.1) then
                    delta=1.q-5
                 endif
                 
                 mq=lv2(i)*1.q0
                 
!!!!  CT li aggiungo qui in quad_precision, 0 stands for value,1 for the analytical derivative
                 
                 call ct_purefsr_quad(p1,p2,p3,p4,q,mq,gvmdq,qph,0)
                 virtual=virtual+real(a_n(i)*cij*gvmdq)
                 
                 call purefsr_quad(p1,p2,p3,p4,q,mq,gvmdq,qph)
                 virtual=virtual+real(a_n(i)*cij*gvmdq)
!!!!!!!!!!!!!TO FIND DERIVATIVE
                 
                 call purefsr_quad(p1,p2,p3,p4,q,mq*(1.q0+1.q0*delta)
     >                ,gvmdq1,qph)

                 derampl=(gvmdq1-gvmdq)/(delta*mq)
                 
                 if(fscoll.eq.1) then
                    call purefsr_quad(p1,p2,p3,p4,q,mq*(1.q0-1.q0*delta)
     >                   ,gvmdq,qph)
                    derampl=(gvmdq1-gvmdq)/(2.q0*delta*mq)
                 endif
                 virtual=virtual+real(a_n(i)**2*lv2(i)*derampl)

!!!   LA DERIVATA DEI CT è fatta qui
                 
                 call ct_purefsr_quad(p1,p2,p3,p4,q,mq,derampl,qph,1)

                 virtual=virtual+real(a_n(i)**2*lv2(i)*derampl)

              enddo
*---------------------------------------------
*------------------FsQED----------------------
*---------------------------------------------
              
           elseif(iffpi.eq.3) then
              if(compttens.eq.3) then
                 virtual=(0.d0,0.d0)
              endif

!--------QUI GENERO S' TRA 4*MP^2 E 9 GeV^2 (ASSUMIAMO VALIDO IL FF FINO A 3 GEV)------           
              
              spmin = 4.d0*mp2
              spmax = 9.d0
              call sp_disp(spmin,spmax,sp,wsp)
              call ffpi(real(sp),ffpisp,ffpifs2)
              
              call purefsr_disp(p1,p2,p3,p4,q,sp,dispsp,qph)

              virtual=virtual-2.d0/pi*real(aimag(ffpisp)/sp*dispsp/wsp)
!--------QUI IL DOPPIO DISPERSIVO
              
              call sp_disp(spmin,spmax,ss,wss)
              call ffpi(real(ss),ffpiss,ffpifs2)
              
              call purefsr_disp(p1,p2,p3,p4,q,ss,dispss,qph)
              
!-------Se |s'-s''|/s'<10^-5 il rapporto in doppia potrebbe diventare instabile--------!
!-------quindi switcho in 4-pla precisione per calcolare il ratio----------------------!
              
              if(abs((sp-ss)/sp).ge.1d-5) then

                 virtual=virtual+1.d0/(pi**2)*
     >                real(aimag(ffpiss)/ss/wss
     >                *aimag(ffpisp)/sp/wsp*
     >                1.d0/(ss-sp)*(dispss*ss-dispsp*sp))

              else
!Amplitude at s'
                 call ct_purefsr_quad(p1,p2,p3,p4,q,sp*1.q0,dispqct,
     >                qph,0)
                 call purefsr_quad(p1,p2,p3,p4,q,sp*1.q0,dispq,qph)
                 dispspq=dispq+dispqct
!Amplitude at s''
                 call ct_purefsr_quad(p1,p2,p3,p4,q,ss*1.q0,dispq1ct,
     >                qph,0)
                 call purefsr_quad(p1,p2,p3,p4,q,ss*1.q0,dispq1,qph)
                 dispssq=dispq1+dispq1ct
                 
                 derampl=(dispssq*ss*1.q0-dispspq*sp*1.q0)
     >                /(ss*1.q0-sp*1.q0)
                 
                 virtual=virtual+1.d0/(pi**2)*
     >                real(aimag(ffpiss)/ss/wss
     >                *aimag(ffpisp)/sp/wsp*derampl)
              endif

              
              
              
           endif
        endif
        
	virtual=4d0*virtual*e6*alpha/(2d0*pi)
        return        
      end
*****************************************************
      subroutine purefsr_disp(p1,p2,p3,p4,q,mq,virtual,qph)
#ifdef COLLIER 
        use collier
#endif
        implicit complex*16 (z)
!     implicit none
        double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
        double precision pt1(0:3),pt2(0:3),pt3(0:3),pt4(0:3),qt(0:3)
        double precision qph(40,0:3)
	double precision virtual1

        complex*16 m1,re,mq,dzphi,dmpi,m12
        complex*16 virtual,deltafsr
	complex*16 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
 	complex*16 e2,e6,x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
        complex*16 ame2,ammu2,me2,me4,mp2,mp4,ampi2,mp6,me6
        real*16 zdeltauv,zphmass2,zmu2dim
        real*16 e,e3,e5,pi2
        common/hparam/ampi
	double precision ame,ammu,convfac,alpha,pi,s,s1,ampi
	common/parameters/ame,ammu,convfac,alpha,pi

        real*16 p1q(0:3),p2q(0:3),p3q(0:3),p4q(0:3),qq(0:3)
        common/momquad/p1q,p2q,p3q,p4q,qq
        external qdot
        real*16 qdot
        external dotq
        real*16 dotq
        real*8 elmatfsrd
        complex*16 elmatifid
	complex*32 elmat,elmat1,elmat2,elmati,elmatf
	character*2 fs
	character*3 ifc
	common/finalstate/fs
	
	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min

	double complex masses(0:4)
	double complex mominv(1:10)
	double complex momvec(0:3,4)

	integer ier
c	double precision routing(5), average
	double complex routing(5), average
	

        external vpolc
        double complex vpolc,alphar
        
        common/ialpharunning/iarun 
        integer iarun 	
        
	integer i

        common/lt/ltint,quadevents,ilt
        integer ltint,quadevents,ilt
        
	common/softphotonmass/phmass
	real*8 phmass

	external dot
	double precision dot

***********CLL
      	double complex cnpi(0:1,0:2,0:2)
	double complex cnpiuv(0:1,0:2,0:2)	
	double complex dnpi(0:1,0:2,0:2,0:2)
	double complex dnpiuv(0:1,0:2,0:2,0:2)    
	         
	double complex d1
	parameter (d1 = (1.d0,0.d0))    
	double complex c1
	parameter (c1 = (1.d0,0.d0))

**********FF

	integer iffpi
	character*10 what_ffpi
	common/ffpioption/iffpi,what_ffpi
	double complex ffpiis
	double precision ffpiis2
	double complex ffpifs
	double precision ffpifs2

        common/radoption/ifisrfsr
        integer ifisrfsr

        common/pavered/p1red,p2red,p3red,p4red,kred
        double precision p1red(0:3),p2red(0:3),p3red(0:3),p4red(0:3)
        double precision kred(0:3)

        common/gvmdchecks/iscoll,fscoll
        integer iscoll,fscoll
        common/checkform/intltgvmd
        integer intltgvmd
**********LOOP FUNCTIONS
        complex*16 Cim(2)
        complex*16 Cijm(4)
        complex*16 C0m

        complex*16 Dim(3)
        complex*16 Dijm(7)
        complex*16 D0m

        complex*16 A0ap,A0am1,B0appm1
        complex*16 B0ax34pp,C0appx34pm1p
        complex*16 B0ax35pm1,B0ax45pm1
        complex*16 DB0apm1p

        complex*16 C0apx350pm1p
        complex*16 C1apx350pm1p,C2apx350pm1p,C00apx350pm1p
        complex*16 C11apx350pm1p,C12apx350pm1p,C22apx350pm1p
        complex*16 C00apx12x35m1pp,C12apx12x35m1pp,C0apx12x35m1pp

        complex*16 C0apx450pm1p
        complex*16 C1apx450pm1p,C2apx450pm1p,C00apx450pm1p
        complex*16 C11apx450pm1p,C12apx450pm1p,C22apx450pm1p

        complex*16 C0apx12x45m1pp,C1apx12x45m1pp,C1ax45x12pm1pp
        complex*16 C00apx12x45m1pp,C11apx12x45m1pp,C11ax45x12pm1pp
        complex*16 C12apx12x45m1pp

        complex*16 C0ax45x12pm1pp,C00ax45x12pm1pp,C12ax45x12pm1pp

        complex*16 C0apx45x12pm1p,C2apx45x12pm1p,C2ax45px12pm1p
        complex*16 C00apx45x12pm1p,C22apx45x12pm1p
        complex*16 C22ax45px12pm1p

        complex*16 C0ax35x12pm1pp,C1apx12x35m1pp,C1ax35x12pm1pp
        complex*16 C00ax35x12pm1pp,C11apx12x35m1pp,C11ax35x12pm1pp
        complex*16 C12ax35x12pm1pp

        complex*16 C0apx35x12pm1p,C2apx35x12pm1p,C2ax35px12pm1p
        complex*16 C00apx35x12pm1p,C22apx35x12pm1p
        complex*16 C22ax35px12pm1p

!!!!!!!!!!BOX

        complex*16 D0appx120x34x35pm1pp,D2app0x12x34x35pm1pp
        complex*16 D3app0x12x34x35pm1pp,D00app0x12x34x35pm1pp
        complex*16 D33app0x12x34x35pm1pp,D0app0x12x34x35pm1pp
        complex*16 D2appx120x34x35pm1pp,D3appx120x34x35pm1pp
        complex*16 D00appx120x34x35pm1pp,D13a0ppx12x35x34ppm1p
        complex*16 D22app0x12x34x35pm1pp,D22appx120x34x35pm1pp
        complex*16 D23app0x12x34x35pm1pp,D23appx120x34x35pm1pp
        complex*16 D33appx120x34x35pm1pp

        complex*16 C1ax12x340ppp,C00a0x12x34ppp
        complex*16 C11ax12x340ppp,C12a0x34x12ppp
        complex*16 C12a0x12x34ppp

!!!!!!!!!!!!!!!!!!!!!!

        complex*16 D0ap0x12px35x34m1ppp
        complex*16 D1ap0x12px35x34m1ppp,D2ap0x12px35x34m1ppp
        complex*16 D3ap0x12px35x34m1ppp
        complex*16 D00ap0x12px35x34m1ppp,D11ap0x12px35x34m1ppp
        complex*16 D12ap0x12px35x34m1ppp
        complex*16 D13ap0x12px35x34m1ppp,D22ap0x12px35x34m1ppp
        complex*16 D23ap0x12px35x34m1ppp,D33ap0x12px35x34m1ppp

        complex*16 D0ap0x12px45x34m1ppp
        complex*16 D1ap0x12px45x34m1ppp,D2ap0x12px45x34m1ppp
        complex*16 D3ap0x12px45x34m1ppp
        complex*16 D00ap0x12px45x34m1ppp,D11ap0x12px45x34m1ppp
        complex*16 D12ap0x12px45x34m1ppp
        complex*16 D13ap0x12px45x34m1ppp,D22ap0x12px45x34m1ppp
        complex*16 D23ap0x12px45x34m1ppp,D33ap0x12px45x34m1ppp        
************FILL COMMON FOR PAVE REDUCTION
        p1red = p1
        p2red = p2
        p3red = p3
        p4red = p4
        kred  = q

        
        re =(1.q0,0.q0)
        e = sqrt(4.q0*pi*alpha*1.q0)
        e3 =  e*e*e
        e5 =  e*e*e*e*e
        pi2 = pi*pi*1.q0
        
	e2 = 4.q0*pi*alpha*1.q0
	e6 = e2*e2*e2*1.q0
	ame2=(ame*1.q0)*(ame*1.q0)*re
	ampi2=(ampi*1.q0)*(ampi*1.q0)*re
        me2=ame2*re
        mp2=ampi2*re
        me4=me2**2*re
        me6=me2**3*re
        mp4=mp2**2*re
        mp6=mp2**3*re

c$$$        print*,'me2',me2
c$$$        print*,'mp2',mp2
        call shellquad(p1,p2,p3,p4,q,p1q,p2q,p3q,p4q,qq)

        p1p2 = qdot(p1q,p2q)*re
        p1p3 = qdot(p1q,p3q)*re
        p1p4 = qdot(p1q,p4q)*re
	kp1  = qdot(p1q,qq)*re
        p2p3 = qdot(p2q,p3q)*re
	p2p4 = qdot(p2q,p4q)*re
        kp2  = qdot(p2q,qq)*re
        kp3  = qdot(p3q,qq)*re
	p3p4 = qdot(p3q,p4q)*re
	kp4  = qdot(p4q,qq)*re

        x12 = qdot(p1q+p2q,p1q+p2q)*re
	x13 = qdot(p1q-p3q,p1q-p3q)*re
	x14 = qdot(p1q-p4q,p1q-p4q)*re
	x15 = qdot(p1q-qq,p1q-qq)*re
        x23 = qdot(p2q-p3q,p2q-p3q)*re
        x24 = qdot(p2q-p4q,p2q-p4q)*re
	x25 = qdot(p2q-qq,p2q-qq)*re
        x34 = qdot(p3q+p4q,p3q+p4q)*re
        x35 = qdot(p3q+qq,p3q+qq)*re
        x45 = qdot(p4q+qq,p4q+qq)*re

        
        elmat1=(0.q0,0.q0)
        elmat2=(0.q0,0.q0)
        
**************
	ffpifs=(1.d0,0.d0)
	ffpifs2=1.d0
        
	ffpiis=(1.d0,0.d0)
	ffpiis2=1.d0
	s=2.d0*(ame**2+p1p2)
	s1=2.d0*(ampi**2+p3p4)
	if(iffpi.gt.0) then

		call ffpi(s,ffpifs,ffpifs2) 
		call ffpi(s1,ffpiis,ffpiis2)
	        
	endif	

******ALPHA RUNNING	
	ffpiis=ffpiis*vpolc(s1)
	ffpiis2=ffpiis*CONJG(ffpiis)
	
	ffpifs=ffpifs*vpolc(s)
	ffpifs2=ffpifs*CONJG(ffpifs)
****************
!        print*,'deltauv',deltauv

        
!        print*,'deltauv',deltauv
        m1      = sqrt(mq)
        m12     = mq
        virtual =(0.q0,0.q0)

        A0ap        = (0.D0, 0.D0)
        A0am1       = (0.D0, 0.D0)
        B0appm1     = (0.D0, 0.D0)
        DB0apm1p    = (0.D0, 0.D0)
        B0ax34pp    = (0.D0, 0.D0)
        C0appx34pm1p= (0.D0, 0.D0)
        B0ax35pm1   = (0.D0, 0.D0)
        B0ax45pm1   = (0.D0, 0.D0)

        ! C 350
        C0apx350pm1p  = (0.D0, 0.D0)
        C1apx350pm1p  = (0.D0, 0.D0)
        C2apx350pm1p  = (0.D0, 0.D0)
        C00apx350pm1p = (0.D0, 0.D0)
        C11apx350pm1p = (0.D0, 0.D0)
        C12apx350pm1p = (0.D0, 0.D0)
        C22apx350pm1p = (0.D0, 0.D0)

        ! C 450
        C0apx450pm1p  = (0.D0, 0.D0)
        C1apx450pm1p  = (0.D0, 0.D0)
        C2apx450pm1p  = (0.D0, 0.D0)
        C00apx450pm1p = (0.D0, 0.D0)
        C11apx450pm1p = (0.D0, 0.D0)
        C12apx450pm1p = (0.D0, 0.D0)
        C22apx450pm1p = (0.D0, 0.D0)

        ! C 12x45
        C0apx12x45m1pp  = (0.D0, 0.D0)
        C1apx12x45m1pp  = (0.D0, 0.D0)
        C1ax45x12pm1pp  = (0.D0, 0.D0)
        C00apx12x45m1pp = (0.D0, 0.D0)
        C11apx12x45m1pp = (0.D0, 0.D0)
        C11ax45x12pm1pp = (0.D0, 0.D0)
        C12apx12x45m1pp = (0.D0, 0.D0)
        C0ax45x12pm1pp  = (0.D0, 0.D0)
        C00ax45x12pm1pp = (0.D0, 0.D0)
        C12ax45x12pm1pp = (0.D0, 0.D0)

        ! C 45x12 p
        C0apx45x12pm1p  = (0.D0, 0.D0)
        C2apx45x12pm1p  = (0.D0, 0.D0)
        C2ax45px12pm1p  = (0.D0, 0.D0)
        C00apx45x12pm1p = (0.D0, 0.D0)
        C22apx45x12pm1p = (0.D0, 0.D0)
        C22ax45px12pm1p = (0.D0, 0.D0)

        ! C 35x12
        C0ax35x12pm1pp  = (0.D0, 0.D0)
        C0apx12x35m1pp  = (0.D0, 0.D0)
        C1apx12x35m1pp  = (0.D0, 0.D0)
        C1ax35x12pm1pp  = (0.D0, 0.D0)
        C12apx12x35m1pp = (0.D0, 0.D0)
        C00ax35x12pm1pp = (0.D0, 0.D0)
        C11apx12x35m1pp = (0.D0, 0.D0)
        C00apx12x35m1pp = (0.D0, 0.D0)
        C11ax35x12pm1pp = (0.D0, 0.D0)
        C12ax35x12pm1pp = (0.D0, 0.D0)
        C0apx35x12pm1p  = (0.D0, 0.D0)
        C2apx35x12pm1p  = (0.D0, 0.D0)
        C2ax35px12pm1p  = (0.D0, 0.D0)
        C00apx35x12pm1p = (0.D0, 0.D0)
        C22apx35x12pm1p = (0.D0, 0.D0)
        C22ax35px12pm1p = (0.D0, 0.D0)

        ! D 35
        D0appx120x34x35pm1pp   = (0.D0, 0.D0)
        D2app0x12x34x35pm1pp   = (0.D0, 0.D0)
        D3app0x12x34x35pm1pp   = (0.D0, 0.D0)
        D00app0x12x34x35pm1pp  = (0.D0, 0.D0)
        D33app0x12x34x35pm1pp  = (0.D0, 0.D0)
        D0app0x12x34x35pm1pp   = (0.D0, 0.D0)
        D2appx120x34x35pm1pp   = (0.D0, 0.D0)
        D3appx120x34x35pm1pp   = (0.D0, 0.D0)
        D00appx120x34x35pm1pp  = (0.D0, 0.D0)
        D13a0ppx12x35x34ppm1p  = (0.D0, 0.D0)
        D22app0x12x34x35pm1pp  = (0.D0, 0.D0)
        D22appx120x34x35pm1pp  = (0.D0, 0.D0)
        D23app0x12x34x35pm1pp  = (0.D0, 0.D0)
        D23appx120x34x35pm1pp  = (0.D0, 0.D0)
        D33appx120x34x35pm1pp  = (0.D0, 0.D0)

        ! C extra
        C1ax12x340ppp   = (0.D0, 0.D0)
        C00a0x12x34ppp  = (0.D0, 0.D0)
        C11ax12x340ppp  = (0.D0, 0.D0)
        C12a0x34x12ppp  = (0.D0, 0.D0)
        C12a0x12x34ppp  = (0.D0, 0.D0)

#ifdef COLLIER

        call A0_cll(A0ap,mp2*c1)
        
        call A0_cll(A0am1,mq*c1)

        call B0_cll(B0ax34pp,x34*c1,mp2*c1,mp2*c1)

        call B0_cll(B0appm1,mp2*c1,mp2*c1,mq*c1)

        call DB0_cll(DB0apm1p,mp2*c1,mq*c1,mp2*c1)

        call B0_cll(B0ax35pm1,x35*c1,mp2*c1,mq*c1)

        call B0_cll(B0ax45pm1,x45*c1,mp2*c1,mq*c1)

        call C0_cll(C0appx34pm1p,mp2*c1,mp2*c1,x34*c1,
     >       mp2*c1,mq*c1,mp2*c1)

        call C_cll(cnpi,cnpiuv,mp2*c1,x35*c1,0.d0*c1,mp2*c1,
     >       mq*c1,mp2*c1,2)
        C0apx350pm1p  = cnpi(0,0,0)
        C1apx350pm1p  = cnpi(0,1,0)
        C2apx350pm1p  = cnpi(0,0,1)
	C00apx350pm1p = cnpi(1,0,0)
	C11apx350pm1p = cnpi(0,2,0)
	C12apx350pm1p = cnpi(0,1,1)
	C22apx350pm1p = cnpi(0,0,2)

        call C_cll(cnpi,cnpiuv,mp2*c1,x45*c1,0.d0*c1,mp2*c1,
     >       mq*c1,mp2*c1,2)
        C0apx450pm1p  = cnpi(0,0,0)
        C1apx450pm1p  = cnpi(0,1,0)
        C2apx450pm1p  = cnpi(0,0,1)
	C00apx450pm1p = cnpi(1,0,0)
	C11apx450pm1p = cnpi(0,2,0)
	C12apx450pm1p = cnpi(0,1,1)
	C22apx450pm1p = cnpi(0,0,2)

        call C_cll(cnpi,cnpiuv,mp2*c1,x12*c1,x45*c1,mq*c1,
     >       mp2*c1,mp2*c1,2)
        C0apx12x45m1pp  = cnpi(0,0,0)
        C1apx12x45m1pp  = cnpi(0,1,0)
        C1ax45x12pm1pp  = cnpi(0,0,1)
        C00apx12x45m1pp = cnpi(1,0,0)
        C11apx12x45m1pp = cnpi(0,2,0)
        C12apx12x45m1pp = cnpi(0,1,1)
        C11ax45x12pm1pp = cnpi(0,0,2)

        C0ax45x12pm1pp  = C0apx12x45m1pp
        C00ax45x12pm1pp = C00apx12x45m1pp
        C12ax45x12pm1pp = C12apx12x45m1pp
        C0apx45x12pm1p  = C0apx12x45m1pp

        C2apx45x12pm1p  = cnpi(0,0,1)
        C00apx45x12pm1p = cnpi(1,0,0)
        C2ax45Px12pm1p  = cnpi(0,1,0)
        C22apx45x12pm1p = cnpi(0,0,2)
        C22ax45px12pm1p = cnpi(0,2,0)

        call C_cll(cnpi,cnpiuv,mp2*c1,x12*c1,x35*c1,mq*c1,
     >       mp2*c1,mp2*c1,2)
        C0ax35x12pm1pp  = cnpi(0,0,0)
        C1apx12x35m1pp  = cnpi(0,1,0)
        C1ax35x12pm1pp  = cnpi(0,0,1)
        C00ax35x12pm1pp = cnpi(1,0,0)
        C11apx12x35m1pp = cnpi(0,2,0)
        C00apx12x35m1pp = cnpi(1,0,0)
        C12apx12x35m1pp = cnpi(0,1,1)
        C12ax35x12pm1pp = cnpi(0,1,1)
        C11ax35x12pm1pp = cnpi(0,0,2)
!!!!!!!!!!!For simmetry of PaVe                                                                                                                                                                               
        C0apx35x12pm1p  = C0ax35x12pm1pp
        c0apx12x35m1pp  = C0ax35x12pm1pp
        C2apx35x12pm1p  = cnpi(0,0,1)
        C00apx35x12pm1p = cnpi(1,0,0)
        C2ax35px12pm1p  = cnpi(0,1,0)
        C22apx35x12pm1p = cnpi(0,0,2)
        C22ax35px12pm1p = cnpi(0,2,0)
        C12ax35x12pm1pp = cnpi(0,1,1)

        call C_cll(cnpi,cnpiuv,x12*c1,x34*c1,0.d0*c1,
     >       mp2*c1,mp2*c1,mp2*c1,2)

        C1ax12x340ppp  = cnpi(0,1,0)
        C11ax12x340ppp = cnpi(0,2,0)
        C12a0x34x12ppp = cnpi(0,1,1)

        call C_cll(cnpi,cnpiuv,0.d0*c1,x12*c1,x34*c1,
     >	     mp2*c1,mp2*c1,mp2*c1,2)
        C00a0x12x34ppp = cnpi(1,0,0)
        C12a0x12x34ppp = cnpi(0,1,1)

        call D_cll(dnpi,dnpiuv,mp2*c1,0.d0*c1,x12*c1,mp2*c1,
     >       x35*c1,x34*c1,mq*c1,mp2*c1,mp2*c1,mp2*c1,2)

        D0ap0x12px35x34m1ppp  = dnpi(0,0,0,0)
        D1ap0x12px35x34m1ppp  = dnpi(0,1,0,0)
        D2ap0x12px35x34m1ppp  = dnpi(0,0,1,0)
        D3ap0x12px35x34m1ppp  = dnpi(0,0,0,1)
        D00ap0x12px35x34m1ppp = dnpi(1,0,0,0)
        D11ap0x12px35x34m1ppp = dnpi(0,2,0,0)
        D12ap0x12px35x34m1ppp = dnpi(0,1,1,0)
        D13ap0x12px35x34m1ppp = dnpi(0,1,0,1)
        D22ap0x12px35x34m1ppp = dnpi(0,0,2,0)
        D23ap0x12px35x34m1ppp = dnpi(0,0,1,1)
        D33ap0x12px35x34m1ppp = dnpi(0,0,0,2)

        call D_cll(dnpi,dnpiuv,mp2*c1,0.d0*c1,x12*c1,mp2*c1,
     >       x45*c1,x34*c1,mq*c1,mp2*c1,mp2*c1,mp2*c1,2)

        D0ap0x12px45x34m1ppp  = dnpi(0,0,0,0)
        D1ap0x12px45x34m1ppp  = dnpi(0,1,0,0)
        D2ap0x12px45x34m1ppp  = dnpi(0,0,1,0)
        D3ap0x12px45x34m1ppp  = dnpi(0,0,0,1)
        D00ap0x12px45x34m1ppp = dnpi(1,0,0,0)
        D11ap0x12px45x34m1ppp = dnpi(0,2,0,0)
        D12ap0x12px45x34m1ppp = dnpi(0,1,1,0)
        D13ap0x12px45x34m1ppp = dnpi(0,1,0,1)
        D22ap0x12px45x34m1ppp = dnpi(0,0,2,0)
        D23ap0x12px45x34m1ppp = dnpi(0,0,1,1)
        D33ap0x12px45x34m1ppp = dnpi(0,0,0,2)

#endif
        
                     
        include 'pipirad/gvmd/purefsrtot_isr.f'

        if(ifisrfsr.eq.11) then
           
           virtual = elmat1*ffpifs*CONJG(ffpiis)

        elseif(ifisrfsr.eq.01) then

           virtual = (0.q0,0.q0)

        endif

        
        include 'pipirad/gvmd/purefsrtot_fsr.f'
                
        virtual = virtual + elmat2*ffpifs2
!        print*,'PURE FSR NO CT',virtual


***************CT_s factorise over born (ISR), so here we have (|M_FSR^0|^2+M_FSR^0 x M_ISR^{0,\dagger})*dzphi****
******************************************************************************************************************
***************REMARK: THE CT_S FOR THE POINT LIKE PART ARE EVALUATED IN CT_PR EXPLOITING THE EXISTING************
******************************************POINTLIKE CALCULATION***************************************************        

        call FSR_pr(p1,p2,p3,p4,q,elmatfsrd,qph)

        call IFI_pr(p1,p2,p3,p4,q,elmatifid,qph)

        dzphi = (2.q0*B0appm1+(4.q0*mp2-mq)*DB0apm1p)

**************TO MATCH THE FACTOR IN LOOPS
        dzphi  = dzphi/(4.q0*e6)
*****************************************
        dmpi=-(2.q0*A0am1-A0ap+(4.q0*mp2-mq)*B0appm1)/mp2
           
         
*************THE IFI GIVES M_ISR*M_FSR^dagger, here i need its c.c.***********************************************
*************IN THE PT.LIKE THIS WAS NOT NEEDED SINCE CT_S WERE REAL*********************************************

        if(ifisrfsr.eq.11) then
           
           virtual = virtual
     >          +dzphi*(elmatfsrd*1.q0+CONJG(elmatifid)*1.q0)
        elseif(ifisrfsr.eq.01) then

           virtual = virtual + dzphi*elmatfsrd*1.q0

        endif


        include 'pipirad/Virtual/CT/CTmpiISR.F'
	include 'pipirad/Virtual/CT/CTmpiFSR.F'
	
	if(ifisrfsr.eq.01) then		
           virtual=virtual+elmat1*ffpifs2
	elseif(ifisrfsr.eq.11) then	
           
           virtual=virtual+elmat2*ffpifs*CONJG(ffpiis)
     >          +elmat1*ffpifs2	
	endif
        
#ifdef LOOPTOOLS
        call clearcache
#endif

        return
        end


      subroutine virtualMixedFSR_pr(p1,p2,p3,p4,q,virtual,qph)
#ifdef COLLIER
        use collier
#endif
	implicit double complex (z)
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision qph(40,0:3),kt(0:3)
	double precision virtual,virtual1
        double complex gvmd,gvmd1,gvmd2
        Complex*32 gvmdq,gvmdq1,gvmdq2,mq,derampl
        complex*32 gvmdq0,deltafsr,dfsr2
        common/maxdiff/diffmc
        real*16 diffm,diffmc
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
	real*8 e2,e6,me2,me4,me6,mp2,mp4,mp6
	double precision ame,ammu,convfac,alpha,pi,ame2,ampi2,ampi
	common/parameters/ame,ammu,convfac,alpha,pi
        common/hparam/ampi
	Complex*16 elmat2,elmat3,elmat4,elmat5,elmat6,elmat
        Complex*32 elmat1
	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min
	common/softphotonmass/phmass
	real*8 phmass

	external dot
	real*8 dot

        double complex vpolc,alphar
        common/indexfsrdisp/iFSRdisp
        integer iFSRdisp
        common/ialpharunning/iarun
        integer iarun
**********FF

	integer iffpi
	character*10 what_ffpi
	common/ffpioption/iffpi,what_ffpi
	double complex ffpiis
	double precision ffpiis2
	double complex ffpifs
	double precision ffpifs2
********FsQED

        common/compotontensor/compttens
        integer compttens
        double precision spmin,spmax
        double complex sp,ss
        double precision wsp,wss
        double complex dispsp,dispss
        double complex ffpisp,ffpiss

**********GVMD
        integer i,j,k
        double complex ct,im,cij
        double precision m_n(1:9),g_n(1:9)
        double complex m1,m2,m3
        double complex a_n(1:9),lv2(1:9),lv2diff
        common/par_gvmd/m_n,g_n,a_n

        common/resonancesgvmd/nres
        integer nres
        
	real*16 delta,delta1
        parameter (delta = 1.q-10)

*********CLL STUFF
	double complex cnpi(0:0,0:1,0:1)
	double complex cnpiuv(0:0,0:1,0:1)

	double complex a1
	parameter (a1 = (1.d0,0.d0))
	double complex c1
	parameter (c1 = (1.d0,0.d0))
	double complex b1
	parameter (b1 = (1.d0,0.d0))



****BUBBLES AND TRIANGLES

	Complex*16 A0ap
	Complex*16 B0ap0p
	Complex*16 C0appx34p0p
	Complex*16 C2appx34p0p

**********GVMD 
        Complex*16 A0am1,B0appm1,B0a0m1m1
        Complex*16 B0ax34pp,C0appx34pm1p,C0appx34pm1pf
        Complex*16 C0appx34pm1pb,C0a0ppm1m1p
        Complex*16 D0a0px34pppm1m1pp

	e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2
	ame2=ame*ame
	ampi2=ampi*ampi
        me2=ame2
        mp2=ampi2
        me4=me2**2
        me6=me2**3
        mp4=mp2**2
        mp6=mp2**3

	p1p2 = dot(p1,p2)
	p1p3 = dot(p1,p3)
	p1p4 = dot(p1,p4)
	kp1  = dot(p1,q)
	p2p3 = dot(p2,p3)
	p2p4 = dot(p2,p4)
	kp2  = dot(p2,q)
	kp3 = dot(p3,q)
	p3p4 = dot(p3,p4)
	kp4  = dot(p4,q)

	x12  = 2d0*ame2+2.d0*p1p2
	x13  = ame2+ampi2-2.d0*p1p3
	x14  = ame2+ampi2-2.d0*p1p4
	x15  = ame2-2.d0*kp1
	x23  = ame2+ampi2-2d0*p2p3
	x24  = ame2+ampi2-2d0*p2p4
	x25 = ame2-2d0*kp2
	x34 = 2d0*ampi2+2d0*p3p4
	x35 = ampi2+2d0*kp3
	x45 = ampi2+2d0*kp4
!        iffpi=0
**************
	ffpifs=(1.d0,0.d0)
	ffpifs2=1.d0

	ffpiis=(1.d0,0.d0)
	ffpiis2=1.d0
	if(iffpi.gt.0) then
		call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2)
		call ffpi(2.d0*(ampi**2+p3p4),ffpiis,ffpiis2)

	endif

******ALPHA RUNNING
	ffpiis=ffpiis*vpolc(2.d0*(ampi**2+p3p4))
	ffpiis2=ffpiis*CONJG(ffpiis)

	ffpifs=ffpifs*vpolc(2.d0*(ame**2+p1p2))
	ffpifs2=ffpifs*CONJG(ffpifs)
****************

#ifdef COLLIER
	call GetMuUV2_cll(mu2dim)
        call GetDeltaUV_cll(deltauv)
!     print*,'mu2dim',mu2dim
        call A0_cll(A0ap,a1*ampi2)

        call B0_cll(B0ap0p,b1*ampi2,b1*0d0,b1*ampi2)

        call C_cll(cnpi,cnpiuv,c1*ampi2,c1*ampi2,c1*x34,
     >       c1*ampi2,c1*0d0,c1*ampi2,1)
		C0appx34p0p=cnpi(0,0,0)
		C2appx34p0p=cnpi(0,0,1)
#endif

****************(VIRTUALE FSR REALE ISR) CROSSATO CON BORN ISR

	include 'pipirad/Virtual_new/mixedfsrisr_isr.f'

**************************
	virtual = real(elmat)*ffpiis2
****************(VIRTUALE FSR REALE ISR) CROSSATO CON BORN FSR

	include 'pipirad/Virtual_new/mixedfsrisr_fsr.f'
**************************

	virtual=virtual+real(elmat*ffpiis*CONJG(ffpifs))
!        print*,'POINTLIKE VIRTUAL',virtual
***** ------------------------------_**********
*******************GVMD***********************
*****-------------------------------**********

        if(iFSRdisp.eq.1) then
           if(iffpi.eq.2) then
!!!!! GIÀ DEFINITI COSÌ NEL FF,OSSIA A_N/CT
              im=(0.d0,1.d0)           
              do i=1,nres

                 lv2(i) = m_n(i)**2 - im*m_n(i)*g_n(i)
              enddo

!!!!!!!!!!!!!!!NOW CALCULATION AS I IMPLEMENTED IT IN THE LATEX
              do i=1,nres
                 
                 cij=-a_n(i)
                 do k=1,nres
                    
                    if(k.ne.i) then
                       cij = cij + 2.d0*a_n(k)*lv2(k)/(lv2(i)-lv2(k))
                    endif
                 enddo
                 mq=lv2(i)*1.q0
!!!!!!!!!!!!!!CENTRAL VALUE
                 call mixedfsr_quad(p1,p2,p3,p4,q,mq,gvmdq,qph)
!!!!!!!!!!!!!TO FIND DERIVATIVE
                 
                 call mixedfsr_quad(p1,p2,p3,p4,q,mq*(1.q0+1.q0*delta)
     >                ,gvmdq1,qph)
                 
                 derampl=(gvmdq1-gvmdq)/(delta*mq)
                 
                 virtual=virtual+real(a_n(i)*cij*gvmdq)
                 virtual=virtual+real(a_n(i)**2*lv2(i)*derampl)
              enddo
*---------------------------------------------
*------------------FsQED----------------------
*---------------------------------------------
              
           elseif(iffpi.eq.3) then
              if(compttens.eq.3) then
                 virtual=(0.d0,0.d0)
              endif

!--------QUI GENERO S' TRA 4*MP^2 E 9 GeV^2 (ASSUMIAMO VALIDO IL FF FINO A 3 GEV)------           
              spmin = 4.d0*mp2
              spmax = 9.d0
              call sp_disp(spmin,spmax,sp,wsp)
              call ffpi(real(sp),ffpisp,ffpifs2)
              
              call mixedfsr_disp(p1,p2,p3,p4,q,sp,dispsp,qph)

              virtual=virtual-2.d0/pi*real(aimag(ffpisp)/sp*dispsp/wsp)

!--------QUI IL DOPPIO DISPERSIVO
              
              call sp_disp(spmin,spmax,ss,wss)
              call ffpi(real(ss),ffpiss,ffpifs2)
              
              call mixedfsr_disp(p1,p2,p3,p4,q,ss,dispss,qph)

              
!-------Se |s'-s''|/s'<10^-2 il rapporto in doppia potrebbe diventare instabile--------!
!-------quindi switcho in 4-pla precisione per calcolare il ratio----------------------!
              if(abs((sp-ss)/sp).ge.1d-5) then

                 virtual=virtual+1.d0/(pi**2)*
     >                real(aimag(ffpiss)/ss/wss
     >                *aimag(ffpisp)/sp/wsp*
     >                1.d0/(ss-sp)*(dispss*ss-dispsp*sp))
              else
                 
                 call mixedfsr_quad(p1,p2,p3,p4,q,sp*1.q0,gvmdq,qph)
                 call mixedfsr_quad(p1,p2,p3,p4,q,ss*1.q0,gvmdq1,qph)

                 derampl=(gvmdq1*ss*1.q0-gvmdq*sp*1.q0)
     >                /(ss*1.q0-sp*1.q0)
                 
                 virtual=virtual+1.d0/(pi**2)*
     >                real(aimag(ffpiss)/ss/wss
     >                *aimag(ffpisp)/sp/wsp*derampl)

              endif
           endif
        endif
!     print*,'post first contrib',virtual
	virtual=4d0*virtual*e6*alpha/(2d0*pi)

	end
     

      subroutine mixedfsr_disp(p1,p2,p3,p4,q,mq,virtual,qph)
#ifdef COLLIER 
        use collier
#endif
        implicit complex*16 (z)
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision pt1(0:3),pt2(0:3),pt3(0:3),pt4(0:3),qt(0:3)
	double precision qph(40,0:3)
	double precision virtual1

        real*16 p1q(0:3),p2q(0:3),p3q(0:3),p4q(0:3),qq(0:3)
        complex*16 m1,re,mq,dzphi,m12
        complex*16 virtual,virtualf,deltafsr
	complex*16 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
 	complex*16 e2,e6,x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
        complex*16 ame2,ammu2,me2,me4,mp2,mp4,ampi2,mp6,me6
        real*16 zdeltauv,zphmass2,zmu2dim
        common/hparam/ampi
	double precision ame,ammu,convfac,alpha,pi,s,s1,ampi
	common/parameters/ame,ammu,convfac,alpha,pi

	external qdot
        external dotq
        real*16 dotq,qdot
        
        real*8 elmatisrd
        complex*16 elmatifid
	complex*32 elmat,elmat1,elmat2,elmati,elmatf
	character*2 fs
	character*3 ifc
	common/finalstate/fs
	
	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min

	double complex masses(0:4)
	double complex mominv(1:10)
	double complex momvec(0:3,4)

	integer ier
c	double precision routing(5), average
	double complex routing(5), average
	

        external vpolc
        double complex vpolc,alphar
        
        common/ialpharunning/iarun 
        integer iarun 	
        
	integer i

        common/lt/ltint,quadevents,ilt
        integer ltint,quadevents,ilt
        
	common/softphotonmass/phmass
	real*8 phmass

	external dot
	double precision dot

**********FF

	integer iffpi
	character*10 what_ffpi
	common/ffpioption/iffpi,what_ffpi
	double complex ffpiis
	double precision ffpiis2
	double complex ffpifs
	double precision ffpifs2

        common/checkform/intltgvmd
        integer intltgvmd
        common/gvmdchecks/iscoll,fscoll
        integer iscoll,fscoll
**********LOOP FUNCTIONS
        Complex*16 A0ap,A0am1,B0appm1
        Complex*16 B0ax34pp,C0appx34pm1p
        Complex*16 DB0apm1p
	double complex c1
	parameter (c1 = (1.d0,0.d0))

        re =(1.q0,0.q0)
        
	e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2
	ame2=(ame*1.q0)*(ame*1.q0)*re
	ampi2=(ampi*1.q0)*(ampi*1.q0)*re
        me2=ame2*re
        mp2=ampi2*re
        me4=me2**2*re
        me6=me2**3*re
        mp4=mp2**2*re
        mp6=mp2**3*re

        pt3=p3
        pt4=p4

        
        call shellquad(p1,p2,p3,p4,q,p1q,p2q,p3q,p4q,qq)
        
        p1p2 = qdot(p1q,p2q)*re
        p1p3 = qdot(p1q,p3q)*re
        p1p4 = qdot(p1q,p4q)*re
	kp1  = qdot(p1q,qq)*re
        p2p3 = qdot(p2q,p3q)*re
	p2p4 = qdot(p2q,p4q)*re
        kp2  = qdot(p2q,qq)*re
        kp3  = qdot(p3q,qq)*re
	p3p4 = qdot(p3q,p4q)*re
	kp4  = qdot(p4q,qq)*re

        x12 = qdot(p1q+p2q,p1q+p2q)*re
	x13 = qdot(p1q-p3q,p1q-p3q)*re
	x14 = qdot(p1q-p4q,p1q-p4q)*re
	x15 = qdot(p1q-qq,p1q-qq)*re
        x23 = qdot(p2q-p3q,p2q-p3q)*re
        x24 = qdot(p2q-p4q,p2q-p4q)*re
	x25 = qdot(p2q-qq,p2q-qq)*re
        x34 = qdot(p3q+p4q,p3q+p4q)*re
        x35 = qdot(p3q+qq,p3q+qq)*re
        x45 = qdot(p4q+qq,p4q+qq)*re

**************
	ffpifs=(1.d0,0.d0)
	ffpifs2=1.d0
        
	ffpiis=(1.d0,0.d0)
	ffpiis2=1.d0
	s=2.d0*(ame**2+p1p2)
	s1=2.d0*(ampi**2+p3p4)
	if(iffpi.gt.0) then

		call ffpi(s,ffpifs,ffpifs2) 
		call ffpi(s1,ffpiis,ffpiis2)
	        
	endif	

******ALPHA RUNNING	
	ffpiis=ffpiis*vpolc(s1)
	ffpiis2=ffpiis*CONJG(ffpiis)
	
	ffpifs=ffpifs*vpolc(s)
	ffpifs2=ffpifs*CONJG(ffpifs)
****************
!        print*,'deltauv',deltauv

        m1       = sqrt(mq)
        m12      = mq
        virtual  = (0.q0,0.q0)
#ifdef COLLIER
        call A0_cll(A0ap,mp2*c1)
        
        call A0_cll(A0am1,mq*c1)

        call B0_cll(B0ax34pp,x34*c1,mp2*c1,mp2*c1)

        call B0_cll(B0appm1,mp2*c1,mp2*c1,mq*c1)

        call DB0_cll(DB0apm1p,mp2*c1,mq*c1,mp2*c1)

        call C0_cll(C0appx34pm1p,mp2*c1,mp2*c1,x34*c1,
     >       mp2*c1,mq*c1,mp2*c1)

#endif

       
        include 'pipirad/gvmd/mixedfsrisrtot_isr.f'
        
        virtual = elmat1*ffpiis2

        include 'pipirad/gvmd/mixedfsrisrtot_fsr.f'
        
        virtual=virtual+elmat2*ffpiis*CONJG(ffpifs)

        
***************CT_s factorise over born (ISR), so here we have (|M_ISR^0|^2+M_ISR^0 x M_FSR^{0,\dagger})*dzphi****
******************************************************************************************************************
***************REMARK: THE CT_S FOR THE POINT LIKE PART ARE EVALUATED IN CT_PR EXPLOITING THE EXISTING************
******************************************POINTLIKE CALCULATION***************************************************        

        call ISR_pr(p1,p2,p3,p4,q,elmatisrd,qph)

        call IFI_pr(p1,p2,p3,p4,q,elmatifid,qph)

        
        dzphi = (2.q0*B0appm1+(4.q0*mp2-mq)*DB0apm1p)

**************TO MATCH THE FACTOR IN LOOPS
        dzphi  = dzphi/(4.q0*e6)
*****************************************
        virtual = virtual +dzphi*(elmatisrd*1.q0+elmatifid*1.q0)

        return
        end


      subroutine virtualIFI_triang_pr(p1,p2,p3,p4,q,virtualc,qph)
#ifdef COLLIER
        use collier
#endif
	implicit double complex (z)
!     implicit none
        character*2 fs
        common/finalstate/fs
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision qph(40,0:3)
	double precision virtual,virtual22
        double complex virtualc
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
	real*8 e2,e6,me2,me4,me6,mp2,mp4,mp6
	double precision ame,ammu,convfac,alpha,pi,ame2,ampi2,ampi
	common/parameters/ame,ammu,convfac,alpha,pi
        common/hparam/ampi
	Complex*16 elmat1,elmat2,elmat3,elmat4,elmat

	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min
	common/softphotonmass/phmass
	real*8 phmass

	external dot
	real*8 dot

        double complex vpolc,alphar

        common/ialpharunning/iarun
        integer iarun
**********FF

	integer iffpi
	character*10 what_ffpi
	common/ffpioption/iffpi,what_ffpi
	double complex ffpiis
	double precision ffpiis2
	double complex ffpifs
	double precision ffpifs2

        common/epssoft/eps
        real*8 eps

      	external aintsv_pions
        real*8 a,aintsv_pions
*********CLL STUFF
	double complex dnpi(0:1,0:2,0:2,0:2)
	double complex dnpiuv(0:1,0:2,0:2,0:2)

	double complex cnpi(0:0,0:1,0:1)
	double complex cnpiuv(0:0,0:1,0:1)

	double complex c1
	parameter (c1 = (1.d0,0.d0))

****TRIANGLES

	Complex*16 C0aex34x15e00
	Complex*16 C1aex34x15e00
	Complex*16 C2aex34x15e00
	Complex*16 C1ax15x34ee00
	Complex*16 C0aex15x340e0


	Complex*16 C0aex34x25e00
	Complex*16 C1aex34x25e00
	Complex*16 C2aex34x25e00
	Complex*16 C1ax25x34ee00
	Complex*16 C0aex25x340e0

	Complex*16 C0aex12ee00
	Complex*16 C1aex12ee00
	Complex*16 C0aeex120e0

****************

        complex*32  box13
        double complex C0ax15ex34m1em2
        double complex D0aex15ex250x34em1em2
        double complex D1aex15ex250x34em1em2
        double complex D3aex15ex250x34em1em2
        double complex D00aex15ex250x34em1em2
        double complex D11aex15ex250x34em1em2
        double complex D12aex15ex250x34em1em2
        double complex D13aex15ex250x34em1em2
        double complex D23aex15ex250x34em1em2
        double complex D33aex15ex250x34em1em2

        double complex m12,m22
        m12=0.d0
        m22=0.d0
        
	e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2
	ame2=ame*ame
	ampi2=ampi*ampi
        me2=ame2
        mp2=ampi2
        me4=me2**2
        me6=me2**3
        mp4=mp2**2
        mp6=mp2**3

	p1p2 = dot(p1,p2)
	p1p3 = dot(p1,p3)
	p1p4 = dot(p1,p4)
	kp1  = dot(p1,q)
	p2p3 = dot(p2,p3)
	p2p4 = dot(p2,p4)
	kp2  = dot(p2,q)
	kp3 = dot(p3,q)
	p3p4 = dot(p3,p4)
	kp4  = dot(p4,q)

	x12  = 2d0*ame2+2.d0*p1p2
	x13  = ame2+ampi2-2.d0*p1p3
	x14  = ame2+ampi2-2.d0*p1p4
	x15  = ame2-2.d0*kp1
	x23  = ame2+ampi2-2d0*p2p3
	x24  = ame2+ampi2-2d0*p2p4
	x25 = ame2-2d0*kp2
	x34 = 2d0*ampi2+2d0*p3p4
	x35 = ampi2+2d0*kp3
	x45 = ampi2+2d0*kp4

**************
	ffpifs=(1.d0,0.d0)
	ffpifs2=1.d0

	ffpiis=(1.d0,0.d0)
	ffpiis2=1.d0
	if(iffpi.gt.0) then
		call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2)
		call ffpi(2.d0*(ampi**2+p3p4),ffpiis,ffpiis2)

	endif

******ALPHA RUNNING
	ffpiis=ffpiis*vpolc(2.d0*(ampi**2+p3p4))
	ffpiis2=ffpiis*CONJG(ffpiis)

	ffpifs=ffpifs*vpolc(2.d0*(ame**2+p1p2))
	ffpifs2=ffpifs*CONJG(ffpifs)
****************

#ifdef COLLIER

        call C_cll(cnpi,cnpiuv,c1*ame2,c1*x34,c1*x15,
     >       c1*ame2,c1*0d0,c1*0d0,1)

		C0aex34x15e00=cnpi(0,0,0)
		C1aex34x15e00=cnpi(0,1,0)
		C2aex34x15e00=cnpi(0,0,1)
		C1ax15x34ee00=cnpi(0,0,1)
		C0aex15x340e0=cnpi(0,0,0)

        call C_cll(cnpi,cnpiuv,c1*ame2,c1*x34,c1*x25,
     >       c1*ame2,c1*0d0,c1*0d0,1)

		C0aex34x25e00=cnpi(0,0,0)
		C1aex34x25e00=cnpi(0,1,0)
		C2aex34x25e00=cnpi(0,0,1)
		C1ax25x34ee00=cnpi(0,0,1)
		C0aex25x340e0=cnpi(0,0,0)

        call C_cll(cnpi,cnpiuv,c1*ame2,c1*x12,c1*ame2,
     >       c1*ame2,c1*0d0,c1*0d0,1)
		C0aex12ee00=cnpi(0,0,0)
		C1aex12ee00=cnpi(0,1,0)
		C0aeex120e0=cnpi(0,0,0)


*******BOX
             
        call C0_cll(C0ax15ex34m1em2,c1*x15,c1*ame2,c1*x34,
     >       c1*m12,c1*ame2,c1*m22)
       
        call D_cll(dnpi,dnpiuv,c1*ame2,c1*x15,c1*ame2,c1*x25,
     >   c1*0.d0,c1*x34,c1*ame2,c1*m12,c1*ame2,c1*m22,2) 
        D0aex15ex250x34em1em2  = dnpi(0,0,0,0)
        D1aex15ex250x34em1em2  = dnpi(0,1,0,0)
        D3aex15ex250x34em1em2  = dnpi(0,0,0,1)
        D00aex15ex250x34em1em2 = dnpi(1,0,0,0)
        D11aex15ex250x34em1em2 = dnpi(0,2,0,0)
        D12aex15ex250x34em1em2 = dnpi(0,1,1,0)
        D13aex15ex250x34em1em2 = dnpi(0,1,0,1)
        D23aex15ex250x34em1em2 = dnpi(0,0,1,1)
        D33aex15ex250x34em1em2 = dnpi(0,0,0,2)

#endif

**************** CROSSATO CON BORN ISR

	include 'pipirad/Virtual_new/triangisr.f'
**************************

	virtualc=elmat1*CONJG(ffpiis)*ffpifs+ 
     >		elmat2*ffpiis2

**************** CROSSATO CON BORN FSR

	include 'pipirad/Virtual_new/triangfsr.f'
**************************

	virtualc=virtualc+elmat1*ffpifs2
     >	 	       +(elmat2)*(ffpiis*CONJG(ffpifs))


	virtualc=4d0*virtualc*e6*alpha/(2d0*pi)
!
!        print*,'VIRTUALE PTLIKE TRIANG',virtual
	end
****************************

      subroutine virtualIFI_Boxes_pr(p1,p2,p3,p4,q,virtualc,qph)
#ifdef COLLIER
        use collier
#endif
	implicit double complex (z)
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision qph(40,0:3)
	double precision virtual
        double complex virtualc
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
	real*8 e2,e6,me2,me4,me6,mp2,mp4,mp6
	double precision ame,ammu,convfac,alpha,pi,ame2,ampi2,ampi
	common/parameters/ame,ammu,convfac,alpha,pi
        common/hparam/ampi
	Complex*16 elmat1,elmat2,elmat3,elmat4,elmat5,elmat6,elmat7
	Complex*16 elmat8,elmat9,elmat10,elmat11,elmat12,elmat13

	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min
	common/softphotonmass/phmass
	real*8 phmass

	external dot
	real*8 dot

        double complex vpolc,alphar

        common/ialpharunning/iarun
        integer iarun

**********FF

	integer iffpi
	character*10 what_ffpi
	common/ffpioption/iffpi,what_ffpi
	double complex ffpiis
	double precision ffpiis2
	double complex ffpifs
	double precision ffpifs2


*********CLL STUFF
	double complex cnpi(0:0,0:1,0:1)
	double complex cnpiuv(0:0,0:1,0:1)

	double complex dnpi(0:1,0:2,0:2,0:2)
	double complex dnpiuv(0:1,0:2,0:2,0:2)

	double complex dnpi1(0:0,0:1,0:1,0:1)
	double complex dnpiuv1(0:0,0:1,0:1,0:1)

	double complex d1
	parameter (d1 = (1.d0,0.d0))
	double complex c1
	parameter (c1 = (1.d0,0.d0))

****TRIANGLES
	Complex*16 C0aex34x15e00
	Complex*16 C0aex15x340e0

	Complex*16 C0aex34x25e00

	Complex*16 C0apx12x35p00
	Complex*16 C1apx12x35p00
	Complex*16 C1ax12x35p00p
	Complex*16 C1apx35x120p0

	Complex*16 C0apx12x45p00
	Complex*16 C1apx12x45p00
	Complex*16 C1ax12x45p00p
	Complex*16 C1apx45x120p0

	Complex*16 C1apx34pp00
	Complex*16 C0appx340p0

	Complex*16 C0aeex120e0

*********BOXES
	Complex*16 D00a0x15x34x25eeee00
	Complex*16 D12a0x15x34x25eeee00
	Complex*16 D22a0x15x34x25eeee00
	Complex*16 D2a0x15x34x25eeee00
	Complex*16 D13a0x15x34x25eeee00
	Complex*16 D23a0x15x34x25eeee00
	Complex*16 D33a0x15x34x25eeee00
	Complex*16 D3a0x15x34x25eeee00


	Complex*16 D0apx14ex12ex350pe0
	Complex*16 D1apx14ex12ex350pe0
	Complex*16 D3apx14ex12ex350pe0

	Complex*16 D0apx24ex12ex350pe0
	Complex*16 D1apx24ex12ex350pe0
	Complex*16 D3apx24ex12ex350pe0

	Complex*16 D0apx13ex12ex450pe0
	Complex*16 D1apx13ex12ex450pe0
	Complex*16 D3apx13ex12ex450pe0

	Complex*16 D0apx23ex12ex450pe0
	Complex*16 D1apx23ex12ex450pe0
	Complex*16 D3apx23ex12ex450pe0

	Complex*16 D0apx23ex34x15p0pe0
	Complex*16 D2apx23ex34x15p0pe0
	Complex*16 D3apx23ex34x15p0pe0

	Complex*16 D0apx24ex34x15p0pe0
	Complex*16 D2apx24ex34x15p0pe0
	Complex*16 D3apx24ex34x15p0pe0

	Complex*16 D0apx13ex34x25p0pe0
	Complex*16 D2apx13ex34x25p0pe0
	Complex*16 D3apx13ex34x25p0pe0

	Complex*16 D0apx14ex34x25p0pe0
	Complex*16 D2apx14ex34x25p0pe0
	Complex*16 D3apx14ex34x25p0pe0

	Complex*16 D3ax45x13ex12ep0pe0
	Complex*16 D3ax45x23ex12ep0pe0
	Complex*16 D3ax35x14ex12ep0pe0
	Complex*16 D3ax35x24ex12ep0pe0
	Complex*16 D3apx23x15x34ep0pe0
	Complex*16 D3apx24x15x34ep0pe0
	Complex*16 D3apx13x25x34ep0pe0
	Complex*16 D3apx14x25x34ep0pe0


****************

	e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2
	ame2=ame*ame
	ampi2=ampi*ampi
        me2=ame2
        mp2=ampi2
        me4=me2**2
        me6=me2**3
        mp4=mp2**2
        mp6=mp2**3

	p1p2 = dot(p1,p2)
	p1p3 = dot(p1,p3)
	p1p4 = dot(p1,p4)
	kp1  = dot(p1,q)
	p2p3 = dot(p2,p3)
	p2p4 = dot(p2,p4)
	kp2  = dot(p2,q)
	kp3 = dot(p3,q)
	p3p4 = dot(p3,p4)
	kp4  = dot(p4,q)


	x12  = 2d0*ame2+2.d0*p1p2
	x13  = ame2+ampi2-2.d0*p1p3
	x14  = ame2+ampi2-2.d0*p1p4
	x15  = ame2-2.d0*kp1
	x23  = ame2+ampi2-2d0*p2p3
	x24  = ame2+ampi2-2d0*p2p4
	x25 = ame2-2d0*kp2
	x34 = 2d0*ampi2+2d0*p3p4
	x35 = ampi2+2d0*kp3
	x45 = ampi2+2d0*kp4

**************
	ffpifs=(1.d0,0.d0)
	ffpifs2=1.d0

	ffpiis=(1.d0,0.d0)
	ffpiis2=1.d0
	if(iffpi.gt.0) then
		call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2)
		call ffpi(2.d0*(ampi**2+p3p4),ffpiis,ffpiis2)

	endif

******ALPHA RUNNING
	ffpiis=ffpiis*vpolc(2.d0*(ampi**2+p3p4))
	ffpiis2=ffpiis*CONJG(ffpiis)

	ffpifs=ffpifs*vpolc(2.d0*(ame**2+p1p2))
	ffpifs2=ffpifs*CONJG(ffpifs)
****************

#ifdef COLLIER

	call C0_cll(C0aex34x25e00,c1*ame2,c1*x34,c1*x25,c1*ame2,
     >		    c1*0d0,c1*0d0)

	call C0_cll(C0aex34x15e00,c1*ame2,c1*x34,c1*x15,c1*ame2,
     >		    c1*0d0,c1*0d0)

	call C0_cll(C0aeex120e0,c1*ame2,c1*ame2,c1*x12,c1*0d0,
     >		    c1*ame2,c1*0d0)

	call C0_cll(C0aex15x340e0,c1*ame2,c1*x15,c1*x34,c1*0d0,
     >		    c1*ame2,c1*0d0)

	call C_cll(cnpi,cnpiuv,c1*ampi2,c1*x12,c1*x35,c1*ampi2,
     >		    c1*0d0,c1*0d0,1)

	     	C0apx12x35p00=cnpi(0,0,0)
		C1apx12x35p00=cnpi(0,1,0)
		C1ax12x35p00p=cnpi(0,0,1)

	call C_cll(cnpi,cnpiuv,c1*ampi2,c1*x35,c1*x12,c1*0d0,
     >		    c1*ampi2,c1*0d0,1)

		C1apx35x120p0=cnpi(0,1,0)

	call C_cll(cnpi,cnpiuv,c1*ampi2,c1*x12,c1*x45,c1*ampi2,
     >		    c1*0d0,c1*0d0,1)

	     	C0apx12x45p00=cnpi(0,0,0)
		C1apx12x45p00=cnpi(0,1,0)
		C1ax12x45p00p=cnpi(0,0,1)

	call C_cll(cnpi,cnpiuv,c1*ampi2,c1*x45,c1*x12,c1*0d0,
     >		    c1*ampi2,c1*0d0,1)

		C1apx45x120p0=cnpi(0,1,0)

	call C_cll(cnpi,cnpiuv,c1*ampi2,c1*x34,c1*ampi2,c1*ampi2,
     >		    c1*0d0,c1*0d0,1)

     		C1apx34pp00=cnpi(0,1,0)
     		C0appx340p0=cnpi(0,0,0)

	call D_cll(dnpi,dnpiuv,d1*0d0,d1*x15,d1*x34,d1*x25,d1*ame2,
     >    d1*ame2,d1*ame2,d1*ame2,d1*0d0,d1*0d0,2)

		D00a0x15x34x25eeee00=dnpi(1,0,0,0)
		D12a0x15x34x25eeee00=dnpi(0,1,1,0)
		D22a0x15x34x25eeee00=dnpi(0,0,2,0)
		D2a0x15x34x25eeee00=dnpi(0,0,1,0)
		D13a0x15x34x25eeee00=dnpi(0,1,0,1)
		D23a0x15x34x25eeee00=dnpi(0,0,1,1)
		D33a0x15x34x25eeee00=dnpi(0,0,0,2)
		D3a0x15x34x25eeee00=dnpi(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*ampi2,d1*x14,d1*ame2,d1*x12,
     >	   d1*ame2,d1*x35,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)

		D0apx14ex12ex350pe0=dnpi1(0,0,0,0)
		D1apx14ex12ex350pe0=dnpi1(0,1,0,0)
		D3apx14ex12ex350pe0=dnpi1(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*ampi2,d1*x24,d1*ame2,d1*x12,
     >	   d1*ame2,d1*x35,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)

		D0apx24ex12ex350pe0=dnpi1(0,0,0,0)
		D1apx24ex12ex350pe0=dnpi1(0,1,0,0)
		D3apx24ex12ex350pe0=dnpi1(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*ampi2,d1*x13,d1*ame2,d1*x12,
     >	   d1*ame2,d1*x45,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)

		D0apx13ex12ex450pe0=dnpi1(0,0,0,0)
		D1apx13ex12ex450pe0=dnpi1(0,1,0,0)
		D3apx13ex12ex450pe0=dnpi1(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*ampi2,d1*x23,d1*ame2,d1*x12,
     >	   d1*ame2,d1*x45,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)

		D0apx23ex12ex450pe0=dnpi1(0,0,0,0)
		D1apx23ex12ex450pe0=dnpi1(0,1,0,0)
		D3apx23ex12ex450pe0=dnpi1(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*ampi2,d1*x23,d1*ame2,d1*x34,
     >	   d1*x15,d1*ampi2,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)

		D0apx23ex34x15p0pe0=dnpi1(0,0,0,0)
		D2apx23ex34x15p0pe0=dnpi1(0,0,1,0)
		D3apx23ex34x15p0pe0=dnpi1(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*ampi2,d1*x24,d1*ame2,d1*x34,
     >	   d1*x15,d1*ampi2,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)

		D0apx24ex34x15p0pe0=dnpi1(0,0,0,0)
		D2apx24ex34x15p0pe0=dnpi1(0,0,1,0)
		D3apx24ex34x15p0pe0=dnpi1(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*ampi2,d1*x13,d1*ame2,d1*x34,
     >	   d1*x25,d1*ampi2,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)

		D0apx13ex34x25p0pe0=dnpi1(0,0,0,0)
		D2apx13ex34x25p0pe0=dnpi1(0,0,1,0)
		D3apx13ex34x25p0pe0=dnpi1(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*ampi2,d1*x14,d1*ame2,d1*x34,
     >	   d1*x25,d1*ampi2,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)

		D0apx14ex34x25p0pe0=dnpi1(0,0,0,0)
		D2apx14ex34x25p0pe0=dnpi1(0,0,1,0)
		D3apx14ex34x25p0pe0=dnpi1(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*x45,d1*x13,d1*ame2,d1*x12,
     >	   d1*ame2,d1*ampi2,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)
		D3ax45x13ex12ep0pe0=dnpi1(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*x45,d1*x23,d1*ame2,d1*x12,
     >	   d1*ame2,d1*ampi2,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)
		D3ax45x23ex12ep0pe0=dnpi1(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*x35,d1*x14,d1*ame2,d1*x12,
     >	   d1*ame2,d1*ampi2,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)
		D3ax35x14ex12ep0pe0=dnpi1(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*x35,d1*x24,d1*ame2,d1*x12,
     >	   d1*ame2,d1*ampi2,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)
		D3ax35x24ex12ep0pe0=dnpi1(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*ampi2,d1*x23,d1*x15,d1*x34,
     >	   d1*ame2,d1*ampi2,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)
		D3apx23x15x34ep0pe0=dnpi1(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*ampi2,d1*x24,d1*x15,d1*x34,
     >	   d1*ame2,d1*ampi2,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)
		D3apx24x15x34ep0pe0=dnpi1(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*ampi2,d1*x13,d1*x25,d1*x34,
     >	   d1*ame2,d1*ampi2,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)
		D3apx13x25x34ep0pe0=dnpi1(0,0,0,1)

	call D_cll(dnpi1,dnpiuv1,d1*ampi2,d1*x14,d1*x25,d1*x34,
     >	   d1*ame2,d1*ampi2,d1*0d0,d1*ampi2,d1*ame2,d1*0d0,1)
		D3apx14x25x34ep0pe0=dnpi1(0,0,0,1)

#endif


**************** CROSSATO CON BORN ISR

	include 'pipirad/Virtual_new/boxisr.f'
!**************************
	virtualc=elmat1*CONJG(ffpiis)*ffpifs+elmat2*ffpiis2

**************** CROSSATO CON BORN FSR

	include 'pipirad/Virtual_new/boxfsr.f'

********************

	virtualc=virtualc+elmat1*ffpifs2+
     >		elmat2*ffpiis*CONJG(ffpifs)	


	virtualc=4d0*virtualc*e6*alpha/(2d0*pi)

!        print*,'VIRTUAL PTLIKE',virtual
	end

*********************

      subroutine virtualIFI_Pent_4pt_pr(p1,p2,p3,p4,q,virtual,qph)
#ifdef COLLIER
        use collier
#endif
	implicit double complex (z)
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision pt1(0:3),pt2(0:3),pt3(0:3),pt4(0:3),qt(0:3)
	double precision qph(40,0:3)
	double complex virtual
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45,p1km1,p2km1
	real*8 e2,e6,me2,me4,me6,mp2,mp4,mp6
	double precision ame,ammu,convfac,alpha,pi,ame2,ampi2,ampi
	common/parameters/ame,ammu,convfac,alpha,pi
        common/hparam/ampi
	Complex*16 elmat1,elmat2,elmati,elmatf
	integer i
	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min
	common/softphotonmass/phmass
	real*8 phmass

        common/instabilities/instabilities
        integer instabilities


	external dot
	real*8 dot

        double complex vpolc,alphar

        common/ialpharunning/iarun
        integer iarun

**********FF

	integer iffpi
	character*10 what_ffpi
	common/ffpioption/iffpi,what_ffpi
	double complex ffpiis,Cffpiis
	double precision ffpiis2
	double complex ffpifs,Cffpifs
	double precision ffpifs2


*********CLL STUFF

	double complex dnpi(0:1,0:2,0:2,0:2)
	double complex dnpiuv(0:1,0:2,0:2,0:2)

	double complex dnpi1(0:0,0:1,0:1,0:1)
	double complex dnpiuv1(0:0,0:1,0:1,0:1)

	double complex d1
	parameter (d1 = (1.d0,0.d0))

	double complex e1
	parameter (e1 = (1.d0,0.d0))

	double complex masses(0:4)
	double complex mominv(1:10)
	double complex momvec(0:3,4)

*************FSR FUNCTIONS
	Complex*16 D00ax12p0px45x3500pp
	Complex*16 D0ax12p0px45x3500pp
	Complex*16 D1ax12p0px45x3500pp
	Complex*16 D12ax12p0px45x3500pp
	Complex*16 D22ax12p0px45x3500pp
	Complex*16 D2ax12p0px45x3500pp
	Complex*16 D23ax12p0px45x3500pp
	Complex*16 D22ax12p0px35x4500pp
	Complex*16 D2ax12p0px35x4500pp
	Complex*16 D12ax12x350x45pp00pp

	Complex*16 D12ax12p0px35x4500pp
	Complex*16 D1ax12p0px35x4500pp
	Complex*16 D12ax12x450x35pp00pp

************ISR FUNCTIONS
	Complex*16 C0a0x25eee0
	Complex*16 C0apx24e0pe
	Complex*16 C0apx23e0pe

	Complex*16 D0a0x13pex24x25eep0
	Complex*16 D1a0x13pex24x25eep0
	Complex*16 D1a0x24px25x13eeep0
	Complex*16 D1apx250x24x13ep0ee

	Complex*16 D00aex34e0x25x15e00e
	Complex*16 D0aex34e0x25x15e00e
	Complex*16 D11aex34e0x25x15e00e
	Complex*16 D1aex34e0x25x15e00e
	Complex*16 D12aex34e0x25x15e00e
	Complex*16 D22aex34e0x25x15e00e
	Complex*16 D2aex34e0x25x15e00e
	Complex*16 D13aex34e0x25x15e00e
	Complex*16 D23aex34e0x25x15e00e
	Complex*16 D3aex34e0x25x15e00e

	Complex*16 D33aex34e0x25x15e00e
	Complex*16 D0a0x14pex23x25eep0
	Complex*16 D1a0x14pex23x25eep0
	Complex*16 D13aex34e0x15x25e00e
	Complex*16 D33aex34e0x15x25e00e
	Complex*16 D23ax34x150x25ee00ee
	Complex*16 D1apx250x23x14ep0ee
	Complex*16 D1a0x23px25x14eeep0
********************
*************GVMD COMMON
        common/gvmdrad/indexnores
        integer indexnores
!        indexnores = 0

	e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2
	ame2=ame*ame
	ampi2=ampi*ampi
        me2=ame2
        mp2=ampi2
        me4=me2**2
        me6=me2**3
        mp4=mp2**2
        mp6=mp2**3

	p1p2 = dot(p1,p2)
	p1p3 = dot(p1,p3)
	p1p4 = dot(p1,p4)
	kp1  = dot(p1,q)
	p2p3 = dot(p2,p3)
	p2p4 = dot(p2,p4)
	kp2  = dot(p2,q)
	kp3 = dot(p3,q)
	p3p4 = dot(p3,p4)
	kp4  = dot(p4,q)
	p1km1=1d0/kp1
	p2km1=1d0/kp2

	x12  = 2d0*ame2+2.d0*p1p2
	x13  = ame2+ampi2-2.d0*p1p3
	x14  = ame2+ampi2-2.d0*p1p4
	x15  = ame2-2.d0*kp1
	x23  = ame2+ampi2-2d0*p2p3
	x24  = ame2+ampi2-2d0*p2p4
	x25 = ame2-2d0*kp2
	x34 = 2d0*ampi2+2d0*p3p4
	x35 = ampi2+2d0*kp3
	x45 = ampi2+2d0*kp4


**************
	ffpifs=(1.d0,0.d0)
	ffpifs2=1.d0

	ffpiis=(1.d0,0.d0)
	ffpiis2=1.d0
	if(iffpi.gt.0) then
		call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2)
		call ffpi(2.d0*(ampi**2+p3p4),ffpiis,ffpiis2)

	endif	

        Cffpiis=CONJG(ffpiis)
        Cffpifs=CONJG(ffpifs)

        if(indexnores.eq.1) then
           ffpiis=(1.d0,0.d0)
           ffpifs=(1.d0,0.d0)
	endif

******ALPHA RUNNING
	ffpiis=ffpiis*vpolc(2.d0*(ampi**2+p3p4))
        Cffpiis=Cffpiis*CONJG(vpolc(2.d0*(ampi**2+p3p4)))

	ffpifs=ffpifs*vpolc(2.d0*(ame**2+p1p2))
        Cffpifs=Cffpifs*CONJG(vpolc(2.d0*(ame**2+p1p2)))
        
        ffpiis2=ffpiis*Cffpiis
        ffpifs2=ffpifs*Cffpifs
****************
c$$$
c$$$	print*,ffpiis
c$$$        print*,ffpifs
c$$$        print*,Cffpiis
c$$$        print*,Cffpifs
c$$$        print*,ffpifs2
c$$$        print*,ffpiis2
!*********PENTAGONS WITH PHOTON ATTACHED TO FS
#ifdef COLLIER

       	call D_cll(dnpi,dnpiuv,d1*x12,d1*ampi2,d1*0d0,d1*ampi2,
     >	   d1*x45,d1*x35,d1*0d0,d1*0d0,d1*ampi2,d1*ampi2,2)
		D00ax12p0px45x3500pp=dnpi(1,0,0,0)
		D0ax12p0px45x3500pp=dnpi(0,0,0,0)
		D1ax12p0px45x3500pp=dnpi(0,1,0,0)
		D12ax12p0px45x3500pp=dnpi(0,1,1,0)
		D22ax12p0px45x3500pp=dnpi(0,0,2,0)
		D2ax12p0px45x3500pp=dnpi(0,0,1,0)
		D23ax12p0px45x3500pp=dnpi(0,0,1,1)
		D22ax12p0px35x4500pp=dnpi(0,0,0,2)
		D2ax12p0px35x4500pp=dnpi(0,0,0,1)
		D12ax12x350x45pp00pp=dnpi(0,1,0,1)

	call D_cll(dnpi,dnpiuv,d1*x12,d1*ampi2,d1*0d0,d1*ampi2,
     >	   d1*x35,d1*x45,d1*0d0,d1*0d0,d1*ampi2,d1*ampi2,2)

		D12ax12p0px35x4500pp=dnpi(0,1,1,0)
		D1ax12p0px35x4500pp=dnpi(0,1,0,0)
		D12ax12x450x35pp00pp=dnpi(0,1,0,1)

#endif

                virtual=0d0

****************DIAGRAMMA 1
 include 'pipirad/Virtual_new/pent1_4pt.f'

 	elmati=elmat1*Cffpiis*ffpifs
 	elmatf=elmat2*ffpifs*Cffpifs

 	virtual=elmati+elmatf

****************DIAGRAMMA 2

 include 'pipirad/Virtual_new/pent2_4pt.f'

 	elmati=elmat1*Cffpiis*ffpifs

 	elmatf=elmat2*ffpifs*Cffpifs

 	virtual=virtual+elmati+elmatf
!        print*,'virtual 12',virtual
        virtual=0.d0
#ifdef COLLIER
	call C0_cll(C0a0x25eee0,d1*0d0,d1*x25,d1*ame2,d1*ame2,
     >	             d1*ame2,d1*0d0)

	call C0_cll(C0apx24e0pe,d1*ampi2,d1*x24,d1*ame2,d1*0d0,
     >	             d1*ampi2,d1*ame2)

     	call C0_cll(C0apx23e0pe,d1*ampi2,d1*x23,d1*ame2,d1*0d0,
     >	             d1*ampi2,d1*ame2)

	call D_cll(dnpi,dnpiuv,d1*ame2,d1*x34,d1*ame2,d1*0d0,
     >	   d1*x25,d1*x15,d1*ame2,d1*0d0,d1*0d0,d1*ame2,2)

		D00aex34e0x25x15e00e=dnpi(1,0,0,0)
		D0aex34e0x25x15e00e=dnpi(0,0,0,0)
		D11aex34e0x25x15e00e=dnpi(0,2,0,0)
		D1aex34e0x25x15e00e=dnpi(0,1,0,0)
		D12aex34e0x25x15e00e=dnpi(0,1,1,0)
		D22aex34e0x25x15e00e=dnpi(0,0,2,0)
		D2aex34e0x25x15e00e=dnpi(0,0,1,0)
		D13aex34e0x25x15e00e=dnpi(0,1,0,1)
		D23aex34e0x25x15e00e=dnpi(0,0,1,1)
		D3aex34e0x25x15e00e=dnpi(0,0,0,1)
		D33aex34e0x25x15e00e=dnpi(0,0,0,2)


         call D_cll(dnpi1,dnpiuv1,d1*ampi2,d1*x25,d1*0d0,d1*x24,
     >         d1*x13,d1*ame2,d1*ampi2,d1*0d0,d1*ame2,d1*ame2,1)

            D0a0x13pex24x25eep0=dnpi1(0,0,0,0)
            D1a0x13pex24x25eep0=dnpi1(0,0,1,0)


            D1apx250x24x13ep0ee=dnpi1(0,1,0,0)
            D1a0x24px25x13eeep0=dnpi1(0,0,0,1)

         call D_cll(dnpi1,dnpiuv1,d1*ampi2,d1*x25,d1*0d0,d1*x23,
     >               d1*x14,d1*ame2,d1*ampi2,d1*0d0,d1*ame2,d1*ame2,1)
           D0a0x14pex23x25eep0=dnpi1(0,0,0,0)
           D1a0x14pex23x25eep0=dnpi1(0,0,1,0)
           D1apx250x23x14ep0ee=dnpi1(0,1,0,0)
           D1a0x23px25x14eeep0=dnpi1(0,0,0,1)


	call D_cll(dnpi,dnpiuv,d1*ame2,d1*x34,d1*ame2,d1*0d0,
     >	  d1*x15,d1*x25,d1*ame2,d1*0d0,d1*0d0,d1*ame2,2)

     		D33aex34e0x15x25e00e=dnpi(0,0,0,2)
		D13aex34e0x15x25e00e=dnpi(0,1,0,1)

	call D_cll(dnpi,dnpiuv,d1*x34,d1*x15,d1*0d0,d1*x25,
     >	   d1*ame2,d1*ame2,d1*0d0,d1*0d0,d1*ame2,d1*ame2,2)

                D23ax34x150x25ee00ee=dnpi(0,0,1,1)

#endif


****************DIAGRAMMA 3
 include 'pipirad/Virtual_new/pent3_4pt.f'

c$$$      print*,elmat1
c$$$      print*,elmat2
 	elmati=elmat1*ffpiis*Cffpiis
 	elmatf=elmat2*Cffpifs*ffpiis
 	virtual=virtual+elmati+elmatf
c$$$        print*,'diag3',elmati,elmatf
****************DIAGRAMMA 4

 include 'pipirad/Virtual_new/pent4_4pt.f'

c$$$      print*,elmat1
c$$$      print*,elmat2
 	elmati=elmat1*ffpiis*Cffpiis
 	elmatf=elmat2*Cffpifs*ffpiis
 	virtual=virtual+elmati+elmatf
c$$$ 	print*,'diag4',elmati,elmatf
!        print*,'virtual fin',virtual
	virtual=4d0*virtual*e6*alpha/(2d0*pi)
c$$$
c$$$        if(instabilities.eq.1) then
c$$$        print*,virtual
c$$$        print*, C0a0x25eee0
c$$$        print*, C0apx24e0pe
c$$$        print*, C0apx23e0pe
c$$$
c$$$        print*, D0a0x13pex24x25eep0
c$$$        print*, D1a0x13pex24x25eep0
c$$$        print*, D1a0x24px25x13eeep0
c$$$        print*, D1apx250x24x13ep0ee
c$$$
c$$$        print*, D00aex34e0x25x15e00e
c$$$        print*, D0aex34e0x25x15e00e
c$$$        print*, D11aex34e0x25x15e00e
c$$$        print*, D1aex34e0x25x15e00e
c$$$        print*, D12aex34e0x25x15e00e
c$$$        print*, D22aex34e0x25x15e00e
c$$$        print*, D2aex34e0x25x15e00e
c$$$        print*, D13aex34e0x25x15e00e
c$$$        print*, D23aex34e0x25x15e00e
c$$$        print*, D3aex34e0x25x15e00e
c$$$
c$$$        print*, D33aex34e0x25x15e00e
c$$$        print*, D0a0x14pex23x25eep0
c$$$        print*, D1a0x14pex23x25eep0
c$$$        print*, D13aex34e0x15x25e00e
c$$$        print*, D33aex34e0x15x25e00e
c$$$        print*, D23ax34x150x25ee00ee
c$$$        print*, D1apx250x23x14ep0ee
c$$$        print*, D1a0x23px25x14eeep0
c$$$           endif
	end

      subroutine virtualIFI_Pent_5pt_pr(p1,p2,p3,p4,q,virtual,qph,i)
#ifdef COLLIER
        use collier
#endif
	implicit none
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision pt1(0:3),pt2(0:3),pt3(0:3),pt4(0:3),qt(0:3)
c	double precision qph(40,0:3),routing1(5),routing2(5)
	double precision qph(40,0:3)
        double complex routing1(5),routing2(5)

	double complex virtual
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
	real*8 e2,e6
	double precision ame,ammu,convfac,alpha,pi,ame2,ampi2,ampi
	common/parameters/ame,ammu,convfac,alpha,pi
        common/hparam/ampi
	Complex*16 elmat1,elmat2,elmati,elmatf
	integer i
	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min
	common/softphotonmass/phmass
	real*8 phmass

        common/instabilities/instabilities
        integer instabilities


	external dot
	real*8 dot
        common/lt/ltint,quadevents,ilt
        integer ltint,quadevents,ilt



        double complex vpolc,alphar

        common/ialpharunning/iarun
        integer iarun

*******CLL STUFF

	double complex e1
	parameter (e1 = (1.d0,0.d0))

	double complex masses(0:4)
	double complex mominv(1:10)
	double complex momvec(0:3,4)

****************FIVE POINT FUNCTIONS
	e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2
	ame2=ame*ame
	ampi2=ampi*ampi
!	instabilities=1
!	print*,"PENT"
	if(i.eq.1) then
		pt1=p1
		pt2=p2
		pt3=p3
		pt4=p4
		qt=q
	elseif(i.eq.2) then
		pt1=p2
		pt2=p1
		pt3=p3
		pt4=p4
		qt=q
	endif

	p1p2 = dot(pt1,pt2)
	p1p3 = dot(pt1,pt3)
	p1p4 = dot(pt1,pt4)
	kp1  = dot(pt1,qt)
	p2p3 = dot(pt2,pt3)
	p2p4 = dot(pt2,pt4)
	kp2  = dot(pt2,qt)
	kp3 = dot(pt3,qt)
	p3p4 = dot(pt3,pt4)
	kp4  = dot(pt4,qt)

********************ROUTING 1

	momvec(:,1)=e1*pt2
	momvec(:,2)=e1*(1d0*pt2-1d0*pt4)
	momvec(:,3)=e1*(1d0*pt3-1d0*pt1)
	momvec(:,4)=-e1*pt1

	mominv(1)=e1*ame2
	mominv(2)=e1*ampi2
	mominv(3)=e1*0d0
	mominv(4)=e1*ampi2
	mominv(5)=e1*ame2
	mominv(6)=e1*(ampi2+ame2-2d0*p2p4)
	mominv(7)=e1*(ampi2+2d0*kp4)
	mominv(8)=e1*(ampi2+2d0*kp3)
	mominv(9)=e1*(ame2+ampi2-2.d0*p1p3)
	mominv(10)=e1*(2d0*ame2+2d0*p1p2)

	masses(0)=e1*ame2
	masses(1)=e1*0d0
	masses(2)=e1*ampi2
	masses(3)=e1*ampi2
	masses(4)=e1*0d0
        if(instabilities.eq.1) then
	call routing_pipi(pt1,pt2,pt3,pt4,qt,
     >		     momvec,mominv,masses,1,1,elmat1,elmat2)
        routing1(1)=elmat1
        routing2(1)=elmat2
        endif

********************ROUTING 2 PER LE FULL TENSOR

	momvec(:,1)=-e1*(1d0*pt2+1d0*pt1)
	momvec(:,2)=-e1*pt2
	momvec(:,3)=-e1*pt4
	momvec(:,4)=-e1*(1d0*p4+1d0*qt)

	mominv(1)=e1*(2d0*ame2+2d0*p1p2)
	mominv(2)=e1*ame2
	mominv(3)=e1*(ampi2+ame2-2d0*p2p4)
	mominv(4)=e1*0d0
	mominv(5)=e1*(ampi2+2d0*kp4)
	mominv(6)=e1*ame2
	mominv(7)=e1*(ampi2+2d0*kp3)
	mominv(8)=e1*(ame2+ampi2-2.d0*p1p3)
	mominv(9)=e1*ampi2
	mominv(10)=e1*ampi2

	masses(0)=e1*0d0
	masses(1)=e1*0d0
	masses(2)=e1*ame2
	masses(3)=e1*ampi2
	masses(4)=e1*ampi2

	if(instabilities.eq.1) then

	call routing_pipi(pt1,pt2,pt3,pt4,qt,
     >		     momvec,mominv,masses,2,1,elmat1,elmat2)
        routing1(2)=elmat1
        routing2(2)=elmat2
	endif
!!	print*,elmat1
!	print*,elmat2
!	write(17,*) real(elmat1)

********************ROUTING 3 PER LE FULL TENSOR

	momvec(:,1)=e1*(pt4)
	momvec(:,2)=-e1*(1d0*pt3+1d0*qt)
	momvec(:,3)=-e1*qt
	momvec(:,4)=e1*(1d0*pt4-1d0*pt2)

	mominv(1)=e1*ampi2
	mominv(2)=e1*(2d0*ame2+2d0*p1p2)
	mominv(3)=e1*ampi2
	mominv(4)=e1*(ame2+ampi2-2.d0*p1p3)
	mominv(5)=e1*(ampi2+ame2-2d0*p2p4)
	mominv(6)=e1*(ampi2+2d0*kp3)
	mominv(7)=e1*(ampi2+2d0*kp4)
	mominv(8)=e1*ame2
	mominv(9)=e1*0d0
	mominv(10)=e1*ame2

	masses(0)=e1*ampi2
	masses(1)=e1*0d0
	masses(2)=e1*0d0
	masses(3)=e1*ampi2
	masses(4)=e1*ame2

	if(instabilities.eq.1) then

	call routing_pipi(pt1,pt2,pt3,pt4,qt,
     >		     momvec,mominv,masses,3,1,elmat1,elmat2)
        routing1(3)=elmat1
        routing2(3)=elmat2
	endif
********************ROUTING 4 PER LE FULL TENSOR

	momvec(:,1)=-e1*pt3
	momvec(:,2)=e1*(pt4+qt)
	momvec(:,3)=e1*qt
	momvec(:,4)=e1*(pt1-pt3)

	mominv(1)=e1*ampi2
	mominv(2)=e1*(2d0*ame2+2d0*p1p2)
	mominv(3)=e1*ampi2
	mominv(4)=e1*(ampi2+ame2-2d0*p2p4)
	mominv(5)=e1*(ame2+ampi2-2.d0*p1p3)
	mominv(6)=e1*(ampi2+2d0*kp4)
	mominv(7)=e1*(ampi2+2d0*kp3)
	mominv(8)=e1*ame2
	mominv(9)=e1*0d0
	mominv(10)=e1*ame2

	masses(0)=e1*ampi2
	masses(1)=e1*0d0
	masses(2)=e1*0d0
	masses(3)=e1*ampi2
	masses(4)=e1*ame2


	call routing_pipi(pt1,pt2,pt3,pt4,qt,
     >		     momvec,mominv,masses,4,1,elmat1,elmat2)
        routing1(4)=elmat1
        routing2(4)=elmat2
********************ROUTING 5 PER LE FULL TENSOR

	momvec(:,1)=e1*(1d0*pt1+1d0*pt2)
	momvec(:,2)=e1*(pt1)
	momvec(:,3)=e1*pt3
	momvec(:,4)=e1*(1d0*pt3+1d0*qt)

	mominv(1)=e1*(2d0*ame2+2d0*p1p2)
	mominv(2)=e1*ame2
	mominv(3)=e1*(ame2+ampi2-2.d0*p1p3)
	mominv(4)=e1*0d0
	mominv(5)=e1*(ampi2+2d0*kp3)
	mominv(6)=e1*ame2
	mominv(7)=e1*(ampi2+2d0*kp4)
	mominv(8)=e1*(ampi2+ame2-2d0*p2p4)
	mominv(9)=e1*ampi2
	mominv(10)=e1*ampi2

	masses(0)=e1*0d0
	masses(1)=e1*0d0
	masses(2)=e1*ame2
	masses(3)=e1*ampi2
	masses(4)=e1*ampi2
	if(instabilities.eq.1) then

	call routing_pipi(pt1,pt2,pt3,pt4,qt,
     >		     momvec,mominv,masses,5,1,elmat1,elmat2)

        routing1(5)=elmat1
        routing2(5)=elmat2
c$$$        call find_closest_two(routing1,elmat1)
c$$$
c$$$        call find_closest_two(routing2,elmat2)
!	print*,'FSR'
!        print*,routing1
!        print*,elmat1
!        print*,'reldiff',elmat1/routing1(4)-1d0
!        print*,routing2
!        print*,elmat2
!        print*,'reldiff',elmat2/routing2(4)-1d0
	endif



        if(i.eq.2) then
           elmat2=-elmat2
        endif
c$$$        print*,'fsr',routing1
c$$$        print*,'fsr',routing2
        virtual=elmat1+elmat2

***********NOW WE EVAUALATE THE TWO ISR PENTAGONS
        if(ltint.eq.0) then
	if(i.eq.1) then
		pt1=p1
		pt2=p2
		pt3=p3
		pt4=p4
		qt=q
	elseif(i.eq.2) then
		pt1=p1
		pt2=p2
		pt3=p4
		pt4=p3
		qt=q
	endif

	p1p2 = dot(pt1,pt2)
	p1p3 = dot(pt1,pt3)
	p1p4 = dot(pt1,pt4)
	kp1  = dot(pt1,qt)
	p2p3 = dot(pt2,pt3)
	p2p4 = dot(pt2,pt4)
	kp2  = dot(pt2,qt)
	kp3 = dot(pt3,qt)
	p3p4 = dot(pt3,pt4)
	kp4  = dot(pt4,qt)

********************QUI PER LE TENSORIALI CON LA CONVENZIONE DI MATHEMATICA SULL'ESPRESSIONE INIZIALE

	momvec(:,1)=-e1*pt3
	momvec(:,2)=-e1*(1d0*pt3+1d0*pt4)
	momvec(:,3)=-e1*(1d0*pt1-1d0*qt)
	momvec(:,4)=-e1*pt1

	mominv(1)=e1*ampi2
	mominv(2)=e1*ampi2
	mominv(3)=e1*(ame2)
	mominv(4)=e1*0d0
	mominv(5)=e1*ame2
	mominv(6)=e1*(2d0*ampi2+2d0*p3p4)
	mominv(7)=e1*(ame2+ampi2-2d0*p2p4)
	mominv(8)=e1*(ame2-2d0*kp2)
	mominv(9)=e1*(ame2-2d0*kp1)
	mominv(10)=e1*(ame2+ampi2-2.d0*p1p3)

	masses(0)=e1*0d0
	masses(1)=e1*ampi2
	masses(2)=e1*0d0
	masses(3)=e1*ame2
	masses(4)=e1*ame2
	if(instabilities.eq.1) then
	call routing_pipi(pt1,pt2,pt3,pt4,qt,
     >		     momvec,mominv,masses,1,3,elmat1,elmat2)
        routing1(1)=elmat1
        routing2(1)=elmat2
	endif

!   	if(i.eq.2) then
!     		elmat1=-elmat1
!	endif
!	virtual=elmat1
********************ROUTING 2 PER LE FULL TENSOR

	momvec(:,1)=e1*pt3
	momvec(:,2)=-e1*pt4
	momvec(:,3)=-e1*(1d0*pt1-1d0*pt3)
	momvec(:,4)=-e1*(1d0*pt4-1d0*pt2)

	mominv(1)=e1*ampi2
	mominv(2)=e1*(2d0*ampi2+2d0*p3p4)
	mominv(3)=e1*(ame2-2d0*kp2)
	mominv(4)=e1*0d0
	mominv(5)=e1*(ame2+ampi2-2d0*p2p4)
	mominv(6)=e1*ampi2
	mominv(7)=e1*ame2
	mominv(8)=e1*ame2
	mominv(9)=e1*(ame2+ampi2-2.d0*p1p3)
	mominv(10)=e1*(ame2-2d0*kp1)

	masses(0)=e1*ampi2
	masses(1)=e1*0d0
	masses(2)=e1*0d0
	masses(3)=e1*ame2
	masses(4)=e1*ame2
	if(instabilities.eq.1) then
	call routing_pipi(pt1,pt2,pt3,pt4,qt,
     >		     momvec,mominv,masses,2,3,elmat1,elmat2)
        routing1(2)=elmat1
        routing2(2)=elmat2
	endif
!   	if(i.eq.2) then
!     		elmat1=-elmat1
!	endif
!	virtual=elmat1

********************ROUTING 3 PER LE FULL TENSOR

	momvec(:,1)=e1*(pt3+pt4)
	momvec(:,2)=e1*pt2
	momvec(:,3)=e1*pt4
	momvec(:,4)=e1*(pt2-qt)

	mominv(1)=e1*(2d0*ampi2+2d0*p3p4)
	mominv(2)=e1*(ame2-2d0*kp1)
	mominv(3)=e1*(ame2+ampi2-2d0*p2p4)
	mominv(4)=e1*(ame2+ampi2-2.d0*p1p3)
	mominv(5)=e1*(ame2-2d0*kp2)
	mominv(6)=e1*ame2
	mominv(7)=e1*ampi2
	mominv(8)=e1*0d0
	mominv(9)=e1*ampi2
	mominv(10)=e1*ame2

	masses(0)=e1*0d0
	masses(1)=e1*0d0
	masses(2)=e1*ame2
	masses(3)=e1*ampi2
	masses(4)=e1*ame2
	if(instabilities.eq.1) then
	call routing_pipi(pt1,pt2,pt3,pt4,qt,
     >		     momvec,mominv,masses,3,3,elmat1,elmat2)
        routing1(3)=elmat1
        routing2(3)=elmat2
	endif
!   	if(i.eq.2) then
!     		elmat1=-elmat1
!	endif
!	virtual=elmat1
!!!
********************ROUTING 4 PER LE FULL TENSOR

	momvec(:,1)=-e1*pt2
	momvec(:,2)=e1*(pt1-qt)
	momvec(:,3)=-e1*qt
	momvec(:,4)=e1*(pt4-pt2)

	mominv(1)=e1*ame2
	mominv(2)=e1*(2d0*ampi2+2d0*p3p4)
	mominv(3)=e1*ame2
	mominv(4)=e1*(ame2+ampi2-2.d0*p1p3)
	mominv(5)=e1*(ame2+ampi2-2d0*p2p4)
	mominv(6)=e1*(ame2-2d0*kp1)
	mominv(7)=e1*(ame2-2d0*kp2)
	mominv(8)=e1*ampi2
	mominv(9)=e1*0d0
	mominv(10)=e1*ampi2

	masses(0)=e1*ame2
	masses(1)=e1*0d0
	masses(2)=e1*0d0
	masses(3)=e1*ame2
	masses(4)=e1*ampi2
	if(instabilities.eq.1) then
	call routing_pipi(pt1,pt2,pt3,pt4,qt,
     >		     momvec,mominv,masses,4,3,elmat1,elmat2)
        routing1(4)=elmat1
        routing2(4)=elmat2
	endif

********************ROUTING 5 PER LE FULL TENSOR

	momvec(:,1)=e1*pt1
	momvec(:,2)=-e1*(pt2-qt)
	momvec(:,3)=e1*qt
	momvec(:,4)=e1*(pt1-pt3)

	mominv(1)=e1*ame2
	mominv(2)=e1*(2d0*ampi2+2d0*p3p4)
	mominv(3)=e1*ame2
	mominv(4)=e1*(ame2+ampi2-2d0*p2p4)
	mominv(5)=e1*(ame2+ampi2-2.d0*p1p3)
	mominv(6)=e1*(ame2-2d0*kp2)
	mominv(7)=e1*(ame2-2d0*kp1)
	mominv(8)=e1*ampi2
	mominv(9)=e1*0d0
	mominv(10)=e1*ampi2

	masses(0)=e1*ame2
	masses(1)=e1*0d0
	masses(2)=e1*0d0
	masses(3)=e1*ame2
	masses(4)=e1*ampi2


	call routing_pipi(pt1,pt2,pt3,pt4,qt,
     >		     momvec,mominv,masses,5,3,elmat1,elmat2)

        routing1(5)=elmat1
        routing2(5)=elmat2

c$$$        print*,'isr',routing1
c$$$        print*,'isr',routing2

c$$$        if(instabilities.eq.1) then
c$$$        call find_closest_two(routing1,elmat1)
c$$$
c$$$        call find_closest_two(routing2,elmat2)
c$$$
c$$$	endif


   	if(i.eq.2) then
     		elmat1=-elmat1
	endif
	virtual=virtual+elmat1+elmat2
        endif
	virtual=4d0*virtual*e6*alpha/(2d0*pi)
	end

c$$$  **************
c$$$

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!ROUTINE FOR FULL PENTAGON !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      subroutine VirtualIFI_Pent(pin1,pin2,p1,p2,q,virtual1,qphmat) 
#ifdef COLLIER
      use collier
#endif
      implicit none
      double complex virtual1,virtual,virtualpent,emtx
      double precision pin1(0:3),pin2(0:3),p1(0:3),p2(0:3),q(0:3)
      double precision qphmat(40,0:3)
      integer iffpi
      character*10 what_ffpi
      real*8 s_max
      common/ffpioption/iffpi,what_ffpi
      double complex ffpis
      double precision ffpis2,s

      external dot
      double precision dot
      common/resonancesgvmd/nres
      integer nres
      integer igvmdwarning
      data igvmdwarning /0/
      
      common/ialpharunning/iarun 
      integer iarun,iarunbck
      
      common/radoption/ifisrfsr
      integer ifisrfsr  

      common/instabilities/instabilities
      integer instabilities
      common/lt/ltint,quadevents,ilt
      integer ltint,quadevents,ilt

      double precision eikfull
      complex*32 quademtx
            

      virtual1=0.d0
      if(ltint.eq.0) then
         call virtualIFI_Pent_4pt_pr(pin1,pin2,p1,p2,q,
     >        virtual,qphmat)
         virtual1=virtual1+virtual
      endif
!       print*,virtual
      if(ltint.eq.1) then

         call virtualPent_pr_quad(pin1,pin2,p1,p2,q,
     >        virtual,qphmat,1)
!          print*,virtual
         virtual1=virtual1+virtual

         call virtualPent_pr_quad(pin1,pin2,p1,p2,q,
     >        virtual,qphmat,2)
!          print*,virtual
         virtual1=virtual1+virtual
         
      endif
      
      call virtualIFI_Pent_5pt_pr(pin1,pin2,p1,p2,q,
     >     virtual,qphmat,1)  
!       print*,virtual
      virtual1=virtual1+virtual
      call virtualIFI_Pent_5pt_pr(pin1,pin2,p1,p2,q,
     >     virtual,qphmat,2)
!       print*,virtual
      if(ltint.eq.1) ltint=0

      virtual1= virtual1+virtual
	end
**************

      subroutine virtual_CT_pr(p1,p2,p3,p4,q,virtual,qph)
#ifdef COLLIER
        use collier
#endif
	implicit none
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision qph(40,0:3)
	double precision virtual
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
	real*8 e2,e6
	double precision ame,ammu,convfac,alpha,pi,ame2,ampi2,ampi
	common/parameters/ame,ammu,convfac,alpha,pi
        common/hparam/ampi
	complex*16 elmat2,born
        double precision elmat1,elmat3
	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min
	common/softphotonmass/phmass
	real*8 phmass

	external dot
	real*8 dot

	double complex dphi
	double complex dmpi

	double complex dpsi
	double complex dme

	common/ialpharunning/iarun
        integer iarun,iarunbck

        common/radoption/ifisrfsr
        integer ifisrfsr

        double complex vpolc,alphar



********FF
	integer iffpi
	character*10 what_ffpi
	common/ffpioption/iffpi,what_ffpi
	double complex ffpiis
	double precision ffpiis2
	double complex ffpifs
	double precision ffpifs2


*********CLL STUFF
	double complex a1
	parameter (a1 = (1.d0,0.d0))
	double complex b1
	parameter (b1 = (1.d0,0.d0))

	Complex*16 A0ae
	Complex*16 A0ap

	Complex*16 B0ae0e
	Complex*16 B0ap0p

	Complex*16 DB0ae0e
	Complex*16 DB0ap0p
****





	e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2
	ame2=ame*ame
	ampi2=ampi*ampi

	p1p2 = dot(p1,p2)
	p1p3 = dot(p1,p3)
	p1p4 = dot(p1,p4)
	kp1  = dot(p1,q)
	p2p3 = dot(p2,p3)
	p2p4 = dot(p2,p4)
	kp2  = dot(p2,q)
	kp3 = dot(p3,q)
	p3p4 = dot(p3,p4)
	kp4  = dot(p4,q)

	x12  = 2d0*ame2+2.d0*p1p2
	x13  = ame2+ampi2-2.d0*p1p3
	x14  = ame2+ampi2-2.d0*p1p4
	x15  = ame2-2.d0*kp1
	x23  = ame2+ampi2-2d0*p2p3
	x24  = ame2+ampi2-2d0*p2p4
	x25 = ame2-2d0*kp2
	x34 = 2d0*ampi2+2d0*p3p4
	x35 = ampi2+2d0*kp3
	x45 = ampi2+2d0*kp4


	ffpifs=(1.d0,0.d0)
	ffpifs2=1.d0

	ffpiis=(1.d0,0.d0)
	ffpiis2=1.d0
	if(iffpi.gt.0) then
		call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2)
		call ffpi(2.d0*(ampi**2+p3p4),ffpiis,ffpiis2)
	ffpifs=CONJG(ffpifs)
	endif

******ALPHA RUNNING
	ffpiis=ffpiis*vpolc(2.d0*(ampi**2+p3p4))
	ffpiis2=ffpiis*CONJG(ffpiis)

	ffpifs=ffpifs*CONJG(vpolc(2.d0*(ame**2+p1p2)))
	ffpifs2=ffpifs*CONJG(ffpifs)
****************

#ifdef COLLIER
!	call SetDeltaUV_cll(1000000d0)
	call A0_cll(A0ap,a1*ampi2)

	call DB0_cll(DB0ap0p,b1*ampi2,b1*0d0,b1*ampi2)

	call B0_cll(B0ap0p,b1*ampi2,b1*0d0,b1*ampi2)

	call A0_cll(A0ae,a1*ame2)

	call DB0_cll(DB0ae0e,b1*ame2,b1*0d0,b1*ame2)

	call B0_cll(B0ae0e,b1*ame2,b1*0d0,b1*ame2)

#endif

	dpsi =(4d0*ame2**2*DB0ae0e-A0ae+ame2)/ame2

	dme=-(ame2*(2d0*B0ae0e-1d0)+A0ae)/ame2

	dphi=2d0*(B0ap0p+2d0*ampi2*DB0ap0p)

	dmpi=(A0ap-4d0*ampi2*B0ap0p)/ampi2

************PROVA

	call ISR_pr(p1,p2,p3,p4,q,elmat1,qph)
	call IFI_pr(p1,p2,p3,p4,q,elmat2,qph)
	call FSR_pr(p1,p2,p3,p4,q,elmat3,qph)

	if(ifisrfsr.eq.10) then

		virtual=real(dpsi*(elmat1))/(4d0*e6)
	elseif(ifisrfsr.eq.01) then

		virtual=real(dphi*(elmat3))/(4d0*e6)
	elseif(ifisrfsr.eq.11) then

	        virtual=real(dpsi*(elmat1+2.d0*elmat2+elmat3))/(4d0*e6)
     >		+real(dphi*(elmat1+2.d0*elmat2+elmat3)/(4d0*e6))
	endif

	include 'pipirad/Virtual/CT/CTmeISR.F'
	include 'pipirad/Virtual/CT/CTmeFSR.F'

	if(ifisrfsr.eq.10) then
		virtual=virtual+real(elmat2)*ffpiis2
	elseif(ifisrfsr.eq.01) then
		virtual=virtual
	elseif(ifisrfsr.eq.11) then

		virtual=virtual+real(elmat2)*ffpiis2
     >          +real(elmat1*ffpiis*ffpifs)
	endif


	include 'pipirad/Virtual/CT/CTmpiISR.F'
	include 'pipirad/Virtual/CT/CTmpiFSR.F'

	if(ifisrfsr.eq.10) then
		virtual=virtual
	elseif(ifisrfsr.eq.01) then
		virtual=virtual+real(elmat1)*ffpifs2
	elseif(ifisrfsr.eq.11) then

		virtual=virtual+real(elmat2*ffpiis*ffpifs)
     >		       +real(elmat1)*ffpifs2
	endif


!QUI HO IL PRIMO CONTROTERMINE, CHE NON È ALTRO CHE (DZPSI+DZPHI)*BORN^2

!	call bk_me(p1,p2,p3,p4,q,elmat1,qph)
!	virtual=real((dphi+dpsi)*(elmat1)/(4d0*e6))
!
!!QUI HO IL CONTROTERMINE DI MASSA DEGLI ELETTRONI
!	include 'pipirad/Virtual/CT/CTmeISR.F'
!	include 'pipirad/Virtual/CT/CTmeFSR.F'
!	virtual=virtual+real(elmat1+elmat2)

!
!
!!QUI HO IL CONTROTERMINE DI MASSA DEI PIONI

!	include 'pipirad/Virtual/CT/CTmpiISR.F'
!	include 'pipirad/Virtual/CT/CTmpiFSR.F'
!	virtual=virtual+real(elmat1+elmat2)
!
!
	virtual=4d0*virtual*e6*alpha/(2d0*pi)



	end
****************************
!-------------------------------------------------
!QUESTO È IL VIRTUALE COMPLETO DEI PIONI RADIATIVI
!-------------------------------------------------

      subroutine Virtual_pr(pin1,pin2,p1,p2,q,virtual1,qphmat)
#ifdef COLLIER
      use collier
#endif
      implicit none
      double precision virtual1,virtual,virtualpent,emtx
      double complex virtualc,virtualtbp
      double precision pin1(0:3),pin2(0:3),p1(0:3),p2(0:3),q(0:3)
      double precision qphmat(40,0:3)
      integer iffpi
      character*10 what_ffpi
      real*8 s_max
      common/ffpioption/iffpi,what_ffpi
      double complex ffpis
      double precision ffpis2,s

      external dot
      double precision dot
      common/resonancesgvmd/nres
      integer nres
      integer igvmdwarning
      data igvmdwarning /0/

	common/ialpharunning/iarun
        integer iarun,iarunbck

        common/radoption/ifisrfsr
        integer ifisrfsr

        common/instabilities/instabilities
        integer instabilities
        common/lt/ltint,quadevents,ilt
        integer ltint,quadevents,ilt

      double precision eikfull
      complex*32 quademtx

      
      common/colliermodedisp/initialisation
      integer initialisation

      common/compotontensor/compttens
      integer compttens

      double precision time1,time2
      
c$$$      if(initialisation.eq.0.and.iffpi.eq.1) then
c$$$         call SetMode_cll(2)
c$$$         initialisation=1
c$$$      endif

        virtual1=0d0

*--------------------------------------------*
*     This first part is not dependent on    *
*     the way the FF is introduced in the    *
*                    loop                    *
*--------------------------------------------*

		call virtualISR_pr(pin1,pin2,p1,p2,q,
     >                                       virtual,qphmat)
     		virtual1=virtual1+virtual

     		call virtualMixedISR_pr(pin1,pin2,p1,p2,q,
     >                                       virtual,qphmat)
                         virtual1=virtual1+virtual
!        print*,'AFTER ISR',virtual1
*--------------------------------------------*
*     This part depends on the way the FF    *
*     is introduced in the loop, but since   *
*     the loop is just on FS legs, the       *
*     virtual GVMD amplitude can be written  *
*     as the FxsQED+finite pieces, so we     *
*     still call the old FxsQED routines     *
*     and we evaluate finite pieces if       *
*                  iffpi=2                   *
*--------------------------------------------*
!      if(compttens.eq.3) then
      		call virtualFSR_pr(pin1,pin2,p1,p2,q,
     >                                       virtual,qphmat)

                     virtual1=virtual1+virtual
     		call virtualMixedFSR_pr(pin1,pin2,p1,p2,q,
     >                                       virtual,qphmat)
                     virtual1=virtual1+virtual

         call virtual_CT_pr(pin1,pin2,p1,p2,q,virtual,qphmat)
         virtual1=virtual1+virtual

!         print*,'AFTER FSR',virtual1
!     return
!      endif
*--------------------------------------------*
*     This part depends on the way the FF    *
*     is introduced in the loop. In this     *
*     case it is not true that the virtual   *
*     GVMD amplitude can be written as the   *
*     FxsQED+finite pieces, so we split      *
*     the routines, the first are used for   *
*     iffpi=1 (FxsQED), the latter for       *
*              iffpi=2  (GVMD)               *
*--------------------------------------------*
!        virtual1=0d0
         if(iffpi.eq.1.or.iffpi.eq.0) then

     		call virtualIFI_triang_pr(pin1,pin2,p1,p2,q,
     >          virtualc,qphmat) 
           virtual1=virtual1+real(virtualc)

     		call virtualIFI_Boxes_pr(pin1,pin2,p1,p2,q,
     >          virtualc,qphmat)
           virtual1 = virtual1 + real(virtualc)

           call VirtualIFI_Pent(pin1,pin2,p1,p2,q,virtualc,qphmat)
           virtual1 = virtual1 + real(virtualc)

        elseif(iffpi.eq.2) then

           call virtualIFI_GVMD(pin1,pin2,p1,p2,q,virtual,qphmat)
                virtual1=virtual1+virtual

        elseif(iffpi.eq.3) then

           call virtualIFI_disp(pin1,pin2,p1,p2,q,virtual,qphmat)
                     virtual1=virtual1+virtual
                endif
!      print*,'after ifi',virtual1
!        return
        end


      !===================================================!
      !  FsQED ROUTINES FOR IFI IN \PI\PI\GAMMA           !   MG 24/04/26
      !===================================================!

!===================================================!                                                                                         
!  ROUTINE FOR IFI WITH MASSES M1,M2 IN \PI\PI\GAMMA!                                                                                         
!===================================================!

      subroutine virtualIFI_ij_Box_disp(p1,p2,p3,p4,m12,m22
     >     ,q,elmatre,elmatim,qph,i)
c     This routines if i=1 returns Re(elmat(l,m)) and Im(elmat(l,m))
c     for l,m=1,2 that stands for IS,FS real photon in virtual-born
c     If i=2 we're evaluating the Im \int ds' D0, for D0_s behaving as
c     1/(s-eta) being eta=M_pp^2,s for IS and FS real photon
c     in virtual diagram. In the latter case Re(elmat(l,m))==0
c     If i=3 we're evaluating the Im(tilde{elmat(l,m)}) where
c     tilde stands for the fact that the ImD0 is subtracted by the
c     Im part behaving as in i=2
c     No FF included in this routine. 
#ifdef COLLIER 
        use collier
#endif
	implicit double complex (z)
!        implicit none
        double complex m12,m22,m1,m2
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
        double precision p1t(0:3),p2t(0:3),p3t(0:3),p4t(0:3),qt(0:3)
        double precision p1oz(0:3),p2oz(0:3),p3oz(0:3),p4oz(0:3)
        double precision qoz(0:3)
        double precision p1tz(0:3),p2tz(0:3),p3tz(0:3),p4tz(0:3)
        double precision qtz(0:3)
	double precision qph(40,0:3)
	double complex virtual,virtualpent,virtualpent1
        complex*32 virtual1,virtualq
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
	real*8 ep2,ep3,e2,e6,me2,me4,me6,mp2,mp4,mp6
        complex*32 me2q
	double precision ame,ammu,convfac,alpha,pi,ame2,ampi2,ampi
	common/parameters/ame,ammu,convfac,alpha,pi
        common/hparam/ampi
	Complex*16 elmat1,elmat2,elmat3,elmat4,elmat

	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min
	common/softphotonmass/phmass
	real*8 phmass
	
	external dot
	real*8 dot
	
        double complex vpolc,alphar
        
        common/ialpharunning/iarun 
        integer iarun 	

        
        complex*16 D0md,re
	complex*16 Dimd(3)
        complex*16 Dijmd(7)
        complex*16 E0md
        complex*16 Eimd(4)
        complex*16 Eijmd(11)
        complex*32 D0mq,rq
        complex*32 Dimq(3)
        complex*32 Dijmq(7)
        complex*16 C0md,C1md,C2md
**********FF

	integer iffpi
	character*10 what_ffpi
	common/ffpioption/iffpi,what_ffpi
	double complex ffpiis
	double precision ffpiis2
	double complex ffpifs
	double precision ffpifs2
        integer i,k
*************GVMD COMMON
        common/gvmdrad/indexnores
        integer indexnores
        
        common/gvmdchecks/iscoll,fscoll
        integer iscoll,fscoll
        common/checkform/intltgvmd
        integer intltgvmd
        
        common/compotontensor/compttens
        integer compttens
*--------FsQED
        real*8 oneoeps,mudimreg2,eulergamma
        common/IRdimreg/oneoeps,mudimreg2,eulergamma
        real*8 elmatre(2,2),elmatim(2,2)
        real*8 elmatres(2,2),elmatims(2,2)
        real*8 elmatresp(2,2),elmatimsp(2,2)
        integer l,m
        double complex setrez
        common/dispisfs/onlyis,onlyfs
        integer onlyis,onlyfs

*********CLL STUFF
        double complex cnpi(0:0,0:1,0:1)
        double complex cnpiuv(0:0,0:1,0:1)
        double complex dnpi(0:1,0:2,0:2,0:2)
	double complex dnpiuv(0:1,0:2,0:2,0:2)
        double complex enpi(0:1,0:2,0:2,0:2,0:2)
        double complex enpiuv(0:1,0:2,0:2,0:2,0:2)
	double complex c1
	parameter (c1 = (1.d0,0.d0))

**************  QUARTIC-VERTEX  PAVE *******************
!!! DIAGRAMS 1-4 IFI_TRIANG AND DIAGRAM 13 IFI_BOXES !!!

	double complex B0aeem1
	double complex B0aeem2
	double complex B0ax12m1m2
	double complex B0ax15em1
	double complex B0ax25em2
	double complex C0aeex12m1em2
	double complex C0aex15x34m2em1
	double complex B0ax34m1m2
	double complex C0aex25x34m1em2

        double complex C0ax15ex34m1em2
        double complex D0aex15ex250x34em1em2
        double complex D1aex15ex250x34em1em2
        double complex D3aex15ex250x34em1em2
        double complex D00aex15ex250x34em1em2
        double complex D11aex15ex250x34em1em2
        double complex D12aex15ex250x34em1em2
        double complex D13aex15ex250x34em1em2
        double complex D23aex15ex250x34em1em2
        double complex D33aex15ex250x34em1em2
        
**************  BOXES  PAVE *******************
!!!!!!!!!! DIAGRAMS 1-12 IFI_BOXES !!!!!!!!!!!!
        double complex C0ax12px45m2m1p
        double complex D0aepx45ex13x12em1pm2
        double complex D1aeex45px12x13m1em2p
        double complex D2aepx45ex13x12em1pm2
        double complex D3aepx45ex13x12em1pm2

        double complex C0ax12px45m1m2p
        double complex D0aex45pex23x12em1pm2
        double complex D1aeepx45x12x23m1em2p
        double complex D1aex45pex23x12em1pm2
        double complex D2aex45pex23x12em1pm2

        double complex C0ax12px35m2m1p
        double complex D0aex12x35x14epem1m2p
        double complex D1aeex35px12x14m1em2p
        double complex D2aex12x35x14epem1m2p
        double complex D3aex12x35x14epem1m2p

        double complex C0ax12px35m1m2p
        double complex D0aex12px24ex35em1m2p
        double complex D1aeepx35x12x24m1em2p
        double complex D1aex12px24ex35em1m2p
        double complex D3aex12px24ex35em1m2p

        double complex C0aex34x15em2m1
        double complex C0appx34m2pm1
        double complex C1apx34ppm1m2
        double complex C1apx34ppm2m1
        double complex D0ax15x34px23epem1m2p
        double complex D1ax15x34px23epem1m2p
        double complex D3ax15x34px23epem1m2p
        double complex D1aex23px34px15m2epm1

        double complex C0appx34m1pm2
        double complex D0ax15x34px24epem1m2p
        double complex D1aex24px34px15m2epm1
        double complex D1ax15x34px24epem1m2p
        double complex D3ax15x34px24epem1m2p

        double complex C0aex34x25em1m2
        double complex D0aex34px13x25pem1m2p
        double complex D1aex13px34px25m1epm2
        double complex D2aex34px13x25pem1m2p
        double complex D3aex34px13x25pem1m2p

        double complex D0aex34px14x25pem1m2p
        double complex D1aex14px34px25m1epm2
        double complex D2aex34px14x25pem1m2p
        double complex D3aex34px14x25pem1m2p
        
        double complex C1apx35x12m1pm2

        double complex C1apx35x12m2pm1

        double complex C1apx45x12m1pm2
        double complex D0aex12x45x13epem1m2p
        double complex D2aex12x45x13epem1m2p
        double complex D3aex12x45x13epem1m2p

        double complex C1apx45x12m2pm1
        double complex D0aex12px23ex45em1m2p
        double complex D1aex12px23ex45em1m2p
        double complex D3aex12px23ex45em1m2p

        double complex m22bck
        double complex C0aepx13em1p
        double complex C0aex13x45m2ep
        external Itot1ovsp_dimreg
        real*8 Itot1ovsp_dimreg

        external Itotspovsp_dimreg
        real*8 Itotspovsp_dimreg
        real*8 betap,cth,tridot
        external tridot
        common/momentadisp/momenta
        double precision momenta(4,0:3)
        
#ifdef COLLIER
        call GetDeltaUV_cll(deltauv)
        call GetMuUV2_cll(mu2dim)
        call GetMuIR2_cll(lambda2)
#endif
!     Initialise all loop functions to (0,0)
        include 'pipirad/gvmd/initial_ifi.f'
!     Get momenta p1,p2,p3,p4 in RF where \vec(p3+p4)=\vec(0) to use dispersive integral on IS radiation
        call new_boost(p3+p4,p1,p1t,1)
        call new_boost(p3+p4,p2,p2t,1)
        call new_boost(p3+p4,p3,p3t,1)
        call new_boost(p3+p4,p4,p4t,1)
        call new_boost(p3+p4,q,qt,1)

        momenta(1,:)=p1
        momenta(2,:)=p2
        momenta(3,:)=p3
        momenta(4,:)=p4
        
        mudimreg2= mu2dim
        eulergamma= 0.57721566490153286060d0
        oneoeps= log(lambda2/mudimreg2) + eulergamma - log(4.d0*pi)
                
        re =(1.d0,0.d0)
        virtual=(0.d0,0.d0)
        e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2
        ep2= 1.d0
        ep3= 1.d0
	ame2=ame*ame
	ampi2=ampi*ampi
        me2=ame2
        mp2=ampi2
        me4=me2**2
        me6=me2**3
        mp4=mp2**2
        mp6=mp2**3

        me2q = (ame*1.q0)*(ame*1.q0)
	p1p2 = dot(p1,p2)
	p1p3 = dot(p1,p3)
	p1p4 = dot(p1,p4)
	kp1  = dot(p1,q)
	p2p3 = dot(p2,p3)
	p2p4 = dot(p2,p4) 
	kp2  = dot(p2,q)
	kp3 = dot(p3,q)
	p3p4 = dot(p3,p4)
	kp4  = dot(p4,q)
					
	x12  = 2d0*ame2+2.d0*p1p2
	x13  = ame2+ampi2-2.d0*p1p3
	x14  = ame2+ampi2-2.d0*p1p4
	x15  = ame2-2.d0*kp1
	x23  = ame2+ampi2-2d0*p2p3
	x24  = ame2+ampi2-2d0*p2p4
	x25 = ame2-2d0*kp2
	x34 = 2d0*ampi2+2d0*p3p4
	x35 = ampi2+2d0*kp3
	x45 = ampi2+2d0*kp4

        m1  = sqrt(m12)
        m2  = sqrt(m22)

**************
	ffpifs=(1.d0,0.d0)
	ffpifs2=1.d0

	ffpiis=(1.d0,0.d0)
	ffpiis2=1.d0
	if(iffpi.gt.0) then
		call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2) 
		call ffpi(2.d0*(ampi**2+p3p4),ffpiis,ffpiis2)

                 endif

******ALPHA RUNNING	
	ffpiis=ffpiis*vpolc(2.d0*(ampi**2+p3p4))
	ffpiis2=ffpiis*CONJG(ffpiis)

	ffpifs=ffpifs*vpolc(2.d0*(ame**2+p1p2))
	ffpifs2=ffpifs*CONJG(ffpifs)
****************
        do l=1,2
           do m=1,2

              elmatre(l,m)=0.d0
              elmatim(l,m)=0.d0
           enddo
        enddo

#ifdef COLLIER

        call B0_cll(B0aeem1,c1*ame2,c1*ame2,c1*m12)

        call B0_cll(B0aeem2,c1*ame2,c1*ame2,c1*m22)

        call B0_cll(B0ax12m1m2,c1*x12,c1*m12,c1*m22)
        
        call B0_cll(B0ax15em1,c1*x15,c1*ame2,c1*m12)

        call B0_cll(B0ax25em2,c1*x25,c1*ame2,c1*m22)

        call B0_cll(B0ax34m1m2,c1*x34,c1*m12,c1*m22)
        
        call C0_cll(C0aeex12m1em2,c1*ame2,c1*ame2,c1*x12,
     >       c1*m12,c1*ame2,c1*m22)

        call C0_cll(C0aex15x34m2em1,c1*ame2,c1*x15,c1*x34,
     >       c1*m22,c1*ame2,c1*m12)

        call C0_cll(C0aex25x34m1em2,c1*ame2,c1*x25,c1*x34,
     >       c1*m12,c1*ame2,c1*m22)

        call C0_cll(C0ax15ex34m1em2,c1*x15,c1*ame2,c1*x34,
     >       c1*m12,c1*ame2,c1*m22)

        
#endif
        if(onlyis.eq.1) then
           call Dij_d(-p1,-q,+p2-q,me2*re,m12,me2*re,
     >          m22,D0md,Dimd,Dijmd)

           D0aex15ex250x34em1em2   = D0md
           D1aex15ex250x34em1em2   = Dimd(1)
           D3aex15ex250x34em1em2   = Dimd(3)
           D00aex15ex250x34em1em2  = Dijmd(1)
           D11aex15ex250x34em1em2  = Dijmd(2)
           D12aex15ex250x34em1em2  = Dijmd(3)
           D13aex15ex250x34em1em2  = Dijmd(4)
           D23aex15ex250x34em1em2  = Dijmd(6)
           D33aex15ex250x34em1em2  = Dijmd(7)
        endif
!     Se i=3 voglio soltanto la parte immaginaria
        if(i.eq.3) then
           B0aeem1     = setrez(B0aeem1)
           B0aeem2     = setrez(B0aeem2)
           B0ax12m1m2  = setrez(B0ax12m1m2)
           B0ax15em1   = setrez(B0ax15em1)
           B0ax25em2   = setrez(B0ax25em2)
           B0ax34m1m2  = setrez(B0ax34m1m2)
           C0aeex12m1em2 = setrez(C0aeex12m1em2)
           C0aex15x34m2em1 = setrez(C0aex15x34m2em1)
           C0aex25x34m1em2 = setrez(C0aex25x34m1em2)
           C0ax15ex34m1em2 = setrez(C0ax15ex34m1em2)
           D0aex15ex250x34em1em2   = setrez(D0md)
           D1aex15ex250x34em1em2   = setrez(Dimd(1))
           D3aex15ex250x34em1em2   = setrez(Dimd(3))
           D00aex15ex250x34em1em2  = setrez(Dijmd(1))
           D11aex15ex250x34em1em2  = setrez(Dijmd(2))
           D12aex15ex250x34em1em2  = setrez(Dijmd(3))
           D13aex15ex250x34em1em2  = setrez(Dijmd(4))
           D23aex15ex250x34em1em2  = setrez(Dijmd(6))
           D33aex15ex250x34em1em2  = setrez(Dijmd(7))
        endif
*-----When ep2=1 & ep3=0 we're looking at virtual diagrams where the------
*-----real photon is IS photon, when ep2=0 & ep3=1 the real is FS--------
        ep2=1.d0
        ep3=0.d0
        
        include 'pipirad/gvmd/triang_isr.f'
        elmatre(1,1) = elmatre(1,1)+real(elmat1)
        elmatim(1,1) = elmatim(1,1)+aimag(elmat1)
        
        include 'pipirad/gvmd/triang_fsr.f'
        elmatre(1,2) = elmatre(1,2)+real(elmat2)
        elmatim(1,2) = elmatim(1,2)+aimag(elmat2)

        ep2=0.d0
        ep3=1.d0

        include 'pipirad/gvmd/triang_isr.f'
        elmatre(2,1) = elmatre(2,1)+real(elmat1)
        elmatim(2,1) = elmatim(2,1)+aimag(elmat1)
        
        include 'pipirad/gvmd/triang_fsr.f'
        elmatre(2,2) = elmatre(2,2)+real(elmat2)
        elmatim(2,2) = elmatim(2,2)+aimag(elmat2)

        
**************BOXES  *******************
!!!!!!!!!!DIAGRAMS 1-12 IFI_BOXES !!!!!!!!!!!!
#ifdef COLLIER

        call C0_cll(C0ax12px45m2m1p,x12*c1,mp2*c1,x45*c1,m22*c1
     >       ,m12*c1,mp2*c1)

        call C0_cll(C0ax12px45m1m2p,x12*c1,mp2*c1,x45*c1,m12*c1
     >       ,m22*c1,mp2*c1)

        call C0_cll(C0ax12px35m2m1p,x12*c1,mp2*c1,x35*c1,m22*c1
     >       ,m12*c1,mp2*c1)

        call C0_cll(C0ax12px35m1m2p,x12*c1,mp2*c1,x35*c1,m12*c1
     >       ,m22*c1,mp2*c1)

        call C0_cll(C0aex34x15em2m1,c1*me2,c1*x34,c1*x15,me2*c1,
     >       c1*m22,c1*m12)

        call C0_cll(C0appx34m2pm1,mp2*c1,mp2*c1,x34*c1,m22*c1,
     >       mp2*c1,m12*c1)

        call C0_cll(C0appx34m1pm2,mp2*c1,mp2*c1,x34*c1,m12*c1,
     >       mp2*c1,m22*c1)

        call C0_cll(C0aex34x25em1m2,c1*me2,c1*x34,c1*x25,me2*c1,
     >       c1*m12,c1*m22)

#endif
        call Ci_d(-p4,+p3,mp2*re,m12,m22,C0md,C1apx34ppm1m2,C2md)
        call Ci_d(-p4,+p3,mp2*re,m22,m12,C0md,C1apx34ppm2m1,C2md)
        call Ci_d(p4,p1+p2,m12,mp2*re,m22,C0md,C1apx35x12m1pm2,
     >       C2md)
        call Ci_d(p4,p1+p2,m22,mp2*re,m12,C0md,C1apx35x12m2pm1,
     >       C2md)
        call Ci_d(p3,p1+p2,m12,mp2*re,m22,C0md,C1apx45x12m1pm2,
     >       C2md)
        call Ci_d(p3,p1+p2,m22,mp2*re,m12,C0md,C1apx45x12m2pm1,
     >       C2md)

        if(i.eq.3) then
           C0ax12px45m2m1p = setrez(C0ax12px45m2m1p)
           C0ax12px45m1m2p = setrez(C0ax12px45m1m2p)
           C0ax12px35m2m1p = setrez(C0ax12px35m2m1p)
           C0ax12px35m1m2p = setrez(C0ax12px35m1m2p)
           C0aex34x15em2m1 = setrez(C0aex34x15em2m1)
           C0appx34m2pm1   = setrez(C0appx34m2pm1)
           C0appx34m1pm2   = setrez(C0appx34m1pm2)
           C0aex34x25em1m2 = setrez(C0aex34x25em1m2)
           C1apx34ppm1m2   = setrez(C1apx34ppm1m2)
           C1apx34ppm2m1   = setrez(C1apx34ppm2m1)
           C1apx35x12m1pm2 = setrez(C1apx35x12m1pm2)
           C1apx35x12m2pm1 = setrez(C1apx35x12m2pm1)
           C1apx45x12m1pm2 = setrez(C1apx45x12m1pm2)
           C1apx45x12m2pm1 = setrez(C1apx45x12m2pm1)
        endif
!     If i=2 i have to take into account only Im D0_s propto 1/(eta-s')
!     that are these one below
        
        if(i.eq.2) then
           do l=1,2
              do m=1,2                 
                 elmatre(l,m)=0.d0
                 elmatim(l,m)=0.d0
              enddo
           enddo
           include 'pipirad/gvmd/initial_ifi.f'        
        endif
        if(onlyfs.eq.1) then
           call Di_disp(-p1,-p1+p3,p2,me2*re,m12,mp2*re,m22,
     >          x12,i,D0md,Dimd)
           D0aepx45ex13x12em1pm2  = D0md
           D2aepx45ex13x12em1pm2  = Dimd(2)
           D3aepx45ex13x12em1pm2  = Dimd(3)

           call Di_disp(p1,p2+p1,p3,m12,me2*re,m22,mp2*re,
     >          x12,i,D0md,Dimd)

           D1aeex45px12x13m1em2p  = Dimd(1)

!!!!!!!Diag 2
           call Di_disp(-p1,p2-p3,p2,me2*re,m22,mp2*re,m12,
     >          x12,i,D0md,Dimd)

           D0aex45pex23x12em1pm2  = D0md
           D1aex45pex23x12em1pm2  = Dimd(1)
           D2aex45pex23x12em1pm2  = Dimd(2)


           call Di_disp(p1,p2+p1,p4+q,m22,me2*re,m12,mp2*re,
     >          x12,i,D0md,Dimd)

           D1aeepx45x12x23m1em2p  = Dimd(1)

!!!!!!Diag 3
           call Di_disp(-p1,p2,p4-p1,me2*re,m12,m22,mp2*re,
     >          x12,i,D0md,Dimd)

           D0aex12x35x14epem1m2p  = D0md
           D2aex12x35x14epem1m2p  = Dimd(2)
           D3aex12x35x14epem1m2p  = Dimd(3)

           call Di_disp(p1,p2+p1,p4,m12,me2*re,m22,mp2*re,
     >          x12,i,D0md,Dimd)

           D1aeex35px12x14m1em2p  = Dimd(1)

!!!!!!Diag 4
           call Di_disp(-p1,p2,p2-p4,me2*re,m22,m12,mp2*re,
     >          x12,i,D0md,Dimd)
           
           D0aex12px24ex35em1m2p = D0md
           D1aex12px24ex35em1m2p = Dimd(1)
           D3aex12px24ex35em1m2p = Dimd(3)
           
           call Di_disp(p1,p1+p2,p3+q,m22,me2*re,m12,mp2*re,
     >          x12,i,D0md,Dimd)
           D1aeepx35x12x24m1em2p = Dimd(1)
        endif
        if(onlyis.eq.1) then
!!!!!!Diag 5
           call Di_disp(-p1t+qt,p2t,p2t-p3t,me2*re,m22,m12,mp2*re,
     >          x34,i,D0md,Dimd)
           D0ax15x34px23epem1m2p = D0md
           D1ax15x34px23epem1m2p = Dimd(1)
           D3ax15x34px23epem1m2p = Dimd(3)

           call Di_disp(-p2t,-p3t,-p3t-p4t,m12,me2*re,mp2*re,m22,
     >          x34,i,D0md,Dimd)
           D1aex23px34px15m2epm1 = Dimd(1)
           
!!!!!!Diag 6
           call Di_disp(-p1t+qt,p2t,p2t-p4t,me2*re,m22,m12,mp2*re,
     >          x34,i,D0md,Dimd)
           D0ax15x34px24epem1m2p = D0md
           D1ax15x34px24epem1m2p = Dimd(1)
           D3ax15x34px24epem1m2p = Dimd(3)

           call Di_disp(-p2t,-p4t,-p3t-p4t,m12,me2*re,mp2*re,m22,
     >          x34,i,D0md,Dimd)
           D1aex24px34px15m2epm1 = Dimd(1)

!!!!!!Diag 7
           call Di_disp(p1t,-p2t+qt,p1t-p3t,me2*re,m12,m22,mp2*re,
     >          x34,i,D0md,Dimd)
           D0aex34px13x25pem1m2p = D0md
           D2aex34px13x25pem1m2p = Dimd(2)
           D3aex34px13x25pem1m2p = Dimd(3)

           call Di_disp(-p1t,-p3t,-p3t-p4t,m12,me2*re,mp2*re,m22,
     >          x34,i,D0md,Dimd)
           D1aex13px34px25m1epm2 = Dimd(1)

!!!!!!Diag 8
           call Di_disp(p1t,-p2t+qt,p1t-p4t,me2*re,m12,m22,mp2*re,
     >          x34,i,D0md,Dimd)
           D0aex34px14x25pem1m2p = D0md
           D2aex34px14x25pem1m2p = Dimd(2)
           D3aex34px14x25pem1m2p = Dimd(3)

           call Di_disp(-p1t,-p4t,-p3t-p4t,m12,me2*re,mp2*re,m22,
     >          x34,i,D0md,Dimd)
           D1aex14px34px25m1epm2 = Dimd(1)
        endif
        
!!!!!!Diag 11
        D0aex12x45x13epem1m2p = D0aepx45ex13x12em1pm2
        D2aex12x45x13epem1m2p = D3aepx45ex13x12em1pm2
        D3aex12x45x13epem1m2p = D2aepx45ex13x12em1pm2

!!!!!!Diag 12
        D0aex12px23ex45em1m2p = D0aex45pex23x12em1pm2
        D1aex12px23ex45em1m2p = D1aex45pex23x12em1pm2
        D3aex12px23ex45em1m2p = D2aex45pex23x12em1pm2

        
*--------Same as before-------*
        ep2=1.d0
        ep3=0.d0
        
        include 'pipirad/gvmd/box_isr.f'
        elmatre(1,1) = elmatre(1,1)+real(elmat1)
        elmatim(1,1) = elmatim(1,1)+aimag(elmat1)
        
        include 'pipirad/gvmd/box_fsr.f'
        elmatre(1,2) = elmatre(1,2)+real(elmat2)
        elmatim(1,2) = elmatim(1,2)+aimag(elmat2)

        ep2=0.d0
        ep3=1.d0

        include 'pipirad/gvmd/box_isr.f'
        elmatre(2,1) = elmatre(2,1)+real(elmat1)
        elmatim(2,1) = elmatim(2,1)+aimag(elmat1)
!	print*,'NUOVA CALL',elmat1
        include 'pipirad/gvmd/box_fsr.f'
        elmatre(2,2) = elmatre(2,2)+real(elmat2)
        elmatim(2,2) = elmatim(2,2)+aimag(elmat2)


        elmatre = elmatre*4.d0*alpha*e6/(2.d0*pi)
        elmatim	= elmatim*4.d0*alpha*e6/(2.d0*pi)

        return
      end


!===================================================!                                                                                                                     
!       ROUTINE FOR FULL FsQED IFI IN \PI\PI\GAMMA  !                                                                                                               
!===================================================!
      subroutine virtualIFI_disp(p1,p2,p3,p4,q,emtx,qph)
      implicit none
      character*2 fs
      common/finalstate/fs
      real*8 ecms,ecmsnom
      common/ecms/ecms,ecmsnom
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision qph(40,0:3)
      double complex virtual,virtual1,bornc
      double precision emtx
      double precision elmat
      double complex ff(2,2),cff(2,2)
      complex*32 virtualq
      external tridot,dot
      double precision tridot,dot,eikfull
      double precision ame,ammu,convfac,alpha,pi,ame2,ampi2,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi

********FsQED
      double precision spmin,spmax
      double complex sp,ss,spp,mz,spbck
      double precision wsp,wspp,wbck
      double complex dispsp,dispss,dispz
      double complex ffpisp,ffpiss,ffpispp
      double complex ipart,rpart
      integer i,j,l,i1,i2
      double precision re(2,2),im(2,2)
      double precision res(2,2),ims(2,2)
      double precision resp(2,2),imsp(2,2)
*------ptlike amplitudes
      double precision rept(2,2),impt(2,2)
      double precision respt(2,2),imspt(2,2)
      double precision resppt(2,2),imsppt(2,2)
*------disp amplitudes
      double precision redisp(2,2),imdisp(2,2)
      double precision resdisp(2,2),imsdisp(2,2)
      double precision respdisp(2,2),imspdisp(2,2)
      double precision imdispIR(2,2),imdispfin(2,2)
*-----Needed for different subtraction
      double precision impents(2,2),impentsp(2,2),imbox(2,2)
*-----Double dispersive
      double precision redoubledisp(2,2),imdoubledisp(2,2)
*-----Sokhotski–Plemelj theorem
      double precision imresidual(2,2)
*-----PV integral
      double precision repv(2,2),impv(2,2)
*-----Common
      common/onlysp/onlysp
      integer onlysp            !     If onlysp=1, the computation of disp=1 in pent is done only with coefficient s',
                                !     otherwise if onlysp=0 then the computation of pent is done with coefficient s' and s (time consuming)
      common/dispisfs/onlyis,onlyfs
      integer onlyis,onlyfs
      double precision elmatre(2,2),elmatim(2,2)
      double precision elmatres(2,2),elmatims(2,2)
      double precision elmatresp(2,2),elmatimsp(2,2)
      double precision eta(2,2)
      integer disp
      double complex setrez
      external setrez
      double complex diff,sppv,ffpipv
      double precision delta,psi(2)
**********FF

      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi
      double complex ffpiis
      double precision ffpiis2
      double complex ffpifs
      double precision ffpifs2
      
      common/compotontensor/compttens
      integer compttens

      real*16 p1q(0:3),p2q(0:3),p3q(0:3),p4q(0:3),qq(0:3)
      common/momquad/p1q,p2q,p3q,p4q,qq

      double precision timein,timefin
      
      i1=1
      i2=2

      onlyis = 1
      onlyfs = 1
      
      if(compttens.eq.1) then
         i2=1
      elseif(compttens.eq.2) then
         i1=2
      endif

!      call getrnd(psi,2)
!!!   Needed for PaVe reduction
      call shellquad(p1,p2,p3,p4,q,p1q,p2q,p3q,p4q,qq)

      if(iffpi.gt.0) then
         call ffpi(2.d0*(ame**2+dot(p1,p2)),ffpifs,ffpifs2) 
         call ffpi(2.d0*(ampi**2+dot(p3,p4)),ffpiis,ffpiis2)
      endif

      
      cff(1,1)=conjg(ffpiis)
      cff(1,2)=conjg(ffpifs)
      cff(2,1)=conjg(ffpiis)
      cff(2,2)=conjg(ffpifs)

      ff(1,1) = ffpiis
      ff(1,2) = ffpiis
      ff(2,1) = ffpifs
      ff(2,2) = ffpifs

      eta(1,1) = dot(p3+p4,p3+p4)
      eta(1,2) = dot(p3+p4,p3+p4)
      eta(2,1) = dot(p1+p2,p1+p2)
      eta(2,2) = dot(p1+p2,p1+p2)


      call bk_me(p1,p2,p3,p4,q,emtx,qph)

      virtual   =  (0.d0,0.d0)
      virtualq  =  (0.d0,0.d0)
      if(iffpi.eq.3) then
         virtual=(0.d0,0.d0)
         rpart=(1.d0,0.d0)
         ipart=(0.d0,1.d0)
         redisp=0.d0
         imdisp=0.d0
         resp=(0.d0,0.d0)
         imsp=(0.d0,0.d0)
         mz=(0.d0,0.d0)
         
*------POINTLIKE NEEDED FOR IR BEHAVIOR OF RE RE---------
*     We need Re(ptlike*ff(eta))*real(C(ff(eta')))
         onlysp = 1
         
         disp=1
         call virtualIFI_ij_Box_disp(p1,p2,p3,p4,mz,mz,q,
     >        re,im,qph,disp)
         
         call VirtualIFI_ij_Pent_disp_full(p1,p2,p3,p4,q,mz,mz
     >     ,res,ims,resp,imsp,qph,disp)

         rept=re+resp
         impt=im+imsp
!!!!! Metto -ptlike per vedere che torni la IR cancellation
!         print*,cff
!         print*,ff
c$$$         do i=i1,i2
c$$$            do j=1,2
c$$$               virtual=virtual-real(
c$$$     >              (rept(i,j)+ipart*impt(i,j))*cff(i,j)*ff(i,j))
c$$$            enddo
c$$$         enddo
c$$$         print*,'PTLIKE',-virtual
c$$$         
c$$$         virtual =(0.d0,0.d0)
*     IR part coming from ReRe (not to be weighted MC)
         do i=i1,i2
            do j=1,2
               virtual=virtual
     >              +rept(i,j)*real(ff(i,j))*real(cff(i,j))
            enddo
         enddo
*---- Generate s' and evaluate F_pi(s')
         spmin = 4.d0*ampi**2
         spmax = 9.d0 
         if(spmax.le.ecmsnom**2*1.5d0) spmax=ecmsnom**2*1.5d0
     
!         call sp_multich(spmin,spmax,aimag(ff(1,1)),sp,wsp)
!         call sp_oosp(spmin,spmax,sp,wsp)
         call sp_pv(spmin,spmax,eta(1,1),sp,wsp)
         call ffpi(real(sp),ffpisp,ffpifs2) 
         spbck=sp
         wbck=wsp
         
*---- DISPERSIVE s'
         onlysp=0
        
         disp=1
         call virtualIFI_ij_Box_disp(p1,p2,p3,p4,mz,sp,q,
     >        re,im,qph,disp)
         

         call VirtualIFI_ij_Pent_disp_full(p1,p2,p3,p4,q,mz,sp
     >        ,res,ims,resp,imsp,qph,disp)

         redisp = redisp +re +resp
         imdisp = imdisp +im +imsp

         impents   = ims
         impentsp =  imsp
         imbox    = im
         
         call virtualIFI_ij_Box_disp(p1,p2,p3,p4,sp,mz,q,
     >        re,im,qph,disp)
         

         call VirtualIFI_ij_Pent_disp_full(p1,p2,p3,p4,q,sp,mz
     >        ,res,ims,resp,imsp,qph,disp)
	
         redisp	= redisp +re +resp
         imdisp = imdisp +im +imsp
         impents  = impents + ims
         impentsp = impentsp + imsp
         imbox  = imbox  + im

*-------------------------------------
*------------REAL PART----------------
*-------------------------------------
         
*     REAL PART (CHECK THAT DOES NOT DEPEND ON LAMBDA)
*     We need amp(s',lambda^2)+amp(lambda^2,s')+eta/(s'-eta)*amp(lambda^2,lambda^2)
*-----------------------------------
*---------------------------------         
*     PV integral
!         virtual=0.d0
!         print*,'SP',sp
!         print*,'virtual pre pv',virtual
         do i=i1,i2
            if(i.eq.1) then
               onlyis = 1
               onlyfs = 0
            elseif(i.eq.2) then
               onlyis =	0
               onlyfs =	1
            endif
            diff = sp -eta(i,1)*rpart
            sppv = eta(i,1)*rpart-diff
            call ffpi(real(sppv),ffpipv,ffpifs2)
!            print*,'SP, SPPV',sp,sppv
!            print*,aimag(ffpipv),aimag(ffpisp)
*     Check if the simmetric point is into the dispersive integration region
            if((real(sppv).le.spmax).and.(real(sppv).ge.spmin)) then
!               print*,'CI E IL PV'
               onlysp = 1

               disp=1
               call virtualIFI_ij_Box_disp(p1,p2,p3,p4,mz,sppv,q,
     >              re,im,qph,disp)               

               call VirtualIFI_ij_Pent_disp_full(p1,p2,p3,p4,q,mz,
     >              sppv,res,ims,resp,imsp,qph,disp)

               repv = +re +resp
               impv = im + imsp
               call virtualIFI_ij_Box_disp(p1,p2,p3,p4,sppv,mz,q,
     >              re,im,qph,disp)
               
               call VirtualIFI_ij_Pent_disp_full(p1,p2,p3,p4,q,sppv,
     >              mz,res,ims,resp,imsp,qph,disp)
               
               repv	= repv +re +resp
               impv     = impv +im +imsp
*     Must take into account that we're generating s' over all range and adding the simmetric, these contribution must be divided by 2
               do j=1,2
                  virtual=virtual-1.d0/pi*aimag(ffpipv)*
     >                 (repv(i,j)+eta(i,j)/(sppv-eta(i,j))*rept(i,j))
     >                 *real(cff(i,j))/wsp/sppv/2.d0
                  virtual=virtual-1.d0/pi*aimag(ffpisp)*
     >                 (redisp(i,j)+eta(i,j)/(sp-eta(i,j))*rept(i,j))
     >                 *real(cff(i,j))/wsp/sp/2.d0
!                  print*,'i,j,virtual',i,j,virtual
               enddo
            else
               do j=1,2
                  virtual=virtual-1.d0/pi*aimag(ffpisp)*
     >                 (redisp(i,j)+eta(i,j)/(sp-eta(i,j))*rept(i,j))
     >                 *real(cff(i,j))/wsp/sp
               enddo
            endif
!            endif
         enddo
         onlyis = 1
         onlyfs = 1
!         print*,'after pv-real',virtual
  
*     Here we have to carefully take into account for the Sokhotski-Plemelj
         do i=i1,i2
*     In the subtracted Real amplitude
*     amp(s',lambda^2)+amp(lambda^2,s')+eta/(s'-eta)*amp(lambda^2,lambda^2)
            
            if(i.eq.1) then
               onlyis = 1
               onlyfs = 0
            elseif(i.eq.2) then
               onlyis =	0
               onlyfs =	1
            endif

            sp=eta(i,1)*(1.d0-1.d0/10.d0**8)*rpart
            onlysp = 1
            
            disp=1
            call virtualIFI_ij_Box_disp(p1,p2,p3,p4,mz,sp,q,
     >           re,im,qph,disp)
            

            call VirtualIFI_ij_Pent_disp_full(p1,p2,p3,p4,q,mz,sp
     >           ,res,ims,resp,imsp,qph,disp)

            imresidual = im + imsp
            
            call virtualIFI_ij_Box_disp(p1,p2,p3,p4,sp,mz,q,
     >           re,im,qph,disp)
            

            call VirtualIFI_ij_Pent_disp_full(p1,p2,p3,p4,q,sp,mz
     >           ,res,ims,resp,imsp,qph,disp)

            imresidual = imresidual + im +imsp

            do j=1,2
*     The dispersive amplitudes has finite limit just on the left, so we have a first term            
*     -1/pi*(pi/2*f(eta-)) where f(s')=Im(f(s'))/s'*Im(amp(s',lambda^2)+amp(lambda^2,s'))*(eta-s')

               virtual=virtual-1.d0/2.d0*aimag(ff(i,j))/eta(i,j)*
     >              imresidual(i,j)*(eta(i,j)-sp)*real(cff(i,j))

            enddo
         enddo
         onlyis = 1
         onlyfs = 1
         sp=spbck
         
         
         
*-------------------------------------         
*---------IMAGINARY PART--------------
*-------------------------------------
*     Here we have to underline that we ALWAYS have to call -CONJG(im(ff(eta')) due to the fact
*     that i*i=-1
         
!         virtual=0.d0
*     We start from the pointlike imaginary part
         do i=i1,i2
            do j=1,2
               virtual=virtual+impt(i,j)*(-aimag(cff(i,j)))
            enddo
         enddo
!         print*,'ptlikeim',virtual
         
!     write(12, *) real(sp), real(virtual)
!         virtual=(0.d0,0.d0)
*     We evaluate the subtracted one. We have already evaluated all the quantities needed

         do i=i1,i2
            do j=1,2
               virtual=virtual-1.d0/pi*(
     >              (aimag(ffpisp)-aimag(ff(i,j)))*imbox(i,j)                 !box
     >              +aimag(ffpisp)*impentsp(i,j)                              !pent 
     >              -aimag(ff(i,j))*impents(i,j)                              !pent
     >              )/real(sp)/wsp*(-aimag(cff(i,j)))
            enddo
         enddo

!         print*,'subim',virtual
*     We evaluate the IR part originating from the imaginary part
         disp=2
         call virtualIFI_ij_Box_disp(p1,p2,p3,p4,mz,sp,q,
     >        re,im,qph,disp)
         
         call VirtualIFI_ij_Pent_disp_full(p1,p2,p3,p4,q,mz,sp
     >        ,res,ims,resp,imsp,qph,disp)

         imdispIR = +im +ims
         
         call virtualIFI_ij_Box_disp(p1,p2,p3,p4,sp,mz,q,
     >        re,im,qph,disp)

         call VirtualIFI_ij_Pent_disp_full(p1,p2,p3,p4,q,sp,mz
     >        ,res,ims,resp,imsp,qph,disp)
	
         imdispIR = imdispIR +im +ims
         
         do i=i1,i2
            do j=1,2
               virtual=virtual-1.d0/pi*(imdispIR(i,j)*aimag(ff(i,j))
     >              *(-aimag(cff(i,j))))
            enddo
         enddo
         
!         print*,'imir',virtual
         
!         write(14, *) real(sp), real(virtual)
!         virtual=(0.d0,0.d0)
*     We evaluate the remaining Dispersive integral over s' of all other PaVe
         disp=3
         call virtualIFI_ij_Box_disp(p1,p2,p3,p4,mz,sp,q,
     >        re,im,qph,disp)
         
         call VirtualIFI_ij_Pent_disp_full(p1,p2,p3,p4,q,mz,sp
     >        ,res,ims,resp,imsp,qph,disp)

         imdispfin = +im +ims
         
         call virtualIFI_ij_Box_disp(p1,p2,p3,p4,sp,mz,q,
     >        re,im,qph,disp)

         call VirtualIFI_ij_Pent_disp_full(p1,p2,p3,p4,q,sp,mz
     >        ,res,ims,resp,imsp,qph,disp)
	
         imdispfin = imdispfin +im +ims

         do i=i1,i2
            do j=1,2
               virtual=virtual-1.d0/pi*(imdispfin(i,j)*aimag(ff(i,j))
     >              *(-aimag(cff(i,j))))/real(sp)/wsp
            enddo
         enddo
!         print*,'after adding piece disp',virtual
         
!     emtx=real(virtual)
!         return
*---- Double dispersive------
!         virtual=0.d0
*---- Generate s'' and evaluate F_pi(s'')
         call cpu_time(timein)
         call sp_disp(spmin,spmax,sp,wsp)
         call ffpi(real(sp),ffpisp,ffpifs2)
         call sp_disp(spmin,spmax,spp,wspp)
         call ffpi(real(spp),ffpispp,ffpifs2) 
         onlysp = 1
         
         disp=1
         call virtualIFI_ij_Box_disp(p1,p2,p3,p4,sp,spp,q,
     >        re,im,qph,disp)
         

         call VirtualIFI_ij_Pent_disp_full(p1,p2,p3,p4,q,sp,spp
     >        ,res,ims,resp,imsp,qph,disp)

         redoubledisp = +re +resp
         imdoubledisp = +im +imsp
                  
         do i=i1,i2
            do j=1,2
               virtual=virtual+1.d0/(pi**2)*aimag(ffpisp)*aimag(ffpispp)
     >              *real((redoubledisp(i,j)+ipart*imdoubledisp(i,j))
     >              *cff(i,j))/real(sp)/real(spp)/wsp/wspp
            enddo
         enddo
         
!     print*,'double disp',virtual
!         stop
         emtx= real(virtual)
      endif

!     print*,'dd',virtual
      return
      end
      

!     Questa routine fa la doppia chiamata ai pent (direct e cross)
      subroutine VirtualIFI_ij_Pent_disp_full(p1,p2,p3,p4,q,m12,m22
     >     ,elmatres,elmatims,elmatresp,elmatimsp,qph,disp)
      implicit none
      double complex m12,m22
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision qph(40,0:3)
      integer disp
      real*8 elmatres1(2,2),elmatims1(2,2)
      real*8 elmatresp1(2,2),elmatimsp1(2,2)
      real*8 elmatres(2,2),elmatims(2,2)
      real*8 elmatresp(2,2),elmatimsp(2,2)
      integer l,m
      
      call VirtualIFI_ij_Pent_disp(p1,p2,p3,p4,q,m12,m22
     >     ,elmatres,elmatims,elmatresp,elmatimsp,qph,1,disp)
      
      call VirtualIFI_ij_Pent_disp(p1,p2,p3,p4,q,m12,m22
     >     ,elmatres1,elmatims1,elmatresp1,elmatimsp1,qph,2,disp)

      do l=1,2
         do m=1,2
            elmatres(l,m)=elmatres(l,m)+elmatres1(l,m)
            elmatresp(l,m)=elmatresp(l,m)+elmatresp1(l,m)
            elmatims(l,m)=elmatims(l,m)+elmatims1(l,m)
            elmatimsp(l,m)=elmatimsp(l,m)+elmatimsp1(l,m)
         enddo
      enddo
      
      return
      end


      subroutine VirtualIFI_ij_Pent_disp(p1,p2,p3,p4,q,m12,m22
     >     ,elmatres,elmatims,elmatresp,elmatimsp,qph,ind,disp) 
#ifdef COLLIER
      use collier
#endif
      implicit double complex (z)
!      implicit none
      double complex m12,m22,m1,m2
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision p1f(0:3),p2f(0:3),p3f(0:3),p4f(0:3),qf(0:3)
      double precision qph(40,0:3),p3t(0:3),p4t(0:3)
      double complex virtual,virtualpent
      integer ind,disp,dispbck,i      
      complex*32 virtual1
      real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
      real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
      real*8 e2,e6,me2,me4,me6,mp2,mp4,mp6
      complex*32 me2q
      double precision ame,ammu,convfac,alpha,pi,ame2,ampi2,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      Complex*16 elmat1,elmat2,elmat3,elmat4,elmat

      real*8 deltauv,lambda2,mu2dim,m2min
      common/loop/deltauv,lambda2,mu2dim,m2min
      common/softphotonmass/phmass
      real*8 phmass
      
      external dot
      real*8 dot
      
      double complex vpolc,alphar
      
      common/ialpharunning/iarun 
      integer iarun
      
      common/onlysp/onlysp
      integer onlysp
      
      complex*16 D0md,re
      complex*32 rq
      complex*16 Dimd(3)
      complex*16 Dijmd(7)
      complex*16 C0md,C1md,C2md
      complex*16 E0md
      complex*16 Eimd(4)
      complex*16 Eijmd(11)
      complex*32 E0mq
      complex*32 Eimq(4)
      complex*32 Eijmq(11)

      complex*32 D0mq
      complex*32 Dimq(3)
**********FF

      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi
      double complex ffpiis
      double precision ffpiis2
      double complex ffpifs
      double precision ffpifs2
      real*8 elmatres(2,2),elmatims(2,2)
      real*8 elmatresp(2,2),elmatimsp(2,2)
      integer l,m
      double complex setrez
      double precision masses(2)
      
*************GVMD COMMON
      common/gvmdrad/indexnores
      integer indexnores
      
      common/gvmdchecks/iscoll,fscoll
      integer iscoll,fscoll

      
      common/compotontensor/compttens
      integer compttens

      double precision ep2,ep3
*********CLL STUFF
      double complex cnpi(0:0,0:1,0:1)
      double complex cnpiuv(0:0,0:1,0:1)
      double complex dnpi(0:1,0:2,0:2,0:2)
      double complex dnpiuv(0:1,0:2,0:2,0:2)
      double complex enpi(0:1,0:2,0:2,0:2,0:2)
      double complex enpiuv(0:1,0:2,0:2,0:2,0:2)
      double complex c1
      parameter (c1 = (1.d0,0.d0))

**************PENT  PAVE *******************
!!!!!!!!!!DIAGRAMS 1 IFI_PENT !!!!!!!!!!!!

      double complex D0ap0px12x35x45m1ppm2
      double complex D1ap0px12x35x45m1ppm2
      double complex D3ap0px12x35x45m1ppm2
      double complex D1ax350x45x12ppm1ppm2
      double complex D00ap0px12x35x45m1ppm2
      double complex D11ap0px12x35x45m1ppm2
      double complex D12ap0px12x35x45m1ppm2
      double complex D13ap0px12x35x45m1ppm2
      double complex D11ax350x45x12ppm1ppm2
      double complex D13ax350x45x12ppm1ppm2

      double complex E3aeep0px12x24x45x35x13m1em2pp
      double complex E0aeex450x35x12x13ppx24m1em2pp
      double complex E1aeex450x35x12x13ppx24m1em2pp
      double complex E2aeex450x35x12x13ppx24m1em2pp
      double complex E3aeex450x35x12x13ppx24m1em2pp
      double complex E13aeep0px12x24x45x35x13m1em2pp
      double complex E23aeep0px12x24x45x35x13m1em2pp
      double complex E33aeep0px12x24x45x35x13m1em2pp
      double complex E00aeex450x35x12x13ppx24m1em2pp
      double complex E11aeex450x35x12x13ppx24m1em2pp
      double complex E12aeex450x35x12x13ppx24m1em2pp
      double complex E22aeex450x35x12x13ppx24m1em2pp
      double complex E13aeex450x35x12x13ppx24m1em2pp
      double complex E23aeex450x35x12x13ppx24m1em2pp
      double complex E33aeex450x35x12x13ppx24m1em2pp
      double complex E34aeex450x35x12x13ppx24m1em2pp

**************PENT  PAVE *******************                                                                                                                                                                                                           
!!!!!!!!!!DIAGRAMS 3 IFI_PENT !!!!!!!!!!!!                                                                                                                                                                                                             

      double complex E4ae0x24px34x15x13epx25m2eepm1
      double complex E3ae0x24px34x15x13epx25m2eepm1
      double complex E1ax150x13px34ex24x25pem2eepm1
      double complex E1ae0x24px34x15x13epx25m2eepm1
      double complex E14ax150x13px34ex24x25pem2eepm1 
      double complex E14ae0x24px34x15x13epx25m2eepm1
      double complex E13ax150x13px34ex24x25pem2eepm1 
      double complex E13ae0x24px34x15x13epx25m2eepm1
      double complex E12ae0x24px34x15x13epx25m2eepm1 
      double complex E11ax150x13px34ex24x25pem2eepm1
      double complex E11ae0x24px34x15x13epx25m2eepm1
      double complex E0ae0x24px34x15x13epx25m2eepm1
      double complex E00ae0x24px34x15x13epx25m2eepm1
      double complex E44ae0x24px34x15x13epx25m2eepm1
      double complex E34ae0x24px34x15x13epx25m2eepm1
      double complex E33ae0x24px34x15x13epx25m2eepm1 

      double complex D0apex15px24x34pm1em2
      double complex D3apex15px24x34pm1em2
      double complex D2apex15px24x34pm1em2
      double complex D1apex15px24x34pm1em2
      
      double complex D0apx25epx13x34pm1em2
      double complex D2apx25epx13x34pm1em2
      double complex D1apx25epx13x34pm1em2      
      double complex D1appex25x34x13m1pm2e
      
      double complex D1aex25ex150x34em1em2
      double complex D0a0ex34ex15x25eem2m1
      double complex D1a0ex34ex25x15eem1m2
      double complex D1a0ex34ex15x25eem2m1
      
      double complex D1aepx130x24x25em1pe
      double complex D1a0x13pex24x25eepm1      
      double complex D1a0epx13x25x24eem1p
      double complex D0a0epx13x25x24eem1p

      double complex C0ax24peepm1
      double complex C0ax13px25epm1
      double complex C0appx34m2pm1

      common/dispisfs/onlyis,onlyfs
      integer onlyis,onlyfs
     
      C0ax24peepm1   = (0.d0,0.d0)
      C0ax13px25epm1 = (0.d0,0.d0)
      C0appx34m2pm1  = (0.d0,0.d0)
      
      masses(1)  = real(m12)
      masses(2)  = real(m22)

      re =(1.d0,0.d0)

      e2 = 4.d0*pi*alpha
      e6 = e2*e2*e2
      ame2=ame*ame
      ampi2=ampi*ampi
      me2=ame2
      mp2=ampi2
      me4=me2**2
      me6=me2**3
      mp4=mp2**2
      mp6=mp2**3

      if(ind.eq.1) then
         p3t=p3
         p4t=p4
      elseif(ind.eq.2) then
         p3t=p4
         p4t=p3
      endif

      call new_boost(p3+p4,p1,p1f,1)
      call new_boost(p3+p4,p2,p2f,1)
      call new_boost(p3+p4,p3t,p3f,1)
      call new_boost(p3+p4,p4t,p4f,1)
      call new_boost(p3+p4,q,qf,1)

      me2q = (ame*1.q0)*(ame*1.q0)
      p1p2 = dot(p1,p2)
      p1p3 = dot(p1,p3t)
      p1p4 = dot(p1,p4t)
      kp1  = dot(p1,q)
      p2p3 = dot(p2,p3t)
      p2p4 = dot(p2,p4t) 
      kp2  = dot(p2,q)
      kp3 = dot(p3t,q)
      p3p4 = dot(p3t,p4t)
      kp4  = dot(p4t,q)
      
      x12  = 2d0*ame2+2.d0*p1p2
      x13  = ame2+ampi2-2.d0*p1p3
      x14  = ame2+ampi2-2.d0*p1p4
      x15  = ame2-2.d0*kp1
      x23  = ame2+ampi2-2d0*p2p3
      x24  = ame2+ampi2-2d0*p2p4
      x25 = ame2-2d0*kp2
      x34 = 2d0*ampi2+2d0*p3p4
      x35 = ampi2+2d0*kp3
      x45 = ampi2+2d0*kp4

      m1  = sqrt(m12)
      m2  = sqrt(m22)

      elmatres=0.d0
      elmatims=0.d0
      elmatresp=0.d0
      elmatimsp=0.d0
!!!!!!!!!!!!!!!PENT 1(FSR)!!!!!!!!!!!!!!!!!
      if(onlyfs.eq.1) then
         if(disp.eq.1) then
            
            call Eij_disp(p1,p1+p2,p3t,p3t+q,m12*re,me2*re,m22*re,mp2*re
     >           ,mp2*re,x12,disp,E0md,Eimd,Eijmd)
            if(fscoll.eq.1) then
               call Eij_dispq(p1,p1+p2,p3t,p3t+q,m12*re,me2*re,m22*re,
     >              mp2*re,mp2*re,x12,disp,E0mq,Eimq,Eijmq)
               E0md  = E0mq
               Eimd  = Eimq
               Eijmd = Eijmq
            endif
            E0aeex450x35x12x13ppx24m1em2pp   = E0md
            E1aeex450x35x12x13ppx24m1em2pp   = Eimd(1)
            E2aeex450x35x12x13ppx24m1em2pp   = Eimd(2)
            E3aeex450x35x12x13ppx24m1em2pp   = Eimd(3)
            E3aeep0px12x24x45x35x13m1em2pp   = Eimd(4)
            E00aeex450x35x12x13ppx24m1em2pp  = Eijmd(1)
            E11aeex450x35x12x13ppx24m1em2pp  = Eijmd(2)
            E12aeex450x35x12x13ppx24m1em2pp  = Eijmd(3)
            E13aeex450x35x12x13ppx24m1em2pp  = Eijmd(4)
            E13aeep0px12x24x45x35x13m1em2pp  = Eijmd(5)
            E22aeex450x35x12x13ppx24m1em2pp  = Eijmd(6)
            E23aeex450x35x12x13ppx24m1em2pp  = Eijmd(7)
            E23aeep0px12x24x45x35x13m1em2pp  = Eijmd(8)
            E33aeex450x35x12x13ppx24m1em2pp  = Eijmd(9)
            E34aeex450x35x12x13ppx24m1em2pp  = Eijmd(10)
            E33aeep0px12x24x45x35x13m1em2pp  = Eijmd(11)

!     Solo se disp=1 mi serve il tutto con i coefficienti s'
            include 'pipirad/gvmd/pent1_isr.f'
            elmatresp(2,1) = elmatresp(2,1)+real(elmat1)
            elmatimsp(2,1) = elmatimsp(2,1)+aimag(elmat1)

            include 'pipirad/gvmd/pent1_fsr.f'

            if(ind.eq.2) then
               elmat2=-elmat2
            endif
            elmatresp(2,2) = elmatresp(2,2)+real(elmat2)
            elmatimsp(2,2) = elmatimsp(2,2)+aimag(elmat2)
         endif
!     print*,'DISP ORIGINALE',disp
         dispbck=disp

         if((disp.eq.1).and.onlysp.eq.1) then
            continue
         else
            if(disp.eq.1) disp=disp+3
            
!     print*,'DISP CAMBIA SOLO SE DISP=1',disp
!     Qui ora valuto la i coefficienti in \eta e non in s' ( sia in penti_isr/fsr.f
!     che nella riduzione delle Eij)
            
            call Eij_disp(p1,p1+p2,p3t,p3t+q,m12*re,me2*re,m22*re,mp2*re
     >           ,mp2*re,x12,disp,E0md,Eimd,Eijmd)
            if(fscoll.eq.1) then
               call Eij_dispq(p1,p1+p2,p3t,p3t+q,m12*re,me2*re,m22*re,
     >              mp2*re,mp2*re,x12,disp,E0mq,Eimq,Eijmq)
               E0md  = E0mq
               Eimd  = Eimq
               Eijmd = Eijmq
            endif

!     print*,E0md
!     print*,Eimd
            E0aeex450x35x12x13ppx24m1em2pp   = E0md
            E1aeex450x35x12x13ppx24m1em2pp   = Eimd(1)
            E2aeex450x35x12x13ppx24m1em2pp   = Eimd(2)
            E3aeex450x35x12x13ppx24m1em2pp   = Eimd(3)
            E3aeep0px12x24x45x35x13m1em2pp   = Eimd(4)
            E00aeex450x35x12x13ppx24m1em2pp  = Eijmd(1)
            E11aeex450x35x12x13ppx24m1em2pp  = Eijmd(2)
            E12aeex450x35x12x13ppx24m1em2pp  = Eijmd(3)
            E13aeex450x35x12x13ppx24m1em2pp  = Eijmd(4)
            E13aeep0px12x24x45x35x13m1em2pp  = Eijmd(5)
            E22aeex450x35x12x13ppx24m1em2pp  = Eijmd(6)
            E23aeex450x35x12x13ppx24m1em2pp  = Eijmd(7)
            E23aeep0px12x24x45x35x13m1em2pp  = Eijmd(8)
            E33aeex450x35x12x13ppx24m1em2pp  = Eijmd(9)
            E34aeex450x35x12x13ppx24m1em2pp  = Eijmd(10)
            E33aeep0px12x24x45x35x13m1em2pp  = Eijmd(11)

            do i=1,2
               if (masses(i).ge.4.d0*ampi**2) then
                  if (i.eq.1) m12 = x12*re
                  if (i.eq.2) m22 = x12*re
               endif
            enddo

            include 'pipirad/gvmd/pent1_isr.f'
            elmatres(2,1) = elmatres(2,1)+real(elmat1)
            elmatims(2,1) = elmatims(2,1)+aimag(elmat1)

            include 'pipirad/gvmd/pent1_fsr.f'
            if(ind.eq.2) then
               elmat2=-elmat2
            endif
            elmatres(2,2) = elmatres(2,2)+real(elmat2)
            elmatims(2,2) = elmatims(2,2)+aimag(elmat2)
            
            m12  =  masses(1)*(1.d0,0.d0)
            m22  =  masses(2)*(1.d0,0.d0)

            disp = dispbck
         endif
      endif
!!!!!!!!!!!!!!!!!PENT 3 (ISR)!!!!!!!!!!!!!!!!!
!!!!  MOMENTI NEL RF FS RIPOSO  !!!!!!!!!!!!!!
!     print*,'DISP DEVE ESSERE TORNATO ORIGINALE',disp
c     Le C0 non danno i poli 1/(eta-s'), quindi per i=2 vanno spente
      if(onlyis.eq.1) then
         if(disp.eq.1) then
            
#ifdef COLLIER

            call C0_cll(C0ax24peepm1,c1*x24,c1*mp2,c1*me2,
     >           c1*me2,c1*mp2,c1*m12)

            call C0_cll(C0ax13px25epm1,c1*x13,c1*mp2,c1*x25,
     >           c1*me2,c1*mp2,c1*m12)

            call C0_cll(C0appx34m2pm1,c1*mp2,c1*mp2,c1*x34,
     >           c1*m22,c1*mp2,c1*m12)
            
#endif      

            call Eij_disp(-p1f,-p1f+qf,-p3f,-p3f-p4f,
     >           m22*re,me2*re,me2*re,mp2*re,m12*re,x34,
     >           disp,E0md,Eimd,Eijmd)
            
            if(iscoll.eq.1) then
               call Eij_dispq(-p1f,-p1f+qf,-p3f,-p3f-p4f,
     >              m22*re,me2*re,me2*re,mp2*re,m12*re,x34,
     >              disp,E0mq,Eimq,Eijmq)   
               E0md  = E0mq
               Eimd  = Eimq
               Eijmd = Eijmq
            endif
            E0ae0x24px34x15x13epx25m2eepm1   = E0md
            E1ae0x24px34x15x13epx25m2eepm1   = Eimd(1)
            E1ax150x13px34ex24x25pem2eepm1   = Eimd(2)
            E3ae0x24px34x15x13epx25m2eepm1   = Eimd(3)
            E4ae0x24px34x15x13epx25m2eepm1   = Eimd(4)
            E00ae0x24px34x15x13epx25m2eepm1  = Eijmd(1)
            E11ae0x24px34x15x13epx25m2eepm1  = Eijmd(2)
            E12ae0x24px34x15x13epx25m2eepm1  = Eijmd(3)
            E13ae0x24px34x15x13epx25m2eepm1  = Eijmd(4)
            E14ae0x24px34x15x13epx25m2eepm1  = Eijmd(5)
            E11ax150x13px34ex24x25pem2eepm1  = Eijmd(6)
            E13ax150x13px34ex24x25pem2eepm1  = Eijmd(7)
            E14ax150x13px34ex24x25pem2eepm1  = Eijmd(8)
            E33ae0x24px34x15x13epx25m2eepm1  = Eijmd(9)
            E34ae0x24px34x15x13epx25m2eepm1  = Eijmd(10)     
            E44ae0x24px34x15x13epx25m2eepm1  = Eijmd(11)

            call Di_disp(-p4f,p2f-p4f,p3f,mp2*re,m12*re,
     >           me2*re,m22*re,x34,disp,D0md,Dimd)
            D0apex15px24x34pm1em2 = D0md
            D1apex15px24x34pm1em2 = Dimd(1)
            D2apex15px24x34pm1em2 = Dimd(2)
            D3apex15px24x34pm1em2 = Dimd(3)
            
            
            call Di_disp(-p4f,-p1f+p3f,p3f,mp2*re,m12*re,
     >           me2*re,m22*re,x34,disp,D0md,Dimd)
            
            D0apx25epx13x34pm1em2 = D0md
            D2apx25epx13x34pm1em2 = Dimd(2)
            D1apx25epx13x34pm1em2 = Dimd(1)
            
            
            call Di_disp(p4f,p3f+p4f,p2f-qf,m12*re,mp2*re,
     >           m22*re,me2*re,x34,disp,D0md,Dimd)
            
            D1appex25x34x13m1pm2e = Dimd(1)  
            
            call Di_disp(-p2f,-qf,p1f-qf,me2*re,m12*re,me2*re,
     >           m22*re,x34,disp,D0md,Dimd)
            D1aex25ex150x34em1em2 = Dimd(1)
            D0a0ex34ex15x25eem2m1 = D0md
            D1a0ex34ex15x25eem2m1 = Dimd(2)
            
            call Di_disp(qf,-p2f+qf,p1f,me2*re,me2*re,m12*re,
     >           m22*re,x34,disp,D0md,Dimd)
            D1a0ex34ex25x15eem1m2 =Dimd(1)
            
            
            call Di_disp(qf,-p2f+qf,p1f-p3f,me2*re,me2*re,
     >           m12*re,mp2*re,x34,disp,D0md,Dimd)
            
            D1a0epx13x25x24eem1p = Dimd(1)
            D0a0epx13x25x24eem1p = D0md
            D1aepx130x24x25em1pe = Dimd(2)
            
            call Di_disp(-qf,p4f-p2f,-p2f,me2*re,me2*re,
     >           mp2*re,m12*re,x34,disp,D0md,Dimd)
            D1a0x13pex24x25eepm1 =Dimd(1)

!     Solo se disp=1 mi serve il tutto con i coefficienti s'

            
            include 'pipirad/gvmd/pent3_isr.f'
            if(ind.eq.2) then
               elmat1=-elmat1
            endif
            
            elmatresp(1,1) = elmatresp(1,1)+real(elmat1)
            elmatimsp(1,1) = elmatimsp(1,1)+aimag(elmat1)

            include 'pipirad/gvmd/pent3_fsr.f'
            
            elmatresp(1,2) = elmatresp(1,2)+real(elmat2)
            elmatimsp(1,2) = elmatimsp(1,2)+aimag(elmat2)

         endif

         dispbck=disp
         if((disp.eq.1).and.onlysp.eq.1) then
            continue
         else

            if(disp.eq.1) disp=disp+3
!     print*,'DISP CAMBIA SE DISP=1',disp
!     Qui ora valuto la i coefficienti in \eta e non in s' ( sia in penti_isr/fsr.f
!     che nella riduzione delle Eij)

            if(disp.ne.2) then
#ifdef COLLIER

               call C0_cll(C0ax24peepm1,c1*x24,c1*mp2,c1*me2,
     >              c1*me2,c1*mp2,c1*m12)

               call C0_cll(C0ax13px25epm1,c1*x13,c1*mp2,c1*x25,
     >              c1*me2,c1*mp2,c1*m12)

               call C0_cll(C0appx34m2pm1,c1*mp2,c1*mp2,c1*x34,
     >              c1*m22,c1*mp2,c1*m12)
               
#endif
               if(disp.eq.3) then
                  C0ax24peepm1=setrez(C0ax24peepm1)
                  C0ax13px25epm1=setrez(C0ax13px25epm1)
                  C0appx34m2pm1=setrez(C0appx34m2pm1)
               endif
               
            endif
            
            call Eij_disp(-p1f,-p1f+qf,-p3f,-p3f-p4f,
     >           m22*re,me2*re,me2*re,mp2*re,m12*re,x34,
     >           disp,E0md,Eimd,Eijmd)
!!!!  Quad if needed
            if(iscoll.eq.1) then
               call Eij_dispq(-p1f,-p1f+qf,-p3f,-p3f-p4f,
     >              m22*re,me2*re,me2*re,mp2*re,m12*re,x34,
     >              disp,E0mq,Eimq,Eijmq)
               E0md  = E0mq
               Eimd  = Eimq
               Eijmd = Eijmq

            endif

            E0ae0x24px34x15x13epx25m2eepm1   = E0md
            E1ae0x24px34x15x13epx25m2eepm1   = Eimd(1)
            E1ax150x13px34ex24x25pem2eepm1   = Eimd(2)
            E3ae0x24px34x15x13epx25m2eepm1   = Eimd(3)
            E4ae0x24px34x15x13epx25m2eepm1   = Eimd(4)
            E00ae0x24px34x15x13epx25m2eepm1  = Eijmd(1)
            E11ae0x24px34x15x13epx25m2eepm1  = Eijmd(2)
            E12ae0x24px34x15x13epx25m2eepm1  = Eijmd(3)
            E13ae0x24px34x15x13epx25m2eepm1  = Eijmd(4)
            E14ae0x24px34x15x13epx25m2eepm1  = Eijmd(5)
            E11ax150x13px34ex24x25pem2eepm1  = Eijmd(6)
            E13ax150x13px34ex24x25pem2eepm1  = Eijmd(7)
            E14ax150x13px34ex24x25pem2eepm1  = Eijmd(8)
            E33ae0x24px34x15x13epx25m2eepm1  = Eijmd(9)
            E34ae0x24px34x15x13epx25m2eepm1  = Eijmd(10)     
            E44ae0x24px34x15x13epx25m2eepm1  = Eijmd(11)

            if(disp.ne.4) disp=disp+3
!     print*,'DISP CAMBIA SE DISP NE 1(4)',disp
            call Di_disp(-p4f,p2f-p4f,p3f,mp2*re,m12*re,
     >           me2*re,m22*re,x34,disp,D0md,Dimd)
            D0apex15px24x34pm1em2 = D0md
            D1apex15px24x34pm1em2 = Dimd(1)
            D2apex15px24x34pm1em2 = Dimd(2)
            D3apex15px24x34pm1em2 = Dimd(3)
            
            
            call Di_disp(-p4f,-p1f+p3f,p3f,mp2*re,m12*re,
     >           me2*re,m22*re,x34,disp,D0md,Dimd)
            
            D0apx25epx13x34pm1em2 = D0md
            D2apx25epx13x34pm1em2 = Dimd(2)
            D1apx25epx13x34pm1em2 = Dimd(1)
            
            
            call Di_disp(p4f,p3f+p4f,p2f-qf,m12*re,mp2*re,
     >           m22*re,me2*re,x34,disp,D0md,Dimd)
            
            D1appex25x34x13m1pm2e = Dimd(1)  
            
            call Di_disp(-p2f,-qf,p1f-qf,me2*re,m12*re,me2*re,
     >           m22*re,x34,disp,D0md,Dimd)
            D1aex25ex150x34em1em2 = Dimd(1)
            D0a0ex34ex15x25eem2m1 = D0md
            D1a0ex34ex15x25eem2m1 = Dimd(2)
            
            call Di_disp(qf,-p2f+qf,p1f,me2*re,me2*re,m12*re,
     >           m22*re,x34,disp,D0md,Dimd)
            D1a0ex34ex25x15eem1m2 =Dimd(1)
            
            
            call Di_disp(qf,-p2f+qf,p1f-p3f,me2*re,me2*re,
     >           m12*re,mp2*re,x34,disp,D0md,Dimd)
            
            D1a0epx13x25x24eem1p = Dimd(1)
            D0a0epx13x25x24eem1p = D0md
            D1aepx130x24x25em1pe = Dimd(2)
            
            call Di_disp(-qf,p4f-p2f,-p2f,me2*re,me2*re,
     >           mp2*re,m12*re,x34,disp,D0md,Dimd)
            D1a0x13pex24x25eepm1 =Dimd(1)

            disp = dispbck
            
            do i=1,2
               if (masses(i).ge.4.d0*ampi**2) then
                  if (i.eq.1) m12 = x34*re
                  if (i.eq.2) m22 = x34*re
               endif
            enddo

            include 'pipirad/gvmd/pent3_isr.f'
            if(ind.eq.2) then
               elmat1=-elmat1
            endif
            
            elmatres(1,1) = elmatres(1,1)+real(elmat1)
            elmatims(1,1) = elmatims(1,1)+aimag(elmat1)

            include 'pipirad/gvmd/pent3_fsr.f'
            
            elmatres(1,2) = elmatres(1,2)+real(elmat2)
            elmatims(1,2) = elmatims(1,2)+aimag(elmat2)

            m12  =  masses(1)*(1.d0,0.d0)
            m22  =  masses(2)*(1.d0,0.d0)

            disp = dispbck
         endif
      endif
      elmatres  = elmatres*4.d0*alpha*e6/(2.d0*pi)
      elmatims	= elmatims*4.d0*alpha*e6/(2.d0*pi)
      elmatresp = elmatresp*4.d0*alpha*e6/(2.d0*pi)
      elmatimsp	= elmatimsp*4.d0*alpha*e6/(2.d0*pi)

      end
      
      

      

      !===================================================!
      !  GVMD ROUTINES FOR IFI IN \PI\PI\GAMMA            !   MG 13/02/26
      !===================================================!

!===================================================!                                                                                                                     
!  ROUTINE FOR IFI WITH MASSES M1,M2 IN \PI\PI\GAMMA!                                                                                                               
!===================================================!

      subroutine virtualIFI_ij_GVMD(p1,p2,p3,p4,m12,m22
     >     ,q,virtual,qph)
#ifdef COLLIER 
        use collier
#endif
	implicit double complex (z)
!        implicit none
        double complex m12,m22,m1,m2
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision qph(40,0:3)
	double complex virtual,virtualpent,virtualpent1
        complex*32 virtual1,virtualq
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
	real*8 ep2,ep3,e2,e6,me2,me4,me6,mp2,mp4,mp6
        complex*32 me2q
	double precision ame,ammu,convfac,alpha,pi,ame2,ampi2,ampi
	common/parameters/ame,ammu,convfac,alpha,pi
        common/hparam/ampi
	Complex*16 elmat1,elmat2,elmat3,elmat4,elmat

	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min
	common/softphotonmass/phmass
	real*8 phmass
	
	external dot
	real*8 dot
	
        double complex vpolc,alphar
        
        common/ialpharunning/iarun 
        integer iarun 	

        
        complex*16 D0md,re
	complex*16 Dimd(3)
        complex*16 Dijmd(7)
        complex*32 D0mq,rq
        complex*32 Dimq(3)
        complex*32 Dijmq(7)
        complex*16 C0md,C1md,C2md
**********FF

	integer iffpi
	character*10 what_ffpi
	common/ffpioption/iffpi,what_ffpi
	double complex ffpiis
	double precision ffpiis2
	double complex ffpifs
	double precision ffpifs2

*************GVMD COMMON
        common/gvmdrad/indexnores
        integer indexnores
        
        common/gvmdchecks/iscoll,fscoll
        integer iscoll,fscoll
        common/checkform/intltgvmd
        integer intltgvmd
        
        common/compotontensor/compttens
        integer compttens
*********CLL STUFF
        double complex cnpi(0:0,0:1,0:1)
        double complex cnpiuv(0:0,0:1,0:1)
        double complex dnpi(0:1,0:2,0:2,0:2)
	double complex dnpiuv(0:1,0:2,0:2,0:2)
        double complex enpi(0:1,0:2,0:2,0:2,0:2)
        double complex enpiuv(0:1,0:2,0:2,0:2,0:2)
	double complex c1
	parameter (c1 = (1.d0,0.d0))

**************  QUARTIC-VERTEX  PAVE *******************
!!! DIAGRAMS 1-4 IFI_TRIANG AND DIAGRAM 13 IFI_BOXES !!!

	double complex B0aeem1
	double complex B0aeem2
	double complex B0ax12m1m2
	double complex B0ax15em1
	double complex B0ax25em2
	double complex C0aeex12m1em2
	double complex C0aex15x34m2em1
	double complex B0ax34m1m2
	double complex C0aex25x34m1em2

        double complex C0ax15ex34m1em2
        double complex D0aex15ex250x34em1em2
        double complex D1aex15ex250x34em1em2
        double complex D3aex15ex250x34em1em2
        double complex D00aex15ex250x34em1em2
        double complex D11aex15ex250x34em1em2
        double complex D12aex15ex250x34em1em2
        double complex D13aex15ex250x34em1em2
        double complex D23aex15ex250x34em1em2
        double complex D33aex15ex250x34em1em2
        
**************  BOXES  PAVE *******************
!!!!!!!!!! DIAGRAMS 1-12 IFI_BOXES !!!!!!!!!!!!
        double complex C0ax12px45m2m1p
        double complex D0aepx45ex13x12em1pm2
        double complex D1aeex45px12x13m1em2p
        double complex D2aepx45ex13x12em1pm2
        double complex D3aepx45ex13x12em1pm2

        double complex C0ax12px45m1m2p
        double complex D0aex45pex23x12em1pm2
        double complex D1aeepx45x12x23m1em2p
        double complex D1aex45pex23x12em1pm2
        double complex D2aex45pex23x12em1pm2

        double complex C0ax12px35m2m1p
        double complex D0aex12x35x14epem1m2p
        double complex D1aeex35px12x14m1em2p
        double complex D2aex12x35x14epem1m2p
        double complex D3aex12x35x14epem1m2p

        double complex C0ax12px35m1m2p
        double complex D0aex12px24ex35em1m2p
        double complex D1aeepx35x12x24m1em2p
        double complex D1aex12px24ex35em1m2p
        double complex D3aex12px24ex35em1m2p

        double complex C0aex34x15em2m1
        double complex C0appx34m2pm1
        double complex C1apx34ppm1m2
        double complex C1apx34ppm2m1
        double complex D0ax15x34px23epem1m2p
        double complex D1ax15x34px23epem1m2p
        double complex D3ax15x34px23epem1m2p
        double complex D1aex23px34px15m2epm1

        double complex C0appx34m1pm2
        double complex D0ax15x34px24epem1m2p
        double complex D1aex24px34px15m2epm1
        double complex D1ax15x34px24epem1m2p
        double complex D3ax15x34px24epem1m2p

        double complex C0aex34x25em1m2
        double complex D0aex34px13x25pem1m2p
        double complex D1aex13px34px25m1epm2
        double complex D2aex34px13x25pem1m2p
        double complex D3aex34px13x25pem1m2p

        double complex D0aex34px14x25pem1m2p
        double complex D1aex14px34px25m1epm2
        double complex D2aex34px14x25pem1m2p
        double complex D3aex34px14x25pem1m2p
        
        double complex C1apx35x12m1pm2

        double complex C1apx35x12m2pm1

        double complex C1apx45x12m1pm2
        double complex D0aex12x45x13epem1m2p
        double complex D2aex12x45x13epem1m2p
        double complex D3aex12x45x13epem1m2p

        double complex C1apx45x12m2pm1
        double complex D0aex12px23ex45em1m2p
        double complex D1aex12px23ex45em1m2p
        double complex D3aex12px23ex45em1m2p

        
        B0aeem1 = (0d0,0d0)
	B0aeem2 = (0d0,0d0)
	B0ax12m1m2 = (0d0,0d0)
	B0ax15em1 = (0d0,0d0)
	B0ax25em2 = (0d0,0d0)
	C0aeex12m1em2 = (0d0,0d0)
	C0aex15x34m2em1 = (0d0,0d0)
	B0ax34m1m2 = (0d0,0d0)
	C0aex25x34m1em2 = (0d0,0d0)
        
        re =(1.d0,0.d0)

        e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2
        ep2= 1.d0
        ep3= 1.d0
	ame2=ame*ame
	ampi2=ampi*ampi
        me2=ame2
        mp2=ampi2
        me4=me2**2
        me6=me2**3
        mp4=mp2**2
        mp6=mp2**3

        me2q = (ame*1.q0)*(ame*1.q0)
	p1p2 = dot(p1,p2)
	p1p3 = dot(p1,p3)
	p1p4 = dot(p1,p4)
	kp1  = dot(p1,q)
	p2p3 = dot(p2,p3)
	p2p4 = dot(p2,p4) 
	kp2  = dot(p2,q)
	kp3 = dot(p3,q)
	p3p4 = dot(p3,p4)
	kp4  = dot(p4,q)
					
	x12  = 2d0*ame2+2.d0*p1p2
	x13  = ame2+ampi2-2.d0*p1p3
	x14  = ame2+ampi2-2.d0*p1p4
	x15  = ame2-2.d0*kp1
	x23  = ame2+ampi2-2d0*p2p3
	x24  = ame2+ampi2-2d0*p2p4
	x25 = ame2-2d0*kp2
	x34 = 2d0*ampi2+2d0*p3p4
	x35 = ampi2+2d0*kp3
	x45 = ampi2+2d0*kp4

        m1  = sqrt(m12)
        m2  = sqrt(m22)

**************
	ffpifs=(1.d0,0.d0)
	ffpifs2=1.d0

	ffpiis=(1.d0,0.d0)
	ffpiis2=1.d0
	if(iffpi.gt.0) then
		call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2) 
		call ffpi(2.d0*(ampi**2+p3p4),ffpiis,ffpiis2)
	         
	endif	

******ALPHA RUNNING	
	ffpiis=ffpiis*vpolc(2.d0*(ampi**2+p3p4))
	ffpiis2=ffpiis*CONJG(ffpiis)
	
	ffpifs=ffpifs*vpolc(2.d0*(ame**2+p1p2))
	ffpifs2=ffpifs*CONJG(ffpifs)
****************

#ifdef COLLIER

	call B0_cll(B0aeem1,c1*ame2,c1*ame2,c1*m12)

        call B0_cll(B0aeem2,c1*ame2,c1*ame2,c1*m22)

        call B0_cll(B0ax12m1m2,c1*x12,c1*m12,c1*m22)
        
        call B0_cll(B0ax15em1,c1*x15,c1*ame2,c1*m12)

        call B0_cll(B0ax25em2,c1*x25,c1*ame2,c1*m22)

        call B0_cll(B0ax34m1m2,c1*x34,c1*m12,c1*m22)
        
        call C0_cll(C0aeex12m1em2,c1*ame2,c1*ame2,c1*x12,
     >       c1*m12,c1*ame2,c1*m22)

        call C0_cll(C0aex15x34m2em1,c1*ame2,c1*x15,c1*x34,
     >       c1*m22,c1*ame2,c1*m12)

        call C0_cll(C0aex25x34m1em2,c1*ame2,c1*x25,c1*x34,
     >       c1*m12,c1*ame2,c1*m22)

        call C0_cll(C0ax15ex34m1em2,c1*x15,c1*ame2,c1*x34,
     >       c1*m12,c1*ame2,c1*m22)

        
#endif

        call Dij_d(-p1,-q,+p2-q,me2*re,m12,me2*re,m22,D0md,Dimd,Dijmd)

       D0aex15ex250x34em1em2   = D0md
       D1aex15ex250x34em1em2   = Dimd(1)
       D3aex15ex250x34em1em2   = Dimd(3)
       D00aex15ex250x34em1em2  = Dijmd(1)
       D11aex15ex250x34em1em2  = Dijmd(2)
       D12aex15ex250x34em1em2  = Dijmd(3)
       D13aex15ex250x34em1em2  = Dijmd(4)
       D23aex15ex250x34em1em2  = Dijmd(6)
       D33aex15ex250x34em1em2  = Dijmd(7)
       
!!!!! QUI NON BISOGNA PIÙ MOLTIPLICARE ESPLICITAMENTE PER NESSUN FF IL VIRTUALE !!!!!!!!!
!!!!!            VIENE IN MODO NATURALE DALL'INTEGRAZIONE DI LOOP               !!!!!!!!!  
!!!!!                   IL BORN INVECE VA TENUTO ESPLICITO                      !!!!!!!!! 
       
********QUARTIC VERTEX
       virtual = (0.d0,0.d0)
**************** CROSSATO CON BORN ISR

	include 'pipirad/gvmd/ifitriang_isr.f'

**************************
	virtual=elmat1*CONJG(ffpiis)

****************CROSSATO CON BORN FSR
	
	include 'pipirad/gvmd/ifitriang_fsr.f'
**************************
	virtual=virtual+elmat2*CONJG(ffpifs)

**************BOXES  *******************
!!!!!!!!!! DIAGRAMS 1-12 IFI_BOXES !!!!!!!!!!!!


#ifdef COLLIER

        call C0_cll(C0ax12px45m2m1p,x12*c1,mp2*c1,x45*c1,m22*c1
     >       ,m12*c1,mp2*c1)

        call C0_cll(C0ax12px45m1m2p,x12*c1,mp2*c1,x45*c1,m12*c1
     >	      ,m22*c1,mp2*c1)

        call C0_cll(C0ax12px35m2m1p,x12*c1,mp2*c1,x35*c1,m22*c1
     >       ,m12*c1,mp2*c1)

        call C0_cll(C0ax12px35m1m2p,x12*c1,mp2*c1,x35*c1,m12*c1
     >        ,m22*c1,mp2*c1)

        call C0_cll(C0aex34x15em2m1,c1*me2,c1*x34,c1*x15,me2*c1,
     >       c1*m22,c1*m12)

        call C0_cll(C0appx34m2pm1,mp2*c1,mp2*c1,x34*c1,m22*c1,
     >       mp2*c1,m12*c1)

        call C0_cll(C0appx34m1pm2,mp2*c1,mp2*c1,x34*c1,m12*c1,
     >       mp2*c1,m22*c1)

        call C0_cll(C0aex34x25em1m2,c1*me2,c1*x34,c1*x25,me2*c1,
     >       c1*m12,c1*m22)

#endif

!!!!!!!!Diag 1        
        call Di_d(-p1,-p1+p3,p2,me2*re,m12,mp2*re,m22,D0md,Dimd)

        D0aepx45ex13x12em1pm2  = D0md
        D2aepx45ex13x12em1pm2  = Dimd(2)
        D3aepx45ex13x12em1pm2  = Dimd(3)

        call Di_d(p1,p2+p1,p3,m12,me2*re,m22,mp2*re,D0md,Dimd)

        D1aeex45px12x13m1em2p  = Dimd(1)

!!!!!!!Diag 2
        call Di_d(-p1,p2-p3,p2,me2*re,m12,mp2*re,m22,D0md,Dimd)

        D0aex45pex23x12em1pm2  = D0md
        D1aex45pex23x12em1pm2  = Dimd(1)
        D2aex45pex23x12em1pm2  = Dimd(2)


        call Di_d(p1,p2+p1,p4+q,m12,me2*re,m22,mp2*re,D0md,Dimd)

        D1aeepx45x12x23m1em2p  = Dimd(1)

!!!!!!Diag 3
        call Di_d(-p1,p2,p4-p1,me2*re,m12,m22,mp2*re,D0md,Dimd)

        D0aex12x35x14epem1m2p  = D0md
        D2aex12x35x14epem1m2p  = Dimd(2)
        D3aex12x35x14epem1m2p  = Dimd(3)

        call Di_d(p1,p2+p1,p4,m12,me2*re,m22,mp2*re,D0md,Dimd)

        D1aeex35px12x14m1em2p  = Dimd(1)

!!!!!!Diag 4
        call Di_d(-p1,p2,p2-p4,me2*re,m12,m22,mp2*re,D0md,Dimd)
        
        D0aex12px24ex35em1m2p = D0md
        D1aex12px24ex35em1m2p = Dimd(1)
        D3aex12px24ex35em1m2p = Dimd(3)

        call Di_d(p1,p1+p2,p3+q,m12,me2*re,m22,mp2*re,D0md,Dimd)
        D1aeepx35x12x24m1em2p = Dimd(1)

!!!!!!Diag 5

        call Ci_d(-p4,+p3,mp2*re,m12,m22,C0md,C1apx34ppm1m2,C2md)
        call Ci_d(-p4,+p3,mp2*re,m22,m12,C0md,C1apx34ppm2m1,C2md)

        call Di_d(-p1+q,p2,p2-p3,me2*re,m12,m22,mp2*re,D0md,Dimd)
        D0ax15x34px23epem1m2p = D0md
        D1ax15x34px23epem1m2p = Dimd(1)
        D3ax15x34px23epem1m2p = Dimd(3)

        call Di_d(-p2,-p3,-p3-p4,m22,me2*re,mp2*re,m12,D0md,Dimd)
        D1aex23px34px15m2epm1 = Dimd(1)

!!!!!!Diag 6
        
        call Di_d(-p1+q,p2,p2-p4,me2*re,m12,m22,mp2*re,D0md,Dimd)
        D0ax15x34px24epem1m2p = D0md
        D1ax15x34px24epem1m2p = Dimd(1)
        D3ax15x34px24epem1m2p = Dimd(3)

        call Di_d(-p2,-p4,-p3-p4,m22,me2*re,mp2*re,m12,D0md,Dimd)
        D1aex24px34px15m2epm1 = Dimd(1)

!!!!!!Diag 7
        
        call Di_d(p1,-p2+q,p1-p3,me2*re,m12,m22,mp2*re,D0md,Dimd)
        D0aex34px13x25pem1m2p = D0md
        D2aex34px13x25pem1m2p = Dimd(2)
        D3aex34px13x25pem1m2p = Dimd(3)

        call Di_d(-p1,-p3,-p3-p4,m12,me2*re,mp2*re,m22,D0md,Dimd)
        D1aex13px34px25m1epm2 = Dimd(1)

!!!!!!Diag 8

        call Di_d(p1,-p2+q,p1-p4,me2*re,m12,m22,mp2*re,D0md,Dimd)
        D0aex34px14x25pem1m2p = D0md
        D2aex34px14x25pem1m2p = Dimd(2)
        D3aex34px14x25pem1m2p = Dimd(3)

        call Di_d(-p1,-p4,-p3-p4,m12,me2*re,mp2*re,m22,D0md,Dimd)
        D1aex14px34px25m1epm2 = Dimd(1)

!!!!!!Diag 9

        call Ci_d(p4,p1+p2,m12,mp2*re,m22,C0md,C1apx35x12m1pm2,C2md)

!!!!!!Diag 10
        call Ci_d(p4,p1+p2,m22,mp2*re,m12,C0md,C1apx35x12m2pm1,C2md)

!!!!!!Diag 11
        call Ci_d(p3,p1+p2,m12,mp2*re,m22,C0md,C1apx45x12m1pm2,C2md)

        D0aex12x45x13epem1m2p = D0aepx45ex13x12em1pm2
        D2aex12x45x13epem1m2p = D3aepx45ex13x12em1pm2
        D3aex12x45x13epem1m2p = D2aepx45ex13x12em1pm2

!!!!!!Diag 12
        call Ci_d(p3,p1+p2,m22,mp2*re,m12,C0md,C1apx45x12m2pm1,C2md)

        D0aex12px23ex45em1m2p = D0aex45pex23x12em1pm2
        D1aex12px23ex45em1m2p = D1aex45pex23x12em1pm2
        D3aex12px23ex45em1m2p = D2aex45pex23x12em1pm2

**************** CROSSATO CON BORN ISR

	include 'pipirad/gvmd/ifiboxes_isr.f'
**************************
	virtual=virtual+elmat1*CONJG(ffpiis)


        
****************CROSSATO CON BORN FSR
	
	include 'pipirad/gvmd/ifiboxes_fsr.f'
**************************
        virtual = virtual + elmat2*CONJG(ffpifs)

        if(compttens.ne.0) then
           if(compttens.eq.1) then
              ep3=0.d0
           elseif(compttens.eq.2) then
              ep2=0.d0
           endif
           virtual=0.d0
           include 'pipirad/gvmd/triang_isr.f'

           virtual=elmat1*CONJG(ffpiis)
           include 'pipirad/gvmd/triang_fsr.f'

           virtual=virtual+elmat2*CONJG(ffpifs)
!           print*,'UP TO TRIANG',virtual*4.d0*e6*alpha/(2d0*pi)
           include	'pipirad/gvmd/box_isr.f'
           virtual=virtual+elmat1*CONJG(ffpiis)

           include	'pipirad/gvmd/box_fsr.f'
           virtual=virtual+elmat2*CONJG(ffpifs)
           
!           print*,'UP TO BOXES',virtual*4.d0*e6*alpha/(2d0*pi)
        endif

        
!        virtual = virtual*4.d0*e6*alpha/(2.d0*pi)
!        return
**************PENT  *******************
!!!!!!!!!!DIAGRAMS 1-4 IFI_PENT!!!!!!!!!
                
*--------------------------------------------------------------------------------------------*        
!!!!!!!!!!HERE WE CALL VirtualIFI_Pent TO GET THE AMPLITUDE WITH M1=M2=0!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!THEN WE CALL TWICE VirtualIFI_ij_Pent WHERE WE EVALUATE DIAGRAMS!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!1 AND 3 OF IFI_PENT WITH THE FIRST CALL AND 2-4 WITH THE SECOND!!!!!!!!!!!!!!!!!!!!!        
*--------------------------------------------------------------------------------------------*
        
c$$$        if((abs(m12).le.1d-17).and.(abs(m22).le.1d-17)) then
c$$$           indexnores=1
c$$$!!!!  CON LA COMMON QUI SOPRA CHIAMO LE ROUTINE CHE HO GIÀ PER I PIONI FATTORIZZATI (M1=M2=0) MA ESCLUDENDO DAL CALCOLO LA MOLTIPLICAZIONE ESPLICITA DEL FF DEL VIRTUALE, CHE NON CI DEVE ESSERE NEL GVMD E RITORNA UN VALORE COMPLEX
c$$$           call VirtualIFI_Pent(p1,p2,p3,p4,q,virtualpent,qph)
c$$$           virtual=virtual+virtualpent/(4.d0*e6*alpha/(2.d0*pi))
c$$$           indexnores=0
c$$$        else
!           compttens=0
           call VirtualIFI_ij_Pent(p1,p2,p3,p4,q,m12,m22
     >          ,virtualpent,qph,1)
!           stop
           virtual=virtual+virtualpent

           call VirtualIFI_ij_Pent(p1,p2,p3,p4,q,m12,m22
     >          ,virtualpent,qph,2)
           virtual=virtual+virtualpent
           
c$$$        endif
        virtual=4d0*virtual*e6*alpha/(2d0*pi)

	end

      subroutine VirtualIFI_ij_Pent(p1,p2,p3,p4,q,m12,m22
     >     ,virtual,qph,ind) 
#ifdef COLLIER
      use collier
#endif
      implicit double complex (z)
!      implicit none
      double complex m12,m22,m1,m2
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision qph(40,0:3),p3t(0:3),p4t(0:3)
      double complex virtual,virtualpent
      integer ind
      complex*32 virtual1
      real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
      real*8 x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
      real*8 e2,e6,me2,me4,me6,mp2,mp4,mp6
      complex*32 me2q
      double precision ame,ammu,convfac,alpha,pi,ame2,ampi2,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      Complex*16 elmat1,elmat2,elmat3,elmat4,elmat

      real*8 deltauv,lambda2,mu2dim,m2min
      common/loop/deltauv,lambda2,mu2dim,m2min
      common/softphotonmass/phmass
      real*8 phmass
      
      external dot
      real*8 dot
      
      double complex vpolc,alphar
      
      common/ialpharunning/iarun 
      integer iarun 	

      
      complex*16 D0md,re
      complex*32 rq
      complex*16 Dimd(3)
      complex*16 Dijmd(7)
      complex*16 C0md,C1md,C2md
      complex*16 E0md
      complex*16 Eimd(4)
      complex*16 Eijmd(11)
      complex*32 E0mq
      complex*32 Eimq(4)
      complex*32 Eijmq(11)

      complex*32 D0mq
      complex*32 Dimq(3)
**********FF

      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi
      double complex ffpiis
      double precision ffpiis2
      double complex ffpifs
      double precision ffpifs2

*************GVMD COMMON
      common/gvmdrad/indexnores
      integer indexnores
      
      common/gvmdchecks/iscoll,fscoll
      integer iscoll,fscoll

      
      common/compotontensor/compttens
      integer compttens

      double precision ep2,ep3
*********CLL STUFF
      double complex cnpi(0:0,0:1,0:1)
      double complex cnpiuv(0:0,0:1,0:1)
      double complex dnpi(0:1,0:2,0:2,0:2)
      double complex dnpiuv(0:1,0:2,0:2,0:2)
      double complex enpi(0:1,0:2,0:2,0:2,0:2)
      double complex enpiuv(0:1,0:2,0:2,0:2,0:2)
      double complex c1
      parameter (c1 = (1.d0,0.d0))

**************PENT  PAVE *******************
!!!!!!!!!!DIAGRAMS 1 IFI_PENT !!!!!!!!!!!!

      double complex D0ap0px12x35x45m1ppm2
      double complex D1ap0px12x35x45m1ppm2
      double complex D3ap0px12x35x45m1ppm2
      double complex D1ax350x45x12ppm1ppm2
      double complex D00ap0px12x35x45m1ppm2
      double complex D11ap0px12x35x45m1ppm2
      double complex D12ap0px12x35x45m1ppm2
      double complex D13ap0px12x35x45m1ppm2
      double complex D11ax350x45x12ppm1ppm2
      double complex D13ax350x45x12ppm1ppm2

      double complex E3aeep0px12x24x45x35x13m1em2pp
      double complex E0aeex450x35x12x13ppx24m1em2pp
      double complex E1aeex450x35x12x13ppx24m1em2pp
      double complex E2aeex450x35x12x13ppx24m1em2pp
      double complex E3aeex450x35x12x13ppx24m1em2pp
      double complex E13aeep0px12x24x45x35x13m1em2pp
      double complex E23aeep0px12x24x45x35x13m1em2pp
      double complex E33aeep0px12x24x45x35x13m1em2pp
      double complex E00aeex450x35x12x13ppx24m1em2pp
      double complex E11aeex450x35x12x13ppx24m1em2pp
      double complex E12aeex450x35x12x13ppx24m1em2pp
      double complex E22aeex450x35x12x13ppx24m1em2pp
      double complex E13aeex450x35x12x13ppx24m1em2pp
      double complex E23aeex450x35x12x13ppx24m1em2pp
      double complex E33aeex450x35x12x13ppx24m1em2pp
      double complex E34aeex450x35x12x13ppx24m1em2pp

**************PENT  PAVE *******************                                                                                                                                                                                                           
!!!!!!!!!!DIAGRAMS 3 IFI_PENT !!!!!!!!!!!!                                                                                                                                                                                                             

      double complex E4ae0x24px34x15x13epx25m2eepm1
      double complex E3ae0x24px34x15x13epx25m2eepm1
      double complex E1ax150x13px34ex24x25pem2eepm1
      double complex E1ae0x24px34x15x13epx25m2eepm1
      double complex E14ax150x13px34ex24x25pem2eepm1 
      double complex E14ae0x24px34x15x13epx25m2eepm1
      double complex E13ax150x13px34ex24x25pem2eepm1 
      double complex E13ae0x24px34x15x13epx25m2eepm1
      double complex E12ae0x24px34x15x13epx25m2eepm1 
      double complex E11ax150x13px34ex24x25pem2eepm1
      double complex E11ae0x24px34x15x13epx25m2eepm1
      double complex E0ae0x24px34x15x13epx25m2eepm1
      double complex E00ae0x24px34x15x13epx25m2eepm1
      double complex E44ae0x24px34x15x13epx25m2eepm1
      double complex E34ae0x24px34x15x13epx25m2eepm1
      double complex E33ae0x24px34x15x13epx25m2eepm1 

      double complex D0apex15px24x34pm1em2
      double complex D3apex15px24x34pm1em2
      double complex D2apex15px24x34pm1em2
      double complex D1apex15px24x34pm1em2
      
      double complex D0apx25epx13x34pm1em2
      double complex D2apx25epx13x34pm1em2
      double complex D1apx25epx13x34pm1em2      
      double complex D1appex25x34x13m1pm2e
      
      double complex D1aex25ex150x34em1em2
      double complex D0a0ex34ex15x25eem2m1
      double complex D1a0ex34ex25x15eem1m2
      double complex D1a0ex34ex15x25eem2m1
      
      double complex D1aepx130x24x25em1pe
      double complex D1a0x13pex24x25eepm1      
      double complex D1a0epx13x25x24eem1p
      double complex D0a0epx13x25x24eem1p

      double complex C0ax24peepm1
      double complex C0ax13px25epm1
      double complex C0appx34m2pm1

      re =(1.d0,0.d0)

      e2 = 4.d0*pi*alpha
      e6 = e2*e2*e2
      ame2=ame*ame
      ampi2=ampi*ampi
      me2=ame2
      mp2=ampi2
      me4=me2**2
      me6=me2**3
      mp4=mp2**2
      mp6=mp2**3

      if(ind.eq.1) then
         p3t=p3
         p4t=p4
      elseif(ind.eq.2) then
         p3t=p4
         p4t=p3
      endif
      me2q = (ame*1.q0)*(ame*1.q0)
      p1p2 = dot(p1,p2)
      p1p3 = dot(p1,p3t)
      p1p4 = dot(p1,p4t)
      kp1  = dot(p1,q)
      p2p3 = dot(p2,p3t)
      p2p4 = dot(p2,p4t) 
      kp2  = dot(p2,q)
      kp3 = dot(p3t,q)
      p3p4 = dot(p3t,p4t)
      kp4  = dot(p4t,q)
      
      x12  = 2d0*ame2+2.d0*p1p2
      x13  = ame2+ampi2-2.d0*p1p3
      x14  = ame2+ampi2-2.d0*p1p4
      x15  = ame2-2.d0*kp1
      x23  = ame2+ampi2-2d0*p2p3
      x24  = ame2+ampi2-2d0*p2p4
      x25 = ame2-2d0*kp2
      x34 = 2d0*ampi2+2d0*p3p4
      x35 = ampi2+2d0*kp3
      x45 = ampi2+2d0*kp4

      m1  = sqrt(m12)
      m2  = sqrt(m22)

**************
      ffpifs=(1.d0,0.d0)
      ffpifs2=1.d0

      ffpiis=(1.d0,0.d0)
      ffpiis2=1.d0
      if(iffpi.gt.0) then
         call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2) 
         call ffpi(2.d0*(ampi**2+p3p4),ffpiis,ffpiis2)
         
      endif	

******ALPHA RUNNING	
      ffpiis=ffpiis*vpolc(2.d0*(ampi**2+p3p4))
      ffpiis2=ffpiis*CONJG(ffpiis)
      
      ffpifs=ffpifs*vpolc(2.d0*(ame**2+p1p2))
      ffpifs2=ffpifs*CONJG(ffpifs)
****************
      ep2=1.d0
      ep3=1.d0
      if(compttens.eq.1) then
         ep3=0.d0
      elseif(compttens.eq.2) then
         ep2=0.d0
      endif

!!!!!!!!!!!!!!!PENT 1(FSR)!!!!!!!!!!!!!!!!!

      if(fscoll.eq.0) then
       call Eij_d(p1,p1+p2,p3t,p3t+q,m12*re,me2*re,m22*re,mp2*re,mp2*re
     >        ,E0md,Eimd,Eijmd)

         E0aeex450x35x12x13ppx24m1em2pp   = E0md
         E1aeex450x35x12x13ppx24m1em2pp   = Eimd(1)
         E2aeex450x35x12x13ppx24m1em2pp   = Eimd(2)
         E3aeex450x35x12x13ppx24m1em2pp   = Eimd(3)
         E3aeep0px12x24x45x35x13m1em2pp   = Eimd(4)
         E00aeex450x35x12x13ppx24m1em2pp  = Eijmd(1)
         E11aeex450x35x12x13ppx24m1em2pp  = Eijmd(2)
         E12aeex450x35x12x13ppx24m1em2pp  = Eijmd(3)
         E13aeex450x35x12x13ppx24m1em2pp  = Eijmd(4)
         E13aeep0px12x24x45x35x13m1em2pp  = Eijmd(5)
         E22aeex450x35x12x13ppx24m1em2pp  = Eijmd(6)
         E23aeex450x35x12x13ppx24m1em2pp  = Eijmd(7)
         E23aeep0px12x24x45x35x13m1em2pp  = Eijmd(8)
         E33aeex450x35x12x13ppx24m1em2pp  = Eijmd(9)
         E34aeex450x35x12x13ppx24m1em2pp  = Eijmd(10)
         E33aeep0px12x24x45x35x13m1em2pp  = Eijmd(11)

!!!SE SONO COLLINEARE A FS, LE E_IJ LE VALUTO CON ROUTINE IN 4-PLA
      
      elseif(fscoll.eq.1) then
        call Eij(p1,p1+p2,p3t,p3t+q,m12*re,me2*re,m22*re,mp2*re,mp2*re
     >      ,E0mq,Eimq,Eijmq)

         E0aeex450x35x12x13ppx24m1em2pp   = E0mq
         E1aeex450x35x12x13ppx24m1em2pp   = Eimq(1)
         E2aeex450x35x12x13ppx24m1em2pp   = Eimq(2)
         E3aeex450x35x12x13ppx24m1em2pp   = Eimq(3)
         E3aeep0px12x24x45x35x13m1em2pp   = Eimq(4)
         E00aeex450x35x12x13ppx24m1em2pp  = Eijmq(1)
         E11aeex450x35x12x13ppx24m1em2pp  = Eijmq(2)
         E12aeex450x35x12x13ppx24m1em2pp  = Eijmq(3)
         E13aeex450x35x12x13ppx24m1em2pp  = Eijmq(4)
         E13aeep0px12x24x45x35x13m1em2pp  = Eijmq(5)
         E22aeex450x35x12x13ppx24m1em2pp  = Eijmq(6)
         E23aeex450x35x12x13ppx24m1em2pp  = Eijmq(7)
         E23aeep0px12x24x45x35x13m1em2pp  = Eijmq(8)
         E33aeex450x35x12x13ppx24m1em2pp  = Eijmq(9)
         E34aeex450x35x12x13ppx24m1em2pp  = Eijmq(10)
         E33aeep0px12x24x45x35x13m1em2pp  = Eijmq(11)
      endif
      
c$$$      call Dij_d(p3t,p3t+q,p1+p2,m12*re,mp2*re,mp2*re,m22*re,
c$$$     >     D0md,Dimd,Dijmd)
c$$$      D0ap0px12x35x45m1ppm2  = D0md
c$$$      D1ap0px12x35x45m1ppm2  = Dimd(1)
c$$$      D1ax350x45x12ppm1ppm2  = Dimd(2)
c$$$      D3ap0px12x35x45m1ppm2  = Dimd(3)
c$$$      D00ap0px12x35x45m1ppm2 = Dijmd(1)
c$$$      D11ap0px12x35x45m1ppm2 = Dijmd(2)
c$$$      D12ap0px12x35x45m1ppm2 = Dijmd(3)
c$$$      D13ap0px12x35x45m1ppm2 = Dijmd(4)
c$$$      D11ax350x45x12ppm1ppm2 = Dijmd(5)
c$$$      D13ax350x45x12ppm1ppm2 = Dijmd(6)

      
      include 'pipirad/gvmd/pent1_isr.f'

      virtual = ep3*elmat1*CONJG(ffpiis)

      include 'pipirad/gvmd/pent1_fsr.f'

      if(ind.eq.2) then
         elmat2=-elmat2
      endif
      
      virtual = virtual+ ep3*elmat2*CONJG(ffpifs)

!!!!!!!!!!!!!!!!!PENT 3 (ISR)!!!!!!!!!!!!!!!!!

#ifdef COLLIER

      call C0_cll(C0ax24peepm1,c1*x24,c1*mp2,c1*me2,
     >     c1*me2,c1*mp2,c1*m12)

      call C0_cll(C0ax13px25epm1,c1*x13,c1*mp2,c1*x25,
     >     c1*me2,c1*mp2,c1*m12)

      call C0_cll(C0appx34m2pm1,c1*mp2,c1*mp2,c1*x34,
     >     c1*m22,c1*mp2,c1*m12)
      
#endif      

      if(iscoll.eq.0) then
         call Eij_d(-p1,-p1+q,-p3t,-p3t-p4t,
     >        m22*re,me2*re,me2*re,mp2*re,m12*re,E0md,Eimd,Eijmd)
         
         E0ae0x24px34x15x13epx25m2eepm1   = E0md
         E1ae0x24px34x15x13epx25m2eepm1   = Eimd(1)
         E1ax150x13px34ex24x25pem2eepm1   = Eimd(2)
         E3ae0x24px34x15x13epx25m2eepm1   = Eimd(3)
         E4ae0x24px34x15x13epx25m2eepm1   = Eimd(4)
         E00ae0x24px34x15x13epx25m2eepm1  = Eijmd(1)
         E11ae0x24px34x15x13epx25m2eepm1  = Eijmd(2)
         E12ae0x24px34x15x13epx25m2eepm1  = Eijmd(3)
         E13ae0x24px34x15x13epx25m2eepm1  = Eijmd(4)
         E14ae0x24px34x15x13epx25m2eepm1  = Eijmd(5)
         E11ax150x13px34ex24x25pem2eepm1  = Eijmd(6)
         E13ax150x13px34ex24x25pem2eepm1  = Eijmd(7)
         E14ax150x13px34ex24x25pem2eepm1  = Eijmd(8)
         E33ae0x24px34x15x13epx25m2eepm1  = Eijmd(9)
         E34ae0x24px34x15x13epx25m2eepm1  = Eijmd(10)     
         E44ae0x24px34x15x13epx25m2eepm1  = Eijmd(11)

      call Di_d(-p4t,p2-p4t,p3t,mp2*re,m12*re,me2*re,m22*re,D0md,Dimd)
         D0apex15px24x34pm1em2 = D0md
         D1apex15px24x34pm1em2 = Dimd(1)
         D2apex15px24x34pm1em2 = Dimd(2)
         D3apex15px24x34pm1em2 = Dimd(3)


       call Di_d(-p4t,-p1+p3t,p3t,mp2*re,m12*re,me2*re,m22*re,D0md,Dimd)

         D0apx25epx13x34pm1em2 = D0md
         D2apx25epx13x34pm1em2 = Dimd(2)
         D1apx25epx13x34pm1em2 = Dimd(1)


       call Di_d(p4t,p3t+p4t,p2-q,m12*re,mp2*re,m22*re,me2*re,D0md,Dimd)
         
         D1appex25x34x13m1pm2e = Dimd(1)  

         call Di_d(-p2,-q,p1-q,me2*re,m12*re,me2*re,m22*re,D0md,Dimd)
         D1aex25ex150x34em1em2 = Dimd(1)
         D0a0ex34ex15x25eem2m1 = D0md
         D1a0ex34ex15x25eem2m1 = Dimd(2)

         call Di_d(q,-p2+q,p1,me2*re,me2*re,m12*re,m22*re,D0md,Dimd)
         D1a0ex34ex25x15eem1m2 =Dimd(1)


         call Di_d(q,-p2+q,p1-p3t,me2*re,me2*re,m12*re,mp2*re,D0md,Dimd)

         D1a0epx13x25x24eem1p = Dimd(1)
         D0a0epx13x25x24eem1p = D0md
         D1aepx130x24x25em1pe = Dimd(2)

         call Di_d(-q,p4t-p2,-p2,me2*re,me2*re,mp2*re,m12*re,D0md,Dimd)
         D1a0x13pex24x25eepm1 =Dimd(1)

      elseif(iscoll.eq.1) then

         call Eij(-p1,-p1+q,-p3t,-p3t-p4t,
     >     m22*re,me2*re,me2*re,mp2*re,m12*re,E0mq,Eimq,Eijmq)

         E0ae0x24px34x15x13epx25m2eepm1   = E0mq
         E1ae0x24px34x15x13epx25m2eepm1   = Eimq(1)
         E1ax150x13px34ex24x25pem2eepm1   = Eimq(2)
         E3ae0x24px34x15x13epx25m2eepm1   = Eimq(3)
         E4ae0x24px34x15x13epx25m2eepm1   = Eimq(4)
         E00ae0x24px34x15x13epx25m2eepm1  = Eijmq(1)
         E11ae0x24px34x15x13epx25m2eepm1  = Eijmq(2)
         E12ae0x24px34x15x13epx25m2eepm1  = Eijmq(3)
         E13ae0x24px34x15x13epx25m2eepm1  = Eijmq(4)
         E14ae0x24px34x15x13epx25m2eepm1  = Eijmq(5)
         E11ax150x13px34ex24x25pem2eepm1  = Eijmq(6)
         E13ax150x13px34ex24x25pem2eepm1  = Eijmq(7)
         E14ax150x13px34ex24x25pem2eepm1  = Eijmq(8)
         E33ae0x24px34x15x13epx25m2eepm1  = Eijmq(9)
         E34ae0x24px34x15x13epx25m2eepm1  = Eijmq(10)     
         E44ae0x24px34x15x13epx25m2eepm1  = Eijmq(11)

         rq= (1.q0,0.q0)
         
         
         call Di(-p4t,p2-p4t,p3t,mp2*rq,m12*rq,me2*rq,m22*rq,D0mq,Dimq)
         D0apex15px24x34pm1em2 = D0mq
         D1apex15px24x34pm1em2 = Dimq(1)
         D2apex15px24x34pm1em2 = Dimq(2)
         D3apex15px24x34pm1em2 = Dimq(3)

         call Di(-p4t,-p1+p3t,p3t,mp2*rq,m12*rq,me2*rq,m22*rq,D0mq,Dimq)

         D0apx25epx13x34pm1em2 = D0mq
         D2apx25epx13x34pm1em2 = Dimq(2)
         D1apx25epx13x34pm1em2 = Dimq(1)

         call Di(p4t,p3t+p4t,p2-q,m12*rq,mp2*rq,m22*rq,me2*rq,D0mq,Dimq)
         
         D1appex25x34x13m1pm2e = Dimq(1)  

         call Di(-p2,-q,p1-q,me2*rq,m12*rq,me2*rq,m22*rq,D0mq,Dimq)
         D1aex25ex150x34em1em2 = Dimq(1)
         D0a0ex34ex15x25eem2m1 = D0mq
         D1a0ex34ex15x25eem2m1 = Dimq(2)

         call Di(q,-p2+q,p1,me2*rq,me2*rq,m12*rq,m22*rq,D0mq,Dimq)
         D1a0ex34ex25x15eem1m2 =Dimq(1)

         call Di(q,-p2+q,p1-p3t,me2*rq,me2*rq,m12*rq,mp2*rq,D0mq,Dimq)

         D1a0epx13x25x24eem1p = Dimq(1)
         D0a0epx13x25x24eem1p = D0mq
         D1aepx130x24x25em1pe = Dimq(2)

         call Di(-q,p4t-p2,-p2,me2*rq,me2*rq,mp2*rq,m12*rq,D0mq,Dimq)
         D1a0x13pex24x25eepm1 =Dimq(1)
      endif

      
      include 'pipirad/gvmd/pent3_isr.f'
      
      if(ind.eq.2) then
         elmat1=-elmat1
      endif
      virtual =  virtual+ ep2*elmat1*CONJG(ffpiis)
      
      include 'pipirad/gvmd/pent3_fsr.f'

      virtual = virtual+ ep2*elmat2*CONJG(ffpifs)

      end

      
!===================================================!                                                                                                                     
!       ROUTINE FOR FULL GVMD IFI IN \PI\PI\GAMMA   !                                                                                                               
!===================================================!

      subroutine virtualIFI_GVMD(p1,p2,p3,p4,q,emtx,qph)
      implicit double precision (b-i,v)
      character*2 fs
      common/finalstate/fs
     
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      double precision qph(40,0:3)
      double complex virtual,virtual1,bornc
      double precision emtx
      double precision elmat
      complex*32 virtualq
      double precision eta
      external tridot,dot
      double precision tridot,dot,eikfull
**********GVMD
      integer i,j,k
      double complex ct,im,cij
      double precision m_n(1:9),g_n(1:9)
      real*16 delta,delta1
      double complex m12,m22
      double complex gvmdij
      double complex a_n(1:9),lv2(1:9),lv2diff
      common/par_gvmd/m_n,g_n,a_n
      double complex math,form
      common/resonancesgvmd/nres
      integer nres
**********FF

      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi
      double complex ffpiis
      double precision ffpiis2
      double complex ffpifs
      double precision ffpifs2
      

      
      common/compotontensor/compttens
      integer compttens
      
      
      virtual   =  (0.d0,0.d0)
      virtualq  =  (0.d0,0.d0)
      if(iffpi.eq.2) then
         m12 = (0.d0,0.d0)
         m22 = (0.d0,0.d0)
         call virtualIFI_ij_GVMD(p1,p2,p3,p4,m12,m22,q,
     >        virtual,qph)
!!!!! GIÀ DEFINITI COSÌ NEL FF,OSSIA A_N/CT
         im=(0.d0,1.d0)           

!         virtual=(0.d0,0.d0)
         do i=1,nres
            lv2(i) = m_n(i)**2 - im*m_n(i)*g_n(i)

            call  virtualIFI_ij_GVMD(p1,p2,p3,p4,m12,lv2(i),q,
     >           gvmdij,qph)
            virtual = virtual -a_n(i)*gvmdij

            call  virtualIFI_ij_GVMD(p1,p2,p3,p4,lv2(i),m22,q,
     >           gvmdij,qph)
            virtual = virtual -a_n(i)*gvmdij

!!!!!!!!!COMPT TENS
c$$$            if(compttens.ne.0) then
c$$$               m12 = (0.d0,0.d0)
c$$$               m22 = (0.d0,0.d0)
c$$$
c$$$               call virtualIFI_ij_GVMD(p1,p2,p3,p4,m12,m22,q,
c$$$     >              gvmdij,qph)
c$$$
c$$$               if(compttens.eq.1) then
c$$$                  eta=dot(p3+p4,p3+p4)
c$$$                  gvmdij=gvmdij*(-1.d0*eta/(lv2(i)-eta))
c$$$               elseif(compttens.eq.2) then
c$$$                  eta=dot(p1+p2,p1+p2)
c$$$                  gvmdij=gvmdij*(-1.d0*eta/(lv2(i)-eta))
c$$$               endif
c$$$               virtual=virtual+a_n(i)*gvmdij
c$$$            endif

         enddo
         
         do i=1,nres
            do j=1,nres
               
               m12 = lv2(i)
               m22 = lv2(j)
               call virtualIFI_ij_GVMD(p1,p2,p3,p4,m12,m22,q,
     >              gvmdij,qph)
               
               virtual = virtual + a_n(i)*a_n(j)*gvmdij
            enddo
         enddo
      else
         print*,'THIS SHOULD NOT HAPPEN, STOPPING THE PROGRAM'
         stop
      endif

      emtx= real(virtual)

      return
      end
	
*-------------------------------------------------------------------------------


      !===================================================!
      !  DISPERSIVE APPROACH - based on arXiv:2207.03495  !   FPU 11.5.24
      !===================================================!

      function virtplussoftpions(p1,p2,p3,p4,eps)

#ifdef COLLIER
      use collier
#endif

      IMPLICIT DOUBLE PRECISION (A-H,L-M,O-Z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),ptmp(0:3)
      character*2 fs
      common/finalstate/fs

      real*8 photonCT,eVertexCT,piVertexCT
      Complex*16 A0e,A0p
      Complex*16 B0ee,B0pp,B0p0p,B0e0e,B00ee,B0see,B0spp,B0s00
      Complex*16 DB0ee,DB0pp,DB00ee,DB0e0e,DB0p0p
      Complex*16 C0eesee,C0ppspp,C0esee00,C0epuep,C0pspp00
      Complex*16 C0eup0ep,C0epte0p,C0ees0e0,C0pps0p0,C0etp0ep
      Complex*16 D0eeppst0e0p,D0ts,D0eeppsu0e0p,D0su
      real*8 deltauv,lambda2,mu2dim,m2min
      common/loop/deltauv,lambda2,mu2dim,m2min
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      common/colliermodedisp/initialisation
      integer initialisation

      complex*16 virtual
      complex*16 virtel,virtph,virtpi
      complex*16 virtboxcross,virtboxsquare,virtboxt
      complex*16 ffi,virtual_box

      complex*16 boxdm
      external boxdmm
      integer ISR,FSR
      common/radiation/ISR,FSR

      real*8 me2,mp2
      real*8 da,dpsi,dphi,de,dmpi,dme
      real*8 cts

      real*8 rsofttot,rsoftISR,rsoftFSR,rsoftIFI,egmin

      real*8 kallen
      external kallen
******FORM FACTOR STUFF

      integer iffpi
      character*10 what_ffpi
      real*8 s_max
      common/ffpioption/iffpi,what_ffpi
      real*8 imffs
      common/smax/s_max

*******Variables needed for the integrals
      real*8 imin,imax,spmin,spmax
      real*8 delta,smdel,spdel
      real*8 Repoledisp
      external Repoledisp

      real*8 intpv

      real*8 redispdisp
      external redispdisp

      real*8 imdispdisp
      external imdispdisp

      real*8 ImFsubtrF_c
      external ImFsubtrF_c

      real*8 impdt,impdu

      real*8 rebox,imbox

      real*8 reboxpp,imboxpp,boxpp
      real*8 reboxpd,imboxpd,boxpd
      real*8 reboxdd,imboxdd,boxdd

      real*8 reir,refi,repp,redd
      real*8 impp,impd,imir,imfi,imdd
      real*8 imbc,imff

      real*8 imfi_t, imfi_u

      double complex c1
      parameter (c1 = (1.d0,0.d0))

******MC INTEGRAL
      real*8 psi(2)
      real*8 xi,csi

      real*8 ff2common

      real*8 boxdm_loopt_integrand
      external boxdm_loopt_integrand
      common/ff2/ff2common

      real*8 MC_1D,mc_2d,mc_pv
      external MC_1D,MC_2d,MC_PV

      real*8 immax

      real*8 box_only_BC_mc
      external box_only_BC_mc

      real*8 im_deltapd_bern
      external im_deltapd_bern
*
      real*8 boxdm_moretti_lambda
      external boxdm_moretti_lambda

      real*8 oneoeps,mudimreg2,eulergamma
      common/IRdimreg/oneoeps,mudimreg2,eulergamma
      complex*16 ilogmoretti,c
*
******DISPERSIVE FSR!!**************
      complex*16 FSRdisp,CTpidisp,FSR_disp,CTpi_disp
      external FSRdisp,CTpidisp
      complex*16 vertexpi,ctdisp
      external vertexpi,ctdisp
      common/indexfsrdisp/iFSRdisp
      integer iFSRdisp
******!! MORETTI LAMBDA CHECKS ********
      real*8 I1moretti,I2moretti,I1ovsptot_lambda
      real*8 Ispovsptot_lambda
      real*8 Itotspovsp_lambda,Itot1ovsp_lambda
      real*8 Itot1ovsp_dimreg,Itotspovsp_dimreg

******RESIDUE added 22.01 ***********
      real*8 residue,delta_disp,resd
      external residue
#ifdef COLLIER
      call GetDeltaUV_cll(deltauv)
      call GetMuUV2_cll(mu2dim)
      call GetMuIR2_cll(lambda2)

#endif
      if(initialisation.eq.0) then
#ifdef COLLIER
         call SetMode_cll(1)
#endif
         initialisation=1
      endif


#ifdef LOOPTOOLS
#include "looptools.h"
      call clearcache
      call ltini
      lambda2=getlambda()
      deltauv=getdelta()
      mu2dim=getmudim()
#endif
*
      mudimreg2= 1.d3
      eulergamma= 0.57721566490153286060d0
      oneoeps= log(lambda2/mudimreg2) + eulergamma - log(4.d0*pi)
*
      do k = 0,3
         ptmp(k) = p1(k) + p2(k)
      enddo
      S = dot(ptmp,ptmp)
      do k = 0,3
         ptmp(k) = p1(k) - p3(k)
      enddo
      T = dot(ptmp,ptmp)
      do k = 0,3
         ptmp(k) = p1(k) - p4(k)
      enddo
      U = dot(ptmp,ptmp)

      p3m = dsqrt(tridot(p3,p3))
      p4m = dsqrt(tridot(p4,p4))

      p1p2 = dot(p1,p2)
      p1p3 = dot(p1,p3)
      p1p4 = dot(p1,p4)
      p2p3 = dot(p2,p3)
      p2p4 = dot(p2,p4)
      p3p4 = dot(p3,p4)

      me2=ame*ame
      mp2=ampi*ampi

****************************************************************

#ifdef LOOPTOOLS
      A0e = A0(me2)
      A0p = A0(mp2)

      B0ee =B0(me2,me2,0d0)
      B0pp =B0(mp2,mp2,0d0)
      B0p0p=B0(mp2,0d0,mp2)
      B0e0e=B0(me2,0d0,me2)
      B00ee=B0(0d0,me2,me2)
      B0see=B0(s,me2,me2)
      B0spp=B0(s,mp2,mp2)
      B0s00=B0(s,0d0,0d0)

      DB0ee =DB0(me2,me2,0d0)
      DB0pp =DB0(mp2,mp2,0d0)
      DB00ee=DB0(0d0,me2,me2)
      DB0e0e=DB0(me2,0d0,me2)
      DB0p0p=DB0(mp2,0d0,mp2)

      C0eesee =C0i(cc0,me2,me2,s,me2,0d0,me2)
      C0ppspp =C0i(cc0,mp2,mp2,s,mp2,0d0,mp2)
      C0esee00=C0i(cc0,me2,s,me2,me2,0d0,0d0)
      C0epuep =C0i(cc0,me2,mp2,u,me2,0d0,mp2)
      C0pspp00=C0i(cc0,mp2,s,mp2,mp2,0d0,0d0)
      C0eup0ep=C0i(cc0,me2,u,mp2,0d0,me2,mp2)
      C0epte0p=C0i(cc0,me2,mp2,t,me2,0d0,mp2)
      C0etp0ep=C0i(cc0,me2,t,mp2,0d0,me2,mp2)
      C0ees0e0=C0i(cc0,me2,me2,s,0d0,me2,0d0)
      C0pps0p0=C0i(cc0,mp2,mp2,s,0d0,mp2,0d0)

      D0eeppst0e0p=D0i(dd0,me2,me2,mp2,mp2,s,t,0d0,me2,0d0,mp2)
      D0eeppsu0e0p=D0i(dd0,me2,me2,mp2,mp2,s,u,0d0,me2,0d0,mp2)

      D0ts=D0eeppst0e0p
      D0su=D0eeppsu0e0p
#endif

#ifdef COLLIER
      call InitEvent_cll(1)

      call A0_cll(A0e,c1*me2)
      call A0_cll(A0p,c1*mp2)

      call B0_cll(B0ee ,c1*me2,c1*me2,c1*0.d0)
      call B0_cll(B0pp ,c1*mp2,c1*mp2,c1*0.d0)
      call B0_cll(B0p0p,c1*mp2,c1*0.d0,c1*mp2)
      call B0_cll(B0e0e,c1*me2,c1*0.d0,c1*me2)
      call B0_cll(B00ee,c1*0.d0,c1*me2,c1*me2)
      call B0_cll(B0see,c1*s,c1*me2,c1*me2)
      call B0_cll(B0spp,c1*s,c1*mp2,c1*mp2)
      call B0_cll(B0s00,c1*s,c1*0.d0,c1*0.d0)

      call DB0_cll(DB0ee ,c1*me2,c1*me2,c1*0d0)
      call DB0_cll(DB0pp ,c1*mp2,c1*mp2,c1*0d0)
      call DB0_cll(DB00ee,c1*0d0,c1*me2,c1*me2)
      call DB0_cll(DB0e0e,c1*me2,c1*0d0,c1*me2)
      call DB0_cll(DB0p0p,c1*mp2,c1*0d0,c1*mp2)

      call C0_cll(C0eesee ,c1*me2,c1*me2,c1*s,c1*me2,c1*0d0,c1*me2)
      call C0_cll(C0ppspp ,c1*mp2,c1*mp2,c1*s,c1*mp2,c1*0d0,c1*mp2)
      call C0_cll(C0esee00,c1*me2,c1*s,c1*me2,c1*me2,c1*0d0,c1*0d0)
      call C0_cll(C0epuep ,c1*me2,c1*mp2,c1*u,c1*me2,c1*0d0,c1*mp2)
      call C0_cll(C0pspp00,c1*mp2,c1*s,c1*mp2,c1*mp2,c1*0d0,c1*0d0)
      call C0_cll(C0eup0ep,c1*me2,c1*u,c1*mp2,c1*0d0,c1*me2,c1*mp2)
      call C0_cll(C0epte0p,c1*me2,c1*mp2,c1*t,c1*me2,c1*0d0,c1*mp2)
      call C0_cll(C0etp0ep,c1*me2,c1*t,c1*mp2,c1*0d0,c1*me2,c1*mp2)
      call C0_cll(C0ees0e0,c1*me2,c1*me2,c1*s,c1*0d0,c1*me2,c1*0d0)
      call C0_cll(C0pps0p0,c1*mp2,c1*mp2,c1*s,c1*0d0,c1*mp2,c1*0d0)

      call D0_cll(D0eeppst0e0p,c1*me2,c1*me2,c1*mp2,c1*mp2,c1*s,c1*t,
     .                         c1*0d0,c1*me2,c1*0d0,c1*mp2)
      call D0_cll(D0eeppsu0e0p,c1*me2,c1*me2,c1*mp2,c1*mp2,c1*s,c1*u,
     .                         c1*0d0,c1*me2,c1*0d0,c1*mp2)

      D0ts=D0eeppst0e0p
      D0su=D0eeppsu0e0p

#endif

****************************************************************

      e=sqrt(4d0*pi*alpha)
      pi2=pi*pi

      call ffpi(s,ffi,ff2)
      ff2common=ff2
      reffs=real(ffi)
      imffs=aimag(ffi)

******COUNTERTERMS***********************************
      dpsi = -(1D0*(4D0 + 2D0*Log(lambda2/ame**2) +
     &      Log(mu2dim/ame**2) + deltauv))
      dphi = -(2D0*(Log(lambda2/ampi**2) -
     &      1D0*Log(mu2dim/ampi**2) - 1D0*deltauv))
      da = -(1.33333333333333D0*
     &    (Log(mu2dim/ame**2) + deltauv))
      !print*,'dpsi',dpsi
      include '/pions/VirtPhdiv.F'
      photonct=-2.d0*da/16.d0/pi/pi*alpha*4.d0*pi
*****************************************************
      cts=0.d0
      virtual=(0.d0,0.d0)

c      loopf=alpha/4d0/pi*2.d0*bornpp   ! to convert from my notation to fulvio

      !! In my version I have two flags to turn on/off the ISR/FSR, not two final states
      !! I translate this implementation to 'pi' and 'pp'

      if(fs.eq.'pp') then
        ISR = 1
        FSR = 1
      elseif (fs.eq.'pi') then
        ISR=1
        FSR=0
      endif

      if(ISR.gt.0.and.FSR.eq.0) then
            include '/pions/VirtEldiv.F'
            virtual=2.d0*virtel
            virtual=virtual*4.d0*pi*alpha
            virtual=virtual/16.d0/pi/pi
            cts=2.d0*dpsi/16.d0/pi/pi*alpha*4.d0*pi
            !print*,virtual,cts

       !!! This can also do FSR only
c      else if(ISR.eq.0.and.FSR.gt.0) then
c            include '/pions/VirtPidiv.F'
c            virtual=2.d0*virtpi
c            virtual=virtual*4.d0*pi*alpha
c            virtual=virtual/16.d0/pi/pi
c            cts=2.d0*dphi/16.d0/pi/pi*alpha*4.d0*pi

      else if(ISR.gt.0.and.FSR.gt.0) then
            include '/pions/VirtEldiv.F'
            include '/pions/VirtPidiv.F'
            include '/pions/VirtBoxSquarediv.F'
            include '/pions/VirtBoxCrossdiv.F'
            include '/pions/VirtBoxTdiv.F'
            virtual_box=virtboxcross+virtboxsquare
     . +      virtboxt
            virtual=virtpi+virtel
            cts=dpsi+dphi
*************************************************************************************
      !!!===================================================================================!!!
      !!!                                 DISPERSIVE APPROACH
      !!!===================================================================================!!!
      !!!   delta_V = (M_L M_0*)/|M_0|^2
      !!!   M_L = virtual corrections
      !!!   M_0 = Born amplitude
      !!!
      !!!   In the dispersive approach delta_V = delta_V(s',s'')
      !!!   where s' and s'' are the photon virtualities in Box diags
      !!!
      !!!   delta_V   = delta(0,0)
      !!!               - 1/pi Int_0^inf ds'/s' Im F(s') [delta(s',0) + delta(0,s')]
      !!!               - 1/pi^2 Int Int ds'/s' ds''/s'' Im F(s') Im F(s'') delta(s',s'')
      !!!
      !!!   M_0 = M_point x F(s)
      !!!==================================================================================!!!

!POLE POLE !
            repp=real(virtual_box)

! POLE DISP !
            redeltaIRt= 2.d0*s*(me2+mp2-t)*real(C0epte0p)
            redeltaIRu=-2.d0*s*(me2+mp2-u)*real(C0epuep)
            redeltaIR=redeltaIRt+redeltaIRu
            reir=2d0/s*(reffs-1.d0)*redeltaIR

!! Generating RND for Dipsersive MC Integration !!
            call getrnd(psi,2)
            xi=psi(1)           ! will be used for x2 (s')
            csi=psi(2)          ! will be used for y2 (s'')

!     ! POLE DISP !
            delta=1d-5

            spmin=4d0*mp2
            spmax=s_max

            refi=MC_PV(repoledisp,
     .           s,t,u,0d0,spmin,spmax,xi)
            refi=-2.d0/pi*refi

            redd=MC_2D(redispdisp,s,t,u,spmin,
     .           spmax,spmin,spmax,xi,csi)

            delta_disp=-1d-7
            resd=residue(s,t,u,delta_disp)

            redd=redd/pi/pi     ! same as fulvio. Checked 2.04.24

            reboxpp=repp
            reboxpd=refi+reir+resd
            reboxdd = redd

****************** F P U the real part checked with Fulvio. 02.04.24 ********

! ******************Imaginary Part*******************

! POLE-POLE !
            impp=aimag(virtual_box)
            immax=(sqrt(s)-sqrt(lambda2))**2

            imff=ImFsubtrF_C(s,t,u,xi)
            imff=-2.d0/pi*imff

!!! Eventually substitute lim s'->s with d/ds'
! h=1d-4
! x2 = spmin + xi*(s-spmin)
! diff=x2-s
! if(abs(diff).le.h) then
!       FFmder=ImFsubtrF_minus_der(s,t,u,x2,0d0,h)
!       FFder=ImFsubtrF_der(s,t,u,x2,0d0,h)

!       FFmder=FFmder*(s-spmin)

!       print*,'FFmder',FFmder
!       FFder=FFder*(s-spmin)
!       print*,'FFder',FFder

!       imff=FFmder+FFder
!       PRINT*,'STOPPO'f
!       stop
! endif

!!! (2)
!     immax=(sqrt(s)-sqrt(lambda2))**2
!     imbc=MC_1D(box_only_BC_mc,s,t,u,0d0,spmin,
!     .            immax,xi)
!     imbc=(-2d0)/pi*imbc*imffs

                  !!! (4)
            imir=-2d0/pi*imffs*
     .           boxdm_moretti_lambda(p1,p2,p3,p4,s,t,u,lambda2)

!     ! DISP DISP !
!!! int ds'/s' ds''/s'' Im F(s')Im F(s'') Im[Box(s',s'')]
            imdd=MC_2D(imdispdisp,s,t,u,spmin,
     .           s,spmin,s,xi,csi)
            imdd=imdd/pi/pi     ! Checked. Does not depend on lambda

            imboxpp=impp
            imboxpd=imff+imbc+imir
            imboxdd = imdd

******************CHECK OK****************
******************TOTAL*******************

            rebox=(reboxpp+reboxpd+reboxdd)*reffs
            imbox=(imboxpp+imboxpd+imboxdd)*imffs

            virtual =(rebox+imbox)/ff2+virtual

            FSR_disp=0d0
            CTpi_disp=0d0
            if(iFSRdisp.eq.1) then
!     here i multiply for the F_\pi(s)* coming from the tree level
               FSR_disp=FSRdisp(s,xi,csi,4d0*mp2,spmax)
               FSR_disp=FSR_disp*dconjg(ffi)/ff2
               virtual = virtual+FSR_disp

               CTpi_disp=CTpidisp(s,xi,csi,4d0*mp2,spmax)
               CTpi_disp=CTpi_disp*dconjg(ffi)/ff2
               cts = cts+CTpi_disp
            endif
c$$$  print*,impp,imdd
c$$$  print*,imff,imbc,imir
c$$$  print*,deltauv,lambda2,mu2dim
c$$$  print*,redd,reboxpp,reboxpd
c$$$  print*,rebox,imbox
c$$$  print*,imir,imdd

!     virtual = virtual+FSR_disp
ccccccccc In my convention, virtual is the correction divided by the born. So I have to multiply by e^2
ccccccccc for the coupling of the amplitude and divide by the loop factor 16pi^2

            virtual=virtual*4.d0*pi*alpha ! e^2 coupling
            virtual=virtual/16.d0/pi/pi ! loop factor
            virtual=2.d0*real(virtual) ! interference
            cts=2.d0*(cts)/16.d0/pi/pi*alpha*4.d0*pi
!print*,'virtual+cts',virtual+cts
      endif

      call realsoft(p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,mp2,me2,
     &                      rsoft,eps)

      egmin  = eps*sqrt(s)/2.d0
      call soft_integral_pions(egmin,p1,p2,p3,p4,rsofttot,rsoftISR,
     .     rsoftFSR,rsoftIFI)

      if(fs.eq.'pi') then
            rsoft = rsoftISR
      elseif(fs.eq.'pp') then
            rsoft = rsoftISR + rsoftFSR + rsoftIFI
      endif

      virtplussoftpions=virtual+rsoft+cts

      return
      END

      function boxdm_loopt_integrand(s,t,u,lambda2)
      implicit none
      real*8 boxdm_loopt_integrand
      real*8 s,t,u,lambda2
      real*8 C1,C2
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 mp2,me2

      real*8 I1lt,I2lt
      external I1lt
      external I2lt
      real*8 I1,I2

      real*8 k,lognum,logden,logd0

      mp2=ampi**2
      me2=ame**2

      C1=4d0*s*(me2+mp2-t)*(me2*(2d0*mp2-t+u-lambda2)
     . +2d0*mp2**2-2d0*mp2*t+t*(t-u+lambda2))
      C1=C1/(4d0*mp2*s-s**2+(t-u)**2)

      C2=-4d0*s*(me2-t)*(me2+mp2-t)
      C2=C2/(4d0*mp2*s-s**2+(t-u)**2)

      I1=I1lt(s,t,lambda2)
      I2=I2lt(s,t,lambda2)

      boxdm_loopt_integrand=I1*C1+I2*C2
      return
      end

      function I1lt(s,t,lambda2)      ! int ds'/s' Im(D0)
      implicit none
      real*8 s,t,lambda2,lambda
      real*8 I1lt
      real*8 logs
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 mp2,me2
      real*8 betap2
      real*8 f,lognum,logden,logd0,k
      mp2=ampi**2
      me2=ame**2
      betap2=1d0-4d0*mp2/s

      lambda=sqrt(lambda2)
      f = (mp2-me2-t)**2 - 4.d0*me2*t
      logs=log(lambda/sqrt(s))+log(1d0)
     .    -log(betap2*s/(4d0*mp2))


      k = me2*me2 + (mp2-t)**2 - 2d0*me2*(mp2+t)
      lognum= me2+mp2-t+sqrt(k)
      logden= me2+mp2-t-sqrt(k)
      logd0=0.5d0*log(lognum/(logden))

      I1lt=-2d0*pi/(sqrt(f)*s)*logs*logd0

      return
      end


      function I2lt(s,t,lambda2)      ! int ds' Im(D0)
      implicit none
      real*8 s,t,lambda2,lambda
      real*8 I2lt
      real*8 logs
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 mp2,me2
      real*8 betap2
      real*8 f,lognum,logden,logd0,k
      mp2=ampi**2
      me2=ame**2

      betap2=1d0-4d0*mp2/s

      lambda=sqrt(lambda2)
      f = (mp2-me2-t)**2 - 4.d0*me2*t
      logs=log(lambda/sqrt(s))+log((1d0)/betap2)
      k = me2*me2 + (mp2-t)**2 - 2d0*me2*(mp2+t)
      lognum= me2+mp2-t+sqrt(k)
      logden= me2+mp2-t-sqrt(k)
      logd0=0.5d0*log(lognum/(logden))

      I2lt=-2d0*pi/(sqrt(f))*logs*logd0

      return
      end

      function MC_1D(func,s,t,u,y2,min,max,xi)
      implicit none
      real*8 func
      external func
      real*8 s,t,u,y2,x2
      real*8 min,max,xi
      real*8 MC_1D
      real*8 funci

      x2 = min + xi*(max-min)

      funci = func(s,t,u,x2,y2)*(max-min)
      MC_1D=funci
      return
      end

      function MC_1D_c(func,s,t,u,y2,min,max,xi)
      implicit none
      complex*16 func
      external func
      real*8 s,t,u,y2,x2
      real*8 min,max,xi
      complex*16 MC_1D_C
      complex*16 funci

      x2 = min + xi*(max-min)

      funci = func(s,t,u,x2,y2)*(max-min)
      MC_1D_C=funci
      return
      end

      function mc_2d(func,s,t,u,min1,max1,min2,max2,xi,csi)
      implicit none
      real*8 s,t,u,spmin,spmax
      real*8 xi,csi
      real*8 mc_2d

      real*8 x2,y2

      real*8 func,funci
      external func

      real*8 max1,max2,min1,min2

      x2=min1+xi*(max1-min1)
      y2=min2+csi*(max2-min2)

      funci=func(s,t,u,x2,y2)
      mc_2d=funci*(max1-min1)*(max2-min2)
      return
      end

      function MC_PV(func,s,t,u,y2,spmin,spmax,xi)
      implicit none
      real*8 x2,x2s,xi,y2
      real*8 s,t,u,spmin,spmax
      real*8 func
      external func
      real*8 MC_PV


      real*8 bornpp
      common/pions/bornpp

      real*8 funcim,funcip
      real*8 funci,funci_rem

      real*8 pvspmin,pvspmax,imin,imax


      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi

      real*8 delta,smdel,spdel
      delta=1d-5
      smdel=s-delta
      spdel=s+delta
! Ucci's patch
c      if(smdel.lt.4d0*ampi**2) then
c         smdel = (s- 4d0*ampi**2)
c         spdel =(s+4d0*ampi**2)
c      endif
!  Carlo's patch
      if(smdel.lt.4d0*ampi**2) then
         delta = (s - 4d0*ampi**2)*0.5d0
         smdel = s - delta
         spdel = s + delta
      endif

      if(abs(spmax-s).lt.abs(spmin-s)) then
            pvspmax=spmax
            pvspmin=s-abs(pvspmax-s)
            imin=spmin
            imax=pvspmin
      else
            pvspmin=spmin
            pvspmax=s+abs(pvspmin-s)
            imin=pvspmax
            imax=spmax
      endif
      x2 = pvspmin + xi*(smdel-pvspmin)
      x2s = s + (s-x2)

      funcim=func(s,t,u,x2,y2)
      funcip=func(s,t,u,x2s,y2)

      funci=(funcip+funcim)*(smdel-pvspmin)

      x2 = imin + xi*(imax-imin)
      funci_rem = func(s,t,u,x2,y2)
      funci_rem=funci_rem*(imax-imin)

      MC_PV=funci+funci_rem

      return
      end

      function ReDispDisp(s,t,u,x2,y2)
      implicit none
      real*8 redispdisp
      complex*16 boxddt,boxddu
      real*8 boxdd

      real*8 s,t,u,x2,y2
      external ffpi
      real*8 imffx2,imffy2,ff2
      complex*16 ffx2,ffy2

      complex*16 BoxDM
      external boxdm

      character*10 filename
      integer is

      call ffpi(x2,ffx2,ff2)
      imffx2=aimag(ffx2)
      call ffpi(y2,ffy2,ff2)
      imffy2=aimag(ffy2)


      boxddt=BoxDM(s,t,u,x2,y2)
      boxddu=-BoxDM(s,u,t,x2,y2)

      boxdd=real(boxddt+boxddu)
      redispdisp=boxdd*imffx2*imffy2/x2/y2
!      print*,boxdd,imffx2,imffy2
      return
      end


      function ImDispDisp(s,t,u,x2,y2)
      implicit none
      real*8 imdispdisp
      complex*16 boxddt,boxddu
      real*8 boxdd

      real*8 s,t,u,x2,y2
      external ffpi
      real*8 imffx2,imffy2,ff2
      complex*16 ffx2,ffy2

      complex*16 BoxDM
      external boxdm

      character*10 filename
      integer is

      call ffpi(x2,ffx2,ff2)
      imffx2=aimag(ffx2)
      call ffpi(y2,ffy2,ff2)
      imffy2=aimag(ffy2)

      boxddt=(0d0,0d0)
      boxddu=(0d0,0d0)
      boxddt=BoxDM(s,t,u,x2,y2)
      boxddu=-BoxDM(s,u,t,x2,y2)
      boxdd=aimag(boxddt+boxddu)

      imdispdisp=boxdd*imffx2*imffy2/x2/y2
!      print*,boxdd,imffx2,imffy2
      return
      end
***
      function Repoledisp(s,t,u,x2,y2)
      implicit none
      real*8 s,t,u,x2,y2
      real*8 Repoledisp

      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi

      complex*16 boxdm
      external boxdm

      real*8 BoxDIR
      external boxdir

      real*8 ff2
      complex*16 ff
      external ffpi

      complex*16 Boxpdt,boxpdu
      real*8 BoxpdtIR,boxpd,BoxpduIR

      real*8 boxpd_t,boxpd_u

      character*10 filename
      integer is


      boxpdt=Boxdm(s,t,u,x2,y2)
      BoxpdtIR=boxdir(s,t,u,x2,y2)

      boxpdu=-Boxdm(s,u,t,x2,y2)
      BoxpduIR=-boxdir(s,u,t,x2,y2)

      boxpd_t=real(boxpdt-BoxpdtIR)
      boxpd_u=real(boxpdu-BoxpduIR)

      call ffpi(x2,ff,ff2)

      boxpd=boxpd_t+boxpd_u
      boxpd=boxpd*aimag(ff)
      repoledisp=boxpd/x2
      return
      end

      function ImFsubtrF_c(s,t,u,x)
      implicit none
      real*8 ImFsubtrF_C
      complex*16 boxddt,boxddu
      real*8 boxdd

      real*8 s,t,u,x
      external ffpi
      real*8 imffsp,imffs,ff2
      complex*16 ffsp,ffs

      complex*16 BoxDM
      external boxdm

      real*8 smsp,sp,norm
      real*8 ame,ammu,convfac,alpha,pi,ampi
      real*8 mp2
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      mp2=ampi**2
      norm = log((s-4.d0*mp2)/(1d-12))
      smsp= (s-4.d0*mp2)*exp(-x*norm)
      sp= s - smsp
      call ffpi(sp,ffsp,ff2)
      imffsp=aimag(ffsp)
      call ffpi(s,ffs,ff2)
      imffs=aimag(ffs)

      boxddt=BoxDM(s,t,u,sp,0d0)
      boxddu=-BoxDM(s,u,t,sp,0d0)

      boxdd=aimag(boxddt+boxddu)*smsp*norm
      ImFsubtrF_c=(imffsp-imffs)*boxdd/sp


      return
      end

      function box_only_BC_mc(s,t,u,x2,y2)
      implicit none
      real*8 box_only_BC_mc
      real*8 s,t,u,x2,y2

      real*8 BoxDM_nod0
      external BoxDM_nod0

      external ffpi
      real*8 ff2p
      complex*16 ffp

      real*8 Boxsubt,boxsubu

      boxsubt=BoxDM_nod0(s,t,u,x2,0d0)
      boxsubu=-BoxDM_nod0(s,u,t,x2,0d0)

      box_only_BC_mc=(boxsubt+boxsubu)/x2

      return
      End

      function BoxDM(s,t,u,x2,y2)

#ifdef COLLIER
      use collier
#endif

      implicit none
      real*8 s,t,u,x2,y2
      complex*16 B0ppy,B0pxp,B0sxy
      complex*16 C0eptexp,C0epteyp,C0eesxey,C0ppsxpy
      complex*16 d0eeppstxeyp
      real*8 deltauv,lambda2,mu2dim,m2min
      common/loop/deltauv,lambda2,mu2dim,m2min
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi

      complex*16 BoxDM,C0eptepy,BoxDIR

      real*8 me2,mp2

      double complex c1
      parameter (c1 = (1.d0,0.d0))

#ifdef LOOPTOOLS
#include "looptools.h"
      call clearcache
      me2=ame**2
      mp2=ampi**2

      B0ppy=B0(mp2,mp2,y2)
      B0pxp=B0(mp2,x2,mp2)
      B0sxy=B0(s,x2,y2)

      C0eptexp=C0i(cc0,me2,mp2,t,me2,x2,mp2)
      C0epteyp=C0i(cc0,me2,mp2,t,me2,y2,mp2)
      C0eesxey=C0i(cc0,me2,me2,s,x2,me2,y2)
      C0ppsxpy=C0i(cc0,mp2,mp2,s,x2,mp2,y2)

      d0eeppstxeyp=D0i(dd0,me2,me2,mp2,mp2,s,t,x2,me2,
     . y2,mp2)

#endif

#ifdef COLLIER
      call InitEvent_cll(1)
      me2=ame**2
      mp2=ampi**2

      call B0_cll(B0ppy,c1*mp2,c1*mp2,c1*y2)
      call B0_cll(B0pxp,c1*mp2,c1*x2,c1*mp2)
      call B0_cll(B0sxy,c1*s,c1*x2,c1*y2)

      call C0_cll(C0eptexp,c1*me2,c1*mp2,c1*t,c1*me2,c1*x2,c1*mp2)
      call C0_cll(C0epteyp,c1*me2,c1*mp2,c1*t,c1*me2,c1*y2,c1*mp2)
      call C0_cll(C0eesxey,c1*me2,c1*me2,c1*s,c1*x2,c1*me2,c1*y2)
      call C0_cll(C0ppsxpy,c1*mp2,c1*mp2,c1*s,c1*x2,c1*mp2,c1*y2)

      call D0_cll(d0eeppstxeyp,c1*me2,c1*me2,c1*mp2,c1*mp2,c1*s,c1*t,
     .                c1*x2,c1*me2,c1*y2,c1*mp2)

#endif
c$$$      print*,me2,mp2
c$$$      print*,B0ppy,B0pxp,B0sxy,C0eptexp,C0epteyp,C0eesxey,C0ppsxpy
c$$$     >     ,d0eeppstxeyp

      BoxDM = -((s*2d0*(C0eptexp*(ame**2 - t*1d0)*
     &         (ame**2 + ampi**2 - t*1d0)*2d0 +
     &        C0epteyp*(ame**2 - t*1d0)*
     &         (ame**2 + ampi**2 - t*1d0)*2d0 -
     &        C0eesxey*2d0*
     &         (t*(t - u*1d0) + ampi**4*2d0 - ampi**2*t*2d0 +
     &           ame**2*(u - t*1d0 + ampi**2*2d0)) -
     &        d0eeppstxeyp*(ame**2 + ampi**2 - t*1d0)*2d0*
     &         (t*(t + x2 + y2 - u*1d0) + ampi**4*2d0 -
     &           ampi**2*t*2d0 +
     &           ame**2*(u - t*1d0 - x2*1d0 - y2*1d0 +
     &              ampi**2*2d0)) +
     &        ((ampi**4 + ame**2*(t + u) - ame**4*1d0 -
     &             t*u*1d0)*(B0ppy + B0pxp - B0sxy*2d0))/
     &         (-(s*1d0) + ampi**2*4d0) -
     &        (C0ppsxpy*1d0*
     &           (t*(u*(x2 + y2) - s**2*1d0 - s*t*3d0) +
     &             ampi**2*t*2d0*(s*5d0 + t*8d0) +
     &             ame**4*
     &              (x2 + y2 - s*3d0 + ampi**2*1.6d1) +
     &             ampi**4*
     &              (s - x2*1d0 - y2*1d0 - t*1.6d1) +
     &             ame**2*
     &              (s**2 - (t + u)*(x2 + y2)*1d0 +
     &                s*t*6d0 + ampi**4*1.6d1 -
     &                ampi**2*2d0*(s*5d0 + t*1.6d1))))/
     &         (-(s*1d0) + ampi**2*4d0)))/
     &    (-(s**2*1d0) + (t - u*1d0)**2 + ampi**2*s*4d0))
      return
      end


      function BoxDM_nod0(s,t,u,x2,y2)

#ifdef COLLIER
      use collier
#endif

      implicit none
      real*8 s,t,u,x2,y2
      complex*16 B0ppy,B0pxp,B0sxy
      complex*16 C0eptexp,C0epteyp,C0eesxey,C0ppsxpy
      complex*16 d0eeppstxeyp
      real*8 deltauv,lambda2,mu2dim,m2min
      common/loop/deltauv,lambda2,mu2dim,m2min
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi

      complex*16 BoxDM,C0eptepy,BoxDIR
      real*8 BoxDM_nod0
      real*8 me2,mp2

      double complex c1
      parameter (c1 = (1.d0,0.d0))

#ifdef LOOPTOOLS
#include "looptools.h"
      call clearcache
      me2=ame**2
      mp2=ampi**2

      B0ppy=B0(mp2,mp2,y2)
      B0pxp=B0(mp2,x2,mp2)
      B0sxy=B0(s,x2,y2)

      C0eptexp=C0i(cc0,me2,mp2,t,me2,x2,mp2)
      C0epteyp=C0i(cc0,me2,mp2,t,me2,y2,mp2)
      C0eesxey=C0i(cc0,me2,me2,s,x2,me2,y2)
      C0ppsxpy=C0i(cc0,mp2,mp2,s,x2,mp2,y2)
#endif

#ifdef COLLIER
      call InitEvent_cll(1)
      me2=ame**2
      mp2=ampi**2

      call B0_cll(B0ppy,c1*mp2,c1*mp2,c1*y2)
      call B0_cll(B0pxp,c1*mp2,c1*x2,c1*mp2)
      call B0_cll(B0sxy,c1*s,c1*x2,c1*y2)

      call C0_cll(C0eptexp,c1*me2,c1*mp2,c1*t,c1*me2,c1*x2,c1*mp2)
      call C0_cll(C0epteyp,c1*me2,c1*mp2,c1*t,c1*me2,c1*y2,c1*mp2)
      call C0_cll(C0eesxey,c1*me2,c1*me2,c1*s,c1*x2,c1*me2,c1*y2)
      call C0_cll(C0ppsxpy,c1*mp2,c1*mp2,c1*s,c1*x2,c1*mp2,c1*y2)
#endif

      d0eeppstxeyp=(0d0,0d0)
      BoxDM = -((s*2d0*(C0eptexp*(ame**2 - t*1d0)*
     &         (ame**2 + ampi**2 - t*1d0)*2d0 +
     &        C0epteyp*(ame**2 - t*1d0)*
     &         (ame**2 + ampi**2 - t*1d0)*2d0 -
     &        C0eesxey*2d0*
     &         (t*(t - u*1d0) + ampi**4*2d0 - ampi**2*t*2d0 +
     &           ame**2*(u - t*1d0 + ampi**2*2d0)) -
     &        d0eeppstxeyp*(ame**2 + ampi**2 - t*1d0)*2d0*
     &         (t*(t + x2 + y2 - u*1d0) + ampi**4*2d0 -
     &           ampi**2*t*2d0 +
     &           ame**2*(u - t*1d0 - x2*1d0 - y2*1d0 +
     &              ampi**2*2d0)) +
     &        ((ampi**4 + ame**2*(t + u) - ame**4*1d0 -
     &             t*u*1d0)*(B0ppy + B0pxp - B0sxy*2d0))/
     &         (-(s*1d0) + ampi**2*4d0) -
     &        (C0ppsxpy*1d0*
     &           (t*(u*(x2 + y2) - s**2*1d0 - s*t*3d0) +
     &             ampi**2*t*2d0*(s*5d0 + t*8d0) +
     &             ame**4*
     &              (x2 + y2 - s*3d0 + ampi**2*1.6d1) +
     &             ampi**4*
     &              (s - x2*1d0 - y2*1d0 - t*1.6d1) +
     &             ame**2*
     &              (s**2 - (t + u)*(x2 + y2)*1d0 +
     &                s*t*6d0 + ampi**4*1.6d1 -
     &                ampi**2*2d0*(s*5d0 + t*1.6d1))))/
     &         (-(s*1d0) + ampi**2*4d0)))/
     &    (-(s**2*1d0) + (t - u*1d0)**2 + ampi**2*s*4d0))
      BoxDM_nod0=aimag(boxdm)
      return
      end

      function BoxDIR(s,t,u,x2,y2)

#ifdef COLLIER
      use collier
#endif

      implicit none
      real*8 BoxDIR
      real*8 s,t,u,y2,x2
      real*8 me2,mp2
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      complex*16 C0eptepy
      double complex c1
      parameter (c1 = (1.d0,0.d0))

#ifdef LOOPTOOLS
#include "looptools.h"
      call clearcache
      me2=ame**2
      mp2=ampi**2
      C0eptepy=C0i(cc0,me2,mp2,t,me2,y2,mp2)
#endif

#ifdef COLLIER
      call InitEvent_cll(1)
      me2=ame**2
      mp2=ampi**2
      call C0_cll(C0eptepy,c1*me2,c1*mp2,c1*t,c1*me2,c1*y2,c1*mp2)
#endif

      BoxDIR=2.d0*s*(me2+mp2-t)*real(C0eptepy)
      BoxDIR=BoxDIR/(s-x2)

      return
      End



      function I1col(s,t,lambda2)      ! int ds'/s' Im(D0)
      implicit none
      real*8 s,t,lambda2,lambda
      real*8 I1col
      real*8 logs,logs2

      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 mp2,me2
      real*8 betap2
      real*8 f,lognum,logden,logd0,k
      mp2=ampi**2
      me2=ame**2
      betap2=1d0-4d0*mp2/s

      lambda=sqrt(lambda2)
      f = (mp2-me2-t)**2 - 4.d0*me2*t
      logs=log(lambda/sqrt(s))+log(1d0+sqrt(-s*t/f))
     .    -log(betap2*s/(4d0*mp2))

      k = me2*me2 + (mp2-t)**2 - 2d0*me2*(mp2+t)
      ! lognum= me2+mp2-t+sqrt(k)
      ! logden= me2+mp2-t-sqrt(k)
      logs=logs*log((mp2-t)/s)

      logs2=log(lambda/sqrt(s))+log(2d0)
     .    -log(betap2*s/(4d0*mp2))

      logs2=logs2*log(s/(ame*ampi))

      I1col=-2d0*pi/(sqrt(f)*s)*(logs+logs2)


      return
      end

      function I2col(s,t,lambda2)      ! int ds'/s' Im(D0)
      implicit none
      real*8 s,t,lambda2,lambda
      real*8 I2col
      real*8 logs,logs2

      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 mp2,me2
      real*8 betap2
      real*8 f,lognum,logden,logd0,k
      mp2=ampi**2
      me2=ame**2
      betap2=1d0-4d0*mp2/s

      lambda=sqrt(lambda2)
      f = (mp2-me2-t)**2 - 4.d0*me2*t
      logs=log(lambda/sqrt(s))+log(1d0+sqrt(-s*t/f))
     .    -log(betap2)

      ! k = me2*me2 + (mp2-t)**2 - 2d0*me2*(mp2+t)
      lognum= me2+mp2-t+sqrt(k)
      ! logden= me2+mp2-t-sqrt(k)
      logs=logs*log(lognum/s)

      logs2=log(lambda/sqrt(s))+log(2d0)
     .    -log(betap2)

      logs2=logs2*log(s/(ame*ampi))

      I2col=-2d0*pi/(sqrt(f))*(logs+logs2)


      return
      end

      function boxdm_col_integrand(s,t,u,lambda2)
      implicit none
      real*8 boxdm_col_integrand
      real*8 s,t,u,lambda2
      real*8 C1,C2

      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 mp2,me2

      real*8 I1col,I2col
      external I1col
      external I2col
      real*8 I1,I2

      mp2=ampi**2
      me2=ame**2


      C1=4d0*s*(me2+mp2-t)*(me2*(2d0*mp2-t+u-lambda2)
     . +2d0*mp2**2-2d0*mp2*t+t*(t-u+lambda2))
      C1=C1/(4d0*mp2*s-s**2+(t-u)**2)

      C2=-4d0*s*(me2-t)*(me2+mp2-t)
      C2=C2/(4d0*mp2*s-s**2+(t-u)**2)

      I1=I1col(s,t,lambda2)
      I2=I2col(s,t,lambda2)

      boxdm_col_integrand=I1*C1+I2*C2
      return
      end


      function boxdm_col_tot(s,t,u,lambda2)
      implicit none
      real*8 boxdm_col_tot
      real*8 s,t,u,lambda2
      real*8 C1,C2

      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 mp2,me2

      real*8 I1col,I2col
      external I1col
      external I2col
      real*8 I1,I2

      real*8 k,lognum,logden,logd0,betap
      real*8 I1t,I1u,I2t,I2u,It,Iu,C1t,C1u,C2t,C2u
      real*8 Ce,cC,cp

      mp2=ampi**2
      me2=ame**2

      betap=sqrt(1d0-4d0*mp2/s)

      Ce=s*(u-t)*pi/(s)*log(ame**2/s)
      Cp=(2d0*mp2-s)*(t-u)*pi/(s*betap)*log((1d0-betap)/(1d0+betap))
      CC=(Ce+Cp)*log(1d0-betap**2)*4d0*s/(4d0*mp2*s-s**2+(t-u)**2)

      C1t=4d0*s*(me2+mp2-t)*(me2*(2d0*mp2-t+u-lambda2)
     . +2d0*mp2**2-2d0*mp2*t+t*(t-u+lambda2))
      C1t=C1t/(4d0*mp2*s-s**2+(t-u)**2)

      C2t=-4d0*s*(me2-t)*(me2+mp2-t)
      C2t=C2t/(4d0*mp2*s-s**2+(t-u)**2)

      I1t=I1col(s,t,lambda2)
      I2t=I2col(s,t,lambda2)

      It=I1t*C1t+I2t*C2t

      C1u=4d0*s*(me2+mp2-u)*(me2*(2d0*mp2-u+t-lambda2)
     . +2d0*mp2**2-2d0*mp2*u+u*(u-t+lambda2))
      C1u=C1u/(4d0*mp2*s-s**2+(t-u)**2)

      C2u=-4d0*s*(me2-u)*(me2+mp2-u)
      C2u=C2u/(4d0*mp2*s-s**2+(t-u)**2)

      I1u=I1col(s,u,lambda2)
      I2u=I2col(s,u,lambda2)

      Iu=I1u*C1u+I2u*C2u
      boxdm_col_tot=-Cc+It-Iu
      return
      end


      function I1colv2(s,t,lambda2)      ! int ds'/s' Im(D0)
      implicit none
      real*8 s,t,lambda2,lambda
      real*8 I1colv2
      real*8 logs,logs2

      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 mp2,me2
      real*8 betap2
      real*8 f,lognum,logden,logd0,k
      mp2=ampi**2
      me2=ame**2
      betap2=1d0-4d0*mp2/s

      lambda=sqrt(lambda2)
      f = (mp2-me2-t)**2 - 4.d0*me2*t
      logs=log(lambda/sqrt(s))+log(1d0+sqrt(-s*t/f))
     .    -log(betap2*s/(4d0*mp2))

      k = me2*me2 + (mp2-t)**2 - 2d0*me2*(mp2+t)
      lognum= me2+mp2-t+sqrt(k)
      logden= me2+mp2-t-sqrt(k)
      logs=logs*0.5d0*log(lognum/logden)


      I1colv2=-2d0*pi/(sqrt(f)*s)*(logs)


      return
      end

      function I2colv2(s,t,lambda2)      ! int ds'/s' Im(D0)
      implicit none
      real*8 s,t,lambda2,lambda
      real*8 I2colv2
      real*8 logs,logs2

      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 mp2,me2
      real*8 betap2
      real*8 f,lognum,logden,logd0,k
      mp2=ampi**2
      me2=ame**2
      betap2=1d0-4d0*mp2/s

      lambda=sqrt(lambda2)
      f = (mp2-me2-t)**2 - 4.d0*me2*t
      logs=log(lambda/sqrt(s))+log(1d0+sqrt(-s*t/f))
     .    -log(betap2)

      k = me2*me2 + (mp2-t)**2 - 2d0*me2*(mp2+t)
      lognum= me2+mp2-t+sqrt(k)
      logden= me2+mp2-t-sqrt(k)

      I2colv2=-2d0*pi/(sqrt(f))*logs*0.5d0*log(lognum/logden)


      return
      end


      function boxdm_col_totv2(s,t,u,lambda2)
      implicit none
      real*8 boxdm_col_totv2
      real*8 s,t,u,lambda2
      real*8 C1,C2

      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 mp2,me2

      real*8 I1colv2,I2colv2
      external I1colv2
      external I2colv2
      real*8 I1,I2

      real*8 k,lognum,logden,logd0,betap
      real*8 I1t,I1u,I2t,I2u,It,Iu,C1t,C1u,C2t,C2u
      real*8 Ce,cC,cp

      mp2=ampi**2
      me2=ame**2

      betap=sqrt(1d0-4d0*mp2/s)

      Ce=s*(u-t)*pi/(s)*log(ame**2/s)
      Cp=(2d0*mp2-s)*(t-u)*pi/(s*betap)*log((1d0-betap)/(1d0+betap))
      CC=(Ce+Cp)*log(1d0-betap**2)*4d0*s/(4d0*mp2*s-s**2+(t-u)**2)

      C1t=4d0*s*(me2+mp2-t)*(me2*(2d0*mp2-t+u-lambda2)
     . +2d0*mp2**2-2d0*mp2*t+t*(t-u+lambda2))
      C1t=C1t/(4d0*mp2*s-s**2+(t-u)**2)

      C2t=-4d0*s*(me2-t)*(me2+mp2-t)
      C2t=C2t/(4d0*mp2*s-s**2+(t-u)**2)

      I1t=I1colv2(s,t,lambda2)
      I2t=I2colv2(s,t,lambda2)

      It=I1t*C1t+I2t*C2t

      C1u=4d0*s*(me2+mp2-u)*(me2*(2d0*mp2-u+t-lambda2)
     . +2d0*mp2**2-2d0*mp2*u+u*(u-t+lambda2))
      C1u=C1u/(4d0*mp2*s-s**2+(t-u)**2)

      C2u=-4d0*s*(me2-u)*(me2+mp2-u)
      C2u=C2u/(4d0*mp2*s-s**2+(t-u)**2)

      I1u=I1colv2(s,u,lambda2)
      I2u=I2colv2(s,u,lambda2)

      Iu=I1u*C1u+I2u*C2u
      boxdm_col_totv2=-Cc+It-Iu
      return
      end


       function im_deltapd_bern(s,z,l2) !!!!! Second part of (3.8) Colangelo
      implicit none
      real*8 im_deltapd_bern
      real*8 l2,s,beta,z
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*16 mp2,me2
      real*8 f1,f2,f3,f4,f5,f6,f60,f61,f62,f63,f7,f70,f71,f72
      real*8 ff2
      complex*16 ff
      external ffpi

      mp2=ampi**2
      beta=sqrt(1d0-4d0*mp2/s)

      f1= 2d0*log(4d0*l2/s)*log((1d0-beta*z)/(1d0+beta*z))
      f2=4d0*log((1d0-beta*z)/2d0)
     .*log(1d0-beta*z+sqrt(1d0-2d0*beta*z+beta**2))
     .-4d0*(log(1d0-beta*z))**2
      f3=-4d0*log((1d0+beta*z)/2d0)
     . *log(1d0+beta*z+sqrt(1d0+2d0*beta*z+beta**2))
     . + 4d0*(log(1d0+beta*z))**2

      f60=2d0*log(1-beta**2)/((1d0-z**2)*beta**2)
      f61=z*(1d0+beta**2)*log((1d0-beta)/(1d0+beta))
      f62=-log((1d0-beta*z)/(1d0+beta*z))
      f63=2d0*beta*z*log((1d0-beta**2*z**2)/(1d0-beta**2))
      f6=f60*(f61+f62+f63)

      f70=-2d0*log((1d0-beta*z)/(1d0+beta*z))/(1d0-z**2)
      f71=(2d0*z**2-1d0)*log((1d0-beta**2)/beta**2)
      f72=2d0*log(beta)
      f7=f70*(f71+f72)
      im_deltapd_bern=f1+f2+f3+f6+f7

      call ffpi(s,ff,ff2)
      im_deltapd_bern=im_deltapd_bern*alpha/pi*(aimag(ff)**2)/ff2
      return
      end
*
*
* Eq. 272 Francesco's notes
*
      function InonLogMoretti(a,b)
      implicit none
      real*8 a,b
      real*8  InonLogMoretti

      real*8 sqta,sqlog

      sqta=sqrt(a)
      sqlog= sqrt((a-a*b)/(a-b))
      InonLogMoretti=1.d0/sqrt((b-1.d0)*(b-a))
     +              * (log(1.d0+sqlog) -log(sqlog-1.d0))
      return
      end

      function ILogMoretti(a,b)
      implicit none
      real*8 a,b
      complex*16 ILogMoretti

      real*8 sqta
      complex*16 pref
      real*8 z

      complex*16 myli2
      external myli2
      complex*16 dilog1,dilog2,dilog3,dilog4,dilog5,dilog6
      complex*16 logquadro1,logquadro2,logquadro3,logquadro4
      complex*16 term3,term41,term4

      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi

      sqta=sqrt(a)
      ILogMoretti=1d0/(sqrt(b-a)*sqrt(b-1d0))

      z=sqrt((b-a)/(b-1))
      dilog1 = -1d0*myli2(dcmplx((1d0-1d0/sqta*z)/(1d0+z)))
      dilog2 = +1d0*myli2(dcmplx((1d0+1d0/sqta*z)/(1d0+z)))
      dilog3 = +1d0*myli2(dcmplx((1d0-z)/(1d0-1d0/sqta*z)))
      dilog4 = -1d0*myli2(dcmplx((1d0-z)/(1d0+1d0/sqta*z)))

      dilog5=-1d0*myli2(dcmplx(1d0-sqta/z))
      logquadro1=-1d0/2d0*(log(sqta/z))**2

      term3=2d0*(dilog5+logquadro1)

      dilog6=-6d0*myli2(dcmplx(1d0/(1d0+(1d0/(sqta/z)))))
      term41=(2d0*pi**2-3d0*log(1d0+1d0/(sqta/z))**2)
      term4=-1d0/3d0*(dilog6+term41)

      logquadro2=1d0/2d0*log((1d0-z/sqta)/(1d0-z))**2
      logquadro3=-1d0/2d0*log((1d0+z/sqta)/(1d0-z))**2
      logquadro4= log(b-1d0)*log((1d0+z/sqta)/(1d0-z/sqta))
      ILogMoretti=ILogMoretti*(dilog1+dilog2+dilog3+dilog4
     .      +term3+term4+logquadro2+logquadro3+logquadro4)

      return
      end
*
* Eq. 291 Francesco's notes (without changing the original terms with C0)
*
      real*8 function boxdm_moretti(p1,p3,p4,s,t,u)
#ifdef COLLIER
      use collier
#endif
      implicit none
      real*8 p1(0:3),p3(0:3),p4(0:3)
      real*8 s,t,u,lambda2
      real*8 C1ovspt,Cspovspt,C1ovspu,Cspovspu

      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 mp2,me2

      real*8 Itot1ovspt,Itot1ovspu
      real*8 Itot1ovsp_dimreg
      external Itot1ovsp_dimreg
      real*8 Itotspovspt,Itotspovspu
      real*8 Itotspovsp_dimreg
      external Itotspovsp_dimreg

      real*8 k,lognum,logden,logd0,betap
      real*8 I1t,I1u,I2t,I2u,It,Iu,C1t,C1u,C2t,C2u
      real*8 Ce,cC,cp

      real*8 I1lt, I2lt
      external i2lt,i1lt
      integer i
      complex*16 D0ltd
      double complex c1,D0eeppst0e0p
      parameter (c1 = (1.d0,0.d0))
      real*8 Itlt,iult

      mp2=ampi**2
      me2=ame**2

      betap=sqrt(1d0-4d0*mp2/s)

      Ce=s*(u-t)*pi/(s)*log(ame**2/s)
      Cp=(2d0*mp2-s)*(t-u)*pi/(s*betap)*log((1d0-betap)/(1d0+betap))
      CC=(Ce+Cp)*log(1d0-betap**2)*4d0*s/(4d0*mp2*s-s**2+(t-u)**2)

      C1ovspt=4d0*s*(me2+mp2-t)*(me2*(2d0*mp2-t+u)
     . +2d0*mp2**2-2d0*mp2*t+t*(t-u))
      C1ovspt=C1ovspt/(4d0*mp2*s-s**2+(t-u)**2)

      Cspovspt=-4d0*s*(me2-t)*(me2+mp2-t)
      Cspovspt=Cspovspt/(4d0*mp2*s-s**2+(t-u)**2)

      Itot1ovspt=Itot1ovsp_dimreg(p1,p3,s)
      Itotspovspt=Itotspovsp_dimreg(p1,p3,s)


      It=Itot1ovspt*C1ovspt+Itotspovspt*Cspovspt


      C1ovspu=4d0*s*(me2+mp2-u)*(me2*(2d0*mp2-u+t)
     . +2d0*mp2**2-2d0*mp2*u+u*(u-t))
      C1ovspu=C1ovspu/(4d0*mp2*s-s**2+(t-u)**2)

      Cspovspu=-4d0*s*(me2-u)*(me2+mp2-u)
      Cspovspu=Cspovspu/(4d0*mp2*s-s**2+(t-u)**2)

      Itot1ovspu=Itot1ovsp_dimreg(p1,p4,s)
      Itotspovspu=Itotspovsp_dimreg(p1,p4,s)


      Iu=Itot1ovspu*C1ovspu+Itotspovspu*Cspovspu

      boxdm_moretti=-Cc+It-Iu
      return
      end
*
* Eq. 289 Francesco'1 notes (term proportional to 1/s')
*
      real*8 function Itot1ovsp_dimreg(p1,p3,s)
      implicit none
      real*8 p1(0:3),p3(0:3)
      real*8 s
      real*8 Itot_nonlog, Itot_log
      external Itot_nonlog, Itot_log
      real*8 A_1ovsp,gamma,B_1ovsp
      real*8 betapi2
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 oneoeps,mudimreg2,eulergamma
      common/IRdimreg/oneoeps,mudimreg2,eulergamma
      real*8 tridot
      external tridot
      betapi2=1d0-4d0*ampi**2/s
      A_1ovsp= 1.d0/8.d0/s
      B_1ovsp= 1.d0/8.d0/s * (2d0*log(betapi2/(1.d0-betapi2))
     +     +log(s/mudimreg2) -2.d0*log(4.d0*pi))
      gamma= log(4.d0*pi) + eulergamma

*
      Itot1ovsp_dimreg= -oneoeps * A_1ovsp + A_1ovsp * gamma + B_1ovsp
      Itot1ovsp_dimreg=Itot1ovsp_dimreg*Itot_nonlog(p1,p3,s)
     .     + A_1ovsp * Itot_log(p1,p3,s)

*
      return
      end
*
* Eq. 290 Francesco'1 notes (term proportional to s'/s')
*
      real*8 function Itotspovsp_dimreg(p1,p3,s)
      implicit none
      real*8 p1(0:3),p3(0:3)
      real*8 s
      real*8 Itot_nonlog, Itot_log
      external Itot_nonlog, Itot_log
      real*8 A_spovsp,gamma,B_spovsp
      real*8 betapi2
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 oneoeps,mudimreg2,eulergamma
      common/IRdimreg/oneoeps,mudimreg2,eulergamma
      real*8 tridot
      external tridot
*
*     Eq. 286 of Francesco's notes
      betapi2=1d0-4d0*ampi**2/s
*
      A_spovsp= 1.d0/8.d0
      B_spovsp= 1d0/8d0*(-2d0*log(4d0*pi )+ 2d0*log(betapi2)
     .  +log(s/mudimreg2))
      gamma= log(4.d0*pi) + eulergamma
*m
      Itotspovsp_dimreg=-oneoeps*A_spovsp + A_spovsp*gamma + B_spovsp
      Itotspovsp_dimreg= Itotspovsp_dimreg * Itot_nonlog(p1,p3,s)
     +                 + A_spovsp * Itot_log(p1,p3,s)

      return
      end
*
* Eq. 267 Francesc's notes
*
      real*8 function Itot_nonlog(p1,p3,s)
      implicit none
      real*8 p1(0:3),p3(0:3)
      real*8 s,p10,p1z,p30,p3y,p3z
      real*8 a1,b1,c1,a2,b2,c2
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 Inonlogmoretti
      external Inonlogmoretti
      real*8 dot,tridot
      external dot,tridot
      real*8 t,betae,betapi,cth
*
      t = dot(p1-p3,p1-p3)
      betae=sqrt(1d0-4d0*ame**2/s)
      betapi=sqrt(1d0-4d0*ampi**2/s)
      cth=(t-ame**2-ampi**2+s/2d0)/(s/2d0*betae*betapi)
      p10= sqrt(s)/2d0
      p1z= p10*betae
      p30= p10
      p3y= betapi*p30*sqrt(1d0-cth**2)
      p3z= betapi*p30*cth

      a1= p10**2/(p3y**2+p3z**2)
      b1= (p3y**2 * p10**2 + (p3z*p30 - p10*p1z)**2)
      b1= b1/(p3y**2 * p1z**2)
      c1= pi * p10 * (p3y**2 + p3z*(p3z-p1z))
      c1= c1/p3y**2/p1z**2/sqrt(p3y**2+p3z**2)

      a2= p10**2/p1z**2
      b2= b1
      c2= pi*p10*(p1z-p3z) / p3y**2 / p1z**2

      !!!!!! checked that a_i are correct
*
      Itot_nonlog= 2.d0*c1*Inonlogmoretti(a1,b1)
     +            +2.d0*c2*Inonlogmoretti(a2,b2)



*
      return
      end
*
* Eq. 268 Francesc's notes
*
      real*8 function Itot_log(p1,p3,s)
      implicit none
      real*8 p1(0:3),p3(0:3)
      real*8 s,p10,p1z,p30,p3y,p3z
      real*8 a1,b1,c1,a2,b2,c2
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      complex*16 Ilogmoretti
      external Ilogmoretti
      real*8 T
      real*8 dot
      external dot
      real*8 betae,betapi,cth
*
      t = dot(p1-p3,p1-p3)
      betae=sqrt(1d0-4d0*ame**2/s)
      betapi=sqrt(1d0-4d0*ampi**2/s)
      cth=(t-ame**2-ampi**2+s/2d0)/(s/2d0*betae*betapi)
      p10= sqrt(s)/2d0
      p1z= p10*betae
      p30= p10
      p3y= betapi*p30*sqrt(1d0-cth**2)
      p3z= betapi*p30*cth

      a1= p30**2/(p3y**2+p3z**2)
      b1= (p3y**2 * p10**2 + (p3z*p30 - p10*p1z)**2)
      b1= b1/(p3y**2 * p1z**2)
      c1= pi * p10 * (p3y**2 + p3z*(p3z-p1z))
      c1= c1/p3y**2/p1z**2/sqrt(p3y**2+p3z**2)

      a2= p10**2/p1z**2
      b2= b1
      c2= pi*p10*(p1z-p3z) / p3y**2 / p1z**2
*
      Itot_log= 2.d0*c1*dble(Ilogmoretti(a1,b1))
     +         +2.d0*c2*dble(Ilogmoretti(a2,b2))
*
      return
      end
*
      function boxdm_lt_tot(s,t,u,lambda2)
      implicit none
      real*8 boxdm_lt_tot
      real*8 s,t,u,lambda2
      real*8 C1,C2

      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 mp2,me2

      real*8 I1lt,I2lt
      external I1lt
      external I2lt
      real*8 I1,I2

      real*8 k,lognum,logden,logd0,betap
      real*8 I1t,I1u,I2t,I2u,It,Iu,C1t,C1u,C2t,C2u
      real*8 Ce,cC,cp

      mp2=ampi**2
      me2=ame**2

      betap=sqrt(1d0-4d0*mp2/s)

      Ce=s*(u-t)*pi/(s)*log(ame**2/s)
      Cp=(2d0*mp2-s)*(t-u)*pi/(s*betap)*log((1d0-betap)/(1d0+betap))
      CC=(Ce+Cp)*log(1d0-betap**2)*4d0*s/(4d0*mp2*s-s**2+(t-u)**2)

      C1t=4d0*s*(me2+mp2-t)*(me2*(2d0*mp2-t+u-lambda2)
     . +2d0*mp2**2-2d0*mp2*t+t*(t-u+lambda2))
      C1t=C1t/(4d0*mp2*s-s**2+(t-u)**2)

      C2t=-4d0*s*(me2-t)*(me2+mp2-t)
      C2t=C2t/(4d0*mp2*s-s**2+(t-u)**2)

      I1t=I1lt(s,t,lambda2)
      I2t=I2lt(s,t,lambda2)

      It=I1t*C1t+I2t*C2t

      C1u=4d0*s*(me2+mp2-u)*(me2*(2d0*mp2-u+t-lambda2)
     . +2d0*mp2**2-2d0*mp2*u+u*(u-t+lambda2))
      C1u=C1u/(4d0*mp2*s-s**2+(t-u)**2)

      C2u=-4d0*s*(me2-u)*(me2+mp2-u)
      C2u=C2u/(4d0*mp2*s-s**2+(t-u)**2)

      I1u=I1lt(s,u,lambda2)
      I2u=I2lt(s,u,lambda2)

      Iu=I1u*C1u+I2u*C2u
      boxdm_lt_tot=-Cc+It-Iu
      return
      End

      !!===============================!!
      !!        DISPERSIVE FSR         !!   FPU 24.07.24
      !!===============================!!

      function FSRdisp(s,xi,csi,min,max)
      implicit none
      real*8 s,t,u,s1,s2
      real*8 xi,csi,min,max

      Complex*16 CTdisp
      external CTdisp

      complex*16 vertexpi
      external vertexpi

      Complex*16 FSRdisp

      external ffpi
      real*8 ff2s1,ff2s2,ff2s
      complex*16 ffs1,ffs2,ffs
      real*8 imffs1,reffs1,imffs2,reffs2

      complex*16 VTpoint,VTpoledisp,VTdispdisp

      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi


      s1 = min + xi*(max-min)
      s2 = min + csi*(max-min)

      call ffpi(s1,ffs1,ff2s1)
      imffs1=aimag(ffs1)
      reffs1=dble(ffs1)

      call ffpi(s2,ffs2,ff2s2)
      imffs2=aimag(ffs2)
      reffs2=dble(ffs2)

      !! Here I use the Eq. (206) of Francesco's notes
!      FSRpoint=CTdisp(0d0)+vertexpi(s,0d0)
      VTpoledisp=vertexpi(s,s1)
      VTpoledisp=VTpoledisp*(-2d0/pi)*imffs1
      VTpoledisp=VTpoledisp*(max-min)

      VTdispdisp=vertexpi(s,s2)*s2-vertexpi(s,s1)*s1
      VTdispdisp =1d0/(pi**2)*(VTdispdisp)/(s2-s1)
      VTdispdisp=VTdispdisp*imffs2*imffs1
      VTdispdisp=VTdispdisp*(max-min)**2

      FSRdisp=VTpoledisp+VTdispdisp
      return
      end

      function CTpidisp(s,xi,csi,min,max)
      implicit none
      real*8 s,t,u,s1,s2
      real*8 xi,csi,min,max

      Complex*16 CTdisp
      external CTdisp

      complex*16 vertexpi
      external vertexpi

      Complex*16 CTpidisp

      external ffpi
      real*8 ff2s1,ff2s2,ff2s
      complex*16 ffs1,ffs2,ffs
      real*8 imffs1,reffs1,imffs2,reffs2

      complex*16 CTpoint,CTpoledisp,CTdispdisp

      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi


      s1 = min + xi*(max-min)
      s2 = min + csi*(max-min)

      call ffpi(s1,ffs1,ff2s1)
      imffs1=aimag(ffs1)
      reffs1=dble(ffs1)

      call ffpi(s2,ffs2,ff2s2)
      imffs2=aimag(ffs2)
      reffs2=dble(ffs2)

      !! Here I use the Eq. (206) of Francesco's notes
!      FSRpoint=CTdisp(0d0)+vertexpi(s,0d0)
      CTpoledisp=CTdisp(s1)
      CTpoledisp=CTpoledisp*(-2d0/pi)*imffs1
      CTpoledisp=CTpoledisp*(max-min)

      CTdispdisp=CTdisp(s2)*s2-CTdisp(s1)*s1
      CTdispdisp =1d0/(pi**2)*(CTdispdisp)/(s2-s1)
      CTdispdisp=CTdispdisp*imffs2*imffs1
      CTdispdisp=CTdispdisp*(max-min)**2

      CTpidisp=CTpoledisp+CTdispdisp
      return
      end


      function CTdisp(sp)
#ifdef COLLIER
      use collier
#endif
      implicit none
      character*2 fs
      common/finalstate/fs

      complex*16 ctdisp

      Complex*16 A0p
      complex*16 B0p0p,B0ppsp,DB0ppsp,DB0p0p
      real*8 deltauv,lambda2,mu2dim,m2min
      common/loop/deltauv,lambda2,mu2dim,m2min
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 oneoeps,mudimreg2,eulergamma
      common/IRdimreg/oneoeps,mudimreg2,eulergamma

      real*8 me2,mp2,sp

      double complex c1
      parameter (c1 = (1.d0,0.d0))

#ifdef COLLIER
      call GetDeltaUV_cll(deltauv)
      call GetMuUV2_cll(mu2dim)
      call GetMuIR2_cll(lambda2)
#endif

#ifdef LOOPTOOLS
#include "looptools.h"
      call clearcache
      call ltini
      lambda2=getlambda()
      deltauv=getdelta()
      mu2dim=getmudim()
#endif

      me2=ame*ame
      mp2=ampi*ampi

****************************************************************

#ifdef LOOPTOOLS
      B0p0p=B0(mp2,lambda2,mp2)
      B0ppsp=B0(mp2,mp2,sp)
      DB0ppsp =DB0(mp2,mp2,sp)
      DB0p0p=DB0(mp2,0d0,mp2)
#endif

#ifdef COLLIER
      call InitEvent_cll(1)
      call A0_cll(A0p,c1*mp2)

      call B0_cll(B0p0p,c1*mp2,c1*0d0,c1*mp2)
      call B0_cll(B0ppsp,c1*mp2,c1*mp2,c1*sp)

      call DB0_cll(DB0ppsp,c1*mp2,c1*mp2,c1*sp)
      call DB0_cll(DB0p0p,c1*mp2,c1*0d0,c1*mp2)
#endif
      ctdisp=(-sp+4d0*mp2)*DB0ppsp + 2d0*B0ppsp

      return
      end

      function Vertexpi(s,sp)
#ifdef COLLIER
      use collier
#endif
      implicit none
      character*2 fs
      common/finalstate/fs

      complex*16 Vertexpi
      complex*16 B0p0p,B0spp,C0ppspspp,B0ppsp
      real*8 deltauv,lambda2,mu2dim,m2min
      common/loop/deltauv,lambda2,mu2dim,m2min
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 oneoeps,mudimreg2,eulergamma
      common/IRdimreg/oneoeps,mudimreg2,eulergamma

      real*8 me2,mp2,sp,s

      double complex c1
      parameter (c1 = (1.d0,0.d0))
#ifdef COLLIER
      call GetDeltaUV_cll(deltauv)
      call GetMuUV2_cll(mu2dim)
      call GetMuIR2_cll(lambda2)
#endif

#ifdef LOOPTOOLS
#include "looptools.h"
      call clearcache
      call ltini
      lambda2=getlambda()
      deltauv=getdelta()
      mu2dim=getmudim()
#endif

      me2=ame*ame
      mp2=ampi*ampi

****************************************************************


#ifdef LOOPTOOLS
      B0ppsp=B0(mp2,mp2,sp)
      B0spp=B0(s,mp2,mp2)

      C0ppspspp =C0i(cc0,mp2,mp2,s,mp2,sp,mp2)

#endif

#ifdef COLLIER
      call InitEvent_cll(1)

      call B0_cll(B0ppsp,c1*mp2,c1*mp2,c1*sp)
      call B0_cll(B0spp,c1*s,c1*mp2,c1*mp2)

      call C0_cll(C0ppspspp ,c1*mp2,c1*mp2,c1*s,c1*mp2,c1*sp,c1*mp2)
#endif
      vertexpi = (1D0*(-((s + sp)*B0ppsp*2D0) +
     &      (-(sp*1D0) - s*2D0 + ampi**2*4D0)*
     &       (-(B0spp*2D0) +
     &         C0ppspspp*
     &          (-(s*1D0) - sp*2D0 + ampi**2*4D0))))/
     &  (-(s*1D0) + ampi**2*4D0)


      return
      end

      !!! checks Mauro Moretti lambda calculatio
      function Itot1ovsp_lambda(p1,p2,p3,p4,s,t,lambda2)
      implicit none
      real*8 p1,p2,p3,p4
      real*8 Itot1ovsp_lambda
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      real*8 I2Moretti,s,T,I1moretti
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi

      real*8 p1p3,p1p2,e1,be1,e3,be3,esoft,lambda2
      complex*16 i13

      real*8 mp2,me2,f

      real*8 num,den
      real*8 dot
      external dot


      mp2=ampi**2
      me2=ame**2


      f= sqrt(t**2+(me2-mp2)**2-2d0*t*(me2+mp2))
      num=mp2+me2-t+F
      den=mp2+me2-t-F
      I2moretti=pi/s*log(s/(4d0*mp2))*1d0/f*log(num/den)

      p1p2=dot(p1,p2)
      me2=ame*ame
      mp2=ampi**2
      S=2.d0*me2+2.d0*p1p2

      p1p3=(me2+mp2-t)/2d0

      e1 = dsqrt(s)/2.d0
      e3=e1

      be1=sqrt(1.d0-4.d0*me2/S)
      be3=sqrt(1.d0-4.d0*mp2/S)
      esoft=sqrt(s)/2d0*be3**2

      call intpipj(pi,p1p3,e1,be1,ame,e3,be3,ampi,esoft,lambda2,i13)
      I1moretti=real(i13)/p1p3/(4d0*s)

      Itot1ovsp_lambda=I1moretti+I2moretti

      return
      end


      function Itotspovsp_lambda(p1,p2,p3,p4,s,t,lambda2)
      implicit none
      real*8 p1,p2,p3,p4
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      real*8 Itotspovsp_lambda
      real*8 s,T

      real*8 p1p3,p1p2,e1,be1,e3,be3,esoft,lambda2
      complex*16 i13
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi

      real*8 me2,mp2
      real*8 dot
      external dot

      p1p2=dot(p1,p2)

      me2=ame*ame
      mp2=ampi**2
      S=2.d0*me2+2.d0*p1p2

      p1p3=(me2+mp2-t)/2d0

      e1 = dsqrt(s)/2.d0
      e3=e1

      be1=sqrt(1.d0-4.d0*me2/S)
      be3=sqrt(1.d0-4.d0*mp2/S)
      esoft=sqrt(s)/2d0*be3**2

      call intpipj(pi,p1p3,e1,be1,ame,e3,be3,ampi,esoft,lambda2,i13)
      Itotspovsp_lambda=real(i13)/p1p3/(4d0)

      return
      end


      real*8 function boxdm_moretti_lambda(p1,p2,p3,p4,s,t,u,lambda2)
#ifdef COLLIER
      use collier
#endif
      implicit none
      real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      real*8 s,t,u,lambda2
      real*8 C1ovspt,Cspovspt,C1ovspu,Cspovspu

      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      real*8 mp2,me2

      real*8 Itot1ovspt,Itot1ovspu
      real*8 Itot1ovsp_lambda
      external Itot1ovsp_lambda
      real*8 Itotspovspt,Itotspovspu
      real*8 Itotspovsp_lambda
      external Itotspovsp_lambda

      real*8 k,lognum,logden,logd0,betap
      real*8 I1t,I1u,I2t,I2u,It,Iu,C1t,C1u,C2t,C2u
      real*8 Ce,cC,cp

      real*8 I1lt, I2lt
      external i2lt,i1lt
      integer i
      complex*16 D0ltd
      double complex c1,D0eeppst0e0p
      parameter (c1 = (1.d0,0.d0))
      real*8 Itlt,iult

      mp2=ampi**2
      me2=ame**2

      betap=sqrt(1d0-4d0*mp2/s)

      Ce=s*(u-t)*pi/(s)*log(ame**2/s)
      Cp=(2d0*mp2-s)*(t-u)*pi/(s*betap)*log((1d0-betap)/(1d0+betap))
      CC=(Ce+Cp)*log(1d0-betap**2)*4d0*s/(4d0*mp2*s-s**2+(t-u)**2)

      C1ovspt=4d0*s*(me2+mp2-t)*(me2*(2d0*mp2-t+u)
     . +2d0*mp2**2-2d0*mp2*t+t*(t-u))
      C1ovspt=C1ovspt/(4d0*mp2*s-s**2+(t-u)**2)

      Cspovspt=-4d0*s*(me2-t)*(me2+mp2-t)
      Cspovspt=Cspovspt/(4d0*mp2*s-s**2+(t-u)**2)

      Itot1ovspt=Itot1ovsp_lambda(p1,p2,p3,p4,s,t,lambda2)
      Itotspovspt=Itotspovsp_lambda(p1,p2,p3,p4,s,t,lambda2)


      It=Itot1ovspt*C1ovspt+Itotspovspt*Cspovspt


      C1ovspu=4d0*s*(me2+mp2-u)*(me2*(2d0*mp2-u+t)
     . +2d0*mp2**2-2d0*mp2*u+u*(u-t))
      C1ovspu=C1ovspu/(4d0*mp2*s-s**2+(t-u)**2)

      Cspovspu=-4d0*s*(me2-u)*(me2+mp2-u)
      Cspovspu=Cspovspu/(4d0*mp2*s-s**2+(t-u)**2)

      Itot1ovspu=Itot1ovsp_lambda(p1,p2,p3,p4,s,u,lambda2)
      Itotspovspu=Itotspovsp_lambda(p1,p2,p3,p4,s,u,lambda2)


      Iu=Itot1ovspu*C1ovspu+Itotspovspu*Cspovspu

      boxdm_moretti_lambda=-Cc+It-Iu
      return
      end

      function residue(s,t,u,delta)
      implicit none
      real*8 s,t,u,x2,y2
      real*8 residue,delta

      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi

      complex*16 boxdm
      external boxdm

      real*8 BoxDIR
      external boxdir

      real*8 ff2
      complex*16 ff
      external ffpi

      complex*16 Boxpdt,boxpdu
      real*8 BoxpdtIR,boxpd,BoxpduIR

      real*8 boxpd_t,boxpd_u

      character*10 filename
      integer is

      x2=s+delta
      y2=0d0
      boxpdt=Boxdm(s,t,u,x2,y2)
      boxpdt=boxpdt*(x2-s)
      boxpdu=-Boxdm(s,u,t,x2,y2)
      boxpdu=boxpdu*(x2-s)

      boxpd_t=aimag(boxpdt)
      boxpd_u=aimag(boxpdu)


      call ffpi(x2,ff,ff2)

      boxpd=boxpd_t+boxpd_u
      boxpd=boxpd*aimag(ff)
      residue=boxpd/x2
      return
      end





*-------------------------------------------------------*
*     ROUTINES TO REDUCE TENSOR SCALAR INTEGRALS        *
*     WHERE PAVE ARE EVALUATED IN DOUBLE PRECISION      *
*                C_I,D_I,D_IJ,E_I,E_IJ                  *
*-------------------------------------------------------*



*-----------ROUTINES TO REDUCE CI,CIJ-------------------*      
      subroutine Ci_d(k1,k2,m1q,m2q,m3q,C0m,C1m,C2m)
c     In this routine we evaluate the numerical decomposition of 
C     C_i as a function of C0,B0,A0. The convention for the loop integral is
c      
c         1             1          {q^mu}      
c     --------     ----------    -------------  
c     (q^2-m1q)  ((q+k1)^2-m2q)  ((q+k2)^2-m3q) 
c

#ifdef COLLIER
      use collier
#endif      

      implicit none
      double precision k1(0:3),k2(0:3),k3(0:3)
      double precision p1LT(0:3),p2LT(0:3),p3LT(0:3),p4LT(0:3)
      complex*16 m1q,m2q,m3q
      complex*16 C0m,C1m,C2m
      complex*16 B0AC,B0AB,B0BC
      complex*16 p1LT2,p2LT2,p3LT2
      complex*16 p1C2,p2C2,p3C2
      integer i
      complex*32 re,im
      complex*16 k1k1,k2k2,k1k2
      
      external dotq
      real*16 dotq

      external dot
      real*8 dot
      double complex c1
      parameter (c1 = (1.d0,0.d0))
      
      external closetomassc
      complex*32 closetomassc

      external closetomasscd
      complex*16 closetomasscd
      re=(1.q0,0.q0)

      k1k1=closetomasscd(dot(k1,k1))
      k2k2=closetomasscd(dot(k2,k2))
      k1k2=closetomasscd(dot(k1,k2))

      
      do i= 0,3
         p1LT(i)= k1(i)
         p2LT(i)= k2(i)-k1(i)
         p3LT(i)= k2(i)
      enddo

*     
      p1LT2= closetomasscd(dot(p1LT,p1LT))
      p2LT2= closetomasscd(dot(p2LT,p2LT))
      p3LT2= closetomasscd(dot(p3LT,p3LT))
      
!!!!!!!!!!!!!!!!!!!!LT non sa nulla di k1k2, lo riscrive partendo da P2LT2-m1^2-m2^2, lo faccio anche qui

      k1k2 = -(p2LT2-k1k1-k2k2)/2.d0

      p1C2 = p1LT2*1.D0
      p2C2 = p2LT2*1.D0
      p3C2 = p3LT2*1.D0
!!!!!!!!!!!!!!!!!!!!
      
      
#ifdef COLLIER      
      call C0_cll(C0m,p1C2*c1,p2C2*c1,p3C2*c1,m1q*c1,m2q*c1,m3q*c1)

      call B0_cll(B0AB,p1C2*c1,m1q*c1,m2q*c1)

      call B0_cll(B0AC,p3C2*c1,m1q*c1,m3q*c1)

      call B0_cll(B0BC,p2C2*c1,m2q*c1,m3q*c1)
#endif


      C1m = (-(B0AB*k1k2) + B0BC*(k1k2 - k2k2) + B0AC*k2k2 +
     >     C0m*(-(k1k1*k2k2) + k2k2*(-m1q + m2q) + k1k2*(k2k2 + m1q -
     >     m3q)))/(-2.q0*k1k2**2 + 2.q0*k1k1*k2k2)
      C2m =(-(B0AB*k1k1) + B0BC*(k1k1 - k1k2) + B0AC*k1k2 +
     >     C0m*k1k2*(-m1q + m2q) + C0m*k1k1*(-k1k2 + k2k2 + m1q - m3q))
     >     /(2.q0*(k1k2**2 - k1k1*k2k2))

      return
      end


*-----------ROUTINES TO REDUCE DI,DIJ-------------------*      
      subroutine Di_d(k1,k2,k3,m1q,m2q,m3q,m4q,D0m,Dim)
c     In this routine we evaluate the numerical decomposition of 
C     D_i as a function of D0,C0,B0,A0. The convention for the loop integral is
c      
c         1             1          {q^mu}       1
c     --------     ----------    -------------   ------------- 
c     (q^2-m1q)  ((q+k1)^2-m2q)  ((q+k2)^2-m3q)  ((q+k3)^2-m4q) 
c

#ifdef COLLIER
      use collier
#endif 
      implicit none
      double precision k1(0:3),k2(0:3),k3(0:3)
      double precision p1LT(0:3),p2LT(0:3),p3LT(0:3),p4LT(0:3)
      complex*16 m1q,m2q,m3q,m4q
      complex*16 C0ABC,C0ABD,C0ACD,C0BCD
      complex*16 D0m
      complex*16 Dim(3)
      complex*16 p1LT2,p2LT2,p3LT2,p4LT2
      complex*16 x23LT,x12LT

      complex*16 p1C2,p2C2,p3C2,p4C2
      complex*16 x23C,x12C
      integer i
      complex*32 re,im
      complex*16 k1k1,k1k2,k1k3
      complex*16 k2k2,k2k3
      complex*16 k3k3

      double complex c1
      parameter (c1 = (1.d0,0.d0))

      external dotq
      real*16 dotq

      external dot
      real*8 dot
      
      external closetomassc
      complex*32 closetomassc

      external closetomasscd
      complex*16 closetomasscd
      
      re=(1.q0,0.q0)

      k1k1=closetomasscd(dot(k1,k1))
      k2k2=closetomasscd(dot(k2,k2))
      k3k3=closetomasscd(dot(k3,k3))

      do i= 0,3
         p1LT(i)= k1(i)
         p2LT(i)= k2(i)-k1(i)
         p3LT(i)= k3(i)-k2(i)
         p4LT(i)= -k3(i)
      enddo

*     
      p1LT2= closetomasscd(dot(p1LT,p1LT))
      p2LT2= closetomasscd(dot(p2LT,p2LT))
      p3LT2= closetomasscd(dot(p3LT,p3LT))
      p4LT2= closetomasscd(dot(p4LT,p4LT))
      
      x12LT= closetomasscd(dot(k2,k2))
      x23LT= closetomasscd(dot(k3-k1,k3-k1))
      
!!!!!!!!!!!!!!!!!!!!LT non sa nulla di k1k2, lo riscrive partendo da P2LT2-m1^2-m2^2, lo faccio anche qui

      k1k2 = -(p2LT2-k1k1-k2k2)/2.d0
      k1k3 = -(x23LT-k1k1-k3k3)/2.d0
      k2k3 = -(p3LT2-k3k3-k2k2)/2.d0
!!!!!!!!!!!!!!!!!!!!

      p1C2 = p1LT2*1.d0
      p2C2 = p2LT2*1.d0
      p3C2 = p3LT2*1.d0
      p4C2 = p4LT2*1.d0

      x12C = x12LT*1.d0
      x23C = x23LT*1.d0
#ifdef COLLIER
      
      call D0_cll(D0m,p1C2*c1,p2C2*c1,p3C2*c1,p4C2*c1,x12C*c1,
     >     x23C*c1,m1q*c1,m2q*c1,m3q*c1,m4q*c1)

      call C0_cll(C0ABC,p1C2,p2C2,x12C,m1q,m2q,m3q)

      call C0_cll(C0ACD,x12C,p3C2,p4C2,m1q,m3q,m4q)

      call C0_cll(C0ABD,p1C2,x23C,p4C2,m1q,m2q,m4q)

      call C0_cll(C0BCD,p2C2,p3C2,x23C,m2q,m3q,m4q)

#endif
      include 'pipirad/decomposition5pt/Di.F'
      return
      end
      
      subroutine Dij_d(k1,k2,k3,m1q,m2q,m3q,m4q,D0m,Dim,Dijm)

c     In this routine we evaluate the numerical decomposition of D_ij
C     D_i as a function of D0,C0,B0,A0. The convention for the loop integral is
c      
c         1             1      {1,q^mu,q^muq^nu}       1
c     --------     ----------    -------------   ------------- 
c     (q^2-m1q)  ((q+k1)^2-m2q)  ((q+k2)^2-m3q)  ((q+k3)^2-m4q) 
c
#ifdef COLLIER
      use collier
#endif 
      implicit none
      double precision k1(0:3),k2(0:3),k3(0:3)
      double precision p1LT(0:3),p2LT(0:3),p3LT(0:3),p4LT(0:3)
      complex*16 m1q,m2q,m3q,m4q
      complex*16 C0ABC,C0ABD,C0ACD,C0BCD
      complex*16 C1ABC,C2ABC
      complex*16 C1ABD,C2ABD
      complex*16 C1ACD,C2ACD
      complex*16 C1BCD,C2BCD
      
      complex*16 D0m
      complex*16 D1m,D2m,D3m
      complex*16 Dim(3),Dijm(7)
      complex*16 bc0,bc1,bc2,bc3,bc4,bc5,bc6
      complex*16 p1LT2,p2LT2,p3LT2,p4LT2
      complex*16 x23LT,X12LT
      integer i
      complex*32 re,im
      complex*16 k1k1,k1k2,k1k3
      complex*16 k2k2,k2k3
      complex*16 k3k3
      
      external dotq
      real*16 dotq

      external dot
      real*8 dot
      
      external closetomassc
      complex*32 closetomassc

      external closetomasscd
      complex*16 closetomasscd
      
      re=(1.q0,0.q0)

      k1k1=closetomasscd(dot(k1,k1))
      k2k2=closetomasscd(dot(k2,k2))
      k3k3=closetomasscd(dot(k3,k3))

      do i= 0,3
         p1LT(i)= k1(i)
         p2LT(i)= k2(i)-k1(i)
         p3LT(i)= k3(i)-k2(i)
         p4LT(i)= -k3(i)
      enddo

*     
      p1LT2= closetomasscd(dot(p1LT,p1LT))
      p2LT2= closetomasscd(dot(p2LT,p2LT))
      p3LT2= closetomasscd(dot(p3LT,p3LT))
      p4LT2= closetomasscd(dot(p4LT,p4LT))
      
      x12LT= closetomasscd(dot(k2,k2))
      x23LT= closetomasscd(dot(k3-k1,k3-k1))
      
!!!!!!!!!!!!!!!!!!!!LT non sa nulla di k1k2, lo riscrive partendo da P2LT2-m1^2-m2^2, lo faccio anche qui

      k1k2 = -(p2LT2-k1k1-k2k2)/2.d0
      k1k3 = -(x23LT-k1k1-k3k3)/2.d0
      k2k3 = -(p3LT2-k3k3-k2k2)/2.d0
!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!Vettoriali e scalare  4 pt      
      call Di_d(k1,k2,k3,m1q,m2q,m3q,m4q,D0m,Dim)
      D1m=Dim(1)
      D2m=Dim(2)
      D3m=Dim(3)
!!!!!!!!!!!!!!!!!!!!!Vettoriali 3 punti
      call Ci_d(k1,k2,m1q,m2q,m3q,C0ABC,C1ABC,C2ABC)
      call Ci_d(k2,k3,m1q,m3q,m4q,C0ACD,C1ACD,C2ACD)
      call Ci_d(k1,k3,m1q,m2q,m4q,C0ABD,C1ABD,C2ABD)
      call Ci_d(k2-k1,k3-k1,m2q,m3q,m4q,C0BCD,C1BCD,C2BCD)

      include 'pipirad/decomposition5pt/Dij.F'
      return
      end
      
*-----------ROUTINES TO REDUCE EI,EIJ-------------------*      

      subroutine Eij_d(k1,k2,k3,k4,m1q,m2q,m3q,m4q,m5q,E0m,Eim,Eijm)
c     In this routine we evaluate the numerical decomposition of E_ij
C     E_i and E_0 as a function of D0,C0,B0,A0. We use the package
c     hexagon.m of https://arxiv.org/pdf/0812.2134 to get the symbolic
c     decomposition. The convention for the loop integral is
c      
c          1             1      {1,q^mu,q^muq^nu}       1              1
c     --------     ----------    -------------   -------------  --------------
c     (q^2-m1q)  ((q+k1)^2-m2q)  ((q+k2)^2-m3q)  ((q+k3)^2-m4q)  ((q+k4)^2-m5q) 
c
#ifdef COLLIER
      use collier
#endif
!      implicit none
      implicit complex*32 (u-w)
      double precision k1(0:3),k2(0:3),k3(0:3),k4(0:3)
      double precision p1LT(0:3),p2LT(0:3),p3LT(0:3),p4LT(0:3),p5LT(0:3)
      complex*16 m1q,m2q,m3q,m4q,m5q
      complex*16 E0m,E1m,E2m,E3m,E4m
      complex*16 E11,E12,E13,E14
      complex*16 E22,E23,E24
      complex*16 E33,E34
      complex*16 E44
      complex*32 den,den1,den2,den3,den4,den5
      complex*16 Eim(4),Eijm(11)
      complex*16 D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q
      complex*16 D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q
      complex*16 D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q
      complex*16 D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q
      complex*16 D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q

      complex*16 C0ap1LT2p2LT2x12LTm1qm2qm3q
      complex*16 C0ap1LT2x23LTx45LTm1qm2qm4q
      complex*16 C0ap1LT2x51LTp5LT2m1qm2qm5q
      complex*16 C0ax12LTp3LT2x45LTm1qm3qm4q
      complex*16 C0ax12LTx34LTp5LT2m1qm3qm5q
      complex*16 C0ax23LTp2LT2p3LT2m4qm2qm3q
      complex*16 C0ax34LTp3LT2p4LT2m5qm3qm4q
      complex*16 C0ax45LTp4LT2p5LT2m1qm4qm5q
      complex*16 C0ax51LTp2LT2x34LTm5qm2qm3q
      complex*16 C0ax51LTx23LTp4LT2m5qm2qm4q      

      complex*32 p1LT2,p2LT2,p3LT2,p4LT2,p5LT2
      complex*32 x12LT,x23LT,x34LT,x45LT,x51LT
      complex*16 p1C2,p2C2,p3C2,p4C2,p5C2
      complex*16 x12C,x23C,x34C,x45C,x51C
      integer i
      complex*32 re,im
      complex*32 k1k1,k1k2,k1k3,k1k4
      complex*32 k2k2,k2k3,k2k4
      complex*32 k3k3,k3k4
      complex*32 k4k4
      
      external dotq
      real*16 dotq

      external closetomassc
      complex*32 closetomassc

      double complex c1
      parameter (c1 = (1.d0,0.d0))

      complex*16 enpi(0:0,0:0,0:0,0:0,0:0)
      complex*16 enpiuv(0:0,0:0,0:0,0:0,0:0)
      
      re=(1.q0,0.q0)

      k1k1=closetomassc(dotq(k1,k1)*re)
      k2k2=closetomassc(dotq(k2,k2)*re)
      k3k3=closetomassc(dotq(k3,k3)*re)
      k4k4=closetomassc(dotq(k4,k4)*re)

      do i= 0,3
         p1LT(i)= k1(i)
         p2LT(i)= k2(i)-k1(i)
         p3LT(i)= k3(i)-k2(i)
         p4LT(i)= k4(i)-k3(i)
         p5LT(i)= -k4(i)
      enddo

*     

      p1LT2= closetomassc(dotq(p1LT,p1LT)*re)
      p2LT2= closetomassc(dotq(p2LT,p2LT)*re)
      p3LT2= closetomassc(dotq(p3LT,p3LT)*re)
      p4LT2= closetomassc(dotq(p4LT,p4LT)*re)
      p5LT2= closetomassc(dotq(p5LT,p5LT)*re)
      
      x12LT= closetomassc(dotq(k2,k2)*re)
      x23LT= closetomassc(dotq(k3-k1,k3-k1)*re)
      x34LT= closetomassc(dotq(k4-k2,k4-k2)*re)
      x45LT= closetomassc(dotq(k3,k3)*re)
      x51LT= closetomassc(dotq(k1-k4,k1-k4)*re)

      p1C2 = p1LT2*1.d0
      p2C2 = p2LT2*1.d0
      p3C2 = p3LT2*1.d0
      p4C2 = p4LT2*1.d0
      p5C2 = p5LT2*1.d0

      
      x12C = x12LT*1.d0
      x23C = x23LT*1.d0
      x34C = x34LT*1.d0
      x45C = x45LT*1.d0
      x51C = x51LT*1.d0
!!!!!!!!!!!!!!!!!!!!LT non sa nulla di kikj
      
      k1k2 = -(p2LT2-k1k1-k2k2)/2.q0
      k1k3 = -(x23LT-k1k1-k3k3)/2.q0
      k1k4 = -(x51LT-k1k1-k4k4)/2.q0
      k2k3 = -(p3LT2-k3k3-k2k2)/2.q0
      k2k4 = -(x34LT-k2k2-k4k4)/2.q0
      k3k4 = -(p4LT2-k3k3-k4k4)/2.q0

      include 'pipirad/decomposition5pt/assignements.f'


#ifdef COLLIER

      call D0_cll(D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q,
     >     c1*x51C,c1*p2C2,c1*p3C2,c1*p4C2,c1*x34C,c1*x23C,
     >     c1*m5q,c1*m2q,c1*m3q,c1*m4q)
      call D0_cll(D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q
     >     ,c1*x12C,c1*p3C2,c1*p4C2,c1*p5C2,c1*x45C,c1*x34C,
     >     c1*m1q,c1*m3q,c1*m4q,c1*m5q)
      call D0_cll(D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q
     >     ,c1*p1C2,c1*x23C,c1*p4C2,c1*p5C2,c1*x45C,c1*x51C,
     >     c1*m1q,c1*m2q,c1*m4q,c1*m5q)
      call D0_cll(D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q,
     >     c1*p1C2,c1*p2C2,c1*x34C,c1*p5C2,c1*x12C,c1*x51C,
     >     c1*m1q,c1*m2q,c1*m3q,c1*m5q)
      call D0_cll(D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q,
     >     c1*p1C2,c1*p2C2,c1*p3C2,c1*x45C,c1*x12C,c1*x23C,
     >     c1*m1q,c1*m2q,c1*m3q,c1*m4q)

      call C0_cll(C0ap1LT2p2LT2x12LTm1qm2qm3q,
     >     c1*p1C2,c1*p2C2,c1*x12C,c1*m1q,c1*m2q,c1*m3q)
      call C0_cll(C0ap1LT2x23LTx45LTm1qm2qm4q,
     >     c1*p1C2,c1*x23C,c1*x45C,c1*m1q,c1*m2q,c1*m4q)
      call C0_cll(C0ap1LT2x51LTp5LT2m1qm2qm5q,
     >     c1*p1C2,c1*x51C,c1*p5C2,c1*m1q,c1*m2q,c1*m5q)
      call C0_cll(C0ax12LTp3LT2x45LTm1qm3qm4q,
     >     c1*x12C,c1*p3C2,c1*x45C,c1*m1q,c1*m3q,c1*m4q)
      call C0_cll(C0ax12LTx34LTp5LT2m1qm3qm5q,
     >     c1*x12C,c1*x34C,c1*p5C2,c1*m1q,c1*m3q,c1*m5q)
      call C0_cll(C0ax23LTp2LT2p3LT2m4qm2qm3q,
     >     c1*x23C,c1*p2C2,c1*p3C2,c1*m4q,c1*m2q,c1*m3q)
      call C0_cll(C0ax34LTp3LT2p4LT2m5qm3qm4q,
     >     c1*x34C,c1*p3C2,c1*p4C2,c1*m5q,c1*m3q,c1*m4q)
      call C0_cll(C0ax45LTp4LT2p5LT2m1qm4qm5q,
     >     c1*x45C,c1*p4C2,c1*p5C2,c1*m1q,c1*m4q,c1*m5q)
      call C0_cll(C0ax51LTp2LT2x34LTm5qm2qm3q,
     >     c1*x51C,c1*p2C2,c1*x34C,c1*m5q,c1*m2q,c1*m3q)
      call C0_cll(C0ax51LTx23LTp4LT2m5qm2qm4q,
     >     c1*x51C,c1*x23C,c1*p4C2,c1*m5q,c1*m2q,c1*m4q)
#endif
!!!!!!!E0!!!!!!!!

      include 'pipirad/decomposition5pt/E0.F'

!!!!!!!Ei!!!!!!!!

      include 'pipirad/decomposition5pt/Ei.F'
      
!!!!!!Eij!!!!!!!!
      
      include 'pipirad/decomposition5pt/Eij.F'
      
      return
      end


*-------------------------------------------------------*
*     D_i IS GIVEN BY A PROPER COMBINATION OF D0,C0     *
*               D_i = a D0+ b s'D0 +c C0                * 
*     THIS ROUTINES RETURNS:                            *
*     -DISP=1                                           *
*              Standard loop decomposition,as D_i rout. *
*     -DISP=2                                           *
*              D_i=a \int ds'/s' ImD0+b \int ds' ImD0   *
*              where it is understood that the Im.      *
*              is the one arising from the cut of the   *
*              2->2 process (namely behaving as 1/s-s') *
*                                                       *      
*     -DISP=3                                           *
*         D_i= a \tilde{ImD0} + bs'\tilde{ImD0} +cIm C0 *
*         where \tilde{ImD0}=ImD0 -ImD0_(22) given      *
*         ImD0_(22) as the previous point               *
*-------------------------------------------------------*

      subroutine Di_disp(k1,k2,k3,m1q,m2q,m3q,m4q,s,disp,D0m,Dim)
c     pi, pf are the momenta that becomes input of the D0 disp. integral 
c     s=|p1+p2|^2 if FS real photon, s=|p3+p4|^2 if IS real photon
#ifdef COLLIER
      use collier
#endif 
      implicit none
      double precision k1(0:3),k2(0:3),k3(0:3)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      double precision pi(0:3),pf(0:3)
      double precision p1LT(0:3),p2LT(0:3),p3LT(0:3),p4LT(0:3)
      complex*16 m1q,m2q,m3q,m4q
      common/hparam/ampi
      double precision ampi,s,sp	
      double precision masses(4),massesbck(4)
      complex*16 C0ABC,C0ABD,C0ACD,C0BCD
      complex*16 D0m
      complex*16 Dim(3),Dim1(3),Dim2(3)
      complex*32 D0mq,Dimq(3)
      complex*16 p1LT2,p2LT2,p3LT2,p4LT2
      complex*16 x23LT,x12LT

      complex*16 p1C2,p2C2,p3C2,p4C2
      complex*16 x23C,x12C
      integer i,k,ifound,disp
      complex*32 re,im
      complex*16 k1k1,k1k2,k1k3
      complex*16 k2k2,k2k3
      complex*16 k3k3

      double complex c1
      parameter (c1 = (1.d0,0.d0))

      external dotq
      real*16 dotq

      external dot
      real*8 dot
      external closetomassc
      complex*32 closetomassc

      
      external closetomasscd
      complex*16 closetomasscd
      complex*16 red
*------------FsQED
      external Itot1ovsp_dimreg
      real*8 Itot1ovsp_dimreg

      external Itotspovsp_dimreg
      real*8 Itotspovsp_dimreg

      
      external Itotspovsp_lambdarad
      real*8 Itotspovsp_lambdarad

      external Itot1ovsp_lambdarad
      real*8 Itot1ovsp_lambdarad

      real*8 deltauv,lambda2,mu2dim,m2min
      common/loop/deltauv,lambda2,mu2dim,m2min

      common/momentadisp/momenta
      double precision momenta(4,0:3)

      external ImD0
      double complex ImD0

      logical cond
      double precision timein,timefin
      p1=momenta(1,:)
      p2=momenta(2,:)
      p3=momenta(3,:)
      p4=momenta(4,:)

      re=(1.q0,0.q0)
      im=(0.q0,1.q0)
      red = (1.d0,0.d0)
      massesbck(1)  = real(m1q)
      massesbck(2)  = real(m2q)
      massesbck(3)  = real(m3q)
      massesbck(4)  = real(m4q)

      k1k1=closetomasscd(dot(k1,k1))
      k2k2=closetomasscd(dot(k2,k2))
      k3k3=closetomasscd(dot(k3,k3))
      
      do i= 0,3
         p1LT(i)= k1(i)
         p2LT(i)= k2(i)-k1(i)
         p3LT(i)= k3(i)-k2(i)
         p4LT(i)= -k3(i)
      enddo

*     
      
      p1LT2= closetomasscd(dot(p1LT,p1LT))
      p2LT2= closetomasscd(dot(p2LT,p2LT))
      p3LT2= closetomasscd(dot(p3LT,p3LT))
      p4LT2= closetomasscd(dot(p4LT,p4LT))
      
      x12LT= closetomasscd(dot(k2,k2))
      x23LT= closetomasscd(dot(k3-k1,k3-k1))
!!!!!!!!!!!!!!!!!!!!LT non sa nulla di k1k2, lo riscrive partendo da P2LT2-m1^2-m2^2, lo faccio anche qui

      k1k2 = -(p2LT2-k1k1-k2k2)/2.d0
      k1k3 = -(x23LT-k1k1-k3k3)/2.d0
      k2k3 = -(p3LT2-k3k3-k2k2)/2.d0
!!!!!!!!!!!!!!!!!!!!

      p1C2 = p1LT2*1.d0
      p2C2 = p2LT2*1.d0
      p3C2 = p3LT2*1.d0
      p4C2 = p4LT2*1.d0

      x12C = x12LT*1.d0
      x23C = x23LT*1.d0

      if(disp.eq.2) then

         D0m   = (1.d0,0.d0)
         C0ABC = (0.d0,0.d0)
         C0ABD = (0.d0,0.d0)
         C0BCD = (0.d0,0.d0)
         C0ACD = (0.d0,0.d0)
!     Se la D0 non è IR, qui restituisci 0
         call D0IRyon(k1,k2,k3,m1q,m2q,m3q,m4q,cond,pi,pf)
         if(.not.cond) then
            D0m=(0.d0,0.d0)
            Dim(1)=(0.d0,0.d0)
            Dim(2)=(0.d0,0.d0)
            Dim(3)=(0.d0,0.d0)
            return
         endif
*------------------------------------------
         
         include 'pipirad/decomposition5pt/Di.F'
*------------QUI L'UNICA M^2 CHE PUÒ ESSERE > 4M_PI^2 È QUELLA DISPERSIVA
         Dim1= Dim
         
         do i=1,4
            if (massesbck(i).ge.4.d0*ampi**2) then
               ifound = i
               if (i.eq.1) m1q = 0.d0
               if (i.eq.2) m2q = 0.d0
               if (i.eq.3) m3q = 0.d0
               if (i.eq.4) m4q = 0.d0
            endif
         enddo
         D0m   = (1.d0,0.d0)
         include 'pipirad/decomposition5pt/Di.F'
         Dim2= Dim
         
!--------QUI ORA RESTITUISCO L'INTEGRALE DISPERSIVO TENENDO CONTO CHE Dim2 \propto s'  
         Dim = ((Dim1-Dim2)*Itotspovsp_lambdarad(pi,pf,s,lambda2)
     >        /massesbck(ifound)   +
     >        Dim2*Itot1ovsp_lambdarad(pi,pf,s,lambda2))*im
!--------LA D0 NEI BOX (CHE È DOVE CHIAMO QUESTA ROUTINE) NON È MAI MOLTIPLICATA PER S'
         D0m = Itot1ovsp_lambdarad(pi,pf,s,lambda2)*im
         m1q  =  massesbck(1)*(1.d0,0.d0)
         m2q  =  massesbck(2)*(1.d0,0.d0)
         m3q  =  massesbck(3)*(1.d0,0.d0)
         m4q  =  massesbck(4)*(1.d0,0.d0)
         return
      elseif((disp.eq.1).or.(disp.eq.3)
     >        .or.(disp.eq.4).or.(disp.eq.6)) then
#ifdef COLLIER
         
         call D0_cll(D0m,p1C2*c1,p2C2*c1,p3C2*c1,p4C2*c1,x12C*c1,
     >        x23C*c1,m1q*c1,m2q*c1,m3q*c1,m4q*c1)

         call C0_cll(C0ABC,p1C2,p2C2,x12C,m1q,m2q,m3q)

         call C0_cll(C0ACD,x12C,p3C2,p4C2,m1q,m3q,m4q)

         call C0_cll(C0ABD,p1C2,x23C,p4C2,m1q,m2q,m4q)

         call C0_cll(C0BCD,p2C2,p3C2,x23C,m2q,m3q,m4q)

#endif
         if((disp.eq.3).or.(disp.eq.6)) then
            do i=1,4
               if (massesbck(i).ge.4.d0*ampi**2) then
                  ifound = i
                  sp = massesbck(i)
               endif
            enddo

            D0m = (D0m - CONJG(D0m))/2.d0
!!!   Here i subtract from the D0 (checking if it is IR)
!!!   the imaginary part the Im coming from the cut
!!!   that originates the behaviour 1/(eta-s') which is integrated analitically
!!!   I do this because in a 2->3 process there are other (finite) cuts contributing to Im

            call D0IRyon(k1,k2,k3,m1q,m2q,m3q,m4q,cond,pi,pf)

            if(cond) then
               D0m = D0m - ImD0(pi,pf,s,sp,lambda2)
            endif
            
            C0ABC = (C0ABC - CONJG(C0ABC))/2.d0
            C0ACD = (C0ACD - CONJG(C0ACD))/2.d0
            C0ABD = (C0ABD - CONJG(C0ABD))/2.d0
            C0BCD = (C0BCD - CONJG(C0BCD))/2.d0

c     Se disp=6 voglio le Pave in s' ma voglio i coefficienti in \eta
            if(disp.eq.6) then
               do i=1,4
                  if (massesbck(i).ge.4.d0*ampi**2) then
                     if (i.eq.1) m1q = s
                     if (i.eq.2) m2q = s
                     if (i.eq.3) m3q = s
                     if (i.eq.4) m4q = s
                  endif
               enddo
            endif
         endif

c     Se disp=4 voglio le Pave in s' ma voglio i coefficienti in \eta
         if(disp.eq.4) then
            do i=1,4
               if (massesbck(i).ge.4.d0*ampi**2) then
                  if (i.eq.1) m1q = s
                  if (i.eq.2) m2q = s
                  if (i.eq.3) m3q = s
                  if (i.eq.4) m4q = s
               endif
            enddo
         endif

         include 'pipirad/decomposition5pt/Di.F'

c$$$         if(disp.eq.1) then
c$$$            do i=1,4
c$$$               if (massesbck(i).ge.4.d0*ampi**2) then
c$$$                  if(abs(massesbck(i)/s-1.d0).le.1d-3) then
c$$$                     call Di(k1,k2,k3,m1q*re,m2q*re,
c$$$     >                    m3q*re,m4q*re,D0mq,Dimq)
c$$$                     D0m=D0mq
c$$$                     Dim=Dimq
c$$$                  endif
c$$$               endif
c$$$            enddo
c$$$         endif

c     go back to original masses value
         m1q  =  massesbck(1)*(1.d0,0.d0)
         m2q  =  massesbck(2)*(1.d0,0.d0)
         m3q  =  massesbck(3)*(1.d0,0.d0)
         m4q  =  massesbck(4)*(1.d0,0.d0)
         
         return
      elseif(disp.eq.5) then
         D0m   = (1.d0,0.d0)
         C0ABC = (0.d0,0.d0)
         C0ABD = (0.d0,0.d0)
         C0BCD = (0.d0,0.d0)
         C0ACD = (0.d0,0.d0)
!     Se la D0 non è IR, qui restituisci 0
         call D0IRyon(k1,k2,k3,m1q,m2q,m3q,m4q,cond,pi,pf)
         if(.not.cond) then
            D0m=(0.d0,0.d0)
            Dim(1)=(0.d0,0.d0)
            Dim(2)=(0.d0,0.d0)
            Dim(3)=(0.d0,0.d0)
            return
         endif

         if(cond) then
            D0m=Itot1ovsp_lambdarad(pi,pf,s,lambda2)*im
         endif
*------------------------------------------
         do i=1,4
            if (massesbck(i).ge.4.d0*ampi**2) then
               ifound = i
               if (i.eq.1) m1q = s
               if (i.eq.2) m2q = s
               if (i.eq.3) m3q = s
               if (i.eq.4) m4q = s
            endif
         enddo

         include 'pipirad/decomposition5pt/Di.F'

         m1q  =  massesbck(1)*(1.d0,0.d0)
         m2q  =  massesbck(2)*(1.d0,0.d0)
         m3q  =  massesbck(3)*(1.d0,0.d0)
         m4q  =  massesbck(4)*(1.d0,0.d0)
         return
      endif

      end

      subroutine Eij_disp(k1,k2,k3,k4,m1q,m2q,m3q,m4q,m5q,
     >     s,disp,E0m,Eim,Eijm)
#ifdef COLLIER
      use collier
#endif      
!     implicit none
      implicit complex*16 (u-w) !cambiato
      double precision k1(0:3),k2(0:3),k3(0:3),k4(0:3)
      double precision p1LT(0:3),p2LT(0:3),p3LT(0:3),p4LT(0:3),p5LT(0:3)
      double precision pi(0:3),pf(0:3)
      complex*16 m1q,m2q,m3q,m4q,m5q
      complex*16 mbck
      double precision masses(5)
      common/hparam/ampi
      double precision ampi,s,sp	
      complex*16 E0m,E1m,E2m,E3m,E4m
      complex*16 E11,E12,E13,E14
      complex*16 E22,E23,E24
      complex*16 E33,E34
      complex*16 E44
      complex*16 den,den1,den2,den3,den4,den5 !! cambiato
      complex*16 Eim(4),Eijm(11)
      complex*32 E0mq,Eimq(4),Eijmq(11)
      complex*16 D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q
      complex*16 D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q
      complex*16 D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q
      complex*16 D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q
      complex*16 D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q

      complex*16 C0ap1LT2p2LT2x12LTm1qm2qm3q
      complex*16 C0ap1LT2x23LTx45LTm1qm2qm4q
      complex*16 C0ap1LT2x51LTp5LT2m1qm2qm5q
      complex*16 C0ax12LTp3LT2x45LTm1qm3qm4q
      complex*16 C0ax12LTx34LTp5LT2m1qm3qm5q
      complex*16 C0ax23LTp2LT2p3LT2m4qm2qm3q
      complex*16 C0ax34LTp3LT2p4LT2m5qm3qm4q
      complex*16 C0ax45LTp4LT2p5LT2m1qm4qm5q
      complex*16 C0ax51LTp2LT2x34LTm5qm2qm3q
      complex*16 C0ax51LTx23LTp4LT2m5qm2qm4q      

      complex*16 p1LT2,p2LT2,p3LT2,p4LT2,p5LT2 !!
      complex*16 x12LT,x23LT,x34LT,x45LT,x51LT !!
      complex*16 p1C2,p2C2,p3C2,p4C2,p5C2 
      complex*16 x12C,x23C,x34C,x45C,x51C
      integer i,ifound
      complex*32 re,im
      complex*16 k1k1,k1k2,k1k3,k1k4 !!
      complex*16 k2k2,k2k3,k2k4 !!
      complex*16 k3k3,k3k4 !!
      complex*16 k4k4 !!
      
      external dotq
      real*16 dotq

      external dot
      real*8 dot
      external closetomassc
      complex*32 closetomassc
      external closetomasscd
      complex*16 closetomasscd
      complex*16 red
      double complex c1
      parameter (c1 = (1.d0,0.d0))

      complex*16 enpi(0:0,0:0,0:0,0:0,0:0)
      complex*16 enpiuv(0:0,0:0,0:0,0:0,0:0)

      
*------------FsQED
      
      common/gvmdchecks/iscoll,fscoll
      integer iscoll,fscoll
      
      external Itot1ovsp_dimreg
      real*8 Itot1ovsp_dimreg

      external Itotspovsp_dimreg
      real*8 Itotspovsp_dimreg

      
      external Itotspovsp_lambdarad
      real*8 Itotspovsp_lambdarad

      external Itot1ovsp_lambdarad
      real*8 Itot1ovsp_lambdarad

      real*8 deltauv,lambda2,mu2dim,m2min
      common/loop/deltauv,lambda2,mu2dim,m2min

      common/momentadisp/momenta
      double precision momenta(4,0:3)

      external ImD0
      double complex ImD0

      logical cond
      integer disp
      double complex setrez

      double precision timein,timefin
      masses(1)  = real(m1q)
      masses(2)  = real(m2q)
      masses(3)  = real(m3q)
      masses(4)  = real(m4q)
      masses(5)  = real(m5q)
      
      D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q = (0.d0, 0.d0)
      D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q = (0.d0, 0.d0)
      D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q = (0.d0, 0.d0)
      D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q = (0.d0, 0.d0)
      D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q = (0.d0, 0.d0)

      C0ap1LT2p2LT2x12LTm1qm2qm3q = (0.d0, 0.d0)
      C0ap1LT2x23LTx45LTm1qm2qm4q = (0.d0, 0.d0)
      C0ap1LT2x51LTp5LT2m1qm2qm5q = (0.d0, 0.d0)
      C0ax12LTp3LT2x45LTm1qm3qm4q = (0.d0, 0.d0)
      C0ax12LTx34LTp5LT2m1qm3qm5q = (0.d0, 0.d0)
      C0ax23LTp2LT2p3LT2m4qm2qm3q = (0.d0, 0.d0)
      C0ax34LTp3LT2p4LT2m5qm3qm4q = (0.d0, 0.d0)
      C0ax45LTp4LT2p5LT2m1qm4qm5q = (0.d0, 0.d0)
      C0ax51LTp2LT2x34LTm5qm2qm3q = (0.d0, 0.d0)
      C0ax51LTx23LTp4LT2m5qm2qm4q = (0.d0, 0.d0)

      
      re=(1.q0,0.q0)
      im=(0.q0,1.q0)
      red = (1.d0,0.d0)

      k1k1=closetomasscd(dot(k1,k1))
      k2k2=closetomasscd(dot(k2,k2))
      k3k3=closetomasscd(dot(k3,k3))
      k4k4=closetomasscd(dot(k4,k4))

      do i= 0,3
         p1LT(i)= k1(i)
         p2LT(i)= k2(i)-k1(i)
         p3LT(i)= k3(i)-k2(i)
         p4LT(i)= k4(i)-k3(i)
         p5LT(i)= -k4(i)
      enddo

*     

      p1LT2= closetomasscd(dot(p1LT,p1LT))
      p2LT2= closetomasscd(dot(p2LT,p2LT))
      p3LT2= closetomasscd(dot(p3LT,p3LT))
      p4LT2= closetomasscd(dot(p4LT,p4LT))
      p5LT2= closetomasscd(dot(p5LT,p5LT))
      
      x12LT= closetomasscd(dot(k2,k2))
      x23LT= closetomasscd(dot(k3-k1,k3-k1))
      x34LT= closetomasscd(dot(k4-k2,k4-k2))
      x45LT= closetomasscd(dot(k3,k3))
      x51LT= closetomasscd(dot(k1-k4,k1-k4))

      p1C2 = p1LT2*1.d0
      p2C2 = p2LT2*1.d0
      p3C2 = p3LT2*1.d0
      p4C2 = p4LT2*1.d0
      p5C2 = p5LT2*1.d0

      
      x12C = x12LT*1.d0
      x23C = x23LT*1.d0
      x34C = x34LT*1.d0
      x45C = x45LT*1.d0
      x51C = x51LT*1.d0
!!!!!!!!!!!!!!!!!!!!LT non sa nulla di kikj
      
      k1k2 = -(p2LT2-k1k1-k2k2)/2.d0
      k1k3 = -(x23LT-k1k1-k3k3)/2.d0
      k1k4 = -(x51LT-k1k1-k4k4)/2.d0
      k2k3 = -(p3LT2-k3k3-k2k2)/2.d0
      k2k4 = -(x34LT-k2k2-k4k4)/2.d0
      k3k4 = -(p4LT2-k3k3-k4k4)/2.d0

      D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q = (0.d0,0.d0)
      D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q = (0.d0,0.d0)
      D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q = (0.d0,0.d0)
      D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q = (0.d0,0.d0)
      D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q = (0.d0,0.d0)

      C0ap1LT2p2LT2x12LTm1qm2qm3q = (0.d0,0.d0)
      C0ap1LT2x23LTx45LTm1qm2qm4q = (0.d0,0.d0)
      C0ap1LT2x51LTp5LT2m1qm2qm5q = (0.d0,0.d0)
      C0ax12LTp3LT2x45LTm1qm3qm4q = (0.d0,0.d0)
      C0ax12LTx34LTp5LT2m1qm3qm5q = (0.d0,0.d0)
      C0ax23LTp2LT2p3LT2m4qm2qm3q = (0.d0,0.d0)
      C0ax34LTp3LT2p4LT2m5qm3qm4q = (0.d0,0.d0)
      C0ax45LTp4LT2p5LT2m1qm4qm5q = (0.d0,0.d0)
      C0ax51LTp2LT2x34LTm5qm2qm3q = (0.d0,0.d0)
      C0ax51LTx23LTp4LT2m5qm2qm4q = (0.d0,0.d0)

!     Classica decomposizione di loop
      if(disp.eq.1) then
#ifdef COLLIER

         call D0_cll(D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q,
     >        c1*x51C,c1*p2C2,c1*p3C2,c1*p4C2,c1*x34C,c1*x23C,
     >        c1*m5q,c1*m2q,c1*m3q,c1*m4q)
         call D0_cll(D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q
     >        ,c1*x12C,c1*p3C2,c1*p4C2,c1*p5C2,c1*x45C,c1*x34C,
     >        c1*m1q,c1*m3q,c1*m4q,c1*m5q)
         call D0_cll(D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q
     >        ,c1*p1C2,c1*x23C,c1*p4C2,c1*p5C2,c1*x45C,c1*x51C,
     >        c1*m1q,c1*m2q,c1*m4q,c1*m5q)
         call D0_cll(D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q,
     >        c1*p1C2,c1*p2C2,c1*x34C,c1*p5C2,c1*x12C,c1*x51C,
     >        c1*m1q,c1*m2q,c1*m3q,c1*m5q)
         call D0_cll(D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q,
     >        c1*p1C2,c1*p2C2,c1*p3C2,c1*x45C,c1*x12C,c1*x23C,
     >        c1*m1q,c1*m2q,c1*m3q,c1*m4q)

         call C0_cll(C0ap1LT2p2LT2x12LTm1qm2qm3q,
     >        c1*p1C2,c1*p2C2,c1*x12C,c1*m1q,c1*m2q,c1*m3q)
         call C0_cll(C0ap1LT2x23LTx45LTm1qm2qm4q,
     >        c1*p1C2,c1*x23C,c1*x45C,c1*m1q,c1*m2q,c1*m4q)
         call C0_cll(C0ap1LT2x51LTp5LT2m1qm2qm5q,
     >        c1*p1C2,c1*x51C,c1*p5C2,c1*m1q,c1*m2q,c1*m5q)
         call C0_cll(C0ax12LTp3LT2x45LTm1qm3qm4q,
     >        c1*x12C,c1*p3C2,c1*x45C,c1*m1q,c1*m3q,c1*m4q)
         call C0_cll(C0ax12LTx34LTp5LT2m1qm3qm5q,
     >        c1*x12C,c1*x34C,c1*p5C2,c1*m1q,c1*m3q,c1*m5q)
         call C0_cll(C0ax23LTp2LT2p3LT2m4qm2qm3q,
     >        c1*x23C,c1*p2C2,c1*p3C2,c1*m4q,c1*m2q,c1*m3q)
         call C0_cll(C0ax34LTp3LT2p4LT2m5qm3qm4q,
     >        c1*x34C,c1*p3C2,c1*p4C2,c1*m5q,c1*m3q,c1*m4q)
         call C0_cll(C0ax45LTp4LT2p5LT2m1qm4qm5q,
     >        c1*x45C,c1*p4C2,c1*p5C2,c1*m1q,c1*m4q,c1*m5q)
         call C0_cll(C0ax51LTp2LT2x34LTm5qm2qm3q,
     >        c1*x51C,c1*p2C2,c1*x34C,c1*m5q,c1*m2q,c1*m3q)
         call C0_cll(C0ax51LTx23LTp4LT2m5qm2qm4q,
     >        c1*x51C,c1*x23C,c1*p4C2,c1*m5q,c1*m2q,c1*m4q)
#endif
         
      elseif(disp.eq.2) then
!     Check quali D0 sono IR e sostituisci D0 con \int ds'/s' Im D0
!     e portano polo 1/(\eta-s')
         call D0IRyon(k1-k4,k2-k4,k3-k4,m5q,m2q,m3q,m4q,cond,pi,pf)
         if(cond) then
            D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q=
     >           Itot1ovsp_lambdarad(pi,pf,s,lambda2)*im
         endif

         call D0IRyon(k2,k3,k4,m1q,m3q,m4q,m5q,cond,pi,pf)
         if(cond) then
            D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q=
     >           Itot1ovsp_lambdarad(pi,pf,s,lambda2)*im
         endif

         call D0IRyon(k1,k3,k4,m1q,m2q,m4q,m5q,cond,pi,pf)
         if(cond) then
            D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q=
     >           Itot1ovsp_lambdarad(pi,pf,s,lambda2)*im
         endif
         
         call D0IRyon(k1,k2,k4,m1q,m2q,m3q,m5q,cond,pi,pf)
         if(cond) then
            D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q=
     >           Itot1ovsp_lambdarad(pi,pf,s,lambda2)*im
         endif

         call D0IRyon(k1,k2,k3,m1q,m2q,m3q,m4q,cond,pi,pf)
         if(cond) then
            D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q=
     >           Itot1ovsp_lambdarad(pi,pf,s,lambda2)*im
         endif

!     Devo restituire il coefficiente d(\eta)*\int ds'/s' Im D0, quindi devo cambiare il valore della massa dispersiva a \eta (qui nella routine s)
         do i=1,5
            if (masses(i).ge.4.d0*ampi**2) then
               ifound = i
               if (i.eq.1) m1q = s
               if (i.eq.2) m2q = s
               if (i.eq.3) m3q = s
               if (i.eq.4) m4q = s
               if (i.eq.5) m5q = s
            endif
         enddo
      elseif(disp.eq.3) then
!     Restituisce la parte Im di tutto ad eccezione delle D0_IR e polari e delle D0_IR (polari) restituisce
!     la parte IM sottratta della parte immaginaria che origina dal taglio che diverge a s'=\eta

         do i=1,5
            if (masses(i).ge.4.d0*ampi**2) then
               ifound = i
               sp = masses(i)
            endif
         enddo

#ifdef COLLIER

         call D0_cll(D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q,
     >        c1*x51C,c1*p2C2,c1*p3C2,c1*p4C2,c1*x34C,c1*x23C,
     >        c1*m5q,c1*m2q,c1*m3q,c1*m4q)
         call D0_cll(D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q
     >        ,c1*x12C,c1*p3C2,c1*p4C2,c1*p5C2,c1*x45C,c1*x34C,
     >        c1*m1q,c1*m3q,c1*m4q,c1*m5q)
         call D0_cll(D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q
     >        ,c1*p1C2,c1*x23C,c1*p4C2,c1*p5C2,c1*x45C,c1*x51C,
     >        c1*m1q,c1*m2q,c1*m4q,c1*m5q)
         call D0_cll(D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q,
     >        c1*p1C2,c1*p2C2,c1*x34C,c1*p5C2,c1*x12C,c1*x51C,
     >        c1*m1q,c1*m2q,c1*m3q,c1*m5q)
         call D0_cll(D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q,
     >        c1*p1C2,c1*p2C2,c1*p3C2,c1*x45C,c1*x12C,c1*x23C,
     >        c1*m1q,c1*m2q,c1*m3q,c1*m4q)

         call C0_cll(C0ap1LT2p2LT2x12LTm1qm2qm3q,
     >        c1*p1C2,c1*p2C2,c1*x12C,c1*m1q,c1*m2q,c1*m3q)
         call C0_cll(C0ap1LT2x23LTx45LTm1qm2qm4q,
     >        c1*p1C2,c1*x23C,c1*x45C,c1*m1q,c1*m2q,c1*m4q)
         call C0_cll(C0ap1LT2x51LTp5LT2m1qm2qm5q,
     >        c1*p1C2,c1*x51C,c1*p5C2,c1*m1q,c1*m2q,c1*m5q)
         call C0_cll(C0ax12LTp3LT2x45LTm1qm3qm4q,
     >        c1*x12C,c1*p3C2,c1*x45C,c1*m1q,c1*m3q,c1*m4q)
         call C0_cll(C0ax12LTx34LTp5LT2m1qm3qm5q,
     >        c1*x12C,c1*x34C,c1*p5C2,c1*m1q,c1*m3q,c1*m5q)
         call C0_cll(C0ax23LTp2LT2p3LT2m4qm2qm3q,
     >        c1*x23C,c1*p2C2,c1*p3C2,c1*m4q,c1*m2q,c1*m3q)
         call C0_cll(C0ax34LTp3LT2p4LT2m5qm3qm4q,
     >        c1*x34C,c1*p3C2,c1*p4C2,c1*m5q,c1*m3q,c1*m4q)
         call C0_cll(C0ax45LTp4LT2p5LT2m1qm4qm5q,
     >        c1*x45C,c1*p4C2,c1*p5C2,c1*m1q,c1*m4q,c1*m5q)
         call C0_cll(C0ax51LTp2LT2x34LTm5qm2qm3q,
     >        c1*x51C,c1*p2C2,c1*x34C,c1*m5q,c1*m2q,c1*m3q)
         call C0_cll(C0ax51LTx23LTp4LT2m5qm2qm4q,
     >        c1*x51C,c1*x23C,c1*p4C2,c1*m5q,c1*m2q,c1*m4q)

#endif
!     Just the imaginary part
         
         D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q = 
     >        setrez(D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q)

         D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q = 
     >        setrez(D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q)

         D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q = 
     >        setrez(D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q)

         D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q = 
     >        setrez(D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q)

         D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q = 
     >        setrez(D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q)

         C0ap1LT2p2LT2x12LTm1qm2qm3q = 
     >        setrez(C0ap1LT2p2LT2x12LTm1qm2qm3q)

         C0ap1LT2x23LTx45LTm1qm2qm4q = 
     >        setrez(C0ap1LT2x23LTx45LTm1qm2qm4q)

         C0ap1LT2x51LTp5LT2m1qm2qm5q = 
     >        setrez(C0ap1LT2x51LTp5LT2m1qm2qm5q)

         C0ax12LTp3LT2x45LTm1qm3qm4q = 
     >        setrez(C0ax12LTp3LT2x45LTm1qm3qm4q)

         C0ax12LTx34LTp5LT2m1qm3qm5q = 
     >        setrez(C0ax12LTx34LTp5LT2m1qm3qm5q)

         C0ax23LTp2LT2p3LT2m4qm2qm3q = 
     >        setrez(C0ax23LTp2LT2p3LT2m4qm2qm3q)

         C0ax34LTp3LT2p4LT2m5qm3qm4q = 
     >        setrez(C0ax34LTp3LT2p4LT2m5qm3qm4q)

         C0ax45LTp4LT2p5LT2m1qm4qm5q = 
     >        setrez(C0ax45LTp4LT2p5LT2m1qm4qm5q)

         C0ax51LTp2LT2x34LTm5qm2qm3q = 
     >        setrez(C0ax51LTp2LT2x34LTm5qm2qm3q)

         C0ax51LTx23LTp4LT2m5qm2qm4q = 
     >        setrez(C0ax51LTx23LTp4LT2m5qm2qm4q)
        
         call D0IRyon(k1-k4,k2-k4,k3-k4,m5q,m2q,m3q,m4q,cond,pi,pf)
         if(cond) then
            D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q=
     >           D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q
     >           -ImD0(pi,pf,s,sp,lambda2)
         endif

         call D0IRyon(k2,k3,k4,m1q,m3q,m4q,m5q,cond,pi,pf)
         if(cond) then
            D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q=
     >           D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q
     >           -ImD0(pi,pf,s,sp,lambda2)
         endif

         call D0IRyon(k1,k3,k4,m1q,m2q,m4q,m5q,cond,pi,pf)
         if(cond) then
            D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q=
     >           D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q
     >           -ImD0(pi,pf,s,sp,lambda2)
         endif
         
         call D0IRyon(k1,k2,k4,m1q,m2q,m3q,m5q,cond,pi,pf)
         if(cond) then
            D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q=
     >           D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q
     >           -ImD0(pi,pf,s,sp,lambda2)
         endif

         call D0IRyon(k1,k2,k3,m1q,m2q,m3q,m4q,cond,pi,pf)
         if(cond) then
            D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q=
     >           D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q
     >           -ImD0(pi,pf,s,sp,lambda2)
         endif


         do i=1,5
            if (masses(i).ge.4.d0*ampi**2) then
               ifound = i
               if (i.eq.1) m1q = s
               if (i.eq.2) m2q = s
               if (i.eq.3) m3q = s
               if (i.eq.4) m4q = s
               if (i.eq.5) m5q = s
            endif
         enddo
      elseif(disp.eq.4) then

#ifdef COLLIER

         call D0_cll(D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q,
     >        c1*x51C,c1*p2C2,c1*p3C2,c1*p4C2,c1*x34C,c1*x23C,
     >        c1*m5q,c1*m2q,c1*m3q,c1*m4q)
         call D0_cll(D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q
     >        ,c1*x12C,c1*p3C2,c1*p4C2,c1*p5C2,c1*x45C,c1*x34C,
     >        c1*m1q,c1*m3q,c1*m4q,c1*m5q)
         call D0_cll(D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q
     >        ,c1*p1C2,c1*x23C,c1*p4C2,c1*p5C2,c1*x45C,c1*x51C,
     >        c1*m1q,c1*m2q,c1*m4q,c1*m5q)
         call D0_cll(D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q,
     >        c1*p1C2,c1*p2C2,c1*x34C,c1*p5C2,c1*x12C,c1*x51C,
     >        c1*m1q,c1*m2q,c1*m3q,c1*m5q)
         call D0_cll(D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q,
     >        c1*p1C2,c1*p2C2,c1*p3C2,c1*x45C,c1*x12C,c1*x23C,
     >        c1*m1q,c1*m2q,c1*m3q,c1*m4q)

         call C0_cll(C0ap1LT2p2LT2x12LTm1qm2qm3q,
     >        c1*p1C2,c1*p2C2,c1*x12C,c1*m1q,c1*m2q,c1*m3q)
         call C0_cll(C0ap1LT2x23LTx45LTm1qm2qm4q,
     >        c1*p1C2,c1*x23C,c1*x45C,c1*m1q,c1*m2q,c1*m4q)
         call C0_cll(C0ap1LT2x51LTp5LT2m1qm2qm5q,
     >        c1*p1C2,c1*x51C,c1*p5C2,c1*m1q,c1*m2q,c1*m5q)
         call C0_cll(C0ax12LTp3LT2x45LTm1qm3qm4q,
     >        c1*x12C,c1*p3C2,c1*x45C,c1*m1q,c1*m3q,c1*m4q)
         call C0_cll(C0ax12LTx34LTp5LT2m1qm3qm5q,
     >        c1*x12C,c1*x34C,c1*p5C2,c1*m1q,c1*m3q,c1*m5q)
         call C0_cll(C0ax23LTp2LT2p3LT2m4qm2qm3q,
     >        c1*x23C,c1*p2C2,c1*p3C2,c1*m4q,c1*m2q,c1*m3q)
         call C0_cll(C0ax34LTp3LT2p4LT2m5qm3qm4q,
     >        c1*x34C,c1*p3C2,c1*p4C2,c1*m5q,c1*m3q,c1*m4q)
         call C0_cll(C0ax45LTp4LT2p5LT2m1qm4qm5q,
     >        c1*x45C,c1*p4C2,c1*p5C2,c1*m1q,c1*m4q,c1*m5q)
         call C0_cll(C0ax51LTp2LT2x34LTm5qm2qm3q,
     >        c1*x51C,c1*p2C2,c1*x34C,c1*m5q,c1*m2q,c1*m3q)
         call C0_cll(C0ax51LTx23LTp4LT2m5qm2qm4q,
     >        c1*x51C,c1*x23C,c1*p4C2,c1*m5q,c1*m2q,c1*m4q)

#endif
         
         do i=1,5
            if (masses(i).ge.4.d0*ampi**2) then
               ifound = i
               if (i.eq.1) m1q = s
               if (i.eq.2) m2q = s
               if (i.eq.3) m3q = s
               if (i.eq.4) m4q = s
               if (i.eq.5) m5q = s
            endif
         enddo

         
      endif


      include 'pipirad/decomposition5pt/assignements.f'

!!!!!!!E0!!!!!!!!

      include 'pipirad/decomposition5pt/E0.F'

!!!!!!!Ei!!!!!!!!

      include 'pipirad/decomposition5pt/Ei.F'
      
!!!!!!Eij!!!!!!!!
      
      include 'pipirad/decomposition5pt/Eij.F'
      
!!!!!!!Ritorno ai valori nominali delle masse perchè per qualche valore
!!!!!!!di disp, modifico i valori
      
      m1q  =  masses(1)*(1.d0,0.d0)
      m2q  =  masses(2)*(1.d0,0.d0)
      m3q  =  masses(3)*(1.d0,0.d0)
      m4q  =  masses(4)*(1.d0,0.d0)
      m5q  =  masses(5)*(1.d0,0.d0)
      
      return
      end
      


*-----------------------------------------
      function lambdad(s,sp,ph2)
      implicit none
      real*8 lambdad
      real*8 s,sp,ph2

      lambdad=(s-sp+ph2)**2-4.d0*s*ph2
      return
      end
*-----------------------------------------      
      function ImD0(p1,p3,s,sp,ph2)
      implicit none
      real*8 p1,p3
      dimension p1(0:3),p3(0:3)
      double complex ImD0
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      
      real*8 mp2,me2,f,ph2,t,s,sp
      real*8 argsqrt,num,den

      double complex im
      
      real*8 dot
      external dot

      external lambdad
      real*8 lambdad

      if(sp.ge.(s+ph2-2.d0*dsqrt(s)*dsqrt(ph2))) then
         ImD0= (0.d0,0.d0)
         return
      endif
      
      im = (0.d0,1.d0)
      mp2=ampi**2
      me2=ame**2

      t = dot(p1-p3,p1-p3)

      argsqrt = lambdad(s,sp,ph2)*((mp2-me2-t)**2-4.d0*me2*t)
     >     -4.d0*ph2*sp*s*t

      num = 2.d0*ph2*s*sp+lambdad(s,sp,ph2)*(me2+mp2-t
     >     +dsqrt(argsqrt/lambdad(s,sp,ph2)))

      den = 2.d0*ph2*s*sp+lambdad(s,sp,ph2)*(me2+mp2-t
     >     -dsqrt(argsqrt/lambdad(s,sp,ph2)))

      ImD0 = pi/dsqrt(argsqrt)*dlog(num/den)

      ImD0 = ImD0*im
      return
      end

      
      real*8 function Itot1ovsp_lambdarad(p1,p3,s,lambda2)
      implicit none
      real*8 p1,p3
      dimension p1(0:3),p3(0:3)
      real*8 I2Moretti,s,T,I1moretti
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi

      real*8 p1p3,p1p2,e1,be1,e3,be3,esoft,lambda2
      complex*16 i13

      real*8 mp2,me2,f

      real*8 num,den
      real*8 dot
      external dot


      mp2=ampi**2
      me2=ame**2

      t = dot(p1-p3,p1-p3)

      f= sqrt(t**2+(me2-mp2)**2-2d0*t*(me2+mp2))
      num=mp2+me2-t+F
      den=mp2+me2-t-F
      I2moretti=pi/s*log(s/(4d0*mp2))*1d0/f*log(num/den)

      p1p3=dot(p1,p3)
      
      e1 = p1(0)
      e3 = p3(0)

      be1=sqrt(1.d0-me2/e1**2)
      be3=sqrt(1.d0-mp2/e3**2)
      esoft=(s-4*mp2)/(2.d0*sqrt(s))

      call intpipj(pi,p1p3,e1,be1,ame,e3,be3,ampi,esoft,lambda2,i13)
      I1moretti=real(i13)/p1p3/(4d0*s)

      Itot1ovsp_lambdarad=I1moretti+I2moretti

      return
      end


      real*8 function Itotspovsp_lambdarad(p1,p3,s,lambda2)
      implicit none
      real*8 p1,p2,p3,p4
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      real*8 s,T

      real*8 p1p3,p1p2,e1,be1,e3,be3,esoft,lambda2
      complex*16 i13
      real*8 ame,ammu,convfac,alpha,pi,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi

      real*8 me2,mp2
      real*8 dot
      external dot

      mp2=ampi**2
      me2=ame**2

      t = dot(p1-p3,p1-p3)

      p1p3=dot(p1,p3)
      
      e1 = p1(0)
      e3 = p3(0)

      be1=sqrt(1.d0-me2/e1**2)
      be3=sqrt(1.d0-mp2/e3**2)
      esoft=(s-4*mp2)/(2.d0*sqrt(s))

      call intpipj(pi,p1p3,e1,be1,ame,e3,be3,ampi,esoft,lambda2,i13)
      Itotspovsp_lambdarad=real(i13)/p1p3/(4d0)

      return
      end


*** the following were in routines.F in Marco's version
*********
      subroutine routing_mumu_old(amme,mmu,pt1,pt2,pt3,pt4,qt,
     >				 momvec,mominv,masses,i,elmat,ifc)

#ifdef COLLIER
        use collier
#endif
        implicit none

	double precision pt1(0:3),pt2(0:3),pt3(0:3),pt4(0:3),qt(0:3)
	double precision pt(0:3)
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	double precision amme,mmu,s,s1,me2,me4,mm2,mm4

        external vpolc
        double complex vpolc,alphar

        common/ialpharunning/iarun
        integer iarun


	double complex tens(1:35)
	double complex tensuv(1:35)

	character*3 ifc

	double complex masses(0:4)
	double complex mominv(1:10)
	double complex momvec(0:3,4)

	double precision E0L
	double precision Evec(0:3)
	double precision Er2(0:9),Er2m(0:3,0:3)
	double precision Er3(0:19),Er3m(0:3,0:3,0:3)

	real*8 elmat,elmat1

	integer i,k,l,m,n,o,p

	double precision p1a
	double precision p2a
	double precision p3a
	double precision p4a
	double precision ka

	double precision gab

	double precision p1ap1b
	double precision p1ap2b
	double precision p1ap3b
	double precision p1ap4b
	double precision kap1b
	double precision p2ap2b
	double precision p2ap3b
	double precision p2ap4b
	double precision kap2b
	double precision p3ap3b
	double precision p3ap4b
	double precision kap3b
	double precision p4ap4b
	double precision kap4b
	double precision kakb

	double precision p1ap1bp2c
	double precision p1ap1bp3c
	double precision p1ap1bp4c
	double precision kap1bp1c
	double precision p1ap2bp2c
	double precision p1ap2bp3c
	double precision p1ap2bp4c
	double precision kap1bp2c
	double precision p1ap3bp3c
	double precision p1ap3bp4c
	double precision kap1bp3c
	double precision p1ap4bp4c
	double precision kap1bp4c
	double precision p2ap2bp3c
	double precision p2ap2bp4c
	double precision kap2bp2c
	double precision p2ap3bp3c
	double precision p2ap3bp4c
	double precision kap2bp3c
	double precision p2ap4bp4c
	double precision kap2bp4c
	double precision p3ap3bp4c
	double precision kap3bp3c
	double precision p3ap4bp4c
	double precision kap3bp4c
	double precision kap4bp4c

	double precision gabkc
	double precision gabp1c
	double precision gabp2c
	double precision gabp3c
	double precision gabp4c

	external dot
	double precision dot



	p1p2 = dot(pt1,pt2)
	p1p3 = dot(pt1,pt3)
	p1p4 = dot(pt1,pt4)
	kp1  = dot(pt1,qt)
	p2p3 = dot(pt2,pt3)
	p2p4 = dot(pt2,pt4)
	kp2  = dot(pt2,qt)
	kp3 = dot(pt3,qt)
	p3p4 = dot(pt3,pt4)
	kp4  = dot(pt4,qt)
	s=dot(pt1+pt2,pt1+pt2)
	s1=dot(pt3+pt4,pt3+pt4)



#ifdef COLLIER

	call Eten_cll(tens,tensuv,momvec,mominv,masses,3)

#endif



************QUI PRENDIAMO I VARI ELEMENTI DI TENS E ISOLIAMO QUELLO SCALARE, VETTORIALE, TENSORIALE DI RANK2 E 3

	E0L=real(tens(1))

	Evec(0)=real(tens(2))
	Evec(1)=real(tens(3))
	Evec(2)=real(tens(4))
	Evec(3)=real(tens(5))

	Er2(0)=real(tens(6))

	Er2(1)=real(tens(7))

	Er2(2)=real(tens(8))

	Er2(3)=real(tens(9))
	Er2(4)=real(tens(10))
	Er2(5)=real(tens(11))
	Er2(6)=real(tens(12))
	Er2(7)=real(tens(13))
	Er2(8)=real(tens(14))
	Er2(9)=real(tens(15))
        if(ifc.eq.'is1') print*,'vettoriali',E0L,Evec
      Er2m(0,0)=Er2(0)
      Er2m(0,1)=Er2(1)
      Er2m(1,0)=Er2(1)
      Er2m(2,0)=Er2(2)
      Er2m(0,2)=Er2(2)
      Er2m(3,0)=Er2(3)
      Er2m(0,3)=Er2(3)
      Er2m(1,1)=Er2(4)
      Er2m(1,2)=Er2(5)
      Er2m(2,1)=Er2(5)
      Er2m(1,3)=Er2(6)
      Er2m(3,1)=Er2(6)
      Er2m(2,2)=Er2(7)
      Er2m(2,3)=Er2(8)
      Er2m(3,2)=Er2(8)
      Er2m(3,3)=Er2(9)

	Er3(0)=real(tens(16))
	Er3(1)=real(tens(17))
	Er3(2)=real(tens(18))
	Er3(3)=real(tens(19))
	Er3(4)=real(tens(20))
	Er3(5)=real(tens(21))
	Er3(6)=real(tens(22))
	Er3(7)=real(tens(23))
	Er3(8)=real(tens(24))
	Er3(9)=real(tens(25))
	Er3(10)=real(tens(26))
	Er3(11)=real(tens(27))
	Er3(12)=real(tens(28))
	Er3(13)=real(tens(29))
	Er3(14)=real(tens(30))
	Er3(15)=real(tens(31))
	Er3(16)=real(tens(32))
	Er3(17)=real(tens(33))
	Er3(18)=real(tens(34))
	Er3(19)=real(tens(35))

	Er3m(0, 0, 0) = Er3(0)
	Er3m(0, 0, 1) = Er3(1)
	Er3m(0, 1, 0) = Er3(1)
	Er3m(1, 0, 0) = Er3(1)
	Er3m(0, 0, 2) = Er3(2)
	Er3m(0, 2, 0) = Er3(2)
	Er3m(2, 0, 0) = Er3(2)
	Er3m(0, 0, 3) = Er3(3)
	Er3m(0, 3, 0) = Er3(3)
	Er3m(3, 0, 0) = Er3(3)
	Er3m(0, 1, 1) = Er3(4)
	Er3m(1, 0, 1) = Er3(4)
	Er3m(1, 1, 0) = Er3(4)
	Er3m(0, 1, 2) = Er3(5)
	Er3m(0, 2, 1) = Er3(5)
	Er3m(1, 0, 2) = Er3(5)
	Er3m(1, 2, 0) = Er3(5)
	Er3m(2, 0, 1) = Er3(5)
	Er3m(2, 1, 0) = Er3(5)
	Er3m(0, 1, 3) = Er3(6)
	Er3m(0, 3, 1) = Er3(6)
	Er3m(1, 0, 3) = Er3(6)
	Er3m(1, 3, 0) = Er3(6)
	Er3m(3, 0, 1) = Er3(6)
	Er3m(3, 1, 0) = Er3(6)
	Er3m(0, 2, 2) = Er3(7)
	Er3m(2, 0, 2) = Er3(7)
	Er3m(2, 2, 0) = Er3(7)
	Er3m(0, 2, 3) = Er3(8)
	Er3m(0, 3, 2) = Er3(8)
	Er3m(2, 0, 3) = Er3(8)
	Er3m(2, 3, 0) = Er3(8)
	Er3m(3, 0, 2) = Er3(8)
	Er3m(3, 2, 0) = Er3(8)
	Er3m(0, 3, 3) = Er3(9)
	Er3m(3, 0, 3) = Er3(9)
	Er3m(3, 3, 0) = Er3(9)
	Er3m(1, 1, 1) = Er3(10)
	Er3m(1, 1, 2) = Er3(11)
	Er3m(1, 2, 1) = Er3(11)
	Er3m(2, 1, 1) = Er3(11)
	Er3m(1, 1, 3) = Er3(12)
	Er3m(1, 3, 1) = Er3(12)
	Er3m(3, 1, 1) = Er3(12)
	Er3m(1, 2, 2) = Er3(13)
	Er3m(2, 1, 2) = Er3(13)
	Er3m(2, 2, 1) = Er3(13)
	Er3m(1, 2, 3) = Er3(14)
	Er3m(1, 3, 2) = Er3(14)
	Er3m(2, 1, 3) = Er3(14)
	Er3m(2, 3, 1) = Er3(14)
	Er3m(3, 1, 2) = Er3(14)
	Er3m(3, 2, 1) = Er3(14)
	Er3m(1, 3, 3) = Er3(15)
	Er3m(3, 1, 3) = Er3(15)
	Er3m(3, 3, 1) = Er3(15)
	Er3m(2, 2, 2) = Er3(16)
	Er3m(2, 2, 3) = Er3(17)
	Er3m(2, 3, 2) = Er3(17)
	Er3m(3, 2, 2) = Er3(17)
	Er3m(2, 3, 3) = Er3(18)
	Er3m(3, 2, 3) = Er3(18)
	Er3m(3, 3, 2) = Er3(18)
	Er3m(3, 3, 3) = Er3(19)

***************vettoriali

	p1a=dot(pt1,Evec)
	p2a=dot(pt2,Evec)
	p3a=dot(pt3,Evec)
	p4a=dot(pt4,Evec)
	ka=dot(qt,Evec)

*********rango 2

	call rank2(Er2,pt3,qt,2,gab)

	call rank2(Er2,pt1,pt1,1,p1ap1b)

	call rank2(Er2,pt1,pt2,1,p1ap2b)

	call rank2(Er2,pt1,pt3,1,p1ap3b)

	call rank2(Er2,pt1,pt4,1,p1ap4b)

	call rank2(Er2,pt1,qt,1,kap1b)

	call rank2(Er2,pt2,pt2,1,p2ap2b)

	call rank2(Er2,pt2,pt3,1,p2ap3b)

	call rank2(Er2,pt2,pt4,1,p2ap4b)

	call rank2(Er2,pt2,qt,1,kap2b)

	call rank2(Er2,pt3,pt3,1,p3ap3b)

	call rank2(Er2,pt3,pt4,1,p3ap4b)

	call rank2(Er2,pt3,qt,1,kap3b)

	call rank2(Er2,pt4,pt4,1,p4ap4b)

	call rank2(Er2,pt4,qt,1,kap4b)

	call rank2(Er2,qt,qt,1,kakb)

**********************rango 3

	call rank3(Er3,pt1,pt1,pt2,1,p1ap1bp2c)

	call rank3(Er3,pt1,pt1,pt3,1,p1ap1bp3c)

	call rank3(Er3,pt1,pt1,pt4,1,p1ap1bp4c)

	call rank3(Er3,pt1,pt1,qt,1,kap1bp1c)

	call rank3(Er3,pt1,pt2,pt2,1,p1ap2bp2c)

	call rank3(Er3,pt1,pt2,pt3,1,p1ap2bp3c)

	call rank3(Er3,pt1,pt2,pt4,1,p1ap2bp4c)

	call rank3(Er3,pt1,pt2,qt,1,kap1bp2c)

	call rank3(Er3,pt1,pt3,pt3,1, p1ap3bp3c)

	call rank3(Er3,pt1,pt3,pt4,1, p1ap3bp4c)

	call rank3(Er3,pt1,pt3,qt,1, kap1bp3c)

	call rank3(Er3,pt1,pt4,pt4,1, p1ap4bp4c)

	call rank3(Er3,pt1,pt4,qt,1, kap1bp4c)

	call rank3(Er3,pt2,pt2,pt3,1, p2ap2bp3c)

	call rank3(Er3,pt2,pt2,pt4,1, p2ap2bp4c)

	call rank3(Er3,pt2,pt2,qt,1, kap2bp2c)

	call rank3(Er3,pt2,pt3,pt3,1, p2ap3bp3c)

	call rank3(Er3,pt2,pt3,pt4,1, p2ap3bp4c)

	call rank3(Er3,pt2,pt3,qt,1, kap2bp3c)

	call rank3(Er3,pt2,pt4,pt4,1, p2ap4bp4c)

	call rank3(Er3,pt2,pt4,qt,1,  kap2bp4c)

	call rank3(Er3,pt3,pt3,pt4,1, p3ap3bp4c)

	call rank3(Er3,pt3,pt3,qt,1, kap3bp3c)

	call rank3(Er3,pt3,pt4,pt4,1, p3ap4bp4c)

	call rank3(Er3,pt3,pt4,qt,1, kap3bp4c)


	call rank3(Er3,pt4,pt4,qt,1, kap4bp4c)

	call rank3(Er3,pt1,pt4,qt,2, gabkc)

	call rank3(Er3,pt1,pt4,pt1,2, gabp1c)

	call rank3(Er3,pt1,pt4,pt2,2, gabp2c)

	call rank3(Er3,pt1,pt4,pt3,2, gabp3c)

	call rank3(Er3,pt1,pt4,pt4,2, gabp4c)


	elmat=0d0
*****SOLO IL ROUTING 5 È USATO, PERÒ LASCIO COMUNQUE IL CODICE COMMENTATO SE MAI SERVISSERO DEI TEST

	if(i.eq.1) then
	    include 'mumurad/Virtual/pent_isr_isr_tensor/Routing1.F'
            if(ifc.eq.'is1'.or.ifc.eq.'fs1') then

        	elmat=elmat-elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    else
		elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    endif


    	    include 'mumurad/Virtual/pent_isr_fsr_tensor/Routing1.F'
	    elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s))

	elseif(i.eq.2) then
	    include 'mumurad/Virtual/pent_isr_isr_tensor/Routing2.F'
            if(ifc.eq.'is1'.or.ifc.eq.'fs1') then
        	elmat=elmat-elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    else
		elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    endif

    	    include 'mumurad/Virtual/pent_isr_fsr_tensor/Routing2.F'
	    elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s))
	elseif(i.eq.3) then
  	    include 'mumurad/Virtual/pent_isr_isr_tensor/Routing3.F'
            if(ifc.eq.'is1'.or.ifc.eq.'fs1') then
        	elmat=elmat-elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    else
		elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    endif

    	    include 'mumurad/Virtual/pent_isr_fsr_tensor/Routing3.F'
	    elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s))
	elseif(i.eq.4) then
    	    include 'mumurad/Virtual/pent_isr_isr_tensor/Routing4.F'
            if(ifc.eq.'is1'.or.ifc.eq.'fs1') then
        	elmat=elmat-elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    else
		elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    endif
    	    include 'mumurad/Virtual/pent_isr_fsr_tensor/Routing4.F'
	    elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s))
!
	elseif(i.eq.5) then
	    include 'mumurad/Virtual/pent_isr_isr_tensor/Routing5.F'

            if(ifc.eq.'is1'.or.ifc.eq.'fs1') then

        	elmat=elmat-elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    else
		elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    endif
    	    include 'mumurad/Virtual/pent_isr_fsr_tensor/Routing5.F'
	    elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s))
	endif

      end
************************
      subroutine routing_mumu(amme,mmu,pt1,pt2,pt3,pt4,qt,
     >				 momvec,mominv,masses,i,elmat,ifc)

#ifdef COLLIER
        use collier
#endif
        implicit none

	double precision pt1(0:3),pt2(0:3),pt3(0:3),pt4(0:3),qt(0:3)
	double precision pt(0:3)
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	double precision amme,mmu,s,s1,me2,me4,mm2,mm4

        external vpolc
        double complex vpolc,alphar

        common/ialpharunning/iarun
        integer iarun


	double complex tens(1:35)
	double complex tensuv(1:35)

	character*3 ifc

	double complex masses(0:4)
	double complex mominv(1:10)
	double complex momvec(0:3,4)

	double complex E0L
	double complex Evec(0:3)
	double complex Er2(0:9),Er2m(0:3,0:3)
	double complex Er3(0:19),Er3m(0:3,0:3,0:3)

	double complex elmat,elmat1

	integer i,k,l,m,n,o,p

	double complex p1a
	double complex p2a
	double complex p3a
	double complex p4a
	double complex ka

        DOUBLE COMPLEX gab

        DOUBLE COMPLEX p1ap1b
        DOUBLE COMPLEX p1ap2b
        DOUBLE COMPLEX p1ap3b
        DOUBLE COMPLEX p1ap4b
        DOUBLE COMPLEX kap1b
        DOUBLE COMPLEX p2ap2b
        DOUBLE COMPLEX p2ap3b
        DOUBLE COMPLEX p2ap4b
        DOUBLE COMPLEX kap2b
        DOUBLE COMPLEX p3ap3b
        DOUBLE COMPLEX p3ap4b
        DOUBLE COMPLEX kap3b
        DOUBLE COMPLEX p4ap4b
        DOUBLE COMPLEX kap4b
        DOUBLE COMPLEX kakb

        DOUBLE COMPLEX p1ap1bp2c
        DOUBLE COMPLEX p1ap1bp3c
        DOUBLE COMPLEX p1ap1bp4c
        DOUBLE COMPLEX kap1bp1c
        DOUBLE COMPLEX p1ap2bp2c
        DOUBLE COMPLEX p1ap2bp3c
        DOUBLE COMPLEX p1ap2bp4c
        DOUBLE COMPLEX kap1bp2c
        DOUBLE COMPLEX p1ap3bp3c
        DOUBLE COMPLEX p1ap3bp4c
        DOUBLE COMPLEX kap1bp3c
        DOUBLE COMPLEX p1ap4bp4c
        DOUBLE COMPLEX kap1bp4c
        DOUBLE COMPLEX p2ap2bp3c
        DOUBLE COMPLEX p2ap2bp4c
        DOUBLE COMPLEX kap2bp2c
        DOUBLE COMPLEX p2ap3bp3c
        DOUBLE COMPLEX p2ap3bp4c
        DOUBLE COMPLEX kap2bp3c
        DOUBLE COMPLEX p2ap4bp4c
        DOUBLE COMPLEX kap2bp4c
        DOUBLE COMPLEX p3ap3bp4c
        DOUBLE COMPLEX kap3bp3c
        DOUBLE COMPLEX p3ap4bp4c
        DOUBLE COMPLEX kap3bp4c
        DOUBLE COMPLEX kap4bp4c

        DOUBLE COMPLEX gabkc
        DOUBLE COMPLEX gabp1c
        DOUBLE COMPLEX gabp2c
        DOUBLE COMPLEX gabp3c
        DOUBLE COMPLEX gabp4c

	external dot
	double precision dot



	p1p2 = dot(pt1,pt2)
	p1p3 = dot(pt1,pt3)
	p1p4 = dot(pt1,pt4)
	kp1  = dot(pt1,qt)
	p2p3 = dot(pt2,pt3)
	p2p4 = dot(pt2,pt4)
	kp2  = dot(pt2,qt)
	kp3 = dot(pt3,qt)
	p3p4 = dot(pt3,pt4)
	kp4  = dot(pt4,qt)
	s=dot(pt1+pt2,pt1+pt2)
	s1=dot(pt3+pt4,pt3+pt4)



#ifdef COLLIER

	call Eten_cll(tens,tensuv,momvec,mominv,masses,3)

#endif



************QUI PRENDIAMO I VARI ELEMENTI DI TENS E ISOLIAMO QUELLO SCALARE, VETTORIALE, TENSORIALE DI RANK2 E 3
        E0L = tens(1)

        Evec(0) = tens(2)
        Evec(1) = tens(3)
        Evec(2) = tens(4)
        Evec(3) = tens(5)

        Er2(0) = tens(6)
        Er2(1) = tens(7)
        Er2(2) = tens(8)
        Er2(3) = tens(9)
        Er2(4) = tens(10)
        Er2(5) = tens(11)
        Er2(6) = tens(12)
        Er2(7) = tens(13)
        Er2(8) = tens(14)
        Er2(9) = tens(15)

      Er2m(0,0)=Er2(0)
      Er2m(0,1)=Er2(1)
      Er2m(1,0)=Er2(1)
      Er2m(2,0)=Er2(2)
      Er2m(0,2)=Er2(2)
      Er2m(3,0)=Er2(3)
      Er2m(0,3)=Er2(3)
      Er2m(1,1)=Er2(4)
      Er2m(1,2)=Er2(5)
      Er2m(2,1)=Er2(5)
      Er2m(1,3)=Er2(6)
      Er2m(3,1)=Er2(6)
      Er2m(2,2)=Er2(7)
      Er2m(2,3)=Er2(8)
      Er2m(3,2)=Er2(8)
      Er2m(3,3)=Er2(9)

      Er3(0) = tens(16)
      Er3(1) = tens(17)
      Er3(2) = tens(18)
      Er3(3) = tens(19)
      Er3(4) = tens(20)
      Er3(5) = tens(21)
      Er3(6) = tens(22)
      Er3(7) = tens(23)
      Er3(8) = tens(24)
      Er3(9) = tens(25)
      Er3(10) = tens(26)
      Er3(11) = tens(27)
      Er3(12) = tens(28)
      Er3(13) = tens(29)
      Er3(14) = tens(30)
      Er3(15) = tens(31)
      Er3(16) = tens(32)
      Er3(17) = tens(33)
      Er3(18) = tens(34)
      Er3(19) = tens(35)

	Er3m(0, 0, 0) = Er3(0)
	Er3m(0, 0, 1) = Er3(1)
	Er3m(0, 1, 0) = Er3(1)
	Er3m(1, 0, 0) = Er3(1)
	Er3m(0, 0, 2) = Er3(2)
	Er3m(0, 2, 0) = Er3(2)
	Er3m(2, 0, 0) = Er3(2)
	Er3m(0, 0, 3) = Er3(3)
	Er3m(0, 3, 0) = Er3(3)
	Er3m(3, 0, 0) = Er3(3)
	Er3m(0, 1, 1) = Er3(4)
	Er3m(1, 0, 1) = Er3(4)
	Er3m(1, 1, 0) = Er3(4)
	Er3m(0, 1, 2) = Er3(5)
	Er3m(0, 2, 1) = Er3(5)
	Er3m(1, 0, 2) = Er3(5)
	Er3m(1, 2, 0) = Er3(5)
	Er3m(2, 0, 1) = Er3(5)
	Er3m(2, 1, 0) = Er3(5)
	Er3m(0, 1, 3) = Er3(6)
	Er3m(0, 3, 1) = Er3(6)
	Er3m(1, 0, 3) = Er3(6)
	Er3m(1, 3, 0) = Er3(6)
	Er3m(3, 0, 1) = Er3(6)
	Er3m(3, 1, 0) = Er3(6)
	Er3m(0, 2, 2) = Er3(7)
	Er3m(2, 0, 2) = Er3(7)
	Er3m(2, 2, 0) = Er3(7)
	Er3m(0, 2, 3) = Er3(8)
	Er3m(0, 3, 2) = Er3(8)
	Er3m(2, 0, 3) = Er3(8)
	Er3m(2, 3, 0) = Er3(8)
	Er3m(3, 0, 2) = Er3(8)
	Er3m(3, 2, 0) = Er3(8)
	Er3m(0, 3, 3) = Er3(9)
	Er3m(3, 0, 3) = Er3(9)
	Er3m(3, 3, 0) = Er3(9)
	Er3m(1, 1, 1) = Er3(10)
	Er3m(1, 1, 2) = Er3(11)
	Er3m(1, 2, 1) = Er3(11)
	Er3m(2, 1, 1) = Er3(11)
	Er3m(1, 1, 3) = Er3(12)
	Er3m(1, 3, 1) = Er3(12)
	Er3m(3, 1, 1) = Er3(12)
	Er3m(1, 2, 2) = Er3(13)
	Er3m(2, 1, 2) = Er3(13)
	Er3m(2, 2, 1) = Er3(13)
	Er3m(1, 2, 3) = Er3(14)
	Er3m(1, 3, 2) = Er3(14)
	Er3m(2, 1, 3) = Er3(14)
	Er3m(2, 3, 1) = Er3(14)
	Er3m(3, 1, 2) = Er3(14)
	Er3m(3, 2, 1) = Er3(14)
	Er3m(1, 3, 3) = Er3(15)
	Er3m(3, 1, 3) = Er3(15)
	Er3m(3, 3, 1) = Er3(15)
	Er3m(2, 2, 2) = Er3(16)
	Er3m(2, 2, 3) = Er3(17)
	Er3m(2, 3, 2) = Er3(17)
	Er3m(3, 2, 2) = Er3(17)
	Er3m(2, 3, 3) = Er3(18)
	Er3m(3, 2, 3) = Er3(18)
	Er3m(3, 3, 2) = Er3(18)
	Er3m(3, 3, 3) = Er3(19)

***************vettoriali
        p1a=dot(pt1,REAL(Evec))+(0.d0,1.d0)*dot(pt1,AIMAG(Evec))
        p2a=dot(pt2,REAL(Evec))+(0.d0,1.d0)*dot(pt2,AIMAG(Evec))
	p3a=dot(pt3,REAL(Evec))+(0.d0,1.d0)*dot(pt3,AIMAG(Evec))
	p4a=dot(pt4,REAL(Evec))+(0.d0,1.d0)*dot(pt4,AIMAG(Evec))
        ka=dot(qt,REAL(Evec))+(0.d0,1.d0)*dot(qt,AIMAG(Evec))

*********rango 2

        call rank2_c(Er2,pt3,qt,2,gab)

        call rank2_c(Er2,pt1,pt1,1,p1ap1b)

        call rank2_c(Er2,pt1,pt2,1,p1ap2b)

        call rank2_c(Er2,pt1,pt3,1,p1ap3b)

        call rank2_c(Er2,pt1,pt4,1,p1ap4b)

        call rank2_c(Er2,pt1,qt,1,kap1b)

        call rank2_c(Er2,pt2,pt2,1,p2ap2b)

        call rank2_c(Er2,pt2,pt3,1,p2ap3b)

        call rank2_c(Er2,pt2,pt4,1,p2ap4b)

        call rank2_c(Er2,pt2,qt,1,kap2b)

        call rank2_c(Er2,pt3,pt3,1,p3ap3b)

        call rank2_c(Er2,pt3,pt4,1,p3ap4b)

        call rank2_c(Er2,pt3,qt,1,kap3b)

        call rank2_c(Er2,pt4,pt4,1,p4ap4b)

        call rank2_c(Er2,pt4,qt,1,kap4b)

        call rank2_c(Er2,qt,qt,1,kakb)

**********************rango 3

        call rank3_c(Er3,pt1,pt1,pt2,1,p1ap1bp2c)

        call rank3_c(Er3,pt1,pt1,pt3,1,p1ap1bp3c)

        call rank3_c(Er3,pt1,pt1,pt4,1,p1ap1bp4c)

        call rank3_c(Er3,pt1,pt1,qt,1,kap1bp1c)

        call rank3_c(Er3,pt1,pt2,pt2,1,p1ap2bp2c)

        call rank3_c(Er3,pt1,pt2,pt3,1,p1ap2bp3c)

        call rank3_c(Er3,pt1,pt2,pt4,1,p1ap2bp4c)

        call rank3_c(Er3,pt1,pt2,qt,1,kap1bp2c)

        call rank3_c(Er3,pt1,pt3,pt3,1,p1ap3bp3c)

        call rank3_c(Er3,pt1,pt3,pt4,1,p1ap3bp4c)

        call rank3_c(Er3,pt1,pt3,qt,1,kap1bp3c)

        call rank3_c(Er3,pt1,pt4,pt4,1,p1ap4bp4c)

        call rank3_c(Er3,pt1,pt4,qt,1,kap1bp4c)

        call rank3_c(Er3,pt2,pt2,pt3,1,p2ap2bp3c)

        call rank3_c(Er3,pt2,pt2,pt4,1,p2ap2bp4c)

        call rank3_c(Er3,pt2,pt2,qt,1,kap2bp2c)

        call rank3_c(Er3,pt2,pt3,pt3,1,p2ap3bp3c)

        call rank3_c(Er3,pt2,pt3,pt4,1,p2ap3bp4c)

        call rank3_c(Er3,pt2,pt3,qt,1,kap2bp3c)

        call rank3_c(Er3,pt2,pt4,pt4,1,p2ap4bp4c)

        call rank3_c(Er3,pt2,pt4,qt,1,kap2bp4c)

        call rank3_c(Er3,pt3,pt3,pt4,1,p3ap3bp4c)

        call rank3_c(Er3,pt3,pt3,qt,1,kap3bp3c)

        call rank3_c(Er3,pt3,pt4,pt4,1,p3ap4bp4c)

        call rank3_c(Er3,pt3,pt4,qt,1,kap3bp4c)

        call rank3_c(Er3,pt4,pt4,qt,1,kap4bp4c)

        call rank3_c(Er3,pt1,pt4,qt,2,gabkc)

        call rank3_c(Er3,pt1,pt4,pt1,2,gabp1c)

        call rank3_c(Er3,pt1,pt4,pt2,2,gabp2c)

        call rank3_c(Er3,pt1,pt4,pt3,2,gabp3c)

        call rank3_c(Er3,pt1,pt4,pt4,2,gabp4c)


	elmat=0d0
*****SOLO IL ROUTING 5 È USATO, PERÒ LASCIO COMUNQUE IL CODICE COMMENTATO SE MAI SERVISSERO DEI TEST

	if(i.eq.1) then
	    include 'mumurad/Virtual/pent_isr_isr_tensor/Routing1.F'
            if(ifc.eq.'is1'.or.ifc.eq.'fs1') then

        	elmat=elmat-elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    else
		elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    endif


    	    include 'mumurad/Virtual/pent_isr_fsr_tensor/Routing1.F'
	    elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s))

	elseif(i.eq.2) then
	    include 'mumurad/Virtual/pent_isr_isr_tensor/Routing2.F'
            if(ifc.eq.'is1'.or.ifc.eq.'fs1') then
        	elmat=elmat-elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    else
		elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    endif

    	    include 'mumurad/Virtual/pent_isr_fsr_tensor/Routing2.F'
	    elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s))
	elseif(i.eq.3) then
  	    include 'mumurad/Virtual/pent_isr_isr_tensor/Routing3.F'
            if(ifc.eq.'is1'.or.ifc.eq.'fs1') then
        	elmat=elmat-elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    else
		elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    endif

    	    include 'mumurad/Virtual/pent_isr_fsr_tensor/Routing3.F'
	    elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s))
	elseif(i.eq.4) then
    	    include 'mumurad/Virtual/pent_isr_isr_tensor/Routing4.F'
            if(ifc.eq.'is1'.or.ifc.eq.'fs1') then
        	elmat=elmat-elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    else
		elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    endif
    	    include 'mumurad/Virtual/pent_isr_fsr_tensor/Routing4.F'
	    elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s))
!
	elseif(i.eq.5) then
	    include 'mumurad/Virtual/pent_isr_isr_tensor/Routing5.F'

            if(ifc.eq.'is1'.or.ifc.eq.'fs1') then

        	elmat=elmat-elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    else
		elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s1))
	    endif
    	    include 'mumurad/Virtual/pent_isr_fsr_tensor/Routing5.F'
	    elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s))
	endif
      end
************************


*************
      subroutine routing_pipi(pt1,pt2,pt3,pt4,qt,
     >		     momvec,mominv,masses,i,j,elmati,elmatf)

#ifdef COLLIER
        use collier
#endif
        implicit none

	double precision pt1(0:3),pt2(0:3),pt3(0:3),pt4(0:3),qt(0:3)
	double precision pt(0:3)
	real*8 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	double precision ame,ammu,convfac,alpha,pi,ame2,ampi2,ampi
	common/parameters/ame,ammu,convfac,alpha,pi
        common/hparam/ampi

	double complex tens(1:15)
	double complex tensuv(1:15)

	character*3 ifc

	double complex masses(0:4)
	double complex mominv(1:10)
	double complex momvec(0:3,4)

	double complex E0L
	double complex Evec(0:3)
	double complex Er2(0:9),Er2m(0:3,0:3)

	double complex elmati,elmatf,elmat1

	integer i,j,k,l,m,n,o,p

	double complex p1a
	double complex p2a
	double complex p3a
	double complex p4a
	double complex ka

	double complex gab

	double complex p1ap1b
	double complex p1ap2b
	double complex p1ap3b
	double complex p1ap4b
	double complex kap1b
	double complex p2ap2b
	double complex p2ap3b
	double complex p2ap4b
	double complex kap2b
	double complex p3ap3b
	double complex p3ap4b
	double complex kap3b
	double complex p4ap4b
	double complex kap4b
	double complex kakb


	external dot
	double precision dot

        double complex vpolc,alphar

        common/ialpharunning/iarun
        integer iarun

**********FF

	integer iffpi
	character*10 what_ffpi
	common/ffpioption/iffpi,what_ffpi
	double complex ffpiis,Cffpiis
	double precision ffpiis2
	double complex ffpifs,Cffpifs
	double precision ffpifs2

**************

*************GVMD COMMON
        common/gvmdrad/indexnores
        integer indexnores

       	p1p2 = dot(pt1,pt2)
	p1p3 = dot(pt1,pt3)
	p1p4 = dot(pt1,pt4)
	kp1  = dot(pt1,qt)
	p2p3 = dot(pt2,pt3)
	p2p4 = dot(pt2,pt4)
	kp2  = dot(pt2,qt)
	kp3 = dot(pt3,qt)
	p3p4 = dot(pt3,pt4)
	kp4  = dot(pt4,qt)

**************
	ffpifs=(1.d0,0.d0)
	ffpifs2=1.d0

	ffpiis=(1.d0,0.d0)
	ffpiis2=1.d0
	if(iffpi.gt.0) then
		call ffpi(2.d0*(ame**2+p1p2),ffpifs,ffpifs2)
		call ffpi(2.d0*(ampi**2+p3p4),ffpiis,ffpiis2)

	endif	

        Cffpiis=CONJG(ffpiis)
        Cffpifs=CONJG(ffpifs)

        if(indexnores.eq.1) then
           ffpiis=(1.d0,0.d0)
           ffpifs=(1.d0,0.d0)
	endif

******ALPHA RUNNING
	ffpiis=ffpiis*vpolc(2.d0*(ampi**2+p3p4))
        Cffpiis=Cffpiis*CONJG(vpolc(2.d0*(ampi**2+p3p4)))

	ffpifs=ffpifs*vpolc(2.d0*(ame**2+p1p2))
        Cffpifs=Cffpifs*CONJG(vpolc(2.d0*(ame**2+p1p2)))
        
        ffpiis2=ffpiis*Cffpiis
        ffpifs2=ffpifs*Cffpifs
****************

#ifdef COLLIER

	call Eten_cll(tens,tensuv,momvec,mominv,masses,2)
!        call SetMode_cll(2)
#endif



************QUI PRENDIAMO I VARI ELEMENTI DI TENS E ISOLIAMO QUELLO SCALARE, VETTORIALE, TENSORIALE DI RANK2 E 3

	E0L=tens(1)

	Evec(0)=tens(2)
	Evec(1)=tens(3)
	Evec(2)=tens(4)
	Evec(3)=tens(5)

	Er2(0)=tens(6)
	Er2(1)=tens(7)
	Er2(2)=tens(8)
	Er2(3)=tens(9)
	Er2(4)=tens(10)
	Er2(5)=tens(11)
	Er2(6)=tens(12)
	Er2(7)=tens(13)
	Er2(8)=tens(14)
	Er2(9)=tens(15)

	Er2m(0,0)=Er2(0)
	Er2m(0,1)=Er2(1)
	Er2m(1,0)=Er2(1)
	Er2m(2,0)=Er2(2)
	Er2m(0,2)=Er2(2)
	Er2m(3,0)=Er2(3)
	Er2m(0,3)=Er2(3)
	Er2m(1,1)=Er2(4)
	Er2m(1,2)=Er2(5)
	Er2m(2,1)=Er2(5)
	Er2m(1,3)=Er2(6)
	Er2m(3,1)=Er2(6)
	Er2m(2,2)=Er2(7)
	Er2m(2,3)=Er2(8)
	Er2m(3,2)=Er2(8)
	Er2m(3,3)=Er2(9)

***************vettoriali

	p1a=dot(pt1,REAL(Evec))+(0.d0,1.d0)*dot(pt1,AIMAG(Evec))
	p2a=dot(pt2,REAL(Evec))+(0.d0,1.d0)*dot(pt2,AIMAG(Evec))
	p3a=dot(pt3,REAL(Evec))+(0.d0,1.d0)*dot(pt3,AIMAG(Evec))
	p4a=dot(pt4,REAL(Evec))+(0.d0,1.d0)*dot(pt4,AIMAG(Evec))
	ka=dot(qt,REAL(Evec))+(0.d0,1.d0)*dot(qt,AIMAG(Evec))


*********rango 2

	call rank2_c(Er2,pt3,qt,2,gab)

	call rank2_c(Er2,pt1,pt1,1,p1ap1b)

	call rank2_c(Er2,pt1,pt2,1,p1ap2b)

	call rank2_c(Er2,pt1,pt3,1,p1ap3b)

	call rank2_c(Er2,pt1,pt4,1,p1ap4b)

	call rank2_c(Er2,pt1,qt,1,kap1b)

	call rank2_c(Er2,pt2,pt2,1,p2ap2b)

	call rank2_c(Er2,pt2,pt3,1,p2ap3b)

	call rank2_c(Er2,pt2,pt4,1,p2ap4b)

	call rank2_c(Er2,pt2,qt,1,kap2b)

	call rank2_c(Er2,pt3,pt3,1,p3ap3b)

	call rank2_c(Er2,pt3,pt4,1,p3ap4b)

	call rank2_c(Er2,pt3,qt,1,kap3b)

	call rank2_c(Er2,pt4,pt4,1,p4ap4b)

	call rank2_c(Er2,pt4,qt,1,kap4b)

	call rank2_c(Er2,qt,qt,1,kakb)




	elmati=(0.d0,0.d0)
	elmatf=(0.d0,0.d0)
	if(i.eq.1) then
	  if(j.eq.1) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/ISR/Pent_Routing1_isr1.F'
	elmati=elmat1*Cffpiis*ffpifs
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/FSR/Pent_Routing1_fsr1.F'
	elmatf=elmat1*ffpifs*Cffpifs
	  elseif(j.eq.3) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing1_isr3.F'
	elmati=elmat1*ffpiis*Cffpiis
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing1_fsr3.F'
	elmatf=elmat1*(ffpiis*Cffpifs)

	  endif
	elseif(i.eq.2) then
	  if(j.eq.1) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/ISR/Pent_Routing2_isr1.F'
	elmati=elmat1*Cffpiis*ffpifs
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/FSR/Pent_Routing2_fsr1.F'
	elmatf=elmat1*ffpifs*Cffpifs
	  elseif(j.eq.3) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing2_isr3.F'
	elmati=elmat1*ffpiis*Cffpiis
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing2_fsr3.F'
	elmatf=elmat1*(ffpiis*Cffpifs)
	  endif
	elseif(i.eq.3) then
	  if(j.eq.1) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/ISR/Pent_Routing3_isr1.F'
	elmati=elmat1*Cffpiis*ffpifs
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/FSR/Pent_Routing3_fsr1.F'
	elmatf=elmat1*ffpifs*Cffpifs
	  elseif(j.eq.3) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing3_isr3.F'
	elmati=elmat1*ffpiis*Cffpiis
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing3_fsr3.F'
	elmatf=elmat1*(ffpiis*Cffpifs)
	  endif
	elseif(i.eq.4) then
	  if(j.eq.1) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/ISR/Pent_Routing4_isr1.F'
	elmati=elmat1*Cffpiis*ffpifs
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/FSR/Pent_Routing4_fsr1.F'
	elmatf=elmat1*ffpifs*Cffpifs
	  elseif(j.eq.3) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing4_isr3.F'
	elmati=elmat1*ffpiis*Cffpiis
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing4_fsr3.F'
	elmatf=elmat1*(ffpiis*Cffpifs)
	  endif
	elseif(i.eq.5) then
	  if(j.eq.1) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/ISR/Pent_Routing5_isr1.F'
	elmati=elmat1*Cffpiis*ffpifs
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/FSR/Pent_Routing5_fsr1.F'
	elmatf=elmat1*ffpifs*Cffpifs
	  elseif(j.eq.3) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing5_isr3.F'
	elmati=elmat1*ffpiis*Cffpiis
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing5_fsr3.F'
	elmatf=elmat1*(ffpiis*Cffpifs)
	  endif
       endif
      end
