*********************************************************************
      subroutine setfilenames(programma,corrections)
#include "shared.F"
      character*100 distname(ndistr),tmp(ndistr),tmp2(ndistr2),
     .     distname2(ndistr2)
      common/filenames/distname,distname2
      character*(*) programma,corrections
      character*2   ec,ac
      character*7   bw

      do k = 1,7
         bw(k:k) = ' '
      enddo
      do k = 1,ndistr
         do i = 1,100
            tmp(k)(i:i) = ' '
         enddo
      enddo
      do k = 1,ndistr2
         do i = 1,100
            tmp2(k)(i:i) = ' '
         enddo
      enddo

      ebeam = ecmsnom/2.d0

      bw='.txt'

******************
! default
      tmp(1)  = 'el_th_'
      tmp(2)  = 'el_en_'
      tmp(3)  = 'acoll_'
      tmp(4)  = 'g1_en_'
      tmp(5)  = 'g1_th_'
      tmp(6)  = 'ep_th_'
      tmp(7)  = 'ep_en_'
      tmp(8)  = 'e+e-m_'
      tmp(9)  = 'g2_en_'
      tmp(10) = 'g2_th_'
      tmp(11) = 'm_jpsi_'
      tmp(12) = 'm_psi2S_'
      tmp(13) = 'm_Y_'
      tmp(14) = 'ctheCRF_'
      tmp(15) = 'ctheCRFf_'
      tmp(16) = 'ctheCRFb_'
      tmp(17) = 'cthe_'
      tmp(18) = 'cthp_'
      tmp(19) = 'q2_34_'
******** 2d distributions
      if (i2d.gt.0) then
         tmp2(1) = '2d_q2_ctheCRF_'
      endif
      do k=1,ndistr2
         if (i2d.gt.0) lung2(k) = len_trim(tmp2(k))
      enddo
******************

      do k=1,ndistr
         lung(k) = len_trim(tmp(k))
      enddo

      if (corrections.eq.'oal') then
         do k=1,ndistr
            tmp(k)(lung(k)+1:) = 'NLO'//bw
         enddo
         do k=1,ndistr2
            if (i2d.gt.0) tmp2(k)(lung2(k)+1:) = 'NLO'//bw
         enddo
      endif
      if (corrections.eq.'oal2') then
         do k=1,ndistr
            tmp(k)(lung(k)+1:) = 'NNLO'//bw
         enddo
         do k=1,ndistr2
            if (i2d.gt.0) tmp2(k)(lung2(k)+1:) = 'NNLO'//bw
         enddo
      endif

      if (corrections.eq.'exp') then
         do k=1,ndistr
            tmp(k)(lung(k)+1:) = 'EXP'//bw
         enddo
         do k=1,ndistr2
            if (i2d.gt.0) tmp2(k)(lung2(k)+1:) = 'NNLO'//bw
         enddo
      endif

      if (corrections.eq.'born') then
         do k=1,ndistr
            tmp(k)(lung(k)+1:) = 'LO'//bw
         enddo
         do k=1,ndistr2
            if (i2d.gt.0) tmp2(k)(lung2(k)+1:) = 'LO'//bw
         enddo
      endif

      if (corrections.eq.'struct') then
         do k=1,ndistr
            tmp(k)(lung(k)+1:) = 'stru_'//bw
         enddo
         do k=1,ndistr2
            if (i2d.gt.0) tmp2(k)(lung2(k)+1:) = 'stru'//bw
         enddo
      endif

      do k = 1,len(programma)
         if (programma(k:k).ne.' ') lpath=k
      enddo

      do k = 1,ndistr
         distname(k)(1:lpath)  = programma
         distname(k)(lpath+1:) = tmp(k)
      enddo

      do k = 1,ndistr2
         if (i2d.gt.0) distname2(k)(1:lpath)  = programma
         if (i2d.gt.0) distname2(k)(lpath+1:) = tmp2(k)
      enddo


      return
      end
*********************************************************************
      subroutine distributions(sd,nc,p3,p4,qph)
#include "shared.F"
      parameter (ibin = 3 * nbin)
      parameter (iperj = ndistr*ibin)
      dimension pin1(0:3),pin2(0:3)
      common/momentainitial/pin1,pin2
      common/reducedtoborn/p1b(0:3),p2b(0:3),iiiref
      character*2 fs
      common/finalstate/fs
      common/idebugging/idebug
      common/distr/distr,s_distr,bmi,bma
      common/distr2/dis2d,edis2d,bmi2,bma2,xaxis2,yaxis2
      common/ionlyfirsttime/ifirst,icount
      data ifirst,icount /1,0/

      common/tmintmax/tmin,tmax,tcur

      common/beamspreadsigmas/esig1,esig2

      common/parameters/ame,ammu,convfac,alpha,pi_i
      common/hparam/ampi

      common/par/pi,radtodeg,ebeam
      common/idarkon/idarkon
      dimension s_distr(ndistr,ibin)
      dimension distr(ndistr,ibin),bmi(ndistr),bma(ndistr),fun(ndistr)

      data ((distr(i,j),i=1,ndistr),j=1,ibin) / iperj * 0.d0/
      data ((s_distr(i,j),i=1,ndistr),j=1,ibin) / iperj * 0.d0/

      integer*8 nc
      integer*8 nccmn
      common/numberofcallsfordistributions/nccmn
      double precision csi(1)

      double precision qa(0:3),qb(0:3),px(0:3),py(0:3)
      double precision pluto(0:3)
      common/massainv/amassainvmin,amassainvmax

      double precision qphtot(40,0:3)
      double precision p1inside(0:3),p2inside(0:3)
      double precision p1rad(0:3),p2rad(0:3)
      integer inside(40),icutgg
      common/cutsggcmn/qphtot,p1inside,p2inside,p1rad,p2rad,
     .     nph,npairs,inside,n_inside,icutgg

#include "strong2020common.F"

      common/todistribute/alford

      ebeam = ecmsnom/2.d0
      nccmn = nc

! default
      if (ifirst.eq.1) then
         pi = 4.d0 * atan(1.d0)
         radtodeg = 180.d0/pi
         bmi(1) = thmine * radtodeg
         bma(1) = thmaxe * radtodeg

         if (fs.eq.'gg') then ! for any fs!
            bmi(1) = 0.d0 * radtodeg - 0.1d0
            bma(1) = pi * radtodeg + 0.1d0
         endif
         bmi(2) = emin  - 0.01d0
         bma(2) = ebeam + 0.01d0!max(esig1,esig2)*2.d0 + 0.000001d0

         bmi(3) =  -.1d0
         bma(3) =  zmax*radtodeg + .01d0
         bmi(4) =  -0.01d0!0.0d0*ebeam
         bma(4) =   bma(2)
         bmi(5) =  -0.1d0!0.d0 * radtodeg - 0.01d0
         bma(5) =   180.1d0! pi * radtodeg + 0.01d0

         bmi(6) =  bmi(1)
         bma(6) =  bma(1)

         bmi(7) =  bmi(2)
         bma(7) =  bma(2)


         bmi(8) = 2.d0*0.511d-3  - 0.00001d0

         bmin = 1.d0
         cccc = cos(pi-zmax)
         bmi(8) = sqrt(2.d0*emin*emin*(1.d0-bmin*bmin*cccc))
         uuu = ame
         if (fs.eq.'mm') uuu = ammu
         if (fs.eq.'pp') uuu = ampi
         if (fs.eq.'pi') uuu = ampi
         bmi(8) = 2.d0*max(emin,uuu)
         bmi(8) = 2.d0*uuu
         bma(8) = ecms

         bmi(9)  =  bmi(4)
         bma(9)  =  bma(4)
         bmi(10) =  bmi(5)
         bma(10) =  bma(5)

c      data mres/
c     .     0.76419892698171199d0, ! rho
c     .     0.78278988240778891d0, ! omega
c     .     1.0198676384707968d0,  ! phi
c     .     3.0969197600196230d0,  ! J/Psi
c     .     3.6861561280010915d0,  ! Psi(2S)
c     .     9.4603013694067908d0,  ! Y
c     .     10.023263041544904d0,  ! Y(2S)
c     .     10.355200529202707d0,  ! Y(3S)
c     .     0.d0,
c     .     0.d0
c     .     /
c      data gres/
c     .     0.149d0,  ! rho
c     .     8.49d-3,  ! omega
c     .     4.266d-3, ! phi
c     .     92.9d-6,  ! J/Psi
c     .     299.d-6,  ! Psi(2S)
c     .     54.02d-6, ! Y
c     .     31.98d-6, ! Y(2S)
c     .     20.32d-6, ! Y(3S)
c     .     0.d0,
c     .     0.d0
c     .     /

         bmi(11) = 3.0969197600196230d0 - 250.d0*92.9d-6
         bma(11) = 3.0969197600196230d0 + 250.d0*92.9d-6
         bmi(12) = 3.6861561280010915d0 - 250.d0*299.d-6
         bma(12) = 3.6861561280010915d0 + 250.d0*299.d-6
         bmi(13) = 9.4603013694067908d0 - 250.d0*54.02d-6
         bma(13) = 9.4603013694067908d0 + 250.d0*54.02d-6

         bmi(14) = cos(thmaxe)
         bma(14) = cos(thmine)

         bmi(15) = cos(thmine)
         bma(15) = 1.d0

         bmi(16) = -1.d0
         bma(16) = cos(thmaxe)

         bmi(17) = -1d0
         bma(17) =  1d0
         bmi(18) = -1d0
         bma(18) =  1d0
         bmi(19) =  0.d0
         if (fs.eq.'mm') then
            bmi(19) = 4.d0*ammu*ammu
         endif
         if (fs.eq.'ee') then
            bmi(19) = 4.d0*ame*ame
         endif
         if (fs.eq.'pp'.or.fs.eq.'pi') then
            bmi(19) = 4.d0*ampi*ampi
         endif
         bma(19) = ecmsnom**2
! filling distr. components when entering the first time
         do k = 1,ndistr
            d = (bma(k) - bmi(k))/nbin
            y = bmi(k)
            do i = 1,nbin
               y = y + d
               distr(k,i*3-2)   = y - d
               s_distr(k,i*3-2) = y - d
            enddo
         enddo

****  2d distributions
         if (i2d.gt.0) then
            bmi2(1,1) = bmi(19)
            bmi2(1,2) = -1.d0 !bmi(14)
            bma2(1,1) = bma(19)
            bma2(1,2) = 1.d0 !bma(14)

            dis2d  = 0.d0
            edis2d = 0.d0
            do k = 1,ndistr2
               dx = (bma2(k,1) - bmi2(k,1))/nbin2
               dy = (bma2(k,2) - bmi2(k,2))/nbin2
               x = bmi2(k,1)
               y = bmi2(k,2)
               do i = 1,nbin2
                  x = x+dx
                  y = y+dy
                  xaxis2(k,i) = x - dx
                  yaxis2(k,i) = y - dy
               enddo
            enddo
         endif
**********************************

         ifirst = 0
      endif
******************
      if (fs.ne.'gg') then
         px = p3
         py = p4
         call duefotonipiuenergetici(qph,qa,qb)
c         qa = qph(1,:)
c         qb = qph(2,:)
      else
         px = p1inside
         py = p2inside
         qa = p1rad
         qb = p2rad
      endif
*** electron angle
      pe      = sqrt(px(1)**2+px(2)**2+px(3)**2)
      fun(17) = px(3)/pe
      fun(1)  = radtodeg*acos(fun(17))
**** electron energy
      fun(2) = px(0)
*** positron angle
      pp      = sqrt(py(1)**2+py(2)**2+py(3)**2)
      fun(18) = py(3)/pp
      fun(6)  = radtodeg*acos(fun(18))
**** positron energy
      fun(7) = py(0)

**** acollinearity
      fun(3)  = abs(180.d0 - fun(1)- fun(6))

***   m.e. photon energy
      fun(4)  = bma(4) + 1.d0
      fun(5)  = bma(5) + 1.d0
      fun(9)  = bma(9) + 1.d0
      fun(10) = bma(10) + 1.d0

      soglia = 0.d0
      if (qa(0).ge.soglia) then
         fun(4)  = qa(0)
         fun(5) = radtodeg*acos(qa(3)/qa(0))
      endif
      if (qb(0).ge.soglia) then
         fun(9)  = qb(0)
         fun(10) = radtodeg*acos(qb(3)/qb(0))
      endif

*** CFR variables
      call new_boost(px+py,px,qa,1)
c      call new_boost(px+py,py,qb,1)
      fun(15) = qa(3)/sqrt(tridot(qa,qa))
      fun(14) = fun(15)
      fun(16) = fun(15)
*****

      fun(19) = abs(dot(px+py,px+py))
      fun(8)  = sqrt(fun(19))
      fun(11) = fun(8)
      fun(12) = fun(8)
      fun(13) = fun(8)
***********************************************************************
      DO k = 1,ndistr
         bmax = bma(k)
         d = (bma(k) - bmi(k))/nbin
         x = bmi(k)
** NEW ***
         if (fun(k).le.bma(k).and.fun(k).ge.bmi(k)) then
            i  = int((fun(k) - bmi(k))/d)
c            if ((i+1).gt.nbin) print*,'WARNING!!',
c     .                            k,i+1,d,fun(k),bmi(k),bma(k)
            i  = min(i+1,nbin)*3

            distr(k,i-1) = distr(k,i-1) + sd
            distr(k,i)   = distr(k,i)   + sd**2
            su  = distr(k,i-1)
            su2 = distr(k,i)
            s_distr(k,i-1) = su / nc
            argument = abs((su2/nc-s_distr(k,i-1)**2)/nc)
            s_distr(k,i) = sqrt(argument)
         endif
*     * NEW ***
      ENDDO
***************************************
      IF (i2d.gt.0) then
*** 2D distributions
         fun2(1,1) = fun(19)
         fun2(1,2) = fun(14)
         DO k = 1,ndistr2
            dx = (bma2(k,1) - bmi2(k,1))/nbin2
            dy = (bma2(k,2) - bmi2(k,2))/nbin2
            if ((fun2(k,1).le.bma2(k,1).and.fun2(k,1).ge.bmi2(k,1)).and.
     .         (fun2(k,2).le.bma2(k,2).and.fun2(k,2).ge.bmi2(k,2))) then
            i = int((fun2(k,1) - bmi2(k,1))/dx)
            j = int((fun2(k,2) - bmi2(k,2))/dy)
            i = min(i+1,nbin2)  ! to avoid roundoffs (?)
            j = min(j+1,nbin2)  ! to avoid roundoffs (?)

            dis2d(k,i,j)  = dis2d(k,i,j)  + sd
            edis2d(k,i,j) = edis2d(k,i,j) + sd*sd
         endif
      ENDDO
      ENDIF

      return
      end
***
      subroutine writedistributions
#include "shared.F"
      common/ifirstwritedistributions/ifirst
      integer*8 nccmn
      common/numberofcallsfordistributions/nccmn

      parameter (ibin = 3 * nbin)
      parameter (iperj = ndistr*ibin)

      common/filenames/distname,distname2
      common/distr/distr,s_distr,bmi,bma
      common/distr2/dis2d,edis2d,bmi2,bma2,xaxis2,yaxis2

      dimension s_distr(ndistr,ibin),bmi(ndistr),bma(ndistr)
      dimension distr(ndistr,ibin)

      character*100 distname(ndistr),distname2(ndistr2)

      data ifirst /1/

      do k = 1,ndistr
         open(12,file=distname(k),status='unknown')
         do i=1, nbin
            sezd  = s_distr(k,i*3-1)/(bma(k)-bmi(k))*nbin
            esezd = s_distr(k,i*3)  /(bma(k)-bmi(k))*nbin
            write(12,*)s_distr(k,i*3-2),sezd,esezd
         enddo
         close(12)
      enddo

*** 2D distributions
      IF (i2d.gt.0) THEN
      do k = 1,ndistr2
         open(12,file=distname2(k),status='unknown')
         do i = 1, nbin2
         do j = 1, nbin2
            av   = 0.d0
            stdd = 0.d0
            stddroot = 0.d0
            if (abs(dis2d(k,i,j)).gt.0.d0) then
               av  = dis2d(k,i,j)/nccmn
               av2 = edis2d(k,i,j)/nccmn
               stdd     = sqrt(abs(av2 -av*av)/nccmn)
            endif
            sezd = av/(bma2(k,1)-bmi2(k,1))*nbin2
            sezd = sezd/(bma2(k,2)-bmi2(k,2))*nbin2
            esezd = stdd/(bma2(k,1)-bmi2(k,1))*nbin2
            esezd = esezd/(bma2(k,2)-bmi2(k,2))*nbin2
            write(12,*)xaxis2(k,i),yaxis2(k,j),sezd,esezd
         enddo
         write(12,*)' '
         enddo
         close(12)
      enddo
      ENDIF
****************************************************

      ifirst = 0
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mostenergeticphoton(q,q1)
      implicit double precision (a-h,o-z)
      dimension q(40,0:3),q1(0:3)
*  LEADING ENERGETIC PHOTON IS EXTRACTED
      q1(0) = 0.d0
      q1(1) = 0.d0
      q1(2) = 0.d0
      q1(3) = 0.d0

ccc      if (q(1,0).lt.1.d-11) return
C for backward compatibility with old ps-based BABAYAGA
      if (q(1,0).lt.1.d-12.and.q(11,0).lt.1.d-12.and.
     .     q(21,0).lt.1.d-12.and.q(31,0).lt.1.d-12) return

      ENPHOT = Q(1,0)
      ENHARD = ENPHOT
      J=1
*
      DO I = 1,40
         if (q(i,0).gt.0.d0) then
            ENPHOT = Q(I,0)
            IF (ENHARD.GE.ENPHOT) THEN
               ENHARD = ENHARD
            ELSE
               J = I
               ENHARD = ENPHOT
            ENDIF
         endif
      ENDDO
*
      ENHARDL = ENHARD
      JL = J
*
      if (j.gt.0) then
         DO I = 0,3
            Q1(I) = Q(J,I)
         ENDDO
      endif
      return
      end
*****************************************
      function get_dphi(p1,p2)
      implicit double precision (a-h,o-z)
      double precision p1(0:3),p2(0:3)
! from ALPGEN
      dphi = (p1(1)*p2(1)+p1(2)*p2(2))/
     .     sqrt(p1(1)**2+p1(2)**2)/sqrt(p2(1)**2+p2(2)**2)
      if (abs(dphi).gt.1.d0) then
         dphi = dphi/abs(dphi)
      endif
      get_dphi = acos(dphi)
      return
      end
***********************

      function getphi_here(ppp)
      implicit double precision (a-h,o-z)
      dimension ppp(0:3)
      parameter (pi = 3.141592653589793238462643383279502884197169399d0)
      pm = sqrt(ppp(1)**2+ppp(2)**2+ppp(3)**2)
      c  = ppp(3)/pm
      s  = sqrt(1.d0-c**2)
      if (s.eq.0.d0) then
         getphi_here = 0.d0
         return
      else
         arg = ppp(1)/pm/s
*  avoiding numerical problems......
         if (abs(arg).ge.1.d0) then
            iarg = arg
            arg  = iarg
         endif
         if (ppp(2).ge.0.d0) getphi_here = acos(arg)
         if (ppp(2).lt.0.d0) getphi_here = 2.d0*pi-acos(arg)
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccc
      subroutine duefotonipiuenergetici(q,q1,q2)
      implicit real*8 (a-h,o-z)
      dimension q(40,0:3),q1(0:3),q2(0:3),ql(40,0:3)
*
*  LEADING ENERGETIC PHOTON IS EXTRACTED
*
      q1 = 0.d0
      q2 = 0.d0

      n = 0
      do k = 1,40
         if (q(k,0).gt.0.d0) then
            n = n+1
            ql(n,:) = q(k,:)
         endif
      enddo
      if (n.eq.0) return

      ENPHOT = QL(1,0)
      ENHARD = ENPHOT
      J=1
*
      DO I = 1,N
         ENPHOT = QL(I,0)
         IF (enphot.ge.enhard) THEN
            J = I
            ENHARD = ENPHOT
         ENDIF
      ENDDO
*
      Q1 = QL(J,:)
      ENHARDL = ENHARD
      JL = J
*
      j1 = j
*
*     NEXT-TO-LEADING ENERGETIC PHOTON IS EXTRACTED
*
      ENHARD = 0.D0
      J = 0
*
      DO I = 1,N
         IF (I.NE.JL) THEN
            ENPHOT = QL(I,0)
            IF (ENPHOT.GE.ENHARD) THEN
               J = I
               ENHARD = ENPHOT
            ENDIF
         ENDIF
      ENDDO
c
      if (j.gt.0) then
         Q2 = QL(J,:)
      endif
      j2=j

      return
      end
************************************************
      subroutine duefotonipiuenergetici_bck(q,q1,q2)
      implicit real*8 (a-h,o-z)
      dimension q(40,0:3),q1(0:3),q2(0:3)
*
*  LEADING ENERGETIC PHOTON IS EXTRACTED
*
        do i = 0,3
           q1(i) = 0.d0
           q2(i) = 0.d0
        enddo

        sum = 0.d0
	do k=1,40
          sum = sum + q(k,0)
        enddo

        if (sum.lt.1.d-11) return

        ENPHOT = Q(1,0)
        ENHARD = ENPHOT
        J=1
*
        DO I = 1,40
           ENPHOT = Q(I,0)
           IF (ENHARD.GE.ENPHOT) THEN
              ENHARD = ENHARD
           ELSE
              J = I
              ENHARD = ENPHOT
           ENDIF
        ENDDO
*
        ENHARDL = ENHARD
        JL = J
*
        if (j.gt.0) then
           DO I = 0,3
              Q1(I) = Q(J,I)
           ENDDO
        endif
        j1 = j
*
*  NEXT-TO-LEADING ENERGETIC PHOTON IS EXTRACTED
*
      ENHARD = 0.D0
      J = 0
*
      DO I = 1,40
         IF (I.NE.JL) THEN
*
            ENPHOT = Q(I,0)
            IF (ENHARD.GE.ENPHOT) THEN
               ENHARD = ENHARD
            ELSE
               J = I
               ENHARD = ENPHOT
            ENDIF
*
         ENDIF
      ENDDO
*
      ENHARD = ENHARD
*
      if (j.gt.0) then
         DO I = 0,3
            Q2(I) = Q(J,I)
         ENDDO
      endif
      j2=j

      return
      end
***************************************************
      subroutine getmasstrack(ng,sqs,p1,p2,m2,m2trk)
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3)
!     should be used only if the number of photons is >= 2!!

      m2trk = m2
      if (ng.lt.2) return

      p1m = sqrt(tridot(p1,p1))
      p2m = sqrt(tridot(p2,p2))

! by Stefan Mueller, from phokhara10
      squared = 0.5d0*((sqs-sqrt(tridot(p1+p2,p1+p2)))**2-p1m**2-p2m**2)
      m2trk = ( squared**2 - p1m**2*p2m**2 ) /
     .     ( 2.d0*squared + p1m**2 + p2m**2 )
      m2trk = abs(m2trk)
c      print*,m2trk,fforgetmasstrack(sqs,m2trk,p1,p2)
***   on gitlab
c      a = (sqs - sqrt(tridot(p1+p2,p1+p2)))**2
c      uu = 0.25d0*a - 0.5d0*(p1m**2+p2m**2)+0.25d0/a*(p1m**2-p2m**2)**2
c      print*,(uu-m2trk)/m2trk

      return
      end
***************************************************
      function fforgetmasstrack(sqs,m2,p1,p2)
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3)
      p1m = sqrt(tridot(p1,p1))
      p2m = sqrt(tridot(p2,p2))
      E1 = sqrt(p1m*p1m + m2)
      E2 = sqrt(p2m*p2m + m2)
c      fforgetmasstrack = sqs*sqs + 2.d0*m2 +
c     .     2d0*E1*E2 - 2d0*tridot(p1,p2)
c     .    -2.d0*sqs*(E1+E2)
c      f1 = fforgetmasstrack
      fforgetmasstrack = (sqs - E1 - E2)**2 - tridot(p1+p2,p1+p2)
      return
      end

