      subroutine phasespaceStrong2020(p1,p2,p3,p4,qph,nphot,m1,m2,esoft,
     .     cth,w,phsp,ie)
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension ap3(0:3),ap4(0:3),pcollorig(0:3)
      dimension p4guess(0:3),p3guess(0:3)
      dimension dir3(0:3),dir4(0:3),ptmp(0:3),vers(0:3)
      double precision kone(0:3),pcoll(0:3),ptld(0:3),koner(0:3)
      double precision lambda
      logical lcond
      character*2 fs
      common/finalstate/fs
      integer emissionpattern(nphot)
      parameter (pi = 3.1415926535897932384626433832795029d0)
      parameter (imaxph = 40)
      double precision qph(imaxph,0:3),oph(imaxph),wten(nphot)
      double precision csi(1),csi2(2),rnddummy
      common/for_debug/pcoll,denphsp,dir3,dir4,br,bq,onm
      common/various/beta,betafs
      common/channelref/iref
      common/iwide/iwider,ntot
      common/regulators/ich,idummy,regulator1,regulator2,wcoll3,wcoll4
      common/radpattern/nph(4)
      common/idebugging/idebug
      common/idarkon/idarkon
#include "strong2020common.F"
      w      = 1.d0
      phsp   = 1.d0
      npart = nphot + 2
      nvar  = 3*nphot + 2

      ie     = 0
      duepigrechi = (2.d0*pi)**(-nvar)

      iwider = 0
      ntot   = nphot

      do j = 0,3
         do k = 1,imaxph
            qph(k,j) = 0.d0
         enddo
         nph(j+1) = 0
      enddo

      if (p1(0)+p2(0).lt.(m1+m2+nphot*esoft)) then
         phsp = 0.d0
         w    = 0.d0
         ie = 1
         return
      endif
*****************
!! I put this outside in order to raise efficiency
!      call get_cos_fer(rnd(1),c,wcos)
!      cth = c
!      w   = w * wcos
      c = cth
*
      if (nphot.eq.0) then
         mass = p1(0)+p2(0)
         sqla = dsqrt(lambda(mass**2,m1**2,m2**2))
         phsp = phsp * sqla/8.d0/mass**2
         s = dsqrt(1.d0 - c**2)
         call getrnd(csi,1)
         call getphi(csi(1),phi,wph)
         w = w*wph
         sphi = sin(phi)
         cphi = cos(phi)
         p1mod = sqla/2.d0/mass
         p3(0) = dsqrt(p1mod**2 + m1**2)
         p3(1) = p1mod * sphi * s
         p3(2) = p1mod * cphi * s
         p3(3) = p1mod * c
         p4(0) = mass - p3(0)
         p4(1) = - p3(1)
         p4(2) = - p3(2)
         p4(3) = - p3(3)
         phsp = phsp * duepigrechi
         return
      endif
*******************
******************
      lcond = fs.eq.'ee'.or.fs.eq.'mm'.or.fs.eq.'mr'
     +     .or.fs.eq.'pp'.or.fs.eq.'pi'.or.fs.eq.'er'
     +     .or.fs.eq.'pr'

      angle = 30.d0 *pi/180.d0
      angle = 0.d0 *pi/180.d0

      cpmax = cos(angle)
      cpmin = cos(pi - angle)
c the following for ng >= 1


      ncharged = 4
      if (fs.eq.'gg'.or.fs.eq.'pi') ncharged = 2

      ir = 1

! cos(th) of electron already generated...
      sth  = dsqrt(1.d0 - cth**2)
      call getrnd(csi,1)
      call getphi(csi(1),phi,wph)
      w    = w*wph
      sphi = dsin(phi)
      cphi = dcos(phi)
      ir   = ir + 2

      call get_pattern(ncharged,nphot,emissionpattern)
      do k = 1,nphot
         nph(emissionpattern(k)) = nph(emissionpattern(k)) + 1
      enddo
      n1 = nph(1)
      n2 = nph(2)
      n3 = nph(3)
      n4 = nph(4)
      nis = n1 + n2
      nfs = n3 + n4

      ei = p1(0)+p2(0)

      dir3(0) = ei/2d0
      dir3(1) = ei/2d0*betafs*cphi*sth
      dir3(2) = ei/2d0*betafs*sphi*sth
      dir3(3) = ei/2d0*betafs*cth
      dir4(0) = dir3(0)
      pcoll(0) = p1(0)+p2(0)
      ptld(0)  = p1(0)+p2(0)
      vers(0)  = 0.d0
      do k = 1,3
         pcoll(k) = p1(k)+p2(k)
         ptld(k)  = p1(k)+p2(k)
         dir4(k)  = -dir3(k)
         vers(k)  = 0.d0
      enddo

      call getrnd(csi,1)
      iref = 3
      if (csi(1).gt.0.5) iref = 4
      wref = 1.d0
      if (nfs.gt.0) then
         if (n4.eq.0) iref = 3
         if (n3.eq.0) iref = 4
      endif
      w = w*wref
      if (fs.eq.'gg'.or.fs.eq.'pi') iref = 3
****
** PHOTONS GENERATION...
      pflatc  = 0.d0
      pflatc3 = 0.d0
      pflatc4 = 0.d0
      pflate  = 0.d0
      expgg = 0.92d0
      ibich = 1
      emaxgamma = (p1(0)+p2(0))/2.d0
      do k = 1,nphot
         wen = 1.d0
         en  = esoft
         if (fs.eq.'gg') then
            call photon_energy(ibich,pflate,expgg,p1(0),
     .           en,wen)
         endif
         wten(k) = wen
         if (fs.eq.'ee'.or.fs.eq.'mm'.or.fs.eq.'mr'
     .   .or.fs.eq.'pp'.or.fs.eq.'pi'.or.fs.eq.'er'
     >	 .or.fs.eq.'pr')  then
         w = w * wen
         endif
         oph(k) = en
         ir = ir + 1
       enddo

      kgiag = -10*nphot ! this recovers the default

#if defined KLOEI || defined B || defined BES3
      call getrnd(csi2,2)
      ir = ir + 2
      kgiag = 1 + nphot*csi2(1)
      omax = emaxgamma
      omin = esoft
#ifdef KLOEI
      omin = egammamin
#endif
#ifdef BES3
      omin = egammamin1
#endif
#ifdef B
      omin = egammamin
#endif
      an  = log(omax/omin)
      en  = omin*exp(an*csi2(2))
      wen = an*en
      oph(kgiag) = en
      w          = w * wen
!      wten(kgiag) = wen ! not used for ee or mm
#endif

** initial state photons
      if (nis.gt.0) then
         if (n1.gt.0) then
            do k = 1,n1
               iaccendi = 1
               call getrnd(csi2,2)
               pfl = 0.d0
#if defined KLOEI || defined KLOEII || defined B || defined BES3
               pfl = 0.05d0
cc               if (k.eq.kgiag) pfl = 0.2d0
#endif
c               if (k.eq.kgiag) pfl = 0.2d0
c               if (k.ne.kgiag) then
                  call collinear(pfl,csi2(1),csi2(2),p1,vers,wverg)
c               else
c                  call collinearstr2020(csi2(1),csi2(2),p1,vers,wverg)
c               endif
               iaccendi = 0
               ir = ir + 2
               w = w * wverg
****
               if (lcond.and.k.ne.kgiag) then
                  cp = vers(3)/vers(0)
                  if (cp.gt.cpmax.or.cp.lt.cpmin) then
                     call photon_energy_new(ibich,pflate,expgg,
     .                    emaxgamma,esoft,ptld,vers,en,wen)
                  else
                     call photon_energy(ibich,pflate,expgg,
     .                    emaxgamma,en,wen)
                  endif
                  w = w * wen
                  oph(k) = en
               endif
***
               qph(k,0) = oph(k)
               qph(k,1) = oph(k)*vers(1)
               qph(k,2) = oph(k)*vers(2)
               qph(k,3) = oph(k)*vers(3)

               pcoll(0) = pcoll(0) - qph(k,0)
               pcoll(1) = pcoll(1) - qph(k,1)
               pcoll(2) = pcoll(2) - qph(k,2)
               pcoll(3) = pcoll(3) - qph(k,3)
            enddo
         endif
         if (n2.gt.0) then
            do k = 1,n2
               iaccendi = 2
               call getrnd(csi2,2)
               pfl = 0.d0
#if defined KLOEI || defined KLOEII || defined B || defined BES3
               pfl = 0.05d0
ccc               if ((k+n1).eq.kgiag) pfl = 0.2d0
#endif
c               if ((k+n1).eq.kgiag) pfl = 0.2d0
c               if ((k+n1).ne.kgiag) then
                  call collinear(pfl,csi2(1),csi2(2),p2,vers,wverg)
c               else
c                  call collinearstr2020(csi2(1),csi2(2),p2,vers,wverg)
c               endif
               iaccendi  = 0
               ir = ir + 2
               w = w * wverg
****
               if (lcond.and.(k+n1).ne.kgiag) then
                  cp = vers(3)/vers(0)
                  if (cp.gt.cpmax.or.cp.lt.cpmin) then
                     call photon_energy_new(ibich,pflate,expgg,
     .                    emaxgamma,esoft,ptld,vers,en,wen)
                  else
                     call photon_energy(ibich,pflate,expgg,
     .                    emaxgamma,en,wen)
                  endif
                  w = w * wen
                  oph(k+n1) = en
               endif
***
               qph(k+n1,0) = oph(k+n1)
               qph(k+n1,1) = oph(k+n1)*vers(1)
               qph(k+n1,2) = oph(k+n1)*vers(2)
               qph(k+n1,3) = oph(k+n1)*vers(3)

               pcoll(0) = pcoll(0) - qph(k+n1,0)
               pcoll(1) = pcoll(1) - qph(k+n1,1)
               pcoll(2) = pcoll(2) - qph(k+n1,2)
               pcoll(3) = pcoll(3) - qph(k+n1,3)
            enddo
         endif
      endif

      if (lcond.and.
     .   (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0)) then
         ie = 1
         phsp = 0.d0
         w  = 0.d0
         return
      endif

** final state photons
c      en3guess = dir3(0)
c      en4guess = dir4(0)
      wcoll3 = 1.d0
      wcoll4 = 1.d0
      if (nfs.gt.0) then
         if (iref.eq.3) then
            if (n3.gt.0) then
               do k = 1,n3
                  iaccendi = 3
                  call getrnd(csi2,2)
                  call collinear(pflatc3,csi2(1),csi2(2),dir3,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond.and.(k+n1+n2).ne.kgiag) then
                     cp = vers(3)/vers(0)
                     if (cp.gt.cpmax.or.cp.lt.cpmin) then
                        call photon_energy_new(ibich,pflate,expgg,
     .                       emaxgamma,esoft,ptld,vers,en,wen)
                     else
                        call photon_energy(ibich,pflate,expgg,
     .                       emaxgamma,en,wen)
                     endif
                     w = w * wen
                     oph(k+n1+n2) = en
                  endif
***
                  wcoll3 = wcoll3 * wverg
                  qph(k+n1+n2,0) = oph(k+n1+n2)
                  qph(k+n1+n2,1) = oph(k+n1+n2)*vers(1)
                  qph(k+n1+n2,2) = oph(k+n1+n2)*vers(2)
                  qph(k+n1+n2,3) = oph(k+n1+n2)*vers(3)

                  pcoll(0) = pcoll(0) - qph(k+n1+n2,0)
                  pcoll(1) = pcoll(1) - qph(k+n1+n2,1)
                  pcoll(2) = pcoll(2) - qph(k+n1+n2,2)
                  pcoll(3) = pcoll(3) - qph(k+n1+n2,3)
c                  en3guess = en3guess - oph(k+n1+n2)
               enddo
            endif
            if (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0) then
               ie = 1
c               print*,'1 -- '
               phsp = 0.d0
               w  = 0.d0
               return
            endif

            do k = 0,3
               pcoll(k) = pcoll(k) - dir3(k)
               pcollorig(k) = pcoll(k)
            enddo

            pcmodmu = 1.d0/sqrt(tridot(pcoll,pcoll))
            pcoll(0) = dir4(0)
            bcoll    = sqrt(1.d0 - m1**2/pcoll(0)**2)
            pcoll(1) = pcoll(1)*pcmodmu*bcoll*pcoll(0)
            pcoll(2) = pcoll(2)*pcmodmu*bcoll*pcoll(0)
            pcoll(3) = pcoll(3)*pcmodmu*bcoll*pcoll(0)

            if (n4.gt.0) then
               do k = 1,n4
                  iaccendi = 4
                  call getrnd(csi2,2)
                  call collinear(pflatc4,csi2(1),csi2(2),pcoll,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond.and.(k+n1+n2+n3).ne.kgiag) then
                     cp = vers(3)/vers(0)
                     if (cp.gt.cpmax.or.cp.lt.cpmin) then
                        call photon_energy_new(ibich,pflate,expgg,
     .                       emaxgamma,esoft,ptld,vers,en,wen)
                     else
                        call photon_energy(ibich,pflate,expgg,
     .                       emaxgamma,en,wen)
                     endif
                     w = w * wen
                     oph(k+n1+n2+n3) = en
                  endif
***
                  wcoll4 = wcoll4 * wverg
                  qph(k+n1+n2+n3,0) = oph(k+n1+n2+n3)
                  qph(k+n1+n2+n3,1) = oph(k+n1+n2+n3)*vers(1)
                  qph(k+n1+n2+n3,2) = oph(k+n1+n2+n3)*vers(2)
                  qph(k+n1+n2+n3,3) = oph(k+n1+n2+n3)*vers(3)
               enddo
            endif
         else                   ! if iref = 4
            if (n4.gt.0) then
               do k = 1,n4
                  iaccendi = 4
                  call getrnd(csi2,2)
                  call collinear(pflatc4,csi2(1),csi2(2),dir4,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond.and.(k+n1+n2+n3).ne.kgiag) then
                     cp = vers(3)/vers(0)
                     if (cp.gt.cpmax.or.cp.lt.cpmin) then
                        call photon_energy_new(ibich,pflate,expgg,
     .                       emaxgamma,esoft,ptld,vers,en,wen)
                     else
                        call photon_energy(ibich,pflate,expgg,
     .                       emaxgamma,en,wen)
                     endif
                     w = w * wen
                     oph(k+n1+n2+n3) = en
                  endif
***
                  wcoll4 = wcoll4 * wverg
                  qph(k+n1+n2+n3,0) = oph(k+n1+n2+n3)
                  qph(k+n1+n2+n3,1) = oph(k+n1+n2+n3)*vers(1)
                  qph(k+n1+n2+n3,2) = oph(k+n1+n2+n3)*vers(2)
                  qph(k+n1+n2+n3,3) = oph(k+n1+n2+n3)*vers(3)

                  pcoll(0) = pcoll(0) - qph(k+n1+n2+n3,0)
                  pcoll(1) = pcoll(1) - qph(k+n1+n2+n3,1)
                  pcoll(2) = pcoll(2) - qph(k+n1+n2+n3,2)
                  pcoll(3) = pcoll(3) - qph(k+n1+n2+n3,3)
c                  en4guess = en4guess - oph(k+n1+n2+n3)
               enddo
            endif
            if (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0) then
               ie = 1
c               print*,'2 -- '
               phsp = 0.d0
               w  = 0.d0
               return
            endif

            do k = 0,3
               pcoll(k) = pcoll(k) - dir4(k)
               pcollorig(k) = pcoll(k)
            enddo

            pcmodmu = 1.d0/sqrt(tridot(pcoll,pcoll))
            pcoll(0) = dir3(0)
            bcoll    = sqrt(1.d0 - m2**2/pcoll(0)**2)
            pcoll(1) = pcoll(1)*pcmodmu*bcoll*pcoll(0)
            pcoll(2) = pcoll(2)*pcmodmu*bcoll*pcoll(0)
            pcoll(3) = pcoll(3)*pcmodmu*bcoll*pcoll(0)

            if (n3.gt.0) then
               do k = 1,n3
                  iaccendi = 3
                  call getrnd(csi2,2)
                  call collinear(pflatc3,csi2(1),csi2(2),pcoll,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond.and.(k+n1+n2).ne.kgiag) then
                     cp = vers(3)/vers(0)
                     if (cp.gt.cpmax.or.cp.lt.cpmin) then
                        call photon_energy_new(ibich,pflate,expgg,
     .                       emaxgamma,esoft,ptld,vers,en,wen)
                     else
                        call photon_energy(ibich,pflate,expgg,
     .                       emaxgamma,en,wen)
                     endif
                     w = w * wen
                     oph(k+n1+n2) = en
                  endif
***
                  wcoll3 = wcoll3 * wverg
                  qph(k+n1+n2,0) = oph(k+n1+n2)
                  qph(k+n1+n2,1) = oph(k+n1+n2)*vers(1)
                  qph(k+n1+n2,2) = oph(k+n1+n2)*vers(2)
                  qph(k+n1+n2,3) = oph(k+n1+n2)*vers(3)
               enddo
            endif
         endif
      endif
***********************************
      do j = 0,3
         kone(j) = 0.d0
         do k = 1,nphot
            kone(j) = kone(j) + qph(k,j)
         enddo
         ptmp(j) = p1(j) + p2(j) - kone(j)
      enddo

***** TEST
      if (fs.eq.'gg') then
         koner = kone-qph(nphot,:) ! kone ridotto dell'ultimo fotone

         pk2 = dot(koner,koner)
         br  = koner(0) - pk2/2.d0/ei
         bq  = koner(0) - tridot(koner,vers)

         onm = ei*0.5 - br
         onm = onm/(1.d0 - bq/ei)

         a  = ei*0.5d0 - koner(0) + pk2/ei*0.5d0
         b  = 1.d0  - dot(koner,vers)/ei
         c  = ei - koner(0) + tridot(koner,dir3)/betafs/dir3(0)
         d  = 1.d0 - tridot(vers,dir3)/betafs/dir3(0)

         ep = (ei*a - esoft*c)/(ei*b-esoft*d)

c         ep = ei*0.5d0

         if (.not.(ep.lt.esoft.or.ep.ge.a/b.or.ep.eq.p1(0))) then
c         if (.FALSE.) then
! qui rimpiazza l'ultimo fotone generato, se non si verifica la condizione
! la condizione non dovrebbe presentarsi quasi mai, c'era beta invece di betafs sopra!!
            call photonenergygg(ep,a,b,esoft,en,wen)

            w = w * wen

            if (isnan(en).or.isnan(wen)) then
               print*,'!!!!!!!',ep,en,wen
               print*,'a b c d dir3',a,b,c,d,dir3
            endif
            wwen = 1.d0
            if (nphot.gt.1) then
               do k = 1,nphot-1
                  wwen = wwen * wten(k)
               enddo
            endif
            w = w * wwen

            qph(nphot,0) = en
            qph(nphot,1) = en*vers(1)
            qph(nphot,2) = en*vers(2)
            qph(nphot,3) = en*vers(3)
            do j = 0,3
               kone(j) = koner(j)+qph(nphot,j)
               ptmp(j) = p1(j) + p2(j) - kone(j)
            enddo
         else
!     altrimenti si tiene quello gia' generato e si va avanti
c            print*,'OOOOPS in phasespace ep',ep,a/b,p1(0),
c     .           tridot(vers,dir3)/betafs/dir3(0),d,
c     .           ep-a/b
            do k = 1,nphot
               w = w * wten(k)
            enddo
         endif
      endif
***** TEST

      if (dot(ptmp,ptmp).lt.(m1+m2)**2.or.ptmp(0).lt.0.d0) then
         ie = 1
c         print*,'4 -- '
         phsp = 0.d0
         w  = 0.d0
         return
      endif
      amkone2 = dot(kone,kone)
***
      if (iref.eq.3) then
         vx  = sth*cphi
         vy  = sth*sphi
         vz  = cth
      else
         vx  = -sth*cphi
         vy  = -sth*sphi
         vz  = -cth
      endif
      o   = kone(0)
      ox  = kone(1)
      oy  = kone(2)
      oz  = kone(3)
      odv = ox*vx + oy*vy + oz*vz
      abig = ei**2 - 2.d0 * ei * o + amkone2 + m1**2 - m2**2
      bbig = -2.d0 * (ei - o)
      cbig = -2.d0 * odv
      arg  = abig**2*bbig**2-m1**2*bbig**2*(bbig**2-cbig**2)
      if (arg.lt.0.d0) then
         phsp = 0.d0
         w    = 0.d0
c         print*,'5 -- '
         ie = 2
         return
      endif
      p1mod_1 = abig*cbig + dsqrt(arg)
      p1mod_1 = p1mod_1 / (bbig**2 - cbig**2)
      p1mod_2 = abig*cbig - dsqrt(arg)
      p1mod_2 = p1mod_2 / (bbig**2 - cbig**2)
      p1jac = 1.d0
      if (p1mod_1.gt.0.d0.and.p1mod_2.gt.0.d0) then
         call getrnd(csi,1)
         if (csi(1).lt.0.5) p1mod = p1mod_1
         if (csi(1).ge.0.5) p1mod = p1mod_2
         p1jac = 2.d0
      endif
      if (p1mod_1.gt.0.d0.and.p1mod_2.le.0.d0) then
         p1mod = p1mod_1
      endif
      if (p1mod_1.le.0.d0.and.p1mod_2.gt.0.d0) then
         p1mod = p1mod_2
      endif
      if (p1mod_1.le.0.d0.and.p1mod_2.le.0.d0) then ! bestiale: con gg ci vuole .le.!!
         phsp = 0.d0
         w    = 0.d0
c         print*,'6 -- '
         ie   = 3
         return
      endif
      w = w * p1jac
      p3(0) = dsqrt(p1mod**2 + m1**2)
      p3(1) = p1mod * vx
      p3(2) = p1mod * vy
      p3(3) = p1mod * vz
      p4(0) =  ei - p3(0) - o
      p4(1) = -ox - p3(1)
      p4(2) = -oy - p3(2)
      p4(3) = -oz - p3(3)
***
      if (fs.eq.'gg'.and.(p3(0).lt.esoft.or.p4(0).lt.esoft)) then ! needed!
         phsp = 0.d0
         w    = 0.d0
         ie = 1
         return
      endif
***
!      phasespace!
      prodomega = 1.d0
      do k = 1,nphot
         prodomega = prodomega * qph(k,0)
      enddo
      den = dabs( p4(0)*p1mod + p3(0)*(p1mod + odv) )
      denphsp = den
      phsp = prodomega * p1mod**2 / denphsp
      phsp = phsp*duepigrechi/2.d0**(nphot+2)
      if (iref.eq.4) call exchange_mom(p3,p4)

***********************
! "regulator"... I try only on FINAL STATE (active only for ee and mm...)
      flatconst = 0.d0 !!!
      den       = 0.d0
      if (n3.gt.0.and.n4.gt.0) then
         prodp3 = 1.d0
         prodp4 = 1.d0
         do k = nis+1,nphot
            do i = 0,3
               ptmp(i) = qph(k,i)
            enddo
            p3q = dot(p3,ptmp)
            p4q = dot(p4,ptmp)
            if (k.le.(n1+n2+n3)) then
               prodp3=prodp3*
     .              (1.d0/p3q*p3(0)*ptmp(0))
            endif
            if (k.gt.(n1+n2+n3)) then
               prodp4=prodp4*
     .              (1.d0/p4q*p4(0)*ptmp(0))
            endif
         enddo
         if (iref.eq.3) anum = prodp3
         if (iref.eq.4) anum = prodp4
         regulator1 = anum/(prodp3+prodp4) * 2.d0
         w = w * regulator1
      endif
*******************
* second "regulator"
      den        = 1.d0
      anum       = 1.d0
      iswitch    = 1
      regulator2 = 1.d0
      if (fs.eq.'gg'.or.fs.eq.'pi') iswitch = 0
      if (fs.ne.'gg') then
         if (n1.gt.0) then
            do k = 1,n1
               do i = 0,3
                  ptmp(i) = qph(k,i)
               enddo
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p1q*p1(0)*ptmp(0))
            enddo
         endif
         if (n2.gt.0) then
            do k = n1+1,n1+n2
               do i = 0,3
                  ptmp(i) = qph(k,i)
               enddo
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p2q*p2(0)*ptmp(0))
            enddo
         endif
         if (n3.gt.0) then
            do k = n1+n2+1,n1+n2+n3
               do i = 0,3
                  ptmp(i) = qph(k,i)
               enddo
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p3q*p3(0)*ptmp(0))
            enddo
         endif
         if (n4.gt.0) then
            do k = n1+n2+n3+1,nphot
               do i = 0,3
                  ptmp(i) = qph(k,i)
               enddo
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p4q*p4(0)*ptmp(0))
            enddo
         endif
         regulator2 = 1.d0* ncharged**nphot * anum/den
      endif
      w = w * regulator2
!-- extra regulator for gg ------------------------
      if (fs.eq.'gg') then
! I try a more symmetryc regulator ! test
         call regulate_gg(nphot,p1,p2,p3,p4,qph,reg)
         w=w/regulator2*reg
      endif
!-------------------------------------------------
#if defined KLOEI || defined B || defined BES3
! OK for nphot = 2      if (nphot.eq.2) then
!         k2 = 1
!         if (k2.eq.kgiag) k2 = 2
!         if (qph(k2,0).lt.0.02d0) w = w*2
!      endif
      k2 = 0
      if (kgiag.gt.0.and.kgiag.le.nis) then
         do k = 1,nphot
            om    = qph(k,0)
            cg    = qph(k,3)/om
            cgmin = cos(thgthigh)
            cgmax = cos(thgtlow)
            lcond = om.gt.omin!.and.abs(cg).lt.cgmax
            if (lcond) then
               k2 = k2 + 1      ! conta quelli dentr
            endif
         enddo
         w = w * nphot/k2      ! sarebbe nphot!/(nphot-1)!/k2 dove k2 conta quelli dentro
      else
         do k = 1,nphot
            om    = qph(k,0)
            lcond = om.gt.omin
            if (lcond) then
               k2 = k2 + 1      ! conta quelli dentro
            endif
         enddo
         w = w * nphot/k2       ! sarebbe nphot!/(nphot-1)!/k2 dove k2 conta quelli dentro
      endif
#endif
!-------------------------------------------------
      if (isnan(w)) then
         print*,''
         print*,'at the end of phasespace w = ',w
         print*,'wph,wen,wverg,p1jac,regulator1,regulator2,reg,nphot,
     .wten,qph'
         print*,wph,wen,wverg,p1jac,regulator1,regulator2,reg,nphot,
     .        wten(1:nphot)
         print*,qph(1:nphot,:)
         print*,''
      endif
      return
      end
**********************************************************************
********************PHASE SPACE CARLO RADIATIVE PROCESS***************
**********************************************************************
******************************************************************************
      subroutine phasespaceCRF(p1,p2,p3,p4,qph,nphot,m1,m2,
     .     esoft,cth,w,phsp,ie)
! derived from phasespaceStrong2020v2
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),Q(0:3),X(0:3)
      dimension dir3(0:3),dir4(0:3),ptmp(0:3),vers(0:3)
      double precision kone(0:3),pcoll(0:3),ptld(0:3),Kfs(0:3),xi(0:3)
      logical lcond,sol1,sol2,li
      double precision regf(4,nphot)
      character*2 fs
      common/finalstate/fs
      parameter (imaxph = 40)
      integer emissionpattern(nphot),epreg(nphot)
      integer empat(imaxph)
      integer*8 pow4(0:imaxph),j8
      integer pow2(0:imaxph)
      common/emissionpattern/empat
      parameter (pi = 3.1415926535897932384626433832795029d0)
      double precision qph(imaxph,0:3),oph(imaxph),wten(nphot)
      double precision csi(1),csi2(2),rnddummy,qbck(nphot,0:3)
      common/for_debug/pcoll,denphsp,dir3,dir4,br,bq,onm
      common/channelref/iref
      common/iwide/iwider,ntot
      common/regulators/ich,idummy,regulator1,regulator2,wcoll3,wcoll4
      common/radpattern/nph(4)
      common/idebugging/idebug
      common/idarkon/idarkon
      common/phspifirstpows/pow4,pow2,ifirst
      data ifirst /0/
      integer ice(10)
      common/counterv2/ic,ice
      data ic,ice/0,10*0/
      double precision lambda
      external lambda
      real strtm,endtm
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/massainv/amassainvmin,amassainvmax
      double precision prob_q2(20)
      common/q2samplingprobabilities/prob_q2
      common/todistribute/alford

      dimension versarr(nphot,0:3),VV(0:3)
      dimension xap(nphot),yap(nphot)
      common/testq2/tq2
      data tq2/-20d0/

#include "strong2020common.F"

      if (ifirst.eq.0) then
         ifirst = 1
         pow2(0) = 1
         pow4(0) = 1
         do k = 1,imaxph
            pow2(k) = pow2(k-1)*2
            pow4(k) = pow4(k-1)*4
         enddo
      endif

      w     = 1.d0

      iwider = 0
      ntot   = nphot

      qph = 0.d0
      nph = 0

      if (p1(0)+p2(0).lt.(m1+m2+nphot*esoft)) then
         phsp = 0.d0
         w    = 0.d0
         ie = 1
         return
      endif

      phsp  = 1.d0
      npart = nphot + 2
      nvar  = 3*nphot + 2

      ie    = 0
      duepigrechi = (2.d0*pi)**(-nvar)

      alford = -10d0

*****************
!! I put this outside in order to raise efficiency
!      call get_cos_fer(rnd(1),c,wcos)
!      cth = c
!      w   = w * wcos
c      c = cth
! mh, no, reputting it here!!
      call get_cos(fs,nphot,c,wcos)
      cth = c
      w   = w * wcos

*
      if (nphot.eq.0) then
         mass = p1(0)+p2(0)
         sqla = dsqrt(lambda(mass**2,m1**2,m2**2))
         phsp = phsp * sqla/8.d0/mass**2
         s = dsqrt(1.d0 - c**2)
         call getrnd(csi,1)
         call getphi(csi(1),phi,wph)
         w = w*wph
         sphi = sin(phi)
         cphi = cos(phi)
         p1mod = sqla/2.d0/mass
         p3(0) = dsqrt(p1mod**2 + m1**2)
         p3(1) = p1mod * sphi * s
         p3(2) = p1mod * cphi * s
         p3(3) = p1mod * c
         p4(0) = mass - p3(0)
         p4(1) = - p3(1)
         p4(2) = - p3(2)
         p4(3) = - p3(3)
         phsp = phsp * duepigrechi
         return
      endif
*******************
******************

      if (fs.eq.'gg') then
         print*,'this phase space should not be used for ''gg'''
         print*,'Stopping...'
         stop
      endif

* the following for ng >= 1
      ncharged = 4
c      ncharged = 2
      if (fs.eq.'pi') ncharged = 2

! cos(th) of electron already generated...
      sth  = sqrt(1.d0 - cth*cth)
      call getrnd(csi,1)
      call getphi(csi(1),phi,wph)
      w    = w*wph
      sphi = sin(phi)
      cphi = cos(phi)

      call get_pattern(ncharged,nphot,emissionpattern)
      do k = 1,nphot
         nph(emissionpattern(k)) = nph(emissionpattern(k)) + 1
      enddo

      n1 = nph(1)
      n2 = nph(2)
      n3 = nph(3)
      n4 = nph(4)
      nis = n1 + n2
      nfs = n3 + n4

      ei = p1(0)+p2(0)
      s  = dot(p1+p2,p1+p2)

      vers    =  0.d0
      vers(0) =  1.d0
      pcoll   =  p1+p2

c      call getrnd(csi,1)
c      iref = 3
c      if (csi(1).gt.0.5d0) iref = 4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      iref = 3  ! no need anymore for iref 3 or 4...??
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      w = w*wref
c      if (fs.eq.'pi') iref = 3

      if (iref.eq.3) then
         vx = sth*cphi
         vy = sth*sphi
         vz = cth
      else
         vx = -sth*cphi
         vy = -sth*sphi
         vz = -cth
      endif

      dir3(0) =  ei*0.5d0
      bb      =  sqrt(1.d0-m1*m1/dir3(0)/dir3(0))
      dir3(1) =  dir3(0)*bb*vx
      dir3(2) =  dir3(0)*bb*vy
      dir3(3) =  dir3(0)*bb*vz
      dir4    = -dir3
      dir4(0) =  dir3(0)

** PHOTONS GENERATION...
      pflatc  = 0.d0
      pflatc3 = 0.d0
      pflatc4 = 0.d0
      pflate  = 0.d0
      expgg = 0.92d0
      ibich = 1
      emaxgamma = (p1(0)+p2(0))*0.5d0

      xxx = lambda(s,0.d0,(m1+m2)**2)/s
      xxx = 0.5d0*sqrt(xxx)
      emaxgamma = xxx

c      eps  = 2.d0*esoft/sqs
      eps  = esoft/Emaxgamma

      do k = 1,nphot
         wen = 1.d0
         en  = esoft
         wten(k) = wen
         oph(k) = en
      enddo

      kgiag = -10*nphot ! this recovers the default
      klast = -10

      if (nis.gt.0) then
         call getrnd(csi,1)
         ndiff = nis*csi(1)+1
         i = 0
         do k = 1,nphot
            iepk = emissionpattern(k)
            if (iepk.le.2) then
               i = i + 1
               empat(i) = k
            endif
         enddo
         klast = empat(ndiff)
      endif
! if here i put klast = -10 I recover the original phase space
c      klast = -10

*** thi is in a common
      empat = 0
      empat(1:nphot) = emissionpattern
*****

**** phsp = dQ^2*d\Phi(P- > Q + K_initialstate)*d\Phi(Q-> p3 + p4 + K_finalstate)

**** initialization
      Q    = p1 + p2
      kone = 0.d0
      Kfs  = 0.d0
      prodomis = 1.d0
      prodomfs = 1.d0
****
** initial state photons
      do k = 1,nphot
         iepk = emissionpattern(k)
         if (iepk.le.2) then
            IF (k.ne.klast) THEN
               call getrnd(csi2,2)
               pfl = 0.d0
               if (iepk.eq.1) then
                  call collinear(pfl,csi2(1),csi2(2),p1,vers,wverg)
               else
                  call collinear(pfl,csi2(1),csi2(2),p2,vers,wverg)
               endif
               w = w * wverg
               call photon_energyv2(pflate,
     .              emaxgamma,esoft,en,wen)
               w = w * wen
               oph(k) = en
               prodomis  = prodomis*en

               qph(k,:) = oph(k)*vers
               pcoll = pcoll - qph(k,:)
               kone = kone + qph(k,:)
            ELSE
               continue
            ENDIF
         endif
      enddo
      if (
     .     (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0)) then
         ie = 1
         phsp = 0.d0
         w  = 0.d0
         ice(1) = ice(1) + 1
         return
      endif
      if (klast.gt.0.and.klast.le.nphot) then
! the first part of the phase space dQ^2*d\Phi(P- > Q + K_initialstate)
         iepk = emissionpattern(klast)
         call getrnd(csi2,2)
         pfl = 0.d0
         if (iepk.eq.1) then
            call collinear(pfl,csi2(1),csi2(2),p1,vers,wverg)
         else
            call collinear(pfl,csi2(1),csi2(2),p2,vers,wverg)
         endif

         w = w * wverg

         sqs  = sqrt(s)
         si   = tridot(kone,vers)
         alam = sqs - kone(0) + si
         ro   = dot(p1 + p2 - kone, p1 + p2 - kone)

         q2min = max((m1+m2)**2,ro - sqs*alam)
         q2min = max((m1+m2)**2,ro - 2.d0*Emaxgamma*alam)

         pmodmin  = sqrt(emin*emin - m1*m1)
         bmin     = pmodmin/emin
c         q2mincut = 2.d0*m1*m1+2.d0*emin*emin*(1.d0-bmin*bmin*cos(zmax)) ! bug
         q2mincut = 2.d0*m1*m1+2.d0*emin*emin*(1.d0+bmin*bmin*cos(zmax))
         q2min    = max(q2min,q2mincut,amassainvmin)

         q2max = ro - eps*2.d0*Emaxgamma*alam
! i.e.         q2max = ro - esoft*2.d0*alam !

***   for nphot = 1 and the photon is from IS, this should ensure
***   the discriminant below (in the second part of phsp) to be > 0
c$$$         al     = vx*vers(1) + vy*vers(2) + vz*vers(3)
c$$$         alford = al
c$$$
c$$$         al2 = al*al
c$$$         x1  = (1.d0+al2) + sqrt(4.d0*al2 + s/m1/m1*(1.d0-al2))
c$$$         x1  = m1*m1*x1/(1.d0 - m1*m1/s*(1.d0-al2))
c$$$
c$$$         q2min_1 = x1
c$$$c         if (nphot.eq.1) q2min = max(q2min,q2min_1)
c$$$
c$$$!   if (-al.lt.0.d0) then ! "-" because al here refers to \vec(Q) anti-parallel to vec(vers)
c$$$         if (al.gt.0d0) then
c$$$!     --- this should ensure at least one positive solution for |p_3|
c$$$            q2min_1 = s*m1 / (sqrt(s) - m1)
c$$$c            if (nphot.eq.1) q2min = max(q2min,q2min_1)
c$$$         endif

*****************
         if (q2max.le.q2min) then
!            print*,q2max,q2min,ro,alam
            ie   = 1
            phsp = 0.d0
            w    = 0.d0
            ice(2) = ice(2) + 1
            return
         endif

**************Q2 sampling *********************
         prob_q2    = 1.d0
         sommatutto = sum(prob_q2)
         prob_q2    = prob_q2/sommatutto

         call sampleq2wR(q2max,q2min,ro,sqs,alam,m1,
     .        q2,wq2)

         w = w * wq2

*****************************************************************
         om = (ro-q2)/alam*0.5d0

         oph(klast)   = om
         qph(klast,:) = oph(klast)*vers

         kone  = kone  + qph(klast,:)
         pcoll = pcoll - qph(klast,:)
         Q     = p1 + p2 - kone

         prodomis = prodomis*om

         phsp = phsp *0.5d0/abs(Q(0) + om + si) * prodomis
      endif

      zm2    = dot(Q,Q)
      bboost = sqrt(1d0-zm2/Q(0)/Q(0))

      if ((zm2.lt.(m1+m2)**2.or.Q(0).lt.0.d0)) then
         ie = 1
         phsp = 0.d0
         w  = 0.d0
         ice(3) = ice(3) + 1
         return
      endif

**************** final state photons and Comoving Reference Frame
****  see appendix A of JHEP 06 (2015) 023 and Hto4l for the version here
c      s34k = dot(Q,Q)
c      if (nfs.gt.0) then
c         do k = 1,nphot
c            iepk = emissionpattern(k)
c            if (iepk.gt.2) then
c               call getrnd(csi2,2)
c               pfl = 0.d0
c               if (iepk.eq.3) then
c                  call collinear(pfl,csi2(1),csi2(2),p1,vers,wverg)
c               else
c                  call collinear(pfl,csi2(1),csi2(2),p2,vers,wverg)
c               endif
c               w = w * wvergc
c               emax = s34k-dot(Kfs,Kfs) - 4d0*m1*m1 -4d0*m1*Kfs(0)
c               emax = emax /(4.d0*m1+2.d0*dot(vers,Kfs))
c               call get_egammasoft(esoft,emax,eg,weg)
c               qph(k,:) = eg*vers
c               Kfs = Kfs + eg*vers
c               prodomfs = prodomfs*eg
c               w = w*weg
c            endif
c         enddo
c      endif
********************************************************
c      E3 = (sqrt(s34k+tridot(Kfs,Kfs)) - Kfs(0))*0.5d0
c      if (E3.lt.m1) then
c         ie = 1
c         w  = 0.d0
c         phsp = 0.d0
c         return
c      endif
c      b3 = sqrt(1.d0-m1*m1/E3/E3)
c      p3(0)   =  E3
c      p3(1:2) =  0d0
c      p3(3)   =  b3*E3
c      p4(0)   =  p3(0)
c      p4(1:3) = -p3(1:3)
c      s34 = dot(p3+p4,p3+p4)
c      phsp = phsp * b3 * prodomfs * s34/s34k *0.5d0
c      phsp = phsp/ (1.d0 + Kfs(0)/sqrt(s34))


!! I'M IN THE REFERENCE WHERE p1 and p2 ARE ALWAYS BACK TO BACK
!! check the notes for conventions and ideas
      qm2     = 4.d0*m1*m1
      VV      = 0.d0
      versarr = 0.d0
      yjac    = 1.d0
      if (nfs.gt.0) then
         ymin = eps !* sqrt((1.d0-bboost)/(1.d0+bboost)) ! without sqrt seems slightly better
         call get_yx(zm2,qm2,nfs,nphot,ymin,yap,xap,yjac,wxy,iexy)
         if (iexy.gt.0) then
            ie   = 1
            phsp = 0.d0
            w    = 0.d0
            return
         endif
         w  = w * wxy

******
         ood3mod =  1.d0/sqrt(tridot(dir3,dir3))
         dir3(0) =  sqrt(zm2)*0.5d0
         bb      =  sqrt(1.d0 - m1*m1/dir3(0)/dir3(0))
         dir3(1) =  dir3(0)*bb*dir3(1)*ood3mod
         dir3(2) =  dir3(0)*bb*dir3(2)*ood3mod
         dir3(3) =  dir3(0)*bb*dir3(3)*ood3mod
         dir4    = -dir3
         dir4(0) = -dir4(0)
******

         ksf = 0
         do k = 1,nphot
            iepk = emissionpattern(k)
            if (iepk.ge.3) then
               ksf = ksf + 1
               call getrnd(csi2,2)
               pfl = 0.d0
               if (iepk.eq.3) then
                  call collinear(pfl,csi2(1),csi2(2),dir3,vers,wverg)
               else
                  call collinear(pfl,csi2(1),csi2(2),dir4,vers,wverg)
               endif
               w  = w * wverg
               versarr(ksf,:) = xap(ksf)*vers
               VV             = VV + versarr(ksf,:) ! this is now V in the notes
            endif
         enddo
      endif

      oodensp = 1.d0 /(1.d0 + 0.25d0 * dot(VV,VV) + VV(0))
      sp      = zm2 * oodensp
      if (sp.lt.qm2) then
         ie   = 1
         phsp = 0.d0
         w    = 0.d0
         return
      endif
      sqspm   = sqrt(sp)*0.5d0
      versarr = versarr * sqspm

      Kfs      = 0.d0
      ksf      = 0
      if (nfs.gt.0) then
         do k = 1,nphot
            iepk = emissionpattern(k)
            if (iepk.ge.3) then
               ksf      =  ksf + 1
               qph(k,:) = versarr(ksf,:)
               Kfs      = Kfs + qph(k,:)
               prodomfs = prodomfs * qph(k,0)
            endif
         enddo
      endif

      E1    = sqspm
      b1    = sqrt( (1.d0-m1/E1)*(1.d0+m1/E1) )
      p1mod = E1 * b1
      p3(0) = E1
      p3(1) = p1mod*vx
      p3(2) = p1mod*vy
      p3(3) = p1mod*vz

      p4(0) =  p3(0)
      p4(1) = -p3(1)
      p4(2) = -p3(2)
      p4(3) = -p3(3)

      yjac  = yjac * (1.d0 + VV(0)) !

      phsp = phsp*b1*0.5d0*prodomfs*oodensp*oodensp
      phsp = phsp * sqspm**nfs * yjac

! this is Q in this frame - I boost back where Q is at rest and then in the lab!
      ptld = p3+p4+Kfs
      call new_boost(ptld,p3,p3,1)
      call new_boost(Q,p3,p3,-1)
      call new_boost(ptld,p4,p4,1)
      call new_boost(Q,p4,p4,-1)
      if (nfs.gt.0) then
         do k = 1,nphot
            iepk = emissionpattern(k)
            if (iepk.ge.3) then
               call new_boost(ptld,qph(k,:),qph(k,:),1)
               call new_boost(Q,qph(k,:),qph(k,:),-1)
               if (qph(k,0).lt.esoft) then
                  ie   = 1
                  phsp = 0.d0
                  w    = 0.d0
                  return
               endif
            endif
         enddo
      endif
*************************************************************
** recalculating kone
! not needed...
c     kone = 0.d0
c      do k = 1,nphot
c         kone = kone + qph(k,:)
c      enddo
**************************
c      print*,p3+p4+kone
c      stop

      if (iref.eq.4) call exchange_mom(p3,p4)

      phsp = phsp * duepigrechi*(0.5d0)**npart

***********************
!     "regulator"... I try only on FINAL STATE (active only for ee and mm...)
! this is not needed anymore with CRF techique!!
c$$$      if (n3.gt.0.and.n4.gt.0) then
c$$$         prodp3 = 1.d0
c$$$         prodp4 = 1.d0
c$$$         do k = 1,nphot
c$$$            iempk = emissionpattern(k)
c$$$            if (iempk.gt.2) then
c$$$               ptmp = qph(k,:)
c$$$               p3q  = dot(p3,ptmp)/p3(0)/ptmp(0)
c$$$               p4q  = dot(p4,ptmp)/p4(0)/ptmp(0)
c$$$               if (iempk.eq.3) then
c$$$                  prodp3 = prodp3/p3q
c$$$               else ! if (iempk.eq.4) then
c$$$                  prodp4 = prodp4/p4q
c$$$               endif
c$$$            endif
c$$$         enddo
c$$$         if (iref.eq.3) anum = prodp3
c$$$         if (iref.eq.4) anum = prodp4
c$$$         regulator1 = anum/(prodp3+prodp4) * 2.d0
c$$$         w = w * regulator1
c$$$      endif
*******************

*******breaking down the product in the denominator
      iswitch    = 1
      if (ncharged.eq.2) iswitch = 0
      regulator2 = 1.d0
      if (fs.eq.'pi') iswitch = 0
      do k = 1,nphot
         regf(1,k) =     1d0 / dot(p1,qph(k,:))*p1(0)
         regf(2,k) =     1d0 / dot(p2,qph(k,:))*p2(0)
         regf(3,k) = iswitch / dot(p3,qph(k,:))*p3(0)
         regf(4,k) = iswitch / dot(p4,qph(k,:))*p4(0)
      enddo
      anum2 = 1.d0
      Q = p3+p4
      do k = 1,nphot
         if (k.ne.klast) then
            if (emissionpattern(k).le.2) then
               anum2 = anum2*regf(emissionpattern(k),k)
            else
               Q = Q + qph(k,:)
               anum2 = anum2*regf(emissionpattern(k),k)
            endif
         endif
      enddo
      if (klast.gt.0.and.klast.le.nphot) then
         ro    = dot(Q+qph(klast,:),Q+qph(klast,:))
         anum2 = anum2*regf(emissionpattern(klast),klast)
         anum2 = anum2*fregen(qph(klast,0),dot(Q,Q),s,ro,m1)
      endif

*** gemini optimized
!     Initialize the total sum
c         call cpu_time(strtm)
      den3 = 0.d0
! Loop over all 2**nphot combinations using 'i' as a bitmask
      do i = 0, 2**nphot - 1
         cis = 1.d0             ! Default value for cis12 (used if na=0)
         cfs = 1.d0             ! Product term for "Group 1"
         na = 0                 ! Counter for "Group 0" (bit = 0)
         nb = 0                 ! Counter for "Group 1" (bit = 1)
         Q  = p3 + p4           ! Initialize vector sum
         X  = p1 + p2
! Inner loop to partition items based on the bits of 'i'
         do k = 1, nphot
! Efficiently get the (k-1)-th bit of i
            ii = ibits(i, k-1, 1)
            select case (ii)
            case (0)            ! Item k is in "Group 0"
               na = na + 1
c     cis = cis * (regf(1,k) + regf(2,k))
            case (1)            ! Item k is in "Group 1"
               nb = nb + 1
               Q = Q + qph(k,:)
               cfs = cfs * (regf(3,k) + regf(4,k))
            end select
! Store the bit value (0 or 1) for this item
            epreg(k) = ii
         enddo
! Calculate the dot product of the summed vector Q
         q2 = dot(q, q)
! Calculate cis12 term
         if (na.gt.0) then
! "Group 0" is not empty. Use the optimized O(N) calculation.
            cis12_prod  = 1.d0
            cis12_sum_b = 0.d0
! Loop once over all photons to find those in "Group 0"
            do k = 1, nphot
               if (epreg(k).eq.0) then
!     Add to the product P = product(A_k)
                  cis12_prod = cis12_prod * (regf(1,k) + regf(2,k))
!     Add to the sum S_B = sum(B_k)
                  ro = dot(Q+qph(k,:),Q+qph(k,:))
                  cis12_sum_b=cis12_sum_b+fregen(qph(k,0),q2,s,ro,m1)
               endif
            enddo
! cis12 = (P * S_B) / na
            cis12 = (cis12_prod * cis12_sum_b) / dble(na)
         else
! "Group 0" is empty, use the default value
            cis12 = 1.d0
         endif
! Add the contribution from this combination to the total sum
         den3 = den3 + cis12 * cfs
      enddo
      regulator2 = 1.d0* ncharged**nphot * anum2/den3
*********************************
      w = w * regulator2
!-------------------------------------------------
      return
      end
**************************************
      function fregen(om,q2,s,ro,m)
! must be dimensionless
      implicit double precision (a-h,m,o-z)
      double precision an(20)
      common/samplingq2cmnwR/an,amro,gro
      double precision prob_q2(20),pr(20)
      common/q2samplingprobabilities/prob_q2
      character*2 fs
      common/finalstate/fs

      double precision mres(10),gres(10),gresee(10),mres2(10),gres2(10)
      double precision rescale(10)
      common/sampledresonances/mres,gres,mres2,gres2,gresee,pr,rescale,
     .     nres
      double precision vbar(10),vbar2(10)
      common/vvbars/vbar,vbar2
! it must reflect more or less the sampling function in sampleq2wR
      fregen = sqrt(1.d0-4.d0*m*m/q2) / q2 * rescale(1) ! GeV^-2

      if (nres.gt.0) then
         do k = 1,nres
            if (pr(k+1).gt.0.d0) then
               ooden  = 1.d0/denBW(q2,mres2(k),gres2(k))
               fregen = fregen + vbar2(k)*mres2(k)*ooden*rescale(k+1)  ! GeV^-2
            endif
         enddo
      endif
      fregen = fregen / (ro-q2) *om*sqrt(ro)*ro ! ciumbia!
      return
      end
*************************************
      subroutine sampleq2wR(q2max,q2min,ro,sqs,alam,m,q2,w)
      implicit double precision (a-h,m,o-z)
      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi
      character*2 fs
      logical lcond
      common/finalstate/fs
      common/ialpharunning/iarun
      double precision prob_q2(20),pr(20)
      common/q2samplingprobabilities/prob_q2
      double precision csi2(2),csi(1)
      double precision an(20),anir(20),anbwq2(20),anbw(20)
      common/samplingq2cmnwRv2/an,anir,anbwq2,anbw,ifirst
      data an,ifirst/20*0.d0,0/

      integer ioffhvp ,inlovp
      common/offhadrVP_nlovp/ioffhvp,inlovp

      double complex ff

      common/warning_ffpi/iffpiwarning,iprintffpiwarn

      
      double precision mres(10),gres(10),xxx(10),gresee(10)
      double precision mres2(10),gres2(10)
      double precision rescale(10)
      common/sampledresonances/mres,gres,mres2,gres2,gresee,pr,rescale,
     .     nres
      double precision vbar(10),vbar2(10)
      common/vvbars/vbar,vbar2
      data vbar,vbar2/20*0.d0/
      data rescale/10*0.d0/
!      m_om  = 0.78248d0
!      g_om  = 0.00855d0
!      m_phi = 1.01947d0
!      g_phi = 0.00425d0
      data mres/
     .     0.76419892698171199d0, ! rho
     .     0.78278988240778891d0, ! omega
     .     1.0198676384707968d0,  ! phi
     .     3.0969197600196230d0,  ! J/Psi
     .     3.6861561280010915d0,  ! Psi(2S)
     .     9.4603013694067908d0,  ! Y
     .     10.023263041544904d0,  ! Y(2S)
     .     10.355200529202707d0,  ! Y(3S)
     .     -1.d0,
     .     -1.d0
     .     /
      data gres/
     .     0.149d0,  ! rho
     .     8.49d-3,  ! omega
     .     4.266d-3, ! phi
     .     92.9d-6,  ! J/Psi
     .     299.d-6,  ! Psi(2S)
     .     54.02d-6, ! Y
     .     31.98d-6, ! Y(2S)
     .     20.32d-6, ! Y(3S)
     .    -1.d0,
     .    -1.d0
     .     /
      data gresee/
     .     7.04d-6,  ! rho
     .     0.60D-6,  ! omega
     .     1.27D-6,  ! phi
     .     5.55D-6,  ! J/Psi
     .     2.36D-6,  ! Psi(2S)
     .     1.34D-6,  ! Y
     .     0.612D-6, ! Y(2S)
     .     0.443D-6, ! Y(3S)
     .    -1.d0,
     .    -1.d0
     .     /

****
*     if AK = 0
*     samples according sqrt(1-4*m*m/q2) / q2 / (ro - q2) + c*BW(q2)/(ro - q2)
*
*     if AK = 1
*     sqrt(1-4*m*m/q2) / q2 / (ro - q2)  --->
*     ---- > sqrt(1-4*m*m/q2) / q2 / (ro - q2)*[ 1+ (q2-q2min)^2/(ro -q2min)^2 ]
*     which adds AP vertex "feature"
****

      lcond = (fs.eq.'pp'.or.fs.eq.'pi'.or.fs.eq.'pr').and.iffpi.gt.0

      if (ifirst.eq.0) then

         if (lcond) then
!     from ffpi_strong2020 in ffpi.F
! the default mass and width for the rho seems to be better...
c     mres(1) = 0.77456d0
c     gres(1) = 0.14832d0
            mres(2) = 1.8668d0
            mres(3) = 2.2645d0
            gres(2) = 0.30334d0
            gres(3) = 0.11327d0
         endif

         do k = 1,10
            mres2(k) = mres(k)*mres(k)
            gres2(k) = gres(k)*gres(k)
         enddo

         do k = 1,10
            rescale(k) = 1.d0
            if (k.gt.1) then
               if (gres(k-1).gt.0.d0) then
                  continue
               else
                  rescale(k) = 0.d0
               endif
            endif
         enddo

         alp = 1.d0/137.03599908d0
         do k = 1,10
            vbar(k)  = 0.d0
            if (mres(k).gt.0.d0) vbar(k) = 3.d0*gresee(k)/mres(k)/alp
            vbar2(k) = vbar(k) *vbar(k)
         enddo

         if (lcond) then                        
            iprintffpiwarn = 0
            call ffpi(mres2(1),ff,ff2)
            vbar2(1) = gres(1)/mres(1)*ff2
            call ffpi(mres2(2),ff,ff2)
            vbar2(2) = gres(2)/mres(2)*(0.0614d0)**2*ff2
            call ffpi(mres2(3),ff,ff2)
            vbar2(3) = gres(3)/mres(3)*(0.0047d0)**2*ff2
            do k = 4,8
               call ffpi(mres2(k),ff,ff2)
               vbar2(k) = 9d0*(gres(k)/mres(k)/alp)**2*ff2
            enddo
            iprintffpiwarn = 1
         endif

         ifirst = 1
      endif

** genero il primo canale come sqrt(1-a/q2)/q2
*     * ora combino il polo coulombiano e quello soffice. combino quelli
*     * che prima erano i primi due canali e genero secondo sqrt(1-a/q2)/q2/(ro-q2)

      nres = 8*iarun*(1-ioffhvp)

      if (lcond) then
         rescale(1) = 0.d0
         nres = 3
         nres = nres + 5*iarun*(1-ioffhvp) ! re-enabling all resonances for pions
      endif

      nc   = 1 + nres

      if (nres.gt.0) then
         do k = 1,nres
!     IR part
            anir(k) = log((ro - q2min)/(ro-q2max))
!     BW x q2 part
            call BWxQ2sampl(-1,q2max,q2min,mres(k),gres(k),anbwq2k,q2,w)
            anbwq2(k) = anbwq2k
!     "pure" BW part
            anbw(k) = atan((q2max - mres2(k))/gres(k)/mres(k))
            xxx(k)  = atan((q2min - mres2(k))/gres(k)/mres(k))
            anbw(k) = anbw(k) - xxx(k)
            anbw(k) = anbw(k)/gres(k)/mres(k)

            an(1+k) = anir(k) + anbwq2(k) + (ro-2.d0*mres2(k)) * anbw(k)

         enddo
      endif
      a     = 4.d0*m*m

      sqrro = sqrt(1.d0 - a/ro)
      sqrma = sqrt(1.d0 - a/q2max)
      sqrmi = sqrt(1.d0 - a/q2min)
!      an(1) =
!     .     + 2.d0*(ro-a)/ro/ro/sqrro*atanh(sqrma/sqrro)-2.d0/ro*sqrma
!      an(1) = an(1)
!     .     - 2.d0*(ro-a)/ro/ro/sqrro*atanh(sqrmi/sqrro)-2.d0/ro*sqrmi
! or
c      an(1) = sqrro*
c     .     log((sqrro+sqrma)*(sqrro-sqrmi)/(sqrro-sqrma)/(sqrro+sqrmi))
c     .     - 2.d0*(sqrma-sqrmi)
c      an(1) = an(1)/ro

      AK = 0.d0
      Z  = 1.d0/(ro-q2min)

      AA = 2d0*sqrro/ro * (1.d0 + AK)

      BB = AK*Z* (Z*(2.d0*ro + a) - 4d0)

      CC = 1.d0 + AK * Z*Z*q2min**2
      CC = -2.d0/ro*CC

      DD = -AK*Z*Z

      aLkma = funLk(sqrro,sqrma)
      aLkmi = funLk(sqrro,sqrmi)
      aL1ma = funLk(1.d0,sqrma)
      aL1mi = funLk(1.d0,sqrmi)

      an(1) = AA * aLkma + BB*aL1ma + CC*sqrma + DD*q2max*sqrma
      an(1) = an(1) - AA * aLkmi - BB*aL1mi - CC*sqrmi - DD*q2min*sqrmi

c      test=fchg(q2max,ro,q2min,AK,Z,a)-fchg(q2min,ro,q2min,AK,Z,a)
c      print*,(an(1)-test)/an(1)

      an1w  = log(q2max/q2min) + log( (ro-q2min) / (ro-q2max) )
      ccc   = log((ro-q2min)/q2min)

      den = 0.d0
      do k = 1,nc
         if (k.eq.1) then
            pr(k) = an(k) * prob_q2(k)*rescale(k)
         else
            pr(k) = an(k) * prob_q2(k)*rescale(k)
     .           *vbar2(k-1)*mres2(k-1)/denBW(ro,mres2(k-1),gres2(k-1))
         endif
         den   = den + pr(k)
      enddo
      pr(1:nc) = pr(1:nc)/den

      call getrnd(csi,1)
      if (csi(1).lt.pr(1)) then
         irig = 1
         do while (irig.eq.1)
            call getrnd(csi,1)
            rrr = exp(an1w*csi(1) - ccc)
            q2 = ro * rrr / (1.d0 + rrr)

            sqr   = sqrt(1.d0-a/q2)
            ratio = sqr * (1.d0 + AK*Z*Z*(q2-q2min)**2)/(1.d0+AK)

            call getrnd(csi,1)
            if (csi(1).lt.ratio) irig = 0
         enddo
      else
         k = 0
         cumul = pr(1)
         do while (csi(1).gt.cumul)
            k = k + 1
            cumul = cumul + pr(1+k)
         enddo

         amax  = max(0.d0,ro-2.d0*mres2(k))
         pir   = anir(k)
         pbwq2 = anbwq2(k)
         pbw   = anbw(k)* amax
         antot = pir + pbwq2 + pbw
         pir   = pir / antot
         pbwq2 = pbwq2 / antot
         pbw   = pbw / antot

         istop = 0
         do while (istop.eq.0)
            call getrnd(csi,1)
            if (csi(1).lt.pir) then
               call getrnd(csi,1)
               q2 = ro - (ro-q2min)*exp(-anir(k) * csi(1))
            elseif (csi(1).lt.(pir+pbwq2)) then
               call BWxQ2sampl(-1,q2max,q2min,mres(k),gres(k),dull,q2,w)
               call BWxQ2sampl( 1,q2max,q2min,mres(k),gres(k),dull,q2,w)
            else
               call getrnd(csi,1)
               q2 = mres2(k)+gres(k)*mres(k)
     .              *tan(gres(k)*mres(k)*anbw(k)*csi(1)+xxx(k))
            endif
            if (amax.gt.0.d0) then
               istop = 1
            else
               ooden  = 1.d0/denBW(q2,mres2(k),gres2(k))
               fwrong = 1.d0/(ro - q2) + q2*ooden
               ftrue  = fwrong + (ro-2d0*mres2(k))*ooden
               call getrnd(csi,1)
               ratio = ftrue/fwrong
               if (csi(1).lt.ratio) istop = 1
            endif
         enddo
         sqr = sqrt(1.d0-a/q2)
      endif

      APnum = (1.d0 + AK*Z*Z*(q2-q2min)**2)
      w     = 1.d0/an(1) * pr(1) * sqr /q2 /(ro-q2) * APnum

      if (nres.gt.0) then
         do k = 1,nres
            if (pr(k+1).gt.0.d0) then
               ooden = 1.d0/denBW(q2,mres2(k),gres2(k))
               w = w + 1d0/an(1+k)*pr(1+k) * (
     .              1.d0 / (ro - q2) +
     .              ( q2 + ro - 2.d0*mres2(k) )*ooden
     .              )
            endif
         enddo
      endif
      w = 1.d0/w
      return
      end
**************************************************************************
      function funLk(k,u)
      implicit none
      double precision funLk,u,k
c      funLK = 0.5d0*log( (k+u) / (k-u) )
      funLK = atanh(u/k)
      return
      end
*********************************************************************
      function denBW(s,m2,g2)
      implicit none
      double precision denBW,s,m2,g2
      denBW = (s-m2)**2 + m2*g2
      return
      end
**************************************************************************
**************************************************************************
      subroutine sampleq2(q2max,q2min,ro,sqs,alam,m,q2,w)
      implicit double precision (a-h,m,o-z)
      double precision csi(1)
      common/samplingq2cmn/an1,an2
      double precision prob_q2(20)
      common/q2samplingprobabilities/prob_q2

      p1 = prob_q2(1)
      p2 = prob_q2(2)

** genero il primo canale come sqrt(1-a/q2)/q2
      a    = 4.d0*m*m
      an1  = -2.d0*sqrt(1.d0-a/q2max) + log(q2max) + 2d0*log(1.d0
     .     + sqrt(1.d0-a/q2max))
      an1  = an1 +2.d0*sqrt(1.d0-a/q2min) - log(q2min)
     .     - 2d0*log(1.d0 + sqrt(1.d0-a/q2min))
      an1w = log(q2max/q2min)

      an2 = log((ro-q2min)/ (ro-q2max) )

      pa = an1*p1/(an1 * p1 + an2*p2)
      pb = 1.d0 - pa

      call getrnd(csi,1)
      if (csi(1).lt.pa) then
         irig = 1
         do while (irig.eq.1)
            call getrnd(csi,1)
            q2 = q2min*exp(an1w*csi(1))
            call getrnd(csi,1)
            sqr = sqrt(1.d0-a/q2)
            if (csi(1).lt.sqr) irig = 0
         enddo
      else
         call getrnd(csi,1)
         q2  = ro + (q2min - ro)*exp(-an2*csi(1))
         sqr = sqrt(1.d0-a/q2)
      endif

      w   = 1.d0/(1.d0/an1*pa*sqr/q2 + 1.d0/an2*pb/(ro-q2))
      return
      end
***************************************************************************
      subroutine get_egammasoft(egmin,egmax,eg,w)
! From Hto4l
      implicit double precision (a-h,o-z)
      double precision rnd(1)

      call getrnd(rnd,1)
!**     as 1/E for egmin < E < egmax
!**   and (depending on the comments with !)
!**    1)   as 1/sqrt(E) + a in 0 < E < Egmin (a to make it continue)
!** or 2)   as 1/egmin in 0 < E < Egmin
      a     = 1.d0/egmin    ! 2)
!     .  - 1.d0/sqrt(egmin) ! 1)

      an0   = a*egmin       ! 2)
!     .  + 2.d0*sqrt(egmin) ! 1)

      anir  = log(egmax/egmin)
      antot = an0 + anir

      csiN = rnd(1) * antot
      if (csiN.le.an0) then
!         x = (-1.d0 + sqrt(1.d0+a*csiN))/a ! 1)
!         eg = x*x                          ! 1)
!         w = antot/(1.d0/x+a)              ! 1)
         eg = csiN/a                        ! 2)
         w  = antot/a                       ! 2)
      else
         eg  = egmin * exp(csiN-an0)
         w   = antot*eg
      endif
      return
      end
****************************************************
      subroutine BWxQ2sampl(imode,q2max,q2min,m,g,an,q2,w)
      implicit double precision (a-h,m,o-z)
      double precision csi(1)
!     from Horace, as described in AcerMC 3.5 paper, pag 32-34, hep-ph/0405247
!     It samples as q2 / ( (q2-m**2)**2 + m**2*g**2)

!     imode = -1, returns only the normalization
!     imode /= -1, actually sample

      common/BWxQ2samplingcmn/am2,amag,ooamag,aa,ooaa,emin,emax,F1min,
     .     F2min, d1,d2,ds

      if (imode.eq.-1) then
         am2    = m*m
         amag   = m*g
         ooamag = 1.d0/amag
         aa     = m/g           ! = am2/amag
         ooaa   = g/m           ! = 1.d0/aa

         emin   = (q2min - am2)*ooamag
         emax   = (q2max - am2)*ooamag

         F1min = aa    * atan(emin)
         F2min = 0.5d0 * log(emin*emin+1.d0)

         d1   = aa * atan(emax) - F1min
         d2   = 0.5d0 * log(emax*emax+1.d0) - F2min
         ds   = d1 + d2
         an   = ds
         return
      endif

      istop = 0
      do while (istop.eq.0)
         call getrnd(csi,1)
         if (csi(1).lt.d2/ds) then
            call getrnd(csi,1)  !! puttana eva!
            x =  d2*csi(1) + F2min
            e =  sqrt(exp(2.d0*x)-1.d0)
            istop = 1
         else
            call getrnd(csi,1)
            x = d1 * csi(1) + F1min
            e = tan(x*ooaa)
            if (e.lt.0.d0) then
               r = (aa+e)/aa
               call getrnd(csi,1)
               if (csi(1).ge.r) e = -e
               if ((e.lt.emax).and.(e.gt.emin)) istop = 1
            else
               istop = 1
            endif
         endif
      enddo
      q2 = amag * e + am2
      w  = ds * ((q2-am2)*(q2-am2) + amag*amag )/q2
      return
      end
**************************************************************************
      subroutine get_yx(s,qm2,n,ntot,ymin,y,x,yjac,w,ie)
! from Hto4l
      implicit double precision (a-h,o-z)
      double precision x(ntot),y(ntot)

!     last notes, 4/6/2015
      ie   = 0
      ymax = 1.d0 - qm2/s
      if (ymin.gt.ymax) then
         w    = 0.d0
         yjac = 0.d0
         x    = 0d0
         y    = 0d0
         ie = 1
         return
      endif
!      ymin = 1d-2
      w    = 1.d0
      yjac = 1.d0
      do k = 1,n
         call get_y(ymin,ymax,yn,wn)
         y(k) = yn
         w    = w * wn
      enddo

      sumy = sum(y(1:n))
      if (sumy.gt.ymax) then
         ie = 1
         w  = 0.d0
         return
      endif

      x    = y    / (1.d0 - sumy)
      yjac = 1.d0 / (1.d0 - sumy)**n

      return
      end
****************************************************
      subroutine get_y(ymin,ymax,y,w)
! From Hto4l
! inspired by AP vertex
      implicit double precision (a-h,o-z)
      double precision rnd(1)
      common/getycmn/a,an0,ifirst
      data ifirst /0/
      if (ifirst.eq.0) then
         a = 2.d0 - 2.d0*ymin +ymin*ymin
         a = a / ymin
         an0   = a*ymin
         ifirst = 0 ! always recalculate because ymin will change...
      endif

      anir = 2.d0 * log(ymax/ymin) - 2.d0 * (ymax-ymin)
     .     + 0.5d0 * (ymax*ymax - ymin*ymin)

      antot = an0 + anir

      ihit = 0
      do while(ihit.eq.0)
         call get_y_wrongAP(ymin,ymax,y,fwrong)
         fy = apfun(ymin,y)
         call getrnd(rnd,1)
         if (rnd(1).le.fy/fwrong) ihit = 1
      enddo

      w = antot/fy
      return
      end
****************************************************
      subroutine get_y_wrongAP(ymin,ymax,y,f)
! From Hto4l
      implicit double precision (a-h,o-z)
      double precision rnd(1)
      a     = 1.d0/ymin
      an0   = a*ymin

      anir  = log(ymax/ymin)
      antot = an0 + anir

      call getrnd(rnd,1)
      csiN = rnd(1) * antot
      if (csiN.le.an0) then
         y = csiN/a
         f = 2.d0*a
      else
         y = ymin * exp(csiN-an0)
         f = 2.d0/y
      endif
      return
      end
****************************************************
****************************************************
      function apfun(ymin,y)
! From Hto4l
      implicit double precision (a-h,o-z)
      common/apfuncmn/a,ifirst
      data ifirst /0/
      if (ifirst.eq.0) then
         a = 2.d0 - 2.d0*ymin +ymin*ymin
         a = a / ymin
         ifirst = 0  ! always recalculate because ymin will change...
      endif
      apfun = a
      if (y.lt.ymin) return
      apfun = 2.d0 - 2.d0*y + y*y
      apfun = apfun/y
      return
      end
****************************************************
***************************************************
      subroutine photon_energyv2(pflat,omax,omin,om,w)
      implicit double precision (a-h,o-z)
      real*8 rnd,csi(1),csib(1)
      character*2 fs
      common/finalstate/fs
      common/epssoft/eps

      double precision omin,anir,anorm,eps

      integer ifirst
      data ifirst /0/

      common/ph_en_v2/anir,anorm,epsi,ifirst

      psing = 1.d0 - pflat

      if (ifirst.eq.0) then
        ifirst = 1
!        omin = eps*omax
        anir  = log(omax/omin)
        epsi  = omin/omax
        anorm = -2.d0*log(epsi)-1.5d0+2.d0*epsi-0.5d0*epsi**2
      endif

      if (pflat.gt.0.d0) then
         call getrnd(csi,1)
      else
         csi(1) = 0.5
      endif

      if (fs.eq.'ee'.or.fs.eq.'mm'
     +     .or.fs.eq.'pp'.or.fs.eq.'pi'
     + 	   .or.fs.eq.'mr'.or.fs.eq.'pr') then
         if (csi(1).lt.psing) then
            call ap_vertexv2(y,x,epsi)
            om = omax*x
            wi = 1.d0/anorm*(1d0 + y**2)/(1.d0-y)
            w = omax*1.d0/(psing*wi + pflat/(omax-omin))
         else
            call getrnd(csi,1)
            om = (omax - omin)*csi(1) + omin
            x = om/omax
            y = 1.d0 - x
            wsing = 1.d0/anorm*(1+y**2)/(1.d0-y)
            w = 1.d0/(psing*wsing + pflat/(omax-omin))
         endif

      elseif (fs.eq.'gg') then
         print*,'SHOULD NOT GO HERE photon_energyv2! stopping...'
         stop
         a = expgg
         b = -a
         upb = 1.d0 + b
         an = 1.d0/upb*(omax - omin)**upb
         call getrnd(csi,1)
         call getrnd(csib,1)
         if (csi(1).lt.0.5d0) then
            om = omax - (omax-omin)*(1.d0-csib(1))**(1.d0/upb)
            wnf = anir*an*om*(omax-om)**a/(anir*om+an*(omax-om)**a)
            w = wnf  * 2.d0
         else
            om = omin*dexp(anir*csib(1))
            wnf = anir*an*om*(omax-om)**a/(anir*om+an*(omax-om)**a)
            w = wnf * 2.d0
         endif
      endif
      return
      end
***************************************************
      subroutine ap_vertexv2(x,omx,eps)
      implicit double precision (a-h,o-z)
*  x generation according to ap splitting function
*  (1+x^2)/(1-x), 0 <= x <= 1-eps
* added 15/10/2005, modified from babayaga, older releases by CMCC
      double precision r(2)

      integer ifirst
      data ifirst /0/
      save ifirst

      double precision alne
      save alne

      if (ifirst.eq.0) then
        ifirst = 1
        alne = dlog(eps)
      endif

      irigenera = 1
      do while(irigenera.eq.1)
         call getrnd(r,2)
         cx  = r(1)
         omx = dexp(cx*alne)
         x   = 1.d0 - omx
         rx  = r(2)*2.d0/omx
         px  = (1.d0 + x**2)/omx
         if (rx.lt.px) irigenera = 0
      enddo
      return
      end
***
**************************************************************************
      function fun_cos(c,fs)
      implicit none
      character*2 fs
      double precision fun_cos,c,c2,upc2
      double precision antot,an1,an2,fmax,fmin,cminh,cmaxh,an0
      integer ifirst
      common/paramsgetcos/antot,an1,an0,an2,fmax,fmin,cminh,cmaxh,ifirst

      if (c.le.cminh) then
         fun_cos = fmin
         return
      elseif (c.ge.cmaxh) then
         fun_cos = fmax
         return
      endif

      select case (fs)
      case ('ee')
!         c2 = c*c
!         upc2 = (1.d0+c)*(1.d0+c)
!         fun_cos = (1.d0+c2)/2.d0 + (4.d0+upc2)/(1.d0-c)/(1.d0-c)
!     .        - upc2/(1.d0 - c)
!     i.e.
         fun_cos = (c*c + 3d0)**2/(1.d0-c)**2*0.5d0
      case ('mm','mr')
         fun_cos = 1.d0 + c*c
      case ('pp','pi','pr')
         fun_cos = 1.d0 - c*c
      end select
      return
      end
**************************************************************************
      subroutine get_cos(fs,nph,c,w)
      implicit double precision (a-h,o-z)
      double precision rnd,csi(1),r(2)
      character*2 fs
      common/coseni/cmax,cmin
      common/enableifirst/ifrstenbld
      double precision anorm, fmax
      integer ifirst,nph
      common/paramsgetcos/antot,an1,an0,an2,fmax,fmin,cminh,cmaxh,ifirst
      common/rescaleminmax/rmin,rmax
      data rmin,rmax /1.d0,1.d0/
      data ifirst /0/

      ilarger = 0
      if (nph.gt.0) ilarger = 1

      cmaxh =  cmax
      cminh =  cmin
c      cmaxh = 0.3d0
c      cminh = -0.3d0

      select case (fs)
      case ('mm','mr')

         rmax = 0.1d0
         rmin = 0.10d0

         fmax  = 1.d0 + cmaxh**2
         fmax  = fmax * rmax
         fmin  = 1.d0 + cminh**2
         fmin  = fmin * rmin
         an1   = fmin * (1.d0 + cminh) * ilarger
         an2   = fmax * (1.d0 - cmaxh) * ilarger
         tetto = 1.d0 + max(cmaxh**2,cminh**2)
         an0   = (cmaxh-cminh) + 1.d0/3.d0*(cmaxh**3-cminh**3)

         antot = an0+an1+an2

         call getrnd(csi,1)
         csi(1) = csi(1) * antot
         if (csi(1).lt.an1) then
            call getrnd(csi,1)
            c = -1d0 + (1.d0+cminh)*csi(1)
         elseif (csi(1).gt.(an0+an1)) then
            call getrnd(csi,1)
            c = cmaxh + (1.d0 - cmaxh)*csi(1)
         else
            iagain = 1
            do while(iagain.eq.1)
               call getrnd(r,2)
               c  = (cmaxh-cminh)*r(1) + cminh
               f  = tetto*r(2)
               fc = (1.d0 + c*c)
               if (f.lt.fc) iagain = 0
            enddo
         endif
         w = antot / fun_cos(c,'mm')

      case ('pp','pi','pr')

         rmax  = 1.d0
         rmin  = 1.d0
         fmax  = 1.d0 - cmaxh**2
         fmax  = fmax * rmax
         fmin  = 1.d0 - cminh**2
         fmin  = fmin * rmin
         tetto = 1.d0
         an1   = fmin * (1.d0 + cminh) * ilarger
         an2   = fmax * (1.d0 - cmaxh) * ilarger
         an0   = cmaxh - cminh - 1.d0/3.d0*(cmaxh**3-cminh**3)

         antot = an0+an1+an2

         call getrnd(csi,1)
         csi(1) = csi(1)*antot
         if (csi(1).lt.an1) then
            call getrnd(csi,1)
            c = -1d0 + (1.d0+cminh)*csi(1)
         elseif (csi(1).gt.(an0+an1)) then
            call getrnd(csi,1)
            c = cmaxh + (1.d0 - cmaxh)*csi(1)
         else
            iagain = 1
            do while(iagain.eq.1)
               call getrnd(r,2)
               c  = (cmaxh-cminh)*r(1) + cminh
               f  = tetto*r(2)
               fc = (1.d0 - c*c)
               if (f.lt.fc) iagain = 0
            enddo
         endif
         w = antot / fun_cos(c,'pp')
      case ('ee')

         rmax = 0.1d0
         rmin = 0.5d0

         fmax  = (cmaxh*cmaxh + 3d0)**2/(1.d0-cmaxh)**2*0.5d0
         fmax  = fmax * rmax
         fmin  = (cminh*cminh + 3d0)**2/(1.d0-cminh)**2*0.5d0
         fmin  = fmin * rmin

         tetto = (max(cminh,cmaxh)**2 + 3d0)**2*0.5d0
         an1   = fmin * (1.d0 + cminh) * ilarger
         an2   = fmax * (1.d0 - cmaxh) * ilarger

         tmax = 1.d0 - cmaxh
         tmin = 1.d0 - cminh
         an0 = !31.d0/6.d0
     .        - 6*tmax + tmax*tmax - tmax**3/6.d0
     .        + 8.d0*(log(tmax) + 1.d0/tmax)
         an0 = an0
     .        ! - 31.d0/6.d0
     .        + 6*tmin - tmin*tmin + tmin**3/6.d0
     .        - 8.d0*(log(tmin) + 1.d0/tmin)

         antot = an0+an1+an2

c         print*,an1/antot,an0/antot,an2/antot

         call getrnd(csi,1)
         csi(1) = csi(1)*antot
         if (csi(1).lt.an1) then
            call getrnd(csi,1)
            c = -1d0 + (1.d0+cminh)*csi(1)
         elseif (csi(1).gt.(an0+an1)) then
            call getrnd(csi,1)
            c = cmaxh + (1.d0 - cmaxh)*csi(1)
         else
            iagain = 1
            do while(iagain.eq.1)
               call getrnd(r,2)
! as 1 / (1-c)**2
               c = ((cmaxh-cminh)*r(1)+(1.d0-cmaxh)*cminh)/
     >              (1.d0-cmaxh+(cmaxh-cminh)*r(1))

               ratio = fun_cos(c,'ee') * (1.d0-c)**2 / tetto

               if (r(2).lt.ratio) iagain = 0
            enddo
         endif
         w = antot / fun_cos(c,'ee')
      case default
         print*,'wrong fs in get_cos. Stopping!'
         stop
      end select

      return
      end
**************************************************************************

***************************************************************************
**********       NEW ROUTINE OPTIMISED FOR A SINGLE PHOTON ****************
***************************************************************************
      subroutine phasespaceStrong2020v2(p1,p2,p3,p4,qph,nphot,m1,m2,
     .     esoft,cth,w,phsp,ie)
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),Q(0:3)
      dimension dir3(0:3),dir4(0:3),ptmp(0:3),vers(0:3)
      double precision kone(0:3),pcoll(0:3),ptld(0:3),koner(0:3),xi(0:3)
      logical lcond,sol1,sol2
      character*2 fs
      common/finalstate/fs
      integer emissionpattern(nphot)
      integer empat(40)
      common/emissionpattern/empat
      parameter (pi = 3.1415926535897932384626433832795029d0)
      parameter (imaxph = 40)
      double precision qph(imaxph,0:3),oph(imaxph),wten(nphot)
      double precision csi(1),csi2(2),rnddummy
      common/for_debug/pcoll,denphsp,dir3,dir4,br,bq,onm
      common/various/beta,betafs
      common/channelref/iref
      common/iwide/iwider,ntot
      common/regulators/ich,idummy,regulator1,regulator2,wcoll3,wcoll4
      common/radpattern/nph(4)
      common/idebugging/idebug
      common/idarkon/idarkon
      double precision lambda
      external lambda
#include "strong2020common.F"
      w     = 1.d0
      phsp  = 1.d0
      npart = nphot + 2
      nvar  = 3*nphot + 2

      ie    = 0
      duepigrechi = (2.d0*pi)**(-nvar)

      iwider = 0
      ntot   = nphot

      do j = 0,3
         do k = 1,imaxph
            qph(k,j) = 0.d0
         enddo
         nph(j+1) = 0
      enddo
      if (p1(0)+p2(0).lt.(m1+m2+nphot*esoft)) then
         phsp = 0.d0
         w    = 0.d0
         ie = 1
         return
      endif

*****************
!! I put this outside in order to raise efficiency
!      call get_cos_fer(rnd(1),c,wcos)
!      cth = c
!      w   = w * wcos
      c = cth
*
      if (nphot.eq.0) then
         mass = p1(0)+p2(0)
         sqla = dsqrt(lambda(mass**2,m1**2,m2**2))
         phsp = phsp * sqla/8.d0/mass**2
         s = dsqrt(1.d0 - c**2)
         call getrnd(csi,1)
         call getphi(csi(1),phi,wph)
         w = w*wph
         sphi = sin(phi)
         cphi = cos(phi)
         p1mod = sqla/2.d0/mass
         p3(0) = dsqrt(p1mod**2 + m1**2)
         p3(1) = p1mod * sphi * s
         p3(2) = p1mod * cphi * s
         p3(3) = p1mod * c
         p4(0) = mass - p3(0)
         p4(1) = - p3(1)
         p4(2) = - p3(2)
         p4(3) = - p3(3)
         phsp = phsp * duepigrechi
         return
      endif
*******************
******************
      lcond = fs.eq.'ee'.or.fs.eq.'mm'
     +     .or.fs.eq.'pp'.or.fs.eq.'pi'
     +	   .or.fs.eq.'er'.or.fs.eq.'mr'
     +     .or.fs.eq.'pr'

      angle = 30.d0 *pi/180.d0
      angle = 0.d0 *pi/180.d0

      cpmax = cos(angle)
      cpmin = cos(pi - angle)
* the following for ng >= 1
      ncharged = 4
      if (fs.eq.'gg'.or.fs.eq.'pi') ncharged = 2
!      ncharged=2
      ir = 1

! cos(th) of electron already generated...
      sth  = sqrt(1.d0 - cth*cth)
      call getrnd(csi,1)
      call getphi(csi(1),phi,wph)
      w    = w*wph
      sphi = dsin(phi)
      cphi = dcos(phi)
      ir   = ir + 2
****  qui sopra ho calcolato il phi del X-

      empat = 0
      call get_pattern(ncharged,nphot,emissionpattern)
      do k = 1,nphot
         nph(emissionpattern(k)) = nph(emissionpattern(k)) + 1
      enddo
      empat(1:nphot) = emissionpattern

c$$$
c$$$      if(nph(3).ne.0.or.nph(4).ne.0) then
c$$$         ie=1
c$$$         phsp=0.d0
c$$$         w=0.d0
c$$$         return
c$$$      endif
c$$$
      n1 = nph(1)
      n2 = nph(2)
      n3 = nph(3)
      n4 = nph(4)
      nis = n1 + n2
      nfs = n3 + n4

      ei = p1(0)+p2(0)
      s  = dot(p1+p2,p1+p2)

      dir3(0) =  ei*0.5d0
      dir3(1) =  dir3(0)*betafs*cphi*sth
      dir3(2) =  dir3(0)*betafs*sphi*sth
      dir3(3) =  dir3(0)*betafs*cth
      vers    =  0.d0
      vers(0) =  1.d0
      ptld    =  p1+p2
      pcoll   =  p1+p2
      dir4    = -dir3
      dir4(0) =  dir3(0)

      call getrnd(csi,1)
      iref = 3
      if (csi(1).gt.0.5d0) iref = 4
      wref = 1.d0
      if (nfs.gt.0) then
         if (n4.eq.0) iref = 3
         if (n3.eq.0) iref = 4
      endif
      w = w*wref
      if (fs.eq.'gg'.or.fs.eq.'pi') iref = 3

****  up to here emissionpattern is used only to count multiplicity on the legs.
****  then n1 are along p1, n2 along p2, etc
****  I rewrite emissionpattern to reflect this, so I can use it below
c      print*, ' '
c      print*,iref,emissionpattern
      do k = 1,nphot
         if (iref.eq.3) then
            if (k.le.n1+n2+n3+n4) emissionpattern(k) = 4
            if (k.le.n1+n2+n3)    emissionpattern(k) = 3
         else
            if (k.le.n1+n2+n4+n3) emissionpattern(k) = 3
            if (k.le.n1+n2+n4)    emissionpattern(k) = 4
         endif
         if (k.le.n1+n2) emissionpattern(k) = 2
         if (k.le.n1)    emissionpattern(k) = 1
      enddo
      empat(1:nphot) = emissionpattern

***** (MARCO) se non capisco male qui ha semplicemente ordinato, cioè empat ha gli indici ordinati

****
** PHOTONS GENERATION...
      pflatc  = 0.d0
      pflatc3 = 0.d0
      pflatc4 = 0.d0
      pflate  = 0.d0
c$$$
c$$$      pflatc  = 1.d0
c$$$      pflatc3 = 1.d0
c$$$      pflatc4 = 1.d0
c$$$      pflate  = 1.d0

      expgg = 0.92d0
      ibich = 1
      emaxgamma = (p1(0)+p2(0))*0.5d0

! this is wrong!! ??????
c      s   = dot(p1+p2,p1+p2)
c      xxx = lambda(s,0.d0,(m1+m2)**2)/s
c      xxx = 0.5d0*sqrt(xxx)
c      emaxgamma = xxx
c      print*,emaxgamma,emaxgammawrong

      do k = 1,nphot
         wen = 1.d0
         en  = esoft
         if (fs.eq.'gg') then
            call photon_energy(ibich,pflate,expgg,p1(0),
     .           en,wen)
         endif
         wten(k) = wen
         if (fs.eq.'ee'.or.fs.eq.'mm'
     .        .or.fs.eq.'pp'.or.fs.eq.'pi'
     -	      .or.fs.eq.'pr'.or.fs.eq.'mr'
     .	      .or.fs.eq.'er') w = w * wen
         oph(k) = en
         ir = ir + 1
      enddo

      kgiag = -10*nphot ! this recovers the default

      klast = -10

!#if defined KLOEI || defined B || defined BES3
!      call getrnd(csi2,2)
!      ir    = ir + 2
!      kgiag = 1 + nphot*csi2(1)
!      omax  = emaxgamma
!      omin  = esoft
!
!#ifdef KLOEI
!      omin = egammamin
!#endif
!#ifdef BES3
!      omin = egammamin1
!#endif
!#ifdef B
!      omin = egammamin
!#endif
!c#ifndef B
!      an  = log(omax/omin)
!      en  = omin*exp(an*csi2(2))
!      wen = an*en
!      oph(kgiag) = en
!      w          = w * wen
!      wten(kgiag) = wen  not used for ee or mm
!#else
!
      if (fs.ne.'gg') then
c         call getrnd(csi,1)
c         klast = 1+nphot*csi(1)
         klast = nphot
      endif
!#endif


** initial state photons
      if (nis.gt.0) then
         if (n1.gt.0) then
            do k = 1,n1
               IF (k.ne.klast) THEN
               iaccendi = 1
               call getrnd(csi2,2)
               pfl = 0.d0

!#if defined KLOEI || defined KLOEII || defined B || defined BES3
!               pfl = 0.05d0
!cc               if (k.eq.kgiag) pfl = 0.2d0
!#endif
c               if (k.eq.kgiag) pfl = 0.2d0
c               if (k.ne.kgiag) then
                  call collinear(pfl,csi2(1),csi2(2),p1,vers,wverg)
c               else
c                  call collinearstr2020(csi2(1),csi2(2),p1,vers,wverg)
c     endif

               iaccendi = 0
               ir = ir + 2
               w = w * wverg
****

               if (lcond.and.k.ne.kgiag) then
                  cp = vers(3)/vers(0)
                  if (cp.gt.cpmax.or.cp.lt.cpmin) then

                     call photon_energy_new(ibich,pflate,expgg,
     .                    emaxgamma,esoft,ptld,vers,en,wen)
                  else
                     call photon_energy(ibich,pflate,expgg,
     .                    emaxgamma,en,wen)
                  endif
                  w = w * wen
                  oph(k) = en
               endif
***
               qph(k,:) = oph(k)*vers
               pcoll = pcoll - qph(k,:)
            ELSE
               continue
            ENDIF
            enddo
         endif
         if (lcond.and.
     .        (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0)) then
            ie = 1
            phsp = 0.d0
            w  = 0.d0
            return
         endif
         if (n2.gt.0) then
            do k = 1,n2
               IF (k+n1.ne.klast) THEN
               iaccendi = 2
               call getrnd(csi2,2)
               pfl = 0.d0
!#if defined KLOEI || defined KLOEII || defined B || defined BES3
!               pfl = 0.05d0
!ccc               if ((k+n1).eq.kgiag) pfl = 0.2d0
!#endif
c               if ((k+n1).eq.kgiag) pfl = 0.2d0
c               if ((k+n1).ne.kgiag) then
                  call collinear(pfl,csi2(1),csi2(2),p2,vers,wverg)
c               else
c                  call collinearstr2020(csi2(1),csi2(2),p2,vers,wverg)
c               endif
               iaccendi  = 0
               ir = ir + 2
               w = w * wverg
****
               if (lcond.and.(k+n1).ne.kgiag) then
                  cp = vers(3)/vers(0)
                  if (cp.gt.cpmax.or.cp.lt.cpmin) then
                     call photon_energy_new(ibich,pflate,expgg,
     .                    emaxgamma,esoft,ptld,vers,en,wen)
                  else
                     call photon_energy(ibich,pflate,expgg,
     .                    emaxgamma,en,wen)
                  endif
                  w = w * wen
                  oph(k+n1) = en
               endif
***
               qph(k+n1,:) = oph(k+n1)*vers
               pcoll = pcoll - qph(k+n1,:)
            ELSE
               continue
            ENDIF
            enddo
         endif
         if (lcond.and.
     .        (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0)) then
            ie = 1
            phsp = 0.d0
            w  = 0.d0
            return
         endif
      endif

** final state photons
      wcoll3 = 1.d0
      wcoll4 = 1.d0
      if (nfs.gt.0) then
         if (iref.eq.3) then
            if (n3.gt.0) then
               do k = 1,n3
                  IF (k+n1+n2.ne.klast) THEN
                  iaccendi = 3
                  call getrnd(csi2,2)
                  call collinear(pflatc3,csi2(1),csi2(2),dir3,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond.and.(k+n1+n2).ne.kgiag) then
                     cp = vers(3)/vers(0)
                     if (cp.gt.cpmax.or.cp.lt.cpmin) then

                        call photon_energy_new(ibich,pflate,expgg,
     .                       emaxgamma,esoft,ptld,vers,en,wen)
                     else
                        call photon_energy(ibich,pflate,expgg,
     .                       emaxgamma,en,wen)
                     endif
                     w = w * wen
                     oph(k+n1+n2) = en
                  endif
***
                  wcoll3 = wcoll3 * wverg
                  qph(k+n1+n2,:) = oph(k+n1+n2)*vers
                  pcoll = pcoll - qph(k+n1+n2,:)

               ELSE
                  continue
               ENDIF
               enddo
            endif
            if (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0) then
               ie = 1
               phsp = 0.d0
               w  = 0.d0
               return
            endif

            pcoll     = pcoll - dir3

            pcmodmu = 1.d0/sqrt(tridot(pcoll,pcoll))
            pcoll(0) = dir4(0)
            bcoll    = sqrt(1.d0 - m1**2/pcoll(0)**2)
            pcoll(1) = pcoll(1)*pcmodmu*bcoll*pcoll(0)
            pcoll(2) = pcoll(2)*pcmodmu*bcoll*pcoll(0)
            pcoll(3) = pcoll(3)*pcmodmu*bcoll*pcoll(0)

            if (n4.gt.0) then
               do k = 1,n4
                  IF (k+n1+n2+n3.ne.klast) THEN
                  iaccendi = 4
                  call getrnd(csi2,2)
                  call collinear(pflatc4,csi2(1),csi2(2),pcoll,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond.and.(k+n1+n2+n3).ne.kgiag) then
                     cp = vers(3)/vers(0)
                     if (cp.gt.cpmax.or.cp.lt.cpmin) then
                        call photon_energy_new(ibich,pflate,expgg,
     .                       emaxgamma,esoft,ptld,vers,en,wen)
                     else
                        call photon_energy(ibich,pflate,expgg,
     .                       emaxgamma,en,wen)
                     endif
                     w = w * wen
                     oph(k+n1+n2+n3) = en
                  endif
***
                  wcoll4 = wcoll4 * wverg
                  qph(k+n1+n2+n3,:) = oph(k+n1+n2+n3)*vers
               ELSE
                  continue
               ENDIF
               enddo
            endif
         elseif (iref.eq.4) then
            if (n4.gt.0) then
               do k = 1,n4
                  IF (k+n1+n2+n3.ne.klast) THEN
                  iaccendi = 4
                  call getrnd(csi2,2)
                  call collinear(pflatc4,csi2(1),csi2(2),dir4,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond.and.(k+n1+n2+n3).ne.kgiag) then
                     cp = vers(3)/vers(0)
                     if (cp.gt.cpmax.or.cp.lt.cpmin) then
                        call photon_energy_new(ibich,pflate,expgg,
     .                       emaxgamma,esoft,ptld,vers,en,wen)
                     else
                        call photon_energy(ibich,pflate,expgg,
     .                       emaxgamma,en,wen)
                     endif
                     w = w * wen
                     oph(k+n1+n2+n3) = en
                  endif
***
                  wcoll4 = wcoll4 * wverg
                  qph(k+n1+n2+n3,:) = oph(k+n1+n2+n3)*vers

                  pcoll = pcoll - qph(k+n1+n2+n3,:)
               ELSE
                  continue
               ENDIF
               enddo
            endif
            if (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0) then
               ie = 1
               phsp = 0.d0
               w  = 0.d0
               return
            endif

            pcoll    = pcoll - dir4

            pcmodmu  = 1.d0/sqrt(tridot(pcoll,pcoll))
            pcoll(0) = dir3(0)
            bcoll    = sqrt(1.d0 - m2**2/pcoll(0)**2)
            pcoll(1) = pcoll(1)*pcmodmu*bcoll*pcoll(0)
            pcoll(2) = pcoll(2)*pcmodmu*bcoll*pcoll(0)
            pcoll(3) = pcoll(3)*pcmodmu*bcoll*pcoll(0)

            if (n3.gt.0) then
               do k = 1,n3
                  IF (k+n1+n2.ne.klast) THEN
                  iaccendi = 3
                  call getrnd(csi2,2)
                  call collinear(pflatc3,csi2(1),csi2(2),pcoll,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond.and.(k+n1+n2).ne.kgiag) then
                     cp = vers(3)/vers(0)
                     if (cp.gt.cpmax.or.cp.lt.cpmin) then
                        call photon_energy_new(ibich,pflate,expgg,
     .                       emaxgamma,esoft,ptld,vers,en,wen)
                     else
                        call photon_energy(ibich,pflate,expgg,
     .                       emaxgamma,en,wen)
                     endif
                     w = w * wen
                     oph(k+n1+n2) = en
                  endif
***
                  wcoll3 = wcoll3 * wverg
                  qph(k+n1+n2,:) = oph(k+n1+n2)*vers
               ELSE
                  continue
               ENDIF
               enddo
            endif
         endif
      endif

**********************************************
      if (klast.gt.0.and.klast.le.nphot) then
         iepklast = emissionpattern(klast)
         call getrnd(csi2,2)
         if (iepklast.eq.1) then
            pf1 = 0.d0
            call collinear(pf1,csi2(1),csi2(2),p1,vers,wverg)
         elseif (iepklast.eq.2) then
            pf2 = 0.d0
            call collinear(pf2,csi2(1),csi2(2),p2,vers,wverg)
         elseif (iepklast.eq.3) then
            if (iref.eq.3) then
               pf3 = 0.d0
               call collinear(pf3,csi2(1),csi2(2),dir3,vers,
     .              wverg)
            else
               pf4 = 0.d0
               call collinear(pf4,csi2(1),csi2(2),pcoll,vers,
     .              wverg)
            endif
         elseif (iepklast.eq.4) then
            if (iref.eq.3) then
               pf3 = 0.d0
               call collinear(pf3,csi2(1),csi2(2),pcoll,vers,
     .              wverg)
            else
               pf4 = 0.d0
               call collinear(pf4,csi2(1),csi2(2),dir4,vers,
     .              wverg)
            endif
         else
            print*,'sth wrong'
         endif

         kone = 0.d0
         do k = 1,nphot
            if (k.ne.klast) kone = kone + qph(k,:)
         enddo

         sqs  = sqrt(s)
         eps  = 2.d0*esoft/sqs
         xi   = (p1 + p2 - kone)/sqs
         xi2  = dot(xi,xi)
         xidv = dot(xi,vers)

c         om = 0.5d0*(s*xi2-q2)/sqs/xidv
c         s*xi2 - s*xidv*eps > q2 >  s*xi2 - s*xidv

         q2max = s*(xi2-xidv*eps)
         q2min = max((m1+m2)**2,s*(xi2-xidv))
         if (q2max.le.q2min) then
            ie = 1
            phsp = 0.d0
            w  = 0.d0
            return
         endif

c     I must generate om as 1/om ie q2 as (s*xi2-q2)
c         an  = log((s*xi2-q2min)/(s*xi2 - q2max))
c         an1 = an
c     or
         an2  = log((s*xi2-q2min)/s/xidv/eps)
cc         call getrnd(csi,1)
cc         q2  = s*xi2 + (q2min - s*xi2)*exp(-an2*csi(1))
c         wq2 = an2*(s*xi2-q2)

         an1 = log(q2max/q2min)
         pa  = an1/(an1+an2)
         pb  = 1.d0 - pa

         call getrnd(csi2,2)
         if (csi2(1).lt.pa) then

            q2 = q2min*exp(an1*csi2(2))
         else
            q2  = s*xi2 + (q2min - s*xi2)*exp(-an2*csi2(2))
         endif
         wq2 = 1.d0/(1.d0/an1*pa/q2 + 1.d0/an2*pb/(s*xi2-q2))
         w   = w * wq2

         om  = 0.5d0*(s*xi2-q2)/sqs/xidv
         en  = om

c         call photon_energy(ibich,pflate,expgg,emaxgamma,en,wen)
c         w = w * wen

         w = w * wverg
         oph(klast)   = en
         qph(klast,:) = oph(klast)*vers
      endif
*********************************************

** recalculating kone
      kone = 0.d0
      do k = 1,nphot
         kone = kone + qph(k,:)
      enddo
      ptmp = p1 + p2 - kone

***** TEST
      if (fs.eq.'gg') then
         koner = kone-qph(nphot,:) ! kone ridotto dell'ultimo fotone

         pk2 = dot(koner,koner)
         br  = koner(0) - pk2/2.d0/ei
         bq  = koner(0) - tridot(koner,vers)

         onm = ei*0.5 - br
         onm = onm/(1.d0 - bq/ei)

         a  = ei*0.5d0 - koner(0) + pk2/ei*0.5d0
         b  = 1.d0  - dot(koner,vers)/ei
         c  = ei - koner(0) + tridot(koner,dir3)/betafs/dir3(0)
         d  = 1.d0 - tridot(vers,dir3)/betafs/dir3(0)

         ep = (ei*a - esoft*c)/(ei*b-esoft*d)

c         ep = ei*0.5d0

         if (.not.(ep.lt.esoft.or.ep.ge.a/b.or.ep.eq.p1(0))) then
c         if (.FALSE.) then
! qui rimpiazza l'ultimo fotone generato, se non si verifica la condizione
! la condizione non dovrebbe presentarsi quasi mai, c'era beta invece di betafs sopra!!
            call photonenergygg(ep,a,b,esoft,en,wen)

            w = w * wen

            if (isnan(en).or.isnan(wen)) then
               print*,'!!!!!!!',ep,en,wen
               print*,'a b c d dir3',a,b,c,d,dir3
            endif
            wwen = 1.d0
            if (nphot.gt.1) then
               do k = 1,nphot-1
                  wwen = wwen * wten(k)
               enddo
            endif
            w = w * wwen

            qph(nphot,0) = en
            qph(nphot,1) = en*vers(1)
            qph(nphot,2) = en*vers(2)
            qph(nphot,3) = en*vers(3)
            do j = 0,3
               kone(j) = koner(j)+qph(nphot,j)
               ptmp(j) = p1(j) + p2(j) - kone(j)
            enddo
         else
!     altrimenti si tiene quello gia' generato e si va avanti
c            print*,'OOOOPS in phasespace ep',ep,a/b,p1(0),
c     .           tridot(vers,dir3)/betafs/dir3(0),d,
c     .           ep-a/b
            do k = 1,nphot
               w = w * wten(k)
            enddo
         endif
      endif
***** TEST

      if (dot(ptmp,ptmp).lt.(m1+m2)**2.or.ptmp(0).lt.0.d0) then
         ie = 1
         phsp = 0.d0
         w  = 0.d0
         return
      endif
      Q      = ptmp
***
      if (iref.eq.3) then
         vx  = sth*cphi
         vy  = sth*sphi
         vz  = cth
      else
         vx  = -sth*cphi
         vy  = -sth*sphi
         vz  = -cth
      endif

      z    = sqrt(tridot(Q,Q))
      ga   = (Q(1)*vx + Q(2)*vy + Q(3)*vz)/z

      q02  = Q(0)*Q(0)
      ap   = q2 + m1*m1 - m2*m2
      bp   = 2*ga*z
      abig = bp*bp - 4.d0*q02
      bbig = ap*bp
      cbig = ap*ap - 4d0*q02*m1*m1

      discr = bbig*bbig - abig*cbig

      if (discr.lt.0.d0) then
         ie = 1
         phsp = 0.d0
         w  = 0.d0
         return
      endif
      discr = sqrt(discr)

      p1mod_1 = -bbig + discr
      p1mod_1 =  p1mod_1 / abig
      p1mod_2 = -bbig - discr
      p1mod_2 =  p1mod_2 / abig
      sol1    = .true.
      sol2    = .true.
      if (p1mod_1.lt.0.d0) sol1 = .false.
      if (p1mod_2.lt.0.d0) sol2 = .false.
      if ((ap+bp*p1mod_1).lt.2.d0*Q(0)*m1) sol1 = .false.
      if ((ap+bp*p1mod_2).lt.2.d0*Q(0)*m1) sol2 = .false.
      if ((Q(0)-sqrt(m1*m1+p1mod_1*p1mod_1)).lt.m2) sol1 = .false.
      if ((Q(0)-sqrt(m1*m1+p1mod_2*p1mod_2)).lt.m2) sol2 = .false.

      p1jac = 1.d0
      if (sol1.and.sol2) then
         call getrnd(csi,1)
         if (csi(1).lt.0.5d0) p1mod = p1mod_1
         if (csi(1).ge.0.5d0) p1mod = p1mod_2
         p1jac = 2.d0
      elseif (sol1.and..not.sol2) then
         p1mod = p1mod_1
      elseif (.not.sol1.and.sol2) then
         p1mod = p1mod_2
      elseif (.not.sol1.and..not.sol2) then
         phsp = 0.d0
         w    = 0.d0
         ie   = 3
         return
      endif

      w     = w * p1jac

      p3(0) = sqrt(p1mod*p1mod + m1*m1)
      p3(1) = p1mod * vx
      p3(2) = p1mod * vy
      p3(3) = p1mod * vz
      p4    = Q - p3
***
      if (fs.eq.'gg'.and.(p3(0).lt.esoft.or.p4(0).lt.esoft)) then ! needed!
         phsp = 0.d0
         w    = 0.d0
         ie = 1
         return
      endif
***
!      phasespace!
      prodomega = 1.d0
      do k = 1,nphot
         prodomega = prodomega * qph(k,0)
      enddo

      phsp = prodomega

!      vers     = qph(klast,:)/qph(klast,0)
!      extrajac = Q(0) + qph(klast,0) - sqs * tridot(xi,vers)
! or
      extrajac = sqs*xidv
      phsp     = phsp/abs(extrajac)

** I put then the d\Phi_2
      den  = p1mod * p4(0) + (p1mod - bp*0.5d0)*p3(0)
! also
!      den = p4(0)*p1mod - p3(0)*tridot(p4,p3)/p1mod
!      print*,(den-aden)/den
! also
!      den = (p3(0)*dot(p3,p4) - p4(0)*m1*m1)/p1mod
!      print*,(den-aden)/den
! also (if m1 = m2)
c      den = (p3(0)*q2-2d0*Q(0)*m1*m1)*0.5d0/p1mod
c      print*,(den-aden)/den

      phsp = phsp  * p1mod*p1mod*0.5d0/abs(den) ! it's already in extrajac!! /Q(0)
*********
      phsp = phsp*duepigrechi*0.5d0**(npart)
      if (iref.eq.4) call exchange_mom(p3,p4)

***********************
! "regulator"... I try only on FINAL STATE (active only for ee and mm...)
      flatconst = 0.d0 !!!
      den       = 0.d0
      if (n3.gt.0.and.n4.gt.0) then
         prodp3 = 1.d0
         prodp4 = 1.d0
         do k = nis+1,nphot
            ptmp = qph(k,:)
            p3q  = dot(p3,ptmp)/p3(0)/ptmp(0)
            p4q  = dot(p4,ptmp)/p4(0)/ptmp(0)
            if (k.le.(n1+n2+n3)) then
               prodp3 = prodp3/p3q
            endif
            if (k.gt.(n1+n2+n3)) then
               prodp4 = prodp4/p4q
            endif
         enddo
         if (iref.eq.3) anum = prodp3
         if (iref.eq.4) anum = prodp4
         regulator1 = anum/(prodp3+prodp4) * 2.d0
         w = w * regulator1
      endif
*******************
* second "regulator"
      den        = 1.d0
      anum       = 1.d0
      iswitch    = 1
      regulator2 = 1.d0
c$$$
c$$$******
c$$$      iswitch=0
c$$$*****

      if (fs.eq.'gg'.or.fs.eq.'pi') iswitch = 0
      if (fs.ne.'gg') then
         if (n1.gt.0) then
            do k = 1,n1
               ptmp = qph(k,:)
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p1q*p1(0)*ptmp(0))
            enddo
         endif
         if (n2.gt.0) then
            do k = n1+1,n1+n2
               ptmp  = qph(k,:)
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p2q*p2(0)*ptmp(0))
            enddo
         endif
         if (n3.gt.0) then
            do k = n1+n2+1,n1+n2+n3
               ptmp = qph(k,:)
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p3q*p3(0)*ptmp(0))
            enddo
         endif
         if (n4.gt.0) then
            do k = n1+n2+n3+1,nphot
               ptmp = qph(k,:)
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p4q*p4(0)*ptmp(0))
            enddo
         endif
         regulator2 = 1.d0* ncharged**nphot * anum/den
      endif

!      regulator2=1.d0

      w = w * regulator2

******* tests
cc      call improvedbornwitheikonal(nphot,p1,p2,p3,p4,qph(1:nphot,:),aa)
*********

!-- extra regulator for gg ------------------------
      if (fs.eq.'gg') then
! I try a more symmetryc regulator ! test
         call regulate_gg(nphot,p1,p2,p3,p4,qph,reg)
         w=w/regulator2*reg
      endif
!-------------------------------------------------
!#if defined KLOEI || defined B || defined BES3
!! OK for nphot = 2      if (nphot.eq.2) then
!!         k2 = 1
!!         if (k2.eq.kgiag) k2 = 2
!!         if (qph(k2,0).lt.0.02d0) w = w*2
!!      endif
!      k2 = 0
!      if (kgiag.gt.0.and.kgiag.le.nis) then
!         do k = 1,nphot
!            om    = qph(k,0)
!            cg    = qph(k,3)/om
!            cgmin = cos(thgthigh)
!            cgmax = cos(thgtlow)
!            lcond = om.gt.omin!.and.abs(cg).lt.cgmax
!            if (lcond) then
!               k2 = k2 + 1      ! conta quelli dentro
!            endif
!         enddo
!         w = w * nphot/k2      ! sarebbe nphot!/(nphot-1)!/k2 dove k2 conta quelli dentro
!      else
!         do k = 1,nphot
!            om    = qph(k,0)
!            lcond = om.gt.omin
!            if (lcond) then
!               k2 = k2 + 1      ! conta quelli dentro
!            endif
!         enddo
!         w = w * nphot/k2       ! sarebbe nphot!/(nphot-1)!/k2 dove k2 conta quelli dentro
!      endif
!#endif
!-------------------------------------------------
      if (isnan(w)) then
         print*,''
         print*,'at the end of phasespace w = ',w
         print*,'wph,wen,wverg,p1jac,regulator1,regulator2,reg,nphot,
     .wten,qph'
         print*,wph,wen,wverg,p1jac,regulator1,regulator2,reg,nphot,
     .        wten(1:nphot)
         print*,qph(1:nphot,:)
         print*,''
      endif

c      if (nfs.gt.0) then
c         if (iref.eq.3) print*,iref,p4/pcoll
c         if (iref.eq.4) print*,iref,p3/pcoll
c      endif

      return
      end
******************************************************************************
***************************************************************************

***************************************************************************
      subroutine phasespaceStrong2020v3(p1,p2,p3,p4,qph,nphot,m1,m2,
     .     esoft,cth,w,phsp,ie)
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),Q(0:3)
      dimension dir3(0:3),dir4(0:3),ptmp(0:3),vers(0:3)
      double precision kone(0:3),pcoll(0:3),ptld(0:3),koner(0:3),xi(0:3)
      logical lcond,sol1,sol2
      character*2 fs
      common/finalstate/fs
      integer emissionpattern(nphot)
      integer empat(40)
      common/emissionpattern/empat
      parameter (pi = 3.1415926535897932384626433832795029d0)
      parameter (imaxph = 40)
      double precision qph(imaxph,0:3),oph(imaxph),wten(nphot)
      double precision csi(1),csi2(2),rnddummy
      common/for_debug/pcoll,denphsp,dir3,dir4,br,bq,onm
      common/various/beta,betafs
      common/channelref/iref
      common/iwide/iwider,ntot
      common/regulators/ich,idummy,regulator1,regulator2,wcoll3,wcoll4
      common/radpattern/nph(4)
      common/idebugging/idebug
      common/idarkon/idarkon
      double precision lambda
      external lambda
#include "strong2020common.F"
      w     = 1.d0
      phsp  = 1.d0
      npart = nphot + 2
      nvar  = 3*nphot + 2

      ie    = 0
      duepigrechi = (2.d0*pi)**(-nvar)

      iwider = 0
      ntot   = nphot

      do j = 0,3
         do k = 1,imaxph
            qph(k,j) = 0.d0
         enddo
         nph(j+1) = 0
      enddo

      if (p1(0)+p2(0).lt.(m1+m2+nphot*esoft)) then
         phsp = 0.d0
         w    = 0.d0
         ie = 1
         return
      endif

*****************
!! I put this outside in order to raise efficiency
!      call get_cos_fer(rnd(1),c,wcos)
!      cth = c
!      w   = w * wcos
      c = cth
*
      if (nphot.eq.0) then
         mass = p1(0)+p2(0)
         sqla = dsqrt(lambda(mass**2,m1**2,m2**2))
         phsp = phsp * sqla/8.d0/mass**2
         s = dsqrt(1.d0 - c**2)
         call getrnd(csi,1)
         call getphi(csi(1),phi,wph)
         w = w*wph
         sphi = sin(phi)
         cphi = cos(phi)
         p1mod = sqla/2.d0/mass
         p3(0) = dsqrt(p1mod**2 + m1**2)
         p3(1) = p1mod * sphi * s
         p3(2) = p1mod * cphi * s
         p3(3) = p1mod * c
         p4(0) = mass - p3(0)
         p4(1) = - p3(1)
         p4(2) = - p3(2)
         p4(3) = - p3(3)
         phsp = phsp * duepigrechi
         return
      endif
*******************
******************
      lcond = fs.eq.'ee'.or.fs.eq.'mm'
     +     .or.fs.eq.'pp'.or.fs.eq.'pi'
     +	   .or.fs.eq.'er'.or.fs.eq.'mr'
     +     .or.fs.eq.'pr'

      angle = 30.d0 *pi/180.d0
      angle = 0.d0 *pi/180.d0

      cpmax = cos(angle)
      cpmin = cos(pi - angle)
* the following for ng >= 1
      ncharged = 4
      if (fs.eq.'gg'.or.fs.eq.'pi') ncharged = 2
      ir = 1

! cos(th) of electron already generated...
      sth  = sqrt(1.d0 - cth*cth)
      call getrnd(csi,1)
      call getphi(csi(1),phi,wph)
      w    = w*wph
      sphi = dsin(phi)
      cphi = dcos(phi)
      ir   = ir + 2
****  qui sopra ho calcolato il phi del X-

      empat = 0
      call get_pattern(ncharged,nphot,emissionpattern)
      do k = 1,nphot
         nph(emissionpattern(k)) = nph(emissionpattern(k)) + 1
      enddo
      empat(1:nphot) = emissionpattern

      n1 = nph(1)
      n2 = nph(2)
      n3 = nph(3)
      n4 = nph(4)
      nis = n1 + n2
      nfs = n3 + n4

      ei = p1(0)+p2(0)
      s  = dot(p1+p2,p1+p2)

      dir3(0) =  ei*0.5d0
      dir3(1) =  dir3(0)*betafs*cphi*sth
      dir3(2) =  dir3(0)*betafs*sphi*sth
      dir3(3) =  dir3(0)*betafs*cth
      vers    =  0.d0
      vers(0) =  1.d0
      ptld    =  p1+p2
      pcoll   =  p1+p2
      dir4    = -dir3
      dir4(0) =  dir3(0)

      call getrnd(csi,1)
      iref = 3
      if (csi(1).gt.0.5d0) iref = 4
      wref = 1.d0
      if (nfs.gt.0) then
         if (n4.eq.0) iref = 3
         if (n3.eq.0) iref = 4
      endif
      w = w*wref
      if (fs.eq.'gg'.or.fs.eq.'pi') iref = 3

****  up to here emissionpattern is used only to count multiplicity on the legs.
****  then n1 are along p1, n2 along p2, etc
****  I rewrite emissionpattern to reflect this, so I can use it below
c      print*, ' '
c      print*,iref,emissionpattern
      do k = 1,nphot
         if (iref.eq.3) then
            if (k.le.n1+n2+n3+n4) emissionpattern(k) = 4
            if (k.le.n1+n2+n3)    emissionpattern(k) = 3
         else
            if (k.le.n1+n2+n4+n3) emissionpattern(k) = 3
            if (k.le.n1+n2+n4)    emissionpattern(k) = 4
         endif
         if (k.le.n1+n2) emissionpattern(k) = 2
         if (k.le.n1)    emissionpattern(k) = 1
      enddo
      empat(1:nphot) = emissionpattern

***** (MARCO) se non capisco male qui ha semplicemente ordinato, cioè empat ha gli indici ordinati

****
** PHOTONS GENERATION...
      pflatc  = 0.d0
      pflatc3 = 0.d0
      pflatc4 = 0.d0
      pflate  = 0.d0
c$$$
c$$$      pflatc  = 1.d0
c$$$      pflatc3 = 1.d0
c$$$      pflatc4 = 1.d0
c$$$      pflate  = 1.d0

      expgg = 0.92d0
      ibich = 1
      emaxgamma = (p1(0)+p2(0))*0.5d0

! this is wrong!! ??????
c      s   = dot(p1+p2,p1+p2)
c      xxx = lambda(s,0.d0,(m1+m2)**2)/s
c      xxx = 0.5d0*sqrt(xxx)
c      emaxgamma = xxx
c      print*,emaxgamma,emaxgammawrong

      do k = 1,nphot
         wen = 1.d0
         en  = esoft
         if (fs.eq.'gg') then
            call photon_energy(ibich,pflate,expgg,p1(0),
     .           en,wen)
         endif
         wten(k) = wen
         if (fs.eq.'ee'.or.fs.eq.'mm'
     .        .or.fs.eq.'pp'.or.fs.eq.'pi'
     -	      .or.fs.eq.'pr'.or.fs.eq.'mr'
     .	      .or.fs.eq.'er') w = w * wen
         oph(k) = en
         ir = ir + 1
      enddo

      kgiag = -10*nphot ! this recovers the default

      klast = -10

!#if defined KLOEI || defined B || defined BES3
!      call getrnd(csi2,2)
!      ir    = ir + 2
!      kgiag = 1 + nphot*csi2(1)
!      omax  = emaxgamma
!      omin  = esoft
!
!#ifdef KLOEI
!      omin = egammamin
!#endif
!#ifdef BES3
!      omin = egammamin1
!#endif
!#ifdef B
!      omin = egammamin
!#endif
!c#ifndef B
!      an  = log(omax/omin)
!      en  = omin*exp(an*csi2(2))
!      wen = an*en
!      oph(kgiag) = en
!      w          = w * wen
!      wten(kgiag) = wen  not used for ee or mm
!#else
!
      if (fs.ne.'gg') then
c         call getrnd(csi,1)
c         klast = 1+nphot*csi(1)
         klast = nphot
      endif
!#endif


** initial state photons
      if (nis.gt.0) then
         nphotel=0
         if (n1.gt.0) then
            do k = 1,n1
               nphotel=nphotel+1

               IF (k.ne.klast) THEN
               iaccendi = 1
               call getrnd(csi2,2)
               pfl = 0.d0
               if(nphotel.eq.1) then
                  pfl=1.d0
               endif
                  call collinear(pfl,csi2(1),csi2(2),p1,vers,wverg)

               iaccendi = 0
               ir = ir + 2
               w = w * wverg
****

               if (lcond.and.k.ne.kgiag) then
                  cp = vers(3)/vers(0)
                  if (cp.gt.cpmax.or.cp.lt.cpmin) then

                     call photon_energy_new(ibich,pflate,expgg,
     .                    emaxgamma,esoft,ptld,vers,en,wen)
                  else
                     call photon_energy(ibich,pflate,expgg,
     .                    emaxgamma,en,wen)
                  endif
                  w = w * wen
                  oph(k) = en
               endif
***
               qph(k,:) = oph(k)*vers
               pcoll = pcoll - qph(k,:)
            ELSE
               continue
            ENDIF
            enddo
         endif
         if (lcond.and.
     .        (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0)) then
            ie = 1
            phsp = 0.d0
            w  = 0.d0
            return
         endif
         if (n2.gt.0) then
            do k = 1,n2
               IF (k+n1.ne.klast) THEN
               iaccendi = 2
               call getrnd(csi2,2)
               pfl = 0.d0
               if(nphotel.eq.1) then
                  pfl=1.d0
               endif

!     #if defined KLOEI || defined KLOEII || defined B || defined BES3
!               pfl = 0.05d0
!ccc               if ((k+n1).eq.kgiag) pfl = 0.2d0
!#endif
c               if ((k+n1).eq.kgiag) pfl = 0.2d0
c               if ((k+n1).ne.kgiag) then
                  call collinear(pfl,csi2(1),csi2(2),p2,vers,wverg)
c               else
c                  call collinearstr2020(csi2(1),csi2(2),p2,vers,wverg)
c               endif
               iaccendi  = 0
               ir = ir + 2
               w = w * wverg
****
               if (lcond.and.(k+n1).ne.kgiag) then
                  cp = vers(3)/vers(0)
                  if (cp.gt.cpmax.or.cp.lt.cpmin) then
                     call photon_energy_new(ibich,pflate,expgg,
     .                    emaxgamma,esoft,ptld,vers,en,wen)
                  else
                     call photon_energy(ibich,pflate,expgg,
     .                    emaxgamma,en,wen)
                  endif
                  w = w * wen
                  oph(k+n1) = en
               endif
***
               qph(k+n1,:) = oph(k+n1)*vers
               pcoll = pcoll - qph(k+n1,:)
            ELSE
               continue
            ENDIF
            enddo
         endif
         if (lcond.and.
     .        (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0)) then
            ie = 1
            phsp = 0.d0
            w  = 0.d0
            return
         endif
      endif

** final state photons
      wcoll3 = 1.d0
      wcoll4 = 1.d0
      if (nfs.gt.0) then
         if (iref.eq.3) then
            if (n3.gt.0) then
               do k = 1,n3
                  IF (k+n1+n2.ne.klast) THEN
                  iaccendi = 3
                  call getrnd(csi2,2)
                  call collinear(pflatc3,csi2(1),csi2(2),dir3,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond.and.(k+n1+n2).ne.kgiag) then
                     cp = vers(3)/vers(0)
                     if (cp.gt.cpmax.or.cp.lt.cpmin) then

                        call photon_energy_new(ibich,pflate,expgg,
     .                       emaxgamma,esoft,ptld,vers,en,wen)
                     else
                        call photon_energy(ibich,pflate,expgg,
     .                       emaxgamma,en,wen)
                     endif
                     w = w * wen
                     oph(k+n1+n2) = en
                  endif
***
                  wcoll3 = wcoll3 * wverg
                  qph(k+n1+n2,:) = oph(k+n1+n2)*vers
                  pcoll = pcoll - qph(k+n1+n2,:)

               ELSE
                  continue
               ENDIF
               enddo
            endif
            if (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0) then
               ie = 1
               phsp = 0.d0
               w  = 0.d0
               return
            endif

            pcoll     = pcoll - dir3

            pcmodmu = 1.d0/sqrt(tridot(pcoll,pcoll))
            pcoll(0) = dir4(0)
            bcoll    = sqrt(1.d0 - m1**2/pcoll(0)**2)
            pcoll(1) = pcoll(1)*pcmodmu*bcoll*pcoll(0)
            pcoll(2) = pcoll(2)*pcmodmu*bcoll*pcoll(0)
            pcoll(3) = pcoll(3)*pcmodmu*bcoll*pcoll(0)

            if (n4.gt.0) then
               do k = 1,n4
                  IF (k+n1+n2+n3.ne.klast) THEN
                  iaccendi = 4
                  call getrnd(csi2,2)
                  call collinear(pflatc4,csi2(1),csi2(2),pcoll,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond.and.(k+n1+n2+n3).ne.kgiag) then
                     cp = vers(3)/vers(0)
                     if (cp.gt.cpmax.or.cp.lt.cpmin) then
                        call photon_energy_new(ibich,pflate,expgg,
     .                       emaxgamma,esoft,ptld,vers,en,wen)
                     else
                        call photon_energy(ibich,pflate,expgg,
     .                       emaxgamma,en,wen)
                     endif
                     w = w * wen
                     oph(k+n1+n2+n3) = en
                  endif
***
                  wcoll4 = wcoll4 * wverg
                  qph(k+n1+n2+n3,:) = oph(k+n1+n2+n3)*vers
               ELSE
                  continue
               ENDIF
               enddo
            endif
         elseif (iref.eq.4) then
            if (n4.gt.0) then
               do k = 1,n4
                  IF (k+n1+n2+n3.ne.klast) THEN
                  iaccendi = 4
                  call getrnd(csi2,2)
                  call collinear(pflatc4,csi2(1),csi2(2),dir4,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond.and.(k+n1+n2+n3).ne.kgiag) then
                     cp = vers(3)/vers(0)
                     if (cp.gt.cpmax.or.cp.lt.cpmin) then
                        call photon_energy_new(ibich,pflate,expgg,
     .                       emaxgamma,esoft,ptld,vers,en,wen)
                     else
                        call photon_energy(ibich,pflate,expgg,
     .                       emaxgamma,en,wen)
                     endif
                     w = w * wen
                     oph(k+n1+n2+n3) = en
                  endif
***
                  wcoll4 = wcoll4 * wverg
                  qph(k+n1+n2+n3,:) = oph(k+n1+n2+n3)*vers

                  pcoll = pcoll - qph(k+n1+n2+n3,:)
               ELSE
                  continue
               ENDIF
               enddo
            endif
            if (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0) then
               ie = 1
               phsp = 0.d0
               w  = 0.d0
               return
            endif

            pcoll    = pcoll - dir4

            pcmodmu  = 1.d0/sqrt(tridot(pcoll,pcoll))
            pcoll(0) = dir3(0)
            bcoll    = sqrt(1.d0 - m2**2/pcoll(0)**2)
            pcoll(1) = pcoll(1)*pcmodmu*bcoll*pcoll(0)
            pcoll(2) = pcoll(2)*pcmodmu*bcoll*pcoll(0)
            pcoll(3) = pcoll(3)*pcmodmu*bcoll*pcoll(0)

            if (n3.gt.0) then
               do k = 1,n3
                  IF (k+n1+n2.ne.klast) THEN
                  iaccendi = 3
                  call getrnd(csi2,2)
                  call collinear(pflatc3,csi2(1),csi2(2),pcoll,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond.and.(k+n1+n2).ne.kgiag) then
                     cp = vers(3)/vers(0)
                     if (cp.gt.cpmax.or.cp.lt.cpmin) then
                        call photon_energy_new(ibich,pflate,expgg,
     .                       emaxgamma,esoft,ptld,vers,en,wen)
                     else
                        call photon_energy(ibich,pflate,expgg,
     .                       emaxgamma,en,wen)
                     endif
                     w = w * wen
                     oph(k+n1+n2) = en
                  endif
***
                  wcoll3 = wcoll3 * wverg
                  qph(k+n1+n2,:) = oph(k+n1+n2)*vers
               ELSE
                  continue
               ENDIF
               enddo
            endif
         endif
      endif

**********************************************
      if (klast.gt.0.and.klast.le.nphot) then
         iepklast = emissionpattern(klast)
         call getrnd(csi2,2)
         if (iepklast.eq.1) then
            pf1 = 0.d0
            call collinear(pf1,csi2(1),csi2(2),p1,vers,wverg)
         elseif (iepklast.eq.2) then
            pf2 = 0.d0
            call collinear(pf2,csi2(1),csi2(2),p2,vers,wverg)
         elseif (iepklast.eq.3) then
            if (iref.eq.3) then
               pf3 = 0.d0
               call collinear(pf3,csi2(1),csi2(2),dir3,vers,
     .              wverg)
            else
               pf4 = 0.d0
               call collinear(pf4,csi2(1),csi2(2),pcoll,vers,
     .              wverg)
            endif
         elseif (iepklast.eq.4) then
            if (iref.eq.3) then
               pf3 = 0.d0
               call collinear(pf3,csi2(1),csi2(2),pcoll,vers,
     .              wverg)
            else
               pf4 = 0.d0
               call collinear(pf4,csi2(1),csi2(2),dir4,vers,
     .              wverg)
            endif
         else
            print*,'sth wrong'
         endif

         kone = 0.d0
         do k = 1,nphot
            if (k.ne.klast) kone = kone + qph(k,:)
         enddo

         sqs  = sqrt(s)
         eps  = 2.d0*esoft/sqs
         xi   = (p1 + p2 - kone)/sqs
         xi2  = dot(xi,xi)
         xidv = dot(xi,vers)

c         om = 0.5d0*(s*xi2-q2)/sqs/xidv
c         s*xi2 - s*xidv*eps > q2 >  s*xi2 - s*xidv

         q2max = s*(xi2-xidv*eps)
         q2min = max((m1+m2)**2,s*(xi2-xidv))
         if (q2max.le.q2min) then
            ie = 1
            phsp = 0.d0
            w  = 0.d0
            return
         endif

c     I must generate om as 1/om ie q2 as (s*xi2-q2)
c         an  = log((s*xi2-q2min)/(s*xi2 - q2max))
c         an1 = an
c     or
         an2  = log((s*xi2-q2min)/s/xidv/eps)
cc         call getrnd(csi,1)
cc         q2  = s*xi2 + (q2min - s*xi2)*exp(-an2*csi(1))
c         wq2 = an2*(s*xi2-q2)

         an1 = log(q2max/q2min)
         pa  = an1/(an1+an2)
         pb  = 1.d0 - pa

         call getrnd(csi2,2)
         if (csi2(1).lt.pa) then
            q2 = q2min*exp(an1*csi2(2))
         else
            q2  = s*xi2 + (q2min - s*xi2)*exp(-an2*csi2(2))
         endif
         wq2 = 1.d0/(1.d0/an1*pa/q2 + 1.d0/an2*pb/(s*xi2-q2))
         w   = w * wq2

         om  = 0.5d0*(s*xi2-q2)/sqs/xidv
         en  = om

c         call photon_energy(ibich,pflate,expgg,emaxgamma,en,wen)
c         w = w * wen

         w = w * wverg
         oph(klast)   = en
         qph(klast,:) = oph(klast)*vers
      endif
*********************************************

** recalculating kone
      kone = 0.d0
      do k = 1,nphot
         kone = kone + qph(k,:)
      enddo
      ptmp = p1 + p2 - kone

***** TEST
      if (fs.eq.'gg') then
         koner = kone-qph(nphot,:) ! kone ridotto dell'ultimo fotone

         pk2 = dot(koner,koner)
         br  = koner(0) - pk2/2.d0/ei
         bq  = koner(0) - tridot(koner,vers)

         onm = ei*0.5 - br
         onm = onm/(1.d0 - bq/ei)

         a  = ei*0.5d0 - koner(0) + pk2/ei*0.5d0
         b  = 1.d0  - dot(koner,vers)/ei
         c  = ei - koner(0) + tridot(koner,dir3)/betafs/dir3(0)
         d  = 1.d0 - tridot(vers,dir3)/betafs/dir3(0)

         ep = (ei*a - esoft*c)/(ei*b-esoft*d)

c         ep = ei*0.5d0

         if (.not.(ep.lt.esoft.or.ep.ge.a/b.or.ep.eq.p1(0))) then
c         if (.FALSE.) then
! qui rimpiazza l'ultimo fotone generato, se non si verifica la condizione
! la condizione non dovrebbe presentarsi quasi mai, c'era beta invece di betafs sopra!!
            call photonenergygg(ep,a,b,esoft,en,wen)

            w = w * wen

            if (isnan(en).or.isnan(wen)) then
               print*,'!!!!!!!',ep,en,wen
               print*,'a b c d dir3',a,b,c,d,dir3
            endif
            wwen = 1.d0
            if (nphot.gt.1) then
               do k = 1,nphot-1
                  wwen = wwen * wten(k)
               enddo
            endif
            w = w * wwen

            qph(nphot,0) = en
            qph(nphot,1) = en*vers(1)
            qph(nphot,2) = en*vers(2)
            qph(nphot,3) = en*vers(3)
            do j = 0,3
               kone(j) = koner(j)+qph(nphot,j)
               ptmp(j) = p1(j) + p2(j) - kone(j)
            enddo
         else
!     altrimenti si tiene quello gia' generato e si va avanti
c            print*,'OOOOPS in phasespace ep',ep,a/b,p1(0),
c     .           tridot(vers,dir3)/betafs/dir3(0),d,
c     .           ep-a/b
            do k = 1,nphot
               w = w * wten(k)
            enddo
         endif
      endif
***** TEST

      if (dot(ptmp,ptmp).lt.(m1+m2)**2.or.ptmp(0).lt.0.d0) then
         ie = 1
         phsp = 0.d0
         w  = 0.d0
         return
      endif
      Q      = ptmp
***
      if (iref.eq.3) then
         vx  = sth*cphi
         vy  = sth*sphi
         vz  = cth
      else
         vx  = -sth*cphi
         vy  = -sth*sphi
         vz  = -cth
      endif

      z    = sqrt(tridot(Q,Q))
      ga   = (Q(1)*vx + Q(2)*vy + Q(3)*vz)/z

      q02  = Q(0)*Q(0)
      ap   = q2 + m1*m1 - m2*m2
      bp   = 2*ga*z
      abig = bp*bp - 4.d0*q02
      bbig = ap*bp
      cbig = ap*ap - 4d0*q02*m1*m1

      discr = bbig*bbig - abig*cbig

      if (discr.lt.0.d0) then
         ie = 1
         phsp = 0.d0
         w  = 0.d0
         return
      endif
      discr = sqrt(discr)

      p1mod_1 = -bbig + discr
      p1mod_1 =  p1mod_1 / abig
      p1mod_2 = -bbig - discr
      p1mod_2 =  p1mod_2 / abig
      sol1    = .true.
      sol2    = .true.
      if (p1mod_1.lt.0.d0) sol1 = .false.
      if (p1mod_2.lt.0.d0) sol2 = .false.
      if ((ap+bp*p1mod_1).lt.2.d0*Q(0)*m1) sol1 = .false.
      if ((ap+bp*p1mod_2).lt.2.d0*Q(0)*m1) sol2 = .false.
      if ((Q(0)-sqrt(m1*m1+p1mod_1*p1mod_1)).lt.m2) sol1 = .false.
      if ((Q(0)-sqrt(m1*m1+p1mod_2*p1mod_2)).lt.m2) sol2 = .false.

      p1jac = 1.d0
      if (sol1.and.sol2) then
         call getrnd(csi,1)
         if (csi(1).lt.0.5d0) p1mod = p1mod_1
         if (csi(1).ge.0.5d0) p1mod = p1mod_2
         p1jac = 2.d0
      elseif (sol1.and..not.sol2) then
         p1mod = p1mod_1
      elseif (.not.sol1.and.sol2) then
         p1mod = p1mod_2
      elseif (.not.sol1.and..not.sol2) then
         phsp = 0.d0
         w    = 0.d0
         ie   = 3
         return
      endif

      w     = w * p1jac

      p3(0) = sqrt(p1mod*p1mod + m1*m1)
      p3(1) = p1mod * vx
      p3(2) = p1mod * vy
      p3(3) = p1mod * vz
      p4    = Q - p3
***
      if (fs.eq.'gg'.and.(p3(0).lt.esoft.or.p4(0).lt.esoft)) then ! needed!
         phsp = 0.d0
         w    = 0.d0
         ie = 1
         return
      endif
***
!      phasespace!
      prodomega = 1.d0
      do k = 1,nphot
         prodomega = prodomega * qph(k,0)
      enddo

      phsp = prodomega

!      vers     = qph(klast,:)/qph(klast,0)
!      extrajac = Q(0) + qph(klast,0) - sqs * tridot(xi,vers)
! or
      extrajac = sqs*xidv
      phsp     = phsp/abs(extrajac)

** I put then the d\Phi_2
      den  = p1mod * p4(0) + (p1mod - bp*0.5d0)*p3(0)
! also
!      den = p4(0)*p1mod - p3(0)*tridot(p4,p3)/p1mod
!      print*,(den-aden)/den
! also
!      den = (p3(0)*dot(p3,p4) - p4(0)*m1*m1)/p1mod
!      print*,(den-aden)/den
! also (if m1 = m2)
c      den = (p3(0)*q2-2d0*Q(0)*m1*m1)*0.5d0/p1mod
c      print*,(den-aden)/den

      phsp = phsp  * p1mod*p1mod*0.5d0/abs(den) ! it's already in extrajac!! /Q(0)
*********
      phsp = phsp*duepigrechi*0.5d0**(npart)
      if (iref.eq.4) call exchange_mom(p3,p4)

***********************
! "regulator"... I try only on FINAL STATE (active only for ee and mm...)
      flatconst = 0.d0 !!!
      den       = 0.d0
      if (n3.gt.0.and.n4.gt.0) then
         prodp3 = 1.d0
         prodp4 = 1.d0
         do k = nis+1,nphot
            ptmp = qph(k,:)
            p3q  = dot(p3,ptmp)/p3(0)/ptmp(0)
            p4q  = dot(p4,ptmp)/p4(0)/ptmp(0)
            if (k.le.(n1+n2+n3)) then
               prodp3 = prodp3/p3q
            endif
            if (k.gt.(n1+n2+n3)) then
               prodp4 = prodp4/p4q
            endif
         enddo
         if (iref.eq.3) anum = prodp3
         if (iref.eq.4) anum = prodp4
         regulator1 = anum/(prodp3+prodp4) * 2.d0
         w = w * regulator1
      endif
*******************
* second "regulator"
      den        = 1.d0
      anum       = 1.d0
      iswitch    = 1
      regulator2 = 1.d0
      if (fs.eq.'gg'.or.fs.eq.'pi') iswitch = 0
      if (fs.ne.'gg') then
         if (n1.gt.0) then
            do k = 1,n1
               ptmp = qph(k,:)
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p1q*p1(0)*ptmp(0))
            enddo
         endif
         if (n2.gt.0) then
            do k = n1+1,n1+n2
               ptmp  = qph(k,:)
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p2q*p2(0)*ptmp(0))
            enddo
         endif
         if (n3.gt.0) then
            do k = n1+n2+1,n1+n2+n3
               ptmp = qph(k,:)
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p3q*p3(0)*ptmp(0))
            enddo
         endif
         if (n4.gt.0) then
            do k = n1+n2+n3+1,nphot
               ptmp = qph(k,:)
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p4q*p4(0)*ptmp(0))
            enddo
         endif
         regulator2 = 1.d0* ncharged**nphot * anum/den
      endif
      w = w * regulator2

******* tests
cc      call improvedbornwitheikonal(nphot,p1,p2,p3,p4,qph(1:nphot,:),aa)
*********

!-- extra regulator for gg ------------------------
      if (fs.eq.'gg') then
! I try a more symmetryc regulator ! test
         call regulate_gg(nphot,p1,p2,p3,p4,qph,reg)
         w=w/regulator2*reg
      endif
!-------------------------------------------------
!#if defined KLOEI || defined B || defined BES3
!! OK for nphot = 2      if (nphot.eq.2) then
!!         k2 = 1
!!         if (k2.eq.kgiag) k2 = 2
!!         if (qph(k2,0).lt.0.02d0) w = w*2
!!      endif
!      k2 = 0
!      if (kgiag.gt.0.and.kgiag.le.nis) then
!         do k = 1,nphot
!            om    = qph(k,0)
!            cg    = qph(k,3)/om
!            cgmin = cos(thgthigh)
!            cgmax = cos(thgtlow)
!            lcond = om.gt.omin!.and.abs(cg).lt.cgmax
!            if (lcond) then
!               k2 = k2 + 1      ! conta quelli dentro
!            endif
!         enddo
!         w = w * nphot/k2      ! sarebbe nphot!/(nphot-1)!/k2 dove k2 conta quelli dentro
!      else
!         do k = 1,nphot
!            om    = qph(k,0)
!            lcond = om.gt.omin
!            if (lcond) then
!               k2 = k2 + 1      ! conta quelli dentro
!            endif
!         enddo
!         w = w * nphot/k2       ! sarebbe nphot!/(nphot-1)!/k2 dove k2 conta quelli dentro
!      endif
!#endif
!-------------------------------------------------
      if (isnan(w)) then
         print*,''
         print*,'at the end of phasespace w = ',w
         print*,'wph,wen,wverg,p1jac,regulator1,regulator2,reg,nphot,
     .wten,qph'
         print*,wph,wen,wverg,p1jac,regulator1,regulator2,reg,nphot,
     .        wten(1:nphot)
         print*,qph(1:nphot,:)
         print*,''
      endif

c      if (nfs.gt.0) then
c         if (iref.eq.3) print*,iref,p4/pcoll
c         if (iref.eq.4) print*,iref,p3/pcoll
c      endif

      return
      end
******************************************************************************
***************************************************************************


      subroutine phasespace(p1,p2,p3,p4,qph,nphot,m1,m2,esoft,
     .     cth,w,phsp,ie)
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension ap3(0:3),ap4(0:3),pcollorig(0:3)
      dimension p4guess(0:3),p3guess(0:3)
      dimension dir3(0:3),dir4(0:3),ptmp(0:3),vers(0:3)
      double precision kone(0:3),pcoll(0:3),ptld(0:3),koner(0:3)
      double precision lambda
      logical lcond
      character*2 fs
      common/finalstate/fs
      integer emissionpattern(nphot)
      parameter (pi = 3.1415926535897932384626433832795029d0)
      parameter (imaxph = 40)
      double precision qph(imaxph,0:3),oph(imaxph),wten(nphot)
      real*8 csi(1),csi2(2),rnddummy
      common/for_debug/pcoll,denphsp,dir3,dir4,br,bq,onm
      common/various/beta,betafs
      common/channelref/iref
      common/iwide/iwider,ntot
      common/regulators/ich,idummy,regulator1,regulator2,wcoll3,wcoll4
      common/radpattern/nph(4)
      common/idebugging/idebug
      common/idarkon/idarkon

      w      = 1.d0
      phsp   = 1.d0
      npart = nphot + 2
      nvar  = 3*nphot + 2

      ie     = 0

      duepigrechi = (2.d0*pi)**(-nvar)

      iwider = 0
      ntot   = nphot

      do j = 0,3
         do k = 1,imaxph
            qph(k,j) = 0.d0
         enddo
         nph(j+1) = 0
      enddo

      if (p1(0)+p2(0).lt.(m1+m2+nphot*esoft)) then
         phsp = 0.d0
         w    = 0.d0
         ie = 1
         return
      endif
*****************
!! I put this outside in order to raise efficiency
!      call get_cos_fer(rnd(1),c,wcos)
!      cth = c
!      w   = w * wcos
      c = cth
*
      if (nphot.eq.0) then
         mass = p1(0)+p2(0)
         sqla = dsqrt(lambda(mass**2,m1**2,m2**2))
         phsp = phsp * sqla/8.d0/mass**2
         s = dsqrt(1.d0 - c**2)
         call getrnd(csi,1)
         call getphi(csi(1),phi,wph)
         w = w*wph
         sphi = sin(phi)
         cphi = cos(phi)
         p1mod = sqla/2.d0/mass
         p3(0) = dsqrt(p1mod**2 + m1**2)
         p3(1) = p1mod * sphi * s
         p3(2) = p1mod * cphi * s
         p3(3) = p1mod * c
         p4(0) = mass - p3(0)
         p4(1) = - p3(1)
         p4(2) = - p3(2)
         p4(3) = - p3(3)
         phsp = phsp * duepigrechi
         return
      endif
*******************
******************
      lcond = fs.eq.'ee'.or.fs.eq.'mm'.or.fs.eq.'mr'


      angle = 30.d0 *pi/180.d0
      angle = 0.d0 *pi/180.d0

      cpmax = cos(angle)
      cpmin = cos(pi - angle)
* the following for ng >= 1
      ncharged = 4
      if (fs.eq.'gg') ncharged = 2
      ir = 1

! cos(th) of electron already generated...
      sth  = dsqrt(1.d0 - cth**2)
      call getrnd(csi,1)
      call getphi(csi(1),phi,wph)
      w    = w*wph
      sphi = dsin(phi)
      cphi = dcos(phi)
      ir   = ir + 2


      call get_pattern(ncharged,nphot,emissionpattern)
      do k = 1,nphot
         nph(emissionpattern(k)) = nph(emissionpattern(k)) + 1
      enddo
      n1 = nph(1)
      n2 = nph(2)
      n3 = nph(3)
      n4 = nph(4)
      nis = n1 + n2
      nfs = n3 + n4

      ei = p1(0)+p2(0)

      dir3(0) = ei/2d0
      dir3(1) = ei/2d0*betafs*cphi*sth
      dir3(2) = ei/2d0*betafs*sphi*sth
      dir3(3) = ei/2d0*betafs*cth
      dir4(0) = dir3(0)
      pcoll(0) = p1(0)+p2(0)
      ptld(0) = p1(0)+p2(0)
      vers(0)  = 0.d0
      do k = 1,3
         pcoll(k) = p1(k)+p2(k)
         ptld(k)  = p1(k)+p2(k)
         dir4(k)  = -dir3(k)
         vers(k)  = 0.d0
      enddo

      call getrnd(csi,1)
      iref = 3
      if (csi(1).gt.0.5) iref = 4
      wref = 1.d0
      if (nfs.gt.0) then
         if (n4.eq.0) iref = 3
         if (n3.eq.0) iref = 4
      endif
      w = w*wref
      if (fs.eq.'gg') iref = 3
****
** PHOTONS GENERATION...
      pflatc  = 0.d0
      pflatc3 = 0.d0
      pflatc4 = 0.d0
      pflate  = 0.d0
      expgg = 0.92d0
      ibich = 1
      do k = 1,nphot
         wen = 1.d0
         en  = esoft
         if (fs.eq.'gg') then
            call photon_energy(ibich,pflate,expgg,p1(0),
     .           en,wen)
         endif
         wten(k) = wen
         if (fs.eq.'ee'.or.fs.eq.'mm'.or.fs.eq.'mr'
     +       ) w = w * wen
         oph(k) = en
         ir = ir + 1
      enddo
** initial state photons
c  si vede su quale gamba ci sono dei fotoni e se ne calcola la
c  direzione con "collinear", poi se ne calcola l'energia, secondo l'AP.
c  con "photon_energy", "photon_energy_dm" calcola l'energia anche con la BW
c  Poi si toglie l'energia del fotone da quella dell'elettrone.

      emaxgamma = (p1(0)+p2(0))/2.d0
      if (nis.gt.0) then
         if (n1.gt.0) then
            do k = 1,n1
               iaccendi = 1
               call getrnd(csi2,2)
               call collinear(0.d0,csi2(1),csi2(2),p1,vers,wverg)
               iaccendi = 0
               ir = ir + 2
               w = w * wverg
****
               if (lcond) then
                  if (idarkon.eq.1) then
                      call photon_energy_dm(emaxgamma,en,wen)
                  else
                      cp = vers(3)/vers(0)
                      if (cp.gt.cpmax.or.cp.lt.cpmin) then
                         call photon_energy_new(ibich,pflate,expgg,
     .                              emaxgamma,esoft,ptld,vers,en,wen)
                      else
                        call photon_energy(ibich,pflate,expgg,
     .                                     emaxgamma,en,wen)
                      endif
                  endif
                  w = w * wen
                  oph(k) = en
               endif
***
               qph(k,0) = oph(k)
               qph(k,1) = oph(k)*vers(1)
               qph(k,2) = oph(k)*vers(2)
               qph(k,3) = oph(k)*vers(3)

               pcoll(0) = pcoll(0) - qph(k,0)
               pcoll(1) = pcoll(1) - qph(k,1)
               pcoll(2) = pcoll(2) - qph(k,2)
               pcoll(3) = pcoll(3) - qph(k,3)
            enddo
         endif
         if (n2.gt.0) then
            do k = 1,n2
               iaccendi = 2
               call getrnd(csi2,2)
               call collinear(0.d0,csi2(1),csi2(2),p2,vers,wverg)
               iaccendi  = 0
               ir = ir + 2
               w = w * wverg
****
               if (lcond) then
                  if (idarkon.eq.1) then
                      call photon_energy_dm(emaxgamma,en,wen)
                  else
                      cp = vers(3)/vers(0)
                      if (cp.gt.cpmax.or.cp.lt.cpmin) then
                         call photon_energy_new(ibich,pflate,expgg,
     .                              emaxgamma,esoft,ptld,vers,en,wen)
                      else
                        call photon_energy(ibich,pflate,expgg,
     .                                     emaxgamma,en,wen)
                      endif
                  endif
                  w = w * wen
                  oph(k+n1) = en
               endif
***
               qph(k+n1,0) = oph(k+n1)
               qph(k+n1,1) = oph(k+n1)*vers(1)
               qph(k+n1,2) = oph(k+n1)*vers(2)
               qph(k+n1,3) = oph(k+n1)*vers(3)

               pcoll(0) = pcoll(0) - qph(k+n1,0)
               pcoll(1) = pcoll(1) - qph(k+n1,1)
               pcoll(2) = pcoll(2) - qph(k+n1,2)
               pcoll(3) = pcoll(3) - qph(k+n1,3)
            enddo
         endif
      endif

      if (lcond.and.
     .   (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0)) then
         ie = 1
         phsp = 0.d0
         w  = 0.d0
         return
      endif

** final state photons
c      en3guess = dir3(0)
c      en4guess = dir4(0)
      wcoll3 = 1.d0
      wcoll4 = 1.d0
      if (nfs.gt.0) then
         if (iref.eq.3) then
            if (n3.gt.0) then
               do k = 1,n3
                  iaccendi = 3
                  call getrnd(csi2,2)
                  call collinear(pflatc3,csi2(1),csi2(2),dir3,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond) then
                    if (idarkon.eq.1) then
                       call photon_energy_dm(emaxgamma,en,wen)
                    else
                       cp = vers(3)/vers(0)
                       if (cp.gt.cpmax.or.cp.lt.cpmin) then
                            call photon_energy_new(ibich,pflate,expgg,
     .                              emaxgamma,esoft,ptld,vers,en,wen)
                       else
                            call photon_energy(ibich,pflate,expgg,
     .                                         emaxgamma,en,wen)
                       endif
                     endif
                     w = w * wen
                     oph(k+n1+n2) = en
                  endif
***
                  wcoll3 = wcoll3 * wverg
                  qph(k+n1+n2,0) = oph(k+n1+n2)
                  qph(k+n1+n2,1) = oph(k+n1+n2)*vers(1)
                  qph(k+n1+n2,2) = oph(k+n1+n2)*vers(2)
                  qph(k+n1+n2,3) = oph(k+n1+n2)*vers(3)

                  pcoll(0) = pcoll(0) - qph(k+n1+n2,0)
                  pcoll(1) = pcoll(1) - qph(k+n1+n2,1)
                  pcoll(2) = pcoll(2) - qph(k+n1+n2,2)
                  pcoll(3) = pcoll(3) - qph(k+n1+n2,3)
c                  en3guess = en3guess - oph(k+n1+n2)
               enddo
            endif
            if (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0) then
               ie = 1
c               print*,'1 -- '
               phsp = 0.d0
               w  = 0.d0
               return
            endif

            do k = 0,3
               pcoll(k) = pcoll(k) - dir3(k)
               pcollorig(k) = pcoll(k)
            enddo

            pcmodmu = 1.d0/sqrt(tridot(pcoll,pcoll))
            pcoll(0) = dir4(0)
            bcoll    = sqrt(1.d0 - m1**2/pcoll(0)**2)
            pcoll(1) = pcoll(1)*pcmodmu*bcoll*pcoll(0)
            pcoll(2) = pcoll(2)*pcmodmu*bcoll*pcoll(0)
            pcoll(3) = pcoll(3)*pcmodmu*bcoll*pcoll(0)

            if (n4.gt.0) then
               do k = 1,n4
                  iaccendi = 4
                  call getrnd(csi2,2)
                  call collinear(pflatc4,csi2(1),csi2(2),pcoll,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond) then
                    if (idarkon.eq.1) then
                       call photon_energy_dm(emaxgamma,en,wen)
                    else
                       cp = vers(3)/vers(0)
                       if (cp.gt.cpmax.or.cp.lt.cpmin) then
                            call photon_energy_new(ibich,pflate,expgg,
     .                              emaxgamma,esoft,ptld,vers,en,wen)
                       else
                            call photon_energy(ibich,pflate,expgg,
     .                                         emaxgamma,en,wen)
                       endif
                     endif
                     w = w * wen
                     oph(k+n1+n2+n3) = en
                  endif
***
                  wcoll4 = wcoll4 * wverg
                  qph(k+n1+n2+n3,0) = oph(k+n1+n2+n3)
                  qph(k+n1+n2+n3,1) = oph(k+n1+n2+n3)*vers(1)
                  qph(k+n1+n2+n3,2) = oph(k+n1+n2+n3)*vers(2)
                  qph(k+n1+n2+n3,3) = oph(k+n1+n2+n3)*vers(3)
               enddo
            endif
         else                   ! if iref = 4
            if (n4.gt.0) then
               do k = 1,n4
                  iaccendi = 4
                  call getrnd(csi2,2)
                  call collinear(pflatc4,csi2(1),csi2(2),dir4,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond) then
                    if (idarkon.eq.1) then
                       call photon_energy_dm(emaxgamma,en,wen)
                    else
                       cp = vers(3)/vers(0)
                       if (cp.gt.cpmax.or.cp.lt.cpmin) then
                            call photon_energy_new(ibich,pflate,expgg,
     .                              emaxgamma,esoft,ptld,vers,en,wen)
                       else
                            call photon_energy(ibich,pflate,expgg,
     .                                         emaxgamma,en,wen)
                       endif
                     endif
                     w = w * wen
                     oph(k+n1+n2+n3) = en
                  endif
***
                  wcoll4 = wcoll4 * wverg
                  qph(k+n1+n2+n3,0) = oph(k+n1+n2+n3)
                  qph(k+n1+n2+n3,1) = oph(k+n1+n2+n3)*vers(1)
                  qph(k+n1+n2+n3,2) = oph(k+n1+n2+n3)*vers(2)
                  qph(k+n1+n2+n3,3) = oph(k+n1+n2+n3)*vers(3)

                  pcoll(0) = pcoll(0) - qph(k+n1+n2+n3,0)
                  pcoll(1) = pcoll(1) - qph(k+n1+n2+n3,1)
                  pcoll(2) = pcoll(2) - qph(k+n1+n2+n3,2)
                  pcoll(3) = pcoll(3) - qph(k+n1+n2+n3,3)
c                  en4guess = en4guess - oph(k+n1+n2+n3)
               enddo
            endif
            if (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0) then
               ie = 1
c               print*,'2 -- '
               phsp = 0.d0
               w  = 0.d0
               return
            endif

            do k = 0,3
               pcoll(k) = pcoll(k) - dir4(k)
               pcollorig(k) = pcoll(k)
            enddo

            pcmodmu = 1.d0/sqrt(tridot(pcoll,pcoll))
            pcoll(0) = dir3(0)
            bcoll    = sqrt(1.d0 - m2**2/pcoll(0)**2)
            pcoll(1) = pcoll(1)*pcmodmu*bcoll*pcoll(0)
            pcoll(2) = pcoll(2)*pcmodmu*bcoll*pcoll(0)
            pcoll(3) = pcoll(3)*pcmodmu*bcoll*pcoll(0)

            if (n3.gt.0) then
               do k = 1,n3
                  iaccendi = 3
                  call getrnd(csi2,2)
                  call collinear(pflatc3,csi2(1),csi2(2),pcoll,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond) then
                    if (idarkon.eq.1) then
                       call photon_energy_dm(emaxgamma,en,wen)
                    else
                       cp = vers(3)/vers(0)
                       if (cp.gt.cpmax.or.cp.lt.cpmin) then
                            call photon_energy_new(ibich,pflate,expgg,
     .                              emaxgamma,esoft,ptld,vers,en,wen)
                       else
                            call photon_energy(ibich,pflate,expgg,
     .                                         emaxgamma,en,wen)
                       endif
                     endif
                     w = w * wen
                     oph(k+n1+n2) = en
                  endif
***
                  wcoll3 = wcoll3 * wverg
                  qph(k+n1+n2,0) = oph(k+n1+n2)
                  qph(k+n1+n2,1) = oph(k+n1+n2)*vers(1)
                  qph(k+n1+n2,2) = oph(k+n1+n2)*vers(2)
                  qph(k+n1+n2,3) = oph(k+n1+n2)*vers(3)
               enddo
            endif
         endif
      endif
***********************************
      do j = 0,3
         kone(j) = 0.d0
         do k = 1,nphot
            kone(j) = kone(j) + qph(k,j)
         enddo
         ptmp(j) = p1(j) + p2(j) - kone(j)
      enddo
***** TEST
      if (fs.eq.'gg') then
         koner = kone-qph(nphot,:) ! kone ridotto dell'ultimo fotone

         pk2 = dot(koner,koner)
         br  = koner(0) - pk2/2.d0/ei
         bq  = koner(0) - tridot(koner,vers)

         onm = ei*0.5 - br
         onm = onm/(1.d0 - bq/ei)

         a  = ei*0.5d0 - koner(0) + pk2/ei*0.5d0
         b  = 1.d0  - dot(koner,vers)/ei
         c  = ei - koner(0) + tridot(koner,dir3)/betafs/dir3(0)
         d  = 1.d0 - tridot(vers,dir3)/betafs/dir3(0)

         ep = (ei*a - esoft*c)/(ei*b-esoft*d)

c         ep = ei*0.5d0

         if (.not.(ep.lt.esoft.or.ep.ge.a/b.or.ep.eq.p1(0))) then
c         if (.FALSE.) then
! qui rimpiazza l'ultimo fotone generato, se non si verifica la condizione
! la condizione non dovrebbe presentarsi quasi mai, c'era beta invece di betafs sopra!!
            call photonenergygg(ep,a,b,esoft,en,wen)

            w = w * wen

            if (isnan(en).or.isnan(wen)) then
               print*,'!!!!!!!',ep,en,wen
               print*,'a b c d dir3',a,b,c,d,dir3
            endif
            wwen = 1.d0
            if (nphot.gt.1) then
               do k = 1,nphot-1
                  wwen = wwen * wten(k)
               enddo
            endif
            w = w * wwen

            qph(nphot,0) = en
            qph(nphot,1) = en*vers(1)
            qph(nphot,2) = en*vers(2)
            qph(nphot,3) = en*vers(3)
            do j = 0,3
               kone(j) = koner(j)+qph(nphot,j)
               ptmp(j) = p1(j) + p2(j) - kone(j)
            enddo
         else
! altrimenti si tiene quello gia' generato e si va avanti
c            print*,'OOOOPS in phasespace ep',ep,a/b,p1(0),
c     .           tridot(vers,dir3)/betafs/dir3(0),d,
c     .           ep-a/b
            do k = 1,nphot
               w = w * wten(k)
            enddo
         endif

      endif
***** TEST

      if (dot(ptmp,ptmp).lt.(m1+m2)**2.or.ptmp(0).lt.0.d0) then
         ie = 1
c         print*,'4 -- '
         phsp = 0.d0
         w  = 0.d0
         return
      endif
      amkone2 = dot(kone,kone)
***
      if (iref.eq.3) then
         vx  = sth*cphi
         vy  = sth*sphi
         vz  = cth
      else
         vx  = -sth*cphi
         vy  = -sth*sphi
         vz  = -cth
      endif
      o   = kone(0)
      ox  = kone(1)
      oy  = kone(2)
      oz  = kone(3)
      odv = ox*vx + oy*vy + oz*vz
      abig = ei**2 - 2.d0 * ei * o + amkone2 + m1**2 - m2**2
      bbig = -2.d0 * (ei - o)
      cbig = -2.d0 * odv
      arg  = abig**2*bbig**2-m1**2*bbig**2*(bbig**2-cbig**2)
      if (arg.lt.0.d0) then
         phsp = 0.d0
         w    = 0.d0
c         print*,'5 -- '
         ie = 2
         return
      endif
      p1mod_1 = abig*cbig + dsqrt(arg)
      p1mod_1 = p1mod_1 / (bbig**2 - cbig**2)
      p1mod_2 = abig*cbig - dsqrt(arg)
      p1mod_2 = p1mod_2 / (bbig**2 - cbig**2)
      p1jac = 1.d0
      if (p1mod_1.gt.0.d0.and.p1mod_2.gt.0.d0) then
         call getrnd(csi,1)
         if (csi(1).lt.0.5) p1mod = p1mod_1
         if (csi(1).ge.0.5) p1mod = p1mod_2
         p1jac = 2.d0
      endif
      if (p1mod_1.gt.0.d0.and.p1mod_2.le.0.d0) then
         p1mod = p1mod_1
      endif
      if (p1mod_1.le.0.d0.and.p1mod_2.gt.0.d0) then
         p1mod = p1mod_2
      endif
      if (p1mod_1.le.0.d0.and.p1mod_2.le.0.d0) then ! bestiale: con gg ci vuole .le.!!
         phsp = 0.d0
         w    = 0.d0
c         print*,'6 -- '
         ie   = 3
         return
      endif
      w = w * p1jac
      p3(0) = dsqrt(p1mod**2 + m1**2)
      p3(1) = p1mod * vx
      p3(2) = p1mod * vy
      p3(3) = p1mod * vz
      p4(0) =  ei - p3(0) - o
      p4(1) = -ox - p3(1)
      p4(2) = -oy - p3(2)
      p4(3) = -oz - p3(3)
***
      if (fs.eq.'gg'.and.(p3(0).lt.esoft.or.p4(0).lt.esoft)) then ! needed!
         phsp = 0.d0
         w    = 0.d0
         ie = 1
         return
      endif
***
!      phasespace!
      prodomega = 1.d0
      do k = 1,nphot
         prodomega = prodomega * qph(k,0)
      enddo
      den = dabs( p4(0)*p1mod + p3(0)*(p1mod + odv) )
      denphsp = den
      phsp = prodomega * p1mod**2 / denphsp
      phsp = phsp*duepigrechi/2.d0**(nphot+2)
      if (iref.eq.4) call exchange_mom(p3,p4)

***********************
! "regulator"... I try only on FINAL STATE (active only for ee and mm...)
      flatconst = 0.d0 !!!
      den       = 0.d0
      if (n3.gt.0.and.n4.gt.0) then
         prodp3 = 1.d0
         prodp4 = 1.d0
         do k = nis+1,nphot
            do i = 0,3
               ptmp(i) = qph(k,i)
            enddo
            p3q = dot(p3,ptmp)
            p4q = dot(p4,ptmp)
            if (k.le.(n1+n2+n3)) then
               prodp3=prodp3*
     .              (1.d0/p3q*p3(0)*ptmp(0))
            endif
            if (k.gt.(n1+n2+n3)) then
               prodp4=prodp4*
     .              (1.d0/p4q*p4(0)*ptmp(0))
            endif
         enddo
         if (iref.eq.3) anum = prodp3
         if (iref.eq.4) anum = prodp4
         regulator1 = anum/(prodp3+prodp4) * 2.d0
         w = w * regulator1
      endif
*******************
* second "regulator"
      den  = 1.d0
      anum = 1.d0
      iswitch = 1
      regulator2 = 1.d0
      if (fs.eq.'gg') iswitch = 0
      if (fs.ne.'gg') then
         if (n1.gt.0) then
            do k = 1,n1
               do i = 0,3
                  ptmp(i) = qph(k,i)
               enddo
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p1q*p1(0)*ptmp(0))
            enddo
         endif
         if (n2.gt.0) then
            do k = n1+1,n1+n2
               do i = 0,3
                  ptmp(i) = qph(k,i)
               enddo
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p2q*p2(0)*ptmp(0))
            enddo
         endif
         if (n3.gt.0) then
            do k = n1+n2+1,n1+n2+n3
               do i = 0,3
                  ptmp(i) = qph(k,i)
               enddo
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p3q*p3(0)*ptmp(0))
            enddo
         endif
         if (n4.gt.0) then
            do k = n1+n2+n3+1,nphot
               do i = 0,3
                  ptmp(i) = qph(k,i)
               enddo
               p1q = dot(p1,ptmp)
               p2q = dot(p2,ptmp)
               p3q = dot(p3,ptmp)
               p4q = dot(p4,ptmp)
               den = den * (
     .              1.d0/p1q*p1(0)*ptmp(0)+
     .              1.d0/p2q*p2(0)*ptmp(0)+
     .              iswitch/p3q*p3(0)*ptmp(0)+
     .              iswitch/p4q*p4(0)*ptmp(0))
               anum = anum*(1.d0/p4q*p4(0)*ptmp(0))
            enddo
         endif
         regulator2 = 1.d0* ncharged**nphot * anum/den
      endif
      w = w * regulator2
!-- extra regulator for gg ------------------------
      if (fs.eq.'gg') then
! I try a more symmetryc regulator ! test
         call regulate_gg(nphot,p1,p2,p3,p4,qph,reg)
         w=w/regulator2*reg
      endif
!-------------------------------------------------
      if (isnan(w)) then
         print*,''
         print*,'at the end of phasespace w = ',w
         print*,'wph,wen,wverg,p1jac,regulator1,regulator2,reg,wten'
         print*,wph,wen,wverg,p1jac,regulator1,regulator2,reg,
     .        wten(1:nphot)
         print*,''
      endif
      return
      end
************************************************************************************
      subroutine phasespaceTMP(p1,p2,p3,p4,qph,nphot,m1,m2,esoft,
     .     cth,w,phsp,ie)
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension ap3(0:3),ap4(0:3),pcollorig(0:3)
      dimension p4guess(0:3),p3guess(0:3),qone(0:3)
      dimension dir3(0:3),dir4(0:3),ptmp(0:3),vers(0:3)
      double precision kone(0:3),pcoll(0:3),ptld(0:3)
      double precision lambda
      logical lcond
      character*2 fs
      common/finalstate/fs
      integer emissionpattern(nphot)
      parameter (pi = 3.1415926535897932384626433832795029d0)
      parameter (imaxph = 40)
      double precision qph(imaxph,0:3),oph(imaxph),wten(nphot)
      real*8 csi(1),csi2(2),rnddummy
      common/for_debug/pcoll,denphsp,dir3,dir4,br,bq,onm
      common/various/beta,betafs
      common/channelref/iref
      common/iwide/iwider,ntot
      common/regulators/ich,idummy,regulator1,regulator2,wcoll3,wcoll4
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/coseni/cmax,cmin

      common/forshatsampling/den,omin,omax,ei
      common/cthcommon/cthc

      common/radpattern/nph(4)
      common/idebugging/idebug
      common/idarkon/idarkon

      common/weightwcos/wcos    ! to be depurated in phasespace if needed

      w      = 1.d0
      phsp   = 1.d0
      npart = nphot + 2
      nvar  = 3*nphot + 2

      ie     = 0

      duepigrechi = (2.d0*pi)**(-nvar)

      iwider = 0
      ntot   = nphot

      do j = 0,3
         do k = 1,imaxph
            qph(k,j) = 0.d0
         enddo
         nph(j+1) = 0
      enddo

      if (p1(0)+p2(0).lt.(m1+m2+nphot*esoft)) then
         phsp = 0.d0
         w    = 0.d0
         ie = 1
         return
      endif
*****************
!! I put this outside in order to raise efficiency
!      call get_cos_fer(rnd(1),c,wcos)
!      cth = c
!      w   = w * wcos
      c    = cth
! cthc in common for the sampling
      cthc = cth
*
      if (nphot.eq.0) then
         mass = p1(0)+p2(0)
         sqla = dsqrt(lambda(mass**2,m1**2,m2**2))
         phsp = phsp * sqla/8.d0/mass**2
         s = dsqrt(1.d0 - c**2)
         call getrnd(csi,1)
         call getphi(csi(1),phi,wph)
         w = w*wph
         sphi = sin(phi)
         cphi = cos(phi)
         p1mod = sqla/2.d0/mass
         p3(0) = dsqrt(p1mod**2 + m1**2)
         p3(1) = p1mod * sphi * s
         p3(2) = p1mod * cphi * s
         p3(3) = p1mod * c
         p4(0) = mass - p3(0)
         p4(1) = - p3(1)
         p4(2) = - p3(2)
         p4(3) = - p3(3)
         phsp = phsp * duepigrechi
         return
      endif
*******************
******************
      lcond = fs.eq.'ee'.or.fs.eq.'mm'.or.fs.eq.'mr'

      angle = 30.d0 *pi/180.d0

      angle = 0.d0 *pi/180.d0

      cpmax = cos(angle)
      cpmin = cos(pi - angle)
* the following for ng >= 1
      ncharged = 4
      if (fs.eq.'gg') ncharged = 2
      ir = 1

      call getrnd(csi,1)
      call getphi(csi(1),phi,wph)
      w    = w*wph
      sphi = dsin(phi)
      cphi = dcos(phi)
      ir   = ir + 2

      call get_patternNEW(ncharged,nphot,emissionpattern,wpat)
      do k = 1,nphot
         nph(emissionpattern(k)) = nph(emissionpattern(k)) + 1
      enddo
      w = w * wpat
      n1 = nph(1)
      n2 = nph(2)
      n3 = nph(3)
      n4 = nph(4)
      nis = n1 + n2
      nfs = n3 + n4

      ei = p1(0)+p2(0)

      ptld(0) = p1(0)+p2(0)
      vers(0)  = 0.d0
      do k = 1,3
         ptld(k)  = p1(k)+p2(k)
         vers(k)  = 0.d0
      enddo

      call getrnd(csi,1)
      iref = 3
      if (csi(1).gt.0.5) iref = 4
      wref = 1.d0
      if (nfs.gt.0) then
         if (n4.eq.0) iref = 3
         if (n3.eq.0) iref = 4
      endif
      w = w*wref
      if (fs.eq.'gg') iref = 3
****
** PHOTONS GENERATION...
      pflatc  = 0.d0
      pflatc3 = 0.d0
      pflatc4 = 0.d0
      pflate  = 0.d0
      expgg = 0.92d0
      ibich = 1
      do k = 1,nphot
         wen = 1.d0
         en  = esoft
         if (fs.eq.'gg') then
            call photon_energy(ibich,pflate,expgg,p1(0),
     .           en,wen)
         endif
         wten(k) = wen
         if (fs.eq.'ee'.or.fs.eq.'mm'.or.fs.eq.'mr'
     +       ) w = w * wen
         oph(k) = en
         ir = ir + 1
      enddo
** initial state photons

      emaxgamma = (p1(0)+p2(0))/2.d0
      omax = emaxgamma
      omin = esoft
      kone = 0
      prodomega = 1.d0
      if (nis.gt.0) then

         call getrnd(csi,1)
         nlast = csi(1)*nis + 1

         if (n1.gt.0) then
            do k = 1,n1
            IF (k.ne.nlast) THEN
               if (omax.lt.omin) then
                  phsp = 0.d0
                  w    = 0.d0
                  ie   = 1
                  return
               endif

               iaccendi = 1
               call getrnd(csi2,2)
               call collinear(0.d0,csi2(1),csi2(2),p1,vers,wverg)
               iaccendi = 0
               ir = ir + 2
               w = w * wverg
****
               if (lcond) then
                  if (idarkon.eq.1) then
                      call photon_energy_dm(emaxgamma,en,wen)
                  else
                     call photon_energy_IS(omin,omax,en,wen)
                  endif
                  w = w * wen
                  oph(k) = en
               endif
***
               qph(k,0) = oph(k)
               qph(k,1) = oph(k)*vers(1)
               qph(k,2) = oph(k)*vers(2)
               qph(k,3) = oph(k)*vers(3)

               kone = kone + qph(k,:)
               omax = omax - en
               prodomega = prodomega * en

            ENDIF
            enddo
         endif
         if (n2.gt.0) then
            do k = n1+1,n1+n2
            IF (k.ne.nlast) THEN
               if (omax.lt.omin) then
                  phsp = 0.d0
                  w    = 0.d0
                  ie   = 1
                  return
               endif

               iaccendi = 2
               call getrnd(csi2,2)
               call collinear(0.d0,csi2(1),csi2(2),p2,vers,wverg)
               iaccendi  = 0
               ir = ir + 2
               w = w * wverg
****
               if (lcond) then
                  if (idarkon.eq.1) then
                      call photon_energy_dm(emaxgamma,en,wen)
                  else
                     call photon_energy_IS(omin,omax,en,wen)
                  endif
                  w = w * wen
                  oph(k) = en
               endif
***
               qph(k,0) = oph(k)
               qph(k,1) = oph(k)*vers(1)
               qph(k,2) = oph(k)*vers(2)
               qph(k,3) = oph(k)*vers(3)

               kone = kone + qph(k,:)
               omax = omax - en

               prodomega = prodomega * en

            ENDIF
            enddo
         endif

*** last IS photon
         if (omax.lt.omin) then
            phsp = 0.d0
            w    = 0.d0
            ie   = 1
            return
         endif

         call getrnd(csi2,2)
         if (nlast.le.n1) then
            iaccendi = 1
            call collinear(0.d0,csi2(1),csi2(2),p1,vers,wverg)
            iaccendi = 0
         else
            iaccendi = 2
            call collinear(0.d0,csi2(1),csi2(2),p2,vers,wverg)
            iaccendi = 0
         endif
         w = w * wverg

         ptld   = p1 + p2 - kone
         ptld2  = dot(ptld,ptld)
         E1pE2  = p1(0)+p2(0)
         v1ltld = dot(vers,kone)
         den    = (E1pE2 - v1ltld)

         if (den.lt.0.d0) then
            phsp = 0.d0
            w    = 0.d0
            ie   = 1
            return
         endif
         smin = ptld2 - 2.d0*omax*den
! default
         smin = max(4.d0*m1*m1,smin)

! this takes into account the acollinearity cut
c         bmin = sqrt(1.d0 - m1*m1/emin/emin)
c         cccc = cos(pi-zmax)
c         smin = 2.d0*m1*m1 + 2.d0*emin*emin*(1.d0-bmin*bmin*cccc)
!!
         smax = ptld2 - 2.d0*omin*den

         if (smin.gt.smax) then
            phsp = 0.d0
            w    = 0.d0
            ie   = 1
            return
         endif

         call get_shat(ptld2,smin,smax,shat,wshat)
         w = w * wshat

         en   = 0.5d0*(ptld2 - shat)/den

         prodomega = prodomega*en*en/(ptld2-shat)

         k = nlast
         oph(k)   = en
         qph(k,0) = oph(k)
         qph(k,1) = oph(k)*vers(1)
         qph(k,2) = oph(k)*vers(2)
         qph(k,3) = oph(k)*vers(3)
 111     continue
         omax = omax - en
         kone = kone + qph(k,:)
      endif

      qone  = p1 + p2 - kone
      pcoll = qone
      qone2 = dot(qone,qone)

      if (lcond.and.
     .   (qone2.lt.(m1+m2+nfs*esoft)**2.or.pcoll(0).lt.0.d0)) then
         ie = 1
         phsp = 0.d0
         w  = 0.d0
         return
      endif


******
ccc      w = w * fakeBW(qone2)
******

*** phase space for inital state photons
      phspis = prodomega
***

***   if I'm really close to J/Psi, better to generate cth flat...
c      sqone2 = sqrt(qone2)
c      if (sqone2.gt.3.097d0.and.sqone2.lt.3.098d0) then
c         w = w/wcos
c         call getrnd(csi,1)
c         cth = (cmax - cmin) * csi(1) + cmin
c         w   = w * (cmax - cmin)
c      endif

! cos(th) and phi of electron already generated...
      sth  = dsqrt(1.d0 - cth**2)
      dir3(0) = ei/2d0
      dir3(1) = ei/2d0*betafs*cphi*sth
      dir3(2) = ei/2d0*betafs*sphi*sth
      dir3(3) = ei/2d0*betafs*cth
      dir4(0) = dir3(0)
      dir4(1)  = -dir3(1)
      dir4(2)  = -dir3(2)
      dir4(3)  = -dir3(3)

*******************************************

** final state photons
c      en3guess = dir3(0)
c      en4guess = dir4(0)
      wcoll3 = 1.d0
      wcoll4 = 1.d0
      if (nfs.gt.0) then
         if (iref.eq.3) then
            if (n3.gt.0) then
               do k = 1,n3
                  iaccendi = 3
                  call getrnd(csi2,2)
                  call collinear(pflatc3,csi2(1),csi2(2),dir3,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond) then
                    if (idarkon.eq.1) then
                       call photon_energy_dm(emaxgamma,en,wen)
                    else
                       cp = vers(3)/vers(0)
                       if (cp.gt.cpmax.or.cp.lt.cpmin) then
                            call photon_energy_new(ibich,pflate,expgg,
     .                              emaxgamma,esoft,ptld,vers,en,wen)
                       else
                            call photon_energy(ibich,pflate,expgg,
     .                                         emaxgamma,en,wen)
                       endif
                     endif
                     w = w * wen
                     oph(k+n1+n2) = en
                  endif
***
                  wcoll3 = wcoll3 * wverg
                  qph(k+n1+n2,0) = oph(k+n1+n2)
                  qph(k+n1+n2,1) = oph(k+n1+n2)*vers(1)
                  qph(k+n1+n2,2) = oph(k+n1+n2)*vers(2)
                  qph(k+n1+n2,3) = oph(k+n1+n2)*vers(3)

                  kone = kone + qph(k+n1+n2,:)

                  pcoll(0) = pcoll(0) - qph(k+n1+n2,0)
                  pcoll(1) = pcoll(1) - qph(k+n1+n2,1)
                  pcoll(2) = pcoll(2) - qph(k+n1+n2,2)
                  pcoll(3) = pcoll(3) - qph(k+n1+n2,3)
c                  en3guess = en3guess - oph(k+n1+n2)
               enddo
            endif
            if (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0) then
               ie = 1
               phsp = 0.d0
               w  = 0.d0
               return
            endif

            do k = 0,3
               pcoll(k) = pcoll(k) - dir3(k)
               pcollorig(k) = pcoll(k)
            enddo

            pcmodmu = 1.d0/sqrt(tridot(pcoll,pcoll))
            pcoll(0) = dir4(0)
            bcoll    = sqrt(1.d0 - m1**2/pcoll(0)**2)
            pcoll(1) = pcoll(1)*pcmodmu*bcoll*pcoll(0)
            pcoll(2) = pcoll(2)*pcmodmu*bcoll*pcoll(0)
            pcoll(3) = pcoll(3)*pcmodmu*bcoll*pcoll(0)

            if (n4.gt.0) then
               do k = 1,n4
                  iaccendi = 4
                  call getrnd(csi2,2)
                  call collinear(pflatc4,csi2(1),csi2(2),pcoll,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond) then
                    if (idarkon.eq.1) then
                       call photon_energy_dm(emaxgamma,en,wen)
                    else
                       cp = vers(3)/vers(0)
                       if (cp.gt.cpmax.or.cp.lt.cpmin) then
                            call photon_energy_new(ibich,pflate,expgg,
     .                              emaxgamma,esoft,ptld,vers,en,wen)
                       else
                            call photon_energy(ibich,pflate,expgg,
     .                                         emaxgamma,en,wen)
                       endif
                     endif
                     w = w * wen
                     oph(k+n1+n2+n3) = en
                  endif
***
                  wcoll4 = wcoll4 * wverg
                  qph(k+n1+n2+n3,0) = oph(k+n1+n2+n3)
                  qph(k+n1+n2+n3,1) = oph(k+n1+n2+n3)*vers(1)
                  qph(k+n1+n2+n3,2) = oph(k+n1+n2+n3)*vers(2)
                  qph(k+n1+n2+n3,3) = oph(k+n1+n2+n3)*vers(3)

                  kone = kone + qph(k+n1+n2+n3,:)

               enddo
            endif
         else                   ! if iref = 4
            if (n4.gt.0) then
               do k = 1,n4
                  iaccendi = 4
                  call getrnd(csi2,2)
                  call collinear(pflatc4,csi2(1),csi2(2),dir4,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond) then
                    if (idarkon.eq.1) then
                       call photon_energy_dm(emaxgamma,en,wen)
                    else
                       cp = vers(3)/vers(0)
                       if (cp.gt.cpmax.or.cp.lt.cpmin) then
                            call photon_energy_new(ibich,pflate,expgg,
     .                              emaxgamma,esoft,ptld,vers,en,wen)
                       else
                            call photon_energy(ibich,pflate,expgg,
     .                                         emaxgamma,en,wen)
                       endif
                     endif
                     w = w * wen
                     oph(k+n1+n2+n3) = en
                  endif
***
                  wcoll4 = wcoll4 * wverg
                  qph(k+n1+n2+n3,0) = oph(k+n1+n2+n3)
                  qph(k+n1+n2+n3,1) = oph(k+n1+n2+n3)*vers(1)
                  qph(k+n1+n2+n3,2) = oph(k+n1+n2+n3)*vers(2)
                  qph(k+n1+n2+n3,3) = oph(k+n1+n2+n3)*vers(3)

                  kone = kone + qph(k+n1+n2+n3,:)

                  pcoll(0) = pcoll(0) - qph(k+n1+n2+n3,0)
                  pcoll(1) = pcoll(1) - qph(k+n1+n2+n3,1)
                  pcoll(2) = pcoll(2) - qph(k+n1+n2+n3,2)
                  pcoll(3) = pcoll(3) - qph(k+n1+n2+n3,3)
c                  en4guess = en4guess - oph(k+n1+n2+n3)
               enddo
            endif
            if (dot(pcoll,pcoll).lt.(m1+m2)**2.or.pcoll(0).lt.0.d0) then
               ie = 1
               phsp = 0.d0
               w  = 0.d0
               return
            endif

            do k = 0,3
               pcoll(k) = pcoll(k) - dir4(k)
               pcollorig(k) = pcoll(k)
            enddo

            pcmodmu = 1.d0/sqrt(tridot(pcoll,pcoll))
            pcoll(0) = dir3(0)
            bcoll    = sqrt(1.d0 - m2**2/pcoll(0)**2)
            pcoll(1) = pcoll(1)*pcmodmu*bcoll*pcoll(0)
            pcoll(2) = pcoll(2)*pcmodmu*bcoll*pcoll(0)
            pcoll(3) = pcoll(3)*pcmodmu*bcoll*pcoll(0)

            if (n3.gt.0) then
               do k = 1,n3
                  iaccendi = 3
                  call getrnd(csi2,2)
                  call collinear(pflatc3,csi2(1),csi2(2),pcoll,vers,
     .                 wverg)
                  iaccendi = 0
                  ir = ir + 2
                  w = w * wverg
****
                  if (lcond) then
                    if (idarkon.eq.1) then
                       call photon_energy_dm(emaxgamma,en,wen)
                    else
                       cp = vers(3)/vers(0)
                       if (cp.gt.cpmax.or.cp.lt.cpmin) then
                            call photon_energy_new(ibich,pflate,expgg,
     .                              emaxgamma,esoft,ptld,vers,en,wen)
                       else
                            call photon_energy(ibich,pflate,expgg,
     .                                         emaxgamma,en,wen)
                       endif
                     endif
                     w = w * wen
                     oph(k+n1+n2) = en
                  endif
***
                  wcoll3 = wcoll3 * wverg
                  qph(k+n1+n2,0) = oph(k+n1+n2)
                  qph(k+n1+n2,1) = oph(k+n1+n2)*vers(1)
                  qph(k+n1+n2,2) = oph(k+n1+n2)*vers(2)
                  qph(k+n1+n2,3) = oph(k+n1+n2)*vers(3)

                  kone = kone + qph(k+n1+n2,:)

               enddo
            endif
         endif
      endif
***********************************
      do j = 0,3
         ptmp(j) = p1(j) + p2(j) - kone(j)
      enddo
***** TEST
      if (fs.eq.'gg') then
         wwen = 1.d0
         if (nphot.gt.1) then
            do k = 1,nphot-1
               wwen = wwen * wten(k)
            enddo
         endif
         w = w * wwen
         do j = 0,3
            kone(j) = kone(j)-qph(nphot,j)
         enddo

         pk2 = dot(kone,kone)
         br  = kone(0) - pk2/2.d0/ei
         bq  = kone(0) - tridot(kone,vers)
         onm = ei/2.d0 - br
         onm = onm/(1.d0 - bq/ei)
         a  = ei/2.d0 - kone(0) + pk2/2.d0/ei
         b  = 1.d0  - dot(kone,vers)/ei
         c  = ei - kone(0) + tridot(kone,dir3)/beta/dir3(0)
         d  = 1.d0 - tridot(vers,dir3)/beta/dir3(0)
         ep = (ei*a - esoft*c)/(ei*b-esoft*d)
         if (ep.lt.esoft.or.ep.ge.a/b.or.ep.eq.p1(0)) then
            w = 0.d0
            phsp = 0.d0
            ie = 1
            return
         endif
         call photonenergygg(ep,a,b,esoft,en,wen)
         w = w * wen
         qph(nphot,0) = en
         qph(nphot,1) = en*vers(1)
         qph(nphot,2) = en*vers(2)
         qph(nphot,3) = en*vers(3)
         do j = 0,3
            kone(j) = kone(j)+qph(nphot,j)
            ptmp(j) = p1(j) + p2(j) - kone(j)
         enddo
      endif
***** TEST
      if (dot(ptmp,ptmp).lt.(m1+m2)**2.or.ptmp(0).lt.0.d0) then
         ie = 1
         phsp = 0.d0
         w  = 0.d0
         return
      endif
      amkone2 = dot(kone,kone)
***
      if (iref.eq.3) then
         vx  = sth*cphi
         vy  = sth*sphi
         vz  = cth
      else
         vx  = -sth*cphi
         vy  = -sth*sphi
         vz  = -cth
      endif
      o   = kone(0)
      ox  = kone(1)
      oy  = kone(2)
      oz  = kone(3)
      odv = ox*vx + oy*vy + oz*vz
      abig = ei**2 - 2.d0 * ei * o + amkone2 + m1**2 - m2**2
      bbig = -2.d0 * (ei - o)
      cbig = -2.d0 * odv
      arg  = abig**2*bbig**2-m1**2*bbig**2*(bbig**2-cbig**2)
      if (arg.lt.0.d0) then
         phsp = 0.d0
         w    = 0.d0
         ie = 2
         return
      endif
      p1mod_1 = abig*cbig + dsqrt(arg)
      p1mod_1 = p1mod_1 / (bbig**2 - cbig**2)
      p1mod_2 = abig*cbig - dsqrt(arg)
      p1mod_2 = p1mod_2 / (bbig**2 - cbig**2)
      p1jac = 1.d0
      if (p1mod_1.gt.0.d0.and.p1mod_2.gt.0.d0) then
         call getrnd(csi,1)
         if (csi(1).lt.0.5) p1mod = p1mod_1
         if (csi(1).ge.0.5) p1mod = p1mod_2
         p1jac = 2.d0
      endif
      if (p1mod_1.gt.0.d0.and.p1mod_2.le.0.d0) then
         p1mod = p1mod_1
      endif
      if (p1mod_1.le.0.d0.and.p1mod_2.gt.0.d0) then
         p1mod = p1mod_2
      endif
      if (p1mod_1.le.0.d0.and.p1mod_2.le.0.d0) then ! bestiale: con gg ci vuole .le.!!
         phsp = 0.d0
         w    = 0.d0
         ie   = 3
         return
      endif
      w = w * p1jac
      p3(0) = dsqrt(p1mod**2 + m1**2)
      p3(1) = p1mod * vx
      p3(2) = p1mod * vy
      p3(3) = p1mod * vz
      p4(0) =  ei - p3(0) - o
      p4(1) = -ox - p3(1)
      p4(2) = -oy - p3(2)
      p4(3) = -oz - p3(3)
***
      if (fs.eq.'gg'.and.(p3(0).lt.esoft.or.p4(0).lt.esoft)) then
         phsp = 0.d0
         w    = 0.d0
         ie = 1
         return
      endif
***
!      phasespace!
      if (nfs.gt.0) then
         do k = nis+1,nis+nfs
            prodomega = prodomega * qph(k,0)
         enddo
      endif
      den = dabs( p4(0)*p1mod + p3(0)*(p1mod + odv) )
      denphsp = den
      phsp = prodomega * p1mod**2 / denphsp
      phsp = phsp*duepigrechi/2.d0**(nphot+2)
      if (iref.eq.4) call exchange_mom(p3,p4)

***********************
! "regulator"... I try only on FINAL STATE (active only for ee and mm...)
      flatconst = 0.d0 !!!
      den       = 0.d0
      if (n3.gt.0.and.n4.gt.0) then
         prodp3 = 1.d0
         prodp4 = 1.d0
         do k = nis+1,nphot
            do i = 0,3
               ptmp(i) = qph(k,i)
            enddo
            p3q = dot(p3,ptmp)
            p4q = dot(p4,ptmp)
            if (k.le.(n1+n2+n3)) then
               prodp3=prodp3*
     .              (1.d0/p3q*p3(0)*ptmp(0))
            endif
            if (k.gt.(n1+n2+n3)) then
               prodp4=prodp4*
     .              (1.d0/p4q*p4(0)*ptmp(0))
            endif
         enddo
         if (iref.eq.3) anum = prodp3
         if (iref.eq.4) anum = prodp4
         regulator1 = anum/(prodp3+prodp4) * 2.d0
         w = w * regulator1
      endif
*******************
* second "regulator"
      den  = 1.d0
      anum = 1.d0

***   regulator(s) ***
      IF (lcond) THEN
      if (n1.gt.0) then
         do k = 1,n1
            anum = anum /dot(p1,qph(k,:)) !* p1(0)*qph(k,0)
         enddo
      endif
      if (n2.gt.0) then
         do k = n1+1,n1+n2
            anum = anum /dot(p2,qph(k,:)) !* p2(0)*qph(k,0)
         enddo
      endif
      if (n3.gt.0) then
         do k = n1+n2+1,n1+n2+n3
            anum = anum /dot(p3,qph(k,:)) !* p3(0)*qph(k,0)
         enddo
      endif
      if (n4.gt.0) then
         do k = n1+n2+n3+1,n1+n2+n3+n4
            anum = anum /dot(p4,qph(k,:)) !* p4(0)*qph(k,0)
         enddo
      endif

cc      anum = anum*fakeBW(qone2)
cc      anum = anum * sdif(qone2,cthc)
cc      anum = anum * spezzasdif(qone2)

      call handregulator(nphot,p1,p2,p3,p4,qph(1:nphot,:),den)

      else

      iswitch = 1
      if (fs.eq.'gg') iswitch = 0
      if (n1.gt.0) then
         do k = 1,n1
            do i = 0,3
               ptmp(i) = qph(k,i)
            enddo
            p1q = dot(p1,ptmp)
            p2q = dot(p2,ptmp)
            p3q = dot(p3,ptmp)
            p4q = dot(p4,ptmp)
            den = den * (
     .           1.d0/p1q*p1(0)*ptmp(0)+
     .           1.d0/p2q*p2(0)*ptmp(0)+
     .           iswitch/p3q*p3(0)*ptmp(0)+
     .           iswitch/p4q*p4(0)*ptmp(0))
            anum = anum*(1.d0/p1q*p1(0)*ptmp(0))
         enddo
      endif
      if (n2.gt.0) then
         do k = n1+1,n1+n2
            do i = 0,3
               ptmp(i) = qph(k,i)
            enddo
            p1q = dot(p1,ptmp)
            p2q = dot(p2,ptmp)
            p3q = dot(p3,ptmp)
            p4q = dot(p4,ptmp)
            den = den * (
     .           1.d0/p1q*p1(0)*ptmp(0)+
     .           1.d0/p2q*p2(0)*ptmp(0)+
     .           iswitch/p3q*p3(0)*ptmp(0)+
     .           iswitch/p4q*p4(0)*ptmp(0))
            anum = anum*(1.d0/p2q*p2(0)*ptmp(0))
         enddo
      endif
      if (n3.gt.0) then
         do k = n1+n2+1,n1+n2+n3
            do i = 0,3
               ptmp(i) = qph(k,i)
            enddo
            p1q = dot(p1,ptmp)
            p2q = dot(p2,ptmp)
            p3q = dot(p3,ptmp)
            p4q = dot(p4,ptmp)
            den = den * (
     .           1.d0/p1q*p1(0)*ptmp(0)+
     .           1.d0/p2q*p2(0)*ptmp(0)+
     .           iswitch/p3q*p3(0)*ptmp(0)+
     .           iswitch/p4q*p4(0)*ptmp(0))
            anum = anum*(1.d0/p3q*p3(0)*ptmp(0))
         enddo
      endif
      if (n4.gt.0) then
         do k = n1+n2+n3+1,nphot
            do i = 0,3
               ptmp(i) = qph(k,i)
            enddo
            p1q = dot(p1,ptmp)
            p2q = dot(p2,ptmp)
            p3q = dot(p3,ptmp)
            p4q = dot(p4,ptmp)
            den = den * (
     .           1.d0/p1q*p1(0)*ptmp(0)+
     .           1.d0/p2q*p2(0)*ptmp(0)+
     .           iswitch/p3q*p3(0)*ptmp(0)+
     .           iswitch/p4q*p4(0)*ptmp(0))
            anum = anum*(1.d0/p4q*p4(0)*ptmp(0))
         enddo
      endif

 321  continue
      ENDIF

      regulator2 = 1.d0 * anum/den
      w = w * regulator2
!-- extra regulator for gg ------------------------
      if (fs.eq.'gg') then
! I try a more symmetryc regulator ! test
         call regulate_gg(nphot,p1,p2,p3,p4,qph,reg)
         w=w/regulator2*reg
      endif
!-------------------------------------------------
      return
      end
************************
      subroutine phasespaceNEW(p1,p2,p3,p4,qph,nphot,m1,m2,esoft,
     .     cth,w,phsp,ie)
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension ap3(0:3),ap4(0:3),pcollorig(0:3)
      dimension p4guess(0:3),p3guess(0:3),direm(0:3)
      dimension dir3(0:3),dir4(0:3),ptmp(0:3),vers(0:3)
      double precision kone(0:3),pcoll(0:3),ptld(0:3)
      double precision qone(0:3),konefs(0:3),pboo1(0:3),pboo2(0:3)
      double precision lambda
      logical lcond
      character*2 fs
      common/finalstate/fs
      integer emissionpattern(nphot)
      parameter (pi = 3.1415926535897932384626433832795029d0)
      parameter (imaxph = 40)
      double precision qph(imaxph,0:3),oph(imaxph),wten(nphot)
      double precision qphfs(imaxph,0:3),ophfs(imaxph)
      real*8 csi(1),csi2(2),rnddummy
      common/for_debug/pcoll,denphsp,dir3,dir4,br,bq,onm
      common/various/beta,betafs
      common/channelref/iref
      common/iwide/iwider,ntot
      common/regulators/ich,idummy,regulator1,regulator2,wcoll3,wcoll4
      common/radpattern/nph(4)
      common/idebugging/idebug
      common/tmintmax/tmin,tmax,tcur
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax

      common/epssoft/eps

      common/weightwcos/wcos ! to be depurated in phasespaceNEW

      w      = 1.d0
      npart = nphot + 2
      nvar  = 3*nphot + 2

      ie     = 0

      duepigrechi = (2.d0*pi)**(-nvar)
      phsp   = 1.d0

      iwider = 0
      ntot   = nphot

      do j = 0,3
         do k = 1,imaxph
            qph(k,j) = 0.d0
         enddo
         nph(j+1) = 0
      enddo

      if (p1(0)+p2(0).lt.(m1+m2+nphot*esoft)) then
         phsp = 0.d0
         w    = 0.d0
         ie = 1
         return
      endif
*****************
!! I put this outside in order to raise efficiency
!      call get_cos_fer(rnd(1),c,wcos)
!      cth = c
!      w   = w * wcos
      c = cth
*
      if (nphot.eq.0) then
         mass = p1(0)+p2(0)
         sqla = dsqrt(lambda(mass**2,m1**2,m2**2))
         phsp = phsp * sqla/8.d0/mass**2
         s = dsqrt(1.d0 - c**2)
         call getrnd(csi,1)
         call getphi(csi(1),phi,wph)
         w = w*wph
         sphi = sin(phi)
         cphi = cos(phi)
         p1mod = sqla/2.d0/mass
         p3(0) = dsqrt(p1mod**2 + m1**2)
         p3(1) = p1mod * sphi * s
         p3(2) = p1mod * cphi * s
         p3(3) = p1mod * c
         p4(0) = mass - p3(0)
         p4(1) = - p3(1)
         p4(2) = - p3(2)
         p4(3) = - p3(3)
         phsp  = phsp*duepigrechi
         return
      endif
*******************
      lcond = fs.eq.'ee'.or.fs.eq.'mm'.or.fs.eq.'mr'
! here only if lcond = .true. !!!

* the following for ng >= 1
      ncharged = 4
      ir = 1

! cos(th) of electron already generated...
      sth  = dsqrt(1.d0 - cth**2)
      call getrnd(csi,1)
      call getphi(csi(1),phi,wph)
      w    = w*wph
      sphi = dsin(phi)
      cphi = dcos(phi)
      ir   = ir + 2
***
      direm(0) = 1.d0
      direm(1) = sphi*sth
      direm(2) = cphi*sth
      direm(3) = cth
***
      call get_pattern(ncharged,nphot,emissionpattern)
      do k = 1,nphot
         nph(emissionpattern(k)) = nph(emissionpattern(k)) + 1
      enddo
      n1 = nph(1)
      n2 = nph(2)
      n3 = nph(3)
      n4 = nph(4)
      nis = n1 + n2
      nfs = n3 + n4
      ei = p1(0)+p2(0)
****
** PHOTONS GENERATION...
      omin = eps * ei * 0.5d0
      omax = ei  * 0.5d0

** initial state photons
      kone = 0
      prodomega = 1.d0
      if (nis.gt.0) then
***
         call getrnd(csi,1)
         nlast = csi(1)*nis + 1
c
cc         nlast = nis + 1
cc
***
         if (n1.gt.0) then
            do k = 1,n1
               if (k.ne.nlast) then
                  if (omax.lt.omin) then
                     phsp = 0.d0
                     w    = 0.d0
                     ie   = 1
                     return
                  endif

                  call getrnd(csi2,2)
                  call collinear_simple(0.d0,csi2(1),csi2(2),p1,
     .                 vers,wverg)
                  w = w * wverg
                  call photon_energy_IS(omin,omax,en,wen)
                  w = w * wen
                  oph(k)  = en
                  qph(k,0)= oph(k)
                  qph(k,1)= oph(k)*vers(1)
                  qph(k,2)= oph(k)*vers(2)
                  qph(k,3)= oph(k)*vers(3)

                  omax = omax - en
                  prodomega = prodomega * 0.5 * en
                  kone = kone + qph(k,:)
               endif
            enddo
         endif
         if (n2.gt.0) then
            do k = n1+1,n1+n2
               if (k.ne.nlast) then
                  if (omax.lt.omin) then
                     phsp = 0.d0
                     w    = 0.d0
                     ie = 1
                     return
                  endif
                  call getrnd(csi2,2)
                  call collinear_simple(0.d0,csi2(1),csi2(2),p2,
     .                 vers,wverg)
                  w = w * wverg
                  call photon_energy_IS(omin,omax,en,wen)
                  w = w * wen
                  oph(k)  = en
                  qph(k,0)= oph(k)
                  qph(k,1)= oph(k)*vers(1)
                  qph(k,2)= oph(k)*vers(2)
                  qph(k,3)= oph(k)*vers(3)

                  omax = omax - en
                  prodomega = prodomega * 0.5 * en
                  kone = kone + qph(k,:)
               endif
            enddo
         endif
*** last IS photon
cc         goto 111

         if (omax.lt.omin) then
            phsp = 0.d0
            w    = 0.d0
            ie = 1
            return
         endif
         call getrnd(csi2,2)
         if (nlast.le.n1) then
            call collinear_simple(0.d0,csi2(1),csi2(2),p1,vers,wverg)
         else
            call collinear_simple(0.d0,csi2(1),csi2(2),p2,vers,wverg)
         endif
         w = w * wverg

         ptld   = p1 + p2 - kone
         ptld2  = dot(ptld,ptld)
         E1pE2  = p1(0)+p2(0)
         v1ltld = dot(vers,kone)
         den    = (E1pE2 - v1ltld)

         if (den.lt.0.d0) then
            phsp = 0.d0
            w    = 0.d0
            ie   = 1
            return
         endif
         smin = ptld2 - 2.d0*omax*den
         smin = max(4.d0*m1*m1,smin)
c         smin = 4.d0*emin**2

         smax = ptld2 - 2.d0*omin*den

         if (smin.gt.smax) then
            phsp = 0.d0
            w    = 0.d0
            ie   = 1
            return
         endif

         call get_shat(ptld2,smin,smax,shat,wshat)
         w = w * wshat

         en   = 0.5d0*(ptld2 - shat)/den

         prodomega = 0.5d0*prodomega*en*en/(ptld2-shat)

         k = nlast
         oph(k)   = en
         qph(k,0) = oph(k)
         qph(k,1) = oph(k)*vers(1)
         qph(k,2) = oph(k)*vers(2)
         qph(k,3) = oph(k)*vers(3)
 111     continue
         omax = omax - en
         kone = kone + qph(k,:)
      endif

      qone  = p1 + p2 - kone
      qone2 = dot(qone,qone)

      if (qone(0).lt.0.d0.or.qone2.lt.(m1+m2+nfs*omin)**2) then
         phsp = 0.d0
         w    = 0.d0
         ie   = 1
         return
      endif

*** phase space for inital state photons
      phspis = prodomega
***
      betaq  = sqrt(1.d0-qone2/qone(0)**2)
      gammaq = 1.d0/sqrt(1.d0-betaq**2)

cc      print*,(1.d0-betaq)*gammaq,sqrt((1.d0-betaq)/(1.d0+betaq))

!! I think now it's better to generate fs momenta in the CRF where p3 is along z, then I boost back where qone is at rest
!! and I rotate everything back

      dir3(0) =  sqrt(qone2)*0.5
      dir3(1) =  0.d0
      dir3(2) =  0.d0
      dir3(3) =  sqrt((1.d0 - m1/dir3(0))*(1.d0 + m1/dir3(0)))*dir3(0)
      dir4(0) =  dir3(0)
      dir4(1) =  0.d0
      dir4(2) =  0.d0
      dir4(3) = -dir3(3) ! m1 = m2

      dir3(1) = dir3(3)*direm(1)
      dir3(2) = dir3(3)*direm(2)
      dir3(3) = dir3(3)*direm(3)
      dir4    = -dir3
      dir4(0) = -dir4(0)

** omin for FS photons, omin it was for IS photons
!  this is really a not trivial task, because ominfs is correct in the frame where q is at rest
!  but not where p3+p4 is at rest and where I generate photons! This affects also Hto4f. Think better.
!  ** Update: I solved this by generating the energy from 0 for FS photons, see the routine
      ominfs = gammaq * omin * ( 1.d0 - betaq)

c      ominfs = omin

!     nota gamma * (1-b) = (1-b)/sqrt((1-b)*(1+b)) = sqrt((1-b)/(1+b))

c      ominfs = ominfs * 0.5d0

c      btest  = ominfs/sqrt(qone2+ominfs**2)
c      bah = (dir3(0)*0.02d0)**2
c      btest = (1.d0-ominfs**2/bah)/(1.d0+ominfs**2/bah)
c      ominfs = ominfs * sqrt((1.d0-btest)/(1.d0+btest))
c      print*,sqrt((1.d0-btest)/(1.d0+btest))

      omaxfs = dir3(0)          ! * gammaq * (1.d0+betaq)

** final state photons
      konefs      = 0.d0
      prodomegafs = 1.d0
      if (n3.gt.0) then
         do k = 1,n3
            if (omaxfs.lt.ominfs) then
               phsp = 0.d0
               w    = 0.d0
               ie = 1
               return
            endif

            call getrnd(csi2,2)
            call collinear_simple(0.d0,csi2(1),csi2(2),dir3,vers,wverg)
            w = w * wverg
            call photon_energy_FS(ominfs,omaxfs,en,wen)
            w = w * wen
            ophfs(k)   = en
            qphfs(k,0) = ophfs(k)
            qphfs(k,1) = ophfs(k)*vers(1)
            qphfs(k,2) = ophfs(k)*vers(2)
            qphfs(k,3) = ophfs(k)*vers(3)

            do i = 0,3
               konefs(i) = konefs(i) + qphfs(k,i)
            enddo
            omaxfs      = omaxfs - en
            prodomegafs = prodomegafs *0.5*en
         enddo
      endif
      if (n4.gt.0) then
         do k = n3+1,n3+n4
            if (omaxfs.lt.ominfs) then
               phsp = 0.d0
               w    = 0.d0
               ie = 1
               return
            endif

            call getrnd(csi2,2)
            call collinear_simple(0.d0,csi2(1),csi2(2),dir4,vers,wverg)
            w = w * wverg
            call photon_energy_FS(ominfs,omaxfs,en,wen)
            w = w * wen
            ophfs(k)   = en
            qphfs(k,0) = ophfs(k)
            qphfs(k,1) = ophfs(k)*vers(1)
            qphfs(k,2) = ophfs(k)*vers(2)
            qphfs(k,3) = ophfs(k)*vers(3)

            do i = 0,3
               konefs(i) = konefs(i) + qphfs(k,i)
            enddo

            omaxfs      = omaxfs - en
            prodomegafs = prodomegafs *0.5*en
         enddo
      endif

      E3 = (sqrt(qone2+tridot(konefs,konefs))-konefs(0))*0.5d0
      if (E3.lt.m1) then
         ie = 1
         w  = 0.d0
         phsp = 0.d0
         return
      endif

      bfs   = sqrt((1.d0-m1/E3)*(1.d0+m1/E3))

      p3(0) = E3
      p3(1) = E3*bfs*direm(1) !0.d0
      p3(2) = E3*bfs*direm(2)  !0.d0
      p3(3) = bfs*E3*direm(3)

      p4    = -p3
      p4(0) = -p4(0)

      sp     = 4.d0*E3*E3 ! = dot(p3+p4,p3+p4)
      phspfs = 0.25d0*bfs *0.5d0*sp/qone2 *
     .         prodomegafs*1.d0/(1.d0+konefs(0)/sqrt(sp))

      phsp = duepigrechi*phspis*phspfs

**** at this stage p3+p4 is at rest, now I go where p3+p4+konefs (i.e. qone) is at rest
      ptld = p3+p4+konefs

      call new_boost(ptld,p3,p3,1)
      call new_boost(ptld,p4,p4,1)
      if (nfs.gt.0) then
         do k = 1,nfs
            dir3 = qphfs(k,:)
            call new_boost(ptld,dir3,dir3,1)
            qphfs(k,:) = dir3
         enddo
      endif

** and now I boost in the lab, where qone has its components
      call new_boost(qone,p3,p3,-1)
      call new_boost(qone,p4,p4,-1)
      if (nfs.gt.0) then
         do k = 1,nfs
            dir3 = qphfs(k,:)
            call new_boost(qone,dir3,dir3,-1)
            qphfs(k,:) = dir3
!! checking in the lab that the FS photons are > omin
            if (qphfs(k,0).lt.omin) then
               ie   = 1
               w    = 0.d0
               phsp = 0.d0
               return
            endif
!!
         enddo
      endif
****
      if (nfs.gt.0) then
         do k = 1,nfs
            qph(nis+k,:) = qphfs(k,:)
            kone = kone + qphfs(k,:)
         enddo
      endif

***   regulator(s) ***
      anum = 1.d0
      if (n1.gt.0) then
         do k = 1,n1
            anum = anum /dot(p1,qph(k,:)) !* p1(0)*qph(k,0)
         enddo
      endif
      if (n2.gt.0) then
         do k = n1+1,n1+n2
            anum = anum /dot(p2,qph(k,:)) !* p2(0)*qph(k,0)
         enddo
      endif
      if (n3.gt.0) then
         do k = n1+n2+1,n1+n2+n3
            anum = anum /dot(p3,qph(k,:)) !* p3(0)*qph(k,0)
         enddo
      endif
      if (n4.gt.0) then
         do k = n1+n2+n3+1,n1+n2+n3+n4
            anum = anum /dot(p4,qph(k,:)) !* p4(0)*qph(k,0)
         enddo
      endif

      call handregulator(nphot,p1,p2,p3,p4,qph(1:nphot,:),den)
c$$$      den = 1.d0
c$$$      do k = 1,nphot
c$$$         den = den * (
c$$$     .         1.d0/dot(p1,qph(k,:)) !* p1(0)*qph(k,0)
c$$$     .       + 1.d0/dot(p2,qph(k,:)) !* p2(0)*qph(k,0)
c$$$     .       + 1.d0/dot(p3,qph(k,:)) !* p3(0)*qph(k,0)
c$$$     .       + 1.d0/dot(p4,qph(k,:)) !* p4(0)*qph(k,0)
c$$$     .        )
c$$$      enddo

      regulator = anum/den * 4**nphot
      w         = w * regulator

      return
      end
*************************************************************************************
**********************************************************
      subroutine handregulator(ng,p1,p2,p3,p4,qph,den)
      implicit double precision (a-h,m,o-z)

      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension qph(ng,0:3),q1(0:3),q2(0:3),q3(0:3)
      dimension p12(0:3),p34(0:3),pz1(0:3),pz2(0:3)

      dimension p14(0:3),p23(0:3),pz12(0:3),pz34(0:3)
      dimension pz14(0:3),pz23(0:3)

      dimension p12q1(0:3),p12q2(0:3)
      dimension p12q1q2(0:3)
      dimension p34q1(0:3),p34q2(0:3)
      dimension p34q1q2(0:3)
      dimension coll1(ng),coll2(ng),coll3(ng),coll4(ng)
      dimension collt(4,ng)
      dimension pt(0:3,4)
      integer   pow4(0:ng)
      integer   pow2(0:ng)

      common/cthcommon/cthc

      pow2(0) = 1
      do k = 1,ng
         pow2(k)    = pow2(k-1)*2
         collt(1,k) = 1.d0/dot(p1,qph(k,:))! * p1(0)*qph(0,k)
         collt(2,k) = 1.d0/dot(p2,qph(k,:))! * p2(0)*qph(0,k)
         collt(3,k) = 1.d0/dot(p3,qph(k,:))! * p3(0)*qph(0,k)
         collt(4,k) = 1.d0/dot(p4,qph(k,:))! * p4(0)*qph(0,k)
      enddo

* OPTIMIZED
      p12 = p1+p2
c      p34 = p3+p4 ! not needed here
      den = 0.d0
      do j = 0,2**ng-1
         pz12 = p12
c         pz34 = p34 ! not needed here
         prod = 1.d0
         do k = 1,ng
            iz = mod(j,pow2(k)) / pow2(k-1) + 1
            if (iz.eq.1) then
c               pz12  = pz12 + qph(k,:)
               pz12  = pz12 - qph(k,:) ! here I need the -...
               prod  = prod * (collt(1,k)+collt(2,k))
            else
c               pz34  = pz34 + qph(k,:) ! not needed here
               prod = prod * (collt(3,k)+collt(4,k))
            endif
         enddo
c     prod = prod*bw(dot(pz12,pz12))*bw(dot(pz34,pz34)) ! to be modified for bhabha
cc         prod = prod * fakeBW(dot(pz12,pz12))
cc         prod = prod * sdif(dot(pz12,pz12),cthc)
cc         prod = prod * spezzasdif(dot(pz12,pz12))
         den  = den + prod
      enddo

      return

      end

***************************************************************************
      !!======================================!!
      !! SAMPLING CHANNELS for e+e- -> X+X- g !!
      !!======================================!!

      subroutine CHN_1g_ISR(p1,p2,p3,m3,p4,m4,qph,
     . w,norm,phsp,icut,sampling,angdip,igen)
***************************************************************************
      !! m_34^2 sampled as 1/[(m_34^2 - M^2)^2 + (G*M)^2]
      !! cos_1 of the photon q is sampled as 1/(1-b^2 cos^2)
      implicit none
      real*8 p1,p2,p3,p4,m3,m4,qph
      integer icut,igen
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),qph(40,0:3)
      real*8 q
      dimension q(0:3)
      real*8 w,phsp
      real*8 M,G,Q2_34_min,Q2_34_max,Q2_34,norm_34,w_34
      character*2 sampling,angdip
      real*8 beta_e
      real*8 ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi

      real*8 w_smpl,norm,w_flat
      real*8 cmin,cmax
      common/coseni/cmax,cmin

      real*8 cos_q,w_cos_q,norm_cos_q
      real*8 phi_q,w_phi_q

      real*8 csi(3),xi(1)
      external getrnd
      real*8 dot
      external dot

      real*8 Q34,m34,s
      dimension Q34(0:3)


      real*8 phsp_2_34_q
      real*8 phi3,w_phi3,cos3,cos3_min,cos3_max,phsp_2_3_4
      real*8 thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
!! parameters of rho resonance
      M = 0.77456d0
      G = 0.14832d0

      !! initialisation of weights and norms
      w= 1d0
      w_smpl = 1d0
      norm = 1d0
      w_flat = 1d0
      phsp = 1d0

      if(igen.eq.0) then
         q(0:3) = qph(1,0:3)
         Q2_34 = dot(p3+p4,p3+p4)
         Q34= p3+p4
         cos_q = qph(1,3)/qph(1,0)

      endif

      s = dot(p1+p2,p1+p2)
      Q2_34_min = (m3+m4)**2
      Q2_34_max = s-2d0*sqrt(s)*egmin

      if(sampling.eq.'BW') then
         call sampling_BW_resonances(M,G,Q2_34_min,Q2_34_max,Q2_34,
     >        w_34,norm_34,igen)
      elseif(sampling.eq.'CO') then
         call sampling_energy_mxx(Q2_34_min,Q2_34_max,Q2_34,w_34,
     >    norm_34,igen)
      else if(sampling.eq.'SC') then
         call sampling_prop(s,Q2_34_min,Q2_34_max,Q2_34,w_34,
     >    norm_34,igen)
      else if(sampling.eq.'FL') then
         if(igen.eq.1) then
            call getrnd(xi,1)
            Q2_34 = Q2_34_min + xi(1)*(Q2_34_max-Q2_34_min)
         endif

         w_34=1.d0
         norm_34=1.d0

         w_flat = w_flat*(Q2_34_max-Q2_34_min)
      endif

      w_smpl = w_smpl * w_34
      norm   = norm * norm_34

      beta_e = sqrt(1d0-ame**2/p1(0)**2)

      if(angdip.eq.'co') then
         call sampling_cos_collsq(beta_e,1d0,-1d0,1d0,
     .        cos_q,w_cos_q,norm_cos_q,igen)

      else if(angdip.eq.'fl') then
         if(igen.eq.1) then
            call getrnd(xi,1)
            cos_q = -1.d0 + 2.d0*xi(1)
         endif
         w_cos_q = 1.d0
         norm_cos_q = 1.d0
         w_flat  = w_flat *  2.d0

      endif

      w_smpl = w_smpl * w_cos_q
      norm = norm*norm_cos_q
      w=w*w_smpl

      if(igen.eq.1) then
         call getrnd(csi,3)

         call getphi(csi(1),phi_q,w_phi_q)
         w_flat = w_flat * w_phi_q

         m34 = sqrt(Q2_34)
         call dLips_2(p1+p2,cos_q,phi_q,q,0d0,Q34,m34,phsp_2_34_q,icut)

         call getphi(csi(2),phi3,w_phi3)
         w_flat = w_flat * w_phi3

         cos3_min = -1d0
         cos3_max = 1d0
         cos3 = cos3_min + (cos3_max - cos3_min)*csi(3)
         w_flat = w_flat *(cos3_max-cos3_min)

         call dLips_2(Q34,cos3,phi3,p3,m3,p4,m4,phsp_2_3_4,icut)

         phsp = phsp * phsp_2_34_q * phsp_2_3_4
         phsp = phsp * (2d0*pi)**7
         phsp = phsp*w_flat
         qph(1,0:3)=q(0:3)
      endif
         if(icut.ne.0) then
         phsp=0d0
      endif
      return
      end


      subroutine CHN_1g_BW_ISR(p1,p2,p3,m3,p4,m4,qph,w,
     . norm,phsp,icut,igen)
***************************************************************************
      !! m_34^2 sampled as 1/[(m_34^2 - M^2)^2 + (G*M)^2]
      !! cos_1 of the photon q is sampled as FLAT
      implicit none
      real*8 p1,p2,p3,p4,m3,m4,q
      integer icut,igen
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      real*8 qph
      dimension qph(40,0:3)
      real*8 w,phsp
      real*8 M,G,Q2_34_min,Q2_34_max,Q2_34,norm_34,w_34

      real*8 beta_e
      real*8 ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi

      real*8 w_smpl,norm,w_flat
      real*8 cmin,cmax
      common/coseni/cmax,cmin

      real*8 cos_q,w_cos_q,norm_cos_q
      real*8 phi_q,w_phi_q

      real*8 csi(4)
      external getrnd
      real*8 dot
      external dot

      real*8 Q34,m34
      dimension Q34(0:3)


      real*8 cos_q_min,cos_q_max


      real*8 phsp_2_34_q
      real*8 phi3,w_phi3,cos3,cos3_min,cos3_max,phsp_2_3_4
      !! parameters of rho resonance
      M = 0.77456d0
      G = 0.14832d0


      if(igen.eq.0) then
         q(0:3) = qph(1,0:3)
         Q2_34 = dot(p3+p4,p3+p4)
      endif

      !! initialisation of weights and norms
      w= 1d0
      w_smpl = 1d0
      norm = 1d0
      w_flat = 1d0
      phsp = 1d0

      Q2_34_min = (m3+m4)**2
      Q2_34_max = dot(p1+p2,p1+p2)
      call sampling_BW_resonances(M,G,Q2_34_min,Q2_34_max,Q2_34,w_34,
     >     norm_34,igen)
      w_smpl = w_smpl * w_34
      norm   = norm * norm_34

   !    beta_e = sqrt(1d0-ame**2/p1(0)**2)
   !    call sampling_cos_collsq(beta_e,1d0,-1d0,1d0,
   !   .                               cos_q,w_cos_q,norm_cos_q,igen)
   !    w_smpl = w_smpl * w_cos_q
   !    norm = norm*norm_cos_q

      w=w*w_smpl

      if(igen.eq.1) then
         call getrnd(csi,4)
         cos_q_min=-1d0
         cos_q_max=1d0
         cos_q = cos_q_min + (cos_q_max-cos_q_min)*csi(4)
         w_flat = w_flat *  (cos_q_max-cos_q_min)
         call getphi(csi(1),phi_q,w_phi_q)
         w_flat = w_flat * w_phi_q

         m34 = sqrt(Q2_34)
         call dLips_2(p1+p2,cos_q,phi_q,q,0d0,Q34,m34,phsp_2_34_q,icut)

         call getphi(csi(2),phi3,w_phi3)
         w_flat = w_flat * w_phi3

         cos3_min = -1d0
         cos3_max = 1d0
         cos3 = cos3_min + (cos3_max - cos3_min)*csi(3)
         w_flat = w_flat *(cos3_max-cos3_min)

         call dLips_2(Q34,cos3,phi3,p3,m3,p4,m4,phsp_2_3_4,icut)

         phsp = phsp * phsp_2_34_q * phsp_2_3_4
         phsp = phsp * (2d0*pi)**7
         phsp = phsp*w_flat

         qph(1,0:3) =q(0:3)
      endif

      if(icut.ne.0) then
         phsp=0d0
      endif
      return
      end



      subroutine CHN_1g_soglia_2m_x_ISR(p1,p2,p3,m3,p4,m4,qph,
     . w,norm,phsp,icut,igen)
***************************************************************************
      !! m_34^2 sampled as 1/[m_34^2]
      !! cos_1 of the photon q is sampled as 1/(1-b^2 cos^2)
      implicit none
      real*8 p1,p2,p3,p4,m3,m4,qph
      integer icut,igen
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),qph(40,0:3)
      real*8 q
      dimension q(0:3)
      real*8 w,phsp
      real*8 M,G,Q2_34_min,Q2_34_max,Q2_34,norm_34,w_34

      real*8 beta_e
      real*8 ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi

      real*8 w_smpl,norm,w_flat
      real*8 cmin,cmax
      common/coseni/cmax,cmin

      real*8 cos_q,w_cos_q,norm_cos_q
      real*8 phi_q,w_phi_q

      real*8 csi(3)
      external getrnd
      real*8 dot
      external dot

      real*8 Q34,m34
      dimension Q34(0:3)


      real*8 phsp_2_34_q
      real*8 phi3,w_phi3,cos3,cos3_min,cos3_max,phsp_2_3_4
      real*8 s
      real*8 thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax

      !! initialisation of weights and norms
      w= 1d0
      w_smpl = 1d0
      norm = 1d0
      w_flat = 1d0
      phsp = 1d0

      if(igen.eq.0) then
         q(0:3) = qph(1,0:3)
         Q2_34 = dot(p3+p4,p3+p4)
         Q34= p3+p4
         cos_q = qph(1,3)/qph(1,0)
      endif

      s = dot(p1+p2,p1+p2)
      Q2_34_min = (m3+m4)**2
      Q2_34_max = s - 2d0* sqrt(s)* egmin

      call sampling_energy_mxx(Q2_34_min,Q2_34_max,
     >     Q2_34,w_34,norm_34,igen)
      w_smpl = w_smpl * w_34
      norm   = norm * norm_34

      beta_e = sqrt(1d0-ame**2/p1(0)**2)
      call sampling_cos_collsq(beta_e,1d0,-1d0,1d0,
     .                               cos_q,w_cos_q,norm_cos_q,igen)
      w_smpl = w_smpl * w_cos_q
      norm = norm*norm_cos_q
      w=w*w_smpl

      if(igen.eq.1) then
         call getrnd(csi,3)

         call getphi(csi(1),phi_q,w_phi_q)
         w_flat = w_flat * w_phi_q

         m34 = sqrt(Q2_34)
         call dLips_2(p1+p2,cos_q,phi_q,q,0d0,Q34,m34,phsp_2_34_q,icut)

         call getphi(csi(2),phi3,w_phi3)
         w_flat = w_flat * w_phi3


         cos3_min = -1d0
         cos3_max = 1d0
         cos3 = cos3_min + (cos3_max - cos3_min)*csi(3)
         w_flat = w_flat *(cos3_max-cos3_min)

         call dLips_2(Q34,cos3,phi3,p3,m3,p4,m4,phsp_2_3_4,icut)

         phsp = phsp * phsp_2_34_q * phsp_2_3_4
         phsp = phsp * (2d0*pi)**7
         phsp = phsp*w_flat
         qph(1,0:3)=q(0:3)
      endif
         if(icut.ne.0) then
         phsp=0d0
      endif
      return
      end



      subroutine CHN_1g_soft_coll_ISR(p1,p2,p3,m3,p4,m4,qph,
     . w,norm,phsp,icut,igen)
***************************************************************************
      !! m_34^2 sampled as 1/[(m_34^2 - M^2)^2 + (G*M)^2]
      !! cos_1 of the photon q is sampled as 1/(1-b^2 cos^2)
      implicit none
      real*8 p1,p2,p3,p4,m3,m4,qph
      integer icut,igen
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),qph(40,0:3)
      real*8 q
      dimension q(0:3)
      real*8 w,phsp
      real*8 M,G,Q2_34_min,Q2_34_max,Q2_34,norm_34,w_34

      real*8 beta_e
      real*8 ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi

      real*8 w_smpl,norm,w_flat
      real*8 cmin,cmax
      common/coseni/cmax,cmin

      real*8 cos_q,w_cos_q,norm_cos_q
      real*8 phi_q,w_phi_q

      real*8 csi(3)
      external getrnd
      real*8 dot
      external dot

      real*8 Q34,m34
      dimension Q34(0:3)


      real*8 phsp_2_34_q
      real*8 phi3,w_phi3,cos3,cos3_min,cos3_max,phsp_2_3_4
      real*8 s
      real*8 thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax

      !! initialisation of weights and norms
      w= 1d0
      w_smpl = 1d0
      norm = 1d0
      w_flat = 1d0
      phsp = 1d0

      if(igen.eq.0) then
         q(0:3) = qph(1,0:3)
         Q2_34 = dot(p3+p4,p3+p4)
         Q34= p3+p4
         cos_q = qph(1,3)/qph(1,0)
      endif

      s = dot(p1+p2,p1+p2)
      Q2_34_min = (m3+m4)**2
      Q2_34_max = s - 2d0* sqrt(s)* egmin
      call sampling_prop(s,Q2_34_min,Q2_34_max,Q2_34,w_34,norm_34,igen)
      w_smpl = w_smpl * w_34
      norm   = norm * norm_34

      beta_e = sqrt(1d0-ame**2/p1(0)**2)
      call sampling_cos_collsq(beta_e,1d0,-1d0,1d0,
     .                               cos_q,w_cos_q,norm_cos_q,igen)
      w_smpl = w_smpl * w_cos_q
      norm = norm*norm_cos_q
      w=w*w_smpl

      if(igen.eq.1) then
         call getrnd(csi,3)

         call getphi(csi(1),phi_q,w_phi_q)
         w_flat = w_flat * w_phi_q

         m34 = sqrt(Q2_34)
         call dLips_2(p1+p2,cos_q,phi_q,q,0d0,Q34,m34,phsp_2_34_q,icut)

         call getphi(csi(2),phi3,w_phi3)
         w_flat = w_flat * w_phi3


         cos3_min = -1d0
         cos3_max = 1d0
         cos3 = cos3_min + (cos3_max - cos3_min)*csi(3)
         w_flat = w_flat *(cos3_max-cos3_min)

         call dLips_2(Q34,cos3,phi3,p3,m3,p4,m4,phsp_2_3_4,icut)

         phsp = phsp * phsp_2_34_q * phsp_2_3_4
         phsp = phsp * (2d0*pi)**7
         phsp = phsp*w_flat
         qph(1,0:3)=q(0:3)
      endif
         if(icut.ne.0) then
         phsp=0d0
      endif
      return
      end

      subroutine CHN_1g_coll_FSR(p1,p2,p3,m3,p4,m4,qph,
     . w,norm,phsp,icut,igen)
***************************************************************************
      !! cos of the photon q is sampled as collinear to p3
      implicit none
      real*8 p1,p2,p3,p4,m3,m4,q
      integer icut,igen
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      real*8 qph(40,0:3)
      real*8 w,phsp
      real*8 M,G,Q2_3g_min,Q2_3g_max,Q2_3g,norm_3g,w_3g

      real*8 beta_e
      real*8 ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi

      real*8 w_smpl,norm,w_flat
      real*8 cmin,cmax
      common/coseni/cmax,cmin

      real*8 cos4,w_cos4,norm_cos4
      real*8 phi4,w_phi4

      real*8 csi(4)
      external getrnd
      real*8 dot
      external dot
      real*8 beta_fs,E_fs

      real*8 Q3g,m3g,phsp_2_3_q
      dimension Q3g(0:3)


      real*8 cos4_min,cos4_max


      real*8 phsp_2_34_q
      real*8 phi3,w_phi3,cos3,cos3_min,cos3_max,phsp_2_3g_4

      real*8 cos_q,cos_q_max,cos_q_min
      real*8 w_cos_q,n_cos_q
      real*8 phi_q,w_phi_q

      real*8 q_cm(0:3),p3_cm(0:3)
      real*8 tridot
      external tridot

      !! initialisation of weights and norms
      w= 1d0
      w_smpl = 1d0
      norm = 1d0
      w_flat = 1d0
      phsp = 1d0

      if(igen.eq.0) then
         q(0:3) = qph(1,0:3)
         Q2_3g = dot(p3+q,p3+q)
         m3g=sqrt(Q2_3g)
         Q3g = p3+q
      endif

      Q2_3g_min = (m3)**2
      Q2_3g_max = dot(p1+p2,p1+p2)

      call sampling_energy_mxx(Q2_3g_min,Q2_3g_max,Q2_3g,w_3g,
     >  norm_3g,igen)

      w_smpl = w_smpl * w_3g
      norm   = norm * norm_3g

      if(igen.eq.1) then
         call getrnd(csi,4)

         ! Q2_3g_min = (m3)**2
         ! Q2_3g_max = dot(p1+p2,p1+p2)
         ! Q2_3g = Q2_3g_min + (Q2_3g_max-Q2_3g_min)*csi(4)
         ! w_Flat = w_flat * (Q2_3g_max-Q2_3g_min)



         cos4_min=-1d0
         cos4_max=1d0
         cos4 = cos4_min + (cos4_max-cos4_min)*csi(1)
         w_flat = w_flat *  (cos4_max-cos4_min)

         call getphi(csi(2),phi4,w_phi4)
         w_flat = w_flat * w_phi4

         m3g = sqrt(Q2_3g)
         call dLips_2(p1+p2,cos4,phi4,Q3g,m3g,p4,m4,phsp_2_3g_4,icut)


         cos3_min = -1d0
         cos3_max = 1d0
         cos3 =cos3_min + (cos3_max - cos3_min)*csi(3)
         w_flat = w_flat *(cos3_max-cos3_min)
      endif


   !    E_fs=(m3g**2+m3**2)/(m3g)
   !    print*,'E_fs GENERATED',E_fs
   !    beta_fs=sqrt(1d0-m3**2/E_fs**2)
   !    cos_q_min=-1d0
   !    cos_q_max=1d0
   !    call sampling_cos_coll(beta_fs,1d0,cos_q_min,cos_q_max,
   !   .              cos_q,w_cos_q,n_cos_q,igen)
   !    print*,'cos_q GENERATED',cos_q
   !    if(w_cos_q.lt.0d0) then
   !       print*,'w_cos_q'
   !       print*,'p3',p3
   !       print*,'p4',p4
   !       print*,'q',Q
   !    endif
   !    w_smpl = w_smpl*w_cos_q
   !    norm = norm*n_cos_q


      w=w*w_smpl

      if(igen.eq.1) then
         call getphi(csi(4),phi_q,w_phi_q)
         w_flat = w_flat * w_phi_q
         ! print*,'decadimento Q3g-> q + p3'
         call dLips_2(Q3g,cos3,phi_q,q,0d0,p3,m3,phsp_2_3_q,icut)
         phsp = phsp * phsp_2_3g_4 * phsp_2_3_q
         phsp = phsp * (2d0*pi)**7
         phsp = phsp*w_flat

         qph(1,0:3) =q(0:3)
      endif
      if(icut.ne.0) then
         phsp=0d0
      endif

      return
      end
****************************************************************************

       subroutine dLips_2(Qlab,cth,phi,p1,m1,p2,m2,w,ierr)
************************************************************
!!    Phase space for Q(sqrt(Q2)) -> p1(m1) + p2(m2) decay
!!    Phase-space jacobian computed as
!!    dPhi_2 = 1/2  1/Pi  beta_1/ 8 Pi dcos(th) dphi
!!    input  | Q2: invariant mass squared
!!    output | p1: momentum of particle 1
!!           | m1: mass of particle 1
!!           | p2: momentum of particle 2
!!           | m2: mass of particle 2
!!           |  w: weight of the phase space

!!  Checks: The volume is identical to the analytic one for all masses.
      implicit none
      integer ierr
      real*8 Q2,p1,m1,p2,m2,w,QQ,Qlab(0:3)
      dimension p1(0:3),p2(0:3),QQ(0:3)
      integer i_no_decay
      common/nodecerr/i_no_decay

      real*8 modp2,modp1
      external dot
      real*8 dot,Q
      real*8 epsilon,condition
      real*8 cth,phi,sth
      real*8 lambda
      real*8 tridot
      external tridot
      external lambda
      real*8 pi
      pi=4d0*atan(1d0)
      condition=(m1+m2)**2
      if(ierr.gt.0) THEN
         w=0d0
         return
      endif
      Q2=dot(Qlab,Qlab)
      if(Q2.lt.(m1+m2)**2) then
            w=0d0
            ierr=1
            return
      endif

      Q=sqrt(Q2)
      QQ(0)=Q
      QQ(1)=0d0
      QQ(2)=0d0
      QQ(3)=0d0

      !! In the reference frame where q is AT REST
      !! First decay Q -> p1 + p2
      !! p1^2 = m1^2   on-shell
      !! p2^2 = m2^2   on-shell

      !! dPhi_2 = d^3 p1 d^3 p2 delta(Q-p1-p2) /((2pi)^3 2)^2

      p2(0)=(Q2+m2**2-m1**2)/(2d0*Q)
      modp2=sqrt(abs((p2(0)**2-m2**2)))
      p1(0)=(Q2-m2**2+m1**2)/(2d0*Q)
      modp1=(sqrt(abs(p1(0)**2-m1**2)))

      sth=sqrt(1d0-cth**2)

      !!! I have to apply a rotation
      p1(1)=sth*cos(phi)*modp1
      p1(2)=sth*sin(phi)*modp1
      p1(3)=cth*modp1


      p2(1)=-sth*cos(phi)*modp1
      p2(2)=-sth*sin(phi)*modp1
      p2(3)=-cth*modp1

      !! I compute the weight of the phase space
      w=sqrt(lambda(Q2,m1**2,m2**2))/Q2
      w=w*1d0/(2d0*((2d0*pi)**3 * 2d0)**2)

      !! I boost the momenta back to the lab frame
      call new_boost(Qlab,p1,p1,-1)
      call new_boost(Qlab,p2,p2,-1)


      return
      end


      subroutine phasespace_sampled(ng,p1,p2,p3,p4,qph,
     .     m3,m4,esoft,w,phsp,ie)
      implicit none
      real*8 p1,p2,p3,p4,m3,m4,q
      real*8 esoft,w,phsp,norm
      integer ie,ng,igen,i,j,nisr,nfsr
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      real*8 qph(40,0:3),q1(0:3),q2(0:3)
      real*8 xi(1)
      integer N
c$$$      common/weightschannel/C_1g,C_2g,C_3g
      real*8 C_1g(10),C_2g(49),C_3g(9),C_3g_FSR(16)
      external CHN_i_1g
      external CHN_i_2g
      external CHN_i_3g
      real*8 dot
      external dot
      common/epssoft/eps
      real*8 eps

      integer isr
      common/chooseisrfsr/isr

      p3=0d0
      p4=0d0
      q=0d0
      phsp=0d0
      w=0d0
      qph=0d0

      if(ng.eq.1) then
         N=10
         C_1g=1d0/10.d0
!     print*,C_1g
         call getrnd(xi,1)
         call multichannel(ng,xi(1),N,C_1g,p1,p2,p3,m3,p4,m4,
     .        qph,w,phsp,ie,CHN_i_1g)
         if(ie.ne.0) then
            phsp=0.d0
         endif
      else if(ng.eq.2) then
         N=49
         C_2g=1.d0/49.d0

         call getrnd(xi,1)
         call multichannel(ng,xi(1),N,C_2g,p1,p2,p3,m3,p4,m4,
     .        qph,w,phsp,ie,CHN_i_2g)

         if(ie.ne.0) then
            phsp=0.d0
         endif

         do i=1,2
            if (qph(i,0).le.eps*dsqrt(dot(p1+p2,p1+p2))/2.d0) then
               ie =1
               phsp =0.d0
               w =0.d0
               return
            endif
            do j=1,3
               if((qph(i,0)-abs(qph(i,j)).lt.-1.d-20)) then
                  ie =1
                  phsp =0.d0
                  w =0.d0
                  return
               endif
            enddo
         enddo

      else if(ng.ge.3) then

         call getrnd(xi,1)
         xi(1)=0.2d0

         if(xi(1).le.0.5d0) then
            N=9
            C_3g=1.d0/9.d0

************ho n phot di isr
            isr=1
            call getrnd(xi,1)
            call multichannel(ng,xi(1),N,C_3g,p1,p2,p3,m3,p4,m4,
     .           qph,w,phsp,ie,CHN_i_3g)

!!!!!!!!!!!!!!!!!ADD CHECK CKKW E VEDI CHE TORNA TUTTO
c$$$
c$$$            call ckkw_crit(p1,p2,p3,p4,qph,ng,nisr,nfsr)
c$$$            if(nisr.ne.ng) then
c$$$               phsp=0.d0
c$$$            endif
c$$$*****************multiply by two due to random choice
c$$$            phsp=phsp*2.d0

         else
            N = 16
            C_3g_FSR =1.d0/16.d0

************ho almeno 1phot di fsr
            isr=0
            call getrnd(xi,1)
            call multichannel(ng,xi(1),N,C_3g_FSR,p1,p2,p3,m3,p4,m4,
     .           qph,w,phsp,ie,CHN_i_3g)

c$$$!!!!!!!!!!!!!!!!!ADD CHECK CKKW E VEDI CHE TORNA TUTTO

            call ckkw_crit(p1,p2,p3,p4,qph,ng,nisr,nfsr)
            if(nisr.eq.ng) then
               phsp=0.d0
            endif

*****************multiply by two due to random choice
            phsp=phsp*2.d0

         endif
         do i=1,ng
            if (qph(i,0).le.eps*dsqrt(dot(p1+p2,p1+p2))/2.d0) then
               ie =1
               phsp =0.d0
               w =0.d0
               return

            endif
         enddo


      endif

      return
      end

************************************************************************************
      !!===========================================!!
      !!   Sampling routines for e+e- -> gg e+e-   !!
      !!===========================================!!

      !! This routine samples cos(theta) as 1/(1-beta^2 cos^2(theta))
      !! input  | en: energy of the leg
      !!        | cmin cmax: max and min cos(theta)
      !!        | igen: 0 computes the weight, 1 generates variables
      !! output | cos: sampled cos(theta)
      !!        | wgt: the weight (1-b^2 cos^2)
      !!        | norm: int_cmin^cmax 1/(1-b^2 cos^2) dcos
      subroutine sampling_cos_collsq(beta_e,beta_q,cmin,cmax,
     .                               cos,wgt,norm,igen)
      implicit none
      real*8 cos,wgt,cmin,cmax
      real*8 beta_e,beta_q
      real*8 b,fac,norm,en,arg
      real*8 xi(1)
      real*8 logmax,logmin
      real*8 ame,ammu,convfac,alpha,pi
      real*16 cos_q,w_q,b_q
      integer igen
      common/parameters/ame,ammu,convfac,alpha,pi

      b = abs(beta_e*beta_q)
      logmax = log((1d0 + b*cmax)/(1d0-b*cmax))
      logmin = log((1d0 - b*cmin)/(1d0+b*cmin))
      norm = 1d0/(2d0*b)*(logmax+logmin)

      if(igen.eq.1) then
         call getrnd(xi,1)
         fac = (1d0+b*cmin)/(1d0-b*cmin)
         arg = 2d0*b*norm*xi(1)+log(fac)
         cos = 1d0/b*((exp(arg)-1d0)/(exp(arg)+1d0))
      else if(igen.eq.0) then
         cos=cos
      endif

      wgt = (1d0-b*b*cos*cos)
      if(abs(1q0-cos*1q0).le.1q-10.or.abs(1q0+cos*1q0).le.1q-10) then
         ! print*,'cosq',cos*1q0
         ! print*,'wgt  ',wgt
         cos_q=cos*1q0
         b_q=b*1q0
         w_q = (1q0-b_q**2*cos_q**2)
         b = abs(beta_e*beta_q)
         wgt=w_q*1D0
         ! print*,'wgt q',wgt
      endif
      return
      end

      !! This routine samples the variable E (energy) as 1/E
      !! input  | emin,emax : min and max energy
      !!        | igen : =0 computes the weight from inputs, =1 generates variables
      !! output | en : sampled energy
      !!        | wgt : weight 1/(1/E)
      !!        | norm : int_emin^emax 1/E
      subroutine sampling_energy(emin,emax,en,wgt,norm,igen)
      implicit none
      real*8 emin,emax,en,xi(1),wgt,norm
      external getrnd
      integer igen
      common/momentainitial/pin1,pin2
      real*8  pin1(0:3),pin2(0:3)
      real*8 BW_min,BW_max,s
      external dot
      real*8 dot

      s =dot(pin1+pin2,pin1+pin2)

      if(igen.eq.1) then
         call getrnd(xi,1)
         en = emin*(emax/emin)**xi(1)
      else if(igen.eq.0) then
         en = en
      endif
      wgt = en/dsqrt(s)
      norm=log(emax/emin)*dsqrt(s)

      return
      end
*************************

      !! This routine samples the variable s (energy) as 1/s dimension GeV^2
      !! input  | emin,emax : min and max energy
      !!        | igen : =0 computes the weight from inputs, =1 generates variables
      !! output | en : sampled energy sqrd
      !!        | wgt : weight 1/(1/s)
      !!        | norm : int_emin^emax 1/s
      subroutine sampling_energy_mxx(emin,emax,en,wgt,norm,igen)
      implicit none
      real*8 emin,emax,en,xi(1),wgt,norm
      external getrnd
      integer igen
      common/momentainitial/pin1,pin2
      real*8  pin1(0:3),pin2(0:3)
      real*8 BW_min,BW_max,s
      external dot
      real*8 dot

      s = dot(pin1+pin2,pin1+pin2)

      if(igen.eq.1) then
         call getrnd(xi,1)
         en = emin*(emax/emin)**xi(1)
      else if(igen.eq.0) then
         en = en
      endif
      wgt = en/s
      norm=log(emax/emin)*s

      return
      end
*************************


!! This routine samples the variable s as 1/(s-m^2) dimension GeV^2
!! input  | emin,emax : min and max energy
!!        | igen : =0 computes the weight from inputs, =1 generates variables
!! output | en : sampled energy sqrd
!!        | wgt : weight 1/(1/s-m^2)
!!        | norm : int_emin^emax 1/(s-m^2)
      subroutine sampling_energy_mxx_sc(m,Q2min,Q2max,Q2,w,norm,igen)
      implicit none
      real*8 s,Q2min,Q2max,Q2,w,norm
      integer igen
      real*8 xi(1),arg_tan,m,emin,emax
      real*8 log_min,log_max
      common/momentainitial/pin1,pin2
      real*8  pin1(0:3),pin2(0:3)
      real*8 BW_min,BW_max,s1
      external dot
      real*8 dot


      s = dot(pin1+pin2,pin1+pin2)

      log_min = log(Q2min-m**2)
      log_max = log(Q2max-m**2)
      norm = (-log_min + log_max)*s

      if(igen.eq.1) then
         call getrnd(xi,1)
         Q2 = m**2 + exp(1d0*norm*xi(1)/s)*(Q2min-m**2)

      endif
      w = (Q2-m**2)/s


      return
      end

**************************
      !!=========================!!
      !!  Breit-Wigner sampling  !!
      !!=========================!!
      !!  Q2 is generated according BW^2 = 1/[(Q2-M^2)^2 + (G M)^2]
      !!  input  |          M,G: mass and width of Breit-Wigner
      !!         | Q2min, Q2max: generation limits on Q2
      !!         |         igen: 0 normalisation, 1 generation
      !!  output |           Q2: virtuality
      !!         |            w: weight
      !!         |         norm: normalisation
      subroutine sampling_BW(M,G,Q2min,Q2max,Q2,w,norm,igen)
      implicit none
      real*8 Q2min,Q2max,M,G,Q2,w,norm
      integer igen
      real*8 xi(1),arg_tan
      common/momentainitial/pin1,pin2
      real*8  pin1(0:3),pin2(0:3)
      real*8 BW_min,BW_max,s
      external dot
      real*8 dot

      s =dot(pin1+pin2,pin1+pin2)

      BW_min = atan((Q2min-M**2)/(G*M))
      BW_max = atan((Q2max-M**2)/(G*M))
!      norm = 1/(G*M)*(BW_max-BW_min)

      norm = (BW_max-BW_min)*s
      if(igen.eq.1) then
         call getrnd(xi,1)
!         arg_tan = G*M*norm*xi(1) + BW_min
         arg_tan = norm*xi(1)/s + BW_min
         Q2 = M**2 + G*M*tan(arg_tan)
      endif

!      w = (Q2-M**2)**2 + (G*M)**2
      w = ((Q2-M**2)**2 + (G*M)**2)/(G*M)/s

      return
      end

*******
      subroutine sampling_BW_resonances(M,G,Q2min,Q2max,Q2,w,norm,igen)
      implicit none
      real*8 Q2min,Q2max,M,G,Q2,w,norm,Q2mini,Q2maxi
      real*8 C_max, C_min,W_i,W_j,N_i,N_j
      integer i_sig,i_succ,ie
      character*2 fs
      common/finalstate/fs
      real*8 xi(1),arg_tan,wgt_sum
      real*8 ame,ammu,convfac,alpha,pi
      integer igen,i,ires,j
      common/parameters/ame,ammu,convfac,alpha,pi
      real*8 BW_min,BW_max,sum,normtot
      integer nres
      parameter (nres=4)
      double precision MA(nres),TG(nres),EG(nres),C(nres)
      common/resonancesparams1/ma,tg,eg
      common/momentainitial/pin1,pin2
      real*8 pin1(0:3),pin2(0:3),sqs
      common/indexres/indres,nreseff
      integer indres,nreseff
      external dot
      real*8 dot
****  UPDATED (from PDG 2024) FOR VP RESONANCES!!
      DATA MA/                  ! GeV
     >     0.78266D0,           ! 1.  omega
     >     1.019460D0,          ! 2.  phi
     >     3.096906D0,          ! 3.  J/Psi
     >     3.686097D0/           ! 4.  Psi(2S)
      DATA TG/                  ! GeV
     >     8.68D-3,             ! 1.  omega
     >     4.249D-3,            ! 2.  phi
     >     92.6D-6,             ! 3.  J/Psi
     >     293.D-6/              ! 4.  Psi(2S)
      DATA EG/                  ! GeV
     >     0.643D-6,            ! 1.  omega
     >     1.256D-6,            ! 2.  phi
     >     5.53D-6,             ! 3.  J/Psi
     >     2.326D-6/             ! 4.  Psi(2S)


      ie=0
      norm=1.d0
      w=1.d0


      if(indres.eq.0) then
         sqs=dsqrt(dot(pin1+pin2,pin1+pin2))
         nreseff=0
         do i=1,nres
            if(sqs.le.MA(i)) then
               goto 118
            else
               nreseff=nreseff+1
            endif
         enddo
 118     indres=1
      endif

      if(fs.eq.'pr') then
         call sampling_BW(M,G,Q2min,Q2max,Q2,w,norm,igen)
         return

      elseif(fs.eq.'mr') then
         call getrnd(xi,1)


         do i=1,nreseff
            C(i)=1.d0/float(nreseff)
         enddo

         normtot=0.d0
         C_min=0d0
         wgt_sum=0d0
         ie=0

         do i=1,nreseff
            w_i=0d0
            N_i=0d0
!! x_i in [C_min,C_min+C(i)]
            C_max = C_min+C(i)
!! I call the generating channel
            if(C_min.lt.xi(1).and.xi(1).le.C_max) then
               i_sig=i
               call sampling_BW(MA(i),TG(i),Q2min,Q2max,Q2,w_i,N_i,igen)

!     print*,'CHN GEN=',i,'ie',ie,'PHSP',PHSP_i,'W_i',W_i
               if(isnan(N_i)) then

                  print*,'!!!!!!!!!!!!!! ABORT !!!!!!!!!!!'
                  stop
               ENDIF

               wgt_sum = wgt_sum + C(i)/W_i
               norm=N_i
               i_succ=I
            endif
            C_min = C_max

         enddo
!! I compute the normalisation calling all the others i!=j channels
         do j=1,nreseff
            if(j.ne.i_succ) then
               w_j=0d0
               N_j=0d0
               call sampling_BW(MA(j),TG(j),Q2min,Q2max,Q2,w_j,N_j,0)
               wgt_sum = wgt_sum + C(j)/W_j
               if(isnan(N_j)) then
                  print*,'!!!!!!!!!!!!!! ABORT !!!!!!!!!!!'
                  stop
               ENDIF
            endif

         enddo

!     ! the weight is given by the normalisation divided by the sum_i C_i/W_i
         W=1.d0/wgt_sum

         return
      endif
      end

      !!===============================!!
      !!  Massive Propagator sampling  !!
      !!===============================!!
      !!  Q2 is generated as 1/(s - Q2)
      !!  input  |            s: virtuality flowing into the propagator
      !!         | Q2min, Q2max: generation limits on Q2
      !!         |         igen: 0 normalisation, 1 generation
      !!  output |           Q2: virtuality
      !!         |            w: weight
      !!         |         norm: normalisation
      subroutine sampling_prop(s,Q2min,Q2max,Q2,w,norm,igen)
      implicit none
      real*8 s,Q2min,Q2max,Q2,w,norm
      integer igen
      real*8 xi(1),arg_tan
      real*8 log_min,log_max
      common/momentainitial/pin1,pin2
      real*8  pin1(0:3),pin2(0:3)
      real*8 BW_min,BW_max,s1
      external dot
      real*8 dot

      s1 =dot(pin1+pin2,pin1+pin2)

      log_min = log(s-Q2min)
      log_max = log(s-Q2max)
      norm = (log_min - log_max)*s1

      if(igen.eq.1) then
         call getrnd(xi,1)
         Q2 = s - exp(-1d0*norm*xi(1)/s1)*(s-Q2min)
      endif

      w = (s - Q2)/s1

      return
      end

      !! This routine samples cos(theta) as 1/(1-beta cos(theta))
      !! input  | en: energy of the leg
      !!        | cmin cmax: max and min cos(theta)
      !!        | igen: 0 computes the weight, 1 generates variables
      !! output | cos: sampled cos(theta)
      !!        | wgt: the weight (1-b cos)
      !!        | norm: int_cmin^cmax 1/(1-b cos) dcos
      subroutine sampling_cos_coll(beta_e,beta_q,cmin,cmax,
     .                               cos,wgt,norm,igen)
      implicit none
      real*8 cos,wgt,cmin,cmax
      real*8 beta_e,beta_q
      real*8 b,fac,norm,en,arg
      real*8 xi(1)
      real*8 logmax,logmin
      real*8 ame,ammu,convfac,alpha,pi
      integer igen
      common/parameters/ame,ammu,convfac,alpha,pi

      real*16 b_q,norm_q,wgt_q

      b = beta_e*beta_q
      logmax = log((1d0 - b*cmax))
      logmin = log((1d0 - b*cmin))
      norm = 1d0/(b)*(logmin-logmax)
      if(igen.eq.1) then
         call getrnd(xi,1)
         arg = -b*norm*xi(1)
         cos = 1d0/b*(1d0-exp(arg)*(1d0-b*cmin))
      else if(igen.eq.0) then
         cos=cos
      endif

      wgt = (1d0-b*cos)

      if((1d0-b).lt.1d-10) then
      if(igen.eq.1) then
         b_q=b*1q0
         logmax = log((1q0 - b_q*cmax*1q0)*1d0)
         logmin= log((1q0 - b_q*cmin*1q0)*1d0)
         norm_q = 1q0/(b_q)*(logmin-logmax)
         arg = -b_q*norm_q*xi(1)*1d0
         cos = 1d0/b*(1d0-exp(arg)*(1d0-b*cmin))
      else if(igen.eq.0) then
         cos=cos
      endif

      wgt_q = (1q0-b_q*cos*1q0)
      wgt=wgt_q*1d0
      endif
      return
      end
************************************************************************************

      !!=============================!!
      !!    MULTICHANNEL SAMPLING    !!
      !!=============================!!
      !! This routine implements the multichannel sampling,
      !! for an arbitrary N number of channels

      !! input  | xi    : random number [0,1]
      !!        | C(N)  : vector of N coefficients, relative C of each channel
      !!        | p1,p2 : initial state momenta
      !! output | p3(m3): final state electron
      !!        | p4(m4): final state positron
      !!        | q1,q2 : photons
      !!        | w     : sampled weight
      !!        | phsp  : phasepsace volume with flat jacobians
      subroutine multichannel(ng,xi,N,C,p1,p2,p3,m3,p4,m4,
     .                        qph,w,phsp,ie,CHN_i)
      implicit none
      real*8 xi,C(N)
      integer i,j,N,i_sig,ng
      common/CHNi/i_sig
      external dot
      real*8 dot
      real*8 C_min,C_max
      real*8 p1,p2,p3,p4,q1,q2
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q1(0:3),q2(0:3)
      real*8 m3,m4,w,phsp,sum,m2
      integer ie,i_succ
      real*8 qph(40,0:3)
      external CHN_i
      logical cond
      integer nres
      parameter (nres=4)
      double precision MA(nres),TG(nres),EG(nres)
      common/resonancesparams1/ma,tg,eg
      !! Variables i=generation channel
      !!           j=normalisation channels
      real*8 W_i,n_i,phsp_i
      real*8 W_j,n_j,phsp_j
      real*8 wgt_sum,norm
      common/channel_multi/i_ch
      integer i_ch


      C_min=0d0
      wgt_sum=0d0
      ie=0
      phsp=0d0
      do i=1,N
         w_i=0d0
         N_i=0d0
         !! x_i in [C_min,C_min+C(i)]
         C_max = C_min+C(i)
         !! I call the generating channel
         if(C_min.lt.xi.and.xi.le.C_max) then
            i_sig=i
            if(ng.le.2) then
               call CHN_i(i,p1,p2,p3,m3,p4,m4,
     .              qph,w_i,N_i,phsp_i,ie,1)

            else
               call CHN_i(ng,i,p1,p2,p3,m3,p4,m4,
     .              qph,w_i,N_i,phsp_i,ie,1)
            endif
! print*,'CHN GEN=',i,'ie',ie,'PHSP',PHSP_i,'W_i',W_i

            if(isnan(N_i)) then
               print*,'!!!!!!!!!!!!!! ABORT !!!!!!!!!!!'
               stop
            ENDIF
            phsp=phsp_i
            wgt_sum = wgt_sum + C(i)/W_i
            norm=N_i
            i_succ=I
         endif
         C_min = C_max

      enddo
            !! I compute the normalisation calling all the others i!=j channels

      if(ie.eq.0) then
         do j=1,N
            if(j.ne.i_succ) then
               phsp_j=0d0
               w_j=0d0
               N_j=0d0
               if(ng.le.2) then

                  call CHN_i(j,p1,p2,p3,m3,p4,m4,
     .                 qph,w_j,N_j,phsp_j,ie,0)
               else
                  call CHN_i(ng,j,p1,p2,p3,m3,p4,m4,
     .                 qph,w_j,N_j,phsp_j,ie,0)
               endif
!     print*,'CHN NORM=',j,'ie',ie,'PHSP',PHSP_j,'W_i',W_j
               wgt_sum = wgt_sum + C(j)/W_j

               if(isnan(N_j)) then
                  print*,'!!!!!!!!!!!!!! ABORT !!!!!!!!!!!'
                  stop
               ENDIF
            endif
         enddo
         endif

      !! the weight is given by the normalisation divided by the sum_i C_i/W_i
c$$$      if(ie.eq.0) then
c$$$         sum=0.d0
c$$$         print*,'wgt_sum',wgt_sum
c$$$         m2=dot(p3+p4,p3+p4)
c$$$         do i=1,4
c$$$          sum = sum+1.d0/16.d0*1.d0/((m2-MA(i)**2)**2+TG(i)**2*MA(i)**2)
c$$$         enddo
c$$$        sum=sum+1.d0/4.d0/(m2)+1.d0/4.d0/(dot(p1+p2,p1+p2)-m2)+1.d0/4.d0
c$$$        print*,'sum',sum
c$$$      endif
c$$$      if(wgt_sum.le.0.3d0) then
c$$$         print*,0.25d0/dot(p3+p4,p3+p4)
c$$$         print*,wgt_sum
c$$$         print*,dsqrt(dot(p3+p4,p3+p4))
c$$$      endif
      i_ch=i_succ

      W=norm/wgt_sum

c$$$      if(W.le.0.d0) then
c$$$
c$$$         print*,norm,wgt_sum,w,i_succ
c$$$      endif
c$$$

      if (ie.eq.1) then
            phsp=0d0
      endif
      return
      end


      !!=================================!!
      !!   i-th CHN for e+e- -> X+X- g   !!
      !!=================================!!
      subroutine CHN_i_1g(i,p1,p2,p3,m3,p4,m4,
     . qph,w,N,phsp,ie,igen)

      implicit none
      real*8 p1,p2,p3,p4,q1,q2
      real*8 m3,m4
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q1(0:3),q2(0:3)
      real*8 w,N,phsp
      integer ie,igen,i
      real*8 qph(40,0:3)
      real*8 tridot
      external tridot
      real*8 norm

      ie=0

      if(i.eq.1) then
         call CHN_1g_ISR(p1,p2,p3,m3,p4,m4,qph,
     .        w,N,phsp,ie,'BW','fl',igen)
      else if(i.eq.2) then
         call CHN_1g_ISR(p1,p2,p3,m3,p4,m4,qph,
     .        w,N,phsp,ie,'BW','co',igen)

      else if(i.eq.3) then
         call CHN_1g_ISR(p1,p2,p3,m3,p4,m4,qph,
     .        w,N,phsp,ie,'SC','fl',igen)

      else if (i.eq.4) then
         call CHN_1g_ISR(p1,p2,p3,m3,p4,m4,qph,
     .        w,N,phsp,ie,'SC','co',igen)

      else if(i.eq.5) then
         call CHN_1g_ISR(p1,p2,p3,m3,p4,m4,qph,
     .        w,N,phsp,ie,'CO','fl',igen)
      else if(i.eq.6) then
         call CHN_1g_ISR(p1,p2,p3,m3,p4,m4,qph,
     .        w,N,phsp,ie,'CO','co',igen)
      else if(i.eq.7) then
         call CHN_1g_ISR(p1,p2,p3,m3,p4,m4,qph,
     .        w,N,phsp,ie,'FL','co',igen)
      else if(i.eq.8) then
         call CHN_1g_ISR(p1,p2,p3,m3,p4,m4,qph,
     .           w,N,phsp,ie,'FL','fl',igen)
      else if(i.eq.9) then

         call CHN_1g_coll_FSR(p1,p2,p3,m3,p4,m4,qph,
     . w,N,phsp,ie,igen)

      else if(i.eq.10) then
         call CHN_1g_coll_FSR(p1,p2,p4,m4,p3,m3,qph,
     . w,N,phsp,ie,igen)

      endif
      return
      end

      !!==================================!!
      !!   i-th CHN for e+e- -> X+X- gg   !!
      !!==================================!!
      subroutine CHN_i_2g(i,p1,p2,p3,m3,p4,m4,
     . qph,w,N,phsp,ie,igen)

      implicit none
      real*8 p1,p2,p3,p4,q1,q2
      real*8 m3,m4
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q1(0:3),q2(0:3)
      real*8 w,N,phsp
      integer ie,igen,i
      real*8 qph(40,0:3)
      real*8 tridot
      external tridot
      real*8 dot
      external dot
      real*8 norm
      common/epssoft/eps
      real*8 eps
      ie=0

      if(igen.eq.0) then
         q1(0:3)=qph(1,0:3)
         q2(0:3)=qph(2,0:3)
      endif




      if(i.eq.1) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,'BW','fl',igen)
      elseif(i.eq.2) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,'BW','co',igen)

      elseif(i.eq.3) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,'CO','fl',igen)

      elseif(i.eq.4) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .           q1,q2,w,N,phsp,ie,'CO','co',igen)

      elseif(i.eq.5) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,'SC','co',igen)
      elseif(i.eq.6) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,'SC','fl',igen)

      elseif(i.eq.7) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,'FL','co',igen)
      elseif(i.eq.8) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,'FL','fl',igen)

c$$$***************RADDOPPIO
      elseif(i.eq.9) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .        q2,q1,w,N,phsp,ie,'BW','fl',igen)
      elseif(i.eq.10) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .        q2,q1,w,N,phsp,ie,'BW','co',igen)

      elseif(i.eq.11) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .        q2,q1,w,N,phsp,ie,'CO','fl',igen)

      elseif(i.eq.12) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .           q2,q1,w,N,phsp,ie,'CO','co',igen)

      elseif(i.eq.13) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .        q2,q1,w,N,phsp,ie,'SC','co',igen)
      elseif(i.eq.14) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .        q2,q1,w,N,phsp,ie,'SC','fl',igen)

      elseif(i.eq.15) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .        q2,q1,w,N,phsp,ie,'FL','co',igen)
      elseif(i.eq.16) then
         call CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     .        q2,q1,w,N,phsp,ie,'FL','fl',igen)

c$$$c$$$******************************
      elseif(i.eq.17) then
         call CHN_2g_flat_ISR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,igen)

      elseif(i.eq.18) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,'BW','fl',igen)
      elseif(i.eq.19) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q1,q2,w,N,phsp,ie,'BW','fl',igen)

      elseif(i.eq.20) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,'BW','co',igen)
      elseif(i.eq.21) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q1,q2,w,N,phsp,ie,'BW','co',igen)

      elseif(i.eq.22) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,'CO','fl',igen)
      elseif(i.eq.23) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q1,q2,w,N,phsp,ie,'CO','fl',igen)

      elseif(i.eq.24) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,'CO','co',igen)
      elseif(i.eq.25) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q1,q2,w,N,phsp,ie,'CO','co',igen)

      elseif(i.eq.26) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,'SC','fl',igen)
      elseif(i.eq.27) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q1,q2,w,N,phsp,ie,'SC','fl',igen)

      elseif(i.eq.28) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,'SC','co',igen)
      elseif(i.eq.29) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q1,q2,w,N,phsp,ie,'SC','co',igen)
      elseif(i.eq.30) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,'FL','co',igen)
      elseif(i.eq.31) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q1,q2,w,N,phsp,ie,'FL','co',igen)

      elseif(i.eq.32) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q1,q2,w,N,phsp,ie,'FL','fl',igen)
      elseif(i.eq.33) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q1,q2,w,N,phsp,ie,'FL','fl',igen)


      elseif(i.eq.40) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q2,q1,w,N,phsp,ie,'BW','fl',igen)
      elseif(i.eq.41) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q2,q1,w,N,phsp,ie,'BW','fl',igen)
      elseif(i.eq.42) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q2,q1,w,N,phsp,ie,'BW','co',igen)
      elseif(i.eq.43) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q2,q1,w,N,phsp,ie,'BW','co',igen)

      elseif(i.eq.44) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q2,q1,w,N,phsp,ie,'CO','fl',igen)
      elseif(i.eq.45) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q2,q1,w,N,phsp,ie,'CO','fl',igen)

      elseif(i.eq.46) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q2,q1,w,N,phsp,ie,'CO','co',igen)
      elseif(i.eq.47) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q2,q1,w,N,phsp,ie,'CO','co',igen)

      elseif(i.eq.48) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q2,q1,w,N,phsp,ie,'SC','fl',igen)
      elseif(i.eq.49) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q2,q1,w,N,phsp,ie,'SC','fl',igen)
      elseif(i.eq.39) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q2,q1,w,N,phsp,ie,'SC','co',igen)
      elseif(i.eq.38) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q2,q1,w,N,phsp,ie,'SC','co',igen)
      elseif(i.eq.37) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q2,q1,w,N,phsp,ie,'FL','co',igen)
      elseif(i.eq.36) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q2,q1,w,N,phsp,ie,'FL','co',igen)
      elseif(i.eq.35) then
         call CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,
     .        q2,q1,w,N,phsp,ie,'FL','fl',igen)
      elseif(i.eq.34) then
         call CHN_2g_ISR_FSR(p1,p2,p4,m4,p3,m3,
     .        q2,q1,w,N,phsp,ie,'FL','fl',igen)

      endif
c$$$
c$$$
c$$$      if(i.eq.1) then
c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$     .        q1,q2,w,N,phsp,ie,'FL',1,igen)
c$$$      elseif(i.eq.2) then
c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$     .        q1,q2,w,N,phsp,ie,'FL',2,igen)
c$$$      elseif(i.eq.3) then
c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$     .        q1,q2,w,N,phsp,ie,'FL',3,igen)
c$$$      elseif(i.eq.4) then
c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$     .        q1,q2,w,N,phsp,ie,'CO',1,igen)
c$$$      elseif(i.eq.5) then
c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$     .        q1,q2,w,N,phsp,ie,'CO',2,igen)
c$$$      elseif(i.eq.6) then
c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$     .        q1,q2,w,N,phsp,ie,'CO',3,igen)
c$$$      elseif(i.eq.7) then
c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$     .        q1,q2,w,N,phsp,ie,'SC',1,igen)
c$$$      elseif(i.eq.8) then
c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$     .        q1,q2,w,N,phsp,ie,'SC',2,igen)
c$$$      elseif(i.eq.9) then
c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$     .        q1,q2,w,N,phsp,ie,'SC',3,igen)
c$$$      elseif(i.eq.10) then
c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$     .        q1,q2,w,N,phsp,ie,'BW',1,igen)
c$$$      elseif(i.eq.11) then
c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$     .        q1,q2,w,N,phsp,ie,'BW',2,igen)
c$$$      elseif(i.eq.12) then
c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$     .        q1,q2,w,N,phsp,ie,'BW',3,igen)
c$$$c$$$c$$$      elseif(i.eq.13) then
c$$$c$$$c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$c$$$c$$$     .        q2,q1,w,N,phsp,ie,'FL',1,igen)
c$$$c$$$c$$$      elseif(i.eq.14) then
c$$$c$$$c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$c$$$c$$$     .        q2,q1,w,N,phsp,ie,'FL',2,igen)
c$$$c$$$c$$$      elseif(i.eq.15) then
c$$$c$$$c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$c$$$c$$$     .        q2,q1,w,N,phsp,ie,'FL',3,igen)
c$$$c$$$c$$$      elseif(i.eq.16) then
c$$$c$$$c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$c$$$c$$$     .        q2,q1,w,N,phsp,ie,'CO',1,igen)
c$$$c$$$c$$$      elseif(i.eq.17) then
c$$$c$$$c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$c$$$c$$$     .        q2,q1,w,N,phsp,ie,'CO',2,igen)
c$$$c$$$c$$$      elseif(i.eq.18) then
c$$$c$$$c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$c$$$c$$$     .        q2,q1,w,N,phsp,ie,'CO',3,igen)
c$$$c$$$c$$$      elseif(i.eq.19) then
c$$$c$$$c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$c$$$c$$$     .        q2,q1,w,N,phsp,ie,'SC',1,igen)
c$$$c$$$c$$$      elseif(i.eq.20) then
c$$$c$$$c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$c$$$c$$$     .        q2,q1,w,N,phsp,ie,'SC',2,igen)
c$$$c$$$c$$$      elseif(i.eq.21) then
c$$$c$$$c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$c$$$c$$$     .        q2,q1,w,N,phsp,ie,'SC',3,igen)
c$$$c$$$c$$$      elseif(i.eq.22) then
c$$$c$$$c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$c$$$c$$$     .        q2,q1,w,N,phsp,ie,'BW',1,igen)
c$$$c$$$c$$$      elseif(i.eq.23) then
c$$$c$$$c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$c$$$c$$$     .        q2,q1,w,N,phsp,ie,'BW',2,igen)
c$$$c$$$c$$$      elseif(i.eq.24) then
c$$$c$$$c$$$         call CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
c$$$c$$$c$$$     .        q2,q1,w,N,phsp,ie,'BW',3,igen)
c$$$c$$$c$$$
c$$$      endif

      qph(1,0:3)=q1(0:3)
      qph(2,0:3)=q2(0:3)

      return
      end

********N_PHOT
      !!==================================!!
      !!   i-th CHN for e+e- -> X+X- (n)g !!
      !!==================================!!
      subroutine CHN_i_3g(ng,i,p1,p2,p3,m3,p4,m4,
     . qph,w,N,phsp,ie,igen)

      implicit none
      real*8 p1,p2,p3,p4,q1,q2
      real*8 m3,m4
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q1(0:3),q2(0:3)
      real*8 w,N,phsp
      integer ie,igen,i,ng,isr
      common/chooseisrfsr/isr
      real*8 qph(40,0:3)
      real*8 tridot
      external tridot
      real*8 dot
      external dot
      real*8 norm
      common/epssoft/eps
      real*8 eps
      ie=0

      if(isr.eq.1) then

         if(i.eq.1) then
            call CHN_ng_soft_ISR(ng,p1,p2,p3,m3,p4,m4,
     .           qph,w,N,phsp,ie,'BW','co',igen)
         elseif(i.eq.2) then
            call CHN_ng_soft_ISR(ng,p1,p2,p3,m3,p4,m4,
     .           qph,w,N,phsp,ie,'BW','fl',igen)

         elseif(i.eq.3) then
            call CHN_ng_soft_ISR(ng,p1,p2,p3,m3,p4,m4,
     .           qph,w,N,phsp,ie,'CO','co',igen)
         elseif(i.eq.4) then
            call CHN_ng_soft_ISR(ng,p1,p2,p3,m3,p4,m4,
     .           qph,w,N,phsp,ie,'CO','fl',igen)

         elseif(i.eq.5) then
            call CHN_ng_soft_ISR(ng,p1,p2,p3,m3,p4,m4,
     .           qph,w,N,phsp,ie,'SC','co',igen)
         elseif(i.eq.6) then
            call CHN_ng_soft_ISR(ng,p1,p2,p3,m3,p4,m4,
     .           qph,w,N,phsp,ie,'SC','fl',igen)

         elseif(i.eq.7) then
            call CHN_ng_soft_ISR(ng,p1,p2,p3,m3,p4,m4,
     .           qph,w,N,phsp,ie,'FL','co',igen)

         elseif(i.eq.8) then
            call CHN_ng_soft_ISR(ng,p1,p2,p3,m3,p4,m4,
     .           qph,w,N,phsp,ie,'FL','fl',igen)


         elseif(i.eq.9) then
            call CHN_ng_flat(ng,p1,p2,p3,m3,p4,m4,
     .           qph,w,N,phsp,ie,igen)
         endif

      elseif(isr.eq.0) then

         if(i.eq.1) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p3,m3,p4,m4,
     .              qph,w,N,phsp,ie,'FL','fl',igen)

         elseif(i.eq.2) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p4,m4,p3,m3,
     .              qph,w,N,phsp,ie,'FL','fl',igen)

         elseif(i.eq.3) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p3,m3,p4,m4,
     .              qph,w,N,phsp,ie,'FL','co',igen)

         elseif(i.eq.4) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p4,m4,p3,m3,
     .              qph,w,N,phsp,ie,'FL','co',igen)

         elseif(i.eq.5) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p3,m3,p4,m4,
     .              qph,w,N,phsp,ie,'BW','fl',igen)

         elseif(i.eq.6) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p4,m4,p3,m3,
     .              qph,w,N,phsp,ie,'BW','fl',igen)

         elseif(i.eq.7) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p3,m3,p4,m4,
     .              qph,w,N,phsp,ie,'BW','co',igen)

         elseif(i.eq.8) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p4,m4,p3,m3,
     .              qph,w,N,phsp,ie,'BW','co',igen)
         elseif(i.eq.9) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p3,m3,p4,m4,
     .              qph,w,N,phsp,ie,'CO','fl',igen)

         elseif(i.eq.10) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p4,m4,p3,m3,
     .              qph,w,N,phsp,ie,'CO','fl',igen)

         elseif(i.eq.11) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p3,m3,p4,m4,
     .              qph,w,N,phsp,ie,'CO','co',igen)

         elseif(i.eq.12) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p4,m4,p3,m3,
     .              qph,w,N,phsp,ie,'CO','co',igen)
         elseif(i.eq.13) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p3,m3,p4,m4,
     .              qph,w,N,phsp,ie,'SC','fl',igen)

         elseif(i.eq.14) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p4,m4,p3,m3,
     .              qph,w,N,phsp,ie,'SC','fl',igen)

         elseif(i.eq.15) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p3,m3,p4,m4,
     .              qph,w,N,phsp,ie,'SC','co',igen)

         elseif(i.eq.16) then
            call CHN_ng_ISR_FSR(ng,p1,p2,p4,m4,p3,m3,
     .              qph,w,N,phsp,ie,'SC','co',igen)

         endif
      endif

      return
      end

**************SOFT COLLINEAR IS CHANNEL
        !! 2->4 phase space as diagram (A) of the notes
      !!      \/
      !!      *
      !! ----*----*~~~~~~~~
      !!          |
      !!          |
      !! ---------*~~~~~~~~
      subroutine CHN_2g_soft_ISR(p1,p2,p3,m3,p4,m4,
     . q1,q2,w,norm,phsp,icut,sampling,angdip,igen)
*******************************************************************
      implicit none
      dimension p1(0:3),p2(0:3)
      real*8 qtmp(0:3)
      dimension p3(0:3),p4(0:3),Q_34(0:3),q1(0:3),q2(0:3)
      real*8 Q_34
      integer icut
      real*8 p1,p2,p3,p4,q1,q2,w,norm
      real*8 Qmin,Qmax,emax
      real*8 cmin_34,cmax_34,cmin_g,cmax_g
      real*8 cmax,cmin
      real*8 qph(40,0:3),xi(1)
      common/coseni/cmax,cmin
      real*8 ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
      real*8 s
      real*8 tridot,dot
      external tridot,dot
      real*8 mod234,phsp,phsp_2
      real*8 m3,m4
      logical cond
      integer igen

      !! Weights: flat and sampled
      real*8 w_flat, w_smpl

      !! Variables of Q34
      real*8 Q2_34,E_34,cos_34,sin_34,phi_34,modQ34
      real*8 beta_e,beta_34

      !! Weights of Q34
      real*8 w_Q2_34,w_E_34,w_cos_34,w_phi_34

      !! Norms of Q34
      real*8 n_Q2_34,n_E_34,n_cos_34

      !! Variables of the e+e- pair
      real*8 cos_3,phi_3,w_phi_3

      !! random numbers
      external getrnd
      external getphi
      integer i
      real*8 csi(5)

      !! Variables for born photons
      real*8 cos_g1,cos_g2,sin_g2,cos_2_34,E_g2
      real*8 w_cos_g2,n_cos_g2,phi_g2,w_phi_g2
      real*8 w_cos_g1,n_cos_g1,phi_g1,w_phi_g1

      !! Rotation
      external rotation_xy
      real*8 jacobian

      real*8 var
      real*8 eps
      common/epssoft/eps

      external dLips_2

      !! CUTS
      real*8 cos_el,cos_ep

      real*16 b1q,b2q,Qmaxq,Qminq

      !! MULTICHANNEL SAMPLING
      real*8 n_smpl
      real*8 w_1,w_2,Q2_13,Q2_24
      real*8 thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      character*2 sampling,angdip
      real*8 cos_12,sin_g1
      real*8 E_2,w_E_2,n_E_2

      !! Parameters of the resonance
      real*8 M,G
      s = dot(p1+p2,p1+p2)

      Q_34(0) = 0d0
      Q_34(1) = 0d0
      Q_34(2) = 0d0
      Q_34(3) = 0d0

      w=1d0
      w_flat=1d0
      w_smpl=1d0
      n_smpl=1d0
      norm=1d0
      phsp=1d0
      !!******************************
      !! RF = Center of Mass  !!
      !!******************************
      !! Invariant mass of the virtual photon sampled as 1/Q_34^2

      if(igen.eq.0) then

         Q_34=p3+p4
         Q2_34=dot(Q_34,Q_34)
         E_34=Q_34(0)

         beta_e = p1(3)/p1(0)

         cos_g2= q2(3)/q2(0)
         cos_g1= q1(3)/q1(0)

         cos_3 = p3(3)/sqrt(abs(tridot(p3,p3)))
         phi_3 = atan2(p3(2),p3(1))

         phi_g2 = atan2(q2(2),q2(1))
         phi_g1 = atan2(q1(2),q1(1))
         E_2 = q2(0)

      endif

      !! Parameters of the rho
      M = 0.77456d0
      G = 0.14832d0

      !! I generate the invariant mass of the p3+p4 pair
      Qmin =(m3+m4)**2
      Qmax = s-2d0*sqrt(s)*egmin

      if(sampling.eq.'BW') then
         call sampling_BW_resonances(M,G,Qmin,Qmax,Q2_34,w_Q2_34,
     >     n_Q2_34,igen)
      else if(sampling.eq.'CO') then
         call sampling_energy_mxx(Qmin,Qmax,Q2_34,w_Q2_34,n_Q2_34,igen)
      else if(sampling.eq.'SC') then
         call sampling_prop(s,Qmin,Qmax,Q2_34,w_Q2_34,n_Q2_34,igen)
      else if(sampling.eq.'FL') then
         if(igen.eq.1) then
            call getrnd(xi,1)
            Q2_34 = Qmin + xi(1)*(Qmax-Qmin)
         endif
         w_Q2_34 = 1.d0
         n_Q2_34 = 1.d0
         w_flat =w_flat*(Qmax-Qmin)
      endif

      w_smpl = w_smpl*w_Q2_34
      n_smpl = n_smpl*n_Q2_34



!! I generate the energy of the photon 2 as 1/E_2
      Emin=eps*sqrt(s)/2d0
      Emax=sqrt(s)/2d0

      call sampling_energy(Emin,Emax,E_2,w_E_2,n_E_2,igen)
      w_smpl = w_smpl*w_E_2
      n_smpl = n_smpl*n_E_2

!! Cos(th_34) of the virtual photon sampled as 1/(1-beta_e^2 beta_34^2 cos^2(th_34))
      if(igen.eq.1) then
         beta_e = p1(3)/p1(0)
      endif


      if(angdip.eq.'co') then
         call sampling_cos_collsq(beta_e,1d0,-1d0,1d0,
     .        cos_g1,w_cos_g1,n_cos_g1,igen)

      else if(angdip.eq.'fl') then
         if(igen.eq.1) then
            call getrnd(xi,1)
            cos_g1 = -1.d0 + 2.d0*xi(1)
         endif

         w_cos_g1 = 1.d0
         n_cos_g1 = 1.d0
         w_flat  = w_flat *  2.d0

      endif

      w_smpl = w_smpl*w_cos_g1
      n_smpl = n_smpl*n_cos_g1

      call sampling_cos_collsq(beta_e,1d0,-1.d0,1.d0,
     .                  cos_g2,w_cos_g2,n_cos_g2,igen)
      w_smpl = w_smpl*w_cos_g2
      n_smpl = n_smpl*n_cos_g2


!! until now I have generated 4 sampled variables
      w=w*w_smpl
      norm=n_smpl

      if(igen.eq.1) then

         call getrnd(csi,4)

         !! Generating the azimuths of 1 and 2 flat
         call getphi(csi(1),phi_g1,w_phi_g1)
         w_flat = w_flat*w_phi_g1

         call getphi(csi(2),phi_g2,w_phi_g2)
         w_flat = w_flat*w_phi_g2

         !! At this point I have the cosines and all the angles. I can calculate the energy
         sin_g1=sqrt(1d0-cos_g1**2)
         sin_g2=sqrt(1d0-cos_g2**2)

         q2(0)=E_2
         q2(1)=E_2*sin_g2*cos(phi_g2)
         q2(2)=E_2*sin_g2*sin(phi_g2)
         q2(3)=E_2*cos_g2

         cos_12 = sin_g2*cos(phi_g2)*sin_g1*cos(phi_g1)
         cos_12 = cos_12 + sin_g1*sin(phi_g1)*sin_g2*sin(phi_g2)
         cos_12 = cos_12 + cos_g1*cos_g2

         q1(0)=(s-2d0*sqrt(s)*E_2-Q2_34)
         q1(0)=q1(0)/(2d0*(sqrt(s)-E_2*(1d0-cos_12)))
         q1(1)=q1(0)*sin_g1*cos(phi_g1)
         q1(2)=q1(0)*sin_g1*sin(phi_g1)
         q1(3)=q1(0)*cos_g1


 !     !*********************************************
      !! RF = Q_34 rest frame of the virtual photon
      !! I compute the decay g(Q_34) -> e-(p3) e+(p4)
      !!*********************************************

      !! At first I generate the cos(th_3) and phi_3 of one of the electrons
      !! NO sampling

         call getphi(csi(3),phi_3,w_phi_3)
         w_flat=w_flat*w_phi_3

         cos_3= -1d0 + 2d0*csi(4)
         w_flat=w_flat*2d0

         !! The trimomentum of q2 is fixed by the delta q2 = -q1-Q_34
         Q_34(0)=sqrt(Q2_34+abs(tridot(q1+q2,q1+q2)))

         do i =1,3
            Q_34(i)=-q2(i)-q1(i)
         enddo

         !! The virtual photon decays into e+e- pair
         call dLips_2(Q_34,cos_3,phi_3,p3,m3,p4,m4,phsp_2,icut)

!        ! 1/(2 (2 pi)^3)^3 coming from the integrals in dq1 dq2 dQ34
         phsp = (1d0/(2d0*(2d0*pi)**3)**3)

         !! 2->2 phase space
         phsp = phsp*phsp_2

         !! Jacobian of the delta^0 * change of variable form |Q_34| to E_34
         jacobian = q2(0)*q1(0)
         jacobian = jacobian/(Q_34(0)+q1(0)+q2(0)*cos_12)

         if(q1(0).le.0.d0) then
            icut=1
            phsp=0.d0
            w=0.d0
            return
         endif
         phsp = phsp*jacobian

         !! (2pi)^4 coming from the overall delta
         !! (2pi)^3 coming from the Q^2_34 integral
         phsp = phsp*(2d0*pi)**4*(2d0*pi)**3
         phsp = phsp*w_flat!*w_smpl*n_smpl

c$$$         if(q2(0).ge.q1(0)) then
c$$$            phsp=0.d0
c$$$            w=0.d0
c$$$            icut=1
c$$$         else
c$$$            phsp=2.d0*phsp
c$$$         endif

      endif


      return
      end

********************


*************FLAT
        !! 2->4 phase space as diagram (A) of the notes
      !!      \/
      !!      *
      !! ----*----*~~~~~~~~
      !!          |
      !!          |
      !! ---------*~~~~~~~~
      subroutine CHN_2g_flat_ISR(p1,p2,p3,m3,p4,m4,
     . q1,q2,w,norm,phsp,icut,igen)
*******************************************************************
      implicit none
      dimension p1(0:3),p2(0:3)
      dimension p3(0:3),p4(0:3),Q_34(0:3),q1(0:3),q2(0:3)
      real*8 Q_34
      integer icut
      real*8 p1,p2,p3,p4,q1,q2,w,norm
      real*8 Qmin,Qmax,emax
      real*8 cmin_34,cmax_34,cmin_g,cmax_g
      real*8 cmax,cmin
      real*8 qph(40,0:3),xi(4)
      common/coseni/cmax,cmin
      real*8 ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
      real*8 s
      real*8 tridot,dot
      external tridot,dot
      real*8 mod234,phsp,phsp_2
      real*8 m3,m4
      logical cond
      integer igen

      !! Weights: flat and sampled
      real*8 w_flat, w_smpl

      !! Variables of Q34
      real*8 Q2_34,E_34,cos_34,sin_34,phi_34,modQ34
      real*8 beta_e,beta_34

      !! Weights of Q34
      real*8 w_Q2_34,w_E_34,w_cos_34,w_phi_34

      !! Norms of Q34
      real*8 n_Q2_34,n_E_34,n_cos_34

      !! Variables of the e+e- pair
      real*8 cos_3,phi_3,w_phi_3

      !! random numbers
      external getrnd
      external getphi
      integer i
      real*8 csi(5)

      !! Variables for born photons
      real*8 cos_g1,cos_g2,sin_g2,cos_2_34,E_g2
      real*8 w_cos_g2,n_cos_g2,phi_g2,w_phi_g2
      real*8 w_cos_g1,n_cos_g1,phi_g1,w_phi_g1

      !! Rotation
      external rotation_xy
      real*8 jacobian

      real*8 var
      real*8 eps
      common/epssoft/eps

      external dLips_2

      !! CUTS
      real*8 cos_el,cos_ep

      real*16 b1q,b2q,Qmaxq,Qminq

      !! MULTICHANNEL SAMPLING
      real*8 n_smpl
      real*8 w_1,w_2,Q2_13,Q2_24
      real*8 thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      character*2 sampling
      real*8 cos_12,sin_g1
      real*8 E_2,w_E_2,n_E_2

      !! Parameters of the resonance
      real*8 M,G
      s = dot(p1+p2,p1+p2)

      Q_34(0) = 0d0
      Q_34(1) = 0d0
      Q_34(2) = 0d0
      Q_34(3) = 0d0

      w=1d0
      w_flat=1d0
      w_smpl=1d0
      n_smpl=1d0
      norm=1d0
      phsp=1d0
      !!******************************
      !! RF = Center of Mass  !!
      !!******************************
      !! Invariant mass of the virtual photon sampled as 1/Q_34^2

      if(igen.eq.0) then

         Q_34=p3+p4
         Q2_34=dot(Q_34,Q_34)
         E_34=Q_34(0)

         beta_e = p1(3)/p1(0)

         cos_g2= q2(3)/q2(0)
         cos_g1= q1(3)/q1(0)

         cos_3 = p3(3)/sqrt(abs(tridot(p3,p3)))
         phi_3 = atan2(p3(2),p3(1))

         phi_g2 = atan2(q2(2),q2(1))
         phi_g1 = atan2(q1(2),q1(1))
         E_2 = q2(0)

      endif

      !! Parameters of the rho
      M = 0.77456d0
      G = 0.14832d0


      if(igen.eq.1) then
         call getrnd(xi,4)

!! I generate the invariant mass of the p3+p4 pair
         Qmin =(m3+m4)**2
         Qmax = s-2d0*sqrt(s)*egmin

         Q2_34 = Qmin+xi(1)*(Qmax-Qmin)

         w_flat = w_flat*(Qmax-Qmin)

!! I generate the energy of the photon 2 as 1/E_2
         Emin=eps*sqrt(s)/2d0
         Emax=sqrt(s)/2d0

         E_2=Emin+xi(2)*(Emax-Emin)

         w_flat = w_flat*(Emax-Emin)

!     ! Cos(th_34) of the virtual photon sampled as 1/(1-beta_e^2 beta_34^2 cos^2(th_34))
         if(igen.eq.1) then
            beta_e = p1(3)/p1(0)
         endif
!! I generate cos_1 and cos_2 of the two photons as collinear
!! Sampling with (1-b^2 cos^2)
         cmin_g = -1d0
         cmax_g = 1d0

         cos_g1=-1.d0+2.d0*xi(3)
         cos_g2=-1.d0+2.d0*xi(4)

!     ! until now I have generated 4 sampled variables
         w_flat=w_flat*4.d0

         call getrnd(csi,4)

         !! Generating the azimuths of 1 and 2 flat
         call getphi(csi(1),phi_g1,w_phi_g1)
         w_flat = w_flat*w_phi_g1

         call getphi(csi(2),phi_g2,w_phi_g2)
         w_flat = w_flat*w_phi_g2

         !! At this point I have the cosines and all the angles. I can calculate the energy
         sin_g1=sqrt(1d0-cos_g1**2)
         sin_g2=sqrt(1d0-cos_g2**2)

         q2(0)=E_2
         q2(1)=E_2*sin_g2*cos(phi_g2)
         q2(2)=E_2*sin_g2*sin(phi_g2)
         q2(3)=E_2*cos_g2

         cos_12 = sin_g2*cos(phi_g2)*sin_g1*cos(phi_g1)
         cos_12 = cos_12 + sin_g1*sin(phi_g1)*sin_g2*sin(phi_g2)
         cos_12 = cos_12 + cos_g1*cos_g2

         q1(0)=(s-2d0*sqrt(s)*E_2-Q2_34)
         q1(0)=q1(0)/(2d0*(sqrt(s)-E_2*(1d0-cos_12)))
         q1(1)=q1(0)*sin_g1*cos(phi_g1)
         q1(2)=q1(0)*sin_g1*sin(phi_g1)
         q1(3)=q1(0)*cos_g1


 !     !*********************************************
      !! RF = Q_34 rest frame of the virtual photon
      !! I compute the decay g(Q_34) -> e-(p3) e+(p4)
      !!*********************************************

      !! At first I generate the cos(th_3) and phi_3 of one of the electrons
      !! NO sampling

         call getphi(csi(3),phi_3,w_phi_3)
         w_flat=w_flat*w_phi_3

         cos_3= -1d0 + 2d0*csi(4)
         w_flat=w_flat*2d0

         !! The trimomentum of q2 is fixed by the delta q2 = -q1-Q_34
         Q_34(0)=sqrt(Q2_34+abs(tridot(q1+q2,q1+q2)))

         do i =1,3
            Q_34(i)=-q2(i)-q1(i)
         enddo

         !! The virtual photon decays into e+e- pair
         call dLips_2(Q_34,cos_3,phi_3,p3,m3,p4,m4,phsp_2,icut)

!        ! 1/(2 (2 pi)^3)^3 coming from the integrals in dq1 dq2 dQ34
         phsp = (1d0/(2d0*(2d0*pi)**3)**3)

         !! 2->2 phase space
         phsp = phsp*phsp_2

         !! Jacobian of the delta^0 * change of variable form |Q_34| to E_34
         jacobian = q2(0)*q1(0)
         jacobian = jacobian/(Q_34(0)+q1(0)+q2(0)*cos_12)

         if(q1(0).le.0.d0) then
            icut=1
            phsp=0.d0
            w=0.d0
            return
         endif
         phsp = phsp*jacobian

         !! (2pi)^4 coming from the overall delta
         !! (2pi)^3 coming from the Q^2_34 integral
         phsp = phsp*(2d0*pi)**4*(2d0*pi)**3
         phsp = phsp*w_flat!*w_smpl*n_smpl


c$$$         if(q2(0).ge.q1(0)) then
c$$$            phsp=0.d0
c$$$            w=0.d0
c$$$            icut=1
c$$$         else
c$$$            phsp=2.d0*phsp
c$$$         endif

      endif

! norm=norm*w_flat
      return
      end

********************


      !! 2->4 phase space as diagram (A) of the notes
      !!      \/
      !!      *
      !! ----*----*~~~~~~~~
      !!          |
      !!          |
      !! ---------*~~~~~~~~
      subroutine CHN_2g_pair_ISR(p1,p2,p3,m3,p4,m4,
     . q1,q2,w,norm,phsp,icut,sampling,igen)
************************************************************************************
      !!  dPhi_4 = dcos_2 dphi_2 dcos_34 dphi_34 dQ_34^2 dE_34  dPhi_2(Q34 -> p3+p4)
      !!   dcos_2 -> 1/(1-b^2 cos_2^2)
      !!   dphi_2 -> flat
      !!  dcos_34 -> 1/[Q_34^2 - 2 E_1 E_34(1-b_1 b_34 cos_34)]
      !!  dphi_34 -> flat
      !!  dQ_34^2 -> sampling = 'BW' BW^2  sampling='coul' 1/Q_34^2
      !!    dE_34 -> flat
      !!   dcos_3 -> flat
      !!   dphi_3 -> flat
************************************************************************************
      implicit none
      dimension p1(0:3),p2(0:3)
      dimension p3(0:3),p4(0:3),Q_34(0:3),q1(0:3),q2(0:3)
      real*8 Q_34
      integer icut
      real*8 p1,p2,p3,p4,q1,q2,w,norm
      real*8 Qmin,Qmax,emax
      real*8 cmin_34,cmax_34,cmin_g,cmax_g
      real*8 cmax,cmin
      real*8 qph(40,0:3)
      common/coseni/cmax,cmin
      real*8 ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
      real*8 s
      real*8 tridot,dot
      external tridot,dot
      real*8 mod234,phsp,phsp_2
      real*8 m3,m4
      logical cond
      integer igen

      !! Weights: flat and sampled
      real*8 w_flat, w_smpl

      !! Variables of Q34
      real*8 Q2_34,E_34,cos_34,sin_34,phi_34,modQ34
      real*8 beta_e,beta_34

      !! Weights of Q34
      real*8 w_Q2_34,w_E_34,w_cos_34,w_phi_34

      !! Norms of Q34
      real*8 n_Q2_34,n_E_34,n_cos_34

      !! Variables of the e+e- pair
      real*8 cos_3,phi_3,w_phi_3

      !! random numbers
      external getrnd
      external getphi
      integer i
      real*8 csi(5)

      !! Variables for born photons
      real*8 cos_g2,sin_g2,cos_2_34,E_g2
      real*8 w_cos_g2,n_cos_g2,phi_g2,w_phi_g2

      !! Rotation
      external rotation_xy
      real*8 jacobian

      real*8 var

      external dLips_2

      !! CUTS
      real*8 cos_el,cos_ep

      real*16 b1q,b2q,Qmaxq,Qminq

      !! MULTICHANNEL SAMLING
      real*8 n_smpl
      real*8 w_1,w_2,Q2_13,Q2_24

      real*8 thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/epssoft/eps
      real*8 egmin_eps
      real*8 eps
      character*2 sampling

      !! Parameters of the resonance
      real*8 M,G
      s = dot(p1+p2,p1+p2)

      Q_34(0) = 0d0
      Q_34(1) = 0d0
      Q_34(2) = 0d0
      Q_34(3) = 0d0

      w=1d0
      w_flat=1d0
      w_smpl=1d0
      n_smpl=1d0
      norm=1d0
      phsp=1d0
      !!******************************
      !! RF = Center of Mass  !!
      !!******************************
      !! Invariant mass of the virtual photon sampled as 1/Q_34^2

      if(igen.eq.0) then
         ! q1(0:3)=qph(1,0:3)
         ! q2(0:3)=qph(2,0:3)
         Q_34=p3+p4
         Q2_34=dot(Q_34,Q_34)
         E_34=Q_34(0)

         beta_e = p1(3)/p1(0)
         beta_34 = sqrt(1d0-Q2_34/(E_34**2))

         phi_34 = atan2(Q_34(2),Q_34(1))
         cos_34 = Q_34(3)/sqrt(abs(tridot(Q_34,Q_34)))
         cos_g2= q2(3)/q2(0)

         cos_3 = p3(3)/sqrt(abs(tridot(p3,p3)))
         phi_3 = atan2(p3(2),p3(1))

         phi_g2 = atan2(q2(2),q2(1))
      endif

      !! Parameters of the rho
      M = 0.77456d0
      G = 0.14832d0

      Qmin =(m3+m4)**2
      Qmax = s-2d0*sqrt(s)*egmin

      if(sampling.eq.'BW') then
         call sampling_BW_resonances(M,G,Qmin,Qmax,Q2_34,w_Q2_34,
     >     n_Q2_34,igen)
      else if(sampling.eq.'CO') then
         call sampling_energy(Qmin,Qmax,Q2_34,w_Q2_34,n_Q2_34,igen)
      endif

      w_smpl= w_smpl*w_Q2_34
      n_smpl = n_smpl*n_Q2_34
      ! print*,'w_Q2_34',w_Q2_34

      !! Energy of the virtual photon sampled as 1/E_34
      if(igen.eq.1) then
         call getrnd(csi,5)
      Emin = sqrt(Q2_34)
      Emax = (s+Q2_34)/(2d0*sqrt(s))
      E_34=Emin +(Emax-Emin)*csi(1)
      w_flat=w_flat*(Emax-Emin)

      ! call sampling_energy(emin,emax,E_34,w_E_34,n_E_34,igen)
      ! w_smpl = w_smpl*w_E_34
      ! n_smpl = n_smpl*n_E_34

      endif
      !! Jacobian because I'm generating the energy

      !! Cos(th_34) of the virtual photon sampled as 1/(1-beta_e^2 beta_34^2 cos^2(th_34))
      if(igen.eq.1) then
         beta_e = p1(3)/p1(0)
         beta_34 = sqrt(1d0-Q2_34/(E_34**2))
      endif
      cmin_34 = -1d0
      cmax_34 = 1d0


      call sampling_cos_invm_quad(Q2_34,p1(0),E_34,beta_e,beta_34,
     . cmin_34,cmax_34,cos_34,w_cos_34,n_cos_34,igen)
      w_smpl = w_smpl*w_cos_34
      n_smpl = n_smpl*n_cos_34
      ! print*,'w_cos_34',w_cos_34



      !! Cos(theta)_2 of photon 2 sampled as 1/(1+beta_e^2 cos^2(theta))
      cmin_g = -1d0
      cmax_g = 1d0
      call sampling_cos_collsq(beta_e,1d0,cmin_g,cmax_g,
     .                  cos_g2,w_cos_g2,n_cos_g2,igen)
      w_smpl = w_smpl*w_cos_g2
      n_smpl = n_smpl*n_cos_g2
      w=w_smpl
      norm=n_smpl

      if(igen.eq.1) then

         call getphi(csi(2),phi_34,w_phi_34)
         w_flat = w_flat*w_phi_34
c
         sin_34 = sqrt(1d0-cos_34**2)
         modQ34 = sqrt(E_34**2-Q2_34)
         Q_34(0) = E_34
         Q_34(1) = modQ34*sin_34*cos(phi_34)
         Q_34(2) = modQ34*sin_34*sin(phi_34)
         Q_34(3) = modQ34*cos_34

      !!*********************************************
      !! RF = Q_34 rest frame of the virtual photon
      !! I compute the decay g(Q_34) -> e-(p3) e+(p4)
      !!*********************************************

      !! At first I generate the cos(th_3) and phi_3 of one of the electrons
      !! NO sampling

         call getphi(csi(3),phi_3,w_phi_3)
         w_flat=w_flat*w_phi_3

         cos_3= -1d0 + 2d0*csi(4)
         w_flat=w_flat*2d0

         !! The virtual photon decays into e+e- pair

         ! cos_g= -1d0 + 2d0*csi(4)
         ! w_flat=w_flat*2d0

         sin_g2=sqrt(1d0-cos_g2**2)
         !! I compute the relative angle between q1 and Q_34
         cos_2_34 = sin_g2*sin_34*cos(phi_34)+ cos_g2*cos_34
         !! Now I can compute the energy of q1
         E_g2 = (s-2d0*sqrt(s)*E_34+Q2_34)
         E_g2 = E_g2 /(2d0*(sqrt(s)-E_34+modQ34*cos_2_34))
         !! Check if the energy is below maximum
         cond = E_g2.gt.0d0.and.E_g2.le.(s-Q2_34)/(2d0*sqrt(s))
         if (.not.cond) then
            icut = 1
            phsp=0d0
            return
         endif

         q2(0) = E_g2
         q2(1) = E_g2*sin_g2
         q2(2) = 0d0
         q2(3) = E_g2*cos_g2

               !! The trimomentum of q2 is fixed by the delta q2 = -q1-Q_34
         do i =1,3
            q1(i)=-q2(i)-Q_34(i)
         enddo

         !! q2 is massless, the energy is given by the momentum
         mod234=sqrt(abs(tridot(q2+Q_34,q2+Q_34)))
         q1(0) = mod234

         !! I have to rotate everything from the x-z plane
         call getphi(csi(5),phi_g2,w_phi_g2)
         w_flat = w_flat*w_phi_g2

         call rotation_xy(phi_g2,q1)
         call rotation_xy(phi_g2,q2)
         call rotation_xy(phi_g2,Q_34)


         !! The virtual photon decays into e+e- pair
         call dLips_2(Q_34,cos_3,phi_3,p3,m3,p4,m4,phsp_2,icut)

         !! 1/(2 (2 pi)^3)^3 coming from the integrals in dq1 dq2 dQ34
         phsp = (1d0/(2d0*(2d0*pi)**3)**3)

         !! 2->2 phase space
         phsp = phsp*phsp_2

         !! Jacobian of the delta^0 * change of variable form |Q_34| to E_34
         jacobian = modQ34*E_g2/(E_g2+ mod234 + modQ34*cos_2_34)
         phsp = phsp*jacobian

         !! (2pi)^4 coming from the overall delta
         !! (2pi)^3 coming from the Q^2_34 integral
         phsp = phsp*(2d0*pi)**4*(2d0*pi)**3
         phsp = phsp*w_flat!*w_smpl*n_smpl

         ! qph(1,0:3)=q1(0:3)
         ! qph(2,0:3)=q2(0:3)

      endif
      ! norm=norm*w_flat

      return
      end


      !! This routine samples cos(theta) as 1/(A-B(1-C cos(theta)))
      !! input  | en: energy of the leg
      !!        | cmin cmax: max and min cos(theta)
      !!        | igen: 0 computes the weight, 1 generates variables
      !! output | cos: sampled cos(theta)
      !!        | wgt: the weight (1-b cos)
      !!        | norm: int_cmin^cmax 1/(1-b cos) dcos
      subroutine sampling_cos_invm(Q2,E1,E2,b1,b2,cmin,cmax,
     .                               cos,wgt,norm,igen)
      implicit none
      real*8 cos,wgt,cmin,cmax
      real*8 Q2,E1,E2,b1,b2
      real*8 bb,fac,norm,en,arg
      real*8 xi(1)
      real*8 logmax,logmin,EE
      real*8 ame,ammu,convfac,alpha,pi
      real*8 A,BBB,C
      integer igen
      common/parameters/ame,ammu,convfac,alpha,pi

      logmax = (Q2-2d0*E1*E2*(1d0-b1*b2*cmax))
      logmin = (Q2-2d0*E1*E2*(1d0-b1*b2*cmin))
      norm=1d0/(b1*b2*E1*E2)*(log(logmax/logmin))
      ! print*,'norm',norm

      A = Q2
      BBB = 2d0*E1*E2
      C = b1*b2
      if(igen.eq.1) then
         call getrnd(xi,1)
         arg = BBB*C*norm*xi(1)
         cos = 1d0/(BBB*C)*(BBB-A-exp(arg)*(A-BBB*(1d0-C*cmin)))
      else if(igen.eq.0) then
         cos=cos
      endif

      wgt = (A-BBB*(1-C*cos))
      return
      end

      subroutine sampling_cos_invm_quad(Q2,E1,E2,b1,b2,cmin,cmax,
     .                               cos,wgt,norm,igen)
      implicit none
      real*16 cos_q,wgt_q
      real*8 cos,wgt,cmin,cmax
      real*8 Q2,E1,E2,b1,b2,norm
      real*16 bb,fac,en,arg
      real*8 xi(1)
      real*16 logmax,logmin,EE
      real*8 ame,ammu,convfac,alpha,pi
      real*16 A,BBB,C,K,Kmin,Kmax
      integer igen
      common/parameters/ame,ammu,convfac,alpha,pi

      A = Q2*1q0
      BBB = 2q0*E1*E2
      C = b1*b2*1q0

      logmax = log(abs(A-BBB*(1q0-C*cmax)))
      logmin = log(abs(A-BBB*(1q0-C*cmin)))
      ! print*,'logmin',logmin


      Kmin = (A-BBB*(1q0-C*cmin))/(abs(A-BBB*(1q0-C*cmin)))
      Kmax = (A-BBB*(1q0-C*cmax))/(abs(A-BBB*(1q0-C*cmax)))
      norm=1q0/(BBB*C)*(-1q0*logmax + 1q0*logmin)

      if(igen.eq.1) then
         call getrnd(xi,1)
         K = norm*BBB*C*xi(1)*1q0 + Kmin*logmin
         cos_q = (-exp(-K)-A+BBB)/(BBB*C)
      else if(igen.eq.0) then
         cos_q=cos
      endif

      wgt_q = abs(A-BBB*(1q0-C*cos_q))
      wgt=wgt_q*1d0
      cos=cos_q*1d0
      return
      end


      !! this routines rotates the momentum p1
      !! by an angle phi on the plane x-y
      subroutine rotation_xy(phi,p1)
      implicit none
      real*8 phi,p1,ptmp
      dimension p1(0:3),ptmp(0:3)
      real*8 c,s
      c=cos(phi)
      s=sin(phi)
      ptmp = p1
      p1(1) = ptmp(1)*c - ptmp(2)*s
      p1(2) = ptmp(1)*s + ptmp(2)*c

      return
      end



!!!  4,5,6,7) photon collinear to the pair direction
!!!  and other photon back to back wrt pair+ph jet
      subroutine CHN_2g_ISR_FSR(p1,p2,p3,m3,p4,m4,q1,q2,w,norm,
     .                       phsp,icut,sampling,angdip,igen)
      implicit none
      dimension p1(0:3),p2(0:3)
      dimension p3(0:3),p4(0:3),q1(0:3),q2(0:3)
      integer icut,i
      character*2 sampling,angdip
      real*8 p1,p2,p3,p4,q1,q2,w,m3,m4,phsp
      real*8 Q234,Q23
      dimension Q234(0:3),Q23(0:3)
      real*8 cmin,cmax
      common/coseni/cmax,cmin
      real*8 ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
      real*8 E_gg_min
      common/pairscut/E_gg_min
      real*8 s
      real*8 tridot,dot
      external tridot,dot
      logical cond
      integer igen
      real*8 norm,xi(1),xi1(1)
      !! Variables for weight
      real*8 w_smpl,w_flat,w_cos_g1,n_cos_g1
      !! Variables for Q234 sampling
      real*8 Q234min,Q234max,Q2_234,w_Q2_234,n_Q2_234
      !! Variables for Q24 sampling
      real*8 Q23min,Q23max,Q2_23,w_Q2_23,n_Q2_23
      !! Routines needed for the generation of 6 angles
      external getrnd,getphi
      real*8 csi(6)
      !! Variables for angles generation
      real*8 cos234,cos234min,cos234max,phi234,w_phi234
      real*8 m234,m23
      real*8 cos4,cos4min,cos4max,phi4,w_phi4
      real*8 cos2,cos2min,cos2max,phi2,w_phi2
      !! 2->2 phasespaces
      real*8 phsp2_234_1,phsp2_24_4,phsp2_2_3
      !! CUTS
      real*8 cos_el,cos_ep,cos_g1,cos_g2
      real*8 w_cos234,n_cos234,beta_e
      real*8 thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/epssoft/eps
      real*8 eps
      real*8 egmin_eps
      real*8 M,G
      real*8 Q_234
      dimension q_234(0:3)
      s = dot(p1+p2,p1+p2)
      w=1d0
      w_smpl=1d0
      w_flat=1d0
      norm=1d0
      phsp=1d0
      !! Initialization of the momenta of the virtual photon-fermions

      beta_e = sqrt(1d0-ame**2/p1(0)**2)
      if(igen.eq.0) then
         Q2_234=dot(q2+p3+p4,q2+p3+p4)
         Q2_23=dot(q2+p3,q2+p3)

         Q234=q2+p3+p4
         Q23=q2+p3
         cos234=Q234(3)/sqrt(abs(tridot(Q234,Q234)))
      endif

!! FIRST DECAY
!! p1+p2 -> q1 + Q234
!! I sample the invariant mass of the first fermion-gamma pair as 1/Q^2

      M = 0.77456d0
      G = 0.14832d0


      Q234min = (m3+m4)**2  !!! In reality it's s-2(sqrt(s)*E_1^max)
      Q234max = s*(1.d0-eps)       !!! The photon E_1^min = 0


      if(sampling.eq.'BW') then
        call sampling_BW_resonances(M,G,Q234min,Q234max,Q2_234,w_Q2_234,
     >     n_Q2_234,igen)

      else if(sampling.eq.'CO') then
         call sampling_energy_mxx(Q234min,Q234max,Q2_234,w_Q2_234,
     >        n_Q2_234,igen)

      else if(sampling.eq.'SC') then
         call sampling_prop(s,Q234min,Q234max,Q2_234,w_Q2_234,
     >        n_Q2_234,igen)

      else if(sampling.eq.'FL') then
         if(igen.eq.1) then
            call getrnd(xi,1)
            Q2_234 = Q234min + xi(1)*(Q234max-Q234min)
         endif

         w_Q2_234=1.d0
         n_Q2_234=1.d0

         w_flat = w_flat*(Q234max-Q234min)
      endif

      w_smpl = w_smpl* w_Q2_234
      norm   = norm *n_Q2_234
      m234=dsqrt(Q2_234)

!! SECOND DECAY
!! Q234 -> p4 + Q23
!! I sample the invariant mass of the second fermion-gamma pair as 1/Q^2
c$$$      Q23min = m3**2
c$$$      Q23max =(Q2_234)
c$$$      call sampling_energy_mxx(Q23min,Q23max,Q2_23,w_Q2_23,n_Q2_23,igen)

      Q23min=m3**2+eps*s/2.d0*(1-sqrt(1.d0-(2.d0*m3)**2/s))

!     Q23min=m3**2*(1.d0+1.d-10)
      Q23max =(Q2_234)
      call sampling_energy_mxx_sc(m3,Q23min,Q23max,Q2_23,w_Q2_23,n_Q2_23
     >     ,igen)

      m23=sqrt(Q2_23)

      w_smpl= w_smpl*w_Q2_23
      norm = norm*n_Q2_23


      if(angdip.eq.'co') then
         call sampling_cos_collsq(beta_e,1d0,-1d0,1d0,
     .        cos234,w_cos234,n_cos234,igen)

      else if(angdip.eq.'fl') then
         if(igen.eq.1) then
            call getrnd(xi1,1)
            cos234 = -1.d0 + 2.d0*xi1(1)
         endif
         w_cos234 = 1.d0
         n_cos234 = 1.d0
         w_flat  = w_flat *  2.d0

      endif

      w_smpl= w_smpl* w_cos234
      norm  = norm  * n_cos234


      w=w_smpl

      if(igen.eq.1) then

         call getrnd(csi,6)

         Q234(0)=0d0
         Q234(1)=0d0
         Q234(2)=0d0
         Q234(3)=0d0

         Q23(0)=0d0
         Q23(1)=0d0
         Q23(2)=0d0
         Q23(3)=0d0

c$$$         cos234min=-1d0
c$$$         cos234max=1d0
c$$$         cos234 = cos234min + (cos234max-cos234min)*csi(1)
c$$$         w_flat = w_flat * (cos234max-cos234min)

         call getphi(csi(2),phi234,w_phi234)
         w_flat = w_flat*w_phi234

         call dLips_2(p1+p2,cos234,phi234,q1,0.d0,Q234,m234,
     .               phsp2_234_1,icut)
c$$$         do i=1,3
c$$$            if(abs(q1(i))-q1(0).le.1.d-14) then
c$$$               phsp=0.d0
c$$$               w=0.d0
c$$$               icut=1
c$$$               return
c$$$            endif
c$$$         enddo
c$$$
         cos4min=-1d0
         cos4max=1d0
         cos4 = cos4min + (cos4max-cos4min)*csi(3)
         w_flat = w_flat*(cos4max-cos4min)

         call getphi(csi(4),phi4,w_phi4)
         w_flat = w_flat*w_phi4

         call dLips_2(Q234,cos4,phi4,Q23,m23,p4,m4,phsp2_24_4,icut)

         cos2min=-1d0
         cos2max=1d0
         cos2 = cos2min + (cos2max-cos2min)*csi(5)
         w_flat = w_flat * (cos2max-cos2min)

         call getphi(csi(6),phi2,w_phi2)
         w_flat = w_flat*w_phi2

         call dLips_2(Q23,cos2,phi2,q2,0d0,p3,m3,phsp2_2_3,icut)

         !! 2pis -> 4 coming from overall delta, 2x3 coming from Q^2 generation
         phsp=(2d0*pi)**4*(2d0*pi)**6
         phsp=phsp*phsp2_234_1*phsp2_24_4*phsp2_2_3
         phsp=phsp*w_flat!*w_smpl*norm

c$$$         if(q2(0).ge.q1(0)) then
c$$$            phsp=0.d0
c$$$            w=0.d0
c$$$            icut=1
c$$$         else
c$$$            phsp=2.d0*phsp
c$$$         endif

      endif


      return
      end


      subroutine CHN_ng_soft_ISR(ng,p1,p2,p3,m3,p4,m4,
     . qph,w,norm,phsp,icut,sampling,angdip,igen)
*******************************************************************
      implicit none
      dimension p1(0:3),p2(0:3),mom_cons(0:3)
      dimension p3(0:3),p4(0:3),Q_34(0:3),q1(0:3),q2(0:3)
      real*8 Q_34,mom_cons
      integer icut,ng
      real*8 p1,p2,p3,p4,q1,q2,w,norm
      real*8 Qmin,Qmax,emax
      real*8 cmin_34,cmax_34,cmin_g,cmax_g,sin_ng,suminc
      real*8 cmax,cmin,sini,sinj,sinng,sumij,sumin,sumen,sumng
      real*8 qph(40,0:3),En(40),qcos(40),qphi(40),cosij(40,40)
      real*8 qval(40,4)
      common/coseni/cmax,cmin
      real*8 ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
      real*8 s
      real*8 tridot,dot
      external tridot,dot
      real*8 mod234,phsp,phsp_2
      real*8 m3,m4
      logical cond
      integer igen

      !! Weights: flat and sampled
      real*8 w_flat, w_smpl

      !! Variables of Q34
      real*8 Q2_34,E_34,cos_34,sin_34,phi_34,modQ34
      real*8 beta_e,beta_34

      !! Weights of Q34
      real*8 w_Q2_34,w_E_34,w_cos_34,w_phi_34

      !! Norms of Q34
      real*8 n_Q2_34,n_E_34,n_cos_34

      !! Variables of the e+e- pair
      real*8 cos_3,phi_3,w_phi_3

      !! random numbers
      external getrnd
      external getphi
      integer i,j
      real*8 csi(40),xi(2),xi1(1)

      !! Variables for born photons
      real*8 cos_g1,cos_g2,sin_g2,cos_2_34,E_g2
      real*8 w_cos_g2,n_cos_g2,phi_g2,w_phi_g2
      real*8 w_cos_g1,n_cos_g1,phi_g1,w_phi_g1

      !! Rotation
      external rotation_xy
      real*8 jacobian

      real*8 var
      real*8 eps
      common/epssoft/eps

      external dLips_2

      !! CUTS
      real*8 cos_el,cos_ep

      real*16 b1q,b2q,Qmaxq,Qminq

      !! MULTICHANNEL SAMPLING
      real*8 n_smpl
      real*8 w_1,w_2,Q2_13,Q2_24
      real*8 thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      character*2 sampling,angdip
      real*8 cos_12,sin_g1
      real*8 E_2,w_E_2,n_E_2

      !! Parameters of the resonance
      real*8 M,G
      s = dot(p1+p2,p1+p2)

      Q_34(0) = 0d0
      Q_34(1) = 0d0
      Q_34(2) = 0d0
      Q_34(3) = 0d0

      w=1d0
      w_flat=1d0
      w_smpl=1d0
      n_smpl=1d0
      norm=1d0
      phsp=1d0
      !!******************************
      !! RF = Center of Mass  !!
      !!******************************
      !! Invariant mass of the virtual photon sampled as 1/Q_34^2
      if(igen.eq.0) then
         Q2_34=dot(p3+p4,p3+p4)

         do i=1,ng-1
            En(i)=qph(i,0)
         enddo

         beta_e = p1(3)/p1(0)
         do i=1,ng
            qcos(i)=qph(i,3)/qph(i,0)
         enddo
      endif

      if(igen.eq.1) then
         beta_e = p1(3)/p1(0)
      endif


      !! Parameters of the rho
      M = 0.77456d0
      G = 0.14832d0

      !! I generate the invariant mass of the p3+p4 pair
      Qmin =(m3+m4)**2
      Qmax = s-2d0*sqrt(s)*egmin

      if(sampling.eq.'BW') then
         call sampling_BW_resonances(M,G,Qmin,Qmax,Q2_34,w_Q2_34,
     >     n_Q2_34,igen)

      else if(sampling.eq.'CO') then
         call sampling_energy_mxx(Qmin,Qmax,Q2_34,w_Q2_34,n_Q2_34,igen)

      else if(sampling.eq.'SC') then
         call sampling_prop(s,Qmin,Qmax,Q2_34,w_Q2_34,n_Q2_34,igen)

      else if(sampling.eq.'FL') then
         if(igen.eq.1) then
            call getrnd(xi,1)
            Q2_34 = Qmin + xi(1)*(Qmax-Qmin)
         endif

         w_Q2_34=1.d0
         n_Q2_34=1.d0

         w_flat = w_flat*(Qmax-Qmin)
      endif

      w_smpl = w_smpl*w_Q2_34
      n_smpl = n_smpl*n_Q2_34

c$$$
c$$$      call getrnd(xi,4)
c$$$      Q2_34=Qmin+xi(1)*(Qmax-Qmin)
c$$$      w_smpl=w_smpl*(Qmax-Qmin)

!! I generate the energy of the photon 2 as 1/E_2
      do i=1,ng-1

         Emin=eps*sqrt(s)/2d0
         Emax=sqrt(s)/2d0

         call sampling_energy(Emin,Emax,En(i),w_E_2,n_E_2,igen)

!         print*,'EN',En(i)
!         print*,'peso',w_E_2
         w_smpl = w_smpl*w_E_2
         n_smpl = n_smpl*n_E_2

      enddo
      !! I generate cos_1 and cos_2 of the two photons as collinear
!! Sampling with (1-b^2 cos^2)

      cmin_g = -1d0
      cmax_g = 1d0
      do i=1,ng-1

         call sampling_cos_collsq(beta_e,1d0,cmin_g,cmax_g,
     .        qcos(i),w_cos_g1,n_cos_g1,igen)


         w_smpl = w_smpl*w_cos_g1
         n_smpl = n_smpl*n_cos_g1

      enddo
!!!until now I have generated 4 sampled variables
      w=w*w_smpl
      norm=n_smpl

*************LAST PHOTON FLAT
      if(angdip.eq.'co') then
         call sampling_cos_collsq(beta_e,1d0,-1d0,1d0,
     .        qcos(ng),w_cos_g1,n_cos_g1,igen)

      else if(angdip.eq.'fl') then
         if(igen.eq.1) then
            call getrnd(xi1,1)
            qcos(ng) = -1.d0 + 2.d0*xi1(1)
         endif
         w_cos_g1 = 1.d0
         n_cos_g1 = 1.d0
         w_flat  = w_flat *  2.d0

      endif

      w = w*w_cos_g1
      norm = norm*n_cos_g1
c$$$
c$$$      if(igen.eq.1) then
c$$$         call getrnd(xi1,1)
c$$$         qcos(ng)=-1.d0+2.d0*xi1(1)
c$$$      endif
c$$$
c$$$      w_flat=w_flat*2.d0


      if(igen.eq.1) then

         call getrnd(csi,40)

         do i=1,ng
!! Generating the azimuths of 1 and 2 flat
            call getphi(csi(i),phi_g1,w_phi_g1)
            w_flat = w_flat*w_phi_g1
            qphi(i)= phi_g1
         enddo

!! At this point I have the cosines and all the angles. I can calculate the momenta of the first n-1 photon
         do i=1,ng-1

            sini=sqrt(1d0-qcos(i)**2)

            qph(i,0)=En(i)
            qph(i,1)=En(i)*sini*cos(qphi(i))
            qph(i,2)=En(i)*sini*sin(qphi(i))
            qph(i,3)=En(i)*qcos(i)

         enddo

         do i=1,ng
            do j=1,ng
              sini=dsqrt(1.d0-qcos(i)**2)
              sinj=dsqrt(1.d0-qcos(j)**2)
              cosij(i,j)= qcos(i)*qcos(j)
              cosij(i,j)= cosij(i,j)+sini*cos(qphi(i))*sinj*cos(qphi(j))
              cosij(i,j)= cosij(i,j)+sini*sin(qphi(i))*sinj*sin(qphi(j))
            enddo
         enddo

******QUI ORA RICAVO K_N E DI CONSEGUENZA Q_34
         sumen=0.d0
         do i=1,ng-1
            sumen=sumen+qph(i,0)
         enddo

         sumij=0.d0
         do i=1,ng-1
            do j=1,ng-1
               sumij=sumij+qph(i,0)*qph(j,0)*(1.d0-cosij(i,j))
            enddo
         enddo

         sumin=0.d0
         do i=1,ng-1
            sumin=sumin+qph(i,0)*(1.d0-cosij(i,ng))
         enddo

         qph(ng,0)=s-Q2_34-2.d0*sqrt(s)*sumen+sumij
         qph(ng,0)=qph(ng,0)/(2.d0*(sqrt(s)-sumin))
         if(qph(ng,0).le.0.d0) then
            icut=1
            phsp=0.d0
            w=0.d0
            return
         endif

         sin_ng=dsqrt(1.d0-qcos(ng)**2)
         qph(ng,1)=qph(ng,0)*sin_ng*cos(qphi(ng))
         qph(ng,2)=qph(ng,0)*sin_ng*sin(qphi(ng))
         qph(ng,3)=qph(ng,0)*qcos(ng)

         Q_34(0)=sqrt(s)
         Q_34(1)=0.d0
         Q_34(2)=0.d0
         Q_34(3)=0.d0

         do i=1,ng
            Q_34(0)=Q_34(0)-qph(i,0)
            do j=1,3
               Q_34(j)=Q_34(j)-qph(i,j)
            enddo
         enddo
         if(Q_34(0).le.0.d0) then
            icut=1
            phsp=0.d0
            w=0.d0
            return
         endif

c$$$  print*,Q_34

         call getrnd(xi,2)

         call getphi(xi(1),phi_3,w_phi_3)
         w_flat=w_flat*w_phi_3

         cos_3= -1d0 + 2d0*xi(2)
         w_flat=w_flat*2d0


         !! The virtual photon decays into e+e- pair
         call dLips_2(Q_34,cos_3,phi_3,p3,m3,p4,m4,phsp_2,icut)

!        ! 1/(2 (2 pi)^3)^3 coming from the integrals in dq1 dq2 dQ34
         phsp = (1d0/(2d0*(2d0*pi)**3)**(float(ng)+1.d0))

         !! 2->2 phase space
         phsp = phsp*phsp_2

         !! Jacobian of the delta^0 * change of variable form |Q_34| to E_34
         suminc = 0.d0
         do i=1,ng-1
            suminc=suminc+qph(i,0)*cosij(i,ng)
         enddo

         jacobian = 1.d0
         do i=1,ng
            jacobian=jacobian*qph(i,0)
         enddo

         jacobian = jacobian/(Q_34(0)+qph(ng,0)+suminc)

         phsp = phsp*jacobian

         !! (2pi)^4 coming from the overall delta
         !! (2pi)^3 coming from the Q^2_34 integral
         phsp = phsp*(2d0*pi)**4*(2d0*pi)**3
         phsp = phsp*w_flat!*w_smpl*n_smpl

      endif

c$$$      call prop_value(p1,p2,p3,p4,qph,ng,qval)
c$$$      cond=.false.
c$$$      do i=1,ng-1
c$$$         if(qval(ng).le.qval(i)) then
c$$$            cond=.true.
c$$$         endif
c$$$      enddo

      cond=.false.
      do i=1,ng-1
         if(qph(ng,0).le.qph(i,0)) then
            cond=.true.
         endif
      enddo

      if(cond) then
         phsp=0.d0
         icut=1
         w=0.d0
      else
         phsp=phsp*float(ng)
      endif

****************************
!      w=w_Q2_34
****************************
! norm=norm*w_flat
      return
      end

********************


      subroutine CHN_ng_flat(ng,p1,p2,p3,m3,p4,m4,
     . qph,w,norm,phsp,icut,igen)
*******************************************************************
      implicit none
      dimension p1(0:3),p2(0:3),mom_cons(0:3)
      dimension p3(0:3),p4(0:3),Q_34(0:3),q1(0:3),q2(0:3)
      real*8 Q_34,mom_cons
      integer icut,ng
      real*8 p1,p2,p3,p4,q1,q2,w,norm
      real*8 Qmin,Qmax,emax
      real*8 cmin_34,cmax_34,cmin_g,cmax_g,sin_ng,suminc
      real*8 cmax,cmin,sini,sinj,sinng,sumij,sumin,sumen,sumng
      real*8 qph(40,0:3),En(40),qcos(40),qphi(40),cosij(40,40)
      common/coseni/cmax,cmin
      real*8 ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
      real*8 s
      real*8 tridot,dot
      external tridot,dot
      real*8 mod234,phsp,phsp_2
      real*8 m3,m4
      logical cond
      integer igen

      !! Weights: flat and sampled
      real*8 w_flat, w_smpl

      !! Variables of Q34
      real*8 Q2_34,E_34,cos_34,sin_34,phi_34,modQ34
      real*8 beta_e,beta_34

      !! Weights of Q34
      real*8 w_Q2_34,w_E_34,w_cos_34,w_phi_34

      !! Norms of Q34
      real*8 n_Q2_34,n_E_34,n_cos_34

      !! Variables of the e+e- pair
      real*8 cos_3,phi_3,w_phi_3

      !! random numbers
      external getrnd
      external getphi
      integer i,j
      real*8 csi(40),xi(2),xi1(40),xi2(1),xi3(40),xi4(40)

      !! Variables for born photons
      real*8 cos_g1,cos_g2,sin_g2,cos_2_34,E_g2
      real*8 w_cos_g2,n_cos_g2,phi_g2,w_phi_g2
      real*8 w_cos_g1,n_cos_g1,phi_g1,w_phi_g1

      !! Rotation
      external rotation_xy
      real*8 jacobian

      real*8 var
      real*8 eps
      common/epssoft/eps

      external dLips_2

      !! CUTS
      real*8 cos_el,cos_ep

      real*16 b1q,b2q,Qmaxq,Qminq

      !! MULTICHANNEL SAMPLING
      real*8 n_smpl
      real*8 w_1,w_2,Q2_13,Q2_24
      real*8 thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      character*2 sampling
      real*8 cos_12,sin_g1
      real*8 E_2,w_E_2,n_E_2

      !! Parameters of the resonance
      real*8 M,G
      s = dot(p1+p2,p1+p2)

      Q_34(0) = 0d0
      Q_34(1) = 0d0
      Q_34(2) = 0d0
      Q_34(3) = 0d0

      w=1d0
      w_flat=1d0
      w_smpl=1d0
      n_smpl=1d0
      norm=1d0
      phsp=1d0
      !!******************************
      !! RF = Center of Mass  !!
      !!******************************
      !! Invariant mass of the virtual photon sampled as 1/Q_34^2
      if(igen.eq.0) then
         Q2_34=dot(p3+p4,p3+p4)

         do i=1,ng-1
            En(i)=qph(i,0)
         enddo

         beta_e = p1(3)/p1(0)
         do i=1,ng
            qcos(i)=qph(i,3)/qph(i,0)
         enddo
      endif

      if(igen.eq.1) then
         beta_e = p1(3)/p1(0)
      endif


      if(igen.eq.1) then

!! I generate the invariant mass of the p3+p4 pair
         Qmin =(m3+m4)**2
         Qmax = s-2d0*sqrt(s)*egmin

         call getrnd(xi2,1)

         Q2_34=Qmin+xi2(1)*(Qmax-Qmin)

         w_flat = w_flat*(Qmax-Qmin)

         call getrnd(xi3,40)

!     ! I generate the energy of the photon 2 as 1/E_2
         do i=1,ng-1

            Emin=eps*sqrt(s)/2d0
            Emax=sqrt(s)/2d0

            EN(i)=Emin+xi3(i)*(Emax-Emin)
!     print*,'peso',w_E_2
            w_flat = w_flat*(Emax-Emin)
         enddo
!! I generate cos_1 and cos_2 of the two photons as collinear
!     ! Sampling with (1-b^2 cos^2)

         cmin_g = -1d0
         cmax_g = 1d0
         call getrnd(xi4,40)

         do i=1,ng
            qcos(i)=-1.d0+2.d0*xi4(i)

            w_flat=w_flat*2.d0

         enddo

         call getrnd(csi,40)

         do i=1,ng
!! Generating the azimuths of 1 and 2 flat
            call getphi(csi(i),phi_g1,w_phi_g1)
            w_flat = w_flat*w_phi_g1
            qphi(i)= phi_g1
         enddo

!! At this point I have the cosines and all the angles. I can calculate the momenta of the first n-1 photon
         do i=1,ng-1

            sini=sqrt(1d0-qcos(i)**2)

            qph(i,0)=En(i)
            qph(i,1)=En(i)*sini*cos(qphi(i))
            qph(i,2)=En(i)*sini*sin(qphi(i))
            qph(i,3)=En(i)*qcos(i)

         enddo

         do i=1,ng
            do j=1,ng
              sini=dsqrt(1.d0-qcos(i)**2)
              sinj=dsqrt(1.d0-qcos(j)**2)
              cosij(i,j)= qcos(i)*qcos(j)
              cosij(i,j)= cosij(i,j)+sini*cos(qphi(i))*sinj*cos(qphi(j))
              cosij(i,j)= cosij(i,j)+sini*sin(qphi(i))*sinj*sin(qphi(j))
            enddo
         enddo

******QUI ORA RICAVO K_N E DI CONSEGUENZA Q_34
         sumen=0.d0
         do i=1,ng-1
            sumen=sumen+qph(i,0)
         enddo

         sumij=0.d0
         do i=1,ng-1
            do j=1,ng-1
               sumij=sumij+qph(i,0)*qph(j,0)*(1.d0-cosij(i,j))
            enddo
         enddo

         sumin=0.d0
         do i=1,ng-1
            sumin=sumin+qph(i,0)*(1.d0-cosij(i,ng))
         enddo

         qph(ng,0)=s-Q2_34-2.d0*sqrt(s)*sumen+sumij
         qph(ng,0)=qph(ng,0)/(2.d0*(sqrt(s)-sumin))
         if(qph(ng,0).le.0.d0) then
            icut=1
            phsp=0.d0
            w=0.d0
            return
         endif

         sin_ng=dsqrt(1.d0-qcos(ng)**2)
         qph(ng,1)=qph(ng,0)*sin_ng*cos(qphi(ng))
         qph(ng,2)=qph(ng,0)*sin_ng*sin(qphi(ng))
         qph(ng,3)=qph(ng,0)*qcos(ng)

         Q_34(0)=sqrt(s)
         Q_34(1)=0.d0
         Q_34(2)=0.d0
         Q_34(3)=0.d0

         do i=1,ng
            Q_34(0)=Q_34(0)-qph(i,0)
            do j=1,3
               Q_34(j)=Q_34(j)-qph(i,j)
            enddo
         enddo
         if(Q_34(0).le.0.d0) then
            icut=1
            phsp=0.d0
            w=0.d0
            return
         endif

c$$$  print*,Q_34

         call getrnd(xi,2)

         call getphi(xi(1),phi_3,w_phi_3)
         w_flat=w_flat*w_phi_3

         cos_3= -1d0 + 2d0*xi(2)
         w_flat=w_flat*2d0


         !! The virtual photon decays into e+e- pair
         call dLips_2(Q_34,cos_3,phi_3,p3,m3,p4,m4,phsp_2,icut)
!        ! 1/(2 (2 pi)^3)^3 coming from the integrals in dq1 dq2 dQ34
         phsp = (1d0/(2d0*(2d0*pi)**3)**(ng+1.d0))
         !! 2->2 phase space
         phsp = phsp*phsp_2
         !! Jacobian of the delta^0 * change of variable form |Q_34| to E_34
         suminc = 0.d0
         do i=1,ng-1
            suminc=suminc+qph(i,0)*cosij(i,ng)
         enddo

         jacobian = 1.d0
         do i=1,ng

            jacobian=jacobian*qph(i,0)
         enddo
         jacobian = jacobian/(Q_34(0)+qph(ng,0)+suminc)

         phsp = phsp*jacobian
         !! (2pi)^4 coming from the overall delta
         !! (2pi)^3 coming from the Q^2_34 integral
         phsp = phsp*(2d0*pi)**4*(2d0*pi)**3
         phsp = phsp*w_flat!*w_smpl*n_smpl

      endif

      cond=.false.
      do i=1,ng-1
         if(qph(ng,0).le.qph(i,0)) then
            cond=.true.
         endif
      enddo

      if(cond) then
         phsp=0.d0
         icut=1
         w=0.d0
         return
      else
         phsp=phsp*float(ng)
      endif
!     print*,norm
!     norm=norm*w_flat
      return
      end

********************

      subroutine CHN_ng_ISR_FSR(ng,p1,p2,p3,m3,p4,m4,
     . qph,w,norm,phsp,icut,sampling,angdip,igen)
*******************************************************************
      implicit none
      dimension p1(0:3),p2(0:3),mom_cons(0:3)
      dimension p3(0:3),p4(0:3),Q_34(0:3),q1(0:3),q2(0:3)
      real*8 Q_34,mom_cons
      integer icut,ng,ngi
      real*8 p1,p2,p3,p4,q1,q2,w,norm
      real*8 Qmin,Qmax,emax
      real*8 cmin_34,cmax_34,cmin_g,cmax_g,sin_ng,suminc
      real*8 cmax,cmin,sini,sinj,sinng,sumij,sumin,sumen,sumng
      real*8 qph(40,0:3),En(40),qcos(40),qphi(40),cosij(40,40)
      real*8 qval(40,4),Q23(0:3)
      common/coseni/cmax,cmin
      real*8 ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
      real*8 s
      real*8 tridot,dot
      external tridot,dot
      real*8 mod234,phsp,phsp_2
      real*8 m3,m4
      logical cond
      integer igen

      real*8 m23,Q2_23, Q23min,Q23Max, W_q2_23,n_q2_23
      real*8 phsp2_2_3,phsp2_24_4
      real*8 cos2,w_cos2,cos4,w_cos4,phi4,w_phi4,phi2,w_phi2
!! Weights: flat and sampled
      real*8 w_flat, w_smpl

      !! Variables of Q34
      real*8 Q2_34,E_34,cos_34,sin_34,phi_34,modQ34
      real*8 beta_e,beta_34

      !! Weights of Q34
      real*8 w_Q2_34,w_E_34,w_cos_34,w_phi_34

      !! Norms of Q34
      real*8 n_Q2_34,n_E_34,n_cos_34

      !! Variables of the e+e- pair
      real*8 cos_3,phi_3,w_phi_3

      !! random numbers
      external getrnd
      external getphi
      integer i,j
      real*8 csi(40),xi(2),xi1(1),eta(4)

      !! Variables for born photons
      real*8 cos_g1,cos_g2,sin_g2,cos_2_34,E_g2
      real*8 w_cos_g2,n_cos_g2,phi_g2,w_phi_g2
      real*8 w_cos_g1,n_cos_g1,phi_g1,w_phi_g1

      !! Rotation
      external rotation_xy
      real*8 jacobian

      real*8 var
      real*8 eps
      common/epssoft/eps

      external dLips_2

      !! CUTS
      real*8 cos_el,cos_ep

      real*16 b1q,b2q,Qmaxq,Qminq

      !! MULTICHANNEL SAMPLING
      real*8 n_smpl
      real*8 w_1,w_2,Q2_13,Q2_24
      real*8 thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      character*2 sampling,angdip
      real*8 cos_12,sin_g1
      real*8 E_2,w_E_2,n_E_2

      !! Parameters of the resonance
      real*8 M,G
      s = dot(p1+p2,p1+p2)

      Q_34(0) = 0d0
      Q_34(1) = 0d0
      Q_34(2) = 0d0
      Q_34(3) = 0d0

      w=1d0
      w_flat=1d0
      w_smpl=1d0
      n_smpl=1d0
      norm=1d0
      phsp=1d0
      !!******************************
      !! RF = Center of Mass  !!
      !!******************************
      !! Invariant mass of the virtual photon sampled as 1/Q_34^2
      if(igen.eq.0) then
         Q2_23=dot(qph(ng,:)+p3,qph(ng,:)+p3)

         Q_34=qph(ng,:)+p3+p4
         Q23=qph(ng,:)+p3

!        cos234=Q234(3)/sqrt(abs(tridot(Q234,Q234)))

         Q2_34=dot(p3+p4+qph(ng,:),p3+p4+qph(ng,:))

         do i=1,ng
            En(i)=qph(i,0)
         enddo

         beta_e = p1(3)/p1(0)
         do i=1,ng
            qcos(i)=qph(i,3)/qph(i,0)
         enddo
      endif

      if(igen.eq.1) then
         beta_e = p1(3)/p1(0)
      endif

      ngi=ng-1
      !! Parameters of the rho
      M = 0.77456d0
      G = 0.14832d0

      !! I generate the invariant mass of the p3+p4 pair
      Qmin =(m3+m4)**2
      Qmax = s*(1.d0-eps)

      if(sampling.eq.'BW') then
         call sampling_BW_resonances(M,G,Qmin,Qmax,Q2_34,w_Q2_34,
     >     n_Q2_34,igen)

      else if(sampling.eq.'CO') then
         call sampling_energy_mxx(Qmin,Qmax,Q2_34,w_Q2_34,n_Q2_34,igen)

      else if(sampling.eq.'SC') then
         call sampling_prop(s,Qmin,Qmax,Q2_34,w_Q2_34,n_Q2_34,igen)

      else if(sampling.eq.'FL') then
         if(igen.eq.1) then
            call getrnd(xi,1)
            Q2_34 = Qmin + xi(1)*(Qmax-Qmin)
         endif

         w_Q2_34=1.d0
         n_Q2_34=1.d0

         w_flat = w_flat*(Qmax-Qmin)
      endif

      w_smpl = w_smpl*w_Q2_34
      n_smpl = n_smpl*n_Q2_34

!! SECOND DECAY
!! Q234 -> p4 + Q23
!! I sample the invariant mass of the second fermion-gamma pair as 1/Q^2
c$$$      Q23min = m3**2
c$$$      Q23max =(Q2_234)
c$$$      call sampling_energy_mxx(Q23min,Q23max,Q2_23,w_Q2_23,n_Q2_23,igen)

      Q23min=m3**2+eps*s/2.d0*(1-sqrt(1.d0-(2.d0*m3)**2/s))

      Q23max =(Q2_34)
      call sampling_energy_mxx_sc(m3,Q23min,Q23max,Q2_23,w_Q2_23,n_Q2_23
     >     ,igen)

      m23=sqrt(Q2_23)


      w_smpl=w_smpl*w_Q2_23
      n_smpl=n_smpl*n_Q2_23

!! I generate the energy of the photon 2 as 1/E_2
      do i=1,ngi-1

         Emin=eps*sqrt(s)/2d0
         Emax=sqrt(s)/2d0

c$$$         call sampling_energy(Emin,Emax,En(i),w_E_2,n_E_2,igen)
c$$$         if(abs(p3(0)/0.49511841186217359-1.d0).le.1.d-8) then
c$$$            print*,'peso en',w_E_2
c$$$         endif
c$$$         if(abs(p4(0)/0.49511841186217359-1.d0).le.1.d-8) then
c$$$            print*,'peso en',w_E_2
c$$$         endif

         w_smpl = w_smpl*w_E_2
         n_smpl = n_smpl*n_E_2

      enddo
      !! I generate cos_1 and cos_2 of the two photons as collinear
!! Sampling with (1-b^2 cos^2)

      cmin_g = -1d0
      cmax_g = 1d0
      do i=1,ngi-1

         call sampling_cos_collsq(beta_e,1d0,cmin_g,cmax_g,
     .        qcos(i),w_cos_g1,n_cos_g1,igen)

c$$$         if(abs(p3(0)/0.49511841186217359-1.d0).le.1.d-8) then
c$$$            print*,'peso ang',w_cos_g1
c$$$         endif
c$$$         if(abs(p4(0)/0.49511841186217359-1.d0).le.1.d-8) then
c$$$            print*,'peso ang',w_cos_g1
c$$$         endif


         w_smpl = w_smpl*w_cos_g1
         n_smpl = n_smpl*n_cos_g1

      enddo
!!!until now I have generated 4 sampled variables
      w=w*w_smpl
      norm=n_smpl

*************LAST PHOTON

      if(angdip.eq.'co') then
         call sampling_cos_collsq(beta_e,1d0,-1d0,1d0,
     .        qcos(ngi),w_cos_g1,n_cos_g1,igen)

      else if(angdip.eq.'fl') then
         if(igen.eq.1) then
            call getrnd(xi1,1)
            qcos(ngi) = -1.d0 + 2.d0*xi1(1)
         endif
         w_cos_g1 = 1.d0
         n_cos_g1 = 1.d0
         w_flat  = w_flat *  2.d0

      endif


c$$$      if(abs(p3(0)/0.49511841186217359-1.d0).le.1.d-8) then
c$$$         print*,'peso 23',w_q2_23
c$$$         print*,'peso n-1 phot',w_cos_g1
c$$$         print*,'peso q2_34',w_q2_34,sampling
c$$$      endif
c$$$      if(abs(p4(0)/0.49511841186217359-1.d0).le.1.d-8) then
c$$$         print*,'peso 23',w_q2_23
c$$$         print*,'peso n-1 phot',w_cos_g1
c$$$         print*,'peso q2_34',w_q2_34,sampling
c$$$      endif


      w = w*w_cos_g1
      norm = norm*n_cos_g1

      if(igen.eq.1) then

         call getrnd(csi,40)

         do i=1,ngi
!! Generating the azimuths of 1 and 2 flat
            call getphi(csi(i),phi_g1,w_phi_g1)
            w_flat = w_flat*w_phi_g1
            qphi(i)= phi_g1
         enddo

!! At this point I have the cosines and all the angles. I can calculate the momenta of the first n-1 photon
         do i=1,ngi-1

            sini=sqrt(1d0-qcos(i)**2)

            qph(i,0)=En(i)
            qph(i,1)=En(i)*sini*cos(qphi(i))
            qph(i,2)=En(i)*sini*sin(qphi(i))
            qph(i,3)=En(i)*qcos(i)

         enddo


         do i=1,ngi
            do j=1,ngi
              sini=dsqrt(1.d0-qcos(i)**2)
              sinj=dsqrt(1.d0-qcos(j)**2)
              cosij(i,j)= qcos(i)*qcos(j)
              cosij(i,j)= cosij(i,j)+sini*cos(qphi(i))*sinj*cos(qphi(j))
              cosij(i,j)= cosij(i,j)+sini*sin(qphi(i))*sinj*sin(qphi(j))
            enddo
         enddo

******QUI ORA RICAVO K_N E DI CONSEGUENZA Q_34
         sumen=0.d0
         do i=1,ngi-1
            sumen=sumen+qph(i,0)
         enddo

         sumij=0.d0
         do i=1,ngi-1
            do j=1,ngi-1
               sumij=sumij+qph(i,0)*qph(j,0)*(1.d0-cosij(i,j))
            enddo
         enddo

         sumin=0.d0
         do i=1,ngi-1
            sumin=sumin+qph(i,0)*(1.d0-cosij(i,ngi))
         enddo

         qph(ngi,0)=s-Q2_34-2.d0*sqrt(s)*sumen+sumij
         qph(ngi,0)=qph(ngi,0)/(2.d0*(sqrt(s)-sumin))

         if(qph(ngi,0).le.0.d0) then
            icut=1
            phsp=0.d0
            w=0.d0
            return
         endif

         sin_ng=dsqrt(1.d0-qcos(ngi)**2)
         qph(ngi,1)=qph(ngi,0)*sin_ng*cos(qphi(ngi))
         qph(ngi,2)=qph(ngi,0)*sin_ng*sin(qphi(ngi))
         qph(ngi,3)=qph(ngi,0)*qcos(ngi)

         Q_34(0)=sqrt(s)
         Q_34(1)=0.d0
         Q_34(2)=0.d0
         Q_34(3)=0.d0

         do i=1,ngi
            Q_34(0)=Q_34(0)-qph(i,0)
            do j=1,3
               Q_34(j)=Q_34(j)-qph(i,j)
            enddo
         enddo
         if(Q_34(0).le.0.d0) then
            icut=1
            phsp=0.d0
            w=0.d0
            return
         endif

***    AT THIS POINT I  HAVE GENERATED N-1 PHOTONS, CHECKED MOM CONSERVATION
***    NOW I HAVE TO DECAY THE Q_34 IN P3,P4, Q(NG)

         call getrnd(eta,4)
         cos4 = -1.d0 + 2.d0*eta(1)
         w_flat = w_flat*2.d0

         call getphi(eta(2),phi4,w_phi4)
         w_flat = w_flat*w_phi4

         call dLips_2(Q_34,cos4,phi4,Q23,m23,p4,m4,phsp2_24_4,icut)

         cos2 = -1.d0 +	2.d0*eta(3)
         w_flat = w_flat * 2.d0

         call getphi(eta(4),phi2,w_phi2)
         w_flat = w_flat*w_phi2

         call dLips_2(Q23,cos2,phi2,qph(ng,:),0d0,p3,m3,phsp2_2_3,icut)


         if(abs(p3(0)/0.49511841186217359-1.d0).le.1.d-8) then
            print*,'IIIIIIIN, peso 23',w_q2_23,sampling,angdip
         endif
         if(abs(p4(0)/0.49511841186217359-1.d0).le.1.d-8) then
            print*,'IIIIIIIN, peso 23',w_q2_23,sampling,angdip
         endif

!        ! 1/(2 (2 pi)^3)^3 coming from the integrals in dq1 dq2 dQ34
         phsp = (1d0/(2d0*(2d0*pi)**3))**(float(ng))

         !! 2->2 phase space
         phsp = phsp*phsp2_24_4*phsp2_2_3

         !! Jacobian of the delta^0 * change of variable form |Q_34| to E_34
         suminc = 0.d0
         do i=1,ngi-1
            suminc=suminc+qph(i,0)*cosij(i,ngi)
         enddo

         jacobian = 1.d0
         do i=1,ngi
            jacobian=jacobian*qph(i,0)
         enddo

         jacobian = jacobian/(Q_34(0)+qph(ngi,0)+suminc)

         phsp = phsp*jacobian

         !! (2pi)^4 coming from the overall delta
         !! (2pi)^3 coming from the Q^2_34 integral
         phsp = phsp*(2d0*pi)**4*(2d0*pi)**6
         phsp = phsp*w_flat!*w_smpl*n_smpl
      endif

      qval=0.d0
      call prop_value(2,p1,p2,p3,p4,qph,ng,qval)
      cond=.false.

      do i=1,ng-1
         if(qval(ng,1).le.qval(i,1)) then
            cond=.true.
         endif
      enddo

      if(cond) then
         phsp=0.d0
         icut=1
         w=0.d0
      else
         phsp=phsp*float(ng)
      endif

c$$$      cond=.false.
c$$$
c$$$
c$$$      do i=1,ng-1
c$$$         if(qval(ng,2).ge.qval(i,2)) then
c$$$            cond=.true.
c$$$         endif
c$$$      enddo
c$$$
c$$$      if(cond) then
c$$$         phsp=0.d0
c$$$         icut=1
c$$$         w=0.d0
c$$$      else
c$$$         phsp=phsp*float(ng)
c$$$      endif


c$$$      call prop_value(1,p1,p2,p3,p4,qph,ng,qval)
c$$$      cond=.false.
c$$$
c$$$      do i=1,ng-1
c$$$         if(qval(ng,1).ge.qval(i,1)) then
c$$$            cond=.true.
c$$$         endif
c$$$c$$$      enddo
c$$$
c$$$      if(cond) then
c$$$         phsp=0.d0
c$$$         icut=1
c$$$         w=0.d0
c$$$      else
c$$$         phsp=phsp*float(ng)
c$$$      endif


      cond=.false.
      do i=1,ng-2
         if(qph(ng-1,0).le.qph(i,0)) then
            cond=.true.
         endif
      enddo

      if(cond) then
         phsp=0.d0
         icut=1
         w=0.d0
      else
         phsp=phsp*float(ng-1)
      endif



****************************
!      w=w_Q2_34
****************************
! norm=norm*w_flat
      return
      end


**************

**************MAURO SAMPLING
      subroutine CHN_2g_multichannel_angular(p1,p2,p3,m3,p4,m4,
     . q1,q2,w,norm,phsp,icut,sampling,coll,igen)
*******************************************************************
      implicit none
      dimension p1(0:3),p2(0:3)
      real*8 qtmp(0:3),pz(0:3)
      dimension p3(0:3),p4(0:3),Q_34(0:3),q1(0:3),q2(0:3)
      real*8 p1z(0:3),p2z(0:3),p3z(0:3),p4z(0:3)
      real*8 Q_34
      integer icut
      real*8 p1,p2,p3,p4,q1,q2,w,norm
      real*8 Qmin,Qmax,emax
      real*8 cmin_34,cmax_34,cmin_g,cmax_g
      real*8 cmax,cmin
      real*8 qph(40,0:3),xi(1),xi1(3),qval(40,4)
      common/coseni/cmax,cmin
      real*8 ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
      real*8 s
      real*8 tridot,dot
      external tridot,dot
      real*8 mod234,phsp,phsp_2
      real*8 m3,m4
      logical cond
      integer igen

      !! Weights: flat and sampled
      real*8 w_flat, w_smpl

      !! Variables of Q34
      real*8 Q2_34,E_34,cos_34,sin_34,phi_34,modQ34
      real*8 beta_e,beta_34,Q34min,Q34max,Q34
      real*8 beta,q1prop,q2prop
      real*8 phi34,w_phi34,cos34
      !! Weights of Q34
      real*8 w_Q2_34,w_E_34,w_cos_34,w_phi_34

      !! Norms of Q34
      real*8 n_Q2_34,n_E_34,n_cos_34

      !! Variables of the e+e- pair
      real*8 cos_3,phi_3,w_phi_3

      real*8 phi_1,phi_2,E_1
      !! random numbers
      external getrnd
      external getphi
      integer i,coll,counter
      common/countercos/counter
      real*8 csi(5)

      !! Variables for born photons
      real*8 cos_g1,cos_g2,sin_g2,cos_2_34,E_g2
      real*8 w_cos_g2,n_cos_g2,phi_g2,w_phi_g2
      real*8 w_cos_g1,n_cos_g1,phi_g1,w_phi_g1

      !! Rotation
      external rotation_xy
      real*8 jacobian

      real*8 var
      real*8 eps
      common/epssoft/eps

      real*8 min1,min2
      external dLips_2

      !! CUTS
      real*8 cos_el,cos_ep

      real*16 b1q,b2q,Qmaxq,Qminq

      !! MULTICHANNEL SAMPLING
      real*8 n_smpl
      real*8 w_1,w_2,Q2_13,Q2_24
      real*8 thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      character*2 sampling,angdip
      real*8 cos_12,sin_g1
      real*8 E_2,w_E_2,n_E_2

      !! Parameters of the resonance
      real*8 M,G
      s = dot(p1+p2,p1+p2)

      Q_34(0) = 0d0
      Q_34(1) = 0d0
      Q_34(2) = 0d0
      Q_34(3) = 0d0
      w=1d0
      w_flat=1d0
      w_smpl=1d0
      n_smpl=1d0
      norm=1d0
      phsp=1d0
      !!******************************
      !! RF = Center of Mass  !!
      !!******************************
      !! Invariant mass of the virtual photon sampled as 1/Q_34^2

      if(igen.eq.0) then

         if(coll.eq.1) then
            pz=p1
         elseif(coll.eq.2) then
            pz=p3
         elseif(coll.eq.3) then
            pz=p4
         endif

         Q_34=p3+p4
         Q2_34=dot(Q_34,Q_34)
         Q34=dsqrt(tridot(Q_34,Q_34))

         cos_g2= tridot(q2,pz)/dsqrt(tridot(pz,pz))/q2(0)

         phi_g2 = atan2(q2(2),q2(1))
         phi_g1 = atan2(q1(2),q1(1))

         E_2 = q2(0)

      endif

      !! Parameters of the rho
      M = 0.77456d0
      G = 0.14832d0

!! I generate the invariant mass of the p3+p4 pair
      Qmin =(m3+m4)**2
      Qmax = s-2d0*sqrt(s)*egmin

      if(sampling.eq.'BW') then
         call sampling_BW_resonances(M,G,Qmin,Qmax,Q2_34,w_Q2_34,
     >     n_Q2_34,igen)
      else if(sampling.eq.'CO') then
         call sampling_energy_mxx(Qmin,Qmax,Q2_34,w_Q2_34,n_Q2_34,igen)
      else if(sampling.eq.'SC') then
         call sampling_prop(s,Qmin,Qmax,Q2_34,w_Q2_34,n_Q2_34,igen)
      else if(sampling.eq.'FL') then
         if(igen.eq.1) then
            call getrnd(xi,1)
            Q2_34 = Qmin + xi(1)*(Qmax-Qmin)
         endif
         w_Q2_34 = 1.d0
         n_Q2_34 = 1.d0
         w_flat =w_flat*(Qmax-Qmin)
      endif

      w_smpl = w_smpl*w_Q2_34
      n_smpl = n_smpl*n_Q2_34

!!!!!!i generate the Q34 spatial part
      Q34min = 0.d0
      Q34max = dsqrt(s-Q2_34)


      beta=p1(3)/p1(0)
      call sampling_cos_collsq(beta,1d0,-1.d0,1.d0,
     .        cos34,w_cos_34,n_cos_34,igen)

      w_smpl=w_smpl*w_cos_34
      n_smpl=n_smpl*n_cos_34

      if(igen.eq.1) then

         call getrnd(xi1,3)
         Q34 = Q34min + xi1(1)*(Q34max-Q34min)
         w_flat = w_flat*(Q34max-Q34min)

c$$$         cos34 = -1.d0+xi1(2)*2.d0
c$$$         w_flat = w_flat*2.d0

         call getphi(xi1(3),phi_34,w_phi34)
         w_flat = w_flat*w_phi34

         Q_34(0) = dsqrt(Q2_34+Q34**2)
         Q_34(1) = Q34*dsqrt(1.d0-cos34**2)*cos(phi_34)
         Q_34(2) = Q34*dsqrt(1.d0-cos34**2)*sin(phi_34)
         Q_34(3) = Q34*cos34

         call getrnd(csi,4)

         call getphi(csi(3),phi_3,w_phi_3)
         w_flat=w_flat*w_phi_3

         cos_3= -1d0 + 2d0*csi(4)
         w_flat=w_flat*2d0

         call dLips_2(Q_34,cos_3,phi_3,p3,m3,p4,m4,phsp_2,icut)

      endif


******MI METTO NEL SISTEMA DI RIFERIMENTO DOVE P_J CON P_J= GAMBA SU CUI CAMPIONO È ORIENTATA LUNGO Z

      if(coll.eq.1) then
         beta=p1(3)/p1(0)
         pz=p1
      elseif(coll.eq.2) then
         beta=dsqrt(1.d0-m3**2/p3(0)**2)
         pz=p3
      elseif(coll.eq.3) then
         beta=dsqrt(1.d0-m4**2/p4(0)**2)

         pz=p4
      endif

      call rot(1,pz,p1,p1)
      call rot(1,pz,p2,p2)
      call rot(1,pz,p3,p3)
      call rot(1,pz,p4,p4)
      call rot(1,pz,Q_34,Q_34)


      if (coll.eq.1) then
         call sampling_cos_collsq(beta,1d0,-1.d0,1.d0,
     .        cos_g2,w_cos_g2,n_cos_g2,igen)
      else
         call sampling_cos_coll(beta,1d0,-1.d0,1.d0,
     .        cos_g2,w_cos_g2,n_cos_g2,igen)
      endif

!      cos_g2=-1.d0+csi(2)*2.d0


      w_smpl = w_smpl*w_cos_g2
      n_smpl = n_smpl*n_cos_g2

!! I generate the energy of the photon 2 as 1/E_2

      Emin=eps*sqrt(s)/2d0
      Emax=sqrt(s)/2.d0


      call sampling_energy(Emin,Emax,E_2,w_E_2,n_E_2,igen)
      w_smpl = w_smpl*w_E_2
      n_smpl = n_smpl*n_E_2

!      E_2 =Emin+ csi(1)*(Emax-Emin)

      w=w*w_smpl
      norm=n_smpl

      if(igen.eq.1) then

**********Set E_1
         E_1 = p1(0)+p2(0) -Q_34(0) -E_2

         if(E_1.le.0.d0) then
            phsp=0.d0
            icut=1
            w=0.d0
            norm=0.d0
            return
         endif
*********Set theta_1
         cos_g1 = (-Q_34(3) - E_2*cos_g2)/E_1

         if(abs(cos_g1).ge.1.d0) then
            phsp=0.d0
            icut=1
            w=0.d0
            norm=0.d0
            return
         endif

*********Set phi1,phi2

         sin_g1=dsqrt(1.d0-cos_g1**2)
         sin_g2=dsqrt(1.d0-cos_g2**2)

         call solvesystem(E_1,E_2,sin_g1,sin_g2,Q_34(1),Q_34(2),
     >        phi_1,phi_2)

         if(isnan(phi_1)) then
            phsp=0.d0
            icut=1
            w=0.d0
            norm=0.d0
            return
         else
            w_flat=w_flat*2.d0
         endif

         q1(0)=E_1
         q1(1)=E_1*sin_g1*cos(phi_1)
         q1(2)=E_1*sin_g1*sin(phi_1)
         q1(3)=E_1*cos_g1

         q2(0)=E_2
         q2(1)=E_2*sin_g2*cos(phi_2)
         q2(2)=E_2*sin_g2*sin(phi_2)
         q2(3)=E_2*cos_g2


!        ! 1/(2 (2 pi)^3)^3 coming from the integrals in dq1 dq2 dQ34
         phsp = (1d0/(2d0*(2d0*pi)**3)**3)

         !! 2->2 phase space
         phsp = phsp*phsp_2

!! Jacobian
         jacobian = Q34**2/Q_34(0)*E_1*E_2
!   Jacobian of delta_z
         jacobian = jacobian/E_1

         jacobian = jacobian/abs(E_1*E_2*sin_g1*sin_g2)
         jacobian = jacobian/(abs(sin(phi_1)*cos(phi_2)
     >        -sin(phi_2)*cos(phi_1)))

         if(q1(0).le.0.d0) then
            icut=1
            phsp=0.d0
            w=0.d0
            norm=0.d0
            return
         endif
         phsp = phsp*jacobian

         !! (2pi)^4 coming from the overall delta
         !! (2pi)^3 coming from the Q^2_34 integral
         phsp = phsp*(2d0*pi)**4*(2d0*pi)**3
         phsp = phsp*w_flat!*w_smpl*n_smpl


         call rot(-1,pz,q1,q1)
         call rot(-1,pz,q2,q2)
      endif

      call rot(-1,pz,p1,p1)
      call rot(-1,pz,p2,p2)
      call rot(-1,pz,p3,p3)
      call rot(-1,pz,p4,p4)
c$$$
c$$$      if(igen.eq.1) then
c$$$         if(coll.eq.1) then
c$$$            continue
c$$$         elseif(coll.eq.2.and.igen.eq.1) then
c$$$            pz=p3
c$$$            print*,p3
c$$$         elseif(coll.eq.3.and.igen.eq.1) then
c$$$            pz=p4
c$$$            print*,p4
c$$$         endif
c$$$         print*,'weight',(1.d0-beta*cos_g2)
c$$$         print*,'weight ex',w_cos_g2
c$$$         print*,coll
c$$$         print*,cos_g2
c$$$         print*,tridot(pz,q2)/dsqrt(tridot(pz,pz))/q2(0)
c$$$         print*,' '
c$$$      endif

      if(igen.eq.1) then

         qph(1,:)=q1
         qph(2,:)=q2

         call prop_value(1,p1,p2,p3,p4,qph,2,qval)
         q1prop=0.d0
         q2prop=0.d0
         do i=1,4
            q1prop=q1prop+qval(1,i)
            q2prop=q2prop+qval(2,i)
         enddo


         if(q1prop.le.q2prop) then
            phsp=0.d0
            icut=1
            w=0.d0
            return
         else
            phsp=phsp*2.d0
         endif
      endif


      return
      end

********************
