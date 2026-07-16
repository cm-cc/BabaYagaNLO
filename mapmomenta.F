      subroutine closest(ng,qph,p1,p2,p3,p4,q1,q2,q3,q4,iclose)
      implicit double precision (a-h,o-z)
      parameter (maxph = 40)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension q1(0:3),q2(0:3),q3(0:3),q4(0:3)
      character*2 fs
      common/finalstate/fs
      dimension qph(maxph,0:3)
      dimension iclose(maxph)
      integer emission(maxph)
      common/radpattern/nph(4)

      logical cond
      common/radoption/ifisrfsr
      integer ifisrfsr
      iclose = 0

      q1 = 0d0
      q2 = 0d0
      q3 = 0d0
      q4 = 0d0

      do k = 1,ng
         q = qph(k,:)


         if(q(0).le.p1(0)) then
            en=q(0)**2.d0
         else
            en=p1(0)**2.d0
         endif
         c1=en*(1-tridot(p1,q)/dsqrt(tridot(p1,p1))/q(0))

         if(q(0).le.p2(0)) then
            en=q(0)**2.d0
         else
            en=p2(0)**2.d0
         endif
         c2=en*(1-tridot(p2,q)/dsqrt(tridot(p2,p2))/q(0))

         if(q(0).le.p3(0)) then
            en=q(0)**2.d0
         else
            en=p3(0)**2.d0
         endif
         c3=en*(1-tridot(p3,q)/dsqrt(tridot(p3,p3))/q(0))

         if(q(0).le.p4(0)) then
            en=q(0)**2.d0
         else
            en=p4(0)**2.d0
         endif
         c4=en*(1-tridot(p4,q)/dsqrt(tridot(p4,p4))/q(0))

         cond=(fs.eq.'mr'.or.fs.eq.'pr')
         if (fs.eq.'ee'.or.fs.eq.'mm'.or.fs.eq.'pp'
     .       .or.fs.eq.'mr'.or.fs.eq.'er'.or.fs.eq.'pr') then
            if (c1.lt.c2.and.c1.lt.c3.and.c1.lt.c4) iclose(k) = 1
            if (c2.lt.c1.and.c2.lt.c3.and.c2.lt.c4) iclose(k) = 2
            if (c3.lt.c2.and.c3.lt.c1.and.c3.lt.c4) iclose(k) = 3
            if (c4.lt.c1.and.c4.lt.c2.and.c4.lt.c3) iclose(k) = 4
         elseif (fs.eq.'gg'.or.fs.eq.'pi') then
            if (c1.lt.c2) iclose(k) = 1
            if (c2.lt.c1) iclose(k) = 2
         endif

         if (cond.and.ifisrfsr.eq.10) then
            if (c1.lt.c2) iclose(k) = 1
            if (c2.lt.c1) iclose(k) = 2
         elseif (cond.and.ifisrfsr.eq.01) then
            if (c3.lt.c4) iclose(k) = 3
            if (c4.lt.c3) iclose(k) = 4
         endif
!!!
ccc   iclose(k) = emission(k)
!!!   ic = iclose(k)
         if (ic.eq.1) then
            q1 = q1 + q
         endif
         if (ic.eq.2) then
            q2 = q2 + q
         endif
         if (ic.eq.3) then
            q3 = q3 + q
         endif
         if (ic.eq.4) then
            q4 = q4 + q
         endif
      enddo

      return
      end
********************************************************************

      subroutine closestOLD(ng,qph,p1,p2,p3,p4,q1,q2,q3,q4,iclose)
      implicit double precision (a-h,o-z)
      parameter (maxph = 40)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension q1(0:3),q2(0:3),q3(0:3),q4(0:3)
      character*2 fs
      common/finalstate/fs
      dimension qph(maxph,0:3)
      dimension iclose(maxph)
      integer emission(maxph)
      common/radpattern/nph(4)

      logical cond
      common/radoption/ifisrfsr
      integer ifisrfsr
      iclose = 0

      q1 = 0d0
      q2 = 0d0
      q3 = 0d0
      q4 = 0d0

c      p1modmu = 1.d0/sqrt(tridot(p1,p1))
c      p2modmu = 1.d0/sqrt(tridot(p2,p2))
c      p3modmu = 1.d0/sqrt(tridot(p3,p3))
c      p4modmu = 1.d0/sqrt(tridot(p4,p4))

      do k = 1,ng
c---
c         if (k.le.nph(1)) emission(k) = 1
c         if (k.gt.nph(1).and.k.le.(nph(1)+nph(2))) emission(k) = 2
c         if (k.gt.(nph(1)+nph(2)).and.k.le.(nph(1)+nph(2)+nph(3)))
c     .        emission(k) = 3
c         if (k.gt.(nph(1)+nph(2)+nph(3))) emission(k) = 4
c---
         q = qph(k,:)

c         c1 = tridot(p1,q)*p1modmu/q(0)
c         c2 = tridot(p2,q)*p2modmu/q(0)
c         c3 = tridot(p3,q)*p3modmu/q(0)
c         c4 = tridot(p4,q)*p4modmu/q(0)

         c1 = dot(p1,q)
         c2 = dot(p2,q)
         c3 = dot(p3,q)
         c4 = dot(p4,q)


         cond=(fs.eq.'mr'.or.fs.eq.'pr')
         if (fs.eq.'ee'.or.fs.eq.'mm'.or.fs.eq.'pp'
     .       .or.fs.eq.'mr'.or.fs.eq.'er'.or.fs.eq.'pr') then
            if (c1.lt.c2.and.c1.lt.c3.and.c1.lt.c4) iclose(k) = 1
            if (c2.lt.c1.and.c2.lt.c3.and.c2.lt.c4) iclose(k) = 2
            if (c3.lt.c2.and.c3.lt.c1.and.c3.lt.c4) iclose(k) = 3
            if (c4.lt.c1.and.c4.lt.c2.and.c4.lt.c3) iclose(k) = 4
         elseif (fs.eq.'gg'.or.fs.eq.'pi') then
            if (c1.lt.c2) iclose(k) = 1
            if (c2.lt.c1) iclose(k) = 2
         endif

         if (cond.and.ifisrfsr.eq.10) then
            if (c1.lt.c2) iclose(k) = 1
            if (c2.lt.c1) iclose(k) = 2
         elseif (cond.and.ifisrfsr.eq.01) then
            if (c3.lt.c4) iclose(k) = 3
            if (c4.lt.c3) iclose(k) = 4
         endif

!!!
ccc         iclose(k) = emission(k)
!!!
         ic = iclose(k)
         if (ic.eq.1) then
            q1 = q1 + q
         endif
         if (ic.eq.2) then
            q2 = q2 + q
         endif
         if (ic.eq.3) then
            q3 = q3 + q
         endif
         if (ic.eq.4) then
            q4 = q4 + q
         endif
      enddo

      return
      end
*
      subroutine closest_BCK(ng,qph,p1,p2,p3,p4,q1,q2,q3,q4,iclose)
      implicit double precision (a-h,o-z)
      parameter (maxph = 40)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension q1(0:3),q2(0:3),q3(0:3),q4(0:3)
      character*2 fs
      common/finalstate/fs
      dimension qph(maxph,0:3)
      dimension iclose(maxph)
      integer emission(maxph)
      common/radpattern/nph(4)

      do k = 1,maxph
         iclose(k) = 0
      enddo

      do k = 0,3
         q1(k) = 0.d0
         q2(k) = 0.d0
         q3(k) = 0.d0
         q4(k) = 0.d0
      enddo

      p1modmu = 1.d0/sqrt(tridot(p1,p1))
      p2modmu = 1.d0/sqrt(tridot(p2,p2))
      p3modmu = 1.d0/sqrt(tridot(p3,p3))
      p4modmu = 1.d0/sqrt(tridot(p4,p4))

      do k = 1,ng
c---
         if (k.le.nph(1)) emission(k) = 1
         if (k.gt.nph(1).and.k.le.(nph(1)+nph(2))) emission(k) = 2
         if (k.gt.(nph(1)+nph(2)).and.k.le.(nph(1)+nph(2)+nph(3)))
     .        emission(k) = 3
         if (k.gt.(nph(1)+nph(2)+nph(3))) emission(k) = 4
c---
         do i = 0,3
            q(i) = qph(k,i)
         enddo
         c1 = tridot(p1,q)*p1modmu/q(0)
         c2 = tridot(p2,q)*p2modmu/q(0)
         c3 = tridot(p3,q)*p3modmu/q(0)
         c4 = tridot(p4,q)*p4modmu/q(0)

         if (fs.eq.'ee'.or.fs.eq.'mm'.or.fs.eq.'pp'
     .       .or.fs.eq.'mr'.or.fs.eq.'er') then
            if (c1.gt.c2.and.c1.gt.c3.and.c1.gt.c4) iclose(k) = 1
            if (c2.gt.c1.and.c2.gt.c3.and.c2.gt.c4) iclose(k) = 2
            if (c3.gt.c2.and.c3.gt.c1.and.c3.gt.c4) iclose(k) = 3
            if (c4.gt.c1.and.c4.gt.c2.and.c4.gt.c3) iclose(k) = 4
         elseif (fs.eq.'gg'.or.fs.eq.'pi') then
            if (c1.gt.c2) iclose(k) = 1
            if (c2.gt.c1) iclose(k) = 2
         endif
!!!
ccc         iclose(k) = emission(k)
!!!
         ic = iclose(k)
         if (ic.eq.1) then
            do j = 0,3
               q1(j) = q1(j) + q(j)
            enddo
         endif
         if (ic.eq.2) then
            do j = 0,3
               q2(j) = q2(j) + q(j)
            enddo
         endif
         if (ic.eq.3) then
            do j = 0,3
               q3(j) = q3(j) + q(j)
            enddo
         endif
         if (ic.eq.4) then
            do j = 0,3
               q4(j) = q4(j) + q(j)
            enddo
         endif
      enddo
      return
      end


*************************
     	subroutine cluster_ang(nchrg,pchrg,nphot,pphot,nsftcll,idxclst,
     .                        idphtme)
c*********************************************************************
      implicit none
c
c     clustered according to angle
c
c            momentum entries p==(px,py,px,p0)
c     INPUT: NCHRG, number of charged particles
c            PCHRG, momentum of charged particles
c            NPHOT, number of final state photons
c            PPHOT, momentum of final state photons
c            NSFTCLL, number of "soft-collinear" photons to be clustered,
c                    the number of matrix element photon will be NPHOT-NSFTCLL
c                 notice that NCHRG,NPHOT,NSFTCLL enter into matirces dimension
c                 should be PARAMETERS in the calling routines
c     OUTPUT: IDXCLST(21,nchrg), IDXCLST(1,J)== # of photons clustered
c                                with particle J,  IDXCLST(L,J)== label of
c                                the photons clustered with the J-th particle
c             IDPHTME(NPHOT-NSFTCLL), labels of unclustered photons
c
      integer nchrg,nphot,nsftcll
      integer idxclst(21,nchrg),idphtme(nphot-nsftcll)
      real*8 pchrg(4,nchrg),pphot(4,nphot)
c
      integer nme,j1,j2,idsftcll(nsftcll),j3,id
      real*8 pclst(4,nchrg),enmepht(nphot-nsftcll),cs,cstmp
      real*8 ppchrg(nchrg)
      logical lg
c
      if(nsftcll.gt.20) then
         write(*,*)'too many soft photons, update array dimension'
         stop
      endif
c
c     finding the hard photons
c
      nme=nphot-nsftcll
      do j1=1,nme
         enmepht(j1)=0.d0
      enddo
c
c     DBG
c      write(*,*)'photon momenta'
c      do j1=1,nphot
c         write(*,*)(pphot(j2,j1),j2=1,4)
c      enddo
c end DBG
c
      do j1=1,nphot
         lg=.true.
         do j2=1,nme
            if( (pphot(4,j1).gt.enmepht(j2)).and.lg) then
               do j3=nme,j2+1,-1
                  enmepht(j3)=enmepht(j3-1)
                  idphtme(j3)=idphtme(j3-1)
               enddo
               enmepht(j2)=pphot(4,j1)
               idphtme(j2)=j1
               lg=.false.
            endif
         enddo
      enddo
      j3=0
      do j1=1,nphot
         lg=.true.
         do j2=1,nme
            if(idphtme(j2).eq.j1) lg=.false.
         enddo
         if(lg) then
            j3=j3+1
            idsftcll(j3)=j1
         endif
      enddo
c     DBG
c      write(*,*)'enmepht',(enmepht(j1),j1=1,nme)
c      write(*,*)'idphtme',(idphtme(j1),j1=1,nme)
c      write(*,*)'idsftcll',(idsftcll(j1),j1=1,nsftcll)
c end DBG
c
c     clustering the soft photons
c
      do j1=1,nchrg
         ppchrg(j1)=sqrt(pchrg(1,j1)**2+pchrg(2,j1)**2+pchrg(3,j1)**2)
         idxclst(1,j1)=0
      enddo
      do j1=1,nsftcll
         id=idsftcll(j1)
c     DBG
c         write(*,*)'id',id
c     end DBG
         cs=0.d0
         do j2=1,nchrg
            cstmp=( pchrg(1,j2)*pphot(1,id)+pchrg(2,j2)*pphot(2,id)+
     .           pchrg(3,j2)*pphot(3,id) )/ppchrg(j2)
c     DBG
c            write(*,*)'cstmp',cstmp
c     end DBG
            if(cstmp.gt.cs) then
               j3=j2
               cs=cstmp
            endif
         enddo
c     DBG
c         write(*,*)'j3',j3
c     end DBG
         idxclst(1,j3)= idxclst(1,j3)+1
         idxclst(idxclst(1,j3)+1,j3)=id
c     DBG
c         write(*,*)'j3',j3
c         write(*,*)'idx',idxclst(1,j3)
c     end DBG
      enddo
c
c     DBG
c      write(*,*)'pchrg'
c      do j1=1,nchrg
c         write(*,*)(pchrg(j2,j1),j2=1,4)
c      enddo
c      write(*,*)'ppchrg',(ppchrg(j1),j1=1,nchrg)
c      write(*,*)'idxclst'
c      do j1=1,nchrg
c         write(*,*)idxclst(1,j1)
c         write(*,*)(idxclst(j2,j1),j2=1,idxclst(1,j1)+1)
c      enddo
c      stop
c     end DBG
      return
      end

**************************************
**************************
      subroutine maptozero2(p1,p2,p3,p4,q1,q2,q3,q4,
     .     p1r,p2r,p3r,p4r,ierror)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension q1(0:3),q2(0:3),q3(0:3),q4(0:3)
      dimension p1r(0:3),p2r(0:3),p3r(0:3),p4r(0:3)
      dimension q(0:3),p(0:3),ptmp(0:3)
      double precision masses(4)
      character*2 fs
      common/momentainitial/pin1(0:3),pin2(0:3)
      common/mapto0_2first/ifirst,iprintwarn
      common/iclosest/iclose(40)
      common/finalstate/fs
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      data ifirst,iprintwarn /0,0/


      ierror = 0

      masses(1) = ame
      masses(2) = ame
      masses(3) = ame
      masses(4) = ame
      if (fs.eq.'gg') then
         masses(3) = 0.d0
         masses(4) = 0.d0
      endif
      if (fs.eq.'mm'.or.fs.eq.'mr') then
         masses(3) = ammu
         masses(4) = ammu
      endif
      if (fs.eq.'pp'.or.fs.eq.'pi') then
         masses(3) = ampi
         masses(4) = ampi
      endif

      do k=0,3
         q(k) = p1(k) + p2(k)
         p(k) = p3(k) + p4(k)
      enddo

      ic = iclose(1)
      shat = dot(q,q)
      sshat = sqrt(shat)    ! invariant

c      if (ic.eq.2.or.ic.eq.4) then
      if (ic.eq.1.or.ic.eq.3) then
c      if (ic.lt.10) then
         do i = 0,3
            ptmp(i) = p1(i) - p3(i) - q3(i)
         enddo
         that = dot(ptmp,ptmp)
         cstar = 1.d0+2.d0*that/shat

         if (cstar.le.-1.d0) then
            cstar = -1.d0 + 0.1d-8
         endif
         sstar = sqrt(1.d0-cstar**2)

         do i = 0,3
            ptmp(i) = p3(i) + q3(i)
         enddo
         ph = getphiangle(ptmp)

         sph = sin(ph)
         cph = cos(ph)

         p3r(0) = 1.d0
         p3r(1) = sph*sstar
         p3r(2) = cph*sstar
         p3r(3) = cstar

         p4r(0) = 1.d0
         p4r(1) = -p3r(1)
         p4r(2) = -p3r(2)
         p4r(3) = -p3r(3)

      else

         do i = 0,3
            ptmp(i) = p2(i) - p4(i) - q4(i)
         enddo
         that = dot(ptmp,ptmp)
         cstar = 1.d0+2.d0*that/shat
         if (cstar.le.-1.d0) cstar = -1.d0 + 0.1d-8

         cstar = -cstar     !!!!!!!!!

         sstar = sqrt(1.d0-cstar**2)

         do i = 0,3
            ptmp(i) = p4(i) + q4(i)
         enddo
         ph = getphiangle(ptmp)
         sph = sin(ph)
         cph = cos(ph)

         p4r(0) = 1.d0
         p4r(1) = sph*sstar
         p4r(2) = cph*sstar
         p4r(3) = cstar

         p3r(0) = 1.d0
         p3r(1) = -p4r(1)
         p3r(2) = -p4r(2)
         p3r(3) = -p4r(3)

      endif
c      sshat = sqrt(dot(q,q)) ! invariant

      e = sshat/2.d0
      ei   = 2.d0*e

      ea = (ei**2 + masses(1)**2 -masses(2)**2)/2.d0/ei
      eb = ei - ea


      if (ea.lt.masses(1).or.eb.lt.masses(2)) then
         iprintwarn = iprintwarn + 1
         if (iprintwarn.le.100.and.iprintwarn.gt.0) then
            print*,' WARNING IN MAPTOZERO'
            print*,' UNPHYSICAL CONFIGURATION, KILLING THIS EVENT!'
            print*,' stopping the warning after 100 prints'
            print*,' but keep killing events'
            print*,'ea,eb,masses',ea,eb,masses
         endif
         if (ea.lt.masses(1)) ea = masses(1)
         if (eb.lt.masses(2)) eb = masses(2)
         ierror = 1
      endif

      beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
      beta2 = sqrt(1.d0 - masses(2)**2/eb**2)

      p1r(0) = ea
      p1r(1) = 0.d0
      p1r(2) = 0.d0
      p1r(3) = ea*beta1

      p2r(0) =  eb
      p2r(1) =  0.d0
      p2r(2) =  0.d0
      p2r(3) = -eb*beta2

      abig = ei**2 + masses(3)**2 - masses(4)**2
      bbig = -2.d0 * ei
      cbig = 0.d0
      arg  = abig**2*bbig**2-masses(3)**2*bbig**2*(bbig**2-cbig**2)

      p1mod = abig*cbig + sqrt(arg)
      p1mod = p1mod / (bbig**2 - cbig**2)

      e1    = sqrt(masses(3)**2 + p1mod**2)
      beta3 = p1mod/e1

      p3mod = sqrt(tridot(p3r,p3r))
      p4mod = sqrt(tridot(p4r,p4r))

      p3r(0) = e1
      p3r(1) = e1*beta3 * p3r(1)/p3mod
      p3r(2) = e1*beta3 * p3r(2)/p3mod
      p3r(3) = e1*beta3 * p3r(3)/p3mod

      p4r(0) = ei - p3r(0)
      p4r(1) = -p3r(1)
      p4r(2) = -p3r(2)
      p4r(3) = -p3r(3)
      return
      end
**************************
      subroutine maptozero(p1,p2,p3,p4,p1r,p2r,p3r,p4r,ierror)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension p1r(0:3),p2r(0:3),p3r(0:3),p4r(0:3)
      dimension q(0:3),p(0:3)
      double precision masses(4)
      character*2 fs
      common/mapto0first/ifirst,iprintwarn
      common/finalstate/fs
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      data ifirst,iprintwarn /0,0/
      ierror = 0
      masses(1) = ame
      masses(2) = ame
      masses(3) = ame
      masses(4) = ame
      if (fs.eq.'gg') then
         masses(3) = 0.d0
         masses(4) = 0.d0
      endif
      if (fs.eq.'mm'.or.fs.eq.'mr') then
         masses(3) = ammu
         masses(4) = ammu
      endif
      if (fs.eq.'pp'.or.fs.eq.'pi') then
         masses(3) = ampi
         masses(4) = ampi
      endif

      do k=0,3
         q(k) = p1(k) + p2(k)
         p(k) = p3(k) + p4(k)
      enddo

      call new_boost(p,p3,p3r,1)
      call new_boost(p,p4,p4r,1)

      sshat = sqrt(dot(q,q)) ! invariant

      e = sshat/2.d0
      ei   = 2.d0*e

      ea = (ei**2 + masses(1)**2 -masses(2)**2)/2.d0/ei
      eb = ei - ea
c      print*,''
c      print*,p1
c      print*,p2
c      print*,p3
c      print*,p4

      if (ea.lt.masses(1).or.eb.lt.masses(2)) then
         iprintwarn = iprintwarn + 1
         if (iprintwarn.le.100.and.iprintwarn.gt.0) then
            print*,' WARNING IN MAPTOZERO'
            print*,' UNPHYSICAL CONFIGURATION, KILLING THIS EVENT!'
            print*,' stopping the warning after 100 prints'
            print*,' but keep killing events'
            print*,'ea,eb,masses',ea,eb,masses
         endif
         if (ea.lt.masses(1)) ea = masses(1)
         if (eb.lt.masses(2)) eb = masses(2)
         ierror = 1
      endif

      beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
      beta2 = sqrt(1.d0 - masses(2)**2/eb**2)

      p1r(0) = ea
      p1r(1) = 0.d0
      p1r(2) = 0.d0
      p1r(3) = ea*beta1

      p2r(0) =  eb
      p2r(1) =  0.d0
      p2r(2) =  0.d0
      p2r(3) = -eb*beta2

      abig = ei**2 + masses(3)**2 - masses(4)**2
      bbig = -2.d0 * ei
      cbig = 0.d0
      arg  = abig**2*bbig**2-masses(3)**2*bbig**2*(bbig**2-cbig**2)

      p1mod = abig*cbig + sqrt(arg)
      p1mod = p1mod / (bbig**2 - cbig**2)

      e1    = sqrt(masses(3)**2 + p1mod**2)
      beta3 = p1mod/e1

      p3mod = sqrt(tridot(p3r,p3r))
      p4mod = sqrt(tridot(p4r,p4r))

      p3r(0) = e1
      p3r(1) = e1*beta3 * p3r(1)/p3mod
      p3r(2) = e1*beta3 * p3r(2)/p3mod
      p3r(3) = e1*beta3 * p3r(3)/p3mod

      p4r(0) = ei - p3r(0)
      p4r(1) = -p3r(1)
      p4r(2) = -p3r(2)
      p4r(3) = -p3r(3)
      return
      end
*******************
**********************
      subroutine maptoone(is,bigq,p3,p4,q,p1r,p2r,p3r,p4r,qr)
      implicit double precision (a-h,o-z)
      dimension bigq(0:3),p3(0:3),p4(0:3)
      dimension p1r(0:3),p2r(0:3),p3r(0:3),p4r(0:3)
      dimension q(0:3),p(0:3),qr(0:3),ptmp(0:3),ptmp2(0:3)
      double precision masses(4)
      character*2 fs
      common/finalstate/fs
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      common/momentainitial/pin1(0:3),pin2(0:3)
      common/counter_mapptoone/icountmto
      data icountmto /0/
      icountmto = icountmto + 1
! is = 1 q is initial state photon
! is = 0 q is final state photon
! This is identical to WP
      masses(1) = ame
      masses(2) = ame
      masses(3) = ame
      masses(4) = ame
      if (fs.eq.'gg') then
         masses(3) = 0.d0
         masses(4) = 0.d0
      endif
      if (fs.eq.'mm'.or.fs.eq.'mr') then
         masses(3) = ammu
         masses(4) = ammu
      endif
      if (fs.eq.'pp'.or.fs.eq.'pi') then
         masses(3) = ampi
         masses(4) = ampi
      endif
      shat = dot(bigq,bigq)
      if (is.eq.1) then
         do k = 0,3
            p(k) = p3(k)+p4(k)
         enddo
         call new_boost(p,p3,p3r,1)
         call new_boost(p,p4,p4r,1)
         e = sqrt(shat)/2.d0

         ei   = 2.d0*e
         abig = ei**2 + masses(3)**2 - masses(4)**2
         bbig = -2.d0 * ei
         cbig = 0.d0
         arg  = abig**2*bbig**2-masses(3)**2*bbig**2*(bbig**2-cbig**2)

         p1mod = abig*cbig + sqrt(arg)
         p1mod = p1mod / (bbig**2 - cbig**2)

         e1    = sqrt(masses(3)**2 + p1mod**2)
         beta3 = p1mod/e1

         p3mod = sqrt(tridot(p3r,p3r))
         p4mod = sqrt(tridot(p4r,p4r))

         p3r(0) = e1
         p3r(1) = e1*beta3 * p3r(1)/p3mod
         p3r(2) = e1*beta3 * p3r(2)/p3mod
         p3r(3) = e1*beta3 * p3r(3)/p3mod

         p4r(0) = ei - p3r(0)
         p4r(1) = -p3r(1)
         p4r(2) = -p3r(2)
         p4r(3) = -p3r(3)
         ! p3r p4r to be intended in the bigq frame

         do k = 0,3
            p(k)    = bigq(k)+q(k)
            ptmp(k) = pin1(k)+pin2(k)-q(k)
         enddo

         scale = sqrt(dot(bigq,bigq)/dot(ptmp,ptmp))

         do k = 0,3
            qr(k) = scale*q(k)
            ptmp(k) = pin1(k)+pin2(k)
         enddo

         e = sqrt(dot(ptmp,ptmp))/2.d0*scale

         ei = 2.d0*e
         ea = (ei**2 + masses(1)**2 -masses(2)**2)/2.d0/ei
         eb = ei - ea
         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)

         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)

         p1r(0) = ea
         p1r(1) = 0.d0
         p1r(2) = 0.d0
         p1r(3) = ea*beta1

         p2r(0) = eb
         p2r(1) = 0.d0
         p2r(2) = 0.d0
         p2r(3) = -eb*beta2

         do k = 0,3
            ptmp(k) = p1r(k)+p2r(k)-qr(k)
         enddo
 !        call new_boost(ptmp,qr,qr,1)
         call new_boost(ptmp,p3r,p3r,-1)
         call new_boost(ptmp,p4r,p4r,-1)
      else

         e=sqrt(shat)/2.d0

         ei = 2.d0*e
         ea = (ei**2 + masses(1)**2 -masses(2)**2)/2.d0/ei
         eb = ei - ea
         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)

         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)

         p1r(0) = ea
         p1r(1) = 0.d0
         p1r(2) = 0.d0
         p1r(3) = ea*beta1

         p2r(0) =  eb
         p2r(1) =  0.d0
         p2r(2) =  0.d0
         p2r(3) = -eb*beta2

         do k = 0,3
            ptmp(k) = p3(k)+p4(k)
         enddo

         call new_boost(ptmp,p3,p3r,1)
         call new_boost(ptmp,p4,p4r,1)
         call new_boost(ptmp,q,qr,1)

         do k = 0,3
            p(k)    = bigq(k)-q(k)
            ptmp(k) = p3(k)+p4(k)+q(k)
         enddo

         scale = sqrt(dot(bigq,bigq)/dot(ptmp,ptmp))

         do k = 0,3
            qr(k)   = scale*qr(k)
            ptmp(k) = p4(k)+p3(k)
         enddo

         e = sqrt(dot(ptmp,ptmp))/2.d0*scale

         ei   =  2.d0*e
         abig =  ei**2 + masses(3)**2 - masses(4)**2
         bbig = -2.d0 * ei
         cbig =  0.d0
         arg  =  abig**2*bbig**2-masses(3)**2*bbig**2*(bbig**2-cbig**2)

         p1mod = abig*cbig + sqrt(arg)
         p1mod = p1mod / (bbig**2 - cbig**2)

         e1    = sqrt(masses(3)**2 + p1mod**2)
         beta3 = p1mod/e1

         p3mod = sqrt(tridot(p3r,p3r))
         p4mod = sqrt(tridot(p4r,p4r))

         p3r(0) = e1
         p3r(1) = e1*beta3 * p3r(1)/p3mod
         p3r(2) = e1*beta3 * p3r(2)/p3mod
         p3r(3) = e1*beta3 * p3r(3)/p3mod

         p4r(0) = ei - p3r(0)

         p4r(1) = -p3r(1)
         p4r(2) = -p3r(2)
         p4r(3) = -p3r(3)
         ! p3r p4r to be intended in the current p frame

         do k = 0,3
            ptmp(k) = p3r(k)+p4r(k)+qr(k)
         enddo

         call new_boost(ptmp,p3r,p3r,1)
         call new_boost(ptmp,p4r,p4r,1)
         call new_boost(ptmp,qr,qr,1)
      endif
      return
      end
************************************************************
************* tentative maptotwo subroutine ****************
************************************************************

      subroutine maptotwo(is,bigq,p3,p4,q1,q2,p1r,p2r,p3r,p4r,q1r,q2r)
      implicit none
      integer is, icountmto, icountmtt
      real*8 bigq(0:3),p3(0:3),p4(0:3),p1r(0:3),p2r(0:3),p3r(0:3)
      real*8 p4r(0:3),q1r(0:3),q2r(0:3),pin1(0:3),pin2(0:3), ptmp(0:3)
      real *8 scale, q1(0:3),q2(0:3),abig,arg,bbig,beta1,beta2,beta3
      real*8 cbig,e,e1,ea,eb,ei,p1mod,p3mod,p4mod,tridot,p(0:3)
      real*8 masses(4), ame, ammu, convfac,alpha,pi,ampi
      character*2 fs
      common/finalstate/fs
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      common/momentainitial/pin1,pin2
      common/counter_mapptoone/icountmto
      data icountmto /0/
      common/counter_mapptotwo/icountmtt
      data icountmtt /0/
      real*8 shat, dot
      external dot
      external tridot
      integer i
      icountmtt = icountmtt + 1
! is = 1  q1,q2 are two initial state photons
! is = 0  q1,q2 are two final state photons
! is = 10 q1 is initial and q2 is final
! is = 11 q2 is initial and q1 is final
      masses(1) = ame
      masses(2) = ame
      masses(3) = ame
      masses(4) = ame
      if (fs.eq.'gg') then
         masses(3) = 0.d0
         masses(4) = 0.d0
      endif
      if (fs.eq.'mm'.or.fs.eq.'mr') then
         masses(3) = ammu
         masses(4) = ammu
      endif
      if (fs.eq.'pp'.or.fs.eq.'pi'.or.fs.eq.'pr') then
         masses(3) = ampi
         masses(4) = ampi
      endif

      shat = dot(bigq,bigq)

      if(is.eq.1) then ! two photons come from the is
         do i=0,3
            p(i) = p3(i) + p4(i)
         enddo
         call new_boost(p,p3,p3r,1)
         call new_boost(p,p4,p4r,1)

         e = sqrt(shat)/2.d0
         ei = 2.d0 * e
         abig = ei**2 + masses(3)**2 - masses(4)**2
         bbig = -2.d0 * ei
         cbig = 0.d0
         arg  = abig**2*bbig**2 - masses(3)**2*bbig**2*(bbig**2-cbig**2)

         p1mod = abig * cbig + sqrt(arg)
         p1mod = p1mod / (bbig**2 - cbig**2)

         e1    = sqrt(masses(3)**2 + p1mod**2)
         beta3 = p1mod / e1

         p3mod = sqrt(tridot(p3r,p3r))
         p4mod = sqrt(tridot(p4r,p4r))

         p3r(0) = e1
         p3r(1) = e1*beta3 * p3r(1)/p3mod
         p3r(2) = e1*beta3 * p3r(2)/p3mod
         p3r(3) = e1*beta3 * p3r(3)/p3mod

         p4r(0) = ei - p3r(0)
         p4r(1) = -p3r(1)
         p4r(2) = -p3r(2)
         p4r(3) = -p3r(3)
! p3r p4r to be intended in the bigq frame

         do i = 0,3
            ptmp(i) = pin1(i) + pin2(i) - q1(i) - q2(i)
         enddo

         scale = sqrt(dot(bigq,bigq)/dot(ptmp,ptmp))

         do i = 0,3
            q1r(i) = scale * q1(i)
            q2r(i) = scale * q2(i)
            ptmp(i) = pin1(i) + pin2(i)
         enddo

         e = sqrt(dot(ptmp,ptmp))/2.d0*scale

         ei = 2.d0*e
         ea = (ei**2 + masses(1)**2 -masses(2)**2)/2.d0/ei
         eb = ei - ea
         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)

         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)

         p1r(0) = ea
         p1r(1) = 0.d0
         p1r(2) = 0.d0
         p1r(3) = ea*beta1

         p2r(0) = eb
         p2r(1) = 0.d0
         p2r(2) = 0.d0
         p2r(3) = -eb*beta2

         do i = 0,3
            ptmp(i) = p1r(i) + p2r(i) - q1r(i) - q2r(i)
         enddo

         call new_boost(ptmp,p3r,p3r,-1)
         call new_boost(ptmp,p4r,p4r,-1)

c$$$         print*, '++ is ++'
c$$$         print*, p1r
c$$$         print*, p2r
c$$$         print*, p3r
c$$$         print*, p4r
c$$$         print*, q1r
c$$$         print*, q2r
c$$$         print*, '-- is --'


!     two photons come from the fs
!     OR
!     one comes from is and one from fs
      elseif(is.eq.0.or.is.eq.10.or.is.eq.11) then

         e=sqrt(shat)/2.d0

         ei = 2.d0*e
         ea = (ei**2 + masses(1)**2 -masses(2)**2)/2.d0/ei
         eb = ei - ea
         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)

         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)

         p1r(0) = ea
         p1r(1) = 0.d0
         p1r(2) = 0.d0
         p1r(3) = ea*beta1

         p2r(0) =  eb
         p2r(1) =  0.d0
         p2r(2) =  0.d0
         p2r(3) = -eb*beta2

         do i = 0,3
            ptmp(i) = p3(i)+p4(i)
         enddo

         call new_boost(ptmp,p3,p3r,1)
         call new_boost(ptmp,p4,p4r,1)
         call new_boost(ptmp,q1,q1r,1)
         call new_boost(ptmp,q2,q2r,1)

         do i = 0,3
            p(i)    = bigq(i)-q1(i)-q2(i)
            ptmp(i) = p3(i)+p4(i)+q1(i)+q2(i)
         enddo

         scale = sqrt(dot(bigq,bigq)/dot(ptmp,ptmp))

         do i = 0,3
            q1r(i)   = scale*q1r(i)
            q2r(i)   = scale*q2r(i)
            ptmp(i)  = p3(i)+p4(i)
         enddo

         e = sqrt(dot(ptmp,ptmp))/2.d0*scale

         ei   =  2.d0*e
         abig =  ei**2 + masses(3)**2 - masses(4)**2
         bbig = -2.d0 * ei
         cbig =  0.d0
         arg  =  abig**2*bbig**2-masses(3)**2*bbig**2*(bbig**2-cbig**2)

         p1mod = abig*cbig + sqrt(arg)
         p1mod = p1mod / (bbig**2 - cbig**2)

         e1    = sqrt(masses(3)**2 + p1mod**2)
         beta3 = p1mod/e1

         p3mod = sqrt(tridot(p3r,p3r))
         p4mod = sqrt(tridot(p4r,p4r))

         p3r(0) = e1
         p3r(1) = e1*beta3 * p3r(1)/p3mod
         p3r(2) = e1*beta3 * p3r(2)/p3mod
         p3r(3) = e1*beta3 * p3r(3)/p3mod

         p4r(0) = ei - p3r(0)
         p4r(1) = -p3r(1)
         p4r(2) = -p3r(2)
         p4r(3) = -p3r(3)

         do i = 0,3
            ptmp(i) = p3r(i)+p4r(i)+q1r(i)+q2r(i)
         enddo

         call new_boost(ptmp,p3r,p3r,1)
         call new_boost(ptmp,p4r,p4r,1)
         call new_boost(ptmp,q1r,q1r,1)
         call new_boost(ptmp,q2r,q2r,1)

c$$$  	 print*, '++ fs ++'
c$$$	 print*, p1r
c$$$         print*, p2r
c$$$         print*, p3r
c$$$         print*, p4r
c$$$         print*, q1r
c$$$         print*, q2r
c$$$         print*, '-- fs --'

      endif
c$$$
c$$$      elseif(is.eq.10) then     ! one photon comes from is and one from fs
c$$$
c$$$         do i=0,3
c$$$            ptmp(i) = pin1(i) + pin2(i) - q1(i)
c$$$         enddo
c$$$         scale = sqrt((dot(bigq,bigq))/dot(ptmp,ptmp))
c$$$
c$$$	 do i=0,3
c$$$            q1r(i) = scale*q1(i)
c$$$            ptmp(i) = pin1(i) + pin2(i) - q1(i)
c$$$         enddo
c$$$
c$$$         e=sqrt(dot(ptmp,ptmp))/2.d0
c$$$
c$$$         ei = 2.d0*e
c$$$         ea = (ei**2 + masses(1)**2 -masses(2)**2)/2.d0/ei
c$$$         eb = ei - ea
c$$$         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
c$$$         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)
c$$$
c$$$         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
c$$$         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)
c$$$
c$$$         p1r(0) = ea
c$$$         p1r(1) = 0.d0
c$$$         p1r(2) = 0.d0
c$$$         p1r(3) = ea*beta1
c$$$
c$$$         p2r(0) =  eb
c$$$         p2r(1) =  0.d0
c$$$         p2r(2) =  0.d0
c$$$         p2r(3) = -eb*beta2
c$$$
c$$$         do i = 0,3
c$$$            ptmp(i) = p3(i)+p4(i) + q1(i)
c$$$         enddo
c$$$         call new_boost(ptmp,p3,p3r,1)
c$$$         call new_boost(ptmp,p4,p4r,1)
c$$$         call new_boost(ptmp,q2,q2r,1)
c$$$
c$$$         do i=0,3
c$$$            ptmp(i) = p3(i) + p4(i) + q2(i)! +q1(i)
c$$$         enddo
c$$$         scale = sqrt((dot(bigq,bigq))/dot(ptmp,ptmp))
c$$$
c$$$         print*, 'scale10', scale
c$$$
c$$$         do i=0,3
c$$$            q2r(i) = q2r(i) * scale
c$$$            ptmp(i) = p3(i) + p4(i) !+q2(i)
c$$$         enddo
c$$$
c$$$         e = sqrt(dot(ptmp,ptmp))/2.d0*scale
c$$$
c$$$         ei   =  2.d0*e
c$$$         abig =  ei**2 + masses(3)**2 - masses(4)**2
c$$$         bbig = -2.d0 * ei
c$$$         cbig =  0.d0
c$$$         arg  =  abig**2*bbig**2-masses(3)**2*bbig**2*(bbig**2-cbig**2)
c$$$
c$$$         p1mod = abig*cbig + sqrt(arg)
c$$$         p1mod = p1mod / (bbig**2 - cbig**2)
c$$$
c$$$         e1    = sqrt(masses(3)**2 + p1mod**2)
c$$$         beta3 = p1mod/e1
c$$$
c$$$         p3mod = sqrt(tridot(p3r,p3r))
c$$$         p4mod = sqrt(tridot(p4r,p4r))
c$$$
c$$$         p3r(0) = e1
c$$$         p3r(1) = e1*beta3 * p3r(1)/p3mod
c$$$         p3r(2) = e1*beta3 * p3r(2)/p3mod
c$$$         p3r(3) = e1*beta3 * p3r(3)/p3mod
c$$$
c$$$         p4r(0) = ei - p3r(0)
c$$$         p4r(1) = -p3r(1)
c$$$         p4r(2) = -p3r(2)
c$$$         p4r(3) = -p3r(3)
c$$$
c$$$         do i = 0,3
c$$$            ptmp(i) = p3r(i)+p4r(i)+q2r(i)
c$$$         enddo
c$$$
c$$$         call new_boost(ptmp,p3r,p3r,1)
c$$$         call new_boost(ptmp,p4r,p4r,1)
c$$$         call new_boost(ptmp,q2r,q2r,1)
c$$$
c$$$      elseif (is.eq.11) then    ! one photon comes from fs and one from is
c$$$
c$$$         do i=0,3
c$$$            ptmp(i) = pin1(i) + pin2(i) - q2(i)
c$$$         enddo
c$$$         scale = sqrt((dot(bigq,bigq))/dot(ptmp,ptmp))
c$$$
c$$$         do i=0,3
c$$$            q2r(i) = scale*q2(i)
c$$$            ptmp(i) = pin1(i) + pin2(i)
c$$$         enddo
c$$$
c$$$         e=sqrt(dot(ptmp,ptmp))/2.d0
c$$$
c$$$         ei = 2.d0*e
c$$$         ea = (ei**2 + masses(1)**2 -masses(2)**2)/2.d0/ei
c$$$         eb = ei - ea
c$$$         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
c$$$         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)
c$$$
c$$$         beta1 = sqrt(1.d0 - masses(1)**2/ea**2)
c$$$         beta2 = sqrt(1.d0 - masses(2)**2/eb**2)
c$$$
c$$$         p1r(0) = ea
c$$$         p1r(1) = 0.d0
c$$$         p1r(2) = 0.d0
c$$$         p1r(3) = ea*beta1
c$$$
c$$$         p2r(0) =  eb
c$$$         p2r(1) =  0.d0
c$$$         p2r(2) =  0.d0
c$$$         p2r(3) = -eb*beta2
c$$$
c$$$         do i = 0,3
c$$$            ptmp(i) = p3(i)+p4(i)
c$$$         enddo
c$$$         call new_boost(ptmp,p3,p3r,1)
c$$$         call new_boost(ptmp,p4,p4r,1)
c$$$         call new_boost(ptmp,q1,q1r,1)
c$$$
c$$$
c$$$         do i=0,3
c$$$            ptmp(i) = p3(i) + p4(i) + q1(i)
c$$$         enddo
c$$$         scale = sqrt((dot(bigq,bigq))/dot(ptmp,ptmp))
c$$$
c$$$         print*, 'scale11',scale
c$$$
c$$$         do i=0,3
c$$$            q1r(i) = q1r(i) * scale
c$$$            ptmp(i) = p3(i) + p4(i)
c$$$         enddo
c$$$
c$$$         e = sqrt(dot(ptmp,ptmp))/2.d0*scale
c$$$
c$$$         ei   =  2.d0*e
c$$$         abig =  ei**2 + masses(3)**2 - masses(4)**2
c$$$         bbig = -2.d0 * ei
c$$$         cbig =  0.d0
c$$$         arg  =  abig**2*bbig**2-masses(3)**2*bbig**2*(bbig**2-cbig**2)
c$$$
c$$$         p1mod = abig*cbig + sqrt(arg)
c$$$         p1mod = p1mod / (bbig**2 - cbig**2)
c$$$
c$$$         e1    = sqrt(masses(3)**2 + p1mod**2)
c$$$         beta3 = p1mod/e1
c$$$
c$$$         p3mod = sqrt(tridot(p3r,p3r))
c$$$         p4mod = sqrt(tridot(p4r,p4r))
c$$$
c$$$         p3r(0) = e1
c$$$         p3r(1) = e1*beta3 * p3r(1)/p3mod
c$$$         p3r(2) = e1*beta3 * p3r(2)/p3mod
c$$$         p3r(3) = e1*beta3 * p3r(3)/p3mod
c$$$
c$$$         p4r(0) = ei - p3r(0)
c$$$         p4r(1) = -p3r(1)
c$$$         p4r(2) = -p3r(2)
c$$$         p4r(3) = -p3r(3)
c$$$
c$$$         do i = 0,3
c$$$            ptmp(i) = p3r(i)+p4r(i)+q1r(i)
c$$$         enddo
c$$$
c$$$         call new_boost(ptmp,p3r,p3r,1)
c$$$         call new_boost(ptmp,p4r,p4r,1)
c$$$         call new_boost(ptmp,q1r,q1r,1)
c$$$
c$$$      endif

      return
      end

************************************************************
********** end of tentative maptotwo subroutine ************
************************************************************


      subroutine maptozero_BCK(p1,p2,p3,p4,p1r,p2r,p3r,p4r)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension p1r(0:3),p2r(0:3),p3r(0:3),p4r(0:3)
      dimension q(0:3),p(0:3)

      do k=0,3
         q(k) = p1(k) + p2(k)
         p(k) = p3(k) + p4(k)
      enddo

      call new_boost(p,p3,p3r,1)
      call new_boost(p,p4,p4r,1)

      amel = sqrt(dot(p3,p3))

      sshat = sqrt(dot(q,q)) ! invariant

      p3mod = sqrt(tridot(p3r,p3r))
      p4mod = sqrt(tridot(p4r,p4r))

      e = sshat/2.d0
      beta = sqrt(1.d0 - amel**2/e**2)

      p1r(0) = e
      p1r(1) = 0.d0
      p1r(2) = 0.d0
      p1r(3) = e*beta

      p2r(0) =  e
      p2r(1) =  0.d0
      p2r(2) =  0.d0
      p2r(3) = -e*beta

      p3r(0) = e
      p3r(1) = e*beta * p3r(1)/p3mod
      p3r(2) = e*beta * p3r(2)/p3mod
      p3r(3) = e*beta * p3r(3)/p3mod

      p4r(0) = e
      p4r(1) = e*beta * p4r(1)/p4mod
      p4r(2) = e*beta * p4r(2)/p4mod
      p4r(3) = e*beta * p4r(3)/p4mod
      return
      end
*******
      subroutine maptoone_BCK(is,bigq,p3,p4,q,p1r,p2r,p3r,p4r,qr)
      implicit double precision (a-h,o-z)
      dimension bigq(0:3),p3(0:3),p4(0:3)
      dimension p1r(0:3),p2r(0:3),p3r(0:3),p4r(0:3)
      dimension q(0:3),p(0:3),qr(0:3),ptmp(0:3)
! is = 1 q is initial state photon
! is = 0 q is final state photon

      shat = dot(bigq,bigq)
      amel = sqrt(dot(p3,p3))

      if (is.eq.1) then
         do k = 0,3
            p(k) = p3(k)+p4(k)
         enddo
         call new_boost(p,p3,p3r,1)
         call new_boost(p,p4,p4r,1)
         e = sqrt(shat)/2.d0
         beta = sqrt(1.d0 - amel**2/e**2)

         p3mod = sqrt(tridot(p3r,p3r))
         p4mod = sqrt(tridot(p4r,p4r))

         p3r(0) = e
         p3r(1) = e*beta * p3r(1)/p3mod
         p3r(2) = e*beta * p3r(2)/p3mod
         p3r(3) = e*beta * p3r(3)/p3mod

         p4r(0) = e
         p4r(1) = e*beta * p4r(1)/p4mod
         p4r(2) = e*beta * p4r(2)/p4mod
         p4r(3) = e*beta * p4r(3)/p4mod
         ! p3r p4r to be intended in the bigq frame

         call new_boost(bigq,p3r,p3r,-1)
         call new_boost(bigq,p4r,p4r,-1)

         do k = 0,3
            p(k) = bigq(k)+q(k)
         enddo

         e = sqrt(dot(p,p))/2.d0
         beta = sqrt(1.d0 - amel**2/e**2)

         p1r(0) = e
         p1r(1) = 0.d0
         p1r(2) = 0.d0
         p1r(3) = e*beta

         p2r(0) = e
         p2r(1) = 0.d0
         p2r(2) = 0.d0
         p2r(3) = -e*beta

         call new_boost(p,q,qr,1)
         call new_boost(p,p3r,p3r,1)
         call new_boost(p,p4r,p4r,1)
      else
         e=sqrt(shat)/2.d0
         beta = sqrt(1.d0 - amel**2/e**2)

         p1r(0) = e
         p1r(1) = 0.d0
         p1r(2) = 0.d0
         p1r(3) = e*beta

         p2r(0) = e
         p2r(1) =  0.d0
         p2r(2) =  0.d0
         p2r(3) = -e*beta

         do k = 0,3
            p(k) = bigq(k)-q(k)
         enddo
         e = sqrt(dot(p,p))/2.d0

         do k = 0,3
            ptmp(k) = p3(k)+p4(k)
         enddo
         call new_boost(ptmp,p3,p3r,1)
         call new_boost(ptmp,p4,p4r,1)

         beta = sqrt(1.d0 - amel**2/e**2)

         p3mod = sqrt(tridot(p3r,p3r))
         p4mod = sqrt(tridot(p4r,p4r))

         p3r(0) = e
         p3r(1) = e*beta * p3r(1)/p3mod
         p3r(2) = e*beta * p3r(2)/p3mod
         p3r(3) = e*beta * p3r(3)/p3mod

         p4r(0) = e
         p4r(1) = e*beta * p4r(1)/p4mod
         p4r(2) = e*beta * p4r(2)/p4mod
         p4r(3) = e*beta * p4r(3)/p4mod
         ! p3r p4r to be intended in the current p frame
         call new_boost(p,p3r,p3r,-1)
         call new_boost(p,p4r,p4r,-1)

         call new_boost(bigq,q,qr,1)
         call new_boost(bigq,p3r,p3r,1)
         call new_boost(bigq,p4r,p4r,1)
      endif
      return
      end



*****************MAP ROUTINE MAURO MORETTI

c***********************************************************************
      subroutine map_mauro(pclin,pclout,mssin,mssout,nout,pin,pout)
c***********************************************************************
      implicit none
c
c     INPUT: PCLIN(4,2), initial state clusters momentum == initial state
c            particle momentum MINUS the momentum of all photons clustered
c            with the particle;
c            PCLOUT(4,NOUT), final state clusters momentum == final state
c            particle momentum PLUS the momentum of all photons clustered
c            with the particle;
c            MSSIN(2), initial state particle masses
c            MSSOUT(NOUT), final state particle masses
c            NOUT number of outcoming particles (after clustering) N.B.
c                 used to dimension various arrays, should be a parameter
c                 in the calling routine
c
c     OUTPUT: PIN(0:3,2)     momentum of initial partons
c             POUT(0:3,NOUT) momentum of final partons
c
c                   convenction for momentum: p(4)==(px,py,pz,p0)
c
      integer nout
      real*8 pclin(4,2),pclout(4,nout),mssin(2),mssout(nout)
      real*8 poutm(0:3,2),pinm(0:3,2)
      real*8 pin(4,2),pout(4,nout),c1,c2,dlt,ppfsq(nout),c0,cms(4)
      real*8 dlt1,dlt2,fct1,fct2,fct,tric

c
      real*8 intrp

      parameter (intrp=1.d0) ! parameter to adjust initial state momenta.
c                              it interpolates linearly between
c                              INTERP=0 == initial partons keep the same
c                              momentum as the initial clusters, center of
c                              of mass energy of initial clusters 'slightly'
c                              affected
c                              INTERP=1 == initial partons keep the same
c                              center of mass energy of initial clusters,
c                              center of mass momentum of initial clusters
c                              slightly affected
c                              Any value of INTERP of order one affects kinematical
c                              invariants at order KT^2, KT being a measure of
c                              photons acollinearity. The approximate dynamics
c                              is accurate up to KT terms and thus all choices are
c                              equivalent. Varying INTERP could provide an handle
c                              to assess the impact of the approximation.
c                              The statement about the accuracy of the mapping
c                              being accurate up to KT^2 terms can be shown
c                              to be true for one emission onto a 2->2 process.
c                              It is REASONABLE/LIKELY it holds for more complicated
c                              setup.
c
      real*8 pis(4),ein,ppis,tmpin(4,2),tmpout(4,nout),pisrv(4),pscnsv
      real*8 pimpcns,ecmin,ecmfin,sq,acc,tmp1,tmp2,fcttot,dltold
      integer j1,j2,j3,flgb
c

      ppis=0.d0
      do j1=1,4
         pis(j1)=pclin(j1,1)+pclin(j1,2) ! initial state momentum
         if(j1.ne.4)then
            pisrv(j1)=-pis(j1)       ! for inverse boost
            ppis=ppis+abs(pis(j1))
         else
            pisrv(j1)=pis(j1)
         endif
      enddo
c
c     DBG
c         write(*,*)'initial state four momentum'
c         write(*,*)(pis(j1),j1=1,4)
c end DBG
      if(ppis.le.(1.d-40)) then    ! no ISR, no boost required
         flgb=1
         do j1=1,4
            pin(j1,1)=pclin(j1,1)
            pin(j1,2)=pclin(j1,2)
            tmpin(j1,1)=pin(j1,1)
            tmpin(j1,2)=pin(j1,2)
         enddo
      else                     ! boost to center of mass and 'adjust' initial state
c                                momenta
         flgb=-1
         sq = pis(4)**2-pis(1)**2-pis(2)**2-pis(3)**2 !invariant mass squared of the
c                                                      initial state clusters
         pscnsv=0.5d0*sqrt(sq-2.d0*(mssin(1)**2+mssin(2)**2)+
     .        (mssin(1)**2-mssin(2)**2)**2/sq) ! partons momentum required to preserve
c                                              center of mass energy
c
c     DBG
c         write(*,*)'in cluster momenta'
c         write(*,*)(pclin(j1,1),j1=1,4)
c         write(*,*)(pclin(j1,2),j1=1,4)
c end DBG
         do j1=1,2

            call boost_mauro(pis,pclin(1,j1),tmpin(1,j1)) !boost in the PIS rest frame
            if(j1.eq.1) then
               pimpcns=sqrt(tmpin(1,j1)**2+tmpin(2,j1)**2+
     .              tmpin(3,j1)**2) ! partons momentum required to preserve clusters
c                                     momentum
            endif
            do j2=1,3
               tmpin(j2,j1)=tmpin(j2,j1)*
     .                  (1.d0+intrp*(pscnsv/pimpcns-1.d0))
            enddo
c     DBG
c            write(*,*)'boosted initial cluster #',j1
c            write(*,*)(tmpin(j2,j1),j2=1,4)
c end DBG
            tmpin(4,j1)=sqrt(tmpin(1,j1)**2+tmpin(2,j1)**2+
     .           tmpin(3,j1)**2+mssin(j1)**2)       !enforce mass shell
c     DBG
c            write(*,*)'boosted initial cluster #',j1,
c     .                           '  after energy adjustment'
c            write(*,*)(tmpin(j2,j1),j2=1,4)
c end DBG
         enddo
c  DBG
c         write(*,*)'in state s', pis(4)**2-pis(1)**2-pis(2)**2
c     .        -pis(3)**2
c         write(*,*)'adjusted in state s',(tmpin(4,1)+tmpin(4,2))**2-
c     .        (tmpin(1,1)+tmpin(1,2))**2-(tmpin(3,1)+tmpin(3,2))**2
c     .        -(tmpin(2,1)+tmpin(2,2))**2
c end DBG
      endif
c
c     Building a set of final state momenta that preserves the angular distribution
c     of final state clusters
c

      ecmin=tmpin(4,1)+tmpin(4,2)
      ecmfin=0.d0
      if(flgb.eq.-1) then
         do j1=1,nout           ! 'adjust' final state momenta
            call boost_mauro(pis,pclout(1,j1),tmpout(1,j1)) ! boost in the PIS rest frame
            tmpout(4,j1)=sqrt(tmpout(1,j1)**2+tmpout(2,j1)**2+
     .           tmpout(3,j1)**2+mssout(j1)**2) !enforce mass shell
            ecmfin=ecmfin+tmpout(4,j1)
         enddo
      else
         do j1=1,nout
            do j2=1,3
               tmpout(j2,j1)=pclout(j2,j1)
            enddo
            tmpout(4,j1)=sqrt(tmpout(1,j1)**2+tmpout(2,j1)**2
     .           +tmpout(3,j1)**2+mssout(j1)**2)
            ecmfin=ecmfin+tmpout(4,j1)
         enddo
      endif
c     DBG
c      write(*,*)'input',(pclout(j1,1),j1=1,4)
c      write(*,*)'output',(tmpout(j1,1),j1=1,4)
c      write(*,*)'input',(pclout(j1,2),j1=1,4)
c      write(*,*)'output',(tmpout(j1,2),j1=1,4)
c      write(*,*)'adjusted in state energy',ecmin
c      write(*,*)'adjusted out state energy',ecmfin
c end DBG
c
c     momentum is now conserved (and 0 since we are in the rest frame), but energy
c     of initial and final state is different.
c     We adjust final state momenta as pf -> pf*fct to ensure energy conservation.
c     With many particles and different masses fct cannot be determined analitically
c     we iterate a second order expansion in fct to achieve the desired accuracy
c
      j1=0
      acc=abs(ecmin-ecmfin)/(ecmin+ecmfin)
      do j2=1,nout
         ppfsq(j2)=tmpout(1,j2)**2+tmpout(2,j2)**2+tmpout(3,j2)**2
      enddo
      fcttot=1.d0
      dltold=0.d0
c     DBG
c      write(*,*)'acc',acc,'  ppfsq',(ppfsq(j2),j2=1,nout)
c end DBG
      do while((j1.le.20).and.acc.ge.1.d-13)
         j1=j1+1
         c1=0.d0
         c2=0.d0
         do j2=1,nout
            c1=c1+ppfsq(j2)/sqrt(ppfsq(j2)+mssout(j2)**2)
            c2=c2+0.5d0*mssout(j2)**2*
     .                     ppfsq(j2)/sqrt(ppfsq(j2)+mssout(j2)**2)**3
         enddo
         c0=-ecmin+ecmfin
         dlt=c1**2-4.d0*c0*c2
c     DBG
c         write(*,*)'dlt',dlt
c end DBG
         if(dlt.lt.0.d0) then
            fct=1.d0-0.5d0*dltold
            dltold=0.5d0*dltold               !should ensure a sort of bisection
         else
            dlt=sqrt(dlt)
            dlt1=0.5d0/c2*(-c1+dlt)
            dlt2=0.5d0/c2*(-c1-dlt)
            fct1=1.d0+dlt1
            fct2=1.d0+dlt2
            tmp1=0.d0
            tmp2=0.d0
c     DBG
c            write(*,*)'fct1, fct2',fct1,fct2
c end DBG
            do j2=1,nout
               tmp1=tmp1+sqrt((fct1**2*ppfsq(j2)+mssout(j2)**2))
               tmp2=tmp1+sqrt((fct2**2*ppfsq(j2)+mssout(j2)**2))
            enddo
c     DBG
c            write(*,*)'tmp1, tmp2, ecmin',tmp1,tmp2,ecmin
c end DBG
            if(abs(tmp1-ecmin).lt.abs(tmp2-ecmin)) then
               fct=fct1
               dltold=dlt1
            else
               fct=fct2
               dltold=dlt2
            endif
         endif
         fcttot=fcttot*fct
c     DBG
c         write(*,*)'fct, fcttot',fct,fcttot
c end DBG
         ecmfin=0.d0
         do j2=1,nout
            ppfsq(j2)=ppfsq(j2)*fct**2
            ecmfin=ecmfin+sqrt(mssout(j2)**2+ppfsq(j2))
         enddo
         acc=abs(ecmin-ecmfin)/(ecmin+ecmfin)
c     DBG
c         write(*,*)'ecmfin,ecmin,acc',ecmfin,ecmin,acc
c end DBG
      enddo
c
      ecmfin=0.d0
      do j1=1,nout
         do j2=1,3
            tmpout(j2,j1)=fcttot*tmpout(j2,j1)
         enddo
         tmpout(4,j1)=sqrt(mssout(j1)**2+tmpout(1,j1)**2
     .        +tmpout(2,j1)**2+tmpout(3,j1)**2)
         ecmfin=ecmfin+tmpout(4,j1)
      enddo
c     DBG
c      write(*,*)'j1,ecmin,ecmfin',j1,ecmin,ecmfin
c end DBG
c
      if(j1.gt.20) then
         write(*,*)'failed prescribed accuracy in subroutin MAP'
         write(*,*)'pclin'
         write(*,*)pclin
         write(*,*)'pclout'
         write(*,*)pclout
      endif
c
c     final refinement to achieve real*8 accuracy on momentum conservation
      pscnsv=0.5d0*sqrt(ecmfin**2-2.d0*(mssin(1)**2+mssin(2)**2)+
     .     (mssin(1)**2-mssin(2)**2)**2/ecmfin**2) ! partons momentum required to
c                                                    match final state
c                                                    center of mass energy
      ppis=sqrt(tmpin(1,1)**2+tmpin(2,1)**2+tmpin(3,1)**2)
      do j1=1,nout
         tmpin(j1,1)=tmpin(j1,1)*pscnsv/ppis
         tmpin(j1,2)=tmpin(j1,2)*pscnsv/ppis
      enddo
      tmpin(4,1)=sqrt(mssin(1)**2+pscnsv**2)
      tmpin(4,2)=sqrt(mssin(2)**2+pscnsv**2)
c     DBG
c      write(*,*)'in state momentum 1',(tmpin(j1,1),j1=1,4)
c      write(*,*)'in state momentum 2',(tmpin(j1,2),j1=1,4)
c      write(*,*)'in state momentum',(tmpin(j1,1)+tmpin(j1,2),j1=1,4)
c      write(*,*)'out state momentum',(tmpout(j1,1)+tmpout(j1,2),j1=1,4)
c end DBG
c     DBG
c         stop
c end DBG
c
      if(flgb.eq.-1) then       !back to the lab frame. NB, if INTRP is not 1
c                                the boost of the fs will be different from the
c                                clusters' one, better check we are not losing
c                                accuracy, maybe PISRV(4) needs to be changed
c                                WARNING !!!!
c
         do j1=1,2
            call boost_mauro(pisrv,tmpin(1,j1),pin(1,j1))
         enddo
         do j1=1,nout
            call boost_mauro(pisrv,tmpout(1,j1),pout(1,j1))
         enddo
      else
         do j1=1,4
            pin(j1,1)=tmpin(j1,1)
            pin(j1,2)=tmpin(j1,2)
            do j2=1,nout
               pout(j1,j2)=tmpout(j1,j2)
            enddo
         enddo
      endif
c

      cms=pin(:,1)+pin(:,2)

      if(abs(cms(1)).ge.1d-20.or.abs(cms(2)).ge.1d-20) then
      call boost_mauro(cms,pin(:,1),pin(:,1))
         call boost_mauro(cms,pin(:,2),pin(:,2))
            do j2=1,nout
               call boost_mauro(cms,pout(:,j2),pout(:,j2))
            enddo
      endif
      return
      end

c***********************************************************************
      subroutine boost_mauro(pboost,pin,pout)
c***********************************************************************
      implicit none
c
c     INPUT: PBOOST(4) (px,py,pz,p0), the routine perform a boost which
c     brings PBOOST in the rest frame (NB it keeps the transverse
c                                      component of PIN untouched);
c            PIN(4), four momentum to be boosted
c     OUTPUt POUT(4), boosted PIN
c
      real*8 pboost(4),pin(4)
      real*8 pout(4)
c
      real*8 mb,pt(3),pl(3),plpb,pi,ppl,chb,shb,pb
      integer j1
c
      pb=pboost(3)**2+pboost(2)**2+pboost(1)**2
      mb=pboost(4)**2-pb
      if (mb.le.0.d0) then
         write(*,*)'mb <=0 in subroutine BOOST'
         stop
      endif
      mb=sqrt(mb)
      pb=sqrt(pb)
c
      pi=sqrt(pin(1)**2+pin(2)**2+pin(3)**2)
      ppl=(pin(1)*pboost(1)+pin(2)*pboost(2)+pin(3)*pboost(3))/(pi*pb)
      do j1=1,3
         pl(j1)=pboost(j1)*ppl/pb*pi
         pt(j1)=pin(j1)-pl(j1)
      enddo
c
      ppl=ppl*pi
      chb=pboost(4)/mb
      shb=-pb/mb
      pout(4)=chb*pin(4)+shb*ppl
      do j1=1,3
         pout(j1)=pt(j1)+pl(j1)*(chb*ppl+shb*pin(4))/ppl
      enddo
c
      return
      end
c***********************************************************************





      subroutine maptoone_rad(pin1cl,pin2cl,p1cl,p2cl,q,closeto,
     >			      ind,pin1br,pin2br,p1br,p2br,qbr)
      implicit double precision (a-h,o-z)
      parameter (nout5=3)
      parameter (nout6=4)
      dimension pclb(4,0:3),pcla(4,0:3)
      dimension q(0:3),p1(0:3),p2(0:3),p3(0:3),p4(0:3),qq(0:3)
      integer ind,is
      dimension pin1cl(0:3),pin2cl(0:3),p1cl(0:3),p2cl(0:3)
      dimension qclm(0:3),qclp(0:3),pmapin(4,2),pmapout(4,3)
      dimension pmapout6(4,4),qa(0:3)
      dimension pin1cl1(4),pin2cl1(4),p1cl1(4),p2cl1(4),qcl(4)
      dimension pclin(4,2),pclout(4,3),pin(4,2),pout(4,3),pout6(4,4)
      dimension pin1br(0:3),pin2br(0:3),p1br(0:3),p2br(0:3),qbr(0:3)
      double precision massesin(2),massesout(3),massesout6(4)
      common/finalstate/fs
      character*2 fs
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      common/momentainitial/pin1(0:3),pin2(0:3)

!#####here i fill some variables i need in the mapping
	massesin(1)=ame
	massesin(2)=ame
	massesout(1)=ame
	massesout(2)=ame
	massesout(3)=0d0
	massesout6(1)=ame
	massesout6(2)=ame
	massesout6(3)=0d0
	massesout6(4)=0d0

	if(fs.eq.'mr'.or.fs.eq.'mm') then
		massesout(1)=ammu
		massesout(2)=ammu
		massesout6(1)=ammu
		massesout6(2)=ammu
	elseif(fs.eq.'pr'.or.fs.eq.'pp') then
		massesout(1)=ampi
		massesout(2)=ampi
		massesout6(1)=ampi
		massesout6(2)=ampi
	endif






!IF QCLM IS ZERO THEN THE MAPPING ROUTINE OF MAURO WORKS PERFECTLY, ALL THE ADDITIONAL PHOTONS W.R.T. THE SIGNAL ONE ARE EXTERNAL, HENCE MUST BE CLUSTERED ON THE EXTERNAL LEGS
      if(ind.eq.1) then
         call conv_mauro(pin1cl,pin2cl,p1cl,p2cl,q,pin,pout)

         call map_mauro(pin,pout,massesin,massesout,
     >	    nout5,pmapin,pmapout)

	call conv_mom(pmapin,pmapout,pin1br,pin2br,p1br,p2br,qbr)
	return
      endif

! ON THE OTHER HAND IF THE PHOTONS ARE INTERNAL W.R.T. THE SIGNAL ONE, CARLO MAPPING WORKS BETTER
      if(ind.eq.2) then
        is=0
        if(closeto.eq.1.or.closeto.eq.2) is=1

	qq=pin1cl+pin2cl-is*(q)

        call maptoone(is,qq,p1cl,p2cl,q,
     .        pin1br,pin2br,p1br,p2br,qbr)

	return
      endif



      return
      end



      subroutine maptotwo_rad(pin1cl,pin2cl,p1cl,p2cl,q,q1,
     >                     pin1br,pin2br,p1br,p2br,qbr,q1br)
      implicit double precision (a-h,o-z)
      parameter (nout5=3)
      parameter (nout6=4)
      dimension pclb(4,0:3),pcla(4,0:3)
      dimension q(0:3),p1(0:3),p2(0:3),p3(0:3),p4(0:3),q1(0:3)
      integer ind,is
      dimension pin1cl(0:3),pin2cl(0:3),p1cl(0:3),p2cl(0:3)
      dimension qclm(0:3),qclp(0:3),pmapin(4,2),pmapout(4,3)
      dimension pmapout6(4,4),qa(0:3)
      dimension pin1cl1(4),pin2cl1(4),p1cl1(4),p2cl1(4),qcl(4)
      dimension pclin(4,2),pclout(4,3),pin(4,2),pout(4,3),pout6(4,4)
      dimension pin1br(0:3),pin2br(0:3),p1br(0:3),p2br(0:3),qbr(0:3)
      dimension q1br(0:3)
      double precision massesin(2),massesout(3),massesout6(4)
      common/finalstate/fs
      character*2 fs
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      common/momentainitial/pin1(0:3),pin2(0:3)

!#####here i fill some variables i need in the mapping
	massesin(1)=ame
	massesin(2)=ame
	massesout(1)=ame
	massesout(2)=ame
	massesout(3)=0d0
	massesout6(1)=ame
	massesout6(2)=ame
	massesout6(3)=0d0
	massesout6(4)=0d0

	if(fs.eq.'mr'.or.fs.eq.'mm') then
		massesout(1)=ammu
		massesout(2)=ammu
		massesout6(1)=ammu
		massesout6(2)=ammu
	elseif(fs.eq.'pr'.or.fs.eq.'pp') then
		massesout(1)=ampi
		massesout(2)=ampi
		massesout6(1)=ampi
		massesout6(2)=ampi
	endif

 	call conv_mom_mauro(pin1cl,pin(:,1))
        call conv_mom_mauro(pin2cl,pin(:,2))
        call conv_mom_mauro(p1cl,pout6(:,1))
        call conv_mom_mauro(p2cl,pout6(:,2))
        call conv_mom_mauro(q,pout6(:,3))
        call conv_mom_mauro(q1,pout6(:,4))

******QUI MAPPO
        call map_mauro(pin,pout6,massesin,massesout6,
     >       nout6,pmapin,pmapout6)

        call conv_mom_mauro1(pmapin(:,1),pin1br)
        call conv_mom_mauro1(pmapin(:,2),pin2br)
        call conv_mom_mauro1(pmapout6(:,1),p1br)
        call conv_mom_mauro1(pmapout6(:,2),p2br)
        call conv_mom_mauro1(pmapout6(:,3),qbr)
        call conv_mom_mauro1(pmapout6(:,4),q1br)





      return
      end

