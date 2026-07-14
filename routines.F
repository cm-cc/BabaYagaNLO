*******************************************************************************
      subroutine initrng(lux,iseed,k1,k2)
      integer iseed,lux,k1,k2,ilux
c      parameter (lux = 4,k1 = 0,k2 = 0) ! for RANLUX

c      call rluxgo(lux,iseed,k1,k2)
c      print*,'But c-ranlux used...'

      ilux = 2
      call drnlx_c_init(ilux,iseed) ! for C version of ranlux
      call srnlx_c_init(ilux,iseed) ! for C version of ranlux
      return
      end
*******************************************************************************
      subroutine getrndstatus(n)
      integer n(104)
      call rnlxdstate_c(n)
      return
      end
*******************************************************************************
      subroutine rndreset(n)
      integer n(104)
      call rnlxdreset_c(n)
      return
      end
************************************************
      subroutine getrnd(csi8,n)
c      subroutine ranlux(csi4,n)
! it can be modified by the user
      integer n
      double precision csi8(n)
      real*4 csi4(n)
c      call ranlux(csi4,n)
c      call srnlx_c(csi4,n) ! for C version of ranlux
c      csi8 = 1.d0 *csi4
c      return
      call drnlx_c(csi8,n)      ! for C version of ranlux

      return
      end
************************************************
      subroutine getdayshuursminutesseconds(sin,d,h,m,s)
      implicit none
      integer sin,d,h,m,s,sh,dins,hins,mins
      parameter (dins = 86400)
      parameter (hins = 3600)
      parameter (mins = 60)
      sh = sin
      d  = sh/dins
      sh = sh - d*dins
      h  = sh / hins
      sh = sh - h*hins
      m  = sh/mins
      s  = sh - m*mins
      return
      end
*****************************
      function check_invariants(p1,p2,p3,p4,q)
      implicit double precision (x)
      real*8 qph(40,0:3),p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      external dot
      real*8 dot
      common/checkinvariants/lim,switchcollier
      integer switchcollier
      double precision lim

      common/lt/ltint,quadevents,ilt
      integer ltint,quadevents,ilt


      x12  =dot(p1,p2)
      x13  =dot(p1,p3)
      x14  =dot(p1,p4)
      x15  =dot(p1,q)
      x23  =dot(p2,p3)
      x24  =dot(p2,p4)
      x25 = dot(p2,q)
      x34 = dot(p3,p4)
      x35 = dot(p3,q)
      x45 = dot(p4,q)

      if((abs(x12).le.lim.or.
     > abs(x13).le.lim.or.
     > abs(x14).le.lim.or.
     > abs(x15).le.lim.or.
     > abs(x23).le.lim.or.
     > abs(x24).le.lim.or.
     > abs(x25).le.lim.or.
     > abs(x34).le.lim.or.
     > abs(x35).le.lim.or.
     > abs(x45).le.lim)) then
         ltint=1
         quadevents=quadevents+1
!         ltint=0
      endif


      end
*******************************************************************************
      subroutine getsignalphot(p1,p2,p3,p4,qph,ng,ind,qind)
      implicit none
      real*8 qph(40,0:3),p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      real*8 enhard,spd,spdi,c1,c2,c3,c4,enphot
      integer k,i,ind,qind,ng
      external dot
      real*8 dot
      external tridot
      real*8 tridot
      real*8 sang1,sang2,sang3,sang4
      common/cutphot/phind
      integer phind

       qind =1

*******THE ME PHOTON IS CHOSEN AS THE HARDEST

        if(ind.eq.1) then
         enhard =0d0

         do i = 1,ng
            if(qph(i,0).gt.0.d0) then
               enphot = qph(i,0)
               if (enphot.ge.enhard) then
                  enhard = enphot
                  qind = i
               else

               endif
            endif
         enddo
         return
	endif
******THE ME PHOTON IS CHOSEN AS THE SIGNAL (IF MORE THAN ONE PASSES THE CUT, WE CHOOSE THE HARDEST(DONE IN CUTS))

	if(ind.eq.2) qind=phind

******THE ME PHOTON IS CHOSE AS THE FURTHEST FROM THE CHARGED LEGS
	if(ind.eq.3) then
	spd=0d0

        do k = 1,ng

		q = qph(k,:)

		c1 = dot(p1,q)
		c2 = dot(p2,q)
		c3 = dot(p3,q)
		c4 = dot(p4,q)

		if (c1.lt.c2.and.c1.lt.c3.and.c1.lt.c4) spdi =c1
		if (c2.lt.c1.and.c2.lt.c3.and.c2.lt.c4) spdi =c2
		if (c3.lt.c2.and.c3.lt.c1.and.c3.lt.c4) spdi =c3
		if (c4.lt.c1.and.c4.lt.c2.and.c4.lt.c3) spdi =c4
		if(spdi.ge.spd) then
		   spd=spdi
		   qind=k
		endif
         enddo
	 return
	 endif
	return
	end
***********
      subroutine cross_product(p, q, r)
      implicit none
      real(8) p(0:3), q(0:3) ! due quadrivettori in input
      real(8) r(0:3) ! risultato: prodotto vettoriale delle parti spaziali

! calcolo del prodotto vettoriale solo sulle componenti 1:3
      r(0) = 0
      r(1) = p(2)*q(3) - p(3)*q(2)
      r(2) = p(3)*q(1) - p(1)*q(3)
      r(3) = p(1)*q(2) - p(2)*q(1)

      end

*******************
      subroutine solvesystem(E1,E2,s1,s2,Qx,Qy,phi1,phi2)
      implicit none
      real*8 E1,E2,s1,s2,Qx,Qy,phi1,phi2
      real*8 A1,A2,A3,B3,xi(1)

      A1 = E1*s1
      A2 = E2*s2
      A3 = Qx
      B3 = Qy

      call getrnd(xi,1)
      if(xi(1).le.0.5d0) then
      phi1=atan2((-A1**2 + A2**2 - A3**2 - B3**2 + (A1**4*A3**2)
     >     /(A1**2*A3**2 + A1**2*B3**2) - (A1**2*A2**2*A3**2)
     >     /(A1**2*A3**2 + A1**2*B3**2) + (A1**2*A3**4)/(A1**2*A3**2 +
     >  A1**2*B3**2) + (A1**2*A3**2*B3**2)/(A1**2*A3**2 + A1**2*B3**2) +
     >     (A1*A3*dsqrt(-(A1**6*B3**2) + 2*A1**4*A2**2*B3**2 -
     >     A1**2*A2**4*B3**2 + 2*A1**4*A3**2*B3**2 + 2*A1**2*A2**2*A3**2
     >     *B3**2 - A1**2*A3**4*B3**2 + 2.d0*A1**4*B3**4 +
     >     2*A1**2*A2**2*B3**4 - 2*A1**2*A3**2*B3**4 - A1**2*B3**6))
     >     /(A1**2*A3**2 + A1**2*B3**2))/(2.d0*A1*B3),(-(A1**3*A3) +
     >     A1*A2**2*A3 - A1*A3**3 - A1*A3*B3**2
     > -dsqrt(-(A1**6*B3**2) + 2.d0*A1**4*A2**2*B3**2 -A1**2*A2**4*B3**2
     >     + 2*A1**4*A3**2*B3**2 + 2*A1**2*A2**2*A3**2*B3**2 -
     >     A1**2*A3**4*B3**2 + 2*A1**4*B3**4 + 2*A1**2*A2**2*B3**4
     >     - 2*A1**2*A3**2*B3**4 - A1**2*B3**6))/(2.d0*(A1**2*A3**2 +
     >     A1**2*B3**2)))


      phi2=atan2(-0.5*(-A1**2 + A2**2
     >     - A3**2 + B3**2 + (A1**4*A3**2)/(A1**2*A3**2 + A1**2*B3**2)
     >     - (A1**2*A2**2*A3**2)/(A1**2*A3**2 + A1**2*B3**2) +
     >     (A1**2*A3**4)/(A1**2*A3**2 + A1**2*B3**2) +
     >     (A1**2*A3**2*B3**2)/(A1**2*A3**2 + A1**2*B3**2) +
     >     (A1*A3*Sqrt(-(A1**6*B3**2) + 2*A1**4*A2**2*B3**2 -
     >     A1**2*A2**4*B3**2 + 2*A1**4*A3**2*B3**2 + 2*A1**2*A2**2*A3**2
     >     *B3**2 - A1**2*A3**4*B3**2 + 2*A1**4*B3**4 +
     >     2*A1**2*A2**2*B3**4 - 2*A1**2*A3**2*B3**4 - A1**2*B3**6))
     >     /(A1**2*A3**2 + A1**2*B3**2))/(A2*B3),
     >     (-A3 + (A1**4*A3)/(2.*(A1**2*A3**2 + A1**2*B3**2)) -
     >     (A1**2*A2**2*A3)/(2.*(A1**2*A3**2 + A1**2*B3**2)) +
     >     (A1**2*A3**3)/(2.*(A1**2*A3**2 + A1**2*B3**2)) +
     >     (A1**2*A3*B3**2)/(2.*(A1**2*A3**2 + A1**2*B3**2)) +
     >     (A1*Sqrt(-(A1**6*B3**2) + 2*A1**4*A2**2*B3**2 -
     >     A1**2*A2**4*B3**2 + 2*A1**4*A3**2*B3**2 +
     >     2*A1**2*A2**2*A3**2*B3**2 - A1**2*A3**4*B3**2 + 2*A1**4*B3**4
     >     + 2*A1**2*A2**2*B3**4 - 2*A1**2*A3**2*B3**4 - A1**2*B3**6))
     >     /(2.*(A1**2*A3**2 + A1**2*B3**2)))/A2)

      else
      phi1 = atan2((-A1**2 + A2**2 - A3**2 - B3**2 +(A1**4*A3**2)/
     >  (A1**2*A3**2 + A1**2*B3**2) - (A1**2*A2**2*A3**2)/(A1**2*A3**2 +
     >     A1**2*B3**2) + (A1**2*A3**4)/(A1**2*A3**2 + A1**2*B3**2)
     >     + (A1**2*A3**2*B3**2)/(A1**2*A3**2 + A1**2*B3**2) -
     >     (A1*A3*dsqrt(-(A1**6*B3**2) + 2*A1**4*A2**2*B3**2 -
     >     A1**2*A2**4*B3**2 + 2*A1**4*A3**2*B3**2 +
     >   2*A1**2*A2**2*A3**2*B3**2 - A1**2*A3**4*B3**2 + 2*A1**4*B3**4 +
     >     2*A1**2*A2**2*B3**4 - 2*A1**2*A3**2*B3**4 - A1**2*B3**6))
     >     /(A1**2*A3**2 + A1**2*B3**2))/(2.d0*A1*B3),(-(A1**3*A3) +
     >     A1*A2**2*A3 - A1*A3**3 - A1*A3*B3**2 +
     >  dsqrt(-(A1**6*B3**2) + 2*A1**4*A2**2*B3**2 - A1**2*A2**4*B3**2 +
     >     2*A1**4*A3**2*B3**2 + 2*A1**2*A2**2*A3**2*B3**2 -
     >     A1**2*A3**4*B3**2 + 2*A1**4*B3**4 + 2*A1**2*A2**2*B3**4 -
     >     2*A1**2*A3**2*B3**4 - A1**2*B3**6))/(2.*(A1**2*A3**2 +
     >     A1**2*B3**2)))

      phi2=atan2(-0.5*(-A1**2 + A2**2 -
     >     A3**2 + B3**2 + (A1**4*A3**2)/(A1**2*A3**2 + A1**2*B3**2) -
     >     (A1**2*A2**2*A3**2)/(A1**2*A3**2 + A1**2*B3**2) +
     >     (A1**2*A3**4)/(A1**2*A3**2 + A1**2*B3**2) +
     >     (A1**2*A3**2*B3**2)/(A1**2*A3**2 + A1**2*B3**2) -
     >     (A1*A3*Sqrt(-(A1**6*B3**2) + 2*A1**4*A2**2*B3**2 -
     >     A1**2*A2**4*B3**2 + 2*A1**4*A3**2*B3**2 +
     >     2*A1**2*A2**2*A3**2*B3**2 - A1**2*A3**4*B3**2 + 2*A1**4*B3**4
     >     + 2*A1**2*A2**2*B3**4 - 2*A1**2*A3**2*B3**4 - A1**2*B3**6))
     >     /(A1**2*A3**2 + A1**2*B3**2))/(A2*B3),(-A3 +
     >     (A1**4*A3)/(2.*(A1**2*A3**2 + A1**2*B3**2)) -
     >     (A1**2*A2**2*A3)/(2.*(A1**2*A3**2 + A1**2*B3**2)) +
     >     (A1**2*A3**3)/(2.*(A1**2*A3**2 + A1**2*B3**2)) +
     >     (A1**2*A3*B3**2)/(2.*(A1**2*A3**2 + A1**2*B3**2)) -
     >     (A1*Sqrt(-(A1**6*B3**2) + 2*A1**4*A2**2*B3**2 -
     >     A1**2*A2**4*B3**2 + 2*A1**4*A3**2*B3**2 +
     >     2*A1**2*A2**2*A3**2*B3**2 - A1**2*A3**4*B3**2 + 2*A1**4*B3**4
     >     + 2*A1**2*A2**2*B3**4 - 2*A1**2*A3**2*B3**4 - A1**2*B3**6))/
     >     (2.*(A1**2*A3**2 + A1**2*B3**2)))/A2)
      endif
      return
      end
*******************************************************************************
      subroutine gettwophot(p1,p2,p3,p4,qph,ng,qind,qind2)
      implicit none
      external tridot
      real*8 tridot,en,eta1,eta2,phi1,phi2
      real*8 qph(40,0:3),qphi(40,0:3)
      real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      real*8 enhard,spd,spdi,c1,c2,c3,c4,enphot,a1,a2,a3,a4
      real*8 b1,b2,b3,b4
      dimension a1(0:3)
      integer k,i,ind,qind,ng,qind2
      external dot
      real*8 dot
        qphi=qph
        spd=0d0

        do k = 1,ng

                q = qphi(k,:)

*********CKKW
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




                if (c1.le.c2.and.c1.le.c3.and.c1.le.c4) spdi =c1
                if (c2.le.c1.and.c2.le.c3.and.c2.le.c4) spdi =c2
                if (c3.le.c2.and.c3.le.c1.and.c3.le.c4) spdi =c3
                if (c4.le.c1.and.c4.le.c2.and.c4.le.c3) spdi =c4
                if(spdi.ge.spd) then
                   spd=spdi
                   qind=k
                endif
             enddo



       qphi(qind,:)=0d0

       spd=0d0

        do k = 1,ng


           if(k.ne.qind) then

              q= qphi(k,:)
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


                if (c1.le.c2.and.c1.le.c3.and.c1.le.c4) spdi =c1
                if (c2.le.c1.and.c2.le.c3.and.c2.le.c4) spdi =c2
                if (c3.le.c2.and.c3.le.c1.and.c3.le.c4) spdi =c3
		if (c4.le.c1.and.c4.le.c2.and.c4.le.c3) spdi =c4

                if(spdi.ge.spd) then
                   spd=spdi
                   qind2=k
                endif
                endif
         enddo

        return
        end
************************************************
*******************************************************************************
      subroutine ckkw_crit(p1,p2,p3,p4,qph,ng,nisr,nfsr)
      implicit none
      external tridot
      real*8 tridot,en,eta1,eta2,phi1,phi2
      real*8 qph(40,0:3),qphi(40,0:3)
      real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      real*8 enhard,spd,spdi,c1,c2,c3,c4,enphot,a1,a2,a3,a4
      real*8 b1,b2,b3,b4
      dimension a1(0:3)
      integer k,i,ind,ng,nisr,nfsr
      external dot
      real*8 dot
        qphi=qph
        spd=0d0
        nisr=0
        nfsr=0
        do k = 1,ng

           q = qphi(k,:)

*********CKKW
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




           if (c1.le.c2.and.c1.le.c3.and.c1.le.c4) nisr=nisr+1
           if (c2.le.c1.and.c2.le.c3.and.c2.le.c4) nisr=nisr+1
           if (c3.le.c2.and.c3.le.c1.and.c3.le.c4) nfsr=nfsr+1
           if (c4.le.c1.and.c4.le.c2.and.c4.le.c3) nfsr=nfsr+1
        enddo



        return
        end
*************************************
      subroutine prop_value(ind,p1,p2,p3,p4,qph,ng,qind)
      implicit none
      external tridot
      real*8 tridot,en,eta1,eta2,phi1,phi2
      real*8 qph(40,0:3),qphi(40,0:3),qind(40,4)
      real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      real*8 enhard,spd,spdi,c1,c2,c3,c4,enphot,a1,a2,a3,a4
      real*8 b1,b2,b3,b4
      dimension a1(0:3)
      integer k,i,ind,ng,nisr,nfsr
      external dot
      real*8 dot
      qphi=qph
      spd=0d0
      qind=0.d0
      nisr=0
      nfsr=0
      if(ind.eq.1) then
         do k = 1,ng

            q = qph(k,:)

            c1 = dot(p1,q)
            c2 = dot(p2,q)
            c3 = dot(p3,q)
            c4 = dot(p4,q)


            if (c1.le.c2.and.c1.le.c3.and.c1.le.c4) qind(k,1)=c1
            if (c2.le.c1.and.c2.le.c3.and.c2.le.c4) qind(k,2)=c2
            if (c3.le.c2.and.c3.le.c1.and.c3.le.c4) qind(k,3)=c3
            if (c4.le.c1.and.c4.le.c2.and.c4.le.c3) qind(k,4)=c4

            if(c1.le.c2) qind(k,1)=c1
            if(c2.le.c1) qind(k,1)=c2

            if(c3.le.c4) qind(k,2)=c3
            if(c4.le.c3) qind(k,2)=c4
         enddo

      elseif(ind.eq.2) then
         do k=1,ng
            q=qph(k,:)

            c1 = acos(tridot(p1,q)/(dsqrt(tridot(p1,p1))*q(0)))
            c2 = acos(tridot(p2,q)/(dsqrt(tridot(p2,p2))*q(0)))
            c3 = acos(tridot(p3,q)/(dsqrt(tridot(p3,p3))*q(0)))
            c4 = acos(tridot(p4,q)/(dsqrt(tridot(p4,p4))*q(0)))

            if (c1.le.c2.and.c1.le.c3.and.c1.le.c4) qind(k,1)=c1
            if (c2.le.c1.and.c2.le.c3.and.c2.le.c4) qind(k,2)=c2
            if (c3.le.c2.and.c3.le.c1.and.c3.le.c4) qind(k,3)=c3
            if (c4.le.c1.and.c4.le.c2.and.c4.le.c3) qind(k,4)=c4

            if(c1.le.c2) qind(k,1)=c1
            if(c2.le.c1) qind(k,1)=c2

            if(c3.le.c4) qind(k,2)=c3
            if(c4.le.c3) qind(k,2)=c4

         enddo
      elseif(ind.eq.3) then
         do k=1,ng
            q = qphi(k,:)

*********CKKW
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

c$$$  if (c1.le.c2.and.c1.le.c3.and.c1.le.c4) qind(k,1)=c1
c$$$  if (c2.le.c1.and.c2.le.c3.and.c2.le.c4) qind(k,2)=c2
c$$$
c$$$  if (c3.le.c2.and.c3.le.c1.and.c3.le.c4) qind(k,3)=c3
c$$$
c$$$  if (c4.le.c1.and.c4.le.c2.and.c4.le.c3) qind(k,4)=c4

            if(c1.le.c2)  qind(k,1)=c1
            if(c2.le.c1)  qind(k,1)=c2
         enddo
      endif


      return
      end
*********
      subroutine check_propagators(p1,p2,qph,ng,p1m,p2m,qph2,error)
      implicit none
      external tridot
      real*8 tridot,en,eta1,eta2,phi1,phi2
      real*8 qph(40,0:3),qph2(40,0:3)
      real*8 p1(0:3),p2(0:3),p1m(0:3),p2m(0:3),q(0:3)
      real*8 enhard,spd,spdtwo,spd1,spd2
      real*8 b1,b2,b3,b4
      common/counterg5/counterg4
      integer counterg4
      integer k,i,ind,qind,ng,qind2,error
      external dot
      real*8 dot
      error=0
      spd=1.d16

      do k=1,ng
         q=qph(k,:)
         spd1=abs(dot(p1-q,p1-q)-dot(p1,p1))
         spd2=abs(dot(p2-q,p2-q)-dot(p2,p2))

         if((spd1.le.spd2).and.spd1.le.spd) spd=spd1
         if((spd2.le.spd1).and.spd2.le.spd) spd=spd2
      enddo

      spdtwo=1.d16
      do k=1,2
         q=qph2(k,:)
         spd1=abs(dot(p1m-q,p1m-q)-dot(p1,p1))
         spd2=abs(dot(p2m-q,p2m-q)-dot(p2,p2))

         if((spd1.le.spd2).and.spd1.le.spdtwo) spdtwo=spd1
	 if((spd2.le.spd1).and.spd2.le.spdtwo) spdtwo=spd2
      enddo
      if(spd/spdtwo.ge.20.d0) then
         error=1
      endif

      return
        end

***********
      subroutine getcmaxmin(cmax,cmin)
! guessed cos_max and min due to beam energy spread...
      implicit double precision (a-h,o-z)
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/ecms/ecms,ecmsnom
      common/beamspreadsigmas/esig1,esig2
      character*2 fs
      common/finalstate/fs

      cmax   = cos(thmin)
      cmin   = cos(thmax)

      if (esig1.lt.1d-12.and.esig2.lt.1d-12) return

      pi = 4.d0*atan(1.d0)

      eb = ecmsnom * 0.5d0

      sig = 3.89d0
      ebmin = min(eb - sig*esig1,eb - sig*esig2)
      ebmax = max(eb + sig*esig1,eb + sig*esig2)

      iwarn = 0
      if (ebmin.lt.0.5d0*eb) iwarn = 1
      if (ebmax.gt.1.5d0*eb) iwarn = 1

      b1 = (ebmax - ebmin)/(ebmax+ebmin)
      b2 = -b1

      g  = 1.d0/sqrt(1.d0-b1*b1)

      c1 = (cmax - b1)/sqrt(1.d0+b1*b1 - 2.d0 * b1 * cmax)
      c2 = (cmax - b2)/sqrt(1.d0+b2*b2 - 2.d0 * b2 * cmax)

      c3 = (cmin - b1)/sqrt(1.d0+b1*b1 - 2.d0 * b1 * cmin)
      c4 = (cmin - b2)/sqrt(1.d0+b2*b2 - 2.d0 * b2 * cmin)

      cmin = min(c1,c2,c3,c4)
      cmax = max(c1,c2,c3,c4)

      cmin = cmin*1.000005d0
      cmax = cmax*1.000005d0

      if (cmax.gt. 1.d0) cmax =  1.d0
      if (cmin.lt.-1.d0) cmin = -1.d0

      if (fs.eq.'ee') then
         if (cmax.gt.0.99999d0.or.cmin.lt.-0.99999d0) then
            do k = 1,10
            print*,'CM angular range too large because of energy spread'
            print*,'The code can not proceed. Stopping here.'
            enddo
            stop
         endif
      endif

      if (iwarn.gt.0) then
         do k = 1,10
            print*,'Energy spread quite large: expect edge effects'
            print*,' Consider generating on a larger angular range'
            print*,' than your volume'
         enddo
      endif

c      print*,cos(thmin),cos(thmax)
c      print*,cmax,cmin
c      print*,acos(cmax)*180.d0/pi,acos(cmax/1.000005d0)*180.d0/pi

      return
      end
****************************************************
      subroutine spreadbeams(s1,s2,ame,am,pb1,pb2,ie)
      implicit double precision (a-h,o-z)
      common/ecms/ecms,ecmsnominal
      dimension pb1(0:3),pb2(0:3)
      common/beamspreadc/ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         ecmsnominal = ecms
         ifirst = 1
      endif

      ie = 0

      ecms = ecmsnominal

      eb1 = ecms * 0.5d0
      eb2 = eb1

      if (s1.gt.1.d-12) eb1 = gau_spread(eb1,s1)
      if (s2.gt.1.d-12) eb2 = gau_spread(eb2,s2)

      if (eb1.lt.ame.or.eb2.lt.ame) then
         ie = 1
         return
      endif

      b1 = (1.d0-ame/eb1)*(1.d0+ame/eb1)
      b1 = sqrt(b1)

      b2 = (1.d0-ame/eb2)*(1.d0+ame/eb2)
      b2 = sqrt(b2)

      pb1(0) =  eb1
      pb1(1) =  0.d0
      pb1(2) =  0.d0
      pb1(3) =  eb1*b1
      pb2(0) =  eb2
      pb2(1) =  0.d0
      pb2(2) =  0.d0
      pb2(3) = -eb2*b2

      ecms = sqrt(dot(pb1+pb2,pb1+pb2))

      if (ecms.lt.2.d0*am) then
         ie = 1
      endif

      return
      end
*********************************************************************
      function gau_spread(e,s)
      implicit double precision (a-h,o-z)
      real*8 csi(2)
      common/gauspreadfun/tpi,x1,x2,ifirst,icalc
      common/gsptest/sum,sum2,ic
      data ifirst /0/
      if (ifirst.eq.0) then
         tpi = 8.d0 * atan(1.d0)
         ifirst = 1
         icalc = 1
         ic = 0
         sum  = 0.d0
         sum2 = 0.d0
      endif

***   using Box-Mueller algorithm to generate 2 normal numbers
      x = 0.d0
      if (icalc.eq.1) then
         call getrnd(csi,2)
         cs1 = csi(1)
         cs2 = csi(2)
         r   = sqrt(-2.d0*log(cs1))
         sf  = sin(tpi*cs2)
         cf  = cos(tpi*cs2)
         x1 = r*sf
         x2 = r*cf

         x     = x1
         icalc = 0
      else
         x     = x2
         icalc = 1
      endif

      gau_spread = e + s*x

c      ic = ic + 1
c      sum  = sum  + gau_spread
c      sum2 = sum2 + gau_spread*gau_spread
c      if (mod(ic,1000000).eq.0)
c     .     print*,sum/ic,sqrt(abs(sum2/ic - (sum/ic)**2))

      return
      end
*********************************************************************
      subroutine printstatus(icode,kl,p1,p2,qph,xs,var,varb,sd,sdm,fm)
      implicit double precision (a-h,o-z)
      integer*8 kl,iwriteout
      character*12 col(20)
      common/colors/col
      dimension p1(0:3),p2(0:3),qph(40,0:3),q(0:3)
      dimension pin1(0:3),pin2(0:3),ptmp(0:3)
      common/ifirstprst/ifirst,icount,procid,statusfile,machine
      data ifirst /1/
      data icount /0/
      character*2 fs
      common/finalstate/fs
      common/ecms/ecms,ecmsnom
      common/nphot_mode/nphotmode
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/zparameters/zm,gz,ve,ae,rv,ra,wm,s2th,gfermi,sqrt2,um
      common/energiabeam/ebeam
      common/parameters/ame,ammu,convfac,alpha,pi
      common/intinput/iwriteout,iseed,nsearch,iverbose
      common/ialpharunning/iarun
      common/idebugging/idebug
      common/momentainitial/pin1,pin2
      common/radpattern/nph(4)

      integer isvec(104)
      common/rlxstatus/isvec

      character*20 statusfile,procid,machine

      if (icount.eq.0) then
         call itoa(getpid(),procid)
         call hostnm(machine)
         statusfile = 'printstatus-'//procid
      endif

      open(543,file=statusfile,status='unknown',access='append')

      icount = icount + 1
      nfot = nph(1)+nph(2)+nph(3)+nph(4)
      amfs = 0.d0
      if (fs.eq.'ee') amfs = ame
      if (fs.eq.'mm'.or.fs.eq.'mr'.or.fs.eq.'mi') amfs = ammu
      write(543,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(543,*)'printing point & event n.: ',icode,', ',kl
      write(543,*)'final state: ',fs
      write(543,'(1x,A,f25.16,A,f25.16)')'cross section: ',xs,' +-',var
      write(543,'(1x,A,f25.16,A,f25.16)')'var now and bfr: ',var,',',
     .     varb
      write(543,'(1x,A,f25.8)')        'diff. xsect: ',sd
      write(543,'(1x,A,f25.8)')        'max xsect:   ',sdm
      write(543,'(1x,A,f25.16)')        'fmax (should be < max. xs): ',
     .     fm
      write(543,*)'n. phot. ',nfot,' and rad. pattern ',nph
c$$$      write(543,'(1x,A,f25.16,f25.16,f25.16,f25.16)')'  p1 ',pin1
c$$$      write(543,'(1x,A,f25.16,f25.16,f25.16,f25.16)')'  p2 ',pin2
c$$$      write(543,'(1x,A,f25.16,f25.16,f25.16,f25.16)')'  p3 ',p1
c$$$      write(543,'(1x,A,f25.16,f25.16,f25.16,f25.16)')'  p4 ',p2
c$$$      if (nfot.gt.0) then
c$$$         do k = 1,nfot
c$$$            q(0) = qph(k,0)
c$$$            q(1) = qph(k,1)
c$$$            q(2) = qph(k,2)
c$$$            q(3) = qph(k,3)
c$$$            write(543,'(1x,i2,A,f25.16,f25.16,f25.16,f25.16)')k,' q ',q
c$$$         enddo
c$$$      endif
      write(543,*)'  p1 ',pin1
      write(543,*)'  p2 ',pin2
      write(543,*)'  p3 ',p1
      write(543,*)'  p4 ',p2
      if (nfot.gt.0) then
         do k = 1,nfot
            q(0) = qph(k,0)
            q(1) = qph(k,1)
            q(2) = qph(k,2)
            q(3) = qph(k,3)
            write(543,*)k,' q ',q
         enddo
      endif

      write(543,'(1x,A,f25.16)')'th. angle p3-p1',angledeg(p1,pin1)
      write(543,'(1x,A,f25.16)')'th. angle p4-p1',angledeg(p2,pin1)
      write(543,'(1x,A,f25.16)')'th. angle p3-p4',angledeg(p1,p2)
      write(543,'(1x,A,f25.16)')'acoll. p3-p4   ',acollinearity(p1,p2)
      if (nfot.gt.0) then
         do k = 1,nfot
            q(0) = qph(k,0)
            q(1) = qph(k,1)
            q(2) = qph(k,2)
            q(3) = qph(k,3)
            write(543,'(1x,i2,A,f25.16)')k,' angle q-p1 ',
     .           angledeg(q,pin1)
            write(543,'(1x,i2,A,f25.16)')k,' angle q-p2 ',
     .           angledeg(q,pin2)
            write(543,'(1x,i2,A,f25.16)')k,' angle q-p3 ',angledeg(q,p1)
            write(543,'(1x,i2,A,f25.16)')k,' angle q-p4 ',angledeg(q,p2)
         enddo
      endif
      write(543,*)'invariants:'
      do k = 0,3
         ptmp(k) = pin1(k)+pin2(k)
      enddo
      s12 = dot(ptmp,ptmp)
      do k = 0,3
         ptmp(k) = p1(k)+p2(k)
      enddo
      s34 = dot(ptmp,ptmp)
      do k = 0,3
         ptmp(k) = pin1(k)-p1(k)
      enddo
      t13 = dot(ptmp,ptmp)
      do k = 0,3
         ptmp(k) = pin2(k)-p2(k)
      enddo
      t24 = dot(ptmp,ptmp)
      do k = 0,3
         ptmp(k) = pin1(k)-p2(k)
      enddo
      u14 = dot(ptmp,ptmp)
      do k = 0,3
         ptmp(k) = pin2(k)-p1(k)
      enddo
      u23 = dot(ptmp,ptmp)
      write(543,'(1x,A,f25.16,f25.16,f25.16)')'s12 s34 t13',s12,s34,t13
      write(543,'(1x,A,f25.16,f25.16,f25.16)')'t24 u14 u23',t24,u14,u23
      if (nfot.gt.0) then
         do k = 1,nfot
            q(0) = qph(k,0)
            q(1) = qph(k,1)
            q(2) = qph(k,2)
            q(3) = qph(k,3)
            write(543,'(1x,i2,A,f25.16)')k,' p1q',dot(pin1,q)
            write(543,'(1x,i2,A,f25.16)')k,' p2q',dot(pin2,q)
            write(543,'(1x,i2,A,f25.16)')k,' p3q',dot(p1,q)
            write(543,'(1x,i2,A,f25.16)')k,' p4q',dot(p2,q)
         enddo
      endif
      write(543,*)'Ranlux status:'
      do kk = 1,104
         write(543,*)'       isvec(',kk,') =',isvec(kk)
      enddo
      write(543,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

      close(543)

      return
      end

      subroutine itoa(int,a)
      implicit integer (h-n)
      implicit character (a-g)
      character*(*) a
      parameter (nbase=10)
      call resetname(a)
      ichar0 = ichar('0')
      idiv = nbase
      i    = int
      ncifre    = 1
      do while((i/idiv).gt.0)
         idiv = idiv * nbase
         ncifre = ncifre+1
      enddo
      idiv = idiv/nbase
      do k =1,ncifre
         icifra = i/idiv
         ia = icifra + ichar0
         a(k:k) = char(ia)
         i = i - icifra*idiv
         idiv = idiv/nbase
      enddo
      return
      end
      subroutine resetname(name)
      integer i
      character*(*) name
      do i = 1,len(name)
         name(i:i) = ' '
      enddo
      return
      end

      subroutine printvector(ip,p1,p2,qph)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),qph(40,0:3),q(0:3)
      common/channelref/iref
      common/ifirstprv/ifirst,icount
      data ifirst /1/
      data icount /0/
      icount = icount + 1
      if (ifirst.eq.1) then
         open(33,file='momenta',status='unknown')
         ifirst = 0
      else
         open(33,file='momenta',status='unknown',access='append')
      endif
c$$$      print*,'set arrow 1 from 0,0,0 to ',p1(1),',',p1(2),',',p1(3)
c$$$      print*,'set arrow 2 from 0,0,0 to ',p2(1),',',p2(2),',',p2(3)
c$$$      print*,'set arrow 3 from 0,0,0 to ',q(1) ,',', q(2),',',q(3)
c$$$      print*,'set'//
c$$$     .' xrange [',min(p1(1),p2(1),q(1)),':',max(p1(1),p2(1),q(1)),'];'//
c$$$     .' set '//
c$$$     .'yrange [',min(p1(2),p2(2),q(2)),':\\',max(p1(2),p2(2),q(2)),'];'
c$$$     .//
c$$$     .' set'//
c$$$     .' zrange [',min(p1(3),p2(3),q(3)),':',max(p1(3),p2(3),q(3)),']'
      write(33,*)'point, count, iref ',ip,icount,iref
c
      write(33,*)'set arrow 1 from 0,0,0 to',p1(3),',',p1(1),',',p1(2),
     .' lt 1'
      write(33,*)'set arrow 2 from 0,0,0 to',p2(3),',',p2(1),',',p2(2),
     .' lt 2'
      n = 0
      do k = 1,40
         if (qph(k,0).gt.0.d0) n = n+1
      enddo
      if (n.gt.0) then
         do k = 3,n+2
            q(0) = qph(k-2,0)
            q(1) = qph(k-2,1)
            q(2) = qph(k-2,2)
            q(3) = qph(k-2,3)
      write(33,*)'set arrow',k
     .,' from 0,0,0 to',q(3),',',q(1),',',q(2),''//
     .' lt 3'
        enddo
      endif

      close(33)

      return
      end

      subroutine randperm(n,irand)
      implicit double precision (a-h,o-z)
      dimension irand(n),ind(n)
      real*8 csi(1)
      do k = 1,n
         ind(k) = k
      enddo
      do j = 0,n-1
         k = n-j
         call getrnd(csi,1)
         ii = k * csi(1) + 1
         if (ii.gt.k) ii = k
         irand(j+1) = ind(ii)
         itmp = ind(k)
         ind(k) = ind(ii)
         ind(ii) = itmp
      enddo
      return
      end
**
      function acollinearity(p1,p2)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      common/acollcommon/ifirst
      data ifirst /0/
      if (ifirst.eq.0) then
c         print*,'CHANGED ACOLLINEARITY DEFINITION!'
c         print*,'CHANGED ACOLLINEARITY DEFINITION!'
c         print*,'CHANGED ACOLLINEARITY DEFINITION!'
         ifirst = 1
      endif
      a1 = acos(p1(3)/sqrt(tridot(p1,p1)))*180.d0/pi
      a2 = acos(p2(3)/sqrt(tridot(p2,p2)))*180.d0/pi

C for me this one is more sensible....
c      a1 = 0.d0
c      a2 = tridot(p1,p2)/sqrt(tridot(p1,p1))/sqrt(tridot(p1,p1))
c      if (a2.lt.-1.d0) a2 = -1.d0
c      if (a2.gt. 1.d0) a2 =  1.d0
c      a2 = acos(a2) * 180.d0/pi

      acollinearity = abs(180.d0 - a1 - a2)
      return
      end

      function acollinearityrad(p1,p2)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      common/acollcommon/ifirst
      data ifirst /0/
      if (ifirst.eq.0) then
c         print*,'CHANGED ACOLLINEARITY DEFINITION!'
c         print*,'CHANGED ACOLLINEARITY DEFINITION!'
c         print*,'CHANGED ACOLLINEARITY DEFINITION!'
         ifirst = 1
      endif

      a1 = acos(p1(3)/sqrt(tridot(p1,p1)))
      a2 = acos(p2(3)/sqrt(tridot(p2,p2)))

C for me this one is more sensible....
c      a1 = 0.d0
c      a2 = tridot(p1,p2)/sqrt(tridot(p1,p1))/sqrt(tridot(p1,p1))
c      if (a2.lt.-1.d0) a2 = -1.d0
c      if (a2.gt. 1.d0) a2 =  1.d0
c      a2 = acos(a2)

      acollinearityrad = abs(pi - a1 - a2)
      return
      end

      function cosine(p1,p2)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3)
      p1modm1 = 1.d0/sqrt(tridot(p1,p1))
      p2modm1 = 1.d0/sqrt(tridot(p2,p2))
      cosine = tridot(p1,p2)*p1modm1*p2modm1
      return
      end

      function angledeg(p1,p2)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      c = cosine(p1,p2)
      if (c.ge. 1.d0) c =  1.d0
      if (c.le.-1.d0) c = -1.d0
      angledeg = acos(c)*180.d0/pi
      return
      end

      function anglerad(p1,p2)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      c = cosine(p1,p2)
      if (c.ge. 1.d0) c =  1.d0
      if (c.le.-1.d0) c = -1.d0
      anglerad = acos(c)
      return
      end

      subroutine rot(idir,vect,pin,pout)
      implicit double precision (a-h,o-z)
      double precision pin(0:3),pout(0:3),pp(0:3),r(3,3),
     >     vers(3),vect(0:3)
* This subroutine rotates the 4-vector pin in the frame where the z-axis is
* directed along the 4-vector vect(0,1,2,3). The rotated vector is stored
* in pout
* idir =  1 ---> direct rotation matrix
* idir = -1 ---> inverse rotation matrix
      pp(0) = pin(0)
      pp(1) = pin(1)
      pp(2) = pin(2)
      pp(3) = pin(3)

      vmo = 1.d0/sqrt(vect(1)**2+vect(2)**2+vect(3)**2)
      vers(1) = vect(1)*vmo
      vers(2) = vect(2)*vmo
      vers(3) = vect(3)*vmo
      vt = sqrt(vers(1)**2+vers(2)**2)

!   BUG - pointed out by CLEO people
!      v1ovt = vers(1)/vt
!      if (vt.eq.0.d0) v1ovt = 0.d0
!      v2ovt = vers(2)/vt
!      if (vt.eq.0.d0) v2ovt = 1.d0

      v1ovt = 0.d0
      v2ovt = 1.d0
      if (vt.gt.0.d0) then
         v1ovt = vers(1)/vt
         v2ovt = vers(2)/vt
      endif

      if (idir.eq.(-1)) then    !! INVERSE rotation matrix
         r(1,1) =  vers(3)*v1ovt
         r(1,2) = -v2ovt
         r(1,3) =  vers(1)
         r(2,1) =  vers(3)*v2ovt
         r(2,2) =  v1ovt
         r(2,3) =  vers(2)
         r(3,1) = -vt
         r(3,2) =  0.d0
         r(3,3) =  vers(3)
      else  ! if (idir.eq.1) !! DIRECT rotation matrix
         r(1,1) =  vers(3)*v1ovt
         r(2,1) = -v2ovt
         r(3,1) =  vers(1)
         r(1,2) =  vers(3)*v2ovt
         r(2,2) =  v1ovt
         r(3,2) =  vers(2)
         r(1,3) = -vt
         r(2,3) =  0.d0
         r(3,3) =  vers(3)
      endif
      pout(0) = pp(0)
      pout(1) = r(1,1)*pp(1) + r(1,2)*pp(2) + r(1,3)*pp(3)
      pout(2) = r(2,1)*pp(1) + r(2,2)*pp(2) + r(2,3)*pp(3)
      pout(3) = r(3,1)*pp(1) + r(3,2)*pp(2) + r(3,3)*pp(3)
      return
      end
*
      subroutine getphi(rnd,ph,w)
      implicit double precision (a-h,o-z)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      real*8 rnd
      ph = 2.d0*pi*rnd
      w  = 2.d0*pi
      return
      end
*
      function getphiangle(ppp)
      implicit double precision (a-h,o-z)
      dimension ppp(0:3)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      pm = sqrt(ppp(1)**2+ppp(2)**2+ppp(3)**2)
      c  = ppp(3)/pm
      s  = sqrt(1.d0-c**2)
      if (s.eq.0.d0) then
         getphiangle = 0.d0
         return
      else
         arg = ppp(1)/pm/s
*  avoiding numerical problems......
         if (abs(arg).ge.1.d0) then
            iarg = arg
            arg  = iarg
         endif
         if (ppp(2).ge.0.d0) getphiangle = acos(arg)
         if (ppp(2).lt.0.d0) getphiangle = 2.d0*pi-acos(arg)
      endif
      return
      end
*
      function factorial_OLD(n)
      integer n,nl
      double precision factorial_OLD
      nl = n
      factorial_OLD = 1.d0
      do while(nl.gt.0)
         factorial_OLD = 1.d0*nl * factorial_OLD
         nl = nl - 1
      enddo
      return
      end
************************************************
      function factorial(n)
      integer n,nl      
      double precision factorial
      double precision factarray(0:15)
      common/factarray/factarray
      data factarray /
     .     1d0,             !  0
     .     1d0,             !  1
     .     2d0,             !  2
     .     6d0,             !  3
     .     24d0,            !  4
     .     120d0,           !  5
     .     720d0,           !  6
     .     5040d0,          !  7
     .     40320d0,         !  8
     .     362880d0,        !  9
     .     3628800d0,       ! 10
     .     39916800d0,      ! 11
     .     479001600d0,     ! 12
     .     6227020800d0,    ! 13
     .     87178291200d0,   ! 14
     .     1307674368000d0/ ! 15
      if (n.le.15) then
         factorial = factarray(n)
         return
      endif      
      nl = n
      factorial = 1.d0
      do while(nl.gt.0)
         factorial = dble(nl) * factorial
         nl = nl - 1
      enddo
      return
      end      
***********
      function eikonalintegral(p1,p2,p3,p4,ch,masses)
      implicit double precision (a-h,l,m,o-z)
      parameter (npart = 4)
      dimension eta(npart),p(0:3),q(0:3),pmat(npart,0:3)
      dimension masses(npart),ch(npart)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
! ch are the charges of the field (not anti-field!!)
! Use this convention: the integer factor in front of the charge
! must be
! -1 --> for incoming particle
! -1 --> for outgoing anti-particle
! +1 --> for outgoing particle
! +1 --> for incoming anti-particle
      eta(1) = -1.d0 * ch(1)
      eta(2) =  1.d0 * ch(2)
      eta(3) =  1.d0 * ch(3)
      eta(4) = -1.d0 * ch(4)
*
      do k = 0,3
         pmat(1,k) = p1(k)
         pmat(2,k) = p2(k)
         pmat(3,k) = p3(k)
         pmat(4,k) = p4(k)
      enddo

! off-diagonal contributions
      softint = 0.d0
      do i = 1,npart-1
         do j = i+1,npart
            etaij = eta(i) * eta(j)
            if (abs(etaij).gt.1.d-3) then
               call rescale_momenta(npart,i,j,pmat,masses,p,q)
               q2 = dot(q,q)
               vl = dot(p,p) - q2
               vl = vl/2.d0
               v  = vl/(p(0) - q(0))
               arglog = 1.d0+2.d0*vl/q2
               terminfra = log(arglog)

               tot = terminfra ! + termfinite(p,q,v)
               tot = - tot *2.d0*dot(p,q)/vl* etaij
               tot = tot * 2.d0 ! this is the double product when
                                ! squaring the eikonal
               softint  = softint + tot
            endif
         enddo
      enddo
      soft_integral = softint
      softint = -4.d0 * (eta(1)**2+eta(2)**2+eta(3)**2+eta(4)**2)
      soft_integral = soft_integral + softint
      soft_integral = soft_integral/16.d0
      eikonalintegral = soft_integral
      return
      end
*     *

*******
      subroutine rescale_momenta(npart,i,j,pmat,masses,p,q)
      implicit double precision (a-h,m,o-z)
      dimension p(0:3),q(0:3),pmat(npart,0:3),masses(npart)
      dimension p1(0:3),p2(0:3)
      do k = 0,3
         p1(k) = pmat(i,k)
         p2(k) = pmat(j,k)
      enddo
      m12 = masses(i)**2
      m22 = masses(j)**2
      p1p2 = dot(p1,p2)
      rho1 = p1p2 + sqrt(p1p2**2 - m12*m22)
      rho1 = rho1 / m12
      rho2 = p1p2 - sqrt(p1p2**2 - m12*m22)
      rho2 = rho2 / m12
      if ( (rho1*p1(0)-p2(0)) .gt. 0.d0) rho = rho1
      if ( (rho2*p1(0)-p2(0)) .gt. 0.d0) rho = rho2
      do k = 0,3
         p(k) = rho * p1(k)
         q(k) = p2(k)
      enddo
      return
      end
***
      function termfinite(p,q,v)
      implicit double precision (a-h,o-z)
      dimension p(0:3),q(0:3)

      u0   = p(0)
      umod = sqrt(tridot(p,p))
      pp = u0 + umod
      pm = u0 - umod

      u0   = q(0)
      umod = sqrt(tridot(q,q))
      qp = u0 + umod
      qm = u0 - umod

! dilog from TOPAZ0
      arg1 = (v-pm)/v
      arg2 = (v-pp)/v
      arg3 = (v-qm)/v
      arg4 = (v-qp)/v

      ddlog1 = ddiloghere(arg1)
      ddlog2 = ddiloghere(arg2)
      ddlog3 = ddiloghere(arg3)
      ddlog3 = ddiloghere(arg4)

c      call tspence(arg1,0.d0,1.d0-arg1,ddlog1,ddim)
c      call tspence(arg2,0.d0,1.d0-arg2,ddlog2,ddim)
c      call tspence(arg3,0.d0,1.d0-arg3,ddlog3,ddim)
c      call tspence(arg4,0.d0,1.d0-arg4,ddlog4,ddim)

      termfin = log(pp/pm)**2/4.d0 + ddlog1 + ddlog2
     >        - log(qp/qm)**2/4.d0 - ddlog3 - ddlog4

      termfinite = termfin
      return
      end
*
      double precision function lambda(x,y,z)
      double precision x,y,z
      lambda = x**2+y**2+z**2 - 2.d0*x*y - 2.d0*x*z - 2.d0*y*z
      return
      end
*
      subroutine exchange_mom(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3),p4(0:3)
      integer i
      do i=0,3
         p4(i)=p1(i)
      enddo
      do i=0,3
         p1(i)=p2(i)
      enddo
      do i=0,3
         p2(i)=p4(i)
      enddo
      return
      end
*
***
      subroutine new_boost(p,q,qq,idir)
      implicit double precision (a-h,o-z)
      integer idir
      double precision p(0:3),q(0:3),qq(0:3),vboost(0:3)
      double precision b,g,p0mu,m_sq
! idir =  1: q is boosted where p is at rest and put in qq
! idir = -1: q, in the frame where p is at rest, is boosted where p is not
!            at rest. The boosted q is put in qq.

      m_sq = p(0)**2 - dot_product(p(1:3), p(1:3))

      if (m_sq .le. 0.d0) then
          qq = q
          return
      endif

      g = p(0)/sqrt(m_sq) ! gamma of the boost

      if (g.eq.1.d0) then
         qq = q
         return
      endif

      p0mu = 1.d0/p(0)
      vboost(0) = 0.d0

      if (idir.eq.1) then
         vboost(1:3) = p(1:3)*p0mu
      else
         vboost(1:3) = -p(1:3)*p0mu
      endif

      call boost(g,vboost,q,qq)
      return
      end
**
      subroutine boost(g,v,p,q)
      implicit double precision (a-h,o-z)
      double precision p(0:3),q(0:3),pp(0:3),v(0:3)
! this is essential if boost is called with the same momenta as arguments!!!!
      pp   = p
      ppdv = tridot(pp,v)
      v2   = v(1)**2+v(2)**2+v(3)**2
*
      q(0)   = g*(pp(0)-ppdv)
c     x      = (g - 1.d0)*ppdv/v2 - g*pp(0)
      x      = g*( g/(1.d0+g) *ppdv - pp(0))
      q(1:3) = pp(1:3) + x*v(1:3)
      return
      end
**************************************************
      function sdot(p,q)
      double precision  p(0:3),q(0:3)
      double precision sdot
      sdot=p(0)*q(0)+p(1)*q(1)+p(2)*q(2)+p(3)*q(3)
      return
      end
**
      function dot(p,q)
      implicit double precision (a-h,o-z)
      double precision p(0:3),q(0:3)
      dot=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
      return
      end

*****
      function dotq(p,q)
      real*16 dotq
      double precision p(0:3),q(0:3)
      dotq=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
      return
      end
*************
      function mdot(p,q)
      implicit double precision (a-h,o-z)
      double precision p(4),q(4)
      mdot=-(p(1)*q(1)+p(2)*q(2)+p(3)*q(3))+p(4)*q(4)
      return
      end
*************
	subroutine conv_mauro(pin1cl,pin2cl,p1cl,p2cl,q,pin,pout)

      implicit double precision (a-h,o-z)
      dimension pin1cl(0:3),pin2cl(0:3),p1cl(0:3),p2cl(0:3),q(0:3)
      dimension pin1cl1(4),pin2cl1(4),p1cl1(4),p2cl1(4),qcl(4)

      dimension pclin(4,2),pclout(4,3),pin(4,2),pout(4,3)


        do ki = 1,3
            pin1cl1(ki) = pin1cl(ki)
	    pin2cl1(ki) = pin2cl(ki)
            p1cl1(ki) = p1cl(ki)
            p2cl1(ki) = p2cl(ki)
            qcl(ki) =q(ki)
        enddo
	pin1cl1(4) = pin1cl(0)
	pin2cl1(4) = pin2cl(0)
	p1cl1(4) = p1cl(0)
	p2cl1(4) = p2cl(0)
	qcl(4) =q(0)

	pin(:,1)=pin1cl1
	pin(:,2)=pin2cl1
	pout(:,1)=p1cl1
	pout(:,2)=p2cl1
	pout(:,3)=qcl
     	 end

****************************
	subroutine conv_mom(pin,pout,pin1,pin2,p1,p2,q)

      implicit double precision (a-h,o-z)
      dimension pin(4,2),pout(4,3)
      dimension pin1(0:3),pin2(0:3),p1(0:3),p2(0:3),q(0:3)

        do ki = 1,3
            pin1(ki) = pin(ki,1)
            pin2(ki) = pin(ki,2)
	    p1(ki)=pout(ki,1)
	    p2(ki)=pout(ki,2)
	    q(ki)=pout(ki,3)
        enddo
	pin1(0) = pin(4,1)
	pin2(0) = pin(4,2)
	p1(0)=pout(4,1)
	p2(0)=pout(4,2)
	q(0)=pout(4,3)

     	 end
******************

      subroutine conv_mom_zero(pin,pout,pin1,pin2,p1,p2)

      implicit double precision (a-h,o-z)
      dimension pin(4,2),pout(4,2)
      dimension pin1(0:3),pin2(0:3),p1(0:3),p2(0:3)

        do ki = 1,3
            pin1(ki) = pin(ki,1)
            pin2(ki) = pin(ki,2)
	    p1(ki)=pout(ki,1)
	    p2(ki)=pout(ki,2)
        enddo
	pin1(0) = pin(4,1)
	pin2(0) = pin(4,2)
	p1(0)=pout(4,1)
	p2(0)=pout(4,2)
     	 end
****************************
	subroutine conv_mom_mauro(pin,pout)

        implicit double precision (a-h,o-z)
        dimension pin(0:3),pout(4)

        do ki = 1,3
            pout(ki) = pin(ki)
        enddo
	pout(4) = pin(0)


     	 end
*************
	subroutine conv_mom_mauro1(pout,pin)

      implicit double precision (a-h,o-z)
      dimension pin(0:3),pout(4)

        do ki = 1,3
            pin(ki) = pout(ki)
        enddo
	pin(0) = pout(4)


     	 end
***************
	subroutine CMUL2(x, A,tmp)
	implicit none
	double precision x(0:3)
	double precision :: A(0:3,0:3)
	double precision :: tmp(0:3)
	INTEGER :: i, j,tip,etai

	    do i = 0, 3
	      tmp(i) = 0d0
	      do j = 0, 3
		if(j.eq.0)then
		  etai=1
		else
		  etai=-1
		endif
		tmp(i) = tmp(i) + etai*A(i,j) *x(j)
	      enddo
	    enddo
	 return
	end



********

*****

      subroutine CMUL2_c(x, A,tmp)
	implicit none
	double precision x(0:3)
	double complex :: A(0:3,0:3)
	double complex :: tmp(0:3)
	INTEGER :: i, j,tip,etai
c QUESTO È PER IMPLEMENTARE IL PRODOTTO SCALARE TRA UNA MATRICE E UN VETTORE, QUINDI IL CASO E^MUNU PMU PNU
	    do i = 0, 3
	      tmp(i) = 0d0
	      do j = 0, 3

		if(j.eq.0)then
		  etai=1
		else
		  etai=-1
		endif
		tmp(i) = tmp(i) + etai*A(i,j)*x(j)
	      enddo
	    enddo
	return
      end


**
      subroutine rank2(E,p,q,tip,sp)
      double precision p(0:3),q(0:3)
      double precision E(0:9)
      double precision :: A(0:3, 0:3)
      double precision :: tmp(0:3)
      double precision  sp
      integer tip

      external dot
      double precision  dot


      A(0,0)=E(0)
      A(0,1)=E(1)
      A(1,0)=E(1)
      A(2,0)=E(2)
      A(0,2)=E(2)
      A(3,0)=E(3)
      A(0,3)=E(3)
      A(1,1)=E(4)
      A(1,2)=E(5)
      A(2,1)=E(5)
      A(1,3)=E(6)
      A(3,1)=E(6)
      A(2,2)=E(7)
      A(2,3)=E(8)
      A(3,2)=E(8)
      A(3,3)=E(9)
c qui implemento il prodotto scalare tra matrice e^munu e due vettori
	if(tip.eq.1) then
	      call CMUL2(q,A, tmp)
!	      print*,tmp
	      sp=dot(tmp,p)
	elseif(tip.eq.2) then
	      sp=A(0,0)-A(1,1)-A(2,2)-A(3,3)
	endif
	return
      end


***************
      subroutine rank2_c(E,p,q,tip,sp)
      double precision p(0:3),q(0:3)
      double complex E(0:9)
      double complex A(0:3,0:3)
      double complex tmp(0:3)
      double complex  sp
      integer tip
      double precision b,c

      external dot
      double precision  dot


      A(0,0)=E(0)
      A(0,1)=E(1)
      A(1,0)=E(1)
      A(2,0)=E(2)
      A(0,2)=E(2)
      A(3,0)=E(3)
      A(0,3)=E(3)
      A(1,1)=E(4)
      A(1,2)=E(5)
      A(2,1)=E(5)
      A(1,3)=E(6)
      A(3,1)=E(6)
      A(2,2)=E(7)
      A(2,3)=E(8)
      A(3,2)=E(8)
      A(3,3)=E(9)
c qui implemento il prodotto scalare tra matrice e^munu e due vettori
	if(tip.eq.1) then
	      call CMUL2_c(q,A, tmp)
!	      print*,tmp
	      sp=dot(REAL(tmp),p)+(0.d0,1.d0)*dot(AIMAG(tmp),p)
	elseif(tip.eq.2) then
	      sp=A(0,0)-A(1,1)-A(2,2)-A(3,3)
	endif
	return
      end
***************

*****

	  subroutine CMUL3_c(x,y,z,A,sp)
	    implicit none
	    double precision x(0:3),y(0:3),z(0:3)
	    double complex :: A(0:3,0:3,0:3)
	    double complex :: sp
	    INTEGER :: i, j,k,tip,etaj,etak,etai
c QEUSTO È PER IMPLEMENTARE IL PRODOTTO SCALARE TRA UNA MATRICE E UN VETTORE, QUINDI IL CASO E^MUNU PMU PNU
		    sp=0d0;
		    do i = 0, 3
		      do j = 0, 3
		        do k=0,3
				if(i.eq.0)then
				  etai=1
				else
				  etai=-1
				endif
				if(j.eq.0)then
				  etaj=1
				else
				  etaj=-1
				endif

				if(k.eq.0)then
				  etak=1
				else
				  etak=-1
				endif
				sp = sp + etak*etaj*etai*A(i,j,k)
     >				      *x(i)*y(j)*z(k)
			enddo
		      enddo
		    enddo
	    return
	  end

**********
      subroutine CMUL3(x,y,z,A,sp)
            implicit none
            double precision x(0:3),y(0:3),z(0:3)
            double precision :: A(0:3,0:3,0:3)
            double precision :: sp
            INTEGER :: i, j,k,tip,etaj,etak,etai
c QEUSTO È PER IMPLEMENTARE IL PRODOTTO SCALARE TRA UNA MATRICE E UN VETTORE, QUINDI IL CASO E^MUNU PMU PNU
                    sp=0d0;
                    do i = 0, 3
                      do j = 0, 3
                        do k=0,3
                                if(i.eq.0)then
                                  etai=1
                                else
                                  etai=-1
                                endif
                                if(j.eq.0)then
                                  etaj=1
                                else
                                  etaj=-1
                                endif

				if(k.eq.0)then
                                  etak=1
                                else
                                  etak=-1
                                endif
                                sp = sp + etak*etaj*etai*A(i,j,k)
     >                                *x(i)*y(j)*z(k)
                        enddo
                      enddo
                    enddo
            return
          end


********

      subroutine rank3(E,p1,p2,p3,tip,sp)
      double precision p1(0:3),p2(0:3),p3(0:3)
      double precision E(0:19)
      double precision :: A(0:3, 0:3,0:3)
      double precision  sp
      integer i,j,k,etai,etaj
      integer tip


      external dot
      double precision  dot

	A(0, 0, 0) = E(0)
	A(0, 0, 1) = E(1)
	A(0, 1, 0) = E(1)
	A(1, 0, 0) = E(1)
	A(0, 0, 2) = E(2)
	A(0, 2, 0) = E(2)
	A(2, 0, 0) = E(2)
	A(0, 0, 3) = E(3)
	A(0, 3, 0) = E(3)
	A(3, 0, 0) = E(3)
	A(0, 1, 1) = E(4)
	A(1, 0, 1) = E(4)
	A(1, 1, 0) = E(4)
	A(0, 1, 2) = E(5)
	A(0, 2, 1) = E(5)
	A(1, 0, 2) = E(5)
	A(1, 2, 0) = E(5)
	A(2, 0, 1) = E(5)
	A(2, 1, 0) = E(5)
	A(0, 1, 3) = E(6)
	A(0, 3, 1) = E(6)
	A(1, 0, 3) = E(6)
	A(1, 3, 0) = E(6)
	A(3, 0, 1) = E(6)
	A(3, 1, 0) = E(6)
	A(0, 2, 2) = E(7)
	A(2, 0, 2) = E(7)
	A(2, 2, 0) = E(7)
	A(0, 2, 3) = E(8)
	A(0, 3, 2) = E(8)
	A(2, 0, 3) = E(8)
	A(2, 3, 0) = E(8)
	A(3, 0, 2) = E(8)
	A(3, 2, 0) = E(8)
	A(0, 3, 3) = E(9)
	A(3, 0, 3) = E(9)
	A(3, 3, 0) = E(9)
	A(1, 1, 1) = E(10)
	A(1, 1, 2) = E(11)
	A(1, 2, 1) = E(11)
	A(2, 1, 1) = E(11)
	A(1, 1, 3) = E(12)
	A(1, 3, 1) = E(12)
	A(3, 1, 1) = E(12)
	A(1, 2, 2) = E(13)
	A(2, 1, 2) = E(13)
	A(2, 2, 1) = E(13)
	A(1, 2, 3) = E(14)
	A(1, 3, 2) = E(14)
	A(2, 1, 3) = E(14)
	A(2, 3, 1) = E(14)
	A(3, 1, 2) = E(14)
	A(3, 2, 1) = E(14)
	A(1, 3, 3) = E(15)
	A(3, 1, 3) = E(15)
	A(3, 3, 1) = E(15)
	A(2, 2, 2) = E(16)
	A(2, 2, 3) = E(17)
	A(2, 3, 2) = E(17)
	A(3, 2, 2) = E(17)
	A(2, 3, 3) = E(18)
	A(3, 2, 3) = E(18)
	A(3, 3, 2) = E(18)
	A(3, 3, 3) = E(19)

	if(tip.eq.1) then
	    call CMUL3(p1,p2,p3,A,sp)
	elseif(tip.eq.2) then
	    sp=0d0;
	    do i=0,3
	     do j=0,3
	      if(j.eq.0) then
	      	etaj=1
	      else
	      	etaj=-1
	      endif
	      if(i.eq.0) then
	      	etai=1
	      else
	      	etai=-1
	      endif
	      sp=sp+A(i,i,j)*etai*etaj*p3(j)
	     enddo
	    enddo
	endif
	return
      end
*********


      subroutine rank3_c(E,p1,p2,p3,tip,sp)
      double precision p1(0:3),p2(0:3),p3(0:3)
      double complex E(0:19)
      double complex :: A(0:3, 0:3,0:3)
      double complex  sp
      integer i,j,k,etai,etaj
      integer tip


      external dot
      double precision  dot

	A(0, 0, 0) = E(0)
	A(0, 0, 1) = E(1)
	A(0, 1, 0) = E(1)
	A(1, 0, 0) = E(1)
	A(0, 0, 2) = E(2)
	A(0, 2, 0) = E(2)
	A(2, 0, 0) = E(2)
	A(0, 0, 3) = E(3)
	A(0, 3, 0) = E(3)
	A(3, 0, 0) = E(3)
	A(0, 1, 1) = E(4)
	A(1, 0, 1) = E(4)
	A(1, 1, 0) = E(4)
	A(0, 1, 2) = E(5)
	A(0, 2, 1) = E(5)
	A(1, 0, 2) = E(5)
	A(1, 2, 0) = E(5)
	A(2, 0, 1) = E(5)
	A(2, 1, 0) = E(5)
	A(0, 1, 3) = E(6)
	A(0, 3, 1) = E(6)
	A(1, 0, 3) = E(6)
	A(1, 3, 0) = E(6)
	A(3, 0, 1) = E(6)
	A(3, 1, 0) = E(6)
	A(0, 2, 2) = E(7)
	A(2, 0, 2) = E(7)
	A(2, 2, 0) = E(7)
	A(0, 2, 3) = E(8)
	A(0, 3, 2) = E(8)
	A(2, 0, 3) = E(8)
	A(2, 3, 0) = E(8)
	A(3, 0, 2) = E(8)
	A(3, 2, 0) = E(8)
	A(0, 3, 3) = E(9)
	A(3, 0, 3) = E(9)
	A(3, 3, 0) = E(9)
	A(1, 1, 1) = E(10)
	A(1, 1, 2) = E(11)
	A(1, 2, 1) = E(11)
	A(2, 1, 1) = E(11)
	A(1, 1, 3) = E(12)
	A(1, 3, 1) = E(12)
	A(3, 1, 1) = E(12)
	A(1, 2, 2) = E(13)
	A(2, 1, 2) = E(13)
	A(2, 2, 1) = E(13)
	A(1, 2, 3) = E(14)
	A(1, 3, 2) = E(14)
	A(2, 1, 3) = E(14)
	A(2, 3, 1) = E(14)
	A(3, 1, 2) = E(14)
	A(3, 2, 1) = E(14)
	A(1, 3, 3) = E(15)
	A(3, 1, 3) = E(15)
	A(3, 3, 1) = E(15)
	A(2, 2, 2) = E(16)
	A(2, 2, 3) = E(17)
	A(2, 3, 2) = E(17)
	A(3, 2, 2) = E(17)
	A(2, 3, 3) = E(18)
	A(3, 2, 3) = E(18)
	A(3, 3, 2) = E(18)
	A(3, 3, 3) = E(19)

	if(tip.eq.1) then
	    call CMUL3_c(p1,p2,p3,A,sp)
	elseif(tip.eq.2) then
	    sp=0d0;
	    do i=0,3
	     do j=0,3
	      if(j.eq.0) then
	      	etaj=1
	      else
	      	etaj=-1
	      endif
	      if(i.eq.0) then
	      	etai=1
	      else
	      	etai=-1
	      endif
	      sp=sp+A(i,i,j)*etai*etaj*p3(j)
	     enddo
	    enddo
	endif
	return
      end

*     *
      subroutine find_closest_two(values,average)
      implicit none
      double precision values(5)
      double precision closest(2)
      double precision :: diff(5, 5)
      double precision :: min_diff
      double precision average
      integer :: i, j, k, idx1, idx2, idx3

! Calculate pairwise differences and store in the diff array
      do i = 1, 5
         do j = i+1, 5
            if (i.ne.j) then
               diff(i, j) = abs(values(i) - values(j))
            else
               diff(i, j) = 1.0D10 ! Set diagonal to a large value to exclude self-comparison
            endif
         enddo
      enddo
!     Find the smallest two differences
      min_diff = 1.0D30
      do i = 1, 5
         do j = i + 1, 5        ! Only compare once for each pair
            if (diff(i, j) < min_diff) then
               min_diff = diff(i, j)
               idx1 = i
               idx2 = j
            endif
         enddo
      enddo

! Store the closest two values
      closest(1) = values(idx1)
      closest(2) = values(idx2)


      average= (closest(1)+closest(2))/2d0

      end

!
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
	double complex ffpiis
	double precision ffpiis2
	double complex ffpifs
	double precision ffpifs2

**************

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

******ALPHA RUNNING
	ffpiis=ffpiis*vpolc(2.d0*(ampi**2+p3p4))
	ffpiis2=ffpiis*CONJG(ffpiis)

	ffpifs=ffpifs*vpolc(2.d0*(ame**2+p1p2))
	ffpifs2=ffpifs*CONJG(ffpifs)
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
	elmati=real(elmat1*CONJG(ffpiis)*ffpifs)
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/FSR/Pent_Routing1_fsr1.F'
	elmatf=elmat1*ffpifs2
	  elseif(j.eq.3) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing1_isr3.F'
	elmati=elmat1*ffpiis2
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing1_fsr3.F'
	elmatf=real(elmat1*(ffpiis*CONJG(ffpifs)))

	  endif
	elseif(i.eq.2) then
	  if(j.eq.1) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/ISR/Pent_Routing2_isr1.F'
	elmati=real(elmat1*CONJG(ffpiis)*ffpifs)
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/FSR/Pent_Routing2_fsr1.F'
	elmatf=elmat1*ffpifs2
	  elseif(j.eq.3) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing2_isr3.F'
	elmati=elmat1*ffpiis2
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing2_fsr3.F'
	elmatf=real(elmat1*(ffpiis*CONJG(ffpifs)))
	  endif
	elseif(i.eq.3) then
	  if(j.eq.1) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/ISR/Pent_Routing3_isr1.F'
	elmati=real(elmat1*CONJG(ffpiis)*ffpifs)
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/FSR/Pent_Routing3_fsr1.F'
	elmatf=elmat1*ffpifs2
	  elseif(j.eq.3) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing3_isr3.F'
	elmati=elmat1*ffpiis2
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing3_fsr3.F'
	elmatf=real(elmat1*(ffpiis*CONJG(ffpifs)))
	  endif
	elseif(i.eq.4) then
	  if(j.eq.1) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/ISR/Pent_Routing4_isr1.F'
	elmati=real(elmat1*CONJG(ffpiis)*ffpifs)
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/FSR/Pent_Routing4_fsr1.F'
	elmatf=elmat1*ffpifs2
	  elseif(j.eq.3) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing4_isr3.F'
	elmati=elmat1*ffpiis2
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing4_fsr3.F'
	elmatf=real(elmat1*(ffpiis*CONJG(ffpifs)))
	  endif
	elseif(i.eq.5) then
	  if(j.eq.1) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/ISR/Pent_Routing5_isr1.F'
	elmati=real(elmat1*CONJG(ffpiis)*ffpifs)
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_1/FSR/Pent_Routing5_fsr1.F'
	elmatf=elmat1*ffpifs2
	  elseif(j.eq.3) then
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing5_isr3.F'
	elmati=elmat1*ffpiis2
  include 'pipirad/Virtual/IFI_PENT/DIAGRAM_3/Pent_Routing5_fsr3.F'
	elmatf=real(elmat1*(ffpiis*CONJG(ffpifs)))
	  endif
       endif
      end


*********************
      function tridot(p,q)
      implicit double precision (a-h,o-z)
      double precision p(0:3),q(0:3)
      tridot=p(1)*q(1)+p(2)*q(2)+p(3)*q(3)
      return
      end
*
*=========== DILOG FROM TOPAZ0 =====================================
*-----SPENCE--------------------------------------------------------
*     COMPUTES  LI_2(X). ACCURACY IS ABOUT 16 DIGITS
      SUBROUTINE TSPENCE(XR,XI,OMXR,CLI2R,CLI2I)
! from TOPAZ0
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      COMMON/TQPARAM/QPI,QPIS,QEPS,QDELTA
*
      DIMENSION B(0:14),BF(0:14)
      DIMENSION CLNX(2),CLNOMX(2),CLNOY(2),CLNZ(2),CLNOMZ(2)
      DIMENSION ADD1(2),ADD2(2),ADD3(2),PAR(2),RES(2),CT(15),SN(15)
      DIMENSION X(2),OMX(2),Y(2),OY(2),OMY(2),Z(2),OMZ(2),T(2),OMT(2)
      common/storedtspence/b,bf,ifirst
      data ifirst /0/
c CARLO
      if (ifirst.eq.0) then
         QDELTA= 9.025809333D0
         QEPS= 1.D-30
         QPI= 3.141592653589793238462643D0
         QPIS= QPI*QPI
      endif
c CARLO
      X(1)= XR
      X(2)= XI
      OMX(1)= OMXR
      OMX(2)= -XI
      IF(XR.LT.0.D0) THEN
          Y(1)= OMXR
          Y(2)= -XI
          SIGN1= -1.D0
          CALL TCQLNX(X,CLNX)
          CALL TCQLNOMX(X,OMX,CLNOMX)
          ADD1(1)= QPIS/6.D0-CLNX(1)*CLNOMX(1)+CLNX(2)*CLNOMX(2)
          ADD1(2)= -CLNX(1)*CLNOMX(2)-CLNX(2)*CLNOMX(1)
      ELSE
          Y(1)= X(1)
          Y(2)= X(2)
          SIGN1= 1.D0
          ADD1(1)= 0.D0
          ADD1(2)= 0.D0
      ENDIF
      OMY(1)= 1.D0-Y(1)
      OMY(2)= -Y(2)
      YM2= Y(1)*Y(1)+Y(2)*Y(2)
      YM= SQRT(YM2)
      IF(YM.GT.1.D0) THEN
          Z(1)= Y(1)/YM2
          Z(2)= -Y(2)/YM2
          SIGN2= -1.D0
          OY(1)= -Y(1)
          OY(2)= -Y(2)
          CALL TCQLNX(OY,CLNOY)
          ADD2(1)= -QPIS/6.D0-0.5D0*((CLNOY(1))**2-(CLNOY(2))**2)
          ADD2(2)= -CLNOY(1)*CLNOY(2)
      ELSE
          Z(1)= Y(1)
          Z(2)= Y(2)
          SIGN2= 1.D0
          ADD2(1)= 0.D0
          ADD2(2)= 0.D0
      ENDIF
      OMZ(1)= 1.D0-Z(1)
      OMZ(2)= -Z(2)
      ZR= Z(1)
      IF(ZR.GT.0.5D0) THEN
          T(1)= 1.D0-Z(1)
          T(2)= -Z(2)
          OMT(1)= 1.D0-T(1)
          OMT(2)= -T(2)
          SIGN3= -1.D0
          CALL TCQLNX(Z,CLNZ)
          CALL TCQLNOMX(Z,OMZ,CLNOMZ)
          ADD3(1)= QPIS/6.D0-CLNZ(1)*CLNOMZ(1)+CLNZ(2)*CLNOMZ(2)
          ADD3(2)= -CLNZ(1)*CLNOMZ(2)-CLNZ(2)*CLNOMZ(1)
      ELSE
          T(1)= Z(1)
          T(2)= Z(2)
          OMT(1)= 1.D0-T(1)
          OMT(2)= -T(2)
          SIGN3= 1.D0
          ADD3(1)= 0.D0
          ADD3(2)= 0.D0
      ENDIF
      CALL TCQLNOMX(T,OMT,PAR)

      if (ifirst.eq.0) then
         B(0)= 1.D0
         B(1)= -1.D0/2.D0
         B(2)= 1.D0/6.D0
         B(4)= -1.D0/30.D0
         B(6)= 1.D0/42.D0
         B(8)= -1.D0/30.D0
         B(10)= 5.D0/66.D0
         B(12)= -691.D0/2730.D0
         B(14)= 7.D0/6.D0
         FACT= 1.D0
         DO N=0,14
            BF(N)= B(N)/FACT
            FACT= FACT*(N+2.D0)
         ENDDO
      endif
      PARR= PAR(1)
      PARI= PAR(2)
      PARM2= PARR*PARR+PARI*PARI
      PARM= SQRT(PARM2)
      CT(1)= PARR/PARM
      SN(1)= PARI/PARM
      DO N=2,15
          CT(N)= CT(1)*CT(N-1)-SN(1)*SN(N-1)
          SN(N)= SN(1)*CT(N-1)+CT(1)*SN(N-1)
      ENDDO
*
      RES(1)= -((((((((BF(14)*CT(15)*PARM2+BF(12)*CT(13))*PARM2+
     #                 BF(10)*CT(11))*PARM2+BF(8)*CT(9))*PARM2+
     #                 BF(6)*CT(7))*PARM2+BF(4)*CT(5))*PARM2+
     #                 BF(2)*CT(3))*(-PARM)+BF(1)*CT(2))*(-PARM)+
     #                 BF(0)*CT(1))*PARM
      RES(2)= -((((((((BF(14)*SN(15)*PARM2+BF(12)*SN(13))*PARM2+
     #                 BF(10)*SN(11))*PARM2+BF(8)*SN(9))*PARM2+
     #                 BF(6)*SN(7))*PARM2+BF(4)*SN(5))*PARM2+
     #                 BF(2)*SN(3))*(-PARM)+BF(1)*SN(2))*(-PARM)+
     #                 BF(0)*SN(1))*PARM
      CLI2R= SIGN1*(SIGN2*(SIGN3*RES(1)+ADD3(1))+ADD2(1))+ADD1(1)
      CLI2I= SIGN1*(SIGN2*(SIGN3*RES(2)+ADD3(2))+ADD2(2))+ADD1(2)
*
      ifirst = 1
      RETURN
      END
*
*-----CQLNX---------------------------------------------
*     COMPUTES  LN(Z)
*
      SUBROUTINE TCQLNX(ARG,RES)
! from TOPAZ0
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      DIMENSION ARG(2),AARG(2),RES(2)
*
      QPI= 3.141592653589793238462643D0
      DO I= 1,2
          AARG(I)= ABS(ARG(I))
      ENDDO
      ZM2= (ARG(1))**2+(ARG(2))**2
      ZM= SQRT(ZM2)
      RES(1)= LOG(ZM)
      IF(ARG(1).EQ.0.D0) THEN
          IF(ARG(2).GT.0.D0) THEN
              TETA= QPI/2.D0
          ELSE
              TETA= -QPI/2.D0
          ENDIF
          RES(2)= TETA
          RETURN
      ELSE IF(ARG(2).EQ.0.D0) THEN
               IF(ARG(1).GT.0.D0) THEN
                   TETA= 0.D0
               ELSE
                   TETA= QPI
               ENDIF
          RES(2)= TETA
          RETURN
      ELSE
          TNTETA= AARG(2)/AARG(1)
          TETA= ATAN(TNTETA)
          SR= ARG(1)/AARG(1)
          SI= ARG(2)/AARG(2)
          IF(SR.GT.0.D0) THEN
              RES(2)= SI*TETA
          ELSE
              RES(2)= SI*(QPI-TETA)
          ENDIF
          RETURN
      ENDIF
      END
*
*-----CQLNOMX---------------------------------------
*     COMPUTES LN(1-X)
*     USUALLY |X| << 1
      SUBROUTINE TCQLNOMX(ARG,OMARG,RES)
! from TOPAZ0
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ARG(2),OMARG(2),RES(2),ARES(2),CT(10),SN(10)
*
      ZR= ARG(1)
      ZI= ARG(2)
      ZM2= ZR*ZR+ZI*ZI
      ZM= SQRT(ZM2)
      IF(ZM.LT.1.D-7) THEN
          CT(1)= ZR/ZM
          SN(1)= ZI/ZM
          DO N=2,10
              CT(N)= CT(1)*CT(N-1)-SN(1)*SN(N-1)
              SN(N)= SN(1)*CT(N-1)+CT(1)*SN(N-1)
          ENDDO
          ARES(1)= CT(10)/10.D0
          ARES(2)= SN(10)/10.D0
          DO K=9,1,-1
              ARES(1)= ARES(1)*ZM+CT(K)/K
              ARES(2)= ARES(2)*ZM+SN(K)/K
          ENDDO
          ARES(1)= -ARES(1)*ZM
          ARES(2)= -ARES(2)*ZM
      ELSE
          CALL TCQLNX(OMARG,ARES)
      ENDIF
      RES(1)= ARES(1)
      RES(2)= ARES(2)
      RETURN
      END


c==================================================================
C this sort subroutine works very well...
c
C From Leonard J. Moss of SLAC:
C Here's a hybrid QuickSort I wrote a number of years ago.  It's
C based on suggestions in Knuth, Volume 3, and performs much better
C than a pure QuickSort on short or partially ordered input arrays.
      SUBROUTINE SORTRX(N,DATA,INDEX)
C===================================================================
C     SORTRX -- SORT, Real*8 input, indeX output
C     Input:  N     INTEGER
C             DATA  REAL
C     Output: INDEX INTEGER (DIMENSION N)
C This routine performs an in-memory sort of the first N elements of
C array DATA, returning into array INDEX the indices of elements of
C DATA arranged in ascending order.  Thus,
C    DATA(INDEX(1)) will be the smallest number in array DATA;
C    DATA(INDEX(N)) will be the largest number in DATA.
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C===================================================================
C SORTRX uses a hybrid QuickSort algorithm, based on several
C suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
C "pivot key" [my term] for dividing each subsequence is chosen to be
C the median of the first, last, and middle values of the subsequence;
C and the QuickSort is cut off when a subsequence has 9 or fewer
C elements, and a straight insertion sort of the entire array is done
C at the end.  The result is comparable to a pure insertion sort for
C very short arrays, and very fast for very large arrays (of order 12
C micro-sec/element on the 3081K for arrays of 10K elements).  It is
C also not subject to the poor performance of the pure QuickSort on
C partially ordered data.
C Created:  15 Jul 1986  Len Moss
C===================================================================
      INTEGER   N,INDEX(N)
      DOUBLE PRECISION      DATA(N)
      INTEGER   LSTK(31),RSTK(31),ISTK
      INTEGER   L,R,I,J,P,INDEXP,INDEXT
      DOUBLE PRECISION      DATAP
C     QuickSort Cutoff
C     Quit QuickSort-ing when a subsequence contains M or fewer
C     elements and finish off at end with straight insertion sort.
C     According to Knuth, V.3, the optimum value of M is around 9.
      INTEGER   M
CARLO
cc      PARAMETER (M=30)
CARLO
      PARAMETER (M=9)

C===================================================================
C
C     Make initial guess for INDEX

      DO 50 I=1,N
         INDEX(I)=I
   50    CONTINUE

C     If array is short, skip QuickSort and go directly to
C     the straight insertion sort.

      IF (N.LE.M) GOTO 900

C===================================================================
C
C     QuickSort
C
C     The "Qn:"s correspond roughly to steps in Algorithm Q,
C     Knuth, V.3, PP.116-117, modified to select the median
C     of the first, last, and middle elements as the "pivot
C     key" (in Knuth's notation, "K").  Also modified to leave
C     data in place and produce an INDEX array.  To simplify
C     comments, let DATA[I]=DATA(INDEX(I)).

C Q1: Initialize
      ISTK=0
      L=1
      R=N

  200 CONTINUE

C Q2: Sort the subsequence DATA[L]..DATA[R].
C
C     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
C     r > R, and L <= m <= R.  (First time through, there is no
C     DATA for l < L or r > R.)

      I=L
      J=R

C Q2.5: Select pivot key
C
C     Let the pivot, P, be the midpoint of this subsequence,
C     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
C     so the corresponding DATA values are in increasing order.
C     The pivot key, DATAP, is then DATA[P].

      P=(L+R)/2
      INDEXP=INDEX(P)
      DATAP=DATA(INDEXP)

      IF (DATA(INDEX(L)) .GT. DATAP) THEN
         INDEX(P)=INDEX(L)
         INDEX(L)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF

      IF (DATAP .GT. DATA(INDEX(R))) THEN
         IF (DATA(INDEX(L)) .GT. DATA(INDEX(R))) THEN
            INDEX(P)=INDEX(L)
            INDEX(L)=INDEX(R)
         ELSE
            INDEX(P)=INDEX(R)
         ENDIF
         INDEX(R)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF

C     Now we swap values between the right and left sides and/or
C     move DATAP until all smaller values are on the left and all
C     larger values are on the right.  Neither the left or right
C     side will be internally ordered yet; however, DATAP will be
C     in its final position.

  300 CONTINUE

C Q3: Search for datum on left >= DATAP
C
C     At this point, DATA[L] <= DATAP.  We can therefore start scanning
C     up from L, looking for a value >= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).

         I=I+1
         IF (DATA(INDEX(I)).LT.DATAP) GOTO 300

  400 CONTINUE

C Q4: Search for datum on right <= DATAP
C
C     At this point, DATA[R] >= DATAP.  We can therefore start scanning
C     down from R, looking for a value <= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).

         J=J-1
         IF (DATA(INDEX(J)).GT.DATAP) GOTO 400

C Q5: Have the two scans collided?

      IF (I.LT.J) THEN

C Q6: No, interchange DATA[I] <--> DATA[J] and continue

         INDEXT=INDEX(I)
         INDEX(I)=INDEX(J)
         INDEX(J)=INDEXT
         GOTO 300
      ELSE

C Q7: Yes, select next subsequence to sort
C
C     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
C     for all L <= l < I and J < r <= R.  If both subsequences are
C     more than M elements long, push the longer one on the stack and
C     go back to QuickSort the shorter; if only one is more than M
C     elements long, go back and QuickSort it; otherwise, pop a
C     subsequence off the stack and QuickSort it.

         IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=J+1
            RSTK(ISTK)=R
            R=I-1
         ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=L
            RSTK(ISTK)=I-1
            L=J+1
         ELSE IF (R-J .GT. M) THEN
            L=J+1
         ELSE IF (I-L .GT. M) THEN
            R=I-1
         ELSE
C Q8: Pop the stack, or terminate QuickSort if empty
            IF (ISTK.LT.1) GOTO 900
            L=LSTK(ISTK)
            R=RSTK(ISTK)
            ISTK=ISTK-1
         ENDIF
         GOTO 200
      ENDIF
  900 CONTINUE
C===================================================================
C Q9: Straight Insertion sort
      DO 950 I=2,N
         IF (DATA(INDEX(I-1)) .GT. DATA(INDEX(I))) THEN
            INDEXP=INDEX(I)
            DATAP=DATA(INDEXP)
            P=I-1
  920       CONTINUE
               INDEX(P+1) = INDEX(P)
               P=P-1
               IF (P.GT.0) THEN
                  IF (DATA(INDEX(P)).GT.DATAP) GOTO 920
               ENDIF
            INDEX(P+1) = INDEXP
         ENDIF
  950    CONTINUE
C===================================================================
C     All done
      END

*********************************************************************
* Function to calculate the width of the U boson knowing its mass and
* its vectorial coupling
*********************************************************************
      double precision function gammafunc(mass,k)
      implicit none
      double precision mass,k

      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi

      double precision me, mmu, r, dr

      me=ame
      mmu=ammu

      gammafunc=0.

      if(mass.gt.2.*me) then
         gammafunc = gammafunc + 1./3.*alpha*k**2*mass*
     >               sqrt(1-4.*me**2/mass**2)*(1+2.*me**2/mass**2)
      endif
      if(mass.gt.2.*mmu) then
         gammafunc = gammafunc + 1./3.*alpha*k**2*mass*
     >               sqrt(1-4.*mmu**2/mass**2)*(1+2.*mmu**2/mass**2)
         call r_intpl(mass,r,dr)
         gammafunc = gammafunc + 1./3.*alpha*k**2*mass*
     >               sqrt(1-4.*mmu**2/mass**2)*(1+2.*mmu**2/mass**2)*r
      endif

      return

      end

**********************************************************************
*
* Sceglie il fotone da calcolare con l'elemento di matrice a caso fra
* quelli che rientrano nei tagli angolari
*********************************************************************
      subroutine sceglifot(ifscelto,ng)
      implicit none
      integer ifscelto,ng,nrilev,nphotscelto
      integer k,i
      real*8 csi(1)
      integer ifnontagliati(40)
      common/fnontagliati/ifnontagliati

      ifscelto = 1
      if (ng.eq.1) return

      nrilev = 0
      do k=1,ng
        nrilev = nrilev + ifnontagliati(k)
      enddo
      if (nrilev.eq.1) then
          nphotscelto = 1
      else
          call getrnd(csi,1)
          nphotscelto = csi(1)*nrilev + 1
      endif
      i = 0
      k = 1
      do while (k.le.ng)
        if (ifnontagliati(k).eq.1) i = i + 1
        if (i.eq.nphotscelto) then
            ifscelto = k
            return
        endif
        k=k+1
      enddo

      return
      end
*
      double precision function trinorma(p)
      implicit none
      double precision p(0:3), tridot

      trinorma=sqrt(tridot(p,p))

      return
      end


*************************************************************
*     FROM nunugpv
*
*.....FUNCTION DOP: D(X)/P(X)/ANORMD(XMIN,XMAX). D(X) IS THE
*                   ELECTRON STRUCTURE FUNCTION
*
*************************************************************

      double precision function dop(x,am,p1,p2)
      implicit none
      double precision x, omx, opx, am
      double precision s
      double precision alx, alomx, alox
      double precision beta, betam, sdelta
      double precision dfunh, pm1
      double precision dot
      integer kk

      double precision ame, ammu, convfac, alpha, pi
      common/parameters/ame,ammu,convfac,alpha,pi

      double precision p1(0:3),p2(0:3),ptmp(0:3)

      save beta, betam, sdelta


      do kk = 0,3
        ptmp(kk)   = p1(kk)+p2(kk)
      enddo
      s = dot(ptmp,ptmp)
*
      beta = 2.d0*alpha/pi*(dlog(s/am**2) - 1.d0)
      betam = beta/2.d0
*
      sdelta = 1+3.d0/8.d0*beta+(9.d0/128.d0-pi**2/48.d0)*beta**2
*
      omx = 1.D0 - x
      opx = 1.D0 + x
*
      alx = dlog(x)
*
      if (omx.LT.1D-28) then
         alomx = 0.D0
         alox = -1.D0
      else
         alomx = dlog(omx)
         alox  = alx/omx
      endif
*
      pm1 = 1.d0/betam*(1.d0 - x)**(1.d0 - betam)

      dfunh = - beta/4.D0*opx
     #        + beta**2/32.D0*(
     #        - 4.d0*opx*alomx + 3.d0*opx*alx - 4.d0*alox
     #        - 5.d0 - x)
      dop = sdelta + dfunh*pm1
*
      return
      end

****************************************************************************
      subroutine emission (p1, p2, w, particle)
      implicit none
      double precision p1(0:3), p2(0:3), ptmp(0:3)
      double precision dir1(1:3), dir2(1:3)
      real*8 csi(2)
      double precision x1,x2
      double precision s
      double precision w
      double precision x1max, x2max, x1min, x2min
      character*1 particle
      integer kk
      double precision dop
      double precision tridot, dot, am
      double precision and1, and2
      double precision ame, ammu, convfac, alpha, pi
      common/parameters/ame,ammu,convfac,alpha,pi
      double precision beta,betam
      double precision amassainv
      character*6   ord
      common/qedORDER/ord

      if (ord.eq.'born') then
        w = 1.d0
        return
      endif

      if (ord.eq.'struct') then

        do kk=0,3
        ptmp(kk)=p1(kk)+p2(kk)
        enddo

        s = dot(ptmp,ptmp)

        if (particle.eq.'e') then
          am = ame
        elseif (particle.eq.'m') then
          am = ammu
        endif

        beta  = 2.d0*alpha/pi * (dlog(s/(am**2))-1.d0)
        betam = beta/2.d0

        x1min = am/p1(0)
        x2min = am/p2(0)
        x1max = 1.d0
        x2max = 1.d0

        call getrnd(csi,2)

        x1=1.d0-((1.D0 - x1min)**(betam)*(1.D0-csi(1)*1.d0)
     >         +(1.D0 - x1max)**(betam)*csi(1))**(1.d0/betam)
*
        x2=1.d0-((1.D0 - x2min)**(betam)*(1.D0-csi(2)*1.d0)
     >         +(1.D0 - x2max)**(betam)*csi(2))**(1.d0/betam)

        w = dop(x1,am,p1,p2) * dop(x2,am,p1,p2)

        p1(0) = x1*p1(0)
        p2(0) = x2*p2(0)

        do kk=1,3
          dir1(kk)=p1(kk)/dsqrt(tridot(p1,p1))
          dir2(kk)=p2(kk)/dsqrt(tridot(p2,p2))
        enddo

        do kk=1,3
          p1(kk) = dir1(kk)*dsqrt(p1(0)**2-am**2)
          p2(kk) = dir2(kk)*dsqrt(p2(0)**2-am**2)
        enddo

        and1 = (1.D0-x1min)**(betam)-(1.D0-x1max)**(betam)
        and2 = (1.D0-x2min)**(betam)-(1.D0-x2max)**(betam)

        w = w * and1 * and2

        return
      endif
      end
