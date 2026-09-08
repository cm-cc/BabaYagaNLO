      subroutine glueALPHA(ng,p1,p2,p3,p4,qph,emtx2)
      implicit double precision (a-h,o-z)
      dimension qph(ng,0:3)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      include 'prm.f'
      logical optspn
      double precision p(4,4),pph(4,nphmx),mss(4)
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      double precision xi(1)
      integer ifirst
      common/gluealphacommon/ifirst
      data ifirst/0/
      character*2 fs
      common/finalstate/fs

      integer jspn_f, jspn_ph
      real*8 r_ph
      parameter (NFERMION_SPIN = 16) ! 2^4 = 16
      integer jspn_ph_max, jspn_ph_loop, nfer_red, jspn_f_red
      integer jspn_f_full, jspn_total, bit_b0, bit_b1
      integer jspn_ph_fixed, jspn_f_loop

      if (ifirst.eq.0) then
         ifirst = 1
         call random_init(.true.,.true.)
      endif

      mss(1) = ame
      mss(3) = ame
      if (fs(1:1).eq.'m') then
         mss(2) = ammu
         mss(4) = ammu
      elseif (fs(1:1).eq.'e') then
         mss(2) = ame
         mss(4) = ame
      elseif (fs(1:1).eq.'p') then
         mss(1) = ame
         mss(2) = ame
         mss(3) = ampi
         mss(4) = ampi   
      endif

      p(1:4,1) = -p1(0:3)
      p(1:4,2) =  p4(0:3)
      p(1:4,3) = -p2(0:3)
      p(1:4,4) =  p3(0:3)

      if (fs(1:1).eq.'p') then
         p(1:4,1) = -p1(0:3)
         p(1:4,2) = -p2(0:3)
         p(1:4,3) =  p4(0:3)
         p(1:4,4) =  p3(0:3)
      endif

      do k = 1,ng
         pph(1:4,k) = qph(k,0:3)
      enddo

      extraw = 1.d0
      optspn = .false.
      emtx2 = 0.d0
      jspn  = 0

!     sum only over fermions, + random photons
c      jspn_ph_max = 2**ng
c      call rndmgnrtr(r_ph)
c      jspn_ph = jspn_ph_max * r_ph
c      do jspn_f = 0, NFERMION_SPIN - 1
c         jspn = jspn_f + jspn_ph * NFERMION_SPIN
c         call eemumu(p, jspn, optspn, mss, ng, pph, emtxalpha)
c         emtx2 = emtx2 + emtxalpha
c      enddo
c      extraw = 1.d0 * jspn_ph_max

! sum over half fermion spin configurations!!
c$$$      jspn_ph_max = 2**ng
c$$$      aa = 0.d0
c$$$      bb = 0.d0
c$$$      cc = 0.d0
c$$$      dd = 0.d0
c$$$      do jspn_f = 0, NFERMION_SPIN- 1
c$$$         do jspn_ph = 0,jspn_ph_max-1
c$$$            jspn = jspn_f + jspn_ph * NFERMION_SPIN
c$$$
c$$$c            print*,mod(jspn,2),mod(jspn/2,2),mod(jspn/4,2),mod(jspn/8,2)
c$$$
c$$$            call eemumu(p, jspn, optspn, mss, ng, pph, emtxalpha)
c$$$            emtx2 = emtx2 + emtxalpha
c$$$c            if (jspn_f.le.3) aa = aa + emtxalpha
c$$$c            if (jspn_f.gt.3.and.jspn_f.le.7) bb = bb + emtxalpha
c$$$c            if (jspn_f.gt.7.and.jspn_f.le.11) cc = cc + emtxalpha
c$$$c            if (jspn_f.gt.11) dd = dd + emtxalpha
c$$$c            print*,jspn_f,jspn_ph,emtxalpha
c$$$         enddo
c$$$      enddo
c$$$
c$$$c      print*,aa*(4.d0*pi*alpha)**(ng+2)
c$$$c      print*,dd*(4.d0*pi*alpha)**(ng+2)
c$$$c      print*,bb*(4.d0*pi*alpha)**(ng+2)
c$$$c      print*,cc*(4.d0*pi*alpha)**(ng+2)
c$$$c      print*,(aa-bb)*(4.d0*pi*alpha)**(ng+2)
c$$$c      print*,NFERMION_SPIN- 1
c$$$c      print*,NFERMION_SPIN/2- 1
c$$$c      extraw = 2.d0
c$$$c      extraw = 1.d0
c$$$

c      p(2:4,:)   = -p(2:4,:)
c      pph(2:4,:) = -pph(2:4,:)

      if (fs(1:1).eq.'m') then
! sum over spin configurations
         do jspn = 0,2**(4+ng)-1
c      do jspn = 2**(4-1+ng),2**(4+ng)-1
c      do jspn = 0,2**(4-1+ng)-1
            call eemumu(p,jspn,optspn,mss,ng,pph,emtxalpha)
            emtx2 = emtx2 + emtxalpha
         enddo
!     sum over spin configurations, but reduced

! random spin configuration
c$$$      extraw = 2**(4+ng)
c$$$      call rndmgnrtr(r)
c$$$      jspn = extraw*r
c$$$      optspn = .false.
c$$$      call eemumu(p,jspn,optspn,mss,ng,pph,emtxalpha)
c$$$      emtx2 = emtx2 + emtxalpha

      elseif(fs(1:1).eq.'e') then
         do jspn = 0,2**(4+ng)-1
c      do jspn = 2**(4-1+ng),2**(4+ng)-1
c      do jspn = 0,2**(4-1+ng)-1
            call eeee(p,jspn,optspn,mss,ng,pph,emtxalpha)
            emtx2 = emtx2 + emtxalpha
         enddo
      elseif(fs(1:1).eq.'p') then
         do jspn = 0,2**(2+ng)-1
            call eepipi(p,jspn,optspn,mss,ng,pph,emtxalpha)
!            print*,emtxalpha
            emtx2 = emtx2 + emtxalpha
         enddo

      endif

      emtx2 = (4.d0*pi*alpha)**(ng+2)*emtx2*extraw

      return
      end
****************************
      subroutine m4v(a, x, y)
      implicit none
      complex*16, intent(in)  :: a(4,4), x(4)
      complex*16, intent(out) :: y(4)
      integer :: i

      do i = 1, 4
         y(i) = a(i,1)*x(1) + a(i,2)*x(2) + a(i,3)*x(3) + a(i,4)*x(4)
      enddo
      end subroutine m4v
      subroutine v4m(x, a, y)
      implicit none
      complex*16, intent(in)  :: x(4), a(4,4)
      complex*16, intent(out) :: y(4)

      y(1) = x(1)*a(1,1) + x(2)*a(2,1) + x(3)*a(3,1) + x(4)*a(4,1)
      y(2) = x(1)*a(1,2) + x(2)*a(2,2) + x(3)*a(3,2) + x(4)*a(4,2)
      y(3) = x(1)*a(1,3) + x(2)*a(2,3) + x(3)*a(3,3) + x(4)*a(4,3)
      y(4) = x(1)*a(1,4) + x(2)*a(2,4) + x(3)*a(3,4) + x(4)*a(4,4)
      end subroutine v4m

! FROM MAURO MORETTI
!     The matrix element is returned by the subroutine EEMUMU, see below
!     The maximum number of photons can be changed in the file PRM.f
!     in the OPTSPN=.TRUE. mode the helicity configuration is somewhat
!     optimized to suppress helicity flipping configuration and uses
!     a random number to perform an hit or miss. The random number is
!     generated in the subroutine RNDMGNRTR which uses fortran built in.
!     you can change random generator there
C***********************************************************************
      subroutine rndmgnrtr(r)
C***********************************************************************
!     It return the random number R, 0<=R<=1, for spin optimization.
!     It uses fortran built in, if you wish change here with your favourite
!      routine
!
      implicit none
      real*8 r
      real r4
      double precision xi(1)
      call random_number(r4)
      r = r4
c      call getrnd(xi,1)
c      r = xi(1)
      return
      end
C***********************************************************************
      function photonpropagatorpipi(q)
C***********************************************************************
!
!     return the photon propagator, modify here to include vacuum polarization
!
      implicit none
!
      complex*16 photonpropagatorpipi
      real*8 q(4)
      real*8 rrlrnzprdct
            
      double complex vpolc
      external vpolc
!
      integer iffpi,iffpibck,vals,u,error

      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi
      complex*16 ffpic
      real*8 ffpi2
!
      real*8 qq
!
      qq=rrlrnzprdct(q,q)
      call ffpi(qq,ffpic,ffpi2)
c      if(qq.lt.3.8d0) then
      photonpropagatorpipi = (1.d0,0.d0)/qq*ffpic
c      else
c         photonpropagator = (1.d0,0.d0)/qq
c      endif
!      write(*,*)'qq',qq,'   prpg',photonpropagator
!
      return
      end

      function photonpropagator(q)
C***********************************************************************
!
!     return the photon propagator, modify here to include vacuum polarization
!
      implicit none
!
      complex*16 photonpropagator
      real*8 q(4)
      double complex vpolc
      external vpolc
!
      real*8 rrlrnzprdct
!
      real*8 qq
!
      qq=rrlrnzprdct(q,q)
c      if(qq.lt.3.8d0) then
      photonpropagator = (1.d0,0.d0)/qq * vpolc(qq)
c      else
c         photonpropagator = (1.d0,0.d0)/qq
c      endif
!      write(*,*)'qq',qq,'   prpg',photonpropagator
!
      return
      end

C***********************************************************************
      subroutine eemumu(p,jspn,optspn,mss,nph,pph,mesq)
C***********************************************************************
!     Input: P(i,j), fermion momentum; i=E,px,py,pz; j fermion label
!            j=1,_4 = e-,mu+,e+,mu-; incoming momenta (e+-) have REVERSED
!            sign (E<0, pj->-pj)
!            JSPN=0,2**(NFER+NPH)-1 label to select helicity configuration
!            OPTSPN, for optimization over helicity configuration. If .TRUE.
!              then a hit or miss is performed and if fails MESQ=0 is returned
!              if .FALSE. all helicities come with weight=1
!            MSS(i) fermion masses
!            NPH, photons number, NPH <=NPHMX which is a parameter in PRM.F
!            PPH(i,j); photon momentum of the j-th photon
!     Outpun: MESQ, squared matrix element. NB alphaem = 1/(4 pi)
!
      implicit none
!
      include 'prm.f'
!
      integer nfer
      parameter (nfer=4)
      integer nph
      real*8 p(4,nfer),pph(4,nphmx),mss(nfer)
      integer jspn
!
      logical OPTSPN
      real*8 mesq
!
      complex*16 gam(4,4,5),iden(4,4) !N.B (gamma_\mu)_{ij} ==
      common/gamma/gam,iden                       ! gam(j,i,mu)
!

      data gam/(1.,0.),(0.,0),(0.,0.),(0.,0.),(0.,0.)
     >     ,(1.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     >     (-1.,0.),(0.,0.),(0.,0.),(0.,0.), (0.,0.),(-1.,0.), !gamma0
     >     (0.,0.),(0.,0.),(0.,0.),(1.,0.),(0.,0.),
     >     (0.,0.),(1.,0.),(0.,0.),(0.,0.),(-1.,0.),
     >     (0.,0.),(0.,0.),(-1.,0.),(0.,0.),(0.,0.),(0.,0.), !gamma1
     >     (0.,0.),(0.,0.),(0.,0.),(0.,-1.),(0.,0.),
     >     (0.,0.),(0.,1.),(0.,0.),(0.,0.),(0.,1.),
     >     (0.,0.),(0.,0.),(0.,-1.),(0.,0.),(0.,0.), (0.,0.), !gamma2
     >     (0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.),
     >     (0.,0.),(0.,0.),(-1.,0.),(-1.,0.),(0.,0.),
     >     (0.,0.),(0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.), !gamma3
     >     (0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.),
     >     (0.,0.),(0.,0.),(1.,0.),(1.,0.),(0.,0.),
     >     (0.,0.),(0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.)/ !gamma5
      data iden/(1.,0.), (0.,0.), (0.,0.), (0.,0.),
     *          (0.,0.), (1.,0.),(0.,0.), (0.,0.), (0.,0.), (0.,0.),
     *          (1.,0.), (0.,0.),(0.,0.), (0.,0.), (0.,0.), (1.,0.)/

      real*8 spn(nfer),epsph(nphmx)
!            SPN(i); +-0.49 fermions; +-0.51 antifermions
!            EPSPH(i); = +-1.d0 polarization of the jth photon
      complex*16 psi(4,0:2**nphmx-1,nfer)
      integer tblmrg(nphmx,0:2**nphmx-1)
      complex*16 curree(4,0:2**nphmx-1),currmumu(4,0:2**nphmx-1),amp
      complex*16 cclrnzprdct,photonpropagator
      real*8 rrlrnzprdct
      logical frst/.true./,mrg
      integer j1,j2,j3,jj
      save frst,tblmrg
      real*8 q(4),s1,s2,pspne,pspnm,prb,r
      complex*16 :: tmpvec(4)
      double precision time,time1
!
!      call cpu_time(time1)
!      print*,'ciao'
      jj = jspn
      spn(1)=0.49d0
      if (mod(jj,2).eq.0) spn(1)=-0.49d0

      jj=jj/2
      spn(2)=0.51d0
      if (mod(jj,2).eq.0) spn(2)=-0.51d0

      jj=jj/2
      spn(3)=0.51d0
      if (mod(jj,2).eq.0) spn(3)=-0.51d0

      jj=jj/2
      spn(4)=0.49d0
      if (mod(jj,2).eq.0) spn(4)=-0.49d0

      s1=sign(1.d0,spn(1))
      s2=sign(1.d0,spn(2))
      s1=s1*sign(1.d0,spn(3))
      s2=s2*sign(1.d0,spn(4))

      jj=jj/2
      do j1=1,nph
         epsph(j1)=1.d0
         if (mod(jj,2).eq.0) epsph(j1)=-1.d0
         jj=jj/2
      enddo
      if(s1.gt.0.d0) then
         pspne=1.d0
      else
         pspne=(mss(1)/min(p(1,1),p(1,3)))**2*1.d1
         pspne=max(0.1d0,pspne)
         pspne=min(1.d0,pspne)
      endif
      if(s2.gt.0.d0) then
         pspnm=1.d0
      else
         pspnm=(mss(2)/min(p(1,2),p(1,4)))**2*1.d1
         pspnm=max(0.1d0,pspnm)
         pspnm=min(1.d0,pspnm)
      endif
      prb=pspne*pspnm
      if (optspn) then
         call rndmgnrtr(r)
         if(r.lt.prb) then
            prb=1.d0/prb
         else
            mesq=0.d0
         endif
      else
         prb=1.d0
      endif
!
!      write(*,*)'spn  ',spn,'  eps',epsph(1:nph)
      call absphotonfld(p,nfer,spn,mss,nph,pph,epsph,psi)
!
*** CHATGPT
c$$$      if(frst) then
c$$$         frst=.false.
c$$$         call mrgng(tblmrg)
c$$$      endif
c$$$      curree=(0.d0,0.d0)
c$$$      currmumu=(0.d0,0.d0)
c$$$      do j1=0,2**nph-1
c$$$         do j2=j1,2**nph-1
c$$$            if(j2.eq.j1.and.j1.ne.0) cycle
c$$$!         do j2=0,2**nph-1
c$$$            mrg=.true.
c$$$            do j3=1,nph
c$$$               mrg=mrg.and.(tblmrg(j3,j1)+tblmrg(j3,j2).le.1)
c$$$            enddo
c$$$            if(mrg) then
c$$$               jj=j1+j2
c$$$               do j3=1,4
c$$$c                 temporary 4-vector for replacing matmul results
c$$$c                  complex*16 :: tmpvec(4)
c$$$
c$$$c                 curree: conjg(psi(1:4,j2,3)) . matmul( psi(1:4,j1,1), gam(:,:,j3) )
c$$$                  call v4m(psi(1:4,j1,1), gam(1:4,1:4,j3), tmpvec)
c$$$                  curree(j3,jj) = dot_product(dconjg(psi(1:4,j2,3)),
c$$$     .                 tmpvec) + curree(j3,jj)
c$$$
c$$$                  if(j1.ne.j2) then
c$$$                     call v4m(psi(1:4,j2,1), gam(1:4,1:4,j3), tmpvec)
c$$$                     curree(j3,jj) = dot_product(dconjg(psi(1:4,j1,3)),
c$$$     .                    tmpvec) + curree(j3,jj)
c$$$                  endif
c$$$
c$$$c                 currmumu: conjg(psi(1:4,j1,4)) . matmul( psi(1:4,j2,2), gam(:,:,j3) )
c$$$                  call v4m(psi(1:4,j2,2), gam(1:4,1:4,j3), tmpvec)
c$$$                  currmumu(j3,jj) = dot_product(dconjg(psi(1:4,j1,4)),
c$$$     .                 tmpvec) + currmumu(j3,jj)
c$$$
c$$$                  if(j1.ne.j2) then
c$$$                     call v4m(psi(1:4,j1,2), gam(1:4,1:4,j3), tmpvec)
c$$$                     currmumu(j3,jj) = dot_product(dconjg(psi(1:4,j2,4)),
c$$$     .                    tmpvec) + currmumu(j3,jj)
c$$$                  endif
c$$$
c$$$               enddo
c$$$            endif
c$$$         enddo
c$$$      enddo
*** END CHATGPT
            
      if(frst) then
         frst=.false.
         call mrgng(tblmrg)
      endif
      
      curree=(0.d0,0.d0)
      currmumu=(0.d0,0.d0)
      do j1=0,2**nph-1
         do j2=j1,2**nph-1
            if(j2.eq.j1.and.j1.ne.0) cycle
!         do j2=0,2**nph-1
            mrg=.true.
            do j3=1,nph
               mrg=mrg.and.(tblmrg(j3,j1)+tblmrg(j3,j2).le.1)
            enddo
            if(mrg) then
               jj=j1+j2
               do j3=1,4
                  curree(j3,jj)=dot_product(dconjg(psi(1:4,j2,3)),
     .                 matmul(psi(1:4,j1,1),gam(1:4,1:4,j3)) )
     .                 +curree(j3,jj)  !  gamma matrix are transposed
                  if(j1.ne.j2) curree(j3,jj)=
     .                 dot_product(dconjg(psi(1:4,j1,3)),
     .                 matmul(psi(1:4,j2,1),gam(1:4,1:4,j3)) )
     .                 +curree(j3,jj)  !  gamma matrix are transposed
                  currmumu(j3,jj)=dot_product(dconjg(psi(1:4,j1,4)),
     .                 matmul(psi(1:4,j2,2),gam(1:4,1:4,j3)) )
     .                 +currmumu(j3,jj)  !  gamma matrix are transposed
                  if(j1.ne.j2) currmumu(j3,jj)=
     .                 dot_product(dconjg(psi(1:4,j2,4)),
     .                 matmul(psi(1:4,j1,2),gam(1:4,1:4,j3)) )
     .                 +currmumu(j3,jj)  !  gamma matrix are transposed
!     per spegnere radiazione dal mu
!                  if(jj.ne.0) currmumu(1:4,jj)=0.d0
               enddo
            endif
         enddo
      enddo
      amp=0.d0
      jj=2**nph-1
      do j1=0,2**nph-1
!     test
!!!         write(*,*)'j1, crre',j1, curree(1:4,j1)
!!!         write(*,*)'j1, crrm',j1, currmumu(1:4,j1)
!!!         write(*,*)'***'
!     test
         q=p(1:4,1)+p(1:4,3)
         do j2=1,nph
            if(tblmrg(j2,j1).eq.1) q=q+pph(1:4,j2)
         enddo
         amp=amp+cclrnzprdct(curree(1:4,j1),currmumu(1:4,jj-j1))
     .        *photonpropagator(q)
c         write(*,*)'xx', cclrnzprdct(curree(1:4,j1),currmumu(1:4,jj-j1))
c     .        **photonpropagator(q),'  amp',amp
      enddo
c      mesq=abs(amp)**2*prb
      mesq=amp*dconjg(amp)*prb
!      call cpu_time(time)
!      print*,time-time1
      
      return
      end
C***********************************************************************
      subroutine mrgng(tblmrg)
C***********************************************************************
      implicit none
c
c     Output: TBMRG(I,J), J represent a set of photons from 1 to NPHMX
c             TBMRG(I,J)=.TRUE. if I-th photon don't belong to the J set
c             NPHABS(J)=.FALSE. it there is only one photon in J set
c             TBMLTPL(I,J), I number of photons, TBMLTPL(0,J) number of
c             photon sets with I photons, TBMLTPL(I,1-n) labels of sets
c             with I photons
c
      include 'prm.f'
c
      integer tblmrg(nphmx,0:2**nphmx-1)
c
      integer j1,jj0,j2
c
      do j1=0,2**nphmx-1
         jj0=j1
         do j2=1,nphmx
            tblmrg(j2,j1)=mod(jj0,2)
            jj0=jj0/2
         enddo
      enddo
c
      return
      end
C***********************************************************************
      subroutine absphotonfld(pfld,nfer,spnfld,mss,nphot,pphot,epsphot
     .     ,psi)
C***********************************************************************
      implicit none
c
c     Input: pfld(i,j), field momentum, i == (E,px,py,pz), initial state fermions E<0,
c                       final state antifermions E>0 , j= j-th fermion enter fields (psi)
c                       before bar fields (psibar) namely initial state particle and final
c                       state antiparticle first
c            nfer, number of fermion fields
c            spnfld, field spn  +- 0.49 particles +- 0.51 antiparticles
c            mss, field mass
c            nphot, photon number
c            pphot(j,k), photon momentum, j=1,_,4 (E,px,py,pz), k= kth photon
c            epsphot(j), helicity of j-th photon
c     Output: psi(i1,i2,i3), field after the absorption of photons in any
c             possible combination. i1 spinor index; i2 label absorbed photons
c             i3 labels absorbing fermion
c
      include 'prm.f'
c
      integer nphot,nfer
      real*8 pfld(4,nfer),spnfld(nfer),pphot(4,nphmx),epsphot(nphmx)
      real*8 mss(nfer)
c
      complex*16 psi(4,0:2**nphmx-1,nfer)
c
      integer tblmrg(nphmx,0:2**nphmx-1)
      logical frst/.true./,tabnw(0:2**nphmx)
      save tblmrg,frst
c
      integer j1,j2,jlow,jhigh,jj,j3,jj0,idabs(2**nphmx),j4,j5
      integer idabsnw(2**nphmx)
      real*8 ptmp(4),kk(4),den,rrlrnzprdct,kksum(4,2**nphmx)
      complex*16 epssl(4,4,nphmx),kksl(4,4,2**nphmx)
      complex*16 psl(4,4),psitmp1(4)
      complex*16 psitmp(4),epstmp(4),frmprp(4,4,2**nphmx,nfer)
c
      if (frst) then
         frst=.false.
         call mrgng(tblmrg)
      endif
c
      do j1=1,nphot
         call photpol(pphot(1,j1),epsphot(j1),epstmp)
!     test
!!!         if(j1.eq.1) epstmp=pphot(1:4,1)
!     end test
         call slashcv(epstmp,epssl(1,1,j1))
!         call slashrv(pphot(1,j1),kksl(1,1,j1))
      enddo
c
      do j4=1,nfer
         ptmp=pfld(1:4,j4)
         call slashrv(ptmp,psl)
         do j1=1,2**nphot-1
            if(j4.eq.1) then
               kk=0.d0
               do j2=1,nphot
                  if(tblmrg(j2,j1).eq.1) then
                     kk=kk+pphot(1:4,j2)
                  endif
               enddo
               kksum(1:4,j1)=kk
               call slashrv(kk,kksl(1,1,j1))
            endif
            den=2.d0*rrlrnzprdct(pfld(1,j4),kksum(1,j1))
            den=den+rrlrnzprdct(kksum(1,j1),kksum(1,j1))
            den=1.d0/den
            frmprp(1:4,1:4,j1,j4)=kksl(1:4,1:4,j1)+psl
            if(j4.le.nfer/2) frmprp(1:4,1:4,j1,j4)=
     .           -frmprp(1:4,1:4,j1,j4)
            do j2=1,4
               frmprp(j2,j2,j1,j4)=frmprp(j2,j2,j1,j4)+mss(j4)
            enddo
            frmprp(1:4,1:4,j1,j4)=frmprp(1:4,1:4,j1,j4)*den
         enddo
      enddo
c
      do j5=1,nfer
         if(j5.le.nfer/2) then
            call fld(spnfld(j5),pfld(1,j5),mss(j5),psi(1,0,j5))
         else
            call fldbr(spnfld(j5),pfld(1,j5),mss(j5),psi(1,0,j5))
         endif
!     test
!!!$$$         ptmp=pfld(1:4,j5)
!!!$$$         call slashrv(ptmp,psl)
!!!$$$         if(j5.le.nfer/2) then
!!!$$$            write(*,*)psi(1:4,0,j5),matmul(psl,psi(1:4,0,j5))
!!!$$$         else
!!!$$$            write(*,*)psi(1:4,0,j5),matmul(psi(1:4,0,j5),psl)
!!!$$$         endif
!     end test
** CHATGPT
c$$$         jhigh=1
c$$$         idabs(1)=0
c$$$         tabnw=.true.
c$$$         do j4=1,nphot
c$$$            jj=0
c$$$            do j1=1,jhigh
c$$$               do j2=1,nphot
c$$$                  if(tblmrg(j2,idabs(j1)).eq.0) then
c$$$                     jj0=2**(j2-1)+idabs(j1)
c$$$                     if(tabnw(jj0)) then
c$$$                        psi(1:4,jj0,j5)=(0.d0,0.d0)
c$$$                        tabnw(jj0)=.false.
c$$$                        jj=jj+1
c$$$                        idabsnw(jj)=jj0
c$$$                     endif
c$$$                     if(j5.gt.nfer/2) then
c$$$c                   psitmp = matmul(epssl(:,:,j2), psi(1:4,idabs(j1),j5))
c$$$c                        complex*16 :: psitmp(4), psitmp1(4)
c$$$                        call m4v(epssl(1:4,1:4,j2),
c$$$     .                       psi(1:4,idabs(j1),j5), psitmp)
c$$$                        call m4v(frmprp(1:4,1:4,jj0,j5), psitmp,
c$$$     .                       psitmp1)
c$$$                     else
c$$$c                        complex*16 :: psitmp(4), psitmp1(4)
c$$$c     psitmp = matmul( psi(1:4,idabs(j1),j5), epssl(:,:,j2) )
c$$$                        call v4m(psi(1:4,idabs(j1),j5), epssl(1:4,1:4,j2), psitmp)
c$$$c     psitmp1 = matmul(psitmp, frmprp(:,:,jj0,j5))
c$$$                        call v4m(psitmp, frmprp(1:4,1:4,jj0,j5),
c$$$     .                       psitmp1)
c$$$                     endif
c$$$                     psi(1:4,jj0,j5)=psi(1:4,jj0,j5)+psitmp1
c$$$                  endif
c$$$               enddo
c$$$            enddo
c$$$            jhigh=jj
c$$$            idabs=idabsnw
c$$$         enddo
c$$$      enddo
** END CHATGPT
         jhigh=1
         idabs(1)=0
         tabnw=.true.
         do j4=1,nphot
            jj=0
            do j1=1,jhigh
               do j2=1,nphot
                  if(tblmrg(j2,idabs(j1)).eq.0) then
                     jj0=2**(j2-1)+idabs(j1)
                     if(tabnw(jj0)) then
                        psi(1:4,jj0,j5)=(0.d0,0.d0)
                        tabnw(jj0)=.false.
                        jj=jj+1
                        idabsnw(jj)=jj0
                     endif
                     if(j5.gt.nfer/2) then
                        psitmp=matmul(epssl(1:4,1:4,j2),
     .                       psi(1:4,idabs(j1),j5))
                        psitmp1=matmul(frmprp(1:4,1:4,jj0,j5),psitmp)
                     else
                        psitmp=matmul(psi(1:4,idabs(j1),j5),
     .                       epssl(1:4,1:4,j2))
                        psitmp1=matmul(psitmp,frmprp(1:4,1:4,jj0,j5))
                     endif
                     psi(1:4,jj0,j5)=psi(1:4,jj0,j5)+psitmp1
                  endif
               enddo
            enddo
            jhigh=jj
            idabs=idabsnw
         enddo
      enddo
c
      return
      end
C***********************************************************************
      subroutine slashcv(cv,cvsl)
C***********************************************************************
c
      implicit none
c
      complex*16 cvsl(4,4),cv(4),test(4,4)
c
      real*8 gmn(4)/1.d0,-1.d0,-1.d0,-1.d0/
      complex*16 gam(4,4,5),iden(4,4)             !N.B (gamma_\mu)_{ij} ==
      common/gamma/gam,iden                       ! gam(j,i,mu)
!                                 keep gamma matrixes here in common elsewhere
      data gam/(1.,0.),(0.,0.),(0.,0.),(0.,0.),
     >         (0.,0.),(1.,0.),(0.,0.),(0.,0.),
     >         (0.,0.),(0.,0.),(-1.,0.),(0.,0.),
     >         (0.,0.),(0.,0.),(0.,0.),(-1.,0.), !gamma0
     >         (0.,0.),(0.,0.),(0.,0.),(1.,0.),
     >         (0.,0.),(0.,0.),(1.,0.),(0.,0.),
     >         (0.,0.),(-1.,0.),(0.,0.),(0.,0.),
     >         (-1.,0.),(0.,0.),(0.,0.),(0.,0.), !gamma1
     >         (0.,0.),(0.,0.),(0.,0.),(0.,-1.),
     >         (0.,0.),(0.,0.),(0.,1.),(0.,0.),
     >         (0.,0.),(0.,1.),(0.,0.),(0.,0.),
     >         (0.,-1.),(0.,0.),(0.,0.), (0.,0.), !gamma2
     >         (0.,0.),(0.,0.),(1.,0.),(0.,0.),
     >         (0.,0.),(0.,0.),(0.,0.),(-1.,0.),
     >         (-1.,0.),(0.,0.),(0.,0.),(0.,0.),
     >         (0.,0.),(1.,0.),(0.,0.),(0.,0.), !gamma3
     >         (0.,0.),(0.,0.),(1.,0.),(0.,0.),
     >         (0.,0.),(0.,0.),(0.,0.),(1.,0.),
     >         (1.,0.),(0.,0.),(0.,0.),(0.,0.),
     >         (0.,0.),(1.,0.),(0.,0.),(0.,0.)/ !gamma5
      data iden/(1.,0.), (0.,0.), (0.,0.), (0.,0.),
     *          (0.,0.), (1.,0.),(0.,0.), (0.,0.), (0.,0.), (0.,0.),
     *          (1.,0.), (0.,0.),(0.,0.), (0.,0.), (0.,0.), (1.,0.)/
c
      integer j1,j2,j3
      double complex v(1:4),v1,v2,v3,v4

c      v(1) =  cv(1)
c      v(2) = -cv(2)
c      v(3) = -cv(3)
c      v(4) = -cv(4)
c      do j1=1,4
c         do j2=1,4
c            cvsl(j1,j2) =
c     &           gam(j1,j2,1) * v(1) +
c     &           gam(j1,j2,2) * v(2) +
c     &           gam(j1,j2,3) * v(3) +
c     &           gam(j1,j2,4) * v(4)
c         enddo
c      enddo
! or
c      cvsl(:,:) = gam(:,:,1)*v(1) + gam(:,:,2)*v(2)
c     .     + gam(:,:,3)*v(3) + gam(:,:,4)*v(4)
c      return

      cvsl(1,1) = cv(1)
      cvsl(2,1) = 0.d0
      cvsl(3,1) = -cv(4)
      cvsl(4,1) = -cv(2) + (0.d0,1.d0)*cv(3)

      cvsl(1,2) = 0.d0
      cvsl(2,2) = cv(1)
      cvsl(3,2) = -cv(2) - (0.d0,1.d0)*cv(3)
      cvsl(4,2) = cv(4)

      cvsl(1,3) =  cv(4)
      cvsl(2,3) =  cv(2) - (0.d0,1.d0)*cv(3)
      cvsl(3,3) = -cv(1)
      cvsl(4,3) = 0.d0

      cvsl(1,4) = cv(2) + (0.d0,1.d0)*cv(3)
      cvsl(2,4) = -cv(4)
      cvsl(3,4) = 0.d0
      cvsl(4,4) = -cv(1)

      return
      end

C***********************************************************************
      subroutine slashrv(rv,cvsl)
C***********************************************************************
c
      implicit none
c
      complex*16 cvsl(4,4),test(4,4)
      real*8 rv(4)
      double precision v(1:4)

c
      real*8 gmn(4)/1.d0,-1.d0,-1.d0,-1.d0/
      complex*16 gam(4,4,5),iden(4,4)             !N.B (gamma_\mu)_{ij} ==
      common/gamma/gam,iden                       ! gam(j,i,mu)
c
      integer j1,j2,j3
c
c      v(1) =  rv(1)
c      v(2) = -rv(2)
c      v(3) = -rv(3)
c      v(4) = -rv(4)
c      do j1=1,4
c         do j2=1,4
c            cvsl(j1,j2) =
c     &           gam(j1,j2,1) * v(1) +
c     &           gam(j1,j2,2) * v(2) +
c     &           gam(j1,j2,3) * v(3) +
c     &           gam(j1,j2,4) * v(4)
c         enddo
c      enddo
! or
c      cvsl(:,:) = gam(:,:,1)*v(1) + gam(:,:,2)*v(2)
c     .     + gam(:,:,3)*v(3) + gam(:,:,4)*v(4)
c      return
      cvsl(1,1) = rv(1)
      cvsl(2,1) = 0.d0
      cvsl(3,1) = -rv(4)
      cvsl(4,1) = -rv(2) + (0.d0,1.d0)*rv(3)

      cvsl(1,2) = 0.d0
      cvsl(2,2) = rv(1)
      cvsl(3,2) = -rv(2) - (0.d0,1.d0)*rv(3)
      cvsl(4,2) = rv(4)

      cvsl(1,3) =  rv(4)
      cvsl(2,3) =  rv(2) - (0.d0,1.d0)*rv(3)
      cvsl(3,3) = -rv(1)
      cvsl(4,3) = 0.d0

      cvsl(1,4) = rv(2) + (0.d0,1.d0)*rv(3)
      cvsl(2,4) = -rv(4)
      cvsl(3,4) = 0.d0
      cvsl(4,4) = -rv(1)
c
      return
      end
C***********************************************************************
      subroutine trnspmtvct(cm44,cv4,cv4out)
c************************************************************************
      implicit none
c
      complex*16 cm44(4,4),cv4(4),cv4out(4)
c
      integer j1,j2
c
      do j1=1,4
         cv4out(j1)=(0.d0,0.d0)
         do j2=1,4
            cv4out(j1)=cm44(j2,j1)*cv4(j2)+cv4out(j1)
         enddo
      enddo
c
      return
      end
c
c************************************************************************
      subroutine mtvct(cm44,cv4,cv4out)
c************************************************************************
      implicit none
c
      complex*16 cm44(4,4),cv4(4),cv4out(4)
c
      integer j1,j2
c
      do j1=1,4
         cv4out(j1)=(0.d0,0.d0)
         do j2=1,4
            cv4out(j1)=cm44(j1,j2)*cv4(j2)+cv4out(j1)
         enddo
      enddo
c
      return
      end
c
c************************************************************************
      subroutine fldbr(hl,p,mss,spnbr)
c************************************************************************
c
      implicit none
c
c     INPUT:
c             real*8 HL spinor helicity
c             real*8 P(4), spinor momentum
c             real*8 MSS, spinor mass
c
c     OUTPUT:
c             complex*16 spnbr(4)
c
      real*8 hl,p(4),mss
c
      complex*16 spnbr(4)
c
      integer j1,j2
      real*8 ptmp(4),hl0
c
      if(p(1).lt.0.d0) then
         do j1=1,4
            ptmp(j1)=-p(j1)
         enddo
      else
         do j1=1,4
            ptmp(j1)=p(j1)
         enddo
      endif
      hl0=hl
      call fermionbarsources(hl0,spnbr,ptmp,mss)
c
      return
      end
c************************************************************************
      subroutine fld(hl,p,mss,spn)
c************************************************************************
c
      implicit none
c
c     INPUT:
c             real*8 HL spinor helicity
c             real*8 P(4), spinor momentum
c             real*8 MSS, spinor mass
c
c     OUTPUT:
c             complex*16 SPN(4)
c
      real*8 hl,p(4),mss
c
      complex*16 spn(4)
c
      real*8 ptmp(4)
      integer j1
c
      if(p(1).lt.0.d0) then
         do j1=1,4
            ptmp(j1)=-p(j1)
         enddo
      else
         do j1=1,4
            ptmp(j1)=p(j1)
         enddo
      endif
      call fermionsources(hl,spn,ptmp,mss)
c
      return
      end
c************************************************************************
      function rrlrnzprdct(rp1,rp2)
c************************************************************************
      implicit none
c
      real*8 rp1(4),rp2(4),rrlrnzprdct
c
      rrlrnzprdct=rp2(1)*rp1(1)-rp2(2)*rp1(2)-rp2(3)*rp1(3)
     .     -rp2(4)*rp1(4)
c
      return
      end
c************************************************************************
      function rrlrnzprdctrvrt(rp1,rp2)
c************************************************************************
      implicit none
c
      real*8 rp1(4),rp2(4),rrlrnzprdctrvrt
c
      rrlrnzprdctrvrt=rp2(4)*rp1(4)-rp2(2)*rp1(2)-rp2(3)*rp1(3)
     .     -rp2(1)*rp1(1)
c
      return
      end
c************************************************************************
      function crlrnzprdct(rp1,cp2)
c************************************************************************
      implicit none
c
      complex*16 crlrnzprdct,cp2(4)
      real*8 rp1(4)
c
      crlrnzprdct=cp2(1)*rp1(1)-cp2(2)*rp1(2)-cp2(3)*rp1(3)
     .     -cp2(4)*rp1(4)
c
      return
      end
c************************************************************************
      function cclrnzprdct(cp1,cp2)
c************************************************************************
      implicit none
c
      complex*16 cclrnzprdct,cp2(4),cp1(4)
c
      cclrnzprdct=cp2(1)*cp1(1)-cp2(2)*cp1(2)-cp2(3)*cp1(3)
     .     -cp2(4)*cp1(4)
c
      return
      end
c************************************************************************
      subroutine eiconal(lfb,lf,lph,pff,pphot,phel,mode,eiccrr)
c************************************************************************
c
c returns eiconal factor
c
c     INPUT:
c            integer LFB,LF,LPH: fbar, f and photon labels, for optimization
c            real*8 PFF(4,2),PPHOT(4), fbar, f and photon momentum (E,px,py,pz)
c            integer MODE, for optimization
c            real*8 PHEL, photon heilicity = +-1
c     OUTPUT:
c            complex*16 EICCRR: eiconal current
c
      implicit none
c
      integer lfb,lf,lph,mode
      real*8 pff(4,2),pphot(4),phel,rrlrnzprdct
      complex*16 eps(4),eiccrr,crlrnzprdct
c
      integer nmax
      parameter (nmax=15)      !up to 100 partons
      logical tab(2,4,4,nmax)
      complex*16 eicstr(2,4,4,nmax)
      save tab,eicstr
c
      integer j1,j2,j3,j4,idbg/0/
      save idbg
c
      if(mode.eq.0) then
         do j1=1,2
            do j2=1,4
               do j3=1,4
                  do j4=1,nmax
                     tab(j1,j2,j3,j4)=.false.
                  enddo
               enddo
            enddo
         enddo
         return
      endif
c
c     DBG
c      idbg=idbg+1
c     end DBG
c
c$$$$      if(phel.gt.0.d0) then
c$$$$         j1=1
c$$$$      else
c$$$$         j1=2
c$$$$      endif
c$$$$c
c$$$$c     DBG
c$$$$      if(tab(j1,lfb,lf,lph)) then
c$$$$c      if(.false.) then
c$$$$c     end DBG
c$$$$         eiccrr=eicstr(j1,lfb,lf,lph)
c$$$$      else
c$$$$c     DBG
c$$$$c         write(*,*)'lbl  ',lph,'  pph',pphot,idbg
c$$$$c     end DBG
c$$$$         call photpol(lph,pphot,phel,mode,eps)
c$$$$         eiccrr=crlrnzprdct(pff(1,1),eps)/rrlrnzprdct(pff(1,1),pphot)
c$$$$         eiccrr=eiccrr-crlrnzprdct(pff(1,2),eps)
c$$$$     .        /rrlrnzprdct(pff(1,2),pphot)
c$$$$c     DBG
c$$$$c         if( tab(j1,lfb,lf,lph) ) then
c$$$$c            if(abs(eiccrr - eicstr(j1,lfb,lf,lph)).gt.1.d-10) then
c$$$$c               write(*,*)'lbl  ',lph,(eiccrr - eicstr(j1,lfb,lf,lph)),
c$$$$c     .              lfb,lf,j1,' eic wr',pphot,idbg
c$$$$c            endif
c$$$$c         else
c$$$$c            write(*,*)'lbl  ',lph,eiccrr,
c$$$$c     .              lfb,lf,j1,' eic nw',pphot,idbg
c$$$$c         endif
c$$$$c     end DBG
c$$$$         eicstr(j1,lfb,lf,lph)=eiccrr
c$$$$         tab(j1,lfb,lf,lph)=.true.
c$$$$      endif
c     DBG
c      write(*,*)'eiccrr',eiccrr,' eps',eps
c      write(*,*)'pff',pff
c     end DBG
c
      return
      end
c************************************************************************
      subroutine photpol(pphot,phel,eps)
c************************************************************************
c
c     INPUT: real*8 PPHOT(4), photon momentum
c            real*8 PHEL, photon helicity
c
c     OUTPUT: complex*16 EPS(4)
c
      implicit none
c
      real*8 pphot(4),phel
      complex*16 eps(4)
c
      integer j1
      real*8 m0
c
      if(phel.gt.0) then
         j1=1
      else
         j1=2
      endif
c
      m0=0.d0
      call sourcemassboson(phel,eps,pphot,m0)
c
      return
      end
c****************************************************************
      subroutine photpolrefp(q,p,lq,lp,mode,ep)
c****************************************************************
!
      implicit none
!     Return polarizations of photon with momentum Q, ep(1:4) the polarization
!     is orthogonal to reference momentum P. Convenction P==(E,px,py,pz)
!     LQ, photon label, LP ref momentum label. MODE=0 init step
!
!     INPUT
      real*8 q(4),p(4)
      integer s,lp,lq,mode
!     OUTPUT
      complex*16 ep(8)
!
      integer nmax
      parameter (nmax=15)       !up to 100 partons
      logical tab(nmax,nmax)
      complex*16 epstr(8,nmax,nmax)
      save tab,epstr
!
      integer j1,j2
      real*8 nprpp(3),r,n(3)
!
      if(mode.eq.0) then
         do j1=1,nmax
            do j2=1,nmax
               tab(j1,j2)=.false.
            enddo
         enddo
         return
      endif
!
      if(tab(lq,lp)) then
c$$$$         do j1=1,8
c$$$$            ep(j1)=epstr(j1,lq,lp)
c$$$$         enddo
c$$$$         return
c$$$$      else
c$$$$         nprpp(1)=q(3)*p(4)-q(4)*p(3)
c$$$$         nprpp(2)=q(4)*p(2)-q(2)*p(4)
c$$$$         nprpp(3)=q(2)*p(3)-q(3)*p(2)
c$$$$         if((abs(nprpp(1))+abs(nprpp(2))+abs(nprpp(3))).lt.
c$$$$     .        (1.d-14*q(1)*p(1))) then
c$$$$            write(*,*)'reference momentum collinear to photon momentum'
c$$$$            write(*,*)'other pol used collinear approx does not work'
c$$$$            r=1.d0
c$$$$            call photpol (lq,q,r,mode,ep(1))
c$$$$            r=-1.d0
c$$$$            call photpol (lq,q,r,mode,ep(5))
c$$$$            return
c$$$$         endif
c$$$$         r=sqrt(nprpp(1)**2+nprpp(3)**2+nprpp(2)**2)
c$$$$         ep(1)=(0.d0,0.d0)
c$$$$         do j1=2,4
c$$$$            ep(j1)=nprpp(j1-1)/r
c$$$$         enddo
c$$$$         r=r*q(1)
c$$$$         n(1)=(q(3)*nprpp(3)-q(4)*nprpp(2))/r
c$$$$         n(2)=(q(4)*nprpp(1)-q(2)*nprpp(3))/r
c$$$$         n(3)=(q(2)*nprpp(2)-q(3)*nprpp(1))/r
c$$$$         ep(5)=(0.d0,0.d0)
c$$$$         do j1=6,8
c$$$$            ep(j1)=n(j1-5)
c$$$$         enddo
c$$$$         tab(lq,lp)=.true.
c$$$$         do j1=1,8
c$$$$            epstr(j1,lq,lp)=ep(j1)
c$$$$         enddo
      endif
!
      return
      end

c
c***********************************************
      function rangen(r)
c***********************************************
      implicit none
      real*8 r,rangen
!
      call random_number(r)
      rangen=r
!
      return
      end
!

C***********************************************************************
          subroutine fermionsources(spinsoura,fieldaux ,mom,mass)
C***********************************************************************
C
C This subroutine return the initialized fermion field.
C Elicity heigenstates
C
          implicit none
C
          real*8  mass         !containing the mass of the external particle
          real*8  mom(4)          !array containing the momenta of the external particle
          real*8 spinsoura    !array containing the spin of the source
          real*8 p,p3p,p3m,mp,coeffp,coeffm
C
          complex*16 fieldaux(4)  !array returning the fermion field configuration
          complex*16 im  !immaginary unity in complex representation
          complex*16 p1p
          data im/(0.,1.)/
          integer spinsour
C
C Static variables
C
          save im
C
C
          p=sqrt(mom(2)**2+mom(3)**2+mom(4)**2)
          p3p=p+mom(4)
          p3m=p-mom(4)
          mp=mass+mom(1)
          p1p=mom(2)+mom(3)*im
C
          if(abs(spinsoura).lt.0.5) spinsoura=-spinsoura
          spinsour=nint(100*spinsoura)
          if(abs(p3m).lt.1.d-10.or.abs(p3p).lt.1.d-10)then
           call fs(spinsour,fieldaux ,mom,mass)
           return
          endif
C
          coeffp=1./Sqrt(2.*p*mp*p3p)
          coeffm=1./Sqrt(2.*p*mp*p3m)
C
C "spin up" ingoing fermion
C
           if (spinsour.eq.49) then
            fieldaux(1)=coeffp*p3p*mp
            fieldaux(2)=coeffp*p1p*mp
            fieldaux(3)=coeffp*p3p*p
            fieldaux(4)=coeffp*p1p*p
           endif
C
C "spin down" ingoing fermion
C
           if (spinsour.eq.-49) then
            fieldaux(1)=coeffm*p3m*mp
            fieldaux(2)=-coeffm*p1p*mp
            fieldaux(3)=-coeffm*p3m*p
            fieldaux(4)=coeffm*p1p*p
           endif
C
C "spin up" outgoing antifermion
C
           if (spinsour.eq.51) then
            fieldaux(1)=-coeffm*p3m*p
            fieldaux(2)=coeffm*p1p*p
            fieldaux(3)=coeffm*p3m*mp
            fieldaux(4)=-coeffm*p1p*mp
           endif
C
C "spin down" outgoing antifermion
C
           if (spinsour.eq.-51) then
            fieldaux(1)=coeffp*p3p*p
            fieldaux(2)=coeffp*p1p*p
            fieldaux(3)=coeffp*p3p*mp
            fieldaux(4)=coeffp*p1p*mp
           endif
C
           return
           end

C***********************************************************************
          subroutine fermionbarsources(spinsoura,fieldaux ,mom,mass)
C***********************************************************************
C
C This subroutine return the initialized fermionbar field.
C Elicity heigenstates
C
          implicit none
C
          real*8  mass         !containing the mass of the external particle
          real*8  mom(4)          !array containing the momenta of the external particle
          real*8 spinsoura    !array containing the spin of the source
          real*8 p,p3p,p3m,mp,coeffp,coeffm
C
          complex*16 fieldaux(4)     !array returning the fermion field configuration
          complex*16 im  !immaginary unity in complex representation
          complex*16 p1p
          data im/(0.,1.)/
          integer spinsour
C
C Static variables
C
          save im
C
C
          p=sqrt(mom(2)**2+mom(3)**2+mom(4)**2)
          p3p=p+mom(4)
          p3m=p-mom(4)
          mp=mass+mom(1)
          p1p=mom(2)-mom(3)*im
C
          if(abs(spinsoura).lt.0.5) spinsoura=-spinsoura
          spinsour=nint(100*spinsoura)
          if(abs(p3m).lt.1.d-10.or.abs(p3p).lt.1.d-10)then
           call fbs(spinsour,fieldaux ,mom,mass)
           return
          endif
C
          coeffp=1./Sqrt(2.*p*mp*p3p)
          coeffm=1./Sqrt(2.*p*mp*p3m)
C
C "spin up" ingoing fermion
C
           if (spinsour.eq.49) then
            fieldaux(1)=coeffp*p3p*mp
            fieldaux(2)=coeffp*p1p*mp
            fieldaux(3)=-coeffp*p3p*p
            fieldaux(4)=-coeffp*p1p*p
           endif
C
C "spin down" ingoing fermion
C
           if (spinsour.eq.-49) then
            fieldaux(1)=coeffm*p3m*mp
            fieldaux(2)=-coeffm*p1p*mp
            fieldaux(3)=coeffm*p3m*p
            fieldaux(4)=-coeffm*p1p*p
           endif
C
C "spin up" outgoing antifermion
C
           if (spinsour.eq.51) then
            fieldaux(1)=-coeffm*p3m*p
            fieldaux(2)=coeffm*p1p*p
            fieldaux(3)=-coeffm*p3m*mp
            fieldaux(4)=coeffm*p1p*mp
           endif
C
C "spin down" outgoing antifermion
C
           if (spinsour.eq.-51) then
            fieldaux(1)=coeffp*p3p*p
            fieldaux(2)=coeffp*p1p*p
            fieldaux(3)=-coeffp*p3p*mp
            fieldaux(4)=-coeffp*p1p*mp
           endif
C
           return
           end

C***********************************************************************
          subroutine fs(spinsour,fieldaux ,mom,mass)
C***********************************************************************
C
C This subroutine return the initialized fermion field.
C
          implicit none
C
          real*8  mass         !containing the mass of the external particle
          real*8  mom(4)       !array containing the momenta of the external particle
          integer spinsour    !array containing the spin of the source
C
          complex*16 fieldaux(4)     !array returning the fermion field configuration
          complex*16 im  !immaginary unity in complex representation
          data im/(0.,1.)/
C
C Static variables
C
          save im
C
C
C "spin up" ingoing fermion
C
           if(mom(4).lt.0) spinsour=-spinsour
           if (spinsour.eq.49) then
            fieldaux(1)=sqrt((mom(1)+mass))
            fieldaux(2)=0.
            fieldaux(3)=mom(4)/sqrt((mass+mom(1)))
            fieldaux(4)=mom(2)/sqrt((mass+mom(1)))
     >                         + mom(3)/sqrt((mass+mom(1)))*im
           endif
C
C "spin down" ingoing fermion
C
           if (spinsour.eq.-49) then
            fieldaux(1)=0.
            fieldaux(2)=sqrt((mom(1)+mass))
            fieldaux(3)= mom(2)/sqrt((mass+mom(1)))
     >                         - mom(3)/sqrt((mass+mom(1)))*im
            fieldaux(4)=-mom(4)/sqrt((mass+mom(1)))
         endif
C
C "spin up" outgoing antifermion
C
           if (spinsour.eq.51) then
            fieldaux(1)=mom(4)/sqrt((mass+mom(1)))
            fieldaux(2)=mom(2)/sqrt((mass+mom(1)))
     >                         + mom(3)/sqrt((mass+mom(1)))*im
            fieldaux(3)=sqrt((mom(1)+mass))
            fieldaux(4)=0.
           endif
C
C "spin down" outgoing antifermion
C
           if (spinsour.eq.-51) then
            fieldaux(1)= mom(2)/sqrt((mass+mom(1)))
     >                         - mom(3)/sqrt((mass+mom(1)))*im
            fieldaux(2)=-mom(4)/sqrt((mass+mom(1)))
            fieldaux(3)=0.
            fieldaux(4)=sqrt((mom(1)+mass))
           endif
C
           return
           end
C***********************************************************************
          subroutine fbs(spinsour,fieldaux,mom,mass)
C***********************************************************************
C
C This subroutine return the initialized fermion field.
C
          implicit none
C
          real*8  mass        !containing the mass of the external particle
          real*8  mom(4)      !array containing the momenta of the external particle
          integer spinsour    !array containing the spin of the source
C
          complex*16 fieldaux(4)   !array returning the fermionbar field configuration
          complex*16 im           !immaginary unity in complex representation
          data im/(0.,1.)/
C
C Static variables
C
          save im
C
C
C "spin up" outgoing fermion
C
           if(mom(4).gt.0) spinsour=-spinsour
           if (spinsour.eq.49) then
            fieldaux(1)=sqrt((mom(1)+mass))
            fieldaux(2)=0.
            fieldaux(3)=-mom(4)/sqrt((mass+mom(1)))
            fieldaux(4)=-mom(2)/sqrt((mass+mom(1)))
     >                         +mom(3)/sqrt((mass+mom(1)))*im
           endif
C
C "spin down" outgoing fermion
C
           if (spinsour.eq.-49) then
            fieldaux(1)=0.
            fieldaux(2)=sqrt((mom(1)+mass))
            fieldaux(3)= -mom(2)/sqrt((mass+mom(1)))
     >                         - mom(3)/sqrt((mass+mom(1)))*im
            fieldaux(4)=mom(4)/sqrt((mass+mom(1)))
                    endif
C
C "spin up" ingoing antifermion
C
           if (spinsour.eq.51) then
            fieldaux(1)=-mom(4)/sqrt((mass+mom(1)))
            fieldaux(2)=-mom(2)/sqrt((mass+mom(1)))
     >                          +mom(3)/sqrt((mass+mom(1)))*im
            fieldaux(3)=sqrt((mom(1)+mass))
            fieldaux(4)=0.
           endif
C
C "spin down" ingoing antifermion
C
           if (spinsour.eq.-51) then
            fieldaux(1)= -mom(2)/sqrt((mass+mom(1)))
     >                         - mom(3)/sqrt((mass+mom(1)))*im
            fieldaux(2)=mom(4)/sqrt((mass+mom(1)))
            fieldaux(3)=0.
            fieldaux(4)=sqrt((mom(1)+mass))
           endif
C
           return
           end

C***************************************************************************
          subroutine sourcemassboson(spinsour,fieldaux,mom,mass)
C***************************************************************************
C
C This subroutine given the four momentum of a massive boson MOM the required
C source polarization SPINSOUR and the boson mass MASS returns in
C FIELDAUX the source term
C
          implicit none
C
          integer n               !loop index
C
          real*8 knorm2           !squared modulus of the three momentum
          real*8 mass             !massive boson mass
          real*8 massa2           !squared mass
          real*8 mom(4)           !four momentum
          real*8 nz(4),nx(4)
          real*8 sour1(4),sour2(4),sour3(4)   !longitudinal (1)
C                                and transverse (2,3) polarizations
          real*8 spinsour   !required polarization
          real*8 scalar3
          real*8 xnorm      !to normalize the sources
          complex*16 fieldaux(4)   !to return the relevan polarization
C
          data nz/0.,0.,0.,1./
          data nx/0.,1.,0.,0./
          data sour1/4*0./
          data sour2/4*0./
          data sour3/4*0./
C
          save nx,nz
C
C
C
          knorm2=mom(2)**2+mom(3)**2+mom(4)**2
C
          if (knorm2.eq.0) then
             do n=1,4
                sour1(n)=0.d0
                sour2(n)=0.d0
                sour3(n)=0.d0
             enddo
             sour1(2)=1.
             sour2(3)=1.
             sour3(4)=1.
             return
          endif
          massa2=mass**2
C
C    pol. longitudinal
C
          do n=2,4
           sour1(n)=mom(n)
          enddo
	  sour1(1)=knorm2/mom(1)
          xnorm=sqrt(scalar3(sour1,sour1))
          if(xnorm.ne.0.) then
           do n=1,4
            sour1(n)=sour1(n)/xnorm
           enddo
          endif
C
C    pol. transverse +
C
          if (mom(4)**2.gt.mom(2)**2) then
           call vector3prod(mom,nx,sour2)
          else
           call vector3prod(mom,nz,sour2)
          endif
          xnorm=Sqrt(scalar3(sour2,sour2))
          do n=1,4
           sour2(n)=sour2(n)/xnorm
          enddo
C
C    pol. transverse -
C
          call vector3prod(mom,sour2,sour3)
          xnorm=sqrt(scalar3(sour3,sour3))
          do n=1,4
           sour3(n)=sour3(n)/xnorm
          enddo
C
          do n=1,4
           fieldaux(n)=0.
           if(nint(spinsour).eq.0)fieldaux(n)=sour1(n)
           if(nint(spinsour).eq.1)fieldaux(n)=sour2(n)
           if(nint(spinsour).eq.-1)fieldaux(n)=sour3(n)
          enddo
C
          return
          end
C*********************************************************************
	  subroutine vector3prod(v1,v2,vfin)
C*********************************************************************
C
          implicit none
C
           real*8 v1(4)
           real*8 v2(4)
           real*8 vfin(4)
C
           vfin(2)=v1(3)*v2(4)-v1(4)*v2(3)
           vfin(3)=v1(4)*v2(2)-v1(2)*v2(4)
           vfin(4)=v1(2)*v2(3)-v1(3)*v2(2)
           vfin(1)=0
C
           return
           end
C***********************************************************************
      function scalar3(v1,v2)
C***********************************************************************
C
      implicit none
C
      real*8 scalar3
      real*8 v1(4)
      real*8 v2(4)
C
      scalar3=abs(v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)-v1(4)*v2(4))
C
      return
      end
****  e+ e- pi pi Ng
      subroutine eepipi(p,jspn,optspn,mss,nph,pph,mesq)
C***********************************************************************
!
!     Input: P(i,j), fermion/pion momentum; i=E,px,py,pz; j particle label
!            j=1,_4 = e-,e+,pi+,pi-; incoming momenta (e+-) have REVERSED
!            sign (E<0, pj->-pj)
!            JSPN=0,2**(NFER+NPH)-1 label to select helicity configuration
!            OPTSPN, for optimization over helicity configuration. If .TRUE.
!              then a hit or miss is performed and if fails MESQ=0 is returned
!              if .FALSE. all helicities come with weight=1
!            MSS(i) fermion masses (i=1,2) or pion masses (i=3,4)
!            NPH, photons number, NPH <=NPHMX which is a parameter in PRM.INC
!            PPH(i,j); photon momentum of the j-th photon
!     Output: MESQ, squared matrix element.
!             NB alphaem = 1/(4 pi) (electric charge = 1)
!
      implicit none
!
      include 'prm.f'
!
      integer nfer
      parameter (nfer=2)
      integer npi
      parameter (npi=2)
      integer nph
      real*8 p(4,nfer+npi),pph(4,nphmx),mss(nfer+npi)
      integer jspn
!
      logical OPTSPN
      real*8 mesq
!
      complex*16 gam(4,4,5),iden(4,4) !N.B (gamma_\mu)_{ij} ==
      common/gamma/gam,iden                       ! gam(j,i,mu)
!
      real*8 spn(nfer),epsph(nphmx)
!            SPN(i); +-0.49 fermions; +-0.51 antifermions
!            EPSPH(i); = +-1.d0 polarization of the jth photon
c      complex*16 psi(4,0:2**nphmx-1,nfer),pi(0:2**nphmx-1,npi)

      complex*16, allocatable :: psi(:,:,:),pi(:,:)

      integer tblmrg(nphmx,0:2**nphmx-1)
      complex*16 curree(4,0:2**nphmx-1),currpipi(4,0:2**nphmx-1),amp
      complex*16 cclrnzprdct,photonpropagatorpipi,eps(4,nphmx)
      real*8 rrlrnzprdct,pin1(4),pin2(4),pout1(4),pout2(4)
      logical frst/.true./,mrg
      integer j1,j2,j3,jj,j5
      save frst,tblmrg
      real*8 q(4),s1,s2,pspne,pspnm,prb,r
!

      
      allocate(psi(4, 0:2**nphmx-1, nfer))
      allocate(pi(0:2**nphmx-1, npi))

      jj=jspn
      spn(1)=0.49d0
      if (mod(jj,2).eq.0) spn(1)=-0.49d0
      s1=sign(1.d0,spn(1))
      jj=jj/2
      spn(2)=0.51d0
      if (mod(jj,2).eq.0) spn(2)=-0.51d0
      s1=sign(1.d0,spn(2))
      jj=jj/2
      do j1=1,nph
         epsph(j1)=1.d0
         if (mod(jj,2).eq.0) epsph(j1)=-1.d0
         jj=jj/2
      enddo
!     test
!      epsph(nph)=-1.d30
!     end test

      if(s1.gt.0.d0) then
         pspne=1.d0
      else
         pspne=(mss(1)/min(p(1,1),p(1,3)))**2*1.d1
         pspne=max(0.1d0,pspne)
         pspne=min(1.d0,pspne)
      endif
      prb=pspne
      if (optspn) then
         call rndmgnrtr(r)
         if(r.lt.prb) then
            prb=1.d0/prb
         else
            mesq=0.d0
         endif
      else
         prb=1.d0
      endif
!
!      write(*,*)'spn  ',spn,'  eps',epsph(1:nph)
      call absphotonfld(p(1:4,1:2),nfer,spn,mss(1:2),nph,pph,epsph,psi)
      call absphotonpion(p(1:4,3:4),npi,mss(3:4),nph,pph,epsph,pi)
!
      if(frst) then
         frst=.false.
         call mrgng(tblmrg)
      endif
      curree=(0.d0,0.d0)
      currpipi=(0.d0,0.d0)
      do j1=1,nph
         call photpol(pph(1:4,j1),epsph(j1),eps(1:4,j1))
      enddo
      
      do j1=0,2**nph-1
         do j2=j1,2**nph-1
            if(j2.eq.j1.and.j1.ne.0) cycle
!         do j2=0,2**nph-1
            mrg=.true.
            do j3=1,nph
               mrg=mrg.and.(tblmrg(j3,j1)+tblmrg(j3,j2).le.1)
            enddo
            if(mrg) then
               jj=j1+j2
               pout1=p(1:4,4)
               pin1=-p(1:4,3)
               pout2=p(1:4,4)
               pin2=-p(1:4,3)
               do j5=1,nph
                  if(tblmrg(j5,j1).eq.1) then
                     pout1=pout1+pph(1:4,j5)
                     pin2=pin2-pph(1:4,j5)
                  endif
                  if(tblmrg(j5,j2).eq.1) then
                     pout2=pout2+pph(1:4,j5)
                     pin1=pin1-pph(1:4,j5)
                  endif
               enddo
               do j3=1,4
                  curree(j3,jj)=dot_product(conjg(psi(1:4,j2,2)),
     .                 matmul(psi(1:4,j1,1),gam(1:4,1:4,j3)) )
     .                 +curree(j3,jj)  !  gamma matrix are transposed
                  if(j1.ne.j2) curree(j3,jj)=
     .                 dot_product(conjg(psi(1:4,j1,2)),
     .                 matmul(psi(1:4,j2,1),gam(1:4,1:4,j3)) )
     .                 +curree(j3,jj)  !  gamma matrix are transposed
!
                  currpipi(j3,jj)=currpipi(j3,jj)
     .                 +(pout1(j3)+pin1(j3))*pi(j1,2)*pi(j2,1)
                  if(j1.ne.j2) currpipi(j3,jj)=currpipi(j3,jj)
     .                 +(pout2(j3)+pin2(j3))*pi(j2,2)*pi(j1,1) ! trilinear
                  do j5=1,nph
!                     cycle
                     if(tblmrg(j5,j1).eq.1) then
                        currpipi(j3,jj)=currpipi(j3,jj)-
     .                       pi(j2,1)*pi(j1-2**(j5-1),2)*eps(j3,j5) ! quadrilinear
!                        if(j1-2**(j5-1).ne.j2)
!     .                       currpipi(j3,jj)=currpipi(j3,jj)-
!     .                       pi(j2,2)*pi(j1-2**(j5-1),1)*eps(j3,j5) ! quadrilinear
                        currpipi(j3,jj)=currpipi(j3,jj)-
     .                       pi(j2,2)*pi(j1-2**(j5-1),1)*eps(j3,j5) ! quadrilinear
                     endif
                     if(tblmrg(j5,j2).eq.1) then
                        currpipi(j3,jj)=currpipi(j3,jj)-
     .                       pi(j2-2**(j5-1),1)*pi(j1,2)*eps(j3,j5) ! quadrilinear
!                        if(j2-2**(j5-1).ne.j1)
!     .                       currpipi(j3,jj)=currpipi(j3,jj)-
!     .                       pi(j2-2**(j5-1),2)*pi(j1,1)*eps(j3,j5) ! quadrilinear
                        currpipi(j3,jj)=currpipi(j3,jj)-
     .                       pi(j2-2**(j5-1),2)*pi(j1,1)*eps(j3,j5) ! quadrilinear
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo
      amp=0.d0
      jj=2**nph-1
      do j1=0,2**nph-1
!     test
!!!         write(*,*)'j1, crre',j1, curree(1:4,j1)
!!!         write(*,*)'j1, crrm',j1, currmumu(1:4,j1)
!!!         write(*,*)'***'
!     test
         q=p(1:4,1)+p(1:4,2)
         do j2=1,nph
            if(tblmrg(j2,j1).eq.1) q=q+pph(1:4,j2)
         enddo
         amp=amp+cclrnzprdct(curree(1:4,j1),currpipi(1:4,jj-j1))
     .        *photonpropagatorpipi(q)
!         write(*,*)'xx', cclrnzprdct(curree(1:4,j1),currmumu(1:4,jj-j1))
!     .        **photonpropagatorpipi(q),'  amp',amp
      enddo
      mesq=abs(amp)**2*prb
!     test
!!!$$$      do j1=0,2**nph-1
!!!$$$         q=p(1:4,1)+p(1:4,3)
!!!$$$         do j2=1,nph
!!!$$$            if(tblmrg(j2,j1).eq.1) q=q+pph(1:4,j2)
!!!$$$         enddo
!!!$$$         write(*,*)'qmu jmu', q(1)*curree(1,j1)-q(2)*curree(2,j1)
!!!$$$     .        -q(3)*curree(3,j1)-q(4)*curree(4,j1)
!!!$$$         write(*,*)'j1,jmu',j1,curree(1:4,j1)
!!!$$$         write(*,*)'******'
!!!$$$      enddo
!
      do j1=0,2**nph-1
         q=p(1:4,3)+p(1:4,4)
         do j2=1,nph
            if(tblmrg(j2,j1).eq.1) q=q+pph(1:4,j2)
         enddo
!!!$$$        write(*,*)'qmu jpi', q(1)*currpipi(1,j1)-q(2)*currpipi(2,j1)
!!!$$$    .        -q(3)*currpipi(3,j1)-q(4)*currpipi(4,j1)
!!!$$$        write(*,*)'j1,jpipi',j1,currpipi(1:4,j1)
!!!$$$        write(*,*)'******'
      enddo
!

      deallocate(psi)
      deallocate(pi)

      return
      end
C***********************************************************************
      subroutine absphotonpion(ppion,npi,mss,nphot,pphot,epsphot,pi)
C***********************************************************************
      implicit none
c
c     Input: ppion(i,j), pion momentum, i == (E,px,py,pz), initial state pions E<0 pj -> -pj,
c                       final state pions E>0 , j= j-th pion  pi+ field before
c                       pi- fields (psibar) namely initial state anti-particle
c                       and final state particle first
c            npi, number of pions
c            mss, pion mass
c            nphot, photon number
c            pphot(j,k), photon momentum, j=1,_,4 (E,px,py,pz), k= kth photon
c            epsphot(j), helicity of j-th photon
c     Output: pi(i2,i3), field after the absorption of photons in any
c             possible combination. i2 label absorbed photons
c             i3 labels absorbing pion
c
      include 'prm.f'
c
      integer nphot,npi
      real*8 ppion(4,npi),pphot(4,nphmx),epsphot(nphmx)
      real*8 mss(npi)
c
      complex*16 pi(0:2**nphmx-1,npi)
c
      integer tblmrg(nphmx,0:2**nphmx-1)
      logical frst/.true./,tabnw(0:2**nphmx)
      save tblmrg,frst
c
      integer j1,j2,jlow,jhigh,jj,j3,jj0,idabs(2**nphmx),j4,j5
      integer idabsnw(2**nphmx),idabsq(2**nphmx),jhighq,jhighqnw
      integer idabsqnw(2**nphmx)
      real*8 pin(4),kk(4),den,rrlrnzprdct,kksum(4,2**nphmx),pout(4)
      real*8 piprp(2**nphmx,npi),ptot(4)
      complex*16 eps(4,nphmx),crlrnzprdct,cclrnzprdct
      complex*16 pitmp
c
      if (frst) then
         frst=.false.
         call mrgng(tblmrg)
      endif
c
      do j1=1,nphot
         call photpol(pphot(1,j1),epsphot(j1),eps(1,j1))
!     test
!!!         if(j1.eq.1) epstmp=pphot(1:4,1)
!     end test
      enddo
c
      do j4=1,npi
         do j1=1,2**nphot-1
            if(j4.eq.1) then
               kk=0.d0
               do j2=1,nphot
                  if(tblmrg(j2,j1).eq.1) then
                     kk=kk+pphot(1:4,j2)
                  endif
               enddo
               kksum(1:4,j1)=kk
            endif
            den=2.d0*rrlrnzprdct(ppion(1,j4),kksum(1,j1))
            den=den+rrlrnzprdct(kksum(1,j1),kksum(1,j1))
            den=1.d0/den
            piprp(j1,j4)=den
         enddo
      enddo
c
      do j5=1,npi
!     end test
         pi(0,j5)=(1.d0,0.d0)
         jhigh=1
         jhighq=0
         idabs(1)=0
         tabnw=.true.
         do j4=1,nphot
            jj=0
            do j1=1,jhigh
               do j2=1,nphot
                  if(tblmrg(j2,idabs(j1)).eq.0) then
                     jj0=2**(j2-1)+idabs(j1)
                     if(tabnw(jj0)) then
                        pi(jj0,j5)=(0.d0,0.d0)
                        tabnw(jj0)=.false.
                        jj=jj+1
                        idabsnw(jj)=jj0
                     endif
                     if(j5.gt.npi/2) then
                        pout=ppion(1:4,j5)+kksum(1:4,jj0)
                        pin=pout-pphot(1:4,j2)
                        ptot=pin+pout
                        pitmp=pi(idabs(j1),j5)*piprp(jj0,j5)*
     .                  crlrnzprdct(ptot,eps(1:4,j2))
                     else
                        pout=-ppion(1:4,j5)-kksum(1:4,jj0)
                        pin=pout+pphot(1:4,j2)
                        ptot=pin+pout
                        pitmp=pi(idabs(j1),j5)*piprp(jj0,j5)*
     .                  crlrnzprdct(ptot,eps(1:4,j2))
                     endif
                     pi(jj0,j5)=pi(jj0,j5)+pitmp
                  endif
               enddo
            enddo
            jhighqnw=jhigh
            jhigh=jj
            idabsqnw=idabs
            idabs=idabsnw
            do j1=1,jhighq
               do j2=1,nphot
                  do j3=j2+1,nphot
                     if((tblmrg(j2,idabsq(j1))).eq.0.and.
     .                    (tblmrg(j3,idabsq(j1)).eq.0)) then
                        jj0=idabsq(j1)+2**(j2-1)+2**(j3-1)
                        pitmp=pi(idabsq(j1),j5)*piprp(jj0,j5)*2.d0
                        pitmp=pitmp*cclrnzprdct(eps(1:4,j2),eps(1:4,j3))
                        pi(jj0,j5)=pi(jj0,j5)-pitmp
                     endif
                  enddo
               enddo
            enddo
            jhighq=jhighqnw
            idabsq=idabsqnw
         enddo
!
      enddo
c
      return
      end
C***********************************************************************


C***********************************************************************
**** bhabha and e+e- -> gg
      subroutine eeee(p,jspn,optspn,mss,nph,pph,mesq)
C***********************************************************************
!
!     Input: P(i,j), fermion momentum; i=E,px,py,pz; j fermion label
!            j=1,_4 = e-,e+, IN; e+,e-; OUT incoming momenta (e+-) have REVERSED
!            sign (E<0, pj->-pj)
!            JSPN=0,2**(NFER+NPH)-1 label to select helicity configuration
!            OPTSPN, for optimization over helicity configuration. If .TRUE.
!              then a hit or miss is performed and if fails MESQ=0 is returned
!              if .FALSE. all helicities come with weight=1
!            MSS(i) fermion masses
!            NPH, photons number, NPH <=NPHMX which is a parameter in PRM.INC
!            PPH(i,j); photon momentum of the j-th photon
!     Outpun: MESQ, squared matrix element. NB alphaem = 1/(4 pi)
!
      implicit none
!
      include 'prm.f'
!
      integer nfer
      parameter (nfer=4)
      integer nph
      real*8 p(4,nfer),pph(4,nphmx),mss(nfer)
      integer jspn
!
      logical OPTSPN
      real*8 mesq
!
      complex*16 gam(4,4,5),iden(4,4) !N.B (gamma_\mu)_{ij} ==
      common/gamma/gam,iden                       ! gam(j,i,mu)
!
      real*8 spn(nfer),epsph(nphmx)
!            SPN(i); +-0.49 fermions; +-0.51 antifermions
!            EPSPH(i); = +-1.d0 polarization of the jth photon
      complex*16 psi(4,0:2**nphmx-1,nfer)
      integer tblmrg(nphmx,0:2**nphmx-1)
      complex*16 curree13(4,0:2**nphmx-1),curree24(4,0:2**nphmx-1),amp
      complex*16 curree14(4,0:2**nphmx-1),curree23(4,0:2**nphmx-1)
      complex*16 cclrnzprdct,photonpropagator
      real*8 rrlrnzprdct
      logical frst/.true./,mrg
      integer j1,j2,j3,jj
      save frst,tblmrg
      real*8 q(4),s1,s2,pspne,pspnm,prb,r
!
      jj=jspn
      spn(1)=0.49d0
      if (mod(jj,2).eq.0) spn(1)=-0.49d0
      s1=sign(1.d0,spn(1))
      jj=jj/2
      spn(2)=0.51d0
      if (mod(jj,2).eq.0) spn(2)=-0.51d0
      s2=sign(1.d0,spn(2))
      jj=jj/2
      spn(3)=0.51d0
      if (mod(jj,2).eq.0) spn(3)=-0.51d0
      jj=jj/2
      s1=s1*sign(1.d0,spn(3))
      spn(4)=0.49d0
      if (mod(jj,2).eq.0) spn(4)=-0.49d0
      s2=s2*sign(1.d0,spn(4))
      jj=jj/2
      do j1=1,nph
         epsph(j1)=1.d0
         if (mod(jj,2).eq.0) epsph(j1)=-1.d0
         jj=jj/2
      enddo
      if(s1.gt.0.d0) then
         pspne=1.d0
      else
         pspne=(mss(1)/min(p(1,1),p(1,3)))**2*1.d1
         pspne=max(0.1d0,pspne)
         pspne=min(1.d0,pspne)
      endif
      if(s2.gt.0.d0) then
         pspnm=1.d0
      else
         pspnm=(mss(2)/min(p(1,2),p(1,4)))**2*1.d1
         pspnm=max(0.1d0,pspnm)
         pspnm=min(1.d0,pspnm)
      endif
      prb=pspne*pspnm
      call rndmgnrtr(r)
      if (optspn) then
         if(r.lt.prb) then
            prb=1.d0/prb
         else
            mesq=0.d0
         endif
      else
         prb=1.d0
      endif
!
!      write(*,*)'spn  ',spn,'  eps',epsph(1:nph)
      call absphotonfld(p,nfer,spn,mss,nph,pph,epsph,psi)
!
      if(frst) then
         frst=.false.
         call mrgng(tblmrg)
      endif
      curree13=(0.d0,0.d0)
      curree24=(0.d0,0.d0)
      curree23=(0.d0,0.d0)
      curree14=(0.d0,0.d0)
      do j1=0,2**nph-1
         do j2=j1,2**nph-1
            if(j2.eq.j1.and.j1.ne.0) cycle
!         do j2=0,2**nph-1
            mrg=.true.
            do j3=1,nph
               mrg=mrg.and.(tblmrg(j3,j1)+tblmrg(j3,j2).le.1)
            enddo
            if(mrg) then
               jj=j1+j2
               do j3=1,4
                  curree13(j3,jj)=dot_product(conjg(psi(1:4,j2,3)),
     .                 matmul(psi(1:4,j1,1),gam(1:4,1:4,j3)) )
     .                 +curree13(j3,jj)  !  gamma matrix are transposed
                  if(j1.ne.j2) curree13(j3,jj)=
     .                 dot_product(conjg(psi(1:4,j1,3)),
     .                 matmul(psi(1:4,j2,1),gam(1:4,1:4,j3)) )
     .                 +curree13(j3,jj)  !  gamma matrix are transposed
                  curree24(j3,jj)=dot_product(conjg(psi(1:4,j1,4)),
     .                 matmul(psi(1:4,j2,2),gam(1:4,1:4,j3)) )
     .                 +curree24(j3,jj)  !  gamma matrix are transposed
                  if(j1.ne.j2) curree24(j3,jj)=
     .                 dot_product(conjg(psi(1:4,j2,4)),
     .                 matmul(psi(1:4,j1,2),gam(1:4,1:4,j3)) )
     .                 +curree24(j3,jj)  !  gamma matrix are transposed
!
                  curree14(j3,jj)=dot_product(conjg(psi(1:4,j2,4)),
     .                 matmul(psi(1:4,j1,1),gam(1:4,1:4,j3)) )
     .                 +curree14(j3,jj)  !  gamma matrix are transposed
                  if(j1.ne.j2) curree14(j3,jj)=
     .                 dot_product(conjg(psi(1:4,j1,4)),
     .                 matmul(psi(1:4,j2,1),gam(1:4,1:4,j3)) )
     .                 +curree14(j3,jj)  !  gamma matrix are transposed
                  curree23(j3,jj)=dot_product(conjg(psi(1:4,j1,3)),
     .                 matmul(psi(1:4,j2,2),gam(1:4,1:4,j3)) )
     .                 +curree23(j3,jj)  !  gamma matrix are transposed
                  if(j1.ne.j2) curree23(j3,jj)=
     .                 dot_product(conjg(psi(1:4,j2,3)),
     .                 matmul(psi(1:4,j1,2),gam(1:4,1:4,j3)) )
     .                 +curree23(j3,jj)  !  gamma matrix are transposed
               enddo
            endif
         enddo
      enddo
      amp=0.d0
      jj=2**nph-1
      do j1=0,2**nph-1
!     test
!!!         write(*,*)'j1, crre',j1, curree(1:4,j1)
!!!         write(*,*)'j1, crrm',j1, currmumu(1:4,j1)
!!!         write(*,*)'***'
!     test
         q=p(1:4,1)+p(1:4,3)
         do j2=1,nph
            if(tblmrg(j2,j1).eq.1) q=q+pph(1:4,j2)
         enddo
         amp=amp+cclrnzprdct(curree13(1:4,j1),curree24(1:4,jj-j1))
     .        *photonpropagator(q)
         q=p(1:4,1)+p(1:4,4)
         do j2=1,nph
            if(tblmrg(j2,j1).eq.1) q=q+pph(1:4,j2)
         enddo
         amp=amp-cclrnzprdct(curree14(1:4,j1),curree23(1:4,jj-j1))
     .        *photonpropagator(q)
!         write(*,*)'xx', cclrnzprdct(curree(1:4,j1),currmumu(1:4,jj-j1))
!     .        **photonpropagator(q),'  amp',amp
      enddo
      mesq=abs(amp)**2*prb
!     test
!!!$$$      do j1=0,2**nph-1
!!!$$$         q=p(1:4,1)+p(1:4,3)
!!!$$$         do j2=1,nph
!!!$$$            if(tblmrg(j2,j1).eq.1) q=q+pph(1:4,j2)
!!!$$$         enddo
!!!$$$         write(*,*)'qmu jmu', q(1)*curree(1,j1)-q(2)*curree(2,j1)
!!!$$$     .        -q(3)*curree(3,j1)-q(4)*curree(4,j1)
!!!$$$         write(*,*)'j1,jmu',j1,curree(1:4,j1)
!!!$$$         write(*,*)'******'
!!!$$$      enddo
!
!!!$$$      do j1=0,2**nph-1
!!!$$$         q=p(1:4,2)+p(1:4,4)
!!!$$$         do j2=1,nph
!!!$$$            if(tblmrg(j2,j1).eq.1) q=q+pph(1:4,j2)
!!!$$$         enddo
!!!$$$         write(*,*)'qmu jmu', q(1)*currmumu(1,j1)-q(2)*currmumu(2,j1)
!!!$$$     .        -q(3)*currmumu(3,j1)-q(4)*currmumu(4,j1)
!!!$$$         write(*,*)'j1,jmu',j1,currmumu(1:4,j1)
!!!$$$         write(*,*)'******'
!!!$$$      enddo
!
      return
      end
C***********************************************************************
      subroutine eegmgm(p,jspn,optspn,mss,nph,pph,mesq)
C***********************************************************************
!
!     Input: P(i,j), fermion momentum; i=E,px,py,pz; j fermion label
!            j=1,2 = e-,e+; incoming momenta (e+-) have REVERSED
!            sign (E<0, pj->-pj)
!            JSPN=0,2**(NFER+NPH)-1 label to select helicity configuration
!            OPTSPN, for optimization over helicity configuration. If .TRUE.
!              then a hit or miss is performed and if fails MESQ=0 is returned
!              if .FALSE. all helicities come with weight=1
!            MSS(i) fermion masses
!            NPH, photons number, NPH <=NPHMX which is a parameter in PRM.INC
!            PPH(i,j); photon momentum of the j-th photon
!     Outpun: MESQ, squared matrix element. NB alphaem = 1/(4 pi) (Qelectron=1)
!
      implicit none
!
      include 'prm.f'
!
      integer nfer
      parameter (nfer=2)
      integer nph
      real*8 p(4,nfer),pph(4,nphmx),mss(nfer)
      integer jspn
!
      logical OPTSPN
      real*8 mesq
!
      complex*16 gam(4,4,5),iden(4,4) !N.B (gamma_\mu)_{ij} ==
      common/gamma/gam,iden                       ! gam(j,i,mu)
!
      real*8 spn(nfer),epsph(nphmx)
!            SPN(i); +-0.49 fermions; +-0.51 antifermions
!            EPSPH(i); = +-1.d0 polarization of the jth photon
      complex*16 psi(4,0:2**nphmx-1,1),psi1(4)
      integer tblmrg(nphmx,0:2**nphmx-1)
      complex*16 curree(4,0:2**nphmx-1),currmumu(4,0:2**nphmx-1),amp
      complex*16 cclrnzprdct,photonpropagator
      real*8 rrlrnzprdct
      logical frst/.true./,mrg
      integer j1,j2,j3,jj
      save frst,tblmrg
      real*8 q(4),s1,s2,pspne,pspnm,prb,r
!
      jj=jspn
      spn(1)=0.49d0
      if (mod(jj,2).eq.0) spn(1)=-0.49d0
      s1=sign(1.d0,spn(1))
      jj=jj/2
      spn(2)=0.51d0
      if (mod(jj,2).eq.0) spn(2)=-0.51d0
      s1=s1*sign(1.d0,spn(2))
      jj=jj/2
      do j1=1,nph
         epsph(j1)=1.d0
         if (mod(jj,2).eq.0) epsph(j1)=-1.d0
         jj=jj/2
      enddo
      if(s1.gt.0.d0) then
         pspne=1.d0
      else
         pspne=(mss(1)/min(p(1,1),p(1,2)))**2*1.d1
         pspne=max(0.1d0,pspne)
         pspne=min(1.d0,pspne)
      endif
      prb=pspne
      call rndmgnrtr(r)
      if (optspn) then
         if(r.lt.prb) then
            prb=1.d0/prb
         else
            mesq=0.d0
         endif
      else
         prb=1.d0
      endif
!
!      write(*,*)'spn  ',spn,'  eps',epsph(1:nph)
      call absphotonfld(p(1:4,2),nfer-1,spn(2),mss(2),nph,pph,epsph,psi)
      call fld(spn(1),p(1:4,1),mss(1),psi1)
      amp=dot_product(conjg(psi1),psi(1:4,2**(nph)-1,1))
      mesq=abs(amp)**2
!
      return
      end



