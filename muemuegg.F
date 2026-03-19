      double precision function muemuegg(p1,p2,p3,p4,k1,k2)
! from MESMER
      implicit double precision (a-h,o-z)
      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
**************************************************************
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,tm1
      double precision p1k,p2k,p3k,p4k,s12,s34,t24,t13,u14,u23
      double precision p1mk2,p2mk2,p3pk2,p4pk2
      common/invariants/p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,tm1
      common/invariants1g/p1k,p2k,p3k,p4k,s12,s34,t24,t13,u14,u23,
     .     p1mk2,p2mk2,p3pk2,p4pk2

      double complex cuno,im,czero
      parameter (cuno  = (1.d0,0.d0))
      parameter (im    = (0.d0,1.d0))
      parameter (czero = (0.d0,0.d0))
      double precision dot
      external dot

      double precision ec,ec2,ec4,mm,mm2,me,me2,mm4,me4,mm6,me6
      common/elchmasses/ec,ec2,ec4,mm,mm2,me,me2,mm4,me4,mm6,me6
****************************************************************************
      double precision k1(0:3),k2(0:3),k1k2,p2mp4(0:3)
      double precision p1k1,p2k1,p3k1,p4k1
      double precision p1k2,p2k2,p3k2,p4k2
      double precision p1k1m1,p2k1m1,p3k1m1,p4k1m1
      double precision p1k2m1,p2k2m1,p3k2m1,p4k2m1
      integer iarun
      common/ialpharunning/iarun

      double precision tm100,tm110,tm102,tm112

      double precision el8
      integer ifirst
      common/ifirstmuemuegg/el8,ifirst
      data ifirst /0/

!      double precision vp13,vp24 already in invariants
      double complex vp13,vp24,vp24d,vp13d
      common/vp1gand1g1L/vp13,vp24,vp13d,vp24d
      double complex vp13k1,vp24k1
      double complex vp13k2,vp24k2
      double complex vp13k1k2,vp24k1k2
      double complex vp13k1d,vp24k1d
      double complex vp13k2d,vp24k2d
      double complex vp13k1k2d,vp24k1k2d
      double complex xxx
      double precision QRe,QRmu
      common/QReQRmu/QRe,QRmu
      data QRe,QRmu /2*1.d0/
      common/vps2gtree/         !already in invariants and a common there vp13,vp24,
     .     vp13k1,vp24k1,
     .     vp13k2,vp24k2,
     .     vp13k1k2,vp24k1k2,
     .     vp13k1d,vp24k1d,
     .     vp13k2d,vp24k2d,
     .     vp13k1k2d,vp24k1k2d

      double complex vpolc

      if (ifirst.eq.0) then

c         QRe  = 1d0
c         QRmu = 1d0

         ec2 = 4.d0*pi*alpha
         ec  = sqrt(ec2)
         ec4 = ec2*ec2
         mm = ammu
         me = ame
         mm2 = mm*mm
         mm4 = mm2*mm2
         mm6 = mm4*mm2
         me2 = me*me
         me4 = me2*me2
         me6 = me4*me2

         el8 = ec4*ec4
         ifirst = 1

      endif

      call loadinvariants(p1,p2,p3,p4)
      p1k1 = dot(p1,k1)
      p2k1 = dot(p2,k1)
      p3k1 = dot(p3,k1)
      p4k1 = dot(p4,k1)
      p1k2 = dot(p1,k2)
      p2k2 = dot(p2,k2)
      p3k2 = dot(p3,k2)
      p4k2 = dot(p4,k2)
      k1k2 = dot(k1,k2)

      p1k1m1 = 1.d0/p1k1
      p2k1m1 = 1.d0/p2k1
      p3k1m1 = 1.d0/p3k1
      p4k1m1 = 1.d0/p4k1
      p1k2m1 = 1.d0/p1k2
      p2k2m1 = 1.d0/p2k2
      p3k2m1 = 1.d0/p3k2
      p4k2m1 = 1.d0/p4k2

      p2mp4 = p2 - p4

      tm100 = 1.d0/dot(p2mp4,p2mp4)
      tm110 = 1.d0/dot(p2mp4-k1,p2mp4-k1)
      tm102 = 1.d0/dot(p2mp4-k2,p2mp4-k2)
      tm112 = 1.d0/dot(p2mp4-k1-k2,p2mp4-k1-k2)

***   VP STUFF
      vp13     = 1.d0
      vp24     = 1.d0
      vp13k1   = 1.d0
      vp24k1   = 1.d0
      vp13k2   = 1.d0
      vp24k2   = 1.d0
      vp13k1k2 = 1.d0
      vp24k1k2 = 1.d0

      if (iarun.eq.1) then
         vp24     = vpolc(1.d0/tm100)
         vp24k1   = vpolc(1.d0/tm110)
         vp24k2   = vpolc(1.d0/tm102)
         vp24k1k2 = vpolc(1.d0/tm112)
      endif

      vp13d     = dconjg(vp13)
      vp24d     = dconjg(vp24)
      vp13k1d   = dconjg(vp13k1)
      vp24k1d   = dconjg(vp24k1)
      vp13k2d   = dconjg(vp13k2)
      vp24k2d   = dconjg(vp24k2)
      vp13k1k2d = dconjg(vp13k1k2)
      vp24k1k2d = dconjg(vp24k1k2)

********
c      print*,'Reenable includes in muemuegg.F!!! stopping'
c      stop
      muemuegg = 0.d0

cc always the matrix element for mu+ used, by crossing symmetry in matrix_model.F
      include 'muemuegg-plus.f'

!! this part is commented in the file above (NOT COMMENTED ANYMORE!!)
c      muemuegg =
c     &     tm112**2*abs(vp24k1k2)**2*Z7_
c     &     + tm102*tm112*vp24k2*vp24k1k2*Z12_
c     &     + tm102**2*abs(vp24k2)**2*Z5_
c     &     + tm110*tm112*vp24k1*vp24k1k2*Z25_
c     &     + tm110*tm102*vp24k1*vp24k2*Z13_
c     &     + tm110**2*abs(vp24k1)**2*Z8_
c     &     + tm100*tm112*vp24*vp24k1k2*Z11_
c     &     + tm100*tm102*vp24*vp24k2*Z9_
c     &     + tm100*tm110*vp24*vp24k1*Z6_
c     &     + tm100**2*abs(vp24)**2*Z10_

ccccccccccccc
ccc not used to speedup compilation time, in matrix_model.F crossing symmetry is used
c      if (Qmu.lt.0.d0) then
c         include 'muemuegg-minus.f'
c      else
c         include 'muemuegg-plus.f'
c     endif
cccccccccccccccccccccccccccc

      muemuegg = el8 * muemuegg
      return
      end
***************************************************************
      double precision function eeeegg(p1,p2,p3,p4,k1,k2)
! derived from MESMER. It must be called with p2,p1,p4,p3,k1,k2!!
      implicit double precision (a-h,o-z)
      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
**************************************************************
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,tm1
      double precision p1k,p2k,p3k,p4k,s12,s34,t24,t13,u14,u23
      double precision p1mk2,p2mk2,p3pk2,p4pk2
      common/invariants/p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,tm1
      common/invariants1g/p1k,p2k,p3k,p4k,s12,s34,t24,t13,u14,u23,
     .     p1mk2,p2mk2,p3pk2,p4pk2

      double complex cuno,im,czero
      parameter (cuno  = (1.d0,0.d0))
      parameter (im    = (0.d0,1.d0))
      parameter (czero = (0.d0,0.d0))
      double precision dot
      external dot

      double precision ec,ec2,ec4,mm,mm2,me,me2,mm4,me4,mm6,me6
      common/elchmasses/ec,ec2,ec4,mm,mm2,me,me2,mm4,me4,mm6,me6
****************************************************************************
      double precision k1(0:3),k2(0:3),k1k2,p2mp4(0:3),p3pp4(0:3)
      double precision p1k1,p2k1,p3k1,p4k1
      double precision p1k2,p2k2,p3k2,p4k2
      double precision p1k1m1,p2k1m1,p3k1m1,p4k1m1
      double precision p1k2m1,p2k2m1,p3k2m1,p4k2m1
      integer iarun
      common/ialpharunning/iarun

      double complex tm100,tm110,tm102,tm112
      double complex sm100,sm110,sm102,sm112
      double complex tm100s,tm110s,tm102s,tm112s
      double complex sm100s,sm110s,sm102s,sm112s

      double precision el8
      integer ifirst
      common/ifirsteeeegg/el8,ifirst
      data ifirst /0/

!      double precision vp13,vp24 already in invariants
      double complex vp13,vp24,vp24d,vp13d
      common/vp1gand1g1L/vp13,vp24,vp13d,vp24d
      double complex vp13k1,vp24k1
      double complex vp13k2,vp24k2
      double complex vp13k1k2,vp24k1k2
      double complex vp13k1d,vp24k1d
      double complex vp13k2d,vp24k2d
      double complex vp13k1k2d,vp24k1k2d
      double complex xxx
      double precision QRe,QRmu
      common/QReQRmu/QRe,QRmu
      data QRe,QRmu /2*1.d0/
      common/vps2gtree/         !already in invariants and a common there vp13,vp24,
     .     vp13k1,vp24k1,
     .     vp13k2,vp24k2,
     .     vp13k1k2,vp24k1k2,
     .     vp13k1d,vp24k1d,
     .     vp13k2d,vp24k2d,
     .     vp13k1k2d,vp24k1k2d

      double complex vpolc

      if (ifirst.eq.0) then

c         QRe  = 1d0
c         QRmu = 1d0

         ec2 = 4.d0*pi*alpha
         ec  = sqrt(ec2)
         ec4 = ec2*ec2
         mm = ammu
         me = ame
         mm2 = mm*mm
         mm4 = mm2*mm2
         mm6 = mm4*mm2
         me2 = me*me
         me4 = me2*me2
         me6 = me4*me2

         el8 = ec4*ec4
         ifirst = 1

      endif

      call loadinvariants(p1,p2,p3,p4)
      p1k1 = dot(p1,k1)
      p2k1 = dot(p2,k1)
      p3k1 = dot(p3,k1)
      p4k1 = dot(p4,k1)
      p1k2 = dot(p1,k2)
      p2k2 = dot(p2,k2)
      p3k2 = dot(p3,k2)
      p4k2 = dot(p4,k2)
      k1k2 = dot(k1,k2)

      p1k1m1 = 1.d0/p1k1
      p2k1m1 = 1.d0/p2k1
      p3k1m1 = 1.d0/p3k1
      p4k1m1 = 1.d0/p4k1
      p1k2m1 = 1.d0/p1k2
      p2k2m1 = 1.d0/p2k2
      p3k2m1 = 1.d0/p3k2
      p4k2m1 = 1.d0/p4k2

      p2mp4 = p2 - p4
      p3pp4 = p3 + p4

      tm100 = 1.d0/dot(p2mp4,p2mp4)
      tm110 = 1.d0/dot(p2mp4-k1,p2mp4-k1)
      tm102 = 1.d0/dot(p2mp4-k2,p2mp4-k2)
      tm112 = 1.d0/dot(p2mp4-k1-k2,p2mp4-k1-k2)
      sm100 = 1.d0/dot(p3pp4,p3pp4)
      sm110 = 1.d0/dot(p3pp4+k1,p3pp4+k1)
      sm102 = 1.d0/dot(p3pp4+k2,p3pp4+k2)
      sm112 = 1.d0/dot(p3pp4+k1+k2,p3pp4+k1+k2)
      if (iarun.eq.1) then
         tm100 = tm100 * vpolc(dble(1.d0/tm100))
         tm110 = tm110 * vpolc(dble(1.d0/tm110))
         tm102 = tm102 * vpolc(dble(1.d0/tm102))
         tm112 = tm112 * vpolc(dble(1.d0/tm112))
         sm100 = sm100 * vpolc(dble(1.d0/sm100))
         sm110 = sm110 * vpolc(dble(1.d0/sm110))
         sm102 = sm102 * vpolc(dble(1.d0/sm102))
         sm112 = sm112 * vpolc(dble(1.d0/sm112))
      endif
      tm100s = dconjg(tm100)
      tm110s = dconjg(tm110)
      tm102s = dconjg(tm102)
      tm112s = dconjg(tm112)
      sm100s = dconjg(sm100)
      sm110s = dconjg(sm110)
      sm102s = dconjg(sm102)
      sm112s = dconjg(sm112)
      

********
c      print*,'Reenable includes in muemuegg.F!!! stopping'
c      stop
      eeeegg = 0.d0

      include 'eeeegg-frm.f'

      eeeegg = el8 * eeeegg
      return
      end
      
