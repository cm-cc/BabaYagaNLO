      subroutine setwarningLT
#ifdef LOOPTOOLS
#include "looptools.h"
      call setwarndigits(40) 
#endif
      return
      end
**********************************
      subroutine virtualPent_mur_quad(p1,p2,p3,p4,q,virtual,qph,ifc)
      implicit real*16 (z)
#ifdef LOOPTOOLS
#include "looptools.h"
#endif
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision pt1(0:3),pt2(0:3),pt3(0:3),pt4(0:3),qt(0:3)
	double precision qph(40,0:3)
	double precision virtual,virtual1
	real*16 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*16 e2,e6,x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
        real*16 ame2,ammu2,me2,me4,mm2,mm4,amme,mmu
	double precision ame,ammu,convfac,alpha,pi,s,s1
	common/parameters/ame,ammu,convfac,alpha,pi

        external dotq
        real*16 dotq
	double complex elmat,elmat1

	character*2 fs
	character*3 ifc
	common/finalstate/fs

	real*8 deltauv,lambda2,mu2dim,m2min
	common/loop/deltauv,lambda2,mu2dim,m2min

	double complex masses(0:4)
	double complex mominv(1:10)
	double complex momvec(0:3,4)

	double precision routing(5), average

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


	double complex e1
	parameter (e1 = (1.d0,0.d0))

        Complex*32 D0aeppx15x24x34e0p0
        Complex*32 E0apx34x150x13peex24x25p00ee
        Complex*32 D1aeppx15x24x34e0p0
        Complex*32 D1aex15ppx34x240e0p
        Complex*32 D1apx15epx24x34p0e0
        Complex*32 E1apx34x150x13peex24x25p00ee
        Complex*32 E2apx34x150x13peex24x25p00ee
        Complex*32 E3apx34x150x13peex24x25p00ee
        Complex*32 E4apx34x150x13peex24x25p00ee
        Complex*32 E00apx34x150x13peex24x25p00ee
        Complex*32 E11apx34x150x13peex24x25p00ee
        Complex*32 E12apx34x150x13peex24x25p00ee
        Complex*32 E13apx34x150x13peex24x25p00ee
        Complex*32 E14apx34x150x13peex24x25p00ee
        Complex*32 E22apx34x150x13peex24x25p00ee
        Complex*32 E23apx34x150x13peex24x25p00ee
        Complex*32 E24apx34x150x13peex24x25p00ee
        Complex*32 E33apx34x150x13peex24x25p00ee
        Complex*32 E34apx34x150x13peex24x25p00ee
        Complex*32 E44apx34x150x13peex24x25p00ee
        Complex*32 E001apx34x150x13peex24x25p00ee
        Complex*32 E002apx34x150x13peex24x25p00ee
        Complex*32 E003apx34x150x13peex24x25p00ee
        Complex*32 E004apx34x150x13peex24x25p00ee
        Complex*32 E111apx34x150x13peex24x25p00ee
        Complex*32 E112apx34x150x13peex24x25p00ee
        Complex*32 E113apx34x150x13peex24x25p00ee
        Complex*32 E114apx34x150x13peex24x25p00ee
        Complex*32 E122apx34x150x13peex24x25p00ee
        Complex*32 E123apx34x150x13peex24x25p00ee
        Complex*32 E124apx34x150x13peex24x25p00ee
        Complex*32 E133apx34x150x13peex24x25p00ee
        Complex*32 E134apx34x150x13peex24x25p00ee
        Complex*32 E144apx34x150x13peex24x25p00ee
        Complex*32 E222apx34x150x13peex24x25p00ee
        Complex*32 E223apx34x150x13peex24x25p00ee
        Complex*32 E224apx34x150x13peex24x25p00ee
        Complex*32 E233apx34x150x13peex24x25p00ee
        Complex*32 E234apx34x150x13peex24x25p00ee
        Complex*32 E244apx34x150x13peex24x25p00ee
        Complex*32 E333apx34x150x13peex24x25p00ee
        Complex*32 E334apx34x150x13peex24x25p00ee
        Complex*32 E344apx34x150x13peex24x25p00ee
        Complex*32 E444apx34x150x13peex24x25p00ee

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
	p1p2 = dotq(pt1,pt2)
	p1p3 = dotq(pt1,pt3)
	p1p4 = dotq(pt1,pt4)
	kp1  = dotq(pt1,qt)
	p2p3 = dotq(pt2,pt3)
	p2p4 = dotq(pt2,pt4)
	kp2  = dotq(pt2,qt)
	kp3 = dotq(pt3,qt)
	p3p4 = dotq(pt3,pt4)
	kp4  = dotq(pt4,qt)

        me2=ame2
        me4=me2*me2
        mm2=ammu2
        mm4=mm2*mm2

        s=dot(pt1+pt2,pt1+pt2)
        s1=dot(pt3+pt4,pt3+pt4)

        x12  = 2d0*ame2+2.d0*p1p2
        x13  = ame2+ammu2-2.d0*p1p3
        x14  = ame2+ammu2-2.d0*p1p4
        x15  = ame2-2.q0*kp1
        x23  = ame2+ammu2-2d0*p2p3
        x24  = ame2+ammu2-2d0*p2p4
        x25 = ame2-2d0*kp2
        x34 = 2d0*ammu2+2d0*p3p4
        x35 = ammu2+2d0*kp3
        x45 = ammu2+2d0*kp4

        if(ilt.eq.0) then
           ilt=1
           zdeltauv=deltauv
           zphmass2=phmass*phmass
#ifdef LOOPTOOLS
           call ltini
           call setdelta(zdeltauv)
           call setmudim(1.q2)
           call setlambda(zphmass2)
           call setwarningLT
#endif
        endif

#ifdef LOOPTOOLS
       D0aeppx15x24x34e0p0 = D0i(dd0, ame2, ammu2, ammu2, x15, x24, x34,
     > ame2, 0.Q0, ammu2, 0.Q0)

      E0apx34x150x13peex24x25p00ee = E0i(ee0, ammu2, x34, x15, 0.Q0, x13
     > ,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      D1aeppx15x24x34e0p0 = D0i(dd1, ame2, ammu2, ammu2, x15, x24, x34,
     > ame2, 0.Q0, ammu2, 0.Q0)

      D1aex15ppx34x240e0p = D0i(dd1, ame2, x15, ammu2, ammu2, x34, x24,
     > 0.Q0, ame2, 0.Q0, ammu2)

      D1apx15epx24x34p0e0 = D0i(dd1, ammu2, x15, ame2, ammu2, x24, x34,
     > ammu2, 0.Q0, ame2, 0.Q0)

      E1apx34x150x13peex24x25p00ee = E0i(ee1, ammu2, x34, x15, 0.Q0, x13
     > ,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E2apx34x150x13peex24x25p00ee = E0i(ee2, ammu2, x34, x15, 0.Q0, x13
     > ,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E3apx34x150x13peex24x25p00ee = E0i(ee3, ammu2, x34, x15, 0.Q0, x13
     > ,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E4apx34x150x13peex24x25p00ee = E0i(ee4, ammu2, x34, x15, 0.Q0, x13
     > ,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E00apx34x150x13peex24x25p00ee = E0i(ee00, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E11apx34x150x13peex24x25p00ee = E0i(ee11, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E12apx34x150x13peex24x25p00ee = E0i(ee12, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E13apx34x150x13peex24x25p00ee = E0i(ee13, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E14apx34x150x13peex24x25p00ee = E0i(ee14, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E22apx34x150x13peex24x25p00ee = E0i(ee22, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E23apx34x150x13peex24x25p00ee = E0i(ee23, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E24apx34x150x13peex24x25p00ee = E0i(ee24, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E33apx34x150x13peex24x25p00ee = E0i(ee33, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E34apx34x150x13peex24x25p00ee = E0i(ee34, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E44apx34x150x13peex24x25p00ee = E0i(ee44, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E001apx34x150x13peex24x25p00ee = E0i(ee001, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E002apx34x150x13peex24x25p00ee = E0i(ee002, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E003apx34x150x13peex24x25p00ee = E0i(ee003, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E004apx34x150x13peex24x25p00ee = E0i(ee004, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E111apx34x150x13peex24x25p00ee = E0i(ee111, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E112apx34x150x13peex24x25p00ee = E0i(ee112, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E113apx34x150x13peex24x25p00ee = E0i(ee113, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E114apx34x150x13peex24x25p00ee = E0i(ee114, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E122apx34x150x13peex24x25p00ee = E0i(ee122, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E123apx34x150x13peex24x25p00ee = E0i(ee123, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E124apx34x150x13peex24x25p00ee = E0i(ee124, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E133apx34x150x13peex24x25p00ee = E0i(ee133, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E134apx34x150x13peex24x25p00ee = E0i(ee134, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E144apx34x150x13peex24x25p00ee = E0i(ee144, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E222apx34x150x13peex24x25p00ee = E0i(ee222, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E223apx34x150x13peex24x25p00ee = E0i(ee223, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E224apx34x150x13peex24x25p00ee = E0i(ee224, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E233apx34x150x13peex24x25p00ee = E0i(ee233, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E234apx34x150x13peex24x25p00ee = E0i(ee234, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E244apx34x150x13peex24x25p00ee = E0i(ee244, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E333apx34x150x13peex24x25p00ee = E0i(ee333, ammu2, x34, x15, 0.Q0,
     > x13, ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E334apx34x150x13peex24x25p00ee = E0i(ee334, ammu2, x34, x15, 0.Q0,
     > x13, ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E344apx34x150x13peex24x25p00ee = E0i(ee344, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)

      E444apx34x150x13peex24x25p00ee = E0i(ee444, ammu2, x34, x15, 0.Q0,
     > x13,ammu2, ame2, ame2, x24, x25, ammu2, 0.Q0, 0.Q0, ame2, ame2)
#endif

	elmat=0.d0
	include 'mumurad/Virtual/pent_LT/isrfsr1.F'
	if(ifc.eq.'is1'.or.ifc.eq.'fs1') then

	elmat=elmat-elmat1*vpolc(s1)*CONJG(vpolc(s1))
	else
	elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s1))
        endif
	include 'mumurad/Virtual/pent_LT/isrisr1.F'
	elmat=elmat+elmat1*vpolc(s1)*CONJG(vpolc(s))

	virtual=4d0*elmat*e6*alpha/(2d0*pi)

#ifdef LOOPTOOLS
        call clearcache
#endif

       end


      subroutine virtualPent_pr_quad(p1,p2,p3,p4,q,virtual,qph,i)
      implicit none
#ifdef LOOPTOOLS
#include "looptools.h"
#endif
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision pt1(0:3),pt2(0:3),pt3(0:3),pt4(0:3),qt(0:3)
	double precision qph(40,0:3)
	double precision virtual1
        double complex virtual
        complex*32 m1,m2
	real*16 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
	real*16 e2,e6,x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
        real*16 ame2,ammu2,me2,me4,mp2,mp4,ampi2,mp6,me6
        real*16 zdeltauv,zphmass2
        common/hparam/ampi
	double precision ame,ammu,convfac,alpha,pi,s,s1,ampi
	common/parameters/ame,ammu,convfac,alpha,pi


        external dotq
        real*16 dotq
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
	double precision routing(5), average


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
	double complex ffpiis,Cffpiis
	double precision ffpiis2
	double complex ffpifs,Cffpifs
	double precision ffpifs2

************ISR FUNCTIONS
	Complex*32 C0a0x25eee0
	Complex*32 C0appx340p0
	Complex*32 C0apx13x250pe
	Complex*32 C0apx24e0pe
	Complex*32 D0a0x13pex24x25eep0
	Complex*32 D0aeppx25x13x34e0p0
	Complex*32 D0aex34e0x25x15e00e
	Complex*32 D0appx15ex34x240p0e
	Complex*32 E0aex34ex24x13x25x15p0pe00ep
	Complex*32 D1a0x13pex24x25eep0
	Complex*32 D1a0x24px25x13eeep0
	Complex*32 D1aeppx25x13x34e0p0
	Complex*32 D1aex13px34px250ep0
	Complex*32 D1aex34e0x25x15e00e
	Complex*32 D1apx15epx24x34p0e0
	Complex*32 D1apx250x24x13ep0ee
	Complex*32 D1apx34ex13px25p00e
	Complex*32 D1apx34x15x24pep00e
	Complex*32 D1ax15x24px34pe0ep0
	Complex*32 E1aex34ex24x13x25x15p0pe00ep
	Complex*32 E1ax13ex34ex24px25x15p0pe00e
	Complex*32 D2aex34e0x25x15e00e
	Complex*32 E2aex34ex24x13x25x15p0pe00ep
	Complex*32 D3aex34e0x25x15e00e
	Complex*32 E3aex34ex24x13x25x15p0pe00ep
	Complex*32 D00aex34e0x25x15e00e
	Complex*32 D11aex34e0x25x15e00e
	Complex*32 E11aex34ex24x13x25x15p0pe00ep
	Complex*32 D12aex34e0x25x15e00e
	Complex*32 E12ax13ex34ex24px25x15p0pe00e
	Complex*32 D13aex34e0x25x15e00e
	Complex*32 E13aex34ex24x13x25x15p0pe00ep
	Complex*32 E13ax13ex34ex24px25x15p0pe00e
	Complex*32 E14ax13ex34ex24px25x15p0pe00e
	Complex*32 E22aex34ex24x13x25x15p0pe00ep
	Complex*32 D23aex34e0x25x15e00e
	Complex*32 E23aex34ex24x13x25x15p0pe00ep
	Complex*32 D22aex34e0x25x15e00e
	Complex*32 E00aex34ex24x13x25x15p0pe00ep
	Complex*32 E11ax13ex34ex24px25x15p0pe00e
	Complex*32 E33aex34ex24x13x25x15p0pe00ep
	Complex*32 E12aex34ex24x13x25x15p0pe00ep
        Complex*32 C0appx34pm1p
********************

*************GVMD COMMON
        common/gvmdrad/indexnores
        integer indexnores
        
        virtual=(0.d0,0.d0)
	e2 = 4.d0*pi*alpha
	e6 = e2*e2*e2
	ame2=ame*ame*1.q0
	ampi2=ampi*ampi*1.q0
        me2=ame2
        mp2=ampi2
        me4=me2**2
        me6=me2**3
        mp4=mp2**2
        mp6=mp2**3

        if(i.eq.1) then
         pt3=p3
         pt4=p4
        elseif(i.eq.2) then
         pt3=p4
         pt4=p3
        endif

	p1p2 = dotq(p1,p2)
	p1p3 = dotq(p1,pt3)
	p1p4 = dotq(p1,pt4)
	kp1  = dotq(p1,q)
	p2p3 = dotq(p2,pt3)
	p2p4 = dotq(p2,pt4)
	kp2  = dotq(p2,q)
	kp3 = dotq(pt3,q)
	p3p4 = dotq(pt3,pt4)
	kp4  = dotq(pt4,q)

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
	s=2.d0*(ame**2+p1p2)
	s1=2.d0*(ampi**2+p3p4)
	if(iffpi.gt.0) then
		call ffpi(s,ffpifs,ffpifs2)
		call ffpi(s1,ffpiis,ffpiis2)
	endif	

        Cffpiis=CONJG(ffpiis)
        Cffpifs=CONJG(ffpifs)

        if(indexnores.eq.1) then
           ffpiis=(1.d0,0.d0)
           ffpifs=(1.d0,0.d0)
	endif

******ALPHA RUNNING
	ffpiis=ffpiis*vpolc(s1)
        Cffpiis=Cffpiis*CONJG(vpolc(s1))

	ffpifs=ffpifs*vpolc(s)
        Cffpifs=Cffpifs*CONJG(vpolc(s))
        
        ffpiis2=ffpiis*Cffpiis
        ffpifs2=ffpifs*Cffpifs
****************



        if(ilt.eq.0) then
           ilt=1
           zdeltauv=deltauv
           zphmass2=phmass*phmass
#ifdef LOOPTOOLS
	 call ltini
	 call setdelta(zdeltauv)
         call setmudim(1.q2)
         call setlambda(zphmass2)
         call setwarningLT
#endif
        endif

#ifdef LOOPTOOLS

      C0a0x25eee0 = C0i(cc0, 0.Q0, x25, ame2, ame2,
     > ame2, 0.Q0)

      C0appx340p0 = C0i(cc0, mp2, mp2, x34, 0.Q0,
     > mp2, 0.Q0)

      C0apx13x250pe = C0i(cc0, mp2, x13, x25, 0.Q0,
     > mp2, ame2)

      C0apx24e0pe = C0i(cc0, mp2, x24, ame2, 0.Q0,
     > mp2, ame2)

      D0a0x13pex24x25eep0 = D0i(dd0, 0.Q0, x13, mp2, ame2,
     > x24, x25, ame2, ame2, mp2, 0.Q0)

      D0aeppx25x13x34e0p0 = D0i(dd0, ame2, mp2, mp2, x25,
     > x13, x34, ame2, 0.Q0, mp2, 0.Q0)

      D0aex34e0x25x15e00e = D0i(dd0, ame2, x34, ame2, 0.Q0,
     > x25, x15, ame2, 0.Q0, 0.Q0, ame2)

      D0appx15ex34x240p0e = D0i(dd0, mp2, mp2, x15, ame2,
     > x34, x24, 0.Q0, mp2, 0.Q0, ame2)

      E0aex34ex24x13x25x15p0pe00ep = E0i(ee0, ame2, x34, ame2,
     > x24, x13, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2, mp2)

      D1a0x13pex24x25eep0 = D0i(dd1, 0.Q0, x13, mp2, ame2,
     > x24, x25, ame2, ame2, mp2, 0.Q0)

      D1a0x24px25x13eeep0 = D0i(dd1, 0.Q0, x24, mp2, x25,
     > x13, ame2, ame2, ame2, mp2, 0.Q0)

      D1aeppx25x13x34e0p0 = D0i(dd1, ame2, mp2, mp2, x25,
     > x13, x34, ame2, 0.Q0, mp2, 0.Q0)

      D1aex13px34px250ep0 = D0i(dd1, ame2, x13, mp2, x34,
     > mp2, x25, 0.Q0, ame2, mp2, 0.Q0)

      D1aex34e0x25x15e00e = D0i(dd1, ame2, x34, ame2, 0.Q0,
     > x25, x15, ame2, 0.Q0, 0.Q0, ame2)

      D1apx15epx24x34p0e0 = D0i(dd1, mp2, x15, ame2, mp2,
     > x24, x34, mp2, 0.Q0, ame2, 0.Q0)

      D1apx250x24x13ep0ee = D0i(dd1, mp2, x25, 0.Q0, x24,
     > x13, ame2, mp2, 0.Q0, ame2, ame2)

      D1apx34ex13px25p00e = D0i(dd1, mp2, x34, ame2, x13,
     > mp2, x25, mp2, 0.Q0, 0.Q0, ame2)

      D1apx34x15x24pep00e = D0i(dd1, mp2, x34, x15, x24,
     > mp2, ame2, mp2, 0.Q0, 0.Q0, ame2)

      D1ax15x24px34pe0ep0 = D0i(dd1, x15, x24, mp2, x34,
     > mp2, ame2, 0.Q0, ame2, mp2, 0.Q0)

      E1aex34ex24x13x25x15p0pe00ep = E0i(ee1, ame2, x34, ame2,
     > x24, x13, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2, mp2)

      E1ax13ex34ex24px25x15p0pe00e = E0i(ee1, x13, ame2, x34,
     > ame2, x24, mp2, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2)

      D2aex34e0x25x15e00e = D0i(dd2, ame2, x34, ame2, 0.Q0,
     > x25, x15, ame2, 0.Q0, 0.Q0, ame2)

      E2aex34ex24x13x25x15p0pe00ep = E0i(ee2, ame2, x34, ame2,
     > x24, x13, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2, mp2)

      D3aex34e0x25x15e00e = D0i(dd3, ame2, x34, ame2, 0.Q0,
     > x25, x15, ame2, 0.Q0, 0.Q0, ame2)

      E3aex34ex24x13x25x15p0pe00ep = E0i(ee3, ame2, x34, ame2,
     > x24, x13, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2, mp2)

      D00aex34e0x25x15e00e = D0i(dd00, ame2, x34, ame2, 0.Q0,
     > x25, x15, ame2, 0.Q0, 0.Q0, ame2)

      D11aex34e0x25x15e00e = D0i(dd11, ame2, x34, ame2, 0.Q0,
     > x25, x15, ame2, 0.Q0, 0.Q0, ame2)

      E11aex34ex24x13x25x15p0pe00ep = E0i(ee11, ame2, x34, ame2,
     > x24, x13, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2, mp2)

      D12aex34e0x25x15e00e = D0i(dd12, ame2, x34, ame2, 0.Q0,
     > x25, x15, ame2, 0.Q0, 0.Q0, ame2)

      E12ax13ex34ex24px25x15p0pe00e = E0i(ee12, x13, ame2, x34,
     > ame2, x24, mp2, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2)

      D13aex34e0x25x15e00e = D0i(dd13, ame2, x34, ame2, 0.Q0,
     > x25, x15, ame2, 0.Q0, 0.Q0, ame2)

      E13aex34ex24x13x25x15p0pe00ep = E0i(ee13, ame2, x34, ame2,
     > x24, x13, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2, mp2)

      E13ax13ex34ex24px25x15p0pe00e = E0i(ee13, x13, ame2, x34,
     > ame2, x24, mp2, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2)

      E14ax13ex34ex24px25x15p0pe00e = E0i(ee14, x13, ame2, x34,
     > ame2, x24, mp2, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2)

      E22aex34ex24x13x25x15p0pe00ep = E0i(ee22, ame2, x34, ame2,
     > x24, x13, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2, mp2)

      D23aex34e0x25x15e00e = D0i(dd23, ame2, x34, ame2, 0.Q0,
     > x25, x15, ame2, 0.Q0, 0.Q0, ame2)

      E23aex34ex24x13x25x15p0pe00ep = E0i(ee23, ame2, x34, ame2,
     > x24, x13, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2, mp2)

      D22aex34e0x25x15e00e = D0i(dd22, ame2, x34, ame2, 0.Q0,
     > x25, x15, ame2, 0.Q0, 0.Q0, ame2)

      E00aex34ex24x13x25x15p0pe00ep = E0i(ee00, ame2, x34, ame2,
     > x24, x13, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2, mp2)

      E11ax13ex34ex24px25x15p0pe00e = E0i(ee11, x13, ame2, x34,
     > ame2, x24, mp2, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2)

      E33aex34ex24x13x25x15p0pe00ep = E0i(ee33, ame2, x34, ame2,
     > x24, x13, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2, mp2)

      E12aex34ex24x13x25x15p0pe00ep = E0i(ee12, ame2, x34, ame2,
     > x24, x13, x25, x15, mp2, 0.Q0, mp2, ame2,
     > 0.Q0, 0.Q0, ame2, mp2)


#endif



      include 'pipirad/Virtual_new/pent3_complete_isr.f'
      include 'pipirad/Virtual_new/pent3_complete_fsr.f'
     	elmati=elmat1*ffpiis*Cffpiis
 	elmatf=elmat2*Cffpifs*ffpiis
        if(i.eq.2) then
           elmati=-elmati
        endif
!        print*,'isr pent',elmati,elmatf
 	virtual=virtual+elmati+elmatf
!        print*,'full interno',virtual
        virtual=4d0*virtual*e6*alpha/(2d0*pi)
!        virtual= C0appx34pm1p

#ifdef LOOPTOOLS
        call clearcache
#endif
       end


      subroutine mixedfsr_quad(p1,p2,p3,p4,q,mq,virtual,qph)
      implicit complex*32 (z)
#ifdef LOOPTOOLS
#include "looptools.h"
#endif
	double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
	double precision pt1(0:3),pt2(0:3),pt3(0:3),pt4(0:3),qt(0:3)
	double precision qph(40,0:3)
	double precision virtual1

        real*16 p1q(0:3),p2q(0:3),p3q(0:3),p4q(0:3),qq(0:3)
        complex*32 m1,re,mq,dzphi,m12
        complex*32 virtual,virtualf,deltafsr
	complex*32 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
 	complex*32 e2,e6,x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
        complex*32 ame2,ammu2,me2,me4,mp2,mp4,ampi2,mp6,me6
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
	double precision routing(5), average
	

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
        Complex*32 A0ap,A0am1,B0appm1
        Complex*32 B0ax34pp,C0appx34pm1p
        Complex*32 DB0apm1p
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
        if(ilt.eq.0) then
           ilt=1
           zdeltauv=deltauv*1.q0
           zphmass2=(phmass*1.q0)*(phmass*1.q0)
           zmu2dim=mu2dim*1.q0
#ifdef LOOPTOOLS
	 call ltini
	 call setdelta(zdeltauv)                                                                                                                                                       
         call setmudim(zmu2dim)
         call setlambda(zphmass2) 
         call setwarningLT
#endif   
        endif

        m1       = sqrt(mq)
        m12      = mq
        virtual  = (0.q0,0.q0)
#ifdef LOOPTOOLS

        A0ap     = A0C(mp2)

        A0am1    = A0C(mq)
        
        B0ax34pp = B0C(x34,mp2,mp2)

        B0appm1  = B0C(mp2,mp2,mq)
        
        DB0apm1p = DB0C(mp2,mq,mp2)

        C0appx34pm1p = C0C(mp2,mp2,x34,mp2,mq,mp2)

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
#ifdef LOOPTOOLS
        call clearcache
#endif

        return
        end


      subroutine purefsr_quad(p1,p2,p3,p4,q,mq,virtual,qph)
      implicit complex*32 (z)
!      implicit none
#ifdef LOOPTOOLS
#include "looptools.h"
#endif
        double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
        double precision pt1(0:3),pt2(0:3),pt3(0:3),pt4(0:3),qt(0:3)
        double precision qph(40,0:3)
	double precision virtual1

        complex*32 m1,re,mq,dzphi,dmpi,m12
        complex*32 virtual,deltafsr
	complex*32 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
 	complex*32 e2,e6,x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
        complex*32 ame2,ammu2,me2,me4,mp2,mp4,ampi2,mp6,me6
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
	double precision routing(5), average
	

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
        complex*32 Cim(2)
        complex*32 Cijm(4)
        complex*32 C0m

        complex*32 Dim(3)
        complex*32 Dijm(7)
        complex*32 D0m

        
        Complex*32 A0ap,A0am1,B0appm1
        Complex*32 B0ax34pp,C0appx34pm1p
        Complex*32 B0ax35pm1,B0ax45pm1
        Complex*32 DB0apm1p

        Complex*32 C0apx350pm1p
        Complex*32 C1apx350pm1p,C2apx350pm1p,C00apx350pm1p
        Complex*32 C11apx350pm1p,C12apx350pm1p,C22apx350pm1p
        Complex*32 C00apx12x35m1pp,C12apx12x35m1pp,C0apx12x35m1pp
        
        Complex*32 C0apx450pm1p
        Complex*32 C1apx450pm1p,C2apx450pm1p,C00apx450pm1p
        Complex*32 C11apx450pm1p,C12apx450pm1p,C22apx450pm1p

       
        Complex*32 C0apx12x45m1pp,C1apx12x45m1pp,C1ax45x12pm1pp
        Complex*32 C00apx12x45m1pp,C11apx12x45m1pp,C11ax45x12pm1pp
        Complex*32 C12apx12x45m1pp

        Complex*32 C0ax45x12pm1pp,C00ax45x12pm1pp,C12ax45x12pm1pp

        
        Complex*32 C0apx45x12pm1p,C2apx45x12pm1p,C2ax45px12pm1p
        Complex*32 C00apx45x12pm1p,C22apx45x12pm1p
        Complex*32 C22ax45px12pm1p

        
        Complex*32 C0ax35x12pm1pp,C1apx12x35m1pp,C1ax35x12pm1pp
        Complex*32 C00ax35x12pm1pp,C11apx12x35m1pp,C11ax35x12pm1pp
        Complex*32 C12ax35x12pm1pp

        Complex*32 C0apx35x12pm1p,C2apx35x12pm1p,C2ax35px12pm1p
        Complex*32 C00apx35x12pm1p,C22apx35x12pm1p
        Complex*32 C22ax35px12pm1p

        
!!!!!!!!!!BOX

        Complex*32 D0appx120x34x35pm1pp,D2app0x12x34x35pm1pp
        Complex*32 D3app0x12x34x35pm1pp,D00app0x12x34x35pm1pp
        Complex*32 D33app0x12x34x35pm1pp,D0app0x12x34x35pm1pp        
        Complex*32 D2appx120x34x35pm1pp,D3appx120x34x35pm1pp
        Complex*32 D00appx120x34x35pm1pp,D13a0ppx12x35x34ppm1p
        Complex*32 D22app0x12x34x35pm1pp,D22appx120x34x35pm1pp
        Complex*32 D23app0x12x34x35pm1pp,D23appx120x34x35pm1pp
        Complex*32 D33appx120x34x35pm1pp

        Complex*32 C1ax12x340ppp,C00a0x12x34ppp
        Complex*32 C11ax12x340ppp,C12a0x34x12ppp
        Complex*32 C12a0x12x34ppp
        
!!!!!!!!!!!!!!!!!!!!!!
        Complex*32 D0ap0x12px35x34m1ppp
        Complex*32 D1ap0x12px35x34m1ppp,D2ap0x12px35x34m1ppp
        Complex*32 D3ap0x12px35x34m1ppp
        Complex*32 D00ap0x12px35x34m1ppp,D11ap0x12px35x34m1ppp
        Complex*32 D12ap0x12px35x34m1ppp
        Complex*32 D13ap0x12px35x34m1ppp,D22ap0x12px35x34m1ppp
        Complex*32 D23ap0x12px35x34m1ppp,D33ap0x12px35x34m1ppp

        Complex*32 D0ap0x12px45x34m1ppp
        Complex*32 D1ap0x12px45x34m1ppp,D2ap0x12px45x34m1ppp
        Complex*32 D3ap0x12px45x34m1ppp
        Complex*32 D00ap0x12px45x34m1ppp,D11ap0x12px45x34m1ppp
        Complex*32 D12ap0x12px45x34m1ppp
        Complex*32 D13ap0x12px45x34m1ppp,D22ap0x12px45x34m1ppp
        Complex*32 D23ap0x12px45x34m1ppp,D33ap0x12px45x34m1ppp

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
        if(ilt.eq.0) then
           ilt=1
           zdeltauv=deltauv*1.q0
           zphmass2=(phmass*1.q0)*(phmass*1.q0)
           zmu2dim=mu2dim*1.q0

#ifdef LOOPTOOLS
	 call ltini
	 call setdelta(zdeltauv)                                                                                                                                                       
         call setmudim(zmu2dim)
         call setlambda(zphmass2) 
         call setwarningLT
#endif   
        endif

        m1      = sqrt(mq)
        m12     = mq
        virtual =(0.q0,0.q0)

	! Inizializzazione a zero (parte A & B)
	A0ap        = (0.Q0, 0.Q0)
	A0am1       = (0.Q0, 0.Q0)
	B0appm1     = (0.Q0, 0.Q0)
        DB0apm1p    = (0.Q0, 0.Q0)
	B0ax34pp    = (0.Q0, 0.Q0)
	C0appx34pm1p= (0.Q0, 0.Q0)
	B0ax35pm1   = (0.Q0, 0.Q0)
	B0ax45pm1   = (0.Q0, 0.Q0)

	! C 350
	C0apx350pm1p  = (0.Q0, 0.Q0)
	C1apx350pm1p  = (0.Q0, 0.Q0)
	C2apx350pm1p  = (0.Q0, 0.Q0)
	C00apx350pm1p = (0.Q0, 0.Q0)
        C11apx350pm1p = (0.Q0, 0.Q0)
	C12apx350pm1p = (0.Q0, 0.Q0)
	C22apx350pm1p = (0.Q0, 0.Q0)

	! C 450
	C0apx450pm1p  = (0.Q0, 0.Q0)
	C1apx450pm1p  = (0.Q0, 0.Q0)
	C2apx450pm1p  = (0.Q0, 0.Q0)
	C00apx450pm1p = (0.Q0, 0.Q0)
	C11apx450pm1p = (0.Q0, 0.Q0)
	C12apx450pm1p = (0.Q0, 0.Q0)
	C22apx450pm1p = (0.Q0, 0.Q0)

	! C 12x45
	C0apx12x45m1pp  = (0.Q0, 0.Q0)
	C1apx12x45m1pp  = (0.Q0, 0.Q0)
	C1ax45x12pm1pp  = (0.Q0, 0.Q0)
	C00apx12x45m1pp = (0.Q0, 0.Q0)
	C11apx12x45m1pp = (0.Q0, 0.Q0)
	C11ax45x12pm1pp = (0.Q0, 0.Q0)
	C12apx12x45m1pp = (0.Q0, 0.Q0)
	C0ax45x12pm1pp  = (0.Q0, 0.Q0)
	C00ax45x12pm1pp = (0.Q0, 0.Q0)
	C12ax45x12pm1pp = (0.Q0, 0.Q0)

	! C 45x12 p
	C0apx45x12pm1p  = (0.Q0, 0.Q0)
	C2apx45x12pm1p  = (0.Q0, 0.Q0)
	C2ax45px12pm1p  = (0.Q0, 0.Q0)
	C00apx45x12pm1p = (0.Q0, 0.Q0)
	C22apx45x12pm1p = (0.Q0, 0.Q0)
	C22ax45px12pm1p = (0.Q0, 0.Q0)

	! C 35x12
	C0ax35x12pm1pp  = (0.Q0, 0.Q0)
        C0apx12x35m1pp  = (0.Q0, 0.Q0)
	C1apx12x35m1pp  = (0.Q0, 0.Q0)
	C1ax35x12pm1pp  = (0.Q0, 0.Q0)
        C12apx12x35m1pp = (0.Q0, 0.Q0)
	C00ax35x12pm1pp = (0.Q0, 0.Q0)
	C11apx12x35m1pp = (0.Q0, 0.Q0)
        C00apx12x35m1pp = (0.Q0, 0.Q0)
        C11ax35x12pm1pp = (0.Q0, 0.Q0)
	C12ax35x12pm1pp = (0.Q0, 0.Q0)
	C0apx35x12pm1p  = (0.Q0, 0.Q0)
	C2apx35x12pm1p  = (0.Q0, 0.Q0)
	C2ax35px12pm1p  = (0.Q0, 0.Q0)
	C00apx35x12pm1p = (0.Q0, 0.Q0)
	C22apx35x12pm1p = (0.Q0, 0.Q0)
	C22ax35px12pm1p = (0.Q0, 0.Q0)

	! D 35
	D0appx120x34x35pm1pp   = (0.Q0, 0.Q0)
	D2app0x12x34x35pm1pp   = (0.Q0, 0.Q0)
	D3app0x12x34x35pm1pp   = (0.Q0, 0.Q0)
	D00app0x12x34x35pm1pp  = (0.Q0, 0.Q0)
	D33app0x12x34x35pm1pp  = (0.Q0, 0.Q0)
	D0app0x12x34x35pm1pp   = (0.Q0, 0.Q0)
	D2appx120x34x35pm1pp   = (0.Q0, 0.Q0)
	D3appx120x34x35pm1pp   = (0.Q0, 0.Q0)
	D00appx120x34x35pm1pp  = (0.Q0, 0.Q0)
	D13a0ppx12x35x34ppm1p  = (0.Q0, 0.Q0)
	D22app0x12x34x35pm1pp  = (0.Q0, 0.Q0)
	D22appx120x34x35pm1pp  = (0.Q0, 0.Q0)
	D23app0x12x34x35pm1pp  = (0.Q0, 0.Q0)
	D23appx120x34x35pm1pp  = (0.Q0, 0.Q0)
	D33appx120x34x35pm1pp  = (0.Q0, 0.Q0)

	! C extra
	C1ax12x340ppp   = (0.Q0, 0.Q0)
	C00a0x12x34ppp  = (0.Q0, 0.Q0)
	C11ax12x340ppp  = (0.Q0, 0.Q0)
	C12a0x34x12ppp  = (0.Q0, 0.Q0)
	C12a0x12x34ppp  = (0.Q0, 0.Q0)


#ifdef LOOPTOOLS

        A0ap     = A0C(mp2)

        A0am1    = A0C(mq)
        
        B0ax34pp = B0C(x34,mp2,mp2)

        B0appm1  = B0C(mp2,mp2,mq)

        DB0apm1p = DB0C(mp2,mq,mp2)
                
        B0ax35pm1 = B0C(x35,mp2,mq)
        
        B0ax45pm1 = B0C(x45,mp2,mq)
        
        C0appx34pm1p = C0C(mp2,mp2,x34,mp2,mq,mp2)        

#endif
        
        call Cij(p3,-q,mp2,mq,mp2,C0apx350pm1p,Cim,Cijm)
        C1apx350pm1p  = Cim(1)
        C2apx350pm1p  = Cim(2)
        C00apx350pm1p = Cijm(1)
        C11apx350pm1p =	Cijm(2)
        C12apx350pm1p =	Cijm(3)
        C22apx350pm1p =	Cijm(4)

        call Cij(p4,-q,mp2,mq,mp2,C0apx450pm1p,Cim,Cijm)
	C1apx450pm1p  = Cim(1)
	C2apx450pm1p  = Cim(2)
	C00apx450pm1p =	Cijm(1)
        C11apx450pm1p = Cijm(2)
        C12apx450pm1p = Cijm(3)
        C22apx450pm1p = Cijm(4)

        
        call Cij(-p3,p4+q,mq,mp2,mp2,C0apx12x45m1pp,Cim,Cijm)
        C1apx12x45m1pp  = Cim(1)
        C1ax45x12pm1pp  = Cim(2)
        C00apx12x45m1pp = Cijm(1)
        C11apx12x45m1pp = Cijm(2)
        C12apx12x45m1pp = Cijm(3)
        C11ax45x12pm1pp = Cijm(4)

        C0ax45x12pm1pp  = C0apx12x45m1pp
        C00ax45x12pm1pp = C00apx12x45m1pp
        C12ax45x12pm1pp = C12apx12x45m1pp

        C0apx45x12pm1p  = C0apx12x45m1pp
        C2apx45x12pm1p  = Cim(2)
        C00apx45x12pm1p = Cijm(1)
        C2ax45Px12pm1p  = Cim(1)
        C22apx45x12pm1p = Cijm(4)
        C22ax45px12pm1p = Cijm(2)

        

        call Cij(-p4,q+p3,mq,mp2,mp2,C0ax35x12pm1pp,Cim,Cijm)
        C1apx12x35m1pp  = Cim(1)
        C1ax35x12pm1pp  = Cim(2)
        C00ax35x12pm1pp = Cijm(1)
        C11apx12x35m1pp = Cijm(2)
        C00apx12x35m1pp = Cijm(1)
        C12apx12x35m1pp = Cijm(3)
        C12ax35x12pm1pp = Cijm(3)
        C11ax35x12pm1pp = Cijm(4)
!!!!!!!!!!!For simmetry of PaVe
        C0apx35x12pm1p  = C0ax35x12pm1pp
        c0apx12x35m1pp  = C0ax35x12pm1pp
        C2apx35x12pm1p  = Cim(2)
        C00apx35x12pm1p = Cijm(1)
        C2ax35px12pm1p  = Cim(1)
        C22apx35x12pm1p = Cijm(4)
        C22ax35px12pm1p = Cijm(2)
        C12ax35x12pm1pp = Cijm(3)

        
        call Cij(p1+p2,q,mp2,mp2,mp2,C0m,Cim,Cijm)
        C1ax12x340ppp  = Cim(1)
        C11ax12x340ppp = Cijm(2)
        C12a0x34x12ppp = Cijm(3)

        call Cij(q,-p3-p4,mp2,mp2,mp2,C0m,Cim,Cijm)
        C00a0x12x34ppp = Cijm(1)
        C12a0x12x34ppp = Cijm(3)

!        fscoll=0
        
!        if(fscoll.eq.1) then
           call Dij(-p3,-p3-q,p4,mq,mp2,mp2,mp2,D0m,Dim,Dijm)
           D0ap0x12px35x34m1ppp  = D0m
           D1ap0x12px35x34m1ppp  = Dim(1)
           D2ap0x12px35x34m1ppp  = Dim(2)
           D3ap0x12px35x34m1ppp  = Dim(3)
           D00ap0x12px35x34m1ppp = Dijm(1)
           D11ap0x12px35x34m1ppp = Dijm(2)
           D12ap0x12px35x34m1ppp = Dijm(3)
           D13ap0x12px35x34m1ppp = Dijm(4)
           D22ap0x12px35x34m1ppp = Dijm(5)
           D23ap0x12px35x34m1ppp = Dijm(6)
           D33ap0x12px35x34m1ppp = Dijm(7)
           
           call Dij(-p4,-p4-q,p3,mq,mp2,mp2,mp2,D0m,Dim,Dijm)
           D0ap0x12px45x34m1ppp  = D0m
           D1ap0x12px45x34m1ppp  = Dim(1)
           D2ap0x12px45x34m1ppp  = Dim(2)
           D3ap0x12px45x34m1ppp  = Dim(3)
           D00ap0x12px45x34m1ppp = Dijm(1)
           D11ap0x12px45x34m1ppp = Dijm(2)
           D12ap0x12px45x34m1ppp = Dijm(3)
           D13ap0x12px45x34m1ppp = Dijm(4)
           D22ap0x12px45x34m1ppp = Dijm(5)
           D23ap0x12px45x34m1ppp = Dijm(6)
           D33ap0x12px45x34m1ppp = Dijm(7)
           
        include 'pipirad/gvmd/purefsrtot_isr.f'

        if(ifisrfsr.eq.11) then
           
           virtual = elmat1*ffpifs*CONJG(ffpiis)

        elseif(ifisrfsr.eq.01) then

           virtual = (0.q0,0.q0)

        endif

        
        include 'pipirad/gvmd/purefsrtot_fsr.f'
                
        virtual = virtual + elmat2*ffpifs2
!        print*,'PURE FSR NO CT',virtual

#ifdef LOOPTOOLS
        call clearcache
#endif

        return
        end

      subroutine ct_purefsr_quad(p1,p2,p3,p4,q,mq,virtual,qph,index)
      implicit complex*32 (z)
!      implicit none
#ifdef LOOPTOOLS
#include "looptools.h"
#endif
        double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
        double precision pt1(0:3),pt2(0:3),pt3(0:3),pt4(0:3),qt(0:3)
        double precision qph(40,0:3)
	double precision virtual1

        complex*32 m1,re,mq,dzphi,dmpi,m12,dzphi1,dmpi1
        complex*32 virtual,deltafsr
	complex*32 p1p2,p1p3,p1p4,kp1,p2p3,p2p4,kp2,p3p4,kp3,kp4
 	complex*32 e2,e6,x12,x13,x14,x15,x23,x24,x25,x34,x35,x45
        complex*32 ame2,ammu2,me2,me4,mp2,mp4,ampi2,mp6,me6
        real*16 zdeltauv,zphmass2,zmu2dim
        real*16 e,e3,e5,pi2
        common/hparam/ampi
	double precision ame,ammu,convfac,alpha,pi,s,s1,ampi
	common/parameters/ame,ammu,convfac,alpha,pi

        integer index
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
	double precision routing(5), average
	

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

        common/radoption/ifisrfsr
        integer ifisrfsr

        common/pavered/p1red,p2red,p3red,p4red,kred
        double precision p1red(0:3),p2red(0:3),p3red(0:3),p4red(0:3)
        double precision kred(0:3)

        common/gvmdchecks/iscoll,fscoll
        integer iscoll,fscoll
        common/checkform/intltgvmd
        integer intltgvm


        complex*32 A0ap
        complex*32 A0am1
        complex*32 B0appm1
        complex*32 DB0apm1p 

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
        pi2 = pi*pi
        
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

        pt3=p3
        pt4=p4

	p1p2 = dotq(p1,p2)*re
	p1p3 = dotq(p1,pt3)*re
	p1p4 = dotq(p1,pt4)*re
	kp1  = dotq(p1,q)*re
	p2p3 = dotq(p2,pt3)*re
	p2p4 = dotq(p2,pt4)*re
	kp2  = dotq(p2,q)*re
	kp3 = dotq(pt3,q)*re
	p3p4 = dotq(pt3,pt4)*re
	kp4  = dotq(pt4,q)*re


      	x12 = dotq(p1+p2,p1+p2)*re
	x13 = dotq(p1-p3,p1-p3)*re
        x14 = dotq(p1-p4,p1-p4)*re
        x15 = dotq(p1-q,p1-q)*re
	x23 = dotq(p2-p3,p2-p3)*re
	x24 = dotq(p2-p4,p2-p4)*re
	x25 = dotq(p2-q,p2-q)*re
	x34 = dotq(p3+p4,p3+p4)*re
        x35 = dotq(p3+q,p3+q)*re
        x45 = dotq(p4+q,p4+q)*re

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

	A0ap        = (0.Q0, 0.Q0)
	A0am1       = (0.Q0, 0.Q0)
	B0appm1     = (0.Q0, 0.Q0)
        DB0apm1p    = (0.Q0, 0.Q0)

        
!        print*,'deltauv',deltauv
        if(ilt.eq.0) then
           ilt=1
           zdeltauv=deltauv*1.q0
           zphmass2=(phmass*1.q0)*(phmass*1.q0)
           zmu2dim=mu2dim*1.q0

#ifdef LOOPTOOLS
	 call ltini
	 call setdelta(zdeltauv)                                                                                                                                                       
         call setmudim(zmu2dim)
         call setlambda(zphmass2) 
         call setwarningLT
#endif   
           
        endif
#ifdef LOOPTOOLS

        A0ap     = A0C(mp2)

        A0am1    = A0C(mq)
        
        B0appm1  = B0C(mp2,mp2,mq)

        DB0apm1p = DB0C(mp2,mq,mp2)
                

#endif
        virtual=(0.q0,0.q0)
***************CT_s factorise over born (ISR), so here we have (|M_FSR^0|^2+M_FSR^0 x M_ISR^{0,\dagger})*dzphi****
******************************************************************************************************************
***************REMARK: THE CT_S FOR THE POINT LIKE PART ARE EVALUATED IN CT_PR EXPLOITING THE EXISTING************
******************************************POINTLIKE CALCULATION***************************************************        

        call FSR_pr(p1,p2,p3,p4,q,elmatfsrd,qph)

        call IFI_pr(p1,p2,p3,p4,q,elmatifid,qph)

        if(index.eq.0) then
           dzphi = (2.q0*B0appm1+(4.q0*mp2-mq)*DB0apm1p)

**************TO MATCH THE FACTOR IN LOOPS
           dzphi  = dzphi/(4.q0*e6)
*****************************************
           dmpi=-(2.q0*A0am1-A0ap+(4.q0*mp2-mq)*B0appm1)/mp2
           
        elseif(index.eq.1) then
           dzphi = (2.q0*B0appm1+(4.q0*mp2-mq)*DB0apm1p)

**************TO MATCH THE FACTOR IN LOOPS
           dzphi  = dzphi/(4.q0*e6)
*****************************************
           
           dmpi=-(2.q0*A0am1-A0ap+(4.q0*mp2-mq)*B0appm1)/mp2

           
           mq=mq*(1.q0+1.q-10)

#ifdef LOOPTOOLS
           
        A0ap     = A0C(mp2)

        A0am1    = A0C(mq)
        
        B0appm1  = B0C(mp2,mp2,mq)

        DB0apm1p = DB0C(mp2,mq,mp2)


#endif
        dzphi1=(2.q0*B0appm1+(4.q0*mp2-mq)*DB0apm1p)/(4.q0*e6)
        dmpi1 = -(2.q0*A0am1-A0ap+(4.q0*mp2-mq)*B0appm1)/mp2

        dmpi  = (dmpi1-dmpi)/(1.q-10*mq)
        dzphi = (dzphi1-dzphi)/(1.q-10*mq)

        mq=mq/(1.q0+1.q-10)
        endif
         
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
        call clearcache
        
        return
        end

*----------------------------------------------*
*     IN THIS ROUTINE WE CHECK IF AN INVARIANT *
*     IS CLOSE TO MASS-SHELL AND WE FORCE IT   *
*     TO BE == M_I^2 (FUNDAMENTAL IF WE WANT   *
*     PAVE TO HAVE THE RIGHT IR BEHAVIOUR      *
*----------------------------------------------*      

      function closetomassc(x)
      implicit none
      complex*32 x,closetomassc
      double precision ame,ammu,convfac,alpha,pi,s,s1,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      integer i,j
!!!!!!!!!!!checckare se l'invariante x è vicina alle masse fisiche, se si vengono messi forzatamente onshel	
      complex*32 ameq, ampiq,ammuq
!!!!!!!!!!!!!!!COMMON TO CHECK IN XIJ REDUCTION
      common/pavered/p1red,p2red,p3red,p4red,kred
      external dotq
      real*16 dotq
      double precision p1red(0:3),p2red(0:3),p3red(0:3),p4red(0:3)
      double precision kred(0:3)
      real*16 mom(5,0:3)
      Complex*32 masses(3),aux
      real*16 p1q(0:3),p2q(0:3),p3q(0:3),p4q(0:3),kq(0:3)
      external qdot
      real*16 qdot,delta
      common/momquad/p1q,p2q,p3q,p4q,kq
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      real*8 thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      mom(1,:)=p1q
      mom(2,:)=p2q
      mom(3,:)=-p3q
      mom(4,:)=-p4q
      mom(5,:)=-kq

      if(abs(x).le.1.d-11) then
         closetomassc=(0.q0,0.q0)
         return
      endif
      
      do i=1,5
         do j=i,5
            aux=qdot(mom(i,:)+mom(j,:),mom(i,:)+mom(j,:))*(1.q0,0.q0)
!            print*,'x',i,j,aux
            if(abs((aux-x)/x).le.1.q-9) then
               closetomassc=aux
               return
            endif
         enddo
      enddo
      
      masses(1)=ame*(1.q0,0.q0)
      masses(2)=ampi*(1.q0,0.q0)
      masses(3)=ammu*(1.q0,0.q0)

*---------------------------------------------------------
*-----    10^-5 is a reasonable  value      *------------- 
*---------------------------------------------------------

*--------When we have (p_i+k)^2 being k a photon we have that----------
*----------((p_i+k)^2-m_x^2)/m_x^2>= E_\gamma/E_X. --------------------

!!!!!!!!!HERE WE CHECK THAT E_\GAMMA/E_IS > 10^-5, IF NOT WE UPDATE!!!!!!!!!
      delta=1.q-4
      if(egmin/mom(1,0).le.delta) then
         delta=egmin/mom(1,0)
      endif
      
      do i=1,3
         if(abs((x-masses(i)**2)/x).le.delta) then
            closetomassc=masses(i)*masses(i)
            return
         endif
      enddo

      closetomassc=x
      return
      end



      function closetomasscd(x)
      implicit none
      double precision x
      double complex closetomasscd
      double precision ame,ammu,convfac,alpha,pi,s,s1,ampi
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      integer i,j
!!!!!!!!!!!!!!!COMMON TO CHECK IN XIJ REDUCTION
      common/momentacheck1/pin1br,pin2br,p1br,p2br,qbr
      double precision pin1br(0:3),pin2br(0:3),p1br(0:3)
      double precision p2br(0:3),qbr(0:3)
      external dot
      real*8 dot,delta
      real*8 mom(5,0:3)
      Complex*16 masses(3),aux
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      real*8 thmin,thmax,emin,zmax,egmin,thgmin,thgmax

      mom(1,:)=pin1br
      mom(2,:)=pin2br
      mom(3,:)=-p1br
      mom(4,:)=-p2br
      mom(5,:)=-qbr
      
      if(abs(x).le.1.d-11) then
         closetomasscd=(0.d0,0.d0)
         return
      endif
      
c$$$      do i=1,5
c$$$         do j=i,5
c$$$            aux=dot(mom(i,:)+mom(j,:),mom(i,:)+mom(j,:))
c$$$!            print*,'x',i,j,aux
c$$$            if(abs((aux-x)/x).le.1.d-9) then
c$$$               closetomasscd=aux
c$$$               return
c$$$            endif
c$$$         enddo
c$$$      enddo
      
      masses(1)=ame
      masses(2)=ampi
      masses(3)=ammu

*---------------------------------------------------------
*-----    10^-5 is a reasonable  value      *------------- 
*---------------------------------------------------------

*--------When we have (p_i+k)^2 being k a photon we have that----------
*----------((p_i+k)^2-m_x^2)/m_x^2>= E_\gamma/E_X. --------------------

!!!!!!!!!HERE WE CHECK THAT E_\GAMMA/E_IS > 10^-5, IF NOT WE UPDATE!!!!!!!!!
      delta=1.d-4
      if(egmin/mom(1,0).le.delta) then
         delta=egmin/mom(1,0)
      endif
      
      do i=1,2
         if(abs((x-masses(i)**2)/x).le.delta) then
            closetomasscd=masses(i)*masses(i)*(1.d0,0.d0)
            return
         endif
      enddo

      closetomasscd=x*(1.d0,0.d0)
      return
      end

      
*-------------------------------------------------------*
*     ROUTINES TO REDUCE TENSOR SCALAR INTEGRALS        *
*     WHERE PAVE ARE EVALUATED IN 4-PLE PRECISION       *
*                   C_I,C_IJ,D_I,D_IJ                   *
*-------------------------------------------------------*


      
*-----------ROUTINES TO REDUCE CI,CIJ-------------------*      
      subroutine Ci(k1,k2,m1q,m2q,m3q,C0m,C1m,C2m)
      implicit complex*32 (z)
#ifdef LOOPTOOLS
#include "looptools.h"
#endif
      double precision k1(0:3),k2(0:3),k3(0:3)
      double precision p1LT(0:3),p2LT(0:3),p3LT(0:3),p4LT(0:3)
      complex*32 m1q,m2q,m3q
      complex*32 C0m,C1m,C2m
      complex*32 B0AC,B0AB,B0BC
      complex*32 p1LT2,p2LT2,p3LT2
      integer i
      complex*32 re,im
      complex*32 k1k1,k2k2,k1k2
      
      external dotq
      real*16 dotq

      external closetomassc
      complex*32 closetomassc

      
      real*8 deltauv,lambda2,mu2dim,m2min
      common/loop/deltauv,lambda2,mu2dim,m2min

      common/lt/ltint,quadevents,ilt
      integer ltint,quadevents,ilt
      common/softphotonmass/phmass
      real*8 phmass
      real*16 zphmass2,zdeltauv
      
      re=(1.q0,0.q0)

      k1k1=closetomassc(dotq(k1,k1)*re)
      k2k2=closetomassc(dotq(k2,k2)*re)
      k1k2=closetomassc(dotq(k1,k2)*re)

      
      do i= 0,3
         p1LT(i)= k1(i)
         p2LT(i)= k2(i)-k1(i)
         p3LT(i)= k2(i)
      enddo

*     
      p1LT2= closetomassc(dotq(p1LT,p1LT)*re)
      p2LT2= closetomassc(dotq(p2LT,p2LT)*re)
      p3LT2= closetomassc(dotq(p3LT,p3LT)*re)
      
!!!!!!!!!!!!!!!!!!!!LT non sa nulla di k1k2, lo riscrive partendo da P2LT2-m1^2-m2^2, lo faccio anche qui

      k1k2 = -(p2LT2-k1k1-k2k2)/2.q0

!!!!!!!!!!!!!!!!!!!!
      
        if(ilt.eq.0) then
           ilt=1
           zdeltauv=deltauv
           zphmass2=phmass*phmass
#ifdef LOOPTOOLS
           call ltini
           call setdelta(zdeltauv)
           call setmudim(1.q2)
           call setlambda(zphmass2)
           call setwarningLT
#endif   
        endif
        

#ifdef LOOPTOOLS
      
      C0m   = C0C(p1LT2,p2LT2,p3LT2,m1q,m2q,m3q)

      B0AB = B0C(p1LT2,m1q,m2q)

      B0AC = B0C(p3LT2,m1q,m3q)

      B0BC = B0C(p2LT2,m2q,m3q)
#endif


      C1m = (-(B0AB*k1k2) + B0BC*(k1k2 - k2k2) + B0AC*k2k2 +
     >     C0m*(-(k1k1*k2k2) + k2k2*(-m1q + m2q) + k1k2*(k2k2 + m1q -
     >     m3q)))/(-2.q0*k1k2**2 + 2.q0*k1k1*k2k2)
      C2m =(-(B0AB*k1k1) + B0BC*(k1k1 - k1k2) + B0AC*k1k2 +
     >     C0m*k1k2*(-m1q + m2q) + C0m*k1k1*(-k1k2 + k2k2 + m1q - m3q))
     >     /(2.q0*(k1k2**2 - k1k1*k2k2))


      call clearcache
      
      return
      end

      subroutine Cij(k1,k2,m1q,m2q,m3q,C0m,Cim,Cijm)
      implicit complex*32 (z)
#ifdef LOOPTOOLS
#include "looptools.h"
#endif
      double precision k1(0:3),k2(0:3),k3(0:3)
      double precision p1LT(0:3),p2LT(0:3),p3LT(0:3),p4LT(0:3)
      complex*32 m1q,m2q,m3q
      complex*32 Cim(2)
      complex*32 Cijm(4)
      complex*32 C0m,C1m,C2m
      complex*32 C00m,C11m,C12m,C22m
      complex*32 B0AC,B0AB,B0BC
      complex*32 A01,A02,A03
      complex*32 p1LT2,p2LT2,p3LT2
      integer i
      complex*32 re,im
      real*16 k1k1,k2k2,k1k2

      external dotq
      real*16 dotq
      
      external closetomassc
      complex*32 closetomassc
      
      real*8 deltauv,lambda2,mu2dim,m2min
      common/loop/deltauv,lambda2,mu2dim,m2min

      
      common/lt/ltint,quadevents,ilt
      integer ltint,quadevents,ilt

      common/softphotonmass/phmass
      real*8 phmass
      
      real*16 zphmass2,zdeltauv
      re=(1.q0,0.q0)
      
      k1k1=closetomassc(dotq(k1,k1)*re)
      k2k2=closetomassc(dotq(k2,k2)*re)
      k1k2=closetomassc(dotq(k1,k2)*re)

      
      do i= 0,3
         p1LT(i)= k1(i)
         p2LT(i)= k2(i)-k1(i)
         p3LT(i)= k2(i)
      enddo

*     
      p1LT2= closetomassc(dotq(p1LT,p1LT)*re)
      p2LT2= closetomassc(dotq(p2LT,p2LT)*re)
      p3LT2= closetomassc(dotq(p3LT,p3LT)*re)
      
!!!!!!!!!!!!!!!!!!!!LT non sa nulla di k1k2, lo riscrive partendo da P2LT2-m1^2-m2^2, lo faccio anche qui

      k1k2 = -(p2LT2-k1k1-k2k2)/2.q0

      
        if(ilt.eq.0) then
           ilt=1
           zdeltauv=deltauv
           zphmass2=phmass*phmass
#ifdef LOOPTOOLS
           call ltini
           call setdelta(zdeltauv)
           call setmudim(1.q2)
           call setlambda(zphmass2)
           call setwarningLT
#endif   
        endif
!!!!!!!!!!!!!!!!!!!!
#ifdef LOOPTOOLS

      C0m   = C0C(p1LT2,p2LT2,p3LT2,m1q,m2q,m3q)

      B0AB = B0C(p1LT2,m1q,m2q)

      B0AC = B0C(p3LT2,m1q,m3q)

      B0BC = B0C(p2LT2,m2q,m3q)

      A01  = A0C(m1q)

      A02  = A0C(m2q)
      
      A03  = A0C(m3q)
#endif

!!!!!!!!!!!!!!!!!!!VETTORIALI      
      C1m = (-(B0AB*k1k2) + B0BC*(k1k2 - k2k2) + B0AC*k2k2 +
     >     C0m*(-(k1k1*k2k2) + k2k2*(-m1q + m2q) + k1k2*(k2k2 + m1q -
     >     m3q)))/(-2.q0*k1k2**2 + 2.q0*k1k1*k2k2)
      C2m =(-(B0AB*k1k1) + B0BC*(k1k1 - k1k2) + B0AC*k1k2 +
     >     C0m*k1k2*(-m1q + m2q) + C0m*k1k1*(-k1k2 + k2k2 + m1q - m3q))
     >     /(2.q0*(k1k2**2 - k1k1*k2k2))

      Cim(1) = C1m
      Cim(2) = C2m

c$$$!!!!!!!!!!!!!!Qui divido per k1k1 o k2k2
c$$$!!!!!!!!!!!!!!Se k1k1 o k2k2 è 0, lo metto a zero numerico in doppia (10^-15)
      if(k1k1.eq.(0.q0,0.q0)) then
         k1k1=1.q-30
      elseif(k2k2.eq.(0.q0,0.q0)) then
         k2k2=1.q-30
      endif
!!!!!!!!!!!!!!!!!!!TENSORIALI
        C00m = (B0AC*k1k2*m3q - B0AC*k1k2*k2k2*1q0 - 
     &    B0AC*k1k2*m1q*1q0 + 
     &    B0AB*k1k2*(k1k1 + m1q - m2q*1q0) - k1k2**2*2q0 + 
     &    k1k1*k2k2*2q0 + C1m*k1k1**2*k2k2*2q0 + 
     &    C1m*k1k1*k1k2*k2k2*2q0 - C2m*k1k1*k1k2*k2k2*2q0 + 
     &    C2m*k1k1*k2k2**2*2q0 + C1m*k1k1*k1k2*m1q*2q0 + 
     &    C1m*k1k1*k2k2*m1q*2q0 + C2m*k1k1*k2k2*m1q*2q0 - 
     &    C2m*k1k2*k2k2*m1q*2q0 - C1m*k1k1*k2k2*m2q*2q0 + 
     &    C2m*k1k2*k2k2*m2q*2q0 - C1m*k1k1*k1k2*m3q*2q0 - 
     &    C2m*k1k1*k2k2*m3q*2q0 + 
     &    B0BC*k1k2*(k2k2 + m2q - m3q*1q0 - k1k2*2q0) - 
     &    B0BC*k1k1*1q0*(k1k2 - k2k2*2q0) - 
     &    C1m*k1k1*k1k2**2*4q0 - C0m*k1k2**2*m1q*4q0 - 
     &    C1m*k1k2**2*m1q*4q0 + C0m*k1k1*k2k2*m1q*4q0 + 
     &    C1m*k1k2**2*m2q*4q0)/
     &  (-(k1k2**2*8q0) + k1k1*k2k2*8q0)

        C11m = -((1.25q-1*(B0AB*k1k1**3*k1k2*k2k2 + 
     &        B0AB*k1k1**2*k1k2*k2k2**2 + 
     &        B0AB*k1k1**2*k1k2*k2k2*m1q + 
     &        B0AB*k1k1*k1k2*k2k2**2*m1q + 
     &        B0BC*k1k1**2*k1k2*k2k2*m2q - 
     &        B0BC*k1k1**3*k1k2*k2k2*1q0 - 
     &        B0AB*k1k1**2*k1k2*k2k2*m2q*1q0 - 
     &        B0AB*k1k1*k1k2*k2k2**2*m2q*1q0 + 
     &        A03*k1k1*k1k2**3*2q0 + 
     &        B0AB*k1k1**2*k1k2**3*2q0 - 
     &        B0BC*k1k1**2*k1k2**3*2q0 - 
     &        A03*k1k1**2*k1k2*k2k2*2q0 - 
     &        A03*k1k1*k1k2**2*k2k2*2q0 - 
     &        k1k1**2*k1k2**2*k2k2*2q0 - 
     &        B0AB*k1k1**2*k1k2**2*k2k2*2q0 + 
     &        B0BC*k1k1**2*k1k2**2*k2k2*2q0 + 
     &        A02*k1k2**3*k2k2*2q0 + 
     &        B0AB*k1k1*k1k2**3*k2k2*2q0 - 
     &        B0BC*k1k1*k1k2**3*k2k2*2q0 - 
     &        A02*k1k1**2*k2k2**2*2q0 + 
     &        A03*k1k1**2*k2k2**2*2q0 + k1k1**3*k2k2**2*2q0 - 
     &        A02*k1k1*k1k2*k2k2**2*2q0 + 
     &        B0BC*k1k1**2*k1k2*k2k2**2*2q0 - 
     &        k1k1*k1k2**2*k2k2**2*2q0 + 
     &        k1k1**2*k2k2**3*2q0 + C2m*k1k1**3*k2k2**3*2q0 + 
     &        C1m*k1k1**2*k1k2*k2k2**3*2q0 + 
     &        C2m*k1k1**2*k2k2**4*2q0 + 
     &        B0AB*k1k1*k1k2**3*m1q*2q0 + 
     &        C1m*k1k1**3*k1k2*k2k2*m1q*2q0 - 
     &        B0AB*k1k1*k1k2**2*k2k2*m1q*2q0 + 
     &        B0AB*k1k2**3*k2k2*m1q*2q0 + 
     &        C2m*k1k1**3*k2k2**2*m1q*2q0 + 
     &        C2m*k1k1**2*k2k2**3*m1q*2q0 - 
     &        B0AB*k1k1*k1k2**3*m2q*2q0 + 
     &        B0BC*k1k1*k1k2**3*m2q*2q0 + 
     &        B0AB*k1k1*k1k2**2*k2k2*m2q*2q0 - 
     &        B0AB*k1k2**3*k2k2*m2q*2q0 + 
     &        B0BC*k1k1**2*k2k2**2*m2q*2q0 + 
     &        A01*k1k2*(k1k1*k2k2 - k1k2**2*1q0)*2q0*
     &         (k1k1 + k2k2 - k1k2*2q0) - 
     &        B0AC*k1k1**2*k1k2*k2k2**2*3q0 - 
     &        B0AC*k1k1*k1k2*k2k2**3*3q0 + 
     &        B0BC*k1k1*k1k2*k2k2**3*3q0 - 
     &        B0AC*k1k1**2*k1k2*k2k2*m1q*3q0 - 
     &        B0AC*k1k1*k1k2*k2k2**2*m1q*3q0 + 
     &        B0BC*k1k1*k1k2*k2k2**2*m2q*3q0 - 
     &        A02*k1k2**4*4q0 - B0AB*k1k1*k1k2**4*4q0 + 
     &        B0BC*k1k1*k1k2**4*4q0 + k1k1*k1k2**3*k2k2*4q0 - 
     &        k1k1**2*k1k2*k2k2**2*4q0 + 
     &        C1m*k1k1*k1k2**3*k2k2**2*4q0 + 
     &        C2m*k1k1*k1k2**2*k2k2**3*4q0 + 
     &        C1m*k1k1**2*k1k2**3*m1q*4q0 - 
     &        B0AB*k1k2**4*m1q*4q0 - 
     &        C0m*k1k1**2*k1k2**2*k2k2*m1q*4q0 + 
     &        C2m*k1k1**2*k1k2**2*k2k2*m1q*4q0 + 
     &        C0m*k1k1**3*k2k2**2*m1q*4q0 - 
     &        C0m*k1k1*k1k2**2*k2k2**2*m1q*4q0 + 
     &        C0m*k1k1**2*k2k2**3*m1q*4q0 + 
     &        B0AB*k1k2**4*m2q*4q0 + 
     &        A02*k1k1*k1k2**2*k2k2*6q0 + 
     &        C1m*k1k1**4*k2k2**2*6q0 - 
     &        C2m*k1k1**3*k1k2*k2k2**2*6q0 + 
     &        B0AC*k1k1*k1k2**2*k2k2**2*6q0 - 
     &        B0BC*k1k1*k1k2**2*k2k2**2*6q0 + 
     &        C1m*k1k1**3*k2k2**3*6q0 + 
     &        B0AC*k1k1*k1k2**2*k2k2*m1q*6q0 + 
     &        C1m*k1k1**3*k2k2**2*m1q*6q0 + 
     &        C1m*k1k1**2*k2k2**3*m1q*6q0 - 
     &        C2m*k1k1*k1k2*k2k2**3*m1q*6q0 - 
     &        C1m*k1k1**3*k2k2**2*m2q*6q0 + 
     &        C2m*k1k1**2*k1k2*k2k2**2*m2q*6q0 - 
     &        C1m*k1k1**2*k2k2**3*m2q*6q0 + 
     &        C2m*k1k1*k1k2*k2k2**3*m2q*6q0 - 
     &        C1m*k1k1*k1k2**4*k2k2*8q0 - 
     &        C2m*k1k1*k1k2**3*k2k2**2*8q0 - 
     &        C1m*k1k1*k1k2**4*m1q*8q0 + 
     &        C0m*k1k1*k1k2**3*k2k2*m1q*8q0 - 
     &        C2m*k1k1*k1k2**3*k2k2*m1q*8q0 - 
     &        C0m*k1k1**2*k1k2*k2k2**2*m1q*8q0 - 
     &        B0BC*k1k1*k1k2**2*k2k2*m2q*8q0 - 
     &        k1k1*m3q*1q0*
     &         ((k1k1 + k2k2 - k1k2*2q0)*
     &            (C2m*k1k1*k2k2**2*2q0 + C1m*k1k2**3*4q0 + 
     &              k1k2*k2k2*
     &               (C1m*k1k1*2q0 - B0AC*3q0 + C2m*k1k2*4q0))
     &             + B0BC*
     &            (k1k2**3*2q0 + k1k1*k2k2**2*2q0 + 
     &              k1k2*k2k2*(k1k1 + k2k2*3q0) - 
     &              k1k2**2*k2k2*8q0)) - 
     &        C1m*k1k1**3*k1k2*k2k2**2*1.q1 - 
     &        C2m*k1k1**2*k1k2*k2k2**3*1.q1 - 
     &        C1m*k1k1**2*k1k2*k2k2**2*m1q*1.q1 - 
     &        C2m*k1k1**2*k1k2*k2k2**2*m1q*1.q1 - 
     &        C1m*k1k1**3*k1k2**2*k2k2*1.2q1 - 
     &        C1m*k1k1*k1k2**2*k2k2**2*m1q*1.2q1 + 
     &        C1m*k1k1**2*k1k2**2*k2k2*m2q*1.2q1 + 
     &        C1m*k1k1**2*k1k2*k2k2**2*m2q*1.2q1 + 
     &        C1m*k1k1*k1k2**2*k2k2**2*m2q*1.2q1 - 
     &        C2m*k1k1*k1k2**2*k2k2**2*m2q*1.2q1 - 
     &        C1m*k1k1**2*k1k2**2*k2k2**2*1.6q1 + 
     &        C2m*k1k1**2*k1k2**2*k2k2**2*1.6q1 - 
     &        C1m*k1k1**2*k1k2**2*k2k2*m1q*1.6q1 + 
     &        C2m*k1k1*k1k2**2*k2k2**2*m1q*1.6q1 - 
     &        C1m*k1k1*k1k2**3*k2k2*m2q*2.4q1 + 
     &        C1m*k1k1**2*k1k2**3*k2k2*2.8q1 + 
     &        C1m*k1k1*k1k2**3*k2k2*m1q*2.8q1))/
     &    (k1k1*(k1k2**2 - k1k1*k2k2*1q0)**2*
     &      (k1k1 + k2k2 - k1k2*2q0)))
        
        C12m = (1.25q-1*(B0BC*k1k2**2*k2k2**2 + 
     &      B0BC*k1k2**2*k2k2*m2q - 
     &      B0AC*k1k1*k1k2**2*k2k2*1q0 - 
     &      B0AC*k1k2**2*k2k2**2*1q0 - 
     &      B0AC*k1k1*k1k2**2*m1q*1q0 - 
     &      B0AC*k1k2**2*k2k2*m1q*1q0 + 
     &      A03*k1k1*k1k2**2*2q0 - A02*k1k2**3*2q0 - 
     &      A03*k1k2**3*2q0 - k1k1*k1k2**3*2q0 - 
     &      A03*k1k1**2*k2k2*2q0 + A02*k1k1*k1k2*k2k2*2q0 + 
     &      A03*k1k1*k1k2*k2k2*2q0 + k1k1**2*k1k2*k2k2*2q0 + 
     &      C1m*k1k1**3*k1k2*k2k2*2q0 + 
     &      A02*k1k2**2*k2k2*2q0 - 
     &      B0BC*k1k1*k1k2**2*k2k2*2q0 + 
     &      C1m*k1k1**2*k1k2**2*k2k2*2q0 - 
     &      C2m*k1k1**2*k1k2**2*k2k2*2q0 - k1k2**3*k2k2*2q0 + 
     &      B0AC*k1k2**3*k2k2*2q0 - B0BC*k1k2**3*k2k2*2q0 - 
     &      A02*k1k1*k2k2**2*2q0 - B0AC*k1k1**2*k2k2**2*2q0 + 
     &      B0BC*k1k1**2*k2k2**2*2q0 + 
     &      k1k1*k1k2*k2k2**2*2q0 + 
     &      C1m*k1k1**2*k1k2*k2k2**2*2q0 - 
     &      B0AC*k1k1*k2k2**3*2q0 + B0BC*k1k1*k2k2**3*2q0 + 
     &      B0AC*k1k2**3*m1q*2q0 - 
     &      B0AC*k1k1**2*k2k2*m1q*2q0 + 
     &      C1m*k1k1**2*k1k2*k2k2*m1q*2q0 + 
     &      C1m*k1k1*k1k2**2*k2k2*m1q*2q0 - 
     &      B0AC*k1k1*k2k2**2*m1q*2q0 + 
     &      C1m*k1k1*k1k2*k2k2**2*m1q*2q0 - 
     &      C2m*k1k2**2*k2k2**2*m1q*2q0 - 
     &      B0BC*k1k1*k1k2*k2k2*m2q*2q0 - 
     &      C1m*k1k1**2*k1k2*k2k2*m2q*2q0 + 
     &      C2m*k1k1*k1k2**2*k2k2*m2q*2q0 + 
     &      B0BC*k1k1*k2k2**2*m2q*2q0 - 
     &      C1m*k1k1*k1k2*k2k2**2*m2q*2q0 + 
     &      C2m*k1k2**2*k2k2**2*m2q*2q0 + 
     &      A01*(k1k1*k2k2 - k1k2**2*1q0)*2q0*
     &       (k1k1 + k2k2 - k1k2*2q0) + 
     &      B0AB*k1k1**2*k1k2**2*3q0 - 
     &      B0BC*k1k1**2*k1k2**2*3q0 + 
     &      B0AB*k1k1*k1k2**2*k2k2*3q0 + 
     &      B0AB*k1k1*k1k2**2*m1q*3q0 + 
     &      B0AB*k1k2**2*k2k2*m1q*3q0 - 
     &      B0AB*k1k1*k1k2**2*m2q*3q0 + 
     &      B0BC*k1k1*k1k2**2*m2q*3q0 - 
     &      B0AB*k1k2**2*k2k2*m2q*3q0 + k1k2**4*4q0 - 
     &      k1k1*k1k2**2*k2k2*4q0 + 
     &      C2m*k1k1*k1k2**3*k2k2*4q0 - 
     &      C2m*k1k1**3*k2k2**2*4q0 + 
     &      B0AC*k1k1*k1k2*k2k2**2*4q0 - 
     &      B0BC*k1k1*k1k2*k2k2**2*4q0 - 
     &      C2m*k1k1**2*k2k2**3*4q0 - 
     &      C0m*k1k1*k1k2**3*m1q*4q0 + 
     &      B0AC*k1k1*k1k2*k2k2*m1q*4q0 + 
     &      C0m*k1k1**2*k1k2*k2k2*m1q*4q0 - 
     &      C0m*k1k2**3*k2k2*m1q*4q0 + 
     &      C2m*k1k2**3*k2k2*m1q*4q0 - 
     &      C2m*k1k1**2*k2k2**2*m1q*4q0 + 
     &      C0m*k1k1*k1k2*k2k2**2*m1q*4q0 - 
     &      C2m*k1k1*k2k2**3*m1q*4q0 - B0BC*k1k2**3*m2q*4q0 + 
     &      C1m*k1k1*k1k2**2*k2k2*m2q*4q0 - 
     &      C2m*k1k2**3*k2k2*m2q*4q0 + 
     &      C2m*k1k1**2*k2k2**2*m2q*4q0 + 
     &      C2m*k1k1*k2k2**3*m2q*4q0 - 
     &      B0AB*k1k1*k1k2**3*6q0 + B0BC*k1k1*k1k2**3*6q0 + 
     &      C1m*k1k1*k1k2**2*k2k2**2*6q0 + 
     &      C2m*k1k1*k1k2*k2k2**3*6q0 + 
     &      C1m*k1k1**2*k1k2**2*m1q*6q0 - 
     &      B0AB*k1k2**3*m1q*6q0 + 
     &      C2m*k1k1**2*k1k2*k2k2*m1q*6q0 + 
     &      B0AB*k1k2**3*m2q*6q0 + 
     &      m3q*(B0AC*(k1k1 + k2k2 - k1k2*2q0)*
     &          (k1k2**2 + k1k1*k2k2*2q0) - 
     &         B0BC*1q0*(k1k1*
     &             (-(k1k2*k2k2*2q0) + k2k2**2*2q0 + 
     &               k1k2**2*3q0) + k1k2**2*(k2k2 - k1k2*4q0))
     &           - k1k1*k1k2*(C1m*k1k2 + C2m*k2k2)*
     &          (k1k1 + k2k2 - k1k2*2q0)*6q0) - 
     &      C1m*k1k1**2*k1k2**3*8q0 + C0m*k1k2**4*m1q*8q0 - 
     &      C0m*k1k1*k1k2**2*k2k2*m1q*8q0 - 
     &      C1m*k1k2**3*k2k2*m1q*8q0 + 
     &      C1m*k1k1*k1k2**3*m2q*8q0 + 
     &      C1m*k1k2**3*k2k2*m2q*8q0 - 
     &      C2m*k1k1*k1k2*k2k2**2*m2q*8q0 + 
     &      C2m*k1k1**2*k1k2*k2k2**2*1.4q1 - 
     &      C2m*k1k1*k1k2**2*k2k2**2*1.4q1 - 
     &      C2m*k1k1*k1k2**2*k2k2*m1q*1.4q1 + 
     &      C2m*k1k1*k1k2*k2k2**2*m1q*1.4q1 + 
     &      C1m*k1k1*k1k2**4*1.6q1 + C1m*k1k2**4*m1q*1.6q1 - 
     &      C1m*k1k2**4*m2q*1.6q1 - 
     &      C1m*k1k1*k1k2**3*k2k2*2.q1 - 
     &      C1m*k1k1*k1k2**3*m1q*2.q1))/
     &  ((k1k2**2 - k1k1*k2k2*1q0)**2*
     &    (k1k1 + k2k2 - k1k2*2q0))
        C22m = -((1.25q-1*(A03*k1k1*k1k2**3*2q0 - 
     &        A03*k1k1**2*k1k2*k2k2*2q0 - 
     &        A02*k1k1*k1k2**2*k2k2*2q0 - 
     &        k1k1**2*k1k2**2*k2k2*2q0 + 
     &        A02*k1k2**3*k2k2*2q0 + 
     &        B0AC*k1k1*k1k2**3*k2k2*2q0 - 
     &        B0BC*k1k1*k1k2**3*k2k2*2q0 + 
     &        A02*k1k1**2*k2k2**2*2q0 - 
     &        A03*k1k1**2*k2k2**2*2q0 + k1k1**3*k2k2**2*2q0 + 
     &        C1m*k1k1**4*k2k2**2*2q0 - 
     &        A02*k1k1*k1k2*k2k2**2*2q0 + 
     &        B0BC*k1k1**2*k1k2*k2k2**2*2q0 + 
     &        C1m*k1k1**3*k1k2*k2k2**2*2q0 - 
     &        k1k1*k1k2**2*k2k2**2*2q0 + 
     &        B0AC*k1k2**3*k2k2**2*2q0 - 
     &        B0BC*k1k2**3*k2k2**2*2q0 + 
     &        k1k1**2*k2k2**3*2q0 + C1m*k1k1**3*k2k2**3*2q0 + 
     &        B0AC*k1k1*k1k2**3*m1q*2q0 + 
     &        B0AC*k1k2**3*k2k2*m1q*2q0 + 
     &        C1m*k1k1**3*k2k2**2*m1q*2q0 + 
     &        C1m*k1k1**2*k1k2*k2k2**2*m1q*2q0 + 
     &        C1m*k1k1**2*k2k2**3*m1q*2q0 - 
     &        B0BC*k1k2**3*k2k2*m2q*2q0 - 
     &        B0BC*k1k1**2*k2k2**2*m2q*2q0 - 
     &        C1m*k1k1**3*k2k2**2*m2q*2q0 - 
     &        C1m*k1k1**2*k2k2**3*m2q*2q0 + 
     &        A01*k1k2*(k1k1*k2k2 - k1k2**2*1q0)*2q0*
     &         (k1k1 + k2k2 - k1k2*2q0) + 
     &        B0AB*k1k1**3*k1k2*k2k2*3q0 - 
     &        B0BC*k1k1**3*k1k2*k2k2*3q0 + 
     &        B0AB*k1k1**2*k1k2*k2k2**2*3q0 + 
     &        B0AB*k1k1**2*k1k2*k2k2*m1q*3q0 + 
     &        B0AB*k1k1*k1k2*k2k2**2*m1q*3q0 - 
     &        B0AB*k1k1**2*k1k2*k2k2*m2q*3q0 + 
     &        B0BC*k1k1**2*k1k2*k2k2*m2q*3q0 - 
     &        B0AB*k1k1*k1k2*k2k2**2*m2q*3q0 - 
     &        A03*k1k2**4*4q0 + k1k1*k1k2**3*k2k2*4q0 + 
     &        C2m*k1k1**2*k1k2**3*k2k2*4q0 - 
     &        B0AC*k1k2**4*k2k2*4q0 + B0BC*k1k2**4*k2k2*4q0 - 
     &        k1k1**2*k1k2*k2k2**2*4q0 + 
     &        C2m*k1k1*k1k2**3*k2k2**2*4q0 - 
     &        B0AC*k1k2**4*m1q*4q0 - 
     &        C0m*k1k1**2*k1k2**2*k2k2*m1q*4q0 + 
     &        C2m*k1k1*k1k2**3*k2k2*m1q*4q0 + 
     &        C0m*k1k1**3*k2k2**2*m1q*4q0 - 
     &        C0m*k1k1*k1k2**2*k2k2**2*m1q*4q0 + 
     &        C2m*k1k2**3*k2k2**2*m1q*4q0 + 
     &        C0m*k1k1**2*k2k2**3*m1q*4q0 - 
     &        B0BC*k1k1*k1k2**2*k2k2*m2q*4q0 - 
     &        C2m*k1k1*k1k2**3*k2k2*m2q*4q0 + 
     &        C1m*k1k1**2*k1k2*k2k2**2*m2q*4q0 - 
     &        C2m*k1k2**3*k2k2**2*m2q*4q0 - 
     &        B0AC*k1k1**2*k1k2*k2k2**2*5q0 - 
     &        B0AC*k1k1*k1k2*k2k2**3*5q0 + 
     &        B0BC*k1k1*k1k2*k2k2**3*5q0 - 
     &        B0AC*k1k1**2*k1k2*k2k2*m1q*5q0 - 
     &        B0AC*k1k1*k1k2*k2k2**2*m1q*5q0 + 
     &        B0BC*k1k1*k1k2*k2k2**2*m2q*5q0 + 
     &        A03*k1k1*k1k2**2*k2k2*6q0 - 
     &        B0AB*k1k1**2*k1k2**2*k2k2*6q0 + 
     &        B0BC*k1k1**2*k1k2**2*k2k2*6q0 + 
     &        C2m*k1k1**3*k2k2**3*6q0 + 
     &        C1m*k1k1**2*k1k2*k2k2**3*6q0 + 
     &        C2m*k1k1**2*k2k2**4*6q0 + 
     &        C1m*k1k1**3*k1k2*k2k2*m1q*6q0 - 
     &        B0AB*k1k1*k1k2**2*k2k2*m1q*6q0 + 
     &        C2m*k1k1**3*k2k2**2*m1q*6q0 + 
     &        C2m*k1k1**2*k2k2**3*m1q*6q0 + 
     &        B0AB*k1k1*k1k2**2*k2k2*m2q*6q0 + 
     &        m3q*(B0AC*k1k2*(k1k1 + k2k2 - k1k2*2q0)*
     &            (-(k1k2**2*2q0) + k1k1*k2k2*5q0) + 
     &           B0BC*k2k2*
     &            (k1k2**3*2q0 + 
     &              k1k1**2*(k2k2*2q0 - k1k2*3q0) + 
     &              k1k1*k1k2*(k1k2*4q0 - k2k2*5q0)) - 
     &           k1k1**2*k2k2*(C1m*k1k2 + C2m*k2k2)*
     &            (k1k1 + k2k2 - k1k2*2q0)*6q0) - 
     &        C1m*k1k1**3*k1k2**2*k2k2*8q0 - 
     &        C2m*k1k1*k1k2**4*k2k2*8q0 + 
     &        C0m*k1k1*k1k2**3*k2k2*m1q*8q0 - 
     &        C2m*k1k2**4*k2k2*m1q*8q0 - 
     &        C0m*k1k1**2*k1k2*k2k2**2*m1q*8q0 - 
     &        C1m*k1k1*k1k2**2*k2k2**2*m1q*8q0 + 
     &        C1m*k1k1**2*k1k2**2*k2k2*m2q*8q0 + 
     &        C2m*k1k2**4*k2k2*m2q*8q0 + 
     &        C1m*k1k1*k1k2**2*k2k2**2*m2q*8q0 - 
     &        C2m*k1k1**3*k1k2*k2k2**2*1.q1 + 
     &        B0AC*k1k1*k1k2**2*k2k2**2*1.q1 - 
     &        B0BC*k1k1*k1k2**2*k2k2**2*1.q1 + 
     &        B0AC*k1k1*k1k2**2*k2k2*m1q*1.q1 - 
     &        C2m*k1k1*k1k2*k2k2**3*m1q*1.q1 + 
     &        C2m*k1k1**2*k1k2*k2k2**2*m2q*1.q1 + 
     &        C2m*k1k1*k1k2*k2k2**3*m2q*1.q1 + 
     &        C1m*k1k1**2*k1k2**3*k2k2*1.6q1 + 
     &        C1m*k1k1*k1k2**3*k2k2*m1q*1.6q1 - 
     &        C1m*k1k1*k1k2**3*k2k2*m2q*1.6q1 - 
     &        C1m*k1k1**2*k1k2**2*k2k2**2*2.q1 + 
     &        C2m*k1k1**2*k1k2**2*k2k2**2*2.q1 - 
     &        C1m*k1k1**2*k1k2**2*k2k2*m1q*2.q1 + 
     &        C2m*k1k1*k1k2**2*k2k2**2*m1q*2.q1 - 
     &        C2m*k1k1*k1k2**2*k2k2**2*m2q*2.q1 - 
     &        C2m*k1k1**2*k1k2*k2k2**3*2.2q1 - 
     &        C2m*k1k1**2*k1k2*k2k2**2*m1q*2.2q1))/
     &    (k2k2*(k1k2**2 - k1k1*k2k2*1q0)**2*
     &      (k1k1 + k2k2 - k1k2*2q0)))
                  
        Cijm(1) = C00m
        Cijm(2) = C11m
        Cijm(3) = C12m
        Cijm(4) = C22m

        call clearcache
      return
      end
      
*-----------ROUTINES TO REDUCE DI,DIJ-------------------*      
      subroutine Di(k1,k2,k3,m1q,m2q,m3q,m4q,D0m,Dim)
      implicit complex*32 (D,z)
#ifdef LOOPTOOLS
#include "looptools.h"
#endif
      double precision k1(0:3),k2(0:3),k3(0:3)
      double precision p1LT(0:3),p2LT(0:3),p3LT(0:3),p4LT(0:3)
      complex*32 m1q,m2q,m3q,m4q
      complex*32 C0ABC,C0ABD,C0ACD,C0BCD
      complex*32 D0m
      complex*32 Dim(3)
      complex*32 p1LT2,p2LT2,p3LT2,p4LT2
      complex*32 x23LT,x12LT
      integer i
      complex*32 re,im
      complex*32 k1k1,k1k2,k1k3
      complex*32 k2k2,k2k3
      complex*32 k3k3
      
      external dotq
      real*16 dotq
      
      real*8 deltauv,lambda2,mu2dim,m2min
      common/loop/deltauv,lambda2,mu2dim,m2min

      external closetomassc
      complex*32 closetomassc
      
      common/lt/ltint,quadevents,ilt
      integer ltint,quadevents,ilt
      common/softphotonmass/phmass
      real*8 phmass
      
      real*16 zphmass2,zdeltauv
      re=(1.q0,0.q0)

      k1k1=closetomassc(dotq(k1,k1)*re)
      k2k2=closetomassc(dotq(k2,k2)*re)
      k3k3=closetomassc(dotq(k3,k3)*re)

      do i= 0,3
         p1LT(i)= k1(i)
         p2LT(i)= k2(i)-k1(i)
         p3LT(i)= k3(i)-k2(i)
         p4LT(i)= -k3(i)
      enddo

      m1q   = closetomassc(m1q)
      m2q   = closetomassc(m2q)
      m3q   = closetomassc(m3q)
      m4q   = closetomassc(m4q)
!      m5q   = closetomassc(m5*re)

*     
      p1LT2= closetomassc(dotq(p1LT,p1LT)*re)
      p2LT2= closetomassc(dotq(p2LT,p2LT)*re)
      p3LT2= closetomassc(dotq(p3LT,p3LT)*re)
      p4LT2= closetomassc(dotq(p4LT,p4LT)*re)
      
      x12LT= closetomassc(dotq(k2,k2)*re)
      x23LT= closetomassc(dotq(k3-k1,k3-k1)*re)

!!!!!!!!!!!!!!!!!!!!LT non sa nulla di k1k2, lo riscrive partendo da P2LT2-m1^2-m2^2, lo faccio anche qui

      k1k2 = -(p2LT2-k1k1-k2k2)/2.q0
      k1k3 = -(x23LT-k1k1-k3k3)/2.q0
      k2k3 = -(p3LT2-k3k3-k2k2)/2.q0
!!!!!!!!!!!!!!!!!!!!

      
        if(ilt.eq.0) then
           ilt=1
           zdeltauv=deltauv
           zphmass2=phmass*phmass
#ifdef LOOPTOOLS
           call ltini
           call setdelta(zdeltauv)
           call setmudim(1.q2)
           call setlambda(zphmass2)
           call setwarningLT
#endif   
        endif
        
#ifdef LOOPTOOLS
      
      D0m   = D0C(p1LT2,p2LT2,p3LT2,p4LT2,x12LT,x23LT,m1q,m2q,m3q,m4q)
      if(m1q.eq.(0.q0,0.q0)) then
         C0ABC = C0C(p1LT2,p2LT2,x12LT,m1q,m2q,m3q)
      elseif(m2q.eq.(0.q0,0.q0)) then
         C0ABC = C0C(p1LT2,x12LT,p2LT2,m2q,m1q,m3q)
      elseif(m3q.eq.(0.q0,0.q0)) then
         C0ABC = C0C(x12LT,p1LT2,p2LT2,m3q,m1q,m2q)
      else
         C0ABC = C0C(p1LT2,p2LT2,x12LT,m1q,m2q,m3q)
      endif

      if(m1q.eq.(0.q0,0.q0)) then
         C0ACD = C0C(x12LT,p3LT2,p4LT2,m1q,m3q,m4q)
      elseif(m3q.eq.(0.q0,0.q0)) then
         C0ACD = C0C(x12LT,p4LT2,p3LT2,m3q,m1q,m4q)
      elseif(m4q.eq.(0.q0,0.q0)) then
         C0ACD = C0C(p4LT2,x12LT,p3LT2,m4q,m1q,m3q)
      else
         C0ACD = C0C(x12LT,p3LT2,p4LT2,m1q,m3q,m4q)
      endif

      if(m1q.eq.(0.q0,0.q0)) then
         C0ABD = C0C(p1LT2,x23LT,p4LT2,m1q,m2q,m4q)
      elseif(m2q.eq.(0.q0,0.q0)) then
         C0ABD = C0C(p1LT2,p4LT2,x23LT,m2q,m1q,m4q)
      elseif(m4q.eq.(0.q0,0.q0)) then
         C0ABD = C0C(p4LT2,p1LT2,x23LT,m4q,m1q,m2q)
      else
         C0ABD = C0C(p1LT2,x23LT,p4LT2,m1q,m2q,m4q)
      endif

      if(m2q.eq.(0.q0,0.q0)) then
         C0BCD = C0C(p2LT2,p3LT2,x23LT,m2q,m3q,m4q)
      elseif(m3q.eq.(0.q0,0.q0)) then
         C0BCD = C0C(p2LT2,x23LT,p3LT2,m3q,m2q,m4q)
      elseif(m4q.eq.(0.q0,0.q0)) then
         C0BCD = C0C(x23LT,p2LT2,p3LT2,m4q,m2q,m3q)
      else
         C0BCD = C0C(p2LT2,p3LT2,x23LT,m2q,m3q,m4q)
      endif
      

#endif

      Dim(1) = -((5.q-1*(-(1q0*(k1k3*k2k3 - k1k2*k3k3*1q0)*
     &           (k2k3*(C0ACD - C0BCD*1q0 - 
     &                D0m*1q0*(k1k1 + m1q - m2q*1q0)) + 
     &             k1k3*(C0BCD - C0ABD*1q0 + 
     &                D0m*(k2k2 + m1q - m3q*1q0)))) + 
     &        (k1k3*k2k2 - k1k2*k2k3*1q0)*
     &         (k3k3*(C0ACD - C0BCD*1q0 - 
     &              D0m*1q0*(k1k1 + m1q - m2q*1q0)) + 
     &           k1k3*(C0BCD - C0ABC*1q0 + 
     &              D0m*(k3k3 + m1q - m4q*1q0)))))/
     &    (k1k3*(k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &        k1k1*(k2k3**2 - k2k2*k3k3*1q0) - 
     &        k1k2*k1k3*k2k3*2q0)))
        Dim(2) = (5.q-1*(C0ABD*k1k3**2 + D0m*k1k1*k1k3*k2k3 + 
     &      C0ACD*k1k2*k3k3 + D0m*k1k2*k1k3*k3k3 + 
     &      D0m*k1k1*k2k2*k3k3 + D0m*k1k2*k1k3*m1q + 
     &      D0m*k1k3*k2k3*m1q + D0m*k1k1*k3k3*m1q + 
     &      D0m*k1k2*k3k3*m2q + D0m*k1k3**2*m3q + 
     &      D0m*k1k1*k2k3*m4q - D0m*k1k3**2*k2k2*1q0 - 
     &      C0ACD*k1k3*k2k3*1q0 - C0ABD*k1k1*k3k3*1q0 - 
     &      D0m*k1k1*k1k2*k3k3*1q0 - D0m*k1k1*k2k3*k3k3*1q0 - 
     &      D0m*k1k3**2*m1q*1q0 - D0m*k1k1*k2k3*m1q*1q0 - 
     &      D0m*k1k2*k3k3*m1q*1q0 - D0m*k1k3*k2k3*m2q*1q0 - 
     &      D0m*k1k1*k3k3*m3q*1q0 - D0m*k1k2*k1k3*m4q*1q0 + 
     &      C0ABC*(k1k1*k2k3 - k1k2*k1k3*1q0) + 
     &      C0BCD*(-(k1k2*k3k3*1q0) + 
     &         k1k3*(k1k2 + k2k3 - k1k3*1q0) + 
     &         k1k1*(k3k3 - k2k3*1q0))))/
     &  (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &    k1k1*(k2k3**2 - k2k2*k3k3*1q0) - k1k2*k1k3*k2k3*2q0)
        Dim(3) = -((5.q-1*(C0ABD*k1k2*k1k3 + D0m*k1k1*k1k3*k2k2 + 
     &        C0ACD*k1k2*k2k3 + D0m*k1k1*k2k2*k2k3 + 
     &        D0m*k1k2**2*k3k3 + D0m*k1k2**2*m1q + 
     &        D0m*k1k3*k2k2*m1q + D0m*k1k1*k2k3*m1q + 
     &        D0m*k1k2*k2k3*m2q + D0m*k1k2*k1k3*m3q + 
     &        D0m*k1k1*k2k2*m4q - C0ACD*k1k3*k2k2*1q0 - 
     &        D0m*k1k2*k1k3*k2k2*1q0 - C0ABD*k1k1*k2k3*1q0 - 
     &        D0m*k1k1*k1k2*k2k3*1q0 - D0m*k1k1*k2k2*k3k3*1q0 - 
     &        D0m*k1k2*k1k3*m1q*1q0 - D0m*k1k1*k2k2*m1q*1q0 - 
     &        D0m*k1k2*k2k3*m1q*1q0 - D0m*k1k3*k2k2*m2q*1q0 - 
     &        D0m*k1k1*k2k3*m3q*1q0 - D0m*k1k2**2*m4q*1q0 + 
     &        C0ABC*(k1k1*k2k2 - k1k2**2*1q0) + 
     &        C0BCD*(k1k2**2 + k1k3*k2k2 + k1k1*k2k3 - 
     &           k1k1*k2k2*1q0 - k1k2*(k1k3 + k2k3)*1q0)))/
     &    (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &      k1k1*(k2k3**2 - k2k2*k3k3*1q0) - 
     &      k1k2*k1k3*k2k3*2q0))

        call clearcache
      return
      end
      
      subroutine Dij(k1,k2,k3,m1q,m2q,m3q,m4q,D0m,Dim,Dijm)
      implicit complex*32 (z)
#ifdef LOOPTOOLS
#include "looptools.h"
#endif
      double precision k1(0:3),k2(0:3),k3(0:3)
      double precision p1LT(0:3),p2LT(0:3),p3LT(0:3),p4LT(0:3)
      complex*32 m1q,m2q,m3q,m4q
      complex*32 C0ABC,C0ABD,C0ACD,C0BCD
      complex*32 C1ABC,C2ABC
      complex*32 C1ABD,C2ABD
      complex*32 C1ACD,C2ACD
      complex*32 C1BCD,C2BCD
      
      complex*32 D0m
      complex*32 D1m,D2m,D3m
      complex*32 Dim(3),Dijm(7)

      complex*32 bc0,bc1,bc2,bc3,bc4,bc5,bc6
      complex*32 d12
      complex*32 p1LT2,p2LT2,p3LT2,p4LT2
      complex*32 x23LT,X12LT
      integer i
      complex*32 re,im
      complex*32 k1k1,k1k2,k1k3
      complex*32 k2k2,k2k3
      complex*32 k3k3
      
      external dotq
      real*16 dotq

      external closetomassc
      complex*32 closetomassc
      
      real*8 deltauv,lambda2,mu2dim,m2min
      common/loop/deltauv,lambda2,mu2dim,m2min
      common/softphotonmass/phmass
      real*8 phmass

      
      common/lt/ltint,quadevents,ilt
      integer ltint,quadevents,ilt
      
      real*16 zphmass2,zdeltauv

      re=(1.q0,0.q0)

      k1k1=closetomassc(dotq(k1,k1)*re)
      k2k2=closetomassc(dotq(k2,k2)*re)
      k3k3=closetomassc(dotq(k3,k3)*re)

      do i= 0,3
         p1LT(i)= k1(i)
         p2LT(i)= k2(i)-k1(i)
         p3LT(i)= k3(i)-k2(i)
         p4LT(i)= -k3(i)
      enddo

*     

c$$$      print*,'MOMENTI RED'
c$$$      print*,k1
c$$$      print*,k2
c$$$      print*,k3
      p1LT2= closetomassc(dotq(p1LT,p1LT)*re)
      p2LT2= closetomassc(dotq(p2LT,p2LT)*re)
      p3LT2= closetomassc(dotq(p3LT,p3LT)*re)
      p4LT2= closetomassc(dotq(p4LT,p4LT)*re)
      
      x12LT= closetomassc(dotq(k2,k2)*re)
      x23LT= closetomassc(dotq(k3-k1,k3-k1)*re)
      
!!!!!!!!!!!!!!!!!!!!LT non sa nulla di k1k2, lo riscrive partendo da P2LT2-m1^2-m2^2, lo faccio anche qui

      k1k2 = -(p2LT2-k1k1-k2k2)/2.q0
      k1k3 = -(x23LT-k1k1-k3k3)/2.q0
      k2k3 = -(p3LT2-k3k3-k2k2)/2.q0
!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!Vettoriali e scalare  4 pt      

        if(ilt.eq.0) then
           ilt=1
           zdeltauv=deltauv
           zphmass2=phmass*phmass
#ifdef LOOPTOOLS
           call ltini
           call setdelta(zdeltauv)
           call setmudim(1.q2)
           call setlambda(zphmass2)
           call setwarningLT
#endif   
        endif
        
      call Di(k1,k2,k3,m1q,m2q,m3q,m4q,D0m,Dim)
      D1m=Dim(1)
      D2m=Dim(2)
      D3m=Dim(3)
!!!!!!!!!!!!!!!!!!!!!Vettoriali 3 punti
      call Ci(k1,k2,m1q,m2q,m3q,C0ABC,C1ABC,C2ABC)
      call Ci(k2,k3,m1q,m3q,m4q,C0ACD,C1ACD,C2ACD)
      call Ci(k1,k3,m1q,m2q,m4q,C0ABD,C1ABD,C2ABD)
      call Ci(k2-k1,k3-k1,m2q,m3q,m4q,C0BCD,C1BCD,C2BCD)

      bc0 = C0BCD + D0m*m1q

      bc1 = (C0BCD*k1k1 + C1ACD*k1k2 - C1BCD*(-k1k1 + k1k2) + C2ACD*k1k3
     >     - C2BCD*(-k1k1 + k1k3) + (D1m*k1k1 + D2m*k1k2 + D3m*k1k3)*
     >     (-k1k1 - m1q + m2q))/2.q0
      
      bc2 = (C0BCD*k1k1 + C1ABD*k1k1 - C1BCD*(-k1k1 + k1k2) + C2ABD*k1k3
     >     - C2BCD*(-k1k1 + k1k3) + (D1m*k1k1 + D2m*k1k2 + D3m*k1k3)*
     >     (-k2k2 - m1q + m3q))/2.q0

!      print*,'B2',bc2
      bc2 = (C0BCD*k1k2 + C1ACD*k2k2 - C1BCD*(-k1k2 + k2k2) + C2ACD*k2k3
     >     - C2BCD*(-k1k2 + k2k3) + (D1m*k1k2 + D2m*k2k2 + D3m*k2k3)*
     >     (-k1k1 - m1q + m2q))/2.q0 +bc2

      bc2 = bc2/2.q0
      
      
      bc3 = (C0BCD*k1k1 + C1ABC*k1k1 + C2ABC*k1k2 - C1BCD*(-k1k1 + k1k2)
     >     - C2BCD*(-k1k1 + k1k3) + (D1m*k1k1 + D2m*k1k2 + D3m*k1k3)
     >     *(-k3k3 - m1q + m4q))/2.q0

      bc3 = (C0BCD*k1k3 + C1ACD*k2k3 - C1BCD*(-k1k3 + k2k3) + C2ACD*k3k3
     >     - C2BCD*(-k1k3 + k3k3) + (D1m*k1k3 + D2m*k2k3 + D3m*k3k3)*
     >     (-k1k1 - m1q + m2q))/2.q0 + bc3

      bc3 = bc3/2.q0
      
      bc4 = (C0BCD*k1k2 + C1ABD*k1k2 - C1BCD*(-k1k2 + k2k2) + C2ABD*k2k3
     >     - C2BCD*(-k1k2 + k2k3) + (D1m*k1k2 + D2m*k2k2 + D3m*k2k3)*
     >     (-k2k2 - m1q + m3q))/2.q0

      
      bc5 = (C0BCD*k1k3 + C1ABD*k1k3 - C1BCD*(-k1k3 + k2k3) + C2ABD*k3k3
     >     - C2BCD*(-k1k3 + k3k3) + (D1m*k1k3 + D2m*k2k3 + D3m*k3k3)*
     >     (-k2k2 - m1q + m3q))/2.q0

      
      bc5 = (C0BCD*k1k2 + C1ABC*k1k2 + C2ABC*k2k2 - C1BCD*(-k1k2 + k2k2)
     >     - C2BCD*(-k1k2 + k2k3) + (D1m*k1k2 + D2m*k2k2 + D3m*k2k3)*
     >     (-k3k3 - m1q + m4q))/2.q0 +bc5

      bc5 =bc5/2.q0
      
      bc6 = (C0BCD*k1k3 + C1ABC*k1k3 + C2ABC*k2k3 - C1BCD*(-k1k3 + k2k3)
     >     - C2BCD*(-k1k3 + k3k3) + (D1m*k1k3 + D2m*k2k3 + D3m*k3k3)*
     >     (-k3k3 - m1q + m4q))/2.q0

      
      Dijm(1) = (bc6*k1k1*k2k2 + bc0*k1k3**2*k2k2 + 
     &    bc0*k1k1*k2k3**2 - bc6*k1k2**2*1.q0 - 
     &    bc4*k1k3**2*1.q0 - bc1*k2k3**2*1.q0 + 
     &    bc5*k1k2*k1k3*2.q0 - bc3*k1k3*k2k2*2.q0 - 
     &    bc5*k1k1*k2k3*2.q0 + bc3*k1k2*k2k3*2.q0 + 
     &    bc2*k1k3*k2k3*2.q0 - bc0*k1k2*k1k3*k2k3*2.q0 + 
     &    k3k3*(bc4*k1k1 + bc0*k1k2**2 + bc1*k2k2 - 
     &       bc0*k1k1*k2k2*1.q0 - bc2*k1k2*2.q0))/
     &  (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &    k1k1*(k2k3**2 - k2k2*k3k3*1.q0) - k1k2*k1k3*k2k3*2.q0)
        Dijm(2) = (bc6*(k1k3**2*k2k2**2 + 
     &       k2k2*k3k3*(k1k1*k2k2 - k1k2**2*1.q0) - 
     &       k1k2*k1k3*k2k2*k2k3*2.q0 + 
     &       k2k3**2*(-(k1k1*k2k2*1.q0) + k1k2**2*2.q0)) + 
     &    k3k3**2*(bc4*(k1k2**2 + k1k1*k2k2) + 
     &       k2k2*(bc0*k1k2**2 - bc0*k1k1*k2k2*1.q0 + 
     &          bc1*k2k2*2.q0 - bc2*k1k2*4.q0)) + 
     &    k2k3*(-(bc5*k1k3**2*k2k2*2.q0) + 
     &       k2k3**3*(-(bc0*k1k1*1.q0) + bc1*2.q0) + 
     &       k2k3**2*2.q0*
     &        (bc5*k1k1 + bc0*k1k2*k1k3 - bc3*k1k2*2.q0 - 
     &          bc2*k1k3*2.q0) + 
     &       k1k3*k2k3*(-(bc0*k1k3*k2k2*1.q0) + bc4*k1k3*2.q0 + 
     &          bc3*k2k2*4.q0)) - 
     &    k3k3*1.q0*(k2k3*2.q0*
     &        (k1k2*(bc5*k1k2 + bc4*k1k3) + 
     &          k2k2*(bc5*k1k1 + bc0*k1k2*k1k3 - 
     &             bc3*k1k2*2.q0 - bc2*k1k3*2.q0)) + 
     &       k2k3**2*(bc4*k1k1 + bc0*k1k2**2 - 
     &          bc0*k1k1*k2k2*2.q0 - bc2*k1k2*4.q0 + 
     &          bc1*k2k2*4.q0) + 
     &       k1k3*k2k2*(bc4*k1k3 - bc0*k1k3*k2k2*1.q0 - 
     &          bc5*k1k2*4.q0 + bc3*k2k2*4.q0)))/
     &  (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &     k1k1*(k2k3**2 - k2k2*k3k3*1.q0) - k1k2*k1k3*k2k3*2.q0
     &     )**2
        Dijm(5) = (k1k3**3*k2k3*2.q0*(bc0*k1k2 - bc2*2.q0) - 
     &    k1k3*k2k3**2*1.q0*
     &     (bc0*k1k1*k1k3 + bc3*k1k1*2.q0 - bc1*k1k3*2.q0) + 
     &    k1k3**3*(-(bc0*k1k3*k2k2*1.q0) + bc4*k1k3*2.q0 + 
     &       bc3*k2k2*2.q0) + 
     &    bc6*(k1k1*(k1k1*(k2k3**2 + k2k2*k3k3) - 
     &          k1k3**2*k2k2*1.q0) - k1k1*k1k2*k1k3*k2k3*2.q0 + 
     &       k1k2**2*(-(k1k1*k3k3*1.q0) + k1k3**2*2.q0)) - 
     &    bc5*(k1k2*k1k3 - k1k1*k2k3*1.q0)*
     &     (k1k3**2 - k1k1*k3k3*1.q0)*4.q0 + 
     &    k3k3**2*((bc1 + bc0*k1k1)*k1k2**2 + 
     &       k1k1*k2k2*(bc1 - bc0*k1k1*1.q0) + 
     &       bc4*k1k1**2*2.q0 - bc2*k1k1*k1k2*4.q0) - 
     &    k3k3*1.q0*(k1k1*k2k3**2*(bc1 - bc0*k1k1*1.q0) + 
     &       bc3*k1k3*(k1k2**2 + k1k1*k2k2)*2.q0 + 
     &       k1k3*k2k3*2.q0*
     &        ((bc1 + bc0*k1k1)*k1k2 - bc2*k1k1*2.q0) - 
     &       bc3*k1k1*k1k2*k2k3*4.q0 + 
     &       k1k3**2*(bc0*k1k2**2 + bc1*k2k2 - 
     &          bc0*k1k1*k2k2*2.q0 + bc4*k1k1*4.q0 - 
     &          bc2*k1k2*4.q0)))/
     &  (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &     k1k1*(k2k3**2 - k2k2*k3k3*1.q0) - k1k2*k1k3*k2k3*2.q0
     &     )**2
        Dijm(7) = ((bc4*k1k1 + bc0*k1k2**2 + bc1*k2k2 - 
     &       bc0*k1k1*k2k2*1.q0 - bc2*k1k2*2.q0)*
     &     (k1k1*k2k2*k3k3 - k1k3**2*k2k2*1.q0 - 
     &       k1k1*k2k3**2*1.q0 - k1k2**2*k3k3*1.q0 + 
     &       k1k2*k1k3*k2k3*2.q0) + 
     &    2.q0*(bc6*(k1k2**2 - k1k1*k2k2*1.q0)**2 + 
     &       k2k3*(k1k2*(bc5*k1k1*k1k2 + bc2*k1k2*k1k3 - 
     &             bc3*k1k2**2*1.q0 - bc4*k1k1*k1k3*1.q0) + 
     &          k2k2*(bc3*k1k1*k1k2 + bc2*k1k1*k1k3 - 
     &             bc5*k1k1**2*1.q0 - bc1*k1k2*k1k3*1.q0))*2.q0
     &        + k2k3**2*(bc4*k1k1**2 + 
     &          k1k2*(bc1*k1k2 - bc2*k1k1*2.q0)) + 
     &       k1k3*(k1k2*k2k2*
     &           (bc5*k1k1 + bc3*k1k2 - bc2*k1k3*1.q0)*2.q0 + 
     &          k2k2**2*(bc1*k1k3 - bc3*k1k1*2.q0) + 
     &          k1k2**2*(bc4*k1k3 - bc5*k1k2*2.q0))))/
     &  (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &     k1k1*(k2k3**2 - k2k2*k3k3*1.q0) - k1k2*k1k3*k2k3*2.q0
     &     )**2
        Dijm(3) = (bc6*(k1k2**3*k3k3 - 
     &       k1k2*(k1k3**2*k2k2 + k1k1*k2k3**2 + 
     &          k1k1*k2k2*k3k3)*1.q0 + k1k1*k1k3*k2k2*k2k3*2.q0)
     &      + k3k3**2*(bc2*k1k1*k2k2 + bc0*k1k1*k1k2*k2k2 - 
     &       bc0*k1k2**3*1.q0 - bc4*k1k1*k1k2*2.q0 - 
     &       bc1*k1k2*k2k2*2.q0 + bc2*k1k2**2*3.q0) + 
     &    k2k3*(k2k3**2*(bc3*k1k1 + bc0*k1k1*k1k3 - 
     &          bc1*k1k3*2.q0) + 
     &       k1k3*k2k3*2.q0*
     &        (bc3*k1k2 - bc0*k1k2*k1k3*1.q0 + bc2*k1k3*2.q0) + 
     &       k1k3**2*(bc0*k1k3*k2k2 - bc4*k1k3*2.q0 - 
     &          bc3*k2k2*3.q0)) + 
     &    bc5*(k1k3**3*k2k2 + k1k2*k1k3**2*k2k3*2.q0 - 
     &       k1k3*1.q0*(k1k1*k2k2*k3k3 + k1k1*k2k3**2*3.q0 + 
     &          k1k2**2*k3k3*3.q0) + k1k1*k1k2*k2k3*k3k3*4.q0)
     &     + k3k3*(bc4*k1k3*(k1k2*k1k3 + k1k1*k2k3)*2.q0 - 
     &       k2k3**2*1.q0*
     &        (bc2*k1k1 + bc0*k1k1*k1k2 - bc1*k1k2*2.q0) - 
     &       k2k3*1.q0*(k2k2*
     &           (bc3*k1k1 + bc0*k1k1*k1k3 - bc1*k1k3*2.q0) + 
     &          k1k2*(bc3*k1k2 - bc0*k1k2*k1k3*1.q0 + 
     &             bc2*k1k3*2.q0)*3.q0) - 
     &       k1k3*k2k2*1.q0*
     &        (bc2*k1k3 + bc0*k1k2*k1k3 - bc3*k1k2*4.q0)))/
     &  (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &     k1k1*(k2k3**2 - k2k2*k3k3*1.q0) - k1k2*k1k3*k2k3*2.q0
     &     )**2
        Dijm(4) = (bc4*k1k3**3*k2k2 + bc2*k1k1*k2k3**3 + 
     &    bc0*k1k1*k1k2*k2k3**3 - bc0*k1k3**3*k2k2**2*1.q0 - 
     &    bc4*k1k1*k1k3*k2k3**2*1.q0 - 
     &    bc3*k1k1*k2k2*k2k3**2*1.q0 - 
     &    bc0*k1k1*k1k3*k2k2*k2k3**2*1.q0 + 
     &    bc5*k1k2**2*k1k3*k2k3*2.q0 + 
     &    bc2*k1k2*k1k3*k2k3**2*2.q0 - 
     &    bc0*k1k2**2*k1k3*k2k3**2*2.q0 + 
     &    bc1*k1k3*k2k2*k2k3**2*2.q0 - bc1*k1k2*k2k3**3*2.q0 - 
     &    bc6*(k1k2**2 - k1k1*k2k2*1.q0)*
     &     (k1k2*k2k3 - k1k3*k2k2*1.q0)*2.q0 - 
     &    bc5*k1k2*k1k3**2*k2k2*3.q0 + 
     &    bc3*k1k3**2*k2k2**2*3.q0 - 
     &    bc2*k1k3**2*k2k2*k2k3*3.q0 + 
     &    bc0*k1k2*k1k3**2*k2k2*k2k3*3.q0 - 
     &    bc5*k1k1*k1k2*k2k3**2*3.q0 + 
     &    bc5*k1k1*k1k3*k2k2*k2k3*4.q0 + 
     &    bc3*k1k2**2*k2k3**2*4.q0 + 
     &    k3k3*(-(bc4*k1k3*(k1k2**2 + k1k1*k2k2)*1.q0) + 
     &       bc5*(k1k2**3 - k1k1*k1k2*k2k2*1.q0) + 
     &       bc4*k1k1*k1k2*k2k3*2.q0 - 
     &       k2k3*1.q0*(bc2*k1k1*k2k2 + bc0*k1k1*k1k2*k2k2 - 
     &          bc0*k1k2**3*1.q0 - bc1*k1k2*k2k2*2.q0 + 
     &          bc2*k1k2**2*3.q0) + 
     &       k2k2*(k2k2*(bc3*k1k1 + bc0*k1k1*k1k3 - 
     &             bc1*k1k3*2.q0) - 
     &          k1k2*1.q0*
     &           (bc3*k1k2 + bc0*k1k2*k1k3 - bc2*k1k3*4.q0)))
     &     - bc3*k1k2*k1k3*k2k2*k2k3*6.q0)/
     &  (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &     k1k1*(k2k3**2 - k2k2*k3k3*1.q0) - k1k2*k1k3*k2k3*2.q0
     &     )**2
        Dijm(6) = -((1.q0*(bc1*k1k3**2*k2k2*k2k3 + 
     &        bc0*k1k1*k1k3**2*k2k2*k2k3 + 
     &        bc0*k1k1**2*k2k3**3 - bc2*k1k3**3*k2k2*1.q0 - 
     &        bc0*k1k2*k1k3**3*k2k2*1.q0 - 
     &        bc1*k1k1*k2k3**3*1.q0 + bc4*k1k2*k1k3**3*2.q0 - 
     &        bc3*k1k2**2*k1k3*k2k3*2.q0 - 
     &        bc4*k1k1*k1k3**2*k2k3*2.q0 - 
     &        bc2*k1k2*k1k3**2*k2k3*2.q0 + 
     &        bc0*k1k2**2*k1k3**2*k2k3*2.q0 + 
     &        bc6*(k1k2**2 - k1k1*k2k2*1.q0)*
     &         (k1k2*k1k3 - k1k1*k2k3*1.q0)*2.q0 + 
     &        bc3*k1k2*k1k3**2*k2k2*3.q0 + 
     &        bc3*k1k1*k1k2*k2k3**2*3.q0 + 
     &        bc2*k1k1*k1k3*k2k3**2*3.q0 - 
     &        bc0*k1k1*k1k2*k1k3*k2k3**2*3.q0 - 
     &        bc3*k1k1*k1k3*k2k2*k2k3*4.q0 + 
     &        k3k3*(bc3*(k1k1*k1k2*k2k2 - k1k2**3*1.q0) + 
     &           k1k3*(bc2*k1k1*k2k2 + bc0*k1k1*k1k2*k2k2 - 
     &              bc0*k1k2**3*1.q0 - bc4*k1k1*k1k2*2.q0 - 
     &              bc1*k1k2*k2k2*2.q0 + bc2*k1k2**2*3.q0) + 
     &           k2k3*((bc1 + bc0*k1k1)*k1k2**2 + 
     &              k1k1*k2k2*(bc1 - bc0*k1k1*1.q0) + 
     &              bc4*k1k1**2*2.q0 - bc2*k1k1*k1k2*4.q0)) + 
     &        bc5*(k1k1*(k1k3**2*k2k2 - 
     &              k1k1*1.q0*(k2k2*k3k3 + k2k3**2*3.q0)) + 
     &           k1k2**2*(k1k1*k3k3 - k1k3**2*4.q0) + 
     &           k1k1*k1k2*k1k3*k2k3*6.q0)))/
     &    (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &       k1k1*(k2k3**2 - k2k2*k3k3*1.q0) - 
     &       k1k2*k1k3*k2k3*2.q0)**2)
        call clearcache
      return
      end

*-----------ROUTINES TO REDUCE EI,EIJ-------------------*      

      
      subroutine Eij(k1,k2,k3,k4,m1,m2,m3,m4,m5,E0m,Eim,Eijm)
c     In this routine we evaluate the numerical decomposition of E_ij
C     E_i and E_0 as a function of D0,C0,B0,A0. We use the package
c     hexagon.m of https://arxiv.org/pdf/0812.2134 to get the symbolic
c     decomposition. The convention for the loop integral is
c      
c          1             1      {1,q^mu,q^muq^nu}       1              1
c     --------     ----------    -------------   -------------  --------------
c     (q^2-m1q)  ((q+k1)^2-m2q)  ((q+k2)^2-m3q)  ((q+k3)^2-m4q)  ((q+k4)^2-m5q) 
c
!      implicit none
      implicit complex*32 (u-z)
#ifdef LOOPTOOLS
#include "looptools.h"
#endif
      double precision k1(0:3),k2(0:3),k3(0:3),k4(0:3)
      double precision p1LT(0:3),p2LT(0:3),p3LT(0:3),p4LT(0:3),p5LT(0:3)
      complex*32 m1q,m2q,m3q,m4q,m5q
      complex*16 m1,m2,m3,m4,m5
      complex*32 E0m,E1m,E2m,E3m,E4m
      complex*32 E11,E12,E13,E14
      complex*32 E22,E23,E24
      complex*32 E33,E34
      complex*32 E44
      complex*32 den,den1,den2,den3,den4,den5
      complex*32 Eim(4),Eijm(11)
      complex*32 D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q
      complex*32 D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q
      complex*32 D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q
      complex*32 D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q
      complex*32 D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q

      complex*32 C0ap1LT2p2LT2x12LTm1qm2qm3q
      complex*32 C0ap1LT2x23LTx45LTm1qm2qm4q
      complex*32 C0ap1LT2x51LTp5LT2m1qm2qm5q
      complex*32 C0ax12LTp3LT2x45LTm1qm3qm4q
      complex*32 C0ax12LTx34LTp5LT2m1qm3qm5q
      complex*32 C0ax23LTp2LT2p3LT2m4qm2qm3q
      complex*32 C0ax34LTp3LT2p4LT2m5qm3qm4q
      complex*32 C0ax45LTp4LT2p5LT2m1qm4qm5q
      complex*32 C0ax51LTp2LT2x34LTm5qm2qm3q
      complex*32 C0ax51LTx23LTp4LT2m5qm2qm4q      

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

      
      common/lt/ltint,quadevents,ilt
      integer ltint,quadevents,ilt
      
      real*8 deltauv,lambda2,mu2dim,m2min
      common/loop/deltauv,lambda2,mu2dim,m2min

      external dotq
      real*16 dotq

      external closetomassc
      complex*32 closetomassc
      
      common/softphotonmass/phmass
      real*8 phmass
      
      real*16 zphmass2,zdeltauv
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
      m1q   = closetomassc(m1*re)
      m2q   = closetomassc(m2*re)
      m3q   = closetomassc(m3*re)
      m4q   = closetomassc(m4*re)
      m5q   = closetomassc(m5*re)

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

!!!!!!!!!!!!!!!!!!!!LT non sa nulla di kikj                                                                               
      k1k2 = -(p2LT2-k1k1-k2k2)/2.q0
      k1k3 = -(x23LT-k1k1-k3k3)/2.q0
      k1k4 = -(x51LT-k1k1-k4k4)/2.q0
      k2k3 = -(p3LT2-k3k3-k2k2)/2.q0
      k2k4 = -(x34LT-k2k2-k4k4)/2.q0
      k3k4 = -(p4LT2-k3k3-k4k4)/2.q0

        if(ilt.eq.0) then
           ilt=1
           zdeltauv=deltauv
           zphmass2=phmass*phmass
#ifdef LOOPTOOLS
           call ltini
           call setdelta(zdeltauv)
           call setmudim(1.q2)
           call setlambda(zphmass2)
           call setwarningLT
#endif   
        endif

        
      include 'pipirad/decomposition5pt/assignements.f'

#ifdef LOOPTOOLS

      D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q = D0C(
     >     x51LT,p2LT2,p3LT2,p4LT2,x34LT,x23LT,
     >     m5q,m2q,m3q,m4q)
      D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q = D0C(
     >     x12LT,p3LT2,p4LT2,p5LT2,x45LT,x34LT,
     >     m1q,m3q,m4q,m5q)
      D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q = D0C(
     >     p1LT2,x23LT,p4LT2,p5LT2,x45LT,x51LT,
     >     m1q,m2q,m4q,m5q)
      D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q = D0C(
     >     p1LT2,p2LT2,x34LT,p5LT2,x12LT,x51LT,
     >     m1q,m2q,m3q,m5q)
      D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q = D0C(
     >     p1LT2,p2LT2,p3LT2,x45LT,x12LT,x23LT,
     >     m1q,m2q,m3q,m4q)

      
      C0ap1LT2p2LT2x12LTm1qm2qm3q = C0C(
     >     p1LT2,p2LT2,x12LT,m1q,m2q,m3q)
      C0ap1LT2x23LTx45LTm1qm2qm4q = C0C(
     >     p1LT2,x23LT,x45LT,m1q,m2q,m4q)
      C0ap1LT2x51LTp5LT2m1qm2qm5q = C0C(
     >     p1LT2,x51LT,p5LT2,m1q,m2q,m5q)
      C0ax12LTp3LT2x45LTm1qm3qm4q = C0C(
     >     x12LT,p3LT2,x45LT,m1q,m3q,m4q)
      C0ax12LTx34LTp5LT2m1qm3qm5q = C0C(
     >     x12LT,x34LT,p5LT2,m1q,m3q,m5q)
      C0ax23LTp2LT2p3LT2m4qm2qm3q = C0C(
     >     x23LT,p2LT2,p3LT2,m4q,m2q,m3q)
      C0ax34LTp3LT2p4LT2m5qm3qm4q = C0C(
     >     x34LT,p3LT2,p4LT2,m5q,m3q,m4q)
      
      C0ax45LTp4LT2p5LT2m1qm4qm5q = C0C(
     >     x45LT,p4LT2,p5LT2,m1q,m4q,m5q)

*------------The routing below gives the correct description when M_pipi \simeq 2*m_\pi and elsewhere gives the right result
      if(m5q.eq.(0.q0,0.q0)) then
         C0ax45LTp4LT2p5LT2m1qm4qm5q = C0C(
     >        p5LT2,p4LT2,x45LT,m5q,m1q,m4q)
      endif
      
      C0ax51LTp2LT2x34LTm5qm2qm3q = C0C(
     >     x51LT,p2LT2,x34LT,m5q,m2q,m3q)
      C0ax51LTx23LTp4LT2m5qm2qm4q = C0C(
     >     x51LT,x23LT,p4LT2,m5q,m2q,m4q)
#endif
!!!!!!!E0!!!!!!!!

      include 'pipirad/decomposition5pt/E0.F'

!!!!!!!Ei!!!!!!!!

      include 'pipirad/decomposition5pt/Ei.F'
      
!!!!!!Eij!!!!!!!!
      
      include 'pipirad/decomposition5pt/Eij.F'
      call clearcache
      return
      
      end
            
      subroutine Eij_dispq(k1,k2,k3,k4,m1,m2,m3,m4,m5,
     >     s,disp,E0m,Eim,Eijm)
      implicit complex*32 (u-z)
!      implicit none
#ifdef LOOPTOOLS
#include "looptools.h"
#endif
      double precision k1(0:3),k2(0:3),k3(0:3),k4(0:3)
      double precision p1LT(0:3),p2LT(0:3),p3LT(0:3),p4LT(0:3),p5LT(0:3)
      double precision pi(0:3),pf(0:3)
      complex*32 m1q,m2q,m3q,m4q,m5q
      complex*16 m1,m2,m3,m4,m5
      complex*32 E0m,E1m,E2m,E3m,E4m
      complex*32 E11,E12,E13,E14
      complex*32 E22,E23,E24
      complex*32 E33,E34
      complex*32 E44
      complex*32 den,den1,den2,den3,den4,den5
      complex*32 Eim(4),Eijm(11)
      complex*32 D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q
      complex*32 D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q
      complex*32 D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q
      complex*32 D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q
      complex*32 D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q

      complex*32 C0ap1LT2p2LT2x12LTm1qm2qm3q
      complex*32 C0ap1LT2x23LTx45LTm1qm2qm4q
      complex*32 C0ap1LT2x51LTp5LT2m1qm2qm5q
      complex*32 C0ax12LTp3LT2x45LTm1qm3qm4q
      complex*32 C0ax12LTx34LTp5LT2m1qm3qm5q
      complex*32 C0ax23LTp2LT2p3LT2m4qm2qm3q
      complex*32 C0ax34LTp3LT2p4LT2m5qm3qm4q
      complex*32 C0ax45LTp4LT2p5LT2m1qm4qm5q
      complex*32 C0ax51LTp2LT2x34LTm5qm2qm3q
      complex*32 C0ax51LTx23LTp4LT2m5qm2qm4q      

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

      
      common/lt/ltint,quadevents,ilt
      integer ltint,quadevents,ilt
      
      real*8 deltauv,lambda2,mu2dim,m2min
      common/loop/deltauv,lambda2,mu2dim,m2min

      external dotq
      real*16 dotq

      external closetomassc
      complex*32 closetomassc
      
      common/softphotonmass/phmass
      real*8 phmass
      
      real*16 zphmass2,zdeltauv

      
*------------FsQED
      common/hparam/ampi
      double precision ampi,s,sp	

      integer ifound
      
      external Itot1ovsp_dimreg
      real*8 Itot1ovsp_dimreg

      external Itotspovsp_dimreg
      real*8 Itotspovsp_dimreg

      
      external Itotspovsp_lambdarad
      real*8 Itotspovsp_lambdarad

      external Itot1ovsp_lambdarad
      real*8 Itot1ovsp_lambdarad

      common/momentadisp/momenta
      double precision momenta(4,0:3)

      external ImD0
      double complex ImD0

      logical cond
      integer disp
      Complex*32 setrezq
      real*16 masses(5)
      
            
      D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q = (0.q0, 0.q0)
      D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q = (0.q0, 0.q0)
      D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q = (0.q0, 0.q0)
      D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q = (0.q0, 0.q0)
      D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q = (0.q0, 0.q0)

      C0ap1LT2p2LT2x12LTm1qm2qm3q = (0.q0, 0.q0)
      C0ap1LT2x23LTx45LTm1qm2qm4q = (0.q0, 0.q0)
      C0ap1LT2x51LTp5LT2m1qm2qm5q = (0.q0, 0.q0)
      C0ax12LTp3LT2x45LTm1qm3qm4q = (0.q0, 0.q0)
      C0ax12LTx34LTp5LT2m1qm3qm5q = (0.q0, 0.q0)
      C0ax23LTp2LT2p3LT2m4qm2qm3q = (0.q0, 0.q0)
      C0ax34LTp3LT2p4LT2m5qm3qm4q = (0.q0, 0.q0)
      C0ax45LTp4LT2p5LT2m1qm4qm5q = (0.q0, 0.q0)
      C0ax51LTp2LT2x34LTm5qm2qm3q = (0.q0, 0.q0)
      C0ax51LTx23LTp4LT2m5qm2qm4q = (0.q0, 0.q0)


      re=(1.q0,0.q0)
      im=(0.q0,1.q0)
      
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
      m1q   = closetomassc(m1*re)
      m2q   = closetomassc(m2*re)
      m3q   = closetomassc(m3*re)
      m4q   = closetomassc(m4*re)
      m5q   = closetomassc(m5*re)

      
      masses(1)  = real(m1q)
      masses(2)  = real(m2q)
      masses(3)  = real(m3q)
      masses(4)  = real(m4q)
      masses(5)  = real(m5q)

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

!!!!!!!!!!!!!!!!!!!!LT non sa nulla di kikj                                                                               
      k1k2 = -(p2LT2-k1k1-k2k2)/2.q0
      k1k3 = -(x23LT-k1k1-k3k3)/2.q0
      k1k4 = -(x51LT-k1k1-k4k4)/2.q0
      k2k3 = -(p3LT2-k3k3-k2k2)/2.q0
      k2k4 = -(x34LT-k2k2-k4k4)/2.q0
      k3k4 = -(p4LT2-k3k3-k4k4)/2.q0

      

!     Classica decomposizione di loop
      if(disp.eq.1) then
         if(ilt.eq.0) then
            ilt=1
            zdeltauv=deltauv
            zphmass2=phmass*phmass
#ifdef LOOPTOOLS
            call ltini
            call setdelta(zdeltauv)
            call setmudim(1.q2)
            call setlambda(zphmass2)
            call setwarningLT
#endif   
         endif

#ifdef LOOPTOOLS

         D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q = D0C(
     >        x51LT,p2LT2,p3LT2,p4LT2,x34LT,x23LT,
     >        m5q,m2q,m3q,m4q)
         D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q = D0C(
     >        x12LT,p3LT2,p4LT2,p5LT2,x45LT,x34LT,
     >        m1q,m3q,m4q,m5q)
         D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q = D0C(
     >        p1LT2,x23LT,p4LT2,p5LT2,x45LT,x51LT,
     >        m1q,m2q,m4q,m5q)
         D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q = D0C(
     >        p1LT2,p2LT2,x34LT,p5LT2,x12LT,x51LT,
     >        m1q,m2q,m3q,m5q)
         D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q = D0C(
     >        p1LT2,p2LT2,p3LT2,x45LT,x12LT,x23LT,
     >        m1q,m2q,m3q,m4q)

         
         C0ap1LT2p2LT2x12LTm1qm2qm3q = C0C(
     >        p1LT2,p2LT2,x12LT,m1q,m2q,m3q)
         C0ap1LT2x23LTx45LTm1qm2qm4q = C0C(
     >        p1LT2,x23LT,x45LT,m1q,m2q,m4q)
         C0ap1LT2x51LTp5LT2m1qm2qm5q = C0C(
     >        p1LT2,x51LT,p5LT2,m1q,m2q,m5q)
         C0ax12LTp3LT2x45LTm1qm3qm4q = C0C(
     >        x12LT,p3LT2,x45LT,m1q,m3q,m4q)
         C0ax12LTx34LTp5LT2m1qm3qm5q = C0C(
     >        x12LT,x34LT,p5LT2,m1q,m3q,m5q)
         C0ax23LTp2LT2p3LT2m4qm2qm3q = C0C(
     >        x23LT,p2LT2,p3LT2,m4q,m2q,m3q)
         C0ax34LTp3LT2p4LT2m5qm3qm4q = C0C(
     >        x34LT,p3LT2,p4LT2,m5q,m3q,m4q)
         
         C0ax45LTp4LT2p5LT2m1qm4qm5q = C0C(
     >        x45LT,p4LT2,p5LT2,m1q,m4q,m5q)

*------------The routing below gives the correct description when M_pipi \simeq 2*m_\pi and elsewhere gives the right result
         if(m5q.eq.(0.q0,0.q0)) then
            C0ax45LTp4LT2p5LT2m1qm4qm5q = C0C(
     >           p5LT2,p4LT2,x45LT,m5q,m1q,m4q)
         endif
         
         C0ax51LTp2LT2x34LTm5qm2qm3q = C0C(
     >        x51LT,p2LT2,x34LT,m5q,m2q,m3q)
         C0ax51LTx23LTp4LT2m5qm2qm4q = C0C(
     >        x51LT,x23LT,p4LT2,m5q,m2q,m4q)
#endif
         
      elseif(disp.eq.2) then
!     Check quali D0 sono IR e sostituisci D0 con \int ds'/s' Im D0
!     e portano polo 1/(\eta-s')
         call D0IRyon(k1-k4,k2-k4,k3-k4,m5,m2,m3,m4,cond,pi,pf)
         if(cond) then
            D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q=
     >           Itot1ovsp_lambdarad(pi,pf,s,lambda2)*im
         endif

         call D0IRyon(k2,k3,k4,m1,m3,m4,m5,cond,pi,pf)
         if(cond) then
            D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q=
     >           Itot1ovsp_lambdarad(pi,pf,s,lambda2)*im
         endif

         call D0IRyon(k1,k3,k4,m1,m2,m4,m5,cond,pi,pf)
         if(cond) then
            D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q=
     >           Itot1ovsp_lambdarad(pi,pf,s,lambda2)*im
         endif
         
         call D0IRyon(k1,k2,k4,m1,m2,m3,m5,cond,pi,pf)
         if(cond) then
            D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q=
     >           Itot1ovsp_lambdarad(pi,pf,s,lambda2)*im
         endif

         call D0IRyon(k1,k2,k3,m1,m2,m3,m4,cond,pi,pf)
         if(cond) then
            D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q=
     >           Itot1ovsp_lambdarad(pi,pf,s,lambda2)*im
         endif

!     Devo restituire il coefficiente d(\eta)*\int ds'/s' Im D0, quindi devo cambiare il valore della massa dispersiva a \eta (qui nella routine s)
         do i=1,5
            if (masses(i).ge.4.d0*ampi**2) then
               ifound = i
               if (i.eq.1) m1q = s*re
               if (i.eq.2) m2q = s*re
               if (i.eq.3) m3q = s*re
               if (i.eq.4) m4q = s*re
               if (i.eq.5) m5q = s*re
            endif
         enddo
      elseif(disp.eq.3) then
!     Restituisce la parte Im di tutto ad eccezione delle D0_IR e polari e delle D0_IR (polari) restituisce
!     la parte IM sottratta della parte immaginaria che origina dal taglio che diverge a s'=\eta

         if(ilt.eq.0) then
            ilt=1
            zdeltauv=deltauv
            zphmass2=phmass*phmass
#ifdef LOOPTOOLS
            call ltini
            call setdelta(zdeltauv)
            call setmudim(1.q2)
            call setlambda(zphmass2)
            call setwarningLT
#endif   
         endif

         do i=1,5
            if (masses(i).ge.4.d0*ampi**2) then
               ifound = i
               sp = masses(i)
            endif
         enddo

#ifdef LOOPTOOLS

         D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q = D0C(
     >        x51LT,p2LT2,p3LT2,p4LT2,x34LT,x23LT,
     >        m5q,m2q,m3q,m4q)
         D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q = D0C(
     >        x12LT,p3LT2,p4LT2,p5LT2,x45LT,x34LT,
     >        m1q,m3q,m4q,m5q)
         D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q = D0C(
     >        p1LT2,x23LT,p4LT2,p5LT2,x45LT,x51LT,
     >        m1q,m2q,m4q,m5q)
         D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q = D0C(
     >        p1LT2,p2LT2,x34LT,p5LT2,x12LT,x51LT,
     >        m1q,m2q,m3q,m5q)
         D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q = D0C(
     >        p1LT2,p2LT2,p3LT2,x45LT,x12LT,x23LT,
     >        m1q,m2q,m3q,m4q)

         
         C0ap1LT2p2LT2x12LTm1qm2qm3q = C0C(
     >        p1LT2,p2LT2,x12LT,m1q,m2q,m3q)
         C0ap1LT2x23LTx45LTm1qm2qm4q = C0C(
     >        p1LT2,x23LT,x45LT,m1q,m2q,m4q)
         C0ap1LT2x51LTp5LT2m1qm2qm5q = C0C(
     >        p1LT2,x51LT,p5LT2,m1q,m2q,m5q)
         C0ax12LTp3LT2x45LTm1qm3qm4q = C0C(
     >        x12LT,p3LT2,x45LT,m1q,m3q,m4q)
         C0ax12LTx34LTp5LT2m1qm3qm5q = C0C(
     >        x12LT,x34LT,p5LT2,m1q,m3q,m5q)
         C0ax23LTp2LT2p3LT2m4qm2qm3q = C0C(
     >        x23LT,p2LT2,p3LT2,m4q,m2q,m3q)
         C0ax34LTp3LT2p4LT2m5qm3qm4q = C0C(
     >        x34LT,p3LT2,p4LT2,m5q,m3q,m4q)
         
         C0ax45LTp4LT2p5LT2m1qm4qm5q = C0C(
     >        x45LT,p4LT2,p5LT2,m1q,m4q,m5q)

*------------The routing below gives the correct description when M_pipi \simeq 2*m_\pi and elsewhere gives the right result
         if(m5q.eq.(0.q0,0.q0)) then
            C0ax45LTp4LT2p5LT2m1qm4qm5q = C0C(
     >           p5LT2,p4LT2,x45LT,m5q,m1q,m4q)
         endif
         
         C0ax51LTp2LT2x34LTm5qm2qm3q = C0C(
     >        x51LT,p2LT2,x34LT,m5q,m2q,m3q)
         C0ax51LTx23LTp4LT2m5qm2qm4q = C0C(
     >        x51LT,x23LT,p4LT2,m5q,m2q,m4q)
#endif
!     Just the imaginary part
         
         D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q = 
     >        setrezq(D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q)

         D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q = 
     >        setrezq(D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q)

         D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q = 
     >        setrezq(D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q)

         D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q = 
     >        setrezq(D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q)

         D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q = 
     >        setrezq(D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q)

         C0ap1LT2p2LT2x12LTm1qm2qm3q = 
     >        setrezq(C0ap1LT2p2LT2x12LTm1qm2qm3q)

         C0ap1LT2x23LTx45LTm1qm2qm4q = 
     >        setrezq(C0ap1LT2x23LTx45LTm1qm2qm4q)

         C0ap1LT2x51LTp5LT2m1qm2qm5q = 
     >        setrezq(C0ap1LT2x51LTp5LT2m1qm2qm5q)

         C0ax12LTp3LT2x45LTm1qm3qm4q = 
     >        setrezq(C0ax12LTp3LT2x45LTm1qm3qm4q)

         C0ax12LTx34LTp5LT2m1qm3qm5q = 
     >        setrezq(C0ax12LTx34LTp5LT2m1qm3qm5q)

         C0ax23LTp2LT2p3LT2m4qm2qm3q = 
     >        setrezq(C0ax23LTp2LT2p3LT2m4qm2qm3q)

         C0ax34LTp3LT2p4LT2m5qm3qm4q = 
     >        setrezq(C0ax34LTp3LT2p4LT2m5qm3qm4q)

         C0ax45LTp4LT2p5LT2m1qm4qm5q = 
     >        setrezq(C0ax45LTp4LT2p5LT2m1qm4qm5q)

         C0ax51LTp2LT2x34LTm5qm2qm3q = 
     >        setrezq(C0ax51LTp2LT2x34LTm5qm2qm3q)

         C0ax51LTx23LTp4LT2m5qm2qm4q = 
     >        setrezq(C0ax51LTx23LTp4LT2m5qm2qm4q)

         call D0IRyon(k1-k4,k2-k4,k3-k4,m5,m2,m3,m4,cond,pi,pf)
         if(cond) then
            D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q=
     >           D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q
     >           -ImD0(pi,pf,s,sp,lambda2)*re
         endif

         call D0IRyon(k2,k3,k4,m1,m3,m4,m5,cond,pi,pf)
         if(cond) then
            D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q=
     >           D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q
     >           -ImD0(pi,pf,s,sp,lambda2)*re
         endif

         call D0IRyon(k1,k3,k4,m1,m2,m4,m5,cond,pi,pf)
         if(cond) then
            D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q=
     >           D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q
     >           -ImD0(pi,pf,s,sp,lambda2)*re
         endif
         
         call D0IRyon(k1,k2,k4,m1,m2,m3,m5,cond,pi,pf)
         if(cond) then
            D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q=
     >           D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q
     >           -ImD0(pi,pf,s,sp,lambda2)*re
         endif

         call D0IRyon(k1,k2,k3,m1,m2,m3,m4,cond,pi,pf)
         if(cond) then
            D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q=
     >           D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q
     >           -ImD0(pi,pf,s,sp,lambda2)*re
         endif
         
         do i=1,5
            if (masses(i).ge.4.d0*ampi**2) then
               ifound = i
               if (i.eq.1) m1q = s*re
               if (i.eq.2) m2q = s*re
               if (i.eq.3) m3q = s*re
               if (i.eq.4) m4q = s*re
               if (i.eq.5) m5q = s*re
            endif
         enddo
      elseif(disp.eq.4) then

#ifdef LOOPTOOLS

         D0ax51LTp2LT2p3LT2p4LT2x34LTx23LTm5qm2qm3qm4q = D0C(
     >        x51LT,p2LT2,p3LT2,p4LT2,x34LT,x23LT,
     >        m5q,m2q,m3q,m4q)
         D0ax12LTp3LT2p4LT2p5LT2x45LTx34LTm1qm3qm4qm5q = D0C(
     >        x12LT,p3LT2,p4LT2,p5LT2,x45LT,x34LT,
     >        m1q,m3q,m4q,m5q)
         D0ap1LT2x23LTp4LT2p5LT2x45LTx51LTm1qm2qm4qm5q = D0C(
     >        p1LT2,x23LT,p4LT2,p5LT2,x45LT,x51LT,
     >        m1q,m2q,m4q,m5q)
         D0ap1LT2p2LT2x34LTp5LT2x12LTx51LTm1qm2qm3qm5q = D0C(
     >        p1LT2,p2LT2,x34LT,p5LT2,x12LT,x51LT,
     >        m1q,m2q,m3q,m5q)
         D0ap1LT2p2LT2p3LT2x45LTx12LTx23LTm1qm2qm3qm4q = D0C(
     >        p1LT2,p2LT2,p3LT2,x45LT,x12LT,x23LT,
     >        m1q,m2q,m3q,m4q)

         
         C0ap1LT2p2LT2x12LTm1qm2qm3q = C0C(
     >        p1LT2,p2LT2,x12LT,m1q,m2q,m3q)
         C0ap1LT2x23LTx45LTm1qm2qm4q = C0C(
     >        p1LT2,x23LT,x45LT,m1q,m2q,m4q)
         C0ap1LT2x51LTp5LT2m1qm2qm5q = C0C(
     >        p1LT2,x51LT,p5LT2,m1q,m2q,m5q)
         C0ax12LTp3LT2x45LTm1qm3qm4q = C0C(
     >        x12LT,p3LT2,x45LT,m1q,m3q,m4q)
         C0ax12LTx34LTp5LT2m1qm3qm5q = C0C(
     >        x12LT,x34LT,p5LT2,m1q,m3q,m5q)
         C0ax23LTp2LT2p3LT2m4qm2qm3q = C0C(
     >        x23LT,p2LT2,p3LT2,m4q,m2q,m3q)
         C0ax34LTp3LT2p4LT2m5qm3qm4q = C0C(
     >        x34LT,p3LT2,p4LT2,m5q,m3q,m4q)
         
         C0ax45LTp4LT2p5LT2m1qm4qm5q = C0C(
     >        x45LT,p4LT2,p5LT2,m1q,m4q,m5q)

*------------The routing below gives the correct description when M_pipi \simeq 2*m_\pi and elsewhere gives the right result
         if(m5q.eq.(0.q0,0.q0)) then
            C0ax45LTp4LT2p5LT2m1qm4qm5q = C0C(
     >           p5LT2,p4LT2,x45LT,m5q,m1q,m4q)
         endif
         
         C0ax51LTp2LT2x34LTm5qm2qm3q = C0C(
     >        x51LT,p2LT2,x34LT,m5q,m2q,m3q)
         C0ax51LTx23LTp4LT2m5qm2qm4q = C0C(
     >        x51LT,x23LT,p4LT2,m5q,m2q,m4q)
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
      
      m1  =  masses(1)*(1.d0,0.d0)
      m2  =  masses(2)*(1.d0,0.d0)
      m3  =  masses(3)*(1.d0,0.d0)
      m4  =  masses(4)*(1.d0,0.d0)
      m5  =  masses(5)*(1.d0,0.d0)
      call clearcache
      return
       end
