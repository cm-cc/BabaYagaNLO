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
	double precision virtual,virtual1
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
	double complex ffpiis
	double precision ffpiis2
	double complex ffpifs
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
********************

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

******ALPHA RUNNING
	ffpiis=ffpiis*vpolc(s1)
	ffpiis2=ffpiis*CONJG(ffpiis)

	ffpifs=ffpifs*vpolc(s)
	ffpifs2=ffpifs*CONJG(ffpifs)
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

        virtual=0.d0
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
     	elmati=real(elmat1)*ffpiis2
 	elmatf=REAL(elmat2*CONJG(ffpifs)*ffpiis)
        if(i.eq.2) then
           elmati=-elmati
        endif
!        print*,'isr pent',elmati,elmatf
 	virtual=virtual+real(elmati+elmatf)
        virtual=4d0*virtual*e6*alpha/(2d0*pi)

#ifdef LOOPTOOLS
        call clearcache
#endif
       end
