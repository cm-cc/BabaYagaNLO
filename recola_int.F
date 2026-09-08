      subroutine recola_init
#ifdef RECOLA
      use recola
      use globals_rcl, only : set_EWloop_weak_rcl,set_EWloop_QED_rcl
      implicit none
      character*100  path
      common/cmnpath/path

      integer icll
      common/colliermode/icll
      common/scaling/scaling

      complex*16 scaling

      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
      double precision phmass
      common/softphotonmass/phmass
      double precision deltauv,mu2dim,mudim,lambda2,m2min,pi216m1,ec6
      common/duvmudim/deltauv,mudim,mu2dim,m2min,pi216m1,ec6
      common/gammamass/lambda2

      double precision eps
      common/epssoft/eps
      double precision ebeam
      common/energiabeam/ebeam
      double precision egmax
      common/egmaxforrecola/egmax
      integer irecola
      common/recolamode/irecola

ccc done in svfactor at the beginning      phmass = 1d-10
      egmax =  eps*ebeam

      if (irecola.eq.0) return
!#default
      deltauv=0d0

      lambda2= phmass*phmass
      m2min= 0.d0

      pi216m1 = pi*pi*16.d0
      pi216m1 = 1.d0/pi216m1

      ec6 = 4.d0*pi*alpha
      ec6 = ec6*ec6*ec6

! the line below, with .false., the amplitudes are 0 with RECOLA2.... (??)
      call set_momenta_correction_rcl(.true.)

c      call set_on_shell_scheme_rcl
      call set_complex_mass_scheme_rcl

      call set_pole_mass_electron_rcl(ame)
      call unset_light_electron_rcl

      call set_pole_mass_muon_rcl(ammu,0.d0)
      call unset_light_muon_rcl

c$$$      call set_pole_mass_Z_rcl(775.26d-3,147.4d-3)
c$$$      call set_pole_mass_Z_rcl(3.5d0,1.0d0)
c$$$      call set_pole_mass_Z_rcl(0.8d0,2.5d-4)

      call set_pole_mass_tau_rcl(1.77686d0,0.d0) ! from amasses in sv.F
      call unset_light_tau_rcl

      call set_pole_mass_top_rcl(172.5d0,0.d0) ! from amasses in sv.F
      call unset_light_top_rcl
c test per Delta independence

      call set_light_fermions_rcl(ame/100d0)


* switching off couplings...

      call switchoff_coupling3_rcl('e-','e+','Z')
      call switchoff_coupling3_rcl('mu-','mu+','Z')
      call switchoff_coupling3_rcl('e-','e+','H')
      call switchoff_coupling3_rcl('mu-','mu+','H')
      call switchoff_coupling3_rcl('p0','p0','H')
      call switchoff_coupling3_rcl('p+','p-','H')
      call switchoff_coupling3_rcl('p0','p0','Z')
      call switchoff_coupling3_rcl('p+','p-','Z')
      call switchoff_coupling3_rcl('p0','p0','A')
      call switchoff_coupling3_rcl('p+','p-','A')
      call switchoff_coupling3_rcl('W+','W-','A')
      call switchoff_coupling3_rcl('e-','e+','H')
      call switchoff_coupling3_rcl('e-','e+','p0')
      call switchoff_coupling3_rcl('mu-','mu+','p0')
      call switchoff_coupling3_rcl('e-','nu_e~','W+')
      call switchoff_coupling3_rcl('e-','nu_e~','W-')
      call switchoff_coupling3_rcl('e-','nu_e~','p+')
      call switchoff_coupling3_rcl('e-','nu_e~','p-')
      call switchoff_coupling3_rcl('e+','nu_e','W+')
      call switchoff_coupling3_rcl('e+','nu_e','W-')
      call switchoff_coupling3_rcl('e+','nu_e','p+')
      call switchoff_coupling3_rcl('e+','nu_e','p-')
      call switchoff_coupling3_rcl('mu-','nu_mu~','W+')
      call switchoff_coupling3_rcl('mu-','nu_mu~','W-')
      call switchoff_coupling3_rcl('mu-','nu_mu~','p+')
      call switchoff_coupling3_rcl('mu-','nu_mu~','p-')
      call switchoff_coupling3_rcl('mu+','nu_mu','W+')
      call switchoff_coupling3_rcl('mu+','nu_mu','W-')
      call switchoff_coupling3_rcl('mu+','nu_mu','p+')
      call switchoff_coupling3_rcl('mu+','nu_mu','p-')

      call use_alpha0_scheme_rcl(alpha)
c      call use_alphaz_scheme_rcl(alpha)

      call use_mass_reg_soft_rcl(phmass)
      call set_delta_uv_rcl(deltauv)
!      call set_print_level_squared_amplitude_rcl(2)
** no...      call switchoff_coupling3_rcl('e-','e+','Z')

      if (irecola.eq.2) then
         call set_EWloop_QED_rcl
!         call set_EWloop_WEAK_rcl
      endif
      ! to run QED corrections only
!     attenzione usando solo EWloop_WEAK, il soft photon va tolto,
!     ed sv deve essere consistente
!      call set_EWloop_WEAK_rcl  ! to run weak corrections only


      icll = 1                  ! collier has now been initialised
      call set_collier_mode_rcl(1)

      call set_output_file_rcl(path(1:len_trim(path))//'recola.out')
c      call set_output_file_rcl('*')
      call set_collier_output_dir_rcl(path(1:len_trim(path))//'cll')

** e+e-\to\gamma\gamma LO
      call define_process_rcl(10, 'e- e+ -> A A', 'LO')
** e+e-\to\gamma\gamma NLO
      call define_process_rcl(11, 'e- e+ -> A A', 'NLO')
** e+e-\to\gamma\gamma\gamma LO
      call define_process_rcl(12, 'e- e+ -> A A A', 'LO')
c$$$
c$$$** e+e-\to 4 \gamma LO
c$$$      call define_process_rcl(13, 'e- e+ -> A A A A', 'LO')
c$$$** e+e-\to 5 \gamma LO
c$$$      call define_process_rcl(14, 'e- e+ -> A A A A A', 'LO')


** e+e-\to e+e- LO
      call define_process_rcl(20, 'e- e+ -> e- e+', 'LO')
** e+e-\to e+e- NLO
      call define_process_rcl(21, 'e- e+ -> e- e+', 'NLO')
** e+e-\to e+e.\gamma LO
      call define_process_rcl(22, 'e- e+ -> e- e+ A', 'LO')
** e+e-\to e+e.\gamma\gamma LO
**      call set_draw_level_branches_rcl(2)
** e+e-\to mu+mu- LO
      call define_process_rcl(30, 'e- e+ -> mu- mu+', 'LO')
*     * e+e-\to mu+mu- NLO
      call define_process_rcl(31, 'e- e+ -> mu- mu+', 'NLO')

******************************************************
*      RADIATIVE PROCESSES - RECOLA DEFINITIONS      *
******************************************************
      call define_process_rcl(90, 'e- e+ -> e- e+ A'  ,   'LO')
      call define_process_rcl(80, 'e- e+ -> mu- mu+ A',   'LO')

      call define_process_rcl(91, 'e- e+ -> e- e+ A'  ,  'NLO')
      call define_process_rcl(81, 'e- e+ -> mu- mu+ A',  'NLO')


      call define_process_rcl(70, 'e- e+ -> e- e+ A A',   'LO')
      call define_process_rcl(99, 'e- e+ -> mu- mu+ A A', 'LO')
      call define_process_rcl(56, 'e- e+ -> mu- mu+ A A A', 'LO')
      call define_process_rcl(57, 'e- e+ -> mu- mu+ A A A A', 'LO')
      call generate_processes_rcl
#endif
      return
      end
***************************************************************
      subroutine exitrecola
#ifdef RECOLA
      use recola
      implicit none
      integer irecola
      common/recolamode/irecola
      if (irecola.eq.0) return
      call reset_recola_rcl
#endif
      end subroutine exitrecola
***************************************************************
      subroutine recola_ampl(ng,p1,p2,p3,p4,p5,born,oneloop,qphmat)
#ifdef RECOLA
      use recola
#endif
      implicit double precision (a-h,o-z)
      dimension p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3)
      parameter (imaxph = 40)
      dimension qphmat(imaxph,0:3)
      dimension pmat(0:3,1:8),A2(2)
      double precision egmax

      common/scaling/scaling
      complex*16 scaling

      double precision phmass
      common/softphotonmass/phmass

      common/egmaxforrecola/egmax
      character*6 ord
      common/qedORDER/ord
      character*2 fs
      common/finalstate/fs
      real*8 rsoft,risr,rfsr,rifi
      dimension p6(0:3),p7(0:3)
      common/p6for4g/p6,p7



      born = 0.d0
      oneloop = 0.d0
#ifdef RECOLA
      pmat(:,1) = p1
      pmat(:,2) = p2
      pmat(:,3) = p3
      pmat(:,4) = p4
      pmat(:,5) = p5
      pmat(:,6) = qphmat(2,:)
      pmat(:,7) = qphmat(3,:)
      pmat(:,8) = qphmat(4,:)


      IF (fs.eq.'gg') THEN ! e+e-\to\gamma\gamma
      if (ng.eq.0) then
         if (ord.eq.'born') then
            call compute_process_rcl(10,pmat(:,1:4),'LO',A2)
            born = A2(1)
            return
         endif
         call compute_process_rcl(11,pmat(:,1:4),'NLO',A2)
         born    = A2(1)
         oneloop = A2(2) +
     >        born*(1.d0+1.d0*soft_integral_gg(egmax,p1,p2,p3,p4))
      elseif(ng.eq.1) then
         call compute_process_rcl(12,pmat(:,1:5),'LO',A2)
         born    = A2(1)
      elseif(ng.eq.2) then
         call compute_process_rcl(13,pmat(:,1:6),'LO',A2)
         born    = A2(1)
      elseif(ng.eq.3) then
         call compute_process_rcl(14,pmat(:,1:7),'LO',A2)
         born    = A2(1)
      else
         print*,'wrong number of extra photons!!'
      endif

      ELSEIF (fs.eq.'ee') THEN ! BHABHA
      if (ng.eq.0) then

         if (ord.eq.'born') then
            call compute_process_rcl(20,pmat(:,1:4),'LO',A2)
            born = A2(1)
            return
         endif
         call compute_process_rcl(21,pmat(:,1:4),'NLO',A2)
         born    = A2(1)
         oneloop = A2(2) +
     >        born*(1.d0+1.d0*soft_integral_ee(egmax,p1,p2,p3,p4))
      elseif(ng.eq.1) then
         call compute_process_rcl(22,pmat(:,1:5),'LO',A2)
         born    = A2(1)
      elseif(ng.eq.2) then
         call compute_process_rcl(23,pmat(:,1:6),'LO',A2)
         born    = A2(1)
      else

         print*,'wrong number of extra photons!!'
      endif
      ELSEIF (fs.eq.'mm') THEN ! mu+mu-
         if (ord.eq.'born') then
            call compute_process_rcl(30,pmat(:,1:4),'LO',A2)
            born = A2(1)
            return
         endif
         call compute_process_rcl(31,pmat(:,1:4),'NLO',A2)
         born    = A2(1)
         oneloop = A2(2) +
     >        born*(1.d0+1.d0*soft_integral_muons(egmax,p1,p2,p3,p4))

      ELSEIF (fs.eq.'er') THEN
      if(ng.eq.1) then
            if(ord.eq.'born') then
		    call compute_process_rcl(90,pmat(:,1:5),'LO',A2)
		    born = A2(1)
                   return
            endif
            call compute_process_rcl(91,pmat(:,1:5),'NLO',A2)
            born    = A2(1)
            oneloop = A2(2) +
     >           born*(1.d0+1.d0*soft_integral_ee(egmax,p1,p2,p3,p4))
	    return
      elseif(ng.eq.2) then
            call compute_process_rcl(70,pmat(:,1:6),'LO',A2)
            born = A2(1)
            return
         endif


      ELSEIF (fs.eq.'mr') THEN
      if(ng.eq.1) then
            if(ord.eq.'born') then
		    call compute_process_rcl(80,pmat(:,1:5),'LO',A2)
		    born = A2(1)
                   return
                endif
            call compute_process_rcl(81,pmat(:,1:5),'NLO',A2)

            born    = A2(1)
            oneloop = A2(2)
     >        +born*(1.d0+1.d0*soft_integral_muons(egmax,p1,p2,p3,p4))

	    return
         elseif(ng.eq.2) then
            call compute_process_rcl(99,pmat(:,1:6),'LO',A2)
            born    = A2(1)
c$$$            oneloop = A2(2) +
c$$$     >           born*(1.d0+1.d0*soft_integral_muons(egmax,p1,p2,p3,p4))

            return
         elseif(ng.eq.3) then
            call compute_process_rcl(56,pmat(:,1:7),'LO',A2)
            born = A2(1)
            return
      elseif(ng.eq.4) then
            call compute_process_rcl(57,pmat(:,1:8),'LO',A2)
            born = A2(1)
            return
      endif


!*****************************************  MUONI RADIATIVI STATO INIZIALE
!      ELSEIF (fs.eq.'mi') THEN
!      if(ng.eq.1) then
!            if(ord.eq.'born') then
!		    call bk_me(p1,p2,p3,p4,p5,born,qphmat)
!                   return
!            endif
!            call bk_me(p1,p2,p3,p4,p5,born,qphmat)
!            call virtualISR_mur(p1,p2,p3,p4,p5,virtual,qphmat)

!c            print*,'oneloop recola', A22/(scaling**2)
!c            print*,'born',born
!c            print*,'virtual',virtual
!            oneloop = virtual
!     >        +born*(1.d0+1.d0*soft_integral_muons(egmax,p1,p2,p3,p4))
!c            print*,'sv factor',oneloop/born
!	    return
!      elseif(ng.eq.2) then
!            call realhard_mur(p1,p2,p3,p4,p5,born,qph)
!            return
!************************************
 !     endif




      ENDIF
#endif
      return
      end
***********************************************************
      subroutine martella(p1,p2,p3,p4,p5,p3m,p4m,p5m)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      dimension p1m(0:3),p2m(0:3),p3m(0:3),p4m(0:3),p5m(0:3)
      real*16 qp1(0:3),qp2(0:3),qp3(0:3),qp4(0:3),qp5(0:3)
      real*16 qtridot,qmod


      p3m(0) = sqrt(tridot(p3,p3))
      p3m(1:3) = p3(1:3)
      p4m(0) = sqrt(tridot(p4,p4))
      p4m(1:3) = p4(1:3)
      p5m(0) = sqrt(tridot(p5,p5))
      p5m(1:3) = p5(1:3)

      return

      qp1 = p1
      qp2 = p2

      qmod = sqrt(qtridot(1.q0*p3,1.q0*p3))
      qp3(1:3) = p3(0)*p3(1:3)/qmod
      qp3(0) = p3(0)

      qmod = sqrt(qtridot(1.q0*p5,1.q0*p5))
      qp5(1:3) = p5(0)*p5(1:3)/qmod
      qp5(0) = p5(0)

      qp4 = qp1+qp2 - qp3 - qp5

      p3m = qp3
      p4m = qp4
      p5m = qp5
      return
      end

      real*16 function qtridot(qp1,qp2)
      real*16 qp1(0:3),qp2(0:3)
      qtridot = qp1(1)*qp2(1)+qp1(2)*qp2(2)+qp1(3)*qp2(3)
      return
      end
