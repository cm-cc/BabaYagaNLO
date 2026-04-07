program WZ_ew_test
  use recola
  implicit none

  integer, parameter :: dp = 8
  integer            :: i,j, active_flavours, n_loops, ai, pol, psp
  real (kind=dp)     :: A2(0:1),B2(0:1),M2, N2,A2u(0:1),B2u(0:1),M2u, N2u
  real (kind=dp)     :: MW, MH, MZ, MT, GW, GZ, GH, GT, MZ_pwg
  real (kind=dp)     :: MWo, MZo, GWo, GZo, muRen
  real (kind=dp)     :: alpha, Gf, pi, muUV, muIR, as_used, DeltaIR2
  real (kind=dp)     :: s35, s34, s45, s345, z, y, split, pref, cos35, vartc, eik
  real (kind=dp),  dimension(0:3,7) :: pr
  real (kind=dp),  dimension(0:3)   :: prw, prp_g, prp_l, prp_v
  real (kind=dp),  dimension(0:3,6) :: pb
  real (kind=dp)     :: alp, delta, pw2, el, ev, eg, pwng, pwnl, ngnl, beta, ang, enl, eng
  real (kind=dp)     :: soft_limit_check, soft_limit_check_acc

  pi =3.1415926535897932385d+00

  ! select the phase-space point
  write(*,*)" "
  write(*,*)" Recola1 test code for u d~ -> W (l v a) Z(l'l') in the limit of a soft photon, for unpolarised or longitudinal W"
  write(*,*)" --------------------------------------------------"
  write(*,*)" choose a phase-space point [insert integer between 1 and 7, especially run 6 or 7]: "
  !read(*,*) psp

  ! input parameters
  MW =   80.385d+00 !   pole values
  GW =   2.0850000000000000d+00 !   pole values
  MZ =   91.1876d+00 !   pole values ! MZ_pwg =   91.2d+00
  GZ =   2.4952d+00 !   pole values
  Gf =   1.1663700000000000d-05
  MH =   125.00000000000000d+00
  GH =   4.07d-03
  GT =   1.3600000000000000d+00
  MT =   173.00000000000000d+00
  muUV  = 1d+02
  muIR   = 1d+02
  DeltaIR2 = 0d0 !pi**2/6d0
  muRen = MZ
  n_loops=2
  active_flavours=5
  alpha  = (sqrt(2d0)/pi)*Gf*(MW**2)*(1d0-(MW**2/MZ**2))  ! -> 7.562339228910d-03

  ! setting Recola
  call set_output_file_rcl('*')
  call set_complex_mass_scheme_rcl
  call set_mu_ir_rcl(muIR)
  call set_mu_uv_rcl(muUV)
  call set_delta_ir_rcl(0d0,DeltaIR2)
  call use_dim_reg_soft_rcl
  call use_gfermi_scheme_rcl(a=alpha)
  call set_pole_mass_W_rcl(MW, GW)
  call set_pole_mass_Z_rcl(MZ, GZ)
  call set_pole_mass_top_rcl(MT, GT)
  call set_pole_mass_h_rcl(MH, GH)
  call set_pole_mass_bottom_rcl(0d0,0d0)
  call set_light_fermions_rcl(1d-3)
  call switchoff_resonant_selfenergies_rcl
  call set_resonant_particle_rcl('W+')
  call set_resonant_particle_rcl('Z')

  ! setting and running alphaS
  call set_alphas_rcl(0.118d+00, MZ, int(active_flavours)) ! useless for pure EW processes

  ! generating processes
  call define_process_rcl(1,"u d~ -> W+[0](e+ nu_e) Z(mu+ mu-)","LO") !  p p > z w, a^4, Born-level (LO)         ! longitudinal
  call define_process_rcl(2,"u d~ -> W+[0](e+ nu_e A*) Z(mu+ mu-)","LO") !  p p > z w a, a^5, Real-level (NLOEW)  ! longitudinal
  call define_process_rcl(3,"u d~ -> W+[4](e+ nu_e) Z(mu+ mu-)","LO") !  p p > z w, a^4, Born-level (LO)         ! unpolarised
  call define_process_rcl(4,"u d~ -> W+[4](e+ nu_e A*) Z(mu+ mu-)","LO") !  p p > z w a, a^5, Real-level (NLOEW)  ! unpolarised
  call unselect_all_gs_powers_BornAmpl_rcl(1)
  call unselect_all_gs_powers_BornAmpl_rcl(2)
  call unselect_all_gs_powers_BornAmpl_rcl(3)
  call unselect_all_gs_powers_BornAmpl_rcl(4)
  call select_gs_power_BornAmpl_rcl(1,0)
  call select_gs_power_BornAmpl_rcl(2,0)
  call select_gs_power_BornAmpl_rcl(3,0)
  call select_gs_power_BornAmpl_rcl(4,0)

  !#############################
  pol = 0 ! longitudinal W boson

  !#############################
  call generate_processes_rcl()

  do psp=1, 7, 1
  ! initialising (n+1)-body momenta (bosons already projected onto mass shell, momenta written in the partonic CM frame)
  select case (psp)
  case(1)
     pr(:,1) = [   580.70311015858420d0,   2.3685309447555420d-012 ,  -3.4536728753924421d-012,   580.70311015858420d0]
     pr(:,2) = [   580.70311015858147d0,   4.0831970959461140d-012 ,  -5.9539129460712052d-012,  -580.70311015858135d0]
     pr(:,3) = [   53.574672231430313d0,  -34.765475875099668d0   ,   39.980529006998971d0   ,  -7.9476092130631208d0 ]
     pr(:,4) = [   287.90146770911457d0,  -127.65125058597040d0   ,   242.41341398879831d0   ,  -88.476833415282925d0 ]
     pr(:,5) = [   238.42904519562697d0,  -156.15784724880268d0   ,   176.02591711059691d0   ,  -38.444932577534360d0 ]
     pr(:,6) = [   379.41200897877366d0,   245.13126211635313d0   ,  -277.26540771176417d0   ,   83.594441065073156d0 ]
     pr(:,7) = [   202.08902620222904d0,   73.443311593519667d0   ,  -181.15445239463006d0   ,   51.274934140807197d0 ]
  case(2)
     pr(:,1) = [   441.03640767982438d0,  -4.0311839678688531d-013 ,-3.3687562652836625d-014 ,   441.03640767982438d0 ]
     pr(:,2) = [   441.03640767982614d0,  -8.6164767286609590d-013 ,-7.2005669291478517d-014 ,  -441.03640767982597d0 ]
     pr(:,3) = [   240.48446218084578d0,   146.88468440082676d0     ,   8.5350321615845939d0 ,   190.22307763486762d0 ]
     pr(:,4) = [   25.446360201901882d0,   18.536688063074269d0     ,   15.584982416442832d0 ,   7.8113229517407898d0 ]
     pr(:,5) = [   174.05497465786095d0,   106.59650225154101d0     ,  -1.8968721390374712d0 ,   137.58169132182513d0 ]
     pr(:,6) = [   340.73952050658244d0,  -185.29274496145118d0     ,  -32.209922021549104d0 ,  -284.13472231156811d0 ]
     pr(:,7) = [   101.34749781245822d0,  -86.725129753990871d0     ,   9.9867795825591514d0 ,  -51.481369596865449d0 ]
  case(3)
     pr(:,1) = [  141.97507438389317d+00,  5.8660379351970675d-15,  -6.1448267515970070d-15,   141.97507438389317d+00 ]
     pr(:,2) = [  141.97507438389303d+00,  1.1873244330141752d-14,  -1.2437531123066544d-14,  -141.97507438389303d+00 ]
     pr(:,3) = [  51.912967524999061d+00,  8.9410831625192753d+00,  -45.323181908934600d+00,   23.681689356592472d+00 ]     ! emitter (l+)
     pr(:,4) = [  86.757389147705155d+00,  16.769230870912509d+00,   18.570733462325382d+00,   83.070845224998465d+00 ]     ! spectator (vl)
     pr(:,5) = [  4.1063550761411965d-02,  2.3076765571936192d-02,  -2.9710440918878939d-02,   1.6461099364221167d-02 ]     ! emissus (photon)
     pr(:,6) = [  77.245517925600126d+00, -24.563796776672682d+00,  -28.643158399386166d+00,  -67.402221064414547d+00 ]
     pr(:,7) = [  67.993210618720568d+00, -1.1695940223310386d+00,   55.425317286914265d+00,  -39.366774616540596d+00 ]
  case(4)
     pr(:,1) = [   85.800045598689465d0,   7.8663857618788702d-016 ,  -1.7306048676133515d-015,   85.800045598689437d0]
     pr(:,2) = [   85.800045598689465d0,   3.2358444843726966d-016 ,  -7.1188578656199329d-016,  -85.800045598689465d0]
     pr(:,3) = [   37.747916268709268d0,   2.5420791416736921d0   ,  -33.745253977684023d0   ,   16.724259333332355d0 ]
     pr(:,4) = [   39.988380993916159d0,  -3.6185741327689476d0   ,   35.210572170708033d0   ,  -18.606239361343622d0 ]
     pr(:,5) = [   2.6633136514805962d0, -0.43429854329065631d0   ,  -1.6869067172084276d0   ,   2.0146885885653933d0 ]
     pr(:,6) = [   45.892233229072446d0,   12.375340721758654d0   ,   43.937908820113918d0   ,   4.7337280537984654d0 ]
     pr(:,7) = [   45.308247054200379d0,  -10.864547187372743d0   ,  -43.716320295929499d0   ,  -4.8664366143525797d0 ]
  case(5)
     pr(:,1) = [ 163.77772872699759d0,  -3.4483804990903599d-015 ,  -2.8628064520750161d-015 ,   163.77772872699757d0 ]
     pr(:,2) = [ 163.77772872699754d0,  -2.4358015314229682d-015 ,  -2.0221748562756719d-015 ,  -163.77772872699754d0 ]
     pr(:,3) = [ 35.166387460939710d0,  -26.020547408876524d0    ,  -23.606455384376495d0    ,   1.5300927360094043d0 ]
     pr(:,4) = [ 125.22057376432562d0,   54.257369830470282d0    ,  -108.30176729485785d0    ,   31.734163196103566d0 ]
     pr(:,5) = [ 0.56158220412461057d0, -0.14720932414785304d0   , -0.53359761450296128d0    , 9.4750053668255313d-02 ]
     pr(:,6) = [ 48.319416702383769d0,   17.799138522285574d0    ,   38.848467079039821d0    ,   22.555560377251787d0 ]
     pr(:,7) = [ 118.28749732222148d0,  -45.888751619731472d0    ,   93.593353214697487d0    ,  -55.914566363032989d0 ]
  case(6)
     pr(:,1) = [   108.21229194929980d0,  -4.7038049030811041d-015 ,   4.7038049030811041d-015 , 108.21229194929977d0 ]
     pr(:,2) = [   108.21229194929974d0,  -2.4016224545199006d-015 ,   2.4016224545199006d-015 ,-108.21229194929974d0 ]
     pr(:,3) = [   68.586147008009917d0,   3.4174682632878302d0,  -39.211386178474697d0      ,  -56.168030640536898d0 ]
     pr(:,4) = [   8.8695543849586329d0,  -1.2461329603724387d0,   3.4579757105855706d0      ,   8.0720847132455322d0 ]
     pr(:,5) = [   26.474659790026994d0,  -21.072950533470298d0,   12.162660751439994d0      ,  -10.435901985787691d0 ]
     pr(:,6) = [   49.037766913388076d0,   47.445086037352795d0,  -6.5401327421904369d0      ,   10.530577309649161d0 ]
     pr(:,7) = [   63.456455802215892d0,  -28.543470806797888d0,   30.130882458639565d0      ,   48.001270603429901d0 ]
  case(7)
     pr(:,1) = [   98.846015976326612d0,  -7.0504138457527723d-015 ,   7.1050682166500815d-16 , 98.846015976326612d0  ]
     pr(:,2) = [   98.846015976326612d0,  -7.2714631719117479d-015 ,   7.3278311034769555d-16 ,-98.846015976326612d0  ]
     pr(:,3) = [   31.321200797531482d0,  -3.3253980410044890d0 ,   21.023559717068760d0 ,  -22.977582207251796d0     ]
     pr(:,4) = [   61.645221570554426d0,   51.327261473323411d0 ,  -23.102367914921334d0 ,   25.138141714443677d0     ]
     pr(:,5) = [   1.1919232189659954d0,  0.88627490891029481d0 , -0.73761824934250397d0 ,  0.30185603188143945d0     ]
     pr(:,6) = [   71.043686350504771d0,  -66.130020263174600d0 ,   8.8572376640992090d0 ,   24.404408028523761d0     ]
     pr(:,7) = [   32.490000015096534d0,   17.241881921945389d0 ,  -6.0408112169041326d0 ,  -26.866823567597081d0     ]
  end select

  ! initialising (n)-body momenta
  prw(:) = pr(:,3) + pr(:,4) + pr(:,5)
  pb(:,1) = [  0d0, 0d0, 0d0, 0d0 ]
  pb(:,2) = [  0d0, 0d0, 0d0, 0d0 ]
  pb(:,3) = [  0d0, 0d0, 0d0, 0d0 ]    ! mapped emitter (l+)
  pb(:,4) = [  0d0, 0d0, 0d0, 0d0 ]    ! mapped spectator (vl)
  pb(:,5) = [  0d0, 0d0, 0d0, 0d0 ]
  pb(:,6) = [  0d0, 0d0, 0d0, 0d0 ]

  ! #########################################################################
  ! now approaching the soft-photon limit starting from the selected real PSP
  ! #########################################################################
  alp   = 1d-00
  write(*,*) " "
  write(*,*) " ------------------------------------------------------------------------------------------"
  write(*,*) " ###   phase-space point Nr. ", psp
  write(*,*) " ###   R/CT = |A(n+1)|^2 / ( K . |A(n)|^2 ) should approach 1 in the soft-photon limit "
  write(*,*) " ###   K is the eikonal kernel in the soft-photon limit "
  write(*,*) " ###   the soft-photon limit is approached keeping the angle between positron and photon fixed "
  write(*,*) " ---------------"
  write(*,*) " energy(A)    ", &
       "             costheta(e+,A)     ", &
       "        energy(e+)    " ,&
       "           R/CT(longit)" , &
       "              R/CT(unpol)"
  write(*,*) " ---------------"


  ! ####################### now loop over descreasing photon energy
  !--------------------------------------------------------------------------
  do ai = 1,8
  !--------------------------------------------------------------------------
  alp   = (2d0)**(-ai)
  delta = pr(0,5)*(1d0-alp)  ! E'(gamma) = alp*pr(0,5)
  pw2   = prw(1)**2+prw(2)**2+prw(3)**2
  el    = pr(0,3)
  ev    = pr(0,4)
  eg    = pr(0,5)
  pwng  = (prw(1) * pr(1,5) + prw(2) * pr(2,5) + prw(3) * pr(3,5)) / &
       sqrt(pr(1,5)**2 + pr(2,5)**2 + pr(3,5)**2)
  pwnl  = (prw(1) * pr(1,3) + prw(2) * pr(2,3) + prw(3) * pr(3,3)) / &
       sqrt(pr(1,3)**2 + pr(2,3)**2 + pr(3,3)**2)
  ngnl  = (pr(1,3) * pr(1,5) + pr(2,3) * pr(2,5) + pr(3,3) * pr(3,5)) / &
       sqrt(pr(1,5)**2 + pr(2,5)**2 + pr(3,5)**2) / &
       sqrt(pr(1,3)**2 + pr(2,3)**2 + pr(3,3)**2)
  beta = (2d0 * pwng * delta + pw2 - 2d0 * pwnl * el - 2d0 * ngnl * delta * el + el**2 - 2d0 * delta * ev &
       - ev**2 - 2d0 * pwng * eg  - 2d0 * delta * eg + 2d0 * ngnl * el * eg + eg**2)/&
       ( 2d0 * delta * ( pwnl - delta + ngnl * delta - el - ev - ngnl * eg ) )
  prp_g(:) = alp * pr(:,5)
  prp_l(:) = (( pr(0,3) + beta * delta) / pr(0,3)) * pr(:,3)
  prp_v(:) = prw(:) - prp_g(:) - prp_l(:)
  pr(:,3) = prp_l(:)
  pr(:,4) = prp_v(:)
  pr(:,5) = prp_g(:)
  ang  = (pr(1,3) * pr(1,5) + pr(2,3) * pr(2,5) + pr(3,3) * pr(3,5)) / &
        sqrt(pr(1,5)**2 + pr(2,5)**2 + pr(3,5)**2) / &
        sqrt(pr(1,3)**2 + pr(2,3)**2 + pr(3,3)**2)
  enl = pr(0,3)
  eng = pr(0,5)
  ! invariants and angles
  s35 = 2d0*(pr(0,3)*pr(0,5)-pr(1,3)*pr(1,5)-pr(2,3)*pr(2,5)-pr(3,3)*pr(3,5))
  s34 = 2d0*(pr(0,3)*pr(0,4)-pr(1,3)*pr(1,4)-pr(2,3)*pr(2,4)-pr(3,3)*pr(3,4))
  s45 = 2d0*(pr(0,4)*pr(0,5)-pr(1,4)*pr(1,5)-pr(2,4)*pr(2,5)-pr(3,4)*pr(3,5))
  s345 = s35+s34+s45
  cos35 = (pr(1,3)*pr(1,5)+pr(2,3)*pr(2,5)+pr(3,3)*pr(3,5))/&
       sqrt(pr(1,3)**2+pr(2,3)**2+pr(3,3)**2)/&
       sqrt(pr(1,5)**2+pr(2,5)**2+pr(3,5)**2)
  ! CS variables, kernel (supposed to mimic the real in soft and collinear limits) and prefactor
  y = s35/s345
  z = s34/(s34+s45)
  vartc = (s35+s45)/s345
  split = (1d0/(1d0-y))*((2d0/(1d0-z*(1d0-y)))-(1d0+z)-(2d0*y/(1d0-z*(1d0-y))**2))
  pref = (-1d0/s35)*(-1d0)*(8d0*pi*alpha) ! (-1/sij)*(casimirfundam)*(8 pi alphaEW)
  eik  = (16d0*pi*alpha)*s345*s45/(s35*(s35+s45)**2)
  ! constructing CS-mapped (n)-body momenta
  do i = 0,3
     pb(i,3) = pr(i,3) + pr(i,5) - pr(i,4) * (s35/(s34+s45))
     pb(i,4) = pr(i,4) * (s345/(s34+s45))
  enddo
  do i = 0,3
     pb(i,1) = pr(i,1)
     pb(i,2) = pr(i,2)
     pb(i,5) = pr(i,6)
     pb(i,6) = pr(i,7)
  enddo
  ! calculate Recola squared amplitudes
  call compute_running_alphas_rcl(0.5d0*(MZ+MW),int(active_flavours),int(n_loops))
  call compute_process_rcl(1,pb, "LO",A2)
  call get_squared_amplitude_rcl(1,0,  "LO",M2)
  call compute_process_rcl(2,pr, "LO",B2)
  call get_squared_amplitude_rcl(2,0,  "LO",N2)
  call compute_process_rcl(3,pb, "LO",A2u)
  call get_squared_amplitude_rcl(3,0,  "LO",M2u)
  call compute_process_rcl(4,pr, "LO",B2u)
  call get_squared_amplitude_rcl(4,0,  "LO",N2u)

  !write(*,*) eng, ang, enl, N2/(M2*pref*split), N2u/(M2u*pref*split) ! this uses a similar kernel CT
  write(*,*) eng, ang, enl, N2/(M2*eik), N2u/(M2u*eik) ! this uses the simple eikonal kernel


  if ((ai .eq. 5) .or. (ai .eq. 6) .or. (ai .eq. 7)) then
      soft_limit_check = abs(N2/(M2*eik) - 1d0)

      ! we check that it gets better, but then numerical accuacy takes over
      if (ai .eq. 5) soft_limit_check_acc = 2d-4
      if (ai .eq. 6) soft_limit_check_acc = 5d-5
      if (ai .eq. 7) soft_limit_check_acc = 3d-4
      if (soft_limit_check .gt. soft_limit_check_acc) then
        call EXIT(1)
      end if
  end if

  if ((ai .eq. 5) .or. (ai .eq. 6) .or. (ai .eq. 7)) then
      soft_limit_check = abs(N2u/(M2u*eik) - 1d0)
      if (ai .eq. 5) soft_limit_check_acc = 2d-4
      if (ai .eq. 6) soft_limit_check_acc = 2d-5
      if (ai .eq. 7) soft_limit_check_acc = 3d-4
      if (soft_limit_check .gt. soft_limit_check_acc) then
        call EXIT(1)
      end if
  end if

  enddo
  write(*,*) " ---------------"
  write(*,*) " "
  !--------------------------------------------------------------------------
  enddo

  call reset_recola_rcl()

end program WZ_ew_test
