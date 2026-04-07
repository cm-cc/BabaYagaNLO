!#####################################################################
!!
!!  File  check_rcl.f90
!!  is part of RECOLA (REcursive Computation of One Loop Amplitudes)
!!
!!  Copyright (C) 2015-2025   Stefano Actis, Ansgar Denner,
!!                            Lars Hofer, Jean-Nicolas Lang,
!!                            Andreas Scharf, Sandro Uccirati
!!
!!  RECOLA is licenced under the GNU GPL version 3,
!!         see COPYING for details.
!!
!#####################################################################

  module check_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl
  use process_definition_rcl
  use process_generation_rcl
  use process_computation_rcl
  use reset_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine check_process_rcl (processIn,delta)

  ! This routine runs recola for the checked process "processIn".
  ! The rules to write the character variable "processIn" are the
  ! same as in define_process_rcl.
  ! The 4-legs processes (EW) checked with the code Pole are:
  ! 'u~ u -> nu_e~ nu_e'
  ! 'u d~ -> nu_e e+'
  ! 'e+ e- -> nu_e~ nu_e'
  ! 'e+ e- -> W+ W-'
  ! 'e+ e- -> Z H'
  ! The 4-legs processes (EW+QCD) checked with the code Pole are:
  ! 'd~ d -> u~ u'
  ! 'u d~ -> W+ H'
  ! The 4-legs processes (QCD) checked with the code OpenLoops are:
  ! 'u d~ -> W+ g'
  ! 'g g -> g g'
  ! 'b b~ -> t t~'
  ! 'g g -> b b~'
  ! The 5-legs processes (EW) checked with the code Pole are:
  ! 'u d~ -> e+ nu_e A'
  ! The 5-legs processes (QCD) checked with the code OpenLoops are:
  ! 'u u~ -> W+ W- g'
  ! 'u u~ -> Z Z g'
  ! 'u u~ -> Z A g'
  ! 'u u~ -> A A g'
  ! 'u d~ -> W+ g g'
  ! 'u d~ -> W+ t t~'
  ! 'u u~ -> Z t t~'
  ! 'd d~ -> Z t t~'
  ! 'g g -> W+ b t~'
  ! 'g g -> Z t t~'
  ! 'u u~ -> Z g g'
  ! 'g g -> g t t~'
  ! 'g g -> g g g'
  ! 'd d~ -> d d~ g'
  ! 'd d~ -> t t~ g'
  ! The 6-legs processes (QCD) checked with the code OpenLoops are:
  ! 'u d~ -> W+ g g g'
  ! 'u u~ -> Z g g g'
  ! 'u u~ -> W+ W- g g'
  ! 'u u~ -> Z Z g g'
  ! 'd d~ -> t t~ b b~'
  ! 'g g -> t t~ b b~'
  ! 'u u~ -> u u~ u u~'
  ! 'g g -> u u~ u u~'
  ! 'g g -> u u~ d d~'
  ! 'g g -> u u~ g g'
  ! 'g g -> t t~ g g'
  ! If a process has not yet been checked, an error message appears
  ! and the program stops.

  character(len=*), intent(in)    :: processIn
  real(dp), intent(out), optional :: delta

  integer               :: i,n,legs,gs0,gs1
  integer, allocatable  :: gs(:,:)
  real(dp)              :: p4a(0:3,4),matrix2Psum(0:1),mat2r,mat2c,del
  real(dp), allocatable :: p(:,:),matrix2P(:,:)
  character             :: order*3,check*10,cpr*99,           &
                           fmt1i*99,fmt2i*99,fmt1*99,fmt2*99, &
                           fmt(0:3)*6,fmtTot*80

  ! variables for all processes

  mass_z  =  91.154892493050440d0; width_z  = 2.4421237356890875d0
  mass_w  =  80.371978311278895d0; width_w  = 2.0842992422669000d0
  mass_h  = 120.d0;                width_h  = 0.d0
  mass_el =   0.d0
  mass_mu =   0.d0;                width_mu = 0.d0
  mass_ta =   0.d0;                width_ta = 0.d0
  mass_u  =   0.d0
  mass_d  =   0.d0
  mass_c  =   0.d0;                width_c  = 0.d0
  mass_s  =   0.d0
  mass_t  = 172.6d0;               width_t  = 0.d0
  mass_b  =   0.d0;                width_b  = 0.d0

  light_el = .true.
  light_mu = .true.
  light_ta = .true.
  light_u  = .true.
  light_d  = .true.
  light_c  = .true.
  light_s  = .true.
  light_t  = .false.
  light_b  = .false.

  reg_soft    = 1; lambda = 1d0

  DeltaUV = 0d0
  DeltaIR = 0d0; DeltaIR2 = 0d0
  muUV = mass_w
  muIR = mass_w

  complex_mass_scheme = 1

  ew_reno_scheme = 1
  gf  = 1.16637d-5
  algf = - 1d0

  masses_in_gf_reno_scheme = 0

  als = 0.12027511899498215d0
  Qren = mass_w
  Nfren = 5

  loopQED  = .true.
  loopWEAK = .true.
  pureQED = .false.
  reguScheme = 2
  check_Pole = .true.
  resIR = .false.
  longitudinal = 0

  resSE = .false.

  dynamic_settings = 0

  momenta_correction = .true.

  colour_optimization = 2
  helicity_optimization = 2
  writeMat  = 0
  writeMat2 = 2
  writeCor  = 0
  draw      = 0

  ifail = 0

  order = 'NLO'

  gs0 = 111
  gs1 = 111

  p4a(0,1) = + 4000d0
  p4a(1,1) =   0d0
  p4a(2,1) =   0d0
  p4a(3,1) = + 4000d0
  p4a(0,2) = + 4000d0
  p4a(1,2) =   0d0
  p4a(2,2) =   0d0
  p4a(3,2) = - 4000d0
  p4a(0,3) = + 4000d0
  p4a(1,3) = + 2144.3833090906010d0
  p4a(2,3) = + 3373.7194807643473d0
  p4a(3,3) = + 140.13239741325378d0
  p4a(0,4) = + 4000d0
  p4a(1,4) = - 2144.3833090906010d0
  p4a(2,4) = - 3373.7194807643473d0
  p4a(3,4) = - 140.13239741325378d0

  ! process

  cpr = adjustl(processIn)//'END'
  n = 99
  do while (index (cpr,'  ').ne.index (cpr,'END')+3)
    n = index (cpr,'  ')
    cpr = cpr(1:n-1)//cpr(n+1:)
  enddo
  n = 99
  do while (index (cpr,' [').ne.0)
    n = index (cpr,' [')
    cpr = cpr(1:n-1)//cpr(n+1:)
  enddo
  n = 99
  do while (index (cpr,' (').ne.0)
    n = index (cpr,' (')
    cpr = cpr(1:n-1)//cpr(n+1:)
  enddo
  n = 99
  do while (index (cpr,'( ').ne.0)
    n = index (cpr,'( ')
    cpr = cpr(1:n)//cpr(n+2:)
  enddo
  n = 99
  do while (index (cpr,' )').ne.0)
    n = index (cpr,' )')
    cpr = cpr(1:n-1)//cpr(n+1:)
  enddo
  n = index (cpr,'END')
  cpr = cpr(1:n-1)

  matrix2Psum = 0d0

  if     (cpr.eq.'u~ u -> nu_e~ nu_e') then

    check = 'Pole'
    legs = 4
    allocate (gs(0:legs,0:1)); gs = 0; gs(0,0:1) = 1
    allocate(p(0:3,legs)); p = p4a
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(0,0) = 9.57814610849876732d-004
    matrix2P(0,1) = &
    - 7.85646307274683107d-004 & ! boxes
    + 2.48717467405677338d-004 & ! triangles
    + 4.01206502134965519d-005 & ! self-energies
    - 7.38975225297194724d-005   ! counterterms

  elseif (cpr.eq.'u d~ -> nu_e e+') then

    check = 'Pole'
    legs = 4
    allocate (gs(0:legs,0:1)); gs = 0; gs(0,0:1) = 1
    allocate(p(0:3,legs)); p = p4a
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(0,0) = 4.05436322743390424d-003
    matrix2P(0,1) = &
    - 3.34685213428006145d-003 & ! boxes
    + 9.54435703234285480d-004 & ! triangles
    + 1.63047377548910415d-004 & ! self-energies
    - 3.30426293431707518d-004   ! counterterms

  elseif (cpr.eq.'d~ d -> u~ u') then

    check = 'Pole'
    legs = 4
    allocate (gs(0:legs,0:1)); gs = 1
    allocate(p(0:3,legs)); p = p4a
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2Psum(0) = 0.41284064171020179d0
    matrix2Psum(1) = -1.4881832162850928d0

  elseif (cpr.eq.'e+ e- -> nu_e~ nu_e') then

    check = 'Pole'
    legs = 4
    allocate (gs(0:legs,0:1)); gs = 1
    allocate(p(0:3,legs)); p = p4a
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2Psum(0) = 7.26179240617077365d-2
    matrix2Psum(1) = -4.70818596988321139d-2

  elseif (cpr.eq.'e+ e- -> W+ W-') then

    check = 'Pole'
    legs = 4
     mass_w = 80.399d0
    width_w =  0.d0
    muUV = mass_w
    muIR = mass_w
    Qren = mass_w
    allocate (gs(0:legs,0:1)); gs = 1
    allocate(p(0:3,legs))
    p(0,1) = + 4000d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 4000d0
    p(0,2) = + 4000d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 4000d0
    p(0,3) = + 4000d0
    p(1,3) = + 2143.9500999571701d0
    p(2,3) = + 3373.0379206689554d0
    p(3,3) = + 140.10408781291710d0
    p(0,4) = + 4000d0
    p(1,4) = - 2143.9500999571701d0
    p(2,4) = - 3373.0379206689554d0
    p(3,4) = - 140.10408781291710d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2Psum(0) = 2.96674776566679364d-2
    matrix2Psum(1) = -3.40773328495682926d-2

  elseif (cpr.eq.'u d~ -> W+ H') then

    check = 'Pole'
    legs = 4
     mass_w = 80.399d0
    width_w =  0.d0
    muUV = mass_w
    muIR = mass_w
    Qren = mass_w
    allocate (gs(0:legs,0:1)); gs = 1
    allocate(p(0:3,legs))
    p(0,1) = + 4000d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 4000d0
    p(0,2) = + 4000d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 4000d0
    p(0,3) = + 3999.5039999500627d0
    p(1,3) = + 2143.6841426646142d0
    p(2,3) = + 3372.6194948701982d0
    p(3,3) = + 140.08670788235273d0
    p(0,4) = + 4000.4960000499377d0
    p(1,4) = - 2143.6841426646142d0
    p(2,4) = - 3372.6194948701982d0
    p(3,4) = - 140.08670788235273d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2Psum(0) = 1.86729421127388080d-3
    matrix2Psum(1) = -3.93591538157348972d-3

  elseif (cpr.eq.'e+ e- -> Z H') then

    check = 'Pole'
    legs = 4
     mass_z = 91.1876d0
    width_z =  0.d0
    allocate (gs(0:legs,0:1)); gs = 1
    allocate(p(0:3,legs))
    p(0,1) = + 4000d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 4000d0
    p(0,2) = + 4000d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 4000d0
    p(0,3) = + 3999.6196986496102d0
    p(1,3) = + 2143.6220891503181d0
    p(2,3) = + 3372.5218671983439d0
    p(3,3) = + 140.08265277350660d0
    p(0,4) = + 4000.3803013503903d0
    p(1,4) = - 2143.6220891503181d0
    p(2,4) = - 3372.5218671983439d0
    p(3,4) = - 140.08265277350660d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2Psum(0) = 2.37149351547372590d-3
    matrix2Psum(1) = -1.63698726254942113d-3

  elseif (cpr.eq.'u d~ -> W+ g') then

    check = 'openLoops'
    legs = 4
    ew_reno_scheme = 3; alZ = 0.0078125d0
     mass_z = 91.1876d0
    width_z =  0.d0
     mass_w = 80.399d0
    width_w =  0.d0
    complex_mass_scheme = 0
    muUV = 100d0
    muIR = 100d0
    Qren = 100d0
    allocate (gs(0:legs,0:1)); gs = 0; gs(1,0) = 1; gs(3,1) = 1
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) = + 503.2319996005000462d0
    p(1,3) = + 302.3147316012482406d0
    p(2,3) = + 163.3466058461748389d0
    p(3,3) = + 358.7507987953716793d0
    p(0,4) = + 496.7680003995000106d0
    p(1,4) = - 302.3147316012482406d0
    p(2,4) = - 163.3466058461748105d0
    p(3,4) = - 358.7507987953716793d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(2,0) = 0.95821667535212d0
    matrix2P(4,1) = -0.15278750780672d0

  elseif (cpr.eq.'g g -> g g') then

    check = 'openLoops'
    legs = 4
    complex_mass_scheme = 0
    muUV = 100d0
    muIR = 100d0
    Qren = 100d0
    allocate (gs(0:legs,0:1)); gs = 1
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) = + 500d0
    p(1,3) = + 304.28160767010678d0
    p(2,3) = + 164.4093477385940218d0
    p(3,3) = + 361.0848509836230278d0
    p(0,4) = + 500d0
    p(1,4) = - 304.28160767010678d0
    p(2,4) = - 164.4093477385940218d0
    p(3,4) = - 361.0848509836230278d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2Psum(0) = 0.24511877831620d3
    matrix2Psum(1) = -0.97387038866635d2

  elseif (cpr.eq.'b b~ -> t t~') then

    check = 'openLoops'
    legs = 4
    gs0 = 2
    gs1 = 4
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
!     mass_b =   5.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    if (mass_b.eq.5d0) then
      p(0,1) = + 500d0
      p(1,1) =   0d0
      p(2,1) =   0d0
      p(3,1) = + 499.9749993749688d0
      p(0,2) = + 500d0
      p(1,2) =   0d0
      p(2,2) =   0d0
      p(3,2) = - 499.9749993749688d0
      p(0,3) =   500.0d0
      p(1,3) =   285.7111940687087d0
      p(2,3) =   154.3753873858148d0
      p(3,3) =   339.0477154521355d0
      p(0,4) =   500.0d0
      p(1,4) = - 285.7111940687087d0
      p(2,4) = - 154.37538738581478d0
      p(3,4) = - 339.0477154521355d0
      allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
      matrix2P(gs0+gs0,0) = + 8.765797797519910d-01
      matrix2P(gs0+gs1,1) = + 8.652605119737087d-01
   elseif (mass_b.eq.0d0) then
      p(0,1) = + 500d0
      p(1,1) =   0d0
      p(2,1) =   0d0
      p(3,1) = + 500d0
      p(0,2) = + 500d0
      p(1,2) =   0d0
      p(2,2) =   0d0
      p(3,2) = - 500d0
      p(0,3) =   500.0d0
      p(1,3) =   285.7111940687087d0
      p(2,3) =   154.3753873858148d0
      p(3,3) =   339.0477154521355d0
      p(0,4) =   500.0d0
      p(1,4) = - 285.7111940687087d0
      p(2,4) = - 154.37538738581478d0
      p(3,4) = - 339.0477154521355d0
      allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
      matrix2P(gs0+gs0,0) = + 8.765497762283978d-01
      matrix2P(gs0+gs1,1) = - 2.855840422350293d-01
    endif

  elseif (cpr.eq.'g g -> b b~') then

    check = 'openLoops'
    legs = 4
    gs0 = 2
    gs1 = 4
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
!     mass_b =   5.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    if (mass_b.eq.5d0) then
      p(0,1) = + 500d0
      p(1,1) =   0d0
      p(2,1) =   0d0
      p(3,1) = + 500d0
      p(0,2) = + 500d0
      p(1,2) =   0d0
      p(2,2) =   0d0
      p(3,2) = - 500d0
      p(0,3) = 500.0d0
      p(1,3) = 304.26639320935226d0
      p(2,3) = 164.40112706568513d0
      p(3,3) = 361.0667962896952d0
      p(0,4) = 500.0d0
      p(1,4) = -304.26639320935226d0
      p(2,4) = -164.4011270656851d0
      p(3,4) = -361.0667962896952d0
      allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
      matrix2P(gs0+gs0,0) = 1.936226265893901d+00
      matrix2P(gs0+gs1,1) = 1.173402037004956d+00
   elseif (mass_b.eq.0d0) then
      p(0,1) = + 500d0
      p(1,1) =   0d0
      p(2,1) =   0d0
      p(3,1) = + 500d0
      p(0,2) = + 500d0
      p(1,2) =   0d0
      p(2,2) =   0d0
      p(3,2) = - 500d0
      p(0,3) = 500.0d0
      p(1,3) = 304.2816076701068d0
      p(2,3) = 164.40934773859402d0
      p(3,3) = 361.084850983623d0
      p(0,4) = 500.0d0
      p(1,4) = -304.2816076701068d0
      p(2,4) = -164.409347738594d0
      p(3,4) = -361.084850983623d0
      allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
      matrix2P(gs0+gs0,0) = + 1.936326942827043d+00
      matrix2P(gs0+gs1,1) = - 1.213562645550502d+00
    endif

  elseif (cpr.eq.'u d~ -> e+ nu_e A*'  ) then

    check = 'Pole'
    legs = 5
    allocate (gs(0:legs,0:1)); gs = 0; gs(0,0:1) = 1
    allocate(p(0:3,legs))
    p(0,1) = + 136.03653095124437d0
    p(1,1) =                    0d0
    p(2,1) =                    0d0
    p(3,1) = + 136.03653095124437d0
    p(0,2) = + 136.03653095124437d0
    p(1,2) =                    0d0
    p(2,2) =                    0d0
    p(3,2) = - 136.03653095124437d0
    p(0,3) = + 34.732854958394526d0
    p(1,3) = - 21.932203254832267d0
    p(2,3) = - 20.060898864526393d0
    p(3,3) = - 17.969697011826309d0
    p(0,4) = + 107.33010513629941d0
    p(1,4) = - 104.90816177806244d0
    p(2,4) = + 20.060898864526393d0
    p(3,4) = - 10.563588294902672d0
    p(0,5) = + 130.01010180779483d0
    p(1,5) = + 126.84036503289471d0
    p(2,5) =                    0d0
    p(3,5) = + 28.533285306728981d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(0,0) = 3.85361389052783732d-008
    matrix2P(0,1) = &
    + 1.10925155882286635d-009 & ! pentagons
    - 7.38908298325202700d-009 & ! boxes
    + 1.17612664353575962d-009 & ! triangles
    + 3.76861917337427360d-009 & ! self-energies
    - 3.81189722853211207d-009   ! counterterms

  elseif (cpr.eq.'u u~ -> W+ W- g') then

    check = 'openLoops'
    legs = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0; gs(1,0) = 1; gs(3,1) = 1
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) = + 438.218295761793172d0
    p(1,3) = + 244.1807322977003309d0
    p(2,3) = + 153.8582132409753171d0
    p(3,3) = + 319.804152149755339d0
    p(0,4) = + 495.5814786225607236d0
    p(1,4) = - 299.1286500603476384d0
    p(2,4) = - 150.133911405357054d0
    p(3,4) = - 356.5373221865085611d0
    p(0,5) = +  66.200225615646147d0
    p(1,5) = +  54.947917762647144d0
    p(2,5) = -   3.7243018356182249d0
    p(3,5) = +  36.7331700367532932d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(2,0) = + 2.478709565073970d-04
    matrix2P(4,1) = - 2.605398472478095d-06

  elseif (cpr.eq.'u u~ -> Z Z g') then

    check = 'openLoops'
    legs = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0; gs(1,0) = 1; gs(3,1) = 1
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) = + 417.4206246099386135d0
    p(1,3) = + 127.8541441091167314d0
    p(2,3) = - 194.2513616124706175d0
    p(3,3) = - 334.4316459335413469d0
    p(0,4) = + 155.8236098065406736d0
    p(1,4) = - 50.091070172690884d0
    p(2,4) = + 90.7786247395575003d0
    p(3,4) = - 72.2214992833172431d0
    p(0,5) = + 426.755765583520656d0
    p(1,5) = - 77.7630739364259114d0
    p(2,5) = + 103.4727368729131598d0
    p(3,5) = + 406.6531452168586611d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(2,0) = + 1.190560462998326d-05
    matrix2P(4,1) = - 8.870018925211744d-08

  elseif (cpr.eq.'u u~ -> Z A* g') then

    check = 'openLoops'
    legs = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0; gs(1,0) = 1; gs(3,1) = 1
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) = + 467.9220320217552853d0
    p(1,3) = - 161.2362660353179535d0
    p(2,3) = - 343.2253414521593413d0
    p(3,3) = - 258.5248172954496795d0
    p(0,4) = + 410.5097303987137138d0
    p(1,4) = +  76.9704063632674007d0
    p(2,4) = + 364.0379035243619228d0
    p(3,4) = + 173.407612560364754d0
    p(0,5) = + 121.5682375795313988d0
    p(1,5) = +  84.2658596720505955d0
    p(2,5) = -  20.8125620722026383d0
    p(3,5) = +  85.1172047350850534d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(2,0) = 4.720966256967035d-06
    matrix2P(4,1) = 1.510513880716169d-07

  elseif (cpr.eq.'u u~ -> A* A* g') then

    check = 'openLoops'
    legs = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0; gs(1,0) = 1; gs(3,1) = 1
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) = + 436.89797548439287d0
    p(1,3) = + 247.64869720720449d0
    p(2,3) = + 156.04337699053843d0
    p(3,3) = + 324.34615498156325d0
    p(0,4) = + 495.96159343298518d0
    p(1,4) = - 303.3770100848065d0
    p(2,4) = - 152.26618094022604d0
    p(3,4) = - 361.60102606943383d0
    p(0,5) = +  67.140431082621404d0
    p(1,5) = +  55.728312877601823d0
    p(2,5) = -   3.7771960503123614d0
    p(3,5) = +  37.254871087870626d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(2,0) = 1.195171665633129d-05
    matrix2P(4,1) = 4.157696725756211d-09

  elseif (cpr.eq.'u d~ -> W+ g g') then

    check = 'openLoops'
    legs = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0; gs(2,0) = 1; gs(4,1) = 1
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) = + 93.3071393593541671d0
    p(1,3) = + 6.3136012593360791d0
    p(2,3) = + 23.425206115008276d0
    p(3,3) = - 40.6647416324440627d0
    p(0,4) = + 437.8263733789333401d0
    p(1,4) = - 300.2349138749669919d0
    p(2,4) = + 191.1319006066975135d0
    p(3,4) = - 254.9892670037541507d0
    p(0,5) = + 468.8664872617123933d0
    p(1,5) = + 293.9213126156308817d0
    p(2,5) = - 214.5571067217057362d0
    p(3,5) = + 295.6540086361981707d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(4,0) = + 3.286801174851247d-04
    matrix2P(6,1) = - 2.090259157260240d-04

  elseif (cpr.eq.'u d~ -> W+ t t~') then

    check = 'openLoops'
    legs = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0; gs(2,0) = 1; gs(4,1) = 1
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) =   302.6058512244932217d0
    p(1,3) =   100.6087122807816741d0
    p(2,3) = - 272.6955006084110664d0
    p(3,3) = -  24.9269523895532146d0
    p(0,4) =   270.476875216737028d0
    p(1,4) =    13.4843528289578565d0
    p(2,4) = -  57.3257929008046219d0
    p(3,4) =   200.2639900833112279d0
    p(0,5) =   426.9172735587698071d0
    p(1,5) = - 114.0930651097395412d0
    p(2,5) =   330.0212935092156954d0
    p(3,5) = - 175.337037693757992d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(4,0) = + 5.636963741941358d-06
    matrix2P(6,1) = - 1.263456886647146d-06

  elseif (cpr.eq.'u u~ -> Z t t~') then

    check = 'openLoops'
    legs = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0; gs(2,0) = 1; gs(4,1) = 1
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) = + 373.51365847353088d0
    p(1,3) = + 205.31393289626283d0
    p(2,3) = + 129.36825347211888d0
    p(3,3) = + 268.90020197977543d0
    p(0,4) = + 445.70368746695988d0
    p(1,4) = - 251.51566631786329d0
    p(2,4) = - 126.23675718259466d0
    p(3,4) = - 299.78647026567d0
    p(0,5) = + 180.78265405950921d0
    p(1,5) = +  46.201733421600331d0
    p(2,5) = -   3.1314962895241925d0
    p(3,5) = +  30.886268285894587d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(4,0) = + 8.756605042513451d-06
    matrix2P(6,1) = - 4.741807687797462d-06

  elseif (cpr.eq.'d d~ -> Z t t~') then

    check = 'openLoops'
    legs = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0; gs(2,0) = 1; gs(4,1) = 1
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) = + 373.51365847353088d0
    p(1,3) = + 205.31393289626283d0
    p(2,3) = + 129.36825347211888d0
    p(3,3) = + 268.90020197977543d0
    p(0,4) = + 445.70368746695988d0
    p(1,4) = - 251.51566631786329d0
    p(2,4) = - 126.23675718259466d0
    p(3,4) = - 299.78647026567d0
    p(0,5) = + 180.78265405950921d0
    p(1,5) = +  46.201733421600331d0
    p(2,5) = -   3.1314962895241925d0
    p(3,5) = +  30.886268285894587d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(4,0) = + 1.228151570891972d-05
    matrix2P(6,1) = - 6.428320975896441d-06

  elseif (cpr.eq.'g g -> W+ b t~') then

    check = 'openLoops'
    legs = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0; gs(2,0) = 1; gs(4,1) = 1
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) = + 387.76679322266511d0
    p(1,3) = + 215.0230581105275d0
    p(2,3) = + 135.48597063818266d0
    p(3,3) = + 281.616269001305d0
    p(0,4) = + 430.62281258884883d0
    p(1,4) = - 263.40963309927599d0
    p(2,4) = - 132.20638849226233d0
    p(3,4) = - 313.96312323944784d0
    p(0,5) = + 181.610394188486d0
    p(1,5) = +  48.38657498874835d0
    p(2,5) = -   3.2795821459202981d0
    p(3,5) = +  32.346854238142846d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(4,0) = + 8.934467532519059d-06
    matrix2P(6,1) = - 4.793069478481384d-06

  elseif (cpr.eq.'g g -> Z t t~') then

    check = 'openLoops'
    legs = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0; gs(2,0) = 1; gs(4,1) = 1
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) = + 373.51365847353088d0
    p(1,3) = + 205.31393289626283d0
    p(2,3) = + 129.36825347211888d0
    p(3,3) = + 268.90020197977543d0
    p(0,4) = + 445.70368746695988d0
    p(1,4) = - 251.51566631786329d0
    p(2,4) = - 126.23675718259466d0
    p(3,4) = - 299.78647026567d0
    p(0,5) = + 180.78265405950921d0
    p(1,5) = +  46.201733421600331d0
    p(2,5) = -   3.1314962895241925d0
    p(3,5) = +  30.886268285894587d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(4,0) = + 4.854341911196551d-06
    matrix2P(6,1) = - 2.435075792140372d-06

  elseif (cpr.eq.'u u~ -> Z g g') then

    check = 'openLoops'
    legs = 5
    gs0 = 2
    gs1 = 4
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) =   442.2491972371382d0
    p(1,3) =   245.2952638366285d0
    p(2,3) =   154.56047926158473d0
    p(3,3) =   321.2638570597092d0
    p(0,4) =   491.24841473399675d0
    p(1,4) = - 300.49398430090025d0
    p(2,4) = - 150.81917832936608d0
    p(3,4) = - 358.16469092540393d0
    p(0,5) =    66.50238802886572d0
    p(1,5) =    55.1987204642716d0
    p(2,5) = -   3.7413009322186235d0
    p(3,5) =    36.90083386569479d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(gs0+gs0,0) = + 4.044494074381270d-4
    matrix2P(gs0+gs1,1) = - 6.871568616190759d-5

  elseif (cpr.eq.'g g -> g t t~') then

    check = 'openLoops'
    legs = 5
    gs0 = 3
    gs1 = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0; gs(gs0,0) = 1; gs(5,1) = 1
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) =   367.60373422843884d0
    p(1,3) =   208.37035412957272d0
    p(2,3) =   131.29410366285185d0
    p(3,3) =   272.90320496831384d0
    p(0,4) =   451.35675954947516d0
    p(1,4) = - 255.25987311474162d0
    p(2,4) = - 128.11599011936772d0
    p(3,4) = - 304.2492639994103d0
    p(0,5) =   181.03950622208615d0
    p(1,5) =    46.88951898516874d0
    p(2,5) = -   3.1781135434841024d0
    p(3,5) =    31.34605903109653d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(gs0+gs0,0) = + 1.497698555532992d-3
    matrix2P(gs0+gs1,1) = - 5.145734157852298d-4

  elseif (cpr.eq.'g g -> g g g') then

    check = 'openLoops'
    legs = 5
    gs0 = 3
    gs1 = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) =   436.8979754843929d0
    p(1,3) =   247.6486972072045d0
    p(2,3) =   156.04337699053843d0
    p(3,3) =   324.34615498156325d0
    p(0,4) =   495.9615934329852d0
    p(1,4) = - 303.3770100848065d0
    p(2,4) = - 152.26618094022604d0
    p(3,4) = - 361.60102606943383d0
    p(0,5) =    67.1404310826214d0
    p(1,5) =    55.72831287760182d0
    p(2,5) = -   3.7771960503123614d0
    p(3,5) =    37.254871087870626d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(gs0+gs0,0) = + 1.586205511886680d0
    matrix2P(gs0+gs1,1) = - 1.106655459389022d0

  elseif (cpr.eq.'d d~ -> d d~ g') then

    check = 'openLoops'
    legs = 5
    gs0 = 3
    gs1 = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) =   436.8979754843929d0
    p(1,3) =   247.6486972072045d0
    p(2,3) =   156.04337699053843d0
    p(3,3) =   324.34615498156325d0
    p(0,4) =   495.9615934329852d0
    p(1,4) = - 303.3770100848065d0
    p(2,4) = - 152.26618094022604d0
    p(3,4) = - 361.60102606943383d0
    p(0,5) =    67.1404310826214d0
    p(1,5) =    55.72831287760182d0
    p(2,5) = -   3.7771960503123614d0
    p(3,5) =    37.254871087870626d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(gs0+gs0,0) = + 1.230186974958428d+00
    matrix2P(gs0+gs1,1) = - 3.031687984375017d-01

  elseif (cpr.eq.'d d~ -> t t~ g') then

    check = 'openLoops'
    legs = 5
    gs0 = 3
    gs1 = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0; gs(gs0,0) = 1; gs(5,1) = 1
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) =   442.9730634176155d0
    p(1,3) =   231.3913515296482d0
    p(2,3) =   145.79962788530543d0
    p(3,3) =   303.0538662669989d0
    p(0,4) =   494.29406027133973d0
    p(1,4) = - 283.46127873151147d0
    p(2,4) = - 142.2703926866928d0
    p(3,4) = - 337.86307410576484d0
    p(0,5) =    62.73287631104489d0
    p(1,5) =    52.06992720186313d0
    p(2,5) = -   3.529235198612626d0
    p(3,5) =    34.80920783876596d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(gs0+gs0,0) = + 2.078818158542508d-03
    matrix2P(gs0+gs1,1) = - 7.053194817113668d-04

  elseif (cpr.eq.'b b~ -> t t~ g') then

    check = 'openLoops'
    legs = 5
    gs0 = 3
    gs1 = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
!     mass_b =   5.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0; gs(gs0,0) = 1; gs(5,1) = 1
    allocate(p(0:3,legs))
    if (mass_b.eq.5d0) then
      p(0,1) = + 500d0
      p(1,1) =   0d0
      p(2,1) =   0d0
      p(3,1) = + 499.9749993749688d0
      p(0,2) = + 500d0
      p(1,2) =   0d0
      p(2,2) =   0d0
      p(3,2) = - 499.9749993749688d0
      p(0,3) =   442.9730634176155d0
      p(1,3) =   231.3913515296482d0
      p(2,3) =   145.79962788530543d0
      p(3,3) =   303.0538662669989d0
      p(0,4) =   494.29406027133973d0
      p(1,4) = - 283.46127873151147d0
      p(2,4) = - 142.2703926866928d0
      p(3,4) = - 337.86307410576484d0
      p(0,5) =    62.73287631104489d0
      p(1,5) =    52.06992720186313d0
      p(2,5) = -   3.529235198612626d0
      p(3,5) =    34.80920783876596d0
      allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
      matrix2P(gs0+gs0,0) = + 2.078770892672683d-03
      matrix2P(gs0+gs1,1) = + 1.940582167713221d-03
   elseif (mass_b.eq.0d0) then
      p(0,1) = + 500d0
      p(1,1) =   0d0
      p(2,1) =   0d0
      p(3,1) = + 500d0
      p(0,2) = + 500d0
      p(1,2) =   0d0
      p(2,2) =   0d0
      p(3,2) = - 500d0
      p(0,3) =   442.9730634176155d0
      p(1,3) =   231.3913515296482d0
      p(2,3) =   145.79962788530543d0
      p(3,3) =   303.0538662669989d0
      p(0,4) =   494.29406027133973d0
      p(1,4) = - 283.46127873151147d0
      p(2,4) = - 142.2703926866928d0
      p(3,4) = - 337.86307410576484d0
      p(0,5) =    62.73287631104489d0
      p(1,5) =    52.06992720186313d0
      p(2,5) = -   3.529235198612626d0
      p(3,5) =    34.80920783876596d0
      allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
      matrix2P(gs0+gs0,0) = + 2.078818158542508d-03
      matrix2P(gs0+gs1,1) = - 7.053194817113668d-04
    endif

  elseif (cpr.eq.'u d~ -> W+ g g g') then

    check = 'openLoops'
    legs = 6
    gs0 = 3
    gs1 = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) =   263.4845819918209d0
    p(1,3) = -  20.526967070412418d0
    p(2,3) =   128.77196215622962d0
    p(3,3) =   214.37479126923805d0
    p(0,4) =   304.4026200608687d0
    p(1,4) = - 300.63757474928536d0
    p(2,4) =    25.834297426690753d0
    p(3,4) = -  40.13219188206375d0
    p(0,5) =    26.69963460475724d0
    p(1,5) =     6.683489195278001d0
    p(2,5) =     5.616637251194153d0
    p(3,5) =    25.23202025587033d0
    p(0,6) =   405.41316334255316d0
    p(1,6) =   314.4810526244197d0
    p(2,6) = - 160.22289683411447d0
    p(3,6) = - 199.47461964304466d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(gs0+gs0,0) = + 1.791371971976059d-06
    matrix2P(gs0+gs1,1) = - 1.917092204236691d-06

  elseif (cpr.eq.'u u~ -> Z g g g') then

    check = 'openLoops'
    legs = 6
    gs0 = 3
    gs1 = 5
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) =   266.1280911257227d0
    p(1,3) = -  20.453291457363225d0
    p(2,3) =   128.30977243171435d0
    p(3,3) =   213.60535494116894d0
    p(0,4) =   303.3100548723616d0
    p(1,4) = - 299.55852310228374d0
    p(2,4) =    25.741572685910604d0
    p(3,4) = -  39.98814898328692d0
    p(0,5) =    26.603803986383504d0
    p(1,5) =     6.659500743302574d0
    p(2,5) =     5.596477955798141d0
    p(3,5) =    25.141457214849993d0
    p(0,6) =   403.9580500155323d0
    p(1,6) =   313.35231381634435d0
    p(2,6) = - 159.64782307342307d0
    p(3,6) = - 198.75866317273204d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(gs0+gs0,0) = + 6.644843472959260d-07
    matrix2P(gs0+gs1,1) = - 6.988208366675967d-07

  elseif (cpr.eq.'u u~ -> W+ W- g g') then

    check = 'openLoops'
    legs = 6
    gs0 = 2
    gs1 = 4
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) =   260.9015189318375d0
    p(1,3) = -  20.304958210187117d0
    p(2,3) =   127.37923246317682d0
    p(3,3) =   212.0562264804183d0
    p(0,4) =   311.6591874642837d0
    p(1,4) = - 297.3860371459931d0
    p(2,4) =    25.554887277751206d0
    p(3,4) = -  39.69814324022019d0
    p(0,5) =    26.41086542484375d0
    p(1,5) =     6.611204097655833d0
    p(2,5) =     5.55589066207754d0
    p(3,5) =    24.959124019472025d0
    p(0,6) =   401.0284281790351d0
    p(1,6) =   311.0797912585243d0
    p(2,6) = - 158.49001040300553d0
    p(3,6) = - 197.31720725967017d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(gs0+gs0,0) = + 1.400987455762900d-06
    matrix2P(gs0+gs1,1) = - 1.075476167625012d-06

  elseif (cpr.eq.'u u~ -> Z Z g g') then

    check = 'openLoops'
    legs = 6
    gs0 = 2
    gs1 = 4
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) =   262.8342977029857d0
    p(1,3) = -  20.16623303040945d0
    p(2,3) =   126.50896684920824d0
    p(3,3) =   210.60744053183606d0
    p(0,4) =   312.6467098449736d0
    p(1,4) = - 295.3542707646298d0
    p(2,4) =    25.380294141675442d0
    p(3,4) = -  39.42692219160577d0
    p(0,5) =    26.23042417417883d0
    p(1,5) =     6.5660357960272195d0
    p(2,5) =     5.517932350469932d0
    p(3,5) =    24.788601188011228d0
    p(0,6) =   398.2885682778618d0
    p(1,6) =   308.95446799901197d0
    p(2,6) = - 157.40719334135358d0
    p(3,6) = - 195.96911952824155d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(gs0+gs0,0) = + 7.525155215596050d-08
    matrix2P(gs0+gs1,1) = - 5.736766868474369d-08

  elseif (cpr.eq.'d d~ -> t t~ b b~') then

    check = 'openLoops'
    legs = 6
    gs0 = 4
    gs1 = 6
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
!     mass_b =   5.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    if (mass_b.eq.5d0) then
      p(0,1) = + 500d0
      p(1,1) =   0d0
      p(2,1) =   0d0
      p(3,1) = + 500d0
      p(0,2) = + 500d0
      p(1,2) =   0d0
      p(2,2) =   0d0
      p(3,2) = - 500d0
      p(0,3) =   284.5880883702362d0
      p(1,3) = -  18.5481459059631d0
      p(2,3) =   116.35820988448233d0
      p(3,3) =   193.70883644830258d0
      p(0,4) =   324.4084490179861d0
      p(1,4) = - 271.65579708569174d0
      p(2,4) =    23.343844046936432d0
      p(3,4) = -  36.263406474090026d0
      p(0,5) =    24.638433047977387d0
      p(1,5) =     6.03919382389566d0
      p(2,5) =     5.075187526664823d0
      p(3,5) =    22.799627027350045d0
      p(0,6) =   366.3650295638003d0
      p(1,6) =   284.16474916775917d0
      p(2,6) = - 144.77724145808355d0
      p(3,6) = - 180.24505700156263d0
      allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
      matrix2P(gs0+gs0,0) = + 3.182318443234540d-09
      matrix2P(gs0+gs1,1) = + 1.897612313201693d-09
   elseif (mass_b.eq.0d0) then
      p(0,1) = + 500d0
      p(1,1) =   0d0
      p(2,1) =   0d0
      p(3,1) = + 500d0
      p(0,2) = + 500d0
      p(1,2) =   0d0
      p(2,2) =   0d0
      p(3,2) = - 500d0
      p(0,3) =   284.7108975623757d0
      p(1,3) = -  18.560754805368916d0
      p(2,3) =   116.4373093788961d0
      p(3,3) =   193.84051835576534d0
      p(0,4) =   324.56700192222627d0
      p(1,4) = - 271.8404667899213d0
      p(2,4) =    23.359713028279646d0
      p(3,4) = -  36.28805808329464d0
      p(0,5) =    24.142162336595277d0
      p(1,5) =     6.043299225470772d0
      p(2,5) =     5.078637603525002d0
      p(3,5) =    22.815126053782954d0
      p(0,6) =   366.5799381788029d0
      p(1,6) =   284.3579223698194d0
      p(2,6) = - 144.8756600107007d0
      p(3,6) = - 180.36758632625367d0
      allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
      matrix2P(gs0+gs0,0) = + 3.227591631303408d-09
      matrix2P(gs0+gs1,1) = - 2.401278721882181d-09
    endif

  elseif (cpr.eq.'g g -> t t~ b b~') then

    check = 'openLoops'
    legs = 6
    gs0 = 4
    gs1 = 6
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
!     mass_b =   5.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    if (mass_b.eq.5d0) then
      p(0,1) = + 500d0
      p(1,1) =   0d0
      p(2,1) =   0d0
      p(3,1) = + 500d0
      p(0,2) = + 500d0
      p(1,2) =   0d0
      p(2,2) =   0d0
      p(3,2) = - 500d0
      p(0,3) =   284.5880883702362d0
      p(1,3) = -  18.5481459059631d0
      p(2,3) =   116.35820988448233d0
      p(3,3) =   193.70883644830258d0
      p(0,4) =   324.4084490179861d0
      p(1,4) = - 271.65579708569174d0
      p(2,4) =    23.343844046936432d0
      p(3,4) = -  36.263406474090026d0
      p(0,5) =    24.638433047977387d0
      p(1,5) =     6.03919382389566d0
      p(2,5) =     5.075187526664823d0
      p(3,5) =    22.799627027350045d0
      p(0,6) =   366.3650295638003d0
      p(1,6) =   284.16474916775917d0
      p(2,6) = - 144.77724145808355d0
      p(3,6) = - 180.24505700156263d0
      allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
      matrix2P(gs0+gs0,0) = + 5.715140567912768d-08
      matrix2P(gs0+gs1,1) = + 3.117569091842899d-08
   elseif (mass_b.eq.0d0) then
      p(0,1) = + 500d0
      p(1,1) =   0d0
      p(2,1) =   0d0
      p(3,1) = + 500d0
      p(0,2) = + 500d0
      p(1,2) =   0d0
      p(2,2) =   0d0
      p(3,2) = - 500d0
      p(0,3) =   284.7108975623757d0
      p(1,3) = -  18.560754805368916d0
      p(2,3) =   116.4373093788961d0
      p(3,3) =   193.84051835576534d0
      p(0,4) =   324.56700192222627d0
      p(1,4) = - 271.8404667899213d0
      p(2,4) =    23.359713028279646d0
      p(3,4) = -  36.28805808329464d0
      p(0,5) =    24.142162336595277d0
      p(1,5) =     6.043299225470772d0
      p(2,5) =     5.078637603525002d0
      p(3,5) =    22.815126053782954d0
      p(0,6) =   366.5799381788029d0
      p(1,6) =   284.3579223698194d0
      p(2,6) = - 144.8756600107007d0
      p(3,6) = - 180.36758632625367d0
      allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
      matrix2P(gs0+gs0,0) = + 7.215798472013058d-08
      matrix2P(gs0+gs1,1) = - 5.103418697049000d-08
    endif

  elseif (cpr.eq.'u u~ -> u u~ u u~') then

    check = 'openLoops'
    legs = 6
    gs0 = 4
    gs1 = 6
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) =   254.11173663453982d0
    p(1,3) = -  20.78819186394802d0
    p(2,3) =   130.41070543048184d0
    p(3,3) =   217.10291035261022d0
    p(0,4) =   308.2764271997574d0
    p(1,4) = - 304.46346817638226d0
    p(2,4) =    26.163062947104724d0
    p(3,4) = -  40.642911439539645d0
    p(0,5) =    27.039412347527474d0
    p(1,5) =     6.768542826393844d0
    p(2,5) =     5.688114196681263d0
    p(3,5) =    25.55312123777459d0
    p(0,6) =   410.57242381817514d0
    p(1,6) =   318.4831172139364d0
    p(2,6) = - 162.2618825742678d0
    p(3,6) = - 202.0131201508452d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(gs0+gs0,0) = + 8.258218711346499d-07
    matrix2P(gs0+gs1,1) = - 3.861723516505142d-07

  elseif (cpr.eq.'g g -> u u~ u u~') then

    check = 'openLoops'
    legs = 6
    gs0 = 4
    gs1 = 6
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) =   254.11173663453982d0
    p(1,3) = -  20.78819186394802d0
    p(2,3) =   130.41070543048184d0
    p(3,3) =   217.10291035261022d0
    p(0,4) =   308.2764271997574d0
    p(1,4) = - 304.46346817638226d0
    p(2,4) =    26.163062947104724d0
    p(3,4) = -  40.642911439539645d0
    p(0,5) =    27.039412347527474d0
    p(1,5) =     6.768542826393844d0
    p(2,5) =     5.688114196681263d0
    p(3,5) =    25.55312123777459d0
    p(0,6) =   410.57242381817514d0
    p(1,6) =   318.4831172139364d0
    p(2,6) = - 162.2618825742678d0
    p(3,6) = - 202.0131201508452d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(gs0+gs0,0) = + 3.711499988284831d-08
    matrix2P(gs0+gs1,1) = - 2.885414903216092d-08

  elseif (cpr.eq.'g g -> u u~ d d~') then

    check = 'openLoops'
    legs = 6
    gs0 = 4
    gs1 = 6
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) =   254.11173663453982d0
    p(1,3) = -  20.78819186394802d0
    p(2,3) =   130.41070543048184d0
    p(3,3) =   217.10291035261022d0
    p(0,4) =   308.2764271997574d0
    p(1,4) = - 304.46346817638226d0
    p(2,4) =    26.163062947104724d0
    p(3,4) = -  40.642911439539645d0
    p(0,5) =    27.039412347527474d0
    p(1,5) =     6.768542826393844d0
    p(2,5) =     5.688114196681263d0
    p(3,5) =    25.55312123777459d0
    p(0,6) =   410.57242381817514d0
    p(1,6) =   318.4831172139364d0
    p(2,6) = - 162.2618825742678d0
    p(3,6) = - 202.0131201508452d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(gs0+gs0,0) = + 8.205424719121766d-08
    matrix2P(gs0+gs1,1) = - 7.420481476358859d-08

  elseif (cpr.eq.'g g -> u u~ g g') then

    check = 'openLoops'
    legs = 6
    gs0 = 4
    gs1 = 6
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) =   254.11173663453982d0
    p(1,3) = -  20.78819186394802d0
    p(2,3) =   130.41070543048184d0
    p(3,3) =   217.10291035261022d0
    p(0,4) =   308.2764271997574d0
    p(1,4) = - 304.46346817638226d0
    p(2,4) =    26.163062947104724d0
    p(3,4) = -  40.642911439539645d0
    p(0,5) =    27.039412347527474d0
    p(1,5) =     6.768542826393844d0
    p(2,5) =     5.688114196681263d0
    p(3,5) =    25.55312123777459d0
    p(0,6) =   410.57242381817514d0
    p(1,6) =   318.4831172139364d0
    p(2,6) = - 162.2618825742678d0
    p(3,6) = - 202.0131201508452d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(gs0+gs0,0) = + 2.124095665038298d-04
    matrix2P(gs0+gs1,1) = - 2.423535768116181d-04

  elseif (cpr.eq.'g g -> t t~ g g') then

    check = 'openLoops'
    legs = 6
    gs0 = 4
    gs1 = 6
    ew_reno_scheme = 3; alZ = 0.0078125d0
    als = 0.125808685692d0; Qren = 100d0
     mass_z =  91.1876d0
    width_z =   0.d0
     mass_w =  80.399d0
    width_w =   0.d0
     mass_h = 125.d0
     mass_t = 172.d0
    complex_mass_scheme = 0
    deltaIR2 = pi**2/6d0
    muUV = 100d0
    muIR = 100d0
    allocate (gs(0:legs,0:1)); gs = 0
    allocate(p(0:3,legs))
    p(0,1) = + 500d0
    p(1,1) =   0d0
    p(2,1) =   0d0
    p(3,1) = + 500d0
    p(0,2) = + 500d0
    p(1,2) =   0d0
    p(2,2) =   0d0
    p(3,2) = - 500d0
    p(0,3) =   284.7108975623757d0
    p(1,3) = -  18.560754805368916d0
    p(2,3) =   116.4373093788961d0
    p(3,3) =   193.84051835576534d0
    p(0,4) =   324.56700192222627d0
    p(1,4) = - 271.8404667899213d0
    p(2,4) =    23.359713028279646d0
    p(3,4) = -  36.28805808329464d0
    p(0,5) =    24.142162336595277d0
    p(1,5) =     6.043299225470772d0
    p(2,5) =     5.078637603525002d0
    p(3,5) =    22.815126053782954d0
    p(0,6) =   366.5799381788029d0
    p(1,6) =   284.3579223698194d0
    p(2,6) = - 144.8756600107007d0
    p(3,6) = - 180.36758632625367d0
    allocate (matrix2P(0:2*legs-2,0:1)); matrix2P = 0d0
    matrix2P(gs0+gs0,0) = + 9.274364276562245d-05
    matrix2P(gs0+gs1,1) = - 9.746105756314310d-05

  else

    if (warnings(381).le.warning_limit) then
      warnings(381) = warnings(381) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 381: No checked process'
      write(nx,*)
      call toomanywarnings(381)
    endif
    call istop (ifail,1)

  endif

  prTot = 0
  call define_process_rcl (1,cpr,order)

  call set_gs_power_rcl (1,gs(0:legs,0:1))

  if (gs0.ne.111) then
    call unselect_all_gs_powers_BornAmpl_rcl (1)
    call select_gs_power_BornAmpl_rcl (1,gs0)
  endif
  if (gs1.ne.111) then
    call unselect_all_gs_powers_LoopAmpl_rcl (1)
    call select_gs_power_LoopAmpl_rcl (1,gs1)
  endif

  call openOutput

  write(nx,*)
  write(nx,'(1x,75("*"))')
  write(nx,*)

  call generate_processes_rcl

  call compute_process_rcl (1,p,'NLO')

  write(nx,*)
  write(nx,'(1x,75("*"))')
  write(nx,*)

  write(nx,'(1x,75("-"))')

  write(nx,'(2x,a)') trim(cpr)

  write(nx,'(1x,75("-"))')

  do i = 1,legs
    fmt(0) = 'f15.10'
    if (abs(p(0,i)).ge.1d4) fmt(0) = 'e15.9'
    fmtTot = '(2x,a,i1,1x,a,1x,"(",'//fmt(0)
    do n = 1,3
      fmt(n) = 'f16.10'
      if (abs(p(n,i)).ge.1d4) fmt(n) = 'e16.9'
      fmtTot = trim(fmtTot)//',",",'//fmt(n)
    enddo
    fmtTot = trim(fmtTot)//',")")'
    write(nx,trim(fmtTot)) 'p',i,'=',p(0:3,i)
  enddo

  fmt1i = '(4x,i2," | ",e21.14,2x,a)'
  fmt2i = '(4x,i2," | ",e21.14,2x,a,5x,a,e21.14)'
  fmt1 = '(3x,a3," | ",e21.14,2x,a)'
  fmt2 = '(3x,a3," | ",e21.14,2x,a,5x,a,e21.14)'

  write(nx,'(1x,75("-"))')
  write(nx,*) '  als |        | A0 |^2      '
  do i = 0,2*legs-4,2
    if (matrix2P(i,0).ne.0d0) then
      mat2r = matrix2(i,0,1)
      mat2c = matrix2P(i,0)
      del = abs(mat2r-mat2c)/abs(mat2r+mat2c)*2
      write(nx,*) ' ========================================'
      write(nx,fmt1i) i/2,mat2r,'recola'
      write(nx,fmt2i) i/2,mat2c,check,'Delta =',del
    endif
  enddo
  if (matrix2Psum(0).ne.0d0) then
    mat2r = sum(matrix2(0:2*legs-4,0,1))
    mat2c =     matrix2Psum(0)
    del = abs(mat2r-mat2c)/abs(mat2r+mat2c)*2
    write(nx,*) ' ========================================'
    write(nx,fmt1) 'all',mat2r,'recola'
    write(nx,fmt2) 'all',mat2c, check,'Delta =',del
  endif
  write(nx,'(1x,75("-"))')

  if (loop(1)) then
    write(nx,*) '  als |   2*Re{ A1 * A0^* }  '
    do i = 0,2*legs-2,2
      if (matrix2P(i,1).ne.0d0) then
        mat2r = matrix2(i,4,1)
        mat2c = matrix2P(i,1)
        del = abs(mat2r-mat2c)/abs(mat2r+mat2c)*2
        write(nx,*) ' ========================================'
        write(nx,fmt1i) i/2,mat2r,'recola'
        write(nx,fmt2i) i/2,mat2c,check,'Delta =',del
      endif
    enddo
    if (matrix2Psum(1).ne.0d0) then
      mat2r = sum(matrix2(0:2*legs-2,4,1))
      mat2c =     matrix2Psum(1)
      del = abs(mat2r-mat2c)/abs(mat2r+mat2c)*2
      write(nx,*) ' ========================================'
      write(nx,fmt1) 'all',mat2r,'recola'
      write(nx,fmt2) 'all',mat2c,check,'Delta =',del
    endif
    write(nx,'(1x,75("-"))')
  endif

  if (present(delta)) delta = del

  deallocate (matrix2P,p,gs)

  call reset_recola_rcl

  end subroutine check_process_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine check_all_processes_rcl

  ! This routine runs recola for all checked processes

  call set_output_file_rcl ('output.rcl')

  ! 4-legs processes checked with Pole
  call check_process_rcl ('u~ u -> nu_e~ nu_e')  ! Delta = 0.11398542008638E-14
  call check_process_rcl ('u d~ -> nu_e e+')     ! Delta = 0.67768053335133E-15
  call check_process_rcl ('d~ d -> u~ u')        ! Delta = 0.74602576650244E-15
  call check_process_rcl ('e+ e- -> nu_e~ nu_e') ! Delta = 0.17685522062960E-14
  call check_process_rcl ('e+ e- -> W+ W-')      ! Delta = 0.56321838985838E-12
  call check_process_rcl ('u d~ -> W+ H')        ! Delta = 0.13222261973147E-14
  call check_process_rcl ('e+ e- -> Z H')        ! Delta = 0.26492623303543E-15

  ! 4-legs processes checked with openLoops
  call check_process_rcl ('u d~ -> W+ g')        ! Delta = 0.11081338630742E-13
  call check_process_rcl ('g g -> g g')          ! Delta = 0.56909352655424E-14
  call check_process_rcl ('b b~ -> t t~')        ! Delta = 0.94572495398593E-11 (mb=0)
                                                ! Delta = 0.94624111514847E-11 (mb=5)
  call check_process_rcl ('g g -> b b~')         ! Delta = 0.94628018748532E-11 (mb=0)
                                                ! Delta = 0.94625201927119E-11 (mb=5)
  ! 5-legs processes checked with Pole
  call check_process_rcl ('u d~ -> e+ nu_e A*')   ! Delta = 0.56249112075080E-14

  ! 5-legs processes checked with openLoops
  call check_process_rcl ('u u~ -> W+ W- g')     ! Delta = 0.99313950048747E-08
  call check_process_rcl ('u u~ -> Z Z g')       ! Delta = 0.27026273673635E-11
  call check_process_rcl ('u u~ -> Z A* g')       ! Delta = 0.76263106415419E-11
  call check_process_rcl ('u u~ -> A* A* g')       ! Delta = 0.10777473344478E-08
  call check_process_rcl ('u d~ -> W+ g g')      ! Delta = 0.10592483036673E-10
  call check_process_rcl ('u d~ -> W+ t t~')     ! Delta = 0.76724965191442E-11
  call check_process_rcl ('u u~ -> Z t t~')      ! Delta = 0.15650738401991E-10
  call check_process_rcl ('d d~ -> Z t t~')      ! Delta = 0.15576036895349E-10
  call check_process_rcl ('g g -> W+ b t~')      ! Delta = 0.15021760468786E-10 (mb=0)
  call check_process_rcl ('g g -> Z t t~')       ! Delta = 0.46821028180503E-10
  call check_process_rcl ('u u~ -> Z g g')       ! Delta = 0.15726020733631E-10
  call check_process_rcl ('g g -> g t t~')       ! Delta = 0.75346041009646E-11
  call check_process_rcl ('g g -> g g g')        ! Delta = 0.12699207342060E-10
  call check_process_rcl ('d d~ -> d d~ g')      ! Delta = 0.12568747092807E-10
  call check_process_rcl ('d d~ -> t t~ g')      ! Delta = 0.11924664738735E-10
  call check_process_rcl ('b b~ -> t t~ g')      ! Delta = 0.11924664738735E-10 (mb=0)
                                                 ! Delta = 0.11737157858432E-10 (mb=5)

  ! 6-legs processes checked with openLoops
  call check_process_rcl ('u d~ -> W+ g g g')    ! Delta = 0.43619656902563E-10
  call check_process_rcl ('u u~ -> Z g g g')     ! Delta = 0.87413423988204E-10
  call check_process_rcl ('u u~ -> W+ W- g g')   ! Delta = 0.17299210026665E-09
  call check_process_rcl ('u u~ -> Z Z g g')     ! Delta = 0.29066251823306E-09
  call check_process_rcl ('d d~ -> t t~ b b~')   ! Delta = 0.20235687442134E-09 (mb=0)
                                                 ! Delta = 0.73955595758043E-10 (mb=5)
  call check_process_rcl ('g g -> t t~ b b~')    ! Delta = 0.58158315425266E-09 (mb=0)
                                                 ! Delta = 0.20179121872741E-07 (mb=5)
  call check_process_rcl ('u u~ -> u u~ u u~')   ! Delta = 0.18829654552338E-09
  call check_process_rcl ('g g -> u u~ u u~')    ! Delta = 0.15243471860270E-09
  call check_process_rcl ('g g -> u u~ d d~')    ! Delta = 0.11050415933361E-10
  call check_process_rcl ('g g -> u u~ g g')     ! Delta = 0.20493088251540E-07
  call check_process_rcl ('g g -> t t~ g g')     ! Delta = 0.59855375919497E-08

  end subroutine check_all_processes_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module check_rcl

!#####################################################################



