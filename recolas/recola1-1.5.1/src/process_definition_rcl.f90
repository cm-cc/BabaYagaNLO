!#####################################################################
!!
!!  File  process_definition_rcl.f90
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

  module process_definition_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine define_process_rcl (npr,processIn,order)

  ! With this subroutine the user can define the processes to be
  ! computed.

  ! "npr" is the process number the user wants to assign to process
  ! described by "processIn". If the subroutine is called more than
  ! once, each process "processIn" must have a different process
  ! number "npr".

  ! "order" is the loop-order at which we want to compute the process
  ! "processIn". It is a character variable which can take only the
  ! values 'LO' and 'NLO'.

  ! "processIn" is a character variable which defines the process.

  ! "processIn" is a list of particles separated by "->". The
  ! incoming (outgoing) particles are on the left (right) side of
  ! "->".
  ! Any number of incoming and outgoing particles is allowed.
  ! The symbols for the particles are:
  ! Scalars:       'H', 'p0', 'p+', 'p-'
  ! Vector bosons: 'g', 'A', 'Z', 'W+', 'W-'
  ! leptons:       'nu_e', 'nu_e~', 'e-', 'e+',
  !                'nu_mu', 'nu_mu~', 'mu-', 'mu+',
  !                'nu_tau', 'nu_tau~', 'tau-', 'tau+'
  ! quarks:        'u', 'u~', 'd', 'd~',
  !                'c', 'c~', 's', 's~',
  !                't', 't~', 'b', 'b~'
  ! All these symbols (particles and "->") must be separated by at
  ! least one blank character.
  ! Examples of allowed values for "processIn":
  ! 'e+ e- -> mu+ mu-'
  ! 'u u~ -> W+ W- g'
  ! 'u d~ -> W+ g g g'
  ! 'u  g  ->  u  g  Z'
  ! 'u    g  -> u        g  tau- tau+'

  ! Additional symbols for helicities are allowed.
  ! The symbols "[-]", "[+]" and "[0]" account for a specific
  ! polarization of a particle.
  ! For fermions only "[-]" and "[+]" are allowed and they represent
  ! left-handed and right-handed fermions respectively.
  ! For massless vector bosons only "[-]" and "[+]" are allowed and
  ! they represent respectively the -1 and +1 transverse
  ! polarizations.
  ! For massive vector bosons "[-]", "[+]" and "[0]" are allowed and
  ! they represent respectively -1, +1 transverse polarizations and
  ! longitudinal polarization.
  ! The symbols "[-]", "[+]", "[0]" must follow the particle in the
  ! character "processIn".
  ! Blank characters can separate "[-],[+],[0]" from the particle.
  ! Examples:
  ! 'e+[+] e- [-] -> Z H':
  ! e+ is right-handed, e- is left-handed, Z is unpolarized.
  ! 'u u~ -> W+[-] W-[+]':
  ! u and u~ are unpolarised, while W+ and W- are transverse

  ! Contributions with specific intermediate states can be selected,
  ! where intermediate states are particles decaying into any number
  ! of other particles.  To this end, in the process declaration the
  ! decaying particle must be followed by round brackets "( ... )"
  ! containing the decay products.
  ! Multiple and nested decays are allowed.
  ! Blank characters can separate "(", ")" and the particles.
  ! Examples:
  ! 'e+ e- -> W+ W-(e- nu_e~)'
  ! 'e+ e- -> Z H ( b~[+] b[-] )'
  ! 'e+ e- -> t( W+(u d~) b) t~(e- nu_e~ b~)'
  ! 'u  u~ -> Z ( mu+(e+ nu_e nu_mu~) mu-(e- nu_e~ nu_mu) ) H'


  integer,          intent(in) :: npr
  character(len=*), intent(in) :: processIn,order

  integer                :: i,j,n,n0,n1,n2,le,leIn,l0,lmax,sign,nd0, &
                            nd,noff,lend(99),pond(99),binnd(99),     &
                            polegs(99),pol(99)
  integer, allocatable   :: inprT(:),legsInT(:),legsOutT(:),         &
                            parT(:,:),helT(:,:),resMaxT(:),          &
                            binResT(:,:),parResT(:,:),powgsT(:,:,:), &
                            refschemeT(:),NalphaT(:,:),              &
                            NLOschemeT(:),NoffPhT(:),                &
                            pa(:),paT(:),he(:),heT(:)
  logical                :: opennd(99), he_aux
  logical, allocatable   :: loopT(:),prexistsT(:)
  integer, allocatable   :: qflowT(:,:),polprojT(:,:)
  integer, parameter     :: nopol=-2
  character              :: cpr*99,che*3
  character, allocatable :: cpa(:)*7,cpaT(:)*7,processT(:)*99

  if (processes_generated) then
    if (warnings(401).le.warning_limit) then
      warnings(401) = warnings(401) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 401: Processes already generated.'
      write(nx,*) '           define_process_rcl can not be called ', &
                         'at this point'
      write(nx,*)
      call toomanywarnings(401)
    endif
  endif

  do i = 1,prTot
    if (inpr(i).eq.npr) then
      if (warnings(402).le.warning_limit) then
        warnings(402) = warnings(402) + 1
        call openOutput
        write(nx,*)
        write(nx,'(1x,a,i3,a)') 'ERROR 402: Index ',npr,' has been already used'
        write(nx,*)             '           for the definition of another process.'
        write(nx,*)             '           Definition of process ', &
                                            trim(adjustl(processIn)), &
                                          ' has failed.'
        write(nx,*)
        call toomanywarnings(402)
      endif
      call istop (ifail,1)
      return
    endif
  enddo

  cpr = adjustl(processIn)

  noff = 0
  n = 99
  do while (scan(cpr,'*').ne.0)
    n = scan(cpr,'*')
    noff = noff + 1
    cpr = cpr(1:n-1)//cpr(n+1:)
  enddo

  n = 99
  do while (index (cpr,' [').ne.0)
    n = index (cpr,' [')
    cpr = cpr(1:n-1)//cpr(n+1:)
  enddo
  n = 99
  do while (scan(cpr,'(').ne.0)
    n = scan(cpr,'(')
    cpr = cpr(1:n-1)//' { '//cpr(n+1:)
  enddo
  n = 99
  do while (scan(cpr,')').ne.0)
    n = scan(cpr,')')
    cpr = cpr(1:n-1)//' } '//cpr(n+1:)
  enddo


  le = 0
  leIn = 0

  pond = 0
  nd0 = 0
  nd = 0

  do while (cpr.ne.'')

    if (le.gt.0) then
      allocate (cpaT(le)); cpaT(:le) = cpa(:le); deallocate (cpa)
      allocate ( paT(le));  paT(:le) =  pa(:le); deallocate ( pa)
      allocate ( heT(le));  heT(:le) =  he(:le); deallocate ( he)
    endif
    le = le + 1
    allocate (cpa(le),pa(le),he(le))
    if (le.gt.1) then
      cpa(:le-1) = cpaT(:le-1); deallocate (cpaT)
       pa(:le-1) =  paT(:le-1); deallocate ( paT)
       he(:le-1) =  heT(:le-1); deallocate ( heT)
    endif

    n0 = scan(cpr,'>')
    n1 = scan(cpr,' ') - 1
    sign = 1
    if (n1.ge.n0) then
      cpr = adjustl(cpr(n0+1:))
      n1 = scan(cpr,' ') - 1
      sign = -1
    else
      leIn = le
    endif

    n2 = scan(cpr(1:n1),'[')
    he_aux = .false.
    if (n2.eq.0) then
      cpa(le) = cpr(1:n1);        he(le) = 111
    else
      cpa(le) = cpr(1:n2-1)
      che     = cpr(n2:n2+2)
      if     (che.eq.'[-]') then; he(le) = -1 * sign
      elseif (che.eq.'[0]') then; he(le) =  0
      elseif (che.eq.'[+]') then; he(le) = +1 * sign
      elseif (che.eq.'[2]') then; he(le) = 2; he_aux = .true.
      elseif (che.eq.'[3]') then; he(le) = 3; he_aux = .true.
      elseif (che.eq.'[T]') then; he(le) = 3; he_aux = .true.
      elseif (che.eq.'[4]') then; he(le) = 4; he_aux = .true.
      else
        he(le) = 111
        if (warnings(397).le.warning_limit) then
          warnings(397) = warnings(397) + 1
          call openOutput
          write(nx,*)
          write(nx,*) 'ERROR 397 (define_process_rcl): '
          write(nx,'(a,i0,2a)') ' In process ',npr, &
                                ' particle '//trim(cpa(le)), &
                                ' has wrong helicity '//trim(che)
          write(nx,*)
          call toomanywarnings(397)
        endif
        call istop (ifail,1)
      endif
    endif

    pa(le) = npar(cpa(le))
    if (sign.eq.-1) pa(le) = anti(pa(le))

    cpr = adjustl(cpr(n1+1:))

    if (cpr(1:1).eq.'{') then
      nd0 = nd0 + 1
      nd  = nd0
      opennd(nd) = .true.
      pond(nd) = pa(le)
      lend(nd) = le
      if (he(le).ne.111) then
        ! todo: info
        pol(nd) = he(le)
        he(le) = 111
        select case (pa(le))
        case (1:16,20:43)
          if (warnings(398).le.warning_limit) then
            warnings(398) = warnings(398) + 1
            call openOutput
            write(nx,*)
            write(nx,*) 'ERROR 398 (define_process_rcl): '
            write(nx,*) ' Polarization of resonant particles implemented', &
                        ' only for massive vector bosons'
            write(nx,*)
            call toomanywarnings(398)
          endif
          call istop (ifail,1)
        end select
      else
        pol(nd) = nopol
      endif
      le = le - 1
      cpr = adjustl(cpr(2:))
    else
      if (he_aux) then
        if (warnings(399).le.warning_limit) then
          warnings(399) = warnings(399) + 1
          call openOutput
          write(nx,*)
          write(nx,*) 'ERROR 399 (define_process_rcl): '
          write(nx,'(a,i0,a)') ' In process ',npr, &
                               ' particle '//cpa(le)//' has wrong helicity '//che
          write(nx,*)
          call toomanywarnings(399)
        endif
        call istop (ifail,1)
      endif
    endif

    do while (cpr(1:1).eq.'}')
      polegs(nd) = le - lend(nd) + 1
      cpr = adjustl(cpr(2:))
      binnd(nd) = 2**(lend(nd)-1)*( 2**polegs(nd) - 1 )
      opennd(nd) = .false.
      nd = nd - 1
      do i = nd0,1,-1
        if (opennd(i)) then
          nd = i
          exit
        endif
      enddo
    enddo

  enddo

  l0 = 0
  if (prTot.gt.0) then
    l0 = maxval(legsIn+legsOut)
    allocate (  processT(            prTot));   processT =  process
    allocate (     inprT(            prTot));      inprT =     inpr
    allocate (   legsInT(            prTot));    legsInT =   legsIn
    allocate (  legsOutT(            prTot));   legsOutT =  legsOut
    allocate (      parT(l0,         prTot));       parT =      par
    allocate (      helT(l0,         prTot));       helT =      hel
    allocate (   resMaxT(            prTot));    resMaxT =   resMax
    allocate (   binResT(l0,         prTot));    binResT =   binRes
    allocate (   parResT(l0,         prTot));    parResT =   parRes
    allocate (     loopT(            prTot));      loopT =     loop
    allocate (    qflowT(l0,         prTot));     qflowT =    qflow
    allocate (polprojT(l0,           prTot));  polprojT = polproj
    allocate (    powgsT(0:l0,0:1,   prTot));     powgsT =    powgs
    allocate ( prexistsT(            prTot));  prexistsT = prexists
    allocate (refschemeT(            prTot)); refschemeT = refscheme
    allocate (   NalphaT(1:4,        prTot));    NalphaT =   Nalpha
    allocate (NLOschemeT(            prTot)); NLOschemeT = NLOscheme
    allocate (   NoffPhT(            prTot));    NoffPhT =   NoffPh
    deallocate (process,inpr,legsIn,legsOut,par,hel,resMax,binRes, &
                parRes,loop,qflow,polproj,powgs,prexists,          &
                refscheme,Nalpha,NLOscheme,NoffPh)
  endif

  prTot = prTot + 1

  lmax = max(l0,le)

  allocate ( process(              prTot))
  allocate (    inpr(              prTot));      inpr = 0
  allocate (  legsIn(              prTot));    legsIn = 0
  allocate ( legsOut(              prTot));   legsOut = 0
  allocate (     par(lmax,         prTot));       par = 0
  allocate (     hel(lmax,         prTot));       hel = 0
  allocate (  resMax(              prTot));    resMax = 0
  allocate (  binRes(lmax,         prTot));    binRes = 0
  allocate (  parRes(lmax,         prTot));    parRes = 0
  allocate (    loop(              prTot));      loop = .false.
  allocate (   qflow(lmax,         prTot));     qflow = 0
  allocate ( polproj(lmax,         prTot));   polproj = nopol
  allocate (   powgs(0:lmax,0:1,   prTot));     powgs = 1
  allocate ( prexists(             prTot));  prexists = .true.
  allocate (refscheme(             prTot)); refscheme = 0
  allocate (   Nalpha(1:4,         prTot));    Nalpha = 0
  allocate (NLOscheme(             prTot)); NLOscheme = 0
  allocate (   NoffPh(             prTot));    NoffPh = 0

  if (prTot.gt.1) then
      process(             :prTot-1) =   processT
         inpr(             :prTot-1) =      inprT
       legsIn(             :prTot-1) =    legsInT
      legsOut(             :prTot-1) =   legsOutT
          par(:l0,         :prTot-1) =       parT
          hel(:l0,         :prTot-1) =       helT
       resMax(             :prTot-1) =    resMaxT
       binRes(:l0,         :prTot-1) =    binResT
       parRes(:l0,         :prTot-1) =    parResT
         loop(             :prTot-1) =      loopT
        qflow(:l0,         :prTot-1) =     qflowT
      polproj(:l0,         :prTot-1) =   polprojT
        powgs(0:l0,0:1,    :prTot-1) =     powgsT
     prexists(             :prTot-1) =  prexistsT
    refscheme(             :prTot-1) = refschemeT
       Nalpha(1:4,         :prTot-1) =    NalphaT
    NLOscheme(             :prTot-1) = NLOschemeT
       NoffPh(             :prTot-1) =    NoffPhT
    deallocate (processT,inprT,legsInT,legsOutT,parT,helT,resMaxT, &
                binResT,parResT,loopT,qflowT,polprojT,powgsT,      &
                prexistsT,refschemeT,NalphaT,NLOschemeT,NoffPhT)
  endif

  process(prTot) = adjustl(processIn)

  inpr(prTot) = npr
  do i = 1,prTot
  do j = 1,prTot
    if (i.ne.j.and.inpr(i).eq.inpr(j)) then
      if (warnings(404).le.warning_limit) then
        warnings(404) = warnings(404) + 1
        call openOutput
        write(nx,*)
        write(nx,*) 'ERROR 404: define_process_rcl called for different '// &
                           'processes with the same process-index'
        write(nx,*)
        call toomanywarnings(404)
      endif
      call istop (ifail,1)
    endif
  enddo
  enddo

  if (leIn.le.0) then
    if (warnings(405).le.warning_limit) then
      warnings(405) = warnings(405) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 405: define_process_rcl called for a process ', &
                         'without incoming particles'
      write(nx,*)
      call toomanywarnings(405)
    endif
    call istop (ifail,1)
  endif
  if (le-leIn.le.0) then
    if (warnings(406).le.warning_limit) then
      warnings(406) = warnings(406) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 406: define_process_rcl called for a process ', &
                         'without outgoing particles'
      write(nx,*)
      call toomanywarnings(406)
    endif
    call istop (ifail,1)
  endif

   legsIn(      prTot) = leIn
  legsOut(      prTot) = le - leIn
      par(:le,  prTot) = pa
      hel(:le,  prTot) = he
   resMax(      prTot) = nd0
   binRes(1:nd0,prTot) = binnd(1:nd0)
   parRes(1:nd0,prTot) = pond(1:nd0)
   NoffPh(      prTot) = noff
   polproj(1:nd0,prTot) = pol(1:nd0)

  if     (order.eq.'LO' ) then; loop(prTot) = .false.
  elseif (order.eq.'NLO') then; loop(prTot) = .true.; loopMax = .true.
  else
    if (warnings(407).le.warning_limit) then
      warnings(407) = warnings(407) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 407: define_process_rcl called at the wrong '// &
                         'loop order '
      write(nx,*) "           (accepted values are order = 'LO','NLO')"
      write(nx,*)
      call toomanywarnings(407)
    endif
    call istop (ifail,1)
  endif

  legsMax = maxval(legsIn+legsOut)

  end subroutine define_process_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_offshell_photons_rcl (npr,n)

  ! This subroutine sets the number of off-shell photons for process
  ! "npr" to "n".
  ! The value for the number of off-shell photons extracted by
  ! define_process_rcl is overwritten.

  integer, intent(in)  :: npr
  integer, intent(out) :: n

  integer :: pr,i,nphotons

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(501).le.warning_limit) then
      warnings(501) = warnings(501) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 501: set_offshell_photons_rcl called '// &
                                     'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(501)
    endif
    call istop (ifail,1)
  endif

  nphotons = 0
  do i = 1, legsIn(pr)+legsOut(pr)
    if (par(i,pr).eq.16) nphotons = nphotons + 1
  enddo
  if ( n.lt.0 .or. n.gt.nphotons ) then
    if (warnings(502).le.warning_limit) then
      warnings(502) = warnings(502) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a)') ' ERROR 502: set_offshell_photons_rcl called '// &
                                  'with a wrong number of off-shell photons'
      write(nx,*)
      call toomanywarnings(502)
    endif
    call istop (ifail,1)
  endif

  NoffPh(pr) = n

  end subroutine set_offshell_photons_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The next subroutines sets the EW renormalization scheme as a mixed
! scheme for a given process "npr"
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_gfermi_mixed_scheme_rcl (npr,Nal0,NalZ,NalMS,NLOs)

  ! This subroutine sets the EW renormalization scheme as a mixed
  ! scheme for process "npr" and takes alpha_GF ("algf") as reference
  ! alpha.
  ! All other alpha's are computed from alpha_GF according to:
  ! al0  = algf * (1 - deltar)
  ! alZ  = algf * (1 - deltar) / (1 - dalZ)
  ! alMS = algf * (1 - deltar) / (1 - dalMS)
  ! where deltar, dalZ, dalMS relate the renormalizations of the
  ! electric charge according to:
  ! e_bare = eGF (1+dZeGF) = e0 (1+dZe0) = eZ (1+dZeZ) = eMS (1+dZeMS)
  ! dZe0 = dZeGF + deltar/2 = dZeZ + dalZ/2 = dZeMS + dalMS/2
  ! The subroutine determines that the process "npr" is renormalized
  ! in the mixed scheme according to:
  ! alpha^N = algf^(N-Nal0-NalZ-NalMS) *
  !           al0^Nal0 * alZ^NalZ * alMS^NalMS
  ! The argument "NLOs", if present, determines which alpha has to be
  ! used for NLO corrections ("alphaNLO") according to:
  ! NLOs = 1 -> alphaNLO = algf
  ! NLOs = 2 -> alphaNLO = al0
  ! NLOs = 3 -> alphaNLO = alZ
  ! NLOs = 4 -> alphaNLO = alMS
  ! If NLOs is not present, alphaNLO is set to algf.

  integer,           intent(in) :: npr,Nal0,NalZ,NalMS
  integer, optional, intent(in) :: NLOs

  integer :: pr,i

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(471).le.warning_limit) then
      warnings(471) = warnings(471) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 471: use_gfermi_mixed_scheme_rcl called '
      write(nx,'(a,i3)') '            with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(2)
    endif
    call istop (ifail,1)
  endif

  if (processes_generated) then
    if (warnings(472).le.warning_limit) then
      warnings(472) = warnings(472) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 472: Processes already generated.'
      write(nx,*) '             The call of use_gfermi_mixed_scheme_rcl ', &
                               'has no effects.'
       write(nx,*)
      call toomanywarnings(472)
    endif
  endif

  ew_reno_scheme = 1

  refscheme(pr) = 1

  if (Nal0.lt.0.or.NalZ.lt.0.or.NalMS.lt.0) then
    if (warnings(473).le.warning_limit) then
      warnings(473) = warnings(473) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 473: use_gfermi_mixed_scheme_rcl called with wrong '
      write(nx,*) '           arguments "Nal0", "NalZ" or "NalMS"'
      write(nx,*)
      call toomanywarnings(473)
    endif
    call istop (ifail,1)
  endif

  Nalpha(2,pr) = Nal0
  Nalpha(3,pr) = NalZ
  Nalpha(4,pr) = NalMS

  if (present(NLOs)) then
    if (NLOs.le.0.or.NLOs.gt.4) then
      if (warnings(474).le.warning_limit) then
        warnings(474) = warnings(474) + 1
        call openOutput
        write(nx,*)
        write(nx,*) &
          'ERROR 474: use_gfermi_mixed_scheme_rcl called with wrong argument "NLOs"'
        write(nx,*)
        call toomanywarnings(474)
      endif
      call istop (ifail,1)
    endif
    NLOscheme(pr) = NLOs
  else
    NLOscheme(pr) = refscheme(pr)
  endif

  end subroutine use_gfermi_mixed_scheme_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_alpha0_mixed_scheme_rcl (npr,Nalgf,NalZ,NalMS,NLOs)

  ! This subroutine sets the EW renormalization scheme as a mixed
  ! scheme for process "npr" and takes alpha(0) ("al0") as reference
  ! alpha.
  ! All other alpha's are computed from alpha(0) according to:
  ! algf = al0 / (1 - deltar)
  ! alZ  = al0 / (1 - dalZ)
  ! alMS = al0 / (1 - dalMS)
  ! where deltar, dalZ, dalMS relate the renormalizations of the
  ! electric charge according to:
  ! e_bare = eGF (1+dZeGF) = e0 (1+dZe0) = eZ (1+dZeZ) = eMS (1+dZeMS)
  ! dZe0 = dZeGF + deltar/2 = dZeZ + dalZ/2 = dZeMS + dalMS/2
  ! The subroutine determines that the process "npr" is renormalized
  ! in the mixed scheme according to:
  ! alpha^N = al0^(N-Nalgf-NalZ-NalMS) *
  !           algf^Nalgf * alZ^NalZ * alMS^NalMS
  ! The argument "NLOs", if present, determines which alpha has to be
  ! used for NLO corrections ("alphaNLO") according to:
  ! NLOs = 1 -> alphaNLO = algf
  ! NLOs = 2 -> alphaNLO = al0
  ! NLOs = 3 -> alphaNLO = alZ
  ! NLOs = 4 -> alphaNLO = alMS
  ! If NLOs is not present, alphaNLO is set to al0.

  integer,           intent(in) :: npr,Nalgf,NalZ,NalMS
  integer, optional, intent(in) :: NLOs

  integer :: pr,i

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(475).le.warning_limit) then
      warnings(475) = warnings(475) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 475: use_alpha0_mixed_scheme_rcl called '
      write(nx,'(a,i3)') '            with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(475)
    endif
    call istop (ifail,1)
  endif

  if (processes_generated) then
    if (warnings(476).le.warning_limit) then
      warnings(476) = warnings(476) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 476: Processes already generated.'
      write(nx,*) '             The call of use_alpha0_mixed_scheme_rcl ', &
                               'has no effects.'
      write(nx,*)
      call toomanywarnings(476)
    endif
  endif

  ew_reno_scheme = 2

  refscheme(pr) = 2

  if (Nalgf.lt.0.or.NalZ.lt.0.or.NalMS.lt.0) then
    if (warnings(477).le.warning_limit) then
      warnings(477) = warnings(477) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 477: use_alpha0_mixed_scheme_rcl called with wrong '
      write(nx,*) '           arguments "Nalgf", "NalZ" or "NalMS"'
      write(nx,*)
      call toomanywarnings(477)
    endif
    call istop (ifail,1)
  endif

  Nalpha(1,pr) = Nalgf
  Nalpha(3,pr) = NalZ
  Nalpha(4,pr) = NalMS

  if (present(NLOs)) then
    if (NLOs.le.0.or.NLOs.gt.4) then
      if (warnings(478).le.warning_limit) then
        warnings(478) = warnings(478) + 1
        call openOutput
        write(nx,*)
        write(nx,*) &
          'ERROR 478: use_alpha0_mixed_scheme_rcl called with wrong argument "NLOs"'
        write(nx,*)
        call toomanywarnings(478)
      endif
      call istop (ifail,1)
    endif
    NLOscheme(pr) = NLOs
  else
    NLOscheme(pr) = refscheme(pr)
  endif

  end subroutine use_alpha0_mixed_scheme_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_alphaz_mixed_scheme_rcl (npr,Nalgf,Nal0,NalMS,NLOs)

  ! This subroutine sets the EW renormalization scheme as a mixed
  ! scheme for process "npr" and takes alpha(M_Z) ("alZ") as reference
  ! alpha.
  ! All other alpha's are computed from alpha(M_Z) according to:
  ! algf = alZ * (1 - dalZ) / (1 - deltar)
  ! al0  = alZ * (1 - dalZ)
  ! alMS = alZ * (1 - dalZ) / (1 - dalMS)
  ! where deltar, dalZ, dalMS relate the renormalizations of the
  ! electric charge according to:
  ! e_bare = eGF (1+dZeGF) = e0 (1+dZe0) = eZ (1+dZeZ) = eMS (1+dZeMS)
  ! dZe0 = dZeGF + deltar/2 = dZeZ + dalZ/2 = dZeMS + dalMS/2
  ! The subroutine determines that the process "npr" is renormalized
  ! in the mixed scheme according to:
  ! alpha^N = alZ^(N-Nalgf-Nal0-NalMS) *
  !           algf^Nalgf * al0^Nal0 * alMS^NalMS
  ! The argument "NLOs", if present, determines which alpha has to be
  ! used for NLO corrections ("alphaNLO") according to:
  ! NLOs = 1 -> alphaNLO = algf
  ! NLOs = 2 -> alphaNLO = al0
  ! NLOs = 3 -> alphaNLO = alZ
  ! NLOs = 4 -> alphaNLO = alMS
  ! If NLOs is not present, alphaNLO is set to alZ.

  integer,           intent(in) :: npr,Nalgf,Nal0,NalMS
  integer, optional, intent(in) :: NLOs

  integer :: pr,i

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(479).le.warning_limit) then
      warnings(479) = warnings(479) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 479: use_alphaz_mixed_scheme_rcl called '
      write(nx,'(a,i3)') '            with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(479)
    endif
    call istop (ifail,1)
  endif

  if (processes_generated) then
    if (warnings(480).le.warning_limit) then
      warnings(480) = warnings(480) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 480: Processes already generated.'
      write(nx,*) '             The call of use_alphaz_mixed_scheme_rcl ', &
                               'has no effects.'
      write(nx,*)
      call toomanywarnings(480)
    endif
  endif

  ew_reno_scheme = 3

  refscheme(pr) = 3

  if (Nalgf.lt.0.or.Nal0.lt.0.or.NalMS.lt.0) then
    if (warnings(481).le.warning_limit) then
      warnings(481) = warnings(481) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 481: use_alphaz_mixed_scheme_rcl called with wrong '
      write(nx,*) '           arguments "Nalgf", "Nal0" or "NalMS"'
      write(nx,*)
      call toomanywarnings(481)
    endif
    call istop (ifail,1)
  endif

  Nalpha(1,pr) = Nalgf
  Nalpha(2,pr) = Nal0
  Nalpha(4,pr) = NalMS

  if (present(NLOs)) then
    if (NLOs.le.0.or.NLOs.gt.4) then
      if (warnings(482).le.warning_limit) then
        warnings(482) = warnings(482) + 1
        call openOutput
        write(nx,*)
        write(nx,*) &
          'ERROR 482: use_alphaz_mixed_scheme_rcl called with wrong argument "NLOs"'
        write(nx,*)
        call toomanywarnings(482)
      endif
      call istop (ifail,1)
    endif
    NLOscheme(pr) = NLOs
  else
    NLOscheme(pr) = refscheme(pr)
  endif

  end subroutine use_alphaz_mixed_scheme_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_alphaMSbar_mixed_scheme_rcl (npr,Nalgf,Nal0,NalZ,NLOs)

  ! This subroutine sets the EW renormalization scheme as a mixed
  ! scheme for process "npr" and takes alphaMSbar ("alMS") as
  ! reference alpha.
  ! All other alpha's are computed from alphaMSbar according to:
  ! algf   = alMS * (1 - dalMS) / (1 - deltar)
  ! al0    = alMS * (1 - dalMS)
  ! alZ    = alMS * (1 - dalMS) / (1 - dalZ)
  ! where deltar, dalZ, dalMS relate the renormalizations of the
  ! electric charge according to:
  ! e_bare = eGF (1+dZeGF) = e0 (1+dZe0) = eZ (1+dZeZ) = eMS (1+dZeMS)
  ! dZe0 = dZeGF + deltar/2 = dZeZ + dalZ/2 = dZeMS + dalMS/2
  ! The subroutine determines that the process "npr" is renormalized
  ! in the mixed scheme according to:
  ! alpha^N = alMS^(N-Nalgf-Nal0-NalZ) *
  !           algf^Nalgf * al0^Nal0 * alZ^NalZ
  ! The argument "NLOs", if present, determines which alpha has to be used
  ! for NLO corrections ("alphaNLO") according to:
  ! NLOs = 1 -> alphaNLO = algf
  ! NLOs = 2 -> alphaNLO = al0
  ! NLOs = 3 -> alphaNLO = alZ
  ! NLOs = 4 -> alphaNLO = alMS
  ! If NLOs is not present, alphaNLO is set to alMS.

  integer,           intent(in) :: npr,Nalgf,Nal0,NalZ
  integer, optional, intent(in) :: NLOs

  integer :: pr,i

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(483).le.warning_limit) then
      warnings(483) = warnings(483) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 483: use_alphaMSbar_mixed_scheme_rcl called '
      write(nx,'(a,i3)') '            with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(483)
    endif
    call istop (ifail,1)
  endif

  if (processes_generated) then
    if (warnings(484).le.warning_limit) then
      warnings(484) = warnings(484) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 484: Processes already generated.'
      write(nx,*) '             The call of use_alphaMSbar_mixed_scheme_rcl ', &
                               'has no effects.'
      write(nx,*)
      call toomanywarnings(484)
    endif
  endif

  ew_reno_scheme = 4

  refscheme(pr) = 4

  if (Nalgf.lt.0.or.Nal0.lt.0.or.NalZ.lt.0) then
    if (warnings(485).le.warning_limit) then
      warnings(485) = warnings(485) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 485: use_alphaMSbar_mixed_scheme_rcl called with wrong '
      write(nx,*) '           arguments "Nalgf", "Nal0" or "NalZ"'
      write(nx,*)
      call toomanywarnings(485)
    endif
    call istop (ifail,1)
  endif

  Nalpha(1,pr) = Nalgf
  Nalpha(2,pr) = Nal0
  Nalpha(3,pr) = NalZ

  if (present(NLOs)) then
    if (NLOs.le.0.or.NLOs.gt.4) then
      if (warnings(486).le.warning_limit) then
        warnings(486) = warnings(486) + 1
        call openOutput
        write(nx,*)
        write(nx,*) &
          'ERROR 486: use_alphaMSbar_mixed_scheme_rcl called with wrong argument "NLOs"'
        write(nx,*)
        call toomanywarnings(486)
      endif
      call istop (ifail,1)
    endif
    NLOscheme(pr) = NLOs
  else
    NLOscheme(pr) = refscheme(pr)
  endif

  end subroutine use_alphaMSbar_mixed_scheme_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_reference_alpha_rcl (npr,a)

  ! This subroutine extracts the actual value of the reference alpha,
  ! setting "a" to the value of the reference alpha for process "npr".

  integer,  intent(in)  :: npr
  real(dp), intent(out) :: a

  integer     :: pr,i
  real(dp)    :: rew2,rez2
  complex(dp) :: mw2,mz2

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(487).le.warning_limit) then
      warnings(487) = warnings(487) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 487: get_reference_alpha_rcl called '// &
                                     'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(487)
    endif
    call istop (ifail,1)
  endif

  if (refscheme(pr).eq.1) then
    if (algf.ge.0d0) then
      a = algf
    elseif (masses_in_gf_reno_scheme.eq.0) then
      rew2 = mass_w**2
      rez2 = mass_z**2
      a = gf*sq2/pi*rew2*( 1d0 - rew2/rez2 )
    elseif (masses_in_gf_reno_scheme.eq.1) then
      mw2 = mass_w**2 - cId0*width_w*mass_w
      mz2 = mass_z**2 - cId0*width_z*mass_z
      a = gf*sq2/pi*abs( mw2*( 1d0 - mw2/mz2 ) )
    endif
  elseif (refscheme(pr).eq.2) then
    a = al0
  elseif (refscheme(pr).eq.3) then
    a = alZ
  elseif (refscheme(pr).eq.4) then
    a = alMS
  endif

  end subroutine get_reference_alpha_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_NLOalpha_rcl (npr,a)

  ! This subroutine extracts the actual value of the NLO alpha,
  ! setting "a" to the value of the NLO alpha for process "npr".

  integer,  intent(in)  :: npr
  real(dp), intent(out) :: a

  integer     :: pr,i
  real(dp)    :: rew2,rez2
  complex(dp) :: mw2,mz2

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(499).le.warning_limit) then
      warnings(499) = warnings(499) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 499: get_NLOalpha_rcl called '// &
                                     'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(499)
    endif
    call istop (ifail,1)
  endif

  if (.not.processes_generated) then
    if (NLOscheme(pr).eq.1) then
      if (algf.ge.0d0) then
        a = algf
      elseif (masses_in_gf_reno_scheme.eq.0) then
        rew2 = mass_w**2
        rez2 = mass_z**2
        a = gf*sq2/pi*rew2*( 1d0 - rew2/rez2 )
      elseif (masses_in_gf_reno_scheme.eq.1) then
        mw2 = mass_w**2 - cId0*width_w*mass_w
        mz2 = mass_z**2 - cId0*width_z*mass_z
        a = gf*sq2/pi*abs( mw2*( 1d0 - mw2/mz2 ) )
      endif
    elseif (NLOscheme(pr).eq.2) then
      a = al0
    elseif (NLOscheme(pr).eq.3) then
      a = alZ
    elseif (NLOscheme(pr).eq.4) then
      a = alMS
    endif
  else
    a = alphai(NLOscheme(pr),pr)
  endif

  end subroutine get_NLOalpha_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_alphaGF_rcl (npr,a)

  ! This subroutine extracts the actual value of alpha_GF for process
  ! "npr", setting "a" equal to that value.

  integer,  intent(in)  :: npr
  real(dp), intent(out) :: a

  integer     :: pr,i
  real(dp)    :: rew2,rez2
  complex(dp) :: mw2,mz2

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(494).le.warning_limit) then
      warnings(494) = warnings(494) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 494: get_alphaGF_rcl called '// &
                                     'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(494)
    endif
    call istop (ifail,1)
  endif

  if (.not.processes_generated) then
    if (algf.ge.0d0) then
      a = algf
    elseif (masses_in_gf_reno_scheme.eq.0) then
      rew2 = mass_w**2
      rez2 = mass_z**2
      a = gf*sq2/pi*rew2*( 1d0 - rew2/rez2 )
    elseif (masses_in_gf_reno_scheme.eq.1) then
      mw2 = mass_w**2 - cId0*width_w*mass_w
      mz2 = mass_z**2 - cId0*width_z*mass_z
      a = gf*sq2/pi*abs( mw2*( 1d0 - mw2/mz2 ) )
    endif
  else
    a = alphai(1,pr)
  endif

  end subroutine get_alphaGF_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_alpha0_rcl (npr,a)

  ! This subroutine extracts the actual value of alpha(0) for process
  ! "npr", setting "a" equal to that value.

  integer,  intent(in)  :: npr
  real(dp), intent(out) :: a

  integer     :: pr,i

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(495).le.warning_limit) then
      warnings(495) = warnings(495) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 495: get_alpha0_rcl called '// &
                                     'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(495)
    endif
    call istop (ifail,1)
  endif

  if (.not.processes_generated) then
    a = al0
  else
    a = alphai(2,pr)
  endif

  end subroutine get_alpha0_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_alphaz_rcl (npr,a)

  ! This subroutine extracts the actual value of alpha(M_Z) for
  ! process "npr", setting "a" equal to that value.

  integer,  intent(in)  :: npr
  real(dp), intent(out) :: a

  integer     :: pr,i

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(496).le.warning_limit) then
      warnings(496) = warnings(496) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 496: get_alphaz_rcl called '// &
                                     'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(496)
    endif
    call istop (ifail,1)
  endif

  if (.not.processes_generated) then
    a = alZ
  else
    a = alphai(3,pr)
  endif

  end subroutine get_alphaz_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_alphaMSbar_rcl (npr,a)

  ! This subroutine extracts the actual value of alphaMSbar for
  ! process "npr", setting "a" equal to that value.

  integer,  intent(in)  :: npr
  real(dp), intent(out) :: a

  integer     :: pr,i

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(497).le.warning_limit) then
      warnings(497) = warnings(497) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 497: get_alphaMSbar_rcl called '// &
                                     'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(497)
    endif
    call istop (ifail,1)
  endif

  if (.not.processes_generated) then
    a = alMS
  else
    a = alphai(4,pr)
  endif

  end subroutine get_alphaMSbar_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_QrenMSbar_rcl (Q)

  ! This subroutine extracts the actual value of "Qren_alMS", setting
  ! "Q" equal to that value.

  real(dp), intent(out) :: Q

  Q = Qren_alMS

  end subroutine get_QrenMSbar_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The next subroutines select/unselect the contributions to the
! amplitude with specific powers of g_s (i.e. the strong coupling:
! g_s^2 = 4*pi*alpha_s).
! If none of them are called all contributions are selected
! (default).
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_gs_power_rcl (npr,gsarray)

  ! This subroutine selects specific powers of g_s, as
  ! specified by the argument "gsarray", for the process with
  ! identifier "npr".
  ! "gsarray(0:,0:1)" is an integer array whose values can be 0 or 1.
  ! The first entry is the power of g_s.
  ! The second entry is the loop order (0 for LO, 1 for NLO)
  ! Example:
  ! process 'd~ d -> u~ u' at NLO
  ! gsarray(0,0) = 1
  ! gsarray(1,0) = 0
  ! gsarray(2,0) = 1
  ! gsarray(0,1) = 0
  ! gsarray(1,1) = 0
  ! gsarray(2,1) = 0
  ! gsarray(3,1) = 0
  ! gsarray(4,1) = 1
  ! Here we select for computation the contributions with g_s^0 and
  ! g_s^2 for the tree amplitude and the g_s^4 contribution for the
  ! loop amplitude. All other contributions are not selected.

  integer, intent(in) :: npr,gsarray(0:,0:)

  integer :: pr,i,legs

  if (processes_generated) then
    if (warnings(408).le.warning_limit) then
      warnings(408) = warnings(408) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 408: Processes already generated.'
      write(nx,'(a,a,i3,a)') &
                 '              The call of set_gs_power_rcl ', &
                               'for process ',npr,' has no effects.'
      write(nx,*)
      call toomanywarnings(408)
    endif
  endif

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(409).le.warning_limit) then
      warnings(409) = warnings(409) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 409: set_gs_power_rcl called with '// &
                                 'undefined process index ',npr
      write(nx,*)
      call toomanywarnings(409)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  do i = 0,legs-2
    if (gsarray(i,0).ne.0.and.gsarray(i,0).ne.1) then
      if (warnings(410).le.warning_limit) then
        warnings(410) = warnings(410) + 1
        call openOutput
        write(nx,*)
        write(nx,*) 'ERROR 410: set_gs_power_rcl called with wrong arguments'
        write(nx,*)
        call toomanywarnings(410)
      endif
      call istop (ifail,1)
    endif
  enddo
  powgs(0:legs-2,0,pr) = gsarray(0:legs-2,0)

  if (loop(pr)) then
    do i = 0,legs
      if (gsarray(i,1).ne.0.and.gsarray(i,1).ne.1) then
        if (warnings(411).le.warning_limit) then
          warnings(411) = warnings(411) + 1
          call openOutput
          write(nx,*)
          write(nx,*) 'ERROR 411: set_gs_power_rcl called with wrong arguments'
          write(nx,*)
          call toomanywarnings(411)
        endif
        call istop (ifail,1)
      endif
    enddo
    powgs(0:legs,1,pr) = gsarray(0:legs,1)
  else
    powgs(0:legs,1,pr) = 0
  endif

  end subroutine set_gs_power_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine select_gs_power_BornAmpl_rcl (npr,gspower)

  ! This subroutine selects the contribution to the Born amplitude
  ! with g_s power "gspower", for the process with process number
  ! "npr" .
  ! The selection of the contributions to the loop amplitude are
  ! uneffected, as well as the selection of the other powers of g_s
  ! for the Born amplitude.

  integer, intent(in)  :: npr,gspower

  integer :: pr,i,legs

  if (processes_generated) then
    if (warnings(412).le.warning_limit) then
      warnings(412) = warnings(412) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 412: Processes already generated.'
      write(nx,'(a,a,a,i3,a)') &
                 '              The call of ', &
                               'select_gs_power_BornAmpl_rcl ', &
                               'for process ',npr,' has no effects.'
      write(nx,*)
      call toomanywarnings(412)
    endif
  endif

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(413).le.warning_limit) then
      warnings(413) = warnings(413) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 413: select_gs_power_BornAmpl_rcl called '
      write(nx,'(a,i3)') '            with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(413)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)
  if (gspower.lt.0.or.gspower.gt.legs-2) then
    if (warnings(414).le.warning_limit) then
      warnings(414) = warnings(414) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 414: select_gs_power_BornAmpl_rcl called with wrong gs power'
      write(nx,*)
      call toomanywarnings(414)
    endif
    call istop (ifail,1)
  endif

  powgs(gspower,0,pr) = 1

  end subroutine select_gs_power_BornAmpl_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine select_gs_power_LoopAmpl_rcl (npr,gspower)

  ! This subroutine selects the contribution to the loop amplitude
  ! with g_s power "gspower" for the process with process number
  ! "npr".
  ! The selection of the contributions to the Born amplitude remains
  ! uneffected, as well as the selection of the other powers of g_s
  ! for the loop amplitude.

  integer, intent(in)  :: npr,gspower

  integer :: pr,i,legs

  if (processes_generated) then
    if (warnings(415).le.warning_limit) then
      warnings(415) = warnings(415) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 415: Processes already generated.'
      write(nx,'(a,a,a,i3,a)') &
                 '              The call of ', &
                               'select_gs_power_LoopAmpl_rcl ', &
                               'for process ',npr,' has no effects.'
      write(nx,*)
      call toomanywarnings(415)
    endif
  endif

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(416).le.warning_limit) then
      warnings(416) = warnings(416) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 416: select_gs_power_LoopAmpl_rcl called '
      write(nx,'(a,i3)') '            with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(416)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)
  if (gspower.lt.0.or.gspower.gt.legs) then
    if (warnings(417).le.warning_limit) then
      warnings(417) = warnings(417) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 417: select_gs_power_LoopAmpl_rcl called with wrong gs power'
      write(nx,*)
      call toomanywarnings(417)
    endif
    call istop (ifail,1)
  endif

  powgs(gspower,1,pr) = 1

  end subroutine select_gs_power_LoopAmpl_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine unselect_gs_power_BornAmpl_rcl (npr,gspower)

  ! This subroutine unselects the contribution to the Born amplitude
  ! with g_s power "gspower" for the process with process number
  ! "npr".
  ! The selection of the contributions to the loop amplitude remains
  ! uneffected, as well as the selection of the other powers of g_s
  ! for the Born amplitude.

  integer, intent(in)  :: npr,gspower

  integer :: pr,i,legs

  if (processes_generated) then
    if (warnings(418).le.warning_limit) then
      warnings(418) = warnings(418) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 418: Processes already generated.'
      write(nx,'(a,a,a,i3,a)') &
                 '              The call of ', &
                               'unselect_gs_power_BornAmpl_rcl ', &
                               'for process ',npr,' has no effects.'
      write(nx,*)
      call toomanywarnings(418)
    endif
  endif

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(419).le.warning_limit) then
      warnings(419) = warnings(419) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 419: unselect_gs_power_BornAmpl_rcl called '
      write(nx,'(a,i3)') '            with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(419)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)
  if (gspower.lt.0.or.gspower.gt.legs-2) then
    if (warnings(420).le.warning_limit) then
      warnings(420) = warnings(420) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 420: unselect_gs_power_BornAmpl_rcl called with wrong gs power'
      write(nx,*)
      call toomanywarnings(420)
    endif
    call istop (ifail,1)
  endif

  powgs(gspower,0,pr) = 0

  end subroutine unselect_gs_power_BornAmpl_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine unselect_gs_power_LoopAmpl_rcl (npr,gspower)

  ! This subroutine unselects the contribution to the loop amplitude
  ! with g_s power "gspower" for the process with process number
  ! "npr".
  ! The selection of the contributions to the Born amplitude remains
  ! uneffected, as well as the selection of the other powers of g_s
  ! for the loop amplitude.

  integer, intent(in)  :: npr,gspower

  integer :: pr,i,legs

  if (processes_generated) then
    if (warnings(421).le.warning_limit) then
      warnings(421) = warnings(421) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 421: Processes already generated.'
      write(nx,'(a,a,a,i3,a)') &
                 '              The call of ', &
                               'unselect_gs_power_LoopAmpl_rcl ', &
                               'for process ',npr,' has no effects.'
      write(nx,*)
      call toomanywarnings(421)
    endif
  endif

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(422).le.warning_limit) then
      warnings(422) = warnings(422) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 422: unselect_gs_power_LoopAmpl_rcl called '
      write(nx,'(a,i3)') '            with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(422)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)
  if (gspower.lt.0.or.gspower.gt.legs) then
    if (warnings(423).le.warning_limit) then
      warnings(423) = warnings(423) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 423: unselect_gs_power_LoopAmpl_rcl called with wrong gs power'
      write(nx,*)
      call toomanywarnings(423)
    endif
    call istop (ifail,1)
  endif

  powgs(gspower,1,pr) = 0

  end subroutine unselect_gs_power_LoopAmpl_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine select_all_gs_powers_BornAmpl_rcl (npr)

  ! This subroutine selects all contribution to the Born amplitude
  ! (with any g_s power) for the process with process number "npr".
  ! The selection of the contributions to the loop amplitude remains
  ! uneffected.

  integer, intent(in)  :: npr

  integer :: pr,i

  if (processes_generated) then
    if (warnings(424).le.warning_limit) then
      warnings(424) = warnings(424) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 424: Processes already generated.'
      write(nx,'(a,a,a,i3,a)') &
                 '              The call of ', &
                               'select_all_gs_powers_BornAmpl_rcl ', &
                               'for process ',npr,' has no effects.'
      write(nx,*)
      call toomanywarnings(424)
    endif
  endif

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(425).le.warning_limit) then
      warnings(425) = warnings(425) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 425: select_all_gs_powers_BornAmpl_rcl called '
      write(nx,'(a,i3)') '            with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(425)
    endif
    call istop (ifail,1)
  endif

  powgs(0:,0,pr) = 1

  end subroutine select_all_gs_powers_BornAmpl_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine select_all_gs_powers_LoopAmpl_rcl (npr)

  ! This subroutine selects all contribution to the loop amplitude
  ! (with any g_s power) for the process with process number "npr".
  ! The selection of the contributions to the Born amplitude remains
  ! uneffected.

  integer, intent(in)  :: npr

  integer :: pr,i

  if (processes_generated) then
    if (warnings(426).le.warning_limit) then
      warnings(426) = warnings(426) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 426: Processes already generated.'
      write(nx,'(a,a,a,i3,a)') &
                 '              The call of ', &
                               'select_all_gs_powers_LoopAmpl_rcl ', &
                               'for process ',npr,' has no effects.'
      write(nx,*)
      call toomanywarnings(426)
    endif
  endif

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(427).le.warning_limit) then
      warnings(427) = warnings(427) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 427: select_all_gs_powers_LoopAmpl_rcl called '
      write(nx,'(a,i3)') '            with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(427)
    endif
    call istop (ifail,1)
  endif

  powgs(0:,1,pr) = 1

  end subroutine select_all_gs_powers_LoopAmpl_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine unselect_all_gs_powers_BornAmpl_rcl (npr)

  ! This subroutine unselects all contribution to the Born amplitude
  ! (with any g_s power) for the process with process number "npr".
  ! The selection of the contributions to the loop amplitude remains
  ! uneffected.

  integer, intent(in)  :: npr

  integer :: pr,i

  if (processes_generated) then
    if (warnings(428).le.warning_limit) then
      warnings(428) = warnings(428) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 428: Processes already generated.'
      write(nx,'(a,a,a,i3,a)') &
                 '              The call of ', &
                               'unselect_all_gs_powers_BornAmpl_rcl ', &
                               'for process ',npr,' has no effects.'
      write(nx,*)
      call toomanywarnings(428)
    endif
  endif

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(429).le.warning_limit) then
      warnings(429) = warnings(429) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 429: unselect_all_gs_powers_BornAmpl_rcl '
      write(nx,'(a,i3)') '            called with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(429)
    endif
    call istop (ifail,1)
  endif

  powgs(0:,0,pr) = 0

  end subroutine unselect_all_gs_powers_BornAmpl_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine unselect_all_gs_powers_LoopAmpl_rcl (npr)

  ! This subroutine unselects all contribution to the loop amplitude
  ! (with any g_s power) for the process with process number "npr".
  ! The selection of the contributions to the Born amplitude remains
  ! uneffected.

  integer, intent(in)  :: npr

  integer :: pr,i

  if (processes_generated) then
    if (warnings(430).le.warning_limit) then
      warnings(430) = warnings(430) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 430: Processes already generated.'
      write(nx,'(a,a,a,i3,a)') &
                 '              The call of ', &
                               'unselect_all_gs_powers_LoopAmpl_rcl ', &
                               'for process ',npr,' has no effects.'
      write(nx,*)
      call toomanywarnings(430)
    endif
  endif

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(431).le.warning_limit) then
      warnings(431) = warnings(431) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 431: unselect_all_gs_powers_LoopAmpl_rcl '
      write(nx,'(a,i3)') '            called with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(431)
    endif
    call istop (ifail,1)
  endif

  powgs(0:,1,pr) = 0

  end subroutine unselect_all_gs_powers_LoopAmpl_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine split_collier_cache_rcl (npr,n)

  ! This subroutine splits the cache of collier for process "npr"
  ! in "n" parts.

  integer, intent(in)  :: npr,n

  integer, allocatable :: nCacheTmp(:)
  logical, allocatable :: cacheOnTmp(:)
  integer              :: pr,i

  if (processes_generated) then
    if (warnings(432).le.warning_limit) then
      warnings(432) = warnings(432) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 432: Processes already generated.'
      write(nx,'(a,a,a,i3,a)') &
                 '              The call of ', &
                               'split_collier_cache_rcl ', &
                               'for process ',npr,' has no effects.'
      write(nx,*)
      call toomanywarnings(432)
    endif
  endif

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(433).le.warning_limit) then
      warnings(433) = warnings(433) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 433: split_collier_cache_rcl called '
      write(nx,'(a,i3)') '            with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(433)
    endif
    call istop (ifail,1)
  endif

  if (n.lt.0) then
    if (warnings(434).le.warning_limit) then
      warnings(434) = warnings(434) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 434: split_collier_cache_rcl called with wrong arguments'
      write(nx,*)
      call toomanywarnings(434)
    endif
    call istop (ifail,1)
  endif

  if (.not.allocated(nCache)) then
    allocate(nCache(prTot))
    nCache = 1
  else if (size(nCache) .ne. prTot) then
    allocate(nCacheTmp(size(nCache)))
    nCacheTmp(:) = nCache(:)
    deallocate(nCache)
    allocate(nCache(prTot))
    nCache = 1
    nCache(1:size(nCacheTmp)) = nCacheTmp(:)
  endif

  if (.not.allocated(cacheOn)) then
    allocate(cacheOn(1:prTot))
    cacheOn = .true.
  else if (size(cacheOn) .ne. prTot) then
    allocate(cacheOnTmp(size(cacheOn)))
    cacheOnTmp(:) = cacheOn(:)
    deallocate(cacheOn)
    allocate(cacheOn(prTot))
    cacheOn = .true.
    cacheOn(1:size(cacheOnTmp)) = cacheOnTmp(:)
  endif

  if (n.eq.0) then
    cacheOn(pr) = .false.
  else
    nCache(pr) = n
  endif

  end subroutine split_collier_cache_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module process_definition_rcl

!#####################################################################



