!#####################################################################
!!
!!  File  globals_rcl.f90
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

  module globals_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!---------------------------------------------------------------------
! Parameters
!---------------------------------------------------------------------

  character(len=9) :: version_rcl = "1.5.0"

  integer,     parameter :: sp = kind (10e0) ! single precision
  integer,     parameter :: dp = kind (23d0) ! double    "

  complex(dp), parameter :: c0d0 = (0d0,0d0)
  complex(dp), parameter :: c1d0 = (1d0,0d0)
  complex(dp), parameter :: cId0 = (0d0,1d0)

  real(dp),    parameter :: sq2 = 1.414213562373095048801689d0
  real(dp),    parameter :: pi  = 3.141592653589793238462643d0

  complex(dp), parameter :: csq2 = (1.414213562373095048801689d0,0d0)

  real(dp),    parameter :: zerocut = 1d-15  ! zero for if-conditions
  real(dp),    parameter :: zerocheck = 1d-7 ! zero for check of phase-space point
  integer,     parameter :: infty = 16       ! terms in series one loop

  integer,     parameter :: nFs = 43        ! number of fields
  integer,     parameter :: nCs =  3        ! number of colours

  integer,     parameter :: nOpenedDef = 99 ! maximal number of
                                            ! output files

!---------------------------------------------------------------------
! Variables for the output file output.rcl
!---------------------------------------------------------------------

  integer       :: RecolaScreen = 1

  character(999) :: outputfile = 'output.rcl'
  integer        :: nx = 934758

  integer       :: nOpened = 0
  character(99) :: nameOpened(nOpenedDef) = ''

!---------------------------------------------------------------------
! Variable for the output directory of collier
! collier_output_dir = 'default': The output directory of collier
!                                 will be the default of collier
! collier_output_dir = 'dir': The output directory of collier will be
!                             'dir', where 'dir' can be a relative or
!                             absolute path.
!---------------------------------------------------------------------

  character(999) :: collier_output_dir = 'default'

!---------------------------------------------------------------------
! Variable to check for late set of input parameters
!---------------------------------------------------------------------

  logical :: processes_generated = .false., &
             changed_lambda      = .false., &
             changed_DeltaUV     = .false., &
             changed_muUV        = .false., &
             changed_DeltaIR     = .false., &
             changed_muIR        = .false.

!---------------------------------------------------------------------
! General variables
!---------------------------------------------------------------------

! compute IR poles
! 0 -> No IR poles computed
! 1 -> IR poles computed via rescaling (3-times)
! 2 -> IR poles computed via I-dip (not implemented yet)
  integer :: compute_ir_poles = 0

! Optimization level for colour algebra
! 0 -> no optimization
! 1 -> optimization of currents just differing by the colour structure
! 2 -> optimization of external U(1)-gluons and of currents just
!      differing by the colour structure
  integer :: colour_optimization = 2

! Optimization level for helicities
! 0 -> no optimization, all helicities configuration of the currents
!      are computed
! 1 -> use helicity conservation for massless fermions, together with
!      a fermion flow for all fermions
! 2 -> use helicity conservation for massless fermions, together with
!      a fermion flow for fermions belonging to a family with a massless
!      and a massive particle
  integer :: helicity_optimization = 2

! Rescaling factor of the couplings of the standard model
  complex(dp), save :: &
      coupling3(nFs,nFs,nFs) = c1d0, &
      coupling4(nFs,nFs,nFs,nFs) =  c1d0

! Resonant particles
  logical :: resPar(nFs) = .false.

! Set the polarization vector of vector bosons to be longitudinal,
! in order to check Ward identities
! longitudinal =   0 -> no vector bosons are longitudinal; this is
!                       the default
! longitudinal = 'i' -> if particle 'i' is a vector boson, its
!                       polarization vector is set equal to its
!                       momentum
! longitudinal = 111 -> For all vector bosons the polarization
!                       vectors are set equal to their momentum
  integer :: longitudinal = 0

! set whether the longitudinal polarization is used exclusively in the loop
! amplitude
  logical :: longitudinal_nlo = .true.

! Separation of QED vs WEAK contributions.
! REMARK:
! The separation is consistently implemented just for processes
! without phi+, phi- W+ W- (neither external nor virtual) at leading
! order). If these particles are present at leading order, both
! variable "qed" and "weak" must be set to ".true."
  logical :: loopQED  = .true.
  logical :: loopWEAK = .true.

! Computation of pure QED contributions (at LO and NLO)
  logical :: pureQED = .false.

! Choose regularization scheme:
! reguScheme = 1 -> four-dimensional helicity (FDH), not checked
! reguScheme = 2 -> 't Hooft-Veltman scheme (HV)
! ATTENTION: FDH scheme not checked
  integer :: reguScheme = 2

! Switch for IR-resonances. For processes with resonances:
! if resIR=F, only normal resonant contributions are selected
! if resIR=T, also IR-resonant contributions are selected.
  logical :: resIR = .false.

! Switch for listings the vertices appeared in the processes under
! consideration in the files vertices2.txt, vertices3.txt,
! vertices4.txt.
  logical :: list_vertices = .false.

! The checks with the code "Pole" need some special treatment.
  logical :: check_Pole = .false.

! Switch for printing the numerical values of the counterterms
  logical :: write_counterterms = .false.

! Init collier after generation process
! 1 -> COLI
! 2 -> DD
! 3 -> COLI+DD
  integer :: collier_mode = 1

! Compute the counterterms using the collier library
  logical :: collier_ct = .true.

!---------------------------------------------------------------------
! Error options
!---------------------------------------------------------------------

! The ifail flag regulates the behaviour of Recola in case of error.
! On input (set by the user calling set_ifail_rcl):
! ifail = -1 -> the code does not stop, when an error occurs
! ifail =  0 -> the code stops, when an error occurs
! On output (which can be got by the user calling get_ifail_rcl):
! ifail = 0 -> no errors occured
! ifail = 1 -> an error occured because of a wrong call of user's
!              subroutines
! ifail = 2 -> an error occured because of an internal bug of Recola
  integer :: ifail = 0

!---------------------------------------------------------------------
! Warnings
!---------------------------------------------------------------------

  integer, dimension(600) :: warnings = 0
  integer, parameter      :: warning_limit = 100

!---------------------------------------------------------------------
! input
!---------------------------------------------------------------------

  real(dp) :: mq(6)

!---------------------------------------------------------------------
! process_definition
!---------------------------------------------------------------------

! number of processes
  integer :: prTot = 0

  integer, allocatable   :: inpr(:),legsIn(:),legsOut(:),par(:,:),   &
                            hel(:,:),resMax(:),binRes(:,:),          &
                            parRes(:,:),powgs(:,:,:)
  logical, allocatable   :: loop(:)
  character, allocatable :: process(:)*99

! all-process variables
  integer :: legsMax
  logical :: loopMax = .false.

! mixed renormalization shemes
  integer, allocatable   :: refscheme(:),Nalpha(:,:),NLOscheme(:),   &
                            NoffPh(:)

! stop and print error if a process does not exist
  logical :: stop_on_nonexisting_process = .false.

  logical :: print_process_generation_summary = .true.

!---------------------------------------------------------------------
! tables
!---------------------------------------------------------------------

  ! binaries_tables
  integer, allocatable      :: levelLeg(:),vectorLeg(:,:),           &
                               firstNumber(:),firstGap(:),           &
                               firstNumbers(:,:),firstGaps(:,:)

  ! tensors_tables
  integer                   :: riTot,gg(0:3,0:3)
  integer, allocatable      :: RtoS(:),riMin(:),riMax(:),ri(:,:),    &
                               RItoR(:),RItoI(:),incRI(:,:)
  logical, allocatable      :: firstRI(:,:)

  ! particles_tables
  integer                   :: parKind(nFs),threeQ(nFs),Qghost(nFs),Nq
  logical                   :: charged(nFs)
  character(8)              :: cpar(nFs)
  character(2)              :: cftype(nFs),cftype2(nFs)
  real(dp)                  :: Nf

  ! masses_tables
  integer                   :: nmasses,nmf(nFs),regf(nFs)
  logical                   :: ffpar(nFs)
  complex(dp)               :: cm2f(nFs),cm2pf(nFs),                 &
                               mw1,mw2,mw3,mw4,mw6,                  &
                               mz1,mz2,mz4,mz6,mh2,mh4,              &
                               st,st2,st3,st4,st6,st8,               &
                               ct,ct2,ct3,ct4,ct6,stct,              &
                               mu1(3),mu2(3),mu4(3),                 &
                               ml1(3),ml2(3),ml4(3),                 &
                               md1(3),md2(3),md4(3)
  complex(dp), allocatable  :: cm2n(:)

  ! couplings_tables
  integer                   :: Nlq0_alMS,Nlq_alMS,                   &
                               Nfren0,Nlq,Nlq0,CMscheme
  integer, allocatable      :: Nlq0R(:,:)
  logical                   :: use_active_qmasses = .false.
  real(dp)                  :: pi2,lam,algf=-1d0,alpha=0d0,          &
                               Qren0_alMS,                           &
                               Ql,Ql2,Ql3,Ql4,                       &
                               Qu,Qu2,Qu3,Qu4,Qd,Qd2,Qd3,Qd4,        &
                               I3n,I3n2,I3n3,I3n4,                   &
                               I3l,I3l2,I3l3,I3l4,                   &
                               I3u,I3u2,I3u3,I3u4,                   &
                               I3d,I3d2,I3d3,I3d4,                   &
                               als0,Qren0,Nc,Nc2,Cf,Ca,              &
                               beta0(6),beta1(6),rb1(6),             &
                               mq2(6),mq2_alMS(6),Qq2_alMS(6)
  real(dp), allocatable     :: als0R(:,:),Qren0R(:,:)
  complex(dp)               :: gpn,gpn2,gpl,gpl2,                    &
                               gpu,gpu2,gpd,gpd2,                    &
                               gmn,gmn2,gml,gml2,                    &
                               gmu,gmu2,gmd,gmd2

  ! counterterms_tables
  complex(dp)               :: dt,dZh,dmh2,dZg,dZgs,dZgs0,dZaa,      &
                               dZaz,dZza,dZe,dZzz,dmz2,dZw,          &
                               dmw2,dmw,dst,dct,dZp,dZp0,            &
                               dgpn,dgpl,dgpu,dgpd,dgmn,dgml,        &
                               dgmu,dgmd
  complex(dp), dimension(3) :: dZnL,dZnR,dZnLc,dZnRc,                &
                               dZuL,dZuLqcd,dZuR,dZuRqcd,            &
                               dZuLc,dZuLcqcd,dZuRc,dZuRcqcd,        &
                               dmu,dmuqcd,                           &
                               dZlR,dZlL,dZlRc,dZlLc,dml,            &
                               dZdL,dZdLqcd,dZdR,dZdRqcd,            &
                               dZdLc,dZdLcqcd,dZdRc,dZdRcqcd,        &
                               dmd,dmdqcd
  complex(dp), allocatable  :: dZgs0R(:)
  real(dp),    allocatable  :: eRefFactor(:),eFactor(:,:),           &
                               alphai(:,:),                          &
                               eLOfactor(:,:),eNLOfactor(:)
  complex(dp), allocatable  :: deltarpr(:),dalZpr(:),dalMSpr(:),     &
                               dZaapr(:),dZepr(:),DdZe(:),           &
                               dgpnpr(:),dgplpr(:),                  &
                               dgpupr(:),dgpdpr(:),                  &
                               dgmnpr(:),dgmlpr(:),                  &
                               dgmupr(:),dgmdpr(:)

!---------------------------------------------------------------------
! process_generation
!---------------------------------------------------------------------

  ! number of external photons (total and on-shell)
  integer, allocatable      :: NPh(:),NonPh(:)

  ! vertices_tables
  logical, allocatable      :: ve2ct(:,:,:), &
                               ve2r2(:,:,:), &
                               ve3tr(:,:,:,:), &
                               ve3ct(:,:,:,:), &
                               ve3r2(:,:,:,:), &
                               ve4tr(:,:,:,:,:), &
                               ve4ct(:,:,:,:,:), &
                               ve4r2(:,:,:,:,:)
  integer                   :: cfMax,csMax,gsMax,gs2Max,             &
                               minlmu(-1:1,nFs),maxlmu(-1:1,nFs),    &
                               fhmin(nFs),fhmax(nFs),ram0,ram1,ram2
  integer, allocatable      :: newleg(:,:),oldleg(:,:),newbin(:,:),  &
                               oldbin(:,:),cfTot(:),                 &
                               csTot(:),pCsTot(:),csIa(:,:,:),       &
                               csIq(:,:,:),nIa(:,:),pIa(:,:,:),      &
                               w0eTot(:),heli(:,:,:),                &
                               lpmax(:),dualheli(:,:,:),cd0sMax(:),  &
                               loopCoef(:,:),w0Tot(:),w1TotMax(:,:), &
                               w1Tot(:),bm0prTot(:),bm1prTot(:),     &
                               w0last(:,:),parw0e(:,:),binw0e(:,:),  &
                               legw0e(:,:),helw0e(:,:),mosm0(:,:),   &
                               binsm0(:,:,:),parsm0(:,:,:),          &
                               gsIncsm0(:,:),gssm0(:,:),cssm0(:,:),  &
                               dasd0(:,:),gssd0(:,:),cssd0(:,:),     &
                               mosm1(:,:),binsm1(:,:,:),parsm1(:,:), &
                               gsIncsm1(:,:),rankInsm1(:,:),         &
                               rankOutsm1(:,:),gssm1(:,:),           &
                               cssm1(:,:),tism1(:,:),dasd1(:,:),     &
                               rankOutsd1(:,:),gssd1(:,:),           &
                               cssd1(:,:),tisd1(:,:),c0EffMax(:),    &
                               c0TOlp(:,:),bm0min(:,:,:,:),          &
                               bm0max(:,:,:,:),bd0min(:,:,:,:),      &
                               bd0max(:,:,:,:),cEffMax(:),cTOt(:,:), &
                               cTOfh(:,:),cTOih1(:,:),               &
                               sbm0(:,:),w0inbm0(:,:,:),             &
                               w0outbm0(:,:),typebm0(:,:),sbd0(:,:), &
                               w0outbd0(:,:),sbm1(:,:),w1inbm1(:,:), &
                               w0inbm1(:,:,:),w1outbm1(:,:),         &
                               typebm1(:,:),sbd1(:,:),w1outbd1(:,:), &
                               modaTot(:),riwMax(:,:),tiTot(:),      &
                               ritiMax(:),legsti(:,:),momsti(:,:,:), &
                               vmti(:,:,:),rankti(:,:),              &
                               colcoef(:,:,:),                       &
                               gsTot(:,:),gs2Tot(:,:),               &
                               gsTotEff(:,:),gs2TotEff(:,:),         &
                               nCache(:),nCacheTot(:),tiCache(:)
  logical, allocatable      :: defp2bin(:,:),defresbin(:,:),         &
                               xsm0(:,:),sesd0(:,:),comp0gs(:,:),    &
                               comp1gs(:,:),ferloopsm1(:,:),         &
                               ferloopsd1(:,:),winitbm0(:,:),        &
                               winitbd0(:,:),winitbm1(:,:),          &
                               winitbd1(:,:),zeroLO(:),cacheOn(:)

  complex (dp), allocatable :: cmONS2(:,:),cmREG2(:,:),cosm0(:,:,:), &
                               cosm1(:,:,:)

  real    (dp), allocatable :: p2bin(:,:),pspbin(:,:),mONS(:,:),     &
                               facsd0(:,:),facsd1(:,:),facIa(:,:,:), &
                               colcoefc(:,:,:,:,:),factor(:)

  logical,      allocatable :: prexists(:)


  type branch_bounds
    integer :: bmin
    integer :: bmax
  end type branch_bounds

  type bbranch
    type(branch_bounds), allocatable :: conf(:,:,:)
  end type bbranch

  type(bbranch), allocatable :: bm1_b(:)
  type(bbranch), allocatable :: bd1_b(:)

!---------------------------------------------------------------------
! Fermion-flow connection
!---------------------------------------------------------------------

  integer, allocatable      :: qflow(:,:)

!---------------------------------------------------------------------
! Vector polarisation projection
!---------------------------------------------------------------------

  logical :: drop_goldstone_transverse = .true.
  logical :: drop_goldstone_longitudinal = .true.
  integer, allocatable      :: polproj(:,:),polprojin(:,:)

!---------------------------------------------------------------------
! Variables for time statistics
!---------------------------------------------------------------------

  ! time statistics for tensor integrals and tensor coefficients
  real(sp)              :: timeGEN
  real(sp), allocatable :: timeTI(:), timeTC(:)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_recola_version_rcl(ret_version)

  character(len=10), intent(out) :: ret_version

  ret_version = version_rcl

  end subroutine get_recola_version_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_EWloop_QED_rcl

  ! This subroutine selects the pure QED corrections in the process
  ! generation, counterterms and rational terms.

  loopQED = .true.
  loopWEAK = .false.

  end subroutine set_EWloop_QED_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_EWloop_WEAK_rcl

  ! This subroutine selects the pure weak corrections in the process
  ! generation, counterterms and rational terms.

  loopQED = .false.
  loopWEAK = .true.

  end subroutine set_EWloop_WEAK_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_EWloop_EW_rcl

  ! This subroutine selects all electroweak corrections, which is the
  ! default, in the process generation, counterterms and rational
  ! terms.

  loopQED = .true.
  loopWEAK = .true.

  end subroutine set_EWloop_EW_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_pure_QED_rcl

  ! This subroutine selects the pure QED contributions to all
  ! processes at LO and NLO.
  ! This automatically implies that the renormalization scheme is
  ! moved to the alpha0 scheme

  pureQED = .true.

  end subroutine set_pure_QED_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine unset_pure_QED_rcl

  ! This subroutine restore the computation of all electroweak (not
  ! just QED) contributions to all processes at LO and NLO, which is
  ! the default

  pureQED = .false.

  end subroutine unset_pure_QED_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function npar (c)

  character(*), intent(in) :: c
  integer                  :: npar

  if     (c.eq.'H')       then; npar = 11
  elseif (c.eq.'p0')      then; npar = 12
  elseif (c.eq.'p+')      then; npar = 13
  elseif (c.eq.'p-')      then; npar = 14
  elseif (c.eq.'g')       then; npar = 15
  elseif (c.eq.'A')       then; npar = 16
  elseif (c.eq.'Z')       then; npar = 17
  elseif (c.eq.'W+')      then; npar = 18
  elseif (c.eq.'W-')      then; npar = 19
  elseif (c.eq.'nu_e')    then; npar = 20
  elseif (c.eq.'nu_mu')   then; npar = 21
  elseif (c.eq.'nu_tau')  then; npar = 22
  elseif (c.eq.'u')       then; npar = 23
  elseif (c.eq.'c')       then; npar = 24
  elseif (c.eq.'t')       then; npar = 25
  elseif (c.eq.'e-')      then; npar = 26
  elseif (c.eq.'mu-')     then; npar = 27
  elseif (c.eq.'tau-')    then; npar = 28
  elseif (c.eq.'d')       then; npar = 29
  elseif (c.eq.'s')       then; npar = 30
  elseif (c.eq.'b')       then; npar = 31
  elseif (c.eq.'nu_e~')   then; npar = 32
  elseif (c.eq.'nu_mu~')  then; npar = 33
  elseif (c.eq.'nu_tau~') then; npar = 34
  elseif (c.eq.'u~')      then; npar = 35
  elseif (c.eq.'c~')      then; npar = 36
  elseif (c.eq.'t~')      then; npar = 37
  elseif (c.eq.'e+')      then; npar = 38
  elseif (c.eq.'mu+')     then; npar = 39
  elseif (c.eq.'tau+')    then; npar = 40
  elseif (c.eq.'d~')      then; npar = 41
  elseif (c.eq.'s~')      then; npar = 42
  elseif (c.eq.'b~')      then; npar = 43
  else
    if (warnings(199).le.warning_limit) then
      warnings(199) = warnings(199) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 199: Particle "',trim(c),'" not defined'
      write(nx,*)
      call toomanywarnings(199)
    endif
    call istop (ifail,1)
  endif

  end function npar

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function anti (p) result (a)

  integer, intent (in) :: p ! read   particle      tag
  integer              :: a ! return anti-particle tag

  select case (p)
  case ( 1: 5);  a = p +  5 ! FP ghost fields
  case ( 6:10);  a = p -  5 !  "   "     "
  case (13,18);  a = p +  1 ! phi+ and W+
  case (14,19);  a = p -  1 ! phi- and W-
  case (20:31);  a = p + 12 ! fermions
  case (32:43);  a = p - 12 ! anti-fermions
  case default;  a = p      ! other particles
  end select

  end function anti

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine istop (i,n)

  integer, intent(inout) :: i
  integer, intent(in)    :: n

  select case (i)
    case (0);     i = n; stop
    case (-1);    i = n
    case default; i = n
  end select

  end subroutine istop

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine openOutput

  integer      :: n,i
  logical      :: nxopen
  character(7) :: status

  if (RecolaScreen.eq.0) then
    write(*,*)
    write(*,'(1x,75("x"))')
    write(*,'(26x,a)') ' _    _   _   _        _ '
    write(*,'(26x,a)') '| )  |_  |   | |  |   |_|'
    write(*,'(26x,a)') '| \  |_  |_  |_|  |_  | |'
    write(*,*)
    write(*,'(16x,a)') 'REcursive Computation of One-Loop Amplitudes'
    write(*,*)
    write(*,'(33x,a)') 'Version ' // trim(version_rcl)
    write(*,*)
    write(*,'(7x,a)') 'by S.Actis, A.Denner, L.Hofer, J.-N.Lang, A.Scharf, S.Uccirati'
    write(*,*)
    write(*,'(1x,75("x"))')
    write(*,*)
    if (trim(outputfile).ne.'*') then
      write(*,*)
      write(*,'(1x,75("x"))')
      write(*,*)
      write(*,*) 'From now on the output of RECOLA is written ', &
                 'to the file ',trim(outputfile)
      write(*,*)
      write(*,'(1x,75("x"))')
      write(*,*)
    endif
    RecolaScreen = 1
  endif

  if (trim(outputfile).eq.'*') then

    nx = 6

  else

    inquire(file=trim(outputfile),number=n)

    if (n.eq.-1) then ! outputfile is not open

      status = 'replace'
      do i = 1,nOpened
        if (trim(nameOpened(i)).eq.trim(outputfile)) then
          status='unknown'; exit
        endif
      enddo

      if (status.eq.'replace') then
        ! "outputfile" is open for the first time
        if (nOpened.ge.nOpenedDef) then
          write(*,*)
          write(*,*) &
            'CODE ERROR: Too many output-files. ', &
                        'Increse nOpenedDef in file globals_rcl.f90'
          write(*,*)
          stop
        else
          nOpened = nOpened + 1
          nameOpened(nOpened) = trim(outputfile)
        endif
      endif

      if (nx.eq.6) nx = 934758

      inquire(nx,opened=nxopen)
      do while (nxopen)
        nx = nx + 1
        inquire(nx,opened=nxopen)
      enddo

      open ( unit=nx,file=trim(outputfile),      &
             status=trim(status),position='append' )

    else              ! outputfile is open

      nx = n

    endif

  endif

  end subroutine openOutput

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine toomanywarnings(iwarn)

  integer :: iwarn

  if (warnings(iwarn).eq.warning_limit+1) then
    call openOutput
    write(nx,*)
    write(nx,'(a,i0,a)') ' Too many WARNINGs of type = ',iwarn,' for this run:'
    write(nx,'(a,i0,a)') ' Further WARNINGs of type = ',iwarn,' will not be printed'
!    write(nx,*) ' Too many WARNINGs for this run:'
!    write(nx,*) ' Further WARNINGs will not be printed'
    write(nx,*)
  endif

  end subroutine toomanywarnings

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine print_warning_summary_rcl(file_handle)

  integer, intent(in), optional :: file_handle
  integer :: i_warn,unit_number

  if(present(file_handle)) then
    unit_number = file_handle
  else
    call openOutput
    unit_number = nx
  end if

  write(unit_number,'(/1x,75("x"))')
  write (unit_number, *) "Summary of Recola errors and warnings"
  write (unit_number, *) "Error No.          #"
  do i_warn = 1, size(warnings)
    if (warnings(i_warn) > 0) &
       write (unit_number, '(i4,7x,i10)') i_warn, warnings(i_warn)
  end do
  write(unit_number,'(1x,75("x")/)')

  end subroutine print_warning_summary_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module globals_rcl

!#####################################################################




