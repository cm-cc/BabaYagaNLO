!#####################################################################
!!
!!  File  process_generation_rcl.f90
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

  module process_generation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl
  use collier_interface_rcl
  use tables_rcl
  use currents_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine generate_processes_rcl

  ! This subroutine generates all defined processes.
  ! It must be called after all processes are defined.

  integer       :: i,pr,gs,ne
  real(sp)      :: timeGENin,timeGENout
  character(5)  :: regel,regmu,regta,regu,regd,regc,regs,regt,regb
  character(30) :: fmt1,fmt2,fmt3

  call cpu_time (timeGENin)

  if (pureQED .and. ew_reno_scheme.le.3) ew_reno_scheme = 2
  do pr = 1,prTot
    if (pureQED .and. refscheme(pr).le.3) refscheme(pr) = 2
  enddo

  !   NPh(pr) = number of external photons for process pr
  ! NonPh(pr) = number of external on-shell photons for process pr
  allocate (NPh(prTot),NonPh(prTot))
  NPh = 0
  do pr = 1,prTot
    do i = 1,legsIn(pr)+legsOut(pr)
      if (par(i,pr).eq.16) NPh(pr) = NPh(pr) + 1
    enddo
    NonPh(pr) = NPh(pr) - NoffPh(pr)
  enddo

  do pr = 1,prTot

    ! If the renormalization scheme of process "pr" is not a mixed
    ! scheme (i.e. is a fixed scheme, i.e. refscheme(pr)=0) and
    ! if the fixed scheme is not the alpha(0) scheme, Nalpha(2,pr)
    ! is set equal to the number of on-shell photons
    if ( refscheme(pr).eq.0 .and. ew_reno_scheme.ne.2 ) &
      Nalpha(2,pr) = NonPh(pr)

    ! If the renormalization scheme of process "pr" is a fixed scheme
    ! (i.e. refscheme(pr)=0), refscheme(pr) is set equal to the
    ! ew_reno_scheme
    if ( refscheme(pr).eq.0 ) refscheme(pr) = ew_reno_scheme

    ! If the renormalization scheme for NLO alpha of process "pr"
    ! has not been fixed (i.e. NLOscheme(pr)=0), NLOscheme(pr) is set
    ! equal to refscheme(pr)
    if ( NLOscheme(pr).eq.0 ) NLOscheme(pr) = refscheme(pr)

  enddo

  ! tables
  call binaries_tables
  if (loopMax) call tensors_tables
  call particles_tables
  call masses_tables
  call couplings_tables

  if (collier_ct) call initialize_collier_ct

  allocate (  eRefFactor(prTot)); eRefFactor = 1d0
  allocate ( eFactor(1:4,prTot)); eFactor    = 1d0
  if (loopMax .or. (sum(Nalpha).gt.0) ) then
    allocate (eNLOfactor(prTot)); eNLOfactor = 1d0
    call counterterms_tables
  endif

  if (regf(26).eq.2) then; regel = 'light'; else; regel = ''; endif
  if (regf(27).eq.2) then; regmu = 'light'; else; regmu = ''; endif
  if (regf(28).eq.2) then; regta = 'light'; else; regta = ''; endif
  if (regf(23).eq.2) then; regu  = 'light'; else; regu  = ''; endif
  if (regf(29).eq.2) then; regd  = 'light'; else; regd  = ''; endif
  if (regf(24).eq.2) then; regc  = 'light'; else; regc  = ''; endif
  if (regf(30).eq.2) then; regs  = 'light'; else; regs  = ''; endif
  if (regf(25).eq.2) then; regt  = 'light'; else; regt  = ''; endif
  if (regf(31).eq.2) then; regb  = 'light'; else; regb  = ''; endif

  call openOutput

  write(nx,*)

  write(nx,'(1x,75("-"))')

  fmt1 = '(2x,a,g21.14,12x,a,g21.14)'
  fmt2 = '(2x,a,g21.14,2x,a)'
  fmt3 = '(2x,a,g21.14,2x,a,5x,a,g21.14)'
  write(nx,'(2x,a)') 'Pole masses and widths [GeV]:'
  write(nx,fmt1) 'M_Z   =', mass_z,         'Width_Z   =', width_z
  write(nx,fmt1) 'M_W   =', mass_w,         'Width_W   =', width_w
  write(nx,fmt1) 'M_H   =', mass_h,         'Width_H   =', width_h
  write(nx,fmt2) 'm_e   =', mass_el, regel
  write(nx,fmt3) 'm_mu  =', mass_mu, regmu, 'Width_mu  =', width_mu
  write(nx,fmt3) 'm_tau =', mass_ta, regta, 'Width_tau =', width_ta
  write(nx,fmt2) 'm_u   =', mass_u,  regu
  write(nx,fmt2) 'm_d   =', mass_d,  regd
  write(nx,fmt3) 'm_c   =', mass_c,  regc,  'Width_c   =', width_c
  write(nx,fmt2) 'm_s   =', mass_s,  regs
  write(nx,fmt3) 'm_t   =', mass_t,  regt,  'Width_t   =', width_t
  write(nx,fmt3) 'm_b   =', mass_b,  regb,  'Width_b   =', width_b

  write(nx,'(1x,75("-"))')

  if     (complex_mass_scheme.eq.0) then
    write(nx,*) ' Renormalization done in the on-shell mass scheme'
  elseif (complex_mass_scheme.eq.1) then
    write(nx,*) ' Renormalization done in the complex-mass scheme'
  endif

  write(nx,'(1x,75("-"))')

  if     (ew_reno_scheme.eq.1) then
    write(nx,'(2x,a,1x,a)') 'fixed EW Renormalization Scheme:','gfermi'
    write(nx,'(2x,a,g21.14,3x,a,g21.14,a)') 'alpha_Gf    =',algf,'(Gf =',gf,' GeV^-2)'
    write(nx,'(2x,a,g21.14,3x,a)')          'alpha(0)    =',al0,'(for on-shell photons)'
  elseif (ew_reno_scheme.eq.2) then
    write(nx,'(2x,a,1x,a)') 'fixed EW Renormalization Scheme:','alpha0'
    write(nx,'(2x,a,g21.14)')      'alpha(0)    =',al0
    write(nx,'(2x,a,g21.14,3x,a)') 'alpha_Gf    =',algf,'(for on-shell photons)'
  elseif (ew_reno_scheme.eq.3) then
    write(nx,'(2x,a,1x,a)') 'fixed EW Renormalization Scheme:','alphaZ'
    write(nx,'(2x,a,g21.14)')      'alpha(M_Z)  =',alZ
    write(nx,'(2x,a,g21.14,3x,a)') 'alpha(0)    =',al0,'(for on-shell photons)'
  elseif (ew_reno_scheme.eq.4) then
    write(nx,'(2x,a,1x,a)') 'fixed EW Renormalization Scheme:','MSbar'
    write(nx,'(2x,a,g21.14)')      'alpha_MSbar =',alMS
    write(nx,'(2x,a,g21.14,3x,a)') 'alpha(0)    =',al0,'(for on-shell photons)'
  endif

  write(nx,'(1x,75("-"))')

  if (Nfren.eq.-1) then
    write(nx,'(2x,a)') &
    'alpha_s Renormalization Scheme: Variable flavours Scheme'
  else
    write(nx,'(2x,a,i1,a)') &
    'alpha_s Renormalization Scheme: ',Nfren,'-flavours Scheme'
  endif
  write(nx,'(2x,a,g21.14,7x,a,g21.14,a)') &
  'alpha_s(Q) =',als,'Q =',Qren,' GeV'
  if (use_active_qmasses) then
    write(nx,'(2x,a)') 'Quark masses in the running of alpha_s [GeV]:'
    write(nx,'(2x,a,g21.14)') 'm_c =',mq(4)
    write(nx,'(2x,a,g21.14)') 'm_b =',mq(5)
    write(nx,'(2x,a,g21.14)') 'm_t =',mq(6)
  end if
  write(nx,'(1x,75("-"))')

  write(nx,'(2x,a,g21.14,7x,a,g21.14,a)') &
  'Delta_UV   =',DeltaUV,'mu_UV =',muUV,' GeV'

  write(nx,'(2x,a,g21.14)') &
  'Delta_IR^2 =',DeltaIR2
  write(nx,'(2x,a,g21.14,7x,a,g21.14,a)') &
  'Delta_IR   =',DeltaIR,'mu_IR =',muIR,' GeV'

  write(nx,'(1x,75("-"))')

  if (reg_soft.eq.1) then
    write(nx,'(2x,a)') &
    'Dimensional regularization for soft singularities'
  else
    write(nx,'(2x,a)') &
    'Mass regularization for soft singularities'
    write(nx,'(2x,a,g21.14,a)') 'Mass regulator = ',lambda,' GeV'
  endif

  write(nx,'(1x,75("-"))')

  write(nx,*)

  call generate_currents

  allocate (eLOfactor(0:maxval(gsTot(:,:)),prTot)); eLOfactor  = 1d0

  do pr = 1,prTot

    do gs = 0,gsTotEff(0,pr)
      if ( sum(Nalpha(:,pr)) .gt. max(0,legsIn(pr)+legsOut(pr)-2-gs) ) then
        comp0gs(gs,pr) = .false.
        if (warnings(461).le.warning_limit) then
          warnings(461) = warnings(461) + 1
          call openOutput
          write(nx,*)
          write(nx,*) 'WARNING 461:'
          write(nx,'(a,i0,a,i0,a)') &
                      ' The Born amplitude for process ',inpr(pr), &
                      ' with ',gs,' powers of gs is set to zero,'
          write(nx,*) 'because the mixed EW renormalization ', &
                      'scheme for this process '
          write(nx,*) 'has been called with too many powers ', &
                      'for the EW coupling constants.'
          write(nx,*)
          call toomanywarnings(461)
        endif
      else
        ne = legsIn(pr) + legsOut(pr) - 2 - gs
        if (ne.gt.0) eLOfactor(gs,pr) = eRefFactor(pr)**ne
        do i = 1,4
          eLOfactor(gs,pr) = eLOfactor(gs,pr) * eFactor(i,pr)
        enddo
      endif
    enddo

    if (loop(pr)) then
      do gs = 0,gsTotEff(1,pr)
        if ( sum(Nalpha(:,pr)) .gt. max(0,legsIn(pr)+legsOut(pr)-gs) ) then
          comp1gs(gs,pr) = .false.
          if (warnings(462).le.warning_limit) then
            warnings(462) = warnings(462) + 1
            call openOutput
            write(nx,*)
            write(nx,*) 'WARNING 462:'
            write(nx,'(a,i0,a,i0,a)') &
                        ' The one-loop amplitude for process ',inpr(pr), &
                        ' with ',gs,' powers of gs is set to zero,'
            write(nx,*) 'because the mixed EW renormalization ', &
                        'scheme for this process '
            write(nx,*) 'has been called with too many powers ', &
                        'of the EW coupling constants.'
            write(nx,*)
            call toomanywarnings(462)
          endif
        endif
        if ( (NLOscheme(pr).ne.refscheme(pr)) .and. &
             (.not.comp0gs(gs,pr)) ) then
          comp1gs(gs,pr) = .false.
          if (warnings(463).le.warning_limit) then
            warnings(463) = warnings(463) + 1
            call openOutput
            write(nx,*)
            write(nx,*) 'WARNING 463:'
            write(nx,*) 'In the mixed EW renormalization scheme with ', &
                        'an NLO alpha different from '
            write(nx,*) 'the reference alpha, all NLO corrections ', &
                        'are intended to be EW corrections '
            write(nx,*) '(just NLO contributions with the same gs ', &
                        'powers as the LO ones are computed). '
            write(nx,'(a,i0,a,i0,a)') &
                        ' Thus the one-loop amplitude for process ',inpr(pr), &
                        ' with ',gs,' powers of gs is set to zero.'
            write(nx,*)
            call toomanywarnings(463)
          endif
        endif
      enddo
    endif

  enddo

  call cpu_time (timeGENout)

  timeGEN = timeGENout - timeGENin

  allocate (timeTI(prTot)); timeTI = 0e0
  allocate (timeTC(prTot)); timeTC = 0e0

  ! collier initialization
  if (loopMax) call initialize_collier

  processes_generated = .true.

  if (print_process_generation_summary) then
    call process_generation_summary
    write(nx,*)
  end if

!+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

  contains

!+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

  subroutine process_generation_summary

    integer np,pr

    call openOutput

    write(nx,'(1x,a)') 'Process generation summary:'
    np = prTot
    if (np .eq. 1) then
      write(nx,'(1x,a)') trim(adjustl(to_str(np))) // ' process defined'
    else
      write(nx,'(1x,a)') trim(adjustl(to_str(np))) // ' processes defined'
    end if

    do pr = 1, np
      write(nx,'(2x,a)') ''
      write(nx,'(2x,a)') "Process " // trim(adjustl(to_str(inpr(pr)))) // &
                         ": "// trim(adjustl(process(pr)))
      write(nx,'(2x,a)') 'Tree currents    = ' // trim(adjustl(to_str(w0Tot(pr))))
      write(nx,'(2x,a)') 'Tree branches    = ' // trim(adjustl(to_str(bm0prTot(pr))))

      if (loop(pr)) then
        write(nx,'(2x,a)') 'Loop currents    = ' // trim(adjustl(to_str(w1Tot(pr))))
        write(nx,'(2x,a)') 'Loop branches    = ' // trim(adjustl(to_str(bm1prTot(pr))))
      endif

      if (loop(pr)) then
        write(nx,'(2x,a)') 'Tensor integrals = ' // trim(adjustl(to_str(tiTot(pr))))
      end if
      write(nx,'(2x,a)') 'Helicities       = ' // trim(adjustl(to_str(cfTot(pr))))
      if (csTot(pr) .ne. pCsTot(pr)) then
        write(nx,'(2x,a)') 'Colourflows      = ' // trim(adjustl(to_str(csTot(pr)))) &
            // ' ('// trim(adjustl(to_str(pCsTot(pr)))) // ' computed)'
      else
        write(nx,'(2x,a)') 'Colourflows      = ' // trim(adjustl(to_str(csTot(pr))))
      end if

    end do
    write(nx,'(1x,75("-"))')

  end subroutine process_generation_summary

!+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

  function to_str(i)
    integer, intent(in) :: i
    character(10)       :: to_str
    write(to_str,'(i10)') i
  end function to_str

  end subroutine generate_processes_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_deltar_rcl (npr,d)

  ! This subroutine extracts the actual value of deltar for process
  ! "npr", setting the value of "d" to deltar

  integer,     intent(in)  :: npr
  complex(dp), intent(out) :: d

  integer :: pr,i

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(488).le.warning_limit) then
      warnings(488) = warnings(488) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i0)') ' ERROR 488: get_deltar_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(488)
    endif
    call istop (ifail,1)
  endif

  if (.not.allocated(deltarpr)) then
    if (warnings(489).le.warning_limit) then
      warnings(489) = warnings(489) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 489: deltar not computed yet.'
      write(nx,*)
      call toomanywarnings(489)
    endif
  else
    d = deltarpr(pr)
  endif

  end subroutine get_deltar_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_dalZ_rcl (npr,d)

  ! This subroutine extracts the actual value of dalZ for process
  ! "npr", setting the value of "d" to dalZ.

  integer,  intent(in)     :: npr
  complex(dp), intent(out) :: d

  integer :: pr,i

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(490).le.warning_limit) then
      warnings(490) = warnings(490) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i0)') ' ERROR 490: get_dalZ_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(490)
    endif
    call istop (ifail,1)
  endif

  if (.not.allocated(dalZpr)) then
    if (warnings(491).le.warning_limit) then
      warnings(491) = warnings(491) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 491: dalZ not computed yet.'
      write(nx,*)
      call toomanywarnings(491)
    endif
  else
    d = dalZpr(pr)
  endif

  end subroutine get_dalZ_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_dalMS_rcl (npr,d)

  ! This subroutine extracts the actual value of dalMS for process
  ! "npr", setting the value of "d" to dalMS.

  integer,     intent(in)  :: npr
  complex(dp), intent(out) :: d

  integer :: pr,i

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(492).le.warning_limit) then
      warnings(492) = warnings(492) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i0)') ' ERROR 492: get_dalMS_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(492)
    endif
    call istop (ifail,1)
  endif

  if (.not.allocated(dalMSpr)) then
    if (warnings(493).le.warning_limit) then
      warnings(493) = warnings(493) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'WARNING 493: dalMS not computed yet.'
      write(nx,*)
      call toomanywarnings(493)
    endif
  else
    d = dalMSpr(pr)
  endif

  end subroutine get_dalMS_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine process_exists_rcl (npr,exists)

  ! returns whether a process with id "npr" exists after generation.

  integer, intent(in)          :: npr
  logical, intent(out)         :: exists

  integer :: i,pr

  if (.not.processes_generated) then
    if (warnings(441).le.warning_limit) then
      warnings(441) = warnings(441) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 441: Call of get_momenta_rcl not allowed:'
      write(nx,*) '           Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(441)
    endif
    call istop (ifail,1)
  endif

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    if (warnings(442).le.warning_limit) then
      warnings(442) = warnings(442) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i0)') ' ERROR 442: get_momenta_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(442)
    endif
    call istop (ifail,1)
  endif

  exists = prexists(pr)

  end subroutine process_exists_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module process_generation_rcl

!#####################################################################



