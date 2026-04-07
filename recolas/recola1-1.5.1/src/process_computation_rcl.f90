!#####################################################################
!!
!!  File  process_computation_rcl.f90
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

  module process_computation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl
  use amplitude_rcl
  use wave_functions_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_resonant_squared_momentum_rcl (npr,res,ps)

  ! This subroutine acts on the resonance number "res" for the process
  ! with process number "npr". It sets the squared momentum of the
  ! denominator of the resonant propagator to "ps".
  ! The resonance-number "res" is defined from the process definition
  ! in the call of define_process_rcl. The first resonant particle
  ! defined there has "res"=1, the second resonant particle has
  ! "res=2", and so on. "ps" is a real number.

  integer,  intent(in)  :: npr,res
  real(dp), intent(in)  :: ps

  integer               :: pr,i,l,e0,e1,e2

  if (.not.processes_generated) then
    if (warnings(1).le.warning_limit) then
      warnings(1) = warnings(1) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 1: Call of set_resonant_squared_momentum_rcl not allowed:'
      write(nx,*) '         Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(1)
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
    if (warnings(2).le.warning_limit) then
      warnings(2) = warnings(2) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 2: set_resonant_squared_momentum_rcl '
      write(nx,'(a,i3)') '          called with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(2)
    endif
    call istop (ifail,1)
  endif

  l  = legsIn(pr) + legsOut(pr)

  if (res.lt.0.or.res.gt.l) then
    if (warnings(3).le.warning_limit) then
      warnings(3) = warnings(3) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 3: set_resonant_squared_momentum_rcl called ', &
                         'with wrong resonance number'
      write(nx,*)
      call toomanywarnings(3)
    endif
    call istop (ifail,1)
  endif

  if (ps.lt.0d0) then
    if (warnings(4).le.warning_limit) then
      warnings(4) = warnings(4) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 4: set_resonant_squared_momentum_rcl called ', &
                         'with wrong squared momentum'
      write(nx,*)
      write(nx,*)
      write(nx,*)
      call toomanywarnings(4)
    endif
    call istop (ifail,1)
  endif

  if (.not.resPar(parRes(res,pr))) then
    if (warnings(5).le.warning_limit) then
      warnings(5) = warnings(5) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 5: set_resonant_squared_momentum_rcl '
      write(nx,*) '         called for a particle not set as resonant'
      write(nx,*)
      call toomanywarnings(5)
    endif
    call istop (ifail,1)
  endif

  e0 = binRes(res,pr)
  e1 = newbin(e0,pr)
  e2 = 2**l - 1 - e1

  defp2bin(e1,pr) = .true.
     p2bin(e1,pr) = ps

  defp2bin(e2,pr) = .true.
     p2bin(e2,pr) = ps

  end subroutine set_resonant_squared_momentum_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_running_alphas_rcl (Q,Nf,lp)

  ! This subroutine can be called to compute the value for alpha_s at
  ! the scale "Q" employing the renormalization-group evolution at
  ! "lp" loops ("lp" can take the value "lp=1" or "lp=2").
  ! The integer argument "Nf" selects the flavour scheme according to
  ! the following rules:
  ! - "Nf = -1":
  !   The variable flavour scheme is selected, where the number of
  !   active flavours contributing to the running of alpha_s is set
  !   to the number of quarks lighter than the scale "Q".
  ! - "Nf = 3,4,5,6":
  !   The fixed Nf-flavour scheme is selected, where the number of
  !   active flavours contributing to the running of alpha_s is set
  !   to "Nf". In this case "Nf" cannot be smaller than the number of
  !   massless quarks (otherwise the code stops).


  real(dp), intent(in) :: Q
  integer,  intent(in) :: Nf,lp

  integer               :: Nl,Nl0,DNl,aDNl,i,n,k
  real(dp)              :: Q0,ia
  real(dp), allocatable :: bL(:)

  if (.not.processes_generated) then
    if (warnings(6).le.warning_limit) then
      warnings(6) = warnings(6) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 6: Call of compute_running_alphas_rcl not allowed:'
      write(nx,*) '         Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(6)
    endif
    call istop (ifail,1)
  endif

  if (Q.le.0d0) then
    if (warnings(7).le.warning_limit) then
      warnings(7) = warnings(7) + 1
      call openOutput
      write(nx,*)
      write(nx,*) &
        'ERROR 7: compute_running_alphas_rcl called with wrong scale'
      write(nx,*)
      call toomanywarnings(7)
    endif
    call istop (ifail,1)
  endif

  Nl0 = Nlq0
  select case (Nf)
  case (-1)
    Nl = 0
    do i = 1,6
      if (mq2(i).lt.Q**2) Nl = Nl + 1
    enddo
  case (3,4,5)
    if (mq2(Nf+1).eq.0d0) then
      if (warnings(8).le.warning_limit) then
        warnings(8) = warnings(8) + 1
        call openOutput
        write(nx,*)
        write(nx,*) 'ERROR 8: compute_running_alphas_rcl called with a wrong number '
        write(nx,*) '         of light flavours Nf (Nf can not be smaller than '
        write(nx,*) '         the number of massless quarks)'
        write(nx,*)
        call toomanywarnings(8)
      endif
      call istop (ifail,1)
    endif
    Nl = Nf
  case (6)
    Nl = Nf
  case default
    if (warnings(9).le.warning_limit) then
      warnings(9) = warnings(9) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 9: compute_running_alphas_rcl called with a wrong number of '
      write(nx,*) '         light flavours Nf (accepted values are Nf = -1,3,4,5,6)'
      write(nx,*)
      call toomanywarnings(9)
    endif
    call istop (ifail,1)
  end select
  DNl = Nl - Nl0
  aDNl = abs(DNl)

  Q0 = Qren0

  allocate ( bL(0:aDNl) )
  select case (DNl)
  case (0)
    bL(0) = beta0(Nl) * log(Q**2/Q0**2)
  case (1:)
    bL(0) = beta0(Nl0) * log(mq2(Nl0+1)/Q0**2)
    do k = 1,aDNl-1; n = Nl0 + k
      bL(k) = beta0(n) * log(mq2(n+1)/mq2(n))
    enddo
    bL(aDNl) = beta0(Nl) * log(Q**2/mq2(Nl))
  case (:-1)
    bL(0) = beta0(Nl0) * log(mq2(Nl0)/Q0**2)
    do k = 1,aDNl-1; n = Nl0 - k
      bL(k) = beta0(n) * log(mq2(n)/mq2(n+1))
    enddo
    bL(aDNl) = beta0(Nl) * log(Q**2/mq2(Nl+1))
  end select

  ia = pi/als0

  select case (lp)
  case (1)
    ia = ia + sum(bL(0:aDNl))
  case (2)
    do k = 0,aDNl; n = Nl0 + k*sign(1,DNl)
      ia = ia + bL(k) + rb1(n) * log( 1d0 + bL(k)/( ia + rb1(n) ) )
    enddo
  case default
    if (warnings(10).le.warning_limit) then
      warnings(10) = warnings(10) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 10: compute_running_alphas_rcl called with a wrong number of '
      write(nx,*) '          loops lp for the running (accepted values are lp = 1,2)'
      write(nx,*)
      call toomanywarnings(10)
    endif
    call istop (ifail,1)
  end select

  als   = pi/ia
  Qren  = Q
  Nfren = Nf
  Nlq   = Nl

  end subroutine compute_running_alphas_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_process_rcl (npr,p,order,A2,momenta_check)

  ! This subroutine computes the structure-dressed helicity
  ! amplitudes and the summed squared amplitude for the process with
  ! process number "npr".
  ! "p" are the external momenta of the process.
  ! "order" is the loop order. Accepted values are 'LO' and 'NLO'.
  ! "A2(0)" is the summed squared  LO amplitude.
  ! "A2(1)" is the summed squared NLO amplitude.
  ! Summed squared amplitude means:
  ! - summed over all computed powers of g_s
  ! - summed over colour and spins for outgoing particles and averaged
  !   for the incoming particles.
  ! If some particles were defined as polarized, no sum/average is
  ! performed on these particles.
  ! "momenta_check" tells whether the phase-space point is good
  ! (momenta_check=.true.) or bad (momenta_check=.false).

  integer,          intent(in)            :: npr
  real(dp),         intent(in)            :: p(0:,:)
  character(len=*), intent(in)            :: order
  real(dp),         intent(out), optional :: A2(0:1)
  logical,          intent(out), optional :: momenta_check

  integer  :: pr,i,legs,long_tmp
  real(dp) :: d1,d2


  if (.not.processes_generated) then
    if (warnings(11).le.warning_limit) then
      warnings(11) = warnings(11) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 11: Call of compute_process_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(11)
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
    if (warnings(12).le.warning_limit) then
      warnings(12) = warnings(12) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 12: compute_process_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(12)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if (size(p,2).ne.legs.or.size(p,1).ne.4) then
    if (warnings(13).le.warning_limit) then
      warnings(13) = warnings(13) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 13: compute_process_rcl called with wrong momenta'
      write(nx,*)
      call toomanywarnings(13)
    endif
    call istop (ifail,1)
  endif

  if (order.ne.'LO'.and.order.ne.'NLO') then
    if (warnings(14).le.warning_limit) then
      warnings(14) = warnings(14) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 14: compute_process_rcl called at the wrong '// &
                         'loop order '
      write(nx,*) "          (accepted values are order = 'LO','NLO')"
      write(nx,*)
      call toomanywarnings(14)
    endif
    call istop (ifail,1)
  endif

  call set_momenta (pr,p)
  if (present(momenta_check)) momenta_check = momcheck

  if (writeMat+writeMat2.ge.1) then
    call print_process_and_momenta (pr)
    call print_rescaling
    call print_parameters_change
  endif

  if (order .eq. 'NLO' .and. dynamic_settings .eq. 1 .and. &
      compute_ir_poles .eq. 1) then

    d1 = DeltaIR
    d2 = DeltaIR2
    call set_delta_ir_rcl(d1 + 100d0, d2)
    call compute_amplitude (pr,order)
    call rescale_amplitude (pr,order)
    call compute_squared_amplitude(pr,order)
    matrix2(0:gs2Tot(1,pr),5,pr) = matrix2(0:gs2Tot(1,pr),4,pr)

    call set_delta_ir_rcl(d1, d2 + 100d0)
    call compute_amplitude (pr,order)
    call rescale_amplitude (pr,order)
    call compute_squared_amplitude(pr,order)
    matrix2(0:gs2Tot(1,pr),6,pr) = matrix2(0:gs2Tot(1,pr),4,pr)

    call set_delta_ir_rcl(d1, d2)
    call compute_amplitude (pr,order)
    call rescale_amplitude (pr,order)

    if (writeMat.ge.1) call print_amplitude (pr,order)
    call compute_squared_amplitude(pr,order)

    ! IR1 pole
    matrix2(0:gs2Tot(1,pr),5,pr) =  &
      (matrix2(0:gs2Tot(1,pr),5,pr) - &
       matrix2(0:gs2Tot(1,pr),4,pr))/100d0
    ! IR2 pole
    matrix2(0:gs2Tot(1,pr),6,pr) =  &
      (matrix2(0:gs2Tot(1,pr),6,pr) - &
       matrix2(0:gs2Tot(1,pr),4,pr))/100d0
  else

    if (order .eq. 'NLO' .and. longitudinal_nlo .and. longitudinal .ne. 0) then
    ! If longitudinal and longitudinal_nlo is selected, the QED Ward identity
    ! is only applied to the one-loop matrix element and the LO matrix
    ! element is recomputed with the full polarization dependence.
      call compute_amplitude (pr,'NLO')
      long_tmp = longitudinal
      longitudinal = 0
      call compute_amplitude (pr,'LO')
      longitudinal = long_tmp
    else
      call compute_amplitude (pr,order)
    end if

    call rescale_amplitude (pr,order)

    if (writeMat.ge.1) call print_amplitude (pr,order)

    call compute_squared_amplitude (pr,order)

  end if

  if (writeMat2.ge.1) call print_squared_amplitude (pr,order)

  if (present(A2)) then
    A2(0) = sum(matrix2(0:gs2Tot(0,pr),0,pr))
    A2(1) = sum(matrix2(0:gs2Tot(1,pr),4,pr))
  endif

  end subroutine compute_process_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine rescale_process_rcl (npr,order,A2)

  ! This subroutine adjusts the results calculated by
  ! compute_process_rcl for a new value of alpha_s, rescaling the
  ! structure-dressed helicity amplitudes and recomputing the summed
  ! squared amplitude for the process with process number "npr".
  ! "order" and "A2" are the same has for compute_process_rcl.

  integer,          intent(in)            :: npr
  character(len=*), intent(in)            :: order
  real(dp),         intent(out), optional :: A2(0:1)

  integer     :: pr,i

  if (.not.processes_generated) then
    if (warnings(15).le.warning_limit) then
      warnings(15) = warnings(15) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 15: Call of rescale_process_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(15)
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
    if (warnings(16).le.warning_limit) then
      warnings(16) = warnings(16) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 16: rescale_process_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(16)
    endif
    call istop (ifail,1)
  endif

  if (order.ne.'LO'.and.order.ne.'NLO') then
    if (warnings(17).le.warning_limit) then
      warnings(17) = warnings(17) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 17: rescale_process_rcl called at the wrong '// &
                         'loop order '
      write(nx,*) "          (accepted values are order = 'LO','NLO')"
      write(nx,*)
      call toomanywarnings(17)
    endif
    call istop (ifail,1)
  endif

  if (writeMat+writeMat2.ge.1) call print_rescaling
  if (writeMat+writeMat2.ge.1) call print_parameters_change

  call rescale_amplitude (pr,order)

  if (writeMat.ge.1) call print_amplitude (pr,order)

  call compute_squared_amplitude (pr,order)

  if (writeMat2.ge.1) call print_squared_amplitude (pr,order)

  if (present(A2)) then
    A2(0) = sum(matrix2(0:gs2Tot(0,pr),0,pr))
    A2(1) = sum(matrix2(0:gs2Tot(1,pr),4,pr))
  endif

  end subroutine rescale_process_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_amplitude_rcl (npr,pow,order,colour,hel,A)

  ! This subroutine extracts a specific contribution to the amplitude
  ! of the process with process number "npr", according to the values
  ! of "pow", "order", "colour" and "hel".
  ! - "pow" is the power of gs of the contribution.
  ! - "order" is the loop-order of the contribution.
  !   It is a character variable with accepted values:
  !   'LO'     -> tree squared amplitude
  !   'NLO'    -> loop squared amplitude
  !   'NLO-D4' -> loop squared amplitude, bare loop contribution
  !   'NLO-CT' -> loop squared amplitude, counterterms contribution
  !   'NLO-R2' -> loop squared amplitude, rational terms
  !               contribution
  ! - "colour(1:l)" describes the colour structure and is a vector of
  !   type integer and length l, where each position in the vector
  !   corresponds to one of the l external particles of process "npr"
  !   (ordered as in the process definition).
  !   For colourless particles, incoming quarks and outgoing
  !   anti-quarks, the corresponding entry in "colour" must be 0.
  !   For all other particles (gluons, outgoing quarks and incoming
  !   anti-quarks), the entries must be, without repetition, the
  !   positions of the gluons, incoming quarks or outgoing
  !   anti-quarks in the process definition.
  ! - "hel" is the helicity of the contribution.
  !   It is a vector whose entries are the helicities of the
  !   particles:
  !   left-handed  fermions/antifermions -> -1
  !   right-handed fermions/antifermions -> +1
  !   logitudinal vector bosons          ->  0
  !   transverse vector bosons           -> -1, +1
  !   scalar particles                   ->  0
  !   Example:
  !   Process 'u g -> W+ d'
  !   "hel" = (/-1,+1,-1,-1/) means
  !   left-handed up-quark
  !   transverse gluon (helicity = +1)
  !   transverse W+ (helicity = -1)
  !   left-handed down-quark

  ! This routine must be called after compute_process_rcl or
  ! rescale_process_rcl for the process with process number "npr".

  integer,              intent(in)  :: npr,pow,colour(1:),hel(1:)
  character(len=*),     intent(in)  :: order
  complex(dp),          intent(out) :: A

  integer               :: pr,i,j,h,o,legs,cs,kf
  integer, allocatable  :: helk(:),cs0(:)
  character(2)          :: t,t2,t2c
  character(8)          :: c

  if (.not.processes_generated) then
    if (warnings(22).le.warning_limit) then
      warnings(22) = warnings(22) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 22: Call of get_amplitude_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(22)
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
    if (warnings(23).le.warning_limit) then
      warnings(23) = warnings(23) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 23: get_amplitude_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(23)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ( size(colour,1).ne.legs) then
    if (warnings(24).le.warning_limit) then
      warnings(24) = warnings(24) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 24: get_amplitude_rcl called with a wrong '// &
                         'colour vector (the number of '
      write(nx,*) '          entries must coincide with the number of '// &
                         'legs of the process)'
      write(nx,*)
      call toomanywarnings(24)
    endif
    call istop (ifail,1)
  endif
  do i = 1,legs
    if ( colour(i).lt.0.or.colour(i).gt.legs ) then
      if (warnings(25).le.warning_limit) then
        warnings(25) = warnings(25) + 1
        call openOutput
        write(nx,*)
        write(nx,*) 'ERROR 25: get_amplitude_rcl called with a wrong '// &
                           'colour vector (the entries can not '
        write(nx,*) '          be negative or bigger than the number '// &
                           'of legs of the process)'
        write(nx,*)
        call toomanywarnings(25)
      endif
      call istop (ifail,1)
    endif
    do j = 1,i-1
      if ( colour(i).eq.colour(j).and.colour(i).ne.0) then
        if (warnings(26).le.warning_limit) then
          warnings(26) = warnings(26) + 1
          call openOutput
          write(nx,*)
          write(nx,*) 'ERROR 26: get_amplitude_rcl called with a wrong '// &
                             'colour vector (the non-zero entries '// &
                             'must differ)'
          write(nx,*)
          call toomanywarnings(26)
        endif
        call istop (ifail,1)
      endif
    enddo
  enddo
  do i = 1,legs
    t2  = cftype2(par(i,pr))
    if (colour(i).ne.0) then; t2c = cftype2(par(colour(i),pr))
    else;                     t2c = '0'
    endif
    if ( ((t2c.ne.'g'.and.t2c.ne.'q' ).and.colour(i).ne.0) .or. &
         ((t2 .eq.'g'.or. t2 .eq.'q~').and.colour(i).eq.0) .or. &
         ((t2 .ne.'g'.and.t2 .ne.'q~').and.colour(i).ne.0)      &
       ) then
      if (warnings(27).le.warning_limit) then
        warnings(27) = warnings(27) + 1
        call openOutput
        write(nx,*)
        write(nx,*) 'ERROR 27: get_amplitude_rcl called with a wrong '// &
                           'colour vector'
        write(nx,*)
        call toomanywarnings(27)
      endif
      call istop (ifail,1)
    endif
  enddo

  allocate (cs0(legs)); cs0 = 0
  do i = 1,legs
    if (colour(i).ne.0) cs0(newleg(i,pr)) = newleg(colour(i),pr)
  enddo
  cs = 0
  do i = 1,csTot(pr)
    if ( sum(abs(csIq(1:legs,i,pr)-cs0(1:legs))).eq.0 ) then
      cs = i
      exit
    endif
  enddo
  deallocate (cs0)

  if (size(hel).ne.legs) then
    if (warnings(28).le.warning_limit) then
      warnings(28) = warnings(28) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 28: get_amplitude_rcl called with a wrong '// &
                         'helicity vector (the number of '
      write(nx,*) '          entries must coincide with the number of '// &
                         'legs of the process)'
      write(nx,*)
      call toomanywarnings(28)
    endif
    call istop (ifail,1)
  endif
  do i = 1,legs
    h = hel(i)
    c = cpar(par(i,pr))
    t = cftype(par(i,pr))
    if ( (h.lt.-1)              .or. &
         (h.gt.+1)              .or. &
         (c.eq.'g' .and.h.eq.0) .or. &
         (c.eq.'A' .and.h.eq.0) .or. &
         (t.eq.'f' .and.h.eq.0) .or. &
         (t.eq.'f~'.and.h.eq.0) .or. &
         (t.eq.'s' .and.h.ne.0)      &
       ) then
      if (warnings(29).le.warning_limit) then
        warnings(29) = warnings(29) + 1
        call openOutput
        write(nx,*)
        write(nx,*) 'ERROR 29: get_amplitude_rcl called with a wrong '// &
                           'helicity vector'
        write(nx,*)
        call toomanywarnings(29)
      endif
      call istop (ifail,1)
    endif
  enddo

  allocate (helk(legs))
  do i = 1,legsIn(pr)
    helk(newleg(i,pr)) = hel(i)
  enddo
  do i = legsIn(pr)+1,legs
    helk(newleg(i,pr)) = - hel(i)
  enddo
  kf = 0
  do i = 1,cfTot(pr)
    if ( sum(abs( heli(1:legs,i,pr) - helk(1:legs) )).eq.0 ) then
      kf = i
      exit
    endif
  enddo

  if     (order.eq.'LO'    ) then; o = 0
  elseif (order.eq.'NLO-D4') then; o = 1
  elseif (order.eq.'NLO-CT') then; o = 2
  elseif (order.eq.'NLO-R2') then; o = 3
  elseif (order.eq.'NLO'   ) then; o = 4
  else
    if (warnings(30).le.warning_limit) then
      warnings(30) = warnings(30) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 30: get_amplitude_rcl called at the '// &
                         'wrong loop order (accepted '
      write(nx,*) "          values are order = 'LO','NLO','NLO-D4',"// &
                         "'NLO-CT','NLO-R2')"
      write(nx,*)
      call toomanywarnings(30)
    endif
    call istop (ifail,1)
  endif

  if ((o.eq.0.and.(pow.lt.0.or.pow.gt.legs-2)) .or. &
      (o.gt.0.and.(pow.lt.0.or.pow.gt.legs  ))      &
     ) then
    if (warnings(31).le.warning_limit) then
      warnings(31) = warnings(31) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 31: get_amplitude_rcl called with wrong gs power'
      write(nx,*)
      call toomanywarnings(31)
    endif
    call istop (ifail,1)
  endif

  A = c0d0
  if (cs.ne.0.and.kf.ne.0) A = matrix(cs,pow,kf,o,pr)

  end subroutine get_amplitude_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_squared_amplitude_rcl (npr,pow,order,A2)

  ! This subroutine extracts the computed value of the summed squared
  ! amplitude with "pow" power of alpha_s at loop-order "order" for
  ! the process with process number "npr".
  ! "order: is a character variable with accepted values:
  ! 'LO'     -> tree squared amplitude
  ! 'NLO'    -> loop squared amplitude
  ! 'NLO-D4' -> loop squared amplitude, bare loop contribution
  ! 'NLO-CT' -> loop squared amplitude, counterterms contribution
  ! 'NLO-R2' -> loop squared amplitude, rational terms contribution
  ! If compute_ir_poles is activated, the IR poles can be obtained with
  ! 'NLO-IR1' -> loop squared amplitude, \epsilon_IR^-1 contribution
  ! 'NLO-IR2' -> loop squared amplitude, \epsilon_IR^-2 contribution

  ! This routine must be called after compute_process_rcl or
  ! rescale_process_rcl for the process with process number "npr".

  integer,          intent(in)  :: npr,pow
  character(len=*), intent(in)  :: order
  real(dp),         intent(out) :: A2

  integer :: pr,i,legs,gspower,o

  if (.not.processes_generated) then
    if (warnings(32).le.warning_limit) then
      warnings(32) = warnings(32) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 32: Call of get_squared_amplitude_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(32)
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
    if (warnings(33).le.warning_limit) then
      warnings(33) = warnings(33) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 33: get_squared_amplitude_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(33)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if     (order.eq.'LO'    )  then; o = 0
  elseif (order.eq.'NLO-D4')  then; o = 1
  elseif (order.eq.'NLO-CT')  then; o = 2
  elseif (order.eq.'NLO-R2')  then; o = 3
  elseif (order.eq.'NLO'   )  then; o = 4
  elseif (order.eq.'NLO-IR1') then; o = 5
  elseif (order.eq.'NLO-IR2') then; o = 6
  else
    if (warnings(34).le.warning_limit) then
      warnings(34) = warnings(34) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 34: get_squared_amplitude_rcl called at '// &
                         'the wrong loop order (accepted '
      write(nx,*) "          values are order = 'LO','NLO','NLO-D4',"// &
                         "'NLO-CT','NLO-R2','NLO-IR1','NLO-IR2')"
      write(nx,*)
      call toomanywarnings(34)
    endif
    call istop (ifail,1)
  endif

  if ((pow.lt.0)                                       .or. &
      (o.eq.0.and.                      pow.gt.legs-2) .or. &
      (o.gt.0.and.(.not.zeroLO(pr)).and.pow.gt.legs-1) .or. &
      (o.gt.0.and.      zeroLO(pr) .and.pow.gt.legs  )      &
     ) then
    if (warnings(35).le.warning_limit) then
      warnings(35) = warnings(35) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 35: get_squared_amplitude_rcl called', &
                        ' with wrong alphas power'
      write(nx,*)
      call toomanywarnings(35)
    endif
    call istop (ifail,1)
  endif

  gspower = 2*pow

  A2 = matrix2(gspower,o,pr)

  end subroutine get_squared_amplitude_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_polarized_squared_amplitude_rcl (npr,pow,order,hel,A2h)

  ! This subroutine extracts the computed value of the squared
  ! amplitude summed over colour with polarization "hel", with "pow"
  ! power of alpha_s at loop-order "order" for the process with
  ! process number "npr".
  ! - "pow" is the power of alpha_s of the contribution.
  ! - "order" is the loop-order of the contribution.
  !   It is a character variable with accepted values:
  !   'LO'     -> tree squared amplitude
  !   'NLO'    -> loop squared amplitude
  !   'NLO-D4' -> loop squared amplitude, bare loop contribution
  !   'NLO-CT' -> loop squared amplitude, counterterms contribution
  !   'NLO-R2' -> loop squared amplitude, rational terms
  !                 contribution
  ! - "hel" is the helicity of the contribution.
  !   It is a vector whose entries are the helicities of the
  !   particles:
  !   left-handed  fermions/antifermions -> -1
  !   right-handed fermions/antifermions -> +1
  !   logitudinal vector bosons          ->  0
  !   transverse vector bosons           -> -1, +1
  !   scalar particles                   ->  0
  !   Example:
  !   Process 'u g -> W+ d'
  !   "hel" = (/-1,+1,-1,-1/) means
  !   left-handed up-quark
  !   transverse gluon (helicity = +1)
  !   transverse W+ (helicity = -1)
  !   left-handed down-quark

  ! This routine must be called after compute_process_rcl or
  ! rescale_process_rcl for the process with process number "npr".

  integer,          intent(in)  :: npr,pow,hel(1:)
  character(len=*), intent(in)  :: order
  real(dp),         intent(out) :: A2h

  integer               :: pr,i,legs,gspower,kf,h,o
  character(2)          :: t
  character(8)          :: c
  integer, allocatable  :: helk(:)

  if (.not.processes_generated) then
    if (warnings(36).le.warning_limit) then
      warnings(36) = warnings(36) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 36: Call of get_polarized_squared_amplitude_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(36)
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
    if (warnings(37).le.warning_limit) then
      warnings(37) = warnings(37) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 37: get_polarized_squared_amplitude_rcl '
      write(nx,'(a,i3)') '           called with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(37)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if (size(hel).ne.legs) then
    if (warnings(38).le.warning_limit) then
      warnings(38) = warnings(38) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 38: get_polarized_squared_amplitude_rcl called '// &
                         'with a wrong '
      write(nx,*) '          helicity vector (the number of entries must '// &
                         'coincide '
      write(nx,*) '          with the number of legs of the process)'
      write(nx,*)
      call toomanywarnings(38)
    endif
    call istop (ifail,1)
  endif
  do i = 1,legs
    h = hel(i)
    c = cpar(par(i,pr))
    t = cftype(par(i,pr))
    if ( (h.lt.-1)              .or. &
         (h.gt.+1)              .or. &
         (c.eq.'g' .and.h.eq.0) .or. &
         (c.eq.'A' .and.h.eq.0) .or. &
         (t.eq.'f' .and.h.eq.0) .or. &
         (t.eq.'f~'.and.h.eq.0) .or. &
         (t.eq.'s' .and.h.ne.0)      &
       ) then
      if (warnings(39).le.warning_limit) then
        warnings(39) = warnings(39) + 1
        call openOutput
        write(nx,*)
        write(nx,*) 'ERROR 39: get_polarized_squared_amplitude_rcl called '// &
                           'with a wrong helicity vector'
        write(nx,*)
        call toomanywarnings(39)
      endif
      call istop (ifail,1)
    endif
  enddo

  allocate (helk(legs))
  do i = 1,legsIn(pr)
    helk(newleg(i,pr)) = hel(i)
  enddo
  do i = legsIn(pr)+1,legs
    helk(newleg(i,pr)) = - hel(i)
  enddo
  kf = 0
  do i = 1,cfTot(pr)
    if ( sum(abs( heli(:,i,pr) - helk(:) )).eq.0 ) then
      kf = i
      exit
    endif
  enddo

  if     (order.eq.'LO'    ) then; o = 0
  elseif (order.eq.'NLO-D4') then; o = 1
  elseif (order.eq.'NLO-CT') then; o = 2
  elseif (order.eq.'NLO-R2') then; o = 3
  elseif (order.eq.'NLO'   ) then; o = 4
  else
    if (warnings(40).le.warning_limit) then
      warnings(40) = warnings(40) + 1
      call openOutput
      write(nx,*)
      write(nx,*) "ERROR 40: get_polarized_squared_amplitude_rcl called "
      write(nx,*) "          at the wrong loop order (accepted values are "
      write(nx,*) "          order = 'LO','NLO','NLO-D4','NLO-CT','NLO-R2')"
      write(nx,*)
      call toomanywarnings(40)
    endif
    call istop (ifail,1)
  endif

  if ((pow.lt.0)                                       .or. &
      (o.eq.0.and.                      pow.gt.legs-2) .or. &
      (o.gt.0.and.(.not.zeroLO(pr)).and.pow.gt.legs-1) .or. &
      (o.gt.0.and.      zeroLO(pr) .and.pow.gt.legs  )      &
     ) then
    if (warnings(41).le.warning_limit) then
      warnings(41) = warnings(41) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR41: get_polarized_squared_amplitude_rcl', &
                        ' called with wrong alphas power'
      write(nx,*)
      call toomanywarnings(41)
    endif
    call istop (ifail,1)
  endif

  gspower = 2*pow

  A2h = 0d0
  if (kf.ne.0) A2h = matrix2h(gspower,kf,o,pr)

  end subroutine get_polarized_squared_amplitude_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_colour_correlation_rcl (npr,p,i1,i2,A2cc, &
                                             momenta_check)

  ! This subroutine computes the colour-correlated summed squared
  ! amplitude, between particle with leg number "i1" and particle
  ! with leg number "i2", for the process with process number "npr".
  ! "p" are the external momenta of the process.
  ! "A2cc" is the colour-correlated summed squared  LO
  ! amplitude.
  ! Summed squared amplitude means:
  ! - summed over all computed powers of g_s
  ! - summed over colour and spins for outgoing particles and averaged
  !   for the incoming particles.
  ! If some particles were defined as polarized, no sum/average is
  ! performed on these particles.
  ! "momenta_check" tells whether the phase-space point is good
  ! (momenta_check=.true.) or bad (momenta_check=.false).

  integer,  intent(in)            :: npr,i1,i2
  real(dp), intent(in)            :: p(0:,:)
  real(dp), intent(out), optional :: A2cc
  logical,  intent(out), optional :: momenta_check

  integer :: pr,i,legs,j1,j2

  if (.not.processes_generated) then
    if (warnings(42).le.warning_limit) then
      warnings(42) = warnings(42) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 42: Call of compute_colour_correlation_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(42)
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
    if (warnings(43).le.warning_limit) then
      warnings(43) = warnings(43) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 43: compute_colour_correlation_rcl called '
      write(nx,'(a,i3)') '           with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(43)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ((i1.le.0.or.i1.gt.legs).or.(i2.le.0.or.i2.gt.legs)) then
    if (warnings(44).le.warning_limit) then
      warnings(44) = warnings(44) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 44: compute_colour_correlation_rcl called ', &
                         'with wrong indices'
      write(nx,*)
      call toomanywarnings(44)
    endif
    call istop (ifail,1)
  endif

  if (size(p,2).ne.legs.or.size(p,1).ne.4) then
    if (warnings(45).le.warning_limit) then
      warnings(45) = warnings(45) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 45: compute_colour_correlation_rcl called ', &
                         'with wrong momenta'
      write(nx,*)
      call toomanywarnings(45)
    endif
    call istop (ifail,1)
  endif

  call set_momenta (pr,p)
  if (present(momenta_check)) momenta_check = momcheck

  if (writeCor.ge.1) then
    call print_process_and_momenta (pr)
    call print_rescaling
    call print_parameters_change
  endif

  call compute_amplitude (pr,'LO')

  call rescale_amplitude (pr,'LO')

  if (writeMat.ge.1) call print_amplitude (pr,'LO')

  call compute_squared_amplitude_cc (pr,i1,i2)

  if (writeCor.ge.1) then
    call print_squared_amplitude_cc (pr,i1,i2)
    write(nx,'(1x,75("x"))')
    write(nx,*)
    write(nx,*)
  endif

  if (present(A2cc)) then
    j1 = newleg(i1,pr)
    j2 = newleg(i2,pr)
    A2cc = sum(matrix2cc(0:gs2Tot(0,pr),j1,j2,pr))
  endif

  end subroutine compute_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine rescale_colour_correlation_rcl (npr,i1,i2,A2cc)

  ! This subroutine adjusts the results calculated by
  ! compute_colour_correlation_rcl for a new value of alpha_s,
  ! rescaling the structure-dressed helicity amplitudes and
  ! recomputing the colour-correlated summed squared amplitude for
  ! the process with process number "npr".
  ! "i1", "i2" and "A2cc" are the same has for
  ! compute_colour_correlation_rcl.

  integer,  intent(in)            :: npr,i1,i2
  real(dp), intent(out), optional :: A2cc

  integer :: pr,i,legs,j1,j2

  if (.not.processes_generated) then
    if (warnings(67).le.warning_limit) then
      warnings(67) = warnings(67) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 67: Call of rescale_colour_correlation_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(67)
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
    if (warnings(68).le.warning_limit) then
      warnings(68) = warnings(68) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 68: rescale_colour_correlation_rcl called '
      write(nx,'(a,i3)') '           with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(68)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ((i1.le.0.or.i1.gt.legs).or.(i2.le.0.or.i2.gt.legs)) then
    if (warnings(69).le.warning_limit) then
      warnings(69) = warnings(69) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 69: rescale_colour_correlation_rcl called ', &
                  'with wrong indices'
      write(nx,*)
      call toomanywarnings(69)
    endif
    call istop (ifail,1)
  endif

  if (writeCor.ge.1) call print_rescaling
  if (writeCor.ge.1) call print_parameters_change

  call rescale_amplitude (pr,'LO')

  if (writeMat.ge.1) call print_amplitude (pr,'LO')

  if (i2.eq.i1) then
    if (.not.allocated(matrix2cc)) then
       allocate (matrix2cc(0:gs2Max,1:legsMax,1:legsMax,prTot))
       matrix2cc = 0d0
     else
      matrix2cc(0:gs2Tot(0,pr),i1,i1,pr) = 0d0
    endif
    do i = 1,legs
      if (i.eq.i1) cycle
      call compute_squared_amplitude_cc (pr,i1,i)
      matrix2cc(0:gs2Tot(0,pr),i1,i1,pr) = &
        + matrix2cc(0:gs2Tot(0,pr),i1,i1,pr) &
        - matrix2cc(0:gs2Tot(0,pr),i1,i,pr)
    enddo
  else
    call compute_squared_amplitude_cc (pr,i1,i2)
  endif

  if (writeCor.ge.1) then
    call print_squared_amplitude_cc (pr,i1,i2)
    write(nx,'(1x,75("x"))')
    write(nx,*)
    write(nx,*)
  endif

  if (present(A2cc)) then
    j1 = newleg(i1,pr)
    j2 = newleg(i2,pr)
    A2cc = sum(matrix2cc(0:gs2Tot(0,pr),j1,j2,pr))
  endif

  end subroutine rescale_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_all_colour_correlations_rcl (npr,p,momenta_check)

  ! This subroutine computes the colour-correlated summed squared
  ! amplitudes for the process with process number "npr".
  ! "p" are the external momenta of the process.
  ! The colour correlation is done for all possible pairs coloured
  ! particles.
  ! Summed squared amplitude means:
  ! - summed over all computed powers of g_s
  ! - summed over colour and spins for outgoing particles and averaged
  !   for the incoming particles.
  ! If some particles were defined as polarized, no sum/average is
  ! performed on these particles.
  ! "momenta_check" tells whether the phase-space point is good
  ! (momenta_check=.true.) or bad (momenta_check=.false).

  integer,  intent(in)            :: npr
  real(dp), intent(in)            :: p(0:,:)
  logical,  intent(out), optional :: momenta_check

  integer :: pr,i,legs,i1,i2

  if (.not.processes_generated) then
    if (warnings(61).le.warning_limit) then
      warnings(61) = warnings(61) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 61: Call of compute_all_colour_correlations_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(61)
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
    if (warnings(62).le.warning_limit) then
      warnings(62) = warnings(62) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 62: compute_all_colour_correlations_rcl '
      write(nx,'(a,i3)') '           called with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(62)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if (size(p,2).ne.legs.or.size(p,1).ne.4) then
    if (warnings(63).le.warning_limit) then
      warnings(63) = warnings(63) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 63: compute_all_colour_correlations_rcl called ', &
                         'with wrong momenta'
      write(nx,*)
      call toomanywarnings(63)
    endif
    call istop (ifail,1)
  endif

  call set_momenta (pr,p)
  if (present(momenta_check)) momenta_check = momcheck

  if (writeCor.ge.1) then
    call print_process_and_momenta (pr)
    call print_rescaling
    call print_parameters_change
  endif

  call compute_amplitude (pr,'LO')

  call rescale_amplitude (pr,'LO')

  if (writeMat.ge.1) call print_amplitude (pr,'LO')

  do i1 = 1,legs
    do i2 = 1,legs
      call compute_squared_amplitude_cc (pr,i1,i2)
      if (writeCor.ge.1) call print_squared_amplitude_cc (pr,i1,i2)
    enddo
  enddo

  if (writeCor.ge.1) then
    call openOutput
    write(nx,'(1x,75("x"))')
    write(nx,*)
    write(nx,*)
  endif

  end subroutine compute_all_colour_correlations_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine rescale_all_colour_correlations_rcl (npr)

  ! This subroutine adjusts the results calculated by
  ! compute_all_colour_correlations_rcl for a new value of alpha_s,
  ! rescaling the structure-dressed helicity amplitudes and
  ! recomputing the colour-correlated summed squared amplitude for
  ! the process with process number "npr".

  integer,  intent(in) :: npr

  integer :: pr,i,legs,i1,i2

  if (.not.processes_generated) then
    if (warnings(70).le.warning_limit) then
      warnings(70) = warnings(70) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 70: Call of rescale_all_colour_correlations_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(70)
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
    if (warnings(71).le.warning_limit) then
      warnings(71) = warnings(71) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 71: rescale_all_colour_correlations_rcl '
      write(nx,'(a,i3)') '           called with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(71)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if (writeCor.ge.1) call print_rescaling
  if (writeCor.ge.1) call print_parameters_change

  call rescale_amplitude (pr,'LO')

  if (writeMat.ge.1) call print_amplitude (pr,'LO')

  do i1 = 1,legs
    if (.not.allocated(matrix2cc))  then
      allocate (matrix2cc(0:gs2Max,1:legsMax,1:legsMax,prTot))
      matrix2cc = 0d0
    else
      matrix2cc(0:gs2Tot(0,pr),i1,i1,pr) = 0d0
    endif
    do i2 = 1,legs
      if (i2.ne.i1) then
        call compute_squared_amplitude_cc (pr,i1,i2)
        if (writeCor.ge.1) call print_squared_amplitude_cc (pr,i1,i2)
        matrix2cc(0:gs2Tot(0,pr),i1,i1,pr) = &
          + matrix2cc(0:gs2Tot(0,pr),i1,i1,pr) &
          - matrix2cc(0:gs2Tot(0,pr),i1,i2,pr)
      endif
    enddo
    if (writeCor.ge.1) call print_squared_amplitude_cc (pr,i1,i1)
  enddo

  if (writeCor.ge.1) then
    call openOutput
    write(nx,'(1x,75("x"))')
    write(nx,*)
    write(nx,*)
  endif

  end subroutine rescale_all_colour_correlations_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_colour_correlation_rcl (npr,pow,i1,i2,A2cc)

  ! This subroutine extracts the computed value of the
  ! colour-correlated summed squared amplitude, with "pow" power of
  ! alpha_s, between particle with leg number "i1" and particle with
  ! leg number "i2", for the process with process number "npr".
  ! This routine must be called after compute_colour_correlation_rcl
  ! or rescale_colour_correlation_rcl (for process with process number
  ! "npr" for particles "i1" and "i2") or after
  ! compute_all_colour_correlations_rcl or
  ! rescale_all_colour_correlations_rcl for process with process
  ! number "npr".

  integer,  intent(in)  :: npr,pow,i1,i2
  real(dp), intent(out) :: A2cc

  integer :: pr,i,legs,j1,j2,gspower

  if (.not.processes_generated) then
    if (warnings(72).le.warning_limit) then
      warnings(72) = warnings(72) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 72: Call of get_colour_correlation_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(72)
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
    if (warnings(73).le.warning_limit) then
      warnings(73) = warnings(73) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 73: get_colour_correlation_rcl called '// &
                                    'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(73)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ((i1.le.0.or.i1.gt.legs).or.(i2.le.0.or.i2.gt.legs)) then
    if (warnings(74).le.warning_limit) then
      warnings(74) = warnings(74) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 74: get_colour_correlation_rcl called ', &
                            'with wrong indices'
      write(nx,*)
      call toomanywarnings(74)
    endif
    call istop (ifail,1)
  endif

  j1 = newleg(i1,pr)
  j2 = newleg(i2,pr)

  if (pow.lt.0.or.pow.gt.legs-2) then
    if (warnings(75).le.warning_limit) then
      warnings(75) = warnings(75) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 75: get_colour_correlation_rcl called ', &
                  '          with wrong alphas power'
      write(nx,*)
      call toomanywarnings(75)
    endif
    call istop (ifail,1)
  endif

  gspower = 2*pow

  A2cc = matrix2cc(gspower,j1,j2,pr)

  end subroutine get_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_spin_colour_correlation_rcl (npr,p,i1,i2,v, &
                                                  A2scc,momenta_check)

  ! This subroutine computes the colour-correlated summed squared
  ! amplitude, between particle with leg number "i1" and particle
  ! with leg number "i2", for the process with process number "npr",
  ! allowing spin-correlation if particle with leg number "i1"
  ! is a gluon.
  ! The spin-correlation is achieved by a user-defined polarization
  ! vector ("v") for the particle with leg number "i1", which
  ! substitutes the usual one.
  ! "p" are the external momenta of the process.
  ! "A2scc" is the spin-colour-correlated summed squared  LO
  ! amplitude.
  ! Summed squared amplitude means:
  ! - summed over all computed powers of g_s
  ! - summed over colour and spins for outgoing particles and averaged
  !   for the incoming particles.
  ! If some particles were defined as polarized, no sum/average is
  ! performed on these particles.
  ! "momenta_check" tells whether the phase-space point is good
  ! (momenta_check=.true.) or bad (momenta_check=.false).

  integer,     intent(in)            :: npr,i1,i2
  real(dp),    intent(in)            :: p(0:,:)
  complex(dp), intent(in)            :: v(0:)
  real(dp),    intent(out), optional :: A2scc
  logical,     intent(out), optional :: momenta_check

  integer :: pr,i,legs,j1,j2

  if (.not.processes_generated) then
    if (warnings(76).le.warning_limit) then
      warnings(76) = warnings(76) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 76: Call of compute_spin_colour_correlation_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(76)
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
    if (warnings(77).le.warning_limit) then
      warnings(77) = warnings(77) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 77: compute_spin_colour_correlation_rcl '
      write(nx,'(a,i3)') '           called with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(77)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ((cpar(par(i1,pr)).ne.'g') .or. &
      (i1.le.0.or.i1.gt.legs)   .or. &
      (i2.le.0.or.i2.gt.legs)        &
    ) then
    if (warnings(78).le.warning_limit) then
      warnings(78) = warnings(78) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 78: compute_spin_colour_correlation_rcl called ', &
                         'with wrong indices'
      write(nx,*)
      call toomanywarnings(78)
    endif
    call istop (ifail,1)
  endif

  if (size(p,2).ne.legs.or.size(p,1).ne.4) then
    if (warnings(79).le.warning_limit) then
      warnings(79) = warnings(79) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 79: compute_spin_colour_correlation_rcl called ', &
                         'with wrong momenta'
      write(nx,*)
      call toomanywarnings(79)
    endif
    call istop (ifail,1)
  endif

  call set_momenta (pr,p)
  if (present(momenta_check)) momenta_check = momcheck

  if (writeCor.ge.1) then
    call print_process_and_momenta (pr)
    call print_rescaling
    call print_parameters_change
  endif

  call compute_amplitude (pr,'LO')

  call rescale_amplitude (pr,'LO')

  if (writeMat.ge.1) call print_amplitude (pr,'LO')

  if (i2.eq.i1) then
    if (.not.allocated(matrix2scc)) then
      allocate (matrix2scc(0:gs2Max,1:legsMax,1:legsMax,prTot))
      matrix2scc = 0d0
    else
      matrix2scc(0:gs2Tot(0,pr),i1,i1,pr) = 0d0
    endif
    do i = 1,legs
      if (i.eq.i1) cycle
      call compute_squared_amplitude_scc (pr,i1,i,v)
      matrix2scc(0:gs2Tot(0,pr),i1,i1,pr) = &
        + matrix2scc(0:gs2Tot(0,pr),i1,i1,pr) &
        - matrix2scc(0:gs2Tot(0,pr),i1,i,pr)
    enddo
  else
    call compute_squared_amplitude_scc (pr,i1,i2,v)
  endif

  if (writeCor.ge.1) call print_squared_amplitude_scc (pr,i1,i2,v)

  if (present(A2scc)) then
    j1 = newleg(i1,pr)
    j2 = newleg(i2,pr)
    A2scc = sum(matrix2scc(0:gs2Tot(0,pr),j1,j2,pr))
  endif

  end subroutine compute_spin_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine rescale_spin_colour_correlation_rcl (npr,i1,i2,v,A2scc)

  ! This subroutine adjusts the results calculated by
  ! compute_spin_colour_correlation_rcl for a new value of alpha_s,
  ! rescaling the structure-dressed helicity amplitudes and
  ! recomputing the spin-colour-correlated summed squared amplitude
  ! for the process with process number "npr".
  ! "i1", "i2", "v" and "A2scc" are the same has for
  ! compute_spin_colour_correlation_rcl.

  integer,     intent(in)            :: npr,i1,i2
  complex(dp), intent(in)            :: v(0:)
  real(dp),    intent(out), optional :: A2scc

  integer :: pr,i,legs,j1,j2

  if (.not.processes_generated) then
    if (warnings(80).le.warning_limit) then
      warnings(80) = warnings(80) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 80: Call of rescale_spin_colour_correlation_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(80)
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
    if (warnings(81).le.warning_limit) then
      warnings(81) = warnings(81) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 81: rescale_spin_colour_correlation_rcl '
      write(nx,'(a,i3)') '           called with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(81)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ((cpar(par(i1,pr)).ne.'g') .or. &
      (i1.le.0.or.i1.gt.legs)   .or. &
      (i2.le.0.or.i2.gt.legs)        &
    ) then
    if (warnings(82).le.warning_limit) then
      warnings(82) = warnings(82) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 82: rescale_spin_colour_correlation_rcl called ', &
                         'with wrong indices'
      write(nx,*)
      call toomanywarnings(82)
    endif
    call istop (ifail,1)
  endif

  if (writeCor.ge.1) call print_rescaling
  if (writeCor.ge.1) call print_parameters_change

  call rescale_amplitude (pr,'LO')

  if (writeMat.ge.1) call print_amplitude (pr,'LO')

  if (i2.eq.i1) then
    if (.not.allocated(matrix2scc)) then
      allocate (matrix2scc(0:gs2Max,1:legsMax,1:legsMax,prTot))
      matrix2scc = 0d0
    else
      matrix2scc(0:gs2Tot(0,pr),i1,i1,pr) = 0d0
    endif
    do i = 1,legs
      if (i.eq.i1) cycle
      call compute_squared_amplitude_scc (pr,i1,i,v)
      matrix2scc(0:gs2Tot(0,pr),i1,i1,pr) = &
        + matrix2scc(0:gs2Tot(0,pr),i1,i1,pr) &
        - matrix2scc(0:gs2Tot(0,pr),i1,i,pr)
    enddo
  else
    call compute_squared_amplitude_scc (pr,i1,i2,v)
  endif

  if (writeCor.ge.1) call print_squared_amplitude_scc (pr,i1,i2,v)

  if (present(A2scc)) then
    j1 = newleg(i1,pr)
    j2 = newleg(i2,pr)
    A2scc = sum(matrix2scc(0:gs2Tot(0,pr),j1,j2,pr))
  endif

  end subroutine rescale_spin_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_spin_colour_correlation_rcl (npr,pow,i1,i2,A2scc)

  ! This subroutine extracts the computed value of the
  ! spin-colour-correlated summed squared amplitude with "pow" power
  ! of alpha_s, between particle with leg number "i1" and particle
  ! with leg number "i2", for the process with process number "npr".
  ! This routine must be called after
  ! compute_spin_colour_correlation_rcl or
  ! rescale_spin_colour_correlation_rcl for process with process
  ! number "npr" for particles "i1" and "i2" with the  desired
  ! polarization vestor.

  integer,  intent(in)  :: npr,pow,i1,i2
  real(dp), intent(out) :: A2scc

  integer               :: pr,i,legs,j1,j2,gspower

  if (.not.processes_generated) then
    if (warnings(83).le.warning_limit) then
      warnings(83) = warnings(83) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 83: Call of get_spin_colour_correlation_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(83)
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
    if (warnings(84).le.warning_limit) then
      warnings(84) = warnings(84) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 84: get_spin_colour_correlation_rcl called '
      write(nx,'(a,i3)') '           with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(84)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ((cpar(par(i1,pr)).ne.'g') .or. &
      (i1.le.0.or.i1.gt.legs)   .or. &
      (i2.le.0.or.i2.gt.legs)        &
     ) then
    if (warnings(85).le.warning_limit) then
      warnings(85) = warnings(85) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 85: get_spin_colour_correlation_rcl called ', &
                         'with wrong indices'
      write(nx,*)
      call toomanywarnings(85)
    endif
    call istop (ifail,1)
  endif

  j1 = newleg(i1,pr)
  j2 = newleg(i2,pr)

  if (pow.lt.0.or.pow.gt.legs-2) then
    if (warnings(86).le.warning_limit) then
      warnings(86) = warnings(86) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 86: get_spin_colour_correlation_rcl called ', &
                         'with wrong alphas power'
      write(nx,*)
      call toomanywarnings(86)
    endif
    call istop (ifail,1)
  endif

  gspower = 2*pow

  A2scc = matrix2scc(gspower,j1,j2,pr)

  end subroutine get_spin_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_spin_correlation_rcl (npr,p,j,v,A2sc, &
                                           momenta_check)

  ! This subroutine computes the summed squared amplitude in the
  ! case where the polarization vector of one of the external photons
  ! or gluons is given by the user.
  ! "p" are the external momenta of the process.
  ! "j" is an integer indicating the leg number for the photon or
  ! gluon with the user-defined polarization vector.
  ! "v" is the user-defined polarization vector.
  ! "A2sc" is the spin-colour-correlated summed squared LO
  ! amplitude.
  ! Summed squared amplitude means:
  ! - summed over all computed powers of g_s
  ! - summed over colour and spins for outgoing particles and averaged
  !   for the incoming particles.
  ! If some particles were defined as polarized, no sum/average is
  ! performed on these particles.
  ! "momenta_check" tells whether the phase-space point is good
  ! (momenta_check=.true.) or bad (momenta_check=.false).

  integer,     intent(in)            :: npr,j
  real(dp),    intent(in)            :: p(0:,:)
  complex(dp), intent(in)            :: v(0:)
  real(dp),    intent(out), optional :: A2sc
  logical,     intent(out), optional :: momenta_check

  integer :: pr,i,legs

  if (.not.processes_generated) then
    if (warnings(87).le.warning_limit) then
      warnings(87) = warnings(87) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 87: Call of compute_spin_correlation_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(87)
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
    if (warnings(88).le.warning_limit) then
      warnings(88) = warnings(88) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 88: compute_spin_correlation_rcl called '
      write(nx,'(a,i3)') '           with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(88)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ( (cpar(par(j,pr)).ne.'A'.and.cpar(par(j,pr)).ne.'g') .or. &
       j.le.0.or.j.gt.legs                                      ) then
    if (warnings(89).le.warning_limit) then
      warnings(89) = warnings(89) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 89: compute_spin_correlation_rcl called ', &
                         'with wrong index'
      write(nx,*)
      call toomanywarnings(89)
    endif
    call istop (ifail,1)
  endif

  if (size(p,2).ne.legs.or.size(p,1).ne.4) then
    if (warnings(90).le.warning_limit) then
      warnings(90) = warnings(90) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 90: compute_spin_correlation_rcl called ', &
                         'with wrong momenta'
      write(nx,*)
      call toomanywarnings(90)
    endif
    call istop (ifail,1)
  endif

  call set_momenta (pr,p)
  if (present(momenta_check)) momenta_check = momcheck

  if (writeCor.ge.1) then
    call print_process_and_momenta (pr)
    call print_rescaling
    call print_parameters_change
  endif

  call compute_amplitude (pr,'LO')

  call rescale_amplitude (pr,'LO')

  if (writeMat.ge.1) call print_amplitude (pr,'LO')

  call compute_squared_amplitude_sc (pr,j,v)

  if (writeCor.ge.1) call print_squared_amplitude_sc (pr,j,v)

  if (present(A2sc)) A2sc = sum(matrix2sc(0:gs2Tot(0,pr),pr))

  end subroutine compute_spin_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine rescale_spin_correlation_rcl (npr,j,v,A2sc)

  ! This subroutine adjusts the results calculated by
  ! compute_spin_correlation_rcl for a new value of alpha_s,
  ! rescaling the structure-dressed helicity amplitudes and
  ! recomputing the spin-correlated summed squared amplitude for
  ! the process with process number "npr".
  ! "j", "v" and "A2sc" are the same has for
  ! compute_spin_correlation_rcl.

  integer,     intent(in)            :: npr,j
  complex(dp), intent(in)            :: v(0:)
  real(dp),    intent(out), optional :: A2sc

  integer :: pr,i,legs

  if (.not.processes_generated) then
    if (warnings(91).le.warning_limit) then
      warnings(91) = warnings(91) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 91: Call of rescale_spin_correlation_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(91)
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
    if (warnings(92).le.warning_limit) then
      warnings(92) = warnings(92) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 92: rescale_spin_correlation_rcl called '
      write(nx,'(a,i3)') '           with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(92)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ( (cpar(par(j,pr)).ne.'A'.and.cpar(par(j,pr)).ne.'g') .or. &
       j.le.0.or.j.gt.legs                                      ) then
    if (warnings(93).le.warning_limit) then
      warnings(93) = warnings(93) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 93: rescale_spin_correlation_rcl called ', &
                         'with wrong index'
      write(nx,*)
      call toomanywarnings(93)
    endif
    call istop (ifail,1)
  endif

  if (writeCor.ge.1) call print_rescaling
  if (writeCor.ge.1) call print_parameters_change

  call rescale_amplitude (pr,'LO')

  if (writeMat.ge.1) call print_amplitude (pr,'LO')

  call compute_squared_amplitude_sc (pr,j,v)

  if (writeCor.ge.1) call print_squared_amplitude_sc (pr,j,v)

  if (present(A2sc)) A2sc = sum(matrix2sc(0:gs2Tot(0,pr),pr))

  end subroutine rescale_spin_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_spin_correlation_rcl (npr,pow,A2sc)

  ! This subroutine extracts the computed value of the
  ! spin-correlated summed squared amplitude with "pow" power
  ! of alpha_s for the process with process number "npr".
  ! This routine must be called after compute_spin_correlation_rcl
  ! or rescale_spin_correlation_rcl for process with process number
  ! "npr" with the desired polarization vestor.

  integer,  intent(in)  :: npr,pow
  real(dp), intent(out) :: A2sc

  integer :: pr,i,legs,gspower

  if (.not.processes_generated) then
    if (warnings(94).le.warning_limit) then
      warnings(94) = warnings(94) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 94: Call of get_spin_correlation_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(94)
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
    if (warnings(95).le.warning_limit) then
      warnings(95) = warnings(95) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 95: get_spin_correlation_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(95)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if (pow.lt.0.or.pow.gt.legs-2) then
    if (warnings(96).le.warning_limit) then
      warnings(96) = warnings(96) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 96: get_spin_correlation_rcl called ', &
                  'with wrong alphas power'
      write(nx,*)
      call toomanywarnings(96)
    endif
    call istop (ifail,1)
  endif

  gspower = 2*pow

  A2sc = matrix2sc(gspower,pr)

  end subroutine get_spin_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_spin_correlation_matrix_rcl (npr,p,j,A2scm,momenta_check)

  ! This subroutine computes the spin-correlation matrix which is the
  ! summed squared amplitude where the polarization vector of a photon
  ! or a gluon has not been inserted (leaving two open Lorentz
  ! indices).
  ! "p" are the external momenta of the process.
  ! "j" is an integer indicating the leg number for the photon or
  ! gluon with wioth no polarization vector.
  ! "A2scm" is the spin_correlation_matrix, summed over all computed
  ! powers of g_s, summed over colour and spins for outgoing particles
  ! and averaged for the incoming particles.
  ! If some particles were defined as polarized, no sum/average is
  ! performed on these particles.
  ! "momenta_check" tells whether the phase-space point is good
  ! (momenta_check=.true.) or bad (momenta_check=.false).

  integer,                    intent(in)  :: npr,j
  real(dp),                   intent(in)  :: p(0:,:)
  real(dp), optional,         intent(out) :: A2scm(0:3,0:3)
  logical, optional,          intent(out) :: momenta_check

  integer      :: pr,legs,i,gs

  if (.not.processes_generated) then
    if (warnings(97).le.warning_limit) then
      warnings(97) = warnings(97) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 97: Call of compute_spin_correlation_matrix_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(97)
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
    if (warnings(98).le.warning_limit) then
      warnings(98) = warnings(98) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 98: compute_spin_correlation_matrix_rcl called '
      write(nx,'(a,i3)') '           with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(98)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ( (cpar(par(j,pr)).ne.'A'.and.cpar(par(j,pr)).ne.'g') .or. &
       j.le.0.or.j.gt.legs                                      ) then
    if (warnings(99).le.warning_limit) then
      warnings(99) = warnings(99) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 99: compute_spin_correlation_matrix_rcl called ', &
                  'with wrong index'
      write(nx,*)
      call toomanywarnings(99)
    endif
    call istop (ifail,1)
  endif

  if (size(p,2).ne.legs.or.size(p,1).ne.4) then
    if (warnings(100).le.warning_limit) then
      warnings(100) = warnings(100) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 100: compute_spin_correlation_matrix_rcl called ', &
                  'with wrong momenta'
      write(nx,*)
      call toomanywarnings(100)
    endif
    call istop (ifail,1)
  endif

  call set_momenta (pr,p)
  if (present(momenta_check)) momenta_check = momcheck

  if (writeCor.ge.1) then
    call print_process_and_momenta (pr)
    call print_rescaling
    call print_parameters_change
  endif

  call compute_amplitude (pr,'LO')

  call rescale_amplitude (pr,'LO')

  if (writeMat.ge.1) call print_amplitude (pr,'LO')

  call compute_squared_amplitude_scm (pr,j)

  if (present(A2scm)) then
    A2scm(:,:) = 0
    do gs = 0, gs2Tot(0,pr)
      A2scm(:,:) = A2scm(:,:) + sum(matrix2scm(:,:,gs,pr))
    end do
  end if

  end subroutine compute_spin_correlation_matrix_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine rescale_spin_correlation_matrix_rcl (npr,j,A2scm)

  ! This subroutine adjusts the results calculated by
  ! compute_spin_correlation_matrix_rcl for a new value of alpha_s,
  ! rescaling the structure-dressed helicity amplitudes and
  ! recomputing the spin-correlated summed squared amplitude for
  ! the process with process number "npr".
  ! "j" and "A2scm" are the same has for
  ! compute_spin_correlation_matrix_rcl.

  integer,     intent(in)            :: npr,j
  real(dp),    intent(out), optional :: A2scm(0:3,0:3)

  integer :: pr,i,legs,gs

  if (.not.processes_generated) then
    if (warnings(101).le.warning_limit) then
      warnings(101) = warnings(101) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 101: Call of rescale_spin_correlation_matrix_rcl not allowed:'
      write(nx,*) '           Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(101)
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
    if (warnings(102).le.warning_limit) then
      warnings(102) = warnings(102) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 102: rescale_spin_correlation_matrix_rcl called '
      write(nx,'(a,i3)') '            with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(102)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ( (cpar(par(j,pr)).ne.'A'.and.cpar(par(j,pr)).ne.'g') .or. &
       j.le.0.or.j.gt.legs                                      ) then
    if (warnings(103).le.warning_limit) then
      warnings(103) = warnings(103) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 103: rescale_spin_correlation_matrix_rcl called ', &
                  'with wrong index'
      write(nx,*)
      call toomanywarnings(103)
    endif
    call istop (ifail,1)
  endif

  if (writeCor.ge.1) call print_rescaling
  if (writeCor.ge.1) call print_parameters_change

  call rescale_amplitude (pr,'LO')

  if (writeMat.ge.1) call print_amplitude (pr,'LO')

  call compute_squared_amplitude_scm (pr,j)

  if (present(A2scm)) then
    A2scm(:,:) = 0
    do gs = 0, gs2Tot(0,pr)
      A2scm(:,:) = A2scm(:,:) + sum(matrix2scm(:,:,gs,pr))
    end do
  end if

  end subroutine rescale_spin_correlation_matrix_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_spin_correlation_matrix_rcl (npr,pow,A2scm)

  ! This subroutine extracts the computed value of the
  ! spin-correlation matrix with "pow" power of alpha_s for the
  ! process with process number "npr".
  ! This routine must be called after compute_spin_correlation_rcl
  ! or rescale_spin_correlation_rcl for process with process number
  ! "npr" for the desired photon or gluon "j".

  integer,  intent(in)  :: npr,pow
  real(dp), intent(out) :: A2scm(0:3,0:3)

  integer :: pr,i,legs,gspower

  if (.not.processes_generated) then
    if (warnings(104).le.warning_limit) then
      warnings(104) = warnings(104) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 104: Call of get_spin_correlation_matrix_rcl not allowed:'
      write(nx,*) '           Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(104)
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
    if (warnings(105).le.warning_limit) then
      warnings(105) = warnings(105) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 105: get_spin_correlation_matrix_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(105)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if (pow.lt.0.or.pow.gt.legs-2) then
    if (warnings(106).le.warning_limit) then
      warnings(106) = warnings(106) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 106: get_spin_correlation_matrix_rcl called ', &
                  'with wrong alphas power'
      write(nx,*)
      call toomanywarnings(106)
    endif
    call istop (ifail,1)
  endif

  gspower = 2*pow

  A2scm(:,:) = matrix2scm(:,:,gspower,pr)

  end subroutine get_spin_correlation_matrix_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_int_colour_correlation_rcl (npr,p,i1,i2, &
                                                 A2ccint,     &
                                                 momenta_check)

  ! This subroutine computes the colour-correlated summed LO-NLO
  ! interference, between particle with leg number "i1" and particle
  ! with leg number "i2", for the process with process number "npr".
  ! "p" are the external momenta of the process.
  ! "A2ccint" is the colour-correlated summed LO-NLO interference.
  ! Summed interference means:
  ! - summed over all computed powers of g_s
  ! - summed over colour and spins for outgoing particles and averaged
  !   for the incoming particles.
  ! If some particles were defined as polarized, no sum/average is
  ! performed on these particles.
  ! "momenta_check" tells whether the phase-space point is good
  ! (momenta_check=.true.) or bad (momenta_check=.false).

  integer,  intent(in)            :: npr,i1,i2
  real(dp), intent(in)            :: p(0:,:)
  real(dp), intent(out), optional :: A2ccint
  logical,  intent(out), optional :: momenta_check

  integer :: pr,i,legs,j1,j2

  if (.not.processes_generated) then
    if (warnings(46).le.warning_limit) then
      warnings(46) = warnings(46) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 46: Call of compute_int_colour_correlation_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(46)
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
    if (warnings(47).le.warning_limit) then
      warnings(47) = warnings(47) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 47: compute_int_colour_correlation_rcl called '
      write(nx,'(a,i3)') '           with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(47)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ((i1.le.0.or.i1.gt.legs).or.(i2.le.0.or.i2.gt.legs)) then
    if (warnings(48).le.warning_limit) then
      warnings(48) = warnings(48) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 48: compute_int_colour_correlation_rcl called ', &
                         'with wrong indices'
      write(nx,*)
      call toomanywarnings(48)
    endif
    call istop (ifail,1)
  endif

  if (size(p,2).ne.legs.or.size(p,1).ne.4) then
    if (warnings(49).le.warning_limit) then
      warnings(49) = warnings(49) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 49: compute_int_colour_correlation_rcl called ', &
                         'with wrong momenta'
      write(nx,*)
      call toomanywarnings(49)
    endif
    call istop (ifail,1)
  endif

  call set_momenta (pr,p)
  if (present(momenta_check)) momenta_check = momcheck

  if (writeCor.ge.1) then
    call print_process_and_momenta (pr)
    call print_rescaling
    call print_parameters_change
  endif

  call compute_amplitude (pr,'NLO')

  call rescale_amplitude (pr,'NLO')

  if (writeMat.ge.1) call print_amplitude (pr,'NLO')

  call compute_squared_amplitude_cc_int (pr,i1,i2)

  if (writeCor.ge.1) then
    call print_squared_amplitude_cc_int (pr,i1,i2)
    write(nx,'(1x,75("x"))')
    write(nx,*)
    write(nx,*)
  endif

  if (present(A2ccint)) then
    j1 = newleg(i1,pr)
    j2 = newleg(i2,pr)
    A2ccint = sum(matrix2ccint(0:gs2Tot(0,pr),j1,j2,pr))
  endif

  end subroutine compute_int_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine rescale_int_colour_correlation_rcl (npr,i1,i2,A2ccint)

  ! This subroutine adjusts the results calculated by
  ! compute_int_colour_correlation_rcl for a new value of alpha_s,
  ! rescaling the structure-dressed helicity amplitudes and
  ! recomputing the colour-correlated summed LO-NLO interference for
  ! the process with process number "npr".
  ! "i1", "i2" and "A2ccint" are the same has for
  ! compute_int_colour_correlation_rcl.

  integer,  intent(in)            :: npr,i1,i2
  real(dp), intent(out), optional :: A2ccint

  integer :: pr,i,legs,j1,j2

  if (.not.processes_generated) then
    if (warnings(54).le.warning_limit) then
      warnings(54) = warnings(54) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 54: Call of rescale_int_colour_correlation_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(54)
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
    if (warnings(55).le.warning_limit) then
      warnings(55) = warnings(55) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 55: rescale_int_colour_correlation_rcl called '
      write(nx,'(a,i3)') '           with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(55)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ((i1.le.0.or.i1.gt.legs).or.(i2.le.0.or.i2.gt.legs)) then
    if (warnings(56).le.warning_limit) then
      warnings(56) = warnings(56) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 56: rescale_int_colour_correlation_rcl called ', &
                         'with wrong indices'
      write(nx,*)
      call toomanywarnings(56)
    endif
    call istop (ifail,1)
  endif

  if (writeCor.ge.1) call print_rescaling
  if (writeCor.ge.1) call print_parameters_change

  call rescale_amplitude (pr,'NLO')

  if (writeMat.ge.1) call print_amplitude (pr,'NLO')

  call compute_squared_amplitude_cc_int (pr,i1,i2)

  if (writeCor.ge.1) then
    call print_squared_amplitude_cc_int (pr,i1,i2)
    write(nx,'(1x,75("x"))')
    write(nx,*)
    write(nx,*)
  endif

  if (present(A2ccint)) then
    j1 = newleg(i1,pr)
    j2 = newleg(i2,pr)
    A2ccint = sum(matrix2ccint(0:gs2Tot(1,pr),j1,j2,pr))
  endif

  end subroutine rescale_int_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_all_int_colour_correlations_rcl (npr,p, &
                                                      momenta_check)

  ! This subroutine computes the colour-correlated summed LO-NLO
  ! interference for the process with process number "npr".
  ! "p" are the external momenta of the process.
  ! The colour correlation is done for all possible pairs coloured
  ! particles.
  ! Summed squared amplitude means:
  ! - summed over all computed powers of g_s
  ! - summed over colour and spins for outgoing particles and averaged
  !   for the incoming particles.
  ! If some particles were defined as polarized, no sum/average is
  ! performed on these particles.
  ! "momenta_check" tells whether the phase-space point is good
  ! (momenta_check=.true.) or bad (momenta_check=.false).

  integer,  intent(in)            :: npr
  real(dp), intent(in)            :: p(0:,:)
  logical,  intent(out), optional :: momenta_check

  integer :: pr,i,legs,i1,i2

  if (.not.processes_generated) then
    if (warnings(64).le.warning_limit) then
      warnings(64) = warnings(64) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 64: Call of compute_all_int_colour_correlations_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(64)
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
    if (warnings(65).le.warning_limit) then
      warnings(65) = warnings(65) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 65: compute_all_int_colour_correlations_rcl '
      write(nx,'(a,i3)') '           called with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(65)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if (size(p,2).ne.legs.or.size(p,1).ne.4) then
    if (warnings(66).le.warning_limit) then
      warnings(66) = warnings(66) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 66: compute_all_int_colour_correlations_rcl called ', &
                         'with wrong momenta'
      write(nx,*)
      call toomanywarnings(66)
    endif
    call istop (ifail,1)
  endif

  call set_momenta (pr,p)
  if (present(momenta_check)) momenta_check = momcheck

  if (writeCor.ge.1) then
    call print_process_and_momenta (pr)
    call print_rescaling
    call print_parameters_change
  endif

  call compute_amplitude (pr,'NLO')

  call rescale_amplitude (pr,'NLO')

  if (writeMat.ge.1) call print_amplitude (pr,'NLO')

  do i1 = 1,legs
  do i2 = 1,legs
    call compute_squared_amplitude_cc_int (pr,i1,i2)
    if (writeCor.ge.1) call print_squared_amplitude_cc_int (pr,i1,i2)
  enddo
  enddo

  if (writeCor.ge.1) then
    call openOutput
    write(nx,'(1x,75("x"))')
    write(nx,*)
    write(nx,*)
  endif

  end subroutine compute_all_int_colour_correlations_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine rescale_all_int_colour_correlations_rcl (npr)

  ! This subroutine adjusts the results calculated by
  ! compute_all_int_colour_correlations_rcl for a new value of
  ! alpha_s, rescaling the structure-dressed helicity amplitudes and
  ! recomputing the colour-correlated summed LO-NLO interference for
  ! the process with process number "npr".

  integer,  intent(in) :: npr

  integer :: pr,i,legs,i1,i2

  if (.not.processes_generated) then
    if (warnings(170).le.warning_limit) then
      warnings(170) = warnings(170) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 170: Call of rescale_all_int_colour_correlations_rcl not allowed:'
      write(nx,*) '           Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(170)
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
    if (warnings(171).le.warning_limit) then
      warnings(171) = warnings(171) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 171: rescale_all_int_colour_correlations_rcl '
      write(nx,'(a,i3)') '            called with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(171)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if (writeCor.ge.1) call print_rescaling
  if (writeCor.ge.1) call print_parameters_change

  call rescale_amplitude (pr,'NLO')

  if (writeMat.ge.1) call print_amplitude (pr,'NLO')

  do i1 = 1,legs
    if (.not.allocated(matrix2ccint)) then
      allocate (matrix2ccint(0:gs2Max,1:legsMax,1:legsMax,prTot))
      matrix2ccint = 0d0
    else
      matrix2ccint(0:gs2Tot(1,pr),i1,i1,pr) = 0d0
    endif
    do i2 = 1,legs
      if (i2.ne.i1) then
        call compute_squared_amplitude_cc_int (pr,i1,i2)
        if (writeCor.ge.1) call print_squared_amplitude_cc_int (pr,i1,i2)
        matrix2ccint(0:gs2Tot(1,pr),i1,i1,pr) = &
          + matrix2ccint(0:gs2Tot(1,pr),i1,i1,pr) &
          - matrix2ccint(0:gs2Tot(1,pr),i1,i2,pr)
      endif
    enddo
    if (writeCor.ge.1) call print_squared_amplitude_cc_int (pr,i1,i1)
  enddo

  if (writeCor.ge.1) then
    call openOutput
    write(nx,'(1x,75("x"))')
    write(nx,*)
    write(nx,*)
  endif

  end subroutine rescale_all_int_colour_correlations_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_int_colour_correlation_rcl (npr,pow,i1,i2,A2ccint)

  ! This subroutine extracts the computed value of the
  ! colour-correlated summed LO-NLO interference, with "pow" power of
  ! alpha_s, between particle with leg number "i1" and particle with
  ! leg number "i2", for the process with process number "npr".
  ! This routine must be called after
  ! compute_int_colour_correlation_rcl or
  ! rescale_int_colour_correlation_rcl (for process with process
  ! number "npr" for particles "i1" and "i2") or after
  ! compute_all_int_colour_correlations_rcl or
  ! rescale_all_int_colour_correlations_rcl for process with process
  ! number "npr".

  integer,  intent(in)  :: npr,pow,i1,i2
  real(dp), intent(out) :: A2ccint

  integer :: pr,i,legs,j1,j2,gspower

  if (.not.processes_generated) then
    if (warnings(50).le.warning_limit) then
      warnings(50) = warnings(50) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 50: Call of get_int_colour_correlation_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(50)
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
    if (warnings(51).le.warning_limit) then
      warnings(51) = warnings(51) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 51: get_int_colour_correlation_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(51)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ((i1.le.0.or.i1.gt.legs).or.(i2.le.0.or.i2.gt.legs)) then
    if (warnings(52).le.warning_limit) then
      warnings(52) = warnings(52) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 52: get_int_colour_correlation_rcl called ', &
                         'with wrong indices'
      write(nx,*)
      call toomanywarnings(52)
    endif
    call istop (ifail,1)
  endif

  j1 = newleg(i1,pr)
  j2 = newleg(i2,pr)

  if (pow.gt.legs-1) then
    if (warnings(53).le.warning_limit) then
      warnings(53) = warnings(53) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 53: get_squared_amplitude_rcl called', &
                        ' with wrong alphas power'
      write(nx,*)
      call toomanywarnings(53)
    endif
    call istop (ifail,1)
  endif

  gspower = 2*pow

  A2ccint = matrix2ccint(gspower,j1,j2,pr)

  end subroutine get_int_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine compute_nlo_colour_correlation_rcl (npr,p,i1,i2, &
                                                 A2ccnlo,     &
                                                 momenta_check)

  ! This subroutine computes the colour-correlated summed squared
  ! NLO amplitude, between particle with leg number "i1" and particle
  ! with leg number "i2", for the process with process number "npr".
  ! "p" are the external momenta of the process.
  ! "A2ccnlo" is the colour-correlated summed squared NLO
  ! amplitude.
  ! Summed squared amplitude means:
  ! - summed over all computed powers of g_s
  ! - summed over colour and spins for outgoing particles and averaged
  !   for the incoming particles.
  ! If some particles were defined as polarized, no sum/average is
  ! performed on these particles.
  ! "momenta_check" tells whether the phase-space point is good
  ! (momenta_check=.true.) or bad (momenta_check=.false).

  integer,  intent(in)            :: npr,i1,i2
  real(dp), intent(in)            :: p(0:,:)
  real(dp), intent(out), optional :: A2ccnlo
  logical,  intent(out), optional :: momenta_check

  integer :: pr,i,legs,j1,j2

  if (.not.processes_generated) then
    if (warnings(57).le.warning_limit) then
      warnings(57) = warnings(57) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 57: Call of compute_nlo_colour_correlation_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(57)
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
    if (warnings(58).le.warning_limit) then
      warnings(58) = warnings(58) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 58: compute_nlo_colour_correlation_rcl called '
      write(nx,'(a,i3)') '        with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(58)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ((i1.le.0.or.i1.gt.legs).or.(i2.le.0.or.i2.gt.legs)) then
    if (warnings(59).le.warning_limit) then
      warnings(59) = warnings(59) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 59: compute_nlo_colour_correlation_rcl called ', &
                         'with wrong indices'
      write(nx,*)
      call toomanywarnings(59)
    endif
    call istop (ifail,1)
  endif

  if (size(p,2).ne.legs.or.size(p,1).ne.4) then
    if (warnings(60).le.warning_limit) then
      warnings(60) = warnings(60) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 60: compute_nlo_colour_correlation_rcl called ', &
                         'with wrong momenta'
      write(nx,*)
      call toomanywarnings(60)
    endif
    call istop (ifail,1)
  endif

  call set_momenta (pr,p)
  if (present(momenta_check)) momenta_check = momcheck

  if (writeCor.ge.1) then
    call print_process_and_momenta (pr)
    call print_rescaling
    call print_parameters_change
  endif

  call compute_amplitude (pr,'NLO')

  call rescale_amplitude (pr,'NLO')

  if (writeMat.ge.1) call print_amplitude (pr,'NLO')

  call compute_squared_amplitude_cc_nlo (pr,i1,i2)

!  if (writeCor.ge.1) then
!    call print_squared_amplitude_cc_nlo (pr,i1,i2)
!    write(nx,'(1x,75("x"))')
!    write(nx,*)
!    write(nx,*)
!  endif

  if (present(A2ccnlo)) then
    j1 = newleg(i1,pr)
    j2 = newleg(i2,pr)
    A2ccnlo = sum(matrix2ccnlo(0:gs2Tot(1,pr),j1,j2,pr))
  endif

  end subroutine compute_nlo_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine rescale_nlo_colour_correlation_rcl (npr,i1,i2,A2ccnlo)

  ! This subroutine adjusts the results calculated by
  ! compute_nlo_colour_correlation_rcl for a new value of alpha_s,
  ! rescaling the structure-dressed helicity amplitudes and
  ! recomputing the colour-correlated summed squared NLO amplitude
  ! for the process with process number "npr".
  ! "i1", "i2" and "A2ccnlo" are the same has for
  ! compute_nlo_colour_correlation_rcl.

  integer,  intent(in)            :: npr,i1,i2
  real(dp), intent(out), optional :: A2ccnlo

  integer :: pr,i,legs,j1,j2

  if (.not.processes_generated) then
    if (warnings(154).le.warning_limit) then
      warnings(154) = warnings(154) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 154: Call of rescale_nlo_colour_correlation_rcl not allowed:'
      write(nx,*) '           Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(154)
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
    if (warnings(155).le.warning_limit) then
      warnings(155) = warnings(155) + 1
      call openOutput
      write(nx,*)
      write(nx,*)         'ERROR 155: rescale_nlo_colour_correlation_rcl called '
      write(nx,'(a,i3)') '            with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(155)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ((i1.le.0.or.i1.gt.legs).or.(i2.le.0.or.i2.gt.legs)) then
    if (warnings(156).le.warning_limit) then
      warnings(156) = warnings(156) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 156: rescale_nlo_colour_correlation_rcl called ', &
                             'with wrong indices'
      write(nx,*)
      call toomanywarnings(156)
    endif
    call istop (ifail,1)
  endif

  if (writeCor.ge.1) call print_rescaling
  if (writeCor.ge.1) call print_parameters_change

  call rescale_amplitude (pr,'NLO')

  if (writeMat.ge.1) call print_amplitude (pr,'NLO')

  call compute_squared_amplitude_cc_nlo (pr,i1,i2)

!  if (writeCor.ge.1) then
!    call print_squared_amplitude_cc_nlo (pr,i1,i2)
!    write(nx,'(1x,75("x"))')
!    write(nx,*)
!    write(nx,*)
!  endif

  if (present(A2ccnlo)) then
    j1 = newleg(i1,pr)
    j2 = newleg(i2,pr)
    A2ccnlo = sum(matrix2ccnlo(0:gs2Tot(1,pr),j1,j2,pr))
  endif

  end subroutine rescale_nlo_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_nlo_colour_correlation_rcl (npr,pow,i1,i2,A2ccnlo)

  ! This subroutine extracts the computed value of the
  ! colour-correlated summed squared NLO amplitude, with "pow" power
  ! of alpha_s, between particle with leg number "i1" and particle
  ! with leg number "i2", for the process with process number "npr".
  ! This routine must be called after
  ! compute_nlo_colour_correlation_rcl or
  ! rescale_nlo_colour_correlation_rcl (for process with process
  ! number "npr" for particles "i1" and "i2").

  integer,  intent(in)  :: npr,pow,i1,i2
  real(dp), intent(out) :: A2ccnlo

  integer :: pr,i,legs,j1,j2,gspower

  if (.not.processes_generated) then
    if (warnings(150).le.warning_limit) then
      warnings(150) = warnings(150) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 150: Call of get_nlo_colour_correlation_rcl not allowed:'
      write(nx,*) '           Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(150)
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
    if (warnings(151).le.warning_limit) then
      warnings(151) = warnings(151) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 151: get_nlo_colour_correlation_rcl called '// &
                                     'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(151)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if ((i1.le.0.or.i1.gt.legs).or.(i2.le.0.or.i2.gt.legs)) then
    if (warnings(152).le.warning_limit) then
      warnings(152) = warnings(152) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 152: get_nlo_colour_correlation_rcl called ', &
                             'with wrong indices'
      write(nx,*)
      call toomanywarnings(152)
    endif
    call istop (ifail,1)
  endif

  j1 = newleg(i1,pr)
  j2 = newleg(i2,pr)

  if (pow.gt.legs-1) then
    if (warnings(153).le.warning_limit) then
      warnings(153) = warnings(153) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 153: get_squared_amplitude_rcl called', &
                            ' with wrong alphas power'
      write(nx,*)
      call toomanywarnings(153)
    endif
    call istop (ifail,1)
  endif

  gspower = 2*pow

  A2ccnlo = matrix2ccnlo(gspower,j1,j2,pr)

  end subroutine get_nlo_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_momenta_rcl (npr,p)

  ! This subroutine extracts the stored momenta of the process with
  ! process number "npr"

  ! This routine must be called after compute_process_rcl,
  ! compute_colour_correlation_rcl,
  ! compute_all_colour_correlations_rcl,
  ! compute_spin_colour_correlation_rcl or
  ! compute_spin_correlation_rcl for process with process number
  ! "npr".

  integer,  intent(in)  :: npr
  real(dp), intent(out) :: p(0:,:)

  integer :: pr,i,legs

  if (.not.processes_generated) then
    if (warnings(107).le.warning_limit) then
      warnings(107) = warnings(107) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 107: Call of get_momenta_rcl not allowed:'
      write(nx,*) '           Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(107)
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
    if (warnings(108).le.warning_limit) then
      warnings(108) = warnings(108) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 108: get_momenta_rcl called '// &
                               'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(108)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  if (size(p,2).ne.legs.or.size(p,1).ne.4) then
    if (warnings(109).le.warning_limit) then
      warnings(109) = warnings(109) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 109: get_momenta_rcl called with wrong momenta'
      write(nx,*)
      call toomanywarnings(109)
    endif
    call istop (ifail,1)
  endif
  p(0:3,1:legs) = momenta(0:3,1:legs)

  end subroutine get_momenta_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_colour_configurations_rcl (npr,cols)

  ! This subroutine extracts all colour configurations of the
  ! process with process number "npr" and write it in the output
  ! variable "cols(1:l,1:csTot)" (csTot is the total number of
  ! colour configuations of the process "npr" and l is its total
  ! number of external legs).
  ! "cols(1:l,n)" describes the n^th colour structure and is a vector
  ! of type integer and length l, where each position in the vector
  ! corresponds to one of the l external particles of process "npr"
  ! (ordered as in the process definition).
  ! For colourless particles, incoming quarks and outgoing
  ! anti-quarks, the corresponding entry in "cols" is 0.
  ! For all other particles (gluons, outgoing quarks and incoming
  ! anti-quarks), the entries are permutations, without repetition,
  ! of the positions of the gluons, incoming quarks or outgoing
  ! anti-quarks in the process definition.
  ! "cols(1:l,n)" can be used as the input variable "colour" in the
  ! subroutine get_amplitude_rcl.

  integer, intent(in)               :: npr
  integer, intent(out), allocatable :: cols(:,:)

  integer :: pr,i,j,legs

  if (.not.processes_generated) then
    if (warnings(18).le.warning_limit) then
      warnings(18) = warnings(18) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 18: Call of get_colour_configurations_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(18)
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
    if (warnings(19).le.warning_limit) then
      warnings(19) = warnings(19) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 19: get_colour_configurations_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(19)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  allocate (cols(1:legs,1:csTot(pr)))

  do i = 1, csTot(pr)
    do j = 1,legs
      if (csIq(j,i,pr).ne.0) then
        cols(oldleg(j,pr),i) = oldleg(csIq(j,i,pr),pr)
      else
        cols(oldleg(j,pr),i) = 0
      endif
    enddo
  enddo

  end subroutine get_colour_configurations_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_helicity_configurations_rcl (npr,hels)

  ! This subroutine extracts all helicity configurations of the
  ! process with process number "npr" and write it in the output
  ! variable "hels(1:l,1:cfTot)" (cfTot is the total number of
  ! helicity configuations of the process "npr" and l is its total
  ! number of external legs).
  ! "hels(1:l,n)" describes the n^th helicity configuration and is a
  ! vector of type integer and length l, where each position in the
  ! vectorcorresponds to one of the l external particles of process
  ! "npr" (ordered as in the process definition).
  ! Its entries are the helicities of the particles:
  !   left-handed  fermions/antifermions -> -1
  !   right-handed fermions/antifermions -> +1
  !   logitudinal vector bosons          ->  0
  !   transverse vector bosons           -> -1, +1
  !   scalar particles                   ->  0
  ! Example:
  !   Process 'u g -> W+ d'
  !   "hels(1:4,n)" = (/-1,+1,-1,-1/) means
  !   left-handed up-quark
  !   transverse gluon (helicity = +1)
  !   transverse W+ (helicity = -1)
  !   left-handed down-quark
  ! "hels(1:l,n)" can be used as the input vaiable "hel" in the
  ! subroutines get_amplitude_rcl and
  ! get_polarized_squared_amplitude_rcl.

  integer, intent(in)               :: npr
  integer, intent(out), allocatable :: hels(:,:)

  integer :: pr,i,j,legs

  if (.not.processes_generated) then
    if (warnings(20).le.warning_limit) then
      warnings(20) = warnings(20) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 20: Call of get_helicity_configurations_rcl not allowed:'
      write(nx,*) '          Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(20)
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
    if (warnings(21).le.warning_limit) then
      warnings(21) = warnings(21) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 21: get_helicity_configurations_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(21)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  allocate (hels(1:legs,1:cfTot(pr)))

  do i = 1, cfTot(pr)
    do j = 1,legsIn(pr)
      hels(j,i) = heli(newleg(j,pr),i,pr)
      hels(j,i) = heli(newleg(j,pr),i,pr)
    enddo
    do j = legsIn(pr)+1,legs
      hels(j,i) = - heli(newleg(j,pr),i,pr)
    enddo
  enddo

  end subroutine get_helicity_configurations_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_TIs_required_accuracy_rcl (acc)

  ! This subroutine sets the required accuracy for TIs to the value
  ! "acc" in collier.

  real(dp), intent(in) :: acc

  call SetReqAcc_cll  (acc)

  end subroutine set_TIs_required_accuracy_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_TIs_required_accuracy_rcl (acc)

  ! This subroutine extracts the value of the required accuracy
  ! for TIs from collier and returns it as value of the output
  ! variable "acc".

  real(dp), intent(out) :: acc

  real(dp) :: cllreqacc

  call GetReqAcc_cll (cllreqacc)

  acc = cllreqacc

  end subroutine get_TIs_required_accuracy_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_TIs_critical_accuracy_rcl (acc)

  ! This subroutine sets the critical accuracy for TIs to the value
  ! "acc" in collier.

  real(dp), intent(in) :: acc

  call SetCritAcc_cll  (acc)

  end subroutine set_TIs_critical_accuracy_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_TIs_critical_accuracy_rcl (acc)

  ! This subroutine extracts the value of the critical accuracy
  ! for TIs from collier and returns it as value of the
  ! output variable "acc".

  real(dp), intent(out) :: acc

  real(dp) :: cllcritacc

  call GetCritAcc_cll (cllcritacc)

  acc = cllcritacc

  end subroutine get_TIs_critical_accuracy_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_TIs_accuracy_flag_rcl (flag)

  ! This subroutine extracts the value of the accuracy flag for TIs
  ! from collier. The output variable "flag" returns global
  ! information on the accuracy of the TIs evaluated in the last call
  ! of compute_process_rcl:
  ! - "flag = 0":
  !   For all TIs, the accuracy is estimated to be better than the
  !   required value.
  ! - "flag = -1":
  !   For at least one TI, the accuracy is estimated to be worse than
  !   the required value, but for all TIs, the accuracy is estimated
  !   to be better than the critical value.
  ! - "flag = -2":
  !   For at least one TI, the accuracy is estimated to be worse than
  !   the critical values.
  ! The value of variable "flag" is determined based on internal
  ! uncertainty estimations performed by collier.

  integer, intent(out) :: flag

  integer  :: cllaccflag

  call GetAccFlag_cll (cllaccflag)

  flag = cllaccflag

  end subroutine get_TIs_accuracy_flag_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module process_computation_rcl

!#####################################################################
