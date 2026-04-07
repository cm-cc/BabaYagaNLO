!#####################################################################
!!
!!  File  wrapper_rcl.f90
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

  module wrapper_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl
  use process_computation_rcl
  use process_definition_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_alphas_masses_nowidtharg_rcl (mc,mb,mt)

  real(dp), intent(in) :: mc,mb,mt

  call set_alphas_masses_rcl (mc,mb,mt)

  end subroutine set_alphas_masses_nowidtharg_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_gfermi_scheme_noarg_rcl

  call use_gfermi_scheme_rcl

  end subroutine use_gfermi_scheme_noarg_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_gfermi_scheme_and_set_alpha_rcl (alpha)

  real(dp), intent(in) :: alpha

  call use_gfermi_scheme_rcl (a=alpha)

  end subroutine use_gfermi_scheme_and_set_alpha_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_gfermi_scheme_and_set_gfermi_rcl (gfermi)

  real(dp), intent(in) :: gfermi

  call use_gfermi_scheme_rcl (g=gfermi)

  end subroutine use_gfermi_scheme_and_set_gfermi_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_gfermi_scheme_real_and_set_gfermi_rcl (gfermi)

  real(dp), intent(in) :: gfermi

  call use_gfermi_scheme_rcl (g=gfermi,massesgf=0)

  end subroutine use_gfermi_scheme_real_and_set_gfermi_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_gfermi_scheme_complex_and_set_gfermi_rcl (gfermi)

  real(dp), intent(in) :: gfermi

  call use_gfermi_scheme_rcl (g=gfermi,massesgf=1)

  end subroutine use_gfermi_scheme_complex_and_set_gfermi_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_alpha0_scheme_noarg_rcl

  call use_alpha0_scheme_rcl

  end subroutine use_alpha0_scheme_noarg_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_alphaz_scheme_noarg_rcl

  call use_alphaz_scheme_rcl

  end subroutine use_alphaz_scheme_noarg_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_alphaMSbar_scheme_noarg_rcl

  call use_alphaMSbar_scheme_rcl

  end subroutine use_alphaMSbar_scheme_noarg_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_set_alphaGF_rcl (alpha)

  real(dp), intent(in) :: alpha

  call set_gfermi_rcl (a=alpha)

  end subroutine wrapper_set_alphaGF_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_set_gfermi_rcl (gfermi)

  real(dp), intent(in) :: gfermi

  call set_gfermi_rcl (g=gfermi)

  end subroutine wrapper_set_gfermi_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_gfermi_mixed_scheme_noNLOs_rcl (npr,Nal0,NalZ,NalMS)

  integer, intent(in) :: npr,Nal0,NalZ,NalMS

  call use_gfermi_mixed_scheme_rcl (npr,Nal0,NalZ,NalMS)

  end subroutine use_gfermi_mixed_scheme_noNLOs_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_alpha0_mixed_scheme_noNLOs_rcl (npr,Nalgf,NalZ,NalMS)

  integer, intent(in) :: npr,Nalgf,NalZ,NalMS

  call use_alpha0_mixed_scheme_rcl (npr,Nalgf,NalZ,NalMS)

  end subroutine use_alpha0_mixed_scheme_noNLOs_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_alphaz_mixed_scheme_noNLOs_rcl (npr,Nalgf,Nal0,NalMS)

  integer, intent(in) :: npr,Nalgf,Nal0,NalMS

  call use_alphaz_mixed_scheme_rcl (npr,Nalgf,Nal0,NalMS)

  end subroutine use_alphaz_mixed_scheme_noNLOs_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine use_alphaMSbar_mixed_scheme_noNLOs_rcl (npr,Nalgf,Nal0,NalZ)

  integer, intent(in) :: npr,Nalgf,Nal0,NalZ

  call use_alphaMSbar_mixed_scheme_rcl (npr,Nalgf,Nal0,NalZ)

  end subroutine use_alphaMSbar_mixed_scheme_noNLOs_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_set_gs_power_rcl (npr,gsarray,gslen)

  integer, intent(in) :: npr,gslen
  integer, intent(in) :: gsarray(0:1,1:gslen)

  call set_gs_power_rcl(npr,transpose(gsarray))

  end subroutine wrapper_set_gs_power_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_set_outgoing_momenta_rcl(npr,pIn,p,legs)

    use globals_rcl, only: dp
    use outgoing_momenta_rcl, only: set_outgoing_momenta_rcl

    integer,  intent(in)  :: npr,legs
    real(dp), intent(in)  :: pIn(0:3,1:2)
    real(dp), intent(out) :: p(0:3,1:legs)

    call set_outgoing_momenta_rcl(npr,pIn,p)
    p(0:3,1) = pIn(0:3,1)
    p(0:3,2) = pIn(0:3,2)

  end subroutine wrapper_set_outgoing_momenta_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_compute_process_rcl (npr,p,legs,order,A2, &
                                          momenta_check)

!  real(dp), intent(out), dimension(0:1), optional :: A2
! we use ALWAYS the A2 argument as mandatory, since the optional form
! conflicts with passing character arrays ("order") from C to Fortran
! in the same function call

  integer,          intent(in)  :: npr,legs
  real(dp),         intent(in)  :: p(0:3,1:legs)
  character(len=*), intent(in)  :: order
  real(dp),         intent(out) :: A2(0:1)
  logical,          intent(out) :: momenta_check

  call compute_process_rcl (npr,p,order,A2,momenta_check)

  end subroutine wrapper_compute_process_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_rescale_process_rcl (npr,order,A2)

  integer,          intent(in)  :: npr
  character(len=*), intent(in)  :: order
  real(dp),         intent(out) :: A2(0:1)

  call rescale_process_rcl (npr,order,A2)

  end subroutine wrapper_rescale_process_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_get_amplitude_rcl (npr,pow,order,colour,hel,legs,A)

  integer,          intent(in)  :: npr,pow,legs
  integer,          intent(in)  :: colour(1:legs),hel(1:legs)
  character(len=*), intent(in)  :: order
  complex(dp),      intent(out) :: A

  call get_amplitude_rcl (npr,pow,order,colour,hel,A)

  end subroutine wrapper_get_amplitude_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_get_polarized_squared_amplitude_rcl (npr,pow,order,hel,legs,A2h)

  integer,          intent(in)  :: npr,pow,legs
  integer,          intent(in)  :: hel(1:legs)
  character(len=*), intent(in)  :: order
  real(dp),         intent(out) :: A2h

  call get_polarized_squared_amplitude_rcl (npr,pow,order,hel,A2h)

  end subroutine wrapper_get_polarized_squared_amplitude_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_compute_colour_correlation_rcl (npr,p,legs,i1,i2, &
                                                     A2cc,momenta_check)

  integer,  intent(in)  :: npr,legs,i1,i2
  real(dp), intent(in)  :: p(0:3,1:legs)
  real(dp), intent(out) :: A2cc
  logical,  intent(out) :: momenta_check

  call compute_colour_correlation_rcl(npr,p,i1,i2,A2cc,momenta_check)

  end subroutine wrapper_compute_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_rescale_colour_correlation_rcl (npr,i1,i2,A2cc)

  integer,  intent(in)  :: npr,i1,i2
  real(dp), intent(out) :: A2cc

  call rescale_colour_correlation_rcl(npr,i1,i2,A2cc)

  end subroutine wrapper_rescale_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_compute_all_colour_correlations_rcl(npr,p,legs, &
                                                         momenta_check)

  integer,  intent(in)  :: npr,legs
  real(dp), intent(in)  :: p(0:3,1:legs)
  logical,  intent(out) :: momenta_check

  call compute_all_colour_correlations_rcl(npr,p,momenta_check)

  end subroutine wrapper_compute_all_colour_correlations_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_compute_spin_colour_correlation_rcl (npr,p,legs, &
                                                          i1,i2,v, &
                                                          A2scc, &
                                                          momenta_check)

  integer,     intent(in)  :: npr,legs,i1,i2
  real(dp),    intent(in)  :: p(0:3,1:legs)
  complex(dp), intent(in)  :: v(0:3)
  real(dp),    intent(out) :: A2scc
  logical,     intent(out) :: momenta_check

  call compute_spin_colour_correlation_rcl(npr,p,i1,i2,v,A2scc,momenta_check)

  end subroutine wrapper_compute_spin_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_rescale_spin_colour_correlation_rcl (npr,i1,i2,v,A2scc)

  integer,     intent(in)  :: npr,i1,i2
  complex(dp), intent(in)  :: v(0:3)
  real(dp),    intent(out) :: A2scc

  call rescale_spin_colour_correlation_rcl(npr,i1,i2,v,A2scc)

  end subroutine wrapper_rescale_spin_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_compute_spin_correlation_rcl (npr,p,legs,j,v,A2sc, &
                                                   momenta_check)

  integer,     intent(in)  :: npr,legs,j
  real(dp),    intent(in)  :: p(0:3,1:legs)
  complex(dp), intent(in)  :: v(0:3)
  real(dp),    intent(out) :: A2sc
  logical,     intent(out) :: momenta_check

  call compute_spin_correlation_rcl (npr,p,j,v,A2sc,momenta_check)

  end subroutine wrapper_compute_spin_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_rescale_spin_correlation_rcl (npr,j,v,A2sc)

  integer,     intent(in)  :: npr,j
  complex(dp), intent(in)  :: v(0:3)
  real(dp),    intent(out) :: A2sc

  call rescale_spin_correlation_rcl (npr,j,v,A2sc)

  end subroutine wrapper_rescale_spin_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_compute_spin_correlation_matrix_rcl(npr,p,legs,j,&
                                                         A2scm,momenta_check)

    use globals_rcl, only: dp
    use process_computation_rcl, only: compute_spin_correlation_matrix_rcl

    integer,          intent(in)  :: npr,legs,j
    real(dp),         intent(in)  :: p(0:3,1:legs)
    real(dp),         intent(out) :: A2scm(0:3,0:3)
    logical,          intent(out) :: momenta_check

    call compute_spin_correlation_matrix_rcl(npr,p,j,A2scm,momenta_check)

  end subroutine wrapper_compute_spin_correlation_matrix_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_rescale_spin_correlation_matrix_rcl(npr,j,A2scm)

    use globals_rcl, only: dp
    use process_computation_rcl, only: rescale_spin_correlation_matrix_rcl

    integer,          intent(in)  :: npr,j
    real(dp),         intent(out) :: A2scm(0:3,0:3)

    call rescale_spin_correlation_matrix_rcl(npr,j,A2scm)

  end subroutine wrapper_rescale_spin_correlation_matrix_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_compute_int_colour_correlation_rcl (npr,p,legs,i1,i2, &
                                                         A2ccint,momenta_check)

  integer,  intent(in)  :: npr,legs,i1,i2
  real(dp), intent(in)  :: p(0:3,1:legs)
  real(dp), intent(out) :: A2ccint
  logical,  intent(out) :: momenta_check

  call compute_int_colour_correlation_rcl(npr,p,i1,i2,A2ccint,momenta_check)

  end subroutine wrapper_compute_int_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_rescale_int_colour_correlation_rcl (npr,i1,i2,A2ccint)

  integer,  intent(in)  :: npr,i1,i2
  real(dp), intent(out) :: A2ccint

  call rescale_int_colour_correlation_rcl(npr,i1,i2,A2ccint)

  end subroutine wrapper_rescale_int_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_compute_all_int_colour_correlations_rcl(npr,p,legs, &
                                                             momenta_check)

  integer,  intent(in)  :: npr,legs
  real(dp), intent(in)  :: p(0:3,1:legs)
  logical,  intent(out) :: momenta_check

  call compute_all_int_colour_correlations_rcl(npr,p,momenta_check)

  end subroutine wrapper_compute_all_int_colour_correlations_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_compute_nlo_colour_correlation_rcl (npr,p,legs,i1,i2, &
                                                         A2ccnlo,momenta_check)

  integer,  intent(in)  :: npr,legs,i1,i2
  real(dp), intent(in)  :: p(0:3,1:legs)
  real(dp), intent(out) :: A2ccnlo
  logical,  intent(out) :: momenta_check

  call compute_nlo_colour_correlation_rcl(npr,p,i1,i2,A2ccnlo,momenta_check)

  end subroutine wrapper_compute_nlo_colour_correlation_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_legs_rcl (npr,legs)

  ! This subroutine extracts the number of legs of the process with
  ! process number "npr"

  integer, intent(in)  :: npr
  integer, intent(out) :: legs

  integer :: pr,i

  if (.not.processes_generated) then
    if (warnings(391).le.warning_limit) then
      warnings(391) = warnings(391) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 391: Call of get_legs_rcl not allowed:'
      write(nx,*) '           Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(391)
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
    if (warnings(392).le.warning_limit) then
      warnings(392) = warnings(392) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 392: get_legs_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(392)
    endif
    call istop (ifail,1)
  endif

  legs = legsIn(pr) + legsOut(pr)

  end subroutine get_legs_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_cs_rcl(npr,cs_count)

  ! This subroutine extracts the number of colour structures of the
  ! process with process number "npr"

  use globals_rcl, only: csTot

  integer, intent(in)  :: npr
  integer, intent(out) :: cs_count

  integer :: pr,i

  if (.not.processes_generated) then
    if (warnings(393).le.warning_limit) then
      warnings(393) = warnings(393) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 393: Call of get_cs_rcl not allowed:'
      write(nx,*) '           Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(393)
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
    if (warnings(394).le.warning_limit) then
      warnings(394) = warnings(394) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 394: get_cs_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(394)
    endif
    call istop (ifail,1)
  endif

  cs_count = csTot(pr)

  end subroutine get_cs_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_cf_rcl(npr,cf_count)

  ! This subroutine extracts the number of helicity structures of the
  ! process with process number "npr"

  use globals_rcl, only: cfTot

  integer, intent(in)  :: npr
  integer, intent(out) :: cf_count

  integer :: pr,i

  if (.not.processes_generated) then
    if (warnings(395).le.warning_limit) then
      warnings(395) = warnings(395) + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR 395: Call of get_cf_rcl not allowed:'
      write(nx,*) '           Processes not generated yet.'
      write(nx,*)
      call toomanywarnings(395)
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
    if (warnings(396).le.warning_limit) then
      warnings(396) = warnings(396) + 1
      call openOutput
      write(nx,*)
      write(nx,'(a,i3)') ' ERROR 396: get_cf_rcl called '// &
                                 'with undefined process index ',npr
      write(nx,*)
      call toomanywarnings(396)
    endif
    call istop (ifail,1)
  endif

  cf_count = cfTot(pr)

  end subroutine get_cf_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_get_momenta_rcl (npr,p,legs)

  integer,  intent(in)  :: npr,legs
  real(dp), intent(out) :: p(0:3,1:legs)

  call get_momenta_rcl (npr,p)

  end subroutine wrapper_get_momenta_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_get_colour_configurations_rcl (npr,legs,csTot,cols)

  integer, intent(in)  :: npr,legs,csTot
  integer, intent(out) :: cols(1:legs,1:csTot)
  integer, allocatable :: icols(:,:)

  call get_colour_configurations_rcl (npr,icols)
  cols(1:size(icols,1), 1:size(icols,2)) = icols(:,:)

  end subroutine wrapper_get_colour_configurations_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_get_helicity_configurations_rcl (npr,legs,cfTot,hels)

  integer, intent(in)  :: npr,legs,cfTot
  integer, intent(out) :: hels(1:legs,1:cfTot)
  integer, allocatable :: ihels(:,:)

  call get_helicity_configurations_rcl (npr,ihels)
  hels(:, :) = ihels(:, :)

  end subroutine wrapper_get_helicity_configurations_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wrapper_get_recola_version_rcl(ret_version,slen)
    use globals_rcl, only: get_recola_version_rcl
    character(len=10), intent(out) :: ret_version
    integer,           intent(out) :: slen

    call get_recola_version_rcl(ret_version)
    slen = len(trim(ret_version))
    ret_version = trim(ret_version)//CHAR(0)

  end subroutine wrapper_get_recola_version_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module wrapper_rcl

!#####################################################################
