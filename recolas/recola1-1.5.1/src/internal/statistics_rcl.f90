!#####################################################################
!!
!!  File  statistics_rcl.f90
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

  module statistics_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine print_generation_statistics_rcl

  ! This subroutine prints the CPU time spent for the generation of all
  ! defined processes.

  call openOutput
  write(nx,*)
  write(nx,*) '   Generation in',timeGEN,'s'
  write(nx,*)

  end subroutine print_generation_statistics_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine print_TI_statistics_rcl (npr,npoints)

  ! This subroutine prints the CPU time spent in the computation of
  ! the tensor integrals for the process with process number "npr".
  ! The subroutine assumes that the process has been computed
  ! "npoints" times (.i.e. the subroutines compute_process_rcl has
  ! been called "npoints" times with process-number "npr").
  ! The average is made on the number of calls "npoints".

  integer, intent(in) :: npr,npoints
  integer             :: pr,i

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    call openOutput
    write(nx,*)
    write(nx,'(a,i3)') ' ERROR: print_TI_statistics_rcl called '// &
                               'with undefined process index ',npr
  endif

  call openOutput
  write(nx,*)
  write(nx,*) &
    'TIs evaluated in',timeTI(npr)/real(npoints,kind=sp)*1e3, &
    'ms per PS point'
  write(nx,*)

  end subroutine print_TI_statistics_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine print_TC_statistics_rcl (npr,npoints)

  ! This subroutine prints the CPU time spent in the computation of
  ! the tensor coefficients for the process with process number
  ! "npr".
  ! The subroutine assumes that the process has been computed
  ! "npoints" times (.i.e. the subroutines compute_process_rcl has
  ! been called "npoints" times with process-number "npr").
  ! The average is made on the number of calls "npoints".

  integer, intent(in) :: npr,npoints
  integer             :: pr,i

  pr = 0
  do i = 1,prTot
    if (inpr(i).eq.npr) then
      pr = i; exit
    endif
  enddo
  if (pr.eq.0) then
    call openOutput
    write(nx,*)
    write(nx,'(a,i3)') ' ERROR: print_TC_statistics_rcl called '// &
                               'with undefined process index ',npr
  endif

  call openOutput
  write(nx,*)
  write(nx,*) &
    'TCs evaluated in', timeTC(npr)/real(npoints,kind=sp)*1e3, &
    'ms per PS point'
  write(nx,*)

  end subroutine print_TC_statistics_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module statistics_rcl

!#####################################################################



