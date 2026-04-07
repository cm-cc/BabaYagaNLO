!#####################################################################
!!
!!  File  twoTOtwo_all.f90
!!  is part of RECOLA (REcursive Computation of One Loop Amplitudes)
!!
!!  Copyright (C) 2015-2017   Stefano Actis, Ansgar Denner,
!!                            Lars Hofer, Jean-Nicolas Lang,
!!                            Andreas Scharf, Sandro Uccirati
!!
!!  RECOLA is licenced under the GNU GPL version 3,
!!         see COPYING for details.
!!
!#####################################################################

program twoTOtwo_all
  use recola, only: set_output_file_rcl
  use check_rcl, only: check_process_rcl
  implicit none
  integer, parameter :: dp = kind (23d0)

  call run_check('u d~ -> e+ nu_e A*',1d-13)
  call run_check('u u~ -> W+ W- g',1d-7)
  call run_check('u u~ -> Z Z g',1d-10)
  call run_check('u u~ -> Z A* g',1d-10)
  call run_check('u u~ -> A* A* g',1d-7)
  call run_check('u d~ -> W+ g g',1d-9)
  call run_check('u d~ -> W+ t t~',1d-10)
  call run_check('u u~ -> Z t t~',1d-9)
  call run_check('d d~ -> Z t t~',1d-9)
  call run_check('g g -> W+ b t~',1d-9)
  call run_check('g g -> Z t t~',1d-9)
  call run_check('u u~ -> Z g g',1d-9)
  call run_check('g g -> g t t~',1d-10)
  call run_check('g g -> g g g',1d-9)
  call run_check('d d~ -> d d~ g',1d-9)
  call run_check('d d~ -> t t~ g',1d-9)
  call run_check('b b~ -> t t~ g',1d-9)

  contains

    subroutine run_check(prstr, threshold)
      character(len=*), intent(in) :: prstr
      real(dp),         intent(in) :: threshold
      real(dp)                     :: delta

      call set_output_file_rcl ('*')
      call check_process_rcl (prstr, delta)
      if (delta .gt. threshold) then
        call EXIT(1)
      end if

    end subroutine run_check

end program twoTOtwo_all
