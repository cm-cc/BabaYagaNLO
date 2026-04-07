!#####################################################################
!!
!!  File  twoTOfour_all.f90
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

program twoTOfour_all
  use recola, only: set_output_file_rcl
  use check_rcl, only: check_process_rcl
  implicit none
  integer, parameter :: dp = kind (23d0)

  call run_check ('u d~ -> W+ g g g',1d-9)
  call run_check ('u u~ -> Z g g g',1d-9)
  call run_check ('u u~ -> W+ W- g g',1d-8)
  call run_check ('u u~ -> Z Z g g',1d-8)
  call run_check ('d d~ -> t t~ b b~',1d-8)
  call run_check ('g g -> t t~ b b~',1d-8)
  call run_check ('u u~ -> u u~ u u~',1d-8)
  call run_check ('g g -> u u~ u u~',1d-8)
  call run_check ('g g -> u u~ d d~',1d-9)
  call run_check ('g g -> u u~ g g',1d-6)
  call run_check ('g g -> t t~ g g',1d-7)

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

end program twoTOfour_all
