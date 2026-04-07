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

  call run_check('u~ u -> nu_e~ nu_e', 1d-13)
  call run_check('u d~ -> nu_e e+', 1d-13)
  call run_check('d~ d -> u~ u', 1d-13)
  call run_check ('e+ e- -> nu_e~ nu_e', 1d-13)
  call run_check ('e+ e- -> W+ W-', 1d-11)
  call run_check ('u d~ -> W+ H', 1d-13)
  call run_check ('e+ e- -> Z H', 1d-14)
  call run_check ('u d~ -> W+ g', 1d-14)

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
