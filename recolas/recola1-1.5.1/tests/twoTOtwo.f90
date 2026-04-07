!#####################################################################
!!
!!  File  twoTOtwo.f90
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

program twoTOtwo
  use recola, only: set_output_file_rcl
  use check_rcl, only: check_process_rcl
  implicit none
  integer, parameter :: dp = kind (23d0)
  real(dp)           :: delta

  call set_output_file_rcl ('*')

  call check_process_rcl ('u~ u -> nu_e~ nu_e', delta)
  if (delta .gt. 1d-14) then
    call EXIT(1)
  end if

  call check_process_rcl ('u d~ -> nu_e e+', delta)
  if (delta .gt. 1d-14) then
    call EXIT(1)
  end if

  call check_process_rcl ('d~ d -> u~ u', delta)
  if (delta .gt. 1d-14) then
    call EXIT(1)
  end if

  call check_process_rcl ('e+ e- -> nu_e~ nu_e', delta)
  if (delta .gt. 1d-14) then
    call EXIT(1)
  end if

end program twoTOtwo
