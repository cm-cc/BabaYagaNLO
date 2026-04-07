!#####################################################################
!!
!!  File  twoTOthree.f90
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

program twoTOthree
  use recola, only: set_output_file_rcl
  use check_rcl, only: check_process_rcl
  implicit none
  integer, parameter :: dp = kind (23d0)
  real(dp)           :: delta

  call set_output_file_rcl ('*')
  call check_process_rcl ('u u~ -> W+ W- g', delta)
  if (delta .gt. 1d-8) then
    call EXIT(1)
  end if

  call set_output_file_rcl ('*')
  call check_process_rcl ('u u~ -> Z Z g', delta)
  if (delta .gt. 1d-11) then
    call EXIT(1)
  end if

  call set_output_file_rcl ('*')
  call check_process_rcl ('u u~ -> Z A* g', delta)
  if (delta .gt. 1d-11) then
    call EXIT(1)
  end if

  call set_output_file_rcl ('*')
  call check_process_rcl ('g g -> g t t~', delta)
  if (delta .gt. 1d-10) then
    call EXIT(1)
  end if

end program twoTOthree
