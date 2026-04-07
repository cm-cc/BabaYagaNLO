!!#####################################################################
!!
!!  File  model_vertices_rcl.f90
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

  module model_vertices_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine vertices_tables
  ! fields incoming

  integer     :: i,f1,f2,f3,f4

  ! vertices_tables
  allocate(ve2ct(0:2,11:nFs,11:nFs),&
           ve2r2(0:2,11:nFs,11:nFs),&
           ve3tr(0:1,nFs,nFs,nFs),  &
           ve3ct(0:3,11:nFs,11:nFs,11:nFs),&
           ve3r2(0:3,11:nFs,11:nFs,11:nFs),&
           ve4tr(0:2,11:19,11:19,11:19,11:19),&
           ve4ct(0:4,11:19,11:19,11:19,11:19),&
           ve4r2(0:4,11:19,11:19,11:19,11:19))

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CT 2-leg vertices
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ve2ct = .false.

  ! 2 scalars
  ve2ct(0,11,11) = .true. ! H H
  ve2ct(0,12,12) = .true. ! p0 p0
  ve2ct(0,13,14) = .true. ! p+ p-

  ! 1 scalar, 1 vector
  ve2ct(0,12,17) = .true. ! p0 Z
  ve2ct(0,13,19) = .true. ! p+ W-
  ve2ct(0,14,18) = .true. ! p- W+

  ! 2 vectors
  ve2ct(2,15,15) = .true. ! g g
  ve2ct(0,16,16) = .true. ! A A
  ve2ct(0,16,17) = .true. ! A Z
  ve2ct(0,17,17) = .true. ! Z Z
  ve2ct(0,18,19) = .true. ! W+ W-

  ! 2 fermions
  do i = 0,2
    f1 = 20 + i; f2 = 32 + i
    ve2ct (0,f1,f2) = .true. ! nu nu~
    f1 = 23 + i; f2 = 35 + i
    ve2ct (0,f1,f2) = .true. ! u u~
    ve2ct (2,f1,f2) = .true. ! u u~
    f1 = 26 + i; f2 = 38 + i
    ve2ct (0,f1,f2) = .true. ! l- l+
    f1 = 29 + i; f2 = 41 + i
    ve2ct (0,f1,f2) = .true. ! d d~
    ve2ct (2,f1,f2) = .true. ! d d~
  enddo

  do f1 = 11,nFs
    do f2 = f1,nFs
      do i = 0,size(ve2ct,1)-1
        if (ve2ct(i,f1,f2)) ve2ct(i,f2,f1) = .true.
      enddo
    enddo
  enddo

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! R2 2-leg vertices
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ve2r2 = ve2ct

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Tree 3-leg vertices
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ve3tr = .false.

  ! 2 FP ghost fields and 1 scalar
  ve3tr (0, 2, 9,13) = .true. ! X_A X+~ p+
  ve3tr (0, 2,10,14) = .true. ! X_A X-~ p-
  ve3tr (0, 3, 8,11) = .true. ! X_Z X_Z~ H
  ve3tr (0, 3, 9,13) = .true. ! X_Z X+~ p+
  ve3tr (0, 3,10,14) = .true. ! X_Z X-~ p-
  ve3tr (0, 4, 8,14) = .true. ! X+ X_Z~ p-
  ve3tr (0, 4, 9,11) = .true. ! X+ X+~ H
  ve3tr (0, 4, 9,12) = .true. ! X+ X+~ p0
  ve3tr (0, 5, 8,13) = .true. ! X- X_Z~ p+
  ve3tr (0, 5,10,11) = .true. ! X- X-~ H
  ve3tr (0, 5,10,12) = .true. ! X- X-~ p0

  ! 2 FP ghost fields and 1 vector
  ve3tr (1, 1, 6,15) = .true. ! X_g X_g~ g
  ve3tr (0, 2, 9,18) = .true. ! X_A X+~ W+
  ve3tr (0, 2,10,19) = .true. ! X_A X-~ W-
  ve3tr (0, 3, 9,18) = .true. ! X_Z X+~ W+
  ve3tr (0, 3,10,19) = .true. ! X_Z X-~ W-
  ve3tr (0, 4, 7,19) = .true. ! X+ X_A~ W-
  ve3tr (0, 4, 8,19) = .true. ! X+ X_Z~ W-
  ve3tr (0, 4, 9,16) = .true. ! X+ X+~ A
  ve3tr (0, 4, 9,17) = .true. ! X+ X+~ Z
  ve3tr (0, 5, 7,18) = .true. ! X- X_A~ W+
  ve3tr (0, 5, 8,18) = .true. ! X- X_Z~ W+
  ve3tr (0, 5,10,16) = .true. ! X- X-~ A
  ve3tr (0, 5,10,17) = .true. ! X- X-~ Z

  ! 3 scalars
  ve3tr (0,11,11,11) = .true. ! H H H
  ve3tr (0,11,12,12) = .true. ! H p0 p0
  ve3tr (0,11,13,14) = .true. ! H p+ p-

  ! 2 scalars and 1 vector
  ve3tr (0,11,12,17) = .true. ! H p0 Z
  ve3tr (0,11,13,19) = .true. ! H p+ W-
  ve3tr (0,11,14,18) = .true. ! H p- W+
  ve3tr (0,12,13,19) = .true. ! p0 p+ W-
  ve3tr (0,12,14,18) = .true. ! p0 p- W+
  ve3tr (0,13,14,16) = .true. ! H p0 A
  ve3tr (0,13,14,17) = .true. ! H p0 Z

  ! 1 scalar and 2 vectors
  ve3tr (0,11,17,17) = .true. ! H Z Z
  ve3tr (0,11,18,19) = .true. ! H W+ W-
  ve3tr (0,13,16,19) = .true. ! p+ A W-
  ve3tr (0,13,17,19) = .true. ! p+ Z W-
  ve3tr (0,14,16,18) = .true. ! p- A W+
  ve3tr (0,14,17,18) = .true. ! p- Z W+

  ! 3 vectors
  ve3tr (1,15,15,15) = .true. ! g g g
  ve3tr (0,16,18,19) = .true. ! A W+ W-
  ve3tr (0,17,18,19) = .true. ! Z W+ W-

  ! 1 scalar and 2 fermions
  do i = 0,2
    f2 = 23 + i; f3 = 35 + i
    ve3tr (0,11,f2,f3) = .true. ! H u u~
    ve3tr (0,12,f2,f3) = .true. ! p0 u u~
    f2 = 26 + i; f3 = 38 + i
    ve3tr (0,11,f2,f3) = .true. ! H l- l+
    ve3tr (0,12,f2,f3) = .true. ! p0 l- l+
    f2 = 29 + i; f3 = 41 + i
    ve3tr (0,11,f2,f3) = .true. ! H d d~
    ve3tr (0,12,f2,f3) = .true. ! p0 d d~
    f2 = 26 + i; f3 = 32 + i
    ve3tr (0,13,f2,f3) = .true. ! p+ l- nu~
    f2 = 29 + i; f3 = 35 + i
    ve3tr (0,13,f2,f3) = .true. ! p+ d u~
    f2 = 20 + i; f3 = 38 + i
    ve3tr (0,14,f2,f3) = .true. ! p- nu l+
    f2 = 23 + i; f3 = 41 + i
    ve3tr (0,14,f2,f3) = .true. ! p- u d~
  enddo

  ! 1 vector and 2 fermions
  do i = 0,2
    f2 = 20 + i; f3 = 32 + i
    ve3tr (0,17,f2,f3) = .true. ! Z nu nu~
    f2 = 23 + i; f3 = 35 + i
    ve3tr (1,15,f2,f3) = .true. ! g u u~
    ve3tr (0,16,f2,f3) = .true. ! A u u~
    ve3tr (0,17,f2,f3) = .true. ! Z u u~
    f2 = 26 + i; f3 = 38 + i
    ve3tr (0,16,f2,f3) = .true. ! A l- l+
    ve3tr (0,17,f2,f3) = .true. ! Z l- l+
    f2 = 29 + i; f3 = 41 + i
    ve3tr (1,15,f2,f3) = .true. ! g d d~
    ve3tr (0,16,f2,f3) = .true. ! A d d~
    ve3tr (0,17,f2,f3) = .true. ! Z d d~
    f2 = 26 + i; f3 = 32 + i
    ve3tr (0,18,f2,f3) = .true. ! W+ l- nu~
    f2 = 29 + i; f3 = 35 + i
    ve3tr (0,18,f2,f3) = .true. ! W+ d u~
    f2 = 20 + i; f3 = 38 + i
    ve3tr (0,19,f2,f3) = .true. ! W- nu l+
    f2 = 23 + i; f3 = 41 + i
    ve3tr (0,19,f2,f3) = .true. ! W- u d~
  enddo

  do f1 = 1,nFs
    do f2 = f1,nFs
      do f3 = f2,nFs
        do i = 0,size(ve3tr,1)-1
          if (ve3tr(i,f1,f2,f3)) then
            ve3tr(i,f1,f3,f2) = .true.
            ve3tr(i,f2,f1,f3) = .true.
            ve3tr(i,f2,f3,f1) = .true.
            ve3tr(i,f3,f1,f2) = .true.
            ve3tr(i,f3,f2,f1) = .true.
          endif
        enddo
      enddo
    enddo
  enddo

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CT 3-leg vertices
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ve3ct = .false.

  ve3ct(0,11:,11:,11:) = ve3tr(0,11:,11:,11:)

  ve3ct(3,11:,11:,11:) = ve3tr(1,11:,11:,11:)

  ! 2 scalars and 1 vector
  ve3ct (0,11,12,16) = .true. ! H p0 A

  ! 1 scalar and 2 vectors
  ve3ct (0,11,16,17) = .true. ! H A Z

  ! 1 scalar and 2 fermions
  do i = 0,2
    f2 = 23 + i; f3 = 35 + i
    ve3ct (2,11,f2,f3) = .true. ! H u u~
    ve3ct (2,12,f2,f3) = .true. ! p0 u u~
    f2 = 29 + i; f3 = 41 + i
    ve3ct (2,11,f2,f3) = .true. ! H d d~
    ve3ct (2,12,f2,f3) = .true. ! p0 d d~
    f2 = 29 + i; f3 = 35 + i
    ve3ct (2,13,f2,f3) = .true. ! p+ d u~
    f2 = 23 + i; f3 = 41 + i
    ve3ct (2,14,f2,f3) = .true. ! p- u d~
  enddo

  ! 1 vector and 2 fermions
  do i = 0,2
    f2 = 20 + i; f3 = 32 + i
    ve3ct (0,16,f2,f3) = .true. ! A nu nu~
    f2 = 23 + i; f3 = 35 + i
    ve3ct (1,15,f2,f3) = .true. ! g u u~
    ve3ct (2,16,f2,f3) = .true. ! A u u~
    ve3ct (2,17,f2,f3) = .true. ! Z u u~
    f2 = 29 + i; f3 = 41 + i
    ve3ct (1,15,f2,f3) = .true. ! g d d~
    ve3ct (2,16,f2,f3) = .true. ! A d d~
    ve3ct (2,17,f2,f3) = .true. ! Z d d~
    f2 = 29 + i; f3 = 35 + i
    ve3ct (2,18,f2,f3) = .true. ! W+ d u~
    f2 = 23 + i; f3 = 41 + i
    ve3ct (2,19,f2,f3) = .true. ! W- u d~
  enddo

  do f1 = 11,nFs
    do f2 = f1,nFs
      do f3 = f2,nFs
        do i = 0,size(ve3ct,1)-1
          if (ve3ct(i,f1,f2,f3)) then
            ve3ct(i,f1,f3,f2) = .true.
            ve3ct(i,f2,f1,f3) = .true.
            ve3ct(i,f2,f3,f1) = .true.
            ve3ct(i,f3,f1,f2) = .true.
            ve3ct(i,f3,f2,f1) = .true.
          endif
        enddo
      enddo
    enddo
  enddo

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! R2 3-leg vertices
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ve3r2 = ve3ct

  ! 1 scalar and 2 vectors
  ve3r2 (2,11,15,15) = .true. ! H g g
  ve3r2 (2,15,11,15) = .true. ! g H g
  ve3r2 (2,15,15,11) = .true. ! g g H
  ve3r2 (0,11,16,16) = .true. ! H A A
  ve3r2 (0,16,11,16) = .true. ! A H A
  ve3r2 (0,16,16,11) = .true. ! A A H

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Tree 4-leg vertices
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ve4tr = .false.

  ! 4 scalars
  ve4tr (0,11,11,11,11) = .true. ! H H H H
  ve4tr (0,11,11,12,12) = .true. ! H H p0 p0
  ve4tr (0,11,11,13,14) = .true. ! H H p+ p-
  ve4tr (0,12,12,12,12) = .true. ! p0 p0 p0 p0
  ve4tr (0,12,12,13,14) = .true. ! p0 p0 p+ p-
  ve4tr (0,13,13,14,14) = .true. ! p+ p+ p- p-

  ! 2 scalars and 2 vectors
  ve4tr (0,11,11,17,17) = .true. ! H H Z Z
  ve4tr (0,11,11,18,19) = .true. ! H H W+ W-
  ve4tr (0,11,13,16,19) = .true. ! H p+ A W-
  ve4tr (0,11,13,17,19) = .true. ! H p+ Z W-
  ve4tr (0,11,14,16,18) = .true. ! H p- A W+
  ve4tr (0,11,14,17,18) = .true. ! H p- Z W+
  ve4tr (0,12,12,17,17) = .true. ! p0 p0 Z Z
  ve4tr (0,12,12,18,19) = .true. ! p0 p0 W+ W-
  ve4tr (0,12,13,16,19) = .true. ! p0 p+ A W-
  ve4tr (0,12,13,17,19) = .true. ! p0 p+ Z W-
  ve4tr (0,12,14,16,18) = .true. ! p0 p- A W+
  ve4tr (0,12,14,17,18) = .true. ! p0 p- Z W+
  ve4tr (0,13,14,16,16) = .true. ! p+ p- A A
  ve4tr (0,13,14,16,17) = .true. ! p+ p- A Z
  ve4tr (0,13,14,17,17) = .true. ! p+ p- Z Z
  ve4tr (0,13,14,18,19) = .true. ! p+ p- W+ W-

  ! 4 vectors
  ve4tr (2,15,15,15,15) = .true. ! g g g g
  ve4tr (0,16,16,18,19) = .true. ! A A W+ W-
  ve4tr (0,16,17,18,19) = .true. ! A Z W+ W-
  ve4tr (0,17,17,18,19) = .true. ! Z Z W+ W-
  ve4tr (0,18,18,19,19) = .true. ! W+ W+ W- W-

  do f1 = 11,19
    do f2 = f1,19
      do f3 = f2,19
        do f4 = f3,19
          do i = 0,size(ve4tr,1)-1
            if (ve4tr(i,f1,f2,f3,f4)) then
              ve4tr(i,f1,f2,f4,f3) = .true.
              ve4tr(i,f1,f3,f2,f4) = .true.
              ve4tr(i,f1,f3,f4,f2) = .true.
              ve4tr(i,f1,f4,f2,f3) = .true.
              ve4tr(i,f1,f4,f3,f2) = .true.
              ve4tr(i,f2,f1,f3,f4) = .true.
              ve4tr(i,f2,f1,f4,f3) = .true.
              ve4tr(i,f2,f3,f1,f4) = .true.
              ve4tr(i,f2,f3,f4,f1) = .true.
              ve4tr(i,f2,f4,f1,f3) = .true.
              ve4tr(i,f2,f4,f3,f1) = .true.
              ve4tr(i,f3,f1,f2,f4) = .true.
              ve4tr(i,f3,f1,f4,f2) = .true.
              ve4tr(i,f3,f2,f1,f4) = .true.
              ve4tr(i,f3,f2,f4,f1) = .true.
              ve4tr(i,f3,f4,f1,f2) = .true.
              ve4tr(i,f3,f4,f2,f1) = .true.
              ve4tr(i,f4,f1,f2,f3) = .true.
              ve4tr(i,f4,f1,f3,f2) = .true.
              ve4tr(i,f4,f2,f1,f3) = .true.
              ve4tr(i,f4,f2,f3,f1) = .true.
              ve4tr(i,f4,f3,f1,f2) = .true.
              ve4tr(i,f4,f3,f2,f1) = .true.
            endif
          enddo
        enddo
      enddo
    enddo
  enddo

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CT 4-leg vertices
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ve4ct = .false.

  ve4ct(0,:,:,:,:) = ve4tr(0,:,:,:,:)

  ve4ct(4,:,:,:,:) = ve4tr(2,:,:,:,:)

  ! 2 scalars and 2 vectors
  ve4ct (0,11,11,16,17) = .true. ! H H A Z
  ve4ct (0,12,12,16,17) = .true. ! p0 p0 A Z

  do f1 = 11,19
    do f2 = f1,19
      do f3 = f2,19
        do f4 = f3,19
          do i = 0,size(ve4ct,1)-1
            if (ve4ct(i,f1,f2,f3,f4)) then
              ve4ct(i,f1,f2,f4,f3) = .true.
              ve4ct(i,f1,f3,f2,f4) = .true.
              ve4ct(i,f1,f3,f4,f2) = .true.
              ve4ct(i,f1,f4,f2,f3) = .true.
              ve4ct(i,f1,f4,f3,f2) = .true.
              ve4ct(i,f2,f1,f3,f4) = .true.
              ve4ct(i,f2,f1,f4,f3) = .true.
              ve4ct(i,f2,f3,f1,f4) = .true.
              ve4ct(i,f2,f3,f4,f1) = .true.
              ve4ct(i,f2,f4,f1,f3) = .true.
              ve4ct(i,f2,f4,f3,f1) = .true.
              ve4ct(i,f3,f1,f2,f4) = .true.
              ve4ct(i,f3,f1,f4,f2) = .true.
              ve4ct(i,f3,f2,f1,f4) = .true.
              ve4ct(i,f3,f2,f4,f1) = .true.
              ve4ct(i,f3,f4,f1,f2) = .true.
              ve4ct(i,f3,f4,f2,f1) = .true.
              ve4ct(i,f4,f1,f2,f3) = .true.
              ve4ct(i,f4,f1,f3,f2) = .true.
              ve4ct(i,f4,f2,f1,f3) = .true.
              ve4ct(i,f4,f2,f3,f1) = .true.
              ve4ct(i,f4,f3,f1,f2) = .true.
              ve4ct(i,f4,f3,f2,f1) = .true.
            endif
          enddo
        enddo
      enddo
    enddo
  enddo

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! R2 4-leg vertices
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ve4r2 = ve4ct

  ! 2 scalars and 2 vectors
  ve4r2 (2,11,11,15,15) = .true. ! H H g g
  ve4r2 (2,12,12,15,15) = .true. ! p0 p0 g g
  ve4r2 (2,13,14,15,15) = .true. ! p+ p- g g
  ve4r2 (0,11,11,16,16) = .true. ! H H A A
  ve4r2 (0,12,12,16,16) = .true. ! p0 p0 A A

  ! 4 vectors
  ve4r2 (3,15,15,15,16) = .true. ! g g g A
  ve4r2 (3,15,15,15,17) = .true. ! g g g Z
  ve4r2 (2,15,15,16,16) = .true. ! g g A A
  ve4r2 (2,15,15,16,17) = .true. ! g g A Z
  ve4r2 (2,15,15,17,17) = .true. ! g g Z Z
  ve4r2 (2,15,15,18,19) = .true. ! g g W+ W-
  ve4r2 (0,16,16,16,16) = .true. ! A A A A
  ve4r2 (0,16,16,16,17) = .true. ! A A A Z
  ve4r2 (0,16,16,17,17) = .true. ! A A Z Z
  ve4r2 (0,16,17,17,17) = .true. ! A Z Z Z
  ve4r2 (0,17,17,17,17) = .true. ! Z Z Z Z

  do f1 = 11,19
    do f2 = f1,19
      do f3 = f2,19
        do f4 = f3,19
          do i = 0,size(ve4r2,1)-1
            if (ve4r2(i,f1,f2,f3,f4)) then
              ve4r2(i,f1,f2,f4,f3) = .true.
              ve4r2(i,f1,f3,f2,f4) = .true.
              ve4r2(i,f1,f3,f4,f2) = .true.
              ve4r2(i,f1,f4,f2,f3) = .true.
              ve4r2(i,f1,f4,f3,f2) = .true.
              ve4r2(i,f2,f1,f3,f4) = .true.
              ve4r2(i,f2,f1,f4,f3) = .true.
              ve4r2(i,f2,f3,f1,f4) = .true.
              ve4r2(i,f2,f3,f4,f1) = .true.
              ve4r2(i,f2,f4,f1,f3) = .true.
              ve4r2(i,f2,f4,f3,f1) = .true.
              ve4r2(i,f3,f1,f2,f4) = .true.
              ve4r2(i,f3,f1,f4,f2) = .true.
              ve4r2(i,f3,f2,f1,f4) = .true.
              ve4r2(i,f3,f2,f4,f1) = .true.
              ve4r2(i,f3,f4,f1,f2) = .true.
              ve4r2(i,f3,f4,f2,f1) = .true.
              ve4r2(i,f4,f1,f2,f3) = .true.
              ve4r2(i,f4,f1,f3,f2) = .true.
              ve4r2(i,f4,f2,f1,f3) = .true.
              ve4r2(i,f4,f2,f3,f1) = .true.
              ve4r2(i,f4,f3,f1,f2) = .true.
              ve4r2(i,f4,f3,f2,f1) = .true.
            endif
          enddo
        enddo
      enddo
    enddo
  enddo

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Checks
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  f = 0
!
!  do f1 = 11,nFs; f(1) = f1
!    do f2 = 11,nFs; f(2) = f2
!      do i = 0,size(ve2ct,1)-1
!        co0 = cou(2,i,f)
!        ve = sum(abs(co0)) .gt. zerocut
!        if (ve2ct(i,f1,f2).neqv.ve) then
!          call openOutput
!          write(nx,*) trim(cpar(f1))//' '//trim(cpar(f2)),' - CT',i,ve2ct(i,f1,f2),ve
!        endif
!      enddo
!      do i = 0,size(ve2r2,1)-1
!        co0 = cou(3,i,f)
!        ve = sum(abs(co0)) .gt. zerocut
!        if (ve2r2(i,f1,f2).neqv.ve) then
!          call openOutput
!          write(nx,*) trim(cpar(f1))//' '//trim(cpar(f2)),' - R2',i,ve2r2(i,f1,f2),ve
!        endif
!      enddo
!    enddo
!  enddo
!
!  do f1 = 1,nFs; f(1) = f1
!    do f2 = 1,nFs; f(2) = f2
!      do f3 = 1,nFs; f(3) = f3
!        do i = 0,size(ve3tr,1)-1
!          co0 = cou(0,i,f)
!          ve = sum(abs(co0)) .gt. zerocut
!          if (ve3tr(i,f1,f2,f3).neqv.ve) then
!            call openOutput
!            write(nx,*) trim(cpar(f1))//' '//trim(cpar(f2))//' '//trim(cpar(f3)),i
!          endif
!        enddo
!      enddo
!    enddo
!  enddo
!
!  do f1 = 11,nFs; f(1) = f1
!    do f2 = 11,nFs; f(2) = f2
!      do f3 = 11,nFs; f(3) = f3
!        do i = 0,size(ve3ct,1)-1
!          co0 = cou(2,i,f)
!          ve = sum(abs(co0)) .gt. zerocut
!          if (ve3ct(i,f1,f2,f3).neqv.ve) then
!            call openOutput
!            write(nx,*) trim(cpar(f1))//' '//trim(cpar(f2))//' '//trim(cpar(f3)),' - CT',i
!          endif
!        enddo
!        do i = 0,size(ve3r2,1)-1
!          co0 = cou(3,i,f)
!          ve = sum(abs(co0)) .gt. zerocut
!          if (ve3r2(i,f1,f2,f3).neqv.ve) then
!            call openOutput
!            write(nx,*) trim(cpar(f1))//' '//trim(cpar(f2))//' '//trim(cpar(f3)),' - R2',i
!          endif
!        enddo
!      enddo
!    enddo
!  enddo
!
!  do f1 = 11,19; f(1) = f1
!    do f2 = 11,19; f(2) = f2
!      do f3 = 11,19; f(3) = f3
!        do f4 = 11,19; f(4) = f4
!          do i = 0,size(ve4tr,1)-1
!            co0 = cou(0,i,f)
!            ve = sum(abs(co0)) .gt. zerocut
!            if (ve4tr(i,f1,f2,f3,f4).neqv.ve) then
!              call openOutput
!              write(nx,*) trim(cpar(f1))//' '//trim(cpar(f2))//' '// &
!                          trim(cpar(f3))//' '//trim(cpar(f4)),i
!            endif
!          enddo
!          do i = 0,size(ve4ct,1)-1
!            co0 = cou(2,i,f)
!            ve = sum(abs(co0)) .gt. zerocut
!            if (ve4ct(i,f1,f2,f3,f4).neqv.ve) then
!              call openOutput
!              write(nx,*) trim(cpar(f1))//' '//trim(cpar(f2))//' '// &
!                          trim(cpar(f3))//' '//trim(cpar(f4)),' - CT',i
!            endif
!          enddo
!          do i = 0,size(ve4r2,1)-1
!            co0 = cou(3,i,f)
!            ve = sum(abs(co0)) .gt. zerocut
!            if (ve4r2(i,f1,f2,f3,f4).neqv.ve) then
!              call openOutput
!              write(nx,*) trim(cpar(f1))//' '//trim(cpar(f2))//' '// &
!                          trim(cpar(f3))//' '//trim(cpar(f4)),' - R2',i
!            endif
!          enddo
!        enddo
!      enddo
!    enddo
!  enddo

  end subroutine vertices_tables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine vert2 ( l,lp,x,f1,f2,cs1,        & ! in
                     v2e,tyv,cs4,coef,colcoef ) ! out
  ! f1,f2 incoming

  integer,                  intent (in)  :: l,lp,x,f1,f2,cs1(-1:)
  integer,     allocatable, intent (out) :: cs4(:,:)
  integer,                  intent (out) :: tyv
  logical,                  intent (out) :: v2e(0:2)
  complex(dp), allocatable, intent (out) :: coef(:,:,:)
  real(dp),    allocatable, intent (out) :: colcoef(:)

  integer       :: i,xlp,f(4)
  logical       :: v2eSum,masslessf2
  complex(dp)   :: co0(1:4),co(1:4,0:2)
  character(99) :: cty


  xlp = x*lp

  select case (xlp)
  case (0); v2e = .false.
  case (2); v2e(0:2) = ve2ct(0:2,f1,f2)
  case (3); v2e(0:2) = ve2r2(0:2,f1,f2)
  end select

  v2eSum = .false.
  do i = 0,size(v2e)-1
    v2eSum = v2eSum .or. v2e(i)
  enddo

  if (.not.v2eSum) return

  f = 0; f(1) = f1; f(2) = f2
  v2eSum = .false.
  do i = 0,size(v2e)-1
    co0 = cou(xlp,i,f)
    co(1:4,i) = co0(1:4)
    select case (dynamic_settings)
    case (0); v2e(i) = sum(abs(co(:,i))) .gt. zerocut
    case (1)
      if (v2e(i).and.(xlp.eq.2)) then; co(:,i) = c1d0
      else; v2e(i) = sum(abs(co(:,i))) .gt. zerocut
      endif
    end select
    v2eSum = v2eSum .or. v2e(i)
  enddo

  if (.not.v2eSum) return

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  tyv = 0

  cty = trim(cftype(f1))//' -> '//trim(cftype(anti(f2)))

  masslessf2 = regf(f2) .le. 2

  select case (cty)
  case   ('s -> s'); tyv = 1001
  case   ('v -> s'); tyv = 1002
  case   ('s -> v'); tyv = 1011
  case   ('v -> v'); tyv = 1012
  case   ('f -> f')
    select case (masslessf2)
    case (.false.);  tyv = 1021 ! massive  on f2
    case (.true.);   tyv = 1022 ! massless on f2
    end select
  case ('f~ -> f~')
    select case (masslessf2)
    case (.false.);  tyv = 1041 ! massive  on f2
    case (.true.);   tyv = 1042 ! massless on f2
    end select
  case default
    if (warnings(310).le.warning_limit) then
      warnings(310) = warnings(310) + 1
      call openOutput
      write(nx,*)
      write(nx,*) &
        'CODE ERROR 310 (model_vertices_rcl.f90): coupling ',cty,' not defined'
      write(nx,*)
      call toomanywarnings(310)
    endif
    call istop (ifail,2)
  end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! colour
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  allocate (cs4(-1:l,1)); cs4(-1:l,1) = cs1(-1:l)

  allocate (colcoef(1)); colcoef(1) = 1d0

  allocate (coef(4,1,0:2)); coef(1:4,1,0:2) = co(1:4,0:2)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (pureQED) then
    do i = 1,2
      select case (f(i))
        case (1:15,17:19); v2e = .false.
      end select
    enddo
  endif

  end subroutine vert2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine vert3 (l,lp,x,last,csL,f1,f2,f3,cs1,cs2,U1g, & ! in
                    v3e,tyv,cs4,U1g4,coef,colcoef         ) ! out
  ! f1,f2,f3 incoming

  integer,                  intent (in)  :: l,lp,x,csL(-1:0),f1,f2, &
                                            f3,cs1(-1:),cs2(-1:)
  logical,                  intent (in)  :: last,U1g(2)
  integer, allocatable,     intent (out) :: cs4(:,:)
  integer,                  intent (out) :: tyv
  logical,                  intent (out) :: v3e(0:3),U1g4
  complex(dp), allocatable, intent (out) :: coef(:,:,:)
  real(dp), allocatable,    intent (out) :: colcoef(:)

  integer       :: i,j,xlp,f(4),n=1,iq1,iq2,iq3,ia1,ia2,ia3,cs(-1:l,6)
  logical       :: v3eSum,masslessf3
  complex(dp)   :: co0(1:4),co(1:2,0:3)
  real(dp)      :: cl(6)
  character(99) :: cty,cty2


  xlp = x*lp

  select case (xlp)
  case (0); v3e(0:1) = ve3tr(0:1,f1,f2,f3); v3e(2:3) = .false.
  case (2); v3e(0:3) = ve3ct(0:3,f1,f2,f3)
  case (3); v3e(0:3) = ve3r2(0:3,f1,f2,f3)
  end select

  v3eSum = .false.
  do i = 0,size(v3e)-1
    v3eSum = v3eSum .or. v3e(i)
  enddo

  if (.not.v3eSum) return

  f = 0; f(1) = f1; f(2) = f2; f(3) = f3
  v3eSum = .false.
  do i = 0,size(v3e)-1
    co0 = cou(xlp,i,f)
    co(1:2,i) = co0(1:2)
    select case (dynamic_settings)
    case (0); v3e(i) = sum(abs(co(:,i))) .gt. zerocut
    case (1)
      if (v3e(i).and.(xlp.eq.2)) then; co(:,i) = c1d0
      else; v3e(i) = sum(abs(co(:,i))) .gt. zerocut
      endif
    end select
    v3eSum = v3eSum .or. v3e(i)
  enddo

  if (.not.v3eSum) return

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  tyv = 0

  cty = trim(cftype(f1))//' '//trim(cftype(f2))//' -> '// &
         trim(cftype(anti(f3)))

  masslessf3 = regf(f3) .le. 2

  select case (cty)
  case   ('s s -> s');              tyv =  1
  case   ('v v -> s');              tyv =  2
  case   ('s v -> s');              tyv =  3
  case   ('v s -> s');              tyv =  4
  case  ('f f~ -> s');              tyv =  5
  case  ('f~ f -> s');              tyv =  6
  case   ('s s -> v');              tyv = 11
  case   ('v v -> v');              tyv = 12
  case   ('s v -> v');              tyv = 13
  case   ('v s -> v');              tyv = 14
  case  ('f f~ -> v');              tyv = 15
  case  ('f~ f -> v');              tyv = 16
  case   ('f s -> f')
    select case (masslessf3)
    case (.false.);                 tyv = 21 ! massive  on f3
    case (.true.);                  tyv = 22 ! massless on f3
    end select
  case   ('s f -> f')
    select case (masslessf3)
    case (.false.);                 tyv = 23 ! massive  on f3
    case (.true.);                  tyv = 24 ! massless on f3
    end select
  case   ('f v -> f')
    select case (masslessf3)
    case (.false.);                 tyv = 31 ! massive  on f3
    case (.true.);                  tyv = 32 ! massless on f3
    end select
  case   ('v f -> f')
    select case (masslessf3)
    case (.false.);                 tyv = 33 ! massive  on f3
    case (.true.);                  tyv = 34 ! massless on f3
    end select
  case ('f~ s -> f~')
    select case (masslessf3)
    case (.false.);                 tyv = 41 ! massive  on f3
    case (.true.);                  tyv = 42 ! massless on f3
    end select
  case ('s f~ -> f~')
    select case (masslessf3)
    case (.false.);                 tyv = 43 ! massive  on f3
    case (.true.);                  tyv = 44 ! massless on f3
    end select
  case ('f~ v -> f~')
    select case (masslessf3)
    case (.false.);                 tyv = 51 ! massive  on f3
    case (.true.);                  tyv = 52 ! massless on f3
    end select
  case ('v f~ -> f~')
    select case (masslessf3)
    case (.false.);                 tyv = 53 ! massive  on f3
    case (.true.);                  tyv = 54 ! massless on f3
    end select
  case  ('x s -> x',  's x -> x');  tyv = 61
  case  ('x v -> x');               tyv = 62
  case  ('v x -> x');               tyv = 63
  case ('x~ s -> x~','s x~ -> x~'); tyv = 71
  case ('x~ v -> x~');              tyv = 72
  case ('v x~ -> x~');              tyv = 73
  case default
    if (warnings(311).le.warning_limit) then
      warnings(311) = warnings(311) + 1
      call openOutput
      write(nx,*)
      write(nx,*) &
        'CODE ERROR 311 (model_vertices_rcl.f90): coupling ',cty,' not defined'
      write(nx,*)
      call toomanywarnings(311)
    endif
    call istop (ifail,2)
  end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! colour
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! cs is the OUTGOING colour structure
! cs(-1) = open iq index
! cs( 0) = open ia index
! cs( i) = ia index in the delta with iq=i (i=1,...,legs)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! iq = incoming quark index
! ia = incoming antiquark index
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! "ia" is the first index in matrices and delta's, "iq" is the second
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cty2 = trim(cftype2(f1))//' '//trim(cftype2(f2))//' '//trim(cftype2(f3))

  n = 1

  cs(:0,1) = 0; cs(1:,1) = cs1(1:) + cs2(1:)

  U1g4 = .false.

  iq1 = cs1(-1); ia1 = cs1(0)
  iq2 = cs2(-1); ia2 = cs2(0)
  iq3 = csL(-1); ia3 = csL(0)

  select case (cty2)
  case ( 'G g G~','g G~ G'  )
    n = 2
    cs(:0,2) = 0; cs(1:,2) = cs1(1:) + cs2(1:)
    !  1) iq: 1 2 3
    !     ia: 3 1 2
    if (last) then; cs(iq1,1) = ia3; else; cs(-1,1) = iq1; endif
                    cs(iq2,1) = ia1
    if (last) then; cs(iq3,1) = ia2; else; cs( 0,1) = ia2; endif
    cl(1) = + 1d0
    !  2) iq: 1 2 3
    !     ia: 2 3 1
                    cs(iq1,2) = ia2
    if (last) then; cs(iq2,2) = ia3; else; cs(-1,2) = iq2; endif
    if (last) then; cs(iq3,2) = ia1; else; cs( 0,2) = ia1; endif
    cl(2) = - 1d0
  case ( 'G~ g G','g G G~','g g g' )
    n = 2
    cs(:0,2) = 0; cs(1:,2) = cs1(1:) + cs2(1:)
    !  1) iq: 1 2 3
    !     ia: 2 3 1
                    cs(iq1,1) = ia2
    if (last) then; cs(iq2,1) = ia3; else; cs(-1,1) = iq2; endif
    if (last) then; cs(iq3,1) = ia1; else; cs( 0,1) = ia1; endif
    cl(1) = + 1d0
    !  2) iq: 1 2 3
    !     ia: 3 1 2
    if (last) then; cs(iq1,2) = ia3; else; cs(-1,2) = iq1; endif
                    cs(iq2,2) = ia1
    if (last) then; cs(iq3,2) = ia2; else; cs( 0,2) = ia2; endif
    cl(2) = - 1d0
  case ( 'g q q~' )
    !  1) iq: 1 2
    !     ia: 3 1
    if (last) then; cs(iq1,1) = ia3; else; cs(-1,1) = iq1; endif
                    cs(iq2,1) = ia1
    cl(1) = + 1d0
    if (U1g(1)) then
      !  2) iq: 1 2
      !     ia: 1 3
      n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:)
                      cs(iq1,n) = ia1
      if (last) then; cs(iq2,n) = ia3; else; cs(-1,n) = iq2; endif
      cl(n) = - 1d0/Nc
    endif
  case ( 'g q~ q' )
    !  1) iq: 1 3
    !     ia: 2 1
                    cs(iq1,1) = ia2
    if (last) then; cs(iq3,1) = ia1; else; cs( 0,1) = ia1; endif
    cl(1) = + 1d0
    if (U1g(1)) then
      !  2) iq: 1 3
      !     ia: 1 2
      n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:)
                      cs(iq1,n) = ia1
      if (last) then; cs(iq3,n) = ia2; else; cs( 0,n) = ia2; endif
      cl(n) = - 1d0/Nc
    endif
  case ( 'q g q~' )
    !  1) iq: 1 2
    !     ia: 2 3
                    cs(iq1,1) = ia2
    if (last) then; cs(iq2,1) = ia3; else; cs(-1,1) = iq2; endif
    cl(1) = + 1d0
    if (U1g(2)) then
      !  2) iq: 1 2
      !     ia: 3 2
      n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:)
      if (last) then; cs(iq1,n) = ia3; else; cs(-1,n) = iq1; endif
                      cs(iq2,n) = ia2
      cl(n) = - 1d0/Nc
    endif
  case ( 'q~ g q' )
    !  1) iq: 2 3
    !     ia: 1 2
                    cs(iq2,1) = ia1
    if (last) then; cs(iq3,1) = ia2; else; cs( 0,1) = ia2; endif
    cl(1) = + 1d0
    if (U1g(2)) then
      !  2) iq: 2 3
      !     ia: 2 1
      n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:)
                      cs(iq2,n) = ia2
      if (last) then; cs(iq3,n) = ia1; else; cs( 0,n) = ia1; endif
      cl(n) = - 1d0/Nc
    endif
  case ( 'q q~ g' )
    !  1) iq: 1 3
    !     ia: 3 2
    if (last) then; cs(iq1,1) = ia3; else; cs(-1,1) = iq1; endif
    if (last) then; cs(iq3,1) = ia2; else; cs( 0,1) = ia2; endif
    cl(1) = + 1d0
    if (last) then;
      !  2) iq: 1 3
      !     ia: 2 3
      n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:)
      cs(iq1,n) = ia2
      cs(iq3,n) = ia3
      cl(n) = - 1d0/Nc
    endif
    U1g4 = .true.
  case ( 'q~ q g' )
    !  1) iq: 2 3
    !     ia: 3 1
    if (last) then; cs(iq2,1) = ia3; else; cs(-1,1) = iq2; endif
    if (last) then; cs(iq3,1) = ia1; else; cs( 0,1) = ia1; endif
    cl(1) = + 1d0
    if (last) then;
      !  2) iq: 2 3
      !     ia: 1 3
      n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:)
      cs(iq2,n) = ia1
      cs(iq3,n) = ia3
      cl(n) = - 1d0/Nc
    endif
    U1g4 = .true.
  case ( 'v q q~','s q q~' )
    !  1) iq: 2
    !     ia: 3
    if (last) then; cs(iq2,1) = ia3; else; cs(-1,1) = iq2; endif
    cl(1) = 1d0
  case ( 'v q~ q','s q~ q' )
    !  1) iq: 3
    !     ia: 2
    if (last) then; cs(iq3,1) = ia2; else; cs( 0,1) = ia2; endif
    cl(1) = 1d0
  case ( 'q v q~','q s q~' )
    !  1) iq: 1
    !     ia: 3
    if (last) then; cs(iq1,1) = ia3; else; cs(-1,1) = iq1; endif
    cl(1) = 1d0
  case ( 'q~ v q','q~ s q' )
    !  1) iq: 3
    !     ia: 1
    if (last) then; cs(iq3,1) = ia1; else; cs( 0,1) = ia1; endif
    cl(1) = 1d0
  case ( 'q q~ v','q q~ s' )
    !  1) iq: 1
    !     ia: 2
    cs(iq1,1) = ia2
    cl(1) = 1d0
  case ( 'q~ q v','q~ q s' )
    !  1) iq: 2
    !     ia: 1
    cs(iq2,1) = ia1
    cl(1) = 1d0
  case ( 'g v g','g s g' )
    !  1) iq: 1 3
    !     ia: 3 1
    if (last) then; cs(iq1,1) = ia3; else; cs(-1,1) = iq1; endif
    if (last) then; cs(iq3,1) = ia1; else; cs( 0,1) = ia1; endif
    cl(1) = + 1d0
    if (last.and.U1g(1)) then;
      !  2) iq: 1 3
      !     ia: 1 3
      n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:)
      cs(iq1,n) = ia1
      cs(iq3,n) = ia3
      cl(n) = - 1d0/Nc
    endif
    U1g4 = U1g(1)
  case ( 'v g g','s g g' )
    !  1) iq: 2 3
    !     ia: 3 2
    if (last) then; cs(iq2,1) = ia3; else; cs(-1,1) = iq2; endif
    if (last) then; cs(iq3,1) = ia2; else; cs( 0,1) = ia2; endif
    cl(1) = + 1d0
    if (last.and.U1g(2)) then;
      !  2) iq: 2 3
      !     ia: 2 3
      n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:)
      cs(iq2,n) = ia2
      cs(iq3,n) = ia3
      cl(n) = - 1d0/Nc
    endif
    U1g4 = U1g(2)
  case ( 'g g v','g g s' )
    !  1) iq: 1 2
    !     ia: 2 1
    cs(iq1,1) = ia2
    cs(iq2,1) = ia1
    cl(1) = + 1d0
    if (U1g(1).and.U1g(2)) then;
      !  2) iq: 1 2
      !     ia: 1 2
      n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:)
      cs(iq1,n) = ia1
      cs(iq2,n) = ia2
      cl(n) = - 1d0/Nc
    endif
  case default
    cl(1) = 1d0
  end select


  if ((last).and.(lp.eq.1)) then
    do i = 1,n
      if (cs(l-1,i).eq.l) then
        cl(i) = Nc*cl(i)
      else
        do j = 1,l
          if (cs(j,i).eq.l) cs(j,i) = cs(l-1,i)
        enddo
      endif
      cs(l-1,i) = 0
      if (cs(l,i).eq.l-1) then
        cl(i) = Nc*cl(i)
      else
        do j = 1,l
          if (cs(j,i).eq.l-1) cs(j,i) = cs(l,i)
        enddo
      endif
      cs(l,i) = 0
    enddo
  endif

  allocate (cs4(-1:l,n)); cs4(-1:l,1:n) = cs(-1:l,1:n)

  allocate (colcoef(n)); colcoef(1:n) = cl(1:n)

  allocate (coef(2,n,0:3))
  do i = 1,n
    coef(1:2,i,0:3) = co(1:2,0:3)
  enddo

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (pureQED) then
    do i = 1,3
      select case (f(i))
        case (1:15,17:19); v3e = .false.
      end select
    enddo
  endif

  end subroutine vert3

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine vert4 ( l,lp,x,last,csL,f1,f2,f3,f4,cs1,cs2,cs3,U1g, & ! in
                     v4e,tyv,cs4,U1g4,coef,colcoef                ) ! out
  ! f1,f2,f3,f4 incoming

  integer,                  intent (in)  :: l,lp,x,csL(-1:0),f1,f2,  &
                                            f3,f4,cs1(-1:),cs2(-1:), &
                                            cs3(-1:)
  logical,                  intent (in)  :: last,U1g(3)
  integer, allocatable,     intent (out) :: cs4(:,:)
  integer,                  intent (out) :: tyv
  logical,                  intent (out) :: v4e(0:4),U1g4
  complex(dp), allocatable, intent (out) :: coef(:,:,:)
  real(dp), allocatable,    intent (out) :: colcoef(:)

  integer       :: i,j,xlp,f(4),n=1,iq1,iq2,iq3,iq4,ia1,ia2,ia3,ia4, &
                   cs(-1:l,24),l1,l2
  logical       :: v4eSum
  complex(dp)   :: co0(1:4),co(1:3,0:4),cA,cB,c(3,24,0:4)
  real(dp)      :: cl(24)
  character(99) :: cty,cty2


  xlp = x*lp

  select case (xlp)
  case (0); v4e(0:2) = ve4tr(0:2,f1,f2,f3,f4); v4e(3:4) = .false.
  case (2); v4e(0:4) = ve4ct(0:4,f1,f2,f3,f4)
  case (3); v4e(0:4) = ve4r2(0:4,f1,f2,f3,f4)
  end select
  v4eSum = .false.
  do i = 0,size(v4e)-1
    v4eSum = v4eSum .or. v4e(i)
  enddo

  if (.not.v4eSum) return

  f = 0; f(1) = f1; f(2) = f2; f(3) = f3; f(4) = f4
  v4eSum = .false.
  do i = 0,size(v4e)-1
    co0 = cou(xlp,i,f)
    co(1:3,i) = co0(1:3)
    select case (dynamic_settings)
    case (0); v4e(i) = sum(abs(co(:,i))) .gt. zerocut
    case (1)
      if (v4e(i).and.(xlp.eq.2)) then; co(:,i) = c1d0
      else; v4e(i) = sum(abs(co(:,i))) .gt. zerocut
      endif
    end select
    v4eSum = v4eSum .or. v4e(i)
  enddo

  if (.not.v4eSum) return

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  tyv = 0

  cty = trim(cftype(f1))//' '//trim(cftype(f2))//' '// &
        trim(cftype(f3))//' -> '//trim(cftype(anti(f4)))

  select case (cty)
  case ( 's s s -> s' ); tyv = -  1
  case ( 's v v -> s' ); tyv = -  2
  case ( 'v s v -> s' ); tyv = -  3
  case ( 'v v s -> s' ); tyv = -  4
  case ( 'v v v -> v' ); tyv = - 11
  case ( 's s v -> v' ); tyv = - 12
  case ( 'v s s -> v' ); tyv = - 13
  case ( 's v s -> v' ); tyv = - 14
  case default
    if (warnings(312).le.warning_limit) then
      warnings(312) = warnings(312) + 1
      call openOutput
      write(nx,*)
      write(nx,*) &
        'CODE ERROR 312 (model_vertices_rcl.f90): coupling ',cty,' not defined'
      write(nx,*)
      call toomanywarnings(312)
    endif
    call istop (ifail,2)
  end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! colour
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! cs is the OUTGOING colour structure
! cs(-1) = open iq index
! cs( 0) = open ia index
! cs( i) = ia index in the delta with iq=i (i=1,...,legs)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! iq = incoming quark index
! ia = incoming antiquark index
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! "ia" is the first index in matrices and delta's, "iq" is the second
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cty2 = trim(cftype2(f1))//' '//trim(cftype2(f2))//' '// &
         trim(cftype2(f3))//' '//trim(cftype2(f4))

  n = 1

  cs(:0,1) = 0; cs(1:,1) = cs1(1:) + cs2(1:) + cs3(1:)

  c(:,1,0:4) = c1d0

  U1g4 = .false.

  select case (cty2)

  case ( 'g g g g' )

    iq1 = cs1(-1); ia1 = cs1(0)
    iq2 = cs2(-1); ia2 = cs2(0)
    iq3 = cs3(-1); ia3 = cs3(0)
    iq4 = csL(-1); ia4 = csL(0)

    n = 6
    do i = 2,n
      cs(:0,i) = 0; cs(1:,i) = cs1(1:) + cs2(1:) + cs3(1:)
    enddo

    !  1) iq: 1 2 3 4
    !     ia: 4 1 2 3
    if (last) then; cs(iq1,1) = ia4; else; cs(-1,1) = iq1; endif
                    cs(iq2,1) = ia1
                    cs(iq3,1) = ia2
    if (last) then; cs(iq4,1) = ia3; else; cs( 0,1) = ia3; endif

    !  2) iq: 1 2 3 4
    !     ia: 2 3 4 1
                    cs(iq1,2) = ia2
                    cs(iq2,2) = ia3
    if (last) then; cs(iq3,2) = ia4; else; cs(-1,2) = iq3; endif
    if (last) then; cs(iq4,2) = ia1; else; cs( 0,2) = ia1; endif

    !  3) iq: 1 2 3 4
    !     ia: 3 4 2 1
                    cs(iq1,3) = ia3
    if (last) then; cs(iq2,3) = ia4; else; cs(-1,3) = iq2; endif
                    cs(iq3,3) = ia2
    if (last) then; cs(iq4,3) = ia1; else; cs( 0,3) = ia1; endif

    !  4) iq: 1 2 3 4
    !     ia: 4 3 1 2
    if (last) then; cs(iq1,4) = ia4; else; cs(-1,4) = iq1; endif
                    cs(iq2,4) = ia3
                    cs(iq3,4) = ia1
    if (last) then; cs(iq4,4) = ia2; else; cs( 0,4) = ia2; endif

    !  5) iq: 1 2 3 4
    !     ia: 3 1 4 2
                    cs(iq1,5) = ia3
                    cs(iq2,5) = ia1
    if (last) then; cs(iq3,5) = ia4; else; cs(-1, 5) = iq3; endif
    if (last) then; cs(iq4,5) = ia2; else; cs( 0, 5) = ia2; endif

    !  6) iq: 1 2 3 4
    !     ia: 2 4 1 3
                    cs(iq1,6) = ia2
    if (last) then; cs(iq2,6) = ia4; else; cs(-1, 6) = iq2; endif
                    cs(iq3,6) = ia1
    if (last) then; cs(iq4,6) = ia3; else; cs( 0, 6) = ia3; endif

    cl(1:6) = 1d0

    if (xlp.eq.3) then

      !  iq: 1 2 3 4
      !  ia: 4 3 2 1
      n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
      if (last) then; cs(iq1,n) = ia4; else; cs(-1,n) = iq1; endif
                      cs(iq2,n) = ia3
                      cs(iq3,n) = ia2
      if (last) then; cs(iq4,n) = ia1; else; cs( 0,n) = ia1; endif
      cl(n) = 1d0

      !  iq: 1 2 3 4
      !  ia: 3 4 1 2
      n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
                      cs(iq1,n) = ia3
      if (last) then; cs(iq2,n) = ia4; else; cs(-1,n) = iq2; endif
                      cs(iq3,n) = ia1
      if (last) then; cs(iq4,n) = ia2; else; cs( 0,n) = ia2; endif
      cl(n) = 1d0

      !  iq: 1 2 3 4
      !  ia: 2 1 4 3
      n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
                      cs(iq1,n) = ia2
                      cs(iq2,n) = ia1
      if (last) then; cs(iq3,n) = ia4; else; cs(-1,n) = iq3; endif
      if (last) then; cs(iq4,n) = ia3; else; cs( 0,n) = ia3; endif
      cl(n) = 1d0

      if (U1g(1)) then

        !  iq: 1 2 3 4
        !  ia: 1 4 2 3
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
                        cs(iq1,n) = ia1
        if (last) then; cs(iq2,n) = ia4; else; cs(-1,n) = iq2; endif
                        cs(iq3,n) = ia2
        if (last) then; cs(iq4,n) = ia3; else; cs( 0,n) = ia3; endif
        cl(n) = - 1d0/2d0 * ( 1d0 - Nf/Nc )

        !  iq: 1 2 3 4
        !  ia: 1 3 4 2
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
                        cs(iq1,n) = ia1
                        cs(iq2,n) = ia3
        if (last) then; cs(iq3,n) = ia4; else; cs(-1,n) = iq3; endif
        if (last) then; cs(iq4,n) = ia2; else; cs( 0,n) = ia2; endif
        cl(n) = - 1d0/2d0 * ( 1d0 - Nf/Nc )

      endif

      if (U1g(2)) then

        !  iq: 1 2 3 4
        !  ia: 3 2 4 1
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
                        cs(iq1,n) = ia3
                        cs(iq2,n) = ia2
        if (last) then; cs(iq3,n) = ia4; else; cs(-1,n) = iq3; endif
        if (last) then; cs(iq4,n) = ia1; else; cs( 0,n) = ia1; endif
        cl(n) = - 1d0/2d0 * ( 1d0 - Nf/Nc )

        !  iq: 1 2 3 4
        !  ia: 4 2 1 3
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
        if (last) then; cs(iq1,n) = ia4; else; cs(-1,n) = iq1; endif
                        cs(iq2,n) = ia2
                        cs(iq3,n) = ia1
        if (last) then; cs(iq4,n) = ia3; else; cs( 0,n) = ia3; endif
        cl(n) = - 1d0/2d0 * ( 1d0 - Nf/Nc )

      endif

      if (U1g(3)) then

        !  iq: 1 2 3 4
        !  ia: 4 1 3 2
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
        if (last) then; cs(iq1,n) = ia4; else; cs(-1,n) = iq1; endif
                        cs(iq2,n) = ia1
                        cs(iq3,n) = ia3
        if (last) then; cs(iq4,n) = ia2; else; cs( 0,n) = ia2; endif
        cl(n) = - 1d0/2d0 * ( 1d0 - Nf/Nc )

        !  iq: 1 2 3 4
        !  ia: 2 4 3 1
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
                        cs(iq1,n) = ia2
        if (last) then; cs(iq2,n) = ia4; else; cs(-1,n) = iq2; endif
                        cs(iq3,n) = ia3
        if (last) then; cs(iq4,n) = ia1; else; cs( 0,n) = ia1; endif
        cl(n) = - 1d0/2d0 * ( 1d0 - Nf/Nc )

      endif

      if (U1g(1).and.U1g(2)) then
        !  iq: 1 2 3 4
        !  ia: 1 2 4 3
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
                        cs(iq1,n) = ia1
                        cs(iq2,n) = ia2
        if (last) then; cs(iq3,n) = ia4; else; cs(-1,n) = iq3; endif
        if (last) then; cs(iq4,n) = ia3; else; cs( 0,n) = ia3; endif
        cl(n) = - Nf/Nc**2
      endif

      if (U1g(1).and.U1g(3)) then
        !  iq: 1 2 3 4
        !  ia: 1 4 3 2
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
                        cs(iq1,n) = ia1
        if (last) then; cs(iq2,n) = ia4; else; cs(-1,n) = iq2; endif
                        cs(iq3,n) = ia3
        if (last) then; cs(iq4,n) = ia2; else; cs( 0,n) = ia2; endif
        cl(n) = - Nf/Nc**2
      endif

      if (U1g(2).and.U1g(3)) then
        !  iq: 1 2 3 4
        !  ia: 4 2 3 1
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
        if (last) then; cs(iq1,n) = ia4; else; cs(-1,n) = iq1; endif
                        cs(iq2,n) = ia2
                        cs(iq3,n) = ia3
        if (last) then; cs(iq4,n) = ia1; else; cs( 0,n) = ia1; endif
        cl(n) = - Nf/Nc**2
      endif

      if (last) then

        !  iq: 1 2 3 4
        !  ia: 3 1 2 4
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
        cs(iq1,n) = ia3
        cs(iq2,n) = ia1
        cs(iq3,n) = ia2
        cs(iq4,n) = ia4
        cl(n) = - 1d0/2d0 * ( 1d0 - Nf/Nc )

        !  iq: 1 2 3 4
        !  ia: 2 3 1 4
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
        cs(iq1,n) = ia2
        cs(iq2,n) = ia3
        cs(iq3,n) = ia1
        cs(iq4,n) = ia4
        cl(n) = - 1d0/2d0 * ( 1d0 - Nf/Nc )

        if (U1g(1)) then
          !  iq: 1 2 3 4
          !  ia: 1 3 2 4
          n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
          cs(iq1,n) = ia1
          cs(iq2,n) = ia3
          cs(iq3,n) = ia2
          cs(iq4,n) = ia4
          cl(n) = - Nf/Nc**2
        endif

        if (U1g(2)) then
          !  iq: 1 2 3 4
          !  ia: 3 2 1 4
          n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
          cs(iq1,n) = ia3
          cs(iq2,n) = ia2
          cs(iq3,n) = ia1
          cs(iq4,n) = ia4
          cl(n) = - Nf/Nc**2
        endif

        if (U1g(3)) then
          !  iq: 1 2 3 4
          !  ia: 2 1 3 4
          n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
          cs(iq1,n) = ia2
          cs(iq2,n) = ia1
          cs(iq3,n) = ia3
          cs(iq4,n) = ia4
          cl(n) = - Nf/Nc**2
        endif

        if (U1g(1).and.U1g(2).and.U1g(3)) then
          !  iq: 1 2 3 4
          !  ia: 1 2 3 4
          n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
          cs(iq1,n) = ia1
          cs(iq2,n) = ia2
          cs(iq3,n) = ia3
          cs(iq4,n) = ia4
          cl(n)    = + 3*Nf/Nc**3
        endif

      endif

      U1g4 = .true.

    endif

    select case (xlp)
    case (0,2)
      cA = + 2d0
      cB = - 1d0
    case (3)
      cA = + Nc*( 3d0 + lam + 5*Nf/(2*Nc) )
      cB = - Nc/2d0*( 5/2d0 + lam + 3*Nf/Nc )
    end select

    c(1,1:2,0:4) = cB
    c(2,1:2,0:4) = cA
    c(3,1:2,0:4) = cB

    c(1,3:4,0:4) = cB
    c(2,3:4,0:4) = cB
    c(3,3:4,0:4) = cA

    c(1,5:6,0:4) = cA
    c(2,5:6,0:4) = cB
    c(3,5:6,0:4) = cB

    c(:,7:n,0:4) = c1d0

  case ( 'g g g v' )

    if (xlp.eq.3) then

      iq1 = cs1(-1); ia1 = cs1(0)
      iq2 = cs2(-1); ia2 = cs2(0)
      iq3 = cs3(-1); ia3 = cs3(0)

      n = 2
      do i = 2,n
        cs(:0,i) = 0; cs(1:,i) = cs1(1:) + cs2(1:) + cs3(1:)
      enddo

      !  1) iq: 1 2 3
      !     ia: 2 3 1
      cs(iq1,1) = ia2
      cs(iq2,1) = ia3
      cs(iq3,1) = ia1

      !  2) iq: 1 2 3
      !     ia: 3 1 2
      cs(iq1,2) = ia3
      cs(iq2,2) = ia1
      cs(iq3,2) = ia2

      cl(1:2) = 1d0

      if (U1g(1)) then
        !  iq: 1 2 3
        !  ia: 1 3 2
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
        cs(iq1,n) = ia1
        cs(iq2,n) = ia3
        cs(iq3,n) = ia2
        cl(n) = - 2d0/Nc
      endif

      if (U1g(2)) then
        !  iq: 1 2 3
        !  ia: 3 2 1
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
        cs(iq1,n) = ia3
        cs(iq2,n) = ia2
        cs(iq3,n) = ia1
        cl(n) = - 2d0/Nc
      endif

      if (U1g(3)) then
        !  iq: 1 2 3
        !  ia: 2 1 3
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
        cs(iq1,n) = ia2
        cs(iq2,n) = ia1
        cs(iq3,n) = ia3
        cl(n) = - 2d0/Nc
      endif

      if (U1g(1).and.U1g(2).and.U1g(3)) then
        !  iq: 1 2 3
        !  ia: 1 2 3
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
        cs(iq1,n) = ia1
        cs(iq2,n) = ia2
        cs(iq3,n) = ia3
        cl(n) = 4d0/Nc**2
      endif

    endif

    c(:,1:n,0:4) = c1d0

  case ( 'g g v g','g v g g','v g g g' )

    if (xlp.eq.3) then

      if (f3.ne.15) then
        l1 = 1; l2 = 2
        iq1 = cs1(-1); ia1 = cs1(0)
        iq2 = cs2(-1); ia2 = cs2(0)
        iq3 = csL(-1); ia3 = csL(0)
      elseif (f2.ne.15) then
        l1 = 1; l2 = 3
        iq1 = cs1(-1); ia1 = cs1(0)
        iq2 = cs3(-1); ia2 = cs3(0)
        iq3 = csL(-1); ia3 = csL(0)
      elseif (f1.ne.15) then
        l1 = 2; l2 = 3
        iq1 = cs2(-1); ia1 = cs2(0)
        iq2 = cs3(-1); ia2 = cs3(0)
        iq3 = csL(-1); ia3 = csL(0)
      endif

      n = 2
      do i = 2,n
        cs(:0,i) = 0; cs(1:,i) = cs1(1:) + cs2(1:) + cs3(1:)
      enddo

      !  1) iq: 1 2 3
      !     ia: 2 3 1
                      cs(iq1,1) = ia2
      if (last) then; cs(iq2,1) = ia3; else; cs(-1,1) = iq2; endif
      if (last) then; cs(iq3,1) = ia1; else; cs( 0,1) = ia1; endif

      !  2) iq: 1 2 3
      !     ia: 3 1 2
      if (last) then; cs(iq1,2) = ia3; else; cs(-1,2) = iq1; endif
                      cs(iq2,2) = ia1
      if (last) then; cs(iq3,2) = ia2; else; cs( 0,2) = ia2; endif

      cl(1:2) = 1d0

      if (U1g(l1)) then
        !  iq: 1 2 3
        !  ia: 1 3 2
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
                        cs(iq1,n) = ia1
        if (last) then; cs(iq2,n) = ia3; else; cs(-1,n) = iq2; endif
        if (last) then; cs(iq3,n) = ia2; else; cs( 0,n) = ia2; endif
        cl(n) = - 2d0/Nc
      endif

      if (U1g(l2)) then
        !  iq: 1 2 3
        !  ia: 3 2 1
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
        if (last) then; cs(iq1,n) = ia3; else; cs(-1,n) = iq1; endif
                        cs(iq2,n) = ia2
        if (last) then; cs(iq3,n) = ia1; else; cs( 0,n) = ia1; endif
        cl(n) = - 2d0/Nc
      endif

      if (last) then

        !  iq: 1 2 3
        !  ia: 2 1 3
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
        cs(iq1,n) = ia2
        cs(iq2,n) = ia1
        cs(iq3,n) = ia3
        cl(n) = - 2d0/Nc

        if (U1g(l1).and.U1g(l2)) then
          !  iq: 1 2 3
          !  ia: 1 2 3
          n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
          cs(iq1,n) = ia1
          cs(iq2,n) = ia2
          cs(iq3,n) = ia3
          cl(n) = 4d0/Nc**2
        endif

      endif

    endif

    c(:,1:n,0:4) = c1d0

    U1g4 = .true.

  case ( 'g g v v','g v g v','v g g v','g g s s','g s g s','s g g s' )

    if (xlp.eq.3) then

      if (f3.ne.15) then
        l1 = 1; l2 = 2
        iq1 = cs1(-1); ia1 = cs1(0)
        iq2 = cs2(-1); ia2 = cs2(0)
      elseif (f2.ne.15) then
        l1 = 1; l2 = 3
        iq1 = cs1(-1); ia1 = cs1(0)
        iq2 = cs3(-1); ia2 = cs3(0)
      elseif (f1.ne.15) then
        l1 = 2; l2 = 3
        iq1 = cs2(-1); ia1 = cs2(0)
        iq2 = cs3(-1); ia2 = cs3(0)
      endif

      !  1) iq: 1 2
      !     ia: 2 1
      cs(iq1,1) = ia2
      cs(iq2,1) = ia1
      cl(1) = 1d0

      if (U1g(l1).and.U1g(l2)) then
        !  iq: 1 2
        !  ia: 1 2
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
        cs(iq1,n) = ia1
        cs(iq2,n) = ia2
        cl(n) = - 1d0/Nc
      endif

    endif

    c(:,1:n,0:4) = c1d0

  case ( 'g v v g','v g v g','v v g g','g s s g','s g s g','s s g g' )

    if (xlp.eq.3) then

      if (f1.eq.15) then
        l1 = 1
        iq1 = cs1(-1); ia1 = cs1(0)
        iq2 = csL(-1); ia2 = csL(0)
      elseif (f2.eq.15) then
        l1 = 2
        iq1 = cs2(-1); ia1 = cs2(0)
        iq2 = csL(-1); ia2 = csL(0)
      elseif (f3.eq.15) then
        l1 = 3
        iq1 = cs3(-1); ia1 = cs3(0)
        iq2 = csL(-1); ia2 = csL(0)
      endif

      !  1) iq: 1 2
      !     ia: 2 1
      if (last) then; cs(iq1,1) = ia2; else; cs(-1,1) = iq1; endif;
      if (last) then; cs(iq2,1) = ia1; else; cs( 0,1) = ia1; endif;
      cl(1) = 1d0

      if (last.and.U1g(l1)) then
        !  2) iq: 1 2
        !     ia: 1 2
        n = n + 1; cs(:0,n) = 0; cs(1:,n) = cs1(1:) + cs2(1:) + cs3(1:)
        cs(iq1,n) = ia1
        cs(iq2,n) = ia2
        cl(n) = - 1d0/Nc
      endif

    endif

    c(:,1:n,0:4) = c1d0

    U1g4 = U1g(l1)

  case default

    cl(1) = 1d0

  end select

  if ((last).and.(lp.eq.1)) then
    do i = 1,n
      if (cs(l-1,i).eq.l) then
        cl(i) = Nc*cl(i)
      else
        do j = 1,l
          if (cs(j,i).eq.l) cs(j,i) = cs(l-1,i)
        enddo
      endif
      cs(l-1,i) = 0
      if (cs(l,i).eq.l-1) then
        cl(i) = Nc*cl(i)
      else
        do j = 1,l
          if (cs(j,i).eq.l-1) cs(j,i) = cs(l,i)
        enddo
      endif
      cs(l,i) = 0
    enddo
  endif

  allocate (cs4(-1:l,n)); cs4(-1:l,1:n) = cs(-1:l,1:n)

  allocate (colcoef(n)); colcoef(1:n) = cl(1:n)

  allocate (coef(3,n,0:4))
  do i = 1,n
    coef(1:3,i,0:4) = c(1:3,i,0:4) * co(1:3,0:4)
  enddo

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (pureQED) then
    do i = 1,4
      select case (f(i))
        case (1:15,17:19); v4e = .false.
      end select
    enddo
  endif

  end subroutine vert4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function cou (xlp,ngs,f) result (co)
  ! f(:) all incoming

  integer, intent (in) :: xlp,ngs,f(:)
  complex(dp)          :: co(1:4)

  integer       :: i,n,ff(size(f,1))
  complex(dp)   :: cP(1:4),cQED(1:4),c0,cA,cB
  real(dp)      :: el,el2,el3,el4,gs,gs2,gs3,gs4
  character(99) :: cv


  ! EW couplings
  el  = 2*sqrt(pi*alpha)
  el2 = el*el
  el3 = el2*el
  el4 = el2*el2

  ! QCD coupling
  gs  = 2*sqrt(pi*als0)
  gs2 = gs*gs
  gs3 = gs2*gs
  gs4 = gs2*gs2

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  n  = 0
  cv = ''
  do i = 1,size(f,1)
    if (f(i).ne.0) then
      n  = n + 1
      ff(n) = f(i)
      cv = trim(cv)//' '//trim(cpar(f(i)))
    endif
  enddo
  cv = trim(adjustl(cv))

  co = c0d0
  cP = c0d0

  cQED = c0d0

  select case (cv)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 2 scalars
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! coupling = i*( co(1)*ps - co(2) )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  case ( 'H H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = dZh
        co(2) = mh2*dZh + dmh2
      case (3)
        co(1) = - el2/(96*pi2*st2)*(                   &
                  + 1/2d0*( 1d0 + 1/(2*ct2) )          & ! bosonic
                  + 1/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  + Nc/mw2*(                           & ! hadronic
                    + mu2(1) + mu2(2) + mu2(3)         & ! hadronic
                    + md2(1) + md2(2) + md2(3)         & ! hadronic
                    )                                  & ! hadronic
                  )
        co(2) = - el2/(16*pi2*st2)*(                             &
                  + mw2/4d0*( 1d0 + 1/(2*ct4) )*( 1d0 - 12*lam ) & ! bosonic
                  + mw2/4d0 + mz2/(8*ct2)                        & ! bosonic
                  + 1/mw2*( ml4(1) + ml4(2) + ml4(3) )           & ! leptonic
                  + Nc/mw2*(                                     & ! hadronic
                    + mu4(1) + mu4(2) + mu4(3)                   & ! hadronic
                    + md4(1) + md4(2) + md4(3)                   & ! hadronic
                    )                                            & ! hadronic
                  )
      end select
    end select

  case ( 'p0 p0' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = dZp0
!        co(2) = - el*dt/(2*st*mw1) + dmz2 ! old (according to 0709.1075)
        co(2) = - el*dt/(2*st*mw1)
      case (3)
        co(1) = - el2/(96*pi2*st2)*(                   &
                  + 1/2d0*( 1d0 + 1/(2*ct2) )          & ! bosonic
                  + 1/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  + Nc/mw2*(                           & ! hadronic
                    + mu2(1) + mu2(2) + mu2(3)         & ! hadronic
                    + md2(1) + md2(2) + md2(3)         & ! hadronic
                    )                                  & ! hadronic
                  )
        co(2) = - el2/(16*pi2*st2)*(                            &
                  + mw2/4d0*( 1d0 + 1/(2*ct4) )*( 1d0 - 4*lam ) & ! bosonic
                  + mw2/4d0 + mh2/(8*ct2)                       & ! bosonic
                  + 1/mw2*( ml4(1) + ml4(2) + ml4(3) )          & ! leptonic
                  + Nc/mw2*(                                    & ! hadronic
                    + mu4(1) + mu4(2) + mu4(3)                  & ! hadronic
                    + md4(1) + md4(2) + md4(3)                  & ! hadronic
                    )                                           & ! hadronic
                  )
      end select
    end select

  case ( 'p+ p-','p- p+' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = dZp
!        co(2) = - el*dt/(2*st*mw1) + dmw2 ! old (according to 0709.1075)
        co(2) = - el*dt/(2*st*mw1)
      case (3)
        co(1) = - el2/(96*pi2*st2)*(                   &
                  + 1/2d0*( 1d0 + 1/(2*ct2) )          & ! bosonic
                  + 1/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  + Nc/mw2*(                           & ! hadronic
                    + mu2(1) + mu2(2) + mu2(3)         & ! hadronic
                    + md2(1) + md2(2) + md2(3)         & ! hadronic
                    )                                  & ! hadronic
                  )
        co(2) = - el2/(32*pi2*st2)*(                                &
                  + mw2/2d0*                                        & ! bosonic
                    ( 3d0 - 4*lam - 2/ct2 + ( 1/2d0 - 2*lam )/ct4 ) & ! bosonic
                  + mw2/(4*ct2) + 1/4d0*( mh2 + mz2 )               & ! bosonic
                  + 1/mw2*( ml4(1) + ml4(2) + ml4(3) )              & ! leptonic
                  + Nc/mw2*(                                        & ! hadronic
                    + ( mu2(1) + md2(1) )**2                        & ! hadronic
                    + ( mu2(2) + md2(2) )**2                        & ! hadronic
                    + ( mu2(3) + md2(3) )**2                        & ! hadronic
                    )                                               & ! hadronic
                  )
      end select
    end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 1 scalar, 1 vector
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! coupling = i*co(1)*p_mu
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! p  = momentum of the vector
! mu = Lorentz index of the vector
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  case ( 'p0 A' )
!
!    select case (ngs)
!    case (0)
!      select case (xlp)
!      case (2)
!        co(1) = - cId0/2d0*dZza * mw1/ct
!      end select
!    end select
!
!  case ( 'A p0' )
!
!    select case (ngs)
!    case (0)
!      select case (xlp)
!      case (2)
!        co(1) = cId0/2d0*dZza * mw1/ct
!      end select
!    end select
!
!  case ( 'p0 Z' )
!    select case (ngs)
!    case (0)
!      select case (xlp)
!      case (2)
!        co(1) = - cId0/2d0*( dZzz + dmz2/mz2 ) * mw1/ct
!      end select
!    end select
!
!  case ( 'Z p0' )
!
!    select case (ngs)
!    case (0)
!      select case (xlp)
!      case (2)
!        co(1) = cId0/2d0*( dZzz + dmz2/mz2 ) * mw1/ct
!      end select
!    end select
!
!  case ( 'p+ W-' )
!
!    select case (ngs)
!    case (0)
!      select case (xlp)
!      case (2)
!!        co(1) = 1/2d0*( dZw + dmw2/mw2 ) * mw1
!        co(1) = 1/2d0*( dZw + dZp + dmw2/mw2 ) * mw1
!      end select
!    end select
!
!  case ( 'W- p+' )
!
!    select case (ngs)
!    case (0)
!      select case (xlp)
!      case (2)
!!        co(1) = - 1/2d0*( dZw + dmw2/mw2 ) * mw1
!        co(1) = - 1/2d0*( dZw + dZp + dmw2/mw2 ) * mw1
!      end select
!    end select
!
!  case ( 'p- W+' )
!
!    select case (ngs)
!    case (0)
!      select case (xlp)
!      case (2)
!!        co(1) = - 1/2d0*( dZw + dmw2/mw2 ) * mw1
!        co(1) = - 1/2d0*( dZw + dZp + dmw2/mw2 ) * mw1
!      end select
!    end select
!
!  case ( 'W+ p-' )
!
!    select case (ngs)
!    case (0)
!      select case (xlp)
!      case (2)
!!        co(1) = 1/2d0*( dZw + dmw2/mw2 ) * mw1
!        co(1) = 1/2d0*( dZw + dZp + dmw2/mw2 ) * mw1
!      end select
!    end select
!
  case ( 'p0 A', 'A p0' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = cId0/2d0*dZza * mw1/ct
      end select
    end select

  case ( 'p0 Z', 'Z p0' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
!        co(1) = cId0/2d0*( dZzz + dmz2/mz2 ) * mw1/ct
        co(1) = cId0/2d0*( dZzz + dZp0 + dmz2/mz2 ) * mw1/ct
      end select
    end select

  case ( 'p+ W-', 'W- p+' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
!        co(1) = - 1/2d0*( dZw + dmw2/mw2 ) * mw1
        co(1) = - 1/2d0*( dZw + dZp + dmw2/mw2 ) * mw1
      end select
    end select

  case ( 'p- W+', 'W+ p-' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
!        co(1) = 1/2d0*( dZw + dmw2/mw2 ) * mw1
        co(1) = 1/2d0*( dZw + dZp + dmw2/mw2 ) * mw1
      end select
    end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 2 vectors
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! coupling = - i*(
!              + ( co(1)*p^2 - co(2) )*g_{mu1,mu2}
!              + co(3)*p_mu1*p_mu2
!              )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! p    = momentum of the 1st vector
! mu_1 = Lorentz index of the 1st vector
! mu_2 = Lorentz index of the 2nd vector
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  case ( 'g g' )

    select case (ngs)
    case (2)
      select case (xlp)
      case (2)
        co(1) = dZg
      case (3)
        co(1) = - gs2*Nc/(48*pi2)*( &
                  + 1/2d0 + lam   & ! ggg+ggg
                  + 1/Nc          & ! gu{u}+gu{u}
                  + 1/Nc          & ! gc{c}+gc{c}
                  + 1/Nc          & ! gt{t}+gt{t}
                  + 1/Nc          & ! gd{d}+gd{d}
                  + 1/Nc          & ! gs{s}+gs{s}
                  + 1/Nc          & ! gb{b}+gb{b}
                  )
        co(2) = - gs2/(8*pi2)*( &
                  + mu2(1)    & ! gu{u}+gu{u}
                  + mu2(2)    & ! gc{c}+gc{c}
                  + mu2(3)    & ! gt{t}+gt{t}
                  + md2(1)    & ! gd{d}+gd{d}
                  + md2(2)    & ! gs{s}+gs{s}
                  + md2(3)    & ! gb{b}+gb{b}
                  )
        co(3) = gs2*Nc/(48*pi2)*lam ! ggg+ggg
      end select
    end select

  case ( 'A A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = dZaa
        co(3) = - dZaa
      case (3)
        co(1) = - el2/(8*pi2)*(           &
                  + 1/6d0*( 1d0 + 2*lam ) & ! bosonic
                  + Ql2                   & ! leptonic
                  + Nc*( Qu2 + Qd2 )      & ! hadronic
                  )
        co(2) = - el2/(4*pi2)*(                        &
                  + mw2/2d0                            & ! bosonic
                  + Ql2*( ml2(1) + ml2(2) + ml2(3) )   & ! leptonic
                  + Nc*(                               & ! hadronic
                    + Qu2*( mu2(1) + mu2(2) + mu2(3) ) & ! hadronic
                    + Qd2*( md2(1) + md2(2) + md2(3) ) & ! hadronic
                    )                                  & ! hadronic
                  )
        co(3) = el2/(24*pi2)*lam ! bosonic
        cQED(1) = - el2/(8*pi2)*(           &
                    + Ql2                   & ! leptonic
                    + Nc*( Qu2 + Qd2 )      & ! hadronic
                    )
        cQED(2) = - el2/(4*pi2)*(                        &
                    + Ql2*( ml2(1) + ml2(2) + ml2(3) )   & ! leptonic
                    + Nc*(                               & ! hadronic
                      + Qu2*( mu2(1) + mu2(2) + mu2(3) ) & ! hadronic
                      + Qd2*( md2(1) + md2(2) + md2(3) ) & ! hadronic
                      )                                  & ! hadronic
                    )
      end select
    end select

  case ( 'A Z','Z A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZaz + dZza )
        co(2) = mz2/2d0*dZza
        co(3) = - 1/2d0*( dZaz + dZza )
      case (3)
        co(1) = + el2/(8*pi2)*(                  &
                  + ct/(6*st)*( 1d0 + 2*lam )    & ! bosonic
                  + Ql/ct*( I3l/(2*st) - Ql*st ) & ! leptonic
                  + Nc/ct*(                      & ! hadronic
                    + Qu*( I3u/(2*st) - Qu*st )  & ! hadronic
                    + Qd*( I3d/(2*st) - Qd*st )  & ! hadronic
                    )                            & ! hadronic
                  )
        co(2) = + el2/(4*pi2)*(                    &
                  + ct/(2*st)*mw2                  & ! bosonic
                  + Ql/ct*( I3l/(2*st) - Ql*st )*  & ! leptonic
                    ( ml2(1) + ml2(2) + ml2(3) )   & ! leptonic
                  + Nc/ct*(                        & ! hadronic
                    + Qu*( I3u/(2*st) - Qu*st )*   & ! hadronic
                      ( mu2(1) + mu2(2) + mu2(3) ) & ! hadronic
                    + Qd*( I3d/(2*st) - Qd*st )*   & ! hadronic
                      ( md2(1) + md2(2) + md2(3) ) & ! hadronic
                    )                              & ! hadronic
                  )
        co(3) = - el2/(24*pi2)*ct/st*lam ! bosonic
      end select
    end select

  case ( 'Z Z' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = dZzz
        co(2) = mz2*dZzz + dmz2
        co(3) = - dZzz
      case (3)
        co(1) = - el2/(8*pi2)*(                       &
                  + ct2/(6*st2)*( 1d0 + 2*lam )       & ! bosonic
                  + 1/ct2*(                           & ! leptonic
                    + I3n2/(2*st2)                    & ! leptonic
                    + I3l2/(2*st2) - Ql*I3l + Ql2*st2 & ! leptonic
                    )                                 & ! leptonic
                  + Nc/ct2*(                          & ! hadronic
                    + I3u2/(2*st2) - Qu*I3u + Qu2*st2 & ! hadronic
                    + I3d2/(2*st2) - Qd*I3d + Qd2*st2 & ! hadronic
                    )                                 & ! hadronic
                  )
        co(2) = - el2/(4*pi2)*(                                &
                  + ct2/(2*st2)*mw2                            & ! bosonic
                  + 1/ct2*( I3l2/(2*st2) - Ql*I3l + Ql2*st2 )* & ! leptonic
                    ( ml2(1) + ml2(2) + ml2(3) )               & ! leptonic
                  + Nc/ct2*(                                   & ! hadronic
                    + ( I3u2/(2*st2) - Qu*I3u + Qu2*st2 )*     & ! hadronic
                      ( mu2(1) + mu2(2) + mu2(3) )             & ! hadronic
                    + ( I3d2/(2*st2) - Qd*I3d + Qd2*st2 )*     & ! hadronic
                      ( md2(1) + md2(2) + md2(3) )             & ! hadronic
                    )                                          & ! hadronic
                  )
        co(3) = el2/(24*pi2)*ct2/st2*lam ! bosonic
      end select
    end select

  case ( 'W+ W-','W- W+' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = dZw
        co(2) = mw2*dZw + dmw2
        co(3) = - dZw
      case (3)
        co(1) = - el2/(8*pi2)*(               &
                  + 1/(6*st2)*( 1d0 + 2*lam ) & ! bosonic
                  + 1/(4*st2)*(               &
                    + 1d0                     & ! leptonic
                    + Nc                      & ! hadronic
                    )                         &
                  )
        co(2) = - el2/(8*pi2*st2)*(                    &
                  + mw2                                & ! bosonic
                  + 1/4d0*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  + 1/4d0*Nc*(                         & ! hadronic
                    + mu2(1) + mu2(2) + mu2(3)         & ! hadronic
                    + md2(1) + md2(2) + md2(3)         & ! hadronic
                    )                                  & ! hadronic
                  )
        co(3) = el2/(24*pi2)/st2*lam ! bosonic
      end select
    end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 2 fermions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! coupling = i*(
!            + co(1)*pslash*omega_-
!            + co(2)*pslash*omega_+
!            + co(3)*omega_-
!            + co(4)*omega_+
!            )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! pslash = p^mu*gamma_mu
! p      = momentum of the 1st fermion
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  case ( 'nu_e nu_e~','nu_e~ nu_e' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZnL(1) + dZnLc(1) )
        co(2) = 1/2d0*( dZnR(1) + dZnRc(1) )
      case (3)
        co(1) = - el2/(16*pi2)*( gmn2 + 1/(2*st2) )*lam
      end select
    end select

  case ( 'nu_mu nu_mu~','nu_mu~ nu_mu' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZnL(2) + dZnLc(2) )
        co(2) = 1/2d0*( dZnR(2) + dZnRc(2) )
      case (3)
        co(1) = - el2/(16*pi2)*( gmn2 + 1/(2*st2) )*lam
      end select
    end select

  case ( 'nu_tau nu_tau~','nu_tau~ nu_tau' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZnL(3) + dZnLc(3) )
        co(2) = 1/2d0*( dZnR(3) + dZnRc(3) )
      case (3)
        co(1) = - el2/(16*pi2)*( gmn2 + 1/(2*st2) )*lam
      end select
    end select

  case ( 'u u~','u~ u' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZuL(1) + dZuLc(1) )
        co(2) = 1/2d0*( dZuR(1) + dZuRc(1) )
        co(3) = - mu1(1)/2d0*( dZuL(1) + dZuRc(1) ) - dmu(1)
        co(4) = - mu1(1)/2d0*( dZuR(1) + dZuLc(1) ) - dmu(1)
      case (3)
        cP(1)   = - el2/(16*pi2)*Qu2*lam
        cP(2)   = - el2/(16*pi2)*Qu2*lam
        cP(3:4) = + el2/(8*pi2)*Qu2*mu1(1)*lam
        co(1)   = - el2/(16*pi2)*( gmu2 + 1/(2*st2) )*lam
        co(2)   = - el2/(16*pi2)*gpu2*lam
        co(3:4) = + el2/(8*pi2)*gpu*gmu*mu1(1)*lam
        cQED = cP
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZuLqcd(1) + dZuLcqcd(1) )
        co(2) = 1/2d0*( dZuRqcd(1) + dZuRcqcd(1) )
        co(3) = - mu1(1)/2d0*( dZuLqcd(1) + dZuRcqcd(1) ) - dmuqcd(1)
        co(4) = - mu1(1)/2d0*( dZuRqcd(1) + dZuLcqcd(1) ) - dmuqcd(1)
      case (3)
        co(1)   = - gs2/(16*pi2)*Cf*lam
        co(2)   = - gs2/(16*pi2)*Cf*lam
        co(3:4) = + gs2/(8*pi2)*Cf*mu1(1)*lam
      end select
    end select

  case ( 'c c~','c~ c' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZuL(2) + dZuLc(2) )
        co(2) = 1/2d0*( dZuR(2) + dZuRc(2) )
        co(3) = - mu1(2)/2d0*( dZuL(2) + dZuRc(2) ) - dmu(2)
        co(4) = - mu1(2)/2d0*( dZuR(2) + dZuLc(2) ) - dmu(2)
      case (3)
        cP(1)   = - el2/(16*pi2)*Qu2*lam
        cP(2)   = - el2/(16*pi2)*Qu2*lam
        cP(3:4) = + el2/(8*pi2)*Qu2*mu1(2)*lam
        co(1)   = - el2/(16*pi2)*( gmu2 + 1/(2*st2) )*lam
        co(2)   = - el2/(16*pi2)*gpu2*lam
        co(3:4) = + el2/(8*pi2)*gpu*gmu*mu1(2)*lam
        cQED = cP
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZuLqcd(2) + dZuLcqcd(2) )
        co(2) = 1/2d0*( dZuRqcd(2) + dZuRcqcd(2) )
        co(3) = - mu1(2)/2d0*( dZuLqcd(2) + dZuRcqcd(2) ) - dmuqcd(2)
        co(4) = - mu1(2)/2d0*( dZuRqcd(2) + dZuLcqcd(2) ) - dmuqcd(2)
      case (3)
        co(1)   = - gs2/(16*pi2)*Cf*lam
        co(2)   = - gs2/(16*pi2)*Cf*lam
        co(3:4) = + gs2/(8*pi2)*Cf*mu1(2)*lam
      end select
    end select

  case ( 't t~','t~ t' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZuL(3) + dZuLc(3) )
        co(2) = 1/2d0*( dZuR(3) + dZuRc(3) )
        co(3) = - mu1(3)/2d0*( dZuL(3) + dZuRc(3) ) - dmu(3)
        co(4) = - mu1(3)/2d0*( dZuR(3) + dZuLc(3) ) - dmu(3)
      case (3)
        cP(1)   = - el2/(16*pi2)*Qu2*lam
        cP(2)   = - el2/(16*pi2)*Qu2*lam
        cP(3:4) = + el2/(8*pi2)*Qu2*mu1(3)*lam
        co(1)   = - el2/(16*pi2)*( gmu2 + 1/(2*st2) )*lam
        co(2)   = - el2/(16*pi2)*gpu2*lam
        co(3:4) = + el2/(8*pi2)*gpu*gmu*mu1(3)*lam
        cQED = cP
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZuLqcd(3) + dZuLcqcd(3) )
        co(2) = 1/2d0*( dZuRqcd(3) + dZuRcqcd(3) )
        co(3) = - mu1(3)/2d0*( dZuLqcd(3) + dZuRcqcd(3) ) - dmuqcd(3)
        co(4) = - mu1(3)/2d0*( dZuRqcd(3) + dZuLcqcd(3) ) - dmuqcd(3)
      case (3)
        co(1)   = - gs2/(16*pi2)*Cf*lam
        co(2)   = - gs2/(16*pi2)*Cf*lam
        co(3:4) = + gs2/(8*pi2)*Cf*mu1(3)*lam
      end select
    end select

  case ( 'e- e+','e+ e-' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1)   = 1/2d0*( dZlL(1) + dZlLc(1) )
        co(2)   = 1/2d0*( dZlR(1) + dZlRc(1) )
        co(3)   = - ml1(1)/2d0*( dZlL(1) + dZlRc(1) ) - dml(1)
        co(4)   = - ml1(1)/2d0*( dZlR(1) + dZlLc(1) ) - dml(1)
      case (3)
        cP(1)   = - el2/(16*pi2)*Ql2*lam
        cP(2)   = - el2/(16*pi2)*Ql2*lam
        cP(3:4) = + el2/(8*pi2)*Ql2*ml1(1)*lam
        co(1)   = - el2/(16*pi2)*( gml2 + 1/(2*st2) )*lam
        co(2)   = - el2/(16*pi2)*gpl2*lam
        co(3:4) = + el2/(8*pi2)*gpl*gml*ml1(1)*lam
        cQED = cP
      end select
    end select

  case ( 'mu- mu+','mu+ mu-' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1)   = 1/2d0*( dZlL(2) + dZlLc(2) )
        co(2)   = 1/2d0*( dZlR(2) + dZlRc(2) )
        co(3)   = - ml1(2)/2d0*( dZlL(2) + dZlRc(2) ) - dml(2)
        co(4)   = - ml1(2)/2d0*( dZlR(2) + dZlLc(2) ) - dml(2)
      case (3)
        cP(1)   = - el2/(16*pi2)*Ql2*lam
        cP(2)   = - el2/(16*pi2)*Ql2*lam
        cP(3:4) = + el2/(8*pi2)*Ql2*ml1(2)*lam
        co(1)   = - el2/(16*pi2)*( gml2 + 1/(2*st2) )*lam
        co(2)   = - el2/(16*pi2)*gpl2*lam
        co(3:4) = + el2/(8*pi2)*gpl*gml*ml1(2)*lam
        cQED = cP
      end select
    end select

  case ( 'tau- tau+','tau+ tau-' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1)   = 1/2d0*( dZlL(3) + dZlLc(3) )
        co(2)   = 1/2d0*( dZlR(3) + dZlRc(3) )
        co(3)   = - ml1(3)/2d0*( dZlL(3) + dZlRc(3) ) - dml(3)
        co(4)   = - ml1(3)/2d0*( dZlR(3) + dZlLc(3) ) - dml(3)
      case (3)
        cP(1)   = - el2/(16*pi2)*Ql2*lam
        cP(2)   = - el2/(16*pi2)*Ql2*lam
        cP(3:4) = + el2/(8*pi2)*Ql2*ml1(3)*lam
        co(1)   = - el2/(16*pi2)*( gml2 + 1/(2*st2) )*lam
        co(2)   = - el2/(16*pi2)*gpl2*lam
        co(3:4) = + el2/(8*pi2)*gpl*gml*ml1(3)*lam
        cQED = cP
      end select
    end select

  case ( 'd d~','d~ d' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZdL(1) + dZdLc(1) )
        co(2) = 1/2d0*( dZdR(1) + dZdRc(1) )
        co(3) = - md1(1)/2d0*( dZdL(1) + dZdRc(1) ) - dmd(1)
        co(4) = - md1(1)/2d0*( dZdR(1) + dZdLc(1) ) - dmd(1)
      case (3)
        cP(1)   = - el2/(16*pi2)*Qd2*lam
        cP(2)   = - el2/(16*pi2)*Qd2*lam
        cP(3:4) = + el2/(8*pi2)*Qd2*md1(1)*lam
        co(1)   = - el2/(16*pi2)*( gmd2 + 1/(2*st2) )*lam
        co(2)   = - el2/(16*pi2)*gpd2*lam
        co(3:4) = + el2/(8*pi2)*gpd*gmd*md1(1)*lam
        cQED = cP
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZdLqcd(1) + dZdLcqcd(1) )
        co(2) = 1/2d0*( dZdRqcd(1) + dZdRcqcd(1) )
        co(3) = - md1(1)/2d0*( dZdLqcd(1) + dZdRcqcd(1) ) - dmdqcd(1)
        co(4) = - md1(1)/2d0*( dZdRqcd(1) + dZdLcqcd(1) ) - dmdqcd(1)
      case (3)
        co(1)   = - gs2/(16*pi2)*Cf*lam
        co(2)   = - gs2/(16*pi2)*Cf*lam
        co(3:4) = + gs2/(8*pi2)*Cf*md1(1)*lam
      end select
    end select

  case ( 's s~','s~ s' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZdL(2) + dZdLc(2) )
        co(2) = 1/2d0*( dZdR(2) + dZdRc(2) )
        co(3) = - md1(2)/2d0*( dZdL(2) + dZdRc(2) ) - dmd(2)
        co(4) = - md1(2)/2d0*( dZdR(2) + dZdLc(2) ) - dmd(2)
      case (3)
        cP(1)   = - el2/(16*pi2)*Qd2*lam
        cP(2)   = - el2/(16*pi2)*Qd2*lam
        cP(3:4) = + el2/(8*pi2)*Qd2*md1(2)*lam
        co(1)   = - el2/(16*pi2)*( gmd2 + 1/(2*st2) )*lam
        co(2)   = - el2/(16*pi2)*gpd2*lam
        co(3:4) = + el2/(8*pi2)*gpd*gmd*md1(2)*lam
        cQED = cP
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZdLqcd(2) + dZdLcqcd(2) )
        co(2) = 1/2d0*( dZdRqcd(2) + dZdRcqcd(2) )
        co(3) = - md1(2)/2d0*( dZdLqcd(2) + dZdRcqcd(2) ) - dmdqcd(2)
        co(4) = - md1(2)/2d0*( dZdRqcd(2) + dZdLcqcd(2) ) - dmdqcd(2)
      case (3)
        co(1)   = - gs2/(16*pi2)*Cf*lam
        co(2)   = - gs2/(16*pi2)*Cf*lam
        co(3:4) = + gs2/(8*pi2)*Cf*md1(2)*lam
      end select
    end select

  case ( 'b b~','b~ b' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZdL(3) + dZdLc(3) )
        co(2) = 1/2d0*( dZdR(3) + dZdRc(3) )
        co(3) = - md1(3)/2d0*( dZdL(3) + dZdRc(3) ) - dmd(3)
        co(4) = - md1(3)/2d0*( dZdR(3) + dZdLc(3) ) - dmd(3)
      case (3)
        cP(1)   = - el2/(16*pi2)*Qd2*lam
        cP(2)   = - el2/(16*pi2)*Qd2*lam
        cP(3:4) = + el2/(8*pi2)*Qd2*md1(3)*lam
        co(1)   = - el2/(16*pi2)*( gmd2 + 1/(2*st2) )*lam
        co(2)   = - el2/(16*pi2)*gpd2*lam
        co(3:4) = + el2/(8*pi2)*gpd*gmd*md1(3)*lam
        cQED = cP
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = 1/2d0*( dZdLqcd(3) + dZdLcqcd(3) )
        co(2) = 1/2d0*( dZdRqcd(3) + dZdRcqcd(3) )
        co(3) = - md1(3)/2d0*( dZdLqcd(3) + dZdRcqcd(3) ) - dmdqcd(3)
        co(4) = - md1(3)/2d0*( dZdRqcd(3) + dZdLcqcd(3) ) - dmdqcd(3)
      case (3)
        co(1)   = - gs2/(16*pi2)*Cf*lam
        co(2)   = - gs2/(16*pi2)*Cf*lam
        co(3:4) = + gs2/(8*pi2)*Cf*md1(3)*lam
      end select
    end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 2 FP ghost fields and 1 scalar
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! coupling = i*co(1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  case ( 'X_A X+~ p+','X_A p+ X+~','X+~ X_A p+', &
         'X+~ p+ X_A','p+ X_A X+~','p+ X+~ X_A', &
         'X_A X-~ p-','X_A p- X-~','X-~ X_A p-', &
         'X-~ p- X_A','p- X_A X-~','p- X-~ X_A'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*mw1
      end select
    end select

  case ( 'X_Z X_Z~ H','X_Z H X_Z~','X_Z~ X_Z H', &
         'X_Z~ H X_Z','H X_Z X_Z~','H X_Z~ X_Z'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*mw1/(2*st*ct2)
      end select
    end select

  case ( 'X_Z X+~ p+','X_Z p+ X+~','X+~ X_Z p+', &
         'X+~ p+ X_Z','p+ X_Z X+~','p+ X+~ X_Z', &
         'X_Z X-~ p-','X_Z p- X-~','X-~ X_Z p-', &
         'X-~ p- X_Z','p- X_Z X-~','p- X-~ X_Z'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*mw1*(st2-ct2)/(2*stct)
      end select
    end select

  case ( 'X+ X_Z~ p-','X+ p- X_Z~','X_Z~ X+ p-', &
         'X_Z~ p- X+','p- X+ X_Z~','p- X_Z~ X+', &
         'X- X_Z~ p+','X- p+ X_Z~','X_Z~ X- p+', &
         'X_Z~ p+ X-','p+ X- X_Z~','p+ X_Z~ X-'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*mw1/(2*stct)
      end select
    end select

  case ( 'X+ X+~ H','X+ H X+~','X+~ X+ H', &
         'X+~ H X+','H X+ X+~','H X+~ X+', &
         'X- X-~ H','X- H X-~','X-~ X- H', &
         'X-~ H X-','H X- X-~','H X-~ X-'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*mw1/(2*st)
      end select
    end select

  case ( 'X+ X+~ p0','X+ p0 X+~','X+~ X+ p0', &
         'X+~ p0 X+','p0 X+ X+~','p0 X+~ X+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - cId0*el*mw1/(2*st)
      end select
    end select

  case ( 'X- X-~ p0','X- p0 X-~','X-~ X- p0', &
         'X-~ p0 X-','p0 X- X-~','p0 X-~ X-'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = cId0*el*mw1/(2*st)
      end select
    end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 2 FP ghost fields and 1 vector
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! coupling = i*co(1)*p_a_mu
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! p_a = momentum of the anti-ghost field
! mu  = Lorentz index of the vector
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  case ( 'X_g X_g~ g','X_g g X_g~','X_g~ X_g g', &
         'X_g~ g X_g','g X_g X_g~','g X_g~ X_g'  )

    select case (ngs)
    case (1)
      select case (xlp)
      case (0)
        co(1) = gs/csq2
      end select
    end select

  case ( 'X_A X+~ W+','X_A W+ X+~','X+~ X_A W+', &
         'X+~ W+ X_A','W+ X_A X+~','W+ X+~ X_A'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el
      end select
    end select

  case ( 'X_A X-~ W-','X_A W- X-~','X-~ X_A W-', &
         'X-~ W- X_A','W- X_A X-~','W- X-~ X_A'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el
      end select
    end select

  case ( 'X_Z X+~ W+','X_Z W+ X+~','X+~ X_Z W+', &
         'X+~ W+ X_Z','W+ X_Z X+~','W+ X+~ X_Z'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*ct/st
      end select
    end select

  case ( 'X_Z X-~ W-','X_Z W- X-~','X-~ X_Z W-', &
         'X-~ W- X_Z','W- X_Z X-~','W- X-~ X_Z'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*ct/st
      end select
    end select

  case ( 'X+ X_A~ W-','X+ W- X_A~','X_A~ X+ W-', &
         'X_A~ W- X+','W- X+ X_A~','W- X_A~ X+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el
      end select
    end select

  case ( 'X+ X_Z~ W-','X+ W- X_Z~','X_Z~ X+ W-', &
         'X_Z~ W- X+','W- X+ X_Z~','W- X_Z~ X+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*ct/st
      end select
    end select

  case ( 'X+ X+~ A','X+ A X+~','X+~ X+ A', &
         'X+~ A X+','A X+ X+~','A X+~ X+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el
      end select
    end select

  case ( 'X+ X+~ Z','X+ Z X+~','X+~ X+ Z', &
         'X+~ Z X+','Z X+ X+~','Z X+~ X+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*ct/st
      end select
    end select

  case ( 'X- X_A~ W+','X- W+ X_A~','X_A~ X- W+', &
         'X_A~ W+ X-','W+ X- X_A~','W+ X_A~ X-'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el
      end select
    end select

  case ( 'X- X_Z~ W+','X- W+ X_Z~','X_Z~ X- W+', &
         'X_Z~ W+ X-','W+ X- X_Z~','W+ X_Z~ X-'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*ct/st
      end select
    end select

  case ( 'X- X-~ A','X- A X-~','X-~ X- A', &
         'X-~ A X-','A X- X-~','A X-~ X-'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el
      end select
    end select

  case ( 'X- X-~ Z','X- Z X-~','X-~ X- Z', &
         'X-~ Z X-','Z X- X-~','Z X-~ X-'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*ct/st
      end select
    end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 3 scalars
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! coupling = i*co(1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  case ( 'H H H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*3/2d0*mh2/(st*mw1)
      case (2)
        co(1) = - el*3/2d0*mh2/(st*mw1)*(    &
                  + dZe - dst/st + dmh2/mh2  &
                  - dmw2/(2*mw2) + 3/2d0*dZh &
                  )                          &
                - el2*3/(4*st2*mw2)*dt
      case (3)
        co(1) = el3/pi2*3/(32*st3*mw1)*(                  &
                + mw2/4d0*( 2d0 + 1/ct4 )*( 1d0 - 4*lam ) & ! bosonic
                + mh2/4d0*( 1d0 + 1/(2*ct2) )             & ! bosonic
                + 1/mw2*( ml4(1) + ml4(2) + ml4(3) )      & ! leptonic
                + Nc/mw2*(                                & ! hadronic
                  + md4(1) + md4(2) + md4(3)              & ! hadronic
                  + mu4(1) + mu4(2) + mu4(3)              & ! hadronic
                  )                                       & ! hadronic
                )
      end select
    end select

  case ( 'H p0 p0','p0 H p0','p0 p0 H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*mh2/(2*st*mw1)
      case (2)
        co(1) = - el*mh2/(2*st*mw1)*(        &
                  + dZe - dst/st + dmh2/mh2  &
                  - dmw2/(2*mw2) + 1/2d0*dZh &
                  + dZp0                     &
                  )                          &
                - el2/(4*st2*mw2)*dt
      case (3)
        co(1) = el3/pi2/(32*st3*mw1)*(                    &
                + mw2/4d0*( 2d0 + 1/ct4 )*( 1d0 - 4*lam ) & ! bosonic
                + mh2/4d0*( 1d0 + 1/(2*ct2) )             & ! bosonic
                + 1/mw2*( ml4(1) + ml4(2) + ml4(3) )      & ! leptonic
                + Nc/mw2*(                                & ! hadronic
                  + md4(1) + md4(2) + md4(3)              & ! hadronic
                  + mu4(1) + mu4(2) + mu4(3)              & ! hadronic
                  )                                       & ! hadronic
                )
      end select
    end select

  case ( 'H p+ p-','H p- p+','p+ H p-','p+ p- H','p- H p+','p- p+ H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*mh2/(2*st*mw1)
      case (2)
        co(1) = - el*mh2/(2*st*mw1)*(      &
                + dZe - dst/st + dmh2/mh2  &
                - dmw2/(2*mw2) + 1/2d0*dZh &
                + dZp                      &
                )                          &
                - el2/(4*st2*mw2)*dt
      case (3)
        co(1) = el3/pi2/(32*st3*mw1)*(                 &
                + mw2/4d0*( 3d0 + st2/ct4 + st2/ct2 )* & ! bosonic
                  ( 1d0 - 4*lam )                      & ! bosonic
                + mh2/4d0*( 1d0 + 1/(2*ct2) )          & ! bosonic
                + 1/mw2*( ml4(1) + ml4(2) + ml4(3) )   & ! leptonic
                + Nc/mw2*(                             & ! hadronic
                  + md4(1) + md4(2) + md4(3)           & ! hadronic
                  + mu4(1) + mu4(2) + mu4(3)           & ! hadronic
                  )                                    & ! hadronic
                )
      end select
    end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 2 scalars and 1 vector
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! coupling = i*co(1)*(p_s1-p_s2)_mu
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! p_s1 = momentum of the 1st scalar
! p_s2 = momentum of the 2nd scalar
! mu   = Lorentz index of the vector
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  case ( 'H p0 A','A H p0','p0 A H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = - cId0*el/(4*stct)*dZza
      case (3)
        co(1) = - cId0*el3/pi2*5/(192*st2)
      end select
    end select

  case ( 'H A p0','A p0 H','p0 H A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = cId0*el/(4*stct)*dZza
      case (3)
        co(1) = cId0*el3/pi2*5/(192*st2)
      end select
    end select

  case ( 'H p0 Z','Z H p0','p0 Z H' )
    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - cId0*el/(2*stct)
      case (2)
        co(1) = - cId0*el/(2*stct)*(             &
                  + dZe + (st2-ct2)/(st*ct2)*dst &
                  + 1/2d0*( dZh + dZzz + dZp0 )  &
                  )
      case (3)
        co(1) = cId0*el3/pi2/(96*stct)*(                   &
                + 1/(8*st2*ct2)*( 1d0 + 2*ct2 + 20*ct4 )   & ! bosonic
                + 1/(mw2*st2)*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                + Nc/(mw2*st2)*(                           & ! hadronic
                  + mu2(1) + mu2(2) + mu2(3)               & ! hadronic
                  + md2(1) + md2(2) + md2(3)               & ! hadronic
                  )                                        & ! hadronic
                )
      end select
    end select

  case ( 'H Z p0','Z p0 H','p0 H Z' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = cId0*el/(2*stct)
      case (2)
        co(1) = cId0*el/(2*stct)*(             &
                + dZe + (st2-ct2)/(st*ct2)*dst &
                + 1/2d0*( dZh + dZzz + dZp0 )  &
                )
      case (3)
        co(1) = - cId0*el3/pi2/(96*stct)*(                   &
                  + 1/(8*st2*ct2)*( 1d0 + 2*ct2 + 20*ct4 )   & ! bosonic
                  + 1/(mw2*st2)*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  + Nc/(mw2*st2)*(                           & ! hadronic
                    + mu2(1) + mu2(2) + mu2(3)               & ! hadronic
                    + md2(1) + md2(2) + md2(3)               & ! hadronic
                    )                                        & ! hadronic
                  )
      end select
    end select

  case ( 'H p+ W-','p+ W- H','W- H p+' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el/(2*st)
      case (2)
        co(1) = el/(2*st)*( dZe - dst/st + 1/2d0*( dZw + dZh + dZp ) )
      case (3)
        co(1) = - el3/pi2/(96*st3)*(                   &
                  + 1/(8*ct2)*( 1d0 + 22*ct2 )         & ! bosonic
                  + 1/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  + Nc/mw2*(                           & ! hadronic
                    + mu2(1) + mu2(2) + mu2(3)         & ! hadronic
                    + md2(1) + md2(2) + md2(3)         & ! hadronic
                    )                                  & ! hadronic
                  )
      end select
    end select

  case ( 'H W- p+','p+ H W-','W- p+ H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el/(2*st)
      case (2)
        co(1) = - el/(2*st)*( dZe - dst/st + 1/2d0*( dZw + dZh + dZp ) )
      case (3)
        co(1) = el3/pi2/(96*st3)*(                   &
                + 1/(8*ct2)*( 1d0 + 22*ct2 )         & ! bosonic
                + 1/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                + Nc/mw2*(                           & ! hadronic
                  + mu2(1) + mu2(2) + mu2(3)         & ! hadronic
                  + md2(1) + md2(2) + md2(3)         & ! hadronic
                  )                                  & ! hadronic
                )
      end select
    end select

  case ( 'H p- W+','p- W+ H','W+ H p-' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el/(2*st)
      case (2)
        co(1) = - el/(2*st)*( dZe - dst/st + 1/2d0*( dZw + dZh + dZp ) )
      case (3)
        co(1) = el3/pi2/(96*st3)*(                   &
                + 1/(8*ct2)*( 1d0 + 22*ct2 )         & ! bosonic
                + 1/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                + Nc/mw2*(                           & ! hadronic
                  + mu2(1) + mu2(2) + mu2(3)         & ! hadronic
                  + md2(1) + md2(2) + md2(3)         & ! hadronic
                  )                                  & ! hadronic
                )
      end select
    end select

  case ( 'H W+ p-','p- H W+','W+ p- H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el/(2*st)
      case (2)
        co(1) = el/(2*st)*( dZe - dst/st + 1/2d0*( dZw + dZh + dZp ) )
      case (3)
        co(1) = - el3/pi2/(96*st3)*(                   &
                  + 1/(8*ct2)*( 1d0 + 22*ct2 )         & ! bosonic
                  + 1/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  + Nc/mw2*(                           & ! hadronic
                    + mu2(1) + mu2(2) + mu2(3)         & ! hadronic
                    + md2(1) + md2(2) + md2(3)         & ! hadronic
                    )                                  & ! hadronic
                  )
      end select
    end select

  case ( 'p0 p+ W-','p+ W- p0','W- p0 p+' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - cId0*el/(2*st)
      case (2)
        co(1) = - cId0*el/(2*st)*( dZe - dst/st + 1/2d0*( dZw + dZp + dZp0 ) )
      case (3)
        co(1) = - cId0*el3/pi2/(48*st3)*(                &
                  - 1/(16*ct2)*( 1d0 + 22*ct2 )          & ! bosonic
                  + I3l/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  - Nc/mw2*(                             & ! hadronic
                    + I3u*( mu2(1) + mu2(2) + mu2(3) )   & ! hadronic
                    - I3d*( md2(1) + md2(2) + md2(3) )   & ! hadronic
                    )                                    & ! hadronic
                  )
      end select
    end select

  case ( 'p0 W- p+','p+ p0 W-','W- p+ p0' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = cId0*el/(2*st)
      case (2)
        co(1) = cId0*el/(2*st)*( dZe - dst/st + 1/2d0*( dZw + dZp + dZp0 ) )
      case (3)
        co(1) = cId0*el3/pi2/(48*st3)*(                &
                - 1/(16*ct2)*( 1d0 + 22*ct2 )          & ! bosonic
                + I3l/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                - Nc/mw2*(                             & ! hadronic
                  + I3u*( mu2(1) + mu2(2) + mu2(3) )   & ! hadronic
                  - I3d*( md2(1) + md2(2) + md2(3) )   & ! hadronic
                  )                                    & ! hadronic
                )
      end select
    end select

  case ( 'p0 p- W+','p- W+ p0','W+ p0 p-' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - cId0*el/(2*st)
      case (2)
        co(1) = - cId0*el/(2*st)*( dZe - dst/st + 1/2d0*( dZw + dZp + dZp0 ) )
      case (3)
        co(1) = - cId0*el3/pi2/(48*st3)*(                &
                  - 1/(16*ct2)*( 1d0 + 22*ct2 )          & ! bosonic
                  + I3l/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  - Nc/mw2*(                             & ! hadronic
                    + I3u*( mu2(1) + mu2(2) + mu2(3) )   & ! hadronic
                    - I3d*( md2(1) + md2(2) + md2(3) )   & ! hadronic
                    )                                    & ! hadronic
                  )
      end select
    end select

  case ( 'p0 W+ p-','p- p0 W+','W+ p- p0' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = cId0*el/(2*st)
      case (2)
        co(1) = cId0*el/(2*st)*( dZe - dst/st + 1/2d0*( dZw + dZp + dZp0 ) )
      case (3)
        co(1) = cId0*el3/pi2/(48*st3)*(                &
                - 1/(16*ct2)*( 1d0 + 22*ct2 )          & ! bosonic
                + I3l/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                - Nc/mw2*(                             & ! hadronic
                  + I3u*( mu2(1) + mu2(2) + mu2(3) )   & ! hadronic
                  - I3d*( md2(1) + md2(2) + md2(3) )   & ! hadronic
                  )                                    & ! hadronic
                )
      end select
    end select

  case ( 'p+ A p-','p- p+ A','A p- p+' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el
      case (2)
        co(1) = - el*( dZe + 1/2d0*dZaa + (st2-ct2)/(4*stct)*dZza + dZp )
      case (3)
        co(1) = el3/pi2/(48*st2)*(                    &
                + 1/(8*ct2)*( 1d0 + 12*ct2 )          & ! bosonic
                - Ql/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                + Nc/mw2*(                            & ! hadronic
                  + mu2(1) + mu2(2) + mu2(3)          & ! hadronic
                  + md2(1) + md2(2) + md2(3)          & ! hadronic
                  )                                   & ! hadronic
                )
      end select
    end select

  case ( 'p+ p- A','p- A p+','A p+ p-' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el
      case (2)
        co(1) = el*( dZe + 1/2d0*dZaa + (st2-ct2)/(4*stct)*dZza + dZp )
      case (3)
        co(1) = - el3/pi2/(48*st2)*(                    &
                  + 1/(8*ct2)*( 1d0 + 12*ct2 )          & ! bosonic
                  - Ql/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  + Nc/mw2*(                            & ! hadronic
                    + mu2(1) + mu2(2) + mu2(3)          & ! hadronic
                    + md2(1) + md2(2) + md2(3)          & ! hadronic
                    )                                   & ! hadronic
                  )
      end select
    end select

  case ( 'p+ Z p-','p- p+ Z','Z p- p+' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*(st2-ct2)/(2*stct)
      case (2)
        co(1) = - el*(st2-ct2)/(2*stct)*( dZe + 1/2d0*dZzz + dZp ) &
                - el/(2*ct3*st2)*dst - el/2d0*dZaz
      case (3)
        co(1) = el3/pi2/(48*stct)*(                                &
                + 1/(16*ct2*st2)*( 1d0 - 24*ct4 )                  & ! bosonic
                - 1/mw2*( Ql + I3n/st2 )*                          & ! leptonic
                  ( ml2(1) + ml2(2) + ml2(3) )                     & ! leptonic
                + Nc/mw2*(                                         & ! hadronic
                  + ( 1d0 + I3d/st2 )*( mu2(1) + mu2(2) + mu2(3) ) & ! hadronic
                  + ( 1d0 - I3u/st2 )*( md2(1) + md2(2) + md2(3) ) & ! hadronic
                  )                                                & ! hadronic
                )
      end select
    end select

  case ( 'p+ p- Z','p- Z p+','Z p+ p-' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*(st2-ct2)/(2*stct)
      case (2)
        co(1) = + el*(st2-ct2)/(2*stct)*( dZe + 1/2d0*dZzz + dZp ) &
                + el/(2*ct3*st2)*dst + el/2d0*dZaz
      case (3)
        co(1) = - el3/pi2/(48*stct)*(                                &
                  + 1/(16*ct2*st2)*( 1d0 - 24*ct4 )                  & ! bosonic
                  - 1/mw2*( Ql + I3n/st2 )*                          & ! leptonic
                    ( ml2(1) + ml2(2) + ml2(3) )                     & ! leptonic
                  + Nc/mw2*(                                         & ! hadronic
                    + ( 1d0 + I3d/st2 )*( mu2(1) + mu2(2) + mu2(3) ) & ! hadronic
                    + ( 1d0 - I3u/st2 )*( md2(1) + md2(2) + md2(3) ) & ! hadronic
                    )                                                & ! hadronic
                  )
      end select
    end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 1 scalar and 2 vectors
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! coupling = i*co(1)*g_{mu_v1,mu_v2}
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mu_v1 = Lorentz index of the 1st vector
! mu_v2 = Lorentz index of the 2nd vector
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  case ( 'H g g','g H g','g g H' )

    select case (ngs)
    case (2)
      select case (xlp)
      case (3)
        co(1) = - el*gs2/(16*pi2*st*mw1)*( &
                  + mu2(1) + mu2(2) + mu2(3)    & ! hadronic
                  + md2(1) + md2(2) + md2(3)    & ! hadronic
                  )
       end select
    end select

  case ( 'H A A','A H A','A A H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (3)
        co(1) = - el3/pi2/(8*st*mw1)*(                 &
                  + mw2/2d0                            & ! bosonic
                  + Ql2*( ml2(1) + ml2(2) + ml2(3) )   & ! leptonic
                  + Nc*(                               & ! hadronic
                    + Qu2*( mu2(1) + mu2(2) + mu2(3) ) & ! hadronic
                    + Qd2*( md2(1) + md2(2) + md2(3) ) & ! hadronic
                    )                                  & ! hadronic
                  )
      end select
    end select

  case ( 'H A Z','H Z A','A H Z','A Z H','Z A H','Z H A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = el*mw1/(2*st*ct2)*dZza
      case (3)
        co(1) = el3/pi2/(8*ct*mw1)*(             &
                + mw2/(4*st2)*( 1d0 + 2*ct2 )    & ! bosonic
                + ( Ql*I3l/(2*st2) - Ql2 )*      & ! leptonic
                  ( ml2(1) + ml2(2) + ml2(3) )   & ! leptonic
                + Nc*(                           & ! hadronic
                  + ( Qu*I3u/(2*st2) - Qu2 )*    & ! hadronic
                    ( mu2(1) + mu2(2) + mu2(3) ) & ! hadronic
                  + ( Qd*I3d/(2*st2) - Qd2 )*    & ! hadronic
                    ( md2(1) + md2(2) + md2(3) ) & ! hadronic
                  )                              & ! hadronic
                )
      end select
    end select

  case ( 'H Z Z','Z H Z','Z Z H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*mw1/(st*ct2)
      case (2)
        co(1) = el*mw1/(st*ct2)*(                 &
                + dZe + dst*(2*st2-ct2)/(st*ct2)  &
                + 1/2d0*( dmw2/mw2 + dZh ) + dZzz &
                )
      case (3)
        co(1) = el3/(8*pi2*mw1)*(                          &
                + mw2/(2*st)*( 1d0 - 2/st2 )               & ! bosonic
                + 1/ct2*( Ql*I3l/st - Ql2*st - I3l2/st3 )* & ! leptonic
                  ( ml2(1) + ml2(2) + ml2(3) )             & ! leptonic
                + Nc/ct2*(                                 & ! hadronic
                  + ( Qu*I3u/st - Qu2*st - I3u2/st3 )*     & ! hadronic
                    ( mu2(1) + mu2(2) + mu2(3) )           & ! hadronic
                  + ( Qd*I3d/st - Qd2*st - I3d2/st3 )*     & ! hadronic
                    ( md2(1) + md2(2) + md2(3) )           & ! hadronic
                  )                                        & ! hadronic
                )
      end select
    end select

  case ( 'H W+ W-','H W- W+','W+ H W-','W+ W- H','W- H W+','W- W+ H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*mw1/st
      case (2)
        co(1) = el*mw1/st*(                                     &
                + dZe - dst/st + 1/2d0*( dmw2/mw2 + dZh ) + dZw &
                )
      case (3)
        co(1) = - el3/pi2/(8*st3*mw1)*(                &
                  + mw2                                & ! bosonic
                  + 1/4d0*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  + 1/4d0*Nc*(                         & ! hadronic
                    + mu2(1) + mu2(2) + mu2(3)         & ! hadronic
                    + md2(1) + md2(2) + md2(3)         & ! hadronic
                    )                                  & ! hadronic
                  )
      end select
    end select

  case ( 'p+ A W-','p+ W- A','W- p+ A','W- A p+','A p+ W-','A W- p+' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*mw1
      case (2)
        co(1) = - el*mw1*(                                    &
                  dZe + 1/2d0*( dmw2/mw2 + dZw + dZaa + dZp ) &
                  )                                           &
                - el*mw1*st/ct/2d0*dZza
      case (3)
        co(1) = el3/pi2/(32*st2*mw1)*(              &
                + mw2                               & ! bosonic
                + Nc*(                              & ! hadronic
                  - Qd*( mu2(1) + mu2(2) + mu2(3) ) & ! hadronic
                  + Qu*( md2(1) + md2(2) + md2(3) ) & ! hadronic
                  )                                 & ! hadronic
                )
      end select
    end select

  case ( 'p+ Z W-','p+ W- Z','W- p+ Z','W- Z p+','Z p+ W-','Z W- p+' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*mw1*st/ct
      case (2)
        co(1) = - el*mw1*st/ct*(                          &
                  + dZe + dst/(st*ct2)                    &
                  + 1/2d0*( dmw2/mw2 + dZw + dZzz + dZp ) &
                  )                                       &
                - el*mw1/2d0*dZaz
      case (3)
        co(1) = el3/pi2/(32*stct*mw1)*(             &
                + mw2                               & ! bosonic
                + Nc*(                              & ! hadronic
                  - Qd*( mu2(1) + mu2(2) + mu2(3) ) & ! hadronic
                  + Qu*( md2(1) + md2(2) + md2(3) ) & ! hadronic
                  )                                 & ! hadronic
                )
      end select
    end select

  case ( 'p- A W+','p- W+ A','A p- W+','A W+ p-','W+ A p-','W+ p- A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*mw1
      case (2)
        co(1) = - el*mw1*(                                    &
                  dZe + 1/2d0*( dmw2/mw2 + dZw + dZaa + dZp ) &
                  )                                           &
                - el*mw1*st/ct/2d0*dZza
      case (3)
        co(1) = el3/pi2/(32*st2*mw1)*(              &
                + mw2                               & ! bosonic
                + Nc*(                              & ! hadronic
                  - Qd*( mu2(1) + mu2(2) + mu2(3) ) & ! hadronic
                  + Qu*( md2(1) + md2(2) + md2(3) ) & ! hadronic
                  )                                 & ! hadronic
                )
      end select
    end select

  case ( 'p- Z W+','p- W+ Z','Z p- W+','Z W+ p-','W+ Z p-','W+ p- Z' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*mw1*st/ct
      case (2)
        co(1) = - el*mw1*st/ct*(                          &
                  + dZe + dst/(st*ct2)                    &
                  + 1/2d0*( dmw2/mw2 + dZw + dZzz + dZp ) &
                  )                                       &
                - el*mw1/2d0*dZaz
      case (3)
        co(1) = el3/pi2/(32*stct*mw1)*(             &
                + mw2                               & ! bosonic
                + Nc*(                              & ! hadronic
                  - Qd*( mu2(1) + mu2(2) + mu2(3) ) & ! hadronic
                  + Qu*( md2(1) + md2(2) + md2(3) ) & ! hadronic
                  )                                 & ! hadronic
                )
      end select
    end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 3 vectors
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! coupling = i*co(1)*(
!            + g_{mu_v1,mu_v2}*( p_v2 - p_v1 )_mu_v3
!            + g_{mu_v2,mu_v3}*( p_v3 - p_v2 )_mu_v1
!            + g_{mu_v3,mu_v1}*( p_v1 - p_v3 )_mu_v2
!            )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  p_v1 = momentum of the 1st vector
!  p_v2 = momentum of the 2nd vector
!  p_v3 = momentum of the 3rd vector
! mu_v1 = Lorentz index of the 1st vector
! mu_v2 = Lorentz index of the 2nd vector
! mu_v3 = Lorentz index of the 3rd vector
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! The vectors are ordered clockwise:
!
!        | v1
!        V
!        |
! --->---0---<---
! v3           v2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  case ( 'g g g' )

    select case (ngs)
    case (1)
      select case (xlp)
      case (0)
        co(1) = - gs/csq2
      end select
    case (3)
      select case (xlp)
      case (2)
        co(1) = - gs/csq2*( dZgs + 3/2d0*dZg )
      case (3)
        co(1) = gs3*Nc/(48*pi2)/csq2*( &
                + 7/4d0 + lam          & ! ggg+ggg+ggg
                + 2/Nc                 & ! gu{u}+gu{u}+gu{u}
                + 2/Nc                 & ! gc{c}+gc{c}+gc{c}
                + 2/Nc                 & ! gt{t}+gt{t}+gt{t}
                + 2/Nc                 & ! gd{d}+gd{d}+gd{d}
                + 2/Nc                 & ! gs{s}+gs{s}+gs{s}
                + 2/Nc                 & ! gb{b}+gb{b}+gb{b}
                )
      end select
    end select

 ! A g g is zero

 ! Z g g cancels when summing over quarks

  case ( 'W+ A W-','A W- W+','W- W+ A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el
      case (2)
        co(1) = el*(                          &
                + dZe + dZw                   &
                + 1/2d0*( dZaa - ct/st*dZza ) &
                )
      case (3)
        co(1) = - el3/pi2*(                     &
                  + 1/(96*st2)*( 7d0 + 4*lam )  & ! bosonic
                  + 3/(48*st2)                  & ! leptonic
                  + 3/(48*st2)*Nc               & ! hadronic
                  )
      end select
    end select

  case ( 'A W+ W-','W- A W+','W+ W- A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el
      case (2)
        co(1) = - el*(                          &
                  + dZe + dZw                   &
                  + 1/2d0*( dZaa - ct/st*dZza ) &
                  )
      case (3)
        co(1) = el3/pi2*(                     &
                + 1/(96*st2)*( 7d0 + 4*lam )  & ! bosonic
                + 3/(48*st2)                  & ! leptonic
                + 3/(48*st2)*Nc               & ! hadronic
                )
      end select
    end select

  case ( 'W+ Z W-','Z W- W+','W- W+ Z' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*ct/st
      case (2)
        co(1) = - el*ct/st*( dZe - dst/(st*ct2) + dZw + 1/2d0*dZzz ) &
                + el/2d0*dZaz
      case (3)
        co(1) = el3/pi2*ct/st*(              &
                + 1/(96*st2)*( 7d0 + 4*lam ) & ! bosonic
                + 3/(48*st2)                 & ! leptonic
                + 3/(48*st2)*Nc              & ! hadronic
                )
      end select
    end select

  case ( 'Z W+ W-','W- Z W+','W+ W- Z' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*ct/st
      case (2)
        co(1) = + el*ct/st*( dZe - dst/(st*ct2) + dZw + 1/2d0*dZzz ) &
                - el/2d0*dZaz
      case (3)
        co(1) = - el3/pi2*ct/st*(              &
                  + 1/(96*st2)*( 7d0 + 4*lam ) & ! bosonic
                  + 3/(48*st2)                 & ! leptonic
                  + 3/(48*st2)*Nc              & ! hadronic
                  )
      end select
    end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 1 scalar and 2 fermions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! coupling = i*( co(1)*omega_- + co(2)*omega_+ )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! omega_- = ( 1 - gamma5 )/2
! omega_+ = ( 1 + gamma5 )/2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  case ( 'H u u~','H u~ u','u H u~','u u~ H','u~ H u','u~ u H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*mu1(1)/(2*st*mw1)
      case (2)
        co(1) = - el*mu1(1)/(2*st*mw1)*(               &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZh + dZuRc(1) + dZuL(1) ) &
                  )                                    &
                - el/(2*st*mw1)*dmu(1)
        co(2) = - el*mu1(1)/(2*st*mw1)*(               &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZh + dZuR(1) + dZuLc(1) ) &
                  )                                    &
                - el/(2*st*mw1)*dmu(1)
      case (3)
        cP(1:2) = + el3/pi2*mu1(1)/(16*mw1*st)*Qu2*( 1d0 + lam )
        co(1:2) = + el3/pi2*mu1(1)/(16*mw1*st)*( &
                    + gpu*gmu*( 1d0 + lam )      &
                    + I3u2/(4*ct2*st2)           &
                    + 1/(8*st2)                  &
                    + md2(1)/(8*mw2*st2)         &
                    )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = - el*mu1(1)/(4*st*mw1)*( dZuRcqcd(1) + dZuLqcd(1) ) &
                - el/(2*st*mw1)*dmuqcd(1)
        co(2) = - el*mu1(1)/(4*st*mw1)*( dZuRqcd(1) + dZuLcqcd(1) ) &
                - el/(2*st*mw1)*dmuqcd(1)
      case (3)
        co(1:2) = el*gs2/(8*pi2)*mu1(1)/(2*st*mw1)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'H c c~','H c~ c','c H c~','c c~ H','c~ H c','c~ c H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*mu1(2)/(2*st*mw1)
      case (2)
        co(1) = - el*mu1(2)/(2*st*mw1)*(               &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZh + dZuRc(2) + dZuL(2) ) &
                  )                                    &
                - el/(2*st*mw1)*dmu(2)
        co(2) = - el*mu1(2)/(2*st*mw1)*(               &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZh + dZuR(2) + dZuLc(2) ) &
                  )                                    &
                - el/(2*st*mw1)*dmu(2)
      case (3)
        cP(1:2) = + el3/pi2*mu1(2)/(16*mw1*st)*Qu2*( 1d0 + lam )
        co(1:2) = + el3/pi2*mu1(2)/(16*mw1*st)*( &
                    + gpu*gmu*( 1d0 + lam )      &
                    + I3u2/(4*ct2*st2)           &
                    + 1/(8*st2)                  &
                    + md2(2)/(8*mw2*st2)         &
                    )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = - el*mu1(2)/(4*st*mw1)*( dZuRcqcd(2) + dZuLqcd(2) ) &
                - el/(2*st*mw1)*dmuqcd(2)
        co(2) = - el*mu1(2)/(4*st*mw1)*( dZuRqcd(2) + dZuLcqcd(2) ) &
                - el/(2*st*mw1)*dmuqcd(2)
      case (3)
        co(1:2) = el*gs2/(8*pi2)*mu1(2)/(2*st*mw1)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'H t t~','H t~ t','t H t~','t t~ H','t~ H t','t~ t H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*mu1(3)/(2*st*mw1)
      case (2)
        co(1) = - el*mu1(3)/(2*st*mw1)*(               &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZh + dZuRc(3) + dZuL(3) ) &
                  )                                    &
                - el/(2*st*mw1)*dmu(3)
        co(2) = - el*mu1(3)/(2*st*mw1)*(               &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZh + dZuR(3) + dZuLc(3) ) &
                  )                                    &
                - el/(2*st*mw1)*dmu(3)
      case (3)
        cP(1:2) = + el3/pi2*mu1(3)/(16*mw1*st)*Qu2*( 1d0 + lam )
        co(1:2) = + el3/pi2*mu1(3)/(16*mw1*st)*( &
                    + gpu*gmu*( 1d0 + lam )      &
                    + I3u2/(4*ct2*st2)           &
                    + 1/(8*st2)                  &
                    + md2(3)/(8*mw2*st2)         &
                    )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = - el*mu1(3)/(4*st*mw1)*( dZuRcqcd(3) + dZuLqcd(3) ) &
                - el/(2*st*mw1)*dmuqcd(3)
        co(2) = - el*mu1(3)/(4*st*mw1)*( dZuRqcd(3) + dZuLcqcd(3) ) &
                - el/(2*st*mw1)*dmuqcd(3)
      case (3)
        co(1:2) = el*gs2/(8*pi2)*mu1(3)/(2*st*mw1)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'H e- e+','H e+ e-','e- H e+','e- e+ H','e+ H e-','e+ e- H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*ml1(1)/(2*st*mw1)
      case (2)
        co(1)   = - el*ml1(1)/(2*st*mw1)*(               &
                    + dZe - dst/st - dmw/mw1             &
                    + 1/2d0*( dZh + dZlRc(1) + dZlL(1) ) &
                    )                                    &
                  - el/(2*st*mw1)*dml(1)
        co(2)   = - el*ml1(1)/(2*st*mw1)*(               &
                    + dZe - dst/st - dmw/mw1             &
                    + 1/2d0*( dZh + dZlR(1) + dZlLc(1) ) &
                    )                                    &
                  - el/(2*st*mw1)*dml(1)
      case (3)
        cP(1:2) = + el3/pi2*ml1(1)/(16*mw1*st)*Ql2*( 1d0 + lam )
        co(1:2) = + el3/pi2*ml1(1)/(16*mw1*st)*( &
                    + gpl*gml*( 1d0 + lam )      &
                    + I3l2/(4*ct2*st2)           &
                    + 1/(8*st2)                  &
                    )
      end select
    end select

  case ( 'H mu- mu+','H mu+ mu-','mu- H mu+', &
         'mu- mu+ H','mu+ H mu-','mu+ mu- H'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*ml1(2)/(2*st*mw1)
      case (2)
        co(1)   = - el*ml1(2)/(2*st*mw1)*(               &
                    + dZe - dst/st - dmw/mw1             &
                    + 1/2d0*( dZh + dZlRc(2) + dZlL(2) ) &
                    )                                    &
                  - el/(2*st*mw1)*dml(2)
        co(2)   = - el*ml1(2)/(2*st*mw1)*(               &
                    + dZe - dst/st - dmw/mw1             &
                    + 1/2d0*( dZh + dZlR(2) + dZlLc(2) ) &
                    )                                    &
                  - el/(2*st*mw1)*dml(2)
      case (3)
        cP(1:2) = + el3/pi2*ml1(2)/(16*mw1*st)*Ql2*( 1d0 + lam )
        co(1:2) = + el3/pi2*ml1(2)/(16*mw1*st)*( &
                    + gpl*gml*( 1d0 + lam )      &
                    + I3l2/(4*ct2*st2)           &
                    + 1/(8*st2)                  &
                    )
      end select
    end select

  case ( 'H tau- tau+','H tau+ tau-','tau- H tau+', &
         'tau- tau+ H','tau+ H tau-','tau+ tau- H'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*ml1(3)/(2*st*mw1)
      case (2)
        co(1)   = - el*ml1(3)/(2*st*mw1)*(               &
                    + dZe - dst/st - dmw/mw1             &
                    + 1/2d0*( dZh + dZlRc(3) + dZlL(3) ) &
                    )                                    &
                  - el/(2*st*mw1)*dml(3)
        co(2)   = - el*ml1(3)/(2*st*mw1)*(               &
                    + dZe - dst/st - dmw/mw1             &
                    + 1/2d0*( dZh + dZlR(3) + dZlLc(3) ) &
                    )                                    &
                  - el/(2*st*mw1)*dml(3)
      case (3)
        cP(1:2) = + el3/pi2*ml1(3)/(16*mw1*st)*Ql2*( 1d0 + lam )
        co(1:2) = + el3/pi2*ml1(3)/(16*mw1*st)*( &
                    + gpl*gml*( 1d0 + lam )      &
                    + I3l2/(4*ct2*st2)           &
                    + 1/(8*st2)                  &
                    )
      end select
    end select

  case ( 'H d d~','H d~ d','d H d~','d d~ H','d~ H d','d~ d H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*md1(1)/(2*st*mw1)
      case (2)
        co(1) = - el*md1(1)/(2*st*mw1)*(               &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZh + dZdRc(1) + dZdL(1) ) &
                  )                                    &
                - el/(2*st*mw1)*dmd(1)
        co(2) = - el*md1(1)/(2*st*mw1)*(               &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZh + dZdR(1) + dZdLc(1) ) &
                  )                                    &
                - el/(2*st*mw1)*dmd(1)
      case (3)
        cP(1:2) = + el3/pi2*md1(1)/(16*mw1*st)*Qd2*( 1d0 + lam )
        co(1:2) = + el3/pi2*md1(1)/(16*mw1*st)*( &
                    + gpd*gmd*( 1d0 + lam )      &
                    + I3d2/(4*ct2*st2)           &
                    + 1/(8*st2)                  &
                    + mu2(1)/(8*mw2*st2)         &
                    )
     end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = - el*md1(1)/(4*st*mw1)*( dZdRcqcd(1) + dZdLqcd(1) ) &
                - el/(2*st*mw1)*dmdqcd(1)
        co(2) = - el*md1(1)/(4*st*mw1)*( dZdRqcd(1) + dZdLcqcd(1) ) &
                - el/(2*st*mw1)*dmdqcd(1)
      case (3)
        co(1:2) = el*gs2/(8*pi2)*md1(1)/(2*st*mw1)*Cf*( 1d0 + lam )
     end select
   end select

  case ( 'H s s~','H s~ s','s H s~','s s~ H','s~ H s','s~ s H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*md1(2)/(2*st*mw1)
      case (2)
        co(1) = - el*md1(2)/(2*st*mw1)*(               &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZh + dZdRc(2) + dZdL(2) ) &
                  )                                    &
                - el/(2*st*mw1)*dmd(2)
        co(2) = - el*md1(2)/(2*st*mw1)*(               &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZh + dZdR(2) + dZdLc(2) ) &
                  )                                    &
                - el/(2*st*mw1)*dmd(2)
      case (3)
        cP(1:2) = + el3/pi2*md1(2)/(16*mw1*st)*Qd2*( 1d0 + lam )
        co(1:2) = + el3/pi2*md1(2)/(16*mw1*st)*( &
                    + gpd*gmd*( 1d0 + lam )      &
                    + I3d2/(4*ct2*st2)           &
                    + 1/(8*st2)                  &
                    + mu2(2)/(8*mw2*st2)         &
                    )
     end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = - el*md1(2)/(4*st*mw1)*( dZdRcqcd(2) + dZdLqcd(2) ) &
                - el/(2*st*mw1)*dmdqcd(2)
        co(2) = - el*md1(2)/(4*st*mw1)*( dZdRqcd(2) + dZdLcqcd(2) ) &
                - el/(2*st*mw1)*dmdqcd(2)
      case (3)
        co(1:2) = el*gs2/(8*pi2)*md1(2)/(2*st*mw1)*Cf*( 1d0 + lam )
     end select
   end select

  case ( 'H b b~','H b~ b','b H b~','b b~ H','b~ H b','b~ b H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*md1(3)/(2*st*mw1)
      case (2)
        co(1) = - el*md1(3)/(2*st*mw1)*(               &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZh + dZdRc(3) + dZdL(3) ) &
                  )                                    &
                - el/(2*st*mw1)*dmd(3)
        co(2) = - el*md1(3)/(2*st*mw1)*(               &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZh + dZdR(3) + dZdLc(3) ) &
                  )                                    &
                - el/(2*st*mw1)*dmd(3)
      case (3)
        cP(1:2) = + el3/pi2*md1(3)/(16*mw1*st)*Qd2*( 1d0 + lam )
        co(1:2) = + el3/pi2*md1(3)/(16*mw1*st)*( &
                    + gpd*gmd*( 1d0 + lam )      &
                    + I3d2/(4*ct2*st2)           &
                    + 1/(8*st2)                  &
                    + mu2(3)/(8*mw2*st2)         &
                    )
     end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = - el*md1(3)/(4*st*mw1)*( dZdRcqcd(3) + dZdLqcd(3) ) &
                - el/(2*st*mw1)*dmdqcd(3)
        co(2) = - el*md1(3)/(4*st*mw1)*( dZdRqcd(3) + dZdLcqcd(3) ) &
                - el/(2*st*mw1)*dmdqcd(3)
      case (3)
        co(1:2) = el*gs2/(8*pi2)*md1(3)/(2*st*mw1)*Cf*( 1d0 + lam )
     end select
   end select

  case ( 'p0 u u~','p0 u~ u','u p0 u~','u u~ p0','u~ p0 u','u~ u p0' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - cId0*el*mu1(1)/(st*mw1)*I3u
        co(2) = - co(1)
      case (2)
        co(1) = - cId0*el*mu1(1)/(st*mw1)*I3u*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZuRc(1) + dZuL(1) + dZp0 ) &
                  )                                     &
                - cId0*el/(st*mw1)*I3u*dmu(1)
        co(2) = + cId0*el*mu1(1)/(st*mw1)*I3u*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZuR(1) + dZuLc(1) + dZp0 ) &
                  )                                     &
                + cId0*el/(st*mw1)*I3u*dmu(1)
      case (3)
        cP(1) = cId0*el3/pi2*mu1(1)/(8*mw1*st)*I3u*Qu2*( 1d0 + lam )
        co(1) = cId0*el3/pi2*mu1(1)/(8*mw1*st)*(          &
                + I3u*gpu*gmu*( 1d0 + lam )               &
                + I3u/(16*ct2*st2)                        &
                + 1/(16*st2)                              &
                - I3u*mu2(1)/(4*mw2*st2)*( I3u2 - 1/4d0 ) &
                - I3d*md2(1)/(8*mw2*st2)                  &
                )
        cP(2) = - cP(1)
        co(2) = - co(1)
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = - cId0*el*mu1(1)/(2*st*mw1)*I3u* &
                  ( dZuRcqcd(1) + dZuLqcd(1) )   &
                - cId0*el/(st*mw1)*I3u*dmuqcd(1)
        co(2) = + cId0*el*mu1(1)/(2*st*mw1)*I3u* &
                  ( dZuRqcd(1) + dZuLcqcd(1) )   &
                + cId0*el/(st*mw1)*I3u*dmuqcd(1)
      case (3)
        co(1) = cId0*el*gs2/pi2*mu1(1)/(8*mw1*st)*I3u*Cf*( 1d0 + lam )
        co(2) = - co(1)
      end select
    end select

  case ( 'p0 c c~','p0 c~ c','c p0 c~','c c~ p0','c~ p0 c','c~ c p0' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - cId0*el*mu1(2)/(st*mw1)*I3u
        co(2) = - co(1)
      case (2)
        co(1) = - cId0*el*mu1(2)/(st*mw1)*I3u*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZuRc(2) + dZuL(2) + dZp0 ) &
                  )                                     &
                - cId0*el/(st*mw1)*I3u*dmu(2)
        co(2) = + cId0*el*mu1(2)/(st*mw1)*I3u*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZuR(2) + dZuLc(2) + dZp0 ) &
                  )                                     &
                + cId0*el/(st*mw1)*I3u*dmu(2)
      case (3)
        cP(1) = cId0*el3/pi2*mu1(2)/(8*mw1*st)*I3u*Qu2*( 1d0 + lam )
        co(1) = cId0*el3/pi2*mu1(2)/(8*mw1*st)*(          &
                + I3u*gpu*gmu*( 1d0 + lam )               &
                + I3u/(16*ct2*st2)                        &
                + 1/(16*st2)                              &
                - I3u*mu2(2)/(4*mw2*st2)*( I3u2 - 1/4d0 ) &
                - I3d*md2(2)/(8*mw2*st2)                  &
                )
        cP(2) = - cP(1)
        co(2) = - co(1)
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = - cId0*el*mu1(2)/(2*st*mw1)*I3u* &
                  ( dZuRcqcd(2) + dZuLqcd(2) )   &
                - cId0*el/(st*mw1)*I3u*dmuqcd(2)
        co(2) = + cId0*el*mu1(2)/(2*st*mw1)*I3u* &
                  ( dZuRqcd(2) + dZuLcqcd(2) )   &
                + cId0*el/(st*mw1)*I3u*dmuqcd(2)
      case (3)
        co(1) = cId0*el*gs2/pi2*mu1(2)/(8*mw1*st)*I3u*Cf*( 1d0 + lam )
        co(2) = - co(1)
      end select
    end select

  case ( 'p0 t t~','p0 t~ t','t p0 t~','t t~ p0','t~ p0 t','t~ t p0' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - cId0*el*mu1(3)/(st*mw1)*I3u
        co(2) = - co(1)
      case (2)
        co(1) = - cId0*el*mu1(3)/(st*mw1)*I3u*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZuRc(3) + dZuL(3) + dZp0 ) &
                  )                                     &
                - cId0*el/(st*mw1)*I3u*dmu(3)
        co(2) = + cId0*el*mu1(3)/(st*mw1)*I3u*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZuR(3) + dZuLc(3) + dZp0 ) &
                  )                                     &
                + cId0*el/(st*mw1)*I3u*dmu(3)
      case (3)
        cP(1) = cId0*el3/pi2*mu1(3)/(8*mw1*st)*I3u*Qu2*( 1d0 + lam )
        co(1) = cId0*el3/pi2*mu1(3)/(8*mw1*st)*(          &
                + I3u*gpu*gmu*( 1d0 + lam )               &
                + I3u/(16*ct2*st2)                        &
                + 1/(16*st2)                              &
                - I3u*mu2(3)/(4*mw2*st2)*( I3u2 - 1/4d0 ) &
                - I3d*md2(3)/(8*mw2*st2)                  &
                )
        cP(2) = - cP(1)
        co(2) = - co(1)
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = - cId0*el*mu1(3)/(2*st*mw1)*I3u* &
                  ( dZuRcqcd(3) + dZuLqcd(3) )   &
                - cId0*el/(st*mw1)*I3u*dmuqcd(3)
        co(2) = + cId0*el*mu1(3)/(2*st*mw1)*I3u* &
                  ( dZuRqcd(3) + dZuLcqcd(3) )   &
                + cId0*el/(st*mw1)*I3u*dmuqcd(3)
      case (3)
        co(1) = cId0*el*gs2/pi2*mu1(3)/(8*mw1*st)*I3u*Cf*( 1d0 + lam )
        co(2) = - co(1)
      end select
    end select

  case ( 'p0 e- e+','p0 e+ e-','e- p0 e+','e- e+ p0','e+ p0 e-','e+ e- p0' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - cId0*el*ml1(1)/(st*mw1)*I3l
        co(2) = - co(1)
      case (2)
        co(1) = - cId0*el*ml1(1)/(st*mw1)*I3l*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZlRc(1) + dZlL(1) + dZp0 ) &
                  )                                     &
                - cId0*el/(st*mw1)*I3l*dml(1)
        co(2) = + cId0*el*ml1(1)/(st*mw1)*I3l*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZlR(1) + dZlLc(1) + dZp0 ) &
                  )                                     &
                + cId0*el/(st*mw1)*I3l*dml(1)
      case (3)
        cP(1) = cId0*el3/pi2*ml1(1)/(8*mw1*st)*I3l*Ql2*( 1d0 + lam )
        co(1) = cId0*el3/pi2*ml1(1)/(8*mw1*st)*(          &
                + I3l*gpl*gml*( 1d0 + lam )               &
                + I3l/(16*ct2*st2)                        &
                - 1/(16*st2)                              &
                - I3l*ml2(1)/(4*mw2*st2)*( I3l2 - 1/4d0 ) &
                )
        cP(2) = - cP(1)
        co(2) = - co(1)
      end select
    end select

  case ( 'p0 mu- mu+','p0 mu+ mu-','mu- p0 mu+', &
         'mu- mu+ p0','mu+ p0 mu-','mu+ mu- p0'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - cId0*el*ml1(2)/(st*mw1)*I3l
        co(2) = - co(1)
      case (2)
        co(1) = - cId0*el*ml1(2)/(st*mw1)*I3l*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZlRc(2) + dZlL(2) + dZp0 ) &
                  )                                     &
                - cId0*el/(st*mw1)*I3l*dml(2)
        co(2) = + cId0*el*ml1(2)/(st*mw1)*I3l*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZlR(2) + dZlLc(2) + dZp0 ) &
                  )                                     &
                + cId0*el/(st*mw1)*I3l*dml(2)
      case (3)
        cP(1) = cId0*el3/pi2*ml1(2)/(8*mw1*st)*I3l*Ql2*( 1d0 + lam )
        co(1) = cId0*el3/pi2*ml1(2)/(8*mw1*st)*(          &
                + I3l*gpl*gml*( 1d0 + lam )               &
                + I3l/(16*ct2*st2)                        &
                - 1/(16*st2)                              &
                - I3l*ml2(2)/(4*mw2*st2)*( I3l2 - 1/4d0 ) &
                )
        cP(2) = - cP(1)
        co(2) = - co(1)
      end select
    end select

  case ( 'p0 tau- tau+','p0 tau+ tau-','tau- p0 tau+', &
         'tau- tau+ p0','tau+ p0 tau-','tau+ tau- p0'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - cId0*el*ml1(3)/(st*mw1)*I3l
        co(2) = - co(1)
      case (2)
        co(1) = - cId0*el*ml1(3)/(st*mw1)*I3l*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZlRc(3) + dZlL(3) + dZp0 ) &
                  )                                     &
                - cId0*el/(st*mw1)*I3l*dml(3)
        co(2) = + cId0*el*ml1(3)/(st*mw1)*I3l*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZlR(3) + dZlLc(3) + dZp0 ) &
                  )                                     &
                + cId0*el/(st*mw1)*I3l*dml(3)
      case (3)
        cP(1) = cId0*el3/pi2*ml1(3)/(8*mw1*st)*I3l*Ql2*( 1d0 + lam )
        co(1) = cId0*el3/pi2*ml1(3)/(8*mw1*st)*(          &
                + I3l*gpl*gml*( 1d0 + lam )               &
                + I3l/(16*ct2*st2)                        &
                - 1/(16*st2)                              &
                - I3l*ml2(3)/(4*mw2*st2)*( I3l2 - 1/4d0 ) &
                )
        cP(2) = - cP(1)
        co(2) = - co(1)
      end select
    end select

  case ( 'p0 d d~','p0 d~ d','d p0 d~','d d~ p0','d~ p0 d','d~ d p0' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - cId0*el*md1(1)/(st*mw1)*I3d
        co(2) = - co(1)
      case (2)
        co(1) = - cId0*el*md1(1)/(st*mw1)*I3d*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZdRc(1) + dZdL(1) + dZp0 ) &
                  )                                     &
                - cId0*el/(st*mw1)*I3d*dmd(1)
        co(2) = + cId0*el*md1(1)/(st*mw1)*I3d*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZdR(1) + dZdLc(1) + dZp0 ) &
                  )                                     &
                + cId0*el/(st*mw1)*I3d*dmd(1)
      case (3)
        cP(1) = cId0*el3/pi2*md1(1)/(8*mw1*st)*I3d*Qd2*( 1d0 + lam )
        co(1) = cId0*el3/pi2*md1(1)/(8*mw1*st)*(          &
                + I3d*gpd*gmd*( 1d0 + lam )               &
                + I3d/(16*ct2*st2)                        &
                - 1/(16*st2)                              &
                - I3d*md2(1)/(4*mw2*st2)*( I3d2 - 1/4d0 ) &
                - I3u*mu2(1)/(8*mw2*st2)                  &
                )
        cP(2) = - cP(1)
        co(2) = - co(1)
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = - cId0*el*md1(1)/(2*st*mw1)*I3d* &
                  ( dZdRcqcd(1) + dZdLqcd(1) )   &
                - cId0*el/(st*mw1)*I3d*dmdqcd(1)
        co(2) = + cId0*el*md1(1)/(2*st*mw1)*I3d* &
                  ( dZdRqcd(1) + dZdLcqcd(1) )   &
                + cId0*el/(st*mw1)*I3d*dmdqcd(1)
      case (3)
        co(1) = cId0*el*gs2/pi2*md1(1)/(8*mw1*st)*I3d*Cf*( 1d0 + lam )
        co(2) = - co(1)
      end select
    end select

  case ( 'p0 s s~','p0 s~ s','s p0 s~','s s~ p0','s~ p0 s','s~ s p0' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - cId0*el*md1(2)/(st*mw1)*I3d
        co(2) = - co(1)
      case (2)
        co(1) = - cId0*el*md1(2)/(st*mw1)*I3d*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZdRc(2) + dZdL(2) + dZp0 ) &
                  )                                     &
                - cId0*el/(st*mw1)*I3d*dmd(2)
        co(2) = + cId0*el*md1(2)/(st*mw1)*I3d*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZdR(2) + dZdLc(2) + dZp0 ) &
                  )                                     &
                + cId0*el/(st*mw1)*I3d*dmd(2)
      case (3)
        cP(1) = cId0*el3/pi2*md1(2)/(8*mw1*st)*I3d*Qd2*( 1d0 + lam )
        co(1) = cId0*el3/pi2*md1(2)/(8*mw1*st)*(          &
                + I3d*gpd*gmd*( 1d0 + lam )               &
                + I3d/(16*ct2*st2)                        &
                - 1/(16*st2)                              &
                - I3d*md2(2)/(4*mw2*st2)*( I3d2 - 1/4d0 ) &
                - I3u*mu2(2)/(8*mw2*st2)                  &
                )
        cP(2) = - cP(1)
        co(2) = - co(1)
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = - cId0*el*md1(2)/(2*st*mw1)*I3d* &
                  ( dZdRcqcd(2) + dZdLqcd(2) )   &
                - cId0*el/(st*mw1)*I3d*dmdqcd(2)
        co(2) = + cId0*el*md1(2)/(2*st*mw1)*I3d* &
                  ( dZdRqcd(2) + dZdLcqcd(2) )   &
                + cId0*el/(st*mw1)*I3d*dmdqcd(2)
      case (3)
        co(1) = cId0*el*gs2/pi2*md1(2)/(8*mw1*st)*I3d*Cf*( 1d0 + lam )
        co(2) = - co(1)
      end select
    end select

  case ( 'p0 b b~','p0 b~ b','b p0 b~','b b~ p0','b~ p0 b','b~ b p0' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - cId0*el*md1(3)/(st*mw1)*I3d
        co(2) = - co(1)
      case (2)
        co(1) = - cId0*el*md1(3)/(st*mw1)*I3d*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZdRc(3) + dZdL(3) + dZp0 ) &
                  )                                     &
                - cId0*el/(st*mw1)*I3d*dmd(3)
        co(2) = + cId0*el*md1(3)/(st*mw1)*I3d*(         &
                  + dZe - dst/st - dmw/mw1              &
                  + 1/2d0*( dZdR(3) + dZdLc(3) + dZp0 ) &
                  )                                     &
                + cId0*el/(st*mw1)*I3d*dmd(3)
      case (3)
        cP(1) = cId0*el3/pi2*md1(3)/(8*mw1*st)*I3d*Qd2*( 1d0 + lam )
        co(1) = cId0*el3/pi2*md1(3)/(8*mw1*st)*(          &
                + I3d*gpd*gmd*( 1d0 + lam )               &
                + I3d/(16*ct2*st2)                        &
                - 1/(16*st2)                              &
                - I3d*md2(3)/(4*mw2*st2)*( I3d2 - 1/4d0 ) &
                - I3u*mu2(3)/(8*mw2*st2)                  &
                )
        cP(2) = - cP(1)
        co(2) = - co(1)
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = - cId0*el*md1(3)/(2*st*mw1)*I3d* &
                  ( dZdRcqcd(3) + dZdLqcd(3) )   &
                - cId0*el/(st*mw1)*I3d*dmdqcd(3)
        co(2) = + cId0*el*md1(3)/(2*st*mw1)*I3d* &
                  ( dZdRqcd(3) + dZdLcqcd(3) )   &
                + cId0*el/(st*mw1)*I3d*dmdqcd(3)
      case (3)
        co(1) = cId0*el*gs2/pi2*md1(3)/(8*mw1*st)*I3d*Cf*( 1d0 + lam )
        co(2) = - co(1)
      end select
    end select

  case ( 'p+ e- nu_e~','p+ nu_e~ e-','e- p+ nu_e~', &
         'e- nu_e~ p+','nu_e~ p+ e-','nu_e~ e- p+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(2) = - el*ml1(1)/(csq2*st*mw1)
      case (2)
        co(2) = - el*ml1(1)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZlR(1) + dZnLc(1) + dZp ) &
                  )                                    &
                - el/(csq2*st*mw1)*dml(1)
      case (3)
        co(2) = + el3/pi2*ml1(1)/(8*csq2*st*mw1)*( &
                  + gpl*gmn*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 - I3n )        &
                  + 3/(16*st2)                     &
                  )
      end select
    end select

  case ( 'p+ mu- nu_mu~','p+ nu_mu~ mu-','mu- p+ nu_mu~', &
         'mu- nu_mu~ p+','nu_mu~ p+ mu-','nu_mu~ mu- p+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(2) = - el*ml1(2)/(csq2*st*mw1)
      case (2)
        co(2) = - el*ml1(2)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZlR(2) + dZnLc(2) + dZp ) &
                  )                                    &
                - el/(csq2*st*mw1)*dml(2)
      case (3)
        co(2) = + el3/pi2*ml1(2)/(8*csq2*st*mw1)*( &
                  + gpl*gmn*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 - I3n )        &
                  + 3/(16*st2)                     &
                  )
      end select
    end select

  case ( 'p+ tau- nu_tau~','p+ nu_tau~ tau-','tau- p+ nu_tau~', &
         'tau- nu_tau~ p+','nu_tau~ p+ tau-','nu_tau~ tau- p+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(2) = - el*ml1(3)/(csq2*st*mw1)
      case (2)
        co(2) = - el*ml1(3)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZlR(3) + dZnLc(3) + dZp ) &
                  )                                    &
                - el/(csq2*st*mw1)*dml(3)
      case (3)
        co(2) = + el3/pi2*ml1(3)/(8*csq2*st*mw1)*( &
                  + gpl*gmn*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 - I3n )        &
                  + 3/(16*st2)                     &
                  )
      end select
    end select

  case ( 'p+ d u~','p+ u~ d','d p+ u~','d u~ p+','u~ p+ d','u~ d p+' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = + el*mu1(1)/(csq2*st*mw1)
        co(2) = - el*md1(1)/(csq2*st*mw1)
      case (2)
        co(1) = + el*mu1(1)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZuRc(1) + dZdL(1) + dZp ) &
                  )                                    &
                + el/(csq2*st*mw1)*dmu(1)
        co(2) = - el*md1(1)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZdR(1) + dZuLc(1) + dZp ) &
                  )                                    &
                - el/(csq2*st*mw1)*dmd(1)
      case (3)
        cP(1) = - el3/pi2*mu1(1)/(8*csq2*st*mw1)*Qu*Qd*( 1d0 + lam )
        cP(2) = + el3/pi2*md1(1)/(8*csq2*st*mw1)*Qu*Qd*( 1d0 + lam )
        co(1) = - el3/pi2*mu1(1)/(8*csq2*st*mw1)*( &
                  + gpu*gmd*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 + I3d )        &
                  + 3/(16*st2)                     &
                  + md2(1)/(8*st2*mw2)             &
                  )
        co(2) = + el3/pi2*md1(1)/(8*csq2*st*mw1)*( &
                  + gmu*gpd*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 - I3u )        &
                  + 3/(16*st2)                     &
                  + mu2(1)/(8*st2*mw2)             &
                  )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = + el*mu1(1)/(2*csq2*st*mw1)*   &
                  ( dZuRcqcd(1) + dZdLqcd(1) ) &
                + el/(csq2*st*mw1)*dmuqcd(1)
        co(2) = - el*md1(1)/(2*csq2*st*mw1)*   &
                  ( dZdRqcd(1) + dZuLcqcd(1) ) &
                - el/(csq2*st*mw1)*dmdqcd(1)
      case (3)
        co(1) = - el*gs2/pi2*mu1(1)/(8*csq2*st*mw1)*Cf*( 1d0 + lam )
        co(2) = + el*gs2/pi2*md1(1)/(8*csq2*st*mw1)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'p+ s c~','p+ c~ s','s p+ c~','s c~ p+','c~ p+ s','c~ s p+' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = + el*mu1(2)/(csq2*st*mw1)
        co(2) = - el*md1(2)/(csq2*st*mw1)
      case (2)
        co(1) = + el*mu1(2)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZuRc(2) + dZdL(2) + dZp ) &
                  )                                    &
                + el/(csq2*st*mw1)*dmu(2)
        co(2) = - el*md1(2)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZdR(2) + dZuLc(2) + dZp ) &
                  )                                    &
                - el/(csq2*st*mw1)*dmd(2)
      case (3)
        cP(1) = - el3/pi2*mu1(2)/(8*csq2*st*mw1)*Qu*Qd*( 1d0 + lam )
        cP(2) = + el3/pi2*md1(2)/(8*csq2*st*mw1)*Qu*Qd*( 1d0 + lam )
        co(1) = - el3/pi2*mu1(2)/(8*csq2*st*mw1)*( &
                  + gpu*gmd*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 + I3d )        &
                  + 3/(16*st2)                     &
                  + md2(2)/(8*st2*mw2)             &
                  )
        co(2) = + el3/pi2*md1(2)/(8*csq2*st*mw1)*( &
                  + gmu*gpd*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 - I3u )        &
                  + 3/(16*st2)                     &
                  + mu2(2)/(8*st2*mw2)             &
                  )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = + el*mu1(2)/(2*csq2*st*mw1)*   &
                  ( dZuRcqcd(2) + dZdLqcd(2) ) &
                + el/(csq2*st*mw1)*dmuqcd(2)
        co(2) = - el*md1(2)/(2*csq2*st*mw1)*   &
                  ( dZdRqcd(2) + dZuLcqcd(2) ) &
                - el/(csq2*st*mw1)*dmdqcd(2)
      case (3)
        co(1) = - el*gs2/pi2*mu1(2)/(8*csq2*st*mw1)*Cf*( 1d0 + lam )
        co(2) = + el*gs2/pi2*md1(2)/(8*csq2*st*mw1)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'p+ b t~','p+ t~ b','b p+ t~','b t~ p+','t~ p+ b','t~ b p+' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = + el*mu1(3)/(csq2*st*mw1)
        co(2) = - el*md1(3)/(csq2*st*mw1)
      case (2)
        co(1) = + el*mu1(3)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZuRc(3) + dZdL(3) + dZp ) &
                  )                                    &
                + el/(csq2*st*mw1)*dmu(3)
        co(2) = - el*md1(3)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZdR(3) + dZuLc(3) + dZp ) &
                  )                                    &
                - el/(csq2*st*mw1)*dmd(3)
      case (3)
        cP(1) = - el3/pi2*mu1(3)/(8*csq2*st*mw1)*Qu*Qd*( 1d0 + lam )
        cP(2) = + el3/pi2*md1(3)/(8*csq2*st*mw1)*Qu*Qd*( 1d0 + lam )
        co(1) = - el3/pi2*mu1(3)/(8*csq2*st*mw1)*( &
                  + gpu*gmd*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 + I3d )        &
                  + 3/(16*st2)                     &
                  + md2(3)/(8*st2*mw2)             &
                  )
        co(2) = + el3/pi2*md1(3)/(8*csq2*st*mw1)*( &
                  + gmu*gpd*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 - I3u )        &
                  + 3/(16*st2)                     &
                  + mu2(3)/(8*st2*mw2)             &
                  )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = + el*mu1(3)/(2*csq2*st*mw1)*   &
                  ( dZuRcqcd(3) + dZdLqcd(3) ) &
                + el/(csq2*st*mw1)*dmuqcd(3)
        co(2) = - el*md1(3)/(2*csq2*st*mw1)*   &
                  ( dZdRqcd(3) + dZuLcqcd(3) ) &
                - el/(csq2*st*mw1)*dmdqcd(3)
      case (3)
        co(1) = - el*gs2/pi2*mu1(3)/(8*csq2*st*mw1)*Cf*( 1d0 + lam )
        co(2) = + el*gs2/pi2*md1(3)/(8*csq2*st*mw1)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'p- nu_e e+','p- e+ nu_e','e+ p- nu_e', &
         'e+ nu_e p-','nu_e p- e+','nu_e e+ p-'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*ml1(1)/(csq2*st*mw1)
      case (2)
        co(1) = - el*ml1(1)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZlRc(1) + dZnL(1) + dZp ) &
                  )                                    &
                - el/(csq2*st*mw1)*dml(1)
      case (3)
        co(1) = + el3/pi2*ml1(1)/(8*csq2*st*mw1)*( &
                  + gmn*gpl*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 - I3n )        &
                  + 3/(16*st2)                     &
                  )
      end select
    end select

  case ( 'p- nu_mu mu+','p- mu+ nu_mu','mu+ p- nu_mu', &
         'mu+ nu_mu p-','nu_mu p- mu+','nu_mu mu+ p-'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*ml1(2)/(csq2*st*mw1)
      case (2)
        co(1) = - el*ml1(2)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZlRc(2) + dZnL(2) + dZp ) &
                  )                                    &
                - el/(csq2*st*mw1)*dml(2)
      case (3)
        co(1) = + el3/pi2*ml1(2)/(8*csq2*st*mw1)*( &
                  + gmn*gpl*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 - I3n )        &
                  + 3/(16*st2)                     &
                  )
      end select
    end select

  case ( 'p- nu_tau tau+','p- tau+ nu_tau','tau+ p- nu_tau', &
         'tau+ nu_tau p-','nu_tau p- tau+','nu_tau tau+ p-'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*ml1(3)/(csq2*st*mw1)
      case (2)
        co(1) = - el*ml1(3)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZlRc(3) + dZnL(3) + dZp ) &
                  )                                    &
                - el/(csq2*st*mw1)*dml(3)
      case (3)
        co(1) = + el3/pi2*ml1(3)/(8*csq2*st*mw1)*( &
                  + gmn*gpl*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 - I3n )        &
                  + 3/(16*st2)                     &
                  )
      end select
    end select

  case ( 'p- d~ u','p- u d~','d~ p- u','d~ u p-','u p- d~','u d~ p-' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*md1(1)/(csq2*st*mw1)
        co(2) = + el*mu1(1)/(csq2*st*mw1)
      case (2)
        co(1) = - el*md1(1)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZdRc(1) + dZuL(1) + dZp ) &
                  )                                    &
                - el/(csq2*st*mw1)*dmd(1)
        co(2) = + el*mu1(1)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZuR(1) + dZdLc(1) + dZp ) &
                  )                                    &
                + el/(csq2*st*mw1)*dmu(1)
      case (3)
        cP(1) = + el3/pi2*md1(1)/(8*csq2*st*mw1)*Qu*Qd*( 1d0 + lam )
        cP(2) = - el3/pi2*mu1(1)/(8*csq2*st*mw1)*Qu*Qd*( 1d0 + lam )
        co(1) = + el3/pi2*md1(1)/(8*csq2*st*mw1)*( &
                  + gmu*gpd*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 - I3u )        &
                  + 3/(16*st2)                     &
                  + mu2(1)/(8*st2*mw2)             &
                  )
        co(2) = - el3/pi2*mu1(1)/(8*csq2*st*mw1)*( &
                  + gpu*gmd*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 + I3d )        &
                  + 3/(16*st2)                     &
                  + md2(1)/(8*st2*mw2)             &
                  )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = - el*md1(1)/(2*csq2*st*mw1)*   &
                  ( dZdRcqcd(1) + dZuLqcd(1) ) &
                - el/(csq2*st*mw1)*dmdqcd(1)
        co(2) = + el*mu1(1)/(2*csq2*st*mw1)*   &
                  ( dZuRqcd(1) + dZdLcqcd(1) ) &
                + el/(csq2*st*mw1)*dmuqcd(1)
      case (3)
        co(1) = + el*gs2/pi2*md1(1)/(8*csq2*st*mw1)*Cf*( 1d0 + lam )
        co(2) = - el*gs2/pi2*mu1(1)/(8*csq2*st*mw1)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'p- s~ c','p- c s~','s~ p- c','s~ c p-','c p- s~','c s~ p-' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*md1(2)/(csq2*st*mw1)
        co(2) = + el*mu1(2)/(csq2*st*mw1)
      case (2)
        co(1) = - el*md1(2)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZdRc(2) + dZuL(2) + dZp ) &
                  )                                    &
                - el/(csq2*st*mw1)*dmd(2)
        co(2) = + el*mu1(2)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZuR(2) + dZdLc(2) + dZp ) &
                  )                                    &
                + el/(csq2*st*mw1)*dmu(2)
      case (3)
        cP(1) = + el3/pi2*md1(2)/(8*csq2*st*mw1)*Qu*Qd*( 1d0 + lam )
        cP(2) = - el3/pi2*mu1(2)/(8*csq2*st*mw1)*Qu*Qd*( 1d0 + lam )
        co(1) = + el3/pi2*md1(2)/(8*csq2*st*mw1)*( &
                  + gmu*gpd*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 - I3u )        &
                  + 3/(16*st2)                     &
                  + mu2(2)/(8*st2*mw2)             &
                  )
        co(2) = - el3/pi2*mu1(2)/(8*csq2*st*mw1)*( &
                  + gpu*gmd*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 + I3d )        &
                  + 3/(16*st2)                     &
                  + md2(2)/(8*st2*mw2)             &
                  )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = - el*md1(2)/(2*csq2*st*mw1)*   &
                  ( dZdRcqcd(2) + dZuLqcd(2) ) &
                - el/(csq2*st*mw1)*dmdqcd(2)
        co(2) = + el*mu1(2)/(2*csq2*st*mw1)*   &
                  ( dZuRqcd(2) + dZdLcqcd(2) ) &
                + el/(csq2*st*mw1)*dmuqcd(2)
      case (3)
        co(1) = + el*gs2/pi2*md1(2)/(8*csq2*st*mw1)*Cf*( 1d0 + lam )
        co(2) = - el*gs2/pi2*mu1(2)/(8*csq2*st*mw1)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'p- b~ t','p- t b~','b~ p- t','b~ t p-','t p- b~','t b~ p-' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el*md1(3)/(csq2*st*mw1)
        co(2) = + el*mu1(3)/(csq2*st*mw1)
      case (2)
        co(1) = - el*md1(3)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZdRc(3) + dZuL(3) + dZp ) &
                  )                                    &
                - el/(csq2*st*mw1)*dmd(3)
        co(2) = + el*mu1(3)/(csq2*st*mw1)*(            &
                  + dZe - dst/st - dmw/mw1             &
                  + 1/2d0*( dZuR(3) + dZdLc(3) + dZp ) &
                  )                                    &
                + el/(csq2*st*mw1)*dmu(3)
      case (3)
        cP(1) = + el3/pi2*md1(3)/(8*csq2*st*mw1)*Qu*Qd*( 1d0 + lam )
        cP(2) = - el3/pi2*mu1(3)/(8*csq2*st*mw1)*Qu*Qd*( 1d0 + lam )
        co(1) = + el3/pi2*md1(3)/(8*csq2*st*mw1)*( &
                  + gmu*gpd*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 - I3u )        &
                  + 3/(16*st2)                     &
                  + mu2(3)/(8*st2*mw2)             &
                  )
        co(2) = - el3/pi2*mu1(3)/(8*csq2*st*mw1)*( &
                  + gpu*gmd*( 1d0 + lam )          &
                  + 1/(8*ct2)*( 1d0 + I3d )        &
                  + 3/(16*st2)                     &
                  + md2(3)/(8*st2*mw2)             &
                  )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = - el*md1(3)/(2*csq2*st*mw1)*   &
                  ( dZdRcqcd(3) + dZuLqcd(3) ) &
                - el/(csq2*st*mw1)*dmdqcd(3)
        co(2) = + el*mu1(3)/(2*csq2*st*mw1)*   &
                  ( dZuRqcd(3) + dZdLcqcd(3) ) &
                + el/(csq2*st*mw1)*dmuqcd(3)
      case (3)
        co(1) = + el*gs2/pi2*md1(3)/(8*csq2*st*mw1)*Cf*( 1d0 + lam )
        co(2) = - el*gs2/pi2*mu1(3)/(8*csq2*st*mw1)*Cf*( 1d0 + lam )
      end select
    end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 1 vector and 2 fermions
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! coupling = i*gamma^mu*( co(1)*omega_- + co(2)*omega_+ )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! omega_- = ( 1 - gamma5 )/2
! omega_+ = ( 1 + gamma5 )/2
! mu      = Lorentz index of the vector
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  case ( 'g u u~','g u~ u','u g u~','u u~ g','u~ g u','u~ u g' )

    select case (ngs)
    case (1)
      select case (xlp)
      case (0)
        co(1:2) = gs/csq2
      case (2)
        co(1)   = gs/csq2/2d0*( dZuL(1) + dZuLc(1) )
        co(2)   = gs/csq2/2d0*( dZuR(1) + dZuRc(1) )
      case (3)
        cP(1)   = - gs*el2/pi2/(16*csq2)*Qu2*( 1d0 + lam )
        cP(2)   = - gs*el2/pi2/(16*csq2)*Qu2*( 1d0 + lam )
        co(1)   = - gs*el2/pi2/(16*csq2)*(                      &
                    + ( gmu2 + 1/(2*st2) )*( 1d0 + lam )        &
                    + mu2(1)/(2*mw2*st2)*( 1/4d0 + I3u2 )       &
                    + md2(1)/(4*mw2*st2)                        &
                    )
        co(2)   = - gs*el2/pi2/(16*csq2)*(                        &
                    + gpu2*( 1d0 + lam )                          &
                    + mu2(1)/(2*st2*mw2)*( 1/2d0 + 1/4d0 + I3u2 ) &
                    )
      end select
    case (3)
      select case (xlp)
      case (2)
        co(1)   = gs/csq2*(                            &
                  + dZgs + 1/2d0*dZg                   &
                  + 1/2d0*( dZuLqcd(1) + dZuLcqcd(1) ) &
                  )
        co(2)   = gs/csq2*(                            &
                  + dZgs + 1/2d0*dZg                   &
                  + 1/2d0*( dZuRqcd(1) + dZuRcqcd(1) ) &
                  )
      case (3)
        co(1:2) = - gs3/pi2/(16*csq2)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'g c c~','g c~ c','c g c~','c c~ g','c~ g c','c~ c g' )

    select case (ngs)
    case (1)
      select case (xlp)
      case (0)
        co(1:2) = gs/csq2
      case (2)
        co(1)   = gs/csq2/2d0*( dZuL(2) + dZuLc(2) )
        co(2)   = gs/csq2/2d0*( dZuR(2) + dZuRc(2) )
      case (3)
        cP(1)   = - gs*el2/pi2/(16*csq2)*Qu2*( 1d0 + lam )
        cP(2)   = - gs*el2/pi2/(16*csq2)*Qu2*( 1d0 + lam )
        co(1)   = - gs*el2/pi2/(16*csq2)*(                      &
                    + ( gmu2 + 1/(2*st2) )*( 1d0 + lam )        &
                    + mu2(2)/(2*mw2*st2)*( 1/4d0 + I3u2 )       &
                    + md2(2)/(4*mw2*st2)                        &
                    )
        co(2)   = - gs*el2/pi2/(16*csq2)*(                        &
                    + gpu2*( 1d0 + lam )                          &
                    + mu2(2)/(2*st2*mw2)*( 1/2d0 + 1/4d0 + I3u2 ) &
                    )
      end select
    case (3)
      select case (xlp)
      case (2)
        co(1)   = gs/csq2*(                            &
                  + dZgs + 1/2d0*dZg                   &
                  + 1/2d0*( dZuLqcd(2) + dZuLcqcd(2) ) &
                  )
        co(2)   = gs/csq2*(                            &
                  + dZgs + 1/2d0*dZg                   &
                  + 1/2d0*( dZuRqcd(2) + dZuRcqcd(2) ) &
                  )
      case (3)
        co(1:2) = - gs3/pi2/(16*csq2)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'g t t~','g t~ t','t g t~','t t~ g','t~ g t','t~ t g' )

    select case (ngs)
    case (1)
      select case (xlp)
      case (0)
        co(1:2) = gs/csq2
      case (2)
        co(1)   = gs/csq2/2d0*( dZuL(3) + dZuLc(3) )
        co(2)   = gs/csq2/2d0*( dZuR(3) + dZuRc(3) )
      case (3)
        cP(1)   = - gs*el2/pi2/(16*csq2)*Qu2*( 1d0 + lam )
        cP(2)   = - gs*el2/pi2/(16*csq2)*Qu2*( 1d0 + lam )
        co(1)   = - gs*el2/pi2/(16*csq2)*(                      &
                    + ( gmu2 + 1/(2*st2) )*( 1d0 + lam )        &
                    + mu2(3)/(2*mw2*st2)*( 1/4d0 + I3u2 )       &
                    + md2(3)/(4*mw2*st2)                        &
                    )
        co(2)   = - gs*el2/pi2/(16*csq2)*(                        &
                    + gpu2*( 1d0 + lam )                          &
                    + mu2(3)/(2*st2*mw2)*( 1/2d0 + 1/4d0 + I3u2 ) &
                    )
      end select
    case (3)
      select case (xlp)
      case (2)
        co(1)   = gs/csq2*(                            &
                  + dZgs + 1/2d0*dZg                   &
                  + 1/2d0*( dZuLqcd(3) + dZuLcqcd(3) ) &
                  )
        co(2)   = gs/csq2*(                            &
                  + dZgs + 1/2d0*dZg                   &
                  + 1/2d0*( dZuRqcd(3) + dZuRcqcd(3) ) &
                  )
      case (3)
        co(1:2) = - gs3/pi2/(16*csq2)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'g d d~','g d~ d','d g d~','d d~ g','d~ g d','d~ d g' )

    select case (ngs)
    case (1)
      select case (xlp)
      case (0)
        co(1:2) = gs/csq2
      case (2)
        co(1)   = gs/csq2/2d0*( dZdL(1) + dZdLc(1) )
        co(2)   = gs/csq2/2d0*( dZdR(1) + dZdRc(1) )
      case (3)
        cP(1)   = - gs*el2/pi2/(16*csq2)*Qd2*( 1d0 + lam )
        cP(2)   = - gs*el2/pi2/(16*csq2)*Qd2*( 1d0 + lam )
        co(1)   = - gs*el2/pi2/(16*csq2)*(                      &
                    + ( gmd2 + 1/(2*st2) )*( 1d0 + lam )        &
                    + md2(1)/(2*mw2*st2)*( 1/4d0 + I3d2 )       &
                    + mu2(1)/(4*mw2*st2)                        &
                    )
        co(2)   = - gs*el2/pi2/(16*csq2)*(                        &
                    + gpd2*( 1d0 + lam )                          &
                    + md2(1)/(2*st2*mw2)*( 1/2d0 + 1/4d0 + I3d2 ) &
                    )
      end select
    case (3)
      select case (xlp)
      case (2)
        co(1)   = gs/csq2*(                            &
                  + dZgs + 1/2d0*dZg                   &
                  + 1/2d0*( dZdLqcd(1) + dZdLcqcd(1) ) &
                  )
        co(2)   = gs/csq2*(                            &
                  + dZgs + 1/2d0*dZg                   &
                  + 1/2d0*( dZdRqcd(1) + dZdRcqcd(1) ) &
                  )
      case (3)
        co(1:2) = - gs3/pi2/(16*csq2)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'g s s~','g s~ s','s g s~','s s~ g','s~ g s','s~ s g' )

    select case (ngs)
    case (1)
      select case (xlp)
      case (0)
        co(1:2) = gs/csq2
      case (2)
        co(1)   = gs/csq2/2d0*( dZdL(2) + dZdLc(2) )
        co(2)   = gs/csq2/2d0*( dZdR(2) + dZdRc(2) )
      case (3)
        cP(1)   = - gs*el2/pi2/(16*csq2)*Qd2*( 1d0 + lam )
        cP(2)   = - gs*el2/pi2/(16*csq2)*Qd2*( 1d0 + lam )
        co(1)   = - gs*el2/pi2/(16*csq2)*(                      &
                    + ( gmd2 + 1/(2*st2) )*( 1d0 + lam )        &
                    + md2(2)/(2*mw2*st2)*( 1/4d0 + I3d2 )       &
                    + mu2(2)/(4*mw2*st2)                        &
                    )
        co(2)   = - gs*el2/pi2/(16*csq2)*(                        &
                    + gpd2*( 1d0 + lam )                          &
                    + md2(2)/(2*st2*mw2)*( 1/2d0 + 1/4d0 + I3d2 ) &
                    )
      end select
    case (3)
      select case (xlp)
      case (2)
        co(1)   = gs/csq2*(                            &
                  + dZgs + 1/2d0*dZg                   &
                  + 1/2d0*( dZdLqcd(2) + dZdLcqcd(2) ) &
                  )
        co(2)   = gs/csq2*(                            &
                  + dZgs + 1/2d0*dZg                   &
                  + 1/2d0*( dZdRqcd(2) + dZdRcqcd(2) ) &
                  )
      case (3)
        co(1:2) = - gs3/pi2/(16*csq2)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'g b b~','g b~ b','b g b~','b b~ g','b~ g b','b~ b g' )

    select case (ngs)
    case (1)
      select case (xlp)
      case (0)
        co(1:2) = gs/csq2
      case (2)
        co(1)   = gs/csq2/2d0*( dZdL(3) + dZdLc(3) )
        co(2)   = gs/csq2/2d0*( dZdR(3) + dZdRc(3) )
      case (3)
        cP(1)   = - gs*el2/pi2/(16*csq2)*Qd2*( 1d0 + lam )
        cP(2)   = - gs*el2/pi2/(16*csq2)*Qd2*( 1d0 + lam )
        co(1)   = - gs*el2/pi2/(16*csq2)*(                      &
                    + ( gmd2 + 1/(2*st2) )*( 1d0 + lam )        &
                    + md2(3)/(2*mw2*st2)*( 1/4d0 + I3d2 )       &
                    + mu2(3)/(4*mw2*st2)                        &
                    )
        co(2)   = - gs*el2/pi2/(16*csq2)*(                        &
                    + gpd2*( 1d0 + lam )                          &
                    + md2(3)/(2*st2*mw2)*( 1/2d0 + 1/4d0 + I3d2 ) &
                    )
      end select
    case (3)
      select case (xlp)
      case (2)
        co(1)   = gs/csq2*(                            &
                  + dZgs + 1/2d0*dZg                   &
                  + 1/2d0*( dZdLqcd(3) + dZdLcqcd(3) ) &
                  )
        co(2)   = gs/csq2*(                            &
                  + dZgs + 1/2d0*dZg                   &
                  + 1/2d0*( dZdRqcd(3) + dZdRcqcd(3) ) &
                  )
      case (3)
        co(1:2) = - gs3/pi2/(16*csq2)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'A nu_e nu_e~','A nu_e~ nu_e','nu_e A nu_e~', &
         'nu_e nu_e~ A','nu_e~ A nu_e','nu_e~ nu_e A'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = el/2d0*gmn*dZza
        co(2) = el/2d0*gpn*dZza
      case (3)
        co(1) = el3/pi2/(32*st2)*(                  &
                + 1d0 + lam                         & ! aww+wnn+wnn
                + Ql*( 1d0 + lam + ml2(1)/(2*mw2) ) & ! all+wnl+wnl
                )
      end select
    end select

  case ( 'A nu_mu nu_mu~','A nu_mu~ nu_mu','nu_mu A nu_mu~', &
         'nu_mu nu_mu~ A','nu_mu~ A nu_mu','nu_mu~ nu_mu A'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = el/2d0*gmn*dZza
        co(2) = el/2d0*gpn*dZza
      case (3)
        co(1) = el3/pi2/(32*st2)*(                  &
                + 1d0 + lam                         & ! aww+wnn+wnn
                + Ql*( 1d0 + lam + ml2(2)/(2*mw2) ) & ! all+wnl+wnl
                )
      end select
    end select

  case ( 'A nu_tau nu_tau~','A nu_tau~ nu_tau','nu_tau A nu_tau~', &
         'nu_tau nu_tau~ A','nu_tau~ A nu_tau','nu_tau~ nu_tau A'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = el/2d0*gmn*dZza
        co(2) = el/2d0*gpn*dZza
      case (3)
        co(1) = el3/pi2/(32*st2)*(                  &
                + 1d0 + lam                         & ! aww+wnn+wnn
                + Ql*( 1d0 + lam + ml2(3)/(2*mw2) ) & ! all+wnl+wnl
                )
      end select
    end select

  case('A u u~','A u~ u','u A u~','u u~ A','u~ A u','u~ u A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*Qu
      case (2)
        co(1)   = - el*Qu*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZuL(1) + dZuLc(1) ) &
                    )                              &
                  + el/2d0*gmu*dZza
        co(2)   = - el*Qu*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZuR(1) + dZuRc(1) ) &
                    )                              &
                  + el/2d0*gpu*dZza
      case (3)
        cP(1)   = el3/(16*pi2)*Qu3*( 1d0 + lam )
        cP(2)   = el3/(16*pi2)*Qu3*( 1d0 + lam )
        co(1)   = el3/(16*pi2)*(                           &
                  + Qu*gmu2*( 1d0 + lam )                  &
                  + 1/(2*st2)*( Qd + 1d0 )*( 1d0 + lam )   &
                  + Qu*mu2(1)/(2*mw2*st2)*( 1/4d0 + I3u2 ) &
                  + Qd*md2(1)/(4*mw2*st2)                  &
                  )
        co(2)   = el3/(16*pi2)*(                                 &
                  + Qu*gpu2*( 1d0 + lam )                        &
                  + mu2(1)/(8*st2*mw2)*( 2*Qd + Qu + 4*Qu*I3u2 ) &
                  )
        cQED = cP
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1)   = - el*Qu/2d0*( dZuLqcd(1) + dZuLcqcd(1) )
        co(2)   = - el*Qu/2d0*( dZuRqcd(1) + dZuRcqcd(1) )
      case (3)
        co(1:2) = el*gs2/(16*pi2)*Qu*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'A c c~','A c~ c','c A c~','c c~ A','c~ A c','c~ c A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*Qu
      case (2)
        co(1)   = - el*Qu*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZuL(2) + dZuLc(2) ) &
                    )                              &
                  + el/2d0*gmu*dZza
        co(2)   = - el*Qu*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZuR(2) + dZuRc(2) ) &
                    )                              &
                  + el/2d0*gpu*dZza
      case (3)
        cP(1)   = el3/(16*pi2)*Qu3*( 1d0 + lam )
        cP(2)   = el3/(16*pi2)*Qu3*( 1d0 + lam )
        co(1)   = el3/(16*pi2)*(                           &
                  + Qu*gmu2*( 1d0 + lam )                  &
                  + 1/(2*st2)*( Qd + 1d0 )*( 1d0 + lam )   &
                  + Qu*mu2(2)/(2*mw2*st2)*( 1/4d0 + I3u2 ) &
                  + Qd*md2(2)/(4*mw2*st2)                  &
                  )
        co(2)   = el3/(16*pi2)*(                                 &
                  + Qu*gpu2*( 1d0 + lam )                        &
                  + mu2(2)/(8*st2*mw2)*( 2*Qd + Qu + 4*Qu*I3u2 ) &
                  )
        cQED = cP
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1)   = - el*Qu/2d0*( dZuLqcd(2) + dZuLcqcd(2) )
        co(2)   = - el*Qu/2d0*( dZuRqcd(2) + dZuRcqcd(2) )
      case (3)
        co(1:2) = el*gs2/(16*pi2)*Qu*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'A t t~','A t~ t','t A t~','t t~ A','t~ A t','t~ t A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*Qu
      case (2)
        co(1)   = - el*Qu*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZuL(3) + dZuLc(3) ) &
                    )                              &
                  + el/2d0*gmu*dZza
        co(2)   = - el*Qu*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZuR(3) + dZuRc(3) ) &
                    )                              &
                  + el/2d0*gpu*dZza
      case (3)
        cP(1)   = el3/(16*pi2)*Qu3*( 1d0 + lam )
        cP(2)   = el3/(16*pi2)*Qu3*( 1d0 + lam )
        co(1)   = el3/(16*pi2)*(                           &
                  + Qu*gmu2*( 1d0 + lam )                  &
                  + 1/(2*st2)*( Qd + 1d0 )*( 1d0 + lam )   &
                  + Qu*mu2(3)/(2*mw2*st2)*( 1/4d0 + I3u2 ) &
                  + Qd*md2(3)/(4*mw2*st2)                  &
                  )
        co(2)   = el3/(16*pi2)*(                                 &
                  + Qu*gpu2*( 1d0 + lam )                        &
                  + mu2(3)/(8*st2*mw2)*( 2*Qd + Qu + 4*Qu*I3u2 ) &
                  )
        cQED = cP
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1)   = - el*Qu/2d0*( dZuLqcd(3) + dZuLcqcd(3) )
        co(2)   = - el*Qu/2d0*( dZuRqcd(3) + dZuRcqcd(3) )
      case (3)
        co(1:2) = el*gs2/(16*pi2)*Qu*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'A e- e+','A e+ e-','e- A e+','e- e+ A','e+ A e-','e+ e- A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*Ql
      case (2)
        co(1)   = - el*Ql*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZlL(1) + dZlLc(1) ) &
                    )                              &
                  + el/2d0*gml*dZza
        co(2)   = - el*Ql*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZlR(1) + dZlRc(1) ) &
                    )                              &
                  + el/2d0*gpl*dZza
      case (3)
        cP(1)   = el3/(16*pi2)*Ql3*( 1d0 + lam )
        cP(2)   = el3/(16*pi2)*Ql3*( 1d0 + lam )
        co(1)   = el3/(16*pi2)*(                           &
                  + Ql*gml2*( 1d0 + lam )                  &
                  - 1/(2*st2)*( 1d0 + lam )                &
                  + Ql*ml2(1)/(2*mw2*st2)*( 1/4d0 + I3l2 ) &
                  )
        co(2)   = el3/(16*pi2)*(                           &
                  + Ql*gpl2*( 1d0 + lam )                  &
                  + Ql*ml2(1)/(2*st2*mw2)*( 1/4d0 + I3l2 ) &
                  )
        cQED = cP
      end select
    end select

  case ( 'A mu- mu+','A mu+ mu-','mu- A mu+', &
         'mu- mu+ A','mu+ A mu-','mu+ mu- A'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*Ql
      case (2)
        co(1)   = - el*Ql*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZlL(2) + dZlLc(2) ) &
                    )                              &
                  + el/2d0*gml*dZza
        co(2)   = - el*Ql*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZlR(2) + dZlRc(2) ) &
                    )                              &
                  + el/2d0*gpl*dZza
      case (3)
        cP(1)   = el3/(16*pi2)*Ql3*( 1d0 + lam )
        cP(2)   = el3/(16*pi2)*Ql3*( 1d0 + lam )
        co(1)   = el3/(16*pi2)*(                           &
                  + Ql*gml2*( 1d0 + lam )                  &
                  - 1/(2*st2)*( 1d0 + lam )                &
                  + Ql*ml2(2)/(2*mw2*st2)*( 1/4d0 + I3l2 ) &
                  )
        co(2)   = el3/(16*pi2)*(                           &
                  + Ql*gpl2*( 1d0 + lam )                  &
                  + Ql*ml2(2)/(2*st2*mw2)*( 1/4d0 + I3l2 ) &
                  )
        cQED = cP
      end select
    end select

  case ( 'A tau- tau+','A tau+ tau-','tau- A tau+', &
         'tau- tau+ A','tau+ A tau-','tau+ tau- A'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*Ql
      case (2)
        co(1)   = - el*Ql*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZlL(3) + dZlLc(3) ) &
                    )                              &
                  + el/2d0*gml*dZza
        co(2)   = - el*Ql*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZlR(3) + dZlRc(3) ) &
                    )                              &
                  + el/2d0*gpl*dZza
      case (3)
        cP(1)   = el3/(16*pi2)*Ql3*( 1d0 + lam )
        cP(2)   = el3/(16*pi2)*Ql3*( 1d0 + lam )
        co(1)   = el3/(16*pi2)*(                           &
                  + Ql*gml2*( 1d0 + lam )                  &
                  - 1/(2*st2)*( 1d0 + lam )                &
                  + Ql*ml2(3)/(2*mw2*st2)*( 1/4d0 + I3l2 ) &
                  )
        co(2)   = el3/(16*pi2)*(                           &
                  + Ql*gpl2*( 1d0 + lam )                  &
                  + Ql*ml2(3)/(2*st2*mw2)*( 1/4d0 + I3l2 ) &
                  )
        cQED = cP
      end select
    end select

  case ( 'A d d~','A d~ d','d A d~','d d~ A','d~ A d','d~ d A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*Qd
      case (2)
        co(1)   = - el*Qd*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZdL(1) + dZdLc(1) ) &
                    )                              &
                  + el/2d0*gmd*dZza
        co(2)   = - el*Qd*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZdR(1) + dZdRc(1) ) &
                    )                              &
                  + el/2d0*gpd*dZza
      case (3)
        cP(1)   = el3/(16*pi2)*Qd3*( 1d0 + lam )
        cP(2)   = el3/(16*pi2)*Qd3*( 1d0 + lam )
        co(1)   = el3/(16*pi2)*(                           &
                  + Qd*gmd2*( 1d0 + lam )                  &
                  + 1/(2*st2)*( Qu - 1d0 )*( 1d0 + lam )   &
                  + Qd*md2(1)/(2*mw2*st2)*( 1/4d0 + I3d2 ) &
                  + Qu*mu2(1)/(4*mw2*st2)                  &
                  )
        co(2)   = el3/(16*pi2)*(                                 &
                  + Qd*gpd2*( 1d0 + lam )                        &
                  + md2(1)/(8*st2*mw2)*( 2*Qu + Qd + 4*Qd*I3d2 ) &
                  )
        cQED = cP
     end select
    case (2)
      select case (xlp)
      case (2)
        co(1)   = - el*Qd/2d0*( dZdLqcd(1) + dZdLcqcd(1) )
        co(2)   = - el*Qd/2d0*( dZdRqcd(1) + dZdRcqcd(1) )
      case (3)
        co(1:2) = el*gs2/(16*pi2)*Qd*Cf*( 1d0 + lam )
     end select
   end select

  case ( 'A s s~','A s~ s','s A s~','s s~ A','s~ A s','s~ s A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*Qd
      case (2)
        co(1)   = - el*Qd*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZdL(2) + dZdLc(2) ) &
                    )                              &
                  + el/2d0*gmd*dZza
        co(2)   = - el*Qd*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZdR(2) + dZdRc(2) ) &
                    )                              &
                  + el/2d0*gpd*dZza
      case (3)
        cP(1)   = el3/(16*pi2)*Qd3*( 1d0 + lam )
        cP(2)   = el3/(16*pi2)*Qd3*( 1d0 + lam )
        co(1)   = el3/(16*pi2)*(                           &
                  + Qd*gmd2*( 1d0 + lam )                  &
                  + 1/(2*st2)*( Qu - 1d0 )*( 1d0 + lam )   &
                  + Qd*md2(2)/(2*mw2*st2)*( 1/4d0 + I3d2 ) &
                  + Qu*mu2(2)/(4*mw2*st2)                  &
                  )
        co(2)   = el3/(16*pi2)*(                                 &
                  + Qd*gpd2*( 1d0 + lam )                        &
                  + md2(2)/(8*st2*mw2)*( 2*Qu + Qd + 4*Qd*I3d2 ) &
                  )
        cQED = cP
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1)   = - el*Qd/2d0*( dZdLqcd(2) + dZdLcqcd(2) )
        co(2)   = - el*Qd/2d0*( dZdRqcd(2) + dZdRcqcd(2) )
      case (3)
        co(1:2) = el*gs2/(16*pi2)*Qd*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'A b b~','A b~ b','b A b~','b b~ A','b~ A b','b~ b A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1:2) = - el*Qd
      case (2)
        co(1)   = - el*Qd*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZdL(3) + dZdLc(3) ) &
                    )                              &
                  + el/2d0*gmd*dZza
        co(2)   = - el*Qd*(                        &
                    + dZe + 1/2d0*dZaa             &
                    + 1/2d0*( dZdR(3) + dZdRc(3) ) &
                    )                              &
                  + el/2d0*gpd*dZza
      case (3)
        cP(1)   = el3/(16*pi2)*Qd3*( 1d0 + lam )
        cP(2)   = el3/(16*pi2)*Qd3*( 1d0 + lam )
        co(1)   = el3/(16*pi2)*(                           &
                  + Qd*gmd2*( 1d0 + lam )                  &
                  + 1/(2*st2)*( Qu - 1d0 )*( 1d0 + lam )   &
                  + Qd*md2(3)/(2*mw2*st2)*( 1/4d0 + I3d2 ) &
                  + Qu*mu2(3)/(4*mw2*st2)                  &
                  )
        co(2)   = el3/(16*pi2)*(                                 &
                  + Qd*gpd2*( 1d0 + lam )                        &
                  + md2(3)/(8*st2*mw2)*( 2*Qu + Qd + 4*Qd*I3d2 ) &
                  )
        cQED = cP
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1)   = - el*Qd/2d0*( dZdLqcd(3) + dZdLcqcd(3) )
        co(2)   = - el*Qd/2d0*( dZdRqcd(3) + dZdRcqcd(3) )
      case (3)
        co(1:2) = el*gs2/(16*pi2)*Qd*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'Z nu_e nu_e~','Z nu_e~ nu_e','nu_e Z nu_e~', &
         'nu_e nu_e~ Z','nu_e~ Z nu_e','nu_e~ nu_e Z'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*gmn
      case (2)
        co(1) = + el*gmn/2d0*( dZzz + dZnL(1) + dZnLc(1) ) + el*dgmn
      case (3)
        co(1) = + el3/(16*pi2)*(                            &
                  - gmn*gmn2*( 1d0 + lam )                  &
                  - gml/(2*st2)*( 1d0 + lam )               &
                  - ct/(2*st3)*( 1d0 + lam )                &
                  + Ql*ml2(1)/(4*stct*mw2)                  &
                  )
      end select
    end select

  case ( 'Z nu_mu nu_mu~','Z nu_mu~ nu_mu','nu_mu Z nu_mu~', &
         'nu_mu nu_mu~ Z','nu_mu~ Z nu_mu','nu_mu~ nu_mu Z'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*gmn
      case (2)
        co(1) = + el*gmn/2d0*( dZzz + dZnL(2) + dZnLc(2) ) + el*dgmn
      case (3)
        co(1) = + el3/(16*pi2)*(                            &
                  - gmn*gmn2*( 1d0 + lam )                  &
                  - gml/(2*st2)*( 1d0 + lam )               &
                  - ct/(2*st3)*( 1d0 + lam )                &
                  + Ql*ml2(2)/(4*stct*mw2)                  &
                  )
      end select
    end select

  case ( 'Z nu_tau nu_tau~','Z nu_tau~ nu_tau','nu_tau Z nu_tau~', &
         'nu_tau nu_tau~ Z','nu_tau~ Z nu_tau','nu_tau~ nu_tau Z'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*gmn
      case (2)
        co(1) = + el*gmn/2d0*( dZzz + dZnL(3) + dZnLc(3) ) + el*dgmn
      case (3)
        co(1) = + el3/(16*pi2)*(                            &
                  - gmn*gmn2*( 1d0 + lam )                  &
                  - gml/(2*st2)*( 1d0 + lam )               &
                  - ct/(2*st3)*( 1d0 + lam )                &
                  + Ql*ml2(3)/(4*stct*mw2)                  &
                  )
      end select
    end select

  case ( 'Z u u~','Z u~ u','u Z u~','u u~ Z','u~ Z u','u~ u Z' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*gmu
        co(2) = el*gpu
      case (2)
        co(1) = + el*gmu/2d0*( dZzz + dZuL(1) + dZuLc(1) ) &
                + el*dgmu - el*Qu/2d0*dZaz
        co(2) = + el*gpu/2d0*( dZzz + dZuR(1) + dZuRc(1) ) &
                + el*dgpu - el*Qu/2d0*dZaz
      case (3)
        cP(1) = - el3/(16*pi2)*gmu*Qu2*( 1d0 + lam )
        cP(2) = - el3/(16*pi2)*gpu*Qu2*( 1d0 + lam )
        co(1) = + el3/(16*pi2)*(                            &
                  - gmu*gmu2*( 1d0 + lam )                  &
                  - gmd/(2*st2)*( 1d0 + lam )               &
                  - ct/(2*st3)*( 1d0 + lam )                &
                  + Qu*mu2(1)/(2*stct*mw2)*( 1/4d0 + I3u2 ) &
                  + Qd*md2(1)/(4*stct*mw2)                  &
                  )
        co(2) = + el3/(16*pi2)*(                     &
                  - gpu*gpu2*( 1d0 + lam )           &
                  - ( gmu + gmd )*mu2(1)/(4*st2*mw2) &
                  )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = el*gmu/2d0*( dZuLqcd(1) + dZuLcqcd(1) )
        co(2) = el*gpu/2d0*( dZuRqcd(1) + dZuRcqcd(1) )
      case (3)
        co(1) = - el*gs2/(16*pi2)*gmu*Cf*( 1d0 + lam )
        co(2) = - el*gs2/(16*pi2)*gpu*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'Z c c~','Z c~ c','c Z c~','c c~ Z','c~ Z c','c~ c Z' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*gmu
        co(2) = el*gpu
      case (2)
        co(1) = + el*gmu/2d0*( dZzz + dZuL(2) + dZuLc(2) ) &
                + el*dgmu - el*Qu/2d0*dZaz
        co(2) = + el*gpu/2d0*( dZzz + dZuR(2) + dZuRc(2) ) &
                + el*dgpu - el*Qu/2d0*dZaz
      case (3)
        cP(1) = - el3/(16*pi2)*gmu*Qu2*( 1d0 + lam )
        cP(2) = - el3/(16*pi2)*gpu*Qu2*( 1d0 + lam )
        co(1) = + el3/(16*pi2)*(                            &
                  - gmu*gmu2*( 1d0 + lam )                  &
                  - gmd/(2*st2)*( 1d0 + lam )               &
                  - ct/(2*st3)*( 1d0 + lam )                &
                  + Qu*mu2(2)/(2*stct*mw2)*( 1/4d0 + I3u2 ) &
                  + Qd*md2(2)/(4*stct*mw2)                  &
                  )
        co(2) = + el3/(16*pi2)*(                     &
                  - gpu*gpu2*( 1d0 + lam )           &
                  - ( gmu + gmd )*mu2(2)/(4*st2*mw2) &
                  )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = el*gmu/2d0*( dZuLqcd(2) + dZuLcqcd(2) )
        co(2) = el*gpu/2d0*( dZuRqcd(2) + dZuRcqcd(2) )
      case (3)
        co(1) = - el*gs2/(16*pi2)*gmu*Cf*( 1d0 + lam )
        co(2) = - el*gs2/(16*pi2)*gpu*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'Z t t~','Z t~ t','t Z t~','t t~ Z','t~ Z t','t~ t Z' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*gmu
        co(2) = el*gpu
      case (2)
        co(1) = + el*gmu/2d0*( dZzz + dZuL(3) + dZuLc(3) ) &
                + el*dgmu - el*Qu/2d0*dZaz
        co(2) = + el*gpu/2d0*( dZzz + dZuR(3) + dZuRc(3) ) &
                + el*dgpu - el*Qu/2d0*dZaz
      case (3)
        cP(1) = - el3/(16*pi2)*gmu*Qu2*( 1d0 + lam )
        cP(2) = - el3/(16*pi2)*gpu*Qu2*( 1d0 + lam )
        co(1) = + el3/(16*pi2)*(                            &
                  - gmu*gmu2*( 1d0 + lam )                  &
                  - gmd/(2*st2)*( 1d0 + lam )               &
                  - ct/(2*st3)*( 1d0 + lam )                &
                  + Qu*mu2(3)/(2*stct*mw2)*( 1/4d0 + I3u2 ) &
                  + Qd*md2(3)/(4*stct*mw2)                  &
                  )
        co(2) = + el3/(16*pi2)*(                     &
                  - gpu*gpu2*( 1d0 + lam )           &
                  - ( gmu + gmd )*mu2(3)/(4*st2*mw2) &
                  )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = el*gmu/2d0*( dZuLqcd(3) + dZuLcqcd(3) )
        co(2) = el*gpu/2d0*( dZuRqcd(3) + dZuRcqcd(3) )
      case (3)
        co(1) = - el*gs2/(16*pi2)*gmu*Cf*( 1d0 + lam )
        co(2) = - el*gs2/(16*pi2)*gpu*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'Z e- e+','Z e+ e-','e- Z e+','e- e+ Z','e+ Z e-','e+ e- Z' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*gml
        co(2) = el*gpl
      case (2)
        co(1) = + el*gml/2d0*( dZzz + dZlL(1) + dZlLc(1) ) &
                + el*dgml - el*Ql/2d0*dZaz
        co(2) = + el*gpl/2d0*( dZzz + dZlR(1) + dZlRc(1) ) &
                + el*dgpl - el*Ql/2d0*dZaz
      case (3)
        cP(1) = - el3/(16*pi2)*gml*Ql2*( 1d0 + lam )
        cP(2) = - el3/(16*pi2)*gpl*Ql2*( 1d0 + lam )
        co(1) = + el3/(16*pi2)*(                            &
                  - gml*gml2*( 1d0 + lam )                  &
                  - gmn/(2*st2)*( 1d0 + lam )               &
                  + ct/(2*st3)*( 1d0 + lam )                &
                  + Ql*ml2(1)/(2*stct*mw2)*( 1/4d0 + I3l2 ) &
                  )
        co(2) = + el3/(16*pi2)*(                     &
                  - gpl*gpl2*( 1d0 + lam )           &
                  - ( gmn + gml )*ml2(1)/(4*st2*mw2) &
                  )
      end select
    end select

  case ( 'Z mu- mu+','Z mu+ mu-','mu- Z mu+', &
         'mu- mu+ Z','mu+ Z mu-','mu+ mu- Z'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*gml
        co(2) = el*gpl
      case (2)
        co(1) = + el*gml/2d0*( dZzz + dZlL(2) + dZlLc(2) ) &
                + el*dgml - el*Ql/2d0*dZaz
        co(2) = + el*gpl/2d0*( dZzz + dZlR(2) + dZlRc(2) ) &
                + el*dgpl - el*Ql/2d0*dZaz
      case (3)
        cP(1) = - el3/(16*pi2)*gml*Ql2*( 1d0 + lam )
        cP(2) = - el3/(16*pi2)*gpl*Ql2*( 1d0 + lam )
        co(1) = + el3/(16*pi2)*(                            &
                  - gml*gml2*( 1d0 + lam )                  &
                  - gmn/(2*st2)*( 1d0 + lam )               &
                  + ct/(2*st3)*( 1d0 + lam )                &
                  + Ql*ml2(2)/(2*stct*mw2)*( 1/4d0 + I3l2 ) &
                  )
        co(2) = + el3/(16*pi2)*(                     &
                  - gpl*gpl2*( 1d0 + lam )           &
                  - ( gmn + gml )*ml2(2)/(4*st2*mw2) &
                  )
      end select
    end select

  case ( 'Z tau- tau+','Z tau+ tau-','tau- Z tau+', &
         'tau- tau+ Z','tau+ Z tau-','tau+ tau- Z'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*gml
        co(2) = el*gpl
      case (2)
        co(1) = + el*gml/2d0*( dZzz + dZlL(3) + dZlLc(3) ) &
                + el*dgml - el*Ql/2d0*dZaz
        co(2) = + el*gpl/2d0*( dZzz + dZlR(3) + dZlRc(3) ) &
                + el*dgpl - el*Ql/2d0*dZaz
      case (3)
        cP(1) = - el3/(16*pi2)*gml*Ql2*( 1d0 + lam )
        cP(2) = - el3/(16*pi2)*gpl*Ql2*( 1d0 + lam )
        co(1) = + el3/(16*pi2)*(                            &
                  - gml*gml2*( 1d0 + lam )                  &
                  - gmn/(2*st2)*( 1d0 + lam )               &
                  + ct/(2*st3)*( 1d0 + lam )                &
                  + Ql*ml2(3)/(2*stct*mw2)*( 1/4d0 + I3l2 ) &
                  )
        co(2) = + el3/(16*pi2)*(                     &
                  - gpl*gpl2*( 1d0 + lam )           &
                  - ( gmn + gml )*ml2(3)/(4*st2*mw2) &
                  )
      end select
    end select

  case ( 'Z d d~','Z d~ d','d Z d~','d d~ Z','d~ Z d','d~ d Z' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*gmd
        co(2) = el*gpd
      case (2)
        co(1) = + el*gmd/2d0*( dZzz + dZdL(1) + dZdLc(1) ) &
                + el*dgmd - el*Qd/2d0*dZaz
        co(2) = + el*gpd/2d0*( dZzz + dZdR(1) + dZdRc(1) ) &
                + el*dgpd - el*Qd/2d0*dZaz
      case (3)
        cP(1) = - el3/(16*pi2)*gmd*Qd2*( 1d0 + lam )
        cP(2) = - el3/(16*pi2)*gpd*Qd2*( 1d0 + lam )
        co(1) = + el3/(16*pi2)*(                            &
                  - gmd*gmd2*( 1d0 + lam )                  &
                  - gmu/(2*st2)*( 1d0 + lam )               &
                  + ct/(2*st3)*( 1d0 + lam )                &
                  + Qd*md2(1)/(2*stct*mw2)*( 1/4d0 + I3d2 ) &
                  + Qu*mu2(1)/(4*stct*mw2)                  &
                  )
        co(2) = + el3/(16*pi2)*(                     &
                  - gpd*gpd2*( 1d0 + lam )           &
                  - ( gmu + gmd )*md2(1)/(4*st2*mw2) &
                  )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = el*gmd/2d0*( dZdLqcd(1) + dZdLcqcd(1) )
        co(2) = el*gpd/2d0*( dZdRqcd(1) + dZdRcqcd(1) )
      case (3)
        co(1) = - el*gs2/(16*pi2)*gmd*Cf*( 1d0 + lam )
        co(2) = - el*gs2/(16*pi2)*gpd*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'Z s s~','Z s~ s','s Z s~','s s~ Z','s~ Z s','s~ s Z' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*gmd
        co(2) = el*gpd
      case (2)
        co(1) = + el*gmd/2d0*( dZzz + dZdL(2) + dZdLc(2) ) &
                + el*dgmd - el*Qd/2d0*dZaz
        co(2) = + el*gpd/2d0*( dZzz + dZdR(2) + dZdRc(2) ) &
                + el*dgpd - el*Qd/2d0*dZaz
      case (3)
        cP(1) = - el3/(16*pi2)*gmd*Qd2*( 1d0 + lam )
        cP(2) = - el3/(16*pi2)*gpd*Qd2*( 1d0 + lam )
        co(1) = + el3/(16*pi2)*(                            &
                  - gmd*gmd2*( 1d0 + lam )                  &
                  - gmu/(2*st2)*( 1d0 + lam )               &
                  + ct/(2*st3)*( 1d0 + lam )                &
                  + Qd*md2(2)/(2*stct*mw2)*( 1/4d0 + I3d2 ) &
                  + Qu*mu2(2)/(4*stct*mw2)                  &
                  )
        co(2) = + el3/(16*pi2)*(                     &
                  - gpd*gpd2*( 1d0 + lam )           &
                  - ( gmu + gmd )*md2(2)/(4*st2*mw2) &
                  )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = el*gmd/2d0*( dZdLqcd(2) + dZdLcqcd(2) )
        co(2) = el*gpd/2d0*( dZdRqcd(2) + dZdRcqcd(2) )
      case (3)
        co(1) = - el*gs2/(16*pi2)*gmd*Cf*( 1d0 + lam )
        co(2) = - el*gs2/(16*pi2)*gpd*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'Z b b~','Z b~ b','b Z b~','b b~ Z','b~ Z b','b~ b Z' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el*gmd
        co(2) = el*gpd
      case (2)
        co(1) = + el*gmd/2d0*( dZzz + dZdL(3) + dZdLc(3) ) &
                + el*dgmd - el*Qd/2d0*dZaz
        co(2) = + el*gpd/2d0*( dZzz + dZdR(3) + dZdRc(3) ) &
                + el*dgpd - el*Qd/2d0*dZaz
      case (3)
        cP(1) = - el3/(16*pi2)*gmd*Qd2*( 1d0 + lam )
        cP(2) = - el3/(16*pi2)*gpd*Qd2*( 1d0 + lam )
        co(1) = + el3/(16*pi2)*(                            &
                  - gmd*gmd2*( 1d0 + lam )                  &
                  - gmu/(2*st2)*( 1d0 + lam )               &
                  + ct/(2*st3)*( 1d0 + lam )                &
                  + Qd*md2(3)/(2*stct*mw2)*( 1/4d0 + I3d2 ) &
                  + Qu*mu2(3)/(4*stct*mw2)                  &
                  )
        co(2) = + el3/(16*pi2)*(                     &
                  - gpd*gpd2*( 1d0 + lam )           &
                  - ( gmu + gmd )*md2(3)/(4*st2*mw2) &
                  )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = el*gmd/2d0*( dZdLqcd(3) + dZdLcqcd(3) )
        co(2) = el*gpd/2d0*( dZdRqcd(3) + dZdRcqcd(3) )
      case (3)
        co(1) = - el*gs2/(16*pi2)*gmd*Cf*( 1d0 + lam )
        co(2) = - el*gs2/(16*pi2)*gpd*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'W+ e- nu_e~','W+ nu_e~ e-','e- W+ nu_e~', &
         'e- nu_e~ W+','nu_e~ W+ e-','nu_e~ e- W+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el/(csq2*st)
      case (2)
        co(1) = el/(csq2*st)*(                       &
                + dZe - dst/st                       &
                + 1/2d0*( dZw + dZnLc(1) + dZlL(1) ) &
                )
      case (3)
        co(1) = - el3/(16*pi2*csq2*st)*( 1d0 + lam )* &
                  ( gmn*gml + 1/st2 )
      end select
    end select

  case ( 'W+ mu- nu_mu~','W+ nu_mu~ mu-','mu- W+ nu_mu~', &
         'mu- nu_mu~ W+','nu_mu~ W+ mu-','nu_mu~ mu- W+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el/(csq2*st)
      case (2)
        co(1) = el/(csq2*st)*(                       &
                + dZe - dst/st                       &
                + 1/2d0*( dZw + dZnLc(2) + dZlL(2) ) &
                )
      case (3)
        co(1) = - el3/(16*pi2*csq2*st)*( 1d0 + lam )* &
                  ( gmn*gml + 1/st2 )
      end select
    end select

  case ( 'W+ tau- nu_tau~','W+ nu_tau~ tau-','tau- W+ nu_tau~', &
         'tau- nu_tau~ W+','nu_tau~ W+ tau-','nu_tau~ tau- W+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el/(csq2*st)
      case (2)
        co(1) = el/(csq2*st)*(                       &
                + dZe - dst/st                       &
                + 1/2d0*( dZw + dZnLc(3) + dZlL(3) ) &
                )
      case (3)
        co(1) = - el3/(16*pi2*csq2*st)*( 1d0 + lam )* &
                  ( gmn*gml + 1/st2 )
      end select
    end select

  case ( 'W+ d u~','W+ u~ d','d W+ u~','d u~ W+','u~ W+ d','u~ d W+' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el/(csq2*st)
      case (2)
        co(1) = el/(csq2*st)*(                       &
                + dZe - dst/st                       &
                + 1/2d0*( dZw + dZuLc(1) + dZdL(1) ) &
                )
      case (3)
        cP(1) = - el3/(16*pi2*csq2*st)*Qu*Qd*( 1d0 + lam )
        co(1) = - el3/(16*pi2*csq2*st)*( 1d0 + lam )* &
                  ( gmu*gmd + 1/st2 )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = + el/(2*csq2*st)*( dZuLcqcd(1) + dZdLqcd(1) )
      case (3)
        co(1) = - el*gs2/(16*pi2*csq2*st)*Cf*( 1d0 + lam )
      end select
    end select


  case ( 'W+ s c~','W+ c~ s','s W+ c~','s c~ W+','c~ W+ s','c~ s W+' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el/(csq2*st)
      case (2)
        co(1) = el/(csq2*st)*(                       &
                + dZe - dst/st                       &
                + 1/2d0*( dZw + dZuLc(2) + dZdL(2) ) &
                )
      case (3)
        cP(1) = - el3/(16*pi2*csq2*st)*Qu*Qd*( 1d0 + lam )
        co(1) = - el3/(16*pi2*csq2*st)*( 1d0 + lam )* &
                  ( gmu*gmd + 1/st2 )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = + el/(2*csq2*st)*( dZuLcqcd(2) + dZdLqcd(2) )
      case (3)
        co(1) = - el*gs2/(16*pi2*csq2*st)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'W+ b t~','W+ t~ b','b W+ t~','b t~ W+','t~ W+ b','t~ b W+' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el/(csq2*st)
      case (2)
        co(1) = el/(csq2*st)*(                       &
                + dZe - dst/st                       &
                + 1/2d0*( dZw + dZuLc(3) + dZdL(3) ) &
                )
      case (3)
        cP(1) = - el3/(16*pi2*csq2*st)*Qu*Qd*( 1d0 + lam )
        co(1) = - el3/(16*pi2*csq2*st)*( 1d0 + lam )* &
                  ( gmu*gmd + 1/st2 )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = + el/(2*csq2*st)*( dZuLcqcd(3) + dZdLqcd(3) )
      case (3)
        co(1) = - el*gs2/(16*pi2*csq2*st)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'W- nu_e e+','W- e+ nu_e','e+ W- nu_e', &
         'e+ nu_e W-','nu_e W- e+','nu_e e+ W-'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el/(csq2*st)
      case (2)
        co(1) = el/(csq2*st)*(                       &
                + dZe - dst/st                       &
                + 1/2d0*( dZw + dZlLc(1) + dZnL(1) ) &
                )
      case (3)
        co(1) = - el3/(16*pi2*csq2*st)*( 1d0 + lam )* &
                  ( gmn*gml + 1/st2 )
      end select
    end select

  case ( 'W- nu_mu mu+','W- mu+ nu_mu','mu+ W- nu_mu', &
         'mu+ nu_mu W-','nu_mu W- mu+','nu_mu mu+ W-'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el/(csq2*st)
      case (2)
        co(1) =  el/(csq2*st)*(                       &
                 + dZe - dst/st                       &
                 + 1/2d0*( dZw + dZlLc(2) + dZnL(2) ) &
                 )
      case (3)
        co(1) = - el3/(16*pi2*csq2*st)*( 1d0 + lam )* &
                  ( gmn*gml + 1/st2 )
      end select
    end select

  case ( 'W- nu_tau tau+','W- tau+ nu_tau','tau+ W- nu_tau', &
         'tau+ nu_tau W-','nu_tau W- tau+','nu_tau tau+ W-'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el/(csq2*st)
      case (2)
        co(1) = el/(csq2*st)*(                       &
                + dZe - dst/st                       &
                + 1/2d0*( dZw + dZlLc(3) + dZnL(3) ) &
                )
      case (3)
        co(1) = - el3/(16*pi2*csq2*st)*( 1d0 + lam )* &
                  ( gmn*gml + 1/st2 )
      end select
    end select

  case ( 'W- d~ u','W- u d~','d~ W- u','d~ u W-','u W- d~','u d~ W-' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el/(csq2*st)
      case (2)
        co(1) = el/(csq2*st)*(                       &
                + dZe - dst/st                       &
                + 1/2d0*( dZw + dZdLc(1) + dZuL(1) ) &
                )
      case (3)
        cP(1) = - el3/(16*pi2*csq2*st)*Qu*Qd*( 1d0 + lam )
        co(1) = - el3/(16*pi2*csq2*st)*( 1d0 + lam )* &
                  ( gmu*gmd + 1/st2 )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = el/(2*csq2*st)*( dZdLcqcd(1) + dZuLqcd(1) )
      case (3)
        co(1) = - el*gs2/(16*pi2*csq2*st)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'W- s~ c','W- c s~','s~ W- c','s~ c W-','c W- s~','c s~ W-' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el/(csq2*st)
      case (2)
        co(1) = el/(csq2*st)*(                       &
                + dZe - dst/st                       &
                + 1/2d0*( dZw + dZdLc(2) + dZuL(2) ) &
                )
      case (3)
        cP(1) = - el3/(16*pi2*csq2*st)*Qu*Qd*( 1d0 + lam )
        co(1) = - el3/(16*pi2*csq2*st)*( 1d0 + lam )* &
                  ( gmu*gmd + 1/st2 )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = el/(2*csq2*st)*( dZdLcqcd(2) + dZuLqcd(2) )
      case (3)
        co(1) = - el*gs2/(16*pi2*csq2*st)*Cf*( 1d0 + lam )
      end select
    end select

  case ( 'W- b~ t','W- t b~','b~ W- t','b~ t W-','t W- b~','t b~ W-' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el/(csq2*st)
      case (2)
        co(1) = el/(csq2*st)*(                       &
                + dZe - dst/st                       &
                + 1/2d0*( dZw + dZdLc(3) + dZuL(3) ) &
                )
      case (3)
        cP(1) = - el3/(16*pi2*csq2*st)*Qu*Qd*( 1d0 + lam )
        co(1) = - el3/(16*pi2*csq2*st)*( 1d0 + lam )* &
                  ( gmu*gmd + 1/st2 )
      end select
    case (2)
      select case (xlp)
      case (2)
        co(1) = el/(2*csq2*st)*( dZdLcqcd(3) + dZuLqcd(3) )
      case (3)
        co(1) = - el*gs2/(16*pi2*csq2*st)*Cf*( 1d0 + lam )
      end select
    end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 4 scalars
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! coupling = i*co(1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  case ( 'H H H H' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el2*mh2*3/(4*mw2*st2)
      case (2)
        co(1) = - el2*mh2*3/(4*mw2*st2)*(                   &
                  + 2*dZe - 2*dst/st + dmh2/mh2             &
                  + el*dt/(2*st*mw1*mh2) - dmw2/mw2 + 2*dZh &
                  )
      case (3)
        co(1) = el4/(64*pi2*st4)*(                           &
                + 1/2d0*( 1d0 + 1/(2*ct4) )*( 1d0 - 12*lam ) & ! bosonic
                + 3/2d0*mh2/mw2*( 1d0 + 1/(2*ct2) )          & ! bosonic
                + 5/mw4*( ml4(1) + ml4(2) + ml4(3) )         & ! leptonic
                + 5/mw4*Nc*(                                 & ! hadronic
                  + md4(1) + md4(2) + md4(3)                 & ! hadronic
                  + mu4(1) + mu4(2) + mu4(3)                 & ! hadronic
                  )                                          & ! hadronic
                )
      end select
    end select

  case ( 'H H p0 p0','H p0 H p0','H p0 p0 H', &
         'p0 p0 H H','p0 H p0 H','p0 H H p0'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el2*mh2/(4*mw2*st2)
      case (2)
        co(1) = - el2*mh2/(4*mw2*st2)*(                          &
                  + 2*dZe - 2*dst/st + dmh2/mh2                  &
                  + el*dt/(2*st*mw1*mh2) - dmw2/mw2 + dZh + dZp0 &
                  )
      case (3)
        co(1) = el4/(192*pi2*st4)*(                          &
                + 1/2d0*( 1d0 + 1/(2*ct4) )*( 1d0 - 12*lam ) & ! bosonic
                + 3/2d0*mh2/mw2*( 1d0 + 1/(2*ct2) )          & ! bosonic
                + 5/mw4*( ml4(1) + ml4(2) + ml4(3) )         & ! leptonic
                + 5/mw4*Nc*(                                 & ! hadronic
                  + md4(1) + md4(2) + md4(3)                 & ! hadronic
                  + mu4(1) + mu4(2) + mu4(3)                 & ! hadronic
                  )                                          & ! hadronic
                )
      end select
    end select

  case ( 'H H p+ p-','H p+ H p-','H p+ p- H', &
         'H H p- p+','H p- H p+','H p- p+ H', &
         'p+ p- H H','p+ H p- H','p+ H H p-', &
         'p- p+ H H','p- H p+ H','p- H H p+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el2*mh2/(4*mw2*st2)
      case (2)
        co(1) = - el2*mh2/(4*mw2*st2)*(                         &
                  + 2*dZe - 2*dst/st + dmh2/mh2                 &
                  + el*dt/(2*st*mw1*mh2) - dmw2/mw2 + dZh + dZp &
                  )
      case (3)
        co(1) = el4/(64*pi2*st4)*(                           &
                + 1/4d0*( 1d0 + st2/(3*ct2) + st2/(3*ct4) )* & ! bosonic
                  ( 1d0 - 12*lam )                           & ! bosonic
                + mh2/(2*mw2)*( 1d0 + 1/(2*ct2) )            & ! bosonic
                + 5/(3*mw4)*( ml4(1) + ml4(2) + ml4(3) )     & ! leptonic
                + 5/(3*mw4)*Nc*(                             & ! hadronic
                  + mu4(1) + mu4(2) + mu4(3)                 & ! hadronic
                  + md4(1) + md4(2) + md4(3)                 & ! hadronic
                  )                                          & ! hadronic
                )
      end select
    end select

  case ( 'p0 p0 p0 p0' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el2*mh2*3/(4*mw2*st2)
      case (2)
        co(1) = - el2*mh2*3/(4*mw2*st2)*(                    &
                  + 2*dZe - 2*dst/st + dmh2/mh2              &
                  + el*dt/(2*st*mw1*mh2) - dmw2/mw2 + 2*dZp0 &
                  )
      case (3)
        co(1) = el4/(64*pi2*st4)*(                           &
                + 1/2d0*( 1d0 + 1/(2*ct4) )*( 1d0 - 12*lam ) & ! bosonic
                + 3/2d0*mh2/mw2*( 1d0 + 1/(2*ct2) )          & ! bosonic
                + 5/mw4*( ml4(1) + ml4(2) + ml4(3) )         & ! leptonic
                + 5/mw4*Nc*(                                 & ! hadronic
                  + md4(1) + md4(2) + md4(3)                 & ! hadronic
                  + mu4(1) + mu4(2) + mu4(3)                 & ! hadronic
                  )                                          & ! hadronic
                )
      end select
    end select

  case ( 'p0 p0 p+ p-','p0 p+ p0 p-','p0 p+ p- p0', &
         'p0 p0 p- p+','p0 p- p0 p+','p0 p- p+ p0', &
         'p+ p- p0 p0','p+ p0 p- p0','p+ p0 p0 p-', &
         'p- p+ p0 p0','p- p0 p+ p0','p- p0 p0 p+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el2*mh2/(4*mw2*st2)
      case (2)
        co(1) = - el2*mh2/(4*mw2*st2)*(                          &
                  + 2*dZe - 2*dst/st + dmh2/mh2                  &
                  + el*dt/(2*st*mw1*mh2) - dmw2/mw2 + dZp + dZp0 &
                  )
      case (3)
        co(1) = el4/(64*pi2*st4)*(                           &
                + 1/4d0*( 1d0 + st2/(3*ct2) + st2/(3*ct4) )* & ! bosonic
                  ( 1d0 - 12*lam )                           & ! bosonic
                + mh2/(2*mw2)*( 1d0 + 1/(2*ct2) )            & ! bosonic
                + 5/(3*mw4)*( ml4(1) + ml4(2) + ml4(3) )     & ! leptonic
                + 5/(3*mw4)*Nc*(                             & ! hadronic
                  + mu4(1) + mu4(2) + mu4(3)                 & ! hadronic
                  + md4(1) + md4(2) + md4(3)                 & ! hadronic
                  )                                          & ! hadronic
                )
      end select
    end select

  case ( 'p+ p+ p- p-','p+ p- p+ p-','p+ p- p- p+', &
         'p- p- p+ p+','p- p+ p- p+','p- p+ p+ p-'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el2*mh2/(2*mw2*st2)
      case (2)
        co(1) = - el2*mh2/(2*mw2*st2)*(                     &
                  + 2*dZe - 2*dst/st + dmh2/mh2             &
                  + el*dt/(2*st*mw1*mh2) - dmw2/mw2 + 2*dZp &
                  )
      case (3)
        co(1) = el4/(32*pi2*st4)*(                        &
                + ( 1d0 + st4 )*( 1/4d0 - 3*lam )         & ! bosonic
                + ( st2 + 2*st6/ct2 )*( 1/6d0 - 2*lam )   & ! bosonic
                + st8/ct4*( 1/12d0 - lam )                & ! bosonic
                + mh2/(2*mw2)*( 1d0 + 1/(2*ct2) )         & ! bosonic
                + 5/(3*mw4)*( ml4(1) + ml4(2) + ml4(3) )  & ! leptonic
                + 5/(3*mw4)*Nc*(                          & ! hadronic
                  + mu4(1) + mu4(2) + mu4(3)              & ! hadronic
                  + md4(1) + md4(2) + md4(3)              & ! hadronic
                  )                                       & ! hadronic
                )
      end select
    end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 2 scalars and 2 vectors
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! coupling = i*co(1)*g_{mu_v1,mu_v2}
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mu_v1 = Lorentz index of the 1st vector
! mu_v2 = Lorentz index of the 2nd vector
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  case ( 'H H g g','H g H g','H g g H',       &
         'g g H H','g H g H','g H H g',       &
         'p0 p0 g g','p0 g p0 g','p0 g g p0', &
         'g g p0 p0','g p0 g p0','g p0 p0 g', &
         'p+ p- g g','p+ g p- g','p+ g g p-', &
         'p- p+ g g','p- g p+ g','p- g g p+', &
         'g g p+ p-','g p+ g p-','g p+ p- g', &
         'g g p- p+','g p- g p+','g p- p+ g'  )

    select case (ngs)
    case (2)
      select case (xlp)
      case (3)
        co(1) = - el2*gs2/(32*st2*mw2*pi2)*( &
                  + mu2(1) + mu2(2) + mu2(3) &
                  + md2(1) + md2(2) + md2(3) &
                  )
      end select
    end select

  case ( 'H H A A','H A H A','H A A H',       &
         'A A H H','A H A H','A H H A',       &
         'p0 p0 A A','p0 A p0 A','p0 A A p0', &
         'A A p0 p0','A p0 A p0','A p0 p0 A'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (3)
        co(1) = el4/(16*pi2*st2)*(                     &
                + 1/12d0                               & ! bosonic
                - Ql2/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                - Nc/mw2*(                             & ! hadronic
                  + Qu2*( mu2(1) + mu2(2) + mu2(3) )   & ! hadronic
                  + Qd2*( md2(1) + md2(2) + md2(3) )   & ! hadronic
                  )                                    & ! hadronic
                )
      end select
    end select

  case ( 'H H A Z','H A H Z','H A Z H',       &
         'H H Z A','H Z H A','H Z A H',       &
         'A Z H H','A H Z H','A H H Z',       &
         'Z A H H','Z H A H','Z H H A',       &
         'p0 p0 A Z','p0 A p0 Z','p0 A Z p0', &
         'p0 p0 Z A','p0 Z p0 A','p0 Z A p0', &
         'A Z p0 p0','A p0 Z p0','A p0 p0 Z', &
         'Z A p0 p0','Z p0 A p0','Z p0 p0 A'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (2)
        co(1) = el2/(4*st2*ct2)*dZza
      case (3)
        co(1) = el4/(16*pi2*st)*(                   &
                + 1/(12*st2*ct)*( 4d0 + st2 )       & ! bosonic
                + Ql/(mw2*ct)*( I3l/(2*st2) - Ql )* & ! leptonic
                  ( ml2(1) + ml2(2) + ml2(3) )      & ! leptonic
                + Nc/(mw2*ct)*(                     & ! hadronic
                  + Qu*( I3u/(2*st2) - Qu )*        & ! hadronic
                    ( mu2(1) + mu2(2) + mu2(3) )    & ! hadronic
                  + Qd*( I3d/(2*st2) - Qd )*        & ! hadronic
                    ( md2(1) + md2(2) + md2(3) )    & ! hadronic
                  )                                 & ! hadronic
                )
      end select
    end select

  case ( 'H H Z Z','H Z H Z','H Z Z H', &
         'Z Z H H','Z H Z H','Z H H Z'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el2/(2*st2*ct2)
      case (2)
        co(1) = el2/(2*st2*ct2)*(                  &
                + 2*dZe + 2*(st2-ct2)/(ct2*st)*dst &
                + dZzz + dZh                       &
                )
      case (3)
        co(1) = - el4/(16*pi2*ct2)*(                                &
                  + 1/(48*st4*ct2)*( 1d0 + 2*ct2 + 40*ct4 - 4*ct6 ) & ! bosonic
                  + 1/mw2*( Ql2 + 4/3d0*I3l2/st4 - Ql*I3l/st2 )*    & ! leptonic
                    ( ml2(1) + ml2(2) + ml2(3) )                    & ! leptonic
                  + Nc/mw2*(                                        & ! hadronic
                    + ( Qu2 + 4/3d0*I3u2/st4 - Qu*I3u/st2 )*        & ! hadronic
                      ( mu2(1) + mu2(2) + mu2(3) )                  & ! hadronic
                    + ( Qd2 + 4/3d0*I3d2/st4 - Qd*I3d/st2 )*        & ! hadronic
                      ( md2(1) + md2(2) + md2(3) )                  & ! hadronic
                    )                                               & ! hadronic
                  )
      end select
    end select

  case ( 'H H W+ W-','H W+ H W-','H W+ W- H', &
         'H H W- W+','H W- H W+','H W- W+ H', &
         'W+ W- H H','W+ H W- H','W+ H H W-', &
         'W- W+ H H','W- H W+ H','W- H H W+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el2/(2*st2)
      case (2)
        co(1) = el2/(2*st2)*( 2*dZe - 2*dst/st + dZw + dZh )
      case (3)
        co(1) = - el4/(48*pi2*st4)*(                   &
                  + 1/(16*ct2)*( 1d0 + 38*ct2 )        & ! bosonic
                  + 1/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  + Nc/mw2*(                           & ! hadronic
                    + mu2(1) + mu2(2) + mu2(3)         & ! hadronic
                    + md2(1) + md2(2) + md2(3)         & ! hadronic
                    )                                  & ! hadronic
                  )
      end select
    end select

  case ( 'H p+ A W-','H p+ W- A','H A p+ W-', &
         'H A W- p+','H W- p+ A','H W- A p+', &
         'p+ H A W-','p+ H W- A','p+ A H W-', &
         'p+ A W- H','p+ W- H A','p+ W- A H', &
         'A H p+ W-','A H W- p+','A p+ H W-', &
         'A p+ W- H','A W- H p+','A W- p+ H', &
         'W- H p+ A','W- H A p+','W- p+ H A', &
         'W- p+ A H','W- A H p+','W- A p+ H', &
         'H p- A W+','H p- W+ A','H A p- W+', &
         'H A W+ p-','H W+ p- A','H W+ A p-', &
         'p- H A W+','p- H W+ A','p- A H W+', &
         'p- A W+ H','p- W+ H A','p- W+ A H', &
         'A H p- W+','A H W+ p-','A p- H W+', &
         'A p- W+ H','A W+ H p-','A W+ p- H', &
         'W+ H p- A','W+ H A p-','W+ p- H A', &
         'W+ p- A H','W+ A H p-','W+ A p- H'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el2/(2*st)
      case (2)
        co(1) = - el2/(2*st)*(                                        &
                  + 2*dZe - dst/st + 1/2d0*( dZw + dZh + dZaa + dZp ) &
                  )                                                   &
                - el2/(4*ct)*dZza
      case (3)
        co(1) = el4/pi2/(192*st3)*(                  &
                + 1/(4*ct2)*( 1d0 + 22*ct2 )         & ! bosonic
                + 1/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                + Nc/mw2*(                           & ! hadronic
                  + 2*( mu2(1) + mu2(2) + mu2(3) )   & ! hadronic
                  + 3*( md2(1) + md2(2) + md2(3) )   & ! hadronic
                  )                                  & ! hadronic
                )
      end select
    end select

  case ( 'H p+ Z W-','H p+ W- Z','H Z p+ W-', &
         'H Z W- p+','H W- p+ Z','H W- Z p+', &
         'p+ H Z W-','p+ H W- Z','p+ Z H W-', &
         'p+ Z W- H','p+ W- H Z','p+ W- Z H', &
         'Z H p+ W-','Z H W- p+','Z p+ H W-', &
         'Z p+ W- H','Z W- H p+','Z W- p+ H', &
         'W- H p+ Z','W- H Z p+','W- p+ H Z', &
         'W- p+ Z H','W- Z H p+','W- Z p+ H', &
         'H p- Z W+','H p- W+ Z','H Z p- W+', &
         'H Z W+ p-','H W+ p- Z','H W+ Z p-', &
         'p- H Z W+','p- H W+ Z','p- Z H W+', &
         'p- Z W+ H','p- W+ H Z','p- W+ Z H', &
         'Z H p- W+','Z H W+ p-','Z p- H W+', &
         'Z p- W+ H','Z W+ H p-','Z W+ p- H', &
         'W+ H p- Z','W+ H Z p-','W+ p- H Z', &
         'W+ p- Z H','W+ Z H p-','W+ Z p- H'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - el2/(2*ct)
      case (2)
        co(1) = - el2/(2*ct)*(                                        &
                  + 2*dZe - dct/ct + 1/2d0*( dZw + dZh + dZzz + dZp ) &
                  )                                                   &
                - el2/(4*st)*dZaz
      case (3)
        co(1) = el4/(pi2*ct)/(192*st2)*(                 &
                + 1/(4*ct2*st2)*( 1d0 + 21*ct2 - 22*ct4) & ! bosonic
                + 1/mw2*( ml2(1) + ml2(2) + ml2(3) )     & ! leptonic
                + Nc/mw2*(                               & ! hadronic
                  + 2*( mu2(1) + mu2(2) + mu2(3) )       & ! hadronic
                  + 3*( md2(1) + md2(2) + md2(3) )       & ! hadronic
                  )                                      & ! hadronic
                )
      end select
    end select

  case ( 'p0 p0 Z Z','p0 Z p0 Z','p0 Z Z p0', &
         'Z Z p0 p0','Z p0 Z p0','Z p0 p0 Z'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el2/(2*st2*ct2)
      case (2)
        co(1) = el2/(2*st2*ct2)*                                   &
                ( 2*dZe + 2*(st2-ct2)/(ct2*st)*dst + dZzz + dZp0 )
      case (3)
        co(1) = - el4/(16*pi2*ct2)*(                                &
                  + 1/(48*st4*ct2)*( 1d0 + 2*ct2 + 40*ct4 - 4*ct6 ) & ! bosonic
                  + 1/mw2*( Ql2 + 4/3d0*I3l2/st4 - Ql*I3l/st2 )*    & ! leptonic
                    ( ml2(1) + ml2(2) + ml2(3) )                    & ! leptonic
                  + Nc/mw2*(                                        & ! hadronic
                    + ( Qu2 + 4/3d0*I3u2/st4 - Qu*I3u/st2 )*        & ! hadronic
                      ( mu2(1) + mu2(2) + mu2(3) )                  & ! hadronic
                    + ( Qd2 + 4/3d0*I3d2/st4 - Qd*I3d/st2 )*        & ! hadronic
                      ( md2(1) + md2(2) + md2(3) )                  & ! hadronic
                    )                                               & ! hadronic
                  )
      end select
    end select

  case ( 'p0 p0 W+ W-','p0 W+ p0 W-','p0 W+ W- p0', &
         'p0 p0 W- W+','p0 W- p0 W+','p0 W- W+ p0', &
         'W+ W- p0 p0','W+ p0 W- p0','W+ p0 p0 W-', &
         'W- W+ p0 p0','W- p0 W+ p0','W- p0 p0 W+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el2/(2*st2)
      case (2)
        co(1) = el2/(2*st2)*( 2*dZe - 2*dst/st + dZw + dZp0 )
      case (3)
        co(1) = - el4/(48*pi2*st4)*(                   &
                  + 1/(16*ct2)*( 1d0 + 38*ct2 )        & ! bosonic
                  + 1/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  + Nc/mw2*(                           & ! hadronic
                    + mu2(1) + mu2(2) + mu2(3)         & ! hadronic
                    + md2(1) + md2(2) + md2(3)         & ! hadronic
                    )                                  & ! hadronic
                  )
      end select
    end select

  case ( 'p0 p+ A W-','p0 p+ W- A','p0 A p+ W-', &
         'p0 A W- p+','p0 W- p+ A','p0 W- A p+', &
         'p+ p0 A W-','p+ p0 W- A','p+ A p0 W-', &
         'p+ A W- p0','p+ W- p0 A','p+ W- A p0', &
         'A p0 p+ W-','A p0 W- p+','A p+ p0 W-', &
         'A p+ W- p0','A W- p0 p+','A W- p+ p0', &
         'W- p0 p+ A','W- p0 A p+','W- p+ p0 A', &
         'W- p+ A p0','W- A p0 p+','W- A p+ p0'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = cId0*el2/(2*st)
      case (2)
        co(1) = + cId0*el2/(2*st)*(                                    &
                  + 2*dZe - dst/st + 1/2d0*( dZw + dZaa + dZp + dZp0 ) &
                  )                                                    &
                + cId0*el2/(4*ct)*dZza
      case (3)
        co(1) = - cId0*el4/pi2/(192*st3)*(             &
                  + 1/(4*ct2)*( 1d0 + 22*ct2 )         & ! bosonic
                  + 1/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  + Nc/mw2*(                           & ! hadronic
                    + 2*( mu2(1) + mu2(2) + mu2(3) )   & ! hadronic
                    + 3*( md2(1) + md2(2) + md2(3) )   & ! hadronic
                    )                                  & ! hadronic
                  )
      end select
    end select

  case ( 'p0 p- A W+','p0 p- W+ A','p0 A p- W+', &
         'p0 A W+ p-','p0 W+ p- A','p0 W+ A p-', &
         'p- p0 A W+','p- p0 W+ A','p- A p0 W+', &
         'p- A W+ p0','p- W+ p0 A','p- W+ A p0', &
         'A p0 p- W+','A p0 W+ p-','A p- p0 W+', &
         'A p- W+ p0','A W+ p0 p-','A W+ p- p0', &
         'W+ p0 p- A','W+ p0 A p-','W+ p- p0 A', &
         'W+ p- A p0','W+ A p0 p-','W+ A p- p0'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - cId0*el2/(2*st)
      case (2)
        co(1) = - cId0*el2/(2*st)*(                                    &
                  + 2*dZe - dst/st + 1/2d0*( dZw + dZaa + dZp + dZp0 ) &
                  )                                                    &
                - cId0*el2/(4*ct)*dZza
      case (3)
        co(1) = + cId0*el4/pi2/(192*st3)*(             &
                  + 1/(4*ct2)*( 1d0 + 22*ct2 )         & ! bosonic
                  + 1/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  + Nc/mw2*(                           & ! hadronic
                    + 2*( mu2(1) + mu2(2) + mu2(3) )   & ! hadronic
                    + 3*( md2(1) + md2(2) + md2(3) )   & ! hadronic
                    )                                  & ! hadronic
                  )
      end select
    end select

  case ( 'p0 p+ Z W-','p0 p+ W- Z','p0 Z p+ W-', &
         'p0 Z W- p+','p0 W- p+ Z','p0 W- Z p+', &
         'p+ p0 Z W-','p+ p0 W- Z','p+ Z p0 W-', &
         'p+ Z W- p0','p+ W- p0 Z','p+ W- Z p0', &
         'Z p0 p+ W-','Z p0 W- p+','Z p+ p0 W-', &
         'Z p+ W- p0','Z W- p0 p+','Z W- p+ p0', &
         'W- p0 p+ Z','W- p0 Z p+','W- p+ p0 Z', &
         'W- p+ Z p0','W- Z p0 p+','W- Z p+ p0'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = + cId0*el2/(2*ct)
      case (2)
        co(1) = + cId0*el2/(2*ct)*(                                    &
                  + 2*dZe - dct/ct + 1/2d0*( dZw + dZzz + dZp + dZp0 ) &
                  )                                                    &
                + cId0*el2/(4*st)*dZaz
      case (3)
        co(1) = - cId0*el4/(pi2*ct)/(192*st2)*(             &
                  + 1/(4*ct2*st2)*( 1d0 + 21*ct2 - 22*ct4 ) & ! bosonic
                  + 1/mw2*( ml2(1) + ml2(2) + ml2(3) )      & ! leptonic
                  + Nc/mw2*(                                & ! hadronic
                    + 2*( mu2(1) + mu2(2) + mu2(3) )        & ! hadronic
                    + 3*( md2(1) + md2(2) + md2(3) )        & ! hadronic
                    )                                       & ! hadronic
                  )
      end select
    end select

  case ( 'p0 p- Z W+','p0 p- W+ Z','p0 Z p- W+', &
         'p0 Z W+ p-','p0 W+ p- Z','p0 W+ Z p-', &
         'p- p0 Z W+','p- p0 W+ Z','p- Z p0 W+', &
         'p- Z W+ p0','p- W+ p0 Z','p- W+ Z p0', &
         'Z p0 p- W+','Z p0 W+ p-','Z p- p0 W+', &
         'Z p- W+ p0','Z W+ p0 p-','Z W+ p- p0', &
         'W+ p0 p- Z','W+ p0 Z p-','W+ p- p0 Z', &
         'W+ p- Z p0','W+ Z p0 p-','W+ Z p- p0'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = - cId0*el2/(2*ct)
      case (2)
        co(1) = - cId0*el2/(2*ct)*(                                    &
                  + 2*dZe - dct/ct + 1/2d0*( dZw + dZzz + dZp + dZp0 ) &
                  )                                                    &
                - cId0*el2/(4*st)*dZaz
      case (3)
        co(1) = + cId0*el4/(pi2*ct)/(192*st2)*(             &
                  + 1/(4*ct2*st2)*( 1d0 + 21*ct2 - 22*ct4 ) & ! bosonic
                  + 1/mw2*( ml2(1) + ml2(2) + ml2(3) )      & ! leptonic
                  + Nc/mw2*(                                & ! hadronic
                    + 2*( mu2(1) + mu2(2) + mu2(3) )        & ! hadronic
                    + 3*( md2(1) + md2(2) + md2(3) )        & ! hadronic
                    )                                       & ! hadronic
                  )
      end select
    end select

  case ( 'p+ p- A A','p+ A p- A','p+ A A p-', &
         'p- p+ A A','p- A p+ A','p- A A p+', &
         'A A p+ p-','A p+ A p-','A p+ p- A', &
         'A A p- p+','A p- A p+','A p- p+ A'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = 2*el2
      case (2)
        co(1) = 2*el2*( 2*dZe + (st2-ct2)/(2*stct)*dZza + dZaa + dZp )
      case (3)
        co(1) = - el4/(12*pi2*st2)*(                   &
                  + 1/(16*ct2)*( 1d0 + 21*ct2 )        & ! bosonic
                  + 1/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  + 5/6d0*Nc/mw2*(                     & ! hadronic
                    + mu2(1) + mu2(2) + mu2(3)         & ! hadronic
                    + md2(1) + md2(2) + md2(3)         & ! hadronic
                    )                                  & ! hadronic
                  )
      end select
    end select

  case ( 'p+ p- A Z','p+ p- Z A','p+ A p- Z', &
         'p+ A Z p-','p+ Z p- A','p+ Z A p-', &
         'p- p+ A Z','p- p+ Z A','p- A p+ Z', &
         'p- A Z p+','p- Z p+ A','p- Z A p+', &
         'A p+ p- Z','A p+ Z p-','A p- p+ Z', &
         'A p- Z p+','A Z p+ p-','A Z p- p+', &
         'Z p+ p- A','Z p+ A p-','Z p- p+ A', &
         'Z p- A p+','Z A p+ p-','Z A p- p+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el2*(st2-ct2)/stct
      case (2)
        co(1) = el2*(                                   &
                + (st2-ct2)/stct*(                      &
                  2*dZe + 1/2d0*dZzz + 1/2d0*dZaa + dZp &
                  )                                     &
                + dst/(st2*ct3)                         &
                + (st2-ct2)**2/(4*st2*ct2)*dZza + dZaz  &
                )
      case (3)
        co(1) = el4/(12*pi2*stct)*(                            &
                + ( 42*ct4 - 10*ct2 - 1d0 )/(32*st2*ct2)       & ! bosonic
                - 1/mw2*( Ql2 + 5*Ql*I3n/(8*st2) )*            & ! leptonic
                  ( ml2(1) + ml2(2) + ml2(3) )                 & ! leptonic
                - Nc/mw2*(                                     & ! hadronic
                  + ( 5/6d0 - Qd*I3d/st2 + 5/8d0*Qu*I3d/st2 )* & ! hadronic
                    ( mu2(1) + mu2(2) + mu2(3) )               & ! hadronic
                  + ( 5/6d0 - Qu*I3u/st2 + 5/8d0*Qd*I3u/st2 )* & ! hadronic
                    ( md2(1) + md2(2) + md2(3) )               & ! hadronic
                  )                                            & ! hadronic
                )
      end select
    end select

  case ( 'p+ p- Z Z','p+ Z p- Z','p+ Z Z p-', &
         'p- p+ Z Z','p- Z p+ Z','p- Z Z p+', &
         'Z Z p+ p-','Z p+ Z p-','Z p+ p- Z', &
         'Z Z p- p+','Z p- Z p+','Z p- p+ Z'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el2*(st2-ct2)**2/(2*st2*ct2)
      case (2)
        co(1) = el2*(                                             &
                + (st2-ct2)**2/(2*st2*ct2)*( 2*dZe + dZzz + dZp ) &
                + (st2-ct2)/(st3*ct4)*dst                         &
                + (st2-ct2)/stct*dZaz                             &
                )
      case (3)
        co(1) = el4/(12*pi2*ct2)*(                                 &
                - 1/(64*st4*ct2)*( 1d0 - 2*ct2 - 44*ct4 + 84*ct6 ) & ! bosonic
                - 1/mw2*( Ql2 + 5/4d0*Ql*I3n/st2 + I3n2/st4 )*     & ! leptonic
                  ( ml2(1) + ml2(2) + ml2(3) )                     & ! leptonic
                - Nc/mw2*(                                         & ! hadronic
                  + (                                              & ! hadronic
                    + 5/6d0                                        & ! hadronic
                    - 2*Qd*I3d/st2 + 5/4d0*Qu*I3d/st2              & ! hadronic
                    + I3d2/st4                                     & ! hadronic
                    )*( mu2(1) + mu2(2) + mu2(3) )                 & ! hadronic
                  + (                                              & ! hadronic
                    + 5/6d0                                        & ! hadronic
                    - 2*Qu*I3u/st2 + 5/4d0*Qd*I3u/st2              & ! hadronic
                    + I3u2/st4                                     & ! hadronic
                    )*( md2(1) + md2(2) + md2(3) )                 & ! hadronic
                  )                                                & ! hadronic
                )
      end select
    end select

  case ( 'p+ p- W+ W-','p+ p- W- W+','p+ W+ p- W-', &
         'p+ W+ W- p-','p+ W- p- W+','p+ W- W+ p-', &
         'p- p+ W+ W-','p- p+ W- W+','p- W+ p+ W-', &
         'p- W+ W- p+','p- W- p+ W+','p- W- W+ p+', &
         'W+ p+ p- W-','W+ p+ W- p-','W+ p- p+ W-', &
         'W+ p- W- p+','W+ W- p+ p-','W+ W- p- p+', &
         'W- p+ p- W+','W- p+ W+ p-','W- p- p+ W+', &
         'W- p- W+ p+','W- W+ p+ p-','W- W+ p- p+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        co(1) = el2/(2*st2)
      case (2)
        co(1) = el2/(2*st2)*( 2*dZe - 2*dst/st + dZw + dZp )
      case (3)
        co(1) = - el4/(48*pi2*st4)*(                   &
                  + 1/(16*ct2)*( 1d0 + 38*ct2 )        & ! bosonic
                  + 1/mw2*( ml2(1) + ml2(2) + ml2(3) ) & ! leptonic
                  + Nc/mw2*(                           & ! hadronic
                    + mu2(1) + mu2(2) + mu2(3)         & ! hadronic
                    + md2(1) + md2(2) + md2(3)         & ! hadronic
                    )                                  & ! hadronic
                  )
      end select
    end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 4 vectors
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! coupling = i*(
!            + co(1)*g_{mu_v1,mu_v4}*g_{mu_v2,mu_v3}
!            + co(2)*g_{mu_v1,mu_v3}*g_{mu_v2,mu_v4}
!            + co(3)*g_{mu_v1,mu_v2}*g_{mu_v3,mu_v4}
!            )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! mu_v1 = Lorentz index of the 1st vector
! mu_v2 = Lorentz index of the 2nd vector
! mu_v3 = Lorentz index of the 3rd vector
! mu_v4 = Lorentz index of the 4th vector
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! The vectors are ordered in the following way:
!
!        | v1
!        V
!        |
! --->---0---<---
! v2     |     v3
!        ^
!        | v4
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  case ( 'g g g g' )
    ! Note that the rational part involves another structure also

    select case (ngs)
    case (2)
      select case (xlp)
      case (0)
        c0 = gs2/2d0
        co(1:3) = c0
      end select
    case (4)
      select case (xlp)
      case (2)
        c0 = gs2*( dZgs + dZg )
        co(1:3) = c0
      case (3)
        c0 = - gs4/(48*pi2)
        co(1:3) = c0
      end select
    end select

  case ( 'g g g A','g g A g','g A g g','A g g g' )
    ! Rational part, the axial part vanishes when summing over quarks

    select case (ngs)
    case (3)
      select case (xlp)
      case (3)
        c0 = - gs3*el*csq2/(32*pi2)*( Qu + Qd )
        co(1:3) = c0
      end select
    end select

  case ( 'g g g Z','g g Z g','g Z g g','Z g g g' )
    ! Rational part, the axial part vanishes when summing over quarks

    select case (ngs)
    case (3)
      select case (xlp)
      case (3)
        c0 = gs3*el*csq2/(64*pi2*stct)*         &
           ( I3u - 2*st2*Qu + I3d - 2*st2*Qd )
        co(1:3) = c0
      end select
    end select

  case ( 'g g A A','g A g A','g A A g', &
         'A g g A','A g A g','A A g g'  )

    select case (ngs)
    case (2)
      select case (xlp)
      case (3)
        c0 = gs2*el2/(8*pi2)*( Qu2 + Qd2 )
        co(1:3) = c0
      end select
    end select

  case ( 'g g A Z','g g Z A','g A g Z', &
         'g A Z g','g Z g A','g Z A g', &
         'A g g Z','A g Z g','A Z g g', &
         'Z g g A','Z g A g','Z A g g'  )

    select case (ngs)
    case (2)
      select case (xlp)
      case (3)
        c0 = - gs2*el2/(16*pi2)/stct*                      &
             ( Qu*I3u - 2*st2*Qu2 + Qd*I3d - 2*st2*Qd2 )
        co(1:3) = c0
      end select
    end select

  case ( 'g g Z Z','g Z g Z','g Z Z g', &
         'Z g g Z','Z g Z g','Z Z g g'  )

    select case (ngs)
    case (2)
      select case (xlp)
      case (3)
        c0 = gs2*el2/(32*pi2*st2*ct2)*(     &
           + ( I3u - 2*st2*Qu )**2 + I3u2 &
           + ( I3d - 2*st2*Qd )**2 + I3d2 &
           )
        co(1:3) = c0
      end select
    end select

  case ( 'g g W+ W-','g g W- W+','g W+ g W-', &
         'g W+ W- g','g W- g W+','g W- W+ g', &
         'W+ g g W-','W+ g W- g','W+ W- g g', &
         'W- g g W+','W- g W+ g','W- W+ g g'  )

    select case (ngs)
    case (2)
      select case (xlp)
      case (3)
        c0 = gs2*el2/(32*pi2*st2)
        co(1:3) = c0
      end select
    end select

  case ( 'A A A A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (3)
        c0 = el4/(12*pi2)*(    &
           - 1d0               & ! bosonic
           + 3*Ql4             & ! leptonic
           + 3*Nc*( Qu4 + Qd4) & ! hadronic
           )
        co(1:3) = c0
        cQED(1:3) = el4/(12*pi2)*(      &
                    + 3*Ql4             & ! leptonic
                    + 3*Nc*( Qu4 + Qd4) & ! hadronic
                    )
      end select
    end select

  case ( 'A A A Z','A A Z A','A Z A A','Z A A A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (3)
        c0 = el4/(12*pi2)*(                     &
           + ct/st                            & ! bosonic
           + 3/ct*( st*Ql4 - Ql3*I3l/(2*st) ) & ! leptonic
           + 3*Nc/ct*(                        & ! hadronic
             + st*Qu4 - Qu3*I3u/(2*st)        & ! hadronic
             + st*Qd4 - Qd3*I3d/(2*st)        & ! hadronic
             )                                & ! hadronic
           )
        co(1:3) = c0
      end select
    end select

  case ( 'A A Z Z','A Z A Z','A Z Z A', &
         'Z A A Z','Z A Z A','Z Z A A'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (3)
        c0 = - el4/(12*pi2)*(                                      &
             + ct2/st2                                           & ! bosonic
             - 3/(2*ct2)*( st2*Ql4 + ( st*Ql2 - Ql*I3l/st )**2 ) & ! leptonic
             - 3/(2*ct2)*Nc*(                                    & ! hadronic
               + st2*Qu4 + ( st*Qu2 - Qu*I3u/st )**2             & ! hadronic
               + st2*Qd4 + ( st*Qd2 - Qd*I3d/st )**2             & ! hadronic
               )                                                 & ! hadronic
             )
        co(1:3) = c0
      end select
    end select

  case ( 'A Z Z Z','Z A Z Z','Z Z A Z','Z Z Z A' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (3)
        c0 = el4/(12*pi2)*(                   &
           + ct3/st3                        & ! bosonic
           + 3/(2*ct3)*Ql*(                 & ! leptonic
             + 2*Ql3*st3 - 3*Ql2*I3l*st     & ! leptonic
             + 3*Ql*I3l2/st - I3l3/st3      & ! leptonic
             )                              & ! leptonic
           + 3/(2*ct3)*Nc*(                 & ! hadronic
             + 2*Qu4*st3 - 3*Qu3*I3u*st     & ! hadronic
             + 3*Qu2*I3u2/st - Qu*I3u3/st3  & ! hadronic
             + 2*Qd4*st3 - 3*Qd3*I3d*st     & ! hadronic
             + 3*Qd2*I3d2/st - Qd*I3d3/st3  & ! hadronic
             )                              & ! hadronic
           )
        co(1:3) = c0
      end select
    end select

  case ( 'Z Z Z Z' )

    select case (ngs)
    case (0)
      select case (xlp)
      case (3)
        c0 = - el4/(12*pi2)*(                 &
             + ct4/st4                      & ! bosonic
             - 3/ct4*(                      & ! leptonic
               + I3n4/(2*st4)               & ! leptonic
               + Ql4*st4 - 2*Ql3*I3l*st2    & ! leptonic
               + 3*Ql2*I3l2 - 2*Ql*I3l3/st2 & ! leptonic
               + I3l4/(2*st4)               & ! leptonic
               )                            & ! leptonic
             - 3*Nc/ct4*(                   & ! hadronic
               + Qu4*st4 - 2*Qu3*I3u*st2    & ! hadronic
               + 3*Qu2*I3u2 - 2*Qu*I3u3/st2 & ! hadronic
               + I3u4/(2*st4)               & ! hadronic
               + Qd4*st4 - 2*Qd3*I3d*st2    & ! hadronic
               + 3*Qd2*I3d2 - 2*Qd*I3d3/st2 & ! hadronic
               + I3d4/(2*st4)               & ! hadronic
               )                            & ! hadronic
             )
        co(1:3) = c0
      end select
    end select

  case ( 'A A W+ W-','A A W- W+','A W+ A W-', &
         'A W+ W- A','A W- A W+','A W- W+ A', &
         'W+ A A W-','W+ A W- A','W+ W- A A', &
         'W- A A W+','W- A W+ A','W- W+ A A'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        c0 = - el2
        cA = + 2*c0
        cB = - c0
      case (2)
        c0 = - el2*( 2*dZe + dZw + dZaa ) + el2*ct/st*dZza
        cA = + 2*c0
        cB = - c0
      case (3)
        cA = el4/(16*pi2*st2)*(       &
             + 1/3d0*( 10d0 + 4*lam ) &
             + 3d0                    &
             + 25/9d0*Nc              &
             )
        cB = - el4/(16*pi2*st2)*(      &
               + 1/3d0*( 7d0 + 2*lam ) &
               + 1d0                   &
               + 11/9d0*Nc             &
               )
      end select
      select case (cv)
      case ( 'A W+ W- A','A W- W+ A','W+ A A W-','W- A A W+' )
        co(1) = cA
        co(2) = cB
        co(3) = cB
      case ( 'A W+ A W-','A W- A W+','W+ A W- A','W- A W+ A' )
        co(1) = cB
        co(2) = cA
        co(3) = cB
      case ( 'A A W+ W-','A A W- W+','W+ W- A A','W- W+ A A' )
        co(1) = cB
        co(2) = cB
        co(3) = cA
      end select
    end select

  case ( 'A Z W+ W-','A Z W- W+','A W+ Z W-', &
         'A W+ W- Z','A W- Z W+','A W- W+ Z', &
         'Z A W+ W-','Z A W- W+','Z W+ A W-', &
         'Z W+ W- A','Z W- A W+','Z W- W+ A', &
         'W+ A Z W-','W+ A W- Z','W+ Z A W-', &
         'W+ Z W- A','W+ W- A Z','W+ W- Z A', &
         'W- A Z W+','W- A W+ Z','W- Z A W+', &
         'W- Z W+ A','W- W+ A Z','W- W+ Z A'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        c0 = el2*ct/st
        cA = + 2*c0
        cB = - c0
      case (2)
        c0 = + el2*ct/st*(                                          &
               + 2*dZe - dst/(ct2*st) + dZw + 1/2d0*( dZzz + dZaa ) &
               )                                                    &
             - el2/2d0*( dZaz + ct2/st2*dZza )
        cA = + 2*c0
        cB = - c0
      case (3)
        cA = + el4/(16*pi2*stct)*(           &
               - ct2/(3*st2)*( 10d0 + 4*lam) & ! bosonic
               + 3*( 1d0 - 11/(12*st2) )     & ! leptonic
               + Nc*( 25/9d0 - 11/(4*st2) )  & ! hadronic
               )
        cB = + el4/(16*pi2*stct)*(           &
               + ct2/(3*st2)*( 7d0 + 2*lam)  & ! bosonic
               - 1d0 + 5/(4*st2)             & ! leptonic
               - Nc*( 11/9d0 - 5/(4*st2) )   & ! hadronic
               )
      end select
      select case (cv)
      case ( 'A W+ W- Z','A W- W+ Z','Z W+ W- A','Z W- W+ A', &
             'W+ A Z W-','W+ Z A W-','W- A Z W+','W- Z A W+' )
        co(1) = cA
        co(2) = cB
        co(3) = cB
      case ( 'A W+ Z W-','A W- Z W+','Z W+ A W-','Z W- A W+', &
             'W+ A W- Z','W+ Z W- A','W- A W+ Z','W- Z W+ A'  )
        co(1) = cB
        co(2) = cA
        co(3) = cB
      case ( 'A Z W+ W-','A Z W- W+','Z A W+ W-','Z A W- W+', &
             'W+ W- A Z','W+ W- Z A','W- W+ A Z','W- W+ Z A'  )
        co(1) = cB
        co(2) = cB
        co(3) = cA
      end select
    end select

  case ( 'Z Z W+ W-','Z Z W- W+','Z W+ Z W-', &
         'Z W+ W- Z','Z W- Z W+','Z W- W+ Z', &
         'W+ Z Z W-','W+ Z W- Z','W+ W- Z Z', &
         'W- Z Z W+','W- Z W+ Z','W- W+ Z Z'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        c0 = - el2*ct2/st2
        cA = + 2*c0
        cB = - c0
      case (2)
        c0 = - el2*ct2/st2*( 2*dZe - 2/(ct2*st)*dst + dZw + dZzz ) &
             + el2*ct/st*dZaz
        cA = + 2*c0
        cB = - c0
      case (3)
        cA = + el4/pi2*(                                   &
               + ct2/(24*st4)*( 5d0 + 2*lam )              &
               + 1/(16*ct2)*(                              &
                 + 3*( 1d0 - 11/(6*st2) + 11/(12*st4) )    &
                 + Nc*( 25/9d0 - 11/(2*st2) + 11/(4*st4) ) &
                 )                                         &
               )
        cB = + el4/pi2*(                                   &
               - ct2/(48*st4)*( 7d0 + 2*lam )              &
               + 1/(16*ct2)*(                              &
                 - 1d0 + 5/(2*st2) - 5/(4*st4)             &
                 + Nc*( - 11/9d0 + 5/(2*st2) - 5/(4*st4) ) &
                 )                                         &
               )
      end select
      select case (cv)
      case ( 'Z W+ W- Z','Z W- W+ Z','W+ Z Z W-','W- Z Z W+' )
        co(1) = cA
        co(2) = cB
        co(3) = cB
      case ( 'Z W+ Z W-','Z W- Z W+','W+ Z W- Z','W- Z W+ Z' )
        co(1) = cB
        co(2) = cA
        co(3) = cB
      case ( 'Z Z W+ W-','Z Z W- W+','W+ W- Z Z','W- W+ Z Z' )
        co(1) = cB
        co(2) = cB
        co(3) = cA
      end select
    end select

  case ( 'W+ W+ W- W-','W+ W- W+ W-','W+ W- W- W+', &
         'W- W+ W+ W-','W- W+ W- W+','W- W- W+ W+'  )

    select case (ngs)
    case (0)
      select case (xlp)
      case (0)
        c0 = el2/st2
        cA = + 2*c0
        cB = - c0
      case (2)
        c0 = el2/st2*2*( dZe - dst/st + dZw )
        cA = + 2*c0
        cB = - c0
      case (3)
        cA = - el4/(8*pi2*st4)*(       &
             + 1/3d0*( 7d0 + 2*lam ) & ! bosonic
             + 5/4d0                 & ! leptonic
             + 5/4d0*Nc              & ! hadronic
             )
        cB = + el4/(16*pi2*st4)*(      &
             + 1/3d0*( 3d0 + 2*lam ) & ! bosonic
             + 3/2d0                 & ! leptonic
             + 3/2d0*Nc              & ! hadronic
             )
      end select
      select case (cv)
      case ( 'W+ W- W- W+','W- W+ W+ W-' )
        co(1) = cA
        co(2) = cB
        co(3) = cB
      case ( 'W+ W- W+ W-','W- W+ W- W+' )
        co(1) = cB
        co(2) = cA
        co(3) = cB
      case ( 'W+ W+ W- W-','W- W- W+ W+' )
        co(1) = cB
        co(2) = cB
        co(3) = cA
      end select
    end select

  end select

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  select case (xlp)
  case (3)
    if (.not.loopWEAK) co = c0d0
    if (loopQED)       co = co + cP
    if (pureQED)       co = cQED
  end select


  if (n.eq.3) co = co * coupling3(ff(1),ff(2),ff(3))
  if (n.eq.4) co = co * coupling4(ff(1),ff(2),ff(3),ff(4))

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  end function cou

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module model_vertices_rcl

!#####################################################################



