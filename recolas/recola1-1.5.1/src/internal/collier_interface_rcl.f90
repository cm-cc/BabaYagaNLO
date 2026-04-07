!#####################################################################
!!
!!  File  collier_interface_rcl.f90
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

  module collier_interface_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl
  use collier

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine initialize_collier

  integer                  :: pr,nPropMax,i,n
  integer, allocatable     :: nProp(:), nCacheTmp(:)
  logical, allocatable     :: cacheOnTmp(:)
  complex(dp), allocatable :: m2(:)

  allocate (nProp(prTot)); nProp = 0
  do pr = 1,prTot
    if (loop(pr)) nProp(pr) = maxval(legsti(:tiTot(pr),pr))
  enddo
  nPropMax = maxval(nProp)

  if (nPropMax.eq.0) return

  if (trim(collier_output_dir) .eq. 'default') then
    call Init_cll(nPropMax,noreset=.true.)
  else
    call Init_cll(nPropMax,                             &
                  folder_name=trim(collier_output_dir), &
                  noreset=.true.)
  endif

  if (.not.allocated(nCache)) then
    allocate(nCache(1:prTot))
    nCache = 1
  else if (size(nCache) .ne. prTot) then
    allocate(nCacheTmp(size(nCache)))
    nCacheTmp(:) = nCache(:)
    deallocate(nCache)
    allocate(nCache(prTot))
    nCache = 1
    nCache(1:size(nCacheTmp)) = nCacheTmp(:)
  endif

  if (.not.allocated(cacheOn)) then
    allocate(cacheOn(1:prTot))
    cacheOn = .true.
  else if (size(cacheOn) .ne. prTot) then
    allocate(cacheOnTmp(size(cacheOn)))
    cacheOnTmp(:) = cacheOn(:)
    deallocate(cacheOn)
    allocate(cacheOn(prTot))
    cacheOn = .true.
    cacheOn(1:size(cacheOnTmp)) = cacheOnTmp(:)
  endif

  allocate(nCacheTot(0:prTot),tiCache(prTot))

  nCacheTot = 0
  do pr = 1,prTot
    if ((.not. prexists(pr)) .or. (.not. loop(pr))) then
      cacheOn(pr) = .false.
      nCache(pr) = 0
    end if

    if (loop(pr)) then
      if (nCache(pr) .gt. 0) then
        tiCache(pr) = ceiling(tiTot(pr)*1d0/nCache(pr))
      else
        tiCache(pr) = 0
      end if
    end if
    nCacheTot(pr) = nCacheTot(pr-1) + nCache(pr)
!    write(nx,*) 'nCache  =',nCache(pr)
!    write(nx,*) 'tiTot   =',tiTot(pr)
!    write(nx,*) 'tiCache =',tiCache(pr)
!    write(nx,*)
  enddo

  call InitCacheSystem_cll(nCacheTot(prTot),nPropMax)

  call SetMode_cll(collier_mode)

  call SwitchOnTenRed_cll
!  call SwitchOffTenRed_cll

  do pr = 1,prTot
    do n = nCacheTot(pr-1)+1,nCacheTot(pr)
      call SetCacheLevel_cll(n,nProp(pr))
      if (.not.cacheOn(pr)) call SwitchOffCache_cll(n)
    enddo
  enddo

  if (ifail.ge.0) call SetErrStop_cll(-11)

  n = 0
  do i = 23,31
    if (regf(i).eq.2) n = n + 1
  enddo
  allocate (m2(n))
  n = 0
  do i = 23,31
    if (regf(i).eq.2) then
      n = n + 1
      if (complex_mass_scheme.eq.1) then
        m2(n) = cm2n(nmf(i))
      else
        m2(n) = real(cm2n(nmf(i)),kind=dp)*c1d0
      endif
    endif
  enddo
  call setminf2_cll (n,m2)
  deallocate (m2)

  call SetDeltaUV_cll (deltaUV)
  call SetDeltaIR_cll (deltaIR,deltaIR2)

  if (reg_soft.eq.1) then
    call SetMuIR2_cll(muIR**2)
  else
    call SetMuIR2_cll(lambda**2)
  endif
  call SetMuUV2_cll(muUV**2)

  deallocate(nProp)

  end subroutine initialize_collier

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine initialize_collier_ct

  integer                  :: i,n
  complex(dp), allocatable :: m2(:)

  if (trim(collier_output_dir) .eq. 'default') then
    call Init_cll(2,noreset=.true.)
  else
    call Init_cll(2,folder_name=trim(collier_output_dir), &
                    noreset=.true.)
  endif

  call SetMode_cll(collier_mode)
  n = 0
  do i = 23,31
    if (regf(i).eq.2) n = n + 1
  enddo
  allocate (m2(n))
  n = 0
  do i = 23,31
    if (regf(i).eq.2) then
      n = n + 1
      if (complex_mass_scheme.eq.1) then
        m2(n) = cm2n(nmf(i))
      else
        m2(n) = real(cm2n(nmf(i)),kind=dp)*c1d0
      endif
    endif
  enddo
  call setminf2_cll (n,m2)
  deallocate (m2)

  call SetDeltaUV_cll (deltaUV)
  call SetDeltaIR_cll (deltaIR,deltaIR2)

  if (reg_soft.eq.1) then
    call SetMuIR2_cll(muIR**2)
  else
    call SetMuIR2_cll(lambda**2)
  endif
  call SetMuUV2_cll(muUV**2)

  end subroutine initialize_collier_ct

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module collier_interface_rcl

!#####################################################################
