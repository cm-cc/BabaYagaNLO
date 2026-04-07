!#####################################################################
!!
!!  File  reset_rcl.f90
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

  module reset_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl
  use amplitude_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine reset_recola_rcl

  ! This subroutine deallocates all global allocatable arrays and sets
  ! internal variables to the initialization value.
  ! It has to be called at the end of the Recola-session.

  integer :: i,n
  logical :: fileopen

  ! internal_variables
  deallocate (timeTI,timeTC,inpr,legsIn,legsOut,par,hel,resMax, &
              binRes,parRes,powgs,loop,process,qflow,           &
              refscheme,Nalpha,NLOscheme,NPh,NonPh,NoffPh,      &
              ve2ct,ve2r2,ve3tr,ve3ct,ve3r2,ve4tr,ve4ct,ve4r2,  &
              polproj,polprojin,prexists)

  if (allocated(momenta))      deallocate(momenta)
  if (allocated(matrixLO))     deallocate(matrixLO)
  if (allocated(matrix))       deallocate(matrix)
  if (allocated(matrix2))      deallocate(matrix2)
  if (allocated(matrix2h))     deallocate(matrix2h)
  if (allocated(matrix2int))   deallocate(matrix2int)
  if (allocated(matrix2cc))    deallocate(matrix2cc)
  if (allocated(matrix2ccint)) deallocate(matrix2ccint)
  if (allocated(matrix2ccnlo)) deallocate(matrix2ccnlo)
  if (allocated(matrix2scc))   deallocate(matrix2scc)
  if (allocated(matrix2sc))    deallocate(matrix2sc)
  if (allocated(matrix2scm))   deallocate(matrix2scm)
  deallocate (factor,newleg,oldleg,newbin,oldbin,defp2bin,p2bin,     &
              defresbin,pspbin,mONS,cmONS2,cmREG2,lpmax,cfTot,csTot, &
              pCsTot,csIa,csIq,nIa,pIa,facIa,w0eTot,heli,dualheli,   &
              cd0sMax,w0Tot,w0last,parw0e,binw0e,legw0e,helw0e,      &
              w1tot,bm0prTot,bm1prTot,                               &
              zeroLO,c0EffMax,modaTot,colcoef,colcoefc,comp0gs,      &
              mosm0,binsm0,parsm0,xsm0,gsIncsm0,cosm0,gssm0,cssm0,   &
              dasd0,sesd0,facsd0,gssd0,cssd0,c0TOlp,bm0min,bm0max,   &
              bd0min,bd0max,sbm0,w0inbm0,w0outbm0,winitbm0,typebm0,  &
              sbd0,w0outbd0,winitbd0,gsTot,gs2Tot,gsTotEff,gs2TotEff)
  if (loopMax) deallocate (loopCoef,cEffMax,tiTot,ritiMax,comp1gs,   &
                           mosm1,binsm1,parsm1,gsIncsm1,cosm1,       &
                           rankInsm1,rankOutsm1,ferloopsm1,gssm1,    &
                           cssm1,tism1,dasd1,facsd1,rankOutsd1,      &
                           ferloopsd1,gssd1,cssd1,tisd1,legsti,      &
                           momsti,vmti,rankti,cTOt,cTOfh,cTOih1,     &
                           sbm1,w1inbm1,w0inbm1,       &
                           w1outbm1,winitbm1,typebm1,  &
                           sbd1,w1outbd1,winitbd1,w1TotMax,riwMax)
  if(allocated(nCache))    deallocate(nCache)
  if(allocated(nCacheTot)) deallocate(nCacheTot)
  if(allocated(tiCache))   deallocate(tiCache)
  if(allocated(CacheOn))   deallocate(CacheOn)

  if (allocated(bm1_b)) deallocate(bm1_b)
  if (allocated(bd1_b)) deallocate(bd1_b)
  ! tables
  deallocate (levelLeg,vectorLeg,firstNumber,firstGap,firstNumbers, &
              firstGaps,cm2n,als0R,Qren0R,Nlq0R)
  deallocate (alphai,eRefFactor,eFactor,eLOfactor)
  if (loopMax) deallocate (RtoS,riMin,riMax,ri,RItoR,RItoI,incRI,   &
                           firstRI)
  if (allocated(dZgs0R))     deallocate(dZgs0R)
  if (allocated(deltarpr))   deallocate(deltarpr)
  if (allocated(dalZpr))     deallocate(dalZpr)
  if (allocated(dalMSpr))    deallocate(dalMSpr)
  if (allocated(dZaapr))     deallocate(dZaapr)
  if (allocated(dZepr))      deallocate(dZepr)
  if (allocated(DdZe))       deallocate(DdZe)
  if (allocated(dgpnpr))     deallocate(dgpnpr)
  if (allocated(dgplpr))     deallocate(dgplpr)
  if (allocated(dgpupr))     deallocate(dgpupr)
  if (allocated(dgpdpr))     deallocate(dgpdpr)
  if (allocated(dgmnpr))     deallocate(dgmnpr)
  if (allocated(dgmlpr))     deallocate(dgmlpr)
  if (allocated(dgmupr))     deallocate(dgmupr)
  if (allocated(dgmdpr))     deallocate(dgmdpr)
  if (allocated(eNLOfactor)) deallocate(eNLOfactor)

! Set internal variables to default values

  loopQED  = .true.
  loopWEAK = .true.
  pureQED = .false.
  reguScheme = 2
  check_Pole = .false.
  resPar = .false.
  resIR = .false.
  longitudinal = 0
  warnings = 0

  prTot = 0
  loopMax = .false.

  ! Set als, Qren, Nfren and Nlq to their initialization values
  als   = als0
  Qren  = Qren0
  Nfren = Nfren0
  Nlq   = Nlq0
  use_active_qmasses = .false.

  ! Set Qren_alMS and Nlq_alMS to their initialization values
  Qren_alMS  = Qren0_alMS
   Nlq_alMS  =  Nlq0_alMS

  ifail = 0

  do i = 1,nOpened
    inquire(file=trim(nameOpened(i)),opened=fileopen,number=n)
    if (fileopen) close(n)
  enddo
  outputfile = 'output.rcl'
  nx = 934758
!  nOpened = 0
!  nameOpened(nOpenedDef) = ''

  collier_output_dir = 'default'

  processes_generated = .false.
  changed_lambda      = .false.
  changed_DeltaUV     = .false.
  changed_muUV        = .false.
  changed_DeltaIR     = .false.
  changed_muIR        = .false.

  end subroutine reset_recola_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module reset_rcl

!#####################################################################



