!#####################################################################
!!
!!  File  draw_current_rcl.f90
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

module draw_current_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl,   only: draw
  use globals_rcl, only: list_vertices

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer :: gr2Max=0,gr3Max=0,gr4Max=0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine texHead(pr)

  integer, intent(in) :: pr
  logical :: is_first_pr
  character(len=4)  :: cpr
  character(len=30) :: file

  is_first_pr = pr.eq.1

  if (draw.ge.1) then

    write(cpr,'(i4)') pr
    file = 'process_'//trim(adjustl(cpr))//'.tex'
    open (unit=99,file=file,status='replace')

    write(99,'(a)') '\documentclass[11pt]{article}'
    write(99,'(a)')
    write(99,'(a)') '\usepackage{axodraw4j}'
    write(99,'(a)')
    write(99,'(a)') '\oddsidemargin -20pt \evensidemargin -20pt'
    write(99,'(a)') '\topmargin -40pt \headheight 00pt \headsep 00pt'
    write(99,'(a)') '\textheight 270mm \textwidth 172mm'
    write(99,'(a)')
    write(99,'(a)') '\begin{document}'
    write(99,'(a)')

  endif

  if (list_vertices.and.is_first_pr) then
    open (unit=2000,file='vertices2_temp.txt',status='replace')
    open (unit=2001,file='vertices2.txt',status='replace')
    open (unit=3000,file='vertices3_temp.txt',status='replace')
    open (unit=3001,file='vertices3.txt',status='replace')
    open (unit=4000,file='vertices4_temp.txt',status='replace')
    open (unit=4001,file='vertices4.txt',status='replace')
  endif

  end subroutine texHead

!---------------------------------------------------------------------

  subroutine texFeet(is_last_pr)
  logical, intent(in) :: is_last_pr
  integer       :: istatus,cty,gr,gr0,grMax,u0,u1
  character(19) :: file
  character(30) :: line

  if (draw.ge.1) then

    write(99,'(a)')
    write(99,'(a)') '\end{document}'

    close (unit=99,status='keep')

  endif

  if (list_vertices.and.is_last_pr) then

    do cty = 2,4
      select case (cty)
      case (2); grMAx=gr2Max; u0=2000; u1=2001; file='vertices2_temp.txt'
      case (3); grMAx=gr3Max; u0=3000; u1=3001; file='vertices3_temp.txt'
      case (4); grMAx=gr4Max; u0=4000; u1=4001; file='vertices4_temp.txt'
      end select
      do gr = 1,grMax
        close (unit=u0,status='keep')
        open  (unit=u0,file=file, &
               status='old',action='read',iostat=istatus)
        do
          read (u0,'(i4,a)',iostat=istatus) gr0,line
          if (istatus.ne.0) exit
          if (gr0.eq.gr) then
            write(u1,'(i4,a)') gr0,trim(line)
          endif
        enddo
        write(u1,'(/)')
      enddo
    enddo

    close (unit=2000,status='delete')
    close (unit=2001,status='keep')
    close (unit=3000,status='delete')
    close (unit=3001,status='keep')
    close (unit=4000,status='delete')
    close (unit=4001,status='keep')

  endif

  end subroutine texFeet

!---------------------------------------------------------------------

  subroutine picture(pr, lp, legs, legsE, e1, e2, e3, e4, &
                     p1, p2, p3, p4, t, park,         &
                     cdsMax, j, wT, csT, xs, newcut, wr)

  use skeleton_rcl, only: hosT, hmaT, binT, gsIncT
  use tables_rcl,   only: nmasses, levelLeg, nmf, firstNumbers
  use globals_rcl,  only: dp, anti

  integer, intent(in)    :: pr, lp, legs, legsE, cdsMax
  integer, intent(in)    :: e1,e2,e3,e4,p1,p2,p3,p4,t
  integer, intent(in)    :: j
  integer, intent(in)    :: csT(-1:,:), park(:), wT(:,:)
  logical, intent(in)    :: xs,newcut
  logical, intent(inout) :: wr(43,0:43,0:43,43)

  integer             :: bin(1:4),par(1:4),parbin(2**(legsE-1)),hm,  &
                         ho(1:size(hosT,1)),gs,e0,p0,height,yoffset, &
                         nprop,mass(legs),row,rowb,rowbl0,rowbl,     &
                         rowbIn(3),rshift,cshift,eProp0,eProp,ee,pp, &
                         l3,column,f1,f2,f3,f4,    &
                         gr2,gr3,gr4,w1,hm0,i,i1,ii,k,kk,l1,l2,n,s,step

  integer,allocatable :: iq(:,:),ia(:,:),iq0(:),ia0(:),gr(:,:,:,:)
  logical             :: axoinv(43)
  real (dp)           :: lenght
  character(len=55)   :: texpar(43)
  character(len=14)   :: axopar(43)
  character(len=9)    :: axoext
  character(len=2)    :: cl
  character(len=6)    :: texmass
  character(len=9)    :: axomass
  character(len=3)    :: txtpar(43)
  character(len=4)    :: ciq,cia
  character(len=999)  :: delta(4)

  if (draw.ge.1) then

    ! Here the legend is written
    if (lp.eq.0.and.newcut) then

      write(99,'(a)') ''

      write(99,'(a)') '\noindent'
      write(99,'(a)') '\makebox[\linewidth]{\rule{\textwidth}{0.4pt}}'
      write(99,'(a)') ''

      write(99,'(a)') '\vspace{.5cm}'
      write(99,'(a)') &
        'Here are listed the branches contributing to the given ', &
        'process.  '
      write(99,'(a)') ''

      write(99,'(a)') '\vspace{.3cm}'
      write(99,'(a)') &
        'On the left side of each figure there are the external ', &
        'currents contributing to the branch, whose binary ', &
        'number and particle are explictly given. ', &
        'Going from left to right, these currents enter the ', &
        'selected branch or a branch computed in a previous ', &
        'step. '
      write(99,'(a)') ''

      write(99,'(a)') '\vspace{.3cm}'
      write(99,'(a)') &
        'In loop branches, there is the first current of the ', &
        'cutted loop line, which is always at the top of the ', &
        'figure. ', &
        'Going downwards, this current enters the selected branch ', &
        'or a previously computed one. ', &
        'All previously computed loop branches are listed ', &
        'downwards and connected by loop propagators, whose mass ', &
        'is explicitly given ($m_0$ indicates a massless ', &
        'propagator). '
      write(99,'(a)') ''

      write(99,'(a)') '\vspace{.3cm}'
      write(99,'(a)') &
        'From the last previously computed branches, going from ', &
        'left to right there is the selected branch and then the ', &
        'out-going current. ', &
        'This one can be an off-shell current (big dot), the last ', &
        'external current (in the tree case) or the second ', &
        'current of the cutted loop line (for loop branches), ', &
        'which virtually close the loop. '
      write(99,'(a)') ''

      write(99,'(a)') '\vspace{.5cm}'
      write(99,'(a)') '\noindent'
      write(99,'(a)') '\begin{tabular}{c p{15.3cm}}'

      write(99,'(a)') '\begin{picture}(6,6)(-3,-3)'
      write(99,'(a)') '\GBoxc(0,0)(6,6){1}'
      write(99,'(a)') '\end{picture}'
      write(99,'(a)') '&'
      write(99,'(a)') 'It indicates the selected branch.'

      write(99,'(a)') '\\[+.3cm]'

      write(99,'(a)') '\begin{picture}(6,6)(-3,-3)'
      write(99,'(a)') '\GBoxc(0,0)(6,6){1}'
      write(99,'(a)') '\Line(-3,-3)(3,3)\Line(-3,3)(3,-3)'
      write(99,'(a)') '\end{picture}'
      write(99,'(a)') '&'
      write(99,'(a)') 'It indicates the selected branch, with special ', &
                      'Feynman rules for a counterterm.'

      write(99,'(a)') '\\[+.3cm]'

      write(99,'(a)') '\begin{picture}(6,6)(-3,-3)'
      write(99,'(a)') '\GBoxc(0,0)(6,6){0}'
      write(99,'(a)') '\Line(-3,-3)(3,3)\Line(-3,3)(3,-3)'
      write(99,'(a)') '\end{picture}'
      write(99,'(a)') '&'
      write(99,'(a)') 'It indicates the selected branch, with special ', &
                      'Feynman rules for a rational term.'

      write(99,'(a)') '\\[+.3cm]'

      write(99,'(a)') '\begin{picture}(6,6)(-3,-3)'
      write(99,'(a)') '\SetWidth{1.5}'
      write(99,'(a)') '\Line(-3,-3)(3,3)\Line(-3,3)(3,-3)'
      write(99,'(a)') '\end{picture}'
      write(99,'(a)') '&'
      write(99,'(a)') 'It indicates the cutted loop line.'

      write(99,'(a)') '\\[+.3cm]'

      write(99,'(a)') '\begin{picture}(6,6)(-3,-3)'
      write(99,'(a)') '\GCirc(0,0){3}{.5}'
      write(99,'(a)') '\end{picture}'
      write(99,'(a)') '&'
      write(99,'(a)') 'It indicates a branch computed in some ', &
                      'previous step.'

      write(99,'(a)') '\\[+.3cm]'

      write(99,'(a)') '\begin{picture}(6,6)(-3,-3)'
      write(99,'(a)') '\GCirc(0,0){1.5}{0}'
      write(99,'(a)') '\end{picture}'
      write(99,'(a)') '&'
      write(99,'(a)') 'It indicates the out-going off-shell ', &
                      'current of the selected branch.'

      write(99,'(a)') '\\[+.3cm]'

      write(99,'(a)') '$g_s^{\#}$'
      write(99,'(a)') '&'
      write(99,'(a)') 'In counterterms/rational terms insertion, ', &
                      'it indicates the power of $g_s$ of the ', &
                      'selected branch.'

      write(99,'(a)') '\end{tabular}'

      write(99,'(a)') ''
      write(99,'(a)') '\vspace{.4cm}'
      write(99,'(a)') &
        'Below each branch the colour structures of the incoming ', &
        'and outgoing currents of the select branch are ', &
        'explicitly written. ', &
        'In front of the colour structures of the incoming ', &
        'currents the binary number of the current is written. ', &
        'If a current has no colour structure, it has been built ', &
        'from colourless external currents. ', &
        'The colour structures are product of $\delta$s. ', &
        'The off-shell currents (i.e. all except the last ones) ', &
        'of coloured particles have always an "open" part, made of ', &
        '$\delta$s containing an $\alpha$ and/or a $\beta$. ', &
        'The other indices of the $\delta$s are the binary numbers ', &
        'of the external currents, which express in compact ', &
        'notation how the colour/anticolour of the extenal currents ', &
        'are contracted. For example ', &
        '$\delta^{2}_{1} \delta^{1}_\beta \delta^\alpha_{2}$ ', &
        'stays for $\delta^{\alpha_2}_{\beta_1} ', &
        '\delta^{\alpha_1}_\beta \delta^\alpha_{\beta_2}$, ', &
        'where $\alpha_i,\beta_i$ go from 1 to 3. '

      write(99,'(a)') ''
      write(99,'(a)') ''
      write(99,'(a)') ''
      write(99,'(a)') ''
      write(99,'(a)') '\vspace{.5cm}'
      write(99,'(a)') '\noindent'
      write(99,'(a)') '\makebox[\linewidth]{\rule{\textwidth}{0.4pt}}'
      write(99,'(a)') ''

      write(99,'(a)') '\clearpage'
      write(99,'(a)') ''

    endif

    bin(1) = e1; par(1) = p1
    bin(2) = e2; par(2) = p2
    bin(3) = e3; par(3) = p3
    bin(4) = e4; par(4) = p4

    w1 = wT(1,j)
    ho = hosT(:,w1)
    hm = hmaT(w1)

    gs = gsIncT(j)

    do i = 1,legsE
      parbin(2**(i-1)) = park(i)
    enddo

    ! define latex symbols for pariticles and their line
    texpar( 1) = '{\scriptstyle{\rm Y}_{\rm g}}'
    texpar( 2) = '{\scriptstyle{\rm X}_\gamma}'
    texpar( 3) = '{\scriptstyle{\rm X}_{\rm Z}}'
    texpar( 4) = '{\scriptstyle{\rm X}^+}'
    texpar( 5) = '{\scriptstyle{\rm X}^-}'
    texpar( 6) = '{\scriptstyle{\bar{\rm Y}}_{\rm g}}'
    texpar( 7) = '{\scriptstyle{\bar{\rm X}}_\gamma}'
    texpar( 8) = '{\scriptstyle{\bar{\rm X}}_{\rm Z}}'
    texpar( 9) = '{\scriptstyle{\bar{\rm X}}^+}'
    texpar(10) = '{\scriptstyle{\bar{\rm X}}^-}'
    texpar(11) = '{\scriptstyle{\rm H}}'
    texpar(12) = '\phi_{_0}'
    texpar(13) = '\phi^+'
    texpar(14) = '\phi^-'
    texpar(15) = '{\rm g}'
    texpar(16) = '\gamma'
    texpar(17) = '{\scriptstyle{\rm Z}}'
    texpar(18) = '{\scriptstyle{\rm W}\!^+}'
    texpar(19) = '{\scriptstyle{\rm W}\!^-}'
    texpar(20) = '\nu_{\rm e}'
    texpar(21) = '\nu_\mu'
    texpar(22) = '\nu_\tau'
    texpar(23) = '{\rm u}'
    texpar(24) = '{\rm c}'
    texpar(25) = '{\rm t}'
    texpar(26) = '{\rm e}^-'
    texpar(27) = '\mu^-'
    texpar(28) = '\tau^-'
    texpar(29) = '{\rm d}'
    texpar(30) = '{\rm s}'
    texpar(31) = '{\rm b}'
    texpar(32) = '\bar{\nu}_{\rm e}'
    texpar(33) = '\bar{\nu}_\mu'
    texpar(34) = '\bar{\nu}_\tau'
    texpar(35) = '\bar{\rm u}'
    texpar(36) = '\bar{\rm c}'
    texpar(37) = '\bar{\rm t}'
    texpar(38) = '{\rm e}^+'
    texpar(39) = '\mu^+'
    texpar(40) = '\tau^+'
    texpar(41) = '\bar{\rm d}'
    texpar(42) = '\bar{\rm s}'
    texpar(43) = '\bar{\rm b}'
    axopar(1:5)   = '\DashArrowLine'; axoinv(1:5)   = .false.
    axopar(6:10)  = '\DashArrowLine'; axoinv(6:10)  = .true.
    axopar(11:14) = '\DashLine'     ; axoinv(11:14) = .false.
    axopar(15)    = '\Gluon'        ; axoinv(15)    = .false.
    axopar(16:19) = '\Photon'       ; axoinv(16:19) = .false.
    axopar(20:31) = '\ArrowLine'    ; axoinv(20:31) = .false.
    axopar(32:43) = '\ArrowLine'    ; axoinv(32:43) = .true.

    ! Here the process is written
    if (lp.eq.0.and.newcut) then
      write(99,'(a)') '\vspace{.5cm}'
      write(99,'(a)') '\begin{center}'
      write(99,'(a,i2,a)') '\begin{tabular}{*{',2*legs,'}{c}}'
      write(99,'(a)') '\hline\hline \\[-.2cm]'
      write(99,'(a)') '\quad'
      do i = 1,legs
        if (i.ne.1) write(99,'(a)') '& $+$ &'
        write(99,'(3a)') '$',trim(texpar(park(i))),'$'
      enddo
      write(99,'(a)') '& $\,\to\quad 0 \quad$\\'
      write(99,'(a)') '\quad'
      do i = 1,legs; e0 = 2**(i-1)
        if (i.ne.1) write(99,'(a)') '& &'
        write(99,'(i3)') e0
      enddo
      write(99,'(a)') '&'
      write(99,'(a)') '\\[+.2cm]\hline\hline'
      write(99,'(a)') '\end{tabular}'
      write(99,'(a)') '\end{center}'
    endif

    if (newcut) then
      write(99,'(a)') ''
      write(99,'(a)') '\vspace{.5cm}'
      write(99,'(a)') '\begin{tabular}{|c|}'
      write(99,'(a)') '\hline\\[-.2cm]'
      select case (lp)
      case (0)
        write(99,'(a)') 'Tree level\\[.2cm]'
      case (1)
        write(99,'(3a)') &
          'Bare loop; \qquad cut-particle: $',trim(texpar(t)),'$\\[.2cm]'
      case (2)
        write(99,'(a)') 'Counterterms\\[.2cm]'
      case (3)
        write(99,'(a)') 'Rational terms\\[.2cm]'
      end select
      write(99,'(a)') '\hline'
      write(99,'(a)') '\end{tabular}'
      write(99,'(a)') '\\[.4cm]'
    endif

    if (e4.lt.2**legs) then
      height = 20*levelLeg(e4) + 5
      yoffset = - 20*levelLeg(e4)
    else
      height = 20*levelLeg(e4) + 10
      yoffset = - 20*levelLeg(e4)
    endif

    if (draw.lt.2) then
      write(99,'(a,i4,a,i4,a)') &
        '\begin{picture}(120,',height,')(-60,',yoffset,')'
    else
      write(99,'(a,i4,a,i4,a)') &
        '\begin{picture}(120,',height+70,')(-60,',yoffset-70,')'
    endif

    l1 = levelLeg(e1)
    l2 = levelLeg(e2)
    l3 = levelLeg(e3)

    if (e4.lt.2**legs) then ! tree branch

      ! The row of the last drawn external current (none at the
      ! moment) is initialized to 20 (later it is increased by -20
      ! each time an external current appears).
      row = 20

      ! Define row of the selected branch (rowb) and of the
      ! previously computed branches (rowbIn(1:3)), if present.
      rowbIn(1) = - 10 * (   l1 - 1 )
      rowbIn(2) = - 10 * ( 2*l1 +   l2 - 1 )
      rowbIn(3) = - 10 * ( 2*l1 + 2*l2 + l3 - 1 )
      if     (e2.eq.0) then; rowb  = rowbIn(1)
      elseif (e3.eq.0) then; rowb  = (rowbIn(1) + rowbIn(2))/2
      else;                  rowb  = rowbIn(2)
      endif

    else ! loop branch

      ! The row of the last drawn external current (none at the
      ! moment) is initialized to 0, because of the presence of the
      ! first loop current. "row" is then increased by -20 each time
      ! an external current appears.
      row  = 0

      ! Write binary number and particle of the first loop current.
      write(99,'(a,i4,a,a,a)') &
        '\Text(3,9)[cc]{\scriptsize $',2**legs, &
        '\;\;',trim(texpar(t)),'$}'

      nprop = maxval(ho)+1
      rowbl = 0
      if (ho(1).gt.0) then
        ! Compute the mass number of each propagator of the history
        hm0 = hm
        step = nmasses + 1
        do ii = legs,1,-1
          mass(ii) = hm0/step**(ii-1)
          hm0 = hm0 - mass(ii)*step**(ii-1)
        enddo
        ! Loop over propagators starting from the second one. When
        ! a propagator is selected, the line and properties of the
        ! previous one is written. This means:
        ! When the 2nd propagator is selected, the 1st is written;
        ! when the 3rd is selected, the 2nd is written; etc. ...
        do i = 2,nprop
          eProp0 = 0
          eProp  = 0
          do k = 1,size(ho)
            if (ho(k).gt.0) then
              if (ho(k)+1.eq.i-1) eProp0 = eProp0 + 2**(k-1)
              if (ho(k)+1.eq.i)   eProp  = eProp  + 2**(k-1)
            endif
          enddo
          rowbl0 = rowbl
          rowbl = row - 10 * ( levelLeg(eProp) + 1 )
          if (i.eq.2) then
            ! Here the line of the first loop current is written, if
            ! the second propagator is selected.
            write(cl,'(i2)') int(abs(rowbl)/4)
            select case (axopar(t))
            case ('\DashArrowLine','\DashLine'); axoext = '{3}'
            case ('\Gluon','\Photon'); axoext = '{1.5}{'//cl//'}'
            case default; axoext = ''
            end select
            if (axoinv(t)) then
              write(99,'(2a,i4,2a)')   &
                trim(axopar(t)),'(0,',rowbl,')(0,0)',axoext
            else
              write(99,'(2a,i4,2a)')   &
                trim(axopar(t)),'(0,0)(0,',rowbl,')',axoext
            endif
          else
            ! Here the mass for the previous propagator is defined
            if     (mass(i-1).eq.1)       then; texmass = 'm_0'
            elseif (mass(i-1).eq.2)       then; texmass = 'm_{_Z}'
            elseif (mass(i-1).eq.3)       then; texmass = 'm_{_W}'
            elseif (mass(i-1).eq.4)       then; texmass = 'm_{_H}'
            elseif (mass(i-1).eq.nmf(23)) then; texmass = 'm_u'
            elseif (mass(i-1).eq.nmf(24)) then; texmass = 'm_c'
            elseif (mass(i-1).eq.nmf(25)) then; texmass = 'm_t'
            elseif (mass(i-1).eq.nmf(26)) then; texmass = 'm_e'
            elseif (mass(i-1).eq.nmf(27)) then; texmass = 'm_\mu'
            elseif (mass(i-1).eq.nmf(28)) then; texmass = 'm_\tau'
            elseif (mass(i-1).eq.nmf(29)) then; texmass = 'm_d'
            elseif (mass(i-1).eq.nmf(30)) then; texmass = 'm_s'
            elseif (mass(i-1).eq.nmf(31)) then; texmass = 'm_b'
            endif
            ! Here the line for the previous propagator is written
            if (mass(i-1).eq.1) then ! massless (m_0)
              axomass = '\DashLine'; axoext = '{1}'
            else                         ! massive
              axomass = '\Line'; axoext = ''
            endif
            ! Here the line and the mass for the previous propagator
            ! is written
            write(99,'(2a,i4,a,i4,3a,i4,3a)') &
              trim(axomass),'(0,',rowbl0,')(0,',rowbl,')', &
              trim(axoext),'\Text(15,',(rowbl+rowbl0)/2, &
              ')[cc]{\scriptsize $',trim(texmass),'$}'
          endif
          ! Loop over the external currents which contributed to the
          ! (previously computed) loop branch, from which the
          ! selected propagator is coming.
          do i1 = 1,levelLeg(eProp)
            ! Binary of the external current is defined
            ee = firstNumbers(eProp,i1); pp = parbin(ee)
            ! "row" is increased by -20
            row = row - 20
            ! Here are written the binary number and the particle of
            ! the external currents which contributed to the loop
            ! branch, from which the selected propagator is coming.
            write(99,'(a,i4,a,i4,a,a,i4,3a)') &
              '\Text(-45,',row,')[rc]{\scriptsize',ee,'}', &
              '\Text(-35,',row,')[lc]{\scriptsize $', &
              trim(texpar(pp)),'$}'
            ! Here are written the lines of the external currents
            ! which contributed to the loop branch, from which the
            ! selected propagator is coming.
            lenght = sqrt( (20**2+(rowbl-row)**2)*1d0 )
            write(cl,'(i2)') int(lenght/4)
            select case (axopar(pp))
            case ('\DashArrowLine','\DashLine'); axoext = '{3}'
            case ('\Gluon','\Photon'); axoext = '{1.5}{'//cl//'}'
            case default; axoext = ''
            end select
            if (axoinv(pp)) then
              write(99,'(2a,i4,a,i4,2a)')                &
                trim(axopar(pp)),'(0,',rowbl,')(-20,',row,')',axoext
            else
              write(99,'(2a,i4,a,i4,2a)')                &
                trim(axopar(pp)),'(-20,',row,')(0,',rowbl,')',axoext
            endif
          enddo
          ! Here is drawn the symbol (grey circle) of the branch from
          ! which the previous propagator is coming.
          if (i.gt.2) &
            write(99,'(a,i4,a)') '\GCirc(0,',rowbl0,'){3}{.5}'
        enddo
      endif
      ! Define row of the selected branch (rowb) and of the
      ! previously computed branches (rowbIn(1:3)), if present.
      rowbIn(1) = rowbl
      rowbIn(2) = row - 10 * ( l2 + 1 )
      rowbIn(3) = row - 10 * ( 2*l2 + l3 + 1 )
      if (e3.eq.0) then; rowb  = (rowbl + rowbIn(2))/2
      else;              rowb  = rowbIn(2)
      endif

    endif

    ! Loop over the in-coming currents of the selected branch
    do i = 1,3
      e0 = bin(i); p0 = par(i)
      if (e0.eq.0) cycle
      if (e0.lt.2**legs) then
        ! Loop over the external currents contributing to the
        ! in-coming currents of the selected branch (not done for the
        ! loop current).
        do i1 = 1,levelLeg(e0)
          ! Binary of the external current is defined
          ee = firstNumbers(e0,i1); pp = parbin(ee)
          ! "row" is increased by -20
          row = row - 20
          ! Here are written the binary number and the particle of
          ! the external currents contributing to the in-coming
          ! currents of the selected branch.
          write(99,'(a,i4,a,i4,a,a,i4,3a)') &
            '\Text(-45,',row,')[rc]{\scriptsize',ee,'}', &
            '\Text(-35,',row,')[lc]{\scriptsize $',trim(texpar(pp)),'$}'
          if (levelLeg(e0).gt.1) then
            ! Here are written the lines of the external currents
            ! contributing to the in-coming currents of the selected
            ! branch.
            lenght = sqrt( (20**2+(rowbIn(i)-row)**2)*1d0 )
            write(cl,'(i2)') int(lenght/4)
            select case (axopar(pp))
            case ('\DashArrowLine','\DashLine'); axoext = '{3}'
            case ('\Gluon','\Photon'); axoext = '{1.5}{'//cl//'}'
            case default; axoext = ''
            end select
            if (axoinv(pp)) then
              write(99,'(2a,i4,a,i4,2a)') &
                trim(axopar(pp)),'(0,',rowbIn(i),')(-20,',row,')',axoext
            else
              write(99,'(2a,i4,a,i4,2a)') &
                trim(axopar(pp)),'(-20,',row,')(0,',rowbIn(i),')',axoext
            endif
          endif
        enddo
      endif
      if (levelLeg(e0).eq.1) then
        ! Here are written the lines of the in-coming currents of the
        ! selected branch, in the case they are external currents.
        if (e0.lt.2**legs) then; column = - 20
        else;                    column =    0
        endif
        lenght = sqrt( ((20-column)**2+(rowb-row)**2)*1d0 )
        write(cl,'(i2)') int(lenght/4)
        select case (axopar(p0))
        case ('\DashArrowLine','\DashLine'); axoext = '{3}'
        case ('\Gluon','\Photon'); axoext = '{1.5}{'//cl//'}'
        case default; axoext = ''
        end select
        if (axoinv(p0)) then
          write(99,'(2a,i4,a,i4,a,i4,2a)') &
            trim(axopar(p0)),'(20,',rowb,')(',column,',',row,')',axoext
        else
          write(99,'(2a,i4,a,i4,a,i4,2a)') &
            trim(axopar(p0)),'(',column,',',row,')(20,',rowb,')',axoext
        endif
      else
        ! Here are written the lines of the in-coming currents of the
        ! selected branch, in the case they are not external currents.
        lenght = sqrt( (20**2+(rowbIn(i)-rowb)**2)*1d0 )
        write(cl,'(i2)') int(lenght/4)
        select case (axopar(p0))
        case ('\DashArrowLine','\DashLine'); axoext = '{3}'
        case ('\Gluon','\Photon'); axoext = '{1.5}{'//cl//'}'
        case default; axoext = ''
        end select
        if (axoinv(p0)) then
          write(99,'(2a,i4,a,i4,2a)') &
            trim(axopar(p0)),'(20,',rowb,')(0,',rowbIn(i),')',axoext
        else
          write(99,'(2a,i4,a,i4,2a)') &
            trim(axopar(p0)),'(0,',rowbIn(i),')(20,',rowb,')',axoext
        endif
        ! Here is written the particle of the in-coming currents of
        ! the selected branch, in the case they are not external
        ! currents.
        if (i.eq.1) then
          if (e2.eq.0) then
            rshift = + 4
            cshift = -12
          else
            rshift = + 1
            cshift = - 6
          endif
        elseif (i.eq.2) then
          if (e3.eq.0) then
            rshift = - 8
            cshift = - 6
          else
            rshift = + 3
            cshift = -17
          endif
        elseif (i.eq.3) then
          rshift = - 8
          cshift = - 6
        endif
        write(99,'(a,i4,a,i4,3a)') &
          '\Text(',20+cshift,',',(rowbIn(i)+rowb)/2+rshift, &
          ')[lb]{\scriptsize $',trim(texpar(p0)),'$}'
        ! Here is drawn the symbol (grey circle) of the branch from
        ! which the in-coming currents of the selected branch are
        ! coming.
        write(99,'(a,i4,a)') '\GCirc(0,',rowbIn(i),'){3}{.5}'
      endif
    enddo

    ! Here are written the line of the out-going current of the
    ! selected branch.
    select case (axopar(p4))
    case ('\DashArrowLine','\DashLine'); axoext = '{3}'
    case ('\Gluon','\Photon'); axoext = '{1.5}{5}'
    case default; axoext = ''
    end select
    if (axoinv(p4)) then
      write(99,'(2a,i4,a,i4,2a)') &
        trim(axopar(p4)),'(40,',rowb,')(20,',rowb,')',axoext
    else
      write(99,'(2a,i4,a,i4,2a)') &
        trim(axopar(p4)),'(20,',rowb,')(40,',rowb,')',axoext
    endif
    ! Here are written the particle of the out-going current of the
    ! selected branch.
    write(99,'(a,i4,3a)') &
      '\Text(46,',rowb,')[lc]{\scriptsize $',trim(texpar(p4)),'$}'
    ! Here is drawn the symbol (white box) of the selected branch
    write(99,'(a,i4,a)') '\GBoxc(20,',rowb,')(6,6){1}'

    if (xs) then
      if (lp.eq.2) then
        ! Here is drawn the symbol (crossed white box) of the selected
        ! branch for counterterm insertion.
        write(99,'(a,i4,a,i4,a,i4,a,i4,a)')      &
          '\Line(17,',rowb-3,')(23,',rowb+3,     &
          ')\Line(17,',rowb+3,')(23,',rowb-3,')'
      elseif (lp.eq.3) then
        ! Here is drawn the symbol (black box) of the selected
        ! branch for rational term insertion.
        write(99,'(a,i4,a)') '\GBoxc(20,',rowb,')(6,6){0}'
      endif
      ! Here is written the power of g_s of the selected branch for
      ! counterterm/rational term insertion.
      write(99,'(a,i4,a,i1,a)') &
        '\Text(28,',rowb+5,')[cb]{\scriptsize $g_s^{\,',gs,'}$}'
    endif

    ! Here are written the colour structures of the currents of the
    ! selected branch
    if (draw.ge.2.and.cdsMax.gt.0) then
      allocate (iq(cdsMax,4)); iq = 0
      allocate (ia(cdsMax,4)); ia = 0
      do i = 1,4
        delta(i) = ''
        if (wT(i,j).eq.0) cycle
        allocate (iq0(cdsMax))
        allocate (ia0(cdsMax))
          ii = 0
          do kk = 1,legs
            if (csT(kk,wT(i,j)).ne.0) then
              ii = ii + 1
              iq0(ii) = kk
              ia0(ii) = csT(kk,wT(i,j))
            endif
          enddo
        iq(1:cdsMax,i) = 2**(iq0(1:cdsMax)-1)
        ia(1:cdsMax,i) = 2**(ia0(1:cdsMax)-1)
        deallocate (iq0)
        deallocate (ia0)
        do n = 1,cdsMax; if (iq(n,i)*ia(n,i).eq.0) exit
          write(ciq,'(i4)') iq(n,i)
          write(cia,'(i4)') ia(n,i)
          delta(i) = trim(adjustl(delta(i)))//'\delta_{j'//trim(adjustl(ciq))//'}'
          delta(i) = trim(adjustl(delta(i)))//'^{i'//trim(adjustl(cia))//'}'
        enddo
        if (csT(-1,wT(i,j)).ne.0) then
          write(ciq,'(i4)') 2**(csT(-1,wT(i,j))-1)
          delta(i) = trim(adjustl(delta(i)))//'\delta_{j'//trim(adjustl(ciq))//'}^i'
        endif
        if (csT(0,wT(i,j)).ne.0) then
          write(cia,'(i4)') 2**(csT(0,wT(i,j))-1)
          delta(i) = trim(adjustl(delta(i)))//'\delta_j^{i'//trim(adjustl(cia))//'}'
        endif
      enddo
      deallocate (iq)
      deallocate (ia)
      s = row - 14
      write(99,'(a,i4,a)') &
        '\Text(-50,',s,')[l]{\scriptsize Incoming Colour Structures:}'
      s = s + 1
      do i = 1,3; if (delta(i).eq.'') cycle
        s = s - 13
        write(99,'(a,i4,a,i4,a)') &
          '\Text(-47,',s,')[lb]{\scriptsize $',binT(wT(i,j)),'$}'
        s = s + 1
        write(99,'(a,i4,a,3a)') &
          '\Text(-30,',s,')[l]{\scriptsize $',trim(adjustl(delta(i))),'$}'
      enddo
      s = s - 12
      write(99,'(a,i4,a)') &
        '\Text(-50,',s,')[l]{\scriptsize Outgoing Colour Structure:}'
      s = s + 1
      if (delta(4).ne.'') then
        s = s - 13
!        write(99,'(a,i4,a,i4,a)') &
!          '\Text(-47,',s,')[lb]{\scriptsize $',binT(wT(4,j)),'$}'
        s = s + 1
        write(99,'(a,i4,3a)') &
          '\Text(-30,',s,')[l]{\scriptsize $',trim(adjustl(delta(4))),'$}'
      endif
    endif

    ! Write the cross for the first and last loop current.
    write(99,'(a)') '\SetWidth{1.5}'
    if (e4.ge.2**legs) then
      ! Write the cross for the first loop current.
      write(99,'(a)') '\Line(-3,-3)(3,3)\Line(-3,3)(3,-3)'
    endif
    if (e4.eq.2**(legs+1)-1) then
      ! Write the cross for the last loop current.
      write(99,'(a,i4,a,i4,a,i4,a,i4,a)') &
        '\Line(37,',rowb-3,')(43,',rowb+3, &
        ')\Line(37,',rowb+3,')(43,',rowb-3,')'
    elseif (e4.lt.2**(legsE-1)-1) then
      ! Write the big dot for the out-going current, if it is not the

      ! last external or loop current.
      write(99,'(a,i4,a)') '\GCirc(40,',rowb,'){1.5}{0}'
    endif

    write(99,'(a)') '\end{picture}'
    write(99,'(a)') '%'

  endif

  if (list_vertices) then

    txtpar( 1) = ' xg'
    txtpar( 2) = ' xA'
    txtpar( 3) = ' xZ'
    txtpar( 4) = ' x+'
    txtpar( 5) = ' x-'
    txtpar( 6) = ' Xg'
    txtpar( 7) = ' XA'
    txtpar( 8) = ' XZ'
    txtpar( 9) = ' X+'
    txtpar(10) = ' X-'
    txtpar(11) = ' H '
    txtpar(12) = ' p0'
    txtpar(13) = ' p+'
    txtpar(14) = ' p-'
    txtpar(15) = ' g '
    txtpar(16) = ' A '
    txtpar(17) = ' Z '
    txtpar(18) = ' W+'
    txtpar(19) = ' W-'
    txtpar(20) = ' ne'
    txtpar(21) = ' nm'
    txtpar(22) = ' nt'
    txtpar(23) = ' u '
    txtpar(24) = ' c '
    txtpar(25) = ' t '
    txtpar(26) = ' e-'
    txtpar(27) = ' m-'
    txtpar(28) = ' t-'
    txtpar(29) = ' d '
    txtpar(30) = ' s '
    txtpar(31) = ' b '
    txtpar(32) = ' Ne'
    txtpar(33) = ' Nm'
    txtpar(34) = ' Nt'
    txtpar(35) = ' U '
    txtpar(36) = ' C '
    txtpar(37) = ' T '
    txtpar(38) = ' e+'
    txtpar(39) = ' m+'
    txtpar(40) = ' t+'
    txtpar(41) = ' D '
    txtpar(42) = ' S '
    txtpar(43) = ' B '

    ! f1, f2, f3, f4 are incoming
    f1 = p1
    f2 = p2
    f3 = p3
    f4 = anti(p4)

    if (e4.ge.2**legs.and.wr(f1,f2,f3,f4)) then
      if (.not. allocated(gr)) allocate(gr(43,43,0:43,43))
      if (e2.eq.0) then

        if     (.not.wr(f4,f2,f3,f1)) then; gr2 = gr(f4,f2,f3,f1)
        else;          gr2Max = gr2Max + 1; gr2 = gr2Max
        endif
        gr(f1,f2,f3,f4) = gr2
        write(2000,'(i4,3x,2a,3x,a,i2,a)') &
          gr2,txtpar(f1),txtpar(f4),'(',pr,')'

      elseif (e3.eq.0) then

        if     (.not.wr(f1,f4,f3,f2)) then; gr3 = gr(f1,f4,f3,f2)
        elseif (.not.wr(f2,f1,f3,f4)) then; gr3 = gr(f2,f1,f3,f4)
        elseif (.not.wr(f2,f4,f3,f1)) then; gr3 = gr(f2,f4,f3,f1)
        elseif (.not.wr(f4,f1,f3,f2)) then; gr3 = gr(f4,f1,f3,f2)
        elseif (.not.wr(f4,f2,f3,f1)) then; gr3 = gr(f4,f2,f3,f1)
        else;          gr3Max = gr3Max + 1; gr3 = gr3Max
        endif
        gr(f1,f2,f3,f4) = gr3
        write(3000,'(i4,3x,3a,3x,a,i2,a)') &
          gr3,txtpar(f1),txtpar(f2),txtpar(4),'(',pr,')'

      else

        if     (.not.wr(f1,f2,f4,f3)) then; gr4 = gr(f1,f2,f4,f3)
        elseif (.not.wr(f1,f3,f2,f4)) then; gr4 = gr(f1,f3,f2,f4)
        elseif (.not.wr(f1,f3,f4,f2)) then; gr4 = gr(f1,f3,f4,f2)
        elseif (.not.wr(f1,f4,f2,f3)) then; gr4 = gr(f1,f4,f2,f3)
        elseif (.not.wr(f1,f4,f3,f2)) then; gr4 = gr(f1,f4,f3,f2)
        elseif (.not.wr(f2,f1,f3,f4)) then; gr4 = gr(f2,f1,f3,f4)
        elseif (.not.wr(f2,f1,f4,f3)) then; gr4 = gr(f2,f1,f4,f3)
        elseif (.not.wr(f2,f3,f1,f4)) then; gr4 = gr(f2,f3,f1,f4)
        elseif (.not.wr(f2,f3,f4,f1)) then; gr4 = gr(f2,f3,f4,f1)
        elseif (.not.wr(f2,f4,f1,f3)) then; gr4 = gr(f2,f4,f1,f3)
        elseif (.not.wr(f2,f4,f3,f1)) then; gr4 = gr(f2,f4,f3,f1)
        elseif (.not.wr(f3,f1,f2,f4)) then; gr4 = gr(f3,f1,f2,f4)
        elseif (.not.wr(f3,f1,f4,f2)) then; gr4 = gr(f3,f1,f4,f2)
        elseif (.not.wr(f3,f2,f1,f4)) then; gr4 = gr(f3,f2,f1,f4)
        elseif (.not.wr(f3,f2,f4,f1)) then; gr4 = gr(f3,f2,f4,f1)
        elseif (.not.wr(f3,f4,f1,f2)) then; gr4 = gr(f3,f4,f1,f2)
        elseif (.not.wr(f3,f4,f2,f1)) then; gr4 = gr(f3,f4,f2,f1)
        elseif (.not.wr(f4,f1,f2,f3)) then; gr4 = gr(f4,f1,f2,f3)
        elseif (.not.wr(f4,f1,f3,f2)) then; gr4 = gr(f4,f1,f3,f2)
        elseif (.not.wr(f4,f2,f1,f3)) then; gr4 = gr(f4,f2,f1,f3)
        elseif (.not.wr(f4,f2,f3,f1)) then; gr4 = gr(f4,f2,f3,f1)
        elseif (.not.wr(f4,f3,f1,f2)) then; gr4 = gr(f4,f3,f1,f2)
        elseif (.not.wr(f4,f3,f2,f1)) then; gr4 = gr(f4,f3,f2,f1)
        else;          gr4Max = gr4Max + 1; gr4 = gr4Max
        endif
        gr(f1,f2,f3,f4) = gr4
        write(4000,'(i4,3x,4a,3x,a,i2,a)') &
          gr4,txtpar(f1),txtpar(f2),txtpar(f3),txtpar(f4),'(',pr,')'

      endif

      wr(f1,f2,f3,f4) = .false.

    endif

  endif
    if (allocated(gr)) deallocate(gr)

  end subroutine picture

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module draw_current_rcl

!#####################################################################
