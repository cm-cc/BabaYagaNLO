! AAFURAMBO.  A NEW MONTE CARLO TREATMENT OF MULTIPARTICLE PHASE SPACE AT
! 1 HIGH ENERGIES.  R. KLEISS, W.J. STIRLING, S.D. ELLIS.
! REF. IN COMP. PHYS. COMMUN. 40 (1986) 359
    subroutine rambo(N,Et,Xm,P,Wt,Lw)
!------------------------------------------------------
!
!               RAMBO
!
!         RA(NdoM)  M(OMENTA)  BO(OSTER)
!
!   A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR
!   AUTHORS:  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING
!
!   N  = NUMBER OF PARTICLES (>1, IN THIS VERSION <101)
!   ET = TOTAL CENTRE-OF-MASS ENERGY
!   XM = PARTICLE MASSES ( DIM=N )
!   P  = PARTICLE MOMENTA ( DIM=(4,N) )
!   WT = WEIGHT OF THE EVENT
!   LW = FLAG FOR EVENT WEIGHTING:
!      LW = 0 WEIGHTED EVENTS
!      LW = 1 UNWEIGHTED EVENTS ( FLAT PHASE SPACE )
!------------------------------------------------------
    implicit none
    real(kind=8) :: a,accu,bq,c,Et,f,f0,g,g0,pm2,rmas,RN
    real(kind=8) :: s,sm2,w,Wt,wt2,wt3,wtm,wtmax,x,x2,xmax,xmt
    real(kind=8) :: Xm(N),P(4,N),q(4,100),xm2(100),r(4),b(3),e(100), &
               v(100),p2(100)

    integer    :: i,iter,k,Lw,N,nm

    real(kind=8), save :: acc=1d-14,z(100),twopi,po2log
    integer, save    :: itmax=6,iwarn(5)=0,ibegin=0

    ! INITIALIZATION STEP: FACTORIALS FOR THE PHASE SPACE WEIGHT
    if (ibegin .eq. 0) then
      ibegin = 1
      twopi = 8.*DATAN(1.D0)
      po2log = DLOG(twopi/4.)
      z(2) = po2log
      do k = 3 , 100
        z(k) = z(k-1) + po2log - 2.*DLOG(DFLOAT(k-2))
      end do
      do k = 3 , 100
        z(k) = (z(k)-DLOG(DFLOAT(k-1)))
      end do
    end if

    ! CHECK ON THE NUMBER OF PARTICLES
    if (N .gt. 1 .and. N .lt. 101) then

      ! CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES
      xmt = 0.
      nm = 0
      do i = 1 , N
        if (Xm(i) .ne. 0d0) nm = nm + 1
        xmt = xmt + DABS(Xm(i))
      end do
      if (xmt .gt. Et) then
        print 99001 , xmt , Et
99001     format (' RAMBO FAILS: TOTAL MASS =',D15.6,' IS NOT',     &
             &' SMALLER THAN TOTAL ENERGY =',D15.6)
        stop
      ! CHECK ON THE WEIGHTING OPTION
      else if (Lw .eq. 1 .or. Lw .eq. 0) then
        ! THE PARAMETER VALUES ARE NOW ACCEPTED
        ! GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE
 20     do i = 1 , N
          c = 2.*RN(1) - 1.
          s = DSQRT(1.-c*c)
          f = twopi*RN(2)
          q(4,i) = -DLOG(RN(3)*RN(4))
          q(3,i) = q(4,i)*c
          q(2,i) = q(4,i)*s*DCOS(f)
          q(1,i) = q(4,i)*s*DSIN(f)
        end do
        ! CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSformatION
        do i = 1 , 4
          r(i) = 0.
        end do
        do i = 1 , N
          do k = 1 , 4
            r(k) = r(k) + q(k,i)
          end do
        end do
        rmas = DSQRT(r(4)**2-r(3)**2-r(2)**2-r(1)**2)
        do k = 1 , 3
          b(k) = -r(k)/rmas
        end do
        g = r(4)/rmas
        a = 1./(1.+g)
        x = Et/rmas
        ! TRANSFORM THE Q'S CONFORMALLY INTO THE P'S
        do i = 1 , N
          bq = b(1)*q(1,i) + b(2)*q(2,i) + b(3)*q(3,i)
          do k = 1 , 3
            P(k,i) = x*(q(k,i)+b(k)*(q(4,i)+a*bq))
          end do
          P(4,i) = x*(g*q(4,i)+bq)
        end do
        ! return FOR UNWEIGHTED MASSLESS MOMENTA
        Wt = 1.D0
        if (nm .eq. 0 .and. Lw .eq. 1) return
        ! CALCULATE WEIGHT AND POSSIBLE WARNINGS
        Wt = po2log
        if (N .ne. 2) Wt = (2.*N-4.)*DLOG(Et) + z(N)
        if (Wt .lt. -180d0) then
          if (iwarn(1) .le. 5) print 99006 , Wt
          iwarn(1) = iwarn(1) + 1
        end if
        if (Wt .gt. 174d0) then
          if (iwarn(2) .le. 5) print 99007 , Wt
          iwarn(2) = iwarn(2) + 1
        end if
        ! return FOR WEIGHTED MASSLESS MOMENTA
        if (nm .ne. 0) then
          ! MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X
          xmax = DSQRT(1.-(xmt/Et)**2)
          do i = 1 , N
            xm2(i) = Xm(i)**2
            p2(i) = P(4,i)**2
          end do
          iter = 0
          x = xmax
          accu = Et*acc
 30       f0 = -Et
          g0 = 0.
          x2 = x*x
          do i = 1 , N
            e(i) = DSQRT(xm2(i)+x2*p2(i))
            f0 = f0 + e(i)
            g0 = g0 + p2(i)/e(i)
          end do
          if (DABS(f0) .gt. accu) then
            iter = iter + 1
            if (iter .le. itmax) then
              x = x - f0/(x*g0)
              GOTO 30
            else
              print 99002 , itmax
99002           format (' RAMBO WARNS:',I3,               &
                   &' ITERATIONS DID NOT GIVE THE',       &
                   &' DESIRED ACCURACY =',D15.6)
            end if
          end if
          do i = 1 , N
            v(i) = x*P(4,i)
            do k = 1 , 3
              P(k,i) = x*P(k,i)
            end do
            P(4,i) = e(i)
          end do
          ! CALCULATE THE MASS-EFFECT WEIGHT FACTOR
          wt2 = 1.
          wt3 = 0.
          do i = 1 , N
            wt2 = wt2*v(i)/e(i)
            wt3 = wt3 + v(i)**2/e(i)
          end do
          wtm = (2.*N-3.)*DLOG(x) + DLOG(wt2/wt3*Et)
          if (Lw .eq. 1) then
            ! UNWEIGHTED MASSIVE MOMENTA REQUIRED: ESTIMATE MAXIMUM WEIGHT
            Wt = DEXP(wtm)
            if (nm .le. 1) then
              ! ONE MASSIVE PARTICLE
              wtmax = xmax**(4*N-6)
            elseif (nm .gt. 2) then
              ! MORE THAN TWO MASSIVE PARTICLES: AN ESTIMATE ONLY
              wtmax = xmax**(2*N-5+nm)
            else
              ! TWO MASSIVE PARTICLES
              sm2 = 0.
              pm2 = 0.
              do i = 1 , N
                if (Xm(i) .ne. 0.d0) then
                  sm2 = sm2 + xm2(i)
                  pm2 = pm2*xm2(i)
                end if
              end do
              wtmax = ((1.-sm2/(Et**2))**2-4.*pm2/Et**4)**(N-1.5)
            end if
            ! DETERMINE WHETHER OR NOT TO ACCEPT THIS EVENT
            w = Wt/wtmax
            if (w .gt. 1.D0) then
              if (iwarn(5) .le. 5) print 99003 , wtmax , w
99003           format (                              &
                  &' RAMBO WARNS: ESTIMATE FOR MAXIMUM WEIGHT ='&
                 & ,D15.6,'    EXCEEDED BY A FACTOR ',D15.6)
              iwarn(5) = iwarn(5) + 1
            end if
            if (w .lt. RN(5)) GOTO 20
            Wt = 1.D0
            return
          else
            ! return FOR  WEIGHTED MASSIVE MOMENTA
            Wt = Wt + wtm
            if (Wt .lt. -180d0) then
              if (iwarn(3) .le. 5) print 99006 , Wt
              iwarn(3) = iwarn(3) + 1
            end if
            if (Wt .gt. 174.D0) then
              if (iwarn(4) .le. 5) print 99007 , Wt
              iwarn(4) = iwarn(4) + 1
            end if
            Wt = DEXP(Wt)
            return
          end if
        else
          Wt = DEXP(Wt)
          return
        end if
      else
        print 99004 , Lw
99004     format (' RAMBO FAILS: LW=',I3,' IS NOT AN ALLOWED OPTION')
        stop
      end if
    else
      print 99005 , N
99005   format (' RAMBO FAILS: # OF PARTICLES =',I5,' IS NOT ALLOWED')
      stop
    end if
99006 format (' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY UNDERFLOW')
99007 format (' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY  OVERFLOW')
    end

    double precision function RN(idmy)
      implicit none
      integer idmy
      real RVEC(1)
      ! idmy = idmy  ! causes a segfault on gfortran 15.X
      call RANLUX (RVEC,1)
      RN = DBLE(RVEC(1))
      return
    end
