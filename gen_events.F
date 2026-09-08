      subroutine gen_events(E1,E2,sdifoutput,p1lab,p2lab,p1cut,
     >     p2cut,qphcut,ng)
#ifdef COLLIER
      use collier
#endif
      implicit double precision (a-h,o-z)
#include "commonmain.F"

      counterbefore=counterg4

ccc   call rluxut(isvec)

****  THE FOLLOWING LINES BETWEEN "*[[ start" AND "* end ]]" MUST
**** BE RECALCULATED FOR EACH EVENT IF ECMS CHANGES EVENT BY EVENT.
**** FURTHERMORE, IF ECMS CHANGES EVENT BY EVENT ifrstenbld MUST
**** BE SET TO 0
*[[ start
        ifrstenbld = 0

****** Beam energy spread no more calculated, taken from input E1 and E2
****** Assuming them to be on the z-axis
!        call spreadbeams(esig1,esig2,ame,am1,p1lab,p2lab,ieb)

        p1lab=0.d0
        p2lab=0.d0

        p1lab(0)=E1
        p1lab(3)=E1*sqrt((1.d0+ame/E1)*(1.d0-ame/E1))

        p2lab(0)=E2
        p2lab(3)=-E2*sqrt((1.d0+ame/E2)*(1.d0-ame/E2))

**************************************************************
        p12lab = p1lab + p2lab
        call new_boost(p12lab,p1lab,pin1,1)
        call new_boost(p12lab,p2lab,pin2,1)
**** to avoid round offs, I recalculate off-shell relations
        et = pin1(0)+pin2(0)
        ebt = 0.5d0*et
        pin1(0) = ebt
        pin1(3) = ebt*sqrt((1.d0+ame/ebt)*(1.d0-ame/ebt))
        pin2(0) =  pin1(0)
        pin2(3) = -pin1(3)
        ecms=sqrt(dot(pin1+pin2,pin1+pin2))
******
        esoft  = eps * ecms/2.d0 ! ecms is in a common, it is changed in spreadbeams

        pin(0) = ecms
        pin(2) = 0.d0
        pin(3) = 0.d0
        pin(1) = 0.d0
        betafs = sqrt(1.d0 - 4.d0*am1*am1/pin(0)**2)
        beta   = sqrt(1.d0 - 4.d0*ame*ame/pin(0)**2)
        do ki = 0,3
           pbeam1(ki) = pin1(ki)
           pbeam2(ki) = pin2(ki)
        enddo
        s = ecms*ecms
* end ]]

** [[end initialization]]

         k = k + 1
** for debugging!!

***

         flux = 8.d0 * (ecms/2.d0)**2

         wcos=1.d0
         cth=0.d0

         call multiplicity(eps,ecms,cth,ng,wnphot)

         if (fs.eq.'gg') then
            if (ng.eq.0) then
               call get_cos_2g_born(cth,wcos)
            else
               call get_cos_2g(cth,wcos)
            endif
         endif

         sdif = wnphot*wcos
***************************************************
** initial state radiation with structure functions
*
         if (ord.eq.'struct') then
            call emission(pbeam1,pbeam2,wisr,'e')

            betafs = sqrt(1.d0 - 4d0*am1**2/(pbeam1(0)+pbeam2(0))**2)
            beta   = sqrt(1.d0 - 4d0*ame**2/(pbeam1(0)+pbeam2(0))**2)
*
            sreduced =   (pbeam1(0)+pbeam2(0))**2
     +           - (pbeam1(3)+pbeam2(3))**2

            pbeam1p(0) = sqrt(sreduced)/2.d0
            pbeam1p(1) = 0d0
            pbeam1p(2) = 0d0
            pbeam1p(3) = beta*pbeam1p(0)
*
            pbeam2p(0) = pbeam1p(0)
            pbeam2p(1) = 0d0
            pbeam2p(2) = 0d0
            pbeam2p(3) = -pbeam1p(3)
*
            sdif = sdif * wisr

         else
            pbeam1p = pbeam1
            pbeam2p = pbeam2
         endif

*
** end of initial state radiation
*******************************************
         fsbck=fs
         if(fs.eq.'er') then

            fs='ee'
         endif
         if (fs.ne.'gg') then
            call phasespaceCRF(pbeam1p,pbeam2p,p1,p2,qph,ng,am1,am2
     .           ,esoft,cth,w,phsp,ie)
         else

            call phasespaceStrong2020(pbeam1p,pbeam2p,p1,p2,qph,ng,
     .           am1,am2,esoft,cth,w,phsp,ie)


         endif

         fs=fsbck

         ie = ie + ieb          ! ieb if beam energy spread fails
         if (ieb.gt.0) phsp = 0.d0
***

         sdif = sdif * phsp * w

         if (phsp /= phsp.or. w /= w) then
            print*,' '
            print*,'NAN found (pos 1)!!! k=',k
            print*,'outfile ',outfile
            print*,'ord,ng,w,phsp ',ord,ng,w,phsp
            print*,'end ----!!!'
            print*,'setting ie>0 and sdif to 0'
            print*,' '
            call printstatus(1,k,p1,p2,qph,xsec,var,varbefore,sdif,
     .           sdifmax,fmax)

            ie = 1
            w = 0.d0
            sdif = 0.d0
            phsp = 0.d0
         endif
*
** final state radiation with structure function
         if (ord.eq.'struct'.and.ie.eq.0) then
            do kk = 0,3
               p1cut(kk) = p1(kk)
               p2cut(kk) = p2(kk)
               pcm(kk)   = pbeam1(kk) + pbeam2(kk)
            enddo
            do ii = 1,ng
               do kk = 0,3
                  qcut(kk)  = qph(ii,kk)
               enddo
               call new_boost(pcm,qcut,qcut,-1)
               do kk = 0,3
                  qphcut(ii,kk)  = qcut(kk)
               enddo
            enddo

            call new_boost(pcm,p1cut,p1cut,-1)
            call new_boost(pcm,p2cut,p2cut,-1)

            if (fs.eq.'mm'.or.fs.eq.'mr') then
                call emission(p1cut,p2cut,wfsr,'m')
            endif
            if (fs.eq.'ee') call emission(p1cut,p2cut,wfsr,'e')
            if (fs.eq.'pp') call emission(p1cut,p2cut,wfsr,'p')
            if (fs.eq.'pi') then
               print*,'call emission with FSR'
               print*,'NON DEVO FINIRE QUI!!!!!!!!'
            endif
         endif
*     * end of final state radiation

*****************************************
         if (ie.ge.1) ie = 1
         if (fs.eq.'gg'.and.ie.eq.0) then
            call mixandsortmomenta(ng,p1,p2,qph)
         endif

         if (ord.ne.'struct'.and.ie.eq.0) then
            call new_boost(p12lab,p1,p1cut,-1)
            call new_boost(p12lab,p2,p2cut,-1)

            qphcut = 0.d0
            if (ng.gt.0) then
               do ii = 1,ng
                  ptmp = qph(ii,:)
                  call new_boost(p12lab,ptmp,ptmp,-1)
                  qphcut(ii,:) = ptmp
               enddo
            endif
         endif
         if (ie.lt.1) then
            if (fs.eq.'gg') then
               call cutsgg(ng,p1cut,p2cut,qphcut,icut)
            else
               call cuts(p1cut,p2cut,qphcut,icut,ng)
            endif
         else
            icut = 1
         endif
         ie = ie + icut

         if (ng.gt.nphmax.and.ie.eq.0) nphmax = ng
         if (icut.eq.0) naccepted = naccepted + 1
         if (icut.eq.0.and.ng.ge.3) naccepted3=naccepted3+1

         irepeat = irepeat + 1

         iswvperr  = 0
         iswvperrt = 0

	 epsbck=eps
         call squared_matrix(model,ng,ecms,p1,p2,pbeam1,pbeam2,qph,
     >        ie,icalc,emtx,prod,deltah)
         irepeat = 0
         bck  = emtx

         emtx = emtx/in_conf_spin ! divided by initial spin conf.
         emtx = emtx*convfac/flux ! divided by the flux

!!!!!!!!per i radiativi e pp (costruito on top del 2->4), analogamente al processo gammagamma aggiungo
!!!!!!!!un fattore di correzione statistico 1/ng
!!!!!!!!considerando che ng include gi  il fotone del radiativo sottostante
	 condfs=(fs.eq.'er'.or.fs.eq.'mr'.or.fs.eq.'pr'.or.fs.eq.'pp')
	 if(condfs.and.ng.ge.3) then

            emtx = emtx / (ng)  !/(ng + 1)
         endif
!!!!!!!!fine correzione dei radiativi

         if (fs.eq.'gg') then
            emtx = emtx/(ng+2)/(ng+1)
         endif

         sv = 1.d0

         if (ie.eq.0) then

***   MODIFY MARCO  FOR MAURO UNWEIGHTING

            modelbck=model
            if(condexp) then
               model=modelunw
            endif
            if(condnlo.and.k.le.nsearch) model='psrad'
            if(condnlo.and.k.gt.nsearch) model='matched'

            call getrndstatus(isvec)
            call svfactor(model,ng,ecms,p1,p2,eps,sv,deltasv,q,qph)

!     To avoid instabilities
            if(abs(deltasv).ge.2.d0) then
                 ltint=1
               fscoll =1
               iscoll =1
!     If dispersive pions we want to evaluate the same s',s'', so we reset isvec to the one before svfactor 
               call rndreset(isvec)
               call svfactor(model,ng,ecms,p1,p2,eps,
     >              sv,deltasv,q,qph)
                 ltint=0
               endif


            if((1d0+deltasv).ge.deltamax) deltamax=1.d0+deltasv

            if((1d0+deltasv).le.deltamin) deltamin=1.d0+deltasv

***********************THIS IS TO CHECK FOR EVENTS WHERE COLLIER CRASHES
#ifdef COLLIER
            call GetErrFlag_cll(errflag)
#endif

            if(errflag.eq.-10) then

#ifdef COLLIER
               call InitErrFlag_cll()
               call setmode_cll(1)
#endif

               call svfactor(model,ng,ecms,p1,p2,eps,sv,deltasv,q,qph)

#ifdef COLLIER
            call GetErrFlag_cll(errflag)
#endif

#ifdef COLLIER
               call setmode_cll(2)
#endif
            endif

****************

            model=modelbck

****************END MODIFICA

            sdif = sdif * sv
         else
            sdif = 0.d0
         endif

         sdif=sdif*emtx

c     -----------------------------    c
c     Radiative phi decays for pr      c
c     -----------------------------    c     
c     This correction is  only for     c
c             radiative pions          c
c     -----------------------------    c
         
!     Occhio al reweighting, deltasv ora è relativo al solo pezzo puro sQED e approssimazioni varie nel loop
         indexcorr = 100*indexbr+10*indexphires+indexdres
         sdifnophi = sdif
!         if(fs.eq.'pr'.and.indexcorr.ne.0.and.ie.eq.0) then
         if((fs.eq.'pp'.or.fs.eq.'pr')
     >        .and.ng.ge.1.and.indexcorr.ge.0.and.ie.eq.0) then
            svphi = 1.d0
!            print*,'passo'
            if(fs.eq.'pp'.and.ord.eq.'alpha') then
               modelphi='born'  !Il sample a 1gamma nel caso di ord=alpha per pp è un born
            else
               modelphi = 'psrad'
            endif
            call svfactor(modelphi,ng,ecms,p1,p2,eps,
     >           svphi,deltasv,q,qph)
            model=modelbck
            call squared_matrix_phi(model,ng,ecms,p1,p2,pbeam1,pbeam2,
     >       qph,ie,icalc,emtxphi,prod,deltah)
            weightpss = wnphot*wcos*phsp * w
            diffphi = emtxphi/in_conf_spin*convfac/flux*weightpss*svphi
            sdif = +sdif +  diffphi
         endif


         if (emtx /= emtx.or. sv /= sv) then
            print*,'NaN found here. Resetting sdif to 0.'
            print*,'emtx sv w phsp cth',emtx,sv,w,phsp,cth
            call printstatus(1,k,p1,p2,qph,xsec,var,varbefore,sdif,
     .           sdifmax,fmax)
            sdif = 0.d0
            phsp = 0.d0
            w    = 0.d0
            ie   = 1
         endif

         if (ord.eq.'struct'.and.ie.eq.0) then
            do kk = 0,3
                p1(kk) = p1cut(kk)
                p2(kk) = p2cut(kk)
                qph(1,kk) = qcut(kk)
            enddo
            sdif = sdif * wfsr
         endif

         if(fs.ne.'mr'.and.fs.ne.'er'
     >                .and.fs.ne.'pr') then
            sdif = sdif/bornme*(bornme + bornmez)
         endif

!     !         emtx = 1.d0 ! uncomment  for phase space integral

         iii = 0
         if (sdif.gt.sdifmax) then
            sdifmax = sdif
            nwhenmax = ng
            iii = 1
         endif
         if (istopsearch.eq.0) then
            if (iverbose.gt.0.and.iii.eq.1) then
               call printstatus(1,k,p1,p2,qph,xsec,var,varbefore,sdif,
     .              sdifmax,fmax)
            endif
         endif

!!!!! unweighting for (non radiative events and radiative events) at NLO
         if (k.gt.nsearch.and.(.not.condexp).and.eventlimiter.eq.'unw')
     >    then
            istopsearch = 1
            factor=1.05d0
            if(condnlo) factor=1.2d0
            if (hitpmiss.lt.1.d0) then
               fmax = factor*sdifmax
               print*,'Starting now also unweighted generation!'
            endif
            hitpmiss = hitpmiss + 1.d0
            call getrnd(csi,1)
            if (fmax*csi(1).lt.sdif) then
               hit = hit + 1.d0
               if (store.eq.'yes'.and.eventlimiter.eq.'unw')
     .              call eventstorage_unweightedneg(k,+1.d0,p1cut,
     >              p2cut,qphcut)
            endif

            if (sdif.lt.-1.d-20) then
               call getrnd(csi,1)
               if (fmax*csi(1).lt.-1.d0*sdif) then
                  hitneg = hitneg + 1.d0
                  if (store.eq.'yes'.and.eventlimiter.eq.'unw') then
                     call eventstorage_unweightedneg(k,-1.d0,p1cut,
     >                    p2cut,qphcut)
                  endif
                  nneg = nneg +1.d0
               endif

            endif
            if (sdif.gt.fmax) then
               nover = nover + 1
               sumover  = sumover  + sdif - fmax
               sum2over = sum2over + (sdif - fmax)**2


               if (iverbose.gt.0) then
                  call printstatus(2,k,p1,p2,qph,xsec,var,varbefore,
     .                 sdif,sdifmax,fmax)
               endif
            endif
         endif

!     ! unweightening for unweighted events...
!!!   HERE I INCREASE FMAX TO 1.2D0*FMAX DUE TO DE FACT THAT THE PRELIMINARY RUN IS DONE IN PSRAD
!!!   and deltasv è al max dell'ordine di 0.2

         
*****************I PUT WEIGHTS TO ZERO IN THE PRE RUN
         
         if (k.le.nsearch.and.condexp) then
            sdif=0.d0
         endif
****************
         if (k.gt.nsearch.and.(condexp)) then
            istopsearch = 1
            if (hitpmiss.lt.1.d0) then
               fmax = 1.2d0*sdifmax
               print*,'Starting now also unweighted generation!'
            endif
            hitpmiss = hitpmiss + 1.d0
            call getrnd(csi,1)
!     Pre-filter to decide wether to evaluate exact SV correction
            wfilt = 2.d0
            if (fmax*csi(1).lt.sdif*wfilt) then

               sdifapprox=sdif
               eventsv=eventsv+1
               q=qph(1,:)

               if(ord.eq.'exprad'.or.ord.eq.'exp') then

                  call getrndstatus(isvec) 
                  call svfactor(model,ng,ecms,p1,p2,
     >                 eps,sv,deltasv,q,qph)

                 
!     To avoid instabilities
                 if(abs(deltasv).ge.2.d0) then
                     ltint=1
                    fscoll =1
                    iscoll =1
!     If dispersive pions we want to evaluate the same s',s'', so we reset isvec to the one before svfactor 
                    call rndreset(isvec)
                    call svfactor(model,ng,ecms,p1,p2,eps,
     >                   sv,deltasv,q,qph)
                     ltint=0
                 endif
                 
!!!   If phi decays active (checked that even with this correction abs(deltasv).le.0.5d0 at KLOE)
                 if(indexcorr.ne.0) then
                    deltasv = sdifnophi*deltasv/sdif
                 endif
                 
                  sdif=sdif*(1d0+deltasv)
               endif
               if (sdif.gt.sdifmax) then
                  sdifmax = sdif
                  nwhenmax = ng
                  iii = 1
               endif
               if((1d0+deltasv).ge.deltamax) deltamax=1d0+deltasv
               if((1d0+deltasv).le.deltamin) deltamin=1d0+deltasv
*--------------UNWEIGHTING

               if(fmax*csi(1).lt.abs(sdif)) then                                                               

!     Peso (+1 o -1)
                  weightunw = sign(1.d0, sdif)                                                                 
                                                                                                               
!     Peso corretto legato ai vari bias possibili introdotti
                  thresh_filt = sdifapprox * wfilt
                  thresh_eff  = min(thresh_filt, fmax)

!     Se l'evento sfora la soglia effettiva, applichiamo la ripesatura
                  if (abs(sdif) .ge. thresh_eff) then
                     weightunw = sign(1.d0, sdif) * abs(sdif)
     >                    / thresh_eff
               endif
                  hit = hit + 1.d0*weightunw                                                                             
                 
c$$$                  if (icalc.eq.1.and.ie.eq.0) then
c$$$                     call distributions(weightunw*fmax,k-nsearch,p1cut
c$$$     >                    ,p2cut,qphcut)
c$$$                  endif

                  if (store.eq.'yes'.and.eventlimiter.eq.'unw') then
                     call eventstorage_unweightedneg(k,weightunw,
     >                    p1cut,p2cut,qphcut)
                  endif
                  
!     1. COUNTER NEGATIVE
               if (sdif.lt.-1.d-20) then
                  nneg  = nneg  + 1
               endif

!     2. COUNTER OVER FMAX
                  if (abs(sdif).gt.fmax) then
                  nover = nover + 1
                  endif

!     3. COUNTER PRE-FILTER
                  if (abs(sdif).ge.sdifapprox*wfilt) then
                     noverdelta   = noverdelta + 1
                  endif

               endif
                              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               corrfact=1.d0/(sdifapprox*wfilt/fmax)
               if(corrfact.ge.1.d0) then
                  sdif=sdif*corrfact
               else
                  sdif=sdif
               endif
            else
               sdif=0.d0
            endif
         endif

         if (isnan(sdif)) then
            print*,'NAN found (pos 2) setting sdif to 0!!!!'
            print*,'pid,hostnm,ord,k,ng,emtx,sv,sdif,phsp,w'
            print*,idproc,hstnm,ord,k,ng,emtx,sv,sdif,phsp,w
            sdif = 0.d0
            w = 0.d0
            phsp = 0.d0
            ie = 1
            call printstatus(1,k,p1,p2,qph,xsec,var,varbefore,sdif,
     .           sdifmax,fmax)
            print*,'end ----!!!'
            if (iverbose.gt.0) then
               call printstatus(5,k,p1,p2,qph,xsec,var,varbefore,
     .              sdif,sdifmax,fmax)
            endif
         endif

         kbck=k
         if (k.gt.nsearch.and.condexp) then
            k=k-nsearch
         endif
         if (idbgprnt.eq.1) then
            print*,sdif,emtx,sv,w,wnphot,wcos,phsp
         endif


         sum  = sum  + sdif

         sum2 = sum2 + sdif**2

         if (ng.eq.0) then
            sumborn = sumborn + sdif/sv
            sum2born = sum2born + (sdif/sv)**2
         elseif (ng.eq.1.and.
     >	      (fs.eq.'mr'.or.fs.eq.'er'.or.fs.eq.'pr'))
     &     then
             sumborn = sumborn + sdif/sv
            sum2born = sum2born + (sdif/sv)**2
         endif

         if (ibranch.eq.1) sumb1 = sumb1 + sdif
         if (ibranch.eq.2) sumb2 = sumb2 + sdif

         sump(ng)  = sump(ng)  + sdif
         sum2p(ng) = sum2p(ng) + sdif**2

         varbefore = var
         xsec = sum/k
         var  = sqrt(abs((sum2/k-xsec**2)/k))




c         tollerate = 1.2d0
         tollerate = 1.8d0
         if (var.gt.tollerate*varbefore
     .        .and.varbefore.gt.0.d0.and.ie.eq.0) then
            if (iverbose.gt.0) then
               call printstatus(3,k,p1,p2,qph,xsec,var,varbefore,sdif,
     .              sdifmax,fmax)
            endif
!---------- Nullifying this event ------------!
            ratio = var/varbefore
c            if (k.gt.5000000.and.ratio.gt.3) then
            if (k.gt.5000000.and.ratio.gt.6.and.ord.ne.'born') then
               print*,' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
               print*,' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
               print*,' >>>>>>>> REJECTING THE EVENT <<<<<<<<'
               print*,' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
               print*,' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
               print*,var,varbefore
               if (iverbose.gt.0) then
                  call printstatus(4,k,p1,p2,qph,xsec,var,varbefore,
     .                 sdif,sdifmax,fmax)
               endif
c               print*,'WELL... actually not rejecting it...'
c               goto 12345
               sum = sum - sdif
               sum2 = sum2 - (sdif)**2
               sump(ng) = sump(ng) - sdif
               sum2p(ng) = sum2p(ng) - (sdif**2)
               xsec = sum/k
               var  = sqrt((sum2/k-xsec**2)/k)
               sdif = 0.d0
12345          continue
            endif
         endif

*[[[  START VP PARAMETRIC UNCERTAINTY
         if (iteubn.eq.0) then
            if (invp.gt.0) then
               do ii = 1,invp
                  sumv(ii)  = sumv(ii)  +  sdif*remtx(ii)
                  sum2v(ii) = sum2v(ii) + (sdif*remtx(ii))**2
                  xsecv(ii) = sumv(ii)/k
                  varv(ii)  = sqrt(abs((sum2v(ii)/k-xsecv(ii)**2)/k))
               enddo
            endif
         else
            if (invpt.gt.0) then
               do ii = 1,nvpst
                  sumvt(ii)  = sumvt(ii)  +  sdif*remtxt(ii)
                  sum2vt(ii) = sum2vt(ii) + (sdif*remtxt(ii))**2
                  xsecvt(ii) = sumvt(ii)/k
                  varvt(ii)  = sqrt(abs((sum2vt(ii)/k-xsecvt(ii)**2)/k))
               enddo
            endif
         endif
*]]]  END VP PARAMETRIC UNCERTAINTY

         if (icalc.eq.1.and.ie.eq.0) then
c            call distributions(sdif,k,p1,p2,qph)
            call distributions(sdif,k,p1cut,p2cut,qphcut)
            if (store.eq.'yes'.and.eventlimiter.eq.'w')
     .           call eventstorage_weighted(k,sdif,p1cut,p2cut,qphcut)
         endif
*



**[[WRITING OUTPUT]]
         writecond = kbck.eq.npoints.or.hit.eq.npoints
         if (iwriteout.gt.0) then
            writecond = writecond.or.mod(k,iwriteout).eq.0
         else
            if (mod(k,3000000).eq.0) then
               call cpu_time(tmptime)
               if ((tmptime-starttime).ge.-iwriteout) then
                  writecond = .true.
               endif
            endif
         endif
         if (writecond) then
            kwritenow = k
            call cpu_time(tmptime)
            dtempo = tmptime-starttime
            starttime = tmptime

            call writedistributions
!! hit or miss cross section...
            hmxsect = 0.d0
            hmerr   = 0.d0
            if (hitpmiss.gt.0.d0)  then
               hmeff   = hit/hitpmiss
               hmxsect = fmax*hmeff
               hmerr   = fmax * sqrt(hmeff*(1-hmeff)/hitpmiss)
               if(.not.condexp) then
                 hmnegeff=hitneg/hitpmiss
                 sezneg    = fmax*hmnegeff
                 errsezneg = fmax * sqrt(hmnegeff*(1-hmnegeff)/hitpmiss)
               endif
            endif
!!
            xsec = sum/k
            xsecb=sumb/k
      	    xsecf=sumf/k

            var  = sqrt((abs(sum2/k-xsec**2))/k)

            xsecborn = sumborn/k
            varborn  = sqrt((abs(sum2born/k-xsecborn**2))/k)

            open(10,file=outfile,status='unknown')
            isecondi = int(tmptime-truestarttime)
            call getdayshuursminutesseconds(isecondi,ida,iho,imi,ise)
            write(10,*)' Process PID:',idproc,' running on ',hstnm
            write(10,*)' started ',ida,' days, ',iho,' hours, ',
     >           imi,' minutes and ',ise,' seconds ago'
            write(10,*)' '
            write(10,*)' Output last time written ~',
     >           dtempo,' s ago'

            idk = kwritenow-kwritelast
            write(10,*)' in which ',idk,' events '//
     >           ' have been generated (',(1.d9*dtempo)/(1.d0*idk),
     >           ' ns/event)'
            kwritelast = kwritenow

            write(10,*)' All events generated at the rate:',
     >           1.d9*tmptime/(1.d0*k),' ns/event'

            write(10,*)' '
            write(10,'(1x,A,A,A)')
     >           'final state = ',fs,' '
            write(10,'(1x,A,f12.4,A)')
     >           'nominal ecms   =',ecmsnom,' GeV'
            write(10,'(1x,A,f12.4,A)')
     >           'thmin for charged FS=',thmin*180.d0/pi,' deg'
         if(fs.eq.'mr'.or.fs.eq.'er'.or.fs.eq.'pr') then
            write(10,'(1x,A,f12.4,A)')
     >           'thmax for charged FS=',thmax*180.d0/pi,' deg'
            write(10,'(1x,A,f12.4,A)')
     >           'thmin for signal photon=',thgmin*180.d0/pi,' deg'
            write(10,'(1x,A,f12.4,A)')
     >           'thmax for signal photon=',thgmax*180.d0/pi,' deg'

            write(10,'(1x,A,f12.4,A)')
     >           'emin for signal photon  = ',egmin,' GeV'
     	 endif
            write(10,'(1x,A,f12.4,A)')
     >           'min invariant XX mass =',amassainvmin,' GeV'
c            write(10,'(1x,A,f12.4,A)')
c     >           'max invariant XX mass =',amassainvmax,' GeV'
            write(10,'(1x,A,f12.4,A)')
     >           'acoll. =',zmax*180.d0/pi,' deg'
            write(10,'(1x,A,f12.4,A)')
     >           'emin   = ',emin,' GeV'

            write(10,'(1x,A)')
     >           'ord    = '//ord
            write(10,'(1x,A)')
     >           'model  = '//model
            write(10,'(1x,A,i5)')
     >           'nphot mode =',nphotmode
            write(10,'(1x,A,i9)')
     >           'seed   =',iseed
            write(10,'(1x,A,i5)')
     >           'iarun  =',iarun
            if (iarun.eq.1) then
               if (iteubn.eq.1) then
                  write(10,*)'(Teubner & al. alpha running)'
               elseif (inskvp.eq.1) then
                  write(10,*)'(Novosibirsk alpha running)'
               else
                  write(10,*)'(Jegerlehner alpha running)'
               endif
            endif
            write(10,'(1x,A,f10.9)')
     >           'eps    = ',eps
            write(10,'(1x,A,i5)')
     >           'darkmod =',idarkon
            if(idarkon.eq.1) then
            write(10,'(1x,A,e14.7,A)')
     >           'U mass  = ',amassU, ' GeV'
            write(10,'(1x,A,e14.7,A)')
     >           'U width = ',gammaU, ' GeV'
             write(10,'(1x,A,e14.7,A)')
     >           'k       = ',gvectU/sqrt(4.*pi*alpha),
     >           ' * electric charge'
            write(10,'(1x,A,f7.3,A)')
     >           'egmin   = ', egmin,  ' GeV'
            write(10,'(1x,A,f7.3,A)')
     >           'thgmin  = ', thgmin*180.d0/pi, ' deg'
            write(10,'(1x,A,f7.3,A)')
     >           'thgmax  = ', thgmax*180.d0/pi, ' deg'

         endif
            write(10,'(1x,A,f7.3)')
     >           'deltamax    = ',deltamax
            write(10,'(1x,A,f7.3)')
     >           'deltamin    = ',deltamin
            write(10,*)' '
            if (eventlimiter.eq.'w') then
               write(10,'(A,f12.0,A)')'~ Generating ',wnpoints,
     .              ' weighted events ~'
            else
               write(10,'(A,i12,A)')'~ Generating ',npoints,
     .              ' unweighted events ~'
            endif
            write(10,*)' '
            write(10,*)'::::::>>>>>> weighted events <<<<<<::::::'
            do i = 0,nphmax
               xsecp(i) = sump(i)/k
               varp(i)  = sqrt((abs(sum2p(i)/k-xsecp(i)**2))/k)
               fractions(i) = xsecp(i)/xsec * 100.d0
               iii = i
               if (fs.eq.'gg') iii = i + 2
               write(10,'(i2,A,f16.8,A,f16.8,A,f8.4,A)')
     :              iii,' photons: ',xsecp(i),' +-',varp(i),
     :              '     (',fractions(i),' %)'
            enddo
            write(10,'(1x,A,f16.8,A,f16.8,A)')
     :           'total:   ',xsec,' +-',var,' nb'
            write(10,*)' '
            
            write(10,'(1x,A,i12,A,f20.0)')
     :           'n. ',k,' of ',wnpoints
            eff = (1.d0*naccepted)/k
            write(10,'(1x,A,f6.2,A)')
     :           'cut points ',100.d0 - eff*100,' %'
            write(10,*)'::::::>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<::::::'
            write(10,*)' '
            write(10,*)'::::::>>>>>> unweighted events <<<<<<::::::'
            write(10,'(1x,A,f6.2,A)')
     :           'hit or miss efficiency ',hmeff*100,' %'


      write(10,'(1x,A,f12.0)')
     :           'unweighted events generated ',hit

            biashit        = 0.d0
            biashitpmiss   = 0.d0
            biasneghit     = 0.d0
            biasneghitmiss = 0.d0
            sezover        = 0.d0
            errsezover     = 0.d0
            sezoverpref    = 0.d0
            errsezoverpref = 0.d0
            if (hit.gt.0.d0) then
               biashit      = 1.d0*nover/hit
               biashitpmiss = 1.d0*nover/hitpmiss
               biasneghit   = 1.d0*nneg/hit
               biasneghitmiss = 1.d0*nneg/hitpmiss
               sezover    = sumover/hitpmiss
               errsezover = (sum2over/hitpmiss - sezover**2)/hitpmiss
               errsezover = sqrt(abs(errsezover))
               sezoverpref    = sumoverpref/hitpmiss
               errsezoverpref = (sumoverpref2/hitpmiss - sezoverpref**2)
     >              /hitpmiss
               errsezoverpref = sqrt(abs(errsezoverpref))
               if(condexp) then
                  sezneg    = sumneg/hitpmiss
                  errsezneg = (sum2neg/hitpmiss - sezneg**2)/hitpmiss
                  errsezneg = sqrt(abs(errsezneg))
               endif
            endif
            write(10,*)' '
            write(10,'(1x,A,f16.8,A,f16.8,A)')
     :           'total (nb):   ',hmxsect,' +-',hmerr,' +'
            write(10,'(1x,A,f16.8,A,f16.8,A)')
     :           '       !!!!   ',sezover,' +-',
     :           errsezover,' (bias over fmax) +'
            write(10,'(1x,A,f16.8,A,f16.8,A)')
     :           '       !!!!   ',-sezneg,' +-',
     :           errsezneg,' (bias negative)'
            write(10,'(1x,A,f16.8,A,f16.8)')
     :           'total + biases: ',hmxsect
     :           +sezover-sezneg,' +-',hmerr+errsezover+errsezneg
            write(10,*)' '
             write(10,'(1x,A,i12)')
     :           'N. accepted events: ',naccepted
             write(10,*)' '
             write(10,'(1x,A,i12)')
     :           'N. points with w > fmax (bias): ',nover
             write(10,'(1x,A,i12)')
     :           'N. points with delta > 0.5 : ',noverdelta
             write(10,'(1x,A,i12)')
     :           'N. points with w < 0:',nneg
             write(10,*)' '

            write(10,'(1x,A,f10.7,A,f10.7,A)')
     :           'bias/hit and bias/(hit+missed):',
     :           biashit*100,' % and ',biashitpmiss*100,' %'
            write(10,'(1x,A,f10.5,A,f10.5,A)')
     :           'biases for w < 0:   ',
     :           biasneghit*100,' % and ',biasneghitmiss*100,' %'
            write(10,*)' '
            write(10,'(1x,A,f14.6,A,f14.6)')
     :           'Upper limits fmax and sdifmax ',fmax,' ',sdifmax
            nwhenprint = nwhenmax
            if (fs.eq.'gg') nwhenprint = nwhenmax+2
            write(10,'(1x,A,i3,A)')
     :           'when there were ',nwhenprint,' photons'
            write(10,*)'::::::>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<::::::'

            if (abs(scal_err_dhad).gt.1.d-5) then
               open(15,file=vpuncout,status='unknown')
               write(15,*)' '
               write(15,*)'VP PARAMETRIC UNCERTAINTIES: '
               IF (iteubn.eq.0) THEN
               ii = 0
               do i1 = -1,1
                  do i2 = -1,1
                     do i3 = -1,1
                        if (i1.ne.0.or.i2.ne.0.or.i3.ne.0) then
                           ii = ii + 1
                           diffv(ii) = xsecv(ii)-xsec
                           perc      = (xsecv(ii)-xsec)/xsec*100
                           percerr   = varv(ii)/xsecv(ii)*100
                           write(15,*)ii,'>',i1,i2,i3,xsecv(ii),'+-',
     :                          varv(ii),' (',perc,'%)'
                        endif
                     enddo
                  enddo
               enddo
               write(15,*)' '
               immloc = minloc(xsecv,MASK=xsecv.gt.0.d0)
               write(15,*)'  Minumum x-sect.:',xsecv(immloc(1)),'+-',
     :              varv(immloc(1)),
     :          ' (',(xsecv(immloc(1))-xsec)/xsec*100,'%)'
               write(15,*)'  Central x-sect.:',xsec,'+-',
     :              var
               immloc = maxloc(xsecv,MASK=xsecv.gt.0.d0)
               write(15,*)'  Maximum x-sect.:',xsecv(immloc(1)),'+-',
     :              varv(immloc(1)),
     :          ' (',(xsecv(immloc(1))-xsec)/xsec*100,'%)'
               ELSE ! for Teubner et al. parameterization
               ii = 0
               do ii = 1,nvpst
                  diffvt(ii) = xsecvt(ii)-xsec
                  perc      = (xsecvt(ii)-xsec)/xsec*100
                  percerr   = varvt(ii)/xsecvt(ii)*100
                  write(15,*)ii,'>',ist(ii,:),xsecvt(ii),'+-',varvt(ii),
     :                 ' (',perc,'%)'
               enddo
               write(15,*)' '
               immloc = minloc(xsecvt,MASK=xsecvt.gt.0.d0)
               write(15,*)'  Minumum x-sect.:',xsecvt(immloc(1)),'+-',
     :              varvt(immloc(1)),
     :          ' (',(xsecvt(immloc(1))-xsec)/xsec*100,'%)'
               write(15,*)'  Central x-sect.:',xsec,'+-',
     :              var
               immloc = maxloc(xsecvt,MASK=xsecvt.gt.0.d0)
               write(15,*)'  Maximum x-sect.:',xsecvt(immloc(1)),'+-',
     :              varvt(immloc(1)),
     :          ' (',(xsecvt(immloc(1))-xsec)/xsec*100,'%)'

               ENDIF
               close(15)
            endif
            write(10,*)' '
            write(10,'(1x,A,f16.8,A,f16.8,A)')
     :           'Born x-sect:   ',xsecborn,' +-',varborn,' nb'

            if (fs.eq.'gg') then
               write(10,*)' '
               write(10,*)' B&K total xsect'
               ak = eps/beta
               am = ame/beta/ecms*2.d0
               v  = 0.5d0*log(4.d0/am/am)

               ss = ecms*ecms

         s0bk = 2.d0*pi*alpha**2/ss*(log(s/ame/ame)-1.d0)*convfac

         s3gbk = 2.d0*alpha**3/ss * convfac
         s3gbk = s3gbk*(-(2.d0*v-1.d0)**2*(2.d0*log(ak)+1.d0)+3.d0)

         ssvbk = 2.d0*alpha**3/ss * convfac
         ssvbk = ssvbk * (2.d0*(2.d0*v-1.d0)**2*log(ak)+4.d0/3.d0*v**3
     .   +3.d0*v*v + v*(2.d0/3.d0*pi*pi-6.d0) -1.d0/12.d0*pi*pi)


         write(10,*)'s0    = ',s0bk
         write(10,*)'s_sv  = ',s0bk+ssvbk
         write(10,*)'s_r** = ',s3gbk
         write(10,*)'s_tot = ',s3gbk+s0bk+ssvbk

         cresult = vpolc(ss)

         write(10,'(A,f15.12)')'\Delta\alpha_lep(s) = ',dallep
         write(10,'(A,f15.12)')'\Delta\alpha_top(s) = ',toppart
         write(10,'(A,f15.12,A,f15.12)')'\Delta\alpha_had(s) = ',dalhad,
     .        ' +- ',edalhad

      endif
            write(10,*)' '
            close(10)

         endif
**[[END WRITING OUTPUT]]
         if (eventlimiter.eq.'unw') wnpoints = npoints
         if (eventlimiter.eq.'unw'.and.hit.ge.npoints) continue

         k=kbck
      return
      end
