
c!--------------------------------------------------------------------
c! acpr_v2_0.f [CPC library]
c! 20.12.2011, changes by Stefano Actis:
c! * no test program
c! * output commented
c!--------------------------------------------------------------------

c!ADFKRANLUXI.  ACCELERATION OF RANLUX.  K.G. HAMILTON, F. JAMES.         ADFK
c!      PROGRAM LUXTST                                                        
c!C         Exercise for the RANLUX Pseudorandom number generator.            
c!C                                                                           
c!      DIMENSION RVEC(1000)                                                  
c!      DIMENSION ISDEXT(25)                                                  
c!C                                                                           
c!C         check that we get the right numbers (machine-indep.)              
c!      WRITE (6,'(/A)')  '  CALL RANLUX(RVEC,100)'                           
c!      CALL RANLUX(RVEC,100)                                                 
c!      WRITE (6,'(A/9X,5F12.8)') ' RANLUX default numbers   1-  5:',         
c!     +    (RVEC(L),L=1,5)                                                   
c!      CALL RANLUX(RVEC,100)                                                 
c!      WRITE (6,'(A/9X,5F12.8)') ' RANLUX default numbers 101-105:',         
c!     +    (RVEC(L),L=1,5)                                                   
c!C                                                                           
c!      WRITE (6,'(/A)')  '  CALL RLUXGO(0,0,0,0)'                            
c!      CALL RLUXGO(0,0,0,0)                                                  
c!      CALL RANLUX(RVEC,100)                                                 
c!      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury level 0,   1-  5:',         
c!     +    (RVEC(L),L=1,5)                                                   
c!      CALL RANLUX(RVEC,100)                                                 
c!      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury level 0, 101-105:',         
c!     +    (RVEC(L),L=1,5)                                                   
c!C                                                                           
c!      WRITE (6,'(/A)')  '   CALL RLUXGO(389,1,0,0)'                         
c!      CALL RLUXGO(389,1,0,0)                                                
c!      CALL RANLUX(RVEC,100)                                                 
c!      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p=389,   1-  5:',           
c!     +    (RVEC(L),L=1,5)                                                   
c!      CALL RANLUX(RVEC,100)                                                 
c!      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p=389, 101-105:',           
c!     +    (RVEC(L),L=1,5)                                                   
c!C                                                                           
c!      WRITE (6,'(/A)')  '  CALL RLUXGO(75,0,0,0)'                           
c!      CALL RLUXGO(75,0,0,0)                                                 
c!      CALL RANLUX(RVEC,100)                                                 
c!      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p= 75,   1-  5:',           
c!     +    (RVEC(L),L=1,5)                                                   
c!      CALL RANLUX(RVEC,100)                                                 
c!      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p= 75, 101-105:',           
c!     +    (RVEC(L),L=1,5)                                                   
c!C                                                                           
c!      WRITE (6,'(/A)')  '  test restarting from the full vector'            
c!      CALL RLUXUT(ISDEXT)                                                   
c!      WRITE (6,'(/A/(1X,5I14))') '  current RANLUX status saved:',ISDEXT    
c!      CALL RANLUX(RVEC,100)                                                 
c!      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 1- 5:',                    
c!     +    (RVEC(L),L=1,5)                                                   
c!      CALL RANLUX(RVEC,100)                                                 
c!      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 101-105:',                 
c!     +    (RVEC(L),L=1,5)                                                   
c!C                                                                           
c!      WRITE (6,'(/A)')   '   previous RANLUX status will be restored'       
c!      CALL RLUXIN(ISDEXT)                                                   
c!      CALL RANLUX(RVEC,100)                                                 
c!      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 1- 5:',                    
c!     +    (RVEC(L),L=1,5)                                                   
c!      CALL RANLUX(RVEC,100)                                                 
c!      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 101-105:',                 
c!     +    (RVEC(L),L=1,5)                                                   
c!C                                                                           
c!      WRITE (6,'(/A)')  '     test the restarting by skipping'              
c!      CALL RLUXGO(4,7674985,0,0)                                            
c!      CALL RLUXAT(I1,I2,I3,I4)                                              
c!      WRITE (6,'(A,4I10)')  '  RLUXAT values =',I1,I2,I3,I4                 
c!      DO 150 LI= 1, 10                                                      
c!  150 CALL RANLUX(RVEC,1000)                                                
c!      CALL RLUXAT(I1,I2,I3,I4)                                              
c!      WRITE (6,'(A,4I10)')  '  RLUXAT values =',I1,I2,I3,I4                 
c!      CALL RANLUX(RVEC,200)                                                 
c!      WRITE (6,'(A,2F10.6)')  '  Next and 200th numbers are:',              
c!     +                             RVEC(1), RVEC(200)                       
c!      CALL RLUXGO(I1,I2,I3,I4)                                              
c!      CALL RANLUX(RVEC,200)                                                 
c!      WRITE (6,'(A,2F10.6)')  '  Next and 200th numbers are:',              
c!     +                             RVEC(1), RVEC(200)                       
c!C                                                                           
c!      WRITE (6,'(/A)') ' The following should provoke an error message'     
c!      CALL RLUXGO(4,11111,31,0)                                             
c!      STOP                                                                  
c!C                                                                           
c!      END                                                                   
c!---- The following is RANLUX2.FOR ----                                      

      SUBROUTINE RANLUX(RVEC,LENV)                                          
c                                                                           
C     Subtract-and-borrow random number generator proposed by               
C     Marsaglia and Zaman, implemented by F. James with the name            
C     RCARRY in 1991, and later improved by Martin Luescher                 
C     in 1993 to produce "Luxury Pseudorandom Numbers".                     
c                                                                           
C     Fortran 77 coded by F. James, 1993                                    
C     Converted to perform the recursion with integer arithmetic,           
C        K. G. Hamilton, 1996                                               
c     (Most of these changes are in lower case letters.)                    
C                                                                           
C     References:                                                           
C     M. Luscher, Computer Physics Communications  79 (1994) 100            
C     F. James, Computer Physics Communications 79 (1994) 111               
C     F. James, Computer Physics Communications 97 (1996) 357               
C                                                                           
C   LUXURY LEVELS.                                                          
C   ------ ------      The available luxury levels are:                     
C                                                                           
C  level 0  (p=24): equivalent to the original RCARRY of Marsaglia          
C           and Zaman, very long period, but fails many tests.              
C  level 1  (p=48): considerable improvement in quality over level 0,       
C           now passes the gap test, but still fails spectral test.         
C  level 2  (p=97): passes all known tests, but theoretically still         
C           defective.                                                      
C  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible             
C           correlations have very small chance of being observed.          
C  level 4  (p=389): highest possible luxury, all 24 bits chaotic.          
C                                                                           
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
C!!!  Calling sequences for RANLUX:                                  ++     
C!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++     
C!!!                   32-bit random floating point numbers between  ++     
C!!!                   zero (not included) and one (also not incl.). ++     
C!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++     
C!!!               one 32-bit integer INT and sets Luxury Level LUX  ++     
C!!!               which is integer between zero and MAXLEV, or if   ++     
C!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++     
C!!!               should be set to zero unless restarting at a break++     
C!!!               point given by output of RLUXAT (see RLUXAT).     ++     
C!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++     
C!!!               which can be used to restart the RANLUX generator ++     
C!!!               at the current point by calling RLUXGO.  K1 and K2++     
C!!!               specify how many numbers were generated since the ++     
C!!!               initialization with LUX and INT.  The restarting  ++     
C!!!               skips over  K1+K2*E9   numbers, so it can be long.++     
C!!!   A more efficient but less convenient way of restarting is by: ++     
C!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++     
C!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++     
C!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++     
C!!!                 32-bit integer seeds, to be used for restarting ++     
C!!!      ISVEC must be dimensioned 25 in the calling program        ++     
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      DIMENSION RVEC(LENV)                                                  
      DIMENSION ISEEDS(24), ISDEXT(25)                                      
      PARAMETER (MAXLEV=4, LXDFLT=3)                                        
      DIMENSION NDSKIP(0:MAXLEV)                                            
      DIMENSION NEXT(24)                                                    
      integer   iuni, icarry, maskhi, masklo                                
      PARAMETER (IGIGA=1000000000, JSDFLT=314159265)                        
      PARAMETER (ITWO24=16777216, ICONS=2147483563)                         
      parameter (masklo=itwo24-1)                                           
      SAVE NOTYET, I24, J24, icarry, ISEEDS, TWOM24, TWOM12, LUXLEV         
      SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED                 
      INTEGER LUXLEV                                                        
      LOGICAL NOTYET                                                        
      DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/      
      DATA I24,J24,icarry/24,10,0/                                          
C                               default                                     
C  Luxury Level   0     1     2   *3*    4                                  
      DATA NDSKIP/0,   24,   73,  199,  365 /                               
Corresponds to p=24    48    97   223   389                                 
C     time factor 1     2     3     6    10   on slow workstation           
C                 1    1.5    2     3     5   on fast mainframe             
c                                                                           
c     MASKLO is a word in which the rightmost 24 bits are set and the       
c     higher ones are clear.  MASKHI is MASKLO's complement: one with       
c     24 clear bits on the right, and all the higher ones set.              
c     If using 32-bit words, MASKLO=z'00FFFFFF', and MASKHI=z'FF000000'.    
c     Thus, ANDing with MASKLO gives a 24-bit reduced result,               
c     while ANDing with MASKHI indicates if there was a "carry."            
c                                                                           
      maskhi = not(masklo)                                                  
C                                                                           
C     NOTYET is .TRUE. if no initialization has been performed yet.         
C     Default Initialization by Multiplicative Congruential                 
c                                                                           
      IF (NOTYET) THEN                                                      
         NOTYET = .FALSE.                                                   
         JSEED = JSDFLT                                                     
         INSEED = JSEED                                                     
c         WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED        
         LUXLEV = LXDFLT                                                    
         NSKIP = NDSKIP(LUXLEV)                                             
         LP = NSKIP + 24                                                    
         IN24 = 0                                                           
         KOUNT = 0                                                          
         MKOUNT = 0                                                         
c         WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',        
c     +       LUXLEV,'      p =',LP                                          
c         WRITE(6,*)
         TWOM24 = 1.                                                        
         DO 25 I= 1, 24                                                     
           TWOM24 = TWOM24 * 0.5                                            
           K = JSEED/53668                                                  
           JSEED = 40014*(JSEED-K*53668)-K*12211                            
           IF (JSEED .LT. 0)  JSEED = JSEED+ICONS                           
           ISEEDS(I) = MOD(JSEED,ITWO24)                                    
           NEXT(I) = I-1                                                    
   25    CONTINUE                                                           
         TWOM12 = TWOM24 * 4096.                                            
         NEXT(1) = 24                                                       
         I24 = 24                                                           
         J24 = 10                                                           
         icarry = 0                                                         
         if (iand(iseeds(24),maskhi).ne.0) icarry = 1                       
      ENDIF                                                                 
C                                                                           
C     The Generator proper: "Subtract-with-borrow",                         
C     as proposed by Marsaglia and Zaman,                                   
C     Florida State University, March, 1989                                 
C                                                                           
      DO 100 IVEC= 1, LENV                                                  
      iuni = iseeds(j24) - iseeds(i24) - icarry                             
      if (iand(iuni,maskhi).ne.0) then                                      
         iuni = iand(iuni,masklo)                                           
         icarry = 1                                                         
      else                                                                  
         icarry = 0                                                         
      endif                                                                 
      iseeds(i24) = iuni                                                    
      I24 = NEXT(I24)                                                       
      J24 = NEXT(J24)                                                       
      uni = float(iuni)*twom24                                              
c                                                                           
C     Small numbers (with less than 12 "significant" bits) are "padded,"    
C     and zero is forbidden in case someone takes a logarithm               
c                                                                           
      if (uni .lt. twom12) then                                             
         uni = uni + (float(iseeds(j24))*twom24)*twom24                     
         if (uni .le. 0.0) uni = twom24*twom24                              
      endif                                                                 
      rvec(ivec) = uni                                                      
c                                                                           
C     Skipping to luxury.  As proposed by Martin Luscher.                   
c                                                                           
      IN24 = IN24 + 1                                                       
      IF (IN24 .EQ. 24)  THEN                                               
         IN24 = 0                                                           
         KOUNT = KOUNT + NSKIP                                              
         DO 90 ISK= 1, NSKIP                                                
         iuni = iseeds(j24) - iseeds(i24) - icarry                          
         if (iand(iuni,maskhi).ne.0) then                                   
           iuni = iand(iuni,masklo)                                         
           icarry = 1                                                       
         else                                                               
           icarry = 0                                                       
         endif                                                              
         iseeds(i24) = iuni                                                 
         I24 = NEXT(I24)                                                    
         J24 = NEXT(J24)                                                    
   90    CONTINUE                                                           
      ENDIF                                                                 
  100 CONTINUE                                                              
      KOUNT = KOUNT + LENV                                                  
      IF (KOUNT .GE. IGIGA)  THEN                                           
         MKOUNT = MKOUNT + 1                                                
         KOUNT = KOUNT - IGIGA                                              
      ENDIF                                                                 
      RETURN                                                                
C                                                                           
C--------------------------------------------------------------------       
C     Entry to input integer seeds from previous run                        
C--------------------------------------------------------------------       
C                                                                           
      ENTRY RLUXIN(ISDEXT)                                                  
      NOTYET = .FALSE.                                                      
      TWOM24 = 1.                                                           
      DO 195 I= 1, 24                                                       
      NEXT(I) = I-1                                                         
      TWOM24 = TWOM24 * 0.5
 195  CONTINUE
      NEXT(1) = 24                                                          
      TWOM12 = TWOM24 * 4096.                                               
      WRITE(6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'     
      WRITE(6,'(5X,5I12)') ISDEXT                                           
      DO 200 I= 1, 24                                                       
      iseeds(i) = isdext(i)                                                 
  200 CONTINUE                                                              
      icarry = 0                                                            
      if (isdext(25) .lt. 0) icarry = 1                                     
      ISD = IABS(ISDEXT(25))                                                
      I24 = MOD(ISD,100)                                                    
      ISD = ISD/100                                                         
      J24 = MOD(ISD,100)                                                    
      ISD = ISD/100                                                         
      IN24 = MOD(ISD,100)                                                   
      ISD = ISD/100                                                         
      LUXLEV = ISD                                                          
      IF (LUXLEV .LE. MAXLEV) THEN                                          
         NSKIP = NDSKIP(LUXLEV)                                             
         WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',      
     +                       LUXLEV                                         
      ELSE  IF (LUXLEV .GE. 24) THEN                                        
         NSKIP = LUXLEV - 24                                                
         WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV      
      ELSE                                                                  
         NSKIP = NDSKIP(MAXLEV)                                             
         WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV        
         LUXLEV = MAXLEV                                                    
      ENDIF                                                                 
      INSEED = -1                                                           
      RETURN                                                                
C                                                                           
C--------------------------------------------------------------------       
C     Entry to ouput seeds as integers                                      
C--------------------------------------------------------------------       
C                                                                           
      ENTRY RLUXUT(ISDEXT)                                                  
      DO 300 I= 1, 24                                                       
        isdext(i) = iseeds(i)                                               
  300 CONTINUE                                                              
      ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV              
      if (icarry .ne. 0.)  isdext(25) = -isdext(25)                         
      RETURN                                                                
C                                                                           
C--------------------------------------------------------------------       
C     Entry to output the "convenient" restart point                        
C--------------------------------------------------------------------       
C                                                                           
      ENTRY RLUXAT(LOUT,INOUT,K1,K2)                                        
      LOUT = LUXLEV                                                         
      INOUT = INSEED                                                        
      K1 = KOUNT                                                            
      K2 = MKOUNT                                                           
      RETURN                                                                
C                                                                           
C--------------------------------------------------------------------       
C     Entry to initialize from one or three integers                        
C--------------------------------------------------------------------       
C                                                                           
      ENTRY RLUXGO(LUX,INS,K1,K2)                                           
c                                                                           
c     LUX either selects a luxury level, or specifies a skipping            
c     number (p) directly.                                                  
c                                                                           
      IF (LUX .LT. 0) THEN                                                  
         LUXLEV = LXDFLT                                                    
      ELSE IF (LUX .LE. MAXLEV) THEN                                        
         LUXLEV = LUX                                                       
      ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN                         
         LUXLEV = MAXLEV                                                    
         WRITE (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX           
      ELSE                                                                  
         LUXLEV = LUX                                                       
         DO 310 ILX= 0, MAXLEV                                              
           IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX                       
  310    CONTINUE                                                           
      ENDIF                                                                 
      IF (LUXLEV .LE. MAXLEV)  THEN                                         
         NSKIP = NDSKIP(LUXLEV)                                             
         WRITE(6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',     
     +        LUXLEV,'     P=', NSKIP+24                                    
      ELSE                                                                  
          NSKIP = LUXLEV - 24                                               
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV     
      ENDIF                                                                 
      IN24 = 0                                                              
c                                                                           
c     INS can supply a seed value.  If INS=0, default seed is used.         
c                                                                           
      IF (INS .LT. 0)  WRITE (6,'(A)')                                      
     +   ' Illegal initialization by RLUXGO, negative input seed'           
      IF (INS .GT. 0)  THEN                                                 
        JSEED = INS                                                         
        WRITE(6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',     
     +      JSEED, K1,K2                                                    
      ELSE                                                                  
        JSEED = JSDFLT                                                      
        WRITE(6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'     
      ENDIF                                                                 
      INSEED = JSEED                                                        
      NOTYET = .FALSE.                                                      
      TWOM24 = 1.                                                           
      DO 325 I= 1, 24                                                       
        TWOM24 = TWOM24 * 0.5                                               
        K = JSEED/53668                                                     
        JSEED = 40014*(JSEED-K*53668)-K*12211                               
        IF (JSEED .LT. 0)  JSEED = JSEED+ICONS                              
        ISEEDS(I) = MOD(JSEED,ITWO24)                                       
        NEXT(I) = I-1                                                       
  325 CONTINUE                                                              
      TWOM12 = TWOM24 * 4096.                                               
      NEXT(1) = 24                                                          
      I24 = 24                                                              
      J24 = 10                                                              
      icarry = 0                                                            
      if (iseeds(24) .eq. 0) icarry = 1                                     
c                                                                           
C     If restarting at a break point, skip K1 + IGIGA*K2                    
C     Note that this is the number of numbers delivered to                  
C     the user PLUS the number skipped (if luxury .GT. 0).                  
c                                                                           
      KOUNT = K1                                                            
      MKOUNT = K2                                                           
      IF (K1+K2 .NE. 0)  THEN                                               
        maskhi = not(masklo)                                                
        DO 500 IOUTER=1,K2+1                                                
          INNER = IGIGA                                                     
          IF (IOUTER .EQ. K2+1)  INNER = K1                                 
          DO 450 ISK= 1, INNER                                              
            iuni = iseeds(j24) - iseeds(i24) - icarry                       
            if (iand(iuni,maskhi).ne.0) then                                
              iuni = iand(iuni,masklo)                                      
              icarry = 1                                                    
            else                                                            
              icarry = 0                                                    
            endif                                                           
            iseeds(i24) = iuni                                              
            I24 = NEXT(I24)                                                 
            J24 = NEXT(J24)                                                 
  450     CONTINUE                                                          
  500   CONTINUE                                                            
C         Get the right value of IN24 by direct calculation                 
        IN24 = MOD(KOUNT, NSKIP+24)                                         
        IF (MKOUNT .GT. 0)  THEN                                            
           IZIP = MOD(IGIGA, NSKIP+24)                                      
           IZIP2 = MKOUNT*IZIP + IN24                                       
           IN24 = MOD(IZIP2, NSKIP+24)                                      
        ENDIF                                                               
C       Now IN24 had better be between zero and 23 inclusive                
        IF (IN24 .GT. 23) THEN                                              
           WRITE (6,'(A/A,3I11,A,I5)')                                      
     +    '  Error in RESTARTING with RLUXGO:','  The values', INS,         
     +     K1, K2, ' cannot occur at luxury level', LUXLEV                  
           IN24 = 0                                                         
        ENDIF                                                               
      ENDIF                                                                 
      RETURN                                                                
      END                                                                   

c!CPC FREE FORMAT DATA 75 CARDS
c!---- The following is LUXTST.OUT ----
c!
c! CALL RANLUX(RVEC,100)
c!RANLUX DEFAULT INITIALIZATION:    314159265
c!RANLUX DEFAULT LUXURY LEVEL =   3      p = 223
c!RANLUX default numbers   1-  5:
c!          0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
c!RANLUX default numbers 101-105:
c!          0.43156743  0.03774416  0.24897110  0.00147784  0.90274453
c!
c! CALL RLUXGO(0,0,0,0)
c!RANLUX LUXURY LEVEL SET BY RLUXGO : 0     P=  24
c!RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED
c!RANLUX luxury level 0,   1-  5:
c!          0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
c!RANLUX luxury level 0, 101-105:
c!          0.41538775  0.05330932  0.58195311  0.91397446  0.67034441
c!
c!  CALL RLUXGO(389,1,0,0)
c!RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
c!RANLUX INITIALIZED BY RLUXGO FROM SEEDS           1           0           0
c!RANLUX luxury p=389,   1-  5:
c!          0.94589490  0.47347850  0.95152789  0.42971975  0.09127384
c!RANLUX luxury p=389, 101-105:
c!          0.02618265  0.03775346  0.97274780  0.13302165  0.43126065
c!
c! CALL RLUXGO(75,0,0,0)
c!RANLUX P-VALUE SET BY RLUXGO TO:   75
c!RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED
c!RANLUX luxury p= 75,   1-  5:
c!          0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
c!RANLUX luxury p= 75, 101-105:
c!          0.25600731  0.23443210  0.59164381  0.59035838  0.07011414
c!
c! test restarting from the full vector
c!
c! current RANLUX status saved:
c!      16156027      16534309      15243811       2751687       6002207
c!       7979506       1301976       4567313       4305996       5872599
c!      12003090       2146823      12606367       4111505       5979640
c!      12739666      10489318      14036909      11729352       8061448
c!       7832659       6069758       3197719       1832730      75080216
c!RANLUX numbers 1- 5:
c!          0.22617835  0..60655993  0.86417443  0.43920082  0.23382509
c!RANLUX numbers 101-105:
c!          0.08107197  0.21466845  0.84856731  0.94078046  0.85626233
c!
c!  previous RANLUX status will be restored
c!FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:
c!        16156027    16534309    15243811     2751687     6002207
c!         7979506     1301976     4567313     4305996     5872599
c!        12003090     2146823    12606367     4111505     5979640
c!        12739666    10489318    14036909    11729352     8061448
c!         7832659     6069758     3197719     1832730    75080216
c!RANLUX P-VALUE SET BY RLUXIN TO:   75
c!RANLUX numbers 1- 5:
c!          0.22617835  0.60655993  0.86417443  0.43920082  0.23382509
c!RANLUX numbers 101-105:
c!          0.08107197  0.21466845  0.84856731  0.94078046  0.85626233
c!
c!    test the restarting by skipping
c!RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
c!RANLUX INITIALIZED BY RLUXGO FROM SEEDS     7674985           0           0
c! RLUXAT values =         4   7674985         0         0
c! RLUXAT values =         4   7674985    161840         0
c! Next and 200th numbers are:  0.019648  0.590586
c!RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
c!RANLUX INITIALIZED BY RLUXGO FROM SEEDS     7674985      161840           0
c! Next and 200th numbers are:  0.019648  0.590586
c!
c!The following should provoke an error message
c!RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
c!RANLUX INITIALIZED BY RLUXGO FROM SEEDS       11111          31           0
c! Error in RESTARTING with RLUXGO:
c! The values      11111         31          0 cannot occur at luxury level    4
c!                                                                            ****
c!                                                                            ****
