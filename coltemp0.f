C******************************************
c this is the old version .................
      SUBROUTINE COLTEMP0(TIME,IPT,NUMNP,AMASS,ADOT,TTTT,
     &                    THICK,TMELT,TBOT,AEFF,DELT,IPLSTRT,IPLT,
     &                    SLOPN,IMELT,BMELT,WTHICK,STHRESHB,SIGMALI)
C******************************************
C THIS SOLVES FOR TIME DEPENDENT VERTICAL TEMPERATURE PROFILE
C WITH SOLVED FOR BASAL TEMPERATURE
C IF BASAL TEMP > PRESSURE MELTING, WILL
C RESOLVE WITH BASAL BC FIXED AT PRESSURE MELTING
C AND CALCULATE BASAL MELTING
C******************************************
      IMPLICIT REAL*8(A-H,O-Z)
#include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      PARAMETER(MMAX=40,BIG=1.E20)
      DIMENSION AMASS(11)
      DIMENSION XXX(MMAX),TTTT(MMAX),COND(MMAX),HEAT(MMAX),RHO(MMAX)
      DIMENSION QQQ(MMAX),SOURCE(MMAX),DELX(MMAX)
      DIMENSION RKX(MMAX),FFFF(MMAX),AAAD(MMAX),AAAU(MMAX),AAAL(MMAX)
      DIMENSION CCCD(MMAX),CCCU(MMAX),CCCL(MMAX)
      DIMENSION WWWW(MMAX),TOLD(MMAX),AT(MMAX),UC(MMAX),DDATA(9)
      DIMENSION U(MMAX),TTTT0(MMAX)
      DIMENSION LMAP(MMAX)
      DIMENSION RJUNK(MMAX),RHOC(MMAX)
      common /heat/ hs
      SIGMAL=-SIGMALI
      SIGADJ=-SIGMALI
C      SIGMAL=-AMASS(11)
      PG = 0.089866*0.3816
      THIRD=1.D0/3.D0
C        PRINT *,'IN COLTEMP0'
C        PRINT *,' IPT = ',IPT
C        PRINT *,' NUMNP = ',NUMNP
C        PRINT *,' AMASS = ',AMASS(11)
C        PRINT *,' TTTT = ',TTTT(1)
C        PRINT *,' THICK = ',THICK
C        PRINT *,' TBOT = ',TBOT
C        PRINT *,' AEFF = ',AEFF
C        PRINT *,' DELT = ',DELT
C        PRINT *,' IPLSTRT = ',IPLSTRT
C        PRINT *,' IPLT = ',IPLT
C        PRINT *,' WTHICK = ',WTHICK
      NPT=MMAX
      NDT=10
      WWWW(1)=0.
C ... QQQ IS FOR HORIZONTAL ADVECTION (MULT BY U (HORIZONTAL VELOCITY AT
C ... EACH DEPTH)
C ... SOURCE IS FOR INTERNAL (SHEAR) HEAT GENERATION
      QQQS=0.0
      WWWW(1)=ADOT
C
      LTOT=0
      DO I=1,NPT/2
        LMAP(I)=I
c        LMAP(I)=1
        LTOT=LTOT+LMAP(I)
      ENDDO
      DO I=NPT/2+1,NPT
        LMAP(I)=LMAP(I-1)-1
c        LMAP(I)=1
        LTOT=LTOT+LMAP(I)
      ENDDO
      XXX(1)=THICK
      XXX(NPT)=0.
      DX=(XXX(NPT)-XXX(1))/(LTOT)
      QQQ(1)=0.0
      TOLD(1)=TTTT(1)
      DO I=2,NPT
        XXX(I)=XXX(I-1)+LMAP(I-1)*DX
        QQQ(I)=0.
        TOLD(I)=TTTT(I)
      ENDDO
C
      TOLER=.001
      FACTOR=1.
      XMAX=-1E30
      XMIN=1E30
9     CONTINUE
      DO I=1,NPT
        XMAX=MAX(XMAX,XXX(I))
        XMIN=MIN(XMIN,XXX(I))
        DEPTH=XXX(1)-XXX(I)
        RHO(I)=DENSITY(DEPTH)
        COND(I)=CONDUCT(TTTT(I),RHO(I))
        HEAT(I)=SPHEAT(TTTT(I))
C*************USING FIXED VALUES OF RHO,COND,HEAT*****************
c         RHO(I)=DENSITY(1000.D0)
c         COND(I)=CONDUCT(-30.D0,RHO(I))
c         HEAT(I)=SPHEAT(-30.D0)
C*************USING FIXED VALUES OF RHO,COND,HEAT*****************
C ..... everything in one term...
c        RKX(I)=COND(I)/RHO(I)/HEAT(I)
c        RHOC(I)=1D0
c ..... ala paterson pg 224
        RKX(I)=COND(I)
        RHOC(I)=RHO(I)*HEAT(I)
        IF(I.GT.1) DELX(I-1)=-(XXX(I)-XXX(I-1))
        AAAD(I)=0.
        AAAU(I)=0.
        AAAL(I)=0.
        CCCD(I)=0.
        CCCU(I)=0.
        CCCL(I)=0.
        FFFF(I)=0.
        U(I)=0.0
        SOURCE(I)=0.0
      ENDDO
C           
C          
      IF(IPLSTRT.EQ.0 .AND. IPLT.EQ.1) THEN
        CALL GRSTRT(600,600)
        IPLSTRT=1
      ENDIF
      TMIN=-60.
      TMAX=60.
      IF(IPLT.EQ.1) THEN
        CALL WINDOW(REAL(TMIN),REAL(TMAX),REAL(XMIN),REAL(XMAX))
        CALL NEWPAG
        CALL LINCLR(1)
        CALL MOVE(0.,REAL(XMIN))
        CALL DRAW(0.,REAL(XMAX))
        CALL MRKCLR(2)
      ENDIF
      WWWW(NPT)=0. 
C --- VERTICAL VELOCITY DISTRIBUTION: LINEAR
      SLOPE=(WWWW(NPT)-WWWW(1))/(XXX(NPT)-XXX(1))
      DO I=2,NPT-1
        WWWW(I)=WWWW(1)+SLOPE*(XXX(I)-XXX(1))
      ENDDO
c      DO I=1,NPT
c        QQQ(I)=QQQS
c      ENDDO
      DO I=1,NPT                                                    
        TTTT0(I)=TTTT(I)                                                
      ENDDO                                                          
C.....BUILD CAPACITANCE MATRIX                                              
      CCCD(1)=THIRD*DELX(1)*RHOC(1)                                            
      CCCD(NPT)=THIRD*DELX(NPT-1)*RHOC(NPT-1)                                       
      DO I=2,NPT-1                                                  
        CCCD(I)=THIRD*(DELX(I)*RHOC(I)+DELX(I-1)*RHOC(I-1))                               
      ENDDO                                                          
      DO I=1,NPT-1                                                  
        CCCU(I)=.5*THIRD*DELX(I)*RHOC(I)
        CCCL(I)=CCCU(I)                                                 
      ENDDO                                                          
C.....END CAPACITANCE MATRIX                                                
C --- LOOP ON ITERATION STEP    
      IF(DELT.GT.0.) THEN
        NDT=1
      ELSE
        NDT=5
      ENDIF
      DO 900 IDT=1,NDT 
C        PRINT *,'ITERATION=',IDT
        IF(IPLT.EQ.1) THEN
          CALL WINDOW(REAL(TMIN),REAL(TMAX),REAL(XMIN),REAL(XMAX))
        ENDIF
C.....COME IN HERE FOR ITERATION ON VERTICAL VELOCITY
700   CONTINUE   
C ... CALCULATE THE INTERNAL HEAT GENERATION DUE TO SHEAR
        DO I=1,NPT
          SOURCE(I)=0.D0
        ENDDO                                             
        DO I=2,NPT
          TMID=0.5*(TTTT(I-1)+TTTT(I))
          ATI=HARDNESS(TMID)
          DEPTH=XXX(1)-0.5*(XXX(I-1)+XXX(I))
          DLX=(XXX(I-1)-XXX(I))
          RHOM=DENSITY(DEPTH)
          HEATM=SPHEAT(TMID)
C          SMID=0.5*DLX*((PG*DEPTH*SLOPN)**4/ATI**3)/RHOM/HEATM
c ....... paterson is 2*epsilon-dot*tau rate of heat/unit volume (times element
c         thickness to get total heat deposited in element)
c         following is 1/2 that, to be deposited at top and bottom of element
c ....... following is in bar/yr
          SMID=DLX*((PG*DEPTH*SLOPN)**4/ATI**3)
c ....... following converstd to cal/yr
          SMID=SMID*2.392D4
C          print *,real(depth),real(slopn),real(ati),
C     &            real(tmid),real(smid)
          SOURCE(I-1)=SOURCE(I-1)+SMID
          SOURCE(I)=SOURCE(I)+SMID
        ENDDO
C ..... THIS DEPOSITS ALL THE HEAT AT THE BED RATHER THAN WITHIN THE COLUMN
        contrib=0.D0
        DO I=1,NPT
          contrib=contrib+SOURCE(I)
          SOURCE(I)=0.0D0
        ENDDO 
        SIGADJ=SIGMAL-CONTRIB
C ........................................................................
c          print *,real(contrib),real(slopn),
c     &            real(xxx(1)-xxx(npt)),
c     &            real(pg*slopn*(xxx(1)-xxx(npt)))
C          print *,'contrib',contrib,sigmal
c        DO I=1,NPT
c          PRINT *,I,SOURCE(I)
c        ENDDO
C ... END INTERNAL HEAT GENERATION SECTION
C ----- FORM LOAD VECTOR 
        FFFF(1)=.5*U(1)*QQQ(1)*DELX(1)+SOURCE(1)                             
        FFFF(NPT)=.5*U(NPT-1)*QQQ(NPT-1)*DELX(NPT-1)+
     &             SOURCE(NPT)-SIGADJ          
        P1=FFFF(1)                                                      
        PN=FFFF(NPT)
        DO I=2,NPT-1    
          FFFF(I)=.5*(U(I)*QQQ(I)*DELX(I)+
     &                U(I-1)*QQQ(I-1)*DELX(I-1))+SOURCE(I)      
        ENDDO  
C         DO I=1,NPT
C           RJUNK(I)=FFFF(I)
C           PRINT *,'RJUNK',I,RJUNK(I)
C         ENDDO
C
C ----- END FORM LOAD VECTOR  
C                                                
C ----- FORM STIFFNESS MATRIX                                                 
C ... OLD FORM
C        AAAD(1)=RKX(1)/DELX(1)+WWWW(1)*.5                               
C        AAAD(NPT)=RKX(NPT-1)/DELX(NPT-1)-WWWW(NPT)*.5                   
C        DO I=2,NPT-1                                                
C          AAAD(I)=RKX(I-1)/DELX(I-1)+RKX(I)/DELX(I)                     
C          AAAD(I)=AAAD(I)+.5*(WWWW(I-1)-WWWW(I))                        
C        ENDDO                                                        
C        DO I=1,NPT-1                                                
C          AAAU(I)=-RKX(I)/DELX(I)                                       
C          AAAL(I)=AAAU(I)-WWWW(I)*.5                                    
C          AAAU(I)=AAAU(I)+WWWW(I)*.5                                    
C        ENDDO                                                        
C ...   NEW FORM
        AAAD(1)=RKX(1)/DELX(1)-RHOC(1)*WWWW(1)*.5                               
        AAAD(NPT)=RKX(NPT-1)/DELX(NPT-1)+RHOC(NPT-1)*WWWW(NPT)*.5                   
        DO I=2,NPT-1   
          AAAD(I)=RKX(I-1)/DELX(I-1)+RKX(I)/DELX(I)                     
          AAAD(I)=AAAD(I)+.5*(RHOC(I-1)*WWWW(I-1)-RHOC(I)*WWWW(I))                        
        ENDDO   
        DO I=1,NPT-1                                                
          AAAU(I)=-RKX(I)/DELX(I)                                       
          AAAL(I)=AAAU(I)-RHOC(I)*WWWW(I)*.5                                    
          AAAU(I)=AAAU(I)+RHOC(I)*WWWW(I)*.5 
        ENDDO 
c save values BEFORE effective matrix is formed
        FFFFNPT=FFFF(NPT)+SIGMAL                                                    
        AAADNPT=AAAD(NPT)
        AAALNPT=AAAL(NPT-1)
c        print *,0d0,cccd(1),cccu(1)
c        do i=2,npt-1
c          print *,cccl(i-1),cccd(i),cccu(i)
c        enddo
c        print *,cccl(npt-1),cccd(npt)
c        pause
c        print *,0d0,aaad(1),aaau(1)
c        do i=2,npt-1
c          print *,aaal(i-1),aaad(i),aaau(i)
c        enddo
c        print *,aaal(npt-1),aaad(npt)
c        pause
C.......FORM EFFECTIVE LOAD VECTOR                                            
        IF(DELT.GT.0.) THEN
C      PRINT *,'FORM EFFECTIVE LOAD VECTOR',DELT
C          FFFF(1)=(FFFF(1)*DELT+                                          
C     &           (CCCD(1))*TTTT(1)+                                       
C     &           (CCCU(1))*TTTT(2)) 
C          FFFF(NPT)=(FFFF(NPT)*DELT+                                      
C     &           (CCCD(NPT))*TTTT(NPT)+                                   
C     &           (CCCL(NPT-1))*TTTT(NPT-1))                               
C          DO I=2,NPT-1                                                
C            FFFF(I)=(FFFF(I)*DELT+                                        
C     &           (CCCD(I))*TTTT(I)   +                                    
C     &           (CCCL(I-1))*TTTT(I-1) +                                  
C     &           (CCCU(I))*TTTT(I+1) )  
          FFFF(1)=FFFF(1)+
     &                   (CCCD(1)*TTTT(1)+                                       
     &                    CCCU(1)*TTTT(2))/DELT
          FFFF(NPT)=FFFF(NPT)+
     &                    (CCCD(NPT)*TTTT(NPT)+                                   
     &                     CCCL(NPT-1)*TTTT(NPT-1))/DELT
          DO I=2,NPT-1                                                
            FFFF(I)=FFFF(I)+
     &                     (CCCD(I)*TTTT(I)     +                                    
     &                      CCCL(I-1)*TTTT(I-1) +                                  
     &                      CCCU(I)*TTTT(I+1))/DELT  
C      PRINT *,'T',I,REAL(TTTT(I)),REAL(TTTT(I-1)),REAL(TTTT(I+1))
C      PRINT *,'C',I,REAL(CCCD(I)),REAL(CCCL(I-1)),REAL(CCCU(I))                                 
          ENDDO    
C         FORM EFFECTIVE STIFFNESS MATRIX                                       
          DO I=1,NPT                                                  
C            AAAD(I)=DELT*AAAD(I)+CCCD(I)                                  
C            AAAU(I)=DELT*AAAU(I)+CCCU(I)                                  
C            AAAL(I)=DELT*AAAL(I)+CCCL(I)     
            AAAD(I)=AAAD(I)+CCCD(I)/DELT                                  
            AAAU(I)=AAAU(I)+CCCU(I)/DELT                                  
            AAAL(I)=AAAL(I)+CCCL(I)/DELT                                  
          ENDDO            
        ENDIF   
c save values AFTER effective matrix is formed                                         
        FFFFNPT=FFFF(NPT)+SIGMAL                                                    
        AAADNPT=AAAD(NPT)
        AAALNPT=AAAL(NPT-1)
C
C       DO I=1,NPT                                                  
C         WRITE(7,*) 'EFF F',I,REAL(FFFF(I)),REAL(RJUNK(I))                                     
C       ENDDO  
C ----- FIX BOUNDARY CONDITION FOR NODE 1                                     
        RK11=AAAD(1)                                                    
        RK12=AAAU(1)                                                    
        RK21=AAAL(1)                                                    
        AAAD(1)=1.                                                      
        AAAU(1)=0.                                                      
        AAAL(1)=0.                                                      
        FFFF(1)=TTTT(1)                                                 
        FFFF(2)=FFFF(2)-RK21*FFFF(1)                                    
C ----- FIX BOUNDARY CONDITION FOR NODE NPT
        RKN1=AAAD(NPT)                                                  
        RKN2=AAAU(NPT-1)                                                
        RK2N=AAAL(NPT-1)    
        IPASS=1     
        ISET=0                                       
C ..... IF THERE IS WATER PRESENT (at least STHRESHB), 
C ..... SET TO FIXED MELTING PT VALUE
c      print *,'wthick',wthick,sthresh
        IF(WTHICK.GT.STHRESHB) THEN
c          PRINT *,IPT,' SET TO FIXED MELTING',
c     &              REAL(WTHICK),REAL(STHRESHB)
          TTTT(NPT)=TMELT
          IPASS=2
        ENDIF
C ----- THE MATRIX SOLUTION
C*****JUMP BACK TO HERE FOR 2ND PASS, FIXED BASAL BC
800   CONTINUE
        IF(IPASS.EQ.2) THEN
          ISET=1
c          IPASS=1
          AAAD(NPT)=1.                                                    
          AAAU(NPT-1)=0.                                                  
          AAAL(NPT-1)=0.                                                  
          FFFF(NPT)=TTTT(NPT)                                             
          FFFF(NPT-1)=FFFF(NPT-1)-RKN2*FFFF(NPT)                          
        ENDIF
C        PRINT *,' BEFORE ',(REAL(FFFF(KKK)),KKK=3,7)
        CALL TRI(NPT,AAAL,AAAD,AAAU,FFFF,TTTT) 
C        PRINT *,' AFTER ',(REAL(FFFF(KKK)),KKK=3,7)
        IF(TTTT(NPT).GT.TMELT .AND. IPASS.EQ.1) THEN
C ..... CHECK FOR BED ABOVE MELTING POINT ...
c          PRINT *,REAL(TTTT(NPT)),REAL(TMELT),' ABOVE MELTING POINT ',
c     &           IPT
          TTTT(NPT)=TMELT
          IPASS=2
          GOTO 800
C        ELSEIF(TTTT(NPT).LT.TMELT .AND. WTHICK.GT.0.0) THEN
C ..... CHECK FOR BED BELOW MELTING POINT WITH WATER AT BED ...
C          PRINT *,REAL(TTTT(NPT)),REAL(TMELT),' BELOW MELTING POINT',
C     &           IPT
C          TTTT(NPT)=TMELT
C          IPASS=2
C          GOTO 800
        ENDIF
C
c ..... experiment to "SMOOTH" temps ...
c        CALL SMOOTH(TTTT,NPT)
C .....
C
C --------------------------
C ..... SCAN THE COLUMN FOR TEMPERATURES ABOVE THE PRESSURE MELTING
C ..... POINT
        ISCAN=0
        DO I=1,NPT-1
          DEPTH=XXX(1)-XXX(I)
          PMPLOC=PMP(DEPTH)
          IF(TTTT(I).GT.PMPLOC) THEN
            TTTT(I)=PMPLOC
            ISCAN=ISCAN+1
          ENDIF
C         TTTT(I)=MIN(TTTT(I),PMPLOC)
        ENDDO
        IF(ISCAN.GT.0) PRINT *,ISCAN,' POLYTHERMAL AT ',IPT
C ..................................................................
        AAAD(1)=RK11                                                    
        AAAU(1)=RK12                                                    
        AAAL(1)=RK21                                                    
        AAAD(NPT)=RKN1                                                  
        AAAU(NPT-1)=RKN2                                                
        AAAL(NPT-1)=RK2N                                                
C ------FOLLOWING FOR RKX FUNCTION OF TEMPERATURE                             
        DO I=1,NPT-1                                                
C --------RKXI= FUNCTION OF TEMPERATURE                                
          DEPTH=XXX(1)-XXX(I)  
          RHO(I)=DENSITY(DEPTH)                                             
          COND(I)=CONDUCT(TTTT(I),RHO(I))
          HEAT(I)=SPHEAT(TTTT(I))
C********************************************
c         RHO(I)=DENSITY(1000.D0)
c         COND(I)=CONDUCT(-30.D0,RHO(I))
c         HEAT(I)=SPHEAT(-30.D0)
C********************************************
C ..... everything in one term...
          RKXI=COND(I)/RHO(I)/HEAT(I)                                     
c ..... ala paterson pg 224
          RKXI=COND(I)
          RKX(I)=((FACTOR-1.)*RKX(I)+RKXI)/FACTOR                      
        ENDDO                                                        
c        SIGMA0=P1-((RK11-CCCD(1))*TTTT(1)+
c     &         RK12*TTTT(2)-CCCU(1)*TTTT0(2))
C        SIGMA0=-SIGMA0/DELT                                             
c        SIGMALC=PN-((RKN1-CCCD(NPT))*TTTT(NPT)+
c     &         RKN2*TTTT(NPT-1)-CCCU(NPT-1)*TTTT0(NPT-1))   
c        DTEMP=(TTTT(NPT-1)-TTTT(NPT))/(XXX(NPT-1)-XXX(NPT))
c        SIGMABC=-RKX(NPT-1)*DTEMP
        SIGMABC=FFFFNPT-AAALNPT*TTTT(NPT-1)-AAADNPT*TTTT(NPT)
C---------------------------------------
C WHAT IS THIS ????
C
c        SIGMALC=-SIGMALC/DELT                                             
C--------------------------------------------
C ..... ONLY CALC BMELT IF TOGGLE IS ON...
        IF(IMELT.EQ.1) THEN
          IF(ISET.EQ.0) THEN
            BMELT=0.0
c            BMELT=-1D-5
          ELSE
            BMELT=-(SIGMAL-SIGMABC)/RKX(NPT-1)
c ... limit freezing to 10 mm/yr ...
            IF(BMELT.LT.0) BMELT=max(-0.01D0,BMELT)
          ENDIF
c ... this is basic offset .02 mm 
          BMELT=BMELT-2D-5
c ... let it play, whatever it makes...
          BMELT=-(SIGMAL-SIGMABC)/RKX(NPT-1)
          PRINT *,SIGMAL,BMELT,SIGMABC
C ...     DISABLED ITERATION ON VERTICAL VELOCITY BASED ON BASAL MELT RATE
          IF(.false.) THEN
            IF(ABS(WWWW(NPT)-BMELT).GT.TOLER) THEN
C             PRINT *,ipt,real(WWWW(NPT)),real(BMELT),
C     &              ' GOING AROUND AGAIN'
              WWWW(NPT)=BMELT
C...........  LOOP BACK TO BEGINNING AFTER RECALCULATING VERTICAL VELOCITY
              SLOPE=(WWWW(NPT)-WWWW(1))/(XXX(NPT)-XXX(1))
              DO I=2,NPT-1
                WWWW(I)=WWWW(1)+SLOPE*(XXX(I)-XXX(1))
              ENDDO
              GOTO 700
            ENDIF
          ENDIF
        ENDIF
C ...   TO HERE
        DO I=1,NPT-1                                                
          DTDX=(TTTT(I+1)-TTTT(I))/DELX(I)                              
          FLUX=-RKX(I)*DTDX 
        ENDDO                        
        IF(IPLT.EQ.1) THEN
          CALL LINCLR(1)                                             
          CALL MOVE(REAL(TOLD(1)),REAL(XXX(1)))                           
          DO I=2,NPT                                                  
            CALL DRAW(REAL(TOLD(I)),REAL(XXX(I)))                         
          ENDDO 
          CALL LINCLR(2)                                             
          CALL MOVE(REAL(TTTT(1)),REAL(XXX(1)))                           
          DO I=2,NPT                                                  
            CALL DRAW(REAL(TTTT(I)),REAL(XXX(I)))                         
          ENDDO                                                        
        ENDIF
        DIFF=0D0
        DO I=1,NPT
          DIFF=DIFF+(TOLD(I)-TTTT(I))**2
          TOLD(I)=TTTT(I)
          TTTT0(I)=TTTT(I)
        ENDDO
C        PRINT *,'DIFF=',DIFF
        CALL VELO(TIME,NPT,5D-3,XXX,TTTT,RHO,AT,U,UC,
     &            BMELT,DDATA,AEFF)                                                         
900   CONTINUE            
C ... TBOT IS DIFFERENCE BETWEEN BASAL TEMP AND PRESSURE MELTING POINT
      TBOT=TTTT(NPT)-TMELT
C....................................................................
C      PRINT *,' TBOT = ',TBOT
C        CALL WAIT(10000)
C      if(adot.lt.0) CALL WAIT(100000)

c ++++++++ EISMINT STUFF FOR SLIDING EXPERIMENT +++++++++++++++++++++++++++++
c      if (IPT.eq.1861) call payneout2(time,IPT,mmax,tttt,xxx,at,wwww)
c      if (IPT.eq.2471) call payneout2(time,IPT,mmax,tttt,xxx,at,wwww)
c      if (IPT.eq.2776) call payneout2(time,IPT,mmax,tttt,xxx,at,wwww)
c      if (IPT.eq.2109) call payneout2(time,IPT,mmax,tttt,xxx,at,wwww)
c      if (IPT.eq.2719) call payneout2(time,IPT,mmax,tttt,xxx,at,wwww)
c ========================================================================
      RETURN
C      IF(IPLT.EQ.1) THEN
C        CALL WAIT(10000)
C        CALL GRSTOP1
C      ENDIF
      END         
