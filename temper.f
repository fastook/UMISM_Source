C----------------------------------------------
      SUBROUTINE TEMPER(MXX, NUMNP, NUMEL, NTYPE, KX, KODE, ADOT, FRACT,
     &            DEPB, HTICE, PSURF, FLOWA, SLDGB, TEMP, ITYPE, TBED,
     &            X, Y, AMASS, SLOPN, TIME, DT, IPLOT, TBAVG,
     &            AFUDGE,GRIDAREA,IMELT,BMELT,WTHICK,VOL,AREA,
     &            GEOFLUX,VEL,IDT,CONTRIB)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM,MMAX=40)
      PARAMETER(NSMAX=MAXTIME)
C USE SURFACE TEMPERATURE TO CALCULATE BED TEMPERATURE, AND HENCE MODIFY
C SLIDING DISTRIBUTION
C
      EXTERNAL ACCUM,PMP
      DIMENSION NTYPE(NMAX),ADOT(NMAX),FRACT(NMAX),X(NMAX),Y(NMAX),
     &          DEPB(NMAX),HTICE(NMAX),FLOWA(NMAX),TBED(NMAX),
     &          SLDGB(NMAX),TEMP(NMAX),LM(4),KX(NMAX,4),AMASS(11),
     &          SLOPN(4,NMAX),KODE(NMAX),ITYPE(NMAX),AFUDGE(NMAX),
     &          ZZZ(MMAX),BMELT(NMAX),WTHICK(NMAX),VEL(NMAX,3),
     &          GEOFLUX(NMAX),BMELTOLD(NMAX),IDT(NMAX),PSURF(NMAX),
     &          CONTRIB(NMAX)
      COMMON /LINE/ NP,NLINE(1000)
      COMMON /TEMPERATURE/ TEMPA(MMAX,NMAX)
      COMMON /TCONST/ TSORIG(NMAX),TFLAG
      LOGICAL TFLAG
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      LOGICAL ISICE,INNLINE
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      COMMON /HEAT/ HTOGG
      LOGICAL HTOGG,HSAVE
      logical file41,file43
      SAVE ISTART,JSTART,BMELTOLD,KSTART
      DATA ISTART /0/, JSTART /0/, KSTART /0/
      DATA CONV /3.D0/
C THRESHOLD VALUE OF WATER THICKNESS (10 MM ) TO PRODUCE SLIDING
      STHRESH=10.D0*1D-3
C THRESHOLD VALUE OF WATER THICKNESS (1 MM ) TO PRODUCE SLIDING
      STHRESH=1.D0*1D-3
C THRESHOLD VALUE OF WATER THICKNESS (0.1 MM ) TO PRODUCE FREEZING
      STHRESHB=0.1D0*1D-3
      IF(MXX.NE.NMAX) THEN
        PRINT *,'CHANGE NMAX IN TEMPER',NMAX,MXX
        STOP
      ENDIF
      !PRINT *,' IN TEMPER ',TIME,DT;pause
C------------------------------------------
      CALL HARD(A0,B0,A1,B1)
      inquire(file='fort.41',exist=file41)
      IF(JSTART.EQ.0. and. file41) THEN
        PRINT *,'READING WATER THICKNESS FILE'
        READ(41,*,END=11) NNN
        IF(NNN.NE.NUMNP) THEN
          PRINT *,'PROBLEMS:'
          PRINT *,'INCOMPATIBLE WITH CURRENT NUMNP=',NUMNP,NNN
          GOTO 11
        ENDIF
        DO I=1,NUMNP
          READ(41,*,END=11) II,WTHICK(I),BMELT(I)
          BMELTOLD(I)=BMELT(I)
          IF(I.NE.II) THEN
            PRINT *,'PROBLEMS:'
            PRINT *,'INCOMPATIBLE WITH CURRENT NMAX=',NMAX
            GOTO 11
          ENDIF
        ENDDO
        JSTART=1
      ENDIF
11    CONTINUE
      IF(JSTART.EQ.0) THEN
        PRINT *,' NONE FOUND, SET TO ZERO ... '
        DO I=1,NUMNP
            WTHICK(I)=0.D0
            BMELT(I)=-1D-5
            BMELTOLD(I)=-1D-5
        ENDDO
        JSTART=1
      ENDIF
      IF(ISTART.EQ.0) THEN
        PRINT *,'READING TEMPERATURE FILE'
        PRINT *,'0:BINARY, 1:ASCII'
        READ(*,*) IREAD
        WRITE(99,*) IREAD
        inquire(file='fort.43',exist=file43)
        if(file43) then
          IF(IREAD.EQ.1) THEN
            DO I=1,NUMNP
              DO J=1,MMAX
                READ(43,*,END=12) II,JJ,TEMPA(J,I)
                IF(I.NE.II .OR. J.NE.JJ) THEN
                  PRINT *,'PROBLEMS:'
                  PRINT *,'INCOMPATIBLE WITH CURRENT MMAX=',MMAX
                  GOTO 12
                ENDIF
              ENDDO
            ENDDO
            ISTART=1
          ELSE
            READ(43) NUMNP1,NPT1
            IF(NUMNP1.NE.NUMNP .OR. NPT1.NE.MMAX) THEN
              PRINT *,'PROBLEMS:'
              PRINT *,'INCOMPATIBLE WITH CURRENT MMAX=',MMAX
              GOTO 12
            ENDIF
            READ(43) ((TEMPA(J,I),J=1,MMAX),I=1,NUMNP)
            ISTART=1
          ENDIF
        endif
      ENDIF
12    CONTINUE
      IF(ISTART.EQ.0) THEN
        PRINT *,'NONE FOUND, STEADY STATE CALCULATION'
        DO I=1,NUMNP
          DO J=1,MMAX
            TEMPA(J,I)=TEMP(I)
          ENDDO
        ENDDO
        DTSAVE=DT
        DT=0.D0
      ENDIF        
C--------------------------------------------          
      IC1=0
      IC2=0
      IFIXED=0
      IMELTING=0
      IMELTED=0
      IFREEZING=0
      IFROZEN=0
      INOICE=0
      IWITHICE=0
      TMIN=1D30
      TMAX=-TMIN
      ISTEP=NUMNP/5
      IF(IPLOT.EQ.7) THEN
        IPLSTRT=0
        IPLT=1
      ELSE
        IPLSTRT=0
        IPLT=0
      ENDIF
      TBAVG=0.D0
      NTAVG=0
      IF(KSTART.EQ.0) THEN
        HSAVE=HTOGG
        HTOGG=.FALSE.
      ENDIF
      DO I=1,NUMNP
C ...... IS I IN NLINE ??
         IF(IPLOT.EQ.7) THEN
           INNLINE=.FALSE.
           DO NN=1,NP
             IF(I.EQ.NLINE(NN)) THEN
               INNLINE=.TRUE.
               GOTO 500
             ENDIF
           ENDDO
500        CONTINUE
           INNLINE=.FALSE.
           IF(INNLINE) THEN
             IPLT=1
           ELSE
             IPLT=0
           ENDIF
         ENDIF
C        IFLAG=0
C        CALL COLUMNZ(MMAX,HTICE(I),DEPB(I),ZZZ)
C        IF(ITYPE(I).NE.0 .AND. KODE(I).EQ.0) THEN
        IF(ITYPE(I).NE.0) THEN
C          IF(MOD(I,ISTEP).EQ.0) PRINT *,I
C ******* REMOVE ACCUM, REPLACE WITH AFUNCT EXPERIMENTAL *********
C          AJUNK=ACCUM(IDT(I),X(I)*.001D0,Y(I)*.001D0,HTICE(I),
C     &                        SLOPN(1,I),
C     &                        PSURF(I),AMASS(9),TEMP(I))
          IF(IDT(I).EQ.24) THEN
            AJUNK=AFUNCT(TIME, 24, AMASS,
     &                   HTICE(I), DEPB(I), PSURF(I),
     &                   SLOPN(1,I),X(I),Y(I),TEMP(I))
          ELSEIF(IDT(I).EQ.20) THEN
            AJUNK=AFUNCT(TIME, 20, AMASS,
     &                   HTICE(I), DEPB(I), PSURF(I),
     &                   SLOPN(1,I),X(I),Y(I),TEMP(I))
          ENDIF
          IF(TFLAG) TEMP(I)=TSORIG(I)
C
C ....... IF SURFACE TEMP GREATER THAN 0.0, SET IT TO ZERO...
          IF(TEMP(I).GT.0.0) TEMP(I)=0.0D0
C
C
          TEMPA(1,I)=TEMP(I)
          IF(HTICE(I).GT.SEALEV) THEN
            THICK=HTICE(I)-DEPB(I)
          ELSE
            THICK=0.D0
          ENDIF
C ....... PRESSURE MELTING POINT AS FUNCTION OF DEPTH (THICKNESS)
C ....... TMELT=-8.7E-4*THICK
          TMELT=PMP(THICK)
C ...............................................................
C ----- A SIMPLE 1D STEADY STATE TEMPERATURE SOLVER
          IF(THICK.LE.1.) THEN
            INOICE=INOICE+1
            IF(DEPB(I).LT.SEALEV) THEN
              TBOT=0.0D0
              TBED(I)=TBOT
              TSURF=TEMP(I)
              TDT=(TBOT-TSURF)/(MMAX-1)
              DO J=1,MMAX
C                TEMPA(J,I)=0.0
                TEMPA(J,I)=TBOT-(MMAX-J+1)*TDT
              ENDDO
            ELSE
              TBOT=TEMP(I)
              TBED(I)=TBOT
              DO J=1,MMAX
                TEMPA(J,I)=TBOT
              ENDDO
            ENDIF        
            TTEMP=MIN(TMELT,TEMP(I))
C ***************THIS IS THE OLD ICE HARDNESS FUNCTION*********
C            IF(TTEMP.GE.-10.) THEN
C              AEFF=A0*EXP(-B0*TTEMP)
C            ELSE
C              AEFF=A1*EXP(-B1*TTEMP)
C            ENDIF
C **********THIS IS THE EISMINT ICE HARDNESS FUNCTION**********
            AEFF=HARDNESS(TTEMP)
            CALL TYPETEST(ITYPE(I),FRACT(I),AFUDGE(I),FLOWA(I),HTICE(I),
     &                    DEPB(I),SLDGB(I),WTHICK(I),STHRESH,
     &                    TBOT,TMELT,IMELTED,IFROZEN,AEFF)
          ELSE  
            IWITHICE=IWITHICE+1
            IF(TEMP(I).GE.0.0 .OR. .NOT.ISICE(HTICE(I),DEPB(I))) THEN
              DO J=1,MMAX
                TEMPA(J,I)=TEMP(I)
              ENDDO
              TBOT=TEMP(I)
C              AEFF=A1
C ****************EISMINT STUFF*************************
               AEFF=HARDNESS(0.D0)
            ELSE
              CALL COLTEMP1(TIME,I,NUMNP,AMASS,ADOT(I),TEMPA(1,I),
     &                      THICK,TMELT,TBOT,AEFF,DT,IPLSTRT,IPLT,
     &                      SLOPN(1,I),IMELT,BMELT(I),WTHICK(I),
     &                      STHRESHB,GEOFLUX(I),FRACT(I),SLDGB(I),
     &                      ITYPE(I),CONTRIB(I))
              CALL MODBMELT(BMELT(I),THICK,SLOPN(1,I))
              BMELT(I) = (CONV*BMELTOLD(I) + BMELT(I))/(CONV+1.D0)
              BMELTOLD(I)=BMELT(I)
C             IFLAG=1
C             BMELT(I)=-0.001
            ENDIF
            TBAVG=TBAVG+TBOT
            NTAVG=NTAVG+1
            CALL TYPETEST(ITYPE(I),FRACT(I),AFUDGE(I),FLOWA(I),HTICE(I),
     &                    DEPB(I),SLDGB(I),WTHICK(I),STHRESH,
     &                    TBOT,TMELT,IMELTED,IFROZEN,AEFF)
C -----------------
            TMIN=MIN(TMIN,TBOT)
            TMAX=MAX(TMAX,TBOT)
            IF(TBOT.GE.0.0) THEN
              IC1=IC1+1
            ELSE
              IC2=IC2+1
            ENDIF
            TBED(I)=TBOT
          ENDIF
        ELSE
          IFIXED=IFIXED+1
        ENDIF
C ... OUTPUT FOR .STF FILE
C        IF(IFLAG.EQ.0) THEN
C          DO J=1,MMAX
C            WRITE(66,*) REAL(X(I)),REAL(Y(I)),REAL(ZZZ(J)),
C     &                  'MISSING'
C          ENDDO
C        ELSE
C          DO J=1,MMAX
C            WRITE(66,*) REAL(X(I)),REAL(Y(I)),REAL(ZZZ(J)),
C     &                 REAL(TEMPA(J,I))
C          ENDDO
C        ENDIF
      ENDDO
      IF(KSTART.EQ.0) THEN
        KSTART=1
        HTOGG=HSAVE
      ENDIF
      IF(NTAVG.GT.0) THEN
        TBAVG=TBAVG/NTAVG
      ENDIF
      IF(ISTART.EQ.0) THEN
        DT=DTSAVE
        ISTART=1
      ENDIF
      PFIXED=DBLE(100*IFIXED)/NUMNP
      PNOICE=DBLE(100*INOICE)/NUMNP
      PWITHICE=DBLE(100*IWITHICE)/NUMNP
      PMELTED=DBLE(100*IC1)/NUMNP
      PFROZEN=DBLE(100*IC2)/NUMNP
      PSLIDE=DBLE(100*IMELTED)/NUMNP
      PNOSLIDE=DBLE(100*IFROZEN)/NUMNP
      AWI=PWITHICE*GRIDAREA/100.D0
      AWS=PMELTED*GRIDAREA/100.D0
      IF(AWI.NE.0.) THEN
        PWS=(AWS/AWI)*100.D0
      ELSE
        PWS=0.0D0
      ENDIF
      IF(IOTOGG) THEN
        WRITE(LIST(IPAGE+1),9) 'NO TEMP CALC',IFIXED,PFIXED
        WRITE(LIST(IPAGE+2),9) 'NO ICE      ',INOICE,PNOICE
        WRITE(LIST(IPAGE+3),9) 'WITH ICE   ',IWITHICE,PWITHICE,
     &                           AWI
        WRITE(LIST(IPAGE+4),9) ' MELTED     ',IC1,PMELTED,
     &                           AWS,PWS
        WRITE(LIST(IPAGE+5),9) ' FROZEN     ',IC2,PFROZEN
        WRITE(LIST(IPAGE+6),9) ' SLIDING    ',IMELTED,PSLIDE
        WRITE(LIST(IPAGE+7),9) ' NO SLIDING ',IFROZEN,PNOSLIDE
9     FORMAT(1X,A12,I6,F9.4,1PE13.6,0PF9.4)
        WRITE(LIST(IPAGE+8),*) REAL(TMIN),REAL(TMAX),REAL(TBAVG)
        IPAGE=IPAGE+8
        WRITE(79,*) 'TIME=',TIME
        WRITE(79,9) 'NO TEMP CALC',IFIXED,PFIXED
        WRITE(79,9) 'NO ICE      ',INOICE,PNOICE
        WRITE(79,9) 'WITH ICE   ',IWITHICE,PWITHICE,
     &                           AWI
        WRITE(79,9) ' MELTED     ',IC1,PMELTED,
     &                           AWS,PWS
        WRITE(79,9) ' FROZEN     ',IC2,PFROZEN
        WRITE(79,9) ' SLIDING    ',IMELTED,PSLIDE
        WRITE(79,9) ' NO SLIDING ',IFROZEN,PNOSLIDE
        WRITE(79,*) 'TMIN,TMAX,TAVG',REAL(TMIN),REAL(TMAX),REAL(TBAVG)
      ENDIF
      ASDF=TMIN
C      CALL WRITEMP(NUMNP)
      IF(IC1+IC2.GT.0) THEN
        RMAF=DBLE(IC1)/DBLE(IC1+IC2)
      ELSE
        RMAF=0.0D0
      ENDIF
C ++++++++ EISMINT STUFF FOR SLIDING EXPERIMENT +++++++++++++++++++++++++++++
C      CALL PAYNEOUT1(TIME,NUMNP,NUMEL,KX,X,Y,HTICE,DEPB,
C     &                    TEMPA,VEL,FLOWA,
C     &                    VOL,AREA,RMAF)
      IF(IPLSTRT.EQ.1) CALL GRSTOP1
      RETURN
      END
C====================================================
      SUBROUTINE BACKSUB(N,DDD,CCC,BBB,XXX)                             
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION BBB(N),CCC(N),DDD(N),XXX(N)                      
      XXX(N)=BBB(N)/DDD(N)                                              
      N1=N-1                                                            
      DO II=1,N1                                                      
        I=N-II                                                            
        XXX(I)=(BBB(I)-CCC(I)*XXX(I+1))/DDD(I)                          
      ENDDO                                                          
      RETURN                                                            
      END                                                               
C====================================================
C THIS IS THE NEW VERSION .................
      SUBROUTINE COLTEMP1(TIME,IPT,NUMNP,AMASS,ADOT,TTTT,
     &                    THICK,TMELT,TBOT,AEFF,DELT,IPLSTRT,IPLT,
     &                    SLOPN,IMELT,BMELT,WTHICK,STHRESHB,SIGMALI,
     &                    FRACT,SLDGB,ITYPE,CONTRIB)
C******************************************
C THIS SOLVES FOR TIME DEPENDENT VERTICAL TEMPERATURE PROFILE
C WITH SOLVED FOR BASAL TEMPERATURE
C IF BASAL TEMP > PRESSURE MELTING, WILL
C RESOLVE WITH BASAL BC FIXED AT PRESSURE MELTING
C AND CALCULATE BASAL MELTING
C******************************************
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      PARAMETER(MMAX=40,BIG=1.D20)
      DIMENSION AMASS(11),DTDT(MMAX)
      DIMENSION XXX(MMAX),TTTT(MMAX),COND(MMAX),HEAT(MMAX),RHO(MMAX)
      DIMENSION QQQ(MMAX),SOURCE(MMAX),DELX(MMAX)
      DIMENSION SAVED(MMAX)
      DIMENSION RKX(MMAX),FFFF(MMAX),AAAD(MMAX),AAAU(MMAX),AAAL(MMAX)
      DIMENSION CCCD(MMAX),CCCU(MMAX),CCCL(MMAX)
      DIMENSION WWWW(MMAX),TOLD(MMAX),AT(MMAX),UC(MMAX),DDATA(9)
      DIMENSION U(MMAX),TTTT0(MMAX)
      DIMENSION LMAP(MMAX)
      DIMENSION RJUNK(MMAX),RHOC(MMAX),TNEW(MMAX)
      COMMON /HEAT/ HTOGG
      LOGICAL HTOGG
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      logical stopflag
      common /stop/ stopflag
      SIGMAL=-SIGMALI
      SIGADJ=-SIGMALI
C      SIGMAL=-AMASS(11)
      PG = 0.089866D0*0.3816d0
      THIRD=1.D0/3.D0
      IF(.FALSE.) THEN
        PRINT *,'IN COLTEMP1 BEFORE'
        PRINT *,' IPT = ',IPT
        PRINT *,' NUMNP = ',NUMNP
        PRINT *,' AMASS = ',AMASS(11)
        PRINT *,' TTTT = ',TTTT(1)
        PRINT *,' THICK = ',THICK
        PRINT *,' TBOT = ',TBOT
        PRINT *,' AEFF = ',AEFF
        PRINT *,' DELT = ',DELT
        PRINT *,' IPLSTRT = ',IPLSTRT
        PRINT *,' IPLT = ',IPLT
        PRINT *,' WTHICK = ',WTHICK
      ENDIF
      NPT=MMAX
      NDT=10
      WWWW(1)=0.D0
C ... QQQ IS FOR HORIZONTAL ADVECTION (MULT BY U (HORIZONTAL VELOCITY AT
C ... EACH DEPTH)
C ... SOURCE IS FOR INTERNAL (SHEAR) HEAT GENERATION
      QQQS=0.0D0
      WWWW(1)=ADOT
C
      LTOT=0
      DO I=1,NPT/2
        LMAP(I)=I
c       LMAP(I)=1
        LTOT=LTOT+LMAP(I)
      ENDDO
      DO I=NPT/2+1,NPT
        LMAP(I)=LMAP(I-1)-1
c       LMAP(I)=1
        LTOT=LTOT+LMAP(I)
      ENDDO
      XXX(1)=THICK
      XXX(NPT)=0.D0
      DX=(XXX(NPT)-XXX(1))/DBLE(LTOT)
      QQQ(1)=0.0D0
      TOLD(1)=TTTT(1)
      DO I=2,NPT
        XXX(I)=XXX(I-1)+LMAP(I-1)*DX
        QQQ(I)=0.D0
        TOLD(I)=TTTT(I)
      ENDDO
C
      TOLER=.001D0
      FACTOR=1.D0
      XMAX=-1D30
      XMIN=1D30
9     CONTINUE
      DO I=1,NPT
        XMAX=MAX(XMAX,XXX(I))
        XMIN=MIN(XMIN,XXX(I))
        DEPTH=XXX(1)-XXX(I)
        RHO(I)=DENSITY(DEPTH)
        COND(I)=CONDUCT(TTTT(I),RHO(I))
        HEAT(I)=SPHEAT(TTTT(I))
C*************USING FIXED VALUES OF RHO,COND,HEAT*****************
C         RHO(I)=DENSITY(1000.D0)
C         COND(I)=CONDUCT(-30.D0,RHO(I))
C         HEAT(I)=SPHEAT(-30.D0)
C*************USING FIXED VALUES OF RHO,COND,HEAT*****************
C ..... EVERYTHING IN ONE TERM...
C        RKX(I)=COND(I)/RHO(I)/HEAT(I)
C        RHOC(I)=1D0
C ..... ALA PATERSON PG 224
        RKX(I)=COND(I)
        RHOC(I)=RHO(I)*HEAT(I)
        IF(I.GT.1) DELX(I-1)=-(XXX(I)-XXX(I-1))
        AAAD(I)=0.D0
        AAAU(I)=0.D0
        AAAL(I)=0.D0
        CCCD(I)=0.D0
        CCCU(I)=0.D0
        CCCL(I)=0.D0
        FFFF(I)=0.D0
        U(I)=0.0D0
        SOURCE(I)=0.0D0
        SAVED(I)=0.0D0
      ENDDO
C           
C          
      IF(IPLSTRT.EQ.0 .AND. IPLT.EQ.1) THEN
        CALL GRSTRT(600,600)
C        CALL WINDOW(REAL(TMIN),REAL(TMAX),REAL(XMIN),REAL(XMAX))
        CALL WINDOW(REAL(TMIN),REAL(TMAX),0.,5000.)
        IPLSTRT=1
      ENDIF
      TMIN=-60.D0
      TMAX=60.D0
      IF(IPLT.EQ.1) THEN
C        CALL WINDOW(REAL(TMIN),REAL(TMAX),REAL(XMIN),REAL(XMAX))
C        CALL WINDOW(REAL(TMIN),REAL(TMAX),0.,5000.)
C        CALL NEWPAG
        CALL LINCLR(1)
        CALL MOVE(0.,REAL(XMIN))
        CALL DRAW(0.,REAL(XMAX))
        CALL MOVE(REAL(TMELT),REAL(XMIN))
        CALL MOVE(REAL(TMELT),REAL(XMIN))
        CALL DRAW(REAL(TMELT),REAL(XMAX))
        CALL MRKCLR(2)
      ENDIF
      WWWW(NPT)=0.D0
C --- VERTICAL VELOCITY DISTRIBUTION: LINEAR
      SLOPE=(WWWW(NPT)-WWWW(1))/(XXX(NPT)-XXX(1))
      DO I=2,NPT-1
        WWWW(I)=WWWW(1)+SLOPE*(XXX(I)-XXX(1))
      ENDDO
C      DO I=1,NPT
C        QQQ(I)=QQQS
C      ENDDO
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
        CCCU(I)=.5D0*THIRD*DELX(I)*RHOC(I)
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
C          CALL WINDOW(REAL(TMIN),REAL(TMAX),REAL(XMIN),REAL(XMAX))
C          CALL WINDOW(REAL(TMIN),REAL(TMAX),0.,5000.)
        ENDIF
C.....COME IN HERE FOR ITERATION ON VERTICAL VELOCITY
700   CONTINUE   
C
C ... HEAT GENERERATION SECTION ...................................
C
C ..... CALCULATE THE INTERNAL HEAT GENERATION DUE TO SHEAR
        DO I=1,NPT
          SOURCE(I)=0.D0
          SAVED(I)=0.D0
        ENDDO
        UFLOW=0.0
        DO I=2,NPT
          TMID=0.5D0*(TTTT(I-1)+TTTT(I))
c         TMID=0.5D0*(TTTT(1)+TTTT(NPT))
          ATI=HARDNESS(TMID)
          DEPTH=XXX(1)-0.5D0*(XXX(I-1)+XXX(I))
          DLX=(XXX(I-1)-XXX(I))
C ....... PATERSON IS 2*EPSILON-DOT*TAU RATE OF HEAT/UNIT VOLUME (TIMES ELEMENT
C         THICKNESS TO GET TOTAL HEAT DEPOSITED IN ELEMENT)
C         FOLLOWING IS 1/2 THAT, TO BE DEPOSITED AT TOP AND BOTTOM OF ELEMENT
C ....... FOLLOWING IS IN BAR/YR
          HTFLOW=DLX*((PG*DEPTH*SLOPN)**4/ATI**3)
          UFLOW=UFLOW+DLX*(PG*DEPTH*SLOPN/ATI)**3
C ....... FOLLOWING CONVERTS TO CAL/YR
          HTFLOW=HTFLOW*2.392D4
C ....... FOLLOWING TURNS DOWN HEAT WHERE SLIDING IS OCCURING
C          HTFLOW=HTFLOW*(1.D0-FRACT)
C
          SOURCE(I-1)=SOURCE(I-1)+HTFLOW
          SOURCE(I)=SOURCE(I)+HTFLOW
          SAVED(I-1)=SAVED(I-1)+HTFLOW
          SAVED(I)=SAVED(I)+HTFLOW
        ENDDO
        do i=1,npt
          saved(i)=saved(i)/(xxx(1)-xxx(npt))
          saved(i)=log10(saved(i))
         enddo
C
C ..... CALCULATE THE INTERNAL HEAT GENERATION DUE TO BASAL SLIDING
        HTSLID=(PG*THICK*SLOPN)**3/SLDGB**2
C ..... FOLLOWING CONVERTS TO CAL/YR
        HTSLID=HTSLID*2.392D4
C ..... SCALE BY WTHICK WITH RLUB FUNCTION (EXPERIMENTAL)
        IF(ITYPE.EQ.8) THEN
          HTSLID=HTSLID*RLUB(WTHICK)
          USLIDE=(PG*THICK*SLOPN/SLDGB)**2*RLUB(WTHICK)
          IF(USLIDE+UFLOW.GT.0.0) THEN
            FRACT=USLIDE/(USLIDE+UFLOW)
          ELSE
            FRACT=0
          ENDIF
C         IF(USLIDE.GT.0) PRINT *,REAL(USLIDE),REAL(UFLOW),REAL(FRACT)
        ENDIF
C ..... REDUCE BY 10 (ALA RAYMOND...)
        HTSLID=HTSLID/10.D0
C
C ..... THIS DEPOSITS ALL THE HEAT AT THE BED RATHER THAN WITHIN THE COLUMN
        HTFLOW=0.D0
        IF(.true.) THEN
          DO I=1,NPT
            HTFLOW=HTFLOW+SOURCE(I)
            SOURCE(I)=0.0D0
          ENDDO 
        ENDIF
C
        IF(.FALSE. .AND. WTHICK.GT.0) PRINT *,REAL(-CONTRIB/SIGMAL),
     &          REAL(-HTSLID/SIGMAL),
     &          REAL(FRACT),
     &          REAL(WTHICK)
C
        IF(ITYPE.NE.8) THEN
C ....... USE HEAT PROPORTIONAL TO FRACT ...
          CONTRIB=(1.D0-FRACT)*HTFLOW+FRACT*HTSLID
        ELSEIF(ITYPE.EQ.8) THEN
C ....... USE HEAT SUM OF FLOW AND SLIDING ...
          CONTRIB=HTFLOW+HTSLID
        ELSEIF(.FALSE.) THEN
C ....... USE HEAT FROM FLOW PROPORTIONAL TO FRACT ...
          CONTRIB=(1.D0-FRACT)*HTFLOW
        ENDIF
        IF(HTOGG) THEN
          SIGADJ=SIGMAL-CONTRIB
        ELSE
C ....... THIS TURNS OFF INTERNAL HEAT GENERATION ...........
          SIGADJ=SIGMAL
          HTSLID=0
          HTFLOW=0
          CONTRIB=0
          DO I=1,NPT
            SOURCE(I)=0.0D0
          ENDDO 
        ENDIF
C ........................................................................
C          PRINT *,REAL(CONTRIB),REAL(SLOPN),
C     &            REAL(XXX(1)-XXX(NPT)),
C     &            REAL(PG*SLOPN*(XXX(1)-XXX(NPT)))
C          PRINT *,'CONTRIB',CONTRIB,SIGMAL
C        DO I=1,NPT
C          PRINT *,I,SOURCE(I)
C        ENDDO
C ... END INTERNAL HEAT GENERATION SECTION
C ----- FORM LOAD VECTOR 
        FFFF(1)=.5D0*U(1)*QQQ(1)*DELX(1)+SOURCE(1)                             
        FFFF(NPT)=.5D0*U(NPT-1)*QQQ(NPT-1)*DELX(NPT-1)+
     &             SOURCE(NPT)-SIGADJ          
        P1=FFFF(1)                                                      
        PN=FFFF(NPT)
        DO I=2,NPT-1    
          FFFF(I)=.5D0*(U(I)*QQQ(I)*DELX(I)+
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
        AAAD(1)=RKX(1)/DELX(1)-RHOC(1)*WWWW(1)*.5D0                               
        AAAD(NPT)=RKX(NPT-1)/DELX(NPT-1)+RHOC(NPT-1)*WWWW(NPT)*.5D0                   
        DO I=2,NPT-1   
          AAAD(I)=RKX(I-1)/DELX(I-1)+RKX(I)/DELX(I)                     
          AAAD(I)=AAAD(I)+.5D0*(RHOC(I-1)*WWWW(I-1)-RHOC(I)*WWWW(I))                        
        ENDDO   
        DO I=1,NPT-1                                                
          AAAU(I)=-RKX(I)/DELX(I)                                       
          AAAL(I)=AAAU(I)-RHOC(I)*WWWW(I)*.5D0                                   
          AAAU(I)=AAAU(I)+RHOC(I)*WWWW(I)*.5D0 
        ENDDO 
C SAVE VALUES BEFORE EFFECTIVE MATRIX IS FORMED
        FFFFNPT=FFFF(NPT)+SIGADJ                                                    
        AAADNPT=AAAD(NPT)
        AAALNPT=AAAL(NPT-1)
        CCCDNPT=CCCD(NPT)
        CCCLNPT=CCCL(NPT-1)
C        PRINT *,0D0,CCCD(1),CCCU(1)
C        DO I=2,NPT-1
C          PRINT *,CCCL(I-1),CCCD(I),CCCU(I)
C        ENDDO
C        PRINT *,CCCL(NPT-1),CCCD(NPT)
C        PAUSE
C        PRINT *,0D0,AAAD(1),AAAU(1)
C        DO I=2,NPT-1
C          PRINT *,AAAL(I-1),AAAD(I),AAAU(I)
C        ENDDO
C        PRINT *,AAAL(NPT-1),AAAD(NPT)
C        PAUSE
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
C SAVE VALUES AFTER EFFECTIVE MATRIX IS FORMED                                         
        IF(.FALSE.) THEN
          FFFFNPT=FFFF(NPT)+SIGMAL                                                    
          AAADNPT=AAAD(NPT)
          AAALNPT=AAAL(NPT-1)
        ENDIF
C
C       DO I=1,NPT                                                  
C         WRITE(7,*) 'EFF F',I,REAL(FFFF(I)),REAL(RJUNK(I))                                     
C       ENDDO  
C ----- FIX BOUNDARY CONDITION FOR NODE 1                                     
        RK11=AAAD(1)                                                    
        RK12=AAAU(1)                                                    
        RK21=AAAL(1)                                                    
        AAAD(1)=1.D0
        AAAU(1)=0.D0
        AAAL(1)=0.D0
        FFFF(1)=TTTT(1)                                                 
        FFFF(2)=FFFF(2)-RK21*FFFF(1)                                    
C ----- FIX BOUNDARY CONDITION FOR NODE NPT
        RKN1=AAAD(NPT)                                                  
        RKN2=AAAU(NPT-1)                                                
        RK2N=AAAL(NPT-1)    
        IPASS=1     
C ..... ISET=0 IS FROZEN BED
        ISET=0       
C                                
C ----- THE MATRIX SOLUTION
C ..... COPY TEMPERATURE INTO TNEW FOR SOLUTION ...
        DO KKK=1,NPT
          TNEW(KKK)=TTTT(KKK)
        ENDDO
C       PRINT *,' TTTT ',(REAL(TNEW(KKK)),KKK=1,NPT)
        CALL TRI(NPT,AAAL,AAAD,AAAU,FFFF,TNEW) 
C       PRINT *,' AFTER ',(REAL(TNEW(KKK)),KKK=1,NPT)
c       if(IPT.eq.1377) then
c         WRITE(LIST(IPAGE+1),*) IPT,xxx(1)-xxx(npt),tnew(npt)
c         IPAGE=IPAGE+1
c       endif
c       if(IPT.eq.1377 .and. mod(nint(TIME),10000).eq.0) then
c         do i=1,npt
c           write(2,*) saved(i),xxx(i)
c           write(3,*) tnew(i),xxx(i)
c         enddo
c         write(2,*) -99999.,0.
c         write(2,19) 'S',TIME
c         write(3,*) -99999.,0.
c         write(3,19) 'T',TIME
c       endif 
C       PAUSE
C
        IF(TNEW(NPT).LT.TMELT .AND. WTHICK.GT.STHRESHB) THEN
C ....... CHECK FOR BED BELOW MELTING POINT WITH WATER PRESENT ...
c         PRINT *,REAL(TNEW(NPT)),REAL(TMELT),' BELOW MELTING POINT ',
c    &           IPT,wthick,sthreshb
          DO KKK=1,NPT
            TNEW(KKK)=TTTT(KKK)
          ENDDO
          TNEW(NPT)=TMELT
C ....... ISET=1 MEANS BED IS FREEZING, AT PMP...
          ISET=1
          AAAD(NPT)=1.D0
          AAAU(NPT-1)=0.D0
          AAAL(NPT-1)=0.D0
          FFFF(NPT)=TNEW(NPT)                                             
          FFFF(NPT-1)=FFFF(NPT-1)-RKN2*FFFF(NPT)                          
          CALL TRI(NPT,AAAL,AAAD,AAAU,FFFF,TNEW) 
        ELSEIF(TNEW(NPT).GT.TMELT .AND. ISET.EQ.0) THEN
C ....... CHECK FOR BED ABOVE MELTING POINT ...
c         PRINT *,REAL(TNEW(NPT)),REAL(TMELT),' ABOVE MELTING POINT ',
c    &           IPT
c         WRITE(LIST(IPAGE+1),*) IPT,real(tnew(npt)),' above PMP'
c         IPAGE=IPAGE+1
          DO KKK=1,NPT
            TNEW(KKK)=TTTT(KKK)
          ENDDO
          TNEW(NPT)=TMELT
C ....... ISET=1 MEANS BED IS MELTING, AT PMP...
          ISET=1
          AAAD(NPT)=1.D0
          AAAU(NPT-1)=0.D0
          AAAL(NPT-1)=0.D0
          FFFF(NPT)=TNEW(NPT)                                             
          FFFF(NPT-1)=FFFF(NPT-1)-RKN2*FFFF(NPT)                          
c         open(1,file='plot-source.data')
c         do i=1,npt
c           write(1,*) saved(i),xxx(i)
c         enddo
c         write(1,*) -99999.,0.
c         write(1,19) 'source',real(slopn)
19    format(1x,a,1pg9.2)
c         close(1)
c         open(1,file='plot-temp.data')
c         do i=1,npt
c           write(1,*) tnew(i),xxx(i)
c           write(2,*) saved(i),xxx(i)
c           write(3,*) tnew(i),xxx(i)
c         enddo
c         write(1,*) -99999.,0.
c         write(1,19) 'temp before',real(slopn)
c         write(2,*) -99999.,0.
c         write(2,19) 'S',TIME
c         write(3,*) -99999.,0.
c         write(3,19) 'temp before',TIME
          CALL TRI(NPT,AAAL,AAAD,AAAU,FFFF,TNEW) 
c         do i=1,npt
c           write(1,*) tnew(i),xxx(i)
c           write(3,*) tnew(i),xxx(i)
c         enddo
c         write(1,*) -99999.,0.
c         write(1,19) 'temp after',real(slopn)
c         write(3,*) -99999.,0.
c         write(3,19) 'temp after',real(slopn)
c         close(1)
c         stopflag=.true.
        ENDIF
C ..... ACCEPT TNEW SOLUTION INTO TTTT ....
        DO KKK=1,NPT
          IF(DELT.GT.0) THEN
            DTDT(KKK)=(TNEW(KKK)-TTTT(KKK))/DELT
          ENDIF
          TTTT(KKK)=TNEW(KKK)
        ENDDO
C
C ..... EXPERIMENT TO "SMOOTH" TEMPS ...
C        CALL SMOOTH(TTTT,NPT)
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
        IF(.FALSE.) THEN
          IF(ISCAN.GT.0..AND.IOTOGG) THEN 
            WRITE(LIST(IPAGE+1),*) ISCAN,' POLYTHERMAL AT ',IPT
            IPAGE=IPAGE+1
          ENDIF
        ENDIF
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
C         RHO(I)=DENSITY(1000.D0)
C         COND(I)=CONDUCT(-30.D0,RHO(I))
C         HEAT(I)=SPHEAT(-30.D0)
C********************************************
C ..... EVERYTHING IN ONE TERM...
          RKXI=COND(I)/RHO(I)/HEAT(I)                                     
C ..... ALA PATERSON PG 224
          RKXI=COND(I)
          RLATENT=8.D4*RHO(I)
          RKX(I)=((FACTOR-1.D0)*RKX(I)+RKXI)/FACTOR                      
        ENDDO                                                        
C        SIGMA0=P1-((RK11-CCCD(1))*TTTT(1)+
C     &         RK12*TTTT(2)-CCCU(1)*TTTT0(2))
C        SIGMA0=-SIGMA0/DELT                                             
C        SIGMALC=PN-((RKN1-CCCD(NPT))*TTTT(NPT)+
C     &         RKN2*TTTT(NPT-1)-CCCU(NPT-1)*TTTT0(NPT-1))   
C        DTEMP=(TTTT(NPT-1)-TTTT(NPT))/(XXX(NPT-1)-XXX(NPT))
C        SIGMABC=-RKX(NPT-1)*DTEMP
        IF(DELT.GT.0) THEN
C ....... THIS IS FULL TIME-DEPENDENT, INCLUDING CAP-MATRIX
          TERM1=-AAALNPT*TTTT(NPT-1)-AAADNPT*TTTT(NPT)
          TERM2=-CCCLNPT*DTDT(NPT-1)-CCCDNPT*DTDT(NPT)
          SIGMABC=FFFFNPT+TERM1+TERM2
        ELSE
C ....... FOLLOWING STEADY-STATE (DELT=0)
          TERM1=-AAALNPT*TTTT(NPT-1)-AAADNPT*TTTT(NPT)
          TERM2=0.D0
          SIGMABC=FFFFNPT+TERM1+TERM2
        ENDIF
C---------------------------------------
C WHAT IS THIS ????
C
C        SIGMALC=-SIGMALC/DELT                                             
C--------------------------------------------
C ..... ONLY CALC BMELT IF TOGGLE IS ON...
        IF(IMELT.EQ.1) THEN
          IF(ISET.EQ.0) THEN
            BMELT=0.0D0
          ELSE
C           IF(ABS((SIGADJ-SIGMABC)/SIGMAL).GT.1E-6) THEN
              BMELT=-(SIGADJ-SIGMABC)/RKX(NPT-1)
              BMELT=-(SIGADJ-SIGMABC)/RLATENT
C           ENDIF
          ENDIF
C ... LET IT PLAY, WHATEVER IT MAKES... MAKE MODIFICATIONS IN MODBMELT
C ....... DISABLED ITERATION ON VERTICAL VELOCITY BASED ON BASAL MELT RATE
          IF(.FALSE.) THEN
            IF(ABS(WWWW(NPT)-BMELT).GT.TOLER) THEN
C             PRINT *,IPT,REAL(WWWW(NPT)),REAL(BMELT),
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
C ....... DISABLEDTO HERE
        ENDIF
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
C      IF(ADOT.LT.0) CALL WAIT(100000)
C ++++++++ EISMINT STUFF FOR SLIDING EXPERIMENT +++++++++++++++++++++++++++++
C      IF (IPT.EQ.1861) CALL PAYNEOUT2(TIME,IPT,MMAX,TTTT,XXX,AT,WWWW)
C      IF (IPT.EQ.2471) CALL PAYNEOUT2(TIME,IPT,MMAX,TTTT,XXX,AT,WWWW)
C      IF (IPT.EQ.2776) CALL PAYNEOUT2(TIME,IPT,MMAX,TTTT,XXX,AT,WWWW)
C      IF (IPT.EQ.2109) CALL PAYNEOUT2(TIME,IPT,MMAX,TTTT,XXX,AT,WWWW)
C      IF (IPT.EQ.2719) CALL PAYNEOUT2(TIME,IPT,MMAX,TTTT,XXX,AT,WWWW)
C ========================================================================
c      IF(IPT.EQ.228 .OR.IPT.EQ.227) THEN
c      IF(BMELT.gt.1.) THEN
      if(.false.) then
7     FORMAT(1X,A,1P3E15.6)
8     FORMAT(1X,A,3I13)
        write(7,*)
        WRITE(7,7) 'IN COLTEMP1 AFTER'
        WRITE(7,8) ' IPT     = ',IPT
        WRITE(7,7) ' NUMNP   = ',NUMNP
        WRITE(7,7) ' AMASS   = ',AMASS(11)
        WRITE(7,7) ' TTTT    = ',TTTT(1)
        WRITE(7,*) ' THICK   = ',THICK
        WRITE(7,7) ' TBOT    = ',(TBOT),(TTTT(NPT)),(TMELT)
        WRITE(7,*) ' AEFF    = ',AEFF
        WRITE(7,7) ' DELT    = ',DELT
        WRITE(7,8) ' IPLSTRT = ',IPLSTRT
        WRITE(7,8) ' IPLT    = ',IPLT
        WRITE(7,7) ' WTHICK  = ',WTHICK
        WRITE(7,7) ' BMELT   = ',BMELT
        WRITE(7,7) ' SIGMAL  = ',SIGMAL,SIGADJ
        WRITE(7,7) ' SIGMABC = ',SIGMABC,SIGMABC/SIGADJ,
     &                          (SIGADJ-SIGMABC)
        WRITE(7,7) ' CONTRIB = ',CONTRIB
        WRITE(7,7) ' HTSLID  = ',HTSLID,PG*THICK*SLOPN
        WRITE(7,*) ' HTFLOW  = ',HTFLOW
        WRITE(7,8) ' ITYPE   = ',ITYPE
        WRITE(7,7) ' FFFFNPT = ',(FFFFNPT),(TERM1),(TERM2)
        write(7,*) ' slope   = ',SLOPN
        write(7,*) ' temps   = '
        write(7,*)   TTTT    
        write(7,*) ' xxxs    = '
        write(7,*)   XXX    
        write(7,*)
        bmelt=0.d0
        pause
      ENDIF
      RETURN
C      IF(IPLT.EQ.1) THEN
C        CALL WAIT(10000)
C        CALL GRSTOP1
C      ENDIF
      END         
C*************************************************
      SUBROUTINE VELO(TIME,NPTS,DH,X,T,RHO,AT,U,UC,
     &                BMELT,DDATA,AEFF)
      PARAMETER(MMAX=40)
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE UOLD,UCOLD,UBAROLD,AVGOLD
      DIMENSION DDATA(9)
      DIMENSION X(NPTS),T(NPTS),RHO(NPTS),AT(NPTS)
      DIMENSION U(NPTS),UC(NPTS),UOLD(MMAX),UCOLD(MMAX)
      DATA UOLD /MMAX*0D0/, UCOLD /MMAX*0D0/, AVG /0D0/,UBAROLD /0D0/
      IF(NPTS.GT.MMAX) THEN
        PRINT *,'INCREASE MMAX IN VELO, MMAX=',MMAX
        STOP
      ENDIF
      THIRD=1D0/3D0
      A0=1.830357D0
      B0=7.738778D-2
      A1=2.633955D0
      B1=4.0990233D-2
      CALL HARD(A0,B0,A1,B1)
      DO I=1,NPTS
C ***************THIS IS THE OLD ICE HARDNESS FUNCTION*********
C        IF(T(I).GE.-10.) THEN
C          AT(I)=A0*EXP(-B0*T(I))
C        ELSE
C          AT(I)=A1*EXP(-B1*T(I))
C        ENDIF
C **********THIS IS THE EISMINT ICE HARDNESS FUNCTION**********
        AT(I)=HARDNESS(T(I))
C
      ENDDO
      TOP=X(1)
      G=9.81D-5*0.3816d0
      AVGAT=0.D0
      AVGT=0.D0
      DO I=NPTS-1,1,-1
        DX=X(I+1)-X(I)
        AVGAT=AVGAT-0.5D0*(AT(I+1)+AT(I))*DX
        AVGT=AVGT-0.5D0*(T(I+1)+T(I))*DX
      ENDDO
      AVGAT=AVGAT/(X(1)-X(NPTS))
      AVGT=AVGT/(X(1)-X(NPTS))
      DDATA(5)=AVGT
C      PRINT *,'AVERAGE A=',AVGAT
      AVGR=0.D0
      DO I=NPTS-1,1,-1
        DX=X(I+1)-X(I)
        AVGR=AVGR-0.5D0*(RHO(I+1)+RHO(I))*DX
      ENDDO
      AVGR=AVGR/(X(1)-X(NPTS))
C      PRINT *,'AVERAGE RHO=',AVGR
      UBAR=.4D0*(AVGR*G*DH/AVGAT)**3*(X(1)-X(NPTS))**4
C      PRINT *,'UBAR=',UBAR
      U(NPTS)=0.D0
      DO I=NPTS-1,1,-1
        SUM=U(I+1)
        V1=2.D0*(RHO(I+1)*G*DH*(TOP-X(I+1))/AT(I+1))**3
        V2=2.D0*(RHO(I)*G*DH*(TOP-X(I))/AT(I))**3
        DX=X(I+1)-X(I)
        U(I)=U(I+1)-0.5D0*(V1+V2)*DX
      ENDDO
      AVG=0.D0
      DO I=NPTS-1,1,-1
        DX=X(I+1)-X(I)
        AVG=AVG-0.5D0*(U(I+1)+U(I))*DX
      ENDDO
      AVG=AVG/(X(1)-X(NPTS))
      DDATA(6)=AVG
      AEFF=(.4D0*(AVGR*G*DH)**3*(X(1)-X(NPTS))**4/AVG)**THIRD
C      PRINT 1002,TIME,REAL(AVGT),REAL(AEFF),REAL(AVG),
C     &                REAL(T(NPTS)),REAL(BMELT)
1002  FORMAT(F10.0,2F10.3,F10.2,F10.3,F10.6)
      DDATA(7)=AEFF
C      PRINT *,'A-EFFECTIVE=',AEFF,3.17E-23/AEFF**3
      C=.5D0*(AVGR*G*DH/AEFF)**3
      THICK=X(1)-X(NPTS)
      DO I=1,NPTS
        UC(I)=C*(THICK**4-(TOP-X(I))**4)
      ENDDO
      YMIN=X(NPTS)
      YMAX=X(1)
      UMIN=1D30
      UMAX=-1D30
C       PRINT 1001,'X','T','A','U','UC'
1001  FORMAT(1X,5A13)
      DO I=1,NPTS
        UMIN=MIN(UMIN,U(I))
        UMAX=MAX(UMAX,U(I))
        UMIN=MIN(UMIN,UC(I))
        UMAX=MAX(UMAX,UC(I))
C         PRINT 1000,X(I),T(I),AT(I),U(I),UC(I)
1000    FORMAT(1X,1P5G13.6)
      ENDDO
C      PRINT *,UMIN,UMAX,YMIN,YMAX
C      CALL GRSTRT(400,400)
C      CALL WINDOW(REAL(0.),REAL(500.),REAL(YMIN),REAL(YMAX))
C      CALL LINCLR(0)
C      CALL MOVE(REAL(UOLD(1)),REAL(X(1)))
C      DO I=2,NPTS
C        CALL DRAW(REAL(UOLD(I)),REAL(X(I)))
C      ENDDO
C      CALL MOVE(REAL(AVGOLD),REAL(YMIN))
C      CALL DRAW(REAL(AVGOLD),REAL(YMAX))
C      CALL LINCLR(0)
C      CALL MOVE(REAL(UCOLD(1)),REAL(X(1)))
C      DO I=2,NPTS
C        CALL DRAW(REAL(UCOLD(I)),REAL(X(I)))
C      ENDDO
C      CALL MOVE(REAL(UBAROLD),REAL(YMIN))
C      CALL DRAW(REAL(UBAROLD),REAL(YMAX))
C      CALL LINCLR(1)
C      CALL MOVE(REAL(U(1)),REAL(X(1)))
C      DO I=2,NPTS
C        CALL DRAW(REAL(U(I)),REAL(X(I)))
C      ENDDO
C      CALL MOVE(REAL(AVG),REAL(YMIN))
C      CALL DRAW(REAL(AVG),REAL(YMAX))
C      CALL LINCLR(2)
C      CALL MOVE(REAL(UC(1)),REAL(X(1)))
C      DO I=2,NPTS
C        CALL DRAW(REAL(UC(I)),REAL(X(I)))
C      ENDDO
C      CALL MOVE(REAL(UBAR),REAL(YMIN))
C      CALL DRAW(REAL(UBAR),REAL(YMAX))
      AVGOLD=AVG
C      UBAROLD=UBAR
      DO I=1,NPTS
        UOLD(I)=U(I)
        UCOLD(I)=UC(I)
      ENDDO
C      CALL GRSTOP1
      END
C*****************************************
      SUBROUTINE HARD(B0P,SLPBP,BB0P,SLPBBP)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NNNN=11)
      SAVE IPASS,B0,SLPB,BB0,SLPBB
      save t,a,b
      REAL*8 T(NNNN),A(NNNN),B(NNNN)
c      DATA T /0.d0,-5.d0,-10.d0,-15.d0,-20.d0,-25.d0,-30.d0,
c    &      -35.d0,-40.d0,-45.d0,-50.d0/
c     DATA A /5.3D-15,1.7D-15,5.2D-16,3.1D-16,1.8D-16,
c    &        1.0D-16,5.4D-17,2.9D-17,1.5D-17,7.7D-18,3.8D-18/
      DATA IPASS /0/
      IF(IPASS.EQ.0) THEN
        IPASS=1
        T(1)=0.d0
        T(2)=-5.d0
        T(3)= -10.d0
        T(4)= -15.d0
        T(5)= -20.d0
        T(6)= -25.d0
        T(7)= -30.d0 
        t(8)=-35.d0
        t(9)= -40.d0
        t(10)= -45.d0
        t(11)= -50.d0 
c
        A(1)=5.3D-15
        A(2)=1.7D-15
        A(3)=5.2D-16
        A(4)=3.1D-16
        A(5)=1.8D-16
        A(6)=1.0D-16
        A(7)=5.4D-17
        A(8)=2.9D-17
        A(9)=1.5D-17
        A(10)=7.7D-18
        A(11)=3.8D-18
        CONVER=1.6D-2/5.2D-16
        DO I=1,NNNN
          TEMP=A(I)*CONVER
          B(I)=1.D0/TEMP
          B(I)=B(I)**(1.D0/3.D0)
C          PRINT *,T(I),A(I),TEMP,B(I)
        ENDDO
        B0=B(1)
        SLPB=-(LOG(B(1))-LOG(B(3)) )/(T(1)-T(3))
        SLPBB=-(LOG(B(3))-LOG(B(11)))/(T(3)-T(11))
        RLNBB0=LOG(B(3))+SLPBB*T(3)
        BB0=EXP(RLNBB0)
      ENDIF
C          PRINT *,'B0,BB0',B0,BB0
C          PRINT *,'R,RR  ',SLPB,SLPBB
      B0P=B0
      SLPBP=SLPB
      BB0P=BB0
      SLPBBP=SLPBB
      END
C*********************************************
      SUBROUTINE WAIT(I)
      DO N=1,I
        X=SIN(REAL(I))
      ENDDO
      END                                                      
C*************************************************
      SUBROUTINE TRI(N,AAA,DDD,CCC,BBB,XXX) 
      PARAMETER(MMAX=40)                            
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION AAA(N),BBB(N),CCC(N),DDD(N),XXX(N)                      
      DIMENSION AAL(MMAX),BBL(MMAX),CCL(MMAX),DDL(MMAX)
      IF(N.GT.MMAX) THEN
        PRINT *,'INCREASE MMAX IN SUB. TRI',MMAX
        STOP
      ENDIF
      DO I=1,N-1
        AAL(I)=AAA(I)
        BBL(I)=BBB(I)
C      PRINT *,I,BBL(I),BBB(I)
        CCL(I)=CCC(I)
        DDL(I)=DDD(I)
      ENDDO
      DDL(N)=DDD(N)
      BBL(N)=BBB(N)
      DO I=2,N                                                        
        XMULT=AAA(I-1)/DDL(I-1)                                         
C        PRINT *,' AI-1,DI-1,XM',
C     &          I,REAL(AAA(I-1)),REAL(DDL(I-1)),REAL(XMULT)
C        PRINT *,' DI,BI-1,CI-1 ',REAL(DDL(I)),REAL(BBL(I-1)),
C     &           REAL(CCL(I-1))
        DDL(I)=DDL(I)-XMULT*CCL(I-1)
C        PRINT *,'BI BEFORE ',I,BBL(I)
        BBL(I)=BBL(I)-XMULT*BBL(I-1)                                    
C        PRINT *,'BI AFTER  ',I,BBL(I)
      ENDDO                                                          
      XXX(N)=BBL(N)/DDL(N)                                              
      N1=N-1                                                            
      DO II=1,N1                                                      
        I=N-II                 
C        PRINT *,I,REAL(XXX(I+1)),REAL(DDL(I)),REAL(BBL(I)),REAL(CCL(I))                                           
        XXX(I)=(BBL(I)-CCL(I)*XXX(I+1))/DDL(I)                          
      ENDDO      
      RETURN                                                            
      END                                                               
C*************************************************
      FUNCTION DENSITY(DEPTH)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(RHOI=910.D0,RHOS=350.D0,CRHO=0.0374D0)
        IF (DEPTH .EQ. 0.) THEN
          DENSITY=RHOS
        ELSEIF ((DEPTH.GT.0.).AND.(DEPTH.LE.55.)) THEN
          DENSITY=RHOI-((RHOI-RHOS)*EXP(-CRHO*DEPTH))
        ELSEIF ((DEPTH.GT.55.) .AND.(DEPTH.LT.65.)) THEN
          DENSITY=7.D0*DEPTH+455.D0
        ELSEIF (DEPTH .GE. 65.) THEN
          DENSITY=RHOI
        ENDIF
      END
C*************************************************
      FUNCTION CONDUCT(TTTT,RHO)
      IMPLICIT REAL*8(A-H,O-Z)
C ... UNITS: CALORIES/M//YR/DEGREE
      IF(TTTT.GT.0.) THEN
        TTT=0.D0
      ELSE
        TTT=TTTT
      ENDIF
      RRR=RHO*.001D0
      CONDUCT=3.1536D06*
     &          (( 4.20102D0-0.00754145D0*TTT)
     &          -(26.48767D0-0.048779D0  *TTT)*RRR
     &          +(57.31865D0-0.141127D0  *TTT)*RRR**2
     &          -(29.55156D0-0.053672D0  *TTT)*RRR**3)
      END
C*************************************************
      FUNCTION SPHEAT(TTTT)
      IMPLICIT REAL*8(A-H,O-Z)
C ... UNITS: CALORIES/KG/DEGREE
      IF(TTTT.GT.0.) THEN
        TTT=0.D0
      ELSE
        TTT=TTTT
      ENDIF
      SPHEAT=1000.D0*(0.494D0+0.00138D0*TTT)                          
      END
C********************************************************
      SUBROUTINE WRITEMP(NUMNP)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      PARAMETER(MMAX=40)
      CHARACTER*1 CHAR
      COMMON /TEMPERATURE/ TEMP(MMAX,NMAX)
      PRINT *,'IN WRITEMP',NUMNP
      PRINT *,'   TO WRITE OUT BACKUP OF TEMPERATURE'
      PRINT *,'   INPUT Y'
      READ(*,1002) CHAR
      WRITE(99,1002) CHAR
1002  FORMAT(A1)
      IF(CHAR.EQ.'y' .OR. CHAR.EQ.'Y') THEN
        PRINT *,'0:BINARY, 1:ASCII'
        READ(*,*) IREAD
        WRITE(99,*) IREAD
        NPT=MMAX
        REWIND 44
        IF(IREAD.EQ.1) THEN
          DO I=1,NUMNP
            DO J=1,NPT
              WRITE(44,*) I,J,TEMP(J,I)
            ENDDO
          ENDDO
        ELSE
          WRITE(44) NUMNP,NPT
          WRITE(44) ((TEMP(J,I),J=1,NPT),I=1,NUMNP)
        ENDIF
      ENDIF
      END
C-------------------------------
      SUBROUTINE SMOOTH(X,N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N)
      DO I=2,N-1
        X(I)=(X(I-1)+X(I)+X(I+1))/3.D0
      ENDDO
      END
C-------------------------------
      FUNCTION HARDNESS(TTT)
      IMPLICIT REAL*8(A-H,O-Z)
c     PARAMETER(FUD=2**(1D0/3D0))
      data fud /1.25d0/
      DATA RRR /8.314D0/
      DATA Q1 /60000.D0/
      DATA A1 /1.14D-5/
      DATA Q2 /139000.D0/
      DATA A2 /5.47D10/
C *****THIS IS 3 FOR GREENLAND EISMINT EXPERIMENT, 1 FOR SLIDING
      DATA EH1 /3.D0/
C      DATA EH1 /1.D0/
C THIS IS 2**(1/3) WHICH FOR SOME REASON MULTIPLIES ENHANCEMENT FACTOR ?
      EH=EH1*FUD
C
      IF(TTT.LT.0.) THEN
        TEMP=TTT+273.15D0
      ELSE
        TEMP=273.15D0
      ENDIF
      IF(TEMP.LE.263.15D0) THEN
        AAA=A1*EXP(-Q1/RRR/TEMP)
      ELSE
        AAA=A2*EXP(-Q2/RRR/TEMP)
      ENDIF
      BBB=EH*AAA*1D15
      BBB=1.D0/BBB
      HARDNESS=BBB**(1.D0/3.D0)
      END
C==========================================
      SUBROUTINE COLUMNZ(NPT,HTICE,DEPB,ZZZ)
      PARAMETER(MMAX=40)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ZZZ(NPT),LMAP(MMAX)
      THICK=HTICE-DEPB
      LTOT=0
      DO I=1,NPT/2
        LMAP(I)=I
C        LMAP(I)=1
        LTOT=LTOT+LMAP(I)
      ENDDO
      DO I=NPT/2+1,NPT
        LMAP(I)=LMAP(I-1)-1
C        LMAP(I)=1
        LTOT=LTOT+LMAP(I)
      ENDDO
      ZZZ(1)=THICK
      ZZZ(NPT)=0.D0
      DX=(ZZZ(NPT)-ZZZ(1))/DBLE(LTOT)
      DO I=2,NPT
        ZZZ(I)=ZZZ(I-1)+LMAP(I-1)*DX
      ENDDO
      DO I=1,NPT
        ZZZ(I)=ZZZ(I)+DEPB
      ENDDO
      END
C==========================================
      FUNCTION PMP(THICK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(SLOPE=0.0074d0*0.089866d0*0.3816d0)
c     PMP=-8.7D-4*THICK
      PMP=-SLOPE*THICK
      END
C==========================================
      SUBROUTINE CHKFLOT(HTICEI,DEPBI,FRACTI)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
      DATA RHOI /0.917D0/, RHOW /1.092D0/
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      IF(DEPBI.LT.SEALEV) THEN
        FLOT=(1.D0-RATDEN)*(DEPBI-SEALEV)
        IF(HTICEI.LE.FLOT) FRACTI=1.D0
      ELSE
        IF(HTICEI.LE.DEPBI) FRACTI=1.D0        
      ENDIF
      END
C==========================================
      LOGICAL FUNCTION ISICE(HTICEI,DEPBI)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
      DATA RHOI /0.917D0/, RHOW /1.092D0/
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      ISICE=.TRUE.
      IF(DEPBI.LT.SEALEV) THEN
        FLOT=(1.D0-RATDEN)*(DEPBI-SEALEV)
        IF(HTICEI.LT.FLOT) ISICE=.FALSE.
      ELSE
        IF(HTICEI.LE.DEPBI) ISICE=.FALSE.
      ENDIF
      END
C-----------------------------------
      SUBROUTINE MODBMELT(BMELTI,THICK,SLOPN)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(.false.) RETURN
C
C ... MODIFY BMELT ......................
      IF(.FALSE.) THEN
C ..... ARBITRARY MULTIPLICATIVE FACTOR 
        BMELTI=BMELTI*0.01
      ENDIF
      IF(.TRUE.) THEN
C ..... THIS IS BASIC OFFSET .02 MM 
        BMELTI=BMELTI-2D-5
      ENDIF
      IF(.false. .AND. BMELTI.LT.0) THEN
C ..... REDUCE FREEZING BY FACTOR OF 10
C       PRINT *,' REDUCE FREEZING BY FACTOR OF 10 ',BMELTI
        BMELTI=BMELTI/10.D0
      ENDIF
      IF(.false. .AND. BMELTI.LT.0) THEN
C ..... LIMIT FREEZE RATE TO 10 CM/YR
C       PRINT *,' LIMIT FREEZE RATE TO 10 CM/YR ',BMELTI
        BMELTI=MAX(-0.1D0,BMELTI)
      ENDIF
      IF(.FALSE. .AND. BMELTI.LT.0) THEN
C ..... NO FREEZING AT ALL
C       PRINT *,' NO FREEZING AT ALL ',BMELTI
        BMELTI=0.D0
      ENDIF
      BLIMIT=0.01D0
      IF(.false. .AND. BMELTI.GT. BLIMIT) THEN
C ..... LIMIT MELT RATE TO 1 CM/YR
C       PRINT *,' LIMIT MELT RATE TO 1 CM/YR ',BMELTI
        BMELTI=MIN(BLIMIT,BMELTI)
      ENDIF
c      BLIMIT=9.1d0
      BLIMIT=0.1d0
      IF(.true. .AND. BMELTI.GT. BLIMIT) THEN
C ..... remove MELT RATEs over 10 CM/YR
c       PRINT *,' REMOVE MELT RATE OVER 10 CM/YR ',BMELTI
        BMELTI=0.d0
      ENDIF
      TLIMIT=10.D0
      IF(.FALSE. .AND. THICK.LT. TLIMIT) THEN
C ..... LIMIT BMELT TO ZERO IN THIN ICE REGIONS
        PRINT *,' MELT RATE TO 0 FOR THIN ICE',BMELTI,THICK
        BMELTI=0.D0
      ENDIF
      TLIMIT=10.D0
      SLIMIT=1E-2
      IF(.false. .AND. SLOPN.GT. SLIMIT .AND. THICK.LT.TLIMIT) THEN
C ..... LIMIT BMELT TO ZERO IN STEEP-THIN ICE REGIONS
        PRINT *,' MELT RATE 0 FOR STEEP/THINICE',BMELTI,SLOPN,THICK
        BMELTI=0.D0
      ENDIF
C ... MODIFY BMELT ......................
C
      END
C-----------------------------------
      SUBROUTINE TYPETEST(ITYPEI,FRACTI,AFUDGEI,FLOWAI,HTICEI,
     &                    DEPBI,SLDGBI,WTHICKI,STHRESH,
     &                    TBOT,TMELT,IMELTED,IFROZEN,AEFF)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL ISICE
      LOGICAL ISICE
C ... FOLLOWING IS AMOUNT OF WATER (M) NECESSARY FOR FULL SLIDING...
      DATA COUPLE /0.1D0/
      IF(ITYPEI.EQ.1) THEN
C ..... MODIFY FRACTION ONLY
        IF(TBOT.GE.TMELT) THEN
C         PRINT *,I,TBOT,TEMPI,THICK
          FRACTI=1.D0
          IMELTED=IMELTED+1
        ELSE
          FRACTI=0.D0
          IFROZEN=IFROZEN+1
        ENDIF
      ELSEIF(ITYPEI.EQ.2) THEN
C ..... MODIFY FLOW ONLY
        AEFF=AEFF*AFUDGEI
        FLOWAI=AEFF
      ELSEIF(ITYPEI.EQ.3) THEN
C ..... MODIFY FRACTION
        IF(TBOT.GE.TMELT) THEN
C         PRINT *,I,TBOT,TEMPI,THICK
          FRACTI=1.D0
          IMELTED=IMELTED+1
        ELSE
          FRACTI=0.D0
          CALL CHKFLOT(HTICEI,DEPBI,FRACTI)
          IF(FRACTI.EQ.0.) THEN
            IFROZEN=IFROZEN+1
          ELSE
            IMELTED=IMELTED+1
          ENDIF
        ENDIF
C ..... MODIFY FLOW
        AEFF=AEFF*AFUDGEI
        FLOWAI=AEFF
      ELSEIF(ITYPEI.EQ.4) THEN
C ..... MODIFY FLOW
        AEFF=AEFF*AFUDGEI
        FLOWAI=AEFF
C ..... MODIFY SLIDING
        SLDGBI=AEFF/200.D0
      ELSEIF(ITYPEI.EQ.5) THEN
C ..... MODIFY FRACTION
        IF(TBOT.GE.TMELT) THEN
C         PRINT *,I,TBOT,TEMPI,THICK
          FRACTI=1.D0
          IMELTED=IMELTED+1
        ELSE
          FRACTI=0.D0
          IFROZEN=IFROZEN+1
        ENDIF
C ..... MODIFY FLOW
        AEFF=AEFF*AFUDGEI
        FLOWAI=AEFF
C ..... MODIFY SLIDING
        SLDGBI=AEFF/200.D0
      ELSEIF(ITYPEI.EQ.6) THEN
C ..... MODIFY AFUDGE TO FIT FLOW
        AFUDGEI=FLOWAI/AEFF
      ELSEIF(ITYPEI.EQ.7) THEN
C ..... MODIFY FRACTION BY WATER THICKNESS
        IF(WTHICKI.GE.STHRESH .AND.
     &    ISICE(HTICEI,DEPBI)) THEN
C          FRACTI=MIN(1.D0,MAX(0.D0,(WTHICKI-STHRESH)/COUPLE))
          FRACTI=MIN(1.D0,MAX(0.D0,(WTHICKI)/COUPLE))
          IMELTED=IMELTED+1
        ELSE
          FRACTI=0.D0
          CALL CHKFLOT(HTICEI,DEPBI,FRACTI)
          IF(FRACTI.EQ.0.) THEN
            IFROZEN=IFROZEN+1
          ELSE
            IMELTED=IMELTED+1
          ENDIF
        ENDIF
C ..... MODIFY FLOW
        AEFF=AEFF*AFUDGEI
        FLOWAI=AEFF
C AS EXPERIMENT, MODIFY SLIDING TO MATCH FRACTION.
CC      IF(ISICE(HTICEI,DEPBI)) THEN
C         SLDGBI=0.133333*(1.5-FRACTI)
CC      ENDIF
      ELSEIF(ITYPEI.EQ.8) THEN
C ..... MODIFY FRACTION BY WATER THICKNESS
        IF(WTHICKI.GE.STHRESH .AND.
     &    ISICE(HTICEI,DEPBI)) THEN
          IMELTED=IMELTED+1
        ELSE
          FRACTI=0.D0
          CALL CHKFLOT(HTICEI,DEPBI,FRACTI)
C          IF(FRACTI.EQ.1) WTHICKI=1.
          IF(FRACTI.EQ.0.) THEN
            IFROZEN=IFROZEN+1
          ELSE
            IMELTED=IMELTED+1
          ENDIF
        ENDIF
C ..... MODIFY FLOW
        AEFF=AEFF*AFUDGEI
        FLOWAI=AEFF
      ENDIF
      END

     
