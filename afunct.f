      SUBROUTINE SCENARIO(TIME)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /EXPER/ TPERIOD,TINIT,TVSTART,TVFINAL,IEXPER
      write(*,*) ' 1-        mass(7) coef      0  10.5 k  -900 -2100'
      write(*,*) ' 2-<  70 k mass(7) coef1     0    20 k     0 -2300'
      write(*,*) '   >  70 k mass(7) coef1    70    90 k -2300     0'
      write(*,*) ' 3-        mass(7) coef      0  20.5 k     0 -2300'
      write(*,*) ' 3-        mass(8) coef    -25 -14.5 k  5600  3800'
      write(*,*) ' 4-<20.5 k mass(7) coef1     0  20.5 k     0 -2300'
      write(*,*) '   <61.5 k mass(7) coef1  20.5    41 k -2300 -1000'
      write(*,*) '   >61.5 k mass(7) coef1  61.5    82 k -2300     0'
      write(*,*) ' 5-<20.5 k mass(7) coef1     0  20.5 k -2000     0'
      write(*,*) '   <61.5 k mass(7) coef1  20.5    41 k     0 -3000'
      write(*,*) '   >61.5 k mass(7) coef1  61.5    82 k     0 -2000'
      write(*,*) ' 6-        mass(7) coef     -3   6.4 k     0 -2000'
      write(*,*) ' 7-        mass(7) coef     -3     7 k -1000 -4000'
c
      write(*,*) ' 8-        tnsl    coef  -9.55 10.55 k     0   -20'
      write(*,*) ' 9-        tnsl    coef      0     1 k    -6   -20'
      write(*,*) '10- northern hemisphere ice sheet cycle (smooth)'
      write(*,*) '   <100 k  tnsl    coef      0    20 k     0   -20'
      write(*,*) '   <120 k  tnsl    coef    100   120 k   -20     0'
      write(*,*) '   <128 k  tnsl    coef    120   128 k     0    -5'
      write(*,*) '   >128 k  tnsl = -5'
      write(*,*) '11- northern hemisphere ice sheet cycle (linear)'
      write(*,*) '   < 40 k  tnsl    coef      0    20 k    -7   -20'
      write(*,*) '   < 41 k  tnsl    coef     40    41 k   -20     0'
      write(*,*) '   < 51 k  tnsl    coef     51    52 k     0    -5'
      write(*,*) '   > 51 k  tnsl = -5'
c      write(*,*) '12         tnsl    coef      0    5.4 k    .5    -4'
c      write(*,*) '           tnsl    coef    5.4   12.3 k    -4     3'
c      write(*,*) '           tnsl = 3'
      write(*,*) '            tnsl    coef      0   10 k    -5    -6'
      write(*,*) '           tnsl    coef    10    20 k   -6    -3'
      write(*,*) '           tnsl = -3'
      write(*,*) '13- 20,000 year cycle of tnsl'
      write(*,*) '   < 10 k  tnsl    coef      0    10 k    -7   -20'
c      write(*,*) '           acom    coef      0    10 k    -5  -9.5'
c      write(*,*) '           mass(7) coef      0    10 k -2000 -2600'
      write(*,*) '   < 20 k  tnsl    coef     10    20 k   -20    -0'
c      write(*,*) '           acom    coef     10    20 k  -9.5    -5'
c      write(*,*) '           mass(7) coef      0    10 k -2600 -2000'
      write(*,*) '   > 20 k  tnsl    =    -7'
c      write(*,*) '           acom    =    -5'
      write(*,*) '           mass(7) = -2000'
      write(*,*) '14- abrupt change'
      write(*,*) '   <  5 k  tnsl    coef      0     5 k    -7   -20'
      write(*,*) '   < 15 k  tnsl = -20'
      write(*,*) '   < 20 k  tnsl    coef     15    20 k   -20     0'
      write(*,*) '   > 20 k  tnsl = 0'
      write(*,*) '15- very abrupt change'
      write(*,*) '   < 30 k  tnsl = -20'
      write(*,*) '   > 30 k  tnsl =   0'
      write(*,*) ' 16- linear deglaciation scenario'
      write(*,*) '           tnsl    coefl     0     5 k   -20     0'
      write(*,*) ' 17- read in TLIST.DATA or user defined'
      write(*,*) ' 18- 20,000 year cycle of tnsl for superglacial'
      write(*,*) '   < 10 k  tnsl    coef      0    10 k    -7   -23'
c      write(*,*) '           acom    coef      0    10 k    -5  -9.5'
c      write(*,*) '           mass(7) coef      0    10 k -2000 -2600'
      write(*,*) '   < 20 k  tnsl    coef     10    20 k   -23    -0'
c      write(*,*) '           acom    coef     10    20 k  -9.5    -5'
c      write(*,*) '           mass(7) coef      0    10 k -2600 -2000'
      write(*,*) '   > 20 k  tnsl    =    -0'
c      write(*,*) '           acom    =    -5'
c      write(*,*) '           mass(7) = -2000'
      write(*,*) ' 19- Abrupt climat change companion to scenario 18'
      write(*,*) '   <  1 k tnsl =   0'
      write(*,*) '   >  1 k tnsl = -23'
      write(*,*) '   > 15 k tnsl =   0'
      write(*,*) ' 20- Anatarctic scenario'
      write(*,*) '   < 10 k  tnsl    coef      0    10 k   -14   -23'
      write(*,*) '   < 20 k  tnsl    coef     10    20 k   -23   -14'
      write(*,*) '   > 20 k  tnsl    =    -0'
      write(*,*) ' 21- scandinavian scenario'
      write(*,*) '     40 k period tnsl   coef 0    20 k   -2   -18'
      write(*,*) ' 22- antarctic scenario'
      write(*,*) '     20 k period tnsl   coef 0    5 k   -20   -8'
      write(*,*) ' 23- user defined scenario'
      write(*,*) '   input period,tstart,tfinal '
      write(*,*) ' 24- eismint scenario'
      
      READ(*,*) IEXPER
      write(99,*) IEXPER
      IF(IEXPER.EQ.23) THEN
        TINIT=TIME
        WRITE(*,*) 'INPUT PERIOD,INITIAL TNSL, AND FINAL TNSL'
        READ(*,*) TPERIOD,TVSTART,TVFINAL
        WRITE(99,*) TPERIOD,TVSTART,TVFINAL
      ENDIF
      END
      FUNCTION AFUNCT(TIME,IDOT,AMASS,HT,BD,PSURF,SLOPN,X,Y,TS)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
      COMMON /SNOW/ SNOLIN,SNO(20)
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      COMMON /TTNSL/ TTIME(NSMAX),TLIST(NSMAX),NTNSL
      COMMON /EXPER/ TPERIOD,TINIT,TVSTART,TVFINAL,IEXPER
      DIMENSION AMASS(11),SLOPN(4)
      data pi /3.14159265358979d0/
      FPC(S,E)=(-7.044897460d+02)*S**2+( 3.286502080d+01)*S+
     1  (-8.363957700d-08)*E**2+( 2.803513780d-04)*E+
     1  (-8.280271290d-02)
      FPX(S,E)=(-1.011654790d+03)*S**2+( 5.757759090d+01)*S+
     1  (-1.945462600d-07)*E**2+( 6.485106420d-04)*E+
     1  (-1.714566350d-01)
      FSX(S,E)=( 1.613718720d+02)*S**2+(-3.064633370d+00)*S+
     1  (-2.274655910d-07)*E**2+( 9.221490470d-04)*E+
     1  (-3.658012150d-01)
      FTM(S,E)=(-2.475505130d+03)*S**2+( 7.388871770d+01)*S+
     1  (-7.129281700d-07)*E**2+( 3.180227710d-03)*E+
     1  (-6.778311130d-01)
      FSM(S,E)=(-1.272729690d+04)*S**2+( 1.028782810d+02)*S+
     1  ( 1.790857030d-07)*E**2+(-1.209567770d-03)*E+
     1  ( 2.244830130d+00)
C CHOOSE COMPENENT OF THE SLOPE, (1) - SQRT(DHDX**2+DHDY**2)
C                                (2) - DHDX
C                                (3) - DHDY
C                                (4) - DHDX*SIN(THETA)+DHDY*COS(THETA)
      SL=SLOPN(1)
C
C FOLLOWING SECTION IS QUITE VARIABLE. IT CONTAINS CALLS TO COEF AND
C COEF1 TO VARY VARIOUS MASS BALANCE PARAMETERS IN A SINUSOIDAL FASHION,
C AND WILL BE DIFFERENT FOR DIFFERENT EXPERIMENTS.
C
C     SNOLIN=0.
c     
      TNSL=AMASS(9)
      IF(IEXPER.EQ.1) THEN
C
        AMASS(7)=COEF(TIME,0.,10500.,-900.,-2100.)
C
      ELSEIF(IEXPER.EQ.2) THEN
C
        IF(TIME.LT.70000.) THEN
          AMASS(7)=COEF1(TIME,0.,20000.,0.,-2300.)
        ELSE
          AMASS(7)=COEF1(TIME,70000.,90000.,-2300.,0.)
        ENDIF
C
      ELSEIF(IEXPER.EQ.3) THEN
C
        DELSNNP=COEF(TIME,0.,20500.,0.,-2300.)
        AMASS(7)=DELSNNP
        DELSNEQ=COEF(TIME,-25000.,-14500.,5600.,3800.)
        AMASS(8)=(DELSNEQ-DELSNNP)*1d-7
C
      ELSEIF(IEXPER.EQ.4) THEN
C
        IF(TIME.LE.20500.) THEN
          AMASS(7)=COEF1(TIME,0.,20500.,0.,-2300.)
        ELSEIF(TIME.LE.61500.) THEN
          AMASS(7)=COEF(TIME,20500.,41000.,-2300.,-1000.)
        ELSE
          AMASS(7)=COEF1(TIME,61500.,82000.,-2300.,0.)
        ENDIF
C
      ELSEIF(IEXPER.EQ.5) THEN
C
        IF(TIME.LE.20500.) THEN
          AMASS(7)=COEF1(TIME,0.,20500.,-2000.,0.)
        ELSEIF(TIME.LE.61500.) THEN
          AMASS(7)=COEF(TIME,20500.,41000.,0.,-3000.)
        ELSE
          AMASS(7)=COEF1(TIME,61500.,82000.,0.,-2000.)
        ENDIF
C
      ELSEIF(IEXPER.EQ.6) THEN
C
        AMASS(7)=COEF(TIME,-3600.,6400.,0.,-2000.)
C
      ELSEIF(IEXPER.EQ.7) THEN
C
        AMASS(7)=COEF(TIME,-3000.,7000.,-1000.,-4000.)
C
c variable stuff, to modify TNSL
      ELSEIF(IEXPER.EQ.8) THEN
        TNSL=COEF(TIME,-9550.,10950.,0.,-20.)

      ELSEIF(IEXPER.EQ.9) THEN
        TNSL=COEF(TIME,0.,1000.,-6.,-20.)

      ELSEIF(IEXPER.EQ.10) THEN
C EXPERIMENT FOR NORTHERN HEMISHHERE ICE SHEET CYCLE (SMOOTH VARIATION)
        IF(TIME.LT.100000.) THEN
          TNSL=COEF(TIME,0.,20000.,-7.,-20.)
        ELSEIF(TIME.LT.120000.) THEN
          TNSL=COEF1(TIME,100000.,120000.,-20.,0.)
        ELSEIF(TIME.LT.128000.) THEN
          TNSL=COEF1(TIME,120000.,128000.,0.,-5.)
        ELSE
          TNSL=-5.d0
        ENDIF
C END OF EXPERIMENT FOR NORTHERN HEMISPHERE ICE SHEET CYCLE

      ELSEIF(IEXPER.EQ.11) THEN
C EXPERIMENT FOR NORTHERN HEMISHHERE ICE SHEET CYCLE (LINEAR VARIATION)
        IF(TIME.LT.40000.) THEN
          TNSL=COEFL(TIME,0.,20000.,-7.,-20.)
        ELSEIF(TIME.LT.41000.) THEN
          TNSL=COEFL(TIME,40000.,41000.,-20.,0.)
        ELSEIF(TIME.LT.51000.) THEN
          TNSL=COEFL(TIME,51000.,52000.,0.,-5.)
        ELSE
          TNSL=-5.d0
        ENDIF
C END OF EXPERIMENT FOR NORTHERN HEMISPHERE ICE SHEET CYCLE

      ELSEIF(IEXPER.EQ.12) THEN
C lake baikal experiment tnsl goes from 0 to -6 and back in 20000
c        IF(TIME.LT.5400.) THEN
c          TNSL=COEF(TIME,0.,5400.,.5,-4.)
c        ELSEIF(TIME.LT.12300.) THEN
c          TNSL=COEF(TIME,5400.,12300.,-4.,3.)
c        ELSE
c          TNSL=3.
c        ENDIF
        IF(TIME.LT.10000.) THEN
          TNSL=COEF(TIME,0.,10000.,-5.,-6.)
        ELSEIF(TIME.LT.20000.) THEN
          TNSL=COEF(TIME,10000.,20000.,-6.,-3.)
        ELSE
          TNSL=-3.d0
        ENDIF
c end lake baikal cycle

      ELSEIF(IEXPER.EQ.13) THEN
C smooth variation experiment, 20,000 yr sin cycle both tnsl and lapse
C rate     
        IF(TIME.LT.10000.) THEN
          TNSL=COEF(TIME,0.,10000.,-7.,-20.)
c          ACOM=COEF(TIME,0.,10000.,-5.,-9.5)
c          AMASS(7)=COEF(TIME,0.,10000.,-2000.,-2600.)
        ELSEIF(TIME.LT.20000.) THEN
          TNSL=COEF(TIME,10000.,20000.,-20.,0.)
c          ACOM=COEF(TIME,10000.,20000.,-9.5,-5.)
c          AMASS(7)=COEF(TIME,10000.,20000,-2600.,-2000.)
        ELSE
          TNSL=0.d0
c          ACOM=-5.d0
c          AMASS(7)=-2000.d0
        ENDIF
c end smooth cycle

      ELSEIF(IEXPER.EQ.14) THEN
c abrupt change experiment
        IF(TIME.LT.5000.) THEN
          TNSL=COEFL(TIME,0.,5000.,-7.,-20.)
        ELSEIF(TIME.LT.15000.) THEN
          TNSL=-20.d0
        ELSEIF(TIME.LT.20000.) THEN
          TNSL=COEFL(TIME,15000.,20000.,-20.,0.)
        ELSE
          TNSL=0.d0
        ENDIF
c end abrupt change experiment

      ELSEIF(IEXPER.EQ.15) THEN
c very abrupt climat change
        IF(TIME.LE.30001.) THEN
          TNSL=-20.d0
        ELSE
          TNSL=0.d0
        ENDIF
c end very abrupt clmate change

      ELSEIF(IEXPER.EQ.16) THEN
c begin linear deglaciation scenario
        TNSL=COEFL(TIME,0.,5000.,-20.,0.)
c end linear deglaciation scenario
      ELSEIF(IEXPER.EQ.17) THEN
c use read in list of tnsl
        IF(NTNSL.GT.0) TNSL=FINDTNSL(TIME)
        AMASS(9)=TNSL
      ELSEIF(IEXPER.EQ.18) THEN
C smooth variation experiment, 20,000 yr sin cycle both tnsl and lapse
C rate     
        IF(TIME.LT.10000.) THEN
          TNSL=COEF(TIME,0.,10000.,-7.,-23.)
c          ACOM=COEF(TIME,0.,10000.,-5.,-9.5)
c          AMASS(7)=COEF(TIME,0.,10000.,-2000.,-2600.)
        ELSEIF(TIME.LT.20000.) THEN
          TNSL=COEF(TIME,10000.,20000.,-23.,-0.)
c          ACOM=COEF(TIME,10000.,20000.,-9.5,-5.)
c          AMASS(7)=COEF(TIME,10000.,20000,-2600.,-2000.)
        ELSE
          TNSL=-0.d0
          ACOM=-5.d0
          AMASS(7)=-2000.d0
        ENDIF
c end smooth cycle

      ELSEIF(IEXPER.EQ.19) THEN
c abrupt climate change companion to scenario 18
        IF(TIME.LT.1000.) THEN
          TNSL=0.d0
c          ACOM=-5
c          AMASS(7)=-2000.
        ELSEIF(TIME.LT.15000.) THEN
          TNSL=-23.d0
c          ACOM=-9.5
c          AMASS(7)=-2600.
        ELSE
          TNSL=0.d0
c          ACOM=-5.
c          AMASS(7)=-2000.
        ENDIF
c end abrupt climate change

      ELSEIF(IEXPER.EQ.20) THEN
C antarrctic variation experiment, 20,000 yr sin cycle tnsl
C      
        IF(TIME.LT.10000.) THEN
          TNSL=COEF(TIME,0.,10000.,-14.,-23.)
        ELSEIF(TIME.LT.20000.) THEN
          TNSL=COEF(TIME,10000.,20000.,-23.,-14.)
        ELSE
          TNSL=-14.d0
        ENDIF
c end smooth cycle
      ELSEIF(IEXPER.EQ.21) THEN
C scandinavian variation experiment, 20,000 yr sin cycle tnsl
C      
        TNSL=COEF(TIME,0.,20000.,-2.,-18.)
c end smooth cycle
      ELSEIF(IEXPER.EQ.22) THEN
C antarctica variation experiment, 20,000 yr sin cycle tnsl
C      
        TNSL=COEF(TIME,-5000.,5000.,-20.,-8.)
c end smooth cycle
      ELSEIF(IEXPER.EQ.23) THEN
C user defined variation experiment.
C      TPERIOD yr sin cycle tnsl
C      TVSTART  initial time
C      TVFINAL final time
C      
        TNSL=COEF(TIME,REAL(TINIT),REAL(TINIT+TPERIOD),
     &                 REAL(TVSTART),REAL(TVFINAL))
        AMASS(9)=TNSL
c end smooth cycle
      ELSEIF(IEXPER.EQ.24) THEN
C eismint variation experiment, 20,000,40000 yr sin cycle acc
c this is also zone 24
C 
c following for uniform, constant accumulation rate
c          afunct=0.3     
c following for uniform accumulation 20000 yr cycle
c          afunct=0.3+0.2*sin(8.*atan(1.)*time/20000.)
c following for uniform accumulation 40000 yr cycle
c          afunct=0.3+0.2*sin(8.*atan(1.)*time/40000.)
c following for ablation zone experiment, steady state
c          distance in km
           dist=sqrt(x**2+y**2)*.001d0
c equilibrium line distance nonuniform accumulation with ablation zone
c following for constant eld
c           eld=450.d0
c following for eld 20000 yr cycle
          eld=450.d0+100.d0*sin(8.d0*atan(1.d0)*time/20000.d0)
c following for eld 40000 yr cycle
c          eld=450.d0+100.d0*sin(8.d0*atan(1.d0)*time/40000.d0)
c
           afunct=min(0.5d0,1d-2*(eld-dist))
c
          return
c end eismint cycle
      ENDIF

        
c end variable stuff of TNSL
C END VARIABLE SECTION
C
C BRANCH FOR VARIOUS MASS BALANCE RELATIONSHIPS
      IF(IDOT.EQ.0) THEN
        AFUNCT=0.0d0
        RETURN
      ENDIF
      GOTO(10,20,30,40,50,60,70,                     ! IDOT 1..7 straight lines
     &     200,210,220,230,240,                      ! IDOT 8..12 
     &     300,310,320,330,340,345,347,348,350,360,  ! IDOT 13..22
     &     370,380,                                  ! IDOT 23..24
     &     400,400,400,400,400,400,400),IDOT         ! IDOT 25..31 data-based
C
10    CONTINUE
c ... mars peak and pole....
c         AMARG=AMASS(1)
c         APEAK=AMASS(2)
c         ADOME=AMASS(3)
c         HE=AMASS(4)
c         HP=AMASS(5)
c         HM=AMASS(6)
          X0=AMASS(1)
          Y0=AMASS(2)
          R0=AMASS(3)
          HP=AMASS(4)
          SL=-HP/R0
          AFUNCT=HP+SL*SQRT((X*0.001-X0)**2+(Y*0.001-Y0)**2)
        return
      GOTO 100
C
20    CONTINUE
c ... mars peak and pole....
          AMARG=AMASS(1)
          APEAK=AMASS(2)
          ADOME=AMASS(3)
          HE=AMASS(4)
          HP=AMASS(5)
          HM=AMASS(6)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
30    CONTINUE
C **** =3 POLAR - B
          AMARG=0.0d0
          APEAK=1.d0
          ADOME=.10d0
          HE=0.d0
          HP=1250.d0
          HM=2500.d0
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
40    CONTINUE
C **** =4 POLAR - C
          AMARG=-1.0d0
          APEAK=1.5d0
          ADOME=.1d0
          HE=100.d0
          HP=1250.d0
          HM=2500.d0
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
50    CONTINUE
C **** IDOT=5 MARITIME - A
          AMARG=-2.0d0
          APEAK=2.d0
          ADOME=.15d0
          HE=1500.d0
          HP=1750.d0
          HM=2000.d0
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
60    CONTINUE
C **** IDOT=6 MARITIME - B
          AMARG=-2.0d0
          APEAK=2.d0
          ADOME=.15d0
          HE=1000.d0
          HP=1500.d0
          HM=2000.d0
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
70    CONTINUE
C **** IDOT=7 MARITIME - C
          AMARG=-3.0d0
          APEAK=2.d0
          ADOME=.15d0
          HE=1000.d0
          HP=1500.d0
          HM=2000.d0
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
100   CONTINUE
      AMASS(1)=AMARG
      AMASS(2)=APEAK
      AMASS(3)=ADOME
      AMASS(4)=HE
      AMASS(5)=HP
      AMASS(6)=HM
      HTT=HT-SNOLIN
      IF(HTT.GT.HM) THEN
        AFUNCT=ADOME
      ELSE
        IF(HTT.LE.0.) THEN
          AFUNCT=AMARG
        ELSE
          IF(HTT.GT.0. .AND. HTT.LE.HE) THEN
            A1=AMARG
            H1=0.d0
            A2=0.d0
            H2=HE
          ENDIF
          IF(HTT.GT.HE .AND. HTT.LE.HP) THEN
            A1=0.d0
            H1=HE
            A2=APEAK
            H2=HP
          ENDIF
          IF(HTT.GT.HP .AND. HTT.LE.HM) THEN
            A1=APEAK
            H1=HP
            A2=ADOME
            H2=HM
          ENDIF
          AFUNCT=(A1*(H2-HTT)-A2*(H1-HTT))/(H2-H1)
        ENDIF
      ENDIF
      IF(BD.LE.-5000.) AFUNCT=-5.d0
      RETURN
C
200   CONTINUE
C **** IDOT=8 SM (SUBPOLAR MARITIME) SLOPE DEPENDENT
          ADOME=.10d0
          SMAX=.01d0
          B0=-5.d0
          HEQ=700.d0
          B0DOT=B0/HEQ
          BNEQ=FSM(SMAX,HEQ)
          HADD=-BNEQ/B0DOT
          IF(HT.LE.HEQ) THEN
            BN=B0-B0DOT*HT
          ELSE IF(HT.LE.(HEQ+HADD)) THEN
            BN=FSM(SL,HT)+B0DOT*(HADD-HT+HEQ)
            IF(BN.LE.0.) BN=ADOME
          ELSE
            BN=FSM(SL,HT)
            IF(BN.LE.0.) BN=ADOME
          ENDIF
          AFUNCT=BN
      GOTO 1000
C
210   CONTINUE
C **** IDOT=9 TM (TEMPERATE MARITIME) SLOPE DEPENDENT
          ADOME=.10d0
          SMAX=.03d0
          B0=-8.d0
          HEQ=1125.d0
          B0DOT=B0/HEQ
          BNEQ=FTM(SMAX,HEQ)
          HADD=-BNEQ/B0DOT
          IF(HT.LE.HEQ) THEN
            BN=B0-B0DOT*HT
          ELSE IF(HT.LE.(HEQ+HADD)) THEN
            BN=FTM(SL,HT)+B0DOT*(HADD-HT+HEQ)
            IF(BN.LE.0.) BN=ADOME
          ELSE
            BN=FTM(SL,HT)
            IF(BN.LE.0.) BN=ADOME
          ENDIF
          AFUNCT=BN
      GOTO 1000
C
220   CONTINUE
C **** IDOT=10 SX (SUBPOLAR MIX) SLOPE DEPENDENT
          ADOME=.10d0
          SMAX=.07d0
          B0=-5.d0
          HEQ=1250.d0
          B0DOT=B0/HEQ
          BNEQ=FSX(SMAX,HEQ)
          HADD=-BNEQ/B0DOT
          IF(HT.LE.HEQ) THEN
            BN=B0-B0DOT*HT
          ELSE IF(HT.LE.(HEQ+HADD)) THEN
            BN=FSX(SL,HT)+B0DOT*(HADD-HT+HEQ)
            IF(BN.LE.0.) BN=ADOME
          ELSE
            BN=FSX(SL,HT)
            IF(BN.LE.0.) BN=ADOME
          ENDIF
          AFUNCT=BN
      GOTO 1000
C
230   CONTINUE
C **** IDOT=11 PX (POLAR MIX) SLOPE DEPENDENT
          ADOME=.05d0
          SMAX=.04d0
          B0=-2.5d0
          HEQ=600.d0
          B0DOT=B0/HEQ
          BNEQ=FPX(SMAX,HEQ)
          HADD=-BNEQ/B0DOT
          IF(HT.LE.HEQ) THEN
            BN=B0-B0DOT*HT
          ELSE IF(HT.LE.(HEQ+HADD)) THEN
            BN=FPX(SL,HT)+B0DOT*(HADD-HT+HEQ)
            IF(BN.LE.0.) BN=ADOME
          ELSE
            BN=FPX(SL,HT)
            IF(BN.LE.0.) BN=ADOME
          ENDIF
          AFUNCT=BN
      GOTO 1000
C
240   CONTINUE
C **** IDOT=12 PC (POLAR CONTINENTAL) SLOPE DEPENDENT
          ADOME=.05d0
          SMAX=.01d0
          B0=-1.d0
          HEQ=300.d0
          B0DOT=B0/HEQ
          BNEQ=FPC(SMAX,HEQ)
          HADD=-BNEQ/B0DOT
          IF(HT.LE.HEQ) THEN
            BN=B0-B0DOT*HT
          ELSE IF(HT.LE.(HEQ+HADD)) THEN
            BN=FPC(SL,HT)+B0DOT*(HADD-HT+HEQ)
            IF(BN.LE.0.) BN=ADOME
          ELSE
            BN=FPC(SL,HT)
            IF(BN.LE.0.) BN=ADOME
          ENDIF
          AFUNCT=BN
      GOTO 1000
C
300   CONTINUE
C **** IDOT=13 SM EXPONENTIAL
          A1=-7.57655d0
          A2=2.56891d0
          RL1=2.773327d-6
          RL2=2.916233d-7
          SNOLIN=SNO(13)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 900
C
310   CONTINUE
C **** IDOT=14 TM EXPONENTIAL
          A1=-11.6877d0
          A2=3.67407d0
          RL1=1.083791d-6
          RL2=3.641119d-8
          SNOLIN=SNO(14)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 900
C
320   CONTINUE
C **** IDOT=15 SX EXPONENTIAL
          A1=-5.685629d0
          A2=0.857527d0
          RL1=1.322827d-6
          RL2=9.461793d-8
          SNOLIN=SNO(15)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 900
C
330   CONTINUE
C **** IDOT=16 PX EXPONENTIAL
          A1=-2.84904d0
          A2=.837780d0
          RL1=8.534590d-6
          RL2=3.405352d-8
          SNOLIN=SNO(16)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 900
C
340   CONTINUE
C **** IDOT=17 PC EXPONENTIAL
          A1=-1.29363d0
          A2=.299798d0
          RL1=1.502293d-5
          RL2=3.661944d-8
          SNOLIN=SNO(17)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 900
C
345   CONTINUE
C **** IDOT=18 ANTARCTIC VALUUES
C **** FOLLOWING IS CATCHALL FOR DIFFERENT MASS BALANCE RELATIONSHIPS
C **** AND MAY HAVE DIFFERENT VALUES FOR DIFFERENT EXPERIMENTS
C **** DOUBLE OF 17
          A1=-2.d0*1.29363d0
          A2=2.d0*.299798d0
          RL1=1.502293d-5
          RL2=3.661944d-8
C **** CTRL VALUES
C         A1=-.838547d0
C         A2=.887967d0
C         RL1=3.471435d-7
C         RL2=3.918019d-8
C **** REGSST VALUES
C         A1=-1.22491d0
C         A2=1.21169d0
C         RL1=3.736239d-7
C         RL2=8.897420d-8
C **** WARM MIN VALUES
C         A1=-15.8931d0
C         A2=15.3500d0
C         RL1=3.919554d-7
C         RL2=2.362680d-7
          SNOLIN=SNO(18)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 900
C
347   CONTINUE
C **** IDOT=19 GROSSWALD'S DIPPING SNOWLINE SCHEME
      DIST=SQRT(X**2+Y**2)
C FOLLOWING FOR 100M/100KM WITH SNOWLINE AT -1000 M AT POLE
      SNSLOPE=.001d0
C FOLLOWING IS FOR SIBERIAN SNOWLINE DATA, 2E-03
      SNSLOPE=2.d-03
      SNSLOPE=AMASS(8)
      SNZERO=AMASS(7)
C FOLLOWING FOR -1 M AT S.L., +1 M AT 2000 M ELEVATION, OR 2M/2000M
      BALGRAD=.001d0
      BALGRAD=AMASS(9)
      SNLINE=SNZERO+SNSLOPE*DIST
      ACC=(HT-SNLINE)*BALGRAD
      IF(ACC.GT.2.) THEN
        ACC=ACC-2.d0*(ACC-2.d0)
        IF(ACC.LT..05) ACC=.05d0
      ELSEIF(ACC.LT.-2.) THEN
        ACC=-2.d0
      ENDIF
      AFUNCT=ACC
      RETURN
348   CONTINUE
C **** IDOT=20 MIKE'S METEROLOGICAL MODEL

      AMASS(9)=TNSL
      AFUNCT=ACCUM(IDOT,X*.001D0,Y*.001D0,HT,SL,PSURF,TNSL,TS)
      if(.true. .and. afunct.gt.0) then
        afunct=afunct*10
      else
        afunct=afunct*1
      endif
      RETURN
C
350   CONTINUE
C **** IDOT=21 WIND/SLOPE METHOD
C WIND/SLOPE METHODA
C FOR NOW JUST AN AUGMENTED METEROLOGICAL ZONE
c     TNSL=AMASS(9)
C     TNSL=COEF(TIME,-300.,800.,-6.,-15.)
      AMASS(9)=TNSL
c     threshangle=25d0
c     threshangle=cos(threshangle*pi/180d0)
      threshangle=0.80d0
      costheta=slopn(4)/slopn(1)
c     print *,180*costheta/pi,slopn(4),slopn(1)
      if(.true.) then
        AFUNCT=ACCUM(IDOT,X*.001D0,Y*.001D0,HT,SL,PSURF,TNSL,TS)
c    &              +SLOPN(4)*WINDIR(2)+0.001
        acc=afunct+abl
        costheta=0.5*(costheta+1.0)
        acc=acc*costheta
        afunct=acc-abl
      elseif(costheta.gt.threshangle) then
        AFUNCT=ACCUM(IDOT,X*.001D0,Y*.001D0,HT,SL,PSURF,TNSL,TS)+
c    &               SLOPN(4)*WINDIR(2)
     &               0
      else
        AFUNCT=ACCUM(IDOT,X*.001D0,Y*.001D0,HT,SL,PSURF,TNSL,TS)
        afunct=-abl
c       afunct=-0.0003
      endif
      if(.true. .and. afunct.gt.0) then
        afunct=afunct*10
      else
        afunct=afunct*1
      endif
      RETURN
C
360   CONTINUE
C **** IDOT=22 HUYBRECHTS METHOD
c     TNSL=AMASS(9)
C     TNSL=COEF(TIME,-300.,800.,-6.,-15.)
      AMASS(9)=TNSL
      AFUNCT=ACCUMH(X*.001D0,Y*.001D0,HT,TNSL,TS)
      RETURN
370   CONTINUE
C **** IDOT=23 WHITE HOLE METHOD (NO ABLATION ...)
c     TNSL=AMASS(9)
C     TNSL=COEF(TIME,-300.,800.,-6.,-15.)
      AMASS(9)=TNSL
      AFUNCT=WOABLATE(X*.001D0,Y*.001D0,HT,SL,0.D0,TNSL,TS)
      RETURN
380   CONTINUE
C **** IDOT=24 MILANKOVICH METHOD
      AMASS(9)=TNSL
      call milank(TIME,X*.001D0,Y*.001D0,HT,SL,0.D0,TNSL,TS,accm)
      AFUNCT=accm
      RETURN
C
400   CONTINUE
C **** IDOT=25 ISLSCP-DATA-BASED MASS BALANCE
C **** IDOT=26 NCEP2-DATA-BASED MASS BALANCE
C **** IDOT=27 GCM-RESULT-BASED MASS BALANCE
C **** IDOT=28 SEALEVEL-DEPENDENT INTERPOLATION 
C      OF NCEP2-DATA-BASED AND GCM-RESULT-BASED MASS BALANCE
C **** IDOT=29 NCEP2-DATA-BASED MASS BALANCE WITH TOPO MOD
C **** IDOT=30 GCM-RESULT-BASED MASS BALANCE WITH TOPO MOD
C **** IDOT=31 SEALEVEL-DEPENDENT INTERPOLATION 
C      OF NCEP2-DATA-BASED AND GCM-RESULT-BASED MASS BALANCE WITH TOPO MOD
C **** 
      AMASS(9)=TNSL
      AFUNCT=ACCUM(IDOT,X*.001D0,Y*.001D0,HT,SL,PSURF,TNSL,TS)
      RETURN
900   CONTINUE
C
C AGAIN A VARIABLE SECTION FOR DIFFERENT EXPERIMENTS, HERE VARYING
C SPECIFIC PARAMETERS IN THE DIFFERENT MASS BALANCE SCHEMES. THIS CAN BE
C USED TO SMOOTHLY SWITCH FROM ONE ZONE TO ANOTHER.
C FOLLOWING FOR TIME VARYING ACCUM CURVES
C COEF IS FUNCTION COEF(TIME,TINIT,TFINAL,VINIT,VFINAL)
C COEF DOES SINUSOIDAL VARIATION REGARDLES OF TIME
C COEF1 FOR TIME<TINIT = VINIT, TIME>TFINAL = VFINAL, IN BETWEEN IS
C COSINE INTERPOLATION
C FOLLOWING ARE INITIAL TIME AND FINAL TIME (HALF CYCLE)
C     TINIT=0.
C     TFINAL=5000.
C FOLLOWING 4 FOR CONSTANT
      A11=A1
      A22=A2
      RL11=RL1
      RL22=RL2
C FOLLOWING 4 FOR TIME VARYING SINUSOIDAL
C     A11=COEF(TIME,TINIT,TFINAL,A1,A1)
C     A22=COEF(TIME,TINIT,TFINAL,A2,2.*A2)
C     RL11=COEF(TIME,TINIT,TFINAL,RL1,RL1*.01)
C     RL22=COEF(TIME,TINIT,TFINAL,RL2,RL2)
C FOLLOWING 4 FOR COSINE INTERPOLATED STEP
C BEGINNING AND END OF TIME STEP
C     TINIT=1000.
C     TFINAL=22000.
C     A11=COEF1(TIME,TINIT,TFINAL,A1,-14.074)
C     A22=COEF1(TIME,TINIT,TFINAL,A2,4.83112)
C     RL11=COEF1(TIME,TINIT,TFINAL,RL1,1.308677E-06)
C     RL22=COEF1(TIME,TINIT,TFINAL,RL2,1.61273E-7)
C     A11=COEF1(TIME,TINIT,TFINAL,A1,A1*18.95)
C     A22=COEF1(TIME,TINIT,TFINAL,A2,17.287*A2)
C     RL11=COEF1(TIME,TINIT,TFINAL,RL1,RL1*1.129)
C     RL22=COEF1(TIME,TINIT,TFINAL,RL2,6.03*RL2)
C     A11=COEF1(TIME,TINIT,TFINAL,A1,A1*1.461)
C     A22=COEF1(TIME,TINIT,TFINAL,A2,1.365*A2)
C     RL11=COEF1(TIME,TINIT,TFINAL,RL1,RL1*1.076)
C     RL22=COEF1(TIME,TINIT,TFINAL,RL2,2.271*RL2)
C END OF VARIABLE SECTION
C
      HTT=HT-SNOLIN
      IF(HTT.LT.0.) HTT=0.d0
      ARG=-RL11*HTT**2
      IF(ARG.GT.-180.) THEN
        TERM1=A11*EXP(ARG)
      ELSE
        TERM1=0.d0
      ENDIF
      ARG=-RL22*HTT**2
      IF(ARG.GT.-180.) THEN
        TERM2=A22*EXP(ARG)
      ELSE
        TERM2=0.d0
      ENDIF
      AFUNCT=TERM1+TERM2
1000  CONTINUE
      RETURN
      END
      FUNCTION COEF(TIME, TINIT, TFINAL, VINIT, VFINAL)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL TINIT, TFINAL, VINIT, VFINAL
C STARTING AT TIME=TINIT SINUSOIDALLY VARY FROM VINIT TO VFINAL WITH
C PERIOD = TFINAL-TINIT. CONTINUE PERIODICITY.
      IF(TIME.LT.TINIT) THEN
        COEF=DBLE(VINIT)
      ELSE
        ZEROL=0.5d0*DBLE(VINIT+VFINAL)
        AMPL=0.5d0*DBLE(VINIT-VFINAL)
        PHASE=(TIME-DBLE(TINIT))/DBLE(TFINAL-TINIT)*3.1415927d0
        COEF=ZEROL+AMPL*COS(PHASE)
      ENDIF
      END
C
      FUNCTION COEF1(TIME, TINIT, TFINAL, VINIT, VFINAL)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL TINIT, TFINAL, VINIT, VFINAL
C STARTING AT TIME TINIT VARY FROM VINIT TO VFINAL SUNUSOIDALLY. FOR
C TIME < TINIT VALUE=VINIT, FOR TIME > VFINAL VALUE= VFINAL.
      IF(TIME.LE.TINIT) THEN
        COEF1=DBLE(VINIT)
      ELSE IF(TIME.GE.TFINAL) THEN
        COEF1=dble(VFINAL)
      ELSE
        ZEROL=0.5d0*dble(VINIT+VFINAL)
        AMPL=0.5d0*dble(VINIT-VFINAL)
        PHASE=(TIME-dble(TINIT))/dble(TFINAL-TINIT)*3.1415927d0
        COEF1=ZEROL+AMPL*COS(PHASE)
      ENDIF
      END
      FUNCTION COEFL(TIME, TINIT, TFINAL, VINIT, VFINAL)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL TINIT, TFINAL, VINIT, VFINAL
C STARTING AT TIME TINIT VARY FROM VINIT TO VFINAL LINEARLY. FOR
C TIME < TINIT VALUE=VINIT, FOR TIME > VFINAL VALUE= VFINAL.
      IF(TIME.LE.TINIT) THEN
        COEFL=dble(VINIT)
      ELSE IF(TIME.GE.TFINAL) THEN
        COEFL=dble(VFINAL)
      ELSE
        COEFL=dble(VINIT)+(TIME-dble(TINIT))*dble(VFINAL-VINIT)/
     &        dble(TFINAL-TINIT)
      ENDIF
      END
C function based on LOCATE in RECIPES, page 90
      FUNCTION FINDTNSL(TIME)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
      COMMON /TTNSL/ TTIME(NSMAX),TLIST(NSMAX),NTNSL
      DATA IPASS /0/
      SAVE IPASS,TSAVE,TNSLSAVE
      IF(IPASS.EQ.0) THEN
        TSAVE=TIME
        IPASS=1
      ELSEIF(TSAVE.EQ.TIME) THEN
C ..... TIME HASNT CHANGED, USE LAST ONE AND RETURN ...
        FINDTNSL=TNSLSAVE
        RETURN
      ENDIF
      TSAVE=TIME
C if time greater than last ttime, tnsl=last tlist
      IF(TIME.GT.TTIME(NTNSL)) THEN
        TNSLSAVE=TLIST(NTNSL)
        FINDTNSL=TNSLSAVE
        RETURN
      ENDIF
C if time less than first ttime, tnsl=first tlist
      IF(TIME.LT.TTIME(1)) THEN
        TNSLSAVE=TLIST(1)
        FINDTNSL=TNSLSAVE
        RETURN
      ENDIF
c if time = one of the ttimes, then tnsl= that tlist
      DO I=1,NTNSL
        IF(TIME.EQ.TTIME(I)) THEN
          TNSLSAVE=TLIST(I)
          FINDTNSL=TNSLSAVE
          RETURN
        ENDIF
      ENDDO
c interpolate
      N=NTNSL
      JL=0
      JU=N+1
10    CONTINUE
        IF(JU-JL.GT.1) THEN
          JM=(JU+JL)/2
          IF((TTIME(N).GT.TTIME(1)).EQV.(TIME.GT.TTIME(JM))) THEN
            JL=JM
          ELSE
            JU=JM
          ENDIF
          GOTO 10
        ENDIF
      J=JL
      SLOPE=(TLIST(J+1)-TLIST(J))/(TTIME(J+1)-TTIME(J))
      TNSL=TLIST(J)+SLOPE*(TIME-TTIME(J))
      TNSLSAVE=TNSL
      FINDTNSL=TNSLSAVE
      END
c=============================================================
      FUNCTION FINDSEAL(TIME)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      DATA IPASS /0/
      SAVE IPASS,TSAVE,SLSAVE
      IF(IPASS.EQ.0) THEN
        TSAVE=TIME
        IPASS=1
      ELSEIF(TSAVE.EQ.TIME) THEN
C ..... TIME HASNT CHANGED, USE LAST ONE AND RETURN ...
        FINDSEAL=SLSAVE
        RETURN
      ENDIF
      TSAVE=TIME
C if time greater than last ttime, tnsl=last tlist
      IF(TIME.GT.STIME(NSEAL)) THEN
        SLSAVE=SLIST(NSEAL)
        FINDSEAL=SLSAVE
        RETURN
      ENDIF
C if time less than first ttime, tnsl=first tlist
      IF(TIME.LT.STIME(1)) THEN
        SLSAVE=SLIST(1)
        FINDSEAL=SLSAVE
        RETURN
      ENDIF
c if time = one of the ttimes, then tnsl= that tlist
      DO I=1,NSEAL
        IF(TIME.EQ.STIME(I)) THEN
          SLSAVE=SLIST(I)
          FINDSEAL=SLSAVE
          RETURN
        ENDIF
      ENDDO
c interpolate
      N=NSEAL
      JL=0
      JU=N+1
10    CONTINUE
        IF(JU-JL.GT.1) THEN
          JM=(JU+JL)/2
          IF((STIME(N).GT.STIME(1)).EQV.(TIME.GT.STIME(JM))) THEN
            JL=JM
          ELSE
            JU=JM
          ENDIF
          GOTO 10
        ENDIF
      J=JL
      SLOPE=(SLIST(J+1)-SLIST(J))/(STIME(J+1)-STIME(J))
      SL=SLIST(J)+SLOPE*(TIME-STIME(J))
      SLSAVE=SL
      FINDSEAL=SLSAVE
      END
