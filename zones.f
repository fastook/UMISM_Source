C      FUNCTION AFUNCT(TIME, IDOT, AMASS, HT, BD, SLOPN, X, Y)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /SNOW/ SNOLIN,SNO(20)
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      COMMON /TTNSL/ TTIME(20000),TLIST(20000),NTNSL
      COMMON /EXPER/ IEXPER
      DIMENSION AMASS(11),SLOPN(4),elev(1000),rmass(1000)
      FPC(S,E)=(-7.044897460E+02)*S**2+( 3.286502080E+01)*S+
     1  (-8.363957700E-08)*E**2+( 2.803513780E-04)*E+
     1  (-8.280271290E-02)
      FPX(S,E)=(-1.011654790E+03)*S**2+( 5.757759090E+01)*S+
     1  (-1.945462600E-07)*E**2+( 6.485106420E-04)*E+
     1  (-1.714566350E-01)
      FSX(S,E)=( 1.613718720E+02)*S**2+(-3.064633370E+00)*S+
     1  (-2.274655910E-07)*E**2+( 9.221490470E-04)*E+
     1  (-3.658012150E-01)
      FTM(S,E)=(-2.475505130E+03)*S**2+( 7.388871770E+01)*S+
     1  (-7.129281700E-07)*E**2+( 3.180227710E-03)*E+
     1  (-6.778311130E-01)
      FSM(S,E)=(-1.272729690E+04)*S**2+( 1.028782810E+02)*S+
     1  ( 1.790857030E-07)*E**2+(-1.209567770E-03)*E+
     1  ( 2.244830130E+00)
      call grstrt(800,800)
      call window(0.,5000.,-5.,5.)
      icolor=1
23    print *,'input zone number'
      read(*,*,end=999) idot
      if(idot.eq.0) goto 999
      call linclr(icolor)
      do i=1,9
        amass(i)=0.
      enddo
      ic=0
      do 1001 ht=0.,5000.,10.
      ic=ic+1
      elev(ic)=ht
C
C BRANCH FOR VARIOUS MASS BALANCE RELATIONSHIPS
      GOTO(10,20,30,40,50,60,70,200,210,220,230,240,
     1     300,310,320,330,340,345,347,348,350),IDOT
C
10    CONTINUE
          AMARG=AMASS(1)
          APEAK=AMASS(2)
          ADOME=AMASS(3)
          HE=AMASS(4)
          HP=AMASS(5)
          HM=AMASS(6)
      GOTO 100
C
20    CONTINUE
C **** IDOT=2 POLAR - A
          AMARG=0.0
          APEAK=.5
          ADOME=.05
          HE=0.
          HP=1250.
          HM=2500.
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
30    CONTINUE
C **** IDOT=3 POLAR - B
          AMARG=0.0
          APEAK=1.
          ADOME=.10
          HE=0.
          HP=1250.
          HM=2500.
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
40    CONTINUE
C **** IDOT=4 POLAR - C
          AMARG=-1.0
          APEAK=1.5
          ADOME=.1
          HE=100.
          HP=1250.
          HM=2500.
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
50    CONTINUE
C **** IDOT=5 MARITIME - A
          AMARG=-2.0
          APEAK=2.
          ADOME=.15
          HE=1500.
          HP=1750.
          HM=2000.
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
60    CONTINUE
C **** IDOT=6 MARITIME - B
          AMARG=-2.0
          APEAK=2.
          ADOME=.15
          HE=1000.
          HP=1500.
          HM=2000.
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
70    CONTINUE
C **** IDOT=7 MARITIME - C
          AMARG=-3.0
          APEAK=2.
          ADOME=.15
          HE=1000.
          HP=1500.
          HM=2000.
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
            H1=0.
            A2=0.
            H2=HE
          ENDIF
          IF(HTT.GT.HE .AND. HTT.LE.HP) THEN
            A1=0.
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
      IF(BD.LE.-5000.) AFUNCT=-5.
      goto 1000
C
200   CONTINUE
C **** IDOT=8 SM (SUBPOLAR MARITIME) SLOPE DEPENDENT
          ADOME=.10
          SMAX=.01
          B0=-5.
          HEQ=700.
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
          ADOME=.10
          SMAX=.03
          B0=-8.
          HEQ=1125.
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
          ADOME=.10
          SMAX=.07
          B0=-5.
          HEQ=1250.
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
          ADOME=.05
          SMAX=.04
          B0=-2.5
          HEQ=600.
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
          ADOME=.05
          SMAX=.01
          B0=-1.
          HEQ=300.
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
          A1=-7.57655
          A2=2.56891
          RL1=2.773327E-6
          RL2=2.916233E-7
          SNOLIN=SNO(13)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 400
C
310   CONTINUE
C **** IDOT=14 TM EXPONENTIAL
          A1=-11.6877
          A2=3.67407
          RL1=1.083791E-6
          RL2=3.641119E-8
          SNOLIN=SNO(14)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 400
C
320   CONTINUE
C **** IDOT=15 SX EXPONENTIAL
          A1=-5.685629
          A2=0.857527
          RL1=1.322827E-6
          RL2=9.461793E-8
          SNOLIN=SNO(15)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 400
C
330   CONTINUE
C **** IDOT=16 PX EXPONENTIAL
          A1=-2.84904
          A2=.837780
          RL1=8.534590E-6
          RL2=3.405352E-8
          SNOLIN=SNO(16)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 400
C
340   CONTINUE
C **** IDOT=17 PC EXPONENTIAL
          A1=-1.29363
          A2=.299798
          RL1=1.502293E-5
          RL2=3.661944E-8
          SNOLIN=SNO(17)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 400
C
345   CONTINUE
C **** IDOT=18 ANTARCTIC VALUUES
C **** FOLLOWING IS CATCHALL FOR DIFFERENT MASS BALANCE RELATIONSHIPS
C **** AND MAY HAVE DIFFERENT VALUES FOR DIFFERENT EXPERIMENTS
C **** DOUBLE OF 17
          A1=-2.*1.29363
          A2=2.*.299798
          RL1=1.502293E-5
          RL2=3.661944E-8
C **** CTRL VALUES
C         A1=-.838547
C         A2=.887967
C         RL1=3.471435E-7
C         RL2=3.918019E-8
C **** REGSST VALUES
C         A1=-1.22491
C         A2=1.21169
C         RL1=3.736239E-7
C         RL2=8.897420E-8
C **** WARM MIN VALUES
C         A1=-15.8931
C         A2=15.3500
C         RL1=3.919554E-7
C         RL2=2.362680E-7
          SNOLIN=SNO(18)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 400
C
347   CONTINUE
C **** IDOT=19 GROSSWALD'S DIPPING SNOWLINE SCHEME
      DIST=SQRT(X**2+Y**2)
C FOLLOWING FOR 100M/100KM WITH SNOWLINE AT -1000 M AT POLE
      SNSLOPE=.001
C FOLLOWING IS FOR SIBERIAN SNOWLINE DATA, 2E-03
      SNSLOPE=2.E-03
      SNSLOPE=AMASS(8)
      SNZERO=AMASS(7)
C FOLLOWING FOR -1 M AT S.L., +1 M AT 2000 M ELEVATION, OR 2M/2000M
      BALGRAD=.001
      BALGRAD=AMASS(9)
      SNLINE=SNZERO+SNSLOPE*DIST
      ACC=(HT-SNLINE)*BALGRAD
      IF(ACC.GT.2.) THEN
        ACC=ACC-2.*(ACC-2.)
        IF(ACC.LT..05) ACC=.05
      ELSEIF(ACC.LT.-2.) THEN
        ACC=-2.
      ENDIF
      AFUNCT=ACC
      goto 1000
348   CONTINUE
C **** IDOT=20 MIKE'S METEROLOGICAL MODEL

      AMASS(9)=TNSL
      AFUNCT=ACCUM(X*.001D0,Y*.001D0,HT,SL,0.D0,TNSL,TS)
      goto 1000
C
350   CONTINUE
C **** IDOT=21 WIND/SLOPE METHOD
C WIND/SLOPE METHODA
C FOR NOW JUST AN AUGMENTED METEROLOGICAL ZONE
      TNSL=AMASS(9)
C     TNSL=COEF(TIME,-300.,800.,-6.,-15.)
      AMASS(9)=TNSL
      AFUNCT=ACCUM(X*.001D0,Y*.001D0,HT,SL,0.D0,TNSL,TS)+
     &             SLOPN(4)*WINDIR(2)
      goto 1000
C
400   CONTINUE
      a11=a1
      rl11=rl1
      a22=a2
      rl22=rl2
C
      HTT=HT-SNOLIN
      IF(HTT.LT.0.) HTT=0.
      ARG=-RL11*HTT**2
      IF(ARG.GT.-180.) THEN
        TERM1=A11*EXP(ARG)
      ELSE
        TERM1=0.
      ENDIF
      ARG=-RL22*HTT**2
      IF(ARG.GT.-180.) THEN
        TERM2=A22*EXP(ARG)
      ELSE
        TERM2=0.
      ENDIF
      AFUNCT=TERM1+TERM2
1000  CONTINUE
      rmass(ic)=afunct
      write(23,*) elev(ic),rmass(ic),idot
1001  continue
      call move(real(elev(1)),real(rmass(1)))
      n=ic
      do ic=1,n
      call draw(real(elev(ic)),real(rmass(ic)))
      enddo
      icolor=icolor+1
      goto 23
999   continue
      END
      FUNCTION ACCUM(X,Y,ELEV1,SLOPE1,SHAPE1,TNSL,TS)                      
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION QI(12),QS(12),TTT(12)                                   
      DATA QI /960.,1036.,1200.,825.,330.,90.,150.,600.,1200.,1020.,    
     &         930.,850./                                               
      DATA QS /-0.667,4.6,11.667,9.167,3.667,1.,1.667,6.667,12.,6.333,  
     &         0.333,-3.333/                                            
C     DATA AAA /-9.14/, BBB /-.68/, CCC /34.461/                        
C     DATA WWW /13.05/, XXX /.664/, ZZZ /2.608/                         
       DATA AAA / -9.62376690     /                                     
       DATA BBB /-0.546917617     /                                     
       DATA CCC /  24.9793854     /                                     
       DATA WWW /  19.1390686     /                                     
       DATA XXX / 0.922791243     /                                     
       DATA ZZZ /-0.738900483     /                                     
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      AAA=ACOM                                                          
C     XXX=0.                                                            
C     WRITE(19,*) X,Y,ELEV1,SLOPE1                                      
C CALCULATE LATITUDE                                                    
      CALL SETRIG                                                       
      CALL RECPOL(X,Y,RLAT,RLONG)                                       
C     PRINT *,'LATITUDE=',RLAT                                          
C ELEVATION (KM)                                                        
      ELEV=ELEV1/1000.                                                  
C SLOPE (M/KM)                                                          
      SLOPE=SLOPE1*1000.                                                
C SHAPE (M/KM/KM)                                                       
      SHAPE=SHAPE1*1000.*1000.                                          
      SHAPE=0.                                                          
C CALCULATE SURFACE MEAN ANNUAL AIR TEMP                                
      TS=AAA*ELEV+BBB*RLAT+CCC+TNSL+TNSLBASE                                
C     PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TS                        
C CALCULATE MEAN ANNUAL TEMP OF FREE ATMOSPHERE-ISOTHEMAL LAYER         
      TF=0.67*(TS+273.0)+88.9                                           
C     PRINT *,'TEMP FREE AT-ISOTHRMAL LAYER=',TF                        
C CALCULATE SATURATION VAPOR PRESSURE                                   
      TERM1=-9.09718*(273.16/TF-1.0)                                    
      TERM2=-3.56654*LOG10(273.16/TF)                                   
      TERM3=0.876793*(1.0-TF/273.16)+0.785835                           
      EXPON=TERM1+TERM2+TERM3                                           
      ES=10.**EXPON                                                     
C     PRINT *,'SATURATION VAPOR PRESSURE=',ES                           
C CALCULATE ACCUMULATION RATE (M/YR)                                    
      TERM1=WWW*ES                                                      
      TERM2=XXX*SLOPE                                                   
      TERM3=ZZZ                                                         
      TERM4=-15.276*SHAPE                                               
      ACC=TERM1+TERM2+TERM3+TERM4                                       
C     WRITE(19,113) TERM1,TERM2,TERM3,ACC                               
113   FORMAT(4G13.6)                                                    
C     PRINT *,'ACCUMULATION=',ACC                                       
C CALCULATE ABLATION                                                    
      QY=0.                                                             
      DO 10 I=1,12                                                      
        QY=QY+QI(I)-QS(I)*RLAT                                          
10    CONTINUE                                                          
      QY=QY/12.                                                         
      PDD=0.                                                            
      DO 20 I=1,12                                                      
C       TTT(I)=TS+0.021*((QI(I)+QS(I)*RLAT)-QY)+8.954                   
        TTT(I)=TS+0.021*((QI(I)-QS(I)*RLAT)-QY)+0.0                     
C       WRITE(19,*) I,TTT(I),QI(I)-QS(I)*RLAT,QY                        
        IF(TTT(I).GT.0.0) PDD=PDD+30.*TTT(I)                            
20    CONTINUE                                                          
      ABL=.6*PDD                                                        
C     PRINT *,'ABLATION=',ABL                                           
C CALCULATE NET ACCUMULATION                                            
      ACCNET=ACC-ABL                                                    
C     ACCNET=ACC                                                        
C     PRINT *,'NET ACCUMULATION/ABLATION=',ACCNET                       
      ACCUM=ACCNET*.01                                                  
C     IF(PDD.GT.0.) THEN                                                
C       WRITE(19,111) ELEV,SLOPE,RLAT,TS,TF-273.,ES,PDD,ACC,ABL,ACCNET  
C     ENDIF                                                             
111   FORMAT(6F7.2,F6.0,3F7.1)                                          
C     WRITE(19,112) (TTT(I),I=1,12)                                     
112   FORMAT(12F6.1)                                                    
      END                                                               

      SUBROUTINE SETRIG
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      PI=4.D0*ATAN(1.D0)
      RADIUS=2.D4/PI
      RADIUS=RADIUS*0.53
      CIRCUM=2.D0*PI*RADIUS
      RKMPDEG=CIRCUM/360.D0
      RADPDEG=PI/180.D0
      DEGPRAD=180.D0/PI
      END
      SUBROUTINE POLREC(RLAT,RLONG,X,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      y= 1000*rkmpdeg*rlat
      x= 1000*rkmpdeg*cos(rlat*radpdeg)*(rlong+127.5)
      END
      SUBROUTINE RECPOL(X,Y,RLAT,RLONG)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      rlat=y*0.001/rkmpdeg
      rlong=-127.5+x*0.001/rkmpdeg/cos(rlat*radpdeg)
      END
