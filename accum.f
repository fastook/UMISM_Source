      FUNCTION ACCUM(IDOT,X,Y,ELEV1,SLOPE1,PSURF1,TNSL,TS)                      
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM,NSMAX=MAXTIME)
      COMMON /TCONST/ TSORIG(NMAX),TFLAG                         
      LOGICAL TFLAG
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      DIMENSION QI(12),QS(12),TTT(12)                                   
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      COMMON /LOCAL/ SIGMA,TMAA,TMS,PI2
      DATA QI /960.d0,1036.d0,1200.d0,825.,330.d0,90.d0,150.d0,
     &         600.d0,1200.d0,1020.d0,    
     &         930.d0,850.d0/                                               
      DATA QS /-0.667d0,4.6d0,11.667d0,9.167d0,3.667d0,
     *         1.d0,1.667d0,6.667d0,12.d0,6.333d0,  
     &         0.333d0,-3.333d0/                                            
C     DATA AAA /-9.14/, BBB /-.68/, CCC /34.461/                        
C     DATA WWW /13.05/, XXX /.664/, ZZZ /2.608/                         
       DATA AAA / -9.62376690d0     /                                     
       DATA BBB /-0.546917617d0     /                                     
       DATA CCC /  24.9793854d0     /                                     
       DATA WWW /  19.1390686d0     /                                     
       DATA XXX / 0.922791243d0     /                                     
       DATA ZZZ /-0.738900483d0     /  
      SHAPE1=0.d0                                   
      AAA=ACOM 
cC THIS IS SPECIAL FOR MISHA'S BAIKAL STUFF
c      BBB=-1
C     XXX=0.                                                            
C     WRITE(19,*) X,Y,ELEV1,SLOPE1                                      
C ... CALCULATE LATITUDE                                                    
C     CALL SETRIG 
      if(IDOT.EQ.25) then
c ..... with the lookup mass balance, there is no "pole"                                                      
        CALL RECPOL(X,Y,RLAT,RLONG)  
c        dist=sqrt((X-XPOLE)**2+(Y-YPOLE)**2)/5d3
c       tnslredu=tnsl*exp(-0.5d0*dist**2)
        tnslredu=tnsl
c        print *,real(X-XPOLE),real(Y-YPOLE),real(dist),tnslredu
        call lookup25(AAA,rlat,rlong,elev1,
     &              tnslredu,TS,acc,abl)
        ACCUM=acc
        if(accum.ne.-999.) return
      endif
      if(IDOT.EQ.26) then
c ..... with the lookup mass balance, there is no "pole"                                                      
        CALL RECPOL(X,Y,RLAT,RLONG)  
c        dist=sqrt((X-XPOLE)**2+(Y-YPOLE)**2)/2d3
c        tnslredu=tnsl*exp(-0.5d0*dist**2)
        tnslredu=tnsl
c        print *,real(X-XPOLE),real(Y-YPOLE),real(dist),tnslredu
        call lookup26(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS,acc,abl)
        ACCUM=acc
        if(accum.ne.-999.) return
      endif
      if(IDOT.EQ.27) then
c ..... with the lookup mass balance, there is no "pole"                                                      
        CALL RECPOL(X,Y,RLAT,RLONG)  
c        dist=sqrt((X-XPOLE)**2+(Y-YPOLE)**2)/2d3
c        tnslredu=tnsl*exp(-0.5d0*dist**2)
        tnslredu=tnsl
c        print *,real(X-XPOLE),real(Y-YPOLE),real(dist),tnslredu
        call lookup27(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS,acc,abl)
        ACCUM=acc
        if(accum.ne.-999.) return
      endif
      if(IDOT.EQ.28) then
c ..... with the lookup mass balance, there is no "pole"
        CALL RECPOL(X,Y,RLAT,RLONG)  
c        dist=sqrt((X-XPOLE)**2+(Y-YPOLE)**2)/2d3
c        tnslredu=tnsl*exp(-0.5d0*dist**2)
        tnslredu=tnsl
c        print *,real(X-XPOLE),real(Y-YPOLE),real(dist),tnslredu
        call lookup26(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS26,acc26,abl26)
        call lookup27(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS27,acc27,abl27)
        slfact=max(0.d0,min(1.d0,-SEALEV/50.d0))
c        print *,slfact,SEALEV
        TS=(1.d0-slfact)*TS26+slfact*TS27
        abl=(1.d0-slfact)*abl26+slfact*abl27
        ACCUM=(1.d0-slfact)*acc26+slfact*acc27
        if(accum.ne.-999.) return
      endif
      if(IDOT.EQ.29) then
c ..... with the lookup mass balance, there is no "pole"                                                      
        CALL RECPOL(X,Y,RLAT,RLONG)  
c        dist=sqrt((X-XPOLE)**2+(Y-YPOLE)**2)/2d3
c        tnslredu=tnsl*exp(-0.5d0*dist**2)
        tnslredu=tnsl
c        print *,real(X-XPOLE),real(Y-YPOLE),real(dist),tnslredu
        call lookup29(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS,acc,abl)
        ACCUM=acc
        if(accum.ne.-999.) return
      endif
      if(IDOT.EQ.30) then
c ..... with the lookup mass balance, there is no "pole"                                                      
        CALL RECPOL(X,Y,RLAT,RLONG)  
c        dist=sqrt((X-XPOLE)**2+(Y-YPOLE)**2)/2d3
c        tnslredu=tnsl*exp(-0.5d0*dist**2)
        tnslredu=tnsl
c        print *,real(X-XPOLE),real(Y-YPOLE),real(dist),tnslredu
        call lookup30(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS,acc,abl)
        ACCUM=acc
        if(accum.ne.-999.) return
      endif
      if(IDOT.EQ.31) then
c ..... with the lookup mass balance, there is no "pole"
        CALL RECPOL(X,Y,RLAT,RLONG)  
c        dist=sqrt((X-XPOLE)**2+(Y-YPOLE)**2)/2d3
c        tnslredu=tnsl*exp(-0.5d0*dist**2)
        tnslredu=tnsl
c        print *,real(X-XPOLE),real(Y-YPOLE),real(dist),tnslredu
        call lookup29(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS29,acc29,abl29)
        call lookup30(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS30,acc30,abl30)
        slfact=max(0.d0,min(1.d0,-SEALEV/50.d0))
c        print *,slfact,SEALEV
        TS=(1.d0-slfact)*TS29+slfact*TS30
        abl=(1.d0-slfact)*abl29+slfact*abl30
        ACCUM=(1.d0-slfact)*acc29+slfact*acc30
        if(accum.ne.-999.) return
      endif
c ... but the fastook-prentice scheme can have a "pole"
      CALL RECPOL(X-XPOLE,Y-YPOLE,RLAT,RLONG)  
c ..........................................................
C      PRINT *,'X,Y,ELEV,SLOPE',X,Y,ELEV1,SLOPE1                                     
C      PRINT *,'LATITUDE=',RLAT                                          
C ... ELEVATION (KM)                                                        
      ELEV=ELEV1/1000.d0
      PSURF=PSURF1/1000.d0
C ... SLOPE (M/KM)                                                          
      SLOPE=SLOPE1*1000.d0
C ... SHAPE (M/KM/KM)                                                       
      SHAPE=SHAPE1*1000.d0*1000.d0
      SHAPE=0.d0
c ***********************************************************                               
c ***********************************************************                               
c *******EISMINT EXPERIMENT FOR GREENLAND *******************
C HUYBRECT'S METHOD FROM A LOOKUP TABLE
C ... PUT .TRUE. FOR EISMINT EXPERIMENT ...
      IF(.FALSE.) THEN
        IF(.TRUE.) THEN
C ....... GREENLAND
          CALL HUYB(ELEV1,RLAT,PDD,TS,TSUMMER,TNSL)
        ELSEIF(.FALSE.) THEN
C ....... ANTARCTICA
          CALL HUYBANT(ELEV1,RLAT,PDD,TS,TSUMMER,TNSL)
          ACC=1.5D2*2.D0**(TS/10.D0)
        ELSEIF(.FALSE.) THEN
C ....... SLIDING EXPERIMENT
          CALL PAYNE(X,Y,TS,ACC)
          PDD=0.0d0
        ENDIF
C ... FOR ICE
        ABL=0.8d0*PDD
C ... FOR SNOW
C      ABL=0.3d0*PDD
C ... WHAT I HAVE BEEN USING
C      ABL=0.6d0*PDD
        ACC=0.0d0
        ACCNET=ACC-ABL                                                    
        ACCUM=ACCNET*.01d0
C        PRINT *,ACC,-ABL,ACC-ABL                                                
        RETURN
      ENDIF
C ***********************************************************                               
C ***********************************************************                               
c ***********************************************************                               
C ... CALCULATE SURFACE MEAN ANNUAL AIR TEMP      

c ... EXPERIMENTAL !!! >>>>>                          
c        tnslredu=0.5*(1.+RLAT/90.)
c        TS=AAA*ELEV+BBB*RLAT+CCC+(TNSL*tnslredu)+TNSLBASE            
c ... EXPERIMENTAL !!! >>>>>                          

      TS=AAA*ELEV+BBB*RLAT+CCC+TNSL+TNSLBASE            
c**************MARS***************************************
      IF(.true.) THEN
        IF(TFLAG) TS=TNSL+TNSLBASE-100
        IF(TFLAG) TS=AAA*ELEV+TNSL+TNSLBASE-100
        ACCUM=accmars(ts,elev*1000,acc,abl)
c       ACCUM=accmars1(ts,elev*1000,acc,abl)
        return
      ENDIF
c**************MARS***************************************
c**************EISMINT***************************************
      IF(.false.) THEN
        Z=MAX(ELEV1,20.d0*(RLAT-65.d0))
        TS=49.13d0-0.007992d0*Z-0.7576d0*RLAT+TNSL
      ENDIF
c**************EISMINT***************************************
C     PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TS                        
C ... CALCULATE MEAN ANNUAL TEMP OF FREE ATMOSPHERE-ISOTHEMAL LAYER         
      TF=0.67d0*(TS+273.0d0)+88.9d0
C     PRINT *,'TEMP FREE AT-ISOTHRMAL LAYER=',TF                        
C ... CALCULATE SATURATION VAPOR PRESSURE                                   
      TERM1=-9.09718d0*(273.16d0/TF-1.0d0)                                    
      TERM2=-3.56654d0*LOG10(273.16d0/TF)                                   
      TERM3=0.876793d0*(1.0d0-TF/273.16d0)+0.785835d0
      EXPON=TERM1+TERM2+TERM3                                           
      ES=10.d0**EXPON                                                     
C     PRINT *,'SATURATION VAPOR PRESSURE=',ES                           
C ... CALCULATE ACCUMULATION RATE (M/YR)                                    
      TERM1=WWW*ES                                                      
      TERM2=XXX*SLOPE 

c **** <<< EXPERIMENTAL >>> **** turn off slope term ...
      TERM2=0
c **** <<< EXPERIMENTAL >>> ****

      TERM3=ZZZ                                                         
      TERM4=-15.276d0*SHAPE                                               
      ACC=max(0.d0,TERM1+TERM2+TERM3+TERM4)
C     WRITE(19,113) TERM1,TERM2,TERM3,ACC 
113   FORMAT(4G13.6)                                                    
C     PRINT *,'ACCUMULATION=',ACC                                       
C ... CALCULATE ABLATION                                                    
      QY=0.d0
      DO I=1,12                                                      
        QY=QY+QI(I)-QS(I)*RLAT                                          
      ENDDO                                                          
      QY=QY/12.d0
      PDD=0.d0
      DO I=1,12                                                      
C       TTT(I)=TS+0.021d0*((QI(I)+QS(I)*RLAT)-QY)+8.954d0
        TTT(I)=TS+0.021d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C       TTT(I)=TS+0.018d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C       WRITE(19,*) I,TTT(I),QI(I)-QS(I)*RLAT,QY                        
        IF(TTT(I).GT.0.0) PDD=PDD+30.d0*TTT(I)                            
      ENDDO
c      print *,ts,pdd
      ABL=.6d0*PDD                                                        
C ... CALCULATE NET ACCUMULATION                                            
      ACCNET=ACC-ABL                                                    
C     ACCNET=ACC                                                        
C     PRINT *,'NET ACCUMULATION/ABLATION=',ACCNET                       
      ACCUM=ACCNET*.01d0                                                  
C     IF(PDD.GT.0.) THEN                                                
C       WRITE(19,111) ELEV,SLOPE,RLAT,TS,TF-273.,ES,PDD,ACC,ABL,ACCNET  
C     ENDIF                                                             
111   FORMAT(6F7.2,F6.0,3F7.1)                                          
C     WRITE(19,112) (TTT(I),I=1,12)                                     
112   FORMAT(12F6.1)        
      END                                                               
c=====================================================================
      FUNCTION ACCUMH(X,Y,ELEV1,TNSL,TMA)                      
      IMPLICIT REAL*8(A-H,O-Z)
      external func
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      COMMON /LOCAL/ SIGMA,TMAA,TMS,PI2
      DATA A1 /34.46D0/, A2 /-0.00914D0/, A3 /-0.68775D0/
      DATA B1 /16.81D0/, B2 /-0.00692D0/, B3 /-0.27937D0/
      DATA C1 /0.78D0/, C2 /2.525D-2/, C3 /2.225D-4/
      sigma=5.d0
      PI2=8.d0*ATAN(1.d0)
      ACOM=A2
c      PRINT *,'X,Y,ELEV=',X,Y,ELEV1                                      
C ... CALCULATE LATITUDE                                                    
C     CALL SETRIG                                                       
      CALL RECPOL(X,Y,RLAT,RLONG)                                       
c      PRINT *,'LATITUDE=',RLAT                                          
C ... CALCULATE SURFACE MEAN ANNUAL AIR TEMP                                
      TMA=A1+A2*ELEV1+A3*RLAT+TNSL+TNSLBASE
      TMAA=TMA                               
c      PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TMA                        
C ... CALCULATE SUMMER ANNUAL TEMP
      tms=B1+B2*ELEV1+B3*RLAT+TNSL+TNSLBASE                                
c      PRINT *,'TEMP SUMMER=',tms  
C ... CALCULATE ACCUMULATION RATE (M/YR)                                    
      ACC=C1+C2*TMA+C3*TMA**2                                       
c      PRINT *,'ACCUMULATION=',ACC                                       
C ... CALCULATE ABLATION 
      CALL QUAD2D(0.D0,365.D0,PDD)
      PDD=PDD/SIGMA/SQRT(PI2)
c      PRINT *,'PDD=',PDD
      ABL=.008d0*PDD                                                        
c      PRINT *,'ABLATION=',ABL                                           
C ... CALCULATE NET ACCUMULATION                                            
      ACCNET=ACC-ABL                                                    
c      PRINT *,'NET ACCUMULATION/ABLATION=',ACCNET                       
      ACCUMH=ACCNET                                                  
      END                                                               
C
      SUBROUTINE QUAD2D(X1,X2,SS)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL H
      CALL QGAUSX(H,X1,X2,SS)
c      PRINT *,'SS IN QUAD2D',SS
      RETURN
      END
C     
      FUNCTION F(YY)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FUNC
      COMMON /XY/ X,Y
      Y=YY
      F=FUNC(X,Y)
      RETURN
      END
C     
      FUNCTION H(XX)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL Y1,Y2,F
      COMMON /XY/ X,Y
      X=XX
      CALL QGAUSY(F,Y1(X),Y2(X),SS)
      H=SS
      RETURN
      END
C
      FUNCTION FUNC(X,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL TD
      COMMON /LOCAL/ SIGMA,TMA,TMS,PI2
      temp=Y*EXP(-0.5d0*(Y-TD(X))**2/SIGMA)
      func=temp
      END
C
      FUNCTION Y1(X)
      IMPLICIT REAL*8(A-H,O-Z)
      Y1=0.d0*X
      END
C
      FUNCTION Y2(X)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL TD
      COMMON /LOCAL/ SIGMA,TMA,TMS,PI2
      Y2=TD(X)+2.5d0*SIGMA
c      Y2=TD(X)
      END
C
      FUNCTION TD(T)
      IMPLICIT REAL*8(A-H,O-Z)
      INTRINSIC COS
      COMMON /LOCAL/ SIGMA,TMA,TMS,PI2
      TD=TMA+(TMA-TMS)*COS(PI2*T/365.d0)
      END
C
      SUBROUTINE QGAUSX(F,A,B,SS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(ITMAX=20)
      DIMENSION X(ITMAX),W(ITMAX)
      EXTERNAL F
      DATA TOL /1D-2/
C      SUMOLD=1d30
      N=ITMAX
C      DO N=ITMAX,ITMAX
        CALL GAULEG(A,B,X,W,N)
        SUM=0.d0
        DO I=1,N
          SUM=SUM+W(I)*F(X(I))
        ENDDO
C        IF(ABS(SUM-SUMOLD).LT.ABS(SUM*TOL) .or. 
C     &          sum-sumold.eq.0.) THEN
C          SS=SUM
C           PRINT *,'qgausx',N,ss
C          RETURN
C        ENDIF
C        PRINT *,'qgausx',I,SUM,SUMOLD
C        SUMOLD=SUM
C      ENDDO
      ss=sum
C      PRINT *,'QGAUSX DIDNT CONVERGE IN ',ITMAX
      END
C
      SUBROUTINE QGAUSY(F,A,B,SS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(ITMAX=8)
      DIMENSION X(ITMAX),W(ITMAX)
      EXTERNAL F
      DATA TOL /1D-2/
C      SUMOLD=1d30
      N=ITMAX
C      DO N=ITMAX,ITMAX
        CALL GAULEG(A,B,X,W,N)
        SUM=0.d0
        DO I=1,N
          SUM=SUM+W(I)*F(X(I))
        ENDDO
C        IF(ABS(SUM-SUMOLD).LT.ABS(SUM*TOL) .or. 
C     &          sum-sumold.eq.0.) THEN
C          SS=SUM
c          PRINT *,'qgausy',N,ss
C          RETURN
C        ENDIF
C        SUMOLD=SUM
C      ENDDO
      ss=sum
C      PRINT *,'QGAUSY DIDNT CONVERGE IN ',ITMAX
      END
C
      SUBROUTINE GAULEG(X1,X2,X,W,N)
C
C ... GIVEN THE LOWER AND UPPER LIMITS OF INTEGRATION X1 AND X2, AND
C ... GIVEN THE NUMBER OF NODAL POINTS DESIRED N, THIS ROUTINE RETURNS
C ... ARRAYS X AND W OF LENGTH N, CONTAINING THE NODES AND WEIGHTS OF
C ... THE GUASS-LEGENDRE N-POINT QUADRATURE FORMULA
C
C ... HIGH PRECISION IS A GOOD IDEA FOR THIS ROUTINE
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X1,X2,X(N),W(N)
C ... INCREASE EPS IF YOU DON'T HAVE THIS PRECISION
      PARAMETER (EPS=3.D-14)
C ... THE ROOTS ARE SYMMETRIC IN THE INTERVAL, SO WE ONLY HAVE TO
C ... FIND HALF OF THEM
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
C ... LOOP OVER THE DESIRED ROOTS
      DO 12 I=1,M
        Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))
C ..... STARTING WITH THE ABOVE APPROXIMATION TO THE ITH ROOT,
C ..... WE ENTER THE MAIN LOOP OF REFINEMENT BY THE NEWTON'S METHOD
1       CONTINUE
          P1=1.D0
          P2=0.D0
C ....... LOOP UP THE RECURRENCE RELATION TO GET THE LEGENDRE
C ....... POLYNOMIAL EVALUATED AT Z
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11        CONTINUE
C ....... P1 IS NOW THE DESIRED LEGENDRE POLYNOMIAL. WE NEXT COMPUTE PP,
C ....... ITS DIRIVATIVE, BY A STANDARD RELATION INVOLVING ALSO P2,
C ....... THE POLYNIMIAL OF ONE LOWER ORDER
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
C ....... NEWTON'S METHOD
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
C ..... SCALE THE ROOT TO ITS DESIRED INTERVAL
        X(I)=XM-XL*Z
C ..... AND PUT IN ITS SYMMETRIC COUNTERPART
        X(N+1-I)=XM+XL*Z
C ..... COMPUTE THE WEIGHT
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
C ..... AND ITS SYMMETRIC COUNTERPART
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END
c================================================
      subroutine huyb(elev1,rlat,pdd,TS,TSUMMER,TNSL)
      IMPLICIT REAL*8(A-H,O-Z)
c following covers greenland
      PARAMETER(JMAX=8,RSTEP=5.D0,RMIN=55.D0,
     &          RMAX=RMIN+(JMAX-1)*RSTEP)
c following for full n. hem
c      PARAMETER(JMAX=19,RSTEP=5.D0,RMIN=0.D0,
c     &          RMAX=RMIN+(JMAX-1)*RSTEP)
c
      PARAMETER(IMAX=26,ESTEP=200.D0,EMIN=0.D0,
     &          EMAX=EMIN+(IMAX-1)*ESTEP)
c high resolution ...
c      PARAMETER(JMAX=31,RSTEP=1.D0,RMIN=55.D0,
c     &          RMAX=RMIN+(JMAX-1)*RSTEP)
c      PARAMETER(IMAX=101,ESTEP=50.D0,EMIN=0.D0,
c     &          EMAX=EMIN+(IMAX-1)*ESTEP)
      DIMENSION RMATRIX(IMAX,JMAX)
      COMMON /LOCAL/ SIGMA,TMA,TMS,PI2
      SAVE IPASS,RMATRIX,TSAVE
      DATA IPASS /0/, TSAVE /-999./
      IF(IPASS.EQ.0 .OR. TNSL.NE.TSAVE) THEN
        IPASS=1
        TSAVE=TNSL
        sigma=5.d0
        PI2=8.d0*ATAN(1.d0)
        print *,'preparing ',real(rmin),real(rmax),real(rstep)
        print *,'          ',real(emin),real(emax),real(estep)
        DO J=1,JMAX
          RL=RMIN+(J-1)*RSTEP
          DO I=1,IMAX
            EL=EMIN+(I-1)*ESTEP
            Z=MAX(EL,20.d0*(RL-65.d0))
c*********************************************************************
c the offset (3.0) in the following is a fudge to get the present config to be
c stable. the temp TS passed out for the surface temperature does not
c include it....
            OFFSET=3.0d0
            TMA=49.13d0-0.007992d0*Z-0.7576d0*RL+TNSL+OFFSET
            TMS=30.78d0-0.006277d0*EL-0.3262d0*RL+TNSL+OFFSET 
C            TMA=49.13d0-0.007992d0*Z-0.7576d0*RL+TNSL
C            TMS=30.78d0-0.006277d0*EL-0.3262d0*RL+TNSL   
c*********************************************************************
C ... THIS IS HUYBRECTS, VERY SLOW...
            CALL QUAD2D(0.D0,365.D0,PDDH)
            PDDH=PDDH/SIGMA/SQRT(PI2)
c           print *,real(tma),real(tms),real(pddh),real(pddh*0.8*.01)
            RMATRIX(I,J)=PDDH
          ENDDO
        ENDDO
      ENDIF
      IF(ELEV1.GE.EMIN .AND. ELEV1.LE.EMAX .AND.
     &   RLAT.GE.RMIN  .AND. RLAT.LE.RMAX) THEN
        I=INT((ELEV1-EMIN)/ESTEP)+1
        I=MIN(I,IMAX)
        I=MAX(I,1)
        J=INT((RLAT-RMIN)/RSTEP)+1
        J=MIN(J,JMAX)
        J=MAX(J,1)
        ELEVI=-999.d0
        RLATI=-999.d0
        IF(I.EQ.IMAX) THEN
          I=I-1
          X=1.D0
        ELSE
          ELEVI=EMIN+(I-1)*ESTEP    
          X=(ELEV1-ELEVI)/ESTEP
        ENDIF
        IF(J.EQ.JMAX) THEN
          J=J-1
          Y=1.D0
        ELSE
          RLATI=RMIN+(J-1)*RSTEP    
          Y=(RLAT-RLATI)/RSTEP
        ENDIF
c       PRINT *,'I,J              ',I,J
c       PRINT *,'ELEV1,RLAT       ',ELEV1,RLAT
c       PRINT *,'ELEVI,RLATI      ',ELEVI,RLATI
c       PRINT *,'X,Y              ',X,Y
c       PRINT *,'(I,J),    (I+1,J)',RMATRIX(I,J),RMATRIX(I+1,J)
c       PRINT *,'(I+1,J+1),(I,J+1)',RMATRIX(I+1,J+1),RMATRIX(I,J+1)
        P1=(1.D0-X)*(1.D0-Y)
        P2=X*(1.D0-Y)
        P3=X*Y
        P4=(1.D0-X)*Y
        PDD=P1*RMATRIX(I,J)+P2*RMATRIX(I+1,J)+
     &      P3*RMATRIX(I+1,J+1)+P4*RMATRIX(I,J+1)
        Z=MAX(ELEV1,20.d0*(RLAT-65.d0))
c*********************************************************************
c these temps are passed out of the rouytine, and dont include the 
c OFFSET degree fudge included in the equations above
        TS=49.13d0-0.007992d0*Z-0.7576d0*RLAT+TNSL
        TSUMMER=30.78d0-0.006277d0*ELEV1-0.3262d0*rlat+TNSL
c*********************************************************************
      ELSE
        PRINT *,'PROBLEMS, OUT OF RANGE',elev1,rlat
        PDD=1000.d0
        Z=MAX(ELEV1,20.d0*(RLAT-65.d0))
        TS=49.13d0-0.007992d0*Z-0.7576d0*RLAT+TNSL
        TSUMMER=30.78d0-0.006277d0*ELEV1-0.3262d0*rlat+TNSL
C        PAUSE
      ENDIF
      END

c================================================
      subroutine huybant(elev1,rlat,pdd,TS,TSUMMER,TNSL)
      IMPLICIT REAL*8(A-H,O-Z)
c following covers greenland
      PARAMETER(JMAX=9,RSTEP=5.D0,RMIN=50.D0,
     &          RMAX=RMIN+(JMAX-1)*RSTEP)
c following for full n. hem
c      PARAMETER(JMAX=19,RSTEP=5.D0,RMIN=0.D0,
c     &          RMAX=RMIN+(JMAX-1)*RSTEP)
      PARAMETER(IMAX=26,ESTEP=200.D0,EMIN=0.D0,
     &          EMAX=EMIN+(IMAX-1)*ESTEP)
c      PARAMETER(JMAX=31,RSTEP=1.D0,RMIN=55.D0,
c     &          RMAX=RMIN+(JMAX-1)*RSTEP)
c      PARAMETER(IMAX=101,ESTEP=50.D0,EMIN=0.D0,
c     &          EMAX=EMIN+(IMAX-1)*ESTEP)
      DIMENSION RMATRIX(IMAX,JMAX)
      COMMON /LOCAL/ SIGMA,TMA,TMS,PI2
      SAVE IPASS,RMATRIX,TSAVE
      DATA IPASS /0/, TSAVE /-999./
      IF(IPASS.EQ.0 .OR. TNSL.NE.TSAVE) THEN
        IPASS=1
        TSAVE=TNSL
        sigma=5.d0
        PI2=8.d0*ATAN(1.d0)
        print *,'preparing ',real(rmin),real(rmax),real(rstep)
        print *,'          ',real(emin),real(emax),real(estep)
        DO J=1,JMAX
          RL=RMIN+(J-1)*RSTEP
          DO I=1,IMAX
            EL=EMIN+(I-1)*ESTEP
            Z=MAX(EL,20.d0*(RL-65.d0))
c*********************************************************************
            TMA=34.46d0-0.00914d0*EL-0.68775d0*RL+TNSL                 
            TMS=16.81d0-0.00692d0*EL-0.27973d0*RL+TNSL
c*********************************************************************
C ... THIS IS HUYBRECTS, VERY SLOW...
            CALL QUAD2D(0.D0,365.D0,PDDH)
            PDDH=PDDH/SIGMA/SQRT(PI2)
            RMATRIX(I,J)=PDDH
          ENDDO
        ENDDO
      ENDIF
      IF(ELEV1.GE.EMIN .AND. ELEV1.LE.EMAX .AND.
     &   RLAT.GE.RMIN  .AND. RLAT.LE.RMAX) THEN
        I=INT((ELEV1-EMIN)/ESTEP)+1
        I=MIN(I,IMAX)
        I=MAX(I,1)
        J=INT((RLAT-RMIN)/RSTEP)+1
        J=MIN(J,JMAX)
        J=MAX(J,1)
        ELEVI=-999.d0
        RLATI=-999.d0
        IF(I.EQ.IMAX) THEN
          I=I-1
          X=1.D0
        ELSE
          ELEVI=EMIN+(I-1)*ESTEP    
          X=(ELEV1-ELEVI)/ESTEP
        ENDIF
        IF(J.EQ.JMAX) THEN
          J=J-1
          Y=1.D0
        ELSE
          RLATI=RMIN+(J-1)*RSTEP    
          Y=(RLAT-RLATI)/RSTEP
        ENDIF
c       PRINT *,'I,J              ',I,J
c       PRINT *,'ELEV1,RLAT       ',ELEV1,RLAT
c       PRINT *,'ELEVI,RLATI      ',ELEVI,RLATI
c       PRINT *,'X,Y              ',X,Y
c       PRINT *,'(I,J),    (I+1,J)',RMATRIX(I,J),RMATRIX(I+1,J)
c       PRINT *,'(I+1,J+1),(I,J+1)',RMATRIX(I+1,J+1),RMATRIX(I,J+1)
        P1=(1.D0-X)*(1.D0-Y)
        P2=X*(1.D0-Y)
        P3=X*Y
        P4=(1.D0-X)*Y
        PDD=P1*RMATRIX(I,J)+P2*RMATRIX(I+1,J)+
     &      P3*RMATRIX(I+1,J+1)+P4*RMATRIX(I,J+1)
        Z=MAX(ELEV1,20.d0*(RLAT-65.d0))
c*********************************************************************
c these temps are passed out of the routine
        TS=34.46d0-0.00914d0*ELEV1-0.68775d0*rlat+TNSL                 
        TSUMMER=16.81d0-0.00692d0*ELEV1-0.27973d0*rlat+TNSL
c*********************************************************************
      ELSE
        PRINT *,'PROBLEMS, OUT OF RANGE',elev1,rlat
        PAUSE
      ENDIF
      END

C===============================================
      FUNCTION WARMING(DT)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      IF(DT.GT.0.0) THEN
        S=0.05D0
      ELSEIF(DT.GT.-10.0) THEN
        S=0.05D0-0.005D0*DT
      ELSE
        S=0.10D0
      ENDIF
      tmp=(1.D0+S)**DT
      WARMING=tmp
      END
C===============================================
      subroutine payne(X,Y,TS,ACC)
      IMPLICIT REAL*8(A-H,O-Z)                                          
c      data bmax /0.5d0/
c      data sb /.01d0/
c      data e /450.d0/
c      data tmin /248.15d0/
c      data st /0.01d0/
c      data xbar,ybar /850.d0,850.d0/
      save ipass,bmax,sb,e,tmin,st,xbar,ybar
      data ipass /0/
      if(ipass.eq.0) then
        ipass=1
        print *,' checking for payne data file... '
        read(89,*,end=999) bmax,sb,e,tmin,st,xbar,ybar
        print *,' found:',bmax,sb,e,tmin,st,xbar,ybar
        goto 1000
999       continue
          bmax=0.5d0
          sb=0.01d0
          e=450.d0
          tmin=248.15d0
          st=0.01d0
          xbar=750.d0
          ybar=750.d0
          print *,' not found:default values'
          print *,bmax,sb,e,tmin,st,xbar,ybar
          rewind 89
          write(89,*) bmax,sb,e,tmin,st,xbar,ybar
1000    continue
      endif
      acc=100.d0*min(bmax,sb*(e-sqrt((x-xbar)**2+(y-ybar)**2)))
      ts=tmin+st*(sqrt((x-xbar)**2+(y-ybar)**2))
      ts=ts-273.15d0
      end
C===============================================
      SUBROUTINE PAYNEOUT1(TIME,NUMNP,NUMEL,KX,X,Y,HTICE,DEPB,
     &                    TEMPA,VEL,FLOWA,
     &                    VOL,AREA,RMAF)
      IMPLICIT REAL*8(A-H,O-Z)   
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM,MMAX=40)
      common /timeout/ itimeout
      dimension x(nmax),y(nmax),depb(nmax),vel(nmax,3)
      dimension flowa(nmax),tempa(mmax,nmax),htice(nmax)
      dimension KX(NMAX,4)
c      return
      itimeout=30000
      itimeout=200000
      print *,'EISMINT DATA DUMP AT',itimeout
      tt=time*.001d0
c ... ume_q_             
      if(nint(time).eq.itimeout) then
c ... plan-form at 200kyr ...			.p.
c ..... x, y, (km) data						x,y
        do i=1,numnp
          xx=x(i)*0.001d0
          yy=y(i)*0.001d0
          thick=htice(i)-depb(i)
c ....... 50 thickness		(m)		.ume_q_p_tk.	htice-depb 
          write(50,1) xx,yy,thick
c ....... 51 basal temperature	(K)		.ume_q_p_tp.	tempa(mmax,i)
          write(51,1) xx,yy,tempa(mmax,i)
        enddo
        do i=1,numel
          xx=0.d0
          yy=0.d0
          thick=0.d0
          do j=1,4
            xx=xx+x(kx(i,j))*0.001d0
            yy=yy+y(kx(i,j))*0.001d0
            thick=thick+htice(kx(i,j))-depb(kx(i,j))
          enddo
          xx=0.25d0*xx
          yy=0.25d0*yy
          thick=0.25d0*thick
c ....... 65 velocity magnitude (m/yr)		.ume_q_p_vmag.	vel(,1)
          write(65,1) xx,yy,vel(i,1)
c ....... 66 horizontal flux(m**2/yr)		.ume_q_p_uvq.	vel(,1)*thick
          write(66,1) xx,yy,vel(i,1)*thick
        enddo
c
      endif
c ... global time series ...			.t.
c ..... time (kyr) data
c ..... 55 volume		(km**3)		.ume_q_t_vo.	vol
        write(55,1) tt,vol
c ..... 56 area 		(km**2)		.ume_q_t_ar.	area
        write(56,1) tt,area
c ..... 57 melt area fraction			.ume_q_t_fr.	rmaf
        write(57,1) tt,rmaf
c
c ... local time series ...			.t.   at stations 1,2,3,4,5
c ..... time (kyr) data
c ..... 58 thickness 		(m)		.ume_q_t_tk_1,2,3,4,5.
        write(58,1) tt,htice(1861)-depb(1861)
     &                ,htice(2471)-depb(2471)
     &                ,htice(2776)-depb(2776)
     &                ,htice(2109)-depb(2109)
     &                ,htice(2719)-depb(2719)
c ..... 59 basal temperature 	(K) 		.ume_q_t_tp_1,2,3,4,5.
        write(59,1) tt,tempa(mmax,1861)
     &                ,tempa(mmax,2471)
     &                ,tempa(mmax,2776)
     &                ,tempa(mmax,2109)
     &                ,tempa(mmax,2719)
1     format(1x,5(2x,1pg13.6))
      end
c==========================================================
      subroutine payneout2(time,ipt,mmax,tttt,xxx,at,wwww)
      IMPLICIT REAL*8(A-H,O-Z) 
      common /timeout/ itimeout
      dimension tttt(mmax),xxx(mmax),at(mmax),wwww(mmax) 
c      return
      if(nint(time).ne.itimeout) return       
c
c ... depth profiles ...			.d.   at stations 1,2,3,4,5
c ... height above bedrock (m), data
      write(60,*) mmax
      write(61,*) mmax
      write(62,*) mmax
      write(63,*) mmax
      write(64,*) mmax
      do i=mmax,1,-1
c ..... 60 temperature 		(K)		.ume_q_d_tp_1,2,3,4,5.
        write(60,1) ipt,xxx(i),tttt(i)
c ..... 61 velocity-x 		(m/yr)		.ume_q_d_uv_1,2,3,4,5.
        write(61,1) ipt,xxx(i),-99999.
c ..... 62 velocity-y 		(m/yr)		.ume_q_d_vv_1,2,3,4,5.
        write(62,1) ipt,xxx(i),-99999.
c ..... 63 velocity-z 		(m/yr)		.ume_q_d_wv_1,2,3,4,5.
        write(63,1) ipt,xxx(i),wwww(i)
c ..... 64 flow coefficient 	(Pa**-2 s**-1)	.ume_q_d_af_1,2,3,4,5.
        atiii=3.168876462d-23/(at(i)**3)
        write(64,1) ipt,xxx(i),atiii
      enddo
1     format(1x,i10,5(1x,1pg13.6))
      end
c==========================================================
      FUNCTION WOABLATE(X,Y,ELEV1,SLOPE1,SHAPE1,TNSL,TS)                      
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION QI(12),QS(12),TTT(12)                                   
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      COMMON /LOCAL/ SIGMA,TMAA,TMS,PI2
      DATA QI /960.d0,1036.d0,1200.d0,825.,330.d0,90.d0,150.d0,
     &         600.d0,1200.d0,1020.d0,    
     &         930.d0,850.d0/                                               
      DATA QS /-0.667d0,4.6d0,11.667d0,9.167d0,3.667d0,
     *         1.d0,1.667d0,6.667d0,12.d0,6.333d0,  
     &         0.333d0,-3.333d0/                                            
C     DATA AAA /-9.14/, BBB /-.68/, CCC /34.461/                        
C     DATA WWW /13.05/, XXX /.664/, ZZZ /2.608/                         
       DATA AAA / -9.62376690d0     /                                     
       DATA BBB /-0.546917617d0     /                                     
       DATA CCC /  24.9793854d0     /                                     
       DATA WWW /  19.1390686d0     /                                     
       DATA XXX / 0.922791243d0     /                                     
       DATA ZZZ /-0.738900483d0     /                                     
      AAA=ACOM 
C ... CALCULATE LATITUDE                                                    
C     CALL SETRIG                                                       
      CALL RECPOL(X-XPOLE,Y-YPOLE,RLAT,RLONG)  
C ... ELEVATION (KM)                                                        
      ELEV=ELEV1/1000.d0
C ... SLOPE (M/KM)                                                          
      SLOPE=SLOPE1*1000.d0
C ... SHAPE (M/KM/KM)                                                       
      SHAPE=SHAPE1*1000.d0*1000.d0
      SHAPE=0.d0
C ... CALCULATE SURFACE MEAN ANNUAL AIR TEMP                                
      TS=AAA*ELEV+BBB*RLAT+CCC+TNSL+TNSLBASE                                
C ... CALCULATE MEAN ANNUAL TEMP OF FREE ATMOSPHERE-ISOTHEMAL LAYER         
      TF=0.67d0*(TS+273.0d0)+88.9d0
C ... CALCULATE SATURATION VAPOR PRESSURE                                   
      TERM1=-9.09718d0*(273.16d0/TF-1.0d0)                                    
      TERM2=-3.56654d0*LOG10(273.16d0/TF)                                   
      TERM3=0.876793d0*(1.0d0-TF/273.16d0)+0.785835d0
      EXPON=TERM1+TERM2+TERM3                                           
      ES=10.d0**EXPON                                                     
C ... CALCULATE ACCUMULATION RATE (M/YR)                                    
      TERM1=WWW*ES                                                      
      TERM2=XXX*SLOPE                                                   

c **** <<< EXPERIMENTAL >>> **** turn off slope term ...
      TERM2=0
c **** <<< EXPERIMENTAL >>> ****

      TERM3=ZZZ                                                         
      TERM4=-15.276d0*SHAPE                                               
      ACC=max(0.d0,TERM1+TERM2+TERM3+TERM4)
C ... CALCULATE ABLATION                                                    
C ... CALCULATE ABLATION                                                    
      QY=0.d0
      DO I=1,12                                                      
        QY=QY+QI(I)-QS(I)*RLAT                                          
      ENDDO                                                          
      QY=QY/12.d0
      PDD=0.d0
      DO I=1,12                                                      
C       TTT(I)=TS+0.021d0*((QI(I)+QS(I)*RLAT)-QY)+8.954d0
        TTT(I)=TS+0.021d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C       TTT(I)=TS+0.018d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C       WRITE(19,*) I,TTT(I),QI(I)-QS(I)*RLAT,QY                        
        IF(TTT(I).GT.0.0) PDD=PDD+30.d0*TTT(I)                            
      ENDDO
      ABL=.6d0*PDD                                                        
c ... LINEAR DECREASE IN ABLATION FROM NOMINAL VALUE AT 200 M (0.2 KM)
C     DOWN TO ZERO AT SEA LEVEL (USING ELEV1, ELEV IN KM)
      FACT=MIN(1.0D0,ELEV/0.2D0)
      ABL=ABL*FACT
C ... CALCULATE NET ACCUMULATION                                            
      ACCNET=ACC-ABL                                                    
      WOABLATE=ACCNET*.01d0                                                  
      END                                                               
c---------------------------------------------------
      subroutine milankold(TIME,x,y,elev1,slope1,shape1,TNSL,tmean,
     &                     accum)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nyear1=-100,nyear2=100,nlat=19,nmonth=12)
      integer insol(nmonth,nlat,nyear1:nyear2)
      integer lat(nlat)
      dimension rmil(nmonth),rmil1(nmonth),rmil2(nmonth)
      character*80 line
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      DIMENSION QI(12),QS(12),TTT(12)                                   
      DATA QI /960.d0,1036.d0,1200.d0,825.,330.d0,90.d0,150.d0,
     &         600.d0,1200.d0,1020.d0,    
     &         930.d0,850.d0/                                               
      DATA QS /-0.667d0,4.6d0,11.667d0,9.167d0,3.667d0,
     *         1.d0,1.667d0,6.667d0,12.d0,6.333d0,  
     &         0.333d0,-3.333d0/                                            
      data ipass /0/, sigma /5.67e-5/,pi /3.1415927/
      data www /  19.1390686e0     /                                     
      data xxx / 0.922791243e0     /                                     
      data zzz /-0.738900483e0     /                                     
      save ipass,insol,iy,lat
      return
C ... CALCULATE LATITUDE                                                    
C     CALL SETRIG                                                       
c     CALL RECPOL(X-XPOLE,Y-YPOLE,RLAT,RLONG)  
      CALL RECPOL(X,Y,RLAT,RLONG)  
      if(ipass.eq.0) then
c ... done once first call ...
        ipass=1
c        open(23,file='milank.dat')
        do iyear=nyear2,nyear1,-1
          read(23,'(a)') line
          do ilat=1,nlat
            read(23,*) iy,lat(ilat),
     &                (insol(im,ilat,iyear),im=1,nmonth)
          enddo
        enddo
      endif
      ryear=TIME/1000.
      if(ryear.lt.nyear1) ryear=nyear1
      if(ryear.gt.nyear2) ryear=nyear2
c ..............................................
C ... CALCULATE LATITUDE (w/ pole shift .......                                                    
C     CALL SETRIG                                                       
      CALL RECPOL(X-XPOLE,Y-YPOLE,RLAT,RLONG)  
c     CALL RECPOL(X,Y,RLAT,RLONG)  
c ... elevation (km)                                                        
      elev=elev1/1000.d0
c ... slope (m/km)                                                          
      slope=slope1*1000.d0
c ... shape (m/km/km)                                                       
      shape=shape1*1000.d0*1000.d0
      shape=0.d0
      if(rlat.gt.90) rlat=90
      if(rlat.lt.-90) rlat=-90
      llat=int(rlat/10)*10
      ilat=1+(90-llat)/10
      rl=1+(90-rlat)/10
      ilat1=rl
      ilat2=rl+1
      if(ilat1.lt.1) ilat1=1
      if(ilat2.gt.nlat) ilat2=nlat
      iyear1=int(abs(ryear))
      iyear2=int(abs(ryear)+1)
      if(ryear.lt.0) then
        itemp=iyear1
        iyear1=-iyear2
        iyear2=-itemp
      endif
      if(iyear1.lt.nyear1) iyear1=nyear1
      if(iyear2.lt.nyear1) iyear2=nyear1
      if(iyear1.gt.nyear2) iyear1=nyear2
      if(iyear2.gt.nyear2) iyear2=nyear2
c ... interpolate latitude ...................
      fractl=(lat(ilat1)-rlat)/10
      do im=1,nmonth
        rmil1(im)=insol(im,ilat1,iyear1)+
     &    fractl*(insol(im,ilat2,iyear1)-insol(im,ilat1,iyear1))
        rmil2(im)=insol(im,ilat1,iyear2)+
     &    fractl*(insol(im,ilat2,iyear2)-insol(im,ilat1,iyear2))
      enddo
      if(.false.) then
        print *,lat(ilat1),rlat,lat(ilat2)
        print *,ilat1,rl,ilat2
        print *,fractl
        print 100,lat(ilat1),(insol(im,ilat1,iyear1),im=1,4)
        print 101,rlat,(rmil1(im),im=1,4)
        print 100,lat(ilat2),(insol(im,ilat2,iyear1),im=1,4)
        print *
        print 100,lat(ilat1),(insol(im,ilat1,iyear2),im=1,4)
        print 101,rlat,(rmil2(im),im=1,4)
        print 100,lat(ilat2),(insol(im,ilat2,iyear2),im=1,4)
100     format(1x,i6,12i6)
101     format(1x,f6.1,12f6.1)
      endif
c ............................................
c ... interpolate year ........................................
      fracty=(ryear-iyear1)
      do im=1,nmonth
        rmil(im)=rmil1(im)+fracty*(rmil2(im)-rmil1(im))
      enddo
      if(.false.) then
        print *,iyear1,ryear,iyear2
        print *,fracty
        print 200,iyear1,(rmil1(im),im=1,4)
        print 201,ryear,(rmil(im),im=1,4)
        print 200,iyear2,(rmil2(im),im=1,4)
200     format(1x,i6,12f6.1)
201     format(1x,f6.1,12f6.1)
c        pause
      endif
c ..............................................................
c ... means and standard deviation for year
      rmean=0.0
      do im=1,nmonth
        rmean=rmean+rmil(im)
      enddo
      rmean=rmean/nmonth
      rplus=0.0
      do im=1,nmonth
        rplus=rplus+(rmil(im)-rmean)**2
      enddo
      rplus=rplus/nmonth
      rplus=sqrt(rplus)
c ...............................................................
c ... convert langleys/day to watt/m**2, multiply by 0.4843 .....
      tmean=rmean*0.4843
      tplus=rplus*0.4843
c ... convert w/m**2 to erg/cm**2/s, multipy by 1000. ...........
      tmean=tmean*1000.
      tplus=tplus*1000.
c ... convert to temperature e=sigma*t**4, ......................
c ............  sigma=5.67e-5 erg/cm**2/deg**4/s ................
c ... tplus is temp with mean+stdev flux ........................
      tplus=((tmean+tplus/16)/sigma)**0.25
c ... tmean is temp with just mean ..............................
      tmean=(tmean/sigma)**0.25
c ... tplus is now amplitude of annual signal ...................
      tplus=tplus-tmean
c ... offset because stephan-boltzman is cold (greenhouse??) ....
      tmean=tmean+25.0-7.65
c ... convert to deg c ..........................................
      tmean=tmean-273.16
c ... lapse rate with elevation .................................
      tmean=tmean+acom*elev+TNSL+TNSLBASE
      TS=tmean
c ... rest is standard climatology ..............................
c     print *,'surface mean annual air temp=',tmean
c ... calculate mean annual temp of free atmosphere-isothemal layer ...         
      tf=0.67d0*(tmean+273.0d0)+88.9d0
c     print *,'temp free at-isothrmal layer=',tf
c ... calculate saturation vapor pressure .............................
      term1=-9.09718d0*(273.16d0/tf-1.0d0)
      term2=-3.56654d0*log10(273.16d0/tf)
      term3=0.876793d0*(1.0d0-tf/273.16d0)+0.785835d0
      expon=term1+term2+term3
      es=10.d0**expon
c     print *,'saturation vapor pressure=',es
c ... calculate accumulation rate (m/yr) ..............................
      term1=www*es
      term2=xxx*slope
      term3=zzz
      term4=-15.276d0*shape
      acc=max(0.d0,term1+term2+term3+term4)
C ... CALCULATE ABLATION THE OLD WAY ...
      if(.true.) then
        QY=0.d0
        DO I=1,12                                                      
          QY=QY+QI(I)-QS(I)*RLAT                                          
        ENDDO                                                          
        QY=QY/12.d0
        PDD=0.d0
        DO I=1,12                                                      
C         TTT(I)=TS+0.021d0*((QI(I)+QS(I)*RLAT)-QY)+8.954d0
          TTT(I)=TS+0.021d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C         TTT(I)=TS+0.018d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C         WRITE(19,*) I,TTT(I),QI(I)-QS(I)*RLAT,QY                        
          IF(TTT(I).GT.0.0) PDD=PDD+30.d0*TTT(I)                            
        ENDDO
      else
c ... sum up positive degree days for ablation calc..............
        pddnew=0
        do i=0,364
          temp=tmean+tplus*sin(2*pi*i/364.)
          if(temp.gt.0.) pddnew=pddnew+temp
        enddo
        pdd=pddnew
      endif
      abl=.6d0*pdd
c      print *,TS,pdd,pddnew
c ... calculate net accumulation ......................................
      accnet=acc-abl
c     print *,'net accumulation/ablation=',accnet
c ... put in m/yr .....................................................
      accum=accnet*.01d0
      abl=abl*.01d0
      end
c------------------------------------------------------------------
      subroutine milank(TIME,x,y,elev1,slope1,shape1,TNSL,tmean,accum)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nyear1=-100,nyear2=100,nlat=19,nmonth=12,amplify=0.2)
      parameter(rmult=4.,rparam=75)
      integer insol(nmonth,nlat,nyear1:nyear2)
      integer lat(nlat)
      dimension rmil(nmonth),rmil1(nmonth),rmil2(nmonth)
      dimension rmilba(nmonth)
      character*80 line
      logical dispon
      common /lapse/ acom,hmax,windir(2),xpole,ypole,abl,tnslbase
      DIMENSION QI(12),QS(12),TTT(12)                                   
      DATA QI /960.d0,1036.d0,1200.d0,825.,330.d0,90.d0,150.d0,
     &         600.d0,1200.d0,1020.d0,    
     &         930.d0,850.d0/                                               
      DATA QS /-0.667d0,4.6d0,11.667d0,9.167d0,3.667d0,
     *         1.d0,1.667d0,6.667d0,12.d0,6.333d0,  
     &         0.333d0,-3.333d0/                                            
      data ipass /0/,pi /3.1415927/, dispon /.false./
      data www /  19.1390686e0     /                                     
      data xxx / 0.922791243e0     /                                     
      data zzz /-0.738900483e0     /                                     
      save ipass,insol,iy,lat
      return
C ... CALCULATE LATITUDE                                                    
C     CALL SETRIG                                                       
c     CALL RECPOL(X-XPOLE,Y-YPOLE,RLAT,RLONG)  
      CALL RECPOL(X,Y,RLAT,RLONG)  
      if(ipass.eq.0) then
c ... done once first call ...
        ipass=1
c        open(23,file='milank.dat')
        do iyear=nyear2,nyear1,-1
          read(23,'(a)') line
          do ilat=1,nlat
            read(23,*) iy,lat(ilat),
     &                (insol(im,ilat,iyear),im=1,nmonth)
          enddo
        enddo
      endif
      ryear=TIME/1000.
      if(ryear.lt.nyear1) ryear=nyear1
      if(ryear.gt.nyear2) ryear=nyear2
c ..............................................
C ... CALCULATE LATITUDE (w/ pole shift .......                                                    
C     CALL SETRIG                                                       
      CALL RECPOL(X-XPOLE,Y-YPOLE,RLAT,RLONG)  
c     CALL RECPOL(X,Y,RLAT,RLONG)  
c ... elevation (km)                                                        
      elev=elev1/1000.d0
c ... slope (m/km)                                                          
      slope=slope1*1000.d0
c ... shape (m/km/km)                                                       
      shape=shape1*1000.d0*1000.d0
      shape=0.d0
      if(rlat.gt.90) rlat=90
      if(rlat.lt.-90) rlat=-90
      llat=int(rlat/10)*10
      ilat=1+(90-llat)/10
      rl=1+(90-rlat)/10
      ilat1=rl
      ilat2=rl+1
      if(.false.) then
        print *,lat(ilat1),rlat,lat(ilat2)
        print *,ilat1,rl,ilat2
      endif
      if(ilat1.lt.1) ilat1=1
      if(ilat2.gt.nlat) ilat2=nlat
      iyear1=int(abs(ryear))
      iyear2=int(abs(ryear)+1)
      if(ryear.lt.0) then
        itemp=iyear1
        iyear1=-iyear2
        iyear2=-itemp
      endif
      if(iyear1.lt.nyear1) iyear1=nyear1
      if(iyear2.lt.nyear1) iyear2=nyear1
      if(iyear1.gt.nyear2) iyear1=nyear2
      if(iyear2.gt.nyear2) iyear2=nyear2
      if(.false.) then
        print *,iyear1,ryear,iyear2
      endif
c ... interpolate latitude ...................
      fractl=(lat(ilat1)-rlat)/10
      do im=1,nmonth
        rmil1(im)=insol(im,ilat1,iyear1)+
     &    fractl*(insol(im,ilat2,iyear1)-insol(im,ilat1,iyear1))
        rmil2(im)=insol(im,ilat1,iyear2)+
     &    fractl*(insol(im,ilat2,iyear2)-insol(im,ilat1,iyear2))
        rmilba(im)=insol(im,ilat1,0)+
     &    fractl*(insol(im,ilat2,0)-insol(im,ilat1,0))
      enddo
      if(.false.) then
        print *,lat(ilat1),rlat,lat(ilat2)
        print *,ilat1,rl,ilat2
        print *,fractl
        print 100,lat(ilat1),(insol(im,ilat1,iyear1),im=1,4)
        print 101,rlat,(rmil1(im),im=1,4)
        print 100,lat(ilat2),(insol(im,ilat2,iyear1),im=1,4)
        print *
        print 100,lat(ilat1),(insol(im,ilat1,iyear2),im=1,4)
        print 101,rlat,(rmil2(im),im=1,4)
        print 100,lat(ilat2),(insol(im,ilat2,iyear2),im=1,4)
100     format(1x,i6,12i6)
101     format(1x,f6.1,12f6.1)
      endif
c ............................................
c ... interpolate year ........................................
      fracty=(ryear-iyear1)
      do im=1,nmonth
        rmil(im)=rmil1(im)+fracty*(rmil2(im)-rmil1(im))
      enddo
      if(.false.) then
        print *,iyear1,ryear,iyear2
        print *,fracty
        print 200,iyear1,(rmil1(im),im=1,4)
        print 201,ryear,(rmil(im),im=1,4)
        print 200,iyear2,(rmil2(im),im=1,4)
200     format(1x,i6,12f6.1)
201     format(1x,f6.1,12f6.1)
c        pause
      endif
c ..............................................................
c ... means and standard deviation for year
      rmeanno=average(nmonth,rmil)
      rmeanba=average(nmonth,rmilba)
c     if(dispon) print *,rlat,ryear,rmeanba,rmeanno,rmeanno-rmeanba
      rplusno=stdev(nmonth,rmil)
      rplusba=stdev(nmonth,rmilba)
      rmean=rmeanno-rmeanba
      rplus=rplusno-rplusba
c ... rmean is EXCESS insolation
c ... rmeanno is PAST insolation
c ... rmeanba is PRESENT insolation
c ... INCREASE PAST INSOLATION
      if(rmean.lt.0) then      
c ..... (Only increase colder)
        rmean=rmean*rmult
        rmeanno=rmeanba+rmean
      else
c ..... (Also increase warmer)
        rmean=rmean*rmult
        rmeanno=rmeanba+rmean
      endif
c ... adjust for latitude: > 45 too cold, <45 too warm
      if(rlat.gt.45) then
        rmeanno=rmeanno+rparam*((abs(rlat-45.))/45.)**0.5
      else
        rmeanno=rmeanno-rparam*((abs(rlat-45.))/45.)**0.5
      endif
c .....................................
c     if(dispon) print *,rlat,ryear,rplusba,rplus,rplus-rplusba
      tmeanno=stepbol(rmeanno)
      tplusno=stepbol(rmeanno+rplusno*amplify)
      tmeanba=stepbol(rmeanba)
      tplusba=stepbol(rmeanba+rplusba*amplify)
      if(dispon) print *,'lat,yr:',rlat,ryear
c     if(dispon) print *,rmil
      if(dispon) print *,'insolno',rmeanno,rplusno,rmeanno+rplusno
      if(dispon) print *,'insolba',rmeanba,rplusba,rmeanba+rplusba
      if(dispon) print *,'tempsno',tmeanno,tplusno-tmeanno,tplusno
      if(dispon) print *,'tempsba',tmeanba,tplusba-tmeanba,tplusba
c     pause
      tplusno=tplusno-tmeanno
      tplusba=tplusba-tmeanba
      rmean=rmeanno
      rplus=rplusno
      tmean=tmeanno
      tplus=tplusno
c ... rest is standard climatology ..............................
c ... lapse rate with elevation .................................
      tmean=tmean+acom*elev+TNSL+TNSLBASE
      if(dispon) print *,'surface mean annual air temp=',tmean,tplus
c ... calculate mean annual temp of free atmosphere-isothemal layer ...         
      tf=0.67d0*(tmean+273.0d0)+88.9d0
c     print *,'temp free at-isothrmal layer=',tf
c ... calculate saturation vapor pressure .............................
      term1=-9.09718d0*(273.16d0/tf-1.0d0)
      term2=-3.56654d0*log10(273.16d0/tf)
      term3=0.876793d0*(1.0d0-tf/273.16d0)+0.785835d0
      expon=term1+term2+term3
      es=10.d0**expon
c     print *,'saturation vapor pressure=',es
c ... calculate accumulation rate (m/yr) ..............................
      term1=www*es
      term2=xxx*slope
      term3=zzz
      term4=-15.276d0*shape
      acc=max(1d-9,term1+term2+term3+term4)
      if(dispon) print *,'accumulation rate=',acc
C ... CALCULATE ABLATION THE OLD WAY ...
c      if(.false.) then
        QY=0.d0
        DO I=1,12                                                      
          QY=QY+QI(I)-QS(I)*RLAT                                          
        ENDDO                                                          
        QY=QY/12.d0
        PDDM=0.d0
        DO I=1,12                                                      
C         TTT(I)=tmean+0.021d0*((QI(I)+QS(I)*RLAT)-QY)+8.954d0
          TTT(I)=tmean+0.021d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C         TTT(I)=tmean+0.018d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C         WRITE(19,*) I,TTT(I),QI(I)-QS(I)*RLAT,QY                        
          IF(TTT(I).GT.0.0) PDDM=PDDM+30.d0*TTT(I)                            
c          if(dispon) print *,i,TTT(I),pddm,max(0.0,TTT(i)*30*0.6)
        ENDDO
c      else
c ..... sum up positive degree days for ablation calc..............
        pdd=0
        do i=0,364,30
          temp=tmean+tplus*sin(2*pi*i/364.)
          if(temp.gt.0.) pdd=pdd+temp*30
c          if(dispon) print *,i,temp,pdd,max(0.0,temp*30*0.6)
        enddo
c      endif
      pdd=pddm
      abl=.6d0*pdd
      if(dispon) print *,'ablation rate=',abl,pdd,PDDM
c      print *,tmean,pdd,pddnew
c ... calculate net accumulation ......................................
      accnet=acc-abl
c     print *,'net accumulation/ablation=',accnet
c ... put in m/yr .....................................................
      accum=accnet*.01d0
      abl=abl*.01d0
c      if(dispon) pause
      end
c--------------------------------------------------
      function stepbol(rmean)
      IMPLICIT REAL*8(A-H,O-Z)
      data sigma /5.67d-5/
      isign=sign(1.d0,rmean)
c ... convert langleys/day to watt/m**2, multiply by 0.4843 .....
      tmean=abs(rmean)*0.4843d0
c ... convert w/m**2 to erg/cm**2/s, multipy by 1000. ...........
      tmean=tmean*1000.d0
c ... convert to temperature e=sigma*t**4, ......................
c ............  sigma=5.67e-5 erg/cm**2/deg**4/s ................
c ... tplus is temp with mean+stdev flux ........................
      albedo=0.12d0 ! dirt
      albedo=0.84d0 ! dry snow
      albedo=0.0d0 ! perfect black body
c ... tmean is temp with just mean ..............................
      tmean=(tmean*(1d0-albedo)/sigma)**0.25d0
c ... offset because stephan-boltzman is cold (greenhouse??) ....
      tmean=tmean+8.d0
c ... convert to deg c ..........................................
      tmean=tmean-273.16d0
      stepbol=isign*tmean
      end
c------------------------------------------------
      function average(n,r)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension r(n)
      sum=0.0d0
      do i=1,n
        sum=sum+r(i)
      enddo
      average=sum/n
      end
c------------------------------------------------
      function stdev(n,r)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension r(n)
      rmean=average(n,r)
      sum=0.0d0
      do i=1,n
        sum=sum+(r(i)-rmean)**2
      enddo
      sum=sum/n
      stdev=sqrt(sum)
      end
c------------------------------------------------
      function accmars0(ts0,elev,abl)
      implicit real*8(a-h,o-z)
c     DATA WWW /  19.1390686d0     /
      DATA WWW /  1.91390686d0     /
      DATA XXX / 0.922791243d0     /
      DATA ZZZ /-0.738900483d0     /
      toffset=60
      ts=ts0+toffset
      tmean=214d0
      pmean=5.6d0
      roverm=192d0
      grav=3.73
      p1=5.6d0
      z1=0d0
      tk=ts+273d0
      TTT=(TK+tmean)*0.5d0
      term=-grav*(elev-z1)/(roverm*TTT)
      pres=p1*exp(term)
C ... CALCULATE MEAN ANNUAL TEMP OF FREE ATMOSPHERE-ISOTHERMAL LAYER
      TF=0.67d0*(TS+273.0d0)+88.9d0
C ... CALCULATE SATURATION VAPOR PRESSURE
      es=satvap(tf)
c     es2=satvap(ts+273.d0)
      es2=es
C ... CALCULATE ACCUMULATION RATE (M/YR)
      TERM1=WWW*ES
      TERM2=XXX*SLOPE
      TERM3=ZZZ
      TERM4=-15.276d0*SHAPE
c
c **** <<< EXPERIMENTAL >>> **** turn off slope term ...
      TERM2=0
      TERM3=0
      TERM4=0
c **** <<< EXPERIMENTAL >>> ****
c
      ACC=max(0.d0,TERM1+TERM2+TERM3+TERM4)
C ... CALCULATE ABLATION
      factor=5
      aaa000=www
      wind=(20000-elev)/20000
      ABL=factor*wind*aaa000*es2/pres
C ... CALCULATE NET ACCUMULATION
      tmp=abl
      abl=acc*0.1d0
      acc=tmp*0.1d0
      ACCNET=ACC-ABL
      if(accnet.lt.0) then
        accnet=accnet*0.1
      endif
      accmars0= ACCNET
c     abl=acc-accnet
      return
      PRINT *,'-----------------------------'
      PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TK
      PRINT *,'Elevation                   =',elev
      PRINT *,'Pressure                    =',pres
      PRINT *,'Wind                        =',wind
      PRINT *,'TEMP FREE AT-ISOTHRMAL LAYER=',TF
      PRINT *,'SATURATION VAPOR PRESSURE   =',real(ES),real(es2)
      PRINT *,'ACCUMULATION                =',ACC*.01
      PRINT *,'ABLATION                    =',-ABL*.01
      PRINT *,'NET                         =',-ACCNET*.01
      end
      function satvap(tf)
      implicit real*8(a-h,o-z)
      TERM1=-9.09718d0*(273.16d0/TF-1.0d0)
      TERM2=-3.56654d0*LOG10(273.16d0/TF)
      TERM3=0.876793d0*(1.0d0-TF/273.16d0)+0.785835d0
      EXPON=TERM1+TERM2+TERM3
      satvap=10.d0**EXPON
      end
c------------------------------------------------
      function accmars1(ts0,elev,acc,abl)
c ... this is the original used in the paper ...
      implicit real*8(a-h,o-z)
c     DATA WWW /  19.1390686d0     /
      DATA WWW /  1.91390686d-0     /
      DATA XXX / 0.922791243d0     /
      DATA ZZZ /-0.738900483d0     /
      toffset=60
      toffset=0
      ts=ts0+toffset
      tmean=214d0
      pmean=5.6d0
      roverm=192d0
      grav=3.73
      p1=5.6d0
      z1=0d0
      tk=ts+273d0
      TTT=(TK+tmean)*0.5d0
      term=-grav*(elev-z1)/(roverm*TTT)
      pres=p1*exp(term)
C ... CALCULATE MEAN ANNUAL TEMP OF FREE ATMOSPHERE-ISOTHERMAL LAYER
      TF=0.67d0*(TS+273.0d0)+88.9d0
C ... CALCULATE SATURATION VAPOR PRESSURE
      es=satvap(tf)
      es2=satvap(ts+273.d0)
C ... CALCULATE ACCUMULATION RATE (M/YR)
      TERM1=WWW*ES
      TERM2=XXX*SLOPE
      TERM3=ZZZ
      TERM4=-15.276d0*SHAPE
c **** <<< EXPERIMENTAL >>> **** turn off slope term ...
      TERM2=0
      TERM3=0
      TERM4=0
c **** <<< EXPERIMENTAL >>> ****
      ACC=max(0.d0,TERM1+TERM2+TERM3+TERM4)
C ... CALCULATE ABLATION
      factor=10
      aaa000=www
      wind=(30000-elev)/30000
c     wind=1
      factor=10
c     pres=p1
c     print *,elev,ts,wind,p1/pres
      ABL=factor*wind*aaa000*es2*p1/pres
C ... CALCULATE NET ACCUMULATION
      acc=acc*0.1d0
      abl=abl*0.1d0
      if(.false.) then
        tmp=abl
        abl=acc
        acc=tmp
      endif
c     acc=es
c     abl=es2
      ACCNET=ACC-ABL
c     if(accnet.lt.0) then
c       accnet=accnet*0.1
c     endif
      accmars1= ACCNET
c     abl=acc-accnet
      return
      PRINT *,'-----------------------------'
      PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TK
      PRINT *,'Elevation                   =',elev
      PRINT *,'Pressure                    =',pres
      PRINT *,'Wind                        =',wind
      PRINT *,'TEMP FREE AT-ISOTHRMAL LAYER=',TF
      PRINT *,'SATURATION VAPOR PRESSURE   =',real(ES),real(es2)
      PRINT *,'ACCUMULATION                =',ACC*.01
      PRINT *,'ABLATION                    =',-ABL*.01
      PRINT *,'NET                         =',-ACCNET*.01
      end
c------------------------------------------------
      function accmars(ts0,elev0,acc,abl)
c ... this is the newest one ...
      implicit real*8(a-h,o-z)
      DATA WWW /  19.1390686d0     /
c     DATA WWW /  1.91390686d-0     /
      DATA XXX / 0.922791243d0     /
      DATA ZZZ /-0.738900483d0     /
      elev=elev0+3000
      abloffset=0.2    !  cm/yr (X10, mm/yr)
      toffset=60
      toffset=2.5
      ts=ts0+toffset
      tmean=214d0
      pmean=5.6d0
      roverm=192d0
      grav=3.74
      p1=5.6d0
      z1= 0
      tk=ts+273d0
      TTT=(TK+tmean)*0.5d0
      term=-grav*(elev-z1)/(roverm*TTT)
      pres=p1*exp(term)
C ... CALCULATE MEAN ANNUAL TEMP OF FREE ATMOSPHERE-ISOTHERMAL LAYER
      TF=0.67d0*(TS+273.0d0)+88.9d0
C ... CALCULATE SATURATION VAPOR PRESSURE
      es=satvap(tf)
      es2=satvap(ts+273.d0)
C ... CALCULATE ABLATION
      aaa000=www
      windf=(30000-elev)/30000
      windf=1
      factor=5
c     pres=p1
c     print *,elev,ts,windf,p1/pres,es,es2
      ABL=factor*windf*aaa000*es2*p1/pres
      ABL=factor*WWW*windf*ES2*p1/pres
c     ABL=ABL+(5e-6*elev)**2
c ... add 0.1 cm/yr (1 mm/yr) to ablation
      ABL=ABL+abloffset
c     ABL=ABL*wind(elev)
C ... CALCULATE ACCUMULATION RATE (M/YR)
      ACC=WWW*windf*ES
C ... CALCULATE NET ACCUMULATION
      acc=acc*0.01d0
      abl=abl*0.01d0
c     print *,elev,ts,p1,pres,acc,abl,acc-abl
      if(.false.) then
        tmp=abl
        abl=acc
        acc=tmp
      endif
c     acc=.1
c     abl=.05
      ACCNET=ACC-ABL
c     if(accnet.lt.0) then
c       accnet=accnet*0.1
c     endif
      accmars= ACCNET
c     print *,accnet,ts0
c     abl=acc-accnet
      return
      PRINT *,'-----------------------------'
      PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TK
      PRINT *,'Elevation                   =',elev
      PRINT *,'Pressure                    =',pres
      PRINT *,'Wind                        =',windf
      PRINT *,'TEMP FREE AT-ISOTHRMAL LAYER=',TF
      PRINT *,'SATURATION VAPOR PRESSURE   =',real(ES),real(es2)
      PRINT *,'ACCUMULATION                =',ACC*.01
      PRINT *,'ABLATION                    =',-ABL*.01
      PRINT *,'NET                         =',-ACCNET*.01
      pause
      end
      function wind(elev)
      implicit real*8(a-h,o-z)
      elevr=elev-10000
      if(elevr.gt.0) then
        wind=1+sign(1.d0,elevr)*(1e-2*elevr)**1
      else
        wind=1
      endif
      end 
