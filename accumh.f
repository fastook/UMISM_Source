      FUNCTION ACCUMH(X,Y,ELEV1,TNSL,TMA)                      
      IMPLICIT REAL*8(A-H,O-Z)
      external func
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      COMMON /LOCAL/ SIGMA,TMAA,TMS,PI2
      DATA A1 /34.46D0/, A2 /-0.00914D0/, A3 /-0.68775D0/
      DATA B1 /16.81D0/, B2 /-0.00692D0/, B3 /-0.27937D0/
      DATA C1 /0.78D0/, C2 /2.525D-2/, C3 /2.225D-4/
      sigma=5.d0
      PI2=8.*ATAN(1.)
      ACOM=A2
c      PRINT *,'X,Y,ELEV=',X,Y,ELEV1                                      
C CALCULATE LATITUDE                                                    
      CALL SETRIG                                                       
      CALL RECPOL(X,Y,RLAT,RLONG)                                       
c      PRINT *,'LATITUDE=',RLAT                                          
C CALCULATE SURFACE MEAN ANNUAL AIR TEMP                                
      TMA=A1+A2*ELEV1+A3*RLAT+TNSL+TNSLBASE
      TMAA=TMA                               
c      PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TMA                        
C CALCULATE SUMMER ANNUAL TEMP
      tms=B1+B2*ELEV1+B3*RLAT+TNSL+TNSLBASE                                
c      PRINT *,'TEMP SUMMER=',tms  
C CALCULATE ACCUMULATION RATE (M/YR)                                    
      ACC=C1+C2*TMA+C3*TMA**2                                       
c      PRINT *,'ACCUMULATION=',ACC                                       
C CALCULATE ABLATION 
      CALL QUAD2D(0.,365.,PDD)
      PDD=PDD/SIGMA/SQRT(PI2)
c      PRINT *,'PDD=',PDD
      ABL=.008*PDD                                                        
c      PRINT *,'ABLATION=',ABL                                           
C CALCULATE NET ACCUMULATION                                            
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
c      write(11,*) 'x,y,td',x,y,td(x),(-0.5*(Y-TD(X))**2/SIGMA)
c      write(11,*) 'func',Y*EXP(-0.5*(Y-TD(X))**2/SIGMA)
c      func=0
      temp=Y*EXP(-0.5*(Y-TD(X))**2/SIGMA)
c      if(temp.gt.0) then
        func=temp
c      endif
      END
C
      FUNCTION Y1(X)
      IMPLICIT REAL*8(A-H,O-Z)
      Y1=0.d0
      END
C
      FUNCTION Y2(X)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL TD
      COMMON /LOCAL/ SIGMA,TMA,TMS,PI2
      Y2=TD(X)+2.5*SIGMA
c      Y2=TD(X)
      END
C
      FUNCTION TD(T)
      IMPLICIT REAL*8(A-H,O-Z)
      INTRINSIC COS
      COMMON /LOCAL/ SIGMA,TMA,TMS,PI2
      TD=TMA+(TMA-TMS)*COS(PI2*T/365.)
      END
C
      SUBROUTINE QGAUSX(F,A,B,SS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(ITMAX=20)
      DIMENSION X(ITMAX),W(ITMAX)
      EXTERNAL F
      DATA TOL /1D-2/
C      SUMOLD=1E30
      N=ITMAX
C      DO N=ITMAX,ITMAX
        CALL GAULEG(A,B,X,W,N)
        SUM=0.
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
C      SUMOLD=1E30
      N=ITMAX
C      DO N=ITMAX,ITMAX
        CALL GAULEG(A,B,X,W,N)
        SUM=0.
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
C GIVEN THE LOWER AND UPPER LIMITS OF INTEGRATION X1 AND X2, AND
C GIVEN THE NUMBER OF NODAL POINTS DESIRED N, THIS ROUTINE RETURNS
C ARRAYS X AND W OF LENGTH N, CONTAINING THE NODES AND WEIGHTS OF
C THE GUASS-LEGENDRE N-POINT QUADRATURE FORMULA
C
C     HIGH PRECISION IS A GOOD IDEA FOR THIS ROUTINE
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X1,X2,X(N),W(N)
C     INCREASE EPS IF YOU DON'T HAVE THIS PRECISION
      PARAMETER (EPS=3.D-14)
C     THE ROOTS ARE SYMMETRIC IN THE INTERVAL, SO WE ONLY HAVE TO
C     FIND HALF OF THEM
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
C     LOOP OVER THE DESIRED ROOTS
      DO 12 I=1,M
        Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))
C       STARTING WITH THE ABOVE APPROXIMATION TO THE ITH ROOT,
C       WE ENTER THE MAIN LOOP OF REFINEMENT BY THE NEWTON'S METHOD
1       CONTINUE
          P1=1.D0
          P2=0.D0
C         LOOP UP THE RECURRENCE RELATION TO GET THE LEGENDRE
C         POLYNOMIAL EVALUATED AT Z
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11        CONTINUE
C         P1 IS NOW THE DESIRED LEGENDRE POLYNOMIAL. WE NEXT COMPUTE PP,
C         ITS DIRIVATIVE, BY A STANDARD RELATION INVOLVING ALSO P2,
C         THE POLYNIMIAL OF ONE LOWER ORDER
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
C         NEWTON'S METHOD
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
C       SCALE THE ROOT TO ITS DESIRED INTERVAL
        X(I)=XM-XL*Z
C       AND PUT IN ITS SYMMETRIC COUNTERPART
        X(N+1-I)=XM+XL*Z
C       COMPUTE THE WEIGHT
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
C       AND ITS SYMMETRIC COUNTERPART
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
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
