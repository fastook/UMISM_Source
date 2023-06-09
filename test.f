      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /LOCAL/ SIGMA,TMA,TMS,PI2
      tnsl=0.0
c 1     print *,'elevation, latitude'
c        read(*,*,end=999) elev1,rlat
      do rlat=40.,90.,10.
      do elev1=0.,3000.,100.
        Z=MAX(ELEV1,20.*(RLAT-65))
        TS=49.13-0.007992*Z-0.7576*RLAT+TNSL
        TSUMMER=30.38-0.006277*ELEV1-0.3262*RLAT+TNSL
C ... THIS IS HUYBRECTS, VERY SLOW...
        sigma=5.d0
        PI2=8.*ATAN(1.)
        TMA=TS
        TMS=TSUMMER
        CALL QUAD2D(0.D0,365.D0,PDDH)
        PDDH=PDDH/SIGMA/SQRT(PI2)
        call huyb(elev1,rlat,pdd,TNSL)
        ABLH=0.8*PDDH
C ... THIS IS SIMPLE SINE, FASTER BUT TWO-TIMES TOO LARGE...
        PDD=0.
        AMP=TSUMMER-TS
        DO IT=1,365
          TD=TS+AMP*COS(PI2*IT/365.)
          IF(TD.GT.0.0) PDD=PDD+TD
        ENDDO
c ... this makes it agree better with huybrects' slower method...
c        PDD=PDD*0.5
        ABL=0.8*PDD
        print *,' elev,z               ',elev1,z
        print *,' annual temp          ',ts
        print *,' summer temp          ',tsummer
        print *,' positive degree days ',pddH,pdd
        print *,' ablation rate        ',ablH,abl
        write(11,*) rlat,elev1,ablH
        write(12,*) rlat,elev1,abl
      enddo
      enddo
c      goto 1
999   continue
      END                                                               
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
      temp=Y*EXP(-0.5*(Y-TD(X))**2/SIGMA)
      func=temp
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
      subroutine huyb(elev1,rlat,pdd,TS,TNSL)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(IMAX=50,JMAX=50,RMIN=0.D0,EMIN=0.D0)
      PARAMETER(RSTEP=5.D0,ESTEP=100.D0)
      DIMENSION RMATRIX(IMAX,JMAX)
      COMMON /LOCAL/ SIGMA,TMA,TMS,PI2
      DATA IPASS /0/
      SAVE IPASS,RMATRIX,RMAX,EMAX
      IF(IPASS.EQ.0) THEN
        IPASS=1
        sigma=5.d0
        PI2=8.*ATAN(1.)
        DO J=1,JMAX
          RL=RMIN+(J-1)*RSTEP
          DO I=1,IMAX
            EL=EMIN+(I-1)*ESTEP
            Z=MAX(EL,20.*(RL-65))
            TMA=49.13-0.007992*Z-0.7576*RL+TNSL
            TMS=30.38-0.006277*ELEV1-0.3262*RL+TNSL
C ... THIS IS HUYBRECTS, VERY SLOW...
            CALL QUAD2D(0.D0,365.D0,PDDH)
            RMATRIX(I,J)=PDDH
          ENDDO
        ENDDO
        RMAX=JMAX*RSTEP
        EMAX=IMAX*ESTEP
      ENDIF
      IF(ELEV1.GE.EMIN .AND. ELEV1.LE.EMAX .AND.
     &   RLAT.GE.RMIN  .AND. RLAT.LE.RMAX) THEN
        I=INT(ELEV1/ESTEP)
        I=MIN(I,IMAX-1)
        I=MAX(I,0)
        J=INT(RLAT/RSTEP)
        J=MIN(J,JMAX-1)
        J=MAX(J,0)
        ELEVI=I*ESTEP    
        RLATI=J*RSTEP    
        X=(ELEV1-ELEVI)/ESTEP
        Y=(RLAT-RLATI)/RSTEP
        P1=(1.D0-X)*(1.D0-Y)
        P2=X*(1.D0-Y)
        P3=X*Y
        P4=(1.D0-X)*Y
        PDD=P1*RMATRIX(I,J)+P2*RMATRIX(I+1,J)+
     &      P2*RMATRIX(I+1,J+1)+P4*RMATRIX(I,J+1)
        Z=MAX(EL,20.*(RL-65))
        TS=49.13-0.007992*Z-0.7576*RLAT+TNSL
      ELSE
        PRINT *,'PROBLEMS, OUT OF RANGE'
        PAUSE
      ENDIF
      END

