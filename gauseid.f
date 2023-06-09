      SUBROUTINE GAUSEID(NMAX,NZ,N,AA,KA,B,X)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION AA(NMAX,NZ),KA(NMAX,NZ+1),B(NMAX),X(NMAX)
      REAL*8 SUM
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      EPS=1D-6
      XMAX=-1d30
C ... A DIFFRENT FIRST GUESS ...
c      DO I=1,N
c        X(I)=B(I)/AA(I,1)
c      ENDDO
C ... --------------------------
      TOL=1d-8
      DO I=1,N
        JMAX=KA(I,NZ+1)
        SUM=B(I)
        DO J=2,JMAX
          SUM=SUM-AA(I,J)*X(KA(I,J))
        ENDDO
        XNEW=SUM/AA(I,1)
        XMAX=MAX(XMAX,ABS(X(I)-XNEW))
        X(I)=XNEW
      ENDDO
      ERRORG=RESID(NMAX,NZ,N,AA,KA,B,X)
c      print 24,0,errorg
      IF(ERRORG.EQ.0.) THEN
C ... RARE BUT POSSIBLE RETURN
        RETURN
      ENDIF
      IF(XMAX.EQ.0.) THEN
C ... RARE BUT POSSIBLE RETURN
        RETURN
      ENDIF
      ITMAX=50
      DO ITER=1,ITMAX
        XMAX=-1d30
        DO I=1,N
          JMAX=KA(I,NZ+1)
          SUM=B(I)
          DO J=2,JMAX
            SUM=SUM-AA(I,J)*X(KA(I,J))
          ENDDO
          XNEW=SUM/AA(I,1)
          XMAX=MAX(XMAX,ABS(X(I)-XNEW))
          X(I)=XNEW
        ENDDO
c       print *,iter,xmax
c       print 23,(x(i),i=1,n)
23    format(5g13.6)
        ERROR=RESID(NMAX,NZ,N,AA,KA,B,X)
        RATIO=ERROR/ERRORG
c        print 24,iter,error,ratio,xmax
24      format(i5,3g13.6)
        IF(ABS(RATIO).LT.EPS) THEN
c          PRINT 25,'A:CONVERGED',iter,error,ratio,xmax
25        format(a,i5,3g13.6)
          RETURN
        ENDIF
        IF(XMAX.LT.TOL) THEN
c          PRINT 25,'B:CONVERGED',iter,error,ratio,xmax
          RETURN
        ENDIF
      ENDDO
      PRINT *,'DIDNOT CONVERGE IN ',itmax,''
c      pause
      END
C=======================================
      FUNCTION RESID(NMAX,NZ,N,AA,KA,B,X)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION AA(NMAX,NZ),KA(NMAX,NZ+1),B(NMAX),X(NMAX)
      REAL*8 SUMSQ,SUM
      SUMSQ=0.D0
      DO I=1,N
        JMAX=KA(I,NZ+1)
        SUM=0.d0
        DO J=1,JMAX
          SUM=SUM+AA(I,J)*X(KA(I,J))
        ENDDO
        SUM=SUM-B(I)
        SUMSQ=SUMSQ+SUM**2
      ENDDO
      RESID=SUMSQ
      END
C=======================================
      SUBROUTINE CONJUG(NMAX,NZ,NUM,ERROR,AA,KA,B,X)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      include "parameter.h"
      PARAMETER(LMAX=MAXNUM)
      DIMENSION AA(NMAX,NZ),KA(NMAX,NZ+1),B(NMAX),X(NMAX)
      DIMENSION R(LMAX),V(LMAX),Z(LMAX)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
100   FORMAT(1X,A,(1X,T10,5G13.6))
      IF(NUM.GT.LMAX) THEN
        WRITE(*,*) 'INCREASE LMAX=',LMAX,' TO',NUM
        STOP
      ENDIF
      M=MAX(1000,1*NUM)
      EPS=0D0
      ITMIN=10
      DELTA=EPS
C ...
c ... output matrix, etc
c      write(93,*) num,nz
c      do i=1,num
c        write(93,*) (aa(i,j),j=1,nz)
c      enddo
c      do i=1,num
c        write(93,*) (ka(i,j),j=1,nz+1)
c      enddo
c      write(93,*) (b(i),i=1,num)
c      write(93,*) (x(i),i=1,num)
c      pause
c ... 
C ... FORM A*X IN R
C     CALL MULTB(MMAX,NWIDE,NUM,A,X,R)
c$doacross local(i,jmax,sum)
      DO I=1,NUM
        JMAX=KA(I,NZ+1)
        SUM=0.D0
        DO J=1,JMAX
          SUM=SUM+AA(I,J)*X(KA(I,J))
        ENDDO
        R(I)=SUM
      ENDDO
C ... FORM R=B-A*X, ALSO V=R, ALSO VV=SQRT(<V,V>), ALSO C=<R,R>
      VV=0.D0
      C=0.D0
      DO I=1,NUM
        R(I)=B(I)-R(I)
        V(I)=R(I)
        VV=VV+V(I)*V(I)
        C=C+R(I)*R(I)
      ENDDO
      VV=SQRT(VV)
c      PRINT 100,' R ',(R(N),N=1,NUM)
c      PRINT 100,' VBASE ',VV
c      PRINT 100,' DBASE ',C
      DBASE=C
      VBASE=VV
C ... BEGIN ITERATION LOOP ...
      DO K=1,M
c        PRINT 100,' V ',(V(N),N=1,NUM)
c        PRINT 100,' VV ',VV
C ..... IF SQRT(<V,V>) < DELTA EXIT
        IF(VV.LE.DELTA .and. k.gt.ITMIN) THEN
          IF(IOTOGG) THEN
            write(list(ipage+1),*) 
     &            '1A:SOLUTION, ITERATION=',K,VV
            ipage=ipage+1
          ENDIF
          RETURN
        ENDIF
        IF(VV/VBASE.LE.ERROR .and. k.gt.ITMIN) THEN
          IF(IOTOGG) THEN
            write(list(ipage+1),*)
     &        '1B:SOLUTION, ITERATION=',K,VV/VBASE
            ipage=ipage+1
          ENDIF
          RETURN
        ENDIF
C ..... FORM Z=A*V
C        CALL MULTB(MMAX,NWIDE,NUM,A,V,Z)
c$doacross local(i,jmax,sum)
        DO I=1,NUM
          JMAX=KA(I,NZ+1)
          SUM=0.D0
          DO J=1,JMAX
            SUM=SUM+AA(I,J)*V(KA(I,J))
          ENDDO
          Z(I)=SUM
        ENDDO
C        PRINT 100,' Z ',(Z(N),N=1,NUM)
C ..... FORM <V,Z> IN T
        T=0.0d0
        DO I=1,NUM
          T=T+V(I)*Z(I)
        ENDDO
C        PRINT 100,' <V,Z> ',T
C ..... PREVENT DIVIDE BY ZERO
        IF(T.EQ.0.0) THEN
          WRITE(*,*) 'DIVIDE BY ZERO IN CONJUG, ITERATION',K
          STOP
        ENDIF
C ..... FORM T=C/<V,Z>
        T=C/T
C        PRINT 100,' C/<V,Z> ',T
C ..... UNDATE X=X+T*V AND R=R-T*Z, FORM D=<R,R>
        D=0.D0
        DO I=1,NUM
          X(I)=X(I)+T*V(I)        
          R(I)=R(I)-T*Z(I) 
          D=D+R(I)*R(I)
        ENDDO
c        PRINT 100,' X ',(X(N),N=1,NUM)
c        PRINT 100,' R ',(R(N),N=1,NUM)
c        PRINT 100,' D ',REAL(K),D
C ..... IF D<EPS EXIT
        IF(D.LE.EPS .and. k.gt.ITMIN) THEN
          IF(IOTOGG) THEN
            write(list(ipage+1),1700) 
     &          '2A:SOLUTION, ITERATION=',K,REAL(D)
            ipage=ipage+1
          ENDIF
          RETURN
        ENDIF
        IF(D/DBASE.LE.ERROR .and. k.gt.ITMIN) THEN
          IF(IOTOGG) THEN
            write(list(ipage+1),1700)
     &            '2B:SOLUTION, ITERATION=',K,
     &             D,D/DBASE
            ipage=ipage+1
1700      format(1x,a,i5,1p2g13.6)
          ENDIF
          RETURN
        ENDIF
C ..... FORM V=R+(D/C)*V, ALSO VV=SQRT(<V,V>)
        CON=D/C
        VV=0.D0
        DO I=1,NUM
          V(I)=R(I)+CON*V(I)
          VV=VV+V(I)*V(I)
        ENDDO
        VV=SQRT(VV)
c        PRINT 100,' D,VV ',REAL(K),D,D/DBASE,VV,VV/VBASE
C        PRINT 100,' V ',(V(N),N=1,NUM)
C ..... REPLACE C WITH D
        C=D
C        PRINT 100,' K,C ',REAL(K),C
C        PAUSE
      ENDDO
      PRINT *,'DIDNOT CONVERGE IN ',M,''
      END
C=======================================
      SUBROUTINE MULTB(MMAX,NWIDE,NUM,A,X,R)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(MMAX,NWIDE),X(MMAX),R(MMAX)
      IB=NWIDE/2
      IB1=IB+1
      DO I=1,NUM
        C=0.D0
        DO K=1,NWIDE
          J=K+I-IB1
          IF(J.GE.1 .AND. J.LE.NUM) THEN
            C=C+A(I,K)*X(J)
          ENDIF
        ENDDO
        R(I)=C
      ENDDO
      END
C=======================================
      FUNCTION When ()
      REAL*8 When
      REAL tdum(2)
*
      When = dble(etime(tdum))
      END
C=======================================
      SUBROUTINE DIAGD(NMAX,NZ,NUM,AA,KA,LDIAG,DMAX,SMAX)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION AA(NMAX,NZ),KA(NMAX,NZ+1)
      LOGICAL LDIAG
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      LDIAG=.TRUE.
      DMAX=0.D0
      DO I=1,NUM
        JMAX=KA(I,NZ+1)
        DIAG=ABS(AA(I,1))
        SUM=0.0d0
        DO J=2,JMAX
          SUM=SUM+ABS(AA(I,J))
c          SUM=SUM+AA(I,J)
        ENDDO
        SUM=ABS(SUM)
        DIFF=ABS(DIAG-SUM)
c        IF(DIFF.GT.1E-6 .AND. DIAG.LT.SUM) THEN
        IF(DIAG.LE.SUM) THEN
          LDIAG=.FALSE.
          DMAX=DIAG
          SMAX=SUM
          IMAX=I
c          IF(IOTOGG) THEN
c            write(list(ipage+1),*) 
c     &          'NOT DIAGONALLY DOMINANT',REAL(DMAX),REAL(SMAX)
c            write(list(ipage+2),*)
c     &          '                       ',IMAX,REAL(DMAX-SMAX)
c            ipage=ipage+2
c          endif
c           print *, 
c     &          'NOT DIAGONALLY DOMINANT',REAL(DMAX),REAL(SMAX)
          RETURN
        ENDIF
        IF(DIAG.GT.DMAX) THEN
          DMAX=DIAG
          SMAX=SUM
          IMAX=I
        ENDIF
      ENDDO
c      IF(LDIAG) THEN
c        IF(IOTOGG) THEN
c          write(list(ipage+1),*)
c     &         'OK, DIAGONALLY DOMINANT',REAL(DMAX),REAL(SMAX)
c          ipage=ipage+1
c         PRINT *,'                       ',IMAX,REAL(DMAX-SMAX)
c        ENDIF
c          print *, 
c     &         'OK, DIAGONALLY DOMINANT',REAL(DMAX),REAL(SMAX)
c      ELSE
c        IF(IOTOGG) THEN
c          write(list(ipage+1),*)
c     &         'NOT DIAGONALLY DOMINANT',DMAX,SMAX
c          write(list(ipage+2),*)
c     &         '                       ',IMAX,DMAX-SMAX
c          ipage=ipage+2
c        ENDIF
c          print *, 
c     &         'NOT DIAGONALLY DOMINANT',DMAX,SMAX
c      ENDIF
      END
c=======================================
      SUBROUTINE ASYMSL(NMAX,NZ,NEQ,AA,KA,B,X,KKK)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION AA(NMAX,NZ),KA(NMAX,NZ+1),B(NMAX),X(NMAX)
      include "parameter.h"
      PARAMETER(MMAX=MAXNUM,MCOL=171)
      DIMENSION A(MCOL,MMAX),Q(MMAX)
      SAVE A,Q
      MBAND=(MCOL-1)/2
      NBAND=MBAND+1
      NCOLS=2*NBAND-1
      IF(KKK.EQ.1) THEN
C ..... COPY THE MATRIX TO BANDED FORM...
c        PRINT *,'COPY THE MATRIX TO BANDED FORM...'
        DO I=1,NEQ
          DO J=1,MCOL
            A(J,I)=0.D0
          ENDDO
          DO J=1,KA(I,NZ+1)
            JJ=KA(I,J)-I+1+MBAND
            IF(JJ.LE.0) PRINT *,JJ,' PROBLEMS...'
            IF(JJ.GT.MCOL) PRINT *,JJ,' PROBLEMS...'
            A(JJ,I)=AA(I,J)
          ENDDO
          Q(I)=B(I)
        ENDDO
c        PRINT *,'DONE COPYING...'
        KMIN=NBAND+1
        DO N=1,NEQ
          IF(A(NBAND,N).EQ.0.0) THEN
            WRITE(*,1001) N,A(NBAND,N)
            STOP
          ENDIF
          IF(A(NBAND,N).NE.1.0) THEN
            C=1.d0/A(NBAND,N)
            DO K=KMIN,NCOLS
              IF(A(K,N).NE.0.0) A(K,N)=C*A(K,N)
            ENDDO
          ENDIF
          DO L=2,NBAND
            JJ=NBAND-L+1
            I=N+L-1
            IF(I.LE.NEQ .AND. A(JJ,I).NE.0.0) THEN
              KI=NBAND+2-L
              KF=NCOLS+1-L
              J=NBAND
              DO K=KI,KF
                J=J+1
                A(K,I)=A(K,I)-A(JJ,I)*A(J,N)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
c        PRINT *,'DONE ELIMINATING...'
        RETURN
      ELSEIF(KKK.EQ.2) THEN
        DO N=1,NEQ
          IF(A(NBAND,N).EQ.0.0) THEN
            WRITE(*,1002) N,A(NBAND,N)
            STOP
          ENDIF
          IF(A(NBAND,N).NE.1.0) Q(N)=Q(N)/A(NBAND,N)
          DO L=2,NBAND
            JJ=NBAND-L+1
            I=N+L-1
            IF(I.LE.NEQ .AND. A(JJ,I).NE.0.0) Q(I)=Q(I)-A(JJ,I)*Q(N)
          ENDDO
        ENDDO
C
C BACK SUBSTITUTION
C
        LL=NBAND+1
        DO M=1,NEQ
          N=NEQ+1-M
          DO L=LL,NCOLS
            IF(A(L,N).NE.0.0) THEN
              K=N+L-NBAND
              Q(N)=Q(N)-A(L,N)*Q(K)
            ENDIF
          ENDDO
          X(N)=Q(N)
        ENDDO
c        PRINT *,'DONE BACK SUBSTITUTION...'
        RETURN
      ENDIF
1001  FORMAT(6X,'SET OF EQUATIONS MAY BE SINGULAR',
     &     /,5X,'DIAGONAL TERM OF EQUATION',I5,' IS EQUAL TO' ,E15.8)
1002  FORMAT(6X,'SET OF EQUATIONS ARE SINGULAR',
     &     /,5X,'DIAGONAL TERM OF EQUATION',I5,'  IS EQUAL TO',E15.8)
      END
C
C
