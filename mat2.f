      IMPLICIT REAL*8(A-H,O-Z)                                          
      PARAMETER(NMAX=29999,NZ=9,NZ1=NZ+1)
      DIMENSION AA(NMAX,NZ),KA(NMAX,NZ1),B(NMAX),X(NMAX)
      READ(1,*) N
      DO I=1,N
        READ(1,*) KA(I,NZ1)
        DO J=1,KA(I,NZ1)
          READ(1,*) I,KA(I,J),AA(I,J)
        ENDDO
        X(I)=0.D0
      ENDDO
      READ(1,*) (B(I),I=1,N)
c     do i=1,n
c       print 23,(aa(i,j),j=1,ka(i,nz1))
c       print 24,(ka(i,j),j=1,ka(i,nz1))
c       print 24,ka(i,nz1)
23    format(5g13.6)
24    format(5i13)
c     enddo
c     pause
      CALL GAUSEID(NMAX,NZ,N,AA,KA,B,X)
      END
      SUBROUTINE GAUSEID(NMAX,NZ,N,AA,KA,B,X)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION AA(NMAX,NZ),KA(NMAX,NZ+1),B(NMAX),X(NMAX)
      EPS=1D-6
      DO I=1,N
        JMAX=KA(I,NZ+1)
        SUM=B(I)
        DO J=2,JMAX
          SUM=SUM-AA(I,J)*X(KA(I,J))
        ENDDO
        X(I)=SUM/AA(I,1)
      ENDDO
      ERRORG=RESID(NMAX,NZ,N,AA,KA,B,X)
      IF(ERRORG.EQ.0.) THEN
C RARE BUT POSSIBLE RETURN
        RETURN
      ENDIF
      ITMAX=100
      DO ITER=1,ITMAX
        DO I=1,N
          JMAX=KA(I,NZ+1)
          SUM=B(I)
          DO J=2,JMAX
            SUM=SUM-AA(I,J)*X(KA(I,J))
          ENDDO
          X(I)=SUM/AA(I,1)
        ENDDO
c     print *,iter
c     print 23,(x(i),i=1,n)
23    format(5g13.6)
        ERROR=RESID(NMAX,NZ,N,AA,KA,B,X)
        RATIO=ERROR/ERRORG
        print 24,iter,error,errorg,ratio
24    format(i5,3g13.6)
        IF(ABS(RATIO).LT.EPS) THEN
          PRINT *,'CONVERGER'
          RETURN
        ENDIF
      ENDDO
      print *,'failed to conveger'
      END
      FUNCTION RESID(NMAX,NZ,N,AA,KA,B,X)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION AA(NMAX,NZ),KA(NMAX,NZ+1),B(NMAX),X(NMAX)
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

      
