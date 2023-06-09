      SUBROUTINE FORMC(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &                 ETA, XI, W, CONST,
     &                 NUMGBC, IBFLUX, BFLUX, NZ, KZ, LM, AADOT,
     &                 ADOT,
     &                 D, B, A, CAP, KA,BOLD,ADIAG)
      IMPLICIT REAL*8(A-H,O-Z)
C FORM STIFFNESS MATRIX
      DIMENSION AADOT(NMAX),ADOT(NMAX)
      DIMENSION IBFLUX(NMAX,2)
      DIMENSION BFLUX(NMAX),KZ(NMAX)
      DIMENSION CONST(NMAX),LM(5),BOLD(NMAX),ADIAG(NMAX)
      DIMENSION KX(NMAX,4),X(NMAX),Y(NMAX),NTYPE(NMAX),D(NMAX),B(NMAX)
      REAL*8 A(NMAX,NZ),CAP(NMAX,NZ)
      INTEGER KA(NMAX,NZ+1)
      DIMENSION P(5),S(5,5),DD(5),CC(5,5)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2)
      DIMENSION PSI(4),DPSI(4,2)
      DIMENSION XY(2,4),XI(2,9),ETA(2,9),W(2,9)
C **********************************************************************
C     FORM CONDUCTIVITY MATRIX FOR COMPLETE BODY
C **********************************************************************
C
C ... ZERO OUT APPROPRIATE ARRAYS ...
c$doacross local(i,j)
      DO I=1,NUMNP
        KZ(I)=1
      ENDDO
      DO I=1,NUMNP
        KA(I,1)=I
      ENDDO
      DO I=1,NUMNP
        D(I)=0.0d0
      ENDDO
      DO I=1,NUMNP
        B(I)=0.0d0
      ENDDO
      DO J=1,NZ
        DO I=1,NUMNP
          A(I,J)=0.0d0
          CAP(I,J)=0.0d0
        ENDDO
      ENDDO
      DO J=1,NZ+1
        DO I=1,NUMNP
          KA(I,J)=0
        ENDDO
      ENDDO
C
C ... BEGIN LOOP OVER ALL THE ELEMENTS ...
      DO 100 N=1,NUMEL
        IF(NTYPE(N).EQ.1) THEN
          NODEN=4
          NINT=9
        ELSE
          NODEN=3
          NINT=4
        ENDIF
        DO I=1,4
          LM(I)=KX(N,I)
        ENDDO
C
C ..... FORM ELEMENT CONDUCTIVITY MATRIX
        DO I=1,4
          DD(I)=0.0d0
          P(I)=0.0d0
          DO J=1,4
            S(I,J)=0.0d0
            CC(I,J)=0.0d0
          ENDDO
        ENDDO
C
        I=LM(1)
        J=LM(2)
        K=LM(3)
        L=LM(4)
        XY(1,1)=X(I)
        XY(1,2)=X(J)
        XY(1,3)=X(K)
C       IF(NTYPE(N).EQ.1) XY(1,4)=X(L)
        XY(2,1)=Y(I)
        XY(2,2)=Y(J)
        XY(2,3)=Y(K)
C       IF(NTYPE(N).EQ.1) XY(2,4)=Y(L)
        IF(NODEN.EQ.4) THEN
          XY(1,4)=X(L)
          XY(2,4)=Y(L)
        ENDIF
C
C ..... FORM ELEMENT MATRIX AND VECTORS
C
C ..... BEGIN INTEGRATION POINT LOOP
        DO L=1,NINT
          CALL FESHAPE(NTYPE(N),XI(NTYPE(N),L),ETA(NTYPE(N),L),PSI,DPSI)
C
C ....... INTEGRATION POINT INTERPOLATION ...
          ADOTN=0.D0
          DO I=1,4
            ADOTN=ADOTN+ADOT(LM(I))*PSI(I)
          ENDDO
C
C ....... CALCULATE DXDS...EQUATION (5.3.6)
          DO I=1,2
            DO J=1,2
              DXDS(I,J)=0.0d0
              DO K=1,4
                DXDS(I,J)=DXDS(I,J)+DPSI(K,J)*XY(I,K)
              ENDDO
            ENDDO
          ENDDO
C
C ....... CALCULATE DSDX...EQUATION (5.2.7)
          DETJ=(DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1))
          IF(DETJ.LE.0.0) THEN
            WRITE(12,1100) DETJ,((XY(MM,NN)/1000,NN=1,4),MM=1,2)
            WRITE(*,1100) DETJ,((XY(MM,NN)/1000,NN=1,4),MM=1,2)
1100        FORMAT(' BAD JACOBIAN AT 161',G13.6,/,4G13.6,/,4G13.6)
            print *,n,lm
            STOP
          ENDIF
          DENOM=1.d0/DETJ
          DSDX(1,1)=DXDS(2,2)*DENOM
          DSDX(2,2)=DXDS(1,1)*DENOM
          DSDX(1,2)=-DXDS(1,2)*DENOM
          DSDX(2,1)=-DXDS(2,1)*DENOM
C
C ....... CALCULATE D(PSI)/DX...EQUATION (5.3.5)
          DO I=1,4
            DPSIX(I)=DPSI(I,1)*DSDX(1,1)+DPSI(I,2)*DSDX(2,1)
            DPSIY(I)=DPSI(I,1)*DSDX(1,2)+DPSI(I,2)*DSDX(2,2)
          ENDDO
C
C ....... ACCUMULATE INTEGRATION POINT VALUES OF INTEGRALS
          FAC=DETJ*W(NTYPE(N),L)
          DO I=1,4
C ......... THIS IS LUMPED CAPACITANCE MATRIX
            DD(I)=DD(I)+PSI(I)*FAC
C ........
C ........  ELEMENT AVERAGE VALUES FROM ELPROP ....
C           P(I)=P(I)+AADOT(N)*PSI(I)*FAC
C ........  INTEGRATION POINT VALUES FROM FEM INTERPOLATIOON ....
            P(I)=P(I)+ADOTN*PSI(I)*FAC
C
            DO J=1,4
              CC(I,J)=CC(I,J)+PSI(I)*PSI(J)*FAC
              IF(CONST(N).GT.1.E-30) THEN
                TERM1=CONST(N)*(DPSIX(I)*DPSIX(J)+DPSIY(I)*DPSIY(J))
                S(I,J)=S(I,J)+TERM1*FAC
              ENDIF
            ENDDO
          ENDDO
        ENDDO
c        write(73,*) 'qqq',n,real(const(n))
c        do i=1,4
c          write(73,*) (real(s(i,j)),j=1,4)
c        enddo
c        write(73,*) '------------------------------'
C
C ..... ADD ELEMENT CONDUCTIVITY TO COMPLETE CONDUCTIVITY MATRIX
C
        DO L=1,4
          I=LM(L)
C ....... THIS (D) IS LUMPED CAPACITANCE MATRIX ...
          D(I)=D(I)+DD(L)
C ....... THIS (B) IS THE LOAD VECTOR .............
          B(I)=B(I)+P(L)
C ....... THIS (A) IS THE STIFFNESS MATRIX ........
          DO M=1,4
            J=LM(M)
            IF(I.EQ.J) THEN
              A(I,1)=A(I,1)+S(L,M)
              CAP(I,1)=CAP(I,1)+CC(L,M)
              KA(I,1)=I
            ELSE
              DO K=2,KZ(I)
                IF(KA(I,K).EQ.J) THEN
                  A(I,K)=A(I,K)+S(L,M)
                  CAP(I,K)=CAP(I,K)+CC(L,M)
                  GOTO 99
                ENDIF
              ENDDO
              KZ(I)=KZ(I)+1
              A(I,KZ(I))=S(L,M)
              CAP(I,KZ(I))=CC(L,M)
              KA(I,KZ(I))=J
            ENDIF
99          CONTINUE
          ENDDO
        ENDDO
100   CONTINUE
C ... END LOOP OVER ALL THE ELEMENTS ...
C
C
C ... BOUNDARY CONDITIONS
C
C
c$doacross local(n)
      DO N=1,NUMNP
        BOLD(N)=B(N)
        ADIAG(N)=A(N,1)
      ENDDO
      IF(NUMGBC.GT.0) THEN
        DO N=1,NUMGBC
          I=IBFLUX(N,1)
          J=IBFLUX(N,2)
          TEMP=BFLUX(N)
          XL=SQRT((X(J)-X(I))**2+(Y(J)-Y(I))**2)
          TEMP=XL*TEMP*.5d0
          B(I)=B(I)+TEMP
          B(J)=B(J)+TEMP
        ENDDO
      ENDIF
c ... remove here, this writes out matrix
c      rewind 73
c      write(73,*) numnp
c      do i=1,numnp
c        write(73,*) kz(i)
c        do j=1,kz(i)
c          write(73,*) i,ka(i,j),a(i,j)
c        enddo
c      enddo
c      write(73,*) (b(i),i=1,numnp)
c ... to here ******
      END
