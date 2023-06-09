      SUBROUTINE SLOPEF(MXX,X,Y,KX,NUMEL,HTICE,SLOPE,WINDIR)
      IMPLICIT REAL*8(A-H,O-Z)
C ... CALCULATES LINEARIZATION CONSTANT FROM CURRENT SOLUTION
      DIMENSION SLOPE(4,MXX),HTICE(MXX)
      DIMENSION KX(MXX,4),X(MXX),Y(MXX)
      DIMENSION LM(5)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2)
      DIMENSION PSI(4),DPSI(4,2)
      DIMENSION XY(2,4),WINDIR(2)
      NODEN=4
      CENTX=0.0d0
      CENTY=0.0d0
      DO J = 1,NUMEL
        SUMX=0.d0
        SUMY=0.d0
        DO I = 1,NODEN
          LM(I) = KX(J,I)
        ENDDO
        I=LM(1)
        JJ=LM(2)
        K=LM(3)
        L=LM(4)
        XY(1,1)=X(I)
        XY(1,2)=X(JJ)
        XY(1,3)=X(K)
        XY(1,4)=X(L)
        XY(2,1)=Y(I)
        XY(2,2)=Y(JJ)
        XY(2,3)=Y(K)
        XY(2,4)=Y(L)
        CALL FESHAPE(1,CENTX,CENTY,PSI,DPSI)
C
C CALCULATE DXDS...EQUATION (5.3.6)
C
        DO I=1,2
          DO L=1,2
            DXDS(I,L)=0.0d0
            DO K=1,NODEN
              DXDS(I,L)=DXDS(I,L)+DPSI(K,L)*XY(I,K)
            ENDDO
          ENDDO
        ENDDO
C
C CALCULATE DSDX...EQUATION (5.2.7)
C
        DETJ=(DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1))
        IF (DETJ.LE.0.0) THEN
          WRITE(12,5544) J,DETJ
          WRITE(*,5544) J,DETJ
          WRITE(12,5545) (JJ,XY(1,JJ),XY(2,JJ),JJ=1,4)
          WRITE(*,5545) (JJ,XY(1,JJ),XY(2,JJ),JJ=1,4)
5545      FORMAT(1X,I5,1X,1PE10.3,E10.3)
          STOP
5544      FORMAT(' BAD JACOBIAN',I5,1PE10.3,/,1X,8E10.3)
        ENDIF
        DSDX(1,1)=DXDS(2,2)/DETJ
        DSDX(2,2)=DXDS(1,1)/DETJ
        DSDX(1,2)=-DXDS(1,2)/DETJ
        DSDX(2,1)=-DXDS(2,1)/DETJ
C
C CALCULATE D(PSI)/DX...EQUATION (5.3.5)
C
        DO I=1,NODEN
          DPSIX(I)=DPSI(I,1)*DSDX(1,1)+DPSI(I,2)*DSDX(2,1)
          DPSIY(I)=DPSI(I,1)*DSDX(1,2)+DPSI(I,2)*DSDX(2,2)
        ENDDO
        DO I = 1,NODEN
          SUMX = SUMX + HTICE(LM(I))*DPSIX(I)
          SUMY = SUMY + HTICE(LM(I))*DPSIY(I)
        ENDDO
C
        DELH = SUMX**2 + SUMY**2
        DELH = SQRT(DELH)
        SLOPE(1,J)=DELH
        SLOPE(2,J)=SUMX
        SLOPE(3,J)=SUMY
        SLOPE(4,J)=SUMX*SIN(WINDIR(1))+SUMY*COS(WINDIR(1))
      ENDDO
      RETURN
      END
