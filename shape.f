      SUBROUTINE FESHAPE(NTYPE,XI,ET,PSI,DPSI)
C ELEMENT SHAPE FUNCTIONS AND DERIVATIVES AT LOCAL COORDINATES (XI,ET)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PSI(4),DPSI(4,2)
      IF(NTYPE.EQ.1) THEN
        PSI(1)=.25d0*(1.d0-XI)*(1.d0-ET)
        PSI(2)=.25d0*(1.d0+XI)*(1.d0-ET)
        PSI(3)=.25d0*(1.d0+XI)*(1.d0+ET)
        PSI(4)=.25d0*(1.d0-XI)*(1.d0+ET)
        DPSI(1,2)=-.25d0*(1.d0-XI)
        DPSI(2,2)=-.25d0*(1.d0+XI)
        DPSI(3,2)=.25d0*(1.d0+XI)
        DPSI(4,2)=.25d0*(1.d0-XI)
        DPSI(1,1)=-.25d0*(1.d0-ET)
        DPSI(2,1)=.25d0*(1.d0-ET)
        DPSI(3,1)=.25d0*(1.d0+ET)
        DPSI(4,1)=-.25d0*(1.d0+ET)
      ELSE
        PSI(1)=1.d0-XI-ET
        PSI(2)=XI
        PSI(3)=ET
        DPSI(1,2)=-1.d0
        DPSI(2,2)=0.d0
        DPSI(3,2)=1.d0
        DPSI(1,1)=-1.d0
        DPSI(2,1)=1.d0
        DPSI(3,1)=0.d0
      ENDIF
      END
