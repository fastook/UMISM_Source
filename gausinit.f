      SUBROUTINE GAUSINIT(XI,ETA,W)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XI(2,9), ETA(2,9), W(2,9)
C     GAUSSIAN QUADRATURE OF ORDER THREE QUADRILATERALS
      XI(1,1)=-SQRT(3.d0/5.d0)
      XI(1,2)=0.d0
      XI(1,3)=-XI(1,1)
      XI(1,4)=XI(1,1)
      XI(1,5)=0.d0
      XI(1,6)=XI(1,3)
      XI(1,7)=XI(1,1)
      XI(1,8)=0.d0
      XI(1,9)=XI(1,3)
      ETA(1,1)=XI(1,1)
      ETA(1,2)=XI(1,1)
      ETA(1,3)=XI(1,1)
      ETA(1,4)=0.d0
      ETA(1,5)=0.d0
      ETA(1,6)=0.d0
      ETA(1,7)=XI(1,3)
      ETA(1,8)=XI(1,3)
      ETA(1,9)=XI(1,3)
      W(1,1)=25.d0/81.d0
      W(1,2)=40.d0/81.d0
      W(1,3)=W(1,1)
      W(1,4)=W(1,2)
      W(1,5)=64.d0/81.d0
      W(1,6)=W(1,2)
      W(1,7)=W(1,1)
      W(1,8)=W(1,2)
      W(1,9)=W(1,1)
C   GAUSSIAN QUADRATURE OF ORDER THREE TRIANGLES
      XI(2,1)=1.d0/3.d0
      XI(2,2)=2.d0/15.d0
      XI(2,3)=XI(2,2)
      XI(2,4)=11.d0/15.d0
      ETA(2,1)=XI(2,1)
      ETA(2,2)=XI(2,4)
      ETA(2,3)=XI(2,2)
      ETA(2,4)=ETA(2,3)
      W(2,1)=-27.d0/96.d0
      W(2,2)=25.d0/96.d0
      W(2,3)=W(2,2)
      W(2,4)=W(2,2)
      END
