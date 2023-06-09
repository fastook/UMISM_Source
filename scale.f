      SUBROUTINE SCALEXY(xx,yy,n,xmin,xmax,ymin,ymax,delx,dely)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION xx(n),yy(n)
      logical zflag
      common /zoomcom/ pxmin,pymin,pxmax,pymax,zflag
      if(zflag) then
        XMIN=PXMIN
        YMIN=PYMIN
        XMAX=PXMAX
        YMAX=PYMAX
      else
        XMIN=1.d30
        YMIN=XMIN
        XMAX=-XMIN
        YMAX=-YMIN
        DO I=1,n
          XMAX=MAX(XMAX,XX(I))
          YMAX=MAX(YMAX,YY(I))
          XMIN=MIN(XMIN,XX(I))
          YMIN=MIN(YMIN,YY(I))
        ENDDO
      endif
      IF(XMAX-XMIN.GT.YMAX-YMIN) THEN
        YMAX=YMIN+XMAX-XMIN
      ELSE
        XMAX=XMIN+YMAX-YMIN
      ENDIF
      DELX=(XMAX-XMIN)/5.d0
      DELY=(YMAX-YMIN)/5.d0
      IF(XMIN.EQ.XMAX) XMAX=XMIN+1.d0
      IF(YMIN.EQ.YMAX) YMAX=YMIN+1.d0
      if(.not.zflag) then
        xmin=xmin*1.d-3
        xmax=xmax*1.d-3
        ymin=ymin*1.d-3
        ymax=ymax*1.d-3
        delx=delx*1.d-3
        dely=dely*1.d-3
      endif
c      print *,'x:',xmin,xmax,delx
c      print *,'y:',ymin,ymax,dely
c      pause
      END
c
      SUBROUTINE SCALE3(AR, RR, NN, II)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AR(NN+2)
      AMIN=1.d30
      AMAX=-1.d30
      DO I=1,NN
        AMAX=MAX(AR(I),AMAX)
        AMIN=MIN(AR(I),AMIN)
      ENDDO
      IF(AMIN.EQ.AMAX) AMAX=AMIN+1.d0
      AR(NN+1)=AMIN
      AR(NN+2)=(AMAX-AMIN)/RR
      AMAX=AMAX+50.d0*AR(NN+2)
      AMIN=AMIN-50.d0*AR(NN+2)
      AR(NN+1)=AMIN
      AR(NN+2)=(AMAX-AMIN)/RR
      RETURN
      END
C
      SUBROUTINE SCALE2(AR, RR, NN, II)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AR(NN+2)
      AMIN=1.d30
      AMAX=-1.d30
      DO I=1,NN
        AMAX=MAX(AR(I),AMAX)
        AMIN=MIN(AR(I),AMIN)
      ENDDO
      IF(AMIN.EQ.AMAX) AMAX=AMIN+1.d0
      AR(NN+1)=AMIN
      AR(NN+2)=(AMAX-AMIN)/RR
      RETURN
      END
C
      SUBROUTINE SORTIX(N,IX)
      DIMENSION IX(N)
      LOGICAL ANYEXC
      ANYEXC=.TRUE.
C WHILE ANY EXCHANGES
30    IF(ANYEXC) THEN
        ANYEXC=.FALSE.
        DO I=1,N-1
          IF(IX(I).LT.IX(I+1)) THEN
            CALL SWAP(IX(I),IX(I+1))
            ANYEXC=.TRUE.
          ENDIF
        ENDDO
        GOTO 30
      ENDIF
      RETURN
      END
C
      SUBROUTINE SWAP(IX1,IX2)
      ITEMP=IX1
      IX1=IX2
      IX2=ITEMP
      RETURN
      END
C
      SUBROUTINE ELMDUP(NN,IX)
      DIMENSION IX(NN)
      NNOUT=NN
      DO N=2,NN
80      IF(IX(N-1).EQ.IX(N) .AND. (IX(N).NE.0)) THEN
          NNOUT=NNOUT-1
          DO J=N,NN-1
            IX(J)=IX(J+1)
          ENDDO
          IX(NN)=0
          GOTO 80
        ENDIF
      ENDDO
      NN=NNOUT
      RETURN
      END
      INTEGER FUNCTION IFIND(X, Y, XC, YC)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XC(5),YC(5),TEST(4)
      XC(5)=XC(1)
      YC(5)=YC(1)
      DO I=1,4
        AX=X-XC(I)
        AY=Y-YC(I)
        BX=XC(I+1)-XC(I)
        BY=YC(I+1)-YC(I)
        TEST(I)=-(AX*BY-AY*BX)
      ENDDO
      IF((TEST(1).GE.0.) .AND.
     &   (TEST(2).GE.0.) .AND.
     &   (TEST(3).GE.0.) .AND.
     &   (TEST(4).GE.0.))THEN
        IFIND=1
      ELSE
        IFIND=0
      ENDIF
      RETURN
      END
