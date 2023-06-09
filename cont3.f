      include "parameter.h"
      PARAMETER(MXX=MAXNUM)
C **************************************************************C
C                                                               C
C   PROGRAM:  CONT3                                             C
C                                                               C
C   DATE:  11 23 87                                             C
C   PROGRAMMER:  FASTOOK                                        C
C                                                               C
C   FUNCTION:                                                   C
C            CONTOURS STAANDARD DATA SET (OUT3**) FOR SURFACE,  C
C            BED, FRACT, FLOWA, THICK ETC. OUTPUT IN STANDARD   C
C            PLOTTER FORM TO OUTC**                             C
C                                                               C
C **************************************************************C
      DIMENSION X(44,40),Y(44,40),Z(44,40)
      DIMENSION XX(MXX),YY(MXX),ZZ(MXX),NNODE(MXX),NBOUND(300)
      DIMENSION XLINE(300),YLINE(300)
      DIMENSION KX(MXX,4),IK(5)
      CHARACTER HED*80
      ISTART=0
C *** FOLLOWING IS LEVEL SPACING
      WRITE(*,*) 'INPUT -1 FOR DIFFERENT COLORS, INTEGER FOR UNIFORM'
      READ(*,*) IUNIFORM
      WRITE(*,*) 'INPUT 0 FOR HTICE'
      WRITE(*,*) '      1 FOR BDROCK'
      WRITE(*,*) '      2 FOR ADOT'
      WRITE(*,*) '      3 FOR FRACT'
      WRITE(*,*) '      4 FOR FLOW CONSTANT'
      WRITE(*,*) '      5 FOR DIFFERENCE'
      WRITE(*,*) '      6 FOR THICKNESS'
      WRITE(*,*) '      7 FOR MARGIN'
      WRITE(*,*) '      8 FOR PSURF'
      WRITE(*,*) '      9 FOR ZONE'
      WRITE(*,*) '     10 FOR SLIDING CONST'
      WRITE(*,*) '     11 FOR TEMPERATURE'
      READ(*,*) IPLOT
  350 FORMAT (I5)
      PRINT *,' THIS IS IPLOT',IPLOT
      RLEVSP=25.
      IF(IPLOT.EQ.0 .OR. IPLOT.EQ.8) RLEVSP=500.
      IF(IPLOT.EQ.1) RLEVSP=500.
      IF(IPLOT.EQ.2) RLEVSP=.2
      IF(IPLOT.EQ.9) RLEVSP=1.
      IF(IPLOT.EQ.3) RLEVSP=.2
      IF(IPLOT.EQ.4) RLEVSP=1.
      IF(IPLOT.EQ.10) RLEVSP=.005
      IF(IPLOT.EQ.11) RLEVSP=5.
      IF(IPLOT.EQ.5) RLEVSP=150.
      IF(IPLOT.EQ.6) RLEVSP=500.
      IF(IPLOT.EQ.7) RLEVSP=1.
      RNINE=-99999.
      ICOLOR=1
      RMINUS=-1.
      RZERO=0.
      RTWO=2.
      IOFF=-1
 1    READ(1,100,END=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &INTER,DT
      IOFF=IOFF+1
      PRINT *,HED
      ICOLOR=1
C     PRINT *,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
C **** THIS IS  PATCH ****
C     NUMGBC=0
100   FORMAT(A80,/,7I6,E15.6)
      NUM1=NUMCOL-1
      XMIN=1.E30
      XMAX=-XMIN
      YMIN=XMIN
      YMAX=XMAX
      DO NUM=1,NUMNP
        READ(1,200) N,KODE,XX(NUM),YY(NUM),HTICE,
     &            ADOT,FRACT,PSURF,BDROCK,FLOWA,SLDGB,TBED
c       xx(num)=-xx(num)
c       yy(num)=-yy(num)
        XMIN=MIN(XMIN,XX(NUM))
        YMIN=MIN(YMIN,YY(NUM))
        XMAX=MAX(XMAX,XX(NUM))
        YMAX=MAX(YMAX,YY(NUM))
C       CHANGE PLOTTED VARIABBLE HERE
        IF(IPLOT.EQ.0) ZZ(N)=HTICE
        IF(IPLOT.EQ.1) ZZ(N)=BDROCK
        IF(IPLOT.EQ.2) ZZ(N)=ADOT
        IF(IPLOT.EQ.9) ZZ(N)=ADOT
        IF(IPLOT.EQ.3) ZZ(N)=FRACT
        IF(IPLOT.EQ.4) ZZ(N)=FLOWA
        IF(IPLOT.EQ.10) ZZ(N)=SLDGB
        IF(IPLOT.EQ.5) ZZ(N)=HTICE-PSURF
        IF(IPLOT.EQ.6 .OR. IPLOT.EQ.7) THEN
          ZZ(N)=HTICE-BDROCK
          IF(BDROCK.LT.0.) THEN
c            FLOT=-BDROCK*1.03/.917
            FLOT=BDROCK*(1.-1.03/.917)
            IF(HTICE.LT.FLOT) ZZ(NUM)=0.
          ENDIF
        ENDIF
        IF(IPLOT.EQ.8) ZZ(N)=PSURF
        IF(IPLOT.EQ.11) ZZ(N)=TBED
      ENDDO
      REWIND 9
200   FORMAT(I6,I4,2E12.5,F10.2,F7.2,F9.2,F10.3,F10.1,F10.2,2F10.3)
      DO I=1,NUMEL
        READ(1,300) NUM,KX(NUM,1),KX(NUM,2),KX(NUM,3),KX(NUM,4)
        NNODE(NUM)=4
        IF(KX(NUM,4).EQ.0) NNODE(NUM)=3
      ENDDO
      DO N=1,NUMGBC
        READ(1,310) I,J,RJUNK
      ENDDO
310   FORMAT(2I6,E13.6)
300   FORMAT(5I6)
      IF(ISTART.EQ.0) THEN
        REWIND 11
        NLINE=0
        READ(11,*,END=9) NLINE
        READ(11,1001) (NBOUND(I),I=1,NLINE)
1001  FORMAT(16I6)
        DO 15 I=1,NLINE
C REVERSE THESE
C       WRITE(7,105) -YY(NBOUND(I))/1000.,XX(NBOUND(I))/1000.,I
        WRITE(7,105) XX(NBOUND(I))/1000.,YY(NBOUND(I))/1000.,I
15    CONTINUE
105   FORMAT(10X,G13.6,2X,G13.6,I13)
C     WRITE(7,105) RNINE,RTWO,ICOLOR
      WRITE(7,105) RNINE,RTWO,-1
      PRINT *,HED
        WRITE(7,*) HED
110   FORMAT(A80)
C       ISTART=1
      GOTO 91
9       CONTINUE
C REVERSE THESE
C       WRITE(7,105) -YMIN*.001,XMIN*.001
C       WRITE(7,105) -YMAX*.001,XMIN*.001
C       WRITE(7,105) -YMAX*.001,XMAX*.001
C       WRITE(7,105) -YMIN*.001,XMAX*.001
        WRITE(7,105) XMIN*.001,YMIN*.001
        WRITE(7,105) XMAX*.001,YMIN*.001
        WRITE(7,105) XMAX*.001,YMAX*.001
        WRITE(7,105) XMIN*.001,YMAX*.001
        WRITE(7,105) RNINE,RTWO,-1
        WRITE(7,*) HED
91      CONTINUE
      ENDIF
      RMAX=-1.D+30
      RMIN=1.D+30
      DO 40 I=1,NUMNP
      RMAX=MAX(ZZ(I),RMAX)
      RMIN=AMIN1(ZZ(I),RMIN)
40    CONTINUE
      WRITE(*,23) RMAX,RMIN
C     WRITE(7,111) RMIN,RMAX,RLEVSP
111   FORMAT('MAX=',G13.6,' MIN=',G13.6,' STEP=',G13.6)
23    FORMAT(2F10.3)
C     RMAX=AINT(RMAX)
C     RMIN=AINT(RMIN)
      WRITE(*,23) RMAX,RMIN
      IF(IPLOT.EQ.0 .OR. IPLOT.EQ.8) THEN
        RMAX=7500.
        RMIN=0.0001
      ENDIF
      IF(IPLOT.EQ.1) THEN
        RMAX=6000.
        RMIN=-1000.0001
        RMIN=-999.99999
      RMIN=-1499.9999
C     RMIN=.0001
C     RMAX=800.
C     RMIN=-500.000
      ENDIF
      IF(IPLOT.EQ.2) THEN
        RMAX=2.0
c        rmax=.1
c        RMIN=-.999
      RMIN=0.00
C     RMAX=-110.
C     RMIN=-120.
      ENDIF
      IF(IPLOT.EQ.9) THEN
      RMAX=-110.
      RMIN=-121.
      ENDIF
      IF(IPLOT.EQ.3) THEN
        RMAX=1.000001
      RMIN=0.000001
      rmin=0
      rmax=1.
      ENDIF
      IF(IPLOT.EQ.4) THEN
        RMAX=5.
        RMIN=1.
      ENDIF
      IF(IPLOT.EQ.10) THEN
        RMAX=.021
        RMIN=-.004
      ENDIF
      IF(IPLOT.EQ.11) THEN
        RMAX=0.0
        RMIN=-50.0001
      ENDIF
      IF(IPLOT.EQ.5) THEN
        RMAX=2100
        RMIN=-1200
      ENDIF
      IF(IPLOT.EQ.6) THEN
        RMAX=6000.
        RMIN=0.001
        rmin=1.0
c     RMIN=-999.999
      ENDIF
      IF(IPLOT.EQ.7) THEN
        RMAX=1.5
        RMIN=1.0000
      ENDIF
      LEVEL=(RMAX-RMIN)/RLEVSP+1
C     LEVEL=20
C     DLEV=(RMAX-RMIN)/(LEVEL-1)
      DLEV=RLEVSP
C **** REMOVE NEXT TWO
C     LEVEL=2
C     DLEV=1.
      ICOLOR=2
      DO 500 LEV=1,LEVEL
      ICOUNT=0
      VAL=RMIN+(LEV-1)*DLEV
      IVAL=0
      PRINT *,VAL
      DO 450 N=1,NUMEL
      DO 415 J=1,NNODE(N)
      IK(J)=KX(N,J)
415   CONTINUE
      IK(NNODE(N)+1)=KX(N,1)
      ICOUNT=0
      DO 425 NN=1,NNODE(N)
      I=IK(NN)
      J=IK(NN+1)
      IF(ZZ(I).LT.VAL .AND. VAL.LE.ZZ(J)) GOTO 430
      IF(ZZ(J).LT.VAL .AND. VAL.LE.ZZ(I)) GOTO 430
      GOTO 425
430   CONTINUE
      ICOUNT=ICOUNT+1
      XINT=XX(J)+(VAL-ZZ(J))*(XX(I)-XX(J))/(ZZ(I)-ZZ(J))
      YINT=YY(J)+(VAL-ZZ(J))*(YY(I)-YY(J))/(ZZ(I)-ZZ(J))
C     IF(ICOUNT.GT.2) THEN
C        IF(IVAL.EQ.1) WRITE(7,105) RNINE,RMINUS
C        IF(IVAL.EQ.0) WRITE(7,105) RNINE,RZERO
C        IVAL=1
C        WRITE(7,97) VAL
C        ICOUNT=0
C     ENDIF
C REVERSE THESE
C     WRITE(7,105) -YINT/1000.,XINT/1000.
      WRITE(7,105) XINT/1000.,YINT/1000.
425   CONTINUE
      IF(ICOUNT.EQ.0) GOTO 450
C REMOVE THIS
C     ICOLOR=2
C FOLLOWING FOR UNIFORM COLOR
      IF(IUNIFORM.NE.-1) ICOLOR=IUNIFORM+IOFF
      IF(IVAL.EQ.1) WRITE(7,105) RNINE,RMINUS,ICOLOR
      IF(IVAL.EQ.0) WRITE(7,105) RNINE,RTWO,ICOLOR
C     IF(IVAL.EQ.1) WRITE(7,105) RNINE,RMINUS,0
C     IF(IVAL.EQ.0) WRITE(7,105) RNINE,RZERO,0
      IVAL=1
      WRITE(7,97) VAL
450   CONTINUE
      ICOLOR=ICOLOR+1
      IF(ICOLOR.GT.14) ICOLOR=2
500   CONTINUE
99    FORMAT('NOT FOUND',F10.3)
98    FORMAT(3G13.4)
97    FORMAT('=',F10.3)
      IF(IPLOT.NE.7) ICOLOR=1
      GOTO 1
999   CONTINUE
      STOP
      END
