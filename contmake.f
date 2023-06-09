      include "parameter.h"
      PARAMETER(MXX=MAXNUM)
C **************************************************************C
C                                                               C
C   PROGRAM:  CONTMAKE                                          C
C                                                               C
C   DATE:  11 30 93                                             C
C   PROGRAMMER:  FASTOOK                                        C
C                                                               C
C   FUNCTION:                                                   C
C            CONTOURS STAANDARD DATA SET (OUT3**) FOR SURFACE,  C
C            BED, FRACT, FLOWA, THICK ETC. OUTPUT IN INFILE     C
C            PLOTTER FORM TO CONTOUR                            C
C                                                               C
C **************************************************************C
      DIMENSION X(44,40),Y(44,40),Z(44,40)
      DIMENSION XX(MXX),YY(MXX),ZZ(MXX),NNODE(MXX),NBOUND(300)
      DIMENSION XLINE(300),YLINE(300)
      DIMENSION KX(MXX,4),IK(5)
      CHARACTER HED*80
      ISTART=0
C *** FOLLOWING IS LEVEL SPACING
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
      READ(*,*) IPLOT
  350 FORMAT (I6)
      PRINT *,' THIS IS IPLOT',IPLOT
 1    READ(1,100,END=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &           INTER,DT
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
      xavg=0.
      yavg=0.
      DO NUM=1,NUMNP
        READ(1,200) N,KODE,XX(NUM),YY(NUM),HTICE,
     &              ADOT,FRACT,PSURF,BDROCK,FLOWA,SLDGB
         xavg=xavg+xx(num)
         yavg=yavg+yy(num)
c        xx(num)=1e-6*xx(num)
c        yy(num)=1e-6*yy(num)
        XMIN=MIN(XMIN,XX(NUM))
        YMIN=MIN(YMIN,YY(NUM))
        XMAX=MAX(XMAX,XX(NUM))
        YMAX=MAX(YMAX,YY(NUM))
C       CHANGE PLOTTED VARIABBLE HERE
        IF(IPLOT.EQ.0) ZZ(N)=1e-3*HTICE
        IF(IPLOT.EQ.1) ZZ(N)=1e-3*BDROCK
        IF(IPLOT.EQ.2) ZZ(N)=ADOT
        IF(IPLOT.EQ.9) ZZ(N)=ADOT
        IF(IPLOT.EQ.3) ZZ(N)=FRACT
        IF(IPLOT.EQ.4) ZZ(N)=FLOWA
        IF(IPLOT.EQ.10) ZZ(N)=SLDGB
        IF(IPLOT.EQ.5) ZZ(N)=1e-3*(HTICE-PSURF)
        IF(IPLOT.EQ.6 .OR. IPLOT.EQ.7) THEN
          ZZ(N)=1e-3*(HTICE-BDROCK)
          IF(BDROCK.LT.0.) THEN
c            FLOT=-BDROCK*1.03/.917
            FLOT=BDROCK*(1.-1.03/.917)
            IF(HTICE.LT.FLOT) ZZ(NUM)=0.
          ENDIF
        ENDIF
        IF(IPLOT.EQ.8) ZZ(N)=1e-3*PSURF
      ENDDO
      xavg=xavg/numnp
      yavg=yavg/numnp
      xmin=1e-3*(xmin-xavg)
      xmax=1e-3*(xmax-xavg)
      ymin=1e-3*(ymin-yavg)
      ymax=1e-3*(ymax-yavg)
200   FORMAT(I6,I4,2E12.5,F10.2,F7.2,F9.2,F10.3,F10.1,F10.2,F10.3)
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
      PRINT *,HED
C REVERSE THESE
      RMAX=-1.D+30
      RMIN=1.D+30
      DO I=1,NUMNP
        RMAX=MAX(ZZ(I),RMAX)
        RMIN=AMIN1(ZZ(I),RMIN)
      ENDDO
      WRITE(*,*) RMAX,RMIN
      WRITE(7,*) XMIN,XMAX,YMIN,YMAX
      WRITE(7,*) NUMCOL,NUMLEV
      WRITE(7,*) (ZZ(I),I=1,NUMNP)
      GOTO 1
999   CONTINUE
      STOP
      END
