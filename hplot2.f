      PARAMETER(NMAX=29999,LINEM=200)
C **************************************************************C
C                                                               C
C   PROGRAM:  HPLOT2                                            C
C                                                               C
C   DATE:  11 23 87                                             C
C   PROGRAMMER:  FASTOOK                                        C
C                                                               C
C   FUNCTION:                                                   C
C            PLOTS CONTOURS (SURFACE, BED VERSUS X (OR Y WITH   C
C            MINOR CHANGE)) FROM FILE LINE* WHICH CONTAINS      C
C            NODE NUMBERS TO BE PLOTTED, OUTPUT O2* (SURF)      C
C            AND O3* (BED) ARE STANDARD FORM FOR X-Y PLOTTER    C
C                                                               C
C **************************************************************C
      DIMENSION X(NMAX),Y(NMAX),HT(NMAX),NLINE(55),
     &          LINE(55,LINEM),ADOT(NMAX),FRACT(NMAX),
     &          PSURF(NMAX),BED(NMAX),FLOWA(NMAX),SLDGB(NMAX),
     &          THICK(NMAX)
      CHARACTER*80 HED
      THOU=1000.
      XMAX=-1.E30
      YMAX=-1.E30
      XMIN=1.E30
      YMIN=1.E30
      RNINE=-99999.
      ICOLOR=1
      RZERO=0.
      RMINUS=2.
      REWIND 11
      I=1
1     CONTINUE
      READ(11,*,END=9) NLINE(I)
      READ(11,*) (LINE(I,J),J=1,NLINE(I))
      write(*,*) (line(i,j),j=1,nline(i))
      I=I+1
      GOTO 1
9     NUMLIN=I-1
10    READ(1,1000,END=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &INTER,DT
      PRINT *,HED
1000  FORMAT(A80,/,7I6,E15.6)
      DO 300 I=1,NUMNP
      READ(1,1003) N,KODE,X(N),Y(N),HT(N),
     &ADOT(N),FRACT(N),PSURF(N),BED(N),FLOWA(N),SLDGB(N),THICK(N)
1003  FORMAT(I6,I4,2E12.5,F10.2,F7.2,F9.2,F10.3,F10.1,F10.2,2F10.3)
200   FORMAT(10X,G13.6,2X,G13.6,I13)
300   CONTINUE                                                                print *,'done with 300'
      DO 600 I=1,NUMEL
      READ (1,6345) JK1,JK2,JK3,JK4,JK5
 6345 FORMAT (5I6)
  600 CONTINUE
      print *,'done with 600'
      DO 610 I=1,NUMGBC
      READ(1,*) IJUNK,JJUNK
610   CONTINUE
      print *,'done reading data'		
      DO 700 I=1,NUMLIN
      DIST=0.0
      WRITE(23,*) 'TITLE'
      DELX=SQRT((X(LINE(1,1))-X(LINE(1,2)))**2+
     &          (Y(LINE(1,1))-Y(LINE(1,2)))**2)
      WRITE(23,*) DELX,1,1,1,1
      WRITE(23,*) 2.,.02,1.
      WRITE(32,*) HED
      WRITE(32,*) NLINE(I)
      WRITE(33,*) NLINE(I)
      DO 650 J=1,NLINE(I)
      L=LINE(I,J)
      IF(J.EQ.1) THEN
        XLAST=X(L)
        YLAST=Y(L)
      ENDIF
      DIST=DIST+SQRT((X(L)-XLAST)**2+(Y(L)-YLAST)**2)
      XLAST=X(L)
      YLAST=Y(L)
C     DIST=X(L)
C     DIST=Y(L)
C     WRITE(23,*) DIST/THOU,HT(L),BED(L),PSURF(L)
C     WRITE(23,*) DIST/THOU,HT(L),BED(L),THICK(L)
      WRITE(2,200) DIST/THOU,HT(L)
      WRITE(32,200) DIST/THOU,HT(L)
c      WRITE(7,200) DIST/THOU,ADOT(L)*1.
      WRITE(7,200) DIST/THOU,PSURF(L)
      XMAX=MAX(DIST,XMAX)
      XMIN=MIN(DIST,XMIN)
      WRITE(3,200) DIST/THOU,BED(L)
      WRITE(33,200) DIST/THOU,BED(L)
C     WRITE(23,123) DIST,HT(L),BED(L),1.,10,100
      WRITE(23,124) 1,J,NINT(100*FRACT(L)),BED(L),1.,ADOT(L)*100,
     &              HT(L),SQRT(X(L)**2+Y(L)**2)
C     WRITE(23,125) DIST,NINT(100*FRACT(L)),BED(L),1.,ADOT(L)*100,
C    &              HT(L),SQRT(X(L)**2+Y(L)**2)
124   FORMAT(3I6,5(1X,F10.2))
125   FORMAT(F10.0,I6,5(1X,F10.2))
123   FORMAT(F11.1,F7.0,F7.0,F7.3,I6,I6)
      YMAX=MAX(HT(L),YMAX)
      YMIN=MIN(HT(L),YMIN)
      YMAX=MAX(BED(L),YMAX)
      YMIN=MIN(BED(L),YMIN)
650   CONTINUE
      WRITE(2,200) RNINE,RMINUS,ICOLOR
C     WRITE(2,200) RNINE,RMINUS,0
      WRITE(7,200) RNINE,RMINUS,ICOLOR
      WRITE(3,200) RNINE,RMINUS,ICOLOR
      ICOLOR=ICOLOR+1
      IF(ICOLOR.GT.8) ICOLOR=1
      WRITE(2,*) 'surface:',HED(1:60)
      WRITE(7,211) 'present:',HED(1:60)
      WRITE(3,*) 'bed:',HED(1:60)
211   FORMAT(A80)
C     WRITE(2,210) I
C     WRITE(3,210) I
210   FORMAT('LINE - ',I3)
700   CONTINUE
      GOTO 10
999   CONTINUE
      WRITE(4,200) XMIN/THOU,YMIN
      WRITE(4,200) XMAX/THOU,YMIN
      WRITE(4,200) XMAX/THOU,YMAX
      WRITE(4,200) XMIN/THOU,YMAX
      WRITE(4,200) XMIN/THOU,YMIN
      WRITE(4,200) RNINE,RMINUS,-1
      WRITE(4,210)
      STOP
      END
