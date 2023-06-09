      include "parameter.h"
      PARAMETER(NMAX=MAXNUM,NLINES=10,LINEM=1000)
      IMPLICIT REAL*8 (A-H,O-Z)
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
      DIMENSION X(NMAX),Y(NMAX),www(NMAX),NLINE(NLINES),
     &          LINE(NLINES,LINEM),thk(NMAX),wrate(nmax)
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
c
c
c
C READ INPUT HEADER
      READ(30,1000,END=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &              INTER,DT
1000  FORMAT (A80,/,7I6,F8.0)
C READ INPUT GRID, THINGS THAT NEVER CHANGE
      READ(31) HED
      READ(31) (KODE,I=1,NUMNP)
      READ(31) (X(I),I=1,NUMNP)
      READ(31) (Y(I),I=1,NUMNP)
      READ(31) (PSURF,I=1,NUMNP)
      READ(31) (AJUNK,I=1,NUMNP)
      READ(31) (KX1,KX2,KX3,KX4,I=1,NUMEL)
      READ(31) (IB1,IB2,BJUNK,I=1,NUMGBC)
      xmin=1e30
      xmax=-xmin
      DO I=1,NUMLIN
        L1=LINE(I,1)
        xy=0.0
        DO J=1,NLINE(I)
          L=LINE(I,J)
          xy=xy+sqrt((x(L1)-x(L))**2+(y(L1)-y(L))**2)
          xmax=max(xmax,xy)
          xmin=min(xmin,xy)
          L1=L
        enddo
      ENDDO

C READ INPUT DIFFER, THINGS THAT DIFFER FROM ONE GRID TO THE NEXT
c
c
c
      call grstrt(800,800)
      ymin=-5000.
      ymax=5000.
      call window(real(xmin),real(xmax),real(ymin),real(ymax))
      iseg=0
10    CONTINUE
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME
      READ(33,END=999) num
      READ(33) (www(I),I=1,NUMNP)
      READ(33) (thk(I),I=1,NUMNP)
      READ(33) (wrate(I),I=1,NUMNP)
      iseg=iseg+1
      call opnseg(iseg)
      CALL MOVE(REAL(XMIN),REAL(YMIN))
      CALL DRAW(REAL(XMIN),REAL(YMAX))
      CALL DRAW(REAL(XMAX),REAL(YMAX))
      CALL DRAW(REAL(XMAX),REAL(YMIN))
      CALL DRAW(REAL(XMIN),REAL(YMIN))
      DO I=1,NUMLIN
        call linclr(i)
        L1=LINE(I,1)
        call move(real(0.),real(www(L1)))
        xy=0.0
        DO J=1,NLINE(I)
          L=LINE(I,J)
          xy=xy+sqrt((x(L1)-x(L))**2+(y(L1)-y(L))**2)
          call draw(real(xy),real(www(L)))
          L1=L
        enddo
      ENDDO
      DO I=1,NUMLIN
        call linclr(2)
        L1=LINE(I,1)
        call move(real(0.),real(thk(L1)))
        xy=0.0
        DO J=1,NLINE(I)
          L=LINE(I,J)
          xy=xy+sqrt((x(L1)-x(L))**2+(y(L1)-y(L))**2)
          call draw(real(xy),real(thk(L)))
          L1=L
        enddo
      ENDDO
      DO I=1,NUMLIN
        call linclr(3)
        L1=LINE(I,1)
        call move(real(0.),real(40*wrate(L1)))
        xy=0.0
        DO J=1,NLINE(I)
          L=LINE(I,J)
          xy=xy+sqrt((x(L1)-x(L))**2+(y(L1)-y(L))**2)
          call draw(real(xy),real(40*wrate(L)))
          L1=L
        enddo
      ENDDO
      CALL CLOSEG
      if(iseg.gt.1) call setvis(iseg-1,.false.)
      DO 700 I=1,NUMLIN
        DIST=0.0
        DELX=SQRT((X(LINE(1,1))-X(LINE(1,2)))**2+
     &            (Y(LINE(1,1))-Y(LINE(1,2)))**2)
        DO 650 J=1,NLINE(I)
          L=LINE(I,J)
          IF(J.EQ.1) THEN
            XLAST=X(L)
            YLAST=Y(L)
          ENDIF
          DIST=DIST+SQRT((X(L)-XLAST)**2+(Y(L)-YLAST)**2)
          XLAST=X(L)
          YLAST=Y(L)
          WRITE(2,200) DIST/THOU,www(L),j
          WRITE(7,200) DIST/THOU,thk(L)
          XMAX=MAX(DIST,XMAX)
          XMIN=MIN(DIST,XMIN)
          YMAX=MAX(www(L),YMAX)
          YMIN=MIN(www(L),YMIN)
          YMAX=MAX(thk(L),YMAX)
          YMIN=MIN(thk(L),YMIN)
650     CONTINUE
        WRITE(2,200) RNINE,RMINUS,ICOLOR
        WRITE(7,200) RNINE,RMINUS,ICOLOR
        ICOLOR=ICOLOR+1
        IF(ICOLOR.GT.8) ICOLOR=1
        WRITE(2,211) HED
        WRITE(7,211) HED
210   FORMAT('LINE - ',I3)
700   CONTINUE
      GOTO 10
999   CONTINUE
      print *,'again??? input 1'
      read(*,*) iagain
      call setvis(iseg,.false.)
      if(iagain.gt.0) then
        call setvis(1,.true.)
        do i=2,iseg
          call setvis(i,.true.)
          do j=1,iagain
            call wait
          enddo
          call setvis(i-1,.false.)
        enddo
        goto 999
      endif
      do i=1,iseg
        call setvis(i,.true.)
        call wait
      enddo
c
124   FORMAT(3I5,5(1X,F10.2))
125   FORMAT(F10.0,I5,5(1X,F10.2))
123   FORMAT(F11.1,F7.0,F7.0,F7.3,I5,I5)
200   FORMAT(10X,G13.6,2X,G13.6,I13)
211   FORMAT(A80)
      call grstop
      STOP
      END
      subroutine wait
      do i=1,10000
        x=sin(1.)
      enddo
      end
