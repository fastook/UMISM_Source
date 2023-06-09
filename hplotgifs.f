      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      PARAMETER(NLINES=10,LINEM=1000)
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
      DIMENSION X(NMAX),Y(NMAX),HT(NMAX),NLINE(NLINES),
     &          LINE(NLINES,LINEM),ADOT(NMAX),FRACT(NMAX),
     &          PSURF(NMAX),BED(NMAX),FLOWA(NMAX),SLDGB(NMAX),
     &          THICK(NMAX)
      CHARACTER*80 HED
      character*4 cbase
      logical iflush
      common /flush/ iflush
      data iflush /.true./
      print *,'input base (1000 or greater)'
      read *,ibase
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
      READ(31) (PSURF(I),I=1,NUMNP)
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
      READ(32) HED
      READ(32) (AJUNK,I=1,NUMNP)
      READ(32) (FRACT(I),I=1,NUMNP)
      READ(32) (FLOWA(I),I=1,NUMNP)
      READ(32) (SLDGB(I),I=1,NUMNP)
c
c
c
      dowhile(.true.)
        ibase=ibase+1
        write(cbase,'(i4)') ibase
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME
        READ(33,END=999) HED
        READ(33) (HT(I),I=1,NUMNP)
        READ(33) (ADOT(I),I=1,NUMNP)
        READ(33) (BED(I),I=1,NUMNP)
        READ(33) (CONST,I=1,NUMEL)
        READ(33) (ACON,I=1,NUMEL)
        READ(34) (FRACT(I),I=1,NUMNP)
        READ(34) (FLOWA(I),I=1,NUMNP)
        READ(34) (SLDGB(I),I=1,NUMNP)
        READ(34) (AFUDGE,I=1,NUMNP)
        do i=1,numnp
          thick(i)=ht(i)-bed(i)
        enddo
        open(2,file='ps.data')
        DO I=1,NUMLIN
          DIST=0.0
          DELX=SQRT((X(LINE(1,1))-X(LINE(1,2)))**2+
     &            (Y(LINE(1,1))-Y(LINE(1,2)))**2)
          DO J=1,NLINE(I)
            L=LINE(I,J)
            IF(J.EQ.1) THEN
              XLAST=X(L)
              YLAST=Y(L)
            ENDIF
            DIST=DIST+SQRT((X(L)-XLAST)**2+(Y(L)-YLAST)**2)
            XLAST=X(L)
            YLAST=Y(L)
            WRITE(2,200) DIST/THOU,psurf(L),j
          enddo
c         WRITE(2,200) RNINE,RMINUS,ICOLOR
c         WRITE(2,211) 'PRESENT'
        enddo
        close(2)
        open(2,file='o2.data')
        open(3,file='o3.data')
        open(4,file='o4.data')
        DO I=1,NUMLIN
          DIST=0.0
          DELX=SQRT((X(LINE(1,1))-X(LINE(1,2)))**2+
     &            (Y(LINE(1,1))-Y(LINE(1,2)))**2)
          DO J=1,NLINE(I)
            L=LINE(I,J)
            IF(J.EQ.1) THEN
              XLAST=X(L)
              YLAST=Y(L)
            ENDIF
            DIST=DIST+SQRT((X(L)-XLAST)**2+(Y(L)-YLAST)**2)
            XLAST=X(L)
            YLAST=Y(L)
            WRITE(2,200) DIST/THOU,HT(L),j
            WRITE(3,200) DIST/THOU,BED(L)
            WRITE(4,200) DIST/THOU,ADOT(L)*1.
          enddo
c         WRITE(2,200) RNINE,RMINUS,ICOLOR
c         WRITE(3,200) RNINE,RMINUS,ICOLOR
c         WRITE(4,200) RNINE,RMINUS,ICOLOR
          ICOLOR=ICOLOR+1
          IF(ICOLOR.GT.8) ICOLOR=1
c         WRITE(2,211) HED
c         WRITE(3,211) HED
c         WRITE(4,211) HED
        enddo
        close(2)
        close(3)
        close(4)
        call strip(80,hed,n1)
c       call system('cat ps.data o2.data o3.data > asdf.data')
c       call system('plotbps-np.e asdf')
        open(9,file='gnu.bat')
        write(9,*) 'set terminal postscript linewidth 3'
        write(9,*) 'set key off'
        write(9,*) 'set grid'
        write(9,*) 'set title "'//HED(:20)//'"'
        write(9,*) 'set xlabel "x"'
        write(9,*) 'set ylabel "y"'
        write(9,*) 'set xrange [0:1000] writeback'
        !write(9,*) 'set yrange [-4000:4000] writeback'
        write(9,*) 'set yrange [1000:6000] writeback'
        write(9,*) 'set output "test.ps"'
        write(9,*) 'plot "o2.data" w l lt 1,\'
        write(9,*) '     "o3.data" w l lt 1'
        close(9)
        call system('chmod +x gnu.bat')
        call system('gnuplot gnu.bat')
        call system('convert -rotate 90 test.ps test.gif')
        call system('convert -rotate 90 test.ps test.pdf')
        !call system('convert -rotate 90 test.ps test.png')
        call system('cp test.gif Slice'//cbase//'.gif')
        call system('cp test.ps Slice'//cbase//'.ps')
        call system('cp test.pdf Slice'//cbase//'.pdf')
        !call system('cp test.png Slice'//cbase//'.png')
        print *,HED(:20)//cbase
      enddo
c---------------------
999   CONTINUE
c
124   FORMAT(3I5,5(1X,F10.2))
125   FORMAT(F10.0,I5,5(1X,F10.2))
123   FORMAT(F11.1,F7.0,F7.0,F7.3,I5,I5)
200   FORMAT(10X,G13.6,2X,G13.6,I13)
210   FORMAT('LINE - ',I3)
211   FORMAT(A80)
      STOP
      END
c---------------------------------------
      SUBROUTINE STRIP(N,JUNK,IC)
      CHARACTER*80 JUNK,JUNKN
      JUNKN=' '
      IC=0
      DO I=1,N
        IF(JUNK(I:I).NE.' ') THEN
          IC=IC+1
          JUNKN(IC:IC)=JUNK(I:I)
        ENDIF
      ENDDO
      JUNK=JUNKN
c      PRINT *,IC,JUNKN
      END
