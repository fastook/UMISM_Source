C READS FROM STANDARD OUT3 AND PLOTS CONTOURS DIRECTLY
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM,NCOLOR=15,NPLOT=10000)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL FOUND
      CHARACTER HED*80,ANSW,JUNK*80
      DIMENSION XOUT(NPLOT),YOUT(NPLOT),IOUT(NPLOT)
      DIMENSION X(44,40),Y(44,40),Z(44,40)
      DIMENSION XX(NMAX),YY(NMAX),ZZ(NMAX),NNODE(NMAX),NBOUND(400)
      DIMENSION PSURF(NMAX),FRACT(NMAX),FLOWA(NMAX),HTICE(NMAX)
      DIMENSION ADOT(NMAX),WWW(NMAX),WRATE(NMAX),THK(NMAX)
      DIMENSION DEPB(NMAX),SLDGB(NMAX),AFUDGE(NMAX),TBED(NMAX)
      DIMENSION WTHICK(NMAX),BMELT(NMAX),WWWBASE(NMAX)
      DIMENSION BSURF(NMAX)
      DIMENSION XLINE(400),YLINE(400)
      DIMENSION KX(NMAX,4),IK(5)
      DIMENSION ICMAP(NCOLOR)
      DIMENSION INARAY(3)
      DIMENSION TIME(NPLOT)
      logical iflush
      common /flush/ iflush
      data iflush /.true./
      DATA ICMAP /1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/
      ipass=0
      IDISP=1
      II=1
      write(*,*) 'input node to plot'
      read(*,*) NUM
21    CONTINUE
C *** FOLLOWING IS LEVEL SPACING
      ISEG=1
3     CONTINUE
      if(iflush) call gflush
      WRITE(*,*) 'INPUT -N FOR DISPLAY STEP'
      WRITE(*,*) 'INPUT  0 FOR HTICE,'
      WRITE(*,*) '       1 FOR BDROCK,'
      WRITE(*,*) '       2 FOR ADOT,'
      WRITE(*,*) '       3 FOR FRACT'
      WRITE(*,*) '       4 FOR FLOW CONSTANT,'
      WRITE(*,*) '       5 FOR DIFFERENCE,'
      WRITE(*,*) '       6 FOR THICKNESS,'
      WRITE(*,*) '       7 FOR MARGIN,'
      WRITE(*,*) '       8 FOR PSURF,'
      WRITE(*,*) '       9 FOR ZONE,'
      WRITE(*,*) '      10 FOR SLIDING CONSTANT,'
      WRITE(*,*) '      11-FOR AFUDGE'
      WRITE(*,*) '      12-FOR TBED'
      WRITE(*,*) '      13-FOR BASAL MELT RATE'
      WRITE(*,*) '      14-FOR WATER THICKNESS '
      WRITE(*,*) '      15-FOR BEDROCK DEPRESSION '
      WRITE(*,*) '      16-FOR RATE OF BEDROCK DEPRESSION '
      WRITE(*,*) '      17-FOR BEDROCK DEPRESSION DIFFERENCE '
      WRITE(*,*) '      18-FOR DIFFERENCE BETWEEN FIRST AND REST '
      WRITE(*,*) '      99-QUIT'
      READ(*,*,END=9999) IPLOT
      IF(IPLOT.EQ.99) GOTO 9999
  350 FORMAT (I6)
      IF(IPLOT.LT.0) THEN
        IDISP=ABS(IPLOT)
        GOTO 3
      ENDIF
      PRINT *,' THIS IS IPLOT',IPLOT
      IF(IPLOT.EQ.0) THEN
        RMAX=5500.
        RMIN=.0001
c       RMIN=-3000.0001
c       RMIN=-2999.999
C       RMIN=.001
C       RMAX=1.
        RLEVSP=500.
      ELSEIF(IPLOT.EQ.1) THEN
        RMAX=6000.
        RMIN=-3000.0001
        RLEVSP=500.
      ELSEIF(IPLOT.EQ.2) THEN
        RMAX=4.
        RMIN=-4.
        RMIN=0.
C       RMIN=-.5
C       RMAX=.5
c        rmax=.2
C       RMIN=0.
C       RMAX=0.
        RLEVSP=.05
        RLEVSP=.2
      ELSEIF(IPLOT.EQ.3) THEN
        RMAX=1.
        RMIN=.01
c       RMAX=.01
        RLEVSP=.15
      ELSEIF(IPLOT.EQ.4) THEN
        RMAX=7.
        RMIN=.0
        IF(IPLOT.EQ.4) RLEVSP=.50
      ELSEIF(IPLOT.EQ.5) THEN
        RMAX=3000.
        RMIN=-650.
        IF(IPLOT.EQ.5) RLEVSP=100.
      ELSEIF(IPLOT.EQ.6) THEN
        RMAX=6000.
        RMIN=0.001
        IF(IPLOT.EQ.6) RLEVSP=500.
c        RMAX=600.
c        RMIN=0.001
c        IF(IPLOT.EQ.6) RLEVSP=50.
      ELSEIF(IPLOT.EQ.7) THEN
        RMAX=1.5
        RMIN=1.0000
        IF(IPLOT.EQ.7) RLEVSP=1.
      ELSEIF(IPLOT.EQ.8) THEN
        RMAX=5000.
        RMIN=.0001
        IF(IPLOT.EQ.8) RLEVSP=500.
      ELSEIF(IPLOT.EQ.9) THEN
        RMAX=-110.
        RMIN=-120.
        IF(IPLOT.EQ.9) RLEVSP=1.
      ELSEIF(IPLOT.EQ.10) THEN
        RMAX=.04
        RMIN=0.
        IF(IPLOT.EQ.10) RLEVSP=.002
      ELSEIF(IPLOT.EQ.11) THEN
        RMAX=3.
        RMIN=0.0000
        IF(IPLOT.EQ.11) RLEVSP=.5
      ELSEIF(IPLOT.EQ.12) THEN
        RMAX=0.
        RMIN=-50.0001
        IF(IPLOT.EQ.12) RLEVSP=5.
      ELSEIF(IPLOT.EQ.13) THEN
        RMAX=10.0
        RMIN=-10.000001
        RMIN=-0.000001
        IF(IPLOT.EQ.13) RLEVSP=1.
      ELSEIF(IPLOT.EQ.14) THEN
        RMAX=1.001
        RMIN=0.001
        IF(IPLOT.EQ.14) RLEVSP=1.
        RMAX=1.0001
        RMIN=0.0001
        IF(IPLOT.EQ.14) RLEVSP=.1
c        RMAX=0.10001
c        RMIN=0.00001
c        IF(IPLOT.EQ.14) RLEVSP=.01
      ELSEIF(IPLOT.EQ.15) THEN
        RMAX=-3000.00
        RMIN=150.00
        IF(IPLOT.EQ.15) RLEVSP=-100.
      ELSEIF(IPLOT.EQ.16) THEN
        RMAX=-127.5
        RMIN=127.5
        IF(IPLOT.EQ.16) RLEVSP=-5.
      ELSEIF(IPLOT.EQ.17) THEN
        RMAX=-1000.00
        RMIN=100.01
        IF(IPLOT.EQ.17) RLEVSP=-50.
      ELSEIF(IPLOT.EQ.18) THEN
        RMAX=3000.
        RMIN=-650.
        IF(IPLOT.EQ.18) RLEVSP=100.
      ENDIF
      ISTART=1
1001  FORMAT(16I6)
C READ INPUT HEADER
      READ(30,1000,END=99) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &              INTER,DT
1000  FORMAT (A80,/,7I6,F8.0)
C READ INPUT GRID, THINGS THAT NEVER CHANGE
      READ(31) HED
      READ(31) (KODE,I=1,NUMNP)
      READ(31) (XX(I),I=1,NUMNP)
      READ(31) (YY(I),I=1,NUMNP)
      READ(31) (PSURF(I),I=1,NUMNP)
      READ(31) (AJUNK,I=1,NUMNP)
      READ(31) (KX(I,1),KX(I,2),KX(I,3),KX(I,4),I=1,NUMEL)
      READ(31) (IB1,IB2,BJUNK,I=1,NUMGBC)
C READ INPUT DIFFER, THINGS THAT DIFFER FROM ONE GRID TO THE NEXT
      READ(32) HED
      READ(32) (AJUNK,I=1,NUMNP)
      READ(32) (FRACT(I),I=1,NUMNP)
      READ(32) (FLOWA(I),I=1,NUMNP)
      READ(32) (SLDGB(I),I=1,NUMNP)
      IREAD=0
1     CONTINUE
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME
        IF(IPLOT.EQ.0 .OR. IPLOT.EQ.1 .OR.IPLOT.EQ.2 .OR.
     &   IPLOT.EQ.5 .OR. IPLOT.EQ.6 .OR.IPLOT.EQ.7 .OR.
     &   IPLOT.EQ.9 .OR. IPLOT.EQ.3 .OR. IPLOT.EQ.4 .OR. 
     &   IPLOT.EQ.11 .OR. IPLOT.EQ.18) THEN
          READ(33,END=999) HED
          READ(33) (HTICE(I),I=1,NUMNP)
          IF(IREAD.EQ.1) THEN
            DO I=1,NUMNP
              BSURF(I)=HTICE(I)
            ENDDO
          ENDIF
          READ(33) (ADOT(I),I=1,NUMNP)
          READ(33) (DEPB(I),I=1,NUMNP)
          READ(33) (CONST,I=1,NUMEL)
          READ(33) (ACON,I=1,NUMEL)
          READ(34) (FRACT(I),I=1,NUMNP)
          READ(34) (FLOWA(I),I=1,NUMNP)
          READ(34) (SLDGB(I),I=1,NUMNP)
          READ(34) (AFUDGE(I),I=1,NUMNP)
        ELSEIF(IPLOT.EQ.12 .OR. IPLOT.EQ.13 .OR. IPLOT.EQ.14) THEN
          READ(36,END=999) HED
          READ(36) (TBED(I),I=1,NUMNP)
          READ(36) (BMELT(I),I=1,NUMNP)
          READ(36) (WTHICK(I),I=1,NUMNP)
        ELSEIF(IPLOT.EQ.15 .OR. IPLOT.EQ.16 .OR. IPLOT.EQ.17) THEN
          READ(37,END=999) HED
          READ(37) (WWW(I),I=1,NUMNP)
          READ(37) (THK(I),I=1,NUMNP)
          READ(37) (WRATE(I),I=1,NUMNP)
          IF(IREAD.EQ.2) THEN
            DO I=1,NUMNP
              WWWBASE(I)=WWW(I)
            ENDDO
          ENDIF
        ENDIF
        IREAD=IREAD+1
        read(hed,'(7x,f10.0)') time(IREAD)
        time(iread)=time(iread)*0.001
c       PRINT *,IREAD,time(IREAD)
C CHANGE PLOTTED VARIABBLE HERE
        IF(IPLOT.EQ.0) ZZ(IREAD)=HTICE(NUM)
        IF(IPLOT.EQ.1) ZZ(IREAD)=DEPB(NUM)
        IF(IPLOT.EQ.2) ZZ(IREAD)=ADOT(NUM)
        IF(IPLOT.EQ.3) ZZ(IREAD)=FRACT(NUM)
        IF(IPLOT.EQ.4) ZZ(IREAD)=FLOWA(NUM)
        IF(IPLOT.EQ.5) ZZ(IREAD)=HTICE(NUM)-PSURF(NUM)
        IF(IPLOT.EQ.6 .OR. IPLOT.EQ.7) THEN
          ZZ(IREAD)=HTICE(NUM)-DEPB(NUM)
          IF(DEPB(NUM).LT.0.) THEN
            FLOT=DEPB(NUM)*(1.-1.03/.917)
            IF(HTICE(NUM).LT.FLOT) ZZ(NUM)=0.
          ENDIF
          IF(ZZ(IREAD).LT.10.) ZZ(IREAD)=0.
        ENDIF
        IF(IPLOT.EQ.8) ZZ(IREAD)=PSURF(NUM)
        IF(IPLOT.EQ.9) ZZ(IREAD)=ADOT(NUM)
        IF(IPLOT.EQ.10) ZZ(IREAD)=SLDGB(NUM)
        IF(IPLOT.EQ.11) ZZ(IREAD)=AFUDGE(NUM)
        IF(IPLOT.EQ.12) ZZ(IREAD)=TBED(NUM)
        IF(IPLOT.EQ.13) ZZ(IREAD)=BMELT(NUM)*1000.
        IF(IPLOT.EQ.14) ZZ(IREAD)=WTHICK(NUM)
        IF(IPLOT.EQ.15) ZZ(IREAD)=WWW(NUM)
        IF(IPLOT.EQ.16) ZZ(IREAD)=WRATE(NUM)
        IF(IPLOT.EQ.17) ZZ(IREAD)=WWW(NUM)-WWWBASE(NUM)
        IF(IPLOT.EQ.18) ZZ(IREAD)=HTICE(NUM)-BSURF(NUM)
      goto 1
999   continue
      tmin=1e30
      tmax=-1e30
      zmin=1e30
      zmax=-1e30
      do i=1,iread
        tmin=min(tmin,time(i))
        tmax=max(tmax,time(i))
        zmin=min(zmin,zz(i))
        zmax=max(zmax,zz(i))
      enddo
      print *,tmin,tmax,zmin,zmax
      if(ipass.eq.0) then
        call grstrt(600,600)
        ipass=1
      endif
      call window(real(tmin),real(tmax),real(zmin),real(zmax))
      call move(real(time(1)),real(zz(1)))
      do i=2,iread
        call draw(real(time(i)),real(zz(i)))
      enddo
      do i=1,iread
        write(9,*) time(i),zz(i)
      enddo
      write(9,*) -99999.,0.
      if(iplot.eq.0) then
        WRITE(9,*) 'HTICE,'
      elseif(iplot.eq.1) then
        WRITE(9,*) 'BDROCK,'
      elseif(iplot.eq.2) then
        WRITE(9,*) 'ADOT,'
      elseif(iplot.eq.3) then
        WRITE(9,*) 'FRACT'
      elseif(iplot.eq.4) then
        WRITE(9,*) 'FLOW CONSTANT,'
      elseif(iplot.eq.5) then
        WRITE(9,*) 'DIFFERENCE,'
      elseif(iplot.eq.6) then
        WRITE(9,*) 'THICKNESS,'
      elseif(iplot.eq.7) then
        WRITE(9,*) 'MARGIN,'
      elseif(iplot.eq.8) then
        WRITE(9,*) 'PSURF,'
      elseif(iplot.eq.9) then
        WRITE(9,*) 'ZONE,'
      elseif(iplot.eq.10) then
        WRITE(9,*) 'SLIDING CONSTANT,'
      elseif(iplot.eq.11) then
        WRITE(9,*) 'AFUDGE'
      elseif(iplot.eq.12) then
        WRITE(9,*) 'TBED'
      elseif(iplot.eq.13) then
        WRITE(9,*) 'BASAL MELT RATE'
      elseif(iplot.eq.14) then
        WRITE(9,*) 'WATER THICKNESS '
      elseif(iplot.eq.15) then
        WRITE(9,*) 'BEDROCK DEPRESSION '
      elseif(iplot.eq.16) then
        WRITE(9,*) 'RATE OF BEDROCK DEPRESSION '
      elseif(iplot.eq.17) then
        WRITE(9,*) 'BEDROCK DEPRESSION DIFFERENCE '
      elseif(iplot.eq.18) then
        WRITE(9,*) ' DIFFERENCE BETWEEN FIRST AND REST '
      endif

200   FORMAT(I6,I4,2E12.5,F10.2,F7.2,F9.2,F10.3,F10.1,F10.2,F10.3)
      REWIND 11
      REWIND 30
      REWIND 31
      REWIND 32
      REWIND 33
      REWIND 34
      REWIND 36
      REWIND 37
      CALL MAKCUR
      GOTO 3
9999  CONTINUE
      if(ipass.eq.1) call grstop1
      STOP
99    continue
      END
      subroutine wait(n)
      do i=1,n
        x=sin(real(i))
      enddo
      end
