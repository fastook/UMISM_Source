      PARAMETER(NMAX=79999,NLINES=10,LINEM=1000)
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
      CHARACTER*80 HED,junk
      logical iflush
      common /flush/ iflush
      data iflush /.true./
      print *,'input time you want to extract'
      read(*,*) jtime
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
c     call grstrt(800,800)
c arsia
      ymin=-1000.
      ymax=20000.
c deuter
c     ymin=-5000.
c     ymax=0
c      ymin=-500.
c      ymax=500.
c     call window(real(xmin),real(xmax),real(ymin),real(ymax))
      iseg=0
      ic=0
      iprint=0
10    CONTINUE
      if(iprint.eq.1) stop
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME
      READ(33,END=999) HED
      PRINT *,HED
      write(junk,'(a)') hed
      print *,'12345678901234567890'
      print *,junk
      read(junk,1234) ttime
1234  format(6x,f20.0)
      itime=nint(ttime)
c     print *,itime
      if(itime.eq.jtime) then
        iprint=1
      else
        iprint=0
      endif
c     print *,'iprint:',iprint,itime,jtime
      READ(33) (HT(I),I=1,NUMNP)
      READ(33) (ADOT(I),I=1,NUMNP)
      READ(33) (BED(I),I=1,NUMNP)
      READ(33) (CONST,I=1,NUMEL)
      READ(33) (ACON,I=1,NUMEL)
      READ(34) (FRACT(I),I=1,NUMNP)
      READ(34) (FLOWA(I),I=1,NUMNP)
      READ(34) (SLDGB(I),I=1,NUMNP)
      READ(34) (AFUDGE,I=1,NUMNP)
      ic=ic+1
c     if(mod(ic,10).ne.0) goto 10
      if(iprint.eq.0) goto 10
      iseg=iseg+1
c     call opnseg(iseg)
c     CALL MOVE(REAL(XMIN),REAL(YMIN))
c     CALL DRAW(REAL(XMIN),REAL(YMAX))
c     CALL DRAW(REAL(XMAX),REAL(YMAX))
c     CALL DRAW(REAL(XMAX),REAL(YMIN))
c     CALL DRAW(REAL(XMIN),REAL(YMIN))
c     DO I=1,NUMLIN
c       call linclr(i)
c       L1=LINE(I,1)
c       call move(real(0.),real(ht(L1)))
c       xy=0.0
c       DO J=1,NLINE(I)
c         L=LINE(I,J)
c         xy=xy+sqrt((x(L1)-x(L))**2+(y(L1)-y(L))**2)
c         call draw(real(xy),real(ht(L)))
c         L1=L
c       enddo
c     ENDDO
c     DO I=1,NUMLIN
c       call linclr(2)
c       L1=LINE(I,1)
c       call move(real(0.),real(psurf(L1)))
c       xy=0.0
c       DO J=1,NLINE(I)
c         L=LINE(I,J)
c         xy=xy+sqrt((x(L1)-x(L))**2+(y(L1)-y(L))**2)
c         call draw(real(xy),real(psurf(L)))
c         L1=L
c       enddo
c     ENDDO
c     DO I=1,NUMLIN
c       call linclr(i)
c       L1=LINE(I,1)
c       call move(real(0.),real(bed(L1)))
c       xy=0.0
c       DO J=1,NLINE(I)
c         L=LINE(I,J)
c         xy=xy+sqrt((x(L1)-x(L))**2+(y(L1)-y(L))**2)
c         call draw(real(xy),real(bed(L)))
c         L1=L
c       enddo
c     ENDDO
c     CALL CLOSEG
c     if(iseg.gt.1) call setvis(iseg-1,.false.)
      do i=1,numnp
        thick(i)=ht(i)-bed(i)
      enddo
      open(2,file='fort.2')
      open(3,file='fort.3')
      DO 700 I=1,NUMLIN
        write(24,201) hed(7:19),ht(line(i,1))
        vol=0.
        do j=1,nline(i)-1
          vol=vol+0.5*(x(line(i,j+1))-x(line(i,j)))*
     &                (ht(line(i,j+1))+ht(line(i,j)))
        enddo
        write(25,201) hed(7:19),vol/1e9
201     format(10x,a13,2x,g13.6)
        DIST=0.0
        WRITE(23,*) 'TITLE'
        DELX=SQRT((X(LINE(1,1))-X(LINE(1,2)))**2+
     &            (Y(LINE(1,1))-Y(LINE(1,2)))**2)
        WRITE(23,*) DELX,1,1,1,1
        WRITE(23,*) 2.,.02,1.
        WRITE(26,*) hed
        write(26,*) nline(i)
        DO 650 J=1,NLINE(I)
          L=LINE(I,J)
          IF(J.EQ.1) THEN
            XLAST=X(L)
            YLAST=Y(L)
          ENDIF
          DIST=DIST+SQRT((X(L)-XLAST)**2+(Y(L)-YLAST)**2)
          write(26,2000) j,real(dist/thou),real(ht(l)),
     &                real(bed(l)),real(adot(l)),
     &                real(flowa(l))
2000  format(1x,i6,5g14.6)
          XLAST=X(L)
          YLAST=Y(L)
C         DIST=X(L)
C         DIST=Y(L)
C         WRITE(23,*) DIST/THOU,HT(L),BED(L),PSURF(L)
C         WRITE(23,*) DIST/THOU,HT(L),BED(L),THICK(L)
          WRITE(2,200) DIST/THOU,HT(L),j
c         WRITE(*,200) DIST/THOU,HT(L),j
          WRITE(9,200) DIST/THOU,ADOT(L)*1.
          WRITE(7,200) DIST/THOU,PSURF(L)
          XMAX=MAX(DIST,XMAX)
          XMIN=MIN(DIST,XMIN)
          WRITE(3,200) DIST/THOU,BED(L)
          WRITE(4,200) DIST/THOU,PSURF(L)
C         WRITE(23,123) DIST,HT(L),BED(L),1.,10,100
          WRITE(23,124) 1,J,NINT(100*FRACT(L)),BED(L),1.,ADOT(L)*100,
     &                  HT(L),SQRT(X(L)**2+Y(L)**2)
C         WRITE(23,125) DIST,NINT(100*FRACT(L)),BED(L),1.,ADOT(L)*100,
C    &                  HT(L),SQRT(X(L)**2+Y(L)**2)
          YMAX=MAX(HT(L),YMAX)
          YMIN=MIN(HT(L),YMIN)
          YMAX=MAX(BED(L),YMAX)
          YMIN=MIN(BED(L),YMIN)
650     CONTINUE
        WRITE(2,200) RNINE,RMINUS,ICOLOR
C       WRITE(2,200) RNINE,RMINUS,0
        WRITE(7,200) RNINE,RMINUS,ICOLOR
        WRITE(9,200) RNINE,RMINUS,ICOLOR
        WRITE(3,200) RNINE,RMINUS,ICOLOR
        ICOLOR=ICOLOR+1
        IF(ICOLOR.GT.8) ICOLOR=1
        WRITE(2,211) HED
        WRITE(7,*) 'Psurf'
        WRITE(9,*) 'Adot'
        WRITE(3,*) 'BED'
C       WRITE(2,210) I
C       WRITE(3,210) I
210   FORMAT('LINE - ',I3)
700   CONTINUE
      close(2);close(3)
      GOTO 10
999   CONTINUE
      WRITE(24,200) RNINE,RMINUS,ICOLOR
      write(24,*) 'dome'
      WRITE(25,200) RNINE,RMINUS,ICOLOR
      write(25,*) 'volume'
9999  CONTINUE
c     if(iflush) call gflush
c     print *,'again??? input 1'
c     read(*,*) iagain
c     call setvis(iseg,.false.)
c     if(iagain.eq.1) then
c       call setvis(1,.true.)
c       do i=2,iseg
c         call setvis(i-1,.false.)
c         call setvis(i,.true.)
c         call wait
c       enddo
c       goto 9999
c     endif
c     do i=1,iseg
c       call setvis(-1,.true.)
c       call wait
c     enddo
c
c     WRITE(4,200) XMIN/THOU,YMIN
c     WRITE(4,200) XMAX/THOU,YMIN
c     WRITE(4,200) XMAX/THOU,YMAX
c     WRITE(4,200) XMIN/THOU,YMAX
c     WRITE(4,200) XMIN/THOU,YMIN
      WRITE(4,200) RNINE,RMINUS,-1
      WRITE(4,210)
124   FORMAT(3I5,5(1X,F10.2))
125   FORMAT(F10.0,I5,5(1X,F10.2))
123   FORMAT(F11.1,F7.0,F7.0,F7.3,I5,I5)
200   FORMAT(10X,G13.6,2X,G13.6,I13)
211   FORMAT(A80)
c     if(iflush) call gflush
c     call grstop1
      STOP
      END
      subroutine wait
      do i=1,10000
        x=sin(1.)
      enddo
      end
