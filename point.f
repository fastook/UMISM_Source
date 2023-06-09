C READS FROM STANDARD OUT3 AND PLOTS CONTOURS DIRECTLY
      PARAMETER(NMAX=79999,NCOLOR=15,NPLOT=10000)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL FOUND
      CHARACTER HED*80,ANSW,JUNK*80
      DIMENSION XOUT(NPLOT),YOUT(NPLOT),IOUT(NPLOT)
      DIMENSION X(44,40),Y(44,40),Z(44,40)
      DIMENSION XX(NMAX),YY(NMAX),ZZ(NMAX),NNODE(NMAX),NBOUND(400)
      DIMENSION PSURF(NMAX),FRACT(NMAX),FLOWA(NMAX),HTICE(NMAX)
      DIMENSION ADOT(NMAX)
      DIMENSION DEPB(NMAX),SLDGB(NMAX),AFUDGE(NMAX),TBED(NMAX)
      DIMENSION WTHICK(NMAX),BMELT(NMAX)
      DIMENSION XLINE(400),YLINE(400)
      DIMENSION KX(NMAX,4),IK(5)
      DIMENSION ICMAP(NCOLOR)
      DIMENSION INARAY(3)
      DATA ICMAP /1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/
      print *,'input node number of point to plot'
      read(*,*) mplot
      print *,'input base time'
      read(*,*) tbase
C READ INPUT HEADER
      READ(30,1000,END=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
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
      IREAD=IREAD+1
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME
      READ(33,END=999) HED
      write(junk,*) hed(8:)
      read(junk,*) time
      print *, hed(:40)
      print *,time,time+tbase
      READ(33) (HTICE(I),I=1,NUMNP)
      READ(33) (ADOT(I),I=1,NUMNP)
      READ(33) (DEPB(I),I=1,NUMNP)
      READ(33) (CONST,I=1,NUMEL)
      READ(33) (ACON,I=1,NUMEL)
      READ(34) (FRACT(I),I=1,NUMNP)
      READ(34) (FLOWA(I),I=1,NUMNP)
      READ(34) (SLDGB(I),I=1,NUMNP)
      READ(34) (AFUDGE(I),I=1,NUMNP)
      READ(36,END=999) HED
      READ(36) (TBED(I),I=1,NUMNP)
      READ(36) (BMELT(I),I=1,NUMNP)
      READ(36) (WTHICK(I),I=1,NUMNP)
      write(21,*) tbase+time,htice(mplot)-depb(mplot)
      write(22,*) tbase+time,depb(mplot)
      write(23,*) tbase+time,htice(mplot)
      if(iread.gt.1) then
        write(24,*) tbase+time,1000*(depb(mplot)-depblast)/(time-tlast)
      endif
      depblast=depb(mplot)
      tlast=time
      GOTO 1
999   CONTINUE
      write(21,*) -99999,0
      write(21,2000) 'load:',mplot
      write(22,*) -99999,0
      write(22,2000) 'bed :',mplot
      write(23,*) -99999,0
      write(23,2000) 'surf :',mplot
      write(24,*) -99999,0
      write(24,2000) 'rate :',mplot
2000  format(1x,a,i6)
      END
