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
      logical iflush
      common /flush/ iflush
      data iflush /.true./
      DATA ICMAP /1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/
      jcount=0
      print *,'input start number,'
      print *,'      display interval, and' 
      print *,'      how many to show'
      read(*,*) jstart,jdisp,jnumber
      if(jstart.lt.0) then
        read(99,*) jstart,jdisp,jnumber
      else
        write(99,*) jstart,jdisp,jnumber
      endif
      IPR=1
      ISTEP=1
      print *,'input base number,step:',ipr,istep
      read(*,*) ipr,istep
      IDISP=1
      II=1
      DO I=1,NPLOT
        READ(12,*,END=21) XIN,YIN,IIN
        IF(XIN.EQ.-99999.) THEN
          READ(12,22) JUNK
22    FORMAT(A80)
        ELSE
          XOUT(II)=XIN
          YOUT(II)=YIN
          IOUT(II)=IIN
          NOUT=II
          II=II+1
        ENDIF
      ENDDO
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
      READ(*,*,END=999) IPLOT
      IF(IPLOT.EQ.99) GOTO 999
      IWBASE=0
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
      ELSEIF(IPLOT.EQ.1) THEN
        RMAX=6000.
        RMIN=-3000.0001
      ELSEIF(IPLOT.EQ.2) THEN
        RMAX=4.
        RMIN=-4.
        RMIN=0.
C       RMIN=-.5
C       RMAX=.5
c        rmax=.2
C       RMIN=0.
C       RMAX=0.
      ELSEIF(IPLOT.EQ.3) THEN
        RMAX=1.
        RMIN=.01
c       RMAX=.01
      ELSEIF(IPLOT.EQ.4) THEN
        RMAX=7.
        RMIN=.0
      ELSEIF(IPLOT.EQ.5) THEN
        RMAX=3000.
        RMIN=-650.
      ELSEIF(IPLOT.EQ.6) THEN
        RMAX=6000.
        RMIN=0.001
      ELSEIF(IPLOT.EQ.7) THEN
        RMAX=1.5
        RMIN=1.0000
      ELSEIF(IPLOT.EQ.8) THEN
        RMAX=5000.
        RMIN=.0001
      ELSEIF(IPLOT.EQ.9) THEN
        RMAX=-110.
        RMIN=-120.
      ELSEIF(IPLOT.EQ.10) THEN
        RMAX=.05
        RMIN=0.
      ELSEIF(IPLOT.EQ.11) THEN
        RMAX=3.
        RMIN=0.0000
      ELSEIF(IPLOT.EQ.12) THEN
        RMAX=0.
        RMIN=-50.0001
      ELSEIF(IPLOT.EQ.13) THEN
        RMAX=10.0
        RMIN=-10.000001
      ELSEIF(IPLOT.EQ.14) THEN
        RMAX=10.001
        RMIN=0.001
      ELSEIF(IPLOT.EQ.15) THEN
        RMAX=-3000.00
        RMIN=150.00
      ELSEIF(IPLOT.EQ.16) THEN
        RMAX=-127.5
        RMIN=127.5
      ELSEIF(IPLOT.EQ.17) THEN
        RMAX=-1000.00
        RMIN=100.01
      ELSEIF(IPLOT.EQ.18) THEN
        RMAX=3000.
        RMIN=-650.
      ENDIF
      NLINE=0
      READ(11,*,END=99) NLINE
      IF(NLINE.GT.0) READ(11,1001) (NBOUND(I),I=1,NLINE)
99    CONTINUE
      ISTART=1
1001  FORMAT(16I6)
      RLEVSP=25.
      IF(IPLOT.EQ.0) RLEVSP=500.
      IF(IPLOT.EQ.1) RLEVSP=500.
      IF(IPLOT.EQ.2) RLEVSP=.05
      IF(IPLOT.EQ.2) RLEVSP=.25
      IF(IPLOT.EQ.3) RLEVSP=.15
      IF(IPLOT.EQ.4) RLEVSP=.50
      IF(IPLOT.EQ.5) RLEVSP=100.
      IF(IPLOT.EQ.6) RLEVSP=500.
      IF(IPLOT.EQ.7) RLEVSP=1.
      IF(IPLOT.EQ.8) RLEVSP=500.
      IF(IPLOT.EQ.9) RLEVSP=1.
      IF(IPLOT.EQ.10) RLEVSP=.010
      IF(IPLOT.EQ.11) RLEVSP=.5
      IF(IPLOT.EQ.12) RLEVSP=5.
      IF(IPLOT.EQ.13) RLEVSP=1.
      IF(IPLOT.EQ.14) RLEVSP=1.
      IF(IPLOT.EQ.15) RLEVSP=-100.
      IF(IPLOT.EQ.16) RLEVSP=-5.
      IF(IPLOT.EQ.17) RLEVSP=-50.
      IF(IPLOT.EQ.18) RLEVSP=100.
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
      DO 9 NUM=1,NUMNP
        XX(NUM)=XX(NUM)*.001
        YY(NUM)=YY(NUM)*.001
C IROT=0 FOR ORDINARY, IROT=1 FOR ROTATED 90 DEG
        IROT=0
        IF(IROT.EQ.1) THEN
          TEMP=XX(NUM)
          XX(NUM)=-YY(NUM)
          YY(NUM)=TEMP
        ENDIF
9     CONTINUE
      IREAD=0
c bypass the first few ...
      do l=1,jstart-1
        IF(IPLOT.EQ.0 .OR. IPLOT.EQ.1 .OR.IPLOT.EQ.2 .OR.
     &    IPLOT.EQ.5 .OR. IPLOT.EQ.6 .OR.IPLOT.EQ.7 .OR.
     &    IPLOT.EQ.9 .OR. IPLOT.EQ.3 .OR. IPLOT.EQ.4 .OR. 
     &    IPLOT.EQ.11) THEN
          READ(33,END=999) HED
          READ(33) (HTICE(I),I=1,NUMNP)
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
          IF(IWBASE.EQ.1) THEN
            IWBASE=IWBASE+1
            DO I=1,NUMNP
              WWWBASE(I)=WWW(I)
            ENDDO
          ELSE
            IWBASE=IWBASE+1
          ENDIF
        ENDIF
        print *,'bypassing: ',hed(1:40)
      enddo
1     CONTINUE
      IREAD=IREAD+1
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
        IF(IWBASE.EQ.1) THEN
          IWBASE=IWBASE+1
          DO I=1,NUMNP
            WWWBASE(I)=WWW(I)
          ENDDO
        ELSE
          IWBASE=IWBASE+1
        ENDIF
      ENDIF
      PRINT *,IREAD,IDISP,MOD(IREAD,IDISP)+1
      IF(MOD(IREAD,IDISP)+1.NE.1) GOTO 1
      jcount=jcount+1
      if(jcount.gt.jnumber*jdisp) then
        call grstop1
        stop
      endif
      if(1+mod(jcount,jdisp).ne.1) goto 1
        print *,'outputting for display',hed
      DO NUM=1,NUMNP
C CHANGE PLOTTED VARIABBLE HERE
        IF(IPLOT.EQ.0) ZZ(NUM)=HTICE(NUM)
        IF(IPLOT.EQ.1) ZZ(NUM)=DEPB(NUM)
        IF(IPLOT.EQ.2) ZZ(NUM)=ADOT(NUM)
        IF(IPLOT.EQ.3) ZZ(NUM)=FRACT(NUM)
        IF(IPLOT.EQ.4) ZZ(NUM)=FLOWA(NUM)
        IF(IPLOT.EQ.5) ZZ(NUM)=HTICE(NUM)-PSURF(NUM)
        IF(IPLOT.EQ.6 .OR. IPLOT.EQ.7) THEN
          ZZ(NUM)=HTICE(NUM)-DEPB(NUM)
          IF(DEPB(NUM).LT.0.) THEN
c            FLOT=-DEPB(NUM)*1.03/.917
            FLOT=DEPB(NUM)*(1.-1.03/.917)
            IF(HTICE(NUM).LT.FLOT) ZZ(NUM)=0.
          ENDIF
          IF(ZZ(NUM).LT.10.) ZZ(NUM)=0.
        ENDIF
        IF(IPLOT.EQ.8) ZZ(NUM)=PSURF(NUM)
        IF(IPLOT.EQ.9) ZZ(NUM)=ADOT(NUM)
        IF(IPLOT.EQ.10) ZZ(NUM)=SLDGB(NUM)
        IF(IPLOT.EQ.11) ZZ(NUM)=AFUDGE(NUM)
        IF(IPLOT.EQ.12) ZZ(NUM)=TBED(NUM)
        IF(IPLOT.EQ.13) ZZ(NUM)=BMELT(NUM)*1000.
        IF(IPLOT.EQ.14) ZZ(NUM)=WTHICK(NUM)
        IF(IPLOT.EQ.15) ZZ(NUM)=WWW(NUM)
        IF(IPLOT.EQ.16) ZZ(NUM)=WRATE(NUM)
        IF(IPLOT.EQ.17) ZZ(NUM)=WWW(NUM)-WWWBASE(NUM)
        IF(IPLOT.EQ.18) ZZ(NUM)=HTICE(NUM)-BSURF(NUM)
      ENDDO
200   FORMAT(I6,I4,2E12.5,F10.2,F7.2,F9.2,F10.3,F10.1,F10.2,F10.3)
      DO I=1,NUMEL
        NNODE(I)=4
        IF(KX(I,4).EQ.0) NNODE(I)=3
      ENDDO
C     RMAX=-1.D+70
      IF(ISTART.EQ.1) THEN
C       ISEG=1
        XMAX=-1.E30
        YMAX=XMAX
        XMIN=-XMAX
        YMIN=XMIN
        IF(NLINE.GT.0) THEN
          DO I=1,NLINE
            XMAX=MAX(XX(NBOUND(I)),XMAX)
            XMIN=MIN(XX(NBOUND(I)),XMIN)
            YMAX=MAX(YY(NBOUND(I)),YMAX)
            YMIN=MIN(YY(NBOUND(I)),YMIN)
          ENDDO
        ELSE
          DO I=1,NUMNP
            XMAX=MAX(XX(I),XMAX)
            XMIN=MIN(XX(I),XMIN)
            YMAX=MAX(YY(I),YMAX)
            YMIN=MIN(YY(I),YMIN)
          ENDDO
        ENDIF
        XDELT=(XMAX-XMIN)*.2
        YDELT=(YMAX-YMIN)*.2
        XMAX=XMAX+XDELT
        YMAX=YMAX+YDELT
        XMIN=XMIN-XDELT
        YMIN=YMIN-YDELT
        IF(ISEG.EQ.1) CALL GRSTRT(600,600)
        IF((XMAX-XMIN).GT.(YMAX-YMIN)) THEN
          YMAX=YMIN+(XMAX-XMIN)*1.
      PRINT *,1.
        ELSE
          XMAX=XMIN+(YMAX-YMIN)*1.
      PRINT *,1.
        ENDIF
        CALL WINDOW(REAL(XMIN),REAL(XMAX),REAL(YMIN),REAL(YMAX))
        WRITE(7,*) XMIN,XMAX,YMIN,YMAX
        ISTART=2
      ENDIF
      IF(ISEG.EQ.1) THEN
        CALL OPNSEG(ISEG)
        CALL MOVE(REAL(XMIN),REAL(YMIN))
        CALL DRAW(REAL(XMIN),REAL(YMAX))
        CALL DRAW(REAL(XMAX),REAL(YMAX))
        CALL DRAW(REAL(XMAX),REAL(YMIN))
        CALL DRAW(REAL(XMIN),REAL(YMIN))
        CALL LINCLR(1)
        CALL DASHPT(3)
        DO I=1,NOUT
          IF(IOUT(I).EQ.1) THEN
            CALL MOVE(REAL(XOUT(I)),REAL(YOUT(I)))
          ELSE
            CALL DRAW(REAL(XOUT(I)),REAL(YOUT(I)))
          ENDIF
        ENDDO
        CALL CLOSEG
        CALL DASHPT(0)
        ISEG=ISEG+1
      ENDIF
      CALL OPNSEG(ISEG)
      CALL MOVE(REAL(XMIN+(XMAX-XMIN)*.05),REAL(YMAX-(YMAX-YMIN)*.05))
      CALL TXTCLR(1)
      CALL TEXT(80,HED)
      CALL LINCLR(8)
      if(nline.gt.0) then
        CALL MOVE(REAL(XX(NBOUND(1))),REAL(YY(NBOUND(1))))
        DO I=2,NLINE
          CALL DRAW(REAL(XX(NBOUND(I))),REAL(YY(NBOUND(I))))
        ENDDO
      endif
C     RMIN=1.D+70
C     DO I=1,NUMNP
C       RMAX=MAX(ZZ(I),RMAX)
C       RMIN=MIN(ZZ(I),RMIN)
C     ENDDO
      LEVEL=(RMAX-RMIN)/RLEVSP+1
      DLEV=RLEVSP
      IF(LEVEL.LE.0) THEN
        LEVEL=2-LEVEL
        DLEV=-DLEV
      ENDIF
      ICOLOR=1
C     IF(IPLOT.EQ.7) ICOLOR=1+MOD(ISEG,NCOLOR)
      CALL LINCLR(ICMAP(ICOLOR))
      IPASS=1
      XPOS=XMAX-(XMAX-XMIN)*.125
      YPOS=YMAX-(YMAX-YMIN)*.1
      DO 500 LEV=1,LEVEL
        FOUND=.FALSE.
        ICOUNT=0
        VAL=RMIN+(LEV-1)*DLEV
        DO 450 N=1,NUMEL
          IPASS=1
          DO 415 J=1,NNODE(N)
            IK(J)=KX(N,J)
415       CONTINUE
          IK(NNODE(N)+1)=KX(N,1)
          ICOUNT=0
          DO 425 NN=1,NNODE(N)
            I=IK(NN)
            J=IK(NN+1)
            IF(ZZ(I).LT.VAL .AND. VAL.LE.ZZ(J)) GOTO 430
            IF(ZZ(J).LT.VAL .AND. VAL.LE.ZZ(I)) GOTO 430
            GOTO 425
430         CONTINUE
            FOUND=.TRUE.
            ICOUNT=ICOUNT+1
            XINT=XX(J)+(VAL-ZZ(J))*(XX(I)-XX(J))/(ZZ(I)-ZZ(J))
            YINT=YY(J)+(VAL-ZZ(J))*(YY(I)-YY(J))/(ZZ(I)-ZZ(J))
            IF(IPASS.EQ.1) THEN
              CALL MOVE(REAL(XINT),REAL(YINT))
              IPASS=2
            ELSE
              CALL DRAW(REAL(XINT),REAL(YINT))
            ENDIF
425       CONTINUE
          IF(ICOUNT.EQ.0) GOTO 450
450     CONTINUE
        IF(FOUND) THEN
          CALL MOVE(REAL(XPOS),REAL(YPOS))
          CALL TXTCLR(ICMAP(ICOLOR))
          CALL RNUMBR(REAL(VAL),3,9)
          YPOS=YPOS-(YMAX-YMIN)*.025
        ENDIF
        ICOLOR=ICOLOR+1
        IF(ICOLOR.GT.NCOLOR) ICOLOR=1
        CALL LINCLR(ICMAP(ICOLOR))
500   CONTINUE
98    FORMAT(3G13.4)
97    FORMAT('=',F10.3)
      IF(IPLOT.NE.7) ICOLOR=1
      CALL CLOSEG
      IF(ISEG.NE.2) CALL SETVIS(ISEG-1,.FALSE.)
C     CALL NEWPAG
      ISEG=ISEG+1
C     CALL GETUTX(19,'INPUT C TO CONTINUE',1,ANSW,IGOT)
      CALL savescreen(IPR,600,600)
      IPR=IPR+ISTEP
      GOTO 1
999   CONTINUE
      CALL GRSTOP1
      STOP
      END
