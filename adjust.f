      SUBROUTINE ADJUST(HED, NUMNP, NUMEL, XX, YY, HTICE, ADOT, ADOTB,
     &              FRACT,BMELT,WTHICK,
     &              TEMP, ITYPE, TBED, PSURF,UNDEPB,
     &              BDROCK, DEPB, FLOWA, SLDGB, THICK, KX, CONST, 
     &              AFUDGE,NNODE, KODE, GEOFLUX, FLUX,
     &              HFIT, NUMCOL, NUMLEV, NUMGBC, NDT, INTER, DT,
     &              IBFLUX, BFLUX, MXX, IDT, SLOPN, AMASS, TIME,
     &              NTSTEP, TTIME, VVOL, AAREA,TTBOT,TTAVG, TTNSL,
     &              TTSEAL, IFIT,
     &              IPLOT, CFACTOR, ACON,ICON, 
     &              ALPHAC,TBASE,
     &              IMELT,TWATER,PWATER,CALV,
     &              NTYPE,AADOT,AFRACT,ABDRCK,PPSURF,
     &              AFLOWA,ASLDGB,PCALV,ADC,WWWMIN,WWW,WRATE,WDIFF,
     &              WWWORIG,ADVANCE,HIRES,ACC,ABLAT,TFRACT)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( MMXX=MAXNUM, NSMAX=MAXTIME)
      EXTERNAL IFIND
      LOGICAL HIRES
      LOGICAL LFOUND
      REAL*4  PX(4),PY(4)
      COMMON /SNOW/ SNOLIN,SNO(20)
      COMMON /LINE/ NP,NLINE(1000)
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      COMMON /EXPER/ TPERIOD,TINIT,TVSTART,TVFINAL,IEXPER
      COMMON /VELOS/ ASCAL,UMAX,VTHRESH,INORM
      COMMON /TCONST/ TSORIG(MMXX),TFLAG
      LOGICAL TFLAG
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      DIMENSION ALPHAC(3)
      DIMENSION NNLINE(1000),DLINE(1000),IFIT(7),TMP(MMXX)
      DIMENSION TTIME(NSMAX),VVOL(NSMAX),AAREA(NSMAX),TTBOT(NSMAX),
     &          TTNSL(NSMAX),TWATER(NSMAX),PWATER(NSMAX),TTAVG(NSMAX),
     &          TTSEAL(NSMAX),
     &          WWWMIN(NSMAX),WWW(3*MXX),WRATE(3*MXX,2),WDIFF(MXX)
      DIMENSION AMASS(11), KODE(MXX), XX(MXX), YY(MXX), HTICE(MXX),
     &          SLOPE(4,MMXX), SLOPN(4,MXX), ADOT(MXX), FRACT(MXX), 
     &          PSURF(MXX), BDROCK(MXX), DEPB(MXX), FLOWA(MXX), 
     &          SLDGB(MXX),UNDEPB(MXX),WWWORIG(MXX),
     &          TEMP(MXX),ADOTB(MXX),TBED(MXX),AFUDGE(MXX),
     &          GEOFLUX(MXX),FLUX(MXX),
     &          KX(MXX,4), CONST(MXX), THICK(MXX), NNODE(MXX),
     &          ZZ(MMXX), ICMAP(16), XCHECK(5), YCHECK(5), IFD(MMXX),
     &          IFP(MMXX), JFP(MMXX), HFIT(MXX), ITYPE(MXX),
     &          IBFLUX(MXX,2), BFLUX(MXX), IDT(MXX),
     &          XA(5), YA(5), IDAT(5),ACON(MXX),VEL(MMXX,3),
     &          BMELT(MMXX),WTHICK(MMXX),CALV(MMXX)
      DIMENSION NTYPE(MMXX), AADOT(MMXX), AFRACT(MMXX), ABDRCK(MMXX),
     &          PPSURF(MMXX), AFLOWA(MMXX), ASLDGB(MMXX), PCALV(MMXX),
     &          ADC(MMXX)
      DIMENSION ACC(MMXX),ABLAT(MMXX)
      REAL*4 XA,YA
c     REAL*8 RAND
      CHARACTER CHAR*1, CHAR2*1,HED*80, JUNK*80, CHAR3*1, CHAR4*1
      DIMENSION ICMAP1(16)
      LOGICAL LDANGLE
      common /flush/ iflush
      logical iflush
      common /heat/ htogg
      logical htogg
      logical zflag
      common /zoomcom/ pxmin,pymin,pxmax,pymax,zflag
      LOGICAL BATCH
      COMMON /BATCHER/ BATCH
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      logical file20,file10,file46
      SAVE ISEED,THRESH,bedthresh,bedthresh2
      DATA ICMAP /0,5,11,4,12,6,13,2,8,7,9,3,10,14,15,1/
      DATA ICMAP1 /0,5,11,4,12,6,13,2,8,7,9,3,10,14,15,1/
C      DATA ICMAP1 /0,4,12,6,13,2,8,7,9,3,10,5,11,1,14,15/                
      DATA ISEED /23/,THRESH /0.0D0/
      data bedthresh /0.d0/, bedthresh2 /0.d0/
c      zflag=.false.
      JUNK=' '
      IGRAP=0
      NP=0
      inquire(file='fort.20',exist=file20)
      if(file20) then
        REWIND 20
        READ(20,1010,END=102) NP
        READ(20,1010) (NLINE(IP),IP=1,NP)
      endif
102   continue
      RHOW=1.092d0
      RHOI=.917d0
      FLOWMIN=1d-6
c      FLOWMIN=2.29d0
      SLDGBMIN=1d-9
      IPASS=0
      ICOUNT=NUMNP
      DO I=1,NUMNP
        IFP(I)=I
        ZZ(I)=HTICE(I)-HFIT(I)
      ENDDO
      XLEN=100.d0
      YLEN=100.d0
      ZLEN=14.d0
      ILINE=0
      GOTO 2100
C
C BEGIN MAIN PLOT
C
2000  CONTINUE
      WRITE(*,*) 'GRAPHICS: '
      WRITE(*,*) '<M>AP OR <P>ROFILE OR <C>ONTOUR OR <3>-D'
      WRITE(*,*) '<S>CALED CONTOUR OR <L>OG CONTOUR '
      WRITE(*,*) '<U>SER DEFINED '
      if(iflush) call gflush
      READ(*,1002) CHAR2
      WRITE(99,1002) CHAR2
      IF(CHAR2.EQ.'m') THEN
        CHAR2='M'
      ELSEIF(CHAR2.EQ.'p') THEN
        CHAR2='P'
      ELSEIF(CHAR2.EQ.'c') THEN
        CHAR2='C'
      ELSEIF(CHAR2.EQ.'s') THEN
        CHAR2='S'
      ELSEIF(CHAR2.EQ.'l') THEN
        CHAR2='L'
      ELSEIF(CHAR2.EQ.'u') THEN
        CHAR2='U'
      ENDIF
      !IF(CHAR2.EQ.'3' .AND. .NOT.BATCH) THEN
      IF(CHAR2.EQ.'3') THEN
        CALL THREED(1,MMXX,NUMNP,NUMEL,XX,YY,ZZ,ZZ,KX,BATCH)
        GOTO 2100
      ELSEIF(CHAR2.EQ.'P') THEN
        IF(NP.LT.1) THEN
          PRINT *,'MUST DEFINE A LINE FIRST: <I>'
          GOTO 2100
        ENDIF
        IF(IPASS.EQ.0) THEN
          CALL GRSTRT(600,600)
          IPASS=1
        ENDIF
        ZMIN=1d30
        ZMAX=-ZMIN
        DMIN=1d30
        DMAX=-ZMIN
        XBASE=XX(NLINE(1))
        YBASE=YY(NLINE(1))
c        DO I=1,NP
c          DLINE(I)=SQRT((XX(NLINE(I))-XBASE)**2+
c     &                 (YY(NLINE(I))-YBASE)**2)
        DLINE(1)=0.0d0
        DMIN=MIN(DMIN,DLINE(1))
        DMAX=MAX(DMAX,DLINE(1))
        ZMIN=MIN(ZMIN,ZZ(NLINE(1)))
        ZMAX=MAX(ZMAX,ZZ(NLINE(1)))
        print *,1,zz(nline(1))
        DO I=2,NP
          DLINE(I)=DLINE(I-1)+
     &            SQRT((XX(NLINE(I))-XX(NLINE(I-1)))**2+
     &                 (YY(NLINE(I))-YY(NLINE(I-1)))**2)
          DMIN=MIN(DMIN,DLINE(I))
          DMAX=MAX(DMAX,DLINE(I))
          ZMIN=MIN(ZMIN,ZZ(NLINE(I)))
          ZMAX=MAX(ZMAX,ZZ(NLINE(I)))
          print *,i,zz(nline(i))
        ENDDO
        IF(IGRAP.EQ.0) THEN
c          CALL NEWPAG
          IGRAP=1
        ENDIF
        IF(CHAR.EQ.'U' .OR. CHAR.EQ.'N' .OR.
     &     CHAR.EQ.'B' .OR. 
     &     CHAR.EQ.'H' .OR. 
     &     CHAR.EQ.'O') THEN
          CALL WINDOW(REAL(DMIN),REAL(DMAX),-1000.,19000.)
          PRINT *,'ZMIN,ZMAX=',-1000.,19000.
c         CALL WINDOW(REAL(DMIN),REAL(DMAX),-6000.,0.)
c         PRINT *,'ZMIN,ZMAX=',-6000.,0.
c         CALL WINDOW(REAL(DMIN),REAL(DMAX),real(zmin),real(zmax))
c         PRINT *,'ZMIN,ZMAX=',zmin,zmax
        elseif(char.eq.'M') then
          CALL WINDOW(REAL(DMIN),REAL(DMAX),-150.,10.)
          PRINT *,'ZMIN,ZMAX=',-150.,10.
        ELSE
          PRINT *,'ZMIN,ZMAX=',ZMIN,ZMAX
          ZPLUS=(ZMAX-ZMIN)/25.d0
          ZMIN=ZMIN-ZPLUS
          ZMAX=ZMAX+ZPLUS
          CALL WINDOW(REAL(DMIN),REAL(DMAX),REAL(ZMIN),REAL(ZMAX))
        ENDIF
        CALL LINCLR(1)
        CALL MOVE(REAL(DLINE(1)),0.)
        CALL DRAW(REAL(DLINE(NP)),0.)
        CALL MOVE(REAL(DLINE(1)),REAL(ZZ(NLINE(1))))
        DO I=2,NP
          CALL DRAW(REAL(DLINE(I)),REAL(ZZ(NLINE(I))))
        ENDDO
c       CALL WINDOW(REAL(DMIN),REAL(DMAX),-1000.,19000.)
        CALL WINDOW(REAL(DMIN),REAL(DMAX),-6000.,0.)
c       CALL WINDOW(REAL(DMIN),REAL(DMAX),real(zmin),real(zmax))
        CALL LINCLR(3)
        CALL MOVE(REAL(DLINE(1)),REAL(DEPB(NLINE(1))))
        DO I=2,NP
          CALL DRAW(REAL(DLINE(I)),REAL(DEPB(NLINE(I))))
        ENDDO
        CALL LINCLR(2)
        CALL MOVE(REAL(DLINE(1)),REAL(PSURF(NLINE(1))))
        DO I=2,NP
          CALL DRAW(REAL(DLINE(I)),REAL(PSURF(NLINE(I))))
        ENDDO
        GOTO 2100
      ELSEIF(CHAR2.EQ.'C' .OR. CHAR2.EQ.'U' .OR.
     &       CHAR2.EQ.'S' .OR.
     &       CHAR2.EQ.'L') THEN
        PRINT *,' CONTOURS'
        IF(IPASS.EQ.0) THEN
          CALL GRSTRT(600,600)
          IPASS=1
        ENDIF
        IF(IGRAP.EQ.0) THEN
c          CALL NEWPAG
          IGRAP=1
        ENDIF
        CALL SCALEXY(xx,yy,NUMNP,xmin,xmax,ymin,ymax,delx,dely)
C
        IF(CHAR.EQ.'L' .OR. CHAR.EQ.'O' .OR.
     &     CHAR.EQ.'Y' .OR. 
     &     CHAR2.EQ.'S' .OR. CHAR2.EQ.'U' .OR.
     &     CHAR2.EQ.'L') THEN
          ZMIN=1.d30
          ZMAX=-ZMIN
          DO I=1,NUMNP
            ZMAX=MAX(ZMAX,ZZ(I))
            ZMIN=MIN(ZMIN,ZZ(I))
c            if(zz(i).ne.0d0) print *,i,zz(i)
          ENDDO
      print *,'min.max=',zmin,zmax
          IF(CHAR2.EQ.'U') THEN
            PRINT *,'INPUT USER MIN,MAX'
            READ *,ZMIN,ZMAX
            WRITE(99,*) ZMIN,ZMAX
          ENDIF
          IF(ZMIN.EQ.ZMAX) ZMAX=ZMIN+1.d0
          DELZ=(ZMAX-ZMIN)/14.d0
          DDELZ=(ZMAX-ZMIN)*1d-9
          PRINT *,ZMIN,ZMAX,DDELZ
c          ZMIN=ZMIN+DELZ/100.d0
          ZMIN=ZMIN+DDELZ
          ZMAX=ZMAX-DELZ/100.d0
        ENDIF
        DELZ=(ZMAX-ZMIN)/14.d0
        IF(CHAR2.EQ.'L' .AND. ZMIN.GT.0.) THEN
          PRINT *,' LOG-PLOT...',ZMIN,ZMAX
          RLOG2=1.d0/LOG(10.d0)
          ZMIN=LOG(ZMIN+0.1)*RLOG2
          ZMAX=LOG(ZMAX+0.1)*RLOG2
          DELZ=(ZMAX-ZMIN)/14.d0
          DO I=1,NUMNP
            IF(ZZ(I).GT.0) THEN
              ZZ(I)=LOG(ZZ(I)+0.1)*RLOG2
            ELSE
              ZZ(I)=-18.d0
            ENDIF
          ENDDO
        ENDIF
        CALL WINDOW(REAL(XMIN-DELX),REAL(XMAX+DELX),
     &              REAL(YMIN-DELY),REAL(YMAX+DELY))
        CALL CONTR(MMXX,NUMEL,XX,YY,KX,
     &                 ZZ,ZMIN,ZMAX,DELZ,
     &                 XMIN,XMAX,YMIN,YMAX)
        WRITE(*,*) ASCAL,UMAX,VTHRESH 
c        IF(CHAR.EQ.'G') THEN
c          CALL DGRAD(MMXX,NUMEL,KX,XX,YY,HTICE,DEPB,.TRUE.)
c        ENDIF
        IF(CHAR4.EQ.'W') THEN
          CALL WVELO(XX,YY, KX, NUMNP, NUMEL,
     &               WTHICK, HTICE, DEPB, .false.)
        ELSE
          CALL DVELO(MMXX,NUMEL,KX,XX,YY,HTICE,DEPB,CONST,
     &               VEL,.TRUE.,DTMIN)
        ENDIF
        CALL DCOAST(1)
        CALL LINCLR(1)                                                    
        CALL MOVE(REAL(XMIN),REAL(YMIN))                                  
        CALL DRAW(REAL(XMIN),REAL(YMAX))                                  
        CALL DRAW(REAL(XMAX),REAL(YMAX))                                  
        CALL DRAW(REAL(XMAX),REAL(YMIN))                                  
        CALL DRAW(REAL(XMIN),REAL(YMIN))                                  
        GOTO 2100
      ENDIF
      IF(IPASS.EQ.0) THEN
        CALL GRSTRT(600,600)
        IPASS=1
c      ELSE
c        CALL NEWPAG
      ENDIF
      IGRAP=0
      CALL WINDOW(0.,100.,0.,100.)
      CALL SCALEXY(xx,yy,NUMNP,xmin,xmax,ymin,ymax,delx,dely)
c      CALL NEWPAG
      CALL WINDOW(REAL(XMIN-DELX),REAL(XMAX+DELX),
     &              REAL(YMIN-DELY),REAL(YMAX+DELY))
c     CALL OPNSEG(1)
      CALL SCALE3(XX,XLEN,NUMNP,1)
      CALL SCALE3(YY,YLEN,NUMNP,1)
      CALL SCALE2(ZZ,ZLEN,NUMNP,1)
C FOLLOWING TO MAKE X AND Y AXIS SAME SCALE
      IF(XX(NUMNP+2).GT.YY(NUMNP+2)) THEN
        YY(NUMNP+2)=XX(NUMNP+2)
      ELSE
        XX(NUMNP+2)=YY(NUMNP+2)
      ENDIF
C END OF SAME SCALE
      XD=XX(NUMNP+2)
      YD=YY(NUMNP+2)
      ZF=ZZ(NUMNP+1)
      ZD=ZZ(NUMNP+2)
      IF(CHAR.EQ.'I') THEN
        ZF=-1300.d0
        ZD=200.d0
c        ZF=-130.d0
c        ZD=20.d0
      ELSEIF(CHAR.EQ.'J') THEN
        ZF=-13.999d0
        ZD=2.d0
      ELSEIF(CHAR.EQ.'Z') THEN
        ZF=7.d0
        ZD=1.d0
      ELSEIF(CHAR.EQ.'A') THEN
        ZF=0.d0
        ZD=.25d0
      ELSEIF(CHAR.EQ.'F') THEN
        ZF=0.d0
        ZD=.5d0
      ELSEIF(CHAR.EQ.'H') THEN
        ZF=-1000.d0
        ZD=1000.d0
      ELSEIF(CHAR.EQ.'U') THEN
        ZF=-1000.d0
        ZD=1000.d0
      ELSEIF(CHAR.EQ.'B') THEN
        ZF=-2000.d0
        ZD=400.d0
C       ZD=100.d0
C       ZF=-1400.d0
      ELSEIF(CHAR.EQ.'D') THEN
        ZF=-.8d0
        ZD=.20d0
c      ELSEIF(CHAR.EQ.'F') THEN
c        ZF=.5d0
c        ZD=.4d0
      ELSEIF(CHAR.EQ.'T') THEN
        ZF=-50.d0
        ZD=100.d0
        ZF=-200.d0
        ZD=250.d0
c following for very thin ice stuff
c        zf=-5.d0
c        zd=10.d0
      ELSEIF(CHAR.EQ.'M') THEN
C        ZF=NINT(ZF)
          ZF=-56.d0
          ZD=4.d0
      ENDIF
      CALL TXTCLR(1)
      ILZ=int(ZLEN+1)
      PXX=80.d0
      PYY=92.5d0
c new ***********
      call linclr(1)
      call move(real(xmin),real(ymin))
      call draw(real(xmax),real(ymin))
      call draw(real(xmax),real(ymax))
      call draw(real(xmin),real(ymax))
      call draw(real(xmin),real(ymin))
      XPOS=XMAX+(XMAX-XMIN)*0.025d0
      YPOS=YMIN+(YMAX-YMIN)*0.90d0
      PXX=xpos
      PYY=ypos
      RNUM=ZF+ZLEN*ZD
      DO I=1,ILZ
        IZ=ILZ-I+2
        CALL MRKCLR(ICMAP(IZ))
        CALL MARKER(REAL(PXX),REAL(PYY),0)
c
        CALL MOVE(REAL(PXX+delx*0.10),REAL(PYY-dely*0.05))
        CALL LINCLR(1)
C       CALL RNUMBR(REAL(RNUM),3,9)
        CALL GNUMBR(REAL(RNUM),3,9)
        RNUM=RNUM-ZD
        YPOS=YPOS-(YMAX-YMIN)*.04d0
        PYY=YPOS
      ENDDO
      RNUM=ZF
      DO I=1,NUMNP
        RZ=(ZZ(I)-ZF)/ZD
        IZ=int(RZ+2)
        IF(IZ.LT.2) IZ=2
        IF(IZ.GT.16) IZ=16
        CALL MRKCLR(ICMAP(IZ))
        xxx=xx(i)*1.d-3
        yyy=yy(i)*1.d-3
        IF(KODE(I).EQ.0) THEN
C         CALL MARKER(XXX,YYY,9)
          CALL MARKER(REAL(XXX),REAL(YYY),0)
        ELSE
          CALL MARKER(REAL(XXX),REAL(YYY),1)
        ENDIF
      ENDDO
      IF(ILINE.NE.0) THEN
        IEC=2
        CALL LINCLR(IEC)
C       CALL DASHPT(1)
        DO I=1,NUMEL
          xxx=xx(KX(I,1))*1.d-3
          yyy=yy(KX(I,1))*1.d-3
          CALL MOVE(REAL(XXX),REAL(YYY))
          DO J=2,NNODE(I)
            xxx=xx(KX(I,J))*1.d-3
            yyy=yy(KX(I,J))*1.d-3
            CALL DRAW(REAL(XXX),REAL(YYY))
          ENDDO
          xxx=xx(KX(I,1))*1.d-3
          yyy=yy(KX(I,1))*1.d-3
          CALL DRAW(REAL(XXX),REAL(YYY))
          IEC=IEC+1
          IF(IEC.GT.14) IEC=2
          CALL LINCLR(IEC)
        ENDDO
        CALL MAKCUR
      ENDIF
      IF(CHAR.EQ.'D') THEN
        CALL MRKCLR(1)
        xxx=xpole
        yyy=ypole
      PRINT *,XXX,XPOLE
      PRINT *,YYY,YPOLE
        CALL MARKER(REAL(XXX),REAL(YYY),1)
      ENDIF
C
C READ AND PLOT OUTLINE IN FILE OUTLINE DATA B
      CALL LINCLR(15)
      IREAD=0
      inquire(file='fort.10',exist=file10)
      if(.not.file10) goto 124
123   READ(10,*,END=124) XOUT,YOUT,IREAD
      IF(XOUT.EQ.-99999.) THEN
        READ(10,1000) JUNK
1000  FORMAT(A80)
        GOTO 123
      ENDIF
      IF(XOUT.EQ.-99999.) GOTO 124
      IF(IREAD.EQ.1) THEN
        CALL MOVE(REAL(XOUT),REAL(YOUT))
      ELSE
        CALL DRAW(REAL(XOUT),REAL(YOUT))
      ENDIF
      IREAD=1
      GOTO 123
124   CONTINUE
      REWIND 10
c ... read the points.GMT file ............
      inquire(file='fort.46',exist=file46)
      if(.not.file46) goto 126
125   READ(46,*,END=126) XOUT,YOUT
        CALL POINT(REAL(XOUT),REAL(YOUT))
        goto 125
126   REWIND 46
c .........................................
C
c     CALL CLOSEG
      CALL MOVE(0.,0.)
      CALL RNUMBR(0.,-1,1)
      CALL MAKCUR
C
C     ******** MAIN MENU ********************************************
C
2100  CONTINUE
      IF(IPASS.EQ.1) CALL MAKCUR
      WRITE(*,*) 'INPUT:'
C
      WRITE(*,*) 'A: GR<A>PH VARIOUS DATA SET PROPERTIES'
      WRITE(*,*) 'B: <B>ACKUP THE DATA SET AS IT CURRENTLY STANDS'
      WRITE(*,*) 'C: <C>HANGE VARIOUS DATA SET PROPERTIES'
      WRITE(*,*) 'D: <D>ISPLAY ASCII DATA ON SELECTED NODES'
      WRITE(*,*) 'E: <E>LEMENT EDGES TOGGLE ON AND OFF'
      WRITE(*,*) 'F: <F>IT  PROPERTIES TO IMPROVE FIT AUTOMATICALLY'
      WRITE(*,*) 'G: <G>RAPH MASS BALANCE'
      IF(WTOGG.eq.0) THEN 
        WRITE(*,*) 'H: <H> TOGGLE BASAL WATER SOLVER: OFF',WTOGG
      elseif(WTOGG.eq.1) then
        WRITE(*,*) 'H: <H> TOGGLE BASAL WATER SOLVER: SIMPLE',WTOGG
      elseif(WTOGG.eq.2) then
        WRITE(*,*) 'H: <H> TOGGLE BASAL WATER SOLVER: JIM''S',WTOGG
      ELSEif(WTOGG.eq.3) then
        WRITE(*,*) 'H: <H> TOGGLE BASAL WATER SOLVER: JESSE''S',WTOGG
      else
        write(*,*) 'problems with WTOGG:',WTOGG
        pause
      ENDIF
      WRITE(*,*) 'I: L<I>NE OF NODES SELECTED WITH CURSOR AND <A> KEY'
      IF(IOTOGG) THEN 
        WRITE(*,*) 'J: <J> RUNNING I/O:TURN OFF'
      ELSE
        WRITE(*,*) 'J: <J> RUNNING I/O:TURN ON'
      ENDIF
      WRITE(*,*) 'K: <K>LIMATE KURVE '
      WRITE(*,*) 'L: NEW E<L>EMENT DEFINITION'
      WRITE(*,*) 'M: <M>OVE NODE'
      WRITE(*,*) 'N: NEW <N>ODE DEFINITION'
      WRITE(*,*) 'O: PL<O>T VARIOUS OUTPUT DURING RUN'
      IF(BTOGG.eq.0) THEN 
        WRITE(*,*) 'P: TOGGLE <P>LATE-BED SOLVER: OLD',BTOGG
      elseif(BTOGG.eq.1) then
        WRITE(*,*) 'P: TOGGLE <P>LATE-BED SOLVER: LOCAL',BTOGG
      elseif(BTOGG.eq.2) then
        WRITE(*,*) 'P: TOGGLE <P>LATE-BED SOLVER: ELASTIC PLATE',BTOGG
      ELSEif(BTOGG.eq.3) then
        WRITE(*,*) 'P: TOGGLE <P>LATE-BED SOLVER: VISCO-ELASTIC',BTOGG
      else
        write(*,*) 'problems with BTOGG:',BTOGG
        pause
      ENDIF
      WRITE(*,*) 'Q: <Q>UIT AND CONTINUE RUN'
      WRITE(*,*) 'R: <R>EMOVE ELEMENTS CONTAINING SELECTED NODES'
      WRITE(*,*) 'S: <S>ELECT NODES WITH COUNTER CLOCKWISE CORNERS'
      IF(ITOGG) THEN
        WRITE(*,*) 'T: <T>URN TNSL-ABLATION FOLLOWER ON',TFRACT
      ELSE
        WRITE(*,*) 'T: <T>URN TNSL-ABLATION FOLLOWER OFF',TFRACT
      ENDIF
      WRITE(*,*) 'U: <U>NZOOM ON GRAPH'
      WRITE(*,*) 'V: <V>OLUME VS TIME DISPLAY'
      IF(CTOGG) THEN 
        WRITE(*,*) 'W: <W> TOGGLE MASS BALANCE SOLVER:TURN OFF'
      ELSE
        WRITE(*,*) 'W: <W> TOGGLE MASS BALANCE SOLVER:TURN ON'
      ENDIF
      WRITE(*,*) 'X: DEFINE FLU<X> ALONG SELECTED BOUNDARY'
      WRITE(*,*) 'Y: MODIF<Y> DATA SET PROPERTIES TO IMPROVE FIT'
      WRITE(*,*) 'Z: <Z>OOM ON GRAPH'
      if(iflush) call gflush
      READ(*,1002) CHAR
      WRITE(99,1002) CHAR
      IF(CHAR.EQ.'a') THEN
        CHAR='A'
      ELSEIF(CHAR.EQ.'b') THEN
        CHAR='B'
      ELSEIF(CHAR.EQ.'c') THEN
        CHAR='C'
      ELSEIF(CHAR.EQ.'d') THEN
        CHAR='D'
      ELSEIF(CHAR.EQ.'e') THEN
        CHAR='E'
      ELSEIF(CHAR.EQ.'f') THEN
        CHAR='F'
      ELSEIF(CHAR.EQ.'g') THEN
        CHAR='G'
      ELSEIF(CHAR.EQ.'h') THEN
        CHAR='H'
      ELSEIF(CHAR.EQ.'i') THEN
        CHAR='I'
      ELSEIF(CHAR.EQ.'j') THEN
        CHAR='J'
      ELSEIF(CHAR.EQ.'k') THEN
        CHAR='K'
      ELSEIF(CHAR.EQ.'l') THEN
        CHAR='L'
      ELSEIF(CHAR.EQ.'m') THEN
        CHAR='M'
      ELSEIF(CHAR.EQ.'n') THEN
        CHAR='N'
      ELSEIF(CHAR.EQ.'o') THEN
        CHAR='O'
      ELSEIF(CHAR.EQ.'p') THEN
        CHAR='P'
      ELSEIF(CHAR.EQ.'q') THEN
        CHAR='Q'
      ELSEIF(CHAR.EQ.'r') THEN
        CHAR='R'
      ELSEIF(CHAR.EQ.'s') THEN
        CHAR='S'
      ELSEIF(CHAR.EQ.'t') THEN
        CHAR='T'
      ELSEIF(CHAR.EQ.'u') THEN
        CHAR='U'
      ELSEIF(CHAR.EQ.'v') THEN
        CHAR='V'
      ELSEIF(CHAR.EQ.'w') THEN
        CHAR='W'
      ELSEIF(CHAR.EQ.'x') THEN
        CHAR='X'
      ELSEIF(CHAR.EQ.'y') THEN
        CHAR='Y'
      ELSEIF(CHAR.EQ.'z') THEN
        CHAR='Z'
      ENDIF
1002  FORMAT(A1)
      IF(CHAR.EQ.'T') GOTO 2810
      IF(CHAR.EQ.'W') GOTO 2820
      IF(CHAR.EQ.'H') GOTO 2821
      IF(CHAR.EQ.'P') GOTO 2822
      IF(CHAR.EQ.'J') GOTO 2823
      IF(CHAR.EQ.'Z') GOTO 2900
      IF(CHAR.EQ.'U') GOTO 2910
      IF(CHAR.EQ.'S') GOTO 3000
      IF(CHAR.EQ.'V') GOTO 3050
      IF(CHAR.EQ.'K') GOTO 3060
      IF(CHAR.EQ.'D') GOTO 3100
      IF(CHAR.EQ.'C') GOTO 3200
      IF(CHAR.EQ.'E') GOTO 3300
      IF(CHAR.EQ.'A') GOTO 3400
      IF(CHAR.EQ.'M') GOTO 3500
      IF(CHAR.EQ.'N') GOTO 3600
      IF(CHAR.EQ.'L') GOTO 3700
      IF(CHAR.EQ.'R') GOTO 3800
      IF(CHAR.EQ.'Q') GOTO 3900
      IF(CHAR.EQ.'B') GOTO 4000
      IF(CHAR.EQ.'I') GOTO 4100
      IF(CHAR.EQ.'Y') GOTO 4200
      IF(CHAR.EQ.'O') GOTO 4300
      IF(CHAR.EQ.'F') GOTO 4400
      IF(CHAR.EQ.'G') GOTO 4450
      IF(CHAR.EQ.'X') GOTO 4460
C
      CALL SETBEL(2)
      CALL RINGBE
      GOTO 2100
C
2810  CONTINUE
C TOGGLE AUGMENTED ABLATION ON(1) OFF(0)
      IF(.NOT.ITOGG) THEN
        ITOGG=.TRUE.
c       TFLAG =ITOGG    !  couple fixed temp flag with augmnted acc toggle
        print *,'current TFRACT:',TFRACT
        READ(*,*) TFRACT
        WRITE(99,*) TFRACT
        DO I=1,NUMNP
          ADOT(I)=TFRACT*ACC(I)-ABLAT(I)
        ENDDO
      ELSE
        ITOGG=.FALSE.
c       TFLAG =ITOGG    !  couple fixed temp flag with augmnted
c       TFRACT=1
        DO I=1,NUMNP
          ADOT(I)=TFRACT*ACC(I)-ABLAT(I)
        ENDDO
      ENDIF
      GOTO 2100
C
2820  CONTINUE
C TOGGLE MASS BALANCE SOLVER ON(.TRUE.) OFF(.FALSE.)
      CTOGG=.NOT. CTOGG
      GOTO 2100
C
2821  CONTINUE
C TOGGLE WATER SOLVER 0,1,2,3
      WTOGG=WTOGG+1
      IF(WTOGG.GT.3) WTOGG=0
      GOTO 2100
C
2822  CONTINUE
C TOGGLE PLATE-BED SOLVER 0,1,2,3
      BTOGG=BTOGG+1
      IF(BTOGG.GT.3) BTOGG=0
      GOTO 2100
2823  CONTINUE
C TOGGLE RUNNING I/O ON(.TRUE.) OFF(.FALSE.)
      IOTOGG=.NOT. IOTOGG
      GOTO 2100
C
2900  CONTINUE
C SELECT A REGION FOR ZOOM
      IF(IPASS.EQ.0) THEN
        WRITE(*,*) 'MUST HAVE GRAPH SHOWING'
        GOTO 2100
      ENDIF
      WRITE(*,*) 'SELECT REGION FOR ZOOM, LOWER LEFT, THEN UPPER RIGHT'
      CALL SETGIN(1)
      CALL SETINK(1)
      CALL SETRUB(1)
      IF(BATCH) THEN
        DO M=1,2
          READ(*,*) XA(M),YA(M),IDAT(M)
        ENDDO
        IGOT=2
      ELSE
        CALL LOCATE(2,XA,YA,IDAT,IGOT)
      ENDIF
      DO M=1,2
        WRITE(99,*) XA(M),YA(M),IDAT(M)
      ENDDO
      IF(IGOT.NE.2) GOTO 2100
      zflag=.true.
      PXMIN=dble(XA(1))
      PYMIN=dble(YA(1))
      PXMAX=dble(XA(2))
      PYMAX=dble(YA(2))
C     CALL ZOOM(REAL(PXMIN),REAL(PXMAX),REAL(PYMIN),REAL(PYMAX))
C     CALL WINDOW(REAL(PXMIN),REAL(PXMAX),REAL(PYMIN),REAL(PYMAX))
      call newpag
      print *,REAL(PXMIN),REAL(PXMAX),REAL(PYMIN),REAL(PYMAX)
C     CALL INUMBR(IERRNM(.TRUE.),10)
      GOTO 2100
C
2910  CONTINUE
C UNZOOM
C     CALL ZOOM(REAL(XMIN),REAL(XMAX),REAL(YMIN),REAL(YMAX))
c      CALL WINDOW(0.,100.,0.,100.)
      zflag=.false.
      if(ipass.ne.0) call newpag
      print *, REAL(XMIN),REAL(XMAX),REAL(YMIN),REAL(YMAX)
      GOTO 2100
C
3000  CONTINUE
C SELECT A REGION FOR MODIFICATION
      IF(IPASS.EQ.0) THEN
        WRITE(*,*) 'MUST HAVE GRAPH SHOWING'
        GOTO 2100
      ENDIF
      WRITE(*,*) 'SELECT A REGION FOR MODIFICATION'
      WRITE(*,*) ' INPUT <A> FOR ALL '
      WRITE(*,*) ' INPUT <R> FOR THOSE IN RANGE ',bedthresh,bedthresh2
      WRITE(*,*) ' INPUT <C> FOR CIRCULAR, PICK CENTER, RADIUS'
      WRITE(*,*) ' INPUT <+> FOR THOSE ABOVE ',bedthresh
      WRITE(*,*) ' INPUT <-> FOR THOSE BELOW ',bedthresh
      WRITE(*,*) ' INPUT <0> FOR THOSE INSIDE SELECTED REGION'
      WRITE(*,*) ' INPUT <1> FOR THOSE OUTSIDE SELECTED REGION'
      WRITE(*,*) ' INPUT <2> FOR THOSE WITH PRESENT ICE '
      WRITE(*,*) ' INPUT <3> FOR THOSE WITHOUT PRESENT ICE '
      WRITE(*,*) ' INPUT <4> FOR THOSE WITH CALCULATED ICE '
      WRITE(*,*) ' INPUT <5> FOR THOSE WITHOUT CALCULATED ICE '
      WRITE(*,*) ' INPUT <6> FOR THOSE BELOW SEA LEVEL W/O ICE '
      WRITE(*,*) ' INPUT <9> REVERSE THOSE SELECTED '
      READ(*,1002) CHAR3
      WRITE(99,1002) CHAR3
      if(CHAR3.eq.'A' .or. CHAR3.eq.'a') then
        ICOUNT=NUMNP
        DO I=1,NUMNP
          IFP(I)=I
          xxx=xx(i)*1.d-3
          yyy=yy(i)*1.d-3
          CALL MRKCLR(1)
          CALL MARKER(REAL(XXX),REAL(YYY),10)
        ENDDO
      elseif(CHAR3.eq.'C' .or. CHAR3.eq.'c') then
        LDANGLE=.false.
        CALL SETGIN(1)
        CALL SETINK(1)
        CALL SETRUB(1)
        IF(BATCH) THEN
          DO M=1,2
            READ(*,*) XA(M),YA(M),IDAT(M)
          ENDDO
          IGOT=2
        ELSE
          CALL LOCATE(2,XA,YA,IDAT,IGOT)
        ENDIF
        DO M=1,2
          WRITE(99,*) XA(M),YA(M),IDAT(M)
        ENDDO
        ICHAR=IDAT(2)
        DO I=1,2
          xxx=dble(xa(i))*1.d3
          yyy=dble(ya(i))*1.d3
          XCHECK(I)=XXX
          YCHECK(I)=YYY
        ENDDO
        radius=(xcheck(1)-xcheck(2))**2+(ycheck(1)-ycheck(2))**2
        ICOUNT=0
        DO I=1,NUMNP
C FIND RETURNS 1 IF FOUND, 0 IF NOT
          dist=(XX(I)-xcheck(1))**2+(YY(I)-ycheck(1))**2
          LFOUND=dist.le.radius
          IF(LFOUND) THEN
            ICOUNT=ICOUNT+1
            IFP(ICOUNT)=I
            xxx=xx(i)*1.d-3
            yyy=yy(i)*1.d-3
            CALL MRKCLR(1)
            CALL MARKER(REAL(XXX),REAL(YYY),10)
          ENDIF
        ENDDO
      elseif(CHAR3.eq.'9') then
        LDANGLE=.true.
        JCOUNT=ICOUNT
        DO J=1,JCOUNT
          JFP(J)=IFP(J)
        ENDDO
        ICOUNT=0
        DO I=1,NUMNP
C FIND RETURNS 1 IF FOUND, 0 IF NOT
          IFOUND=0
          DO J=1,JCOUNT
            IF(I.EQ.JFP(J)) THEN
              IFOUND=0
              GOTO 3001
            ENDIF
          ENDDO
          IFOUND=1
3001      CONTINUE
          IF(IFOUND.EQ.1) THEN
            ICOUNT=ICOUNT+1
            IFP(ICOUNT)=I
            xxx=xx(i)*1.d-3
            yyy=yy(i)*1.d-3
            CALL MRKCLR(1)
            CALL MARKER(REAL(XXX),REAL(YYY),10)
          ENDIF
        ENDDO
      elseif(CHAR3.eq.'R' .or. CHAR3.eq.'r') then
        WRITE(*,*) 'THRESHHOLDS:',bedthresh,bedthresh2
        read(*,*) bedthresh,bedthresh2
        write(99,*) bedthresh,bedthresh2
        LDANGLE=.true.
        ICOUNT=0
        DO I=1,NUMNP
C FIND RETURNS 1 IF FOUND, 0 IF NOT
          IFOUND=0
          IF(DEPB(i).gt.bedthresh .and. 
     &       DEPB(i).lt.bedthresh2) IFOUND=1
          IF(IFOUND.EQ.1) THEN
            ICOUNT=ICOUNT+1
            IFP(ICOUNT)=I
            xxx=xx(i)*1.d-3
            yyy=yy(i)*1.d-3
            CALL MRKCLR(1)
            CALL MARKER(REAL(XXX),REAL(YYY),10)
          ENDIF
        ENDDO
      elseif(CHAR3.eq.'-') then
        WRITE(*,*) 'THRESHHOLD:',bedthresh
        read(*,*) bedthresh
        write(99,*) bedthresh
        LDANGLE=.true.
        ICOUNT=0
        DO I=1,NUMNP
C FIND RETURNS 1 IF FOUND, 0 IF NOT
          IFOUND=0
          IF(DEPB(i).lt.bedthresh) IFOUND=1
          IF(IFOUND.EQ.1) THEN
            ICOUNT=ICOUNT+1
            IFP(ICOUNT)=I
            xxx=xx(i)*1.d-3
            yyy=yy(i)*1.d-3
            CALL MRKCLR(1)
            CALL MARKER(REAL(XXX),REAL(YYY),10)
          ENDIF
        ENDDO
      elseif(CHAR3.eq.'+') then
        WRITE(*,*) 'THRESHHOLD:',bedthresh
        read(*,*) bedthresh
        write(99,*) bedthresh
        LDANGLE=.true.
        ICOUNT=0
        DO I=1,NUMNP
C FIND RETURNS 1 IF FOUND, 0 IF NOT
          IFOUND=0
          !IF(DEPB(i).ge.bedthresh) IFOUND=1
          IF(BDROCK(i).ge.bedthresh) IFOUND=1
          IF(IFOUND.EQ.1) THEN
            ICOUNT=ICOUNT+1
            IFP(ICOUNT)=I
            xxx=xx(i)*1.d-3
            yyy=yy(i)*1.d-3
            CALL MRKCLR(1)
            CALL MARKER(REAL(XXX),REAL(YYY),10)
          ENDIF
        ENDDO
      elseif(CHAR3.eq.'6') then
        LDANGLE=.true.
        ICOUNT=0
        DO I=1,NUMNP
C FIND RETURNS 1 IF FOUND, 0 IF NOT
          IFOUND=0
          IF(DEPB(I).LT.SEALEV) THEN
              FLOT=(1.d0-RATDEN)*(DEPB(I)-SEALEV)
              IF(FLOT.GE.PSURF(I)) IFOUND=1
          ENDIF
          IF(IFOUND.EQ.1) THEN
            ICOUNT=ICOUNT+1
            IFP(ICOUNT)=I
            xxx=xx(i)*1.d-3
            yyy=yy(i)*1.d-3
            CALL MRKCLR(1)
            CALL MARKER(REAL(XXX),REAL(YYY),10)
          ENDIF
        ENDDO
      elseif(CHAR3.eq.'2') then
c ..... THOSE WITH PRESENT ICE '
        LDANGLE=.true.
        ICOUNT=0
        DO I=1,NUMNP
C FIND RETURNS 1 IF FOUND, 0 IF NOT
          IFOUND=0
          IF(BDROCK(I).LT.SEALEV) THEN
            FLOT=(1.d0-RATDEN)*(BDROCK(I)-SEALEV)
            IF(FLOT.LT.PSURF(I)) IFOUND=1
          ELSE
            IF(PSURF(I).GT.BDROCK(I)) IFOUND=1
          ENDIF
          IF(IFOUND.EQ.1) THEN
            ICOUNT=ICOUNT+1
            IFP(ICOUNT)=I
            xxx=xx(i)*1.d-3
            yyy=yy(i)*1.d-3
            CALL MRKCLR(1)
            CALL MARKER(REAL(XXX),REAL(YYY),10)
          ENDIF
        ENDDO
      elseif(CHAR3.eq.'3') then
c ..... THOSE WTHOUT PRESENT ICE '
        LDANGLE=.true.
        ICOUNT=0
        DO I=1,NUMNP
C FIND RETURNS 1 IF FOUND, 0 IF NOT
          IFOUND=0
          IF(BDROCK(I).LT.SEALEV) THEN
            FLOT=(1.d0-RATDEN)*(BDROCK(I)-SEALEV)
            IF(FLOT.GE.PSURF(I)) IFOUND=1
          ELSE
            IF(PSURF(I).LE.BDROCK(I)) IFOUND=1
          ENDIF
          IF(IFOUND.EQ.1) THEN
            ICOUNT=ICOUNT+1
            IFP(ICOUNT)=I
            xxx=xx(i)*1.d-3
            yyy=yy(i)*1.d-3
            CALL MRKCLR(1)
            CALL MARKER(REAL(XXX),REAL(YYY),10)
          ENDIF
        ENDDO
      elseif(CHAR3.eq.'4') then
c ..... THOSE WITH CALCULATED ICE '
        LDANGLE=.true.
        ICOUNT=0
        DO I=1,NUMNP
C FIND RETURNS 1 IF FOUND, 0 IF NOT
          IFOUND=0
          IF(DEPB(I).LT.SEALEV) THEN
            FLOT=(1.d0-RATDEN)*(DEPB(I)-SEALEV)
            IF(FLOT.LT.HTICE(I)) IFOUND=1
          ELSE
            IF(HTICE(I).GT.DEPB(I)) IFOUND=1
          ENDIF
          IF(IFOUND.EQ.1) THEN
            ICOUNT=ICOUNT+1
            IFP(ICOUNT)=I
            xxx=xx(i)*1.d-3
            yyy=yy(i)*1.d-3
            CALL MRKCLR(1)
            CALL MARKER(REAL(XXX),REAL(YYY),10)
          ENDIF
        ENDDO
      elseif(CHAR3.eq.'5') then
c ..... THOSE WITHOUT CALCULATED ICE '
        LDANGLE=.true.
        ICOUNT=0
        DO I=1,NUMNP
C FIND RETURNS 1 IF FOUND, 0 IF NOT
          IFOUND=0
          IF(DEPB(I).LT.SEALEV) THEN
            FLOT=(1.d0-RATDEN)*(DEPB(I)-SEALEV)
            IF(FLOT.GE.HTICE(I)) IFOUND=1
          ELSE
            IF(HTICE(I).LE.DEPB(I)) IFOUND=1
          ENDIF
          IF(IFOUND.EQ.1) THEN
            ICOUNT=ICOUNT+1
            IFP(ICOUNT)=I
            xxx=xx(i)*1.d-3
            yyy=yy(i)*1.d-3
            CALL MRKCLR(1)
            CALL MARKER(REAL(XXX),REAL(YYY),10)
          ENDIF
        ENDDO
      else
        LDANGLE=.false.
        CALL SETGIN(1)
        CALL SETINK(1)
        CALL SETRUB(1)
        IF(BATCH) THEN
          DO M=1,4
            READ(*,*) XA(M),YA(M),IDAT(M)
          ENDDO
          IGOT=4
        ELSE
          CALL LOCATE(4,XA,YA,IDAT,IGOT)
        ENDIF
        DO M=1,4
          WRITE(99,*) XA(M),YA(M),IDAT(M)
        ENDDO
        ICHAR=IDAT(4)
        DO I=1,4
          xxx=dble(xa(i))*1.d3
          yyy=dble(ya(i))*1.d3
          XCHECK(I)=XXX
          YCHECK(I)=YYY
        ENDDO
        ICOUNT=0
        DO I=1,NUMNP
C FIND RETURNS 1 IF FOUND, 0 IF NOT
          IFOUND=IFIND(XX(I),YY(I),XCHECK,YCHECK)
          IF(IFOUND.EQ.1 .and. char3.eq.'0') THEN
            ICOUNT=ICOUNT+1
            IFP(ICOUNT)=I
            xxx=xx(i)*1.d-3
            yyy=yy(i)*1.d-3
            CALL MRKCLR(1)
            CALL MARKER(REAL(XXX),REAL(YYY),10)
          ENDIF
          IF(IFOUND.EQ.0 .and. char3.eq.'1') THEN
            ICOUNT=ICOUNT+1
            IFP(ICOUNT)=I
            xxx=xx(i)*1.d-3
            yyy=yy(i)*1.d-3
            CALL MRKCLR(1)
            CALL MARKER(REAL(XXX),REAL(YYY),10)
          ENDIF
        ENDDO
      endif
      GOTO 2100
C
3050  CONTINUE
C PLOT VOLUME VS TIME
      XMIN=1.d30
      XMAX=-XMIN
      VMIN=XMIN
      VMAX=-VMIN
      AMIN=XMIN
      AMAX=-AMIN
      TMIN=XMIN
      TMAX=-TMIN
      TAMIN=XMIN
      TAMAX=-TMIN
      TTMIN=XMIN
      TTMAX=-TTMIN
      SLMIN=XMIN
      SLMAX=-TTMIN
      WTMIN=XMIN
      WTMAX=-WTMIN
      PWMIN=XMIN
      PWMAX=-PWMIN
      WWMIN=XMIN
      WWMAX=-PWMIN
      DO I=2,NTSTEP
        XMAX=MAX(XMAX,TTIME(I))
        VMAX=MAX(VMAX,VVOL(I))
        AMAX=MAX(AMAX,AAREA(I))
        TMAX=MAX(TMAX,TTBOT(I))
        TAMAX=MAX(TAMAX,TTAVG(I))
        TTMAX=MAX(TTMAX,TTNSL(I))
        SLMAX=MAX(SLMAX,TTSEAL(I))
        WTMAX=MAX(WTMAX,TWATER(I))
        PWMAX=MAX(PWMAX,PWATER(I))
        WWMAX=MAX(WWMAX,WWWMIN(I))
        XMIN=MIN(XMIN,TTIME(I))
        VMIN=MIN(VMIN,VVOL(I))
        AMIN=MIN(AMIN,AAREA(I))
        TMIN=MIN(TMIN,TTBOT(I))
        TAMIN=MIN(TAMIN,TTAVG(I))
        TTMIN=MIN(TTMIN,TTNSL(I))
        SLMIN=MIN(SLMIN,TTSEAL(I))
        WTMIN=MIN(WTMIN,TWATER(I))
        PWMIN=MIN(PWMIN,PWATER(I))
        WWMIN=MIN(WWMIN,WWWMIN(I))
        IF(NTSTEP-I.LE.10) PRINT 111,TTIME(I),
     &            VVOL(I),AAREA(I),TTBOT(I),
     &            TTNSL(I),TTSEAL(I),TWATER(I),
     &            PWATER(I),
     &            WWWMIN(I)
111     FORMAT(1X,1P9G10.3)
      ENDDO
      PRINT *,XMIN,XMAX
      PRINT *,VMIN,VMAX
      PRINT *,AMIN,AMAX
      PRINT *,TMIN,TMAX
      PRINT *,TTMIN,TTMAX
      PRINT *,SLMIN,SLMAX
      PRINT *,WTMIN,WTMAX
      PRINT *,PWMIN,PWMAX
      PRINT *,WWMIN,WWMAX
      DX=(XMAX-XMIN)/20.d0
      DA=(AMAX-AMIN)/20.d0
      DTT=(TMAX-TMIN)/20.d0
      DTTA=(TAMAX-TAMIN)/20.d0
      DTTT=(TTMAX-TTMIN)/20.d0
      DTSL=(SLMAX-SLMIN)/20.d0
      DWT=(WTMAX-WTMIN)/20.d0
      DPW=(PWMAX-PWMIN)/20.d0
      DV=(VMAX-VMIN)/20.d0
      DWW=(WWMAX-WWMIN)/20.d0
      IF(DX.EQ.0.) DX=1.d0
      IF(DA.EQ.0.) DA=1.d0
      IF(DTT.LT.1E-6) DTT=1.d0
      IF(DTTA.LT.1E-6) DTTA=1.d0
      IF(DWT.EQ.0.) DWT=1.d0
      IF(DPW.EQ.0.) DPW=1.d0
      IF(DTTT.EQ.0.) DTTT=1.d0
      IF(DTSL.EQ.0.) DTSL=1.d0
      IF(DV.EQ.0.) DV=1.d0
      IF(DWW.EQ.0.) DWW=1.d0
      XMAX=XMAX+DX
      XMIN=XMIN-DX
      AMAX=AMAX+DA
      AMIN=AMIN-DA
      VMAX=VMAX+DV
      VMIN=VMIN-DV
      TMAX=TMAX+DTT
      TMIN=TMIN-DTT
      TAMAX=TAMAX+DTTA
      TAMIN=TAMIN-DTTA
      TTMAX=TTMAX+DTTT
      TTMIN=TTMIN-DTTT
      SLMAX=SLMAX+DTSL
      SLMIN=SLMIN-DTSL
      WTMAX=WTMAX+DWT
      WTMIN=WTMIN-DWT
      PWMAX=PWMAX+DPW
      PWMIN=PWMIN-DPW
      WWMAX=WWMAX+DWW
      WWMIN=WWMIN-DWW
      DA=(AMAX-AMIN)
      DTT=(TMAX-TMIN)
      DTTA=(TAMAX-TAMIN)
      DTTT=(TTMAX-TTMIN)
      DTSL=(SLMAX-SLMIN)
      DWT=(WTMAX-WTMIN)
      DPW=(PWMAX-PWMIN)
      DWW=(WWMAX-WWMIN)
      DV=(VMAX-VMIN)
      IF(DX.EQ.0.) DX=1.d0
      IF(DA.EQ.0.) DA=1.d0
      IF(DTT.EQ.0.) DTT=1.d0
      IF(DTTA.EQ.0.) DTTA=1.d0
      IF(DTTT.EQ.0.) DTTT=1.d0
      IF(DTSL.EQ.0.) DTSL=1.d0
      IF(DWT.EQ.0.) DWT=1.d0
      IF(DPW.EQ.0.) DPW=1.d0
      IF(DWW.EQ.0.) DWW=1.d0
      IF(DV.EQ.0.) DV=1.d0
      IF(IPASS.EQ.0) CALL GRSTRT(600,600)
      CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &            REAL(VMIN),REAL(VMAX+5*dv))
      IGRAP=0
      CALL NEWPAG
      CALL LINCLR(1)
      CALL MOVE(REAL(TTIME(NTSTEP)),REAL(VVOL(2)))
      CALL DRAW(REAL(TTIME(2)),REAL(VVOL(2)))
      DO I=2,NTSTEP
        CALL DRAW(REAL(TTIME(I)),REAL(VVOL(I)))
      ENDDO
      CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &            REAL(AMIN-da),REAL(AMAX+4*da))
      CALL LINCLR(2)
      CALL MOVE(REAL(TTIME(NTSTEP)),REAL(AAREA(2)))
      CALL DRAW(REAL(TTIME(2)),REAL(AAREA(2)))
      DO I=2,NTSTEP
        CALL DRAW(REAL(TTIME(I)),REAL(AAREA(I)))
      ENDDO
      CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &            REAL(TAMIN-2*DTTA),REAL(TAMAX+3*DTTA))
      CALL LINCLR(8)
      CALL MOVE(REAL(TTIME(NTSTEP)),REAL(TTAVG(2)))
      CALL DRAW(REAL(TTIME(2)),REAL(TTAVG(2)))
      DO I=2,NTSTEP
        CALL DRAW(REAL(TTIME(I)),REAL(TTAVG(I)))
      ENDDO
      CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &            REAL(TMIN-2*dtt),REAL(TMAX+3*dtt))
      CALL LINCLR(3)
      CALL MOVE(REAL(TTIME(NTSTEP)),REAL(TTBOT(2)))
      CALL DRAW(REAL(TTIME(2)),REAL(TTBOT(2)))
      DO I=2,NTSTEP
        CALL DRAW(REAL(TTIME(I)),REAL(TTBOT(I)))
      ENDDO
      CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &            REAL(TTMIN-3*dttt),REAL(TTMAX+2*dttt))
      CALL LINCLR(4)
      CALL MOVE(REAL(TTIME(NTSTEP)),REAL(TTNSL(2)))
      CALL DRAW(REAL(TTIME(2)),REAL(TTNSL(2)))
      DO I=2,NTSTEP
        CALL DRAW(REAL(TTIME(I)),REAL(TTNSL(I)))
      ENDDO
      CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &            REAL(SLMIN-3*dtsl),REAL(SLMAX+2*dtsl))
      CALL LINCLR(11)
      CALL MOVE(REAL(TTIME(NTSTEP)),REAL(TTSEAL(2)))
      CALL DRAW(REAL(TTIME(2)),REAL(TTSEAL(2)))
      DO I=2,NTSTEP
        CALL DRAW(REAL(TTIME(I)),REAL(TTSEAL(I)))
      ENDDO
      CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &            REAL(WTMIN-4*dwt),REAL(WTMAX+1*dwt))
      CALL LINCLR(5)
      CALL MOVE(REAL(TTIME(NTSTEP)),REAL(TWATER(2)))
      CALL DRAW(REAL(TTIME(2)),REAL(TWATER(2)))
      DO I=2,NTSTEP
        CALL DRAW(REAL(TTIME(I)),REAL(TWATER(I)))
      ENDDO
      CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &            REAL(PWMIN-4*DPW),REAL(PWMAX+1*DPW))
      CALL LINCLR(6)
      CALL MOVE(REAL(TTIME(NTSTEP)),REAL(PWATER(2)))
      CALL DRAW(REAL(TTIME(2)),REAL(PWATER(2)))
      DO I=2,NTSTEP
        CALL DRAW(REAL(TTIME(I)),REAL(PWATER(I)))
      ENDDO
      CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &            REAL(WWMIN-5*DWW),REAL(WWMAX+0*DWW))
      CALL LINCLR(7)
      CALL MOVE(REAL(TTIME(NTSTEP)),REAL(WWWMIN(2)))
      CALL DRAW(REAL(TTIME(2)),REAL(WWWMIN(2)))
      DO I=2,NTSTEP
        CALL DRAW(REAL(TTIME(I)),REAL(WWWMIN(I)))
      ENDDO
C      CALL GRSTOP
C      IPASS=0
      GOTO 2100
C
3060  CONTINUE
C PLOT TLIST FILE VS TIME
      CALL DTLIST(IPASS,TIME,AMASS(9))
      GOTO 2100
C
3100  CONTINUE
C DISPLAY ALL FOR CHOSEN REGION
      WRITE(*,*) ' DISPLAY ALL FOR CHOSEN REGION'
      WRITE(*,1050) ICOUNT
      IF(ICOUNT.GT.100) THEN
        WRITE(*,*) 'THIS IS ALOT, DO YOU REALLY WANT THIS?'
        if(iflush) call gflush
        READ(*,1002) CHAR
        WRITE(99,1002) CHAR
        IF(CHAR.NE.'Y') GOTO 2100
      ENDIF
1050  FORMAT(' NODAL VALUES FOR ',I5,' NODES',/,
     &'   N  K      X         Y   SURFACE BED  THICK   ACC   FRACT',
     &' PSURF FLOW SLIDE')
      DO I=1,ICOUNT
        J=IFP(I)
        THIK=HTICE(J)-DEPB(J)
        WRITE(*,1100) J, KODE(J), XX(J), YY(J), HTICE(J), DEPB(J),
     &                THIK, ADOT(J), FRACT(J), PSURF(J), FLOWA(J),
     &                SLDGB(J),slopn(1,j)
1100  FORMAT(1X,I7,I2,1P,2G10.3,0P,3F6.0,2F7.2,F6.0,F6.2,F6.3,1pg12.5)
      ENDDO
      GOTO 2100
C
3200  CONTINUE
C CHANGE VARIOUS PHYSICAL PARAMETERS
      WRITE(*,*) 'CHANGE VARIOUS PHYSICAL PARAMETERS OF SELECTED NODES'
      WRITE(*,*) 'INPUT:'
      WRITE(*,*) 'A: <A>FUDGE'
      WRITE(*,*) 'B: <B>ED'
      WRITE(*,*) 'C: <C>LEAR LINEARIZATION CONSTANT'
      WRITE(*,*) 'D: A<D>OT NET MASS BALANCE'
      WRITE(*,*) 'E: INITIAL TIM<E>'
      WRITE(*,*) 'F: <F>LOW LAW CONSTANT'
      WRITE(*,*) 'G: <G>EOTHERMAL GRADIENT, EFFECTIVE FLOW FUDGE FACTOR'
      WRITE(*,*) 'H: <H>TICE INITIAL ICE SURFACE'
      WRITE(*,*) 'I: W<I>ND DIRECTION AND MAGNITUDE FOR ZONE 21'
      WRITE(*,*) 'J: <J> CONSTANTS IN WATER THICKNESS CALCULATION'
      WRITE(*,*) 'K: <K>ODE FIXED OR FREE'
      WRITE(*,*) 'L: <L>APSE RATE FOR ZONE 20,21'
      WRITE(*,*) 'M: EXPERI<M>ENT'
      WRITE(*,*) 'N: S<N>OW LINE ELEVATION FOR ZONES 13..18'
      WRITE(*,*) 'O: P<O>LE POSITION'
      WRITE(*,*) 'P: <P>ERCENT OF VELOCITY DUE TO SLIDING'
      WRITE(*,*) 'Q: TOGGLE ELEMENT(1)/NODE(0) ACON <Q>?',ICON
      WRITE(*,*) 'R: WATE<R> THICKNESS '
      WRITE(*,*) 'S: <S>LIDING LAW CONSTANT'
      WRITE(*,*) 'T: <T>IME STEP NUMBER, OUTPUT STEP, TIME STEP'
      WRITE(*,*) 'U: PS<U>RF PRESENT SURFACE'
      WRITE(*,*) 'V: CAL<V>ING FACTOR'
      WRITE(*,*) 'W: SNO<W>LINE GRADIENTS FOR ZONES 20 AND 21'
      WRITE(*,*) 'X: FLU<X> MODIFICATION'
      WRITE(*,*) 'Y: CHANGE TEMPERATURE CALCULATION T<Y>PE'
      IF(IMELT.EQ.1) THEN
        WRITE(*,*) 'Z: <Z> TOGGLE BASAL MELT CALCULATION OFF'
      ELSE
        WRITE(*,*) 'Z: <Z> TOGGLE BASAL MELT CALCULATION ON'
      ENDIF
      WRITE(*,*) '1: CHANGE TEMPERATURE NEAR SEA LEVEL BASE'
      WRITE(*,*) '2: CHANGE ADVANCE RATE PARAMETER (0-1):',ADVANCE
      IF(htogg) THEN
        WRITE(*,*) '3: <3> INTERNAL HEAT GENERATION OFF'
      ELSE
        WRITE(*,*) '3: <3> INTERNAL HEAT GENERATION ON'
      ENDIF
      WRITE(*,*) '4: CHANGE SEA LEVEL:',SEALEV
      WRITE(*,*) '5: <5> SEALEVEL CONSTANT:                  0',SLTOGG
      WRITE(*,*) '       SEALEVEL FROM EXTERNAL FILE:        1'
      WRITE(*,*) '       SEALEVEL FROM INTERNAL CALCULATION: 2'
      IF(HIRES) THEN
        WRITE(*,*) '6: <6> SET LO-RES ON',HIRES
      ELSE
        WRITE(*,*) '6: <6> SET HI-RES ON',HIRES
      ENDIF
      if(iflush) call gflush
      READ(*,1002) CHAR
      WRITE(99,1002) CHAR
      IF(CHAR.EQ.'k') THEN
        CHAR='K'
      ELSEIF(CHAR.EQ.'j') THEN
        CHAR='J'
      ELSEIF(CHAR.EQ.'z') THEN
        CHAR='Z'
      ELSEIF(CHAR.EQ.'h') THEN
        CHAR='H'
      ELSEIF(CHAR.EQ.'a') THEN
        CHAR='A'
      ELSEIF(CHAR.EQ.'y') THEN
        CHAR='Y'
      ELSEIF(CHAR.EQ.'g') THEN
        CHAR='G'
      ELSEIF(CHAR.EQ.'i') THEN
        CHAR='I'
      ELSEIF(CHAR.EQ.'d') THEN
        CHAR='D'
      ELSEIF(CHAR.EQ.'u') THEN
        CHAR='U'
      ELSEIF(CHAR.EQ.'p') THEN
        CHAR='P'
      ELSEIF(CHAR.EQ.'f') THEN
        CHAR='F'
      ELSEIF(CHAR.EQ.'s') THEN
        CHAR='S'
      ELSEIF(CHAR.EQ.'b') THEN
        CHAR='B'
      ELSEIF(CHAR.EQ.'c') THEN
        CHAR='C'
      ELSEIF(CHAR.EQ.'t') THEN
        CHAR='T'
      ELSEIF(CHAR.EQ.'v') THEN
        CHAR='V'
      ELSEIF(CHAR.EQ.'m') THEN
        CHAR='M'
      ELSEIF(CHAR.EQ.'n') THEN
        CHAR='N'
      ELSEIF(CHAR.EQ.'w') THEN
        CHAR='W'
      ELSEIF(CHAR.EQ.'x') THEN
        CHAR='X'
      ELSEIF(CHAR.EQ.'l') THEN
        CHAR='L'
      ELSEIF(CHAR.EQ.'e') THEN
        CHAR='E'
      ELSEIF(CHAR.EQ.'o') THEN
        CHAR='O'
      ELSEIF(CHAR.EQ.'q') THEN
        CHAR='Q'
      ELSEIF(CHAR.EQ.'r') THEN
        CHAR='R'
      ENDIF
      IF(CHAR.EQ.'K') THEN
        WRITE(*,*) 'INPUT NEW KODE FOR SELECTED NODES'
        WRITE(*,*) '-3 FOR FIX ALL WITH NO CURRENT ICE'
        WRITE(*,*) '-2 FOR FIX ALL GRID BOUNDARY, FREE ALL INTERIOR'
        WRITE(*,*) '-1 FOR FIX ALL BELOW FLOTATION LINE'
        WRITE(*,*) '   AND SET AFUDGE TO 0.01'
        WRITE(*,*) ' 0 FOR CALCULATED (FREE)'
        WRITE(*,*) ' 1 FOR SPECIFIED (FIXED)'
        if(iflush) call gflush
        READ(*,*) KODEN
        WRITE(99,*) KODEN
        IF(KODEN.GE.0) THEN
          DO I=1,ICOUNT
            KODE(IFP(I))=KODEN
          ENDDO
        ELSEIF(KODEN.EQ.-3) THEN
          DO I=1,ICOUNT
            IF(PSURF(IFP(I)).LE.BDROCK(IFP(I))) THEN
              KODE(IFP(I))=1
              AFUDGE(IFP(I))=0.01d0
            ELSE
              KODE(IFP(I))=0
            ENDIF
          ENDDO
        ELSEIF(KODEN.EQ.-2) THEN
          N=1
          DO J=1,NUMLEV
            DO I=1,NUMCOL
              IF(I.EQ.1 .OR. I.EQ.NUMCOL) THEN
                KODE(N)=1
              ELSEIF(J.EQ.1 .OR. J.EQ.NUMLEV) THEN
                KODE(N)=1
              ELSE
                KODE(N)=0
              ENDIF
              N=N+1
            ENDDO
          ENDDO
        ELSEIF(KODEN.EQ.-1) THEN
          DO I=1,ICOUNT
            IF(DEPB(IFP(I)).LT.SEALEV) THEN
              FLOT=(1.d0-RATDEN)*(DEPB(IFP(I))-SEALEV)
              IF(FLOT+25.GE.PSURF(IFP(I))) THEN
                KODE(IFP(I))=1
                AFUDGE(IFP(I))=0.01d0
              ENDIF
            ENDIF
          ENDDO
        ENDIF  
        GOTO 2100
      ELSEIF(CHAR.EQ.'1') THEN
        WRITE(*,*) 'INPUT NEW TNSLBASE',TNSLBASE
        if(iflush) call gflush
        READ(*,*) TNSLBASE
        WRITE(99,*) TNSLBASE
        DO I=1,NUMNP
          ADOT(I)=AFUNCT(TIME, IDT(I), AMASS, HTICE(I),
     &                   BDROCK(I), PSURF(I), SLOPN(1,I),
     &                   XX(I),YY(I),TEMP(I))
        ENDDO
        GOTO 2100
      ELSEIF(CHAR.EQ.'3') THEN
        htogg=.not.htogg
        GOTO 2100
      ELSEIF(CHAR.EQ.'5') THEN
        WRITE(*,*) 'SEALEVEL CONSTANT:                   0',SLTOGG
        WRITE(*,*) 'SEALEVEL FROM EXTERNAL FILE:         1'
        WRITE(*,*) 'SEA LEVEL FROM INTERNAL CALCULATION: 2'
        WRITE(*,*) 'INPUT SLTOGG',SLTOGG
        if(iflush) call gflush
        READ(*,*) SLTOGG
        WRITE(99,*) SLTOGG
        GOTO 2100
      ELSEIF(CHAR.EQ.'2') THEN
        WRITE(*,*) 'INPUT NEW ADVANCE RATE PARAMETER',ADVANCE
        if(iflush) call gflush
        READ(*,*) ADVANCE
        WRITE(99,*) ADVANCE
        GOTO 2100
      ELSEIF(CHAR.EQ.'6') THEN
        HIRES=.NOT.HIRES
        GOTO 2100
      ELSEIF(CHAR.EQ.'4') THEN
        SEASAVE=SEALEV
        WRITE(*,*) 'INPUT NEW SEA LEVEL',SEALEV
        if(iflush) call gflush
        READ(*,*) SEALEV
        WRITE(99,*) SEALEV
        DO I=1,NUMNP
          if(abs(htice(i)-seasave).lt.1e-5) then
            htice(i)=max(sealev,depb(i))
          endif
        ENDDO

        GOTO 2100
      ELSEIF(CHAR.EQ.'X') THEN
        RFLUX=0.d0
        IFLUX=0
        DO I=1,NUMGBC
          IF(BFLUX(I).NE.0.0) THEN
            RFLUX=RFLUX+BFLUX(I)
            IFLUX=IFLUX+1
          ENDIF
        ENDDO
        IF(IFLUX.GT.0) RFLUX=RFLUX/dble(IFLUX)
        WRITE(*,*) 'INPUT MULTIPLYING FACTOR FOR FLUXES'
        WRITE(*,*) 'CURRENT AVERAGE FLUX IS',RFLUX
        WRITE(*,*) '-999 TO SMOOTH...'
        if(iflush) call gflush
        READ(*,*) RFLUX
        WRITE(99,*) RFLUX
        IF(RFLUX.NE.-999.) THEN
          DO I=1,NUMGBC
            BFLUX(I)=BFLUX(I)*RFLUX
          ENDDO
        ELSE
          DO I=2,NUMGBC-1
            TMP(I)=(BFLUX(I-1)+BFLUX(I)+BFLUX(I+1))/3.d0
          ENDDO
          TMP(1)=(BFLUX(NUMGBC)+BFLUX(1)+BFLUX(2))/3.d0
          TMP(NUMGBC)=(BFLUX(NUMGBC-1)+BFLUX(NUMGBC)+BFLUX(1))/3.d0
          DO I=1,NUMGBC
            BFLUX(I)=TMP(I)
          ENDDO
        ENDIF
        GOTO 2100
      ELSEIF(CHAR.EQ.'W') THEN
        WRITE(*,*) 'INPUT POLE SNOWLINE ELEVATION'
        WRITE(*,*) ' AND LATITUDE GRADIENT, TNSL'
        WRITE(*,*) AMASS(7),AMASS(8),AMASS(9)
        if(iflush) call gflush
        READ(*,*) AMASS(7),AMASS(8),AMASS(9)
        WRITE(99,*) AMASS(7),AMASS(8),AMASS(9)
        WRITE(73,*) TIME,AMASS(9)
C**********EISMINT ala huybrechts ********************
        WRM=WARMING(AMASS(9)+TNSLBASE)
C*****************************************************
        DO I=1,ICOUNT
          ADOTI=AFUNCT(TIME, IDT(IFP(I)), AMASS, HTICE(IFP(I)),
     &                      BDROCK(IFP(I)), PSURF(IFP(I)),
     &                      SLOPN(1,IFP(I)),
     &                      XX(IFP(I)),YY(IFP(I)),TEMP(IFP(I)))
          IF(IDT(IFP(I)).GT.0) THEN
            ADOT(IFP(I))=ADOTI
          ELSE
            AJUNK=ACCUM(IDT(IFP(I)),XX(IFP(I))*.001D0,
     &                  YY(IFP(I))*.001D0,
     &                  HTICE(IFP(I)),SLOPN(1,IFP(I)),PSURF(I),
     &                  AMASS(9),TEMP(IFP(I)))
            IF(TFLAG) TEMP(I)=TSORIG(I)
            IF(ITOGG) THEN
C             ADOT(I)=ADOTB(I)-ABL*.01
C*****************************************************
C**********EISMINT ala huybrechts ********************
              ADOT(I)=WRM*ADOTB(I)-ABL*.01d0
              IF(TFLAG) ADOT(I)=AMARS(WRM,ADOTB(I),ABL)
C*****************************************************
            ELSE
              ADOT(IFP(I))=ADOTB(IFP(I))
            ENDIF
          ENDIF
        ENDDO
        GOTO 2100
      ELSEIF(CHAR.EQ.'R') THEN
        WRITE(*,*) 'INPUT WATER THICKNESS '
        WTHICKN=0.0d0
        DO I=1,ICOUNT
          WTHICKN=WTHICKN+WTHICK(IFP(I))
        ENDDO
        WTHICKN=WTHICKN/ICOUNT
        WRITE(*,*) 'CURRENT WATER THICKNESS:',REAL(WTHICKN)
        WRITE(*,*) 'INPUT NEW WATER THICKNESS FOR SELECTED NODES'
        if(iflush) call gflush
        READ(*,*) WTHICKN
        WRITE(99,*) WTHICKN
        DO I=1,ICOUNT
          WTHICK(IFP(I))=WTHICKN
        ENDDO
        GOTO 2100
      ELSEIF(CHAR.EQ.'G') THEN
        WRITE(*,*) 'INPUT GEOTHERMAL GRADIENT'
        GEOFLUXN=0.0d0
        DO I=1,ICOUNT
          GEOFLUXN=GEOFLUXN+GEOFLUX(IFP(I))
        ENDDO
        GEOFLUXN=GEOFLUXN/ICOUNT
        WRITE(*,*) 'CURRENT GEOTHERMAL GRADIENT:',REAL(GEOFLUXN)
        WRITE(*,*) 'INPUT NEW GEOTHERMAL GRADIENT FOR SELECTED NODES'
        if(iflush) call gflush
        READ(*,*) GEOFLUXN
        WRITE(99,*) GEOFLUXN
        DO I=1,ICOUNT
          GEOFLUX(IFP(I))=GEOFLUXN
        ENDDO
        GOTO 2100
      ELSEIF(CHAR.EQ.'H') THEN
        WRITE(*,*) 'INPUT NEW HTICE FOR SELECTED NODES'
        WRITE(*,*) '      NEGATIVE SETS TO PRESENT SURFACE'
        WRITE(*,*) '      ZERO SETS TO PRESENT BED'
        if(iflush) call gflush
        READ(*,*) HTICEN
        WRITE(99,*) HTICEN
        DO I=1,ICOUNT
          IF(HTICEN.LT.0.) THEN
            HTICE(IFP(I))=MAX(HTICEN,PSURF(IFP(I)))
          ELSE
            HTICE(IFP(I))=MAX(HTICEN,DEPB(IFP(I)))
          ENDIF
          IF(HTICE(IFP(I)).LT.0.) HTICE(IFP(I))=0.d0
        ENDDO
        GOTO 2100
      ELSEIF(CHAR.EQ.'Z') THEN
        WRITE(*,*) 'INPUT NEW BMELT FOR SELECTED NODES'
        PRINT *,' TO TURN  ON BASAL MELT CALCULATION INPUT  1'
        PRINT *,' TO TURN OFF BASAL MELT CALCULATION INPUT -1'
        if(iflush) call gflush
        READ(*,*) IMELT
        WRITE(99,*) IMELT
        IF(IMELT.EQ.1) THEN
          PRINT *,' TOGGLE BMELT CALC ON...'
        ELSE
          PRINT *,' TOGGLE BMELT CALC OFF: SPECIFY BMELT...'
          if(iflush) call gflush
          READ(*,*) BMELTN
          WRITE(99,*) BMELTN
          DO I=1,ICOUNT
            BMELT(IFP(I))=BMELTN
            PRINT *,I,IFP(I),BMELTN
          ENDDO
        ENDIF
        GOTO 2100
      ELSEIF(CHAR.EQ.'D') THEN
        GOTO 3230
      ELSEIF(CHAR.EQ.'U') THEN
        WRITE(*,*) 'INPUT NEW PRESENT SURFACE FOR SELECTED NODES'
        WRITE(*,*) 'IF NEGATIVE SET ALL BELOW FLOTATION LINE TO 0'
        if(iflush) call gflush
        READ(*,*) PSURFN
        WRITE(99,*) PSURFN
        IF(PSURFN.GE.0) THEN
          DO I=1,ICOUNT
            PSURF(IFP(I))=PSURFN
          ENDDO
        ELSE
          DO I=1,ICOUNT
            FLOT=(1.d0-RATDEN)*(BDROCK(IFP(I))-SEALEV)
            IF(PSURF(IFP(I)).LT.FLOT) THEN
              PSURF(IFP(I))=SEALEV
              AFUDGE(IFP(I))=0.1d0
            ENDIF
          ENDDO
        ENDIF
        GOTO 2100
      ELSEIF(CHAR.EQ.'P') THEN
        WRITE(*,*) 'INPUT NEW FRACT'
        WRITE(*,*) ' (1 FOR ALL SLIDING, 0 FOR NONE)'
        WRITE(*,*) '      NEGATIVE FOR RANDOM...'
        if(iflush) call gflush
        READ(*,*) FRACTN
        WRITE(99,*) FRACTN
        IF(FRACTN.GE.0.) THEN
          DO I=1,ICOUNT
            FRACT(IFP(I))=FRACTN
          ENDDO
        ELSE
          WRITE(*,*) 'INPUT SEED AND THRESHOLD (0-NONE, 1-ALL)',
     &              ISEED,THRESH
          if(iflush) call gflush
          READ(*,*) ISEED,THRESH
          WRITE(99,*) ISEED,THRESH
          CALL SRAND(ISEED)
          DO I=1,ICOUNT
            XRAND=dble(RAND())
            IF(XRAND.GT.THRESH) THEN
              FRACT(IFP(I))=1.d0
            ELSE
              FRACT(IFP(I))=0.d0
            ENDIF
          ENDDO
        ENDIF
        GOTO 2100
      ELSEIF(CHAR.EQ.'F') THEN
        WRITE(*,*) 'INPUT NEW FLOWA FOR SELECTED NODES,'
        WRITE(*,*) '       <0 USE AS MULTIPLIER'
        WRITE(*,*) '   **** NOTE ****'
        WRITE(*,*) '  THIS CHANGES ALL ELEMENT VLAUES OF FLOW'
        WRITE(*,*) '  ***************************************'
        if(iflush) call gflush
        READ(*,*) FLOWAN
        WRITE(99,*) FLOWAN
        IF(FLOWAN.GT.0.) THEN
          DO I=1,ICOUNT
            FLOWA(IFP(I))=FLOWAN
          ENDDO
        ELSE
          DO I=1,ICOUNT
            FLOWA(IFP(I))=-FLOWAN*FLOWA(IFP(I))
          ENDDO
        ENDIF
        IF(FLOWAN.GT.0.) THEN
          DO I=1,NUMEL
            ACON(I)=FLOWAN
          ENDDO
        ELSE
          DO I=1,NUMEL
            ACON(I)=-FLOWAN*ACON(I)
          ENDDO
        ENDIF
        GOTO 2100
      ELSEIF(CHAR.EQ.'S') THEN
        WRITE(*,*) 'INPUT NEW SLDGB FOR SELECTED NODES,'
        WRITE(*,*) '       <0 USE AS MULTIPLIER'
        if(iflush) call gflush
        READ(*,*) SLDGBN
        WRITE(99,*) SLDGBN
        IF(SLDGBN.GT.0.) THEN
          DO I=1,ICOUNT
            SLDGB(IFP(I))=SLDGBN
          ENDDO
        ELSE
          DO I=1,ICOUNT
            SLDGB(IFP(I))=-SLDGBN*SLDGB(IFP(I))
          ENDDO
        ENDIF  
        GOTO 2100
      ELSEIF(CHAR.EQ.'B') THEN
        WRITE(*,*) 'INPUT NEW BED FOR SELECTED NODES'
        if(iflush) call gflush
        READ(*,*) BEDN
        WRITE(99,*) BEDN
        DO I=1,ICOUNT
          BDROCK(IFP(I))=BEDN
        ENDDO
        GOTO 2100
      ELSEIF(CHAR.EQ.'C') THEN
        WRITE(*,*) 'SET ALL CONST=0.00'
        DO I=1,NUMEL
C         CONST(I)=1.92d7
          CONST(I)=0.0d0
        ENDDO
        GOTO 2100
      ELSEIF(CHAR.EQ.'T') THEN
        WRITE(*,*) 'CURRENT NUMBER OF STEPS,OUTPUT INTERVAL,TIME STEP'
        WRITE(*,*) NDT,INTER,DT
        if(iflush) call gflush
        READ(*,*) NDT,INTER,DT
        WRITE(99,*) NDT,INTER,DT
        GOTO 2100
      ELSEIF(CHAR.EQ.'N') THEN
        WRITE(*,*) 'WHICH SNOW LINE DO YOU WANT TO CHANGE'
        WRITE(*,*) ' (-1 CHANGES ALL)'
        WRITE(*,*) '13=>SM,14=>TM,15=>SX,16=>PX,17=>PC,18=>??'
        if(iflush) call gflush
        READ(*,*) IDOT
        WRITE(99,*) IDOT
        IF(IDOT.GE.0) THEN
          WRITE(*,*) 'CURRENT SNOW LINE ',SNO(IDOT)
          if(iflush) call gflush
          READ(*,*) SNO(IDOT)
          WRITE(99,*) SNO(IDOT)
        ELSE
          WRITE(*,*) (I,SNO(I),I=13,18)
          if(iflush) call gflush
          READ(*,*) SNOCON
          WRITE(99,*) SNOCON
          DO I=13,18
            SNO(I)=SNOCON
          ENDDO
        ENDIF
        DO I=1,NUMNP
          ADOTI=AFUNCT(TIME, IDT(I), AMASS,             
     &                   HTICE(I), BDROCK(I), PSURF(I),                
     &                   SLOPN(1,I),XX(I),YY(I),TEMP(I))            
          IF(IDT(I).GT.0) ADOT(I)=ADOTI           
        ENDDO
        GOTO 2100
      ELSEIF(CHAR.EQ.'L') THEN
        WRITE(*,*) 'PPRESENT LAPSE RATE',ACOM                             
        if(iflush) call gflush
        READ(*,*) ACOM                                                    
        WRITE(99,*) ACOM                                                    
        DO I=1,NUMNP                                                 
          ADOTI=AFUNCT(TIME, IDT(I), AMASS,             
     &                 HTICE(I), BDROCK(I), PSURF(I),                
     &                 SLOPN(1,I),XX(I),YY(I),TEMP(I))            
          IF(IDT(I).GT.0) ADOT(I)=ADOTI            
        ENDDO                                                          
        GOTO 2100                                                         
      ELSEIF(CHAR.EQ.'V') THEN
        WRITE(*,*) 'PPRESENT CALVING FACTOR',CFACTOR
        WRITE(*,*) '  (POSITIVE FOR GLOBAL VALUE) '
        WRITE(*,*) '  (NEGATIVE FOR LOCAL VALUES) '
        if(iflush) call gflush
        READ(*,*) CFACTOR                                                    
        WRITE(99,*) CFACTOR
        IF(CFACTOR.LT.0.) THEN
          CAVG=0.0d0
          DO I=1,ICOUNT
            CAVG=CAVG+CALV(IFP(I))
          ENDDO
          IF(ICOUNT.GT.0) CAVG=CAVG/ICOUNT
          WRITE(*,*) 'LOCAL CALVING FACTOR FOR SPECIFIED NODES',CAVG
          WRITE(*,*) '      (NEGATIVE TO MULTIPLY... '
          if(iflush) call gflush
          READ(*,*) CAVG                              
          write(99,*) CAVG                              
          IF(CAVG.GE.0.0) THEN
            DO I=1,ICOUNT
              CALV(IFP(I))=CAVG
            ENDDO
          ELSE
            DO I=1,ICOUNT
              CALV(IFP(I))=ABS(CAVG)*CALV(IFP(I))
            ENDDO
          ENDIF
        ENDIF
        DO I=1,NUMNP                                                 
          ADOT(I)=ADOTB(I)
          ADOTI=AFUNCT(TIME, IDT(I), AMASS,             
     &                 HTICE(I), BDROCK(I), PSURF(I),                
     &                 SLOPN(1,I),XX(I),YY(I),TEMP(I))            
          IF(IDT(I).GT.0) ADOT(I)=ADOTI   
        ENDDO                                                          
        CALL CALVING(MMXX, NUMEL, NTYPE, KX, HTICE, DEPB, ADOT, 
     &           AADOT,CFACTOR,ADC,CALV,PCALV,AMASS(9))
        CALL ELPROP(MMXX, NUMEL, NTYPE, KX, ADOT, 
     &             AADOT, FRACT,
     &             AFRACT, BDROCK, ABDRCK, PSURF, PPSURF, FLOWA,ACON,
     &             ICON, AFLOWA, SLDGB, ASLDGB,CALV,PCALV)
        GOTO 2100                                                         
      ELSEIF(CHAR.EQ.'J') THEN
        WRITE(*,'(1x,a/1x,a/1x,1p3e13.6)') 
     &         'PPRESENT WATER THICKNESS CONSTANTS',
     &         'advection, diffusion, and leakage ',ALPHAC
        if(iflush) call gflush
        READ(*,*) ALPHAC
        WRITE(99,*) ALPHAC
        GOTO 2100                                                         
      ELSEIF(CHAR.EQ.'I') THEN
        WINDEG=WINDIR(1)*180.d0/3.14159d0
        WINDSP=WINDIR(2)                                                        
        WRITE(*,*) 'PRESENT WIND DIRECTION AND SPEED'
        WRITE(*,*) 'FROM:0-SOUTH,90-WEST,180-NORTH,270-EAST',
     &              WINDEG,WINDSP
        if(iflush) call gflush
        READ(*,*) WINDEG,WINDSP
        WRITE(99,*) WINDEG,WINDSP
        WINDIR(1)=WINDEG/180.d0*3.14159d0
        WINDIR(2)=WINDSP
        CALL SLOPEF(MXX,XX,YY,KX,NUMEL,HTICE,SLOPE,WINDIR)
        CALL NODESL(MXX, NUMNP, NUMEL, KX, SLOPE, SLOPN)
        DO I=1,NUMNP                                                 
          ADOTI=AFUNCT(TIME, IDT(I), AMASS,             
     &                   HTICE(I), BDROCK(I), PSURF(I),                
     &                   SLOPN(1,I),XX(I),YY(I),TEMP(I))            
          IF(IDT(I).GT.0) ADOT(I)=ADOTI           
        ENDDO                                                          
        GOTO 2100
      ELSEIF(CHAR.EQ.'E') THEN
        WRITE(*,*) 'INPUT INITIAL TIME',TIME
        if(iflush) call gflush
        READ(*,*) TIME
        WRITE(99,*) TIME
        TBASE=TIME
        WRITE(73,*) TIME,AMASS(9)
        DO I=1,NUMNP                                                 
          ADOTI=AFUNCT(TIME, IDT(I), AMASS,             
     &                   HTICE(I), BDROCK(I), PSURF(I),                
     &                   SLOPN(1,I),XX(I),YY(I),TEMP(I))            
          IF(IDT(I).GT.0) ADOT(I)=ADOTI           
        ENDDO
        GOTO 2100                                                          
      ELSEIF(CHAR.EQ.'O') THEN
        WRITE(*,*) 'CURRENT POLE POSITION',XPOLE,YPOLE
        if(iflush) call gflush
        READ(*,*) XPOLE,YPOLE
        WRITE(99,*) XPOLE,YPOLE
        DO I=1,NUMNP                                                 
          ADOTI=AFUNCT(TIME, IDT(I), AMASS,             
     &                   HTICE(I), BDROCK(I), PSURF(I),                
     &                   SLOPN(1,I),XX(I),YY(I),TEMP(I))            
          IF(IDT(I).GT.0) ADOT(I)=ADOTI           
        ENDDO
        GOTO 2100                                                          
      ELSEIF(CHAR.EQ.'Y') THEN
        WRITE(*,*) 'INPUT: 0 - NO TEMPERATURE CALCULATION'
        WRITE(*,*) '       1 - MODIFY FRACTION ONLY'
        WRITE(*,*) '       2 - MODIFY FLOW     ONLY'
        WRITE(*,*) '       3 - MODIFY FLOW AND FRACTION BY TEMP'
        WRITE(*,*) '       4 - MODIFY FLOW AND SLIDING'
        WRITE(*,*) '       5 - MODIFY FLOW, SLIDING, AND FRACTION'
        WRITE(*,*) '       6 - MODIFY AFUDGE TO FIT FLOW'
        WRITE(*,*) '       7 - MODIFY FLOW BY TEMP '
        WRITE(*,*) '                  AND FRACTION BY WATER'
        WRITE(*,*) '       8 - MODIFY FLOW BY TEMP '
        WRITE(*,*) '                  AND LUBRICATION FROM WATER'
        WRITE(*,*) '                  AND FRACT FROM FLOW AND SLIDING V''S'

        if(iflush) call gflush
        READ(*,*) ITY
        WRITE(99,*) ITY
        DO I=1,ICOUNT
          ITYPE(IFP(I))=ITY
        ENDDO
        GOTO 2100                                                          
      ELSEIF(CHAR.EQ.'M') THEN
        CALL SCENARIO(TIME)
        GOTO 2100                                                          
      ELSEIF(CHAR.EQ.'A') THEN
        WRITE(*,*) 'INPUT NEW AFUDGE FOR SELECTED NODES'
        if(iflush) call gflush
        READ(*,*) AFUDGEN
        WRITE(99,*) AFUDGEN
        AMASS(10)=AFUDGEN
        DO I=1,ICOUNT
          AFUDGE(IFP(I))=AFUDGEN
        ENDDO
        GOTO 2100
      ELSEIF(CHAR.EQ.'Q') THEN
        WRITE(*,*) ' TOGGLE ACON '
        IF(ICON.EQ.1) THEN
          ICON=0
        ELSEIF(ICON.EQ.0) THEN
          ICON=1
        ENDIF
        GOTO 2100
      ELSE
C
        GOTO 3200
      ENDIF
C
3230  CONTINUE
      WRITE(*,*) 'INPUT NEW ADOT FOR SELECTED NODES (-999 USES FIT)'
      if(iflush) call gflush
      READ(*,*) ADOTN
      WRITE(99,*) ADOTN
      IDOT=0
      IF(ADOTN.EQ.-999) THEN
3231    CONTINUE
        WRITE(*,*) '0 - MODIFY'
        WRITE(*,*) '1 - MARS, CIRCULAR DEFINED'
        WRITE(*,*) '2 - MARS, USER DEFINED ELA, PEAK, ETC.'
        WRITE(*,*) 'LINEAR POLAR <= 3-B, 4-C'
        WRITE(*,*) 'LINEAR MARITIME <= 5-A, 6-B, 7-C'
        WRITE(*,*) 'SLOPE DEPENDENT 8=>SM,9=>TM,10=>SX,11=>PX,12=>PC'
      WRITE(*,*) 'EXPONENTIAL 13=>SM,14=>TM,15=>SX,16=>PX,17=>PC,18=ANT'
        WRITE(*,*) '19=>GROSSWALD,20=>METEOROLOGICAL'
        WRITE(*,*) '21=>WIND/SLOPE METHOD'
        WRITE(*,*) '22=>HUYBRECHTS METHOD'
        WRITE(*,*) '23=>WHITE HOLE'
        WRITE(*,*) '24=>MILANKOVICH'
        WRITE(*,*) '25=>ISLSCP DATA-BASED METHOD'
        WRITE(*,*) '26=>NCEP2 DATA-BASED METHOD'
        WRITE(*,*) '27=>LGM-GCM RESULT-BASED METHOD'
        WRITE(*,*) '28=>S.L-DEPENDENT INTERPOLATION OF '
        WRITE(*,*) '  NCEP2 DATA-BASED AND LGM-GCM RESULT-BASED METHOD'
        WRITE(*,*) '29=>NCEP2 DATA-BASED METHOD TOPO-MOD'
        WRITE(*,*) '30=>LGM-GCM RESULT-BASED METHOD TOPO-MOD'
        WRITE(*,*) '31=>S.L-DEPENDENT INTERPOLATION OF '
        WRITE(*,*) '  NCEP2 DATA-BASED AND LGM-GCM RESULT-BASED METHOD'
        WRITE(*,*) '  WITH TOPO-MOD'
  
        if(iflush) call gflush
        READ(*,*) IDOT
        WRITE(99,*) IDOT
        IF(IDOT.LT.0 .OR. IDOT.GT.31) GOTO 3231
        IF(IDOT.EQ.1) THEN
            WRITE(*,*) 'INPUT X,Y,RADIUS,PEAK'
            WRITE(*,*) (real(AMASS(I)),I=1,4)
            if(iflush) call gflush
            READ(*,*) (AMASS(I),I=1,4)
            WRITE(99,*) (AMASS(I),I=1,4)
        ENDIF
        IF(IDOT.EQ.2) THEN
3232      CONTINUE
331       FORMAT(3F10.4,3F10.0)
c          READ(*,*) (AMASS(I),I=1,6)
c          WRITE(99,*) (AMASS(I),I=1,6)
            if(idot.eq.2) WRITE(*,*) 'INPUT AMARG,APEAK,ADOME,HE,HP,HM'
            WRITE(*,*) (real(AMASS(I)),I=1,6)
            WRITE(*,*) '      1 TO CHANGE ELA'
            WRITE(*,*) '      2 TO CHANGE PEAK AND DOME ACCUMULATIONS'
            WRITE(*,*) '      3 TO CHANGE ALL CONSTANTS'
            if(iflush) call gflush
            READ(*,*) ICH
            WRITE(99,*) ICH
            IF(ICH.EQ.1 .OR. ICH.EQ.2) THEN
              IF(ICH.EQ.1) THEN
                WRITE(*,*) ' ELA AT ',AMASS(4),' CHANGE BY?'
                if(iflush) call gflush
                READ(*,*) CHNGELA
                WRITE(99,*) CHNGELA
                BALGRAD=AMASS(1)/AMASS(4)
                DO I=4,6
                  AMASS(I)=AMASS(I)+CHNGELA
                ENDDO
                AMASS(1)=BALGRAD*AMASS(4)
              ELSEIF(ICH.EQ.2) THEN
                WRITE(*,*) 'PEAK AND DOME',AMASS(2),AMASS(3)
                if(iflush) call gflush
                READ(*,*) AMASS(2),AMASS(3)
                WRITE(99,*) AMASS(2),AMASS(3)
              ENDIF
              WRITE(*,331) (AMASS(I),I=1,6)
              DO I=1,NUMNP
                ADOT(I)=AFUNCT(TIME, IDT(I), AMASS, HTICE(I),
     &                          BDROCK(I), PSURF(I), SLOPN(1,I),
     &                          XX(I),YY(I),TEMP(I))
              ENDDO
              goto 2100
            ELSE
              if(iflush) call gflush
              READ(*,*) (AMASS(I),I=1,6)
              WRITE(99,*) (AMASS(I),I=1,6)
            ENDIF
        ENDIF
        IF(IDOT.EQ.0) THEN
          WRITE(*,*) 'INPUT -1 TO COPY, 0 TO MULTIPLY, 1 TO ADD'
          if(iflush) call gflush
          READ(*,*) IACTION
          WRITE(99,*) IACTION
          IF(IACTION.EQ.-1) THEN
            DO I=1,ICOUNT
              IDT(IFP(I))=0
              ADOTB(IFP(I))=ADOT(IFP(I))
            ENDDO
            GOTO 2100
          ENDIF
          WRITE(*,*) 'INPUT -1 ABLATION, 0 ALL, +1 ACCUMULATION'
          WRITE(*,*) 'TO LIMIT MODIFICATION'
          if(iflush) call gflush
          READ(*,*) ISCREEN
          WRITE(99,*) ISCREEN
          IF(IACTION.EQ.0) THEN
            WRITE(*,*) 'MULTIPLICATIVE FACTOR'
            if(iflush) call gflush
            READ(*,*) AMULT
            WRITE(99,*) AMULT
          ELSE
            WRITE(*,*) 'ADDITIVE FACTOR'
            if(iflush) call gflush
            READ(*,*) AMULT
            WRITE(99,*) AMULT
          ENDIF
        ENDIF
        IF(IDOT.EQ.19) THEN
          WRITE(*,*) 'INPUT POLE SNOWLINE, SNOWLINE GRADIENT, BALANCE'
          WRITE(*,*) 'GRADIENT',(AMASS(II),II=7,9)
          if(iflush) call gflush
          READ(*,*) AMASS(7),AMASS(8),AMASS(9)
          WRITE(99,*) AMASS(7),AMASS(8),AMASS(9)
        ENDIF
      ENDIF
      DO I=1,ICOUNT
        IF(ADOTN.NE.-999.) THEN
          IF(BDROCK(IFP(I)).GT.-9999.) THEN
            ADOT(IFP(I))=ADOTN
            ADOTB(IFP(I))=ADOTN
          ENDIF
          IDT(IFP(I))=0
        ELSE
          IF(IDOT.NE.0) THEN
            ADOT(IFP(I))=AFUNCT(TIME, IDOT, AMASS, HTICE(IFP(I)),
     &                          BDROCK(IFP(I)), PSURF(IFP(I)),
     &                          SLOPN(1,IFP(I)),
     &                          XX(IFP(I)),YY(IFP(I)),TEMP(IFP(I)))
            IDT(IFP(I))=IDOT
          ELSE
            IF(IDT(IFP(I)).EQ.0) THEN
              IF(BDROCK(IFP(I)).GT.-9999.) THEN
                IF(IACTION.EQ.0) THEN
                  IF(ISCREEN.EQ.-1 .AND. ADOT(IFP(I)).LT.0.) THEN
                    ADOT(IFP(I))=AMULT*ADOT(IFP(I))
                    ADOTB(IFP(I))=ADOT(IFP(I))
                  ELSEIF(ISCREEN.EQ.1 .AND. ADOT(IFP(I)).GT.0.) THEN
                    ADOT(IFP(I))=AMULT*ADOT(IFP(I))
                    ADOTB(IFP(I))=ADOT(IFP(I))
                  ELSEIF(ISCREEN.EQ.0) THEN
                    ADOT(IFP(I))=AMULT*ADOT(IFP(I))
                    ADOTB(IFP(I))=ADOT(IFP(I))
                  ENDIF
                ELSE
                  IF(ISCREEN.EQ.-1 .AND. ADOT(IFP(I)).LT.0.) THEN
                    ADOT(IFP(I))=AMULT+ADOT(IFP(I))
                    ADOTB(IFP(I))=ADOT(IFP(I))
                  ELSEIF(ISCREEN.EQ.1 .AND. ADOT(IFP(I)).GT.0.) THEN
                    ADOT(IFP(I))=AMULT+ADOT(IFP(I))
                    ADOTB(IFP(I))=ADOT(IFP(I))
                  ELSEIF(ISCREEN.EQ.0) THEN
                    ADOT(IFP(I))=AMULT+ADOT(IFP(I))
                    ADOTB(IFP(I))=ADOT(IFP(I))
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
C IF U REMOVE THE IF(IDT.. ABOVE, TURN ON THE NEXT LINE. THIS ALLOWS
C MULTIPLUCATION OF ZONES TO CONVERT FROM ZONE TO REGULAR
c            IDT(IFP(I))=0
          ENDIF
        ENDIF
      ENDDO
      GOTO 2100
C
C
3300  CONTINUE
C TOGGLE ELEMENT DISPLAY ON OR OFF
      WRITE(*,*) ' TOGGLE ELEMENT DISPLAY ON OR OFF'
      IF(ILINE.EQ.0) THEN
        ILINE=1
        WRITE(*,*) 'NOW ON'
      ELSE
        ILINE=0
        WRITE(*,*) 'NOW OFF'
      ENDIF
      GOTO 2100
C
3400  CONTINUE
C REDRAW PICTURE WITH NEW MAIN VARIABLE
c     IF(IPASS.EQ.1) CALL DELSEG(1)
      WRITE(*,*) 'INPUT:'
      WRITE(*,*) 'A: <A>FUDGE '
      WRITE(*,*) 'B: <B>ED ELEVATION'
      WRITE(*,*) 'C: ELEMENT <C>ONST AND A<C>ON'
      WRITE(*,*) 'D: A<D>OT NET MASS BALANCE' 
      WRITE(*,*) 'E: <E>RASE'
      WRITE(*,*) 'F: <F>LOW LAW CONSTANT'
      WRITE(*,*) 'G: BASAL WATER PRESSURE <G>RADIENT'
      WRITE(*,*) 'H: <H>TICE ICE SURFACE ELEVATION'
      WRITE(*,*) 'I: D<I>FF. BETWEEN CALCULATED AND PRESENT SURFACE'
      WRITE(*,*) 'J: <J> BASAL MELT RATE' 
      WRITE(*,*) 'K: <K>ODE FOR FIXED (1) OR FREE (0) NODES' 
      WRITE(*,*) 'L: S<L>OPE OF THE SURFACE'
      WRITE(*,*) 'M: TE<M>PERATURE OF THE SURFACE'
      WRITE(*,*) 'N: S<N>OWLINE ELVATION'
      WRITE(*,*) 'O: FL<O>ATATION LINE ELEVATION'
      WRITE(*,*) 'P: <P>ERCENT OF VELOCITY DUE TO SLIDING'
      WRITE(*,*) 'R; CALVING FACTO<R>'
      WRITE(*,*) 'S: <S>LIDING LAW CONSTANT'
      WRITE(*,*) 'T: <T>HICKNESS OF ICE'
      WRITE(*,*) 'U: PS<U>RF PRESENT SURFACE ELEVATION'
      WRITE(*,*) 'V: <V>ELOCITY <V>ECTORS'
      WRITE(*,*) 'W: BASAL <W>ATER LAYER THICKNESS' 
      WRITE(*,*) 'X: GEOTHERMAL FLU<X> ' 
      WRITE(*,*) 'Y: T<Y>PE OF TEMPERATURE CALCULATION'
      WRITE(*,*) 'Z: CLIMATE <Z>ONE USED'
      WRITE(*,*) '1: <1> BOUNDARY FLUX DISPLAY'
      WRITE(*,*) '2: <2> BEDROCK DEPRESSION '
      WRITE(*,*) '4: <4> BEDROCK UPLIFT RATE '
      WRITE(*,*) '5: <5> BEDROCK DEPRESSION DIFFERENCE'
      WRITE(*,*) '6: <6> INTERNAL FLUX'
      WRITE(*,*) '7: <7> INTERNAL FLUX VELOCITY'
      WRITE(*,*) '3: <3>-D SURFACE'
      if(iflush) call gflush
      READ(*,1002) CHAR
      WRITE(99,1002) CHAR
      IF(CHAR.EQ.'k') THEN
        CHAR='K'
      ELSEIF(CHAR.EQ.'w') THEN
        CHAR='W'
      ELSEIF(CHAR.EQ.'y') THEN
        CHAR='Y'
      ELSEIF(CHAR.EQ.'a') THEN
        CHAR='A'
      ELSEIF(CHAR.EQ.'h') THEN
        CHAR='H'
      ELSEIF(CHAR.EQ.'t') THEN
        CHAR='T'
      ELSEIF(CHAR.EQ.'c') THEN
        CHAR='C'
      ELSEIF(CHAR.EQ.'d') THEN
        CHAR='D'
      ELSEIF(CHAR.EQ.'g') THEN
        CHAR='G'
      ELSEIF(CHAR.EQ.'z') THEN
        CHAR='Z'
      ELSEIF(CHAR.EQ.'u') THEN
        CHAR='U'
      ELSEIF(CHAR.EQ.'p') THEN
        CHAR='P'
      ELSEIF(CHAR.EQ.'r') THEN
        CHAR='R'
      ELSEIF(CHAR.EQ.'e') THEN
        CHAR='E'
      ELSEIF(CHAR.EQ.'f') THEN
        CHAR='F'
      ELSEIF(CHAR.EQ.'s') THEN
        CHAR='S'
      ELSEIF(CHAR.EQ.'b') THEN
        CHAR='B'
      ELSEIF(CHAR.EQ.'l') THEN
        CHAR='L'
      ELSEIF(CHAR.EQ.'i') THEN
        CHAR='I'
      ELSEIF(CHAR.EQ.'j') THEN
        CHAR='J'
      ELSEIF(CHAR.EQ.'u') THEN
        CHAR='U'
      ELSEIF(CHAR.EQ.'m') THEN
        CHAR='M'
      ELSEIF(CHAR.EQ.'n') THEN
        CHAR='N'
      ELSEIF(CHAR.EQ.'o') THEN
        CHAR='O'
      ELSEIF(CHAR.EQ.'v') THEN
        CHAR='V'
      ELSEIF(CHAR.EQ.'x') THEN
        CHAR='X'
      ENDIF
      IF(CHAR.EQ.'E') THEN
        IF(IPASS.EQ.1) CALL NEWPAG
        GOTO 2100
      ELSEIF(CHAR.EQ.'1') THEN
        CALL FLUXES(IPASS,NUMNP,NUMEL,XX,YY,HTICE,DEPB,KX,
     &              NUMGBC,IBFLUX,BFLUX)
        GOTO 2100
      !ELSEIF(CHAR.EQ.'3' .AND. .NOT.BATCH) THEN
      ELSEIF(CHAR.EQ.'3') THEN
        CALL THREED(2,MMXX,NUMNP,NUMEL,XX,YY,HTICE,BDROCK,KX,BATCH)
        GOTO 2100
      ELSEIF(CHAR.EQ.'N') THEN
        WRITE(*,*) 'SNOWLINE ELVATION (PATIENCE, THIS TAKES AWHILE... '
        DO I=1,NUMNP
          ZZ(I)=7000.d0
          AA=0.d0
          BB=4096.d0
          FA=ACCUM(IDT(I),XX(I)*.001D0,YY(I)*.001D0,AA,
     &                  SLOPN(1,I),PSURF(I),AMASS(9),TJUNK)
          FB=ACCUM(IDT(I),XX(I)*.001D0,YY(I)*.001D0,BB,
     &                  SLOPN(1,I),PSURF(I),AMASS(9),TJUNK)
          IF(FA.GT.0) THEN
            ZZ(I)=0.d0
          ELSEIF(FB.LT.0) THEN
            ZZ(I)=BB
          ELSEIF(FA*FB.GT.0) THEN
            ZZ(I)=0.d0
          ELSE
            DO K=1,13
              CC=(AA+BB)*0.5d0
              FC=ACCUM(IDT(I),XX(I)*.001D0,YY(I)*.001D0,CC,
     &                 SLOPN(1,I),PSURF(I),AMASS(9),TJUNK)
              IF(FA*FC.LT.0) THEN
                BB=CC
                FB=FC
              ELSE
                AA=CC
                FA=FC
              ENDIF
            ENDDO
            ZZ(I)=CC
          ENDIF
        ENDDO
        GOTO 2000
      ELSEIF(CHAR.EQ.'6') THEN
        WRITE(*,*) ' INTERNAL FLUX '
        DO I=1,NUMNP
          ZZ(I)=FLUX(I)*1d-9
        ENDDO
        ZMIN=0.0d0
        ZMAX=1.0D+03
        GOTO 2000
      ELSEIF(CHAR.EQ.'7') THEN
        WRITE(*,*) ' INTERNAL FLUX VELOCITY '
        DO I=1,NUMNP
          if(flux(i).gt.0 .AND. THICK(I).GT.0) then
            ZZ(I)=FLUX(I)/THICK(I)
            zz(i)=log10(zz(i))
          else
            zz(i)=0.0
          endif
        ENDDO
        ZMIN=0.0d0
        ZMAX=1.0D+04
        GOTO 2000
      ELSEIF(CHAR.EQ.'W') THEN
        WRITE(*,*) ' WTHICK (m) '
        DO I=1,NUMNP
          ZZ(I)=WTHICK(I)
        ENDDO
        ZMIN=-0.09999999D0
        ZMAX=1.30000001D0
        GOTO 2000
      ELSEIF(CHAR.EQ.'2') THEN
        WRITE(*,*) ' BED DEPRESSION: WWW '
        DO I=1,NUMNP
          ZZ(I)=WWW((I-1)*3+1)                       
        ENDDO
        ZMIN=-650.D0
        ZMAX=50.D0
        GOTO 2000
      ELSEIF(CHAR.EQ.'4') THEN
        WRITE(*,*) ' UPLIFT RATE: '
        DO I=1,NUMNP
          ZZ(I)=1000*WRATE((I-1)*3+1,1)
        ENDDO
        ZMIN=-98.D0
        ZMAX=98.D0
        GOTO 2000
      ELSEIF(CHAR.EQ.'5') THEN
        WRITE(*,*) ' DEPRESSION DIFFERENCE '
        DO I=1,NUMNP
          ZZ(I)=WDIFF(I)
        ENDDO
        ZMIN=-98.D0
        ZMAX=98.D0
        GOTO 2000
      ELSEIF(CHAR.EQ.'X') THEN
        WRITE(*,*) ' GEOTHERMAL FLUX '
        DO I=1,NUMNP
          ZZ(I)=GEOFLUX(I)
        ENDDO
        ZMIN=-0.99d0
        ZMAX=5.01d0
        GOTO 2000
      ELSEIF(CHAR.EQ.'R') THEN
        WRITE(*,*) ' CALVING FACTOR '
        DO I=1,NUMNP
          ZZ(I)=CALV(I)
        ENDDO
        ZMIN=0.0d0
        ZMAX=0.5d0
        GOTO 2000
      ELSEIF(CHAR.EQ.'J') THEN
        WRITE(*,*) ' BMELT (mm/yr)'
        DO I=1,NUMNP
          ZZ(I)=BMELT(I)*1000.d0
        ENDDO
        ZMIN=-1.999d0
        ZMAX=12.001d0
        zmin=-2.d0
        zmax=12.d0
        GOTO 2000
      ELSEIF(CHAR.EQ.'K') THEN
        WRITE(*,*) 'KODE '
        DO I=1,NUMNP
          ZZ(I)=KODE(I)
        ENDDO
        GOTO 2000
      ELSEIF(CHAR.EQ.'A') THEN
        WRITE(*,*) ' AFUDGE '
        DO I=1,NUMNP
          ZZ(I)=AFUDGE(I)
        ENDDO
        ZMIN=0.0d0
        ZMAX=3.5d0
        GOTO 2000
      ELSEIF(CHAR.EQ.'O') THEN
        WRITE(*,*) ' FLOATION HEIGHT '
        DO I=1,NUMNP
          ZZ(I)=SEALEV
          IF(DEPB(I).LT.SEALEV) THEN
c            ZZ(I)=-DEPB(I)*RATDEN
            FLOT=(1.d0-RATDEN)*(DEPB(I)-SEALEV)
            ZZ(I)=FLOT
c            FLOT=(1.d0-RATDEN)*bed(I)
c            ZZ(I)=FLOT-psurf(i)
c            ZZ(I)=max(-100.d0,FLOT-psurf(i))
          ENDIF
        ENDDO
C SELF SCALING
        GOTO 2000
      ELSEIF(CHAR.EQ.'G') THEN
        WRITE(*,*) ' BASAL WATER PRESSURE POTENTIAL '
        DO I=1,NUMNP
C          ZZ(I)=RHOI*HTICE(I)+(RHOW-RHOI)*DEPB(I)
          ZZ(I)=HTICE(I)+(RHOW-RHOI)*DEPB(I)/RHOI
        ENDDO
        ZMIN=-300.0d0
        ZMAX=3900.d0
        GOTO 2000
      ELSEIF(CHAR.EQ.'H') THEN
        WRITE(*,*) 'HTICE '
        DO I=1,NUMNP
          ZZ(I)=HTICE(I)
        ENDDO
        ZMIN=-1000.0d0
        ZMAX=19000.d0
        GOTO 2000
      ELSEIF(CHAR.EQ.'T') THEN
        WRITE(*,*) 'THICKNESS'
        DO I=1,NUMNP
          ZZ(I)=HTICE(I)-DEPB(I)
          IF(DEPB(I).LT.SEALEV) THEN
c            FLOT=-DEPB(I)*RATDEN
            FLOT=(1.d0-RATDEN)*(DEPB(I)-SEALEV)
            IF(HTICE(I).LT.FLOT) THEN
              ZZ(I)=0.0
            ENDIF
          ENDIF
        ENDDO
        ZMIN=1.0d0
        ZMAX=3501.0d0
        GOTO 2000
      ELSEIF(CHAR.EQ.'D') THEN
        WRITE(*,*) 'ADOT '
        DO I=1,NUMNP
          ZZ(I)=ADOT(I)
        ENDDO
        ZMIN=-0.8d0
        ZMAX=2.0d0
        GOTO 2000
      ELSEIF(CHAR.EQ.'Z') THEN
        WRITE(*,*) 'ADOT ZONE '
        DO I=1,NUMNP
          ZZ(I)=IDT(I)
        ENDDO
        ZMIN=11.0d0
        ZMAX=25.0d0
        GOTO 2000
      ELSEIF(CHAR.EQ.'U') THEN
        WRITE(*,*) 'PSURF '
        PRINT *,'<T>HICKENSS, OR <S>URFACE'
        if(iflush) call gflush
        READ(*,1002) CHAR2
        WRITE(99,1002) CHAR2
        IF(CHAR2.EQ.'t') THEN
          CHAR2='T'
        ELSEIF(CHAR2.EQ.'s') THEN
          CHAR2='S'
        ENDIF
        IF(CHAR2.EQ.'S') THEN
          print *,' surface '
          DO I=1,NUMNP
            ZZ(I)=PSURF(I)
          ENDDO
          ZMIN=-1000.0d0
          ZMAX=19000.0d0
        ELSE
          print *,' thickness '
          DO I=1,NUMNP
            ZZ(I)=PSURF(I)-BDROCK(I)
            IF(DEPB(I).LT.SEALEV) THEN
c             FLOT=-BDROCK(I)*RATDEN
              FLOT=(1.d0-RATDEN)*(BDROCK(I)-SEALEV)
              IF(PSURF(I).LT.FLOT) THEN
                ZZ(I)=SEALEV
              ENDIF
            ENDIF
          ENDDO
          ZMIN=1.0d0
          ZMAX=3501.0d0
        ENDIF
        GOTO 2000
      ELSEIF(CHAR.EQ.'M') THEN
        WRITE(*,*) 'TEMPERATURE '
        PRINT *,'<B>OTTOM'
        PRINT *,'<T>OP'
        print *,'<D>IFFERENCE'
        print *,'<O>RIGINAL or <P>ROFILES'
        if(iflush) call gflush
        READ(*,1002) CHAR2
        WRITE(99,1002) CHAR2
        IF(CHAR2.EQ.'b') THEN
          CHAR2='B'
        ELSEIF(CHAR2.EQ.'t') THEN
          CHAR2='T'
        ELSEIF(CHAR2.EQ.'d') THEN
          CHAR2='D'
        ELSEIF(CHAR2.EQ.'o') THEN
          CHAR2='O'
        ELSEIF(CHAR2.EQ.'p') THEN
          CHAR2='P'
        ENDIF
        IF(CHAR2.EQ.'T') THEN
          DO I=1,NUMNP
            ZZ(I)=TEMP(I)
          ENDDO
        ELSEIF(CHAR2.EQ.'B') THEN
          DO I=1,NUMNP
            ZZ(I)=TBED(I)
          ENDDO
        ELSEIF(CHAR2.EQ.'D') THEN
          DO I=1,NUMNP
            ZZ(I)=-(TBED(I)-TEMP(I))
          ENDDO
        ELSEIF(CHAR2.EQ.'O') THEN
          DO I=1,NUMNP
            ZZ(I)=TSORIG(I)
          ENDDO
        ELSEIF(CHAR2.EQ.'P') THEN
        CALL TPROF(NUMNP,HTICE,DEPB,.true.,0,TIME,
     &             BMELT,WTHICK)
          goto 2100
        ENDIF
        ZMIN=-56.0d0
        ZMAX=0.0d0
        GOTO 2000
      ELSEIF(CHAR.EQ.'P') THEN
        WRITE(*,*) ' FRACT '
        DO I=1,NUMNP
           ZZ(I)=FRACT(I)
        ENDDO
        ZMIN=0.0d0
        ZMAX=1.0d0
        GOTO 2000
      ELSEIF(CHAR.EQ.'Y') THEN
        WRITE(*,*) ' ITYPE '
        DO I=1,NUMNP
           ZZ(I)=ITYPE(I)
        ENDDO
C SELF SCALING
        GOTO 2000
      ELSEIF(CHAR.EQ.'F') THEN
        WRITE(*,*) 'FLOWA '
        DO I=1,NUMNP
          ZZ(I)=FLOWA(I)
        ENDDO
        ZMIN=0.0d0
        ZMAX=7.0d0
        GOTO 2000
      ELSEIF(CHAR.EQ.'S') THEN
        WRITE(*,*) 'SLDGB '
        DO I=1,NUMNP
          ZZ(I)=SLDGB(I)
        ENDDO
        ZMIN=0.0d0
        ZMAX=0.1d0
        GOTO 2000
      ELSEIF(CHAR.EQ.'B') THEN
        WRITE(*,*) 'BED: <U>NDEPRESSED, <D>EPRESSED, OR INPUT <B>ED? '
        if(iflush) call gflush
        READ(*,1002) CHAR2
        WRITE(99,1002) CHAR2
        IF(CHAR2.EQ.'d') THEN
          CHAR2='D'
        ELSEIF(CHAR2.EQ.'u') THEN
          CHAR2='U'
        ELSEIF(CHAR2.EQ.'b') THEN
          CHAR2='B'
        ENDIF
        IF(CHAR2.EQ.'U') THEN
          DO I=1,NUMNP
            ZZ(I)=BDROCK(I)
            ZZ(I)=UNDEPB(I)
          ENDDO
        ELSEIF(CHAR2.EQ.'D') THEN
          DO I=1,NUMNP
            ZZ(I)=DEPB(I)
          ENDDO
        ELSEIF(CHAR2.EQ.'B') THEN
          DO I=1,NUMNP
            ZZ(I)=BDROCK(I)
          ENDDO
        ENDIF
        ZMIN=-2000.0d0
        ZMAX=3600.0d0
        GOTO 2000
      ELSEIF(CHAR.EQ.'L') THEN
        WRITE(*,*) 'RAW SLOPE '
        DO I=1,NUMNP
c         following is wind/slope
c         ZZ(I)=SLOPN(4,I)/slopn(1,i)
c         following is raw slope
          ZZ(I)=SLOPN(1,I)
        ENDDO
C SELF SCALING
        GOTO 2000
      ELSEIF(CHAR.EQ.'I') THEN
        WRITE(*,*) 'DIFFERENCE '
        DO I=1,NUMNP
C         ZZ(I)=HTICE(I)-HFIT(I)
          ZZ(I)=HTICE(I)-PSURF(I)
        ENDDO
        ZMIN=-1300.0d0
        ZMAX=1500.0d0
        GOTO 2000
      ELSEIF(CHAR.EQ.'U') THEN
        WRITE(*,*) 'PRESENT SURFACE '
        DO I=1,NUMNP
          ZZ(I)=PSURF(I)
        ENDDO
        ZMIN=-1000.0d0
        ZMAX=19000.0d0
        GOTO 2000

      ELSEIF(CHAR.EQ.'C') THEN
        PRINT *,'CONSTANT AND ACON FOR ELEMENTS'
        PRINT *,'C:<C>ONST                 (SCALED TO FIT)'
        PRINT *,'V:<V>ELOCITY IN ACON SLOT (SCALED TO FIT)'
        PRINT *,'A:<A>CON                  (NOT SCALED)'
        if(iflush) call gflush
        READ(*,1002) CHAR2
        WRITE(99,1002) CHAR2
        IF(CHAR2.EQ.'c') CHAR2='C'
        IF(CHAR2.EQ.'a') CHAR2='A'
        IF(CHAR2.EQ.'v') CHAR2='V'
        IF(IPASS.EQ.0) THEN
          CALL GRSTRT(600,600)
          IPASS=1
        ENDIF
        XMIN=1.d30
        YMIN=XMIN
        XMAX=-XMIN
        YMAX=-YMIN
        DO I=1,NUMNP
          XMAX=MAX(XMAX,XX(I))
          YMAX=MAX(YMAX,YY(I))
          XMIN=MIN(XMIN,XX(I))
          YMIN=MIN(YMIN,YY(I))
        ENDDO
        XWID=XMAX-XMIN
        YWID=YMAX-YMIN
        IF(XWID.GT.YWID) THEN
          YWID=XWID
          YMAX=YMIN+YWID
        ELSE
          XWID=YWID
          XMAX=XMIN+XWID
        ENDIF
        XMIN=XMIN-XWID*.05d0
        XMAX=XMAX+XWID*.25d0
        YMIN=YMIN-YWID*.05d0
        YMAX=YMAX+YWID*.25d0
        VMIN=1d30
        VMAX=-VMIN
        IF(CHAR2.EQ.'C') THEN
          DO I=1,NUMEL
            VMAX=MAX(VMAX,CONST(I))
            VMIN=MIN(VMIN,CONST(I))
          ENDDO
        ELSEIF(CHAR2.EQ.'V') THEN
          DO I=1,NUMEL
            VMAX=MAX(VMAX,ACON(I))
            VMIN=MIN(VMIN,ACON(I))
          ENDDO
          PRINT *,'CURRENT MIN,MAX:',VMIN,VMAX
          if(iflush) call gflush
          READ(*,*) VMIN,VMAX
          WRITE(99,*) VMIN,VMAX
        ELSEIF(CHAR2.EQ.'A') THEN
          VMIN=0.d0
          VMAX=7.d0
        ENDIF
        NCOLOR=15
        VDELT=(VMAX-VMIN)/(NCOLOR-1)
c        CALL NEWPAG
        CALL WINDOW(REAL(XMIN),REAL(XMAX),REAL(YMIN),REAL(YMAX))
        XPOS1=XMAX-XWID*.25d0
        XPOS2=XMAX-XWID*.15d0
        YPOS=YMIN+YWID*.15d0
        YDELT=YWID/(NCOLOR-1)
        DO I=1,NCOLOR
          CALL LINCLR(ICMAP1(I))
          CALL MOVE(REAL(XPOS1),REAL(YPOS-YDELT/2))
          CALL DRAW(REAL(XPOS2),REAL(YPOS-YDELT/2))
          CALL MOVE(REAL(XPOS2),REAL(YPOS))
          RNUM=VMIN+(I-1)*VDELT
          CALL LINCLR(1)
          CALL GNUMBR(REAL(RNUM),2,6)
          YPOS=YPOS+YDELT
        ENDDO
        DO I=1,NUMEL
          IF(CHAR2.EQ.'C') THEN
            ZPLOT=CONST(I)
          ELSE
            ZPLOT=ACON(I)
          ENDIF
          RZ=(ZPLOT-VMIN)/VDELT
          ICOLOR=int(RZ+2)
          IF(ICOLOR.LT.2) ICOLOR=2
          IF(ICOLOR.GT.16) ICOLOR=16
          PRINT *,ZPLOT,ICOLOR
          DO L=1,4
            PX(L)=real(XX(KX(I,L)))
            PY(L)=real(YY(KX(I,L)))
          ENDDO
          CALL FILPAN(ICMAP1(ICOLOR),.FALSE.)
          CALL PANEL(4,PX,PY)
        ENDDO
        GOTO 2100
      
      ELSEIF(CHAR.EQ.'V') THEN
        PRINT *,'VELOCITY w/ CONTOURS or MAP/PROFILES'
        WRITE(*,*) 'INPUT I FOR <I>CE VELOCITIES, W FOR <W>ATER'
        WRITE(*,*) '     OR <M> FOR MAP/PROFILES...'
        READ(*,1002) CHAR4
        WRITE(99,1002) CHAR4
        IF(CHAR4.EQ.'w') THEN
          CHAR4='W'
        ENDIF
        IF(CHAR4.EQ.'i') THEN
          CHAR4='I'
        ENDIF
        IF(CHAR4.EQ.'m') THEN
          CHAR4='M'
        ENDIF
        IF(CHAR4.NE.'M') THEN
          IF(IPASS.EQ.0) THEN
            CALL GRSTRT(600,600)
            IPASS=1
          ENDIF
          CALL SCALEXY(xx,yy,NUMNP,xmin,xmax,ymin,ymax,delx,dely)
c         CALL NEWPAG
          CALL WINDOW(REAL(XMIN-DELX),REAL(XMAX+DELX),
     &              REAL(YMIN-DELY),REAL(YMAX+DELY))
          ZMAX=4200.00001D0
          ZMIN=0.00001D0
          DELZ=(ZMAX-ZMIN)/14.d0
          CALL CONTR(MMXX,NUMEL,XX,YY,KX,
     &                 HTICE,ZMIN,ZMAX,DELZ,
     &                 XMIN,XMAX,YMIN,YMAX)
          WRITE(*,*) 'INPUT SCALE FACTOR '
          WRITE(*,*) 'FOR ARROWS, UMAX, AND THRESHOLD'
          WRITE(*,*) 'INPUT 1 FOR NORMALIZED VECTORS, 0 FOR ABSOLUTE'  
          WRITE(*,*) ASCAL,UMAX,VTHRESH,INORM                                       
          if(iflush) call gflush
          READ(*,*) ASCAL,UMAX,VTHRESH,INORM                                       
          WRITE(99,*) ASCAL,UMAX,VTHRESH,INORM                  
        endif                     
        IF(CHAR4.EQ.'W') THEN
          CALL WVELO(XX,YY, KX, NUMNP, NUMEL,
     &               WTHICK, HTICE, DEPB, .true.)
        ELSEIF(CHAR4.EQ.'I') THEN
          CALL DVELO(MMXX,NUMEL,KX,XX,YY,HTICE,DEPB,CONST,
     &               VEL,.TRUE.,DTMIN)
        ELSEIF(CHAR4.EQ.'M') THEN
          CALL DVELO(MMXX,NUMEL,KX,XX,YY,HTICE,DEPB,CONST,
     &               VEL,.false.,DTMIN)
          DO I=1,NUMNP
            ZZ(I)=0
            TMP(I)=0
          ENDDO
          DO I=1,NUMEL
            DO J=1,4
              ZZ(KX(I,J))=ZZ(KX(I,J))+VEL(I,1)
              TMP(KX(I,J))=TMP(KX(I,J))+1
            ENDDO
          ENDDO
          DO I=1,NUMNP
            IF(TMP(I).EQ.0) THEN
              ZZ(I)=0
            ELSE
              ZZ(I)=ZZ(I)/TMP(I)
            ENDIF
          ENDDO
          ZMIN=0.0D0
          ZMAX=1.0D+03
          GOTO 2000
        ELSE
          GOTO 2100
        ENDIF
        CALL DCOAST(1)
        CALL LINCLR(1)                                                    
        CALL MOVE(REAL(XMIN),REAL(YMIN))                                  
        CALL DRAW(REAL(XMIN),REAL(YMAX))                                  
        CALL DRAW(REAL(XMAX),REAL(YMAX))                                  
        CALL DRAW(REAL(XMAX),REAL(YMIN))                                  
        CALL DRAW(REAL(XMIN),REAL(YMIN))                                  
        GOTO 2100
      ELSE
        CALL SETBEL(2)
        CALL RINGBE
        GOTO 2100
      ENDIF
C
C
C
3500  CONTINUE
C CHANGE XY COORDINATE
      print *,'not active'
      IF(IPASS.EQ.0 .or. .true.) THEN
        WRITE(*,*) 'MUST HAVE GRAPH SHOWING'
c        CALL SETBEL(2)
c        CALL RINGBE
        GOTO 2100
      ENDIF
      IF(BATCH) THEN
        DO M=1,1
          READ(*,*) XA(M),YA(M),IDAT(M)
        ENDDO
        IGOT=1
      ELSE
        CALL LOCATE(1,XA,YA,IDAT,IGOT)
      ENDIF
      DO M=1,1
        WRITE(99,*) XA(M),YA(M),IDAT(M)
      ENDDO
      ICHAR=IDAT(1)
      xxx=dble(xa(1))*1.d3
      yyy=dble(ya(1))*1.d3
      DMIN=1.d30
      DO I=1,NUMNP
        DIST=(XXX-XX(I))**2+(YYY-YY(I))**2
        IF(DIST.LT.DMIN) THEN
          DMIN=DIST
          IFOUND=I
        ENDIF
      ENDDO
      xxx=xx(ifound)*1.d-3
      yyy=yy(ifound)*1.d-3
      CALL MRKCLR(1)
C     CALL MARKER(XXX,YYY,9)
      CALL MARKER(REAL(XXX),REAL(YYY),0)
      IF(BATCH) THEN
        DO M=1,1
          READ(*,*) XA(M),YA(M),IDAT(M)
        ENDDO
        IGOT=1
      ELSE
        CALL LOCATE(1,XA,YA,IDAT,IGOT)
      ENDIF
      DO M=1,1
        WRITE(99,*) XA(M),YA(M),IDAT(M)
      ENDDO
      ICHAR=IDAT(1)
      CALL MRKCLR(2)
      CALL MARKER(REAL(XA(1)),REAL(YA(1)),10)
      CALL MRKCLR(0)
C     CALL MARKER(REAL(XXX),REAL(YYY),9)
      CALL MARKER(REAL(XXX),REAL(YYY),0)
      CALL MRKCLR(1)
      xxx=dble(xa(1))*1.d3
      yyy=dble(ya(1))*1.d3
      XX(IFOUND)=XXX
      YY(IFOUND)=YYY
      GOTO 2100
C
3600  CONTINUE
C ADD NEW NODE
      print *,'not active'
      IF(IPASS.EQ.0 .or. .true.) THEN
        WRITE(*,*) 'MUST HAVE GRAPH SHOWING'
        GOTO 2100
      ENDIF
      NUMNP=NUMNP+1
      IF(BATCH) THEN
        DO M=1,1
          READ(*,*) XA(M),YA(M),IDAT(M)
        ENDDO
        IGOT=1
      ELSE
        CALL LOCATE(1,XA,YA,IDAT,IGOT)
      ENDIF
      DO M=1,1
        WRITE(99,*) XA(M),YA(M),IDAT(M)
      ENDDO
      ICHAR=IDAT(1)
C     CALL MARKER(XA(1),YA(1),9)
      CALL MARKER(REAL(XA(1)),REAL(YA(1)),0)
      xx(numnp)=dble(xa(1))*1.d3
      yy(numnp)=dble(ya(1))*1.d3
      KODE(NUMNP)=0
      HTICE(NUMNP)=100.d0
      ADOT(NUMNP)=1.0d0
      FRACT(NUMNP)=0.0d0
      PSURF(NUMNP)=0.0d0
      BDROCK(NUMNP)=0.0d0
      FLOWA(NUMNP)=2.0d0
      SLDGB(NUMNP)=.02d0
      GOTO 2100
C
3700  CONTINUE
C DEFINE NEW ELEMENT CONNECTIVITY
      print *,'not active'
      IF(IPASS.EQ.0 .or. .true.) THEN
        WRITE(*,*) 'MUST HAVE GRAPH SHOWING'
        GOTO 2100
      ENDIF
      NUMEL=NUMEL+1
      WRITE(*,*) 'INPUT NUMBER OF NODES FOR NEW ELEMENT ',NUMEL
      if(iflush) call gflush
      READ(*,*) NN
      WRITE(99,*) NN
      print *,'NOT WORKING'
      if(.false.) then
        NNODE(NUMEL)=NN
        CONST(NUMEL)=1.d6
        KX(NUMEL,4)=0
        WRITE(*,*) 'PICK ENDPOINTS OF NEW ELEMENT BOUNDARY '
        CALL LOCATE(NN+1,XA,YA,IDAT,IGOT)
        ICHAR=IDAT(NN+1)
        DO I=1,NN
          xxx=dble(xa(1))*1.d3
          yyy=dble(ya(1))*1.d3
          XCHECK(I)=XXX
          YCHECK(I)=YYY
        ENDDO
        DO I=1,NN
          DMIN=1.d30
          KX(NUMEL,I)=0
          DO JJ=1,NUMNP
            DIST=(XCHECK(I)-XX(JJ))**2+(YCHECK(I)-YY(JJ))**2
            IF(DIST.LT.DMIN) THEN
              DMIN=DIST
              KX(NUMEL,I)=JJ
            ENDIF
          ENDDO
        ENDDO
      endif
      GOTO 2100
C
3800  CONTINUE
      IFOUND=0
      DO IC=1,ICOUNT
        DO 3810 I=1,NUMEL
          IF(KX(I,4).EQ.0) THEN
            NN=3
          ELSE
            NN=4
          ENDIF
          DO N=1,NN
            IF(KX(I,N).EQ.IFP(IC)) THEN
C             FOUND ELEMENT CONTAINING NODE, DELETE
              IFOUND=IFOUND+1
              IFD(IFOUND)=I
              GOTO 3810
            ENDIF
          ENDDO
3810    CONTINUE
        xxx=xx(IFP(IC))*1.d-3
        yyy=yy(IFP(IC))*1.d-3
        CALL MRKCLR(2)
        CALL point(REAL(XXX),REAL(YYY))
c        CALL MARKER(REAL(XXX),REAL(YYY),10)
        if(iflush) call gflush
      ENDDO
      print *,'sorting'
      CALL SORTIX(IFOUND,IFD)
      CALL ELMDUP(IFOUND,IFD)
      print *,'sorting done'
      DO  I=1,IFOUND
        DO K=1,4
c          xxx=xx(IFD(IC))*1.d-3
c          yyy=yy(IFD(IC))*1.d-3
          xxx=xx(IFD(I))*1.d-3
          yyy=yy(IFD(I))*1.d-3
          kode(KX(IFD(I),K))=1
          CALL MRKCLR(3)
          CALL point(REAL(XXX),REAL(YYY))
c          CALL MARKER(REAL(XXX),REAL(YYY),10)
          if(iflush) call gflush
        ENDDO
        DO J=IFD(I),NUMEL-1
          DO K=1,4
            KX(J,K)=KX(J+1,K)
            CONST(J)=CONST(J+1)
          ENDDO
        ENDDO
        NUMEL=NUMEL-1
      ENDDO
      print *,'sorting again'
      CALL SORTIX(ICOUNT,IFP)
      print *,'sorting done'
      DO J=1,ICOUNT
        xxx=xx(IFP(j))*1.d-3
        yyy=yy(IFP(j))*1.d-3
        CALL MRKCLR(4)
        CALL point(REAL(XXX),REAL(YYY))
c        CALL MARKER(REAL(XXX),REAL(YYY),10)
        if(iflush) call gflush
        DO N=1,NUMEL
          DO K=1,4
            IF(KX(N,K).GE.IFP(J)) KX(N,K)=KX(N,K)-1
          ENDDO
        ENDDO
      ENDDO
      DO I=1,ICOUNT
        xxx=xx(IFP(i))*1.d-3
        yyy=yy(IFP(i))*1.d-3
        CALL MRKCLR(5)
        CALL point(REAL(XXX),REAL(YYY))
c        CALL MARKER(REAL(XXX),REAL(YYY),10)
        if(iflush) call gflush
        DO J=IFP(I),NUMNP-1
          XX(J)=XX(J+1)
          YY(J)=YY(J+1)
          HTICE(J)=HTICE(J+1)
          ADOT(J)=ADOT(J+1)
          ADOTB(J)=ADOTB(J+1)
          BDROCK(J)=BDROCK(J+1)
          FRACT(J)=FRACT(J+1)
          PSURF(J)=PSURF(J+1)
          FLOWA(J)=FLOWA(J+1)
          SLDGB(J)=SLDGB(J+1)
          THICK(J)=THICK(J+1)
          CALV(J)=CALV(J+1)
          GEOFLUX(J)=GEOFLUX(J+1)
          ITYPE(J)=ITYPE(J+1)
          KODE(J)=KODE(J+1)
          IDT(J)=IDT(J+1)
          TEMP(J)=TEMP(J+1)
        ENDDO
        NUMNP=NUMNP-1
      ENDDO
      if(LDANGLE) then
        LDANGLE=.false.
        print *,'removing danglers...'
        ICOUNT=0
        DO I=1,NUMNP
          if(mod(i,1000).eq.0) print *,'checking:',i,numnp
          IFOUND=0
          DO N=1,NUMEL
            DO K=1,4
              IF(KX(N,K).eq.i) ifound=ifound+1
            ENDDO
            if(ifound.ne.0) goto 3790
          ENDDO
          IF(IFOUND.eq.0) THEN
            print *,'removing ',I,ICOUNT+1
            ICOUNT=ICOUNT+1
            IFP(ICOUNT)=I
            xxx=xx(I)*1.d-3
            yyy=yy(I)*1.d-3
            CALL MRKCLR(1)
            CALL point(REAL(XXX),REAL(YYY))
c           CALL MARKER(REAL(XXX),REAL(YYY),10)
          ENDIF
3790      continue
        ENDDO
        if(icount.gt.0) goto 3800
      endif
      GOTO 2100
C
3900  CONTINUE
C QUIT
      GOTO 999
C
4000  CONTINUE
C BACKUP
      CALL WRITPROF(NUMNP,XX,YY,HTICE,BDROCK,ADOT)
      CALL WRITEMP(NUMNP)
      CALL WRITEH2O(NUMNP,WTHICK,BMELT)
      if(BTOGG.eq.1) then
c ..... use following for new point load visco-elastic plate backup .........
        CALL WRITEDEPS(NUMNP,3*NUMNP)
      elseif(BTOGG.eq.3) then
c ..... use following for new visco-elastic plate backup .........
        CALL WRITEDEPN(NUMNP,3*NUMNP)
      elseif(BTOGG.eq.2) then
c ..... use following for old elastic plate backup .........
        CALL WRITEDEPO(NUMNP,3*NUMNP)
      endif
      REWIND 8
      AMASS(11)=SEALEV
      WJUNK=WINDIR(1)*180.d0/3.14159d0    
      WRITE(8,*) AMASS,ACOM,WJUNK,WINDIR(2),
     &             XPOLE,YPOLE,CFACTOR,ALPHAC,TBASE,
     &             TNSLBASE,CTOGG,WTOGG,ITOGG,BTOGG
C
      REWIND 21
100   FORMAT(A80,/,7I6,F8.0)
200   FORMAT(I6,I4,1P,2E12.5,0P,F10.2,F7.2,F9.2,F10.3,F10.1,
     &          F10.3,2F10.5,I5,F10.3,F10.3,1PE10.3)
300   FORMAT(5I6,1P2E17.10)
310   FORMAT(2I6,1PE13.6)
      WRITE(21,100) HED, NUMNP, NUMEL, NUMCOL, NUMLEV, NUMGBC, NDT,
     &              INTER, DT
      DO NUM=1,NUMNP
        IF(IDT(NUM).EQ.0) THEN
          IF(BTOGG.eq.2 .or. BTOGG.eq.3) THEN
c         WRITE(21,200) NUM, KODE(NUM), XX(NUM), YY(NUM), HTICE(NUM),
          WRITE(21,*) NUM, KODE(NUM), XX(NUM), YY(NUM), HTICE(NUM),
c    &                  ADOTB(NUM),FRACT(NUM),PSURF(NUM),DEPB(NUM),
     &                  ACC(NUM),FRACT(NUM),PSURF(NUM),DEPB(NUM),
c    &                  FLOWA(NUM), SLDGB(NUM), TEMP(NUM),ITYPE(NUM),
     &                  3.0, SLDGB(NUM), TEMP(NUM),ITYPE(NUM),
c    &                  AFUDGE(NUM),GEOFLUX(NUM),CALV(NUM)
     &                  AFUDGE(NUM),GEOFLUX(NUM),ABLAT(NUM),0.,0.
          ELSE
c         WRITE(21,200) NUM, KODE(NUM), XX(NUM), YY(NUM), HTICE(NUM),
          WRITE(21,*) NUM, KODE(NUM), XX(NUM), YY(NUM), HTICE(NUM),
c    &                  ADOTB(NUM),FRACT(NUM),PSURF(NUM),BDROCK(NUM),
     &                  ACC(NUM),FRACT(NUM),PSURF(NUM),BDROCK(NUM),
c    &                  FLOWA(NUM), SLDGB(NUM), TEMP(NUM),ITYPE(NUM),
     &                  3.0, SLDGB(NUM), TEMP(NUM),ITYPE(NUM),
c    &                  AFUDGE(NUM),GEOFLUX(NUM),CALV(NUM)
     &                  AFUDGE(NUM),GEOFLUX(NUM),ABLAT(NUM),0.,0.
          ENDIF
        ELSE
          ADDT=-100-IDT(NUM)
          IF(BTOGG.eq.2 .or. BTOGG.eq.3) THEN
c         WRITE(21,200) NUM, KODE(NUM), XX(NUM), YY(NUM), HTICE(NUM),
          WRITE(21,*) NUM, KODE(NUM), XX(NUM), YY(NUM), HTICE(NUM),
     &                  ADDT, FRACT(NUM), PSURF(NUM), DEPB(NUM),
c    &                  FLOWA(NUM), SLDGB(NUM), TEMP(NUM),ITYPE(NUM),
     &                  3.0, SLDGB(NUM), TEMP(NUM),ITYPE(NUM),
     &                  AFUDGE(NUM),GEOFLUX(NUM),CALV(NUM),0.,0.
          ELSE
c         WRITE(21,200) NUM, KODE(NUM), XX(NUM), YY(NUM), HTICE(NUM),
          WRITE(21,*) NUM, KODE(NUM), XX(NUM), YY(NUM), HTICE(NUM),
     &                  ADDT, FRACT(NUM), PSURF(NUM), BDROCK(NUM),
c    &                  FLOWA(NUM), SLDGB(NUM), TEMP(NUM),ITYPE(NUM),
     &                  3.0, SLDGB(NUM), TEMP(NUM),ITYPE(NUM),
     &                  AFUDGE(NUM),GEOFLUX(NUM),CALV(NUM),0.,0.
          ENDIF
        ENDIF
C       THICK(NUM)=HTICE(NUM)-BDROCK(NUM)
        THICK(NUM)=HTICE(NUM)-DEPB(NUM)
        ZZ(NUM)=HTICE(NUM)
      ENDDO
      WRITE(*,*) 'OUTPUT VELOCITIES:(1) OR FLOW CONST:(0)'
      if(iflush) call gflush
      READ(*,*) IOUT
      WRITE(99,*) IOUT
      IF(IOUT.EQ.1) THEN
        CALL DVELO(MMXX,NUMEL,KX,XX,YY,HTICE,DEPB,CONST,
     &             VEL,.FALSE.,DTMIN)
        DO I=1,NUMEL
c         WRITE(21,300) I, (KX(I,II),II=1,4),CONST(I),VEL(I,1)
          WRITE(21,*) I, (KX(I,II),II=1,4),CONST(I),VEL(I,1)
        ENDDO
      ELSE
        DO I=1,NUMEL
c         WRITE(21,300) I, (KX(I,II),II=1,4),CONST(I),ACON(I)
          WRITE(21,*) I, (KX(I,II),II=1,4),CONST(I),ACON(I)
        ENDDO
      ENDIF
      DO N=1,NUMGBC
c       WRITE(21,310) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
        WRITE(21,*) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
      ENDDO
      GOTO 2100
C
C DEFINE LINE
4100  CONTINUE
      WRITE(*,*) '<D>RAW,<2>-POINT OR <N>-POINT LINE?'
      if(iflush) call gflush
      READ(*,1002) CHAR2
      WRITE(99,1002) CHAR2
      IF(CHAR2.EQ.'n') CHAR2='N'
      IF(CHAR2.EQ.'d') CHAR2='D'
      IF(CHAR2.EQ.'D') THEN
        NP=0
        inquire(file='fort.20',exist=file20)
        if(file20) then
          REWIND 20
          READ(20,1010,END=103) NP
          READ(20,1010) (NLINE(IP),IP=1,NP)
        endif
103     CONTINUE
        CALL MRKCLR(1)
        DO I=1,NP
          xxx=xx(NLINE(i))*1.d-3
          yyy=yy(NLINE(i))*1.d-3
          CALL MARKER(REAL(XXX),REAL(YYY),10)
        ENDDO
        GOTO 2100
      ELSEIF(CHAR2.EQ.'N') THEN
        NP=0
4105    CONTINUE
        IF(IPASS.EQ.0) THEN
          WRITE(*,*) 'MUST HAVE GRAPH SHOWING'
          GOTO 2100
        ENDIF
        IF(BATCH) THEN
          DO M=1,1
            READ(*,*) XA(M),YA(M),IDAT(M)
          ENDDO
          IGOT=1
        ELSE
          CALL LOCATE(1,XA,YA,IDAT,IGOT)
        ENDIF
        DO M=1,1
          WRITE(99,*) XA(M),YA(M),IDAT(M)
        ENDDO
        ICHAR=IDAT(1)
        IF(ICHAR.EQ.113 .OR. ICHAR.EQ.81) GOTO 4120
        xxx=dble(xa(1))*1.d3
        yyy=dble(ya(1))*1.d3
        DMIN=1.d30
        DO I=1,NUMNP
          DIST=(XXX-XX(I))**2+(YYY-YY(I))**2
          IF(DIST.LT.DMIN) THEN
            DMIN=DIST
            IFOUND=I
          ENDIF
        ENDDO
        CALL MRKCLR(1)
        xxx=xx(ifound)*1.d-3
        yyy=yy(ifound)*1.d-3
        CALL MARKER(REAL(XXX),REAL(YYY),10)
        NLINE(NP+1)=IFOUND
        NP=NP+1
        GOTO 4105
      ELSEIF(CHAR2.EQ.'2') THEN
        IF(IPASS.EQ.0) THEN
          WRITE(*,*) 'MUST HAVE GRAPH SHOWING'
          GOTO 2100
        ENDIF
        IF(BATCH) THEN
          DO M=1,2
            READ(*,*) XA(M),YA(M),IDAT(M)
          ENDDO
          IGOT=2
        ELSE
          CALL LOCATE(2,XA,YA,IDAT,IGOT)
        ENDIF
        DO M=1,2
          WRITE(99,*) XA(M),YA(M),IDAT(M)
        ENDDO
        xx1=dble(xa(1))*1.d3
        yy1=dble(ya(1))*1.d3
        xx2=dble(xa(2))*1.d3
        yy2=dble(ya(2))*1.d3
        DX=(XX2-XX1)/(100-1)
        DY=(YY2-YY1)/(100-1)
        NNP=0
        DO I=1,100
          XXC=XX1+(I-1)*DX
          YYC=YY1+(I-1)*DY
          DISTMIN=1d30
          N=0
          DO J=1,NUMNP
            DIST=SQRT((XX(J)-XXC)**2+(YY(J)-YYC)**2)
            IF(DIST.LT.DISTMIN) THEN
              N=J
              DISTMIN=DIST
            ENDIF
          ENDDO
          IF(N.EQ.0) PRINT *,'PROBLEMS...'
          NNP=NNP+1
          NNLINE(NNP)=N
          CALL MRKCLR(1)
          xxx=xx(n)*1.d-3
          yyy=yy(n)*1.d-3
          CALL MARKER(REAL(XXX),REAL(YYY),0)
        ENDDO
        NP=1
        NLINE(1)=NNLINE(1)
        DO I=2,NNP
          IF(NNLINE(I).NE.NLINE(NP)) THEN
            NP=NP+1
            NLINE(NP)=NNLINE(I)
          ENDIF
        ENDDO
        CALL MRKCLR(1)
        DO I=1,NP
          xxx=xx(NLINE(I))*1.d-3
          yyy=yy(NLINE(I))*1.d-3
          CALL MARKER(REAL(XXX),REAL(YYY),0)
        ENDDO
        REWIND 20
      ENDIF
4120  CONTINUE
      WRITE(20,1010) NP
      WRITE(20,1010) (NLINE(IP),IP=1,NP)
1010  FORMAT(13I6)
      GOTO 2100
C
4200  CONTINUE
C FOLLOWING IS REFERENCE SUURFACE FOR ADJUSTMENTS
      RSURF=-500.d0
C CHANGE VARIOUS PHYSICAL PARAMETERS
      WRITE(*,*) ' MODIFY VARIOUS PHYSICAL PARAMETERS'
      WRITE(*,*) 'A<D>OT, <P>ERCENT, <S>LIDING, <F>LOW'
      if(iflush) call gflush
      READ(*,1002) CHAR
      WRITE(99,1002) CHAR
      IF(CHAR.EQ.'d') THEN
        CHAR='D'
      ELSEIF(CHAR.EQ.'p') THEN
        CHAR='P'
      ELSEIF(CHAR.EQ.'s') THEN
        CHAR='S'
      ELSEIF(CHAR.EQ.'f') THEN
        CHAR='F'
      ENDIF
      WRITE(*,*) 'INPUT RATIO FACTOR'
      if(iflush) call gflush
      READ(*,*) RFACT
      WRITE(99,*) RFACT
      IF(CHAR.EQ.'D') GOTO 4210
      IF(CHAR.EQ.'P') GOTO 4210
      IF(CHAR.EQ.'F') GOTO 4210
      IF(CHAR.EQ.'S') GOTO 4210
      GOTO 4200
C
4210  CONTINUE
      DO I=1,ICOUNT
        J=IFP(I)
        IF(KODE(J).EQ.1) THEN
          RATIOD=0.d0
          RATIOF=1.d0
        ELSE
          RATIOD=1.d0-(HTICE(J)-RSURF)/(PSURF(J)-RSURF)
          RATIOD=RATIOD*RFACT
          RATIOF=RFACT*(PSURF(J)-RSURF)/(HTICE(J)-RSURF)
        ENDIF
        IF(CHAR.EQ.'D') THEN
          IF(CHAR.EQ.'D') THEN
            ADOT(J)=ADOT(J)+RATIOD
            ADOTB(J)=ADOTB(J)+RATIOD
          ENDIF
        ELSEIF(CHAR.EQ.'P') THEN
c          FRACT(J)=FRACT(J)-FRACT(J)*RATIOD
          FRACT(J)=FRACT(J)-0.1d0*RATIOD
          IF(FRACT(J).GT.1.) FRACT(J)=1.d0
          IF(FRACT(J).LT.0.) FRACT(J)=0.d0
        ELSEIF(CHAR.EQ.'F') THEN
          FLOWA(J)=FLOWA(J)*RATIOF
          IF(FLOWA(J).GT.10.) FLOWA(J)=10.d0
          IF(FLOWA(J).LT.FLOWMIN) FLOWA(J)=FLOWMIN
        ELSEIF(CHAR.EQ.'S') THEN
          SLDGB(J)=SLDGB(J)*RATIOF
          IF(SLDGB(J).GT..20) SLDGB(J)=.20d0
          IF(SLDGB(J).LT.SLDGBMIN) SLDGB(J)=SLDGBMIN
        ENDIF
      ENDDO
      GOTO 2100
C
4300  CONTINUE
      WRITE(*,*) ' DISPLAY FOLLOWING DURING RUN:'
      WRITE(*,*) 'F: <F>RACT'
      WRITE(*,*) 'T: <T>HICKNESS'
      WRITE(*,*) 'S: <S>URFACE'
      WRITE(*,*) 'V: <V>OLUME'
      WRITE(*,*) 'P: <P>ROFILE'
      WRITE(*,*) 'M: TE<M>PERATURES'
      WRITE(*,*) 'L: VE<L>OCITY VECTORS'
      WRITE(*,*) 'W: <W>ATER THICKNESS'
      WRITE(*,*) 'X: BOUNDARY FLU<X>ES'
      WRITE(*,*) '<N>OTHING'
      if(iflush) call gflush
      READ(*,1002) CHAR
      WRITE(99,1002) CHAR
      IF(CHAR.EQ.'t') THEN
        CHAR='T'
      ELSEIF(CHAR.EQ.'s') THEN
        CHAR='S'
      ELSEIF(CHAR.EQ.'p') THEN
        CHAR='P'
      ELSEIF(CHAR.EQ.'v') THEN
        CHAR='V'
      ELSEIF(CHAR.EQ.'n') THEN
        CHAR='N'
      ELSEIF(CHAR.EQ.'f') THEN
        CHAR='F'
      ELSEIF(CHAR.EQ.'m') THEN
        CHAR='M'
      ELSEIF(CHAR.EQ.'l') THEN
        CHAR='L'
      ELSEIF(CHAR.EQ.'w') THEN
        CHAR='W'
      ELSEIF(CHAR.EQ.'x') THEN
        CHAR='X'
      ENDIF
      IF(CHAR.EQ.'X') THEN
        IPLOT=9
      ELSEIF(CHAR.EQ.'W') THEN
        IPLOT=8
      ELSEIF(CHAR.EQ.'M') THEN
        IPLOT=7
      ELSEIF(CHAR.EQ.'L') THEN
        IPLOT=6
      ELSEIF(CHAR.EQ.'P') THEN
        IPLOT=5
      ELSEIF(CHAR.EQ.'F') THEN
        IPLOT=4
      ELSEIF(CHAR.EQ.'T') THEN
        IPLOT=3
      ELSEIF(CHAR.EQ.'S') THEN
        IPLOT=2
      ELSEIF(CHAR.EQ.'V') THEN
        IPLOT=1
      ELSEIF(CHAR.EQ.'N') THEN
        IPLOT=0
      ELSE
        IPLOT=0
      ENDIF
      GOTO 2100
4400  CONTINUE
      DO I=1,6
        IFIT(I)=0
      ENDDO
      WRITE(*,*) 'FIT ADOT 1-YES, 0-NO'
      if(iflush) call gflush
      READ(*,*) IFIT(1)
      WRITE(99,*) IFIT(1)
      WRITE(*,*) 'FIT FRACT 1-YES, 0-NO'
      if(iflush) call gflush
      READ(*,*) IFIT(2)
      WRITE(99,*) IFIT(2)
      WRITE(*,*) 'FIT FLOWA 1-YES, 0-NO'
      if(iflush) call gflush
      READ(*,*) IFIT(3)
      WRITE(99,*) IFIT(3)
      WRITE(*,*) 'FIT SLDGB 1-YES, 0-NO'
      if(iflush) call gflush
      READ(*,*) IFIT(4)
      WRITE(99,*) IFIT(4)
      WRITE(*,*) 'FIT AFUDGE 1-YES, 0-NO'
      if(iflush) call gflush
      READ(*,*) IFIT(5)
      WRITE(99,*) IFIT(5)
      WRITE(*,*) 'FIT AFUDGE/SLDGB 1-YES, 0-NO'
      if(iflush) call gflush
      READ(*,*) IFIT(6)
      WRITE(99,*) IFIT(6)
      WRITE(*,*) 'HOW OFTEN DO I APPLY MODIFY (STEPS)',IFIT(6)
      if(iflush) call gflush
      READ(*,*) IFIT(7)
      WRITE(99,*) IFIT(7)
      GOTO 2100
999   CONTINUE
c     IF(IPASS.EQ.1) CALL DELSEG(-1)
      ISKIP=0
      IF(ISKIP.EQ.0) then
        CALL GRSTOP1
        RETURN
      endif
4450  CONTINUE
C PLOT MASS BALANCE VS ELEVATION
      XMIN=1.d30
      YMIN=XMIN
      XMAX=-XMIN
      YMAX=-YMIN
      DO I=1,NUMNP
        if(htice(i).gt.SEALEV+1e-6) then
          XMAX=MAX(XMAX,ADOT(I))
          YMAX=MAX(YMAX,HTICE(I))
          XMIN=MIN(XMIN,ADOT(I))
          YMIN=MIN(YMIN,HTICE(I))
        endif
      ENDDO
c     YMIN=0.d0
c     YMAX=5000.d0
      IF(IPASS.EQ.0) CALL GRSTRT(600,600)
      DELX=(XMAX-XMIN)/5.d0
      DELY=(YMAX-YMIN)/5.d0
      CALL WINDOW(REAL(XMIN-DELX),REAL(XMAX+DELX),
     &            REAL(YMIN-DELY),REAL(YMAX+DELY))
      IGRAP=0
      CALL NEWPAG
      CALL LINCLR(1)
      CALL XAXIS('ADOT',5.,REAL(XMIN),REAL(XMAX),REAL(YMIN),REAL(YMAX))
      CALL YAXIS('ELEV',10.,REAL(XMIN),REAL(XMAX),REAL(YMIN),REAL(YMAX))
      CALL MOVE(0.,REAL(YMIN))
      CALL DRAW(0.,REAL(YMAX))
      DO I=1,NUMNP
        THIK=HTICE(I)-DEPB(I)
        IF(DEPB(I).LT.SEALEV) THEN
          FLOT=(1.d0-RATDEN)*(DEPB(I)-SEALEV)
          IF(HTICE(I).LT.FLOT) THEN
            THIK=0.0d0
          ENDIF
        ENDIF
        IF(THIK.GT.0.0) THEN
          CALL LINCLR(1)
        ELSE
          CALL LINCLR(2)
        ENDIF
        CALL MOVE(REAL(ADOT(I)),REAL(HTICE(I)))
        CALL DRAW(REAL(ADOT(I)),REAL(HTICE(I)))
        CALL marker(REAL(ADOT(I)),REAL(HTICE(I)),4)
        WRITE(19,202) ADOT(I),HTICE(I)
        !WRITE(29,202) TEMP(I),HTICE(I)
202     FORMAT(10X,G13.6,2X,G13.6,I13)
      ENDDO
      WRITE(19,202) -99999.,-1.
      WRITE(19,*) 'MASS BALANCE'
      WRITE(29,202) -99999.,-1.
      WRITE(29,*) 'TEMPERATURE'
C      CALL GRSTOP
C      IPASS=0
      GOTO 2100
4460  CONTINUE
      PRINT *,'WHOLE BOUNDARY:0, OR CURRENT SELECTED:1'
      if(iflush) call gflush
      READ(*,*) IWHOLE
      WRITE(99,*) IWHOLE
      IF(IWHOLE.EQ.1) THEN
        NUMGBC=ICOUNT-1
        DO I=1,NUMGBC
          J1=IFP(I)
          J2=IFP(I+1)
          IBFLUX(I,1)=J1
          IBFLUX(I,2)=J2
          BFLUX(I)=1.d0
        ENDDO
      ELSEIF(IWHOLE.EQ.0) THEN
        NUMGBC=0
        DO I=1,NUMCOL-1
          NUMGBC=NUMGBC+1
          J2=I
          J1=I+1
          IBFLUX(NUMGBC,1)=J1
          IBFLUX(NUMGBC,2)=J2
          BFLUX(NUMGBC)=1.d0
        ENDDO
        DO I=(NUMLEV-1)*NUMCOL+1,NUMCOL*NUMLEV-1
          NUMGBC=NUMGBC+1
          J1=I
          J2=I+1
          IBFLUX(NUMGBC,1)=J1
          IBFLUX(NUMGBC,2)=J2
          BFLUX(NUMGBC)=1.d0
        ENDDO
        DO I=1,NUMCOL*(NUMLEV-1),NUMCOL
          NUMGBC=NUMGBC+1
          J1=I
          J2=I+NUMCOL
          IBFLUX(NUMGBC,1)=J1
          IBFLUX(NUMGBC,2)=J2
          BFLUX(NUMGBC)=1.d0
        ENDDO
        DO I=NUMCOL,NUMCOL*(NUMLEV-1),NUMCOL
          NUMGBC=NUMGBC+1
          J2=I
          J1=I+NUMCOL
          IBFLUX(NUMGBC,1)=J1
          IBFLUX(NUMGBC,2)=J2
          BFLUX(NUMGBC)=1.d0
        ENDDO
      ENDIF
      GOTO 2100
9999  CALL GRSTOP1
      RETURN
      END
C==========================================================
      SUBROUTINE GNUMBR1(X,IDIGIT,IWIDE)
      CHARACTER*20 JUNK
      IF(ABS(X).GE..1 .AND. ABS(X).LE.999999.) THEN
        WRITE(JUNK,100) X
100     FORMAT(G13.6)
      ELSEIF(X.EQ.0.) THEN
        WRITE(JUNK,101) X
101     FORMAT(2X,F7.5)
      ELSE
        WRITE(JUNK,102) X
102     FORMAT(1PE9.2)
      ENDIF
C     WRITE(JUNK,103) X
103   FORMAT(F6.0)
      CALL TEXT(20,JUNK)
      END
C==========================================================
      SUBROUTINE XAXIS(LABEL,XLEN,XMIN,XMAX,YMIN,YMAX)
      CHARACTER*20 LABEL
      CALL TXFONT(0)
      XTIC=(XMAX-XMIN)/100.
      YTIC=(YMAX-YMIN)/100.
      CALL MOVE(XMIN,YMIN)
      DO X=0.,XLEN,1.
        XP=XMIN+X*(XMAX-XMIN)/XLEN
        CALL DRAW(XP,YMIN)
        CALL DRAW(XP,YMIN-3.*YTIC)
        CALL DRAW(XP,YMIN+3.*YTIC)
        CALL MOVE(XP-7.*XTIC,YMIN-7.*YTIC)
        CALL  GNUMBR(XP,2,8)
        CALL MOVE(XP,YMIN)
      ENDDO
      CALL MOVE(XMIN,YMIN)
      DO X=0.,XLEN,.1
        XP=XMIN+X*(XMAX-XMIN)/XLEN
        CALL MOVE(XP,YMIN)
        CALL DRAW(XP,YMIN-YTIC)
      ENDDO
      CALL MOVE(XMIN+30.*XTIC,YMIN-15.*YTIC)
      CALL TEXT(20,LABEL)
      END
C==========================================================
      SUBROUTINE YAXIS(LABEL,YLEN,XMIN,XMAX,YMIN,YMAX)
      CHARACTER*20 LABEL
      CALL TXFONT(0)
      CALL TXANGL(90.)
      XTIC=(XMAX-XMIN)/100.
      YTIC=(YMAX-YMIN)/100.
      CALL MOVE(XMIN,YMIN)
      DO Y=0.,YLEN,1.
        YP=YMIN+Y*(YMAX-YMIN)/YLEN
        CALL DRAW(XMIN,YP)
        CALL DRAW(XMIN-3.*XTIC,YP)
        CALL DRAW(XMAX,YP)
        CALL MOVE(XMIN,YP+1.*YTIC)
        CALL  GNUMBR(YP,2,8)
        CALL MOVE(XMIN,YP)
      ENDDO
      CALL MOVE(XMIN,YMIN)
      DO Y=0.,YLEN,.1
        YP=YMIN+Y*(YMAX-YMIN)/YLEN
        CALL MOVE(XMIN,YP)
        CALL DRAW(XMIN-XTIC,YP)
      ENDDO
      CALL MOVE(XMAX,YMIN)
      CALL DRAW(XMAX,YMAX)
      CALL MOVE(XMIN-15.*XTIC,YMIN+23.*YTIC)
C     CALL TXANGL(90.)
      CALL TEXT(20,LABEL)
      CALL TXANGL(0.)
      END
C==========================================================
      subroutine txangl(a)
      end
C==========================================================
      subroutine txfont(i)
      end
C==========================================================
      SUBROUTINE DTLIST(IPASS,TIME,TNSL)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
      COMMON /TTNSL/ TTIME(NSMAX),TLIST(NSMAX),NTNSL
      LOGICAL BATCH
      COMMON /BATCHER/ BATCH
      common /flush/ iflush
      logical iflush
      XMIN=.001d0*TTIME(1)
      XMAX=.001d0*TTIME(NTNSL)
      TMIN=1d30
      TMAX=-1d30
      DO I=1,NTNSL
        TMIN=MIN(TMIN,TLIST(I))
        TMAX=MAX(TMAX,TLIST(I))
      ENDDO
      IF(.NOT.BATCH) THEN
        IF(IPASS.EQ.0) THEN
          CALL GRSTRT(600,600)
          IPASS=1
        ENDIF
        CALL WINDOW(REAL(XMIN),REAL(XMAX),REAL(TMIN),REAL(TMAX))
        CALL MOVE(REAL(.001d0*TTIME(1)),REAL(TLIST(1)))
        DO I=2,NTNSL
          CALL DRAW(REAL(.001d0*TTIME(I)),REAL(TLIST(I)))
        ENDDO
        CALL LINCLR(2)
        CALL MOVE(REAL(0.001d0*TIME),REAL(TMIN))
        CALL DRAW(REAL(0.001d0*TIME),REAL(TMAX))
        CALL MOVE(REAL(XMIN),REAL(TNSL))
        CALL DRAW(REAL(XMAX),REAL(TNSL))
      ENDIF
      if(.false.) then
        WRITE(*,*) 'scale and offset',tmin,tmax
        if(iflush) call gflush
        READ(*,*) scale,offset
        WRITE(99,*) scale,offset
        TMIN=1d30
        TMAX=-1d30
        DO I=1,NTNSL
          tlist(i)=scale*tlist(i)+offset
          TMIN=MIN(TMIN,TLIST(I))
          TMAX=MAX(TMAX,TLIST(I))
        ENDDO
      else
        WRITE(*,*) 'new min and max',tmin,tmax
        if(iflush) call gflush
        READ(*,*) tnmin,tnmax
        WRITE(99,*) tnmin,tnmax
        scale=(tnmax-tnmin)/(tmax-tmin)
        offset=tnmin-scale*tmin
        print *,scale,offset
        TMIN=1d30
        TMAX=-1d30
        DO I=1,NTNSL
          tlist(i)=scale*tlist(i)+offset
          TMIN=MIN(TMIN,TLIST(I))
          TMAX=MAX(TMAX,TLIST(I))
        ENDDO
      endif
      WRITE(*,*) 'after scale and offset',tmin,tmax
      END
C===========================================
      FUNCTION IZSET(NCOLOR,ZZZ,ZF,ZD,IOFF)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      IMAX=NCOLOR+IOFF
      IMIN=IOFF+1
C       IC=IMAX-NINT((ZZZ-ZF)/ZD-.5d0)
        ZIC=(ZZZ-ZF)/ZD+.5d0
        IC=NINT(ZIC)
C       WRITE(*,*) ZZZ,ZIC,IC
      IF(IC.LT.IMIN) IC=IMIN
      IF(IC.GT.IMAX) IC=IMAX
      IZSET=IC
      RETURN
      END

      SUBROUTINE FLUXES(IPASS,NUMNP,NUMEL,XX,YY,HTICE,DEPB,KX,
     &                  NUMGBC,IBFLUX,BFLUX)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM)
      DIMENSION XX(NMAX),YY(NMAX),KX(NMAX,4)
      DIMENSION HTICE(NMAX),DEPB(NMAX)
      DIMENSION IBFLUX(NMAX,2),BFLUX(NMAX)
      common /flush/ iflush
      logical iflush
      SAVE FACTOR
      DATA FACTOR /1.D1/
      IF(IPASS.EQ.0) THEN
        CALL GRSTRT(600,600)
        IPASS=1
      ENDIF
      CALL LINCLR(1)
      if(ipass.eq.1) then
        PRINT *,'CURRENT FACTOR',FACTOR
        if(iflush) call gflush
        READ(*,*) FACTOR
        WRITE(99,*) FACTOR
      endif
      XMIN=1.d30
      YMIN=XMIN
      XMAX=-XMIN
      YMAX=-YMIN
      DO I=1,NUMNP
        XMAX=MAX(XMAX,XX(I))
        YMAX=MAX(YMAX,YY(I))
        XMIN=MIN(XMIN,XX(I))
        YMIN=MIN(YMIN,YY(I))
      ENDDO
      IF(XMAX-XMIN.GT.YMAX-YMIN) THEN
        YMAX=YMIN+XMAX-XMIN
      ELSE
        XMAX=XMIN+YMAX-YMIN
      ENDIF
      DELX=(XMAX-XMIN)/5.d0
      DELY=(YMAX-YMIN)/5.d0
      CALL WINDOW(REAL(XMIN-DELX),REAL(XMAX+DELX),
     &            REAL(YMIN-DELY),REAL(YMAX+DELY))
      DO I=1,NUMGBC
        XP=0.5d0*(XX(IBFLUX(I,1))+XX(IBFLUX(I,2)))
        YP=0.5d0*(YY(IBFLUX(I,1))+YY(IBFLUX(I,2)))
        XD=XX(IBFLUX(I,2))-XX(IBFLUX(I,1))
        YD=YY(IBFLUX(I,2))-YY(IBFLUX(I,1))
        DD=SQRT(XD**2+YD**2)
        IF(DD.GT.0.) THEN
          XD=XD/DD
          YD=YD/DD
        ELSE
          PRINT *,'ERROR, TWO BC NODES ARE THE SAME...'
          STOP
        ENDIF
c        PRINT *,I,REAL(XP),REAL(YP),REAL(BFLUX(I)),
c     &             REAL(BFLUX(I)/(HTICE(I)-DEPB(I)))
        CALL MOVE(REAL(XP),REAL(YP))
        XP1=XP+YD*FACTOR*BFLUX(I)
        YP1=YP-XD*FACTOR*BFLUX(I)
        CALL DRAW(REAL(XP1),REAL(YP1))
        XP=XX(IBFLUX(I,1))
        YP=YY(IBFLUX(I,1))
        CALL MOVE(REAL(XP),REAL(YP))
        XP=XX(IBFLUX(I,2))
        YP=YY(IBFLUX(I,2))
        CALL DRAW(REAL(XP),REAL(YP))
      ENDDO
      END
C===========================================
      SUBROUTINE WRITPROF(NUMNP,XX,YY,HTICE,BDROCK,ADOT)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM,NSMAX=MAXTIME)
      DIMENSION XX(NUMNP),YY(NUMNP),HTICE(NUMNP),BDROCK(NUMNP)
      DIMENSION ADOT(NUMNP)
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      COMMON /LINE/ NP,NLINE(1000)
      logical file20
1010  FORMAT(13I6)
      NP=0
      inquire(file='fort.20',exist=file20)
      if(file20) then
        REWIND 20
        READ(20,1010,END=102) NP
        READ(20,1010) (NLINE(IP),IP=1,NP)
      endif
102   CONTINUE
      RHOI=0.917
      RHOW=1.092
      RATIO=RHOI/(RHOW-RHOI)
      if(np.gt.0) then
        write(74,*) 'test'
c       WRITE(74,*) 0.0
c       write(74,*) '    LAT      LONG     SURF BED ADOT'
        j1=nline(1)
        DO I=1,NP
          J=NLINE(I)
          surf=htice(j)
          if(j.gt.0 .and. j.le.numnp) then
            IF(BDROCK(J).LT.SEALEV) THEN
              FLOT=(1.d0-RATDEN)*(BDROCK(J)-SEALEV)
              if(surf.le.flot) then
                if(surf.eq.0) surf=1
                BASE=-RATIO*surf
              ELSE
                BASE=BDROCK(J)
              endif
            else
              BASE=BDROCK(J)
            endif
            CALL RECPOL(0.001d0*XX(J),0.001d0*YY(J),RLAT,RLONG)
            place=sqrt((xx(j)-xx(j1))**2+(yy(j)-yy(j1))**2)
            WRITE(74,*) REAL(place),REAL(surf),
     &                  REAL(BASE),REAL(BDROCK(j)),REAL(ADOT(J))
            WRITE(*,*) REAL(place),REAL(surf),
     &                  REAL(BASE),REAL(BDROCK(j)),REAL(ADOT(J))
c           WRITE(74,*) REAL(RLAT),REAL(RLONG),REAL(HTICE(J)),
c    &                  REAL(BDROCK(J)),REAL(ADOT(J))
          endif
        ENDDO
      endif
      END

