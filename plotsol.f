      SUBROUTINE PLOTSOL(NUMNP,XX,YY,HTICE,DEPB,KODE,FRACT,
     &                   PSURF,WTHICK,TWATER,PWATER,AFUDGE,
     &                   NTSTEP,TTIME,VOL,AREA,TTBOT,TTAVG,TTNSL,
     &                   TTSEAL,
     &                   IPLOT,
     &                   NUMEL,KX,CONST,VEL,WWWMIN,DTMIN)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM, NSMAX=MAXTIME )
      COMMON /PICTURE/ XLL,XD,YLL,YD
      COMMON /LINE/ NP,NLINE(1000)
      COMMON /VELOS/ ASCAL,UMAX,VTHRESH,INORM
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      DIMENSION XX(NMAX),YY(NMAX),ZZ(NMAX),PSURF(NMAX)
      DIMENSION HTICE(NMAX),DEPB(NMAX),KODE(NMAX),FRACT(NMAX)
      DIMENSION WTHICK(NMAX),AFUDGE(NMAX)
      DIMENSION KX(NMAX,4),CONST(NMAX),VEL(NMAX,3)
      DIMENSION VOL(NSMAX),TTNSL(NSMAX),TWATER(NSMAX),TTSEAL(NSMAX),
     &          AREA(NSMAX),TTIME(NSMAX),TTBOT(NSMAX),
     &          PWATER(NSMAX),WWWMIN(NSMAX),TTAVG(NSMAX)
      DIMENSION ICMAP(16),DIST(1000)
      DATA ICMAP /0,5,11,4,12,6,13,2,8,7,9,3,10,14,15,1/
c
      IF(IPLOT.EQ.7) THEN
        RETURN
      ELSEIF(IPLOT.EQ.8) THEN
        CALL SCALEXY(xx,yy,NUMNP,xmin,xmax,ymin,ymax,delx,dely)
        CALL NEWPAG
        CALL WINDOW(REAL(XMIN-DELX),REAL(XMAX+DELX),
     &              REAL(YMIN-DELY),REAL(YMAX+DELY))
c        CALL CONTR(NMAX,NUMEL,XX,YY,KX,
c     &                 WTHICK,-0.999999999D0,10.000000001D0,1.D0,
c     &                 XMIN,XMAX,YMIN,YMAX)
        CALL CONTR(NMAX,NUMEL,XX,YY,KX,
     &                 WTHICK,-0.0999999999D0,1.000000001D0,0.1D0,
     &                 XMIN,XMAX,YMIN,YMAX)
C        CALL DVELO(NMAX,NUMEL,KX,XX,YY,HTICE,DEPB,CONST,
C     &             VEL,.TRUE.,DTMIN)
         call WVELO(XX, YY, KX, NUMNP, NUMEL,
     &              WTHICK, HTICE, DEPB, .false.)
        CALL DCOAST(15)
        CALL LINCLR(1)                                                    
        CALL MOVE(REAL(XMIN),REAL(YMIN))                                  
        CALL DRAW(REAL(XMIN),REAL(YMAX))                                  
        CALL DRAW(REAL(XMAX),REAL(YMAX))                                  
        CALL DRAW(REAL(XMAX),REAL(YMIN))                                  
        CALL DRAW(REAL(XMIN),REAL(YMIN))                                  
        RETURN
      ELSEIF(IPLOT.EQ.6) THEN
        CALL SCALEXY(xx,yy,NUMNP,xmin,xmax,ymin,ymax,delx,dely)
        CALL NEWPAG
        CALL WINDOW(REAL(XMIN-DELX),REAL(XMAX+DELX),
     &              REAL(YMIN-DELY),REAL(YMAX+DELY))
        CALL CONTR(NMAX,NUMEL,XX,YY,KX,
     &                 HTICE,-999.999D0,19000.D0,1000.D0,
c    &                 HTICE,-5500.D0,0.D0,250.D0,
     &                 XMIN,XMAX,YMIN,YMAX)
        CALL DVELO(NMAX,NUMEL,KX,XX,YY,HTICE,DEPB,CONST,
     &             VEL,.TRUE.,DTMIN)
        CALL DCOAST(15)
        CALL LINCLR(1)                                                    
        CALL MOVE(REAL(XMIN),REAL(YMIN))                                  
        CALL DRAW(REAL(XMIN),REAL(YMAX))                                  
        CALL DRAW(REAL(XMAX),REAL(YMAX))                                  
        CALL DRAW(REAL(XMAX),REAL(YMIN))                                  
        CALL DRAW(REAL(XMIN),REAL(YMIN))                                  
        RETURN
      ELSEIF(IPLOT.EQ.5) THEN
        ZMIN=1d30
        ZMAX=-ZMIN
        DMIN=1d30
        DMAX=-ZMIN
        ZMIN=-1000.d0
        XBASE=XX(NLINE(1))
        YBASE=YY(NLINE(1))
c        DO I=1,NP
c          DIST(I)=SQRT((XX(NLINE(I))-XBASE)**2+
c     &                 (YY(NLINE(I))-YBASE)**2)
        DIST(1)=0.0d0
        DMIN=MIN(DMIN,DIST(1))
        DMAX=MAX(DMAX,DIST(1))
        ZMIN=MIN(ZMIN,HTICE(NLINE(1)))
        ZMAX=MAX(ZMAX,HTICE(NLINE(1)))
        DO I=2,NP
          DIST(I)=DIST(I-1)+
     &            SQRT((XX(NLINE(I))-XX(NLINE(I-1)))**2+
     &                 (YY(NLINE(I))-YY(NLINE(I-1)))**2)
          DMIN=MIN(DMIN,DIST(I))
          DMAX=MAX(DMAX,DIST(I))
          ZMIN=MIN(ZMIN,HTICE(NLINE(I)))
          ZMAX=MAX(ZMAX,HTICE(NLINE(I)))
        ENDDO
C        DO I=1,NP
C          ZMIN=MIN(ZMIN,DEPB(NLINE(I)))
C          ZMAX=MAX(ZMAX,DEPB(NLINE(I)))
C        ENDDO
        DO I=1,NP
          ZMIN=MIN(ZMIN,PSURF(NLINE(I)))
          ZMAX=MAX(ZMAX,PSURF(NLINE(I)))
        ENDDO
c set specific vertical range 
        ZMAX=19000.d0
        ZMIN=-1000.d0
c set specific vertical range 
        ZMAX=0.d0
        ZMIN=-6000.d0
c set specific vertical range 
        ZMAX=10000.d0
        ZMIN=-1000.d0
        CALL NEWPAG
        CALL WINDOW(REAL(DMIN),REAL(DMAX),REAL(ZMIN),REAL(ZMAX))
        CALL LINCLR(1)
        CALL MOVE(0.,0.)
        CALL DRAW(REAL(NP),0.)
        CALL MOVE(REAL(DIST(1)),REAL(HTICE(NLINE(1))))
        DO I=2,NP
          CALL DRAW(REAL(DIST(I)),REAL(HTICE(NLINE(I))))
        ENDDO
        CALL LINCLR(3)
        CALL MOVE(REAL(DIST(1)),REAL(DEPB(NLINE(1))))
        DO I=2,NP
          CALL DRAW(REAL(DIST(I)),REAL(DEPB(NLINE(I))))
        ENDDO
        CALL LINCLR(2)
        CALL MOVE(REAL(DIST(1)),REAL(PSURF(NLINE(1))))
        DO I=2,NP
          CALL DRAW(REAL(DIST(I)),REAL(PSURF(NLINE(I))))
        ENDDO
        CALL LINCLR(4)
        CALL MOVE(REAL(DIST(1)),
     &            REAL(-250+500*AFUDGE(NLINE(1))))
        DO I=2,NP
          CALL DRAW(REAL(DIST(I)),
     &              REAL(-250+500*AFUDGE(NLINE(I))))
        ENDDO
        call gflush()
        RETURN
      ELSEIF(IPLOT.EQ.4) THEN
        DO I=1,NUMNP
          ZZ(I)=FRACT(I)
        ENDDO
      ELSEIF(IPLOT.EQ.3) THEN
        RHOI=.917d0
        RHOW=1.092d0
        DO I=1,NUMNP
          ZZ(I)=HTICE(I)-DEPB(I)
          IF(DEPB(I).LT.SEALEV) THEN
            FLOT=(1.d0-RATDEN)*(DEPB(I)-SEALEV)
            IF(HTICE(I).LT.FLOT) THEN
              ZZ(I)=0.0
            ENDIF
          ENDIF
        ENDDO
      ELSEIF(IPLOT.EQ.2) THEN
        DO I=1,NUMNP
          ZZ(I)=HTICE(I)
        ENDDO
       ELSEIF(IPLOT.EQ.1) THEN
         continue
      ELSEIF(IPLOT.EQ.0) THEN
        RETURN
      ELSE
        RETURN
      ENDIF
      XLEN=100.d0
      YLEN=100.d0
      ZLEN=14.d0
      IF(IPLOT.EQ.4 .OR. IPLOT.EQ.3 .OR. IPLOT.EQ.2) THEN
        CALL SCALE3(XX,XLEN,NUMNP,1)
        CALL SCALE3(YY,YLEN,NUMNP,1)
c       CALL SCALE2(ZZ,ZLEN,NUMNP,1)
C ... FOLLOWING TO MAKE X AND Y AXIS SAME SCALE
        IF(XX(NUMNP+2).GT.YY(NUMNP+2)) THEN
          YY(NUMNP+2)=XX(NUMNP+2)
        ELSE
          XX(NUMNP+2)=YY(NUMNP+2)
        ENDIF
C ... END OF SAME SCALE
        XLL=XX(NUMNP+1)
        YLL=YY(NUMNP+1)
        XD=XX(NUMNP+2)
        YD=YY(NUMNP+2)
        XUR=XLL+XLEN*XD
        YUR=YLL+YLEN*YD
        ZF=ZZ(NUMNP+1)
        ZD=ZZ(NUMNP+2)
      ELSEIF(IPLOT.EQ.1) THEN
         continue
c        CALL NEWPAG
c        XLEN=100.
c        YLEN=100.
c        CALL SCALE2(TTIME,XLEN,NTSTEP,1)
c        XLL=TTIME(NTSTEP+1)
c        XD=TTIME(NTSTEP+2)
c        CALL SCALE2(VOL,YLEN,NTSTEP,1)
c        CALL SCALE2(AREA,YLEN,NTSTEP,1)
c        CALL SCALEB(NTSTEP,VOL,AREA,YLEN,YLL,YD)
      ELSE
        RETURN
      ENDIF
      IF(IPLOT.EQ.4) THEN
        ZF=0.d0
        ZD=.1d0
        DO I=1,NUMNP
          RZ=(ZZ(I)-ZF)/ZD
          IZ=INT(RZ+2)
          IF(IZ.LT.2) IZ=2
          IF(IZ.GT.16) IZ=16
          XXX=(XX(I)-XLL)/XD
          YYY=(YY(I)-YLL)/YD
          CALL MRKCLR(ICMAP(IZ))
c          CALL MARKER(REAL(XXX),REAL(YYY),1)
          CALL POINT(REAL(XXX),REAL(YYY))
        ENDDO
      ELSEIF(IPLOT.EQ.3) THEN
        ZF=-100.d0
        ZD=200.d0
        DO I=1,NUMNP
          RZ=(ZZ(I)-ZF)/ZD
          IZ=INT(RZ+2)
          IF(IZ.LT.2) IZ=2
          IF(IZ.GT.16) IZ=16
          XXX=(XX(I)-XLL)/XD
          YYY=(YY(I)-YLL)/YD
          IF(ZZ(I).GT.1.) THEN
            CALL MRKCLR(ICMAP(IZ))
c            CALL MARKER(REAL(XXX),REAL(YYY),1)
            CALL POINT(REAL(XXX),REAL(YYY))
          ELSE
            CALL MRKCLR(0)
C            CALL MARKER(REAL(XXX),REAL(YYY),1)
            CALL POINT(REAL(XXX),REAL(YYY))
          ENDIF
        ENDDO
      ELSEIF(IPLOT.EQ.2) THEN
        ZF=-1000.d0
        ZD=500.d0
        DO I=1,NUMNP
          RZ=(ZZ(I)-ZF)/ZD
          IZ=INT(RZ+2)
          IF(IZ.LT.2) IZ=2
          IF(IZ.GT.16) IZ=16
          XXX=(XX(I)-XLL)/XD
          YYY=(YY(I)-YLL)/YD
          CALL MRKCLR(ICMAP(IZ))
C          CALL MARKER(REAL(XXX),REAL(YYY),1)
          CALL POINT(REAL(XXX),REAL(YYY))
        ENDDO
      ELSEIF(IPLOT.EQ.1) THEN
c          XXX=(TTIME(1)-XLL)/XD
c          YYY=(VOL(1)-YLL)/YD
c          CALL LINCLR(1)
c          CALL MOVE(REAL(XXX),REAL(YYY))
c          DO I=1,NTSTEP
c            XXX=(TTIME(I)-XLL)/XD
c            YYY=(VOL(I)-YLL)/YD
c            CALL DRAW(REAL(XXX),REAL(YYY))
c          ENDDO
c          XXX=(TTIME(1)-XLL)/XD
c          YYY=(AREA(1)-YLL)/YD
c          CALL LINCLR(2)
c          CALL MOVE(REAL(XXX),REAL(YYY))
c          DO I=1,NTSTEP
c            XXX=(TTIME(I)-XLL)/XD
c            YYY=(AREA(I)-YLL)/YD
c            CALL DRAW(REAL(XXX),REAL(YYY))
c          ENDDO
c-----------------------------------------
C ........ PLOT VOLUME VS TIME
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
             VMAX=MAX(VMAX,VOL(I))
             AMAX=MAX(AMAX,AREA(I))
             TMAX=MAX(TMAX,TTBOT(I))
             TAMAX=MAX(TAMAX,TTAVG(I))
             TTMAX=MAX(TTMAX,TTNSL(I))
             SLMAX=MAX(SLMAX,TTSEAL(I))
             WTMAX=MAX(WTMAX,TWATER(I))
             PWMAX=MAX(PWMAX,PWATER(I))
             WWMAX=MAX(WWMAX,WWWMIN(I))
             XMIN=MIN(XMIN,TTIME(I))
             VMIN=MIN(VMIN,VOL(I))
             AMIN=MIN(AMIN,AREA(I))
             TMIN=MIN(TMIN,TTBOT(I))
             TAMIN=MIN(TAMIN,TTAVG(I))
             TTMIN=MIN(TTMIN,TTNSL(I))
             SLMIN=MIN(SLMIN,TTSEAL(I))
             WTMIN=MIN(WTMIN,TWATER(I))
             PWMIN=MIN(PWMIN,PWATER(I))
             WWMIN=MIN(WWMIN,WWWMIN(I))
           ENDDO
           DX=(XMAX-XMIN)/20.d0
           DA=(AMAX-AMIN)/20.d0
           DTT=(TMAX-TMIN)/20.d0
           DTTA=(TAMAX-TAMIN)/20.d0
           DTTT=(TTMAX-TTMIN)/20.d0
           DTSL=(SLMAX-SLMIN)/20.d0
           DWT=(WTMAX-WTMIN)/20.d0
           DPW=(PWMAX-PWMIN)/20.d0
           DWW=(WWMAX-WWMIN)/20.d0
           DV=(VMAX-VMIN)/20.d0
           IF(DX.EQ.0.) DX=1.d0
           IF(DA.EQ.0.) DA=1.d0
           IF(DTT.LT.1E-6) DTT=1.d0
           IF(DTTA.LT.1E-6) DTTA=1.d0
           IF(DTTT.EQ.0.) DTTT=1.d0
           IF(DTSL.EQ.0.) DTSL=1.d0
           IF(DWT.EQ.0.) DWT=1.d0
           IF(DPW.EQ.0.) DPW=1.d0
           IF(DWW.EQ.0.) DWW=1.d0
           IF(DV.EQ.0.) DV=1.d0
c           PRINT *,XMIN,XMAX
c           PRINT *,VMIN,VMAX
c           PRINT *,AMIN,AMAX
c           PRINT *,TMIN,TMAX
c           PRINT *,TAMIN,TAMAX
c           PRINT *,TTMIN,TTMAX
c           PRINT *,WTMIN,WTMAX
c           PRINT *,PWMIN,PWMAX
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
           CALL NEWPAG
           CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &                 REAL(VMIN),REAL(VMAX+5*DV))
           CALL LINCLR(1)
           CALL MOVE(REAL(TTIME(NTSTEP)),REAL(VOL(2)))
           CALL DRAW(REAL(TTIME(2)),REAL(VOL(2)))
           DO I=2,NTSTEP
             CALL DRAW(REAL(TTIME(I)),REAL(VOL(I)))
           ENDDO
           CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &                 REAL(AMIN-DA),REAL(AMAX+4*DA))
           CALL LINCLR(2)
           CALL MOVE(REAL(TTIME(NTSTEP)),REAL(AREA(2)))
           CALL DRAW(REAL(TTIME(2)),REAL(AREA(2)))
           DO I=2,NTSTEP
             CALL DRAW(REAL(TTIME(I)),REAL(AREA(I)))
           ENDDO

           CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &                 REAL(TAMIN-2*DTTA),REAL(TAMAX+3*DTTA))
           CALL LINCLR(8)
           CALL MOVE(REAL(TTIME(NTSTEP)),REAL(TTAVG(2)))
           CALL DRAW(REAL(TTIME(2)),REAL(TTAVG(2)))
           DO I=2,NTSTEP
             CALL DRAW(REAL(TTIME(I)),REAL(TTAVG(I)))
           ENDDO

           CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &                 REAL(TMIN-2*DTT),REAL(TMAX+3*DTT))
           CALL LINCLR(3)
           CALL MOVE(REAL(TTIME(NTSTEP)),REAL(TTBOT(2)))
           CALL DRAW(REAL(TTIME(2)),REAL(TTBOT(2)))
           DO I=2,NTSTEP
             CALL DRAW(REAL(TTIME(I)),REAL(TTBOT(I)))
           ENDDO

           CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &                 REAL(TTMIN-3*DTTT),REAL(TTMAX+2*DTTT))
           CALL LINCLR(4)
           CALL MOVE(REAL(TTIME(NTSTEP)),REAL(TTNSL(2)))
           CALL DRAW(REAL(TTIME(2)),REAL(TTNSL(2)))
           DO I=2,NTSTEP
             CALL DRAW(REAL(TTIME(I)),REAL(TTNSL(I)))
           ENDDO
           CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &                 REAL(SLMIN-3*DTSL),REAL(SLMAX+2*DTSL))
           CALL LINCLR(11)
           CALL MOVE(REAL(TTIME(NTSTEP)),REAL(TTSEAL(2)))
           CALL DRAW(REAL(TTIME(2)),REAL(TTSEAL(2)))
           DO I=2,NTSTEP
             CALL DRAW(REAL(TTIME(I)),REAL(TTSEAL(I)))
           ENDDO

           CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &                 REAL(WTMIN-4*DWT),REAL(WTMAX+1*DWT))
           CALL LINCLR(5)
           CALL MOVE(REAL(TTIME(NTSTEP)),REAL(TWATER(2)))
           CALL DRAW(REAL(TTIME(2)),REAL(TWATER(2)))
           DO I=2,NTSTEP
             CALL DRAW(REAL(TTIME(I)),REAL(TWATER(I)))
           ENDDO
           CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &                 REAL(PWMIN-4*DPW),REAL(PWMAX+1*DPW))
           CALL LINCLR(6)
           CALL MOVE(REAL(TTIME(NTSTEP)),REAL(PWATER(2)))
           CALL DRAW(REAL(TTIME(2)),REAL(PWATER(2)))
           DO I=2,NTSTEP
             CALL DRAW(REAL(TTIME(I)),REAL(PWATER(I)))
           ENDDO

           CALL WINDOW(REAL(XMIN),REAL(XMAX),
     &                 REAL(WWMIN-5*DWW),REAL(WWMAX+0*DWW))
           CALL LINCLR(7)
           CALL MOVE(REAL(TTIME(NTSTEP)),REAL(WWWMIN(2)))
           CALL DRAW(REAL(TTIME(2)),REAL(WWWMIN(2)))
           DO I=2,NTSTEP
             CALL DRAW(REAL(TTIME(I)),REAL(WWWMIN(I)))
           ENDDO
c-------------------------------
        RETURN
      ELSEIF(IPLOT.EQ.0) THEN
        RETURN
      ELSE
        RETURN
      ENDIF  
      END
      SUBROUTINE DRAWOUT
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(MXX=MAXNUM)
      COMMON /PICTURE/ XLL,XD,YLL,YD
      DIMENSION XP(MXX),YP(MXX),IR(MXX)
      DIMENSION XPOINT(MXX),YPOINT(MXX)
      SAVE II,XP,YP,IR,XPOINT,YPOINT,NPOINT,IDONE
      DATA IDONE /0/
C
      CALL LINCLR(15)

      IF(IDONE.EQ.0) THEN
        IDONE=1
        II=0
C ..... READ AND PLOT OUTLINE IN FILE OUTLINE DATA B
        REWIND 10
        IREAD=0
123     CONTINUE
        READ(10,*,END=124) XOUT,YOUT,IREAD
        II=II+1
        IF(XOUT.EQ.-99999.) THEN
          II=II-1
          READ(10,1000) JUNK
1000  FORMAT(A80)
          GOTO 123
        ENDIF
        IF(XOUT.EQ.-99999.) GOTO 124
        XOUT=(XOUT*1000.d0-XLL)/XD
        XP(II)=XOUT
        YOUT=(YOUT*1000.d0-YLL)/YD
        YP(II)=YOUT
        IR(II)=IREAD
        IF(IREAD.EQ.1) THEN
          CALL MOVE(REAL(XOUT),REAL(YOUT))
        ELSE
          CALL DRAW(REAL(XOUT),REAL(YOUT))
        ENDIF
        IREAD=1
        GOTO 123
124     REWIND 10
        REWIND 46
        NPOINT=0
        DO I=1,MXX
          READ(46,*,END=125) XPOINT(I),YPOINT(I)
          XPOINT(I)=(XPOINT(I)*1000.d0-XLL)/XD
          YPOINT(I)=(YPOINT(I)*1000.d0-YLL)/YD
          NPOINT=I
        ENDDO
        PRINT *,'TOO MANY IN DRAWOUT '
        STOP
125     REWIND 46
      ELSE
        DO I=1,II
        IF(IR(I).EQ.1) THEN
          CALL MOVE(REAL(XP(I)),REAL(YP(I)))
        ELSE
          CALL DRAW(REAL(XP(I)),REAL(YP(I)))
        ENDIF
        ENDDO 
        DO I=1,NPOINT
          CALL POINT(REAL(XPOINT(I)),REAL(YPOINT(I)))
        ENDDO
      ENDIF
c
      END
      SUBROUTINE SCALEB(NTSTEP,VOL,AREA,YLEN,YLL,YD)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VOL(NTSTEP),AREA(NTSTEP)
      RMIN=1d30
      RMAX=-RMIN
      DO I=1,NTSTEP
        RMIN=MIN(RMIN,VOL(I))
        RMIN=MIN(RMIN,AREA(I))
        RMAX=MAX(RMAX,VOL(I))
        RMAX=MAX(RMAX,AREA(I))
      ENDDO
      IF(RMIN.EQ.RMAX) RMAX=RMIN+1.d0
      YLL=RMIN
      YD=(RMAX-RMIN)/YLEN
      END
C====================================
      SUBROUTINE CONTR(NMAX,NUMEL,XX,YY,KX,
     &                 ZZ,RMIN,RMAX,RLEVSP,
     &                 XMIN,XMAX,YMIN,YMAX)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION ZZ(NMAX),XX(NMAX),YY(NMAX),KX(NMAX,4),IK(5)
      LOGICAL FOUND
      DIMENSION ICMAP(16)
c      DATA ICMAP /15,4,12,6,13,2,8,7,9,3,10,5,11,14,1,15/                
      DATA ICMAP /5,11,4,12,6,13,2,8,7,9,3,10,14,15,1,0/
c      DATA ICMAP /1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,1/                
c draw contours
      LEVEL=INT((RMAX-RMIN)/RLEVSP+1)
      DLEV=RLEVSP
      IF(LEVEL.LE.0) THEN
        LEVEL=2-LEVEL
        DLEV=-DLEV
      ENDIF
      LCONT=1
      XPOS=XMAX+(XMAX-XMIN)*0.025d0
      YPOS=YMIN+(YMAX-YMIN)*.25d0
      ICO=1
c      CALL LINCLR(ICO)
      CALL LINCLR(ICMAP(ICO))
      DO LEV=1,LEVEL
        FOUND=.FALSE.
        ICOUNT=0
        VAL=RMIN+(LEV-1)*DLEV
        DO N=1,NUMEL
          LCONT=1
          NNODE=4                                                         
          DO J=1,NNODE
            IK(J)=KX(N,J)
          ENDDO
          IK(NNODE+1)=KX(N,1)
          ICOUNT=0
          DO NN=1,NNODE
            I=IK(NN)
            J=IK(NN+1)
            IF(ZZ(I).LT.VAL .AND. VAL.LE.ZZ(J) .OR. 
     &         ZZ(J).LT.VAL .AND. VAL.LE.ZZ(I)) THEN
              FOUND=.TRUE.
              ICOUNT=ICOUNT+1
              XINT=XX(J)+(VAL-ZZ(J))*(XX(I)-XX(J))/(ZZ(I)-ZZ(J))
              YINT=YY(J)+(VAL-ZZ(J))*(YY(I)-YY(J))/(ZZ(I)-ZZ(J))
              XINT=XINT*.001d0
              YINT=YINT*.001d0
              IF(LCONT.EQ.1) THEN
                CALL MOVE(REAL(XINT),REAL(YINT))
                LCONT=2
              ELSE
                CALL DRAW(REAL(XINT),REAL(YINT))
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        IF(FOUND) THEN
          CALL MOVE(REAL(XPOS),REAL(YPOS))
c          CALL TXTCLR(ICO)
          CALL TXTCLR(ICMAP(ICO))
          CALL RNUMBR(REAL(VAL),3,9)
          YPOS=YPOS+(YMAX-YMIN)*.04d0
        ENDIF
        ICO=ICO+1
        IF(ICO.GT.15) ICO=1
c        CALL LINCLR(ICO)
        CALL LINCLR(ICMAP(ICO))
      ENDDO
      END
C==========================================================
      SUBROUTINE DCOAST(ICOLOR)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(MXX=MAXNUM)
      DIMENSION XP(MXX),YP(MXX),IR(MXX)
      DIMENSION XPOINT(MXX),YPOINT(MXX)
      SAVE II,XP,YP,IR,XPOINT,YPOINT,NPOINT,IDONE
      DATA IDONE /0/
C
      CALL LINCLR(ICOLOR)

      IF(IDONE.EQ.0) THEN
        IDONE=1
        II=0
C ..... READ AND PLOT OUTLINE IN FILE OUTLINE DATA B
        REWIND 10
        IREAD=0
123     CONTINUE
        READ(10,*,END=124) XOUT,YOUT,IREAD
        II=II+1
        IF(XOUT.EQ.-99999.) THEN
          II=II-1
          READ(10,1000) JUNK
1000  FORMAT(A80)
          GOTO 123
        ENDIF
        IF(XOUT.EQ.-99999.) GOTO 124
        XP(II)=XOUT
        YP(II)=YOUT
        IR(II)=IREAD
        IF(IREAD.EQ.1) THEN
          CALL MOVE(REAL(XOUT),REAL(YOUT))
        ELSE
          CALL DRAW(REAL(XOUT),REAL(YOUT))
        ENDIF
        IREAD=1
        GOTO 123
124     REWIND 10
        REWIND 46
        NPOINT=0
        DO I=1,MXX
          READ(46,*,END=125) XPOINT(I),YPOINT(I)
          NPOINT=I
        ENDDO
        PRINT *,'TOO MANY IN DRAWOUT '
        STOP
125     REWIND 46
      ELSE
        DO I=1,II
        IF(IR(I).EQ.1) THEN
          CALL MOVE(REAL(XP(I)),REAL(YP(I)))
        ELSE
          CALL DRAW(REAL(XP(I)),REAL(YP(I)))
        ENDIF
        ENDDO 
        DO I=1,NPOINT
          CALL POINT(REAL(XPOINT(I)),REAL(YPOINT(I)))
        ENDDO
      ENDIF
c
      END
C==========================================================
      SUBROUTINE DCOAST1(ICOLOR)
      CHARACTER*80 JUNK
C READ AND PLOT OUTLINE IN FILE OUTLINE DATA B
      CALL LINCLR(ICOLOR)
      IREAD=0
123   READ(10,*,END=124) XOUT,YOUT,IREAD
        IF(XOUT.EQ.-99999.) THEN
          READ(10,1000) JUNK
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
125   READ(46,*,END=126) XOUT,YOUT
        CALL POINT(REAL(XOUT),REAL(YOUT))
        goto 125
126   REWIND 46
1000  FORMAT(A80)
      END
C==========================================================
      SUBROUTINE DVELO(NMAX,NUMEL,KX,XX,YY,HTICE,DEPB,CONST,VEL,IGO,
     &                 DTMIN)
      IMPLICIT REAL*8(A-H,O-Z) 
      LOGICAL IGO                                         
      DIMENSION HTICE(NMAX),DEPB(NMAX),XX(NMAX),YY(NMAX)
      DIMENSION KX(NMAX,4),CONST(NMAX),VEL(NMAX,3)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2),LM(4)             
      DIMENSION PSI(4),DPSI(4,2)                                        
      DIMENSION XY(2,4),XARO(2),YARO(2)
      DIMENSION ICMAP(16)
      COMMON /VELOS/ ASCAL,UMAX,VTHRESH,INORM
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      DATA ICMAP /0,4,12,6,13,2,8,7,9,3,10,5,11,1,14,15/   
      DTMIN=1.d30
      VMAX=-1.d30
      HMIN=10.d0
      UDELTA=UMAX/13.d0
      DO 600 J=1,NUMEL  
        VEL(J,1)=0D0
        VEL(J,2)=0D0
        VEL(J,3)=0D0
        NNODE=4                                                         
        CENTX=0.0D00                                                    
        CENTY=0.0D00                                                    
        HH=0.0D00                                                         
        SUMX=0.0D00                                                       
        SUMY=0.0D00                                                       
        DO I=1,NNODE                                                  
          LM(I)=KX(J,I)                                                     
        ENDDO
        I=LM(1)                                                           
        JJ=LM(2)                                                          
        K=LM(3)                                                           
        L=LM(4)                                                           
        XY(1,1)=XX(I)                                                     
        XY(1,2)=XX(JJ)                                                    
        XY(1,3)=XX(K)                                                     
        XY(1,4)=XX(L)                                   
        XY(2,1)=YY(I)                                                     
        XY(2,2)=YY(JJ)                                                    
        XY(2,3)=YY(K)                                                     
        XY(2,4)=YY(L)                                   
        XCENT=(XX(I)+XX(JJ)+XX(K)+XX(L))/4000.d0
        YCENT=(YY(I)+YY(JJ)+YY(K)+YY(L))/4000.d0
        DISTMAX=0.D0
        XC=XCENT*1000.D0
        YC=YCENT*1000.D0
        DO I=1,4
          DIST=SQRT( (XC-XX(LM(I)))**2+(YC-YY(LM(I)))**2)
          DISTMAX=MAX(DISTMAX,DIST)
        ENDDO
        DISTMAX=DISTMAX*2D0
C                                                                       
        CALL FESHAPE(1,CENTX,CENTY,PSI,DPSI)                         
      
C CALCULATE DXDS...EQUATION (5.3.6)                                     
        DO I=1,2                                                      
          DO L=1,2                                                      
            DXDS(I,L)=0.0d0
            DO K=1,NNODE                                                  
              DXDS(I,L)=DXDS(I,L)+DPSI(K,L)*XY(I,K)
            ENDDO
          ENDDO
        ENDDO
   
C CALCULATE DSDX...EQUATION (5.2.7)                                     
        DETJ=(DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1))                    
        IF (DETJ.LE.0.0) THEN
          WRITE(*,*) ' BAD JACOBIAN',DETJ,X                                             
          STOP                                                              
        ENDIF
        DSDX(1,1)=DXDS(2,2)/DETJ                                          
        DSDX(2,2)=DXDS(1,1)/DETJ                                          
        DSDX(1,2)=-DXDS(1,2)/DETJ                                         
        DSDX(2,1)=-DXDS(2,1)/DETJ                                         
C CALCULATE D(PSI)/DX...EQUATION (5.3.5)                                
        DO I=1,NNODE                                                  
          DPSIX(I)=DPSI(I,1)*DSDX(1,1)+DPSI(I,2)*DSDX(2,1)                  
          DPSIY(I)=DPSI(I,1)*DSDX(1,2)+DPSI(I,2)*DSDX(2,2)                  
        ENDDO                                                          
        DO I=1,NNODE                                                  
          SUMX=SUMX + HTICE(LM(I))*DPSIX(I)                                 
          SUMY=SUMY + HTICE(LM(I))*DPSIY(I)                                 
          HH=HH + (HTICE(LM(I))-DEPB(LM(I)))*PSI(I)                                       
        ENDDO                                                          
C                                                                       
        DELH=SUMX**2 + SUMY**2                                            
        DELH=SQRT(DELH)                                                   
        IF(HH.GT.HMIN) THEN                                               
          UX=-CONST(J)*SUMX/HH                                            
          UY=-CONST(J)*SUMY/HH                                            
        ELSE                                                              
          UX=0.d0
          UY=0.d0
        ENDIF                                                             
        UMAG=SQRT(UX**2+UY**2)                                        
        IF(UMAG.GT.VTHRESH) THEN
          VEL(J,1)=0.0d0
          VEL(J,2)=0.0d0
          VEL(J,3)=0.0d0
        ELSE                                       
          VEL(J,1)=UMAG
          VEL(J,2)=UX
          VEL(J,3)=UY
        ENDIF
        IF(UMAG.GT.0.) THEN
          DTELEM=DISTMAX/UMAG
c         DTELEM=DISTMAX/1d4
        ELSE
          DTELEM=1d30
        ENDIF
        DTMIN=MIN(DTMIN,DTELEM)                                               
        VMAX=MAX(VMAX,UMAG)                                               
        IF(UMAG.GT.VTHRESH .OR. .NOT.IGO) GOTO 600                                       
        IF(UMAG.EQ.0.) THEN                                               
          ICOLOR=1                                                        
        ELSE                                                              
          ICOLOR=1+NINT(UMAG/UDELTA)                                    
        ENDIF                                                             
        IF(ICOLOR.LT.2) ICOLOR=2                                          
        IF(ICOLOR.GT.14) ICOLOR=14                                        
        XARO(1)=XCENT                                                     
        YARO(1)=YCENT                                                     
        IF(UMAG.NE.0) THEN
          IF(INORM.EQ.1) THEN
            XARO(2)=XCENT+ASCAL*UX/UMAG
            YARO(2)=YCENT+ASCAL*UY/UMAG
          ELSE
            XARO(2)=XCENT+ASCAL*UX
            YARO(2)=YCENT+ASCAL*UY
          ENDIF
        ELSE
          XARO(2)=XCENT                                                     
          YARO(2)=YCENT                                                     
        ENDIF
        CALL LINCLR(ICMAP(ICOLOR))                                      
        ICSAV=ICOLOR                                                    
        CALL MOVE(REAL(XARO(1)),REAL(YARO(1)))
        CALL DRAW(REAL(XARO(2)),REAL(YARO(2)))
600   CONTINUE
      IF(IOTOGG) then
        WRITE(list(ipage+1),*) 'UMAX=',VMAX,' DTMIN=',DTMIN
        ipage=ipage+1
      endif
      END
C==========================================================
      SUBROUTINE DGRAD(NMAX,NUMEL,KX,XX,YY,HTICE,DEPB,IGO)
      IMPLICIT REAL*8(A-H,O-Z) 
      LOGICAL IGO                                         
      DIMENSION HTICE(NMAX),DEPB(NMAX),XX(NMAX),YY(NMAX)
      DIMENSION KX(NMAX,4)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2),LM(4)             
      DIMENSION PSI(4),DPSI(4,2)                                        
      DIMENSION XY(2,4),XARO(2),YARO(2)
      DIMENSION ICMAP(16)
      SAVE ASCAL,UMAX,VTHRESH,INORM
      DATA ICMAP /0,4,12,6,13,2,8,7,9,3,10,5,11,1,14,15/   
C      DATA ASCAL /50.D0/, UMAX /0.01D0/, VTHRESH /1.D0/, INORM /1/
      DATA ASCAL /10000.D0/, UMAX /0.01D0/, VTHRESH /1.D0/, INORM /0/
      WRITE(*,*) 'INPUT SCALE FACTOR FOR ARROWS, UMAX, AND THRESHOLD'
      WRITE(*,*) 'INPUT 1 FOR NORMALIZED VECTORS, 0 FOR ABSOLUTE'  
      WRITE(*,*) ASCAL,UMAX,VTHRESH,INORM                                       
      READ(*,*) ASCAL,UMAX,VTHRESH,INORM                                       
      WRITE(99,*) ASCAL,UMAX,VTHRESH,INORM                                       
      VMAX=-1.d30
      HMIN=10.d0
      UDELTA=UMAX/13.d0
      DO 600 J=1,NUMEL  
        NNODE=4                                                         
        CENTX=0.0D00                                                    
        CENTY=0.0D00                                                    
        SUMHX=0.0D00                                                       
        SUMHY=0.0D00                                                       
        SUMBX=0.0D00                                                       
        SUMBY=0.0D00                                                       
        HH=0.0D00
        DO I=1,NNODE                                                  
          LM(I)=KX(J,I)                                                     
        ENDDO
        I=LM(1)                                                           
        JJ=LM(2)                                                          
        K=LM(3)                                                           
        L=LM(4)                                                           
        XY(1,1)=XX(I)                                                     
        XY(1,2)=XX(JJ)                                                    
        XY(1,3)=XX(K)                                                     
        XY(1,4)=XX(L)                                   
        XY(2,1)=YY(I)                                                     
        XY(2,2)=YY(JJ)                                                    
        XY(2,3)=YY(K)                                                     
        XY(2,4)=YY(L)                                   
        XCENT=(XX(I)+XX(JJ)+XX(K)+XX(L))/4000.d0
        YCENT=(YY(I)+YY(JJ)+YY(K)+YY(L))/4000.d0
C                                                                       
        CALL FESHAPE(1,CENTX,CENTY,PSI,DPSI)                         
      
C CALCULATE DXDS...EQUATION (5.3.6)                                     
        DO I=1,2                                                      
          DO L=1,2                                                      
            DXDS(I,L)=0.0d0
            DO K=1,NNODE                                                  
              DXDS(I,L)=DXDS(I,L)+DPSI(K,L)*XY(I,K)
            ENDDO
          ENDDO
        ENDDO
   
C CALCULATE DSDX...EQUATION (5.2.7)                                     
        DETJ=(DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1))                    
        IF (DETJ.LE.0.0) THEN
          WRITE(*,*) ' BAD JACOBIAN',DETJ,X                                             
          STOP                                                              
        ENDIF
        DSDX(1,1)=DXDS(2,2)/DETJ                                          
        DSDX(2,2)=DXDS(1,1)/DETJ                                          
        DSDX(1,2)=-DXDS(1,2)/DETJ                                         
        DSDX(2,1)=-DXDS(2,1)/DETJ                                         
C CALCULATE D(PSI)/DX...EQUATION (5.3.5)                                
        DO I=1,NNODE                                                  
          DPSIX(I)=DPSI(I,1)*DSDX(1,1)+DPSI(I,2)*DSDX(2,1)                  
          DPSIY(I)=DPSI(I,1)*DSDX(1,2)+DPSI(I,2)*DSDX(2,2)                  
        ENDDO                                                          
        DO I=1,NNODE                                                  
          SUMHX=SUMHX + HTICE(LM(I))*DPSIX(I)                                 
          SUMHY=SUMHY + HTICE(LM(I))*DPSIY(I)                                 
          SUMBX=SUMBX + DEPB(LM(I))*DPSIX(I)                                 
          SUMBY=SUMBY + DEPB(LM(I))*DPSIY(I)                                 
          HH=HH + (HTICE(LM(I))-DEPB(LM(I)))*PSI(I)                                       
        ENDDO                                                          
C                                                                       
        DELH=SUMHX**2 + SUMHY**2                                            
        DELH=SQRT(DELH)                                                   
        DELB=SUMBX**2 + SUMBY**2                                            
        DELB=SQRT(DELB)                                                   
        IF(HH.GT.HMIN) THEN                                               
          UX=-(SUMHX+0.09d0*SUMBX)                                            
          UY=-(SUMHY+0.09d0*SUMBY)                                           
        ELSE                                                              
          UX=0.d0
          UY=0.d0
        ENDIF                                                             
        UMAG=SQRT(UX**2+UY**2)                                        
        VMAX=MAX(VMAX,UMAG)                                               
        IF(UMAG.GT.VTHRESH .OR. .NOT.IGO) GOTO 600                                       
        IF(UMAG.EQ.0.) THEN                                               
          ICOLOR=1                                                        
        ELSE                                                              
          ICOLOR=1+NINT(UMAG/UDELTA)                                    
        ENDIF                                                             
        IF(ICOLOR.LT.2) ICOLOR=2                                          
        IF(ICOLOR.GT.14) ICOLOR=14                                        
        XARO(1)=XCENT                                                     
        YARO(1)=YCENT                                                     
        IF(UMAG.NE.0) THEN
          IF(INORM.EQ.1) THEN
            XARO(2)=XCENT+ASCAL*UX/UMAG
            YARO(2)=YCENT+ASCAL*UY/UMAG
          ELSE
            XARO(2)=XCENT+ASCAL*UX
            YARO(2)=YCENT+ASCAL*UY
          ENDIF
        ELSE
          XARO(2)=XCENT                                                     
          YARO(2)=YCENT                                                     
        ENDIF
        CALL LINCLR(ICMAP(ICOLOR))                                      
        ICSAV=ICOLOR                                                    
        CALL MOVE(REAL(XARO(1)),REAL(YARO(1)))
        CALL DRAW(REAL(XARO(2)),REAL(YARO(2)))
600   CONTINUE
      WRITE(*,*) 'GMAX=',VMAX                                           
      END
C===========================================
      SUBROUTINE TPROF(NUMNP,HTICE,DEPB,clear,IPLSTRT,TIME,
     &                 BMELT,WTHICK)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM,MMAX=40,LMAX=281)
      COMMON /TEMPERATURE/ TEMPA(MMAX,NMAX)
      COMMON /LINE/ NP,NLINE(1000)
      common /flush/ iflush
      logical iflush
      DIMENSION THICK(NMAX),depb(nmax),HTICE(nmax),BMELT(NMAX)
      DIMENSION WTHICK(NMAX)
      DIMENSION LMAP(MMAX),xxx(mmax),tttt(mmax)
      dimension ilab(8)
      data ilab /1,30,61,92,123,154,185,215/
      logical clear
      dimension xp(lmax),yp(lmax,mmax),tp(lmax,mmax)
      dimension bm(lmax),wt(lmax)
      dimension icol(nmax,mmax),xy(3,4),ic(4),x2(2)
      real*4 px(4),py(4)
      character*14 junk
      save itype
      data itype /1/
      do i=1,numnp
        THICK(i)=HTICE(i)-DEPB(i)
c        print *,i,htice(i),depb(i)
      enddo
      if(IPLSTRT.eq.0 ) then
        call grstrt(600,600)
        PRINT *,'TYPE:0-profile,1:color contour',itype
        if(iflush) call gflush
        READ(*,*) itype
        WRITE(99,*) itype
      endif
c old way profiles ...
      if(itype.eq.0) then
        if(.false.) then
          xmin=0.0d0
          xmax=dble(np+1)
          ymin=1d30
          ymax=-ymin
          do nn=1,np
            n=nline(nn)
            ymin=min(ymin,depb(n))
            ymax=max(ymax,depb(n)+thick(n))
c           print *,nn,n,real(depb(n)),real(thick(n))
          enddo
c         print *,xmin,xmax
c         print *,ymin,ymax
          call newpag       
          call window(real(-xmax),real(xmax),
     &                real(1.5d0*ymin),real(1.5d0*ymax))
          do nn=1,np
            call linclr(1)
            n=nline(1)
            call move(real(xmin),0.0)
            call draw(real(xmax),0.0)
            call move(1.0,real(depb(n)))
            do mm=2,np
              n=nline(mm)
              call draw(real(mm),real(depb(n)))
            enddo
            n=nline(1)
            call move(1.0,real(depb(n)+thick(n)))
            do mm=2,np
              n=nline(mm)
              call draw(real(mm),real(depb(n)+thick(n)))
            enddo
            call linclr(2)  
            n=nline(nn)
            LTOT=0
            DO I=1,MMAX/2
              LMAP(I)=I
              LTOT=LTOT+LMAP(I)
            ENDDO
            DO I=MMAX/2+1,MMAX
              LMAP(I)=LMAP(I-1)-1
              LTOT=LTOT+LMAP(I)
            ENDDO
            TMELT=PMP(THICK(n))
            TBASE=nn
            XXX(1)=THICK(n)+depb(n)
            XXX(MMAX)=depb(n)
            call move(real(TBASE),real(xxx(1)))
            call draw(real(TBASE),real(xxx(MMAX)))
            DX=(XXX(MMAX)-XXX(1))/dble(LTOT)
            DO I=2,MMAX
              XXX(I)=XXX(I-1)+LMAP(I-1)*DX
            ENDDO
            DO I=1,MMAX
              XXX(I)=XXX(I)
            ENDDO
            DO I=1,MMAX
              XMAX=MAX(XMAX,XXX(I))
              XMIN=MIN(XMIN,XXX(I))
              TTTT(i)=dscale*TEMPA(i,n)+TBASE
            enddo
            CALL LINCLR(1)                                             
            CALL MOVE(REAL(1*(TTTT(1)-TMELT)),REAL(XXX(1)))                           
            DO I=2,MMAX                                                  
              CALL DRAW(REAL(1*(TTTT(I)-TMELT)),REAL(XXX(I)))                         
            ENDDO  
            call linclr(3)
            call draw(real(TBASE),real(xxx(MMAX)))
            if(iflush) call gflush
c           pause   
c           call newpag       
          enddo
        endif
        if(.true.) then
          xmin=0.0
          xmax=dble(np+1)
          ymin=1d30
          ymax=-ymin
          do nn=1,np
            n=nline(nn)
            ymin=min(ymin,depb(n))
            ymax=max(ymax,depb(n)+thick(n))
c           print *,nn,n,real(depb(n)),real(thick(n))
          enddo
c         print *,xmin,xmax
c         print *,ymin,ymax
          delx=(xmax-xmin)/10.
          dely=(ymax-ymin)/10.
          call newpag       
          call window(real(-xmax-delx),real(xmax+delx),
     &                real(ymin-dely),real(ymax+dely))
          dmin=1e30
          dmax=-1e30
          do nn=1,np
            n=nline(nn)
            LTOT=0
            DO I=1,MMAX/2
              LMAP(I)=I
              LTOT=LTOT+LMAP(I)
            ENDDO
            DO I=MMAX/2+1,MMAX
              LMAP(I)=LMAP(I-1)-1
              LTOT=LTOT+LMAP(I)
            ENDDO
            TMELT=PMP(THICK(n))
            TBASE=nn
            XXX(1)=THICK(n)+depb(n)
            XXX(MMAX)=depb(n)
c            dmin=min(dmin,TBASE)
c            dmax=max(dmax,TBASE)
            DX=(XXX(MMAX)-XXX(1))/dble(LTOT)
            DO I=2,MMAX
              XXX(I)=XXX(I-1)+LMAP(I-1)*DX
            ENDDO
            DO I=1,MMAX
              XXX(I)=XXX(I)
            ENDDO
            DO I=1,MMAX
              TTTT(i)=TEMPA(i,n)+TBASE
              TTTT(i)=TEMPA(i,n)-TMELT
              dmin=min(dmin,TTTT(I))
              dmax=max(dmax,TTTT(I))
            enddo
c            DO I=1,MMAX                                                  
c              dmin=min(dmin,TTTT(I)-TMELT)
c              dmax=max(dmax,TTTT(I)-TMELT)
c            ENDDO  
c           pause   
          enddo
          dmax=0
          dmin=-40
          if(dmax-dmin.ne.0) then
            dscale=xmax/(dmax-dmin)
          else
            dscale=1
          endif
c          write(7,*) 'dscale',dscale,dmin,dmax,xmax
          do nn=1,np
            call linclr(1)
            n=nline(1)
            call move(real(xmin),0.0)
            call draw(real(xmax),0.0)
            call move(1.0,real(depb(n)))
            do mm=2,np
              n=nline(mm)
              call draw(real(mm),real(depb(n)))
            enddo
            n=nline(1)
            call move(1.0,real(depb(n)+thick(n)))
            do mm=2,np
              n=nline(mm)
              call draw(real(mm),real(depb(n)+thick(n)))
            enddo
            call linclr(2)  
            n=nline(nn)
            LTOT=0
            DO I=1,MMAX/2
              LMAP(I)=I
              LTOT=LTOT+LMAP(I)
            ENDDO
            DO I=MMAX/2+1,MMAX
              LMAP(I)=LMAP(I-1)-1
              LTOT=LTOT+LMAP(I)
            ENDDO
            TMELT=PMP(THICK(n))
            TBASE=nn
            XXX(1)=THICK(n)+depb(n)
            XXX(MMAX)=depb(n)
            call move(real(TBASE),real(xxx(1)))
            call draw(real(TBASE),real(xxx(MMAX)))
            DX=(XXX(MMAX)-XXX(1))/dble(LTOT)
            DO I=2,MMAX
              XXX(I)=XXX(I-1)+LMAP(I-1)*DX
            ENDDO
            DO I=1,MMAX
              XXX(I)=XXX(I)
            ENDDO
            DO I=1,MMAX
              XMAX=MAX(XMAX,XXX(I))
              XMIN=MIN(XMIN,XXX(I))
              TTTT(i)=TEMPA(i,n)-TMELT
              TTTT(i)=dscale*TTTT(i)+TBASE
           enddo
            CALL LINCLR(1)                                             
            CALL MOVE(REAL(1*(TTTT(1))),REAL(XXX(1)))                           
            DO I=2,MMAX                                                  
              CALL DRAW(REAL(1*(TTTT(I))),REAL(XXX(I)))                         
            ENDDO  
            call linclr(3)
            call draw(real(TBASE),real(xxx(MMAX)))
            if(iflush) call gflush
c           pause   
c           call newpag       
          enddo
        endif
c new way color contours ...
      else
        ncolor=215
        xmax=-1e30
        xmin=1e30
        ymax=-1e30
        ymin=1e30
        tmax=-1e30
        tmin=1e30
        imax=0
        imin=ncolor
        do nn=1,np
          n=nline(nn)
          LTOT=0
          DO I=1,MMAX/2
            LMAP(I)=I
            LTOT=LTOT+LMAP(I)
          ENDDO
          DO I=MMAX/2+1,MMAX
            LMAP(I)=LMAP(I-1)-1
            LTOT=LTOT+LMAP(I)
          ENDDO
          TMELT=PMP(THICK(n))
          TBASE=0.0
          xp(nn)=nn
          xp(nn)=nn
          yp(nn,1)=THICK(n)+depb(n)
          yp(nn,MMAX)=depb(n)
          dy=(yp(nn,MMAX)-yp(nn,1))/dble(LTOT)
          DO I=2,MMAX
            yp(nn,I)=yp(nn,I-1)+LMAP(I-1)*dy
          ENDDO
          xmax=MAX(xmax,xp(nn))
          xmin=MIN(xmin,xp(nn))
          DO I=1,MMAX
            tp(nn,i)=TEMPA(i,n)-TBASE
            ymax=MAX(ymax,yp(nn,I))
            ymin=MIN(ymin,yp(nn,I))
            tmax=MAX(tmax,tp(nn,I))
            tmin=MIN(tmin,tp(nn,I))
          enddo
          tmax=0.0
          tmin=-49.0
          tmin=-110.0
          ymin=-1000.
          ymax=19000.
          tdelt=(tmax-tmin)/(ncolor-1)
          do j=1,mmax
            if(tp(nn,j).eq.-99999.) then
              icol(nn,j)=-99999
            else
              icol(nn,j)=izset(ncolor,tp(nn,j),tmin,tdelt,0)
              imin=min(imin,icol(nn,j))
              imax=max(imax,icol(nn,j))
            endif
          enddo
        enddo
        ddx=(xmax-xmin)/10
        call newpag       
        call window(real(xmin),real(xmax),
     &              real(ymin),real(ymax))

          
        dx=(xmax-xmin-2*ddx)/ncolor
        dy=0.05*(ymax-ymin)
        ybox=ymax-dy
        call filpan(1,.false.)
        do i=1,ncolor
          xbox=xmin+(i-1)*dx+ddx
            px(1)=xbox
            py(1)=ybox

            px(2)=xbox+dx
            py(2)=ybox

            px(3)=xbox+dx
            py(3)=ybox+dy

            px(4)=xbox
            py(4)=ybox+dy

            call linclr1(i)
            call panel(4,px,py)
        enddo
        ybox=ymax-2*dy
        do ii=1,8
          i=ilab(ii)
          xbox=xmin+(i-1)*dx+ddx/2
          rnum=tmin+(i-1)*tdelt
          rnum=real(nint(rnum))
          icc=izset(ncolor,rnum,tmin,tdelt,0)    
          write(junk,'(f8.3)') rnum
c          write(*,'(i4,f8.3,i4)') i,rnum,izset(ncolor,rnum,zmin,zdelt,0)
          call linclr1(i)
          call move(real(xbox),real(ybox))
          call text(14,junk)
        enddo

        do nn=1,np-1
          px(1)=xp(nn)
          px(2)=xp(nn)
          px(3)=xp(nn+1)
          px(4)=xp(nn+1)
          do ll=1,MMAX-1
            py(1)=yp(nn,ll)
            py(2)=yp(nn,ll+1)
            py(3)=yp(nn+1,ll+1)
            py(4)=yp(nn+1,ll)
            icc=nint((icol(nn,ll)+icol(nn,ll+1)+
     &          icol(nn+1,ll+1)+icol(nn+1,ll))/4.)
            call linclr1(icc)
            call panel(4,px,py)
          enddo
        enddo
        write(junk,'(a5,f8.1)') 'TIME=',TIME
        call linclr(1)
        ybox=ymax-3*dy
        xbox=xmin+ddx/2
        call move(real(xbox),real(ybox))
        call text(14,junk)

        call linclr(1)
        n=nline(1)
        call move(real(xmin),0.0)
        call draw(real(xmax),0.0)
        call move(1.0,real(depb(n)))
        do mm=2,np
          n=nline(mm)
          call draw(real(mm),real(depb(n)))
        enddo
        n=nline(1)
        call move(1.0,real(depb(n)+thick(n)))
        do mm=2,np
          n=nline(mm)
          call draw(real(mm),real(depb(n)+thick(n)))
        enddo
        pllev=-3000.
        call move(real(xmin),real(pllev))
        call draw(real(xmax),real(pllev))
        if(.true.) then
          n=nline(1)
          do mm=1,np
            n=nline(mm)
            call linclr(1)
            call move(real(mm-0.1),real(pllev))
            down=pllev+bmelt(n)*1d5
            call draw(real(mm-0.1),real(down))
            call linclr(2)
            call move(real(mm+0.1),real(pllev))
            down=pllev+wthick(n)*1d4
            call draw(real(mm+0.1),real(down))
            bm(mm)=bmelt(n)
            wt(mm)=wthick(n)
          enddo
        else
          n=nline(1)
          down=real(pllev)+bmelt(n)*1d6
          call move(real(1),real(down))
          do mm=1,np
            n=nline(mm)
            down=real(pllev)+bmelt(n)*1d6
            call draw(real(mm),real(down))
          enddo
        endif
!       if(.false.) then
!         write(37) junk
!         write(37) np,mmax
!         write(37) (xp(i),i=1,np)
!         write(37) ((yp(i,j),i=1,np),j=1,mmax)
!         write(37) ((tp(i,j),i=1,np),j=1,mmax)
c         write(37) (bm(i),i=1,np)
c         write(37) (wt(i),i=1,np)
!       else
!         do i=1,np
!           write(37,*) xp(i),yp(i,j),tp(i,j)
!         enddo
!       endif
      endif

c      if(clear) pause
c      if(IPLSTRT.eq.0 ) call grstop1
      end                                           
