      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM,NZ=9,NZ1=NZ+1,NSMAX=MAXTIME,N3=NMAX*3,
     &           N23=N3*2)
C ***************************************************************C
C                                                                C
C   PROGRAM:  MAP5                                               C
C                                                                C
C   DATE:  7-31-95                                               C
C   PROGRAMMER:  J.  FASTOOK                                     C
C                                                                C
C   FUNCTION:                                                    C
C           THIS IS A PROGRAM TO MODEL THE FLOW OF A GLACIER     C
C           WITH MATERIAL PROPERTIES READ FOR EACH NODAL POINT   C
C           AND THE AVERAGE USED FOR THE ELEMENT.                C
C           IT CALCULATES AN ELEMENT MATRIX ALA BECKER, 2-D      C
C           PROGRAM WITH LINEAR SHAPE FUNCTIONS.                 C
C           DOES TIME DEPENDENT CASE USING A LUMPED CAPACITANCE  C
C           MATRIX AND A BACKWARD DIFFERENCE SCHEME. ALLOWING    C
C           FOR VERY QUICK OPERATION.                            C
C ***************************************************************C
C     PROGRAM FOR STEADY AND UNSTEADY STATE FLOW ANALYSIS
C     USING FINITE ELEMENTS
C
C ***************************************************************C
      DIMENSION ZERO(NMAX)
      data zero /nmax*0.d0/
      DIMENSION WWW(N3),THICKL(NMAX),WRATE(N3,2)
      DIMENSION WWWORIG(NMAX),WDIFF(NMAX)
      CHARACTER HED*80,SCRTCH*80,IADJ*2
      REAL*8 A(NMAX,NZ),CAP(NMAX,NZ)
      COMMON /MAT/  A,KA(NMAX,NZ1),NUMNP
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      COMMON /VELOS/ ASCAL,UMAX,VTHRESH,INORM
      DIMENSION IFIT(7), AMASS(11), B(NMAX), X(NMAX), Y(NMAX), D(NMAX),
     &          BOLD(NMAX),FLUX(NMAX),ADIAG(NMAX),ITKODE(NMAX),
     &          KODE(NMAX), CONST(NMAX), ACON(NMAX), LM(5), ITYPE(NMAX),
     &          KX(NMAX,4), IDT(NMAX), ADOTB(NMAX), 
     &          ADOT(NMAX), BDROCK(NMAX), FLOWA(NMAX), SLDGB(NMAX),
     &          PSURF(NMAX), PPSURF(NMAX), FRACT(NMAX),
     &          CNEW(NMAX), QHOLD(NMAX), HTICE(NMAX), THICK(NMAX),
     &          HFIT(NMAX), IBFLUX(NMAX,2), BFLUX(NMAX), DEPB(NMAX),
     &          QQQ(NMAX),TEMP(NMAX),TBED(NMAX),VEL(NMAX,3),
     &          UNDEPB(NMAX), SLOPE(4,NMAX), SLOPN(4,NMAX), KZ(NMAX),
     &          WTHICK(NMAX),WTELEM(NMAX),GEOFLUX(NMAX),
     &          CALV(NMAX),PCALV(NMAX),
     &          WVELX(NMAX),WVELY(NMAX),CONTRIB(NMAX)
      DIMENSION AFUDGE(NMAX),ADP(NMAX),ADM(NMAX),ADC(NMAX),BMELT(NMAX)
      DIMENSION TTIME(NSMAX),VVOL(NSMAX),AAREA(NSMAX),TTBOT(NSMAX),
     &          TTNSL(NSMAX),TWATER(NSMAX),PWATER(NSMAX),WWWMIN(NSMAX),
     &          TTSEAL(NSMAX)
      DIMENSION TTAVG(NSMAX)
      DIMENSION NTYPE(NMAX), NNODE(NMAX),
     &          XI(2,9), ETA(2,9), W(2,9),
     &          AADOT(NMAX), AFRACT(NMAX), AFLOWA(NMAX),
     &          ABDRCK(NMAX), ASLDGB(NMAX)
      DIMENSION ALPHAC(3)
C     DIMENSION DPSIX(9), DPSIY(9), DXDS(2,2), DSDX(2,2),
C    &          PSI(4), DPSI(4,2), CNST(NMAX),
C    &          XY(2,4)
      DIMENSION XL(NMAX),YL(NMAX),FLUXL(NMAX)
      DIMENSION ACC(NMAX),ABLAT(NMAX)
      COMMON /COORD/ CLAT(NMAX),CLONG(NMAX)
      LOGICAL HIRES
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      LOGICAL ZFLAG
      COMMON /ZOOMCOM/ PXMIN,PYMIN,PXMAX,PYMAX,ZFLAG
      DATA ZFLAG /.FALSE./
      DATA PXMIN,PYMIN,PXMAX,PYMAX /4*0.D0/
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      LOGICAL IFLUSH
      COMMON /FLUSH/ IFLUSH
      COMMON /HEAT/ HTOGG
      LOGICAL HTOGG
      DATA HTOGG /.TRUE./
      COMMON /IOLIST/ LIST
C ... TIMER STUFF, FOR SGI ONLY ...
      REAL*4 TB(2)
c     REAL*4 TB(2),ETIME,DTIME
c     EXTERNAL ETIME,DTIME
C ... TIMER STUFF, FOR SGI ONLY ...
      REAL*8 TTT(NMAX),TNEW(NMAX)
      EXTERNAL WARMING
      DATA HFIT /NMAX*0.0/
      DATA ASCAL /1.D0/, UMAX /100.D0/, VTHRESH /1000.D0/
      DATA INORM /1/
      DATA WWW /N3*0.D0/
      DATA WRATE /N23*0.D0/
      DATA IFLUSH /.true./
      DATA RHOI /0.917D0/, RHOW /1.092D0/, RHOR /4.0D0/
      DATA ADVANCE /0.95D0/
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      logical file09
      character*41 string
      data IPLOT /0/
      COMMON /TCONST/ TSORIG(NMAX),TFLAG
      LOGICAL TFLAG
      logical stopflag
      common /stop/ stopflag
      stopflag=.false.
      open(2,file='plot-junk.data')
      open(3,file='plot-temps.data')
      string='(1x,39(/))'
123   FORMAT(A25,T30,1PG13.6,G13.6)
      HIRES=.true.
      RATDEN=RHOW/RHOI
      IPASS=0
      ICON=0
      CALL SETRIG
      BIG=1D20
      ITOGG=.FALSE.
      TFRACT=1
      TFLAG =.true.
c     TFLAG =.false.
c     TFLAG =ITOGG    !  couple fixed temp flag with augmnted acc toggle
      SLTOGG=0
      CTOGG=.TRUE.
      WTOGG=0
      BTOGG=1
      IOTOGG=.TRUE.
      XPOLE=0.0D0
      YPOLE=0.0D0
      ACOM=-9.6237669D0
      WINDIR(1)=90.D0/180.D0*3.14159D0
      WINDIR(2)=100.D0
C     WINDIR(1)=315.D0/180.D0*3.14159D0
C     WINDIR(2)=50.D0
C FOLLOWING IS DEFAULT SNOWLINE ELEVATION AT POLE, GRADIENT, TNSL
      AMASS(1)=-0.1d0
      AMASS(2)=0.1d0
      AMASS(3)=0.01d0
      AMASS(4)=5000.D0
      AMASS(5)=6000.D0
      AMASS(6)=7000.D0
      AMASS(7)=-3000.D0
      AMASS(8)=.56D-3
      AMASS(8)=.001D0
      AMASS(9)=-14.D0
      AMASS(10)=.4D0
      AMASS(11)=-10000.D0                     ! SEALEV
      CFACTOR=0.D0
C A(1): ADVECTION COEFFICIENT, A(2): DIFFUSION, A(3): LOSS TERM
      ALPHAC(1)=1.D0
      ALPHAC(2)=1.D0
      ALPHAC(3)=1D-2
      TBASE=0.0D0
      TNSLBASE=14.D0
C*******************
C ... CHECK FO DEFAULTS FILE. IF NONE EXIST, CREATE ONE ...
      inquire(file='fort.9',exist=file09)
      if(.not.file09) goto 100
      READ(9,*,END=100) AMASS,ACOM,WJUNK,WINDIR(2),
     &                  XPOLE,YPOLE,CFACTOR,ALPHAC,TBASE,
     &                  TNSLBASE,CTOGG,WTOGG,ITOGG,BTOGG
c     TFLAG =ITOGG    !  couple fixed temp flag with augmnted acc toggle
C                                T   (0,1,2,3)   F  (1,2,3)
        WINDIR(1)=WJUNK/180.D0*3.14159D0
        PRINT *,'**** DEFAULTS FILE FOUND ****'
        IF(AMASS(11).EQ.1.2D0) THEN
          PRINT *,' OLD FILE, CHANGE ',AMASS(11)
          PAUSE
        ENDIF
        GOTO 101
100   PRINT *,'NO DEFAULTS FILE FOUND'
        REWIND 9
        WJUNK=WINDIR(1)*180.D0/3.14159D0   
        WRITE(9,*) AMASS,ACOM,WJUNK,WINDIR(2),
     &             XPOLE,YPOLE,CFACTOR,ALPHAC,TBASE,
     &             TNSLBASE,CTOGG,WTOGG,ITOGG,BTOGG
101   CONTINUE
      PRINT *,' TNSL       = ',AMASS(9)
      PRINT *,' FUDGE      = ',AMASS(10)
      PRINT *,' SEA LEV    = ',AMASS(11)
      PRINT *,' ACOM       = ',ACOM
      PRINT *,' WINDIR     = ',WJUNK,WINDIR(2)
      PRINT *,' XPOLE,YPOLE= ',XPOLE,YPOLE
      PRINT *,' TIME       = ',TBASE
      PRINT *,' TNSLBASE   = ',TNSLBASE
      SEALEV=AMASS(11)
C*******************
      DO I=1,6
        IFIT(I)=0
      ENDDO
      IFIT(7)=1
      NTSTEP=0
C
C ... FOLLOWING IS SEALEVEL REFERENCED TO PRESENT=0.
      SEALEV=-10000.D0
      SEALEV=AMASS(11)
C
C ... FOLLOWING SETS RATE OF CONVERGENCE, UP TO 5 WORKS WELL
      CONV=1.0D0
      TIME=TBASE
      PG = 0.089866D0*0.3816d0
      NCOL=NZ
C
C ... INITIALIZE INTEGRATION POINTS AND WEIGHTS
C ... GAUSSIAN QUADRATURE OF ORDER THREE QUADRILATERALS
      CALL GAUSINIT(XI,ETA,W)
C
      DO I=1,4
        LM(I)=0
      ENDDO
C
C ... FOLLOWING READN FOR SPLIT DATA SETS
      CALL READN(NMAX,HED,NUMNP,NUMEL,NUMGBC,NDT,INTER,DT,
     &           KODE,X,Y,HTICE,ADOT,ADOTB,FRACT,PSURF,RHOI,RHOW,
     &           RHOR,DEPB,BDROCK,
     &           UNDEPB,FLOWA,ACON,SLDGB,TEMP,ITYPE,AFUDGE,GEOFLUX,
     &           THICK,KX,CONST,IBFLUX,BFLUX,QHOLD,NTYPE,NNODE,NCOL,
     &           AADOT,AFRACT,ABDRCK,PPSURF,AFLOWA,ASLDGB,IDT,AMASS,
     &           NUMCOL,NUMLEV,CALV,PCALV,WWW,WRATE,THICKL,
     &           WWWORIG,TSORIG,TIME,ACC,ABLAT)
      CALL READCOORD(NUMNP)
C ..... CALCULATES VOLUMES (FLOTATION AND TOTAL) AND AREAL EXTENT
        CALL VOLUME(NMAX, TIME, NUMNP, NUMEL, X, Y, KX, PSURF, BDROCK,
     &             DEPB, ADOT, RHOI, RHOW, VOL, AREA, AMASS,
     &             GRIDAREA,.FALSE.,SEABASE)
      ACOMSAVE=ACOM
      IF(.FALSE.) THEN
        WRITE(7,123) ' TIME BEFORE TEMPER ',ETIME(TB),DTIME(TB)
        CALL RESETTEMP(NUMNP,AMASS(9),HTICE,TEMP,AMASS,TIME)
        CALL TEMPER(NMAX, NUMNP, NUMEL, NTYPE, KX, KODE, ADOT, FRACT,
     &              DEPB, HTICE,PSURF,FLOWA, SLDGB, TEMP, ITYPE, TBED,
     &              X, Y, AMASS, SLOPN, TIME, DT, IPLOT, TBAVG,
     &              AFUDGE,GRIDAREA,IMELT,BMELT,WTHICK,VOL,AREA,
     &              GEOFLUX,VEL,IDT,CONTRIB)
        WRITE(7,123) ' TIME AFTER TEMPER ',ETIME(TB),DTIME(TB)
      ENDIF
      II=1
      DO I=1,NUMNP
        WDIFF(I)=WWW(II)-WWWORIG(I)
        II=II+3
      ENDDO
C......EXPERIMENTAL SECTION FOR BED DEPRESSION ..................!
C ... THIS GENERATES EQUILIBRUIUM DEPRESSION, DELT=0.0           !
      IF(.FALSE.) THEN                                           !
C     IF(BTOGG.EQ.2 .OR. BTOGG.EQ.3) THEN                                             !
        DO I=1,NUMNP                                             !
          THICKL(I)=HTICE(I)-BDROCK(I)                           !
          FLOT=(1.D0-RATDEN)*(BDROCK(I)-SEALEV)                        !
          IF(HTICE(I).LT.FLOT) THEN                              !
            THICKL(I)=SEALEV                                       !
          ENDIF                                                  !
        ENDDO                                                    !
C        CALL PLATE(NTSTEP,NMAX,N3,NUMNP,NUMEL,X,Y,KX,THICKL,    !
C     &             0.D0,WWW,WRATE,WMIN,TIME,WWWORIG)            !
        CALL EPLATE(NTSTEP,NUMNP,NUMEL,X,Y,KX,THICKL,KODE,       !
     &             0.D0,WWW,WRATE,WMIN,TIME,WWWORIG,FNET,        !
     &                 WWW)                                      !
        DO I=1,NUMNP                                             !
C         HTICE(I)=HTICE(I)+WWW((I-1)*3+1)                       !
          DEPB(I)=UNDEPB(I)+WWW((I-1)*3+1)                       !
          IF(HTICE(I).LT.DEPB(I)) HTICE(I)=DEPB(I)               !
        ENDDO                                                    !
C        CALL POUTSTF(NUMNP*3,WWW,WRATE(1,1),NUMNP,THICKL,NUMCOL, !
C     &               NUMLEV,TIME,WDIFF)                                !
      ENDIF                                                      !
C................................................................!
      DTLOCAL=DT
      DO I=1,NUMNP
C THIS WILL NEED TO BE READ AND WRITTEN EVENTUALLY...
        IMELT=1
C        BMELT(I)=-.0001
        BMELT(I)=0D0
        WTHICK(I)=-1D-6
        WTHICK(I)=0D0
c------experimental-----------
c       HTICE(I)=DEPB(I)
        QQQ(I)=HTICE(I)
        TTT(I)=HTICE(I)
        ADP(I)=0.0D0
        ADM(I)=0.0D0
        ADC(I)=0.0D0
      ENDDO
C
C ... SET LINEARIZATION CONSTANT USING INITIAL CONFIGURATION
       WRITE(7,123) ' TIME BEFORE NCONST ',ETIME(TB),DTIME(TB)
       CALL NCONST(NMAX, X, Y, KX, NTYPE, NUMEL,
     &             AFRACT, ASLDGB, LM, AFLOWA, BDROCK, DEPB,
     &             UNDEPB, PG, QQQ, CNEW, SLOPE, RHOI, WINDIR,
     &             WTHICK,ITYPE)
       WRITE(7,123) ' TIME AFTER NCONST ',ETIME(TB),DTIME(TB)
C$DOACROSS LOCAL(I)
      IF(.false.) THEN
C ..... THIS DISCARDS ANY VALUES THAT HAVE BEEN READ IN...
        DO I=1,NUMEL
          CONST(I)=CNEW(I)
        ENDDO
      ENDIF
C ....... PRINT OUT INITIAL CONFIGURATION .................
          IF(.TRUE.) THEN
            WRITE(*,*)
            WRITE(*,*) 'DUMPING INITIAL CONFIGURATION . . .'
            WRITE(*,*)
            WRITE(SCRTCH,*) 'TIME=',NINT(TIME)
            WRITE(34) SCRTCH
            WRITE(34) (HTICE(I),I=1,NUMNP)
            WRITE(34) (ADOT(I),I=1,NUMNP)
            WRITE(34) (DEPB(I),I=1,NUMNP)
            WRITE(34) (CONST(I),I=1,NUMEL)
            WRITE(34) (ACON(I),I=1,NUMEL)
            WRITE(35) (FRACT(I),I=1,NUMNP)
            WRITE(35) (FLOWA(I),I=1,NUMNP)
            WRITE(35) (SLDGB(I),I=1,NUMNP)
            WRITE(35) (AFUDGE(I),I=1,NUMNP)
            WRITE(36) SCRTCH
            WRITE(36) (TBED(I),I=1,NUMNP)
            WRITE(36) (BMELT(I),I=1,NUMNP)
            WRITE(36) (WTHICK(I),I=1,NUMNP)
            WRITE(36) (CONTRIB(I),I=1,NUMNP)
          ENDIF
C
      WRITE(*,*) 'TIME STEP=',DT
      WRITE(*,*) 'INPUT 1 TO ENTER ADJUST, 0 TO BYPASS'
      IF(IFLUSH) CALL GFLUSH
      READ(*,4000) IADJ
      WRITE(99,4000) IADJ
      IF(IADJ.EQ.'1') THEN
C
C ..... CALCULATE SLOPES IN CASE NEEDED BY ADJUST
        CALL NODESL(NMAX, NUMNP, NUMEL, KX, SLOPE, SLOPN)
C
C ..... ENTER INTERACTIVE DATA SET MANIPULATOR
        CALL ADJUST(HED, NUMNP, NUMEL, X, Y, HTICE, ADOT, 
     &              ADOTB,FRACT,BMELT,WTHICK, 
     &              TEMP, ITYPE, TBED, PSURF, UNDEPB,
     &              BDROCK, DEPB, FLOWA, SLDGB, THICK, KX, CONST, 
     &              AFUDGE,NNODE, KODE, GEOFLUX, FLUX,
     &              HFIT, NUMCOL, NUMLEV, NUMGBC, NDT, INTER, DT,
     &              IBFLUX, BFLUX, NMAX, IDT, SLOPN, AMASS, TIME,
     &              NTSTEP, TTIME, VVOL, AAREA,TTBOT,TTAVG, TTNSL,
     &              TTSEAL, IFIT,
     &              IPLOT, CFACTOR, ACON,ICON, 
     &              ALPHAC,TBASE,
     &              IMELT,TWATER,PWATER,CALV,
     &              NTYPE,AADOT,AFRACT,ABDRCK,PPSURF,
     &              AFLOWA,ASLDGB,PCALV,ADC,WWWMIN,WWW,WRATE,WDIFF,
     &              WWWORIG,ADVANCE,HIRES,ACC,ABLAT,TFRACT)
        DTLOCAL=DT
C
      ENDIF
C$DOACROSS LOCAL(I)
c      DO I=1,NUMNP
c        QQQ(I)=HTICE(I)
c        TTT(I)=HTICE(I)
c      ENDDO
C
C ... LOAD NODAL PROPERTIES INTO ELEMENT PROPERTIES
      CALL CALVING(NMAX, NUMEL, NTYPE, KX, HTICE, DEPB, ADOT, 
     &           AADOT,CFACTOR,ADC,CALV,PCALV,AMASS(9))
      CALL ELPROP(NMAX, NUMEL, NTYPE, KX, ADOT, 
     &             AADOT, FRACT,
     &             AFRACT, BDROCK, ABDRCK, PSURF, PPSURF, FLOWA,ACON,
     &             ICON, AFLOWA, SLDGB, ASLDGB,CALV,PCALV)
C
C
C
C
C
C .....................*** MAIN LOOP ***............................
C
C
C
215   CONTINUE
C
C 
C ..... FOLLOWING CALCULATES NODE SLOPE FROM ELEMENT SLOPE
        CALL NODESL(NMAX, NUMNP, NUMEL, KX, SLOPE, SLOPN)
C
C ..... FOLLOWING ADJUSTS ADOT FOR FITTED ACCUMULATION
        CALL ACCSET(NUMNP,AMASS,ACOMSAVE,TIME,
     &                  THICK,IDT,HTICE,BDROCK,PSURF,SLOPN,X,Y,TEMP,
     &                  ADOT,ADOTB,ADP,ADM,TJUNK,TFRACT,ACC,ABLAT)
        CALL RESETTEMP(NUMNP,AMASS(9),HTICE,TEMP,AMASS,TIME)

C---------------------------------------------------------------------
      IF(.FALSE.) THEN
C**********EISMINT ALA HUYBRECHTS ********************
        WRM=WARMING(AMASS(9)+TNSLBASE)
C*****************************************************
C$DOACROSS LOCAL(I,ADOTI,AJUNK), SHARED(TIME)
        DO I=1,NUMNP
C          IF(THICK(I).GT.0) THEN
C            ACOM=ACOMSAVE
C          ELSE
C            ACOM=ACOMSAVE/3.
C          ENDIF
          ADOTI=AFUNCT(TIME, IDT(I), AMASS,
     &                   HTICE(I), BDROCK(I), PSURF(I),
     &                   SLOPN(1,I),X(I),Y(I),TEMP(I))
          IF(IDT(I).GT.0) THEN
            ADOT(I)=ADOTI
          ELSE
            IF(ITOGG) THEN
              AJUNK=ACCUM(IDT(I),X(I)*.001D0,Y(I)*.001D0,HTICE(I),
     &                    SLOPN(1,I),
     &                    PSURF(I),AMASS(9),TEMP(I))
C             ADOT(I)=ADOTB(I)-ABL*.01
C*****************************************************
C**********EISMINT ALA HUYBRECHTS ********************
              ADP(I)=WRM*ADOTB(I)
              ADM(I)=-ABL*0.01D0
C STRAIGHT
              ADOT(I)=WRM*ADOTB(I)-ABL*.01D0
C SLOW DOWN ACCUM RATE CHANGE...-------------------------
C              ADOT(I)=0.5*(ADOT(I)+WRM*ADOTB(I)-ABL*.01)
C--------------------------------------------------------
C*****************************************************
            ELSE
              ADOT(I)=ADOTB(I)
            ENDIF
          ENDIF
C ......  FOLLOWING ZEROS ACCUM FOR STEEP SLOPE, NO ICE...
          CALL STEEPSET(I,HTICE(I),BDROCK(I),SLOPN(1,I),ADOT(I))
        ENDDO
      ENDIF
C---------------------------------------------------------------------

C
        IF(IPASS.EQ.0) THEN
          CALL DVELO(NMAX,NUMEL,KX,X,Y,HTICE,DEPB,CONST,
     &             VEL,.FALSE.,DTMIN)
          IPASS=1
          WRITE(7,123) ' TIME BEFORE TEMPER ',ETIME(TB),DTIME(TB)
          CALL RESETTEMP(NUMNP,AMASS(9),HTICE,TEMP,AMASS,TIME)
          CALL TEMPER(NMAX, NUMNP, NUMEL, NTYPE, KX, KODE, ADOT, FRACT,
     &                DEPB, HTICE,PSURF,FLOWA, SLDGB, TEMP, ITYPE, TBED,
     &                X, Y, AMASS, SLOPN, TIME, DT, IPLOT, TBAVG,
     &                AFUDGE,GRIDAREA,IMELT,BMELT,WTHICK,VOL,AREA,
     &                GEOFLUX,VEL,IDT,CONTRIB)
          WRITE(7,123) ' TIME AFTER TEMPER ',ETIME(TB),DTIME(TB)
        ENDIF
C ..... LOADS NEW NODAL MATERIAL PROPERTIES INTO ELEMENT MATERIAL PROPERTIES
        CALL CALVING(NMAX, NUMEL, NTYPE, KX, HTICE, DEPB, ADOT, 
     &           AADOT,CFACTOR,ADC,CALV,PCALV,AMASS(9))
        CALL ELPROP(NMAX, NUMEL, NTYPE, KX, ADOT, 
     &             AADOT, FRACT,
     &             AFRACT, BDROCK, ABDRCK, PSURF, PPSURF, FLOWA,ACON,
     &             ICON, AFLOWA, SLDGB, ASLDGB,CALV,PCALV)
C
C ..... CALCULATES VOLUMES (FLOTATION AND TOTAL) AND AREAL EXTENT
        CALL VOLUME(NMAX, TIME, NUMNP, NUMEL, X, Y, KX, HTICE, BDROCK,
     &             DEPB, ADOT, RHOI, RHOW, VOL, AREA, AMASS,
     &             GRIDAREA,.FALSE.,sealow)
C
C        NTSTEP=NTSTEP+1
C        TTIME(NTSTEP)=TIME
C        VVOL(NTSTEP)=VOL*1.E-15
C        WWWMIN(NTSTEP)=WMIN
C        WWWMIN(NTSTEP)=WAVG
C        AAREA(NTSTEP)=AREA*1.E-12
C        TTBOT(NTSTEP)=TBAVG
C        TTAVG(NTSTEP)=TJUNK
C        TTNSL(NTSTEP)=AMASS(9)
C        TWATER(NTSTEP)=TOTALW*1E-9
C        PWATER(NTSTEP)=TOTALP
C
C
        LL=0
        LF=0
        IF(IPLOT.GT.0 .AND. IPLOT.LT.10) THEN
           CALL GRSTRT(600,1)
           CALL WINDOW(0.,100.,0.,100.)
        ENDIF
        IDONE=0
C
C
C
C ..... ****************************************************************
C ..... LOOP ON NUMBER OF TIME STEPS *******************************
C ..... ****************************************************************
        DO 450 L=1,NDT
          IF(DT.GT.0.0) THEN
            TIME=TIME+DT
          ELSE
            TIME=TIME+1
          ENDIF
          if(.not.ITOGG) then     ! HERE IS WHERE TFRACT GETS SET
            if(.false.) then
              TFRACT=coef(time,0.e0,1.e5,0.8e0,1.e0)
            elseif(.false.) then
              TFRACT=amass(9)
            elseif(.true.) then
              VOLMAX=1.7d6        ! npole cap volume km^3
              VOLMAX=5d6 * 1.6    ! current caps volume km^3

!             VOLMAX=VOLMAX*2     ! times 2
              VOLMAX=VOLMAX*5     ! times 5
!             VOLMAX=VOLMAX*10    ! times 10
!             VOLMAX=VOLMAX*15    ! times 15
!             VOLMAX=VOLMAX*20    ! times 20
! for KAT SCANLON'S RUNS, THIS (40) is 1X
!             VOLMAX = 0.125 * VOLMAX  ! 0.125x
!             VOLMAX = 0.25 * VOLMAX   ! 0.25x
!             VOLMAX = 0.5 * VOLMAX    ! 0.5x
!             VOLMAX = 1.0 * VOLMAX    ! 1.0x
              VOLMAX = 2.0 * VOLMAX    ! 2.0x
!             VOLMAX = 5.0 * VOLMAX    ! 5.0x
! for KAT SCANLON'S RUNS, THIS (40) is 1X

! special setting, VOLMAX in GEL and then converted to m^3
              VOLMAX = 34 ! current inventory GEL
              VOLMAX = VOLMAX * 16 ! what was 2X
c  1X case-------------------------------------------
              VOLMAX = 1*34  ! current inventory GEL
c  5X case-------------------------------------------
c             VOLMAX = 5*34  ! current inventory GEL
c 16X case-------------------------------------------
c             VOLMAX = 16*34 ! current inventory GEL
c end cases------------------------------------------
              VOLMAX = VOLMAX/1000 ! GEL in km
              VOLMAX = VOLMAX * 144.8d6 ! km^3
! special setting, VOLMAX in GEL and then converted to m^3

              VOLMAX=VOLMAX*1d9   ! m^3
              if(.true.) then ! use this when setting flist
                vratio=(vol/VOLMAX)
                TFRACT=(1-min(1.d0,vratio))
                write(55,*) time,tfract
              else             ! and this when reading flist
                read(54,*) tjunk,tfract
                write(55,*) tjunk,tfract
              endif
c             pause
            endif
          endif
          IF(SLTOGG.EQ.1 .and. NSEAL.GT.0) then
            SEALEV=FINDSEAL(TIME)
            DO I=1,NUMNP
              if(abs(htice(i)-SEALEV).lt.1e-5) then
                htice(i)=max(sealev,depb(i))
              endif
            ENDDO
          ELSEIF(SLTOGG.EQ.2) then
            SEALEV=SEALOW-SEABASE
            WRITE(77,*) TIME,SEALEV
            DO I=1,NUMNP
              if(abs(htice(i)-SEALEV).lt.1e-5) then
                htice(i)=max(sealev,depb(i))
              endif
            ENDDO
          endif
          IPAGE=1
          WRITE(LIST(IPAGE),*) 'TIME=',TIME
          WRITE(LIST(IPAGE+1),*) 'FLIST TFRACT=',REAL(TFRACT),REAL(VOL),
     &                            real(v ratio)
          IPAGE=IPAGE+1
C
C ..... FORM STIFFNESS,CAPACITANCE AND LOAD
C ..... DUMP DATA FOR FORMC
C        REWIND 73
C        WRITE(73) NUMNP,NUMEL,NUMGBC
C        WRITE(73),ETA,XI,W,LM
C        WRITE(73) (X(I),I=1,NUMNP)
C        WRITE(73) (Y(I),I=1,NUMNP)
C        WRITE(73) (T(I),I=1,NUMNP)
C        WRITE(73) (KODE(I),I=1,NUMNP)
C        WRITE(73) (KZ(I),I=1,NUMNP)
C        WRITE(73) (NTYPE(I),I=1,NUMEL)
C        WRITE(73) (CONST(I),I=1,NUMEL)
C        WRITE(73) (AADOT(I),I=1,NUMEL)
C        WRITE(73) ((KX(I,J),J=1,4),I=1,NUMEL)
C        WRITE(73) ((IBFLUX(I,J),J=1,2),I=1,NUMGBC)
C        WRITE(73) (BFLUX(I),I=1,NUMGBC)
C        PAUSE
C .....................................
C .... SKIP THE CONTINUITY SOLVER CALCULATION
C .... (DO IT, BUT SEE CTOGG BELOW WHERE SOLUTION IS ACCEPTED...)
        IF(CTOGG .OR. .FALSE.) THEN
C ......BEGIN EXPERIMENTAL...............
C ..... SET SURFACE, TTT, TO FLOTATION HEIGHT WHEREVER ZERO ...
          IF(.TRUE.) THEN
            DO JK=1,NUMNP
              IF(DEPB(JK).LT.SEALEV .AND. KODE(JK).EQ.0) THEN
                FLOT=(1.D0-RATDEN)*(DEPB(JK)-SEALEV)
                IF(ADOT(JK).GT.0) THEN
                  FLOT=(FLOT-ADOT(JK)*DT*1)*ADVANCE
                ELSE
                  FLOT=ADVANCE*FLOT
                ENDIF
                IF(TTT(JK).LT.FLOT) THEN
                  TTT(JK)=MIN(100.D0,FLOT)
C                  TTT(JK)=MAX(0.D0,TTT(JK))
                  TTT(JK)=MAX(SEALEV,TTT(JK))
                ENDIF
              ENDIF
            ENDDO
          ENDIF
C ......END EXPERIMENTAL...............
          WRITE(7,123) ' TIME BEFORE FORMC ',ETIME(TB),DTIME(TB)
          CALL FORMC(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL, 
     &                 ETA, XI, W, CONST,
     &                 NUMGBC, IBFLUX, BFLUX, NZ, KZ, LM, AADOT, 
     &                 ADOT,D, B, A, CAP,KA,BOLD,ADIAG)
          WRITE(7,123) ' TIME AFTER FORMC ',ETIME(TB),DTIME(TB)
C
          IF(.false.) THEN
C ......... BEGIN FORWARD DIFFERENCE SCHEME ..................
c ......... TNEW(N+1)=TTT(N)+DT*(D**-1)*(B-A*TTT(N))
            IF(DT.GT.0) THEN
c ..........  first the d-inverse, also mult by dt
              DO I=1,NUMNP
                IF (D(I).NE.0.) THEN
                  D(I)=DT/D(I)
                endif
              enddo
c ..........  form A*TTT store in FLUX (MATMULT RETURNS NEGATIVE)
              CALL MATMULT(NMAX,NZ,NZ1,NUMNP,A,KA,KZ,ADIAG,TTT,
     &                     ZERO,FLUX)
c ..........  dot with d-invers and change sign, keep in TNEW ...
              DO I=1,NUMNP
                TNEW(I)=d(i)*(b(i)-FLUX(i))
              enddo
c ..........  and add to the old vector, if NOT fixed...
              DO I=1,NUMNP
                IF(kode(i).eq.0) then
                  TNEW(I)=TTT(i)+TNEW(I)
                else
                  TNEW(i)=TTT(i)
                endif
              enddo
              DO I=1,NUMNP
                if(real(tnew(i)-ttt(i)).gt.0.and..false.) then
                  write(73,*) i
                  write(73,*) 'tnew   :',real(tnew(i))
                  write(73,*) 'change :',real(tnew(i)-ttt(i))
                  write(73,*) 'b      :',real(b(i))
                  write(73,*) 'k*t    :',real(flux(i))
                  write(73,*) 'b-k*t  :',real(b(i)-flux(i))
                  write(73,*) 'd      :',real(d(i))
                  write(73,*) 'change :',real(d(i)*(b(i)-FLUX(i)))
                endif
              enddo
            ELSE
              PRINT *,'ASKING FOR STEADY STATE, not ready...'
              STOP
            ENDIF
C ......... END FORWARD DIFFERENCE SCHEME ..................
          ELSE
C ......... BEGIN BACKWARD DIFFERENCE SCHEME ..................
C ......... CALCULATE EFFECTIVE LOAD  AND STIFFNESS MATRIX
            IF(DT.GT.0.0) THEN
              DT2=1.0D0/DT
              DO I=1,NUMNP
                IF (KODE(I).EQ.0) THEN
c ... LUMPED CAPACITANCE ...
                  if(.true.) then
                    IF (D(I).NE.0.) THEN
                      D(I)=DT2*D(I)
c      print *,i,real(a(i,1)),real(d(i)),
c     &          real(b(i)),real(d(i)*TTT(i))
                      A(I,1)=A(I,1)+D(I)
                      B(I)=B(I)+D(I)*TTT(I)
C                     BOLD(I)=BOLD(I)+D(I)*TTT(I)
                    ENDIF
                  else
c ... FULL CAPACITANCE ...
                    vect=0.d0
                    do j=1,kz(i)
                      cap(i,j)=cap(i,j)*dt2
                      A(I,j)=A(I,j)+cap(I,j)
                      vect=vect+cap(i,j)*ttt(ka(i,j))
                    enddo
                    b(i)=b(i)+vect
                  endif
                ENDIF
              ENDDO
            ENDIF





C
C ... FIXED BOUNDARY CONDITIONS BY PENALTY METHOD
C
C$DOACROSS LOCAL(N)
            DO N=1,NUMNP
              IF(KODE(N).EQ.1) THEN
                A(N,1)=BIG
                B(N)=TTT(N)*BIG
              ENDIF
              KA(N,NZ+1)=KZ(N)
            ENDDO
C
C ......... REMOVE HERE, THIS WRITES OUT MATRIX
            if(.false.) then
              REWIND 73
              WRITE(73,*) NUMNP
              DO I=1,NUMNP
                WRITE(73,*) KA(I,10)
                DO J=1,KA(I,10)
                  WRITE(73,*) I,KA(I,J),A(I,J)
                ENDDO
              ENDDO
              WRITE(73,*) (TTT(I),I=1,NUMNP)
              WRITE(73,*) (B(I),I=1,NUMNP)
C             PAUSE
            endif
C ......... TO HERE ******
C$DOACROSS LOCAL(JK)
            DO JK=1,NUMNP
              TNEW(JK)=TTT(JK)
            ENDDO
C
C ......... CONJUGATE-GRADIENT ITERATIVE SOLVER
C           SEC1=WHEN()
            WRITE(7,123) ' TIME BEFORE CONJUG ',ETIME(TB),DTIME(TB)
            CALL CONJUG(NMAX,NZ,NUMNP,1.D-10,A,KA,B,TNEW)
c              WRITE(73,*) (tnew(I),I=1,NUMNP)
            WRITE(7,123) ' TIME AFTER CONJUG ',ETIME(TB),DTIME(TB)
C           CALL CONJUG(NMAX,NZ,NUMNP,0.D-6,A,KA,B,TNEW)
C           SEC2=WHEN()
C           PRINT *,SEC2-SEC1
C ......... GAUSS-SEIDEL ITERATIVE SOLVER
C           SEC1=WHEN()
C           CALL GAUSEID(NMAX,NZ,NUMNP,A,KA,B,TNEW)
C           SEC2=WHEN()
C           PRINT *,SEC2-SEC1
C
C ......... END BACKWARD DIFFERENCE SCHEME ..................
          ENDIF
C ....... WITH CTOGG OFF, WILL KEEP ORIGINAL SOLUTION, BUT DO ALL REST
          IF(CTOGG) THEN
C$DOACROSS LOCAL(JK)
            DO JK=1,NUMNP
              if(tnew(jk)-ttt(jk).ge.0.0 .and. HIRES) then
c                print *,'----',jk,real(ttt(jk)),
c    &                  real(tnew(jk)),
c    &                  real(tnew(jk))-real(ttt(jk)),
c    &                  real(adot(jk))
                IF(DEPB(JK).LT.SEALEV) THEN
                  FLOT=(1.D0-RATDEN)*(DEPB(JK)-SEALEV)
                  if(tnew(jk).lt.sealev) then
                    tnew(jk)=SEALEV+tnew(jk)-ttt(jk)
                  endif
                ELSE
                  if(tnew(jk).lt.DEPB(JK)) then
                    tnew(jk)=DEPB(JK)+tnew(jk)-ttt(jk)
                  endif
                ENDIF
              endif
              QQQ(JK)=TNEW(JK)
            ENDDO
          ENDIF
        ENDIF
C ..... TO HERE SKIPPING CONTINUITY SOLVER
C
          LL=LL+1
          LF=LF+1
          HMAX=-1.D30
          DIFF=0.D0
          NDIFF=0
C ....... MEASURES DIFFERENCE
C ....... BETWEEN PRESENT SOLUTION AND PRESENT SURFACE.
C ....... ALSO SETS ITKODE, THE BC-KODE FOR THE WATER THICKNESS CALCULATION
          DO JK=1,NUMNP
            ITKODE(JK)=KODE(JK)
            ITKODE(JK)=0
            IF(DEPB(JK).LT.SEALEV) THEN
              FLOT=(1.D0-RATDEN)*(DEPB(JK)-SEALEV)
              THIK=QQQ(JK)-FLOT
              SURF=SEALEV
C             SURF=PSURF(JK)
              IF(THIK.LT.0. .AND. .NOT.HIRES) THEN
                QQQ(JK)=SURF
              ENDIF
            ELSE
              THIK=QQQ(JK)-DEPB(JK)
              SURF=DEPB(JK)
            ENDIF
            IF(THIK.LT.0.) THEN
c ..........  let it get negative ...
              IF(.NOT.HIRES) QQQ(JK)=SURF
              IF(IMELT.EQ.1) BMELT(JK)=-0.1D0
              IF(IMELT.EQ.1) BMELT(JK)=0D0
C              WTHICK(JK)=-1D-6
C              WTHICK(JK)=0D0
C              ITKODE(JK)=1
            ENDIF
            IF(KODE(JK).NE.1) THEN
              DIFF=DIFF+(QQQ(JK)-PSURF(JK))**2
              NDIFF=NDIFF+1
            ELSE
              ITKODE(JK)=KODE(JK)
              ITKODE(JK)=0
            ENDIF
            IF(QQQ(JK).GT.HMAX) THEN
              HMAX=QQQ(JK)
              NNMAX=JK
            ENDIF 
C ... UPDATE BASAL WATER THICKNESS (ONLY FOR QUICK AND DIRTY WATER MOVER)
C            IF(BMELT(JK).NE.0.0) THEN
C              WTHICK(JK)=MAX(0.0,WTHICK(JK)+BMELT(JK)*DT)
C       PRINT *,JK,WTHICK(JK),BMELT(JK)
C            ENDIF
          ENDDO
C .....   FOLLOWING ACCEPT SOLUTION AND
C .....   CHECKS TO MAKE SURE THE SURFACE IS NOT BELOW   !
C .....   THE BED OR THE FLOTATION LINE                            !
          DO I=1,NUMNP                                             !
C ......... ACCEPT SOLUTION ...
            HTICE(I)=QQQ(I)
            IF(DEPB(I).LT.SEALEV) THEN                                 !
              FLOT=(1.D0-RATDEN)*(DEPB(I)-SEALEV)                          !
              THIK=HTICE(I)-FLOT                                   !
              SURF=SEALEV                                              !
            ELSE                                                   !
              THIK=HTICE(I)-DEPB(I)                                !
              SURF=DEPB(I)                                         !
            ENDIF                                                  !
            IF(THIK.LT.0.) THEN                                    !
              HTICE(I)=SURF                                        !
            ENDIF                                                  !
c      if(htice(i).ne.qqq(i)) print *,i,htice(i),qqq(i)
          ENDDO                                                    !
C ....... QUICK AND DIRTY WATER MOVER, AVERAGE EACH ELEMENT
          IF(.FALSE.) THEN
            DO I=2,NUMCOL-1
              DO J=2,NUMLEV-1
                IP=(J-1)*NUMCOL+I
                WTELEM(IP)=0.25d0*(WTHICK(IP+1)+WTHICK(IP-1)+
     &                           WTHICK(IP+NUMCOL)+WTHICK(IP-NUMCOL))
              ENDDO
            ENDDO
            DO I=2,NUMCOL-1
              DO J=2,NUMLEV-1
                IP=(J-1)*NUMCOL+I
                WTHICK(IP)=WTELEM(IP)
              ENDDO
            ENDDO
          ELSEIF(.TRUE.) THEN
            IF(WTOGG.gt.0) THEN
C ......... FEM WATER MOVER, CONSERVATION EQUATION ...
              WRITE(7,123) ' TIME BEFORE WMOVER ',ETIME(TB),DTIME(TB)
              IF(.FALSE.) THEN
                REWIND(67)
                WRITE(67) NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL, DT,
     &                   DTLOCAL,
     &                   ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &                   NZ, KZ, LM, TNEW,
     &                   BMELT, D, B, A, KA, ALPHAC, TOTALW, TOTALP
              ENDIF
              IF(WTOGG.eq.2) THEN
                CALL WMOVER(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL, DT,
     &                   DTLOCAL,
     &                   ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &                   NZ, KZ, LM, TNEW,
     &                   BMELT, D, B, A, KA, ALPHAC, TOTALW, TOTALP,
     &                   IPLOT)
              ELSEIF(WTOGG.eq.1) THEN
                CALL WMOVERSIMP(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL, DT,
     &                 ETA, XI, W, WTHICK, HTICE, DEPB, LM,
     &                 BMELT, ALPHAC, TOTALW, TOTALP,IPLOT)
              ELSEIF(WTOGG.eq.3) THEN
                CALL WMOVERJ(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL, DT,
     &                   DTLOCAL,
     &                   ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &                   NZ, KZ, LM, TNEW,
     &                   BMELT, D, B, A, KA, ALPHAC, TOTALW, TOTALP,
     &                   IPLOT,WVELX,WVELY)
              ELSE
                PRINT *,'WTOGG:',WTOGG,' NOT 1,2, OR 3'
                STOP
              ENDIF
c ... fix nodes with no ice, below flotation line to have water
              if(.false.) then
                DO I=1,NUMNP
                  FLOT=(1.D0-RATDEN)*(BDROCK(I)-SEALEV)
                  IF(HTICE(I).LT.FLOT) THEN
                    wthick(I)=1.d0
                  ENDIF
                ENDDO
              endif
c ... fix nodes with no ice, below flotation line to have water

              WRITE(7,123) ' TIME AFTER WMOVER ',ETIME(TB),DTIME(TB)
            ENDIF
          ENDIF
C ............................................................
C ....... SECTION TO CALCULATE FLUX A*X=B, THEN B-BOLD=FLUX
          CALL MATMULT(NMAX,NZ,NZ1,NUMNP,A,KA,KZ,ADIAG,QQQ,BOLD,FLUX)
          CALL FLUXUP(NMAX,NUMNP,NUMGBC,IBFLUX,BFLUX,FLUX,
     &         X,Y,KODE,TIME)
C
C ....... OUTPUT TIME STEP INFO TO SCREEN
          IF(IOTOGG) THEN
            WRITE(LIST(IPAGE+1),*) 'MAXSURF=',HMAX,
     &              'NODE',NNMAX,' DIF=',
     &              REAL(SQRT(DIFF/DBLE(NDIFF)))
            IPAGE=IPAGE+1
          ENDIF
C
C ....... DERIVE SLOPE IN CASE NEEDED BY MASS BALANCE PARAMETERIZATION
          WRITE(7,123) ' TIME BEFORE NODESL ',ETIME(TB),DTIME(TB)
          CALL NODESL(NMAX, NUMNP, NUMEL, KX, SLOPE, SLOPN)
          WRITE(7,123) ' TIME AFTER NODESL ',ETIME(TB),DTIME(TB)
C
C ....... MODIFY TEMPERATURE FIELD...
          CALL DVELO(NMAX,NUMEL,KX,X,Y,HTICE,DEPB,CONST,
     &             VEL,.FALSE.,DTMIN)
          WRITE(7,123) ' TIME BEFORE TEMPER ',ETIME(TB),DTIME(TB)
          CALL RESETTEMP(NUMNP,AMASS(9),HTICE,TEMP,AMASS,TIME)
          CALL TEMPER(NMAX, NUMNP, NUMEL, NTYPE, KX, KODE, ADOT, FRACT,
     &            DEPB, HTICE,PSURF,FLOWA, SLDGB, TEMP, ITYPE, TBED,
     &            X, Y, AMASS, SLOPN, TIME, DT, IPLOT,TBAVG,
     &            AFUDGE,GRIDAREA,IMELT,BMELT,WTHICK,VOL,AREA,
     &            GEOFLUX,VEL,IDT,CONTRIB)
          WRITE(7,123) ' TIME AFTER TEMPER ',ETIME(TB),DTIME(TB)
C ....... UPDATE ACCUMULTION RATES FOR NEW COBFIGURATION (NEW POSITION)
          !print *,'here before ACCSET';pause
          CALL ACCSET(NUMNP,AMASS,ACOMSAVE,TIME,
     &                  THICK,IDT,HTICE,BDROCK,PSURF,SLOPN,X,Y,TEMP,
     &                  ADOT,ADOTB,ADP,ADM,TJUNK,TFRACT,ACC,ABLAT)
        CALL RESETTEMP(NUMNP,AMASS(9),HTICE,TEMP,AMASS,TIME)

C---------------------------------------------------------------------
      IF(.FALSE.) THEN
C**********EISMINT ALA HUYBRECHTS ********************
          WRM=WARMING(AMASS(9)+TNSLBASE)
C*****************************************************
C$DOACROSS LOCAL(I,ADOTI,AJUNK), SHARED(TIME)
          DO I=1,NUMNP
C            IF(THICK(I).GT.0) THEN
C              ACOM=ACOMSAVE
C            ELSE
C              ACOM=ACOMSAVE/3.
C            ENDIF
            ADOTI=AFUNCT(TIME, IDT(I), AMASS,
     &                     HTICE(I), BDROCK(I),PSURF(I),
     &                     SLOPN(1,I),X(I),Y(I),TEMP(I))
            IF(IDT(I).GT.0) THEN
              ADOT(I)=ADOTI
            ELSE
              IF(ITOGG) THEN
                AJUNK=ACCUM(IDT(I),X(I)*.001D0,Y(I)*.001D0,QQQ(I),
     &                      SLOPN(1,I),
     &                      PSURF(I),AMASS(9),TEMP(I))
C               ADOT(I)=ADOTB(I)-ABL*.01
C*****************************************************
C**********EISMINT ALA HUYBRECHTS ********************
                ADP(I)=WRM*ADOTB(I)
                ADM(I)=-ABL*0.01D0
C STRAIGHT
                ADOT(I)=WRM*ADOTB(I)-ABL*.01D0
C SLOW DOWN ACCUM RATE CHANGE...-------------------------
C                ADOT(I)=0.5D0*(ADOT(I)+WRM*ADOTB(I)-ABL*.01)
C--------------------------------------------------------
C*****************************************************
              ELSE
                ADOT(I)=ADOTB(I)
              ENDIF
            ENDIF
            CALL STEEPSET(I,HTICE(I),BDROCK(I),SLOPN(1,I),ADOT(I))
          ENDDO
      ENDIF
C---------------------------------------------------------------------

          CALL CALVING(NMAX, NUMEL, NTYPE, KX, HTICE, DEPB, ADOT, 
     &           AADOT,CFACTOR,ADC,CALV,PCALV,AMASS(9))
          CALL ELPROP(NMAX, NUMEL, NTYPE, KX, ADOT, 
     &             AADOT, FRACT,
     &             AFRACT, BDROCK, ABDRCK, PSURF, PPSURF, FLOWA,ACON,
     &             ICON, AFLOWA, SLDGB, ASLDGB,CALV,PCALV)
C
C ....... OBTAIN NEW LINEARIZATION CONSTANT FROM LATEST SOLUTION 
          WRITE(7,123) ' TIME BEFORE NCONST ',ETIME(TB),DTIME(TB)
          CALL NCONST(NMAX, X, Y, KX, NTYPE, NUMEL,
     &             AFRACT, ASLDGB, LM, AFLOWA, BDROCK, DEPB, UNDEPB,
     &             PG, QQQ, CNEW, SLOPE, RHOI, WINDIR,
     &             WTHICK,ITYPE)
          WRITE(7,123) ' TIME AFTER NCONST ',ETIME(TB),DTIME(TB)
C
C ....... UPDATE ACCUMULTION RATES FOR NEW COBFIGURATION
C          CALL ACCSET(NUMNP,AMASS,ACOMSAVE,TIME,
C     &                  THICK,IDT,HTICE,BDROCK,PSURF,SLOPN,X,Y,TEMP,
C     &                  ADOT,ADOTB,ADP,ADM,TJUNK,TFRACT,ACC,ABLAT)
C          CALL RESETTEMP(NUMNP,AMASS(9),HTICE,TEMP,AMASS,TIME)

C-------------------------------------------
C     IF(.FALSE.) THEN
C**********EISMINT ALA HUYBRECHTS ********************
C          WRM=WARMING(AMASS(9)+TNSLBASE)
C*****************************************************
C          DO I=1,NUMNP
C            ADOTI=AFUNCT(TIME, IDT(I), AMASS,
C     &                     HTICE(I), BDROCK(I),PSURF(I),
C     &                     SLOPN(1,I),X(I),Y(I),TEMP(I))
C            IF(IDT(I).GT.0) THEN
C              ADOT(I)=ADOTI
C            ELSE
C              IF(ITOGG) THEN
C                AJUNK=ACCUM(IDT(I),X(I)*.001D0,Y(I)*.001D0,QQQ(I),
C     &                      SLOPN(1,I),
C     &                      PSURF(I),AMASS(9),TEMP(I))
C                ADOT(I)=ADOTB(I)-ABL*.01
C*****************************************************
C**********EISMINT ALA HUYBRECHTS ********************
C                ADP(I)=WRM*ADOTB(I)
C                ADM(I)=-ABL*0.01
C                ADOT(I)=WRM*ADOTB(I)-ABL*.01
C*****************************************************
C              ELSE
C                ADOT(I)=ADOTB(I)
C              ENDIF
C            ENDIF
C          ENDDO
C     ENDIF
C-------------------------------------------

C
C ....... CALCULATE VOLUMES (FLOTATION AND TOTAL) AND AREA
          WRITE(7,123) ' TIME BEFORE VOLUME ',ETIME(TB),DTIME(TB)
          CALL VOLUME(NMAX, TIME, NUMNP, NUMEL, X, Y, KX, HTICE, BDROCK,
     &             DEPB, ADOT, RHOI, RHOW, VOL, AREA, AMASS,
     &             GRIDAREA,.TRUE.,sealow)
          WRITE(7,123) ' TIME AFTER VOLUME ',ETIME(TB),DTIME(TB)
C
          IF(.TRUE.) THEN
            CALL DUMPBC(TIME,NUMNP,X,Y,QQQ,FLUX,KODE,
     &                  NUMGBC,IBFLUX,BFLUX)
          ENDIF
          NTSTEP=NTSTEP+1
          TTIME(NTSTEP)=TIME
          VVOL(NTSTEP)=VOL*1.D-15
          AAREA(NTSTEP)=AREA*1.D-12
          TTBOT(NTSTEP)=TBAVG
          TTAVG(NTSTEP)=TJUNK
c         WRITE(17,*) TJUNK,'       POLE, SURF=0 TEMP'
          WRITE(17,*) -999,'       junk'
          WRITE(17,*) TBAVG,'       AVERAG BASAL TEMP'
          TTNSL(NTSTEP)=AMASS(9)
          TTSEAL(NTSTEP)=SEALEV
          TWATER(NTSTEP)=TOTALW*1D-9
          PWATER(NTSTEP)=TOTALP
          if(TWATER(NTSTEP).gt.0) then
            WRITE(17,*) TWATER(NTSTEP),'       TOTAL   WATER NONZERO'
          else
            WRITE(17,*) TWATER(NTSTEP),'       TOTAL   WATER'
          endif
          WRITE(17,*) PWATER(NTSTEP),'       PERCENT WATER'
          DO I=1,NUMNP
            TTT(I)=QQQ(I)
          ENDDO
C
C ....... LOAD NEW LINEARIZATION CONSTANT INTO OLD
          IF(DT.EQ.0.) CONV=5.D0
          DO I=1,NUMEL
C            CONST(I)=CNEW(I)
            CONST(I) = (CONV*CONST(I) + CNEW(I))/(CONV+1.D0)
          ENDDO
C
C$DOACROSS LOCAL(I)
c          DO N=1,NUMNP
c            HTICE(N)=QQQ(N)
c          ENDDO
C......EXPERIMENTAL SECTION FOR BED DEPRESSION ......................!
          WRITE(7,123) ' TIME BEFORE PLATE ',ETIME(TB),DTIME(TB)
          IF(BTOGG.EQ.1) THEN
            DO I=1,NUMNP                                             !
              THICKL(I)=HTICE(I)-DEPB(I)                             !
              IF(DEPB(I).LT.SEALEV) THEN
                FLOT=(1.D0-RATDEN)*(DEPB(I)-SEALEV)                          !
                IF(HTICE(I).LT.FLOT) THEN                              !
                  THICKL(I)=0.D0                                       !
                  THICKL(I)=FLOT-DEPB(I)                                       !
                ENDIF
              ENDIF
              IF(DEPB(I).LE.-9999.) THICKL(I)=0.D0                                       !
              IF(THICKL(I).LT.0.) THICKL(I)=0.D0                                       !
            ENDDO
            CALL SPLATE(NTSTEP,NUMNP,THICKL,UNDEPB,DEPB,RHOI,RHOR,RHOW,
     &                  DT,TIME,WRATE,WWW,WWWORIG)
            WMIN=1D30
            WMAX=-1D30
            II=1
            DO I=1,NUMNP
              WDIFF(I)=WWW(II)-WWWORIG(I)
              WMIN=MIN(WMIN,WDIFF(I))
              WMAX=MAX(WMAX,WDIFF(I))
              WAVG=WAVG+WWW(II)
              II=II+3
            ENDDO
            WMIN=WMIN
            WAVG=WAVG/NUMNP
            IF(IOTOGG) THEN
              WRITE(LIST(IPAGE+1),*) 'DEPRESSION DIFFERENCE=',
     &                                REAL(WMIN),REAL(WMAX),REAL(WAVG)
              IPAGE=IPAGE+1
            ENDIF
            DO I=1,NUMNP                                             !
C              HTICE(I)=HTICE(I)+WWW((I-1)*3+1)                       !
              DEPB(I)=UNDEPB(I)+WWW((I-1)*3+1)                       !
              IF(HTICE(I).LT.DEPB(I)) HTICE(I)=DEPB(I)               !
            ENDDO                                                    !
            IF(LL.GE.INTER) THEN                                     !
              CALL POUTSTF(3*NUMNP,WWW,WRATE(1,1),NUMNP,THICKL,      !
     &                     NUMCOL,NUMLEV,TIME,WDIFF)                 !
            ENDIF
C .......   FOLLOWING CHECKS TO MAKE SURE THE SURFACE IS NOT BELOW   !
C .......   THE BED OR THE FLOTATION LINE                            !
            DO I=1,NUMNP                                             !
              IF(DEPB(I).LT.SEALEV) THEN                                 !
                FLOT=(1.D0-RATDEN)*(DEPB(I)-SEALEV)                          !
                THIK=HTICE(I)-FLOT                                   !
                SURF=SEALEV                                              !
              ELSE                                                   !
                THIK=HTICE(I)-DEPB(I)                                !
                SURF=DEPB(I)                                         !
              ENDIF                                                  !
              IF(THIK.LT.0.) THEN                                    !
                HTICE(I)=SURF                                        !
              ENDIF                                                  !
            ENDDO                                                    !
          ELSEIF(BTOGG.EQ.2 .OR. BTOGG.EQ.3) THEN                                             !
            DO I=1,NUMNP                                             !
              THICKL(I)=HTICE(I)-DEPB(I)                             !
              IF(DEPB(I).LT.SEALEV) THEN
                FLOT=(1.D0-RATDEN)*(DEPB(I)-SEALEV)                          !
                IF(HTICE(I).LT.FLOT) THEN                              !
                 THICKL(I)=0.D0                                       !
                 THICKL(I)=FLOT-DEPB(I)                                         !
                ENDIF
              ENDIF
              IF(THICKL(I).LT.0.) THICKL(I)=0.D0                                       !
              IF(DEPB(I).LE.-9999.) THICKL(I)=0.D0                                       !
            ENDDO
            IF(BTOGG.EQ.3) THEN
              CALL VEPLATE(NTSTEP,NUMNP,NUMEL,X,Y,KX,THICKL,KODE,
     &                   DT,WWW,WRATE,TIME,WWWORIG,
     &                   WMIN,WMINE,WMINV,FNETE,FNETV)
            ELSEIF(BTOGG.EQ.2) THEN
              CALL OPLATE(NTSTEP,NMAX,N3,NUMNP,NUMEL,X,Y,KX,THICKL,
     &                 DT,WWW,WRATE,WMIN,TIME,WWWORIG)
            ELSE
              PRINT *,' PROBLEMS WITH BTOGG ', BTOGG
              PAUSE
            ENDIF
            WMIN=1D30
            WMAX=-WMIN
            WAVG=0.D0
            II=1
            DO I=1,NUMNP
              WDIFF(I)=WWW(II)-WWWORIG(I)
              WMIN=MIN(WMIN,WDIFF(I))
              WMAX=MAX(WMAX,WDIFF(I))
              WAVG=WAVG+WWW(II)
              II=II+3
            ENDDO
            WAVG=WAVG/NUMNP
            IF(IOTOGG) THEN
              WRITE(LIST(IPAGE+1),*) 'DEPRESSION DIFFERENCE=',WMIN,WAVG
              IPAGE=IPAGE+1
            ENDIF
            DO I=1,NUMNP                                             !
C              HTICE(I)=HTICE(I)+WWW((I-1)*3+1)                       !
              DEPB(I)=UNDEPB(I)+WWW((I-1)*3+1)                       !
              IF(HTICE(I).LT.DEPB(I)) HTICE(I)=DEPB(I)               !
            ENDDO                                                    !
            IF(LL.GE.INTER) THEN                                     !
              CALL POUTSTF(NUMNP*3,WWW,WRATE(1,1),NUMNP,THICKL,      !
     &                     NUMCOL,NUMLEV,TIME,WDIFF)                 !
            ENDIF                                                    !
C .......   FOLLOWING CHECKS TO MAKE SURE THE SURFACE IS NOT BELOW   !
C .......   THE BED OR THE FLOTATION LINE                            !
            DO I=1,NUMNP                                             !
              IF(DEPB(I).LT.SEALEV) THEN                                 !
                FLOT=(1.D0-RATDEN)*(DEPB(I)-SEALEV)                          !
                THIK=HTICE(I)-FLOT                                   !
                SURF=SEALEV                                              !
              ELSE                                                   !
                THIK=HTICE(I)-DEPB(I)                                !
                SURF=DEPB(I)                                         !
              ENDIF                                                  !
              IF(THIK.LT.0.) THEN                                    !
                HTICE(I)=SURF                                        !
              ENDIF                                                  !
            ENDDO                                                    !
          ENDIF                                                      !
          WRITE(7,123) ' TIME AFTER  PLATE ',ETIME(TB),DTIME(TB)
C....................................................................!
          WWWMIN(NTSTEP)=WMIN
          WWWMIN(NTSTEP)=WAVG
          WRITE(17,*) WWWMIN(NTSTEP),'       BED DEPRESSION'
C ....... EISMINT PRINT OUT ........
          WRITE(7,123) ' TIME BEFORE PRINTOUT ',ETIME(TB),DTIME(TB)
          CALL DVELO(NMAX,NUMEL,KX,X,Y,HTICE,DEPB,CONST,
     &             VEL,.FALSE.,DTMIN)
          IF(.FALSE.) THEN
            IF(DTMIN.LT.1.) THEN
              DT=DTMIN
            ELSEIF(DTMIN.LT.10.) THEN
              DT=NINT(DTMIN)
            ELSEIF(DTMIN.LT.100.) THEN
              DT=10.D0*NINT(DTMIN/10.D0)
            ELSEIF(DTMIN.LT.1000.) THEN
              DT=100.D0*NINT(DTMIN/100.D0)
            ELSEIF(DTMIN.LT.10000.) THEN
              DT=1000.D0*NINT(DTMIN/1000.D0)
            ELSE
              DT=1000.D0
            ENDIF
          ENDIF
          CALL EISMINT(TIME,X,Y,HTICE,ADOT,ADOTB,DEPB,THICK,TBED,VOL,
     &                 AREA,GRIDAREA,VEL,SLOPN,ADP,ADM,ADC,BMELT,
     &                 NUMNP,NUMEL,KX)
C ....... PRINT OUT SOLUTION FOR APPROPRIATE TIME STEPS
          IF(LL.GE.INTER) THEN
            WRITE(LIST(IPAGE+1),*) '. . . . . . . . . . . . . . . . .'
            WRITE(LIST(IPAGE+2),*) '. DUMPING SOLUTION  ',real(TIME)
            WRITE(LIST(IPAGE+3),*) '. . . . . . . . . . . . . . . . .'
            IPAGE=IPAGE+3
            WRITE(SCRTCH,*) 'TIME=',NINT(TIME)
            WRITE(34) SCRTCH
            WRITE(34) (HTICE(I),I=1,NUMNP)
            WRITE(34) (ADOT(I),I=1,NUMNP)
            WRITE(34) (DEPB(I),I=1,NUMNP)
            WRITE(34) (CONST(I),I=1,NUMEL)
            WRITE(34) (ACON(I),I=1,NUMEL)
            WRITE(35) (FRACT(I),I=1,NUMNP)
            WRITE(35) (FLOWA(I),I=1,NUMNP)
            WRITE(35) (SLDGB(I),I=1,NUMNP)
            WRITE(35) (AFUDGE(I),I=1,NUMNP)
            WRITE(36) SCRTCH
            WRITE(36) (TBED(I),I=1,NUMNP)
            WRITE(36) (BMELT(I),I=1,NUMNP)
            WRITE(36) (WTHICK(I),I=1,NUMNP)
            WRITE(36) (CONTRIB(I),I=1,NUMNP)
            LL=0
          ENDIF
          WRITE(7,123) ' TIME AFTER PRINTOUT ',ETIME(TB),DTIME(TB)
          IF(LF.GE.IFIT(7)) THEN
c            DO N=1,NUMNP
c              HTICE(N)=QQQ(N)
c            ENDDO
            LF=0
            CALL MODIFY(NMAX,IFIT,NUMNP,KODE,HTICE,DEPB,PSURF,
     &                  ADOT,FRACT,FLOWA,SLDGB,BDROCK,AFUDGE,IDT,
     &                  NUMEL,KX)
            CALL UNLOAD(NUMNP,BDROCK,UNDEPB,PSURF,RHOI,RHOR)
          ENDIF
          IF(IPLOT.EQ.7) THEN
            CALL TPROF(NUMNP,HTICE,DEPB,.FALSE.,1,TIME,
     &                 BMELT,WTHICK)
          ELSEIF(IPLOT.EQ.9) THEN
            IPLOTLOC=6
            CALL PLOTSOL(NUMNP,X,Y,HTICE,DEPB,KODE,FRACT,
     &                   PSURF,WTHICK,TWATER,PWATER,AFUDGE,
     &                   NTSTEP,TTIME,VVOL,
     &                   AAREA,TTBOT,TTAVG,TTNSL,TTSEAL,
     &                   IPLOTLOC,
     &                   NUMEL,KX,CONST,VEL,WWWMIN,DTMIN)
            CALL FLUXES(2,NUMNP,NUMEL,X,Y,HTICE,DEPB,KX,
     &                  NUMGBC,IBFLUX,BFLUX)
          ELSEIF(IPLOT.NE.0) THEN
            CALL PLOTSOL(NUMNP,X,Y,HTICE,DEPB,KODE,FRACT,
     &                   PSURF,WTHICK,TWATER,PWATER,AFUDGE,
     &                   NTSTEP,TTIME,VVOL,
     &                   AAREA,TTBOT,TTAVG,TTNSL,TTSEAL,
     &                   IPLOT,
     &                   NUMEL,KX,CONST,VEL,WWWMIN,DTMIN)
            IF(IPLOT.EQ.2.OR.IPLOT.EQ.3.OR.IPLOT.EQ.4) THEN
              CALL DRAWOUT
            ENDIF
          ENDIF
          WRITE(7,123) ' >>>>>>TIME ITERATION END ',ETIME(TB),DTIME(TB)
c         IF(L.GT.1) THEN
c           DO J=1,NPAGE
c             PRINT *,J,'-----------------------'
c           ENDDO
c         print string
c         ENDIF
          IF(L.LT.NDT) THEN
            WRITE(*,'(A72)') (LIST(J),J=1,NPAGE)
          ELSE
            WRITE(*,'(A72)') (LIST(J),J=1,NPAGE-4)
          ENDIF
          DO J=1,NPAGE
            WRITE(LIST(J),*) j,'--------------'
          ENDDO
c         call sleep(1)
          IF(IFLUSH) CALL GFLUSH
          if(stopflag) then
            WRITE(LIST(IPAGE+1),*) '. . . . . . . . . . . . . . . . .'
            WRITE(LIST(IPAGE+2),*) '. DUMPING SOLUTION  ',real(TIME)
            WRITE(LIST(IPAGE+3),*) '. . . . . . . . . . . . . . . . .'
            IPAGE=IPAGE+3
            WRITE(SCRTCH,*) 'TIME=',NINT(TIME)
            WRITE(34) SCRTCH
            WRITE(34) (HTICE(I),I=1,NUMNP)
            WRITE(34) (ADOT(I),I=1,NUMNP)
            WRITE(34) (DEPB(I),I=1,NUMNP)
            WRITE(34) (CONST(I),I=1,NUMEL)
            WRITE(34) (ACON(I),I=1,NUMEL)
            WRITE(35) (FRACT(I),I=1,NUMNP)
            WRITE(35) (FLOWA(I),I=1,NUMNP)
            WRITE(35) (SLDGB(I),I=1,NUMNP)
            WRITE(35) (AFUDGE(I),I=1,NUMNP)
            WRITE(36) SCRTCH
            WRITE(36) (TBED(I),I=1,NUMNP)
            WRITE(36) (BMELT(I),I=1,NUMNP)
            WRITE(36) (WTHICK(I),I=1,NUMNP)
            WRITE(36) (CONTRIB(I),I=1,NUMNP)
            !goto 451
          endif
450     CONTINUE
451     CONTINUE
        stopflag=.false.
C
C ..... ****************************************************************
C ..... **** END OF TIME STEP SECTION **********************************
C ..... ****************************************************************
C
         IF(IPLOT.NE.0) CALL GRSTOP1
        HMAX=-1.D30
        DIFF=0.D0
        NDIFF=0
C ..... FOLLOWING CHECKS TO MAKE SURE THE SURFACE IS NOT BELOW 
C ..... THE BED OR THE FLOTATION LINE. ALSO MEASURES DIFFERENCE
C ..... BETWEEN PRESENT SOLUTION AND PRESENT SURFACE.
        DO JK=1,NUMNP
          HTICE(JK)=QQQ(JK)
          IF(DEPB(JK).LT.SEALEV) THEN
            FLOT=(1.D0-RATDEN)*(DEPB(JK)-SEALEV)
            THIK=HTICE(JK)-FLOT
            SURF=SEALEV
          ELSE
            THIK=HTICE(JK)-DEPB(JK)
            SURF=DEPB(JK)
          ENDIF
          IF(THIK.LT.0.) THEN
C***********THIS IS A TEST************************
C EISMINT EXPERIMENT, LET SURFACE STAY BELOW FLOTATION LINE
C BUT DONT LET IT GET NEGATIVE...
C            HTICE(JK)=SURF-10.
C***********THIS IS A TEST************************
            HTICE(JK)=SURF
            IF(IMELT.EQ.1) BMELT(JK)=-0.1D0
            IF(IMELT.EQ.1) BMELT(JK)=0D0
C            WTHICK(JK)=-1D-6
C            WTHICK(JK)=0D0
          ENDIF
          IF(KODE(JK).NE.1) THEN
            DIFF=DIFF+(HTICE(JK)-PSURF(JK))**2
            NDIFF=NDIFF+1
          ENDIF
          IF(HTICE(JK).GT.HMAX) THEN
            HMAX=HTICE(JK)
            NNMAX=JK
          ENDIF
        ENDDO
C        DO N=1,NUMNP
C          IF(KODE(N).NE.1) THEN
C            DIFF=DIFF+(QQQ(N)-PSURF(N))**2
C            NDIFF=NDIFF+1
C          ENDIF
C          HTICE(N)=QQQ(N)
C          IF(HTICE(N).GT.HMAX) THEN
C            HMAX=HTICE(N)
C            NNMAX=N
C          ENDIF
C
C ....... BEGIN FLOTATION PART
C
C ....... DEAL WITH BED ABOVE SEA LEVEL AND SURFACE BELOW BED
C          IF(HTICE(N).LE.UNDEPB(N)) HTICE(N)=UNDEPB(N)
C
C ....... DEAL WITH BED BELOW SEA LEVEL 
C ........AND SURFACE BELOW FLOTATION HEIGHT
C          IF(DEPB(N).LE.SEALEV) THEN
C            FLOT=(1.D0-RATDEN)*(DEPB(N)-SEALEV)
CCC           IF(HTICE(N).LE.FLOT) HTICE(N)=FLOT
C            IF(HTICE(N).LE.FLOT) HTICE(N)=SEALEV
C            IF(BDROCK(N).LE.-9999.) HTICE(N)=SEALEV
C          ELSE
C            IF(HTICE(N).LE.DEPB(N)) HTICE(N)=DEPB(N)   
C          ENDIF
C ....... END FLOTATION PART ****
C        ENDDO
C
C ..... OUTPUT TIME STEP INFO TO SCREEN
        WRITE(*,*) 'NMAX SURF=', HMAX, 
     &              ' AT NODE',NNMAX,' DIFF=',
     &              SQRT(DIFF/DBLE(NDIFF))
C
C
C
        WRITE(*,*) ' END OF RUN'
C
C ..... END OF RUN, ENTER ADJUST TO GO ROUND AGAIN.
        WRITE(*,*) 'INPUT 1 TO ENTER ADJUST -9 TO SKIP AND END'
        IF(IFLUSH) CALL GFLUSH
        READ(*,4000,END=999) IADJ
        WRITE(99,4000) IADJ
        IF(IADJ.EQ.'-9') GOTO 999
C
C ..... ENTER INTERACTIVE DATA SET MANIPULATOR
        CALL ADJUST(HED, NUMNP, NUMEL, X, Y, HTICE, ADOT,
     &              ADOTB,FRACT,BMELT,WTHICK,
     &              TEMP, ITYPE, TBED, PSURF, UNDEPB,
     &              BDROCK, DEPB, FLOWA, SLDGB, THICK, KX, CONST, 
     &              AFUDGE,NNODE, KODE, GEOFLUX, FLUX,
     &              HFIT, NUMCOL, NUMLEV, NUMGBC, NDT, INTER, DT,
     &              IBFLUX, BFLUX, NMAX, IDT, SLOPN, AMASS, TIME,
     &              NTSTEP, TTIME, VVOL, AAREA,TTBOT,TTAVG, TTNSL,
     &              TTSEAL, IFIT,
     &              IPLOT, CFACTOR, ACON,ICON, 
     &              ALPHAC,TBASE,
     &              IMELT,TWATER,PWATER,CALV,
     &              NTYPE,AADOT,AFRACT,ABDRCK,PPSURF,
     &              AFLOWA,ASLDGB,PCALV,ADC,WWWMIN,WWW,WRATE,WDIFF,
     &              WWWORIG,ADVANCE,HIRES,ACC,ABLAT,TFRACT)
        DTLOCAL=DT
C
c        DO I=1,NUMNP 
c          QQQ(I)=HTICE(I)
c          TTT(I)=HTICE(I)
c        ENDDO
C
C ..... LOAD NODAL PROPERTIES INTO ELEMENT PROPERTIES
        CALL CALVING(NMAX, NUMEL, NTYPE, KX, HTICE, DEPB, ADOT, 
     &           AADOT,CFACTOR,ADC,CALV,PCALV,AMASS(9))
        CALL ELPROP(NMAX, NUMEL, NTYPE, KX, ADOT, 
     &             AADOT, FRACT,
     &             AFRACT, BDROCK, ABDRCK, PSURF, PPSURF, FLOWA,ACON,
     &             ICON, AFLOWA, SLDGB, ASLDGB,CALV,PCALV)
C
        WRITE(*,*) 'INPUT 1 TO CONTINUE WITH NEW SET, 0 TO STOP '
        IF(IFLUSH) CALL GFLUSH
        READ(*,4000,END=999) IADJ
        WRITE(99,4000) IADJ
        IF(IADJ.EQ.'0') GOTO 999
        ITER=1
C ..... GOTO BEGINNING OF MAIN LOOP
C
C
C
      GOTO 215
C
C
C
999   CONTINUE
C
C ... VERBOSE WRITER OFF
C     CALL WRITER(HED, NUMNP, NUMEL, X, Y, HTICE, ADOT, FRACT, PSURF,
C    &            RHOI, BDROCK, FLOWA, SLDGB, THICK, KX, CONST, NNODE,
C    &            KODE, HFIT, NUMCOL, NUMLEV, NUMGBC, NDT, INTER, DT,
C    &            IBFLUX, BFLUX, AADOT, AFRACT, AFLOWA, ABDRCK,
C    &            ASLDGB)
C
      WRITE(18,2000) -99999.,2.,0,HED
C
C ... FORMAT STATEMENTS
2000  FORMAT(10X,G13.6,2X,G13.6,I13,/,A80)
4000  FORMAT(A2)
C
      STOP
      END
C===========================================================
      SUBROUTINE EISMINT(TIME,X,Y,HTICE,ADOT,ADOTB,DEPB,THICK,TBED,VOL,
     &                 AREA,GRIDAREA,VEL,SLOPN,ADP,ADM,ADC,BMELT,
     &                 NUMNP,NUMEL,KX)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM,MMAX=40)
      DIMENSION HTICE(NMAX),ADOT(NMAX),DEPB(NMAX),THICK(NMAX)
      DIMENSION X(NMAX),Y(NMAX),TBED(NMAX),ADOTB(NMAX)
      DIMENSION VEL(NMAX,3),SLOPN(4,NMAX),ADP(NMAX),ADM(NMAX)
      DIMENSION ADC(NMAX),BMELT(NMAX),KX(NMAX,4)
      DIMENSION TOUT(100)
      COMMON /TEMPERATURE/ TEMP(MMAX,NMAX)
      COMMON /LINE/ NP,NLINE(1000)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      logical file97,file98,file20
      SAVE IPASS,IGRIP,NX,NY,IRES,NTIMES,TOUT
      DATA IPASS /0/
C ... EISMINT OFF
      !RETURN
      CALL V2DFIELD(TIME,X,Y,HTICE,DEPB,SLOPN,ADOT,IRES)
      RETURN
C ...............
      IF(IPASS.EQ.0) THEN
        IPASS=1
        NTIMES=0
        inquire(file='fort.97',exist=file97)
        if(file97) then
          READ(97,*,END=102) NTIMES
          READ(97,*) (TOUT(IP),IP=1,NTIMES)
        else
          IF(IOTOGG) THEN
            WRITE(LIST(IPAGE+1),*) 'NO GRIP TIMES DEFINED'
            IPAGE=IPAGE+1
          ENDIF
        ENDIF
102     CONTINUE
        IGRIP=1
        NX=0
        NY=0
        IRES=0
        inquire(file='fort.98',exist=file98)
        if(file98) then
          READ(98,*,END=108) IGRIP,NX,NY,IRES
        else
          WRITE(LIST(IPAGE+1),*) 'NO GRIP POINT DEFINED'
          IPAGE=IPAGE+1
        ENDIF
108     CONTINUE
        NP=0
        inquire(file='fort.20',exist=file20)
        if(file20) then
          REWIND 20
          READ(20,1010,END=104) NP
1010    FORMAT(13I6)
          READ(20,1010) (NLINE(IP),IP=1,NP)
        else
          PRINT *,'NO GRIP LINE DEFINED'
        endif
104     CONTINUE
      ENDIF
C      WRITE(*,*) (TOUT(IP),IP=1,NTIMES)
C      WRITE(*,*) IGRIP
C      WRITE(*,1010) (NLINE(IP),IP=1,NP)
      DO IP=1,NTIMES
        IF(TIME.EQ.TOUT(IP)) THEN
C
          PRINT *,'EISMINT: HORIZONTAL 2D FIELD HF',TIME
          CALL H2DFIELD(TIME,NX,NY,IRES,HTICE,THICK,DEPB,TBED,
     &                  TEMP,ADP,ADM,VEL,NUMNP,NUMEL,KX)
C
          PRINT *,'EISMINT: VERTICAL   2D FIELD VF'
          CALL V2DFIELD(TIME,X,Y,HTICE,DEPB,SLOPN,ADOT,IRES)

          GOTO 106
        ENDIF
      ENDDO
106   CONTINUE
      IF(NTIMES.LT.0) THEN
        IF(MOD(REAL(TIME),REAL(ABS(NTIMES))).EQ.0.) THEN
          PRINT *,'EISMINT: VERTICAL   2D FIELD VF'
          CALL V2DFIELD(TIME,X,Y,HTICE,DEPB,SLOPN,ADOT,IRES)
        ENDIF
      ENDIF
      IF(TIME.GT.-130000.) THEN
        IF(TIME.LT.0.) THEN
          IF(MOD(TIME,100.D0).EQ.0.) THEN
            PRINT *,'EISMINT: (-) TIME DEPENDENT VARIABLE TV,GR',TIME
            HTMAX=FINDMAX(HTICE,NX,NY,IMAX,JMAX)
            AMELT=FMELT(HTICE,DEPB,TBED,NX,NY)
            AMELT=AMELT*GRIDAREA
C            PRINT *,AMELT/1E12
C ORIGINAL SPECS
C            WRITE(96,1003) TIME,AREA,VOL,HTMAX,IMAX,JMAX,AMELT,
C     &                       HTICE(IGRIP),THICK(IGRIP),TEMP(MMAX,IGRIP),
C     &                       TEMP(1,IGRIP),ADOT(IGRIP)
C 1003        FORMAT(1X,1P4G14.6,2I6,6G14.6)
C NEW SPECS ...
            RNINE=999.999D0
            WRITE(96,1003) NINT(TIME),AREA,VOL,RNINE,RNINE,
     &                     RNINE,RNINE,AMELT,
     &                     HTMAX,IMAX,JMAX
            WRITE(93,1004) NINT(TIME),HTICE(IGRIP),THICK(IGRIP),
     &                     TEMP(MMAX,IGRIP),TEMP(1,IGRIP),
     &                     ADOT(IGRIP)
1003        FORMAT(1X,I7,8(1X,E12.5),2I4)
1004        FORMAT(1X,I7,4(1X,F9.3),1X,E12.5)
          ENDIF
        ELSE
          IF(MOD(TIME,10.D0).EQ.0.) THEN
            IF(IOTOGG) THEN
              WRITE(LIST(IPAGE+1),*)
     &          'EISMINT: (+) TIME DEPENDENT VARIABLE TV,GR',TIME
              IPAGE=IPAGE+1
            ENDIF
            HTMAX=FINDMAX(HTICE,NX,NY,IMAX,JMAX)
            AMELT=FMELT(HTICE,DEPB,TBED,NX,NY)
            AMELT=AMELT*GRIDAREA
C            PRINT *,AMELT/1E12
C ORIGINAL SPECS
C            WRITE(96,1003) TIME,AREA,VOL,HTMAX,IMAX,JMAX,AMELT,
C     &                       HTICE(IGRIP),THICK(IGRIP),TEMP(MMAX,IGRIP),
C     &                       TEMP(1,IGRIP),ADOT(IGRIP)
C 1003        FORMAT(1X,1P4G14.6,2I6,6G14.6)
C NEW SPECS ...
            WRITE(96,1003) NINT(TIME),AREA,VOL,999.999,999.999,
     &                     999.999,999.999,AMELT,
     &                     HTMAX,IMAX,JMAX
            WRITE(93,1004) NINT(TIME),HTICE(IGRIP),THICK(IGRIP),
     &                     TEMP(MMAX,IGRIP),ADOT(IGRIP)
          ENDIF
        ENDIF        
      ENDIF        
      END
C=================================================
      FUNCTION FINDMAX(R,NX,NY,IMAX,JMAX)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM)
      DIMENSION R(NMAX)
      VAL=-1D30
      IC=0
      DO IX=1,NX
        DO IY=1,NY
          IC=IC+1
          IF(R(IC).GT.VAL) THEN
            VAL=R(IC)
            IMAX=IX
            JMAX=IY
          ENDIF
        ENDDO
      ENDDO
      FINDMAX=VAL
      END
C=================================================
      FUNCTION FMELT(HTICE,DEPB,TBED,NX,NY)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM,NSMAX=MAXTIME)
      DIMENSION HTICE(NMAX),DEPB(NMAX),TBED(NMAX)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      DATA RHOI /0.917D0/, RHOW /1.092D0/, THRESH /1.D0/
      IVAL=0
      IC=0
      DO IX=1,NX
        DO IY=1,NY
          IC=IC+1
          IF(DEPB(IC).LT.SEALEV) THEN
            FLOT=(1.D0-RATDEN)*(DEPB(IC)-SEALEV)
          ELSE
            FLOT=DEPB(IC)
          ENDIF
          IF(HTICE(IC).GT.FLOT+THRESH) THEN
            IF(TBED(IC).GE.0.0) THEN
              IVAL=IVAL+1
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      VAL=DBLE(IVAL)/DBLE((NX-1)*(NY-1))
      IF(IOTOGG) THEN
        WRITE(LIST(IPAGE+1),*) 'FMELT=',VAL
        IPAGE=IPAGE+1
      ENDIF
      FMELT=VAL
      END
C===========================================================
      SUBROUTINE H2DFIELD(TIME,NX,NY,IRES,HTICE,THICK,DEPB,TBED,
     &                    TEMP,ADP,ADM,VEL,NUMNP,NUMEL,KX)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM )
      DIMENSION HTICE(NMAX),THICK(NMAX),DEPB(NMAX),TBED(NMAX)
      DIMENSION TEMP(NMAX)
      DIMENSION ADP(NMAX),ADM(NMAX),VEL(NMAX,3),KX(NMAX,4)
      DIMENSION VELNODX(NMAX),VELNODY(NMAX),INUM(NMAX)
      DO I=1,NUMNP
        VELNODX(I)=0.0D0
        VELNODY(I)=0.0D0
        INUM(I)=0
      ENDDO
      DO I=1,NUMEL
        DO J=1,4
          LM=KX(I,J)
          VELNODX(LM)=VELNODX(LM)+VEL(I,2)
          VELNODY(LM)=VELNODY(LM)+VEL(I,3)
          INUM(LM)=INUM(LM)+1
        ENDDO
      ENDDO
      DO I=1,NUMNP
        IF(INUM(I).EQ.0) THEN
          VELNODX(I)=0.0D0
          VELNODY(I)=0.0D0
        ELSE
          VELNODX(I)=VELNODX(I)/INUM(I)
          VELNODY(I)=VELNODY(I)/INUM(I)
        ENDIF
      ENDDO
      WRITE(94,1000) 'TIME=',TIME
1000  FORMAT(1X,A,F15.0)
      WRITE(94,1001) NX,NY,IRES
      WRITE(94,*)
1001  FORMAT(3I7)
      IC=0
      DO IX=1,NX
        DO IY=1,NY
          IC=IC+1
          WRITE(94,1002) HTICE(IC),THICK(IC),DEPB(IC),
     &                   TBED(IC),TEMP(IC),
     &                   ADM(IC),ADP(IC),
     &                   VELNODX(IC),VELNODY(IC)
1002      FORMAT(1X,1P9G14.6)
        ENDDO
      ENDDO
      END
C===========================================================
      SUBROUTINE V2DFIELD(TIME,X,Y,HTICE,DEPB,SLOPN,ADOT,IRES)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM,MMAX=40,NSMAX=MAXTIME)
      DIMENSION X(NMAX),Y(NMAX),HTICE(NMAX),DEPB(NMAX)
      DIMENSION ADOT(NMAX),SLOPN(4,NMAX)
      DIMENSION LMAP(MMAX),XXX(MMAX),RHO(MMAX),TTTT(MMAX)
      DIMENSION TLEV(MMAX),XLEV(MMAX)
      DIMENSION AT(MMAX),UX(MMAX),UY(MMAX),UZ(MMAX),UC(NMAX)
      DIMENSION DDATA(9)
      COMMON /TEMPERATURE/ TEMP(MMAX,NMAX)
      COMMON /LINE/ NP,NLINE(1000)
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      DATA RHOI /0.917D0/, RHOW /1.092D0/, THRESH /1.D0/
      !WRITE(95,1000) 'TIME=',TIME
1000  FORMAT(1X,A,F15.0)
      XBASE=X(NLINE(1))
      YBASE=Y(NLINE(1))
      DO K=1,NP
        I=NLINE(K)
        DIST=SQRT((X(I)-XBASE)**2+(Y(I)-YBASE)**2)
c       DIST=(K-1)*IRES
        IF(DEPB(I).LT.SEALEV) THEN
          FLOT=(1.D0-RATDEN)*(DEPB(I)-SEALEV)
        ELSE
          FLOT=DEPB(I)
        ENDIF
        IF(.true. .or. HTICE(I).GT.FLOT+THRESH) THEN
          LTOT=0
          DO II=1,MMAX/2
            LMAP(II)=II
C           LMAP(II)=1
            LTOT=LTOT+LMAP(II)
          ENDDO
          DO II=MMAX/2+1,MMAX
            LMAP(II)=LMAP(II-1)-1
C           LMAP(II)=1
            LTOT=LTOT+LMAP(II)
          ENDDO
C          XXX(1)=THICK(I)
C          XXX(MMAX)=0.
          XXX(1)=HTICE(I)
          XXX(MMAX)=DEPB(I)
          DX=(XXX(MMAX)-XXX(1))/DBLE(LTOT)
          DO II=2,MMAX
           XXX(II)=XXX(II-1)+LMAP(II-1)*DX
          ENDDO
          DO II=1,MMAX
            DEPTH=XXX(1)-XXX(II)
            RHO(II)=DENSITY(DEPTH)
          ENDDO
          DO II=1,MMAX
            XLEV(II)=XXX(II)
            TLEV(II)=TEMP(II,I)
            TTTT(II)=TEMP(II,I)
          ENDDO
C ....... X VELOCITIES .........................................
          CALL VELO(TIME,MMAX,SLOPN(2,I),XXX,TTTT,RHO,AT,UX,UC,
     &              BMELT,DDATA,AEFF)                                                         
C ....... Y VELOCITIES .........................................
          CALL VELO(TIME,MMAX,SLOPN(3,I),XXX,TTTT,RHO,AT,UY,UC,
     &              BMELT,DDATA,AEFF)  
C ....... Z VELOCITIES .........................................
          UZ(1)=ADOT(I)
          UZ(MMAX)=0.D0
C ......  VERTICAL VELOCITY DISTRIBUTION: LINEAR
          SLOPE=(UZ(MMAX)-UZ(1))/(XXX(MMAX)-XXX(1))
          DO II=2,MMAX-1
            UZ(II)=UZ(1)+SLOPE*(XXX(II)-XXX(1))
          ENDDO
        ELSE
          DO II=1,MMAX
            XXX(II)=DEPB(I)
            UX(II)=-99999.0D0
            UY(II)=-99999.0D0
            UZ(II)=-99999.0D0
            TTTT(II)=-99999.0D0
          ENDDO
        ENDIF
!       DO J=2,MMAX
!         WRITE(38,*) 0.001*DIST,XLEV(J),TLEV(J)
!       ENDDO
        DO J=1,MMAX
C          WRITE(95,1002) X(I),XXX(J),UX(J),UY(J),UZ(J),TEMP(J,I)
!          WRITE(95,1002) DIST,XXX(J),UX(J),UY(J),UZ(J),
!    &                   TTTT(J)
1002    FORMAT(1X,1P8G14.6)
        ENDDO
      ENDDO
      END
C=======================================================
      SUBROUTINE MATMULT(NMAX,NZ,NZ1,NUMNP,A,KA,KZ,ADIAG,
     &                   QQQ,BOLD,FLUX)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(NMAX,NZ)
      DIMENSION KA(NMAX,NZ1),KZ(NMAX),QQQ(NMAX),BOLD(NMAX)
      DIMENSION FLUX(NMAX),ADIAG(NMAX)
C$DOACROSS LOCAL(K,I,TEMP)
      DO K=1,NUMNP
        TEMP=ADIAG(K)*QQQ(KA(K,1))
        DO I=2,KZ(K)
          TEMP=TEMP+A(K,I)*QQQ(KA(K,I))
c          if(abs(a(k,i)).gt.1e-5) then
c            print *,k,ka(k,i),real(A(K,I)),real(QQQ(KA(K,I))),
c     &            real(temp)
c          endif
        ENDDO
        FLUX(K)=TEMP-BOLD(K)
C        PRINT *,K,REAL(FLUX(K))
      ENDDO
c      pause
      END
C=======================================================
      SUBROUTINE FIND4(XI,YI,NUMNPL,XL,YL,QL,QI)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM )
      REAL*4 XL(NMAX),YL(NMAX),QL(NMAX)
      DIMENSION XMIN(4),YMIN(4),DISTMIN(4),IMIN(4)
      DO I=1,4
        DISTMIN(I)=1D30
        IMIN(I)=-999
      ENDDO
      DO N=1,NUMNPL
        DIST=(XI-XL(N))**2+(YI-YL(N))**2
        IF(XL(N).GT.XI .AND. YL(N).GT.YI) THEN
C QUAD-1
          IQUAD=1
        ELSEIF(XL(N).LE.XI .AND. YL(N).GT.YI) THEN
C QUAD-2
          IQUAD=2
        ELSEIF(XL(N).LE.XI .AND. YL(N).LE.YI) THEN
C QUAD-3
          IQUAD=3
        ELSEIF(XL(N).GT.XI .AND. YL(N).LE.YI) THEN
C QUAD-4
          IQUAD=4
        ENDIF
        IF(DIST.LT.DISTMIN(IQUAD)) THEN
          DISTMIN(IQUAD)=DIST
          IMIN(IQUAD)=N
        ENDIF
      ENDDO
      DO I=1,4
        IF(IMIN(I).EQ.-999) THEN
C          PRINT *,'PROBLEMS IN FIND4',I
          QI=-999.
          RETURN
        ENDIF
      ENDDO
      QI=0
      WTOT=0.
      DO I=1,4
        IF(DISTMIN(I).EQ.0) THEN
          QI=QL(IMIN(I))
          GOTO 100
        ENDIF
        WEIGHT=1./DISTMIN(I)
        QI=QI+QL(IMIN(I))*WEIGHT
        WTOT=WTOT+WEIGHT
      ENDDO
      QI=QI/WTOT
100   CONTINUE
C      PRINT *,REAL(XI/1000),REAL(YI/1000),REAL(QI)
C      DO I=1,4
C        PRINT *,IMIN(I),XL(IMIN(I))/1000,YL(IMIN(I))/1000,
C     &          SQRT(REAL(DISTMIN(I)))/1000,
C     &  REAL(QL(IMIN(I)))
C      ENDDO
      END
C=======================================================
      SUBROUTINE FINDMIN(NUMNP,XI,YI,NUMNPL,XL,YL,IMIN)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 XL(NUMNPL),YL(NUMNPL)
      DISTMIN=1E30
      IMIN=-999
      DO N=1,NUMNPL
        DIST=(XI-XL(N))**2+(YI-YL(N))**2
        IF(DIST.LT.DISTMIN) THEN
          DISTMIN=DIST
          IMIN=N
        ENDIF
      ENDDO
      IF(IMIN.EQ.-999) THEN
        PRINT *,'PROBLEMS IN FINDER'
        STOP
      ENDIF
C      PRINT *,REAL(XI),REAL(YI),REAL(XL(IMIN)),REAL(YL(IMIN))
C      PRINT *,IMIN,SQRT(DISTMIN)
      END
C=======================================================
      SUBROUTINE FLUXUP(NMAX,NUMNP,NUMGBC,IBFLUX,BFLUX,FLUX,
     &                  X,Y,KODE,TIME)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IBFLUX(NMAX,2),BFLUX(NMAX),FLUX(NMAX),KODE(NMAX)
      DIMENSION X(NMAX),Y(NMAX)
      DO N=1,NUMGBC
        IF(KODE(IBFLUX(N,1)).EQ.1 .AND. 
     &     KODE(IBFLUX(N,2)).EQ.1) THEN
          I=IBFLUX(N,1)
          J=IBFLUX(N,2)
          XLENGTH=SQRT((X(J)-X(I))**2+(Y(J)-Y(I))**2)
          BFLUXN=0.5D0*(FLUX(I)+FLUX(J))/XLENGTH
C          PRINT *,N,(IBFLUX(N,J),J=1,2),REAL(BFLUX(N)),REAL(BFLUXN)
          BFLUX(N)=BFLUXN
C A STRANGE FIX APPLIED HERE FOR REASONS UNKNOWN...
C          XLENGTH=SQRT(XL)
          BFLUX(N)=-XLENGTH*BFLUX(N)/(400.D0*10000.D0)
        ENDIF
      ENDDO
      END
C=======================================================
      function PDDCOUNT(temp,ampli)
      IMPLICIT REAL*8(A-H,O-Z)
      pi = 4d0 * atan(1d0)
      PDDCOUNT = 0
      if(temp+0.5*ampli .lt. 273.16) return
      ndays = 687
      const = 2*pi/ndays
      do i = 1, ndays/2
        t=i
        tmp = temp + 0.5 * ampli * sin(const*t)
        if(tmp.gt.273.16) then
          PDDCOUNT = PDDCOUNT + (tmp - 273.16)
          !PDDCOUNT = PDDCOUNT + 1
        endif
        !print *,temp,tmp,PDDCOUNT
      enddo
      !if(PDDCOUNT.gt.0) print *,temp,PDDCOUNT
      end
C=======================================================
      SUBROUTINE ACCSET(NUMNP,AMASS,ACOMSAVE,TIME,
     &                  THICK,IDT,HTICE,BDROCK,PSURF,SLOPN,X,Y,TEMP,
     &                  ADOT,ADOTB,ADP,ADM,TJUNK,TFRACT,ACC,ABLAT)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM)
      DIMENSION AMASS(11),THICK(NMAX),IDT(NMAX),HTICE(NMAX),BDROCK(NMAX)
      DIMENSION SLOPN(4,NMAX),X(NMAX),Y(NMAX),TEMP(NMAX),ADOT(NMAX)
      DIMENSION ADOTB(NMAX),ADP(NMAX),ADM(NMAX),PSURF(NMAX)
      DIMENSION ACC(NMAX),ABLAT(NMAX)
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      COMMON /TCONST/ TSORIG(NMAX),TFLAG
      LOGICAL TFLAG
      EXTERNAL WARMING
      !print *,'in ACCSET';pause
C**********EISMINT ALA HUYBRECHTS ********************
      WRM=WARMING(AMASS(9)+TNSLBASE)
C*****************************************************
C$DOACROSS LOCAL(I,ADOTI,AJUNK), SHARED(TIME)
      TAVGSURF=0.
      TMINSURF=1E30
      TMAXSURF=-1E30
      DO I=1,NUMNP
C*****************************************************
c-------experimental for icy highland--------
C*****************************************************
        if(.false.) then
          ablat(I) = .05 * 0
          if(htice(I).gt.0) then
            acc(I)=min(htice(I)/22000.,ablat(I)+.01)
          else
            acc(I)=0.0
          endif
        elseif(.false.) then
          hlow=-3000;hela=1000;abl=-0.005;accmax=0.01
          hlow=-3000;hela=3000;abl=-0.005;accmax=0.01
          azzz=-abl/(hela-hlow)
          bzzz=hlow*abl/(hela-hlow)
          ablat(I)=-abl
          acc(I) = max(0.d0,azzz*htice(I) + bzzz)
          acc(I) = min(acc(I), ablat(I)+accmax)

          ablat(I) = 0.01
          acc(I) =   0.005

        elseif(.false.) then ! turn on for icy highland experiment
          TNSL = AMASS(9)
          hlow=-3000;hela=1000;high=2000;abl=0.005;accmax=0.01
          ablat(I)=abl
          if(.true.) then
            hlow = hlow - 1000 * (TNSL + 50) / ACOM
            hela = hela - 1000 * (TNSL + 50) / ACOM
            high = high - 1000 * (TNSL + 50) / ACOM
          endif
          if(.false.) then
            PDD  = PDDCOUNT(temp(i)+273.16d0,40.d0)
            !ablat(I) = ablat(I) + 0.00144 * PDD
            ablat(I) = ablat(I) + 0.00144 * PDD *0.75 ! reduced sun
          endif
          if(htice(i).lt.hlow) then
            acc(I)=0.0
          elseif(htice(i).lt.hela) then
            acc(I)=abl*(htice(I)-hlow)/(hela-hlow)
          elseif(htice(i).lt.high) then
            acc(I)=abl+accmax*(htice(I)-hela)/(high-hela)
          else
            acc(I)=abl+accmax
          endif
        endif
C*****************************************************
c-------experimental for icy highland--------
C*****************************************************
        ADOT(I)=TFRACT*ACC(I)-ABLAT(I)
        !ADOT(I) = ADOTB(I)
        !print *,i,real(adot(i)),real(acc(i)),real(ablat(i))
      ENDDO
C*****************************************************
c************* nothing is used beyond here ***********
C*****************************************************
      RETURN
C*****************************************************
C*****************************************************
      DO I=1,NUMNP
C ...... FOLLOWING MAKES LAPSE RATE DIFFERENT FOR BARE GROUND...
C        IF(BDROCK(I).GT.0) THEN
C          THL=HTICE(I)-BDROCK(I)
C          IF(THL.GT.1) THEN
C            ACOM=ACOMSAVE
C          ELSE
C            ACOM=ACOMSAVE/2.
C          ENDIF
C        ENDIF
        ADOTI=AFUNCT(TIME, IDT(I), AMASS,
     &               HTICE(I), BDROCK(I),PSURF(I),
     &               SLOPN(1,I),X(I),Y(I),TEMP(I))
        TAVGSURF=TAVGSURF+TEMP(I)
        TMINSURF=MIN(TMINSURF,TEMP(I))
        TMAXSURF=MAX(TMAXSURF,TEMP(I))
        IF(IDT(I).GT.0) THEN
          ADOT(I)=ADOTI
        ELSE
          IF(ITOGG) THEN
            AJUNK=ACCUM(IDT(I),X(I)*.001D0,Y(I)*.001D0,HTICE(I),
     &                  SLOPN(1,I),
     &                  PSURF(I),AMASS(9),TEMP(I))
            IF(TFLAG) TEMP(I)=TSORIG(I)
C           ADOT(I)=ADOTB(I)-ABL*.01
C*****************************************************
C**********EISMINT ALA HUYBRECHTS ********************
            ADP(I)=WRM*ADOTB(I)
            ADM(I)=-ABL*0.01D0
C STRAIGHT
c           ADOT(I)=WRM*ADOTB(I)-ABL*.01D0
            ADOT(I)=AMARS(WRM,ADOTB(I),ABL)
C SLOW DOWN ACCUM RATE CHANGE...-------------------------
C           ADOT(I)=0.5D0*(ADOT(I)+WRM*ADOTB(I)-ABL*.01D0)
C--------------------------------------------------------
C****************************************** ***********
c         ELSEIF(TFLAG) THEN  
c           print *,adotb(i),amass(9)
c           ADOT(I)=ADOTB(I)+AMASS(9)
          ELSE
            ADOT(I)=ADOTB(I)
          ENDIF
        ENDIF
C ......  FOLLOWING ZEROS ACCUM FOR STEEP SLOPE, NO ICE...
C        CALL STEEPSET(I,HTICE(I),BDROCK(I),SLOPN(1,I),ADOT(I))
        ENDDO
      TAVGSURF=TAVGSURF/NUMNP
      CALL MILANK(TIME,0.D0,0.D0,0.D0,0.D0,0.D0,TNSL,TJUNK,RJUNK)
      avg=0
      do i=1,numnp
        avg=avg+adot(i)
      enddo
      TJUNK=avg/numnp
      IF(IOTOGG) THEN
        WRITE(LIST(IPAGE+1),*) 'TAVGSURF=',
     &        REAL(TMINSURF),REAL(TAVGSURF),REAL(TMAXSURF),
     &        REAL(TJUNK)
        IPAGE=IPAGE+1
      ENDIF
      TJUNK=TAVGSURF
      END
C=======================================================
      SUBROUTINE STEEPSET(I,HTI,BEDI,SLOPNI,ADOTI)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(BEDI.LT.0) RETURN
      IF(HTI-BEDI.LT.1. .AND. ABS(SLOPNI).GT. 0.005) THEN
        ADOTI=-10.
C        PRINT *,'SETTING:',I,' -1.',REAL(HTI-BEDI),REAL(SLOPNI)
      ENDIF
      END
C=======================================================
      SUBROUTINE DUMPBC1(TIME,NUMNP,X,Y,QQQ,FLUX,KODE,
     &                  NUMGBC,IBFLUX,BFLUX)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM)
      DIMENSION X(NMAX),Y(NMAX),QQQ(NMAX),KODE(NMAX),FLUX(NMAX)
      REAL*4 XL(NMAX),YL(NMAX),QL(NMAX),FLUXL(NMAX)
      DIMENSION IBFLUX(NMAX,2),BFLUX(NMAX)
      LOGICAL DUMPOUT,DUMPIN
      logical file69
      DATA IPASS /0/
      SAVE IPASS,DUMPOUT,DUMPIN,XL,YL
C ... FLUX IS BUGGERED UP SO IT IS ALL TURNED OFF. FIXED BC WORKS...
      IF(IPASS.EQ.0) THEN
        IPASS=1
        IDUMP=0
        PRINT *,' DO YOU WANT TO 1: WRITE CONFIG FOR LATER USE.'
        PRINT *,'                2: READ CONFIG PREVIOUSLY WRITTEN'
        PRINT *,'                3: WRITE/READ CONFIG'
        PRINT *,'                0: DO NOTHING WITH CONFIG'
        READ(*,*) IDUMP
        WRITE(99,*) IDUMP,'     DUMP:1-WRITE,2-READ,3-W/R,0-NOTHING'
        IF(IDUMP.EQ.1) THEN
          DUMPOUT=.TRUE.
          DUMPIN=.FALSE.
          WRITE(70) TIME,NUMNP
          WRITE(70) (REAL(X(I)),I=1,NUMNP)
          WRITE(70) (REAL(Y(I)),I=1,NUMNP)
        ELSEIF(IDUMP.EQ.2) THEN
          DUMPOUT=.FALSE.
          DUMPIN=.TRUE.
          inquire(file='fort.69',exist=file69)
          if(file69) then
            READ(69,END=999) TIMEL,NUMNPL
c           print *, TIMEL,NUMNPL
            READ(69) (XL(I),I=1,NUMNPL)
            READ(69) (YL(I),I=1,NUMNPL)
          else
            goto 999
          endif
        ELSEIF(IDUMP.EQ.3) THEN
          DUMPOUT=.TRUE.
          DUMPIN=.TRUE.
          inquire(file='fort.69',exist=file69)
          if(file69) then
            READ(69,END=999) TIMEL,NUMNPL
c           print *, TIMEL,NUMNPL
            READ(69) (XL(I),I=1,NUMNPL)
            READ(69) (YL(I),I=1,NUMNPL)
          else
            goto 999
          endif
          WRITE(70) TIME,NUMNP
          WRITE(70) (REAL(X(I)),I=1,NUMNP)
          WRITE(70) (REAL(Y(I)),I=1,NUMNP)
        ELSE
          DUMPOUT=.FALSE.
          DUMPIN=.FALSE.
        ENDIF
      ENDIF
      IF(DUMPOUT) THEN
        WRITE(70) TIME,NUMNP
        WRITE(70) (REAL(QQQ(I)),I=1,NUMNP)
C       WRITE(70) (REAL(FLUX(I)),I=1,NUMNP)
      ENDIF
      IF(DUMPIN) THEN
        READ(69,END=999) TIMEL,NUMNPL
c          print *, TIMEL,NUMNPL
        READ(69) (QL(I),I=1,NUMNPL)
C       READ(69) (FLUXL(I),I=1,NUMNPL)
        IF(TIME.NE.TIMEL) THEN
          PRINT *,'PROBLEMS IN DUMPBC',TIME,TIMEL
          STOP
        ENDIF
        DO I=1,NUMNP
          IF(KODE(I).EQ.1) THEN
C ......... FIND XL,YL NEAREST TO X,Y .........
C ......... SIMPLE FIND THE CLOSEST
            IF(.FALSE.) THEN
              CALL FINDMIN(NUMNP,X(I),Y(I),NUMNPL,XL,YL,II)
              QQQ(I)=QL(II)
            ELSE
C ........... MORE COMPLEX, FIND THE CLOSEST 4 (1 IN EACH QUADRANT) 
C             AND INTERPOLATE W/
C             DISTANCE-WEIGHTED MEAN
              XI=X(I)
              YI=Y(I)
              CALL FIND4(XI,YI,NUMNPL,XL,YL,QL,QI)
              IF(QI.NE.-999.) THEN
                QQQ(I)=QI
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        IF(.FALSE.) THEN
          DO N=1,NUMGBC
            I=IBFLUX(N,1)
            J=IBFLUX(N,2)
C ....  ... FIND XL,YL NEAREST TO X,Y .........
            CALL FINDMIN(NUMNP,X(I),Y(I),NUMNPL,XL,YL,II)
            CALL FINDMIN(NUMNP,X(J),Y(J),NUMNPL,XL,YL,JJ)
            XLENGTH=SQRT((X(J)-X(I))**2+(Y(J)-Y(I))**2)
            BFLUXN=0.5D0*(FLUXL(II)+FLUXL(JJ))/XLENGTH
C           PRINT *,N,(IBFLUX(N,J),J=1,2),REAL(BFLUX(N)),REAL(BFLUXN)
            BFLUX(N)=-BFLUXN*2
C A STRANGE FIX APPLIED HERE FOR REASONS UNKNOWN...
C           XLENGTH=SQRT(XL)
C           BFLUX(N)=XLENGTH*BFLUX(N)/(400.D0*10000.D0)
          ENDDO
        ENDIF
      ENDIF
999   CONTINUE
      END
C=======================================================
      SUBROUTINE DUMPBC(TIME,NUMNP,X,Y,QQQ,FLUX,KODE,
     &                  NUMGBC,IBFLUX,BFLUX)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM)
      DIMENSION X(NMAX),Y(NMAX),QQQ(NMAX),KODE(NMAX),FLUX(NMAX)
      REAL*4 XL(NMAX),YL(NMAX),Q1(NMAX),Q2(NMAX)
C     REAL*4 FLUX1(NMAX),FLUX2(NMAX)
      DIMENSION IBFLUX(NMAX,2),BFLUX(NMAX)
      LOGICAL DUMPOUT,DUMPIN
      DATA IPASS /0/
      SAVE IPASS,DUMPOUT,DUMPIN,XL,YL,Q1,Q2,TIME1,TIME2,NUMNPL
C ... FLUX IS BUGGERED UP SO IT IS ALL TURNED OFF. FIXED BC WORKS...
      IF(IPASS.EQ.0) THEN
        IPASS=1
        IDUMP=0
        PRINT *,' DO YOU WANT TO 1: WRITE CONFIG FOR LATER USE.'
        PRINT *,'                2: READ CONFIG PREVIOUSLY WRITTEN'
        PRINT *,'                3: READ/WRITE CONFIG'
        PRINT *,'                0: DO NOTHING WITH CONFIG'
        READ(*,*) IDUMP
        WRITE(99,*) IDUMP,'     DUMP:1-WRITE,2-READ,3-R/W,0-NOTHING'
        IF(IDUMP.EQ.1) THEN
          DUMPOUT=.TRUE.
          DUMPIN=.FALSE.
          WRITE(70) TIME,NUMNP
          WRITE(70) (REAL(X(I)),I=1,NUMNP)
          WRITE(70) (REAL(Y(I)),I=1,NUMNP)
        ELSEIF(IDUMP.EQ.2) THEN
          DUMPOUT=.FALSE.
          DUMPIN=.TRUE.
          READ(69,END=999) TIMEL,NUMNPL
c          print *, TIMEL,NUMNPL
          READ(69) (XL(I),I=1,NUMNPL)
          READ(69) (YL(I),I=1,NUMNPL)
          READ(69,END=999) TIME1,NUMNPL
c          print *, TIME1,NUMNPL
          READ(69) (Q1(I),I=1,NUMNPL)
C         READ(69) (FLUX1(I),I=1,NUMNPL)
          READ(69,END=999) TIME2,NUMNPL
c          print *, TIME2,NUMNPL
          READ(69) (Q2(I),I=1,NUMNPL)
C         READ(69) (FLUX2(I),I=1,NUMNPL)
C          WRITE(7,*) 'READ 1ST TWO:',REAL(TIME1),REAL(TIME2),REAL(TIME)
          IF(TIME.LT.TIMEL) THEN
C            PRINT *,'PROBLEMS IN DUMPBC',REAL(TIME),
C     &               REAL(TIMEL),REAL(TIME2)
            GOTO 999
          ENDIF
        ELSEIF(IDUMP.EQ.3) THEN
          DUMPOUT=.TRUE.
          DUMPIN=.TRUE.
          WRITE(70) TIME,NUMNP
          WRITE(70) (REAL(X(I)),I=1,NUMNP)
          WRITE(70) (REAL(Y(I)),I=1,NUMNP)
          READ(69,END=999) TIMEL,NUMNPL
c          print *, TIMEL,NUMNPL
          READ(69) (XL(I),I=1,NUMNPL)
          READ(69) (YL(I),I=1,NUMNPL)
          READ(69,END=999) TIME1,NUMNPL
c          print *, TIME1,NUMNPL
          READ(69) (Q1(I),I=1,NUMNPL)
C         READ(69) (FLUX1(I),I=1,NUMNPL)
          READ(69,END=999) TIME2,NUMNPL
c          print *, TIME2,NUMNPL
          READ(69) (Q2(I),I=1,NUMNPL)
C         READ(69) (FLUX2(I),I=1,NUMNPL)
C          WRITE(7,*) 'READ 1ST TWO:',REAL(TIME1),REAL(TIME2),REAL(TIME)
          IF(TIME.LT.TIMEL) THEN
C            PRINT *,'PROBLEMS IN DUMPBC',REAL(TIME),
C     &               REAL(TIMEL),REAL(TIME2)
            GOTO 999
          ENDIF
        ELSE
          DUMPOUT=.FALSE.
          DUMPIN=.FALSE.
        ENDIF
      ENDIF
      IF(DUMPOUT) THEN
        WRITE(70) TIME,NUMNP
        WRITE(70) (REAL(QQQ(I)),I=1,NUMNP)
C       WRITE(70) (REAL(FLUX(I)),I=1,NUMNP)
      ENDIF
      IF(DUMPIN) THEN
        DOWHILE(TIME.GE.TIME2)
C          WRITE(7,*) TIME,'>=',TIME2
          TIME1=TIME2
          DO I=1,NUMNPL
            Q1(I)=Q2(I)
C           FLUX1(I)=FLUX2(I)
          ENDDO
          READ(69,END=999) TIME2,NUMNPL
c          print *, TIME2,NUMNPL
          READ(69) (Q2(I),I=1,NUMNPL)
C         READ(69) (FLUX2(I),I=1,NUMNPL)
C          WRITE(7,*) 'READ ANOTHER:',REAL(TIME),
C     &                REAL(TIME1),REAL(TIME2)
        ENDDO
        FACTOR=(TIME-TIME1)/(TIME2-TIME1)
C        WRITE(7,*) ' FACTOR=',REAL(TIME1),REAL(TIME),
C     &               REAL(TIME2),REAL(FACTOR)
        IF(FACTOR.LT.0) GOTO 999
        DO I=1,NUMNP
          IF(KODE(I).EQ.1) THEN
C ......... FIND XL,YL NEAREST TO X,Y .........
C ......... SIMPLE FIND THE CLOSEST
            IF(.FALSE.) THEN
              CALL FINDMIN(NUMNP,X(I),Y(I),NUMNPL,XL,YL,II)
              QQQ(I)=((1.D0-FACTOR)*Q1(II)+FACTOR*Q2(II))
            ELSE
C ........... MORE COMPLEX, FIND THE CLOSEST 4 (1 IN EACH QUADRANT) 
C             AND INTERPOLATE W/
C             DISTANCE-WEIGHTED MEAN or WITH FEM INTERPOLATION
              XI=X(I)
              YI=Y(I)
              CALL FIND4A(XI,YI,NUMNPL,XL,YL,FACTOR,Q1,Q2,QI)
              IF(QI.NE.-999.) THEN
                QQQ(I)=QI
              ENDIF
            ENDIF
C            WRITE(7,*) I,II,REAL(Q1(II)),REAL(QQQ(I)),REAL(Q2(II))
          ENDIF
        ENDDO
      ENDIF
999   CONTINUE
      END
C=======================================================
      SUBROUTINE FIND4A(XI,YI,NUMNPL,XL,YL,FACTOR,Q1,Q2,QI)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM )
      REAL*4 XL(NMAX),YL(NMAX),Q1(NMAX),Q2(NMAX)
      REAL*4 XMIN(4),YMIN(4),DISTMIN(4),value,xtest,ytest
      DIMENSION IMIN(4)
      DO I=1,4
        DISTMIN(I)=1E30
        IMIN(I)=-999
      ENDDO
      DO N=1,NUMNPL
        DIST=(XI-XL(N))**2+(YI-YL(N))**2
        IF(XL(N).GT.XI .AND. YL(N).GT.YI) THEN
C QUAD-1
          IQUAD=1
        ELSEIF(XL(N).LE.XI .AND. YL(N).GT.YI) THEN
C QUAD-2
          IQUAD=2
        ELSEIF(XL(N).LE.XI .AND. YL(N).LE.YI) THEN
C QUAD-3
          IQUAD=3
        ELSEIF(XL(N).GT.XI .AND. YL(N).LE.YI) THEN
C QUAD-4
          IQUAD=4
        ENDIF
        IF(DIST.LT.DISTMIN(IQUAD)) THEN
          DISTMIN(IQUAD)=DIST
          IMIN(IQUAD)=N
        ENDIF
      ENDDO
      DO I=1,4
        IF(IMIN(I).EQ.-999) THEN
C          PRINT *,'PROBLEMS IN FIND4',I
          QI=-999.
          RETURN
        ENDIF
      ENDDO
      if(.false.) then
        QI=0
        WTOT=0.
        DO I=1,4
          IF(DISTMIN(I).EQ.0) THEN
            II=IMIN(I)
            QINT=((1.D0-FACTOR)*Q1(II)+FACTOR*Q2(II))
            QI=QINT
            GOTO 100
          ENDIF
          WEIGHT=1./DISTMIN(I)
          II=IMIN(I)
          QINT=((1.D0-FACTOR)*Q1(II)+FACTOR*Q2(II))
          QI=QI+QINT*WEIGHT
          WTOT=WTOT+WEIGHT
        ENDDO
        QI=QI/WTOT
100     CONTINUE
        PRINT *,REAL(XI/1000),REAL(YI/1000),REAL(QI)
        DO I=1,4
          PRINT *,IMIN(I),XL(IMIN(I))/1000,YL(IMIN(I))/1000,
     &          SQRT(REAL(DISTMIN(I)))/1000,
     &          REAL(Q1(IMIN(I))),REAL(Q2(IMIN(I)))
        ENDDO
      ELSE
        do i=1,4
          xmin(i)=xl(imin(i))
          ymin(i)=yl(imin(i))
          QINT=((1.D0-FACTOR)*Q1(imin(i))+FACTOR*Q2(imin(i)))
          distmin(i)=qint
        enddo
        xtest=xi
        ytest=yi
        call interp1(xmin,ymin,distmin,xtest,ytest,value)
        qi=value
      ENDIF
      END
c======================================================
      SUBROUTINE RESETTEMP(NUMNP,TNSL,HTICE,TEMP,AMASS,TIME)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM)
      COMMON /TCONST/ TSORIG(NMAX),TFLAG
      COMMON /COORD/ CLAT(NMAX),CLONG(NMAX)
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      LOGICAL TFLAG
      dimension TEMP(numnp),HTICE(numnp)
c***************************************A*****************
      dimension AMASS(11)
      COMMON /EXPER/ TPERIOD,TINIT,TVSTART,TVFINAL,IEXPER
      !print *,TFLAG;pause
      IF(IEXPER.EQ.23) THEN
C user defined variation experiment.
C      TPERIOD yr sin cycle tnsl
C      TVSTART  initial time
C      TVFINAL final time
C
        TNSL=COEF(TIME,REAL(TINIT),REAL(TINIT+TPERIOD),
     &                 REAL(TVSTART),REAL(TVFINAL))
        !AMASS(9) = TNSL
      ENDIF
      RLAPSE = ACOM
c********************************************************
      if(TFLAG) THEN
        !print *,'here in reset'
        do i=1,numnp
          temp(i)=tsorig(i)
        enddo
      else
        !print *,'here in reset else'
        radpdeg=4d0*atan(1d0)/180d0
        RLAPSE=ACOM;TBASE=TNSL;AMPLI=50
        if(.true.) then
c---------HIGH--------
          RLAPSE = -2.4
          TBASE  = -50
c****************************
          TBASE  = TNSL
c****************************
          AMPLI  =  20
        elseif(.true.) then
c---------MED--------
          RLAPSE = -2
          TBASE  = -100
          AMPLI  =  50
        elseif(.true.) then
c---------LOW--------
          RLAPSE = -0.5
          TBASE  = -115
          AMPLI  =  53
        endif
        do i=1,numnp
          temp(i)=TBASE+AMPLI*cos(CLAT(I)*radpdeg)
          temp(i)=temp(i)+RLAPSE*htice(i)/1000
          !write(29,*) htice(i),temp(i)+273.16,clat(i)
          !write(29,*) temp(i)+273.16,htice(i),clat(i)
        enddo
        !stop
      endif
      ACOM = RLAPSE
      end
c======================================================
      function AMARS(WRM,ADOTBI,ABL)
      IMPLICIT REAL*8(A-H,O-Z)
c     ADOTI=WRM*ADOTBI-ABL*1
c     ADOTI=sqrt(wrm)*ADOTBI-ABL*1
c     ADOTI=ADOTBI-wrm**0.333333*ABL
c     ADOTI=wrm**0.33333*(ADOTBI-ABL)
      ADOTI=(ADOTBI-ABL)
c     print *,adoti,adotbi,wrm,abl,idt(i)
      AMARS=ADOTI
      end
