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
c 1X case-------------------------------------------
c             VOLMAX = 2*34 ! current inventory GEL
c 1X case-------------------------------------------
              VOLMAX = VOLMAX/1000 ! GEL in km
              VOLMAX = VOLMAX * 144.8d6 ! km^3
! special setting, VOLMAX in GEL and then converted to m^3

              VOLMAX=VOLMAX*1d9   ! m^3
              if(.false.) then ! use this when setting flist
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
     &                            real(vratio)
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
      FUNCTION ACCUM(IDOT,X,Y,ELEV1,SLOPE1,PSURF1,TNSL,TS)                      
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM,NSMAX=MAXTIME)
      COMMON /TCONST/ TSORIG(NMAX),TFLAG                         
      LOGICAL TFLAG
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      DIMENSION QI(12),QS(12),TTT(12)                                   
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      COMMON /LOCAL/ SIGMA,TMAA,TMS,PI2
      DATA QI /960.d0,1036.d0,1200.d0,825.,330.d0,90.d0,150.d0,
     &         600.d0,1200.d0,1020.d0,    
     &         930.d0,850.d0/                                               
      DATA QS /-0.667d0,4.6d0,11.667d0,9.167d0,3.667d0,
     *         1.d0,1.667d0,6.667d0,12.d0,6.333d0,  
     &         0.333d0,-3.333d0/                                            
C     DATA AAA /-9.14/, BBB /-.68/, CCC /34.461/                        
C     DATA WWW /13.05/, XXX /.664/, ZZZ /2.608/                         
       DATA AAA / -9.62376690d0     /                                     
       DATA BBB /-0.546917617d0     /                                     
       DATA CCC /  24.9793854d0     /                                     
       DATA WWW /  19.1390686d0     /                                     
       DATA XXX / 0.922791243d0     /                                     
       DATA ZZZ /-0.738900483d0     /  
      SHAPE1=0.d0                                   
      AAA=ACOM 
cC THIS IS SPECIAL FOR MISHA'S BAIKAL STUFF
c      BBB=-1
C     XXX=0.                                                            
C     WRITE(19,*) X,Y,ELEV1,SLOPE1                                      
C ... CALCULATE LATITUDE                                                    
C     CALL SETRIG 
      if(IDOT.EQ.25) then
c ..... with the lookup mass balance, there is no "pole"                                                      
        CALL RECPOL(X,Y,RLAT,RLONG)  
c        dist=sqrt((X-XPOLE)**2+(Y-YPOLE)**2)/5d3
c       tnslredu=tnsl*exp(-0.5d0*dist**2)
        tnslredu=tnsl
c        print *,real(X-XPOLE),real(Y-YPOLE),real(dist),tnslredu
        call lookup25(AAA,rlat,rlong,elev1,
     &              tnslredu,TS,acc,abl)
        ACCUM=acc
        if(accum.ne.-999.) return
      endif
      if(IDOT.EQ.26) then
c ..... with the lookup mass balance, there is no "pole"                                                      
        CALL RECPOL(X,Y,RLAT,RLONG)  
c        dist=sqrt((X-XPOLE)**2+(Y-YPOLE)**2)/2d3
c        tnslredu=tnsl*exp(-0.5d0*dist**2)
        tnslredu=tnsl
c        print *,real(X-XPOLE),real(Y-YPOLE),real(dist),tnslredu
        call lookup26(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS,acc,abl)
        ACCUM=acc
        if(accum.ne.-999.) return
      endif
      if(IDOT.EQ.27) then
c ..... with the lookup mass balance, there is no "pole"                                                      
        CALL RECPOL(X,Y,RLAT,RLONG)  
c        dist=sqrt((X-XPOLE)**2+(Y-YPOLE)**2)/2d3
c        tnslredu=tnsl*exp(-0.5d0*dist**2)
        tnslredu=tnsl
c        print *,real(X-XPOLE),real(Y-YPOLE),real(dist),tnslredu
        call lookup27(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS,acc,abl)
        ACCUM=acc
        if(accum.ne.-999.) return
      endif
      if(IDOT.EQ.28) then
c ..... with the lookup mass balance, there is no "pole"
        CALL RECPOL(X,Y,RLAT,RLONG)  
c        dist=sqrt((X-XPOLE)**2+(Y-YPOLE)**2)/2d3
c        tnslredu=tnsl*exp(-0.5d0*dist**2)
        tnslredu=tnsl
c        print *,real(X-XPOLE),real(Y-YPOLE),real(dist),tnslredu
        call lookup26(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS26,acc26,abl26)
        call lookup27(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS27,acc27,abl27)
        slfact=max(0.d0,min(1.d0,-SEALEV/50.d0))
c        print *,slfact,SEALEV
        TS=(1.d0-slfact)*TS26+slfact*TS27
        abl=(1.d0-slfact)*abl26+slfact*abl27
        ACCUM=(1.d0-slfact)*acc26+slfact*acc27
        if(accum.ne.-999.) return
      endif
      if(IDOT.EQ.29) then
c ..... with the lookup mass balance, there is no "pole"                                                      
        CALL RECPOL(X,Y,RLAT,RLONG)  
c        dist=sqrt((X-XPOLE)**2+(Y-YPOLE)**2)/2d3
c        tnslredu=tnsl*exp(-0.5d0*dist**2)
        tnslredu=tnsl
c        print *,real(X-XPOLE),real(Y-YPOLE),real(dist),tnslredu
        call lookup29(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS,acc,abl)
        ACCUM=acc
        if(accum.ne.-999.) return
      endif
      if(IDOT.EQ.30) then
c ..... with the lookup mass balance, there is no "pole"                                                      
        CALL RECPOL(X,Y,RLAT,RLONG)  
c        dist=sqrt((X-XPOLE)**2+(Y-YPOLE)**2)/2d3
c        tnslredu=tnsl*exp(-0.5d0*dist**2)
        tnslredu=tnsl
c        print *,real(X-XPOLE),real(Y-YPOLE),real(dist),tnslredu
        call lookup30(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS,acc,abl)
        ACCUM=acc
        if(accum.ne.-999.) return
      endif
      if(IDOT.EQ.31) then
c ..... with the lookup mass balance, there is no "pole"
        CALL RECPOL(X,Y,RLAT,RLONG)  
c        dist=sqrt((X-XPOLE)**2+(Y-YPOLE)**2)/2d3
c        tnslredu=tnsl*exp(-0.5d0*dist**2)
        tnslredu=tnsl
c        print *,real(X-XPOLE),real(Y-YPOLE),real(dist),tnslredu
        call lookup29(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS29,acc29,abl29)
        call lookup30(AAA,rlat,rlong,elev1,psurf1,
     &              tnslredu,TS30,acc30,abl30)
        slfact=max(0.d0,min(1.d0,-SEALEV/50.d0))
c        print *,slfact,SEALEV
        TS=(1.d0-slfact)*TS29+slfact*TS30
        abl=(1.d0-slfact)*abl29+slfact*abl30
        ACCUM=(1.d0-slfact)*acc29+slfact*acc30
        if(accum.ne.-999.) return
      endif
c ... but the fastook-prentice scheme can have a "pole"
      CALL RECPOL(X-XPOLE,Y-YPOLE,RLAT,RLONG)  
c ..........................................................
C      PRINT *,'X,Y,ELEV,SLOPE',X,Y,ELEV1,SLOPE1                                     
C      PRINT *,'LATITUDE=',RLAT                                          
C ... ELEVATION (KM)                                                        
      ELEV=ELEV1/1000.d0
      PSURF=PSURF1/1000.d0
C ... SLOPE (M/KM)                                                          
      SLOPE=SLOPE1*1000.d0
C ... SHAPE (M/KM/KM)                                                       
      SHAPE=SHAPE1*1000.d0*1000.d0
      SHAPE=0.d0
c ***********************************************************                               
c ***********************************************************                               
c *******EISMINT EXPERIMENT FOR GREENLAND *******************
C HUYBRECT'S METHOD FROM A LOOKUP TABLE
C ... PUT .TRUE. FOR EISMINT EXPERIMENT ...
      IF(.FALSE.) THEN
        IF(.TRUE.) THEN
C ....... GREENLAND
          CALL HUYB(ELEV1,RLAT,PDD,TS,TSUMMER,TNSL)
        ELSEIF(.FALSE.) THEN
C ....... ANTARCTICA
          CALL HUYBANT(ELEV1,RLAT,PDD,TS,TSUMMER,TNSL)
          ACC=1.5D2*2.D0**(TS/10.D0)
        ELSEIF(.FALSE.) THEN
C ....... SLIDING EXPERIMENT
          CALL PAYNE(X,Y,TS,ACC)
          PDD=0.0d0
        ENDIF
C ... FOR ICE
        ABL=0.8d0*PDD
C ... FOR SNOW
C      ABL=0.3d0*PDD
C ... WHAT I HAVE BEEN USING
C      ABL=0.6d0*PDD
        ACC=0.0d0
        ACCNET=ACC-ABL                                                    
        ACCUM=ACCNET*.01d0
C        PRINT *,ACC,-ABL,ACC-ABL                                                
        RETURN
      ENDIF
C ***********************************************************                               
C ***********************************************************                               
c ***********************************************************                               
C ... CALCULATE SURFACE MEAN ANNUAL AIR TEMP      

c ... EXPERIMENTAL !!! >>>>>                          
c        tnslredu=0.5*(1.+RLAT/90.)
c        TS=AAA*ELEV+BBB*RLAT+CCC+(TNSL*tnslredu)+TNSLBASE            
c ... EXPERIMENTAL !!! >>>>>                          

      TS=AAA*ELEV+BBB*RLAT+CCC+TNSL+TNSLBASE            
c**************MARS***************************************
      IF(.true.) THEN
        IF(TFLAG) TS=TNSL+TNSLBASE-100
        IF(TFLAG) TS=AAA*ELEV+TNSL+TNSLBASE-100
        ACCUM=accmars(ts,elev*1000,acc,abl)
c       ACCUM=accmars1(ts,elev*1000,acc,abl)
        return
      ENDIF
c**************MARS***************************************
c**************EISMINT***************************************
      IF(.false.) THEN
        Z=MAX(ELEV1,20.d0*(RLAT-65.d0))
        TS=49.13d0-0.007992d0*Z-0.7576d0*RLAT+TNSL
      ENDIF
c**************EISMINT***************************************
C     PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TS                        
C ... CALCULATE MEAN ANNUAL TEMP OF FREE ATMOSPHERE-ISOTHEMAL LAYER         
      TF=0.67d0*(TS+273.0d0)+88.9d0
C     PRINT *,'TEMP FREE AT-ISOTHRMAL LAYER=',TF                        
C ... CALCULATE SATURATION VAPOR PRESSURE                                   
      TERM1=-9.09718d0*(273.16d0/TF-1.0d0)                                    
      TERM2=-3.56654d0*LOG10(273.16d0/TF)                                   
      TERM3=0.876793d0*(1.0d0-TF/273.16d0)+0.785835d0
      EXPON=TERM1+TERM2+TERM3                                           
      ES=10.d0**EXPON                                                     
C     PRINT *,'SATURATION VAPOR PRESSURE=',ES                           
C ... CALCULATE ACCUMULATION RATE (M/YR)                                    
      TERM1=WWW*ES                                                      
      TERM2=XXX*SLOPE 

c **** <<< EXPERIMENTAL >>> **** turn off slope term ...
      TERM2=0
c **** <<< EXPERIMENTAL >>> ****

      TERM3=ZZZ                                                         
      TERM4=-15.276d0*SHAPE                                               
      ACC=max(0.d0,TERM1+TERM2+TERM3+TERM4)
C     WRITE(19,113) TERM1,TERM2,TERM3,ACC 
113   FORMAT(4G13.6)                                                    
C     PRINT *,'ACCUMULATION=',ACC                                       
C ... CALCULATE ABLATION                                                    
      QY=0.d0
      DO I=1,12                                                      
        QY=QY+QI(I)-QS(I)*RLAT                                          
      ENDDO                                                          
      QY=QY/12.d0
      PDD=0.d0
      DO I=1,12                                                      
C       TTT(I)=TS+0.021d0*((QI(I)+QS(I)*RLAT)-QY)+8.954d0
        TTT(I)=TS+0.021d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C       TTT(I)=TS+0.018d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C       WRITE(19,*) I,TTT(I),QI(I)-QS(I)*RLAT,QY                        
        IF(TTT(I).GT.0.0) PDD=PDD+30.d0*TTT(I)                            
      ENDDO
c      print *,ts,pdd
      ABL=.6d0*PDD                                                        
C ... CALCULATE NET ACCUMULATION                                            
      ACCNET=ACC-ABL                                                    
C     ACCNET=ACC                                                        
C     PRINT *,'NET ACCUMULATION/ABLATION=',ACCNET                       
      ACCUM=ACCNET*.01d0                                                  
C     IF(PDD.GT.0.) THEN                                                
C       WRITE(19,111) ELEV,SLOPE,RLAT,TS,TF-273.,ES,PDD,ACC,ABL,ACCNET  
C     ENDIF                                                             
111   FORMAT(6F7.2,F6.0,3F7.1)                                          
C     WRITE(19,112) (TTT(I),I=1,12)                                     
112   FORMAT(12F6.1)        
      END                                                               
c=====================================================================
      FUNCTION ACCUMH(X,Y,ELEV1,TNSL,TMA)                      
      IMPLICIT REAL*8(A-H,O-Z)
      external func
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      COMMON /LOCAL/ SIGMA,TMAA,TMS,PI2
      DATA A1 /34.46D0/, A2 /-0.00914D0/, A3 /-0.68775D0/
      DATA B1 /16.81D0/, B2 /-0.00692D0/, B3 /-0.27937D0/
      DATA C1 /0.78D0/, C2 /2.525D-2/, C3 /2.225D-4/
      sigma=5.d0
      PI2=8.d0*ATAN(1.d0)
      ACOM=A2
c      PRINT *,'X,Y,ELEV=',X,Y,ELEV1                                      
C ... CALCULATE LATITUDE                                                    
C     CALL SETRIG                                                       
      CALL RECPOL(X,Y,RLAT,RLONG)                                       
c      PRINT *,'LATITUDE=',RLAT                                          
C ... CALCULATE SURFACE MEAN ANNUAL AIR TEMP                                
      TMA=A1+A2*ELEV1+A3*RLAT+TNSL+TNSLBASE
      TMAA=TMA                               
c      PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TMA                        
C ... CALCULATE SUMMER ANNUAL TEMP
      tms=B1+B2*ELEV1+B3*RLAT+TNSL+TNSLBASE                                
c      PRINT *,'TEMP SUMMER=',tms  
C ... CALCULATE ACCUMULATION RATE (M/YR)                                    
      ACC=C1+C2*TMA+C3*TMA**2                                       
c      PRINT *,'ACCUMULATION=',ACC                                       
C ... CALCULATE ABLATION 
      CALL QUAD2D(0.D0,365.D0,PDD)
      PDD=PDD/SIGMA/SQRT(PI2)
c      PRINT *,'PDD=',PDD
      ABL=.008d0*PDD                                                        
c      PRINT *,'ABLATION=',ABL                                           
C ... CALCULATE NET ACCUMULATION                                            
      ACCNET=ACC-ABL                                                    
c      PRINT *,'NET ACCUMULATION/ABLATION=',ACCNET                       
      ACCUMH=ACCNET                                                  
      END                                                               
C
      SUBROUTINE QUAD2D(X1,X2,SS)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL H
      CALL QGAUSX(H,X1,X2,SS)
c      PRINT *,'SS IN QUAD2D',SS
      RETURN
      END
C     
      FUNCTION F(YY)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FUNC
      COMMON /XY/ X,Y
      Y=YY
      F=FUNC(X,Y)
      RETURN
      END
C     
      FUNCTION H(XX)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL Y1,Y2,F
      COMMON /XY/ X,Y
      X=XX
      CALL QGAUSY(F,Y1(X),Y2(X),SS)
      H=SS
      RETURN
      END
C
      FUNCTION FUNC(X,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL TD
      COMMON /LOCAL/ SIGMA,TMA,TMS,PI2
      temp=Y*EXP(-0.5d0*(Y-TD(X))**2/SIGMA)
      func=temp
      END
C
      FUNCTION Y1(X)
      IMPLICIT REAL*8(A-H,O-Z)
      Y1=0.d0*X
      END
C
      FUNCTION Y2(X)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL TD
      COMMON /LOCAL/ SIGMA,TMA,TMS,PI2
      Y2=TD(X)+2.5d0*SIGMA
c      Y2=TD(X)
      END
C
      FUNCTION TD(T)
      IMPLICIT REAL*8(A-H,O-Z)
      INTRINSIC COS
      COMMON /LOCAL/ SIGMA,TMA,TMS,PI2
      TD=TMA+(TMA-TMS)*COS(PI2*T/365.d0)
      END
C
      SUBROUTINE QGAUSX(F,A,B,SS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(ITMAX=20)
      DIMENSION X(ITMAX),W(ITMAX)
      EXTERNAL F
      DATA TOL /1D-2/
C      SUMOLD=1d30
      N=ITMAX
C      DO N=ITMAX,ITMAX
        CALL GAULEG(A,B,X,W,N)
        SUM=0.d0
        DO I=1,N
          SUM=SUM+W(I)*F(X(I))
        ENDDO
C        IF(ABS(SUM-SUMOLD).LT.ABS(SUM*TOL) .or. 
C     &          sum-sumold.eq.0.) THEN
C          SS=SUM
C           PRINT *,'qgausx',N,ss
C          RETURN
C        ENDIF
C        PRINT *,'qgausx',I,SUM,SUMOLD
C        SUMOLD=SUM
C      ENDDO
      ss=sum
C      PRINT *,'QGAUSX DIDNT CONVERGE IN ',ITMAX
      END
C
      SUBROUTINE QGAUSY(F,A,B,SS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(ITMAX=8)
      DIMENSION X(ITMAX),W(ITMAX)
      EXTERNAL F
      DATA TOL /1D-2/
C      SUMOLD=1d30
      N=ITMAX
C      DO N=ITMAX,ITMAX
        CALL GAULEG(A,B,X,W,N)
        SUM=0.d0
        DO I=1,N
          SUM=SUM+W(I)*F(X(I))
        ENDDO
C        IF(ABS(SUM-SUMOLD).LT.ABS(SUM*TOL) .or. 
C     &          sum-sumold.eq.0.) THEN
C          SS=SUM
c          PRINT *,'qgausy',N,ss
C          RETURN
C        ENDIF
C        SUMOLD=SUM
C      ENDDO
      ss=sum
C      PRINT *,'QGAUSY DIDNT CONVERGE IN ',ITMAX
      END
C
      SUBROUTINE GAULEG(X1,X2,X,W,N)
C
C ... GIVEN THE LOWER AND UPPER LIMITS OF INTEGRATION X1 AND X2, AND
C ... GIVEN THE NUMBER OF NODAL POINTS DESIRED N, THIS ROUTINE RETURNS
C ... ARRAYS X AND W OF LENGTH N, CONTAINING THE NODES AND WEIGHTS OF
C ... THE GUASS-LEGENDRE N-POINT QUADRATURE FORMULA
C
C ... HIGH PRECISION IS A GOOD IDEA FOR THIS ROUTINE
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X1,X2,X(N),W(N)
C ... INCREASE EPS IF YOU DON'T HAVE THIS PRECISION
      PARAMETER (EPS=3.D-14)
C ... THE ROOTS ARE SYMMETRIC IN THE INTERVAL, SO WE ONLY HAVE TO
C ... FIND HALF OF THEM
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
C ... LOOP OVER THE DESIRED ROOTS
      DO 12 I=1,M
        Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))
C ..... STARTING WITH THE ABOVE APPROXIMATION TO THE ITH ROOT,
C ..... WE ENTER THE MAIN LOOP OF REFINEMENT BY THE NEWTON'S METHOD
1       CONTINUE
          P1=1.D0
          P2=0.D0
C ....... LOOP UP THE RECURRENCE RELATION TO GET THE LEGENDRE
C ....... POLYNOMIAL EVALUATED AT Z
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11        CONTINUE
C ....... P1 IS NOW THE DESIRED LEGENDRE POLYNOMIAL. WE NEXT COMPUTE PP,
C ....... ITS DIRIVATIVE, BY A STANDARD RELATION INVOLVING ALSO P2,
C ....... THE POLYNIMIAL OF ONE LOWER ORDER
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
C ....... NEWTON'S METHOD
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
C ..... SCALE THE ROOT TO ITS DESIRED INTERVAL
        X(I)=XM-XL*Z
C ..... AND PUT IN ITS SYMMETRIC COUNTERPART
        X(N+1-I)=XM+XL*Z
C ..... COMPUTE THE WEIGHT
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
C ..... AND ITS SYMMETRIC COUNTERPART
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END
c================================================
      subroutine huyb(elev1,rlat,pdd,TS,TSUMMER,TNSL)
      IMPLICIT REAL*8(A-H,O-Z)
c following covers greenland
      PARAMETER(JMAX=8,RSTEP=5.D0,RMIN=55.D0,
     &          RMAX=RMIN+(JMAX-1)*RSTEP)
c following for full n. hem
c      PARAMETER(JMAX=19,RSTEP=5.D0,RMIN=0.D0,
c     &          RMAX=RMIN+(JMAX-1)*RSTEP)
c
      PARAMETER(IMAX=26,ESTEP=200.D0,EMIN=0.D0,
     &          EMAX=EMIN+(IMAX-1)*ESTEP)
c high resolution ...
c      PARAMETER(JMAX=31,RSTEP=1.D0,RMIN=55.D0,
c     &          RMAX=RMIN+(JMAX-1)*RSTEP)
c      PARAMETER(IMAX=101,ESTEP=50.D0,EMIN=0.D0,
c     &          EMAX=EMIN+(IMAX-1)*ESTEP)
      DIMENSION RMATRIX(IMAX,JMAX)
      COMMON /LOCAL/ SIGMA,TMA,TMS,PI2
      SAVE IPASS,RMATRIX,TSAVE
      DATA IPASS /0/, TSAVE /-999./
      IF(IPASS.EQ.0 .OR. TNSL.NE.TSAVE) THEN
        IPASS=1
        TSAVE=TNSL
        sigma=5.d0
        PI2=8.d0*ATAN(1.d0)
        print *,'preparing ',real(rmin),real(rmax),real(rstep)
        print *,'          ',real(emin),real(emax),real(estep)
        DO J=1,JMAX
          RL=RMIN+(J-1)*RSTEP
          DO I=1,IMAX
            EL=EMIN+(I-1)*ESTEP
            Z=MAX(EL,20.d0*(RL-65.d0))
c*********************************************************************
c the offset (3.0) in the following is a fudge to get the present config to be
c stable. the temp TS passed out for the surface temperature does not
c include it....
            OFFSET=3.0d0
            TMA=49.13d0-0.007992d0*Z-0.7576d0*RL+TNSL+OFFSET
            TMS=30.78d0-0.006277d0*EL-0.3262d0*RL+TNSL+OFFSET 
C            TMA=49.13d0-0.007992d0*Z-0.7576d0*RL+TNSL
C            TMS=30.78d0-0.006277d0*EL-0.3262d0*RL+TNSL   
c*********************************************************************
C ... THIS IS HUYBRECTS, VERY SLOW...
            CALL QUAD2D(0.D0,365.D0,PDDH)
            PDDH=PDDH/SIGMA/SQRT(PI2)
c           print *,real(tma),real(tms),real(pddh),real(pddh*0.8*.01)
            RMATRIX(I,J)=PDDH
          ENDDO
        ENDDO
      ENDIF
      IF(ELEV1.GE.EMIN .AND. ELEV1.LE.EMAX .AND.
     &   RLAT.GE.RMIN  .AND. RLAT.LE.RMAX) THEN
        I=INT((ELEV1-EMIN)/ESTEP)+1
        I=MIN(I,IMAX)
        I=MAX(I,1)
        J=INT((RLAT-RMIN)/RSTEP)+1
        J=MIN(J,JMAX)
        J=MAX(J,1)
        ELEVI=-999.d0
        RLATI=-999.d0
        IF(I.EQ.IMAX) THEN
          I=I-1
          X=1.D0
        ELSE
          ELEVI=EMIN+(I-1)*ESTEP    
          X=(ELEV1-ELEVI)/ESTEP
        ENDIF
        IF(J.EQ.JMAX) THEN
          J=J-1
          Y=1.D0
        ELSE
          RLATI=RMIN+(J-1)*RSTEP    
          Y=(RLAT-RLATI)/RSTEP
        ENDIF
c       PRINT *,'I,J              ',I,J
c       PRINT *,'ELEV1,RLAT       ',ELEV1,RLAT
c       PRINT *,'ELEVI,RLATI      ',ELEVI,RLATI
c       PRINT *,'X,Y              ',X,Y
c       PRINT *,'(I,J),    (I+1,J)',RMATRIX(I,J),RMATRIX(I+1,J)
c       PRINT *,'(I+1,J+1),(I,J+1)',RMATRIX(I+1,J+1),RMATRIX(I,J+1)
        P1=(1.D0-X)*(1.D0-Y)
        P2=X*(1.D0-Y)
        P3=X*Y
        P4=(1.D0-X)*Y
        PDD=P1*RMATRIX(I,J)+P2*RMATRIX(I+1,J)+
     &      P3*RMATRIX(I+1,J+1)+P4*RMATRIX(I,J+1)
        Z=MAX(ELEV1,20.d0*(RLAT-65.d0))
c*********************************************************************
c these temps are passed out of the rouytine, and dont include the 
c OFFSET degree fudge included in the equations above
        TS=49.13d0-0.007992d0*Z-0.7576d0*RLAT+TNSL
        TSUMMER=30.78d0-0.006277d0*ELEV1-0.3262d0*rlat+TNSL
c*********************************************************************
      ELSE
        PRINT *,'PROBLEMS, OUT OF RANGE',elev1,rlat
        PDD=1000.d0
        Z=MAX(ELEV1,20.d0*(RLAT-65.d0))
        TS=49.13d0-0.007992d0*Z-0.7576d0*RLAT+TNSL
        TSUMMER=30.78d0-0.006277d0*ELEV1-0.3262d0*rlat+TNSL
C        PAUSE
      ENDIF
      END

c================================================
      subroutine huybant(elev1,rlat,pdd,TS,TSUMMER,TNSL)
      IMPLICIT REAL*8(A-H,O-Z)
c following covers greenland
      PARAMETER(JMAX=9,RSTEP=5.D0,RMIN=50.D0,
     &          RMAX=RMIN+(JMAX-1)*RSTEP)
c following for full n. hem
c      PARAMETER(JMAX=19,RSTEP=5.D0,RMIN=0.D0,
c     &          RMAX=RMIN+(JMAX-1)*RSTEP)
      PARAMETER(IMAX=26,ESTEP=200.D0,EMIN=0.D0,
     &          EMAX=EMIN+(IMAX-1)*ESTEP)
c      PARAMETER(JMAX=31,RSTEP=1.D0,RMIN=55.D0,
c     &          RMAX=RMIN+(JMAX-1)*RSTEP)
c      PARAMETER(IMAX=101,ESTEP=50.D0,EMIN=0.D0,
c     &          EMAX=EMIN+(IMAX-1)*ESTEP)
      DIMENSION RMATRIX(IMAX,JMAX)
      COMMON /LOCAL/ SIGMA,TMA,TMS,PI2
      SAVE IPASS,RMATRIX,TSAVE
      DATA IPASS /0/, TSAVE /-999./
      IF(IPASS.EQ.0 .OR. TNSL.NE.TSAVE) THEN
        IPASS=1
        TSAVE=TNSL
        sigma=5.d0
        PI2=8.d0*ATAN(1.d0)
        print *,'preparing ',real(rmin),real(rmax),real(rstep)
        print *,'          ',real(emin),real(emax),real(estep)
        DO J=1,JMAX
          RL=RMIN+(J-1)*RSTEP
          DO I=1,IMAX
            EL=EMIN+(I-1)*ESTEP
            Z=MAX(EL,20.d0*(RL-65.d0))
c*********************************************************************
            TMA=34.46d0-0.00914d0*EL-0.68775d0*RL+TNSL                 
            TMS=16.81d0-0.00692d0*EL-0.27973d0*RL+TNSL
c*********************************************************************
C ... THIS IS HUYBRECTS, VERY SLOW...
            CALL QUAD2D(0.D0,365.D0,PDDH)
            PDDH=PDDH/SIGMA/SQRT(PI2)
            RMATRIX(I,J)=PDDH
          ENDDO
        ENDDO
      ENDIF
      IF(ELEV1.GE.EMIN .AND. ELEV1.LE.EMAX .AND.
     &   RLAT.GE.RMIN  .AND. RLAT.LE.RMAX) THEN
        I=INT((ELEV1-EMIN)/ESTEP)+1
        I=MIN(I,IMAX)
        I=MAX(I,1)
        J=INT((RLAT-RMIN)/RSTEP)+1
        J=MIN(J,JMAX)
        J=MAX(J,1)
        ELEVI=-999.d0
        RLATI=-999.d0
        IF(I.EQ.IMAX) THEN
          I=I-1
          X=1.D0
        ELSE
          ELEVI=EMIN+(I-1)*ESTEP    
          X=(ELEV1-ELEVI)/ESTEP
        ENDIF
        IF(J.EQ.JMAX) THEN
          J=J-1
          Y=1.D0
        ELSE
          RLATI=RMIN+(J-1)*RSTEP    
          Y=(RLAT-RLATI)/RSTEP
        ENDIF
c       PRINT *,'I,J              ',I,J
c       PRINT *,'ELEV1,RLAT       ',ELEV1,RLAT
c       PRINT *,'ELEVI,RLATI      ',ELEVI,RLATI
c       PRINT *,'X,Y              ',X,Y
c       PRINT *,'(I,J),    (I+1,J)',RMATRIX(I,J),RMATRIX(I+1,J)
c       PRINT *,'(I+1,J+1),(I,J+1)',RMATRIX(I+1,J+1),RMATRIX(I,J+1)
        P1=(1.D0-X)*(1.D0-Y)
        P2=X*(1.D0-Y)
        P3=X*Y
        P4=(1.D0-X)*Y
        PDD=P1*RMATRIX(I,J)+P2*RMATRIX(I+1,J)+
     &      P3*RMATRIX(I+1,J+1)+P4*RMATRIX(I,J+1)
        Z=MAX(ELEV1,20.d0*(RLAT-65.d0))
c*********************************************************************
c these temps are passed out of the routine
        TS=34.46d0-0.00914d0*ELEV1-0.68775d0*rlat+TNSL                 
        TSUMMER=16.81d0-0.00692d0*ELEV1-0.27973d0*rlat+TNSL
c*********************************************************************
      ELSE
        PRINT *,'PROBLEMS, OUT OF RANGE',elev1,rlat
        PAUSE
      ENDIF
      END

C===============================================
      FUNCTION WARMING(DT)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      IF(DT.GT.0.0) THEN
        S=0.05D0
      ELSEIF(DT.GT.-10.0) THEN
        S=0.05D0-0.005D0*DT
      ELSE
        S=0.10D0
      ENDIF
      tmp=(1.D0+S)**DT
      WARMING=tmp
      END
C===============================================
      subroutine payne(X,Y,TS,ACC)
      IMPLICIT REAL*8(A-H,O-Z)                                          
c      data bmax /0.5d0/
c      data sb /.01d0/
c      data e /450.d0/
c      data tmin /248.15d0/
c      data st /0.01d0/
c      data xbar,ybar /850.d0,850.d0/
      save ipass,bmax,sb,e,tmin,st,xbar,ybar
      data ipass /0/
      if(ipass.eq.0) then
        ipass=1
        print *,' checking for payne data file... '
        read(89,*,end=999) bmax,sb,e,tmin,st,xbar,ybar
        print *,' found:',bmax,sb,e,tmin,st,xbar,ybar
        goto 1000
999       continue
          bmax=0.5d0
          sb=0.01d0
          e=450.d0
          tmin=248.15d0
          st=0.01d0
          xbar=750.d0
          ybar=750.d0
          print *,' not found:default values'
          print *,bmax,sb,e,tmin,st,xbar,ybar
          rewind 89
          write(89,*) bmax,sb,e,tmin,st,xbar,ybar
1000    continue
      endif
      acc=100.d0*min(bmax,sb*(e-sqrt((x-xbar)**2+(y-ybar)**2)))
      ts=tmin+st*(sqrt((x-xbar)**2+(y-ybar)**2))
      ts=ts-273.15d0
      end
C===============================================
      SUBROUTINE PAYNEOUT1(TIME,NUMNP,NUMEL,KX,X,Y,HTICE,DEPB,
     &                    TEMPA,VEL,FLOWA,
     &                    VOL,AREA,RMAF)
      IMPLICIT REAL*8(A-H,O-Z)   
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM,MMAX=40)
      common /timeout/ itimeout
      dimension x(nmax),y(nmax),depb(nmax),vel(nmax,3)
      dimension flowa(nmax),tempa(mmax,nmax),htice(nmax)
      dimension KX(NMAX,4)
c      return
      itimeout=30000
      itimeout=200000
      print *,'EISMINT DATA DUMP AT',itimeout
      tt=time*.001d0
c ... ume_q_             
      if(nint(time).eq.itimeout) then
c ... plan-form at 200kyr ...			.p.
c ..... x, y, (km) data						x,y
        do i=1,numnp
          xx=x(i)*0.001d0
          yy=y(i)*0.001d0
          thick=htice(i)-depb(i)
c ....... 50 thickness		(m)		.ume_q_p_tk.	htice-depb 
          write(50,1) xx,yy,thick
c ....... 51 basal temperature	(K)		.ume_q_p_tp.	tempa(mmax,i)
          write(51,1) xx,yy,tempa(mmax,i)
        enddo
        do i=1,numel
          xx=0.d0
          yy=0.d0
          thick=0.d0
          do j=1,4
            xx=xx+x(kx(i,j))*0.001d0
            yy=yy+y(kx(i,j))*0.001d0
            thick=thick+htice(kx(i,j))-depb(kx(i,j))
          enddo
          xx=0.25d0*xx
          yy=0.25d0*yy
          thick=0.25d0*thick
c ....... 65 velocity magnitude (m/yr)		.ume_q_p_vmag.	vel(,1)
          write(65,1) xx,yy,vel(i,1)
c ....... 66 horizontal flux(m**2/yr)		.ume_q_p_uvq.	vel(,1)*thick
          write(66,1) xx,yy,vel(i,1)*thick
        enddo
c
      endif
c ... global time series ...			.t.
c ..... time (kyr) data
c ..... 55 volume		(km**3)		.ume_q_t_vo.	vol
        write(55,1) tt,vol
c ..... 56 area 		(km**2)		.ume_q_t_ar.	area
        write(56,1) tt,area
c ..... 57 melt area fraction			.ume_q_t_fr.	rmaf
        write(57,1) tt,rmaf
c
c ... local time series ...			.t.   at stations 1,2,3,4,5
c ..... time (kyr) data
c ..... 58 thickness 		(m)		.ume_q_t_tk_1,2,3,4,5.
        write(58,1) tt,htice(1861)-depb(1861)
     &                ,htice(2471)-depb(2471)
     &                ,htice(2776)-depb(2776)
     &                ,htice(2109)-depb(2109)
     &                ,htice(2719)-depb(2719)
c ..... 59 basal temperature 	(K) 		.ume_q_t_tp_1,2,3,4,5.
        write(59,1) tt,tempa(mmax,1861)
     &                ,tempa(mmax,2471)
     &                ,tempa(mmax,2776)
     &                ,tempa(mmax,2109)
     &                ,tempa(mmax,2719)
1     format(1x,5(2x,1pg13.6))
      end
c==========================================================
      subroutine payneout2(time,ipt,mmax,tttt,xxx,at,wwww)
      IMPLICIT REAL*8(A-H,O-Z) 
      common /timeout/ itimeout
      dimension tttt(mmax),xxx(mmax),at(mmax),wwww(mmax) 
c      return
      if(nint(time).ne.itimeout) return       
c
c ... depth profiles ...			.d.   at stations 1,2,3,4,5
c ... height above bedrock (m), data
      write(60,*) mmax
      write(61,*) mmax
      write(62,*) mmax
      write(63,*) mmax
      write(64,*) mmax
      do i=mmax,1,-1
c ..... 60 temperature 		(K)		.ume_q_d_tp_1,2,3,4,5.
        write(60,1) ipt,xxx(i),tttt(i)
c ..... 61 velocity-x 		(m/yr)		.ume_q_d_uv_1,2,3,4,5.
        write(61,1) ipt,xxx(i),-99999.
c ..... 62 velocity-y 		(m/yr)		.ume_q_d_vv_1,2,3,4,5.
        write(62,1) ipt,xxx(i),-99999.
c ..... 63 velocity-z 		(m/yr)		.ume_q_d_wv_1,2,3,4,5.
        write(63,1) ipt,xxx(i),wwww(i)
c ..... 64 flow coefficient 	(Pa**-2 s**-1)	.ume_q_d_af_1,2,3,4,5.
        atiii=3.168876462d-23/(at(i)**3)
        write(64,1) ipt,xxx(i),atiii
      enddo
1     format(1x,i10,5(1x,1pg13.6))
      end
c==========================================================
      FUNCTION WOABLATE(X,Y,ELEV1,SLOPE1,SHAPE1,TNSL,TS)                      
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION QI(12),QS(12),TTT(12)                                   
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      COMMON /LOCAL/ SIGMA,TMAA,TMS,PI2
      DATA QI /960.d0,1036.d0,1200.d0,825.,330.d0,90.d0,150.d0,
     &         600.d0,1200.d0,1020.d0,    
     &         930.d0,850.d0/                                               
      DATA QS /-0.667d0,4.6d0,11.667d0,9.167d0,3.667d0,
     *         1.d0,1.667d0,6.667d0,12.d0,6.333d0,  
     &         0.333d0,-3.333d0/                                            
C     DATA AAA /-9.14/, BBB /-.68/, CCC /34.461/                        
C     DATA WWW /13.05/, XXX /.664/, ZZZ /2.608/                         
       DATA AAA / -9.62376690d0     /                                     
       DATA BBB /-0.546917617d0     /                                     
       DATA CCC /  24.9793854d0     /                                     
       DATA WWW /  19.1390686d0     /                                     
       DATA XXX / 0.922791243d0     /                                     
       DATA ZZZ /-0.738900483d0     /                                     
      AAA=ACOM 
C ... CALCULATE LATITUDE                                                    
C     CALL SETRIG                                                       
      CALL RECPOL(X-XPOLE,Y-YPOLE,RLAT,RLONG)  
C ... ELEVATION (KM)                                                        
      ELEV=ELEV1/1000.d0
C ... SLOPE (M/KM)                                                          
      SLOPE=SLOPE1*1000.d0
C ... SHAPE (M/KM/KM)                                                       
      SHAPE=SHAPE1*1000.d0*1000.d0
      SHAPE=0.d0
C ... CALCULATE SURFACE MEAN ANNUAL AIR TEMP                                
      TS=AAA*ELEV+BBB*RLAT+CCC+TNSL+TNSLBASE                                
C ... CALCULATE MEAN ANNUAL TEMP OF FREE ATMOSPHERE-ISOTHEMAL LAYER         
      TF=0.67d0*(TS+273.0d0)+88.9d0
C ... CALCULATE SATURATION VAPOR PRESSURE                                   
      TERM1=-9.09718d0*(273.16d0/TF-1.0d0)                                    
      TERM2=-3.56654d0*LOG10(273.16d0/TF)                                   
      TERM3=0.876793d0*(1.0d0-TF/273.16d0)+0.785835d0
      EXPON=TERM1+TERM2+TERM3                                           
      ES=10.d0**EXPON                                                     
C ... CALCULATE ACCUMULATION RATE (M/YR)                                    
      TERM1=WWW*ES                                                      
      TERM2=XXX*SLOPE                                                   

c **** <<< EXPERIMENTAL >>> **** turn off slope term ...
      TERM2=0
c **** <<< EXPERIMENTAL >>> ****

      TERM3=ZZZ                                                         
      TERM4=-15.276d0*SHAPE                                               
      ACC=max(0.d0,TERM1+TERM2+TERM3+TERM4)
C ... CALCULATE ABLATION                                                    
C ... CALCULATE ABLATION                                                    
      QY=0.d0
      DO I=1,12                                                      
        QY=QY+QI(I)-QS(I)*RLAT                                          
      ENDDO                                                          
      QY=QY/12.d0
      PDD=0.d0
      DO I=1,12                                                      
C       TTT(I)=TS+0.021d0*((QI(I)+QS(I)*RLAT)-QY)+8.954d0
        TTT(I)=TS+0.021d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C       TTT(I)=TS+0.018d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C       WRITE(19,*) I,TTT(I),QI(I)-QS(I)*RLAT,QY                        
        IF(TTT(I).GT.0.0) PDD=PDD+30.d0*TTT(I)                            
      ENDDO
      ABL=.6d0*PDD                                                        
c ... LINEAR DECREASE IN ABLATION FROM NOMINAL VALUE AT 200 M (0.2 KM)
C     DOWN TO ZERO AT SEA LEVEL (USING ELEV1, ELEV IN KM)
      FACT=MIN(1.0D0,ELEV/0.2D0)
      ABL=ABL*FACT
C ... CALCULATE NET ACCUMULATION                                            
      ACCNET=ACC-ABL                                                    
      WOABLATE=ACCNET*.01d0                                                  
      END                                                               
c---------------------------------------------------
      subroutine milankold(TIME,x,y,elev1,slope1,shape1,TNSL,tmean,
     &                     accum)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nyear1=-100,nyear2=100,nlat=19,nmonth=12)
      integer insol(nmonth,nlat,nyear1:nyear2)
      integer lat(nlat)
      dimension rmil(nmonth),rmil1(nmonth),rmil2(nmonth)
      character*80 line
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      DIMENSION QI(12),QS(12),TTT(12)                                   
      DATA QI /960.d0,1036.d0,1200.d0,825.,330.d0,90.d0,150.d0,
     &         600.d0,1200.d0,1020.d0,    
     &         930.d0,850.d0/                                               
      DATA QS /-0.667d0,4.6d0,11.667d0,9.167d0,3.667d0,
     *         1.d0,1.667d0,6.667d0,12.d0,6.333d0,  
     &         0.333d0,-3.333d0/                                            
      data ipass /0/, sigma /5.67e-5/,pi /3.1415927/
      data www /  19.1390686e0     /                                     
      data xxx / 0.922791243e0     /                                     
      data zzz /-0.738900483e0     /                                     
      save ipass,insol,iy,lat
      return
C ... CALCULATE LATITUDE                                                    
C     CALL SETRIG                                                       
c     CALL RECPOL(X-XPOLE,Y-YPOLE,RLAT,RLONG)  
      CALL RECPOL(X,Y,RLAT,RLONG)  
      if(ipass.eq.0) then
c ... done once first call ...
        ipass=1
c        open(23,file='milank.dat')
        do iyear=nyear2,nyear1,-1
          read(23,'(a)') line
          do ilat=1,nlat
            read(23,*) iy,lat(ilat),
     &                (insol(im,ilat,iyear),im=1,nmonth)
          enddo
        enddo
      endif
      ryear=TIME/1000.
      if(ryear.lt.nyear1) ryear=nyear1
      if(ryear.gt.nyear2) ryear=nyear2
c ..............................................
C ... CALCULATE LATITUDE (w/ pole shift .......                                                    
C     CALL SETRIG                                                       
      CALL RECPOL(X-XPOLE,Y-YPOLE,RLAT,RLONG)  
c     CALL RECPOL(X,Y,RLAT,RLONG)  
c ... elevation (km)                                                        
      elev=elev1/1000.d0
c ... slope (m/km)                                                          
      slope=slope1*1000.d0
c ... shape (m/km/km)                                                       
      shape=shape1*1000.d0*1000.d0
      shape=0.d0
      if(rlat.gt.90) rlat=90
      if(rlat.lt.-90) rlat=-90
      llat=int(rlat/10)*10
      ilat=1+(90-llat)/10
      rl=1+(90-rlat)/10
      ilat1=rl
      ilat2=rl+1
      if(ilat1.lt.1) ilat1=1
      if(ilat2.gt.nlat) ilat2=nlat
      iyear1=int(abs(ryear))
      iyear2=int(abs(ryear)+1)
      if(ryear.lt.0) then
        itemp=iyear1
        iyear1=-iyear2
        iyear2=-itemp
      endif
      if(iyear1.lt.nyear1) iyear1=nyear1
      if(iyear2.lt.nyear1) iyear2=nyear1
      if(iyear1.gt.nyear2) iyear1=nyear2
      if(iyear2.gt.nyear2) iyear2=nyear2
c ... interpolate latitude ...................
      fractl=(lat(ilat1)-rlat)/10
      do im=1,nmonth
        rmil1(im)=insol(im,ilat1,iyear1)+
     &    fractl*(insol(im,ilat2,iyear1)-insol(im,ilat1,iyear1))
        rmil2(im)=insol(im,ilat1,iyear2)+
     &    fractl*(insol(im,ilat2,iyear2)-insol(im,ilat1,iyear2))
      enddo
      if(.false.) then
        print *,lat(ilat1),rlat,lat(ilat2)
        print *,ilat1,rl,ilat2
        print *,fractl
        print 100,lat(ilat1),(insol(im,ilat1,iyear1),im=1,4)
        print 101,rlat,(rmil1(im),im=1,4)
        print 100,lat(ilat2),(insol(im,ilat2,iyear1),im=1,4)
        print *
        print 100,lat(ilat1),(insol(im,ilat1,iyear2),im=1,4)
        print 101,rlat,(rmil2(im),im=1,4)
        print 100,lat(ilat2),(insol(im,ilat2,iyear2),im=1,4)
100     format(1x,i6,12i6)
101     format(1x,f6.1,12f6.1)
      endif
c ............................................
c ... interpolate year ........................................
      fracty=(ryear-iyear1)
      do im=1,nmonth
        rmil(im)=rmil1(im)+fracty*(rmil2(im)-rmil1(im))
      enddo
      if(.false.) then
        print *,iyear1,ryear,iyear2
        print *,fracty
        print 200,iyear1,(rmil1(im),im=1,4)
        print 201,ryear,(rmil(im),im=1,4)
        print 200,iyear2,(rmil2(im),im=1,4)
200     format(1x,i6,12f6.1)
201     format(1x,f6.1,12f6.1)
c        pause
      endif
c ..............................................................
c ... means and standard deviation for year
      rmean=0.0
      do im=1,nmonth
        rmean=rmean+rmil(im)
      enddo
      rmean=rmean/nmonth
      rplus=0.0
      do im=1,nmonth
        rplus=rplus+(rmil(im)-rmean)**2
      enddo
      rplus=rplus/nmonth
      rplus=sqrt(rplus)
c ...............................................................
c ... convert langleys/day to watt/m**2, multiply by 0.4843 .....
      tmean=rmean*0.4843
      tplus=rplus*0.4843
c ... convert w/m**2 to erg/cm**2/s, multipy by 1000. ...........
      tmean=tmean*1000.
      tplus=tplus*1000.
c ... convert to temperature e=sigma*t**4, ......................
c ............  sigma=5.67e-5 erg/cm**2/deg**4/s ................
c ... tplus is temp with mean+stdev flux ........................
      tplus=((tmean+tplus/16)/sigma)**0.25
c ... tmean is temp with just mean ..............................
      tmean=(tmean/sigma)**0.25
c ... tplus is now amplitude of annual signal ...................
      tplus=tplus-tmean
c ... offset because stephan-boltzman is cold (greenhouse??) ....
      tmean=tmean+25.0-7.65
c ... convert to deg c ..........................................
      tmean=tmean-273.16
c ... lapse rate with elevation .................................
      tmean=tmean+acom*elev+TNSL+TNSLBASE
      TS=tmean
c ... rest is standard climatology ..............................
c     print *,'surface mean annual air temp=',tmean
c ... calculate mean annual temp of free atmosphere-isothemal layer ...         
      tf=0.67d0*(tmean+273.0d0)+88.9d0
c     print *,'temp free at-isothrmal layer=',tf
c ... calculate saturation vapor pressure .............................
      term1=-9.09718d0*(273.16d0/tf-1.0d0)
      term2=-3.56654d0*log10(273.16d0/tf)
      term3=0.876793d0*(1.0d0-tf/273.16d0)+0.785835d0
      expon=term1+term2+term3
      es=10.d0**expon
c     print *,'saturation vapor pressure=',es
c ... calculate accumulation rate (m/yr) ..............................
      term1=www*es
      term2=xxx*slope
      term3=zzz
      term4=-15.276d0*shape
      acc=max(0.d0,term1+term2+term3+term4)
C ... CALCULATE ABLATION THE OLD WAY ...
      if(.true.) then
        QY=0.d0
        DO I=1,12                                                      
          QY=QY+QI(I)-QS(I)*RLAT                                          
        ENDDO                                                          
        QY=QY/12.d0
        PDD=0.d0
        DO I=1,12                                                      
C         TTT(I)=TS+0.021d0*((QI(I)+QS(I)*RLAT)-QY)+8.954d0
          TTT(I)=TS+0.021d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C         TTT(I)=TS+0.018d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C         WRITE(19,*) I,TTT(I),QI(I)-QS(I)*RLAT,QY                        
          IF(TTT(I).GT.0.0) PDD=PDD+30.d0*TTT(I)                            
        ENDDO
      else
c ... sum up positive degree days for ablation calc..............
        pddnew=0
        do i=0,364
          temp=tmean+tplus*sin(2*pi*i/364.)
          if(temp.gt.0.) pddnew=pddnew+temp
        enddo
        pdd=pddnew
      endif
      abl=.6d0*pdd
c      print *,TS,pdd,pddnew
c ... calculate net accumulation ......................................
      accnet=acc-abl
c     print *,'net accumulation/ablation=',accnet
c ... put in m/yr .....................................................
      accum=accnet*.01d0
      abl=abl*.01d0
      end
c------------------------------------------------------------------
      subroutine milank(TIME,x,y,elev1,slope1,shape1,TNSL,tmean,accum)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nyear1=-100,nyear2=100,nlat=19,nmonth=12,amplify=0.2)
      parameter(rmult=4.,rparam=75)
      integer insol(nmonth,nlat,nyear1:nyear2)
      integer lat(nlat)
      dimension rmil(nmonth),rmil1(nmonth),rmil2(nmonth)
      dimension rmilba(nmonth)
      character*80 line
      logical dispon
      common /lapse/ acom,hmax,windir(2),xpole,ypole,abl,tnslbase
      DIMENSION QI(12),QS(12),TTT(12)                                   
      DATA QI /960.d0,1036.d0,1200.d0,825.,330.d0,90.d0,150.d0,
     &         600.d0,1200.d0,1020.d0,    
     &         930.d0,850.d0/                                               
      DATA QS /-0.667d0,4.6d0,11.667d0,9.167d0,3.667d0,
     *         1.d0,1.667d0,6.667d0,12.d0,6.333d0,  
     &         0.333d0,-3.333d0/                                            
      data ipass /0/,pi /3.1415927/, dispon /.false./
      data www /  19.1390686e0     /                                     
      data xxx / 0.922791243e0     /                                     
      data zzz /-0.738900483e0     /                                     
      save ipass,insol,iy,lat
      return
C ... CALCULATE LATITUDE                                                    
C     CALL SETRIG                                                       
c     CALL RECPOL(X-XPOLE,Y-YPOLE,RLAT,RLONG)  
      CALL RECPOL(X,Y,RLAT,RLONG)  
      if(ipass.eq.0) then
c ... done once first call ...
        ipass=1
c        open(23,file='milank.dat')
        do iyear=nyear2,nyear1,-1
          read(23,'(a)') line
          do ilat=1,nlat
            read(23,*) iy,lat(ilat),
     &                (insol(im,ilat,iyear),im=1,nmonth)
          enddo
        enddo
      endif
      ryear=TIME/1000.
      if(ryear.lt.nyear1) ryear=nyear1
      if(ryear.gt.nyear2) ryear=nyear2
c ..............................................
C ... CALCULATE LATITUDE (w/ pole shift .......                                                    
C     CALL SETRIG                                                       
      CALL RECPOL(X-XPOLE,Y-YPOLE,RLAT,RLONG)  
c     CALL RECPOL(X,Y,RLAT,RLONG)  
c ... elevation (km)                                                        
      elev=elev1/1000.d0
c ... slope (m/km)                                                          
      slope=slope1*1000.d0
c ... shape (m/km/km)                                                       
      shape=shape1*1000.d0*1000.d0
      shape=0.d0
      if(rlat.gt.90) rlat=90
      if(rlat.lt.-90) rlat=-90
      llat=int(rlat/10)*10
      ilat=1+(90-llat)/10
      rl=1+(90-rlat)/10
      ilat1=rl
      ilat2=rl+1
      if(.false.) then
        print *,lat(ilat1),rlat,lat(ilat2)
        print *,ilat1,rl,ilat2
      endif
      if(ilat1.lt.1) ilat1=1
      if(ilat2.gt.nlat) ilat2=nlat
      iyear1=int(abs(ryear))
      iyear2=int(abs(ryear)+1)
      if(ryear.lt.0) then
        itemp=iyear1
        iyear1=-iyear2
        iyear2=-itemp
      endif
      if(iyear1.lt.nyear1) iyear1=nyear1
      if(iyear2.lt.nyear1) iyear2=nyear1
      if(iyear1.gt.nyear2) iyear1=nyear2
      if(iyear2.gt.nyear2) iyear2=nyear2
      if(.false.) then
        print *,iyear1,ryear,iyear2
      endif
c ... interpolate latitude ...................
      fractl=(lat(ilat1)-rlat)/10
      do im=1,nmonth
        rmil1(im)=insol(im,ilat1,iyear1)+
     &    fractl*(insol(im,ilat2,iyear1)-insol(im,ilat1,iyear1))
        rmil2(im)=insol(im,ilat1,iyear2)+
     &    fractl*(insol(im,ilat2,iyear2)-insol(im,ilat1,iyear2))
        rmilba(im)=insol(im,ilat1,0)+
     &    fractl*(insol(im,ilat2,0)-insol(im,ilat1,0))
      enddo
      if(.false.) then
        print *,lat(ilat1),rlat,lat(ilat2)
        print *,ilat1,rl,ilat2
        print *,fractl
        print 100,lat(ilat1),(insol(im,ilat1,iyear1),im=1,4)
        print 101,rlat,(rmil1(im),im=1,4)
        print 100,lat(ilat2),(insol(im,ilat2,iyear1),im=1,4)
        print *
        print 100,lat(ilat1),(insol(im,ilat1,iyear2),im=1,4)
        print 101,rlat,(rmil2(im),im=1,4)
        print 100,lat(ilat2),(insol(im,ilat2,iyear2),im=1,4)
100     format(1x,i6,12i6)
101     format(1x,f6.1,12f6.1)
      endif
c ............................................
c ... interpolate year ........................................
      fracty=(ryear-iyear1)
      do im=1,nmonth
        rmil(im)=rmil1(im)+fracty*(rmil2(im)-rmil1(im))
      enddo
      if(.false.) then
        print *,iyear1,ryear,iyear2
        print *,fracty
        print 200,iyear1,(rmil1(im),im=1,4)
        print 201,ryear,(rmil(im),im=1,4)
        print 200,iyear2,(rmil2(im),im=1,4)
200     format(1x,i6,12f6.1)
201     format(1x,f6.1,12f6.1)
c        pause
      endif
c ..............................................................
c ... means and standard deviation for year
      rmeanno=average(nmonth,rmil)
      rmeanba=average(nmonth,rmilba)
c     if(dispon) print *,rlat,ryear,rmeanba,rmeanno,rmeanno-rmeanba
      rplusno=stdev(nmonth,rmil)
      rplusba=stdev(nmonth,rmilba)
      rmean=rmeanno-rmeanba
      rplus=rplusno-rplusba
c ... rmean is EXCESS insolation
c ... rmeanno is PAST insolation
c ... rmeanba is PRESENT insolation
c ... INCREASE PAST INSOLATION
      if(rmean.lt.0) then      
c ..... (Only increase colder)
        rmean=rmean*rmult
        rmeanno=rmeanba+rmean
      else
c ..... (Also increase warmer)
        rmean=rmean*rmult
        rmeanno=rmeanba+rmean
      endif
c ... adjust for latitude: > 45 too cold, <45 too warm
      if(rlat.gt.45) then
        rmeanno=rmeanno+rparam*((abs(rlat-45.))/45.)**0.5
      else
        rmeanno=rmeanno-rparam*((abs(rlat-45.))/45.)**0.5
      endif
c .....................................
c     if(dispon) print *,rlat,ryear,rplusba,rplus,rplus-rplusba
      tmeanno=stepbol(rmeanno)
      tplusno=stepbol(rmeanno+rplusno*amplify)
      tmeanba=stepbol(rmeanba)
      tplusba=stepbol(rmeanba+rplusba*amplify)
      if(dispon) print *,'lat,yr:',rlat,ryear
c     if(dispon) print *,rmil
      if(dispon) print *,'insolno',rmeanno,rplusno,rmeanno+rplusno
      if(dispon) print *,'insolba',rmeanba,rplusba,rmeanba+rplusba
      if(dispon) print *,'tempsno',tmeanno,tplusno-tmeanno,tplusno
      if(dispon) print *,'tempsba',tmeanba,tplusba-tmeanba,tplusba
c     pause
      tplusno=tplusno-tmeanno
      tplusba=tplusba-tmeanba
      rmean=rmeanno
      rplus=rplusno
      tmean=tmeanno
      tplus=tplusno
c ... rest is standard climatology ..............................
c ... lapse rate with elevation .................................
      tmean=tmean+acom*elev+TNSL+TNSLBASE
      if(dispon) print *,'surface mean annual air temp=',tmean,tplus
c ... calculate mean annual temp of free atmosphere-isothemal layer ...         
      tf=0.67d0*(tmean+273.0d0)+88.9d0
c     print *,'temp free at-isothrmal layer=',tf
c ... calculate saturation vapor pressure .............................
      term1=-9.09718d0*(273.16d0/tf-1.0d0)
      term2=-3.56654d0*log10(273.16d0/tf)
      term3=0.876793d0*(1.0d0-tf/273.16d0)+0.785835d0
      expon=term1+term2+term3
      es=10.d0**expon
c     print *,'saturation vapor pressure=',es
c ... calculate accumulation rate (m/yr) ..............................
      term1=www*es
      term2=xxx*slope
      term3=zzz
      term4=-15.276d0*shape
      acc=max(1d-9,term1+term2+term3+term4)
      if(dispon) print *,'accumulation rate=',acc
C ... CALCULATE ABLATION THE OLD WAY ...
c      if(.false.) then
        QY=0.d0
        DO I=1,12                                                      
          QY=QY+QI(I)-QS(I)*RLAT                                          
        ENDDO                                                          
        QY=QY/12.d0
        PDDM=0.d0
        DO I=1,12                                                      
C         TTT(I)=tmean+0.021d0*((QI(I)+QS(I)*RLAT)-QY)+8.954d0
          TTT(I)=tmean+0.021d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C         TTT(I)=tmean+0.018d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C         WRITE(19,*) I,TTT(I),QI(I)-QS(I)*RLAT,QY                        
          IF(TTT(I).GT.0.0) PDDM=PDDM+30.d0*TTT(I)                            
c          if(dispon) print *,i,TTT(I),pddm,max(0.0,TTT(i)*30*0.6)
        ENDDO
c      else
c ..... sum up positive degree days for ablation calc..............
        pdd=0
        do i=0,364,30
          temp=tmean+tplus*sin(2*pi*i/364.)
          if(temp.gt.0.) pdd=pdd+temp*30
c          if(dispon) print *,i,temp,pdd,max(0.0,temp*30*0.6)
        enddo
c      endif
      pdd=pddm
      abl=.6d0*pdd
      if(dispon) print *,'ablation rate=',abl,pdd,PDDM
c      print *,tmean,pdd,pddnew
c ... calculate net accumulation ......................................
      accnet=acc-abl
c     print *,'net accumulation/ablation=',accnet
c ... put in m/yr .....................................................
      accum=accnet*.01d0
      abl=abl*.01d0
c      if(dispon) pause
      end
c--------------------------------------------------
      function stepbol(rmean)
      IMPLICIT REAL*8(A-H,O-Z)
      data sigma /5.67d-5/
      isign=sign(1.d0,rmean)
c ... convert langleys/day to watt/m**2, multiply by 0.4843 .....
      tmean=abs(rmean)*0.4843d0
c ... convert w/m**2 to erg/cm**2/s, multipy by 1000. ...........
      tmean=tmean*1000.d0
c ... convert to temperature e=sigma*t**4, ......................
c ............  sigma=5.67e-5 erg/cm**2/deg**4/s ................
c ... tplus is temp with mean+stdev flux ........................
      albedo=0.12d0 ! dirt
      albedo=0.84d0 ! dry snow
      albedo=0.0d0 ! perfect black body
c ... tmean is temp with just mean ..............................
      tmean=(tmean*(1d0-albedo)/sigma)**0.25d0
c ... offset because stephan-boltzman is cold (greenhouse??) ....
      tmean=tmean+8.d0
c ... convert to deg c ..........................................
      tmean=tmean-273.16d0
      stepbol=isign*tmean
      end
c------------------------------------------------
      function average(n,r)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension r(n)
      sum=0.0d0
      do i=1,n
        sum=sum+r(i)
      enddo
      average=sum/n
      end
c------------------------------------------------
      function stdev(n,r)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension r(n)
      rmean=average(n,r)
      sum=0.0d0
      do i=1,n
        sum=sum+(r(i)-rmean)**2
      enddo
      sum=sum/n
      stdev=sqrt(sum)
      end
c------------------------------------------------
      function accmars0(ts0,elev,abl)
      implicit real*8(a-h,o-z)
c     DATA WWW /  19.1390686d0     /
      DATA WWW /  1.91390686d0     /
      DATA XXX / 0.922791243d0     /
      DATA ZZZ /-0.738900483d0     /
      toffset=60
      ts=ts0+toffset
      tmean=214d0
      pmean=5.6d0
      roverm=192d0
      grav=3.73
      p1=5.6d0
      z1=0d0
      tk=ts+273d0
      TTT=(TK+tmean)*0.5d0
      term=-grav*(elev-z1)/(roverm*TTT)
      pres=p1*exp(term)
C ... CALCULATE MEAN ANNUAL TEMP OF FREE ATMOSPHERE-ISOTHERMAL LAYER
      TF=0.67d0*(TS+273.0d0)+88.9d0
C ... CALCULATE SATURATION VAPOR PRESSURE
      es=satvap(tf)
c     es2=satvap(ts+273.d0)
      es2=es
C ... CALCULATE ACCUMULATION RATE (M/YR)
      TERM1=WWW*ES
      TERM2=XXX*SLOPE
      TERM3=ZZZ
      TERM4=-15.276d0*SHAPE
c
c **** <<< EXPERIMENTAL >>> **** turn off slope term ...
      TERM2=0
      TERM3=0
      TERM4=0
c **** <<< EXPERIMENTAL >>> ****
c
      ACC=max(0.d0,TERM1+TERM2+TERM3+TERM4)
C ... CALCULATE ABLATION
      factor=5
      aaa000=www
      wind=(20000-elev)/20000
      ABL=factor*wind*aaa000*es2/pres
C ... CALCULATE NET ACCUMULATION
      tmp=abl
      abl=acc*0.1d0
      acc=tmp*0.1d0
      ACCNET=ACC-ABL
      if(accnet.lt.0) then
        accnet=accnet*0.1
      endif
      accmars0= ACCNET
c     abl=acc-accnet
      return
      PRINT *,'-----------------------------'
      PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TK
      PRINT *,'Elevation                   =',elev
      PRINT *,'Pressure                    =',pres
      PRINT *,'Wind                        =',wind
      PRINT *,'TEMP FREE AT-ISOTHRMAL LAYER=',TF
      PRINT *,'SATURATION VAPOR PRESSURE   =',real(ES),real(es2)
      PRINT *,'ACCUMULATION                =',ACC*.01
      PRINT *,'ABLATION                    =',-ABL*.01
      PRINT *,'NET                         =',-ACCNET*.01
      end
      function satvap(tf)
      implicit real*8(a-h,o-z)
      TERM1=-9.09718d0*(273.16d0/TF-1.0d0)
      TERM2=-3.56654d0*LOG10(273.16d0/TF)
      TERM3=0.876793d0*(1.0d0-TF/273.16d0)+0.785835d0
      EXPON=TERM1+TERM2+TERM3
      satvap=10.d0**EXPON
      end
c------------------------------------------------
      function accmars1(ts0,elev,acc,abl)
c ... this is the original used in the paper ...
      implicit real*8(a-h,o-z)
c     DATA WWW /  19.1390686d0     /
      DATA WWW /  1.91390686d-0     /
      DATA XXX / 0.922791243d0     /
      DATA ZZZ /-0.738900483d0     /
      toffset=60
      toffset=0
      ts=ts0+toffset
      tmean=214d0
      pmean=5.6d0
      roverm=192d0
      grav=3.73
      p1=5.6d0
      z1=0d0
      tk=ts+273d0
      TTT=(TK+tmean)*0.5d0
      term=-grav*(elev-z1)/(roverm*TTT)
      pres=p1*exp(term)
C ... CALCULATE MEAN ANNUAL TEMP OF FREE ATMOSPHERE-ISOTHERMAL LAYER
      TF=0.67d0*(TS+273.0d0)+88.9d0
C ... CALCULATE SATURATION VAPOR PRESSURE
      es=satvap(tf)
      es2=satvap(ts+273.d0)
C ... CALCULATE ACCUMULATION RATE (M/YR)
      TERM1=WWW*ES
      TERM2=XXX*SLOPE
      TERM3=ZZZ
      TERM4=-15.276d0*SHAPE
c **** <<< EXPERIMENTAL >>> **** turn off slope term ...
      TERM2=0
      TERM3=0
      TERM4=0
c **** <<< EXPERIMENTAL >>> ****
      ACC=max(0.d0,TERM1+TERM2+TERM3+TERM4)
C ... CALCULATE ABLATION
      factor=10
      aaa000=www
      wind=(30000-elev)/30000
c     wind=1
      factor=10
c     pres=p1
c     print *,elev,ts,wind,p1/pres
      ABL=factor*wind*aaa000*es2*p1/pres
C ... CALCULATE NET ACCUMULATION
      acc=acc*0.1d0
      abl=abl*0.1d0
      if(.false.) then
        tmp=abl
        abl=acc
        acc=tmp
      endif
c     acc=es
c     abl=es2
      ACCNET=ACC-ABL
c     if(accnet.lt.0) then
c       accnet=accnet*0.1
c     endif
      accmars1= ACCNET
c     abl=acc-accnet
      return
      PRINT *,'-----------------------------'
      PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TK
      PRINT *,'Elevation                   =',elev
      PRINT *,'Pressure                    =',pres
      PRINT *,'Wind                        =',wind
      PRINT *,'TEMP FREE AT-ISOTHRMAL LAYER=',TF
      PRINT *,'SATURATION VAPOR PRESSURE   =',real(ES),real(es2)
      PRINT *,'ACCUMULATION                =',ACC*.01
      PRINT *,'ABLATION                    =',-ABL*.01
      PRINT *,'NET                         =',-ACCNET*.01
      end
c------------------------------------------------
      function accmars(ts0,elev0,acc,abl)
c ... this is the newest one ...
      implicit real*8(a-h,o-z)
      DATA WWW /  19.1390686d0     /
c     DATA WWW /  1.91390686d-0     /
      DATA XXX / 0.922791243d0     /
      DATA ZZZ /-0.738900483d0     /
      elev=elev0+3000
      abloffset=0.2    !  cm/yr (X10, mm/yr)
      toffset=60
      toffset=2.5
      ts=ts0+toffset
      tmean=214d0
      pmean=5.6d0
      roverm=192d0
      grav=3.74
      p1=5.6d0
      z1= 0
      tk=ts+273d0
      TTT=(TK+tmean)*0.5d0
      term=-grav*(elev-z1)/(roverm*TTT)
      pres=p1*exp(term)
C ... CALCULATE MEAN ANNUAL TEMP OF FREE ATMOSPHERE-ISOTHERMAL LAYER
      TF=0.67d0*(TS+273.0d0)+88.9d0
C ... CALCULATE SATURATION VAPOR PRESSURE
      es=satvap(tf)
      es2=satvap(ts+273.d0)
C ... CALCULATE ABLATION
      aaa000=www
      windf=(30000-elev)/30000
      windf=1
      factor=5
c     pres=p1
c     print *,elev,ts,windf,p1/pres,es,es2
      ABL=factor*windf*aaa000*es2*p1/pres
      ABL=factor*WWW*windf*ES2*p1/pres
c     ABL=ABL+(5e-6*elev)**2
c ... add 0.1 cm/yr (1 mm/yr) to ablation
      ABL=ABL+abloffset
c     ABL=ABL*wind(elev)
C ... CALCULATE ACCUMULATION RATE (M/YR)
      ACC=WWW*windf*ES
C ... CALCULATE NET ACCUMULATION
      acc=acc*0.01d0
      abl=abl*0.01d0
c     print *,elev,ts,p1,pres,acc,abl,acc-abl
      if(.false.) then
        tmp=abl
        abl=acc
        acc=tmp
      endif
c     acc=.1
c     abl=.05
      ACCNET=ACC-ABL
c     if(accnet.lt.0) then
c       accnet=accnet*0.1
c     endif
      accmars= ACCNET
c     print *,accnet,ts0
c     abl=acc-accnet
      return
      PRINT *,'-----------------------------'
      PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TK
      PRINT *,'Elevation                   =',elev
      PRINT *,'Pressure                    =',pres
      PRINT *,'Wind                        =',windf
      PRINT *,'TEMP FREE AT-ISOTHRMAL LAYER=',TF
      PRINT *,'SATURATION VAPOR PRESSURE   =',real(ES),real(es2)
      PRINT *,'ACCUMULATION                =',ACC*.01
      PRINT *,'ABLATION                    =',-ABL*.01
      PRINT *,'NET                         =',-ACCNET*.01
      pause
      end
      function wind(elev)
      implicit real*8(a-h,o-z)
      elevr=elev-10000
      if(elevr.gt.0) then
        wind=1+sign(1.d0,elevr)*(1e-2*elevr)**1
      else
        wind=1
      endif
      end 
      SUBROUTINE GAUSINIT(XI,ETA,W)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XI(2,9), ETA(2,9), W(2,9)
C     GAUSSIAN QUADRATURE OF ORDER THREE QUADRILATERALS
      XI(1,1)=-SQRT(3.d0/5.d0)
      XI(1,2)=0.d0
      XI(1,3)=-XI(1,1)
      XI(1,4)=XI(1,1)
      XI(1,5)=0.d0
      XI(1,6)=XI(1,3)
      XI(1,7)=XI(1,1)
      XI(1,8)=0.d0
      XI(1,9)=XI(1,3)
      ETA(1,1)=XI(1,1)
      ETA(1,2)=XI(1,1)
      ETA(1,3)=XI(1,1)
      ETA(1,4)=0.d0
      ETA(1,5)=0.d0
      ETA(1,6)=0.d0
      ETA(1,7)=XI(1,3)
      ETA(1,8)=XI(1,3)
      ETA(1,9)=XI(1,3)
      W(1,1)=25.d0/81.d0
      W(1,2)=40.d0/81.d0
      W(1,3)=W(1,1)
      W(1,4)=W(1,2)
      W(1,5)=64.d0/81.d0
      W(1,6)=W(1,2)
      W(1,7)=W(1,1)
      W(1,8)=W(1,2)
      W(1,9)=W(1,1)
C   GAUSSIAN QUADRATURE OF ORDER THREE TRIANGLES
      XI(2,1)=1.d0/3.d0
      XI(2,2)=2.d0/15.d0
      XI(2,3)=XI(2,2)
      XI(2,4)=11.d0/15.d0
      ETA(2,1)=XI(2,1)
      ETA(2,2)=XI(2,4)
      ETA(2,3)=XI(2,2)
      ETA(2,4)=ETA(2,3)
      W(2,1)=-27.d0/96.d0
      W(2,2)=25.d0/96.d0
      W(2,3)=W(2,2)
      W(2,4)=W(2,2)
      END
      SUBROUTINE NCONST(MXX,X,Y,KX,NTYPE,NUMEL,AFRACT,ASLDGB,
     &    LM,AFLOWA,BDROCK,DEPB,UNDEPB,PG,Q,CNEW,SLOPE,RHOI,WINDIR,
     &    WTHICK,ITYPE)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
C CALCULATES LINEARIZATION CONSTANT FROM CURRENT SOLUTION
      DIMENSION AFRACT(MXX),ASLDGB(MXX),AFLOWA(MXX),WTHICK(MXX)
      DIMENSION BDROCK(MXX),DEPB(MXX),UNDEPB(MXX),SLOPE(4,MXX)
      DIMENSION KX(MXX,4),X(MXX),Y(MXX),NTYPE(MXX),ITYPE(MXX)
      DIMENSION LM(5),CNEW(MXX)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2)
      DIMENSION PSI(4),DPSI(4,2)
      DIMENSION XY(2,4),WINDIR(2)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      REAL*8 Q(MXX)
      SAVE IPASS
      DATA IPASS /0/
      LOGICAL BATCH
      COMMON /BATCHER/ BATCH
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
c ... though the below originally was just a reminder, it is now
c     used to distinguish between interactive runs (a 1 is entered)
c     and batch runs (a 0 is entered). as such it reads a 1 value, but
c     writes a zero. there is then a 0 in the script file used fo
c     batches
      IF(IPASS.EQ.0) THEN
        IPASS=1
        print *,'USING THE 0.4 FACTOR IN THE TERM EQUATION (nconst.f)'
        IF(.FALSE.) THEN
          print *,'USING THE EISMINT EQUATION FOR SURFACE '
          print *,'      TEMP (accum.f)'
          print *,'USING THE EISMINT HARDNESS FUNCTION'
          print *,'           FOR FLOW CONST (temper.f)'
        ENDIF
        read(*,*) iwait
        write(99,*) 0
        BATCH=iwait.EQ.0
      ENDIF
      RHOR=4.0d0
      IF(BTOGG.ne.0) THEN
        FDEP1=1.d0
        FDEP2=0.d0
      ELSE
        FDEP1=RHOR/(RHOR-RHOI)
        FDEP2=RHOI/(RHOI-RHOR)
c experimental, turn off bed depression for BTOGG=0 case
        FDEP1=1.d0
        FDEP2=0.d0
      ENDIF
      DO J = 1,NUMEL
        IF(NTYPE(J).EQ.1) THEN
          NODEN=4
          CENTX=0.0d0
          CENTY=0.0d0
        ELSE
          NODEN=3
          CENTX=1.D0/3.D0
          CENTY=1.D0/3.D0
        ENDIF
        SUMHH=0.d0
        SUMX=0.d0
        SUMY=0.d0
        WTHIK=0.d0
        LITYPE=0
        DO I = 1,NODEN
          LM(I) = KX(J,I)
          LITYPE=MAX(LITYPE,ITYPE(LM(I)))
        ENDDO
        I=LM(1)
        JJ=LM(2)
        K=LM(3)
        L=LM(4)
        XY(1,1)=X(I)
        XY(1,2)=X(JJ)
        XY(1,3)=X(K)
        IF(NTYPE(J).EQ.1) XY(1,4)=X(L)
        XY(2,1)=Y(I)
        XY(2,2)=Y(JJ)
        XY(2,3)=Y(K)
        IF(NTYPE(J).EQ.1) XY(2,4)=Y(L)
        CALL FESHAPE(NTYPE(J),CENTX,CENTY,PSI,DPSI)
C
C ..... CALCULATE DXDS...EQUATION (5.3.6)
C
        DO I=1,2
          DO L=1,2
            DXDS(I,L)=0.0d0
            DO K=1,NODEN
              DXDS(I,L)=DXDS(I,L)+DPSI(K,L)*XY(I,K)
            ENDDO
          ENDDO
        ENDDO
C
C ..... CALCULATE DSDX...EQUATION (5.2.7)
C
        DETJ=(DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1))
        IF (DETJ.LE.0.0) THEN
          WRITE(12,5544) J,DETJ
          WRITE(*,5544) J,DETJ
          write(*,*) lm
          WRITE(12,5545) (JJ,XY(1,JJ),XY(2,JJ),JJ=1,4)
          WRITE(*,5545) (JJ,XY(1,JJ),XY(2,JJ),JJ=1,4)
5545      FORMAT(1X,I5,1X,1PE10.3,E10.3)
c          STOP
5544      FORMAT(' BAD JACOBIAN',I5,1PE10.3,/,1X,8E10.3)
        ENDIF
        DSDX(1,1)=DXDS(2,2)/DETJ
        DSDX(2,2)=DXDS(1,1)/DETJ
        DSDX(1,2)=-DXDS(1,2)/DETJ
        DSDX(2,1)=-DXDS(2,1)/DETJ
C
C ..... CALCULATE D(PSI)/DX...EQUATION (5.3.5)
C
        DO I=1,NODEN
          DPSIX(I)=DPSI(I,1)*DSDX(1,1)+DPSI(I,2)*DSDX(2,1)
          DPSIY(I)=DPSI(I,1)*DSDX(1,2)+DPSI(I,2)*DSDX(2,2)
        ENDDO
        IF(BTOGG.ne.0) THEN
          DO I = 1,NODEN
            SUMX = SUMX + Q(LM(I))*DPSIX(I)
            SUMY = SUMY + Q(LM(I))*DPSIY(I)
            WTHIK = WTHIK + WTHICK(LM(I))*PSI(I)
            THIK=Q(LM(I))-DEPB(LM(I))
            IF(DEPB(LM(I)).LT.SEALEV) THEN
              FLOT=(1.d0-RATDEN)*(DEPB(LM(I))-SEALEV)
              IF(Q(LM(I)).le.FLOT) THIK=0.d0
            ENDIF
c           if(thik.lt.0.) print *,'negative thickness at ',lm(i)
c            IF(THIK.GT.0.) SUMHH=SUMHH+THIK
            IF(THIK.GT.0.) SUMHH=SUMHH+THIK*PSI(I)
          ENDDO

        ELSE
          DO I = 1,NODEN
            SUMX = SUMX + Q(LM(I))*DPSIX(I)
            SUMY = SUMY + Q(LM(I))*DPSIY(I)
            WTHIK = WTHIK + WTHICK(LM(I))*PSI(I)
            IF(BDROCK(LM(I)).LE.-9999.) THEN
              DEPB(LM(I))=0.d0
            ELSE
              DEPB(LM(I))=FDEP2*Q(LM(I))+FDEP1*UNDEPB(LM(I))
            ENDIF
            THIK=(Q(LM(I))-UNDEPB(LM(I)))*FDEP1
            IF(DEPB(LM(I)).LT.SEALEV) THEN
              FLOT=(1.d0-RATDEN)*(DEPB(LM(I))-SEALEV)
              IF(Q(LM(I)).le.FLOT) THIK=0.d0
            ENDIF
c           if(thik.lt.0.) print *,'negative thickness at ',lm(i)
c            IF(THIK.GT.0.) SUMHH=SUMHH+THIK
            IF(THIK.GT.0.) SUMHH=SUMHH+THIK*PSI(I)
          ENDDO
        ENDIF
C
        DELH = SUMX**2 + SUMY**2
        DELH = SQRT(DELH)
        SLOPE(1,J)=DELH
        SLOPE(2,J)=SUMX
        SLOPE(3,J)=SUMY
        SLOPE(4,J)=SUMX*SIN(WINDIR(1))+SUMY*COS(WINDIR(1))
c        HH = SUMHH/DBLE(NODEN)
        HH = SUMHH
        if(LITYPE.eq.7) then
          TERM1 = AFRACT(J)*((PG/ASLDGB(J))**2)*(HH**3)*DELH
          TERM2 = (1.d0-AFRACT(J))*0.4d0*((PG/AFLOWA(J))**3)*
     &            (HH**5)*(DELH**2)
        else
c ....... experimental, eliminate use of fract ... 
c ....... link sliding directly to water thickness ...
c
          TERM1 = ((PG/ASLDGB(J))**2)*(HH**3)*DELH*RLUB(wthik)
          TERM2 = 0.4d0*((PG/AFLOWA(J))**3)*(HH**5)*(DELH**2)
          if(.false. .and. wthik.gt.-10) then
            us=0
            uf=0
            ut=0
            if(hh.gt.0) us=term1*delh/hh
            if(hh.gt.0) uf=term2*delh/hh
            if(hh.gt.0) ut=(term1+term2)*delh/hh
            if(us.gt.0) print *,real(us),real(uf),real(us+uf),
     &                          real(us/ut),real(wthik)
          endif
        endif
        CNEW(J) = TERM1 + TERM2
c      if(lm(1).eq.199.or.
c     &   lm(2).eq.199.or.
c     &   lm(3).eq.199.or.
c     &   lm(4).eq.199) then
c         print *,(lm(m),m=1,4)
c         print *,'s',(real(q(lm(m))),m=1,4)
c         print *,'b',(real(depb(lm(m))),m=1,4)
c         print *,'t',(real(q(lm(m))-depb(lm(m))),m=1,4)
c         print *,j,real(cnew(j)),real(hh),real(delh)
c         print *,real(sumx),real(sumy)
c       endif
c         if(hh.gt.0) then
c           print *,j,real(cnew(j)),real(hh),real(delh)
c            print *,j,real(term1),real(term2),real(term1/term2)
c          print *,(real(q(lm(m))),m=1,4)
c          print *,(real(depb(lm(m))),m=1,4)
c          print *,(real(q(lm(m))-depb(lm(m))),m=1,4)
c         endif
      ENDDO
      RETURN
      END
C==========================================
      FUNCTION RLUB(WTHICK)
c ... this function is used to calculate multiplying factor for sliding
c ... it appears in nconst.f with non-linear constant, and  in temper.f
c ... where heat generated at the bed due to sliding is calculated.
      IMPLICIT REAL*8(A-H,O-Z)
      RLUB=WTHICK*1d0
      END
C==========================================
      SUBROUTINE READN(MXX,HED,NUMNP,NUMEL,NUMGBC,NDT,INTER,DT,
     &                 KODE,X,Y,HTICE,ADOT,ADOTB,FRACT,PSURF,RHOI,
     &                 RHOW,RHOR,DEPB,
     &                 BDROCK,UNDEPB,FLOWA,ACON,SLDGB,TEMP,ITYPE,
     &                 AFUDGE,GEOFLUX,THICK,KX,CONST,IBFLUX,BFLUX,
     &                 QHOLD,NTYPE,NNODE,NCOL,AADOT,AFRACT,ABDRCK,
     &                 PPSURF,AFLOWA,ASLDGB,IDT,AMASS,NUMCOL,NUMLEV,
     &                 CALV,PCALV,WWW,WRATE,THICKL,WWWORIG,TSORIG,
     &                 TIME,ACC,ABLAT)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
C
C ... READ FROM SPLIT DATA SETS, PRODUCED BY UNPACK, REASSEMBLED BY PACK
C
      COMMON /SNOW/ SNOLIN,SNO(20)
      COMMON /TTNSL/ TTIME(NSMAX),TLIST(NSMAX),NTNSL
      COMMON /EXPER/ TPERIOD,TINIT,TVSTART,TVFINAL,IEXPER
      LOGICAL CTOGG,ITOGG,IOTOGG
      LOGICAL file67,file76
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      CHARACTER HED*80
      DIMENSION AMASS(11)
      DIMENSION KODE(MXX),X(MXX),Y(MXX),HTICE(MXX),ADOT(MXX),FRACT(MXX),
     &          PSURF(MXX),BDROCK(MXX),FLOWA(MXX),SLDGB(MXX),
     &          THICK(MXX),ADOTB(MXX),ITYPE(MXX),
     &          KX(MXX,4),CONST(MXX),IBFLUX(MXX,2),BFLUX(MXX),
     &          QHOLD(MXX),IDT(MXX),GEOFLUX(MXX),
     &          NNODE(MXX),NTYPE(MXX),AADOT(MXX),AFRACT(MXX),
     &          ABDRCK(MXX),AFUDGE(MXX),
     &          PPSURF(MXX),AFLOWA(MXX),ASLDGB(MXX),UNDEPB(MXX),
     &          LM(4),TEMP(MXX),ZERO(4),ACON(MXX),
     &          CALV(MXX),PCALV(MXX),WWW(3*MXX),WRATE(3*MXX,2),
     &          THICKL(MXX),DEPB(MXX),WWWORIG(MXX),TSORIG(MXX)
      DIMENSION ACC(MXX),ABLAT(MXX)
C ... TIMER STUFF, FOR SGI ONLY ...
      REAL*4 TB(2)
c     REAL*4 TB(2),ETIME,DTIME
c     EXTERNAL ETIME,DTIME
      DATA ZERO /4*0D0/
123   FORMAT(A25,T30,1PG13.6,G13.6)

      CALL SCENARIO(0.d0)
      DO I=1,NSMAX
        TTIME(I)=0.d0
        TLIST(I)=0.d0
        STIME(I)=0.d0
        SLIST(I)=0.d0
      ENDDO
      NSEAL=0
      DO I=1,NSMAX
        READ(75,*,END=889) STIME(I),SLIST(I)
        NSEAL=I
      ENDDO
      PRINT *,'INCREASE NUMBER OF STIME,STLIST',NSMAX
      STOP
889   CONTINUE
      IF(NSEAL.GT.0) THEN
        PRINT *,'THERE IS A SEA LEVEL LIST FOR THIS DATA SET'
        PRINT *,'THERE ARE',NSEAL,' TIMES IN THE LIST'
        PRINT *,'STARTING AT',STIME(1),' AND ENDING AT',STIME(NSEAL)
      ENDIF
      NTNSL=0
      DO I=1,NSMAX
        READ(72,*,END=888) TTIME(I),TLIST(I)
        NTNSL=I
      ENDDO
      PRINT *,'INCREASE NUMBER OF TTIME,TTLIST',NSMAX
      STOP
888   CONTINUE
      IF(NTNSL.GT.0) THEN
        PRINT *,'THERE IS A TNSL LIST FOR THIS DATA SET'
        PRINT *,'THERE ARE',NTNSL,' TIMES IN THE LIST'
        PRINT *,'STARTING AT',TTIME(1),' AND ENDING AT',TTIME(NTNSL)
      ENDIF
      DO I=1,20
        SNO(I)=0.d0
      ENDDO
C
C ... READ INPUT HEADER
      READ(30,1000,END=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &                        INTER,DT
1000  FORMAT (A80,/,7I6,F8.0)
C
C ... READ INPUT GRID, THINGS THAT NEVER CHANGE
      READ(31) HED
      READ(31) (KODE(I),I=1,NUMNP)
      READ(31) (X(I),I=1,NUMNP)
      READ(31) (Y(I),I=1,NUMNP)
      READ(31) (PSURF(I),I=1,NUMNP)
      READ(31) (BDROCK(I),I=1,NUMNP)
      READ(31) (KX(I,1),KX(I,2),KX(I,3),KX(I,4),I=1,NUMEL)
      READ(31) (IBFLUX(I,1),IBFLUX(I,2),BFLUX(I),I=1,NUMGBC)
C
C ... READ INPUT DIFFER, THINGS THAT DIFFER FROM ONE GRID TO THE NEXT
      READ(32) HED
      READ(32) (ADOT(I),I=1,NUMNP)
      READ(32) (FRACT(I),I=1,NUMNP)
      READ(32) (FLOWA(I),I=1,NUMNP)
      READ(32) (SLDGB(I),I=1,NUMNP)
      READ(32) (TEMP(I),I=1,NUMNP)
      READ(32) (ITYPE(I),I=1,NUMNP)
      READ(32) (AFUDGE(I),I=1,NUMNP)
      READ(32) (GEOFLUX(I),I=1,NUMNP)
      IF(GEOFLUX(1).GT.0. .AND. GEOFLUX(1).LE.10.) THEN
        PRINT *,'THIS IS OLD DATA SET, CHANGE GEOTHERMAL'
c        PAUSE
      ENDIF
      READ(32) (CALV(I),I=1,NUMNP)
      DO I=1,NUMNP
        TSORIG(I)=TEMP(I)
        IF(TSORIG(I).EQ.-999.) TSORIG(I)=0.0d0
        IF(AFUDGE(I).EQ.0.) AFUDGE(I)=1.D0
        IF(FLOWA(I).LE.0.) FLOWA(I)=1D-6
        IF(SLDGB(I).LE.0.) SLDGB(I)=1D-9
        IF(CALV(I).LT.0.) CALV(I)=0.01d0
c ************EISMINT GEOTHERMAL FLUX VALUES ********************
c ... 42 milliWatts/m**2 sliding experiment
c        IF(GEOFLUX(I).EQ.0.) GEOFLUX(I)=3.17D5
c ... 50 milliWatts/m**2 greenland experiment 
c        IF(GEOFLUX(I).EQ.0.) GEOFLUX(I)=3.77D5
c ... 54.6 milliWatts/m**2 ANTARCTIC VALUE (HUYBRECHTS) 
        IF(GEOFLUX(I).EQ.0.) GEOFLUX(I)=4.121D5
c ... 25 milliWatts/m**2 greenland experiment 
         IF(GEOFLUX(I).EQ.0.) GEOFLUX(I)=3.77D5/2
c        IF(GEOFLUX(I).LE.10.) THEN
c          PRINT *,'THIS IS OLD DATA SET, CHANGE GEOTHERMAL'
c          GEOFLUX(I)=4.121D5
c        ENDIF
      ENDDO
C
C ... READ INPUT TIME, THINGS THAT CHANGE WITH TIME
      READ(33) HED
      READ(33) (HTICE(I),I=1,NUMNP)
      READ(33) (ACC(I),I=1,NUMNP)
      READ(33) (ABLAT(I),I=1,NUMNP)
      READ(33) (CONST(I),I=1,NUMEL)
      READ(33) (ACON(I),I=1,NUMEL)
      DO I=1,NUMEL
        IF(ACON(I).LT.1E-30) ACON(I)=0.D0
      ENDDO
      IF(NUMNP.GT.MXX) THEN
        PRINT *,'NUMNP=',NUMNP,' MXX=',MXX,' INCREASE MXX'
        STOP
      ENDIF
      DO N=1,NUMNP
        ADOTB(N)=ADOT(N)
        IDT(N)=0
c ... turn off for madeleine/forget Mass Balance
        IF(.false. .and. ADOT(N).LT.-100.) THEN
         IDT(N)=-INT(100+ADOT(N))
         ADOT(N)=AFUNCT(TIME,IDT(N),AMASS,HTICE(N),BDROCK(N),PSURF(N),
     &                  ZERO,X(N),Y(N),TEMP(N))
        ENDIF
c       IF(PSURF(N).LT.BDROCK(N)-SEALEV) PSURF(N)=BDROCK(N)
        THICK(N)=HTICE(N)-BDROCK(N)
      ENDDO
C ... UNLOAD THE BED ...
      IF(BTOGG.ne.0) THEN
        WRITE(7,123) ' TIME BEFORE 1ST PLATE ',ETIME(TB),DTIME(TB)
        print *,'unloading bed with present surface'
C----------------------------------------------------------------!
        LDEP=0                                                   1
        DO I=1,NUMNP                                             !
          THICKL(I)=PSURF(I)-BDROCK(I)                           !
          if(bdrock(i).lt.SEALEV) then
             FLOT=(1.D0-RATDEN)*(BDROCK(I)-SEALEV)                     !
             IF(PSURF(I).LT.FLOT) THEN                           !
               THICKL(I)=0.D0                                    !
               THICKL(I)=FLOT-BDROCK(I)                          !
             ENDIF                                               !
          endif
          IF(BDROCK(I).le.-9999.) THICKL(I)=0.D0                                       !
          IF(THICKL(I) .GT. 0.D0) LDEP=1                         !
        ENDDO                                                    !
        IF(LDEP.EQ.1) THEN                                       !
          if(BTOGG.eq.1) then
            CALL SPLATE(-1,NUMNP,THICKL,UNDEPB,DEPB,RHOI,RHOR,RHOW,
     &                  0.d0,TIME,wrate,www,wwworig)
          else
            CALL OPLATE(-1,MXX,3*MXX,NUMNP,NUMEL,X,Y,KX,THICKL,   !
     &                 0.d0,WWW,WRATE,WMIN,0.d0,WWWORIG)          !
          endif
          IF(.FALSE.) THEN                                       !
            CALL POUTSTF(NUMNP*3,WWW,WRATE(1,1),NUMNP,THICKL,    !
     &                   NUMCOL,NUMLEV,0.d0,WWWORIG)             !
          ENDIF                                                  !
c         sum=0.                                                 !
          DO I=1,NUMNP      
            if(.false.) then
              UNDEPB(I)=BDROCK(I)-WWW((I-1)*3+1)                   !
            else
              UNDEPB(I)=BDROCK(I)-WWWORIG(I)                       !
            endif
            if(.true.) then
              DEPB(I)=UNDEPB(I)+WWW((I-1)*3+1)                     !
            else
              DEPB(I)=UNDEPB(I)+WWWORIG(I)                         !
            endif

c            DEPB(I)=BDROCK(I)
            UNDEPBI=RHOI*(PSURF(I)-BDROCK(I))/RHOR+BDROCK(I)
c           print *,real(WWW((I-1)*3+1)),REAL(UNDEPB(I)),REAL(depb(i))                                     !
c           PRINT *,I,REAL(UNDEPB(I)),REAL(UNDEPBI),
c    &            REAL(bdrock(I)-(UNDEPBI))/real(WWW((I-1)*3+1))
c           sum=sum+(bdrock(I)-(UNDEPBI))/(WWW((I-1)*3+1))
          ENDDO                  
c         print *,sum/numnp
        ENDIF                                                    !
        WRITE(7,123) ' TIME AFTER  1ST PLATE ',ETIME(TB),DTIME(TB)
C----------------------------------------------------------------!
      ELSEIF(BTOGG.eq.0) then
C ... OLD WAY ...
        DO I=1,NUMNP
          IF(BDROCK(I).GT.-9999.) THEN
            UNDEPB(I)=RHOI*(PSURF(I)-BDROCK(I))/RHOR+BDROCK(I)
            UNDEPB(I)=BDROCK(I)
          ENDIF
          DEPB(I)=UNDEPB(I)
        ENDDO
      else
        print *,' problems with BTOGG',BTOGG
        pause
      ENDIF
      NELMAX=-10000
      NCLMAX=-10000
      DO N=1,NUMEL
        NTYPE(N)=1
        IF(KX(N,4).EQ.0) NTYPE(N)=2
        NODEN=4
        NNODE(N)=4
        IF(NTYPE(N).EQ.2) NODEN=3
        IF(NTYPE(N).EQ.2) NNODE(N)=3
        IMAX=-10000
        IMIN=10000
        DO I=1,NODEN
          IF(KX(N,I).GT.IMAX) IMAX=KX(N,I)
          IF(KX(N,I).LT.IMIN) IMIN=KX(N,I)
        ENDDO
C **** FOR MAP1 WITH ASYMSL ****
        NCOL1=2*(IMAX-IMIN)+1
        IF(NCOL1.GT.NCLMAX) THEN
          NCLMAX=NCOL1
          NELMAX=N
        ENDIF
      ENDDO
      WRITE(*,*) 'MAX NCOL=',NCLMAX,' IN ELEMENT',NELMAX
c     IF(NCLMAX.GT.NCOL) THEN
c       WRITE(*,*) 'CURRENTLY USING NCOL=',NCOL,' INCREASE'
c       STOP
c     ENDIF
C
      SUMC = 0.0d0
      DO N = 1,NUMEL
        QHOLD(N) = 0.0d0
c       IF (CONST(N).EQ.0D0 .AND. DT.EQ.0D0) CONST(N)=1.92d7
        IF (CONST(N).EQ.0D0) CONST(N)=1.92d7
c       CONST(N)=1.92d7
        SUMC = SUMC + CONST(N)
        IF(NTYPE(N).EQ.1) THEN
          NODEN=4
          NINT=9
        ELSE
          NODEN=3
          NINT=4
        ENDIF
        AAADOT=0.d0
        AAFRCT=0.d0
        AABDRK=0.d0
        APSURF=0.d0
        AAFLOW=0.d0
        AASLDG=0.d0
        ACALV=0.d0
        DO I=1,4
          LM(I)=KX(N,I)
          AAADOT=AAADOT+ADOT(LM(I))
          AAFRCT=AAFRCT+FRACT(LM(I))
          AABDRK=AABDRK+BDROCK(LM(I))
          APSURF=APSURF+PSURF(LM(I))
          AAFLOW=AAFLOW+FLOWA(LM(I))
          AASLDG=AASLDG+SLDGB(LM(I))
          ACALV=ACALV+CALV(LM(I))
        ENDDO
        DENOM=1.d0/DBLE(NODEN)
        AAADOT=AAADOT*DENOM
        AADOT(N)=AAADOT
        AAFRCT=AAFRCT*DENOM
        AFRACT(N)=AAFRCT
        AABDRK=AABDRK*DENOM
        ABDRCK(N)=AABDRK
        APSURF=APSURF*DENOM
        PPSURF(N)=APSURF
        ACALV=ACALV*DENOM
        PCALV(N)=ACALV
        AAFLOW=AAFLOW*DENOM
        IF(ACON(N).EQ.0.) THEN
          AFLOWA(N)=AAFLOW
        ELSE
          AFLOWA(N)=ACON(N)
        ENDIF
        AASLDG=AASLDG*DENOM
        ASLDGB(N)=AASLDG
      ENDDO
C
      SUMC = SUMC/NUMEL
      WRITE(*,*) SUMC
      if(.false.) then
        CALL WRITPROF(NUMNP,X,Y,HTICE,BDROCK,ADOT)
        pause
      endif

      RETURN
999   CONTINUE
      PRINT *,'NO DATA SET'
      STOP
      END
c--------------------------------------------
      SUBROUTINE READCOORD(NUMNP)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      common /coord/ CLAT(MAXNUM),CLONG(MAXNUM)
      READ(53,end=100) HED
      print *,'reading coords'
      READ(53) (CLAT(I),I=1,NUMNP)
      READ(53) (CLONG(I),I=1,NUMNP)
      return
100   continue
      print *,'no coord file available'
      END
      SUBROUTINE FESHAPE(NTYPE,XI,ET,PSI,DPSI)
C ELEMENT SHAPE FUNCTIONS AND DERIVATIVES AT LOCAL COORDINATES (XI,ET)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PSI(4),DPSI(4,2)
      IF(NTYPE.EQ.1) THEN
        PSI(1)=.25d0*(1.d0-XI)*(1.d0-ET)
        PSI(2)=.25d0*(1.d0+XI)*(1.d0-ET)
        PSI(3)=.25d0*(1.d0+XI)*(1.d0+ET)
        PSI(4)=.25d0*(1.d0-XI)*(1.d0+ET)
        DPSI(1,2)=-.25d0*(1.d0-XI)
        DPSI(2,2)=-.25d0*(1.d0+XI)
        DPSI(3,2)=.25d0*(1.d0+XI)
        DPSI(4,2)=.25d0*(1.d0-XI)
        DPSI(1,1)=-.25d0*(1.d0-ET)
        DPSI(2,1)=.25d0*(1.d0-ET)
        DPSI(3,1)=.25d0*(1.d0+ET)
        DPSI(4,1)=-.25d0*(1.d0+ET)
      ELSE
        PSI(1)=1.d0-XI-ET
        PSI(2)=XI
        PSI(3)=ET
        DPSI(1,2)=-1.d0
        DPSI(2,2)=0.d0
        DPSI(3,2)=1.d0
        DPSI(1,1)=-1.d0
        DPSI(2,1)=1.d0
        DPSI(3,1)=0.d0
      ENDIF
      END
      SUBROUTINE ELPROP(MXX, NUMEL, NTYPE, KX, ADOT, 
     &           AADOT, FRACT,
     &           AFRACT, BDROCK, ABDRCK, PSURF, PPSURF, FLOWA, ACON,
     &           ICON, AFLOWA, SLDGB, ASLDGB,CALV,PCALV)
      IMPLICIT REAL*8(A-H,O-Z)
C ... LOAD NODAL PROPERTIES INTO ELEMENT PROPERTIES
      DIMENSION NTYPE(MXX),ADOT(MXX),AADOT(MXX),FRACT(MXX),AFRACT(MXX),
     &          BDROCK(MXX),ABDRCK(MXX),PSURF(MXX),PPSURF(MXX),
     &          FLOWA(MXX),AFLOWA(MXX),SLDGB(MXX),ASLDGB(MXX),
     &          LM(4),KX(MXX,4),ACON(MXX),CALV(MXX),PCALV(MXX)
C      REWIND 90
c      PRINT *,'IN ELPROP'
      DO N = 1,NUMEL
        IF(NTYPE(N).EQ.1) THEN
          NODEN=4
        ELSE
          NODEN=3
        ENDIF
        AAADOT=0.d0
        AAFRCT=0.d0
        AABDRK=0.d0
        APSURF=0.d0
        AAFLOW=0.d0
        AASLDG=0.d0
        ADNSTY=0.d0
        ACALV=0.d0
        DO I=1,4
          LM(I)=KX(N,I)
          AAADOT=AAADOT+ADOT(LM(I))
          AAFRCT=AAFRCT+FRACT(LM(I))
          AABDRK=AABDRK+BDROCK(LM(I))
          APSURF=APSURF+PSURF(LM(I))
          AAFLOW=AAFLOW+FLOWA(LM(I))
          AASLDG=AASLDG+SLDGB(LM(I))
          ACALV=ACALV+CALV(LM(I))
        ENDDO
        DENOM=1.d0/DBLE(NODEN)
        AAADOT=AAADOT*DENOM
        AADOT(N)=AAADOT
        AAFRCT=AAFRCT*DENOM
        AFRACT(N)=AAFRCT
c        WRITE(90,*) N,AFRACT(N)
        AABDRK=AABDRK*DENOM
        ABDRCK(N)=AABDRK
        ADNSTY=ADNSTY*DENOM
        PPSURF(N)=APSURF
        ACALV=ACALV*DENOM
c ..... calving is AVERAGE of 4 nodal values...
c       PCALV(N)=ACALV
c ..... calving is MINIMUM of 4 nodal values...
        PCALV(N)=MIN(CALV(LM(1)),CALV(LM(2)),CALV(LM(3)),CALV(LM(4)))
c ..... calving is MAXIMUM of 4 nodal values...
c       PCALV(N)=MAX(CALV(LM(1)),CALV(LM(2)),CALV(LM(3)),CALV(LM(4)))
        AAFLOW=AAFLOW*DENOM
        IF(ACON(N).EQ.0. .OR. ICON.EQ.0) THEN
          AFLOWA(N)=AAFLOW
        ELSE
          AFLOWA(N)=ACON(N)
        ENDIF
        AASLDG=AASLDG*DENOM
        ASLDGB(N)=AASLDG
      ENDDO
      RETURN
      END
c==========================================================
      SUBROUTINE CALVING(MXX, NUMEL, NTYPE, KX, HTICE, BED, ADOT, 
     &           AADOT,CFACTOR,ADC,CALV,PCALV,TNSL)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
C ... LOAD NODAL PROPERTIES INTO ELEMENT PROPERTIES
      DIMENSION NTYPE(MXX),ADOT(MXX),AADOT(MXX),
     &          HTICE(MXX),BED(MXX),ADC(MXX),
     &          LM(4),KX(MXX,4),CALV(MXX),PCALV(MXX)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      DATA RHOI /0.917D0/, RHOW /1.092D0/
      IF(CFACTOR.EQ.0.) RETURN
      NCALV=0
      acalv=0.d0
      if(.true.) then
        rmult=1.d0
      else
c ... make calving depend on tnsl X 0 at -10, X 1 at 0
        ttemp=4.5d0
        ttemp=6.0d0
        rmult=(ttemp+tnsl)/ttemp
        rmult=min(rmult,1.d0)
        rmult=max(rmult,0.d0)
      endif
      if(.false.) then
c ..... old way ......
        DO N = 1,NUMEL
          IF(CFACTOR.LT.0) THEN
            CFACTLOC=PCALV(N)
          ELSE
            CFACTLOC=CFACTOR
            CFACTLOC=CFACTOR*rmult
          ENDIF
          IF(NTYPE(N).EQ.1) THEN
            NODEN=4
          ELSE
            NODEN=3
          ENDIF
          AABED=0.d0
          AAHTICE=0.d0
          DO I=1,4
            LM(I)=KX(N,I)
            AABED=AABED+BED(LM(I))
            AAHTICE=AAHTICE+HTICE(LM(I))
          ENDDO
          DENOM=1.d0/DBLE(NODEN)
          AABED=AABED*DENOM
          AAHTICE=AAHTICE*DENOM
C ....... THIS IS FOR CALVING ...
          IF(AABED.LT.SEALEV) THEN
            FLOT=(1.d0-RATDEN)*(AABED-SEALEV)
            THIK=AAHTICE-FLOT
            IF(THIK.LT.0.) THEN
              NCALV=NCALV+1
              DC=CFACTLOC*THIK
              AADOT(N)=MIN(DC,AADOT(N))
              DO I=1,4
                ADC(LM(I))=DC
                ADOT(LM(I))=AADOT(N)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      else
c ..... new way ......
        DENOM=0.25d0
        DO N = 1,NUMEL
          IF(CFACTOR.LT.0) THEN
            CFACTLOC=PCALV(N)
c            CFACTLOC=PCALV(N)*rmult
          ELSE
            CFACTLOC=CFACTOR
c            CFACTLOC=CFACTOR*rmult
          ENDIF
          AABED=0.d0
          AAHTICE=0.d0
          icalv=0
          DO I=1,4
            LM(I)=KX(N,I)
            AABED=AABED+BED(LM(I))
            AAHTICE=AAHTICE+HTICE(LM(I))
            if(BED(LM(I)).lt.SEALEV) then
              FLOT=(1.d0-RATDEN)*(BED(LM(I))-SEALEV)
              if(HTICE(LM(I)).lt.FLOT) icalv=icalv+1
            endif
          ENDDO

c ... old way
          if(icalv.gt.0 .and. icalv.lt.5 .and. .false.) then
            AABED=AABED*DENOM
            AAHTICE=AAHTICE*DENOM
C ......... THIS IS FOR CALVING ...
            IF(AABED.LT.SEALEV) THEN
              FLOT=(1.d0-RATDEN)*(AABED-SEALEV)
              THIK=AAHTICE-FLOT
              IF(THIK.LT.0.) THEN
                DC=CFACTLOC*THIK
                if(icalv.lt.4) then
                  NCALV=NCALV+1
                  acalv=acalv+dc
                endif
c              AADOT(N)=MIN(DC,AADOT(N))
              AADOT(N)=DC+AADOT(N)
                DO I=1,4
                  ADC(LM(I))=DC
                  ADOT(LM(I))=AADOT(N)
                ENDDO
              ENDIF
            ENDIF
          endif

c ... new way 28 feb 2003 ...
          if(icalv.gt.0 .and. icalv.lt.5) then
            AABED=AABED*DENOM
            AAHTICE=AAHTICE*DENOM
C ......... THIS IS FOR CALVING ...
            IF(AABED.LT.SEALEV) THEN
              FLOT=(1.d0-RATDEN)*(AABED-SEALEV)
c ... following is weertman spreading of unconstrained ice shelf (Z**3)
                DC=(1.1e-7*((1.d0-RATDEN)*(FLOT-AABED))**3)*CFACTLOC 
c ... fudge to make CFACTLOC roughly 1.0 ...
                DC=DC*10    
                dcrmult=dc*rmult**2
                if(icalv.lt.4) then
                  NCALV=NCALV+1
                  acalv=acalv+dcrmult
                endif
                DO I=1,4
                  if(BED(LM(I)).lt.SEALEV) then
                    FLOT=(1.d0-RATDEN)*(BED(LM(I))-SEALEV)
                    if(HTICE(LM(I)).lt.FLOT) then
                      ADC(LM(I))=dcrmult
                      ADOT(LM(I))=max(-100.d0,ADOT(LM(I))+dcrmult)
c                      ADOT(LM(I))=dcrmult
                    endif
                  endif
                ENDDO
            ENDIF
          endif


        ENDDO
      endif
      IF(IOTOGG) then
        write(list(ipage+1),*) NCALV,'= NUMBER CALVING...',real(acalv),
     &                         real(rmult)
        ipage=ipage+1
      endif
      END
      SUBROUTINE NODESL(MXX,NUMNP,NUMEL,KX,SLOPE,SLOPN)
C ... CALCULATES NODE SLOPES FROM ELEMENT SLOPES 
C ... FOR ACCUMULATION PARAMETERIZATIONS
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(MMXX=MAXNUM)
      DIMENSION KX(MXX,4),SLOPE(4,MXX),SLOPN(4,MXX),ICOUNT(MMXX)
      if(.false.) then ! average of element slopes
        DO I=1,NUMNP
          DO L=1,4
            SLOPN(L,I)=0.d0
          ENDDO
          ICOUNT(I)=0
        ENDDO
        DO I=1,NUMEL
          DO J=1,4
            DO L=1,4
              SLOPN(L,KX(I,J))=SLOPN(L,KX(I,J))+SLOPE(L,I)
            ENDDO
            ICOUNT(KX(I,J))=ICOUNT(KX(I,J))+1
          ENDDO
        ENDDO
        DO I=1,NUMNP
          IF(ICOUNT(I).GT.0) THEN
            DO L=1,4
              SLOPN(L,I)=SLOPN(L,I)/DBLE(ICOUNT(I))
            ENDDO
c           SLOPN(1,I)=SQRT(SLOPN(2,I)**2+SLOPN(3,I)**2)
          ENDIF
        ENDDO
      else                     ! minimum of element slopes
        DO I=1,NUMNP
          DO L=1,4
            SLOPN(L,I)=1d30
          ENDDO
          ICOUNT(I)=0
        ENDDO
        DO I=1,NUMEL
          DO J=1,4
            DO L=1,4
              SLOPN(L,KX(I,J))=min(SLOPN(L,KX(I,J)),SLOPE(L,I))
            ENDDO
            ICOUNT(KX(I,J))=ICOUNT(KX(I,J))+1
          ENDDO
        ENDDO
        DO I=1,NUMNP
          IF(ICOUNT(I).eq.0) THEN
            print *,'problems in node',i,' with nodesl'
          ENDIF
        ENDDO
      endif
      RETURN
      END
      SUBROUTINE SCALEXY(xx,yy,n,xmin,xmax,ymin,ymax,delx,dely)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION xx(n),yy(n)
      logical zflag
      common /zoomcom/ pxmin,pymin,pxmax,pymax,zflag
      if(zflag) then
        XMIN=PXMIN
        YMIN=PYMIN
        XMAX=PXMAX
        YMAX=PYMAX
      else
        XMIN=1.d30
        YMIN=XMIN
        XMAX=-XMIN
        YMAX=-YMIN
        DO I=1,n
          XMAX=MAX(XMAX,XX(I))
          YMAX=MAX(YMAX,YY(I))
          XMIN=MIN(XMIN,XX(I))
          YMIN=MIN(YMIN,YY(I))
        ENDDO
      endif
      IF(XMAX-XMIN.GT.YMAX-YMIN) THEN
        YMAX=YMIN+XMAX-XMIN
      ELSE
        XMAX=XMIN+YMAX-YMIN
      ENDIF
      DELX=(XMAX-XMIN)/5.d0
      DELY=(YMAX-YMIN)/5.d0
      IF(XMIN.EQ.XMAX) XMAX=XMIN+1.d0
      IF(YMIN.EQ.YMAX) YMAX=YMIN+1.d0
      if(.not.zflag) then
        xmin=xmin*1.d-3
        xmax=xmax*1.d-3
        ymin=ymin*1.d-3
        ymax=ymax*1.d-3
        delx=delx*1.d-3
        dely=dely*1.d-3
      endif
c      print *,'x:',xmin,xmax,delx
c      print *,'y:',ymin,ymax,dely
c      pause
      END
c
      SUBROUTINE SCALE3(AR, RR, NN, II)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AR(NN+2)
      AMIN=1.d30
      AMAX=-1.d30
      DO I=1,NN
        AMAX=MAX(AR(I),AMAX)
        AMIN=MIN(AR(I),AMIN)
      ENDDO
      IF(AMIN.EQ.AMAX) AMAX=AMIN+1.d0
      AR(NN+1)=AMIN
      AR(NN+2)=(AMAX-AMIN)/RR
      AMAX=AMAX+50.d0*AR(NN+2)
      AMIN=AMIN-50.d0*AR(NN+2)
      AR(NN+1)=AMIN
      AR(NN+2)=(AMAX-AMIN)/RR
      RETURN
      END
C
      SUBROUTINE SCALE2(AR, RR, NN, II)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AR(NN+2)
      AMIN=1.d30
      AMAX=-1.d30
      DO I=1,NN
        AMAX=MAX(AR(I),AMAX)
        AMIN=MIN(AR(I),AMIN)
      ENDDO
      IF(AMIN.EQ.AMAX) AMAX=AMIN+1.d0
      AR(NN+1)=AMIN
      AR(NN+2)=(AMAX-AMIN)/RR
      RETURN
      END
C
      SUBROUTINE SORTIX(N,IX)
      DIMENSION IX(N)
      LOGICAL ANYEXC
      ANYEXC=.TRUE.
C WHILE ANY EXCHANGES
30    IF(ANYEXC) THEN
        ANYEXC=.FALSE.
        DO I=1,N-1
          IF(IX(I).LT.IX(I+1)) THEN
            CALL SWAP(IX(I),IX(I+1))
            ANYEXC=.TRUE.
          ENDIF
        ENDDO
        GOTO 30
      ENDIF
      RETURN
      END
C
      SUBROUTINE SWAP(IX1,IX2)
      ITEMP=IX1
      IX1=IX2
      IX2=ITEMP
      RETURN
      END
C
      SUBROUTINE ELMDUP(NN,IX)
      DIMENSION IX(NN)
      NNOUT=NN
      DO N=2,NN
80      IF(IX(N-1).EQ.IX(N) .AND. (IX(N).NE.0)) THEN
          NNOUT=NNOUT-1
          DO J=N,NN-1
            IX(J)=IX(J+1)
          ENDDO
          IX(NN)=0
          GOTO 80
        ENDIF
      ENDDO
      NN=NNOUT
      RETURN
      END
      INTEGER FUNCTION IFIND(X, Y, XC, YC)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XC(5),YC(5),TEST(4)
      XC(5)=XC(1)
      YC(5)=YC(1)
      DO I=1,4
        AX=X-XC(I)
        AY=Y-YC(I)
        BX=XC(I+1)-XC(I)
        BY=YC(I+1)-YC(I)
        TEST(I)=-(AX*BY-AY*BX)
      ENDDO
      IF((TEST(1).GE.0.) .AND.
     &   (TEST(2).GE.0.) .AND.
     &   (TEST(3).GE.0.) .AND.
     &   (TEST(4).GE.0.))THEN
        IFIND=1
      ELSE
        IFIND=0
      ENDIF
      RETURN
      END
      SUBROUTINE VOLUME(MXX,TIME,NUMNP,NUMEL,X,Y,KX,Q,BDROCK,DEPB,
     &                  ADOT,RHOI,RHOW,VOL,AREATOT,AMASS,
     &                  GRIDAREA,iprt,sealow)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
C ... CALCULATES VOLUMES (FLOTATION AND TOTAL) AND AREA
      DIMENSION BDROCK(MXX),KX(MXX,4),X(MXX),Y(MXX),LM(4)
      DIMENSION DEPB(MXX),ADOT(MXX)
      REAL*8 Q(MXX)
      DIMENSION AMASS(11)
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      logical iprt
      RNET=0.0d0
      VOL=0.0d0
      VOL1=0.d0
      AREATOT=0.0d0
      GRIDAREA=0.0d0
      THICKMAX=0.0d0
      ITHMAX=0
      DO I=1,NUMEL
        SUMH=0.d0
        SUMT=0.d0
        SUMA=0.d0
        IC=0
        DO J=1,4
          LM(J)=KX(I,J)
          IF(BDROCK(LM(J)).LE.-9999.) GOTO 100
          IF(DEPB(LM(J)).LT.SEALEV) THEN
            FLOT=(1.d0-RATDEN)*(DEPB(LM(J))-SEALEV)
          ELSE
            FLOT=DEPB(LM(J))
          ENDIF
          HEIGHT=Q(LM(J))-FLOT
          THICK=Q(LM(J))-DEPB(LM(J))
          IF(HEIGHT.GT.1.) SUMH=SUMH+HEIGHT
          IF(THICK.GT.1. .and. Q(LM(J)).GT.SEALEV+1E-6) THEN
            SUMT=SUMT+THICK
            SUMA=SUMA+ADOT(LM(J))
            IC=IC+1
          ENDIF
        ENDDO
        HEIGHT=SUMH*0.25d0
        THICK=SUMT*0.25d0
        ADOTAVG=SUMA*0.25d0
        AREA=0.5d0*((X(LM(2))-X(LM(1)))*(Y(LM(3))-Y(LM(2)))-
     &            (X(LM(3))-X(LM(2)))*(Y(LM(2))-Y(LM(1)))+
     &            (X(LM(4))-X(LM(3)))*(Y(LM(1))-Y(LM(4)))-
     &            (X(LM(1))-X(LM(4)))*(Y(LM(4))-Y(LM(3))))
        GRIDAREA=GRIDAREA+AREA
C ... AREA WITH ICE ON IT (AREA*IC/4.)
        AREA=AREA*DBLE(IC)*0.25d0
c        if(ic.ne.4 .and. ic.ne.0) then
c          print *,area,ic
c          print *,thick,height
c          print *,Q(LM(1))-DEPB(LM(1)),Q(LM(1)),DEPB(LM(1))
c          print *,Q(LM(2))-DEPB(LM(2)),Q(LM(2)),DEPB(LM(2))
c          print *,Q(LM(3))-DEPB(LM(3)),Q(LM(3)),DEPB(LM(3))
c          print *,Q(LM(4))-DEPB(LM(4)),Q(LM(4)),DEPB(LM(4))
c        endif
C ... ************************************************************
C ... ************************************************************
C ... THIS AREA IS UNCORRECTED FOR THE EFFECTS OF CURVATURE AND 
C ... WILL PROVIDE AN OVERESTIMATE OF THE VOLUME BY AS MUCH AS 15%
C ... IN MID LATITUDES. THE FOLLOWING IS A PATCH TO CORRECT THIS
C ... AND SHOULD BE REMOVED IF THE PROGRAM EVER MOVES TO 
C ... SPHERICAL COORDINATES
        if(.false.) then
          rlatmin=1d30
          rlatmax=-rlatmin
          do k=1,4
            call recpol(0.001d0*x(lm(k)),
     &                  0.001d0*y(lm(k)),rlat,rlong)
            rlatmin=min(rlatmin,rlat)
            rlatmax=max(rlatmax,rlat)
          enddo
          phi1=(90.d0-rlatmin)*radpdeg
          phi2=(90.d0-rlatmax)*radpdeg
          ratio=2*(cos(phi1)-cos(phi2))/(phi2**2-phi1**2)
          area=area*ratio
        endif
c ... ************************************************************
c ...  ************************************************************
C
C ... IF THICKNESS LT 1 METER, NEGLECT
        IF(HEIGHT.GT.1.) THEN
          RNET=RNET+AREA*ADOTAVG
          VOL=VOL+AREA*HEIGHT
          VOL1=VOL1+AREA*THICK
          AREATOT=AREATOT+AREA
          IF(THICKMAX.LT.THICK) THEN
            THICKMAX=THICK
            ITHMAX=I
          ENDIF
        ENDIF
100     CONTINUE
      ENDDO
      IF(AREATOT.GT.0.) THEN
        AVGHGT=VOL1/AREATOT
        stuff=RNET/AREATOT
      ELSE
        AVGHGT=0.d0
        stuff=0
      ENDIF
      IF(IOTOGG) THEN
        WRITE(list(ipage+1),1000) VOL*1.d-15,AREATOT*1.d-12,AVGHGT
        WRITE(list(ipage+2),1001) VOL1*1.d-15,AMASS(7),AMASS(8)
c       WRITE(list(ipage+1),1000) VOL*1.d-15/0.4,AREATOT*1.d-12,AVGHGT
c       WRITE(list(ipage+2),1001) VOL1*1.d-15/0.4,AMASS(7),AMASS(8)
        WRITE(list(ipage+3),1002) AMASS(9),ACOM,THICKMAX
        WRITE(list(ipage+4),1003) RNET*1d-15,stuff,SEALEV
        sealow=-VOL*1.d-15/0.4
        ipage=ipage+4
      ENDIF
c      if(ithmax.ne.0) then
c      write(*,*) kx(ithmax,1),q(kx(ithmax,1)),depb(kx(ithmax,1))
c      write(*,*) kx(ithmax,2),q(kx(ithmax,2)),depb(kx(ithmax,2))
c      write(*,*) kx(ithmax,3),q(kx(ithmax,3)),depb(kx(ithmax,3))
c      write(*,*) kx(ithmax,4),q(kx(ithmax,4)),depb(kx(ithmax,4))
c      endif
      if(iprt) then
        WRITE(18,2000) TIME,VOL1*1.d-15
        WRITE(17,1004) TIME,VOL1*1d-15,VOL*1d-15,AREATOT*1d-12,AMASS(7),
     &              AMASS(9),AVGHGT,HMAX,ACOM,AMASS(8)
      endif
2000  FORMAT(10X,G13.6,2X,G13.6)
1000  FORMAT(' Vflot  =',1PG13.6,' MKM**3, A    =',G13.6,
     &       ' MKM**2, <H>=',G13.6)
1001  FORMAT(' Vtot   =',1PG13.6,' MKM**3, DELSN=',G13.6,
     &       ' M    , SLSLP=',G13.6)
c1000  FORMAT(' Vflot  =',1PG13.6,' m-s.l., A    =',G13.6,
c     &       ' MKM**2, <H>=',G13.6)
c1001  FORMAT(' Vtot   =',1PG13.6,' m-s.l., DELSN=',G13.6,
c     &       ' M    , SLSLP=',G13.6)
1002  FORMAT(' TNSL   =',1PG13.6,' LAPSE RATE=',G13.6,
     &       ' MAX THICK=',G13.6)
1003  format(1X,'NET MASS BALANCE=',1p2g13.6,' SEALEV',g13.6)
1004  format(1p10g13.6)
c 1004  FORMAT(1X,1P,G13.6,0P,3F7.3,F8.0,F6.1,F7.1,F7.1,F7.3,F8.5)               
      RETURN
      END
      SUBROUTINE SCENARIO(TIME)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /EXPER/ TPERIOD,TINIT,TVSTART,TVFINAL,IEXPER
      write(*,*) ' 1-        mass(7) coef      0  10.5 k  -900 -2100'
      write(*,*) ' 2-<  70 k mass(7) coef1     0    20 k     0 -2300'
      write(*,*) '   >  70 k mass(7) coef1    70    90 k -2300     0'
      write(*,*) ' 3-        mass(7) coef      0  20.5 k     0 -2300'
      write(*,*) ' 3-        mass(8) coef    -25 -14.5 k  5600  3800'
      write(*,*) ' 4-<20.5 k mass(7) coef1     0  20.5 k     0 -2300'
      write(*,*) '   <61.5 k mass(7) coef1  20.5    41 k -2300 -1000'
      write(*,*) '   >61.5 k mass(7) coef1  61.5    82 k -2300     0'
      write(*,*) ' 5-<20.5 k mass(7) coef1     0  20.5 k -2000     0'
      write(*,*) '   <61.5 k mass(7) coef1  20.5    41 k     0 -3000'
      write(*,*) '   >61.5 k mass(7) coef1  61.5    82 k     0 -2000'
      write(*,*) ' 6-        mass(7) coef     -3   6.4 k     0 -2000'
      write(*,*) ' 7-        mass(7) coef     -3     7 k -1000 -4000'
c
      write(*,*) ' 8-        tnsl    coef  -9.55 10.55 k     0   -20'
      write(*,*) ' 9-        tnsl    coef      0     1 k    -6   -20'
      write(*,*) '10- northern hemisphere ice sheet cycle (smooth)'
      write(*,*) '   <100 k  tnsl    coef      0    20 k     0   -20'
      write(*,*) '   <120 k  tnsl    coef    100   120 k   -20     0'
      write(*,*) '   <128 k  tnsl    coef    120   128 k     0    -5'
      write(*,*) '   >128 k  tnsl = -5'
      write(*,*) '11- northern hemisphere ice sheet cycle (linear)'
      write(*,*) '   < 40 k  tnsl    coef      0    20 k    -7   -20'
      write(*,*) '   < 41 k  tnsl    coef     40    41 k   -20     0'
      write(*,*) '   < 51 k  tnsl    coef     51    52 k     0    -5'
      write(*,*) '   > 51 k  tnsl = -5'
c      write(*,*) '12         tnsl    coef      0    5.4 k    .5    -4'
c      write(*,*) '           tnsl    coef    5.4   12.3 k    -4     3'
c      write(*,*) '           tnsl = 3'
      write(*,*) '            tnsl    coef      0   10 k    -5    -6'
      write(*,*) '           tnsl    coef    10    20 k   -6    -3'
      write(*,*) '           tnsl = -3'
      write(*,*) '13- 20,000 year cycle of tnsl'
      write(*,*) '   < 10 k  tnsl    coef      0    10 k    -7   -20'
c      write(*,*) '           acom    coef      0    10 k    -5  -9.5'
c      write(*,*) '           mass(7) coef      0    10 k -2000 -2600'
      write(*,*) '   < 20 k  tnsl    coef     10    20 k   -20    -0'
c      write(*,*) '           acom    coef     10    20 k  -9.5    -5'
c      write(*,*) '           mass(7) coef      0    10 k -2600 -2000'
      write(*,*) '   > 20 k  tnsl    =    -7'
c      write(*,*) '           acom    =    -5'
      write(*,*) '           mass(7) = -2000'
      write(*,*) '14- abrupt change'
      write(*,*) '   <  5 k  tnsl    coef      0     5 k    -7   -20'
      write(*,*) '   < 15 k  tnsl = -20'
      write(*,*) '   < 20 k  tnsl    coef     15    20 k   -20     0'
      write(*,*) '   > 20 k  tnsl = 0'
      write(*,*) '15- very abrupt change'
      write(*,*) '   < 30 k  tnsl = -20'
      write(*,*) '   > 30 k  tnsl =   0'
      write(*,*) ' 16- linear deglaciation scenario'
      write(*,*) '           tnsl    coefl     0     5 k   -20     0'
      write(*,*) ' 17- read in TLIST.DATA or user defined'
      write(*,*) ' 18- 20,000 year cycle of tnsl for superglacial'
      write(*,*) '   < 10 k  tnsl    coef      0    10 k    -7   -23'
c      write(*,*) '           acom    coef      0    10 k    -5  -9.5'
c      write(*,*) '           mass(7) coef      0    10 k -2000 -2600'
      write(*,*) '   < 20 k  tnsl    coef     10    20 k   -23    -0'
c      write(*,*) '           acom    coef     10    20 k  -9.5    -5'
c      write(*,*) '           mass(7) coef      0    10 k -2600 -2000'
      write(*,*) '   > 20 k  tnsl    =    -0'
c      write(*,*) '           acom    =    -5'
c      write(*,*) '           mass(7) = -2000'
      write(*,*) ' 19- Abrupt climat change companion to scenario 18'
      write(*,*) '   <  1 k tnsl =   0'
      write(*,*) '   >  1 k tnsl = -23'
      write(*,*) '   > 15 k tnsl =   0'
      write(*,*) ' 20- Anatarctic scenario'
      write(*,*) '   < 10 k  tnsl    coef      0    10 k   -14   -23'
      write(*,*) '   < 20 k  tnsl    coef     10    20 k   -23   -14'
      write(*,*) '   > 20 k  tnsl    =    -0'
      write(*,*) ' 21- scandinavian scenario'
      write(*,*) '     40 k period tnsl   coef 0    20 k   -2   -18'
      write(*,*) ' 22- antarctic scenario'
      write(*,*) '     20 k period tnsl   coef 0    5 k   -20   -8'
      write(*,*) ' 23- user defined scenario'
      write(*,*) '   input period,tstart,tfinal '
      write(*,*) ' 24- eismint scenario'
      
      READ(*,*) IEXPER
      write(99,*) IEXPER
      IF(IEXPER.EQ.23) THEN
        TINIT=TIME
        WRITE(*,*) 'INPUT PERIOD,INITIAL TNSL, AND FINAL TNSL'
        READ(*,*) TPERIOD,TVSTART,TVFINAL
        WRITE(99,*) TPERIOD,TVSTART,TVFINAL
      ENDIF
      END
      FUNCTION AFUNCT(TIME,IDOT,AMASS,HT,BD,PSURF,SLOPN,X,Y,TS)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
      COMMON /SNOW/ SNOLIN,SNO(20)
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      COMMON /TTNSL/ TTIME(NSMAX),TLIST(NSMAX),NTNSL
      COMMON /EXPER/ TPERIOD,TINIT,TVSTART,TVFINAL,IEXPER
      DIMENSION AMASS(11),SLOPN(4)
      data pi /3.14159265358979d0/
      FPC(S,E)=(-7.044897460d+02)*S**2+( 3.286502080d+01)*S+
     1  (-8.363957700d-08)*E**2+( 2.803513780d-04)*E+
     1  (-8.280271290d-02)
      FPX(S,E)=(-1.011654790d+03)*S**2+( 5.757759090d+01)*S+
     1  (-1.945462600d-07)*E**2+( 6.485106420d-04)*E+
     1  (-1.714566350d-01)
      FSX(S,E)=( 1.613718720d+02)*S**2+(-3.064633370d+00)*S+
     1  (-2.274655910d-07)*E**2+( 9.221490470d-04)*E+
     1  (-3.658012150d-01)
      FTM(S,E)=(-2.475505130d+03)*S**2+( 7.388871770d+01)*S+
     1  (-7.129281700d-07)*E**2+( 3.180227710d-03)*E+
     1  (-6.778311130d-01)
      FSM(S,E)=(-1.272729690d+04)*S**2+( 1.028782810d+02)*S+
     1  ( 1.790857030d-07)*E**2+(-1.209567770d-03)*E+
     1  ( 2.244830130d+00)
C CHOOSE COMPENENT OF THE SLOPE, (1) - SQRT(DHDX**2+DHDY**2)
C                                (2) - DHDX
C                                (3) - DHDY
C                                (4) - DHDX*SIN(THETA)+DHDY*COS(THETA)
      SL=SLOPN(1)
C
C FOLLOWING SECTION IS QUITE VARIABLE. IT CONTAINS CALLS TO COEF AND
C COEF1 TO VARY VARIOUS MASS BALANCE PARAMETERS IN A SINUSOIDAL FASHION,
C AND WILL BE DIFFERENT FOR DIFFERENT EXPERIMENTS.
C
C     SNOLIN=0.
c     
      TNSL=AMASS(9)
      IF(IEXPER.EQ.1) THEN
C
        AMASS(7)=COEF(TIME,0.,10500.,-900.,-2100.)
C
      ELSEIF(IEXPER.EQ.2) THEN
C
        IF(TIME.LT.70000.) THEN
          AMASS(7)=COEF1(TIME,0.,20000.,0.,-2300.)
        ELSE
          AMASS(7)=COEF1(TIME,70000.,90000.,-2300.,0.)
        ENDIF
C
      ELSEIF(IEXPER.EQ.3) THEN
C
        DELSNNP=COEF(TIME,0.,20500.,0.,-2300.)
        AMASS(7)=DELSNNP
        DELSNEQ=COEF(TIME,-25000.,-14500.,5600.,3800.)
        AMASS(8)=(DELSNEQ-DELSNNP)*1d-7
C
      ELSEIF(IEXPER.EQ.4) THEN
C
        IF(TIME.LE.20500.) THEN
          AMASS(7)=COEF1(TIME,0.,20500.,0.,-2300.)
        ELSEIF(TIME.LE.61500.) THEN
          AMASS(7)=COEF(TIME,20500.,41000.,-2300.,-1000.)
        ELSE
          AMASS(7)=COEF1(TIME,61500.,82000.,-2300.,0.)
        ENDIF
C
      ELSEIF(IEXPER.EQ.5) THEN
C
        IF(TIME.LE.20500.) THEN
          AMASS(7)=COEF1(TIME,0.,20500.,-2000.,0.)
        ELSEIF(TIME.LE.61500.) THEN
          AMASS(7)=COEF(TIME,20500.,41000.,0.,-3000.)
        ELSE
          AMASS(7)=COEF1(TIME,61500.,82000.,0.,-2000.)
        ENDIF
C
      ELSEIF(IEXPER.EQ.6) THEN
C
        AMASS(7)=COEF(TIME,-3600.,6400.,0.,-2000.)
C
      ELSEIF(IEXPER.EQ.7) THEN
C
        AMASS(7)=COEF(TIME,-3000.,7000.,-1000.,-4000.)
C
c variable stuff, to modify TNSL
      ELSEIF(IEXPER.EQ.8) THEN
        TNSL=COEF(TIME,-9550.,10950.,0.,-20.)

      ELSEIF(IEXPER.EQ.9) THEN
        TNSL=COEF(TIME,0.,1000.,-6.,-20.)

      ELSEIF(IEXPER.EQ.10) THEN
C EXPERIMENT FOR NORTHERN HEMISHHERE ICE SHEET CYCLE (SMOOTH VARIATION)
        IF(TIME.LT.100000.) THEN
          TNSL=COEF(TIME,0.,20000.,-7.,-20.)
        ELSEIF(TIME.LT.120000.) THEN
          TNSL=COEF1(TIME,100000.,120000.,-20.,0.)
        ELSEIF(TIME.LT.128000.) THEN
          TNSL=COEF1(TIME,120000.,128000.,0.,-5.)
        ELSE
          TNSL=-5.d0
        ENDIF
C END OF EXPERIMENT FOR NORTHERN HEMISPHERE ICE SHEET CYCLE

      ELSEIF(IEXPER.EQ.11) THEN
C EXPERIMENT FOR NORTHERN HEMISHHERE ICE SHEET CYCLE (LINEAR VARIATION)
        IF(TIME.LT.40000.) THEN
          TNSL=COEFL(TIME,0.,20000.,-7.,-20.)
        ELSEIF(TIME.LT.41000.) THEN
          TNSL=COEFL(TIME,40000.,41000.,-20.,0.)
        ELSEIF(TIME.LT.51000.) THEN
          TNSL=COEFL(TIME,51000.,52000.,0.,-5.)
        ELSE
          TNSL=-5.d0
        ENDIF
C END OF EXPERIMENT FOR NORTHERN HEMISPHERE ICE SHEET CYCLE

      ELSEIF(IEXPER.EQ.12) THEN
C lake baikal experiment tnsl goes from 0 to -6 and back in 20000
c        IF(TIME.LT.5400.) THEN
c          TNSL=COEF(TIME,0.,5400.,.5,-4.)
c        ELSEIF(TIME.LT.12300.) THEN
c          TNSL=COEF(TIME,5400.,12300.,-4.,3.)
c        ELSE
c          TNSL=3.
c        ENDIF
        IF(TIME.LT.10000.) THEN
          TNSL=COEF(TIME,0.,10000.,-5.,-6.)
        ELSEIF(TIME.LT.20000.) THEN
          TNSL=COEF(TIME,10000.,20000.,-6.,-3.)
        ELSE
          TNSL=-3.d0
        ENDIF
c end lake baikal cycle

      ELSEIF(IEXPER.EQ.13) THEN
C smooth variation experiment, 20,000 yr sin cycle both tnsl and lapse
C rate     
        IF(TIME.LT.10000.) THEN
          TNSL=COEF(TIME,0.,10000.,-7.,-20.)
c          ACOM=COEF(TIME,0.,10000.,-5.,-9.5)
c          AMASS(7)=COEF(TIME,0.,10000.,-2000.,-2600.)
        ELSEIF(TIME.LT.20000.) THEN
          TNSL=COEF(TIME,10000.,20000.,-20.,0.)
c          ACOM=COEF(TIME,10000.,20000.,-9.5,-5.)
c          AMASS(7)=COEF(TIME,10000.,20000,-2600.,-2000.)
        ELSE
          TNSL=0.d0
c          ACOM=-5.d0
c          AMASS(7)=-2000.d0
        ENDIF
c end smooth cycle

      ELSEIF(IEXPER.EQ.14) THEN
c abrupt change experiment
        IF(TIME.LT.5000.) THEN
          TNSL=COEFL(TIME,0.,5000.,-7.,-20.)
        ELSEIF(TIME.LT.15000.) THEN
          TNSL=-20.d0
        ELSEIF(TIME.LT.20000.) THEN
          TNSL=COEFL(TIME,15000.,20000.,-20.,0.)
        ELSE
          TNSL=0.d0
        ENDIF
c end abrupt change experiment

      ELSEIF(IEXPER.EQ.15) THEN
c very abrupt climat change
        IF(TIME.LE.30001.) THEN
          TNSL=-20.d0
        ELSE
          TNSL=0.d0
        ENDIF
c end very abrupt clmate change

      ELSEIF(IEXPER.EQ.16) THEN
c begin linear deglaciation scenario
        TNSL=COEFL(TIME,0.,5000.,-20.,0.)
c end linear deglaciation scenario
      ELSEIF(IEXPER.EQ.17) THEN
c use read in list of tnsl
        IF(NTNSL.GT.0) TNSL=FINDTNSL(TIME)
        AMASS(9)=TNSL
      ELSEIF(IEXPER.EQ.18) THEN
C smooth variation experiment, 20,000 yr sin cycle both tnsl and lapse
C rate     
        IF(TIME.LT.10000.) THEN
          TNSL=COEF(TIME,0.,10000.,-7.,-23.)
c          ACOM=COEF(TIME,0.,10000.,-5.,-9.5)
c          AMASS(7)=COEF(TIME,0.,10000.,-2000.,-2600.)
        ELSEIF(TIME.LT.20000.) THEN
          TNSL=COEF(TIME,10000.,20000.,-23.,-0.)
c          ACOM=COEF(TIME,10000.,20000.,-9.5,-5.)
c          AMASS(7)=COEF(TIME,10000.,20000,-2600.,-2000.)
        ELSE
          TNSL=-0.d0
          ACOM=-5.d0
          AMASS(7)=-2000.d0
        ENDIF
c end smooth cycle

      ELSEIF(IEXPER.EQ.19) THEN
c abrupt climate change companion to scenario 18
        IF(TIME.LT.1000.) THEN
          TNSL=0.d0
c          ACOM=-5
c          AMASS(7)=-2000.
        ELSEIF(TIME.LT.15000.) THEN
          TNSL=-23.d0
c          ACOM=-9.5
c          AMASS(7)=-2600.
        ELSE
          TNSL=0.d0
c          ACOM=-5.
c          AMASS(7)=-2000.
        ENDIF
c end abrupt climate change

      ELSEIF(IEXPER.EQ.20) THEN
C antarrctic variation experiment, 20,000 yr sin cycle tnsl
C      
        IF(TIME.LT.10000.) THEN
          TNSL=COEF(TIME,0.,10000.,-14.,-23.)
        ELSEIF(TIME.LT.20000.) THEN
          TNSL=COEF(TIME,10000.,20000.,-23.,-14.)
        ELSE
          TNSL=-14.d0
        ENDIF
c end smooth cycle
      ELSEIF(IEXPER.EQ.21) THEN
C scandinavian variation experiment, 20,000 yr sin cycle tnsl
C      
        TNSL=COEF(TIME,0.,20000.,-2.,-18.)
c end smooth cycle
      ELSEIF(IEXPER.EQ.22) THEN
C antarctica variation experiment, 20,000 yr sin cycle tnsl
C      
        TNSL=COEF(TIME,-5000.,5000.,-20.,-8.)
c end smooth cycle
      ELSEIF(IEXPER.EQ.23) THEN
C user defined variation experiment.
C      TPERIOD yr sin cycle tnsl
C      TVSTART  initial time
C      TVFINAL final time
C      
        TNSL=COEF(TIME,REAL(TINIT),REAL(TINIT+TPERIOD),
     &                 REAL(TVSTART),REAL(TVFINAL))
        AMASS(9)=TNSL
c end smooth cycle
      ELSEIF(IEXPER.EQ.24) THEN
C eismint variation experiment, 20,000,40000 yr sin cycle acc
c this is also zone 24
C 
c following for uniform, constant accumulation rate
c          afunct=0.3     
c following for uniform accumulation 20000 yr cycle
c          afunct=0.3+0.2*sin(8.*atan(1.)*time/20000.)
c following for uniform accumulation 40000 yr cycle
c          afunct=0.3+0.2*sin(8.*atan(1.)*time/40000.)
c following for ablation zone experiment, steady state
c          distance in km
           dist=sqrt(x**2+y**2)*.001d0
c equilibrium line distance nonuniform accumulation with ablation zone
c following for constant eld
c           eld=450.d0
c following for eld 20000 yr cycle
          eld=450.d0+100.d0*sin(8.d0*atan(1.d0)*time/20000.d0)
c following for eld 40000 yr cycle
c          eld=450.d0+100.d0*sin(8.d0*atan(1.d0)*time/40000.d0)
c
           afunct=min(0.5d0,1d-2*(eld-dist))
c
          return
c end eismint cycle
      ENDIF

        
c end variable stuff of TNSL
C END VARIABLE SECTION
C
C BRANCH FOR VARIOUS MASS BALANCE RELATIONSHIPS
      IF(IDOT.EQ.0) THEN
        AFUNCT=0.0d0
        RETURN
      ENDIF
      GOTO(10,20,30,40,50,60,70,                     ! IDOT 1..7 straight lines
     &     200,210,220,230,240,                      ! IDOT 8..12 
     &     300,310,320,330,340,345,347,348,350,360,  ! IDOT 13..22
     &     370,380,                                  ! IDOT 23..24
     &     400,400,400,400,400,400,400),IDOT         ! IDOT 25..31 data-based
C
10    CONTINUE
c ... mars peak and pole....
c         AMARG=AMASS(1)
c         APEAK=AMASS(2)
c         ADOME=AMASS(3)
c         HE=AMASS(4)
c         HP=AMASS(5)
c         HM=AMASS(6)
          X0=AMASS(1)
          Y0=AMASS(2)
          R0=AMASS(3)
          HP=AMASS(4)
          SL=-HP/R0
          AFUNCT=HP+SL*SQRT((X*0.001-X0)**2+(Y*0.001-Y0)**2)
        return
      GOTO 100
C
20    CONTINUE
c ... mars peak and pole....
          AMARG=AMASS(1)
          APEAK=AMASS(2)
          ADOME=AMASS(3)
          HE=AMASS(4)
          HP=AMASS(5)
          HM=AMASS(6)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
30    CONTINUE
C **** =3 POLAR - B
          AMARG=0.0d0
          APEAK=1.d0
          ADOME=.10d0
          HE=0.d0
          HP=1250.d0
          HM=2500.d0
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
40    CONTINUE
C **** =4 POLAR - C
          AMARG=-1.0d0
          APEAK=1.5d0
          ADOME=.1d0
          HE=100.d0
          HP=1250.d0
          HM=2500.d0
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
50    CONTINUE
C **** IDOT=5 MARITIME - A
          AMARG=-2.0d0
          APEAK=2.d0
          ADOME=.15d0
          HE=1500.d0
          HP=1750.d0
          HM=2000.d0
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
60    CONTINUE
C **** IDOT=6 MARITIME - B
          AMARG=-2.0d0
          APEAK=2.d0
          ADOME=.15d0
          HE=1000.d0
          HP=1500.d0
          HM=2000.d0
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
70    CONTINUE
C **** IDOT=7 MARITIME - C
          AMARG=-3.0d0
          APEAK=2.d0
          ADOME=.15d0
          HE=1000.d0
          HP=1500.d0
          HM=2000.d0
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 100
C
100   CONTINUE
      AMASS(1)=AMARG
      AMASS(2)=APEAK
      AMASS(3)=ADOME
      AMASS(4)=HE
      AMASS(5)=HP
      AMASS(6)=HM
      HTT=HT-SNOLIN
      IF(HTT.GT.HM) THEN
        AFUNCT=ADOME
      ELSE
        IF(HTT.LE.0.) THEN
          AFUNCT=AMARG
        ELSE
          IF(HTT.GT.0. .AND. HTT.LE.HE) THEN
            A1=AMARG
            H1=0.d0
            A2=0.d0
            H2=HE
          ENDIF
          IF(HTT.GT.HE .AND. HTT.LE.HP) THEN
            A1=0.d0
            H1=HE
            A2=APEAK
            H2=HP
          ENDIF
          IF(HTT.GT.HP .AND. HTT.LE.HM) THEN
            A1=APEAK
            H1=HP
            A2=ADOME
            H2=HM
          ENDIF
          AFUNCT=(A1*(H2-HTT)-A2*(H1-HTT))/(H2-H1)
        ENDIF
      ENDIF
      IF(BD.LE.-5000.) AFUNCT=-5.d0
      RETURN
C
200   CONTINUE
C **** IDOT=8 SM (SUBPOLAR MARITIME) SLOPE DEPENDENT
          ADOME=.10d0
          SMAX=.01d0
          B0=-5.d0
          HEQ=700.d0
          B0DOT=B0/HEQ
          BNEQ=FSM(SMAX,HEQ)
          HADD=-BNEQ/B0DOT
          IF(HT.LE.HEQ) THEN
            BN=B0-B0DOT*HT
          ELSE IF(HT.LE.(HEQ+HADD)) THEN
            BN=FSM(SL,HT)+B0DOT*(HADD-HT+HEQ)
            IF(BN.LE.0.) BN=ADOME
          ELSE
            BN=FSM(SL,HT)
            IF(BN.LE.0.) BN=ADOME
          ENDIF
          AFUNCT=BN
      GOTO 1000
C
210   CONTINUE
C **** IDOT=9 TM (TEMPERATE MARITIME) SLOPE DEPENDENT
          ADOME=.10d0
          SMAX=.03d0
          B0=-8.d0
          HEQ=1125.d0
          B0DOT=B0/HEQ
          BNEQ=FTM(SMAX,HEQ)
          HADD=-BNEQ/B0DOT
          IF(HT.LE.HEQ) THEN
            BN=B0-B0DOT*HT
          ELSE IF(HT.LE.(HEQ+HADD)) THEN
            BN=FTM(SL,HT)+B0DOT*(HADD-HT+HEQ)
            IF(BN.LE.0.) BN=ADOME
          ELSE
            BN=FTM(SL,HT)
            IF(BN.LE.0.) BN=ADOME
          ENDIF
          AFUNCT=BN
      GOTO 1000
C
220   CONTINUE
C **** IDOT=10 SX (SUBPOLAR MIX) SLOPE DEPENDENT
          ADOME=.10d0
          SMAX=.07d0
          B0=-5.d0
          HEQ=1250.d0
          B0DOT=B0/HEQ
          BNEQ=FSX(SMAX,HEQ)
          HADD=-BNEQ/B0DOT
          IF(HT.LE.HEQ) THEN
            BN=B0-B0DOT*HT
          ELSE IF(HT.LE.(HEQ+HADD)) THEN
            BN=FSX(SL,HT)+B0DOT*(HADD-HT+HEQ)
            IF(BN.LE.0.) BN=ADOME
          ELSE
            BN=FSX(SL,HT)
            IF(BN.LE.0.) BN=ADOME
          ENDIF
          AFUNCT=BN
      GOTO 1000
C
230   CONTINUE
C **** IDOT=11 PX (POLAR MIX) SLOPE DEPENDENT
          ADOME=.05d0
          SMAX=.04d0
          B0=-2.5d0
          HEQ=600.d0
          B0DOT=B0/HEQ
          BNEQ=FPX(SMAX,HEQ)
          HADD=-BNEQ/B0DOT
          IF(HT.LE.HEQ) THEN
            BN=B0-B0DOT*HT
          ELSE IF(HT.LE.(HEQ+HADD)) THEN
            BN=FPX(SL,HT)+B0DOT*(HADD-HT+HEQ)
            IF(BN.LE.0.) BN=ADOME
          ELSE
            BN=FPX(SL,HT)
            IF(BN.LE.0.) BN=ADOME
          ENDIF
          AFUNCT=BN
      GOTO 1000
C
240   CONTINUE
C **** IDOT=12 PC (POLAR CONTINENTAL) SLOPE DEPENDENT
          ADOME=.05d0
          SMAX=.01d0
          B0=-1.d0
          HEQ=300.d0
          B0DOT=B0/HEQ
          BNEQ=FPC(SMAX,HEQ)
          HADD=-BNEQ/B0DOT
          IF(HT.LE.HEQ) THEN
            BN=B0-B0DOT*HT
          ELSE IF(HT.LE.(HEQ+HADD)) THEN
            BN=FPC(SL,HT)+B0DOT*(HADD-HT+HEQ)
            IF(BN.LE.0.) BN=ADOME
          ELSE
            BN=FPC(SL,HT)
            IF(BN.LE.0.) BN=ADOME
          ENDIF
          AFUNCT=BN
      GOTO 1000
C
300   CONTINUE
C **** IDOT=13 SM EXPONENTIAL
          A1=-7.57655d0
          A2=2.56891d0
          RL1=2.773327d-6
          RL2=2.916233d-7
          SNOLIN=SNO(13)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 900
C
310   CONTINUE
C **** IDOT=14 TM EXPONENTIAL
          A1=-11.6877d0
          A2=3.67407d0
          RL1=1.083791d-6
          RL2=3.641119d-8
          SNOLIN=SNO(14)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 900
C
320   CONTINUE
C **** IDOT=15 SX EXPONENTIAL
          A1=-5.685629d0
          A2=0.857527d0
          RL1=1.322827d-6
          RL2=9.461793d-8
          SNOLIN=SNO(15)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 900
C
330   CONTINUE
C **** IDOT=16 PX EXPONENTIAL
          A1=-2.84904d0
          A2=.837780d0
          RL1=8.534590d-6
          RL2=3.405352d-8
          SNOLIN=SNO(16)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 900
C
340   CONTINUE
C **** IDOT=17 PC EXPONENTIAL
          A1=-1.29363d0
          A2=.299798d0
          RL1=1.502293d-5
          RL2=3.661944d-8
          SNOLIN=SNO(17)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 900
C
345   CONTINUE
C **** IDOT=18 ANTARCTIC VALUUES
C **** FOLLOWING IS CATCHALL FOR DIFFERENT MASS BALANCE RELATIONSHIPS
C **** AND MAY HAVE DIFFERENT VALUES FOR DIFFERENT EXPERIMENTS
C **** DOUBLE OF 17
          A1=-2.d0*1.29363d0
          A2=2.d0*.299798d0
          RL1=1.502293d-5
          RL2=3.661944d-8
C **** CTRL VALUES
C         A1=-.838547d0
C         A2=.887967d0
C         RL1=3.471435d-7
C         RL2=3.918019d-8
C **** REGSST VALUES
C         A1=-1.22491d0
C         A2=1.21169d0
C         RL1=3.736239d-7
C         RL2=8.897420d-8
C **** WARM MIN VALUES
C         A1=-15.8931d0
C         A2=15.3500d0
C         RL1=3.919554d-7
C         RL2=2.362680d-7
          SNOLIN=SNO(18)
          SNOLIN=AMASS(7)+AMASS(8)*SQRT(X**2+Y**2)
      GOTO 900
C
347   CONTINUE
C **** IDOT=19 GROSSWALD'S DIPPING SNOWLINE SCHEME
      DIST=SQRT(X**2+Y**2)
C FOLLOWING FOR 100M/100KM WITH SNOWLINE AT -1000 M AT POLE
      SNSLOPE=.001d0
C FOLLOWING IS FOR SIBERIAN SNOWLINE DATA, 2E-03
      SNSLOPE=2.d-03
      SNSLOPE=AMASS(8)
      SNZERO=AMASS(7)
C FOLLOWING FOR -1 M AT S.L., +1 M AT 2000 M ELEVATION, OR 2M/2000M
      BALGRAD=.001d0
      BALGRAD=AMASS(9)
      SNLINE=SNZERO+SNSLOPE*DIST
      ACC=(HT-SNLINE)*BALGRAD
      IF(ACC.GT.2.) THEN
        ACC=ACC-2.d0*(ACC-2.d0)
        IF(ACC.LT..05) ACC=.05d0
      ELSEIF(ACC.LT.-2.) THEN
        ACC=-2.d0
      ENDIF
      AFUNCT=ACC
      RETURN
348   CONTINUE
C **** IDOT=20 MIKE'S METEROLOGICAL MODEL

      AMASS(9)=TNSL
      AFUNCT=ACCUM(IDOT,X*.001D0,Y*.001D0,HT,SL,PSURF,TNSL,TS)
      if(.true. .and. afunct.gt.0) then
        afunct=afunct*10
      else
        afunct=afunct*1
      endif
      RETURN
C
350   CONTINUE
C **** IDOT=21 WIND/SLOPE METHOD
C WIND/SLOPE METHODA
C FOR NOW JUST AN AUGMENTED METEROLOGICAL ZONE
c     TNSL=AMASS(9)
C     TNSL=COEF(TIME,-300.,800.,-6.,-15.)
      AMASS(9)=TNSL
c     threshangle=25d0
c     threshangle=cos(threshangle*pi/180d0)
      threshangle=0.80d0
      costheta=slopn(4)/slopn(1)
c     print *,180*costheta/pi,slopn(4),slopn(1)
      if(.true.) then
        AFUNCT=ACCUM(IDOT,X*.001D0,Y*.001D0,HT,SL,PSURF,TNSL,TS)
c    &              +SLOPN(4)*WINDIR(2)+0.001
        acc=afunct+abl
        costheta=0.5*(costheta+1.0)
        acc=acc*costheta
        afunct=acc-abl
      elseif(costheta.gt.threshangle) then
        AFUNCT=ACCUM(IDOT,X*.001D0,Y*.001D0,HT,SL,PSURF,TNSL,TS)+
c    &               SLOPN(4)*WINDIR(2)
     &               0
      else
        AFUNCT=ACCUM(IDOT,X*.001D0,Y*.001D0,HT,SL,PSURF,TNSL,TS)
        afunct=-abl
c       afunct=-0.0003
      endif
      if(.true. .and. afunct.gt.0) then
        afunct=afunct*10
      else
        afunct=afunct*1
      endif
      RETURN
C
360   CONTINUE
C **** IDOT=22 HUYBRECHTS METHOD
c     TNSL=AMASS(9)
C     TNSL=COEF(TIME,-300.,800.,-6.,-15.)
      AMASS(9)=TNSL
      AFUNCT=ACCUMH(X*.001D0,Y*.001D0,HT,TNSL,TS)
      RETURN
370   CONTINUE
C **** IDOT=23 WHITE HOLE METHOD (NO ABLATION ...)
c     TNSL=AMASS(9)
C     TNSL=COEF(TIME,-300.,800.,-6.,-15.)
      AMASS(9)=TNSL
      AFUNCT=WOABLATE(X*.001D0,Y*.001D0,HT,SL,0.D0,TNSL,TS)
      RETURN
380   CONTINUE
C **** IDOT=24 MILANKOVICH METHOD
      AMASS(9)=TNSL
      call milank(TIME,X*.001D0,Y*.001D0,HT,SL,0.D0,TNSL,TS,accm)
      AFUNCT=accm
      RETURN
C
400   CONTINUE
C **** IDOT=25 ISLSCP-DATA-BASED MASS BALANCE
C **** IDOT=26 NCEP2-DATA-BASED MASS BALANCE
C **** IDOT=27 GCM-RESULT-BASED MASS BALANCE
C **** IDOT=28 SEALEVEL-DEPENDENT INTERPOLATION 
C      OF NCEP2-DATA-BASED AND GCM-RESULT-BASED MASS BALANCE
C **** IDOT=29 NCEP2-DATA-BASED MASS BALANCE WITH TOPO MOD
C **** IDOT=30 GCM-RESULT-BASED MASS BALANCE WITH TOPO MOD
C **** IDOT=31 SEALEVEL-DEPENDENT INTERPOLATION 
C      OF NCEP2-DATA-BASED AND GCM-RESULT-BASED MASS BALANCE WITH TOPO MOD
C **** 
      AMASS(9)=TNSL
      AFUNCT=ACCUM(IDOT,X*.001D0,Y*.001D0,HT,SL,PSURF,TNSL,TS)
      RETURN
900   CONTINUE
C
C AGAIN A VARIABLE SECTION FOR DIFFERENT EXPERIMENTS, HERE VARYING
C SPECIFIC PARAMETERS IN THE DIFFERENT MASS BALANCE SCHEMES. THIS CAN BE
C USED TO SMOOTHLY SWITCH FROM ONE ZONE TO ANOTHER.
C FOLLOWING FOR TIME VARYING ACCUM CURVES
C COEF IS FUNCTION COEF(TIME,TINIT,TFINAL,VINIT,VFINAL)
C COEF DOES SINUSOIDAL VARIATION REGARDLES OF TIME
C COEF1 FOR TIME<TINIT = VINIT, TIME>TFINAL = VFINAL, IN BETWEEN IS
C COSINE INTERPOLATION
C FOLLOWING ARE INITIAL TIME AND FINAL TIME (HALF CYCLE)
C     TINIT=0.
C     TFINAL=5000.
C FOLLOWING 4 FOR CONSTANT
      A11=A1
      A22=A2
      RL11=RL1
      RL22=RL2
C FOLLOWING 4 FOR TIME VARYING SINUSOIDAL
C     A11=COEF(TIME,TINIT,TFINAL,A1,A1)
C     A22=COEF(TIME,TINIT,TFINAL,A2,2.*A2)
C     RL11=COEF(TIME,TINIT,TFINAL,RL1,RL1*.01)
C     RL22=COEF(TIME,TINIT,TFINAL,RL2,RL2)
C FOLLOWING 4 FOR COSINE INTERPOLATED STEP
C BEGINNING AND END OF TIME STEP
C     TINIT=1000.
C     TFINAL=22000.
C     A11=COEF1(TIME,TINIT,TFINAL,A1,-14.074)
C     A22=COEF1(TIME,TINIT,TFINAL,A2,4.83112)
C     RL11=COEF1(TIME,TINIT,TFINAL,RL1,1.308677E-06)
C     RL22=COEF1(TIME,TINIT,TFINAL,RL2,1.61273E-7)
C     A11=COEF1(TIME,TINIT,TFINAL,A1,A1*18.95)
C     A22=COEF1(TIME,TINIT,TFINAL,A2,17.287*A2)
C     RL11=COEF1(TIME,TINIT,TFINAL,RL1,RL1*1.129)
C     RL22=COEF1(TIME,TINIT,TFINAL,RL2,6.03*RL2)
C     A11=COEF1(TIME,TINIT,TFINAL,A1,A1*1.461)
C     A22=COEF1(TIME,TINIT,TFINAL,A2,1.365*A2)
C     RL11=COEF1(TIME,TINIT,TFINAL,RL1,RL1*1.076)
C     RL22=COEF1(TIME,TINIT,TFINAL,RL2,2.271*RL2)
C END OF VARIABLE SECTION
C
      HTT=HT-SNOLIN
      IF(HTT.LT.0.) HTT=0.d0
      ARG=-RL11*HTT**2
      IF(ARG.GT.-180.) THEN
        TERM1=A11*EXP(ARG)
      ELSE
        TERM1=0.d0
      ENDIF
      ARG=-RL22*HTT**2
      IF(ARG.GT.-180.) THEN
        TERM2=A22*EXP(ARG)
      ELSE
        TERM2=0.d0
      ENDIF
      AFUNCT=TERM1+TERM2
1000  CONTINUE
      RETURN
      END
      FUNCTION COEF(TIME, TINIT, TFINAL, VINIT, VFINAL)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL TINIT, TFINAL, VINIT, VFINAL
C STARTING AT TIME=TINIT SINUSOIDALLY VARY FROM VINIT TO VFINAL WITH
C PERIOD = TFINAL-TINIT. CONTINUE PERIODICITY.
      IF(TIME.LT.TINIT) THEN
        COEF=DBLE(VINIT)
      ELSE
        ZEROL=0.5d0*DBLE(VINIT+VFINAL)
        AMPL=0.5d0*DBLE(VINIT-VFINAL)
        PHASE=(TIME-DBLE(TINIT))/DBLE(TFINAL-TINIT)*3.1415927d0
        COEF=ZEROL+AMPL*COS(PHASE)
      ENDIF
      END
C
      FUNCTION COEF1(TIME, TINIT, TFINAL, VINIT, VFINAL)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL TINIT, TFINAL, VINIT, VFINAL
C STARTING AT TIME TINIT VARY FROM VINIT TO VFINAL SUNUSOIDALLY. FOR
C TIME < TINIT VALUE=VINIT, FOR TIME > VFINAL VALUE= VFINAL.
      IF(TIME.LE.TINIT) THEN
        COEF1=DBLE(VINIT)
      ELSE IF(TIME.GE.TFINAL) THEN
        COEF1=dble(VFINAL)
      ELSE
        ZEROL=0.5d0*dble(VINIT+VFINAL)
        AMPL=0.5d0*dble(VINIT-VFINAL)
        PHASE=(TIME-dble(TINIT))/dble(TFINAL-TINIT)*3.1415927d0
        COEF1=ZEROL+AMPL*COS(PHASE)
      ENDIF
      END
      FUNCTION COEFL(TIME, TINIT, TFINAL, VINIT, VFINAL)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL TINIT, TFINAL, VINIT, VFINAL
C STARTING AT TIME TINIT VARY FROM VINIT TO VFINAL LINEARLY. FOR
C TIME < TINIT VALUE=VINIT, FOR TIME > VFINAL VALUE= VFINAL.
      IF(TIME.LE.TINIT) THEN
        COEFL=dble(VINIT)
      ELSE IF(TIME.GE.TFINAL) THEN
        COEFL=dble(VFINAL)
      ELSE
        COEFL=dble(VINIT)+(TIME-dble(TINIT))*dble(VFINAL-VINIT)/
     &        dble(TFINAL-TINIT)
      ENDIF
      END
C function based on LOCATE in RECIPES, page 90
      FUNCTION FINDTNSL(TIME)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
      COMMON /TTNSL/ TTIME(NSMAX),TLIST(NSMAX),NTNSL
      DATA IPASS /0/
      SAVE IPASS,TSAVE,TNSLSAVE
      IF(IPASS.EQ.0) THEN
        TSAVE=TIME
        IPASS=1
      ELSEIF(TSAVE.EQ.TIME) THEN
C ..... TIME HASNT CHANGED, USE LAST ONE AND RETURN ...
        FINDTNSL=TNSLSAVE
        RETURN
      ENDIF
      TSAVE=TIME
C if time greater than last ttime, tnsl=last tlist
      IF(TIME.GT.TTIME(NTNSL)) THEN
        TNSLSAVE=TLIST(NTNSL)
        FINDTNSL=TNSLSAVE
        RETURN
      ENDIF
C if time less than first ttime, tnsl=first tlist
      IF(TIME.LT.TTIME(1)) THEN
        TNSLSAVE=TLIST(1)
        FINDTNSL=TNSLSAVE
        RETURN
      ENDIF
c if time = one of the ttimes, then tnsl= that tlist
      DO I=1,NTNSL
        IF(TIME.EQ.TTIME(I)) THEN
          TNSLSAVE=TLIST(I)
          FINDTNSL=TNSLSAVE
          RETURN
        ENDIF
      ENDDO
c interpolate
      N=NTNSL
      JL=0
      JU=N+1
10    CONTINUE
        IF(JU-JL.GT.1) THEN
          JM=(JU+JL)/2
          IF((TTIME(N).GT.TTIME(1)).EQV.(TIME.GT.TTIME(JM))) THEN
            JL=JM
          ELSE
            JU=JM
          ENDIF
          GOTO 10
        ENDIF
      J=JL
      SLOPE=(TLIST(J+1)-TLIST(J))/(TTIME(J+1)-TTIME(J))
      TNSL=TLIST(J)+SLOPE*(TIME-TTIME(J))
      TNSLSAVE=TNSL
      FINDTNSL=TNSLSAVE
      END
c=============================================================
      FUNCTION FINDSEAL(TIME)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      DATA IPASS /0/
      SAVE IPASS,TSAVE,SLSAVE
      IF(IPASS.EQ.0) THEN
        TSAVE=TIME
        IPASS=1
      ELSEIF(TSAVE.EQ.TIME) THEN
C ..... TIME HASNT CHANGED, USE LAST ONE AND RETURN ...
        FINDSEAL=SLSAVE
        RETURN
      ENDIF
      TSAVE=TIME
C if time greater than last ttime, tnsl=last tlist
      IF(TIME.GT.STIME(NSEAL)) THEN
        SLSAVE=SLIST(NSEAL)
        FINDSEAL=SLSAVE
        RETURN
      ENDIF
C if time less than first ttime, tnsl=first tlist
      IF(TIME.LT.STIME(1)) THEN
        SLSAVE=SLIST(1)
        FINDSEAL=SLSAVE
        RETURN
      ENDIF
c if time = one of the ttimes, then tnsl= that tlist
      DO I=1,NSEAL
        IF(TIME.EQ.STIME(I)) THEN
          SLSAVE=SLIST(I)
          FINDSEAL=SLSAVE
          RETURN
        ENDIF
      ENDDO
c interpolate
      N=NSEAL
      JL=0
      JU=N+1
10    CONTINUE
        IF(JU-JL.GT.1) THEN
          JM=(JU+JL)/2
          IF((STIME(N).GT.STIME(1)).EQV.(TIME.GT.STIME(JM))) THEN
            JL=JM
          ELSE
            JU=JM
          ENDIF
          GOTO 10
        ENDIF
      J=JL
      SLOPE=(SLIST(J+1)-SLIST(J))/(STIME(J+1)-STIME(J))
      SL=SLIST(J)+SLOPE*(TIME-STIME(J))
      SLSAVE=SL
      FINDSEAL=SLSAVE
      END
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
      SUBROUTINE SETRIG
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      PI=4.D0*ATAN(1.D0)
      RADIUS=2.D4/PI
      RADIUS=RADIUS*0.53
      CIRCUM=2.D0*PI*RADIUS
      RKMPDEG=CIRCUM/360.D0
      RADPDEG=PI/180.D0
      DEGPRAD=180.D0/PI
      END
      SUBROUTINE POLREC(RLAT,RLONG,X,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      y= 1000*rkmpdeg*rlat
      x= 1000*rkmpdeg*cos(rlat*radpdeg)*(rlong+127.5)
      END
      SUBROUTINE RECPOL(X,Y,RLAT,RLONG)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      rlat=y*0.001/rkmpdeg
      rlong=-127.5+x*0.001/rkmpdeg/cos(rlat*radpdeg)
      END
      SUBROUTINE ZOOM(X1,X2,Y1,Y2)
      PRINT *,'DOESNT WORK'
      END
      FUNCTION IERRNM(TF)
      LOGICAL TF
      PRINT *,'DOESNT WORK'
      IERRNM=0
      END
      SUBROUTINE INUMBR(I,J)
      PRINT *,'DOESNT WORK'
      END
      SUBROUTINE FORMC(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &                 ETA, XI, W, CONST,
     &                 NUMGBC, IBFLUX, BFLUX, NZ, KZ, LM, AADOT,
     &                 ADOT,
     &                 D, B, A, CAP, KA,BOLD,ADIAG)
      IMPLICIT REAL*8(A-H,O-Z)
C FORM STIFFNESS MATRIX
      DIMENSION AADOT(NMAX),ADOT(NMAX)
      DIMENSION IBFLUX(NMAX,2)
      DIMENSION BFLUX(NMAX),KZ(NMAX)
      DIMENSION CONST(NMAX),LM(5),BOLD(NMAX),ADIAG(NMAX)
      DIMENSION KX(NMAX,4),X(NMAX),Y(NMAX),NTYPE(NMAX),D(NMAX),B(NMAX)
      REAL*8 A(NMAX,NZ),CAP(NMAX,NZ)
      INTEGER KA(NMAX,NZ+1)
      DIMENSION P(5),S(5,5),DD(5),CC(5,5)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2)
      DIMENSION PSI(4),DPSI(4,2)
      DIMENSION XY(2,4),XI(2,9),ETA(2,9),W(2,9)
C **********************************************************************
C     FORM CONDUCTIVITY MATRIX FOR COMPLETE BODY
C **********************************************************************
C
C ... ZERO OUT APPROPRIATE ARRAYS ...
c$doacross local(i,j)
      DO I=1,NUMNP
        KZ(I)=1
      ENDDO
      DO I=1,NUMNP
        KA(I,1)=I
      ENDDO
      DO I=1,NUMNP
        D(I)=0.0d0
      ENDDO
      DO I=1,NUMNP
        B(I)=0.0d0
      ENDDO
      DO J=1,NZ
        DO I=1,NUMNP
          A(I,J)=0.0d0
          CAP(I,J)=0.0d0
        ENDDO
      ENDDO
      DO J=1,NZ+1
        DO I=1,NUMNP
          KA(I,J)=0
        ENDDO
      ENDDO
C
C ... BEGIN LOOP OVER ALL THE ELEMENTS ...
      DO 100 N=1,NUMEL
        IF(NTYPE(N).EQ.1) THEN
          NODEN=4
          NINT=9
        ELSE
          NODEN=3
          NINT=4
        ENDIF
        DO I=1,4
          LM(I)=KX(N,I)
        ENDDO
C
C ..... FORM ELEMENT CONDUCTIVITY MATRIX
        DO I=1,4
          DD(I)=0.0d0
          P(I)=0.0d0
          DO J=1,4
            S(I,J)=0.0d0
            CC(I,J)=0.0d0
          ENDDO
        ENDDO
C
        I=LM(1)
        J=LM(2)
        K=LM(3)
        L=LM(4)
        XY(1,1)=X(I)
        XY(1,2)=X(J)
        XY(1,3)=X(K)
C       IF(NTYPE(N).EQ.1) XY(1,4)=X(L)
        XY(2,1)=Y(I)
        XY(2,2)=Y(J)
        XY(2,3)=Y(K)
C       IF(NTYPE(N).EQ.1) XY(2,4)=Y(L)
        IF(NODEN.EQ.4) THEN
          XY(1,4)=X(L)
          XY(2,4)=Y(L)
        ENDIF
C
C ..... FORM ELEMENT MATRIX AND VECTORS
C
C ..... BEGIN INTEGRATION POINT LOOP
        DO L=1,NINT
          CALL FESHAPE(NTYPE(N),XI(NTYPE(N),L),ETA(NTYPE(N),L),PSI,DPSI)
C
C ....... INTEGRATION POINT INTERPOLATION ...
          ADOTN=0.D0
          DO I=1,4
            ADOTN=ADOTN+ADOT(LM(I))*PSI(I)
          ENDDO
C
C ....... CALCULATE DXDS...EQUATION (5.3.6)
          DO I=1,2
            DO J=1,2
              DXDS(I,J)=0.0d0
              DO K=1,4
                DXDS(I,J)=DXDS(I,J)+DPSI(K,J)*XY(I,K)
              ENDDO
            ENDDO
          ENDDO
C
C ....... CALCULATE DSDX...EQUATION (5.2.7)
          DETJ=(DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1))
          IF(DETJ.LE.0.0) THEN
            WRITE(12,1100) DETJ,((XY(MM,NN)/1000,NN=1,4),MM=1,2)
            WRITE(*,1100) DETJ,((XY(MM,NN)/1000,NN=1,4),MM=1,2)
1100        FORMAT(' BAD JACOBIAN AT 161',G13.6,/,4G13.6,/,4G13.6)
            print *,n,lm
            STOP
          ENDIF
          DENOM=1.d0/DETJ
          DSDX(1,1)=DXDS(2,2)*DENOM
          DSDX(2,2)=DXDS(1,1)*DENOM
          DSDX(1,2)=-DXDS(1,2)*DENOM
          DSDX(2,1)=-DXDS(2,1)*DENOM
C
C ....... CALCULATE D(PSI)/DX...EQUATION (5.3.5)
          DO I=1,4
            DPSIX(I)=DPSI(I,1)*DSDX(1,1)+DPSI(I,2)*DSDX(2,1)
            DPSIY(I)=DPSI(I,1)*DSDX(1,2)+DPSI(I,2)*DSDX(2,2)
          ENDDO
C
C ....... ACCUMULATE INTEGRATION POINT VALUES OF INTEGRALS
          FAC=DETJ*W(NTYPE(N),L)
          DO I=1,4
C ......... THIS IS LUMPED CAPACITANCE MATRIX
            DD(I)=DD(I)+PSI(I)*FAC
C ........
C ........  ELEMENT AVERAGE VALUES FROM ELPROP ....
C           P(I)=P(I)+AADOT(N)*PSI(I)*FAC
C ........  INTEGRATION POINT VALUES FROM FEM INTERPOLATIOON ....
            P(I)=P(I)+ADOTN*PSI(I)*FAC
C
            DO J=1,4
              CC(I,J)=CC(I,J)+PSI(I)*PSI(J)*FAC
              IF(CONST(N).GT.1.E-30) THEN
                TERM1=CONST(N)*(DPSIX(I)*DPSIX(J)+DPSIY(I)*DPSIY(J))
                S(I,J)=S(I,J)+TERM1*FAC
              ENDIF
            ENDDO
          ENDDO
        ENDDO
c        write(73,*) 'qqq',n,real(const(n))
c        do i=1,4
c          write(73,*) (real(s(i,j)),j=1,4)
c        enddo
c        write(73,*) '------------------------------'
C
C ..... ADD ELEMENT CONDUCTIVITY TO COMPLETE CONDUCTIVITY MATRIX
C
        DO L=1,4
          I=LM(L)
C ....... THIS (D) IS LUMPED CAPACITANCE MATRIX ...
          D(I)=D(I)+DD(L)
C ....... THIS (B) IS THE LOAD VECTOR .............
          B(I)=B(I)+P(L)
C ....... THIS (A) IS THE STIFFNESS MATRIX ........
          DO M=1,4
            J=LM(M)
            IF(I.EQ.J) THEN
              A(I,1)=A(I,1)+S(L,M)
              CAP(I,1)=CAP(I,1)+CC(L,M)
              KA(I,1)=I
            ELSE
              DO K=2,KZ(I)
                IF(KA(I,K).EQ.J) THEN
                  A(I,K)=A(I,K)+S(L,M)
                  CAP(I,K)=CAP(I,K)+CC(L,M)
                  GOTO 99
                ENDIF
              ENDDO
              KZ(I)=KZ(I)+1
              A(I,KZ(I))=S(L,M)
              CAP(I,KZ(I))=CC(L,M)
              KA(I,KZ(I))=J
            ENDIF
99          CONTINUE
          ENDDO
        ENDDO
100   CONTINUE
C ... END LOOP OVER ALL THE ELEMENTS ...
C
C
C ... BOUNDARY CONDITIONS
C
C
c$doacross local(n)
      DO N=1,NUMNP
        BOLD(N)=B(N)
        ADIAG(N)=A(N,1)
      ENDDO
      IF(NUMGBC.GT.0) THEN
        DO N=1,NUMGBC
          I=IBFLUX(N,1)
          J=IBFLUX(N,2)
          TEMP=BFLUX(N)
          XL=SQRT((X(J)-X(I))**2+(Y(J)-Y(I))**2)
          TEMP=XL*TEMP*.5d0
          B(I)=B(I)+TEMP
          B(J)=B(J)+TEMP
        ENDDO
      ENDIF
c ... remove here, this writes out matrix
c      rewind 73
c      write(73,*) numnp
c      do i=1,numnp
c        write(73,*) kz(i)
c        do j=1,kz(i)
c          write(73,*) i,ka(i,j),a(i,j)
c        enddo
c      enddo
c      write(73,*) (b(i),i=1,numnp)
c ... to here ******
      END
      SUBROUTINE GAUSEID(NMAX,NZ,N,AA,KA,B,X)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION AA(NMAX,NZ),KA(NMAX,NZ+1),B(NMAX),X(NMAX)
      REAL*8 SUM
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      EPS=1D-6
      XMAX=-1d30
C ... A DIFFRENT FIRST GUESS ...
c      DO I=1,N
c        X(I)=B(I)/AA(I,1)
c      ENDDO
C ... --------------------------
      TOL=1d-8
      DO I=1,N
        JMAX=KA(I,NZ+1)
        SUM=B(I)
        DO J=2,JMAX
          SUM=SUM-AA(I,J)*X(KA(I,J))
        ENDDO
        XNEW=SUM/AA(I,1)
        XMAX=MAX(XMAX,ABS(X(I)-XNEW))
        X(I)=XNEW
      ENDDO
      ERRORG=RESID(NMAX,NZ,N,AA,KA,B,X)
c      print 24,0,errorg
      IF(ERRORG.EQ.0.) THEN
C ... RARE BUT POSSIBLE RETURN
        RETURN
      ENDIF
      IF(XMAX.EQ.0.) THEN
C ... RARE BUT POSSIBLE RETURN
        RETURN
      ENDIF
      ITMAX=50
      DO ITER=1,ITMAX
        XMAX=-1d30
        DO I=1,N
          JMAX=KA(I,NZ+1)
          SUM=B(I)
          DO J=2,JMAX
            SUM=SUM-AA(I,J)*X(KA(I,J))
          ENDDO
          XNEW=SUM/AA(I,1)
          XMAX=MAX(XMAX,ABS(X(I)-XNEW))
          X(I)=XNEW
        ENDDO
c       print *,iter,xmax
c       print 23,(x(i),i=1,n)
23    format(5g13.6)
        ERROR=RESID(NMAX,NZ,N,AA,KA,B,X)
        RATIO=ERROR/ERRORG
c        print 24,iter,error,ratio,xmax
24      format(i5,3g13.6)
        IF(ABS(RATIO).LT.EPS) THEN
c          PRINT 25,'A:CONVERGED',iter,error,ratio,xmax
25        format(a,i5,3g13.6)
          RETURN
        ENDIF
        IF(XMAX.LT.TOL) THEN
c          PRINT 25,'B:CONVERGED',iter,error,ratio,xmax
          RETURN
        ENDIF
      ENDDO
      PRINT *,'DIDNOT CONVERGE IN ',itmax,''
c      pause
      END
C=======================================
      FUNCTION RESID(NMAX,NZ,N,AA,KA,B,X)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION AA(NMAX,NZ),KA(NMAX,NZ+1),B(NMAX),X(NMAX)
      REAL*8 SUMSQ,SUM
      SUMSQ=0.D0
      DO I=1,N
        JMAX=KA(I,NZ+1)
        SUM=0.d0
        DO J=1,JMAX
          SUM=SUM+AA(I,J)*X(KA(I,J))
        ENDDO
        SUM=SUM-B(I)
        SUMSQ=SUMSQ+SUM**2
      ENDDO
      RESID=SUMSQ
      END
C=======================================
      SUBROUTINE CONJUG(NMAX,NZ,NUM,ERROR,AA,KA,B,X)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      include "parameter.h"
      PARAMETER(LMAX=MAXNUM)
      DIMENSION AA(NMAX,NZ),KA(NMAX,NZ+1),B(NMAX),X(NMAX)
      DIMENSION R(LMAX),V(LMAX),Z(LMAX)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
100   FORMAT(1X,A,(1X,T10,5G13.6))
      IF(NUM.GT.LMAX) THEN
        WRITE(*,*) 'INCREASE LMAX=',LMAX,' TO',NUM
        STOP
      ENDIF
      M=MAX(1000,1*NUM)
      EPS=0D0
      ITMIN=10
      DELTA=EPS
C ...
c ... output matrix, etc
c      write(93,*) num,nz
c      do i=1,num
c        write(93,*) (aa(i,j),j=1,nz)
c      enddo
c      do i=1,num
c        write(93,*) (ka(i,j),j=1,nz+1)
c      enddo
c      write(93,*) (b(i),i=1,num)
c      write(93,*) (x(i),i=1,num)
c      pause
c ... 
C ... FORM A*X IN R
C     CALL MULTB(MMAX,NWIDE,NUM,A,X,R)
c$doacross local(i,jmax,sum)
      DO I=1,NUM
        JMAX=KA(I,NZ+1)
        SUM=0.D0
        DO J=1,JMAX
          SUM=SUM+AA(I,J)*X(KA(I,J))
        ENDDO
        R(I)=SUM
      ENDDO
C ... FORM R=B-A*X, ALSO V=R, ALSO VV=SQRT(<V,V>), ALSO C=<R,R>
      VV=0.D0
      C=0.D0
      DO I=1,NUM
        R(I)=B(I)-R(I)
        V(I)=R(I)
        VV=VV+V(I)*V(I)
        C=C+R(I)*R(I)
      ENDDO
      VV=SQRT(VV)
c      PRINT 100,' R ',(R(N),N=1,NUM)
c      PRINT 100,' VBASE ',VV
c      PRINT 100,' DBASE ',C
      DBASE=C
      VBASE=VV
C ... BEGIN ITERATION LOOP ...
      DO K=1,M
c        PRINT 100,' V ',(V(N),N=1,NUM)
c        PRINT 100,' VV ',VV
C ..... IF SQRT(<V,V>) < DELTA EXIT
        IF(VV.LE.DELTA .and. k.gt.ITMIN) THEN
          IF(IOTOGG) THEN
            write(list(ipage+1),*) 
     &            '1A:SOLUTION, ITERATION=',K,VV
            ipage=ipage+1
          ENDIF
          RETURN
        ENDIF
        IF(VV/VBASE.LE.ERROR .and. k.gt.ITMIN) THEN
          IF(IOTOGG) THEN
            write(list(ipage+1),*)
     &        '1B:SOLUTION, ITERATION=',K,VV/VBASE
            ipage=ipage+1
          ENDIF
          RETURN
        ENDIF
C ..... FORM Z=A*V
C        CALL MULTB(MMAX,NWIDE,NUM,A,V,Z)
c$doacross local(i,jmax,sum)
        DO I=1,NUM
          JMAX=KA(I,NZ+1)
          SUM=0.D0
          DO J=1,JMAX
            SUM=SUM+AA(I,J)*V(KA(I,J))
          ENDDO
          Z(I)=SUM
        ENDDO
C        PRINT 100,' Z ',(Z(N),N=1,NUM)
C ..... FORM <V,Z> IN T
        T=0.0d0
        DO I=1,NUM
          T=T+V(I)*Z(I)
        ENDDO
C        PRINT 100,' <V,Z> ',T
C ..... PREVENT DIVIDE BY ZERO
        IF(T.EQ.0.0) THEN
          WRITE(*,*) 'DIVIDE BY ZERO IN CONJUG, ITERATION',K
          STOP
        ENDIF
C ..... FORM T=C/<V,Z>
        T=C/T
C        PRINT 100,' C/<V,Z> ',T
C ..... UNDATE X=X+T*V AND R=R-T*Z, FORM D=<R,R>
        D=0.D0
        DO I=1,NUM
          X(I)=X(I)+T*V(I)        
          R(I)=R(I)-T*Z(I) 
          D=D+R(I)*R(I)
        ENDDO
c        PRINT 100,' X ',(X(N),N=1,NUM)
c        PRINT 100,' R ',(R(N),N=1,NUM)
c        PRINT 100,' D ',REAL(K),D
C ..... IF D<EPS EXIT
        IF(D.LE.EPS .and. k.gt.ITMIN) THEN
          IF(IOTOGG) THEN
            write(list(ipage+1),1700) 
     &          '2A:SOLUTION, ITERATION=',K,REAL(D)
            ipage=ipage+1
          ENDIF
          RETURN
        ENDIF
        IF(D/DBASE.LE.ERROR .and. k.gt.ITMIN) THEN
          IF(IOTOGG) THEN
            write(list(ipage+1),1700)
     &            '2B:SOLUTION, ITERATION=',K,
     &             D,D/DBASE
            ipage=ipage+1
1700      format(1x,a,i5,1p2g13.6)
          ENDIF
          RETURN
        ENDIF
C ..... FORM V=R+(D/C)*V, ALSO VV=SQRT(<V,V>)
        CON=D/C
        VV=0.D0
        DO I=1,NUM
          V(I)=R(I)+CON*V(I)
          VV=VV+V(I)*V(I)
        ENDDO
        VV=SQRT(VV)
c        PRINT 100,' D,VV ',REAL(K),D,D/DBASE,VV,VV/VBASE
C        PRINT 100,' V ',(V(N),N=1,NUM)
C ..... REPLACE C WITH D
        C=D
C        PRINT 100,' K,C ',REAL(K),C
C        PAUSE
      ENDDO
      PRINT *,'DIDNOT CONVERGE IN ',M,''
      END
C=======================================
      SUBROUTINE MULTB(MMAX,NWIDE,NUM,A,X,R)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(MMAX,NWIDE),X(MMAX),R(MMAX)
      IB=NWIDE/2
      IB1=IB+1
      DO I=1,NUM
        C=0.D0
        DO K=1,NWIDE
          J=K+I-IB1
          IF(J.GE.1 .AND. J.LE.NUM) THEN
            C=C+A(I,K)*X(J)
          ENDIF
        ENDDO
        R(I)=C
      ENDDO
      END
C=======================================
      FUNCTION When ()
      REAL*8 When
      REAL tdum(2)
*
      When = dble(etime(tdum))
      END
C=======================================
      SUBROUTINE DIAGD(NMAX,NZ,NUM,AA,KA,LDIAG,DMAX,SMAX)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION AA(NMAX,NZ),KA(NMAX,NZ+1)
      LOGICAL LDIAG
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      LDIAG=.TRUE.
      DMAX=0.D0
      DO I=1,NUM
        JMAX=KA(I,NZ+1)
        DIAG=ABS(AA(I,1))
        SUM=0.0d0
        DO J=2,JMAX
          SUM=SUM+ABS(AA(I,J))
c          SUM=SUM+AA(I,J)
        ENDDO
        SUM=ABS(SUM)
        DIFF=ABS(DIAG-SUM)
c        IF(DIFF.GT.1E-6 .AND. DIAG.LT.SUM) THEN
        IF(DIAG.LE.SUM) THEN
          LDIAG=.FALSE.
          DMAX=DIAG
          SMAX=SUM
          IMAX=I
c          IF(IOTOGG) THEN
c            write(list(ipage+1),*) 
c     &          'NOT DIAGONALLY DOMINANT',REAL(DMAX),REAL(SMAX)
c            write(list(ipage+2),*)
c     &          '                       ',IMAX,REAL(DMAX-SMAX)
c            ipage=ipage+2
c          endif
c           print *, 
c     &          'NOT DIAGONALLY DOMINANT',REAL(DMAX),REAL(SMAX)
          RETURN
        ENDIF
        IF(DIAG.GT.DMAX) THEN
          DMAX=DIAG
          SMAX=SUM
          IMAX=I
        ENDIF
      ENDDO
c      IF(LDIAG) THEN
c        IF(IOTOGG) THEN
c          write(list(ipage+1),*)
c     &         'OK, DIAGONALLY DOMINANT',REAL(DMAX),REAL(SMAX)
c          ipage=ipage+1
c         PRINT *,'                       ',IMAX,REAL(DMAX-SMAX)
c        ENDIF
c          print *, 
c     &         'OK, DIAGONALLY DOMINANT',REAL(DMAX),REAL(SMAX)
c      ELSE
c        IF(IOTOGG) THEN
c          write(list(ipage+1),*)
c     &         'NOT DIAGONALLY DOMINANT',DMAX,SMAX
c          write(list(ipage+2),*)
c     &         '                       ',IMAX,DMAX-SMAX
c          ipage=ipage+2
c        ENDIF
c          print *, 
c     &         'NOT DIAGONALLY DOMINANT',DMAX,SMAX
c      ENDIF
      END
c=======================================
      SUBROUTINE ASYMSL(NMAX,NZ,NEQ,AA,KA,B,X,KKK)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION AA(NMAX,NZ),KA(NMAX,NZ+1),B(NMAX),X(NMAX)
      include "parameter.h"
      PARAMETER(MMAX=MAXNUM,MCOL=171)
      DIMENSION A(MCOL,MMAX),Q(MMAX)
      SAVE A,Q
      MBAND=(MCOL-1)/2
      NBAND=MBAND+1
      NCOLS=2*NBAND-1
      IF(KKK.EQ.1) THEN
C ..... COPY THE MATRIX TO BANDED FORM...
c        PRINT *,'COPY THE MATRIX TO BANDED FORM...'
        DO I=1,NEQ
          DO J=1,MCOL
            A(J,I)=0.D0
          ENDDO
          DO J=1,KA(I,NZ+1)
            JJ=KA(I,J)-I+1+MBAND
            IF(JJ.LE.0) PRINT *,JJ,' PROBLEMS...'
            IF(JJ.GT.MCOL) PRINT *,JJ,' PROBLEMS...'
            A(JJ,I)=AA(I,J)
          ENDDO
          Q(I)=B(I)
        ENDDO
c        PRINT *,'DONE COPYING...'
        KMIN=NBAND+1
        DO N=1,NEQ
          IF(A(NBAND,N).EQ.0.0) THEN
            WRITE(*,1001) N,A(NBAND,N)
            STOP
          ENDIF
          IF(A(NBAND,N).NE.1.0) THEN
            C=1.d0/A(NBAND,N)
            DO K=KMIN,NCOLS
              IF(A(K,N).NE.0.0) A(K,N)=C*A(K,N)
            ENDDO
          ENDIF
          DO L=2,NBAND
            JJ=NBAND-L+1
            I=N+L-1
            IF(I.LE.NEQ .AND. A(JJ,I).NE.0.0) THEN
              KI=NBAND+2-L
              KF=NCOLS+1-L
              J=NBAND
              DO K=KI,KF
                J=J+1
                A(K,I)=A(K,I)-A(JJ,I)*A(J,N)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
c        PRINT *,'DONE ELIMINATING...'
        RETURN
      ELSEIF(KKK.EQ.2) THEN
        DO N=1,NEQ
          IF(A(NBAND,N).EQ.0.0) THEN
            WRITE(*,1002) N,A(NBAND,N)
            STOP
          ENDIF
          IF(A(NBAND,N).NE.1.0) Q(N)=Q(N)/A(NBAND,N)
          DO L=2,NBAND
            JJ=NBAND-L+1
            I=N+L-1
            IF(I.LE.NEQ .AND. A(JJ,I).NE.0.0) Q(I)=Q(I)-A(JJ,I)*Q(N)
          ENDDO
        ENDDO
C
C BACK SUBSTITUTION
C
        LL=NBAND+1
        DO M=1,NEQ
          N=NEQ+1-M
          DO L=LL,NCOLS
            IF(A(L,N).NE.0.0) THEN
              K=N+L-NBAND
              Q(N)=Q(N)-A(L,N)*Q(K)
            ENDIF
          ENDDO
          X(N)=Q(N)
        ENDDO
c        PRINT *,'DONE BACK SUBSTITUTION...'
        RETURN
      ENDIF
1001  FORMAT(6X,'SET OF EQUATIONS MAY BE SINGULAR',
     &     /,5X,'DIAGONAL TERM OF EQUATION',I5,' IS EQUAL TO' ,E15.8)
1002  FORMAT(6X,'SET OF EQUATIONS ARE SINGULAR',
     &     /,5X,'DIAGONAL TERM OF EQUATION',I5,'  IS EQUAL TO',E15.8)
      END
C
C
      SUBROUTINE MODIFY(NMAX,IFIT,NUMNP,KODE,HTICE,DEPB,PSURF,
     &                  ADOT,FRACT,FLOWA,SLDGB,BDROCK,AFUDGE,IDT,
     &                  NUMEL,KX)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IFIT(7),KODE(NMAX),HTICE(NMAX),PSURF(NMAX)
      DIMENSION ADOT(NMAX),FRACT(NMAX),FLOWA(NMAX),SLDGB(NMAX)
      DIMENSION BDROCK(NMAX),DEPB(NMAX),AFUDGE(NMAX),IDT(NMAX)
      DIMENSION KX(NMAX,4)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      FLOWMIN=1d-6
c      FLOWMIN=2.29d0
      SLDGBMIN=1d-9
      RSURF=-500.d0
      RTHICK=1000.d0
      RFACT=1.d0
      RFACTE=1.d0
      AFMAX=10.d0
      AFMIN=0.01d0
      IF(IOTOGG) then
        IF(IFIT(1).EQ.1) then
          WRITE(list(ipage+1),*) 'MODIFYING ADOT:NODE'
          ipage=ipage+1
        ENDIF
        IF(IFIT(2).EQ.1) then
          WRITE(list(ipage+1),*) 'MODIFYING FRACT:NODE'
          ipage=ipage+1
        endif
        IF(IFIT(3).EQ.1) then
          WRITE(list(ipage+1),*) 'MODIFYING FLOW:NODE'
          ipage=ipage+1
        endif
        IF(IFIT(4).EQ.1) then
          WRITE(list(ipage+1),*) 'MODIFYING SLDGB:NODE'
          ipage=ipage+1
        endif
        IF(IFIT(5).EQ.1) then
          WRITE(list(ipage+1),*) 'MODIFYING AFUDGE:NODE'
          ipage=ipage+1
        endif
        IF(IFIT(6).EQ.1) then
          WRITE(list(ipage+1),*) 'MODIFYING AFUDGE/SLDGB:NODE'
          ipage=ipage+1
        endif
        IF(IFIT(1).EQ.-1) then
          WRITE(list(ipage+1),*) 'MODIFYING ADOT:ELEMENT'
          ipage=ipage+1
        endif
        IF(IFIT(2).EQ.-1) then
          WRITE(list(ipage+1),*) 'MODIFYING FRACT:ELEMENT'
          ipage=ipage+1
        endif
        IF(IFIT(3).EQ.-1) then
          WRITE(list(ipage+1),*) 'MODIFYING FLOW:ELEMENT'
          ipage=ipage+1
        endif
        IF(IFIT(4).EQ.-1) then
          WRITE(list(ipage+1),*) 'MODIFYING SLDGB:ELEMENT'
          ipage=ipage+1
        endif
        IF(IFIT(5).EQ.-1) then
          WRITE(list(ipage+1),*) 'MODIFYING AFUDGE:ELEMENT'
          ipage=ipage+1
        endif
        IF(IFIT(6).EQ.-1) then
          WRITE(list(ipage+1),*) 'MODIFYING AFUDGE/SLDGB:ELEMENT'
          ipage=ipage+1
        endif
      endif
      IF(IFIT(1).EQ.1 .OR. 
     &   IFIT(2).EQ.1 .OR. 
     &   IFIT(3).EQ.1 .OR.
     &   IFIT(4).EQ.1 .OR.
     &   IFIT(5).EQ.1 .OR.
     &   IFIT(6).EQ.1) THEN
        DO J=1,NUMNP
          THICKC=HTICE(J)-DEPB(J)
          THICKA=PSURF(J)-BDROCK(J)
          IF(KODE(J).EQ.1) THEN
            RATIOD=0.d0
            RATIOF=1.d0
          ELSE
            RATIOD=1.d0-(HTICE(J)-RSURF)/(PSURF(J)-RSURF)
            RATIOD=RATIOD*RFACT
            RATIOF=RFACT*(PSURF(J)-RSURF)/(HTICE(J)-RSURF)
          ENDIF
          IF(IFIT(1).EQ.1) THEN
            IF(IDT(J).EQ.0) THEN
              ADOT(J)=ADOT(J)+RATIOD
            ENDIF
          ENDIF
          IF(IFIT(2).EQ.1) THEN
c            FRACT(J)=FRACT(J)-FRACT(J)*RATIOD
            FRACT(J)=FRACT(J)-0.1d0*RATIOD
            IF(FRACT(J).GT.1.) FRACT(J)=1.d0
            IF(FRACT(J).LT.0.) FRACT(J)=0.d0
          ENDIF
          IF(IFIT(3).EQ.1) THEN
            FLOWA(J)=FLOWA(J)*RATIOF
            IF(FLOWA(J).GT.10.) FLOWA(J)=10.d0
            IF(FLOWA(J).LT.FLOWMIN) FLOWA(J)=FLOWMIN
          ENDIF
          IF(IFIT(4).EQ.1) THEN
            SLDGB(J)=SLDGB(J)*RATIOF
            IF(SLDGB(J).GT..20) SLDGB(J)=.20d0
            IF(SLDGB(J).LT.SLDGBMIN) SLDGB(J)=SLDGBMIN
          ENDIF
          IF(IFIT(5).EQ.1) THEN
c           IF(KODE(J).EQ.1) THEN
c             RATIOD=0.
c           ELSE
c             RATIOD=1.-(THICKC-RTHICK)/(THICKA-RTHICK)
c             RATIOD=RATIOD*RFACT
c           ENDIF
C            BDROCK(J)=BDROCK(J)+RATIOD*100.
C            IF(BDROCK(J).LT.-9999.) BDROCK(J)=-9999.
C            IF(BDROCK(J).GT.PSURF(J)) BDROCK(J)=PSURF(J)
            AFUDGE(J)=AFUDGE(J)*RATIOF
            IF(AFUDGE(J).GT.AFMAX) AFUDGE(J)=AFMAX
            IF(AFUDGE(J).LT.AFMIN) AFUDGE(J)=AFMIN
C      PRINT 23,J,HTICE(J),PSURF(J),RATIOD*100.,BDROCK(J)
23    FORMAT(I5,4G13.6)
          ENDIF
          IF(IFIT(6).EQ.1) THEN
            IF(FRACT(J).LT.0.5) THEN
              AFUDGE(J)=AFUDGE(J)*RATIOF
              IF(AFUDGE(J).GT.AFMAX) AFUDGE(J)=AFMAX
              IF(AFUDGE(J).LT.AFMIN) AFUDGE(J)=AFMIN
            ELSE
              SLDGB(J)=SLDGB(J)*RATIOF
              IF(SLDGB(J).GT..20) SLDGB(J)=.20d0
              IF(SLDGB(J).LT.SLDGBMIN) SLDGB(J)=SLDGBMIN
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      IF(IFIT(1).EQ.-1 .OR. 
     &   IFIT(2).EQ.-1 .OR. 
     &   IFIT(3).EQ.-1 .OR.
     &   IFIT(4).EQ.-1 .OR.
     &   IFIT(5).EQ.-1 .OR.
     &   IFIT(6).EQ.-1) THEN
        DO J=1,NUMEL
          HT=0.d0
          DB=0.d0
          PS=0.d0
          BD=0.d0
          DO K=1,4
            HT=HT+HTICE(KX(J,K))
            DB=DB+DEPB(KX(J,K))
            PS=PS+PSURF(KX(J,K))
            BD=BD+BDROCK(KX(J,K))
          ENDDO
          HT=0.25d0*HT
          DB=0.25d0*DB
          PS=0.25d0*PS
          BD=0.25d0*BD
          THICKC=HT-DB
          THICKA=PS-BD
          RATIOD=1.d0-(HT-RSURF)/(PS-RSURF)
          RATIOD=RATIOD*RFACTE
          RATIOF=RFACTE*(PS-RSURF)/(HT-RSURF)
C         PRINT *,REAL(HT),REAL(PS),REAL(RATIOD),REAL(RATIOF)
          IF(IFIT(1).EQ.-1) THEN
            DO K=1,4
              IF(IDT(KX(J,K)).EQ.0) THEN
                ADOT(KX(J,K))=ADOT(KX(J,K))+RATIOD
              ENDIF
            ENDDO
          ENDIF
          IF(IFIT(2).EQ.-1) THEN
            DO K=1,4
c             FRACT(KX(J,K))=FRACT(KX(J,K))-FRACT(KX(J,K))*RATIOD
              FRACT(KX(J,K))=FRACT(KX(J,K))-0.1d0*RATIOD
              IF(FRACT(KX(J,K)).GT.1.) FRACT(KX(J,K))=1.d0
              IF(FRACT(KX(J,K)).LT.0.) FRACT(KX(J,K))=0.d0
            ENDDO
          ENDIF
          IF(IFIT(3).EQ.-1) THEN
            DO K=1,4
              FLOWA(KX(J,K))=
     &              0.5d0*(FLOWA(KX(J,K))+FLOWA(KX(J,K))*RATIOF)
              IF(FLOWA(KX(J,K)).GT.10.) FLOWA(KX(J,K))=10.d0
              IF(FLOWA(KX(J,K)).LT.FLOWMIN) FLOWA(KX(J,K))=FLOWMIN
            ENDDO
          ENDIF
          IF(IFIT(4).EQ.-1) THEN
            DO K=1,4
              SLDGB(KX(J,K))=
     &              0.5d0*(SLDGB(KX(J,K))+SLDGB(KX(J,K))*RATIOF)
              IF(SLDGB(KX(J,K)).GT..20) SLDGB(KX(J,K))=.20d0
              IF(SLDGB(KX(J,K)).LT.SLDGBMIN) SLDGB(KX(J,K))=SLDGBMIN
            ENDDO
          ENDIF
          IF(IFIT(5).EQ.-1) THEN
            DO K=1,4
             AFUDGE(KX(J,K))=
     &              0.5d0*(AFUDGE(KX(J,K))+AFUDGE(KX(J,K))*RATIOF)
              IF(AFUDGE(KX(J,K)).GT.AFMAX) AFUDGE(KX(J,K))=AFMAX
              IF(AFUDGE(KX(J,K)).LT.AFMIN) AFUDGE(KX(J,K))=AFMIN
            ENDDO
          ENDIF
          IF(IFIT(6).EQ.-1) THEN
            DO K=1,4
              IF(FRACT(KX(J,K)).LT.0.5) THEN
                AFUDGE(KX(J,K))=
     &              0.5d0*(AFUDGE(KX(J,K))+AFUDGE(KX(J,K))*RATIOF)
                IF(AFUDGE(KX(J,K)).GT.AFMAX) AFUDGE(KX(J,K))=AFMAX
                IF(AFUDGE(KX(J,K)).LT.AFMIN) AFUDGE(KX(J,K))=AFMIN
              ELSE
                SLDGB(KX(J,K))=
     &              0.5d0*(SLDGB(KX(J,K))+SLDGB(KX(J,K))*RATIOF)
                IF(SLDGB(KX(J,K)).GT..20) SLDGB(KX(J,K))=.20d0
                IF(SLDGB(KX(J,K)).LT.SLDGBMIN) SLDGB(KX(J,K))=SLDGBMIN
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      IF(IFIT(1).EQ.2 .OR. 
     &   IFIT(2).EQ.2 .OR. 
     &   IFIT(3).EQ.2 .OR.
     &   IFIT(4).EQ.2 .OR.
     &   IFIT(5).EQ.2) THEN
        THICKDIFF=0.0
        THICKAVG=0.0
        DO J=1,NUMNP
          THICKC=HTICE(J)-DEPB(J)
          THICKA=PSURF(J)-BDROCK(J)
          THICKAVG=THICKAVG+THICKA+THICKC
          THICKDIFF=THICKDIFF+(THICKC-THICKA)
        ENDDO
        THICKAVG=0.5d0*THICKAVG/NUMNP
        THICKDIFF=THICKDIFF/NUMNP
        ratio=1.d0-1*THICKDIFF/THICKAVG
        IF(IFIT(5).EQ.2) then
          IF(IOTOGG) then
            WRITE(list(ipage+1),*) 'MODIFYING AFUDGE:GLOBAL',
     &                             real(thickdiff),real(ratio)
            ipage=ipage+1
          ENDIF
          avg1=0
          avg2=0
          do i=1,numnp
            avg1=avg1+afudge(i)
            afudge(i)=afudge(i)*ratio
            IF(AFUDGE(i).GT.AFMAX) AFUDGE(i)=AFMAX
            IF(AFUDGE(i).LT.AFMIN) AFUDGE(i)=AFMIN
            avg2=avg2+afudge(i)
          enddo
          avg1=avg1/numnp
          avg2=avg2/numnp
          IF(IOTOGG) then
            WRITE(list(ipage+1),*) '                       ',
     &                             real(avg1),real(avg2)
            ipage=ipage+1
          ENDIF
        ENDIF
        IF(IFIT(4).EQ.2) then
          IF(IOTOGG) then
            WRITE(list(ipage+1),*) 'MODIFYING SLDGB:GLOBAL',
     &                             real(thickdiff),real(ratio)
            ipage=ipage+1
          ENDIF
          avg1=0
          avg2=0
          do i=1,numnp
            avg1=avg1+SLDGB(i)
            SLDGB(i)=SLDGB(i)*ratio
            IF(SLDGB(J).GT..20) SLDGB(J)=.20d0
            IF(SLDGB(J).LT.SLDGBMIN) SLDGB(J)=SLDGBMIN
            avg2=avg2+SLDGB(i)
          enddo
          avg1=avg1/numnp
          avg2=avg2/numnp
          IF(IOTOGG) then
            WRITE(list(ipage+1),*) '                       ',
     &                             real(avg1),real(avg2)
            ipage=ipage+1
          ENDIF
        ENDIF
      ENDIF
      END

      SUBROUTINE GNUMBR(X,IDIGIT,IWIDE)
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
      SUBROUTINE SLOPEF(MXX,X,Y,KX,NUMEL,HTICE,SLOPE,WINDIR)
      IMPLICIT REAL*8(A-H,O-Z)
C ... CALCULATES LINEARIZATION CONSTANT FROM CURRENT SOLUTION
      DIMENSION SLOPE(4,MXX),HTICE(MXX)
      DIMENSION KX(MXX,4),X(MXX),Y(MXX)
      DIMENSION LM(5)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2)
      DIMENSION PSI(4),DPSI(4,2)
      DIMENSION XY(2,4),WINDIR(2)
      NODEN=4
      CENTX=0.0d0
      CENTY=0.0d0
      DO J = 1,NUMEL
        SUMX=0.d0
        SUMY=0.d0
        DO I = 1,NODEN
          LM(I) = KX(J,I)
        ENDDO
        I=LM(1)
        JJ=LM(2)
        K=LM(3)
        L=LM(4)
        XY(1,1)=X(I)
        XY(1,2)=X(JJ)
        XY(1,3)=X(K)
        XY(1,4)=X(L)
        XY(2,1)=Y(I)
        XY(2,2)=Y(JJ)
        XY(2,3)=Y(K)
        XY(2,4)=Y(L)
        CALL FESHAPE(1,CENTX,CENTY,PSI,DPSI)
C
C CALCULATE DXDS...EQUATION (5.3.6)
C
        DO I=1,2
          DO L=1,2
            DXDS(I,L)=0.0d0
            DO K=1,NODEN
              DXDS(I,L)=DXDS(I,L)+DPSI(K,L)*XY(I,K)
            ENDDO
          ENDDO
        ENDDO
C
C CALCULATE DSDX...EQUATION (5.2.7)
C
        DETJ=(DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1))
        IF (DETJ.LE.0.0) THEN
          WRITE(12,5544) J,DETJ
          WRITE(*,5544) J,DETJ
          WRITE(12,5545) (JJ,XY(1,JJ),XY(2,JJ),JJ=1,4)
          WRITE(*,5545) (JJ,XY(1,JJ),XY(2,JJ),JJ=1,4)
5545      FORMAT(1X,I5,1X,1PE10.3,E10.3)
          STOP
5544      FORMAT(' BAD JACOBIAN',I5,1PE10.3,/,1X,8E10.3)
        ENDIF
        DSDX(1,1)=DXDS(2,2)/DETJ
        DSDX(2,2)=DXDS(1,1)/DETJ
        DSDX(1,2)=-DXDS(1,2)/DETJ
        DSDX(2,1)=-DXDS(2,1)/DETJ
C
C CALCULATE D(PSI)/DX...EQUATION (5.3.5)
C
        DO I=1,NODEN
          DPSIX(I)=DPSI(I,1)*DSDX(1,1)+DPSI(I,2)*DSDX(2,1)
          DPSIY(I)=DPSI(I,1)*DSDX(1,2)+DPSI(I,2)*DSDX(2,2)
        ENDDO
        DO I = 1,NODEN
          SUMX = SUMX + HTICE(LM(I))*DPSIX(I)
          SUMY = SUMY + HTICE(LM(I))*DPSIY(I)
        ENDDO
C
        DELH = SUMX**2 + SUMY**2
        DELH = SQRT(DELH)
        SLOPE(1,J)=DELH
        SLOPE(2,J)=SUMX
        SLOPE(3,J)=SUMY
        SLOPE(4,J)=SUMX*SIN(WINDIR(1))+SUMY*COS(WINDIR(1))
      ENDDO
      RETURN
      END
C----------------------------------------------
      SUBROUTINE TEMPER(MXX, NUMNP, NUMEL, NTYPE, KX, KODE, ADOT, FRACT,
     &            DEPB, HTICE, PSURF, FLOWA, SLDGB, TEMP, ITYPE, TBED,
     &            X, Y, AMASS, SLOPN, TIME, DT, IPLOT, TBAVG,
     &            AFUDGE,GRIDAREA,IMELT,BMELT,WTHICK,VOL,AREA,
     &            GEOFLUX,VEL,IDT,CONTRIB)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM,MMAX=40)
      PARAMETER(NSMAX=MAXTIME)
C USE SURFACE TEMPERATURE TO CALCULATE BED TEMPERATURE, AND HENCE MODIFY
C SLIDING DISTRIBUTION
C
      EXTERNAL ACCUM,PMP
      DIMENSION NTYPE(NMAX),ADOT(NMAX),FRACT(NMAX),X(NMAX),Y(NMAX),
     &          DEPB(NMAX),HTICE(NMAX),FLOWA(NMAX),TBED(NMAX),
     &          SLDGB(NMAX),TEMP(NMAX),LM(4),KX(NMAX,4),AMASS(11),
     &          SLOPN(4,NMAX),KODE(NMAX),ITYPE(NMAX),AFUDGE(NMAX),
     &          ZZZ(MMAX),BMELT(NMAX),WTHICK(NMAX),VEL(NMAX,3),
     &          GEOFLUX(NMAX),BMELTOLD(NMAX),IDT(NMAX),PSURF(NMAX),
     &          CONTRIB(NMAX)
      COMMON /LINE/ NP,NLINE(1000)
      COMMON /TEMPERATURE/ TEMPA(MMAX,NMAX)
      COMMON /TCONST/ TSORIG(NMAX),TFLAG
      LOGICAL TFLAG
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      LOGICAL ISICE,INNLINE
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      COMMON /HEAT/ HTOGG
      LOGICAL HTOGG,HSAVE
      logical file41,file43
      SAVE ISTART,JSTART,BMELTOLD,KSTART
      DATA ISTART /0/, JSTART /0/, KSTART /0/
      DATA CONV /3.D0/
C THRESHOLD VALUE OF WATER THICKNESS (10 MM ) TO PRODUCE SLIDING
      STHRESH=10.D0*1D-3
C THRESHOLD VALUE OF WATER THICKNESS (1 MM ) TO PRODUCE SLIDING
      STHRESH=1.D0*1D-3
C THRESHOLD VALUE OF WATER THICKNESS (0.1 MM ) TO PRODUCE FREEZING
      STHRESHB=0.1D0*1D-3
      IF(MXX.NE.NMAX) THEN
        PRINT *,'CHANGE NMAX IN TEMPER',NMAX,MXX
        STOP
      ENDIF
      !PRINT *,' IN TEMPER ',TIME,DT;pause
C------------------------------------------
      CALL HARD(A0,B0,A1,B1)
      inquire(file='fort.41',exist=file41)
      IF(JSTART.EQ.0. and. file41) THEN
        PRINT *,'READING WATER THICKNESS FILE'
        READ(41,*,END=11) NNN
        IF(NNN.NE.NUMNP) THEN
          PRINT *,'PROBLEMS:'
          PRINT *,'INCOMPATIBLE WITH CURRENT NUMNP=',NUMNP,NNN
          GOTO 11
        ENDIF
        DO I=1,NUMNP
          READ(41,*,END=11) II,WTHICK(I),BMELT(I)
          BMELTOLD(I)=BMELT(I)
          IF(I.NE.II) THEN
            PRINT *,'PROBLEMS:'
            PRINT *,'INCOMPATIBLE WITH CURRENT NMAX=',NMAX
            GOTO 11
          ENDIF
        ENDDO
        JSTART=1
      ENDIF
11    CONTINUE
      IF(JSTART.EQ.0) THEN
        PRINT *,' NONE FOUND, SET TO ZERO ... '
        DO I=1,NUMNP
            WTHICK(I)=0.D0
            BMELT(I)=-1D-5
            BMELTOLD(I)=-1D-5
        ENDDO
        JSTART=1
      ENDIF
      IF(ISTART.EQ.0) THEN
        PRINT *,'READING TEMPERATURE FILE'
        PRINT *,'0:BINARY, 1:ASCII'
        READ(*,*) IREAD
        WRITE(99,*) IREAD
        inquire(file='fort.43',exist=file43)
        if(file43) then
          IF(IREAD.EQ.1) THEN
            DO I=1,NUMNP
              DO J=1,MMAX
                READ(43,*,END=12) II,JJ,TEMPA(J,I)
                IF(I.NE.II .OR. J.NE.JJ) THEN
                  PRINT *,'PROBLEMS:'
                  PRINT *,'INCOMPATIBLE WITH CURRENT MMAX=',MMAX
                  GOTO 12
                ENDIF
              ENDDO
            ENDDO
            ISTART=1
          ELSE
            READ(43) NUMNP1,NPT1
            IF(NUMNP1.NE.NUMNP .OR. NPT1.NE.MMAX) THEN
              PRINT *,'PROBLEMS:'
              PRINT *,'INCOMPATIBLE WITH CURRENT MMAX=',MMAX
              GOTO 12
            ENDIF
            READ(43) ((TEMPA(J,I),J=1,MMAX),I=1,NUMNP)
            ISTART=1
          ENDIF
        endif
      ENDIF
12    CONTINUE
      IF(ISTART.EQ.0) THEN
        PRINT *,'NONE FOUND, STEADY STATE CALCULATION'
        DO I=1,NUMNP
          DO J=1,MMAX
            TEMPA(J,I)=TEMP(I)
          ENDDO
        ENDDO
        DTSAVE=DT
        DT=0.D0
      ENDIF        
C--------------------------------------------          
      IC1=0
      IC2=0
      IFIXED=0
      IMELTING=0
      IMELTED=0
      IFREEZING=0
      IFROZEN=0
      INOICE=0
      IWITHICE=0
      TMIN=1D30
      TMAX=-TMIN
      ISTEP=NUMNP/5
      IF(IPLOT.EQ.7) THEN
        IPLSTRT=0
        IPLT=1
      ELSE
        IPLSTRT=0
        IPLT=0
      ENDIF
      TBAVG=0.D0
      NTAVG=0
      IF(KSTART.EQ.0) THEN
        HSAVE=HTOGG
        HTOGG=.FALSE.
      ENDIF
      DO I=1,NUMNP
C ...... IS I IN NLINE ??
         IF(IPLOT.EQ.7) THEN
           INNLINE=.FALSE.
           DO NN=1,NP
             IF(I.EQ.NLINE(NN)) THEN
               INNLINE=.TRUE.
               GOTO 500
             ENDIF
           ENDDO
500        CONTINUE
           INNLINE=.FALSE.
           IF(INNLINE) THEN
             IPLT=1
           ELSE
             IPLT=0
           ENDIF
         ENDIF
C        IFLAG=0
C        CALL COLUMNZ(MMAX,HTICE(I),DEPB(I),ZZZ)
C        IF(ITYPE(I).NE.0 .AND. KODE(I).EQ.0) THEN
        IF(ITYPE(I).NE.0) THEN
C          IF(MOD(I,ISTEP).EQ.0) PRINT *,I
C ******* REMOVE ACCUM, REPLACE WITH AFUNCT EXPERIMENTAL *********
C          AJUNK=ACCUM(IDT(I),X(I)*.001D0,Y(I)*.001D0,HTICE(I),
C     &                        SLOPN(1,I),
C     &                        PSURF(I),AMASS(9),TEMP(I))
          IF(IDT(I).EQ.24) THEN
            AJUNK=AFUNCT(TIME, 24, AMASS,
     &                   HTICE(I), DEPB(I), PSURF(I),
     &                   SLOPN(1,I),X(I),Y(I),TEMP(I))
          ELSEIF(IDT(I).EQ.20) THEN
            AJUNK=AFUNCT(TIME, 20, AMASS,
     &                   HTICE(I), DEPB(I), PSURF(I),
     &                   SLOPN(1,I),X(I),Y(I),TEMP(I))
          ENDIF
          IF(TFLAG) TEMP(I)=TSORIG(I)
C
C ....... IF SURFACE TEMP GREATER THAN 0.0, SET IT TO ZERO...
          IF(TEMP(I).GT.0.0) TEMP(I)=0.0D0
C
C
          TEMPA(1,I)=TEMP(I)
          IF(HTICE(I).GT.SEALEV) THEN
            THICK=HTICE(I)-DEPB(I)
          ELSE
            THICK=0.D0
          ENDIF
C ....... PRESSURE MELTING POINT AS FUNCTION OF DEPTH (THICKNESS)
C ....... TMELT=-8.7E-4*THICK
          TMELT=PMP(THICK)
C ...............................................................
C ----- A SIMPLE 1D STEADY STATE TEMPERATURE SOLVER
          IF(THICK.LE.1.) THEN
            INOICE=INOICE+1
            IF(DEPB(I).LT.SEALEV) THEN
              TBOT=0.0D0
              TBED(I)=TBOT
              TSURF=TEMP(I)
              TDT=(TBOT-TSURF)/(MMAX-1)
              DO J=1,MMAX
C                TEMPA(J,I)=0.0
                TEMPA(J,I)=TBOT-(MMAX-J+1)*TDT
              ENDDO
            ELSE
              TBOT=TEMP(I)
              TBED(I)=TBOT
              DO J=1,MMAX
                TEMPA(J,I)=TBOT
              ENDDO
            ENDIF        
            TTEMP=MIN(TMELT,TEMP(I))
C ***************THIS IS THE OLD ICE HARDNESS FUNCTION*********
C            IF(TTEMP.GE.-10.) THEN
C              AEFF=A0*EXP(-B0*TTEMP)
C            ELSE
C              AEFF=A1*EXP(-B1*TTEMP)
C            ENDIF
C **********THIS IS THE EISMINT ICE HARDNESS FUNCTION**********
            AEFF=HARDNESS(TTEMP)
            CALL TYPETEST(ITYPE(I),FRACT(I),AFUDGE(I),FLOWA(I),HTICE(I),
     &                    DEPB(I),SLDGB(I),WTHICK(I),STHRESH,
     &                    TBOT,TMELT,IMELTED,IFROZEN,AEFF)
          ELSE  
            IWITHICE=IWITHICE+1
            IF(TEMP(I).GE.0.0 .OR. .NOT.ISICE(HTICE(I),DEPB(I))) THEN
              DO J=1,MMAX
                TEMPA(J,I)=TEMP(I)
              ENDDO
              TBOT=TEMP(I)
C              AEFF=A1
C ****************EISMINT STUFF*************************
               AEFF=HARDNESS(0.D0)
            ELSE
              CALL COLTEMP1(TIME,I,NUMNP,AMASS,ADOT(I),TEMPA(1,I),
     &                      THICK,TMELT,TBOT,AEFF,DT,IPLSTRT,IPLT,
     &                      SLOPN(1,I),IMELT,BMELT(I),WTHICK(I),
     &                      STHRESHB,GEOFLUX(I),FRACT(I),SLDGB(I),
     &                      ITYPE(I),CONTRIB(I))
              CALL MODBMELT(BMELT(I),THICK,SLOPN(1,I))
              BMELT(I) = (CONV*BMELTOLD(I) + BMELT(I))/(CONV+1.D0)
              BMELTOLD(I)=BMELT(I)
C             IFLAG=1
C             BMELT(I)=-0.001
            ENDIF
            TBAVG=TBAVG+TBOT
            NTAVG=NTAVG+1
            CALL TYPETEST(ITYPE(I),FRACT(I),AFUDGE(I),FLOWA(I),HTICE(I),
     &                    DEPB(I),SLDGB(I),WTHICK(I),STHRESH,
     &                    TBOT,TMELT,IMELTED,IFROZEN,AEFF)
C -----------------
            TMIN=MIN(TMIN,TBOT)
            TMAX=MAX(TMAX,TBOT)
            IF(TBOT.GE.0.0) THEN
              IC1=IC1+1
            ELSE
              IC2=IC2+1
            ENDIF
            TBED(I)=TBOT
          ENDIF
        ELSE
          IFIXED=IFIXED+1
        ENDIF
C ... OUTPUT FOR .STF FILE
C        IF(IFLAG.EQ.0) THEN
C          DO J=1,MMAX
C            WRITE(66,*) REAL(X(I)),REAL(Y(I)),REAL(ZZZ(J)),
C     &                  'MISSING'
C          ENDDO
C        ELSE
C          DO J=1,MMAX
C            WRITE(66,*) REAL(X(I)),REAL(Y(I)),REAL(ZZZ(J)),
C     &                 REAL(TEMPA(J,I))
C          ENDDO
C        ENDIF
      ENDDO
      IF(KSTART.EQ.0) THEN
        KSTART=1
        HTOGG=HSAVE
      ENDIF
      IF(NTAVG.GT.0) THEN
        TBAVG=TBAVG/NTAVG
      ENDIF
      IF(ISTART.EQ.0) THEN
        DT=DTSAVE
        ISTART=1
      ENDIF
      PFIXED=DBLE(100*IFIXED)/NUMNP
      PNOICE=DBLE(100*INOICE)/NUMNP
      PWITHICE=DBLE(100*IWITHICE)/NUMNP
      PMELTED=DBLE(100*IC1)/NUMNP
      PFROZEN=DBLE(100*IC2)/NUMNP
      PSLIDE=DBLE(100*IMELTED)/NUMNP
      PNOSLIDE=DBLE(100*IFROZEN)/NUMNP
      AWI=PWITHICE*GRIDAREA/100.D0
      AWS=PMELTED*GRIDAREA/100.D0
      IF(AWI.NE.0.) THEN
        PWS=(AWS/AWI)*100.D0
      ELSE
        PWS=0.0D0
      ENDIF
      IF(IOTOGG) THEN
        WRITE(LIST(IPAGE+1),9) 'NO TEMP CALC',IFIXED,PFIXED
        WRITE(LIST(IPAGE+2),9) 'NO ICE      ',INOICE,PNOICE
        WRITE(LIST(IPAGE+3),9) 'WITH ICE   ',IWITHICE,PWITHICE,
     &                           AWI
        WRITE(LIST(IPAGE+4),9) ' MELTED     ',IC1,PMELTED,
     &                           AWS,PWS
        WRITE(LIST(IPAGE+5),9) ' FROZEN     ',IC2,PFROZEN
        WRITE(LIST(IPAGE+6),9) ' SLIDING    ',IMELTED,PSLIDE
        WRITE(LIST(IPAGE+7),9) ' NO SLIDING ',IFROZEN,PNOSLIDE
9     FORMAT(1X,A12,I6,F9.4,1PE13.6,0PF9.4)
        WRITE(LIST(IPAGE+8),*) REAL(TMIN),REAL(TMAX),REAL(TBAVG)
        IPAGE=IPAGE+8
        WRITE(79,*) 'TIME=',TIME
        WRITE(79,9) 'NO TEMP CALC',IFIXED,PFIXED
        WRITE(79,9) 'NO ICE      ',INOICE,PNOICE
        WRITE(79,9) 'WITH ICE   ',IWITHICE,PWITHICE,
     &                           AWI
        WRITE(79,9) ' MELTED     ',IC1,PMELTED,
     &                           AWS,PWS
        WRITE(79,9) ' FROZEN     ',IC2,PFROZEN
        WRITE(79,9) ' SLIDING    ',IMELTED,PSLIDE
        WRITE(79,9) ' NO SLIDING ',IFROZEN,PNOSLIDE
        WRITE(79,*) 'TMIN,TMAX,TAVG',REAL(TMIN),REAL(TMAX),REAL(TBAVG)
      ENDIF
      ASDF=TMIN
C      CALL WRITEMP(NUMNP)
      IF(IC1+IC2.GT.0) THEN
        RMAF=DBLE(IC1)/DBLE(IC1+IC2)
      ELSE
        RMAF=0.0D0
      ENDIF
C ++++++++ EISMINT STUFF FOR SLIDING EXPERIMENT +++++++++++++++++++++++++++++
C      CALL PAYNEOUT1(TIME,NUMNP,NUMEL,KX,X,Y,HTICE,DEPB,
C     &                    TEMPA,VEL,FLOWA,
C     &                    VOL,AREA,RMAF)
      IF(IPLSTRT.EQ.1) CALL GRSTOP1
      RETURN
      END
C====================================================
      SUBROUTINE BACKSUB(N,DDD,CCC,BBB,XXX)                             
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION BBB(N),CCC(N),DDD(N),XXX(N)                      
      XXX(N)=BBB(N)/DDD(N)                                              
      N1=N-1                                                            
      DO II=1,N1                                                      
        I=N-II                                                            
        XXX(I)=(BBB(I)-CCC(I)*XXX(I+1))/DDD(I)                          
      ENDDO                                                          
      RETURN                                                            
      END                                                               
C====================================================
C THIS IS THE NEW VERSION .................
      SUBROUTINE COLTEMP1(TIME,IPT,NUMNP,AMASS,ADOT,TTTT,
     &                    THICK,TMELT,TBOT,AEFF,DELT,IPLSTRT,IPLT,
     &                    SLOPN,IMELT,BMELT,WTHICK,STHRESHB,SIGMALI,
     &                    FRACT,SLDGB,ITYPE,CONTRIB)
C******************************************
C THIS SOLVES FOR TIME DEPENDENT VERTICAL TEMPERATURE PROFILE
C WITH SOLVED FOR BASAL TEMPERATURE
C IF BASAL TEMP > PRESSURE MELTING, WILL
C RESOLVE WITH BASAL BC FIXED AT PRESSURE MELTING
C AND CALCULATE BASAL MELTING
C******************************************
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      PARAMETER(MMAX=40,BIG=1.D20)
      DIMENSION AMASS(11),DTDT(MMAX)
      DIMENSION XXX(MMAX),TTTT(MMAX),COND(MMAX),HEAT(MMAX),RHO(MMAX)
      DIMENSION QQQ(MMAX),SOURCE(MMAX),DELX(MMAX)
      DIMENSION SAVED(MMAX)
      DIMENSION RKX(MMAX),FFFF(MMAX),AAAD(MMAX),AAAU(MMAX),AAAL(MMAX)
      DIMENSION CCCD(MMAX),CCCU(MMAX),CCCL(MMAX)
      DIMENSION WWWW(MMAX),TOLD(MMAX),AT(MMAX),UC(MMAX),DDATA(9)
      DIMENSION U(MMAX),TTTT0(MMAX)
      DIMENSION LMAP(MMAX)
      DIMENSION RJUNK(MMAX),RHOC(MMAX),TNEW(MMAX)
      COMMON /HEAT/ HTOGG
      LOGICAL HTOGG
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      logical stopflag
      common /stop/ stopflag
      SIGMAL=-SIGMALI
      SIGADJ=-SIGMALI
C      SIGMAL=-AMASS(11)
      PG = 0.089866D0*0.3816d0
      THIRD=1.D0/3.D0
      IF(.FALSE.) THEN
        PRINT *,'IN COLTEMP1 BEFORE'
        PRINT *,' IPT = ',IPT
        PRINT *,' NUMNP = ',NUMNP
        PRINT *,' AMASS = ',AMASS(11)
        PRINT *,' TTTT = ',TTTT(1)
        PRINT *,' THICK = ',THICK
        PRINT *,' TBOT = ',TBOT
        PRINT *,' AEFF = ',AEFF
        PRINT *,' DELT = ',DELT
        PRINT *,' IPLSTRT = ',IPLSTRT
        PRINT *,' IPLT = ',IPLT
        PRINT *,' WTHICK = ',WTHICK
      ENDIF
      NPT=MMAX
      NDT=10
      WWWW(1)=0.D0
C ... QQQ IS FOR HORIZONTAL ADVECTION (MULT BY U (HORIZONTAL VELOCITY AT
C ... EACH DEPTH)
C ... SOURCE IS FOR INTERNAL (SHEAR) HEAT GENERATION
      QQQS=0.0D0
      WWWW(1)=ADOT
C
      LTOT=0
      DO I=1,NPT/2
        LMAP(I)=I
c       LMAP(I)=1
        LTOT=LTOT+LMAP(I)
      ENDDO
      DO I=NPT/2+1,NPT
        LMAP(I)=LMAP(I-1)-1
c       LMAP(I)=1
        LTOT=LTOT+LMAP(I)
      ENDDO
      XXX(1)=THICK
      XXX(NPT)=0.D0
      DX=(XXX(NPT)-XXX(1))/DBLE(LTOT)
      QQQ(1)=0.0D0
      TOLD(1)=TTTT(1)
      DO I=2,NPT
        XXX(I)=XXX(I-1)+LMAP(I-1)*DX
        QQQ(I)=0.D0
        TOLD(I)=TTTT(I)
      ENDDO
C
      TOLER=.001D0
      FACTOR=1.D0
      XMAX=-1D30
      XMIN=1D30
9     CONTINUE
      DO I=1,NPT
        XMAX=MAX(XMAX,XXX(I))
        XMIN=MIN(XMIN,XXX(I))
        DEPTH=XXX(1)-XXX(I)
        RHO(I)=DENSITY(DEPTH)
        COND(I)=CONDUCT(TTTT(I),RHO(I))
        HEAT(I)=SPHEAT(TTTT(I))
C*************USING FIXED VALUES OF RHO,COND,HEAT*****************
C         RHO(I)=DENSITY(1000.D0)
C         COND(I)=CONDUCT(-30.D0,RHO(I))
C         HEAT(I)=SPHEAT(-30.D0)
C*************USING FIXED VALUES OF RHO,COND,HEAT*****************
C ..... EVERYTHING IN ONE TERM...
C        RKX(I)=COND(I)/RHO(I)/HEAT(I)
C        RHOC(I)=1D0
C ..... ALA PATERSON PG 224
        RKX(I)=COND(I)
        RHOC(I)=RHO(I)*HEAT(I)
        IF(I.GT.1) DELX(I-1)=-(XXX(I)-XXX(I-1))
        AAAD(I)=0.D0
        AAAU(I)=0.D0
        AAAL(I)=0.D0
        CCCD(I)=0.D0
        CCCU(I)=0.D0
        CCCL(I)=0.D0
        FFFF(I)=0.D0
        U(I)=0.0D0
        SOURCE(I)=0.0D0
        SAVED(I)=0.0D0
      ENDDO
C           
C          
      IF(IPLSTRT.EQ.0 .AND. IPLT.EQ.1) THEN
        CALL GRSTRT(600,600)
C        CALL WINDOW(REAL(TMIN),REAL(TMAX),REAL(XMIN),REAL(XMAX))
        CALL WINDOW(REAL(TMIN),REAL(TMAX),0.,5000.)
        IPLSTRT=1
      ENDIF
      TMIN=-60.D0
      TMAX=60.D0
      IF(IPLT.EQ.1) THEN
C        CALL WINDOW(REAL(TMIN),REAL(TMAX),REAL(XMIN),REAL(XMAX))
C        CALL WINDOW(REAL(TMIN),REAL(TMAX),0.,5000.)
C        CALL NEWPAG
        CALL LINCLR(1)
        CALL MOVE(0.,REAL(XMIN))
        CALL DRAW(0.,REAL(XMAX))
        CALL MOVE(REAL(TMELT),REAL(XMIN))
        CALL MOVE(REAL(TMELT),REAL(XMIN))
        CALL DRAW(REAL(TMELT),REAL(XMAX))
        CALL MRKCLR(2)
      ENDIF
      WWWW(NPT)=0.D0
C --- VERTICAL VELOCITY DISTRIBUTION: LINEAR
      SLOPE=(WWWW(NPT)-WWWW(1))/(XXX(NPT)-XXX(1))
      DO I=2,NPT-1
        WWWW(I)=WWWW(1)+SLOPE*(XXX(I)-XXX(1))
      ENDDO
C      DO I=1,NPT
C        QQQ(I)=QQQS
C      ENDDO
      DO I=1,NPT                                                    
        TTTT0(I)=TTTT(I)                                                
      ENDDO                                                          
C.....BUILD CAPACITANCE MATRIX                                              
      CCCD(1)=THIRD*DELX(1)*RHOC(1)                                            
      CCCD(NPT)=THIRD*DELX(NPT-1)*RHOC(NPT-1)                                       
      DO I=2,NPT-1                                                  
        CCCD(I)=THIRD*(DELX(I)*RHOC(I)+DELX(I-1)*RHOC(I-1))                               
      ENDDO                                                          
      DO I=1,NPT-1                                                  
        CCCU(I)=.5D0*THIRD*DELX(I)*RHOC(I)
        CCCL(I)=CCCU(I)                                                 
      ENDDO                                                          
C.....END CAPACITANCE MATRIX                                                
C --- LOOP ON ITERATION STEP    
      IF(DELT.GT.0.) THEN
        NDT=1
      ELSE
        NDT=5
      ENDIF
      DO 900 IDT=1,NDT 
C        PRINT *,'ITERATION=',IDT
        IF(IPLT.EQ.1) THEN
C          CALL WINDOW(REAL(TMIN),REAL(TMAX),REAL(XMIN),REAL(XMAX))
C          CALL WINDOW(REAL(TMIN),REAL(TMAX),0.,5000.)
        ENDIF
C.....COME IN HERE FOR ITERATION ON VERTICAL VELOCITY
700   CONTINUE   
C
C ... HEAT GENERERATION SECTION ...................................
C
C ..... CALCULATE THE INTERNAL HEAT GENERATION DUE TO SHEAR
        DO I=1,NPT
          SOURCE(I)=0.D0
          SAVED(I)=0.D0
        ENDDO
        UFLOW=0.0
        DO I=2,NPT
          TMID=0.5D0*(TTTT(I-1)+TTTT(I))
c         TMID=0.5D0*(TTTT(1)+TTTT(NPT))
          ATI=HARDNESS(TMID)
          DEPTH=XXX(1)-0.5D0*(XXX(I-1)+XXX(I))
          DLX=(XXX(I-1)-XXX(I))
C ....... PATERSON IS 2*EPSILON-DOT*TAU RATE OF HEAT/UNIT VOLUME (TIMES ELEMENT
C         THICKNESS TO GET TOTAL HEAT DEPOSITED IN ELEMENT)
C         FOLLOWING IS 1/2 THAT, TO BE DEPOSITED AT TOP AND BOTTOM OF ELEMENT
C ....... FOLLOWING IS IN BAR/YR
          HTFLOW=DLX*((PG*DEPTH*SLOPN)**4/ATI**3)
          UFLOW=UFLOW+DLX*(PG*DEPTH*SLOPN/ATI)**3
C ....... FOLLOWING CONVERTS TO CAL/YR
          HTFLOW=HTFLOW*2.392D4
C ....... FOLLOWING TURNS DOWN HEAT WHERE SLIDING IS OCCURING
C          HTFLOW=HTFLOW*(1.D0-FRACT)
C
          SOURCE(I-1)=SOURCE(I-1)+HTFLOW
          SOURCE(I)=SOURCE(I)+HTFLOW
          SAVED(I-1)=SAVED(I-1)+HTFLOW
          SAVED(I)=SAVED(I)+HTFLOW
        ENDDO
        do i=1,npt
          saved(i)=saved(i)/(xxx(1)-xxx(npt))
          saved(i)=log10(saved(i))
         enddo
C
C ..... CALCULATE THE INTERNAL HEAT GENERATION DUE TO BASAL SLIDING
        HTSLID=(PG*THICK*SLOPN)**3/SLDGB**2
C ..... FOLLOWING CONVERTS TO CAL/YR
        HTSLID=HTSLID*2.392D4
C ..... SCALE BY WTHICK WITH RLUB FUNCTION (EXPERIMENTAL)
        IF(ITYPE.EQ.8) THEN
          HTSLID=HTSLID*RLUB(WTHICK)
          USLIDE=(PG*THICK*SLOPN/SLDGB)**2*RLUB(WTHICK)
          IF(USLIDE+UFLOW.GT.0.0) THEN
            FRACT=USLIDE/(USLIDE+UFLOW)
          ELSE
            FRACT=0
          ENDIF
C         IF(USLIDE.GT.0) PRINT *,REAL(USLIDE),REAL(UFLOW),REAL(FRACT)
        ENDIF
C ..... REDUCE BY 10 (ALA RAYMOND...)
        HTSLID=HTSLID/10.D0
C
C ..... THIS DEPOSITS ALL THE HEAT AT THE BED RATHER THAN WITHIN THE COLUMN
        HTFLOW=0.D0
        IF(.true.) THEN
          DO I=1,NPT
            HTFLOW=HTFLOW+SOURCE(I)
            SOURCE(I)=0.0D0
          ENDDO 
        ENDIF
C
        IF(.FALSE. .AND. WTHICK.GT.0) PRINT *,REAL(-CONTRIB/SIGMAL),
     &          REAL(-HTSLID/SIGMAL),
     &          REAL(FRACT),
     &          REAL(WTHICK)
C
        IF(ITYPE.NE.8) THEN
C ....... USE HEAT PROPORTIONAL TO FRACT ...
          CONTRIB=(1.D0-FRACT)*HTFLOW+FRACT*HTSLID
        ELSEIF(ITYPE.EQ.8) THEN
C ....... USE HEAT SUM OF FLOW AND SLIDING ...
          CONTRIB=HTFLOW+HTSLID
        ELSEIF(.FALSE.) THEN
C ....... USE HEAT FROM FLOW PROPORTIONAL TO FRACT ...
          CONTRIB=(1.D0-FRACT)*HTFLOW
        ENDIF
        IF(HTOGG) THEN
          SIGADJ=SIGMAL-CONTRIB
        ELSE
C ....... THIS TURNS OFF INTERNAL HEAT GENERATION ...........
          SIGADJ=SIGMAL
          HTSLID=0
          HTFLOW=0
          CONTRIB=0
          DO I=1,NPT
            SOURCE(I)=0.0D0
          ENDDO 
        ENDIF
C ........................................................................
C          PRINT *,REAL(CONTRIB),REAL(SLOPN),
C     &            REAL(XXX(1)-XXX(NPT)),
C     &            REAL(PG*SLOPN*(XXX(1)-XXX(NPT)))
C          PRINT *,'CONTRIB',CONTRIB,SIGMAL
C        DO I=1,NPT
C          PRINT *,I,SOURCE(I)
C        ENDDO
C ... END INTERNAL HEAT GENERATION SECTION
C ----- FORM LOAD VECTOR 
        FFFF(1)=.5D0*U(1)*QQQ(1)*DELX(1)+SOURCE(1)                             
        FFFF(NPT)=.5D0*U(NPT-1)*QQQ(NPT-1)*DELX(NPT-1)+
     &             SOURCE(NPT)-SIGADJ          
        P1=FFFF(1)                                                      
        PN=FFFF(NPT)
        DO I=2,NPT-1    
          FFFF(I)=.5D0*(U(I)*QQQ(I)*DELX(I)+
     &                U(I-1)*QQQ(I-1)*DELX(I-1))+SOURCE(I)      
        ENDDO  
C         DO I=1,NPT
C           RJUNK(I)=FFFF(I)
C           PRINT *,'RJUNK',I,RJUNK(I)
C         ENDDO
C
C ----- END FORM LOAD VECTOR  
C                                                
C ----- FORM STIFFNESS MATRIX                                                 
C ... OLD FORM
C        AAAD(1)=RKX(1)/DELX(1)+WWWW(1)*.5                               
C        AAAD(NPT)=RKX(NPT-1)/DELX(NPT-1)-WWWW(NPT)*.5                   
C        DO I=2,NPT-1                                                
C          AAAD(I)=RKX(I-1)/DELX(I-1)+RKX(I)/DELX(I)                     
C          AAAD(I)=AAAD(I)+.5*(WWWW(I-1)-WWWW(I))                        
C        ENDDO                                                        
C        DO I=1,NPT-1                                                
C          AAAU(I)=-RKX(I)/DELX(I)                                       
C          AAAL(I)=AAAU(I)-WWWW(I)*.5                                    
C          AAAU(I)=AAAU(I)+WWWW(I)*.5                                    
C        ENDDO                                                        
C ...   NEW FORM
        AAAD(1)=RKX(1)/DELX(1)-RHOC(1)*WWWW(1)*.5D0                               
        AAAD(NPT)=RKX(NPT-1)/DELX(NPT-1)+RHOC(NPT-1)*WWWW(NPT)*.5D0                   
        DO I=2,NPT-1   
          AAAD(I)=RKX(I-1)/DELX(I-1)+RKX(I)/DELX(I)                     
          AAAD(I)=AAAD(I)+.5D0*(RHOC(I-1)*WWWW(I-1)-RHOC(I)*WWWW(I))                        
        ENDDO   
        DO I=1,NPT-1                                                
          AAAU(I)=-RKX(I)/DELX(I)                                       
          AAAL(I)=AAAU(I)-RHOC(I)*WWWW(I)*.5D0                                   
          AAAU(I)=AAAU(I)+RHOC(I)*WWWW(I)*.5D0 
        ENDDO 
C SAVE VALUES BEFORE EFFECTIVE MATRIX IS FORMED
        FFFFNPT=FFFF(NPT)+SIGADJ                                                    
        AAADNPT=AAAD(NPT)
        AAALNPT=AAAL(NPT-1)
        CCCDNPT=CCCD(NPT)
        CCCLNPT=CCCL(NPT-1)
C        PRINT *,0D0,CCCD(1),CCCU(1)
C        DO I=2,NPT-1
C          PRINT *,CCCL(I-1),CCCD(I),CCCU(I)
C        ENDDO
C        PRINT *,CCCL(NPT-1),CCCD(NPT)
C        PAUSE
C        PRINT *,0D0,AAAD(1),AAAU(1)
C        DO I=2,NPT-1
C          PRINT *,AAAL(I-1),AAAD(I),AAAU(I)
C        ENDDO
C        PRINT *,AAAL(NPT-1),AAAD(NPT)
C        PAUSE
C.......FORM EFFECTIVE LOAD VECTOR                                            
        IF(DELT.GT.0.) THEN
C      PRINT *,'FORM EFFECTIVE LOAD VECTOR',DELT
C          FFFF(1)=(FFFF(1)*DELT+                                          
C     &           (CCCD(1))*TTTT(1)+                                       
C     &           (CCCU(1))*TTTT(2)) 
C          FFFF(NPT)=(FFFF(NPT)*DELT+                                      
C     &           (CCCD(NPT))*TTTT(NPT)+                                   
C     &           (CCCL(NPT-1))*TTTT(NPT-1))                               
C          DO I=2,NPT-1                                                
C            FFFF(I)=(FFFF(I)*DELT+                                        
C     &           (CCCD(I))*TTTT(I)   +                                    
C     &           (CCCL(I-1))*TTTT(I-1) +                                  
C     &           (CCCU(I))*TTTT(I+1) )  
          FFFF(1)=FFFF(1)+
     &                   (CCCD(1)*TTTT(1)+                                       
     &                    CCCU(1)*TTTT(2))/DELT
          FFFF(NPT)=FFFF(NPT)+
     &                    (CCCD(NPT)*TTTT(NPT)+                                   
     &                     CCCL(NPT-1)*TTTT(NPT-1))/DELT
          DO I=2,NPT-1                                                
            FFFF(I)=FFFF(I)+
     &                     (CCCD(I)*TTTT(I)     +                                    
     &                      CCCL(I-1)*TTTT(I-1) +                                  
     &                      CCCU(I)*TTTT(I+1))/DELT  
C      PRINT *,'T',I,REAL(TTTT(I)),REAL(TTTT(I-1)),REAL(TTTT(I+1))
C      PRINT *,'C',I,REAL(CCCD(I)),REAL(CCCL(I-1)),REAL(CCCU(I))                                 
          ENDDO    
C         FORM EFFECTIVE STIFFNESS MATRIX                                       
          DO I=1,NPT                                                  
C            AAAD(I)=DELT*AAAD(I)+CCCD(I)                                  
C            AAAU(I)=DELT*AAAU(I)+CCCU(I)                                  
C            AAAL(I)=DELT*AAAL(I)+CCCL(I)     
            AAAD(I)=AAAD(I)+CCCD(I)/DELT                                  
            AAAU(I)=AAAU(I)+CCCU(I)/DELT                                  
            AAAL(I)=AAAL(I)+CCCL(I)/DELT                                  
          ENDDO            
        ENDIF   
C SAVE VALUES AFTER EFFECTIVE MATRIX IS FORMED                                         
        IF(.FALSE.) THEN
          FFFFNPT=FFFF(NPT)+SIGMAL                                                    
          AAADNPT=AAAD(NPT)
          AAALNPT=AAAL(NPT-1)
        ENDIF
C
C       DO I=1,NPT                                                  
C         WRITE(7,*) 'EFF F',I,REAL(FFFF(I)),REAL(RJUNK(I))                                     
C       ENDDO  
C ----- FIX BOUNDARY CONDITION FOR NODE 1                                     
        RK11=AAAD(1)                                                    
        RK12=AAAU(1)                                                    
        RK21=AAAL(1)                                                    
        AAAD(1)=1.D0
        AAAU(1)=0.D0
        AAAL(1)=0.D0
        FFFF(1)=TTTT(1)                                                 
        FFFF(2)=FFFF(2)-RK21*FFFF(1)                                    
C ----- FIX BOUNDARY CONDITION FOR NODE NPT
        RKN1=AAAD(NPT)                                                  
        RKN2=AAAU(NPT-1)                                                
        RK2N=AAAL(NPT-1)    
        IPASS=1     
C ..... ISET=0 IS FROZEN BED
        ISET=0       
C                                
C ----- THE MATRIX SOLUTION
C ..... COPY TEMPERATURE INTO TNEW FOR SOLUTION ...
        DO KKK=1,NPT
          TNEW(KKK)=TTTT(KKK)
        ENDDO
C       PRINT *,' TTTT ',(REAL(TNEW(KKK)),KKK=1,NPT)
        CALL TRI(NPT,AAAL,AAAD,AAAU,FFFF,TNEW) 
C       PRINT *,' AFTER ',(REAL(TNEW(KKK)),KKK=1,NPT)
c       if(IPT.eq.1377) then
c         WRITE(LIST(IPAGE+1),*) IPT,xxx(1)-xxx(npt),tnew(npt)
c         IPAGE=IPAGE+1
c       endif
c       if(IPT.eq.1377 .and. mod(nint(TIME),10000).eq.0) then
c         do i=1,npt
c           write(2,*) saved(i),xxx(i)
c           write(3,*) tnew(i),xxx(i)
c         enddo
c         write(2,*) -99999.,0.
c         write(2,19) 'S',TIME
c         write(3,*) -99999.,0.
c         write(3,19) 'T',TIME
c       endif 
C       PAUSE
C
        IF(TNEW(NPT).LT.TMELT .AND. WTHICK.GT.STHRESHB) THEN
C ....... CHECK FOR BED BELOW MELTING POINT WITH WATER PRESENT ...
c         PRINT *,REAL(TNEW(NPT)),REAL(TMELT),' BELOW MELTING POINT ',
c    &           IPT,wthick,sthreshb
          DO KKK=1,NPT
            TNEW(KKK)=TTTT(KKK)
          ENDDO
          TNEW(NPT)=TMELT
C ....... ISET=1 MEANS BED IS FREEZING, AT PMP...
          ISET=1
          AAAD(NPT)=1.D0
          AAAU(NPT-1)=0.D0
          AAAL(NPT-1)=0.D0
          FFFF(NPT)=TNEW(NPT)                                             
          FFFF(NPT-1)=FFFF(NPT-1)-RKN2*FFFF(NPT)                          
          CALL TRI(NPT,AAAL,AAAD,AAAU,FFFF,TNEW) 
        ELSEIF(TNEW(NPT).GT.TMELT .AND. ISET.EQ.0) THEN
C ....... CHECK FOR BED ABOVE MELTING POINT ...
c         PRINT *,REAL(TNEW(NPT)),REAL(TMELT),' ABOVE MELTING POINT ',
c    &           IPT
c         WRITE(LIST(IPAGE+1),*) IPT,real(tnew(npt)),' above PMP'
c         IPAGE=IPAGE+1
          DO KKK=1,NPT
            TNEW(KKK)=TTTT(KKK)
          ENDDO
          TNEW(NPT)=TMELT
C ....... ISET=1 MEANS BED IS MELTING, AT PMP...
          ISET=1
          AAAD(NPT)=1.D0
          AAAU(NPT-1)=0.D0
          AAAL(NPT-1)=0.D0
          FFFF(NPT)=TNEW(NPT)                                             
          FFFF(NPT-1)=FFFF(NPT-1)-RKN2*FFFF(NPT)                          
c         open(1,file='plot-source.data')
c         do i=1,npt
c           write(1,*) saved(i),xxx(i)
c         enddo
c         write(1,*) -99999.,0.
c         write(1,19) 'source',real(slopn)
19    format(1x,a,1pg9.2)
c         close(1)
c         open(1,file='plot-temp.data')
c         do i=1,npt
c           write(1,*) tnew(i),xxx(i)
c           write(2,*) saved(i),xxx(i)
c           write(3,*) tnew(i),xxx(i)
c         enddo
c         write(1,*) -99999.,0.
c         write(1,19) 'temp before',real(slopn)
c         write(2,*) -99999.,0.
c         write(2,19) 'S',TIME
c         write(3,*) -99999.,0.
c         write(3,19) 'temp before',TIME
          CALL TRI(NPT,AAAL,AAAD,AAAU,FFFF,TNEW) 
c         do i=1,npt
c           write(1,*) tnew(i),xxx(i)
c           write(3,*) tnew(i),xxx(i)
c         enddo
c         write(1,*) -99999.,0.
c         write(1,19) 'temp after',real(slopn)
c         write(3,*) -99999.,0.
c         write(3,19) 'temp after',real(slopn)
c         close(1)
c         stopflag=.true.
        ENDIF
C ..... ACCEPT TNEW SOLUTION INTO TTTT ....
        DO KKK=1,NPT
          IF(DELT.GT.0) THEN
            DTDT(KKK)=(TNEW(KKK)-TTTT(KKK))/DELT
          ENDIF
          TTTT(KKK)=TNEW(KKK)
        ENDDO
C
C ..... EXPERIMENT TO "SMOOTH" TEMPS ...
C        CALL SMOOTH(TTTT,NPT)
C .....
C
C --------------------------
C ..... SCAN THE COLUMN FOR TEMPERATURES ABOVE THE PRESSURE MELTING
C ..... POINT
        ISCAN=0
        DO I=1,NPT-1
          DEPTH=XXX(1)-XXX(I)
          PMPLOC=PMP(DEPTH)
          IF(TTTT(I).GT.PMPLOC) THEN
            TTTT(I)=PMPLOC
            ISCAN=ISCAN+1
          ENDIF
C         TTTT(I)=MIN(TTTT(I),PMPLOC)
        ENDDO
        IF(.FALSE.) THEN
          IF(ISCAN.GT.0..AND.IOTOGG) THEN 
            WRITE(LIST(IPAGE+1),*) ISCAN,' POLYTHERMAL AT ',IPT
            IPAGE=IPAGE+1
          ENDIF
        ENDIF
C ..................................................................
        AAAD(1)=RK11                                                    
        AAAU(1)=RK12                                                    
        AAAL(1)=RK21                                                    
        AAAD(NPT)=RKN1                                                  
        AAAU(NPT-1)=RKN2                                                
        AAAL(NPT-1)=RK2N                                                
C ------FOLLOWING FOR RKX FUNCTION OF TEMPERATURE                             
        DO I=1,NPT-1                                                
C --------RKXI= FUNCTION OF TEMPERATURE                                
          DEPTH=XXX(1)-XXX(I)  
          RHO(I)=DENSITY(DEPTH)                                             
          COND(I)=CONDUCT(TTTT(I),RHO(I))
          HEAT(I)=SPHEAT(TTTT(I))
C********************************************
C         RHO(I)=DENSITY(1000.D0)
C         COND(I)=CONDUCT(-30.D0,RHO(I))
C         HEAT(I)=SPHEAT(-30.D0)
C********************************************
C ..... EVERYTHING IN ONE TERM...
          RKXI=COND(I)/RHO(I)/HEAT(I)                                     
C ..... ALA PATERSON PG 224
          RKXI=COND(I)
          RLATENT=8.D4*RHO(I)
          RKX(I)=((FACTOR-1.D0)*RKX(I)+RKXI)/FACTOR                      
        ENDDO                                                        
C        SIGMA0=P1-((RK11-CCCD(1))*TTTT(1)+
C     &         RK12*TTTT(2)-CCCU(1)*TTTT0(2))
C        SIGMA0=-SIGMA0/DELT                                             
C        SIGMALC=PN-((RKN1-CCCD(NPT))*TTTT(NPT)+
C     &         RKN2*TTTT(NPT-1)-CCCU(NPT-1)*TTTT0(NPT-1))   
C        DTEMP=(TTTT(NPT-1)-TTTT(NPT))/(XXX(NPT-1)-XXX(NPT))
C        SIGMABC=-RKX(NPT-1)*DTEMP
        IF(DELT.GT.0) THEN
C ....... THIS IS FULL TIME-DEPENDENT, INCLUDING CAP-MATRIX
          TERM1=-AAALNPT*TTTT(NPT-1)-AAADNPT*TTTT(NPT)
          TERM2=-CCCLNPT*DTDT(NPT-1)-CCCDNPT*DTDT(NPT)
          SIGMABC=FFFFNPT+TERM1+TERM2
        ELSE
C ....... FOLLOWING STEADY-STATE (DELT=0)
          TERM1=-AAALNPT*TTTT(NPT-1)-AAADNPT*TTTT(NPT)
          TERM2=0.D0
          SIGMABC=FFFFNPT+TERM1+TERM2
        ENDIF
C---------------------------------------
C WHAT IS THIS ????
C
C        SIGMALC=-SIGMALC/DELT                                             
C--------------------------------------------
C ..... ONLY CALC BMELT IF TOGGLE IS ON...
        IF(IMELT.EQ.1) THEN
          IF(ISET.EQ.0) THEN
            BMELT=0.0D0
          ELSE
C           IF(ABS((SIGADJ-SIGMABC)/SIGMAL).GT.1E-6) THEN
              BMELT=-(SIGADJ-SIGMABC)/RKX(NPT-1)
              BMELT=-(SIGADJ-SIGMABC)/RLATENT
C           ENDIF
          ENDIF
C ... LET IT PLAY, WHATEVER IT MAKES... MAKE MODIFICATIONS IN MODBMELT
C ....... DISABLED ITERATION ON VERTICAL VELOCITY BASED ON BASAL MELT RATE
          IF(.FALSE.) THEN
            IF(ABS(WWWW(NPT)-BMELT).GT.TOLER) THEN
C             PRINT *,IPT,REAL(WWWW(NPT)),REAL(BMELT),
C     &              ' GOING AROUND AGAIN'
              WWWW(NPT)=BMELT
C...........  LOOP BACK TO BEGINNING AFTER RECALCULATING VERTICAL VELOCITY
              SLOPE=(WWWW(NPT)-WWWW(1))/(XXX(NPT)-XXX(1))
              DO I=2,NPT-1
                WWWW(I)=WWWW(1)+SLOPE*(XXX(I)-XXX(1))
              ENDDO
              GOTO 700
            ENDIF
          ENDIF
C ....... DISABLEDTO HERE
        ENDIF
        DO I=1,NPT-1                                                
          DTDX=(TTTT(I+1)-TTTT(I))/DELX(I)                              
          FLUX=-RKX(I)*DTDX 
        ENDDO                        
        IF(IPLT.EQ.1) THEN
          CALL LINCLR(1)                                             
          CALL MOVE(REAL(TOLD(1)),REAL(XXX(1)))                           
          DO I=2,NPT                                                  
            CALL DRAW(REAL(TOLD(I)),REAL(XXX(I)))                         
          ENDDO 
          CALL LINCLR(2)                                             
          CALL MOVE(REAL(TTTT(1)),REAL(XXX(1)))                           
          DO I=2,NPT                                                  
            CALL DRAW(REAL(TTTT(I)),REAL(XXX(I)))                         
          ENDDO                                                        
        ENDIF
        DIFF=0D0
        DO I=1,NPT
          DIFF=DIFF+(TOLD(I)-TTTT(I))**2
          TOLD(I)=TTTT(I)
          TTTT0(I)=TTTT(I)
        ENDDO
C        PRINT *,'DIFF=',DIFF
        CALL VELO(TIME,NPT,5D-3,XXX,TTTT,RHO,AT,U,UC,
     &            BMELT,DDATA,AEFF)                                                         
900   CONTINUE            
C ... TBOT IS DIFFERENCE BETWEEN BASAL TEMP AND PRESSURE MELTING POINT
      TBOT=TTTT(NPT)-TMELT
C....................................................................
C      PRINT *,' TBOT = ',TBOT
C        CALL WAIT(10000)
C      IF(ADOT.LT.0) CALL WAIT(100000)
C ++++++++ EISMINT STUFF FOR SLIDING EXPERIMENT +++++++++++++++++++++++++++++
C      IF (IPT.EQ.1861) CALL PAYNEOUT2(TIME,IPT,MMAX,TTTT,XXX,AT,WWWW)
C      IF (IPT.EQ.2471) CALL PAYNEOUT2(TIME,IPT,MMAX,TTTT,XXX,AT,WWWW)
C      IF (IPT.EQ.2776) CALL PAYNEOUT2(TIME,IPT,MMAX,TTTT,XXX,AT,WWWW)
C      IF (IPT.EQ.2109) CALL PAYNEOUT2(TIME,IPT,MMAX,TTTT,XXX,AT,WWWW)
C      IF (IPT.EQ.2719) CALL PAYNEOUT2(TIME,IPT,MMAX,TTTT,XXX,AT,WWWW)
C ========================================================================
c      IF(IPT.EQ.228 .OR.IPT.EQ.227) THEN
c      IF(BMELT.gt.1.) THEN
      if(.false.) then
7     FORMAT(1X,A,1P3E15.6)
8     FORMAT(1X,A,3I13)
        write(7,*)
        WRITE(7,7) 'IN COLTEMP1 AFTER'
        WRITE(7,8) ' IPT     = ',IPT
        WRITE(7,7) ' NUMNP   = ',NUMNP
        WRITE(7,7) ' AMASS   = ',AMASS(11)
        WRITE(7,7) ' TTTT    = ',TTTT(1)
        WRITE(7,*) ' THICK   = ',THICK
        WRITE(7,7) ' TBOT    = ',(TBOT),(TTTT(NPT)),(TMELT)
        WRITE(7,*) ' AEFF    = ',AEFF
        WRITE(7,7) ' DELT    = ',DELT
        WRITE(7,8) ' IPLSTRT = ',IPLSTRT
        WRITE(7,8) ' IPLT    = ',IPLT
        WRITE(7,7) ' WTHICK  = ',WTHICK
        WRITE(7,7) ' BMELT   = ',BMELT
        WRITE(7,7) ' SIGMAL  = ',SIGMAL,SIGADJ
        WRITE(7,7) ' SIGMABC = ',SIGMABC,SIGMABC/SIGADJ,
     &                          (SIGADJ-SIGMABC)
        WRITE(7,7) ' CONTRIB = ',CONTRIB
        WRITE(7,7) ' HTSLID  = ',HTSLID,PG*THICK*SLOPN
        WRITE(7,*) ' HTFLOW  = ',HTFLOW
        WRITE(7,8) ' ITYPE   = ',ITYPE
        WRITE(7,7) ' FFFFNPT = ',(FFFFNPT),(TERM1),(TERM2)
        write(7,*) ' slope   = ',SLOPN
        write(7,*) ' temps   = '
        write(7,*)   TTTT    
        write(7,*) ' xxxs    = '
        write(7,*)   XXX    
        write(7,*)
        bmelt=0.d0
        pause
      ENDIF
      RETURN
C      IF(IPLT.EQ.1) THEN
C        CALL WAIT(10000)
C        CALL GRSTOP1
C      ENDIF
      END         
C*************************************************
      SUBROUTINE VELO(TIME,NPTS,DH,X,T,RHO,AT,U,UC,
     &                BMELT,DDATA,AEFF)
      PARAMETER(MMAX=40)
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE UOLD,UCOLD,UBAROLD,AVGOLD
      DIMENSION DDATA(9)
      DIMENSION X(NPTS),T(NPTS),RHO(NPTS),AT(NPTS)
      DIMENSION U(NPTS),UC(NPTS),UOLD(MMAX),UCOLD(MMAX)
      DATA UOLD /MMAX*0D0/, UCOLD /MMAX*0D0/, AVG /0D0/,UBAROLD /0D0/
      IF(NPTS.GT.MMAX) THEN
        PRINT *,'INCREASE MMAX IN VELO, MMAX=',MMAX
        STOP
      ENDIF
      THIRD=1D0/3D0
      A0=1.830357D0
      B0=7.738778D-2
      A1=2.633955D0
      B1=4.0990233D-2
      CALL HARD(A0,B0,A1,B1)
      DO I=1,NPTS
C ***************THIS IS THE OLD ICE HARDNESS FUNCTION*********
C        IF(T(I).GE.-10.) THEN
C          AT(I)=A0*EXP(-B0*T(I))
C        ELSE
C          AT(I)=A1*EXP(-B1*T(I))
C        ENDIF
C **********THIS IS THE EISMINT ICE HARDNESS FUNCTION**********
        AT(I)=HARDNESS(T(I))
C
      ENDDO
      TOP=X(1)
      G=9.81D-5*0.3816d0
      AVGAT=0.D0
      AVGT=0.D0
      DO I=NPTS-1,1,-1
        DX=X(I+1)-X(I)
        AVGAT=AVGAT-0.5D0*(AT(I+1)+AT(I))*DX
        AVGT=AVGT-0.5D0*(T(I+1)+T(I))*DX
      ENDDO
      AVGAT=AVGAT/(X(1)-X(NPTS))
      AVGT=AVGT/(X(1)-X(NPTS))
      DDATA(5)=AVGT
C      PRINT *,'AVERAGE A=',AVGAT
      AVGR=0.D0
      DO I=NPTS-1,1,-1
        DX=X(I+1)-X(I)
        AVGR=AVGR-0.5D0*(RHO(I+1)+RHO(I))*DX
      ENDDO
      AVGR=AVGR/(X(1)-X(NPTS))
C      PRINT *,'AVERAGE RHO=',AVGR
      UBAR=.4D0*(AVGR*G*DH/AVGAT)**3*(X(1)-X(NPTS))**4
C      PRINT *,'UBAR=',UBAR
      U(NPTS)=0.D0
      DO I=NPTS-1,1,-1
        SUM=U(I+1)
        V1=2.D0*(RHO(I+1)*G*DH*(TOP-X(I+1))/AT(I+1))**3
        V2=2.D0*(RHO(I)*G*DH*(TOP-X(I))/AT(I))**3
        DX=X(I+1)-X(I)
        U(I)=U(I+1)-0.5D0*(V1+V2)*DX
      ENDDO
      AVG=0.D0
      DO I=NPTS-1,1,-1
        DX=X(I+1)-X(I)
        AVG=AVG-0.5D0*(U(I+1)+U(I))*DX
      ENDDO
      AVG=AVG/(X(1)-X(NPTS))
      DDATA(6)=AVG
      AEFF=(.4D0*(AVGR*G*DH)**3*(X(1)-X(NPTS))**4/AVG)**THIRD
C      PRINT 1002,TIME,REAL(AVGT),REAL(AEFF),REAL(AVG),
C     &                REAL(T(NPTS)),REAL(BMELT)
1002  FORMAT(F10.0,2F10.3,F10.2,F10.3,F10.6)
      DDATA(7)=AEFF
C      PRINT *,'A-EFFECTIVE=',AEFF,3.17E-23/AEFF**3
      C=.5D0*(AVGR*G*DH/AEFF)**3
      THICK=X(1)-X(NPTS)
      DO I=1,NPTS
        UC(I)=C*(THICK**4-(TOP-X(I))**4)
      ENDDO
      YMIN=X(NPTS)
      YMAX=X(1)
      UMIN=1D30
      UMAX=-1D30
C       PRINT 1001,'X','T','A','U','UC'
1001  FORMAT(1X,5A13)
      DO I=1,NPTS
        UMIN=MIN(UMIN,U(I))
        UMAX=MAX(UMAX,U(I))
        UMIN=MIN(UMIN,UC(I))
        UMAX=MAX(UMAX,UC(I))
C         PRINT 1000,X(I),T(I),AT(I),U(I),UC(I)
1000    FORMAT(1X,1P5G13.6)
      ENDDO
C      PRINT *,UMIN,UMAX,YMIN,YMAX
C      CALL GRSTRT(400,400)
C      CALL WINDOW(REAL(0.),REAL(500.),REAL(YMIN),REAL(YMAX))
C      CALL LINCLR(0)
C      CALL MOVE(REAL(UOLD(1)),REAL(X(1)))
C      DO I=2,NPTS
C        CALL DRAW(REAL(UOLD(I)),REAL(X(I)))
C      ENDDO
C      CALL MOVE(REAL(AVGOLD),REAL(YMIN))
C      CALL DRAW(REAL(AVGOLD),REAL(YMAX))
C      CALL LINCLR(0)
C      CALL MOVE(REAL(UCOLD(1)),REAL(X(1)))
C      DO I=2,NPTS
C        CALL DRAW(REAL(UCOLD(I)),REAL(X(I)))
C      ENDDO
C      CALL MOVE(REAL(UBAROLD),REAL(YMIN))
C      CALL DRAW(REAL(UBAROLD),REAL(YMAX))
C      CALL LINCLR(1)
C      CALL MOVE(REAL(U(1)),REAL(X(1)))
C      DO I=2,NPTS
C        CALL DRAW(REAL(U(I)),REAL(X(I)))
C      ENDDO
C      CALL MOVE(REAL(AVG),REAL(YMIN))
C      CALL DRAW(REAL(AVG),REAL(YMAX))
C      CALL LINCLR(2)
C      CALL MOVE(REAL(UC(1)),REAL(X(1)))
C      DO I=2,NPTS
C        CALL DRAW(REAL(UC(I)),REAL(X(I)))
C      ENDDO
C      CALL MOVE(REAL(UBAR),REAL(YMIN))
C      CALL DRAW(REAL(UBAR),REAL(YMAX))
      AVGOLD=AVG
C      UBAROLD=UBAR
      DO I=1,NPTS
        UOLD(I)=U(I)
        UCOLD(I)=UC(I)
      ENDDO
C      CALL GRSTOP1
      END
C*****************************************
      SUBROUTINE HARD(B0P,SLPBP,BB0P,SLPBBP)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NNNN=11)
      SAVE IPASS,B0,SLPB,BB0,SLPBB
      save t,a,b
      REAL*8 T(NNNN),A(NNNN),B(NNNN)
c      DATA T /0.d0,-5.d0,-10.d0,-15.d0,-20.d0,-25.d0,-30.d0,
c    &      -35.d0,-40.d0,-45.d0,-50.d0/
c     DATA A /5.3D-15,1.7D-15,5.2D-16,3.1D-16,1.8D-16,
c    &        1.0D-16,5.4D-17,2.9D-17,1.5D-17,7.7D-18,3.8D-18/
      DATA IPASS /0/
      IF(IPASS.EQ.0) THEN
        IPASS=1
        T(1)=0.d0
        T(2)=-5.d0
        T(3)= -10.d0
        T(4)= -15.d0
        T(5)= -20.d0
        T(6)= -25.d0
        T(7)= -30.d0 
        t(8)=-35.d0
        t(9)= -40.d0
        t(10)= -45.d0
        t(11)= -50.d0 
c
        A(1)=5.3D-15
        A(2)=1.7D-15
        A(3)=5.2D-16
        A(4)=3.1D-16
        A(5)=1.8D-16
        A(6)=1.0D-16
        A(7)=5.4D-17
        A(8)=2.9D-17
        A(9)=1.5D-17
        A(10)=7.7D-18
        A(11)=3.8D-18
        CONVER=1.6D-2/5.2D-16
        DO I=1,NNNN
          TEMP=A(I)*CONVER
          B(I)=1.D0/TEMP
          B(I)=B(I)**(1.D0/3.D0)
C          PRINT *,T(I),A(I),TEMP,B(I)
        ENDDO
        B0=B(1)
        SLPB=-(LOG(B(1))-LOG(B(3)) )/(T(1)-T(3))
        SLPBB=-(LOG(B(3))-LOG(B(11)))/(T(3)-T(11))
        RLNBB0=LOG(B(3))+SLPBB*T(3)
        BB0=EXP(RLNBB0)
      ENDIF
C          PRINT *,'B0,BB0',B0,BB0
C          PRINT *,'R,RR  ',SLPB,SLPBB
      B0P=B0
      SLPBP=SLPB
      BB0P=BB0
      SLPBBP=SLPBB
      END
C*********************************************
      SUBROUTINE WAIT(I)
      DO N=1,I
        X=SIN(REAL(I))
      ENDDO
      END                                                      
C*************************************************
      SUBROUTINE TRI(N,AAA,DDD,CCC,BBB,XXX) 
      PARAMETER(MMAX=40)                            
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION AAA(N),BBB(N),CCC(N),DDD(N),XXX(N)                      
      DIMENSION AAL(MMAX),BBL(MMAX),CCL(MMAX),DDL(MMAX)
      IF(N.GT.MMAX) THEN
        PRINT *,'INCREASE MMAX IN SUB. TRI',MMAX
        STOP
      ENDIF
      DO I=1,N-1
        AAL(I)=AAA(I)
        BBL(I)=BBB(I)
C      PRINT *,I,BBL(I),BBB(I)
        CCL(I)=CCC(I)
        DDL(I)=DDD(I)
      ENDDO
      DDL(N)=DDD(N)
      BBL(N)=BBB(N)
      DO I=2,N                                                        
        XMULT=AAA(I-1)/DDL(I-1)                                         
C        PRINT *,' AI-1,DI-1,XM',
C     &          I,REAL(AAA(I-1)),REAL(DDL(I-1)),REAL(XMULT)
C        PRINT *,' DI,BI-1,CI-1 ',REAL(DDL(I)),REAL(BBL(I-1)),
C     &           REAL(CCL(I-1))
        DDL(I)=DDL(I)-XMULT*CCL(I-1)
C        PRINT *,'BI BEFORE ',I,BBL(I)
        BBL(I)=BBL(I)-XMULT*BBL(I-1)                                    
C        PRINT *,'BI AFTER  ',I,BBL(I)
      ENDDO                                                          
      XXX(N)=BBL(N)/DDL(N)                                              
      N1=N-1                                                            
      DO II=1,N1                                                      
        I=N-II                 
C        PRINT *,I,REAL(XXX(I+1)),REAL(DDL(I)),REAL(BBL(I)),REAL(CCL(I))                                           
        XXX(I)=(BBL(I)-CCL(I)*XXX(I+1))/DDL(I)                          
      ENDDO      
      RETURN                                                            
      END                                                               
C*************************************************
      FUNCTION DENSITY(DEPTH)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(RHOI=910.D0,RHOS=350.D0,CRHO=0.0374D0)
        IF (DEPTH .EQ. 0.) THEN
          DENSITY=RHOS
        ELSEIF ((DEPTH.GT.0.).AND.(DEPTH.LE.55.)) THEN
          DENSITY=RHOI-((RHOI-RHOS)*EXP(-CRHO*DEPTH))
        ELSEIF ((DEPTH.GT.55.) .AND.(DEPTH.LT.65.)) THEN
          DENSITY=7.D0*DEPTH+455.D0
        ELSEIF (DEPTH .GE. 65.) THEN
          DENSITY=RHOI
        ENDIF
      END
C*************************************************
      FUNCTION CONDUCT(TTTT,RHO)
      IMPLICIT REAL*8(A-H,O-Z)
C ... UNITS: CALORIES/M//YR/DEGREE
      IF(TTTT.GT.0.) THEN
        TTT=0.D0
      ELSE
        TTT=TTTT
      ENDIF
      RRR=RHO*.001D0
      CONDUCT=3.1536D06*
     &          (( 4.20102D0-0.00754145D0*TTT)
     &          -(26.48767D0-0.048779D0  *TTT)*RRR
     &          +(57.31865D0-0.141127D0  *TTT)*RRR**2
     &          -(29.55156D0-0.053672D0  *TTT)*RRR**3)
      END
C*************************************************
      FUNCTION SPHEAT(TTTT)
      IMPLICIT REAL*8(A-H,O-Z)
C ... UNITS: CALORIES/KG/DEGREE
      IF(TTTT.GT.0.) THEN
        TTT=0.D0
      ELSE
        TTT=TTTT
      ENDIF
      SPHEAT=1000.D0*(0.494D0+0.00138D0*TTT)                          
      END
C********************************************************
      SUBROUTINE WRITEMP(NUMNP)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      PARAMETER(MMAX=40)
      CHARACTER*1 CHAR
      COMMON /TEMPERATURE/ TEMP(MMAX,NMAX)
      PRINT *,'IN WRITEMP',NUMNP
      PRINT *,'   TO WRITE OUT BACKUP OF TEMPERATURE'
      PRINT *,'   INPUT Y'
      READ(*,1002) CHAR
      WRITE(99,1002) CHAR
1002  FORMAT(A1)
      IF(CHAR.EQ.'y' .OR. CHAR.EQ.'Y') THEN
        PRINT *,'0:BINARY, 1:ASCII'
        READ(*,*) IREAD
        WRITE(99,*) IREAD
        NPT=MMAX
        REWIND 44
        IF(IREAD.EQ.1) THEN
          DO I=1,NUMNP
            DO J=1,NPT
              WRITE(44,*) I,J,TEMP(J,I)
            ENDDO
          ENDDO
        ELSE
          WRITE(44) NUMNP,NPT
          WRITE(44) ((TEMP(J,I),J=1,NPT),I=1,NUMNP)
        ENDIF
      ENDIF
      END
C-------------------------------
      SUBROUTINE SMOOTH(X,N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N)
      DO I=2,N-1
        X(I)=(X(I-1)+X(I)+X(I+1))/3.D0
      ENDDO
      END
C-------------------------------
      FUNCTION HARDNESS(TTT)
      IMPLICIT REAL*8(A-H,O-Z)
c     PARAMETER(FUD=2**(1D0/3D0))
      data fud /1.25d0/
      DATA RRR /8.314D0/
      DATA Q1 /60000.D0/
      DATA A1 /1.14D-5/
      DATA Q2 /139000.D0/
      DATA A2 /5.47D10/
C *****THIS IS 3 FOR GREENLAND EISMINT EXPERIMENT, 1 FOR SLIDING
      DATA EH1 /3.D0/
C      DATA EH1 /1.D0/
C THIS IS 2**(1/3) WHICH FOR SOME REASON MULTIPLIES ENHANCEMENT FACTOR ?
      EH=EH1*FUD
C
      IF(TTT.LT.0.) THEN
        TEMP=TTT+273.15D0
      ELSE
        TEMP=273.15D0
      ENDIF
      IF(TEMP.LE.263.15D0) THEN
        AAA=A1*EXP(-Q1/RRR/TEMP)
      ELSE
        AAA=A2*EXP(-Q2/RRR/TEMP)
      ENDIF
      BBB=EH*AAA*1D15
      BBB=1.D0/BBB
      HARDNESS=BBB**(1.D0/3.D0)
      END
C==========================================
      SUBROUTINE COLUMNZ(NPT,HTICE,DEPB,ZZZ)
      PARAMETER(MMAX=40)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ZZZ(NPT),LMAP(MMAX)
      THICK=HTICE-DEPB
      LTOT=0
      DO I=1,NPT/2
        LMAP(I)=I
C        LMAP(I)=1
        LTOT=LTOT+LMAP(I)
      ENDDO
      DO I=NPT/2+1,NPT
        LMAP(I)=LMAP(I-1)-1
C        LMAP(I)=1
        LTOT=LTOT+LMAP(I)
      ENDDO
      ZZZ(1)=THICK
      ZZZ(NPT)=0.D0
      DX=(ZZZ(NPT)-ZZZ(1))/DBLE(LTOT)
      DO I=2,NPT
        ZZZ(I)=ZZZ(I-1)+LMAP(I-1)*DX
      ENDDO
      DO I=1,NPT
        ZZZ(I)=ZZZ(I)+DEPB
      ENDDO
      END
C==========================================
      FUNCTION PMP(THICK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(SLOPE=0.0074d0*0.089866d0*0.3816d0)
c     PMP=-8.7D-4*THICK
      PMP=-SLOPE*THICK
      END
C==========================================
      SUBROUTINE CHKFLOT(HTICEI,DEPBI,FRACTI)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
      DATA RHOI /0.917D0/, RHOW /1.092D0/
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      IF(DEPBI.LT.SEALEV) THEN
        FLOT=(1.D0-RATDEN)*(DEPBI-SEALEV)
        IF(HTICEI.LE.FLOT) FRACTI=1.D0
      ELSE
        IF(HTICEI.LE.DEPBI) FRACTI=1.D0        
      ENDIF
      END
C==========================================
      LOGICAL FUNCTION ISICE(HTICEI,DEPBI)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
      DATA RHOI /0.917D0/, RHOW /1.092D0/
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      ISICE=.TRUE.
      IF(DEPBI.LT.SEALEV) THEN
        FLOT=(1.D0-RATDEN)*(DEPBI-SEALEV)
        IF(HTICEI.LT.FLOT) ISICE=.FALSE.
      ELSE
        IF(HTICEI.LE.DEPBI) ISICE=.FALSE.
      ENDIF
      END
C-----------------------------------
      SUBROUTINE MODBMELT(BMELTI,THICK,SLOPN)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(.false.) RETURN
C
C ... MODIFY BMELT ......................
      IF(.FALSE.) THEN
C ..... ARBITRARY MULTIPLICATIVE FACTOR 
        BMELTI=BMELTI*0.01
      ENDIF
      IF(.TRUE.) THEN
C ..... THIS IS BASIC OFFSET .02 MM 
        BMELTI=BMELTI-2D-5
      ENDIF
      IF(.false. .AND. BMELTI.LT.0) THEN
C ..... REDUCE FREEZING BY FACTOR OF 10
C       PRINT *,' REDUCE FREEZING BY FACTOR OF 10 ',BMELTI
        BMELTI=BMELTI/10.D0
      ENDIF
      IF(.false. .AND. BMELTI.LT.0) THEN
C ..... LIMIT FREEZE RATE TO 10 CM/YR
C       PRINT *,' LIMIT FREEZE RATE TO 10 CM/YR ',BMELTI
        BMELTI=MAX(-0.1D0,BMELTI)
      ENDIF
      IF(.FALSE. .AND. BMELTI.LT.0) THEN
C ..... NO FREEZING AT ALL
C       PRINT *,' NO FREEZING AT ALL ',BMELTI
        BMELTI=0.D0
      ENDIF
      BLIMIT=0.01D0
      IF(.false. .AND. BMELTI.GT. BLIMIT) THEN
C ..... LIMIT MELT RATE TO 1 CM/YR
C       PRINT *,' LIMIT MELT RATE TO 1 CM/YR ',BMELTI
        BMELTI=MIN(BLIMIT,BMELTI)
      ENDIF
c      BLIMIT=9.1d0
      BLIMIT=0.1d0
      IF(.true. .AND. BMELTI.GT. BLIMIT) THEN
C ..... remove MELT RATEs over 10 CM/YR
c       PRINT *,' REMOVE MELT RATE OVER 10 CM/YR ',BMELTI
        BMELTI=0.d0
      ENDIF
      TLIMIT=10.D0
      IF(.FALSE. .AND. THICK.LT. TLIMIT) THEN
C ..... LIMIT BMELT TO ZERO IN THIN ICE REGIONS
        PRINT *,' MELT RATE TO 0 FOR THIN ICE',BMELTI,THICK
        BMELTI=0.D0
      ENDIF
      TLIMIT=10.D0
      SLIMIT=1E-2
      IF(.false. .AND. SLOPN.GT. SLIMIT .AND. THICK.LT.TLIMIT) THEN
C ..... LIMIT BMELT TO ZERO IN STEEP-THIN ICE REGIONS
        PRINT *,' MELT RATE 0 FOR STEEP/THINICE',BMELTI,SLOPN,THICK
        BMELTI=0.D0
      ENDIF
C ... MODIFY BMELT ......................
C
      END
C-----------------------------------
      SUBROUTINE TYPETEST(ITYPEI,FRACTI,AFUDGEI,FLOWAI,HTICEI,
     &                    DEPBI,SLDGBI,WTHICKI,STHRESH,
     &                    TBOT,TMELT,IMELTED,IFROZEN,AEFF)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL ISICE
      LOGICAL ISICE
C ... FOLLOWING IS AMOUNT OF WATER (M) NECESSARY FOR FULL SLIDING...
      DATA COUPLE /0.1D0/
      IF(ITYPEI.EQ.1) THEN
C ..... MODIFY FRACTION ONLY
        IF(TBOT.GE.TMELT) THEN
C         PRINT *,I,TBOT,TEMPI,THICK
          FRACTI=1.D0
          IMELTED=IMELTED+1
        ELSE
          FRACTI=0.D0
          IFROZEN=IFROZEN+1
        ENDIF
      ELSEIF(ITYPEI.EQ.2) THEN
C ..... MODIFY FLOW ONLY
        AEFF=AEFF*AFUDGEI
        FLOWAI=AEFF
      ELSEIF(ITYPEI.EQ.3) THEN
C ..... MODIFY FRACTION
        IF(TBOT.GE.TMELT) THEN
C         PRINT *,I,TBOT,TEMPI,THICK
          FRACTI=1.D0
          IMELTED=IMELTED+1
        ELSE
          FRACTI=0.D0
          CALL CHKFLOT(HTICEI,DEPBI,FRACTI)
          IF(FRACTI.EQ.0.) THEN
            IFROZEN=IFROZEN+1
          ELSE
            IMELTED=IMELTED+1
          ENDIF
        ENDIF
C ..... MODIFY FLOW
        AEFF=AEFF*AFUDGEI
        FLOWAI=AEFF
      ELSEIF(ITYPEI.EQ.4) THEN
C ..... MODIFY FLOW
        AEFF=AEFF*AFUDGEI
        FLOWAI=AEFF
C ..... MODIFY SLIDING
        SLDGBI=AEFF/200.D0
      ELSEIF(ITYPEI.EQ.5) THEN
C ..... MODIFY FRACTION
        IF(TBOT.GE.TMELT) THEN
C         PRINT *,I,TBOT,TEMPI,THICK
          FRACTI=1.D0
          IMELTED=IMELTED+1
        ELSE
          FRACTI=0.D0
          IFROZEN=IFROZEN+1
        ENDIF
C ..... MODIFY FLOW
        AEFF=AEFF*AFUDGEI
        FLOWAI=AEFF
C ..... MODIFY SLIDING
        SLDGBI=AEFF/200.D0
      ELSEIF(ITYPEI.EQ.6) THEN
C ..... MODIFY AFUDGE TO FIT FLOW
        AFUDGEI=FLOWAI/AEFF
      ELSEIF(ITYPEI.EQ.7) THEN
C ..... MODIFY FRACTION BY WATER THICKNESS
        IF(WTHICKI.GE.STHRESH .AND.
     &    ISICE(HTICEI,DEPBI)) THEN
C          FRACTI=MIN(1.D0,MAX(0.D0,(WTHICKI-STHRESH)/COUPLE))
          FRACTI=MIN(1.D0,MAX(0.D0,(WTHICKI)/COUPLE))
          IMELTED=IMELTED+1
        ELSE
          FRACTI=0.D0
          CALL CHKFLOT(HTICEI,DEPBI,FRACTI)
          IF(FRACTI.EQ.0.) THEN
            IFROZEN=IFROZEN+1
          ELSE
            IMELTED=IMELTED+1
          ENDIF
        ENDIF
C ..... MODIFY FLOW
        AEFF=AEFF*AFUDGEI
        FLOWAI=AEFF
C AS EXPERIMENT, MODIFY SLIDING TO MATCH FRACTION.
CC      IF(ISICE(HTICEI,DEPBI)) THEN
C         SLDGBI=0.133333*(1.5-FRACTI)
CC      ENDIF
      ELSEIF(ITYPEI.EQ.8) THEN
C ..... MODIFY FRACTION BY WATER THICKNESS
        IF(WTHICKI.GE.STHRESH .AND.
     &    ISICE(HTICEI,DEPBI)) THEN
          IMELTED=IMELTED+1
        ELSE
          FRACTI=0.D0
          CALL CHKFLOT(HTICEI,DEPBI,FRACTI)
C          IF(FRACTI.EQ.1) WTHICKI=1.
          IF(FRACTI.EQ.0.) THEN
            IFROZEN=IFROZEN+1
          ELSE
            IMELTED=IMELTED+1
          ENDIF
        ENDIF
C ..... MODIFY FLOW
        AEFF=AEFF*AFUDGEI
        FLOWAI=AEFF
      ENDIF
      END

     
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

      subroutine threed(ipass,npts,nelem,x,y,htice,bdrock,kx)
      include "parameter.h"
      parameter(nmax=MAXNUM)
      !include "fshort.h"
      real*8 x,y,htice,bdrock
      dimension x(nmax),y(nmax),htice(nmax),bdrock(nmax),kx(nmax,4)
      real xyzbed(3,nmax),xyzsrf(3,nmax),xyzgrd(3,nmax)
      integer icxyzb(nmax),icxyzs(nmax),kxp(nmax,4)
      common /plotdata/ numnp,numel,xyzbed,xyzsrf,xyzgrd,
     &                  kxp,icxyzb,icxyzs
      return
      end

      SUBROUTINE WMOVER0(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL, DT,
     &                 DTLOCAL,
     &                 ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &                 NZ, KZ, LM, TNEW,
     &                 BMELT, D, B, A, KA, ALPHAC, TOTALW, TOTALP,IPLOT)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NNMAX=MAXNUM)
      DIMENSION TESTVEC1(NNMAX),TESTVEC2(nnmax),TESTVEC3(nnmax)
      DIMENSION TESTVEC4(NNMAX),svec(nnmax),yyy(nnmax),phi(nnmax)
      DIMENSION ALPHAC(3)
      DIMENSION BMELT(NMAX),TNEW(NMAX)
      DIMENSION ITKODE(NMAX)
      DIMENSION KZ(NMAX)
      DIMENSION LM(5),HTICE(NMAX),DEPB(NMAX)
      DIMENSION KX(NMAX,4),X(NMAX),Y(NMAX),NTYPE(NMAX),D(NMAX),B(NMAX)
      REAL*8 A(NMAX,NZ),WTHICK(NMAX)
      INTEGER KA(NMAX,NZ+1)
      DIMENSION XI(2,9),ETA(2,9),W(2,9)
      REAL*4 TB(2)
c     REAL*4 TB(2),ETIME,DTIME
c     EXTERNAL ETIME,DTIME
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      SAVE NSTEP,ISTART,WSAVE,TSAVE,PSAVE
      DATA ISTART /0/, BIG /1D30/
      DATA RHOI /0.917D0/, RHOW /1.092D0/, GRAV /3.74d0/
      dtlocal=0.1
      NSTEP=dt/dtlocal
      TLOCAL=0
      do i=1,numnp
        phi(i)=rhoi*grav*(htice(i)-depb(i))+rhow*grav*depb(i)
      enddo
C ... CALL CHECKER THAT PUTS WATER ON/OFF ICE-FREE NODES
c ... (LAST ARG IS VALUE TO PUT ON ICE-FREE NODES, ZERO IT BEFORE
c     (turn on/off internally)
      CALL CHECKER(NMAX, NUMNP, NUMEL, NTYPE, KX, HTICE, DEPB, 
     &             WTHICK,0.0d0)
c ... form matrices
      H=DTLOCAL
      H2=0.5d0*H
      H6=H/6.d0
      do n=1,nstep
        TLOCAL=TLOCAL+DTLOCAL
c ..... save water thickness for convergence test ...
        do i=1,numnp
          svec(i)=wthick(i)
        enddo
c ..... FIRST RK STEP ...............................
        CALL FORMNT(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &             ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &             NZ, KZ, LM,
     &             BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &             AREAWET, AREADRY, AREATOT)
        IF(ISTART.EQ.0) THEN
          ISTART=1
          WSAVE=TOTALW
          IF(AREAWET.EQ.0.0) THEN
            TSAVE=0D0
          ELSE
            TSAVE=100.D0*TOTALW/AREAWET
          ENDIF
          IF(AREATOT.EQ.0.0) THEN
            PSAVE=0D0
          ELSE
            PSAVE=100.D0*AREAWET/AREATOT
          ENDIF
        ENDIF
        call formrhs(numnp,ka,kz,nz,a,b,d,phi,testvec1,rat1)
c ..... update water thickness for first RK step ........
        do i=1,numnp
          yyy(i)=max(0.d0,wthick(i)*(1-h2*alphac(3)))
          yyy(i)=yyy(i)+testvec1(i)*h2
          yyy(i)=max(0.d0,yyy(i))
        enddo
C ..... SECOND RK STEP .......................................
        CALL FORMNT(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &             ETA, XI, W, yyy, ITKODE, HTICE, DEPB,
     &             NZ, KZ, LM,
     &             BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &             AREAWET, AREADRY, AREATOT)
        call formrhs(numnp,ka,kz,nz,a,b,d,phi,testvec2,rat2)
c ..... update water thickness for SECOND RK step ........
        do i=1,numnp
          yyy(i)=max(0.d0,wthick(i)*(1-h2*alphac(3)))
          yyy(i)=yyy(i)+testvec2(i)*h2
          yyy(i)=max(0.d0,yyy(i))
        enddo
C ..... THIRD RK STEP .......................................
        CALL FORMNT(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &             ETA, XI, W, yyy, ITKODE, HTICE, DEPB,
     &             NZ, KZ, LM,
     &             BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &             AREAWET, AREADRY, AREATOT)
        call formrhs(numnp,ka,kz,nz,a,b,d,phi,testvec3,rat3)
c ..... update water thickness for THIRD RK step ........
        do i=1,numnp
          yyy(i)=max(0.d0,wthick(i)*(1-H*alphac(3)))
          yyy(i)=yyy(i)+testvec3(i)*H
          yyy(i)=max(0.d0,yyy(i))
        enddo
C ..... FOURTH RK STEP .......................................
        CALL FORMNT(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &             ETA, XI, W, yyy, ITKODE, HTICE, DEPB,
     &             NZ, KZ, LM,
     &             BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &             AREAWET, AREADRY, AREATOT)
        call formrhs(numnp,ka,kz,nz,a,b,d,phi,testvec4,rat4)
c ..... update water thickness for FOURTH RK step ........
        do i=1,numnp
          yyy(i)=max(0.d0,wthick(i)*(1-H*alphac(3)))
          yyy(i)=yyy(i)+
     &           h6*(testvec1(i)+
     &           2.d0*(testvec2(i)+testvec3(i))+testvec4(i))
          wthick(i)=max(0.d0,yyy(i))
        enddo
c ..... form DIFF for convergence test
        diff=0.d0
        do i=1,numnp
          diff=diff+abs(svec(i)-wthick(i))
        enddo
        diff=diff/numnp
        if(N.eq.1) then
          dsave=diff
        endif
        ratiod=diff/dsave
        if(mod(n,nstep/10).eq.0) then
          print *,n,real(diff),real(rat4),real(ratiod),
     &              REAL(1D-9*(TOTALW))
        endif
        if(ratiod.lt.1e-2) goto 10
      enddo
c      print *,'didnt converge to steady state'
10    continue
      RATIOW=(TOTALW-WSAVE)/DT
      TOTALT=100.D0*TOTALW/AREAWET
      TOTALP=100.D0*AREAWET/AREATOT
      RATIOT=(TOTALT-TSAVE)/DT
      RATIOP=(TOTALP-PSAVE)/DT
      IF(IOTOGG) THEN
        WRITE(LIST(IPAGE+1),*) ISTEP,
     &       '********************************************'
        WRITE(LIST(IPAGE+2),*) 
     &         'TOTAL WATER (KM**3)=',REAL(1D-9*TOTALW),
     &         REAL(1D-9*RATIOW),REAL(1D-9*RATIOW*DT)
C       WRITE(LIST(IPAGE+1),*) 'WET AREA    (KM**2)=',REAL(1D-6*AREAWET)
C       WRITE(LIST(IPAGE+1),*) 'DRY AREA    (KM**2)=',REAL(1D-6*AREADRY)
C       WRITE(LIST(IPAGE+1),*) 'TOTAL AREA  (KM**2)=',REAL(1D-6*AREATOT)
        WRITE(LIST(IPAGE+3),*) 'PERCENT WET        =',REAL(TOTALP),
     &         REAL(RATIOP),REAL(RATIOP*DT)
        WRITE(LIST(IPAGE+4),*) 'AVG THICK   ( CM  )=',REAL(TOTALT),
     &         REAL(RATIOT),REAL(RATIOT*DT)
        WRITE(LIST(IPAGE+5),*) 'used',n,' steps, ratio',real(ratiod),
     &                         real(TLOCAL)
        WRITE(LIST(IPAGE+6),*) 
     &       '********************************************'
        IPAGE=IPAGE+6
      ENDIF
      WSAVE=TOTALW
      PSAVE=TOTALP
      TSAVE=TOTALT
C ... CALL CHECKER THAT PUTS WATER ON/OFF ICE-FREE NODES
c ... (LAST ARG IS VALUE TO PUT ON ICE-FREE NODES, 0.1 IT AFTER
c     (turn on/off internally)
      CALL CHECKER(NMAX, NUMNP, NUMEL, NTYPE, KX, HTICE, DEPB, 
     &             WTHICK,1d-5)






      end
C-----------------------------------------------------------------
      SUBROUTINE FORMNT(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &                 ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &                 NZ, KZ, LM,
     &                 BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &                 AREAWET, AREADRY, AREATOT)
      IMPLICIT REAL*8(A-H,O-Z)
C FORM STIFFNESS MATRIX
      DIMENSION ALPHAC(3)
      DIMENSION BMELT(NMAX)
      DIMENSION ITKODE(NMAX)
      DIMENSION KZ(NMAX)
      DIMENSION LM(5),HTICE(NMAX),DEPB(NMAX)
      DIMENSION KX(NMAX,4),X(NMAX),Y(NMAX),NTYPE(NMAX),D(NMAX),B(NMAX)
      REAL*8 A(NMAX,NZ),WTHICK(NMAX)
      INTEGER KA(NMAX,NZ+1)
      DIMENSION P(5),S(5,5),DD(5)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2)
      DIMENSION PSI(4),DPSI(4,2)
      DIMENSION DWSIX(9),DWSIY(9)
      DIMENSION WSI(4),DWSI(4,2)
      DIMENSION XY(2,4),XI(2,9),ETA(2,9),W(2,9)
      DIMENSION ALPHA(2),BETA(2)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      LOGICAL UPWIND
      DATA BIG /1.D30/
      DATA ASCALE /1.D0/
C **********************************************************************
C     FORM CONDUCTIVITY MATRIX FOR COMPLETE BODY
C **********************************************************************
C
C ... ZERO OUT APPROPRIATE ARRAYS ...
      DO I=1,NUMNP
        KZ(I)=1
        KA(I,1)=I
        D(I)=0.0D0
        B(I)=0.0D0
        KA(I,NZ+1)=0
        DO J=1,NZ
          KA(I,J)=0
          A(I,J)=0.0D0
        ENDDO
      ENDDO
C
C ... BEGIN LOOP OVER ALL THE ELEMENTS ...
      TOTALW=0.0D0
      AREAWET=0.0D0
      AREADRY=0.0D0
      AREATOT=0.0D0
      avg=0.
      iavg=0
      DO 100 N=1,NUMEL
        IF(NTYPE(N).EQ.1) THEN
          NODEN=4
          NINT=9
        ELSE
          NODEN=3
          NINT=4
          PRINT *,'TRIANGLES NOT ALLOWED!'
          STOP
        ENDIF
        DO I=1,NODEN
          LM(I)=KX(N,I)
        ENDDO
C
C ..... FORM ELEMENT CONDUCTIVITY MATRIX
        DO I=1,NODEN
          DD(I)=0.0D0
          P(I)=0.0D0
          DO J=1,NODEN
            S(I,J)=0.0D0
          ENDDO
        ENDDO
C
        I=LM(1)
        J=LM(2)
        K=LM(3)
        L=LM(4)
        XY(1,1)=X(I)
        XY(1,2)=X(J)
        XY(1,3)=X(K)
        XY(2,1)=Y(I)
        XY(2,2)=Y(J)
        XY(2,3)=Y(K)
        IF(NODEN.EQ.4) THEN
          XY(1,4)=X(L)
          XY(2,4)=Y(L)
        ENDIF
C ..... USE FOLLOWING TO GENERATE CENTROID VALUES FOR GRADIENTS AND
C       OTHER MATERIAL PROPERTIES 
        IF(.TRUE.) THEN
          CALL FESHAPE(NTYPE(N),0D0,0D0,PSI,DPSI)
          CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
          TWTHIK=0.D0
          BMELTN=0.D0
          DGDX=0.D0
          DGDY=0.D0
          SURFX=0
          SURFY=0
          BEDX=0
          BEDY=0
          DO I=1,NODEN
            SURFX=SURFX+HTICE(LM(I))*DPSIX(I)
            SURFY=SURFY+HTICE(LM(I))*DPSIY(I)
          ENDDO
          SURFXY=SQRT(SURFX**2+SURFY**2)
          DO I=1,NODEN
            BEDX=BEDX+DEPB(LM(I))*DPSIX(I)
            BEDY=BEDY+DEPB(LM(I))*DPSIY(I)
            BMELTN=BMELTN+BMELT(LM(I))*PSI(I)
            TWTHIK=TWTHIK+WTHICK(LM(I))*PSI(I)
            IF(WTHICK(LM(I)).GT.0.001) THEN
C       ....  COEFFICIENT COMES FROM ALLEY'S PAPER
              RNPRES=5.0D-5*(HTICE(LM(I))-DEPB(LM(I)))*SURFXY/
     &              WTHICK(LM(I))
            ELSE
              RNPRES=0.D0
            ENDIF
c effective pressure turned off below...
            RNPRES=0.D0
            GGG=-(HTICE(LM(I))+0.09D0*DEPB(LM(I)))-RNPRES
            DGDX=DGDX+GGG*DPSIX(I)
            DGDY=DGDY+GGG*DPSIY(I)
          ENDDO
          DGXY=SQRT(DGDX**2+DGDY**2)
          SURFXY=SQRT(SURFX**2+SURFY**2)
          BEDXY=SQRT(BEDX**2+BEDY**2)
        ENDIF
C........... TO HERE
c
C ..... FORM ELEMENT MATRIX AND VECTORS
C
C ..... BEGIN INTEGRATION POINT LOOP
        DO L=1,NINT
          CALL FESHAPE(NTYPE(N),XI(NTYPE(N),L),ETA(NTYPE(N),L),PSI,DPSI)
C
C ......  GENERATE DEIRVATIVES, JUST NON-UPWINDING, COPY INTO UPWINDING
          CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
          DO I=1,NODEN
            WSI(I)=PSI(I)
            DWSI(I,1)=DPSI(I,1)
            DWSI(I,2)=DPSI(I,2)
            DWSIX(I)=DPSIX(I)
            DWSIY(I)=DPSIY(I)
          ENDDO
C ....... ACCUMULATE INTEGRATION POINT VALUES OF INTEGRALS
          FAC=DETJ*W(NTYPE(N),L)
          TOTALW=TOTALW+TWTHIK*FAC
          IF(TWTHIK.GT.0.0) THEN
            AREAWET=AREAWET+FAC
          ELSE
            AREADRY=AREADRY+FAC
          ENDIF
          AREATOT=AREATOT+FAC
C ....... LATEST NEW WAY ..................
          GRADG=SQRT(DGDX**2+DGDY**2)
          IF(TWTHIK.LE.0.01 .or. .true.) THEN
            IPP=2
            IQQ=1
            ACONST=ASCALE*ALPHAC(1)*TWTHIK**(IPP+1)*GRADG**(IQQ-1)
          ELSE
            RPP=.5D0
            RQQ=.7D0
            ACONST=ASCALE*ALPHAC(1)*TWTHIK**(RPP+1.d0)*GRADG**(RQQ-1)
          ENDIF
          avg=avg+aconst*gradg
          iavg=iavg+1
          DO I=1,NODEN
C ......... THIS IS LUMPED CAPACITANCE MATRIX
            DD(I)=DD(I)+WSI(I)*FAC
C ........
C ........  INTEGRATION POINT VALUES FROM FEM INTERPOLATIOON ....
C ......... THIS IS LOAD VECTOR
            P(I)=P(I)+BMELTN*WSI(I)*FAC
C
            DO J=1,NODEN
              TERM1=-ACONST*(dwsix(i)*dpsix(j)+dwsiy(i)*dpsiy(j))
              TERM2=0.0
              TERM3=0.0
              S(I,J)=S(I,J)+(TERM1+TERM2+TERM3)*FAC
            ENDDO
          ENDDO
        ENDDO
C       IF(S(1,1).NE.0) THEN
C         PRINT *,N
C         DO II=1,4
C         PRINT *,(REAL(S(II,JJ)),JJ=1,4)
C         ENDDO
C       ENDIF
C
C ..... ADD ELEMENT CONDUCTIVITY TO COMPLETE CONDUCTIVITY MATRIX
C
        DO L=1,NODEN
          I=LM(L)
C ....... THIS IS LUMPED CAPACITANCE MATRIX
          D(I)=D(I)+DD(L)
C
C ....... THIS IS LOAD VECTOR
          B(I)=B(I)+P(L)
          DO M=1,NODEN
            J=LM(M)
            IF(I.EQ.J) THEN
              A(I,1)=A(I,1)+S(L,M)
              KA(I,1)=I
            ELSE
              DO K=2,KZ(I)
                IF(KA(I,K).EQ.J) THEN
                  A(I,K)=A(I,K)+S(L,M)
                  GOTO 99
                ENDIF
              ENDDO
              KZ(I)=KZ(I)+1
              A(I,KZ(I))=S(L,M)
              KA(I,KZ(I))=J
            ENDIF
99          CONTINUE
          ENDDO
          KA(I,NZ+1)=KZ(I)
        ENDDO
100   CONTINUE
C ... END LOOP OVER ALL THE ELEMENTS ...
C      IF(NFIX.EQ.0) THEN
C        PRINT *,'NO NODES ARE FIXED....'
C        PAUSE
C      ELSE
C        PRINT *,NFIX,' NODES ARE FIXED....',BIG
C      ENDIF
C ... REMOVE HERE, THIS WRITES OUT MATRIX
C      REWIND 73
C      WRITE(73,*) NUMNP
C      DO I=1,NUMNP
C        WRITE(73,*) KA(I,10)
C        DO J=1,KA(I,10)
C          WRITE(73,*) I,KA(I,J),A(I,J)
C        ENDDO
C      ENDDO
C      WRITE(73,*) (B(I),I=1,NUMNP)
C ... TO HERE ******
C      IF(IOTOGG) THEN
C        WRITE(LIST(IPAGE+1),*) 
C     &     'TOTAL WATER (KM**3)=',REAL(1E-9*TOTALW)
C        IPAGE=IPAGE+1
C       PRINT *,'WET AREA    (KM**2)=',REAL(1E-6*AREAWET)
C       PRINT *,'DRY AREA    (KM**2)=',REAL(1E-6*AREADRY)
C       PRINT *,'TOTAL AREA  (KM**2)=',REAL(1E-6*AREATOT)
C       PRINT *,'PERCENT WET        =',REAL(100.*AREAWET/AREATOT)
C       PRINT *,'AVG THICK   ( CM  )=',REAL(100.*TOTALW/AREAWET)
C      ENDIF
      END


      subroutine formrhs(numnp,ka,kz,nz,a,b,d,phi,testvec,rat)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      DIMENSION TESTVEC(NMAX)
      DIMENSION KZ(NMAX)
      DIMENSION phi(NMAX)
      DIMENSION KX(NMAX,4),D(NMAX),B(NMAX)
      REAL*8 A(NMAX,NZ)
      INTEGER KA(NMAX,NZ+1)
c ..... stiffness matrix * pressure potential
      do i=1,numnp
        testvec(i)=0.d0
        do jj=1,kz(i)
          j=ka(i,jj)
          testvec(i)=testvec(i)+a(i,jj)*phi(j)
        enddo
      enddo
c ... add load vector
      rat=0
      irat=0
      do i=1,numnp
        if(b(i).gt.0) then
          rat=rat+testvec(i)/b(i)
          irat=irat+1
        endif
        testvec(i)=testvec(i)+b(i)
      enddo
      if(irat.gt.0) then
        rat=rat/irat
      endif
c ... scale with 1/d (M^-1)
      do i=1,numnp
        testvec(i)=testvec(i)/d(i)
      enddo
      end
      SUBROUTINE WMOVER(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL, DT,
     &                 DTLOCAL,
     &                 ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &                 NZ, KZ, LM, TNEW,
     &                 BMELT, D, B, A, KA, ALPHAC, TOTALW, TOTALP,IPLOT)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL LDIAG
      DIMENSION ALPHAC(3)
      DIMENSION BMELT(NMAX),TNEW(NMAX)
      DIMENSION ITKODE(NMAX)
      DIMENSION KZ(NMAX)
      DIMENSION LM(5),HTICE(NMAX),DEPB(NMAX)
      DIMENSION KX(NMAX,4),X(NMAX),Y(NMAX),NTYPE(NMAX),D(NMAX),B(NMAX)
      REAL*8 A(NMAX,NZ),WTHICK(NMAX)
      INTEGER KA(NMAX,NZ+1)
      DIMENSION XI(2,9),ETA(2,9),W(2,9)
      REAL*4 TB(2)
c     REAL*4 TB(2),ETIME,DTIME
c     EXTERNAL ETIME,DTIME
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
123   FORMAT(A25,T30,1PG13.6,G13.6)
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      SAVE NSTEP,ISTART,WSAVE,TSAVE,PSAVE
      DATA ISTART /0/, BIG /1D30/
C      IF(DTLOCAL.EQ.DT) NSTEP=1
      DTLOCAL=DT/10.d0
      TLOCAL=0
C ... CALL CHECKER THAT PUTS WATER ON/OFF ICE-FREE NODES
c ... (LAST ARG IS VALUE TO PUT ON ICE-FREE NODES, ZERO IT BEFORE
c     (turn on/off internally)
      CALL CHECKER(NMAX, NUMNP, NUMEL, NTYPE, KX, HTICE, DEPB, 
     &             WTHICK,0.0d0)
      CALL FORMT(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &           ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &           NZ, KZ, LM,
     &           BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &           AREAWET, AREADRY, AREATOT,
     &           DTLOCAL)
      IF(ISTART.EQ.0) THEN
        ISTART=1
        WSAVE=TOTALW
        IF(AREAWET.EQ.0.0) THEN
          TSAVE=0D0
        ELSE
          TSAVE=100.D0*TOTALW/AREAWET
        ENDIF
        IF(AREATOT.EQ.0.0) THEN
          PSAVE=0D0
        ELSE
          PSAVE=100.D0*AREAWET/AREATOT
        ENDIF
      ENDIF
C
      DTSAVE=0.0D0
      IC=0
1000  CONTINUE
      IC=IC+1
      CALL EFFECTIVE(NMAX,NZ,NUMNP,DTLOCAL,DTSAVE,A,B,D,
     &               WTHICK,ITKODE)
      LDIAG=.TRUE.
      CALL DIAGD(NMAX,NZ,NUMNP,A,KA,LDIAG,DMAX,SMAX)
      IF(SMAX.EQ.0) THEN
        RATIO=1.0D0
      ELSE
        RATIO=DMAX/SMAX
      ENDIF
C     PRINT *,'FIRST',DMAX,SMAX,RATIO
      IF(.NOT.LDIAG) THEN
        DTLOCAL=DTLOCAL*0.5D0
        NSTEP=INT(DT/DTLOCAL)
C       PRINT *,'GO AROUND AGAIN WITH SMALLER DT=',DTLOCAL,NSTEP
        IF(IC.LT.100) GOTO 1000
      ELSEIF(RATIO.GT.10) THEN
        DTLOCAL=DTLOCAL*2D0
        NSTEP=INT(DT/DTLOCAL)
        IF(DTLOCAL.LE.DT) THEN
C         PRINT *,'GO AROUND AGAIN WITH LARGER DT=',DTLOCAL,NSTEP
          IF(IC.LT.100) GOTO 1000
        ELSE
          DTLOCAL=DT
          NSTEP=1
        ENDIF
      ENDIF
      NSTEP=INT(DT/DTLOCAL)
      CALL EFFECTIVE(NMAX,NZ,NUMNP,DTLOCAL,DTSAVE,A,B,D,
     &                 WTHICK,ITKODE)
C
C      DO ISTEP=1,NSTEP
      ISTEP=0
      DOWHILE(TLOCAL.LT.DT)
        ISTEP=ISTEP+1
        TLOCAL=TLOCAL+DTLOCAL
        IF(NSTEP.GT.100) PRINT *,'DT,TLOCAL',REAL(DT),REAL(TLOCAL),
     &                          REAL(DTLOCAL),NSTEP
        DO JK=1,NUMNP
          TNEW(JK)=WTHICK(JK)
        ENDDO
C
C ..... CONJUGATE-GRADIENT ITERATIVE SOLVER
        ICOUNT=0
1100    CONTINUE
C        IF(LDIAG) THEN
C        IF(.TRUE.) THEN
C          IF(IOTOGG) THEN
C            WRITE(LIST(IPAGE+1),123) ' USING GAUSEID (FAST) '
C            IPAGE=IPAGE+1
C          ENDIF
C          WRITE(*,123) ' TIME BEFORE GAUSEID ',ETIME(TB),DTIME(TB)
          CALL GAUSEID(NMAX,NZ,NUMNP,A,KA,B,TNEW)
C          WRITE(*,123) ' TIME  AFTER GAUSEID ',ETIME(TB),DTIME(TB)
C        ELSE
C          WRITE(*,123) ' USING ASYMSL (SLOW) '
C          WRITE(*,123) ' TIME BEFORE ASYMSL ',ETIME(TB),DTIME(TB)
C          CALL ASYMSL(NMAX,NZ,NUMNP,A,KA,B,TNEW,1)
C          CALL ASYMSL(NMAX,NZ,NUMNP,A,KA,B,TNEW,2)
C          WRITE(*,123) ' TIME  AFTER ASYMSL ',ETIME(TB),DTIME(TB)
C        ENDIF
        IF(.FALSE.) THEN
          IFIX=0
          DO JK=1,NUMNP
            IF(DABS(TNEW(JK)).LT.1E-20) THEN
C             PRINT *,JK,TNEW(JK)
              TNEW(JK)=0D0
            ENDIF
            IF(TNEW(JK).LT.0D0) THEN
              IFIX=IFIX+1
              A(JK,1)=BIG
              B(JK)=0D0
            ENDIF
          ENDDO
C         PAUSE
          IF(IFIX.GT.0) THEN
            DO JK=1,NUMNP
              TNEW(JK)=WTHICK(JK)
            ENDDO
C           PRINT *,' GOING AROUND AGAIN !',IFIX
            ICOUNT=ICOUNT+1
            IF(ICOUNT.GT.10) THEN
              PRINT *,'PROBLEMS, 10 PASSES'
            ELSE
              GOTO 1100
            ENDIF
          ENDIF    
        ENDIF    
C
C ..... ACCEPT ONLY THICKNESSES GREATER THAN OR EQUAL TO ZERO
C ..... ACCEPT ONLY THICKNESSES    LESS THAN OR EQUAL TO 10.0
        DO JK=1,NUMNP
           WTHICK(JK)=TNEW(JK)
           WTHICK(JK)=MAX(0D0,TNEW(JK))
           WTHICK(JK)=MIN(10D0,WTHICK(JK))
        ENDDO
C ....  CALL EDGE DETECTOR THAT ALLOWS LEAKAGE OUT EDGES
c       (turn on/off internally)
        CALL EDGENODE(NMAX, NUMNP, NUMEL, NTYPE, KX, HTICE, DEPB, 
     &                WTHICK)
        IF(.FALSE. .AND. IPLOT.EQ.8) THEN
C          CALL NEWPAG
C          CALL CONTR(NMAX,NUMEL,X,Y,KX,
C     &                 HTICE,-2999.999D0,5500.D0,250.D0,
C     &                 -1000.D0,1000.D0,-1000.D0,1000.D0)
          CALL CONTR(NMAX,NUMEL,X,Y,KX,
     &             WTHICK,-0.0999999999D0,1.000000001D0,0.1D0,
     &            -1000.D0,1000.D0,-1000.D0,1000.D0)
        ENDIF

C        IF(ISTEP.NE.NSTEP) THEN
        IF(.TRUE.) THEN
          CALL FORMT(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &           ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &           NZ, KZ, LM,
     &           BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &           AREAWET, AREADRY, AREATOT,
     &           DTLOCAL)
          DTSAVE=0.0D0
          IC=0
1001      CONTINUE
          IC=IC+1
          CALL EFFECTIVE(NMAX,NZ,NUMNP,DTLOCAL,DTSAVE,A,B,D,
     &                     WTHICK,ITKODE)
          LDIAG=.TRUE.
          CALL DIAGD(NMAX,NZ,NUMNP,A,KA,LDIAG,DMAX,SMAX)
          IF(SMAX.EQ.0) THEN
            RATIO=1.0D0
          ELSE
            RATIO=DMAX/SMAX
          ENDIF
C         PRINT *,'SECOND',REAL(DMAX),REAL(SMAX),REAL(RATIO),LDIAG
          IF(.NOT.LDIAG) THEN
            DTLOCAL=DTLOCAL*0.5D0
            NSTEP=INT(DT/DTLOCAL)
C           PRINT *,'2ND:GO AROUND AGAIN WITH SMALLER DT=',DTLOCAL,NSTEP
            IF(IC.LT.100) GOTO 1001
          ELSEIF(RATIO.GT.10) THEN
            DTLOCAL=DTLOCAL*2D0
            NSTEP=INT(DT/DTLOCAL)
            IF(DTLOCAL.LE.(DT-TLOCAL)) THEN
C             PRINT *,'2ND:GO AROUND AGAIN WITH LARGER DT=',DTLOCAL,NSTEP
              IF(IC.LT.100) GOTO 1001
            ELSE
              DTLOCAL=DT-TLOCAL
              NSTEP=1
            ENDIF
          ENDIF
          NSTEP=INT(DT/DTLOCAL)
          CALL EFFECTIVE(NMAX,NZ,NUMNP,DTLOCAL,DTSAVE,A,B,D,
     &                     WTHICK,ITKODE)
        ENDIF
      ENDDO
      RATIOW=(TOTALW-WSAVE)/DT
      TOTALT=100.D0*TOTALW/AREAWET
      TOTALP=100.D0*AREAWET/AREATOT
      RATIOT=(TOTALT-TSAVE)/DT
      RATIOP=(TOTALP-PSAVE)/DT
      IF(IOTOGG) THEN
        WRITE(LIST(IPAGE+1),*) ISTEP,
     &       '********************************************'
        WRITE(LIST(IPAGE+2),*) 
     &         'TOTAL WATER (KM**3)=',REAL(1D-9*TOTALW),
     &         REAL(1D-9*RATIOW),REAL(1D-9*RATIOW*DT)
C       WRITE(LIST(IPAGE+1),*) 'WET AREA    (KM**2)=',REAL(1D-6*AREAWET)
C       WRITE(LIST(IPAGE+1),*) 'DRY AREA    (KM**2)=',REAL(1D-6*AREADRY)
C       WRITE(LIST(IPAGE+1),*) 'TOTAL AREA  (KM**2)=',REAL(1D-6*AREATOT)
        WRITE(LIST(IPAGE+3),*) 'PERCENT WET        =',REAL(TOTALP),
     &         REAL(RATIOP),REAL(RATIOP*DT)
        WRITE(LIST(IPAGE+4),*) 'AVG THICK   ( CM  )=',REAL(TOTALT),
     &         REAL(RATIOT),REAL(RATIOT*DT)
        WRITE(LIST(IPAGE+5),*) 
     &       '********************************************'
        IPAGE=IPAGE+5
      ENDIF
      WSAVE=TOTALW
      PSAVE=TOTALP
      TSAVE=TOTALT
C ... CALL CHECKER THAT PUTS WATER ON/OFF ICE-FREE NODES
c ... (LAST ARG IS VALUE TO PUT ON ICE-FREE NODES, 0.1 IT AFTER
c     (turn on/off internally)
      CALL CHECKER(NMAX, NUMNP, NUMEL, NTYPE, KX, HTICE, DEPB, 
     &             WTHICK,0.0d0)
      END
C-----------------------------------------------------------------
      SUBROUTINE FORMT(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &                 ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &                 NZ, KZ, LM,
     &                 BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &                 AREAWET, AREADRY, AREATOT,
     &                 DTLOCAL)
      IMPLICIT REAL*8(A-H,O-Z)
C FORM STIFFNESS MATRIX
      DIMENSION ALPHAC(3)
      DIMENSION BMELT(NMAX)
      DIMENSION ITKODE(NMAX)
      DIMENSION KZ(NMAX)
      DIMENSION LM(5),HTICE(NMAX),DEPB(NMAX)
      DIMENSION KX(NMAX,4),X(NMAX),Y(NMAX),NTYPE(NMAX),D(NMAX),B(NMAX)
      REAL*8 A(NMAX,NZ),WTHICK(NMAX)
      INTEGER KA(NMAX,NZ+1)
      DIMENSION P(5),S(5,5),DD(5)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2)
      DIMENSION PSI(4),DPSI(4,2)
      DIMENSION DWSIX(9),DWSIY(9)
      DIMENSION WSI(4),DWSI(4,2)
      DIMENSION XY(2,4),XI(2,9),ETA(2,9),W(2,9)
      DIMENSION ALPHA(2),BETA(2)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      LOGICAL UPWIND
      DATA BIG /1.D30/
      DATA ASCALE /1.D3/
C **********************************************************************
C     FORM CONDUCTIVITY MATRIX FOR COMPLETE BODY
C **********************************************************************
C
C ... ZERO OUT APPROPRIATE ARRAYS ...
      DO I=1,NUMNP
        KZ(I)=1
        KA(I,1)=I
        D(I)=0.0D0
        B(I)=0.0D0
        KA(I,NZ+1)=0
        DO J=1,NZ
          KA(I,J)=0
          A(I,J)=0.0D0
        ENDDO
      ENDDO
C
C ... BEGIN LOOP OVER ALL THE ELEMENTS ...
      TOTALW=0.0D0
      AREAWET=0.0D0
      AREADRY=0.0D0
      AREATOT=0.0D0
      DO 100 N=1,NUMEL
        IF(NTYPE(N).EQ.1) THEN
          NODEN=4
          NINT=9
        ELSE
          NODEN=3
          NINT=4
          PRINT *,'TRIANGLES NOT ALLOWED!'
          STOP
        ENDIF
        DO I=1,NODEN
          LM(I)=KX(N,I)
        ENDDO
C
C ..... FORM ELEMENT CONDUCTIVITY MATRIX
        DO I=1,NODEN
          DD(I)=0.0D0
          P(I)=0.0D0
          DO J=1,NODEN
            S(I,J)=0.0D0
          ENDDO
        ENDDO
C
        I=LM(1)
        J=LM(2)
        K=LM(3)
        L=LM(4)
        XY(1,1)=X(I)
        XY(1,2)=X(J)
        XY(1,3)=X(K)
        XY(2,1)=Y(I)
        XY(2,2)=Y(J)
        XY(2,3)=Y(K)
        IF(NODEN.EQ.4) THEN
          XY(1,4)=X(L)
          XY(2,4)=Y(L)
        ENDIF
C ..... USE FOLLOWING TO GENERATE CENTROID VALUES FOR GRADIENTS AND
C       OTHER MATERIAL PROPERTIES 
        IF(.TRUE.) THEN
          CALL FESHAPE(NTYPE(N),0D0,0D0,PSI,DPSI)
          CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
          TWTHIK=0.D0
          BMELTN=0.D0
          DGDX=0.D0
          DGDY=0.D0
          SURFX=0
          SURFY=0
          BEDX=0
          BEDY=0
          DO I=1,NODEN
            SURFX=SURFX+HTICE(LM(I))*DPSIX(I)
            SURFY=SURFY+HTICE(LM(I))*DPSIY(I)
          ENDDO
          SURFXY=SQRT(SURFX**2+SURFY**2)
          DO I=1,NODEN
            BEDX=BEDX+DEPB(LM(I))*DPSIX(I)
            BEDY=BEDY+DEPB(LM(I))*DPSIY(I)
            BMELTN=BMELTN+BMELT(LM(I))*PSI(I)
            TWTHIK=TWTHIK+WTHICK(LM(I))*PSI(I)
            IF(WTHICK(LM(I)).GT.0.001) THEN
C       ....  COEFFICIENT COMES FROM ALLEY'S PAPER
              RNPRES=5.0D-5*(HTICE(LM(I))-DEPB(LM(I)))*SURFXY/
     &              WTHICK(LM(I))
            ELSE
              RNPRES=0.D0
            ENDIF
            GGG=-(HTICE(LM(I))+0.09D0*DEPB(LM(I)))-RNPRES
            DGDX=DGDX+GGG*DPSIX(I)
            DGDY=DGDY+GGG*DPSIY(I)
          ENDDO
          DGXY=SQRT(DGDX**2+DGDY**2)
          SURFXY=SQRT(SURFX**2+SURFY**2)
          BEDXY=SQRT(BEDX**2+BEDY**2)
          IF(.FALSE.) THEN
            ANG=ATAN(DGDX/DGDY)*180.D0/3.14159D0
            ANGSURF=ATAN(SURFX/SURFY)*180.D0/3.14159D0
            ANGBED=ATAN(BEDX/BEDY)*180.D0/3.14159D0
            PRINT *,N,REAL(DGDX),REAL(DGDY),REAL(ANG),
     &                REAL(DGXY)
            PRINT *,N,REAL(SURFX),REAL(SURFY),REAL(ANGSURF),
     &                REAL(11*SURFXY)
            PRINT *,N,REAL(BEDX),REAL(BEDY),REAL(ANGBED),REAL(BEDXY)
          ENDIF
        ENDIF
C........... TO HERE
C TURN ON/OFF UPWINDING HERE (TRUE:ON)
      UPWIND=.FALSE.
      IF(UPWIND) THEN
C ..... VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
C SET ALPHA AND BETA FOR UPWINDING SHAPE FUNCTION
        ALPHA(1)=-1.D0
        ALPHA(2)=-1.D0
        BETA(1)=1.D0
        BETA(2)=1.D0
C        XCENT=(X(I)+X(J)+X(K)+X(L))/4.D0
C        YCENT=(Y(I)+Y(J)+Y(K)+Y(L))/4.D0
C ... THE FOLLOWING IS JUST TO GET THE DIRECTION OF THE GRADIENT...
        CALL FESHAPE(NTYPE(N),0D0,0D0,PSI,DPSI)
        CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
        DGDX=0.D0
        DGDY=0.D0
        SURFX=0
        SURFY=0
        BEDX=0
        BEDY=0
        DO I=1,NODEN
          GGG=-(HTICE(LM(I))+0.09D0*DEPB(LM(I)))
          DGDX=DGDX+GGG*DPSIX(I)
          DGDY=DGDY+GGG*DPSIY(I)
          SURFX=SURFX+HTICE(LM(I))*DPSIX(I)
          SURFY=SURFY+HTICE(LM(I))*DPSIY(I)
          BEDX=BEDX+DEPB(LM(I))*DPSIX(I)
          BEDY=BEDY+DEPB(LM(I))*DPSIY(I)
        ENDDO
        SURFXY=SQRT(SURFX**2+SURFY**2)
        BEDXY=SQRT(BEDX**2+BEDY**2)
C        IF(BEDXY.GE.11*SURFXY) THEN
C          PRINT *,NBEDXY,SURFXY
C          PAUSE
C        ENDIF
   
C ... THE PRECEEDING IS JUST TO GET THE DIRECTION OF THE GRADIENT...
C OPTIMIZE ALPHA AND BETA
        I=LM(1)
        J=LM(2)
        K=LM(3)
        L=LM(4)
        H12=SQRT((X(I)-X(J))**2+(Y(I)-Y(J))**2)
        H23=SQRT((X(J)-X(K))**2+(Y(J)-Y(K))**2)
        H43=SQRT((X(K)-X(L))**2+(Y(K)-Y(L))**2)
        H14=SQRT((X(L)-X(I))**2+(Y(L)-Y(I))**2)
C        PRINT *,'H',REAL(H12),REAL(H23),REAL(H43),REAL(H14)
C
        V12=-(DGDX*(X(J)-X(I))+DGDY*(Y(J)-Y(I)))/H12
        V23=-(DGDX*(X(K)-X(J))+DGDY*(Y(K)-Y(J)))/H23
        V43=(DGDX*(X(L)-X(K))+DGDY*(Y(L)-Y(K)))/H43
        V14=(DGDX*(X(I)-X(L))+DGDY*(Y(I)-Y(L)))/H14
C        PRINT *,'V',REAL(V12),REAL(V23),REAL(V43),REAL(V14)
C
        IF(ALPHAC(2).NE.0.0) THEN
          G12=1.0D0*ALPHAC(1)*V12*H12/ALPHAC(2)
          G23=1.0D0*ALPHAC(1)*V23*H23/ALPHAC(2)
          G43=1.0D0*ALPHAC(1)*V43*H43/ALPHAC(2)
          G14=1.0D0*ALPHAC(1)*V14*H14/ALPHAC(2)
        ELSE
          G12=0D0
          G23=0D0
          G43=0D0
          G14=0D0
        ENDIF
C        PRINT *,'G',REAL(G12),REAL(G23),REAL(G43),REAL(G14)
C
        IPR=0
        RMOD=0.1D0
        IF(G12.NE.0.) THEN
          ALPHA(1)=((1.D0/TANH(G12))-1.D0/G12)*RMOD
          IPR=1
        ELSE
          ALPHA(1)=0.D0
        ENDIF
        IF(G43.NE.0.) THEN
          ALPHA(2)=((1.D0/TANH(G43))-1.D0/G43)*RMOD
          IPR=1
        ELSE
          ALPHA(2)=0.D0
        ENDIF
        IF(G23.NE.0.) THEN
          BETA(1)=((1.D0/TANH(G23))-1.D0/G23)*RMOD
          IPR=1
        ELSE
          BETA(1)=0.D0
        ENDIF
        IF(G14.NE.0.) THEN
          BETA(2)=((1.D0/TANH(G14))-1.D0/G14)*RMOD
          IPR=1
        ELSE
          BETA(2)=0.D0
        ENDIF
C        PRINT *,'A',REAL(ALPHA(1)),REAL(ALPHA(2)),
C     &              REAL(BETA(1)),REAL(BETA(2))
        IF(IPR.EQ.2) THEN
          PRINT *,N,REAL(ALPHA(1)),REAL(ALPHA(2)),REAL(DGDX),REAL(DGDY)
          PRINT *,N,REAL(BETA(1)),REAL(BETA(2))
          IPR=0
        ENDIF
      ENDIF
C ..... ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C ..... FORM ELEMENT MATRIX AND VECTORS
C
C ..... BEGIN INTEGRATION POINT LOOP
        DO L=1,NINT
C SHAPE IS W/O UPWINDING, SHAPEUP IS WITH ...
          IF(UPWIND) THEN
            CALL SHAPEUP(ALPHA,BETA,
     &                   XI(NTYPE(N),L),ETA(NTYPE(N),L),
     &                   PSI,DPSI,WSI,DWSI)
C
C ........ GENERATE DEIRVATIVES, BOTH UPWINDING AND NON-UPWINDING
            CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
            CALL DERIVE(XY,DXDS,DWSI,DETJ,DSDX,DWSIX,DWSIY)
C
          ELSE
            CALL FESHAPE(NTYPE(N),XI(NTYPE(N),L),ETA(NTYPE(N),L),
     &                   PSI,DPSI)
C
C ........  GENERATE DEIRVATIVES, JUST NON-UPWINDING, COPY INTO UPWINDING
            CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
            DO I=1,NODEN
              WSI(I)=PSI(I)
              DWSI(I,1)=DPSI(I,1)
              DWSI(I,2)=DPSI(I,2)
              DWSIX(I)=DPSIX(I)
              DWSIY(I)=DPSIY(I)
            ENDDO
          ENDIF
C2345
C
C
C ....... ACCUMULATE INTEGRATION POINT VALUES OF INTEGRALS
          FAC=DETJ*W(NTYPE(N),L)
          TOTALW=TOTALW+TWTHIK*FAC
          IF(TWTHIK.GT.0.0) THEN
            AREAWET=AREAWET+FAC
          ELSE
            AREADRY=AREADRY+FAC
          ENDIF
          AREATOT=AREATOT+FAC
C ....... DELETE FOLLOWING TO USE ELEMENT CENTROID VALUES FOR GRADIENTS
C         AND OTHER MATERIAL PROPERTIES
          IF(.FALSE.) THEN
            DGDX=0.D0
            DGDY=0.D0
            TWTHIK=0.D0
            BMELTN=0.D0
            DO I=1,NODEN
              GGG=-(HTICE(LM(I))+0.09D0*DEPB(LM(I)))
              DGDX=DGDX+GGG*DPSIX(I)
              DGDY=DGDY+GGG*DPSIY(I)
              BMELTN=BMELTN+BMELT(LM(I))*PSI(I)
              TWTHIK=TWTHIK+WTHICK(LM(I))*PSI(I)
            ENDDO
C .....................................
          ENDIF
          IF(.FALSE.) THEN
C ......... OLD WAY ..................
            DO I=1,NODEN
C ........... THIS IS LUMPED CAPACITANCE MATRIX
C             DD(I)=DD(I)+PSI(I)*FAC
              DD(I)=DD(I)+WSI(I)*FAC
C ..........
C ..........  INTEGRATION POINT VALUES FROM FEM INTERPOLATIOON ....
C             P(I)=P(I)+BMELTN*PSI(I)*FAC
              P(I)=P(I)+BMELTN*WSI(I)*FAC
C
              DO J=1,NODEN
                TERM1=ALPHAC(2)*(DWSIX(I)*DPSIX(J)+DWSIY(I)*DPSIY(J))
C               TERM2=ALPHAC(1)*WSI(I)*(DGDX*DPSIX(J)+DGDY*DPSIY(J))
                TERM2=ALPHAC(1)*WSI(J)*(DGDX*DPSIX(I)+DGDY*DPSIY(I))
                TERM3=ALPHAC(3)*WSI(I)*PSI(J)*TWTHIK
                S(I,J)=S(I,J)+(TERM1+TERM2+TERM3)*FAC
              ENDDO
            ENDDO
C ......... END OLD WAY ..................
          ELSEIF(.FALSE.) THEN
C ......... FIRST NEW WAY ..................
            IPP=2
            IQQ=1
            GRADG=SQRT(DGDX**2+DGDY**2)
            IF(GRADG.EQ.0.D0) THEN
              ACONST=0.D0
              BCONST=0.D0
            ELSE
              ACONST=TWTHIK**IPP*GRADG**(IQQ-1)
              BCONST=TWTHIK**IPP*GRADG**(IQQ-1)
            ENDIF
C ......... NEED TO DEAL WITH POSSIBILITY THAT CONSTANT IS ZERO ...
            IF(ACONST.EQ.0) ACONST=1D-16
            IF(BCONST.EQ.0) BCONST=1D-16
C ...............................................................
            DO I=1,NODEN
C ........... THIS IS LUMPED CAPACITANCE MATRIX
              DD(I)=DD(I)+WSI(I)*FAC
              TERM1=ACONST*ALPHAC(2)*(DGDX*DWSIX(I)+DGDY*DWSIY(I))
              TERM2=BMELTN*WSI(I)
              P(I)=P(I)+(TERM1+TERM2)*FAC
              DO J=1,NODEN
                TERM2=BCONST*ALPHAC(1)*WSI(I)*
     &                (DGDX*DPSIX(J)+DGDY*DPSIY(J))
                TERM3=ALPHAC(3)*WSI(I)*PSI(J)*TWTHIK
                S(I,J)=S(I,J)+(TERM2+TERM3)*FAC
              ENDDO
            ENDDO
C ......... END FIRST NEW WAY ..................
          ELSE
C ......... LATEST NEW WAY ..................
            GRADG=SQRT(DGDX**2+DGDY**2)
            IPP=3
            IQQ=3
            IF(.FALSE. .AND. TWTHIK.GT.0.0) THEN
              IPP=2
              IQQ=1
              IF(GRADG.EQ.0.D0) THEN
                ACONST1=0.D0
                BCONST1=0.D0
                ACONST2=0.D0
                BCONST2=0.D0
              ELSE
                ACONST1=ASCALE*TWTHIK**IPP*GRADG**(IQQ-1)
                BCONST1=ASCALE*TWTHIK**IPP*GRADG**(IQQ-1)
                RPP=.5D0
                RQQ=.7D0
                ACONST2=ASCALE*TWTHIK**RPP*GRADG**(RQQ-1)/1D5
                BCONST2=ASCALE*TWTHIK**RPP*GRADG**(RQQ-1)/1D5
                PRINT *,'A',ACONST1,ACONST2,ACONST2/ACONST1
                PRINT *,'B',BCONST1,BCONST2,BCONST2/BCONST1
              ENDIF
            ENDIF
            IF(TWTHIK.LE.0.01) THEN
              IPP=2
              IQQ=1
              IF(GRADG.EQ.0.D0) THEN
                ACONST=0.D0
                BCONST=0.D0
              ELSE
                ACONST=ASCALE*TWTHIK**IPP*GRADG**(IQQ-1)
                BCONST=ASCALE*TWTHIK**IPP*GRADG**(IQQ-1)
              ENDIF
            ELSE
              RPP=.5D0
              RQQ=.7D0
              IF(GRADG.EQ.0.D0) THEN
                ACONST=0.D0
                BCONST=0.D0
              ELSE
                ACONST=ASCALE*TWTHIK**RPP*GRADG**(RQQ-1)/1D5
                BCONST=ASCALE*TWTHIK**RPP*GRADG**(RQQ-1)/1D5
              ENDIF
            ENDIF
C ......... NEED TO DEAL WITH POSSIBILITY THAT CONSTANT IS ZERO ...
C            IF(ACONST.EQ.0) THEN
C              ACONST=1E-16
C              BCONST=1E-16
C            ELSE
C              PRINT *,N,ACONST,TWTHIK
C            ENDIF
            DO I=1,NODEN
C ........... THIS IS LUMPED CAPACITANCE MATRIX
              DD(I)=DD(I)+WSI(I)*FAC
C ..........
C ..........  INTEGRATION POINT VALUES FROM FEM INTERPOLATIOON ....
              P(I)=P(I)+BMELTN*WSI(I)*FAC
C
              DO J=1,NODEN
                TERM1=-ACONST*ALPHAC(2)*
     &                (DGDX*(WSI(I)*DPSIX(J)+DWSIX(I)*PSI(J))+
     &                 DGDY*(WSI(I)*DPSIY(J)+DWSIY(I)*PSI(J)))
                TERM2=BCONST*ALPHAC(1)*WSI(I)*
     &                (DGDX*DPSIX(J)+DGDY*DPSIY(J))
                TERM3=ALPHAC(3)*WSI(I)*PSI(J)*TWTHIK
                S(I,J)=S(I,J)+(TERM1+TERM2+TERM3)*FAC
              ENDDO
            ENDDO
C ......... END LATEST NEW WAY ..................
          ENDIF
        ENDDO
C       IF(S(1,1).NE.0) THEN
C         PRINT *,N
C         DO II=1,4
C         PRINT *,(REAL(S(II,JJ)),JJ=1,4)
C         ENDDO
C       ENDIF
C
C ..... ADD ELEMENT CONDUCTIVITY TO COMPLETE CONDUCTIVITY MATRIX
C
        DO L=1,NODEN
          I=LM(L)
          D(I)=D(I)+DD(L)
C
C ....... THIS (D) IS LUMPED CAPACITANCE MATRIX
          B(I)=B(I)+P(L)
          DO M=1,NODEN
            J=LM(M)
            IF(I.EQ.J) THEN
              A(I,1)=A(I,1)+S(L,M)
              KA(I,1)=I
            ELSE
              DO K=2,KZ(I)
                IF(KA(I,K).EQ.J) THEN
                  A(I,K)=A(I,K)+S(L,M)
                  GOTO 99
                ENDIF
              ENDDO
              KZ(I)=KZ(I)+1
              A(I,KZ(I))=S(L,M)
              KA(I,KZ(I))=J
            ENDIF
99          CONTINUE
          ENDDO
          KA(I,NZ+1)=KZ(I)
        ENDDO
100   CONTINUE
C ... END LOOP OVER ALL THE ELEMENTS ...
C
C
C ... BOUNDARY CONDITIONS
C
C
C ... FIXED BOUNDARY CONDITIONS BY PENALTY METHOD
C
      NFIX=0
      DO N=1,NUMNP
        IF(ITKODE(N).EQ.1) THEN
          NFIX=NFIX+1
          A(N,1)=BIG
          B(N)=WTHICK(N)*BIG
        ENDIF
        KA(N,NZ+1)=KZ(N)
      ENDDO
C      IF(NFIX.EQ.0) THEN
C        PRINT *,'NO NODES ARE FIXED....'
C        PAUSE
C      ELSE
C        PRINT *,NFIX,' NODES ARE FIXED....',BIG
C      ENDIF

C ... REMOVE HERE, THIS WRITES OUT MATRIX
C      REWIND 73
C      WRITE(73,*) NUMNP
C      DO I=1,NUMNP
C        WRITE(73,*) KA(I,10)
C        DO J=1,KA(I,10)
C          WRITE(73,*) I,KA(I,J),A(I,J)
C        ENDDO
C      ENDDO
C      WRITE(73,*) (B(I),I=1,NUMNP)
C ... TO HERE ******
C      IF(IOTOGG) THEN
C        WRITE(LIST(IPAGE+1),*) 
C     &     'TOTAL WATER (KM**3)=',REAL(1E-9*TOTALW)
C        IPAGE=IPAGE+1
C       PRINT *,'WET AREA    (KM**2)=',REAL(1E-6*AREAWET)
C       PRINT *,'DRY AREA    (KM**2)=',REAL(1E-6*AREADRY)
C       PRINT *,'TOTAL AREA  (KM**2)=',REAL(1E-6*AREATOT)
C       PRINT *,'PERCENT WET        =',REAL(100.*AREAWET/AREATOT)
C       PRINT *,'AVG THICK   ( CM  )=',REAL(100.*TOTALW/AREAWET)
C      ENDIF
      END
C---------------------------------------------
      SUBROUTINE WRITEH2O(NUMNP,WTHICK,BMELT)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CHAR
      DIMENSION WTHICK(NUMNP),BMELT(NUMNP)
      PRINT *,'IN WRITEH2O',NUMNP
      PRINT *,'   TO WRITE OUT BACKUP OF WATER THICKNESS '
      PRINT *,'   INPUT Y'
      READ(*,1002) CHAR
      WRITE(99,1002) CHAR
1002  FORMAT(A1)
      IF(CHAR.EQ.'Y' .OR. CHAR.EQ.'y') THEN
        REWIND 42
        WRITE(42,*) NUMNP
        DO I=1,NUMNP
          WRITE(42,*) I,WTHICK(I),BMELT(I)
        ENDDO
      ENDIF
      END
C=================================================================
      SUBROUTINE SHAPEUP(AA,BB,XXI,ET,PSI,DPSI,WSI,DWSI)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AA(2),BB(2)
      DIMENSION PSI(4),DPSI(4,2)
      DIMENSION WSI(4),DWSI(4,2)
      A12=AA(1)
      A21=-A12
      A43=AA(2)
      A34=-A43
      B14=BB(2)
      B41=-B14
      B23=BB(1)
      B32=-B23
C CHANGE THESE 4 IF PROBLEMS OCCUR
C        A21=A12
C        A34=A43
C        B41=B14
C        B32=B23
      PSI(1)=.25D0*(1.D0-XXI)*(1.D0-ET)
      PSI(2)=.25D0*(1.D0+XXI)*(1.D0-ET)
      PSI(3)=.25D0*(1.D0+XXI)*(1.D0+ET)
      PSI(4)=.25D0*(1.D0-XXI)*(1.D0+ET)
      DPSI(1,2)=-.25D0*(1.D0-XXI)
      DPSI(2,2)=-.25D0*(1.D0+XXI)
      DPSI(3,2)=.25D0*(1.D0+XXI)
      DPSI(4,2)=.25D0*(1.D0-XXI)
      DPSI(1,1)=-.25D0*(1.D0-ET)
      DPSI(2,1)=.25D0*(1.D0-ET)
      DPSI(3,1)=.25D0*(1.D0+ET)
      DPSI(4,1)=-.25D0*(1.D0+ET)
      WSI(1)=PSI(1)*(1.D0+1.5D0*(A12*(1.D0+XXI)+B14*(1.D0+ET))+
     &        9.D0*A12*B14*PSI(3))
      WSI(2)=PSI(2)*(1.D0+1.5D0*(A21*(1.D0-XXI)+B23*(1.D0+ET))+
     &        9.D0*A21*B23*PSI(4))
      WSI(3)=PSI(3)*(1.D0+1.5D0*(A34*(1.D0-XXI)+B32*(1.D0-ET))+
     &       9.D0*A34*B32
     +       *PSI(1))
      WSI(4)=PSI(4)*(1.D0+1.5D0*(A43*(1.D0+XXI)+B41*(1.D0-ET))+
     &        9.D0*A43*B41
     +       *PSI(2))
      DWSI(1,1)=DPSI(1,1)*WSI(1)/PSI(1)+
     &PSI(1)*(1.5D0*B14+9.D0*A12*B14*DPSI(3,1))
      DWSI(2,1)=DPSI(2,1)*WSI(2)/PSI(2)+
     &PSI(2)*(1.5D0*B23+9.D0*A21*B23*DPSI(4,1))
      DWSI(3,1)=DPSI(3,1)*WSI(3)/PSI(3)+
     &PSI(3)*(-1.5D0*B32+9.D0*A34*B32*DPSI(1,1))
      DWSI(4,1)=DPSI(4,1)*WSI(4)/PSI(4)+
     &PSI(4)*(-1.5D0*B41+9.D0*A43*B41*DPSI(2,1))
      DWSI(1,2)=DPSI(1,2)*WSI(1)/PSI(1)+
     &PSI(1)*(1.5D0*A12+9.D0*A12*B14*DPSI(3,2))
      DWSI(2,2)=DPSI(2,2)*WSI(2)/PSI(2)+
     &PSI(2)*(-1.5D0*A21+9.D0*A21*B23*DPSI(4,2))
      DWSI(3,2)=DPSI(3,2)*WSI(3)/PSI(3)+
     &PSI(3)*(-1.5D0*A34+9.D0*A34*B32*DPSI(1,2))
      DWSI(4,2)=DPSI(4,2)*WSI(4)/PSI(4)+
     &PSI(4)*(1.5D0*A43+9.D0*A43*B41*DPSI(2,2))
C REMOVE FOLLOWING TO TURN ON UPWINDING
C     DO 100 I=1,4
C       WSI(I)=PSI(I)
C       DWSI(I,1)=DPSI(I,1)
C       DWSI(I,2)=DPSI(I,2)
100   CONTINUE
      RETURN
      END
C=================================================================
      SUBROUTINE DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DPSI(4,2)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2)
      DIMENSION XY(2,4)
C
C ....CALCULATE DXDS...EQUATION (5.3.6)
      DO I=1,2
        DO J=1,2
          DXDS(I,J)=0.0D0
          DO K=1,4
            DXDS(I,J)=DXDS(I,J)+DPSI(K,J)*XY(I,K)
          ENDDO
        ENDDO
      ENDDO
C
C .......   CALCULATE DSDX...EQUATION (5.2.7)
      DETJ=(DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1))
      IF(DETJ.LE.0.0) THEN
        WRITE(12,*) DETJ,((XY(MM,NN),NN=1,4),MM=1,2)
        WRITE(*,*) 'IN DERIVE, BAD JACOBIAN...'
        WRITE(*,*) DETJ,((XY(MM,NN),NN=1,4),MM=1,2)
        STOP
      ENDIF
      DENOM=1.D0/DETJ
      DSDX(1,1)=DXDS(2,2)*DENOM
      DSDX(2,2)=DXDS(1,1)*DENOM
      DSDX(1,2)=-DXDS(1,2)*DENOM
      DSDX(2,1)=-DXDS(2,1)*DENOM
C
C .......   CALCULATE D(PSI)/DX...EQUATION (5.3.5)
      DO I=1,4
        DPSIX(I)=DPSI(I,1)*DSDX(1,1)+DPSI(I,2)*DSDX(2,1)
        DPSIY(I)=DPSI(I,1)*DSDX(1,2)+DPSI(I,2)*DSDX(2,2)
      ENDDO
      END


C==========================================================
      SUBROUTINE EDGENODE(MXX, NUMNP, NUMEL, NTYPE, KX, HTICE, BED, 
     &                    WTHICK)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
C
C ... LOAD NODAL PROPERTIES INTO ELEMENT PROPERTIES
      DIMENSION NTYPE(MXX),
     &          HTICE(MXX),BED(MXX),
     &          LM(4),KX(MXX,4),WTHICK(NUMNP)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      DATA RHOI /0.917D0/, RHOW /1.092D0/
c......................................
      IF(.TRUE.) RETURN
c......................................
      NEDGE1=0
      NEDGE2=0
      NEDGE3=0
      NEDGE4=0
      AEDGE1=0.0D0
      AEDGE2=0.0D0
      AEDGE3=0.0D0
      AEDGE4=0.0D0
C ... NEW WAY ......
      DENOM=0.25D0
      DO N = 1,NUMEL
        AABED=0.D0
        AAHTICE=0.D0
        IEDGE1=0
        IEDGE2=0
        DO I=1,4
          LM(I)=KX(N,I)
          AABED=AABED+BED(LM(I))
          AAHTICE=AAHTICE+HTICE(LM(I))
          IF(BED(LM(I)).LT.SEALEV) THEN
            FLOT=(1.D0-RATDEN)*(BED(LM(I))-SEALEV)
            IF(HTICE(LM(I)).LT.FLOT) THEN
              IEDGE1=IEDGE1+1
            ENDIF
          ELSE
            IF(HTICE(LM(I))-.01D0.LE.BED(LM(I))) THEN
              IEDGE2=IEDGE2+1
            ENDIF
          ENDIF
        ENDDO
        IEDGE=IEDGE1+IEDGE2
        IF(IEDGE.EQ.4) THEN
          NEDGE4=NEDGE4+1
          RLEAK=1.0D0
          RKEEP=1.0D0-RLEAK
          DO I=1,4
            AEDGE4=AEDGE4+RLEAK*WTHICK(LM(I))
            WTHICK(LM(I))=RKEEP*WTHICK(LM(I))
          ENDDO
        ELSEIF(IEDGE.EQ.3) THEN
          NEDGE3=NEDGE3+1
          RLEAK=1.0D0
          RKEEP=1.0D0-RLEAK
          DO I=1,4
            AEDGE3=AEDGE3+RLEAK*WTHICK(LM(I))
            WTHICK(LM(I))=RKEEP*WTHICK(LM(I))
          ENDDO      
        ELSEIF(IEDGE.EQ.2) THEN
          NEDGE2=NEDGE2+1
          RLEAK=1.0D0
          RKEEP=1.0D0-RLEAK
          DO I=1,4
            AEDGE2=AEDGE2+RLEAK*WTHICK(LM(I))
            WTHICK(LM(I))=RKEEP*WTHICK(LM(I))
          ENDDO      
        ELSEIF(IEDGE.EQ.1) THEN
          NEDGE1=NEDGE1+1
          RLEAK=1.0D0
          RKEEP=1.0D0-RLEAK
          DO I=1,4
            AEDGE1=AEDGE1+RLEAK*WTHICK(LM(I))
            WTHICK(LM(I))=RKEEP*WTHICK(LM(I))
          ENDDO      
        ELSEIF(IEDGE.GT.0 .AND. IEDGE.LT.5) THEN
C ...... THIS STUFF BEYOND HERE NEVER HAPPENS....
          AABED=AABED*DENOM
          AAHTICE=AAHTICE*DENOM
          IF(AABED.LT.SEALEV) THEN
            FLOT=(1.D0-RATDEN)*(AABED-SEALEV)
            THIK=AAHTICE-FLOT
            IF(THIK.LT.0.) THEN
              IF(IEDGE.LT.5) THEN
                NEDGE=NEDGE+1
                DO I=1,4
                  AEDGE=AEDGE+WTHICK(LM(I))
                  WTHICK(LM(I))=0.0D0
                ENDDO
              ENDIF
            ENDIF
          ELSE
            THIK=AAHTICE-AABED
            IF(THIK.LT.0.) THEN
              IF(IEDGE.LT.5) THEN
                NEDGE=NEDGE+1
                DO I=1,4
                  AEDGE=AEDGE+WTHICK(LM(I))
                  WTHICK(LM(I))=0.0D0
                ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      NEDGE=NEDGE1+NEDGE2+NEDGE3+NEDGE4
      AEDGE=AEDGE1+AEDGE2+AEDGE3+AEDGE4
      IF(IOTOGG) THEN
        WRITE(LIST(IPAGE+1),*) NEDGE,'= NUMBER EDGE...',AEDGE
        IPAGE=IPAGE+1
      ENDIF
      END


      SUBROUTINE EFFECTIVE(NMAX,NZ,NUMNP,DTLOCAL,DTSAVE,A,B,D,
     &                     WTHICK,ITKODE)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(NMAX,NZ)
      DIMENSION D(NMAX),B(NMAX),WTHICK(NMAX),ITKODE(NMAX)
C ... IF SAVED TIME, THEN RESTORE MATRIX ...
      IF(DTSAVE.GT.0.0) THEN
        DT2=1.0D0/DTSAVE
        DO I=1,NUMNP
          IF (ITKODE(I).EQ.0) THEN
            IF (D(I).NE.0.) THEN
              A(I,1)=A(I,1)-DT2*D(I)
              B(I)=B(I)-DT2*D(I)*WTHICK(I)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C ... CALCULATE EFFECTIVE LOAD  AND STIFFNESS MATRIX
      IF(DTLOCAL.GT.0.0) THEN
        DT2=1.0D0/DTLOCAL
        DO I=1,NUMNP
          IF (ITKODE(I).EQ.0) THEN
            IF (D(I).NE.0.) THEN
              A(I,1)=A(I,1)+DT2*D(I)
C ... IS SIGN CORRECT HERE ?????????????????????????????
              B(I)=B(I)+DT2*D(I)*WTHICK(I)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      DTSAVE=DTLOCAL
      END


      SUBROUTINE WVELO(X, Y, KX, NUMNP, NUMEL,
     &                 WTHICK, HTICE, DEPB, LPRT)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM, NSMAX=MAXTIME)
      DIMENSION ALPHAC(3)
      DIMENSION LM(5),HTICE(NMAX),DEPB(NMAX)
      DIMENSION KX(NMAX,4),X(NMAX),Y(NMAX)
      DIMENSION WTHICK(NMAX)
      DIMENSION WVELX(NMAX),WVELY(NMAX)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2)
      DIMENSION PSI(4),DPSI(4,2)
      REAL*4 TB(2)
c     REAL*4 TB(2),ETIME,DTIME
c     EXTERNAL ETIME,DTIME
      DIMENSION XY(2,4),XARO(2),YARO(2)
      LOGICAL LPRT
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      COMMON /VELOS/ ASCAL,UMAX,VTHRESH,INORM
      LOGICAL IFLUSH
      COMMON /FLUSH/ IFLUSH
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      DIMENSION ICMAP(16)
      DATA ICMAP /0,4,12,6,13,2,8,7,9,3,10,5,11,1,14,15/   
      ASCALE=1D3
C      VTHRESH=1000.
      UMAXL=0.1
C      INORM=1
C      ASCAL=10
      UDELTA=UMAXL/13.D0
      NODEN=4
      UMAXL=0.0
      DO N=1,NUMEL
        DO I=1,NODEN
          LM(I)=KX(N,I)
          XY(1,I)=X(LM(I))
          XY(2,I)=Y(LM(I))
        ENDDO
C
        CALL FESHAPE(1,0D0,0D0,PSI,DPSI)
        CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
        TWTHIK=0.D0
        DGDX=0.D0
        DGDY=0.D0
        XCENT=0
        YCENT=0
        DO I=1,NODEN
          GGG=-(HTICE(LM(I))+0.09D0*DEPB(LM(I)))
          DGDX=DGDX+GGG*DPSIX(I)
          DGDY=DGDY+GGG*DPSIY(I)
          XCENT=XCENT+X(LM(I))*PSI(I)
          YCENT=YCENT+Y(LM(I))*PSI(I)
          TWTHIK=TWTHIK+WTHICK(LM(I))*PSI(I)
        ENDDO
        XCENT=XCENT/1000.D0
        YCENT=YCENT/1000.D0
        DGXY=SQRT(DGDX**2+DGDY**2)
        IPP=3
        IQQ=3
        IF(.FALSE. .AND. TWTHIK.GT.0.0) THEN
          IPP=2
          IQQ=1
          IF(GRADG.EQ.0.D0) THEN
            ACONST1=0.D0
            BCONST1=0.D0
            ACONST2=0.D0
            BCONST3=0.D0
          ELSE
            ACONST1=ASCALE*TWTHIK**IPP*DGXY**(IQQ-1)
            BCONST1=ASCALE*TWTHIK**IPP*DGXY**(IQQ-1)
            RPP=.5D0
            RQQ=.7D0
            ACONST2=ASCALE*TWTHIK**RPP*DGXY**(RQQ-1)/1D5
            BCONST2=ASCALE*TWTHIK**RPP*DGXY**(RQQ-1)/1D5
            PRINT *,'A',ACONST1,ACONST2,ACONST2/ACONST1
            PRINT *,'B',BCONST1,BCONST2,BCONST2/BCONST1
          ENDIF
        ENDIF
        IF(TWTHIK.LE.0.01) THEN
          IPP=2
          IQQ=1
          IF(GRADG.EQ.0.D0) THEN
            ACONST=0.D0
            BCONST=0.D0
          ELSE
            ACONST=ASCALE*TWTHIK**IPP*DGXY**(IQQ-1)
            BCONST=ASCALE*TWTHIK**IPP*DGXY**(IQQ-1)
          ENDIF
        ELSE
          RPP=.5D0
          RQQ=.7D0
          IF(GRADG.EQ.0.D0) THEN
            ACONST=0.D0
            BCONST=0.D0
          ELSE
            ACONST=ASCALE*TWTHIK**RPP*DGXY**(RQQ-1)/1D5
            BCONST=ASCALE*TWTHIK**RPP*DGXY**(RQQ-1)/1D5
          ENDIF
        ENDIF
        WVELX(N)=ACONST*DGDX
        WVELY(N)=ACONST*DGDY
        UMAG=SQRT(WVELX(N)**2+WVELY(N)**2)    
        UMAXL=MAX(UMAXL,UMAG)
      ENDDO
      IF(LPRT) PRINT *,'MAX VELOCITY:',UMAXL
      UDELTA=UMAX/13.D0
      DO N=1,NUMEL
        DO I=1,NODEN
          LM(I)=KX(N,I)
          XY(1,I)=X(LM(I))
          XY(2,I)=Y(LM(I))
        ENDDO
C
        CALL FESHAPE(1,0D0,0D0,PSI,DPSI)
        CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
        XCENT=0
        YCENT=0
        DO I=1,NODEN
          XCENT=XCENT+X(LM(I))*PSI(I)
          YCENT=YCENT+Y(LM(I))*PSI(I)
        ENDDO
        XCENT=XCENT/1000.D0
        YCENT=YCENT/1000.D0
        UX=WVELX(N)
        UY=WVELY(N)
        UMAG=SQRT(WVELX(N)**2+WVELY(N)**2)    
        IF(UMAG.LE.VTHRESH) THEN
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
          CALL MOVE(REAL(XARO(1)),REAL(YARO(1)))
          CALL DRAW(REAL(XARO(2)),REAL(YARO(2)))
C          CALL POINT(REAL(XARO(1)),REAL(YARO(1)))
C          IF(UMAG.NE.0) PRINT *,ICOLOR,UMAG,UMAX
        ENDIF
      ENDDO
      IF(IFLUSH) CALL GFLUSH
      END


C==========================================================
      SUBROUTINE CHECKER(MXX, NUMNP, NUMEL, NTYPE, KX, HTICE, BED, 
     &                    WTHICK,amount)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
C ... LOAD NODAL PROPERTIES INTO ELEMENT PROPERTIES
      DIMENSION NTYPE(MXX),
     &          HTICE(MXX),BED(MXX),
     &          LM(4),KX(MXX,4),WTHICK(NUMNP)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      DATA RHOI /0.917D0/, RHOW /1.092D0/
c...................................
      IF(.true.) RETURN
c...................................
      WTADD=0.D0
      NADD=0
      DO N = 1,NUMNP
        IF(BED(N).LT.SEALEV) THEN
          FLOT=(1.D0-RATDEN)*(BED(N)-SEALEV)
          IF(HTICE(N).LT.FLOT) THEN
            WTHICK(N)=amount
            WTADD=WTADD+WTHICK(N)
            NADD=NADD+1
          ENDIF
        ELSE
          IF(HTICE(N)-.01D0.LE.BED(N)) THEN
            WTHICK(N)=amount
            WTADD=WTADD+WTHICK(N)
            NADD=NADD+1
          ENDIF
        ENDIF
      ENDDO
c      IF(IOTOGG) THEN
c        WRITE(LIST(IPAGE+1),*) WTADD,
c     &                 ' WATER ADDED UNDER ICE-FREE NODES',NADD
c        IPAGE=IPAGE+1
c      ENDIF
      END
      SUBROUTINE WMOVERJ(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL, DT,
     &                 DTLOCAL,
     &                 ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &                 NZ, KZ, LM, TNEW,
     &                 BMELT, D, B, A, KA, ALPHAC, TOTALW, TOTALP,IPLOT,
     &                 WVELX,WVELY)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL LDIAG
      DIMENSION ALPHAC(3)
      DIMENSION BMELT(NMAX),TNEW(NMAX)
      DIMENSION ITKODE(NMAX)
      DIMENSION KZ(NMAX)
      DIMENSION LM(5),HTICE(NMAX),DEPB(NMAX)
      DIMENSION KX(NMAX,4),X(NMAX),Y(NMAX),NTYPE(NMAX),D(NMAX),B(NMAX)
      REAL*8 A(NMAX,NZ),WTHICK(NMAX),WVELX(NMAX),WVELY(NMAX)
      INTEGER KA(NMAX,NZ+1)
      DIMENSION XI(2,9),ETA(2,9),W(2,9)
      REAL*4 TB(2)
c     REAL*4 TB(2),ETIME,DTIME
c     EXTERNAL ETIME,DTIME
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      SAVE NSTEP,ISTART,WSAVE,TSAVE,PSAVE
      DATA ISTART /0/, BIG /1D30/
C      IF(DTLOCAL.EQ.DT) NSTEP=1
      DTLOCAL=DT
      TLOCAL=0
      CALL FORMTJ(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &           ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &           NZ, KZ, LM,
     &           BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &           AREAWET, AREADRY, AREATOT,
     &           DTLOCAL,WVELX,WVELY)
      IF(ISTART.EQ.0) THEN
        ISTART=1
        WSAVE=TOTALW
        IF(AREAWET.EQ.0.0) THEN
          TSAVE=0D0
        ELSE
          TSAVE=100.D0*TOTALW/AREAWET
        ENDIF
        IF(AREATOT.EQ.0.0) THEN
          PSAVE=0D0
        ELSE
          PSAVE=100.D0*AREAWET/AREATOT
        ENDIF
      ENDIF
C
      DTSAVE=0.0D0
      IC=0
1000  CONTINUE
      IC=IC+1
      CALL EFFECTIVE(NMAX,NZ,NUMNP,DTLOCAL,DTSAVE,A,B,D,
     &                     WTHICK,ITKODE)
      LDIAG=.TRUE.
      CALL DIAGD(NMAX,NZ,NUMNP,A,KA,LDIAG,DMAX,SMAX)
      IF(SMAX.EQ.0) THEN
        RATIO=1.0D0
      ELSE
        RATIO=DMAX/SMAX
      ENDIF
C     PRINT *,'FIRST',DMAX,SMAX,RATIO
      IF(.NOT.LDIAG) THEN
        DTLOCAL=DTLOCAL*0.5D0
        NSTEP=INT(DT/DTLOCAL)
C       PRINT *,'GO AROUND AGAIN WITH SMALLER DT=',DTLOCAL,NSTEP
        IF(IC.LT.100) GOTO 1000
      ELSEIF(RATIO.GT.10) THEN
        DTLOCAL=DTLOCAL*2D0
        NSTEP=INT(DT/DTLOCAL)
        IF(DTLOCAL.LE.DT) THEN
C         PRINT *,'GO AROUND AGAIN WITH LARGER DT=',DTLOCAL,NSTEP
          IF(IC.LT.100) GOTO 1000
        ELSE
          DTLOCAL=DT
          NSTEP=1
        ENDIF
      ENDIF
      NSTEP=INT(DT/DTLOCAL)
      CALL EFFECTIVE(NMAX,NZ,NUMNP,DTLOCAL,DTSAVE,A,B,D,
     &                     WTHICK,ITKODE)
C
C      DO ISTEP=1,NSTEP
      ISTEP=0
      DOWHILE(TLOCAL.LT.DT)
        ISTEP=ISTEP+1
        TLOCAL=TLOCAL+DTLOCAL
        IF(NSTEP.GT.100) PRINT *,'DT,TLOCAL',REAL(DT),REAL(TLOCAL),
     &                          REAL(DTLOCAL),NSTEP
        DO JK=1,NUMNP
          TNEW(JK)=WTHICK(JK)
        ENDDO
C
C ..... CONJUGATE-GRADIENT ITERATIVE SOLVER
        ICOUNT=0
1100    CONTINUE
123     FORMAT(A25,1PG13.6,G13.6)
C        IF(LDIAG) THEN
C        IF(.TRUE.) THEN
C          IF(IOTOGG) THEN
C            WRITE(LIST(IPAGE+1),123) ' USING GAUSEID (FAST) '
C            IPAGE=IPAGE+1
C          ENDIF
C          WRITE(*,123) ' TIME BEFORE GAUSEID ',ETIME(TB),DTIME(TB)
          CALL GAUSEID(NMAX,NZ,NUMNP,A,KA,B,TNEW)
C          WRITE(*,123) ' TIME  AFTER GAUSEID ',ETIME(TB),DTIME(TB)
C        ELSE
C          WRITE(*,123) ' USING ASYMSL (SLOW) '
C          WRITE(*,123) ' TIME BEFORE ASYMSL ',ETIME(TB),DTIME(TB)
C          CALL ASYMSL(NMAX,NZ,NUMNP,A,KA,B,TNEW,1)
C          CALL ASYMSL(NMAX,NZ,NUMNP,A,KA,B,TNEW,2)
C          WRITE(*,123) ' TIME  AFTER ASYMSL ',ETIME(TB),DTIME(TB)
C        ENDIF
        IF(.FALSE.) THEN
          IFIX=0
          DO JK=1,NUMNP
            IF(DABS(TNEW(JK)).LT.1E-20) THEN
C             PRINT *,JK,TNEW(JK)
              TNEW(JK)=0D0
            ENDIF
            IF(TNEW(JK).LT.0D0) THEN
              IFIX=IFIX+1
              A(JK,1)=BIG
              B(JK)=0D0
            ENDIF
          ENDDO
C         PAUSE
          IF(IFIX.GT.0) THEN
            DO JK=1,NUMNP
              TNEW(JK)=WTHICK(JK)
            ENDDO
C           PRINT *,' GOING AROUND AGAIN !',IFIX
            ICOUNT=ICOUNT+1
            IF(ICOUNT.GT.10) THEN
              PRINT *,'PROBLEMS, 10 PASSES'
            ELSE
              GOTO 1100
            ENDIF
          ENDIF    
        ENDIF    
C
C ..... ACCEPT ONLY THICKNESSES GREATER THAN OR EQUAL TO ZERO
C ..... ACCEPT ONLY THICKNESSES    LESS THAN OR EQUAL TO 10.0
        DO JK=1,NUMNP
           WTHICK(JK)=TNEW(JK)
           WTHICK(JK)=MAX(0D0,TNEW(JK))
           WTHICK(JK)=MIN(10D0,WTHICK(JK))
        ENDDO
C ....  CALL EDGE DETECTOR THAT ALLOWS LEAKAGE OUT EDGES
        CALL EDGENODE(NMAX, NUMNP, NUMEL, NTYPE, KX, HTICE, DEPB, 
     &                WTHICK)
        IF(.FALSE. .AND. IPLOT.EQ.8) THEN
C          CALL NEWPAG
C          CALL CONTR(NMAX,NUMEL,X,Y,KX,
C     &                 HTICE,-2999.999D0,5500.D0,250.D0,
C     &                 -1000.D0,1000.D0,-1000.D0,1000.D0)
          CALL CONTR(NMAX,NUMEL,X,Y,KX,
     &             WTHICK,-0.0999999999D0,1.000000001D0,0.1D0,
     &            -1000.D0,1000.D0,-1000.D0,1000.D0)
        ENDIF

C        IF(ISTEP.NE.NSTEP) THEN
        IF(.TRUE.) THEN
          CALL FORMTJ(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &           ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &           NZ, KZ, LM,
     &           BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &           AREAWET, AREADRY, AREATOT,
     &           DTLOCAL,WVELX,WVELY)
          DTSAVE=0.0D0
          IC=0
1001      CONTINUE
          IC=IC+1
          CALL EFFECTIVE(NMAX,NZ,NUMNP,DTLOCAL,DTSAVE,A,B,D,
     &                     WTHICK,ITKODE)
          LDIAG=.TRUE.
          CALL DIAGD(NMAX,NZ,NUMNP,A,KA,LDIAG,DMAX,SMAX)
          IF(SMAX.EQ.0) THEN
            RATIO=1.0D0
          ELSE
            RATIO=DMAX/SMAX
          ENDIF
C         PRINT *,'SECOND',REAL(DMAX),REAL(SMAX),REAL(RATIO),LDIAG
          IF(.NOT.LDIAG) THEN
            DTLOCAL=DTLOCAL*0.5D0
            NSTEP=INT(DT/DTLOCAL)
C           PRINT *,'2ND:GO AROUND AGAIN WITH SMALLER DT=',DTLOCAL,NSTEP
            IF(IC.LT.100) GOTO 1001
          ELSEIF(RATIO.GT.10) THEN
            DTLOCAL=DTLOCAL*2D0
            NSTEP=INT(DT/DTLOCAL)
            IF(DTLOCAL.LE.(DT-TLOCAL)) THEN
C             PRINT *,'2ND:GO AROUND AGAIN WITH LARGER DT=',DTLOCAL,NSTEP
              IF(IC.LT.100) GOTO 1001
            ELSE
              DTLOCAL=DT-TLOCAL
              NSTEP=1
            ENDIF
          ENDIF
          NSTEP=INT(DT/DTLOCAL)
          CALL EFFECTIVE(NMAX,NZ,NUMNP,DTLOCAL,DTSAVE,A,B,D,
     &                     WTHICK,ITKODE)
        ENDIF
      ENDDO
      RATIOW=(TOTALW-WSAVE)/DT
      TOTALT=100.D0*TOTALW/AREAWET
      TOTALP=100.D0*AREAWET/AREATOT
      RATIOT=(TOTALT-TSAVE)/DT
      RATIOP=(TOTALP-PSAVE)/DT
      IF(IOTOGG) THEN
        WRITE(LIST(IPAGE+1),*) ISTEP,
     &       '********************************************'
        WRITE(LIST(IPAGE+2),*) 
     &         'TOTAL WATER (KM**3)=',REAL(1D-9*TOTALW),
     &         REAL(1D-9*RATIOW),REAL(1D-9*RATIOW*DT)
C       WRITE(LIST(IPAGE+1),*) 'WET AREA    (KM**2)=',REAL(1D-6*AREAWET)
C       WRITE(LIST(IPAGE+1),*) 'DRY AREA    (KM**2)=',REAL(1D-6*AREADRY)
C       WRITE(LIST(IPAGE+1),*) 'TOTAL AREA  (KM**2)=',REAL(1D-6*AREATOT)
        WRITE(LIST(IPAGE+3),*) 'PERCENT WET        =',REAL(TOTALP),
     &         REAL(RATIOP),REAL(RATIOP*DT)
        WRITE(LIST(IPAGE+4),*) 'AVG THICK   ( CM  )=',REAL(TOTALT),
     &         REAL(RATIOT),REAL(RATIOT*DT)
        WRITE(LIST(IPAGE+5),*) 
     &       '********************************************'
        IPAGE=IPAGE+5
      ENDIF
      WSAVE=TOTALW
      PSAVE=TOTALP
      TSAVE=TOTALT
      END
C-----------------------------------------------------------------
      SUBROUTINE FORMTJ(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &                 ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &                 NZ, KZ, LM,
     &                 BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &                 AREAWET, AREADRY, AREATOT,
     &                 DTLOCAL,WVELX,WVELY)
      IMPLICIT REAL*8(A-H,O-Z)
C FORM STIFFNESS MATRIX
      DIMENSION ALPHAC(3)
      DIMENSION BMELT(NMAX)
      DIMENSION ITKODE(NMAX)
      DIMENSION KZ(NMAX)
      DIMENSION LM(5),HTICE(NMAX),DEPB(NMAX)
      DIMENSION KX(NMAX,4),X(NMAX),Y(NMAX),NTYPE(NMAX),D(NMAX),B(NMAX)
      REAL*8 A(NMAX,NZ),WTHICK(NMAX),WVELX(NMAX),WVELY(NMAX),NPRES
      INTEGER KA(NMAX,NZ+1)
      DIMENSION P(4),S(4,4),DD(4)
      DIMENSION P2(4)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2)
      DIMENSION PSI(4),DPSI(4,2)
      DIMENSION DWSIX(9),DWSIY(9)
      DIMENSION WSI(4),DWSI(4,2)
      DIMENSION XY(2,4),XI(2,9),ETA(2,9),W(2,9)
      DIMENSION ALPHA(2),BETA(2)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      LOGICAL UPWIND
      DATA BIG /1.D30/
      DATA ASCALE /1.D0/
C **********************************************************************
C     FORM CONDUCTIVITY MATRIX FOR COMPLETE BODY
C **********************************************************************
C
C ... ZERO OUT APPROPRIATE ARRAYS ...
      DO I=1,NUMNP
        KZ(I)=1
        KA(I,1)=I
        D(I)=0.0D0
        B(I)=0.0D0
        KA(I,NZ+1)=0
        DO J=1,NZ
          KA(I,J)=0
          A(I,J)=0.0D0
        ENDDO
      ENDDO
C
C ... BEGIN LOOP OVER ALL THE ELEMENTS ...
      TOTALW=0.0D0
      AREAWET=0.0D0
      AREADRY=0.0D0
      AREATOT=0.0D0
      DO 100 N=1,NUMEL
        IF(NTYPE(N).EQ.1) THEN
          NODEN=4
          NINT=9
        ELSE
          NODEN=3
          NINT=4
          PRINT *,'TRIANGLES NOT ALLOWED!'
          STOP
        ENDIF
        DO I=1,NODEN
          LM(I)=KX(N,I)
        ENDDO
C
C ..... FORM ELEMENT CONDUCTIVITY MATRIX
        DO I=1,NODEN
          DD(I)=0.0D0
          P(I)=0.0D0
          P2(I)=0.0D0
          DO J=1,NODEN
            S(I,J)=0.0D0
          ENDDO
        ENDDO
C
        I=LM(1)
        J=LM(2)
        K=LM(3)
        L=LM(4)
        XY(1,1)=X(I)
        XY(1,2)=X(J)
        XY(1,3)=X(K)
        XY(2,1)=Y(I)
        XY(2,2)=Y(J)
        XY(2,3)=Y(K)
        IF(NODEN.EQ.4) THEN
          XY(1,4)=X(L)
          XY(2,4)=Y(L)
        ENDIF
C ..... USE FOLLOWING TO GENERATE CENTROID VALUES FOR GRADIENTS AND
C       OTHER MATERIAL PROPERTIES 
        IF(.TRUE.) THEN
          CALL FESHAPE(NTYPE(N),0D0,0D0,PSI,DPSI)
          CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
          TWTHIK=0.D0
          BMELTN=0.D0
          DGDX=0.D0
          DGDY=0.D0
          SURFX=0
          SURFY=0
          BEDX=0
          BEDY=0
          HEADDX=0
          HEADDY=0
          DO I=1,NODEN
            SURFX=SURFX+HTICE(LM(I))*DPSIX(I)
            SURFY=SURFY+HTICE(LM(I))*DPSIY(I)
          ENDDO
          SURFXY=SQRT(SURFX**2+SURFY**2)
          DO I=1,NODEN
            IF(WTHICK(LM(I)).GT.000001) THEN
C             COEFFICIENT COMES FROM ALLEY'S PAPER
              NPRES=5.0D-5*(HTICE(LM(I))-DEPB(LM(I)))*SURFXY/
     &              WTHICK(LM(I))
            ELSE
              NPRES=0.D0
            ENDIF
C
C GGG IS THE POTENTIAL ENERGY, OR THE PRESSURE IN THE WATER SYSTEM. IT IS 
C SIMPLY THE ICE OVERBURDEN PRESSURE, RHO(ICE)*G*THICKNESS PLUS THE PUTENTIAL 
C ENERGY STORED IN THE WATER AT POSITION Z ABOVE SEA LEVEL. THICKNESS IS THE
C 
C INTRODUCE THE HEAD. THIS IS THE FLOATATION HEIGHT PLUS THE BED
C ELEVAION.
C COMPUTE IT'S DERIVATIVES.

C           GGG=-(HTICE(LM(I))+0.09D0*DEPB(LM(I)))-NPRES
            GGG=-(HTICE(LM(I))+0.092D0*DEPB(LM(I)))
C
            HEAD=-(0.82D0*HTICE(LM(I))+0.092D0*DEPB(LM(I)))
C    
            DGDX=DGDX+GGG*DPSIX(I)
            DGDY=DGDY+GGG*DPSIY(I)
C    
            HEADDX=HEADDX+HEAD*DPSIX(I)
            HEADDY=HEADDY+HEAD*DPSIY(I)
C    
            BEDX=BEDX+DEPB(LM(I))*DPSIX(I)
            BEDY=BEDY+DEPB(LM(I))*DPSIY(I)
            BMELTN=BMELTN+BMELT(LM(I))*PSI(I)
            TWTHIK=TWTHIK+WTHICK(LM(I))*PSI(I)
          ENDDO
          DGXY=SQRT(DGDX**2+DGDY**2)
          SURFXY=SQRT(SURFX**2+SURFY**2)
          BEDXY=SQRT(BEDX**2+BEDY**2)

        ENDIF
C........... TO HERE


C ..... ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C ..... FORM ELEMENT MATRIX AND VECTORS
C
C ..... BEGIN INTEGRATION POINT LOOP
        DO L=1,NINT
            CALL FESHAPE(NTYPE(N),XI(NTYPE(N),L),ETA(NTYPE(N),L),
     &                   PSI,DPSI)
C
C ........  GENERATE DEIRVATIVES, JUST NON-UPWINDING, COPY INTO UPWINDING
            CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
            DO I=1,NODEN
              WSI(I)=PSI(I)
              DWSI(I,1)=DPSI(I,1)
              DWSI(I,2)=DPSI(I,2)
              DWSIX(I)=DPSIX(I)
              DWSIY(I)=DPSIY(I)
            ENDDO

C2345
C
C
C ....... ACCUMULATE INTEGRATION POINT VALUES OF INTEGRALS
          FAC=DETJ*W(NTYPE(N),L)
          TOTALW=TOTALW+TWTHIK*FAC
          IF(TWTHIK.GT.0.0) THEN
            AREAWET=AREAWET+FAC
          ELSE
            AREADRY=AREADRY+FAC
          ENDIF
          AREATOT=AREATOT+FAC
C ....... DELETE FOLLOWING TO USE ELEMENT CENTROID VALUES FOR GRADIENTS
C         AND OTHER MATERIAL PROPERTIES
          IF(.FALSE.) THEN
            CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
            TWTHIK=0.D0
            BMELTN=0.D0
            DGDX=0.D0
            DGDY=0.D0
            SURFX=0
            SURFY=0
            BEDX=0
            BEDY=0
            HEADDX=0
            HEADDY=0
            DO I=1,NODEN
              SURFX=SURFX+HTICE(LM(I))*DPSIX(I)
              SURFY=SURFY+HTICE(LM(I))*DPSIY(I)
            ENDDO
            SURFXY=SQRT(SURFX**2+SURFY**2)
            DO I=1,NODEN
              IF(WTHICK(LM(I)).GT.000001) THEN
C               COEFFICIENT COMES FROM ALLEY'S PAPER
                NPRES=5.0D-5*(HTICE(LM(I))-DEPB(LM(I)))*SURFXY/
     &                WTHICK(LM(I))
              ELSE
                NPRES=0.D0
              ENDIF
C
C GGG IS THE POTENTIAL ENERGY, OR THE PRESSURE IN THE WATER SYSTEM. IT IS 
C SIMPLY THE ICE OVERBURDEN PRESSURE, RHO(ICE)*G*THICKNESS PLUS THE PUTENTIAL 
C ENERGY STORED IN THE WATER AT POSITION Z ABOVE SEA LEVEL. THICKNESS IS THE
C 
C INTRODUCE THE HEAD. THIS IS THE FLOATATION HEIGHT PLUS THE BED
C ELEVAION.
C COMPUTE IT'S DERIVATIVES.

C             GGG=-(HTICE(LM(I))+0.09D0*DEPB(LM(I)))-NPRES
              GGG=-(HTICE(LM(I))+0.092D0*DEPB(LM(I)))
C
              HEAD=-(0.82D0*HTICE(LM(I))+0.092D0*DEPB(LM(I)))
            
              DGDX=DGDX+GGG*DPSIX(I)
              DGDY=DGDY+GGG*DPSIY(I)
C
              HEADDX=HEADDX+HEAD*DPSIX(I)
              HEADDY=HEADDY+HEAD*DPSIY(I)
C
              BEDX=BEDX+DEPB(LM(I))*DPSIX(I)
              BEDY=BEDY+DEPB(LM(I))*DPSIY(I)
              BMELTN=BMELTN+BMELT(LM(I))*PSI(I)
              TWTHIK=TWTHIK+WTHICK(LM(I))*PSI(I)
            ENDDO
            DGXY=SQRT(DGDX**2+DGDY**2)
            SURFXY=SQRT(SURFX**2+SURFY**2)
            BEDXY=SQRT(BEDX**2+BEDY**2)

C .....................................
          ENDIF


C ......... LATEST NEW WAY ..................
C EXPONENTS FOR LAMINAR FLOW
          RPP=2.D0/3.D0
          RQQ=1.D0/2.D0
C EXPONENTS FOR TURBULENT FLOW
C FOLLOWING HOOKE'S SUGGESTION WATER THICKER THAN 1 CM IS TURBULENT
C         IF (TWTHIK.GE..01) THEN		
C           RPP=.5D0
C           RQQ=.7D0
C         ENDIF
            GRADG=SQRT(DGDX**2+DGDY**2)
            ACONST=ASCALE*TWTHIK**RPP*GRADG**(RQQ-1.D0)

C ......... NEED TO DEAL WITH POSSIBILITY THAT CONSTANT IS ZERO ...
            IF(ACONST.EQ.0) THEN
              ACONST=1D-16
            ELSE
C              PRINT *,N,ACONST,TWTHIK
            ENDIF
C.....GET WATER VELOCITIES
            IF (TWTHIK.GT.0) THEN
              DO I=1,NODEN
                WVELX(LM(I))=ACONST*DGDX
                WVELY(LM(I))=ACONST*DGDY
              ENDDO
            ELSE
              DO I=1,NODEN
                WVELX(LM(I))=0
                WVELY(LM(I))=0
              ENDDO
            ENDIF

            DO I=1,NODEN
C ........... THIS IS LUMPED CAPACITANCE MATRIX
              DD(I)=DD(I)+WSI(I)*FAC
C ..........
C ..........  INTEGRATION POINT VALUES FROM FEM INTERPOLATIOON ....
              RHS1=BMELTN*WSI(I)
              RHS2=ALPHAC(3)*6*(HEADDX*DPSIX(I)+HEADDY*DPSIY(I))
              P(I)=P(I)+(RHS1)*FAC
              P2(I)=P2(I)+(RHS2)*FAC
C
              DO J=1,NODEN
C......TERM 1 DOES NOT REALLY EXIST IN ANY ANALYSIS THAT I HAVE DONE, BUT
C MAY BE INCLUDED TO INDICATE DIFFUSION THROUGH THE AQUIFER?
C IT IS NICE BECASUE IT IMPARTS A NUMBERICAL STABILITY.
                TERM1=-ACONST*ALPHAC(2)*
     &                (DGDX*(WSI(I)*DPSIX(J)+DWSIX(I)*PSI(J))+
     &                 DGDY*(WSI(I)*DPSIY(J)+DWSIY(I)*PSI(J)))
C.... THE ORGIN AND NATURE OF TERM TO ARE DETAILED IN THE BASAL WATER PAPER
C THAT I HAVE WRITTEN. THIS IS REALLY THE RESULT OF MY RESEARCH SO FAR.
C...	ALPHAC(1) IS 1/(2**.666667*N) WHERE N IS MANNING ROUGHNESS COEFFICIENT
C...	REASONABLE VALUES OF ALPHAC(2) ARE BETWEEN 13 AND 16
                TERM2=ACONST*ALPHAC(1)*WSI(I)*
     &                (DGDX*DPSIX(J)+DGDY*DPSIY(J))
C....	THIS TERM IS LEAKAGE INTO THE AQUIFER, IT SHOULD DEPEND ON THE PRESSURE
C	IN THE SYSTEM.
C               TERM3=ALPHAC(3)*WSI(I)*PSI(J)*TWTHIK
C		TERM3=-ALPHAC(3)*6*(HEADDX*DPSIX(J)+HEADDY*DPSIY(J))	
		TERM3=0.D0
                S(I,J)=S(I,J)+(TERM1+TERM2+TERM3)*FAC
              ENDDO
            ENDDO
C ......... END LATEST NEW WAY ..................

        ENDDO
C      PRINT *,'P ',(REAL(P(II)),II=1,4)
C      PRINT *,'P2',(REAL(P2(II)),II=1,4)
C       IF(S(1,1).NE.0) THEN
C         PRINT *,N
C         DO II=1,4
C         PRINT *,(REAL(S(II,JJ)),JJ=1,4)
C         ENDDO
C       ENDIF
C
C ..... ADD ELEMENT CONDUCTIVITY TO COMPLETE CONDUCTIVITY MATRIX
C
        DO L=1,NODEN
          I=LM(L)
          D(I)=D(I)+DD(L)
C
C ....... THIS (D) IS LUMPED CAPACITANCE MATRIX
          B(I)=B(I)+P(L)+P2(L)
          DO M=1,NODEN
            J=LM(M)
            IF(I.EQ.J) THEN
              A(I,1)=A(I,1)+S(L,M)
              KA(I,1)=I
            ELSE
              DO K=2,KZ(I)
                IF(KA(I,K).EQ.J) THEN
                  A(I,K)=A(I,K)+S(L,M)
                  GOTO 99
                ENDIF
              ENDDO
              KZ(I)=KZ(I)+1
              A(I,KZ(I))=S(L,M)
              KA(I,KZ(I))=J
            ENDIF
99          CONTINUE
          ENDDO
          KA(I,NZ+1)=KZ(I)
        ENDDO
100   CONTINUE
C ... END LOOP OVER ALL THE ELEMENTS ...
C
C
C ... BOUNDARY CONDITIONS
C
C
C ... FIXED BOUNDARY CONDITIONS BY PENALTY METHOD
C
      NFIX=0
      DO N=1,NUMNP
        IF(ITKODE(N).EQ.1) THEN
          NFIX=NFIX+1
          A(N,1)=BIG
          B(N)=WTHICK(N)*BIG
        ENDIF
        KA(N,NZ+1)=KZ(N)
      ENDDO
C      IF(NFIX.EQ.0) THEN
C        PRINT *,'NO NODES ARE FIXED....'
C        PAUSE
C      ELSE
C        PRINT *,NFIX,' NODES ARE FIXED....',BIG
C      ENDIF

C ... REMOVE HERE, THIS WRITES OUT MATRIX
C     REWIND 73
C     WRITE(73,*) NUMNP
C     DO I=1,NUMNP
C       WRITE(73,*) KA(I,10)
C       DO J=1,KA(I,10)
C         WRITE(73,*) I,KA(I,J),A(I,J)
C       ENDDO
C     ENDDO
C     WRITE(73,*) (B(I),I=1,NUMNP)
C ... TO HERE ******
C      IF(IOTOGG) THEN
C        WRITE(LIST(IPAGE+1),*) 
C     &     'TOTAL WATER (KM**3)=',REAL(1E-9*TOTALW)
C        IPAGE=IPAGE+1
C       PRINT *,'WET AREA    (KM**2)=',REAL(1E-6*AREAWET)
C       PRINT *,'DRY AREA    (KM**2)=',REAL(1E-6*AREADRY)
C       PRINT *,'TOTAL AREA  (KM**2)=',REAL(1E-6*AREATOT)
C       PRINT *,'PERCENT WET        =',REAL(100.*AREAWET/AREATOT)
C       PRINT *,'AVG THICK   ( CM  )=',REAL(100.*TOTALW/AREAWET)
C      ENDIF
      END

      SUBROUTINE OPLATE(ITIME,NMAX,N3,NUMNP,NUMEL,X,Y,KX,THICK,
     &                 DELT,WWW,WRATE,WMIN,TIME,WWWORIG)
c-----------------------------------------------------------------------
c 4th order plate solver with one-time generation of stiffness and
c capacitance matrix. uses my sparse matrix storage for the static matrice
c and ITPACK sprse storage for the time-dependent modified matrices, and JCG
c iterative matrix solver. j fastook march 11, 1999
c-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)                                          
      include "parameter.h"
      PARAMETER(MMAX=MAXNUM,NZ=27, NZ1=NZ+1, M3=3*MMAX)
      parameter(nit=m3*nz,nmax1=m3+1,nmax3=m3*3)
      parameter(itmax=m3/10,ncg=4*itmax,nw=4*m3+ncg)
c ....arrays for ITPACK sparse storage...................
      dimension GK(nit),ja(nit),ia(nmax1)
      dimension iwksp(nmax3),wksp(nw),iwork(nit)
      dimension iparm(12),rparm(12)
      DIMENSION THICK(NMAX),WWW(N3),X(NMAX),Y(NMAX),KX(NMAX,4)
      DIMENSION WRATE(N3,2),WWWORIG(NMAX)
      DIMENSION KKX(MMAX,12),LTEMP(MMAX,3)
      DIMENSION XI(2,4),W(4)
      DIMENSION EK(12,12),EC(12,12),EF(12)
      DIMENSION GF(M3)
      DIMENSION GK0(M3,NZ),GC0(M3,NZ)
      DIMENSION KA(M3,NZ1),KZ(M3)
C     DIMENSION WWWSAVE(MMAX)
      CHARACTER*80 HED
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      common /oelastic/ wwwo(M3),worate(M3,2),wwwoorig(MMAX)
      common /elastic/ wwwe(M3),werate(M3,2),wwweorig(MMAX)
      common /viscous/ wwwv(M3),wvrate(M3,2),wwwvorig(MMAX)
      logical file40

C ... TIMER STUFF, FOR SGI ONLY ...!
      REAL*4 TB(2)
c     REAL*4 TB(2),ETIME,DTIME     !
c     EXTERNAL ETIME,DTIME         !
C .................................!
      SAVE IPASS,KKX,WSAVE,NNODE,XI,W,GK0,GC0,KA,KZ
C     SAVE WWWSAVE,iparm,rparm,IREAD
      DATA IPASS /0/,IREAD /0/
1000  FORMAT(1X,A,T25,1PG13.6,G13.6)
1001  FORMAT(1X,6(1X,1PG12.5))
c      WRITE(7,1000) ' TIME BEGINNING PLATE ',ETIME(TB),DTIME(TB)
      IF(IPASS.EQ.0) THEN
C ....  STUFF YOU ONLY DO THE FIRST PASS ..............
        IF(NMAX.NE.MMAX .OR.N3.NE.M3) THEN
          PRINT *,'PROBLEMS WITH NMAX,MMAX:',NMAX,MMAX
          STOP
        ENDIF
        CALL PCONNECT(NMAX,NUMNP,NUMEL,KX,KKX,LTEMP)
        WSAVE=0.D0
        CALL PSETINT(XI,W)
        CALL OSETMAT
        NNODE=3*NUMNP
        PRINT *,'READING BEDROCK DEPRESSION FILE'
        inquire(file='fort.40',exist=file40)
        iok=-1
        if(file40) then
          READ(40,*,IOSTAT=IOK) NNN
        endif
        IF(IOK.EQ.0) THEN
          IF(NNN.NE.NNODE) THEN
            PRINT *,'PROBLEMS:'
            PRINT *,'incompatible with current NNODE=',NNODE,NNN
            IOK=1
          ENDIF
          do i=1,nnode
            read(40,*) ii,wwwe(i),wwwv(i)
            werate(i,2)=wwwe(i)
            wvrate(i,2)=wwwv(i)
            if(i.ne.ii) then
              print *,'problems:reading wwwe,wwwv'
              print *,'incompatible with current (www)'
              iok=1
            endif
          enddo
          do i=1,nnode
            read(40,*) ii,werate(I,1),wvrate(I,1)
            if(i.ne.ii) then
              print *,'problems:reading wrates(1)'
              print *,'incompatible with current (www)'
              iok=1
            endif
          enddo
          do i=1,nnode
            read(40,*) ii,werate(I,2),wvrate(I,2)
            if(i.ne.ii) then
              print *,'problems:reading wrates(2)'
              print *,'incompatible with current (www)'
              iok=1
            endif
          enddo
          do i=1,numnp
            read(40,*) ii,wwweorig(i),wwwvorig(i)
            if(i.ne.ii) then
              print *,'problems:reading wwworigs'
              print *,'incompatible with current (wwworig)'
              iok=1
            endif
          enddo
          do i=1,numnp*3
            www(i)=wwwe(i)+wwwv(i)
c           if(mod(i,30).eq.1) print *,i,wwwe(i),wwwv(i)
            wrate(i,1)=werate(i,1)+wvrate(i,1)
            wrate(i,2)=werate(i,2)+wvrate(i,2)
          enddo
          do i=1,numnp
            wwworig(i)=wwweorig(i)+wwwvorig(i)
          enddo
          PRINT *,' BEDROCK DEPRESSION FILE FOUND'
          PRINT *,'    AND READ SUCCESSFULLY '
          if(itime.lt.0) then
            print *,'abandoning unloading'
            rewind 40
            return
          endif
          IREAD=1
        ENDIF
        IF(IOK.NE.0) THEN
          PRINT *,' NONE FOUND, SET TO ZERO ... '
c$doacross local(i)
          DO I=1,NNODE
            WWW(I)=0.D0
            WRATE(I,2)=0.D0
          ENDDO
          WRITE(HED,*) ' TIME=',NINT(TIME-DELT)
          WRITE(88) HED
          WRITE(88) (WWW(I),I=1,NNODE,3)
          WRITE(88) (THICK(I),I=1,NUMNP)
          WRITE(88) (1000*WRATE(I,1),I=1,NNODE,3)
        ENDIF
        CALL OFORMGK(NMAX,N3,NZ,NUMEL,NNODE,
     &              X,Y,KX,KKX,EK,EC,XI,W,
     &              GK0,GC0,KA,KZ)
      ENDIF
C ......................................................
C ....FORM STIFFNESS AND LOAD ..............
      CALL OFORMGF(NMAX,N3,NUMEL,NNODE,
     &            X,Y,THICK,KX,KKX,EF,XI,W,
     &            GF)
C ....TIME DEPENDENT CASE AND VARYING LOAD .......
      CALL OTIMEDEP(N3,NZ,NNODE,GK0,GC0,GK,GF,KA,KZ,
     &              WWW,DELT,nit,nmax1,ia,ja,iwork)
C ......................................................
C.....DUMP MATRIX FOR EXAMINATION ......................
      CALL ODUMPMAT(N3,NZ,NNODE,KZ,KA,GK0,GF,WWW,.false.)
      CALL NDUMPMAT(NIT,GK,JA,NMAX1,IA,NNODE,GF,.false.)
C.......................................................
C     IF(.TRUE.) THEN
C.......SOLVE EQUATIONS WITH JORDAN CONJUGATE-GRADIENT ITPACK ....!
c        WRITE(7,1000) ' TIME BEFORE JCG ',ETIME(TB),DTIME(TB)!
        call dfault(iparm,rparm)
        if(delt.ne.0.0) then
          iparm(1)=150  ! MAX number of iteration
        else
          iparm(1)=1000 ! MAX number of iteration
        endif
        rparm(1)=1d-6   ! ZETA, stopping criteria
c       iparm(2)=2      ! LEVEL of output (-1:none)
c       iparm(4)=7      ! OUTPUT unit number
        iparm(5)=1      ! NONsymmetric matrix (1)
        iparm(6)=0      ! NON-ADAPTIVE (0) (1 doesnt work)
c       iparm(10)=1     ! REMOVES large diagonal entries (doesnt work)
C        DO I=1,12
C          PRINT *,I,IPARM(I),RPARM(I)
C        ENDDO
        CALL JCG(NNODE,ia,ja,GK,GF,WWW,
     &         iwksp,nw,wksp,iparm,rparm,ier)
C         DO I=1,12
C           PRINT *,I,IPARM(I),RPARM(I)
C         ENDDO
        IF(IOTOGG) THEN
          write(list(ipage+1),*) ' relative error=',rparm(1),
     &            ' in iterations = ',iparm(1)
          write(list(ipage+2),*) rparm(11),rparm(12)
          ipage=ipage+2
        ENDIF
        if(ier.ne.0) then
          CALL NDUMPMAT(NIT,GK,JA,NMAX1,IA,NNODE,GF,.TRUE.)
          print *,'JCG error:',ier
c          print '(1x,i3,1pg13.6)',(i,rparm(i),i=1,12)
          pause
        endif
c        WRITE(7,1000) ' TIME AFTER JCG ',ETIME(TB),DTIME(TB) !
C ................................................................!
C     ELSE
C.......SOLVE EQUATIONS WITH CONJUGATE-GRADIENT ..................!
C       WRITE(7,1000) ' TIME BEFORE CONJUG ',ETIME(TB),DTIME(TB)  !
C       CALL CONJUG(N3,NZ,NNODE,1.D-6,GK,KA,GF,WWW)              !
C       WRITE(7,1000) ' TIME AFTER CONJUG ',ETIME(TB),DTIME(TB)   !
C ................................................................!
C     ENDIF
      IF(DELT.GT.0.) THEN
        DO I=1,NNODE
          WRATE(I,1)=(WWW(I)-WRATE(I,2))/DELT
        ENDDO
      ENDIF
c$doacross local(i)
      DO I=1,NNODE
        WRATE(I,2)=WWW(I)
      ENDDO
      WMIN=1D30
      WMAX=-1D30
      WTOT=0.D0
      DO I=1,NNODE,3
        WMIN=MIN(WMIN,WWW(I))
        WMAX=MAX(WMAX,WWW(I))
        WTOT=WTOT+WWW(I)
c        WRITE(*,*) (WWW(I+J),J=0,2)
        wwwo(i)=www(i)
        worate(i,1)=wrate(i,1)
        worate(i,2)=wrate(i,2)
      ENDDO
      THTOT=0.D0
      DO I=1,NUMNP
        THTOT=THTOT+THICK(I)
      ENDDO
      IF(IREAD.EQ.0 .AND. IPASS.EQ.0) THEN
        II=1
        DO I=1,NUMNP
          WWWORIG(I)=WWW(II)
          WWWOORIG(I)=WWW(II)
          II=II+3
        ENDDO
      ENDIF
      IF(IOTOGG) THEN
        write(list(ipage+1),*) 
     &       '*******************************************'
        WRITE(list(ipage+2),1001) TIME,REAL(THTOT),
     &              REAL(-THTOT/WTOT/ROCKICE),
     &              REAL(WMIN),REAL(WMAX),REAL(WSAVE-WMIN)
        write(list(ipage+3),*) 
     &       '*******************************************'
        ipage=ipage+3
      ENDIF
      WRITE(92,*) TIME,WMIN
      WSAVE=WMIN
      IPASS=1
      END
C=============================================================
      SUBROUTINE OELEMEK(XY,N,EKB,EC,NL,XI,W)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XY(2,N),EKB(12,12),EKS(12,12),EKH(12,12)
      DIMENSION EC(12,12)
      DIMENSION DPSIX(4),DPSIY(4)
      DIMENSION PSI(4),DPSI(4,2),XS(2)
      DIMENSION BB(3,12),BS(2,12),DB(3,3),DS(2,2)
      DIMENSION BTDB(3,12)
      DIMENSION XI(2,4),W(4)
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
C.....INITIALIZE ELEMENT ARRAYS
c$doacross local(i,j)
      DO I=1,12
        DO J=1,12
          EKB(I,J)=0.0D0
          EKS(I,J)=0.0D0
          EKH(I,J)=0.0D0
          EC(I,J)=0.0D0
        ENDDO
      ENDDO
C ... FORM D-MATRICES ...............
      TTT3=TTT**3/12.D0
      DB(1,1)=TTT3*(2.D0*RMU+RLAMBDA)
      DB(1,2)=TTT3*RLAMBDA
      DB(1,3)=0
      DB(2,1)=TTT3*RLAMBDA
      DB(2,2)=TTT3*(2.D0*RMU+RLAMBDA)
      DB(2,3)=0
      DB(3,1)=0
      DB(3,2)=0
      DB(3,3)=TTT3*RMU
      DS(1,1)=TTT*RMU
      DS(1,2)=0
      DS(2,1)=0
      DS(2,2)=TTT*RMU
C.....BEGIN 2X2 INTEGRATION LOOP
      DO L=1,NL
        XS(1)=XI(1,L)
        XS(2)=XI(2,L)
        CALL PGENSHAPE(N,XS,XY,PSI,DPSI,DETJ,DPSIX,DPSIY)
        CALL PLOADB(PSI,DPSIX,DPSIY,BB,BS)
C
C.......ACCUMULATE INTEGRATION POINT VALUE OF INTEGRALS
        FAC=DETJ*W(L)
C  .... RATIOS OF DENSITY OF ROCK AND ICE FOR OVERBURDEN LOAD
        XB=ROCKICE*RHOG
C ..... FORM CAPACITANCE MATRIX .....................
        DO I=1,4
          IP1=3*I-2
          IP2=3*I-1
          IP3=3*I
          DO J=1,4
            JP1=3*J-2
            JP2=3*J-1
            JP3=3*J
            EKH(IP1,JP1)=EKH(IP1,JP1)+FAC*XB*PSI(I)*PSI(J)
            TERM=FAC*CTIME*PSI(I)*PSI(J)
            EC(IP1,JP1)=EC(IP1,JP1)+TERM
            EC(IP1,JP2)=EC(IP1,JP2)+TERM
            EC(IP1,JP3)=EC(IP1,JP3)+TERM
            EC(IP2,JP1)=EC(IP2,JP1)+TERM
            EC(IP2,JP2)=EC(IP2,JP2)+TERM
            EC(IP2,JP3)=EC(IP2,JP3)+TERM
            EC(IP3,JP1)=EC(IP3,JP1)+TERM
            EC(IP3,JP2)=EC(IP3,JP2)+TERM
            EC(IP3,JP3)=EC(IP3,JP3)+TERM
          ENDDO
        ENDDO
C ..... FORM KB STIFFNESS MATRIX ..............
C ..... OK, NOW HERE GOES THE MATRIX MULTIPLICATION THAT 
C ..... GENERATES THE BT D B ALA BOOK...
C ..... FIRST BBT*DB (12X3)*(3X3)
c$doacross local(i,j,sum)
        DO I=1,12
          DO J=1,3
            SUM=0.D0
            DO K=1,3
              SUM=SUM+BB(K,I)*DB(J,K)
            ENDDO
            BTDB(J,I)=SUM
          ENDDO
        ENDDO
C THEN (BBT*DB)*BB (12X12)*(3X12)
c$doacross local(i,j,sum)
        DO I=1,12
          DO J=1,12
            SUM=0.D0
            DO K=1,3
              SUM=SUM+BTDB(K,I)*BB(K,J)
            ENDDO
            EKB(I,J)=EKB(I,J)+FAC*SUM
          ENDDO
        ENDDO
      ENDDO
C END OF 2X2 INTEGRATION LOOP
C
C.....BEGIN 1X1 REDUCED INTEGRATION LOOP
      DO L=1,1
        XS(1)=0.D0
        XS(2)=0.D0
        CALL PGENSHAPE(N,XS,XY,PSI,DPSI,DETJ,DPSIX,DPSIY)
        CALL PLOADB(PSI,DPSIX,DPSIY,BB,BS)
C
C.......ACCUMULATE INTEGRATION POINT VALUE OF INTEGRALS
        FAC=DETJ*4.D0
C ..... FORM KS STIFFNESS MATRIX ........................
C ..... OK, NOW HERE GOES THE MATRIX MULTIPLICATION THAT 
C ..... GENERATES BST*DS*BS ALA BOOK...
C ..... FIRST BST*DS (12X2)*(2X2)
c$doacross local(i,j,sum)
        DO I=1,12
          DO J=1,2
            SUM=0.D0
            DO K=1,2
              SUM=SUM+BS(K,I)*DS(J,K)
            ENDDO
            BTDB(J,I)=SUM
          ENDDO
        ENDDO
C THEN (BBT*DB)*BB (12X12)*(2X12)
c$doacross local(i,j,sum)
        DO I=1,12
          DO J=1,12
            SUM=0.D0
            DO K=1,2
              SUM=SUM+BTDB(K,I)*BS(K,J)
            ENDDO
            EKS(I,J)=EKS(I,J)+FAC*SUM
          ENDDO
        ENDDO
      ENDDO
C END OF 1X1 INTEGRATION LOOP
C
C ... COMBINE KB AND KS STIFFNESS MATRICES .......
      DO I=1,12
        DO J=1,12
C          PRINT '(2I3,1P3G13.6)',I,J,EKB(I,J),EKS(I,J),EKH(I,J)
          EKB(I,J)=EKB(I,J)+EKS(I,J)+EKH(I,J)
        ENDDO
      ENDDO
C      PAUSE
C      IF(.FALSE.) THEN
C        DO I=1,12
C          PRINT 1003,(EKB(I,J)/EKB(I,I),J=1,12)
C        ENDDO
C        PRINT *,'---------------------------'
C        DO I=1,12
C          PRINT 1003,(EC(I,J)/EC(I,I),J=1,12)
C        ENDDO
C        PAUSE
C      ENDIF
1003  FORMAT(1X,1P6G13.6)
      RETURN
1000  FORMAT(1X,'BAD JACOBIAN',E10.3)
1001  FORMAT(A,/,(2I5,3X,1PE13.6))
1002  FORMAT(A,/,(1P2E13.6))
      END
C=============================================================
      SUBROUTINE OELEMGF(ETHICK,XY,N,EF,NL,XI,W)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XY(2,N)
      DIMENSION EF(12)
      DIMENSION DPSIX(4),DPSIY(4)
      DIMENSION PSI(4),DPSI(4,2),XS(2)
      DIMENSION XI(2,4),W(4)
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
C.....INITIALIZE ELEMENT ARRAYS
c$doacross local(i)
      DO I=1,12
        EF(I)=0.0D0
      ENDDO
C.....BEGIN 2X2 INTEGRATION LOOP
      DO L=1,NL
        XS(1)=XI(1,L)
        XS(2)=XI(2,L)
        CALL PGENSHAPE(N,XS,XY,PSI,DPSI,DETJ,DPSIX,DPSIY)
C
C.......ACCUMULATE INTEGRATION POINT VALUE OF INTEGRALS
        FAC=DETJ*W(L)
        XF=-ETHICK
C ..... FORM CAPACITANCE MATRIX .....................
        DO I=1,4
          IP1=3*I-2
          EF(IP1)=EF(IP1)+XF*PSI(I)*FAC
        ENDDO
      ENDDO
      RETURN
1000  FORMAT(1X,'BAD JACOBIAN',E10.3)
1001  FORMAT(A,/,(2I5,3X,1PE13.6))
1002  FORMAT(A,/,(1P2E13.6))
      END
C=============================================================
      SUBROUTINE OPASSMBGK(N3,NZ,GK0,GC0,KA,KZ,EK,EC,N,NODE)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GK0(N3,NZ),GC0(N3,NZ)
      DIMENSION KA(N3,NZ+1),KZ(N3)
      DIMENSION EK(12,12),EC(12,12),NODE(N)
C........................
C||||||||||||||||||||||||
C     PRINT *,'IN PASSMB'
      DO L=1,N
        I=NODE(L)
        DO M=1,N
          J=NODE(M)
C.........ASSEMBLE GLOBAL STIFFNESS MATRIX GK
            IF(I.EQ.J) THEN
              GK0(I,1)=GK0(I,1)+EK(L,M)
              GC0(I,1)=GC0(I,1)+EC(L,M)
              KA(I,1)=I
            ELSE
              DO K=2,KZ(I)
                IF(KA(I,K).EQ.J) THEN
                  GK0(I,K)=GK0(I,K)+EK(L,M)
c FIX FIX FIX in other versions...
c                 GC0(I,K)=GC0(I,K)+EK(L,M)
                  GC0(I,K)=GC0(I,K)+EC(L,M)
c FIX FIX FIX in other versions...
                  GOTO 99
                ENDIF
              ENDDO
              KZ(I)=KZ(I)+1
              GK0(I,KZ(I))=EK(L,M)
              GC0(I,KZ(I))=EC(L,M)
              KA(I,KZ(I))=J
            ENDIF
99          CONTINUE
        ENDDO
      ENDDO
C||||||||||||||||||||||||
C........................
      END
C=============================================================
      SUBROUTINE OPASSMBGF(N3,GF,EF,N,NODE)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GF(N3)
      DIMENSION EF(12),NODE(N)
C........................
C||||||||||||||||||||||||
C     PRINT *,'IN PASSMBGF'
      DO L=1,N
        I=NODE(L)
C.......ASSEMBLE GLOBAL VECTOR GF
        GF(I)=GF(I)+EF(L)
      ENDDO
C||||||||||||||||||||||||
C........................
      END
C=============================================================
      SUBROUTINE OFORMGK(NMAX,N3,NZ,NUMEL,NNODE,
     &                  X,Y,KX,KKX,EK,EC,XI,W,
     &                  GK0,GC0,KA,KZ)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NMAX),Y(NMAX),KX(NMAX,4),KKX(NMAX,12)
      DIMENSION GK0(N3,NZ),GC0(N3,NZ)
      DIMENSION KA(N3,NZ+1),KZ(N3)
      DIMENSION XI(2,4),W(4)
      DIMENSION EK(12,12),EC(12,12)
      DIMENSION LM(4),XY(2,4),LLM(12)
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
c$doacross local(i,j)
      DO I=1,NNODE
        KZ(I)=1
        KA(I,1)=I
        KA(I,NZ+1)=0
        DO J=1,NZ
          KA(I,J)=0
          GK0(I,J)=0.0D0
          GC0(I,J)=0.0D0
        ENDDO
      ENDDO
      DO NEL=1,NUMEL
        DO L=1,4
          LM(L)=KX(NEL,L)
          XY(1,L)=X(LM(L))
          XY(2,L)=Y(LM(L))
        ENDDO
c$doacross local(L)
        DO L=1,12
          LLM(L)=KKX(NEL,L)
        ENDDO
        CALL OELEMEK(XY,4,EK,EC,4,XI,W)
C..................................
        IF(.false.) THEN
          write(7,*) NEL
          DO I=1,12
            write(7,*) (EK(I,J),J=1,12)
          ENDDO
c          PAUSE
        ENDIF
C..................................
        CALL OPASSMBGK(N3,NZ,GK0,GC0,KA,KZ,EK,EC,12,LLM)
C..................................
C        DO I=1,NNODE
C          PRINT 1000,(GK0(I,J),J=1,NNODE)
C        ENDDO
C        PAUSE
C..................................
      ENDDO
c$doacross local(i)
      DO I=1,NNODE
C        PRINT *,KZ(I),KA(I,NZ+1)
        KA(I,NZ+1)=KZ(I)
      ENDDO
C|||||||||||||||||||||||||||||||||||
C...................................
1000  FORMAT(1X,10F8.3)
      END
C=============================================================
      SUBROUTINE OFORMGF(NMAX,N3,NUMEL,NNODE,
     &                  X,Y,THICK,KX,KKX,EF,XI,W,
     &                  GF)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NMAX),Y(NMAX),KX(NMAX,4),KKX(NMAX,12)
      DIMENSION GF(N3)
      DIMENSION THICK(NMAX)
      DIMENSION XI(2,4),W(4)
      DIMENSION EF(12)
      DIMENSION LM(4),XY(2,4),LLM(12)
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
c$doacross local(i)
      DO I=1,NNODE
        GF(I)=0.D0
      ENDDO
      DO NEL=1,NUMEL
        ETHICK=0.D0
        DO L=1,4
          LM(L)=KX(NEL,L)
          XY(1,L)=X(LM(L))
          XY(2,L)=Y(LM(L))
          ETHICK=ETHICK+THICK(LM(L))
        ENDDO
c$doacross local(L)
        DO L=1,12
          LLM(L)=KKX(NEL,L)
        ENDDO
        ETHICK=RHOG*ETHICK/DBLE(4)
        CALL OELEMGF(ETHICK,XY,4,EF,4,XI,W)
C..................................
        IF(.false.) THEN
          write(7,*) NEL
          DO I=1,12
            write(7,*) EF(I)
          ENDDO
c          PAUSE
        ENDIF
C..................................
        CALL OPASSMBGF(N3,GF,EF,12,LLM)
C..................................
C        DO I=1,NNODE
C          PRINT 1000,GF(I)
C        ENDDO
C        PAUSE
C..................................
      ENDDO
C|||||||||||||||||||||||||||||||||||
C...................................
1000  FORMAT(1X,10F8.3)
      END
C=============================================================
      SUBROUTINE OTIMEDEP(N3,NZ,NNODE,GK0,GC0,GK,GF,KA,KZ,
     &                    WWW,DELT,nit,nmax1,ia,ja,iwork)
      IMPLICIT REAL*8(A-H,O-Z)
c ....arrays for ITPACK sparse storage...................
      dimension GK(nit),ja(nit),ia(nmax1)
      dimension iwork(nit)
      DIMENSION GF(N3),WWW(N3)
      DIMENSION GK0(N3,NZ),GC0(N3,NZ)
      DIMENSION KA(N3,NZ+1),KZ(N3)
      IF(DELT.GT.0.D0) THEN
        call sbini(NNODE,nit,ia,ja,GK,iwork)
c       call sbini(NNODE,NNODE,ia,ja,GK,iwork)
C ..... FORM MODIFIED STIFFNESS AND LOAD .......
        DELT1=1.D0/DELT    
        DO I = 1,NNODE
          DO J = 1,KZ(I)
            JG=KA(I,J)
C           PRINT *,I,J,GF(I),GC0(I,J)*WWW(JG)*DELT1
c            if(.false.) then ! full matrix...
              GF(I)=GF(I)+GC0(I,J)*WWW(JG)*DELT1
              GKIJ=GK0(I,J)+GC0(I,J)*DELT1
              call sbsij(NNODE,nit,ia,ja,GK,iwork,I,JG,GKIJ,
     &                   0,0,6,ier)
c            else ! lumped
c              GF(I)=GF(I)+GC0(I,J)*WWW(I)*DELT1
c              GKIJ=GK0(I,1)+GC0(I,J)*DELT1
c              call sbsij(NNODE,nit,ia,ja,GK,iwork,I,I,GKIJ,
c     &                   1,0,6,ier)
c            endif
          ENDDO 
        ENDDO
        call sbend(NNODE,nit,ia,ja,GK,iwork)
      else
        call sbini(NNODE,nit,ia,ja,GK,iwork)
c       call sbini(NNODE,NNODE,ia,ja,GK,iwork)
        DO I = 1,NNODE
          DO J = 1,KZ(I)
            JG=KA(I,J)
            GKIJ=GK0(I,J)
            call sbsij(NNODE,nit,ia,ja,GK,iwork,I,JG,GKIJ,0,1,6,ier)
          ENDDO 
        ENDDO
        call sbend(NNODE,nit,ia,ja,GK,iwork)
      ENDIF
      END
C=======================================
      SUBROUTINE OSETMAT
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
      logical ipass
      logical file39
      data ipass /.true./
      save ipass,e,rnu,thcrust,relax
      if(ipass) then
        ipass=.false.
        print *,'reading from input..defdep'
        inquire(file='fort.39',exist=file39)
        iok=-1
        if(file39) then
          read(39,*,IOSTAT=iok) e,rnu,thcrust,relax
        endif
        if(iok.ne.0) then
          print *,'not found, defaults used'
          e=1.D11
          rnu=0.5D0
          thcrust=110.D0
          relax=2000.D0                  
          rewind(39)
          write(39,*) e,rnu,thcrust,relax
          print *,'  youngs modulus: ', e
          print *,'  poissons ratio: ', rnu
          print *,' crust thickness: ', thcrust
          print *,' relaxation time: ', relax
        else
          print *,'defdep found, values used:'
          print *,'  youngs modulus: ', e
          print *,'  poissons ratio: ', rnu
          print *,' crust thickness: ', thcrust
          print *,' relaxation time: ', relax
        endif
      endif


c ... MATERIAL PROPERTIES AND THICKNESS OF CRUST
C ... YOUNGS MODULUS (NT/M**2)
c      E=5.D10
c      E=1.D11
C ... POISONS RATIO (DIMENSIONLESS)
c      RNU=0.5D0
C ... LAMES COEFFICIENTS
      RLAMBDA=RNU*E/(1.D0-RNU**2)
      RMU=0.5D0*E/(1.D0+RNU)
C      RMU=10000.D0
C      RLAMBDA=10000.D0
C ... THICKNESS OF CRUST (90-130 KM) (110,000 M)
c      TTT=110.D0*1000.D0
      TTT=thcrust*1000.D0
C ... RHOICE*GRAV (NT/M**3) TIMES THICKNESS, YIELDS PRESSURE, NT/M**2
      GRAV=9.8d0*0.3816d0
      RHOR=4000.D0
      RHOW=1092.D0
      RHOI=917.D0
      RHOG=RHOI*GRAV
      ROCKICE=RHOR/RHOI
C ... RELAXATION TIME CONSTANT ......!
c      RELAX=6000.D0                  !
c      RELAX=2000.D0                  !
      CTIME=3.16516D8*relax/6000.d0  !
      END
C---------------------------------------------
      SUBROUTINE WRITEDEPO(NUMNP,NNODE)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      parameter(nmax=maxnum,n3=3*nmax)
      CHARACTER*1 CHAR
      common /oelastic/ www(n3),wrate(n3,2),wwworig(nmax)
      PRINT *,'IN WRITEDEP',NNODE
      PRINT *,'   TO WRITE OUT BACKUP OF BEDROCK DEPRESSION '
      PRINT *,'   INPUT Y'
      READ(*,1002) CHAR
      WRITE(99,1002) CHAR
1002  FORMAT(A1)
      IF(CHAR.EQ.'Y' .OR. CHAR.EQ.'y') THEN
        REWIND 45
        WRITE(45,*) NNODE
        DO I=1,NNODE
          WRITE(45,*) I,WWW(I),0.d0
        ENDDO
        DO I=1,NNODE
          WRITE(45,*) I,wrate(I,1),0.d0
        ENDDO
        DO I=1,NNODE
          WRITE(45,*) I,wrate(I,2),0.d0
        ENDDO
        DO I=1,NUMNP
          WRITE(45,*) I,WWWORIG(I),0.d0
        ENDDO
      ENDIF
      END
      subroutine eplate(itime,numnp,numel,x,y,kx,thick,kode,
     &                 delt,www,wrate,wmin,time,wwworig,fnet,
     &                 wwwtot)
c-----------------------------------------------------------------------
c ... elastic solver .............................................
c 4th order plate solver with one-time generation of stiffness and
c capacitance matrix. uses my sparse matrix storage for the static matrice
c and itpack sprse storage for the time-dependent modified matrices, and jcg
c iterative matrix solver. j fastook july 2000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)                                          
      include "parameter.h"
      parameter(nmax=maxnum,nz=27, nz1=nz+1, n3=3*nmax)
      parameter(nit=n3*nz,nmax1=n3+1,nmax3=n3*3)
      parameter(itmax=n3/10,ncg=4*itmax,nw=4*n3+ncg)
c ....arrays for itpack sparse storage...................
      dimension gk(nit),ja(nit),ia(nmax1)
      dimension iwksp(nmax3),wksp(nw),iwork(nit)
      dimension iparm(12),rparm(12)
      dimension thick(nmax),x(nmax),y(nmax),kx(nmax,4)
      dimension wwwtot(n3),kode(nmax)
      dimension www(n3),wrate(n3,2),wwworig(nmax)
      dimension kkx(nmax,12),ltemp(nmax,3)
      dimension xi(2,4),w(4)
      dimension ek(12,12),ec(12,12),ef(12)
      dimension gf(n3)
      dimension gk0(n3,nz),gc0(n3,nz)
      dimension ka(n3,nz1),kz(n3)
c     dimension wwwsave(nmax)
      character*80 hed
      common /pmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      common /iolist/ list

c ... timer stuff, for sgi only ...!
      real*4 tb(2)
c     real*4 tb(2),etime,dtime     !
c     external etime,dtime         !
c .................................!
      save ipass,kkx,wsave,nnode,xi,w,gk0,gc0,ka,kz
c     save wwwsave,iparm,rparm,iread
      data ipass /0/,iread /0/
1000  format(1x,a,t25,1pg13.6,g13.6)
1001  format(6g12.5)
c      write(7,1000) ' time beginning plate ',etime(tb),dtime(tb)
      if(ipass.eq.0) then
c ....  stuff you only do the first pass ..............
c        if(nmax.ne.nmax .or.n3.ne.m3) then
c          print *,'problems with nmax,mmax:',nmax,mmax
c          stop
c        endif
        call pconnect(nmax,numnp,numel,kx,kkx,ltemp)
        wsave=0.d0
        call psetint(xi,w)
        call psetmat
        nnode=3*numnp
        call pformgk(nmax,n3,nz,numel,nnode,
     &              x,y,kx,kkx,ek,ec,xi,w,
     &              gk0,gc0,ka,kz)
      endif
c ......................................................
c ....form stiffness and load ..............
      call pformgf(nmax,n3,numel,nnode,
     &              x,y,thick,wwwtot,kx,kkx,ef,xi,w,
     &              gf,fnet)
c ....time dependent case and varying load ...........
c ....and load into sparse storage
      call ptimedep(n3,nz,nnode,gk0,gc0,gk,gf,ka,kz,
     &              www,delt,nit,nmax1,ia,ja,iwork,nmax,kode)
c ......................................................
c.....dump matrix for examination ......................
      call odumpmat(n3,nz,nnode,kz,ka,gk0,gf,www,.false.)
      call ndumpmat(nit,gk,ja,nmax1,ia,nnode,gf,.false.)
c.......................................................
c     if(.true.) then
c.......solve equations with jordan conjugate-gradient itpack ....!
c        write(7,1000) ' time before jcg ',etime(tb),dtime(tb)!
        call dfault(iparm,rparm)
        if(delt.ne.0.0) then
          iparm(1)=150  ! max number of iteration
        else
          iparm(1)=1000 ! max number of iteration
        endif
        rparm(1)=1d-6   ! zeta, stopping criteria
c       iparm(2)=2      ! level of output (-1:none)
c       iparm(4)=7      ! output unit number
        iparm(5)=1      ! nonsymmetric matrix (1)
        iparm(6)=0      ! non-adaptive (0) (1 doesnt work)
c       iparm(10)=1     ! removes large diagonal entries (doesnt work)
c        do i=1,12
c          print *,i,iparm(i),rparm(i)
c        enddo
        call jcg(nnode,ia,ja,gk,gf,www,
     &         iwksp,nw,wksp,iparm,rparm,ier)
c         do i=1,12
c           print *,i,iparm(i),rparm(i)
c         enddo
c        if(iotogg) then
c          write(list(ipage+1),*) ' relative error=',rparm(1),
c     &            ' in iterations = ',iparm(1)
c          write(list(ipage+2),*) rparm(11),rparm(12)
c          ipage=ipage+2
c        endif
        if(ier.ne.0) then
          call ndumpmat(nit,gk,ja,nmax1,ia,nnode,gf,.true.)
          print *,'jcg error:',ier
c          print '(1x,i3,1pg13.6)',(i,rparm(i),i=1,12)
          pause
        endif
c        write(7,1000) ' time after jcg ',etime(tb),dtime(tb) !
c ................................................................!
c     else
c.......solve equations with conjugate-gradient ..................!
c       write(7,1000) ' time before conjug ',etime(tb),dtime(tb)  !
c       call conjug(n3,nz,nnode,1.d-6,gk,ka,gf,www)              !
c       write(7,1000) ' time after conjug ',etime(tb),dtime(tb)   !
c ................................................................!
c     endif
      if(delt.gt.0.) then
        do i=1,nnode
          wrate(i,1)=(www(i)-wrate(i,2))/delt
        enddo
      endif
c$doacross local(i)
      do i=1,nnode
        wrate(i,2)=www(i)
      enddo
      wmin=1d30
      wmax=-1d30
      wtot=0.d0
      do i=1,nnode,3
        wmin=min(wmin,www(i))
        wmax=max(wmax,www(i))
        wtot=wtot+www(i)
c        write(*,*) (www(i+j),j=0,2)
      enddo
      thtot=0.d0
      do i=1,numnp
        thtot=thtot+thick(i)
      enddo
      if(iread.eq.0 .and. ipass.eq.0) then
        ii=1
        do i=1,numnp
          wwworig(i)=www(ii)
          ii=ii+3
        enddo
      endif
      if(iotogg) then
        write(list(ipage+1),*) 
     &       '*********** elastic **********************'
        write(list(ipage+2),1001) time,thtot,
     &              -thtot/wtot/rockice,
     &              wmin,wmax,-1000*(wsave-wmin)/delt
c        write(list(ipage+3),*) 
c     &       '*******************************************'
c        ipage=ipage+3
        ipage=ipage+2
      endif
      write(92,*) time,wmin
      wsave=wmin
      ipass=1
      end
c=============================================================
      subroutine psetint(xi,w)
      implicit real*8(a-h,o-z)
      dimension xi(2,4),w(4)
c ...............................
      xi(1,1)=1.d0/sqrt(3.d0)
      xi(2,1)=xi(1,1)
      xi(1,2)=-xi(1,1)
      xi(2,2)=xi(1,1)
      xi(1,3)=xi(1,1)
      xi(2,3)=-xi(1,1)
      xi(1,4)=-xi(1,1)
      xi(2,4)=-xi(1,1)
      w(1)=1.d0
      w(2)=1.d0
      w(3)=1.d0
      w(4)=1.d0
c ...............................
      end
c=============================================================
      subroutine pelemek(xy,n,ekb,ec,nl,xi,w)
      implicit real*8(a-h,o-z)
      dimension xy(2,n),ekb(12,12),eks(12,12),ekh(12,12)
      dimension ec(12,12)
      dimension dpsix(4),dpsiy(4)
      dimension psi(4),dpsi(4,2),xs(2)
      dimension bb(3,12),bs(2,12),db(3,3),ds(2,2)
      dimension btdb(3,12)
      dimension xi(2,4),w(4)
      common /pmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c.....initialize element arrays
c$doacross local(i,j)
      do i=1,12
        do j=1,12
          ekb(i,j)=0.0d0
          eks(i,j)=0.0d0
          ekh(i,j)=0.0d0
          ec(i,j)=0.0d0
        enddo
      enddo
c ... form d-matrices ...............
      ttt3=ttt**3/12.d0
      db(1,1)=ttt3*(2.d0*rmu+rlambda)
      db(1,2)=ttt3*rlambda
      db(1,3)=0
      db(2,1)=ttt3*rlambda
      db(2,2)=ttt3*(2.d0*rmu+rlambda)
      db(2,3)=0
      db(3,1)=0
      db(3,2)=0
      db(3,3)=ttt3*rmu
      ds(1,1)=ttt*rmu
      ds(1,2)=0
      ds(2,1)=0
      ds(2,2)=ttt*rmu
c.....begin 2x2 integration loop
      do l=1,nl
        xs(1)=xi(1,l)
        xs(2)=xi(2,l)
        call pgenshape(n,xs,xy,psi,dpsi,detj,dpsix,dpsiy)
        call ploadb(psi,dpsix,dpsiy,bb,bs)
c
c.......accumulate integration point value of integrals
        fac=detj*w(l)
c  .... ratios of density of rock and ice for overburden load
        xb=rockice*rhog
c ..... form capacitance matrix .....................
        do i=1,4
          ip1=3*i-2
          ip2=3*i-1
          ip3=3*i
          do j=1,4
            jp1=3*j-2
            jp2=3*j-1
            jp3=3*j
            ekh(ip1,jp1)=ekh(ip1,jp1)+fac*xb*psi(i)*psi(j)
            term=fac*ctime*psi(i)*psi(j)
            ec(ip1,jp1)=ec(ip1,jp1)+term
            ec(ip1,jp2)=ec(ip1,jp2)+term
            ec(ip1,jp3)=ec(ip1,jp3)+term
            ec(ip2,jp1)=ec(ip2,jp1)+term
            ec(ip2,jp2)=ec(ip2,jp2)+term
            ec(ip2,jp3)=ec(ip2,jp3)+term
            ec(ip3,jp1)=ec(ip3,jp1)+term
            ec(ip3,jp2)=ec(ip3,jp2)+term
            ec(ip3,jp3)=ec(ip3,jp3)+term
          enddo
        enddo
c ..... form kb stiffness matrix ..............
c ..... ok, now here goes the matrix multiplication that 
c ..... generates the bt d b ala book...
c ..... first bbt*db (12x3)*(3x3)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,3
            sum=0.d0
            do k=1,3
              sum=sum+bb(k,i)*db(j,k)
            enddo
            btdb(j,i)=sum
          enddo
        enddo
c then (bbt*db)*bb (12x12)*(3x12)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,12
            sum=0.d0
            do k=1,3
              sum=sum+btdb(k,i)*bb(k,j)
            enddo
            ekb(i,j)=ekb(i,j)+fac*sum
          enddo
        enddo
      enddo
c end of 2x2 integration loop
c
c.....begin 1x1 reduced integration loop
      do l=1,1
        xs(1)=0.d0
        xs(2)=0.d0
        call pgenshape(n,xs,xy,psi,dpsi,detj,dpsix,dpsiy)
        call ploadb(psi,dpsix,dpsiy,bb,bs)
c
c.......accumulate integration point value of integrals
        fac=detj*4.d0
c ..... form ks stiffness matrix ........................
c ..... ok, now here goes the matrix multiplication that 
c ..... generates bst*ds*bs ala book...
c ..... first bst*ds (12x2)*(2x2)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,2
            sum=0.d0
            do k=1,2
              sum=sum+bs(k,i)*ds(j,k)
            enddo
            btdb(j,i)=sum
          enddo
        enddo
c then (bbt*db)*bb (12x12)*(2x12)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,12
            sum=0.d0
            do k=1,2
              sum=sum+btdb(k,i)*bs(k,j)
            enddo
            eks(i,j)=eks(i,j)+fac*sum
          enddo
        enddo
      enddo
c end of 1x1 integration loop
c
c ... combine kb and ks stiffness matrices .......
      do i=1,12
        do j=1,12
c          print '(2i3,1p3g13.6)',i,j,ekb(i,j),eks(i,j),ekh(i,j)
c          ekb(i,j)=ekb(i,j)+eks(i,j)+ekh(i,j)
          ekb(i,j)=ekb(i,j)+eks(i,j)
        enddo
      enddo
c      pause
c      if(.false.) then
c        do i=1,12
c          print 1003,(ekb(i,j)/ekb(i,i),j=1,12)
c        enddo
c        print *,'---------------------------'
c        do i=1,12
c          print 1003,(ec(i,j)/ec(i,i),j=1,12)
c        enddo
c        pause
c      endif
1003  format(1x,1p6g13.6)
      return
1000  format(1x,'bad jacobian',e10.3)
1001  format(a,/,(2i5,3x,1pe13.6))
1002  format(a,/,(1p2e13.6))
      end
c=============================================================
      subroutine pelemgf(ethick,xy,n,ef,nl,xi,w)
      implicit real*8(a-h,o-z)
      dimension xy(2,n)
      dimension ef(12)
      dimension dpsix(4),dpsiy(4)
      dimension psi(4),dpsi(4,2),xs(2)
      dimension xi(2,4),w(4)
c      common /pmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c.....initialize element arrays
c$doacross local(i)
      do i=1,12
        ef(i)=0.0d0
      enddo
c.....begin 2x2 integration loop
      do l=1,nl
        xs(1)=xi(1,l)
        xs(2)=xi(2,l)
        call pgenshape(n,xs,xy,psi,dpsi,detj,dpsix,dpsiy)
c
c.......accumulate integration point value of integrals
        fac=detj*w(l)
        xf=-ethick
c ..... form capacitance matrix .....................
        do i=1,4
          ip1=3*i-2
          ef(ip1)=ef(ip1)+xf*psi(i)*fac
        enddo
      enddo
      return
1000  format(1x,'bad jacobian',e10.3)
1001  format(a,/,(2i5,3x,1pe13.6))
1002  format(a,/,(1p2e13.6))
      end
c=============================================================
      subroutine passmbgk(n3,nz,gk0,gc0,ka,kz,ek,ec,n,node)
      implicit real*8(a-h,o-z)
      dimension gk0(n3,nz),gc0(n3,nz)
      dimension ka(n3,nz+1),kz(n3)
      dimension ek(12,12),ec(12,12),node(n)
c........................
c||||||||||||||||||||||||
c     print *,'in passmb'
      do l=1,n
        i=node(l)
        do m=1,n
          j=node(m)
c.........assemble global stiffness matrix gk
            if(i.eq.j) then
              gk0(i,1)=gk0(i,1)+ek(l,m)
              gc0(i,1)=gc0(i,1)+ec(l,m)
              ka(i,1)=i
            else
              do k=2,kz(i)
                if(ka(i,k).eq.j) then
                  gk0(i,k)=gk0(i,k)+ek(l,m)
c fix fix fix in other versions...
c                 gc0(i,k)=gc0(i,k)+ek(l,m)
                  gc0(i,k)=gc0(i,k)+ec(l,m)
c fix fix fix in other versions...
                  goto 99
                endif
              enddo
              kz(i)=kz(i)+1
              gk0(i,kz(i))=ek(l,m)
              gc0(i,kz(i))=ec(l,m)
              ka(i,kz(i))=j
            endif
99          continue
        enddo
      enddo
c||||||||||||||||||||||||
c........................
      end
c=============================================================
      subroutine aplybc(n3,gf,nmax,kode,nnode)
      implicit real*8(a-h,o-z)
      dimension gf(n3),kode(nmax)
      data big /1d30/
      do i=1,nnode
        n=1+(i-1)/3
        if(kode(n).eq.1) gf(i)=0.d0
      enddo
      end
c=============================================================
      subroutine passmbgf(n3,gf,ef,n,node)
      implicit real*8(a-h,o-z)
      dimension gf(n3)
      dimension ef(12),node(n)
c........................
c||||||||||||||||||||||||
c     print *,'in passmbgf'
      do l=1,n
        i=node(l)
c.......assemble global vector gf
        gf(i)=gf(i)+ef(l)
      enddo
c||||||||||||||||||||||||
c........................
      end
c=============================================================
      subroutine pformgk(nmax,n3,nz,numel,nnode,
     &                  x,y,kx,kkx,ek,ec,xi,w,
     &                  gk0,gc0,ka,kz)
      implicit real*8(a-h,o-z)
      dimension x(nmax),y(nmax),kx(nmax,4),kkx(nmax,12)
      dimension gk0(n3,nz),gc0(n3,nz)
      dimension ka(n3,nz+1),kz(n3)
      dimension xi(2,4),w(4)
      dimension ek(12,12),ec(12,12)
      dimension lm(4),xy(2,4),llm(12)
c      common /pmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c$doacross local(i,j)
      do i=1,nnode
        kz(i)=1
        ka(i,1)=i
        ka(i,nz+1)=0
        do j=1,nz
          ka(i,j)=0
          gk0(i,j)=0.0d0
          gc0(i,j)=0.0d0
        enddo
      enddo
      do nel=1,numel
        do l=1,4
          lm(l)=kx(nel,l)
          xy(1,l)=x(lm(l))
          xy(2,l)=y(lm(l))
        enddo
c$doacross local(l)
        do l=1,12
          llm(l)=kkx(nel,l)
        enddo
        call pelemek(xy,4,ek,ec,4,xi,w)
c..................................
        if(.false.) then
          write(7,*) nel
          do i=1,12
            write(7,*) (ek(i,j),j=1,12)
          enddo
c          pause
        endif
c..................................
        call passmbgk(n3,nz,gk0,gc0,ka,kz,ek,ec,12,llm)
c..................................
c        do i=1,nnode
c          print 1000,(gk0(i,j),j=1,nnode)
c        enddo
c        pause
c..................................
      enddo
c$doacross local(i)
      do i=1,nnode
c        print *,kz(i),ka(i,nz+1)
        ka(i,nz+1)=kz(i)
      enddo
c|||||||||||||||||||||||||||||||||||
c...................................
1000  format(1x,10f8.3)
      end
c=============================================================
      subroutine pformgf(nmax,n3,numel,nnode,
     &              x,y,thick,www,kx,kkx,ef,xi,w,
     &              gf,fnet)
      implicit real*8(a-h,o-z)
      dimension x(nmax),y(nmax),kx(nmax,4),kkx(nmax,12)
      dimension gf(n3),www(n3)
      dimension thick(nmax)
      dimension xi(2,4),w(4)
      dimension ef(12)
      dimension lm(4),xy(2,4),llm(12)
      common /pmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c$doacross local(i)
      do i=1,nnode
        gf(i)=0.d0
      enddo
      fnet=0.d0
      do nel=1,numel
        ethick=0.d0
        ewww=0.d0
        do l=1,4
          lm(l)=kx(nel,l)
          xy(1,l)=x(lm(l))
          xy(2,l)=y(lm(l))
          ethick=ethick+thick(lm(l))
          ewww=ewww+www((lm(l)-1)*3+1)
        enddo
c$doacross local(l)
        do l=1,12
          llm(l)=kkx(nel,l)
        enddo
        ethick=ethick/dble(4)
        ewww=ewww/dble(4)
c         print '(i5,5g13.6)',nel,real(ethick),real(ewww),
c     &          real(rhog*ethick),real(rockice*rhog*ewww),
c     &          real(rhog*ethick+rockice*rhog*ewww)
        ethick=rhog*ethick
        ewww=rockice*rhog*ewww
c       print *,nel,real(ethick),real(ewww),real(ethick+ewww)
        ethick=ethick+ewww
c ...... EXPERIMENTAL ..........
c        ethick=max(0.d0,ethick)
c ...... EXPERIMENTAL ..........
        fnet=fnet+ethick
        call pelemgf(ethick,xy,4,ef,4,xi,w)
c..................................
        if(.false.) then
          write(7,*) nel
          do i=1,12
            write(7,*) ef(i)
          enddo
c          pause
        endif
c..................................
        call passmbgf(n3,gf,ef,12,llm)
c..................................
c        do i=1,nnode
c          print 1000,gf(i)
c        enddo
c        pause
c..................................
      enddo
c|||||||||||||||||||||||||||||||||||
c...................................
1000  format(1x,10f8.3)
      end
c=============================================================
      subroutine pshape(xs,psi,dpsi)
      implicit real*8(a-h,o-z)
      dimension xs(2),psi(4),dpsi(4,2)
c...................................
c|||||||||||||||||||||||||||||||||||
      eta=xs(1)
      rnu=xs(2)
      px0=1.d0-eta
      px1=1.d0+eta
      py0=1.d0-rnu
      py1=1.d0+rnu
      psi(1)=px0*py0/4.d0
      psi(2)=px1*py0/4.d0
      psi(3)=px1*py1/4.d0
      psi(4)=px0*py1/4.d0
      dpsi(1,1)=-py0/4.d0
      dpsi(2,1)=-dpsi(1,1)
      dpsi(3,1)=py1/4.d0
      dpsi(4,1)=-dpsi(3,1)
      dpsi(1,2)=-px0/4.d0
      dpsi(2,2)=-px1/4.d0
      dpsi(3,2)=-dpsi(2,2)
      dpsi(4,2)=-dpsi(1,2)
c|||||||||||||||||||||||||||||||||||
c...................................
      end
c=============================================================
      subroutine ploadb(psi,dpsix,dpsiy,bb,bs)
      implicit real*8(a-h,o-z)
      dimension psi(4),dpsix(4),dpsiy(4)
      dimension bb(3,12),bs(2,12)
      do ia=1,4
        ip1=3*ia-2
        ip2=ip1+1
        ip3=ip1+2
        bb(1,ip1)=0.d0
        bb(2,ip1)=0.d0
        bb(3,ip1)=0.d0
        bb(1,ip2)=dpsix(ia)
        bb(2,ip2)=0.d0
        bb(3,ip2)=dpsiy(ia)
        bb(1,ip3)=0.d0
        bb(2,ip3)=dpsiy(ia)
        bb(3,ip3)=dpsix(ia)
        bs(1,ip1)=dpsix(ia)
        bs(2,ip1)=dpsiy(ia)
        bs(1,ip2)=-psi(ia)
        bs(2,ip2)=0.d0
        bs(1,ip3)=0.d0
        bs(2,ip3)=-psi(ia)
      enddo
      end
c=============================================================
      subroutine pconnect(nmax,numnp,numel,kx,kkx,ltemp)
      implicit real*8(a-h,o-z)
      dimension kx(nmax,4),kkx(nmax,12),ltemp(nmax,3)
      ic=1
      do i=1,numnp
        ltemp(i,1)=ic
        ltemp(i,2)=ic+1
        ltemp(i,3)=ic+2
        ic=ic+3
c        print *,i,(ltemp(i,j),j=1,3)
      enddo
c$doacross local(i,j,ip1,ip2,ip3)
      do i=1,numel
        do j=1,4
          ip1=3*j-2
          ip2=ip1+1
          ip3=ip1+2
          kkx(i,ip1)=ltemp(kx(i,j),1)
          kkx(i,ip2)=ltemp(kx(i,j),2)
          kkx(i,ip3)=ltemp(kx(i,j),3)
        enddo
c        print 1000,i,(kx(i,j),j=1,4)
c        print 1000,i,(kkx(i,j),j=1,12)
c        pause
      enddo
1000  format(1x,13i5)
      end
c=============================================================
      subroutine poutstf(nnode,www,wrate,numnp,thick,numcol,numlev,
     &           time,wdiff)
      implicit real*8(a-h,o-z)
      include "parameter.h"
      parameter(nmax=maxnum,n3=3*nmax)
      dimension www(nnode),thick(numnp),tmp(n3),wrate(nnode)
      dimension wdiff(numnp)
      character*80 hed
c ... turn it off ...
c      return
c
c ... depression ...
      do i=1,nnode
        tmp(i)=www(i)
        if(www(i).gt.0) tmp(i)=tmp(i)*1
      enddo
      if(.false.) then
        write(81,1000)
        write(81,1001) numcol,numlev
        write(81,1002) 0,numcol*1000,0,numlev*1000
        write(81,1003) 'depression'
        write(81,1008) time
        write(81,1004)
        write(81,1005)
        write(81,1006) (tmp(i),i=1,nnode,3)
        write(81,1007)
      endif
c ... load ...
      if(.false.) then
        write(80,1000)
        write(80,1001) numcol,numlev
        write(80,1002) 0,numcol*1000,0,numlev*1000
        write(80,1003) 'load'
        write(80,1008) time
        write(80,1004)
        write(80,1005)
        write(80,1006) (thick(i),i=1,numnp)
        write(80,1007)
      endif
c ... rate ... mm/yr
c$doacross local(i)
      do i=1,nnode
        tmp(i)=wrate(i)*1000.d0
      enddo
      if(.false.) then
        write(82,1000)
        write(82,1001) numcol,numlev
        write(82,1002) 0,numcol*1000,0,numlev*1000
        write(82,1003) 'rate'
        write(82,1008) time
        write(82,1004)
        write(82,1005)
        write(82,1006) (tmp(i),i=1,nnode,3)
        write(82,1007)
      endif
c ... depression difference ...
      if(.false.) then
        write(83,1000)
        write(83,1001) numcol,numlev
        write(83,1002) 0,numcol*1000,0,numlev*1000
        write(83,1003) 'difference'
        write(83,1008) time
        write(83,1004)
        write(83,1005)
        write(83,1006) (wdiff(i),i=1,numnp)
        write(83,1007)
      endif
      write(hed,*) ' time=',nint(time)
      write(88) hed
      write(88) (www(i),i=1,nnode,3)
      write(88) (thick(i),i=1,numnp)
      write(88) (1000*wrate(i),i=1,nnode,3)
c ... 
1000  format(1x,'rank 2')
1001  format(1x,'dimensions',2i7)
1002  format(1x,'bounds',4i13)
1003  format(1x,'name ',a)
1004  format(1x,'scalar')
1005  format(1x,'data')
1006  format(1x,1p5g14.6)      
1007  format(1x,'end')
1008  format(1x,'time',g13.6)
1033  format(1x,'name load')
      end
c=============================================================
      subroutine ptimedep(n3,nz,nnode,gk0,gc0,gk,gf,ka,kz,
     &                    www,delt,nit,nmax1,ia,ja,iwork,nmax,kode)
      implicit real*8(a-h,o-z)
c ....arrays for itpack sparse storage...................
      dimension gk(nit),ja(nit),ia(nmax1)
      dimension iwork(nit)
      dimension gf(n3),www(n3),kode(nmax)
      dimension gk0(n3,nz),gc0(n3,nz)
      dimension ka(n3,nz+1),kz(n3)
      data big /1d30/
      if(delt.gt.0.d0) then
        call sbini(nnode,nit,ia,ja,gk,iwork)
c ..... form modified stiffness and load .......
        delt1=1.d0/delt    
        do i = 1,nnode
          n=1+(i-1)/3
          if(kode(n).eq.1) then
            call sbsij(nnode,nit,ia,ja,gk,iwork,i,i,big,
     &                   0,0,6,ier)
            gf(i)=0.d0
          else
            do j = 1,kz(i)
              jg=ka(i,j)
c             print *,i,j,gf(i),gc0(i,j)*www(jg)*delt1
c             if(.false.) then ! full matrix...
                gf(i)=gf(i)+gc0(i,j)*www(jg)*delt1
                gkij=gk0(i,j)+gc0(i,j)*delt1
                call sbsij(nnode,nit,ia,ja,gk,iwork,i,jg,gkij,
     &                     0,0,6,ier)
c              else ! lumped
c                gf(i)=gf(i)+gc0(i,j)*www(i)*delt1
c                gkij=gk0(i,1)+gc0(i,j)*delt1
c                call sbsij(nnode,nit,ia,ja,gk,iwork,i,i,gkij,
c     &                     1,0,6,ier)
c              endif
            enddo 
          endif
        enddo
        call sbend(nnode,nit,ia,ja,gk,iwork)
      else
        call sbini(nnode,nit,ia,ja,gk,iwork)
        do i = 1,nnode
c          n=1+(i-1)/3
c          if(kode(n).eq.1) then
c            call sbsij(nnode,nit,ia,ja,gk,iwork,i,i,big,
c     &                   0,0,6,ier)
c            gf(i)=0.d0
c          else
            do j = 1,kz(i)
              jg=ka(i,j)
              gkij=gk0(i,j)
              call sbsij(nnode,nit,ia,ja,gk,iwork,i,jg,gkij,0,1,6,ier)
            enddo 
c          endif
        enddo
        call sbend(nnode,nit,ia,ja,gk,iwork)
      endif
      end
c=============================================================
      subroutine pgauseid(nmax,nz,n,eps,aa,ka,b,x)
      implicit real*8(a-h,o-z)                                          
      dimension aa(nmax,nz),ka(nmax,nz+1),b(nmax),x(nmax)
      real*8 sum
c      data eps /1d-6/, itmax /100/, tol /0.d0/
      data itmax /100/, tol /0.d0/
c ... a diffrent first guess ...
c     if(.false.) then
c       do i=1,n
c         x(i)=b(i)/aa(i,1)
c       enddo
c     endif
c ... ..........................
      xmax=-1d30
      do i=1,n
        jmax=ka(i,nz+1)
        sum=b(i)
        do j=2,jmax
          sum=sum-aa(i,j)*x(ka(i,j))
        enddo
        xnew=sum/aa(i,1)
        xmax=max(xmax,abs(x(i)-xnew))
        x(i)=xnew
      enddo
      errorg=presid(nmax,nz,n,aa,ka,b,x)
c     print 1001,   0,errorg
      write(7,1001) 0,errorg
c ... rare but possible return
      if(errorg.eq.0.) return
      do iter=1,itmax
      xmax=-1d30
        do i=1,n
          jmax=ka(i,nz+1)
          sum=b(i)
          do j=2,jmax
            sum=sum-aa(i,j)*x(ka(i,j))
          enddo
          xnew=sum/aa(i,1)
          xmax=max(xmax,abs(x(i)-xnew))
          x(i)=xnew
        enddo
c       if(.true.) then
c         print 1001,iter,xmax
c         print 1000,(x(i),i=1,n)
c         write(7,1001) iter,xmax
c         write(7,1000) (x(i),i=1,n)
c       endif
        error=presid(nmax,nz,n,aa,ka,b,x)
        ratio=error/errorg
c       print 1001,   iter,error,ratio,xmax
        write(7,1001) iter,error,ratio,xmax
c ..... a normal return when converged ...
        if(abs(ratio).lt.eps) then
c         print 1002,   'a:converged',iter,error,ratio,xmax
          write(7,1002) 'a:converged',iter,error,ratio,xmax
          return
        endif
c ..... another normal return when converged ...
        if(xmax.lt.tol) then
c         print 1002,   'b:converged',iter,error,ratio,xmax
          write(7,1002) 'b:converged',iter,error,ratio,xmax
          return
        endif
      enddo
      print *,   'didnot converge in ',itmax,''
      write(7,*) 'didnot converge in ',itmax
c     pause
1000  format(1x,'pgauseid:',(1x,t10,5g14.6))
1001  format(1x,'pgauseid:',t21,i5,3g13.6)
1002  format(1x,'pgauseid:',a,i5,3g13.6)
      end
c=======================================
      function presid(nmax,nz,n,aa,ka,b,x)
      implicit real*8(a-h,o-z)                                          
      dimension aa(nmax,nz),ka(nmax,nz+1),b(nmax),x(nmax)
      real*8 sumsq,sum
      sumsq=0.d0
      do i=1,n
        jmax=ka(i,nz+1)
        sum=0.d0
        do j=1,jmax
          sum=sum+aa(i,j)*x(ka(i,j))
        enddo
        sum=sum-b(i)
        sumsq=sumsq+sum**2
      enddo
      presid=sumsq
      end
c=======================================
      subroutine pgenshape(n,xs,xy,psi,dpsi,detj,dpsix,dpsiy)
      implicit real*8(a-h,o-z)                                          
      dimension xy(2,n)
      dimension dpsix(4),dpsiy(4),dxds(2,2),dsdx(2,2)
      dimension psi(4),dpsi(4,2),xs(2)
        call pshape(xs,psi,dpsi)
c.......calculate dsds...equation(5.3.6)
        do i=1,2
          do j=1,2
            dxds(i,j)=0.0d0
            do k=1,4
              dxds(i,j)=dxds(i,j)+dpsi(k,j)*xy(i,k)
            enddo
          enddo
        enddo
c.......calculate dsdx...equation(5.2.7)
        detj=dxds(1,1)*dxds(2,2)-dxds(1,2)*dxds(2,1)
        if(detj.le.0.0) then
          write(*,1000) detj
          write(*,1001) 'dxds',((i,j,dxds(i,j),i=1,2),j=1,2)
          write(*,1001) 'dpsi',((i,j,dpsi(j,i),i=1,2),j=1,4)
          write(*,1002) 'xy',((xy(i,j),i=1,2),j=1,n)
          stop
        endif
        dsdx(1,1)=dxds(2,2)/detj    
        dsdx(2,2)=dxds(1,1)/detj    
        dsdx(1,2)=-dxds(1,2)/detj    
        dsdx(2,1)=-dxds(2,1)/detj  
c.......calculate d(psi)/dx...equation(5.3.5)
        do i=1,n
          dpsix(i)=dpsi(i,1)*dsdx(1,1)+dpsi(i,2)*dsdx(2,1)  
          dpsiy(i)=dpsi(i,1)*dsdx(1,2)+dpsi(i,2)*dsdx(2,2) 
        enddo
1000  format(1x,'bad jacobian',e10.3)
1001  format(a,/,(2i5,3x,1pe13.6))
1002  format(a,/,(1p2e13.6))
      end
c=======================================
      subroutine ndumpmat(nit,gk,ja,nmax1,ia,nnode,gf,disp)
      implicit real*8(a-h,o-z)                                          
      logical disp
      dimension gk(nit),ja(nit),ia(nmax1),gf(nnode)
      if(disp) then
        do l=1,nnode
          sum=0.0d0
          do j=ia(l)+1,ia(l+1)-1
            sum=sum+gk(j)
          enddo
          gdiag=gk(ia(l))
          write(13,30) 'row',l,' gf=',gf(l),gdiag,sum,gdiag/sum
          write(13,10) (gk(j),j=ia(l),ia(l+1)-1)
          write(13,10) (gk(j)/gdiag,j=ia(l),ia(l+1)-1)
          write(13,20) (ja(j),j=ia(l),ia(l+1)-1)
        enddo
        pause
      endif
10    format(5x,1p6g13.6)
20    format(1x,6i13)
30    format(1x,a,i6,a,1p4g13.6)
      end
c=======================================
      subroutine odumpmat(n3,nz,nnode,kz,ka,gk,gf,www,disp)
      implicit real*8(a-h,o-z)                                          
      logical disp
      dimension gk(n3,nz),gf(n3),www(n3)
      dimension ka(n3,nz+1),kz(n3)
        if(.false.) then
          do i=1,nnode
            write(7,*) 'diagonal='
            write(7,*) kz(i),i,gk(i,1)
            sum=0
            do j=2,kz(i)
              sum=sum+abs(gk(i,j))
              jj=ka(i,j)
              write(7,*) i,jj,gk(i,j),gk(i,j)/gk(i,1)
            enddo
            if(sum/gk(i,1).gt.1.1) then
              write(7,*) 'eq:',i,' not diag dom.',gk(i,1)/sum
c              pause
            endif
c            pause
            write(7,*) i,gf(i)
          enddo
c          pause
        endif
        if(disp) then
          write(13,*) nnode
          do i=1,nnode
            write(13,*) (gk(i,j),j=1,nz)
            write(13,*) (ka(i,j),j=1,nz+1)
            write(13,*) gf(i),www(i)
          enddo
        endif
      end
c
c=======================================================================
c
c
c
      subroutine vplate(itime,numnp,numel,x,y,kx,thick,kode,
     &                 delt,www,wrate,wmin,time,wwworig,fnet,
     &                 wwwtot)
c-----------------------------------------------------------------------
c 4th order viscous plate solver with one-time generation of stiffness
c matrix. uses my sparse matrix storage for the static matrice
c and itpack sprse storage for the time-dependent modified matrices, and c             ipage=ipage+2

c iterative matrix solver. j fastook july, 1999
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)                                          
      include "parameter.h"
      parameter(nmax=maxnum,nz=27, nz1=nz+1, n3=3*nmax)
      parameter(nit=n3*nz,nmax1=n3+1,nmax3=n3*3)
      parameter(itmax=n3/10,ncg=4*itmax,nw=4*n3+ncg)
c ....arrays for itpack sparse storage...................
      dimension gk(nit),ja(nit),ia(nmax1)
      dimension iwksp(nmax3),wksp(nw),iwork(nit)
      dimension iparm(12),rparm(12)
      dimension thick(nmax),www(n3),x(nmax),y(nmax),kx(nmax,4)
      dimension wwwtot(n3),kode(nmax)
      dimension wrate(n3,2),wwworig(nmax)
      dimension kkx(nmax,12),ltemp(nmax,3)
      dimension xi(2,4),w(4)
      dimension ek(12,12),ef(12)
      dimension gf(n3)
      dimension gk0(n3,nz)
      dimension ka(n3,nz1),kz(n3)
      dimension wwwtmp(n3)
      character*80 hed
      common /vmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      common /iolist/ list

c ... timer stuff, for sgi only ...!
      real*4 tb(2)
c     real*4 tb(2),etime,dtime     !
c     external etime,dtime         !
c .................................!
      save ipass,kkx,wsave,nnode,xi,w,gk0,ka,kz
      save gk,ia,ja
c     save wwwsave,iparm,rparm,iread
      data ipass /0/,iread /0/
1000  format(1x,a,t25,1pg13.6,g13.6)
1001  format(6g12.5)
c      write(7,1000) ' time beginning plate ',etime(tb),dtime(tb)
      if(ipass.eq.0) then
c ....  stuff you only do the first pass ..............
c        if(nmax.ne.mmax .or.n3.ne.m3) then
c          print *,'problems with nmax,mmax:',nmax,mmax
c          stop
c        endif
        call pconnect(nmax,numnp,numel,kx,kkx,ltemp)
        wsave=0.d0
        call psetint(xi,w)
        call vsetmat
        nnode=3*numnp
        call vformgk(nmax,n3,nz,numel,nnode,
     &              x,y,kx,kkx,ek,xi,w,
     &              gk0,ka,kz)
        call vloadgk(n3,nz,nnode,gk0,gk,ka,kz,
     &                    nit,nmax1,ia,ja,iwork,nmax,kode)
      endif
      do i=1,nnode
        wwwtmp(i)=www(i)
      enddo
c ... 2-step method ...
      do ii=1,2
c ......................................................
c ......form stiffness and load ..............
        if(.true.) then
c ......  use own solution ....
          call vformgf(nmax,n3,numel,nnode,
     &              x,y,thick,wwwtmp,kx,kkx,ef,xi,w,
     &              gf,fnet,kode)
        else
c ......  use other solution ....
          call vformgf(nmax,n3,numel,nnode,
     &              x,y,thick,wwwtot,kx,kkx,ef,xi,w,
     &              gf,fnet,kode)
        endif
c ......................................................
c.......dump matrix for examination ......................
        call odumpmat(n3,nz,nnode,kz,ka,gk0,gf,wwwtmp,.false.)
        call ndumpmat(nit,gk,ja,nmax1,ia,nnode,gf,.false.)
c.......................................................
c       if(.true.) then
c.........solve equations with jordan conjugate-gradient itpack ....!
c         write(7,1000) ' time before jcg ',etime(tb),dtime(tb)!
          call dfault(iparm,rparm)
          if(delt.ne.0.0) then
            iparm(1)=1000  ! max number of iteration
          else
            iparm(1)=1000 ! max number of iteration
          endif
          rparm(1)=1d-3   ! zeta, stopping criteria
c         iparm(2)=2      ! level of output (-1:none)
c         iparm(4)=7      ! output unit number
          iparm(5)=1      ! nonsymmetric matrix (1)
          iparm(6)=0      ! non-adaptive (0) (1 doesnt work)
c         iparm(10)=1     ! removes large diagonal entries (doesnt work)
c         do i=1,12
c           print *,i,iparm(i),rparm(i)
c         enddo
          do i=1,nnode
            wrate(i,ii)=0.d0
          enddo
          call jcg(nnode,ia,ja,gk,gf,wrate(1,ii),
     &            iwksp,nw,wksp,iparm,rparm,ier)

c ...................................................
c ....... experimental...............................
c ....... remove rigid body motion ..................
c ....... this could be done better by imposing bc ..
c ....... along the edge of the grid ................
          if(.false.) then
            wratebase=wrate(1,ii)                     !
            do i=1,nnode                              !
              wrate(i,ii)=wrate(i,ii)-wratebase       !
            enddo                                     !
          endif
c ....... experimental, remove rigid body motion ....
c ...................................................

c         do i=1,12
c           print *,i,iparm(i),rparm(i)
c         enddo
c           if(iotogg) then
c             write(list(ipage+1),*) ' relative error=',rparm(1),
c     &              ' in iterations = ',iparm(1)
c             write(list(ipage+2),*) rparm(11),rparm(12)
c             ipage=ipage+2
c             ipage=ipage+1
c           endif
          if(ier.ne.0) then
            call ndumpmat(nit,gk,ja,nmax1,ia,nnode,gf,.true.)
            print *,'jcg error:',ier
            print '(1x,i3,i10,1pg13.6)',(i,iparm(i),rparm(i),i=1,12)
            pause
          endif
c         write(7,1000) ' time after jcg ',etime(tb),dtime(tb) !
c ................................................................!
c       else
c.........solve equations with conjugate-gradient ..................!
c         write(7,1000) ' time before conjug ',etime(tb),dtime(tb)  !
c         call conjug(n3,nz,nnode,1.d-6,gk,ka,gf,www)              !
c         write(7,1000) ' time after conjug ',etime(tb),dtime(tb)   !
c ................................................................!
c       endif
        if(ii.eq.1) then
          do i=1,nnode
            wwwtmp(i)=www(i)+wrate(i,ii)*delt
          enddo
        endif
      enddo
      do i=1,nnode
        www(i)=www(i)+0.5d0*(wrate(i,1)+wrate(i,2))*delt
      enddo
      wmin=1d30
      wmax=-1d30
      wtot=0.d0
      do i=1,nnode,3
        wmin=min(wmin,www(i))
        wmax=max(wmax,www(i))
        wtot=wtot+www(i)
c        write(*,*) (real(www(i+j)),j=0,2)
c        write(*,*) (real(wrate(i+j,1)),j=0,2)
      enddo
      thtot=0.d0
      do i=1,numnp
        thtot=thtot+thick(i)
      enddo
      if(iread.eq.0 .and. ipass.eq.0) then
        ii=1
        do i=1,numnp
          wwworig(i)=www(ii)
          ii=ii+3
        enddo
      endif
      if(iotogg) then
        write(list(ipage+1),*) 
     &       '*********** viscous **********************'
        write(list(ipage+2),1001) time,thtot,
     &              -thtot/wtot/rockice,
     &              wmin,wmax,-1000*(wsave-wmin)/delt
        write(list(ipage+3),*) 
     &       '*******************************************'
        ipage=ipage+3
      endif
      write(92,*) time,wmin
      wsave=wmin
      ipass=1
      end
c=============================================================
      subroutine velemek(xy,n,ekb,nl,xi,w)
      implicit real*8(a-h,o-z)
      dimension xy(2,n),ekb(12,12),eks(12,12)
      dimension dpsix(4),dpsiy(4)
      dimension psi(4),dpsi(4,2),xs(2)
      dimension bb(3,12),bs(2,12),db(3,3),ds(2,2)
      dimension btdb(3,12)
      dimension xi(2,4),w(4)
      common /vmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c.....initialize element arrays
c$doacross local(i,j)
      do i=1,12
        do j=1,12
          ekb(i,j)=0.0d0
          eks(i,j)=0.0d0
        enddo
      enddo
c ... form d-matrices ...............
c ... what is this ctime1? why is it 10 times ctime?
      ctime1=ctime*10
      ctime1=1.0d0
c ..................................................
      ttt3=ttt**3/12.d0
      db(1,1)=ttt3*(2.d0*rmu+rlambda)*ctime1
      db(1,2)=ttt3*rlambda*ctime1
      db(1,3)=0
      db(2,1)=ttt3*rlambda*ctime1
      db(2,2)=ttt3*(2.d0*rmu+rlambda)*ctime1
      db(2,3)=0
      db(3,1)=0
      db(3,2)=0
      db(3,3)=ttt3*rmu*ctime1
      ds(1,1)=ttt*rmu*ctime1
      ds(1,2)=0
      ds(2,1)=0
      ds(2,2)=ttt*rmu*ctime1
c.....begin 2x2 integration loop
      do l=1,nl
        xs(1)=xi(1,l)
        xs(2)=xi(2,l)
        call pgenshape(n,xs,xy,psi,dpsi,detj,dpsix,dpsiy)
        call ploadb(psi,dpsix,dpsiy,bb,bs)
c
c.......accumulate integration point value of integrals
        fac=detj*w(l)
c ..... form kb stiffness matrix ..............
c ..... ok, now here goes the matrix multiplication that 
c ..... generates the bt d b ala book...
c ..... first bbt*db (12x3)*(3x3)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,3
            sum=0.d0
            do k=1,3
              sum=sum+bb(k,i)*db(j,k)
            enddo
            btdb(j,i)=sum
          enddo
        enddo
c then (bbt*db)*bb (12x12)*(3x12)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,12
            sum=0.d0
            do k=1,3
              sum=sum+btdb(k,i)*bb(k,j)
            enddo
            ekb(i,j)=ekb(i,j)+fac*sum
          enddo
        enddo
      enddo
c end of 2x2 integration loop
c
c.....begin 1x1 reduced integration loop
      do l=1,1
        xs(1)=0.d0
        xs(2)=0.d0
        call pgenshape(n,xs,xy,psi,dpsi,detj,dpsix,dpsiy)
        call ploadb(psi,dpsix,dpsiy,bb,bs)
c
c.......accumulate integration point value of integrals
        fac=detj*4.d0
c ..... form ks stiffness matrix ........................
c ..... ok, now here goes the matrix multiplication that 
c ..... generates bst*ds*bs ala book...
c ..... first bst*ds (12x2)*(2x2)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,2
            sum=0.d0
            do k=1,2
              sum=sum+bs(k,i)*ds(j,k)
            enddo
            btdb(j,i)=sum
          enddo
        enddo
c then (bbt*db)*bb (12x12)*(2x12)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,12
            sum=0.d0
            do k=1,2
              sum=sum+btdb(k,i)*bs(k,j)
            enddo
            eks(i,j)=eks(i,j)+fac*sum
          enddo
        enddo
      enddo
c end of 1x1 integration loop
c
c ... combine kb and ks stiffness matrices .......
      do i=1,12
        do j=1,12
c          print '(2i3,1p3g13.6)',i,j,ekb(i,j),eks(i,j)
          ekb(i,j)=ekb(i,j)+eks(i,j)
        enddo
      enddo
c      pause
c      if(.false.) then
c        do i=1,12
c          print 1003,(ekb(i,j)/ekb(i,i),j=1,12)
c        enddo
c      endif
1003  format(1x,1p6g13.6)
      return
1000  format(1x,'bad jacobian',e10.3)
1001  format(a,/,(2i5,3x,1pe13.6))
1002  format(a,/,(1p2e13.6))
      end
c=============================================================
      subroutine vassmbgk(n3,nz,gk0,ka,kz,ek,n,node)
      implicit real*8(a-h,o-z)
      dimension gk0(n3,nz)
      dimension ka(n3,nz+1),kz(n3)
      dimension ek(12,12),node(n)
c........................
c||||||||||||||||||||||||
c     print *,'in passmb'
      do l=1,n
        i=node(l)
        do m=1,n
          j=node(m)
c.........assemble global stiffness matrix gk
            if(i.eq.j) then
              gk0(i,1)=gk0(i,1)+ek(l,m)
              ka(i,1)=i
            else
              do k=2,kz(i)
                if(ka(i,k).eq.j) then
                  gk0(i,k)=gk0(i,k)+ek(l,m)
                  goto 99
                endif
              enddo
              kz(i)=kz(i)+1
              gk0(i,kz(i))=ek(l,m)
              ka(i,kz(i))=j
            endif
99          continue
        enddo
      enddo
c||||||||||||||||||||||||
c........................
      end
c=============================================================
      subroutine vformgk(nmax,n3,nz,numel,nnode,
     &                  x,y,kx,kkx,ek,xi,w,
     &                  gk0,ka,kz)
      implicit real*8(a-h,o-z)
      dimension x(nmax),y(nmax),kx(nmax,4),kkx(nmax,12)
      dimension gk0(n3,nz)
      dimension ka(n3,nz+1),kz(n3)
      dimension xi(2,4),w(4)
      dimension ek(12,12)
      dimension lm(4),xy(2,4),llm(12)
c      common /vmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c$doacross local(i,j)
      do i=1,nnode
        kz(i)=1
        ka(i,1)=i
        ka(i,nz+1)=0
        do j=1,nz
          ka(i,j)=0
          gk0(i,j)=0.0d0
        enddo
      enddo
      do nel=1,numel
        do l=1,4
          lm(l)=kx(nel,l)
          xy(1,l)=x(lm(l))
          xy(2,l)=y(lm(l))
        enddo
c$doacross local(l)
        do l=1,12
          llm(l)=kkx(nel,l)
        enddo
        call velemek(xy,4,ek,4,xi,w)
c..................................
        if(.false.) then
          write(7,*) nel
          do i=1,12
            write(7,*) (ek(i,j),j=1,12)
          enddo
c          pause
        endif
c..................................
        call vassmbgk(n3,nz,gk0,ka,kz,ek,12,llm)
c..................................
c        do i=1,nnode
c          print 1000,(gk0(i,j),j=1,nnode)
c        enddo
c        pause
c..................................
      enddo
c$doacross local(i)
      do i=1,nnode
c        print *,kz(i),ka(i,nz+1)
        ka(i,nz+1)=kz(i)
      enddo
c|||||||||||||||||||||||||||||||||||
c...................................
1000  format(1x,10f8.3)
      end
c=============================================================
      subroutine vformgf(nmax,n3,numel,nnode,
     &                  x,y,thick,www,kx,kkx,ef,xi,w,
     &                  gf,fnet,kode)
      implicit real*8(a-h,o-z)
      dimension x(nmax),y(nmax),kx(nmax,4),kkx(nmax,12)
      dimension gf(n3),www(n3)
      dimension thick(nmax),kode(nmax)
      dimension xi(2,4),w(4)
      dimension ef(12)
      dimension lm(4),xy(2,4),llm(12)
      common /vmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c$doacross local(i)
      do i=1,nnode
        gf(i)=0.d0
      enddo
      fnet=0.d0
      do nel=1,numel
        ethick=0.d0
        ewww=0.d0
        do l=1,4
          lm(l)=kx(nel,l)
          xy(1,l)=x(lm(l))
          xy(2,l)=y(lm(l))
          ethick=ethick+thick(lm(l))
          ewww=ewww+www((lm(l)-1)*3+1)
        enddo
c$doacross local(l)
        do l=1,12
          llm(l)=kkx(nel,l)
        enddo
        ethick=ethick/dble(4)
        ewww=ewww/dble(4)
c        print *,nel,real(ethick),real(ewww),
c     &          real(rhog*ethick+rockice*rhog*ewww)
        ethick=rhog*ethick
        ewww=rockice*rhog*ewww
c        print *,nel,real(ethick),real(ewww),real(ethick+ewww)
        ethick=ethick+ewww
        fnet=fnet+ethick
        call pelemgf(ethick,xy,4,ef,4,xi,w)
c..................................
        if(.false.) then
          write(7,*) nel
          do i=1,12
            write(7,*) ef(i)
          enddo
c          pause
        endif
c..................................
        call passmbgf(n3,gf,ef,12,llm)
c        call aplybc(n3,gf,nmax,kode,nnode)
c..................................
c        do i=1,nnode
c          print 1000,gf(i)
c        enddo
c        pause
c..................................
      enddo
c|||||||||||||||||||||||||||||||||||
c...................................
1000  format(1x,10f8.3)
      end
c=============================================================
      subroutine vloadgk(n3,nz,nnode,gk0,gk,ka,kz,
     &                    nit,nmax1,ia,ja,iwork,nmax,kode)
      implicit real*8(a-h,o-z)
c ....arrays for itpack sparse storage...................
      dimension gk(nit),ja(nit),ia(nmax1)
      dimension iwork(nit)
      dimension gk0(n3,nz)
      dimension ka(n3,nz+1),kz(n3),kode(nmax)
      data big /1d30/
      call sbini(nnode,nit,ia,ja,gk,iwork)
      do i = 1,nnode
c        n=1+(i-1)/3
c        if(kode(n).eq.1) then
c          call sbsij(nnode,nit,ia,ja,gk,iwork,i,i,big,0,1,6,ier)
c        else          
          do j = 1,kz(i)
            jg=ka(i,j)
            gkij=gk0(i,j)
            call sbsij(nnode,nit,ia,ja,gk,iwork,i,jg,gkij,0,1,6,ier)
          enddo 
c        endif
      enddo
      call sbend(nnode,nit,ia,ja,gk,iwork)
      end

c=======================================
      subroutine psetmat
c ... elastic material properties ...
      implicit real*8(a-h,o-z)                                          
      common /pmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c ... material properties and thickness of crust
c ... youngs modulus (nt/m**2)
      e=5.d10
      e=1.d11
c      e=1d20
c ... poisons ratio (dimensionless)
      rnu=0.5d0
c ... lames coefficients
      rlambda=rnu*e/(1.d0-rnu**2)
      rmu=0.5d0*e/(1.d0+rnu)
c ... thickness of crust (90-130 km) (110,000 m)
      ttt=110.d0*1000.d0
c      ttt=30*1000
c ... rhoice*grav (nt/m**3) times thickness, yields pressure, nt/m**2
      grav=9.8d0*0.3816d0
      rhor=5000.d0
      rhor=4000.d0
      rhow=1092.d0
      rhoi=917.d0
      rhog=rhoi*grav
      rockice=2.0d0*rhor/rhoi
c      rockice=1.0d0*rhor/rhoi
c ... relaxation time constant ......!
      relax=6000.d0                  !
      relax=500.d0                  !
      ctime=3.16516d8*relax/6000.d0  !
      end
c=======================================
      subroutine vsetmat
c ... viscous material properties ...
      implicit real*8(a-h,o-z)                                          
      common /vmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c ... material properties and thickness of crust
c ... youngs modulus (nt/m**2)
      e=5.d10
      e=1.d11
      e=1.d17
c ... poisons ratio (dimensionless)
      rnu=0.5d0
c ... lames coefficients
      rlambda=rnu*e/(1.d0-rnu**2)
      rmu=0.5d0*e/(1.d0+rnu)
c ... thickness of crust (90-130 km) (110,000 m)
      ttt=110.d0*1000.d0
       ttt=40*1000
c ... rhoice*grav (nt/m**3) times thickness, yields pressure, nt/m**2
      grav=9.8d0*0.3816d0
      rhor=5000.d0
      rhor=4000.d0
      rhow=1092.d0
      rhoi=917.d0
      rhog=rhoi*grav
      rockice=2.0d0*rhor/rhoi
c      rockice=1.0d0*rhor/rhoi
c ... relaxation time constant ......! NOT USED ...
      relax=6000.d0                  !
      ctime=3.16516d8*relax/6000.d0  !
      end
c=============================================================
      subroutine veplate(itime,numnp,numel,x,y,kx,thick,kode,
     &                 delt,www,wrate,time,wwworig,
     &                 wmin,wmine,wminv,fnete,fnetv)
c-----------------------------------------------------------------------
c 4th order visco-elastic plate solver with one-time generation of stiffness
c matrix. uses my sparse matrix storage for the static matrice
c and itpack sprse storage for the time-dependent modified matrices, and c             ipage=ipage+2

c iterative matrix solver. j fastook july, 1999
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)                                          
      include "parameter.h"
      parameter(nmax=maxnum,nz=27, nz1=nz+1, n3=3*nmax)
      parameter(nit=n3*nz,nmax1=n3+1,nmax3=n3*3)
      parameter(itmax=n3/10,ncg=4*itmax,nw=4*n3+ncg)
c ....arrays for itpack sparse storage...................
c      dimension gk(nit),ja(nit),ia(nmax1)
c      dimension iwksp(nmax3),wksp(nw),iwork(nit)
c      dimension iparm(12),rparm(12)
      dimension thick(nmax),x(nmax),y(nmax),kx(nmax,4)
      dimension www(n3),wrate(n3,2),wwworig(nmax)
      common /elastic/ wwwe(n3),werate(n3,2),wwweorig(nmax)
      common /viscous/ wwwv(n3),wvrate(n3,2),wwwvorig(nmax)
      dimension www0(n3)
      dimension kode(nmax)
c      dimension kkx(nmax,12),ltemp(nmax,3)
c      dimension xi(2,4),w(4)
c      dimension ek(12,12),ef(12)
c      dimension gf(n3)
c      dimension gk0(n3,nz)
c      dimension ka(n3,nz1),kz(n3)
c      dimension wwwtmp(n3)
      character*80 hed
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      common /line/ np,nline(1000)
      parameter(npage=39)
      character*80 list(1000)
      common /iolist/ list
      data www0 /n3*0.d0/
      logical file40
c ... timer stuff, for sgi only ...!
      real*4 tb(2)
c     real*4 tb(2),etime,dtime     !
c     external etime,dtime         !
c .................................!
      save ipass
c      save kkx,wsave,nnode,xi,w,gk0,ka,kz
c      save wwwsave,iparm,rparm,iread
c      save wwwe,werate,wwweorig
c      save wwwv,wvrate,wwwvorig
      data ipass /0/,iread /0/
      nnode=3*numnp
      if(ipass.eq.0) then
        ipass=1
        print *,'reading bedrock depression file'
        inquire(file='fort.40',exist=file40)
        iok=-1
        if(file40) then
          rewind 40
          read(40,*,iostat=iok) nnn
        endif
        if(iok.eq.0) then
          if(nnn.ne.nnode) then
            print *,'problems:'
            print *,'incompatible with current nnode=',nnode,nnn
            iok=1
          endif
          do i=1,nnode
            read(40,*) ii,wwwe(i),wwwv(i)
            werate(i,2)=wwwe(i)
            wvrate(i,2)=wwwv(i)
            if(i.ne.ii) then
              print *,'problems:reading wwwe,wwwv'
              print *,'incompatible with current (www)'
              iok=1
            endif
          enddo
          do i=1,nnode
            read(40,*) ii,werate(I,1),wvrate(I,1)
            if(i.ne.ii) then
              print *,'problems:reading wrates(1)'
              print *,'incompatible with current (www)'
              iok=1
            endif
          enddo
          do i=1,nnode
            read(40,*) ii,werate(I,2),wvrate(I,2)
            if(i.ne.ii) then
              print *,'problems:reading wrates(2)'
              print *,'incompatible with current (www)'
              iok=1
            endif
          enddo
          do i=1,numnp
            read(40,*) ii,wwweorig(i),wwwvorig(i)
            if(i.ne.ii) then
              print *,'problems:reading wwworigs'
              print *,'incompatible with current (wwworig)'
              iok=1
            endif
          enddo
          do i=1,numnp*3
            www(i)=wwwe(i)+wwwv(i)
c           if(mod(i,30).eq.1) print *,i,wwwe(i),wwwv(i)
            wrate(i,1)=werate(i,1)+wvrate(i,1)
            wrate(i,2)=werate(i,2)+wvrate(i,2)
          enddo
          do i=1,numnp
            wwworig(i)=wwweorig(i)+wwwvorig(i)
          enddo
          print *,' bedrock depression file found'
          print *,'    and read successfully '
          if(itime.lt.0) then
            print *,'abandoning unloading'
            rewind 40
            return
          endif
          iread=1
        endif
        if(iok.ne.0) then
          print *,' none found, set to zero ... (or part or www '
c$doacross local(i)
          if(.false.) then
            do i=1,nnode
              www(i)=0.d0
              wwwe(i)=0.d0
              wwwv(i)=0.d0
              wrate(i,2)=0.d0
              werate(i,2)=0.d0
              wvrate(i,2)=0.d0
            enddo
          elseif(.false.) then
            do i=1,nnode
              wwwe(i)=1.d0*www(i)
              wwwv(i)=0.d0*www(i)
              werate(i,2)=1.d0*wrate(i,2)
              wvrate(i,2)=0.d0*wrate(i,2)
            enddo
          else
            do i=1,numnp*3
              www(i)=wwwe(i)+wwwv(i)
c             if(mod(i,30).eq.1) print *,i,wwwe(i),wwwv(i)
              wrate(i,1)=werate(i,1)+wvrate(i,1)
              wrate(i,2)=werate(i,2)+wvrate(i,2)
            enddo
          endif
        endif
        write(hed,*) ' time=',nint(time-delt)
c        write(88) hed
c        write(88) (www(i),i=1,nnode,3)
c        write(88) (thick(i),i=1,numnp)
c        write(88) (1000*wrate(i,1),i=1,nnode,3)
      endif
      if(.true.) then
c ..... elastic solution .............................................
        call eplate(itime,numnp,numel,x,y,kx,thick,kode,
     &                 delt,wwwe,werate,wmine,time,wwweorig,fnete,wwwe)
      endif
      if(.true.) then
c ..... viscous solution .............................................
        call vplate(itime,numnp,numel,x,y,kx,thick,kode,
     &                 delt,wwwv,wvrate,wminv,time,wwwvorig,fnetv,wwwv)
      endif
      write(list(ipage+1),*) ' fnete,fnetv=',fnete*1d-6,fnetv*1d-6
      ipage=ipage+1
      do i=1,numnp*3
        www(i)=wwwe(i)+wwwv(i)
c        if(mod(i,30).eq.1) print *,i,wwwe(i),wwwv(i)
        wrate(i,1)=werate(i,1)+wvrate(i,1)
        wrate(i,2)=werate(i,2)+wvrate(i,2)
      enddo
c      do i=1,numnp
c        wwworig(i)=wwweorig(i)+wwwvorig(i)
c      enddo
c      pause
      wmin=1d30
      wmax=-1d30
      wmine=1d30
      wmaxe=-1d30
      wminv=1d30
      wmaxv=-1d30
      do i=1,numnp*3,3
        wmin=min(wmin,www(i))
        wmax=max(wmax,www(i))
        wmine=min(wmine,wwwe(i))
        wmaxe=max(wmaxe,wwwe(i))
        wminv=min(wminv,wwwv(i))
        wmaxv=max(wmaxv,wwwv(i))
      enddo
      if(.false.) then
        print *,'total  :',wmin,wmax
        print *,'elastic:',wmine,wmaxe
        print *,'viscous:',wminv,wmaxv
        do i=1,np
          ii=nline(i)
          print *,ii,real(wwwe((ii-1)*3+1)),
     &               real(wwwv((ii-1)*3+1)),
     &               real(www((ii-1)*3+1)),real(thick(ii))
        enddo
      endif
      end
C---------------------------------------------
      SUBROUTINE WRITEDEPN(NUMNP,NNODE)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      parameter(nmax=maxnum,n3=3*nmax)
      CHARACTER*1 CHAR
      common /elastic/ wwwe(n3),werate(n3,2),wwweorig(nmax)
      common /viscous/ wwwv(n3),wvrate(n3,2),wwwvorig(nmax)
      PRINT *,'IN WRITEDEPN',NUMNP,NNODE
      PRINT *,'   TO WRITE OUT BACKUP OF BEDROCK DEPRESSION '
      PRINT *,'   INPUT Y'
      READ(*,1002) CHAR
      WRITE(99,1002) CHAR
1002  FORMAT(A1)
      IF(CHAR.EQ.'Y' .OR. CHAR.EQ.'y') THEN
        REWIND 45
        WRITE(45,*) NNODE
        DO I=1,NNODE
          WRITE(45,*) I,wwwe(I),wwwv(I)
        ENDDO
        DO I=1,NNODE
          WRITE(45,*) I,werate(I,1),wvrate(I,1)
        ENDDO
        DO I=1,NNODE
          WRITE(45,*) I,werate(I,2),wvrate(I,2)
        ENDDO
        DO I=1,NUMNP
          WRITE(45,*) I,wwweorig(I),wwwvorig(i)
        ENDDO
      ENDIF
      END
C===================================================
      SUBROUTINE UNLOAD(NUMNP,BDROCK,UNDEPB,PSURF,RHOI,RHOR)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      DIMENSION BDROCK(NUMNP),UNDEPB(NUMNP),PSURF(NUMNP)
      RETURN
      IF(BTOGG.ne.0) RETURN
      DO N=1,NUMNP
        IF(BDROCK(N).GT.-9999.) THEN
          UNDEPB(N)=RHOI*(PSURF(N)-BDROCK(N))/RHOR+BDROCK(N)
        ENDIF
      ENDDO
      END
C===================================================
      SUBROUTINE SPLATE(ITIME,NUMNP,THICK,UNDEPB,DEPB,RHOI,RHOR,RHOW,
     &                  DELT,TIME,WRATE,WWW,WWWORIG)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM,N3=NMAX*3)
      PARAMETER(RELAX=3000.D0,RK=1.D0/RELAX)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      DIMENSION UNDEPB(NMAX),DEPB(NMAX),THICK(NMAX)
      DIMENSION WWW(N3),WRATE(N3,2),WWWORIG(NMAX)
      CHARACTER*80 HED
      COMMON /SELASTIC/ WWWS(N3),WSRATE(N3,2),WWWSORIG(NMAX)
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      logical file40
      DATA IPASS /0/, IREAD /0/
      SAVE IPASS,IREAD,RATIO,SIXTH,ROCKICE
      NNODE=3*NUMNP
      IF(IPASS.EQ.0) THEN
C ....  STUFF YOU ONLY DO THE FIRST PASS ..............
        SIXTH=1D0/6D0
        RATIO=RHOI/RHOR
        ROCKICE=RHOR/RHOI
        PRINT *,'READING BEDROCK DEPRESSION FILE'
        inquire(file='fort.40',exist=file40)
        iok=-1
        if(file40) then
          READ(40,*,IOSTAT=IOK) NNN
        endif
        IF(IOK.EQ.0) THEN
          IF(NNN.NE.NNODE) THEN
            PRINT *,'PROBLEMS:'
            PRINT *,'INCOMPATIBLE WITH CURRENT NNODE=',NNODE,NNN
            IOK=1
          ENDIF
          DO I=1,NNODE
            READ(40,*) II,WWW(I),WJUNK
            IF(I.NE.II) THEN
              PRINT *,'PROBLEMS:READING WWW'
              PRINT *,'INCOMPATIBLE WITH CURRENT (WWW)'
              IOK=1
            ENDIF
          ENDDO
          DO I=1,NNODE
            READ(40,*) II,WRATE(I,1),WJUNK
            IF(I.NE.II) THEN
              PRINT *,'PROBLEMS:READING WRATES(1)'
              PRINT *,'INCOMPATIBLE WITH CURRENT (WWW)'
              IOK=1
            ENDIF
          ENDDO
          DO I=1,NNODE
            READ(40,*) II,WRATE(I,2),WJUNK
            IF(I.NE.II) THEN
              PRINT *,'PROBLEMS:READING WRATES(2)'
              PRINT *,'INCOMPATIBLE WITH CURRENT (WWW)'
              IOK=1
            ENDIF
          ENDDO
          DO I=1,NUMNP
            READ(40,*) II,WWWORIG(I),WJUNK
            WWWSORIG(I)=WWWORIG(i)
            IF(I.NE.II) THEN
              PRINT *,'PROBLEMS:READING WWWORIGS'
              PRINT *,'INCOMPATIBLE WITH CURRENT (WWWORIG)'
              IOK=1
            ENDIF
          ENDDO
          PRINT *,' BEDROCK DEPRESSION FILE FOUND'
          PRINT *,'    AND READ SUCCESSFULLY '
          IF(ITIME.LT.0) THEN
            PRINT *,'ABANDONING UNLOADING'
            REWIND 40
            RETURN
          ENDIF
          IREAD=1
        ENDIF
        IF(IOK.NE.0) THEN
          PRINT *,' NONE FOUND, SET TO ZERO ... AND UNLOADED ...'
C$DOACROSS LOCAL(I)
          DO I=1,NNODE
            WWW(I)=0.D0
            WRATE(I,1)=0.D0
            WRATE(I,2)=0.D0
          ENDDO
          IF(DELT.LE.0.0) THEN
            II=1
            DO I=1,NUMNP
              WWW(II)=-RATIO*THICK(I)
              WRATE(II,2)=WWW(II)
              II=II+3
            ENDDO
          ENDIF
          WRITE(HED,*) ' TIME=',NINT(TIME-DELT)
          WRITE(88) HED
          WRITE(88) (WWW(I),I=1,NNODE,3)
          WRITE(88) (THICK(I),I=1,NUMNP)
          WRITE(88) (1000*WRATE(I,1),I=1,NNODE,3)
        ENDIF
      ENDIF
C ... END OF STUFF DONE ONLY ONCE ...
      II=1
      DO I=1,NUMNP
        RHS1=-RK*(WWW(II)           +RATIO*THICK(I))*DELT
        RHS2=-RK*(WWW(II)+0.5D0*RHS1+RATIO*THICK(I))*DELT
        RHS3=-RK*(WWW(II)+0.5D0*RHS2+RATIO*THICK(I))*DELT
        RHS4=-RK*(WWW(II)+RHS3      +RATIO*THICK(I))*DELT
        SUM=(RHS1+RHS2+RHS2+RHS3+RHS3+RHS4)*SIXTH
        WWW(II)=WWW(II)+SUM
        WRATE(II,1)=SUM/DELT
c        DEPB(I)=UNDEPB(I)+WWW(II)
C        IF(THICK(I).NE.0.) THEN
C          PRINT 100,I,THICK(I),WWW(II),WRATE(II,1),DEPB(I),UNDEPB(I)
C        ENDIF
        II=II+3
      ENDDO
100   FORMAT(1X,I5,5G13.6)
      IF(DELT.GT.0. .and. .false.) THEN
        DO I=1,NNODE
          WRATE(I,1)=(WWW(I)-WRATE(I,2))/DELT
        ENDDO
      ENDIF
      DO I=1,NNODE
        WRATE(I,2)=WWW(I)
      ENDDO
      WMIN=1D30
      WMAX=-1D30
      WTOT=0.D0
      DO I=1,NNODE,3
        WMIN=MIN(WMIN,WWW(I))
        WMAX=MAX(WMAX,WWW(I))
        WTOT=WTOT+WWW(I)
C        WRITE(*,*) (WWW(I+J),J=0,2)
        WWWS(I)=WWW(I)
        WSRATE(I,1)=WRATE(I,1)
        WSRATE(I,2)=WRATE(I,2)
      ENDDO
      THTOT=0.D0
      DO I=1,NUMNP
        THTOT=THTOT+THICK(I)
      ENDDO
      IF(IREAD.EQ.0 .AND. IPASS.EQ.0) THEN
        II=1
        DO I=1,NUMNP
          WWWORIG(I)=WWW(II)
          WWWSORIG(I)=WWW(II)
          II=II+3
        ENDDO
      ENDIF
      IF(IOTOGG) THEN
        WRITE(LIST(IPAGE+1),*) 
     &       '*******************************************'
        WRITE(LIST(IPAGE+2),1001) TIME,REAL(THTOT),
     &              REAL(-THTOT/WTOT/ROCKICE),
     &              REAL(WMIN),REAL(WMAX),REAL(WSAVE-WMIN)
        WRITE(LIST(IPAGE+3),*) 
     &       '*******************************************'
        IPAGE=IPAGE+3
      ENDIF
      WRITE(92,*) TIME,WMIN
      WSAVE=WMIN
      IPASS=1
1001  FORMAT(1X,6(1X,1PG12.5))
      END
        

C---------------------------------------------
      SUBROUTINE WRITEDEPS(NUMNP,NNODE)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM,N3=3*NMAX)
      CHARACTER*1 CHAR
      COMMON /SELASTIC/ WWW(N3),WRATE(N3,2),WWWORIG(NMAX)
      PRINT *,'IN WRITEDEP',NNODE
      PRINT *,'   TO WRITE OUT BACKUP OF BEDROCK DEPRESSION '
      PRINT *,'   INPUT Y'
      READ(*,1002) CHAR
      WRITE(99,1002) CHAR
1002  FORMAT(A1)
      IF(CHAR.EQ.'Y' .OR. CHAR.EQ.'y') THEN
        REWIND 45
        WRITE(45,*) NNODE
        DO I=1,NNODE
          WRITE(45,*) I,WWW(I),0.D0
        ENDDO
        DO I=1,NNODE
          WRITE(45,*) I,WRATE(I,1),0.D0
        ENDDO
        DO I=1,NNODE
          WRITE(45,*) I,WRATE(I,2),0.D0
        ENDDO
        DO I=1,NUMNP
          WRITE(45,*) I,WWWORIG(I),0.D0
        ENDDO
      ENDIF
      END
c==========================================================
      SUBROUTINE JCG (NN,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM,IERR)
C       
C     ITPACK 2C MAIN SUBROUTINE  JCG  (JACOBI CONJUGATE GRADIENT)   
C     EACH OF THE MAIN SUBROUTINES:   
C           JCG, JSI, SOR, SSORCG, SSORSI, RSCG, RSSI     
C     CAN BE USED INDEPENDENTLY OF THE OTHERS   
C       
C ... FUNCTION:   
C       
C          THIS SUBROUTINE, JCG, DRIVES THE JACOBI CONJUGATE
C          GRADIENT ALGORITHM.
C       
C ... PARAMETER LIST:       
C       
C          N      INPUT INTEGER.  DIMENSION OF THE MATRIX. (= NN)   
C          IA,JA  INPUT INTEGER VECTORS.  THE TWO INTEGER ARRAYS OF 
C                 THE SPARSE MATRIX REPRESENTATION.       
C          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
C                 MATRIX REPRESENTATION.
C          RHS    INPUT D.P. VECTOR.  CONTAINS THE RIGHT HAND SIDE  
C                 OF THE MATRIX PROBLEM.
C          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U CONTAINS THE 
C                 INITIAL GUESS TO THE SOLUTION. ON OUTPUT, IT CONTAINS 
C                 THE LATEST ESTIMATE TO THE SOLUTION.    
C          IWKSP  INTEGER VECTOR WORKSPACE OF LENGTH 3*N  
C          NW     INPUT INTEGER.  LENGTH OF AVAILABLE WKSP.  ON OUTPUT, 
C                 IPARM(8) IS AMOUNT USED.      
C          WKSP   D.P. VECTOR USED FOR WORKING SPACE.  JACOBI CONJUGATE 
C                 GRADIENT NEEDS THIS TO BE IN LENGTH AT LEAST      
C                 4*N + 2*ITMAX,  IF ISYM = 0  (SYMMETRIC STORAGE)  
C                 4*N + 4*ITMAX,  IF ISYM = 1  (NONSYMMETRIC STORAGE) 
C                 HERE ITMAX = IPARM(1) AND ISYM = IPARM(5) 
C                 (ITMAX IS THE MAXIMUM ALLOWABLE NUMBER OF ITERATIONS) 
C          IPARM  INTEGER VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY
C                 SOME INTEGER PARAMETERS WHICH AFFECT THE METHOD.  
C          RPARM  D.P. VECTOR OF LENGTH 12. ALLOWS USER TO SPECIFY SOME 
C                 D.P. PARAMETERS WHICH AFFECT THE METHOD.
C          IER    OUTPUT INTEGER.  ERROR FLAG. (= IERR)   
C       
C ... JCG SUBPROGRAM REFERENCES:      
C       
C          FROM ITPACK    BISRCH, CHGCON, DETERM, DFAULT, ECHALL,   
C                         ECHOUT, EIGVNS, EIGVSS, EQRT1S, ITERM, TIMER, 
C                         ITJCG, IVFILL, PARCON, PERMAT,  
C                         PERROR1, PERVEC, PJAC, PMULT, PRBNDX,      
C                         PSTOP, QSORT, DAXPY, SBELM, SCAL, DCOPY,  
C                         DDOT, SUM3, UNSCAL, VEVMW, VFILL, VOUT,   
C                         WEVMW, ZBRENT 
C          SYSTEM         DABS, DLOG10, DBLE(AMAX0), DMAX1, MOD, DSQRT
C       
C     VERSION:  ITPACK 2C (MARCH 1982)
C       
C     CODE WRITTEN BY:  DAVID KINCAID, ROGER GRIMES, JOHN RESPESS   
C                       CENTER FOR NUMERICAL ANALYSIS     
C                       UNIVERSITY OF TEXAS     
C                       AUSTIN, TX  78712       
C                       (512) 471-1242
C       
C     FOR ADDITIONAL DETAILS ON THE   
C          (A) SUBROUTINE SEE TOMS ARTICLE 1982 
C          (B) ALGORITHM  SEE CNA REPORT 150    
C       
C     BASED ON THEORY BY:  DAVID YOUNG, DAVID KINCAID, LOU HAGEMAN  
C       
C     REFERENCE THE BOOK:  APPLIED ITERATIVE METHODS      
C                          L. HAGEMAN, D. YOUNG 
C                          ACADEMIC PRESS, 1981 
C       
C     **************************************************  
C     *               IMPORTANT NOTE                   *  
C     *                                                *  
C     *      WHEN INSTALLING ITPACK ROUTINES ON A      *  
C     *  DIFFERENT COMPUTER, RESET SOME OF THE VALUES  *  
C     *  IN  SUBROUTNE DFAULT.   MOST IMPORTANT ARE    *  
C     *                                                *  
C     *   DRELPR      MACHINE RELATIVE PRECISION       *  
C     *   RPARM(1)    STOPPING CRITERION               *  
C     *                                                *  
C     *   ALSO CHANGE SYSTEM-DEPENDENT ROUTINE         *  
C     *   SECOND USED IN TIMER                         *  
C     *                                                *  
C     **************************************************  
C       
C     SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(*),JA(*),IWKSP(*),IPARM(12),NN,NW,IERR   
      DOUBLE PRECISION A(*),RHS(NN),U(NN),WKSP(NW),RPARM(12)
C       
C     SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IB1,IB2,IB3,IB4,IB5,IDGTS,IER,IERPER,ITMAX1,LOOP,N,NB,N3
      DOUBLE PRECISION DIGIT1,DIGIT2,TEMP,TIME1,TIME2,TOL 
C       
C **** BEGIN: ITPACK COMMON 
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C **** END  : ITPACK COMMON 
C       
C ... VARIABLES IN COMMON BLOCK - ITCOM1
C       
C     IN     - ITERATION NUMBER       
C     IS     - ITERATION NUMBER WHEN PARAMETERS LAST CHANGED
C     ISYM   - SYMMETRIC/NONSYMMETRIC STORAGE FORMAT SWITCH 
C     ITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED       
C     LEVEL  - LEVEL OF OUTPUT CONTROL SWITCH   
C     NOUT   - OUTPUT UNIT NUMBER     
C       
C ... VARIABLES IN COMMON BLOCK - ITCOM2
C       
C     ADAPT  - FULLY ADAPTIVE PROCEDURE SWITCH  
C     BETADT - SWITCH FOR ADAPTIVE DETERMINATION OF BETA  
C     CASEII - ADAPTIVE PROCEDURE CASE SWITCH   
C     HALT   - STOPPING TEST SWITCH   
C     PARTAD - PARTIALLY ADAPTIVE PROCEDURE SWITCH
C       
C ... VARIABLES IN COMMON BLOCK - ITCOM3
C       
C     BDELNM - TWO NORM OF B TIMES DELTA-SUPER-N
C     BETAB  - ESTIMATE FOR THE SPECTRAL RADIUS OF LU MATRIX
C     CME    - ESTIMATE OF LARGEST EIGENVALUE   
C     DELNNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION N      
C     DELSNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION S      
C     FF     - ADAPTIVE PROCEDURE DAMPING FACTOR
C     GAMMA  - ACCELERATION PARAMETER 
C     OMEGA  - OVERRELAXATION PARAMETER FOR SOR AND SSOR  
C     QA     - PSEUDO-RESIDUAL RATIO  
C     QT     - VIRTUAL SPECTRAL RADIUS
C     RHO    - ACCELERATION PARAMETER 
C     RRR    - ADAPTIVE PARAMETER     
C     SIGE   - PARAMETER SIGMA-SUB-E  
C     SME    - ESTIMATE OF SMALLEST EIGENVALUE  
C     SPECR  - SPECTRAL RADIUS ESTIMATE FOR SSOR
C     DRELPR - MACHINE RELATIVE PRECISION       
C     STPTST - STOPPING PARAMETER     
C     UDNM   - TWO NORM OF U
C     ZETA   - STOPPING CRITERION     
C       
C ... INITIALIZE COMMON BLOCKS
C       
      LEVEL = IPARM(2)      
      NOUT = IPARM(4)       
      IF (LEVEL.GE.1) WRITE (NOUT,10) 
   10 FORMAT ('0'///1X,'BEGINNING OF ITPACK SOLUTION MODULE  JCG')  
      IER = 0     
      IF (IPARM(1).LE.0) RETURN       
      N = NN      
      IF (IPARM(11).EQ.0) TIMJ1 = TIMER(DUMMY)  
      IF (LEVEL.GE.3) GO TO 20
      CALL ECHOUT (IPARM,RPARM,1)     
      GO TO 30    
   20 CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,1) 
   30 TEMP = 5.0D2*DRELPR   
      IF (ZETA.GE.TEMP) GO TO 50      
      IF (LEVEL.GE.1) WRITE (NOUT,40) ZETA,DRELPR,TEMP    
   40 FORMAT ('0','*** W A R N I N G ************'/'0',   
     *   '    IN ITPACK ROUTINE JCG'/' ','    RPARM(1) =',D10.3,    
     *   ' (ZETA)'/' ','    A VALUE THIS SMALL MAY HINDER CONVERGENCE '/
     *   ' ','    SINCE MACHINE PRECISION DRELPR =',D10.3/' ',      
     *   '    ZETA RESET TO ',D10.3)  
      ZETA = TEMP 
   50 CONTINUE    
      TIME1 = RPARM(9)      
      TIME2 = RPARM(10)     
      DIGIT1 = RPARM(11)    
      DIGIT2 = RPARM(12)    
C       
C ... VERIFY N    
C       
      IF (N.GT.0) GO TO 70  
      IER = 11    
      IF (LEVEL.GE.0) WRITE (NOUT,60) N 
   60 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    CALLED FROM ITPACK ROUTINE JCG '/' ',       
     *   '    INVALID MATRIX DIMENSION, N =',I8)
      GO TO 370   
   70 CONTINUE    
C       
C ... REMOVE ROWS AND COLUMNS IF REQUESTED      
C       
      IF (IPARM(10).EQ.0) GO TO 90    
      TOL = RPARM(8)
      CALL IVFILL (N,IWKSP,0) 
      CALL VFILL (N,WKSP,0.0D0)       
      CALL SBELM (N,IA,JA,A,RHS,IWKSP,WKSP,TOL,ISYM,LEVEL,NOUT,IER) 
      IF (IER.EQ.0) GO TO 90
      IF (LEVEL.GE.0) WRITE (NOUT,80) IER,TOL   
   80 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    CALLED FROM ITPACK ROUTINE JCG '/' ',       
     *   '    ERROR DETECTED IN SUBROUTINE  SBELM '/' ',  
     *   '    WHICH REMOVES ROWS AND COLUMNS OF SYSTEM '/' ',       
     *   '    WHEN DIAGONAL ENTRY TOO LARGE  '/' ','    IER = ',I5,5X,
     *   ' RPARM(8) = ',D10.3,' (TOL)') 
      GO TO 370   
C       
C ... INITIALIZE WKSP BASE ADDRESSES. 
C       
   90 IB1 = 1     
      IB2 = IB1+N 
      IB3 = IB2+N 
      IB4 = IB3+N 
      IB5 = IB4+N 
      IPARM(8) = 4*N+2*ITMAX
      IF (ISYM.NE.0) IPARM(8) = IPARM(8)+2*ITMAX
      IF (NW.GE.IPARM(8)) GO TO 110   
      IER = 12    
      IF (LEVEL.GE.0) WRITE (NOUT,100) NW,IPARM(8)
  100 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    CALLED FROM ITPACK ROUTINE JCG '/' ',       
     *   '    NOT ENOUGH WORKSPACE AT ',I10/' ','    SET IPARM(8) =',I10
     *   ,' (NW)')
      GO TO 370   
C       
C ... PERMUTE TO  RED-BLACK SYSTEM IF REQUESTED 
C       
  110 NB = IPARM(9) 
      IF (NB.LT.0) GO TO 170
      N3 = 3*N    
      CALL IVFILL (N3,IWKSP,0)
      CALL PRBNDX (N,NB,IA,JA,IWKSP,IWKSP(IB2),LEVEL,NOUT,IER)      
      IF (IER.EQ.0) GO TO 130 
      IF (LEVEL.GE.0) WRITE (NOUT,120) IER,NB   
  120 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    CALLED FROM ITPACK ROUTINE JCG '/' ',       
     *   '    ERROR DETECTED IN SUBROUTINE  PRBNDX'/' ',  
     *   '    WHICH COMPUTES THE RED-BLACK INDEXING'/' ','    IER = ',I5
     *   ,' IPARM(9) = ',I5,' (NB)')  
      GO TO 370   
C       
C ... PERMUTE MATRIX AND RHS
C       
  130 IF (LEVEL.GE.2) WRITE (NOUT,140) NB       
  140 FORMAT (/10X,'ORDER OF BLACK SUBSYSTEM = ',I5,' (NB)')
      CALL PERMAT (N,IA,JA,A,IWKSP,IWKSP(IB3),ISYM,LEVEL,NOUT,IER)  
      IF (IER.EQ.0) GO TO 160 
      IF (LEVEL.GE.0) WRITE (NOUT,150) IER      
  150 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    CALLED FROM ITPACK ROUTINE JCG '/' ',       
     *   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  
     *   '    WHICH DOES THE RED-BLACK PERMUTATION'/' ','    IER = ',I5)
      GO TO 370   
  160 CALL PERVEC (N,RHS,IWKSP)       
      CALL PERVEC (N,U,IWKSP) 
C       
C ... SCALE LINEAR SYSTEM, U, AND RHS BY THE SQUARE ROOT OF THE     
C ... DIAGONAL ELEMENTS.    
C       
  170 CONTINUE    
      CALL VFILL (IPARM(8),WKSP,0.0D0)
      CALL SCAL (N,IA,JA,A,RHS,U,WKSP,LEVEL,NOUT,IER)     
      IF (IER.EQ.0) GO TO 190 
      IF (LEVEL.GE.0) WRITE (NOUT,180) IER      
  180 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    CALLED FROM ITPACK ROUTINE JCG '/' ',       
     *   '    ERROR DETECTED IN SUBROUTINE  SCAL  '/' ',  
     *   '    WHICH SCALES THE SYSTEM   '/' ','    IER = ',I5)      
      GO TO 370   
  190 IF (LEVEL.LE.2) GO TO 220       
      WRITE (NOUT,200)      
  200 FORMAT (///1X,'IN THE FOLLOWING, RHO AND GAMMA ARE',
     *   ' ACCELERATION PARAMETERS')  
      IF (ADAPT) WRITE (NOUT,210)     
  210 FORMAT (1X,'CME IS THE ESTIMATE OF THE LARGEST EIGENVALUE OF',
     *   ' THE JACOBI MATRIX')
  220 IF (IPARM(11).NE.0) GO TO 230   
      TIMI1 = TIMER(DUMMY)  
C       
C ... COMPUTE INITIAL PSEUDO-RESIDUAL 
C       
  230 CONTINUE    
      CALL DCOPY (N,RHS,1,WKSP(IB2),1)
      CALL PJAC (N,IA,JA,A,U,WKSP(IB2)) 
      CALL VEVMW (N,WKSP(IB2),U)      
C       
C ... ITERATION SEQUENCE    
C       
      ITMAX1 = ITMAX+1      
      DO 250 LOOP = 1,ITMAX1
         IN = LOOP-1
         IF (MOD(IN,2).EQ.1) GO TO 240
C       
C ... CODE FOR THE EVEN ITERATIONS.   
C       
C     U           = U(IN)             WKSP(IB2) = DEL(IN) 
C     WKSP(IB1)   = U(IN-1)           WKSP(IB3) = DEL(IN-1) 
C       
         CALL ITJCG (N,IA,JA,A,U,WKSP(IB1),WKSP(IB2),WKSP(IB3),WKSP(IB4)
     *      ,WKSP(IB5))     
C       
         IF (HALT) GO TO 280
         GO TO 250
C       
C ... CODE FOR THE ODD ITERATIONS.    
C       
C     U           = U(IN-1)           WKSP(IB2) = DEL(IN-1) 
C     WKSP(IB1)   = U(IN)             WKSP(IB3) = DEL(IN) 
C       
  240    CALL ITJCG (N,IA,JA,A,WKSP(IB1),U,WKSP(IB3),WKSP(IB2),WKSP(IB4)
     *      ,WKSP(IB5))     
C       
         IF (HALT) GO TO 280
  250 CONTINUE    
C       
C ... ITMAX HAS BEEN REACHED
C       
      IF (IPARM(11).NE.0) GO TO 260   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  260 IER = 13    
      IF (LEVEL.GE.1) WRITE (NOUT,270) ITMAX    
  270 FORMAT ('0','*** W A R N I N G ************'/'0',   
     *   '    IN ITPACK ROUTINE JCG'/' ','    FAILURE TO CONVERGE IN',I5
     *   ,' ITERATIONS')    
      IF (IPARM(3).EQ.0) RPARM(1) = STPTST      
      GO TO 310   
C       
C ... METHOD HAS CONVERGED  
C       
  280 IF (IPARM(11).NE.0) GO TO 290   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  290 IF (LEVEL.GE.1) WRITE (NOUT,300) IN       
  300 FORMAT (/1X,'JCG  HAS CONVERGED IN ',I5,' ITERATIONS')
C       
C ... PUT SOLUTION INTO U IF NOT ALREADY THERE. 
C       
  310 CONTINUE    
      IF (MOD(IN,2).EQ.1) CALL DCOPY (N,WKSP(IB1),1,U,1)  
C       
C ... UNSCALE THE MATRIX, SOLUTION, AND RHS VECTORS.      
C       
      CALL UNSCAL (N,IA,JA,A,RHS,U,WKSP)
C       
C ... UN-PERMUTE MATRIX,RHS, AND SOLUTION       
C       
      IF (IPARM(9).LT.0) GO TO 340    
      CALL PERMAT (N,IA,JA,A,IWKSP(IB2),IWKSP(IB3),ISYM,LEVEL,NOUT, 
     *   IERPER)  
      IF (IERPER.EQ.0) GO TO 330      
      IF (LEVEL.GE.0) WRITE (NOUT,320) IERPER   
  320 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    CALLED FROM ITPACK ROUTINE JCG '/' ',       
     *   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  
     *   '    WHICH UNDOES THE RED-BLACK PERMUTATION   '/' ',       
     *   '    IER = ',I5)   
      IF (IER.EQ.0) IER = IERPER      
      GO TO 370   
  330 CALL PERVEC (N,RHS,IWKSP(IB2))  
      CALL PERVEC (N,U,IWKSP(IB2))    
C       
C ... OPTIONAL ERROR ANALYSIS 
C       
  340 IDGTS = IPARM(12)     
      IF (IDGTS.LT.0) GO TO 350       
      IF (IPARM(2).LE.0) IDGTS = 0    
      CALL PERROR1 (N,IA,JA,A,RHS,U,WKSP,DIGIT1,DIGIT2,IDGTS)
C       
C ... SET RETURN PARAMETERS IN IPARM AND RPARM  
C       
  350 IPARM(8) = IPARM(8)-2*(ITMAX-IN)
      IF (IPARM(11).NE.0) GO TO 360   
      TIMJ2 = TIMER(DUMMY)  
      TIME2 = DBLE(TIMJ2-TIMJ1)       
  360 IF (ISYM.NE.0) IPARM(8) = IPARM(8)-2*(ITMAX-IN)     
      IF (IPARM(3).NE.0) GO TO 370    
      IPARM(1) = IN 
      IPARM(9) = NB 
      RPARM(2) = CME
      RPARM(3) = SME
      RPARM(9) = TIME1      
      RPARM(10) = TIME2     
      RPARM(11) = DIGIT1    
      RPARM(12) = DIGIT2    
C       
  370 CONTINUE    
      IERR = IER  
      IF (LEVEL.GE.3) CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,2)     
C       
      RETURN      
      END 
c==================================================================
      SUBROUTINE SBINI (N,NZ,IA,JA,A,IWORK)     
C       
C***********************************************************************
C       
C     SBINI IS THE FIRST OF A SUITE OF THREE SUBROUTINES TO AID     
C     THE USER TO CONSTRUCT THE IA, JA, A DATA STRUCTURE USED       
C     IN ITPACK.  
C       
C     SBINI INITIALIZES THE ARRAYS IA, JA, IWORK, AND A.  THE OTHER 
C     SUBROUTINES IN THE SUITE ARE SBSIJ ( WHICH BUILDS A LINKED    
C     LIST REPRESENTATION OF THE MATRIX STRUCTURE ) AND SBEND ( WHICH 
C     RESTRUCTURE THE LINKED LIST FORM INTO THE FINAL FORM ).       
C       
C ... PARAMETERS  
C       
C ...... INPUT    
C       
C     N          THE ORDER OF THE LINEAR SYSTEM 
C       
C     NZ         THE MAXIMUM NUMBER OF NONZEROES ALLOWED IN THE     
C                LINEAR SYSTEM.       
C       
C ...... OUTPUT   
C       
C     IA         INTEGER ARRAY OF LENGTH N+1.  SBINI SETS THIS ARRAY
C                TO -I FOR I = 1 THRU N.  IA(N+1) IS SET TO NZ.     
C       
C     JA         INTEGER ARRAY OF LENGTH NZ.  INITIALIZED TO ZERO HERE. 
C       
C     A          D.P. ARRAY OF LENGTH NZ.  INITIALIZED TO ZERO HERE.
C       
C     IWORK       INTEGER ARRAY OF LENGTH NZ.  INITIALIZED TO ZERO HERE.
C       
C***********************************************************************
C       
      INTEGER N,NZ,IA(NZ+1),JA(NZ),IWORK(NZ),I     
      DOUBLE PRECISION A(NZ)
C       
C***********************************************************************
C       
      DO 10 I = 1,N 
         IA(I) = -I 
   10 CONTINUE    
      IA(N+1) = NZ
C       
      CALL IVFILL (NZ,JA,0) 
      CALL IVFILL (NZ,IWORK,0)
      CALL VFILL (NZ,A,0.D0)
C       
      RETURN      
      END 
C======================================================================
      SUBROUTINE SBSIJ (N,NZ,IA,JA,A,IWORK,II,JJ,VALL,MODE,LEVELL,NOUTT,
     *   IERR)    
C       
C***********************************************************************
C       
C     SBSIJ IS THE SECOND OF A SUITE OF THREE SUBROUTINES TO AID IN 
C     THE CONSTRUCTION OF THE IA, JA, A DATA STRUCTURE USED IN      
C     ITPACK.     
C       
C     SBSIJ TAKES THE INDIVIDUAL ENTRIES OF THE SPARSE MATRIX AS    
C     GIVEN TO IT AT EACH CALL VIA  (I,J,VAL) AND INSERTS IT INTO   
C     A LINKED LIST REPRESENTATION OF THE SPARSE MATRIX.  
C       
C     EACH ROW OF THE SPARSE MATRIX IS ASSOCIATED WITH A CIRCULAR   
C     LINKED LIST BEGINNING AT IA(I).  THE LAST ENTERED ELEMENT IN  
C     EACH LIST POINTS BACK TO IA(I) WITH THE VALUE -I.  THE LINKS  
C     ARE STORED IN THE ARRAY IWORK, WHILE JA AND A STORE THE COLUMN
C     NUMBER AND VALUE IN PARALLEL TO IWORK.  THE LINKED LISTED ARE 
C     STORED BEGINNING AT ENTRY NZ AND WORKING BACKWARDS TOWARDS 1. 
C       
C ... PARAMETERS  
C       
C ...... INPUT    
C       
C     N       THE ORDER OF THE LINEAR SYSTEM    
C       
C     NZ      THE LENGTH OF THE ARRAYS  JA, A, AND IWORK  
C       
C     I, J    THE ROW AND COLUMN NUMBERS OF THE ENTRY OF THE SPARSE 
C             LINEAR SYSTEM TO BE ENTERED IN THE DATA STRUCTURE(=II,JJ) 
C       
C     VAL     THE NONZERO VALUE ASSOCIATED WITH (I,J)  (= VALL)     
C       
C     MODE    IF THE (I,J) ENTRY HAS ALREADY BEEN SET, MODE SPECIFIES 
C             THE WAY IN WHICH THE ENTRY IS TO BE TREATED.
C             IF   MODE .LT. 0  LET THE VALUE REMAIN AS IS
C                       .EQ. 0  RESET IT TO THE NEW VALUE 
C                       .GT. 0  ADD THE NEW VALUE TO THE OLD VALUE  
C       
C     NOUT  OUTPUT FILE NUMBER (= NOUTT)
C       
C     LEVEL   OUTPUT FILE SWITCH (= LEVELL)     
C ... INPUT/OUTPUT
C       
C     IA      INTEGER ARRAY OF LENGTH N+1.  THE FIRST N ENTRIES     
C             POINT TO THE BEGINNING OF THE LINKED LIST FOR EACH    
C             ROW.  IA(N+1) POINTS TO THE NEXT ENTRY AVAILABLE FOR  
C             STORING THE CURRENT ENTRY INTO THE LINKED LIST.       
C       
C     JA      INTEGER ARRAY OF LENGTH NZ.  JA STORES THE COLUMN     
C             NUMBERS OF THE NONZERO ENTRIES.   
C       
C     A       D.P. ARRAY OF LENGTH NZ.  A STORES THE VALUE OF THE   
C             NONZERO ENTRIES.
C       
C     IWORK   INTEGER ARRAY OF LENGTH NZ. IWORK STORES THE LINKS.   
C       
C     IER     ERROR FLAG.(= IERR)  POSSIBLE RETURNS ARE   
C             IER =    0   SUCCESSFUL COMPLETION
C                 =  700   ENTRY WAS ALREADY SET,  VALUE HANDLED    
C                          AS SPECIFIED BY MODE.
C                 =  701   IMPROPER VALUE OF EITHER I OR J INDEX    
C                 =  702   NO ROOM REMAINING, NZ TOO SMALL. 
C       
C***********************************************************************
C       
      INTEGER N,NZ,IA(*),JA(NZ),IWORK(NZ),II,JJ,MODE,LEVELL,NOUTT,
     &        IERR
      DOUBLE PRECISION A(NZ),VALL     
C       
      INTEGER LINK,NEXT,NPL1,I,J,LEVEL,NOUT,IER 
      DOUBLE PRECISION VAL,TEMP       
C       
C***********************************************************************
C       
C ... CHECK THE VALIDITY OF THE (I,J) ENTRY     
C       
      I = II      
      J = JJ      
      VAL = VALL  
      LEVEL = LEVELL
      NOUT = NOUTT
      IER = 0     
      IF (I.LE.0.OR.I.GT.N) IER = 701 
      IF (J.LE.0.OR.J.GT.N) IER = 701 
      IF (IER.EQ.0) GO TO 20
      IF (LEVEL.GE.0) WRITE (NOUT,10) IER,I,J   
   10 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', 
     *   '    ( ',I10,' , ',I10,' )'/' ',       
     *   '    IMPROPER VALUE FOR I OR J ')      
      GO TO 130   
C       
C ... TRAVERSE THE LINK LIST POINTED TO BY IA(I) UNTIL EITHER       
C ... THE J ENTRY OR THE END OF THE LIST HAS BEEN FOUND.  
C       
   20 NPL1 = N+1  
      LINK = IA(I)
C       
C ...... SPECIAL CASE FOR THE FIRST ENTRY IN THE ROW      
C       
      IF (LINK.GT.0) GO TO 30 
      NEXT = IA(NPL1)       
      IF (NEXT.LT.1) GO TO 110
C       
      IA(I) = NEXT
      JA(NEXT) = J
      A(NEXT) = VAL 
      IWORK(NEXT) = -I      
      IA(NPL1) = NEXT-1     
      GO TO 130   
C       
C ... FOLLOW THE LINK LIST UNTIL J OR THE END OF THE LIST IS FOUND  
C       
   30 IF (JA(LINK).EQ.J) GO TO 40     
      IF (IWORK(LINK).LE.0) GO TO 100 
      LINK = IWORK(LINK)    
      GO TO 30    
C       
C:      
C ... ENTRY (I,J) ALREADY HAS BEEN SET.  RESET VALUE DEPENDING ON MODE
C       
   40 IER = 700   
      IF (MODE.GE.0) GO TO 60 
      IF (LEVEL.GE.1) WRITE (NOUT,50) IER,I,J,A(LINK)     
   50 FORMAT ('0','*** W A R N I N G ************'/'0',   
     *   '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', 
     *   '    ( ',I10,' , ',I10,' )'/' ',       
     *   '    ENTRY ALREADY SET AND IS LEFT AS ',D15.8)   
      GO TO 130   
   60 IF (MODE.GE.1) GO TO 80 
      IF (LEVEL.GE.1) WRITE (NOUT,70) IER,I,J,A(LINK),VAL 
   70 FORMAT ('0','*** W A R N I N G ************'/'0',   
     *   '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', 
     *   '    ( ',I10,' , ',I10,' )'/' ',       
     *   '    ENTRY ALREADY SET - CURRENT VALUE OF',D15.8/' ',      
     *   '                                RESET TO',D15.8)
      A(LINK) = VAL 
      GO TO 130   
   80 TEMP = A(LINK)+VAL    
      IF (LEVEL.GE.1) WRITE (NOUT,90) IER,I,J,A(LINK),TEMP
   90 FORMAT ('0','*** W A R N I N G ************'/'0',   
     *   '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', 
     *   '    ( ',I10,' , ',I10,' )'/' ',       
     *   '    ENTRY ALREADY SET - CURRENT VALUE OF',D15.8/' ',      
     *   '                                RESET TO',D15.8)
      A(LINK) = TEMP
      GO TO 130   
C       
C ... ENTRY (I,J) HAS NOT BEEN SET.  ENTER IT INTO THE LINKED LIST  
C       
  100 NEXT = IA(NPL1)       
      IF (NEXT.LT.1) GO TO 110
C       
      IWORK(LINK) = NEXT    
      JA(NEXT) = J
      A(NEXT) = VAL 
      IWORK(NEXT) = -I      
      IA(NPL1) = NEXT-1     
      GO TO 130   
C       
C***********************************************************************
C       
C ... ERROR TRAP FOR NO ROOM REMAINING
C       
  110 IER = 702   
      IF (LEVEL.GE.0) WRITE (NOUT,120) IER      
  120 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', 
     *   '    NZ TOO SMALL - NO ROOM FOR NEW ENTRY')      
C       
  130 CONTINUE    
      IERR = IER  
      RETURN      
      END 
c=======================
      SUBROUTINE SBEND (N,NZ,IA,JA,A,IWORK)     
C       
C***********************************************************************
C       
C     SBEND IS THE THIRD OF A SUITE OF SUBROUTINES TO AID THE       
C     USER TO CONSTRUCT THE  IA, JA, A DATA STRUCTURE USED IN       
C     ITPACK.     
C       
C     SBEND RESTRUCTURES THE LINKED LIST DATA STRUCTURE BUILT BY    
C     SBINI AND SBSIJ INTO THE FINAL DATA STRUCTURE REQUIRE BY      
C     ITPACK.  THE RESTRUCTURING CAN TAKE PLACE IN THE MINIMUM      
C     AMOUNT OF MEMORY REQUIRED TO HOLD THE NONZERO STRUCTURE OF    
C     THE SPARSE MATRIX BUT WILL RUN QUICKER IF MORE STORAGE
C     IS ALLOWED. 
C       
C     SBEND IS BASED ON SUBROUTINE BUILD OF THE SPARSE MATRIX       
C     PACKAGE SPARSPAK DEVELOPED BY ALAN GEORGE AND JOSEPH LUI      
C     OF THE UNIVERSITY OF WATERLOO, WATERLOO, ONTARIO.   
C       
C ... PARAMETERS  
C       
C ...... INPUT    
C       
C     N       THE ORDER OF THE LINEAR SYSTEM    
C       
C     NZ      THE LENGTH OF THE ARRAYS JA, IWORK, AND A.  
C       
C ...... INPUT/OUTPUT       
C       
C     IA      INTEGER ARRAY OF LENGTH N+1.  THE FIRST N ENTRIES     
C             POINT TO THE BEGINNING OF THE LINKED LIST FOR EACH    
C             ROW.  IA(N+1)-1 IS THE TOP OF THE LINKED LISTS
C             CONTAINED IN JA, IWORK, AND A.  ON OUTPUT IA WILL     
C             POINT TO THE FIRST ENTRY OF EACH ROW IN THE FINAL     
C             DATA STRUCTURE. 
C       
C     JA      INTEGER ARRAY OF LENGTH NZ.  ON INPUT JA STORES THE   
C             COLUMN NUMBERS OF THE NONZERO ENTRIES AS INDICATED    
C             BY THE LINKED LISTS.  ON OUTPUT JA STORES THE 
C             COLUMN NUMBERS IN ROW ORDERED FORM. 
C       
C     A       D.P. ARRAY OF LENGTH NZ.  ON INPUT A STORES THE       
C             VALUE OF THE NOZERO ENTRIES AS INDICATED BY THE       
C             LINKED LISTS.  ON OUTPUT A STORES THE VALUES IN       
C             ROW ORDERED FORM.       
C       
C     IWORK    INTEGER ARRAY OF LENGTH NZ.  ON INPUT IWORK STORES THE 
C             THE LINKS OF THE LINKED LISTS.  ON OUTPUT IT IS       
C             DESTROYED.    
C       
C***********************************************************************
C       
      INTEGER N,NZ,IA(*),JA(NZ),IWORK(NZ)       
      DOUBLE PRECISION A(NZ)
C       
      INTEGER MAXTOP,NEXT,TOP,IDEG,NULINK,JAJ,HLINK,OHLINK,L,I,LINK,
     *   MHLINK   
      DOUBLE PRECISION VAL  
C       
C***********************************************************************
C       
C ... INITIALIZATION
C       
C ...... THE VARIABLES NEXT AND TOP RESPECTIVELY POINT TO THE       
C        NEXT AVAILABLE ENTRY FOR THE FINAL DATA STRUCTURE AND      
C        THE TOP OF THE REMAINDER OF THE LINKED LISTS.    
C       
      NEXT = 1    
      TOP = IA(N+1)+1       
      MAXTOP = NZ-IA(N+1)+1 
C       
C***********************************************************************
C       
C ... CONVERT EACH ROW INTO FINAL FORM
C       
      DO 90 I = 1,N 
         IDEG = 0 
         NULINK = IA(I)     
C       
C ... LOOP OVER EACH NODE IN THE LINKED LIST OF ROW I     
C       
   10    LINK = NULINK      
         IF (LINK.LE.0) GO TO 80      
         NULINK = IWORK(LINK) 
         JAJ = JA(LINK)     
         VAL = A(LINK)      
C       
C ... CHECK TO SEE IF A COLLISION BETWEEN THE LINKED LISTS
C     AND THE FINAL FORM HAS OCCURRED.
C       
         IF (NEXT.GE.TOP.AND.LINK.NE.TOP) GO TO 20
C       
C ... COLLISION HAS NOT OCCURRED.  FREE THE SPACE FOR THE TRIPLE    
C     (JA(LINK), A(LINK), IWORK(LINK))
C       
         JA(LINK) = 0       
         A(LINK) = 0.0D0    
         IWORK(LINK) = 0    
C       
C ... SPECIAL CASE TO MOVE  TOP  DOWN IF LINK .EQ. TOP    
C       
         IF (LINK.EQ.TOP) GO TO 60    
         GO TO 70 
C       
C***********************************************************************
C       
C ... COLLISION HAS OCCURRED.  CLEAR OFF SOME SPACE FOR THE CURRENT 
C     ENTRY BY MOVING THE TRIPLE ( JA(TOP),A(TOP),IWORK(TOP) )      
C     DOWNWARDS TO THE FREED TRIPLE ( JA(LINK),A(LINK),IWORK(LINK) ). 
C     THEN ADJUST THE LINK FIELDS.    
C       
C ...... PATCH UP THE LINKED LIST FOR THE CURRENT ROW I.  THEN      
C        TRAVERSE THE LINKED LIST CONTAINING TOP UNTIL THE POINTER  
C        POINTER BACK TO IA IS FOUND. 
C       
   20    IA(I) = LINK       
         HLINK = TOP
C       
   30    HLINK = IWORK(HLINK) 
         IF (HLINK.GT.0) GO TO 30     
C       
C ...... NOW FOLLOW THE LINKED LIST BACK TO TOP KEEPING TRACK       
C        OF THE OLD LINK.   
C       
C ......... SPECIAL CASE IF IA(-HLINK) = TOP    
C       
         MHLINK = -HLINK    
         IF (IA(MHLINK).NE.TOP) GO TO 40
C       
         IWORK(LINK) = IWORK(TOP)     
         JA(LINK) = JA(TOP) 
         A(LINK) = A(TOP)   
         IA(MHLINK) = LINK  
         IF (NULINK.EQ.TOP) NULINK = LINK       
         GO TO 60 
C       
C ......... USUAL CASE.     
C       
   40    HLINK = IA(MHLINK) 
   50    OHLINK = HLINK     
         HLINK = IWORK(OHLINK)
         IF (HLINK.NE.TOP) GO TO 50   
C       
         IWORK(LINK) = IWORK(TOP)     
         JA(LINK) = JA(TOP) 
         A(LINK) = A(TOP)   
         IF (OHLINK.NE.LINK) IWORK(OHLINK) = LINK 
         IF (NULINK.EQ.TOP) NULINK = LINK       
C       
C ... COLLAPSE TOP OF LINK LIST BY AS MUCH AS POSSIBLE    
C       
   60    TOP = TOP+1
         IF (TOP.GE.MAXTOP) GO TO 70  
         IF (IWORK(TOP).NE.0) GO TO 70
         GO TO 60 
C       
C***********************************************************************
C       
C ... PUT THE CURRENT TRIPLE INTO THE FINAL DATA STRUCTURE
C       
   70    JA(NEXT) = JAJ     
         A(NEXT) = VAL      
         NEXT = NEXT+1      
         IDEG = IDEG+1      
         GO TO 10 
C       
C ... FINAL STRUCTURE FOR ROW I IS COMPLETE.  LINKED LIST IS
C     DESTROYED AND WILL BE RECAPTURED AS NECESSARY BY THE
C     LOOP ON LABEL 60      
C       
   80    IA(I) = IDEG       
C       
   90 CONTINUE    
C       
C***********************************************************************
C       
C ... FINALIZE THE DATA STRUCTURE BY BUILDING THE FINAL VERSION OF  
C     IA. 
C       
      L = IA(1)+1 
      IA(1) = 1   
      DO 100 I = 1,N
         IDEG = IA(I+1)     
         IA(I+1) = L
         L = L+IDEG 
  100 CONTINUE    
C       
C ... FINAL IA, JA, A DATA STRUCTURE BUILT.     
C       
      RETURN      
      END 
c==================================
      SUBROUTINE DFAULT (IPARM,RPARM) 
C       
C ... THIS SUBROUTINE SETS THE DEFAULT VALUES OF IPARM AND RPARM.   
C       
C ... PARAMETER LIST:       
C       
C          IPARM  
C           AND   
C          RPARM  ARRAYS SPECIFYING OPTIONS AND TOLERANCES
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IPARM(12)     
      DOUBLE PRECISION RPARM(12)      
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
C     DRELPR  - COMPUTER PRECISION (APPROX.)    
C     IF INSTALLER OF PACKAGE DOES NOT KNOW DRELPR VALUE, 
C     AN APPROXIMATE VALUE CAN BE DETERMINED FROM A SIMPLE
C     FORTRAN PROGRAM SUCH AS 
C       
C     DOUBLE PRECISION DRELPR, TEMP   
C     DRELPR = 1.0D0
C   2 DRELPR = 0.5D0*DRELPR 
C     TEMP = DRELPR + 1.0D0 
C     IF(TEMP .GT. 1.0D0)  GO TO 2    
C     WRITE(6,3) DRELPR     
C   3 FORMAT(5X,D15.8)      
C     STOP
C     END 
C       
C     SOME VALUES ARE:      
C       
C     DRELPR = 1.26D-29  FOR CDC CYBER 170/750  (APPROX.) 2**-96    
C            = 2.22D-16  FOR DEC 10             (APPROX.) 2**-52    
C            = 7.11D-15  FOR VAX 11/780         (APPROX.) 2**-47    
C            = 1.14D-13  FOR IBM 370/158        (APPROX.) 2**-43    
C       
C             *** SHOULD BE CHANGED FOR OTHER MACHINES ***
C       
C     TO FACILITATE CONVERGENCE, RPARM(1) SHOULD BE SET TO
C          500.*DRELPR OR LARGER      
C       
      DRELPR = 7.11D-15
C       
      IPARM(1) = 100
      IPARM(2) = 0
      IPARM(3) = 0
      IPARM(4) = 6
      IPARM(5) = 0
      IPARM(6) = 1
      IPARM(7) = 1
      IPARM(8) = 0
      IPARM(9) = -1 
      IPARM(10) = 0 
      IPARM(11) = 0 
      IPARM(12) = 0 
C       
      RPARM(1) = 0.5D-5     
      RPARM(2) = 0.D0       
      RPARM(3) = 0.D0       
      RPARM(4) = .75D0      
      RPARM(5) = 1.D0       
      RPARM(6) = 0.D0       
      RPARM(7) = .25D0      
      RPARM(8) = 1.D2*DRELPR
      RPARM(9) = 0.D0       
      RPARM(10) = 0.D0      
      RPARM(11) = 0.D0      
      RPARM(12) = 0.D0      
C       
      RETURN      
      END 
c===============================================
      FUNCTION TIMER (TIMDMY) 
C
C ... TIMER IS A ROUTINE TO RETURN THE EXECUTION TIME IN
C ... SECONDS.
C
C ... PARAMETERS -- 
C
C          TIMDMY   DUMMY ARGUMENT
C
C
C *********************************************
C **                                         **
C **   THIS ROUTINE IS NOT PORTABLE.         **
C **                                         **
C *********************************************
C
      REAL TIMDMY
C
C ... CRAY Y-MP.
C
C     TIMER = SECOND ()
C
C ... UNIX ETIME FACILITY.
C
c     REAL ETIME
c     EXTERNAL ETIME
      DIMENSION TARRAY(2)
      REAL TIMER
      TOTAL = ETIME (TARRAY)
      TIMER = TOTAL
C
C ... IBM RISC SYSTEM/6000.
C
C     TIMER = FLOAT(MCLOCK())/100.0
C
      RETURN
      END 
c====================================
      SUBROUTINE ECHOUT (IPARM,RPARM,IMTHD)     
C       
C     THIS ROUTINE INITIALIZES THE ITPACK COMMON BLOCKS FROM THE    
C     INFORMATION CONTAINED IN IPARM AND RPARM. 
C       
C ... PARAMETER LIST:       
C       
C          IPARM  
C           AND   
C          RPARM  ARRAYS OF PARAMETERS SPECIFYING OPTIONS AND       
C                    TOLERANCES       
C          IMTHD  INDICATOR OF METHOD 
C                    IMTHD = 1,  JCG  
C                    IMTHD = 2,  JSI  
C                    IMTHD = 3,  SOR  
C                    IMTHD = 4,  SSORCG 
C                    IMTHD = 5,  SSORSI 
C                    IMTHD = 6,  RSCG 
C                    IMTHD = 7,  RSSI 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IPARM(12),IMTHD 
      DOUBLE PRECISION RPARM(12)      
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
C ... INITIALIZE ITPACK COMMON
C       
      ZETA = RPARM(1)       
      CME = RPARM(2)
      SME = RPARM(3)
      FF = RPARM(4) 
      OMEGA = RPARM(5)      
      SPECR = RPARM(6)      
      BETAB = RPARM(7)      
      ITMAX = IPARM(1)      
      LEVEL = IPARM(2)      
      ISYM = IPARM(5)       
C       
      ADAPT = .FALSE.       
      PARTAD = .FALSE.      
      BETADT = .FALSE.      
      IF (IPARM(6).EQ.1.OR.IPARM(6).EQ.3) ADAPT = .TRUE.  
      IF (IPARM(6).EQ.1) BETADT = .TRUE.
      IF (IPARM(6).EQ.2) PARTAD = .TRUE.
C       
      CASEII = .FALSE.      
      IF (IPARM(7).EQ.2) CASEII = .TRUE.
      IF (CASEII) SME = -CME
      IF (.NOT.CASEII.AND.SME.EQ.0.D0) SME = -1.D0
      SPR = SME   
C       
C ... SET REST OF COMMON VARIABLES TO ZERO      
C       
      IN = 0      
      IS = 0      
      HALT = .FALSE.
      BDELNM = 0.D0 
      DELNNM = 0.D0 
      DELSNM = 0.D0 
      GAMMA = 0.D0
      QA = 0.D0   
      QT = 0.D0   
      RHO = 0.D0  
      RRR = 0.D0  
      SIGE = 0.D0 
      STPTST = 0.D0 
      UDNM = 0.D0 
      IF (LEVEL.LE.2) RETURN
C       
C ... THIS SECTION OF ECHOUT ECHOES THE INPUT VALUES FOR THE INITIAL
C     ITERATIVE PARAMETERS  
C       
      WRITE (NOUT,10) ISYM,ITMAX,ZETA,ADAPT,CASEII
   10 FORMAT (///30X,'INITIAL ITERATIVE PARAMETERS',3X,   
     *   'RELEVANT SWITCHES'/35X,'ISYM   =',I15,8X,'IPARM(5)'/35X,  
     *   'ITMAX  =',I15,8X,'IPARM(1)'/35X,'ZETA   =',D15.8,8X,'RPARM(1)'
     *   /35X,'ADAPT  =',L15,8X,'IPARM(6)'/35X,'CASEII =',L15,8X,   
     *   'IPARM(7)')
      GO TO (80,20,100,60,40,80,20), IMTHD      
C       
C ... JSI, RSSI   
C       
   20 WRITE (NOUT,30) FF,CME,SME      
   30 FORMAT (35X,'FF     =',D15.8,8X,'RPARM(4)'/35X,'CME    =',D15.8,8X
     *   ,'RPARM(2)'/35X,'SME    =',D15.8,8X,'RPARM(3)'///) 
      RETURN      
C       
C ... SSORSI      
C       
   40 WRITE (NOUT,50) PARTAD,FF,CME,OMEGA,SPECR,BETAB,BETADT
   50 FORMAT (35X,'PARTAD =',L15,8X,'IPARM(6)'/35X,'FF     =',D15.8,8X, 
     *   'RPARM(4)'/35X,'CME    =',D15.8,8X,'RPARM(2)'/35X,'OMEGA  =',
     *   D15.8,8X,'RPARM(5)'/35X,'SPECR  =',D15.8,8X,'RPARM(6)'/35X,
     *   'BETAB  =',D15.8,8X,'RPARM(7)'/35X,'BETADT =',L15,8X,'IPARM(6)'
     *   ///)     
      RETURN      
C       
C ... SSORCG      
C       
   60 WRITE (NOUT,70) PARTAD,CME,OMEGA,SPECR,BETAB,BETADT 
   70 FORMAT (35X,'PARTAD =',L15,8X,'IPARM(6)'/35X,'CME    =',D15.8,8X, 
     *   'RPARM(2)'/35X,'OMEGA  =',D15.8,8X,'RPARM(5)'/35X,'SPECR  =',
     *   D15.8,8X,'RPARM(6)'/35X,'BETAB  =',D15.8,8X,'RPARM(7)'/35X,
     *   'BETADT =',L15,8X,'IPARM(6)'///)       
      RETURN      
C       
C ... JCG, RSCG   
C       
   80 IF (ADAPT) RETURN     
      WRITE (NOUT,90) CME   
   90 FORMAT (35X,'CME    =',D15.8,8X,'RPARM(2)'///)      
C       
  100 CONTINUE    
      RETURN      
      END 
c===================================================
      SUBROUTINE ECHALL (NN,IA,JA,A,RHS,IPARM,RPARM,ICALL)
C       
C ... THIS ROUTINE INITIALIZES THE ITPACK COMMON BLOCKS FROM THE    
C ... INFORMATION CONTAINED IN IPARM AND RPARM. ECHALL ALSO PRINTS THE
C ... VALUES OF ALL THE PARAMETERS IN IPARM AND RPARM.    
C       
C ... PARAMETER LIST:       
C       
C          IPARM  
C           AND   
C          RPARM  ARRAYS OF PARAMETERS SPECIFYING OPTIONS AND       
C                    TOLERANCES       
C          ICALL  INDICATOR OF WHICH PARAMETERS ARE BEING PRINTED   
C                    ICALL = 1,  INITIAL PARAMETERS       
C                    ICALL = 2,  FINAL PARAMETERS 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(*),JA(*),IPARM(12),NN,ICALL    
      DOUBLE PRECISION A(*),RHS(NN),RPARM(12)   
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,N,NP1,NZRO  
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
      IF (ICALL.NE.1) GO TO 100       
      N = NN      
      NP1 = N+1   
      NZRO = IA(NP1)-1      
C       
C ... INITIALIZE ITPACK COMMON
C       
      ZETA = RPARM(1)       
      CME = RPARM(2)
      SME = RPARM(3)
      FF = RPARM(4) 
      OMEGA = RPARM(5)      
      SPECR = RPARM(6)      
      BETAB = RPARM(7)      
      ITMAX = IPARM(1)      
      LEVEL = IPARM(2)      
      ISYM = IPARM(5)       
C       
      ADAPT = .FALSE.       
      PARTAD = .FALSE.      
      BETADT = .FALSE.      
      IF (IPARM(6).EQ.1.OR.IPARM(6).EQ.3) ADAPT = .TRUE.  
      IF (IPARM(6).EQ.1) BETADT = .TRUE.
      IF (IPARM(6).EQ.2) PARTAD = .TRUE.
C       
      CASEII = .FALSE.      
      IF (IPARM(7).EQ.2) CASEII = .TRUE.
      IF (CASEII) SME = -CME
      IF (.NOT.CASEII.AND.SME.EQ.0.D0) SME = -1.D0
      SPR = SME   
C       
C ... SET REST OF COMMON VARIABLES TO ZERO      
C       
      IN = 0      
      IS = 0      
      HALT = .FALSE.
      BDELNM = 0.D0 
      DELNNM = 0.D0 
      DELSNM = 0.D0 
      GAMMA = 0.D0
      QA = 0.D0   
      QT = 0.D0   
      RHO = 0.D0  
      RRR = 0.D0  
      SIGE = 0.D0 
      STPTST = 0.D0 
      UDNM = 0.D0 
C       
      IF (LEVEL.LE.4) GO TO 80
C       
C     THIS SECTION OF ECHALL CAUSES PRINTING OF THE LINEAR SYSTEM AND 
C     THE ITERATIVE PARAMETERS
C       
      WRITE (NOUT,10)       
   10 FORMAT (///30X,'THE LINEAR SYSTEM IS AS FOLLOWS')   
      WRITE (NOUT,20)       
   20 FORMAT (/2X,'IA ARRAY') 
      WRITE (NOUT,30) (IA(I),I=1,NP1) 
   30 FORMAT (2X,10(2X,I8)) 
      WRITE (NOUT,40)       
   40 FORMAT (/2X,'JA ARRAY') 
      WRITE (NOUT,30) (JA(I),I=1,NZRO)
      WRITE (NOUT,50)       
   50 FORMAT (/2X,' A ARRAY') 
      WRITE (NOUT,60) (A(I),I=1,NZRO) 
   60 FORMAT (2X,5(2X,D20.13))
      WRITE (NOUT,70)       
   70 FORMAT (/2X,'RHS ARRAY')
      WRITE (NOUT,60) (RHS(I),I=1,N)  
   80 WRITE (NOUT,90)       
   90 FORMAT (///30X,'INITIAL ITERATIVE PARAMETERS')      
      GO TO 120   
  100 WRITE (NOUT,110)      
  110 FORMAT (///30X,'FINAL ITERATIVE PARAMETERS')
  120 WRITE (NOUT,130) IPARM(1),LEVEL,IPARM(3),NOUT,ISYM,IPARM(6)   
  130 FORMAT (35X,'IPARM(1)  =',I15,4X,'(ITMAX)'/35X,'IPARM(2)  =',I15, 
     *   4X,'(LEVEL) '/35X,'IPARM(3)  =',I15,4X,'(IRESET)'/35X,     
     *   'IPARM(4)  =',I15,4X,'(NOUT)  '/35X,'IPARM(5)  =',I15,4X,  
     *   '(ISYM)  '/35X,'IPARM(6)  =',I15,4X,'(IADAPT)')  
      WRITE (NOUT,140) IPARM(7),IPARM(8),IPARM(9),IPARM(10),IPARM(11),
     *   IPARM(12)
  140 FORMAT (35X,'IPARM(7)  =',I15,4X,'(ICASE)'/35X,'IPARM(8)  =',I15, 
     *   4X,'(NWKSP)'/35X,'IPARM(9)  =',I15,4X,'(NB)    '/35X,      
     *   'IPARM(10) =',I15,4X,'(IREMOVE)'/35X,'IPARM(11) =',I15,4X, 
     *   '(ITIME)'/35X,'IPARM(12) =',I15,4X,'(IDGTS)')    
      WRITE (NOUT,150) ZETA,CME,SME,FF,OMEGA,SPECR
  150 FORMAT (35X,'RPARM(1)  =',D15.8,4X,'(ZETA)  '/35X,'RPARM(2)  =',
     *   D15.8,4X,'(CME)   '/35X,'RPARM(3)  =',D15.8,4X,'(SME)   '/35X, 
     *   'RPARM(4)  =',D15.8,4X,'(FF)    '/35X,'RPARM(5)  =',D15.8,4X,
     *   '(OMEGA) '/35X,'RPARM(6)  =',D15.8,4X,'(SPECR) ')
      WRITE (NOUT,160) BETAB,RPARM(8),RPARM(9),RPARM(10),RPARM(11), 
     *   RPARM(12)
  160 FORMAT (35X,'RPARM(7)  =',D15.8,4X,'(BETAB) '/35X,'RPARM(8)  =',
     *   D15.8,4X,'(TOL)'/35X,'RPARM(9)  =',D15.8,4X,'(TIME1)'/35X, 
     *   'RPARM(10) =',D15.8,4X,'(TIME2)'/35X,'RPARM(11) =',D15.8,4X, 
     *   '(DIGIT1)'/35X,'RPARM(12) =',D15.8,4X,'(DIGIT2)')
C       
      RETURN      
      END 
c==================================
      SUBROUTINE IVFILL (N,IV,IVAL)   
C       
C     FILLS AN INTEGER VECTOR, IV, WITH AN INTEGER VALUE, IVAL.     
C       
C ... PARAMETER LIST:       
C       
C          N      INTEGER LENGTH OF VECTOR IV   
C          IV     INTEGER VECTOR      
C          IVAL   INTEGER CONSTANT THAT FILLS FIRST N LOCATIONS OF IV 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N,IVAL,IV(N)  
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,M,MP1       
C       
      IF (N.LE.0) RETURN    
C       
C     CLEAN UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 10  
C       
      M = MOD(N,10) 
      IF (M.EQ.0) GO TO 20  
      DO 10 I = 1,M 
         IV(I) = IVAL       
   10 CONTINUE    
      IF (N.LT.10) RETURN   
C       
   20 MP1 = M+1   
      DO 30 I = MP1,N,10    
         IV(I) = IVAL       
         IV(I+1) = IVAL     
         IV(I+2) = IVAL     
         IV(I+3) = IVAL     
         IV(I+4) = IVAL     
         IV(I+5) = IVAL     
         IV(I+6) = IVAL     
         IV(I+7) = IVAL     
         IV(I+8) = IVAL     
         IV(I+9) = IVAL     
   30 CONTINUE    
C       
      RETURN      
      END 
c========================================
      SUBROUTINE VFILL (N,V,VAL)      
C       
C     FILLS A VECTOR, V, WITH A CONSTANT VALUE, VAL.      
C       
C ... PARAMETER LIST:       
C       
C          N      INTEGER LENGTH OF VECTOR V    
C          V      D.P. VECTOR 
C          VAL    D.P. CONSTANT THAT FILLS FIRST N LOCATIONS OF V   
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N   
      DOUBLE PRECISION V(N),VAL       
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,M,MP1       
C       
      IF (N.LE.0) RETURN    
C       
C     CLEAN UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 10  
C       
      M = MOD(N,10) 
      IF (M.EQ.0) GO TO 20  
      DO 10 I = 1,M 
         V(I) = VAL 
   10 CONTINUE    
      IF (N.LT.10) RETURN   
C       
   20 MP1 = M+1   
      DO 30 I = MP1,N,10    
         V(I) = VAL 
         V(I+1) = VAL       
         V(I+2) = VAL       
         V(I+3) = VAL       
         V(I+4) = VAL       
         V(I+5) = VAL       
         V(I+6) = VAL       
         V(I+7) = VAL       
         V(I+8) = VAL       
         V(I+9) = VAL       
   30 CONTINUE    
C       
      RETURN      
      END 
c===============================
      SUBROUTINE SBELM (NN,IA,JA,A,RHS,IW,RW,TOL,ISYM,LEVEL,NOUT,IER) 
C       
C ... SBELM IS DESIGNED TO REMOVE ROWS AND COLUMNS OF THE MATRIX    
C ... WHERE DABS(A(I,J))/A(I,I) .LE. TOL FOR J = 1 TO N AND A(I,I)  
C ... .GT. 0. THIS IS TO TAKE CARE OF MATRICES ARISING    
C ... FROM FINITE ELEMENT DISCRETIZATIONS OF PDE^S WITH DIRICHLET   
C ... BOUNDARY CONDITIONS.  ANY SUCH ROWS AND CORRESPONDING COLUMNS 
C ... ARE THEN SET TO THE IDENTITY AFTER CORRECTING RHS.  
C       
C ... PARAMETER LIST:       
C       
C          N      DIMENSION OF MATRIX (= NN)    
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
C          IW,RW  WORK ARRAYS OF LENGTH N       
C          TOL    TOLERANCE FACTOR    
C          ISYM   FLAG FOR TYPE OF STORAGE FOR SYSTEM     
C                 (0: SYMMETRIC, 1:NONSYMMETRIC)
C          LEVEL  PRINTING SWITCH FOR ERROR CONDITION     
C          NOUT OUTPUT TAPE NUMBER    
C          IER    ERROR FLAG: NONZERO VALUE ON RETURN MEANS 
C                    101 : DIAGONAL ENTRY NOT POSITIVE    
C                    102 : THERE IS NO DIAGONAL ENTRY IN ROW
C       
C********************************************************************** 
C       
C     UPDATE.  SBELM HAS BEEN REWRITTEN TO SPEED UP THE LOCATION OF 
C              OF ROWS WHICH ARE TO BE ELIMINATED.  THIS IS DONE BY 
C              FIRST STORING THE LARGEST ELEMENT OF EACH ROW IN     
C              THE ARRAY RW.  THE DIAGONAL ENTRY IS THEN COMPARED   
C              WITH THE CORRESPONDING ELEMENT IN RW.  IF IT IS      
C              DECIDED TO ELIMINATE THE ROW THEN IT IS MARKED FOR   
C              ELIMINATION. 
C       
C              WHEN A ROW IS TO BE ELIMINATED ITS DIAGONAL ENTRY    
C              IS STORED IN  RW  AND  IW IS MARKED BY A NONZERO     
C              (WHICH IS THIS ROW NUMBER)       
C       
C              ROWS WHICH HAVE ONLY DIAGONAL ENTRIES ARE NOT
C              ALTERED.     
C       
C*********************************************************************
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER NN,IA(*),JA(*),IW(NN),ISYM,LEVEL,NOUT,IER   
      DOUBLE PRECISION A(*),RHS(NN),RW(NN),TOL  
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IBGN,ICNT,IEND,JJ,JJDI,KK,N       
      DOUBLE PRECISION DI   
C       
      N = NN      
C       
C        IF (N .GE. 1) GO TO 10       
C           IER = 100       
C           RETURN
C 10     CONTINUE 
C       
C ... STORE THE LARGEST (DABSOLUTE VALUE) OFF DIAGONAL ENTRY FOR    
C ... ROW II IN RW(II).     
C       
      IER = 0     
      ICNT = 0    
      DO 10 II = 1,N
         RW(II) = 0.0D0     
         IW(II) = 0 
   10 CONTINUE    
      DO 20 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 140  
         DO 20 JJ = IBGN,IEND 
            KK = JA(JJ)     
            IF (KK.EQ.II) GO TO 20    
            RW(II) = DMAX1(RW(II),DABS(A(JJ)))  
            IF (ISYM.NE.0) GO TO 20   
            RW(KK) = DMAX1(RW(KK),DABS(A(JJ)))  
   20 CONTINUE    
C       
C ... FOR II = 1 TO N FIND THE DIAGONAL ENTRY IN ROW II   
C       
      DO 80 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         DO 40 JJ = IBGN,IEND 
            IF (JA(JJ).NE.II) GO TO 40
            DI = A(JJ)      
            JJDI = JJ       
            IF (DI.GT.0.D0) GO TO 50  
            IER = 101       
            IF (LEVEL.GE.0) WRITE (NOUT,30) II,DI 
   30       FORMAT ('0','*** F A T A L     E R R O R ************'/'0', 
     *         '    IN ITPACK ROUTINE SBELM   '/' ',      
     *         '    DIAGONAL ELEMENT',I10,' NOT POSITIVE  '/' ',    
     *         '    CURRENT VALUE = ',D15.8)    
            RETURN
   40    CONTINUE 
         GO TO 140
   50    CONTINUE 
C       
C ... CHECK THE SIZE OF THE LARGEST OFF DIAGONAL ELEMENT  
C ... ( STORED IN RW(II) ) AGAINST THE DIAGONAL ELEMENT DII.
C       
         IF (RW(II).NE.0.0D0) GO TO 60
         IF (1.0D0/DI.LE.TOL) GO TO 70
         GO TO 80 
   60    IF (RW(II)/DI.GT.TOL) GO TO 80 
C       
C ... THE OFF DIAGONAL ELEMENTS ARE SMALL COMPARED TO THE DIAGONAL  
C ... THEREFORE MARK IT FOR ELIMINATION AND PERFORM INITIAL 
C ... PROCESSING  
C       
   70    ICNT = ICNT+1      
         IW(II) = II
         RW(II) = DI
         A(JJDI) = 1.0D0    
         RHS(II) = RHS(II)/DI 
C       
   80 CONTINUE    
C       
C ... ELIMINATE THE ROWS AND COLUMNS INDICATED BY THE NONZERO       
C ... ENTRIES IN IW.  THERE ARE ICNT OF THEM    
C       
      IF (ICNT.EQ.0) GO TO 130
C       
C ... THE ELIMINATION IS AS FOLLOWS:  
C       
C     FOR II = 1 TO N DO    
C        IF ( IW(II) .NE. 0 ) THEN    
C           SET DIAGONAL VALUE TO 1.0  ( ALREADY DONE )   
C           SET RHS(II) = RHS(II) / RW(II)   ( ALREADY DONE )       
C           FIND NONZERO OFFDIAGONAL ENTRIES  KK
C           IF ( IW(KK) .EQ. 0 ) FIX UP RHS(KK)  WHEN USING SYMMETRIC ST
C           SET A(II,KK) = 0.0
C        ELSE ( I.E.  IW(II) .EQ. 0  )
C           FIND NONZERO OFFDIAGONAL ENTRIES   KK 
C           IF ( IW(KK) .NE. 0 ) FIX UP RHS(II) 
C                                AND SET A(II,KK) = 0.0   
C        END IF   
C     END DO      
C       
      DO 120 II = 1,N       
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IW(II).EQ.0) GO TO 100   
C       
C ... THE II-TH ROW IS TO BE ELIMINATED 
C       
         DO 90 JJ = IBGN,IEND 
            KK = JA(JJ)     
            IF (KK.EQ.II) GO TO 90    
            IF ((IW(KK).EQ.0).AND.(ISYM.EQ.0)) RHS(KK) = RHS(KK)-A(JJ)* 
     *         RHS(II)      
            A(JJ) = 0.0D0   
   90    CONTINUE 
         GO TO 120
C       
C ... THE II-TH ROW IS KEPT.  CHECK THE OFF-DIAGONAL ENTRIES
C       
  100    DO 110 JJ = IBGN,IEND
            KK = JA(JJ)     
            IF (KK.EQ.II.OR.IW(KK).EQ.0) GO TO 110
            RHS(II) = RHS(II)-A(JJ)*RHS(KK)     
            A(JJ) = 0.0D0   
  110    CONTINUE 
C       
  120 CONTINUE    
C       
  130 RETURN      
C       
C ... ERROR TRAPS -- NO DIAGONAL ENTRY IN ROW II (ROW MAY BE EMPTY).
C       
  140 CONTINUE    
      IER = 102   
      IF (LEVEL.GE.0) WRITE (NOUT,150) II       
  150 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE SBELM   '/' ',  
     *   '    NO DIAGONAL ENTRY IN ROW  ',I10)  
C       
      RETURN      
      END 
c=====================================
      SUBROUTINE PRBNDX (NN,NBLACK,IA,JA,P,IP,LEVEL,NOUT,IER)       
C       
C**************************************************************     
C       
C     THIS SUBROUTINE COMPUTES THE RED-BLACK PERMUTATION  
C     VECTORS P ( AND ITS INVERSE IP ) IF POSSIBLE.       
C       
C     THE ALGORITHM IS TO MARK THE FIRST NODE AS RED (ARBITRARY).   
C     ALL OF ITS ADJACENT NODES ARE MARKED BLACK AND PLACED IN      
C     A STACK.  THE REMAINDER OF THE CODE PULLS THE FIRST NODE      
C     OFF THE TOP OF THE STACK AND TRIES TO TYPE ITS ADJACENT NODES.
C     THE TYPING OF THE ADJACENT POINT IS A FIVE WAY CASE STATEMENT 
C     WHICH IS WELL COMMENTED BELOW (SEE DO LOOP 100).    
C       
C     THE ARRAY P IS USED BOTH TO KEEP TRACK OF THE COLOR OF A NODE 
C     (RED NODE IS POSITIVE, BLACK IS NEGATIVE) BUT ALSO THE FATHER 
C     NODE THAT CAUSED THE COLOR MARKING OF THAT POINT.  SINCE      
C     COMPLETE INFORMATION ON THE ADJACENCY STRUCTURE IS HARD TO COME 
C     BY THIS FORMS A LINK TO ENABLE THE COLOR CHANGE OF A PARTIAL  
C     TREE WHEN A RECOVERABLE COLOR CONFLICT OCCURS.      
C       
C     THE ARRAY IP IS USED AS A STACK TO POINT TO THE SET OF NODES  
C     LEFT TO BE TYPED THAT ARE KNOWN TO BE ADJACENT TO THE CURRENT 
C     FATHER NODE.
C       
C*********************************************************************
C       
C     INPUT PARAMETERS      
C       
C        N      NUMBER OF NODES.  (INTEGER, SCALAR) (= NN)
C       
C        IA,JA  ADJACENCY STRUCTURE ARRAYS.  CAN BE EITHER THE      
C               SYMMETRIC OR NONSYMMETRIC FORM.  IT IS ASSUMED      
C               THAT FOR EVERY ROW WHERE ONLY ONE ELEMENT IS
C               STORED THAT ELEMENT CORRESPONDS TO THE DIAGONAL     
C               ENTRY.  THE DIAGONAL DOES NOT HAVE TO BE THE FIRST  
C               ENTRY STORED.  (INTEGER, ARRAYS)
C        LEVEL  SWITCH FOR PRINTING   
C        NOUT OUTPUT TAPE NUMBER      
C       
C     OUTPUT PARAMETERS     
C       
C        NBLACK NUMBER OF BLACK NODES.  NUMBER OF RED NODES IS      
C               N - NBLACK.  (INTEGER, SCALAR)  
C       
C        P, IP  PERMUTATION AND INVERSE PERMUTATION VECTORS.
C               (INTEGER, ARRAYS EACH OF LENGTH N)
C       
C        IER    ERROR FLAG. (INTEGER, SCALAR)   
C       
C               IER = 0, NORMAL RETURN.  INDEXING PERFORMED 
C                        SUCCESSFULLY 
C               IER =201, RED-BLACK INDEXING NOT POSSIBLE.
C       
C******************************************************************** 
C       
      INTEGER NN,NBLACK,IA(*),JA(*),P(NN),IP(NN),IER      
C       
      INTEGER FIRST,NEXT,LAST,I,OLD,YOUNG,IBGN,IEND,J,K,CURTYP,NXTTYP,
     *   TYPE,NRED,N
C       
C-----------------------------------------------------------------------
C       
      N = NN      
      IER = 0     
C       
C        IF ( N .LE. 0 ) GO TO 8000   
C       
      DO 10 I = 1,N 
         P(I) = 0 
         IP(I) = 0
   10 CONTINUE    
C       
C ... HANDLE THE FIRST SET OF POINTS UNTIL SOME ADJACENT POINTS     
C ... ARE FOUND   
C       
      FIRST = 1   
C       
   20 P(FIRST) = FIRST      
      IF (IA(FIRST+1)-IA(FIRST).GT.1) GO TO 40  
C       
C ... SEARCH FOR NEXT ENTRY THAT HAS NOT BEEN MARKED      
C       
      IF (FIRST.EQ.N) GO TO 130       
      IBGN = FIRST+1
      DO 30 I = IBGN,N      
         IF (P(I).NE.0) GO TO 30      
         FIRST = I
         GO TO 20 
   30 CONTINUE    
      GO TO 130   
C       
C ... FIRST SET OF ADJACENT POINTS FOUND
C       
   40 NEXT = 1    
      LAST = 1    
      IP(1) = FIRST 
C       
C ... LOOP OVER LABELED POINTS INDICATED IN THE STACK STORED IN     
C ... THE ARRAY IP
C       
   50 K = IP(NEXT)
      CURTYP = P(K) 
      NXTTYP = -CURTYP      
      IBGN = IA(K)
      IEND = IA(K+1)-1      
      IF (IBGN.GT.IEND) GO TO 110     
      DO 100 I = IBGN,IEND  
         J = JA(I)
         TYPE = P(J)
         IF (J.EQ.K) GO TO 100
C       
C================================================================== 
C       
C     THE FOLLOWING IS A FIVE WAY CASE STATEMENT DEALING WITH THE   
C     LABELING OF THE ADJACENT NODE.  
C       
C ... CASE I.  IF THE ADJACENT NODE HAS ALREADY BEEN LABELED WITH   
C              LABEL EQUAL TO NXTTYP, THEN SKIP TO THE NEXT ADJACENT
C              NODE.
C       
         IF (TYPE.EQ.NXTTYP) GO TO 100
C       
C ... CASE II.  IF THE ADJACENT NODE HAS NOT BEEN LABELED YET LABEL 
C               IT WITH NXTTYP AND ENTER IT IN THE STACK  
C       
         IF (TYPE.NE.0) GO TO 60      
         LAST = LAST+1      
         IP(LAST) = J       
         P(J) = NXTTYP      
         GO TO 100
C       
C ... CASE III.  IF THE ADJACENT NODE HAS ALREADY BEEN LABELED WITH 
C                OPPOSITE COLOR AND THE SAME FATHER SEED, THEN THERE
C                IS AN IRRECOVERABLE COLOR CONFLICT.      
C       
   60    IF (TYPE.EQ.CURTYP) GO TO 160
C       
C ... CASE IV.  IF THE ADJACENT NODE HAS THE RIGHT COLOR AND A DIFFERENT
C               FATHER NODE, THEN CHANGE ALL NODES OF THE YOUNGEST FATHE
C               NODE TO POINT TO THE OLDEST FATHER SEED AND RETAIN THE
C               SAME COLORS.
C       
         IF (TYPE*NXTTYP.LT.1) GO TO 80 
         OLD = MIN0(IABS(TYPE),IABS(NXTTYP))    
         YOUNG = MAX0(IABS(TYPE),IABS(NXTTYP))  
         DO 70 J = YOUNG,N  
            IF (IABS(P(J)).EQ.YOUNG) P(J) = ISIGN(OLD,P(J)) 
   70    CONTINUE 
         CURTYP = P(K)      
         NXTTYP = -CURTYP   
         GO TO 100
C       
C ... CASE V.  IF THE ADJACENT NODE HAS THE WRONG COLOR AND A DIFFERENT 
C              FATHER NODE, THEN CHANGE ALL NODES OF THE YOUNGEST FATHER
C              NODE TO POINT TO THE OLDEST FATHER NODE ALONG WITH   
C              CHANGING THEIR COLORS.  SINCE UNTIL THIS TIME THE    
C              YOUNGEST FATHER NODE TREE HAS BEEN INDEPENDENT NO OTHER
C              COLOR CONFLICTS WILL ARISE FROM THIS CHANGE. 
C       
   80    OLD = MIN0(IABS(TYPE),IABS(NXTTYP))    
         YOUNG = MAX0(IABS(TYPE),IABS(NXTTYP))  
         DO 90 J = YOUNG,N  
            IF (IABS(P(J)).EQ.YOUNG) P(J) = ISIGN(OLD,-P(J))
   90    CONTINUE 
         CURTYP = P(K)      
         NXTTYP = -CURTYP   
C       
C ... END OF CASE STATEMENT 
C       
C================================================================== 
C       
  100 CONTINUE    
C       
C ... ADVANCE TO NEXT NODE IN THE STACK 
C       
  110 NEXT = NEXT+1 
      IF (NEXT.LE.LAST) GO TO 50      
C       
C ... ALL NODES IN THE STACK HAVE BEEN REMOVED  
C       
C ... CHECK FOR NODES NOT LABELED.  IF ANY ARE FOUND      
C ... START THE LABELING PROCESS AGAIN AT THE FIRST       
C ... NODE FOUND THAT IS NOT LABELED. 
C       
      IBGN = FIRST+1
      DO 120 I = IBGN,N     
         IF (P(I).NE.0) GO TO 120     
         FIRST = I
         GO TO 20 
  120 CONTINUE    
C       
C===================================================================
C       
C ... ALL NODES ARE NOW TYPED EITHER RED OR BLACK 
C       
C ... GENERATE PERMUTATION VECTORS    
C       
  130 NRED = 0    
      NBLACK = 0  
      DO 150 I = 1,N
         IF (P(I).LT.0) GO TO 140     
C       
C       RED POINT 
C       
         NRED = NRED+1      
         IP(NRED) = I       
         P(I) = NRED
         GO TO 150
C       
C     BLACK POINT 
C       
  140    NBLACK = NBLACK+1  
         J = N-NBLACK+1     
         IP(J) = I
         P(I) = J 
C       
  150 CONTINUE    
C       
C ... SUCCESSFUL RED-BLACK ORDERING COMPLETED   
C       
      GO TO 180   
C       
C ........ ERROR TRAPS      
C       
C ...... N .LE. 0 
C       
C8000    IER = 200
C        GO TO 9000 
C       
C ...... TYPE CONFLICT      
C       
  160 IER = 201   
      IF (LEVEL.GE.0) WRITE (NOUT,170)
  170 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE PRBNDX  '/' ',  
     *   '    RED-BLACK INDEXING NOT POSSIBLE') 
C       
C ... RETURN      
C       
  180 CONTINUE    
      RETURN      
      END 
c===================================
      SUBROUTINE PERMAT (NN,IA,JA,A,P,NEWIA,ISYM,LEVEL,NOUT,IERR)   
C       
C*********************************************************************
C       
C ... SUBROUTINE PERMAT TAKES THE SPARSE MATRIX REPRESENTATION      
C     OF THE MATRIX STORED IN THE ARRAYS IA, JA, AND A AND
C     PERMUTES BOTH ROWS AND COLUMNS OVERWRITING THE PREVIOUS       
C     STRUCTURE.  
C       
C ... PARAMETER LIST:       
C       
C         N      ORDER OF SYSTEM (= NN) 
C         IA,JA  INTEGER ARRAYS OF THE SPARSE MATRIX REPRESENTATION 
C         A      D.P. ARRAY OF THE SPARSE MATRIX REPRESENTATION     
C         P      PERMUTATION VECTOR   
C         NEWIA  INTEGER WORK VECTOR OF LENGTH N
C         ISYM   SYMMETRIC/NONSYMMETRIC STORAGE SWITCH    
C         LEVEL  SWITCH CONTROLLING LEVEL OF OUTPUT       
C         NOUT OUTPUT UNIT NUMBER     
C         IER    OUTPUT ERROR FLAG (= IERR)     
C       
C                   IER =   0  NORMAL RETURN    
C                   IER = 301  NO ENTRY IN ITH ROW OF ORIGINAL      
C                              MATRIX. IF LEVEL IS GREATER THAN     
C                              0, I WILL BE PRINTED       
C                   IER = 302  THERE IS NO ENTRY IN THE ITH ROW     
C                              OF THE PERMUTED MATRIX     
C                   IER = 303  ERROR RETURN FROM QSORT IN 
C                              SORTING THE ITH ROW OF THE 
C                              PERMUTED MATRIX  
C ... IT IS ASSUMED THAT THE I-TH ENTRY OF THE PERMUTATION VECTOR   
C     P INDICATES THE ROW THE I-TH ROW GETS MAPPED INTO.  (I.E.     
C     IF ( P(I) = J ) ROW I GETS MAPPED INTO ROW J.)      
C       
C ... THE ARRAY NEWIA IS AN INTEGER WORK VECTOR OF LENGTH N WHICH   
C     KEEPS TRACK OF WHERE THE ROWS BEGIN IN THE PERMUTED STRUCTURE.
C       
C ... PERMAT IS CAPABLE OF PERMUTING BOTH THE SYMMETRIC AND NON-    
C     SYMMETRIC FORM OF IA, JA, AND A.  IF ( ISYM .EQ. 0 ) SYMMETRIC
C     FORM IS ASSUMED.      
C       
C ... TWO EXTERNAL MODULES ARE USED BY PERMAT.  THE FIRST IS INTEGER
C     FUNCTION BISRCH WHICH USES A BISECTION SEARCH ( ORDER LOG-BASE-2
C     OF N+1 ) THROUGH THE ARRAY IA TO FIND THE ROW INDEX OF AN ARBI- 
C     TRARY ENTRY EXTRACTED FROM THE ARRAY JA. THE SECOND IS SUBROUTINE 
C     QSORT WHICH PERFORMS A QUICK SORT TO PLACE THE ENTRIES IN     
C     THE PERMUTED ROWS IN COLUMN ORDER.
C       
C*********************************************************************
C       
      INTEGER NN,IA(*),JA(*),P(NN),NEWIA(NN),ISYM,IERR    
      DOUBLE PRECISION A(*) 
C       
C ... INTERNAL VARIABLES    
C       
      INTEGER BISRCH,I,IBGN,IEND,IP,IPP,J,JAJ,JP,IER,K,N,NELS,NEXT,NPL1 
C       
      DOUBLE PRECISION SAVE,TEMP      
C       
C*********************************************************************
C       
C ... PREPROCESSING PHASE   
C       
C ...... DETERMINE THE NUMBER OF NONZEROES IN THE ROWS OF THE PERMUTED
C        MATRIX AND STORE THAT IN NEWIA.  THEN SWEEP THRU NEWIA TO MAKE 
C        NEWIA(I) POINT TO THE BEGINNING OF EACH ROW IN THE PERMUTED
C        DATA STRUCTURE.  ALSO NEGATE ALL THE ENTRIES IN JA TO INDICATE 
C        THAT THOSE ENTRIES HAVE NOT BEEN MOVED YET.      
C       
      N = NN      
      IER = 0     
      NPL1 = N+1  
      NELS = IA(NPL1)-1     
      DO 10 I = 1,N 
         NEWIA(I) = 0       
   10 CONTINUE    
      DO 30 I = 1,N 
         IP = P(I)
         IBGN = IA(I)       
         IEND = IA(I+1)-1   
         IF (IBGN.GT.IEND) GO TO 90   
         DO 20 J = IBGN,IEND
            IPP = IP
            JAJ = JA(J)     
            JP = P(JAJ)     
            IF (ISYM.EQ.0.AND.IP.GT.JP) IPP = JP
            NEWIA(IPP) = NEWIA(IPP)+1 
            JA(J) = -JAJ    
   20    CONTINUE 
   30 CONTINUE    
      IBGN = 1    
      DO 40 I = 1,N 
         K = IBGN+NEWIA(I)  
         NEWIA(I) = IBGN    
         IBGN = K 
   40 CONTINUE    
C       
C ...... PREPROCESSING NOW FINISHED.  
C       
C ...... NOW PERMUTE JA AND A.  THIS PERMUTATION WILL PERFORM THE   
C        FOLLOWING STEPS    
C       
C           1.  FIND THE FIRST ENTRY IN JA NOT PERMUTED WHICH IS    
C               INDICATED BY AN NEGATIVE VALUE IN JA      
C           2.  COMPUTE WHICH ROW THE CURRENT ENTRY IS IN.  THIS    
C               IS COMPUTED BY A BISECTION SEARCH THRU THE ARRAY    
C               IA. 
C           3.  USING THE PERMUTATION ARRAY P AND THE ARRAY NEWIA   
C               COMPUTE WHERE THE CURRENT ENTRY IS TO BE PLACED.    
C           4.  THEN PICK UP THE ENTRY WHERE THE CURRENT ENTRY WILL 
C               GO.  PUT THE CURRENT ENTRY IN PLACE.  THEN MAKE THE 
C               DISPLACED ENTRY THE CURRENT ENTRY AND LOOP TO STEP 2. 
C           5.  THIS PROCESS WILL END WHEN THE NEXT ENTRY HAS ALREADY 
C               BEEN MOVED.  THEN LOOP TO STEP 1. 
C       
      DO 70 J = 1,NELS      
         IF (JA(J).GT.0) GO TO 70     
         JAJ = -JA(J)       
         SAVE = A(J)
         NEXT = J 
         JA(J) = JAJ
C       
   50    JP = P(JAJ)
         I = BISRCH(NPL1,IA,NEXT)     
         IP = P(I)
         IPP = IP 
         IF (ISYM.NE.0.OR.IP.LE.JP) GO TO 60    
         IPP = JP 
         JP = IP  
   60    NEXT = NEWIA(IPP)  
C       
         TEMP = SAVE
         SAVE = A(NEXT)     
         A(NEXT) = TEMP     
C       
         JAJ = -JA(NEXT)    
         JA(NEXT) = JP      
         NEWIA(IPP) = NEWIA(IPP)+1    
         IF (JAJ.GT.0) GO TO 50       
C       
   70 CONTINUE    
C       
C ...... THE MATRIX IS NOW PERMUTED BUT THE ROWS MAY NOT BE IN      
C        ORDER.  THE REMAINDER OF THIS SUBROUTINE PERFORMS
C        A QUICK SORT ON EACH ROW TO SORT THE ENTRIES IN  
C        COLUMN ORDER.  THE IA ARRAY IS ALSO CORRECTED FROM 
C        INFORMATION STORED IN THE NEWIA ARRAY.  NEWIA(I) NOW       
C        POINTS TO THE FIRST ENTRY OF ROW I+1.  
C       
      IA(1) = 1   
      DO 80 I = 1,N 
         IA(I+1) = NEWIA(I) 
         K = IA(I+1)-IA(I)  
         IF (K.EQ.1) GO TO 80 
         IF (K.LT.1) GO TO 110
C       
         IBGN = IA(I)       
         CALL QSORT (K,JA(IBGN),A(IBGN),IER)    
         IF (IER.NE.0) GO TO 130      
C       
   80 CONTINUE    
C       
C ...... END OF MATRIX PERMUTATION    
C       
      GO TO 150   
C       
C ... ERROR TRAPS 
C       
C ...... NO ENTRY IN ROW I IN THE ORIGINAL SYSTEM 
C       
   90 IER = 301   
      IF (LEVEL.GE.0) WRITE (NOUT,100) I
  100 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE PERMAT  '/' ','    NO ENTRY IN ROW ',I10
     *   ,' OF ORIGINAL MATRIX ')     
      GO TO 150   
C       
C ...... NO ENTRY IN ROW I IN THE PERMUTED SYSTEM 
C       
  110 IER = 302   
      IF (LEVEL.GE.0) WRITE (NOUT,120) I
  120 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE PRBNDX  '/' ','    NO ENTRY IN ROW ',I10
     *   ,' OF PERMUTED MATRIX ')     
      GO TO 150   
C       
C ...... ERROR RETURN FROM SUBROUTINE QSORT     
C       
  130 IER = 303   
      IF (LEVEL.GE.0) WRITE (NOUT,140) I
  140 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE QSORT   '/' ',  
     *   '    ERROR IN SORTING PERMUTED ROW ',I12/' ',    
     *   '    CALLED FROM ITPACK ROUTINE PRBNDX   ')      
C       
  150 CONTINUE    
      IERR = IER  
      RETURN      
      END 
c===================================
      SUBROUTINE PERVEC (N,V,P)       
C       
C     THIS SUBROUTINE PERMUTES A D.P. VECTOR AS DICTATED BY THE     
C     PERMUTATION VECTOR, P.  IF P(I) = J, THEN V(J) GETS V(I).     
C       
C ... PARAMETER LIST:       
C       
C          V      D.P. VECTOR OF LENGTH N       
C          P     INTEGER PERMUTATION VECTOR     
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N,P(N)
      DOUBLE PRECISION V(N) 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER II,NEXT,NOW   
      DOUBLE PRECISION SAVE,TEMP      
C       
      IF (N.LE.0) RETURN    
C       
      DO 20 II = 1,N
         IF (P(II).LT.0) GO TO 20     
C       
         NEXT = P(II)       
         SAVE = V(II)       
C       
   10    CONTINUE 
         IF (P(NEXT).LT.0) GO TO 20   
         TEMP = SAVE
         SAVE = V(NEXT)     
         V(NEXT) = TEMP     
C       
         NOW = NEXT 
         NEXT = P(NOW)      
         P(NOW) = -NEXT     
         GO TO 10 
C       
   20 CONTINUE    
C       
      DO 30 II = 1,N
         P(II) = -P(II)     
   30 CONTINUE    
C       
      RETURN      
      END 
c=================================
      SUBROUTINE SCAL (NN,IA,JA,A,RHS,U,D,LEVEL,NOUT,IER) 
C       
C ... ORIGINAL MATRIX IS SCALED TO A UNIT DIAGONAL MATRIX.  RHS     
C ... AND U ARE SCALED ACCORDINGLY.  THE MATRIX IS THEN SPLIT AND   
C ... IA, JA, AND A RESHUFFLED.       
C       
C ... PARAMETER LIST:       
C       
C          N      DIMENSION OF MATRIX (= NN)    
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
C          U      LATEST ESTIMATE OF SOLUTION   
C          D      OUTPUT VECTOR CONTAINING THE SQUARE ROOTS 
C                    OF THE DIAGONAL ENTRIES    
C          LEVEL  PRINTING SWITCH FOR ERROR CONDITION     
C          NOUT OUTPUT TAPE NUMBER    
C          IER    ERROR FLAG: ON RETURN NONZERO VALUES MEAN 
C                    401 : THE ITH DIAGONAL ELEMENT IS .LE. 0.      
C                    402 : NO DIAGONAL ELEMENT IN ROW I   
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(*),JA(*),NN,LEVEL,NOUT,IER     
      DOUBLE PRECISION A(*),RHS(NN),U(NN),D(NN) 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,IBGN,IEND,II,IM1,J,JADD,JAJJ,JJ,JJPI,N,NP1
      DOUBLE PRECISION DI   
C       
C ... EXTRACT SQUARE ROOT OF THE DIAGONAL OUT OF A AND SCALE U AND RHS
C       
      N = NN      
      IER = 0     
      DO 80 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 50   
         DO 40 JJ = IBGN,IEND 
            IF (JA(JJ).NE.II) GO TO 40
            DI = A(JJ)      
            IF (DI.GT.0.D0) GO TO 70  
            IF (DI.EQ.0.D0) GO TO 20  
            IER = 401       
            IF (LEVEL.GE.0) WRITE (NOUT,10) II  
   10       FORMAT ('0','*** F A T A L     E R R O R ************'/'0', 
     *         '    IN ITPACK ROUTINE SCAL    '/' ',      
     *         '    DIAGONAL ENTRY IN ROW ',I10,' NEGATIVE')
            RETURN
   20       IER = 401       
            IF (LEVEL.GE.0) WRITE (NOUT,30)     
   30       FORMAT ('0','*** F A T A L     E R R O R ************'/'0', 
     *         '    IN ITPACK ROUTINE SCAL    '/' ',      
     *         '    DIAGONAL ENTRY IN ROW ',I10,' IS ZERO') 
            RETURN
   40    CONTINUE 
   50    IER = 402
         IF (LEVEL.GE.0) WRITE (NOUT,60) II     
   60    FORMAT ('0','*** F A T A L     E R R O R ************'/'0',
     *      '    IN ITPACK ROUTINE SCAL    '/' ', 
     *      '    NO DIAGONAL ENTRY IN ROW',I10) 
         RETURN   
C       
   70    CONTINUE 
         DI = DSQRT(DABS(DI)) 
         RHS(II) = RHS(II)/DI 
         U(II) = U(II)*DI   
         D(II) = DI 
   80 CONTINUE    
C       
C ... SHIFT MATRIX TO ELIMINATE DIAGONAL ENTRIES
C       
      IF (N.EQ.1) GO TO 110 
      NP1 = N+1   
      DO 100 I = 1,N
         IM1 = I-1
         II = NP1-I 
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         JADD = IBGN+IEND   
         DO 90 J = IBGN,IEND
            JJ = JADD-J     
            JJPI = JJ+IM1   
            IF (JA(JJ).EQ.II) IM1 = I 
            A(JJPI) = A(JJ) 
            JA(JJPI) = JA(JJ) 
   90    CONTINUE 
         IA(II+1) = IA(II+1)+I-1      
  100 CONTINUE    
  110 IA(1) = IA(1)+N       
C       
C ... SCALE SHIFTED MATRIX AND STORE D ARRAY IN FIRST N ENTRIES OF A
C       
      DO 140 II = 1,N       
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         DI = D(II) 
         IF (IBGN.GT.IEND) GO TO 130  
         DO 120 JJ = IBGN,IEND
            JAJJ = JA(JJ)   
            A(JJ) = A(JJ)/(DI*D(JAJJ))
  120    CONTINUE 
  130    CONTINUE 
         A(II) = DI 
  140 CONTINUE    
C       
      RETURN      
      END 
c================================
      SUBROUTINE DCOPY (N,DX,INCX,DY,INCY)      
C       
C     COPY DOUBLE PRECISION DX TO DOUBLE PRECISION DY.    
C       
      DOUBLE PRECISION DX(*),DY(*)    
      IF (N.LE.0) RETURN    
      IF (INCX.EQ.INCY) IF (INCX-1) 10 , 30 , 70
   10 CONTINUE    
C       
C        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.      
C       
      IX = 1      
      IY = 1      
      IF (INCX.LT.0) IX = (-N+1)*INCX+1 
      IF (INCY.LT.0) IY = (-N+1)*INCY+1 
      DO 20 I = 1,N 
         DY(IY) = DX(IX)    
         IX = IX+INCX       
         IY = IY+INCY       
   20 CONTINUE    
      RETURN      
C       
C        CODE FOR BOTH INCREMENTS EQUAL TO 1    
C       
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7. 
C       
   30 M = N-(N/7)*7 
      IF (M.EQ.0) GO TO 50  
      DO 40 I = 1,M 
         DY(I) = DX(I)      
   40 CONTINUE    
      IF (N.LT.7) RETURN    
   50 MP1 = M+1   
      DO 60 I = MP1,N,7     
         DY(I) = DX(I)      
         DY(I+1) = DX(I+1)  
         DY(I+2) = DX(I+2)  
         DY(I+3) = DX(I+3)  
         DY(I+4) = DX(I+4)  
         DY(I+5) = DX(I+5)  
         DY(I+6) = DX(I+6)  
   60 CONTINUE    
      RETURN      
C       
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.    
C       
   70 CONTINUE    
      NS = N*INCX 
      DO 80 I = 1,NS,INCX   
         DY(I) = DX(I)      
   80 CONTINUE    
      RETURN      
      END 
c==========================================
      SUBROUTINE PJAC (NN,IA,JA,A,U,RHS)
C       
C     ... THIS SUBROUTINE PERFORMS ONE JACOBI ITERATION.  
C       
C ... PARAMETER LIST:       
C       
C          N      DIMENSION OF MATRIX (= NN)    
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          U      ESTIMATE OF SOLUTION OF A MATRIX PROBLEM
C          RHS    ON INPUT: CONTAINS THE RIGHT HAND SIDE OF 
C                    A MATRIX PROBLEM 
C                 ON OUTPUT: CONTAINS A*U + RHS 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(*),JA(*),NN
      DOUBLE PRECISION A(*),U(NN),RHS(NN)       
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IBGN,IEND,II,JAJJ,JJ,N  
      DOUBLE PRECISION RHSII,UII      
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
      N = NN      
      IF (ISYM.EQ.0) GO TO 30 
C       
C     *************** NON - SYMMETRIC SECTION ****************      
C       
      DO 20 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 20   
         RHSII = RHS(II)    
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHSII = RHSII-A(JJ)*U(JAJJ) 
   10    CONTINUE 
         RHS(II) = RHSII    
   20 CONTINUE    
      RETURN      
C       
C     ************** SYMMETRIC SECTION **********************       
C       
   30 DO 50 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 50   
         RHSII = RHS(II)    
         UII = U(II)
         DO 40 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHSII = RHSII-A(JJ)*U(JAJJ) 
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UII     
   40    CONTINUE 
         RHS(II) = RHSII    
   50 CONTINUE    
      RETURN      
C       
      END 
c========================================
      SUBROUTINE VEVMW (N,V,W)
C       
C ... VEVMW COMPUTES V = V - W
C       
C ... PARAMETER LIST:       
C       
C          N      INTEGER LENGTH OF VECTORS V AND W       
C          V      D.P. VECTOR 
C          W      D.P. VECTOR SUCH THAT   V(I) = V(I) - W(I)
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N   
      DOUBLE PRECISION V(N),W(N)      
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,M,MP1       
C       
      IF (N.LE.0) RETURN    
      M = MOD(N,4)
C       
      IF (M.EQ.0) GO TO 20  
      DO 10 I = 1,M 
         V(I) = V(I)-W(I)   
   10 CONTINUE    
      IF (N.LT.4) RETURN    
C       
   20 MP1 = M+1   
      DO 30 I = MP1,N,4     
         V(I) = V(I)-W(I)   
         V(I+1) = V(I+1)-W(I+1)       
         V(I+2) = V(I+2)-W(I+2)       
         V(I+3) = V(I+3)-W(I+3)       
   30 CONTINUE    
      RETURN      
C       
      END 
c==============================
      SUBROUTINE ITJCG (NN,IA,JA,A,U,U1,D,D1,DTWD,TRI)    
C       
C ... FUNCTION:   
C       
C          THIS SUBROUTINE, ITJCG, PERFORMS ONE ITERATION OF THE    
C          JACOBI CONJUGATE GRADIENT ALGORITHM.  IT IS CALLED BY JCG. 
C       
C ... PARAMETER LIST:       
C       
C          N      INPUT INTEGER.  DIMENSION OF THE MATRIX. (= NN)   
C          IA,JA  INPUT INTEGER VECTORS.  CONTAINS INFORMATION DEFINING 
C                 THE SPARSE MATRIX REPRESENTATION.       
C          A      INPUT D.P. VECTOR. CONTAINS THE NONZERO VALUES OF THE 
C                 LINEAR SYSTEM.      
C          U      INPUT D.P. VECTOR.  CONTAINS THE VALUE OF THE     
C                 SOLUTION VECTOR AT THE END OF IN ITERATIONS.      
C          U1     INPUT/OUTPUT D.P. VECTOR.  ON INPUT, IT CONTAINS  
C                 THE VALUE OF THE SOLUTION AT THE END OF THE IN-1  
C                 ITERATION.  ON OUTPUT, IT WILL CONTAIN THE NEWEST 
C                 ESTIMATE FOR THE SOLUTION VECTOR.       
C          D      INPUT D.P. VECTOR.  CONTAINS THE PSEUDO-RESIDUAL  
C                 VECTOR AFTER IN ITERATIONS.   
C          D1     INPUT/OUTPUT D.P. VECTOR.  ON INPUT, D1 CONTAINS  
C                 THE PSEUDO-RESIDUAL VECTOR AFTER IN-1 ITERATIONS.  ON 
C                 OUTPUT, IT WILL CONTAIN THE NEWEST PSEUDO-RESIDUAL
C                 VECTOR.   
C          DTWD   D.P. ARRAY.  USED IN THE COMPUTATIONS OF THE      
C                 ACCELERATION PARAMETER GAMMA AND THE NEW PSEUDO-  
C                 RESIDUAL. 
C          TRI    D.P. ARRAY.  STORES THE TRIDIAGONAL MATRIX ASSOCIATED 
C                 WITH THE EIGENVALUES OF THE CONJUGATE GRADIENT    
C                 POLYNOMIAL. 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(*),JA(*),NN
      DOUBLE PRECISION A(*),U(NN),U1(NN),D(NN),D1(NN),DTWD(NN),TRI(2,1) 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER N   
      DOUBLE PRECISION CON,C1,C2,C3,C4,DNRM,DTNRM,GAMOLD,RHOOLD,RHOTMP
      LOGICAL Q1  
C       
C ... SPECIFICATIONS FOR FUNCTION SUBPROGRAMS   
C       
      DOUBLE PRECISION DDOT 
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN SUBROUTINE JCG   
C       
C ... COMPUTE NEW ESTIMATE FOR CME IF ADAPT = .TRUE.      
C       
      IF (ADAPT) CALL CHGCON (TRI,GAMOLD,RHOOLD,1)
C       
C ... TEST FOR STOPPING     
C       
      N = NN      
      DELNNM = DDOT(N,D,1,D,1)
      DNRM = DELNNM 
      CON = CME   
      CALL PSTOP (N,U,DNRM,CON,1,Q1)  
      IF (HALT) GO TO 30    
C       
C ... COMPUTE RHO AND GAMMA - ACCELERATION PARAMETERS     
C       
      CALL VFILL (N,DTWD,0.D0)
      CALL PJAC (N,IA,JA,A,D,DTWD)    
      DTNRM = DDOT(N,D,1,DTWD,1)      
      IF (ISYM.EQ.0) GO TO 10 
      RHOTMP = DDOT(N,DTWD,1,D1,1)    
      CALL PARCON (DTNRM,C1,C2,C3,C4,GAMOLD,RHOTMP,1)     
      RHOOLD = RHOTMP       
      GO TO 20    
   10 CALL PARCON (DTNRM,C1,C2,C3,C4,GAMOLD,RHOOLD,1)     
C       
C ... COMPUTE U(IN+1) AND D(IN+1)     
C       
   20 CALL SUM3 (N,C1,D,C2,U,C3,U1)   
      CALL SUM3 (N,C1,DTWD,C4,D,C3,D1)
C       
C ... OUTPUT INTERMEDIATE INFORMATION 
C       
   30 CALL ITERM (N,A,U,DTWD,1)       
C       
      RETURN      
      END 
c==================================
      SUBROUTINE UNSCAL (N,IA,JA,A,RHS,U,D)     
C       
C ... THIS SUBROUTINE REVERSES THE PROCESS OF SCAL.       
C       
C ... PARAMETER LIST:       
C       
C          N      DIMENSION OF MATRIX 
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
C          U      LATEST ESTIMATE OF SOLUTION   
C          D      VECTOR CONTAINING THE SQUARE ROOTS      
C                    OF THE DIAGONAL ENTRIES    
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(*),JA(*),N 
      DOUBLE PRECISION A(*),RHS(N),U(N),D(N)    
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IBGN,IEND,II,INEW,IS,JAJJ,JJ,JJPI 
      DOUBLE PRECISION DI   
C       
C ... EXTRACT DIAGONAL FROM SCALED A AND UNSCALE U AND RHS
C       
      DO 10 II = 1,N
         DI = A(II) 
         U(II) = U(II)/DI   
         RHS(II) = RHS(II)*DI 
         D(II) = DI 
   10 CONTINUE    
C       
C ... UNSCALE A   
C       
      DO 30 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 30   
         DI = D(II) 
         DO 20 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            A(JJ) = A(JJ)*DI*D(JAJJ)  
   20    CONTINUE 
   30 CONTINUE    
C       
C ... INSERT DIAGONAL BACK INTO A     
C       
      DO 60 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IS = N-II
         INEW = IBGN-IS-1   
         A(INEW) = D(II)**2 
         JA(INEW) = II      
         IF (IS.EQ.0.OR.IBGN.GT.IEND) GO TO 50  
         DO 40 JJ = IBGN,IEND 
            JJPI = JJ-IS    
            A(JJPI) = A(JJ) 
            JA(JJPI) = JA(JJ) 
   40    CONTINUE 
   50    CONTINUE 
         IA(II) = INEW      
   60 CONTINUE    
C       
      RETURN      
      END 
c=============================
      INTEGER FUNCTION BISRCH (N,K,L) 
C       
C ... BISRCH IS AN INTEGER FUNCTION WHICH USES A BISECTION SEARCH   
C     TO FIND THE ENTRY J IN THE ARRAY K SUCH THAT THE VALUE L IS   
C     GREATER THAN OR EQUAL TO K(J) AND STRICTLY LESS THAN K(J+1).  
C       
C ... PARAMETER LIST:       
C       
C          N      INTEGER LENGTH OF VECTOR K    
C          K      INTEGER VECTOR      
C          L      INTEGER CONSTANT SUCH THAT  K(J) .GE. L .LT. K(J+1) 
C                 WITH J RETURNED AS VALUE OF INTEGER FUNCTION BISRCH 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N,L,K(N)      
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER JLEFT,JMID,JRIGHT       
C       
      JLEFT = 1   
      JRIGHT = N  
      IF (N.EQ.2) GO TO 40  
      JMID = (N+1)/2
C       
   10 IF (L.GE.K(JMID)) GO TO 20      
C       
C ...... L .GE. K(LEFT)  AND  L .LT. K(JMID)    
C       
      JRIGHT = JMID 
      GO TO 30    
C       
C ...... L .GE. K(JMID)  AND  L .LT. K(JRIGHT)  
C       
   20 JLEFT = JMID
C       
C ...... TEST FOR CONVERGENCE 
C       
   30 IF (JRIGHT-JLEFT.EQ.1) GO TO 40 
      JMID = JLEFT+(JRIGHT-JLEFT+1)/2 
      GO TO 10    
C       
C ...... BISECTION SEARCH FINISHED    
C       
   40 BISRCH = JLEFT
C       
      RETURN      
      END 
c============================
      SUBROUTINE CHGCON (TRI,GAMOLD,RHOOLD,IBMTH) 
C       
C     COMPUTES THE NEW ESTIMATE FOR THE LARGEST EIGENVALUE FOR      
C     CONJUGATE GRADIENT ACCELERATION.
C       
C ... PARAMETER LIST:       
C       
C          TRI    TRIDIAGONAL MATRIX ASSOCIATED WITH THE EIGENVALUES
C                    OF THE CONJUGATE GRADIENT POLYNOMIAL 
C          GAMOLD 
C            AND  
C          RHOOLD PREVIOUS VALUES OF ACCELERATION PARAMETERS
C          IBMTH  INDICATOR OF BASIC METHOD BEING ACCELERATED BY CG 
C                      IBMTH = 1,  JACOBI       
C                            = 2,  REDUCED SYSTEM 
C                            = 3,  SSOR 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IBMTH 
      DOUBLE PRECISION TRI(2,1),GAMOLD,RHOOLD   
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IB2,IB3,IER,IP
      DOUBLE PRECISION CMOLD,END,START,EIGVSS,EIGVNS      
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
      GO TO (10,20,30), IBMTH 
C       
C ... JACOBI CONJUGATE GRADIENT       
C       
   10 START = CME 
      IP = IN     
      GO TO 40    
C       
C ... REDUCED SYSTEM CG     
C       
   20 START = CME**2
      IP = IN     
      GO TO 40    
C       
C ... SSOR CG     
C       
   30 IF (ADAPT) START = SPR
      IF (.NOT.ADAPT) START = SPECR   
      IP = IN-IS  
C       
C ... DEFINE THE MATRIX     
C       
   40 IF (IP.GE.2) GO TO 60 
      IF (IP.EQ.1) GO TO 50 
C       
C ... IP = 0      
C       
      END = 0.D0  
      CMOLD = 0.D0
      GO TO 110   
C       
C ... IP = 1      
C       
   50 END = 1.D0-1.D0/GAMMA 
      TRI(1,1) = END
      TRI(2,1) = 0.D0       
      GO TO 110   
C       
C ... IP > 1      
C       
   60 IF ((IP.GT.2).AND.(DABS(START-CMOLD).LE.ZETA*START)) GO TO 120
      CMOLD = START 
C       
C ... COMPUTE THE LARGEST EIGENVALUE  
C       
      TRI(1,IP) = 1.D0-1.D0/GAMMA     
      TRI(2,IP) = (RHO-1.D0)/(RHO*RHOOLD*GAMMA*GAMOLD)    
      IF (ISYM.NE.0) GO TO 80 
      END = EIGVSS(IP,TRI,START,ZETA,ITMAX,IER) 
      IF (IER.EQ.0) GO TO 100 
      IF (LEVEL.GE.2) WRITE (NOUT,70) IER       
   70 FORMAT (/10X,'DIFFICULTY IN COMPUTATION OF MAXIMUM EIGENVALUE'/15X
     *   ,'OF ITERATION MATRIX'/10X,'SUBROUTINE ZBRENT RETURNED IER =', 
     *   I5)      
      GO TO 100   
   80 IB2 = 1+IP  
      IB3 = IB2+IP/2+1      
      END = EIGVNS(IP,TRI,TRI(1,IB2),TRI(1,IB3),IER)      
      IF (IER.EQ.0) GO TO 100 
      IF (LEVEL.GE.2) WRITE (NOUT,90) IER       
   90 FORMAT (/10X,'DIFFICULTY IN COMPUTATION OF MAXIMUM EIGENVALUE'/15X
     *   ,'OF ITERATION MATRIX'/10X,'SUBROUTINE EQRT1S RETURNED IER =', 
     *   I5)      
  100 CONTINUE    
      IF (IER.NE.0) GO TO 130 
C       
C ... SET SPECTRAL RADIUS FOR THE VARIOUS METHODS 
C       
  110 IF (IBMTH.EQ.1) CME = END       
      IF (IBMTH.EQ.2) CME = DSQRT(DABS(END))    
      IF (IBMTH.EQ.3.AND.ADAPT) SPR = END       
      IF (IBMTH.EQ.3.AND..NOT.ADAPT) SPECR = END
      RETURN      
C       
C ... RELATIVE CHANGE IN CME IS LESS THAN ZETA.  THEREFORE STOP     
C     CHANGING.   
C       
  120 ADAPT = .FALSE.       
      PARTAD = .FALSE.      
      RETURN      
C       
C ... ESTIMATE FOR CME > 1.D0.  THEREFORE NEED TO STOP ADAPTIVE     
C     PROCEDURE AND KEEP OLD VALUE OF CME.      
C       
  130 ADAPT = .FALSE.       
      PARTAD = .FALSE.      
      IF (LEVEL.GE.2) WRITE (NOUT,140) IN,START 
  140 FORMAT (/10X,'ESTIMATE OF MAXIMUM EIGENVALUE OF JACOBI   '/15X, 
     *   'MATRIX (CME) NOT ACCURATE'/10X,       
     *   'ADAPTIVE PROCEDURE TURNED OFF AT ITERATION ',I5/10X,      
     *   'FINAL ESTIMATE OF MAXIMUM EIGENVALUE =',D15.7/) 
C       
      RETURN      
      END 
c==================================
      DOUBLE PRECISION FUNCTION DDOT (N,DX,INCX,DY,INCY)  
C       
C     RETURNS THE DOT PRODUCT OF DOUBLE PRECISION DX AND DY.
C       
      DOUBLE PRECISION DX(*),DY(*)    
      DDOT = 0.D0 
      IF (N.LE.0) RETURN    
      IF (INCX.EQ.INCY) IF (INCX-1) 10 , 30 , 70
   10 CONTINUE    
C       
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.     
C       
      IX = 1      
      IY = 1      
      IF (INCX.LT.0) IX = (-N+1)*INCX+1 
      IF (INCY.LT.0) IY = (-N+1)*INCY+1 
      DO 20 I = 1,N 
         DDOT = DDOT+DX(IX)*DY(IY)    
         IX = IX+INCX       
         IY = IY+INCY       
   20 CONTINUE    
      RETURN      
C       
C        CODE FOR BOTH INCREMENTS EQUAL TO 1.   
C       
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5. 
C       
   30 M = N-(N/5)*5 
      IF (M.EQ.0) GO TO 50  
      DO 40 I = 1,M 
         DDOT = DDOT+DX(I)*DY(I)      
   40 CONTINUE    
      IF (N.LT.5) RETURN    
   50 MP1 = M+1   
      DO 60 I = MP1,N,5     
         DDOT = DDOT+DX(I)*DY(I)+DX(I+1)*DY(I+1)+DX(I+2)*DY(I+2)+DX(I+3)
     *      *DY(I+3)+DX(I+4)*DY(I+4)  
   60 CONTINUE    
      RETURN      
C       
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.       
C       
   70 CONTINUE    
      NS = N*INCX 
      DO 80 I = 1,NS,INCX   
         DDOT = DDOT+DX(I)*DY(I)      
   80 CONTINUE    
      RETURN      
      END 
c=====================================
      SUBROUTINE PSTOP (N,U,DNRM,CCON,IFLAG,Q1) 
C       
C     THIS SUBROUTINE PERFORMS A TEST TO SEE IF THE ITERATIVE       
C     METHOD HAS CONVERGED TO A SOLUTION INSIDE THE ERROR 
C     TOLERANCE, ZETA.      
C       
C ... PARAMETER LIST:       
C       
C          N      ORDER OF SYSTEM     
C          U      PRESENT SOLUTION ESTIMATE     
C          DNRM   INNER PRODUCT OF PSEUDO-RESIDUALS AT PRECEDING    
C                    ITERATION
C          CON    STOPPING TEST PARAMETER (= CCON)
C          IFLAG  STOPPING TEST INTEGER FLAG    
C                    IFLAG = 0,  SOR ITERATION ZERO       
C                    IFLAG = 1,  NON-RS METHOD  
C                    IFLAG = 2,  RS METHOD      
C          Q1     STOPPING TEST LOGICAL FLAG    
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N,IFLAG       
      DOUBLE PRECISION U(N),DNRM,CCON 
      LOGICAL Q1  
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      DOUBLE PRECISION CON,TL,TR,UOLD 
C       
C ... SPECIFICATIONS FOR ARGUMENT SUBROUTINES   
C       
      DOUBLE PRECISION DDOT 
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
      CON = CCON  
      HALT = .FALSE.
C       
C     SPECIAL PROCEDURE FOR ZEROTH ITERATION    
C       
      IF (IN.GE.1) GO TO 10 
      Q1 = .FALSE.
      UDNM = 1.D0 
      STPTST = 1.D3 
      IF (IFLAG.LE.0) RETURN
C       
C ... TEST IF UDNM NEEDS TO BE RECOMPUTED       
C       
   10 CONTINUE    
      IF (Q1) GO TO 20      
      IF ((IN.GT.5).AND.(MOD(IN,5).NE.0)) GO TO 20
      UOLD = UDNM 
      UDNM = DDOT(N,U,1,U,1)
      IF (UDNM.EQ.0.D0) UDNM = 1.D0   
      IF ((IN.GT.5).AND.(DABS(UDNM-UOLD).LE.UDNM*ZETA)) Q1 = .TRUE. 
C       
C ... COMPUTE STOPPING TEST 
C       
   20 TR = DSQRT(UDNM)      
      TL = 1.D0   
      IF (CON.EQ.1.D0) GO TO 40       
      IF (IFLAG.EQ.2) GO TO 30
      TL = DSQRT(DNRM)      
      TR = TR*(1.D0-CON)    
      GO TO 40    
   30 TL = DSQRT(2.D0*DNRM) 
      TR = TR*(1.D0-CON*CON)
   40 STPTST = TL/TR
      IF (TL.GE.TR*ZETA) RETURN       
      HALT = .TRUE. 
C       
      RETURN      
      END 
c========================================
      SUBROUTINE PARCON (DTNRM,C1,C2,C3,C4,GAMOLD,RHOTMP,IBMTH)     
C       
C     COMPUTES ACCELERATION PARAMETERS FOR CONJUGATE GRADIENT       
C     ACCELERATED METHODS.  
C       
C ... PARAMETER LIST:       
C       
C          DTNRM  INNER PRODUCT OF RESIDUALS    
C          C1     OUTPUT: RHO*GAMMA   
C          C2     OUTPUT: RHO 
C          C3     OUTPUT: 1-RHO       
C          C4     OUTPUT: RHO*(1-GAMMA) 
C          GAMOLD OUTPUT: VALUE OF GAMMA AT PRECEDING ITERATION     
C          RHOTMP LAST ESTIMATE FOR VALUE OF RHO
C          IBMTH  INDICATOR OF BASIC METHOD BEING ACCELERATED BY CG 
C                      IBMTH = 1,   JACOBI      
C                            = 2,   REDUCED SYSTEM
C                            = 3,   SSOR
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IBMTH 
      DOUBLE PRECISION DTNRM,C1,C2,C3,C4,GAMOLD,RHOTMP    
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IP  
      DOUBLE PRECISION RHOOLD 
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
      IP = IN-IS  
C       
C ... SET RHOOLD AND GAMOLD 
C       
      RHOOLD = RHO
      GAMOLD = GAMMA
C       
C ... COMPUTE GAMMA (IN+1)  
C       
C ... FOR JACOBI OR REDUCED SYSTEM CG 
C       
      IF (IBMTH.LE.2) GAMMA = 1.D0/(1.D0-DTNRM/DELNNM)    
C       
C ... FOR SSOR CG 
C       
      IF (IBMTH.EQ.3) GAMMA = DELNNM/DTNRM      
C       
C ... COMPUTE RHO (IN+1)    
C       
      RHO = 1.D0  
      IF (IP.EQ.0) GO TO 20 
      IF (ISYM.EQ.0) GO TO 10 
      RHO = 1.D0/(1.D0-GAMMA*RHOTMP/DELSNM)     
      GO TO 20    
   10 RHO = 1.D0/(1.D0-GAMMA*DELNNM/(GAMOLD*DELSNM*RHOOLD)) 
C       
C ... COMPUTE CONSTANTS C1, C2, C3, AND C4      
C       
   20 DELSNM = DELNNM       
      RHOTMP = RHOOLD       
      C1 = RHO*GAMMA
      C2 = RHO    
      C3 = 1.D0-RHO 
      C4 = RHO*(1.D0-GAMMA) 
C       
      RETURN      
      END 
c====================================
      SUBROUTINE SUM3 (N,C1,X1,C2,X2,C3,X3)     
C       
C ... COMPUTES X3 = C1*X1 + C2*X2 + C3*X3       
C       
C ... PARAMETER LIST:       
C       
C          N        INTEGER LENGTH OF VECTORS X1, X2, X3  
C          C1,C2,C3 D.P. CONSTANTS    
C          X1,X2,X3 D.P. VECTORS SUCH THAT      
C                   X3(I) = C1*X1(I) + C2*X2(I) + C3*X3(I)
C                   X3(I) = C1*X1(I) + C2*X2(I)  IF C3 = 0. 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I   
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
      DOUBLE PRECISION X1(N),X2(N),X3(N),C1,C2,C3 
C       
      IF (N.LE.0) RETURN    
      IF (DABS(C3).EQ.0.D0) GO TO 20  
C       
      DO 10 I = 1,N 
         X3(I) = C1*X1(I)+C2*X2(I)+C3*X3(I)     
   10 CONTINUE    
      RETURN      
C       
C ... COMPUTE X3 = C1*X1 + C2*X2      
C       
   20 DO 30 I = 1,N 
         X3(I) = C1*X1(I)+C2*X2(I)    
   30 CONTINUE    
C       
      RETURN      
      END 
c================================
      SUBROUTINE ITERM (NN,A,U,WK,IMTHDD)       
C       
C     THIS ROUTINE PRODUCES THE ITERATION SUMMARY LINE AT THE END   
C     OF EACH ITERATION. IF LEVEL = 5, THE LATEST APPROXIMATION     
C     TO THE SOLUTION WILL BE PRINTED.
C       
C ... PARAMETER LIST:       
C       
C          NN     ORDER OF SYSTEM OR, FOR REDUCED SYSTEM  
C                    ROUTINES, ORDER OF BLACK SUBSYSTEM   
C          A      ITERATION MATRIX    
C          U      SOLUTION ESTIMATE   
C          WK     WORK ARRAY OF LENGTH NN       
C          IMTHD  INDICATOR OF METHOD (=IMTHDD) 
C                    IMTHD = 1,  JCG  
C                    IMTHD = 2,  JSI  
C                    IMTHD = 3,  SOR  
C                    IMTHD = 4,  SSORCG 
C                    IMTHD = 5,  SSORSI 
C                    IMTHD = 6,  RSCG 
C                    IMTHD = 7,  RSSI 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER NN,IMTHD      
      DOUBLE PRECISION A(*),U(NN),WK(NN)
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,IMTHDD,IP,N 
      DOUBLE PRECISION QTFF 
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
      N = NN      
      IMTHD = IMTHDD
C       
C ... PRINT VARIOUS PARAMETERS AFTER EACH ITERATION       
C       
      IF (LEVEL.LT.2) RETURN
      GO TO (10,110,170,210,50,10,110), IMTHD   
   10 IF (IN.GT.0) GO TO 30 
C       
C ... PRINT HEADER FOR JCG AND RSCG   
C       
      WRITE (NOUT,20)       
   20 FORMAT (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'//  
     *   ' NUMBER OF',5X,'CONVERGENCE',7X,'CME ',11X,'RHO',12X,'GAMMA'/ 
     *   ' ITERATIONS',4X,'TEST '//)  
C       
C ... PRINT SUMMARY LINE    
C       
   30 WRITE (NOUT,40) IN,STPTST,CME,RHO,GAMMA   
   40 FORMAT (4X,I5,3X,4D15.7)
      IF (LEVEL.GE.4) GO TO 250       
C       
      RETURN      
C       
   50 IF (IN.GT.0) GO TO 70 
C       
C ... PRINT HEADER FOR SSOR-SI
C       
      WRITE (NOUT,60)       
   60 FORMAT (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'//  
     *   ' NUMBER OF',4X,'CONVERGENCE',7X,'PARAMETER CHANGE TEST',10X,
     *   'RHO',12X,'GAMMA'/' ITERATIONS',3X,'TEST ',11X,'LHS(QA)',7X, 
     *   'RHS(QT**FF)'//)   
C       
C ... PRINT SUMMARY LINE    
C       
   70 IP = IN-IS  
      IF (IMTHD.EQ.7) IP = 2*IP       
      IF (IP.LT.3) GO TO 90 
      QTFF = QT**FF 
      WRITE (NOUT,80) IN,STPTST,QA,QTFF,RHO,GAMMA 
   80 FORMAT (4X,I5,3X,5D15.7)
      IF (LEVEL.GE.4) GO TO 250       
      RETURN      
C       
   90 WRITE (NOUT,100) IN,STPTST,RHO,GAMMA      
  100 FORMAT (4X,I5,3X,D15.7,30X,2D15.7)
      IF (LEVEL.GE.4) GO TO 250       
      RETURN      
C       
  110 IF (IN.GT.0) GO TO 130
C       
C ... PRINT HEADER FOR J-SI AND RS-SI 
C       
      WRITE (NOUT,120)      
  120 FORMAT (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'//  
     *   ' NUMBER OF',4X,'CONVERGENCE',7X,'PARAMETER CHANGE TEST',10X,
     *   'RHO'/' ITERATIONS',3X,'TEST ',11X,'LHS(QA)',7X,'RHS(QT**FF)'//
     *   )
C       
C ... PRINT SUMMARY LINE    
C       
  130 IP = IN-IS  
      IF (IMTHD.EQ.7) IP = 2*IP       
      IF (IP.LT.3) GO TO 150
      QTFF = QT**FF 
      WRITE (NOUT,140) IN,STPTST,QA,QTFF,RHO    
  140 FORMAT (4X,I5,3X,5D15.7)
      IF (LEVEL.GE.4) GO TO 250       
      RETURN      
C       
  150 WRITE (NOUT,160) IN,STPTST,RHO  
  160 FORMAT (4X,I5,3X,D15.7,30X,D15.7) 
      IF (LEVEL.GE.4) GO TO 250       
      RETURN      
C       
C ... PRINT VARIOUS PARAMETERS AFTER EACH ITERATION FOR SOR.
C       
  170 IF (IN.GT.0) GO TO 190
C       
C ... PRINT HEADER FOR SOR  
C       
      WRITE (NOUT,180)      
  180 FORMAT (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'//  
     *   ' NUMBER OF',4X,'CONVERGENCE',6X,'CME ',9X,'OMEGA',7X,     
     *   'SPECTRAL'/' ITERATIONS',3X,'TEST',38X,'RADIUS'//) 
C       
C ... PRINT SUMMARY LINE FOR SOR      
C       
  190 CONTINUE    
      WRITE (NOUT,200) IN,STPTST,CME,OMEGA,SPECR
  200 FORMAT (4X,I5,3X,4D14.7)
      IF (LEVEL.GE.4) GO TO 250       
C       
      RETURN      
C       
C ... PRINT VARIOUS PARAMETERS AFTER EACH ITERATION FOR SSOR-CG.    
C       
  210 IF (IN.GT.0) GO TO 230
C       
C ... PRINT HEADER FOR SSOR-CG
C       
      WRITE (NOUT,220)      
  220 FORMAT (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'//  
     *   ' NUMBER OF',4X,'CONVERGENCE',3X,' SPECTRAL',6X,'S-PRIME',9X,
     *   'RHO',10X,'GAMMA'/' ITERATIONS',3X,'TEST ',10X,'RADIUS'//) 
C       
C ... PRINT SUMMARY LINE FOR SSOR-CG  
C       
  230 CONTINUE    
      WRITE (NOUT,240) IN,STPTST,SPECR,SPR,RHO,GAMMA      
  240 FORMAT (4X,I5,3X,5D14.7)
      IF (LEVEL.GE.4) GO TO 250       
      RETURN      
C       
  250 IF (IMTHD.GT.5) GO TO 270       
      WRITE (NOUT,260) IN   
  260 FORMAT ('0',2X,'ESTIMATE OF SOLUTION AT ITERATION ',I5)       
      GO TO 290   
  270 WRITE (NOUT,280) IN   
  280 FORMAT ('0',2X,'ESTIMATE OF SOLUTION AT BLACK POINTS ',       
     *   'AT ITERATION ',I5)
  290 DO 300 I = 1,N
         WK(I) = U(I)/A(I)  
  300 CONTINUE    
      WRITE (NOUT,310) (WK(I),I=1,N)  
  310 FORMAT (2X,5(2X,D20.13))
      WRITE (NOUT,320)      
  320 FORMAT (//) 
C       
      RETURN      
      END 
c=======================================
      DOUBLE PRECISION FUNCTION EIGVSS (N,TRI,START,ZETA,ITMAX,IER) 
C       
C     COMPUTES THE LARGEST EIGENVALUE OF A SYMMETRIC TRIDIAGONAL MATRIX 
C     FOR CONJUGATE GRADIENT ACCELERATION.      
C     MODIFIED IMSL ROUTINE ZBRENT USED.
C       
C ... PARAMETER LIST:       
C       
C          N      ORDER OF TRIDIAGONAL SYSTEM   
C          TRI    SYMMETRIC TRIDIAGONAL MATRIX OF ORDER N 
C          START  INITIAL LOWER BOUND OF INTERVAL CONTAINING ROOT   
C          ZETA   STOPPING CRITERIA   
C          IER    ERROR FLAG: ON RETURN, IER=0 INDICATES THAT       
C                    THE LARGEST EIGENVALUE OF TRI WAS FOUND.       
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N,ITMAX,IER   
      DOUBLE PRECISION TRI(2,1),START,ZETA,A,B,EPS
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER MAXFN,NSIG,ITMP 
C       
      EIGVSS = 0.D0 
      ITMP = IFIX(SNGL(-DLOG10(DABS(ZETA))))    
      NSIG = MAX0(ITMP,4)   
      MAXFN = MAX0(ITMAX,50)
C       
C     EPS = DMIN1(ZETA,0.5D-4)
C       
      EPS = 0.0D0 
      A = START   
      B = 1.0D0   
      CALL ZBRENT (N,TRI,EPS,NSIG,A,B,MAXFN,IER)
      EIGVSS = B  
C       
      RETURN      
      END 
c====================================
      DOUBLE PRECISION FUNCTION EIGVNS (N,TRI,D,E2,IER)   
C       
C     COMPUTES THE LARGEST EIGENVALUE OF A SYMMETRIC TRIDIAGONAL MATRIX 
C     FOR CONJUGATE GRADIENT ACCELERATION.      
C       
C ... PARAMETER LIST:       
C       
C          N      ORDER OF TRIDIAGONAL SYSTEM   
C          TRI    SYMMETRIC TRIDIAGONAL MATRIX OF ORDER N 
C          D      ARRAY FOR EQRT1S (NEGATIVE DIAGONAL ELEMENTS)     
C          E2     ARRAY FOR EQRT1S (SUPER DIAGONAL ELEMENTS)
C          IER    ERROR FLAG: ON RETURN, IER=0 INDICATES THAT       
C                    THE LARGEST EIGENVALUE OF TRI WAS FOUND.       
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N,IER 
      DOUBLE PRECISION TRI(2,1),D(N),E2(N)      
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I   
C       
      EIGVNS = 0.D0 
C       
      D(1) = -TRI(1,1)      
      DO 10 I = 2,N 
         D(I) = -TRI(1,I)   
         E2(I) = DABS(TRI(2,I))       
   10 CONTINUE    
C       
      CALL EQRT1S (D,E2,N,1,0,IER)    
      EIGVNS = -D(1)
C       
      RETURN      
      END 
c=============================
      SUBROUTINE ZBRENT (N,TRI,EPS,NSIG,AA,BB,MAXFNN,IER) 
C       
C   MODIFIED IMSL ROUTINE NAME   - ZBRENT       
C       
C-----------------------------------------------------------------------
C       
C   COMPUTER            - CDC/SINGLE  
C       
C   LATEST REVISION     - JANUARY 1, 1978       
C       
C   PURPOSE             - ZERO OF A FUNCTION WHICH CHANGES SIGN IN A
C                           GIVEN INTERVAL (BRENT ALGORITHM)
C       
C   USAGE               - CALL ZBRENT (F,EPS,NSIG,A,B,MAXFN,IER)    
C       
C   ARGUMENTS    TRI    - A TRIDIAGONAL MATRIX OF ORDER N 
C                EPS    - FIRST CONVERGENCE CRITERION (INPUT).  A ROOT, 
C                           B, IS ACCEPTED IF DABS(F(B)) IS LESS THAN OR
C                           EQUAL TO EPS.  EPS MAY BE SET TO ZERO.  
C                NSIG   - SECOND CONVERGENCE CRITERION (INPUT).  A ROOT,
C                           B, IS ACCEPTED IF THE CURRENT APPROXIMATION 
C                           AGREES WITH THE TRUE SOLUTION TO NSIG   
C                           SIGNIFICANT DIGITS. 
C                A,B    - ON INPUT, THE USER MUST SUPPLY TWO POINTS, A
C                           AND B, SUCH THAT F(A) AND F(B) ARE OPPOSITE 
C                           IN SIGN. (= AA, BB) 
C                           ON OUTPUT, BOTH A AND B ARE ALTERED.  B 
C                           WILL CONTAIN THE BEST APPROXIMATION TO THE
C                           ROOT OF F. SEE REMARK 1.      
C                MAXFN  - ON INPUT, MAXFN SHOULD CONTAIN AN UPPER BOUND 
C                           ON THE NUMBER OF FUNCTION EVALUATIONS   
C                           REQUIRED FOR CONVERGENCE.  ON OUTPUT, MAXFN 
C                           WILL CONTAIN THE ACTUAL NUMBER OF FUNCTION
C                           EVALUATIONS USED. (= MAXFNN)  
C                IER    - ERROR PARAMETER. (OUTPUT)       
C                         TERMINAL ERROR
C                           IER = 501 INDICATES THE ALGORITHM FAILED TO 
C                             CONVERGE IN MAXFN EVALUATIONS.
C                           IER = 502 INDICATES F(A) AND F(B) HAVE THE
C                             SAME SIGN.
C       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32 
C                       - SINGLE/H36,H48,H60    
C       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND       
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL  
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C       
C   REMARKS  1.  LET F(X) BE THE CHARACTERISTIC FUNCTION OF THE MATRIX
C                TRI EVALUATED AT X. FUNCTION DETERM EVALUATES F(X).
C                ON EXIT FROM ZBRENT, WHEN IER=0, A AND B SATISFY THE 
C                FOLLOWING, 
C                F(A)*F(B) .LE.0,     
C                DABS(F(B)) .LE. DABS(F(A)), AND
C                EITHER DABS(F(B)) .LE. EPS OR  
C                DABS(A-B) .LE. MAX(DABS(B),0.1)*10.0**(-NSIG).     
C                THE PRESENCE OF 0.1 IN THIS ERROR CRITERION CAUSES 
C                LEADING ZEROES TO THE RIGHT OF THE DECIMAL POINT TO BE 
C                COUNTED AS SIGNIFICANT DIGITS. SCALING MAY BE REQUIRED 
C                IN ORDER TO ACCURATELY DETERMINE A ZERO OF SMALL   
C                MAGNITUDE. 
C            2.  ZBRENT IS GUARANTEED TO REACH CONVERGENCE WITHIN   
C                K = (DLOG((B-A)/D)+1.0)**2 FUNCTION EVALUATIONS WHERE
C                  D=MIN(OVER X IN (A,B) OF     
C                    MAX(DABS(X),0.1)*10.0**(-NSIG)).     
C                THIS IS AN UPPER BOUND ON THE NUMBER OF EVALUATIONS. 
C                RARELY DOES THE ACTUAL NUMBER OF EVALUATIONS USED BY 
C                ZBRENT EXCEED DSQRT(K). D CAN BE COMPUTED AS FOLLOWS,
C                  P = DBLE(AMIN1(DABS(A),DABS(B)))       
C                  P = DMAX1(0.1,P)   
C                  IF ((A-0.1)*(B-0.1).LT.0.0) P = 0.1    
C                  D = P*10.0**(-NSIG)
C       
C   COPYRIGHT           - 1977 BY IMSL, INC. ALL RIGHTS RESERVED.   
C       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.    
C       
C-----------------------------------------------------------------------
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
C                                  SPECIFICATIONS FOR ARGUMENTS     
C       
      INTEGER NSIG,MAXFNN,IER 
      DOUBLE PRECISION TRI(2,1),EPS,AA,BB       
C       
C                                  SPECIFICATIONS FOR LOCAL VARIABLES 
C       
      INTEGER IC,MAXFN      
      DOUBLE PRECISION ZERO,HALF,ONE,THREE,TEN,A,B,T,FA,FB,C,FC,D,E,TOL,
     *   RM,S,P,Q,R,RONE,TEMP,DETERM  
      DATA ZERO / 0.D0 / ,HALF / 5.D-1 / ,ONE / 1.D0 / ,THREE / 3.D0 / ,
     *   TEN / 10.D0 /      
C       
C                                  FIRST EXECUTABLE STATEMENT       
C       
      A = AA      
      B = BB      
      MAXFN = MAXFNN
      IER = 0     
      T = TEN**(-NSIG)      
      IC = 2      
      FA = DETERM(N,TRI,A)  
      FB = DETERM(N,TRI,B)  
      S = B       
C       
C                                  TEST FOR SAME SIGN     
C       
      IF (FA*FB.GT.ZERO) GO TO 110    
   10 C = A       
      FC = FA     
      D = B-C     
      E = D       
   20 IF (DABS(FC).GE.DABS(FB)) GO TO 30
      A = B       
      B = C       
      C = A       
      FA = FB     
      FB = FC     
      FC = FA     
   30 CONTINUE    
      TOL = T*DMAX1(DABS(B),0.1D0)    
      RM = (C-B)*HALF       
C       
C                                  TEST FOR FIRST CONVERGENCE CRITERIA
C       
      IF (DABS(FB).LE.EPS) GO TO 80   
C       
C                                  TEST FOR SECOND CONVERGENCE CRITERIA 
C       
      IF (DABS(C-B).LE.TOL) GO TO 80  
C       
C                                  CHECK EVALUATION COUNTER 
C       
      IF (IC.GE.MAXFN) GO TO 90       
C       
C                                  IS BISECTION FORCED    
C       
      IF (DABS(E).LT.TOL) GO TO 60    
      IF (DABS(FA).LE.DABS(FB)) GO TO 60
      S = FB/FA   
      IF (A.NE.C) GO TO 40  
C       
C                                  LINEAR INTERPOLATION   
C       
      P = (C-B)*S 
      Q = ONE-S   
      GO TO 50    
C       
C                                  INVERSE QUADRATIC INTERPOLATION  
C       
   40 Q = FA/FC   
      R = FB/FC   
      RONE = R-ONE
      P = S*((C-B)*Q*(Q-R)-(B-A)*RONE)
      Q = (Q-ONE)*RONE*(S-ONE)
   50 IF (P.GT.ZERO) Q = -Q 
      IF (P.LT.ZERO) P = -P 
      S = E       
      E = D       
C       
C                                  IF DABS(P/Q).GE.75*DABS(C-B) THEN
C                                     FORCE BISECTION     
C       
      IF (P+P.GE.THREE*RM*Q) GO TO 60 
C       
C                                  IF DABS(P/Q).GE..5*DABS(S) THEN FORCE
C                                     BISECTION. S = THE VALUE OF P/Q 
C                                     ON THE STEP BEFORE THE LAST ONE 
C       
      IF (P+P.GE.DABS(S*Q)) GO TO 60  
      D = P/Q     
      GO TO 70    
C       
C                                  BISECTION    
C       
   60 E = RM      
      D = E       
C       
C                                  INCREMENT B  
C       
   70 A = B       
      FA = FB     
      TEMP = D    
      IF (DABS(TEMP).LE.HALF*TOL) TEMP = DSIGN(HALF*TOL,RM) 
      B = B+TEMP  
      S = B       
      FB = DETERM(N,TRI,S)  
      IC = IC+1   
      IF (FB*FC.LE.ZERO) GO TO 20     
      GO TO 10    
C       
C                                  CONVERGENCE OF B       
C       
   80 A = C       
      MAXFN = IC  
      GO TO 130   
C       
C                                  MAXFN EVALUATIONS      
C       
   90 IER = 501   
      A = C       
      MAXFN = IC  
      IF (LEVEL.GE.1) WRITE (NOUT,100) MAXFN    
  100 FORMAT ('0','*** W A R N I N G ************'/'0',   
     *   '    IN ITPACK ROUTINE ZBRENT  '/' ',  
     *   '    ALGORITHM FAILED TO CONVERGE   '/' ','    IN',I6,     
     *   ' ITERATIONS ')    
      GO TO 130   
C       
C                                  TERMINAL ERROR - F(A) AND F(B) HAVE
C                                  THE SAME SIGN
C       
  110 IER = 502   
      MAXFN = IC  
      IF (LEVEL.GE.1) WRITE (NOUT,120)
  120 FORMAT ('0','*** W A R N I N G ************'/'0',   
     *   '    IN ITPACK ROUTINE ZBRENT  '/' ',  
     *   '    F(A) AND F(B) HAVE SAME SIGN   ') 
  130 CONTINUE    
      AA = A      
      BB = B      
      MAXFNN = MAXFN
      RETURN      
      END 
c=========================================
      SUBROUTINE EQRT1S (D,E2,NN,M,ISW,IERR)    
C       
C   MODIFIED IMSL ROUTINE NAME   - EQRT1S       
C       
C-----------------------------------------------------------------------
C       
C   COMPUTER            - CDC/SINGLE  
C       
C   LATEST REVISION     - JUNE 1, 1980
C       
C   PURPOSE             - SMALLEST OR LARGEST M EIGENVALUES OF A    
C                           SYMMETRIC TRIDIAGONAL MATRIX  
C       
C   USAGE               - CALL EQRT1S (D,E2,N,M,ISW,IER)  
C       
C   ARGUMENTS    D      - INPUT VECTOR OF LENGTH N CONTAINING       
C                           THE DIAGONAL ELEMENTS OF THE MATRIX.  THE 
C                           COMPUTED EIGENVALUES REPLACE THE FIRST M
C                           COMPONENTS OF THE VECTOR D IN NON-      
C                           DECREASING SEQUENCE, WHILE THE REMAINING
C                           COMPONENTS ARE LOST.
C                E2     - INPUT VECTOR OF LENGTH N CONTAINING       
C                           THE SQUARES OF THE OFF-DIAGONAL ELEMENTS
C                           OF THE MATRIX.  INPUT E2 IS DESTROYED.  
C                N      - INPUT SCALAR CONTAINING THE ORDER OF THE  
C                           MATRIX. (= NN)      
C                M      - INPUT SCALAR CONTAINING THE NUMBER OF     
C                           SMALLEST EIGENVALUES DESIRED (M IS      
C                           LESS THAN OR EQUAL TO N).     
C                ISW    - INPUT SCALAR MEANING AS FOLLOWS - 
C                           ISW=1 MEANS THAT THE MATRIX IS KNOWN TO BE
C                             POSITIVE DEFINITE.
C                           ISW=0 MEANS THAT THE MATRIX IS NOT KNOWN
C                             TO BE POSITIVE DEFINITE.    
C                IER    - ERROR PARAMETER. (OUTPUT) (= IERR)
C                           WARNING ERROR       
C                             IER = 601 INDICATES THAT SUCCESSIVE   
C                               ITERATES TO THE K-TH EIGENVALUE WERE NOT
C                               MONOTONE INCREASING. THE VALUE K IS 
C                               STORED IN E2(1).
C                           TERMINAL ERROR      
C                             IER = 602 INDICATES THAT ISW=1 BUT MATRIX 
C                               IS NOT POSITIVE DEFINITE  
C       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32 
C                       - SINGLE/H36,H48,H60    
C       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND       
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL  
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C       
C   REMARKS      AS WRITTEN, THE ROUTINE COMPUTES THE M SMALLEST    
C                EIGENVALUES. TO COMPUTE THE M LARGEST EIGENVALUES, 
C                REVERSE THE SIGN OF EACH ELEMENT OF D BEFORE AND   
C                AFTER CALLING THE ROUTINE. IN THIS CASE, ISW MUST  
C                EQUAL ZERO.
C       
C   COPYRIGHT           - 1980 BY IMSL, INC. ALL RIGHTS RESERVED.   
C       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.    
C       
C-----------------------------------------------------------------------
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
C                                  SPECIFICATIONS FOR ARGUMENTS     
C       
      INTEGER NN,M,ISW,IERR 
      DOUBLE PRECISION D(NN),E2(NN)   
C       
C                                  SPECIFICATIONS FOR LOCAL VARIABLES 
C       
      INTEGER II,I,JJ,J,K1,K,N,IER    
      DOUBLE PRECISION DELTA,DLAM,EP,ERR,F,P,QP,Q,R,S,TOT 
C       
C                                  DRELPR = MACHINE PRECISION       
C                                  FIRST EXECUTABLE STATEMENT       
C       
      N = NN      
      IER = 0     
      DLAM = 0.0D0
      ERR = 0.0D0 
      S = 0.0D0   
C       
C                                  LOOK FOR SMALL SUB-DIAGONAL ENTRIES
C                                  DEFINE INITIAL SHIFT FROM LOWER  
C                                  GERSCHGORIN BOUND.     
C       
      TOT = D(1)  
      Q = 0.0D0   
      J = 0       
      DO 30 I = 1,N 
         P = Q    
         IF (I.EQ.1) GO TO 10 
         IF (P.GT.DRELPR*(DABS(D(I))+DABS(D(I-1)))) GO TO 20
   10    E2(I) = 0.0D0      
C       
C                                  COUNT IF E2(I) HAS UNDERFLOWED   
C       
   20    IF (E2(I).EQ.0.D0) J = J+1   
         Q = 0.0D0
         IF (I.NE.N) Q = DSQRT(DABS(E2(I+1)))   
         TOT = DMIN1(D(I)-P-Q,TOT)    
   30 CONTINUE    
      IF (ISW.EQ.1.AND.TOT.LT.0.0D0) GO TO 50   
      DO 40 I = 1,N 
         D(I) = D(I)-TOT    
   40 CONTINUE    
      GO TO 60    
   50 TOT = 0.0D0 
   60 DO 200 K = 1,M
C       
C                                  NEXT QR TRANSFORMATION 
C       
   70    TOT = TOT+S
         DELTA = D(N)-S     
         I = N    
         F = DABS(DRELPR*TOT) 
         IF (DLAM.LT.F) DLAM = F      
         IF (DELTA.GT.DLAM) GO TO 90  
         IF (DELTA.GE.(-DLAM)) GO TO 170
         IER = 602
         IF (LEVEL.GE.1) WRITE (NOUT,80)
   80    FORMAT ('0','*** W A R N I N G ************'/' ',
     *      '    IN ITPACK ROUTINE EQRT1S  '/' ', 
     *      '    PARAMETER ISW = 1 BUT MATRIX   '/' ',    
     *      '    NOT POSITIVE DEFINITE')
         GO TO 210
C       
C                                  REPLACE SMALL SUB-DIAGONAL SQUARES 
C                                  BY ZERO TO REDUCE THE INCIDENCE OF 
C                                  UNDERFLOWS   
C       
   90    IF (K.EQ.N) GO TO 110
         K1 = K+1 
         DO 100 J = K1,N    
            IF (E2(J).LE.(DRELPR*(D(J)+D(J-1)))**2) E2(J) = 0.0D0   
  100    CONTINUE 
  110    F = E2(N)/DELTA    
         QP = DELTA+F       
         P = 1.0D0
         IF (K.EQ.N) GO TO 140
         K1 = N-K 
         DO 130 II = 1,K1   
            I = N-II
            Q = D(I)-S-F    
            R = Q/QP
            P = P*R+1.0D0   
            EP = F*R
            D(I+1) = QP+EP  
            DELTA = Q-EP    
            IF (DELTA.GT.DLAM) GO TO 120
            IF (DELTA.GE.(-DLAM)) GO TO 170     
            IER = 602       
            IF (LEVEL.GE.0) WRITE (NOUT,80)     
            GO TO 210       
  120       F = E2(I)/Q     
            QP = DELTA+F    
            E2(I+1) = QP*EP 
  130    CONTINUE 
  140    D(K) = QP
         S = QP/P 
         IF (TOT+S.GT.TOT) GO TO 70   
         IER = 601
         E2(1) = K
         IF (LEVEL.GE.1) WRITE (NOUT,150) K     
  150    FORMAT ('0','*** W A R N I N G ************'/'0',
     *      '    IN ITPACK ROUTINE EQRT1S  '/' ', 
     *      '    SUCCESSIVE ITERATES TO THE',I10/' ',     
     *      '    EIGENVALUE WERE NOT MONOTONE INCREASING ') 
C       
C                                  SET ERROR -- IRREGULAR END       
C                                  DEFLATE MINIMUM DIAGONAL ELEMENT 
C       
         S = 0.0D0
         DELTA = QP 
         DO 160 J = K,N     
            IF (D(J).GT.DELTA) GO TO 160
            I = J 
            DELTA = D(J)    
  160    CONTINUE 
C       
C                                  CONVERGENCE  
C       
  170    IF (I.LT.N) E2(I+1) = E2(I)*F/QP       
         IF (I.EQ.K) GO TO 190
         K1 = I-K 
         DO 180 JJ = 1,K1   
            J = I-JJ
            D(J+1) = D(J)-S 
            E2(J+1) = E2(J) 
  180    CONTINUE 
  190    D(K) = TOT 
         ERR = ERR+DABS(DELTA)
         E2(K) = ERR
  200 CONTINUE    
      IF (IER.EQ.0) GO TO 220 
  210 CONTINUE    
  220 IERR = IER  
      RETURN      
      END 
c==========================
      DOUBLE PRECISION FUNCTION DETERM (N,TRI,XLMDA)      
C       
C     THIS SUBROUTINE COMPUTES THE DETERMINANT OF A SYMMETRIC       
C     TRIDIAGONAL MATRIX GIVEN BY TRI. DET(TRI - XLMDA*I) = 0       
C       
C ... PARAMETER LIST
C       
C          N      ORDER OF TRIDIAGONAL SYSTEM   
C          TRI    SYMMETRIC TRIDIAGONAL MATRIX OF ORDER N 
C          XLMDA  ARGUMENT FOR CHARACTERISTIC EQUATION    
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N   
      DOUBLE PRECISION TRI(2,1),XLMDA 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER ICNT,L,NM1    
      DOUBLE PRECISION D1,D2,D3       
C       
      NM1 = N-1   
      D2 = TRI(1,N)-XLMDA   
      D1 = D2*(TRI(1,NM1)-XLMDA)-TRI(2,N)       
      IF (N.EQ.2) GO TO 20  
C       
C ... BEGINNING OF LOOP     
C       
      DO 10 ICNT = 2,NM1    
         L = NM1-ICNT+2     
         D3 = D2  
         D2 = D1  
         D1 = (TRI(1,L-1)-XLMDA)*D2-D3*TRI(2,L) 
   10 CONTINUE    
C       
C ... DETERMINANT COMPUTED  
C       
   20 DETERM = D1 
C       
      RETURN      
      END 
c===============================
      SUBROUTINE PERROR1 (NN,IA,JA,A,RHS,U,W,DIGTT1,DIGTT2,IDGTTS)   
C       
C     PERROR1 COMPUTES THE RESIDUAL, R = RHS - A*U.  THE USER
C     ALSO HAS THE OPTION OF PRINTING THE RESIDUAL AND/OR THE       
C     UNKNOWN VECTOR DEPENDING ON IDGTS.
C       
C ... PARAMETER LIST:       
C       
C          N      DIMENSION OF MATRIX (= NN)    
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
C          U      LATEST ESTIMATE OF SOLUTION   
C          W      WORKSPACE VECTOR    
C          DIGIT1 OUTPUT: MEASURE OF ACCURACY OF STOPPING TEST (= DIGTT1
C          DIGIT2 OUTPUT: MEASURE OF ACCURACY OF SOLUTION (= DIGTT2)
C          IDGTS   PARAMETER CONTROLING LEVEL OF OUTPUT (= IDGTTS)  
C                    IF IDGTS < 1 OR IDGTS > 4, THEN NO OUTPUT.     
C                            = 1, THEN NUMBER OF DIGITS IS PRINTED, PRO-
C                                 VIDED LEVEL .GE. 1      
C                            = 2, THEN SOLUTION VECTOR IS PRINTED, PRO- 
C                                 VIDED LEVEL .GE. 1      
C                            = 3, THEN RESIDUAL VECTOR IS PRINTED, PRO- 
C                                 VIDED LEVEL .GE. 1      
C                            = 4, THEN BOTH VECTORS ARE PRINTED, PRO- 
C                                 VIDED LEVEL .GE. 1      
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(*),JA(*),NN,IDGTTS   
      DOUBLE PRECISION A(*),RHS(NN),U(NN),W(NN),DIGTT1,DIGTT2       
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IDGTS,N       
      DOUBLE PRECISION BNRM,DIGIT1,DIGIT2,RNRM,TEMP       
C       
C ... SPECIFICATIONS FOR FUNCTION SUBPROGRAMS   
C       
      DOUBLE PRECISION DDOT 
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
      N = NN      
      IDGTS = IDGTTS
      DIGIT1 = 0.D0 
      DIGIT2 = 0.D0 
      IF (N.LE.0) GO TO 40  
C       
      DIGIT1 = -DLOG10(DABS(DRELPR))  
      IF (STPTST.GT.0.D0) DIGIT1 = -DLOG10(DABS(STPTST))  
      BNRM = DDOT(N,RHS,1,RHS,1)      
      IF (BNRM.EQ.0.D0) GO TO 10      
      CALL PMULT (N,IA,JA,A,U,W)      
      CALL WEVMW (N,RHS,W)  
      RNRM = DDOT(N,W,1,W,1)
      TEMP = RNRM/BNRM      
      IF (TEMP.EQ.0.D0) GO TO 10      
      DIGIT2 = -DLOG10(DABS(TEMP))/2.D0 
      GO TO 20    
C       
   10 DIGIT2 = -DLOG10(DABS(DRELPR))  
C       
   20 IF ((IDGTS.LT.1).OR.(LEVEL.LE.0)) GO TO 40
      WRITE (NOUT,30) DIGIT1,DIGIT2   
   30 FORMAT (/6X,'APPROX. NO. OF DIGITS (EST. REL. ERROR) =',F5.1,2X,
     *   '(DIGIT1)'/3X,'APPROX. NO. OF DIGITS (EST. REL. RESIDUAL) =',
     *   F5.1,2X,'(DIGIT2)')
C       
      IF (IDGTS.LE.1.OR.IDGTS.GT.4) GO TO 40    
      IF (IDGTS.NE.3) CALL VOUT (N,U,2,NOUT)    
      IF (IDGTS.GE.3) CALL VOUT (N,W,1,NOUT)    
C       
   40 CONTINUE    
      DIGTT1 = DIGIT1       
      DIGTT2 = DIGIT2       
      RETURN      
      END 
c=======================
      SUBROUTINE QSORT (NN,KEY,DATA,ERROR)      
C       
C     ==================================================================
C       
C     Q U I C K S O R T     
C       
C         IN THE STYLE OF THE CACM PAPER BY BOB SEDGEWICK, OCTOBER 1978 
C       
C     INPUT:      
C         N    -- NUMBER OF ELEMENTS TO BE SORTED (= NN)  
C         KEY  -- AN ARRAY OF LENGTH  N  CONTAINING THE VALUES      
C                 WHICH ARE TO BE SORTED
C         DATA -- A SECOND ARRAY OF LENGTH  N  CONTAINING DATA      
C                 ASSOCIATED WITH THE INDIVIDUAL KEYS.    
C       
C     OUTPUT:     
C         KEY  -- WILL BE ARRANGED SO THAT VALUES ARE IN INCREASING 
C                 ORDER     
C         DATA -- REARRANGED TO CORRESPOND TO REARRANGED KEYS       
C         ERROR -- WILL BE ZERO UNLESS YOUR INPUT FILE WAS OF TRULY 
C                  ENORMOUS LENGTH, IN WHICH CASE IT WILL BE EQUAL TO 1.
C       
C     ==================================================================
C       
      INTEGER NN,ERROR,KEY(NN)
      DOUBLE PRECISION DATA(NN)       
C       
C     ------------------------
C       
      INTEGER TOP,LEFT,RIGHT,I,J,TINY,V,K,IP1,JM1,LLEN,RLEN,N       
      LOGICAL DONE
      DOUBLE PRECISION D    
      INTEGER STKLEN,STACK(30)
C       
      DATA TINY,STKLEN / 9,30 /       
C       
C     -----------------------------------       
C       
C     ... PROGRAM IS A DIRECT TRANSLATION INTO FORTRAN OF SEDGEWICK^S 
C         PROGRAM 2, WHICH IS NON-RECURSIVE, IGNORES FILES OF LENGTH
C         LESS THAN 'TINY' DURING PARTITIONING, AND USES MEDIAN OF THREE
C         PARTITIONING.     
C       
      N = NN      
      IF (N.EQ.1) RETURN    
      IF (N.LE.0) GO TO 240 
C       
      ERROR = 0   
      TOP = 1     
      LEFT = 1    
      RIGHT = N   
      DONE = (N.LE.TINY)    
C       
      IF (DONE) GO TO 150   
      CALL IVFILL (STKLEN,STACK,0)    
C       
C     ===========================================================   
C     QUICKSORT -- PARTITION THE FILE UNTIL NO SUBFILE REMAINS OF   
C     LENGTH GREATER THAN 'TINY'      
C     ===========================================================   
C       
C     ... WHILE NOT DONE DO ...       
C       
   10 IF (DONE) GO TO 150   
C       
C         ... FIND MEDIAN OF LEFT, RIGHT AND MIDDLE ELEMENTS OF CURRENT 
C             SUBFILE, WHICH IS  KEY(LEFT), ..., KEY(RIGHT) 
C       
      LFRH2 = (LEFT+RIGHT)/2
      K = KEY(LFRH2)
      D = DATA(LFRH2)       
      KEY(LFRH2) = KEY(LEFT)
      DATA(LFRH2) = DATA(LEFT)
      KEY(LEFT) = K 
      DATA(LEFT) = D
C       
      IF (KEY(LEFT+1).LE.KEY(RIGHT)) GO TO 20   
      K = KEY(LEFT+1)       
      D = DATA(LEFT+1)      
      KEY(LEFT+1) = KEY(RIGHT)
      DATA(LEFT+1) = DATA(RIGHT)      
      KEY(RIGHT) = K
      DATA(RIGHT) = D       
C       
   20 IF (KEY(LEFT).LE.KEY(RIGHT)) GO TO 30     
      K = KEY(LEFT) 
      D = DATA(LEFT)
      KEY(LEFT) = KEY(RIGHT)
      DATA(LEFT) = DATA(RIGHT)
      KEY(RIGHT) = K
      DATA(RIGHT) = D       
C       
   30 IF (KEY(LEFT+1).LE.KEY(LEFT)) GO TO 40    
      K = KEY(LEFT+1)       
      D = DATA(LEFT+1)      
      KEY(LEFT+1) = KEY(LEFT) 
      DATA(LEFT+1) = DATA(LEFT)       
      KEY(LEFT) = K 
      DATA(LEFT) = D
C       
   40 V = KEY(LEFT) 
C       
C         ... V IS NOW THE MEDIAN VALUE OF THE THREE KEYS.  NOW MOVE
C             FROM THE LEFT AND RIGHT ENDS SIMULTANEOUSLY, EXCHANGING 
C             KEYS AND DATA UNTIL ALL KEYS LESS THAN  V  ARE PACKED TO
C             THE LEFT, ALL KEYS LARGER THAN  V  ARE PACKED TO THE  
C             RIGHT.
C       
      I = LEFT+1  
      J = RIGHT   
C       
C         LOOP    
C             REPEAT I = I+1 UNTIL KEY(I) >= V; 
C             REPEAT J = J-1 UNTIL KEY(J) <= V; 
C         EXIT IF J < I;    
C             << EXCHANGE KEYS I AND J >>       
C         END     
C       
   50 CONTINUE    
   60 I = I+1     
      IF (KEY(I).LT.V) GO TO 60       
C       
   70 J = J-1     
      IF (KEY(J).GT.V) GO TO 70       
C       
      IF (J.LT.I) GO TO 80  
      K = KEY(I)  
      D = DATA(I) 
      KEY(I) = KEY(J)       
      DATA(I) = DATA(J)     
      KEY(J) = K  
      DATA(J) = D 
      GO TO 50    
C       
   80 K = KEY(LEFT) 
      D = DATA(LEFT)
      KEY(LEFT) = KEY(J)    
      DATA(LEFT) = DATA(J)  
      KEY(J) = K  
      DATA(J) = D 
C       
C         ... WE HAVE NOW PARTITIONED THE FILE INTO TWO SUBFILES,   
C             ONE IS (LEFT ... J-1)  AND THE OTHER IS (I...RIGHT).  
C             PROCESS THE SMALLER NEXT.  STACK THE LARGER ONE.      
C       
      LLEN = J-LEFT 
      RLEN = RIGHT-I+1      
      IF (MAX0(LLEN,RLEN).GT.TINY) GO TO 100    
C       
C             ... BOTH SUBFILES ARE TINY, SO UNSTACK NEXT LARGER FILE 
C       
      IF (TOP.EQ.1) GO TO 90
      TOP = TOP-2 
      LEFT = STACK(TOP)     
      RIGHT = STACK(TOP+1)  
      GO TO 10    
C       
   90 DONE = .TRUE. 
C       
      GO TO 10    
C       
C             ... ELSE ONE OR BOTH SUBFILES ARE LARGE     
C       
  100 IF (MIN0(LLEN,RLEN).GT.TINY) GO TO 120    
C       
C             ... ONE SUBFILE IS SMALL, ONE LARGE.  IGNORE THE SMALL ONE
C       
      IF (LLEN.GT.RLEN) GO TO 110     
      LEFT = I    
      GO TO 10    
C       
  110 RIGHT = J-1 
C       
      GO TO 10    
C       
C         ... ELSE BOTH ARE LARGER THAN TINY.  ONE MUST BE STACKED. 
C       
  120 IF (TOP.GE.STKLEN) GO TO 240    
      IF (LLEN.GT.RLEN) GO TO 130     
      STACK(TOP) = I
      STACK(TOP+1) = RIGHT  
      RIGHT = J-1 
      GO TO 140   
C       
  130 STACK(TOP) = LEFT     
      STACK(TOP+1) = J-1    
      LEFT = I    
C       
  140 TOP = TOP+2 
C       
      GO TO 10    
C       
C     ------------------------------------------------------------  
C     INSERTION SORT THE ENTIRE FILE, WHICH CONSISTS OF A LIST      
C     OF 'TINY' SUBFILES, LOCALLY OUT OF ORDER, GLOBALLY IN ORDER.  
C     ------------------------------------------------------------  
C       
C     ... FIRST, FIND LARGEST ELEMENT IN 'KEY'  
C       
  150 I = N-1     
      LEFT = MAX0(0,N-TINY) 
      K = KEY(N)  
      J = N       
C       
  160 IF (I.LE.LEFT) GO TO 180
      IF (KEY(I).LE.K) GO TO 170      
      K = KEY(I)  
      J = I       
C       
  170 I = I-1     
      GO TO 160   
C       
  180 IF (J.EQ.N) GO TO 190 
C       
C     ... LARGEST ELEMENT WILL BE IN  KEY(N)    
C       
      KEY(J) = KEY(N)       
      KEY(N) = K  
      D = DATA(N) 
      DATA(N) = DATA(J)     
      DATA(J) = D 
C       
C     ... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...       
C       
  190 I = N-1     
      IP1 = N     
C       
  200 IF (KEY(I).LE.KEY(IP1)) GO TO 220 
C       
C             ... OUT OF ORDER ... MOVE UP TO CORRECT PLACE 
C       
      K = KEY(I)  
      D = DATA(I) 
      J = IP1     
      JM1 = I     
C       
C             ... REPEAT ... UNTIL 'CORRECT PLACE FOR K FOUND'      
C       
  210 KEY(JM1) = KEY(J)     
      DATA(JM1) = DATA(J)   
      JM1 = J     
      J = J+1     
      IF (KEY(J).LT.K) GO TO 210      
C       
      KEY(JM1) = K
      DATA(JM1) = D 
C       
  220 IP1 = I     
      I = I-1     
      IF (I.GT.0) GO TO 200 
C       
  230 RETURN      
C       
  240 ERROR = 1   
      GO TO 230   
C       
      END 
c=================================
      SUBROUTINE PMULT (NN,IA,JA,A,U,W) 
C       
C     ... THIS SUBROUTINE PERFORMS ONE MATRIX-VECTOR MULTIPLICATION.
C       
C ... PARAMETER LIST:       
C       
C          N      DIMENSION OF MATRIX (= NN)    
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          U      LATEST ESTIMATE OF SOLUTION   
C          W      ON RETURN W CONTAINS A*U      
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(*),JA(*),NN
      DOUBLE PRECISION A(*),U(NN),W(NN) 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IBGN,IEND,II,JJ,N       
      DOUBLE PRECISION SUM,UII,WII    
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
      N = NN      
      IF (N.LE.0) RETURN    
      IF (ISYM.EQ.0) GO TO 40 
C       
C     *************** NON - SYMMETRIC SECTION **********************
C       
      DO 30 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = 0.0D0
         IF (IBGN.GT.IEND) GO TO 20   
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM+A(JJ)*U(JAJJ)   
   10    CONTINUE 
   20    W(II) = SUM
   30 CONTINUE    
      RETURN      
C       
C     ***************** SYMMETRIC SECTION **************************
C       
   40 CALL VFILL (N,W,0.D0) 
      DO 70 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         UII = U(II)
         WII = W(II)
         IF (IBGN.GT.IEND) GO TO 60   
         DO 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            WII = WII+A(JJ)*U(JAJJ)   
            W(JAJJ) = W(JAJJ)+A(JJ)*UII 
   50    CONTINUE 
   60    W(II) = WII
   70 CONTINUE    
      RETURN      
C       
      END 
c============================
      SUBROUTINE WEVMW (N,V,W)
C       
C ... WEVMW COMPUTES W = V - W
C       
C ... PARAMETER LIST:       
C       
C          N      INTEGER LENGTH OF VECTORS V AND W       
C          V      D.P. VECTOR 
C          W      D.P. VECTOR SUCH THAT   W(I) = V(I) - W(I)
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N   
      DOUBLE PRECISION V(N),W(N)      
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,M,MP1       
C       
      IF (N.LE.0) RETURN    
      M = MOD(N,4)
      IF (M.EQ.0) GO TO 20  
      DO 10 I = 1,M 
         W(I) = V(I)-W(I)   
   10 CONTINUE    
      IF (N.LT.4) RETURN    
C       
   20 MP1 = M+1   
      DO 30 I = MP1,N,4     
         W(I) = V(I)-W(I)   
         W(I+1) = V(I+1)-W(I+1)       
         W(I+2) = V(I+2)-W(I+2)       
         W(I+3) = V(I+3)-W(I+3)       
   30 CONTINUE    
C       
      RETURN      
      END 
c==============================
      SUBROUTINE VOUT (N,V,ISWT,NOUTT)
C       
C     THIS SUBROUTINE EFFECTS PRINTING OF RESIDUAL AND SOLUTION     
C     VECTORS - CALLED FROM PERROR1    
C       
C ... PARAMETER LIST:       
C       
C          V      VECTOR OF LENGTH N  
C          ISWT   LABELLING INFORMATION 
C          NOUT OUTPUT DEVICE NUMBER (= NOUTT)  
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N,ISWT,NOUTT  
      DOUBLE PRECISION V(N) 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,J,JM1,K,KUPPER,NOUT   
C       
      NOUT = NOUTT
C       
C        IF (N .LE. 0) RETURN 
C       
      KUPPER = MIN0(N,8)    
      IF (ISWT.EQ.1) WRITE (NOUT,10)  
   10 FORMAT (//5X,'RESIDUAL VECTOR') 
      IF (ISWT.EQ.2) WRITE (NOUT,20)  
   20 FORMAT (//5X,'SOLUTION VECTOR') 
      WRITE (NOUT,30) (I,I=1,KUPPER)  
   30 FORMAT (10X,8I15)     
      WRITE (NOUT,40)       
   40 FORMAT (10X,120('-')/)
C       
      DO 60 J = 1,N,8       
         KUPPER = MIN0(J+7,N) 
         JM1 = J-1
         WRITE (NOUT,50) JM1,(V(K),K=J,KUPPER)  
   50    FORMAT (4X,I5,'+  ',8D15.5)  
   60 CONTINUE    
C       
      RETURN      
      END 
c====================================
      SUBROUTINE SBAGN (N,NZ,IA,JA,A,IWORK,LEVELL,NOUTT,IERR)       
C       
C ... THE ROUTINES SBINI, SBSIJ, AND SBEND CREATE A SPARSE
C     MATRIX STRUCTURE BY MEANS OF A LINKED LIST WHICH IS 
C     DESTROYED BY SBEND. SBAGN CREATES A NEW LINKED LIST 
C     SO THAT ELEMENTS MAY BE ADDED TO THE MATRIX AFTER SBEND       
C     HAS BEEN CALLED. SBAGN SHOULD BE CALLED WITH THE APPRO-       
C     PRIATE PARAMETERS, AND THEN SBSIJ AND SBEND CAN BE CALLED     
C     TO ADD THE ELEMENTS AND COMPLETE THE SPARSE MATRIX STRUC-     
C     TURE.       
C       
C ... PARAMETER LIST:       
C       
C           N       ORDER OF THE SYSTEM 
C           NZ      MAXIMUM NUMBER OF NON-ZERO ELEMENTS   
C                   IN THE SYSTEM     
C           IA, JA  INTEGER ARRAYS OF THE SPARSE
C                   MATRIX STRUCTURE  
C           A       D.P. ARRAY OF THE SPARSE MATRIX       
C                   STRUCTURE 
C           IWORK   WORK ARRAY OF DIMENSION NZ  
C           LEVEL   OUTPUT LEVEL CONTROL (= LEVELL)       
C           NOUT  OUTPUT FILE NUMBER (= NOUTT)  
C           IER     ERROR FLAG (= IERR). POSSIBLE RETURNS ARE       
C                      IER = 0, SUCCESSFUL COMPLETION     
C                          = 703, NZ TOO SMALL - NO MORE  
C                                 ELEMENTS CAN BE ADDED   
C       
C ... SPECIFICTIONS FOR ARGUMENTS     
C       
      INTEGER NZ,IA(*),JA(*),IWORK(NZ),N,LEVELL,NOUTT,IERR
      DOUBLE PRECISION A(NZ)
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,IER,J,LEVEL,NOUT,NADD,NADDP1,NOW,NP1,NTO,NTN
C       
C ... INITIALIZE LOCAL VARIABLES AND MAKE ERROR CHECK     
C       
      NOW = IA(N+1)-1       
      NADD = NZ-NOW 
      IER = 0     
      LEVEL = LEVELL
      NOUT = NOUTT
      IF (NADD.LE.0) IER = 703
      IF (IER.EQ.0) GO TO 20
      IF (LEVEL.GE.0) WRITE (NOUT,10) IER       
   10 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE SBAGN   '/' ','    IER = ',I10/' ', 
     *   '    NZ TOO SMALL - NO ROOM FOR NEW ENTRY')      
      GO TO 90    
C       
C ... SHIFT ELEMENTS OF A AND JA DOWN AND ADD ZERO FILL   
C       
   20 NTO = NOW   
      NTN = NZ    
      DO 30 I = 1,NOW       
         JA(NTN) = JA(NTO)  
         A(NTN) = A(NTO)    
         NTO = NTO-1
         NTN = NTN-1
   30 CONTINUE    
      DO 40 I = 1,NADD      
         JA(I) = 0
         A(I) = 0.D0
   40 CONTINUE    
C       
C ... UPDATE IA TO REFLECT DOWNWARD SHIFT IN A AND JA     
C       
      NP1 = N+1   
      DO 50 I = 1,NP1       
         IA(I) = IA(I)+NADD 
   50 CONTINUE    
C       
C ... CREATE LINKED LIST    
C       
      NADDP1 = NADD+1       
      DO 60 I = NADDP1,NZ   
         IWORK(I) = I+1     
   60 CONTINUE    
      DO 70 I = 1,NADD      
         IWORK(I) = 0       
   70 CONTINUE    
      DO 80 I = 1,N 
         J = IA(I+1)-1      
         IWORK(J) = -I      
   80 CONTINUE    
C       
C ... INDICATE IN LAST POSITION OF IA HOW MANY SPACES     
C     ARE LEFT IN A AND JA FOR ADDITION OF ELEMENTS       
C       
      IA(N+1) = NADD
      RETURN      
C       
C ... ERROR RETURN
C       
   90 IERR = IER  
      RETURN      
      END 
c=================================
      subroutine lookup25(rlapse,rlat1,rlong1,elev1,tnsl,
     &                  tsurf,acc,abl)
      real*8 rlapse,rlat1,rlong1,elev1,tnsl,tsurf,acc,abl,WARMING
      parameter(nlat=80,nlong=360,ntot=nlat*nlong,
     &          nmon=12,nyr=2,mmon=nmon*nyr,
     &          mtot=nlong*nlat*nmon*nyr)
      dimension temp(nlong,nlat,nmon,nyr)
      dimension tempstd(nlong,nlat,nmon,nyr)
      dimension precip(nlong,nlat,nmon,nyr)
      dimension fill(nlat,nmon,nyr)
      dimension rmb(nlong,nlat)
      dimension snow(nlong,nlat)
      dimension rain(nlong,nlat)
      dimension tmean(nlong,nlat),nmean(nlong,nlat)
      dimension pmean(nlong,nlat)
      dimension pdd(nlong,nlat)
      dimension ampli(nlong,nlat),ttmon(mmon)
      dimension mday(12)
      dimension topo(nlong,nlat)
      data mday /31,28,31,30,31,30,31,31,30,31,30,31/
      data nmean /ntot*0/
      data tmean /ntot*0.0/
      data pmean /ntot*0.0/
      data snow /ntot*0.0/
      data rain /ntot*0.0/
      data pdd /ntot*0.0/
      data temp /mtot*0.0/
      data zero /273.16/
      data ipass /0/
      data fudge /0./
      data factor /-1./
      save topo,temp,precip,tempstd,ipass
c
      rlat=rlat1
      rlong=rlong1
      elev=elev1*0.001
      toffset=tnsl
      call ll2ji(rlong,rlat,kklong,kklat)
      n4=4
      n12=12
      n2=2
      if(ipass.eq.0) then
        print *,'reading climate data'
        ipass=1
        open(49,file='../ext-data/topo.bindmp',
     &       form='unformatted')
        read(49) topo
        do i=1,nlat
          do j=1,nlong
            if(topo(j,i).lt.0.0) topo(j,i)=0.0
          enddo
        enddo
        close(49)
        open(49,file='../ext-data/temperature.bindmp',
     &       form='unformatted')
        read(49) temp
        do iyr=1,n2
          do imon=1,n12
            do i=1,nlat
              fill(i,imon,iyr)=0.0
              nfill=0
              do j=1,nlong
                if(temp(j,i,imon,iyr).ne.0.0) then
                  fill(i,imon,iyr)=fill(i,imon,iyr)+
     &                             temp(j,i,imon,iyr)-
     &                             rlapse*topo(j,i)
                  nfill=nfill+1
                endif
              enddo
              if(nfill.ne.0) then
                fill(i,imon,iyr)=fill(i,imon,iyr)/nfill
              else
                fill(i,imon,iyr)=0
              endif
            enddo
          enddo
        enddo
        do iyr=1,n2
          do imon=1,n12
            do i=1,nlat
              do j=1,nlong
                if(temp(j,i,imon,iyr).eq.0.0) then
                  temp(j,i,imon,iyr)=fill(i,imon,iyr)
                endif
              enddo
            enddo
          enddo
        enddo
        close(49)
        open(49,file='../ext-data/temp2mstd.bindmp',
     &       form='unformatted')
        read(49) tempstd
        do iyr=1,n2
          do imon=1,n12
            do i=1,nlat
              fill(i,imon,iyr)=0.0
              nfill=0
              do j=1,nlong
                if(tempstd(j,i,imon,iyr).ne.0.0) then
                  fill(i,imon,iyr)=fill(i,imon,iyr)+
     &                             tempstd(j,i,imon,iyr)
                  nfill=nfill+1
                endif
              enddo
              if(nfill.ne.0) then
                fill(i,imon,iyr)=fill(i,imon,iyr)/nfill
              else
                fill(i,imon,iyr)=0
              endif
            enddo
          enddo
        enddo
        do iyr=1,n2
          do imon=1,n12
            do i=1,nlat
              do j=1,nlong
                if(tempstd(j,i,imon,iyr).eq.0.0) then
                  tempstd(j,i,imon,iyr)=fill(i,imon,iyr)
                endif
              enddo
            enddo
          enddo
        enddo
        close(49)
        open(49,file='../ext-data/precip.bindmp',
     &       form='unformatted')
        read(49) precip
        do iyr=1,n2
          do imon=1,n12
            do i=1,nlat
              fill(i,imon,iyr)=0.0
              nfill=0
              do j=1,nlong
                if(precip(j,i,imon,iyr).ne.0.0) then
                  fill(i,imon,iyr)=fill(i,imon,iyr)+
     &                             precip(j,i,imon,iyr)
                  nfill=nfill+1
                endif
              enddo
              if(nfill.ne.0) then
                fill(i,imon,iyr)=fill(i,imon,iyr)/nfill
              else
                fill(i,imon,iyr)=0
              endif
            enddo
          enddo
        enddo
        do iyr=1,n2
          do imon=1,n12
            do i=1,nlat
              do j=1,nlong
                if(precip(j,i,imon,iyr).eq.0.0) then
                  precip(j,i,imon,iyr)=fill(i,imon,iyr)
                endif
              enddo
            enddo
          enddo
        enddo
        close(49)
      endif
      do i=kklat,kklat
        do j=kklong,kklong
          tmean(j,i)=0.0
          pmean(j,i)=0.0
          snow(j,i)=0.0
          rain(j,i)=0.0
          pdd(j,i)=0.0
          nmean(j,i)=0
        enddo
      enddo
      do iyr=1,n2
        do imon=1,n12
          do i=kklat,kklat
            do j=kklong,kklong
              if(precip(j,i,imon,iyr).ne.0.0) then
                pmean(j,i)=pmean(j,i)+precip(j,i,imon,iyr)
                nmean(j,i)=nmean(j,i)+1
              endif
            enddo
          enddo
        enddo
      enddo
      do i=kklat,kklat
        do j=kklong,kklong
          if(pmean(j,i).ne.0) then
            pmean(j,i)=0.001*pmean(j,i)/n2
          else
            pmean(j,i)=-999.
          endif
        enddo
      enddo
c
      do i=kklat,kklat
        do j=kklong,kklong
          pdd(j,i)=0.0
          snow(j,i)=0.0
          rain(j,i)=0.0
          tmean(j,i)=0.0
          nmean(j,i)=0
        enddo
      enddo
       do iyr=1,n2
        do imon=1,n12
          do i=kklat,kklat
            do j=kklong,kklong
              if(temp(j,i,imon,iyr).ne.0.0) then
c ... use the following to REDUCE preciptation as climate cools ...
                WRM=WARMING(TNSL+rlapse*(elev-topo(j,i)))
c               WRM=1.0
                temploc=temp(j,i,imon,iyr)+rlapse*(elev-topo(j,i))+
     &                  toffset
                fudge=tempstd(j,i,imon,iyr)*factor
                tmean(j,i)=tmean(j,i)+temploc
                nmean(j,i)=nmean(j,i)+1
                if(temploc.gt.(zero+fudge)) then
                  pdd(j,i)=pdd(j,i)+
     &            (temploc-(zero+fudge))*mday(imon)
                  rain(j,i)=rain(j,i)+
     &                      WRM*precip(j,i,imon,iyr)*0.001/n2
                else
                  snow(j,i)=snow(j,i)+
     &                      WRM*precip(j,i,imon,iyr)*0.001/n2
                endif
              endif
            enddo
          enddo
        enddo
      enddo
c
      if(.false.) then
        do i=kklat,kklat
          do j=kklong,kklong
            ampli(j,i)=-999.
          enddo
        enddo
        do i=kklat,kklat
          do j=kklong,kklong
            jmon=0
            do iyr=1,n2
              do imon=1,n12
                if(temp(j,i,imon,iyr).ne.0.0) then
                  temploc=temp(j,i,imon,iyr)+rlapse*(elev-topo(j,i))
                  jmon=jmon+1
                  ttmon(jmon)=temploc
                  endif
              enddo
            enddo
            if(jmon.gt.1) then
              call fourier(jmon,ttmon,cccmax,ncmax)
              ampli(j,i)=cccmax
            else
              ampli(j,i)=-999.
            endif
          enddo
        enddo
      endif
c
      do i=kklat,kklat
        do j=kklong,kklong
          if(nmean(j,i).ne.0) then
            tmean(j,i)=tmean(j,i)/nmean(j,i)
          else
            tmean(j,i)=-999.
          endif
          if(pdd(j,i).ne.0.0) then
            pdd(j,i)=pdd(j,i)*0.6e-3/n2
          else
            pdd(j,i)=-999.
          endif
          if(pmean(j,i).eq.-999.) then
            rmb(j,i)=-999.
          elseif(pdd(j,i).eq.-999.) then
            rmb(j,i)=snow(j,i)
          else
            rmb(j,i)=snow(j,i)-pdd(j,i)
          endif
          if(rain(j,i).eq.0) rain(j,i)=-999.
          if(snow(j,i).eq.0) snow(j,i)=-999.
        enddo
      enddo
      if(tmean(kklong,kklat).ne.-999.) then
        tsurf=tmean(kklong,kklat)-zero
      else
        tsurf=-999.
      endif
      acc=rmb(kklong,kklat)
      abl=pdd(kklong,kklat)
      if(tsurf.ne.-999. .and. .false.) then
        print *,' topo         (km)  ',topo(kklong,kklat),elev
        print *,' mean temp    (degC)',tmean(kklong,kklat)-zero
        print *,' mean precip  (m)   ',pmean(kklong,kklat)
        print *,' rain         (m)   ',-rain(kklong,kklat),
     &          100*rain(kklong,kklat)/pmean(kklong,kklat)
        print *,' snow         (m)   ',snow(kklong,kklat),
     &          100*snow(kklong,kklat)/pmean(kklong,kklat)
        print *,' pdd-ablation (m)   ',-pdd(kklong,kklat)
        print *,' mass balance (m)   ',rmb(kklong,kklat)
        print *,' seasonal amplitude ',ampli(kklong,kklat)
      endif
      end
C===========================================
      subroutine fourier(npts,f,cccmax,ncmax)
      parameter(m=4)
c     parameter(pi=4.*atan(1.))
      data pi /3.1415927/
      dimension f(npts),aaa(0:m),bbb(m),ccc(m),phi(m)
      rl=real(npts/2)
      sum=0.
      do i=1,npts
        sum=sum+f(i)
      enddo
      aaa(0)=sum/rl
c      print *,0,aaa(0)/2
      aaamax=-1e30
      bbbmax=-1e30
      namax=0
      nbmax=0
      fmin=1e30
      fmax=-fmin
      do i=1,npts
        fmin=min(fmin,f(i))
        fmax=max(fmax,f(i))
      enddo
      do j=1,m
        suma=0.0
        sumb=0.0
        do i=1,npts
          suma=suma+f(i)*cos(j*pi*i/rl)

          sumb=sumb+f(i)*sin(j*pi*i/rl)
        enddo
        aaa(j)=suma/rl
        bbb(j)=sumb/rl
        if(abs(aaa(j)).gt.aaamax) then
          aaamax=abs(aaa(j))
          namax=j
        endif
        if(abs(bbb(j)).gt.bbbmax) then
          bbbmax=abs(bbb(j))
          nbmax=j
        endif
c        print *,j,aaa(j),bbb(j)
      enddo
      cccmax=-1e30
      ncmax=0
      do n=1,m
        ccc(n)=sqrt(aaa(n)**2+bbb(n)**2)
        if(aaa(n).lt.0) ccc(n)=-ccc(n)
        if(abs(ccc(n)).gt.cccmax) then
          cccmax=abs(ccc(n))
          ncmax=n
        endif
        phi(n)=atan(aaa(n)/bbb(n))
      enddo
      if(.false.) then
        call grstrt(600,600)
        call window(0.,real(24+1),fmin,fmax)
        call linclr(1)
        call move(1.,f(1))
        do jj=1,npts
          call draw(real(jj),f(jj))
          call point(real(jj),f(jj))
        enddo
        call linclr(2)
        do jj=1,npts
          sum=aaa(0)/2.
            do n=namax,namax
              sum=sum+aaa(n)*cos(n*pi*jj/rl)
            enddo
            do n=nbmax,nbmax
              sum=sum+bbb(n)*sin(n*pi*jj/rl)
            enddo
            call point(real(jj),sum)
        enddo
        if(.false.) then
          call newpag
        endif
        call grstop1
      endif
      end
c===================================      
      subroutine ji2ll(j,i,rlong,rlat)
      rlong=j-180-0.5
      if(rlong.lt.0) rlong=rlong+360.
      rlat=90-i+0.5
      end
c=========================================
      subroutine ll2ji(rlong,rlat,j,i)
      j=nint(rlong+180.5)
      if(j.gt.360) j=j-360
      i=nint(89.5-rlat)
      if(j.le.0) then
        j=1
      elseif(j.gt.360) then
        j=360
      endif
      if(i.le.0) then
        i=1
      elseif(i.gt.80) then
        j=80
      endif
      end
c=============================================
      subroutine interp1(x,y,v,xtest,ytest,val)
      implicit real*8(u)
      dimension x(4),y(4),v(4)
c this is machine generated code, from macsyma. dont change anything....
c (it's about 20 times faster than inverting the matrix.................
c-----------------------------------------------------------------------
      U0=X(1)
      U1=Y(1)
      U2=U0*U1
      U3=X(2)
      U4=-(U1*U3)
      U5=U4+U2
      U6=Y(2)
      U7=X(3)
      U8=U5*U6*U7
      U9=-(U0*U1*U3)
      U10=U0*U3*U6
      U11=U1*U3
      U12=-(U0*U6)+U11
      U13=Y(3)
      U14=(U12*U7+U10+U9)*U13
      U15=-(U0*U1)
      U16=U11+U15
      U17=-(U3*U6)
      U18=U6-U1
      U19=(U18*U7+U17+U2)*U13+U16*U6
      U20=X(4)
      U21=U0*U1*U3
      U22=-(U0*U3*U6)
      U23=U3*U6+U15
      U24=U23*U7
      U25=-U3+U0
      U26=U25*U7*U13
      U27=U0*U6
      U28=-U6+U1
      U29=U3-U0
      U30=U29*U13+U28*U7+U27+U4
      U31=Y(4)
      U32=V(4)
      U33=V(3)
      U34=V(2)
      U35=U0*U1*U34
      U36=V(1)
      U37=-(U36*U3*U6)
      U38=-(U1*U34)
      U39=U36*U6
      U40=(U39+U38)*U7
      U41=-(U0*U1*U34)
      U42=U36*U3*U6
      U43=U0*U34
      U44=-(U36*U3)
      U45=U44+U43
      U46=U1*U34
      U47=-(U36*U6)
      U48=(U47+U46)*U7
      U49=-(U0*U34)
      U50=U36*U3
      U51=U50+U49
      U52=-U34+U36
      U53=U34-U36
      VAL=((((U53*U7+U25*U33+U50+U49)*U31+(U52*U13+U18*U33+U47+U46)*U20+
     . U30*U32+U45*U13+U40+U12*U33)*XTEST+(U52*U7+U29*U33+U44+U43)*U20*
     . U31+(U53*U7*U13+(U17+U2)*U33+U42+U41)*U20+(U26+U24+U22+U21)*U32+
     . U51*U7*U13+(U37+U35)*U7+(U10+U9)*U33)*YTEST+(((U53*U13+U28*U33+
     . U39+U38)*U20+U52*U7*U13+U23*U33+U37+U35)*U31+U19*U32+(U48+U42+U41
     . )*U13+U5*U6*U33)*XTEST+((U51*U13+U48+(U27+U4)*U33)*U20+U45*U7*U13
     . +(U42+U41)*U7+(U22+U21)*U33)*U31+((U40+U37+U35)*U13+U16*U6*U33)*
     . U20+(U14+U8)*U32)/((U30*U20+U26+U24+U22+U21)*U31+U19*U20+U14+U8)
c-----------------------------------------------------------------------
      end

c=============================================
      subroutine interp0(xxx,yyy,value,xtest,ytest,val)
      dimension rmat(4,4),rhs(4)
      dimension xxx(4),yyy(4),value(4)
      do i=1,4
        rhs(i)=value(i)
      enddo
      do i=1,4
        rmat(i,1)=1.0            
        rmat(i,2)=xxx(i)            
        rmat(i,3)=yyy(i)            
        rmat(i,4)=xxx(i)*yyy(i)
      enddo  

      call GAUSSJ(rmat,4,4,rhs,1,1)

      val=rhs(1)+rhs(2)*xtest+
     &         rhs(3)*ytest+rhs(4)*xtest*ytest
      end
c=============================================

      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      PARAMETER (NMAX=50)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) PAUSE 'Singular matrix.'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
c==============================================================
      subroutine lookup26(rlapse,rlat1,rlong1,elev1,psurf1,
     &                    tnsl,tsurf,acc,abl)
      real*8 rlapse,rlat1,rlong1,elev1,psurf1,
     &       tnsl,tsurf,acc,abl,WARMING
      parameter(nlong=192,nlat=94,nmon=12)
      dimension rlong(nlong), rlat(nlat)
      dimension pmean(nlong,nlat,nmon)
      dimension ptot(nlong,nlat)
      dimension snow(nlong,nlat)
      dimension rain(nlong,nlat)
      dimension pdd(nlong,nlat)
      dimension rmb(nlong,nlat)
      dimension tmean(nlong,nlat,nmon)
      dimension tannual(nlong,nlat)
      dimension nannual(nlong,nlat)
      dimension topo(nlong,nlat)
      dimension mday(12),dist(4),value(4),xxx(4),yyy(4)
      data mday /31,28,31,30,31,30,31,31,30,31,30,31/
      character*8 cmon(12)
      data cmon /'january','february','march','april','may','june',
     &           'july','august','sepember','october','november',
     &           'december'/
      data ipass /0/
      data zero /0.0/
      logical quiet,close
      parameter(close=.false.,quiet=.true.)
      data fudge /-5./
      data reject /-999./
      data correct /0./
      save rlong,rlat,pmean,tmean,topo,ipass
c ... read rlong,rlat,pmean,tmean topo .....................
c      fudge=factor
      factor=fudge
      if(ipass.eq.0) then
        print *,' reading ncep2 climate data '
        open(12,file='../ext-data/lookup26.bin',form='unformatted')
        read(12) nlongt,nlatt,nmont
        if(nlongt.ne.nlong.or.nlatt.ne.nlat.or.nmont.ne.nmon) then
          print *,' error, stopping '
          stop
        endif
        read(12) rlong,rlat,pmean,tmean,topo
        close(12)
        ipass=1
      endif
      elev=real((elev1-psurf1)*0.001)
      toffset=real(tnsl)+correct
      if(close) then
        call closest(nlong,nlat,rlong,rlat,
     &             real(rlong1),real(rlat1),kklong,kklat)
      else
        call bracket(nlong,nlat,rlong,rlat,
     &             real(rlong1),real(rlat1),kklong,kklat,dist)
        xxx(1)=rlong(kklong  )
        xxx(2)=rlong(kklong+1)
        xxx(3)=rlong(kklong+1)
        xxx(4)=rlong(kklong  )
        yyy(1)=rlat(kklat  )
        yyy(2)=rlat(kklat  )
        yyy(3)=rlat(kklat+1)
        yyy(4)=rlat(kklat+1)
      endif
c      print *,real(rlong1),real(rlat1),rlong(kklong),rlat(kklat)
      if(quiet) then
        lat1=kklat
        long1=kklong
        if(close) then
          lat2=lat1
          long2=long1
        else
          lat2=lat1+1
          long2=long1+1
        endif
      else
        lat1=1
        lat2=nlat
        long1=1
        long2=nlong
      endif
c ... zero out arrays .................................
      do i=lat1,lat2
        do j=long1,long2
          ptot(j,i)=0.0
          pdd(j,i)=0.0
          snow(j,i)=0.0
          rain(j,i)=0.0
          rmb(j,i)=0.0
          tannual(j,i)=0.0
          nannual(j,i)=0
        enddo
      enddo
c ... calculate tannual, mean annual temperature (screen10##)
      do imon=1,nmon
        do i=lat1,lat2
          do j=long1,long2
            tempji=tmean(j,i,imon)+toffset
            tannual(j,i)=tannual(j,i)+tempji+
     &                   rlapse*elev
            nannual(j,i)=nannual(j,i)+1
          enddo
        enddo
      enddo
      do i=lat1,lat2
        do j=long1,long2
          tannual(j,i)=tannual(j,i)/nannual(j,i)
        enddo
      enddo
c ... calculate pdd .......................................
      do imon=1,nmon
        do i=lat1,lat2
          do j=long1,long2
c           if(.not.quiet) elev=topo(j,i)
            WRM=REDUCE(real(TNSL+rlapse*elev))
c           WRM=REDUCE2(real(TNSL),real(rlapse*elev))
c           WRM=1.
            temploc=tmean(j,i,imon)+rlapse*elev+toffset
            if(temploc.gt.(zero+fudge)) then
              pdd(j,i)=pdd(j,i)+
     &            (temploc-(zero+fudge))*mday(imon)
              rain(j,i)=rain(j,i)+
     &                  WRM*pmean(j,i,imon)*1
            else
              snow(j,i)=snow(j,i)+
     &                  WRM*pmean(j,i,imon)*1
            endif
          enddo
        enddo
      enddo
      do i=lat1,lat2
        do j=long1,long2
          if(pdd(j,i).ne.0.0) then
            pdd(j,i)=pdd(j,i)*0.6e-3
          else
            pdd(j,i)=reject
          endif
          if(pdd(j,i).eq.reject) then
            rmb(j,i)=snow(j,i)
          else
            rmb(j,i)=snow(j,i)-pdd(j,i)
          endif
        enddo
      enddo
      if(tannual(kklong,kklat).ne.reject) then
        if(close) then
          tsurf=tannual(kklong,kklat)-zero
        else
          value(1)=tannual(kklong  ,kklat  )
          value(2)=tannual(kklong+1,kklat  )
          value(3)=tannual(kklong+1,kklat+1)
          value(4)=tannual(kklong  ,kklat+1)
          tsurf=wmean(dist,value)
          call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
          tsurf=val
        endif
      else
        tsurf=reject
      endif
      if(close) then
        acc=rmb(kklong,kklat)
      else
        value(1)=rmb(kklong  ,kklat  )
        value(2)=rmb(kklong+1,kklat  )
        value(3)=rmb(kklong+1,kklat+1)
        value(4)=rmb(kklong  ,kklat+1)
c        acc=wmean(dist,value)
        call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
        acc=val
c ... EXPERIMENTAL arbitrary negative offset 0.1
c       acc=acc-0.1
      endif
      if(.false.) then
        print *,'tannual:',tannual(kklong,kklat)
        do i=1,12
          print *,'tmon:',i,tmean(kklong,kklat,i),pmean(kklong,kklat,i)
        enddo
        print *,'rm,sn,rn,pd    ',rmb(kklong,kklat),
     &                      snow(kklong,kklat),
     &                      rain(kklong,kklat),
     &                      pdd(kklong,kklat)
      endif
      if(close) then
        abl=pdd(kklong,kklat)
      else
        value(1)=pdd(kklong  ,kklat  )
        value(2)=pdd(kklong+1,kklat  )
        value(3)=pdd(kklong+1,kklat+1)
        value(4)=pdd(kklong  ,kklat+1)
        abl=wmean(dist,value)
        call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
        abl=val
      endif
      end
c=========================================
      subroutine lookup27(rlapse,rlat1,rlong1,elev1,psurf1,
     &                    tnsl,tsurf,acc,abl)
      real*8 rlapse,rlat1,rlong1,elev1,psurf1,
     &       tnsl,tsurf,acc,abl,WARMING
      parameter(nlong=128,nlat=64,nmon=12)
      dimension rlong(nlong), rlat(nlat)
      dimension pmean(nlong,nlat,nmon)
      dimension ptot(nlong,nlat)
      dimension snow(nlong,nlat)
      dimension rain(nlong,nlat)
      dimension pdd(nlong,nlat)
      dimension rmb(nlong,nlat)
      dimension tmean(nlong,nlat,nmon)
      dimension tannual(nlong,nlat)
      dimension nannual(nlong,nlat)
      dimension topo(nlong,nlat)
      dimension mday(12),dist(4),value(4),xxx(4),yyy(4)
      data mday /31,28,31,30,31,30,31,31,30,31,30,31/
      character*8 cmon(12)
      data cmon /'january','february','march','april','may','june',
     &           'july','august','sepember','october','november',
     &           'december'/
      data ipass /0/
      data zero /0.0/
      logical quiet,close
      parameter(close=.false.,quiet=.true.)
      data fudge /-5./
      data reject /-999./
      data correct /14./
      save rlong,rlat,pmean,tmean,topo,ipass
c ... read rlong,rlat,pmean,tmean topo .....................
c      fudge=factor
      factor=fudge
      if(ipass.eq.0) then
        print *,' reading ncep2 climate data '
        open(12,file='../ext-data/lookup27.bin',form='unformatted')
        read(12) nlongt,nlatt,nmont
        if(nlongt.ne.nlong.or.nlatt.ne.nlat.or.nmont.ne.nmon) then
          print *,' error, stopping '
          stop
        endif
        read(12) rlong,rlat,pmean,tmean,topo
        close(12)
        ipass=1
      endif
      elev=real((elev1-psurf1)*0.001)
      toffset=real(tnsl)+correct
      if(close) then
        call closest(nlong,nlat,rlong,rlat,
     &             real(rlong1),real(rlat1),kklong,kklat)
      else
        call bracket(nlong,nlat,rlong,rlat,
     &             real(rlong1),real(rlat1),kklong,kklat,dist)
        xxx(1)=rlong(kklong  )
        xxx(2)=rlong(kklong+1)
        xxx(3)=rlong(kklong+1)
        xxx(4)=rlong(kklong  )
        yyy(1)=rlat(kklat  )
        yyy(2)=rlat(kklat  )
        yyy(3)=rlat(kklat+1)
        yyy(4)=rlat(kklat+1)
      endif
c      print *,real(rlong1),real(rlat1),rlong(kklong),rlat(kklat)
      if(quiet) then
        lat1=kklat
        long1=kklong
        if(close) then
          lat2=lat1
          long2=long1
        else
          lat2=lat1+1
          long2=long1+1
        endif
      else
        lat1=1
        lat2=nlat
        long1=1
        long2=nlong
      endif
c ... zero out arrays .................................
      do i=lat1,lat2
        do j=long1,long2
          ptot(j,i)=0.0
          pdd(j,i)=0.0
          snow(j,i)=0.0
          rain(j,i)=0.0
          rmb(j,i)=0.0
          tannual(j,i)=0.0
          nannual(j,i)=0
        enddo
      enddo
c ... calculate tannual, mean annual temperature (screen10##)
      do imon=1,nmon
        do i=lat1,lat2
          do j=long1,long2
            tempji=tmean(j,i,imon)
            tannual(j,i)=tannual(j,i)+tempji+toffset+
     &                   rlapse*elev
            nannual(j,i)=nannual(j,i)+1
          enddo
        enddo
      enddo
      do i=lat1,lat2
        do j=long1,long2
          tannual(j,i)=tannual(j,i)/nannual(j,i)
        enddo
      enddo
c ... calculate pdd .......................................
      do imon=1,nmon
        do i=lat1,lat2
          do j=long1,long2
c           if(.not.quiet) elev=topo(j,i)
            WRM=REDUCE(real(TNSL+rlapse*elev))
c           WRM=REDUCE2(real(TNSL),real(rlapse*elev))
c           WRM=1.
            temploc=tmean(j,i,imon)+rlapse*elev+toffset
            if(temploc.gt.(zero+fudge)) then
              pdd(j,i)=pdd(j,i)+
     &            (temploc-(zero+fudge))*mday(imon)
              rain(j,i)=rain(j,i)+
     &                  WRM*pmean(j,i,imon)*1
            else
              snow(j,i)=snow(j,i)+
     &                  WRM*pmean(j,i,imon)*1
            endif
          enddo
        enddo
      enddo
      do i=lat1,lat2
        do j=long1,long2
          if(pdd(j,i).ne.0.0) then
            pdd(j,i)=pdd(j,i)*0.6e-3
          else
            pdd(j,i)=reject
          endif
          if(pdd(j,i).eq.reject) then
            rmb(j,i)=snow(j,i)
          else
            rmb(j,i)=snow(j,i)-pdd(j,i)
          endif
        enddo
      enddo
      if(tannual(kklong,kklat).ne.reject) then
        if(close) then
          tsurf=tannual(kklong,kklat)-zero
        else
          value(1)=tannual(kklong  ,kklat  )
          value(2)=tannual(kklong+1,kklat  )
          value(3)=tannual(kklong+1,kklat+1)
          value(4)=tannual(kklong  ,kklat+1)
          tsurf=wmean(dist,value)
          call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
          tsurf=val
        endif
      else
        tsurf=reject
      endif
      if(close) then
        acc=rmb(kklong,kklat)
      else
        value(1)=rmb(kklong  ,kklat  )
        value(2)=rmb(kklong+1,kklat  )
        value(3)=rmb(kklong+1,kklat+1)
        value(4)=rmb(kklong  ,kklat+1)
        acc=wmean(dist,value)
        call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
        acc=val
c ... EXPERIMENTAL arbitrary negative offset 0.1
c        acc=acc-0.1
      endif
      if(.false.) then
        print *,'tannual:',tannual(kklong,kklat)
        do i=1,12
          print *,'tmon:',i,tmean(kklong,kklat,i),pmean(kklong,kklat,i)
        enddo
        print *,'rm,sn,rn,pd    ',rmb(kklong,kklat),
     &                      snow(kklong,kklat),
     &                      rain(kklong,kklat),
     &                      pdd(kklong,kklat)
      endif
      if(close) then
        abl=pdd(kklong,kklat)
      else
        value(1)=pdd(kklong  ,kklat  )
        value(2)=pdd(kklong+1,kklat  )
        value(3)=pdd(kklong+1,kklat+1)
        value(4)=pdd(kklong  ,kklat+1)
        abl=wmean(dist,value)
        call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
        abl=val
      endif
      end
c===========================================================
      subroutine lookup29(rlapse,rlat1,rlong1,elev1,psurf1,
     &                    tnsl,tsurf,acc,abl)
      real*8 rlapse,rlat1,rlong1,elev1,psurf1,
     &       tnsl,tsurf,acc,abl,WARMING
      parameter(nlong=192,nlat=94,nmon=12)
      dimension rlong(nlong), rlat(nlat)
      dimension pmean(nlong,nlat,nmon)
      dimension ptot(nlong,nlat)
      dimension snow(nlong,nlat)
      dimension rain(nlong,nlat)
      dimension pdd(nlong,nlat)
      dimension rmb(nlong,nlat)
      dimension tmean(nlong,nlat,nmon)
      dimension tannual(nlong,nlat)
      dimension nannual(nlong,nlat)
      dimension topo(nlong,nlat)
      dimension mday(12),dist(4),value(4),xxx(4),yyy(4)
      data mday /31,28,31,30,31,30,31,31,30,31,30,31/
      character*8 cmon(12)
      data cmon /'january','february','march','april','may','june',
     &           'july','august','sepember','october','november',
     &           'december'/
      data ipass /0/
      data zero /0.0/
      logical quiet,close
      parameter(close=.false.,quiet=.true.)
      data fudge /-5./
      data reject /-999./
      data correct /0./
      save rlong,rlat,pmean,tmean,topo,ipass
c ... read rlong,rlat,pmean,tmean topo .....................
      if(ipass.eq.0) then
        print *,' reading ncep2 climate data '
        open(12,file='../ext-data/lookup26.bin',form='unformatted')
        read(12) nlongt,nlatt,nmont
        if(nlongt.ne.nlong.or.nlatt.ne.nlat.or.nmont.ne.nmon) then
          print *,' error, stopping '
          stop
        endif
        read(12) rlong,rlat,pmean,tmean,topo
        close(12)
        ipass=1
      endif
      elev=real((elev1-psurf1)*0.001)
      toffset=real(tnsl)+correct
      if(close) then
        call closest(nlong,nlat,rlong,rlat,
     &             real(rlong1),real(rlat1),kklong,kklat)
      else
        call bracket(nlong,nlat,rlong,rlat,
     &             real(rlong1),real(rlat1),kklong,kklat,dist)
        xxx(1)=rlong(kklong  )
        xxx(2)=rlong(kklong+1)
        xxx(3)=rlong(kklong+1)
        xxx(4)=rlong(kklong  )
        yyy(1)=rlat(kklat  )
        yyy(2)=rlat(kklat  )
        yyy(3)=rlat(kklat+1)
        yyy(4)=rlat(kklat+1)
      endif
c      print *,real(rlong1),real(rlat1),rlong(kklong),rlat(kklat)
      if(quiet) then
        lat1=kklat
        long1=kklong
        if(close) then
          lat2=lat1
          long2=long1
        else
          lat2=lat1+1
          long2=long1+1
        endif
      else
        lat1=1
        lat2=nlat
        long1=1
        long2=nlong
      endif
c ... zero out arrays .................................
      do i=lat1,lat2
        do j=long1,long2
          ptot(j,i)=0.0
          pdd(j,i)=0.0
          snow(j,i)=0.0
          rain(j,i)=0.0
          rmb(j,i)=0.0
          tannual(j,i)=0.0
          nannual(j,i)=0
        enddo
      enddo
c ... calculate tannual, mean annual temperature (screen10##)
      do imon=1,nmon
        do i=lat1,lat2
          do j=long1,long2
            ediff=elev+(psurf1*0.001-max(topo(j,i),0.0))*1
            tempji=tmean(j,i,imon)+toffset
            tannual(j,i)=tannual(j,i)+tempji+
     &                   rlapse*ediff
            nannual(j,i)=nannual(j,i)+1
          enddo
        enddo
      enddo
      do i=lat1,lat2
        do j=long1,long2
          tannual(j,i)=tannual(j,i)/nannual(j,i)
        enddo
      enddo
c ... calculate pdd .......................................
      do imon=1,nmon
        do i=lat1,lat2
          do j=long1,long2
            ediff=elev+(psurf1*0.001-max(topo(j,i),0.0))*1
c           WRM=WARMING(dble(TNSL+rlapse*ediff))
c           if(.not.quiet) elev=topo(j,i)
            WRM=REDUCE(real(TNSL+rlapse*ediff))
c           WRM=REDUCE2(real(TNSL),real(rlapse*ediff))
c           WRM=1.
            temploc=tmean(j,i,imon)+rlapse*ediff+toffset
            if(temploc.gt.(zero+fudge)) then
              pdd(j,i)=pdd(j,i)+
     &            (temploc-(zero+fudge))*mday(imon)
              rain(j,i)=rain(j,i)+
     &                  WRM*pmean(j,i,imon)*1
            else
              snow(j,i)=snow(j,i)+
     &                  WRM*pmean(j,i,imon)*1
            endif
          enddo
        enddo
      enddo
      do i=lat1,lat2
        do j=long1,long2
          if(pdd(j,i).ne.0.0) then
            pdd(j,i)=pdd(j,i)*0.6e-3
          else
            pdd(j,i)=reject
          endif
          if(pdd(j,i).eq.reject) then
            rmb(j,i)=snow(j,i)
          else
            rmb(j,i)=snow(j,i)-pdd(j,i)
          endif
        enddo
      enddo
      if(tannual(kklong,kklat).ne.reject) then
        if(close) then
          tsurf=tannual(kklong,kklat)-zero
        else
          value(1)=tannual(kklong  ,kklat  )
          value(2)=tannual(kklong+1,kklat  )
          value(3)=tannual(kklong+1,kklat+1)
          value(4)=tannual(kklong  ,kklat+1)
          tsurf=wmean(dist,value)
          call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
          tsurf=val
        endif
      else
        tsurf=reject
      endif
      if(close) then
        acc=rmb(kklong,kklat)
      else
        value(1)=rmb(kklong  ,kklat  )
        value(2)=rmb(kklong+1,kklat  )
        value(3)=rmb(kklong+1,kklat+1)
        value(4)=rmb(kklong  ,kklat+1)
c        acc=wmean(dist,value)
        call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
        acc=val
c ... EXPERIMENTAL arbitrary negative offset 0.1
c       acc=acc-0.1
      endif
      if(.false.) then
        print *,'tannual:',tannual(kklong,kklat)
        do i=1,12
          print *,'tmon:',i,tmean(kklong,kklat,i),pmean(kklong,kklat,i)
        enddo
        print *,'rm,sn,rn,pd    ',rmb(kklong,kklat),
     &                      snow(kklong,kklat),
     &                      rain(kklong,kklat),
     &                      pdd(kklong,kklat)
      endif
      if(close) then
        abl=pdd(kklong,kklat)
      else
        value(1)=pdd(kklong  ,kklat  )
        value(2)=pdd(kklong+1,kklat  )
        value(3)=pdd(kklong+1,kklat+1)
        value(4)=pdd(kklong  ,kklat+1)
        abl=wmean(dist,value)
        call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
        abl=val
      endif
      end
c=========================================
      subroutine lookup30(rlapse,rlat1,rlong1,elev1,psurf1,
     &                    tnsl,tsurf,acc,abl)
      real*8 rlapse,rlat1,rlong1,elev1,psurf1,
     &       tnsl,tsurf,acc,abl,WARMING
      parameter(nlong=128,nlat=64,nmon=12)
      dimension rlong(nlong), rlat(nlat)
      dimension pmean(nlong,nlat,nmon)
      dimension ptot(nlong,nlat)
      dimension snow(nlong,nlat)
      dimension rain(nlong,nlat)
      dimension pdd(nlong,nlat)
      dimension rmb(nlong,nlat)
      dimension tmean(nlong,nlat,nmon)
      dimension tannual(nlong,nlat)
      dimension nannual(nlong,nlat)
      dimension topo(nlong,nlat)
      dimension mday(12),dist(4),value(4),xxx(4),yyy(4)
      data mday /31,28,31,30,31,30,31,31,30,31,30,31/
      character*8 cmon(12)
      data cmon /'january','february','march','april','may','june',
     &           'july','august','sepember','october','november',
     &           'december'/
      data ipass /0/
      data zero /0.0/
      logical quiet,close
      parameter(close=.false.,quiet=.true.)
      data fudge /-5./
      data reject /-999./
      data correct /14./
      save rlong,rlat,pmean,tmean,topo,ipass
c ... read rlong,rlat,pmean,tmean topo .....................
      if(ipass.eq.0) then
        print *,' reading ncep2 climate data '
        open(12,file='../ext-data/lookup27.bin',form='unformatted')
        read(12) nlongt,nlatt,nmont
        if(nlongt.ne.nlong.or.nlatt.ne.nlat.or.nmont.ne.nmon) then
          print *,' error, stopping '
          stop
        endif
        read(12) rlong,rlat,pmean,tmean,topo
        close(12)
        ipass=1
      endif
      elev=real((elev1-psurf1)*0.001)
      toffset=real(tnsl)+correct
      if(close) then
        call closest(nlong,nlat,rlong,rlat,
     &             real(rlong1),real(rlat1),kklong,kklat)
      else
        call bracket(nlong,nlat,rlong,rlat,
     &             real(rlong1),real(rlat1),kklong,kklat,dist)
        xxx(1)=rlong(kklong  )
        xxx(2)=rlong(kklong+1)
        xxx(3)=rlong(kklong+1)
        xxx(4)=rlong(kklong  )
        yyy(1)=rlat(kklat  )
        yyy(2)=rlat(kklat  )
        yyy(3)=rlat(kklat+1)
        yyy(4)=rlat(kklat+1)
      endif
c      print *,real(rlong1),real(rlat1),rlong(kklong),rlat(kklat)
      if(quiet) then
        lat1=kklat
        long1=kklong
        if(close) then
          lat2=lat1
          long2=long1
        else
          lat2=lat1+1
          long2=long1+1
        endif
      else
        lat1=1
        lat2=nlat
        long1=1
        long2=nlong
      endif
c ... zero out arrays .................................
      do i=lat1,lat2
        do j=long1,long2
          ptot(j,i)=0.0
          pdd(j,i)=0.0
          snow(j,i)=0.0
          rain(j,i)=0.0
          rmb(j,i)=0.0
          tannual(j,i)=0.0
          nannual(j,i)=0
        enddo
      enddo
c ... calculate tannual, mean annual temperature (screen10##)
      do imon=1,nmon
        do i=lat1,lat2
          do j=long1,long2
            ediff=elev+(psurf1*0.001-max(topo(j,i),0.0))*1
c           print*,psurf1*0.001,max(topo(j,i),0.0)
            tempji=tmean(j,i,imon)+toffset
            tannual(j,i)=tannual(j,i)+tempji+
     &                   rlapse*ediff
            nannual(j,i)=nannual(j,i)+1
          enddo
        enddo
      enddo
      do i=lat1,lat2
        do j=long1,long2
          tannual(j,i)=tannual(j,i)/nannual(j,i)
        enddo
      enddo
c ... calculate pdd .......................................
      do imon=1,nmon
        do i=lat1,lat2
          do j=long1,long2
            ediff=elev+(psurf1*0.001-max(topo(j,i),0.0))*1
c           WRM=WARMING(dble(TNSL+rlapse*ediff))
c           if(.not.quiet) elev=topo(j,i)
            WRM=REDUCE(real(TNSL+rlapse*ediff))
c           WRM=REDUCE2(real(TNSL),real(rlapse*ediff))
c           WRM=1.
            temploc=tmean(j,i,imon)+rlapse*ediff+toffset
            if(temploc.gt.(zero+fudge)) then
              pdd(j,i)=pdd(j,i)+
     &            (temploc-(zero+fudge))*mday(imon)
              rain(j,i)=rain(j,i)+
     &                  WRM*pmean(j,i,imon)*1
            else
              snow(j,i)=snow(j,i)+
     &                  WRM*pmean(j,i,imon)*1
            endif
          enddo
        enddo
      enddo
      do i=lat1,lat2
        do j=long1,long2
          if(pdd(j,i).ne.0.0) then
            pdd(j,i)=pdd(j,i)*0.6e-3
          else
            pdd(j,i)=reject
          endif
          if(pdd(j,i).eq.reject) then
            rmb(j,i)=snow(j,i)
          else
            rmb(j,i)=snow(j,i)-pdd(j,i)
          endif
        enddo
      enddo
      if(tannual(kklong,kklat).ne.reject) then
        if(close) then
          tsurf=tannual(kklong,kklat)-zero
        else
          value(1)=tannual(kklong  ,kklat  )
          value(2)=tannual(kklong+1,kklat  )
          value(3)=tannual(kklong+1,kklat+1)
          value(4)=tannual(kklong  ,kklat+1)
          tsurf=wmean(dist,value)
          call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
          tsurf=val
        endif
      else
        tsurf=reject
      endif
      if(close) then
        acc=rmb(kklong,kklat)
      else
        value(1)=rmb(kklong  ,kklat  )
        value(2)=rmb(kklong+1,kklat  )
        value(3)=rmb(kklong+1,kklat+1)
        value(4)=rmb(kklong  ,kklat+1)
        acc=wmean(dist,value)
        call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
        acc=val
c ... EXPERIMENTAL arbitrary negative offset 0.1
c        acc=acc-0.1
      endif
      if(.false.) then
        print *,'tannual:',tannual(kklong,kklat)
        do i=1,12
          print *,'tmon:',i,tmean(kklong,kklat,i),pmean(kklong,kklat,i)
        enddo
        print *,'rm,sn,rn,pd    ',rmb(kklong,kklat),
     &                      snow(kklong,kklat),
     &                      rain(kklong,kklat),
     &                      pdd(kklong,kklat)
      endif
      if(close) then
        abl=pdd(kklong,kklat)
      else
        value(1)=pdd(kklong  ,kklat  )
        value(2)=pdd(kklong+1,kklat  )
        value(3)=pdd(kklong+1,kklat+1)
        value(4)=pdd(kklong  ,kklat+1)
        abl=wmean(dist,value)
        call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
        abl=val
      endif
      end






c=========================================
      real function wmean(dist,value)
      dimension dist(4),value(4)
      do i=1,4
        if(dist(i).eq.0) then
          wmean=value(i)
          return
        endif
      enddo
      total=0.
      wmean=0.0
      do i=1,4
        denom=1./dist(i)
        total=total+denom
        wmean=wmean+denom*value(i)
      enddo
      wmean=wmean/total
      end
c=========================================
      subroutine closest(nlong,nlat,rlong,rlat,
     &             rlong1,rlat1,j,i)
      dimension rlong(nlong),rlat(nlat)
      j=0
      dmin=1e30
      do n=1,nlong
        dist=abs(rlong(n)-rlong1)
        if(dist.lt.dmin) then
          dmin=dist
          j=n
        endif
      enddo
      i=0
      dmin=1e30
      do n=1,nlat
        dist=abs(rlat(n)-rlat1)
        if(dist.lt.dmin) then
          dmin=dist
          i=n
        endif
      enddo
      print *,j,i
      end
c=========================================
      subroutine bracket(nlong,nlat,rlong,rlat,
     &             rlong1,rlat1,j,i,dist)
      dimension rlong(nlong),rlat(nlat),dist(4)
      call finder(rlong,nlong,rlong1,j)
      call finder(rlat,nlat,rlat1,i)
c      print *,j,i
c      print *,rlong(j),rlong1,rlong(j+1)
c      print *,rlat(i),rlat1,rlat(i+1)
      if(j.lt.1) j=1
      if(j.ge.nlong) j=nlong-1
      if(i.lt.1) i=1
      if(i.ge.nlat) i=nlat-1
      dist(1)=(rlong(j  )-rlong1)**2+(rlat(i  )-rlat1)**2
      dist(2)=(rlong(j+1)-rlong1)**2+(rlat(i  )-rlat1)**2
      dist(3)=(rlong(j+1)-rlong1)**2+(rlat(i+1)-rlat1)**2
      dist(4)=(rlong(j  )-rlong1)**2+(rlat(i+1)-rlat1)**2
c      print *,dist
c      pause
      end
c=========================================
      SUBROUTINE finder(XX,N,X,J)
      DIMENSION XX(N)
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
      J=JL
      RETURN
      END
C===============================================
      FUNCTION REDUCE(DT)
      const=0.007D0
c      const=0.03D0
      IF(DT.GT.0.0) THEN
        S=const
      ELSEIF(DT.GT.-10.0) THEN
        S=const-const*DT/10.
      ELSE
        S=const*2
      ENDIF
      REDUCE=(1.0+S)**DT
      END

C===============================================
      FUNCTION REDUCE2(tt,dt)
c
      TF=0.67d0*(tt+273.16d0)+88.9d0
      TERM1=-9.09718d0*(273.16d0/TF-1.0d0)                                    
      TERM2=-3.56654d0*LOG10(273.16d0/TF)                                   
      TERM3=0.876793d0*(1.0d0-TF/273.16d0)+0.785835d0
      EXPON=TERM1+TERM2+TERM3                                           
      ES0=10.d0**EXPON                                                     
c
      TF=0.67d0*(TT+DT+273.16d0)+88.9d0
      TERM1=-9.09718d0*(273.16d0/TF-1.0d0)                                    
      TERM2=-3.56654d0*LOG10(273.16d0/TF)                                   
      TERM3=0.876793d0*(1.0d0-TF/273.16d0)+0.785835d0
      EXPON=TERM1+TERM2+TERM3                                           
      ES=10.d0**EXPON    
c
      REDUCE2=ES/ES0
      end                                                 
      SUBROUTINE WMOVERSIMP(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL, DT,
     &                 ETA, XI, W, WTHICK, HTICE, DEPB, LM,
     &                 BMELT, ALPHAC, TOTALW, TOTALP,IPLOT)
c
c ... a purely local water solver, simple 1/e leakage, with source=bmelt
c
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION BMELT(NMAX),ALPHAC(3)
      DIMENSION LM(5),HTICE(NMAX),DEPB(NMAX)
      DIMENSION KX(NMAX,4),X(NMAX),Y(NMAX),NTYPE(NMAX)
      REAL*8 WTHICK(NMAX)

      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2)
      DIMENSION PSI(4),DPSI(4,2)
      DIMENSION XY(2,4),XI(2,9),ETA(2,9),W(2,9)

      REAL*4 TB(2)
c     REAL*4 TB(2),ETIME,DTIME
c     EXTERNAL ETIME,DTIME
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
123   FORMAT(A25,T30,1PG13.6,G13.6)
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      SAVE ISTART,WSAVE,TSAVE,PSAVE
      DATA ISTART /0/
      AREAWET=0
      AREATOT=0
C ... CALL CHECKER THAT PUTS WATER ON/OFF ICE-FREE NODES
c ... (LAST ARG IS VALUE TO PUT ON ICE-FREE NODES, ZERO IT BEFORE
c     (turn on/off internally)
      CALL CHECKER(NMAX, NUMNP, NUMEL, NTYPE, KX, HTICE, DEPB, 
     &             WTHICK,0.0d0)
      IF(ISTART.EQ.0) THEN
        ISTART=1
        WSAVE=TOTALW
        IF(AREAWET.EQ.0.0) THEN
          TSAVE=0D0
        ELSE
          TSAVE=100.D0*TOTALW/AREAWET
        ENDIF
        IF(AREATOT.EQ.0.0) THEN
          PSAVE=0D0
        ELSE
          PSAVE=100.D0*AREAWET/AREATOT
        ENDIF
      ENDIF
C
      const=ALPHAC(3)
      do i=1,numnp
c ..... forward difference .....
c        wthick(i)=wthick(i)+bmelt(i)*dt-const*wthick(i)*dt
c ..... backward difference .....
        wthick(i)=(wthick(i)+bmelt(i)*dt)/(1.+const*dt)
      enddo



C
C ... ACCEPT ONLY THICKNESSES GREATER THAN OR EQUAL TO ZERO
C ... ACCEPT ONLY THICKNESSES    LESS THAN OR EQUAL TO 10.0
      DO JK=1,NUMNP
         WTHICK(JK)=MAX(0D0,WTHICK(JK))
         WTHICK(JK)=MIN(10D0,WTHICK(JK))
      ENDDO
C ..  CALL EDGE DETECTOR THAT ALLOWS LEAKAGE OUT EDGES
c     (turn on/off internally)
      CALL EDGENODE(NMAX, NUMNP, NUMEL, NTYPE, KX, HTICE, DEPB, 
     &                WTHICK)
      IF(.TRUE. .AND. IPLOT.EQ.8) THEN
C       CALL NEWPAG
C       CALL CONTR(NMAX,NUMEL,X,Y,KX,
C    &                 HTICE,-2999.999D0,5500.D0,250.D0,
C    &                 -1000.D0,1000.D0,-1000.D0,1000.D0)
        CALL CONTR(NMAX,NUMEL,X,Y,KX,
     &             WTHICK,-0.0999999999D0,1.000000001D0,0.1D0,
     &            -1000.D0,1000.D0,-1000.D0,1000.D0)
      ENDIF

C ... BEGIN LOOP OVER ALL THE ELEMENTS ...
      TOTALW=0.0D0
      AREAWET=0.0D0
      AREADRY=0.0D0
      AREATOT=0.0D0
      DO 100 N=1,NUMEL
        DO I=1,4
          LM(I)=KX(N,I)
          XY(1,I)=X(LM(I))
          XY(2,I)=Y(LM(I))
        ENDDO
        CALL FESHAPE(1,0D0,0D0,PSI,DPSI)
        CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
        TWTHIK=0.D0
        DO I=1,4
          TWTHIK=TWTHIK+WTHICK(LM(I))*PSI(I)
        ENDDO
C
C ..... ACCUMULATE INTEGRATION POINT VALUES OF INTEGRALS
        DO L=1,9
          CALL FESHAPE(1,XI(1,L),ETA(1,L),PSI,DPSI)
          CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
          FAC=DETJ*W(1,L)
          TOTALW=TOTALW+TWTHIK*FAC
          IF(TWTHIK.GT.0.0) THEN
            AREAWET=AREAWET+FAC
          ELSE
            AREADRY=AREADRY+FAC
          ENDIF
          AREATOT=AREATOT+FAC
        ENDDO
C ..... ADD ELEMENT CONDUCTIVITY TO COMPLETE CONDUCTIVITY MATRIX
C
100   CONTINUE
C ... END LOOP OVER ALL THE ELEMENTS ...
C
C
      RATIOW=(TOTALW-WSAVE)/DT
      TOTALT=100.D0*TOTALW/AREAWET
      TOTALP=100.D0*AREAWET/AREATOT
      RATIOT=(TOTALT-TSAVE)/DT
      RATIOP=(TOTALP-PSAVE)/DT
      IF(IOTOGG) THEN
        WRITE(LIST(IPAGE+1),*)
     &       '********************************************'
        WRITE(LIST(IPAGE+2),*) 
     &         'TOTAL WATER (KM**3)=',REAL(1D-9*TOTALW),
     &         REAL(1D-9*RATIOW),REAL(1D-9*RATIOW*DT)
C       WRITE(LIST(IPAGE+1),*) 'WET AREA    (KM**2)=',REAL(1D-6*AREAWET)
C       WRITE(LIST(IPAGE+1),*) 'DRY AREA    (KM**2)=',REAL(1D-6*AREADRY)
C       WRITE(LIST(IPAGE+1),*) 'TOTAL AREA  (KM**2)=',REAL(1D-6*AREATOT)
        WRITE(LIST(IPAGE+3),*) 'PERCENT WET        =',REAL(TOTALP),
     &         REAL(RATIOP),REAL(RATIOP*DT)
        WRITE(LIST(IPAGE+4),*) 'AVG THICK   ( CM  )=',REAL(TOTALT),
     &         REAL(RATIOT),REAL(RATIOT*DT)
        WRITE(LIST(IPAGE+5),*) 
     &       '********************************************'
        IPAGE=IPAGE+5
      ENDIF
      WSAVE=TOTALW
      PSAVE=TOTALP
      TSAVE=TOTALT
C ... CALL CHECKER THAT PUTS WATER ON/OFF ICE-FREE NODES
c ... (LAST ARG IS VALUE TO PUT ON ICE-FREE NODES, 0.1 IT AFTER
c     (turn on/off internally)
      CALL CHECKER(NMAX, NUMNP, NUMEL, NTYPE, KX, HTICE, DEPB, 
     &             WTHICK,0.0d0)
      END




