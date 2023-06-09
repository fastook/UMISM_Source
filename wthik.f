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
