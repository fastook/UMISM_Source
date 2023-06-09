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

