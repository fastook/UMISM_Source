C CONVERTS MULTIPLE OUTPUT OF MAP5 BACK INTO SINGLE DATA SET FORMAT
C INPUT . . .
C UNIT 30 <= INPUT&1 HEADER B
C     HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
C UNIT 31 <= INPUT&1 GRID B (THINGS THAT NEVER CHANGE
C     HED
C     (KODE(I),I=1,NUMNP)
C     (X(I),I=1,NUMNP)
C     (Y(I),I=1,NUMNP)
C     (PSURF(I),I=1,NUMNP)
C     (BDROCK(I),I=1,NUMNP)
C     (KX(I,1),KX(I,2),KX(I,3),KX(I,4),I=1,NUMEL)
C     (IBFLUX(I,1),IBFLUX(I,2),BFLUX(I),I=1,NUMGBC)
C UNIT 32 <= INPUT&1 DIFF B (THINGS THAT DIFFER FROM GRID TO GRID
C     HED
C     (ADOT(I),I=1,NUMNP)
C     (FRACT(I),I=1,NUMNP)
C     (FLOWA(I),I=1,NUMNP)
C     (SLDGB(I),I=1,NUMNP)
C UNIT 33 <= OUT&2 TIME B (THINGS THAT CHANGE WITH TIME AT TIME STEPS
C     HED
C     (HTICE(I),I=1,NUMNP)
C     (ADOT(I),I=1,NUMNP)
C     (BDROCK(I),I=1,NUMNP)
C     (CONST(I),I=1,NUMEL)
C     (ACON(I),I=1,NUMEL)
C OUTPUT . . .
C UNIT 1 >= OUT3&3 DATA B  (SINGLE DATA SET IN STANDARD FORMAT
C UNITS 31,32,33 ARE BINARY, UNIT 1 (RECFM F LRECL 130
C RUN BY EXEC PACK EXEC A:
C           FI 1 DISK OUT3&3 DATA B (LRECL 130 RECFM F
C           FI 30 DISK INPUT&1 HEAD B
C           FI 31 DISK INPUT&1 GRID B
C           FI 32 DISK INPUT&1 DIFF B
C           FI 33 DISK OUT&2 TIME B
C           GL TXT VSF2FORT CMSLIB
C           LOAD PACK (START NOMAP
      PARAMETER(NMAX=29999)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER HED*80,HED1*80,JUNK*80
      DIMENSION KODE(NMAX),X(NMAX),Y(NMAX),HTICE(NMAX),ADOT(NMAX),
     &          FRACT(NMAX),PSURF(NMAX),BDROCK(NMAX),FLOWA(NMAX),
     &          SLDGB(NMAX),KX(NMAX,4),CONST(NMAX),BED(NMAX),
     &          IBFLUX(NMAX,2),BFLUX(NMAX),PFRACT(NMAX),
     &          PICE(NMAX),TBED(NMAX),PTEMP(NMAX),ATIME(NMAX),
     &          ACON(NMAX)
      PRINT *,'INPUT 0 FOR ICE'
      PRINT *,'      1 FOR FROZEN'
      PRINT *,'      2 FOR TEMPERATURE'
      PRINT *,'        AND THRESHOLD'
      READ(*,*) ITYPE,THRESH
      RHOW=1.092
      RHOI=0.917
C READ INPUT HEADER
      READ(30,1000) HED1,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      PRINT *,HED1
1000  FORMAT (A80,/,7I6,F8.0)
      IF(NUMNP.GT.NMAX) THEN
        PRINT *,'NUMNP=',NUMNP,' NMAX=',NMAX,' INCREASE NMAX'
        STOP
      ENDIF
      DO I=1,NUMNP
        ATIME(I)=0.
      ENDDO
C READ INPUT GRID, THINGS THAT NEVER CHANGE
      READ(31) HED
      PRINT *,HED
      READ(31) (KODE(I),I=1,NUMNP)
      READ(31) (X(I),I=1,NUMNP)
      READ(31) (Y(I),I=1,NUMNP)
      READ(31) (PSURF(I),I=1,NUMNP)
      READ(31) (BED(I),I=1,NUMNP)
      READ(31) (KX(I,1),KX(I,2),KX(I,3),KX(I,4),I=1,NUMEL)
      READ(31) (IBFLUX(I,1),IBFLUX(I,2),BFLUX(I),I=1,NUMGBC)
C READ INPUT DIFFER, THINGS THAT DIFFER FROM ONE GRID TO THE NEXT
      READ(32) HED
      PRINT *,HED
      READ(32) (ADOT(I),I=1,NUMNP)
      READ(32) (FRACT(I),I=1,NUMNP)
      READ(32) (FLOWA(I),I=1,NUMNP)
      READ(32) (SLDGB(I),I=1,NUMNP)
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME
      ICOUNT=0
85    READ(33,END=999) HED
        ICOUNT=ICOUNT+1
        WRITE(JUNK,*) HED
        READ(JUNK,1010) TIME
        PRINT *,HED
        READ(33) (HTICE(I),I=1,NUMNP)
        READ(33) (ADOT(I),I=1,NUMNP)
        READ(33) (BDROCK(I),I=1,NUMNP)
        READ(33) (CONST(I),I=1,NUMEL)
        READ(33) (ACON(I),I=1,NUMEL)
        READ(34) (FRACT(I),I=1,NUMNP)
        READ(34) (FLOWA(I),I=1,NUMNP)
        READ(34) (SLDGB(I),I=1,NUMNP)
        READ(34) (AFUDGE,I=1,NUMNP)
        READ(36,END=999) HEDT
        READ(36) (TBED(I),I=1,NUMNP)
        READ(36) (BMELT,I=1,NUMNP)
        READ(36) (WTHICK,I=1,NUMNP)
C--------------------------------------------
        IF(ICOUNT.EQ.1) THEN
          TBASE=TIME
          TOLD=TIME
        ELSE
          DT=TIME-TOLD
          DO I=1,NUMNP
            IF(ITYPE.EQ.1) THEN
              IF(BED(I).LT. 0.) THEN
                FLOT=BDROCK(I)*(1-RHOW/RHOI)
                IF(HTICE(I)-FLOT.GT.0.0) THEN
                  ATIME(I)=ATIME(I)+DT
                ENDIF
              ELSE
                IF(HTICE(I)-BDROCK(I).GT.0.0) THEN
                  ATIME(I)=ATIME(I)+DT
                ENDIF
              ENDIF
              IF(FRACT(I).LE. THRESH) THEN
                PFRACT(I)=PFRACT(I)+DT
              ENDIF
            ELSEIF(ITYPE.EQ.0) THEN
              IF(BED(I).LT. 0.) THEN
                FLOT=BDROCK(I)*(1-RHOW/RHOI)
                IF(HTICE(I)-FLOT.GT.THRESH) THEN
                  PICE(I)=PICE(I)+DT
                ENDIF
              ELSE
                IF(HTICE(I)-BDROCK(I).GT.THRESH) THEN
                  PICE(I)=PICE(I)+DT
                ENDIF
              ENDIF
            ELSEIF(ITYPE.EQ.2) THEN
c              IF(BED(I).LT. 0.) THEN
c                FLOT=BDROCK(I)*(1-RHOW/RHOI)
c                IF(HTICE(I)-FLOT.GT.0.0) THEN
c                  ATIME(I)=ATIME(I)+DT
c                  IF(TBED(I).LT.THRESH) THEN
c                    PTEMP(I)=PTEMP(I)+DT
c                  ENDIF
c                ENDIF
c              ELSE
c                IF(HTICE(I)-BDROCK(I).GT.0.0) THEN
c                  ATIME(I)=ATIME(I)+DT
c                  IF(TBED(I).LT.THRESH) THEN
c                    PTEMP(I)=PTEMP(I)+DT
c                  ENDIF
c                ENDIF
c              ENDIF
              IF(TBED(I).LT.THRESH) THEN
                PTEMP(I)=PTEMP(I)+DT
              ENDIF
            ENDIF
          ENDDO
          TOLD=TIME
        ENDIF
      GOTO 85
999   CONTINUE
      DO I=1,NUMNP
        PFRACT(I)=PFRACT(I)/(TIME-TBASE)
        PTEMP(I)=PTEMP(I)/(TIME-TBASE)
        PICE(I)=PICE(I)/(TIME-TBASE)
c        IF(ATIME(I).NE.0.0) THEN
c          PTEMP(I)=PTEMP(I)/ATIME(I)
c          PFRACT(I)=PFRACT(I)/ATIME(I)
c        ELSE
c          PTEMP(I)=0.0
c          PFRACT(I)=0.0
c        ENDIF
      ENDDO
      WRITE(1,1000) HED1,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      DO N=1,NUMNP
        IF(ITYPE.EQ.0) THEN
          WRITE(1,1001) N,KODE(N),X(N),Y(N),HTICE(N),
     &                  ADOT(N),PICE(N),PSURF(N),
     &                  BDROCK(N),FLOWA(N),SLDGB(N)
        ELSEIF(ITYPE.EQ.1) THEN
          WRITE(1,1001) N,KODE(N),X(N),Y(N),HTICE(N),
     &                  ADOT(N),PFRACT(N),PSURF(N),
     &                  BDROCK(N),FLOWA(N),SLDGB(N)
        ELSEIF(ITYPE.EQ.2) THEN
          WRITE(1,1001) N,KODE(N),X(N),Y(N),HTICE(N),
     &                  ADOT(N),PTEMP(N),PSURF(N),
     &                  BDROCK(N),FLOWA(N),SLDGB(N)
        ENDIF
      ENDDO
      DO N=1,NUMEL
        WRITE(1,1002) N,KX(N,1),KX(N,2),KX(N,3),KX(N,4),CONST(N),
     &                ACON(N)
      ENDDO
      IF(NUMGBC.GT.0) THEN
        DO N=1,NUMGBC
          WRITE(1,1007) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
        ENDDO
      ENDIF
1001  FORMAT(I6,I4,1P2E12.5,0PF10.2,F7.2,F9.2,F10.3,F10.1,F10.5,2F10.5,
     &       I5,F10.3)
1002  FORMAT(5I6,1P2E17.10)
1007  FORMAT(2I6,E13.6)
1010  FORMAT(6X,F20.0)
      END
