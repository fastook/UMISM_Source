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
       include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER HED*80,HEDT*80,junk*80
      DIMENSION KODE(NMAX),X(NMAX),Y(NMAX),HTICE(NMAX),ADOT(NMAX),
     &          FRACT(NMAX),PSURF(NMAX),BDROCK(NMAX),FLOWA(NMAX),
     &          SLDGB(NMAX),KX(NMAX,4),CONST(NMAX),IBFLUX(NMAX,2),
     &          BFLUX(NMAX),TBED(NMAX),AFUDGE(NMAX),ACON(NMAX),
     &          ITYPE(NMAX),GEOFLUX(NMAX),CALV(NMAX)
C .......... FORMAT STATEMENTS ..............
1000  FORMAT (A80,/,7I6,F8.0)
1001  FORMAT(I6,I4,1P2E12.5,0PF10.2,F7.2,F9.2,F10.3,F10.1,F10.5,2F10.5,
     &       I5,F10.3,F10.3,1PE10.3)
1002  FORMAT(5I6,1P2E17.10)
1007  FORMAT(2I6,E13.6)
C .......... FORMAT STATEMENTS ..............
C READ INPUT HEADER
      READ(30,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      PRINT *,HED
      IF(NUMNP.GT.NMAX) THEN
        PRINT *,'NUMNP=',NUMNP,' NMAX=',NMAX,' INCREASE NMAX'
        STOP
      ENDIF
C READ INPUT GRID, THINGS THAT NEVER CHANGE
      READ(31) HED
      PRINT *,HED
      READ(31) (KODE(I),I=1,NUMNP)
      READ(31) (X(I),I=1,NUMNP)
      READ(31) (Y(I),I=1,NUMNP)
      READ(31) (PSURF(I),I=1,NUMNP)
      READ(31) (BDROCK(I),I=1,NUMNP)
      READ(31) (KX(I,1),KX(I,2),KX(I,3),KX(I,4),I=1,NUMEL)
      READ(31) (IBFLUX(I,1),IBFLUX(I,2),BFLUX(I),I=1,NUMGBC)
C READ INPUT DIFFER, THINGS THAT DIFFER FROM ONE GRID TO THE NEXT
      READ(32) HED
      PRINT *,HED
      READ(32) (ADOT(I),I=1,NUMNP)
      READ(32) (FRACT(I),I=1,NUMNP)
      READ(32) (FLOWA(I),I=1,NUMNP)
      READ(32) (SLDGB(I),I=1,NUMNP)
      READ(32) (TBED(I),I=1,NUMNP)
      READ(32) (ITYPE(I),I=1,NUMNP)
      READ(32) (AFUDGE(I),I=1,NUMNP)
      READ(32) (GEOFLUX(I),I=1,NUMNP)
      READ(32) (CALV(I),I=1,NUMNP)
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME
85    READ(33,END=999) HED
      PRINT *,HED
      READ(33) (HTICE(I),I=1,NUMNP)
      READ(33) (ADOT(I),I=1,NUMNP)
      READ(33) (BDROCK(I),I=1,NUMNP)
      READ(33) (CONST(I),I=1,NUMEL)
      READ(33) (ACON(I),I=1,NUMEL)
      READ(34) (FRACT(I),I=1,NUMNP)
      READ(34) (FLOWA(I),I=1,NUMNP)
      READ(34) (SLDGB(I),I=1,NUMNP)
      READ(34) (AFUDGE(I),I=1,NUMNP)
      READ(36,END=999) HEDT
      READ(36) (TBED(I),I=1,NUMNP)
      READ(36) (BMELT,I=1,NUMNP)
      READ(36) (WTHICK,I=1,NUMNP)
      print *,HEDT
      PRINT *,HED
      PRINT *,'INPUT 1 TO PRINT, 0 TO SKIP, -1 TO STOP'
      READ(*,*) IPRINT
      IF(IPRINT.LT.0) STOP
c     iprint=1
      IF(IPRINT.EQ.1) THEN
        print *,'writing out ',hedt
        WRITE(1,1000) HEDT,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
        DO N=1,NUMNP
          WRITE(1,1001) N,KODE(N),X(N),Y(N),HTICE(N),
     &    ADOT(N),FRACT(N),PSURF(N),BDROCK(N),FLOWA(N),SLDGB(N),
     &    TBED(N),ITYPE(n),AFUDGE(N),GEOFLUX(N),CALV(N)
        ENDDO
        DO N=1,NUMEL
          WRITE(1,1002) N,KX(N,1),KX(N,2),KX(N,3),KX(N,4),CONST(N),
     &                 ACON(N)
        ENDDO
        IF(NUMGBC.GT.0) THEN
          DO N=1,NUMGBC
            WRITE(1,1007) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
          ENDDO
        ENDIF
      endif
      GOTO 85
999   CONTINUE
      END
