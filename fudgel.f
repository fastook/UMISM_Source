      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER HED*80
      DIMENSION KODE(NMAX),X(NMAX),Y(NMAX),HTICE(NMAX),ADOT(NMAX),
     &          FRACT(NMAX),PSURF(NMAX),BDROCK(NMAX),FLOWA(NMAX),
     &          SLDGB(NMAX),KX(NMAX,4),CONST(NMAX),IBFLUX(NMAX,2),
     &          BFLUX(NMAX),temp(nmax),itype(nmax),AFUDGE(NMAX)
1000  FORMAT (A80,/,7I6,F8.1)
1001  FORMAT(I6,I4,2E12.5,F10.2,F7.2,F9.2,F10.3,F10.1,F10.5,2F10.5,
     &       i5,F10.3)
1002  FORMAT(5I6,2E17.10)
1003  FORMAT(5I6,1P2E17.10)
1007  FORMAT(2I6,E13.6)
      READ(1,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      WRITE(11,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      PRINT *,HED
      print *,' input upper limit on fudge'
      read(*,*) fudgel
      IF(NUMNP.GT.NMAX) THEN
        PRINT *,'NUMNP=',NUMNP,' NMAX=',NMAX,' INCREASE NMAX'
        STOP
      ENDIF
      DO N=1,NUMNP
        READ(1,1001) NUM,KODE(N),X(N),Y(N),HTICE(N),
     &                   ADOT(N),FRACT(N),PSURF(N),
     &                   BDROCK(N),FLOWA(N),SLDGB(N),
     &                   temp(n),itype(n),AFUDGE(N)
c put modifications here
      if(afudge(n).gt.fudgel) afudge(n)=fudgel

c to here
       WRITE(11,1001) NUM,KODE(N),X(N),Y(N),HTICE(N),
     &                   ADOT(N),FRACT(N),PSURF(N),
     &                   BDROCK(N),FLOWA(N),SLDGB(N),
     &                   temp(n),itype(n),AFUDGE(N)
      enddo
      DO N=1,NUMEL
        READ(1,1002) NUM,KX(NUM,1),KX(NUM,2),KX(NUM,3),KX(NUM,4),
     &             CONST(N),vel
        WRITE(11,1003) NUM,KX(NUM,1),KX(NUM,2),KX(NUM,3),KX(NUM,4),
     &               CONST(N),vel
      enddo
      IF(NUMGBC.GT.0) THEN
        DO N=1,NUMGBC
          READ(1,1007) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
          WRITE(11,1007) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
        enddo
      ENDIF
      END
