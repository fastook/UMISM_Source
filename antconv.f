      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER HED*80,ice*3
      DIMENSION KODE(NMAX),X(NMAX),Y(NMAX),HTICE(NMAX),ADOT(NMAX),
     &          FRACT(NMAX),PSURF(NMAX),BDROCK(NMAX),FLOWA(NMAX),
     &          SLDGB(NMAX),KX(NMAX,4),CONST(NMAX),IBFLUX(NMAX,2),
     &          BFLUX(NMAX),temp(nmax),itype(nmax),AFUDGE(NMAX)
      READ(1,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      WRITE(11,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      PRINT *,HED
      call setrig
1000  FORMAT (A80,/,7I6,F8.1)
      IF(NUMNP.GT.NMAX) THEN
        PRINT *,'NUMNP=',NUMNP,' NMAX=',NMAX,' INCREASE NMAX'
        STOP
      ENDIF
      DO 100 N=1,NUMNP
        READ(1,1001) NUM,KODE(N),X(N),Y(N),HTICE(N),
     &                   ADOT(N),FRACT(N),PSURF(N),
     &                   BDROCK(N),FLOWA(N),SLDGB(N),
     &                   temp(n),itype(n),AFUDGE(N)
1001  FORMAT(I6,I4,2E12.5,F10.2,F7.2,F9.2,F10.3,F10.1,F10.5,2F10.5,
     &       i5,F10.3)
c put modifications here
c from here
c
       if(kode(n).eq.1) then
         adot(n)=-5.
         kode(n)=0
       endif
       if(psurf(n).eq.0.) then
          afudge(n)=0.1
          adot(n)=-5.
          kode(n)=1
          fract(n)=1.
          itype(n)=0
       else
c          if(adot(n).ne.-5.) adot(n)=-120.
          itype(n)=2
           if(kode(n).eq.0) then
            kode(n)=0
            fract(n)=0.01
          else
            fract(n)=1.
          endif
        endif
       if(kode(n).eq.1 .and. adot(n).eq.-120.) kode(n)=0
       if(psurf(n).eq.0. .and. bdrock(n).gt.-1000) kode(n)=0
       if(kode(n).eq.0) itype(n)=2
123   continue
      RHOW=1.092
      RHOI=0.917
      IF(bdrock(n).LT.0.) THEN
        FLOT=(1.-RHOW/RHOI)*bdrock(n)
        if(psurf(n).lt.flot) psurf(n)=0.
      ENDIF



c
c to here is standard antarctic correction

c to here
       WRITE(11,1001) NUM,KODE(N),X(N),Y(N),HTICE(N),
     &                   ADOT(N),FRACT(N),PSURF(N),
     &                   BDROCK(N),FLOWA(N),SLDGB(N),
     &                   temp(n),itype(n),AFUDGE(N)
100   CONTINUE
      DO 90 N=1,NUMEL
      READ(1,1002) NUM,KX(NUM,1),KX(NUM,2),KX(NUM,3),KX(NUM,4),CONST(N)
      WRITE(11,1003) NUM,KX(NUM,1),KX(NUM,2),KX(NUM,3),KX(NUM,4),
     &               CONST(N)
 1002 FORMAT(5I6,E17.10)
 1003 FORMAT(5I6,1PE17.10)
 90   CONTINUE
      IF(NUMGBC.GT.0) THEN
        DO 210 N=1,NUMGBC
          READ(1,1007) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
          WRITE(11,1007) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
 1007 FORMAT(2I6,E13.6)
210     CONTINUE
      ENDIF
      END

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
