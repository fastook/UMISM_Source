      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER HED*80,ice*3
      DIMENSION KODE(NMAX),X(NMAX),Y(NMAX),HTICE(NMAX),ADOT(NMAX),
     &          FRACT(NMAX),PSURF(NMAX),BDROCK(NMAX),FLOWA(NMAX),
     &          SLDGB(NMAX),KX(NMAX,4),CONST(NMAX),IBFLUX(NMAX,2),
     &          BFLUX(NMAX),temp(nmax),itype(nmax),AFUDGE(NMAX),
     &          vel(nmax),geoflux(nmax),calv(nmax)
C
C--------------- READ IN OLD FROM UNIT 1 -------------------------
C
      RHOW=1.092
      RHOI=0.917
      READ(1,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      PRINT *,HED
      call setrig
1000  FORMAT (A80,/,7I6,F8.1)
      IF(NUMNP.GT.NMAX) THEN
        PRINT *,'NUMNP=',NUMNP,' NMAX=',NMAX,' INCREASE NMAX'
        STOP
      ENDIF
      DO N=1,NUMNP
        READ(1,1001) NUM,KODE(N),X(N),Y(N),HTICE(N),
     &                   ADOT(N),FRACT(N),PSURF(N),
     &                   BDROCK(N),FLOWA(N),SLDGB(N),
     &                   temp(n),itype(n),AFUDGE(N),GEOFLUX(N),
     &                   calv(n)
1001  FORMAT(I6,I4,1P2E12.5,0PF10.2,F7.2,F9.2,F10.3,F10.1,F10.5,2F10.5,
     &       I5,F10.3,F10.3,1PE10.3)
      enddo
      DO N=1,NUMEL
        READ(1,1002) NUM,KX(NUM,1),KX(NUM,2),KX(NUM,3),KX(NUM,4),
     &               CONST(N),vel(n)
      enddo
 1002 FORMAT(5I6,2E17.10)
      IF(NUMGBC.GT.0) THEN
        DO N=1,NUMGBC
          READ(1,1007) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
 1007 FORMAT(2I6,E13.6)
        enddo
      ENDIF
C --------------- END OF OLD READ ---------------------
C
C ---------------- MODIFICATIONS -------------------------
C
      do N=1,numnp
c put modifications here
c      read(2,2000) nn,ibed,ice
c2000  format(i3,i6,1x,a3)
c      if(ibed.eq.-5000) ibed=4999.
c      if(ice.eq.'ice') then
c        print *,'an ice point'
c        htice(n)=real(ibed)
c        psurf(n)=real(ibed)
c        bdrock(n)=real(ibed)-100.
c        print *,n,psurf(n),bdrock(n)
c      else
c        htice(n)=real(ibed)
c        psurf(n)=real(ibed)
c        bdrock(n)=real(ibed)
c      endif
c     psurf(n)=htice(n)
c     htice(n)=psurf(n)
c     if(adot(n).lt.0.) then
c       kode(n)=1
c     else
c       kode(n)=0
c     endif
c      fract(n)=0.
c      call RECPOL(.001*X(n),.001*Y(n),RLAT,RLONG)
c      if(rlat.gt.50.) fract(n)=1.
c 
c      call RECPOL(.001*X(n),.001*Y(n),RLAT,RLONG)
c      if(rlat.gt.72.) then
c        adot(n)=.1
c      else
c        adot(n)=-1.0
c      endif
c
c      if(bdrock(n).lt.-0. .or. bdrock(n).gt.500.) fract(n)=1.
c       if(bdrock(n).lt.0.) then
c         fract(n)=1.
c       else
c         fract(n)=0.
c       endif
c        if(bdrock(n).lt.0.) adot(n)=-1.
c      if(bdrock(n).lt.-500 .and. psurf(n).eq.0.0) then
c        adot(n)=-10.
c      else
c        adot(n)=-120.
c      endif
c      if(bdrock(n).gt.0 .and. bdrock(n).lt.500.) then
c        fract(n)=0.
c      elseif(bdrock(n).gt.500 .and. psurf(n).gt.bdrock(n)) then
c        fract(n)=0.
c      elseif(bdrock(n).gt.500 .and. psurf(n).eq.bdrock(n)) then
c        fract(n)=1.
c      elseif(bdrock(n).gt.0 .and. psurf(n).gt.bdrock(n)) then
c        fract(n)=0.
c      elseif(bdrock(n).gt.0 .and. psurf(n).eq.bdrock(n)) then
c        fract(n)=0.
c      elseif(bdrock(n).le.0. .and. psurf(n).eq.0.) then
c        fract(n)=1.
c      elseif(bdrock(n).le.0. .and. psurf(n).ne.0.) then
c        fract(n)=0.
c      else
c        print *,n,'didnt find',bdrock(n),psurf(n)
c      endif
c      if(bdrock(n).lt.-0. .or. bdrock(n).gt.500.) then
c        if(bdrock(n).lt.0 .and. psurf(n).eq.0) then
c          fract(n)=0.
c        elseif(bdrock(n).gt.0 .and. psurf(n).gt.bdrock(n)) then
c          fract(n)=0.
c        else
c          fract(n)=1.
c        endif
c      else
c        fract(n)=0.
c      endif
c      psurf(n)=htice(n)
c      if(bdrock(n).gt.100) then
c        fract(n)=0.
c      else
c        fract(n)=1.
c      endif
c
       if(.false.) then
         IF(bdrock(n).LT.0.) THEN
           if(psurf(n).lt.100.) htice(n)=0.
           if(psurf(n).lt.100.) psurf(n)=0.
           FLOT=(1.-RHOW/RHOI)*bdrock(n)
           if(psurf(n).lt.flot) then
             psurf(n)=0.
             htice(n)=0.
             afudge(n)=.1
             kode(n)=1
           else
             htice(n)=psurf(n)
             afudge(n)=1.
             kode(n)=0
           endif
         ENDIF
       endif 
c       if(psurf(n).eq.0.) then
c          kode(n)=1
c       endif
c from here standard antarctic conversion ...
c
       if(.false.) then
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
            if(adot(n).ne.-5.) adot(n)=-120.
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
        endif
c
c to here is standard antarctic correction
c switches all above threshold elevation to type=2, below, type=3
      if(.true.) then
        if(htice(n).gt.2000.) then
          itype(n)=2
        else
          itype(n)=3
        endif
      endif
      if(.true.) then
           FLOT=(1.-RHOW/RHOI)*bdrock(n)
           if(psurf(n).lt.flot) then
             calv(n)=0.1
           else
             calv(n)=0.001
           endif
      endif
      enddo
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
C
C ---------------- MODIFICATIONS TO HERE -------------------------
C
c
c -------------- NOW WRITE NEW SET OUT TO 11 -----------------------
c
      WRITE(11,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      do N=1,numnp
       WRITE(11,1001) n,KODE(N),X(N),Y(N),HTICE(N),
     &                   ADOT(N),FRACT(N),PSURF(N),
     &                   BDROCK(N),FLOWA(N),SLDGB(N),
     &                   temp(n),itype(n),AFUDGE(N),
     &                   geoflux(n),calv(n)
      enddo
      do n=1,numel
        WRITE(11,1003) n,KX(n,1),KX(n,2),KX(n,3),KX(n,4),
     &                 CONST(N),vel(n)
      enddo
 1003 FORMAT(5I6,1P2E17.10)
      IF(NUMGBC.GT.0) THEN
        DO N=1,NUMGBC
          WRITE(11,1007) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
        enddo
      endif
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
