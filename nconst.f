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
