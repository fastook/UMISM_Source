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
