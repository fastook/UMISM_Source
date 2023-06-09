      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      IMPLICIT REAL*4(A-H,O-Z)
      CHARACTER HED*80
      DIMENSION KODE(NMAX),X(NMAX),Y(NMAX),HTICE(NMAX),ADOT(NMAX),
     &          FRACT(NMAX),PSURF(NMAX),BDROCK(NMAX),FLOWA(NMAX),
     &          SLDGB(NMAX),KX(NMAX,4),CONST(NMAX),AFUDGE(NMAX),
     &          IBFLUX(NMAX,2),BFLUX(NMAX),temp(nmax),itype(nmax),
     &          ACON(NMAX),icmap(16)
      DATA ICMAP /0,5,11,4,12,6,13,2,8,7,9,3,10,14,15,1/
      READ(1,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      PRINT *,HED
1000  FORMAT (A80,/,7I6,F8.1)
      IF(NUMNP.GT.NMAX) THEN
        PRINT *,'NUMNP=',NUMNP,' NMAX=',NMAX,' INCREASE NMAX'
        STOP
      ENDIF
      print *,numnp,' nodes,',numel,' elements'
      amin=1e30
      amax=-amin
      xmin=1e30
      xmax=-xmin
      ymin=1e30
      ymax=-xmin
      DO N=1,NUMNP
        READ(1,1001) NUM,KODE(N),X(N),Y(N),HTICE(N),
     &                   ADOT(N),FRACT(N),PSURF(N),
     &                   BDROCK(N),FLOWA(N),SLDGB(N),
     &                   temp(n),itype(n),AFUDGE(N)
        xmin=min(xmin,x(n))
        xmax=max(xmax,x(n))
        ymin=min(ymin,y(n))
        ymax=max(ymax,y(n))
        amin=min(amin,flowa(n))
        amax=max(amax,flowa(n))
      enddo
1001  FORMAT(I6,I4,1P2E12.5,0PF10.2,F7.2,F9.2,F10.3,F10.1,F10.5,2F10.5,
     &       I5,F10.3)
      DO N=1,NUMEL
        READ(1,1002) NUM,(KX(NUM,i),i=1,4),CONST(N),ACON(N)
        amin=min(amin,acon(n))
        amax=max(amax,acon(n))
      enddo
 1002 FORMAT(5I6,1P2E17.10)
      IF(NUMGBC.GT.0) THEN
        DO N=1,NUMGBC
          READ(1,1007) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
 1007 FORMAT(2I6,E13.6)
        enddo
      ENDIF
      xd=(xmax-xmin)/20.
      print *,xmin,xmax      
      print *,ymin,ymax      
      print *,amin,amax      
      call grstrt(800,800)
      call window(xmin-xd,xmax+xd,ymin-xd,ymax+xd)
      zf=amin
      zd=1.
      do k=1,30
      DO I=1,NUMNP
        RZ=(flowa(I)-ZF)/ZD
        IZ=RZ+2
        IF(IZ.LT.2) IZ=2
        IF(IZ.GT.16) IZ=16
        CALL MRKCLR(ICMAP(IZ))
        XXX=x(I)
        YYY=y(I)
        CALL BMARKER(REAL(XXX),REAL(YYY),1)
      ENDDO
c      call newpag
      DO I=1,numel
        RZ=(acon(I)-ZF)/ZD
        IZ=RZ+2
        IF(IZ.LT.2) IZ=2
        IF(IZ.GT.16) IZ=16
        CALL MRKCLR(ICMAP(IZ))
        XXX=0.
        YYY=0.
        do j=1,4
          xxx=xxx+x(kx(i,j))
          yyy=yyy+y(kx(i,j))
        enddo
        xxx=0.25*xxx
        yyy=0.25*yyy
        CALL BMARKER(REAL(XXX),REAL(YYY),1)
      ENDDO
c      call newpag
      enddo


      call grstop1
      END
