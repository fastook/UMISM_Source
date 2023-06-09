C PROGRAM TO CONVERT OLD DATA SETS (SINGLE COMPLETE) TO NEW MULTIPLE
C DATA SETS.
C INPUT CONSISTS OF:
C   UNIT 1 <= INPUT&1 DATA B (RECFM F LRECL 130
C TITLE  CARD 1 (A80) HED,7I6,F8.1)
C        CARD 2 (7I6,F8.1)
C              NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
C NUMNP  CARD 3'S (I6,I4,2E12.5,F10.2,F7.2,F9.2,F10.3,F10.1,F10.2,F10.3)
C              N,KODE(N),X(N),Y(N),HTICE(N),
C              ADOT(N),FRACT(N),PSURF(N),BDROCK(N),FLOWA(N),SLDGB(N)
C NUMEL  CARD 4'S (5I6,1PE17.10)
C              N,KX(N,1),KX(N,2),KX(N,3),KX(N,4),CONST(N)
C NUMGBC CARD 5'S (2I6,E13.6)
C              IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
C
C OUTPUT CONSISTS OF:
C   UNIT 30 => INPUT&1 HEAD B (INPUT HEADER
C              HED           = TITLE 80 CHARACTERS LONG
C              NUMNP         = NUMBER OF NODAL POINTS
C              NUMEL         = NUMBER OF ELEMENTS
C              NUMCOL,NUMLEV = NUMBER OF COLUMNS AND ROWS IF APPLICABLE
C              NUMGBC        = NUMBER OF FLUX BOUNDARY CONDITIONS
C              NDT           = NUMBER OF TIME STEPS
C              INTER         = HOW OFTEN TO OUTPUT SOLUTION
C              DT            = TIME STEP
C   UNIT 31 => INPUT&1 GRID B (THINGS THAT NEVER CHANGE
C              HED                    = TITLE
C              (KODE(I),I=1,NUMNP)    = NODE FIXED OR FREE
C              (X(I),I=1,NUMNP)       = X COORDINATES
C              (Y(I),I=1,NUMNP)       = Y COORDINATES
C              (PSURF(I),I=1,NUMNP)   = PRESENT SURFACE
C              (BDROCK(I),I=1,NUMNP)  = PRESENT BEDROCK
C              (KX(I,1),KX(I,2),KX(I,3),KX(I,4),I=1,NUMEL)
C                                     = CONNECTIVITY ARRAY
C              (IBFLUX(I,1),IBFLUX(I,2),BFLUX(I),I=1,NUMGBC)
C                                     = EDGE NODES AND DEFINED FLUX
C   UNIT 32 => INPUT&1 DIFF B (THINGS THAT DIFFER FROM ONE GRID TO
C                                 ANOTHER
C              HED                  = TITLE
C              (ADOT(I),I=1,NUMNP)  = MASS BALANCE, EITHER RATE OR CODE
C              (FRACT(I),I=1,NUMNP) = FRACTION OF BED MELTED
C              (FLOWA(I),I=1,NUMNP) = FLOW LAW CONSTANT
C              (SLDGB(I),I=1,NUMNP) = SLIDING LAW CONSTANT
C   UNIT 33 => INPUT&1 TIME B (THINGS THAT CHANGE WITH TIME
C              HED                  = TITLE
C              (HTICE(I),I=1,NUMNP) = INITIAL ICE SURFACE
C              (ADOT(I),I=1,NUMNP)  = MASS BALANCE, EITHER RATE OR CODE
C              (BDROCK(I),I=1,NUMNP)= BEDROCK (MAY BE DEPRESSED)
C              (CONST(I),I=1,NUMEL) = LINEARIZATION CONSTANT
C              (AFUDGE(I),I=1,NUMNP)= ICE HARDNESS FUDGE FACTOR
C UNITS 31,32,33 ARE BINARY AND ARE MACHINE READABLE ONLY
C RUN BY UNPACK EXEC A AS FOLLOWS:
C          FI 1 DISK INPUT&1 DATA B (LRECL 130 RECFM F
C          FI 30 DISK INPUT&2 HEAD B
C          FI 31 DISK INPUT&2 GRID B
C          FI 32 DISK INPUT&2 DIFF B
C          FI 33 DISK INPUT&2 TIME B
C          GL TXT VSF2FORT CMSLIB
C          LOAD UNPACK (START NOMAP
      include "parameter.h" 
      PARAMETER(NMAX=MAXNUM)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER HED*80
      DIMENSION KODE(NMAX),X(NMAX),Y(NMAX),HTICE(NMAX),ADOT(NMAX),
     &          FRACT(NMAX),PSURF(NMAX),BDROCK(NMAX),FLOWA(NMAX),
     &          SLDGB(NMAX),KX(NMAX,4),CONST(NMAX),AFUDGE(NMAX),
     &          IBFLUX(NMAX,2),BFLUX(NMAX),temp(nmax),itype(nmax),
     &          ACON(NMAX),GEOFLUX(NMAX),calv(nmax),
     &          ACC(NMAX),ABLAT(NMAX),CLAT(nmax),CLONG(NMAX)
      dimension lm(5),xy(2,4)
      DIMENSION PSI(4),DPSI(4,2)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2)
      print *,'input MAT in Kelvin'
      read *, rmat
      rmat = rmat -273
      print *,'input geothermal flux in mW/m^2'
      read *,flux
      flux = flux *3.77e5/50
c     READ(1,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      READ(1,'(a)') HED
      READ(1,*) NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
1000  FORMAT (A80,/,7I6,F8.1)
      IF(NUMNP.GT.NMAX) THEN
        PRINT *,'NUMNP=',NUMNP,' NMAX=',NMAX,' INCREASE NMAX'
        STOP
      ENDIF
      print *,numnp,' nodes,',numel,' elements'
      DO N=1,NUMNP
c       READ(1,1001) NUM,KODE(N),X(N),Y(N),HTICE(N),
        READ(1,*) NUM,KODE(N),X(N),Y(N),HTICE(N),
     &                   ACC(N),FRACT(N),PSURF(N),
     &                   BDROCK(N),FLOWA(N),SLDGB(N),
     &                   TEMP(N),ITYPE(N),AFUDGE(N),GEOFLUX(N),
     &                   ABLAT(N),CLAT(N),CLONG(N)
        TEMP(N) = rmat
        GEOFLUX(N) = flux
        ADOT(N)=ACC(N)-ABLAT(N)
        CALV(N)=0.
1001  FORMAT(I6,I4,1P2E12.5,0PF10.2,F7.2,F9.2,F10.3,F10.1,F10.5,2F10.5,
     &       I5,F10.3,F10.3,1PE10.3)
      enddo
      DO N=1,NUMEL
c       READ(1,1002) NUM,(KX(NUM,i),i=1,4),CONST(N),ACON(N)
        READ(1,*) NUM,(KX(NUM,i),i=1,4),CONST(N),ACON(N)
 1002 FORMAT(5I6,1P2E17.10)
      enddo
      IF(NUMGBC.GT.0) THEN
        DO N=1,NUMGBC
c         READ(1,1007) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
          READ(1,*) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
 1007 FORMAT(2I6,E13.6)
        enddo
      ENDIF
C WRITE INPUT HEADER
      WRITE(30,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
C WRITE INPUT GRID, THINGS THAT NEVER CHANGE
      WRITE(31) HED
      WRITE(31) (KODE(I),I=1,NUMNP)
      WRITE(31) (X(I),I=1,NUMNP)
      WRITE(31) (Y(I),I=1,NUMNP)
      WRITE(31) (PSURF(I),I=1,NUMNP)
      WRITE(31) (BDROCK(I),I=1,NUMNP)
      WRITE(31) (KX(I,1),KX(I,2),KX(I,3),KX(I,4),I=1,NUMEL)
      WRITE(31) (IBFLUX(I,1),IBFLUX(I,2),BFLUX(I),I=1,NUMGBC)
C WRITE INPUT DIFFER, THINGS THAT DIFFER FROM ONE GRID TO THE NEXT
      WRITE(32) HED
      WRITE(32) (ADOT(I),I=1,NUMNP)
      WRITE(32) (FRACT(I),I=1,NUMNP)
      WRITE(32) (FLOWA(I),I=1,NUMNP)
      WRITE(32) (SLDGB(I),I=1,NUMNP)
      WRITE(32) (temp(I),I=1,NUMNP)
      WRITE(32) (itype(I),I=1,NUMNP)
      WRITE(32) (AFUDGE(I),I=1,NUMNP)
      WRITE(32) (GEOFLUX(I),I=1,NUMNP)
      WRITE(32) (calv(I),I=1,NUMNP)
C WRITE INPUT TIME, THINGS THAT CHANGE WITH TIME
      WRITE(33) HED
      WRITE(33) (HTICE(I),I=1,NUMNP)
      WRITE(33) (ACC(I),I=1,NUMNP)
      WRITE(33) (ABLAT(I),I=1,NUMNP)
      WRITE(33) (CONST(I),I=1,NUMEL)
      WRITE(33) (ACON(I),I=1,NUMEL)
C WRITE INPUT COORD, lat and long of grid point
      if(CLAT(1).ne.0.) then
        WRITE(34) HED
        WRITE(34) (CLAT(I),I=1,NUMNP)
        WRITE(34) (CLONG(I),I=1,NUMNP)
      endif
      xmin=1e30
      xmax=-1e30
      ymin=1e30
      ymax=-1e30
      do i=1,numnp
        xmin=min(xmin,x(i))
        xmax=max(xmax,x(i))
        ymin=min(ymin,y(i))
        ymax=max(ymax,y(i))
      enddo
      call setrig
      open(98,file='last-rn2.data')
      call recpol(xmin,ymin,rlat,rlong)
      print *,rlat,rlong
      write(98,*) rlat,rlong
      call recpol(xmax,ymin,rlat,rlong)
      print *,rlat,rlong
      write(98,*) rlat,rlong
      call recpol(xmax,ymax,rlat,rlong)
      print *,rlat,rlong
      write(98,*) rlat,rlong
      call recpol(xmin,ymax,rlat,rlong)
      print *,rlat,rlong
      write(98,*) rlat,rlong
      amin=1e30
      amax=-1e30
      dmin=1e30
      dmax=-1e30
      areatot=0
      do i = 1, numel
        do j=1,4
          lm(j)=kx(i,j)
        enddo
        lm(5)=lm(1)
        do m=1,4
          xy(1,m)=x(lm(m))
          xy(2,m)=y(lm(m))
        enddo
C ..... USE FOLLOWING TO GENERATE CENTROID VALUES FOR GRADIENTS AND
C       OTHER MATERIAL PROPERTIES
        CALL FESHAPE(1,0D0,0D0,PSI,DPSI)
        CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
        area=4*abs(detj)
        areatot=areatot+area
        amin=min(amin,area)
        amax=max(amax,area)
        do j=1,4
          dist=sqrt( (x(lm(j+1))-x(lm(j)))**2 +
     &               (y(lm(j+1))-y(lm(j)))**2 )
          dmin=min(dmin,dist)
          dmax=max(dmax,dist)
        enddo
      enddo
      print *,'area-based grid spacing',real(sqrt(amin)*0.001),
     &                       real(sqrt(amax)*0.001)
      print *,'           grid spacing',real(dmin*0.001),
     &                                  real(dmax*0.001)
      print *,real(xmin*.001),real(xmax*.001),
     &        real(ymin*.001),real(ymax*.001)
      print *,'total area=',areatot/1e12
      print *,'resolution=',sqrt(areatot/numel)/1000
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
c     IF(DETJ.LE.0.0) THEN
c       WRITE(12,*) DETJ,((XY(MM,NN),NN=1,4),MM=1,2)
c       WRITE(*,*) 'IN DERIVE, BAD JACOBIAN...'
c       WRITE(*,*) DETJ,((XY(MM,NN),NN=1,4),MM=1,2)
c       STOP
c     ENDIF
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
c-------------------------------------------------
      SUBROUTINE FESHAPE(NTYPE,XI,ET,PSI,DPSI)
C ELEMENT SHAPE FUNCTIONS AND DERIVATIVES AT LOCAL COORDINATES (XI,ET)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PSI(4),DPSI(4,2)
      IF(NTYPE.EQ.1) THEN
        PSI(1)=.25d0*(1.d0-XI)*(1.d0-ET)
        PSI(2)=.25d0*(1.d0+XI)*(1.d0-ET)
        PSI(3)=.25d0*(1.d0+XI)*(1.d0+ET)
        PSI(4)=.25d0*(1.d0-XI)*(1.d0+ET)
        DPSI(1,2)=-.25d0*(1.d0-XI)
        DPSI(2,2)=-.25d0*(1.d0+XI)
        DPSI(3,2)=.25d0*(1.d0+XI)
        DPSI(4,2)=.25d0*(1.d0-XI)
        DPSI(1,1)=-.25d0*(1.d0-ET)
        DPSI(2,1)=.25d0*(1.d0-ET)
        DPSI(3,1)=.25d0*(1.d0+ET)
        DPSI(4,1)=-.25d0*(1.d0+ET)
      ELSE
        PSI(1)=1.d0-XI-ET
        PSI(2)=XI
        PSI(3)=ET
        DPSI(1,2)=-1.d0
        DPSI(2,2)=0.d0
        DPSI(3,2)=1.d0
        DPSI(1,1)=-1.d0
        DPSI(2,1)=1.d0
        DPSI(3,1)=0.d0
      ENDIF
      END
