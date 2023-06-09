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
     &          ACON(NMAX),GEOFLUX(NMAX),calv(nmax)
      dimension lm(5)
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
     &                   ADOT(N),FRACT(N),PSURF(N),
     &                   BDROCK(N),FLOWA(N),SLDGB(N),
     &                   temp(n),itype(n),AFUDGE(N),GEOFLUX(N),
     &                   calv(n)
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
      WRITE(33) (ADOT(I),I=1,NUMNP)
      WRITE(33) (BDROCK(I),I=1,NUMNP)
      WRITE(33) (CONST(I),I=1,NUMEL)
      WRITE(33) (ACON(I),I=1,NUMEL)
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

      dmin=1e30
      dmax=-1e30
      do i = 1, numel
        do j=1,4
          lm(j)=kx(i,j)
        enddo
        lm(5)=lm(1)
        do j=1,4
          dist=sqrt( (x(lm(j+1))-x(lm(j)))**2 +
     &               (y(lm(j+1))-y(lm(j)))**2 )
          dmin=min(dmin,dist)
          dmax=max(dmax,dist)
        enddo
      enddo
      print *,'grid spacing',real(dmin*0.001),real(dmax*0.001)
      print *,real(xmin*.001),real(xmax*.001),
     &        real(ymin*.001),real(ymax*.001)
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
