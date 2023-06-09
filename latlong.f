C READS FROM BINARY DATA SET, OUTPUTS:
C NODE,LATITUDE,LONGITUDE,SURFACE ELEVATION, AND THICKNESS IN ASCII
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
C UNIT 1 >= LATLONG&3 DATA B  (SINGLE DATA SET IN ASCII
C           NODE,LAT,LONG,ELEV,THICK
C UNITS 31,32,33 ARE BINARY, UNIT 
C RUN BY EXEC LATLONG.E:
C rm fort.*
C ln -s latlong$3.data fort.1
C ln -s input$1.head fort.30
C ln -s input$1.grid fort.31
C ln -s input$1.diff fort.32
C ln -s out$2.time fort.33
C latlong.x
C
      PARAMETER(NMAX=29999)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER HED*80,HEDT*80,JUNK*80
      DIMENSION KODE(NMAX),X(NMAX),Y(NMAX),HTICE(NMAX),ADOT(NMAX),
     &FRACT(NMAX),PSURF(NMAX),BDROCK(NMAX),FLOWA(NMAX),SLDGB(NMAX),
     &KX(NMAX,4),CONST(NMAX),IBFLUX(NMAX,2),BFLUX(NMAX),ACON(NMAX)
      print *,'impose offset for old grids (subtract 5 degree from'
      print *,'               longitudes'
      print *,'1 for offset, 0 for none'
      read(*,*) ioff
      print *,'input 0 for northern hemisphere (READN1 from global map)'
      print *,'      1 for southern hemisphere (GSMO5 from budds map)'
      read(*,*) ihem
      RHOI=0.917D0
      RHOW=1.092D0
      RATIO=RHOW/RHOI
      CALL SETRIG
C READ INPUT HEADER
      READ(30,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      PRINT *,HED
1000  FORMAT (A80,/,7I6,F8.0)
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
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME
85    READ(33,END=999) HEDT
      PRINT *,HEDT
      WRITE(JUNK,'(A80)') HEDT
      READ(JUNK,'(A6,F20.0)') HEDT,TIME
      READ(33) (HTICE(I),I=1,NUMNP)
      READ(33) (ADOT(I),I=1,NUMNP)
      READ(33) (BDROCK(I),I=1,NUMNP)
      READ(33) (CONST(I),I=1,NUMEL)
      READ(33) (ACON(I),I=1,NUMEL)
      PRINT *,'SHOULD I PRINT THIS? 1=YES, 0=NO, -1 TO QUIT'
      READ(*,*) IPRINT
      IF(IPRINT.LT.0) STOP
c     iprint=1
      IF(IPRINT.EQ.1) THEN
        WRITE(1,1010) HED,TIME
        DO 100 N=1,NUMNP
          IF(BDROCK(N).GT.0.) THEN
            FLOT=BDROCK(N)
          ELSE
            FLOT=(1.D0-RHOW/RHOI)*(BDROCK(N)-SEALEV)
          ENDIF
          CALL RECPOL(X(N)*.001D0,Y(N)*.001D0,RLAT,RLONG)
          if(ihem.eq.0) then
	      if(ioff.eq.1) then
		rlong=rlong-5.
		if(rlong.lt.0.) rlong=rlong+360.
	      endif          
          else
	      if(ioff.eq.1) then
		rlong=rlong-5.
		if(rlong.lt.0.) rlong=rlong+360.
	      endif
              rlong=-rlong+90.
5             if(rlong.gt.360.) then
                rlong=rlong-360.
                goto 5
              endif
10            if(rlong.lt.0.) then
                rlong=rlong+360.
                goto 10
              endif
          endif
          if(ihem.eq.0) then
            IF(HTICE(N).GT.FLOT+.1d0) THEN
              WRITE(1,1020) N,RLAT,RLONG,HTICE(N),BDROCK(N)
            ELSE
              WRITE(1,1030) N,RLAT,RLONG,HTICE(N),BDROCK(N)
            ENDIF
          else
            IF(HTICE(N).GT.FLOT+.1d0) THEN
              WRITE(1,1020) N,-RLAT,RLONG,HTICE(N),BDROCK(N)
            ELSE
              WRITE(1,1030) N,-RLAT,RLONG,HTICE(N),BDROCK(N)
            ENDIF
          endif
100     CONTINUE
      ENDIF
      GOTO 85
999   CONTINUE
1010  FORMAT(A80,/,
     &       ' AT MODEL TIME ',F9.1,/,
     &       ' node, latitude,  longitude, elevation,   bedrock')
1020  FORMAT(I6,4f11.4,' ice')
1030  FORMAT(I6,4f11.4)
      END
