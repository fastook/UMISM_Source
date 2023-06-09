c latlong*.data file viewer
      PARAMETER(NMAX=29999)
      DIMENSION RLAT(NMAX),RLONG(NMAX),SURF(NMAX),BED(NMAX)
      DIMENSION XOUT(NMAX),YOUT(NMAX)
      CHARACTER JUNK1*80,JUNK2*80,JUNK3*80,ICE(NMAX)*3
      ncolors=10
      READ(1,100) JUNK1
      READ(1,100) JUNK2
      READ(1,100) JUNK3
100   FORMAT(A80)
      XMIN=1E30
      XMAX=-XMIN
      YMIN=XMIN
      YMAX=XMAX
      ZMIN=XMIN
      ZMAX=XMAX
      DO I=1,NMAX
        READ(1,200,END=998) N,RLAT(I),RLONG(I),SURF(I),BED(I),ICE(I)
200     FORMAT(I5,F11.4,F11.4,F11.4,F11.4,1X,A3)
        XMAX=MAX(XMAX,RLONG(I))
        YMAX=MAX(YMAX,RLAT(I))
        ZMAX=MAX(ZMAX,SURF(I))
        XMIN=MIN(XMIN,RLONG(I))
        YMIN=MIN(YMIN,RLAT(I))
        ZMIN=MIN(ZMIN,SURF(I))
        NPTS=I
      ENDDO
998   CONTINUE
      ZDELT=(ZMAX-ZMIN)/NCOLORS
      print *,xmin,xmax,ymin,ymax,zmin,zmax,zdelt
      XBORD=(XMAX-XMIN)/20.
      YBORD=(YMAX-YMIN)/20.
      XMAX=XMAX+XBORD
      XMIN=XMIN-XBORD
      YMAX=YMAX+YBORD
      YMIN=YMIN-YBORD
      CALL GRSTRT(800,800)
      CALL WINDOW(XMIN,XMAX,YMIN,YMAX)
      DO I=1,NPTS
        IF(ICE(I).EQ.'ICE' .OR. ICE(I).EQ.'ice') THEN
          ic=izset(ncolors,surf(i),zmin,zdelt,0)
          call linclr(ic)
          CALL MOVE(RLONG(I),RLAT(I))
          CALL DRAW(RLONG(I),RLAT(I))
        else
          call linclr(0)
          CALL MOVE(RLONG(I),RLAT(I))
          CALL DRAW(RLONG(I),RLAT(I))
        ENDIF
      ENDDO
      READ(2,100) JUNK1
      READ(2,100) JUNK2
      READ(2,100) JUNK3
      READ(2,*,END=999) XOUT(1),YOUT(1),isave
      call linclr(1)
      CALL MOVE(XOUT(1),YOUT(1))
      DO I=1,NMAX
        READ(2,*,END=999) XOUT(I),YOUT(I),iout
        if(xout(i).lt.0.) xout(i)=xout(i)+360.
        if(isave.eq.0) then
          CALL DRAW(XOUT(I),YOUT(I))
        else
          CALL move(XOUT(I),YOUT(I))
        endif
        isave=iout
      ENDDO
999   continue
      CALL GRSTOP
      END
      
      FUNCTION IZSET(NCOLOR,ZZZ,ZF,ZD,IOFF)
      IMPLICIT REAL*4(A-H,O-Z)                                          
      IMAX=NCOLOR+IOFF
      IMIN=IOFF+1
        ZIC=(ZZZ-ZF)/ZD+.5
        IC=NINT(ZIC)
      IF(IC.LT.IMIN) IC=IMIN
      IF(IC.GT.IMAX) IC=IMAX
      IZSET=IC
      RETURN
      END

