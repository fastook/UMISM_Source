      IMPLICIT REAL*8(A-H,O-Z)
C **************************************************************C
C                                                               C
C   PROGRAM:  CONT3                                             C
C                                                               C
C   DATE:  11 23 87                                             C
C   PROGRAMMER:  FASTOOK                                        C
C                                                               C
C   FUNCTION:                                                   C
C            CONTOURS STAANDARD DATA SET (OUT3**) FOR SURFACE,  C
C            BED, FRACT, FLOWA, THICK ETC. OUTPUT IN STANDARD   C
C            PLOTTER FORM TO OUTC**                             C
C                                                               C
C **************************************************************C
      CHARACTER HED*80
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      CALL SETRIG
      WRITE(*,*) 'INPUT 0 FOR NORTHERN HEMISPHERE, 1 FOR SOUTH'
      READ(*,*) IHEMI
      WRITE(*,*) 'INPUT 0 FOR LAT-LONG, 1 FOR X-Y'
      READ(*,*) IXY
C
      READ(1,100,END=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &                    INTER,DT
      PRINT *,HED
100   FORMAT(A80,/,7I5,E15.6)
      NUM1=NUMCOL-1
      OFFSET=5./60.
      HMIN=1.E30
      HMAX=-HMIN
      AMIN=1.E30
      AMAX=-HMIN
      FRMIN=1.E30
      FRMAX=-HMIN
      PMIN=1.E30
      PMAX=-HMIN
      BMIN=1.E30
      BMAX=-HMIN
      FLMIN=1.E30
      FLMAX=-HMIN
      TMIN=1.E30
      TMAX=-HMIN
      THMIN=1.E30
      THMAX=-HMIN
      DMIN=1.E30
      DMAX=-HMIN
      DO NUM=1,NUMNP
        READ(1,200) N,KODE,XX,YY,HTICE,
     &              ADOT,FRACT,PSURF,BDROCK,FLOWA,SLDGB,
     &              TBED
        XX=XX*.001
        YY=YY*.001
        THICK=HTICE-BDROCK
c        IF(BDROCK.LT.0.) THEN
c          FLOT=-BDROCK*1.03/.917
c          FLOT=BDROCK*(1.-1.03/.917)
c          IF(HTICE.LT.FLOT) THICK=0.
c        ENDIF
        DIFF=HTICE-PSURF
C
        HMIN=MIN(HMIN,HTICE)
        HMAX=MAX(HMAX,HTICE)
        AMIN=MIN(AMIN,ADOT)
        AMAX=MAX(AMAX,ADOT)
        FRMIN=MIN(FRMIN,FRACT)
        FRMAX=MAX(FRMAX,FRACT)
        PMIN=MIN(PMIN,PSURF)
        PMAX=MAX(PMAX,PSURF)
        BMIN=MIN(BMIN,BDROCK)
        BMAX=MAX(BMAX,BDROCK)
        FLMIN=MIN(FLMIN,FLOWA)
        FLMAX=MAX(FLMAX,FLOWA)
        SMIN=MIN(SMIN,SLDGB)
        SMAX=MAX(SMAX,SLDGB)
        TMIN=MIN(TMIN,TBED)
        TMAX=MAX(TMAX,TBED)
        THMIN=MIN(THMIN,THICK)
        THMAX=MAX(THMAX,THICK)
        DMIN=MIN(DMIN,DIFF)
        DMAX=MAX(DMAX,DIFF)
        CALL RECPOL(XX,YY,RLAT,RLONG)
        IF(IHEMI.EQ.1) THEN
          RLAT=-RLAT
          RLONG=90-RLONG
          IF(RLONG.LT.0) RLONG=RLONG+360
        ENDIF
        RLONG=RLONG+OFFSET
        IF(IXY.EQ.0) THEN
          WRITE(10,*) RLONG,RLAT,HTICE
          WRITE(11,*) RLONG,RLAT,BDROCK
          WRITE(12,*) RLONG,RLAT,ADOT
          WRITE(13,*) RLONG,RLAT,FRACT
          WRITE(14,*) RLONG,RLAT,FLOWA
          WRITE(15,*) RLONG,RLAT,DIFF
          WRITE(16,*) RLONG,RLAT,THICK
          WRITE(17,*) RLONG,RLAT,SLDGB
          WRITE(18,*) RLONG,RLAT,TBED
        ELSE
          WRITE(20,*) XX,YY,HTICE
          WRITE(21,*) XX,YY,BDROCK
          WRITE(22,*) XX,YY,ADOT
          WRITE(23,*) XX,YY,FRACT
          WRITE(24,*) XX,YY,FLOWA
          WRITE(25,*) XX,YY,DIFF
          WRITE(26,*) XX,YY,THICK
          WRITE(27,*) XX,YY,SLDGB
          WRITE(28,*) XX,YY,TBED
        ENDIF
      ENDDO
      WRITE(30,101) HMIN,250,250,250,HMAX,50,50,50
      WRITE(31,101) BMIN,250,250,250,BMAX,50,50,50
      WRITE(32,101) AMIN,250,250,250,AMAX,50,50,50
      WRITE(33,101) FRMIN,250,250,250,FRMAX,50,50,50
      WRITE(34,101) FLMIN,250,250,250,FLMAX,50,50,50
      WRITE(35,101) DMIN,250,250,250,DMAX,50,50,50
      WRITE(36,101) THMIN,250,250,250,THMAX,50,50,50
      WRITE(37,101) SMIN,250,250,250,SMAX,50,50,50
      WRITE(38,101) TMIN,250,250,250,TMAX,50,50,50
101   FORMAT(1X,2(1PE13.6,3I4))
200   FORMAT(2I4,2E12.5,F10.2,F7.2,F9.2,F10.3,F10.1,F10.2,2F10.3)
999   continue
      END
c***************************************************

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
