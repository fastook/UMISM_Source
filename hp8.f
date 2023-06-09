      PARAMETER(MAX=6000,NUMTOT=30)
C **************************************************************C
C                                                               C
C   PROGRAM:  HPLOT6                                            C
C                                                               C
C   DATE:  12 13 87                                             C
C   PROGRAMMER:  FASTOOK                                        C
C                                                               C
C   FUNCTION:                                                   C
C            PLOTS CONTOURS (SURFACE, BED VERSUS X (OR Y WITH   C
C            MINOR CHANGE)) FROM FILE LINE* WHICH CONTAINS      C
C            NODE NUMBERS TO BE PLOTTED, OUTPUT O2* (SURF)      C
C            AND O3* (BED) ARE STANDARD FORM FOR X-Y PLOTTER    C
C                                                               C
C **************************************************************C
      CHARACTER*80 HED(numtot)
      DIMENSION PXARAY(10),PYARAY(10)
      DIMENSION ZZ(MAX,NUMTOT)
      DIMENSION X(MAX),Y(MAX)
      DIMENSION KX(MAX,5),INDEX(MAX),DIST(MAX)
      DIMENSION XX(MAX),YY(MAX),LM(5)
      DIMENSION ICMAP(125)
C FOLLOWING IS ORIGINAL COLOR SCHEME, WITH UP TO 23 COLORS
      DATA ICMAP /50,174,170,165,161,155,150,137,106,80,65,70,92,
     1            98,99,74,69,63,58,53,79,84,114,144,101*0/
      DATA NCOLOR/240/
C FOLLOWING IS ALTERNATE SCHEME, MORE MAGENTA, LESS DARKS, ALSO 23 COLOR
C     DATA ICMAP /0,-1,-7,165,161,155,-8,157,163,134,129,-6,104,109,
C    1            89,69,124,99,-5,97,-9,120,145,147,101*0/
C     DATA NCOLOR/23/
C FOLLOWING IS FOR 15 COLORS, IN TERMS OF COLOR INDICES 2-15. CAN BE
C REMAPPED INTO SEQUNTIALS WITH F4, SHFT F4
C     DATA ICMAP /0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14
C    1            ,-15,109*0/
C     DATA NCOLOR/15/
C FOLLOWING IS FOR GREY SCALE, FOR BLACK AND WHITE, 14 COLORS
C     DATA ICMAP /0,1,2,3,4,5,7,8,9,10,11,12,13,14,16,-14,-15,108*0/
C    DATA NCOLOR/14/
C TEST OF ALL GREEN, 24 COLORS
C     DATA ICMAP /0,-1,147,146,145,122,121,120,96,71,70,117,97,
C    1            95,72,115,90,91,66,65,86,85,61,60,55,100*0/
C     DATA NCOLOR/24/
      ISKIP=0
      IXO=0
      IYO=0
      IXTOT=3130-IXO
      IYTOT=3130-IYO
      ZLEN=NCOLOR
      XLEN=6.0
      YLEN=6.0
      XLEN=100.
      YLEN=100.
      XO=2.
      YO=0.
      XTOT=8.
      YTOT=6.
      IPASS=0
      HED(1)=' '
C     REWIND 1
      ISKIP=0
      IF(IPASS.EQ.0) THEN
        WRITE(*,1010) HED(1)
1010  FORMAT(1X,'0-ICE SURFACE',/,
     1 ' 1-BED',/,
     2 ' 2-THICKNESS',/,
     3 ' 3-VELOCITY (NOT WORKING)',/,
     4 ' 4-TAU (NOT WORKING)',/,
     5 ' 5-ADOT',/,
     6 ' 6-FLOWA',/,
     7 ' 7-SLDGB',/,
     8 ' 8-FRACT',/,
     9 ' 9-PSURF',/,
     9 ' 10-SKIP THIS ONE',/,
     9 A80)
        READ(*,*,END=1) IPLOT
        WRITE(*,*) 'INPUT RSCALE ANGLE AND 0 TO AUTO SCALE,1 TO DEFINE'
        READ(*,*) RSCALE,ANGLE,ISCALE
        IF(ISCALE.EQ.1) THEN
          WRITE(*,*) 'MINIMUM AND DELTA VALUES'
          READ(*,*) ZMINN,ZDELTT
        ENDIF
      ENDIF
      XMAX=-1.E30
      YMAX=-1.E30
      XMIN=1.E30
      YMIN=1.E30
       DO 900 NUMBER=1,NUMTOT
10    READ(1,1000,END=1) HED(NUMBER),NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &INTER,DT
1000  FORMAT(A80,/,7I6,E15.6)
      NUMREAD=NUMBER
      write(*,*) 'reading data'
      DO 300 I=1,NUMNP
      READ(1,1003) N,KODE,X(N),Y(N),HT,
     &ADOT,FRACT,PSURF,BED,FLOWA,SLDGB
      IF(IPLOT.EQ.0) ZZ(I,NUMBER)=HT
      IF(IPLOT.EQ.1) ZZ(I,NUMBER)=BED
      IF(IPLOT.EQ.2) ZZ(I,NUMBER)=HT-BED
      IF(IPLOT.EQ.5) ZZ(I,NUMBER)=ADOT
      IF(IPLOT.EQ.6) ZZ(I,NUMBER)=FLOWA
      IF(IPLOT.EQ.7) ZZ(I,NUMBER)=SLDGB
      IF(IPLOT.EQ.8) ZZ(I,NUMBER)=FRACT
      IF(IPLOT.EQ.9) ZZ(I,NUMBER)=PSURF
      XMAX=AMAX1(X(N),XMAX)
      XMIN=AMIN1(X(N),XMIN)
      YMAX=AMAX1(Y(N),YMAX)
      YMIN=AMIN1(Y(N),YMIN)
c 1003  FORMAT(I6,I4,2E12.5,F10.2,F7.2,F9.2,F10.3,F10.1,F10.2,F10.3)
1003  FORMAT(I6,I4,2E12.5,F10.2,F7.2,F9.2,F10.3,F10.1,F10.2,F10.3)
1001  FORMAT(I6,I4,F9.0,F9.3,F9.4,F7.2,1PE10.3,
     &0PF9.3,F7.4)
200   FORMAT(10X,G13.6,2X,G13.6,I13)
300   CONTINUE
      DO 600 I=1,NUMEL
      READ(1,6345) NUM,(KX(NUM,J),J=1,4)
      KX(NUM,5)=KX(NUM,1)
 6345 FORMAT (5I6)
  600 CONTINUE
      DO 610 I=1,NUMGBC
      READ(1,*) IJUNK,JJUNK
610   CONTINUE
      write(*,*) 'done reading data'
      write(*,*) hed(number)
900   CONTINUE
      write(*,*) 'all done with reading data'
1     CONTINUE
      do 9999 number=1,numread
      IF(IPASS.EQ.0) THEN
        CALL GRSTRT(800,-1)
        CALL WINDOW(0.,XLEN*1.2,0.,YLEN)
        IPASS=1
        RSCL=1./RSCALE
        CALL LINCLR(0)
        XORG=0.
        YORG=0.
        OFFSET=-2000.
        ICOLOR=1
        COSS=COS(3.14159*ANGLE/180.)
        SINS=SIN(3.14159*ANGLE/180.)
        RNINE=-99999.
        ICOLOR=1
        RZERO=0.
        RMINUS=-1.
        IF(ANGLE.GT.0. .AND. ANGLE.LE.90.) THEN
          XORG=XMIN
          YORG=YMIN
        ENDIF
        IF(ANGLE.GT.90. .AND. ANGLE.LE.180.) THEN
          XORG=XMAX
          YORG=YMIN
        ENDIF
        IF(ANGLE.GT.180. .AND. ANGLE.LE.270.) THEN
          XORG=XMIN
          YORG=YMAX
        ENDIF
        IF(ANGLE.GT.270. .AND. ANGLE.LE.360.) THEN
          XORG=XMAX
          YORG=YMAX
        ENDIF
        CALL SCALE2(ZZ(1,NUMBER),ZLEN,NUMNP,1)
        ZF=ZZ(NUMNP+1,NUMBER)
        ZD=ZZ(NUMNP+2,NUMBER)
C REMOVE NEXT 3 FOR AUTO SCALING
        IF(ISCALE.EQ.1) THEN
          ZF=ZMINN
C         ZMAX=1300.
          ZD=(ZMAX-ZF)/ZLEN
          ZD=ZDELTT
        ENDIF
      ENDIF
c     do 1235 lk=1,2
      if(lk.eq.2) call newpag
      DO 580 I=1,NUMEL
      L1=KX(I,1)
      L2=KX(I,2)
      L3=KX(I,3)
      L4=KX(I,4)
      XCENT=.25*(X(L1)+X(L2)+X(L3)+X(L4))
      YCENT=.25*(Y(L1)+Y(L2)+Y(L3)+Y(L4))
      DIST(I)=SQRT((XORG-XCENT)**2+(YORG-YCENT)**2)
580   CONTINUE
      CALL INDEXX(NUMEL,DIST,INDEX)
      DO 590 I=1,NUMEL
      DO 590 J=1,4
      LM(J)=KX(INDEX(NUMEL-I+1),J)
      XX(LM(J))=SIGN(1.,180.-ANGLE)*(X(LM(J))-Y(LM(J))*COSS)*RSCL
      YY(LM(J))=(RSCALE*ZZ(LM(J),NUMBER)+Y(LM(J))*SINS)*RSCL
      IQ=LM(J)
590   CONTINUE
      CALL SCALE1(XX,XLEN,NUMNP,1)
      CALL SCALE1(YY,YLEN,NUMNP,1)
      IF(XX(NUMNP+2).GT.YY(NUMNP+2)) THEN
        YY(NUMNP+2)=XX(NUMNP+2)
      ELSE
        XX(NUMNP+2)=YY(NUMNP+2)
      ENDIF
      XF=XX(NUMNP+1)
      XD=XX(NUMNP+2)
      YF=YY(NUMNP+1)
      YD=YY(NUMNP+2)
      DO 595 I=1,NUMEL
      ZZZ=0.
      DO 593 J=1,4
      LM(J)=KX(INDEX(NUMEL-I+1),J)
      ZZZ=ZZZ+ZZ(LM(J),NUMBER)
593   CONTINUE
      ZZZ=.25*ZZZ
      IZ=IZSET(NCOLOR,ZZZ,ZF,ZD,1)
          IPAT=ICMAP(IZ)
          IPAT=IZ
      CALL FILPAN(IPAT,.TRUE.)
      CALL MAKCUR
C     IF(IPAT.LE.50 .OR. IPAT.GE.174) WRITE(23,*) ZZZ,IPAT
          DO 201 K=1,4
C MUST SCALE PX,PY TO FIT IN WORLD COORDINATES
            PXARAY(K)=(XX(LM(K))-XF)/XD
            PYARAY(K)=(YY(LM(K))-YF)/YD
201       CONTINUE
      PXARAY(5)=PXARAY(1)
      PYARAY(5)=PYARAY(1)
      CALL PANEL(5,PXARAY,PYARAY)
595   CONTINUE
      CALL LINCLR(1)
      IPL=0
2323  READ(3,*,END=9998) XP,YP
      IF(XP.EQ.-99999.) GOTO 9998
      XP=SIGN(1.,180.-ANGLE)*(XP*1000.-YP*1000.*COSS)*RSCL
      YP=(RSCALE*0.+YP*1000.*SINS)*RSCL
      XP=(XP-XF)/XD
      YP=(YP-YF)/YD
      IF(IPL.EQ.0) THEN
        CALL MOVE(XP,YP)
      ELSE
        CALL DRAW(XP,YP)
      ENDIF
      IPL=1
      GOTO 2323
9998  CONTINUE
      CALL LINCLR(0)
C     ISKIP=0
c     IF(ISKIP.EQ.1) GOTO 1234
      ISKIP=1
C     CALL LABEL1(3030,2800,IPLOT,HED(NUMBER))
      CALL TXTCLR(0)
      CALL MOVE(10.,95.)
      CALL TEXT(80,HED(NUMBER))
      CALL MOVE(10.,92.5)
      IF(IPLOT.EQ.0) CALL TEXT(17,'SURFACE ELEVATION')
      IF(IPLOT.EQ.1) CALL TEXT(17,'BEDROCK ELEVATION')
      IF(IPLOT.EQ.2) CALL TEXT(17,'THICKNESS        ')
      IF(IPLOT.EQ.3) CALL TEXT(17,'VELOCITY         ')
      IF(IPLOT.EQ.4) CALL TEXT(17,'BASAL STRESS     ')
      IF(IPLOT.EQ.5) CALL TEXT(17,'ACCUMULATION RATE')
      IF(IPLOT.EQ.6) CALL TEXT(17,'FLOW LAW CONSTANT')
      IF(IPLOT.EQ.7) CALL TEXT(20,'SLIDING LAW CONSTANT')
      IF(IPLOT.EQ.8) CALL TEXT(17,'FRACTION MELTED  ')
      IF(IPLOT.EQ.9) CALL TEXT(17,'PRESENT SURFACE  ')
      XXX=105.
      YYY=92.5
      IYY=2500
      IXX=3130
      RNUM=ZF+ZLEN*ZD
      DO 345 I=2,NCOLOR+1
      CALL FILPAN(I,.FALSE.)
      PXARAY(1)=XXX
      PYARAY(1)=YYY
      PXARAY(2)=XXX
      PYARAY(2)=YYY+.3
      PXARAY(3)=XXX+4.
      PYARAY(3)=YYY+.3
      PXARAY(4)=XXX+4.
      PYARAY(4)=YYY
      CALL PANEL(4,PXARAY,PYARAY)
C     CALL WRNUM(IXX+150,IYY+75,RNUM,0)
      CALL MOVE(XXX+4.,YYY+.35)
      CALL TXTCLR(0)
      CALL RNUMBR(RNUM,0,6)
      RNUM=RNUM-ZD
      IYY=IYY-100
      YYY=YYY-.3
345   CONTINUE
      RNUM=ZF
C     CALL WRNUM(IXX+150,IYY+75,RNUM,0)
      CALL MOVE(XXX+6.,YYY+3.)
      CALL RNUMBR(RNUM,0,6)
1234  CONTINUE
      CALL MAKCUR
1235  continue
9999  CONTINUE
      CALL GRSTOP
      write(*,*) 'go around again (1), quit (0)'
      read(*,*) igo
      ipass=0
      if(igo.eq.1) goto 1
      STOP
      END
      SUBROUTINE INDEXX(N,ARRIN,INDX)
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1) THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1) THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR) THEN
          IF(J.LT.IR) THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1))) J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J))) THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GOTO 20
        ENDIF
        INDX(I)=INDXT
      GOTO 10
      END
      FUNCTION IZSET(NCOLOR,ZZZ,ZF,ZD,IOFF)
      IMAX=NCOLOR+IOFF
      IMIN=IOFF+1
        IC=IMAX-NINT((ZZZ-ZF)/ZD-.5)
      IF(IC.LT.IMIN) IC=IMIN
      IF(IC.GT.IMAX) IC=IMAX
      IZSET=IC
      RETURN
      END
      SUBROUTINE SCALE1(AR,RR,NN,II)
      DIMENSION AR(NN)
      AMIN=1.E30
      AMAX=-1.E30
      DO 100 I=1,NN
      AMAX=AMAX1(AR(I),AMAX)
      AMIN=AMIN1(AR(I),AMIN)
100   CONTINUE
      IF(AMIN.EQ.AMAX) AMAX=AMIN+1.
      AR(NN+1)=AMIN
      AR(NN+2)=(AMAX-AMIN)/RR
      AMAX=AMAX+AR(NN+2)
      AMIN=AMIN-AR(NN+2)
      AR(NN+1)=AMIN
      AR(NN+2)=(AMAX-AMIN)/RR
      RETURN
      END
      SUBROUTINE SCALE2(AR,RR,NN,II)
      DIMENSION AR(NN)
      AMIN=1.E30
      AMAX=-1.E30
      DO 100 I=1,NN
      AMAX=AMAX1(AR(I),AMAX)
      AMIN=AMIN1(AR(I),AMIN)
100   CONTINUE
      IF(AMIN.EQ.AMAX) AMAX=AMIN+1.
      AR(NN+1)=AMIN
      AR(NN+2)=(AMAX-AMIN)/RR
      RETURN
      END
