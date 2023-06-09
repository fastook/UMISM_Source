C **************************************************************C       VEL00010
C                                                               C       VEL00020
C   PROGRAM:  VELO1                                             C       VEL00030
C                                                               C       VEL00040
C   DATE:  11 23 87                                             C       VEL00050
C   PROGRAMMER:  FASTOOK                                        C       VEL00060
C                                                               C       VEL00070
C   FUNCTION:                                                   C       VEL00080
C            PLOTS VELOCITY VECTORS COLOR-CODED FOR MAGNITUDE   C       VEL00090
C            FROM THE STANDARD OLD FORMAT OUTPUT SET OUT3**     C       VEL00100
C            AFTER ADJUSTING SCALE (PAIR, FIRST LENGHT OF       C       VEL00110
C            ARROW, SECOND MAX VELOCITY) OUTPUTS STANDARD       C       VEL00120
C            PLOTTER DATA SET TO OUTV** DATA B                  C       VEL00130
C                                                               C       VEL00140
C **************************************************************C       VEL00150
      PARAMETER(MXX=29999,HMIN=10.)                                     VEL00160
      CHARACTER HED*80                                                  VEL00170
      DIMENSION ICMAP(16)                                               VEL00180
      DIMENSION A(80),B(80),XARO(9),YARO(9),NBOUND(MXX)                 VEL00190
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2),LM(4)             VEL00200
      DIMENSION PSI(4),DPSI(4,2),CNST(100)                              VEL00210
      DIMENSION XY(2,4),XI(9),ETA(9),W(9)                               VEL00220
      DIMENSION XX(MXX),YY(MXX),ZZ(MXX),THICK(MXX),NTYPE(MXX)           VEL00230
      DIMENSION XLINE(MXX),YLINE(MXX)                                   VEL00240
      DIMENSION KX(MXX,4),IK(5),CONST(MXX)                              VEL00250
      DATA ICMAP /0,4,12,6,13,2,8,7,9,3,10,5,11,1,14,15/                VEL00260
C     DATA ICMAP /0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/                VEL00270
      ASCAL=10.                                                         VEL00280
      IARO=0                                                            VEL00290
      WRITE(*,*) 'INPUT SCALE FACTOR FOR ARROWS, UMAX, AND THRESHOLD'   VEL00300
      READ(*,*) ASCAL,UMAX,THRESH                                       VEL00310
      UDELTA=UMAX/13.                                                   VEL00320
      XLEN=5.0                                                          VEL00330
      YLEN=5.0                                                          VEL00340
      CALL GRSTRT(500,1)                                                VEL00350
      RNINE=-99999.                                                     VEL00360
      RMINUS=-1.                                                        VEL00370
      RZERO=0.                                                          VEL00380
      RTWO=2.                                                           VEL00390
      RAROW=20.                                                         VEL00400
      ICOLOR=1                                                          VEL00410
      ICSAV=ICOLOR                                                      VEL00420
 1    READ(1,100,END=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,     VEL00430
     &INTER,DT                                                          VEL00440
C **** THIS IS  PATCH ****                                              VEL00450
C     NUMGBC=0                                                          VEL00460
100   FORMAT(A80,/,7I6,E15.6)                                           VEL00470
      DO 10 NUM=1,NUMNP                                                 VEL00480
      READ(1,200) N,KODE,XX(NUM),YY(NUM),HTICE,                         VEL00490
     &ADOT,FRACT,DENS,BED,FLOWA,SLDGB                                   VEL00500
      ZZ(NUM)=HTICE                                                     VEL00510
      THICK(NUM)=HTICE-BED                                              VEL00520
10    CONTINUE                                                          VEL00530
200   FORMAT(I6,I4,2E12.5,F10.2,F7.2,F9.2,F10.3,F10.1,F10.2,F10.3)        VEL00540
      DO 12 I=1,NUMEL                                                   VEL00550
      READ(1,300) NUM, (KX(NUM,II),II=1,4),CONST(NUM)                   VEL00560
      NTYPE(NUM)=1                                                      VEL00570
      IF(KX(NUM,4).EQ.0) NTYPE(NUM)=2                                   VEL00580
12    CONTINUE                                                          VEL00590
      DO 13 N=1,NUMGBC                                                  VEL00600
13    READ(1,310) I,J,RJUNK                                             VEL00610
310   FORMAT(2I6,E13.6)                                                 VEL00620
300   FORMAT(5I6,1PE17.10)                                              VEL00630
      NLINE=0                                                           VEL00640
      READ(11,*,END=99) NLINE                                           VEL00650
      READ(11,*) (NBOUND(I),I=1,NLINE)                                  VEL00660
      REWIND 11                                                         VEL00670
99    CONTINUE                                                          VEL00680
      XMAX=-1.E30                                                       VEL00690
      XMIN=1.E30                                                        VEL00700
      YMAX=-1.E30                                                       VEL00710
      YMIN=1.E30                                                        VEL00720
      IF(NLINE.GT.0) THEN                                               VEL00730
        DO 15 I=1,NLINE                                                 VEL00740
        XLINE(I)=XX(NBOUND(I))/1000.                                    VEL00750
        YLINE(I)=YY(NBOUND(I))/1000.                                    VEL00760
        XMAX=MAX(XMAX,XLINE(I))                                         VEL00770
        XMIN=MIN(XMIN,XLINE(I))                                         VEL00780
        YMAX=MAX(YMAX,YLINE(I))                                         VEL00790
        YMIN=MIN(YMIN,YLINE(I))                                         VEL00800
        WRITE(7,105) XX(NBOUND(I))/1000.,YY(NBOUND(I))/1000.,I          VEL00810
15      CONTINUE                                                        VEL00820
105     FORMAT(10X,G13.6,2X,G13.6,I13)                                  VEL00830
        WRITE(7,105) RNINE,RTWO                                         VEL00840
        WRITE(7,111)                                                    VEL00850
111     FORMAT('BOUNDARY')                                              VEL00860
      ELSE                                                              VEL00870
        DO 151 I=1,NUMNP                                                VEL00880
        XMAX=MAX(XMAX,XX(I)*.001)                                       VEL00890
        XMIN=MIN(XMIN,XX(I)*.001)                                       VEL00900
        YMAX=MAX(YMAX,YY(I)*.001)                                       VEL00910
        YMIN=MIN(YMIN,YY(I)*.001)                                       VEL00920
151     CONTINUE                                                        VEL00930
      ENDIF                                                             VEL00940
      XBORD=(XMAX-XMIN)/10.                                             VEL00950
      XMIN=XMIN-XBORD                                                   VEL00960
      XMAX=XMAX+XBORD                                                   VEL00970
      YBORD=(YMAX-YMIN)/10.                                             VEL00980
      YMIN=YMIN-YBORD                                                   VEL00990
      YMAX=YMAX+YBORD                                                   VEL01000
23    CONTINUE                                                          VEL01010
      CALL WINDOW(XMIN,XMAX,YMIN,YMAX)                                  VEL01020
C     REWIND 9                                                          VEL01030
      WRITE(9,*) 'ELEMENT,X,Y,UX,UY,UMAG'                               VEL01040
      CALL LINCLR(1)                                                    VEL01050
      IF(NLINE.GT.0) THEN                                               VEL01060
        CALL MOVE(XLINE(1),YLINE(1))                                    VEL01070
        DO 73 I=2,NLINE                                                 VEL01080
          CALL DRAW(XLINE(I),YLINE(I))                                  VEL01090
73      CONTINUE                                                        VEL01100
      ENDIF                                                             VEL01110
      VMAX=-1.E30                                                       VEL01120
  557 DO 600 J=1,NUMEL                                                  VEL01130
      IF(NTYPE(J).EQ.1) THEN                                            VEL01140
        NNODE=4                                                         VEL01150
        CENTX=0.0D00                                                    VEL01160
        CENTY=0.0D00                                                    VEL01170
      ELSE                                                              VEL01180
        NNODE=3                                                         VEL01190
        CENTX=1.D0/3.D0                                                 VEL01200
        CENTY=1.D0/3.D0                                                 VEL01210
      ENDIF                                                             VEL01220
      HH=0.0D00                                                         VEL01230
      SUMX=0.0D00                                                       VEL01240
      SUMY=0.0D00                                                       VEL01250
      DO 560 I=1,NNODE                                                  VEL01260
  560 LM(I)=KX(J,I)                                                     VEL01270
      I=LM(1)                                                           VEL01280
      JJ=LM(2)                                                          VEL01290
      K=LM(3)                                                           VEL01300
      L=LM(4)                                                           VEL01310
      XY(1,1)=XX(I)                                                     VEL01320
      XY(1,2)=XX(JJ)                                                    VEL01330
      XY(1,3)=XX(K)                                                     VEL01340
      IF(NTYPE(J).EQ.1) XY(1,4)=XX(L)                                   VEL01350
      XY(2,1)=YY(I)                                                     VEL01360
      XY(2,2)=YY(JJ)                                                    VEL01370
      XY(2,3)=YY(K)                                                     VEL01380
      IF(NTYPE(J).EQ.1) XY(2,4)=YY(L)                                   VEL01390
      IF(NTYPE(J).EQ.1) THEN                                            VEL01400
        XCENT=(XX(I)+XX(JJ)+XX(K)+XX(L))/4000.                          VEL01410
        YCENT=(YY(I)+YY(JJ)+YY(K)+YY(L))/4000.                          VEL01420
      ELSE                                                              VEL01430
        XCENT=(XX(I)+XX(JJ)+XX(K))/3000.                                VEL01440
        YCENT=(YY(I)+YY(JJ)+YY(K))/3000.                                VEL01450
      ENDIF                                                             VEL01460
C                                                                       VEL01470
      CALL SHAPE(NTYPE(J),CENTX,CENTY,PSI,DPSI)                         VEL01480
      
C CALCULATE DXDS...EQUATION (5.3.6)                                     VEL01490
      DO 565 I=1,2                                                      VEL01500
      DO 565 L=1,2                                                      VEL01510
      DXDS(I,L)=0.0                                                     VEL01520
      DO 565 K=1,NNODE                                                  VEL01530
 565  DXDS(I,L)=DXDS(I,L)+DPSI(K,L)*XY(I,K)                             VEL01540
C CALCULATE DSDX...EQUATION (5.2.7)                                     VEL01550
      DETJ=(DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1))                    VEL01560
      IF (DETJ.LE.0.0) GOTO 750                                         VEL01570
      DSDX(1,1)=DXDS(2,2)/DETJ                                          VEL01580
      DSDX(2,2)=DXDS(1,1)/DETJ                                          VEL01590
      DSDX(1,2)=-DXDS(1,2)/DETJ                                         VEL01600
      DSDX(2,1)=-DXDS(2,1)/DETJ                                         VEL01610
C CALCULATE D(PSI)/DX...EQUATION (5.3.5)                                VEL01620
      DO 570 I=1,NNODE                                                  VEL01630
      DPSIX(I)=DPSI(I,1)*DSDX(1,1)+DPSI(I,2)*DSDX(2,1)                  VEL01640
      DPSIY(I)=DPSI(I,1)*DSDX(1,2)+DPSI(I,2)*DSDX(2,2)                  VEL01650
 570  CONTINUE                                                          VEL01660
      DO 580 I=1,NNODE                                                  VEL01670
      SUMX=SUMX + ZZ(LM(I))*DPSIX(I)                                    VEL01680
      SUMY=SUMY + ZZ(LM(I))*DPSIY(I)                                    VEL01690
      HH=HH + THICK(LM(I))*PSI(I)                                       VEL01700
  580 CONTINUE                                                          VEL01710
C                                                                       VEL01720
      DELH=SUMX**2 + SUMY**2                                            VEL01730
      DELH=SQRT(DELH)                                                   VEL01740
      IF(HH.GT.HMIN) THEN                                               VEL01750
        UX=-CONST(J)*SUMX/HH                                            VEL01760
        UY=-CONST(J)*SUMY/HH                                            VEL01770
      ELSE                                                              VEL01800
        UX=0.                                                           VEL01810
        UY=0.                                                           VEL01820
      ENDIF                                                             VEL01830
      UMAG=(UX**2+UY**2)**.5                                            VEL01840
      VMAX=MAX(VMAX,UMAG)                                               VEL01850
2222  format(4g13.6)
      IF(UMAG.GT.THRESH) GOTO 600                                       VEL01860
      IF(UMAG.EQ.0.) THEN                                               VEL01870
        ICOLOR=1                                                        VEL01880
      ELSE                                                              VEL01890
        ICOLOR=1+NINT((UMAG)/UDELTA)                                    VEL01900
      ENDIF                                                             VEL01910
      IF(ICOLOR.LT.2) ICOLOR=2                                          VEL01920
      IF(ICOLOR.GT.14) ICOLOR=14                                        VEL01930
C     ICOLOR=ICMAP(ICOLOR)                                              VEL01940
      WRITE(9,1001) J,XCENT,YCENT,UX,UY,UMAG,CONST(J)
1001  FORMAT(I6,6G13.6)                                                 VEL01960
      XARO(1)=XCENT                                                     VEL01970
      XARO(2)=XCENT+ASCAL*UX                                            VEL01980
      YARO(1)=YCENT                                                     VEL01990
      YARO(2)=YCENT+ASCAL*UY                                            VEL02000
202   FORMAT(10X,G13.6,2X,G13.6,I13)                                    VEL02010
      WRITE(7,202) XARO(1),YARO(1)                                      VEL02020
      WRITE(7,202) XARO(2),YARO(2)                                      VEL02030
      WRITE(7,202) RNINE,RAROW,ICMAP(ICOLOR)-1                          VEL02040
      WRITE(7,*) ' ELEMENT',J                                           VEL02050
      IF(ICOLOR.NE.ICSAV) THEN                                          VEL02060
        CALL LINCLR(ICMAP(ICOLOR))                                      VEL02070
        ICSAV=ICOLOR                                                    VEL02080
      ENDIF                                                             VEL02090
      CALL MOVE(XARO(1),YARO(1))                                        VEL02100
      CALL DRAW(XARO(2),YARO(2))                                        VEL02110
  600 CONTINUE                                                          VEL02120
      XARO(1)=XMIN+(XMAX-XMIN)/XLEN                                     VEL02130
      XARO(2)=XARO(1)+(XMAX-XMIN)/XLEN                                  VEL02140
      YARO(1)=YMAX+(YMAX-YMIN)/YLEN                                     VEL02150
      YARO(2)=YARO(1)                                                   VEL02160
      CALL MOVE(XARO(1),YARO(1))                                        VEL02170
      CALL DRAW(XARO(2),YARO(2))                                        VEL02180
      ZSCAL=(XMAX-XMIN)/ASCAL/XLEN                                      VEL02190
C     WRITE(7,202) XARO(1),YARO(1)                                      VEL02200
C     WRITE(7,202) XARO(2),YARO(2)                                      VEL02210
      WRITE(7,202) XMIN,YMIN                                            VEL02220
      WRITE(7,202) XMIN,YMIN                                            VEL02230
      WRITE(7,202) RNINE,RAROW+1.                                       VEL02240
      WRITE(7,2323) ZSCAL                                               VEL02250
2323  FORMAT(G10.3)                                                     VEL02260
      CALL GRSTOP                                                       VEL02270
      WRITE(*,*) 'UMAX=',VMAX                                           VEL02280
      READ(*,*,END=24) ASCAL1,UMAX,THRESH                               VEL02290
      UDELTA=UMAX/13.                                                   VEL02300
      IF(ASCAL1.EQ.0.) GOTO 24                                          VEL02310
      ASCAL=ASCAL1                                                      VEL02320
      REWIND 7                                                          VEL02330
      CALL GRSTRT(500,1)                                                VEL02340
      GOTO 23                                                           VEL02350
24    CONTINUE                                                          VEL02360
      GOTO 1                                                            VEL02370
999   CONTINUE                                                          VEL02380
      CALL GRSTOP                                                       VEL02390
      STOP                                                              VEL02400
 750  WRITE (12,5544) DETJ,X                                            VEL02410
      WRITE (*,5544) DETJ,X                                             VEL02420
 5544 FORMAT (' BAD JACOBIAN at 750',E14.6,E14.6)                       VEL02430
      STOP                                                              VEL02440
      END                                                               VEL02450
      SUBROUTINE SHAPE(NTYPE,XI,ET,PSI,DPSI)                            VEL02460
C     IMPLICIT REAL*8(A-H,O-Z)                                          VEL02470
      DIMENSION PSI(4),DPSI(4,2)                                        VEL02480
      IF(NTYPE.EQ.2) GOTO 100                                           VEL02490
      PSI(1)=.25*(1.-XI)*(1.-ET)                                        VEL02500
      PSI(2)=.25*(1.+XI)*(1.-ET)                                        VEL02510
      PSI(3)=.25*(1.+XI)*(1.+ET)                                        VEL02520
      PSI(4)=.25*(1.-XI)*(1.+ET)                                        VEL02530
      DPSI(1,2)=-.25*(1.-XI)                                            VEL02540
      DPSI(2,2)=-.25*(1.+XI)                                            VEL02550
      DPSI(3,2)=.25*(1.+XI)                                             VEL02560
      DPSI(4,2)=.25*(1.-XI)                                             VEL02570
      DPSI(1,1)=-.25*(1.-ET)                                            VEL02580
      DPSI(2,1)=.25*(1.-ET)                                             VEL02590
      DPSI(3,1)=.25*(1.+ET)                                             VEL02600
      DPSI(4,1)=-.25*(1.+ET)                                            VEL02610
      RETURN                                                            VEL02620
100   CONTINUE                                                          VEL02630
      PSI(1)=1.-XI-ET                                                   VEL02640
      PSI(2)=XI                                                         VEL02650
      PSI(3)=ET                                                         VEL02660
      DPSI(1,2)=-1.                                                     VEL02670
      DPSI(2,2)=0.                                                      VEL02680
      DPSI(3,2)=1.                                                      VEL02690
      DPSI(1,1)=-1.                                                     VEL02700
      DPSI(2,1)=1.                                                      VEL02710
      DPSI(3,1)=0.                                                      VEL02720
      RETURN                                                            VEL02730
      END                                                               VEL02740
