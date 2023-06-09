C CONTOUR.F TAKES THE OUTPUT OF CONT3.F AND ASSEMBLES THE SEGMENTED
C CONTOURS INTO INDIVIDUAL LINES, FOR EITHER PLOTTING WITH PLOTB.F 
C OR MANIPULATION WITH MKPLOTC.F TO CREATE A .NAMELIST FILE THAT CAN
C BE USED WITH PLOTC TO GENERATE POSTCRIPT FILES FOR PRINTING ON BUDDY.
      PARAMETER(NMAX=1000,NLMAX=2000)
      DIMENSION X(4),Y(4),XX(NLMAX,NMAX),YY(NLMAX,NMAX),RLINE(NLMAX)
      DIMENSION NPOINTS(NLMAX),IUSED(NLMAX)
      CHARACTER*80 TITLE1,LINE1,LINE2,LINES(NLMAX)
      IFINISH=0
      IC=1
      LINE1=' '
      XMAX=-1E30
      XMIN=-XMAX
      YMAX=-1E30
      YMIN=-XMAX
      DO I=1,4
        READ(1,*) X(I),Y(I)
        WRITE(11,202) X(I),Y(I)
        XMAX=MAX(XMAX,X(I))
        XMIN=MIN(XMIN,X(I))
        YMAX=MAX(YMAX,Y(I))
        YMIN=MIN(YMIN,Y(I))
      ENDDO
      WRITE(11,202) -99999.,21.,0
      WRITE(11,*) 'OUTLINE'
      PRINT *,XMIN,XMAX,YMIN,YMAX
      CALL GRSTRT(800,800)
      CALL WINDOW(XMIN,XMAX,YMIN,YMAX)
      READ(1,100) TITLE1
      READ(1,100) TITLE1
100   FORMAT(A80)
      DO NCONT=1,NLMAX
        CALL LINCLR(NCONT)
        DO I=1,5
          READ(1,*,END=99) XP,YP
          IF(XP.EQ.-99999.) THEN
            READ(1,100) LINE1
            GOTO 180
          ELSE
            XX(1,I)=XP
            YY(1,I)=YP
            NPOINTS(1)=I
          ENDIF
        ENDDO
        PRINT *,'MORE THAN 4 IN FIRST SEGMENT'
        STOP
180     CONTINUE
        DO NSEG=2,NLMAX
          DO I=1,5
            READ(1,*,END=195) XP,YP
            IF(XP.EQ.-99999.) THEN
              READ(1,100) LINE2
              IF(LINE1.NE.LINE2) THEN
                RLINE(NCONT)=PARSE(LINE1)
                LINES(NCONT)=LINE1
                LINE1=LINE2
                NNSEG=NSEG
                DO K=1,I+1
                  BACKSPACE 1
                ENDDO
                IFINISH=0
                GOTO 200
              ENDIF
              GOTO 190
            ELSE
              XX(NSEG,I)=XP
              YY(NSEG,I)=YP
              NPOINTS(NSEG)=I
              CALL MOVE(XP,YP)
              CALL DRAW(XP,YP)
            ENDIF
          ENDDO 
          PRINT *,'MORE THAN 4 POINTS ON A CONTOUR'
          STOP
190       CONTINUE
          NNSEG=NSEG
        ENDDO
        PRINT *,'NLMAX TOO SMALL',NLMAX
        STOP
195     CONTINUE
        RLINE(NCONT)=PARSE(LINE2)
        LINES(NCONT)=LINE2
        PRINT *,'END OF FILE ENCOUNTERED',NCONT,NSEG,NNSEG
        IFINISH=1
200     CONTINUE
        NSEG=NNSEG
        PRINT *,NSEG,' IS THE NUMBER OF SEGMENTS FOR ',LINES(NCONT)
        DO ISEG=1,NSEG
          IUSED(ISEG)=0
        ENDDO
        DO ICSEG=1,NSEG
          NNSEG=0
          IF(IUSED(ICSEG).EQ.0) THEN
            IUSED(ICSEG)=1
            ICHECK=ICSEG
            JCHECK=0
210       CONTINUE
            DO ISEG=1,NSEG
              TOL=1E-03
              IF(IUSED(ISEG).EQ.0) THEN
                DX1=ABS(XX(ICHECK,NPOINTS(ICHECK))-XX(ISEG,1))
                DY1=ABS(YY(ICHECK,NPOINTS(ICHECK))-YY(ISEG,1))
                DX2=ABS(XX(ICHECK,NPOINTS(ICHECK))-
     &                  XX(ISEG,NPOINTS(ISEG)))
                DY2=ABS(YY(ICHECK,NPOINTS(ICHECK))-
     &                  YY(ISEG,NPOINTS(ISEG)))
                DX3=ABS(XX(ICHECK,1)-XX(ISEG,1))
                DY3=ABS(YY(ICHECK,1)-YY(ISEG,1))
                DX4=ABS(XX(ICHECK,1)-XX(ISEG,NPOINTS(ISEG)))
                DY4=ABS(YY(ICHECK,1)-YY(ISEG,NPOINTS(ISEG)))
                IF(DX1.LT.TOL .AND. DY1.LT.TOL .AND. 
     &             DX2.LT.TOL .AND. DY2.LT.TOL) THEN
                   PRINT *,'DUPLICATE 1'
                   IUSED(ISEG)=1
                ELSEIF(DX3.LT.TOL .AND. DY3.LT.TOL .AND. 
     &             DX4.LT.TOL .AND. DY4.LT.TOL) THEN
                   PRINT *,'DUPLICATE 2'
                   IUSED(ISEG)=1
                ELSEIF(DX1.LT.TOL .AND. DY1.LT.TOL) THEN
                  IUSED(ISEG)=1
                  NNSEG=NNSEG+1
                  IF(JCHECK.EQ.0) THEN
                    WRITE(11,202) XX(ICHECK,1),YY(ICHECK,1),1
                    CALL MOVE(XX(ICHECK,1),YY(ICHECK,1))
                    DO NP=2,NPOINTS(ICHECK)
                      WRITE(11,202) XX(ICHECK,NP),YY(ICHECK,NP),1
                      CALL DRAW(XX(ICHECK,NP),YY(ICHECK,NP))
                    ENDDO
                  ENDIF
                  WRITE(11,202) XX(ISEG,2),YY(ISEG,2),11
                  CALL DRAW(XX(ISEG,2),YY(ISEG,2))
                  DO NP=3,NPOINTS(ISEG)
                    WRITE(11,202) XX(ISEG,NP),YY(ISEG,NP),11
                    CALL DRAW(XX(ISEG,NP),YY(ISEG,NP))
                  ENDDO
C                  WRITE(11,*) 'CASE 1'
                  ICHECK=ISEG
                  JCHECK=1
                  GOTO 210
                ELSEIF(DX2.LT.TOL .AND. DY2.LT.TOL) THEN
                  IUSED(ISEG)=1
                  NNSEG=NNSEG+1
                  IF(JCHECK.EQ.0) THEN
                    WRITE(11,202) XX(ICHECK,1),
     &                            YY(ICHECK,1),2
                    CALL MOVE(XX(ICHECK,1),
     &                        YY(ICHECK,1))
                    DO NP=2,NPOINTS(ICHECK)
                      WRITE(11,202) XX(ICHECK,NP),YY(ICHECK,NP),2
                      CALL DRAW(XX(ICHECK,NP),YY(ICHECK,NP))
                    ENDDO
                  ENDIF
                  WRITE(11,202) XX(ISEG,NPOINTS(ISEG)-1),
     &                          YY(ISEG,NPOINTS(ISEG)-1),22
                  CALL DRAW(XX(ISEG,NPOINTS(ISEG)-1),
     &                      YY(ISEG,NPOINTS(ISEG)-1))
                  DO NP=NPOINTS(ISEG)-2,1,-1
                    WRITE(11,202) XX(ISEG,NP),YY(ISEG,NP),22
                    CALL DRAW(XX(ISEG,NP),YY(ISEG,NP))
                  ENDDO

                  ICHECK=ISEG
                  JCHECK=1
C                  WRITE(11,*) 'CASE 2'
                  GOTO 210
                ELSEIF(DX3.LT.TOL .AND. DY3.LT.TOL) THEN
                  IUSED(ISEG)=1
                  NNSEG=NNSEG+1
                  IF(JCHECK.EQ.0) THEN
                    WRITE(11,202) XX(ICHECK,1),YY(ICHECK,1),3
                    CALL MOVE(XX(ICHECK,1),YY(ICHECK,1))
                    DO NP=2,NPOINTS(ICHECK)
                      WRITE(11,202) XX(ICHECK,NP),YY(ICHECK,NP),3
                      CALL DRAW(XX(ICHECK,NP),YY(ICHECK,NP))
                    ENDDO
                  ENDIF
                  WRITE(11,202) XX(ISEG,2),YY(ISEG,2),33
                  CALL DRAW(XX(ISEG,2),YY(ISEG,2))
                  DO NP=3,NPOINTS(ISEG)
                    WRITE(11,202) XX(ISEG,NP),YY(ISEG,NP),33
                    CALL DRAW(XX(ISEG,NP),YY(ISEG,NP))
                  ENDDO
                  ICHECK=ISEG
                  JCHECK=1
C                  WRITE(11,*) 'CASE 3'
                  GOTO 210
               ELSEIF(DX4.LT.TOL .AND. DY4.LT.TOL) THEN
                  NNSEG=NNSEG+1
                  IUSED(ISEG)=1
                  IF(JCHECK.EQ.0) THEN
                    WRITE(11,202) XX(ICHECK,NPOINTS(ICHECK)),
     &                            YY(ICHECK,NPOINTS(ICHECK)),4
                    CALL MOVE(XX(ICHECK,NPOINTS(ICHECK)),
     &                        YY(ICHECK,NPOINTS(ICHECK)))
                    DO NP=NPOINTS(ICHECK)-1,1,-1
                      WRITE(11,202) XX(ICHECK,NP),YY(ICHECK,NP),4
                      CALL DRAW(XX(ICHECK,NP),YY(ICHECK,NP))
                    ENDDO
                  ENDIF
                  WRITE(11,202) XX(ISEG,NPOINTS(ISEG)-1),
     &                          YY(ISEG,NPOINTS(ISEG)-1),44
                  CALL DRAW(XX(ISEG,NPOINTS(ISEG)-1),
     &                      YY(ISEG,NPOINTS(ISEG)-1))
                  DO NP=NPOINTS(ISEG)-2,1,-1
                    WRITE(11,202) XX(ISEG,NP),YY(ISEG,NP),44
                    CALL DRAW(XX(ISEG,NP),YY(ISEG,NP))
                  ENDDO
                  ICHECK=ISEG
                  JCHECK=1
C                  WRITE(11,*) 'CASE 4'
                  GOTO 210
                ENDIF
              ENDIF
            ENDDO
            PRINT *,'NOT FOUND, STARTING OVER'
            PRINT *,'NNSEG=',NNSEG
            IF(NNSEG.GT.0) THEN
              WRITE(11,202) -99999.,21.,IC
              WRITE(11,*) LINES(NCONT)
            ENDIF
          ENDIF
        ENDDO
C      WRITE(11,202) -99999.,21.,IC
      IC=IC+1
       IF(IC.GT.10) IC=1
202   FORMAT(10X,G13.6,2X,G13.6,I13)
C      WRITE(11,*) LINES(NCONT)
      IF(IFINISH.EQ.1) GOTO 99
      ENDDO
99    CONTINUE
      CALL GRSTOP1
      END
      REAL FUNCTION PARSE(LINE)
      CHARACTER*80 LINE,JUNK
      WRITE(JUNK,*) LINE
      READ(JUNK,1000) VAL
1000  FORMAT(1X,F20.0)
      PARSE=VAL
      END
