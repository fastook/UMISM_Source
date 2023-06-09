      IMPLICIT REAL*8(A-H,O-Z)  
      include "parameter.h"                                       
      PARAMETER(MXX=MAXNUM,HMIN=10.,NPLOT=101)                           
      CHARACTER HED*80                                                  
      DIMENSION XX(MXX),YY(MXX),HTICE(MXX,NPLOT)
      DIMENSION KX(MXX,4),CONST(MXX),ACON(MXX)
      DIMENSION PSURF(MXX),FRACT(MXX),FLOWA(MXX),SLDGB(MXX)             
      DIMENSION ADOT(MXX,NPLOT),DEPB(MXX,NPLOT),KODE(MXX),AJUNK(MXX)
      logical iflush
      common /flush/ iflush
      data iflush /.true./
C READ INPUT HEADER                   
      PRINT *,'INPUT START, STEP, HOW MANY?'
      READ(*,*) ISTART,ISKIP,NPLT
      ISKIP=ISKIP-1
      NPLT=MIN(NPLT,NPLOT)
      READ(30,1000,END=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,   
     &              INTER,DT                                            
1000  FORMAT (A80,/,7I6,F8.0)                                           
C READ INPUT GRID, THINGS THAT NEVER CHANGE                             
      READ(31) HED                                                      
      READ(31) (KODE(I),I=1,NUMNP)                                         
      READ(31) (XX(I),I=1,NUMNP)                                        
      READ(31) (YY(I),I=1,NUMNP)                                        
      READ(31) (PSURF(I),I=1,NUMNP)                                     
      READ(31) (AJUNK(I),I=1,NUMNP)                                        
      READ(31) (KX(I,1),KX(I,2),KX(I,3),KX(I,4),I=1,NUMEL)              
      READ(31) (IB1,IB2,BJUNK,I=1,NUMGBC)                               
C READ INPUT DIFFER, THINGS THAT DIFFER FROM ONE GRID TO THE NEXT       
      READ(32) HED                                                      
      READ(32) (AJUNK(I),I=1,NUMNP)                                        
      READ(32) (FRACT(I),I=1,NUMNP)                                     
      READ(32) (FLOWA(I),I=1,NUMNP)                                     
      READ(32) (SLDGB(I),I=1,NUMNP)                                     
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME                         
      XMIN=1E30
      XMAX=-XMIN
      YMIN=1E30
      YMAX=-YMAX
      II=0
      DO K=1,ISTART
        READ(33,END=99) HED                                              
        II=II+1
        WRITE(*,*) 'BYPASS ',II,' ',HED(1:40)                                              
        READ(33) (RJUNK,I=1,NUMNP)                                     
        READ(33) (RJUNK,I=1,NUMNP)                                      
        READ(33) (RJUNK,I=1,NUMNP)                                      
        READ(33) (RJUNK,I=1,NUMEL)                                     
        READ(33) (RJUNK,I=1,NUMEL)                                     
      ENDDO
      DO J=1,NPLT
        READ(33,END=99) HED                                              
        II=II+1
        WRITE(*,*) II,' ',HED(1:40)
        NUM=J
        READ(33) (HTICE(I,J),I=1,NUMNP)                                     
        READ(33) (ADOT(I,J),I=1,NUMNP)                                      
        READ(33) (DEPB(I,J),I=1,NUMNP)                                      
        READ(33) (CONST(I),I=1,NUMEL)                                     
        READ(33) (ACON(I),I=1,NUMEL)                                     
        DO I=1,NUMNP
          XMAX=MAX(XMAX,HTICE(I,J))
          XMIN=MIN(XMIN,HTICE(I,J))
          YMAX=MAX(YMAX,ADOT(I,J))
          YMIN=MIN(YMIN,ADOT(I,J))
        ENDDO
        DO K=1,ISKIP
          READ(33,END=99) HED                                              
          II=II+1
          WRITE(*,*) ' SKIP ',II,' ',HED(1:40)                                              
          READ(33) (RJUNK,I=1,NUMNP)                                     
          READ(33) (RJUNK,I=1,NUMNP)                                      
          READ(33) (RJUNK,I=1,NUMNP)                                      
          READ(33) (RJUNK,I=1,NUMEL)                                     
          READ(33) (RJUNK,I=1,NUMEL)                                     
        ENDDO
      ENDDO
99    CONTINUE
      PRINT *,'DONE READING DATA',XMIN,XMAX,YMIN,YMAX,NUM
      YMIN=-1.
      YMAX=2.
      XMIN=0.
      XMAX=4000.
      CALL GRSTRT(800,1)
      XBORD=(XMAX-XMIN)/5.D0
      YBORD=(YMAX-YMIN)/5.D0
      XMIN1=XMIN-XBORD
      XMAX1=XMAX+XBORD
      YMIN1=YMIN-YBORD
      YMAX1=YMAX+YBORD
C     CALL WINDOW(REAL(XMIN1),REAL(XMAX1),REAL(YMIN1),REAL(YMAX1))
      CALL WINDOW(REAL(YMIN1),REAL(YMAX1),REAL(XMIN1),REAL(XMAX1))
      CALL OPNSEG(1)
      CALL LINCLR(1)
C     CALL MOVE(REAL(XMIN),0.)
C     CALL DRAW(REAL(XMAX),0.)
C     CALL MOVE(0.,REAL(YMIN))
C     CALL DRAW(0.,REAL(YMAX))
C     CALL XAXIS('X',5.,REAL(XMIN),REAL(XMAX),REAL(YMAX),REAL(YMIN))
C     CALL XAXIS('X',5.,REAL(XMIN),REAL(XMAX),REAL(YMIN),REAL(YMAX))
C     CALL YAXIS('Y',5.,REAL(XMIN),REAL(XMAX),REAL(YMIN),REAL(YMAX))
C     CALL YAXIS('Y',5.,REAL(XMAX),REAL(XMIN),REAL(YMIN),REAL(YMAX))
      CALL MOVE(REAL(YMIN),0.)
      CALL DRAW(REAL(YMAX),0.)
      CALL MOVE(0.,REAL(XMIN))
      CALL DRAW(0.,REAL(XMAX))
      CALL XAXIS('X',5.,REAL(YMIN),REAL(YMAX),REAL(XMAX),REAL(XMIN))
      CALL XAXIS('X',5.,REAL(YMIN),REAL(YMAX),REAL(XMIN),REAL(XMAX))
      CALL YAXIS('Y',5.,REAL(YMIN),REAL(YMAX),REAL(XMIN),REAL(XMAX))
      CALL YAXIS('Y',5.,REAL(YMAX),REAL(YMIN),REAL(XMIN),REAL(XMAX))
      CALL CLOSEG
      CALL SETVIS(1,.TRUE.)
      DO ISEG=2,NUM
        J=ISEG
        IF(ISEG.GT.2) CALL SETVIS(ISEG-1,.FALSE.)
        CALL OPNSEG(ISEG)
        DO I=1,NUMNP
C         CALL MOVE(REAL(HTICE(I,J)),REAL(ADOT(I,J)))
          CALL MOVE(REAL(ADOT(I,J)),REAL(HTICE(I,J)))
          IF(HTICE(I,J).GT.DEPB(I,J)+1.) THEN
            CALL LINCLR(3)
          ELSE
            CALL LINCLR(2)
          ENDIF
C         CALL DRAW(REAL(HTICE(I,J)),REAL(ADOT(I,J)))
          CALL DRAW(REAL(ADOT(I,J)),REAL(HTICE(I,J)))
          CALL POINT(REAL(ADOT(I,J)),REAL(HTICE(I,J)))
        ENDDO
        CALL LSTSQ(NUM,HTICE(1,J),ADOT(1,J),AAA,BBB)
        YP1=AAA*XMIN+BBB
        YP2=AAA*XMAX+BBB
        XP1=(YMIN-BBB)/AAA
        XP2=(YMAX-BBB)/AAA
        CALL LINCLR(1)
        CALL MOVE(REAL(YP1),REAL(XMIN))
C       CALL DRAW(REAL(YP2),REAL(XMAX))
C       CALL MOVE(REAL(YMIN),REAL(XP1))
        CALL DRAW(REAL(YMAX),REAL(XP2))
        CALL CLOSEG
        CALL SETVIS(ISEG,.TRUE.)
      ENDDO
999   CONTINUE
      PRINT *,'0 TO PAUSE, 1 TO SCAN THRU'
      READ(*,*) ISKIP
      CALL SETVIS(NUM,.FALSE.)
      CALL SETVIS(1,.TRUE.)
      IF(ISKIP.EQ.0) READ(*,*) IJUNK
      CALL SETVIS(2,.TRUE.)
      IF(ISKIP.EQ.0) READ(*,*) IJUNK
      DO I=3,NUM
        CALL SETVIS(I,.TRUE.)
      IF(ISKIP.EQ.0) READ(*,*) IJUNK
        CALL WAIT(IAGAIN)
        CALL SETVIS(I-1,.FALSE.)
      IF(ISKIP.EQ.0) READ(*,*) IJUNK
      ENDDO
      if(iflush) call gflush
      WRITE(*,*) 'INPUT 1 TO SEE AGAIN, 0 TO END'
      READ(*,*) IAGAIN
      IF(IAGAIN.GT.0) GOTO 999
      CALL GRSTOP1
      END
      SUBROUTINE XAXIS(LABEL,XLEN,XMIN,XMAX,YMIN,YMAX)
      CHARACTER*20 LABEL
      CALL TXFONT(0)
      XTIC=(XMAX-XMIN)/100.
      YTIC=(YMAX-YMIN)/100.
      CALL MOVE(XMIN,YMIN)
      DO 10 X=0.,XLEN,1.
        XP=XMIN+X*(XMAX-XMIN)/XLEN
        CALL DRAW(XP,YMIN)
        CALL DRAW(XP,YMIN-3.*YTIC)
        CALL DRAW(XP,YMIN+3.*YTIC)
        CALL MOVE(XP-7.*XTIC,YMIN-7.*YTIC)
        CALL GNUMBR(XP,2,8)
        CALL MOVE(XP,YMIN)
10    CONTINUE
      CALL MOVE(XMIN,YMIN)
      DO 11 X=0.,XLEN,.1
        XP=XMIN+X*(XMAX-XMIN)/XLEN
        CALL MOVE(XP,YMIN)
        CALL DRAW(XP,YMIN-YTIC)
11    CONTINUE
      CALL MOVE(XMIN+30.*XTIC,YMIN-15.*YTIC)
      CALL TEXT(20,LABEL)
      END
      SUBROUTINE YAXIS(LABEL,YLEN,XMIN,XMAX,YMIN,YMAX)
      CHARACTER*20 LABEL
      CALL TXFONT(0)
      CALL TXANGL(90.)
      XTIC=(XMAX-XMIN)/100.
      YTIC=(YMAX-YMIN)/100.
      CALL MOVE(XMIN,YMIN)
      DO 10 Y=0.,YLEN,1.
        YP=YMIN+Y*(YMAX-YMIN)/YLEN
        CALL DRAW(XMIN,YP)
        CALL DRAW(XMIN-3.*XTIC,YP)
        CALL DRAW(XMIN+3.*XTIC,YP)
        CALL MOVE(XMIN-7.*XTIC,YP-8.*YTIC)
        CALL  GNUMBR(YP,2,8)
        CALL MOVE(XMIN,YP)
10    CONTINUE
      CALL MOVE(XMIN,YMIN)
      DO 11 Y=0.,YLEN,.1
        YP=YMIN+Y*(YMAX-YMIN)/YLEN
        CALL MOVE(XMIN,YP)
        CALL DRAW(XMIN-XTIC,YP)
11    CONTINUE
      CALL MOVE(XMIN-15.*XTIC,YMIN+23.*YTIC)
C     CALL TXANGL(90.)
      CALL TEXT(20,LABEL)
      CALL TXANGL(0.)
      END
      SUBROUTINE GNUMBR(X,IDIGIT,IWIDE)
      CHARACTER*20 JUNK
      IF(ABS(X).GE..1 .AND. ABS(X).LE.999999.) THEN
        WRITE(JUNK,100) X
100     FORMAT(G13.6)
      ELSEIF(X.EQ.0.) THEN
        WRITE(JUNK,101) X
101     FORMAT(2X,F7.5)
      ELSE
        WRITE(JUNK,102) X
102     FORMAT(1PE9.2)
      ENDIF
C     WRITE(JUNK,103) X
103   FORMAT(F6.0)
      CALL TEXT(20,JUNK)
      END
      SUBROUTINE ZOOM(A,B,C,D)
      PRINT *,'NOT WORKING'
      END
      SUBROUTINE TXANGL(A)
      END
      SUBROUTINE TXFONT(I)
      END
      SUBROUTINE LSTSQ(N,X,Y,AAA,BBB)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION X(N),Y(N)
      EXTERNAL AVG
      XBAR=AVG(N,X)
      YBAR=AVG(N,Y)
      SUMTOP=0.0
      SUMBOT=0.0
      DO I=1,N
        SUMTOP=SUMTOP+(X(I)-XBAR)*(Y(I)-YBAR)
        SUMBOT=SUMBOT+(X(I)-XBAR)**2
      ENDDO
      AAA=SUMTOP/SUMBOT
      BBB=YBAR-AAA*XBAR
      END
      FUNCTION AVG(N,X)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION X(N)
      SUM=0.0
      DO I=1,N
        SUM=SUM+X(I)
      ENDDO
      AVG=SUM/N
      END
      subroutine wait(n)
      sum=0
      do i=1,n
        sum=sum+sin(real(i))
      enddo
      end
