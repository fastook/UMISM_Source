C **************************************************************C       
C                                                               C       
C   PROGRAM:  VELO4a                                            C       
C                                                               C       
C   DATE:     1995                                              C       
C   PROGRAMMER:  FASTOOK                                        C       
C                                                               C       
C   FUNCTION:                                                   C       
C            PLOTS VELOCITY VECTORS COLOR-CODED FOR MAGNITUDE   C       
C            FROM THE NEW SPLIT FORMAT OUTPUT SET **            C       
C            AFTER ADJUSTING SCALE (PAIR, FIRST LENGHT OF       C       
C            ARROW, SECOND MAX VELOCITY) OUTPUTS STANDARD       C       
C            PLOTTER DATA SET TO OUTV** DATA B                  C       
C                                                               C       
C **************************************************************C       
      IMPLICIT REAL*8(A-H,O-Z)                                          
      PARAMETER(MXX=79999,HMIN=10.,NPLOT=10000)                           
      CHARACTER HED*80                                                  
      DIMENSION ICMAP(16),INARAY(2)
      DIMENSION XARO(9),YARO(9),NBOUND(MXX)                 
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2),LM(4)             
      DIMENSION PSI(4),DPSI(4,2)                                        
      DIMENSION XY(2,4),XI(9),ETA(9),W(9)                               
      DIMENSION XX(MXX),YY(MXX),HTICE(MXX),THICK(MXX),NTYPE(MXX)        
      DIMENSION XLINE(MXX),YLINE(MXX)                                   
      DIMENSION KX(MXX,4),IK(5),CONST(MXX)                              
      DIMENSION PSURF(MXX),FRACT(MXX),FLOWA(MXX),SLDGB(MXX)             
      DIMENSION ADOT(MXX),DEPB(MXX),KODE(MXX),AJUNK(MXX)
      DATA ICMAP /0,4,12,6,13,2,8,7,9,3,10,5,11,1,14,15/                
C     DATA ICMAP /0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/                
      DIMENSION XOUT(NPLOT),YOUT(NPLOT),IOUT(NPLOT)    
      LOGICAL LPLOT
      logical iflush
      common /flush/ iflush
      data iflush /.false./
      LPLOT=.true.                 
      CALL SETRIG
      print *,'input 1 for normalized vectors, 0 for absolute lengths'
      if(iflush) call gflush
      read(*,*) INORM
      XLEN=5.0                                                          
      YLEN=5.0                                                          
      RNINE=-99999.                                                     
      RMINUS=-1.                                                        
      RZERO=0.                                                          
      RTWO=2.                                                           
      RAROW=20.                                                         
      ICOLOR=1                                                          
      ICSAV=ICOLOR                                                      
      II=1                                                              
      DO I=1,NPLOT                                                      
        READ(12,*,END=21) XIN,YIN,IIN                                   
        IF(XIN.EQ.-99999.) THEN                                         
          READ(12,22) JUNK                                              
22        FORMAT(A80)                                                       
        ELSE                                                            
          XOUT(II)=XIN                                                  
          YOUT(II)=YIN                                                  
          IOUT(II)=IIN                                                  
          NOUT=II                                                       
          II=II+1                                                       
        ENDIF                                                           
      ENDDO                                                             
21    CONTINUE                                                          
      WRITE(*,*) 'INPUT SCALE FACTOR FOR ARROWS, UMAX, AND THRESHOLD'   
      if(iflush) call gflush
      READ(*,*) ASCAL,UMAX,THRESH                                       
      ASCAL1=ASCAL
      UDELTA=UMAX/13.                                                   
      IF(LPLOT) CALL GRSTRT(800,1) 
C                                               
C ASCII FORMATTED READ...
      READ(1,100,END=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &                    INTER,DT
100   FORMAT(A80,/,7I6,E15.6)
      DO NUM=1,NUMNP
        READ(1,200) N,KODE(NUM),XX(NUM),YY(NUM),HTICE(NUM),
     &              ADOT(NUM),FRACT(NUM),DENS,DEPB(NUM),
     &              FLOWA(NUM),SLDGB(NUM)
      ENDDO
200   FORMAT(I6,I4,2E12.5,F10.2,F7.2,F9.2,F10.3,F10.1,F10.2,F10.3)
      DO I=1,NUMEL
        READ(1,300) NUM, (KX(NUM,II),II=1,4),CONST(NUM)
        NTYPE(NUM)=1
        IF(KX(NUM,4).EQ.0) NTYPE(NUM)=2
      ENDDO
300   FORMAT(5I6,1PE17.10)
      DO N=1,NUMGBC
        READ(1,310) I,J,RJUNK
      ENDDO
310   FORMAT(2I6,E13.6)
C
      NLINE=0                                                           
      READ(11,*,END=99) NLINE                                           
      READ(11,*) (NBOUND(I),I=1,NLINE)                                  
      REWIND 11                                                         
99    CONTINUE                                                          
      XMAX=-1.E30                                                       
      XMIN=1.E30                                                        
      YMAX=-1.E30                                                       
      YMIN=1.E30                                                        
      IF(NLINE.GT.0) THEN                                               
        DO 15 I=1,NLINE                                                 
        XLINE(I)=XX(NBOUND(I))/1000.                                    
        YLINE(I)=YY(NBOUND(I))/1000.                                    
        XMAX=MAX(XMAX,XLINE(I))
        XMIN=MIN(XMIN,XLINE(I))
        YMAX=MAX(YMAX,YLINE(I))
        YMIN=MIN(YMIN,YLINE(I))
15      CONTINUE                                                        
105     FORMAT(10X,G13.6,2X,G13.6,I13)                                  
      ELSE                                                              
        DO 151 I=1,NUMNP                                                
          XMAX=MAX(XMAX,XX(I)*.001)
          XMIN=MIN(XMIN,XX(I)*.001)
          YMAX=MAX(YMAX,YY(I)*.001)
          YMIN=MIN(YMIN,YY(I)*.001)
151     CONTINUE                                                        
      ENDIF                                                             
      XBORD=(XMAX-XMIN)/10.                                             
      XMIN=XMIN-XBORD                                                   
      XMAX=XMAX+XBORD                                                   
      YBORD=(YMAX-YMIN)/10.                                             
      YMIN=YMIN-YBORD                                                   
      YMAX=YMAX+YBORD                                                   
c     WRITE(7,105) XMIN,YMIN
c     WRITE(7,105) XMAX,YMIN
c     WRITE(7,105) XMAX,YMAX
c     WRITE(7,105) XMIN,YMIN
c     WRITE(7,105) XMIN,YMIN
c     WRITE(7,105) RNINE,RTWO                                         
c     WRITE(7,111)                                                    
111     FORMAT('BOUNDARY')                                              
      IPASS=0
      ISEG=1
      REWIND 1
1     CONTINUE                                                          
C ASCII FORMATTED READ...
      READ(1,100,END=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &                    INTER,DT
      DO NUM=1,NUMNP
        READ(1,200) N,KODE(NUM),XX(NUM),YY(NUM),HTICE(NUM),
     &              ADOT(NUM),FRACT(NUM),DENS,DEPB(NUM),
     &              FLOWA(NUM),SLDGB(NUM)
      ENDDO
      DO I=1,NUMEL
        READ(1,300) NUM, (KX(NUM,II),II=1,4),CONST(NUM)
        NTYPE(NUM)=1
        IF(KX(NUM,4).EQ.0) NTYPE(NUM)=2
      ENDDO
      DO N=1,NUMGBC
        READ(1,310) I,J,RJUNK
      ENDDO
C

      DO NUM=1,NUMNP                                                    
        THICK(NUM)=HTICE(NUM)-DEPB(NUM)                                     
      ENDDO                                                             
23    CONTINUE                                                          
      IF(LPLOT) CALL OPNSEG(ISEG)
      IF(IPASS.EQ.0) THEN
        IF(LPLOT) THEN
          CALL WINDOW(REAL(XMIN),REAL(XMAX),REAL(YMIN),REAL(YMAX))
        ENDIF
        IPASS=1
      ENDIF
      IF(LPLOT) THEN
        CALL LINCLR(1)                                                    
        CALL MOVE(REAL(XMIN),REAL(YMIN))                                  
        CALL DRAW(REAL(XMIN),REAL(YMAX))                                  
        CALL DRAW(REAL(XMAX),REAL(YMAX))                                  
        CALL DRAW(REAL(XMAX),REAL(YMIN))                                  
        CALL DRAW(REAL(XMIN),REAL(YMIN))                                  
        CALL DASHPT(3)                                                    
        DO I=1,NOUT                                                       
          IF(IOUT(I).EQ.1) THEN                                           
            CALL MOVE(REAL(XOUT(I)),REAL(YOUT(I)))                        
          ELSE                                                            
            CALL DRAW(REAL(XOUT(I)),REAL(YOUT(I)))                        
          ENDIF                                                           
        ENDDO                                                             
        CALL LINCLR(1)                                                    
        IF(NLINE.GT.0) THEN                                               
          CALL MOVE(REAL(XLINE(1)),REAL(YLINE(1)))
          DO 73 I=2,NLINE                                                 
            CALL DRAW(REAL(XLINE(I)),REAL(YLINE(I)))
73        CONTINUE                                                        
        ENDIF 
      ENDIF                                                            
c      REWIND 9                                                          
c      WRITE(9,*) HED
c      WRITE(9,*) 'ELEMENT,X,Y,UX,UY,UMAG'                               
      VMAX=-1.E30                                                       
      DO 600 J=1,NUMEL                                                  
        IF(NTYPE(J).EQ.1) THEN                                            
          NNODE=4                                                         
          CENTX=0.0D00                                                    
          CENTY=0.0D00                                                    
        ELSE                                                              
          NNODE=3                                                         
          CENTX=1.D0/3.D0                                                 
          CENTY=1.D0/3.D0                                                 
        ENDIF                                                             
        HH=0.0D00                                                         
        SUMX=0.0D00                                                       
        SUMY=0.0D00                                                       
        DO I=1,NNODE                                                  
          LM(I)=KX(J,I)                                                     
        ENDDO
        I=LM(1)                                                           
        JJ=LM(2)                                                          
        K=LM(3)                                                           
        L=LM(4)                                                           
        XY(1,1)=XX(I)                                                     
        XY(1,2)=XX(JJ)                                                    
        XY(1,3)=XX(K)                                                     
        IF(NTYPE(J).EQ.1) XY(1,4)=XX(L)                                   
        XY(2,1)=YY(I)                                                     
        XY(2,2)=YY(JJ)                                                    
        XY(2,3)=YY(K)                                                     
        IF(NTYPE(J).EQ.1) XY(2,4)=YY(L)                                   
        IF(NTYPE(J).EQ.1) THEN                                            
          XCENT=(XX(I)+XX(JJ)+XX(K)+XX(L))/4000.                          
          YCENT=(YY(I)+YY(JJ)+YY(K)+YY(L))/4000.                          
        ELSE                                                              
          XCENT=(XX(I)+XX(JJ)+XX(K))/3000.                                
          YCENT=(YY(I)+YY(JJ)+YY(K))/3000.                                
        ENDIF                                                             
C                                                                       
        CALL SHAPE(NTYPE(J),CENTX,CENTY,PSI,DPSI)                         
      
C CALCULATE DXDS...EQUATION (5.3.6)                                     
        DO I=1,2                                                      
          DO L=1,2                                                      
            DXDS(I,L)=0.0                                                     
            DO K=1,NNODE                                                  
              DXDS(I,L)=DXDS(I,L)+DPSI(K,L)*XY(I,K)
            ENDDO
          ENDDO
        ENDDO
   
C CALCULATE DSDX...EQUATION (5.2.7)                                     
        DETJ=(DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1))                    
        IF (DETJ.LE.0.0) GOTO 750                                         
c       IF (DETJ.GE.0.0) GOTO 750                                         
        DSDX(1,1)=DXDS(2,2)/DETJ                                          
        DSDX(2,2)=DXDS(1,1)/DETJ                                          
        DSDX(1,2)=-DXDS(1,2)/DETJ                                         
        DSDX(2,1)=-DXDS(2,1)/DETJ                                         
C CALCULATE D(PSI)/DX...EQUATION (5.3.5)                                
        DO I=1,NNODE                                                  
          DPSIX(I)=DPSI(I,1)*DSDX(1,1)+DPSI(I,2)*DSDX(2,1)                  
          DPSIY(I)=DPSI(I,1)*DSDX(1,2)+DPSI(I,2)*DSDX(2,2)                  
        ENDDO                                                          
        DO I=1,NNODE                                                  
          SUMX=SUMX + HTICE(LM(I))*DPSIX(I)                                 
          SUMY=SUMY + HTICE(LM(I))*DPSIY(I)                                 
          HH=HH + THICK(LM(I))*PSI(I)                                       
        ENDDO                                                          
C                                                                       
        DELH=SUMX**2 + SUMY**2                                            
        DELH=SQRT(DELH)                                                   
        IF(HH.GT.HMIN) THEN                                               
          UX=-CONST(J)*SUMX/HH                                            
          UY=-CONST(J)*SUMY/HH                                            
        ELSE                                                              
          UX=0.                                                           
          UY=0.                                                           
        ENDIF                                                             
        UMAG=SQRT(UX**2+UY**2)                                        
        VMAX=MAX(VMAX,UMAG)                                               
        IF(UMAG.GT.THRESH) GOTO 600                                       
        IF(UMAG.EQ.0.) THEN                                               
          ICOLOR=1                                                        
        ELSE                                                              
          ICOLOR=1+NINT((UMAG)/UDELTA)                                    
        ENDIF                                                             
        IF(ICOLOR.LT.2) ICOLOR=2                                          
        IF(ICOLOR.GT.14) ICOLOR=14                                        
C       ICOLOR=ICMAP(ICOLOR)                                              
c        WRITE(9,1001) J,XCENT,YCENT,UX,UY,UMAG,CONST(J)
1001  FORMAT(I6,6G13.6)                                                 
        XARO(1)=XCENT                                                     
        YARO(1)=YCENT                                                     
        IF(UMAG.NE.0) THEN
          IF(INORM.EQ.1) THEN
            XARO(2)=XCENT+ASCAL*UX/UMAG
            YARO(2)=YCENT+ASCAL*UY/UMAG
          ELSE
            XARO(2)=XCENT+ASCAL*UX
            YARO(2)=YCENT+ASCAL*UY
          ENDIF
        ELSE
          XARO(2)=XCENT                                                     
          YARO(2)=YCENT                                                     
        ENDIF
        WRITE(10,'(A1)') '>'
        CALL RECPOL(XARO(1),YARO(1),RLAT,RLONG)
        WRITE(10,*) REAL(RLONG),REAL(RLAT)
        CALL RECPOL(XARO(2),YARO(2),RLAT,RLONG)
        WRITE(10,*) REAL(RLONG),REAL(RLAT)
        WRITE(13,'(A1)') '>'
        WRITE(13,*) REAL(XARO(1)),REAL(YARO(1))
        WRITE(13,*) REAL(XARO(2)),REAL(YARO(2))
202   FORMAT(10X,G13.6,2X,G13.6,I13)                                    
c       WRITE(7,202) XARO(1),YARO(1)                                      
c       WRITE(7,202) XARO(2),YARO(2)                                      
c       WRITE(7,202) RNINE,RAROW,ICMAP(ICOLOR)-1                          
c       WRITE(7,*) ' ELEMENT',J                                           
c       IF(ICOLOR.NE.ICSAV) THEN                                          
          IF(LPLOT) CALL LINCLR(ICMAP(ICOLOR))                                      
          ICSAV=ICOLOR                                                    
c       ENDIF                                                             
        IF(LPLOT) CALL MOVE(REAL(XARO(1)),REAL(YARO(1)))
        IF(LPLOT) CALL DRAW(REAL(XARO(2)),REAL(YARO(2)))
600   CONTINUE                                                          
      XARO(1)=XMIN+(XMAX-XMIN)/XLEN                                     
      XARO(2)=XARO(1)+(XMAX-XMIN)/XLEN                                  
      YARO(1)=YMAX+(YMAX-YMIN)/YLEN                                     
      YARO(2)=YARO(1)                                                   
      IF(LPLOT) CALL MOVE(REAL(XARO(1)),REAL(YARO(1)))
      IF(LPLOT) CALL DRAW(REAL(XARO(2)),REAL(YARO(2)))
      IF(LPLOT) CALL CLOSEG
      IF(ISEG.NE.1 .AND.LPLOT) CALL SETVIS(ISEG-1,.FALSE.)
      ISEG=ISEG+1
      ZSCAL=(XMAX-XMIN)/ASCAL/XLEN                                      
C     WRITE(7,202) XARO(1),YARO(1)                                      
C     WRITE(7,202) XARO(2),YARO(2)                                      
c     WRITE(7,202) XMIN,YMIN                                            
c     WRITE(7,202) XMIN,YMIN                                            
c     WRITE(7,202) RNINE,RAROW+1.                                       
c     WRITE(7,2323) ZSCAL                                               
2323  FORMAT(G10.3)                                                     
c     CALL GRSTOP1                                                      
      WRITE(*,*) 'UMAX=',VMAX                                           
      IF(ASCAL1.GT.0.) then
        IF(LPLOT) call gflush()
        READ(*,*,END=24) ASCAL1,UMAX,THRESH 
      endif
      UDELTA=UMAX/13.                                                   
      IF(ASCAL1.LE.0.) GOTO 24                                          
      ASCAL=ASCAL1                                                      
      REWIND 7                                                          
      REWIND 10                                                          
c     WRITE(7,105) XMIN,YMIN
c     WRITE(7,105) XMAX,YMIN
c     WRITE(7,105) XMAX,YMAX
c     WRITE(7,105) XMIN,YMIN
c     WRITE(7,105) XMIN,YMIN
c     WRITE(7,105) RNINE,RTWO                                         
c     WRITE(7,111)                                                    
c     CALL GRSTRT(800,1)                                                
      GOTO 23                                                           
24    CONTINUE                                                          
      GOTO 1                                                            
999   CONTINUE                                                          
      NSEG=ISEG-1
900   CONTINUE
c     OPEN(*)
      IF(LPLOT) CALL SETVIS(NSEG,.FALSE.)
      CALL GETUIN(-1,'PAUSE 1, NO PAUSE 0',1,INARAY,IGOT)
      IF(IGOT.NE.1) THEN
        IPAUSE=0
      ELSE
        IPAUSE=INARAY(1)
      ENDIF
      CALL SETVIS(1,.TRUE.)
      DO 910 I=2,NSEG
        IF(LPLOT) CALL SETVIS(I,.TRUE.)
        IF(IPAUSE.EQ.1) CALL GETUIN(-1,'CONTINUE?',1,INARAY,IGOT)
        IF(LPLOT) CALL SETVIS(I-1,.FALSE.)
C       CALL NEWPAG
        IF(IPAUSE.EQ.1) CALL GETUIN(-1,'CONTINUE?',1,INARAY,IGOT)
910   CONTINUE
      CALL GETUIN(-1,'AGAIN? 1, FINISH? 0',1,INARAY,IGOT)
      IF(LPLOT) CALL SETVIS(NSEG,.FALSE.)
      IF(IGOT.EQ.1) THEN
        IF(INARAY(1).EQ.1) GOTO 900
      ENDIF
911   continue
      IF(LPLOT) call setvis(-1,.true.)
      print *,'input 1 to do again'
      IF(LPLOT .and. iflush) call gflush
      read(*,*) igot
      if(igot.eq.1) then
        IF(LPLOT) call newpag
        goto 911
      endif
      IF(LPLOT) CALL GRSTOP1                                                      
      STOP                                                              
 750  CONTINUE                                                          
      WRITE (*,5544) DETJ,X                                             
 5544 FORMAT (' BAD JACOBIAN at 750',E14.6,E14.6)                       
      STOP                                                              
      END                                                               
C******************************************************
      SUBROUTINE SHAPE(NTYPE,XI,ET,PSI,DPSI)                            
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION PSI(4),DPSI(4,2)                                        
      IF(NTYPE.EQ.2) GOTO 100                                           
      PSI(1)=.25*(1.-XI)*(1.-ET)                                        
      PSI(2)=.25*(1.+XI)*(1.-ET)                                        
      PSI(3)=.25*(1.+XI)*(1.+ET)                                        
      PSI(4)=.25*(1.-XI)*(1.+ET)                                        
      DPSI(1,2)=-.25*(1.-XI)                                            
      DPSI(2,2)=-.25*(1.+XI)                                            
      DPSI(3,2)=.25*(1.+XI)                                             
      DPSI(4,2)=.25*(1.-XI)                                             
      DPSI(1,1)=-.25*(1.-ET)                                            
      DPSI(2,1)=.25*(1.-ET)                                             
      DPSI(3,1)=.25*(1.+ET)                                             
      DPSI(4,1)=-.25*(1.+ET)                                            
      RETURN                                                            
100   CONTINUE                                                          
      PSI(1)=1.-XI-ET                                                   
      PSI(2)=XI                                                         
      PSI(3)=ET                                                         
      DPSI(1,2)=-1.                                                     
      DPSI(2,2)=0.                                                      
      DPSI(3,2)=1.                                                      
      DPSI(1,1)=-1.                                                     
      DPSI(2,1)=1.                                                      
      DPSI(3,1)=0.                                                      
      RETURN                                                            
      END                                                               
C******************************************************
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
