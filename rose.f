C **************************************************************C       
C                                                               C       
C   PROGRAM:  VELO5                                             C       
C                                                               C       
C   DATE:  12 04 91                                             C       
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
      LOGICAL LDRAW                                  
      PARAMETER(NMAX=29999,HMIN=10.,NPLOT=100000)                           
      CHARACTER HED*80     
      REAL*4 XA,YA                                             
      DIMENSION ICMAP(16),INARAY(2)
      DIMENSION XARO(9),YARO(9),NBOUND(NMAX)                 
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2),LM(4)             
      DIMENSION PSI(4),DPSI(4,2)                                        
      DIMENSION XY(2,4),XA(1),YA(1),IDAT(1)
      DIMENSION XX(NMAX),YY(NMAX),HTICE(NMAX),THICK(NMAX),NTYPE(NMAX)        
      DIMENSION XLINE(NMAX),YLINE(NMAX)                                   
      DIMENSION KX(NMAX,4),CONST(NMAX)
      DIMENSION PSURF(NMAX),FRACT(NMAX),FLOWA(NMAX),SLDGB(NMAX)             
      DIMENSION ADOT(NMAX),DEPB(NMAX),KODE(NMAX),AJUNK(NMAX)
      DATA ICMAP /0,4,12,6,13,2,8,7,9,3,10,5,11,1,14,15/                
C     DATA ICMAP /0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/                
      DIMENSION XOUT(NPLOT),YOUT(NPLOT),IOUT(NPLOT)                     
      DATA RHOI /0.917D0/, RHOW /1.092D0/, THRESH /1.D0/
      PG = 0.089866*0.3816
      ipass=0
      LDRAW=.false.
      PI=4.D0*ATAN(1.D0)
      IVELO=0
c     PRINT *,'TOTAL VELOCITY (0) OR JUST SLIDING (1)'
c     READ(*,*) IVELO
      CALL SETRIG
C READ INPUT HEADER                                                     
      READ(30,1000,END=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,   
     &              INTER,DT                                            
1000  FORMAT (A80,/,7I6,F8.0)                                           
C READ INPUT GRID, THINGS THAT NEVER CHANGE                             
      READ(31) HED                                                      
      READ(31) (KODE(I),I=1,NUMNP)                                         
      READ(31) (XX(I),I=1,NUMNP)                                        
      READ(31) (YY(I),I=1,NUMNP)                                        
      do i=1,numnp
        call rotate(xx(i),yy(i))
      enddo
      READ(31) (PSURF(I),I=1,NUMNP)                                     
      READ(31) (AJUNK(I),I=1,NUMNP)                                        
      READ(31) (KX(I,1),KX(I,2),KX(I,3),KX(I,4),I=1,NUMEL)              
      DO NUM=1,NUMEL                                                      
        NTYPE(NUM)=1                                                    
        IF(KX(NUM,4).EQ.0) NTYPE(NUM)=2                                 
      ENDDO                                                             
      READ(31) (IB1,IB2,BJUNK,I=1,NUMGBC)                      
C READ INPUT DIFFER, THINGS THAT DIFFER FROM ONE GRID TO THE NEXT       
      READ(32) HED                                                      
      READ(32) (AJUNK(I),I=1,NUMNP)                                        
      READ(32) (FRACT(I),I=1,NUMNP)                                     
      READ(32) (FLOWA(I),I=1,NUMNP)                                     
      READ(32) (SLDGB(I),I=1,NUMNP)                                     
      READ(33,END=999) HED                                              
      READ(33) (HTICE(I),I=1,NUMNP)                                     
      READ(33) (ADOT(I),I=1,NUMNP)                                      
      READ(33) (DEPB(I),I=1,NUMNP)  
      rewind 33                                    
      NLINE=0                                                           
      READ(11,*,END=99) NLINE                                           
      READ(11,*) (NBOUND(I),I=1,NLINE)                                  
      REWIND 11                                                         
99    CONTINUE                                                          
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
      do i=1,NOUT
        call rotate(XOUT(i),YOUT(i))
      enddo
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
      XWID=XMAX-XMIN
      YWID=YMAX-YMIN
      IF(XWID.GT.YWID) THEN
        YWID=XWID
        YMAX=YMIN+YWID
      ELSE
        XWID=YWID
        XMAX=XMIN+XWID
      ENDIF
C 105     FORMAT(10X,G13.6,2X,G13.6,I13)                                  
C     WRITE(7,105) XMIN,YMIN
C     WRITE(7,105) XMAX,YMIN
C     WRITE(7,105) XMAX,YMAX
C     WRITE(7,105) XMIN,YMIN
C     WRITE(7,105) XMIN,YMIN
C     WRITE(7,105) RNINE,RTWO                                         
C     WRITE(7,111)                                                    
c 111     FORMAT('BOUNDARY')                                              
      ivect=1
c      PRINT *,'INPUT 1 FOR GMT OUTPUT VECTORS, 0 FOR NONE'
c      CALL GFLUSH()
c      READ(*,*) IVECT
      inorm=0
c      PRINT *,'INPUT 1 FOR NORMALIZED VECTORS, 0 FOR ABSOLUTE LENGTHS'
c      CALL GFLUSH()
c      READ(*,*) INORM
      CALL GRSTRT(800,1)                                                
      IF(IPASS.EQ.0) THEN
        CALL WINDOW(REAL(XMIN),REAL(XMAX),REAL(YMIN),REAL(YMAX))
        IPASS=1
      ENDIF
      CALL LINCLR(1)                                                    
      CALL MOVE(REAL(XMIN),REAL(YMIN))                                  
      CALL DRAW(REAL(XMIN),REAL(YMAX))                                  
      CALL DRAW(REAL(XMAX),REAL(YMAX))                                  
      CALL DRAW(REAL(XMAX),REAL(YMIN))                                  
      CALL DRAW(REAL(XMIN),REAL(YMIN))                                  
      CALL DASHPT(3)                                                    
      DO I=1,NOUT                                                       
c       CALL LINCLR(1)                                                    
c       CALL point(REAL(XOUT(I)),REAL(YOUT(I)))                        
        IF(IOUT(I).EQ.1) THEN                                           
          CALL MOVE(REAL(XOUT(I)),REAL(YOUT(I)))                        
        ELSE                                                            
          CALL DRAW(REAL(XOUT(I)),REAL(YOUT(I)))                        
        ENDIF                                                           
      ENDDO                                                             
      CALL CONTR(NMAX,NUMEL,XX,YY,KX,NTYPE,
     &           DEPB,-1.d0,1.d0,1.d0)
231   continue
      read(97,*,iostat=iok) RLATlast,RLONGlast
      rewind 97
      if(iok.eq.0) then
        CALL POLREC(RLATlast,RLONGlast,xlast,ylast)
        xlast=xlast/1000.
        ylast=ylast/1000.
        call rotate(xlast,ylast)
        call linclr(2)
        call point(real(xlast),real(ylast))
        call marker(real(xlast),real(ylast),9)
        CALL GFLUSH()
      else
        print *,'no last point'
      endif
c     call grstop
c     stop
      call gflush
      call gflush
      call gflush
      call gflush
      call gflush
      call gflush
      call gflush
      call gflush
      call gflush
      call gflush
      WRITE(*,*) 'INPUT POSITION OF DESIRED POINT with cursor'
      write(*,*) ' ESC to do manually'
      igot=0
      if(.true.) then
        CALL LOCATE(1,XA,YA,IDAT,IGOT)
      else
        igot=0
      endif
      if(igot.eq.1) then
        XCK=XA(1)
        YCK=YA(1)
        call rotate(yck,xck)
          call recpol(XCK,YCK,RLATCK,RLONGCK)
          print *,RLATCK,RLONGCK
          write(97,*) RLATCK,RLONGCK
        rewind 97
      else
        read(97,*,iostat=iok) RLATCK,RLONGCK
        rewind 97
        if(iok.ne.0) then
          write(*,*) 'INPUT LAT AND LONG OF DESIRED POINT'
          read(*,*) RLATCK,RLONGCK
          write(97,*) RLATCK,RLONGCK
        endif
        write(*,*) RLATCK,RLONGCK
        call polrec(RLATCK,RLONGCK,XCK,YCK)
      endif
c        PRINT *,'INPUT LAT AND LONG OF DESIRED POINT'
        CALL GFLUSH()
c        READ(*,*) RLATCK,RLONGCK
        CALL POLREC(RLATCK,RLONGCK,XCK,YCK)
        xck=xck/1000.
        yck=yck/1000.
        call rotate(xck,yck)
        call point(real(XCK),real(YCK))
        CALL GFLUSH()
        print *,'ok? 1, retry? 0'
        read(*,*) itry
      if(itry.eq.0) goto 231
      print *,'checking point',xck,yck
      write(96,103) '4.2  7.1 20 0 1 1 latitude:',
     &              RLATCK,' longitude:',RLONGCK
      write(98,103) '0.2  5.7 10 0 1 1 latitude:',
     &              RLATCK,' longitude:',RLONGCK
103   format(1x,2(a,f10.3))
      CALL GFLUSH()
      XLEN=5.0                                                          
      YLEN=5.0                                                          
      RNINE=-99999.                                                     
      RMINUS=-1.                                                        
      RZERO=0.                                                          
      RTWO=2.                                                           
      RAROW=20.                                                         
      ICOLOR=1                                                          
      ICSAV=ICOLOR                                                      
      CALL GFLUSH()
      WRITE(*,*) 'INPUT SCALE FACTOR FOR ARROWS, UMAX, AND THRESHOLD'   
      CALL GFLUSH()
      READ(*,*) ASCAL,UMAX,THRESH                                       
      ASCAL1=ASCAL
      UDELTA=UMAX/13.  
      XARO(1)=XCK                                                     
      YARO(1)=YCK                                                     
      IF(INORM.EQ.1) THEN
        XARO(2)=XCK+ASCAL*1.
        YARO(2)=YCK+ASCAL*0.0
      ELSE
        XARO(2)=XCK+ASCAL*UMAX
        YARO(2)=YCK+ASCAL*0.0
      ENDIF
      CALL MOVE(REAL(XARO(1)),REAL(YARO(1)))
      CALL DRAW(REAL(XARO(2)),REAL(YARO(2)))
      IPASS=0
      ISEG=1
C FIND THE ELEMENT IN WHICH OUR POINT RESIDES
      DMIN=1E30
      IFOUND=0
      DO I=1,NUMEL
        DO J=1,4                                                  
          LM(J)=KX(I,J)                               
          XCENT=XCENT+XX(LM(J))                      
          YCENT=YCENT+YY(LM(J))                      
        ENDDO
        XCENT=XCENT/4000.
        YCENT=YCENT/4000.
        DIST=SQRT((XCK-XCENT)**2+(YCK-YCENT)**2)
        IF(DMIN.gt.DIST) THEN
          DMIN=DIST
          IFOUND=I
          if(.false.) then
            XMAX=-1.E30                                                       
            XMIN=1.E30                                                        
            YMAX=-1.E30                                                       
            YMIN=1.E30                                                        
            DO J=1,4
              XMAX=MAX(XMAX,XX(LM(J))*.001)
              XMIN=MIN(XMIN,XX(LM(J))*.001)
              YMAX=MAX(YMAX,YY(LM(J))*.001)
              YMIN=MIN(YMIN,YY(LM(J))*.001)
            ENDDO
          endif
        ENDIF
      ENDDO   
      if(.false.) then
        XBORD=(XMAX-XMIN)/10.                                             
        XMIN=XMIN-XBORD                                                   
        XMAX=XMAX+XBORD                                                   
        YBORD=(YMAX-YMIN)/10.                                             
        YMIN=YMIN-YBORD                                                   
        YMAX=YMAX+YBORD                                                   
        XWID=XMAX-XMIN
        YWID=YMAX-YMIN
        IF(XWID.GT.YWID) THEN
          YWID=XWID
          YMAX=YMIN+YWID
        ELSE
          XWID=YWID
          XMAX=XMIN+XWID
        ENDIF
      endif
      PRINT *,IFOUND 
      PRINT *,xx(IFOUND)/1000.,yy(ifound)/1000.
      print *,xck,yck
      print *,xmin,xmax
      print *,ymin,ymax    
C END OF ELEMENT FINDER ...................... 
1     CONTINUE                                                          
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME                         
      READ(33,END=999) HED                                              
      READ(33) (HTICE(I),I=1,NUMNP)                                     
      READ(33) (ADOT(I),I=1,NUMNP)                                      
      READ(33) (DEPB(I),I=1,NUMNP)                                      
      READ(33) (CONST(I),I=1,NUMEL)                                     
      READ(33) (ACON,I=1,NUMEL)                                     
      READ(34) (FRACT(I),I=1,NUMNP)
      READ(34) (FLOWA(I),I=1,NUMNP)
      READ(34) (SLDGB(I),I=1,NUMNP)
      READ(34) (AFUDGE,I=1,NUMNP)
      if(.false.) then
        DO NUM=1,NUMNP                                                    
          THICK(NUM)=HTICE(NUM)-DEPB(NUM)                                     
        ENDDO                                                             
      else
        DO NUM=1,NUMNP                                                    
          IF(DEPB(NUM).LT.0.) THEN
            FLOT=(1.-RHOW/RHOI)*DEPB(NUM)
            THICK(NUM)=HTICE(NUM)-FLOT
            SURF=0.
c           SURF=PSURF(NUM)
          ELSE
            THICK(NUM)=HTICE(NUM)-DEPB(NUM)
            SURF=DEPB(NUM)
          ENDIF
        ENDDO                                                             
      endif
23    CONTINUE                                                          
      if(ldraw) CALL OPNSEG(ISEG)
      IF(IPASS.EQ.0) THEN
        CALL WINDOW(REAL(XMIN),REAL(XMAX),REAL(YMIN),REAL(YMAX))
        IPASS=1
      ENDIF
      if(.false.) then
        CALL CONTR(NMAX,NUMEL,XX,YY,KX,NTYPE,
     &             HTICE,-2999.999D0,5500.D0,500.D0)
      endif
      if(.false.) then
        CALL CONTR(NMAX,NUMEL,XX,YY,KX,NTYPE,
     &             thick,1d0,1.5d0,1d0)
      endif
      if(.false.) then
        CALL CONTR(NMAX,NUMEL,XX,YY,KX,NTYPE,
     &             FRACT,-0.001D0,0.9D0,0.5D0)
      endif
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
C     REWIND 9                                                          
C     WRITE(9,*) HED
C     WRITE(9,*) 'ELEMENT,X,Y,UX,UY,UMAG'                               
      VMAX=-1.E30                                                       
      if(.true.) then
        CALL MOVE(REAL(XMIN+(XMAX-XMIN)*.05),REAL(YMAX-(YMAX-YMIN)*.05))
        CALL TXTCLR(1)
        CALL TEXT(80,HED)                                               
      endif
      read(hed,'(6x,f12.0)') time
      DO 600 J=IFOUND,IFOUND                                                  
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
        AFRACT=0.0D00                                                       
        AFLOWA=0.0D00                                                       
        ASLDGB=0.0D00                                                       
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
        CALL linclr(1)
        CALL point(real(xcent),real(ycent))
C                                                                       
        CALL SHAPEFEM(NTYPE(J),CENTX,CENTY,PSI,DPSI)                         
      
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
C       IF (DETJ.GE.0.0) GOTO 750                                         
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
          AFLOWA=AFLOWA + FLOWA(LM(I))*PSI(I)                                       
          ASLDGB=ASLDGB + SLDGB(LM(I))*PSI(I)                                       
          AFRACT=AFRACT + FRACT(LM(I))*PSI(I)                                       
        ENDDO                                                          
C                                                                       
        DELH=SUMX**2 + SUMY**2                                            
        DELH=SQRT(DELH)                                                   
        IF(HH.GT.HMIN) THEN                                               
          IF(IVELO.EQ.1) THEN
C ....... JUST THE SLIDING VELOCITY, BASED ON FRACT ....
            TERM1 = AFRACT*((PG/ASLDGB)**2)*(HH**3)*DELH
C
C .....     the .2 in here is wrong, but its been that way for a long time
C .....     the reduction in FLOW that results in the same flux is .79. humm..
c           TERM2=(1.-AFRACT)*(0.2)*((PG/AFLOWA)**3)*(HH**5)*(DELH**2)
            TERM2=(1.-AFRACT)*(0.4)*((PG/AFLOWA)**3)*(HH**5)*(DELH**2)
            CNEW = TERM1
            UX=-CNEW*SUMX/HH                                            
            UY=-CNEW*SUMY/HH                                            
          else
            UX=-CONST(J)*SUMX/HH                                            
            UY=-CONST(J)*SUMY/HH                                            
          endif
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
C        WRITE(9,1001) J,XCENT,YCENT,UX,UY,UMAG,CONST(J)
C 1001  FORMAT(I6,6G13.6)                                                 
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
        IF(IVECT.EQ.1) THEN
          WRITE(10,'(A1)') '>'
          CALL RECPOL(XARO(1),YARO(1),RLAT,RLONG)
          WRITE(10,*) REAL(RLONG),REAL(RLAT)
          CALL RECPOL(XARO(2),YARO(2),RLAT,RLONG)
          WRITE(10,*) REAL(RLONG),REAL(RLAT)
          WRITE(13,'(A1)') '>'
          WRITE(13,*) REAL(XARO(1)),REAL(YARO(1))
          WRITE(13,*) REAL(XARO(2)),REAL(YARO(2))
          WRITE(14,*) REAL(XARO(1)),REAL(YARO(1)),REAL(UMAG)
          theta=angle(xcent,ycent,ux,uy)
          write(99,*) real(time),real(theta),real(UMAG)
        ENDIF
c 202   FORMAT(10X,G13.6,2X,G13.6,I13)                                    
C       WRITE(7,202) XARO(1),YARO(1)                                      
C       WRITE(7,202) XARO(2),YARO(2)                                      
C       WRITE(7,202) RNINE,RAROW,ICMAP(ICOLOR)-1                          
C       WRITE(7,*) ' ELEMENT',J                                           
C       IF(ICOLOR.NE.ICSAV) THEN                                          
          CALL LINCLR(ICMAP(ICOLOR))                                      
          ICSAV=ICOLOR                                                    
C       ENDIF                                                             
        CALL MOVE(REAL(XARO(1)),REAL(YARO(1)))
        CALL DRAW(REAL(XARO(2)),REAL(YARO(2)))
600   CONTINUE                                                          
      XARO(1)=XMIN+(XMAX-XMIN)/XLEN                                     
      XARO(2)=XARO(1)+(XMAX-XMIN)/XLEN                                  
      YARO(1)=YMAX+(YMAX-YMIN)/YLEN                                     
      YARO(2)=YARO(1)                                                   
      CALL MOVE(REAL(XARO(1)),REAL(YARO(1)))
      CALL DRAW(REAL(XARO(2)),REAL(YARO(2)))
      if(ldraw) CALL CLOSEG
      IF(ISEG.NE.1) then
        if(ldraw) CALL SETVIS(ISEG-1,.FALSE.)
      endif
      ISEG=ISEG+1
      ZSCAL=(XMAX-XMIN)/ASCAL/XLEN                                      
C     WRITE(7,202) XARO(1),YARO(1)                                      
C     WRITE(7,202) XARO(2),YARO(2)                                      
C     WRITE(7,202) XMIN,YMIN                                            
C     WRITE(7,202) XMIN,YMIN                                            
C     WRITE(7,202) RNINE,RAROW+1.                                       
C     WRITE(7,2323) ZSCAL                                               
C 2323  FORMAT(G10.3)                                                     
C     CALL GRSTOP1                                                      
      if(vmax.gt.0) WRITE(*,2324) HED,' UMAX=',VMAX
2324  FORMAT(A20,A,4G13.6)
      IF(ASCAL1.GT.0.) THEN
        CALL GFLUSH()
        READ(*,*,END=24) ASCAL1,UMAX,THRESH 
      ENDIF
      UDELTA=UMAX/13.                                                   
      IF(ASCAL1.LE.0.) GOTO 24                                          
      ASCAL=ASCAL1                                                      
      REWIND 7                                                          
      REWIND 10                                                          
C     WRITE(7,105) XMIN,YMIN
C     WRITE(7,105) XMAX,YMIN
C     WRITE(7,105) XMAX,YMAX
C     WRITE(7,105) XMIN,YMIN
C     WRITE(7,105) XMIN,YMIN
C     WRITE(7,105) RNINE,RTWO                                         
C     WRITE(7,111)                                                    
C     CALL GRSTRT(800,1)                                                
      GOTO 23                                                           
24    CONTINUE                                                          
      GOTO 1                                                            
999   CONTINUE                                                          
      NSEG=ISEG-1
900   CONTINUE
      if(.false.) then
        if(ldraw) CALL SETVIS(NSEG,.FALSE.)
        CALL GFLUSH()
        CALL GETUIN(-1,'NO PAUSE 1, PAUSE 0',1,INARAY,IGOT)
        IF(IGOT.NE.1) THEN
          IPAUSE=1
        ELSE
          IPAUSE=INARAY(1)
        ENDIF
        CALL GFLUSH()
        if(ldraw) CALL SETVIS(1,.TRUE.)
        DO 910 I=2,NSEG
          if(ldraw) CALL SETVIS(I,.TRUE.)
          IF(IPAUSE.EQ.0) CALL GETUIN(-1,'CONTINUE?',1,INARAY,IGOT)
          if(ldraw) CALL SETVIS(I-1,.FALSE.)
C         CALL NEWPAG
          IF(IPAUSE.EQ.0) CALL GETUIN(-1,'CONTINUE?',1,INARAY,IGOT)
910     CONTINUE
        CALL GFLUSH()
        CALL GETUIN(-1,'AGAIN? 1, FINISH? 0',1,INARAY,IGOT)
        if(ldraw) CALL SETVIS(NSEG,.FALSE.)
        IF(IGOT.EQ.1) THEN
          IF(INARAY(1).EQ.1) GOTO 900
        ENDIF
911     CONTINUE
        if(ldraw) CALL SETVIS(-1,.TRUE.)
        PRINT *,'INPUT 1 TO DO AGAIN, 0 to finish, -1 to save picture'
        CALL GFLUSH()
        READ(*,*) IGOT
        IF(IGOT.EQ.1) THEN
          CALL NEWPAG
          GOTO 911
        ENDIF
c       if(igot.eq.-1) call savescreen(1,800,800)
      endif
      CALL GRSTOP1                                                      
      STOP                                                              
 750  CONTINUE                                                          
      WRITE (*,5544) DETJ,X                                             
 5544 FORMAT (' BAD JACOBIAN AT 750',E14.6,E14.6)                       
      STOP                                                              
      END                                                               
C******************************************************
      SUBROUTINE SHAPEFEM(NTYPE,XI,ET,PSI,DPSI)                            
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
C====================================
      SUBROUTINE CONTR(NMAX,NUMEL,XX,YY,KX,NTYPE,
     &                 ZZ,RMIN,RMAX,RLEVSP)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION ZZ(NMAX),XX(NMAX),YY(NMAX),KX(NMAX,4),IK(5)
      DIMENSION NTYPE(NMAX)
      LOGICAL FOUND
      DIMENSION ICMAP(16)
      DATA ICMAP /0,4,12,6,13,2,8,7,9,3,10,5,11,1,14,15/                
C     DATA ICMAP /0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/                
C DRAW CONTOURS
      LEVEL=(RMAX-RMIN)/RLEVSP+1
      DLEV=RLEVSP
      IF(LEVEL.LE.0) THEN
        LEVEL=2-LEVEL
        DLEV=-DLEV
      ENDIF
      ICO=1
      CALL LINCLR(ICO)
      LCONT=1
      XPOS=XMAX-(XMAX-XMIN)*.125
      YPOS=YMAX-(YMAX-YMIN)*.1
      ICO=1
      CALL LINCLR(ICO)
      DO LEV=1,LEVEL
        FOUND=.FALSE.
        ICOUNT=0
        VAL=RMIN+(LEV-1)*DLEV
        DO N=1,NUMEL
          LCONT=1
          IF(NTYPE(N).EQ.1) THEN                                            
            NNODE=4                                                         
          ELSE                                                              
            NNODE=3                                                         
          ENDIF                                                             
          DO J=1,NNODE
            IK(J)=KX(N,J)
          ENDDO
          IK(NNODE+1)=KX(N,1)
          ICOUNT=0
          DO NN=1,NNODE
            I=IK(NN)
            J=IK(NN+1)
            IF(ZZ(I).LT.VAL .AND. VAL.LE.ZZ(J) .OR. 
     &         ZZ(J).LT.VAL .AND. VAL.LE.ZZ(I)) THEN
              FOUND=.TRUE.
              ICOUNT=ICOUNT+1
              XINT=XX(J)+(VAL-ZZ(J))*(XX(I)-XX(J))/(ZZ(I)-ZZ(J))
              YINT=YY(J)+(VAL-ZZ(J))*(YY(I)-YY(J))/(ZZ(I)-ZZ(J))
              XINT=XINT*.001
              YINT=YINT*.001
              IF(LCONT.EQ.1) THEN
                CALL MOVE(REAL(XINT),REAL(YINT))
                LCONT=2
              ELSE
                CALL DRAW(REAL(XINT),REAL(YINT))
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        IF(FOUND) THEN
          CALL MOVE(REAL(XPOS),REAL(YPOS))
          CALL TXTCLR(ICMAP(ICO))
          CALL RNUMBR(REAL(VAL),3,9)
          YPOS=YPOS-(YMAX-YMIN)*.025
        ENDIF
        ICO=ICO+1
        IF(ICO.GT.6) ICO=1
        CALL LINCLR(ICO)
      ENDDO
      END
      SUBROUTINE rotate(xx,yy)
      IMPLICIT REAL*8(A-H,O-Z)  
      return
      tt=xx
      xx=yy
      yy=-tt
      END
      REAL*8 FUNCTION angle(xcent,ycent,ux,uy)
      IMPLICIT REAL*8(A-H,O-Z)  
      PI=4.D0*ATAN(1.D0)
      dot=(xcent*ux+ycent*uy)
      cross=(xcent*uy-ycent*ux)
      xmag=sqrt(xcent**2+ycent**2)
      umag=sqrt(UX**2+UY**2)
      if(UMAG.gt.1) then
        theta=dot/xmag/umag
        theta=180.-acos(theta)*180./pi
        theta=sign(theta,cross)
        if(theta.lt.0) theta=theta+360.
      else
        theta=-999.
      endif
      angle=theta
      END
