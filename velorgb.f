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
      include "parameter.h"                                     
      PARAMETER(NMAX=MAXNUM,HMIN=10.,NPLOT=20000)                           
      CHARACTER HED*80                                                  
      REAL*4  PX(4),PY(4)
      DIMENSION ICMAP1(16)
c      DATA ICMAP1 /0,4,12,6,13,2,8,7,9,3,10,5,11,1,14,15/                
      DATA ICMAP1 /0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/                
      DIMENSION ICMAP(16),INARAY(2)
      DIMENSION XARO(9),YARO(9),NBOUND(NMAX)                 
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2),LM(4)             
      DIMENSION PSI(4),DPSI(4,2)                                        
      DIMENSION XY(2,4)
      DIMENSION XX(NMAX),YY(NMAX),HTICE(NMAX),THICK(NMAX),NTYPE(NMAX)
      DIMENSION TBED(NMAX),WTHICK(NMAX),BMELT(NMAX)      
      DIMENSION XLINE(NMAX),YLINE(NMAX)                                   
      DIMENSION KX(NMAX,4),CONST(NMAX)
      DIMENSION PSURF(NMAX),FRACT(NMAX),FLOWA(NMAX),SLDGB(NMAX)             
      DIMENSION ADOT(NMAX),DEPB(NMAX),KODE(NMAX),AJUNK(NMAX)
      DATA ICMAP /0,4,12,6,13,2,8,7,9,3,10,5,11,1,14,14/                
C     DATA ICMAP /0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/                
      DIMENSION XOUT(NPLOT),YOUT(NPLOT),IOUT(NPLOT) 
      logical iflush
      common /flush/ iflush
      data iflush /.false./
      COMMON /MAXES/ XMAX,XMIN,YMAX,YMIN                    
      CALL SETRIG
      icount=0
      print *,'input start number,'
      print *,'      display interval, and' 
      print *,'      how many to show'
      read(*,*) istart,idisp,inumber
      PRINT *,'INPUT 1 FOR NORMALIZED VECTORS, 0 FOR ABSOLUTE LENGTHS'
      READ(*,*) INORM
      XLEN=5.0                                                          
      YLEN=5.0                                                          
      ICOLOR=1                                                          
      II=1                                                              
      DO I=1,NPLOT                                                      
        READ(12,*,END=21) XIN,YIN,IIN                                   
        IF(XIN.EQ.-99999.) THEN                                         
          READ(12,100) JUNK                                              
100       FORMAT(A80)                                                       
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
      READ(*,*) ASCAL,UMAX,THRESH                                       
      ASCAL1=ASCAL
      UDELTA=UMAX/13.                                                   
      CALL GRSTRT(600,600)                                                
C READ INPUT HEADER                                                     
      READ(30,1000,END=991) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,   
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
        DO I=1,NLINE                                                 
          XLINE(I)=XX(NBOUND(I))/1000.                                    
          YLINE(I)=YY(NBOUND(I))/1000.                                    
          XMAX=MAX(XMAX,XLINE(I))
          XMIN=MIN(XMIN,XLINE(I))
          YMAX=MAX(YMAX,YLINE(I))
          YMIN=MIN(YMIN,YLINE(I))
        ENDDO                                                        
      ELSE                                                              
        DO I=1,NUMNP                                                
          XMAX=MAX(XMAX,XX(I)*.001)
          XMIN=MIN(XMIN,XX(I)*.001)
          YMAX=MAX(YMAX,YY(I)*.001)
          YMIN=MIN(YMIN,YY(I)*.001)
        ENDDO                                                        
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
      IPASS=0
      IPR=1
      idstep=1
      print *,'input base number,step:',ipr,idstep
      read(*,*) ipr,idstep
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME                         
c bypass the first few ...
      do l=1,istart-1
        READ(33,END=992) HED           
        print *,'bypassing',hed                                   
        READ(33) (HTICE(I),I=1,NUMNP)                                     
        READ(33) (ADOT(I),I=1,NUMNP)                                      
        READ(33) (DEPB(I),I=1,NUMNP)                                      
        READ(33) (CONST(I),I=1,NUMEL)                                     
        READ(33) (ACON,I=1,NUMEL)                                     
        READ(34) (FRACT(I),I=1,NUMNP)
        READ(34) (FLOWA(I),I=1,NUMNP)
        READ(34) (SLDGB(I),I=1,NUMNP)
        READ(34) (AFUDGE,I=1,NUMNP)
        READ(36,END=993) HED
        READ(36) (TBED(I),I=1,NUMNP)
        READ(36) (BMELT(I),I=1,NUMNP)
        READ(36) (WTHICK(I),I=1,NUMNP)
      enddo
1     CONTINUE                                                          
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME                         
c ... now start drawing them...
      READ(33,END=994) HED                                              
      READ(33) (HTICE(I),I=1,NUMNP)                                     
      READ(33) (ADOT(I),I=1,NUMNP)                                      
      READ(33) (DEPB(I),I=1,NUMNP)                                      
      READ(33) (CONST(I),I=1,NUMEL)                                     
      READ(33) (ACON,I=1,NUMEL)                                     
      READ(34) (FRACT(I),I=1,NUMNP)
      READ(34) (FLOWA(I),I=1,NUMNP)
      READ(34) (SLDGB(I),I=1,NUMNP)
      READ(34) (AFUDGE,I=1,NUMNP)
      READ(36,END=995) HED
      READ(36) (TBED(I),I=1,NUMNP)
      READ(36) (BMELT(I),I=1,NUMNP)
      READ(36) (WTHICK(I),I=1,NUMNP)
      icount=icount+1
      if(icount.gt.inumber*idisp) stop
      if(1+mod(icount,idisp).ne.1) goto 1
        print *,'outputting for display',hed
c--------------------------------------      
      DO NUM=1,NUMNP                                                    
        THICK(NUM)=HTICE(NUM)-DEPB(NUM)                                     
      ENDDO                                                             
23    CONTINUE                                                          
      IF(IPASS.EQ.0) THEN
        CALL WINDOW(REAL(XMIN),REAL(XMAX),REAL(YMIN),REAL(YMAX))
        IPASS=1
      ENDIF
      print *,REAL(XMIN),REAL(XMAX),REAL(YMIN),REAL(YMAX)
c
      do i=1,numel
        favg=0
        wavg=0
        DO L=1,4
          PX(L)=XX(KX(I,L))/1000.
          PY(L)=YY(KX(I,L))/1000.
          favg=favg+FRACT(KX(I,L))
          wavg=wavg+WTHICK(KX(I,L))
        ENDDO
        favg=favg*0.25
        wavg=wavg*0.25
        favg=4*wavg+0.49
        if(.true.) then
          if(favg.gt.0.1) then
            icmin= 0
            aaa=(icmin-145)/2.0
            bbb=icmin-aaa
            CALL FILPAN(ICMAP1(15),.false.)
            iccc=aaa*favg+bbb
            iccc=izset(200,favg,0.d0,1d-2,0)
            call linclrg(iccc)
            CALL PANEL(4,PX,PY)  
          endif
        else
          if(favg.gt.0.5) then
            CALL FILPAN(ICMAP1(15),.false.)
            CALL PANEL(4,PX,PY)
          endif
        endif
      enddo
      call offmap
c
      CALL MOVE(REAL(XMIN+(XMAX-XMIN)*.05),REAL(YMAX-(YMAX-YMIN)*.05))
      CALL TXTCLR(1)
      CALL TEXT(80,HED)                                               
      if(.false.) then
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
        call savescreen(ipr+1000,600,600)                                                         
c       goto 1234
      endif
      CALL CONTR(NMAX,NUMEL,XX,YY,KX,NTYPE,
     &                 HTICE,-2999.999D0,5500.D0,250.D0,1)
c      CALL CONTR(NMAX,NUMEL,XX,YY,KX,NTYPE,
c     &                 FRACT,-0.001D0,0.9D0,0.5D0,0)
c      call minmax(numnp,wthick,wmin,wmax)
c      wlev=(wmax-wmin)/15.
c      print *,wmin,wmax,wlev
c      CALL CONTR(NMAX,NUMEL,XX,YY,KX,NTYPE,
c     &                 WTHICK,0.0d0,0.5d0,0.05d0,0)
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
73      CONTINUE                                                        
      ENDIF                                                             
      VMAX=-1.E30                                                       
      CALL MOVE(REAL(XMIN+(XMAX-XMIN)*.05),REAL(YMAX-(YMAX-YMIN)*.05))
      CALL TXTCLR(1)
      CALL TEXT(80,HED)                                               
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
        icolor=izset(13,umag,0.d0,udelta,1)
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
        CALL LINCLR(ICMAP(ICOLOR))                                      
        CALL MOVE(REAL(XARO(1)),REAL(YARO(1)))
        CALL DRAW(REAL(XARO(2)),REAL(YARO(2)))
600   CONTINUE                                                          
      XARO(1)=XMIN+(XMAX-XMIN)*.025                                     
      XARO(2)=XARO(1)+(XMAX-XMIN)* .12                                  
      YARO(1)=YMIN+(YMAX-YMIN)*.025                                      
      YARO(2)=YARO(1)  
      YPOS=YARO(2)
      XPOS=XARO(1)
      j=1
      CALL TXTCLR(1)
      CALL MOVE(REAL(xmin),REAL(YPOS))
      CALL RNUMBR(REAL(UDELTA*(j-1)),3,9)
      YPOS=YPOS+(YMAX-YMIN)*.025
      do j=2,14
        CALL LINCLR(ICMAP(j))                                      
c        print *,real(xmin),real(xmax)
c        print *,real(ymin),real(ymax)
c        print *,REAL(XARO(1)),REAL(YPOS)
c        print *,REAL(XARO(2)),REAL(YPOS)
        CALL MOVE(REAL(XARO(1)),REAL(YPOS))
        CALL DRAW(REAL(XARO(2)),REAL(YPOS))
        YPOS=YPOS+(YMAX-YMIN)*.01
        CALL TXTCLR(1)                                      
        CALL MOVE(REAL(XPOS),REAL(YPOS))
        CALL RNUMBR(REAL(UDELTA*(j-1)),3,9)
        YPOS=YPOS+(YMAX-YMIN)*.025
      enddo
      if(iflush) call gflush
      ZSCAL=(XMAX-XMIN)/ASCAL/XLEN                                      
      WRITE(*,2324) HED,' UMAX=',VMAX                                           
2324  FORMAT(A20,A,G13.6)
      IF(ASCAL1.GT.0.) READ(*,*) ASCAL1,UMAX,THRESH 
      UDELTA=UMAX/13.                                                   
      IF(ASCAL1.GT.0.) THEN
        ASCAL=ASCAL1                                                      
        GOTO 23                                                           
      ENDIF 
      call savescreen(ipr,600,600)                                                         
1234  continue
      CALL NEWPAG
      IPR=IPR+idstep
      GOTO 1                                                            
 750  CONTINUE                                                          
      WRITE (*,5544) DETJ,X                                             
 5544 FORMAT (' BAD JACOBIAN AT 750',E14.6,E14.6)                       
      CALL GRSTOP1                                                      
      STOP                                                              
991   CONTINUE                                                          
      print *,'end-of-file unit 30 encountered...'
      CALL GRSTOP1          
      stop                                            
992   CONTINUE                                                          
      print *,'end-of-file unit 33 encountered...'
      CALL GRSTOP1          
      stop                                            
993   CONTINUE                                                          
      print *,'end-of-file unit 36 encountered...'
      CALL GRSTOP1          
      stop                                            
994   CONTINUE                                                          
      print *,'end-of-file unit 33 encountered...(second)'
      CALL GRSTOP1          
      stop                                            
995   CONTINUE                                                          
      print *,'end-of-file unit 36 encountered...(second)'
      CALL GRSTOP1          
      stop                                            
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
     &                 ZZ,RMIN,RMAX,RLEVSP,ILABEL)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION ZZ(NMAX),XX(NMAX),YY(NMAX),KX(NMAX,4),IK(5)
      DIMENSION NTYPE(NMAX)
      LOGICAL FOUND
      DIMENSION ICMAP(16)
      COMMON /MAXES/ XMAX,XMIN,YMAX,YMIN                    
      DATA ICMAP /0,4,12,6,13,2,8,7,9,3,10,5,11,1,14,14/                
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
      XPOS=XMAX-(XMAX-XMIN)*.13
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
        IF(FOUND .AND. ILABEL.EQ.1) THEN
          CALL MOVE(REAL(XPOS),REAL(YPOS))
          CALL TXTCLR(ICO)
          CALL RNUMBR(REAL(VAL),3,9)
          YPOS=YPOS-(YMAX-YMIN)*.025
        ENDIF
        ICO=ICO+1
        IF(ICO.GT.6) ICO=1
        CALL LINCLR(ICO)
      ENDDO
      END
      subroutine minmax(n,x,xmin,xmax)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      dimension x(n)
      xmin=1e30
      xmax=-1e30
      do i=1,n
        xmin=min(xmin,x(i))
        xmax=max(xmax,x(i))
      enddo
      end

C===========================================
      FUNCTION IZSET(NCOLOR,ZZZ,ZF,ZD,IOFF)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      IMAX=NCOLOR+IOFF
      IMIN=IOFF+1
C       IC=IMAX-NINT((ZZZ-ZF)/ZD-.5d0)
        ZIC=(ZZZ-ZF)/ZD+.5d0
        IC=NINT(ZIC)
C       WRITE(*,*) ZZZ,ZIC,IC
      IF(IC.LT.IMIN) IC=IMIN
      IF(IC.GT.IMAX) IC=IMAX
      IZSET=IC
      RETURN
      END

