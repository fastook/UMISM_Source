C **************************************************************C       
C                                                               C       
C   PROGRAM:  VELOEXTRACT                                       C       
C                                                               C       
C   DATE:  04 18 23                                             C       
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
      PARAMETER(NMAX=29999,HMIN=1.)                           
      CHARACTER HED*80                                                  
      REAL*4  PX(4),PY(4)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2),LM(4)             
      DIMENSION PSI(4),DPSI(4,2),XY(2,4)
      DIMENSION XX(NMAX),YY(NMAX),HTICE(NMAX),THICK(NMAX),NTYPE(NMAX)        
      DIMENSION KX(NMAX,4),CONST(NMAX)
      DIMENSION PSURF(NMAX),FRACT(NMAX),FLOWA(NMAX),SLDGB(NMAX)             
      DIMENSION ADOT(NMAX),DEPB(NMAX),KODE(NMAX),AJUNK(NMAX)
C READ INPUT HEADER                                                     
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
        DO NUM=1,NUMNP                                                    
          THICK(NUM)=HTICE(NUM)-DEPB(NUM)                                     
        ENDDO                                                             
        VAVG=0; NUMVELO = 0
        HAVG=0; NUMTHICK = 0
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
c         IF (DETJ.GE.0.0) GOTO 750                                         
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
              HAVG = HAVG + HH
              NUMTHICK = NUMTHICK + 1
          ELSE                                                              
            UX=0.                                                           
            UY=0.                                                           
          ENDIF                                                             
          UMAG=SQRT(UX**2+UY**2)
          IF(UMAG.GT.0) THEN
            VAVG = VAVG + UMAG
            NUMVELO = NUMVELO + 1
          ENDIF
          VMAX=MAX(VMAX,UMAG)                                               
c         WRITE(9,1001) J,XCENT,YCENT,UX,UY,UMAG,CONST(J)
1001  FORMAT(I6,6G13.6)                                                 
202   FORMAT(10X,G13.6,2X,G13.6,I13)                                    
600     CONTINUE                                                          
        IF(NUMVELO.GT.0) VAVG = VAVG / NUMVELO
        IF(NUMTHICK.GT.0) HAVG = HAVG / NUMTHICK
        read(HED,2325) ITIME
2325    format(10x,i10)
        TIME = ITIME/1000.
        WRITE(*,2324) HED(10:20),' UMAX=',VMAX,VAVG,NUMVELO
        write(11,*) TIME,VMAX
        write(12,*) TIME,VAVG
        write(13,*) TIME,NUMVELO
        write(14,*) TIME,HAVG
2324  format(a,a,2g13.6,i10)
c       WRITE(9,*) HED
c       WRITE(9,*) 'ELEMENT,X,Y,UX,UY,UMAG'                               
        !GOTO 23                                                           
24      CONTINUE                                                          
      GOTO 1                                                            
999   CONTINUE                                                          
      write(11,*) -99999,0
      write(11,*) 'VMAX'
      write(12,*) -99999,0
      write(12,*) 'VAVG'
      write(13,*) -99999,0
      write(13,*) 'NUM NONZERO'
      write(14,*) -99999,0
      write(14,*) 'HAVG'
      STOP                                                              
 750  CONTINUE                                                          
      WRITE (*,5544) DETJ,X                                             
 5544 FORMAT (' BAD JACOBIAN at 750',E14.6,E14.6)                       
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
