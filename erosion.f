      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER HED*80,HED1*80,JUNK*80
      DIMENSION KODE(NMAX),X(NMAX),Y(NMAX),HTICE(NMAX),ADOT(NMAX),
     &          FRACT(NMAX),PSURF(NMAX),BDROCK(NMAX),FLOWA(NMAX),
     &          SLDGB(NMAX),KX(NMAX,4),CONST(NMAX),BED(NMAX),
     &          IBFLUX(NMAX,2),BFLUX(NMAX),NTYPE(NMAX),
     &          TBED(NMAX),PTEMP(NMAX),ATIME(NMAX),THICK(NMAX),
     &          ACON(NMAX),ptemp1(nmax),ptemp2(nmax),ttime(nmax)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2),LM(4)             
      DIMENSION PSI(4),DPSI(4,2)                                        
      DIMENSION XY(2,4),XA(1),YA(1),IDAT(1)
      data ntype /nmax*1/
      data ttime /nmax*0.d0/
      HMIN=0.1
      print *,'input start time, stop time'
      read(*,*) tstart, tstop
      IVELO=0
c      PRINT *,'TOTAL VELOCITY (0) OR JUST SLIDING (1)'
c      READ(*,*) IVELO
      THRESH=1000.
      PG = 0.089866
      RHOW=1.092
      RHOI=0.917
C READ INPUT HEADER
      READ(30,1000) HED1,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      PRINT *,HED1
1000  FORMAT (A80,/,7I6,F8.0)
      IF(NUMNP.GT.NMAX) THEN
        PRINT *,'NUMNP=',NUMNP,' NMAX=',NMAX,' INCREASE NMAX'
        STOP
      ENDIF
      DO I=1,NUMNP
        ATIME(I)=0.
      ENDDO
C READ INPUT GRID, THINGS THAT NEVER CHANGE
      READ(31) HED
      PRINT *,HED
      READ(31) (KODE(I),I=1,NUMNP)
      READ(31) (X(I),I=1,NUMNP)
      READ(31) (Y(I),I=1,NUMNP)
      READ(31) (PSURF(I),I=1,NUMNP)
      READ(31) (BED(I),I=1,NUMNP)
      READ(31) (KX(I,1),KX(I,2),KX(I,3),KX(I,4),I=1,NUMEL)
      READ(31) (IBFLUX(I,1),IBFLUX(I,2),BFLUX(I),I=1,NUMGBC)
C READ INPUT DIFFER, THINGS THAT DIFFER FROM ONE GRID TO THE NEXT
      READ(32) HED
      PRINT *,HED
      READ(32) (ADOT(I),I=1,NUMNP)
      READ(32) (FRACT(I),I=1,NUMNP)
      READ(32) (FLOWA(I),I=1,NUMNP)
      READ(32) (SLDGB(I),I=1,NUMNP)
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME
      ICOUNT=0
c ... read and discard first output as times might be wrong ...
      READ(33,END=999) HED
      WRITE(JUNK,*) HED
      READ(JUNK,1010) TIME
      PRINT *,HED(1:20),TIME
      READ(33) (HTICE(I),I=1,NUMNP)
      READ(33) (ADOT(I),I=1,NUMNP)
      READ(33) (BDROCK(I),I=1,NUMNP)
      READ(33) (CONST(I),I=1,NUMEL)
      READ(33) (ACON(I),I=1,NUMEL)
      READ(34) (FRACT(I),I=1,NUMNP)
      READ(34) (FLOWA(I),I=1,NUMNP)
      READ(34) (SLDGB(I),I=1,NUMNP)
      READ(34) (AFUDGE,I=1,NUMNP)
      if(.false.) then
        DO NUM=1,NUMNP                                                    
          THICK(NUM)=HTICE(NUM)-BDROCK(NUM)                                     
        ENDDO                                                             
      else
        DO NUM=1,NUMNP                                                    
          IF(BDROCK(NUM).LT.0.) THEN
            FLOT=(1.-RHOW/RHOI)*BDROCK(NUM)
            THICK(NUM)=HTICE(NUM)-FLOT
            SURF=0.
c           SURF=PSURF(NUM)
          ELSE
            THICK(NUM)=HTICE(NUM)-BDROCK(NUM)
            SURF=BDROCK(NUM)
          ENDIF
        ENDDO                                                             
      endif
      READ(36,END=999) HEDT
      READ(36) (TBED(I),I=1,NUMNP)
      READ(36) (BMELT,I=1,NUMNP)
      READ(36) (WTHICK,I=1,NUMNP)
c ... now start reading for real ...
85    READ(33,END=999) HED
        WRITE(JUNK,*) HED
        READ(JUNK,1010) TIME
        READ(33) (HTICE(I),I=1,NUMNP)
        READ(33) (ADOT(I),I=1,NUMNP)
        READ(33) (BDROCK(I),I=1,NUMNP)
        READ(33) (CONST(I),I=1,NUMEL)
        READ(33) (ACON(I),I=1,NUMEL)
        READ(34) (FRACT(I),I=1,NUMNP)
        READ(34) (FLOWA(I),I=1,NUMNP)
        READ(34) (SLDGB(I),I=1,NUMNP)
        READ(34) (AFUDGE,I=1,NUMNP)
        if(time.lt.tstart) goto 85
        if(time.gt.tstop) goto 999
        PRINT *,HED(1:20),TIME
        ICOUNT=ICOUNT+1
        if(.false.) then
          DO NUM=1,NUMNP                                                    
            THICK(NUM)=HTICE(NUM)-BDROCK(NUM)                                     
          ENDDO                                                             
        else
          DO NUM=1,NUMNP                                                    
            IF(BDROCK(NUM).LT.0.) THEN
              FLOT=(1.-RHOW/RHOI)*BDROCK(NUM)
              THICK(NUM)=HTICE(NUM)-FLOT
              SURF=0.
c             SURF=PSURF(NUM)
            ELSE
              THICK(NUM)=HTICE(NUM)-BDROCK(NUM)
              SURF=BDROCK(NUM)
            ENDIF
          ENDDO                                                             
        endif
        READ(36,END=999) HEDT
        READ(36) (TBED(I),I=1,NUMNP)
        READ(36) (BMELT,I=1,NUMNP)
        READ(36) (WTHICK,I=1,NUMNP)
C--------------------------------------------
        IF(ICOUNT.EQ.1) THEN
          TBASE=TIME
          TOLD=TIME
        ELSE
          DT=TIME-TOLD
          DO J=1,NUMEL                                                  
           ATIME(j)=ATIME(j)+DT
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
            XY(1,1)=X(I)                                                     
            XY(1,2)=X(JJ)                                                    
            XY(1,3)=X(K)                                                     
            IF(NTYPE(J).EQ.1) XY(1,4)=X(L)                                   
            XY(2,1)=Y(I)                                                     
            XY(2,2)=Y(JJ)                                                    
            XY(2,3)=Y(K)                                                     
            IF(NTYPE(J).EQ.1) XY(2,4)=Y(L)                                   
            IF(NTYPE(J).EQ.1) THEN                                            
              XCENT=(X(I)+X(JJ)+X(K)+X(L))/4000.                          
              YCENT=(Y(I)+Y(JJ)+Y(K)+Y(L))/4000.                          
            ELSE                                                              
              XCENT=(X(I)+X(JJ)+X(K))/3000.                                
              YCENT=(Y(I)+Y(JJ)+Y(K))/3000.                                
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
c           IF (DETJ.GE.0.0) GOTO 750                                         
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
              ttime(j)=ttime(j)+dt
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
            IF(UMAG.LT.THRESH .and. umag.gt.0) then
              ptemp1(j)=ptemp1(j)+umag*dt
c             PRINT *,j,REAL(PTEMP1(j)),REAL(ATIME(j)),real(dt),
c    &                   real(umag)
             endif
          ENDDO
          TOLD=TIME
        ENDIF
      GOTO 85
999   CONTINUE
      DO I=1,numel
        IF(ATIME(I).NE.0.0) THEN
          ptemp2(i)=ptemp1(i)/1000
          PTEMP1(I)=PTEMP1(I)/ATIME(I)
          ttime(I)=100*ttime(I)/ATIME(I)
          print *,i,ttime(i)
        ELSE
          PTEMP1(I)=0.0
          ttime(I)=0.0
        ENDIF
      ENDDO
      WRITE(1,1000) HED1,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      DO N=1,NUMNP
c         WRITE(1,1001) N,KODE(N),X(N),Y(N),HTICE(N),
          WRITE(1,*) N,KODE(N),X(N),Y(N),HTICE(N),
     &                  max(ADOT(N),-999.),FRACT(N),PSURF(N),
     &                  BDROCK(N),FLOWA(N),SLDGB(N),0.,0,0.,0.
      ENDDO
      DO N=1,NUMEL
        WRITE(1,1002) N,KX(N,1),KX(N,2),KX(N,3),KX(N,4),PTEMP1(N),
c    &                PTEMP2(N)
     &                ttime(N)
      ENDDO
      IF(NUMGBC.GT.0) THEN
        DO N=1,NUMGBC
          WRITE(1,1007) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
        ENDDO
      ENDIF
1001  FORMAT(I6,I4,1P2E12.5,0PF10.2,F7.2,F9.2,F10.3,F10.1,F10.5,2F10.5,
     &       I5,F10.3)
1002  FORMAT(5I6,1P2E17.10)
1007  FORMAT(2I6,E13.6)
1010  FORMAT(7X,F20.0)
      stop
 750  CONTINUE   
      WRITE(*,*) j,ntype(j)                                                       
      WRITE (*,5544) DETJ                                             
      WRITE (*,5544) xy                                             
 5544 FORMAT (' BAD JACOBIAN AT 750',E14.6,E14.6)                       
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
