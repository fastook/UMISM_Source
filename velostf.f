      PARAMETER(NMAX=29999)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER HED*80
      DIMENSION KODE(NMAX),X(NMAX),Y(NMAX),HTICE(NMAX),ADOT(NMAX),
     &          FRACT(NMAX),PSURF(NMAX),BDROCK(NMAX),FLOWA(NMAX),
     &          SLDGB(NMAX),
     &          KX(NMAX,4),CONST(NMAX),IBFLUX(NMAX,2),BFLUX(NMAX),
     &          TBED(NMAX),
     &          ACON(NMAX),NTYPE(NMAX),thick(nmax),
     &          BMELT(NMAX),WTHICK(NMAX)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2),LM(4)             
      DIMENSION PSI(4),DPSI(4,2)                                        
      DIMENSION XY(2,4)
      HMIN=10.
      icnt=0
      icount=0
      print *,'input display interval and how many to show'
      read(*,*) idisp,inumber
      print *,'input 1 for curvilinear, 0 for uniform'
      read(*,*) ibound
      if(ibound.eq.1) then
        write(*,*) 'surface display'
        WRITE(*,*) 'INPUT  1 FOR CALCULATED SURFACE,'
        WRITE(*,*) '       2 FOR MASS BALANCE,'
        WRITE(*,*) '       3 FOR FLOW CONSTANT,'
        WRITE(*,*) '       4 FOR DIFFERENCE,'
        WRITE(*,*) '       5 FOR THICKNESS,'
        WRITE(*,*) '       6 FOR PRESENT SURFACE,'
        READ(*,*) ISPLOT
        write(*,*) 'basal display'
        WRITE(*,*) '       1 FOR BEDROCK,'
        WRITE(*,*) '       2 FOR SLIDING FRACTION,'
        WRITE(*,*) '       3 FOR SLIDING CONSTANT,'
        WRITE(*,*) '       4 FOR BASAL TEMPERATURE,'
        WRITE(*,*) '       5 FOR BASAL MELT RATE,'
        WRITE(*,*) '       6 FOR WATER THICKNESS,'
        READ(*,*) IBPLOT

      endif
C READ INPUT HEADER
      READ(30,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      PRINT *,HED
1000  FORMAT (A80,/,7I6,F8.0)
      IF(NUMNP.GT.NMAX) THEN
        PRINT *,'NUMNP=',NUMNP,' NMAX=',NMAX,' INCREASE NMAX'
        STOP
      ENDIF
C READ INPUT GRID, THINGS THAT NEVER CHANGE
      READ(31) HED
      PRINT *,HED
      READ(31) (KODE(I),I=1,NUMNP)
      READ(31) (X(I),I=1,NUMNP)
      READ(31) (Y(I),I=1,NUMNP)
      READ(31) (PSURF(I),I=1,NUMNP)
      READ(31) (BDROCK(I),I=1,NUMNP)
      READ(31) (KX(I,1),KX(I,2),KX(I,3),KX(I,4),I=1,NUMEL)
      READ(31) (IBFLUX(I,1),IBFLUX(I,2),BFLUX(I),I=1,NUMGBC)
      DO NUM=1,NUMEL                                                      
        NTYPE(NUM)=1                                                    
        IF(KX(NUM,4).EQ.0) NTYPE(NUM)=2                                 
      ENDDO                                                             
C READ INPUT DIFFER, THINGS THAT DIFFER FROM ONE GRID TO THE NEXT
      READ(32) HED
      PRINT *,HED
      READ(32) (ADOT(I),I=1,NUMNP)
      READ(32) (FRACT(I),I=1,NUMNP)
      READ(32) (FLOWA(I),I=1,NUMNP)
      READ(32) (SLDGB(I),I=1,NUMNP)
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME
10    continue
      READ(33,END=999) HED
      PRINT *,HED
      READ(33) (HTICE(I),I=1,NUMNP)
      READ(33) (ADOT(I),I=1,NUMNP)
      READ(33) (BDROCK(I),I=1,NUMNP)
      READ(33) (CONST(I),I=1,NUMEL)
      READ(33) (ACON(I),I=1,NUMEL)
      READ(34) (FRACT(I),I=1,NUMNP)
      READ(34) (FLOWA(I),I=1,NUMNP)
      READ(34) (SLDGB(I),I=1,NUMNP)
      READ(34) (AFUDGE,I=1,NUMNP)
      DO NUM=1,NUMNP                                                    
        THICK(NUM)=HTICE(NUM)-BDROCK(NUM)                                     
      ENDDO                                                             
      READ(36,END=999) HED
      READ(36) (TBED(I),I=1,NUMNP)
      READ(36) (BMELT(I),I=1,NUMNP)
      READ(36) (WTHICK(I),I=1,NUMNP)
      icount=icount+1
      if(icount.gt.inumber*idisp) stop
      if(1+mod(icount,idisp).eq.1) then
        print *,'outputting for display'
c--------------------------------------      
      xmin=1e30
      xmax=-xmin
      ymin=1e30
      ymax=-ymin
      do i=1,numnp
        xmin=min(xmin,x(i))
        xmax=max(xmax,x(i))
        ymin=min(ymin,y(i))
        ymax=max(ymax,y(i))
      enddo
      if(ibound.eq.0) then
        facth=1.
        factx=1.e-5
        facty=1.e-5
      else
        facth=1./1000.
        factx=1./1000./100.
        facty=1./1000./100.
      endif
      write(11,*) 'RANK 2'
      write(11,*) 'DIMENSIONS ',NUMCOL,NUMLEV
      if(ibound.eq.0) then
        write(11,*) 'BOUNDS ',real(factx*xmin),real(factx*xmax)
     &                       ,real(facty*ymin),real(facty*ymax)
      endif
      write(11,*) 'NAME BEDROCK'
      write(11,*) 'TIME',icnt
c      write(11,*) 'TIME',hed(6:20)
      write(11,*) 'SCALAR'
      WRITE(11,*) 'DATA'
      if(ibplot.eq.1 .or. ibound.eq.0) then
        WRITE(11,*) (real(facth*BDROCK(I)),I=1,NUMNP)
      elseif(ibplot.eq.2) then
        WRITE(11,*) (real(FRACT(I)),I=1,NUMNP)
      elseif(ibplot.eq.3) then
        WRITE(11,*) (real(SLDGB(I)),I=1,NUMNP)
      elseif(ibplot.eq.4) then
        WRITE(11,*) (real(TBED(I)),I=1,NUMNP)
      elseif(ibplot.eq.5) then
        WRITE(11,*) (real(1000.*BMELT(I)),I=1,NUMNP)
      elseif(ibplot.eq.6) then
        WRITE(11,*) (real(WTHICK(I)),I=1,NUMNP)
      endif
      if(ibound.eq.1) then
        write(11,*) 'INTERLACED'
        write(11,*) 'VECTOR 3'
        WRITE(11,*) 'GRID'
        do i=1,numnp
          WRITE(11,*) factx*X(I),facty*Y(I),facth*bdrock(i)
c          WRITE(11,*) 1e-5*X(I),1e-5*Y(I),facth*bdrock(i)
        enddo
      endif
      WRITE(11,*) 'END'
      write(12,*) 'RANK 2'
      write(12,*) 'DIMENSIONS ',NUMCOL,NUMLEV
      if(ibound.eq.0) then
        write(12,*) 'BOUNDS ',real(factx*xmin),real(factx*xmax)
     &                       ,real(facty*ymin),real(facty*ymax)
      endif
      write(12,*) 'NAME SURF'
      write(12,*) 'TIME',icnt
      write(12,*) 'SCALAR'
      WRITE(12,*) 'DATA'
      if(isplot.eq.1 .or. ibound.eq.0) then
        WRITE(12,*) (real(facth*htice(I)),I=1,NUMNP)
      elseif(isplot.eq.2) then
        WRITE(12,*) (real(adot(I)),I=1,NUMNP)
      elseif(isplot.eq.3) then
        WRITE(12,*) (real(flowa(I)),I=1,NUMNP)
      elseif(isplot.eq.4) then
        WRITE(12,*) (real(facth*(htice(I)-psurf(i))),I=1,NUMNP)
      elseif(isplot.eq.5) then
        WRITE(12,*) (real(facth*(htice(I)-bdrock(i))),I=1,NUMNP)
      elseif(isplot.eq.6) then
        WRITE(12,*) (real(facth*psurf(I)),I=1,NUMNP)
      endif
      if(ibound.eq.1) then
        write(12,*) 'INTERLACED'
        write(12,*) 'VECTOR 3'
        WRITE(12,*) 'GRID'
        do i=1,numnp
          WRITE(12,*) factx*X(I),facty*Y(I),facth*htice(i)
c          WRITE(12,*) 1e-5*X(I),1e-5*Y(I),facth*htice(i)
        enddo
      endif
      write(13,*) 'RANK 2'
      write(13,*) 'DIMENSIONS ',NUMCOL-1,NUMLEV-1
      if(ibound.eq.0) then
        write(13,*) 'BOUNDS ',real(factx*xmin),real(factx*xmax)
     &                       ,real(facty*ymin),real(facty*ymax)
      endif
      write(13,*) 'NAME VELO'
      write(13,*) 'TIME',icnt
      write(13,*) 'INTERLACED'
      write(13,*) 'VECTOR 2'
      WRITE(13,*) 'DATA'
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
        XY(1,1)=X(I)                                                     
        XY(1,2)=X(JJ)                                                    
        XY(1,3)=X(K)                                                     
        IF(NTYPE(J).EQ.1) XY(1,4)=X(L)                                   
        XY(2,1)=Y(I)                                                     
        XY(2,2)=Y(JJ)                                                    
        XY(2,3)=Y(K)                                                     
        IF(NTYPE(J).EQ.1) XY(2,4)=Y(L)                                   
        IF(NTYPE(J).EQ.1) THEN                                            
          XCENT=(X(I)+X(JJ)+X(K)+X(L))/4.                          
          YCENT=(Y(I)+Y(JJ)+Y(K)+Y(L))/4.                          
        ELSE                                                              
          XCENT=(X(I)+X(JJ)+X(K))/3.                                
          YCENT=(Y(I)+Y(JJ)+Y(K))/3.                                
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
        IF (DETJ.LE.0.0) THEN                                         
          WRITE (*,5544) DETJ,X                                             
5544      FORMAT (' BAD JACOBIAN at 750',E14.6,E14.6)                       
          STOP
        ENDIF
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

        if(UMAG.gt.0. .and. UMAG.lt.1e5) then                                               
          WRITE(13,2000) real(ux*1),real(uy*1)
        else
          WRITE(13,*) ' Missing Missing '
        endif
600   CONTINUE                                                          
      if(ibound.eq.1) then
        write(13,*) 'INTERLACED'
        write(13,*) 'VECTOR 3'
        WRITE(13,*) 'GRID'
        DO J=1,NUMEL                                                  
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
          XY(1,1)=X(I)                                                     
          XY(1,2)=X(JJ)                                                    
          XY(1,3)=X(K)                                                     
          IF(NTYPE(J).EQ.1) XY(1,4)=X(L)                                   
          XY(2,1)=Y(I)                                                     
          XY(2,2)=Y(JJ)                                                    
          XY(2,3)=Y(K)                                                     
          IF(NTYPE(J).EQ.1) XY(2,4)=Y(L)                                   
          IF(NTYPE(J).EQ.1) THEN                                            
            XCENT=(X(I)+X(JJ)+X(K)+X(L))/4.                          
            YCENT=(Y(I)+Y(JJ)+Y(K)+Y(L))/4.                          
            HCENT=(HTICE(I)+HTICE(JJ)+HTICE(K)+HTICE(L))/4.                          
          ELSE                                                              
            XCENT=(X(I)+X(JJ)+X(K))/3.                                
            YCENT=(Y(I)+Y(JJ)+Y(K))/3.                                
            HCENT=(HTICE(I)+HTICE(JJ)+HTICE(K))/3.                                
          ENDIF                                                             
          WRITE(13,*) factx*XCENT,facty*YCENT,facth*HTICE(i)
c           WRITE(13,2000) real(1e-5*XCENT),real(1e-5*YCENT),
c     &                 real(facth*HTICE(i))
2000  format(3(1x,f13.6))
        enddo
      endif
      WRITE(12,*) 'END'
      WRITE(13,*) 'END'
      icnt=icnt+1
c------------------------
      endif
      goto 10
999   END
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
