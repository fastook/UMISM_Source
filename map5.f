      PARAMETER( MXX=29999, NZ=9 )
      PARAMETER( NAUX1=3*MXX+1000, NAUX2=MXX*(NZ-1)*1.5+(MXX+6) )       
C     IMPLICIT REAL*8(A-H,O-Z)                                          
C ***************************************************************C      
C                                                                C      
C   PROGRAM:  MAP5                                               C      
C                                                                C      
C   DATE:  8-24-90                                               C      
C   PROGRAMMER:  J.  FASTOOK                                     C      
C                                                                C      
C   FUNCTION:                                                    C      
C           THIS IS A PROGRAM TO MODEL THE FLOW OF A GLACIER     C      
C           WITH MATERIAL PROPERTIES READ FOR EACH NODAL POINT   C      
C           AND THE AVERAGE USED FOR THE ELEMENT.                C      
C           IT CALCULATES AN ELEMENT MATRIX ALA BECKER, 2-D      C      
C           PROGRAM WITH LINEAR SHAPE FUNCTIONS.                 C      
C           DOES TIME DEPENDENT CASE USING A LUMPED CAPACITANCE  C      
C           MATRIX AND A BACKWARD DIFFERENCE SCHEME. ALLOWING    C      
C           FOR VERY QUICK OPERATION, HOWEVER, MATERIAL PROPS    C      
C           ARE HELD CONSTANT THRUOUT THE TIME ITERATIONS        C      
C           SO THAT THE FINAL EQUILIBRIUM STATE DIFFERS SLIGHLTY C      
C           FROM ONE WHICH RECALCULATED THE STIFFNESS MAATRIX AT C      
C           EACH TIME STEP                                       C      
C           THIS IS HANDLED BY OUTPUTING A COMPLETE DATA SET AT  C      
C           THE END OF THE TIME STEP, (INDEED AT THE END OF      C      
C           EVERY PROGRAM RUN SO THAT BY TIME STEPPING ONE STEP  C      
C           AT A TIME AND USING THE NEW DATA SET IT IS AS IF THE C      
C           COMPLETE MATRIX WERE REDEVELOPED AT EACH STEP.       C      
C                                                                C      
C ***************************************************************C      
C     PROGRAM FOR STEADY AND UNSTEADY STATE FLOW ANALYSIS               
C     USING FINITE ELEMENTS                                             
C RUN BY MAP5 EXEC A:                                                   
C       FI * CLEAR                                                      
C       FI 30 DISK INPUT&1 HEAD B                                       
C       FI 31 DISK INPUT&1 GRID B                                       
C       FI 32 DISK INPUT&1 DIFF B                                       
C       FI 33 DISK INPUT&1 TIME B                                       
C       FI 34 DISK OUT&2 TIME B                                         
C       FI 10 DISK OUTLINE DATA B (LRECL 80 RECFM F                     
C       FI 7 DISK OUT1&2 DATA B (LRECL 130 RECFM F                      
C       FI 12 DISK OUT2&2 DATA B (LRECL 130 RECFM F                     
C       FI 11 DISK OUT3&2 DATA B (LRECL 130                             
C       FI 21 DISK OUT4&2 DATA B (LRECL 130                             
C       FI 20 DISK LINEM DATA B (LRECL 80                               
C       FI 26 DISK OUTPUT&2 DATA B (LRECL 130                           
C       FI 13 DISK MAT&2 DATA B                                         
C       FI 18 DISK VOL&2 DATA B                                         
C       FI 17 DISK VOLT&2 DATA B                                        
C       GL TXT  VSF2FORT CMSLIB IGLLIB ESSL                             
C       LOAD MAP5 (START NOMAP                                          
      CHARACTER*80 HED,SCRTCH                                           
      DIMENSION AMASS(9), B(MXX), X(MXX), Y(MXX), D(MXX),               
     &          KODE(MXX), CONST(MXX), LM(5), IX(3), E(3,3),            
     &          KX(MXX,4), P(5), S(5,5), DD(5), IDT(MXX),               
     &          ADOT(MXX), BDROCK(MXX), FLOWA(MXX), SLDGB(MXX),         
     &          PSURF(MXX), PPSURF(MXX), FRACT(MXX),                    
     &          CNEW(MXX), QHOLD(MXX), HTICE(MXX), THICK(MXX),          
     &          HFIT(MXX), IBFLUX(MXX,2), BFLUX(MXX), DEPB(MXX),        
     &          UNDEPB(MXX), SLOPE(MXX), SLOPN(MXX), KZ(MXX)            
      DIMENSION TTIME(2000),VVOL(2000)                                  
      DIMENSION NTYPE(MXX), NNODE(MXX),                                 
     &          DPSIX(9), DPSIY(9), DXDS(2,2), DSDX(2,2),               
     &          PSI(4), DPSI(4,2), CNST(MXX),                           
     &          XY(2,4), XI(2,9), ETA(2,9), W(2,9),                     
     &          AADOT(MXX), AFRACT(MXX), AFLOWA(MXX),                   
     &          ABDRCK(MXX), ASLDGB(MXX)                                
      REAL*8 A(MXX,NZ),Q(MXX),T(MXX),AUX1(NAUX1),AUX2(NAUX2)            
      REAL*8 RPARM(3),XGUESS(MXX)                                       
      INTEGER IPARM(4),KA(MXX,NZ)                                       
      COMMON /LAPSE/ ACOM,HMAX                                          
      DATA HFIT /MXX*0.0/                                               
      DATA IPARM /500,1,0,0/, RPARM /1.D-10,0.D0,0.D0/                  
      ACOM=-9.6237669                                                   
      DO 1 I=1,MXX                                                      
        B(I)=0.                                                         
        X(I)=0                                                          
        Y(I)=0.                                                         
        D(I)=0.                                                         
        KODE(I)=0                                                       
        CONST(I)=0.                                                     
        IDT(I)=0                                                        
        ADOT(I)=0.                                                      
        BDROCK(I)=0.                                                    
        FLOWA(I)=0.                                                     
        SLDGB(I)=0.                                                     
        PSURF(I)=0.                                                     
        PPSURF(I)=0.                                                    
        FRACT(I)=0.                                                     
        CNEW(I)=0.                                                      
        QHOLD(I)=0.                                                     
        HTICE(I)=0.                                                     
        THICK(I)=0.                                                     
        HFIT(I)=0.                                                      
        BFLUX(I)=0.                                                     
        DEPB(I)=0.                                                      
        UNDEPB(I)=0.                                                    
        SLOPE(I)=0.                                                     
        SLOPN(I)=0.                                                     
        KZ(I)=0                                                         
        NTYPE(I)=0                                                      
        NNODE(I)=0                                                      
        CNST(I)=0.                                                      
        AADOT(I)=0.                                                     
        AFRACT(I)=0.                                                    
        AFLOWA(I)=0.                                                    
        ABDRCK(I)=0.                                                    
        ASLDGB(I)=0.                                                    
1     CONTINUE                                                          
      CALL XUFLOW(0)                                                    
C FOLLOWING IS DEFAULT SNOWLINE ELEVATION AT POLE AND GRADIENT          
      AMASS(7)=-3000.                                                   
      AMASS(8)=.56E-3                                                   
      AMASS(8)=.001                                                     
      AMASS(9)=-14.                                                     
C **** FOLLOWING IS SEALEVEL REFERENCED TO PRESENT=0.                   
      NTSTEP=0                                                          
      SEALEV=0.                                                         
C **** FOLLOWING SETS RATE OF CONVERGENCE, UP TO 5 WORKS WELL           
      CONV=1.                                                           
      TIME=0.0                                                          
      RHOW=1.092                                                        
      PG = 0.089866*0.3816
      NUMCOL=0                                                          
      NUMLEV=0                                                          
      NCOL=MXX                                                          
C                                                                       
C     INITIALIZE INTEGRATION POINTS AND WEIGHTS                         
C     GAUSSIAN QUADRATURE OF ORDER THREE QUADRILATERALS                 
      CALL GAUSINIT(XI,ETA,W)                                           
C                                                                       
      DO 10 I=1,4                                                       
        LM(I)=0                                                         
10    CONTINUE                                                          
C                                                                       
C **** FOLLOWING READN FOR NEW SPLIT DATA SETS                          
      CALL READN(MXX, IDEP, HED, NUMNP, NUMEL, NUMGBC, NDT,             
     &           INTER, DT, KODE, X, Y, HTICE, ADOT, FRACT,             
     &           PSURF, RHOI, BDROCK, UNDEPB, FLOWA, SLDGB,             
     &           THICK, KX, CONST, IBFLUX, BFLUX, QHOLD,                
     &           NTYPE, NNODE, NCOL, DEPB, AADOT, AFRACT,               
     &           ABDRCK, PPSURF, AFLOWA, ASLDGB, IDT, AMASS)            
      DO 90 I=1,NUMNP                                                   
        Q(I)=HTICE(I)                                                   
        T(I)=HTICE(I)                                                   
90    CONTINUE                                                          
C                                                                       
C SET LINEARIZATION CONSTANT USING INITIAL CONFIGURATION                
      CALL NCONST(MXX, IDEP, X, Y, KX, NTYPE, NUMNP, NUMEL,             
     &            AFRACT, ASLDGB, LM, AFLOWA, BDROCK, DEPB,             
     &            UNDEPB, KODE, PG, Q, CNEW, SLOPE, RHOI)               
C                                                                       
      WRITE(*,*) 'TIME STEP=',DT                                        
      WRITE(*,*) 'INPUT 1 TO CALL ADJUST, 0 TO BYPASS'                  
      READ(*,*) IADJ                                                    
      IF(IADJ.EQ.1) THEN                                                
C                                                                       
C CALCULATE SLOPES IN CASE NEEDED BY ADJUST                             
        CALL NODESL(MXX, NUMNP, NUMEL, KX, SLOPE, SLOPN)                
C                                                                       
C ENTER INTERACTIVE DATA SET MANIPULATOR                                
        CALL ADJUST(HED, NUMNP, NUMEL, X, Y, HTICE, ADOT, FRACT,        
     &         PSURF, BDROCK, DEPB, FLOWA, SLDGB, THICK, KX, CONST,     
     &              NNODE, KODE, HFIT, NUMCOL, NUMLEV, NUMGBC, NDT,     
     &              INTER, DT, IBFLUX, BFLUX, MXX, IDT, SLOPN, AMASS,   
     &              TIME, NTSTEP,TTIME,VVOL)                            
C                                                                       
      ENDIF                                                             
      DO 91 I=1,NUMNP                                                   
        Q(I)=HTICE(I)                                                   
91    CONTINUE                                                          
C                                                                       
C LOAD NODAL PROPERTIES INTO ELEMENT PROPERTIES                         
      CALL ELPROP(MXX, NUMEL, NTYPE, KX, ADOT, AADOT, FRACT, AFRACT,    
     &            BDROCK, ABDRCK, PSURF, PPSURF, FLOWA, AFLOWA,         
     &            SLDGB, ASLDGB)                                        
C                                                                       
C                                                                       
C                                                                       
C                                                                       
C                                                                       
C                      *** MAIN LOOP ***                                
C                                                                       
C                                                                       
C                                                                       
215   CONTINUE                                                          
C                                                                       
C       FOLLOWING CALCULATES NODE SLOPE FROM ELEMENT SLOPE              
        CALL NODESL(MXX, NUMNP, NUMEL, KX, SLOPE, SLOPN)                
C                                                                       
C       FOLLOWING ADJUSTS ADOT FOR FITTEED ACCUMULATION                 
        DO 101 I=1,NUMNP                                                
          IF(IDT(I).GT.0) ADOT(I)=AFUNCT(TIME, IDT(I), AMASS,           
     &                                 REAL(Q(I)), BDROCK(I),           
     &                                 SLOPN(I),X(I),Y(I))              
101     CONTINUE                                                        
C                                                                       
C LOADS NEW NODAL MATERIAL PROPERTIES INTO ELEMENT MATERIAL PROPERTIES  
        CALL ELPROP(MXX, NUMEL, NTYPE, KX, ADOT, AADOT, FRACT, AFRACT,  
     &            BDROCK, ABDRCK, PSURF, PPSURF, FLOWA, AFLOWA,         
     &            SLDGB, ASLDGB)                                        
C                                                                       
C CALCULATES VOLUMES (FLOTATION AND TOTAL) AND AREAL EXTENT             
        CALL VOLUME(MXX, TIME, NUMNP, NUMEL, X, Y, KX, Q, BDROCK,       
     &             DEPB, SEALEV, RHOI, RHOW, VOL, AMASS)                
C                                                                       
      NTSTEP=NTSTEP+1                                                   
      TTIME(NTSTEP)=TIME                                                
      VVOL(NTSTEP)=VOL*1.E-15                                           
C                                                                       
C       FORM STIFFNESS,CAPACITANCE AND LOAD                             
        CALL FORMC(MXX, X, Y, KX, NTYPE, NUMNP, NUMEL, NCOL, ETA, XI,   
     &            W, CONST, ADOT, FRACT, BDROCK, PSURF, RHOI, FLOWA,    
     &           SLDGB, T, KODE, NUMGBC, IBFLUX, BFLUX, NZ, KZ, LM,     
     &           AADOT, AFRACT, ABDRCK, AFLOWA, ASLDGB, D, B, A, KA)    
C                                                                       
C       FORM EFFECTIVE CONDUCTIVITY MATRIX FOR TIME INCREMENT           
          DT2=1.0/DT                                                    
          DO 320 N=1,NUMNP                                              
            IF (KODE(N).EQ.0) THEN                                      
              IF (D(N).NE.0.) THEN                                      
                D(N)=DT2*D(N)                                           
                A(N,1)=A(N,1)+D(N)                                      
              ENDIF                                                     
            ENDIF                                                       
320       CONTINUE                                                      
          LL=0                                                          
C                                                                       
C                                                                       
C                                                                       
C         LOOP ON NUMBER OF TIME STEPS *******************************  
          DO 450 L=1,NDT                                                
C                                                                       
C           CALCULATE EFFECTIVE LOAD MATRIX                             
            DO 360 I=1,NUMNP                                            
              Q(I)=B(I)                                                 
              IF (KODE(I).EQ.0) Q(I)=B(I)+D(I)*T(I)                     
360         CONTINUE                                                    
C                                                                       
C FOR ITERATIVE MATRIX SOLVER USE AS INITIAL GUESS OF SOLUTION...       
            DO 1212 LK=1,NUMNP                                          
              XGUESS(LK)=T(LK)                                          
1212        CONTINUE                                                    
C                                                                       
C ESSL CONJUGATE GRADIENT MATRIX SOLVER                                 
C SOLUTION RETURNS IN XGUESS                                            
            CALL DSMCG(NUMNP, NZ, A, KA, MXX, Q, XGUESS, IPARM, RPARM,  
     &             AUX1, NAUX1, AUX2, NAUX2)                            
C                                                                       
            DO 5521 JK=1,NUMNP                                          
              Q(JK)=XGUESS(JK)                                          
5521        CONTINUE                                                    
            TIME=TIME+DT                                                
            WRITE(*,*) 'TIME=',TIME                                     
            LL=LL+1                                                     
            HMAX=-1.E30                                                 
            DIFF=0.                                                     
            NDIFF=0                                                     
            DO 5511 JK=1,NUMNP                                          
              IF(KODE(JK).NE.1) THEN                                    
                DIFF=DIFF+(Q(JK)-PSURF(JK))**2                          
                NDIFF=NDIFF+1                                           
                IF(Q(JK).GT.HMAX) THEN                                  
                  HMAX=Q(JK)                                            
                  NMAX=JK                                               
                ENDIF                                                   
                IF (Q(JK).LE.UNDEPB(JK)) Q(JK)=UNDEPB(JK)               
C                                                                       
C   ****        TO DEAL WITH BED BELOW SEA LEVEL AND                    
C   ****        SURFACE BELOW FLOTATION HEIGHT                          
                IF(UNDEPB(JK).LE.0.) THEN                               
                  FLOT=(1.-RHOW/RHOI)*UNDEPB(JK)                        
C                 IF(Q(JK).LE.FLOT) Q(JK)=FLOT                          
                  IF(Q(JK).LE.FLOT) Q(JK)=0.                            
                  IF(BDROCK(JK).LE.-9999.) Q(JK)=0.                     
                ENDIF                                                   
C   ****        END FLOTATION PART ****                                 
              ENDIF                                                     
5511        CONTINUE                                                    
C                                                                       
C THIS IS SPECIAL TO DEAL WITH TERRYS TEMP DEP STUFF, PLEASE REMOVE     
C       DO 5512 JK=1,NUMNP                                              
C         IF(Q(JK).GT.1250.) THEN                                       
C           FRACT(JK)=FRACT(JK)+.05                                     
C           IF(FRACT(JK).GT.1.) FRACT(JK)=1.                            
C         ELSE                                                          
C           FRACT(JK)=FRACT(JK)-.05                                     
C           IF(FRACT(JK).LT.0.) FRACT(JK)=0.                            
C         ENDIF                                                         
C 5512  CONTINUE                                                        
C       CALL ELPROP(MXX, NUMEL, NTYPE, KX, ADOT, AADOT, FRACT, AFRACT,  
C      &            BDROCK, ABDRCK, PSURF, PPSURF, FLOWA, AFLOWA,       
C      &            SLDGB, ASLDGB)                                      
C END TERRYS SPECIAL STUFF                                              
C                                                                       
C OUTPUT TIME STEP INFO TO SCREEN                                       
            WRITE(*,*) 'MXX SURF=',HMAX,' AT NODE',NMAX,' DIFF=',       
     &              SQRT(DIFF/REAL(NDIFF))                              
C                                                                       
C OBTAIN NEW LINEARIZATION CONSTANT FROM LATEST SOLUTION                
            CALL NCONST(MXX, IDEP, X, Y, KX, NTYPE, NUMNP, NUMEL,       
     &             AFRACT, ASLDGB, LM, AFLOWA, BDROCK, DEPB, UNDEPB,    
     &             KODE, PG, Q, CNEW, SLOPE, RHOI)                      
C                                                                       
C DERIVE SLOPE IN CASE NEEDED BY MASS BALANCE PARAMETERIZATION          
            CALL NODESL(MXX, NUMNP, NUMEL, KX, SLOPE, SLOPN)            
            DO 102 I=1,NUMNP                                            
              IF(IDT(I).GT.0) ADOT(I)=AFUNCT(TIME, IDT(I), AMASS,       
     &                            REAL(Q(I)), BDROCK(I), SLOPN(I),      
     &                            X(I),Y(I))                            
102         CONTINUE                                                    
C                                                                       
C LOAD NEW NODE PROPERTIES INTO ELEMENT PROPERTIES                      
            CALL ELPROP(MXX, NUMEL, NTYPE, KX, ADOT, AADOT, FRACT,      
     &             AFRACT, BDROCK, ABDRCK, PSURF, PPSURF, FLOWA,        
     &             AFLOWA, SLDGB, ASLDGB)                               
C                                                                       
C CALCULATE VOLUMES (FLOTATION AND TOTAL) AND AREA                      
            CALL VOLUME(MXX, TIME, NUMNP, NUMEL, X, Y, KX, Q, BDROCK,   
     &            DEPB, SEALEV, RHOI, RHOW, VOL, AMASS)                 
C                                                                       
      NTSTEP=NTSTEP+1                                                   
      TTIME(NTSTEP)=TIME                                                
      VVOL(NTSTEP)=VOL*1.E-15                                           
            DO 380 I=1,NUMNP                                            
              T(I)=Q(I)                                                 
380         CONTINUE                                                    
C                                                                       
C LOAD NEW LINEARIZATION CONSTANT INTO OLD                              
            DO 381 I=1,NUMEL                                            
C             CONST(I)=CNEW(I)                                          
              CONST(I) = (CONV*CONST(I) + CNEW(I))/(CONV+1)             
381         CONTINUE                                                    
C                                                                       
C           FORM STIFF,CAP,LOAD FOR NEXT PASS                           
            CALL FORMC(MXX, X, Y, KX, NTYPE, NUMNP, NUMEL, NCOL, ETA,   
     &            XI, W, CONST, ADOT, FRACT, BDROCK, PSURF, RHOI,       
     &            FLOWA, SLDGB, T, KODE, NUMGBC, IBFLUX, BFLUX, NZ,     
     &            KZ, LM, AADOT, AFRACT, ABDRCK, AFLOWA, ASLDGB, D, B,  
     &            A, KA)                                                
C                                                                       
            DO 3201 N=1,NUMNP                                           
              IF (KODE(N).EQ.0) THEN                                    
                IF (D(N).NE.0.0) THEN                                   
                  D(N)=DT2*D(N)                                         
                  A(N,1)=A(N,1)+D(N)                                    
                ENDIF                                                   
              ENDIF                                                     
3201        CONTINUE                                                    
C                                                                       
C PRINT OUT SOLUTION FOR APPROPRIATE TIME STEPS                         
            IF(LL.GE.INTER) THEN                                        
              DO 731 N=1,NUMNP                                          
                HTICE(N)=Q(N)                                           
731           CONTINUE                                                  
              WRITE(SCRTCH,*) 'TIME=',TIME                              
              WRITE(34) SCRTCH                                          
              WRITE(34) (HTICE(I),I=1,NUMNP)                            
              WRITE(34) (ADOT(I),I=1,NUMNP)                             
              WRITE(34) (DEPB(I),I=1,NUMNP)                             
              WRITE(34) (CONST(I),I=1,NUMEL)                            
              LL=0                                                      
            ENDIF                                                       
450       CONTINUE                                                      
C                                                                       
C     END OF TIME STEP SECTION ***************************************  
C                                                                       
C                                                                       
C                                                                       
      HMAX=-1.E30                                                       
      DIFF=0.                                                           
      NDIFF=0                                                           
      DO 730 N=1,NUMNP                                                  
        IF(KODE(N).NE.1) THEN                                           
          DIFF=DIFF+(Q(N)-PSURF(N))**2                                  
          NDIFF=NDIFF+1                                                 
        ENDIF                                                           
        HTICE(N)=Q(N)                                                   
        IF(HTICE(N).GT.HMAX) THEN                                       
          HMAX=HTICE(N)                                                 
          NMAX=N                                                        
        ENDIF                                                           
        IF (HTICE(N).LE.UNDEPB(N)) HTICE(N)=UNDEPB(N)                   
C ****  TO DEAL WITH BED BELOW SEA LEVEL                                
C ****  AND SURFACE BELOW FLOTATION HEIGHT                              
        IF(UNDEPB(N).LE.0.) THEN                                        
          FLOT=(1.-RHOW/RHOI)*UNDEPB(N)                                 
C         IF(HTICE(N).LE.FLOT) HTICE(N)=FLOT                            
          IF(HTICE(N).LE.FLOT) HTICE(N)=0.                              
          IF(BDROCK(N).LE.-9999.) HTICE(N)=0.                           
        ENDIF                                                           
C ****  END FLOTATION PART ****                                         
730   CONTINUE                                                          
C                                                                       
C OUTPUT TIME STEP INFO TO SCREEN                                       
      WRITE(*,*) 'MXX SURF=', HMAX, ' AT NODE',NMAX,' DIFF=',           
     &            SQRT(DIFF/REAL(NDIFF))                                
C                                                                       
C     FORMAT STATEMENTS                                                 
C                                                                       
1001  FORMAT(I6,I4,1P2E12.5,0PF10.2,F7.2,F9.2,F10.3,F10.1,F10.5,2F10.5,
     &       I5,F10.3)
1002  FORMAT(5I5,1PE17.10)                                              
1003  FORMAT(I10,4F10.0)                                                
1007  FORMAT(2I5,E13.6)                                                 
2002  FORMAT(2I2,1PE10.3,2E10.3,/,10X,6E10.3)                           
2004  FORMAT(5I3,7E14.4)                                                
2005  FORMAT(' TIME=', E12.5,/,(I6,E14.6,I6,E14.6,I6,E14.6,I6,E14.6,    
     & I6,E14.6,I6,E14.6))                                              
2007  FORMAT(2I5,2E15.6)                                                
2009  FORMAT('  M ',14X,'K',14X,'C',14X,'D',14X,'Q',/,(I4,4E15.6))      
2010  FORMAT(' AXISYMMETRIC SOLID BODY')                                
2011  FORMAT('TWO DIMENSIONAL PLANE BODY ')                             
2020  FORMAT(' CARD NO. ',I4, ' OUT OF ORDER')                          
2021  FORMAT(' BAD CARD NO. ',I4)                                       
3005  FORMAT(2F10.0,G13.6)                                              
C                                                                       
C                                                                       
      WRITE(*,8989)                                                     
8989  FORMAT (' END OF RUN')                                            
C     END OF RUN, ENTER ADJUST TO GO ROUND AGAIN.                       
      WRITE(*,*) 'INPUT 1 TO CALL ADJUST -9 TO SKIP AND END'            
      READ(*,*,END=999) IADJ                                            
      IF(IADJ.EQ.-9) GOTO 999                                           
C                                                                       
C ENTER INTERACTIVE DATA SET MANIPULATOR                                
      CALL ADJUST(HED, NUMNP, NUMEL, X, Y, HTICE, ADOT, FRACT, PSURF,   
     &        BDROCK, DEPB, FLOWA, SLDGB, THICK, KX, CONST, NNODE, KODE,
     &            HFIT, NUMCOL, NUMLEV, NUMGBC, NDT, INTER, DT,         
     &            IBFLUX, BFLUX, MXX, IDT, SLOPN, AMASS, TIME,          
     &            NTSTEP, TTIME, VVOL)                                  
C                                                                       
      DO 92 I=1,NUMNP                                                   
        Q(I)=HTICE(I)                                                   
92    CONTINUE                                                          
C                                                                       
C LOAD NODAL PROPERTIES INTO ELEMENT PROPERTIES                         
      CALL ELPROP(MXX, NUMEL, NTYPE, KX, ADOT, AADOT, FRACT, AFRACT,    
     &            BDROCK, ABDRCK, PSURF, PPSURF, FLOWA, AFLOWA, SLDGB,  
     &            ASLDGB)                                               
C                                                                       
      WRITE(*,*) 'INPUT 1 TO CONTINUE WITH NEW SET '                    
      READ(*,*,END=999) IADJ                                            
      IF(IADJ.EQ.-9) GOTO 999                                           
      ITER=1                                                            
C     GOTO BEGINNING OF MAIN LOOP                                       
C                                                                       
C                                                                       
C                                                                       
      GOTO 215                                                          
C                                                                       
C                                                                       
C                                                                       
999   CONTINUE                                                          
C                                                                       
C VERBOSE WRITER OFF                                                    
C     CALL WRITER(HED, NUMNP, NUMEL, X, Y, HTICE, ADOT, FRACT, PSURF,   
C    &            RHOI, BDROCK, FLOWA, SLDGB, THICK, KX, CONST, NNODE,  
C    &            KODE, HFIT, NUMCOL, NUMLEV, NUMGBC, NDT, INTER, DT,   
C    &            IBFLUX, BFLUX, AADOT, AFRACT, AFLOWA, ABDRCK,         
C    &            ASLDGB)                                               
C                                                                       
      WRITE(18,2000) -99999.,2.,0,HED                                   
2000  FORMAT(10X,G13.6,2X,G13.6,I13,/,A80)                              
      STOP                                                              
      END                                                               
C                                                                       
      SUBROUTINE SHAPE(NTYPE,XI,ET,PSI,DPSI)                            
C ELEMENT SHAPE FUNCTIONS AND DERIVATIVES AT LOCAL COORDINATES (XI,ET)  
C     IMPLICIT REAL*8(A-H,O-Z)                                          
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
C                                                                       
      SUBROUTINE NCONST(MXX,IDEP,X,Y,KX,NTYPE,NUMNP,NUMEL,AFRACT,ASLDGB,
     &    LM,AFLOWA,BDROCK,DEPB,UNDEPB,KODE,PG,Q,CNEW,SLOPE,RHOI)       
C CALCULATES LINEARIZATION CONSTANT FROM CURRENT SOLUTION               
      DIMENSION KODE(MXX),AFRACT(MXX),ASLDGB(MXX),AFLOWA(MXX)           
      DIMENSION BDROCK(MXX),DEPB(MXX),UNDEPB(MXX),SLOPE(MXX)            
      DIMENSION KX(MXX,4),X(MXX),Y(MXX),NTYPE(MXX)                      
      DIMENSION LM(5),CNEW(MXX)                                         
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2)                   
      DIMENSION PSI(4),DPSI(4,2)                                        
      DIMENSION XY(2,4)                                                 
      REAL*8 Q(MXX)                                                     
      RHOR=4.0                                                          
C                                                                       
C IDEP=1 DOES UNDEPRESSION AND DEPRESSION, IDEP=0 RIGID BED             
C                                                                       
      IF(IDEP.EQ.1) THEN                                                
        FDEP1=RHOR/(RHOR-RHOI)                                          
        FDEP2=RHOI/(RHOI-RHOR)                                          
      ELSE                                                              
        FDEP1=1.                                                        
        FDEP2=0.                                                        
      ENDIF                                                             
      DO 600 J = 1,NUMEL                                                
        IF(NTYPE(J).EQ.1) THEN                                          
          NODEN=4                                                       
          CENTX=0.0                                                     
          CENTY=0.0                                                     
        ELSE                                                            
          NODEN=3                                                       
          CENTX=1.D0/3.D0                                               
          CENTY=1.D0/3.D0                                               
        ENDIF                                                           
        SUMHH=0.                                                        
        SUMX=0.                                                         
        SUMY=0.                                                         
        DO 560 I = 1,NODEN                                              
          LM(I) = KX(J,I)                                               
560     CONTINUE                                                        
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
        CALL SHAPE(NTYPE(J),CENTX,CENTY,PSI,DPSI)                       
C                                                                       
C CALCULATE DXDS...EQUATION (5.3.6)                                     
C                                                                       
        DO 565 I=1,2                                                    
        DO 565 L=1,2                                                    
          DXDS(I,L)=0.0                                                 
          DO 565 K=1,NODEN                                              
            DXDS(I,L)=DXDS(I,L)+DPSI(K,L)*XY(I,K)                       
565     CONTINUE                                                        
C                                                                       
C CALCULATE DSDX...EQUATION (5.2.7)                                     
C                                                                       
        DETJ=(DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1))                  
        IF (DETJ.LE.0.0) THEN                                           
          WRITE(12,5544) J,DETJ                                         
          WRITE(*,5544) J,DETJ                                          
          WRITE(12,5545) (JJ,XY(1,JJ),XY(2,JJ),JJ=1,4)                  
          WRITE(*,5545) (JJ,XY(1,JJ),XY(2,JJ),JJ=1,4)                   
5545      FORMAT(1X,I5,1X,1PE10.3,E10.3)                                
          STOP                                                          
5544      FORMAT(' BAD JACOBIAN',I5,1PE10.3,/,1X,8E10.3)                
        ENDIF                                                           
        DSDX(1,1)=DXDS(2,2)/DETJ                                        
        DSDX(2,2)=DXDS(1,1)/DETJ                                        
        DSDX(1,2)=-DXDS(1,2)/DETJ                                       
        DSDX(2,1)=-DXDS(2,1)/DETJ                                       
C                                                                       
C CALCULATE D(PSI)/DX...EQUATION (5.3.5)                                
C                                                                       
        DO 570 I=1,NODEN                                                
          DPSIX(I)=DPSI(I,1)*DSDX(1,1)+DPSI(I,2)*DSDX(2,1)              
          DPSIY(I)=DPSI(I,1)*DSDX(1,2)+DPSI(I,2)*DSDX(2,2)              
570     CONTINUE                                                        
        DO 580 I = 1,NODEN                                              
          SUMX = SUMX + Q(LM(I))*DPSIX(I)                               
          SUMY = SUMY + Q(LM(I))*DPSIY(I)                               
          IF(BDROCK(LM(I)).LE.-9999.) THEN                              
            DEPB(LM(I))=0.                                              
          ELSE                                                          
            DEPB(LM(I))=FDEP2*Q(LM(I))+FDEP1*UNDEPB(LM(I))              
          ENDIF                                                         
          THIK=(Q(LM(I))-UNDEPB(LM(I)))*FDEP1                           
          IF(THIK.GT.0.) SUMHH=SUMHH+THIK                               
580     CONTINUE                                                        
C                                                                       
        DELH = SUMX**2 + SUMY**2                                        
        DELH = SQRT(DELH)                                               
        SLOPE(J)=DELH                                                   
        HH = SUMHH/FLOAT(NODEN)                                         
        TERM1 = AFRACT(J)*((PG/ASLDGB(J))**2)*(HH**3)*DELH              
        TERM2=(1.-AFRACT(J))*(.2)*((PG/AFLOWA(J))**3)*(HH**5)*(DELH**2) 
        CNEW(J) = TERM1 + TERM2                                         
  600 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE ELPROP(MXX,NUMEL,NTYPE,KX,ADOT,AADOT,FRACT,AFRACT,     
     &BDROCK,ABDRCK,PSURF,PPSURF,FLOWA,AFLOWA,SLDGB,ASLDGB)             
C LOAD NODAL PROPERTIES INTO ELEMENT PROPERTIES                         
      DIMENSION NTYPE(MXX),ADOT(MXX),AADOT(MXX),FRACT(MXX),AFRACT(MXX), 
     &BDROCK(MXX),ABDRCK(MXX),PSURF(MXX),PPSURF(MXX),FLOWA(MXX),        
     &AFLOWA(MXX),SLDGB(MXX),ASLDGB(MXX),LM(4),KX(MXX,4)                
      DO 100 N = 1,NUMEL                                                
        IF(NTYPE(N).EQ.1) THEN                                          
          NODEN=4                                                       
          NINT=9                                                        
        ELSE                                                            
          NODEN=3                                                       
          NINT=4                                                        
        ENDIF                                                           
        AAADOT=0.                                                       
        AAFRCT=0.                                                       
        AABDRK=0.                                                       
        APSURF=0.                                                       
        AAFLOW=0.                                                       
        AASLDG=0.                                                       
        ADNSTY=0.                                                       
        DO 90 I=1,4                                                     
          LM(I)=KX(N,I)                                                 
          AAADOT=AAADOT+ADOT(LM(I))                                     
          AAFRCT=AAFRCT+FRACT(LM(I))                                    
          AABDRK=AABDRK+BDROCK(LM(I))                                   
          APSURF=APSURF+PSURF(LM(I))                                    
          AAFLOW=AAFLOW+FLOWA(LM(I))                                    
          AASLDG=AASLDG+SLDGB(LM(I))                                    
90      CONTINUE                                                        
        DENOM=1./FLOAT(NODEN)                                           
        AAADOT=AAADOT*DENOM                                             
        AADOT(N)=AAADOT                                                 
        AAFRCT=AAFRCT*DENOM                                             
        AFRACT(N)=AAFRCT                                                
        AABDRK=AABDRK*DENOM                                             
        ABDRCK(N)=AABDRK                                                
        ADNSTY=ADNSTY*DENOM                                             
        PPSURF(N)=APSURF                                                
        AAFLOW=AAFLOW*DENOM                                             
        AFLOWA(N)=AAFLOW                                                
        AASLDG=AASLDG*DENOM                                             
        ASLDGB(N)=AASLDG                                                
100   CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE FORMC(MXX,X,Y,KX,NTYPE,NUMNP,NUMEL,NCOL,ETA,XI,W,      
     &     CONST,ADOT,FRACT,BDROCK,PSURF,RHOI,FLOWA,SLDGB,              
     &     T,KODE,NUMGBC,IBFLUX,BFLUX,                                  
     &        NZ,KZ,LM,AADOT,AFRACT,ABDRCK,AFLOWA,ASLDGB,D,B,A,KA)      
C FORM STIFFNESS MATRIX                                                 
      DIMENSION AADOT(MXX),AFRACT(MXX),AFLOWA(MXX)                      
      DIMENSION ABDRCK(MXX),ASLDGB(MXX),KODE(MXX),IBFLUX(MXX,2)         
      DIMENSION BFLUX(MXX),KZ(MXX)                                      
      DIMENSION ADOT(MXX),FRACT(MXX),PSURF(MXX),FLOWA(MXX)              
      DIMENSION BDROCK(MXX),SLDGB(MXX),CONST(MXX),LM(5)                 
      DIMENSION KX(MXX,4),X(MXX),Y(MXX),NTYPE(MXX),D(MXX),B(MXX)        
      REAL*8 A(MXX,NZ),T(MXX)                                           
      INTEGER KA(MXX,NZ)                                                
      DIMENSION P(5),S(5,5),DD(5)                                       
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2)                   
      DIMENSION PSI(4),DPSI(4,2)                                        
      DIMENSION XY(2,4),XI(2,9),ETA(2,9),W(2,9)                         
      BIG=1.E20                                                         
C **********************************************************************
C     FORM CONDUCTIVITY MATRIX FOR COMPLETE BODY                        
C **********************************************************************
      DO 110 I=1,NUMNP                                                  
        KZ(I)=1                                                         
        KA(I,1)=I                                                       
        D(I)=0.0                                                        
        B(I)=0.0                                                        
        DO 110 J=1,NZ                                                   
          KA(I,J)=0                                                     
          A(I,J)=0.0                                                    
110   CONTINUE                                                          
      DO 160 N=1,NUMEL                                                  
        IF(NTYPE(N).EQ.1) THEN                                          
          NODEN=4                                                       
          NINT=9                                                        
        ELSE                                                            
          NODEN=3                                                       
          NINT=4                                                        
        ENDIF                                                           
        DO 125 I=1,4                                                    
          LM(I)=KX(N,I)                                                 
125     CONTINUE                                                        
C                                                                       
C                                                                       
C     2. FORM ELEMENT CONDUCTIVITY MATRIX                               
C                                                                       
        DO 130 I=1,4                                                    
          DD(I)=0.0                                                     
          P(I)=0.0                                                      
          DO 130 J=1,4                                                  
            S(I,J)=0.0                                                  
130     CONTINUE                                                        
C                                                                       
        I=LM(1)                                                         
        J=LM(2)                                                         
        K=LM(3)                                                         
        L=LM(4)                                                         
        XY(1,1)=X(I)                                                    
        XY(1,2)=X(J)                                                    
        XY(1,3)=X(K)                                                    
C       IF(NTYPE(N).EQ.1) XY(1,4)=X(L)                                  
        XY(2,1)=Y(I)                                                    
        XY(2,2)=Y(J)                                                    
        XY(2,3)=Y(K)                                                    
C       IF(NTYPE(N).EQ.1) XY(2,4)=Y(L)                                  
        IF(NODEN.EQ.4) THEN                                             
          XY(1,4)=X(L)                                                  
          XY(2,4)=Y(L)                                                  
        ENDIF                                                           
C                                                                       
C FORM ELEMENT MATRIX AND VECTORS                                       
C                                                                       
C BEGIN INTEGRATION POINT LOOP                                          
        DO 139 L=1,NINT                                                 
          CALL SHAPE(NTYPE(N),XI(NTYPE(N),L),ETA(NTYPE(N),L),PSI,DPSI)  
C CALCULATE DXDS...EQUATION (5.3.6)                                     
          DO 133 I=1,2                                                  
          DO 133 J=1,2                                                  
            DXDS(I,J)=0.0                                               
            DO 133 K=1,4                                                
              DXDS(I,J)=DXDS(I,J)+DPSI(K,J)*XY(I,K)                     
133       CONTINUE                                                      
C CALCULATE DSDX...EQUATION (5.2.7)                                     
          DETJ=(DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1))                
          IF(DETJ.LE.0.0) THEN                                          
            WRITE(12,1100) DETJ,((XY(MM,NN),NN=1,4),MM=1,2)             
1100        FORMAT(' BAD JACOBIAN AT 161',G13.6,/,4G13.6,/,4G13.6)      
            STOP                                                        
          ENDIF                                                         
          DENOM=1./DETJ                                                 
          DSDX(1,1)=DXDS(2,2)*DENOM                                     
          DSDX(2,2)=DXDS(1,1)*DENOM                                     
          DSDX(1,2)=-DXDS(1,2)*DENOM                                    
          DSDX(2,1)=-DXDS(2,1)*DENOM                                    
C CALCULATE D(PSI)/DX...EQUATION (5.3.5)                                
          DO 136 I=1,4                                                  
            DPSIX(I)=DPSI(I,1)*DSDX(1,1)+DPSI(I,2)*DSDX(2,1)            
            DPSIY(I)=DPSI(I,1)*DSDX(1,2)+DPSI(I,2)*DSDX(2,2)            
136       CONTINUE                                                      
C ACCUMULATE INTEGRATION POINT VALUES OF INTEGRALS                      
          FAC=DETJ*W(NTYPE(N),L)                                        
          DO 138 I=1,4                                                  
            DD(I)=DD(I)+PSI(I)*FAC                                      
C THIS IS LUMPED CAPACITANCE MATRIX                                     
            P(I)=P(I)+AADOT(N)*PSI(I)*FAC                               
            DO 138 J=1,4                                                
C                                                                       
              IF(CONST(N).GT.1.E-30) THEN                               
                TERM1=CONST(N)*(DPSIX(I)*DPSIX(J)+DPSIY(I)*DPSIY(J))    
                S(I,J)=S(I,J)+TERM1*FAC                                 
              ENDIF                                                     
138       CONTINUE                                                      
139     CONTINUE                                                        
C                                                                       
C     3. ADD ELEMENT CONDUCTIVITY TO COMPLETE CONDUCTIVITY MATRIX       
C                                                                       
        DO 155 L=1,4                                                    
          I=LM(L)                                                       
          D(I)=D(I)+DD(L)                                               
C THIS (D) IS LUMPED CAPACITANCE MATRIX                                 
          B(I)=B(I)+P(L)                                                
          DO 155 M=1,4                                                  
            J=LM(M)                                                     
            IF(I.EQ.J) THEN                                             
              A(I,1)=A(I,1)+S(L,M)                                      
              KA(I,1)=I                                                 
            ELSE                                                        
              DO 154 K=2,KZ(I)                                          
                IF(KA(I,K).EQ.J) THEN                                   
                  A(I,K)=A(I,K)+S(L,M)                                  
                  GOTO 155                                              
                ENDIF                                                   
154           CONTINUE                                                  
            KZ(I)=KZ(I)+1                                               
            A(I,KZ(I))=S(L,M)                                           
            KA(I,KZ(I))=J                                               
          ENDIF                                                         
155     CONTINUE                                                        
160   CONTINUE                                                          
C                                                                       
C                                                                       
C     BOUNDARY CONDITIONS                                               
C                                                                       
C                                                                       
      IF(NUMGBC.GT.0) THEN                                              
        DO 2101 N=1,NUMGBC                                              
          I=IBFLUX(N,1)                                                 
          J=IBFLUX(N,2)                                                 
          TEMP=BFLUX(N)                                                 
          XL=SQRT((X(J)-X(I))**2+(Y(J)-Y(I))**2)                        
          TEMP=XL*TEMP*.5                                               
          B(I)=B(I)+TEMP                                                
          B(J)=B(J)+TEMP                                                
2101    CONTINUE                                                        
      ENDIF                                                             
C                                                                       
C     2. FIXED BOUNDARY CONDITIONS                                      
C            BY PENALTY METHOD                                          
C                                                                       
      DO 300 N=1,NUMNP                                                  
        IF(KODE(N).EQ.1) THEN                                           
          A(N,1)=BIG                                                    
          B(N)=T(N)*BIG                                                 
        ENDIF                                                           
 300  CONTINUE                                                          
      END                                                               
C                                                                       
      SUBROUTINE VOLUME(MXX,TIME,NUMNP,NUMEL,X,Y,KX,Q,BDROCK,DEPB,      
     &SEALEV,RHOI,RHOW,VOL,AMASS)                                       
C CALCULATES VOLUMES (FLOTATION AND TOTAL) AND AREA                     
      DIMENSION BDROCK(MXX),KX(MXX,4),X(MXX),Y(MXX),LM(4)               
      DIMENSION DEPB(MXX)                                               
      REAL*8 Q(MXX)                                                     
      DIMENSION AMASS(9)                                                
      COMMON /LAPSE/ ACOM,HMAX                                          
      VOL=0.0                                                           
      VOL1=0.                                                           
      AREATOT=0.0                                                       
      DO 100 I=1,NUMEL                                                  
        SUMH=0.                                                         
        SUMT=0.                                                         
        DO 90 J=1,4                                                     
        LM(J)=KX(I,J)                                                   
        IF(BDROCK(LM(J)).LE.-9999.) GOTO 100                            
        IF(DEPB(LM(J)).LT.SEALEV) THEN                                  
          FLOT=(1.-RHOW/RHOI)*(DEPB(LM(J))-SEALEV)                      
        ELSE                                                            
          FLOT=DEPB(LM(J))                                              
        ENDIF                                                           
        HEIGHT=Q(LM(J))-FLOT                                            
        THICK=Q(LM(J))-DEPB(LM(J))                                      
        IF(HEIGHT.GT.1.) SUMH=SUMH+HEIGHT                               
        IF(THICK.GT.1.) SUMT=SUMT+THICK                                 
90    CONTINUE                                                          
      HEIGHT=SUMH*.25                                                   
      THICK=SUMT*.25                                                    
      AREA=0.5*((X(LM(2))-X(LM(1)))*(Y(LM(3))-Y(LM(2)))-                
     &          (X(LM(3))-X(LM(2)))*(Y(LM(2))-Y(LM(1)))+                
     &          (X(LM(4))-X(LM(3)))*(Y(LM(1))-Y(LM(4)))-                
     &          (X(LM(1))-X(LM(4)))*(Y(LM(4))-Y(LM(3))))                
C IF THICKNESS LT 1 METER, NEGLECT                                      
      IF(HEIGHT.GT.1.) THEN                                             
        VOL=VOL+AREA*HEIGHT                                             
        VOL1=VOL1+AREA*THICK                                            
        AREATOT=AREATOT+AREA                                            
      ENDIF                                                             
100   CONTINUE                                                          
      IF(AREATOT.GT.0.) THEN                                            
        AVGHGT=VOL/AREATOT                                              
      ELSE                                                              
        AVGHGT=0.                                                       
      ENDIF                                                             
      WRITE(*,1000) VOL*1.E-15,AREATOT*1.E-12,AVGHGT                    
      WRITE(*,1010) VOL1*1.E-15,AMASS(7),AMASS(8)                       
C     WRITE(*,1010) VOL1*1.E-15,ACOM,AMASS(9)                           
C     WRITE(18,2000) TIME,VOL*1.E-15                                    
      WRITE(18,2000) TIME,VOL1*1.E-15                                   
2000  FORMAT(10X,G13.6,2X,G13.6)                                        
1000  FORMAT(' V=',G13.6,' MEG KM**3, A=',G13.6,' MEG KM**2, <H>=',G13.6
     &)                                                                 
1010  FORMAT(' V=',G13.6,' MEG KM**3, L=',G13.6,' DEG / M  ,TSNL=',G13.6
     &)                                                                 
C     WRITE(17,1004) TIME,VOL1*1E-15,VOL*1E-15,AREATOT*1E-12,AMASS(7),  
      WRITE(17,1004) TIME,VOL*1E-15,VOL1*1E-15,AREATOT*1E-12,AMASS(7),  
     &              AMASS(9),AVGHGT,HMAX                                
1004  FORMAT(1X,1P,G13.6,0P,3F7.3,F7.0,F6.1,F7.1,F7.1)                  
1001  FORMAT(' TIME=',G13.6)                                            
1002  FORMAT(' TOTAL VOLUME=',G13.6,' MEG KM**3')                       
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE NODESL(MXX,NUMNP,NUMEL,KX,SLOPE,SLOPN)                 
C CALCULATES SLOES FOR ACCUMULATION PARAMETERIZATIONS                   
      PARAMETER(MMXX=29999)                                              
      DIMENSION KX(MXX,4),SLOPE(MXX),SLOPN(MXX),ICOUNT(MMXX)            
      DO 100 I=1,NUMNP                                                  
        SLOPN(I)=0.                                                     
        ICOUNT(I)=0                                                     
100   CONTINUE                                                          
      DO 200 I=1,NUMEL                                                  
        DO 200 J=1,4                                                    
          SLOPN(KX(I,J))=SLOPN(KX(I,J))+SLOPE(I)                        
          ICOUNT(KX(I,J))=ICOUNT(KX(I,J))+1                             
200   CONTINUE                                                          
      DO 300 I=1,NUMNP                                                  
        IF(ICOUNT(I).GT.0) THEN                                         
          SLOPN(I)=SLOPN(I)/REAL(ICOUNT(I))                             
        ENDIF                                                           
300   CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE GAUSINIT(XI,ETA,W)                                     
      DIMENSION XI(2,9), ETA(2,9), W(2,9)                               
C     GAUSSIAN QUADRATURE OF ORDER THREE QUADRILATERALS                 
      XI(1,1)=-SQRT(3./5.)                                              
      XI(1,2)=0.                                                        
      XI(1,3)=-XI(1,1)                                                  
      XI(1,4)=XI(1,1)                                                   
      XI(1,5)=0.                                                        
      XI(1,6)=XI(1,3)                                                   
      XI(1,7)=XI(1,1)                                                   
      XI(1,8)=0.                                                        
      XI(1,9)=XI(1,3)                                                   
      ETA(1,1)=XI(1,1)                                                  
      ETA(1,2)=XI(1,1)                                                  
      ETA(1,3)=XI(1,1)                                                  
      ETA(1,4)=0.                                                       
      ETA(1,5)=0.                                                       
      ETA(1,6)=0.                                                       
      ETA(1,7)=XI(1,3)                                                  
      ETA(1,8)=XI(1,3)                                                  
      ETA(1,9)=XI(1,3)                                                  
      W(1,1)=25./81.                                                    
      W(1,2)=40./81.                                                    
      W(1,3)=W(1,1)                                                     
      W(1,4)=W(1,2)                                                     
      W(1,5)=64./81.                                                    
      W(1,6)=W(1,2)                                                     
      W(1,7)=W(1,1)                                                     
      W(1,8)=W(1,2)                                                     
      W(1,9)=W(1,1)                                                     
C   GAUSSIAN QUADRATURE OF ORDER THREE TRIANGLES                        
      XI(2,1)=1./3.                                                     
      XI(2,2)=2./15.                                                    
      XI(2,3)=XI(2,2)                                                   
      XI(2,4)=11./15.                                                   
      ETA(2,1)=XI(2,1)                                                  
      ETA(2,2)=XI(2,4)                                                  
      ETA(2,3)=XI(2,2)                                                  
      ETA(2,4)=ETA(2,3)                                                 
      W(2,1)=-27./96.                                                   
      W(2,2)=25./96.                                                    
      W(2,3)=W(2,2)                                                     
      W(2,4)=W(2,2)                                                     
      END                                                               
