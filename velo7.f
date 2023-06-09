C **************************************************************C       
C                                                               C       
C   PROGRAM:  VELO7                                             C       
C                                                               C       
C   DATE:  04 14 92                                             C       
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
      PARAMETER(MXX=29999,HMIN=10.,NPLOT=10000)                           
      CHARACTER HED*80,junk*80                                                  
      LOGICAL FOUND
      REAL*4  PX(4),PY(4)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2),LM(4)             
      DIMENSION PSI(4),DPSI(4,2)                                        
      DIMENSION XY(2,4),ZZ(MXX),ZZZ(MXX)
      DIMENSION XX(MXX),YY(MXX),HTICE(MXX),THICK(MXX),NTYPE(MXX)        
      DIMENSION KX(MXX,4),CONST(MXX)
      DIMENSION PSURF(MXX),FRACT(MXX),FLOWA(MXX),SLDGB(MXX)             
      DIMENSION ADOT(MXX),DEPB(MXX),KODE(MXX),AJUNK(MXX)
      DIMENSION UXELEM(MXX),UYELEM(MXX),UMAGE(MXX),UMAGN(MXX)
      DIMENSION JCOUNT(MXX),IK(5)
      DIMENSION ICMAP1(16)
      DATA ICMAP1 /0,4,12,6,13,2,8,7,9,3,10,5,11,1,14,15/                
      DIMENSION ICMAP2(3)
      DATA ICMAP2 /1,14,15/
C     DIMENSION ICMAP1(11)
C     DATA ICMAP1 /11,10,9,8,7,6,5,4,3,2,1/
C     DATA ICMAP1 /0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/                
      NCOLOR=12
      RLOGP1=LOG(.1)
      print *,'input time you want to extract'
      read(*,*) jtime
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
      DO NUM=1,NUMNP
        ZZZ(NUM)=PSURF(NUM)
      ENDDO
      READ(31) (IB1,IB2,BJUNK,I=1,NUMGBC)                               
C READ INPUT DIFFER, THINGS THAT DIFFER FROM ONE GRID TO THE NEXT       
      READ(32) HED                                                      
      READ(32) (AJUNK(I),I=1,NUMNP)                                        
      READ(32) (FRACT(I),I=1,NUMNP)                                     
      READ(32) (FLOWA(I),I=1,NUMNP)                                     
      READ(32) (SLDGB(I),I=1,NUMNP)                                     
      XMAX=-1.E30                                                       
      XMIN=1.E30                                                        
      YMAX=-1.E30                                                       
      YMIN=1.E30                                                        
      DO I=1,NUMNP                                                
        XMAX=MAX(XMAX,XX(I))
        XMIN=MIN(XMIN,XX(I))
        YMAX=MAX(YMAX,YY(I))
        YMIN=MIN(YMIN,YY(I))
      ENDDO                                                        
      XWID=XMAX-XMIN
      YWID=YMAX-YMIN
      IF(XWID.GT.YWID) THEN
        YWID=XWID
        YMAX=YMIN+YWID
      ELSE
        XWID=YWID
        XMAX=XMIN+XWID
      ENDIF
      XMIN=XMIN-XWID*.05
      XMAX=XMAX+XWID*.25
      YMIN=YMIN-YWID*.05
      YMAX=YMAX+YWID*.25
      VMIN=RLOGP1
      VMAX=LOG(3000.)
      VDELT=(VMAX-VMIN)/(NCOLOR-1)
      CALL GRSTRT(800,1)                                                
      CALL WINDOW(REAL(XMIN),REAL(XMAX),REAL(YMIN),REAL(YMAX))
      XPOS1=XMAX-XWID*.25
      XPOS2=XMAX-XWID*.15
      YPOS=YMIN+YWID*.15
      YDELT=YWID/(NCOLOR-1)
      CALL OPNSEG(1)
      DO I=1,NCOLOR
        CALL LINCLR(ICMAP1(I))
        CALL MOVE(REAL(XPOS1),REAL(YPOS))
        CALL DRAW(REAL(XPOS2),REAL(YPOS))
        RNUM=VMIN+(I-1)*VDELT
        RNUM=EXP(RNUM)
        CALL LINCLR(1)
        CALL GNUMBR(REAL(RNUM),2,6)
        YPOS=YPOS+YDELT
      ENDDO
      CALL LINCLR(1)                                                    
      CALL MOVE(REAL(XMIN),REAL(YMIN))                                  
      CALL DRAW(REAL(XMIN),REAL(YMAX))                                  
      CALL DRAW(REAL(XMAX),REAL(YMAX))                                  
      CALL DRAW(REAL(XMAX),REAL(YMIN))                                  
      CALL DRAW(REAL(XMIN),REAL(YMIN))                                  
      CALL LINCLR(1)                                                    
      CALL CLOSEG
      ISEG=1
1     CONTINUE                                                          
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME                         
      READ(33,END=99) HED                                              
      write(junk,*) hed
      read(junk,2000) ttime
2000  format(5x,f20.0)
      itime=nint(ttime)
      print *,itime
      READ(33) (HTICE(I),I=1,NUMNP)                                     
      READ(33) (ADOT(I),I=1,NUMNP)                                      
      READ(33) (DEPB(I),I=1,NUMNP)                                      
      READ(33) (CONST(I),I=1,NUMEL)                                     
      READ(33) (ACON,I=1,NUMEL)                                     
      ISEG=ISEG+1
      CALL OPNSEG(ISEG)
      DO NUM=1,NUMNP                                                    
        THICK(NUM)=HTICE(NUM)-DEPB(NUM)                                     
        ZZ(NUM)=HTICE(NUM)
      ENDDO                                                             
23    CONTINUE                                                          
      DO I=1,NUMNP
        JCOUNT(I)=0
        UMAGN(I)=0.
      ENDDO
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
        II=LM(1)                                                           
        JJ=LM(2)                                                          
        KK=LM(3)                                                           
        LL=LM(4)                                                           
        XY(1,1)=XX(II)                                                     
        XY(1,2)=XX(JJ)                                                    
        XY(1,3)=XX(KK)                                                     
        IF(NTYPE(J).EQ.1) XY(1,4)=XX(LL)                                   
        XY(2,1)=YY(II)                                                     
        XY(2,2)=YY(JJ)                                                    
        XY(2,3)=YY(KK)                                                     
        IF(NTYPE(J).EQ.1) XY(2,4)=YY(LL)                                   
        IF(NTYPE(J).EQ.1) THEN                                            
          XCENT=(XX(II)+XX(JJ)+XX(KK)+XX(LL))/4000.                          
          YCENT=(YY(II)+YY(JJ)+YY(KK)+YY(LL))/4000.                          
        ELSE                                                              
          XCENT=(XX(II)+XX(JJ)+XX(KK))/3000.                                
          YCENT=(YY(II)+YY(JJ)+YY(KK))/3000.                                
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
          WRITE(*,*) 'DET WRONG, ELEMENT ',J
          STOP
        ENDIF
C       IF (DETJ.GE.0.0) GOTO 750                                         
        DENOM=1./DETJ
        DSDX(1,1)=DXDS(2,2)*DENOM
        DSDX(2,2)=DXDS(1,1)*DENOM
        DSDX(1,2)=-DXDS(1,2)*DENOM
        DSDX(2,1)=-DXDS(2,1)*DENOM
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
          UXELEM(J)=-CONST(J)*SUMX/HH 
          UYELEM(J)=-CONST(J)*SUMY/HH
        ELSE                                                              
          UXELEM(J)=0.
          UYELEM(J)=0.
        ENDIF                                                             
        UMAGE(J)=SQRT(UXELEM(J)**2+UYELEM(J)**2)
c       DO N=1,NNODE
c         JCOUNT(LM(N))=JCOUNT(LM(N))+1
c         UMAGN(LM(N))=UMAGN(LM(N))+UMAGE(J)
c       ENDDO
600   CONTINUE                                                          
c     DO I=1,NUMNP
c       IF(JCOUNT(I).NE.0) THEN
c         UMAGN(I)=UMAGN(I)/REAL(JCOUNT(I))
c       ELSE
c         UMAGN(I)=0.
c       ENDIF
c     ENDDO
      DO I=1,NUMEL
        IF(UMAGE(I).LT. .1) THEN
          ZPLOT=RLOGP1
        ELSE
          ZPLOT=LOG(UMAGE(I))
        ENDIF
        ICOLOR=IZSET(NCOLOR,ZPLOT,VMIN,VDELT,0)
        DO L=1,4
          PX(L)=XX(KX(I,L))
          PY(L)=YY(KX(I,L))
        ENDDO
        CALL FILPAN(ICMAP1(ICOLOR),.FALSE.)
        CALL PANEL(4,PX,PY)
      ENDDO


C DRAW CONTOURS
      RMAX=5500.
      RMIN=0.1
      RLEVSP=500./3.
      LEVEL=(RMAX-RMIN)/RLEVSP+1
      DLEV=RLEVSP
      IF(LEVEL.LE.0) THEN
        LEVEL=2-LEVEL
        DLEV=-DLEV
      ENDIF
C DRAW CONTOURS OF PRESENT SURFACE
      ISKIP=1
      IF(ISKIP.EQ.1) GOTO 2300
      ICO=1
      CALL LINCLR(ICMAP2(ICO))
      LCONT=1
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
            IF(ZZZ(I).LT.VAL .AND. VAL.LE.ZZZ(J) .OR. 
     &         ZZZ(J).LT.VAL .AND. VAL.LE.ZZZ(I)) THEN
              FOUND=.TRUE.
              ICOUNT=ICOUNT+1
              XINT=XX(J)+(VAL-ZZZ(J))*(XX(I)-XX(J))/(ZZZ(I)-ZZZ(J))
              YINT=YY(J)+(VAL-ZZZ(J))*(YY(I)-YY(J))/(ZZZ(I)-ZZZ(J))
              XINT=XINT*1.
              YINT=YINT*1.
              IF(LCONT.EQ.1) THEN
                CALL MOVE(REAL(XINT),REAL(YINT))
                LCONT=2
              ELSE
                CALL DRAW(REAL(XINT),REAL(YINT))
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        ICO=ICO+1
        IF(ICO.GT.3) ICO=1
        CALL LINCLR(ICMAP2(ICO))
      ENDDO
C DRAW CONTOURS OF CALCULATED SURFACE
2300  CONTINUE
      ICO=1
      CALL LINCLR(ICMAP2(ICO))
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
              IF(LCONT.EQ.1) THEN
                CALL MOVE(REAL(XINT),REAL(YINT))
                LCONT=2
              ELSE
                CALL DRAW(REAL(XINT),REAL(YINT))
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        ICO=ICO+1
        IF(ICO.GT.3) ICO=1
        CALL LINCLR(ICMAP2(ICO))
      ENDDO
      CALL CLOSEG
      if(itime.eq.jtime) then
        pause
      endif


      GOTO 1
99    CONTINUE
      DO I=1,ISEG
        CALL SETVIS(I,.TRUE.)
      ENDDO
      write(*,*) 'input 1 to go around again, 0 to quit'
      read(*,*) igo
      if(igo.eq.1) goto 99
      CALL GRSTOP1
999   CONTINUE
      END                                                               
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
      FUNCTION IZSET(NCOLOR,ZZZ,ZF,ZD,IOFF)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      IMAX=NCOLOR+IOFF
      IMIN=IOFF+1
C       IC=IMAX-NINT((ZZZ-ZF)/ZD-.5)
        ZIC=(ZZZ-ZF)/ZD+.5
        IC=NINT(ZIC)
C       WRITE(*,*) ZZZ,ZIC,IC
      IF(IC.LT.IMIN) IC=IMIN
      IF(IC.GT.IMAX) IC=IMAX
      IZSET=IC
      RETURN
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
     
