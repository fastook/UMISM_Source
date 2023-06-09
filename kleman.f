      include "parameter.h"
c      PARAMETER(NMAX=MAXNUM)
      PARAMETER(NMAX=79999)
      IMPLICIT REAL*8(A-H,O-Z)
C **************************************************************C
C                                                               C
C   PROGRAM:  CONT3                                             C
C                                                               C
C   DATE:  11 23 87                                             C
C   PROGRAMMER:  FASTOOK                                        C
C                                                               C
C   FUNCTION:                                                   C
C            CONTOURS STAANDARD DATA SET (OUT3**) FOR SURFACE,  C
C            BED, FRACT, FLOWA, THICK ETC. OUTPUT IN STANDARD   C
C            PLOTTER FORM TO OUTC**                             C
C                                                               C
C **************************************************************C
C      DIMENSION X(44,40),Y(44,40),Z(44,40)
      DIMENSION XX(NMAX),YY(NMAX),ZZ(NMAX),CONST(NMAX),ACON(NMAX)
      DIMENSION KX(NMAX,4),HTICE(NMAX),BDROCK(NMAX),THIK(NMAX)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2),LM(4)             
      DIMENSION PSI(4),DPSI(4,2)                                        
      DIMENSION XY(2,4)
      CHARACTER HED*80,JUNK1*80,JUNK2*80,JUNK3*80,LABEL*26
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      CALL SETRIG
      IHEMI=0
c     WRITE(*,*) 'INPUT 0 FOR NORTHERN HEMISPHERE, 1 FOR SOUTH'
c     READ(*,*) IHEMI
      Istrad=0
c      write(*,*) 'input 1 for straddling greenwich, 0 for not'
c      read(*,*) Istrad
C
 1    READ(1,100,END=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &                    INTER,DT
      PRINT *,HED
100   FORMAT(A80,/,7I6,E15.6)
      NUM1=NUMCOL-1
      IF(IHEMI.EQ.0) THEN
        IOFF=0
        JOFF=0
c       PRINT *,'OFFSET: IN MINUTES LONGITUDE/LATITUDE'
c       PRINT *,'        FOR SCOTLAND USE 10,-2.5'
c       READ(*,*) IOFF,JOFF
      ELSE
        IOFF=0
        JOFF=0
      ENDIF
      OFFSETI=IOFF/60.
      OFFSETJ=JOFF/60.
      write(9,*) 'The following contains 2 sets:'
      write(9,'(1x,i6,a)') numnp,
     &                     ' nodal point values of surf,bed,thick'
      write(9,'(1x,i6,a)') numel,
     &                     ' element centroid values of ux,uy,umag'
      write(9,*) '------------------------------------------------------
     &--------------------'
      write(9,*) '     LONG      LAT         X          Y       SURF
     &   BED     THICK'
      DO NUM=1,NUMNP
        READ(1,200) N,KODE,XX(NUM),YY(NUM),HTICE(NUM),
     &              ADOT,FRACT,PSURF,BDROCK(NUM),FLOWA,SLDGB,
     &              TBED,itype,AFUDGE
        CALL RECPOL(.001*XX(NUM),.001*YY(NUM),RLAT,RLONG)
        IF(IHEMI.EQ.1) THEN
          RLAT=-RLAT
          RLONG=90-RLONG
          IF(RLONG.LT.0) RLONG=RLONG+360
        ENDIF
        RLONG=RLONG+OFFSETI
        RLAT=RLAT+OFFSETJ
        if(Istrad.eq.1) then
          if(rlong.gt.180.) rlong=rlong-360.
        endif
C CHANGE PLOTTED VARIABBLE HERE
        THIK(N)=max(0.,HTICE(N)-BDROCK(N))
        IF(BDROCK(N).LT.0.) THEN
          FLOT=BDROCK(N)*(1.-1.03/.917)
          IF(HTICE(N).LT.FLOT) THIK(NUM)=0.
        ENDIF
        WRITE(9,1000) RLONG,RLAT,XX(NUM),YY(NUM),HTICE(NUM),
     &             BDROCK(NUM),THIK(NUM)
1000  format(1x,2f9.4,2f11.0,3f10.3)
      ENDDO
      NNODE=4                                                         
      CENTX=0.0D00                                                    
      CENTY=0.0D00                                                    
      write(9,*) '     LONG      LAT         X          Y       UX      
     &    UY        UMAG'
      DO I=1,NUMEL
        READ(1,300) NUM,KX(NUM,1),KX(NUM,2),KX(NUM,3),KX(NUM,4),
     &              CONST(NUM),ACON(NUM)
        XPLOT=0.
        YPLOT=0.
        HH=0.0D00
        BD=0.0D00
        SUMX=0.0D00
        SUMY=0.0D00 
        DO ii=1,NNODE 
          LM(ii)=KX(NUM,ii)
          XY(1,ii)=XX(LM(ii))
          XY(2,ii)=YY(LM(ii))
        ENDDO
C                                                                       
        CALL SHAPEFEM(1,CENTX,CENTY,PSI,DPSI)                         
C CALCULATE DXDS...EQUATION (5.3.6)                                     
        DO ii=1,2                                                      
          DO ll=1,2
                               
            DXDS(ii,ll)=0.0
            DO K=1,NNODE  
              DXDS(ii,ll)=DXDS(ii,ll)+DPSI(K,ll)*XY(ii,K)
            ENDDO
          ENDDO
        ENDDO
   
C CALCULATE DSDX...EQUATION (5.2.7)                                     
        DETJ=(DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1))
        IF (DETJ.LE.0.0) then                                         
C       IF (DETJ.GE.0.0) then                                         
          WRITE (*,5544) DETJ,X                      
5544      FORMAT (' BAD JACOBIAN AT 750',E14.6,E14.6)
          STOP                                      
        endif
        DSDX(1,1)=DXDS(2,2)/DETJ                   
        DSDX(2,2)=DXDS(1,1)/DETJ                  
        DSDX(1,2)=-DXDS(1,2)/DETJ               
        DSDX(2,1)=-DXDS(2,1)/DETJ              
C CALCULATE D(PSI)/DX...EQUATION (5.3.5)                                
        DO ii=1,NNODE                                                  
          DPSIX(ii)=DPSI(ii,1)*DSDX(1,1)+DPSI(ii,2)*DSDX(2,1)
          DPSIY(ii)=DPSI(ii,1)*DSDX(1,2)+DPSI(ii,2)*DSDX(2,2)
        ENDDO                                                          
        DO ii=1,NNODE                                                  
          XPLOT=XPLOT+XX(LM(ii))*PSI(ii)
          YPLOT=YPLOT+YY(LM(ii))*PSI(ii)
          SUMX=SUMX + HTICE(LM(ii))*DPSIX(ii)              
          SUMY=SUMY + HTICE(LM(ii))*DPSIY(ii)             
          HH=HH + HTICE(LM(ii))*PSI(ii)                  
          BD=BD + BDROCK(LM(ii))*PSI(ii)               
        ENDDO                                                          
C                                                                      
        THICK=max(0.,HH-BD)
        IF(BD.LT.0.) THEN
          FLOT=BD*(1.-1.03/.917)
          IF(HH.LT.FLOT) THICK=0.
        ENDIF
        CALL RECPOL(.001*XPLOT,.001*YPLOT,RLAT,RLONG)
        IF(IHEMI.EQ.1) THEN
          RLAT=-RLAT
          RLONG=90-RLONG
          IF(RLONG.LT.0) RLONG=RLONG+360
        ENDIF
        RLONG=RLONG+OFFSET
        if(Istrad.eq.1) then
          if(rlong.gt.180) rlong=rlong-360.
        endif
        IF(THICK.GT.0.) THEN
          UX=-CONST(NUM)*SUMX/THICK
          UY=-CONST(NUM)*SUMY/THICK
        ELSE                      
          UX=0.                 
          UY=0.                
        ENDIF                 
        UMAG=SQRT(UX**2+UY**2)                                        
        WRITE(9,2000) RLONG,RLAT,XPLOT,YPLOT,UX,UY,UMAG
      ENDDO
2000  format(1x,2f9.4,2f11.0,1P3E12.4)
      DO N=1,NUMGBC
        READ(1,310) I,J,RJUNK
      ENDDO
999   continue
200   FORMAT(I6,I4,1P2E12.5,0PF10.2,F7.2,F9.2,F10.3,F10.1,F10.5,2F10.5,
     &       I5,F10.3)
300   FORMAT(5I6,1P2E17.10)
310   FORMAT(2I6,E13.6)

      END
c***************************************************
      
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
C*************************
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
