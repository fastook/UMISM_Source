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
      DIMENSION KX(NMAX,4),HTICE(NMAX),BDROCK(NMAX),GEOFLUX(NMAX)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2),LM(4)             
      DIMENSION PSI(4),DPSI(4,2)                                        
      DIMENSION XY(2,4)
      CHARACTER HED*80,JUNK1*80,JUNK2*80,JUNK3*80,LABEL*26
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      DATA RHOI /0.917D0/, RHOW /1.092D0/, GRAV /3.74d0/
      rhog=RHOI*GRAV/100
      CALL SETRIG
      WRITE(*,*) 'INPUT 0 FOR NORTHERN HEMISPHERE, 1 FOR SOUTH'
      READ(*,*) IHEMI
      Istrad=1
c      write(*,*) 'input 1 for straddling greenwich, 0 for not'
c      read(*,*) Istrad
C
      WRITE(*,*) 'INPUT 0 FOR HTICE'
      WRITE(*,*) '      1 FOR BDROCK'
      WRITE(*,*) '      2 FOR ADOT'
      WRITE(*,*) '      3 FOR FRACT'
      WRITE(*,*) '      4 FOR FLOW CONSTANT'
      WRITE(*,*) '      5 FOR DIFFERENCE'
      WRITE(*,*) '      6 FOR THICKNESS'
      WRITE(*,*) '      7 FOR MARGIN'
      WRITE(*,*) '      8 FOR PSURF'
      WRITE(*,*) '      9 FOR ZONE'
      WRITE(*,*) '     10 FOR SLIDING CONST'
      WRITE(*,*) '     11 FOR TEMPERATURE'
      WRITE(*,*) '     12 FOR AFUDGE '
      WRITE(*,*) '     13 FOR ELEMENT ACON '
      WRITE(*,*) '     14 FOR ELEMENT CONSTANT '
      WRITE(*,*) '     15 FOR ELEMENT VELOCITY '
      WRITE(*,*) '     16 FOR ELEMENT GRADIENT '
      WRITE(*,*) '     17 FOR BASAL WATER THICKNESS '
      WRITE(*,*) '     18 FOR BASAL MELT RATE '
      WRITE(*,*) '     19 FOR HTICE ON BEDROCK WITH BATHYMETRY '
      WRITE(*,*) '     20 FOR BASAL WATER POTENTIAL '
      WRITE(*,*) '     21 FOR GEOTHERMAL HEAT FLUX'
      WRITE(*,*) '     22 FOR DRIVING STRESS'
      READ(*,*) IPLOT
  350 FORMAT (I6)
      PRINT *,' THIS IS IPLOT',IPLOT
      CALL SCALPL(IPLOT,RLEVSP,RMAX,RMIN,LABEL)
 1    READ(1,100,END=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &                    INTER,DT
      PRINT *,HED
100   FORMAT(A80,/,7I6,E15.6)
      NUM1=NUMCOL-1
      IF(IHEMI.EQ.0) THEN
        PRINT *,'OFFSET: IN MINUTES LONGITUDE/LATITUDE'
        PRINT *,'        FOR SCOTLAND USE 10,-2.5'
        READ(*,*) IOFF,JOFF
      ELSE
        IOFF=0
        JOFF=0
      ENDIF
      OFFSETI=IOFF/60.
      OFFSETJ=JOFF/60.
      XMIN=1.E30
      XMAX=-XMIN
      YMIN=XMIN
      YMAX=XMAX
      RLATMIN=1.E30
      RLATMAX=-XMIN
      RLONGMIN=XMIN
      RLONGMAX=XMAX
      READ(2,*,IOSTAT=IOK) NN
      write(*,*) NN,' NN'
      IF(IOK.NE.0) THEN
        PRINT *,' NO DATA SET FOR WTHICK, BMELT '
      ELSEIF(NN.NE.NUMNP) THEN
        PRINT *,' PROBLEMS, DATA SETS DONT MATCH',nn,numnp
      ENDIF
      xcent=0.0
      ycent=0.0
      DO NUM=1,NUMNP
        IF(IOK.EQ.0) then
          READ(2,*) N,WTHICK,BMELT
        endif
c       READ(1,200) N,KODE,XX(NUM),YY(NUM),HTICE(NUM),
        READ(1,*) N,KODE,XX(NUM),YY(NUM),HTICE(NUM),
     &              ADOT,FRACT,PSURF,BDROCK(NUM),FLOWA,SLDGB,
     &              TBED,itype,AFUDGE,GEOFLUX(NUM)
        xcent=xcent+xx(num)
        ycent=ycent+yy(num)
        CALL RECPOL(.001*XX(NUM),.001*YY(NUM),RLAT,RLONG)
        IF(IHEMI.EQ.1) THEN
          RLAT=-RLAT
          RLONG=90-RLONG
          IF(RLONG.LT.0) RLONG=RLONG+360
        ENDIF
        RLONG=RLONG+OFFSETI
        RLAT=RLAT+OFFSETJ
        if(Istrad.eq.1) then
c         if(rlong.gt.180.) rlong=rlong-360.
          if(rlong.gt.180.) rlong=rlong-355.
        endif
        RLATMIN=MIN(RLATMIN,RLAT)
        RLATMAX=MAX(RLATMAX,RLAT)
        RLONGMIN=MIN(RLONGMIN,RLONG)
        RLONGMAX=MAX(RLONGMAX,RLONG)
        XMIN=MIN(XMIN,XX(NUM))
        YMIN=MIN(YMIN,YY(NUM))
        XMAX=MAX(XMAX,XX(NUM))
        YMAX=MAX(YMAX,YY(NUM))
C CHANGE PLOTTED VARIABBLE HERE
        IF(IPLOT.EQ.0) ZZ(N)=HTICE(N)
        IF(IPLOT.EQ.1) ZZ(N)=BDROCK(N)
c       IF(IPLOT.EQ.2) ZZ(N)=ADOT        *1e6 ! micrometers/yr
c       IF(IPLOT.EQ.2) ZZ(N)=ADOT        *1e3 ! millimeters/yr
        IF(IPLOT.EQ.2) ZZ(N)=ADOT        *1e2 ! centimeters/yr
        IF(IPLOT.EQ.3) ZZ(N)=FRACT
        IF(IPLOT.EQ.4) ZZ(N)=FLOWA
        IF(IPLOT.EQ.5) ZZ(N)=HTICE(N)-PSURF
        IF(IPLOT.EQ.6 .OR. IPLOT.EQ.7) THEN
          ZZ(N)=max(0.,HTICE(N)-BDROCK(N))
c         IF(BDROCK(N).LT.0.) THEN
c           FLOT=BDROCK(N)*(1.-1.03/.917)
c           IF(HTICE(N).LT.FLOT) ZZ(NUM)=0.
c         ENDIF
        ENDIF
        IF(IPLOT.EQ.8) ZZ(N)=PSURF
        IF(IPLOT.EQ.9) ZZ(N)=ADOT
        IF(IPLOT.EQ.10) ZZ(N)=SLDGB
        IF(IPLOT.EQ.11) ZZ(N)=TBED
        IF(IPLOT.EQ.12) ZZ(N)=AFUDGE
        IF(IPLOT.EQ.17) ZZ(N)=WTHICK
        IF(IPLOT.EQ.18) ZZ(N)=BMELT*1000.
        IF(IPLOT.EQ.19) THEN
          ZZ(N)=HTICE(N)
c         IF(BDROCK(N).LT.0.) THEN
c           FLOT=BDROCK(N)*(1.-1.03/.917)
c           IF(HTICE(N).LT.FLOT) ZZ(NUM)=BDROCK(N)
c         ENDIF
        ENDIF
        IF(IPLOT.EQ.20) ZZ(N)=HTICE(N)+(1.03-.917)*BDROCK(N)/.917
        IF(IPLOT.EQ.21) ZZ(N)=GEOFLUX(N)*13.25e-5
        IF(IPLOT.GE.0 .AND. IPLOT.LE.12) THEN
          WRITE(9,123) real(RLONG),real(RLAT),real(ZZ(NUM))
123   format(1x,1p3e14.6)
          if(adot.eq.-123) then
            WRITE(19,*) real(RLONG),real(RLAT),real(ZZ(NUM))
          endif
          WRITE(10,*) real(XX(NUM))*.001,real(YY(NUM))*.001,
     &                real(ZZ(NUM))
        ENDIF
        IF(IPLOT.GE.17 .AND. IPLOT.LE.21) THEN
          WRITE(9,123) real(RLONG),real(RLAT),real(ZZ(NUM))
          WRITE(10,*) real(XX(NUM))*.001,real(YY(NUM))*.001,
     &                real(ZZ(NUM))
        ENDIF
      ENDDO
      xcent=xcent/num
      ycent=ycent/num
      CALL RECPOL(0.001*xcent,0.001*ycent,rlatc,rlongc)
      print *,xcent,ycent
      print *,rlongc,rlatc
c      pause
c      PRINT *,' RLATMIN ',RLATMIN
c      PRINT *,' RLATMAX ',RLATMAX
c      PRINT *,' RLONGMIN ',RLONGMIN
c      PRINT *,' RLONGMAX ',RLONGMAX
      DISTMAX=0
      DISTMIN=1E30
      VMAX=0
      VMIN=1E30
      DO I=1,NUMEL
        READ(1,300) NUM,KX(NUM,1),KX(NUM,2),KX(NUM,3),KX(NUM,4),
     &              CONST(NUM),ACON(NUM)
        CONST(NUM)=CONST(NUM)
        if(.false.) then
          XPLOT=0.
          YPLOT=0.
          HH=0.
          BD=0.
          DO J=1,4
            XPLOT=XPLOT+XX(KX(NUM,J))
            YPLOT=YPLOT+YY(KX(NUM,J))
            HH=HH+HTICE(KX(NUM,J))
            BD=BD+BDROCK(KX(NUM,J))
          ENDDO
          XPLOT=0.25*XPLOT
          YPLOT=0.25*YPLOT
          HH=0.25*HH
          BD=0.25*BD
          H21=HTICE(KX(NUM,2))-HTICE(KX(NUM,1))
          X21=XX(KX(NUM,2))-XX(KX(NUM,1))
          H34=HTICE(KX(NUM,3))-HTICE(KX(NUM,4))
          X34=XX(KX(NUM,3))-XX(KX(NUM,4))
          H41=HTICE(KX(NUM,4))-HTICE(KX(NUM,1))
          Y41=YY(KX(NUM,4))-YY(KX(NUM,1))
          H32=HTICE(KX(NUM,3))-HTICE(KX(NUM,4))
          Y32=YY(KX(NUM,3))-YY(KX(NUM,2))
          DHDX=(H21/X21+H34/X34)*0.5
          DHDY=(H41/Y41+H32/Y32)*0.5
          DH=SQRT(DHDX**2+DHDY**2)
        else

          NNODE=4                                                         
          XPLOT=0.
          YPLOT=0.
          CENTX=0.0D00                                                    
          CENTY=0.0D00                                                    
          HH=0.0D00                                                         
          DH=0.0D00                                                         
          BD=0.0D00                                                         
          SUMX=0.0D00                                                       
          SUMY=0.0D00                                                       
          DO ii=1,NNODE                                                  
            LM(ii)=KX(NUM,ii)                                                     
          ENDDO
          ii=LM(1)                                                           
          JJ=LM(2)                                                          
          kk=LM(3)                                                           
          ll=LM(4)                                                           
          XY(1,1)=XX(ii)                                                     
          XY(1,2)=XX(JJ)                                                    
          XY(1,3)=XX(kk)                                                     
          XY(1,4)=XX(ll)                                   
          XY(2,1)=YY(ii)                                                     
          XY(2,2)=YY(JJ)                                                    
          XY(2,3)=YY(kk)                                                     
          XY(2,4)=YY(ll)                                   
          XCENT=(XX(ii)+XX(JJ)+XX(kk)+XX(ll))/4000.                          
          YCENT=(YY(ii)+YY(JJ)+YY(kk)+YY(ll))/4000.                          
C                                                                       
          CALL SHAPE(1,CENTX,CENTY,PSI,DPSI)                         
      
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
C         IF (DETJ.GE.0.0) then                                         
            WRITE (*,5544) DETJ,X                                             
 5544       FORMAT (' BAD JACOBIAN AT 750',E14.6,E14.6)                       
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
          enddo
          CALL RECPOL(.001*XPLOT,.001*YPLOT,RLAT,RLONG)
C                                                                      
          DH=SUMX**2 + SUMY**2                                            
          DH=SQRT(DH)                                                   
        endif


        THICK=HH-BD
c        IF(BD.LT.0.) THEN
c          FLOT=BD*(1.-1.03/.917)
c          IF(HH.LT.FLOT) THICK=0.
c        ENDIF
        IF(IHEMI.EQ.1) THEN
          RLAT=-RLAT
          RLONG=90-RLONG
          IF(RLONG.LT.0) RLONG=RLONG+360
        ENDIF
        RLONG=RLONG+OFFSETI
        RLAT=RLAT+OFFSETJ
        if(Istrad.eq.1) then
c         if(rlong.gt.180.) rlong=rlong-360.
          if(rlong.gt.180.) rlong=rlong-354.5
        endif
        IF(IPLOT.EQ.13) THEN
          WRITE(9,123) RLONG,RLAT,ACON(NUM)
          WRITE(10,*) XPLOT*.001,YPLOT*.001,ACON(NUM)
        ELSEIF(IPLOT.EQ.14) THEN
          WRITE(9,123) RLONG,RLAT,CONST(NUM)
          WRITE(10,*) XPLOT*.001,YPLOT*.001,CONST(NUM)
        ELSEIF(IPLOT.EQ.15) THEN
          IF(THICK.GT.0.) THEN
            UX=-CONST(NUM)*SUMX/THICK
            UY=-CONST(NUM)*SUMY/THICK
          ELSE                                                              
            UX=0.                                                           
            UY=0.                                                           
          ENDIF                                                             
          UMAG=SQRT(UX**2+UY**2)                                        
          VEL=100*10*UMAG ! mm/yr
c         VEL=UMAG        ! m/yr
          vel=log10(1e3*umag+1e-6)
c         if(vel.gt.-5) vel=log10(1e3*umag)
          if(.false.) then
            if(umag.gt.0) then
              vel=log10(1e3*umag)
            else
              vel=0
            endif
          endif
          VMIN=MIN(VMIN,VEL)
          VMAX=MAX(VMAX,VEL)
          WRITE(9,123) RLONG,RLAT,VEL
          WRITE(10,*) XPLOT*.001,YPLOT*.001,VEL
        ELSEIF(IPLOT.EQ.16) THEN
          VMIN=MIN(VMIN,DH*1000.)
          VMAX=MAX(VMAX,DH*1000.)
          WRITE(9,123) RLONG,RLAT,DH*1000.
          WRITE(10,*) XPLOT*.001,YPLOT*.001,DH*1000.
        ELSEIF(IPLOT.EQ.22) THEN
          alpha=sqrt(sumx**2+sumy**2)
          WRITE(9,123) RLONG,RLAT,rhog*thick*alpha
          WRITE(10,*) XPLOT*.001,YPLOT*.001,rhog*thick*alpha
        ENDIF
        DIST1=(XX(KX(NUM,1))-XX(KX(NUM,3)))**2+
     &        (YY(KX(NUM,1))-YY(KX(NUM,3)))**2
        DIST2=(XX(KX(NUM,2))-XX(KX(NUM,4)))**2+
     &        (YY(KX(NUM,2))-YY(KX(NUM,4)))**2
        DIST3=(XX(KX(NUM,2))-XX(KX(NUM,1)))**2+
     &        (YY(KX(NUM,2))-YY(KX(NUM,1)))**2
        DIST4=(XX(KX(NUM,3))-XX(KX(NUM,2)))**2+
     &        (YY(KX(NUM,3))-YY(KX(NUM,2)))**2
        DIST5=(XX(KX(NUM,4))-XX(KX(NUM,3)))**2+
     &        (YY(KX(NUM,4))-YY(KX(NUM,3)))**2
        DIST6=(XX(KX(NUM,1))-XX(KX(NUM,4)))**2+
     &        (YY(KX(NUM,1))-YY(KX(NUM,4)))**2
        DISTMAX=MAX(DISTMAX,DIST1,DIST2,DIST3,DIST4,DIST5,DIST6)
        DISTMIN=MIN(DISTMIN,DIST1,DIST2,DIST3,DIST4,DIST5,DIST6)
      ENDDO
      XYDMAX=SQRT(DISTMAX)*0.001
      XYDMIN=SQRT(DISTMIN)*0.001
      PRINT *,XYDMAX,XYDMIN
      DISTMAX=SQRT(DISTMAX)*0.001/RKMPDEG
      DISTMIN=SQRT(DISTMIN)*0.001/RKMPDEG
      PRINT *,'DISTMAX',DISTMAX,DISTMIN
      DO N=1,NUMGBC
        READ(1,310) I,J,RJUNK
      ENDDO
C
      iskip=1
      if(iskip.eq.0) then
cVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
      CALL RECPOL(.001*XMIN,.001*YMIN,RLAT,RLONG)
      IF(IHEMI.EQ.1) THEN
        RLAT=-RLAT
        RLONG=90-RLONG
        IF(RLONG.LT.0) RLONG=RLONG+360
      ENDIF
      RLONG=RLONG+OFFSETI
        RLAT=RLAT+OFFSETJ
      RLATMIN=MIN(RLATMIN,RLAT)
      RLATMAX=MAX(RLATMAX,RLAT)
      RLONGMIN=MIN(RLONGMIN,RLONG)
      RLONGMAX=MAX(RLONGMAX,RLONG)
      CALL RECPOL(.001*XMAX,.001*YMIN,RLAT,RLONG)
      IF(IHEMI.EQ.1) THEN
        RLAT=-RLAT
        RLONG=90-RLONG
        IF(RLONG.LT.0) RLONG=RLONG+360
      ENDIF
      RLONG=RLONG+OFFSETI
        RLAT=RLAT+OFFSETJ
      RLATMIN=MIN(RLATMIN,RLAT)
      RLATMAX=MAX(RLATMAX,RLAT)
      RLONGMIN=MIN(RLONGMIN,RLONG)
      RLONGMAX=MAX(RLONGMAX,RLONG)
      CALL RECPOL(.001*XMIN,.001*YMAX,RLAT,RLONG)
      IF(IHEMI.EQ.1) THEN
        RLAT=-RLAT
        RLONG=90-RLONG
        IF(RLONG.LT.0) RLONG=RLONG+360
      ENDIF
      RLONG=RLONG+OFFSETI
        RLAT=RLAT+OFFSETJ
      RLATMIN=MIN(RLATMIN,RLAT)
      RLATMAX=MAX(RLATMAX,RLAT)
      RLONGMIN=MIN(RLONGMIN,RLONG)
      RLONGMAX=MAX(RLONGMAX,RLONG)
      CALL RECPOL(.001*XMAX,.001*YMIN,RLAT,RLONG)
      IF(IHEMI.EQ.1) THEN
        RLAT=-RLAT
        RLONG=90-RLONG
        IF(RLONG.LT.0) RLONG=RLONG+360
      ENDIF
      RLONG=RLONG+OFFSETI
        RLAT=RLAT+OFFSETJ
      RLATMIN=MIN(RLATMIN,RLAT)
      RLATMAX=MAX(RLATMAX,RLAT)
      RLONGMIN=MIN(RLONGMIN,RLONG)
      RLONGMAX=MAX(RLONGMAX,RLONG)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      endif
C
      LATMIN=NINT(RLATMIN)
      LATMAX=NINT(RLATMAX)
      LONGMIN=NINT(RLONGMIN)
      LONGMAX=NINT(RLONGMAX)
      IF(LATMIN.GT.LATMAX) CALL ISWAP(LATMIN,LATMAX)
      IF(LONGMIN.GT.LONGMAX) CALL ISWAP(LONGMIN,LONGMAX)
C gmt.e file, unit 11
      WRITE(JUNK1,900) RLONGMIN,RLONGMAX,RLATMIN,RLATMAX
900   FORMAT('-R',F7.3,'/',F7.3,'/',F7.3,'/',F7.3)
901   FORMAT(' ',F7.3,'/',F7.3,'/',F7.3,'/',F7.3)
      CALL STRIP(80,JUNK1,N)
      WRITE(11,1000) JUNK1(1:N)
1000  FORMAT('pscoast -V -U ',A,' -Jm.1 ',
     &       '-B10g5/10g5',
     &       ' -G220 -P -X1 -Y1 -W5 -K  > test.ps')
C
      WRITE(JUNK1,905) 2*DISTMAX,2*DISTMAX
905   FORMAT('-I',F10.3,'/',F10.3)
      WRITE(JUNK2,906) 4.*DISTMAX
906   FORMAT('-S',F10.3)
      CALL STRIP(80,JUNK1,N1)
      CALL STRIP(80,JUNK2,N2)
      WRITE(11,1010) JUNK1(1:N1)
1010  FORMAT('#blockmean gmt.xyz -V ',A,' -R > temp.xyz')
C
      WRITE(11,1020) JUNK1(1:N1)
1020  FORMAT('#surface temp.xyz -V -A.707 -Gout.grd ',A,' -R -T.75')
C
      WRITE(11,1025) JUNK1(1:N1),JUNK2(1:N2)
c 1025  FORMAT('nearneighbor gmt.xyz -E0 -V -Gout.grd ',A,' -R ',A)
1025  FORMAT('nearneighbor gmt.xyz -V -Gout.grd ',A,' -R ',A)
C
      WRITE(JUNK1,910) '-C',RLEVSP
      IF(IPLOT.EQ.3) THEN
        WRITE(JUNK2,910) '-A',RLEVSP
      ELSE
c        WRITE(JUNK2,910) '-A',RLEVSP*2.
        WRITE(JUNK2,910) '-A',RLEVSP
      ENDIF
      WRITE(JUNK3,920) RMIN,RMAX
910   FORMAT(A,F10.3)
920   FORMAT('-L',F10.3,'/',F10.3)
      CALL STRIP(80,JUNK1,N1)
      CALL STRIP(80,JUNK2,N2)
      CALL STRIP(80,JUNK3,N3)
c      WRITE(*,1030) JUNK1(1:N1),JUNK2(1:N2),JUNK3(1:N3)
      WRITE(11,1030) JUNK1(1:N1),JUNK2(1:N2),JUNK3(1:N3)
      IF(IPLOT.EQ.6) THEN
        WRITE(JUNK1,910) '-C',1.
        WRITE(JUNK2,910) '-A',1.
        WRITE(JUNK3,920) 1.,1.
        CALL STRIP(80,JUNK1,N1)
        CALL STRIP(80,JUNK2,N2)
        CALL STRIP(80,JUNK3,N3)
c        WRITE(*,1030) JUNK1(1:N1),JUNK2(1:N2),JUNK3(1:N3)
        WRITE(11,1030) JUNK1(1:N1),JUNK2(1:N2),JUNK3(1:N3)
      ENDIF
1030  FORMAT('grdcontour -O -V  -S5 out.grd -Jm -R ',A,
     &       ' -B -P ',A,
     &       ' -Wa -T ',A,
     &       ' -K  >> test.ps')
      WRITE(11,1040)
1040  FORMAT('pstext -P -R0/11/0/8.5 -Jx1 -O -K << END >> test.ps')
      WRITE(11,1050) 6.5,LABEL
1050  FORMAT('0 ',f4.1,' 20 0 1 1 ',A)
      WRITE(11,1060)
1060  FORMAT('END')
      WRITE(11,1041)
1041  FORMAT('pstext -P -R0/11/0/8.5 -Jx1 -O  << END >> test.ps')
      WRITE(11,1051) 6.2,HED
1051  FORMAT('0.2 ',f4.1,' 10 0 1 1 ',A)
      WRITE(11,1061)
1061  FORMAT('END')
      WRITE(11,1070)
c  1070  FORMAT('xpsview -wp -skipc test.ps &')
1070  FORMAT('display test.ps &')
c
C img.e file, unit 12
      WRITE(*,900) RLONGMIN,RLONGMAX,RLATMIN,RLATMAX
      RLONGMIN=1*nint(RLONGMIN/1.)
      RLONGMAX=1*nint(RLONGMAX/1.)
      RLATMIN=1*nint(RLATMIN/1.)
      RLATMAX=1*nint(RLATMAX/1.)
      WRITE(*,900) RLONGMIN,RLONGMAX,RLATMIN,RLATMAX
      DISTMIN=0.01*nint(DISTMIN/0.01)
      DISTMAX=0.01*nint(DISTMAX/0.01)
      DISTX=(XMAX-XMIN)
      DISTY=(YMAX-YMIN)
      PRINT *,DISTX,DISTY,DISTMIN,DISTMAX
      nincx=nint(-(RLONGMIN-RLONGMAX)/DISTMIN)
      xwid=nincx*DISTMIN
      nincy=nint(-(RLATMIN-RLATMAX)/DISTMIN)
      ywid=nincy*DISTMIN
      print *,xwid,nincx,ywid,nincy
      RLONGMAX=RLONGMIN+xwid
      RLATMAX=RLATMIN+ywid
      WRITE(*,900) RLONGMIN,RLONGMAX,RLATMIN,RLATMAX
      WRITE(JUNK1,900) RLONGMIN,RLONGMAX,RLATMIN,RLATMAX
      CALL STRIP(80,JUNK1,N1)
      WRITE(JUNK2,900) RLONGMIN,RLATMIN,RLONGMAX,RLATMAX
      CALL STRIP(80,JUNK2,N2)
902   FORMAT(' ',F7.3,'/',F7.3,'/1:',I10)
      if(IHEMI.eq.0) then
c        WRITE(JUNK3,901) (RLONGMIN+RLONGMAX)*.5,90.,2.,60.
        WRITE(JUNK3,902) (RLONGMIN+RLONGMAX)*.5,90.,
     &                    nint(10*MAX(DISTx,DISTy))
      else
c        WRITE(JUNK3,901) (RLONGMIN+RLONGMAX)*.5,-90.,2.,-60.
        WRITE(JUNK3,902) (RLONGMIN+RLONGMAX)*.5,-90.,
     &                    nint(10*MAX(DISTx,DISTy))
      endif
      CALL STRIP(80,JUNK3,N3)
c psbasemap
      WRITE(22,1080) JUNK1(1:N1)
1080  FORMAT('psbasemap -V -U ',A,' -Jm.08 ',
     &       '-B10g5/10g5',
     &       ' -P -X1 -Y2 -K  > test.ps')
      WRITE(12,1081) JUNK1(1:N1),JUNK3(1:N3)
      WRITE(12,1082) JUNK2(1:N2),JUNK3(1:N3)
1081  FORMAT('#psbasemap -V -U ',A,' -Js',A,
     &       ' -B10g5/10g5',
     &       ' -P -X1 -Y3 -K  > test.ps')
1082  FORMAT('psbasemap -V -U ',A,'r -Js',A,
     &       ' -B10g5/10g5',
     &       ' -P -X1 -Y3 -K  > test.ps')
C
      WRITE(JUNK1,905) DISTMIN,DISTMAX
      WRITE(JUNK2,906) 4.*DISTMAX
      CALL STRIP(80,JUNK1,N1)
      CALL STRIP(80,JUNK2,N2)
C nearneighbor
c      WRITE(JUNK1,905) DISTMIN/2.,DISTMAX/2.
      WRITE(JUNK1,905) DISTMIN,DISTMIN
      CALL STRIP(80,JUNK1,N1)
      WRITE(22,1025) JUNK1(1:N1),JUNK2(1:N2)
      WRITE(12,1025) JUNK1(1:N1),JUNK2(1:N2)
      if(iplot.eq.19) WRITE(12,1026) JUNK1(1:N1),JUNK2(1:N2)
1026  FORMAT('#nearneighbor /d/people/shamis/topo/latlong.nhem -V',
     &       ' -Getopo5.grd ',A,' -R ',A)

C grdmath
      WRITE(22,1090)
      if(iplot.eq.19) then
        WRITE(12,1089)
        WRITE(12,1090)
      else
        WRITE(12,1091)
      endif
1089  FORMAT('cp out.grd tmp.grd')
1091  FORMAT('grdmath out.grd 1 MUL = tmp.grd')
1090  FORMAT('#grdmath out.grd etopo5.grd AND = tmp.grd')
C grd2cpt or makecpt
      WRITE(22,1092)
      WRITE(12,1092)
c 1092  FORMAT('grd2cpt tmp.grd -Z -V -Crelief > topo.cpt')
1092  FORMAT('makecpt -Z -T-2000/2000/250 -V -Crainbow > topo.cpt')
      WRITE(22,1093)
      WRITE(12,1093)
c 1093  FORMAT('grd2cpt tmp.grd -Z -V -L-2000/2000 -S-2000/2000/250 -V 
c      &-Cno_green > topo.cpt')
c 1093  FORMAT('grd2cpt tmp.grd -Z -V -Cno_green > topo.cpt')
1093  FORMAT('grd2cpt tmp.grd -Z -V -Crainbow > topo.cpt')
      WRITE(22,1095)
      WRITE(12,1095)
1095  FORMAT('grdgradient tmp.grd -A45 -Gtmpt.grd -Nt -V')
c grdimage
      WRITE(22,1100)
1100  FORMAT('grdimage tmp.grd -Itmpt.grd -V -R -Jm -O ',
     &       '-K -P -U -V -Ctopo.cpt',
     &       '  >> test.ps')
      WRITE(12,1101)
1101  FORMAT('grdimage tmp.grd -Itmpt.grd -V -R -Js -O ',
     &       '-K -P -U -V -Ctopo.cpt',
     &       '  >> test.ps')
c psxy for vectors
      WRITE(22,1105)
1105  FORMAT('psxy -R -Jm -O -K -M velo.vector -W0.1 -V >> test.ps')
      WRITE(12,1106)
1106  FORMAT('psxy -R -Js -O -K -M velo.vector -W0.1 -V >> test.ps')
c pscoast
      WRITE(22,1110)
1110  FORMAT('pscoast -Na -Lfx5.4/1.8/50/500 -V -U -R -Jm -P -W2 -K -O 
     &>> test.ps')
      WRITE(12,1111)
1111  FORMAT('pscoast -Na -Lfx5.4/1.8/50/500 -V -U -R -Js -P -W2 -K -O 
     &>> test.ps')
C grdcontour
      WRITE(JUNK1,910) '-C',RLEVSP
      IF(IPLOT.EQ.3) THEN
        WRITE(JUNK2,910) '-A',RLEVSP
      ELSE
c        WRITE(JUNK2,910) '-A',RLEVSP*2.
        WRITE(JUNK2,910) '-A',RLEVSP
      ENDIF
      WRITE(JUNK3,920) RMIN,RMAX
      CALL STRIP(80,JUNK1,N1)
      CALL STRIP(80,JUNK2,N2)
      CALL STRIP(80,JUNK3,N3)
c      WRITE(*,1031) JUNK1(1:N1),JUNK2(1:N2),JUNK3(1:N3)
      WRITE(22,1031) JUNK1(1:N1),JUNK2(1:N2),JUNK3(1:N3)
      WRITE(12,1032) JUNK1(1:N1),JUNK2(1:N2),JUNK3(1:N3)
      IF(IPLOT.EQ.6) THEN
        WRITE(JUNK1,910) '-C',1.
        WRITE(JUNK2,910) '-A',1.
        WRITE(JUNK3,920) 1.,1.
        CALL STRIP(80,JUNK1,N1)
        CALL STRIP(80,JUNK2,N2)
        CALL STRIP(80,JUNK3,N3)
c        WRITE(*,1031) JUNK1(1:N1),JUNK2(1:N2),JUNK3(1:N3)
        WRITE(22,1031) JUNK1(1:N1),JUNK2(1:N2),JUNK3(1:N3)
        WRITE(12,1032) JUNK1(1:N1),JUNK2(1:N2),JUNK3(1:N3)
      ENDIF
1031  FORMAT('grdcontour -O -V  -S5 out.grd -Jm -R ',A,
     &       ' -B -P ',A,
     &       ' -Wa -T ',A,
     &       ' -K  >> test.ps')
1032  FORMAT('grdcontour -O -V  -S5 tmp.grd -Js -R ',A,
     &       ' -B -P ',A,
     &       ' -Wa -T ',A,
     &       ' -K  >> test.ps')
c psscale
      WRITE(JUNK1,*) LABEL
      CALL STRIP(80,JUNK1,N1)
      WRITE(22,1120)
      WRITE(12,1120)
1120  FORMAT('psscale -E -V -Ctopo.cpt -D5/4/3.3/0.3',
     &       ' -O -K -L >> test.ps')
c
C pstext's
      WRITE(22,1040)
      WRITE(22,1050) 6.5,LABEL
      WRITE(22,1060)
      WRITE(22,1041)
      WRITE(22,1051) 6.2,HED
      WRITE(22,1061)
      WRITE(22,1070)
      WRITE(12,1040)
      WRITE(12,1050) 6.5,LABEL
      WRITE(12,1060)
      WRITE(12,1041)
      WRITE(12,1051) 6.2,HED
      WRITE(12,1061)
      WRITE(12,1070)
C GMT.e file, unit 14
c .... make x,y scales the same for southern hemisphere grids ...
      xxmin= -2713.6*1000
      xxmax=  2761.4*1000
      yymin=(-2304.0+5)*1000
      yymax= (2246+5)*1000
      xmin=max(xmin,xxmin)
      xmax=min(xmax,xxmax)
      ymin=max(ymin,yymin)
      ymax=min(ymax,yymax)

      xmid=0.5*(xmax+xmin)
      ymid=0.5*(ymax+ymin)
      if(xmax-xmin.lt.ymax-ymin) then
        ymax=ymid+0.5*(xmax-xmin)
        ymin=ymid-0.5*(xmax-xmin)
      else
        xmax=xmid+0.5*(ymax-ymin)
        xmin=xmid-0.5*(ymax-ymin)
      endif
      xmin1=0.001*xmin
      xmax1=0.001*xmax
      ymin1=0.001*ymin
      ymax1=0.001*ymax
c      xmin1=100*nint(xmin1/100.)
c      xmax1=100*nint(xmax1/100.)
c      ymin1=100*nint(ymin1/100.)
c      ymax1=100*nint(ymax1/100.)
c      xydmin=10*nint(xydmin/10.)
      if(.false.) then
        nincx=nint((xmax1-xmin1)/xydmin)
        xwid=nincx*xydmin
        nincy=nint((ymax1-ymin1)/xydmin)
        ywid=nincy*xydmin
        xmax1=xmin1+xwid
        ymax1=ymin1+ywid
      endif
c      xydmin=(xmax1-xmin1)/ninc
      WRITE(JUNK1,2900) xmin1,xmax1,ymin1,ymax1
      if(ihemi.eq.1) then
        WRITE(16,'(a)') 'extract.e << END'
        WRITE(16,*) real(xmin1),real(xmax1),real(ymin1),real(ymax1)
        WRITE(16,'(a)') 'END'
        WRITE(16,*)
      endif
2900  FORMAT('-R',F7.1,'/',F7.1,'/',F7.1,'/',F7.1)
      CALL STRIP(80,JUNK1,N1)
c psbasemap
      WRITE(14,2080) JUNK1(1:N1)
2080  FORMAT('psbasemap -V -U ',A,' -JX5/5 ',
     &       '-B200g200/200g200',
     &       ' -P -X1 -Y2 -K  > test.ps')
      if(ihemi.eq.1) WRITE(16,2180) JUNK1(1:N1)
2180  FORMAT('psbasemap -V -U ',A,' -JX6/6 ',
     &       '-B500g50/500g50',
     &       ' -P -X1 -Y2 -K  > test.ps')
C
      WRITE(JUNK1,905) XYDMIN/1,XYDMIN/1
      WRITE(JUNK2,906) XYDMIN*2
      CALL STRIP(80,JUNK1,N1)
      CALL STRIP(80,JUNK2,N2)
C nearneighbor
      WRITE(14,2025) JUNK1(1:N1),JUNK2(1:N2)
      if(ihemi.eq.1) then
        WRITE(16,2025) JUNK1(1:N1),JUNK2(1:N2)
        WRITE(16,*)
      ENDIF
c 2025  FORMAT('nearneighbor gmt.XYZ -E0 -V -Gout.grd ',A,' -R ',A)
2025  FORMAT('nearneighbor gmt.XYZ -V -Gout.grd ',A,' -R ',A)
C grd2cpt or makecpt
      WRITE(14,2092)
      if(ihemi.eq.1) WRITE(16,2092)
2092  FORMAT('makecpt -Z -T-3/3/0.500 -V -Cno_green > topo.cpt')
      WRITE(14,2093)
      if(ihemi.eq.1) WRITE(16,2093)
c 2093  FORMAT('grd2cpt out.grd -Z -V -L-2000/2000 -S-2000/2000/250 -V 
c      &-Cno_green > topo.cpt')
2093  FORMAT('grd2cpt out.grd -Z -V -Crainbow > topo.cpt')
      WRITE(14,2026)
      if(ihemi.eq.1) WRITE(16,2026)
2026  FORMAT('grdgradient out.grd -A45 -Goutt.grd -Nt -V')
c grdimage
      if(iplot.eq.0 .or. iplot.eq.1 .or. 
     &   iplot.eq.6 .or. iplot.eq.8 .or. iplot.eq.19 .or.
     &   iplot.eq.20 .or. iplot.eq.21) then
        WRITE(14,2101)
        if(ihemi.eq.1) WRITE(16,2101)
      else
        WRITE(14,2102)
        if(ihemi.eq.1) WRITE(16,2102)
      endif
2101  FORMAT('grdimage out.grd -Ioutt.grd -V -R -JX -O -K -P -U ',
     &       '-V -Ctopo.cpt  >> test.ps')
2102  FORMAT('grdimage out.grd -V -R -JX -O -K -P -U ',
     &       '-V -Ctopo.cpt  >> test.ps')
c psxy for vectors
      WRITE(14,2106)
      if(ihemi.eq.1) then
        WRITE(16,2116)
        WRITE(16,2106)
2116  FORMAT(/'psimage asdf.ras -V -P -W6 -V -O -K -C0/0 >> test.ps'/)
        WRITE(16,2117)
2117  FORMAT(/'grdcontour -O -V  -S5 out.grd -JX -R -B -P  -A- ',
     &        ' -W+c5 -Ctopo.cpt -T -K  >> test.ps'/)
      endif
2106  FORMAT('psxy -R -JX -O -K -M velo.xy -W0.1 -V >> test.ps')
C grdcontour
      WRITE(JUNK1,910) '-C',RLEVSP
      IF(IPLOT.EQ.3) THEN
        WRITE(JUNK2,910) '-A',RLEVSP
      ELSE
c        WRITE(JUNK2,910) '-A',RLEVSP*2.
        WRITE(JUNK2,910) '-A',RLEVSP
      ENDIF
      WRITE(JUNK3,920) RMIN,RMAX
      CALL STRIP(80,JUNK1,N1)
      CALL STRIP(80,JUNK2,N2)
      CALL STRIP(80,JUNK3,N3)
      IF(IPLOT.EQ.15) THEN
        WRITE(14,2031) '-C1','-A1','-L1/2'
        WRITE(14,2031) '-C10','-A10','-L10/11'
        WRITE(14,2031) '-C100','-A100','-L100/101'
        WRITE(14,2031) '-C200','-A200','-L200/201'
        WRITE(14,2031) '-C500','-A500','-L500/501'
        WRITE(14,2031) '-C1000','-A1000','-L1000/1001'
        WRITE(14,2031) '-C1500','-A1500','-L1500/1501'
        WRITE(14,2031) '-C5000','-A5000','-L5000/5001'
      ELSE
        WRITE(14,2031) JUNK1(1:N1),JUNK2(1:N2),JUNK3(1:N3)
      ENDIF
      IF(IPLOT.EQ.6) THEN
        WRITE(JUNK1,910) '-C',1.
        WRITE(JUNK2,910) '-A',1.
        WRITE(JUNK3,920) 1.,1.
        CALL STRIP(80,JUNK1,N1)
        CALL STRIP(80,JUNK2,N2)
        CALL STRIP(80,JUNK3,N3)
        WRITE(14,2031) JUNK1(1:N1),JUNK2(1:N2),JUNK3(1:N3)
      ENDIF
2031  FORMAT('grdcontour -O -V  -S5 out.grd -JX -R ',A,
     &       ' -B -P ',A,
     &       ' -Wa -T ',A,
     &       ' -K  >> test.ps')
c psxy
      WRITE(14,2110)
      if(ihemi.eq.1) WRITE(16,2110)
2110  FORMAT('psxy -V -O -R -JX -B -K -M ',
     &       'outline.GMT -W10 -P >> test.ps')
      WRITE(14,2111)
      WRITE(14,2112)
      WRITE(14,2113)
      if(ihemi.eq.1) then
        WRITE(16,2111)
        WRITE(16,2112)
        WRITE(16,2113)
      endif
2111  FORMAT('psxy -V -O -R -JX -B -K -M ',
     &       'lake2.GMT -W5/0/255/0 -SC0.025 -P >> test.ps')
2112  FORMAT('psxy -V -O -R -JX -B -K -M ',
     &       'lake1.GMT -W5/0/255/0 -SC0.025 -P >> test.ps')
2113  FORMAT('psxy -V -O -R -JX -B -K -M ',
     &       'points.GMT -W5/0/0/0 -SA0.025 -P >> test.ps')
c psscale
      WRITE(JUNK1,*) LABEL
      CALL STRIP(80,JUNK1,N1)
      WRITE(14,2120)
2120  FORMAT('psscale -E -V -Ctopo.cpt -D6/4/3.3/0.3',
     &       ' -O -K -L >> test.ps')
      if(ihemi.eq.1) WRITE(16,2220)
2220  FORMAT('psscale -E -V -Ctopo.cpt -D6.1/3/6/0.3',
     &       ' -O -K -L >> test.ps')
C pstext's
      WRITE(14,1040)
      WRITE(14,1050) 6.5,LABEL
      WRITE(14,1060)
      WRITE(14,1041)
      WRITE(14,1051) 6.2,HED
      WRITE(14,1061)
      WRITE(14,1070)
      if(ihemi.eq.1) then
        WRITE(16,1040)
        WRITE(16,1050) 6.5,LABEL
        WRITE(16,1060)
        WRITE(16,1041)
        WRITE(16,1051) 6.2,HED
        WRITE(16,1061)
        WRITE(16,1070)
      endif
C topo.cpt
      ZMIN=1E30
      ZMAX=-ZMIN
      IF(IPLOT.GE.0 .AND. IPLOT.LE.12) THEN
        DO I=1,NUMNP
          ZMIN=MIN(ZMIN,ZZ(I))
          ZMAX=MAX(ZMAX,ZZ(I))
        ENDDO
      ELSEIF(IPLOT.GE.17 .OR. IPLOT.LE.20) THEN
        DO I=1,NUMNP
          ZMIN=MIN(ZMIN,ZZ(I))
          ZMAX=MAX(ZMAX,ZZ(I))
        ENDDO
      ELSEIF(IPLOT.EQ.13) THEN
        DO I=1,NUMEL
          ZMIN=MIN(ZMIN,ACON(I))
          ZMAX=MAX(ZMAX,ACON(I))
        ENDDO
      ELSEIF(IPLOT.EQ.14) THEN
        DO I=1,NUMEL
          ZMIN=MIN(ZMIN,CONST(I))
          ZMAX=MAX(ZMAX,CONST(I))
        ENDDO
      ELSEIF(IPLOT.EQ.15) THEN
        ZMIN=VMIN
        ZMAX=VMAX
      ELSEIF(IPLOT.EQ.16) THEN
        ZMIN=VMIN
        ZMAX=VMAX
      ENDIF       
      COLMIN=250
      COLMAX=50
      IBAND=1
      IF(IBAND.EQ.0) THEN
c following for color bands
c        ZMIN=RMIN
c        ZMAX=RMAX
        ZMIN=-1500.
        ZMAX=2500.
C
        NLEV=10
c        NLEV=NINT((ZMAX-ZMIN)/100.)
        COLSTEP=(COLMIN-COLMAX)/(NLEV-1)
c
        ZSTEP=(ZMAX-ZMIN)/(NLEV)
        DO I=1,NLEV
          ZVAL=ZMIN+(I-1)*ZSTEP
          IF(ZVAL.GE.0.) THEN
            CB=0
          ELSE
            CB=250
          ENDIF
          CG=COLMIN-(I-1)*COLSTEP
          CR=COLMAX+(I-1)*COLSTEP
          WRITE(13,1130) ZVAL,NINT(CR),NINT(CG),NINT(CB),
     &                   ZVAL+ZSTEP,NINT(CR),NINT(CG),NINT(CB)
          WRITE(*,1130) ZVAL,NINT(CR),NINT(CG),NINT(CB),
     &                   ZVAL+ZSTEP,NINT(CR),NINT(CG),NINT(CB)
        ENDDO
      ELSE
c following for graded colors
        ICOL1=COLMIN
        ICOL2=COLMAX
c8888888888888888888888888888888888
        zmin=zmin-1.
        zmax=zmax+1.
c8888888888888888888888888888888888
        WRITE(13,1130) ZMIN,ICOL1,ICOL1,ICOL1,ZMAX,ICOL2,ICOL2,ICOL2
        WRITE(*,*) ZMIN,ICOL1,ICOL1,ICOL1,ZMAX,ICOL2,ICOL2,ICOL2
      ENDIF
1130  FORMAT(1X,F10.3,3I4,F10.3,3I4)
200   FORMAT(I6,I4,1P2E12.5,0PF10.2,F7.2,F9.2,F10.3,F10.1,F10.5,2F10.5,
     &       I5,F10.3,F10.3,1PE10.3)
300   FORMAT(5I6,1P2E17.10)
310   FORMAT(2I6,E13.6)
999   continue
      END
c***************************************************
      SUBROUTINE ISWAP(I,J)
      ITEMP=I
      I=J
      J=ITEMP
      END
C******************************
      SUBROUTINE STRIP(N,JUNK,IC)
      CHARACTER*80 JUNK,JUNKN
      JUNKN=' '
      IC=0
      DO I=1,N
        IF(JUNK(I:I).NE.' ') THEN
          IC=IC+1
          JUNKN(IC:IC)=JUNK(I:I)
        ENDIF
      ENDDO
      JUNK=JUNKN
c      PRINT *,IC,JUNKN
      END
C*****************************************
      SUBROUTINE SCALPL(IPLOT,RLEVSP,RMAX,RMIN,LABEL)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER LABEL*26
      IF(IPLOT.EQ.0 .OR. IPLOT.EQ.8) THEN
        RLEVSP=500.
        RMAX=7500.
        RMIN=0.0
        IF(IPLOT.EQ.0) LABEL=' SURFACE '
        IF(IPLOT.EQ.8) LABEL=' PRESENT SURFACE '
      ELSEIF(IPLOT.EQ.1) THEN
        RLEVSP=100.
        RMAX=6000.
        RMIN=-1500.
        LABEL=' BEDROCK '
      ELSEIF(IPLOT.EQ.2) THEN
        RLEVSP=.2
        RMAX=2.0
        RMIN=-2.0
        LABEL=' MASS BALANCE '
      ELSEIF(IPLOT.EQ.3) THEN
        RLEVSP=.5
        RMAX=1.0
        RMIN=0.0
        LABEL=' PERCENT '
      ELSEIF(IPLOT.EQ.4 ) THEN
        RLEVSP=1.
        RMAX=7.
        RMIN=0.
        LABEL=' FLOW CONSTANT '
      ELSEIF(IPLOT.EQ.5) THEN
        RLEVSP=150.
        RMAX=2100.
        RMIN=-1200.
        LABEL=' DIFFERENCE '
      ELSEIF(IPLOT.EQ.6) THEN
        RLEVSP=500.
        RMAX=6000.
        RMIN=0.0
        LABEL=' THICKNESS '
      ELSEIF(IPLOT.EQ.7) THEN
        RLEVSP=1.
        RMAX=1.5
        RMIN=1.0000
        LABEL=' MARGIN '
      ELSEIF(IPLOT.EQ.9) THEN
        RLEVSP=1.
        RMAX=-110.
        RMIN=-121.
        LABEL=' ZONE '
      ELSEIF(IPLOT.EQ.10) THEN
        RLEVSP=.005
        RMAX=.021
        RMIN=-.004
        LABEL=' SLIDING CONSTANT '
      ELSEIF(IPLOT.EQ.11) THEN
        RLEVSP=5.
        RMAX=0.0
        RMIN=-50.
        LABEL=' TEMPERATURE '
      ELSEIF(IPLOT.EQ.12) THEN
        RLEVSP=.1
        RMAX=3.0
        RMIN=0.
        LABEL=' A FUDGE '
      ELSEIF(IPLOT.EQ.13 ) THEN
        RLEVSP=1.
        RMAX=7.
        RMIN=0.
        LABEL=' ELEMENT FLOW CONSTANT '
      ELSEIF(IPLOT.EQ.14 ) THEN
        RLEVSP=100.
        RMAX=1000.
        RMIN=0.
        LABEL=' NONLINEAR CONSTANT (*1E6)'
      ELSEIF(IPLOT.EQ.15 ) THEN
        RLEVSP=10.
        RMAX=10000.
        RMIN=0.
        LABEL=' VELOCITY '
      ELSEIF(IPLOT.EQ.16 ) THEN
        RLEVSP=1.
        RMAX=100.
        RMIN=0.
        LABEL=' SURFACE GRADIENT (M/KM)'
      ELSEIF(IPLOT.EQ.17 ) THEN
        RLEVSP=1.
        RMAX=20.
        RMIN=0.
        LABEL=' BASAL WATER THICKNESS (M)'
      ELSEIF(IPLOT.EQ.18 ) THEN
        RLEVSP=1.
        RMAX=5.
        RMIN=-5.
        LABEL=' BASAL MELT RATE (MM/YR)'
      ELSEIF(IPLOT.EQ.19) THEN
        RLEVSP=500.
        RMAX=6000.
        RMIN=-2000.
        LABEL=' ICE SURFACE AND BEDROCK '
      ELSEIF(IPLOT.EQ.20) THEN
        RLEVSP=500.
        RMAX=6000.
        RMIN=-2000.
        LABEL=' BASAL WATER POTENTIAL '
      ELSE
        PRINT *,'PROBLEMS...'
        LABEL=' PROBLEMS '
      ENDIF
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
      rlat=y/rkmpdeg
      rlong=-127.5+x/rkmpdeg/cos(rlat*radpdeg)
      if(rlong.lt.0) rlong=360+rlong
      END
