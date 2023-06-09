      implicit real*8(a-h,o-z)
      if(.false.) then
        ts=0
        elev=0
1       print *,'input ts,elev',ts,elev
          read(*,*,end=99) ts,elev
          acc=accearth(ts,elev,a1,a2)
          print *,ts,acc
        goto 1
      elseif(.false.) then
        open(11,file='plot.xyz')
        open(12,file='plot.data')
        rlapse=3.
        do ts0=-90.,0.,2.
          do elev=0.,5000.,100.
            ts=ts0-rlapse*elev/1000
            acc=accearth(ts,elev,a1,a2) ! m/yr
            acc=acc*1         ! m/yr
            if(acc.lt.0) acc=0
            write(11,*) ts0,elev,acc
            write(12,*) acc,elev
          enddo
          write(12,*) -99999.,0.
          write(12,*) 'temp:',ts
        enddo
        close(11)
        close(12)
c       call system('plotbps.e plot')
c       call system('plotxyz.e plot.xyz')
        call system('mb-earth.gmt')
      else
        open(11,file='plota-e.data')
        open(12,file='plot-e.data')
        rlapse=6.
        tbase=0
        do elev=0.,5000.,100.
          ts=tbase-rlapse*elev/1000
          acc=accearth(ts,elev,a1,a2)
          write(11,*) 1*acc,elev
          write(12,*) 1*acc,elev
        enddo
        write(11,*) -99999.,0.
        write(11,*) 'm/yr,lapse:',rlapse
        write(12,*) -99999.,0.
        write(12,*) 'm/yr,lapse:',rlapse
        do elev=0.,5000.,100.
          ts=tbase-rlapse*elev/1000
          acc=accearth(ts,elev,a1,a2)
          write(12,*) 1*a1,elev
        enddo
        write(12,*) -99999.,0.
        write(12,*) 'acc'
        do elev=1000.,5000.,100.
          ts=tbase-rlapse*elev/1000
          acc=accearth(ts,elev,a1,a2)
          write(12,*) 1*a2,elev
        enddo
        write(12,*) -99999.,0.
        write(12,*) 'abl'
        close(11)
        close(12)
        call system('plotbps.e plot-e')
c       call system('plotbps.e plota-e')
      endif
99    end
      function accmars0(ts0,elev,acc,abl)
      implicit real*8(a-h,o-z)
c     DATA WWW /  19.1390686d0     /
      DATA WWW /  1.91390686d0     /
      DATA XXX / 0.922791243d0     /
      DATA ZZZ /-0.738900483d0     /
      toffset=60
      ts=ts0+toffset
      tmean=214d0
      pmean=5.6d0
      roverm=192d0
      grav=3.73
      p1=5.6d0
      z1=5000
      tk=ts+273d0
      TTT=(TK+tmean)*0.5d0
      term=-grav*(elev-z1)/(roverm*TTT)
      pres=p1*exp(term)
C ... CALCULATE MEAN ANNUAL TEMP OF FREE ATMOSPHERE-ISOTHERMAL LAYER
      TF=0.67d0*(TS+273.0d0)+88.9d0
C ... CALCULATE SATURATION VAPOR PRESSURE
      es=satvap(tf)
c     es2=satvap(ts+273.d0)
      es2=es
C ... CALCULATE ACCUMULATION RATE (M/YR)
      TERM1=WWW*ES
      TERM2=XXX*SLOPE
      TERM3=ZZZ
      TERM4=-15.276d0*SHAPE
c **** <<< EXPERIMENTAL >>> **** turn off slope term ...
      TERM2=0
      TERM3=0
      TERM4=0
c **** <<< EXPERIMENTAL >>> ****
      ACC=max(0.d0,TERM1+TERM2+TERM3+TERM4)
C ... CALCULATE ABLATION
      factor=5
      aaa000=www
      ebase=z1
      wind=(ebase-elev)**4/ebase**4
      wind=(1-elev/ebase)**2
      print *,elev,wind,p1/pres
      wind=1
      factor=3
c     pres=1
      ABL=factor*wind*aaa000*es2*(p1/pres)
C ... CALCULATE NET ACCUMULATION
      ACCNET=ACC-ABL
      acc=acc*0.01
      abl=abl*0.01
      accmars=-ACCNET*.01d0
      accmars= ACCNET*.01d0
      return
      PRINT *,'-----------------------------'
      PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TK
      PRINT *,'Elevation                   =',elev
      PRINT *,'Pressure                    =',pres
      PRINT *,'Wind                        =',wind
      PRINT *,'TEMP FREE AT-ISOTHRMAL LAYER=',TF
      PRINT *,'SATURATION VAPOR PRESSURE   =',real(ES),real(es2)
      PRINT *,'ACCUMULATION                =',ACC*.01
      PRINT *,'ABLATION                    =',-ABL*.01
      PRINT *,'NET                         =',-ACCNET*.01
      end
      function satvap(tf)
      implicit real*8(a-h,o-z)
      TERM1=-9.09718d0*(273.16d0/TF-1.0d0)
      TERM2=-3.56654d0*LOG10(273.16d0/TF)
      TERM3=0.876793d0*(1.0d0-TF/273.16d0)+0.785835d0
      EXPON=TERM1+TERM2+TERM3
      satvap=10.d0**EXPON
      end

c------------------------------------------------
      function accmars(ts0,elev,acc,abl)
      implicit real*8(a-h,o-z)
c     DATA WWW /  19.1390686d0     /
      DATA WWW /  1.91390686d-2     /
      DATA XXX / 0.922791243d0     /
      DATA ZZZ /-0.738900483d0     /
      toffset=60
      toffset=0
      ts=ts0+toffset
      tmean=214d0
      pmean=5.6d0
      roverm=192d0
      grav=3.73
      p1=5.6d0
      z1=0d0
      tk=ts+273d0
      TTT=(TK+tmean)*0.5d0
      term=-grav*(elev-z1)/(roverm*TTT)
      pres=p1*exp(term)
C ... CALCULATE MEAN ANNUAL TEMP OF FREE ATMOSPHERE-ISOTHERMAL LAYER
      TF=0.67d0*(TS+273.0d0)+88.9d0
C ... CALCULATE SATURATION VAPOR PRESSURE
      es=satvap(tf)
      es2=satvap(ts+273.d0)
C ... CALCULATE ACCUMULATION RATE (M/YR)
      TERM1=WWW*ES
      TERM2=XXX*SLOPE
      TERM3=ZZZ
      TERM4=-15.276d0*SHAPE
c **** <<< EXPERIMENTAL >>> **** turn off slope term ...
      TERM2=0
      TERM3=0
      TERM4=0
c **** <<< EXPERIMENTAL >>> ****
      ACC=max(0.d0,TERM1+TERM2+TERM3+TERM4)
C ... CALCULATE ABLATION
      factor=10
      aaa000=www
      wind=(30000-elev)/30000
c     wind=1
      factor=10
c     pres=p1
      print *,elev,ts,wind,p1/pres
      ABL=factor*wind*aaa000*es2*p1/pres
C ... CALCULATE NET ACCUMULATION
      acc=acc*0.1d0
      abl=abl*0.1d0
      if(.false.) then
        tmp=abl
        abl=acc
        acc=tmp
      endif
c     acc=es
c     abl=es2
      ACCNET=ACC-ABL
c     if(accnet.lt.0) then
c       accnet=accnet*0.1
c     endif
      accmars= ACCNET
c     abl=acc-accnet
      return
      PRINT *,'-----------------------------'
      PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TK
      PRINT *,'Elevation                   =',elev
      PRINT *,'Pressure                    =',pres
      PRINT *,'Wind                        =',wind
      PRINT *,'TEMP FREE AT-ISOTHRMAL LAYER=',TF
      PRINT *,'SATURATION VAPOR PRESSURE   =',real(ES),real(es2)
      PRINT *,'ACCUMULATION                =',ACC*.01
      PRINT *,'ABLATION                    =',-ABL*.01
      PRINT *,'NET                         =',-ACCNET*.01
      end
      function accearth(ts0,elev0,acc,abl)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION QI(12),QS(12),TTT(12)                                   
      DATA QI /960.d0,1036.d0,1200.d0,825.,330.d0,90.d0,150.d0,
     &         600.d0,1200.d0,1020.d0,    
     &         930.d0,850.d0/                                               
      DATA QS /-0.667d0,4.6d0,11.667d0,9.167d0,3.667d0,
     *         1.d0,1.667d0,6.667d0,12.d0,6.333d0,  
     &         0.333d0,-3.333d0/                                            
       DATA AAA / -9.62376690d0     /                                     
       DATA BBB /-0.546917617d0     /                                     
       DATA CCC /  24.9793854d0     /                                     
       DATA WWW /  19.1390686d0     /                                     
       DATA XXX / 0.922791243d0     /                                     
       DATA ZZZ /-0.738900483d0     /  
      rlat=45
      SHAPE1=0.d0                                   
C ... ELEVATION (KM)                                                        
      ELEV=ELEV0/1000.d0
C ... SLOPE (M/KM)                                                          
c     SLOPE=SLOPE1*1000.d0
      SLOPE=0.d0
C ... SHAPE (M/KM/KM)                                                       
c     SHAPE=SHAPE1*1000.d0*1000.d0
      SHAPE=0.d0
      TS=ts0
C     PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TS                        
C ... CALCULATE MEAN ANNUAL TEMP OF FREE ATMOSPHERE-ISOTHEMAL LAYER         
      TF=0.67d0*(TS+273.0d0)+88.9d0
C     PRINT *,'TEMP FREE AT-ISOTHRMAL LAYER=',TF                        
C ... CALCULATE SATURATION VAPOR PRESSURE                                   
      TERM1=-9.09718d0*(273.16d0/TF-1.0d0)                                    
      TERM2=-3.56654d0*LOG10(273.16d0/TF)                                   
      TERM3=0.876793d0*(1.0d0-TF/273.16d0)+0.785835d0
      EXPON=TERM1+TERM2+TERM3                                           
      ES=10.d0**EXPON                                                     
C     PRINT *,'SATURATION VAPOR PRESSURE=',ES                           
C ... CALCULATE ACCUMULATION RATE (M/YR)                                    
      TERM1=WWW*ES                                                      
      TERM2=XXX*SLOPE 
c **** <<< EXPERIMENTAL >>> **** turn off slope term ...
      TERM2=0
c **** <<< EXPERIMENTAL >>> ****
      TERM3=ZZZ                                                         
      TERM4=-15.276d0*SHAPE                                               
      ACC=max(0.d0,TERM1+TERM2+TERM3+TERM4)
C     WRITE(19,113) TERM1,TERM2,TERM3,ACC 
113   FORMAT(4G13.6)                                                    
C     PRINT *,'ACCUMULATION=',ACC                                       
C ... CALCULATE ABLATION                                                    
      QY=0.d0
      DO I=1,12                                                      
        QY=QY+QI(I)-QS(I)*RLAT                                          
      ENDDO                                                          
      QY=QY/12.d0
      PDD=0.d0
      DO I=1,12                                                      
C       TTT(I)=TS+0.021d0*((QI(I)+QS(I)*RLAT)-QY)+8.954d0
        TTT(I)=TS+0.021d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C       TTT(I)=TS+0.018d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C       WRITE(19,*) I,TTT(I),QI(I)-QS(I)*RLAT,QY                        
        IF(TTT(I).GT.0.0) PDD=PDD+30.d0*TTT(I)                            
      ENDDO
c      print *,ts,pdd
      ABL=.6d0*PDD                                                        
C ... CALCULATE NET ACCUMULATION                                            
      ACCNET=ACC-ABL                                                    
      accearth=ACCNET*.01d0                                                  
      ACC=ACC*.01d0                                                  
      ABL=ABL*.01d0                                                  
      END                                                               
