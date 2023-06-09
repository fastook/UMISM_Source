      implicit real*8(a-h,o-z)
      character label*20
      ebase=5000
      if(.false.) then
        ts=0
        elev=0
1       print *,'input ts,elev',ts,elev
          read(*,*,end=99) ts,elev
          acc=accmars(ts,elev,a1,a2)
          print *,ts,acc,a1,a2
        goto 1
      elseif(.false.) then
        rlapse=5.
        do rlapse=0.,5.,1.0
          open(11,file='plot.xyz')
          open(12,file='plot.data')
          do ts0=-140.,-40.,5.
            do elev=-2000.,20000.,100.
              ts=ts0-rlapse*(elev-ebase)/1000
              acc=accmars(ts,elev,a1,a2) ! m/yr
              acc=acc*1000         ! mm/yr
c             if(acc.lt.0) acc=0
              write(11,*) ts0,elev,acc
              write(12,*) acc,elev
            enddo
            write(12,*) -99999.,0.
            write(12,*) 'temp:',ts0
          enddo
          close(11)
          close(12)
c         call system('plotbps.e plot')
c         call system('plotxyz.e plot.xyz')
          write(label,*) rlapse
          call system('mb-mars.gmt '//label)
        enddo
      else
        rlapse=3.
        tbase=-80
        open(11,file='plota-m.data')
        open(12,file='plot-m.data')
        do elev=-2000.,20000.,100.
          ts=tbase-rlapse*(elev-ebase)/1000
          acc=accmars(ts,elev,a1,a2)
          if(.true. .or. acc.gt.0) then
            write(11,*) 1000*acc,elev
            write(12,*) 1000*acc,elev
            write(12,*) 0*acc,elev
            write(12,*) 1000*acc,elev
c           write(*,*) 'acc,ratio acc/abl',acc,a1/a2
c           print *,elev,ts,1000*acc
          endif
        enddo
        write(11,*) -99999.,0.
        write(11,*) rlapse,tbase
        write(12,*) -99999.,0.
        write(12,*) rlapse,tbase
        do elev=-2000.,20000.,100.
          ts=tbase-rlapse*(elev-ebase)/1000
          acc=accmars(ts,elev,a1,a2)
          if(.true. .or. acc.gt.0) then
            write(12,*) 1000*a1,elev
          endif
        enddo
        write(12,*) -99999.,0.
        write(12,*) 'acc'
        do elev=-2000.,20000.,100.
          ts=tbase-rlapse*(elev-ebase)/1000
          acc=accmars(ts,elev,a1,a2)
          if(.true. .or. acc.gt.0) then
            write(12,*) 1000*a2,elev
          endif
        enddo
        write(12,*) -99999.,0.
        write(12,*) 'abl'
        close(11)
        close(12)
        call system('plotbps.e plot-m')
        call system('plotbps.e plota-m')
      endif
99    end
c------------------------------------------------
      function satvap(tf)
      implicit real*8(a-h,o-z)
      TERM1=-9.09718d0*(273.16d0/TF-1.0d0)
      TERM2=-3.56654d0*LOG10(273.16d0/TF)
      TERM3=0.876793d0*(1.0d0-TF/273.16d0)+0.785835d0
      EXPON=TERM1+TERM2+TERM3
      satvap=10.d0**EXPON
      end
c------------------------------------------------
      function accmars(ts0,elev0,acc,abl)
c ... this is the newest one ...
      implicit real*8(a-h,o-z)
      DATA WWW /  19.1390686d0     /
c     DATA WWW /  1.91390686d-0     /
      DATA XXX / 0.922791243d0     /
      DATA ZZZ /-0.738900483d0     /
      elev=elev0+3000
      abloffset=3*0.1
      toffset=60
      toffset=5
      ts=ts0+toffset
      tmean=214d0
      pmean=5.6d0
      roverm=192d0
      grav=3.74
      p1=5.6d0
      z1= 0
      tk=ts+273d0
      TTT=(TK+tmean)*0.5d0
      term=-grav*(elev-z1)/(roverm*TTT)
      pres=p1*exp(term)
c     print *,elev,p1/pres
      
C ... CALCULATE MEAN ANNUAL TEMP OF FREE ATMOSPHERE-ISOTHERMAL LAYER
      TF=0.67d0*(TS+273.0d0)+88.9d0
c     TF=(TS+273.0d0)+40
C ... CALCULATE SATURATION VAPOR PRESSURE
      es=satvap(tf)
      es2=satvap(ts+273.d0)
c     es=satvap(tk+30)
c     es2=satvap(tk+30)
      print *,elev,ts+273,es2,tf,es,tf-(ts+273)
C ... CALCULATE ABLATION
      aaa000=www
      windf=(30000-elev)/30000
      windf=1
      factor=5
c     pres=p1
c     print *,elev,ts,windf,p1/pres,es,es2
c     print *,elev,ts,p1/pres,es,es2
      ABL=factor*windf*aaa000*es2*p1/pres
      ABL=factor*WWW*windf*ES2*p1/pres
c     ABL=ABL+(5e-6*elev)**2
c ... add 0.3 cm/yr (3 mm/yr) to ablation
      ABL=ABL+abloffset
c     ABL=ABL*wind(elev)
C ... CALCULATE ACCUMULATION RATE (M/YR)
      ACC=WWW*windf*ES
C ... CALCULATE NET ACCUMULATION
      acc=acc*0.01d0
      abl=abl*0.01d0
c     print *,elev,ts,p1,pres,acc,abl,acc-abl
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
      if(accmars.gt.0) then
        accmars=10*accmars
      else
        accmars=0.1*accmars
      endif
c     abl=acc-accnet
      return
      PRINT *,'-----------------------------'
      PRINT *,'SURFACE MEAN ANNUAL AIR TEMP=',TK
      PRINT *,'Elevation                   =',elev
      PRINT *,'Pressure                    =',pres
      PRINT *,'Wind                        =',windf
      PRINT *,'TEMP FREE AT-ISOTHRMAL LAYER=',TF
      PRINT *,'SATURATION VAPOR PRESSURE   =',real(ES),real(es2)
      PRINT *,'ACCUMULATION                =',ACC*.01
      PRINT *,'ABLATION                    =',-ABL*.01
      PRINT *,'NET                         =',-ACCNET*.01
      end
      function wind(elev)
      implicit real*8(a-h,o-z)
      elevr=elev-10000
      if(elevr.gt.0) then
        wind=1+sign(1.d0,elevr)*(1e-2*elevr)**1
      else
        wind=1
      endif
      end 
