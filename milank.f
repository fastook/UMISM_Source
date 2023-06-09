      subroutine milank(TIME,x,y,elev1,slope1,shape1,TNSL,tmean,accum)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nyear1=-100,nyear2=100,nlat=19,nmonth=12,amplify=0.1)
      integer insol(nmonth,nlat,nyear1:nyear2)
      integer lat(nlat)
      real rmil(nmonth),rmil1(nmonth),rmil2(nmonth)
      real rmilba(nmonth)
      character*80 line
      common /lapse/ acom,hmax,windir(2),xpole,ypole,abl,tnslbase
      DIMENSION QI(12),QS(12),TTT(12)                                   
      DATA QI /960.d0,1036.d0,1200.d0,825.,330.d0,90.d0,150.d0,
     &         600.d0,1200.d0,1020.d0,    
     &         930.d0,850.d0/                                               
      DATA QS /-0.667d0,4.6d0,11.667d0,9.167d0,3.667d0,
     *         1.d0,1.667d0,6.667d0,12.d0,6.333d0,  
     &         0.333d0,-3.333d0/                                            
      data ipass /0/,pi /3.1415927/
      data www /  19.1390686e0     /                                     
      data xxx / 0.922791243e0     /                                     
      data zzz /-0.738900483e0     /                                     
      save ipass,insol,iy,lat
C ... CALCULATE LATITUDE                                                    
C     CALL SETRIG                                                       
c     CALL RECPOL(X-XPOLE,Y-YPOLE,RLAT,RLONG)  
      CALL RECPOL(X,Y,RLAT,RLONG)  
      if(ipass.eq.0) then
c ... done once first call ...
        ipass=1
c        open(23,file='milank.dat')
        do iyear=nyear2,nyear1,-1
          read(23,'(a)') line
          do ilat=1,nlat
            read(23,*) iy,lat(ilat),
     &                (insol(im,ilat,iyear),im=1,nmonth)
          enddo
        enddo
      endif
      ryear=TIME/1000.
      if(ryear.lt.nyear1) ryear=nyear1
      if(ryear.gt.nyear2) ryear=nyear2
c ..............................................
C ... CALCULATE LATITUDE (w/ pole shift .......                                                    
C     CALL SETRIG                                                       
      CALL RECPOL(X-XPOLE,Y-YPOLE,RLAT,RLONG)  
c     CALL RECPOL(X,Y,RLAT,RLONG)  
c ... elevation (km)                                                        
      elev=elev1/1000.d0
c ... slope (m/km)                                                          
      slope=slope1*1000.d0
c ... shape (m/km/km)                                                       
      shape=shape1*1000.d0*1000.d0
      shape=0.d0
      if(rlat.gt.90) rlat=90
      if(rlat.lt.-90) rlat=-90
      llat=int(rlat/10)*10
      ilat=1+(90-llat)/10
      rl=1+(90-rlat)/10
      ilat1=rl
      ilat2=rl+1
      if(.false.) then
        print *,lat(ilat1),rlat,lat(ilat2)
        print *,ilat1,rl,ilat2
      endif
      if(ilat1.lt.1) ilat1=1
      if(ilat2.gt.nlat) ilat2=nlat
      iyear1=int(abs(ryear))
      iyear2=int(abs(ryear)+1)
      if(ryear.lt.0) then
        itemp=iyear1
        iyear1=-iyear2
        iyear2=-itemp
      endif
      if(iyear1.lt.nyear1) iyear1=nyear1
      if(iyear2.lt.nyear1) iyear2=nyear1
      if(iyear1.gt.nyear2) iyear1=nyear2
      if(iyear2.gt.nyear2) iyear2=nyear2
      if(.false.) then
        print *,iyear1,ryear,iyear2
      endif
c ... interpolate latitude ...................
      fractl=(lat(ilat1)-rlat)/10
      do im=1,nmonth
        rmil1(im)=insol(im,ilat1,iyear1)+
     &    fractl*(insol(im,ilat2,iyear1)-insol(im,ilat1,iyear1))
        rmil2(im)=insol(im,ilat1,iyear2)+
     &    fractl*(insol(im,ilat2,iyear2)-insol(im,ilat1,iyear2))
        rmilba(im)=insol(im,ilat1,0)+
     &    fractl*(insol(im,ilat2,0)-insol(im,ilat1,0))
      enddo
      if(.false.) then
        print *,lat(ilat1),rlat,lat(ilat2)
        print *,ilat1,rl,ilat2
        print *,fractl
        print 100,lat(ilat1),(insol(im,ilat1,iyear1),im=1,4)
        print 101,rlat,(rmil1(im),im=1,4)
        print 100,lat(ilat2),(insol(im,ilat2,iyear1),im=1,4)
        print *
        print 100,lat(ilat1),(insol(im,ilat1,iyear2),im=1,4)
        print 101,rlat,(rmil2(im),im=1,4)
        print 100,lat(ilat2),(insol(im,ilat2,iyear2),im=1,4)
100     format(1x,i6,12i6)
101     format(1x,f6.1,12f6.1)
      endif
c ............................................
c ... interpolate year ........................................
      fracty=(ryear-iyear1)
      do im=1,nmonth
        rmil(im)=rmil1(im)+fracty*(rmil2(im)-rmil1(im))
      enddo
      if(.false.) then
        print *,iyear1,ryear,iyear2
        print *,fracty
        print 200,iyear1,(rmil1(im),im=1,4)
        print 201,ryear,(rmil(im),im=1,4)
        print 200,iyear2,(rmil2(im),im=1,4)
200     format(1x,i6,12f6.1)
201     format(1x,f6.1,12f6.1)
c        pause
      endif
c ..............................................................
c ... means and standard deviation for year
      rmeanno=average(nmonth,rmil)
      rmeanba=average(nmonth,rmilba)
c     if(dispon) print *,rlat,ryear,rmeanba,rmeanno,rmeanno-rmeanba
      rplusno=stdev(nmonth,rmil)
      rplusba=stdev(nmonth,rmilba)
      rmean=rmeanno-rmeanba
      rplus=rplusno-rplusba
c ... rmean is EXCESS insolation
c ... rmeanno is PAST insolation
c ... rmeanba is PRESENT insolation
c ... INCREASE PAST INSOLATION
c ... (Only increase colder)
      if(rmean.lt.0) then      
        rmean=rmean*1
        rmeanno=rmeanba+rmean
      endif
c .....................................
c     if(dispon) print *,rlat,ryear,rplusba,rplus,rplus-rplusba
      tmeanno=stepbol(rmeanno)
      tplusno=stepbol(rmeanno+rplusno*amplify)
      tmeanba=stepbol(rmeanba)
      tplusba=stepbol(rmeanba+rplusba*amplify)
      if(dispon) print *,'lat,yr:',rlat,ryear
c     if(dispon) print *,rmil
      if(dispon) print *,'insolno',rmeanno,rplusno,rmeanno+rplusno
      if(dispon) print *,'insolba',rmeanba,rplusba,rmeanba+rplusba
      if(dispon) print *,'tempsno',tmeanno,tplusno-tmeanno,tplusno
      if(dispon) print *,'tempsba',tmeanba,tplusba-tmeanba,tplusba
c     pause
      tplusno=tplusno-tmeanno
      tplusba=tplusba-tmeanba
      rmean=rmeanno
      rplus=rplusno
      tmean=tmeanno
      tplus=tplusno
c ... rest is standard climatology ..............................
c ... lapse rate with elevation .................................
      tmean=tmean+acom*elev+TNSL+TNSLBASE
      if(dispon) print *,'surface mean annual air temp=',tmean,tplus
c ... calculate mean annual temp of free atmosphere-isothemal layer ...         
      tf=0.67d0*(tmean+273.0d0)+88.9d0
c     print *,'temp free at-isothrmal layer=',tf
c ... calculate saturation vapor pressure .............................
      term1=-9.09718d0*(273.16d0/tf-1.0d0)
      term2=-3.56654d0*log10(273.16d0/tf)
      term3=0.876793d0*(1.0d0-tf/273.16d0)+0.785835d0
      expon=term1+term2+term3
      es=10.d0**expon
c     print *,'saturation vapor pressure=',es
c ... calculate accumulation rate (m/yr) ..............................
      term1=www*es
      term2=xxx*slope
      term3=zzz
      term4=-15.276d0*shape
      acc=max(1e-9,term1+term2+term3+term4)
      if(dispon) print *,'accumulation rate=',acc
C ... CALCULATE ABLATION THE OLD WAY ...
c      if(.false.) then
        QY=0.d0
        DO I=1,12                                                      
          QY=QY+QI(I)-QS(I)*RLAT                                          
        ENDDO                                                          
        QY=QY/12.d0
        PDDM=0.d0
        DO I=1,12                                                      
C         TTT(I)=tmean+0.021d0*((QI(I)+QS(I)*RLAT)-QY)+8.954d0
          TTT(I)=tmean+0.021d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C         TTT(I)=tmean+0.018d0*((QI(I)-QS(I)*RLAT)-QY)+0.0d0
C         WRITE(19,*) I,TTT(I),QI(I)-QS(I)*RLAT,QY                        
          IF(TTT(I).GT.0.0) PDDM=PDDM+30.d0*TTT(I)                            
c          if(dispon) print *,i,TTT(I),pddm,max(0.0,TTT(i)*30*0.6)
        ENDDO
c      else
c ..... sum up positive degree days for ablation calc..............
        pdd=0
        do i=0,364,30
          temp=tmean+tplus*sin(2*pi*i/364.)
          if(temp.gt.0.) pdd=pdd+temp*30
c          if(dispon) print *,i,temp,pdd,max(0.0,temp*30*0.6)
        enddo
c      endif
      abl=.6d0*pdd
      if(dispon) print *,'ablation rate=',abl,pdd,PDDM
c      print *,tmean,pdd,pddnew
c ... calculate net accumulation ......................................
      accnet=acc-abl
c     print *,'net accumulation/ablation=',accnet
c ... put in m/yr .....................................................
      accum=accnet*.01d0
      abl=abl*.01d0
c      if(dispon) pause
      end
c--------------------------------------------------
      real function stepbol(rmean)
      data sigma /5.67e-5/
      isign=sign(1.,rmean)
c ... convert langleys/day to watt/m**2, multiply by 0.4843 .....
      tmean=abs(rmean)*0.4843
c ... convert w/m**2 to erg/cm**2/s, multipy by 1000. ...........
      tmean=tmean*1000.
c ... convert to temperature e=sigma*t**4, ......................
c ............  sigma=5.67e-5 erg/cm**2/deg**4/s ................
c ... tplus is temp with mean+stdev flux ........................
      albedo=0.12 ! dirt
      albedo=0.84 ! dry snow
      albedo=0.0 ! perfect black body
c ... tmean is temp with just mean ..............................
      tmean=(tmean*(1-albedo)/sigma)**0.25
c ... offset because stephan-boltzman is cold (greenhouse??) ....
c      tmean=tmean+25.0-7.65
c ... convert to deg c ..........................................
      tmean=tmean-273.16
      stepbol=isign*tmean
      end
c------------------------------------------------
      real function average(n,r)
      real r(n)
      sum=0.0
      do i=1,n
        sum=sum+r(i)
      enddo
      average=sum/n
      end
c------------------------------------------------
      real function stdev(n,r)
      real r(n)
      rmean=average(n,r)
      sum=0.0
      do i=1,n
        sum=sum+(r(i)-rmean)**2
      enddo
      sum=sum/n
      stdev=sqrt(sum)
      end
