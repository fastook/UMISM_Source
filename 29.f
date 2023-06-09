      subroutine lookup29(rlapse,rlat1,rlong1,elev1,psurf1,
     &                    tnsl,tsurf,acc,abl)
      real*8 rlapse,rlat1,rlong1,elev1,psurf1,
     &       tnsl,tsurf,acc,abl,WARMING
      parameter(nlong=192,nlat=94,nmon=12)
      dimension rlong(nlong), rlat(nlat)
      dimension pmean(nlong,nlat,nmon)
      dimension ptot(nlong,nlat)
      dimension snow(nlong,nlat)
      dimension rain(nlong,nlat)
      dimension pdd(nlong,nlat)
      dimension rmb(nlong,nlat)
      dimension tmean(nlong,nlat,nmon)
      dimension tannual(nlong,nlat)
      dimension nannual(nlong,nlat)
      dimension topo(nlong,nlat)
      dimension mday(12),dist(4),value(4),xxx(4),yyy(4)
      data mday /31,28,31,30,31,30,31,31,30,31,30,31/
      character*8 cmon(12)
      data cmon /'january','february','march','april','may','june',
     &           'july','august','sepember','october','november',
     &           'december'/
      data ipass /0/
      data zero /0.0/
      logical quiet,close
      parameter(close=.false.,quiet=.true.)
      data fudge /-5./
      data reject /-999./
      data correct /0./
      save rlong,rlat,pmean,tmean,topo,ipass
c ... read rlong,rlat,pmean,tmean topo .....................
      if(ipass.eq.0) then
        print *,' reading ncep2 climate data '
        open(12,file='../ext-data/lookup26.bin',form='unformatted')
        read(12) nlongt,nlatt,nmont
        if(nlongt.ne.nlong.or.nlatt.ne.nlat.or.nmont.ne.nmon) then
          print *,' error, stopping '
          stop
        endif
        read(12) rlong,rlat,pmean,tmean,topo
        close(12)
        ipass=1
      endif
      elev=real((elev1-psurf1)*0.001)
      toffset=real(tnsl)+correct
      if(close) then
        call closest(nlong,nlat,rlong,rlat,
     &             real(rlong1),real(rlat1),kklong,kklat)
      else
        call bracket(nlong,nlat,rlong,rlat,
     &             real(rlong1),real(rlat1),kklong,kklat,dist)
        xxx(1)=rlong(kklong  )
        xxx(2)=rlong(kklong+1)
        xxx(3)=rlong(kklong+1)
        xxx(4)=rlong(kklong  )
        yyy(1)=rlat(kklat  )
        yyy(2)=rlat(kklat  )
        yyy(3)=rlat(kklat+1)
        yyy(4)=rlat(kklat+1)
      endif
c      print *,real(rlong1),real(rlat1),rlong(kklong),rlat(kklat)
      if(quiet) then
        lat1=kklat
        long1=kklong
        if(close) then
          lat2=lat1
          long2=long1
        else
          lat2=lat1+1
          long2=long1+1
        endif
      else
        lat1=1
        lat2=nlat
        long1=1
        long2=nlong
      endif
c ... zero out arrays .................................
      do i=lat1,lat2
        do j=long1,long2
          ptot(j,i)=0.0
          pdd(j,i)=0.0
          snow(j,i)=0.0
          rain(j,i)=0.0
          rmb(j,i)=0.0
          tannual(j,i)=0.0
          nannual(j,i)=0
        enddo
      enddo
c ... calculate tannual, mean annual temperature (screen10##)
      do imon=1,nmon
        do i=lat1,lat2
          do j=long1,long2
            ediff=elev+(psurf1*0.001-max(topo(j,i),0.0))*0
            tempji=tmean(j,i,imon)+toffset
            tannual(j,i)=tannual(j,i)+tempji+
     &                   rlapse*ediff
            nannual(j,i)=nannual(j,i)+1
          enddo
        enddo
      enddo
      do i=lat1,lat2
        do j=long1,long2
          tannual(j,i)=tannual(j,i)/nannual(j,i)
        enddo
      enddo
c ... calculate pdd .......................................
      do imon=1,nmon
        do i=lat1,lat2
          do j=long1,long2
            ediff=elev+(psurf1*0.001-max(topo(j,i),0.0))*0
            WRM=WARMING(dble(TNSL+rlapse*ediff))
c           WRM=1.
            temploc=tmean(j,i,imon)+rlapse*ediff+toffset
            if(temploc.gt.(zero+fudge)) then
              pdd(j,i)=pdd(j,i)+
     &            (temploc-(zero+fudge))*mday(imon)
              rain(j,i)=rain(j,i)+
     &                  WRM*pmean(j,i,imon)*1
            else
              snow(j,i)=snow(j,i)+
     &                  WRM*pmean(j,i,imon)*1
            endif
          enddo
        enddo
      enddo
      do i=lat1,lat2
        do j=long1,long2
          if(pdd(j,i).ne.0.0) then
            pdd(j,i)=pdd(j,i)*0.6e-3
          else
            pdd(j,i)=reject
          endif
          if(pdd(j,i).eq.reject) then
            rmb(j,i)=snow(j,i)
          else
            rmb(j,i)=snow(j,i)-pdd(j,i)
          endif
        enddo
      enddo
      if(tannual(kklong,kklat).ne.reject) then
        if(close) then
          tsurf=tannual(kklong,kklat)-zero
        else
          value(1)=tannual(kklong  ,kklat  )
          value(2)=tannual(kklong+1,kklat  )
          value(3)=tannual(kklong+1,kklat+1)
          value(4)=tannual(kklong  ,kklat+1)
          tsurf=wmean(dist,value)
          call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
          tsurf=val
        endif
      else
        tsurf=reject
      endif
      if(close) then
        acc=rmb(kklong,kklat)
      else
        value(1)=rmb(kklong  ,kklat  )
        value(2)=rmb(kklong+1,kklat  )
        value(3)=rmb(kklong+1,kklat+1)
        value(4)=rmb(kklong  ,kklat+1)
c        acc=wmean(dist,value)
        call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
        acc=val
c ... EXPERIMENTAL arbitrary negative offset 0.1
c       acc=acc-0.1
      endif
      if(.false.) then
        print *,'tannual:',tannual(kklong,kklat)
        do i=1,12
          print *,'tmon:',i,tmean(kklong,kklat,i),pmean(kklong,kklat,i)
        enddo
        print *,'rm,sn,rn,pd    ',rmb(kklong,kklat),
     &                      snow(kklong,kklat),
     &                      rain(kklong,kklat),
     &                      pdd(kklong,kklat)
      endif
      if(close) then
        abl=pdd(kklong,kklat)
      else
        value(1)=pdd(kklong  ,kklat  )
        value(2)=pdd(kklong+1,kklat  )
        value(3)=pdd(kklong+1,kklat+1)
        value(4)=pdd(kklong  ,kklat+1)
        abl=wmean(dist,value)
        call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
        abl=val
      endif
      end
