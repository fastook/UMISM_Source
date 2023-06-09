      subroutine lookup25(rlapse,rlat1,rlong1,elev1,tnsl,
     &                  tsurf,acc,abl)
      real*8 rlapse,rlat1,rlong1,elev1,tnsl,tsurf,acc,abl,WARMING
      parameter(nlat=80,nlong=360,ntot=nlat*nlong,
     &          nmon=12,nyr=2,mmon=nmon*nyr,
     &          mtot=nlong*nlat*nmon*nyr)
      dimension temp(nlong,nlat,nmon,nyr)
      dimension tempstd(nlong,nlat,nmon,nyr)
      dimension precip(nlong,nlat,nmon,nyr)
      dimension fill(nlat,nmon,nyr)
      dimension rmb(nlong,nlat)
      dimension snow(nlong,nlat)
      dimension rain(nlong,nlat)
      dimension tmean(nlong,nlat),nmean(nlong,nlat)
      dimension pmean(nlong,nlat)
      dimension pdd(nlong,nlat)
      dimension ampli(nlong,nlat),ttmon(mmon)
      dimension mday(12)
      dimension topo(nlong,nlat)
      data mday /31,28,31,30,31,30,31,31,30,31,30,31/
      data nmean /ntot*0/
      data tmean /ntot*0.0/
      data pmean /ntot*0.0/
      data snow /ntot*0.0/
      data rain /ntot*0.0/
      data pdd /ntot*0.0/
      data temp /mtot*0.0/
      data zero /273.16/
      data ipass /0/
      data fudge /0./
      data factor /-1./
      save topo,temp,precip,tempstd,ipass
c
      rlat=rlat1
      rlong=rlong1
      elev=elev1*0.001
      toffset=tnsl
      call ll2ji(rlong,rlat,kklong,kklat)
      n4=4
      n12=12
      n2=2
      if(ipass.eq.0) then
        print *,'reading climate data'
        ipass=1
        open(49,file='../ext-data/topo.bindmp',
     &       form='unformatted')
        read(49) topo
        do i=1,nlat
          do j=1,nlong
            if(topo(j,i).lt.0.0) topo(j,i)=0.0
          enddo
        enddo
        close(49)
        open(49,file='../ext-data/temperature.bindmp',
     &       form='unformatted')
        read(49) temp
        do iyr=1,n2
          do imon=1,n12
            do i=1,nlat
              fill(i,imon,iyr)=0.0
              nfill=0
              do j=1,nlong
                if(temp(j,i,imon,iyr).ne.0.0) then
                  fill(i,imon,iyr)=fill(i,imon,iyr)+
     &                             temp(j,i,imon,iyr)-
     &                             rlapse*topo(j,i)
                  nfill=nfill+1
                endif
              enddo
              if(nfill.ne.0) then
                fill(i,imon,iyr)=fill(i,imon,iyr)/nfill
              else
                fill(i,imon,iyr)=0
              endif
            enddo
          enddo
        enddo
        do iyr=1,n2
          do imon=1,n12
            do i=1,nlat
              do j=1,nlong
                if(temp(j,i,imon,iyr).eq.0.0) then
                  temp(j,i,imon,iyr)=fill(i,imon,iyr)
                endif
              enddo
            enddo
          enddo
        enddo
        close(49)
        open(49,file='../ext-data/temp2mstd.bindmp',
     &       form='unformatted')
        read(49) tempstd
        do iyr=1,n2
          do imon=1,n12
            do i=1,nlat
              fill(i,imon,iyr)=0.0
              nfill=0
              do j=1,nlong
                if(tempstd(j,i,imon,iyr).ne.0.0) then
                  fill(i,imon,iyr)=fill(i,imon,iyr)+
     &                             tempstd(j,i,imon,iyr)
                  nfill=nfill+1
                endif
              enddo
              if(nfill.ne.0) then
                fill(i,imon,iyr)=fill(i,imon,iyr)/nfill
              else
                fill(i,imon,iyr)=0
              endif
            enddo
          enddo
        enddo
        do iyr=1,n2
          do imon=1,n12
            do i=1,nlat
              do j=1,nlong
                if(tempstd(j,i,imon,iyr).eq.0.0) then
                  tempstd(j,i,imon,iyr)=fill(i,imon,iyr)
                endif
              enddo
            enddo
          enddo
        enddo
        close(49)
        open(49,file='../ext-data/precip.bindmp',
     &       form='unformatted')
        read(49) precip
        do iyr=1,n2
          do imon=1,n12
            do i=1,nlat
              fill(i,imon,iyr)=0.0
              nfill=0
              do j=1,nlong
                if(precip(j,i,imon,iyr).ne.0.0) then
                  fill(i,imon,iyr)=fill(i,imon,iyr)+
     &                             precip(j,i,imon,iyr)
                  nfill=nfill+1
                endif
              enddo
              if(nfill.ne.0) then
                fill(i,imon,iyr)=fill(i,imon,iyr)/nfill
              else
                fill(i,imon,iyr)=0
              endif
            enddo
          enddo
        enddo
        do iyr=1,n2
          do imon=1,n12
            do i=1,nlat
              do j=1,nlong
                if(precip(j,i,imon,iyr).eq.0.0) then
                  precip(j,i,imon,iyr)=fill(i,imon,iyr)
                endif
              enddo
            enddo
          enddo
        enddo
        close(49)
      endif
      do i=kklat,kklat
        do j=kklong,kklong
          tmean(j,i)=0.0
          pmean(j,i)=0.0
          snow(j,i)=0.0
          rain(j,i)=0.0
          pdd(j,i)=0.0
          nmean(j,i)=0
        enddo
      enddo
      do iyr=1,n2
        do imon=1,n12
          do i=kklat,kklat
            do j=kklong,kklong
              if(precip(j,i,imon,iyr).ne.0.0) then
                pmean(j,i)=pmean(j,i)+precip(j,i,imon,iyr)
                nmean(j,i)=nmean(j,i)+1
              endif
            enddo
          enddo
        enddo
      enddo
      do i=kklat,kklat
        do j=kklong,kklong
          if(pmean(j,i).ne.0) then
            pmean(j,i)=0.001*pmean(j,i)/n2
          else
            pmean(j,i)=-999.
          endif
        enddo
      enddo
c
      do i=kklat,kklat
        do j=kklong,kklong
          pdd(j,i)=0.0
          snow(j,i)=0.0
          rain(j,i)=0.0
          tmean(j,i)=0.0
          nmean(j,i)=0
        enddo
      enddo
       do iyr=1,n2
        do imon=1,n12
          do i=kklat,kklat
            do j=kklong,kklong
              if(temp(j,i,imon,iyr).ne.0.0) then
c ... use the following to REDUCE preciptation as climate cools ...
                WRM=WARMING(TNSL+rlapse*(elev-topo(j,i)))
c               WRM=1.0
                temploc=temp(j,i,imon,iyr)+rlapse*(elev-topo(j,i))+
     &                  toffset
                fudge=tempstd(j,i,imon,iyr)*factor
                tmean(j,i)=tmean(j,i)+temploc
                nmean(j,i)=nmean(j,i)+1
                if(temploc.gt.(zero+fudge)) then
                  pdd(j,i)=pdd(j,i)+
     &            (temploc-(zero+fudge))*mday(imon)
                  rain(j,i)=rain(j,i)+
     &                      WRM*precip(j,i,imon,iyr)*0.001/n2
                else
                  snow(j,i)=snow(j,i)+
     &                      WRM*precip(j,i,imon,iyr)*0.001/n2
                endif
              endif
            enddo
          enddo
        enddo
      enddo
c
      if(.false.) then
        do i=kklat,kklat
          do j=kklong,kklong
            ampli(j,i)=-999.
          enddo
        enddo
        do i=kklat,kklat
          do j=kklong,kklong
            jmon=0
            do iyr=1,n2
              do imon=1,n12
                if(temp(j,i,imon,iyr).ne.0.0) then
                  temploc=temp(j,i,imon,iyr)+rlapse*(elev-topo(j,i))
                  jmon=jmon+1
                  ttmon(jmon)=temploc
                  endif
              enddo
            enddo
            if(jmon.gt.1) then
              call fourier(jmon,ttmon,cccmax,ncmax)
              ampli(j,i)=cccmax
            else
              ampli(j,i)=-999.
            endif
          enddo
        enddo
      endif
c
      do i=kklat,kklat
        do j=kklong,kklong
          if(nmean(j,i).ne.0) then
            tmean(j,i)=tmean(j,i)/nmean(j,i)
          else
            tmean(j,i)=-999.
          endif
          if(pdd(j,i).ne.0.0) then
            pdd(j,i)=pdd(j,i)*0.6e-3/n2
          else
            pdd(j,i)=-999.
          endif
          if(pmean(j,i).eq.-999.) then
            rmb(j,i)=-999.
          elseif(pdd(j,i).eq.-999.) then
            rmb(j,i)=snow(j,i)
          else
            rmb(j,i)=snow(j,i)-pdd(j,i)
          endif
          if(rain(j,i).eq.0) rain(j,i)=-999.
          if(snow(j,i).eq.0) snow(j,i)=-999.
        enddo
      enddo
      if(tmean(kklong,kklat).ne.-999.) then
        tsurf=tmean(kklong,kklat)-zero
      else
        tsurf=-999.
      endif
      acc=rmb(kklong,kklat)
      abl=pdd(kklong,kklat)
      if(tsurf.ne.-999. .and. .false.) then
        print *,' topo         (km)  ',topo(kklong,kklat),elev
        print *,' mean temp    (degC)',tmean(kklong,kklat)-zero
        print *,' mean precip  (m)   ',pmean(kklong,kklat)
        print *,' rain         (m)   ',-rain(kklong,kklat),
     &          100*rain(kklong,kklat)/pmean(kklong,kklat)
        print *,' snow         (m)   ',snow(kklong,kklat),
     &          100*snow(kklong,kklat)/pmean(kklong,kklat)
        print *,' pdd-ablation (m)   ',-pdd(kklong,kklat)
        print *,' mass balance (m)   ',rmb(kklong,kklat)
        print *,' seasonal amplitude ',ampli(kklong,kklat)
      endif
      end
C===========================================
      subroutine fourier(npts,f,cccmax,ncmax)
      parameter(m=4)
c     parameter(pi=4.*atan(1.))
      data pi /3.1415927/
      dimension f(npts),aaa(0:m),bbb(m),ccc(m),phi(m)
      rl=real(npts/2)
      sum=0.
      do i=1,npts
        sum=sum+f(i)
      enddo
      aaa(0)=sum/rl
c      print *,0,aaa(0)/2
      aaamax=-1e30
      bbbmax=-1e30
      namax=0
      nbmax=0
      fmin=1e30
      fmax=-fmin
      do i=1,npts
        fmin=min(fmin,f(i))
        fmax=max(fmax,f(i))
      enddo
      do j=1,m
        suma=0.0
        sumb=0.0
        do i=1,npts
          suma=suma+f(i)*cos(j*pi*i/rl)

          sumb=sumb+f(i)*sin(j*pi*i/rl)
        enddo
        aaa(j)=suma/rl
        bbb(j)=sumb/rl
        if(abs(aaa(j)).gt.aaamax) then
          aaamax=abs(aaa(j))
          namax=j
        endif
        if(abs(bbb(j)).gt.bbbmax) then
          bbbmax=abs(bbb(j))
          nbmax=j
        endif
c        print *,j,aaa(j),bbb(j)
      enddo
      cccmax=-1e30
      ncmax=0
      do n=1,m
        ccc(n)=sqrt(aaa(n)**2+bbb(n)**2)
        if(aaa(n).lt.0) ccc(n)=-ccc(n)
        if(abs(ccc(n)).gt.cccmax) then
          cccmax=abs(ccc(n))
          ncmax=n
        endif
        phi(n)=atan(aaa(n)/bbb(n))
      enddo
      if(.false.) then
        call grstrt(600,600)
        call window(0.,real(24+1),fmin,fmax)
        call linclr(1)
        call move(1.,f(1))
        do jj=1,npts
          call draw(real(jj),f(jj))
          call point(real(jj),f(jj))
        enddo
        call linclr(2)
        do jj=1,npts
          sum=aaa(0)/2.
            do n=namax,namax
              sum=sum+aaa(n)*cos(n*pi*jj/rl)
            enddo
            do n=nbmax,nbmax
              sum=sum+bbb(n)*sin(n*pi*jj/rl)
            enddo
            call point(real(jj),sum)
        enddo
        if(.false.) then
          call newpag
        endif
        call grstop1
      endif
      end
c===================================      
      subroutine ji2ll(j,i,rlong,rlat)
      rlong=j-180-0.5
      if(rlong.lt.0) rlong=rlong+360.
      rlat=90-i+0.5
      end
c=========================================
      subroutine ll2ji(rlong,rlat,j,i)
      j=nint(rlong+180.5)
      if(j.gt.360) j=j-360
      i=nint(89.5-rlat)
      if(j.le.0) then
        j=1
      elseif(j.gt.360) then
        j=360
      endif
      if(i.le.0) then
        i=1
      elseif(i.gt.80) then
        j=80
      endif
      end
c=============================================
      subroutine interp1(x,y,v,xtest,ytest,val)
      implicit real*8(u)
      dimension x(4),y(4),v(4)
c this is machine generated code, from macsyma. dont change anything....
c (it's about 20 times faster than inverting the matrix.................
c-----------------------------------------------------------------------
      U0=X(1)
      U1=Y(1)
      U2=U0*U1
      U3=X(2)
      U4=-(U1*U3)
      U5=U4+U2
      U6=Y(2)
      U7=X(3)
      U8=U5*U6*U7
      U9=-(U0*U1*U3)
      U10=U0*U3*U6
      U11=U1*U3
      U12=-(U0*U6)+U11
      U13=Y(3)
      U14=(U12*U7+U10+U9)*U13
      U15=-(U0*U1)
      U16=U11+U15
      U17=-(U3*U6)
      U18=U6-U1
      U19=(U18*U7+U17+U2)*U13+U16*U6
      U20=X(4)
      U21=U0*U1*U3
      U22=-(U0*U3*U6)
      U23=U3*U6+U15
      U24=U23*U7
      U25=-U3+U0
      U26=U25*U7*U13
      U27=U0*U6
      U28=-U6+U1
      U29=U3-U0
      U30=U29*U13+U28*U7+U27+U4
      U31=Y(4)
      U32=V(4)
      U33=V(3)
      U34=V(2)
      U35=U0*U1*U34
      U36=V(1)
      U37=-(U36*U3*U6)
      U38=-(U1*U34)
      U39=U36*U6
      U40=(U39+U38)*U7
      U41=-(U0*U1*U34)
      U42=U36*U3*U6
      U43=U0*U34
      U44=-(U36*U3)
      U45=U44+U43
      U46=U1*U34
      U47=-(U36*U6)
      U48=(U47+U46)*U7
      U49=-(U0*U34)
      U50=U36*U3
      U51=U50+U49
      U52=-U34+U36
      U53=U34-U36
      VAL=((((U53*U7+U25*U33+U50+U49)*U31+(U52*U13+U18*U33+U47+U46)*U20+
     . U30*U32+U45*U13+U40+U12*U33)*XTEST+(U52*U7+U29*U33+U44+U43)*U20*
     . U31+(U53*U7*U13+(U17+U2)*U33+U42+U41)*U20+(U26+U24+U22+U21)*U32+
     . U51*U7*U13+(U37+U35)*U7+(U10+U9)*U33)*YTEST+(((U53*U13+U28*U33+
     . U39+U38)*U20+U52*U7*U13+U23*U33+U37+U35)*U31+U19*U32+(U48+U42+U41
     . )*U13+U5*U6*U33)*XTEST+((U51*U13+U48+(U27+U4)*U33)*U20+U45*U7*U13
     . +(U42+U41)*U7+(U22+U21)*U33)*U31+((U40+U37+U35)*U13+U16*U6*U33)*
     . U20+(U14+U8)*U32)/((U30*U20+U26+U24+U22+U21)*U31+U19*U20+U14+U8)
c-----------------------------------------------------------------------
      end

c=============================================
      subroutine interp0(xxx,yyy,value,xtest,ytest,val)
      dimension rmat(4,4),rhs(4)
      dimension xxx(4),yyy(4),value(4)
      do i=1,4
        rhs(i)=value(i)
      enddo
      do i=1,4
        rmat(i,1)=1.0            
        rmat(i,2)=xxx(i)            
        rmat(i,3)=yyy(i)            
        rmat(i,4)=xxx(i)*yyy(i)
      enddo  

      call GAUSSJ(rmat,4,4,rhs,1,1)

      val=rhs(1)+rhs(2)*xtest+
     &         rhs(3)*ytest+rhs(4)*xtest*ytest
      end
c=============================================

      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      PARAMETER (NMAX=50)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) PAUSE 'Singular matrix.'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
c==============================================================
      subroutine lookup26(rlapse,rlat1,rlong1,elev1,psurf1,
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
c      fudge=factor
      factor=fudge
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
            tempji=tmean(j,i,imon)+toffset
            tannual(j,i)=tannual(j,i)+tempji+
     &                   rlapse*elev
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
c           if(.not.quiet) elev=topo(j,i)
            WRM=REDUCE(real(TNSL+rlapse*elev))
c           WRM=REDUCE2(real(TNSL),real(rlapse*elev))
c           WRM=1.
            temploc=tmean(j,i,imon)+rlapse*elev+toffset
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
c=========================================
      subroutine lookup27(rlapse,rlat1,rlong1,elev1,psurf1,
     &                    tnsl,tsurf,acc,abl)
      real*8 rlapse,rlat1,rlong1,elev1,psurf1,
     &       tnsl,tsurf,acc,abl,WARMING
      parameter(nlong=128,nlat=64,nmon=12)
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
      data correct /14./
      save rlong,rlat,pmean,tmean,topo,ipass
c ... read rlong,rlat,pmean,tmean topo .....................
c      fudge=factor
      factor=fudge
      if(ipass.eq.0) then
        print *,' reading ncep2 climate data '
        open(12,file='../ext-data/lookup27.bin',form='unformatted')
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
            tempji=tmean(j,i,imon)
            tannual(j,i)=tannual(j,i)+tempji+toffset+
     &                   rlapse*elev
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
c           if(.not.quiet) elev=topo(j,i)
            WRM=REDUCE(real(TNSL+rlapse*elev))
c           WRM=REDUCE2(real(TNSL),real(rlapse*elev))
c           WRM=1.
            temploc=tmean(j,i,imon)+rlapse*elev+toffset
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
        acc=wmean(dist,value)
        call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
        acc=val
c ... EXPERIMENTAL arbitrary negative offset 0.1
c        acc=acc-0.1
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
c===========================================================
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
            ediff=elev+(psurf1*0.001-max(topo(j,i),0.0))*1
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
            ediff=elev+(psurf1*0.001-max(topo(j,i),0.0))*1
c           WRM=WARMING(dble(TNSL+rlapse*ediff))
c           if(.not.quiet) elev=topo(j,i)
            WRM=REDUCE(real(TNSL+rlapse*ediff))
c           WRM=REDUCE2(real(TNSL),real(rlapse*ediff))
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
c=========================================
      subroutine lookup30(rlapse,rlat1,rlong1,elev1,psurf1,
     &                    tnsl,tsurf,acc,abl)
      real*8 rlapse,rlat1,rlong1,elev1,psurf1,
     &       tnsl,tsurf,acc,abl,WARMING
      parameter(nlong=128,nlat=64,nmon=12)
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
      data correct /14./
      save rlong,rlat,pmean,tmean,topo,ipass
c ... read rlong,rlat,pmean,tmean topo .....................
      if(ipass.eq.0) then
        print *,' reading ncep2 climate data '
        open(12,file='../ext-data/lookup27.bin',form='unformatted')
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
            ediff=elev+(psurf1*0.001-max(topo(j,i),0.0))*1
c           print*,psurf1*0.001,max(topo(j,i),0.0)
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
            ediff=elev+(psurf1*0.001-max(topo(j,i),0.0))*1
c           WRM=WARMING(dble(TNSL+rlapse*ediff))
c           if(.not.quiet) elev=topo(j,i)
            WRM=REDUCE(real(TNSL+rlapse*ediff))
c           WRM=REDUCE2(real(TNSL),real(rlapse*ediff))
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
        acc=wmean(dist,value)
        call interp1(xxx,yyy,value,real(rlong1),real(rlat1),val)
        acc=val
c ... EXPERIMENTAL arbitrary negative offset 0.1
c        acc=acc-0.1
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






c=========================================
      real function wmean(dist,value)
      dimension dist(4),value(4)
      do i=1,4
        if(dist(i).eq.0) then
          wmean=value(i)
          return
        endif
      enddo
      total=0.
      wmean=0.0
      do i=1,4
        denom=1./dist(i)
        total=total+denom
        wmean=wmean+denom*value(i)
      enddo
      wmean=wmean/total
      end
c=========================================
      subroutine closest(nlong,nlat,rlong,rlat,
     &             rlong1,rlat1,j,i)
      dimension rlong(nlong),rlat(nlat)
      j=0
      dmin=1e30
      do n=1,nlong
        dist=abs(rlong(n)-rlong1)
        if(dist.lt.dmin) then
          dmin=dist
          j=n
        endif
      enddo
      i=0
      dmin=1e30
      do n=1,nlat
        dist=abs(rlat(n)-rlat1)
        if(dist.lt.dmin) then
          dmin=dist
          i=n
        endif
      enddo
      print *,j,i
      end
c=========================================
      subroutine bracket(nlong,nlat,rlong,rlat,
     &             rlong1,rlat1,j,i,dist)
      dimension rlong(nlong),rlat(nlat),dist(4)
      call finder(rlong,nlong,rlong1,j)
      call finder(rlat,nlat,rlat1,i)
c      print *,j,i
c      print *,rlong(j),rlong1,rlong(j+1)
c      print *,rlat(i),rlat1,rlat(i+1)
      if(j.lt.1) j=1
      if(j.ge.nlong) j=nlong-1
      if(i.lt.1) i=1
      if(i.ge.nlat) i=nlat-1
      dist(1)=(rlong(j  )-rlong1)**2+(rlat(i  )-rlat1)**2
      dist(2)=(rlong(j+1)-rlong1)**2+(rlat(i  )-rlat1)**2
      dist(3)=(rlong(j+1)-rlong1)**2+(rlat(i+1)-rlat1)**2
      dist(4)=(rlong(j  )-rlong1)**2+(rlat(i+1)-rlat1)**2
c      print *,dist
c      pause
      end
c=========================================
      SUBROUTINE finder(XX,N,X,J)
      DIMENSION XX(N)
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
      J=JL
      RETURN
      END
C===============================================
      FUNCTION REDUCE(DT)
      const=0.007D0
c      const=0.03D0
      IF(DT.GT.0.0) THEN
        S=const
      ELSEIF(DT.GT.-10.0) THEN
        S=const-const*DT/10.
      ELSE
        S=const*2
      ENDIF
      REDUCE=(1.0+S)**DT
      END

C===============================================
      FUNCTION REDUCE2(tt,dt)
c
      TF=0.67d0*(tt+273.16d0)+88.9d0
      TERM1=-9.09718d0*(273.16d0/TF-1.0d0)                                    
      TERM2=-3.56654d0*LOG10(273.16d0/TF)                                   
      TERM3=0.876793d0*(1.0d0-TF/273.16d0)+0.785835d0
      EXPON=TERM1+TERM2+TERM3                                           
      ES0=10.d0**EXPON                                                     
c
      TF=0.67d0*(TT+DT+273.16d0)+88.9d0
      TERM1=-9.09718d0*(273.16d0/TF-1.0d0)                                    
      TERM2=-3.56654d0*LOG10(273.16d0/TF)                                   
      TERM3=0.876793d0*(1.0d0-TF/273.16d0)+0.785835d0
      EXPON=TERM1+TERM2+TERM3                                           
      ES=10.d0**EXPON    
c
      REDUCE2=ES/ES0
      end                                                 
