#include "parameter.h"
      parameter( nmax=maxnum,nz=9,nz1=nz+1,nsmax=maxtime,n3=nmax*3,
     &           n23=n3*2)
      implicit real*8(a-h,o-z)
c ***************************************************************c
c                                                                c
c   program:  map5                                               c
c                                                                c
c   date:  7-31-95                                               c
c   programmer:  j.  fastook                                     c
c                                                                c
c   function:                                                    c
c           this is a program to model the flow of a glacier     c
c           with material properties read for each nodal point   c
c           and the average used for the element.                c
c           it calculates an element matrix ala becker, 2-d      c
c           program with linear shape functions.                 c
c           does time dependent case using a lumped capacitance  c
c           matrix and a backward difference scheme. allowing    c
c           for very quick operation.                            c
c ***************************************************************c
c     program for steady and unsteady state flow analysis
c     using finite elements
c
c ***************************************************************c
      dimension www(n3),thickl(nmax),wrate(n3,2)
      dimension wwworig(nmax),wdiff(nmax)
      character hed*80,scrtch*80,iadj*2
      real*8 a(nmax,nz)
      common /mat/  a,ka(nmax,nz1),numnp
      common /lapse/ acom,hmax,windir(2),xpole,ypole,abl,tnslbase
      common /velos/ ascal,umax,vthresh,inorm
      dimension ifit(6), amass(11), b(nmax), x(nmax), y(nmax), d(nmax),
     &          bold(nmax),flux(nmax),adiag(nmax),itkode(nmax),
     &          kode(nmax), const(nmax), acon(nmax), lm(5), itype(nmax),
     &          kx(nmax,4), idt(nmax), adotb(nmax), 
     &          adot(nmax), bdrock(nmax), flowa(nmax), sldgb(nmax),
     &          psurf(nmax), ppsurf(nmax), fract(nmax),
     &          cnew(nmax), qhold(nmax), htice(nmax), thick(nmax),
     &          hfit(nmax), ibflux(nmax,2), bflux(nmax), depb(nmax),
     &          q(nmax),temp(nmax),tbed(nmax),vel(nmax,3),
     &          undepb(nmax), slope(4,nmax), slopn(4,nmax), kz(nmax),
     &          wthick(nmax),wtelem(nmax),geoflux(nmax),
     &          calv(nmax),pcalv(nmax),tsorig(nmax),htice0(nmax)
      dimension afudge(nmax),adp(nmax),adm(nmax),adc(nmax),bmelt(nmax)
      dimension ttime(nsmax),vvol(nsmax),aarea(nsmax),ttbot(nsmax),
     &          ttnsl(nsmax),twater(nsmax),pwater(nsmax),wwwmin(nsmax)
      dimension ntype(nmax), nnode(nmax),
     &          xi(2,9), eta(2,9), w(2,9),
     &          aadot(nmax), afract(nmax), aflowa(nmax),
     &          abdrck(nmax), asldgb(nmax)
      dimension alphac(3)
c     dimension dpsix(9), dpsiy(9), dxds(2,2), dsdx(2,2),
c    &          psi(4), dpsi(4,2), cnst(nmax),
c    &          xy(2,4)
      logical ctogg,wtogg,itogg,btogg,iotogg
      common /toggles/ ctogg,wtogg,itogg,btogg,iotogg,ipage
      parameter(npage=39)
      character*80 list(1000)
      common /iolist/ list
      common /line/ np,nline(1000)
      logical iflush
      common /flush/ iflush
      data iflush /.true./
c ... timer stuff, for sgi only ...
      real*4 tb(2),etime,dtime
      external etime,dtime
c ... timer stuff, for sgi only ...
      real*8 t(nmax),tnew(nmax)
      external warming
      data hfit /nmax*0.0/
      data ascal /1.d0/, umax /100.d0/, vthresh /1000.d0/
      data inorm /1/
      data www /n3*0.d0/
      data wrate /n23*0.d0/
123   format(a25,t30,1pg13.6,g13.6)
1010  format(13i6)
      rewind 20
      read(20,1010,end=102) np
        read(20,1010) (nline(ip),ip=1,np)
      write(*,1010) np
        write(*,1010) (nline(ip),ip=1,np)
        goto 103
102   continue
        np=0
103   continue
      rewind 20
      iplot=5
      ipass=0
      ntstep=0
      time=0.d0
      icon=0
      call setrig
      big=1d20
      itogg=.false.
      ctogg=.true.
      wtogg=.false.
      btogg=.true.
      iotogg=.true.
      xpole=0.0d0
      ypole=0.0d0
      acom=-9.6237669d0
      windir(1)=90.d0/180.d0*3.14159d0
      windir(2)=100.d0
c     windir(1)=315.d0/180.d0*3.14159d0
c     windir(2)=50.d0
c following is default snowline elevation at pole, gradient, tnsl
      amass(1)=-5.d0
      amass(2)=1.d0
      amass(3)=.1d0
      amass(4)=500.d0
      amass(5)=1000.d0
      amass(6)=1500.d0
      amass(7)=-3000.d0
      amass(8)=.56d-3
      amass(8)=.001d0
      amass(9)=-14.d0
      amass(10)=.4d0
      amass(11)=1.2d0
      cfactor=0.d0
c a(1): advection coefficient, a(2): diffusion, a(3): loss term
      alphac(1)=1d0
      alphac(2)=1d0
      alphac(3)=1d-3
      tbase=0.0d0
      tnslbase=14.d0
c*******************
c ... check fo defaults file. if none exist, create one ...
      read(9,*,end=100) amass,acom,wjunk,windir(2),
     &                  xpole,ypole,cfactor,alphac,tbase,
     &                  tnslbase,ctogg,wtogg,itogg,btogg
c                                t     f     f     t
        windir(1)=wjunk/180.d0*3.14159d0
        print *,'**** defaults file found ****'
        goto 101
100   print *,'no defaults file found'
        rewind 9
        wjunk=windir(1)*180.d0/3.14159d0      
        write(9,*) amass,acom,wjunk,windir(2),
     &             xpole,ypole,cfactor,alphac,tbase,
     &             tnslbase,ctogg,wtogg,itogg,btogg
101   continue
      print *,' tnsl       = ',amass(9)
      print *,' fudge      = ',amass(10)
      print *,' geo grad   = ',amass(11)
      print *,' acom       = ',acom
      print *,' windir     = ',wjunk,windir(2)
      print *,' xpole,ypole= ',xpole,ypole
      print *,' time       = ',tbase
      print *,' tnslbase   = ',tnslbase
c*******************
      do i=1,5
        ifit(i)=0
      enddo
      ifit(6)=1
      ntstep=0
c
c ... following is sealevel referenced to present=0.
      sealev=0.d0
c
c ... following sets rate of convergence, up to 5 works well
      conv=1.0d0
      time=tbase
      rhow=1.092d0
      pg = 0.089866d0*0.3816d0
      ncol=nz
c
c ... initialize integration points and weights
c ... gaussian quadrature of order three quadrilaterals
      call gausinit(xi,eta,w)
c
      do i=1,4
        lm(i)=0
      enddo
c
c ... following readn for split data sets
      call readn(nmax,hed,numnp,numel,numgbc,ndt,inter,dt,
     &           kode,x,y,htice,adot,adotb,fract,psurf,rhoi,rhow,
     &           depb,bdrock,
     &           undepb,flowa,acon,sldgb,temp,itype,afudge,geoflux,
     &           thick,kx,const,ibflux,bflux,qhold,ntype,nnode,ncol,
     &           aadot,afract,abdrck,ppsurf,aflowa,asldgb,idt,amass,
     &           numcol,numlev,calv,pcalv,www,wrate,thickl,
     &           wwworig,tsorig)
      acomsave=acom
      ii=1
      do i=1,numnp
        htice0(i)=htice(i)
        wdiff(i)=www(ii)-wwworig(i)
        ii=ii+3
      enddo
      do i=1,numnp
        thickl(i)=htice(i)-bdrock(i)
        flot=(1.d0-rhow/rhoi)*bdrock(i)
        if(htice(i).lt.flot) then
          thickl(i)=0.d0
        endif
      enddo
1     continue
      write(*,*) 'time step=',dt
      write(*,*) 'input 1 to enter adjust and continue, 0 to bypass'
      if(iflush) call gflush
      read(*,4000) iadj
4000  format(a2)
      write(99,4000) iadj
      if(iadj.eq.'1') then
c
c ..... enter interactive data set manipulator
        call adjust(hed, numnp, numel, x, y, htice, adot, 
     &              adotb,fract,bmelt,wthick, 
     &              temp, itype, tbed, psurf, undepb,
     &              bdrock, depb, flowa, sldgb, thick, kx, const, 
     &              afudge,nnode, kode, geoflux, flux,
     &              hfit, numcol, numlev, numgbc, ndt, inter, dt,
     &              ibflux, bflux, nmax, idt, slopn, amass, time,
     &              ntstep, ttime, vvol, aarea,ttbot, ttnsl, ifit,
     &              iplot, cfactor, acon,icon, 
     &              alphac,tbase,
     &              imelt,twater,pwater,calv,
     &              ntype,aadot,afract,abdrck,ppsurf,
     &              aflowa,asldgb,pcalv,adc,wwwmin,www,wrate,wdiff,
     &              wwworig,tsorig)
        print *,' update changes ??? 0-yes, 1 no '
        read(*,*) iupdate
        write(99,*) iupdate
        if(iupdate.ne.1) then
          do i=1,numnp
            htice0(i)=htice(i)
            thickl(i)=htice(i)-bdrock(i)
            flot=(1.d0-rhow/rhoi)*bdrock(i)
            if(htice(i).lt.flot) then
              thickl(i)=0.d0
            endif
            thickl(i)=max(0.d0,thickl(i))
          enddo
        endif 
c
        write(*,*) 'input 1 to continue, 0 to quit'
        if(iflush) call gflush
        read(*,4000) iadj
        write(99,4000) iadj
        if(iadj.eq.'0') stop
      endif
c......experimental section for bed depression ..................
      if(.true.) then
        call grstrt(800,800)
        do l=1,ndt
          time=time+dt
          ntstep=ntstep+1
          if(.false.) then
            call eplate(ntstep,numnp,numel,x,y,kx,thickl,kode,
     &           dt,www,wrate,wwwmin(ntstep),time,wwworig,fnete,www)
          elseif(.false.) then
            call vplate(ntstep,numnp,numel,x,y,kx,thickl,kode,
     &           dt,www,wrate,wwwmin(ntstep),time,wwworig,fnetv,www)
          else
            call veplate(ntstep,numnp,numel,x,y,kx,thickl,kode,
     &           dt,www,wrate,time,wwworig,
     &           wmin,wmine,wminv,fnete,fnetv)
          endif
          print *,l,real(time),': fnete=',real(fnete/1d6),
     &              real(fnetv/1d6)
          vvol(ntstep)=fnete
          aarea(ntstep)=fnetv
          wwwmin(ntstep)=wmin
          ttnsl(ntstep)=0.d0
          twater(ntstep)=wmine
          pwater(ntstep)=wminv
          ttime(ntstep)=time
          if(.true.) then
            do i=1,numnp
              htice(i)=htice0(i)+www((i-1)*3+1)
              depb(i)=undepb(i)+www((i-1)*3+1)
              if(htice(i).lt.depb(i)) then
                htice(i)=depb(i)
c                thickl(i)=0.d0
              else
c                thickl(i)=htice(i)-depb(i)
              endif
c             thickl(i)=htice(i)-depb(i)
             if(depb(i).lt.0) then
               flot=(1.d0-rhow/rhoi)*depb(i)
c               if(htice(i).lt.flot) then
c                 thickl(i)=0.d0
c               endif
             endif
            enddo
          endif
          if(.true.) then
            write(scrtch,*) 'time=',nint(time)
            write(34) scrtch
            write(34) (htice(i),i=1,numnp)
            write(34) (adot(i),i=1,numnp)
            write(34) (depb(i),i=1,numnp)
            write(34) (const(i),i=1,numel)
            write(34) (acon(i),i=1,numel)
            write(35) (fract(i),i=1,numnp)
            write(35) (flowa(i),i=1,numnp)
            write(35) (sldgb(i),i=1,numnp)
            write(35) (afudge(i),i=1,numnp)
            write(36) scrtch
            write(36) (tbed(i),i=1,numnp)
            write(36) (bmelt(i),i=1,numnp)
            write(36) (wthick(i),i=1,numnp)
          endif
          call poutstf(numnp*3,www,wrate(1,1),numnp,thickl,numcol,
     &               numlev,time,wdiff)
c          if(iplot.ne.0) then
          if(.true. .and. ntstep.gt.2) then
            call plotsol(numnp,x,y,htice,depb,kode,fract,
     &                   psurf,wthick,twater,pwater,afudge,
     &                   ntstep,ttime,vvol,
     &                   aarea,ttbot,ttnsl,
     &                   iplot,
     &                   numel,kx,const,vel,wwwmin,dtmin)
            if(iplot.eq.2.or.iplot.eq.3.or.iplot.eq.4) then
              call drawout(idone)
            endif
          endif
        enddo
        call grstop1
      endif
      goto 1
c................................................................
c
      stop
      end
