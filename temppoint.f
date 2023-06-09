      parameter(nmax=141,mmax=40)
      dimension x(nmax,mmax),y(nmax,mmax),t(nmax,mmax)
      dimension vx(nmax,mmax),vy(nmax,mmax),vz(nmax,mmax)
      dimension zz(nmax,mmax),zsave(mmax),ysave(mmax)
      dimension dsave(nmax),tsave(nmax)
      dimension dcore(nmax),tcore(nmax)
      dimension dcores(nmax),tcores(nmax)
      dimension icol(nmax,mmax),xy(3,4),ic(4),x2(2)
      character*14 out(nmax,mmax),junk
      data ipass /0/
      read(99,10) junk
10    format(a14)
      dbase=3234.17
      tbase=-31.7731
      bbase=327.165
      do i=1,nmax
        read(99,*,end=20) dcore(i),tcore(i)
        ncore=i
      enddo
20    continue
      fact=(dbase-bbase)/(dcore(1)-dcore(ncore))
      do i=1,ncore
        dcore(i)=fact*dcore(i)
      enddo
      shiftd=dbase-dcore(1)
      shiftt=tbase-tcore(1)
      do i=1,ncore
        dcores(i)=dcore(i)+shiftd
        tcores(i)=tcore(i)+shiftt
c        write(*,*) dcore(i),dcores(i),tcore(i),tcores(i)
      enddo
      print *,'input: 1-temperature '
      print *,'input: 2-x-velocity '
      print *,'input: 3-y-velocity '
      print *,'input: 4-z-velocity '
      read(*,*) iplot
      print *,'input first, display interval and how many to show'
      read(*,*) ifirst,istep,inumber
      print *,'input 0 for no dump, 1 for screen dump'
      read(*,*) idump
      read(98,*) igrip,npts,mpts,ires
      icc=0
      do j=1,mpts
        do i=1,npts
          icc=icc+1
          if(icc.eq.igrip) then
            igrip1=i
            goto 23
          endif
        enddo
      enddo
23    continue
      print *,'GRIP core at column',igrip1
      do igrip=igrip1-0,igrip1+0
      rewind 95
      icount=0
      ipr=0
1     read(95,100,end=999) time
100   format(7x,f15.0)
      xmin=1e30
      xmax=-1e30
      ymin=1e30
      ymax=-1e30
      do i=1,npts
        do j=1,mmax
          read(95,*) x(i,j),y(i,j),vx(i,j),vy(i,j),vz(i,j),
     &               t(i,j)
          if(iplot.eq.1) then
            zz(i,j)=t(i,j)
            zmax=0.0
            zmin=-60.
          elseif(iplot.eq.2) then
            zz(i,j)=vx(i,j)
            zmax=10.
            zmin=-10.
          elseif(iplot.eq.3) then
            zz(i,j)=vy(i,j)
            zmax=10.
            zmin=-10.
          elseif(iplot.eq.4) then
            zz(i,j)=vz(i,j)
            zmax=1.
            zmin=-0.
          endif
          xmin=min(xmin,x(i,j))
          xmax=max(xmax,x(i,j))
          ymin=min(ymin,y(i,j))
          ymax=max(ymax,y(i,j))
        enddo
      enddo
      if(ipass.eq.0) then
        ipass=1
        call grstrt(400,400)
        tmin=-60.
        tmax=0.
        hmin=0.
        hmax=4000.
        call window(tmin,tmax+5.,hmin,hmax)
        tlab=tmin+(tmax-tmin)/20
        hlab=hmax-(hmax-hmin)/20.
      endif
      icount=icount+1
      if(icount.gt.ifirst-1+inumber*istep) goto 999
      if(icount.ge.ifirst .and. 
     *   1+mod(icount-ifirst,istep).eq.1) then
        ipr=ipr+1
        print *,'number displayed',ipr,time
        if(ipr.gt.1) then
          if(mod(ipr,10).eq.0) then
            call linclr(14)
          else
            call linclr(0)
          endif
          call linclr(0)
          call move(zsave(1),ysave(1))
c          call marker(zsave(1),ysave(1),10)
          do i=2,mmax
            call draw(zsave(i),ysave(i))
c            call marker(zsave(i),ysave(i),10)
          enddo
          call linclr(3)
          do i=1,ncore
            call marker(tcores(i),dcores(i),10)
          enddo
          call move(tlab,hlab)
          call linclr(0)
          call text(14,junk)
          write(junk,'(1pg13.6)') time
          call linclr(1)
          call text(14,junk)
        endif
        call linclr(14)
        do tline=tmin,0.0,5.
          call move(tline,hmin)
          call draw(tline,hmax)
        enddo
        call linclr(2)
c        print *,zz(igrip,1),y(igrip,1),y(igrip,mmax)
        call move(zz(igrip,1),y(igrip,1))
c        call marker(zz(igrip,1),y(igrip,1),10)
        zsave(1)=zz(igrip,1)
        ysave(1)=y(igrip,1)
        do i=2,mmax
          call draw(zz(igrip,i),y(igrip,i))
c          call marker(zz(igrip,i),y(igrip,i),10)
          zsave(i)=zz(igrip,i)
          ysave(i)=y(igrip,i)
        enddo
        if(idump.eq.1) call savescreen(ipr,400,400)
      else
        print *,'number read     ',icount,time
      endif
c
      goto 1
999   continue
      enddo

      call grstop
      end


      subroutine wait(n)
      do i=1,n
        x=sin(real(i))
      enddo
      end
