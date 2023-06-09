c mkplotc takes the output in contour form (made perhaps by contour.f,
c which reassembles the segmented contours that come out of cont3.f) and
c builds a .namelist file that can be used with plotc to draw postscript
c files. 
c
c it requires as input a .data file with the contours assembled.
c                      an outline.data with the continenet outlines
      parameter(nmax=1000,nplots=10000)
      dimension x(nmax,nplots),y(nmax,nplots),npoints(nplots)
      dimension ltype(nplots),ictype(nplots),istype(nplots)
      dimension istypv(nplots)
      character ptitle*80,xtitle*20,ytitle*20
      character*20 legend(nplots),junk
      character*1 comma,equal,quote
      data comma /','/,equal /'='/, quote /''''/
      ioutline=0
      do nplot=1,nplots
        do np=1,nmax
          read(2,*,end=1010) xp,yp,ip
c          write(*,202) xp,yp,ip
          ioutline=1
          if(xp.eq.-99999.) then
            ltype(nplot)=ip+1
            read(2,100) junk
            goto 1000
          endif
          x(np,nplot)=xp
          y(np,nplot)=yp
c          x(np,nplot)=xp-6.
c          y(np,nplot)=yp-7000.
          npoints(nplot)=np
        enddo
        write(*,*) 'too many points',nmax
        stop
1000    continue
      enddo
      write(*,*) 'too many plots',nplots
1010  continue
      read(1,202) x1,y1
      read(1,202) x2,y2
      read(1,202) x3,y3
      read(1,202) x4,y4
      read(1,202) x5,x6
      read(1,100) junk
      xmin=min(x1,x2,x3,x4)
      xmax=max(x1,x2,x3,x4)
      ymin=min(y1,y2,y3,y4)
      ymax=max(y1,y2,y3,y4)
      print *,'case 2 below'
      print *,'xmin,xmax',xmin,xmax
      print *,'ymin,ymax',ymin,ymax
      dx=xmax-xmin
      dy=ymax-ymin
      xmid=(xmax+xmin)/2.
      ymid=(ymax+ymin)/2.
      if(dx.gt.dy) then
        xmin1=xmid-dx/2.
        xmax1=xmid+dx/2.
        ymin1=ymid-dx/2.
        ymax1=ymid+dx/2.
      else
        xmin1=xmid-dy/2.
        xmax1=xmid+dy/2.
        ymin1=ymid-dy/2.
        ymax1=ymid+dy/2.
      endif
      print *,'case 3 below'
      print *,'xmin,xmax',xmin1,xmax1
      print *,'ymin,ymax',ymin1,ymax1
      nsolid=nplot-1
      do nplot=nsolid+1,nplots
        do np=1,nmax
          read(1,202,end=2010) xp,yp,ip
          if(xp.eq.-99999.) then
            ltype(nplot)=ip+1
            read(1,100) legend(nplot-nsolid)
            goto 2000
          endif
          x(np,nplot)=xp
          y(np,nplot)=yp
          npoints(nplot)=np
        enddo
        write(*,*) 'too many points',nmax
        stop
2000    continue
      enddo
      write(*,*) 'too many plots',nplots
2010  continue
      nplot=nplot-1
      write(*,*) '-1 -> autoscale, equal axis'
      write(*,*) '+1 -> autoscale with plotc'
      write(*,*) ' 0 -> manual scaling'
      write(*,*) ' 2 -> use values above'
      write(*,*) ' 3 -> use values above, equal axis'
      read(*,*) iscale
c
c      iscale=0
c
      write(*,*) 'input plot title'
      read(*,200) ptitle
      write(*,*) 'input x-axis title'
c      read(*,100) xtitle
c
      xtitle='X (KM)'
c
      write(*,*) 'input y-axis title'
c      read(*,100) ytitle
c
      ytitle='Y (KM)'
c
      if(iscale.eq.0) then
        write(*,*) 'input xmax'
        read(*,*) xmax
        write(*,*) 'input xmin'
        read(*,*) xmin
        write(*,*) 'input ymax'
        read(*,*) ymax
        write(*,*) 'input ymin'
        read(*,*) ymin
c
c        xmin=500.
c        xmax=625.
c        ymin=625.
c        ymax=750.
      endif
      if(iscale.eq.3) then
        xmax=xmax1
        ymax=ymax1
        xmin=xmin1
        ymin=ymin1
      endif
      do i=1,nplot
        ictype(i)=mod(ltype(i)-1,5)+1
        istype(i)=mod(ltype(i)-1,10)+1
        istypv(i)=mod(ltype(i)-1,10)+1
      enddo
      write(11,*) '$ptype'
      write(11,*) '  irot=',1,comma
      write(11,*) '  ncurv=',-nplot,comma
      write(11,*) '  nsymb=',0,comma
      write(11,*) '  itype=2',comma
c      write(11,*) '  ictype=',(ictype(i),comma,i=1,nplot)
      write(11,1001) '  ictype=',(1,comma,i=1,nplot)
      write(11,1001) '  istype=',(istype(i),comma,i=1,nplot)
      write(11,1001) '  istypv=',(istypv(i),comma,i=1,nplot)
1001  format(a,12(i3,a))
      write(11,*) '  isize=5',comma
      write(11,*) '  ifreq=50',comma
      write(11,*) '  ititle= 1,1,1',comma
      write(11,*) '  ilegnd=-1',comma
      write(11,*) '  xleg=8.1',comma
      write(11,*) '  yleg=6.75',comma
      if(ioutline.eq.1) then
        write(11,*) '  aleg(1,2)=',quote,'outline',quote,comma
        write(11,*) '            ',quote,legend(1),quote,comma
      else
        write(11,*) '  aleg(1,2)=',quote,legend(1),quote,comma
      endif
      do i=2,nplot
        if(legend(i).ne.legend(i-1)) then
          write(11,*) '            ',quote,legend(i),quote,comma
        endif
      enddo
      write(11,*) '  ilegbx=1',comma
      write(11,*) '  modec=2',comma
      write(11,*) '  iclip=1',comma
      write(11,*) '  idebug=1,1,1,1,1,1,1,1,1,1,1',comma
      write(11,*) '$end'
      write(11,*) '$axes'
      if(iscale.eq.0 .or. iscale.eq.2 .or. iscale.eq.3) then
        write(11,*) '  iscale=1',comma
        write(11,*) '  xmin=',xmin,comma
        write(11,*) '  xmax=',xmax,comma
        write(11,*) '  ymin=',ymin,comma
        write(11,*) '  ymax=',ymax,comma
      elseif(iscale.eq.1) then
        write(11,*) '  iscale=0',comma
      elseif(iscale.eq.-1) then
        write(11,*) '  iscale=4',comma
      endif
      write(11,*) '  ixaxis=2',comma
      write(11,*) '  iyaxis=2',comma
      write(11,*) '$end'
      write(11,*) xtitle
      write(11,*) ytitle
      write(11,*) ptitle
      do np=1,nplot
        write(11,*) npoints(np)
c        write(*,*) npoints(np)
        do i=1,npoints(np)
          write(11,*) x(i,np),y(i,np)
c          write(*,*) x(i,np),y(i,np)
       enddo
      enddo
100   format(a20)
200   format(a80)
202   format(10x,g13.6,2x,g13.6,i13)
      end

      

