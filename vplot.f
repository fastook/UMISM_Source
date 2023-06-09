      parameter(nmax=29999,nvmax=1000)
      character*80 junk
      dimension x1(nmax),y1(nmax),vx1(nmax),vy1(nmax),vmag1(nmax)
      dimension x2(nmax),y2(nmax),vx2(nmax),vy2(nmax),vmag2(nmax)
      real*8 rand
      call srand(23)
      read(1,1000) junk
      read(2,1000) junk
1000  format(a80)
      do i=1,nmax
        read(1,*,end=98) n,x1(i),y1(i),vx1(i),vy1(i),vmag1(i)
        npt1=i
      enddo
      print *,'increase nmax',nmax
      stop
98    continue
      do i=1,nmax
        read(2,*,end=99) n,x2(i),y2(i),vx2(i),vy2(i),vmag2(i)
        npt2=i
      enddo
      print *,'increase nmax',nmax
      stop
99    continue
      istep=nmax
      percent=2.
      if(npt1.gt.nvmax) then
        print *,npt1
        percent=1.-real(nvmax)/real(npt1)
        istep=nvmax/(npt1-nvmax)
        print *,'too many, eliminating',100*percent,'percent',istep
      endif
      xmin=1e30
      ymin=xmin
      xmax=-xmin
      ymax=xmax
      do i=1,npt1
        xmin=min(xmin,x1(i))
        xmax=max(xmax,x1(i))
        ymin=min(ymin,y1(i))
        ymax=max(ymax,y1(i))
      enddo
      do i=1,npt2
        xmin=min(xmin,x2(i))
        xmax=max(xmax,x2(i))
        ymin=min(ymin,y2(i))
        ymax=max(ymax,y2(i))
      enddo
      print *,xmin,xmax,ymin,ymax
      call grstrt(800,800)
      call window(xmin,xmax,ymin,ymax)
1     continue
      print *,'input scale1,scale2',scale1,scale2
      read(*,*,end=999) scale1,scale2
      call newpag
      call linclr(1)
      nc=0
      call srand(23)
      do i=1,npt1
        if(rand().gt.percent) nc=nc+1
      enddo
      print *,npt1,nc
      write(11,1001) nc,xmin,xmax,ymin,ymax
1001  format(1x,'$ptype'/
     &       3x,'irot=1,'/
     &       3x,'ncurv=',i5,','/
     &       3x,'ititle=1,1,1,'/
     &       3x,'nsymb=0,itype=0,modec=2,iclip=1,'/
     &       3x,'idebug=1,1,1,1,1,1,1,1,1,1,1,'/
     &       1x,'$end'/
     &       1x,'$axes'/
     &       3x,'iscale=1,'/
     &       3x,'xmin=',g14.7,','/
     &       3x,'xmax=',g14.7,','/
     &       3x,'ymin=',g14.7,','/
     &       3x,'ymax=',g14.7,','/
     &       3x,'ixaxis=2,iyaxis=2,'/
     &       1x,'$end'/
     &       1x,'X (KM)'/
     &       1x,'Y (KM)'/
     &       1x,'TITLE')
      call srand(23)
      do i=1,npt1
        if(rand().gt.percent) then
          call move(x1(i),y1(i))
          xp=x1(i)+scale1*vx1(i)
          yp=y1(i)+scale1*vy1(i)
          call draw(xp,yp)
          write(11,*) 2
          write(11,*) x1(i),y1(i)
          write(11,*) xp,yp
        endif
      enddo
      call linclr(2)
      write(12,1001) nc,xmin,xmax,ymin,ymax
      call srand(23)
      do i=1,npt2
        if(rand().gt.percent) then
          call move(x2(i),y2(i))
          xp=x2(i)+scale2*vx2(i)
          yp=y2(i)+scale2*vy2(i)
          call draw(xp,yp)
          write(12,*) 2
          write(12,*) x2(i),y2(i)
          write(12,*) xp,yp
        endif
      enddo
      goto 1
999   continue
      call grstop1
      end



