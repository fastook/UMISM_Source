      read(1,*) xmin,xmax,ymin,ymax
      call grstrt(800,800)
      call window(xmin,xmax,ymin,ymax)
      read(1,*) rc
      ncurves=int(rc)
      do ic=1,ncurves
        read(1,*,end=999) rn
        n=int(rn)
        read(1,*) x,y
        call move(x,y)
        do i=2,n
          read(1,*) x,y
          call draw(x,y)
        enddo
      enddo
999   continue
      call grstop
      end
