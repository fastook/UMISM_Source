      parameter(nmax=10000)
      dimension x1(nmax),y1(nmax)
      dimension x2(nmax),y2(nmax)
      character junk*25
100   format(a25,f17.0,f13.0)
      do i=1,nmax
        read(1,100,end=10) junk,y1(i),x1(i)
        n1=i
      enddo
10    continue
      do i=1,nmax
        read(2,100,end=20) junk,y2(i),x2(i)
        n2=i
      enddo
20    continue
      ymax=-1e30
      do i=1,n1
        ymax=max(ymax,y1(i))
      enddo
      do i=1,n2
        ymax=max(ymax,y2(i))
      enddo
      call grstrt(600,600)
      call window(0.0,real(max(n1,n2)),0.0,ymax)
      call linclr(2)
      call move(0.0,y1(1))
      do i=2,n1
        call draw(real(i),y1(i))
      enddo
      call linclr(1)
      call move(0.0,y2(1))
      do i=2,n2
        call draw(real(i),y2(i))
      enddo
      call linclr(2)
      call move(0.0,x1(1))
      do i=2,n1
        call draw(real(i),10*x1(i))
      enddo
      call linclr(1)
      call move(0.0,x2(1))
      do i=2,n2
        call draw(real(i),10*x2(i))
      enddo
      call grstop
      end

