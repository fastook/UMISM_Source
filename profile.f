      parameter(nmax=100000)
      dimension ttime(nmax),dtime(nmax)
      character junk*25
      dtmin=1e30
      dtmax=-1e30
      do i=1,nmax
        read(1,100,end=99) junk,ttime(i),dtime(i)
        print *,junk
100     format(a25,g13.6,g13.6)
        dtmin=min(dtmin,dtime(i))
        dtmax=max(dtmax,dtime(i))
        n=i
      enddo
      print *,'increase nmax=',nmax
99    continue
      ttmin=0.0
      ttmax=ttime(n)
      call grstrt(800,800)
      call linclr(1)
      call window(0.,real(n+1),ttmin,ttmax)
      call move(1.,ttime(1))
      do i=2,n
        call draw(real(i),ttime(i))
      enddo
      call window(0.,real(n+1),dtmin,dtmax)
      call linclr(2)
      call move(1.,dtime(1))
      do i=2,n
        call draw(real(i),dtime(i))
      enddo
      call grstop
      end
