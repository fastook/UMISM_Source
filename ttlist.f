      parameter(nmax=1000)
      dimension ttime(nmax),tlist(nmax)
      print *,'input tnsl offset'
      read(*,*) offset
      do i=1,nmax
        read(1,*,end=999) ttime(i),tlist(i)
        ttime(i)=ttime(i)*1e6
        npts=i
      enddo
      print *,'increase nmax',nmax
      stop
999   continue
      print *,'what is time zero in years?'
      read(*,*) tzero
      do i=1,npts
        write(11,*) tzero-ttime(i),tlist(i)+offset
      enddo
      end
      
