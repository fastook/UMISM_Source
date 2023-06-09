      parameter(nmax=9999)
      dimension age(nmax),delox(nmax)
      character*80 junk
      read(1,100) junk
      print *,'input time offset, temp offset'
      read(*,*) timeset,tempset
100   format(a80)
      do i=1,nmax
        read(1,*,end=999) age(i),depth,ox,ox4,delox(i)
        npts=i
      enddo
999   continue
      do i=npts,1,-1
	  time=age(npts)-age(i)-timeset
	  tnsl=delox(i)-tempset
	  write(11,*) time*1000.,tnsl
      enddo
      end
