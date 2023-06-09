      parameter(nmax=10000)
      real time(nmax),theta(nmax),vmag(nmax)
      real cumm1(0:359),cumm2(0:359)
      tstart=-120000.
      tstop=0
      print *,'input start/stop time'
      read(*,*) tstart,tstop
      open(1,file='rose.data')
      open(11,file='rose1.d')
      open(12,file='rose2.d')
      open(13,file='rose.time')
      write(13,100) tstart,tstop
100   format(1x,'0.2 5.5 10 0 1 1 time:',f8.0,'->',f8.0)
      read(1,*,end=99) timei,thetai,vmagi
      do i=0,359
        cumm1(i)=0.0
        cumm2(i)=0.0
      enddo
      do i=1,nmax
        read(1,*,end=99) time(i),theta(i),vmag(i)
        n=i
      enddo
      print *,'increase nmax:',nmax
      stop
99    continue
      do i=2,n
        ttime=(time(i)+time(i-1))/2.
        if(tstart.lt.ttime .and. ttime.lt.tstop) then
          dt=time(i)-time(i-1)
          if(theta(i).ne.-999.) then
            iang=int(theta(i))
              cumm1(iang)=cumm1(iang)+dt
              cumm2(iang)=cumm2(iang)+dt*vmag(i)
          endif
        endif
      enddo
        do i=0,359
        if(cumm1(i).gt.0.) then
          write(11,'(1x,i4,f16.1)') i,cumm1(i)
        endif
        if(cumm2(i).gt.0.) then
          write(12,'(1x,i4,1pe16.6)') i,cumm2(i)/1000.
        endif
      enddo
      end

        
      
