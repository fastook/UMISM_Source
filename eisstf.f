      parameter(nmax=100,mmax=40)
      dimension x(nmax,mmax),y(nmax,mmax),t(nmax,mmax)
      dimension vx(nmax,mmax),vy(nmax,mmax),vz(nmax,mmax)
      character*14 out(nmax,mmax),junk
      print *,'input first, display interval and how many to show'
      read(*,*) ifirst,istep,inumber
      read(98,*) njunk1,npts,njunk2,ires
      icount=0
      ipr=0
1     read(95,100,end=999) time
      print *,time
100   format(7x,f15.0)
      do i=1,npts
        do j=1,mmax
          read(95,*) x(i,j),y(i,j),vx(i,j),vy(i,j),vz(i,j),
     &               t(i,j)
        enddo
      enddo
      icount=icount+1
      if(icount.gt.ifirst-1+inumber*istep) stop
      if(icount.ge.ifirst .and. 
     *   1+mod(icount-ifirst,istep).eq.1) then
        ipr=ipr+1
        print *,'outputting for display',icount,time,ipr
      write(11,*) ' RANK 2 '
      write(11,*) ' DIMENSIONS ',mmax,npts
      write(11,*) ' TIME ',time
c
      write(11,*) ' NAME temperature '
      write(11,*) ' SCALAR '
      write(11,*) ' DATA '
      do i=1,npts
        do j=1,mmax
          if(t(i,j).eq.-99999.) then
            write(junk,*) 'Missing'
          else
            write(junk,1000) t(i,j)
          endif
          out(i,j)=junk
1000  format(1pg14.6)
        enddo
      enddo
      write(11,*) ((out(i,j),j=1,mmax),i=1,npts)
      write(11,*) ' INTERLACED '
      write(11,*) ' VECTOR 2 '
      write(11,*) ' GRID '
      do i=1,npts
        do j=1,mmax
          write(11,*) x(i,j),y(i,j)
        enddo
      enddo
      write(11,*) ' END '
c
      write(12,*) ' RANK 2 '
      write(12,*) ' DIMENSIONS ',mmax,npts
      write(12,*) ' TIME ',time
      write(12,*) ' NAME x-velocity '
      write(12,*) ' SCALAR '
      write(12,*) ' DATA '
      do i=1,npts
        do j=1,mmax
          if(t(i,j).eq.-99999.) then
            write(junk,*) 'Missing'
          else
            write(junk,1000) vx(i,j)
          endif
          out(i,j)=junk
        enddo
      enddo
      write(12,*) ((out(i,j),j=1,mmax),i=1,npts)
      write(12,*) ' INTERLACED '
      write(12,*) ' VECTOR 2 '
      write(12,*) ' GRID '
      do i=1,npts
        do j=1,mmax
          write(12,*) x(i,j),y(i,j)
        enddo
      enddo
      write(12,*) ' END '
c
      write(13,*) ' RANK 2 '
      write(13,*) ' DIMENSIONS ',mmax,npts
      write(13,*) ' TIME ',time
      write(13,*) ' NAME y-velocity '
      write(13,*) ' SCALAR '
      write(13,*) ' DATA '
      do i=1,npts
        do j=1,mmax
          if(t(i,j).eq.-99999.) then
            write(junk,*) 'Missing'
          else
            write(junk,1000) vy(i,j)
          endif
          out(i,j)=junk
        enddo
      enddo
      write(13,*) ((out(i,j),j=1,mmax),i=1,npts)
      write(13,*) ' INTERLACED '
      write(13,*) ' VECTOR 2 '
      write(13,*) ' GRID '
      do i=1,npts
        do j=1,mmax
          write(13,*) x(i,j),y(i,j)
        enddo
      enddo
      write(13,*) ' END '
c
      write(14,*) ' RANK 2 '
      write(14,*) ' DIMENSIONS ',mmax,npts
      write(14,*) ' TIME ',time
      write(14,*) ' NAME z-velocity '
      write(14,*) ' SCALAR '
      write(14,*) ' DATA '
      do i=1,npts
        do j=1,mmax
          if(t(i,j).eq.-99999.) then
            write(junk,*) ' Missing'
          else
            write(junk,1000) vz(i,j)
          endif
          out(i,j)=junk
        enddo
      enddo
      write(14,*) ((out(i,j),j=1,mmax),i=1,npts)
      write(14,*) ' INTERLACED '
      write(14,*) ' VECTOR 2 '
      write(14,*) ' GRID '
      do i=1,npts
        do j=1,mmax
          write(14,*) x(i,j),y(i,j)
        enddo
      enddo
      write(14,*) ' END '
      endif
 



      goto 1
999   end
