      parameter(ndim=281,nmax=ndim*ndim,nlevels=40)
      dimension temp(nlevels,ndim,ndim)
      character*80 hed
      read(10,100) hed
100   format(a80)
      read(10,*) nnode, nelem, ncol,nrow
      do i=1,nrow
        do j=1,ncol
          do k=1,nlevels
            read(11,*) ii,jj,temp(k,j,i)
          enddo
        enddo
      enddo
      write(20,*) 'RANK 3'
c     write(20,*) 'DIMENSIONS',nlevels,ncol,nrow
      write(20,*) 'DIMENSIONS',ncol,nrow,nlevels
      write(20,*) 'BOUNDS -1. 1. -1. 1. -1. 1. '
      write(20,*) 'NAME TEMPERATURE'
      write(20,*) 'TIME 0'
      write(20,*) 'SCALAR'
      write(20,*) 'DATA'
c     write(20,*) (((temp(k,j,i),k=1,nlevels),j=1,ncol),i=1,nrow)
      write(20,*) (((temp(k,j,i),j=1,ncol),i=1,nrow),k=nlevels,1,-1)
      write(20,*) 'END'
      end
      
