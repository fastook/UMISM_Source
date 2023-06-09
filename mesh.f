      parameter(nmax=10000)
      xs=1e30
      nx=0
      do i=1,nmax
        read(1,*,end=99) x,y,z
        nx=nx+1
        npt=i
        if(x.ne.xs) then
          nx=1
          xs=x
        endif
      enddo
99    continue
      ny=nx
      nx=npt/ny
      !nt=nx;nx=ny;ny=nt
      print *,npt,nx,ny
      do ix=1,nx-1
        ioff=(ix-1)*ny
        do iy=1,ny-1 
          kx1=iy+ioff
          kx2=iy+ny+ioff
          kx3=iy+ny+ioff+1
          kx4=iy+ioff+1
          write(11,*) kx1-1,kx2-1,kx4-1
          write(11,*) kx2-1,kx3-1,kx4-1
c         write(*,*) kx1-1,kx2-1,kx4-1
c         write(*,*) kx2-1,kx3-1,kx4-1
        enddo
      enddo
      end

        
