      subroutine threed(ichk,nmax,numnp,numel,x,y,htice,bdrock,kx,batch)
      implicit real*8(a-h,o-z)
      logical batch
      dimension x(nmax),y(nmax),htice(nmax),bdrock(nmax)
      dimension kx(nmax,4)
      integer ibool(12)
      real*4 obmat(4,4),ltmat(4,4)
c     data ibool /0,0,1,0,1,0,1,0,0,0,1,1/
      data ibool /0,0,1,0,1,0,0,1,0,0,1,1/
c     data obmat /1.,0.,0.,0.,
c    &            0.,1.,0.,0.,
c    &            0.,0.,1.,0.,
c    &            0.,0.,0.,1./
c     data ltmat /1.,0.,0.,0.,
c    &            0.,1.,0.,0.,
c    &            0.,0.,1.,0.,
c    &            0.,0.,0.,1./

      data obmat /0.999840, -0.017844, -0.000444, 0.000000,
     &            0.012326, 0.708199, -0.705904, 0.000000,
     &            0.012911, 0.705787, 0.708306, 0.000000,
     &            0.000000, 0.000000, 0.000000, 1.000000/
 

      data ltmat / 0.298508,0.274015,-0.914227,0.000000,
     &             0.377915, 0.845669, 0.376861, 0.000000,
     &             0.876399, -0.457996, 0.148885, 0.000000,
     &             0.000000, 0.000000, 0.000000, 1.000000/

      integer icolors(3,256)
c     save ibool,obmat,ltmat
      if(ichk.eq.1) then
        ibool(12)=0
      else
        ibool(12)=1
      endif
      iunit=91
      open(iunit,file='color.map')
      do i=1,256
        read(iunit,*,end=999) (icolors(j,i),j=1,3)
        ncolors=i
      enddo
999   continue
      close(iunit)
      open(iunit,file='tmp.color')
      write(iunit,*) ncolors
      do i=1,ncolors
        write(iunit,*) (icolors(j,i)/255.,j=1,3)
      enddo
      close(iunit)
      open(iunit,file='tmp.data',form='formatted')
      write(iunit,*) numnp,numel
      do i=1,numnp
        write(iunit,100) real(x(i))/1000
      enddo
      do i=1,numnp
        write(iunit,100) real(y(i))/1000
      enddo
      do i=1,numnp
        write(iunit,100) real(htice(i)),real(bdrock(i))
      enddo
      do i=1,numel
        do n=1,4
          write(iunit,200) kx(i,n)
        enddo
      enddo
      close(iunit)
      if(.true.) then
        open(iunit,file='backedup.data')
        write(iunit,*) (ibool(i),i=1,12)
c       write(*,*) (ibool(i),i=1,12)
        do j=1,4
          write(iunit,*) (obmat(i,j),i=1,4)
c         write(*,*) (obmat(i,j),i=1,4)
        enddo
        do j=1,4
          write(iunit,*) (ltmat(i,j),i=1,4)
c         write(*,*) (ltmat(i,j),i=1,4)
        enddo
        close(iunit)
100     format(1x,g15.8)
200     format(1x,i10)
      endif
      if(.not.batch) call system('threedc.x')
      if(.true.) then
        open(iunit,file='backup.data')
c       write(*,*) 'reading backup...'
        read(iunit,*) (ibool(i),i=1,12)
c       write(*,*) (ibool(i),i=1,12)
        do j=1,4
          read(iunit,*) (obmat(i,j),i=1,4)
c         write(*,*) (obmat(i,j),i=1,4)
        enddo
        do j=1,4
          read(iunit,*) (ltmat(i,j),i=1,4)
c         write(*,*) (ltmat(i,j),i=1,4)
        enddo
      endif
      end

