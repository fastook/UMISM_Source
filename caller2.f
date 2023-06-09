      PARAMETER(NMAX=79999)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER HED*80
      DIMENSION KODE(NMAX),X(NMAX),Y(NMAX),HTICE(NMAX),ADOT(NMAX),
     &          FRACT(NMAX),PSURF(NMAX),BDROCK(NMAX),FLOWA(NMAX),
     &          SLDGB(NMAX),KX(NMAX,4),CONST(NMAX),AFUDGE(NMAX),
     &          IBFLUX(NMAX,2),BFLUX(NMAX),temp(nmax),itype(nmax),
     &          ACON(NMAX),GEOFLUX(NMAX),calv(nmax)
      dimension lm(5)
      open(1,file='inpute.data')
      READ(1,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      PRINT *,HED
1000  FORMAT (A80,/,7I6,F8.1)
      IF(NUMNP.GT.NMAX) THEN
        PRINT *,'NUMNP=',NUMNP,' NMAX=',NMAX,' INCREASE NMAX'
        STOP
      ENDIF
      print *,numnp,' nodes,',numel,' elements'
      DO N=1,NUMNP
        READ(1,1001) NUM,KODE(N),X(N),Y(N),HTICE(N),
     &                   ADOT(N),FRACT(N),PSURF(N),
     &                   BDROCK(N),FLOWA(N),SLDGB(N),
     &                   temp(n),itype(n),AFUDGE(N),GEOFLUX(N),
     &                   calv(n)
1001  FORMAT(I6,I4,1P2E12.5,0PF10.2,F7.2,F9.2,F10.3,F10.1,F10.5,2F10.5,
     &       I5,F10.3,F10.3,1PE10.3)
      enddo
      DO N=1,NUMEL
        READ(1,1002) NUM,(KX(NUM,i),i=1,4),CONST(N),ACON(N)
 1002 FORMAT(5I6,1P2E17.10)
      enddo
      IF(NUMGBC.GT.0) THEN
        DO N=1,NUMGBC
          READ(1,1007) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
 1007 FORMAT(2I6,E13.6)
        enddo
      ENDIF
      call threed(3,nmax,numnp,numel,x,y,htice,bdrock,kx)
      print *,'called'
      END

      subroutine threed(ichk,nmax,numnp,numel,x,y,htice,bdrock,kx)
      implicit real*8(a-h,o-z)
      dimension x(nmax),y(nmax),htice(nmax),bdrock(nmax)
      dimension kx(nmax,4)
      integer ibool(11)
      real*4 obmat(4,4),ltmat(4,4)
c     data ibool /0,0,1,0,1,0,1,0,0,0,1/
      data ibool /0,0,1,0,1,0,0,1,0,0,1/
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
      iunit=11
      open(iunit,file='color.map.rain')
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
      open(iunit,file='tmp2.data',form='formatted')
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
        write(iunit,*) (ibool(i),i=1,11)
c       write(*,*) (ibool(i),i=1,11)
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
      call system('threed2.x')
      if(.true.) then
        open(iunit,file='backup.data')
c       write(*,*) 'reading backup...'
        read(iunit,*) (ibool(i),i=1,11)
c       write(*,*) (ibool(i),i=1,11)
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


      SUBROUTINE SETRIG
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      PI=4.D0*ATAN(1.D0)
      RADIUS=2.D4/PI
      RADIUS=RADIUS*0.53
      CIRCUM=2.D0*PI*RADIUS
      RKMPDEG=CIRCUM/360.D0
      RADPDEG=PI/180.D0
      DEGPRAD=180.D0/PI
      END
      SUBROUTINE POLREC(RLAT,RLONG,X,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      y= 1000*rkmpdeg*rlat
      x= 1000*rkmpdeg*cos(rlat*radpdeg)*(rlong+127.5)
      END
      SUBROUTINE RECPOL(X,Y,RLAT,RLONG)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      rlat=y*0.001/rkmpdeg
      rlong=-127.5+x*0.001/rkmpdeg/cos(rlat*radpdeg)
      END
