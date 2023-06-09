      parameter(nmax=141,mmax=40)
c             ------------ #include "fshort.h"
      dimension x(nmax,mmax),y(nmax,mmax),t(nmax,mmax)
      dimension vx(nmax,mmax),vy(nmax,mmax),vz(nmax,mmax)
      dimension zz(nmax,mmax)
      dimension icol(nmax,mmax),xy(3,4),ic(4),x2(2)
      character*14 out(nmax,mmax),junk
      data ipass /0/
      print *,'input: 1-temperature '
      print *,'input: 2-x-velocity '
      print *,'input: 3-y-velocity '
      print *,'input: 4-z-velocity '
      print *,'input: 5-total-velocity '
      print *,'input: 6-ice hardness '
      read(*,*) iplot
      print *,'input -1 to print, 0 not to...'
      read(*,*) iprint
      print *,'input first, display interval and how many to show'
      read(*,*) ifirst,istep,inumber
      read(98,*) njunk1,npts,njunk2,ires
      icount=0
      ipr=0
      xmin=1e30
      xmax=-1e30
      ymin=-1000.
      ymax=5000.
      zmin=1e30
      zmax=-1e30
1     read(95,100,end=999) time
100   format(7x,f15.0)
      do i=1,npts
        do j=1,mmax
          read(95,*) x(i,j),y(i,j),vx(i,j),vy(i,j),vz(i,j),
     &               t(i,j)
          if(iplot.eq.1) then
            zz(i,j)=t(i,j)
            zmax=0.0
            zmin=-50.
            ispace=2
          elseif(iplot.eq.2) then
            zz(i,j)=vx(i,j)
            zmax=50.
            zmin=-50.
            ispace=5
          elseif(iplot.eq.3) then
            zz(i,j)=vy(i,j)
            zmax=50.
            zmin=-50.
            ispace=5
          elseif(iplot.eq.4) then
            zz(i,j)=vz(i,j)
            zmax=1.
            zmin=-0.
            ispace=1
          elseif(iplot.eq.5) then
            if(vx(i,j).eq.-99999. .or. 
     &         vy(i,j).eq.-99999. .or.
     &         vz(i,j).eq.-99999.) then
              zz(i,j)=-99999.
            else
              zz(i,j)=sqrt(vx(i,j)**2+vy(i,j)**2+vz(i,j)**2)
            endif
             zmax=50.
             zmin=0.
            ispace=5
          elseif(iplot.eq.6) then
            zz(i,j)=hardness(t(i,j))
            zmax=18.75
            zmin=0.
            ispace=1
          endif
          if(ipass.eq.0) then
            xmin=min(xmin,x(i,j))
            xmax=max(xmax,x(i,j))
            ymin=min(ymin,y(i,j))
            ymax=max(ymax,y(i,j))
          endif
        enddo
      enddo
      if(ipass.eq.0) then
        ipass=1
        dx=(xmax-xmin)/20.
        call initialize(ncolor,xmin-dx,xmax+5*dx,ymin,2.*ymax)
        zdelt=(zmax-zmin)/ncolor
        dx=(xmax-xmin)/ncolor
        xtime=xmin+dx
        ytime=1.7*ymax
      endif
      icount=icount+1
      print *,'number read',icount,time
      if(icount.gt.ifirst-1+inumber*istep) goto 999
      if(icount.ge.ifirst .and. 
     *   1+mod(icount-ifirst,istep).eq.1) then
        ipr=ipr+1
        if(iprint.ne.-1) call makeob(ipr)
        print *,'number displayed',ipr,time
        write(junk,'(1pg13.6)') time
        call RGBcol(255,255,255)
        call cmov2(xtime,ytime)
        call charst(junk,14)

        imin=2000
        imax=0
        do i=1,npts
          do j=1,mmax
            if(zz(i,j).eq.-99999.) then
              icol(i,j)=-99999
            else
              icol(i,j)=izset(ncolor,zz(i,j),zmin,zdelt,0)
              imin=min(imin,icol(i,j))
              imax=max(imax,icol(i,j))
            endif
          enddo
        enddo
        ncol=1+imax-imin
        dx=(xmax-xmin)/ncolor
        ii=1
        dy=0.1*ymax
        ybox=2.*ymax-dy
        do i=1,ncolor-1
          xbox=xmin+(ii-1)*dx
          ii=ii+1
          call bgnpol
            x2(1)=xbox
            x2(2)=ybox
            call setcolor(i)
            call v2f(x2)
            x2(1)=xbox+dx
            x2(2)=ybox
            call setcolor(i+1)
            call v2f(x2)
            x2(1)=xbox+dx
            x2(2)=ybox+dy
            call setcolor(i+1)
            call v2f(x2)
            x2(1)=xbox
            x2(2)=ybox+dy
            call setcolor(i)
            call v2f(x2)
          call endpol
        enddo
        dy=(2*ymax-0.9*ymin)/(ncolor+10)
        dx=0.1*xmax
        xbox=xmax+dx
c        print *,imin,imax
        imin=10*nint(imin/10.)
        imax=10*nint(imax/10.)
        if(imin.eq.imax) imax=imin+10
        do i=imin+1,imax+1,10
c        do i=1,ncolor+1,10
          ybox=ymin+(i-0)*dy+5*dy
          rnum=zmin+(i-1)*zdelt
          icc=izset(ncolor,rnum,zmin,zdelt,0)    
          write(junk,'(f8.3)') rnum
c          write(*,'(i4,f8.3,i4)') i,rnum,izset(ncolor,rnum,zmin,zdelt,0)
          call setcolor(icc)
          call cmov2(xbox,ybox)
          call charst(junk,14)
        enddo
        do i=1,npts-1
          do j=1,mmax-1
              xy(1,1)=x(i,j)
              xy(2,1)=y(i,j)
              xy(3,1)=zz(i,j)
              ic(1)=icol(i,j)
              xy(1,2)=x(i+1,j)
              xy(2,2)=y(i+1,j)
              xy(3,2)=zz(i+1,j)
              ic(2)=icol(i+1,j)
              xy(1,3)=x(i+1,j+1)
              xy(2,3)=y(i+1,j+1)
              xy(3,3)=zz(i+1,j+1)
              ic(3)=icol(i+1,j+1)
              xy(1,4)=x(i,j+1)
              xy(2,4)=y(i,j+1)
              xy(3,4)=zz(i,j+1)
              ic(4)=icol(i,j+1)
              call bgnpol
              do n=1,4
                 if(ic(n).ne.-99999) then
                   call setcolor(ic(n))
                   call v2f(xy(1,n))
                 endif
              enddo
              call endpol
              tmin=min(xy(3,1),xy(3,2),xy(3,3),xy(3,4))
              tmax=max(xy(3,1),xy(3,2),xy(3,3),xy(3,4))
              lev1=tmin
              lev2=tmax
              num=0
              ips=0
              if(lev1.ne.-99999 .and. lev2.ne.-99999) then
                if(2*(lev1/2).ne.lev1) lev1=lev1-1
                if(2*(lev2/2).ne.lev2) lev2=lev2+1
                call setcolor(0)
                do lev=lev1,lev2,ispace
                  val=lev
                  do n=1,4
                    ii=n
                    jj=n+1
                    if(jj.gt.4) jj=1
                    IF(xy(3,ii).LT.VAL .AND. VAL.LE.xy(3,jj)) then
                      num=num+1
                      x2(1)=xy(1,jj)+(VAL-xy(3,jj))*(xy(1,ii)-xy(1,jj))/
     &                     (xy(3,ii)-xy(3,jj))
                      x2(2)=xy(2,jj)+(VAL-xy(3,jj))*(xy(2,ii)-xy(2,jj))/
     &                     (xy(3,ii)-xy(3,jj))
                      IF(ips.EQ.0) THEN
                        call bgnlin
                        CALL v2f(x2)
                        ips=1
                      ELSE
                        CALL v2f(x2)
                        call endlin
                        ips=0
                      ENDIF
                    endif
                    IF(xy(3,jj).LT.VAL .AND. VAL.LE.xy(3,ii)) then
                      num=num+1
                      x2(1)=xy(1,jj)+(VAL-xy(3,jj))*(xy(1,ii)-xy(1,jj))/
     &                     (xy(3,ii)-xy(3,jj))
                      x2(2)=xy(2,jj)+(VAL-xy(3,jj))*(xy(2,ii)-xy(2,jj))/
     &                     (xy(3,ii)-xy(3,jj))
                      IF(ips.EQ.0) THEN
                        call setcolor(0)
                        call bgnlin
                        CALL v2f(x2)
                        ips=1
                      ELSE
                        CALL v2f(x2)
                        call endlin
                        ips=0
                      ENDIF
                    endif
                  enddo
                enddo
              endif
          enddo
        enddo
        call swapbu
        call setcolor(0)
        call clear
        call rclear(xmin,xmax,ymin,ymax)
        if(iprint.ne.-1) call closeo
        if(iprint.ne.-1) call callob(ipr)
        if(iprint.eq.-1) call savescreen(ipr,850,425)

      endif
c      call closeo
c      call callob(ipr)

c
      goto 1
999   continue
      pause 'at 999'
      if(iprint.ne.-1) then
        iagain=1
123     continue
        call setcolor(0)
        call clear
        do i=1,ipr
          call setcolor(0)
          call clear
          call rclear(xmin,xmax,ymin,ymax)
          call callob(i)
          if(iagain.eq.-1) call savescreen(i,850,425)
        enddo
c this is for double buffering
c       call swapbu
        print *,'input 1 to do again,-1 with screendump'
        read(*,*) iagain
        if(abs(iagain).eq.1) goto 123
      endif
      call winclo(iwid)
      call gexit
      end
c==============================================================================
        subroutine initialize(ncolorr,xmin,xmax,ymin,ymax)

c #include <gl/fgl.h>
c #include <gl/fdevice.h>
#include "fshort.h"
      real xy1(2)
c     CALL PREFPO(430,430+850,599,599+425)
      isize=850
      jsize=425
      ihsize=getgde(0)-20
      ivsize=getgde(1)-40
      CALL PREFPO(1+ihsize-isize,ihsize,1+ivsize-jsize,ivsize)
        call foregr
        call keepas (4, 2)
        igid = winope ('My Window', 9)

        call double
        call RGBmod
        call zbuffe (.TRUE.)

        call qdevic (ESCKEY)
        call qdevic (KEYBD)
        call qdevic (MKEY)
        call qdevic (LEFTMO)
        call qdevic (MIDDLE)
        call qdevic (RIGHTM)
        call qdevic (MOUSEX)
        call qdevic (MOUSEY)
        call shadem (GOURAU)
c        call shadem (FLAT)
        call turnonmap(ncolors)
        call gconfi
        ncolorr=ncolors
c        call zclear
        call ortho2(xmin,xmax,ymin,ymax)
        call setcolor(0)
        call clear
        call rclear(xmin,xmax,ymin,ymax)
        call swapbu
        call setcolor(0)
        call clear
        call rclear(xmin,xmax,ymin,ymax)
        call swapbu
        return
        end
c==============================================================================
      subroutine rclear(xmin,xmax,ymin,ymax)
      real xy(2)
        call bgnpol
          xy(1)=xmin
          xy(2)=ymin
          call v2f(xy)
          xy(1)=xmax
          call v2f(xy)
          xy(2)=ymax
          call v2f(xy)
          xy(1)=xmin
          call v2f(xy)
        call endpol
      end
c----------------------------------
      FUNCTION IZSET(NCOLOR,ZZZ,ZF,ZD,IOFF)
      IMPLICIT REAL*4(A-H,O-Z)                                          
      IMAX=NCOLOR+IOFF
      IMIN=IOFF+1
C       IC=IMAX-NINT((ZZZ-ZF)/ZD-.5)
        ZIC=(ZZZ-ZF)/ZD+.5
        IC=NINT(ZIC)
C       WRITE(*,*) ZZZ,ZIC,IC
      IF(IC.LT.IMIN) IC=IMIN
      IF(IC.GT.IMAX) IC=IMAX
      IZSET=IC
      RETURN
      END
C===========================================
      FUNCTION HARDNESS(TTT)
      IMPLICIT REAL*4(A-H,O-Z)
      DATA RRR /8.314D0/
      DATA Q1 /60000.D0/
      DATA A1 /1.14D-5/
      DATA Q2 /139000.D0/
      DATA A2 /5.47D10/
      DATA EH /3.D0/
C
      IF(TTT.LT.0.) THEN
        TEMP=TTT+273.15D0
      ELSE
        TEMP=273.15D0
      ENDIF
      IF(TEMP.LE.263.15D0) THEN
        AAA=A1*EXP(-Q1/RRR/TEMP)
      ELSE
        AAA=A2*EXP(-Q2/RRR/TEMP)
      ENDIF
      BBB=EH*AAA*1D15
      BBB=1./BBB
      HARDNESS=BBB**(1.D0/3.D0)
      END
c-------------------------------------------------------------
      subroutine setcolor(i)
      common /colors/ ncolor,ired(0:300),igreen(0:300),iblue(0:300)
c       print *,i,ired(i),igreen(i),iblue(i)
      call RGBcol(ired(i),igreen(i),iblue(i))
      end
c-------------------------------------------------------------
      subroutine turnonmap(ncolors)
      common /colors/ ncolor,ired(0:300),igreen(0:300),iblue(0:300)
c read in color.map
      ncolor=0
      ired(0)=0
      igreen(0)=0
      iblue(0)=0
      do i=1,250
        read(72,*,end=9) ired(i),igreen(i),iblue(i)
        ncolor=i
      enddo      
9     continue
      if(ncolor.gt.0) then
        print *,'color map read',ncolor
        ncolors=ncolor
      else
        print *,'no color map found, link fort.72 to color.map'
        print *,'using internally generated colormap'
        ncolor=250
        ncolors=ncolor
c define color.map for red
        slope=real(0-240)/(40-1)
        do i=1,40
          ired(i)=240+(i-1)*slope
        enddo
        do i=40,140
          ired(i)=0
        enddo
        slope=real(255-0)/(180-140)
        do i=140,180
          ired(i)=0+(i-140)*slope
        enddo
        do i=180,230
          ired(i)=255
        enddo
        slope=real(200-255)/(250-230)
        do i=230,250
          ired(i)=255+(i-230)*slope
        enddo
c define color.map for green
        do i=1,40
          igreen(i)=0
        enddo
        slope=real(255-0)/(90-40)
        do i=40,90
          igreen(i)=0+(i-40)*slope
        enddo
        do i=90,190
          igreen(i)=255
        enddo
        slope=real(0-255)/(220-190)
        do i=190,220
          igreen(i)=255+(i-190)*slope
        enddo
        do i=220,250
          igreen(i)=0
        enddo
c define color.map for blue
        do i=1,90
          iblue(i)=255
        enddo
        slope=real(0-255)/(130-90)
        do i=90,130
          iblue(i)=255+(i-90)*slope
        enddo
        do i=130,250
          iblue(i)=0
        enddo
      endif
c      print *,' colors loaded:', ncolor
c      do i=0,ncolor
c        print *,i,ired(i),igreen(i),iblue(i)
c      enddo
      end
