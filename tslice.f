      program tslice
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM,MMAX=40,LMAX=281)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      dimension ilab(8)
      data ilab /1,30,61,92,123,154,185,215/
      dimension xp(lmax),yp(lmax,mmax),tp(lmax,mmax)
      dimension icol(nmax,mmax),xy(3,4),ic(4),x2(2)
      real*4 px(4),py(4)
      character*14 junk,hed
      IPR=1
      call grstrt(600,600)
1     continue
        read(37,end=999) hed
        read(37) np,llmax
        print *, np,llmax; pause
        read(37) (xp(i),i=1,np)
        read(37) ((yp(i,j),i=1,np),j=1,llmax)
        read(37) ((tp(i,j),i=1,np),j=1,llmax)
        ncolor=215
        xmax=-1e30
        xmin=1e30
        ymax=-1e30
        ymin=1e30
        tmax=-1e30
        tmin=1e30
        imax=0
        imin=ncolor
        do nn=1,np
          print *,nn,abs(yp(nn,llmax)-yp(nn,1))
          if(abs(yp(nn,llmax)-yp(nn,1)).gt.1) then
            print *,'dumping column',nn,xp(nn)
            DO I=1,llmax
              write(11,*) xp(nn),yp(nn,I),tp(nn,I)
            enddo
          endif
        enddo
        stop
        do nn=1,np
          xmax=MAX(xmax,xp(nn))
          xmin=MIN(xmin,xp(nn))
          DO I=1,MMAX
            ymax=MAX(ymax,yp(nn,I))
            ymin=MIN(ymin,yp(nn,I))
            tmax=MAX(tmax,tp(nn,I))
            tmin=MIN(tmin,tp(nn,I))
            write(11,*) xp(nn),yp(nn,I),tp(nn,I)
          enddo
c         tmax=0.0
c         tmin=-49.0
c         ymin=-6000.
c         ymax=6000.
          tdelt=(tmax-tmin)/(ncolor-1)
          do j=1,llmax
            if(tp(nn,j).eq.-99999.) then
              icol(nn,j)=-99999
            else
              icol(nn,j)=izset(ncolor,tp(nn,j),tmin,tdelt,0)
              imin=min(imin,icol(nn,j))
              imax=max(imax,icol(nn,j))
            endif
          enddo
        enddo
        ddx=(xmax-xmin)/10
        call newpag       
        call window(real(xmin),real(xmax),
     &              real(ymin),real(ymax))
        call linclr(1)
        call move(real(xmin),0.0)
        call draw(real(xmax),0.0)
        call move(1.0,real(yp(1,llmax)))
        do mm=2,np
          call draw(real(mm),real(yp(mm,llmax)))
        enddo
        call move(1.0,real(yp(1,1)))
        do mm=2,np
          call draw(real(mm),real(yp(mm,1)))
        enddo
          
        dx=(xmax-xmin-2*ddx)/ncolor
        dy=0.05*(ymax-ymin)
        ybox=ymax-dy
        call filpan(1,.false.)
        do i=1,ncolor
          xbox=xmin+(i-1)*dx+ddx
            px(1)=xbox
            py(1)=ybox

            px(2)=xbox+dx
            py(2)=ybox

            px(3)=xbox+dx
            py(3)=ybox+dy

            px(4)=xbox
            py(4)=ybox+dy

            call linclr1(i)
            call panel(4,px,py)
        enddo
        ybox=ymax-2*dy
        do ii=1,8
          i=ilab(ii)
          xbox=xmin+(i-1)*dx+ddx/2
          rnum=tmin+(i-1)*tdelt
          rnum=real(nint(rnum))
          icc=izset(ncolor,rnum,tmin,tdelt,0)    
          write(junk,'(f8.3)') rnum
c          write(*,'(i4,f8.3,i4)') i,rnum,izset(ncolor,rnum,zmin,zdelt,0)
          call linclr1(i)
          call move(real(xbox),real(ybox))
          call text(14,junk)
        enddo

        do nn=1,np-1
          px(1)=xp(nn)
          px(2)=xp(nn)
          px(3)=xp(nn+1)
          px(4)=xp(nn+1)
          do ll=1,MMAX-1
            py(1)=yp(nn,ll)
            py(2)=yp(nn,ll+1)
            py(3)=yp(nn+1,ll+1)
            py(4)=yp(nn+1,ll)
            icc=nint((icol(nn,ll)+icol(nn,ll+1)+
     &          icol(nn+1,ll+1)+icol(nn+1,ll))/4.)
            call linclr1(icc)
            call panel(4,px,py)
          enddo
        enddo
        call linclr(1)
        ybox=ymax-3*dy
        xbox=xmin+ddx/2
        call move(real(xbox),real(ybox))
        call text(14,hed)
        CALL savescreen(IPR,600,600)
        IPR=IPR+1
      goto 1
999   continue

      call grstop 
      end                                           
C===========================================
      FUNCTION IZSET(NCOLOR,ZZZ,ZF,ZD,IOFF)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      IMAX=NCOLOR+IOFF
      IMIN=IOFF+1
C       IC=IMAX-NINT((ZZZ-ZF)/ZD-.5d0)
        ZIC=(ZZZ-ZF)/ZD+.5d0
        IC=NINT(ZIC)
C       WRITE(*,*) ZZZ,ZIC,IC
      IF(IC.LT.IMIN) IC=IMIN
      IF(IC.GT.IMAX) IC=IMAX
      IZSET=IC
      RETURN
      END
