      subroutine eplate(itime,numnp,numel,x,y,kx,thick,kode,
     &                 delt,www,wrate,wmin,time,wwworig,fnet,
     &                 wwwtot)
c-----------------------------------------------------------------------
c ... elastic solver .............................................
c 4th order plate solver with one-time generation of stiffness and
c capacitance matrix. uses my sparse matrix storage for the static matrice
c and itpack sprse storage for the time-dependent modified matrices, and jcg
c iterative matrix solver. j fastook july 2000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)                                          
      include "parameter.h"
      parameter(nmax=maxnum,nz=27, nz1=nz+1, n3=3*nmax)
      parameter(nit=n3*nz,nmax1=n3+1,nmax3=n3*3)
      parameter(itmax=n3/10,ncg=4*itmax,nw=4*n3+ncg)
c ....arrays for itpack sparse storage...................
      dimension gk(nit),ja(nit),ia(nmax1)
      dimension iwksp(nmax3),wksp(nw),iwork(nit)
      dimension iparm(12),rparm(12)
      dimension thick(nmax),x(nmax),y(nmax),kx(nmax,4)
      dimension wwwtot(n3),kode(nmax)
      dimension www(n3),wrate(n3,2),wwworig(nmax)
      dimension kkx(nmax,12),ltemp(nmax,3)
      dimension xi(2,4),w(4)
      dimension ek(12,12),ec(12,12),ef(12)
      dimension gf(n3)
      dimension gk0(n3,nz),gc0(n3,nz)
      dimension ka(n3,nz1),kz(n3)
c     dimension wwwsave(nmax)
      character*80 hed
      common /pmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      common /iolist/ list

c ... timer stuff, for sgi only ...!
      real*4 tb(2)
c     real*4 tb(2),etime,dtime     !
c     external etime,dtime         !
c .................................!
      save ipass,kkx,wsave,nnode,xi,w,gk0,gc0,ka,kz
c     save wwwsave,iparm,rparm,iread
      data ipass /0/,iread /0/
1000  format(1x,a,t25,1pg13.6,g13.6)
1001  format(6g12.5)
c      write(7,1000) ' time beginning plate ',etime(tb),dtime(tb)
      if(ipass.eq.0) then
c ....  stuff you only do the first pass ..............
c        if(nmax.ne.nmax .or.n3.ne.m3) then
c          print *,'problems with nmax,mmax:',nmax,mmax
c          stop
c        endif
        call pconnect(nmax,numnp,numel,kx,kkx,ltemp)
        wsave=0.d0
        call psetint(xi,w)
        call psetmat
        nnode=3*numnp
        call pformgk(nmax,n3,nz,numel,nnode,
     &              x,y,kx,kkx,ek,ec,xi,w,
     &              gk0,gc0,ka,kz)
      endif
c ......................................................
c ....form stiffness and load ..............
      call pformgf(nmax,n3,numel,nnode,
     &              x,y,thick,wwwtot,kx,kkx,ef,xi,w,
     &              gf,fnet)
c ....time dependent case and varying load ...........
c ....and load into sparse storage
      call ptimedep(n3,nz,nnode,gk0,gc0,gk,gf,ka,kz,
     &              www,delt,nit,nmax1,ia,ja,iwork,nmax,kode)
c ......................................................
c.....dump matrix for examination ......................
      call odumpmat(n3,nz,nnode,kz,ka,gk0,gf,www,.false.)
      call ndumpmat(nit,gk,ja,nmax1,ia,nnode,gf,.false.)
c.......................................................
c     if(.true.) then
c.......solve equations with jordan conjugate-gradient itpack ....!
c        write(7,1000) ' time before jcg ',etime(tb),dtime(tb)!
        call dfault(iparm,rparm)
        if(delt.ne.0.0) then
          iparm(1)=150  ! max number of iteration
        else
          iparm(1)=1000 ! max number of iteration
        endif
        rparm(1)=1d-6   ! zeta, stopping criteria
c       iparm(2)=2      ! level of output (-1:none)
c       iparm(4)=7      ! output unit number
        iparm(5)=1      ! nonsymmetric matrix (1)
        iparm(6)=0      ! non-adaptive (0) (1 doesnt work)
c       iparm(10)=1     ! removes large diagonal entries (doesnt work)
c        do i=1,12
c          print *,i,iparm(i),rparm(i)
c        enddo
        call jcg(nnode,ia,ja,gk,gf,www,
     &         iwksp,nw,wksp,iparm,rparm,ier)
c         do i=1,12
c           print *,i,iparm(i),rparm(i)
c         enddo
c        if(iotogg) then
c          write(list(ipage+1),*) ' relative error=',rparm(1),
c     &            ' in iterations = ',iparm(1)
c          write(list(ipage+2),*) rparm(11),rparm(12)
c          ipage=ipage+2
c        endif
        if(ier.ne.0) then
          call ndumpmat(nit,gk,ja,nmax1,ia,nnode,gf,.true.)
          print *,'jcg error:',ier
c          print '(1x,i3,1pg13.6)',(i,rparm(i),i=1,12)
          pause
        endif
c        write(7,1000) ' time after jcg ',etime(tb),dtime(tb) !
c ................................................................!
c     else
c.......solve equations with conjugate-gradient ..................!
c       write(7,1000) ' time before conjug ',etime(tb),dtime(tb)  !
c       call conjug(n3,nz,nnode,1.d-6,gk,ka,gf,www)              !
c       write(7,1000) ' time after conjug ',etime(tb),dtime(tb)   !
c ................................................................!
c     endif
      if(delt.gt.0.) then
        do i=1,nnode
          wrate(i,1)=(www(i)-wrate(i,2))/delt
        enddo
      endif
c$doacross local(i)
      do i=1,nnode
        wrate(i,2)=www(i)
      enddo
      wmin=1d30
      wmax=-1d30
      wtot=0.d0
      do i=1,nnode,3
        wmin=min(wmin,www(i))
        wmax=max(wmax,www(i))
        wtot=wtot+www(i)
c        write(*,*) (www(i+j),j=0,2)
      enddo
      thtot=0.d0
      do i=1,numnp
        thtot=thtot+thick(i)
      enddo
      if(iread.eq.0 .and. ipass.eq.0) then
        ii=1
        do i=1,numnp
          wwworig(i)=www(ii)
          ii=ii+3
        enddo
      endif
      if(iotogg) then
        write(list(ipage+1),*) 
     &       '*********** elastic **********************'
        write(list(ipage+2),1001) time,thtot,
     &              -thtot/wtot/rockice,
     &              wmin,wmax,-1000*(wsave-wmin)/delt
c        write(list(ipage+3),*) 
c     &       '*******************************************'
c        ipage=ipage+3
        ipage=ipage+2
      endif
      write(92,*) time,wmin
      wsave=wmin
      ipass=1
      end
c=============================================================
      subroutine psetint(xi,w)
      implicit real*8(a-h,o-z)
      dimension xi(2,4),w(4)
c ...............................
      xi(1,1)=1.d0/sqrt(3.d0)
      xi(2,1)=xi(1,1)
      xi(1,2)=-xi(1,1)
      xi(2,2)=xi(1,1)
      xi(1,3)=xi(1,1)
      xi(2,3)=-xi(1,1)
      xi(1,4)=-xi(1,1)
      xi(2,4)=-xi(1,1)
      w(1)=1.d0
      w(2)=1.d0
      w(3)=1.d0
      w(4)=1.d0
c ...............................
      end
c=============================================================
      subroutine pelemek(xy,n,ekb,ec,nl,xi,w)
      implicit real*8(a-h,o-z)
      dimension xy(2,n),ekb(12,12),eks(12,12),ekh(12,12)
      dimension ec(12,12)
      dimension dpsix(4),dpsiy(4)
      dimension psi(4),dpsi(4,2),xs(2)
      dimension bb(3,12),bs(2,12),db(3,3),ds(2,2)
      dimension btdb(3,12)
      dimension xi(2,4),w(4)
      common /pmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c.....initialize element arrays
c$doacross local(i,j)
      do i=1,12
        do j=1,12
          ekb(i,j)=0.0d0
          eks(i,j)=0.0d0
          ekh(i,j)=0.0d0
          ec(i,j)=0.0d0
        enddo
      enddo
c ... form d-matrices ...............
      ttt3=ttt**3/12.d0
      db(1,1)=ttt3*(2.d0*rmu+rlambda)
      db(1,2)=ttt3*rlambda
      db(1,3)=0
      db(2,1)=ttt3*rlambda
      db(2,2)=ttt3*(2.d0*rmu+rlambda)
      db(2,3)=0
      db(3,1)=0
      db(3,2)=0
      db(3,3)=ttt3*rmu
      ds(1,1)=ttt*rmu
      ds(1,2)=0
      ds(2,1)=0
      ds(2,2)=ttt*rmu
c.....begin 2x2 integration loop
      do l=1,nl
        xs(1)=xi(1,l)
        xs(2)=xi(2,l)
        call pgenshape(n,xs,xy,psi,dpsi,detj,dpsix,dpsiy)
        call ploadb(psi,dpsix,dpsiy,bb,bs)
c
c.......accumulate integration point value of integrals
        fac=detj*w(l)
c  .... ratios of density of rock and ice for overburden load
        xb=rockice*rhog
c ..... form capacitance matrix .....................
        do i=1,4
          ip1=3*i-2
          ip2=3*i-1
          ip3=3*i
          do j=1,4
            jp1=3*j-2
            jp2=3*j-1
            jp3=3*j
            ekh(ip1,jp1)=ekh(ip1,jp1)+fac*xb*psi(i)*psi(j)
            term=fac*ctime*psi(i)*psi(j)
            ec(ip1,jp1)=ec(ip1,jp1)+term
            ec(ip1,jp2)=ec(ip1,jp2)+term
            ec(ip1,jp3)=ec(ip1,jp3)+term
            ec(ip2,jp1)=ec(ip2,jp1)+term
            ec(ip2,jp2)=ec(ip2,jp2)+term
            ec(ip2,jp3)=ec(ip2,jp3)+term
            ec(ip3,jp1)=ec(ip3,jp1)+term
            ec(ip3,jp2)=ec(ip3,jp2)+term
            ec(ip3,jp3)=ec(ip3,jp3)+term
          enddo
        enddo
c ..... form kb stiffness matrix ..............
c ..... ok, now here goes the matrix multiplication that 
c ..... generates the bt d b ala book...
c ..... first bbt*db (12x3)*(3x3)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,3
            sum=0.d0
            do k=1,3
              sum=sum+bb(k,i)*db(j,k)
            enddo
            btdb(j,i)=sum
          enddo
        enddo
c then (bbt*db)*bb (12x12)*(3x12)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,12
            sum=0.d0
            do k=1,3
              sum=sum+btdb(k,i)*bb(k,j)
            enddo
            ekb(i,j)=ekb(i,j)+fac*sum
          enddo
        enddo
      enddo
c end of 2x2 integration loop
c
c.....begin 1x1 reduced integration loop
      do l=1,1
        xs(1)=0.d0
        xs(2)=0.d0
        call pgenshape(n,xs,xy,psi,dpsi,detj,dpsix,dpsiy)
        call ploadb(psi,dpsix,dpsiy,bb,bs)
c
c.......accumulate integration point value of integrals
        fac=detj*4.d0
c ..... form ks stiffness matrix ........................
c ..... ok, now here goes the matrix multiplication that 
c ..... generates bst*ds*bs ala book...
c ..... first bst*ds (12x2)*(2x2)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,2
            sum=0.d0
            do k=1,2
              sum=sum+bs(k,i)*ds(j,k)
            enddo
            btdb(j,i)=sum
          enddo
        enddo
c then (bbt*db)*bb (12x12)*(2x12)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,12
            sum=0.d0
            do k=1,2
              sum=sum+btdb(k,i)*bs(k,j)
            enddo
            eks(i,j)=eks(i,j)+fac*sum
          enddo
        enddo
      enddo
c end of 1x1 integration loop
c
c ... combine kb and ks stiffness matrices .......
      do i=1,12
        do j=1,12
c          print '(2i3,1p3g13.6)',i,j,ekb(i,j),eks(i,j),ekh(i,j)
c          ekb(i,j)=ekb(i,j)+eks(i,j)+ekh(i,j)
          ekb(i,j)=ekb(i,j)+eks(i,j)
        enddo
      enddo
c      pause
c      if(.false.) then
c        do i=1,12
c          print 1003,(ekb(i,j)/ekb(i,i),j=1,12)
c        enddo
c        print *,'---------------------------'
c        do i=1,12
c          print 1003,(ec(i,j)/ec(i,i),j=1,12)
c        enddo
c        pause
c      endif
1003  format(1x,1p6g13.6)
      return
1000  format(1x,'bad jacobian',e10.3)
1001  format(a,/,(2i5,3x,1pe13.6))
1002  format(a,/,(1p2e13.6))
      end
c=============================================================
      subroutine pelemgf(ethick,xy,n,ef,nl,xi,w)
      implicit real*8(a-h,o-z)
      dimension xy(2,n)
      dimension ef(12)
      dimension dpsix(4),dpsiy(4)
      dimension psi(4),dpsi(4,2),xs(2)
      dimension xi(2,4),w(4)
c      common /pmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c.....initialize element arrays
c$doacross local(i)
      do i=1,12
        ef(i)=0.0d0
      enddo
c.....begin 2x2 integration loop
      do l=1,nl
        xs(1)=xi(1,l)
        xs(2)=xi(2,l)
        call pgenshape(n,xs,xy,psi,dpsi,detj,dpsix,dpsiy)
c
c.......accumulate integration point value of integrals
        fac=detj*w(l)
        xf=-ethick
c ..... form capacitance matrix .....................
        do i=1,4
          ip1=3*i-2
          ef(ip1)=ef(ip1)+xf*psi(i)*fac
        enddo
      enddo
      return
1000  format(1x,'bad jacobian',e10.3)
1001  format(a,/,(2i5,3x,1pe13.6))
1002  format(a,/,(1p2e13.6))
      end
c=============================================================
      subroutine passmbgk(n3,nz,gk0,gc0,ka,kz,ek,ec,n,node)
      implicit real*8(a-h,o-z)
      dimension gk0(n3,nz),gc0(n3,nz)
      dimension ka(n3,nz+1),kz(n3)
      dimension ek(12,12),ec(12,12),node(n)
c........................
c||||||||||||||||||||||||
c     print *,'in passmb'
      do l=1,n
        i=node(l)
        do m=1,n
          j=node(m)
c.........assemble global stiffness matrix gk
            if(i.eq.j) then
              gk0(i,1)=gk0(i,1)+ek(l,m)
              gc0(i,1)=gc0(i,1)+ec(l,m)
              ka(i,1)=i
            else
              do k=2,kz(i)
                if(ka(i,k).eq.j) then
                  gk0(i,k)=gk0(i,k)+ek(l,m)
c fix fix fix in other versions...
c                 gc0(i,k)=gc0(i,k)+ek(l,m)
                  gc0(i,k)=gc0(i,k)+ec(l,m)
c fix fix fix in other versions...
                  goto 99
                endif
              enddo
              kz(i)=kz(i)+1
              gk0(i,kz(i))=ek(l,m)
              gc0(i,kz(i))=ec(l,m)
              ka(i,kz(i))=j
            endif
99          continue
        enddo
      enddo
c||||||||||||||||||||||||
c........................
      end
c=============================================================
      subroutine aplybc(n3,gf,nmax,kode,nnode)
      implicit real*8(a-h,o-z)
      dimension gf(n3),kode(nmax)
      data big /1d30/
      do i=1,nnode
        n=1+(i-1)/3
        if(kode(n).eq.1) gf(i)=0.d0
      enddo
      end
c=============================================================
      subroutine passmbgf(n3,gf,ef,n,node)
      implicit real*8(a-h,o-z)
      dimension gf(n3)
      dimension ef(12),node(n)
c........................
c||||||||||||||||||||||||
c     print *,'in passmbgf'
      do l=1,n
        i=node(l)
c.......assemble global vector gf
        gf(i)=gf(i)+ef(l)
      enddo
c||||||||||||||||||||||||
c........................
      end
c=============================================================
      subroutine pformgk(nmax,n3,nz,numel,nnode,
     &                  x,y,kx,kkx,ek,ec,xi,w,
     &                  gk0,gc0,ka,kz)
      implicit real*8(a-h,o-z)
      dimension x(nmax),y(nmax),kx(nmax,4),kkx(nmax,12)
      dimension gk0(n3,nz),gc0(n3,nz)
      dimension ka(n3,nz+1),kz(n3)
      dimension xi(2,4),w(4)
      dimension ek(12,12),ec(12,12)
      dimension lm(4),xy(2,4),llm(12)
c      common /pmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c$doacross local(i,j)
      do i=1,nnode
        kz(i)=1
        ka(i,1)=i
        ka(i,nz+1)=0
        do j=1,nz
          ka(i,j)=0
          gk0(i,j)=0.0d0
          gc0(i,j)=0.0d0
        enddo
      enddo
      do nel=1,numel
        do l=1,4
          lm(l)=kx(nel,l)
          xy(1,l)=x(lm(l))
          xy(2,l)=y(lm(l))
        enddo
c$doacross local(l)
        do l=1,12
          llm(l)=kkx(nel,l)
        enddo
        call pelemek(xy,4,ek,ec,4,xi,w)
c..................................
        if(.false.) then
          write(7,*) nel
          do i=1,12
            write(7,*) (ek(i,j),j=1,12)
          enddo
c          pause
        endif
c..................................
        call passmbgk(n3,nz,gk0,gc0,ka,kz,ek,ec,12,llm)
c..................................
c        do i=1,nnode
c          print 1000,(gk0(i,j),j=1,nnode)
c        enddo
c        pause
c..................................
      enddo
c$doacross local(i)
      do i=1,nnode
c        print *,kz(i),ka(i,nz+1)
        ka(i,nz+1)=kz(i)
      enddo
c|||||||||||||||||||||||||||||||||||
c...................................
1000  format(1x,10f8.3)
      end
c=============================================================
      subroutine pformgf(nmax,n3,numel,nnode,
     &              x,y,thick,www,kx,kkx,ef,xi,w,
     &              gf,fnet)
      implicit real*8(a-h,o-z)
      dimension x(nmax),y(nmax),kx(nmax,4),kkx(nmax,12)
      dimension gf(n3),www(n3)
      dimension thick(nmax)
      dimension xi(2,4),w(4)
      dimension ef(12)
      dimension lm(4),xy(2,4),llm(12)
      common /pmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c$doacross local(i)
      do i=1,nnode
        gf(i)=0.d0
      enddo
      fnet=0.d0
      do nel=1,numel
        ethick=0.d0
        ewww=0.d0
        do l=1,4
          lm(l)=kx(nel,l)
          xy(1,l)=x(lm(l))
          xy(2,l)=y(lm(l))
          ethick=ethick+thick(lm(l))
          ewww=ewww+www((lm(l)-1)*3+1)
        enddo
c$doacross local(l)
        do l=1,12
          llm(l)=kkx(nel,l)
        enddo
        ethick=ethick/dble(4)
        ewww=ewww/dble(4)
c         print '(i5,5g13.6)',nel,real(ethick),real(ewww),
c     &          real(rhog*ethick),real(rockice*rhog*ewww),
c     &          real(rhog*ethick+rockice*rhog*ewww)
        ethick=rhog*ethick
        ewww=rockice*rhog*ewww
c       print *,nel,real(ethick),real(ewww),real(ethick+ewww)
        ethick=ethick+ewww
c ...... EXPERIMENTAL ..........
c        ethick=max(0.d0,ethick)
c ...... EXPERIMENTAL ..........
        fnet=fnet+ethick
        call pelemgf(ethick,xy,4,ef,4,xi,w)
c..................................
        if(.false.) then
          write(7,*) nel
          do i=1,12
            write(7,*) ef(i)
          enddo
c          pause
        endif
c..................................
        call passmbgf(n3,gf,ef,12,llm)
c..................................
c        do i=1,nnode
c          print 1000,gf(i)
c        enddo
c        pause
c..................................
      enddo
c|||||||||||||||||||||||||||||||||||
c...................................
1000  format(1x,10f8.3)
      end
c=============================================================
      subroutine pshape(xs,psi,dpsi)
      implicit real*8(a-h,o-z)
      dimension xs(2),psi(4),dpsi(4,2)
c...................................
c|||||||||||||||||||||||||||||||||||
      eta=xs(1)
      rnu=xs(2)
      px0=1.d0-eta
      px1=1.d0+eta
      py0=1.d0-rnu
      py1=1.d0+rnu
      psi(1)=px0*py0/4.d0
      psi(2)=px1*py0/4.d0
      psi(3)=px1*py1/4.d0
      psi(4)=px0*py1/4.d0
      dpsi(1,1)=-py0/4.d0
      dpsi(2,1)=-dpsi(1,1)
      dpsi(3,1)=py1/4.d0
      dpsi(4,1)=-dpsi(3,1)
      dpsi(1,2)=-px0/4.d0
      dpsi(2,2)=-px1/4.d0
      dpsi(3,2)=-dpsi(2,2)
      dpsi(4,2)=-dpsi(1,2)
c|||||||||||||||||||||||||||||||||||
c...................................
      end
c=============================================================
      subroutine ploadb(psi,dpsix,dpsiy,bb,bs)
      implicit real*8(a-h,o-z)
      dimension psi(4),dpsix(4),dpsiy(4)
      dimension bb(3,12),bs(2,12)
      do ia=1,4
        ip1=3*ia-2
        ip2=ip1+1
        ip3=ip1+2
        bb(1,ip1)=0.d0
        bb(2,ip1)=0.d0
        bb(3,ip1)=0.d0
        bb(1,ip2)=dpsix(ia)
        bb(2,ip2)=0.d0
        bb(3,ip2)=dpsiy(ia)
        bb(1,ip3)=0.d0
        bb(2,ip3)=dpsiy(ia)
        bb(3,ip3)=dpsix(ia)
        bs(1,ip1)=dpsix(ia)
        bs(2,ip1)=dpsiy(ia)
        bs(1,ip2)=-psi(ia)
        bs(2,ip2)=0.d0
        bs(1,ip3)=0.d0
        bs(2,ip3)=-psi(ia)
      enddo
      end
c=============================================================
      subroutine pconnect(nmax,numnp,numel,kx,kkx,ltemp)
      implicit real*8(a-h,o-z)
      dimension kx(nmax,4),kkx(nmax,12),ltemp(nmax,3)
      ic=1
      do i=1,numnp
        ltemp(i,1)=ic
        ltemp(i,2)=ic+1
        ltemp(i,3)=ic+2
        ic=ic+3
c        print *,i,(ltemp(i,j),j=1,3)
      enddo
c$doacross local(i,j,ip1,ip2,ip3)
      do i=1,numel
        do j=1,4
          ip1=3*j-2
          ip2=ip1+1
          ip3=ip1+2
          kkx(i,ip1)=ltemp(kx(i,j),1)
          kkx(i,ip2)=ltemp(kx(i,j),2)
          kkx(i,ip3)=ltemp(kx(i,j),3)
        enddo
c        print 1000,i,(kx(i,j),j=1,4)
c        print 1000,i,(kkx(i,j),j=1,12)
c        pause
      enddo
1000  format(1x,13i5)
      end
c=============================================================
      subroutine poutstf(nnode,www,wrate,numnp,thick,numcol,numlev,
     &           time,wdiff)
      implicit real*8(a-h,o-z)
      include "parameter.h"
      parameter(nmax=maxnum,n3=3*nmax)
      dimension www(nnode),thick(numnp),tmp(n3),wrate(nnode)
      dimension wdiff(numnp)
      character*80 hed
c ... turn it off ...
c      return
c
c ... depression ...
      do i=1,nnode
        tmp(i)=www(i)
        if(www(i).gt.0) tmp(i)=tmp(i)*1
      enddo
      if(.false.) then
        write(81,1000)
        write(81,1001) numcol,numlev
        write(81,1002) 0,numcol*1000,0,numlev*1000
        write(81,1003) 'depression'
        write(81,1008) time
        write(81,1004)
        write(81,1005)
        write(81,1006) (tmp(i),i=1,nnode,3)
        write(81,1007)
      endif
c ... load ...
      if(.false.) then
        write(80,1000)
        write(80,1001) numcol,numlev
        write(80,1002) 0,numcol*1000,0,numlev*1000
        write(80,1003) 'load'
        write(80,1008) time
        write(80,1004)
        write(80,1005)
        write(80,1006) (thick(i),i=1,numnp)
        write(80,1007)
      endif
c ... rate ... mm/yr
c$doacross local(i)
      do i=1,nnode
        tmp(i)=wrate(i)*1000.d0
      enddo
      if(.false.) then
        write(82,1000)
        write(82,1001) numcol,numlev
        write(82,1002) 0,numcol*1000,0,numlev*1000
        write(82,1003) 'rate'
        write(82,1008) time
        write(82,1004)
        write(82,1005)
        write(82,1006) (tmp(i),i=1,nnode,3)
        write(82,1007)
      endif
c ... depression difference ...
      if(.false.) then
        write(83,1000)
        write(83,1001) numcol,numlev
        write(83,1002) 0,numcol*1000,0,numlev*1000
        write(83,1003) 'difference'
        write(83,1008) time
        write(83,1004)
        write(83,1005)
        write(83,1006) (wdiff(i),i=1,numnp)
        write(83,1007)
      endif
      write(hed,*) ' time=',nint(time)
      write(88) hed
      write(88) (www(i),i=1,nnode,3)
      write(88) (thick(i),i=1,numnp)
      write(88) (1000*wrate(i),i=1,nnode,3)
c ... 
1000  format(1x,'rank 2')
1001  format(1x,'dimensions',2i7)
1002  format(1x,'bounds',4i13)
1003  format(1x,'name ',a)
1004  format(1x,'scalar')
1005  format(1x,'data')
1006  format(1x,1p5g14.6)      
1007  format(1x,'end')
1008  format(1x,'time',g13.6)
1033  format(1x,'name load')
      end
c=============================================================
      subroutine ptimedep(n3,nz,nnode,gk0,gc0,gk,gf,ka,kz,
     &                    www,delt,nit,nmax1,ia,ja,iwork,nmax,kode)
      implicit real*8(a-h,o-z)
c ....arrays for itpack sparse storage...................
      dimension gk(nit),ja(nit),ia(nmax1)
      dimension iwork(nit)
      dimension gf(n3),www(n3),kode(nmax)
      dimension gk0(n3,nz),gc0(n3,nz)
      dimension ka(n3,nz+1),kz(n3)
      data big /1d30/
      if(delt.gt.0.d0) then
        call sbini(nnode,nit,ia,ja,gk,iwork)
c ..... form modified stiffness and load .......
        delt1=1.d0/delt    
        do i = 1,nnode
          n=1+(i-1)/3
          if(kode(n).eq.1) then
            call sbsij(nnode,nit,ia,ja,gk,iwork,i,i,big,
     &                   0,0,6,ier)
            gf(i)=0.d0
          else
            do j = 1,kz(i)
              jg=ka(i,j)
c             print *,i,j,gf(i),gc0(i,j)*www(jg)*delt1
c             if(.false.) then ! full matrix...
                gf(i)=gf(i)+gc0(i,j)*www(jg)*delt1
                gkij=gk0(i,j)+gc0(i,j)*delt1
                call sbsij(nnode,nit,ia,ja,gk,iwork,i,jg,gkij,
     &                     0,0,6,ier)
c              else ! lumped
c                gf(i)=gf(i)+gc0(i,j)*www(i)*delt1
c                gkij=gk0(i,1)+gc0(i,j)*delt1
c                call sbsij(nnode,nit,ia,ja,gk,iwork,i,i,gkij,
c     &                     1,0,6,ier)
c              endif
            enddo 
          endif
        enddo
        call sbend(nnode,nit,ia,ja,gk,iwork)
      else
        call sbini(nnode,nit,ia,ja,gk,iwork)
        do i = 1,nnode
c          n=1+(i-1)/3
c          if(kode(n).eq.1) then
c            call sbsij(nnode,nit,ia,ja,gk,iwork,i,i,big,
c     &                   0,0,6,ier)
c            gf(i)=0.d0
c          else
            do j = 1,kz(i)
              jg=ka(i,j)
              gkij=gk0(i,j)
              call sbsij(nnode,nit,ia,ja,gk,iwork,i,jg,gkij,0,1,6,ier)
            enddo 
c          endif
        enddo
        call sbend(nnode,nit,ia,ja,gk,iwork)
      endif
      end
c=============================================================
      subroutine pgauseid(nmax,nz,n,eps,aa,ka,b,x)
      implicit real*8(a-h,o-z)                                          
      dimension aa(nmax,nz),ka(nmax,nz+1),b(nmax),x(nmax)
      real*8 sum
c      data eps /1d-6/, itmax /100/, tol /0.d0/
      data itmax /100/, tol /0.d0/
c ... a diffrent first guess ...
c     if(.false.) then
c       do i=1,n
c         x(i)=b(i)/aa(i,1)
c       enddo
c     endif
c ... ..........................
      xmax=-1d30
      do i=1,n
        jmax=ka(i,nz+1)
        sum=b(i)
        do j=2,jmax
          sum=sum-aa(i,j)*x(ka(i,j))
        enddo
        xnew=sum/aa(i,1)
        xmax=max(xmax,abs(x(i)-xnew))
        x(i)=xnew
      enddo
      errorg=presid(nmax,nz,n,aa,ka,b,x)
c     print 1001,   0,errorg
      write(7,1001) 0,errorg
c ... rare but possible return
      if(errorg.eq.0.) return
      do iter=1,itmax
      xmax=-1d30
        do i=1,n
          jmax=ka(i,nz+1)
          sum=b(i)
          do j=2,jmax
            sum=sum-aa(i,j)*x(ka(i,j))
          enddo
          xnew=sum/aa(i,1)
          xmax=max(xmax,abs(x(i)-xnew))
          x(i)=xnew
        enddo
c       if(.true.) then
c         print 1001,iter,xmax
c         print 1000,(x(i),i=1,n)
c         write(7,1001) iter,xmax
c         write(7,1000) (x(i),i=1,n)
c       endif
        error=presid(nmax,nz,n,aa,ka,b,x)
        ratio=error/errorg
c       print 1001,   iter,error,ratio,xmax
        write(7,1001) iter,error,ratio,xmax
c ..... a normal return when converged ...
        if(abs(ratio).lt.eps) then
c         print 1002,   'a:converged',iter,error,ratio,xmax
          write(7,1002) 'a:converged',iter,error,ratio,xmax
          return
        endif
c ..... another normal return when converged ...
        if(xmax.lt.tol) then
c         print 1002,   'b:converged',iter,error,ratio,xmax
          write(7,1002) 'b:converged',iter,error,ratio,xmax
          return
        endif
      enddo
      print *,   'didnot converge in ',itmax,''
      write(7,*) 'didnot converge in ',itmax
c     pause
1000  format(1x,'pgauseid:',(1x,t10,5g14.6))
1001  format(1x,'pgauseid:',t21,i5,3g13.6)
1002  format(1x,'pgauseid:',a,i5,3g13.6)
      end
c=======================================
      function presid(nmax,nz,n,aa,ka,b,x)
      implicit real*8(a-h,o-z)                                          
      dimension aa(nmax,nz),ka(nmax,nz+1),b(nmax),x(nmax)
      real*8 sumsq,sum
      sumsq=0.d0
      do i=1,n
        jmax=ka(i,nz+1)
        sum=0.d0
        do j=1,jmax
          sum=sum+aa(i,j)*x(ka(i,j))
        enddo
        sum=sum-b(i)
        sumsq=sumsq+sum**2
      enddo
      presid=sumsq
      end
c=======================================
      subroutine pgenshape(n,xs,xy,psi,dpsi,detj,dpsix,dpsiy)
      implicit real*8(a-h,o-z)                                          
      dimension xy(2,n)
      dimension dpsix(4),dpsiy(4),dxds(2,2),dsdx(2,2)
      dimension psi(4),dpsi(4,2),xs(2)
        call pshape(xs,psi,dpsi)
c.......calculate dsds...equation(5.3.6)
        do i=1,2
          do j=1,2
            dxds(i,j)=0.0d0
            do k=1,4
              dxds(i,j)=dxds(i,j)+dpsi(k,j)*xy(i,k)
            enddo
          enddo
        enddo
c.......calculate dsdx...equation(5.2.7)
        detj=dxds(1,1)*dxds(2,2)-dxds(1,2)*dxds(2,1)
        if(detj.le.0.0) then
          write(*,1000) detj
          write(*,1001) 'dxds',((i,j,dxds(i,j),i=1,2),j=1,2)
          write(*,1001) 'dpsi',((i,j,dpsi(j,i),i=1,2),j=1,4)
          write(*,1002) 'xy',((xy(i,j),i=1,2),j=1,n)
          stop
        endif
        dsdx(1,1)=dxds(2,2)/detj    
        dsdx(2,2)=dxds(1,1)/detj    
        dsdx(1,2)=-dxds(1,2)/detj    
        dsdx(2,1)=-dxds(2,1)/detj  
c.......calculate d(psi)/dx...equation(5.3.5)
        do i=1,n
          dpsix(i)=dpsi(i,1)*dsdx(1,1)+dpsi(i,2)*dsdx(2,1)  
          dpsiy(i)=dpsi(i,1)*dsdx(1,2)+dpsi(i,2)*dsdx(2,2) 
        enddo
1000  format(1x,'bad jacobian',e10.3)
1001  format(a,/,(2i5,3x,1pe13.6))
1002  format(a,/,(1p2e13.6))
      end
c=======================================
      subroutine ndumpmat(nit,gk,ja,nmax1,ia,nnode,gf,disp)
      implicit real*8(a-h,o-z)                                          
      logical disp
      dimension gk(nit),ja(nit),ia(nmax1),gf(nnode)
      if(disp) then
        do l=1,nnode
          sum=0.0d0
          do j=ia(l)+1,ia(l+1)-1
            sum=sum+gk(j)
          enddo
          gdiag=gk(ia(l))
          write(13,30) 'row',l,' gf=',gf(l),gdiag,sum,gdiag/sum
          write(13,10) (gk(j),j=ia(l),ia(l+1)-1)
          write(13,10) (gk(j)/gdiag,j=ia(l),ia(l+1)-1)
          write(13,20) (ja(j),j=ia(l),ia(l+1)-1)
        enddo
        pause
      endif
10    format(5x,1p6g13.6)
20    format(1x,6i13)
30    format(1x,a,i6,a,1p4g13.6)
      end
c=======================================
      subroutine odumpmat(n3,nz,nnode,kz,ka,gk,gf,www,disp)
      implicit real*8(a-h,o-z)                                          
      logical disp
      dimension gk(n3,nz),gf(n3),www(n3)
      dimension ka(n3,nz+1),kz(n3)
        if(.false.) then
          do i=1,nnode
            write(7,*) 'diagonal='
            write(7,*) kz(i),i,gk(i,1)
            sum=0
            do j=2,kz(i)
              sum=sum+abs(gk(i,j))
              jj=ka(i,j)
              write(7,*) i,jj,gk(i,j),gk(i,j)/gk(i,1)
            enddo
            if(sum/gk(i,1).gt.1.1) then
              write(7,*) 'eq:',i,' not diag dom.',gk(i,1)/sum
c              pause
            endif
c            pause
            write(7,*) i,gf(i)
          enddo
c          pause
        endif
        if(disp) then
          write(13,*) nnode
          do i=1,nnode
            write(13,*) (gk(i,j),j=1,nz)
            write(13,*) (ka(i,j),j=1,nz+1)
            write(13,*) gf(i),www(i)
          enddo
        endif
      end
c
c=======================================================================
c
c
c
      subroutine vplate(itime,numnp,numel,x,y,kx,thick,kode,
     &                 delt,www,wrate,wmin,time,wwworig,fnet,
     &                 wwwtot)
c-----------------------------------------------------------------------
c 4th order viscous plate solver with one-time generation of stiffness
c matrix. uses my sparse matrix storage for the static matrice
c and itpack sprse storage for the time-dependent modified matrices, and c             ipage=ipage+2

c iterative matrix solver. j fastook july, 1999
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)                                          
      include "parameter.h"
      parameter(nmax=maxnum,nz=27, nz1=nz+1, n3=3*nmax)
      parameter(nit=n3*nz,nmax1=n3+1,nmax3=n3*3)
      parameter(itmax=n3/10,ncg=4*itmax,nw=4*n3+ncg)
c ....arrays for itpack sparse storage...................
      dimension gk(nit),ja(nit),ia(nmax1)
      dimension iwksp(nmax3),wksp(nw),iwork(nit)
      dimension iparm(12),rparm(12)
      dimension thick(nmax),www(n3),x(nmax),y(nmax),kx(nmax,4)
      dimension wwwtot(n3),kode(nmax)
      dimension wrate(n3,2),wwworig(nmax)
      dimension kkx(nmax,12),ltemp(nmax,3)
      dimension xi(2,4),w(4)
      dimension ek(12,12),ef(12)
      dimension gf(n3)
      dimension gk0(n3,nz)
      dimension ka(n3,nz1),kz(n3)
      dimension wwwtmp(n3)
      character*80 hed
      common /vmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      common /iolist/ list

c ... timer stuff, for sgi only ...!
      real*4 tb(2)
c     real*4 tb(2),etime,dtime     !
c     external etime,dtime         !
c .................................!
      save ipass,kkx,wsave,nnode,xi,w,gk0,ka,kz
      save gk,ia,ja
c     save wwwsave,iparm,rparm,iread
      data ipass /0/,iread /0/
1000  format(1x,a,t25,1pg13.6,g13.6)
1001  format(6g12.5)
c      write(7,1000) ' time beginning plate ',etime(tb),dtime(tb)
      if(ipass.eq.0) then
c ....  stuff you only do the first pass ..............
c        if(nmax.ne.mmax .or.n3.ne.m3) then
c          print *,'problems with nmax,mmax:',nmax,mmax
c          stop
c        endif
        call pconnect(nmax,numnp,numel,kx,kkx,ltemp)
        wsave=0.d0
        call psetint(xi,w)
        call vsetmat
        nnode=3*numnp
        call vformgk(nmax,n3,nz,numel,nnode,
     &              x,y,kx,kkx,ek,xi,w,
     &              gk0,ka,kz)
        call vloadgk(n3,nz,nnode,gk0,gk,ka,kz,
     &                    nit,nmax1,ia,ja,iwork,nmax,kode)
      endif
      do i=1,nnode
        wwwtmp(i)=www(i)
      enddo
c ... 2-step method ...
      do ii=1,2
c ......................................................
c ......form stiffness and load ..............
        if(.true.) then
c ......  use own solution ....
          call vformgf(nmax,n3,numel,nnode,
     &              x,y,thick,wwwtmp,kx,kkx,ef,xi,w,
     &              gf,fnet,kode)
        else
c ......  use other solution ....
          call vformgf(nmax,n3,numel,nnode,
     &              x,y,thick,wwwtot,kx,kkx,ef,xi,w,
     &              gf,fnet,kode)
        endif
c ......................................................
c.......dump matrix for examination ......................
        call odumpmat(n3,nz,nnode,kz,ka,gk0,gf,wwwtmp,.false.)
        call ndumpmat(nit,gk,ja,nmax1,ia,nnode,gf,.false.)
c.......................................................
c       if(.true.) then
c.........solve equations with jordan conjugate-gradient itpack ....!
c         write(7,1000) ' time before jcg ',etime(tb),dtime(tb)!
          call dfault(iparm,rparm)
          if(delt.ne.0.0) then
            iparm(1)=1000  ! max number of iteration
          else
            iparm(1)=1000 ! max number of iteration
          endif
          rparm(1)=1d-3   ! zeta, stopping criteria
c         iparm(2)=2      ! level of output (-1:none)
c         iparm(4)=7      ! output unit number
          iparm(5)=1      ! nonsymmetric matrix (1)
          iparm(6)=0      ! non-adaptive (0) (1 doesnt work)
c         iparm(10)=1     ! removes large diagonal entries (doesnt work)
c         do i=1,12
c           print *,i,iparm(i),rparm(i)
c         enddo
          do i=1,nnode
            wrate(i,ii)=0.d0
          enddo
          call jcg(nnode,ia,ja,gk,gf,wrate(1,ii),
     &            iwksp,nw,wksp,iparm,rparm,ier)

c ...................................................
c ....... experimental...............................
c ....... remove rigid body motion ..................
c ....... this could be done better by imposing bc ..
c ....... along the edge of the grid ................
          if(.false.) then
            wratebase=wrate(1,ii)                     !
            do i=1,nnode                              !
              wrate(i,ii)=wrate(i,ii)-wratebase       !
            enddo                                     !
          endif
c ....... experimental, remove rigid body motion ....
c ...................................................

c         do i=1,12
c           print *,i,iparm(i),rparm(i)
c         enddo
c           if(iotogg) then
c             write(list(ipage+1),*) ' relative error=',rparm(1),
c     &              ' in iterations = ',iparm(1)
c             write(list(ipage+2),*) rparm(11),rparm(12)
c             ipage=ipage+2
c             ipage=ipage+1
c           endif
          if(ier.ne.0) then
            call ndumpmat(nit,gk,ja,nmax1,ia,nnode,gf,.true.)
            print *,'jcg error:',ier
            print '(1x,i3,i10,1pg13.6)',(i,iparm(i),rparm(i),i=1,12)
            pause
          endif
c         write(7,1000) ' time after jcg ',etime(tb),dtime(tb) !
c ................................................................!
c       else
c.........solve equations with conjugate-gradient ..................!
c         write(7,1000) ' time before conjug ',etime(tb),dtime(tb)  !
c         call conjug(n3,nz,nnode,1.d-6,gk,ka,gf,www)              !
c         write(7,1000) ' time after conjug ',etime(tb),dtime(tb)   !
c ................................................................!
c       endif
        if(ii.eq.1) then
          do i=1,nnode
            wwwtmp(i)=www(i)+wrate(i,ii)*delt
          enddo
        endif
      enddo
      do i=1,nnode
        www(i)=www(i)+0.5d0*(wrate(i,1)+wrate(i,2))*delt
      enddo
      wmin=1d30
      wmax=-1d30
      wtot=0.d0
      do i=1,nnode,3
        wmin=min(wmin,www(i))
        wmax=max(wmax,www(i))
        wtot=wtot+www(i)
c        write(*,*) (real(www(i+j)),j=0,2)
c        write(*,*) (real(wrate(i+j,1)),j=0,2)
      enddo
      thtot=0.d0
      do i=1,numnp
        thtot=thtot+thick(i)
      enddo
      if(iread.eq.0 .and. ipass.eq.0) then
        ii=1
        do i=1,numnp
          wwworig(i)=www(ii)
          ii=ii+3
        enddo
      endif
      if(iotogg) then
        write(list(ipage+1),*) 
     &       '*********** viscous **********************'
        write(list(ipage+2),1001) time,thtot,
     &              -thtot/wtot/rockice,
     &              wmin,wmax,-1000*(wsave-wmin)/delt
        write(list(ipage+3),*) 
     &       '*******************************************'
        ipage=ipage+3
      endif
      write(92,*) time,wmin
      wsave=wmin
      ipass=1
      end
c=============================================================
      subroutine velemek(xy,n,ekb,nl,xi,w)
      implicit real*8(a-h,o-z)
      dimension xy(2,n),ekb(12,12),eks(12,12)
      dimension dpsix(4),dpsiy(4)
      dimension psi(4),dpsi(4,2),xs(2)
      dimension bb(3,12),bs(2,12),db(3,3),ds(2,2)
      dimension btdb(3,12)
      dimension xi(2,4),w(4)
      common /vmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c.....initialize element arrays
c$doacross local(i,j)
      do i=1,12
        do j=1,12
          ekb(i,j)=0.0d0
          eks(i,j)=0.0d0
        enddo
      enddo
c ... form d-matrices ...............
c ... what is this ctime1? why is it 10 times ctime?
      ctime1=ctime*10
      ctime1=1.0d0
c ..................................................
      ttt3=ttt**3/12.d0
      db(1,1)=ttt3*(2.d0*rmu+rlambda)*ctime1
      db(1,2)=ttt3*rlambda*ctime1
      db(1,3)=0
      db(2,1)=ttt3*rlambda*ctime1
      db(2,2)=ttt3*(2.d0*rmu+rlambda)*ctime1
      db(2,3)=0
      db(3,1)=0
      db(3,2)=0
      db(3,3)=ttt3*rmu*ctime1
      ds(1,1)=ttt*rmu*ctime1
      ds(1,2)=0
      ds(2,1)=0
      ds(2,2)=ttt*rmu*ctime1
c.....begin 2x2 integration loop
      do l=1,nl
        xs(1)=xi(1,l)
        xs(2)=xi(2,l)
        call pgenshape(n,xs,xy,psi,dpsi,detj,dpsix,dpsiy)
        call ploadb(psi,dpsix,dpsiy,bb,bs)
c
c.......accumulate integration point value of integrals
        fac=detj*w(l)
c ..... form kb stiffness matrix ..............
c ..... ok, now here goes the matrix multiplication that 
c ..... generates the bt d b ala book...
c ..... first bbt*db (12x3)*(3x3)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,3
            sum=0.d0
            do k=1,3
              sum=sum+bb(k,i)*db(j,k)
            enddo
            btdb(j,i)=sum
          enddo
        enddo
c then (bbt*db)*bb (12x12)*(3x12)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,12
            sum=0.d0
            do k=1,3
              sum=sum+btdb(k,i)*bb(k,j)
            enddo
            ekb(i,j)=ekb(i,j)+fac*sum
          enddo
        enddo
      enddo
c end of 2x2 integration loop
c
c.....begin 1x1 reduced integration loop
      do l=1,1
        xs(1)=0.d0
        xs(2)=0.d0
        call pgenshape(n,xs,xy,psi,dpsi,detj,dpsix,dpsiy)
        call ploadb(psi,dpsix,dpsiy,bb,bs)
c
c.......accumulate integration point value of integrals
        fac=detj*4.d0
c ..... form ks stiffness matrix ........................
c ..... ok, now here goes the matrix multiplication that 
c ..... generates bst*ds*bs ala book...
c ..... first bst*ds (12x2)*(2x2)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,2
            sum=0.d0
            do k=1,2
              sum=sum+bs(k,i)*ds(j,k)
            enddo
            btdb(j,i)=sum
          enddo
        enddo
c then (bbt*db)*bb (12x12)*(2x12)
c$doacross local(i,j,sum)
        do i=1,12
          do j=1,12
            sum=0.d0
            do k=1,2
              sum=sum+btdb(k,i)*bs(k,j)
            enddo
            eks(i,j)=eks(i,j)+fac*sum
          enddo
        enddo
      enddo
c end of 1x1 integration loop
c
c ... combine kb and ks stiffness matrices .......
      do i=1,12
        do j=1,12
c          print '(2i3,1p3g13.6)',i,j,ekb(i,j),eks(i,j)
          ekb(i,j)=ekb(i,j)+eks(i,j)
        enddo
      enddo
c      pause
c      if(.false.) then
c        do i=1,12
c          print 1003,(ekb(i,j)/ekb(i,i),j=1,12)
c        enddo
c      endif
1003  format(1x,1p6g13.6)
      return
1000  format(1x,'bad jacobian',e10.3)
1001  format(a,/,(2i5,3x,1pe13.6))
1002  format(a,/,(1p2e13.6))
      end
c=============================================================
      subroutine vassmbgk(n3,nz,gk0,ka,kz,ek,n,node)
      implicit real*8(a-h,o-z)
      dimension gk0(n3,nz)
      dimension ka(n3,nz+1),kz(n3)
      dimension ek(12,12),node(n)
c........................
c||||||||||||||||||||||||
c     print *,'in passmb'
      do l=1,n
        i=node(l)
        do m=1,n
          j=node(m)
c.........assemble global stiffness matrix gk
            if(i.eq.j) then
              gk0(i,1)=gk0(i,1)+ek(l,m)
              ka(i,1)=i
            else
              do k=2,kz(i)
                if(ka(i,k).eq.j) then
                  gk0(i,k)=gk0(i,k)+ek(l,m)
                  goto 99
                endif
              enddo
              kz(i)=kz(i)+1
              gk0(i,kz(i))=ek(l,m)
              ka(i,kz(i))=j
            endif
99          continue
        enddo
      enddo
c||||||||||||||||||||||||
c........................
      end
c=============================================================
      subroutine vformgk(nmax,n3,nz,numel,nnode,
     &                  x,y,kx,kkx,ek,xi,w,
     &                  gk0,ka,kz)
      implicit real*8(a-h,o-z)
      dimension x(nmax),y(nmax),kx(nmax,4),kkx(nmax,12)
      dimension gk0(n3,nz)
      dimension ka(n3,nz+1),kz(n3)
      dimension xi(2,4),w(4)
      dimension ek(12,12)
      dimension lm(4),xy(2,4),llm(12)
c      common /vmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c$doacross local(i,j)
      do i=1,nnode
        kz(i)=1
        ka(i,1)=i
        ka(i,nz+1)=0
        do j=1,nz
          ka(i,j)=0
          gk0(i,j)=0.0d0
        enddo
      enddo
      do nel=1,numel
        do l=1,4
          lm(l)=kx(nel,l)
          xy(1,l)=x(lm(l))
          xy(2,l)=y(lm(l))
        enddo
c$doacross local(l)
        do l=1,12
          llm(l)=kkx(nel,l)
        enddo
        call velemek(xy,4,ek,4,xi,w)
c..................................
        if(.false.) then
          write(7,*) nel
          do i=1,12
            write(7,*) (ek(i,j),j=1,12)
          enddo
c          pause
        endif
c..................................
        call vassmbgk(n3,nz,gk0,ka,kz,ek,12,llm)
c..................................
c        do i=1,nnode
c          print 1000,(gk0(i,j),j=1,nnode)
c        enddo
c        pause
c..................................
      enddo
c$doacross local(i)
      do i=1,nnode
c        print *,kz(i),ka(i,nz+1)
        ka(i,nz+1)=kz(i)
      enddo
c|||||||||||||||||||||||||||||||||||
c...................................
1000  format(1x,10f8.3)
      end
c=============================================================
      subroutine vformgf(nmax,n3,numel,nnode,
     &                  x,y,thick,www,kx,kkx,ef,xi,w,
     &                  gf,fnet,kode)
      implicit real*8(a-h,o-z)
      dimension x(nmax),y(nmax),kx(nmax,4),kkx(nmax,12)
      dimension gf(n3),www(n3)
      dimension thick(nmax),kode(nmax)
      dimension xi(2,4),w(4)
      dimension ef(12)
      dimension lm(4),xy(2,4),llm(12)
      common /vmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c$doacross local(i)
      do i=1,nnode
        gf(i)=0.d0
      enddo
      fnet=0.d0
      do nel=1,numel
        ethick=0.d0
        ewww=0.d0
        do l=1,4
          lm(l)=kx(nel,l)
          xy(1,l)=x(lm(l))
          xy(2,l)=y(lm(l))
          ethick=ethick+thick(lm(l))
          ewww=ewww+www((lm(l)-1)*3+1)
        enddo
c$doacross local(l)
        do l=1,12
          llm(l)=kkx(nel,l)
        enddo
        ethick=ethick/dble(4)
        ewww=ewww/dble(4)
c        print *,nel,real(ethick),real(ewww),
c     &          real(rhog*ethick+rockice*rhog*ewww)
        ethick=rhog*ethick
        ewww=rockice*rhog*ewww
c        print *,nel,real(ethick),real(ewww),real(ethick+ewww)
        ethick=ethick+ewww
        fnet=fnet+ethick
        call pelemgf(ethick,xy,4,ef,4,xi,w)
c..................................
        if(.false.) then
          write(7,*) nel
          do i=1,12
            write(7,*) ef(i)
          enddo
c          pause
        endif
c..................................
        call passmbgf(n3,gf,ef,12,llm)
c        call aplybc(n3,gf,nmax,kode,nnode)
c..................................
c        do i=1,nnode
c          print 1000,gf(i)
c        enddo
c        pause
c..................................
      enddo
c|||||||||||||||||||||||||||||||||||
c...................................
1000  format(1x,10f8.3)
      end
c=============================================================
      subroutine vloadgk(n3,nz,nnode,gk0,gk,ka,kz,
     &                    nit,nmax1,ia,ja,iwork,nmax,kode)
      implicit real*8(a-h,o-z)
c ....arrays for itpack sparse storage...................
      dimension gk(nit),ja(nit),ia(nmax1)
      dimension iwork(nit)
      dimension gk0(n3,nz)
      dimension ka(n3,nz+1),kz(n3),kode(nmax)
      data big /1d30/
      call sbini(nnode,nit,ia,ja,gk,iwork)
      do i = 1,nnode
c        n=1+(i-1)/3
c        if(kode(n).eq.1) then
c          call sbsij(nnode,nit,ia,ja,gk,iwork,i,i,big,0,1,6,ier)
c        else          
          do j = 1,kz(i)
            jg=ka(i,j)
            gkij=gk0(i,j)
            call sbsij(nnode,nit,ia,ja,gk,iwork,i,jg,gkij,0,1,6,ier)
          enddo 
c        endif
      enddo
      call sbend(nnode,nit,ia,ja,gk,iwork)
      end

c=======================================
      subroutine psetmat
c ... elastic material properties ...
      implicit real*8(a-h,o-z)                                          
      common /pmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c ... material properties and thickness of crust
c ... youngs modulus (nt/m**2)
      e=5.d10
      e=1.d11
c      e=1d20
c ... poisons ratio (dimensionless)
      rnu=0.5d0
c ... lames coefficients
      rlambda=rnu*e/(1.d0-rnu**2)
      rmu=0.5d0*e/(1.d0+rnu)
c ... thickness of crust (90-130 km) (110,000 m)
      ttt=110.d0*1000.d0
c      ttt=30*1000
c ... rhoice*grav (nt/m**3) times thickness, yields pressure, nt/m**2
      grav=9.8d0*0.3816d0
      rhor=5000.d0
      rhor=4000.d0
      rhow=1092.d0
      rhoi=917.d0
      rhog=rhoi*grav
      rockice=2.0d0*rhor/rhoi
c      rockice=1.0d0*rhor/rhoi
c ... relaxation time constant ......!
      relax=6000.d0                  !
      relax=500.d0                  !
      ctime=3.16516d8*relax/6000.d0  !
      end
c=======================================
      subroutine vsetmat
c ... viscous material properties ...
      implicit real*8(a-h,o-z)                                          
      common /vmatprop/ rmu,rlambda,ttt,rhog,ctime,rockice
c ... material properties and thickness of crust
c ... youngs modulus (nt/m**2)
      e=5.d10
      e=1.d11
      e=1.d17
c ... poisons ratio (dimensionless)
      rnu=0.5d0
c ... lames coefficients
      rlambda=rnu*e/(1.d0-rnu**2)
      rmu=0.5d0*e/(1.d0+rnu)
c ... thickness of crust (90-130 km) (110,000 m)
      ttt=110.d0*1000.d0
       ttt=40*1000
c ... rhoice*grav (nt/m**3) times thickness, yields pressure, nt/m**2
      grav=9.8d0*0.3816d0
      rhor=5000.d0
      rhor=4000.d0
      rhow=1092.d0
      rhoi=917.d0
      rhog=rhoi*grav
      rockice=2.0d0*rhor/rhoi
c      rockice=1.0d0*rhor/rhoi
c ... relaxation time constant ......! NOT USED ...
      relax=6000.d0                  !
      ctime=3.16516d8*relax/6000.d0  !
      end
c=============================================================
      subroutine veplate(itime,numnp,numel,x,y,kx,thick,kode,
     &                 delt,www,wrate,time,wwworig,
     &                 wmin,wmine,wminv,fnete,fnetv)
c-----------------------------------------------------------------------
c 4th order visco-elastic plate solver with one-time generation of stiffness
c matrix. uses my sparse matrix storage for the static matrice
c and itpack sprse storage for the time-dependent modified matrices, and c             ipage=ipage+2

c iterative matrix solver. j fastook july, 1999
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)                                          
      include "parameter.h"
      parameter(nmax=maxnum,nz=27, nz1=nz+1, n3=3*nmax)
      parameter(nit=n3*nz,nmax1=n3+1,nmax3=n3*3)
      parameter(itmax=n3/10,ncg=4*itmax,nw=4*n3+ncg)
c ....arrays for itpack sparse storage...................
c      dimension gk(nit),ja(nit),ia(nmax1)
c      dimension iwksp(nmax3),wksp(nw),iwork(nit)
c      dimension iparm(12),rparm(12)
      dimension thick(nmax),x(nmax),y(nmax),kx(nmax,4)
      dimension www(n3),wrate(n3,2),wwworig(nmax)
      common /elastic/ wwwe(n3),werate(n3,2),wwweorig(nmax)
      common /viscous/ wwwv(n3),wvrate(n3,2),wwwvorig(nmax)
      dimension www0(n3)
      dimension kode(nmax)
c      dimension kkx(nmax,12),ltemp(nmax,3)
c      dimension xi(2,4),w(4)
c      dimension ek(12,12),ef(12)
c      dimension gf(n3)
c      dimension gk0(n3,nz)
c      dimension ka(n3,nz1),kz(n3)
c      dimension wwwtmp(n3)
      character*80 hed
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      common /line/ np,nline(1000)
      parameter(npage=39)
      character*80 list(1000)
      common /iolist/ list
      data www0 /n3*0.d0/
      logical file40
c ... timer stuff, for sgi only ...!
      real*4 tb(2)
c     real*4 tb(2),etime,dtime     !
c     external etime,dtime         !
c .................................!
      save ipass
c      save kkx,wsave,nnode,xi,w,gk0,ka,kz
c      save wwwsave,iparm,rparm,iread
c      save wwwe,werate,wwweorig
c      save wwwv,wvrate,wwwvorig
      data ipass /0/,iread /0/
      nnode=3*numnp
      if(ipass.eq.0) then
        ipass=1
        print *,'reading bedrock depression file'
        inquire(file='fort.40',exist=file40)
        iok=-1
        if(file40) then
          rewind 40
          read(40,*,iostat=iok) nnn
        endif
        if(iok.eq.0) then
          if(nnn.ne.nnode) then
            print *,'problems:'
            print *,'incompatible with current nnode=',nnode,nnn
            iok=1
          endif
          do i=1,nnode
            read(40,*) ii,wwwe(i),wwwv(i)
            werate(i,2)=wwwe(i)
            wvrate(i,2)=wwwv(i)
            if(i.ne.ii) then
              print *,'problems:reading wwwe,wwwv'
              print *,'incompatible with current (www)'
              iok=1
            endif
          enddo
          do i=1,nnode
            read(40,*) ii,werate(I,1),wvrate(I,1)
            if(i.ne.ii) then
              print *,'problems:reading wrates(1)'
              print *,'incompatible with current (www)'
              iok=1
            endif
          enddo
          do i=1,nnode
            read(40,*) ii,werate(I,2),wvrate(I,2)
            if(i.ne.ii) then
              print *,'problems:reading wrates(2)'
              print *,'incompatible with current (www)'
              iok=1
            endif
          enddo
          do i=1,numnp
            read(40,*) ii,wwweorig(i),wwwvorig(i)
            if(i.ne.ii) then
              print *,'problems:reading wwworigs'
              print *,'incompatible with current (wwworig)'
              iok=1
            endif
          enddo
          do i=1,numnp*3
            www(i)=wwwe(i)+wwwv(i)
c           if(mod(i,30).eq.1) print *,i,wwwe(i),wwwv(i)
            wrate(i,1)=werate(i,1)+wvrate(i,1)
            wrate(i,2)=werate(i,2)+wvrate(i,2)
          enddo
          do i=1,numnp
            wwworig(i)=wwweorig(i)+wwwvorig(i)
          enddo
          print *,' bedrock depression file found'
          print *,'    and read successfully '
          if(itime.lt.0) then
            print *,'abandoning unloading'
            rewind 40
            return
          endif
          iread=1
        endif
        if(iok.ne.0) then
          print *,' none found, set to zero ... (or part or www '
c$doacross local(i)
          if(.false.) then
            do i=1,nnode
              www(i)=0.d0
              wwwe(i)=0.d0
              wwwv(i)=0.d0
              wrate(i,2)=0.d0
              werate(i,2)=0.d0
              wvrate(i,2)=0.d0
            enddo
          elseif(.false.) then
            do i=1,nnode
              wwwe(i)=1.d0*www(i)
              wwwv(i)=0.d0*www(i)
              werate(i,2)=1.d0*wrate(i,2)
              wvrate(i,2)=0.d0*wrate(i,2)
            enddo
          else
            do i=1,numnp*3
              www(i)=wwwe(i)+wwwv(i)
c             if(mod(i,30).eq.1) print *,i,wwwe(i),wwwv(i)
              wrate(i,1)=werate(i,1)+wvrate(i,1)
              wrate(i,2)=werate(i,2)+wvrate(i,2)
            enddo
          endif
        endif
        write(hed,*) ' time=',nint(time-delt)
c        write(88) hed
c        write(88) (www(i),i=1,nnode,3)
c        write(88) (thick(i),i=1,numnp)
c        write(88) (1000*wrate(i,1),i=1,nnode,3)
      endif
      if(.true.) then
c ..... elastic solution .............................................
        call eplate(itime,numnp,numel,x,y,kx,thick,kode,
     &                 delt,wwwe,werate,wmine,time,wwweorig,fnete,wwwe)
      endif
      if(.true.) then
c ..... viscous solution .............................................
        call vplate(itime,numnp,numel,x,y,kx,thick,kode,
     &                 delt,wwwv,wvrate,wminv,time,wwwvorig,fnetv,wwwv)
      endif
      write(list(ipage+1),*) ' fnete,fnetv=',fnete*1d-6,fnetv*1d-6
      ipage=ipage+1
      do i=1,numnp*3
        www(i)=wwwe(i)+wwwv(i)
c        if(mod(i,30).eq.1) print *,i,wwwe(i),wwwv(i)
        wrate(i,1)=werate(i,1)+wvrate(i,1)
        wrate(i,2)=werate(i,2)+wvrate(i,2)
      enddo
c      do i=1,numnp
c        wwworig(i)=wwweorig(i)+wwwvorig(i)
c      enddo
c      pause
      wmin=1d30
      wmax=-1d30
      wmine=1d30
      wmaxe=-1d30
      wminv=1d30
      wmaxv=-1d30
      do i=1,numnp*3,3
        wmin=min(wmin,www(i))
        wmax=max(wmax,www(i))
        wmine=min(wmine,wwwe(i))
        wmaxe=max(wmaxe,wwwe(i))
        wminv=min(wminv,wwwv(i))
        wmaxv=max(wmaxv,wwwv(i))
      enddo
      if(.false.) then
        print *,'total  :',wmin,wmax
        print *,'elastic:',wmine,wmaxe
        print *,'viscous:',wminv,wmaxv
        do i=1,np
          ii=nline(i)
          print *,ii,real(wwwe((ii-1)*3+1)),
     &               real(wwwv((ii-1)*3+1)),
     &               real(www((ii-1)*3+1)),real(thick(ii))
        enddo
      endif
      end
C---------------------------------------------
      SUBROUTINE WRITEDEPN(NUMNP,NNODE)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      parameter(nmax=maxnum,n3=3*nmax)
      CHARACTER*1 CHAR
      common /elastic/ wwwe(n3),werate(n3,2),wwweorig(nmax)
      common /viscous/ wwwv(n3),wvrate(n3,2),wwwvorig(nmax)
      PRINT *,'IN WRITEDEPN',NUMNP,NNODE
      PRINT *,'   TO WRITE OUT BACKUP OF BEDROCK DEPRESSION '
      PRINT *,'   INPUT Y'
      READ(*,1002) CHAR
      WRITE(99,1002) CHAR
1002  FORMAT(A1)
      IF(CHAR.EQ.'Y' .OR. CHAR.EQ.'y') THEN
        REWIND 45
        WRITE(45,*) NNODE
        DO I=1,NNODE
          WRITE(45,*) I,wwwe(I),wwwv(I)
        ENDDO
        DO I=1,NNODE
          WRITE(45,*) I,werate(I,1),wvrate(I,1)
        ENDDO
        DO I=1,NNODE
          WRITE(45,*) I,werate(I,2),wvrate(I,2)
        ENDDO
        DO I=1,NUMNP
          WRITE(45,*) I,wwweorig(I),wwwvorig(i)
        ENDDO
      ENDIF
      END
