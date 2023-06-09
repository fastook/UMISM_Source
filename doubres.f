      include "parameter.h"
      parameter(nmax=MAXNUM)
      implicit real*8(a-h,o-z)
      logical iodd
      character hed*80
      dimension kode(nmax),x(nmax),y(nmax),htice(nmax),adot(nmax),
     &          fract(nmax),psurf(nmax),bdrock(nmax),flowa(nmax),
     &          sldgb(nmax),kx(nmax,4),const(nmax),afudge(nmax),
     &          ibflux(nmax,2),bflux(nmax),temp(nmax),itype(nmax)
      dimension kode2(nmax),x2(nmax),y2(nmax),htice2(nmax),adot2(nmax),
     &          fract2(nmax),psurf2(nmax),bdrock2(nmax),flowa2(nmax),
     &          sldgb2(nmax),const2(nmax),afudge2(nmax),
     &          ibflux2(nmax,2),bflux2(nmax),temp2(nmax),itype2(nmax)
      external iodd
c ... read input file for resolution doubling
      read(1,1000) hed,numnp,numel,numcol,numlev,numgbc,ndt,inter,dt
      print *,hed
      if(numnp.gt.nmax) then
        print *,'numnp=',numnp,' nmax=',nmax,' increase nmax'
        stop
      endif
      print *,numnp,' nodes,',numel,' elements'
      do n=1,numnp
        read(1,1001) num,kode(n),x(n),y(n),htice(n),
     &                   adot(n),fract(n),psurf(n),
     &                   bdrock(n),flowa(n),sldgb(n),
     &                   temp(n),itype(n),afudge(n)
      enddo
      do n=1,numel
        read(1,1002) num,(kx(num,i),i=1,4),const(n)
      enddo
      if(numgbc.gt.0) then
        do n=1,numgbc
          read(1,1007) ibflux(n,1),ibflux(n,2),bflux(n)
        enddo
      endif
c
      numcol2=numcol*2-1
      numlev2=numlev*2-1
      numnp2=numcol2*numlev2
      if(numnp2.gt.nmax) then
        print *,'new node number',numnp2,' >',nmax
        stop
      endif
      numel2=(numcol2-1)*(numlev2-1)
      numgbc2=numgbc*2
      do i=1,numnp
        iy=nint(real(i+1)/real(numcol))
        ix=mod(i-1,numcol)+1
        call k2ij(i,numcol,ix,iy)
        j=ix+(iy-1)*numcol
        call ij2k(ix,iy,numcol,j)
      enddo
      do iy=1,numlev2
        if(iodd(iy)) then
          iiy=iy/2+1
          do jx=1,numcol2
            if(iodd(jx)) then
              jjx=jx/2+1
              call ij2k(jx,iy,numcol2,k2)
              call ij2k(jjx,iiy,numcol,k)
              x2(k2)=x(k)
              y2(k2)=y(k)
              htice2(k2)=htice(k)
              adot2(k2)=adot(k)
              fract2(k2)=fract(k)
              psurf2(k2)=psurf(k)
              bdrock2(k2)=bdrock(k)
              flowa2(k2)=flowa(k)
              sldgb2(k2)=sldgb(k)
              afudge2(k2)=afudge(k)
              temp2(k2)=temp(k)
              kode2(k2)=kode(k)
              itype2(k2)=itype(k)
            else
              jjx=jx/2+1
              call ij2k(jx,iy,numcol2,k2)
              call ij2k(jjx,iiy,numcol,k)
              call ij2k(jjx-1,iiy,numcol,l)
              x2(k2)=(x(k)+x(l))/2.
              y2(k2)=(y(k)+y(l))/2.
c              y2(k2)=y(k)
              htice2(k2)=(htice(k)+htice(l))/2.
              adot2(k2)=(adot(k)+adot(l))/2.
              fract2(k2)=(fract(k)+fract(l))/2.
              psurf2(k2)=(psurf(k)+psurf(l))/2.
              bdrock2(k2)=(bdrock(k)+bdrock(l))/2.
              flowa2(k2)=(flowa(k)+flowa(l))/2.
              sldgb2(k2)=(sldgb(k)+sldgb(l))/2.
              afudge2(k2)=(afudge(k)+afudge(l))/2.
              temp2(k2)=(temp(k)+temp(l))/2.
              kode2(k2)=min(kode(k),kode(l))
              itype2(k2)=max(itype(k),itype(l))
            endif
          enddo
        else
          iiy=iy/2+1
          do jx=1,numcol2
            if(iodd(jx)) then
              jjx=jx/2+1
              call ij2k(jx,iy,numcol2,k2)
              call ij2k(jjx,iiy,numcol,k)
              call ij2k(jjx,iiy-1,numcol,l)
c              x2(k2)=x(k)
              x2(k2)=(x(k)+x(l))/2.
              y2(k2)=(y(k)+y(l))/2.
              htice2(k2)=(htice(k)+htice(l))/2.
              adot2(k2)=(adot(k)+adot(l))/2.
              fract2(k2)=(fract(k)+fract(l))/2.
              psurf2(k2)=(psurf(k)+psurf(l))/2.
              bdrock2(k2)=(bdrock(k)+bdrock(l))/2.
              flowa2(k2)=(flowa(k)+flowa(l))/2.
              sldgb2(k2)=(sldgb(k)+sldgb(l))/2.
              afudge2(k2)=(afudge(k)+afudge(l))/2.
              temp2(k2)=(temp(k)+temp(l))/2.
              kode2(k2)=min(kode(k),kode(l))
              itype2(k2)=max(itype(k),itype(l))
            else
              jjx=jx/2+1
              call ij2k(jx,iy,numcol2,k2)
              call ij2k(jjx,iiy,numcol,k)
              call ij2k(jjx,iiy-1,numcol,l)
              call ij2k(jjx-1,iiy,numcol,m)
              call ij2k(jjx-1,iiy-1,numcol,n)
              x2(k2)=(x(k)+x(l)+x(m)+x(n))/4.
              y2(k2)=(y(k)+y(l)+y(m)+y(n))/4.
              htice2(k2)=(htice(k)+htice(l)+htice(m)+htice(n))/4.
              adot2(k2)=(adot(k)+adot(l)+adot(m)+adot(n))/4.
              fract2(k2)=(fract(k)+fract(l)+fract(m)+fract(n))/4.
              psurf2(k2)=(psurf(k)+psurf(l)+psurf(m)+psurf(n))/4.
              bdrock2(k2)=(bdrock(k)+bdrock(l)+bdrock(m)+bdrock(n))/4.
              flowa2(k2)=(flowa(k)+flowa(l)+flowa(m)+flowa(n))/4.
              sldgb2(k2)=(sldgb(k)+sldgb(l)+sldgb(m)+sldgb(n))/4.
              afudge2(k2)=(afudge(k)+afudge(l)+afudge(m)+afudge(n))/4.
              temp2(k2)=(temp(k)+temp(l)+temp(m)+temp(n))/4.
              kode2(k2)=min(kode(k),kode(l),kode(m),kode(n))
              itype2(k2)=max(itype(k),itype(l),itype(m),itype(n))
            endif
          enddo
        endif
      enddo      
      do n=1,numel
        call k2ij(n,numcol-1,i,j)
        i1=2*i-1
        j1=2*j-1
        i2=i1+1
        j2=j1
        i3=i1
        j3=j1+1
        i4=i1+1
        j4=j1+1
        call ij2k(i1,j1,numcol2-1,k1)
        call ij2k(i2,j2,numcol2-1,k2)
        call ij2k(i3,j3,numcol2-1,k3)
        call ij2k(i4,j4,numcol2-1,k4)
        const2(k1)=const(n)
        const2(k2)=const(n)
        const2(k3)=const(n)
        const2(k4)=const(n)
      enddo
c ... write output file for resolution doubling
      write(11,1000) hed,numnp2,numel2,numcol2,numlev2,numgbc2,ndt,inter,dt
      print *,'new',numnp2,' nodes,',numel2,' elements'
      do n=1,numnp2
        write(11,1001) n,kode2(n),x2(n),y2(n),htice2(n),
     &                   adot2(n),fract2(n),psurf2(n),
     &                   bdrock2(n),flowa2(n),sldgb2(n),
     &                   temp2(n),itype2(n),afudge2(n)
      enddo
      n=1
      do i=1,numlev2-1
        istart=1+(i-1)*numcol2
        do j=1,numcol2-1
          kxx1=istart+j-1
          kxx2=kxx1+1
          kxx3=kxx2+numcol2
          kxx4=kxx3-1
          write(11,1002) n,kxx1,kxx2,kxx3,kxx4,const2(n)
          n=n+1
        enddo
      enddo
      if(numgbc2.gt.0) then
        m=1
        do n=1,numgbc
          call k2ij(ibflux(n,1),numcol,i1,j1)
          ii1=i1*2-1
          jj1=j1*2-1
          call k2ij(ibflux(n,2),numcol,i3,j3)
          ii3=i3*2-1
          jj3=j3*2-1
          ii2=(ii1+ii3)/2
          jj2=(jj1+jj3)/2
          call ij2k(ii1,jj1,numcol2,k1)
          call ij2k(ii2,jj2,numcol2,k2)
          call ij2k(ii3,jj3,numcol2,k3)
          ibflux2(m,1)=k1
          ibflux2(m,2)=k2
          bflux2(m)=bflux(n)
          ibflux2(m+1,1)=k2
          ibflux2(m+1,2)=k3
          bflux2(m+1)=bflux(n)
          m=m+2
        enddo
        do n=1,numgbc2
          write(11,1007) ibflux2(n,1),ibflux2(n,2),bflux2(n)
        enddo
      endif
1000  format (a80,/,7I6,f8.1)
1001  FORMAT(I6,I4,1P,2E12.5,0P,F10.2,F7.2,F9.4,F10.3,F10.1,
     &          F10.5,2F10.5,I5,F10.3)
1002  format(5i6,1pe17.10)
1007  format(2i6,1pe13.6)

      end
c=======================================
      subroutine ij2k(i,j,numcol,k)
      k=i+(j-1)*numcol
      end
c=======================================
      subroutine k2ij(k,numcol,i,j)
      i=mod(k-1,numcol)+1
      j=1+(k-1)/numcol
      end
c=======================================
      logical function iodd(i)
      if(2*(i/2).eq.i) then
        iodd=.false.
      else
        iodd=.true.
      endif
      end
