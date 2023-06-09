      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      parameter(nmax=MAXNUM)
      real*4 rplot(nmax),xplot(nmax),eplot(nmax),tplot(nmax)
      real*4 aplot(nmax),dplot(nmax)
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      COMMON /LOCAL/ SIGMA,TMAA,TMS,PI2
      acom=-9.62376690d0
      acom=-6
      tnsl=40.
      xpole=0.0
      ypole=0.0
      call setrig
      n=0
      emin=0.
      do rl=85.,50.,-1.d-1
        call POLREC(rl,0.d0,X,Y)
        do elev=emin,10000.d0,1.d-1
          acc=ACCUM(x/1000.,y/1000.,elev,0.d0,0.d0,tnsl,TS) 
          if(acc.gt.0) then
            n=n+1
            xplot(n)=x
            rplot(n)=rl
            eplot(n)=elev
            emin=elev
            tplot(n)=ts
            aplot(n)=abl
            goto 99
          endif
        enddo
99      continue
      enddo
      call grstrt(800,800)
      call linclr(1)
      xmin=1e30
      xmax=-1e30
      ymin=1e30
      ymax=-1e30
      do i=1,n
        xmin=min(xmin,rplot(i))
        xmax=max(xmax,rplot(i))
        ymin=min(ymin,eplot(i))
        ymax=max(ymax,eplot(i))
      enddo
      print *,'elevation',real(xmin),real(xmax),real(ymin),real(ymax)
      call window(real(xmax),real(xmin),real(ymin),real(ymax))
      call move(real(rplot(1)),real(eplot(1)))
      do i=2,n
        call draw(real(rplot(i)),real(eplot(i)))
        dplot(i)=1e5*(eplot(i)-eplot(i-1))/(xplot(i)-xplot(i-1))
c        print *,rplot(i),dplot(i)
      enddo
      call linclr(2)
      xmin=1e30
      xmax=-1e30
      ymin=1e30
      ymax=-1e30
      do i=1,n
        xmin=min(xmin,rplot(i))
        xmax=max(xmax,rplot(i))
        ymin=min(ymin,tplot(i))
        ymax=max(ymax,tplot(i))
      enddo
      print *,'ts       ',real(xmin),real(xmax),real(ymin),real(ymax)
      call window(real(xmax),real(xmin),real(ymin),real(ymax))
      call move(real(rplot(1)),real(tplot(1)))
      do i=2,n
        call draw(real(rplot(i)),real(tplot(i)))
      enddo
      call linclr(3)
      xmin=1e30
      xmax=-1e30
      ymin=1e30
      ymax=-1e30
      do i=2,n
        xmin=min(xmin,rplot(i))
        xmax=max(xmax,rplot(i))
        ymin=min(ymin,dplot(i))
        ymax=max(ymax,dplot(i))
      enddo
      print *,'gradient ',real(xmin),real(xmax),real(ymin),real(ymax)
      call window(real(xmax),real(xmin),real(ymin),real(ymax))
      call move(real(rplot(2)),real(dplot(2)))
      do i=3,n
        call draw(real(rplot(i)),real(dplot(i)))
      enddo
      call grstop
      do i=1,n
        write(*,100) rplot(i),xplot(i)/1000.,eplot(i),dplot(i)
      enddo
100   format(1x,5g13.6)
      end

