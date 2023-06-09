      subroutine threed(ipass,npts,nelem,x,y,htice,bdrock,kx)
      include "parameter.h"
      parameter(nmax=MAXNUM)
      !include "fshort.h"
      real*8 x,y,htice,bdrock
      dimension x(nmax),y(nmax),htice(nmax),bdrock(nmax),kx(nmax,4)
      real xyzbed(3,nmax),xyzsrf(3,nmax),xyzgrd(3,nmax)
      integer icxyzb(nmax),icxyzs(nmax),kxp(nmax,4)
      common /plotdata/ numnp,numel,xyzbed,xyzsrf,xyzgrd,
     &                  kxp,icxyzb,icxyzs
      return
      end

