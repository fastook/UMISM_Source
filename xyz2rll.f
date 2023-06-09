      implicit real*8(a-h,o-z)
      open(1,file='gmt.XYZ')
      open(11,file='gmt.rll')
      rlat0=42.5;rlon0=67.5
      rlat =-45;rlon =270
      call polrec(0,rlat0,rlon0,0.0d0,0.0d0,
     &              rlat,rlon,xdble,ydble)
      print *,xdble,ydble
      call polrec(1,rlat0,rlon0,0.0d0,0.0d0,
     &              rlat1,rlon1,xdble,ydble)
      print *,rlat,rlon
      print *,rlat1,rlon1
      do i=1,10000000
        read(1,*,end=99) x,y,z
        x=x*1000d0;y=y*1000d0
        call polrec(1,rlat0,rlon0,0.0d0,0.0d0,
     &              rlat1,rlon1,x,y)
        write(11,*) rlon1,rlat1,z
      enddo
      print *,'increase nmax',nmax
      stop
99    continue
      end
      include "proj-st.f"
