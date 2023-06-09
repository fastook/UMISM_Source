      dimension lev(4,10)
c original
c     data lev    /  1, 240,   0, 255,
c    &              40,   0,   0, 255,
c    &              90,   0, 255, 255,
c    &             130,   0, 255,   0,
c    &             140,   0, 255,   0,
c    &             180, 255, 255,   0,
c    &             190, 255, 255,   0,
c    &             220, 255,   0,   0,
c    &             230, 255,   0,   0,
c    &             250, 200,   0,   0/
c simple
      data lev    /  1, 240,   0, 255,
     &               2,   0,   0, 255,
     &               3,   0, 255, 255,
     &               4,   0, 255,   0,
     &               5, 255, 255,   0,
     &               6, 255, 127,   0,
     &               7, 255,   0,   0,
     &               8, 255,   0,   0,
     &             230, 255,   0,   0,
     &             250, 200,   0,   0/
1     print *,'input level spacing and min, max values'
      read(*,*,end=999) rlevsp,zmin,zmax
      NLEV=NINT((zmax-zmin)/RLEVSP)+1
c      print *,nlev
      COLMIN=250
      COLMAX=50
      COLSTEP=(COLMAX-COLMIN)/(NLEV-1)
      ZSTEP=(ZMAX-ZMIN)/(NLEV-1)
      rewind 13
      DO I=1,NLEV
        ZVAL=ZMIN+(I-1)*ZSTEP
        ICOL=COLMIN+(I-1)*COLSTEP
c        print *,zval,icol,zval+zstep
        WRITE(13,1130) ZVAL,ICOL,ICOL,ICOL,ZVAL+ZSTEP,ICOL,ICOL,ICOL
      ENDDO
1130  FORMAT(1X,F10.3,3I4,F10.3,3I4)
      rewind 13
      slope=(zmax-zmin)/7
      do i=1,7
        first=zmin+slope*(lev(1,i)-1)
        rlast=zmin+slope*(lev(1,i+1)-1)
c        if(rlast.ge.10.) then
c          print 100,first,(lev(j,i),j=2,4),rlast,(lev(j,i+1),j=2,4)
c          write(13,100) first,(lev(j,i),j=2,4),rlast,(lev(j,i+1),j=2,4)
c        else
          print 101,first,(lev(j,i),j=2,4),rlast,(lev(j,i+1),j=2,4)
          write(13,101) first,(lev(j,i),j=2,4),rlast,(lev(j,i+1),j=2,4)
c        endif
      enddo
100   format(1x,f10.0,3i5,f10.0,3i5)
101   format(1x,f10.2,3i5,f10.2,3i5)
c      goto 1
999   end
