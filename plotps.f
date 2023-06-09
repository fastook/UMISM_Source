      PARAMETER(NPLMAX=50,NMAX=1000)
      DIMENSION X(NMAX,NPLMAX),Y(NMAX,NPLMAX),NPTS(NPLMAX)
      dimension ictype(nplmax),istype(nplmax),istypv(nplmax)
      CHARACTER*20 TITLE(NPLMAX)
      character ptitle*80,xtitle*20,ytitle*20
      character*1 comma,equal,quote
      data comma /','/,equal /'='/, quote /''''/
      DO NP=1,NPLMAX
        DO N=1,NMAX
          read(1,202,end=999) xp,yp,ip
          IF(XP.EQ.-99999.) THEN
            NPLOT=NP
            read(1,100) title(nplot)
            GOTO 10
          ENDIF
          X(N,NP)=XP
          Y(N,NP)=YP
          NPTS(NP)=N
        ENDDO
        PRINT *,'TOO MANY POINTS IN PLOT',NP,' NMAX=',NMAX
        STOP
10      CONTINUE
      ENDDO
      PRINT *,'TOO MANY PLOTS ',NPLMAX
      STOP
999   CONTINUE
      xmin=1e30
      ymin=xmin
      xmax=-xmin
      ymax=xmax
      DO NP=1,NPLOT
        DO N=1,NPTS(NP)
          xmax=max(xmax,x(n,np))
          xmin=min(xmin,x(n,np))
          ymax=max(ymax,y(n,np))
          ymin=min(ymin,y(n,np))
        ENDDO
      ENDDO
      write(*,*) 'xmin,xmax=',xmin,xmax
      write(*,*) 'ymin,ymax=',ymin,ymax
      write(*,*) '-1 -> autoscale, equal axis'
      write(*,*) '+1 -> autoscale with plotc'
      write(*,*) ' 0 -> manual scaling'
      write(*,*) ' 2 -> use values above'
      read(*,*) iscale
      write(*,*) 'input plot title'
      read(*,200) ptitle
      write(*,*) 'input x-axis title'
      read(*,100) xtitle
      write(*,*) 'input y-axis title'
      read(*,100) ytitle
      if(iscale.eq.0) then
        write(*,*) 'input xmax'
        read(*,*) xmax
        write(*,*) 'input xmin'
        read(*,*) xmin
        write(*,*) 'input ymax'
        read(*,*) ymax
        write(*,*) 'input ymin'
        read(*,*) ymin
      endif
      do i=1,nplot
        if(i.lt.11) then
          ictype(i)=1
        elseif(i.lt.21) then
          ictype(i)=2
        elseif(i.lt.21) then
          ictype(i)=3
        elseif(i.lt.21) then
          ictype(i)=4
        else
          ictype(i)=5
        endif
        istype(i)=mod(i,10)+1
        istypv(i)=mod(i,10)+1
      enddo
      write(11,*) '$ptype'
      write(11,*) '  irot=',1,comma
      write(11,*) '  ncurv=',-nplot,comma
      write(11,*) '  nsymb=',0,comma
      write(11,*) '  itype=2',comma
      write(11,*) '  ictype=',(ictype(i),comma,i=1,nplot)
      write(11,*) '  istype=',((i),comma,i=1,nplot)
      write(11,*) '  istypv=',(istypv(i),comma,i=1,nplot)
      write(11,*) '  isize=5',comma
      write(11,*) '  ifreq=1',comma
      write(11,*) '  ititle= 1,1,1',comma
      write(11,*) '  ilegnd=-1',comma
      write(11,*) '  xleg=8.1',comma
      write(11,*) '  yleg=6.75',comma
      write(11,*) '  aleg(1,2)=',quote,title(1),quote,comma
      do i=2,nplot
        write(11,*) '            ',quote,title(i),quote,comma
      enddo
      write(11,*) '  ilegbx=1',comma
      write(11,*) '  modec=2',comma
      write(11,*) '  iclip=1',comma
      write(11,*) '  idebug=1,1,1,1,1,1,1,1,1,1,1',comma
      write(11,*) '$end'
      write(11,*) '$axes'
      if(iscale.eq.0 .or. iscale.eq.2) then
        write(11,*) '  iscale=1',comma
        write(11,*) '  xmin=',xmin,comma
        write(11,*) '  xmax=',xmax,comma
        write(11,*) '  ymin=',ymin,comma
        write(11,*) '  ymax=',ymax,comma
      elseif(iscale.eq.1) then
        write(11,*) '  iscale=0',comma
      elseif(iscale.eq.-1) then
        write(11,*) '  iscale=4',comma
      endif
      write(11,*) '  ixaxis=2',comma
      write(11,*) '  iyaxis=2',comma
      write(11,*) '$end'
      write(11,*) xtitle
      write(11,*) ytitle
      write(11,*) ptitle
      do np=1,nplot
        write(11,*) npts(np)
        do i=1,npts(np)
          write(11,*) x(i,np),y(i,np)
       enddo
      enddo
100   format(a20)
200   format(a80)
202   format(10x,g13.6,2x,g13.6,i13)
      END


              
