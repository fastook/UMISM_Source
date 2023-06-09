      parameter(nmax=10001,icol=9,icolmax=9)                                        
      dimension t(nmax),a(icolmax,nmax),amin(icolmax),amax(icolmax)
      dimension aa(icolmax,nmax)                
      logical iflush
      common /flush/ iflush
      data iflush /.true./
      do 100 i=1,nmax                                                 
        read(1,*,end=101) t(i),(a(j,i),j=1,icol) 
        do j=icol+1,icolmax
          read(1,*) a(j,i)
        enddo
        t(i)=t(i)/1000. 
        if(mod(i,1).eq.0) then                  
          write(21,202) t(i),a(5,i)                                   
          write(22,202) t(i),a(1,i)                                       
          write(23,202) t(i),a(2,i)                                       
          write(24,202) t(i),a(3,i)                                       
          write(25,202) t(i),a(6,i)                                       
        
          if(i.eq.2) then
            a(7,1)=a(7,2)
            write(26,202) t(1),a(7,i)                                       
            write(26,202) t(i),a(7,i)
          elseif(i.gt.2) then
            write(26,202) t(i),a(7,i)
          endif            
c          write(27,202) t(i),a(10,i)                                       
c          write(28,202) t(i),a(11,i)                                       
c          write(29,202) t(i),a(12,i)                                       
c          write(30,202) t(i),a(13,i)                                       
        endif                           
        if(i.gt.10 .and. mod(i,3000).eq.0) then                  
          write(21,202) -99999.,21.,0                                       
          write(22,202) -99999.,21.,0                                       
          write(23,202) -99999.,21.,0                                       
          write(24,202) -99999.,21.,0                                       
          write(25,202) -99999.,21.,0                                       
          write(26,202) -99999.,21.,0                                       
c          write(27,202) -99999.,21.,0                                       
c          write(28,202) -99999.,21.,0                                       
c          write(29,202) -99999.,21.,0                                       
c          write(30,202) -99999.,21.,0                                       
          write(21,*) 'tnsl'                                                
          write(22,*) ' total volume'                                        
          write(23,*) ' flotation volume'                                        
          write(24,*) ' areal extent'                                        
          write(25,*) ' average thickness'                                        
          write(26,*) ' dome elevation'                                        
c          write(27,*) ' avg basal temp  '                                         
c          write(28,*) ' total water  '                                         
c          write(29,*) ' percent water  '                                         
c          write(30,*) ' bed depression  '                                         
        endif
        npt=i                                                           
100   continue  
      print *,'incomplete read, nmax=',nmax                                                        
101   continue                                                          
      write(21,202) -99999.,21.,0                                       
      write(22,202) -99999.,21.,0                                       
      write(23,202) -99999.,21.,0                                       
      write(24,202) -99999.,21.,0                                       
      write(25,202) -99999.,21.,0                                       
      write(26,202) -99999.,21.,0                                       
c      write(27,202) -99999.,21.,0                                       
c      write(28,202) -99999.,21.,0                                       
c      write(29,202) -99999.,21.,0                                       
c      write(30,202) -99999.,21.,0                                       
      write(21,*) 'tnsl'                                                
      write(22,*) ' total volume'                                        
      write(23,*) ' flotation volume'                                        
      write(24,*) ' areal extent'                                        
      write(25,*) ' average thickness'                                        
      write(26,*) ' dome elevation'                                        
c      write(27,*) ' avg basal temp  '                                         
c      write(28,*) ' total water  '                                         
c      write(29,*) ' percent water  '                                         
c      write(30,*) ' bed depression  '                                         
      do j=1,icolmax                                                       
        amax(j)=-1.e30                                                  
        amin(j)=1.e30                                                   
      enddo                                                             
      do i=1,npt                                                        
        do j=1,icolmax
          aa(j,i)=a(j,i)                                                     
          amax(j)=max(amax(j),a(j,i))                                   
          amin(j)=min(amin(j),a(j,i))                                   
        enddo    
      enddo      
      print 1000,(j,amin(j),amax(j),j=1,icolmax)                           
1000  format(i5,2g13.6)                                                 
c      do j=1,icolmax
c        delta=(amax(j)-amin(j))/20.
c        amax(j)=amax(j)+delta
c        amin(j)=amin(j)-delta                                                       
c      enddo                                             
      do i=1,npt                                                        
        do j=1,icolmax                                                     
          if(amax(j).ne.amin(j)) then
            a(j,i)=(a(j,i)-amin(j))/(amax(j)-amin(j))
          else
            a(j,i)=.5
          endif
        enddo                                                           
      enddo                                                             
      call grstrt(800,1)                                               
      call window(t(1),t(npt),0.,1.)                                    
      icolor=2
1     continue
      if(iflush) call gflush
      write(*,*) ' -1-to plot one vs another'
      write(*,*) '  0-quit'
      write(*,*) '  1-vol'
      write(*,*) '  2-volt'
      write(*,*) '  3-area'
      write(*,*) '  4-snowline depression at pole'
      write(*,*) '  5-tnsl'
      write(*,*) '  6-avg ht'
      write(*,*) '  7-dome'
      write(*,*) '  8-lapse'
      write(*,*) '  9-snowline slope'
      write(*,*) ' 10-avg basal temp '
      write(*,*) ' 11-total water '
      write(*,*) ' 12-percent water '
      write(*,*) ' 13-bed depression '
      write(*,*) ' anything else to clear'
      read(*,*,end=999) j
      if(j.eq.0) goto 999
      if(j.eq.-1) then
        write(*,*) 'input components to plot vs ...'
        read(*,*) iplot,jplot
        do i=1,npt
          write(99,202) aa(iplot,i),aa(jplot,i)
        enddo
        write(99,202) -99999.,2.,0
        write(99,*) 'versus'
      endif                                               
      if(j.gt.icolmax .or. j.lt.1) then                                       
        call newpag                                                     
        icolor=2
        rewind 11                                                       
        goto 1                                                          
      endif                                                             
        call linclr(icolor)
        call move(t(1),a(j,1))                                          
        offset=0.0
        write(11,202) t(1),a(j,1)+offset
202     FORMAT(10X,G13.6,2X,G13.6,I13)                                  
        do i=2,npt                                                      
          call draw(t(i),a(j,i))                                        
          if(mod(i,3).eq.0) then                  
            write(11,202) t(i),a(j,i)+offset
          endif
        enddo                                                           
        write(11,202) -99999.,21.,icolor
        icolor=icolor+1
        if(j.eq.1) then                                                 
          write(11,2001) 'vflt ' ,amin(j),amax(j)
2001      format(a,2f7.2)
        elseif(j.eq.2) then                                             
          write(11,2002) 'vtot ' ,amin(j),amax(j)
2002      format(a,2f7.2)
        elseif(j.eq.3) then                                             
          write(11,2003) 'area ' ,amin(j),amax(j)
2003      format(a,2f7.2)
        elseif(j.eq.4) then                                             
          write(11,2004) 'delt ' ,amin(j),amax(j)
2004      format(a,2f7.0)
        elseif(j.eq.5) then                                             
          write(11,2005) 'tnsl ' ,amin(j),amax(j)
2005      format(a,2f7.2)
        elseif(j.eq.6) then                                             
          write(11,2006) 'thck ' ,amin(j),amax(j)
2006      format(a,2f7.0)
        elseif(j.eq.7) then                                             
          write(11,2007) 'dome ' ,amin(j),amax(j)
2007      format(a,2f7.0)
        elseif(j.eq.8) then                                             
          write(11,2008) 'lapse' ,amin(j),amax(j)
2008      format(a,2f7.3)
        elseif(j.eq.9) then                                             
          write(11,2009) 'snowline slope' ,amin(j),amax(j)
2009      format(a,2f8.5)
        elseif(j.eq.10) then                                             
          write(11,2010) 'avg basal temp' ,amin(j),amax(j)
2010      format(a,2f8.3)
        elseif(j.eq.11) then                                             
          write(11,2011) 'total water' ,amin(j),amax(j)
2011      format(a,2f8.3)
        elseif(j.eq.12) then                                             
          write(11,2012) 'percent water' ,amin(j),amax(j)
2012      format(a,2f8.3)
        elseif(j.eq.13) then                                             
          write(11,2013) 'bed depression' ,amin(j),amax(j)
2013      format(a,2f8.3)
        else                                                            
          rewind 11                                                     
        endif                                                           
      goto 1                                                            
999   continue                                                          
                                                                        
      call grstop1                                                      
                                                                        
      end                                                               
