      implicit real*8(a-h,o-z)
1     print *,'input adot,tsurf,thick'
      read(*,*,end=999) adot,tsurf,thick
      call coltemp(adot,tsurf,thick,tbot)
      print *,tbot
      goto 1
999   end
c**********************************
      subroutine coltemp(adot,tsurf,thick,tbot)
c******************************************
c this solves for steady state vertical temperature profile
c with solved for basal temperature
c if basal temp > pressure melting, will
c resolve with basal bc fixed at pressure melting
c and calculate basal melting
c******************************************
      implicit real*8(a-h,o-z)
      parameter(nmax=100,big=1.e20)
      dimension xxx(nmax),tttt(nmax),cond(nmax),heat(nmax),rho(nmax)
      dimension qqq(nmax),source(nmax),delx(nmax)
      dimension rkx(nmax),ffff(nmax),aaad(nmax),aaau(nmax),aaal(nmax)
      dimension cccd(nmax),cccu(nmax),cccl(nmax)
      dimension wwww(nmax),told(nmax)
      dimension u(nmax)
      dimension lmap(nmax)
      print *,' adot = ',adot
      print *,' tsurf= ',tsurf
      print *,' thick= ',thick
      npt=50
      ndt=10
      wwww(1)=0.
      sigmal=1.2046
      qqqs=0.
      wwww(1)=-adot
      tttt(1)=tsurf
c
      ltot=0
      do i=1,npt/2
        lmap(i)=i
        ltot=ltot+lmap(i)
      enddo
      do i=npt/2+1,npt
        lmap(i)=lmap(i-1)-1
        ltot=ltot+lmap(i)
      enddo
      xxx(1)=thick
      xxx(npt)=0.
      dx=(xxx(npt)-xxx(1))/(ltot)
      do i=2,npt
        tttt(i)=tttt(1)
        xxx(i)=xxx(i-1)+lmap(i-1)*dx
        qqq(i)=0.
        source(i)=0.
        told(i)=0.
      enddo
c
      toler=.001
      factor=1.
      xmax=-1e30
      xmin=1e30
9     continue
      do i=1,npt
        xmax=max(xmax,xxx(i))
        xmin=min(xmin,xxx(i))
        depth=xxx(1)-xxx(i)
        rho(i)=density(depth)
        cond(i)=conduct(tttt(i),rho(i))
        heat(i)=spheat(tttt(i))
        rkx(i)=cond(i)/rho(i)/heat(i)
        if(i.gt.1) delx(i-1)=xxx(i)-xxx(i-1)
        aaad(i)=0.
        aaau(i)=0.
        aaal(i)=0.
        cccd(i)=0.
        cccu(i)=0.
        cccl(i)=0.
        ffff(i)=0.
      enddo
c           
c          
      call grstrt(600,600)
      tmin=-60.
      tmax=60.
      call window(real(tmin),real(tmax),real(xmin),real(xmax))
      call linclr(1)
      call move(0.,real(xmin))
      call draw(0.,real(xmax))
      call mrkclr(2)
      wwww(npt)=0. 
c --- vertical velocity distribution: linear
      slope=(wwww(npt)-wwww(1))/(xxx(npt)-xxx(1))
      do i=2,npt-1
        wwww(i)=wwww(1)+slope*(xxx(i)-xxx(1))
      enddo
      do i=1,npt
        qqq(i)=qqqs
      enddo
c --- loop on iteration step                                                     
      do 900 idt=1,ndt             
        call window(real(tmin),real(tmax),real(xmin),real(xmax))
c ----- form load vector 
        ffff(1)=.5*u(1)*qqq(1)*delx(1)+source(1)                             
        ffff(npt)=.5*u(npt-1)*qqq(npt-1)*delx(npt-1)+
     &             source(npt)-sigmal          
        p1=ffff(1)                                                      
        pn=ffff(npt)                                                    
        do i=2,npt-1                                                
          ffff(i)=.5*(u(i)*qqq(i)*delx(i)+
     &                u(i-1)*qqq(i-1)*delx(i-1))+source(i)      
        enddo                                                        
c ----- end form load vector  
c                                                
c ----- form stiffness matrix                                                 
        aaad(1)=rkx(1)/delx(1)+wwww(1)*.5                               
        aaad(npt)=rkx(npt-1)/delx(npt-1)-wwww(npt)*.5                   
        do i=2,npt-1                                                
          aaad(i)=rkx(i-1)/delx(i-1)+rkx(i)/delx(i)                     
          aaad(i)=aaad(i)+.5*(wwww(i-1)-wwww(i))                        
        enddo                                                        
        do i=1,npt-1                                                
          aaau(i)=-rkx(i)/delx(i)                                       
          aaal(i)=aaau(i)-wwww(i)*.5                                    
          aaau(i)=aaau(i)+wwww(i)*.5                                    
        enddo                                                        
c ----- fix boundary condition for node 1                                     
        rk11=aaad(1)                                                    
        rk12=aaau(1)                                                    
        rk21=aaal(1)                                                    
        aaad(1)=1.                                                      
        aaau(1)=0.                                                      
        aaal(1)=0.                                                      
        ffff(1)=tttt(1)                                                 
        ffff(2)=ffff(2)-rk21*ffff(1)                                    
c ----- fix boundary condition for node npt
        rkn1=aaad(npt)                                                  
        rkn2=aaau(npt-1)                                                
        rk2n=aaal(npt-1)    
        ipass=1                                            
c ----- the matrix solution
        call tri(npt,aaal,aaad,aaau,ffff,tttt) 
c --------------------------
        aaad(1)=rk11                                                    
        aaau(1)=rk12                                                    
        aaal(1)=rk21                                                    
        aaad(npt)=rkn1                                                  
        aaau(npt-1)=rkn2                                                
        aaal(npt-1)=rk2n                                                
c ------following for rkx function of temperature                             
        do i=1,npt-1                                                
c --------rkxi= function of temperature                                
          depth=xxx(1)-xxx(i)  
          rho(i)=density(depth)                                             
          cond(i)=conduct(tttt(i),rho(i))
          heat(i)=spheat(tttt(i))
          rkxi=cond(i)/rho(i)/heat(i)                                     
          rkx(i)=((factor-1.)*rkx(i)+rkxi)/factor                      
        enddo                                                        
          call linclr(1)                                             
          call move(real(told(1)),real(xxx(1)))                           
          do i=2,npt                                                  
            call draw(real(told(i)),real(xxx(i)))                         
          enddo                                                        
          call linclr(2)                                             
          call move(real(tttt(1)),real(xxx(1)))                           
          do i=2,npt                                                  
            call draw(real(tttt(i)),real(xxx(i)))                         
          enddo                                                        
          diff=0d0
          do i=1,npt
            diff=diff+(told(i)-tttt(i))**2
            told(i)=tttt(i)
          enddo
          diff=sqrt(diff)
        print 1103,idt,real(tttt(npt)),real(diff)
1103    format(i10,f10.3,f10.6)
        if(diff.lt.toler) then
          call wait(100000)
          call grstop1
          tbot=tttt(npt)
          print *,' tbot = ',tbot
          return
        endif
900   continue                                                          
      print *,'failed to converge'
      call wait(10000)
      call grstop1
      tbot=-999.
      end         
c*********************************************
      subroutine wait(i)
      do n=1,i
        x=sin(real(i))
      enddo
      end                                                      
c*************************************************
      subroutine tri(n,aaa,ddd,ccc,bbb,xxx) 
      parameter(nmax=100)                            
      implicit real*8(a-h,o-z)                                          
      dimension aaa(n),bbb(n),ccc(n),ddd(n),xxx(n)                      
      dimension aal(nmax),bbl(nmax),ccl(nmax),ddl(nmax)
      if(n.gt.nmax) then
        print *,'increase nmax in sub. tri',nmax
        stop
      endif
      do i=1,n-1
        aal(i)=aaa(i)
        bbl(i)=bbb(i)
        ccl(i)=ccc(i)
        ddl(i)=ddd(i)
      enddo
      ddl(n)=ddd(n)
      bbl(n)=bbb(n)
      do i=2,n                                                        
        xmult=aaa(i-1)/ddl(i-1)                                         
        ddl(i)=ddl(i)-xmult*ccl(i-1)                                    
        bbl(i)=bbl(i)-xmult*bbl(i-1)                                    
      enddo                                                          
      xxx(n)=bbl(n)/ddl(n)                                              
      n1=n-1                                                            
      do ii=1,n1                                                      
        i=n-ii                                                            
        xxx(i)=(bbl(i)-ccl(i)*xxx(i+1))/ddl(i)                          
      enddo                                                          
      return                                                            
      end                                                               
c*************************************************
      function density(depth)
      implicit real*8(a-h,o-z)
      parameter(rhoi=910.d0,rhos=350.d0,crho=0.0374d0)
        if (depth .eq. 0.) then
          density=rhos
        elseif ((depth.gt.0.).and.(depth.le.55.)) then
          density=rhoi-((rhoi-rhos)*exp(-crho*depth))
        elseif ((depth.gt.55.) .and.(depth.lt.65.)) then
          density=7.d0*depth+455.d0
        elseif (depth .ge. 65.) then
          density=rhoi
        endif
      end
c*************************************************
      function conduct(tttt,rho)
      implicit real*8(a-h,o-z)
      if(tttt.gt.0.) then
        ttt=0.
      else
        ttt=tttt
      endif
      rrr=rho*.001
      conduct=3.1536d06*
     &          (( 4.20102d0-0.00754145d0*ttt)
     &          -(26.48767d0-0.048779d0  *ttt)*rrr
     &          +(57.31865d0-0.141127d0  *ttt)*rrr**2
     &          -(29.55156d0-0.053672d0  *ttt)*rrr**3)
      end
c*************************************************
      function spheat(tttt)
      implicit real*8(a-h,o-z)
      if(tttt.gt.0.) then
        ttt=0.
      else
        ttt=tttt
      endif
      spheat=1000.d0*(0.494d0+0.00138d0*ttt)                          
      end

