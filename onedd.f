      implicit real*8(a-h,o-z)
1     print *,'input adot,tsurf,thick'
      read(*,*,end=999) adot,tsurf,thick
      call coltemp(adot,tsurf,thick,tbot,aeff)
      print *,tbot,aeff
      goto 1
999   end
c**********************************
      subroutine coltemp(adot,tsurf,thick,tbot,aeff)
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
      dimension wwww(nmax),told(nmax),at(nmax),uc(nmax),ddata(9)
      dimension u(nmax),tttt0(nmax)
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
      do i=1,npt                                                    
        tttt0(i)=tttt(i)                                                
      enddo                                                          
c --- loop on iteration step                                                     
      do 900 idt=1,ndt 
        time=idt            
        call window(real(tmin),real(tmax),real(xmin),real(xmax))
c.....come in here for iteration on vertical velocity
700   continue                                                
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
c*****jump back to here for 2nd pass, fixed basal bc
800   continue
        if(ipass.eq.2) then
          ipass=1
          aaad(npt)=1.                                                    
          aaau(npt-1)=0.                                                  
          aaal(npt-1)=0.                                                  
          ffff(npt)=tttt(npt)                                             
          ffff(npt-1)=ffff(npt-1)-rkn2*ffff(npt)                          
        endif
        call tri(npt,aaal,aaad,aaau,ffff,tttt) 
        if(tttt(npt).gt.0.0) then
c          print *,tttt(npt),' going around again'
          tttt(npt)=0.
          ipass=2
          goto 800
        endif
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
        sigma0=p1-((rk11-cccd(1))*tttt(1)+
     &         rk12*tttt(2)-cccu(1)*tttt0(2))
        sigma0=-sigma0/delt                                             
        sigmalc=pn-((rkn1-cccd(npt))*tttt(npt)+
     &         rkn2*tttt(npt-1)-cccu(npt-1)*tttt0(npt-1))   
        dtemp=(tttt(npt-1)-tttt(npt))/(xxx(npt-1)-xxx(npt))
        sigmabc=-rkx(npt-1)*dtemp
        sigmalc=-sigmalc/delt                                             
        bmelt=-(sigmal-sigmabc)/rkx(npt-1)
c       print *,sigmal,bmelt,sigmabc
        if(wwwn.eq.-99999.) then
          if(abs(wwww(npt)-bmelt).gt.toler) then
            print *,wwww(npt),bmelt,' going around again'
            wwww(npt)=bmelt
c...........loop back to beginning after recalculating vertical velocity
            slope=(wwww(npt)-wwww(1))/(xxx(npt)-xxx(1))
            do i=2,npt-1
              wwww(i)=wwww(1)+slope*(xxx(i)-xxx(1))
c             write(7,*) 'wwww',i,wwww(i)
            enddo
            goto 700
          endif
        endif
        do i=1,npt-1                                                
          dtdx=(tttt(i+1)-tttt(i))/delx(i)                              
          flux=-rkx(i)*dtdx                                             
          write(11,202) .5*(xxx(i+1)+xxx(i)),flux                       
        enddo                                                        
        write(11,202) -99999.,-1.,0                                        
202     format(10x,g13.6,2x,g13.6,i13)                                  
        write(11,*) 'sigma0=',sigma0,sigmalc                             
        write(12,202) time,sigma0                                       
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
          call wait(500000)
          call grstop1
          tbot=tttt(npt)
          print *,' tbot = ',tbot
          return
        endif
        call velo(time,npt,5d-3,xxx,tttt,rho,at,u,uc,
     &            bmelt,ddata,aeff)                                                         
900   continue                                                          
      do i=1,npt                                                    
        tttt0(i)=tttt(i)                                                
      enddo                                                          
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
c*************************************************
      subroutine velo(time,npts,dh,x,t,rho,at,u,uc,
     &                bmelt,ddata,aeff)
      parameter(nmax=100)
      implicit real*8(a-h,o-z)
      save uold,ucold,ubarold,avgold
      dimension ddata(9)
      dimension x(npts),t(npts),rho(npts),at(npts)
      dimension u(npts),uc(npts),uold(nmax),ucold(nmax)
      data uold /nmax*0d0/, ucold /nmax*0d0/, avg /0d0/,ubarold /0d0/
      if(npts.gt.nmax) then
        print *,'increase nmax in velo, nmax=',nmax
        stop
      endif
      a0=1.830357
      b0=7.738778e-2
      a1=2.633955
      b1=4.0990233e-2
      call hard(a0,b0,a1,b1)
      do i=1,npts
        if(t(i).ge.-10.) then
          at(i)=a0*exp(-b0*t(i))
        else
          at(i)=a1*exp(-b1*t(i))
        endif
      enddo
      top=x(1)
      g=9.81e-5
      avgat=0.
      avgt=0.
      do i=npts-1,1,-1
        dx=x(i+1)-x(i)
        avgat=avgat-0.5*(at(i+1)+at(i))*dx
        avgt=avgt-0.5*(t(i+1)+t(i))*dx
      enddo
      avgat=avgat/(x(1)-x(npts))
      avgt=avgt/(x(1)-x(npts))
      ddata(5)=avgt
c      print *,'average a=',avgat
      avgr=0.
      do i=npts-1,1,-1
        dx=x(i+1)-x(i)
        avgr=avgr-0.5*(rho(i+1)+rho(i))*dx
      enddo
      avgr=avgr/(x(1)-x(npts))
c      print *,'average rho=',avgr
      ubar=.4*(avgr*g*dh/avgat)**3*(x(1)-x(npts))**4
c      print *,'ubar=',ubar
      u(npts)=0.
      do i=npts-1,1,-1
        sum=u(i+1)
        v1=2.*(rho(i+1)*g*dh*(top-x(i+1))/at(i+1))**3
        v2=2.*(rho(i)*g*dh*(top-x(i))/at(i))**3
        dx=x(i+1)-x(i)
        u(i)=u(i+1)-0.5*(v1+v2)*dx
      enddo
      avg=0.
      do i=npts-1,1,-1
        dx=x(i+1)-x(i)
        avg=avg-0.5*(u(i+1)+u(i))*dx
      enddo
      avg=avg/(x(1)-x(npts))
      ddata(6)=avg
      write(13,202) avgt,avg
      write(14,202) time,avg
      write(15,202) time,avgt
      write(16,202) time,t(npts)
202   format(10x,g13.6,2x,g13.6,i13)
      aeff=(.4*(910.*g*dh)**3*(x(1)-x(npts))**4/avg)**(1./3.)
      print 1002,time,real(avgt),real(aeff),real(avg),
     &                real(t(npts)),real(bmelt)
1002  format(f10.0,2f10.3,f10.2,f10.3,f10.6)
      ddata(7)=aeff
c      print *,'a-effective=',aeff
      c=.5*(910.*g*dh/aeff)**3
      thick=x(1)-x(npts)
      do i=1,npts
        uc(i)=c*(thick**4-(top-x(i))**4)
      enddo
      ymin=x(npts)
      ymax=x(1)
      umin=1e30
      umax=-1e30
c       print 1001,'x','t','a','u','uc'
1001  format(1x,5a13)
      do i=1,npts
        umin=min(umin,u(i))
        umax=max(umax,u(i))
        umin=min(umin,uc(i))
        umax=max(umax,uc(i))
c         print 1000,x(i),t(i),at(i),u(i),uc(i)
1000    format(1x,1p5g13.6)
      enddo
c      print *,umin,umax,ymin,ymax
c      call grstrt(400,400)
      call window(real(0.),real(500.),real(ymin),real(ymax))
      call linclr(0)
      call move(real(uold(1)),real(x(1)))
      do i=2,npts
        call draw(real(uold(i)),real(x(i)))
      enddo
      call move(real(avgold),real(ymin))
      call draw(real(avgold),real(ymax))
      call linclr(0)
      call move(real(ucold(1)),real(x(1)))
      do i=2,npts
        call draw(real(ucold(i)),real(x(i)))
      enddo
c      call move(real(ubarold),real(ymin))
c      call draw(real(ubarold),real(ymax))
      call linclr(1)
      call move(real(u(1)),real(x(1)))
      do i=2,npts
        call draw(real(u(i)),real(x(i)))
      enddo
      call move(real(avg),real(ymin))
      call draw(real(avg),real(ymax))
      call linclr(2)
      call move(real(uc(1)),real(x(1)))
      do i=2,npts
        call draw(real(uc(i)),real(x(i)))
      enddo
c      call move(real(ubar),real(ymin))
c      call draw(real(ubar),real(ymax))
      avgold=avg
c      ubarold=ubar
      do i=1,npts
        uold(i)=u(i)
        ucold(i)=uc(i)
      enddo
c      call grstop1
      end
c*****************************************
      subroutine hard(b0p,slpbp,bb0p,slpbbp)
      implicit real*8(a-h,o-z)
      parameter(nmax=11)
      save ipass,b0,slpb,bb0,slpbb
      real*8 t(nmax),a(nmax),b(nmax)
      data t /0.0,    -5.0,   -10.0,   -15.0,   -20.0,   -25.0,   -30.0,
     &      -35.0,   -40.0,   -45.0,   -50.0/
      data ipass /0/
      if(ipass.eq.0) then
        ipass=1
        a(1)=5.3e-15
        a(2)=1.7e-15
        a(3)=5.2e-16
        a(4)=3.1e-16
        a(5)=1.8e-16
        a(6)=1.0e-16
        a(7)=5.4e-17
        a(8)=2.9e-17
        a(9)=1.5e-17
        a(10)=7.7e-18
        a(11)=3.8e-18
        conver=1.6e-2/5.2e-16
        do i=1,nmax
          temp=a(i)*conver
          b(i)=1./temp
          b(i)=b(i)**(1./3.)
c          print *,t(i),a(i),temp,b(i)
        enddo
        b0=b(1)
        slpb=-(log(b(1))-log(b(3)) )/(t(1)-t(3))
        slpbb=-(log(b(3))-log(b(11)))/(t(3)-t(11))
        rlnbb0=log(b(3))+slpbb*t(3)
        bb0=exp(rlnbb0)
      endif
c         print *,'b0,bb0',b0,bb0
c         print *,'r,rr  ',slpb,slpbb
      b0p=b0
      slpbp=slpb
      bb0p=bb0
      slpbbp=slpbb
      end
