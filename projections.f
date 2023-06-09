c======================================================
      SUBROUTINE SETRIG(POLE,RLONG0)
      implicit real*8(a-h,o-z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD                               
      if(POLE.gt.0) then
        call setproj( 90.d0,rlong0); return
      else
        call setproj(-90.d0,rlong0); return
      endif
c---------------------------------------------
      RADIUS=6370.D0                                                    
      PI=4.D0*ATAN(1.D0)                                                
      RADIUS=2.0D4/PI                                              
      CIRCUM=2.D0*PI*RADIUS                                             
      RKMPDEG=CIRCUM/360.D0                                             
      RADPDEG=PI/180.D0                                                 
      DEGPRAD=180.D0/PI                                                 
      END                                                               
c======================================================
      SUBROUTINE POLREC(RLAT,RLONG,X,Y)                                 
      implicit real*8(a-h,o-z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD                               
      real*8 dlat,dlong,xee,ynn
      dlat=dble(rlat);dlong=dble(rlong)
      call ll2xyB(dlat,dlong,xee,ynn)
      x=real(xee);y=real(ynn)
c     x=real(0.001*xee);y=real(0.001*ynn)
c     print *,rlat,rlong,x,y
      return
c---------------------------------------------
      X=1000.D0*(90.D0-RLAT)*RKMPDEG*COS(RLONG*RADPDEG)                 
      Y=1000.D0*(90.D0-RLAT)*RKMPDEG*SIN(RLONG*RADPDEG)                 
      END                                                               
c======================================================
      SUBROUTINE RECPOL(X,Y,RLAT,RLONG)                                 
      implicit real*8(a-h,o-z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD                               
      real*8 dlat,dlong,xee,ynn
c     xee=1000.d0*dble(x);ynn=1000.d0*dble(y)
      xee=1.d0*dble(x);ynn=1.d0*dble(y)
      call xy2llB(xee,ynn,dlat,dlong)
      rlat=real(dlat);rlong=real(dlong) 
c     print *,rlat,rlong
c     pause
      return
c---------------------------------------------
      IF(X.EQ.0.D0) THEN                                                
        IF(Y.GE.0.D0) THEN                                              
          THETA=90.D0                                                   
        ELSE                                                            
          THETA=270.D0                                                  
        ENDIF                                                           
      ELSE                                                              
        THETA=DEGPRAD*ATAN(Y/X)                                         
        IF(X.LT.0.D0) THEN                                              
          THETA=180.D0+THETA                                            
        ELSEIF(Y.LT.0.D0) THEN                                          
          THETA=360.D0+THETA                                            
        ELSEIF(THETA.GT.360.D0) THEN                                    
          THETA=THETA-360.D0                                            
        ENDIF                                                           
      ENDIF                                                             
      R=.001D0*SQRT(X**2+Y**2)                                          
      RLAT=90.D0-R/RKMPDEG                                              
      RLONG=THETA                                                       
      END                                                               
c======================================================
      subroutine setproj(rlat0,rlong00)
      implicit real*8(a-h,o-z)
      common /projcons/ aaa,eee,pi,radpdeg,rlong0,rlatf,phi0,
     &                   rk0,fe,fn,io
c ... EARTH ...
c     parameter(aaa=6377397.155d0,loverf=299.15281d0,eee=0.08169683d0)
c ...  MARS ...
c     parameter(aaa=3396200.d0,loverf=169.81d0,eee=0.1083660d0)
c     parameter(k0=0.9999079d0)
c     parameter(k0=1.d0)

      io=0
      phi0=rlat0
c     aaa=6378137.d0; eee=0.081819191d0;rk0=0.972769012891835d0 !EARTH
      aaa=3396200.d0; eee=0.1083660d0;rk0=1.0d0                 ! MARS
      rlong0= rlong00; rlatf= 71.d0
      pi=4.d0*atan(1.d0); radpdeg=pi/180.d0
      fe=0d6; fn=0d6
      end
c======================================================
      subroutine ll2xyB(rlat,rlong,xee,ynn)
      implicit real*8(a-h,o-z)
      common /projcons/ aaa,eee,pi,radpdeg,rlong0,rlatf,phi0,
     &                   rk0,fe,fn,io
      logical npole
      npole =(phi0.gt.0d0); rlamb0=rlong0*radpdeg
      phif  =sign(rlatf*radpdeg,phi0)  
      phi   =rlat  *radpdeg; rlamb =rlong *radpdeg
      term1=1d0+eee*sin(phif);term2=1d0-eee*sin(phif)
      denom1=(term1/term2)**(eee/2)
      term1=1d0+eee*sin(phi);term2=1d0-eee*sin(phi)
      denom2=(term1/term2)**(eee/2)
      if(npole) then
        tf= tan(pi/4d0-phif/2d0)*denom1
        ttt=tan(pi/4d0- phi/2d0)*denom2
      else
        tf= tan(pi/4d0+phif/2d0)/denom1
        ttt=tan(pi/4d0+ phi/2d0)/denom2
      endif;
      rmf=cos(phif)/sqrt(1d0-eee**2*sin(phif)**2)
      denom3=sqrt((1+eee)**(1+eee)*(1-eee)**(1-eee))
      rk0=rmf*denom3/(2d0*tf)
      rho=2d0*aaa*rk0*ttt/denom3
      theta=rlamb-rlamb0
      if(npole) then
        de=rho*sin(theta); dn=-rho*cos(theta)
      else
        de=rho*sin(theta); dn=rho*cos(theta)
      endif
      xee=fe+de; ynn=fn+dn
      end
c======================================================
      subroutine xy2llB(xee,ynn,rlat,rlong)
      implicit real*8(a-h,o-z)
      common /projcons/ aaa,eee,pi,radpdeg,rlong0,rlatf,phi0,
     &                   rk0,fe,fn,io
      logical npole
      phif =sign(rlatf*radpdeg,phi0)  
      npole=(phi0.gt.0d0); rlamb0=rlong0*radpdeg
      term1=1d0+eee*sin(phif);term2=1d0-eee*sin(phif)
      denom1=(term1/term2)**(eee/2)
      if(npole) then
        tf=tan(pi/4d0-phif/2d0)*denom1
      else
        tf=tan(pi/4d0+phif/2d0)/denom1
      endif
      rmf=cos(phif)/sqrt(1d0-eee**2*sin(phif)**2)
      denom3=sqrt((1+eee)**(1+eee)*(1-eee)**(1-eee))
      rk0=rmf*denom3/(2d0*tf)
      rho=sqrt((xee-fe)**2+(ynn-fn)**2)  
      ttt=rho*denom3/(2d0*aaa*rk0)
      if(npole) then
        chi=pi/2d0-2d0*atan(ttt)
      else
        chi=2d0*atan(ttt)-pi/2d0
      endif  
      term1=eee**2/2d0+5d0*eee**4/24d0+eee**6/12d0+13d0*eee**8/360d0
      term2=7d0*eee**4/48d0+29d0*eee**6/240d0+811d0*eee**8/11520d0
      term3=7d0*eee**6+81d0*eee**8/1120d0
      term4=4279d0*eee**8/161280d0
      phi=chi+term1*sin(2d0*chi)+term2*sin(4d0*chi)+
     &        term3*sin(6d0*chi)+term4*sin(8d0*chi)
      if(xee.eq.fe) then
        rlamb=rlamb0
      elseif(npole) then
        rlamb=rlamb0+atan2((xee-fe),(fn-ynn))
      else
c       rlamb=rlamb0+atan2((xee-fe),(ynn-fn))+pi/2.d0
        rlamb=rlamb0+atan2((xee-fe),(ynn-fn))
      endif
      rlat=phi/radpdeg; rlong=rlamb/radpdeg
      if(rlong.lt.0) rlong=rlong+360d0
      if(rlong.gt.360) rlong=rlong-360d0
      end
