      subroutine polrec(ll,rlat0,rlong0,FE,FN,rlat,rlong,x,y)
c-----------------------------------------------------------------------
c OGP Surveying and Positioning Guidance Note number 7, part 2 – January 2009
c-----------------------------------------------------------------------
c EPSG dataset coordinate operation method code 9809
c The Stereographic projection may be imagined to be a projection of the 
c earth's surface onto a plane in contact with the earth at a single 
c tangent point from a projection point at the opposite end of the 
c diameter through that tangent point.
c This projection is best known in its polar form and is frequently used 
c for mapping polar areas where it complements the Universal Transverse 
c Mercator used for lower latitudes. Its spherical form has also been 
c widely used by the US Geological Survey for planetary mapping and the 
c mapping at small scale of continental hydrocarbon provinces. In its 
c transverse or oblique ellipsoidal forms it is useful for mapping limited 
c areas centred on the point where the plane of the projection is regarded 
c as tangential to the ellipsoid., e.g. the Netherlands. The tangent point 
c is the origin of the projected coordinate system and the meridian through 
c it is regarded as the central meridian. In order to reduce the scale error 
c at the extremities of the projection area it is usual to introduce a scale 
c factor of less than unity at the origin such that a unitary scale factor 
c applies on a near circle centred at the origin and some distance from it.
c The coordinate conversion from geographical to projected coordinates is 
c executed via the distance and azimuth of the point from the centre point 
c or origin. For a sphere the formulas are relatively simple. For the 
c ellipsoid the parameters defining the conformal sphere at the tangent 
c point as origin are first derived. The conformal latitudes and longitudes 
c are substituted for the geodetic latitudes and longitudes of the spherical 
c formulas for the origin and the point.
c-----------------------------------------------------------------------
      implicit none
      integer ipass,ll,m
      real*8 rlat0,rlong0,FE,FN,rlat,rlong,x,y
      real*8 radpdeg,aaa,lovef,eee,k0,loverf
      real*8 phi0,lambda0,phi,lambda,clambda,clambda0,chi,chi0
      real*8 sinphi0,term,rho0,vvv0,RRR,n,c,S1,S2,w1,w2,w,Sa,Sb
      real*8 sinphi,B,sinchi0
      real*8 g,h,i,j,pi,psi,psii,phii,change
      logical io
      data io /.false./
      save pi,radpdeg,ipass
      data ipass /0/
c ... EARTH ...
c     parameter(aaa=6377397.155d0,loverf=299.15281d0,eee=0.08169683d0)
c ...  MARS ...
      parameter(aaa=3396200.d0,loverf=169.81d0,eee=0.1083660d0)
c     parameter(k0=0.9999079d0)
      parameter(k0=1.d0)
      if(ipass.eq.0) then
        ipass=1;pi=4.d0*atan(1.d0);radpdeg=pi/180.d0
      endif
      phi0=radpdeg*rlat0
      lambda0=radpdeg*rlong0
      sinphi0=sin(phi0)
      term=1d0-(eee*sinphi0)**2
      vvv0=aaa/sqrt(term)
      rho0=vvv0*(1d0-eee*eee)/term
      RRR = sqrt(rho0*vvv0)
      n = sqrt(1d0+((eee*eee*cos(phi0)**4)/(1d0-eee*eee)))
      S1 = (1d0+sinphi0)/(1d0-sinphi0) 
      S2 = (1d0-eee*sinphi0)/(1d0+eee*sinphi0)
      w1 = (S1*(S2)**eee)**n 
      sinchi0 =(w1-1d0)/(w1+1d0) 
      c=(n+sinphi0)*(1d0-sinchi0)/((n-sinphi0)*(1d0+sinchi0))
      w2 =c*w1
      chi0 =asin((w2-1d0)/(w2+1d0))
      clambda0 =lambda0
      if(ll.eq.0) then
c ...   phi = latitude, lambda = longitude of point to (x,y)
        phi=radpdeg*rlat
        lambda=radpdeg*rlong
        sinphi=sin(phi)
        clambda=n*(lambda-clambda0 )+clambda0
        Sa=(1d0+sinphi)/(1d0-sinphi) 
        Sb=(1d0-eee*sinphi)/(1d0+eee*sinphi)
        w=c*(Sa*(Sb)**eee)**n 
        chi=asin((w-1d0)/(w+1d0))
        B = (1d0+sin(chi)*sin(chi0)+
     &      cos(chi)*cos(chi0)*cos(clambda-clambda0))
        x=FE+2d0*RRR*k0*cos(chi)*sin(clambda-clambda0)/B 
        y=FN+2d0*RRR*k0*(sin(chi)*cos(chi0)-
     &      cos(chi)*sin(chi0)*cos(clambda-clambda0))/B 
        if(io) then
          print *,'phi0,lambda0',phi0,lambda0
          print *,'phi ,lambda ',phi ,lambda 
          print *,'rho0        ',6374588.71-rho0
          print *,'vvv0        ',6390710.613-vvv0
          print *,'RRR         ',6382644.571-vvv0
          print *,'n           ',1.0004758571-n
          print *,'S1          ',8.509582274-S1
          print *,'S2          ',0.878790173-S2
          print *,'w1          ',8.428769183-w1
          print *,'sinchi0     ',0.787883237-sinchi0
          print *,'c           ',1.007576465-c
          print *,'w2          ',8.492629457-w2
          print *,'chi0        ',0.909684757-chi0
          print *,'clambda     ',0.104724841-clambda
          print *,'chi         ',0.924394997-chi
          print *,'B           ',1.999870665-B
          print *,'x           ',196105.283-x
          print *,'y           ',557057.739-y
        endif
        return
      else
        g=2d0*RRR*k0*tan(pi/4d0-chi0/2d0) 
        h=4d0*RRR*k0*tan(chi0)+g 
        i = atan((x-FE)/(h+(y-FN))) 
        j = atan((x-FE)/(g-(y-FN)))-i
        chi=chi0+2d0*atan(((y-FN)-(x-FE)*tan(j/2d0))/(2d0*RRR*k0))
        clambda=j+2d0*i+clambda0
        lambda=(clambda-clambda0)/n+clambda0
        psi=0.5*log((1d0+sin(chi))/(c*(1d0-sin(chi))))/n
        phii=2d0*atan(exp(psi))-pi/2d0
        if(io) then
          print *,'g           ',4379954.188,g
          print *,'h           ',37197327.960,h
          print *,'i           ',0.001102255,i
          print *,'j           ',0.008488122,j
          print *,'chi         ',0.924394767,chi
          print *,'clambda     ',0.10472467,clambda
          print *,'lambda      ',0.104719584,lambda
          print *,'psi         ',1.089495123,psi
          print *,'phii        ',0.921804948,phii
        endif
        do m=2,10
          psii=log((tan(phii/2d0+pi/4d0))*((1d0-eee*sin(phii))/
     &             (1d0+eee*sin(phii)))**(eee/2d0))
          change=-(psii-psi)*cos(phii)*(1d0-(eee*sin(phii))**2)/
     &            (1d0-eee*eee)
          phii=phii+change
          if(io) print *,m,phii,psii
          if(abs(change).lt.1d-8) then
            rlat=phii/radpdeg
            rlong=lambda/radpdeg
            return
          endif
        enddo
        print *,'problems, did not converge'
        stop
      endif
      end
