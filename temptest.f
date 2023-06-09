      PARAMETER(NMAX=29999,MMAX=40)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension amass(11),tempa(mmax)
      common /heat/ hs
      dt=10.
      thick=1850.
      SLOPN=0d0
      adot=-0.5
      ttop=-20.
1     continue
      dt=10.
      SLOPN=0d0
      ttop=-20.
      print *,'input dt,THICK,heat',thick,adot
      read(*,*,end=999) thick,adot
      SLOPN=SLOPNN
      do i=1,mmax
         tempa(i)=ttop
      enddo
      geoflux=3.80e5
      fract=0.
      numnp=1
      i=1
      aeff=3.
      iplt=1
      IPLSTRT=0
      tbot=0.0
      wthick=0.D0
      wthresh=0d0
      TMELT=-8.7E-4*THICK
      CALL COLTEMP1(0.d0,I,NUMNP,AMASS,ADOT,TEMPA,
     &                      THICK,TMELT,TBOT,AEFF,0d0,IPLSTRT,IPLT,
     &                      SLOPN,1,BMELT,wthick,wthresh,geoflux,
     &                      fract)
      print *,'   steady state: tbot=',real(tbot),real(bmelt*1000)
      pause
      total=1000.
      nsteps=int(total/dt)
      istep=nint(100./dt)
      tempa(1)=ttop
      write(9,*) 0.,tempa(MMAX)
      do i=1,nsteps
c        tempa(1)=ttop-i*10./nsteps
C        adot=-i*3./nsteps
        CALL COLTEMP1(0.d0,I,NUMNP,AMASS,ADOT,TEMPA,
     &                      THICK,TMELT,TBOT,AEFF,DT,IPLSTRT,IPLT,
     &                      SLOPN,1,BMELT,wthick,wthresh,geoflux,
     &                      fract)
        
        if(mod(i,istep).eq.0) then
          print *,' time dependent:tbot=',
     &              real(tbot),real(i*dt),real(bmelt*1000)
          write(9,*) i*dt,tempa(MMAX)
        endif
        call wait(1)
      enddo
c       write(9,*) -99999.,0.
c       write(9,*) ' steady  ',thick
      pause
      total2=30000.
      nsteps=int(total2/dt)
      tempa(1)=ttop+10.
      do i=1,nsteps
        CALL COLTEMP1(0.d0,I,NUMNP,AMASS,ADOT,TEMPA,
     &                      THICK,TMELT,TBOT,AEFF,DT,IPLSTRT,IPLT,
     &                      SLOPNN,1,BMELT,wthick,wthresh,geoflux,
     &                      fract)
        wthick=max(0d0,wthick+bmelt*dt)
        if(mod(i,istep).eq.0) then
          print *,' time dep.:',
     &              real(tbot),real(i*dt),real(bmelt*1000),real(wthick)
          write(11,*) total+i*dt,real(wthick)
        endif
        call wait(1)
      enddo
      write(11,*) -99999.,0.
      write(11,*) ' time ',thick
      write(9,*) -99999.,0.
      write(9,*) ' time ',thick
      goto 1
999   continue
      call grstop1
      end



