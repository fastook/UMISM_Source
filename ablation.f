      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /LOCAL/ SIGMA,TMA,TMS,PI2
      tnsl=0.0
      ncol=0
      do rlat=40.,90.,10.
        ncol=ncol+1
      enddo
      nrow=0
      do elev1=0.,3000.,100.
        nrow=nrow+1
      enddo
1     print *,'elevation, latitude'
      read(*,*,end=999) elev1,rlat
      ic=0
c      do rlat=60.,80.,10.
c      do elev1=0.,3000.,10.
        ic=ic+1
        Z=MAX(ELEV1,20.*(RLAT-65))
        TS=49.13-0.007992*Z-0.7576*RLAT+TNSL
        TSUMMER=30.38-0.006277*ELEV1-0.3262*RLAT+TNSL
C ... THIS IS HUYBRECTS, VERY SLOW...
        sigma=5.d0
        PI2=8.*ATAN(1.)
        TMA=TS
        TMS=TSUMMER
c*********************************************************************
c the 3.5 in the following is a fudge to get the present config to be
c stable. the temp TS passed out for the surface temperature does not
c include it....
        TMA=49.13-0.007992*Z-0.7576*RLAT+TNSL+3.5                  
        TMS=30.38-0.006277*ELEV1-0.3262*RLAT+TNSL+3.5   
        CALL QUAD2D(0.D0,365.D0,PDDH)
        PDDH=PDDH/SIGMA/SQRT(PI2)
        call huyb(elev1,rlat,pdd,tss,TSSUMER,TNSL)
        ABLH=0.8*PDDH
C ... THIS IS SIMPLE SINE, FASTER BUT TWO-TIMES TOO LARGE...
c        PDD=0.
c        AMP=TSUMMER-TS
c        DO IT=1,365
c          TD=TS+AMP*COS(PI2*IT/365.)
c          IF(TD.GT.0.0) PDD=PDD+TD
c        ENDDO
c ... this makes it agree better with huybrects' slower method...
c        PDD=PDD*0.5
        ABL=0.8*PDD
        print *,'---------------------------------'
        print *,' latitude             ',rlat
        print *,' elev,z               ',elev1,z
        print *,' annual temp          ',ts,TSS
        print *,' summer temp          ',tsummer,TSSUMER
        print *,' positive degree days ',pddH,pdd
        print *,' ablation rate        ',ablH,abl
        write(11,*) rlat,elev1,ablH
        write(12,*) rlat,elev1,abl
        write(13,*) real(ic),real(ablh-abl)
c      enddo
c      enddo
      write(13,*) -99999.,0.
      write(13,*) 'test'
      goto 1
999   continue
      END                                                               
