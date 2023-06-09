      IMPLICIT REAL*8(A-H,O-Z)
      temp1 = 0d0
      temp2 = 226d0-273.16d0
      rho = 917d0
      print *,'input geothermal in mW/m2'
      read *,geomWpm2
      c=3.77d5/50d0 ! umism units/mWm2
      print 1,'mW/m2',geomWpm2
      flux = geomWpm2*c
      print 1,'umism units',flux
      dtdz1 = 1000d0*flux/conduct(temp1,rho)
      dtdz2 = 1000d0*flux/conduct(temp2,rho)
      thick1 = -temp2/dtdz1
      thick2 = -temp2/dtdz2
      print 1,'temp',temp1,'gradient',dtdz1,'thickness',thick1
      print 1,'temp',temp2,'gradient',dtdz2,'thickness',thick2
1     format(1x,3(a10,1pg13.6))
      end
C*************************************************
      FUNCTION CONDUCT(TTTT,RHO)
      IMPLICIT REAL*8(A-H,O-Z)
C ... UNITS: CALORIES/M//YR/DEGREE
      IF(TTTT.GT.0.) THEN
        TTT=0.D0
      ELSE
        TTT=TTTT
      ENDIF
      RRR=RHO*.001D0
      CONDUCT=3.1536D06*
     &          (( 4.20102D0-0.00754145D0*TTT)
     &          -(26.48767D0-0.048779D0  *TTT)*RRR
     &          +(57.31865D0-0.141127D0  *TTT)*RRR**2
     &          -(29.55156D0-0.053672D0  *TTT)*RRR**3)
      END

      
