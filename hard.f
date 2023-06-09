      implicit real*8(a-h,o-z)
      dimension ttt(12),aaa(12),aa1(12),aa2(12),bbb(12)
      data ttt /0.,-2.,-5.,-10.,-15.,-20.,
     &         -25.,-30.,-35.,-40.,-45.,-50./
c 1 bar=10**6 dyne/cm**2
c 1 dyne = force(cgs) = gr-cm/s**2
c 1 bar = 10**6 gr/cm/s**2
c 1 Pa = 1 Nt/m**2
c 1 Nt = kg-m/s**2
c 1 Pa = kg/m/s**2
c 1 bar = 10**5 Pa = 100 kPa
c a1,a2 in units of 1/Pa**3/yr
c huybrect flow law espdot=aaa * tau**n
c mine is           espdot=(tau/b)**n
      data aaa /680.d0,240.d0,160.d0,49.d0,29.d0,17.d0,9.4d0,
     &          5.1d0,2.7d0,1.4d0,.73d0,.36d0/
      conv=365.d0*24.d0*60.d0*60.d0
      rrr=8.314d0
      q1=60000.d0
      a1=1.14d-5
      q2=139000.d0
      a2=5.47d10
      do i=1,12
        aaa(i)=aaa(i)*1d-17
        aa1(i)=aaa(i)*conv*1d-9
        temp=ttt(i)+273.15d0
        if(temp.le.263.15d0) then
          aa2(i)=a1*exp(-q1/rrr/temp)
        else
          aa2(i)=a2*exp(-q2/rrr/temp)
        endif
        bbb(i)=aa2(i)*1d15
        bbb(i)=1.d0/bbb(i)
        bbb(i)=bbb(i)**(1.d0/3.d0)
      enddo
      do i=1,12
        print 100,ttt(i),bbb(i),hardness(ttt(i))
100     format(1x,f6.0,3x,1p5e13.6)
      enddo
      end
      function hardness(ttt)
      implicit real*8(a-h,o-z)
      data rrr /8.314d0/
      data q1 /60000.d0/
      data a1 /1.14d-5/
      data q2 /139000.d0/
      data a2 /5.47d10/
c
      temp=ttt+273.15d0
      if(temp.le.263.15d0) then
        aaa=a1*exp(-q1/rrr/temp)
      else
        aaa=a2*exp(-q2/rrr/temp)
      endif
      bbb=aaa*1d15
      bbb=1./bbb
      hardness=bbb**(1.d0/3.d0)
      end
