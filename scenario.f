      subroutine scenario
      common /exper/ iexper
      print *,' 1-        mass(7) coef      0  10.5 k  -900 -2100'
      print *,'--------------------------------------------------'
      print *,' 2-<  70 k mass(7) coef1     0    20 k     0 -2300'
      print *,'   >  70 k mass(7) coef1    70    90 k -2300     0'
      print *,'--------------------------------------------------'
      print *,' 3-        mass(7) coef      0  20.5 k     0 -2300'
      print *,' 3-        mass(8) coef    -25 -14.5 k  5600  3800'
      print *,'--------------------------------------------------'
      print *,' 4-<20.5 k mass(7) coef1     0  20.5 k     0 -2300'
      print *,'   <61.5 k mass(7) coef1  20.5    41 k -2300 -1000'
      print *,'   >61.5 k mass(7) coef1  61.5    82 k -2300     0'
      print *,'--------------------------------------------------'
      print *,' 5-<20.5 k mass(7) coef1     0  20.5 k -2000     0'
      print *,'   <61.5 k mass(7) coef1  20.5    41 k     0 -3000'
      print *,'   >61.5 k mass(7) coef1  61.5    82 k     0 -2000'
      print *,'--------------------------------------------------'
      print *,' 6-        mass(7) coef     -3   6.4 k     0 -2000'
      print *,'--------------------------------------------------'
      print *,' 7-        mass(7) coef     -3     7 k -1000 -4000'
c
      print *,'--------------------------------------------------'
      print *,' 8-        tnsl    coef  -9.55 10.55 k     0   -20'
      print *,'--------------------------------------------------'
      print *,' 9-        tnsl    coef      0     1 k    -6   -14'
      print *,'--------------------------------------------------'
      print *,'10- experiment for northern hemisphere ice sheet cycle (smooth)'
      print *,'   <100 k  tnsl    coef      0    20 k     0   -20'
      print *,'   <120 k  tnsl    coef    100   120 k   -20     0'
      print *,'   <128 k  tnsl    coef    120   128 k     0    -5'
      print *,'   >128 k  tnsl = -5'
      print *,'--------------------------------------------------'
      print *,'11- experiment for northern hemisphere ice sheet cycle (linear)'
      print *,'   < 40 k  tnsl    coef      0    20 k    -7   -20'
      print *,'   < 41 k  tnsl    coef     40    41 k   -20     0'
      print *,'   < 51 k  tnsl    coef     51    52 k     0    -5'
      print *,'   > 51 k  tnsl = -5'
      print *,'--------------------------------------------------'
      print *,'12 < 10 k  tnsl    coef      0    10 k     0    -6'
      print *,'   < 20 k  tnsl    coef     10    20 k    -6     0'
      print *,'   > 20 k  tnsl = 0'
      print *,'--------------------------------------------------'
      print *,'13- 20,000 year cycle of both tnsl and lapse rate'
      print *,'   < 10 k  tnsl    coef      0    10 k    -7   -20'
      print *,'           acom    coef      0    10 k    -5  -9.5'
      print *,'   < 20 k  tnsl    coef     10    20 k   -20    -5'
      print *,'           acom    coef     10    20 k  -9.5    -5'
      print *,'   > 20 k  tnsl = -5'
      print *,'           acom = -5'
      print *,'--------------------------------------------------'
      print *,'14- abrupt change'
      print *,'   <  5 k  tnsl    coef      0     5 k    -7   -20'
      print *,'   < 15 k  tnsl = -20'
      print *,'   < 20 k  tnsl    coef     15    20 k   -20     0'
      print *,'   > 20 k  tnsl = 0'
      print *,'--------------------------------------------------'
      print *,'15- very abrupt change'
      print *,'   < 30 k  tnsl = -20'
      print *,'   > 30 k  tnsl =   0'
      print *,'--------------------------------------------------'
      print *,' 16- linear deglaciation scenario'
      print *,'           tnsl    coefl     0     5 k   -20     0'
      print *,'--------------------------------------------------'
      print *,' 17- read in TLIST.DATA or user defined'
      read(*,*) iexper
      end
