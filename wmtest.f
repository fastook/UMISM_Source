      include "parameter.h"
      PARAMETER( NMAX=MAXNUM,NZ=9,NZ1=NZ+1,NSMAX=MAXTIME,N3=NMAX*3,
     &           N23=N3*2)
      IMPLICIT REAL*8(A-H,O-Z)
c
      DIMENSION WWW(N3),THICKL(NMAX),WRATE(N3,2)
      DIMENSION WWWORIG(NMAX),WDIFF(NMAX)
      CHARACTER HED*80,SCRTCH*80,IADJ*2
      REAL*8 A(NMAX,NZ)
      COMMON /MAT/  A,KA(NMAX,NZ1),NUMNP
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL,TNSLBASE
      COMMON /VELOS/ ASCAL,UMAX,VTHRESH,INORM
      DIMENSION IFIT(6), AMASS(11), B(NMAX), X(NMAX), Y(NMAX), D(NMAX),
     &          BOLD(NMAX),FLUX(NMAX),ADIAG(NMAX),ITKODE(NMAX),
     &          KODE(NMAX), CONST(NMAX), ACON(NMAX), LM(5), ITYPE(NMAX),
     &          KX(NMAX,4), IDT(NMAX), ADOTB(NMAX), 
     &          ADOT(NMAX), BDROCK(NMAX), FLOWA(NMAX), SLDGB(NMAX),
     &          PSURF(NMAX), PPSURF(NMAX), FRACT(NMAX),
     &          CNEW(NMAX), QHOLD(NMAX), HTICE(NMAX), THICK(NMAX),
     &          HFIT(NMAX), IBFLUX(NMAX,2), BFLUX(NMAX), DEPB(NMAX),
     &          Q(NMAX),TEMP(NMAX),TBED(NMAX),VEL(NMAX,3),
     &          UNDEPB(NMAX), SLOPE(4,NMAX), SLOPN(4,NMAX), KZ(NMAX),
     &          WTHICK(NMAX),WTELEM(NMAX),GEOFLUX(NMAX),
     &          CALV(NMAX),PCALV(NMAX)
      DIMENSION AFUDGE(NMAX),ADP(NMAX),ADM(NMAX),ADC(NMAX),BMELT(NMAX)
      DIMENSION TTIME(NSMAX),VVOL(NSMAX),AAREA(NSMAX),TTBOT(NSMAX),
     &          TTNSL(NSMAX),TWATER(NSMAX),PWATER(NSMAX),WWWMIN(NSMAX)
      DIMENSION NTYPE(NMAX), NNODE(NMAX),
     &          XI(2,9), ETA(2,9), W(2,9),
     &          AADOT(NMAX), AFRACT(NMAX), AFLOWA(NMAX),
     &          ABDRCK(NMAX), ASLDGB(NMAX)
      DIMENSION ALPHAC(3)
c     DIMENSION DPSIX(9), DPSIY(9), DXDS(2,2), DSDX(2,2),
c    &          PSI(4), DPSI(4,2), CNST(NMAX),
c    &          XY(2,4)
      LOGICAL CTOGG,WTOGG,ITOGG,IOTOGG
      INTEGER BTOGG
      COMMON /TOGGLES/ CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
C ... TIMER STUFF, FOR SGI ONLY ...
      REAL*4 TB(2)
c     REAL*4 TB(2),ETIME,DTIME
c     EXTERNAL ETIME,DTIME
C ... TIMER STUFF, FOR SGI ONLY ...
      REAL*8 T(NMAX),TNEW(NMAX)
      EXTERNAL WARMING
      DATA HFIT /NMAX*0.0/
      DATA ASCAL /1.D0/, UMAX /100.D0/, VTHRESH /1000.D0/
      DATA INORM /1/
      DATA WWW /N3*0.D0/
      DATA WRATE /N23*0.D0/
      IOTOGG=.false.
      WRITE(*,123) ' TIME BEFORE READ ',ETIME(TB),DTIME(TB)
      read(67) MMAX, X, Y, KX, NTYPE, NUMNP, NUMEL, DT,
     &                   DTLOCAL,
     &                   ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &                   MZ, KZ, LM, TNEW,
     &                   BMELT, D, B, A, KA, ALPHAC, TOTALW, TOTALP
      WRITE(*,123) ' TIME AFTER READ ',ETIME(TB),DTIME(TB)
c--------------------------
      XMIN=1.E30
      YMIN=XMIN
      XMAX=-XMIN
      YMAX=-YMIN
      DO I=1,NUMNP
        XMAX=MAX(XMAX,x(I)*0.001)
        YMAX=MAX(YMAX,y(I)*0.001)
        XMIN=MIN(XMIN,x(I)*0.001)
        YMIN=MIN(YMIN,y(I)*0.001)
      ENDDO
      IF(XMAX-XMIN.GT.YMAX-YMIN) THEN
        YMAX=YMIN+XMAX-XMIN
      ELSE
        XMAX=XMIN+YMAX-YMIN
      ENDIF
      DELX=(XMAX-XMIN)/5.
      DELY=(YMAX-YMIN)/5.
      time=0
      dt=10.0
      dtlocal=dt
      CALL GRSTRT(800,800)
      CALL WINDOW(REAL(XMIN-DELX),REAL(XMAX+DELX),
     &              REAL(YMIN-DELY),REAL(YMAX+DELY))
c--------------------------
      CALL CONTR(NMAX,NUMEL,x,y,KX,
     &                 WTHICK,-0.0999999999D0,1.000000001D0,0.1D0,
     &                 -1000.d0,1000.d0,-1000.d0,1000.d0)
      DO L=1,100
        TIME=TIME+DT
        ipage=1
        WRITE(*,*) 'TIME=',TIME
c        WRITE(*,123) ' TIME BEFORE WMOVER ',ETIME(TB),DTIME(TB)
        CALL WMOVER(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL, DT,
     &                   DTLOCAL,
     &                   ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &                   NZ, KZ, LM, TNEW,
     &                   BMELT, D, B, A, KA, ALPHAC, TOTALW, TOTALP)
         dtlocal=dt
c        WRITE(*,123) ' TIME AFTER WMOVER ',ETIME(TB),DTIME(TB)
123   FORMAT(A25,T30,1PG13.6,G13.6)
        call newpag
        CALL CONTR(NMAX,NUMEL,x,y,KX,
     &                 HTICE,-2999.999D0,5500.D0,250.D0,
     &                 -1000.d0,1000.d0,-1000.d0,1000.d0)
        CALL CONTR(NMAX,NUMEL,x,y,KX,
     &                 WTHICK,-0.0999999999D0,1.000000001D0,0.1D0,
     &                 XMIN,XMAX,YMIN,YMAX)
        CALL LINCLR(1)                                                    
        CALL MOVE(REAL(XMIN),REAL(YMIN))                                  
        CALL DRAW(REAL(XMIN),REAL(YMAX))                                  
        CALL DRAW(REAL(XMAX),REAL(YMAX))                                  
        CALL DRAW(REAL(XMAX),REAL(YMIN))                                  
        CALL DRAW(REAL(XMIN),REAL(YMIN))                                  
      enddo
      call grstop
      end
