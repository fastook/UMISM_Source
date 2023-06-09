      SUBROUTINE WMOVER0(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL, DT,
     &                 DTLOCAL,
     &                 ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &                 NZ, KZ, LM, TNEW,
     &                 BMELT, D, B, A, KA, ALPHAC, TOTALW, TOTALP,IPLOT)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NNMAX=MAXNUM)
      DIMENSION TESTVEC1(NNMAX),TESTVEC2(nnmax),TESTVEC3(nnmax)
      DIMENSION TESTVEC4(NNMAX),svec(nnmax),yyy(nnmax),phi(nnmax)
      DIMENSION ALPHAC(3)
      DIMENSION BMELT(NMAX),TNEW(NMAX)
      DIMENSION ITKODE(NMAX)
      DIMENSION KZ(NMAX)
      DIMENSION LM(5),HTICE(NMAX),DEPB(NMAX)
      DIMENSION KX(NMAX,4),X(NMAX),Y(NMAX),NTYPE(NMAX),D(NMAX),B(NMAX)
      REAL*8 A(NMAX,NZ),WTHICK(NMAX)
      INTEGER KA(NMAX,NZ+1)
      DIMENSION XI(2,9),ETA(2,9),W(2,9)
      REAL*4 TB(2)
c     REAL*4 TB(2),ETIME,DTIME
c     EXTERNAL ETIME,DTIME
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      SAVE NSTEP,ISTART,WSAVE,TSAVE,PSAVE
      DATA ISTART /0/, BIG /1D30/
      DATA RHOI /0.917D0/, RHOW /1.092D0/, GRAV /3.74d0/
      dtlocal=0.1
      NSTEP=dt/dtlocal
      TLOCAL=0
      do i=1,numnp
        phi(i)=rhoi*grav*(htice(i)-depb(i))+rhow*grav*depb(i)
      enddo
C ... CALL CHECKER THAT PUTS WATER ON/OFF ICE-FREE NODES
c ... (LAST ARG IS VALUE TO PUT ON ICE-FREE NODES, ZERO IT BEFORE
c     (turn on/off internally)
      CALL CHECKER(NMAX, NUMNP, NUMEL, NTYPE, KX, HTICE, DEPB, 
     &             WTHICK,0.0d0)
c ... form matrices
      H=DTLOCAL
      H2=0.5d0*H
      H6=H/6.d0
      do n=1,nstep
        TLOCAL=TLOCAL+DTLOCAL
c ..... save water thickness for convergence test ...
        do i=1,numnp
          svec(i)=wthick(i)
        enddo
c ..... FIRST RK STEP ...............................
        CALL FORMNT(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &             ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &             NZ, KZ, LM,
     &             BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &             AREAWET, AREADRY, AREATOT)
        IF(ISTART.EQ.0) THEN
          ISTART=1
          WSAVE=TOTALW
          IF(AREAWET.EQ.0.0) THEN
            TSAVE=0D0
          ELSE
            TSAVE=100.D0*TOTALW/AREAWET
          ENDIF
          IF(AREATOT.EQ.0.0) THEN
            PSAVE=0D0
          ELSE
            PSAVE=100.D0*AREAWET/AREATOT
          ENDIF
        ENDIF
        call formrhs(numnp,ka,kz,nz,a,b,d,phi,testvec1,rat1)
c ..... update water thickness for first RK step ........
        do i=1,numnp
          yyy(i)=max(0.d0,wthick(i)*(1-h2*alphac(3)))
          yyy(i)=yyy(i)+testvec1(i)*h2
          yyy(i)=max(0.d0,yyy(i))
        enddo
C ..... SECOND RK STEP .......................................
        CALL FORMNT(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &             ETA, XI, W, yyy, ITKODE, HTICE, DEPB,
     &             NZ, KZ, LM,
     &             BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &             AREAWET, AREADRY, AREATOT)
        call formrhs(numnp,ka,kz,nz,a,b,d,phi,testvec2,rat2)
c ..... update water thickness for SECOND RK step ........
        do i=1,numnp
          yyy(i)=max(0.d0,wthick(i)*(1-h2*alphac(3)))
          yyy(i)=yyy(i)+testvec2(i)*h2
          yyy(i)=max(0.d0,yyy(i))
        enddo
C ..... THIRD RK STEP .......................................
        CALL FORMNT(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &             ETA, XI, W, yyy, ITKODE, HTICE, DEPB,
     &             NZ, KZ, LM,
     &             BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &             AREAWET, AREADRY, AREATOT)
        call formrhs(numnp,ka,kz,nz,a,b,d,phi,testvec3,rat3)
c ..... update water thickness for THIRD RK step ........
        do i=1,numnp
          yyy(i)=max(0.d0,wthick(i)*(1-H*alphac(3)))
          yyy(i)=yyy(i)+testvec3(i)*H
          yyy(i)=max(0.d0,yyy(i))
        enddo
C ..... FOURTH RK STEP .......................................
        CALL FORMNT(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &             ETA, XI, W, yyy, ITKODE, HTICE, DEPB,
     &             NZ, KZ, LM,
     &             BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &             AREAWET, AREADRY, AREATOT)
        call formrhs(numnp,ka,kz,nz,a,b,d,phi,testvec4,rat4)
c ..... update water thickness for FOURTH RK step ........
        do i=1,numnp
          yyy(i)=max(0.d0,wthick(i)*(1-H*alphac(3)))
          yyy(i)=yyy(i)+
     &           h6*(testvec1(i)+
     &           2.d0*(testvec2(i)+testvec3(i))+testvec4(i))
          wthick(i)=max(0.d0,yyy(i))
        enddo
c ..... form DIFF for convergence test
        diff=0.d0
        do i=1,numnp
          diff=diff+abs(svec(i)-wthick(i))
        enddo
        diff=diff/numnp
        if(N.eq.1) then
          dsave=diff
        endif
        ratiod=diff/dsave
        if(mod(n,nstep/10).eq.0) then
          print *,n,real(diff),real(rat4),real(ratiod),
     &              REAL(1D-9*(TOTALW))
        endif
        if(ratiod.lt.1e-2) goto 10
      enddo
c      print *,'didnt converge to steady state'
10    continue
      RATIOW=(TOTALW-WSAVE)/DT
      TOTALT=100.D0*TOTALW/AREAWET
      TOTALP=100.D0*AREAWET/AREATOT
      RATIOT=(TOTALT-TSAVE)/DT
      RATIOP=(TOTALP-PSAVE)/DT
      IF(IOTOGG) THEN
        WRITE(LIST(IPAGE+1),*) ISTEP,
     &       '********************************************'
        WRITE(LIST(IPAGE+2),*) 
     &         'TOTAL WATER (KM**3)=',REAL(1D-9*TOTALW),
     &         REAL(1D-9*RATIOW),REAL(1D-9*RATIOW*DT)
C       WRITE(LIST(IPAGE+1),*) 'WET AREA    (KM**2)=',REAL(1D-6*AREAWET)
C       WRITE(LIST(IPAGE+1),*) 'DRY AREA    (KM**2)=',REAL(1D-6*AREADRY)
C       WRITE(LIST(IPAGE+1),*) 'TOTAL AREA  (KM**2)=',REAL(1D-6*AREATOT)
        WRITE(LIST(IPAGE+3),*) 'PERCENT WET        =',REAL(TOTALP),
     &         REAL(RATIOP),REAL(RATIOP*DT)
        WRITE(LIST(IPAGE+4),*) 'AVG THICK   ( CM  )=',REAL(TOTALT),
     &         REAL(RATIOT),REAL(RATIOT*DT)
        WRITE(LIST(IPAGE+5),*) 'used',n,' steps, ratio',real(ratiod),
     &                         real(TLOCAL)
        WRITE(LIST(IPAGE+6),*) 
     &       '********************************************'
        IPAGE=IPAGE+6
      ENDIF
      WSAVE=TOTALW
      PSAVE=TOTALP
      TSAVE=TOTALT
C ... CALL CHECKER THAT PUTS WATER ON/OFF ICE-FREE NODES
c ... (LAST ARG IS VALUE TO PUT ON ICE-FREE NODES, 0.1 IT AFTER
c     (turn on/off internally)
      CALL CHECKER(NMAX, NUMNP, NUMEL, NTYPE, KX, HTICE, DEPB, 
     &             WTHICK,1d-5)






      end
C-----------------------------------------------------------------
      SUBROUTINE FORMNT(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL,
     &                 ETA, XI, W, WTHICK, ITKODE, HTICE, DEPB,
     &                 NZ, KZ, LM,
     &                 BMELT, D, B, A, KA, ALPHAC, TOTALW,
     &                 AREAWET, AREADRY, AREATOT)
      IMPLICIT REAL*8(A-H,O-Z)
C FORM STIFFNESS MATRIX
      DIMENSION ALPHAC(3)
      DIMENSION BMELT(NMAX)
      DIMENSION ITKODE(NMAX)
      DIMENSION KZ(NMAX)
      DIMENSION LM(5),HTICE(NMAX),DEPB(NMAX)
      DIMENSION KX(NMAX,4),X(NMAX),Y(NMAX),NTYPE(NMAX),D(NMAX),B(NMAX)
      REAL*8 A(NMAX,NZ),WTHICK(NMAX)
      INTEGER KA(NMAX,NZ+1)
      DIMENSION P(5),S(5,5),DD(5)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2)
      DIMENSION PSI(4),DPSI(4,2)
      DIMENSION DWSIX(9),DWSIY(9)
      DIMENSION WSI(4),DWSI(4,2)
      DIMENSION XY(2,4),XI(2,9),ETA(2,9),W(2,9)
      DIMENSION ALPHA(2),BETA(2)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      LOGICAL UPWIND
      DATA BIG /1.D30/
      DATA ASCALE /1.D0/
C **********************************************************************
C     FORM CONDUCTIVITY MATRIX FOR COMPLETE BODY
C **********************************************************************
C
C ... ZERO OUT APPROPRIATE ARRAYS ...
      DO I=1,NUMNP
        KZ(I)=1
        KA(I,1)=I
        D(I)=0.0D0
        B(I)=0.0D0
        KA(I,NZ+1)=0
        DO J=1,NZ
          KA(I,J)=0
          A(I,J)=0.0D0
        ENDDO
      ENDDO
C
C ... BEGIN LOOP OVER ALL THE ELEMENTS ...
      TOTALW=0.0D0
      AREAWET=0.0D0
      AREADRY=0.0D0
      AREATOT=0.0D0
      avg=0.
      iavg=0
      DO 100 N=1,NUMEL
        IF(NTYPE(N).EQ.1) THEN
          NODEN=4
          NINT=9
        ELSE
          NODEN=3
          NINT=4
          PRINT *,'TRIANGLES NOT ALLOWED!'
          STOP
        ENDIF
        DO I=1,NODEN
          LM(I)=KX(N,I)
        ENDDO
C
C ..... FORM ELEMENT CONDUCTIVITY MATRIX
        DO I=1,NODEN
          DD(I)=0.0D0
          P(I)=0.0D0
          DO J=1,NODEN
            S(I,J)=0.0D0
          ENDDO
        ENDDO
C
        I=LM(1)
        J=LM(2)
        K=LM(3)
        L=LM(4)
        XY(1,1)=X(I)
        XY(1,2)=X(J)
        XY(1,3)=X(K)
        XY(2,1)=Y(I)
        XY(2,2)=Y(J)
        XY(2,3)=Y(K)
        IF(NODEN.EQ.4) THEN
          XY(1,4)=X(L)
          XY(2,4)=Y(L)
        ENDIF
C ..... USE FOLLOWING TO GENERATE CENTROID VALUES FOR GRADIENTS AND
C       OTHER MATERIAL PROPERTIES 
        IF(.TRUE.) THEN
          CALL FESHAPE(NTYPE(N),0D0,0D0,PSI,DPSI)
          CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
          TWTHIK=0.D0
          BMELTN=0.D0
          DGDX=0.D0
          DGDY=0.D0
          SURFX=0
          SURFY=0
          BEDX=0
          BEDY=0
          DO I=1,NODEN
            SURFX=SURFX+HTICE(LM(I))*DPSIX(I)
            SURFY=SURFY+HTICE(LM(I))*DPSIY(I)
          ENDDO
          SURFXY=SQRT(SURFX**2+SURFY**2)
          DO I=1,NODEN
            BEDX=BEDX+DEPB(LM(I))*DPSIX(I)
            BEDY=BEDY+DEPB(LM(I))*DPSIY(I)
            BMELTN=BMELTN+BMELT(LM(I))*PSI(I)
            TWTHIK=TWTHIK+WTHICK(LM(I))*PSI(I)
            IF(WTHICK(LM(I)).GT.0.001) THEN
C       ....  COEFFICIENT COMES FROM ALLEY'S PAPER
              RNPRES=5.0D-5*(HTICE(LM(I))-DEPB(LM(I)))*SURFXY/
     &              WTHICK(LM(I))
            ELSE
              RNPRES=0.D0
            ENDIF
c effective pressure turned off below...
            RNPRES=0.D0
            GGG=-(HTICE(LM(I))+0.09D0*DEPB(LM(I)))-RNPRES
            DGDX=DGDX+GGG*DPSIX(I)
            DGDY=DGDY+GGG*DPSIY(I)
          ENDDO
          DGXY=SQRT(DGDX**2+DGDY**2)
          SURFXY=SQRT(SURFX**2+SURFY**2)
          BEDXY=SQRT(BEDX**2+BEDY**2)
        ENDIF
C........... TO HERE
c
C ..... FORM ELEMENT MATRIX AND VECTORS
C
C ..... BEGIN INTEGRATION POINT LOOP
        DO L=1,NINT
          CALL FESHAPE(NTYPE(N),XI(NTYPE(N),L),ETA(NTYPE(N),L),PSI,DPSI)
C
C ......  GENERATE DEIRVATIVES, JUST NON-UPWINDING, COPY INTO UPWINDING
          CALL DERIVE(XY,DXDS,DPSI,DETJ,DSDX,DPSIX,DPSIY)
          DO I=1,NODEN
            WSI(I)=PSI(I)
            DWSI(I,1)=DPSI(I,1)
            DWSI(I,2)=DPSI(I,2)
            DWSIX(I)=DPSIX(I)
            DWSIY(I)=DPSIY(I)
          ENDDO
C ....... ACCUMULATE INTEGRATION POINT VALUES OF INTEGRALS
          FAC=DETJ*W(NTYPE(N),L)
          TOTALW=TOTALW+TWTHIK*FAC
          IF(TWTHIK.GT.0.0) THEN
            AREAWET=AREAWET+FAC
          ELSE
            AREADRY=AREADRY+FAC
          ENDIF
          AREATOT=AREATOT+FAC
C ....... LATEST NEW WAY ..................
          GRADG=SQRT(DGDX**2+DGDY**2)
          IF(TWTHIK.LE.0.01 .or. .true.) THEN
            IPP=2
            IQQ=1
            ACONST=ASCALE*ALPHAC(1)*TWTHIK**(IPP+1)*GRADG**(IQQ-1)
          ELSE
            RPP=.5D0
            RQQ=.7D0
            ACONST=ASCALE*ALPHAC(1)*TWTHIK**(RPP+1.d0)*GRADG**(RQQ-1)
          ENDIF
          avg=avg+aconst*gradg
          iavg=iavg+1
          DO I=1,NODEN
C ......... THIS IS LUMPED CAPACITANCE MATRIX
            DD(I)=DD(I)+WSI(I)*FAC
C ........
C ........  INTEGRATION POINT VALUES FROM FEM INTERPOLATIOON ....
C ......... THIS IS LOAD VECTOR
            P(I)=P(I)+BMELTN*WSI(I)*FAC
C
            DO J=1,NODEN
              TERM1=-ACONST*(dwsix(i)*dpsix(j)+dwsiy(i)*dpsiy(j))
              TERM2=0.0
              TERM3=0.0
              S(I,J)=S(I,J)+(TERM1+TERM2+TERM3)*FAC
            ENDDO
          ENDDO
        ENDDO
C       IF(S(1,1).NE.0) THEN
C         PRINT *,N
C         DO II=1,4
C         PRINT *,(REAL(S(II,JJ)),JJ=1,4)
C         ENDDO
C       ENDIF
C
C ..... ADD ELEMENT CONDUCTIVITY TO COMPLETE CONDUCTIVITY MATRIX
C
        DO L=1,NODEN
          I=LM(L)
C ....... THIS IS LUMPED CAPACITANCE MATRIX
          D(I)=D(I)+DD(L)
C
C ....... THIS IS LOAD VECTOR
          B(I)=B(I)+P(L)
          DO M=1,NODEN
            J=LM(M)
            IF(I.EQ.J) THEN
              A(I,1)=A(I,1)+S(L,M)
              KA(I,1)=I
            ELSE
              DO K=2,KZ(I)
                IF(KA(I,K).EQ.J) THEN
                  A(I,K)=A(I,K)+S(L,M)
                  GOTO 99
                ENDIF
              ENDDO
              KZ(I)=KZ(I)+1
              A(I,KZ(I))=S(L,M)
              KA(I,KZ(I))=J
            ENDIF
99          CONTINUE
          ENDDO
          KA(I,NZ+1)=KZ(I)
        ENDDO
100   CONTINUE
C ... END LOOP OVER ALL THE ELEMENTS ...
C      IF(NFIX.EQ.0) THEN
C        PRINT *,'NO NODES ARE FIXED....'
C        PAUSE
C      ELSE
C        PRINT *,NFIX,' NODES ARE FIXED....',BIG
C      ENDIF
C ... REMOVE HERE, THIS WRITES OUT MATRIX
C      REWIND 73
C      WRITE(73,*) NUMNP
C      DO I=1,NUMNP
C        WRITE(73,*) KA(I,10)
C        DO J=1,KA(I,10)
C          WRITE(73,*) I,KA(I,J),A(I,J)
C        ENDDO
C      ENDDO
C      WRITE(73,*) (B(I),I=1,NUMNP)
C ... TO HERE ******
C      IF(IOTOGG) THEN
C        WRITE(LIST(IPAGE+1),*) 
C     &     'TOTAL WATER (KM**3)=',REAL(1E-9*TOTALW)
C        IPAGE=IPAGE+1
C       PRINT *,'WET AREA    (KM**2)=',REAL(1E-6*AREAWET)
C       PRINT *,'DRY AREA    (KM**2)=',REAL(1E-6*AREADRY)
C       PRINT *,'TOTAL AREA  (KM**2)=',REAL(1E-6*AREATOT)
C       PRINT *,'PERCENT WET        =',REAL(100.*AREAWET/AREATOT)
C       PRINT *,'AVG THICK   ( CM  )=',REAL(100.*TOTALW/AREAWET)
C      ENDIF
      END


      subroutine formrhs(numnp,ka,kz,nz,a,b,d,phi,testvec,rat)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      DIMENSION TESTVEC(NMAX)
      DIMENSION KZ(NMAX)
      DIMENSION phi(NMAX)
      DIMENSION KX(NMAX,4),D(NMAX),B(NMAX)
      REAL*8 A(NMAX,NZ)
      INTEGER KA(NMAX,NZ+1)
c ..... stiffness matrix * pressure potential
      do i=1,numnp
        testvec(i)=0.d0
        do jj=1,kz(i)
          j=ka(i,jj)
          testvec(i)=testvec(i)+a(i,jj)*phi(j)
        enddo
      enddo
c ... add load vector
      rat=0
      irat=0
      do i=1,numnp
        if(b(i).gt.0) then
          rat=rat+testvec(i)/b(i)
          irat=irat+1
        endif
        testvec(i)=testvec(i)+b(i)
      enddo
      if(irat.gt.0) then
        rat=rat/irat
      endif
c ... scale with 1/d (M^-1)
      do i=1,numnp
        testvec(i)=testvec(i)/d(i)
      enddo
      end
