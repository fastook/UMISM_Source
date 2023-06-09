      SUBROUTINE READN(MXX,HED,NUMNP,NUMEL,NUMGBC,NDT,INTER,DT,
     &                 KODE,X,Y,HTICE,ADOT,ADOTB,FRACT,PSURF,RHOI,
     &                 RHOW,RHOR,DEPB,
     &                 BDROCK,UNDEPB,FLOWA,ACON,SLDGB,TEMP,ITYPE,
     &                 AFUDGE,GEOFLUX,THICK,KX,CONST,IBFLUX,BFLUX,
     &                 QHOLD,NTYPE,NNODE,NCOL,AADOT,AFRACT,ABDRCK,
     &                 PPSURF,AFLOWA,ASLDGB,IDT,AMASS,NUMCOL,NUMLEV,
     &                 CALV,PCALV,WWW,WRATE,THICKL,WWWORIG,TSORIG,
     &                 TIME,ACC,ABLAT)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
C
C ... READ FROM SPLIT DATA SETS, PRODUCED BY UNPACK, REASSEMBLED BY PACK
C
      COMMON /SNOW/ SNOLIN,SNO(20)
      COMMON /TTNSL/ TTIME(NSMAX),TLIST(NSMAX),NTNSL
      COMMON /EXPER/ TPERIOD,TINIT,TVSTART,TVFINAL,IEXPER
      LOGICAL CTOGG,ITOGG,IOTOGG
      LOGICAL file67,file76
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      CHARACTER HED*80
      DIMENSION AMASS(11)
      DIMENSION KODE(MXX),X(MXX),Y(MXX),HTICE(MXX),ADOT(MXX),FRACT(MXX),
     &          PSURF(MXX),BDROCK(MXX),FLOWA(MXX),SLDGB(MXX),
     &          THICK(MXX),ADOTB(MXX),ITYPE(MXX),
     &          KX(MXX,4),CONST(MXX),IBFLUX(MXX,2),BFLUX(MXX),
     &          QHOLD(MXX),IDT(MXX),GEOFLUX(MXX),
     &          NNODE(MXX),NTYPE(MXX),AADOT(MXX),AFRACT(MXX),
     &          ABDRCK(MXX),AFUDGE(MXX),
     &          PPSURF(MXX),AFLOWA(MXX),ASLDGB(MXX),UNDEPB(MXX),
     &          LM(4),TEMP(MXX),ZERO(4),ACON(MXX),
     &          CALV(MXX),PCALV(MXX),WWW(3*MXX),WRATE(3*MXX,2),
     &          THICKL(MXX),DEPB(MXX),WWWORIG(MXX),TSORIG(MXX)
      DIMENSION ACC(MXX),ABLAT(MXX)
C ... TIMER STUFF, FOR SGI ONLY ...
      REAL*4 TB(2)
c     REAL*4 TB(2),ETIME,DTIME
c     EXTERNAL ETIME,DTIME
      DATA ZERO /4*0D0/
123   FORMAT(A25,T30,1PG13.6,G13.6)

      CALL SCENARIO(0.d0)
      DO I=1,NSMAX
        TTIME(I)=0.d0
        TLIST(I)=0.d0
        STIME(I)=0.d0
        SLIST(I)=0.d0
      ENDDO
      NSEAL=0
      DO I=1,NSMAX
        READ(75,*,END=889) STIME(I),SLIST(I)
        NSEAL=I
      ENDDO
      PRINT *,'INCREASE NUMBER OF STIME,STLIST',NSMAX
      STOP
889   CONTINUE
      IF(NSEAL.GT.0) THEN
        PRINT *,'THERE IS A SEA LEVEL LIST FOR THIS DATA SET'
        PRINT *,'THERE ARE',NSEAL,' TIMES IN THE LIST'
        PRINT *,'STARTING AT',STIME(1),' AND ENDING AT',STIME(NSEAL)
      ENDIF
      NTNSL=0
      DO I=1,NSMAX
        READ(72,*,END=888) TTIME(I),TLIST(I)
        NTNSL=I
      ENDDO
      PRINT *,'INCREASE NUMBER OF TTIME,TTLIST',NSMAX
      STOP
888   CONTINUE
      IF(NTNSL.GT.0) THEN
        PRINT *,'THERE IS A TNSL LIST FOR THIS DATA SET'
        PRINT *,'THERE ARE',NTNSL,' TIMES IN THE LIST'
        PRINT *,'STARTING AT',TTIME(1),' AND ENDING AT',TTIME(NTNSL)
      ENDIF
      DO I=1,20
        SNO(I)=0.d0
      ENDDO
C
C ... READ INPUT HEADER
      READ(30,1000,END=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &                        INTER,DT
1000  FORMAT (A80,/,7I6,F8.0)
C
C ... READ INPUT GRID, THINGS THAT NEVER CHANGE
      READ(31) HED
      READ(31) (KODE(I),I=1,NUMNP)
      READ(31) (X(I),I=1,NUMNP)
      READ(31) (Y(I),I=1,NUMNP)
      READ(31) (PSURF(I),I=1,NUMNP)
      READ(31) (BDROCK(I),I=1,NUMNP)
      READ(31) (KX(I,1),KX(I,2),KX(I,3),KX(I,4),I=1,NUMEL)
      READ(31) (IBFLUX(I,1),IBFLUX(I,2),BFLUX(I),I=1,NUMGBC)
C
C ... READ INPUT DIFFER, THINGS THAT DIFFER FROM ONE GRID TO THE NEXT
      READ(32) HED
      READ(32) (ADOT(I),I=1,NUMNP)
      READ(32) (FRACT(I),I=1,NUMNP)
      READ(32) (FLOWA(I),I=1,NUMNP)
      READ(32) (SLDGB(I),I=1,NUMNP)
      READ(32) (TEMP(I),I=1,NUMNP)
      READ(32) (ITYPE(I),I=1,NUMNP)
      READ(32) (AFUDGE(I),I=1,NUMNP)
      READ(32) (GEOFLUX(I),I=1,NUMNP)
      IF(GEOFLUX(1).GT.0. .AND. GEOFLUX(1).LE.10.) THEN
        PRINT *,'THIS IS OLD DATA SET, CHANGE GEOTHERMAL'
c        PAUSE
      ENDIF
      READ(32) (CALV(I),I=1,NUMNP)
      DO I=1,NUMNP
        TSORIG(I)=TEMP(I)
        IF(TSORIG(I).EQ.-999.) TSORIG(I)=0.0d0
        IF(AFUDGE(I).EQ.0.) AFUDGE(I)=1.D0
        IF(FLOWA(I).LE.0.) FLOWA(I)=1D-6
        IF(SLDGB(I).LE.0.) SLDGB(I)=1D-9
        IF(CALV(I).LT.0.) CALV(I)=0.01d0
c ************EISMINT GEOTHERMAL FLUX VALUES ********************
c ... 42 milliWatts/m**2 sliding experiment
c        IF(GEOFLUX(I).EQ.0.) GEOFLUX(I)=3.17D5
c ... 50 milliWatts/m**2 greenland experiment 
c        IF(GEOFLUX(I).EQ.0.) GEOFLUX(I)=3.77D5
c ... 54.6 milliWatts/m**2 ANTARCTIC VALUE (HUYBRECHTS) 
        IF(GEOFLUX(I).EQ.0.) GEOFLUX(I)=4.121D5
c ... 25 milliWatts/m**2 greenland experiment 
         IF(GEOFLUX(I).EQ.0.) GEOFLUX(I)=3.77D5/2
c        IF(GEOFLUX(I).LE.10.) THEN
c          PRINT *,'THIS IS OLD DATA SET, CHANGE GEOTHERMAL'
c          GEOFLUX(I)=4.121D5
c        ENDIF
      ENDDO
C
C ... READ INPUT TIME, THINGS THAT CHANGE WITH TIME
      READ(33) HED
      READ(33) (HTICE(I),I=1,NUMNP)
      READ(33) (ACC(I),I=1,NUMNP)
      READ(33) (ABLAT(I),I=1,NUMNP)
      READ(33) (CONST(I),I=1,NUMEL)
      READ(33) (ACON(I),I=1,NUMEL)
      DO I=1,NUMEL
        IF(ACON(I).LT.1E-30) ACON(I)=0.D0
      ENDDO
      IF(NUMNP.GT.MXX) THEN
        PRINT *,'NUMNP=',NUMNP,' MXX=',MXX,' INCREASE MXX'
        STOP
      ENDIF
      DO N=1,NUMNP
        ADOTB(N)=ADOT(N)
        IDT(N)=0
c ... turn off for madeleine/forget Mass Balance
        IF(.false. .and. ADOT(N).LT.-100.) THEN
         IDT(N)=-INT(100+ADOT(N))
         ADOT(N)=AFUNCT(TIME,IDT(N),AMASS,HTICE(N),BDROCK(N),PSURF(N),
     &                  ZERO,X(N),Y(N),TEMP(N))
        ENDIF
c       IF(PSURF(N).LT.BDROCK(N)-SEALEV) PSURF(N)=BDROCK(N)
        THICK(N)=HTICE(N)-BDROCK(N)
      ENDDO
C ... UNLOAD THE BED ...
      IF(BTOGG.ne.0) THEN
        WRITE(7,123) ' TIME BEFORE 1ST PLATE ',ETIME(TB),DTIME(TB)
        print *,'unloading bed with present surface'
C----------------------------------------------------------------!
        LDEP=0                                                   1
        DO I=1,NUMNP                                             !
          THICKL(I)=PSURF(I)-BDROCK(I)                           !
          if(bdrock(i).lt.SEALEV) then
             FLOT=(1.D0-RATDEN)*(BDROCK(I)-SEALEV)                     !
             IF(PSURF(I).LT.FLOT) THEN                           !
               THICKL(I)=0.D0                                    !
               THICKL(I)=FLOT-BDROCK(I)                          !
             ENDIF                                               !
          endif
          IF(BDROCK(I).le.-9999.) THICKL(I)=0.D0                                       !
          IF(THICKL(I) .GT. 0.D0) LDEP=1                         !
        ENDDO                                                    !
        IF(LDEP.EQ.1) THEN                                       !
          if(BTOGG.eq.1) then
            CALL SPLATE(-1,NUMNP,THICKL,UNDEPB,DEPB,RHOI,RHOR,RHOW,
     &                  0.d0,TIME,wrate,www,wwworig)
          else
            CALL OPLATE(-1,MXX,3*MXX,NUMNP,NUMEL,X,Y,KX,THICKL,   !
     &                 0.d0,WWW,WRATE,WMIN,0.d0,WWWORIG)          !
          endif
          IF(.FALSE.) THEN                                       !
            CALL POUTSTF(NUMNP*3,WWW,WRATE(1,1),NUMNP,THICKL,    !
     &                   NUMCOL,NUMLEV,0.d0,WWWORIG)             !
          ENDIF                                                  !
c         sum=0.                                                 !
          DO I=1,NUMNP      
            if(.false.) then
              UNDEPB(I)=BDROCK(I)-WWW((I-1)*3+1)                   !
            else
              UNDEPB(I)=BDROCK(I)-WWWORIG(I)                       !
            endif
            if(.true.) then
              DEPB(I)=UNDEPB(I)+WWW((I-1)*3+1)                     !
            else
              DEPB(I)=UNDEPB(I)+WWWORIG(I)                         !
            endif

c            DEPB(I)=BDROCK(I)
            UNDEPBI=RHOI*(PSURF(I)-BDROCK(I))/RHOR+BDROCK(I)
c           print *,real(WWW((I-1)*3+1)),REAL(UNDEPB(I)),REAL(depb(i))                                     !
c           PRINT *,I,REAL(UNDEPB(I)),REAL(UNDEPBI),
c    &            REAL(bdrock(I)-(UNDEPBI))/real(WWW((I-1)*3+1))
c           sum=sum+(bdrock(I)-(UNDEPBI))/(WWW((I-1)*3+1))
          ENDDO                  
c         print *,sum/numnp
        ENDIF                                                    !
        WRITE(7,123) ' TIME AFTER  1ST PLATE ',ETIME(TB),DTIME(TB)
C----------------------------------------------------------------!
      ELSEIF(BTOGG.eq.0) then
C ... OLD WAY ...
        DO I=1,NUMNP
          IF(BDROCK(I).GT.-9999.) THEN
            UNDEPB(I)=RHOI*(PSURF(I)-BDROCK(I))/RHOR+BDROCK(I)
            UNDEPB(I)=BDROCK(I)
          ENDIF
          DEPB(I)=UNDEPB(I)
        ENDDO
      else
        print *,' problems with BTOGG',BTOGG
        pause
      ENDIF
      NELMAX=-10000
      NCLMAX=-10000
      DO N=1,NUMEL
        NTYPE(N)=1
        IF(KX(N,4).EQ.0) NTYPE(N)=2
        NODEN=4
        NNODE(N)=4
        IF(NTYPE(N).EQ.2) NODEN=3
        IF(NTYPE(N).EQ.2) NNODE(N)=3
        IMAX=-10000
        IMIN=10000
        DO I=1,NODEN
          IF(KX(N,I).GT.IMAX) IMAX=KX(N,I)
          IF(KX(N,I).LT.IMIN) IMIN=KX(N,I)
        ENDDO
C **** FOR MAP1 WITH ASYMSL ****
        NCOL1=2*(IMAX-IMIN)+1
        IF(NCOL1.GT.NCLMAX) THEN
          NCLMAX=NCOL1
          NELMAX=N
        ENDIF
      ENDDO
      WRITE(*,*) 'MAX NCOL=',NCLMAX,' IN ELEMENT',NELMAX
c     IF(NCLMAX.GT.NCOL) THEN
c       WRITE(*,*) 'CURRENTLY USING NCOL=',NCOL,' INCREASE'
c       STOP
c     ENDIF
C
      SUMC = 0.0d0
      DO N = 1,NUMEL
        QHOLD(N) = 0.0d0
c       IF (CONST(N).EQ.0D0 .AND. DT.EQ.0D0) CONST(N)=1.92d7
        IF (CONST(N).EQ.0D0) CONST(N)=1.92d7
c       CONST(N)=1.92d7
        SUMC = SUMC + CONST(N)
        IF(NTYPE(N).EQ.1) THEN
          NODEN=4
          NINT=9
        ELSE
          NODEN=3
          NINT=4
        ENDIF
        AAADOT=0.d0
        AAFRCT=0.d0
        AABDRK=0.d0
        APSURF=0.d0
        AAFLOW=0.d0
        AASLDG=0.d0
        ACALV=0.d0
        DO I=1,4
          LM(I)=KX(N,I)
          AAADOT=AAADOT+ADOT(LM(I))
          AAFRCT=AAFRCT+FRACT(LM(I))
          AABDRK=AABDRK+BDROCK(LM(I))
          APSURF=APSURF+PSURF(LM(I))
          AAFLOW=AAFLOW+FLOWA(LM(I))
          AASLDG=AASLDG+SLDGB(LM(I))
          ACALV=ACALV+CALV(LM(I))
        ENDDO
        DENOM=1.d0/DBLE(NODEN)
        AAADOT=AAADOT*DENOM
        AADOT(N)=AAADOT
        AAFRCT=AAFRCT*DENOM
        AFRACT(N)=AAFRCT
        AABDRK=AABDRK*DENOM
        ABDRCK(N)=AABDRK
        APSURF=APSURF*DENOM
        PPSURF(N)=APSURF
        ACALV=ACALV*DENOM
        PCALV(N)=ACALV
        AAFLOW=AAFLOW*DENOM
        IF(ACON(N).EQ.0.) THEN
          AFLOWA(N)=AAFLOW
        ELSE
          AFLOWA(N)=ACON(N)
        ENDIF
        AASLDG=AASLDG*DENOM
        ASLDGB(N)=AASLDG
      ENDDO
C
      SUMC = SUMC/NUMEL
      WRITE(*,*) SUMC
      if(.false.) then
        CALL WRITPROF(NUMNP,X,Y,HTICE,BDROCK,ADOT)
        pause
      endif

      RETURN
999   CONTINUE
      PRINT *,'NO DATA SET'
      STOP
      END
c--------------------------------------------
      SUBROUTINE READCOORD(NUMNP)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      common /coord/ CLAT(MAXNUM),CLONG(MAXNUM)
      READ(53,end=100) HED
      print *,'reading coords'
      READ(53) (CLAT(I),I=1,NUMNP)
      READ(53) (CLONG(I),I=1,NUMNP)
      return
100   continue
      print *,'no coord file available'
      END
