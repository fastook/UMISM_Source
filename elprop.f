      SUBROUTINE ELPROP(MXX, NUMEL, NTYPE, KX, ADOT, 
     &           AADOT, FRACT,
     &           AFRACT, BDROCK, ABDRCK, PSURF, PPSURF, FLOWA, ACON,
     &           ICON, AFLOWA, SLDGB, ASLDGB,CALV,PCALV)
      IMPLICIT REAL*8(A-H,O-Z)
C ... LOAD NODAL PROPERTIES INTO ELEMENT PROPERTIES
      DIMENSION NTYPE(MXX),ADOT(MXX),AADOT(MXX),FRACT(MXX),AFRACT(MXX),
     &          BDROCK(MXX),ABDRCK(MXX),PSURF(MXX),PPSURF(MXX),
     &          FLOWA(MXX),AFLOWA(MXX),SLDGB(MXX),ASLDGB(MXX),
     &          LM(4),KX(MXX,4),ACON(MXX),CALV(MXX),PCALV(MXX)
C      REWIND 90
c      PRINT *,'IN ELPROP'
      DO N = 1,NUMEL
        IF(NTYPE(N).EQ.1) THEN
          NODEN=4
        ELSE
          NODEN=3
        ENDIF
        AAADOT=0.d0
        AAFRCT=0.d0
        AABDRK=0.d0
        APSURF=0.d0
        AAFLOW=0.d0
        AASLDG=0.d0
        ADNSTY=0.d0
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
c        WRITE(90,*) N,AFRACT(N)
        AABDRK=AABDRK*DENOM
        ABDRCK(N)=AABDRK
        ADNSTY=ADNSTY*DENOM
        PPSURF(N)=APSURF
        ACALV=ACALV*DENOM
c ..... calving is AVERAGE of 4 nodal values...
c       PCALV(N)=ACALV
c ..... calving is MINIMUM of 4 nodal values...
        PCALV(N)=MIN(CALV(LM(1)),CALV(LM(2)),CALV(LM(3)),CALV(LM(4)))
c ..... calving is MAXIMUM of 4 nodal values...
c       PCALV(N)=MAX(CALV(LM(1)),CALV(LM(2)),CALV(LM(3)),CALV(LM(4)))
        AAFLOW=AAFLOW*DENOM
        IF(ACON(N).EQ.0. .OR. ICON.EQ.0) THEN
          AFLOWA(N)=AAFLOW
        ELSE
          AFLOWA(N)=ACON(N)
        ENDIF
        AASLDG=AASLDG*DENOM
        ASLDGB(N)=AASLDG
      ENDDO
      RETURN
      END
c==========================================================
      SUBROUTINE CALVING(MXX, NUMEL, NTYPE, KX, HTICE, BED, ADOT, 
     &           AADOT,CFACTOR,ADC,CALV,PCALV,TNSL)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NSMAX=MAXTIME)
C ... LOAD NODAL PROPERTIES INTO ELEMENT PROPERTIES
      DIMENSION NTYPE(MXX),ADOT(MXX),AADOT(MXX),
     &          HTICE(MXX),BED(MXX),ADC(MXX),
     &          LM(4),KX(MXX,4),CALV(MXX),PCALV(MXX)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      COMMON /SEALV/ SEALEV,RATDEN,STIME(NSMAX),SLIST(NSMAX),NSEAL
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      DATA RHOI /0.917D0/, RHOW /1.092D0/
      IF(CFACTOR.EQ.0.) RETURN
      NCALV=0
      acalv=0.d0
      if(.true.) then
        rmult=1.d0
      else
c ... make calving depend on tnsl X 0 at -10, X 1 at 0
        ttemp=4.5d0
        ttemp=6.0d0
        rmult=(ttemp+tnsl)/ttemp
        rmult=min(rmult,1.d0)
        rmult=max(rmult,0.d0)
      endif
      if(.false.) then
c ..... old way ......
        DO N = 1,NUMEL
          IF(CFACTOR.LT.0) THEN
            CFACTLOC=PCALV(N)
          ELSE
            CFACTLOC=CFACTOR
            CFACTLOC=CFACTOR*rmult
          ENDIF
          IF(NTYPE(N).EQ.1) THEN
            NODEN=4
          ELSE
            NODEN=3
          ENDIF
          AABED=0.d0
          AAHTICE=0.d0
          DO I=1,4
            LM(I)=KX(N,I)
            AABED=AABED+BED(LM(I))
            AAHTICE=AAHTICE+HTICE(LM(I))
          ENDDO
          DENOM=1.d0/DBLE(NODEN)
          AABED=AABED*DENOM
          AAHTICE=AAHTICE*DENOM
C ....... THIS IS FOR CALVING ...
          IF(AABED.LT.SEALEV) THEN
            FLOT=(1.d0-RATDEN)*(AABED-SEALEV)
            THIK=AAHTICE-FLOT
            IF(THIK.LT.0.) THEN
              NCALV=NCALV+1
              DC=CFACTLOC*THIK
              AADOT(N)=MIN(DC,AADOT(N))
              DO I=1,4
                ADC(LM(I))=DC
                ADOT(LM(I))=AADOT(N)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      else
c ..... new way ......
        DENOM=0.25d0
        DO N = 1,NUMEL
          IF(CFACTOR.LT.0) THEN
            CFACTLOC=PCALV(N)
c            CFACTLOC=PCALV(N)*rmult
          ELSE
            CFACTLOC=CFACTOR
c            CFACTLOC=CFACTOR*rmult
          ENDIF
          AABED=0.d0
          AAHTICE=0.d0
          icalv=0
          DO I=1,4
            LM(I)=KX(N,I)
            AABED=AABED+BED(LM(I))
            AAHTICE=AAHTICE+HTICE(LM(I))
            if(BED(LM(I)).lt.SEALEV) then
              FLOT=(1.d0-RATDEN)*(BED(LM(I))-SEALEV)
              if(HTICE(LM(I)).lt.FLOT) icalv=icalv+1
            endif
          ENDDO

c ... old way
          if(icalv.gt.0 .and. icalv.lt.5 .and. .false.) then
            AABED=AABED*DENOM
            AAHTICE=AAHTICE*DENOM
C ......... THIS IS FOR CALVING ...
            IF(AABED.LT.SEALEV) THEN
              FLOT=(1.d0-RATDEN)*(AABED-SEALEV)
              THIK=AAHTICE-FLOT
              IF(THIK.LT.0.) THEN
                DC=CFACTLOC*THIK
                if(icalv.lt.4) then
                  NCALV=NCALV+1
                  acalv=acalv+dc
                endif
c              AADOT(N)=MIN(DC,AADOT(N))
              AADOT(N)=DC+AADOT(N)
                DO I=1,4
                  ADC(LM(I))=DC
                  ADOT(LM(I))=AADOT(N)
                ENDDO
              ENDIF
            ENDIF
          endif

c ... new way 28 feb 2003 ...
          if(icalv.gt.0 .and. icalv.lt.5) then
            AABED=AABED*DENOM
            AAHTICE=AAHTICE*DENOM
C ......... THIS IS FOR CALVING ...
            IF(AABED.LT.SEALEV) THEN
              FLOT=(1.d0-RATDEN)*(AABED-SEALEV)
c ... following is weertman spreading of unconstrained ice shelf (Z**3)
                DC=(1.1e-7*((1.d0-RATDEN)*(FLOT-AABED))**3)*CFACTLOC 
c ... fudge to make CFACTLOC roughly 1.0 ...
                DC=DC*10    
                dcrmult=dc*rmult**2
                if(icalv.lt.4) then
                  NCALV=NCALV+1
                  acalv=acalv+dcrmult
                endif
                DO I=1,4
                  if(BED(LM(I)).lt.SEALEV) then
                    FLOT=(1.d0-RATDEN)*(BED(LM(I))-SEALEV)
                    if(HTICE(LM(I)).lt.FLOT) then
                      ADC(LM(I))=dcrmult
                      ADOT(LM(I))=max(-100.d0,ADOT(LM(I))+dcrmult)
c                      ADOT(LM(I))=dcrmult
                    endif
                  endif
                ENDDO
            ENDIF
          endif


        ENDDO
      endif
      IF(IOTOGG) then
        write(list(ipage+1),*) NCALV,'= NUMBER CALVING...',real(acalv),
     &                         real(rmult)
        ipage=ipage+1
      endif
      END
