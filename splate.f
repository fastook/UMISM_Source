C===================================================
      SUBROUTINE UNLOAD(NUMNP,BDROCK,UNDEPB,PSURF,RHOI,RHOR)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      DIMENSION BDROCK(NUMNP),UNDEPB(NUMNP),PSURF(NUMNP)
      RETURN
      IF(BTOGG.ne.0) RETURN
      DO N=1,NUMNP
        IF(BDROCK(N).GT.-9999.) THEN
          UNDEPB(N)=RHOI*(PSURF(N)-BDROCK(N))/RHOR+BDROCK(N)
        ENDIF
      ENDDO
      END
C===================================================
      SUBROUTINE SPLATE(ITIME,NUMNP,THICK,UNDEPB,DEPB,RHOI,RHOR,RHOW,
     &                  DELT,TIME,WRATE,WWW,WWWORIG)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM,N3=NMAX*3)
      PARAMETER(RELAX=3000.D0,RK=1.D0/RELAX)
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      DIMENSION UNDEPB(NMAX),DEPB(NMAX),THICK(NMAX)
      DIMENSION WWW(N3),WRATE(N3,2),WWWORIG(NMAX)
      CHARACTER*80 HED
      COMMON /SELASTIC/ WWWS(N3),WSRATE(N3,2),WWWSORIG(NMAX)
      PARAMETER(NPAGE=39)
      CHARACTER*80 LIST(1000)
      COMMON /IOLIST/ LIST
      logical file40
      DATA IPASS /0/, IREAD /0/
      SAVE IPASS,IREAD,RATIO,SIXTH,ROCKICE
      NNODE=3*NUMNP
      IF(IPASS.EQ.0) THEN
C ....  STUFF YOU ONLY DO THE FIRST PASS ..............
        SIXTH=1D0/6D0
        RATIO=RHOI/RHOR
        ROCKICE=RHOR/RHOI
        PRINT *,'READING BEDROCK DEPRESSION FILE'
        inquire(file='fort.40',exist=file40)
        iok=-1
        if(file40) then
          READ(40,*,IOSTAT=IOK) NNN
        endif
        IF(IOK.EQ.0) THEN
          IF(NNN.NE.NNODE) THEN
            PRINT *,'PROBLEMS:'
            PRINT *,'INCOMPATIBLE WITH CURRENT NNODE=',NNODE,NNN
            IOK=1
          ENDIF
          DO I=1,NNODE
            READ(40,*) II,WWW(I),WJUNK
            IF(I.NE.II) THEN
              PRINT *,'PROBLEMS:READING WWW'
              PRINT *,'INCOMPATIBLE WITH CURRENT (WWW)'
              IOK=1
            ENDIF
          ENDDO
          DO I=1,NNODE
            READ(40,*) II,WRATE(I,1),WJUNK
            IF(I.NE.II) THEN
              PRINT *,'PROBLEMS:READING WRATES(1)'
              PRINT *,'INCOMPATIBLE WITH CURRENT (WWW)'
              IOK=1
            ENDIF
          ENDDO
          DO I=1,NNODE
            READ(40,*) II,WRATE(I,2),WJUNK
            IF(I.NE.II) THEN
              PRINT *,'PROBLEMS:READING WRATES(2)'
              PRINT *,'INCOMPATIBLE WITH CURRENT (WWW)'
              IOK=1
            ENDIF
          ENDDO
          DO I=1,NUMNP
            READ(40,*) II,WWWORIG(I),WJUNK
            WWWSORIG(I)=WWWORIG(i)
            IF(I.NE.II) THEN
              PRINT *,'PROBLEMS:READING WWWORIGS'
              PRINT *,'INCOMPATIBLE WITH CURRENT (WWWORIG)'
              IOK=1
            ENDIF
          ENDDO
          PRINT *,' BEDROCK DEPRESSION FILE FOUND'
          PRINT *,'    AND READ SUCCESSFULLY '
          IF(ITIME.LT.0) THEN
            PRINT *,'ABANDONING UNLOADING'
            REWIND 40
            RETURN
          ENDIF
          IREAD=1
        ENDIF
        IF(IOK.NE.0) THEN
          PRINT *,' NONE FOUND, SET TO ZERO ... AND UNLOADED ...'
C$DOACROSS LOCAL(I)
          DO I=1,NNODE
            WWW(I)=0.D0
            WRATE(I,1)=0.D0
            WRATE(I,2)=0.D0
          ENDDO
          IF(DELT.LE.0.0) THEN
            II=1
            DO I=1,NUMNP
              WWW(II)=-RATIO*THICK(I)
              WRATE(II,2)=WWW(II)
              II=II+3
            ENDDO
          ENDIF
          WRITE(HED,*) ' TIME=',NINT(TIME-DELT)
          WRITE(88) HED
          WRITE(88) (WWW(I),I=1,NNODE,3)
          WRITE(88) (THICK(I),I=1,NUMNP)
          WRITE(88) (1000*WRATE(I,1),I=1,NNODE,3)
        ENDIF
      ENDIF
C ... END OF STUFF DONE ONLY ONCE ...
      II=1
      DO I=1,NUMNP
        RHS1=-RK*(WWW(II)           +RATIO*THICK(I))*DELT
        RHS2=-RK*(WWW(II)+0.5D0*RHS1+RATIO*THICK(I))*DELT
        RHS3=-RK*(WWW(II)+0.5D0*RHS2+RATIO*THICK(I))*DELT
        RHS4=-RK*(WWW(II)+RHS3      +RATIO*THICK(I))*DELT
        SUM=(RHS1+RHS2+RHS2+RHS3+RHS3+RHS4)*SIXTH
        WWW(II)=WWW(II)+SUM
        WRATE(II,1)=SUM/DELT
c        DEPB(I)=UNDEPB(I)+WWW(II)
C        IF(THICK(I).NE.0.) THEN
C          PRINT 100,I,THICK(I),WWW(II),WRATE(II,1),DEPB(I),UNDEPB(I)
C        ENDIF
        II=II+3
      ENDDO
100   FORMAT(1X,I5,5G13.6)
      IF(DELT.GT.0. .and. .false.) THEN
        DO I=1,NNODE
          WRATE(I,1)=(WWW(I)-WRATE(I,2))/DELT
        ENDDO
      ENDIF
      DO I=1,NNODE
        WRATE(I,2)=WWW(I)
      ENDDO
      WMIN=1D30
      WMAX=-1D30
      WTOT=0.D0
      DO I=1,NNODE,3
        WMIN=MIN(WMIN,WWW(I))
        WMAX=MAX(WMAX,WWW(I))
        WTOT=WTOT+WWW(I)
C        WRITE(*,*) (WWW(I+J),J=0,2)
        WWWS(I)=WWW(I)
        WSRATE(I,1)=WRATE(I,1)
        WSRATE(I,2)=WRATE(I,2)
      ENDDO
      THTOT=0.D0
      DO I=1,NUMNP
        THTOT=THTOT+THICK(I)
      ENDDO
      IF(IREAD.EQ.0 .AND. IPASS.EQ.0) THEN
        II=1
        DO I=1,NUMNP
          WWWORIG(I)=WWW(II)
          WWWSORIG(I)=WWW(II)
          II=II+3
        ENDDO
      ENDIF
      IF(IOTOGG) THEN
        WRITE(LIST(IPAGE+1),*) 
     &       '*******************************************'
        WRITE(LIST(IPAGE+2),1001) TIME,REAL(THTOT),
     &              REAL(-THTOT/WTOT/ROCKICE),
     &              REAL(WMIN),REAL(WMAX),REAL(WSAVE-WMIN)
        WRITE(LIST(IPAGE+3),*) 
     &       '*******************************************'
        IPAGE=IPAGE+3
      ENDIF
      WRITE(92,*) TIME,WMIN
      WSAVE=WMIN
      IPASS=1
1001  FORMAT(1X,6(1X,1PG12.5))
      END
        

C---------------------------------------------
      SUBROUTINE WRITEDEPS(NUMNP,NNODE)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM,N3=3*NMAX)
      CHARACTER*1 CHAR
      COMMON /SELASTIC/ WWW(N3),WRATE(N3,2),WWWORIG(NMAX)
      PRINT *,'IN WRITEDEP',NNODE
      PRINT *,'   TO WRITE OUT BACKUP OF BEDROCK DEPRESSION '
      PRINT *,'   INPUT Y'
      READ(*,1002) CHAR
      WRITE(99,1002) CHAR
1002  FORMAT(A1)
      IF(CHAR.EQ.'Y' .OR. CHAR.EQ.'y') THEN
        REWIND 45
        WRITE(45,*) NNODE
        DO I=1,NNODE
          WRITE(45,*) I,WWW(I),0.D0
        ENDDO
        DO I=1,NNODE
          WRITE(45,*) I,WRATE(I,1),0.D0
        ENDDO
        DO I=1,NNODE
          WRITE(45,*) I,WRATE(I,2),0.D0
        ENDDO
        DO I=1,NUMNP
          WRITE(45,*) I,WWWORIG(I),0.D0
        ENDDO
      ENDIF
      END
