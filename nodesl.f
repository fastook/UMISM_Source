      SUBROUTINE NODESL(MXX,NUMNP,NUMEL,KX,SLOPE,SLOPN)
C ... CALCULATES NODE SLOPES FROM ELEMENT SLOPES 
C ... FOR ACCUMULATION PARAMETERIZATIONS
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER(MMXX=MAXNUM)
      DIMENSION KX(MXX,4),SLOPE(4,MXX),SLOPN(4,MXX),ICOUNT(MMXX)
      if(.false.) then ! average of element slopes
        DO I=1,NUMNP
          DO L=1,4
            SLOPN(L,I)=0.d0
          ENDDO
          ICOUNT(I)=0
        ENDDO
        DO I=1,NUMEL
          DO J=1,4
            DO L=1,4
              SLOPN(L,KX(I,J))=SLOPN(L,KX(I,J))+SLOPE(L,I)
            ENDDO
            ICOUNT(KX(I,J))=ICOUNT(KX(I,J))+1
          ENDDO
        ENDDO
        DO I=1,NUMNP
          IF(ICOUNT(I).GT.0) THEN
            DO L=1,4
              SLOPN(L,I)=SLOPN(L,I)/DBLE(ICOUNT(I))
            ENDDO
c           SLOPN(1,I)=SQRT(SLOPN(2,I)**2+SLOPN(3,I)**2)
          ENDIF
        ENDDO
      else                     ! minimum of element slopes
        DO I=1,NUMNP
          DO L=1,4
            SLOPN(L,I)=1d30
          ENDDO
          ICOUNT(I)=0
        ENDDO
        DO I=1,NUMEL
          DO J=1,4
            DO L=1,4
              SLOPN(L,KX(I,J))=min(SLOPN(L,KX(I,J)),SLOPE(L,I))
            ENDDO
            ICOUNT(KX(I,J))=ICOUNT(KX(I,J))+1
          ENDDO
        ENDDO
        DO I=1,NUMNP
          IF(ICOUNT(I).eq.0) THEN
            print *,'problems in node',i,' with nodesl'
          ENDIF
        ENDDO
      endif
      RETURN
      END
