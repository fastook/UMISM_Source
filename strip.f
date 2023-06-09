      PARAMETER(NMAX=29999)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER HED*80
      DIMENSION KODE(NMAX),X(NMAX),Y(NMAX),HTICE(NMAX),ADOT(NMAX),
     &FRACT(NMAX),PSURF(NMAX),bed(NMAX),FLOWA(NMAX),SLDGB(NMAX),
     &KX(NMAX,4),CONST(NMAX),IBFLUX(NMAX,2),BFLUX(NMAX),temp(nmax),
     &       ifd(nmax),ifp(nmax)
      READ(1,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      PRINT *,HED
1000  FORMAT (A80,/,7I6,F8.1)
      IF(NUMNP.GT.NMAX) THEN
        PRINT *,'NUMNP=',NUMNP,' NMAX=',NMAX,' INCREASE NMAX'
        STOP
      ENDIF
      print *,numnp,' nodes,',numel,' elements'
      DO N=1,NUMNP
        READ(1,1001) NUM,KODE(N),X(N),Y(N),HTICE(N),
     &               ADOT(N),FRACT(N),PSURF(N),bed(N),
     &               FLOWA(N),SLDGB(N),temp(n)
1001  FORMAT(I6,I4,1P2E12.5,0PF10.2,F7.2,F9.2,F10.3,F10.1,F10.5,2F10.5,
     &       I5,F10.3)
      enddo
      DO N=1,NUMEL
        READ(1,1002) NUM,KX(NUM,1),KX(NUM,2),KX(NUM,3),KX(NUM,4),
     &               CONST(N)
 1002 FORMAT(5I6,1PE17.10)
      enddo
      IF(NUMGBC.GT.0) THEN
        DO N=1,NUMGBC
          READ(1,1007) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
 1007 FORMAT(2I6,E13.6)
        enddo
      ENDIF
      IFOUND=0
      icount=0
      do i=1,numnp
         if(kode(i).eq.1) then
c        if(adot(i).ne.-115.) then
c        if(bed(i).le.-1200) then
          icount=icount+1
          ifp(icount)=i
        endif
      enddo
      print *,'removing',icount
      istep=icount/20
      DO IC=1,ICOUNT
        if(mod(ic,istep).eq.0) print *,'first pass',ic
        DO I=1,NUMEL
          IF(KX(I,4).EQ.0) THEN
            NN=3
          ELSE
            NN=4
          ENDIF
          DO N=1,NN
            IF(KX(I,N).EQ.IFP(IC)) THEN
C             FOUND ELEMENT CONTAINING NODE, DELETE
              IFOUND=IFOUND+1
              IFD(IFOUND)=I
              GOTO 3810
            ENDIF
          enddo
3810      CONTINUE
        enddo
      enddo
      CALL SORTIX(IFOUND,IFD)
      CALL ELMDUP(IFOUND,IFD)
      DO I=1,IFOUND
        DO J=IFD(I),NUMEL-1
          DO K=1,4
            KX(J,K)=KX(J+1,K)
            CONST(J)=CONST(J+1)
          enddo
        enddo
        NUMEL=NUMEL-1
      enddo
      CALL SORTIX(ICOUNT,IFP)
      DO J=1,ICOUNT
        if(mod(j,istep).eq.0) print *,'second pass',j
        DO N=1,NUMEL
          DO K=1,4
            IF(KX(N,K).GE.IFP(J)) KX(N,K)=KX(N,K)-1
          enddo
        enddo
      enddo
      DO I=1,ICOUNT
        if(mod(i,istep).eq.0) print *,'third pass',i
        DO J=IFP(I),NUMNP-1
          x(J)=x(J+1)
          y(J)=y(J+1)
          HTICE(J)=HTICE(J+1)
          ADOT(J)=ADOT(J+1)
          BED(J)=BED(J+1)
          FRACT(J)=FRACT(J+1)
          PSURF(J)=PSURF(J+1)
          FLOWA(J)=FLOWA(J+1)
          SLDGB(J)=SLDGB(J+1)
          KODE(J)=KODE(J+1)
          ADOT(J)=ADOT(J+1)
          TEMP(J)=TEMP(J+1)
        enddo
        NUMNP=NUMNP-1
      enddo
      write(11,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      PRINT *,HED
      print *,numnp,' nodes,',numel,' elements'
      DO N=1,NUMNP
        write(11,1001) n,KODE(N),X(N),Y(N),HTICE(N),
     &               ADOT(N),FRACT(N),PSURF(N),bed(N),
     &               FLOWA(N),SLDGB(N),temp(n)
      enddo
      DO N=1,NUMEL
        write(11,1002) n,KX(n,1),KX(n,2),KX(n,3),KX(n,4),
     &               CONST(N)
      enddo
      IF(NUMGBC.GT.0) THEN
        DO N=1,NUMGBC
          write(11,1007) IBFLUX(N,1),IBFLUX(N,2),BFLUX(N)
        enddo
      ENDIF
      END
c====================================================================
      SUBROUTINE SORTIX(N,IX)
      DIMENSION IX(N)
      LOGICAL ANYEXC
      ANYEXC=.TRUE.
C WHILE ANY EXCHANGES
30    IF(ANYEXC) THEN
        ANYEXC=.FALSE.
        iswap=0
        DO 40 I=1,N-1
          IF(IX(I).LT.IX(I+1)) THEN
            iswap=iswap+1
            CALL SWAP(IX(I),IX(I+1))
            ANYEXC=.TRUE.
          ENDIF
40      CONTINUE
        GOTO 30
      ENDIF
      RETURN
      END
c======================================================
      SUBROUTINE ELMDUP(NN,IX)
      DIMENSION IX(NN)
      NNOUT=NN
      DO 100 N=2,NN
80      IF(IX(N-1).EQ.IX(N) .AND. (IX(N).NE.0)) THEN
          NNOUT=NNOUT-1
          DO 90 J=N,NN-1
            IX(J)=IX(J+1)
90        CONTINUE
          IX(NN)=0
          GOTO 80
        ENDIF
100   CONTINUE
      NN=NNOUT
      RETURN
      END
c==================================================
      SUBROUTINE SWAP(IX1,IX2)
      ITEMP=IX1
      IX1=IX2
      IX2=ITEMP
      RETURN
      END
