      SUBROUTINE PLATE(ITIME,NMAX,N3,NUMNP,NUMEL,X,Y,KX,THICK,
     &                 DELT,WWW,WRATE,WMIN,TIME)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      PARAMETER(MMAX=12999,NZ=27, NZ1=NZ+1, M3=3*MMAX)
c
      parameter(nit=m3*nz,nmax1=m3+1,nmax3=m3*3)
      parameter(itmax=m3/10,ncg=4*itmax,nw=4*m3+ncg)
c ....arrays for ITPACK sparse storage...................
      dimension GK(nit),ja(nit),ia(nmax1)
      dimension iwksp(nmax3),wksp(nw),iwork(nit)
      dimension iparm(12),rparm(12)
c
      DIMENSION THICK(NMAX),WWW(N3),X(NMAX),Y(NMAX),KX(NMAX,4)
      DIMENSION WRATE(N3,2)
      DIMENSION KKX(MMAX,12),LTEMP(MMAX,3)
      DIMENSION XI(2,4),W(4)
      DIMENSION EK(12,12),EC(12,12),EF(12)
      DIMENSION GF(M3)
      CHARACTER*80 HED
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
C ... TIMER STUFF, FOR SGI ONLY ...!
      REAL*4 TB(2)
c     REAL*4 TB(2),ETIME,DTIME     !
c     EXTERNAL ETIME,DTIME         !
C .................................!
      SAVE IPASS,KKX,WSAVE,NNODE,XI,W
      DATA IPASS /0/
1000  FORMAT(1X,A,T25,1PG13.6,G13.6)
1001  FORMAT(1X,I5,1P5(1X,G13.6))
      IF(IPASS.EQ.0) THEN
C ....  STUFF YOU ONLY DO THE FIRST PASS ..............
        IF(NMAX.NE.MMAX .OR.N3.NE.M3) THEN
          PRINT *,'PROBLEMS WITH NMAX,MMAX:',NMAX,MMAX
          STOP
        ENDIF
        CALL PCONNECT(NMAX,NUMNP,NUMEL,KX,KKX,LTEMP)
        WSAVE=0.D0
        CALL PSETINT(XI,W)
        CALL PSETMAT
        NNODE=3*NUMNP
        PRINT *,'READING BEDROCK DEPRESSION FILE'
        READ(40,*,END=11) NNN
        IF(NNN.NE.NNODE) THEN
          PRINT *,'PROBLEMS:'
          PRINT *,'incompatible with current NNODE=',NNODE,NNN
          GOTO 11
        ENDIF
        DO I=1,NNODE
          READ(40,*,END=11) II,WWW(I)
          WRATE(I,2)=WWW(I)
          IF(I.NE.II) THEN
            PRINT *,'PROBLEMS:'
            PRINT *,'incompatible with current NMAX=',NMAX
            GOTO 11
          ENDIF
        ENDDO
        IPASS=1
      PRINT *,' BEDROCK DEPRESSION FILE FOUND AND READ SUCCESSFULLY '
      ENDIF
11    CONTINUE
      IF(IPASS.EQ.0) THEN
        PRINT *,' NONE FOUND, SET TO ZERO ... '
        DO I=1,NUMNP
          WWW(I)=0.D0
          WRATE(I,2)=WWW(I)
        ENDDO
        IPASS=1
        WRITE(HED,*) ' TIME=',NINT(TIME-DELT)
        WRITE(88) HED
        WRITE(88) (WWW(I),I=1,NNODE,3)
        WRITE(88) (THICK(I),I=1,NUMNP)
        WRITE(88) (1000*WRATE(I,1),I=1,NNODE,3)
C ......................................................
      ENDIF
C ....FORM STIFFNESS AND LOAD ..............
      WRITE(7,1000) ' TIME BEFORE PFORMKF ',ETIME(TB),DTIME(TB)
      call sbini(NNODE,nit,ia,ja,GK,iwork)
      CALL PFORMKF(NMAX,N3,NZ,NUMEL,NNODE,
     &            X,Y,THICK,KX,KKX,EK,EC,EF,XI,W,
     &            GK,GF,DELT,
     &            nit,ia,nmax1,ja,iwork)
C ....TIME DEPENDENT CASE AND VARYING LOAD .......
      CALL PTIMEDEP(N3,NZ,NNODE,GK,GC,GF,KA,KZ,WWW,DELT)
C ......................................................
C.....DUMP MATRIX FOR EXAMINATION ......................
      CALL PDUMPMAT(N3,NZ,NNODE,KZ,KA,GK,GF,WWW,.true.)
C.......................................................
C     IF(.TRUE.) THEN
C.......SOLVE EQUATIONS WITH GAUSS-SEIDEL ........................!
        WRITE(7,1000) ' TIME BEFORE PGAUSEID ',ETIME(TB),DTIME(TB)!
        CALL PGAUSEID(N3,NZ,NNODE,1.D-6,GK,KA,GF,WWW)                   !
        WRITE(7,1000) ' TIME AFTER PGAUSEID ',ETIME(TB),DTIME(TB) !
C ................................................................!
C     ELSE
C.......SOLVE EQUATIONS WITH CONJUGATE-GRADIENT ..................!
C       WRITE(7,1000) ' TIME BEFORE CONJUG ',ETIME(TB),DTIME(TB)  !
C       CALL CONJUG(N3,NZ,NNODE,1.D-6,GK,KA,GF,WWW)              !
C       WRITE(7,1000) ' TIME AFTER CONJUG ',ETIME(TB),DTIME(TB)   !
C ................................................................!
C     ENDIF
      IF(DELT.GT.0.) THEN
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
c        WRITE(*,*) (WWW(I+J),J=0,2)
      ENDDO
      THTOT=0.D0
      DO I=1,NUMNP
        THTOT=THTOT+THICK(I)
      ENDDO
      write(*,*) '*******************************************'
      WRITE(*,1001) ITIME,REAL(THTOT),REAL(-THTOT/WTOT/ROCKICE),
     &              REAL(WMIN),REAL(WMAX),REAL(WSAVE-WMIN)
      write(*,*) '*******************************************'
      WRITE(92,*) ITIME*DELT,WMIN
      WSAVE=WMIN
      END
C=============================================================
      SUBROUTINE PSETINT(XI,W)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XI(2,4),W(4)
C ...............................
      XI(1,1)=1.D0/SQRT(3.D0)
      XI(2,1)=XI(1,1)
      XI(1,2)=-XI(1,1)
      XI(2,2)=XI(1,1)
      XI(1,3)=XI(1,1)
      XI(2,3)=-XI(1,1)
      XI(1,4)=-XI(1,1)
      XI(2,4)=-XI(1,1)
      W(1)=1.D0
      W(2)=1.D0
      W(3)=1.D0
      W(4)=1.D0
C ...............................
      END
C=============================================================
      SUBROUTINE PELEM(nel,ETHICK,XY,N,EKB,EC,EF,NL,XI,W,DELT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XY(2,N),EKB(12,12),EKS(12,12),EKH(12,12)
      DIMENSION EF(12),EC(12,12)
      DIMENSION DPSIX(4),DPSIY(4)
      DIMENSION PSI(4),DPSI(4,2),XS(2)
      DIMENSION BB(3,12),BS(2,12),DB(3,3),DS(2,2)
      DIMENSION BTDB(3,12)
      DIMENSION XI(2,4),W(4)
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
C.....INITIALIZE ELEMENT ARRAYS
      DO I=1,12
        EF(I)=0.0D0
        DO J=1,12
          EKB(I,J)=0.0D0
          EKS(I,J)=0.0D0
          EKH(I,J)=0.0D0
          EC(I,J)=0.0D0
        ENDDO
      ENDDO
C ... FORM D-MATRICES ...............
      TTT3=TTT**3/12.D0
      DB(1,1)=TTT3*(2.D0*RMU+RLAMBDA)
      DB(1,2)=TTT3*RLAMBDA
      DB(1,3)=0
      DB(2,1)=TTT3*RLAMBDA
      DB(2,2)=TTT3*(2.D0*RMU+RLAMBDA)
      DB(2,3)=0
      DB(3,1)=0
      DB(3,2)=0
      DB(3,3)=TTT3*RMU
      DS(1,1)=TTT*RMU
      DS(1,2)=0
      DS(2,1)=0
      DS(2,2)=TTT*RMU
C.....BEGIN 2X2 INTEGRATION LOOP
      DO L=1,NL
        XS(1)=XI(1,L)
        XS(2)=XI(2,L)
        CALL PGENSHAPE(N,XS,XY,PSI,DPSI,DETJ,DPSIX,DPSIY)
        CALL PLOADB(PSI,DPSIX,DPSIY,BB,BS)
C
C.......ACCUMULATE INTEGRATION POINT VALUE OF INTEGRALS
        FAC=DETJ*W(L)
        XF=-ETHICK
C  .... RATIOS OF DENSITY OF ROCK AND ICE FOR OVERBURDEN LOAD
        XB=ROCKICE*RHOG
C ..... FORM CAPACITANCE MATRIX .....................
        DO I=1,4
          IP1=3*I-2
          IP2=3*I-1
          IP3=3*I
          EF(IP1)=EF(IP1)+XF*PSI(I)*FAC
          DO J=1,4
            JP1=3*J-2
            JP2=3*J-1
            JP3=3*J
            EKH(IP1,JP1)=EKH(IP1,JP1)+FAC*XB*PSI(I)*PSI(J)
            IF(DELT.GT.0) THEN
              TERM=FAC*CTIME*PSI(I)*PSI(J)
              EC(IP1,JP1)=EC(IP1,JP1)+TERM
              EC(IP1,JP2)=EC(IP1,JP2)+TERM
              EC(IP1,JP3)=EC(IP1,JP3)+TERM
              EC(IP2,JP1)=EC(IP2,JP1)+TERM
              EC(IP2,JP2)=EC(IP2,JP2)+TERM
              EC(IP2,JP3)=EC(IP2,JP3)+TERM
              EC(IP3,JP1)=EC(IP3,JP1)+TERM
              EC(IP3,JP2)=EC(IP3,JP2)+TERM
              EC(IP3,JP3)=EC(IP3,JP3)+TERM
            ENDIF
          ENDDO
        ENDDO
C ..... FORM KB STIFFNESS MATRIX ..............
C ..... OK, NOW HERE GOES THE MATRIX MULTIPLICATION THAT 
C ..... GENERATES THE BT D B ALA BOOK...
C ..... FIRST BBT*DB (12X3)*(3X3)
        DO I=1,12
          DO J=1,3
            SUM=0.D0
            DO K=1,3
              SUM=SUM+BB(K,I)*DB(J,K)
            ENDDO
            BTDB(J,I)=SUM
          ENDDO
        ENDDO
C THEN (BBT*DB)*BB (12X12)*(3X12)
        DO I=1,12
          DO J=1,12
            SUM=0.D0
            DO K=1,3
              SUM=SUM+BTDB(K,I)*BB(K,J)
            ENDDO
            EKB(I,J)=EKB(I,J)+FAC*SUM
          ENDDO
        ENDDO
      ENDDO
C END OF 2X2 INTEGRATION LOOP
C
C.....BEGIN 1X1 REDUCED INTEGRATION LOOP
      DO L=1,1
        XS(1)=0.D0
        XS(2)=0.D0
        CALL PGENSHAPE(N,XS,XY,PSI,DPSI,DETJ,DPSIX,DPSIY)
        CALL PLOADB(PSI,DPSIX,DPSIY,BB,BS)
C
C.......ACCUMULATE INTEGRATION POINT VALUE OF INTEGRALS
        FAC=DETJ*4.D0
C ..... FORM KS STIFFNESS MATRIX ........................
C ..... OK, NOW HERE GOES THE MATRIX MULTIPLICATION THAT 
C ..... GENERATES BST*DS*BS ALA BOOK...
C ..... FIRST BST*DS (12X2)*(2X2)
        DO I=1,12
          DO J=1,2
            SUM=0.D0
            DO K=1,2
              SUM=SUM+BS(K,I)*DS(J,K)
            ENDDO
            BTDB(J,I)=SUM
          ENDDO
        ENDDO
C THEN (BBT*DB)*BB (12X12)*(2X12)
        DO I=1,12
          DO J=1,12
            SUM=0.D0
            DO K=1,2
              SUM=SUM+BTDB(K,I)*BS(K,J)
            ENDDO
            EKS(I,J)=EKS(I,J)+FAC*SUM
          ENDDO
        ENDDO
      ENDDO
C END OF 1X1 INTEGRATION LOOP
C
C ... COMBINE KB AND KS STIFFNESS MATRICES .......
      DO I=1,12
        DO J=1,12
          EKB(I,J)=EKB(I,J)+EKS(I,J)+EKH(I,J)
        ENDDO
      ENDDO
C      PAUSE
C.....CALCULATE LOWER SYMMETRIC PART OF EK
C      DO I=1,4
C        DO J=1,I
C          EK(I,J)=EK(J,I)
C        ENDDO
C      ENDDO
      RETURN
1000  FORMAT(1X,'BAD JACOBIAN',E10.3)
1001  FORMAT(A,/,(2I5,3X,1PE13.6))
1002  FORMAT(A,/,(1P2E13.6))
      END
C=============================================================
      SUBROUTINE PASSMB(N3,NZ,nnode,GK,GF,EK,EC,EF,N,NODE,
     &                  nit,ia,nmax1,ja,iwork)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GK(nit),GF(N3),ia(nmax1),ja(nit),iwork(nit)
      DIMENSION EK(12,12),EC(12,12),EF(12),NODE(N)
C........................
C||||||||||||||||||||||||
C     PRINT *,'IN PASSMB'
      DO L=1,N
        I=NODE(L)
C.......ASSEMBLE GLOBAL VECTOR GF
        GF(I)=GF(I)+EF(L)
        DO M=1,N
          J=NODE(M)
C.........ASSEMBLE GLOBAL STIFFNESS MATRIX GK
            call sbsij(nnode,nit,ia,ja,GK,iwork,i,j,ek(l,m),1,1,6,ier)
            call sbsij(nnode,nit,ia,ja,GK,iwork,i,j,ec(l,m),1,1,6,ier)
        ENDDO
      ENDDO
C||||||||||||||||||||||||
C........................
      END
C=============================================================
      SUBROUTINE PFORMKF(NMAX,N3,NZ,NUMEL,NNODE,
     &                  X,Y,THICK,KX,KKX,EK,EC,EF,XI,W,
     &                  GK,GF,DELT,
     &                  nit,ia,nmax1,ja,iwork)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NMAX),Y(NMAX),KX(NMAX,4),KKX(NMAX,12)
      DIMENSION GK(nit),GF(N3),ia(nmax1)ja(nit),iwork(nit)
      DIMENSION THICK(NMAX)
      DIMENSION XI(2,4),W(4)
      DIMENSION EK(12,12),EC(12,12),EF(12)
      DIMENSION LM(4),XY(2,4),LLM(12)
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
      DATA N /4/,NL /4/
      DO NEL=1,NUMEL
        ETHICK=0.D0
        DO L=1,4
          LM(L)=KX(NEL,L)
          XY(1,L)=X(LM(L))
          XY(2,L)=Y(LM(L))
          ETHICK=ETHICK+THICK(LM(L))
        ENDDO
        DO L=1,12
          LLM(L)=KKX(NEL,L)
        ENDDO
        ETHICK=RHOG*ETHICK/DBLE(4)
        CALL PELEM(nel,ETHICK,XY,4,EK,EC,EF,NL,XI,W,DELT)
C..................................
        IF(.false.) THEN
          write(7,*) NEL
          DO I=1,12
            write(7,*) (EK(I,J),J=1,12),EF(I)
          ENDDO
c          PAUSE
        ENDIF
C..................................
        CALL PASSMB(N3,NZ,nnode,GK,GF,EK,EC,EF,12,LLM,
     &                  nit,ia,nmax1,ja,iwork)
C..................................
C        DO I=1,NNODE
C          PRINT 1000,(GK(I,J),J=1,NNODE),GF(I)
C        ENDDO
C        PAUSE
C..................................
      ENDDO
C|||||||||||||||||||||||||||||||||||
C...................................
1000  FORMAT(1X,10F8.3)
      END
C=============================================================
      SUBROUTINE PSHAPE(XS,PSI,DPSI)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XS(2),PSI(4),DPSI(4,2)
C...................................
C|||||||||||||||||||||||||||||||||||
      ETA=XS(1)
      RNU=XS(2)
      PX0=1.D0-ETA
      PX1=1.D0+ETA
      PY0=1.D0-RNU
      PY1=1.D0+RNU
      PSI(1)=PX0*PY0/4.D0
      PSI(2)=PX1*PY0/4.D0
      PSI(3)=PX1*PY1/4.D0
      PSI(4)=PX0*PY1/4.D0
      DPSI(1,1)=-PY0/4.D0
      DPSI(2,1)=-DPSI(1,1)
      DPSI(3,1)=PY1/4.D0
      DPSI(4,1)=-DPSI(3,1)
      DPSI(1,2)=-PX0/4.D0
      DPSI(2,2)=-PX1/4.D0
      DPSI(3,2)=-DPSI(2,2)
      DPSI(4,2)=-DPSI(1,2)
C|||||||||||||||||||||||||||||||||||
C...................................
      END
C=============================================================
      SUBROUTINE PLOADB(PSI,DPSIX,DPSIY,BB,BS)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PSI(4),DPSIX(4),DPSIY(4)
      DIMENSION BB(3,12),BS(2,12)
      DO IA=1,4
        IP1=3*IA-2
        IP2=IP1+1
        IP3=IP1+2
        BB(1,IP1)=0.D0
        BB(2,IP1)=0.D0
        BB(3,IP1)=0.D0
        BB(1,IP2)=DPSIX(IA)
        BB(2,IP2)=0.D0
        BB(3,IP2)=DPSIY(IA)
        BB(1,IP3)=0.D0
        BB(2,IP3)=DPSIY(IA)
        BB(3,IP3)=DPSIX(IA)
        BS(1,IP1)=DPSIX(IA)
        BS(2,IP1)=DPSIY(IA)
        BS(1,IP2)=-PSI(IA)
        BS(2,IP2)=0.D0
        BS(1,IP3)=0.D0
        BS(2,IP3)=-PSI(IA)
      ENDDO
      END
C=============================================================
      SUBROUTINE PCONNECT(NMAX,NUMNP,NUMEL,KX,KKX,LTEMP)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION KX(NMAX,4),KKX(NMAX,12),LTEMP(NMAX,3)
      IC=1
      DO I=1,NUMNP
        LTEMP(I,1)=IC
        LTEMP(I,2)=IC+1
        LTEMP(I,3)=IC+2
        IC=IC+3
C        PRINT *,I,(LTEMP(I,J),J=1,3)
      ENDDO
      DO I=1,NUMEL
        DO J=1,4
          IP1=3*J-2
          IP2=IP1+1
          IP3=IP1+2
          KKX(I,IP1)=LTEMP(KX(I,J),1)
          KKX(I,IP2)=LTEMP(KX(I,J),2)
          KKX(I,IP3)=LTEMP(KX(I,J),3)
        ENDDO
C        PRINT 1000,I,(KX(I,J),J=1,4)
C        PRINT 1000,I,(KKX(I,J),J=1,12)
C        PAUSE
      ENDDO
1000  FORMAT(1X,13I5)
      END
C=============================================================
      SUBROUTINE POUTSTF(NNODE,WWW,WRATE,NUMNP,THICK,NUMCOL,NUMLEV,TIME)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=10404,N3=3*NMAX)
      DIMENSION WWW(NNODE),THICK(NUMNP),TMP(N3),WRATE(NNODE)
      CHARACTER*80 HED
C ... DEPRESSION ...
      DO I=1,NNODE
        TMP(I)=WWW(I)
        IF(WWW(I).GT.0) TMP(I)=TMP(I)*1
      ENDDO
      WRITE(81,1000)
      WRITE(81,1001) NUMCOL,NUMLEV
      WRITE(81,1002) 0,NUMCOL*1000,0,NUMLEV*1000
      WRITE(81,1003) 'DEPRESSION'
      WRITE(81,1008) TIME
      WRITE(81,1004)
      WRITE(81,1005)
      WRITE(81,1006) (TMP(I),I=1,NNODE,3)
      WRITE(81,1007)
C ... LOAD ...
      WRITE(80,1000)
      WRITE(80,1001) NUMCOL,NUMLEV
      WRITE(80,1002) 0,NUMCOL*1000,0,NUMLEV*1000
      WRITE(80,1003) 'LOAD'
      WRITE(80,1008) TIME
      WRITE(80,1004)
      WRITE(80,1005)
      WRITE(80,1006) (THICK(I),I=1,NUMNP)
      WRITE(80,1007)
C ... RATE ... mm/yr
      DO I=1,NNODE
        TMP(I)=WRATE(I)*1000.D0
      ENDDO
      WRITE(82,1000)
      WRITE(82,1001) NUMCOL,NUMLEV
      WRITE(82,1002) 0,NUMCOL*1000,0,NUMLEV*1000
      WRITE(82,1003) 'RATE'
      WRITE(82,1008) TIME
      WRITE(82,1004)
      WRITE(82,1005)
      WRITE(82,1006) (TMP(I),I=1,NNODE,3)
      WRITE(82,1007)
      WRITE(HED,*) ' TIME=',NINT(TIME)
      WRITE(88) HED
      WRITE(88) (WWW(I),I=1,NNODE,3)
      WRITE(88) (THICK(I),I=1,NUMNP)
      WRITE(88) (1000*WRATE(I),I=1,NNODE,3)
C ... 
1000  FORMAT(1X,'RANK 2')
1001  FORMAT(1X,'DIMENSIONS',2I7)
1002  FORMAT(1X,'BOUNDS',4I13)
1003  FORMAT(1X,'NAME ',A)
1004  FORMAT(1X,'SCALAR')
1005  FORMAT(1X,'DATA')
1006  FORMAT(1X,1P5G14.6)      
1007  FORMAT(1X,'END')
1008  FORMAT(1X,'TIME',G13.6)
1033  FORMAT(1X,'NAME LOAD')
      END
C=============================================================
      SUBROUTINE PTIMEDEP(N3,NZ,NNODE,GK,GC,GF,KA,KZ,WWW,DELT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GK(N3,NZ),GC(N3,NZ),GF(N3),WWW(N3)
      DIMENSION KA(N3,NZ+1),KZ(N3)
      IF(DELT.GT.0.D0) THEN
C ..... FORM MODIFIED STIFFNESS AND LOAD .......
        DELT1=1.D0/DELT    
        DO I = 1,NNODE
          DO J = 1,KZ(I)
            JG=KA(I,J)
C           PRINT *,I,J,GF(I),GC(I,J)*WWW(JG)*DELT1
            GF(I)=GF(I)+GC(I,J)*WWW(JG)*DELT1
            GK(I,J)=GK(I,J)+GC(I,J)*DELT1
          ENDDO 
        ENDDO
      ENDIF
      END
C=============================================================
      SUBROUTINE PGAUSEID(NMAX,NZ,N,EPS,AA,KA,B,X)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION AA(NMAX,NZ),KA(NMAX,NZ+1),B(NMAX),X(NMAX)
      REAL*8 SUM
C      DATA EPS /1D-6/, ITMAX /100/, TOL /0.D0/
      DATA ITMAX /100/, TOL /0.D0/
C ... A DIFFRENT FIRST GUESS ...
C     IF(.FALSE.) THEN
C       DO I=1,N
C         X(I)=B(I)/AA(I,1)
C       ENDDO
C     ENDIF
C ... ..........................
      XMAX=-1D30
      DO I=1,N
        JMAX=KA(I,NZ+1)
        SUM=B(I)
        DO J=2,JMAX
          SUM=SUM-AA(I,J)*X(KA(I,J))
        ENDDO
        XNEW=SUM/AA(I,1)
        XMAX=MAX(XMAX,ABS(X(I)-XNEW))
        X(I)=XNEW
      ENDDO
      ERRORG=PRESID(NMAX,NZ,N,AA,KA,B,X)
C     PRINT 1001,   0,ERRORG
      WRITE(7,1001) 0,ERRORG
C ... RARE BUT POSSIBLE RETURN
      IF(ERRORG.EQ.0.) RETURN
      DO ITER=1,ITMAX
      XMAX=-1D30
        DO I=1,N
          JMAX=KA(I,NZ+1)
          SUM=B(I)
          DO J=2,JMAX
            SUM=SUM-AA(I,J)*X(KA(I,J))
          ENDDO
          XNEW=SUM/AA(I,1)
          XMAX=MAX(XMAX,ABS(X(I)-XNEW))
          X(I)=XNEW
        ENDDO
C       IF(.TRUE.) THEN
C         PRINT 1001,ITER,XMAX
C         PRINT 1000,(X(I),I=1,N)
C         WRITE(7,1001) ITER,XMAX
C         WRITE(7,1000) (X(I),I=1,N)
C       ENDIF
        ERROR=PRESID(NMAX,NZ,N,AA,KA,B,X)
        RATIO=ERROR/ERRORG
C       PRINT 1001,   ITER,ERROR,RATIO,XMAX
        WRITE(7,1001) ITER,ERROR,RATIO,XMAX
C ..... A NORMAL RETURN WHEN CONVERGED ...
        IF(ABS(RATIO).LT.EPS) THEN
C         PRINT 1002,   'A:CONVERGED',ITER,ERROR,RATIO,XMAX
          WRITE(7,1002) 'A:CONVERGED',ITER,ERROR,RATIO,XMAX
          RETURN
        ENDIF
C ..... ANOTHER NORMAL RETURN WHEN CONVERGED ...
        IF(XMAX.LT.TOL) THEN
C         PRINT 1002,   'B:CONVERGED',ITER,ERROR,RATIO,XMAX
          WRITE(7,1002) 'B:CONVERGED',ITER,ERROR,RATIO,XMAX
          RETURN
        ENDIF
      ENDDO
      PRINT *,   'DIDNOT CONVERGE IN ',ITMAX,''
      WRITE(7,*) 'DIDNOT CONVERGE IN ',ITMAX
C     PAUSE
1000  FORMAT(1X,'PGAUSEID:',(1X,T10,5G14.6))
1001  FORMAT(1X,'PGAUSEID:',T21,I5,3G13.6)
1002  FORMAT(1X,'PGAUSEID:',A,I5,3G13.6)
      END
C=======================================
      FUNCTION PRESID(NMAX,NZ,N,AA,KA,B,X)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION AA(NMAX,NZ),KA(NMAX,NZ+1),B(NMAX),X(NMAX)
      REAL*8 SUMSQ,SUM
      SUMSQ=0.D0
      DO I=1,N
        JMAX=KA(I,NZ+1)
        SUM=0.D0
        DO J=1,JMAX
          SUM=SUM+AA(I,J)*X(KA(I,J))
        ENDDO
        SUM=SUM-B(I)
        SUMSQ=SUMSQ+SUM**2
      ENDDO
      PRESID=SUMSQ
      END
C=======================================
      SUBROUTINE PGENSHAPE(N,XS,XY,PSI,DPSI,DETJ,DPSIX,DPSIY)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION XY(2,N)
      DIMENSION DPSIX(4),DPSIY(4),DXDS(2,2),DSDX(2,2)
      DIMENSION PSI(4),DPSI(4,2),XS(2)
        CALL PSHAPE(XS,PSI,DPSI)
C.......CALCULATE DSDS...EQUATION(5.3.6)
        DO I=1,2
          DO J=1,2
            DXDS(I,J)=0.0D0
            DO K=1,4
              DXDS(I,J)=DXDS(I,J)+DPSI(K,J)*XY(I,K)
            ENDDO
          ENDDO
        ENDDO
C.......CALCULATE DSDX...EQUATION(5.2.7)
        DETJ=DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1)
        IF(DETJ.LE.0.0) THEN
          WRITE(*,1000) DETJ
          WRITE(*,1001) 'DXDS',((I,J,DXDS(I,J),I=1,2),J=1,2)
          WRITE(*,1001) 'DPSI',((I,J,DPSI(J,I),I=1,2),J=1,4)
          WRITE(*,1002) 'XY',((XY(I,J),I=1,2),J=1,N)
          STOP
        ENDIF
        DSDX(1,1)=DXDS(2,2)/DETJ    
        DSDX(2,2)=DXDS(1,1)/DETJ    
        DSDX(1,2)=-DXDS(1,2)/DETJ    
        DSDX(2,1)=-DXDS(2,1)/DETJ  
C.......CALCULATE D(PSI)/DX...EQUATION(5.3.5)
        DO I=1,N
          DPSIX(I)=DPSI(I,1)*DSDX(1,1)+DPSI(I,2)*DSDX(2,1)  
          DPSIY(I)=DPSI(I,1)*DSDX(1,2)+DPSI(I,2)*DSDX(2,2) 
        ENDDO
1000  FORMAT(1X,'BAD JACOBIAN',E10.3)
1001  FORMAT(A,/,(2I5,3X,1PE13.6))
1002  FORMAT(A,/,(1P2E13.6))
      END
C=======================================
      SUBROUTINE PDUMPMAT(N3,NZ,NNODE,KZ,KA,GK,GF,WWW,DISP)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      LOGICAL DISP
      DIMENSION GK(N3,NZ),GF(N3),WWW(N3)
      DIMENSION KA(N3,NZ+1),KZ(N3)
        IF(.false.) THEN
          DO I=1,NNODE
            write(7,*) 'DIAGONAL='
            write(7,*) KZ(I),I,GK(I,1)
            SUM=0
            DO J=2,KZ(I)
              SUM=SUM+ABS(GK(I,J))
              JJ=KA(I,J)
              write(7,*) I,JJ,GK(I,J),GK(I,J)/GK(I,1)
            ENDDO
            IF(SUM/GK(I,1).GT.1.1) THEN
              write(7,*) 'EQ:',I,' NOT DIAG DOM.',GK(I,1)/SUM
c              PAUSE
            ENDIF
c            PAUSE
            write(7,*) I,GF(I)
          ENDDO
c          PAUSE
        ENDIF
        IF(DISP) THEN
          write(13,*) nnode
          DO I=1,NNODE
            write(13,*) (gk(i,j),J=1,nz)
            write(13,*) (ka(i,j),J=1,nz+1)
            write(13,*) gf(i),WWW(I)
          ENDDO
        ENDIF
      END
C=======================================
C=======================================
      SUBROUTINE PSETMAT
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
c ... MATERIAL PROPERTIES AND THICKNESS OF CRUST
C ... YOUNGS MODULUS (NT/M**2)
      E=5.D10
      E=1.D11
C ... POISONS RATIO (DIMENSIONLESS)
      RNU=0.5D0
C ... LAMES COEFFICIENTS
      RLAMBDA=RNU*E/(1.D0-RNU**2)
      RMU=0.5D0*E/(1.D0+RNU)
C      RMU=10000.D0
C      RLAMBDA=10000.D0
C ... THICKNESS OF CRUST (90-130 KM) (110,000 M)
      TTT=110.D0*1000.D0
C ... RHOICE*GRAV (NT/M**3) TIMES THICKNESS, YIELDS PRESSURE, NT/M**2
      GRAV=9.8d0*0.3816d0
      RHOR=4000.D0
      RHOW=1092.D0
      RHOI=917.D0
      RHOG=RHOI*GRAV
      ROCKICE=RHOR/RHOI
C ... RELAXATION TIME CONSTANT ......!
      RELAX=6000.D0                  !
      CTIME=3.16516D8*relax/6000.d0  !
      END
