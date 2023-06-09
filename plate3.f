      SUBROUTINE PLATE(ITIME,NMAX,N3,NUMNP,NUMEL,X,Y,KX,THICK,
     &                 DELT,WWW,WRATE,WMIN,TIME,WWWORIG)
c-----------------------------------------------------------------------
c 4th order plate solver with one-time generation of stiffness and
c capacitance matrix. uses my sparse matrix storage for the static matrice
c and ITPACK sprse storage for the time-dependent modified matrices, and JCG
c iterative matrix solver. j fastook march 11, 1999
c-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)                                          
      PARAMETER(MMAX=12999,NZ=27, NZ1=NZ+1, M3=3*MMAX)
      parameter(nit=m3*nz,nmax1=m3+1,nmax3=m3*3)
      parameter(itmax=m3/10,ncg=4*itmax,nw=4*m3+ncg)
c ....arrays for ITPACK sparse storage...................
      dimension GK(nit),ja(nit),ia(nmax1)
      dimension iwksp(nmax3),wksp(nw),iwork(nit)
      dimension iparm(12),rparm(12)
      DIMENSION THICK(NMAX),WWW(N3),X(NMAX),Y(NMAX),KX(NMAX,4)
      DIMENSION WRATE(N3,2),WWWORIG(NMAX)
      DIMENSION KKX(MMAX,12),LTEMP(MMAX,3)
      DIMENSION XI(2,4),W(4)
      DIMENSION EK(12,12),EC(12,12),EF(12)
      DIMENSION GF(M3)
      DIMENSION GK0(M3,NZ),GC0(M3,NZ)
      DIMENSION KA(M3,NZ1),KZ(M3)
C     DIMENSION WWWSAVE(MMAX)
      CHARACTER*80 HED
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
C ... TIMER STUFF, FOR SGI ONLY ...!
      REAL*4 TB(2)
c     REAL*4 TB(2),ETIME,DTIME     !
c     EXTERNAL ETIME,DTIME         !
C .................................!
      SAVE IPASS,KKX,WSAVE,NNODE,XI,W,GK0,GC0,KA,KZ
C     SAVE WWWSAVE,iparm,rparm,IREAD
      DATA IPASS /0/,IREAD /0/
1000  FORMAT(1X,A,T25,1PG13.6,G13.6)
1001  FORMAT(1X,1P6(1X,G12.5))
c      WRITE(7,1000) ' TIME BEGINNING PLATE ',ETIME(TB),DTIME(TB)
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
        READ(40,*,IOSTAT=IOK) NNN
        IF(IOK.EQ.0) THEN
          IF(NNN.NE.NNODE) THEN
            PRINT *,'PROBLEMS:'
            PRINT *,'incompatible with current NNODE=',NNODE,NNN
            IOK=1
          ENDIF
          DO I=1,NNODE
            READ(40,*) II,WWW(I)
            WRATE(I,2)=WWW(I)
            IF(I.NE.II) THEN
              PRINT *,'PROBLEMS:'
              PRINT *,'incompatible with current (WWW)'
              IOK=1
            ENDIF
          ENDDO
          DO I=1,NUMNP
            READ(40,*) II,WWWORIG(I)
            IF(I.NE.II) THEN
              PRINT *,'PROBLEMS:'
              PRINT *,'incompatible with current (WWWORIG)'
              IOK=1
            ENDIF
          ENDDO
          PRINT *,' BEDROCK DEPRESSION FILE FOUND'
          PRINT *,'    AND READ SUCCESSFULLY '
          if(itime.lt.0) then
            print *,'abandoning unloading'
            rewind 40
            return
          endif
          IREAD=1
        ENDIF
        IF(IOK.NE.0) THEN
          PRINT *,' NONE FOUND, SET TO ZERO ... '
          DO I=1,NUMNP
            WWW(I)=0.D0
            WRATE(I,2)=WWW(I)
          ENDDO
          WRITE(HED,*) ' TIME=',NINT(TIME-DELT)
          WRITE(88) HED
          WRITE(88) (WWW(I),I=1,NNODE,3)
          WRITE(88) (THICK(I),I=1,NUMNP)
          WRITE(88) (1000*WRATE(I,1),I=1,NNODE,3)
        ENDIF
        CALL PFORMGK(NMAX,N3,NZ,NUMEL,NNODE,
     &              X,Y,KX,KKX,EK,EC,XI,W,
     &              GK0,GC0,KA,KZ)
      ENDIF
C ......................................................
C ....FORM STIFFNESS AND LOAD ..............
      CALL PFORMGF(NMAX,N3,NUMEL,NNODE,
     &            X,Y,THICK,KX,KKX,EF,XI,W,
     &            GF)
C ....TIME DEPENDENT CASE AND VARYING LOAD .......
      CALL PTIMEDEP(N3,NZ,NNODE,GK0,GC0,GK,GF,KA,KZ,
     &              WWW,DELT,nit,nmax1,ia,ja,iwork)
C ......................................................
C.....DUMP MATRIX FOR EXAMINATION ......................
      CALL PDUMPMAT(N3,NZ,NNODE,KZ,KA,GK,GF,WWW,.false.)
C.......................................................
C     IF(.TRUE.) THEN
C.......SOLVE EQUATIONS WITH JORDAN CONJUGATE-GRADIENT ITPACK ....!
c        WRITE(7,1000) ' TIME BEFORE JCG ',ETIME(TB),DTIME(TB)!
        call dfault(iparm,rparm)
        if(delt.ne.0.0) then
          iparm(1)=100  ! MAX number of iteration
        else
          iparm(1)=1000 ! MAX number of iteration
        endif
        rparm(1)=1d-6   ! ZETA, stopping criteria
c       iparm(2)=-1     ! LEVEL of output (-1:none)
c       iparm(4)=7      ! OUTPUT unit number
        iparm(5)=1      ! NONsymmetric matrix (1)
        iparm(6)=0      ! NON-ADAPTIVE (0) (1 doesnt work)
c       iparm(10)=1     ! REMOVES large diagonal entries (doesnt work)
C        DO I=1,12
C          PRINT *,I,IPARM(I),RPARM(I)
C        ENDDO
        CALL JCG(NNODE,ia,ja,GK,GF,WWW,
     &         iwksp,nw,wksp,iparm,rparm,ier)
C        DO I=1,12
C          PRINT *,I,IPARM(I),RPARM(I)
C        ENDDO
        if(ier.ne.0) then
          print *,'JCG error:',ier
          pause
        endif
        print *,' relative error=',rparm(1),' in iterations = ',iparm(1)
        print *,rparm(11),rparm(12)
c        WRITE(7,1000) ' TIME AFTER JCG ',ETIME(TB),DTIME(TB) !
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
      IF(IREAD.EQ.0 .AND. IPASS.EQ.0) THEN
        II=1
        DO I=1,NUMNP
          WWWORIG(I)=WWW(II)
          II=II+3
        ENDDO
      ENDIF
      write(*,*) '*******************************************'
      WRITE(*,1001) TIME,REAL(THTOT),REAL(-THTOT/WTOT/ROCKICE),
     &              REAL(WMIN),REAL(WMAX),REAL(WSAVE-WMIN)
      write(*,*) '*******************************************'
      WRITE(92,*) TIME,WMIN
      WSAVE=WMIN
      IPASS=1
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
      SUBROUTINE PELEMEK(XY,N,EKB,EC,NL,XI,W)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XY(2,N),EKB(12,12),EKS(12,12),EKH(12,12)
      DIMENSION EC(12,12)
      DIMENSION DPSIX(4),DPSIY(4)
      DIMENSION PSI(4),DPSI(4,2),XS(2)
      DIMENSION BB(3,12),BS(2,12),DB(3,3),DS(2,2)
      DIMENSION BTDB(3,12)
      DIMENSION XI(2,4),W(4)
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
C.....INITIALIZE ELEMENT ARRAYS
      DO I=1,12
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
C  .... RATIOS OF DENSITY OF ROCK AND ICE FOR OVERBURDEN LOAD
        XB=ROCKICE*RHOG
C ..... FORM CAPACITANCE MATRIX .....................
        DO I=1,4
          IP1=3*I-2
          IP2=3*I-1
          IP3=3*I
          DO J=1,4
            JP1=3*J-2
            JP2=3*J-1
            JP3=3*J
            EKH(IP1,JP1)=EKH(IP1,JP1)+FAC*XB*PSI(I)*PSI(J)
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
      SUBROUTINE PELEMGF(ETHICK,XY,N,EF,NL,XI,W)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XY(2,N)
      DIMENSION EF(12)
      DIMENSION DPSIX(4),DPSIY(4)
      DIMENSION PSI(4),DPSI(4,2),XS(2)
      DIMENSION XI(2,4),W(4)
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
C.....INITIALIZE ELEMENT ARRAYS
      DO I=1,12
        EF(I)=0.0D0
      ENDDO
C.....BEGIN 2X2 INTEGRATION LOOP
      DO L=1,NL
        XS(1)=XI(1,L)
        XS(2)=XI(2,L)
        CALL PGENSHAPE(N,XS,XY,PSI,DPSI,DETJ,DPSIX,DPSIY)
C
C.......ACCUMULATE INTEGRATION POINT VALUE OF INTEGRALS
        FAC=DETJ*W(L)
        XF=-ETHICK
C ..... FORM CAPACITANCE MATRIX .....................
        DO I=1,4
          IP1=3*I-2
          EF(IP1)=EF(IP1)+XF*PSI(I)*FAC
        ENDDO
      ENDDO
      RETURN
1000  FORMAT(1X,'BAD JACOBIAN',E10.3)
1001  FORMAT(A,/,(2I5,3X,1PE13.6))
1002  FORMAT(A,/,(1P2E13.6))
      END
C=============================================================
      SUBROUTINE PASSMBGK(N3,NZ,GK0,GC0,KA,KZ,EK,EC,N,NODE)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GK0(N3,NZ),GC0(N3,NZ)
      DIMENSION KA(N3,NZ+1),KZ(N3)
      DIMENSION EK(12,12),EC(12,12),NODE(N)
C........................
C||||||||||||||||||||||||
C     PRINT *,'IN PASSMB'
      DO L=1,N
        I=NODE(L)
        DO M=1,N
          J=NODE(M)
C.........ASSEMBLE GLOBAL STIFFNESS MATRIX GK
            IF(I.EQ.J) THEN
              GK0(I,1)=GK0(I,1)+EK(L,M)
              GC0(I,1)=GC0(I,1)+EC(L,M)
              KA(I,1)=I
            ELSE
              DO K=2,KZ(I)
                IF(KA(I,K).EQ.J) THEN
                  GK0(I,K)=GK0(I,K)+EK(L,M)
                  GC0(I,K)=GC0(I,K)+EK(L,M)
                  GOTO 99
                ENDIF
              ENDDO
              KZ(I)=KZ(I)+1
              GK0(I,KZ(I))=EK(L,M)
              GC0(I,KZ(I))=EC(L,M)
              KA(I,KZ(I))=J
            ENDIF
99          CONTINUE
        ENDDO
      ENDDO
C||||||||||||||||||||||||
C........................
      END
C=============================================================
      SUBROUTINE PASSMBGF(N3,GF,EF,N,NODE)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GF(N3)
      DIMENSION EF(12),NODE(N)
C........................
C||||||||||||||||||||||||
C     PRINT *,'IN PASSMBGF'
      DO L=1,N
        I=NODE(L)
C.......ASSEMBLE GLOBAL VECTOR GF
        GF(I)=GF(I)+EF(L)
      ENDDO
C||||||||||||||||||||||||
C........................
      END
C=============================================================
      SUBROUTINE PFORMGK(NMAX,N3,NZ,NUMEL,NNODE,
     &                  X,Y,KX,KKX,EK,EC,XI,W,
     &                  GK0,GC0,KA,KZ)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NMAX),Y(NMAX),KX(NMAX,4),KKX(NMAX,12)
      DIMENSION GK0(N3,NZ),GC0(N3,NZ)
      DIMENSION KA(N3,NZ+1),KZ(N3)
      DIMENSION XI(2,4),W(4)
      DIMENSION EK(12,12),EC(12,12)
      DIMENSION LM(4),XY(2,4),LLM(12)
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
      DO I=1,NNODE
        KZ(I)=1
        KA(I,1)=I
        KA(I,NZ+1)=0
        DO J=1,NZ
          KA(I,J)=0
          GK0(I,J)=0.0D0
          GC0(I,J)=0.0D0
        ENDDO
      ENDDO
      DO NEL=1,NUMEL
        DO L=1,4
          LM(L)=KX(NEL,L)
          XY(1,L)=X(LM(L))
          XY(2,L)=Y(LM(L))
        ENDDO
        DO L=1,12
          LLM(L)=KKX(NEL,L)
        ENDDO
        CALL PELEMEK(XY,4,EK,EC,4,XI,W)
C..................................
        IF(.false.) THEN
          write(7,*) NEL
          DO I=1,12
            write(7,*) (EK(I,J),J=1,12)
          ENDDO
c          PAUSE
        ENDIF
C..................................
        CALL PASSMBGK(N3,NZ,GK0,GC0,KA,KZ,EK,EC,12,LLM)
C..................................
C        DO I=1,NNODE
C          PRINT 1000,(GK0(I,J),J=1,NNODE)
C        ENDDO
C        PAUSE
C..................................
      ENDDO
      DO I=1,NNODE
C        PRINT *,KZ(I),KA(I,NZ+1)
        KA(I,NZ+1)=KZ(I)
      ENDDO
C|||||||||||||||||||||||||||||||||||
C...................................
1000  FORMAT(1X,10F8.3)
      END
C=============================================================
      SUBROUTINE PFORMGF(NMAX,N3,NUMEL,NNODE,
     &                  X,Y,THICK,KX,KKX,EF,XI,W,
     &                  GF)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NMAX),Y(NMAX),KX(NMAX,4),KKX(NMAX,12)
      DIMENSION GF(N3)
      DIMENSION THICK(NMAX)
      DIMENSION XI(2,4),W(4)
      DIMENSION EF(12)
      DIMENSION LM(4),XY(2,4),LLM(12)
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
      DO I=1,NNODE
        GF(I)=0.D0
      ENDDO
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
        CALL PELEMGF(ETHICK,XY,4,EF,4,XI,W)
C..................................
        IF(.false.) THEN
          write(7,*) NEL
          DO I=1,12
            write(7,*) EF(I)
          ENDDO
c          PAUSE
        ENDIF
C..................................
        CALL PASSMBGF(N3,GF,EF,12,LLM)
C..................................
C        DO I=1,NNODE
C          PRINT 1000,GF(I)
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
      SUBROUTINE PTIMEDEP(N3,NZ,NNODE,GK0,GC0,GK,GF,KA,KZ,
     &                    WWW,DELT,nit,nmax1,ia,ja,iwork)
      IMPLICIT REAL*8(A-H,O-Z)
c ....arrays for ITPACK sparse storage...................
      dimension GK(nit),ja(nit),ia(nmax1)
      dimension iwork(nit)
      DIMENSION GF(N3),WWW(N3)
      DIMENSION GK0(N3,NZ),GC0(N3,NZ)
      DIMENSION KA(N3,NZ+1),KZ(N3)
      IF(DELT.GT.0.D0) THEN
        call sbini(NNODE,nit,ia,ja,GK,iwork)
C ..... FORM MODIFIED STIFFNESS AND LOAD .......
        DELT1=1.D0/DELT    
        DO I = 1,NNODE
          DO J = 1,KZ(I)
            JG=KA(I,J)
C           PRINT *,I,J,GF(I),GC0(I,J)*WWW(JG)*DELT1
            GF(I)=GF(I)+GC0(I,J)*WWW(JG)*DELT1
            GKIJ=GK0(I,J)+GC0(I,J)*DELT1
            call sbsij(NNODE,nit,ia,ja,GK,iwork,I,JG,GKIJ,0,1,6,ier)
          ENDDO 
        ENDDO
        call sbend(NNODE,nit,ia,ja,GK,iwork)
      else
        call sbini(NNODE,nit,ia,ja,GK,iwork)
        DO I = 1,NNODE
          DO J = 1,KZ(I)
            JG=KA(I,J)
            GKIJ=GK0(I,J)
            call sbsij(NNODE,nit,ia,ja,GK,iwork,I,JG,GKIJ,0,1,6,ier)
          ENDDO 
        ENDDO
        call sbend(NNODE,nit,ia,ja,GK,iwork)
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
