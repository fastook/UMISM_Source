      SUBROUTINE OPLATE(ITIME,NMAX,N3,NUMNP,NUMEL,X,Y,KX,THICK,
     &                 DELT,WWW,WRATE,WMIN,TIME,WWWORIG)
c-----------------------------------------------------------------------
c 4th order plate solver with one-time generation of stiffness and
c capacitance matrix. uses my sparse matrix storage for the static matrice
c and ITPACK sprse storage for the time-dependent modified matrices, and JCG
c iterative matrix solver. j fastook march 11, 1999
c-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)                                          
      include "parameter.h"
      PARAMETER(MMAX=MAXNUM,NZ=27, NZ1=NZ+1, M3=3*MMAX)
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
      LOGICAL CTOGG,ITOGG,IOTOGG
      INTEGER BTOGG,SLTOGG,WTOGG
      COMMON /TOGGLES/ SLTOGG,CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(1000)
      COMMON /IOLIST/ LIST
      common /oelastic/ wwwo(M3),worate(M3,2),wwwoorig(MMAX)
      common /elastic/ wwwe(M3),werate(M3,2),wwweorig(MMAX)
      common /viscous/ wwwv(M3),wvrate(M3,2),wwwvorig(MMAX)
      logical file40

C ... TIMER STUFF, FOR SGI ONLY ...!
      REAL*4 TB(2)
c     REAL*4 TB(2),ETIME,DTIME     !
c     EXTERNAL ETIME,DTIME         !
C .................................!
      SAVE IPASS,KKX,WSAVE,NNODE,XI,W,GK0,GC0,KA,KZ
C     SAVE WWWSAVE,iparm,rparm,IREAD
      DATA IPASS /0/,IREAD /0/
1000  FORMAT(1X,A,T25,1PG13.6,G13.6)
1001  FORMAT(1X,6(1X,1PG12.5))
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
        CALL OSETMAT
        NNODE=3*NUMNP
        PRINT *,'READING BEDROCK DEPRESSION FILE'
        inquire(file='fort.40',exist=file40)
        iok=-1
        if(file40) then
          READ(40,*,IOSTAT=IOK) NNN
        endif
        IF(IOK.EQ.0) THEN
          IF(NNN.NE.NNODE) THEN
            PRINT *,'PROBLEMS:'
            PRINT *,'incompatible with current NNODE=',NNODE,NNN
            IOK=1
          ENDIF
          do i=1,nnode
            read(40,*) ii,wwwe(i),wwwv(i)
            werate(i,2)=wwwe(i)
            wvrate(i,2)=wwwv(i)
            if(i.ne.ii) then
              print *,'problems:reading wwwe,wwwv'
              print *,'incompatible with current (www)'
              iok=1
            endif
          enddo
          do i=1,nnode
            read(40,*) ii,werate(I,1),wvrate(I,1)
            if(i.ne.ii) then
              print *,'problems:reading wrates(1)'
              print *,'incompatible with current (www)'
              iok=1
            endif
          enddo
          do i=1,nnode
            read(40,*) ii,werate(I,2),wvrate(I,2)
            if(i.ne.ii) then
              print *,'problems:reading wrates(2)'
              print *,'incompatible with current (www)'
              iok=1
            endif
          enddo
          do i=1,numnp
            read(40,*) ii,wwweorig(i),wwwvorig(i)
            if(i.ne.ii) then
              print *,'problems:reading wwworigs'
              print *,'incompatible with current (wwworig)'
              iok=1
            endif
          enddo
          do i=1,numnp*3
            www(i)=wwwe(i)+wwwv(i)
c           if(mod(i,30).eq.1) print *,i,wwwe(i),wwwv(i)
            wrate(i,1)=werate(i,1)+wvrate(i,1)
            wrate(i,2)=werate(i,2)+wvrate(i,2)
          enddo
          do i=1,numnp
            wwworig(i)=wwweorig(i)+wwwvorig(i)
          enddo
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
c$doacross local(i)
          DO I=1,NNODE
            WWW(I)=0.D0
            WRATE(I,2)=0.D0
          ENDDO
          WRITE(HED,*) ' TIME=',NINT(TIME-DELT)
          WRITE(88) HED
          WRITE(88) (WWW(I),I=1,NNODE,3)
          WRITE(88) (THICK(I),I=1,NUMNP)
          WRITE(88) (1000*WRATE(I,1),I=1,NNODE,3)
        ENDIF
        CALL OFORMGK(NMAX,N3,NZ,NUMEL,NNODE,
     &              X,Y,KX,KKX,EK,EC,XI,W,
     &              GK0,GC0,KA,KZ)
      ENDIF
C ......................................................
C ....FORM STIFFNESS AND LOAD ..............
      CALL OFORMGF(NMAX,N3,NUMEL,NNODE,
     &            X,Y,THICK,KX,KKX,EF,XI,W,
     &            GF)
C ....TIME DEPENDENT CASE AND VARYING LOAD .......
      CALL OTIMEDEP(N3,NZ,NNODE,GK0,GC0,GK,GF,KA,KZ,
     &              WWW,DELT,nit,nmax1,ia,ja,iwork)
C ......................................................
C.....DUMP MATRIX FOR EXAMINATION ......................
      CALL ODUMPMAT(N3,NZ,NNODE,KZ,KA,GK0,GF,WWW,.false.)
      CALL NDUMPMAT(NIT,GK,JA,NMAX1,IA,NNODE,GF,.false.)
C.......................................................
C     IF(.TRUE.) THEN
C.......SOLVE EQUATIONS WITH JORDAN CONJUGATE-GRADIENT ITPACK ....!
c        WRITE(7,1000) ' TIME BEFORE JCG ',ETIME(TB),DTIME(TB)!
        call dfault(iparm,rparm)
        if(delt.ne.0.0) then
          iparm(1)=150  ! MAX number of iteration
        else
          iparm(1)=1000 ! MAX number of iteration
        endif
        rparm(1)=1d-6   ! ZETA, stopping criteria
c       iparm(2)=2      ! LEVEL of output (-1:none)
c       iparm(4)=7      ! OUTPUT unit number
        iparm(5)=1      ! NONsymmetric matrix (1)
        iparm(6)=0      ! NON-ADAPTIVE (0) (1 doesnt work)
c       iparm(10)=1     ! REMOVES large diagonal entries (doesnt work)
C        DO I=1,12
C          PRINT *,I,IPARM(I),RPARM(I)
C        ENDDO
        CALL JCG(NNODE,ia,ja,GK,GF,WWW,
     &         iwksp,nw,wksp,iparm,rparm,ier)
C         DO I=1,12
C           PRINT *,I,IPARM(I),RPARM(I)
C         ENDDO
        IF(IOTOGG) THEN
          write(list(ipage+1),*) ' relative error=',rparm(1),
     &            ' in iterations = ',iparm(1)
          write(list(ipage+2),*) rparm(11),rparm(12)
          ipage=ipage+2
        ENDIF
        if(ier.ne.0) then
          CALL NDUMPMAT(NIT,GK,JA,NMAX1,IA,NNODE,GF,.TRUE.)
          print *,'JCG error:',ier
c          print '(1x,i3,1pg13.6)',(i,rparm(i),i=1,12)
          pause
        endif
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
c$doacross local(i)
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
        wwwo(i)=www(i)
        worate(i,1)=wrate(i,1)
        worate(i,2)=wrate(i,2)
      ENDDO
      THTOT=0.D0
      DO I=1,NUMNP
        THTOT=THTOT+THICK(I)
      ENDDO
      IF(IREAD.EQ.0 .AND. IPASS.EQ.0) THEN
        II=1
        DO I=1,NUMNP
          WWWORIG(I)=WWW(II)
          WWWOORIG(I)=WWW(II)
          II=II+3
        ENDDO
      ENDIF
      IF(IOTOGG) THEN
        write(list(ipage+1),*) 
     &       '*******************************************'
        WRITE(list(ipage+2),1001) TIME,REAL(THTOT),
     &              REAL(-THTOT/WTOT/ROCKICE),
     &              REAL(WMIN),REAL(WMAX),REAL(WSAVE-WMIN)
        write(list(ipage+3),*) 
     &       '*******************************************'
        ipage=ipage+3
      ENDIF
      WRITE(92,*) TIME,WMIN
      WSAVE=WMIN
      IPASS=1
      END
C=============================================================
      SUBROUTINE OELEMEK(XY,N,EKB,EC,NL,XI,W)
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
c$doacross local(i,j)
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
c$doacross local(i,j,sum)
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
c$doacross local(i,j,sum)
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
c$doacross local(i,j,sum)
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
c$doacross local(i,j,sum)
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
C          PRINT '(2I3,1P3G13.6)',I,J,EKB(I,J),EKS(I,J),EKH(I,J)
          EKB(I,J)=EKB(I,J)+EKS(I,J)+EKH(I,J)
        ENDDO
      ENDDO
C      PAUSE
C      IF(.FALSE.) THEN
C        DO I=1,12
C          PRINT 1003,(EKB(I,J)/EKB(I,I),J=1,12)
C        ENDDO
C        PRINT *,'---------------------------'
C        DO I=1,12
C          PRINT 1003,(EC(I,J)/EC(I,I),J=1,12)
C        ENDDO
C        PAUSE
C      ENDIF
1003  FORMAT(1X,1P6G13.6)
      RETURN
1000  FORMAT(1X,'BAD JACOBIAN',E10.3)
1001  FORMAT(A,/,(2I5,3X,1PE13.6))
1002  FORMAT(A,/,(1P2E13.6))
      END
C=============================================================
      SUBROUTINE OELEMGF(ETHICK,XY,N,EF,NL,XI,W)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XY(2,N)
      DIMENSION EF(12)
      DIMENSION DPSIX(4),DPSIY(4)
      DIMENSION PSI(4),DPSI(4,2),XS(2)
      DIMENSION XI(2,4),W(4)
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
C.....INITIALIZE ELEMENT ARRAYS
c$doacross local(i)
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
      SUBROUTINE OPASSMBGK(N3,NZ,GK0,GC0,KA,KZ,EK,EC,N,NODE)
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
c FIX FIX FIX in other versions...
c                 GC0(I,K)=GC0(I,K)+EK(L,M)
                  GC0(I,K)=GC0(I,K)+EC(L,M)
c FIX FIX FIX in other versions...
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
      SUBROUTINE OPASSMBGF(N3,GF,EF,N,NODE)
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
      SUBROUTINE OFORMGK(NMAX,N3,NZ,NUMEL,NNODE,
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
c$doacross local(i,j)
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
c$doacross local(L)
        DO L=1,12
          LLM(L)=KKX(NEL,L)
        ENDDO
        CALL OELEMEK(XY,4,EK,EC,4,XI,W)
C..................................
        IF(.false.) THEN
          write(7,*) NEL
          DO I=1,12
            write(7,*) (EK(I,J),J=1,12)
          ENDDO
c          PAUSE
        ENDIF
C..................................
        CALL OPASSMBGK(N3,NZ,GK0,GC0,KA,KZ,EK,EC,12,LLM)
C..................................
C        DO I=1,NNODE
C          PRINT 1000,(GK0(I,J),J=1,NNODE)
C        ENDDO
C        PAUSE
C..................................
      ENDDO
c$doacross local(i)
      DO I=1,NNODE
C        PRINT *,KZ(I),KA(I,NZ+1)
        KA(I,NZ+1)=KZ(I)
      ENDDO
C|||||||||||||||||||||||||||||||||||
C...................................
1000  FORMAT(1X,10F8.3)
      END
C=============================================================
      SUBROUTINE OFORMGF(NMAX,N3,NUMEL,NNODE,
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
c$doacross local(i)
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
c$doacross local(L)
        DO L=1,12
          LLM(L)=KKX(NEL,L)
        ENDDO
        ETHICK=RHOG*ETHICK/DBLE(4)
        CALL OELEMGF(ETHICK,XY,4,EF,4,XI,W)
C..................................
        IF(.false.) THEN
          write(7,*) NEL
          DO I=1,12
            write(7,*) EF(I)
          ENDDO
c          PAUSE
        ENDIF
C..................................
        CALL OPASSMBGF(N3,GF,EF,12,LLM)
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
      SUBROUTINE OTIMEDEP(N3,NZ,NNODE,GK0,GC0,GK,GF,KA,KZ,
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
c       call sbini(NNODE,NNODE,ia,ja,GK,iwork)
C ..... FORM MODIFIED STIFFNESS AND LOAD .......
        DELT1=1.D0/DELT    
        DO I = 1,NNODE
          DO J = 1,KZ(I)
            JG=KA(I,J)
C           PRINT *,I,J,GF(I),GC0(I,J)*WWW(JG)*DELT1
c            if(.false.) then ! full matrix...
              GF(I)=GF(I)+GC0(I,J)*WWW(JG)*DELT1
              GKIJ=GK0(I,J)+GC0(I,J)*DELT1
              call sbsij(NNODE,nit,ia,ja,GK,iwork,I,JG,GKIJ,
     &                   0,0,6,ier)
c            else ! lumped
c              GF(I)=GF(I)+GC0(I,J)*WWW(I)*DELT1
c              GKIJ=GK0(I,1)+GC0(I,J)*DELT1
c              call sbsij(NNODE,nit,ia,ja,GK,iwork,I,I,GKIJ,
c     &                   1,0,6,ier)
c            endif
          ENDDO 
        ENDDO
        call sbend(NNODE,nit,ia,ja,GK,iwork)
      else
        call sbini(NNODE,nit,ia,ja,GK,iwork)
c       call sbini(NNODE,NNODE,ia,ja,GK,iwork)
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
C=======================================
      SUBROUTINE OSETMAT
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
      logical ipass
      logical file39
      data ipass /.true./
      save ipass,e,rnu,thcrust,relax
      if(ipass) then
        ipass=.false.
        print *,'reading from input..defdep'
        inquire(file='fort.39',exist=file39)
        iok=-1
        if(file39) then
          read(39,*,IOSTAT=iok) e,rnu,thcrust,relax
        endif
        if(iok.ne.0) then
          print *,'not found, defaults used'
          e=1.D11
          rnu=0.5D0
          thcrust=110.D0
          relax=2000.D0                  
          rewind(39)
          write(39,*) e,rnu,thcrust,relax
          print *,'  youngs modulus: ', e
          print *,'  poissons ratio: ', rnu
          print *,' crust thickness: ', thcrust
          print *,' relaxation time: ', relax
        else
          print *,'defdep found, values used:'
          print *,'  youngs modulus: ', e
          print *,'  poissons ratio: ', rnu
          print *,' crust thickness: ', thcrust
          print *,' relaxation time: ', relax
        endif
      endif


c ... MATERIAL PROPERTIES AND THICKNESS OF CRUST
C ... YOUNGS MODULUS (NT/M**2)
c      E=5.D10
c      E=1.D11
C ... POISONS RATIO (DIMENSIONLESS)
c      RNU=0.5D0
C ... LAMES COEFFICIENTS
      RLAMBDA=RNU*E/(1.D0-RNU**2)
      RMU=0.5D0*E/(1.D0+RNU)
C      RMU=10000.D0
C      RLAMBDA=10000.D0
C ... THICKNESS OF CRUST (90-130 KM) (110,000 M)
c      TTT=110.D0*1000.D0
      TTT=thcrust*1000.D0
C ... RHOICE*GRAV (NT/M**3) TIMES THICKNESS, YIELDS PRESSURE, NT/M**2
      GRAV=9.8d0*0.3816d0
      RHOR=4000.D0
      RHOW=1092.D0
      RHOI=917.D0
      RHOG=RHOI*GRAV
      ROCKICE=RHOR/RHOI
C ... RELAXATION TIME CONSTANT ......!
c      RELAX=6000.D0                  !
c      RELAX=2000.D0                  !
      CTIME=3.16516D8*relax/6000.d0  !
      END
C---------------------------------------------
      SUBROUTINE WRITEDEPO(NUMNP,NNODE)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      parameter(nmax=maxnum,n3=3*nmax)
      CHARACTER*1 CHAR
      common /oelastic/ www(n3),wrate(n3,2),wwworig(nmax)
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
          WRITE(45,*) I,WWW(I),0.d0
        ENDDO
        DO I=1,NNODE
          WRITE(45,*) I,wrate(I,1),0.d0
        ENDDO
        DO I=1,NNODE
          WRITE(45,*) I,wrate(I,2),0.d0
        ENDDO
        DO I=1,NUMNP
          WRITE(45,*) I,WWWORIG(I),0.d0
        ENDDO
      ENDIF
      END
