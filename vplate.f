      SUBROUTINE VPLATE(ITIME,NMAX,N3,NUMNP,NUMEL,X,Y,KX,THICK,
     &                 DELT,WWW,WRATE,WMIN,TIME,WWWORIG)
c-----------------------------------------------------------------------
c 4th order viscous plate solver with one-time generation of stiffness
c matrix. uses my sparse matrix storage for the static matrice
c and ITPACK sprse storage for the time-dependent modified matrices, and c             ipage=ipage+2

c iterative matrix solver. j fastook july, 1999
c-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)                                          
#include "parameter.h"
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
      DIMENSION EK(12,12),EF(12)
      DIMENSION GF(M3)
      DIMENSION GK0(M3,NZ)
      DIMENSION KA(M3,NZ1),KZ(M3)
      DIMENSION WWWTMP(M3)
      CHARACTER*80 HED
      COMMON /PMATPROP/ RMU,RLAMBDA,TTT,RHOG,CTIME,ROCKICE
      LOGICAL CTOGG,WTOGG,ITOGG,IOTOGG
      INTEGER BTOGG
      COMMON /TOGGLES/ CTOGG,WTOGG,ITOGG,BTOGG,IOTOGG,IPAGE
      parameter(npage=39)
      character*80 list(npage)
      COMMON /IOLIST/ LIST

C ... TIMER STUFF, FOR SGI ONLY ...!
      REAL*4 TB(2)
c     REAL*4 TB(2),ETIME,DTIME     !
c     EXTERNAL ETIME,DTIME         !
C .................................!
      SAVE IPASS,KKX,WSAVE,NNODE,XI,W,GK0,KA,KZ
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
c$doacross local(i)
c           DO I=1,NUMNP
c             WWW(I)=0.D0
c             WRATE(I,1)=0.D0
c             WRATE(I,2)=0.D0
c           ENDDO
          WRITE(HED,*) ' TIME=',NINT(TIME-DELT)
          WRITE(88) HED
          WRITE(88) (WWW(I),I=1,NNODE,3)
          WRITE(88) (THICK(I),I=1,NUMNP)
          WRITE(88) (1000*WRATE(I,1),I=1,NNODE,3)
        ENDIF
        CALL VFORMGK(NMAX,N3,NZ,NUMEL,NNODE,
     &              X,Y,KX,KKX,EK,XI,W,
     &              GK0,KA,KZ)
        CALL VLOADGK(N3,NZ,NNODE,GK0,GK,KA,KZ,
     &                    nit,nmax1,ia,ja,iwork)
      ENDIF
      DO I=1,NNODE
        WWWTMP(I)=WWW(I)
      ENDDO
C ... 2-STEP METHOD ...
      DO II=1,2
C ......................................................
C ......FORM STIFFNESS AND LOAD ..............
        CALL VFORMGF(NMAX,N3,NUMEL,NNODE,
     &              X,Y,THICK,WWWTMP,KX,KKX,EF,XI,W,
     &              GF)
C ......................................................
C.......DUMP MATRIX FOR EXAMINATION ......................
        CALL ODUMPMAT(N3,NZ,NNODE,KZ,KA,GK0,GF,WWWTMP,.false.)
        CALL NDUMPMAT(NIT,GK,JA,NMAX1,IA,NNODE,GF,.false.)
C.......................................................
C       IF(.TRUE.) THEN
C.........SOLVE EQUATIONS WITH JORDAN CONJUGATE-GRADIENT ITPACK ....!
c         WRITE(7,1000) ' TIME BEFORE JCG ',ETIME(TB),DTIME(TB)!
          call dfault(iparm,rparm)
          if(delt.ne.0.0) then
            iparm(1)=1000  ! MAX number of iteration
          else
            iparm(1)=1000 ! MAX number of iteration
          endif
          rparm(1)=1d-3   ! ZETA, stopping criteria
c         iparm(2)=2      ! LEVEL of output (-1:none)
c         iparm(4)=7      ! OUTPUT unit number
          iparm(5)=1      ! NONsymmetric matrix (1)
          iparm(6)=0      ! NON-ADAPTIVE (0) (1 doesnt work)
c         iparm(10)=1     ! REMOVES large diagonal entries (doesnt work)
C         DO I=1,12
C           PRINT *,I,IPARM(I),RPARM(I)
C         ENDDO
          do i=1,nnode
            wrate(i,II)=0.d0
          enddo
          CALL JCG(NNODE,ia,ja,GK,GF,WRATE(1,II),
     &            iwksp,nw,wksp,iparm,rparm,ier)

C ...................................................
c ....... EXPERIMENTAL...............................
c ....... REMOVE RIGID BODY MOTION ..................
C ....... THIS COULD BE DONE BETTER BY IMPOSING BC ..
C ....... ALONG THE EDGE OF THE GRID ................
          WRATEBASE=WRATE(1,II)                     !
          do i=1,nnode                              !
            wrate(i,II)=WRATE(I,II)-WRATEBASE       !
          enddo                                     !
c ....... EXPERIMENTAL, REMOVE RIGID BODY MOTION ....
C ...................................................

C         DO I=1,12
C           PRINT *,I,IPARM(I),RPARM(I)
C         ENDDO
           IF(IOTOGG) THEN
             write(list(ipage+1),*) ' relative error=',rparm(1),
     &              ' in iterations = ',iparm(1)
c             write(list(ipage+2),*) rparm(11),rparm(12)
c             ipage=ipage+2
             ipage=ipage+1
           ENDIF
          if(ier.ne.0) then
            CALL NDUMPMAT(NIT,GK,JA,NMAX1,IA,NNODE,GF,.TRUE.)
            print *,'JCG error:',ier
            print '(1x,i3,i10,1pg13.6)',(i,iparm(i),rparm(i),i=1,12)
            pause
          endif
c         WRITE(7,1000) ' TIME AFTER JCG ',ETIME(TB),DTIME(TB) !
C ................................................................!
C       ELSE
C.........SOLVE EQUATIONS WITH CONJUGATE-GRADIENT ..................!
C         WRITE(7,1000) ' TIME BEFORE CONJUG ',ETIME(TB),DTIME(TB)  !
C         CALL CONJUG(N3,NZ,NNODE,1.D-6,GK,KA,GF,WWW)              !
C         WRITE(7,1000) ' TIME AFTER CONJUG ',ETIME(TB),DTIME(TB)   !
C ................................................................!
C       ENDIF
        IF(II.EQ.1) THEN
          DO I=1,NNODE
            WWWTMP(I)=WWW(I)+WRATE(I,II)*DELT
          ENDDO
        ENDIF
      ENDDO
      DO I=1,NNODE
        WWW(I)=WWW(I)+0.5D0*(WRATE(I,1)+WRATE(I,2))*DELT
      ENDDO
      WMIN=1D30
      WMAX=-1D30
      WTOT=0.D0
      DO I=1,NNODE,3
        WMIN=MIN(WMIN,WWW(I))
        WMAX=MAX(WMAX,WWW(I))
        WTOT=WTOT+WWW(I)
C        WRITE(*,*) (real(WWW(I+J)),J=0,2)
c        WRITE(*,*) (real(wrate(I+J,1)),J=0,2)
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
      SUBROUTINE VELEMEK(XY,N,EKB,NL,XI,W)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XY(2,N),EKB(12,12),EKS(12,12)
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
        ENDDO
      ENDDO
C ... FORM D-MATRICES ...............
      CTIME1=CTIME*10
      TTT3=TTT**3/12.D0
      DB(1,1)=TTT3*(2.D0*RMU+RLAMBDA)*CTIME1
      DB(1,2)=TTT3*RLAMBDA*CTIME1
      DB(1,3)=0
      DB(2,1)=TTT3*RLAMBDA*CTIME1
      DB(2,2)=TTT3*(2.D0*RMU+RLAMBDA)*CTIME1
      DB(2,3)=0
      DB(3,1)=0
      DB(3,2)=0
      DB(3,3)=TTT3*RMU*CTIME1
      DS(1,1)=TTT*RMU*CTIME1
      DS(1,2)=0
      DS(2,1)=0
      DS(2,2)=TTT*RMU*CTIME1
C.....BEGIN 2X2 INTEGRATION LOOP
      DO L=1,NL
        XS(1)=XI(1,L)
        XS(2)=XI(2,L)
        CALL PGENSHAPE(N,XS,XY,PSI,DPSI,DETJ,DPSIX,DPSIY)
        CALL PLOADB(PSI,DPSIX,DPSIY,BB,BS)
C
C.......ACCUMULATE INTEGRATION POINT VALUE OF INTEGRALS
        FAC=DETJ*W(L)
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
C          PRINT '(2I3,1P3G13.6)',I,J,EKB(I,J),EKS(I,J)
          EKB(I,J)=EKB(I,J)+EKS(I,J)
        ENDDO
      ENDDO
C      PAUSE
C      IF(.FALSE.) THEN
C        DO I=1,12
C          PRINT 1003,(EKB(I,J)/EKB(I,I),J=1,12)
C        ENDDO
C      ENDIF
1003  FORMAT(1X,1P6G13.6)
      RETURN
1000  FORMAT(1X,'BAD JACOBIAN',E10.3)
1001  FORMAT(A,/,(2I5,3X,1PE13.6))
1002  FORMAT(A,/,(1P2E13.6))
      END
C=============================================================
      SUBROUTINE VASSMBGK(N3,NZ,GK0,KA,KZ,EK,N,NODE)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GK0(N3,NZ)
      DIMENSION KA(N3,NZ+1),KZ(N3)
      DIMENSION EK(12,12),NODE(N)
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
              KA(I,1)=I
            ELSE
              DO K=2,KZ(I)
                IF(KA(I,K).EQ.J) THEN
                  GK0(I,K)=GK0(I,K)+EK(L,M)
                  GOTO 99
                ENDIF
              ENDDO
              KZ(I)=KZ(I)+1
              GK0(I,KZ(I))=EK(L,M)
              KA(I,KZ(I))=J
            ENDIF
99          CONTINUE
        ENDDO
      ENDDO
C||||||||||||||||||||||||
C........................
      END
C=============================================================
      SUBROUTINE VFORMGK(NMAX,N3,NZ,NUMEL,NNODE,
     &                  X,Y,KX,KKX,EK,XI,W,
     &                  GK0,KA,KZ)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NMAX),Y(NMAX),KX(NMAX,4),KKX(NMAX,12)
      DIMENSION GK0(N3,NZ)
      DIMENSION KA(N3,NZ+1),KZ(N3)
      DIMENSION XI(2,4),W(4)
      DIMENSION EK(12,12)
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
        CALL VELEMEK(XY,4,EK,4,XI,W)
C..................................
        IF(.false.) THEN
          write(7,*) NEL
          DO I=1,12
            write(7,*) (EK(I,J),J=1,12)
          ENDDO
c          PAUSE
        ENDIF
C..................................
        CALL VASSMBGK(N3,NZ,GK0,KA,KZ,EK,12,LLM)
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
      SUBROUTINE VFORMGF(NMAX,N3,NUMEL,NNODE,
     &                  X,Y,THICK,WWW,KX,KKX,EF,XI,W,
     &                  GF)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NMAX),Y(NMAX),KX(NMAX,4),KKX(NMAX,12)
      DIMENSION GF(N3),WWW(N3)
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
        EWWW=0.D0
        DO L=1,4
          LM(L)=KX(NEL,L)
          XY(1,L)=X(LM(L))
          XY(2,L)=Y(LM(L))
          ETHICK=ETHICK+THICK(LM(L))
          EWWW=EWWW+WWW((LM(L)-1)*3+1)
        ENDDO
c$doacross local(L)
        DO L=1,12
          LLM(L)=KKX(NEL,L)
        ENDDO
        ETHICK=ETHICK/DBLE(4)
        EWWW=EWWW/DBLE(4)
C        print *,nel,real(ethick),real(ewww)
        ETHICK=RHOG*ETHICK
        EWWW=ROCKICE*RHOG*EWWW
C        print *,nel,real(ethick),real(ewww),real(ethick+ewww)
        ETHICK=ETHICK+EWWW
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
      SUBROUTINE VLOADGK(N3,NZ,NNODE,GK0,GK,KA,KZ,
     &                    nit,nmax1,ia,ja,iwork)
      IMPLICIT REAL*8(A-H,O-Z)
c ....arrays for ITPACK sparse storage...................
      dimension GK(nit),ja(nit),ia(nmax1)
      dimension iwork(nit)
      DIMENSION GK0(N3,NZ)
      DIMENSION KA(N3,NZ+1),KZ(N3)
      call sbini(NNODE,nit,ia,ja,GK,iwork)
      DO I = 1,NNODE
        DO J = 1,KZ(I)
          JG=KA(I,J)
          GKIJ=GK0(I,J)
          call sbsij(NNODE,nit,ia,ja,GK,iwork,I,JG,GKIJ,0,1,6,ier)
        ENDDO 
      ENDDO
      call sbend(NNODE,nit,ia,ja,GK,iwork)
      END
