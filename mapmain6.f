      PARAMETER( NMAX=29999)
      IMPLICIT REAL*8(A-H,O-Z)
C ***************************************************************C
C                                                                C
C   PROGRAM:  MAP6                                               C
C                                                                C
C   DIFFERS FROM MAP5 IN THAT IT USES A WAVEFRONT SOLVER.        C
C      UNFORTUNATELY IT IS MUCH SLOWER THAN THE ITERATIVE        C
C      SOLVER USED IN MAP5, SO WE WONT USE IT...                 C
C      THE ONE THING IT DOES IS SHOW THAT THE ITERATIVE          C
C      SOLVER IS OK, SINCE THIS DIRECT METHOD GIVE THE SAME      C
C      ANSWERS AS THE ITERATIVE SOLVER (TO 3-4 PLACES)           C
C                                                                C
C   DATE:  7-31-95                                               C
C   PROGRAMMER:  J.  FASTOOK                                     C
C                                                                C
C   FUNCTION:                                                    C
C           THIS IS A PROGRAM TO MODEL THE FLOW OF A GLACIER     C
C           WITH MATERIAL PROPERTIES READ FOR EACH NODAL POINT   C
C           AND THE AVERAGE USED FOR THE ELEMENT.                C
C           IT CALCULATES AN ELEMENT MATRIX ALA BECKER, 2-D      C
C           PROGRAM WITH LINEAR SHAPE FUNCTIONS.                 C
C           DOES TIME DEPENDENT CASE USING A LUMPED CAPACITANCE  C
C           MATRIX AND A BACKWARD DIFFERENCE SCHEME. ALLOWING    C
C           FOR VERY QUICK OPERATION.                            C
C ***************************************************************C
C     PROGRAM FOR STEADY AND UNSTEADY STATE FLOW ANALYSIS
C     USING FINITE ELEMENTS
C
C ***************************************************************C
      CHARACTER HED*80,SCRTCH*80,IADJ*2
      COMMON /LAPSE/ ACOM,HMAX,WINDIR(2),XPOLE,YPOLE,ABL
      DIMENSION IFIT(6), AMASS(11),  X(NMAX), Y(NMAX), 
     &          KODE(NMAX), CONST(NMAX), LM(5), ITYPE(NMAX),
     &          KX(NMAX,4), IDT(NMAX), ADOTB(NMAX), 
     &          ADOT(NMAX), BDROCK(NMAX), FLOWA(NMAX), SLDGB(NMAX),
     &          PSURF(NMAX), PPSURF(NMAX), FRACT(NMAX),
     &          CNEW(NMAX), QHOLD(NMAX), HTICE(NMAX), THICK(NMAX),
     &          HFIT(NMAX), IBFLUX(NMAX,2), BFLUX(NMAX), DEPB(NMAX),
     &          Q(NMAX),TEMP(NMAX),TBED(NMAX),
     &          UNDEPB(NMAX), SLOPE(4,NMAX), SLOPN(4,NMAX)
      DIMENSION AFUDGE(NMAX)
      DIMENSION TTIME(10000),VVOL(10000),AAREA(10000),TTBOT(10000),
     &          TTNSL(10000)
      DIMENSION NTYPE(NMAX), NNODE(NMAX),
     &          AADOT(NMAX), AFRACT(NMAX), AFLOWA(NMAX),
     &          ABDRCK(NMAX), ASLDGB(NMAX)
      REAL*8 T(NMAX),TNEW(NMAX)
      DATA HFIT /NMAX*0.0/
      CALL SETRIG
      ITOGG=0
      XPOLE=0.0D0
      YPOLE=0.0D0
      ACOM=-9.6237669
      WINDIR(1)=90.D0/180.*3.14159
      WINDIR(2)=100.D0
c     WINDIR(1)=315.D0/180.*3.14159
c     WINDIR(2)=50.D0
C FOLLOWING IS DEFAULT SNOWLINE ELEVATION AT POLE, GRADIENT, TNSL
      AMASS(1)=-5.
      AMASS(2)=1.
      AMASS(3)=.1
      AMASS(4)=500.
      AMASS(5)=1000.
      AMASS(6)=1500.
      AMASS(7)=-3000.
      AMASS(8)=.56E-3
      AMASS(8)=.001
      AMASS(9)=-14.
      AMASS(10)=.4
      AMASS(11)=1.2
C*******************
C ... CHECK FO DEFAULTS FILE. IF NONE EXIST, CREATE ONE ...
      READ(9,*,END=100) AMASS,ACOM,WJUNK,WINDIR(2),
     &                   XPOLE,YPOLE
        WINDIR(1)=WJUNK/180.*3.14159
        PRINT *,'**** DEFAULTS FILE FOUND ****'
        GOTO 101
100   PRINT *,'NO DEFAULTS FILE FOUND'
        REWIND 9
        WJUNK=WINDIR(1)*180./3.14159      
        WRITE(9,*) AMASS,ACOM,WJUNK,WINDIR(2),
     &             XPOLE,YPOLE
101   CONTINUE
      PRINT *,' TNSL       = ',AMASS(9)
      PRINT *,' FUDGE      = ',AMASS(10)
      PRINT *,' GEO GRAD   = ',AMASS(11)
      PRINT *,' ACOM       = ',ACOM
      PRINT *,' WINDIR     = ',WJUNK,WINDIR(2)
      PRINT *,' XPOLE,YPOLE= ',XPOLE,YPOLE
C*******************
      DO I=1,5
        IFIT(I)=0
      ENDDO
      IFIT(6)=1
      NTSTEP=0
C
C ... FOLLOWING IS SEALEVEL REFERENCED TO PRESENT=0.
      SEALEV=0.
C
C ... FOLLOWING SETS RATE OF CONVERGENCE, UP TO 5 WORKS WELL
      CONV=1.0
      TIME=0.0
      RHOW=1.092
      PG = 0.089866*0.3816
      NCOL=9
C
      DO I=1,4
        LM(I)=0
      ENDDO
C
C ... FOLLOWING READN FOR SPLIT DATA SETS
      CALL READN(NMAX, IDEP, HED, NUMNP, NUMEL, NUMGBC, NDT,
     &           INTER, DT, KODE, X, Y, HTICE, ADOT, ADOTB, FRACT,
     &           PSURF, RHOI, BDROCK, UNDEPB, FLOWA, SLDGB,TEMP,
     &           ITYPE,AFUDGE,
     &           THICK, KX, CONST, IBFLUX, BFLUX, QHOLD,
     &           NTYPE, NNODE, NCOL, AADOT, AFRACT,
     &           ABDRCK, PPSURF, AFLOWA, ASLDGB, IDT, AMASS,
     &           NUMCOL,NUMLEV)
      DO I=1,NUMNP
        Q(I)=HTICE(I)
        T(I)=HTICE(I)
      ENDDO
C
C ... SET LINEARIZATION CONSTANT USING INITIAL CONFIGURATION
      CALL NCONST(NMAX, IDEP, X, Y, KX, NTYPE, NUMEL,
     &            AFRACT, ASLDGB, LM, AFLOWA, BDROCK, DEPB,
     &            UNDEPB, PG, Q, CNEW, SLOPE, RHOI, WINDIR)
C
      WRITE(*,*) 'TIME STEP=',DT
      WRITE(*,*) 'INPUT 1 TO CALL ADJUST, 0 TO BYPASS'
      READ(*,4000) IADJ
      WRITE(99,4000) IADJ
      IF(IADJ.EQ.'1') THEN
C
C ..... CALCULATE SLOPES IN CASE NEEDED BY ADJUST
        CALL NODESL(NMAX, NUMNP, NUMEL, KX, SLOPE, SLOPN)
C
C ..... ENTER INTERACTIVE DATA SET MANIPULATOR
        CALL ADJUST(HED, NUMNP, NUMEL, X, Y, HTICE, ADOT, ADOTB,FRACT,
     &              TEMP,ITYPE,TBED,
     &         PSURF, BDROCK, DEPB, FLOWA, SLDGB, THICK, KX, CONST,
     &              AFUDGE,
     &              NNODE, KODE, HFIT, NUMCOL, NUMLEV, NUMGBC, NDT,
     &              INTER, DT, IBFLUX, BFLUX, NMAX, IDT, SLOPN, AMASS,
     &              TIME, NTSTEP,TTIME,VVOL,AAREA,TTBOT,TTNSL,IFIT,
     &              ITOGG,IPLOT)
C
      ENDIF
      DO I=1,NUMNP
        Q(I)=HTICE(I)
      ENDDO
C
C ... LOAD NODAL PROPERTIES INTO ELEMENT PROPERTIES
      CALL ELPROP(NMAX, NUMEL, NTYPE, KX, ADOT, AADOT, FRACT, AFRACT,
     &            BDROCK, ABDRCK, PSURF, PPSURF, FLOWA, AFLOWA,
     &            SLDGB, ASLDGB)
C
C
C
C
C
C .....................*** MAIN LOOP ***............................
C
C
C
215   CONTINUE
C
C 
C ..... FOLLOWING CALCULATES NODE SLOPE FROM ELEMENT SLOPE
        CALL NODESL(NMAX, NUMNP, NUMEL, KX, SLOPE, SLOPN)
C
C ..... FOLLOWING ADJUSTS ADOT FOR FITTED ACCUMULATION
        DO I=1,NUMNP
          IF(IDT(I).GT.0) THEN
            ADOT(I)=AFUNCT(TIME, IDT(I), AMASS,
     &                     Q(I), BDROCK(I),
     &                     SLOPN(1,I),X(I),Y(I),TEMP(I))
          ELSE
            IF(ITOGG.EQ.1) THEN
              AJUNK=ACCUM(X(I)*.001D0,Y(I)*.001D0,Q(I),SLOPN(1,I),
     &                    0.D0,AMASS(9),TEMP(I))
              ADOT(I)=ADOTB(I)-ABL*.01
            ELSE
              ADOT(I)=ADOTB(I)
            ENDIF
          ENDIF
        ENDDO
C
C       CALL TEMPER(NMAX, NUMNP, NTYPE, KX, KODE, ADOT, FRACT,
C     &              DEPB, HTICE, FLOWA, SLDGB, TEMP, ITYPE, TBED,
C     &              X, Y, AMASS, SLOPN, TIME, DT, IPLOT, TBAVG,
C     &              AFUDGE)
C ..... LOADS NEW NODAL MATERIAL PROPERTIES INTO ELEMENT MATERIAL PROPERTIES
        CALL ELPROP(NMAX, NUMEL, NTYPE, KX, ADOT, AADOT, FRACT, AFRACT,
     &            BDROCK, ABDRCK, PSURF, PPSURF, FLOWA, AFLOWA,
     &            SLDGB, ASLDGB)
C
C ..... CALCULATES VOLUMES (FLOTATION AND TOTAL) AND AREAL EXTENT
        CALL VOLUME(NMAX, TIME, NUMNP, NUMEL, X, Y, KX, Q, BDROCK,
     &             DEPB, ADOT, SEALEV, RHOI, RHOW, VOL, AREA, AMASS)
C
        NTSTEP=NTSTEP+1
        TTIME(NTSTEP)=TIME
        VVOL(NTSTEP)=VOL*1.E-15
        AAREA(NTSTEP)=AREA*1.E-12
        TTBOT(NTSTEP)=TBAVG
        TTNSL(NTSTEP)=AMASS(9)
C
C
        LL=0
        LF=0
        IF(IPLOT.GT.0 .AND. IPLOT.LT.6) THEN
          CALL GRSTRT(500,1)
          CALL WINDOW(0.,100.,0.,100.)
        ENDIF
        IDONE=0
C
C
C
C ..... ****************************************************************
C ..... LOOP ON NUMBER OF TIME STEPS *******************************
C ..... ****************************************************************
        DO 450 L=1,NDT
C
C ..... FORM STIFFNESS,CAPACITANCE AND LOAD
c        CALL FORMC(NMAX, X, Y, KX, NTYPE, NUMNP, NUMEL, 
c     &                 ETA, XI, W, CONST, T, KODE, 
c     &                 NUMGBC, IBFLUX, BFLUX, NZ, KZ, LM, AADOT, 
c     &                 D, B, A, KA)
C
C ......... CALCULATE EFFECTIVE LOAD  AND STIFFNESS MATRIX
c          DT2=1.0/DT
c          DO I=1,NUMNP
c            IF (KODE(I).EQ.0) THEN
c              IF (D(I).NE.0.) THEN
c                D(I)=DT2*D(I)
c                A(I,1)=A(I,1)+D(I)
c                B(I)=B(I)+D(I)*T(I)
c              ENDIF
c            ENDIF
c          ENDDO
C
          DO JK=1,NUMNP
            TNEW(JK)=T(JK)
          ENDDO
c
C ....... GAUSS-SEIDEL ITERATIVE SOLVER
c          CALL GAUSEID(NMAX,NZ,NUMNP,A,KA,B,TNEW)
c
c
          call formwave(numnp, DT, X, Y, KX, NUMNP, NUMEL, 
     &                 CONST, KODE, 
     &                 NUMGBC, IBFLUX, BFLUX, AADOT, 
     &                 TNEW)
c
          DO JK=1,NUMNP
            Q(JK)=TNEW(JK)
          ENDDO
c
          TIME=TIME+DT
          WRITE(*,*) 'TIME=',TIME
          LL=LL+1
          LF=LF+1
          HMAX=-1.E30
          DIFF=0.
          NDIFF=0
C ....... FOLLOWING CHECKS TO MAKE SURE THE SURFACE IS NOT BELOW 
C ....... THE BED OR THE FLOTATION LINE. ALSO MEASURES DIFFERENCE
C ....... BETWEEN PRESENT SOLUTION AND PRESENT SURFACE.
          DO JK=1,NUMNP
            IF(DEPB(JK).LT.0.) THEN
              FLOT=(1.-RHOW/RHOI)*DEPB(JK)
              THIK=Q(JK)-FLOT
              SURF=0.
            ELSE
              THIK=Q(JK)-UNDEPB(JK)
              SURF=UNDEPB(JK)
            ENDIF
            IF(THIK.LT.0.) Q(JK)=SURF
            IF(KODE(JK).NE.1) THEN
              DIFF=DIFF+(Q(JK)-PSURF(JK))**2
              NDIFF=NDIFF+1
            ENDIF
            IF(Q(JK).GT.HMAX) THEN
              HMAX=Q(JK)
              NNMAX=JK
            ENDIF
          ENDDO
C
C ....... OUTPUT TIME STEP INFO TO SCREEN
          WRITE(*,*) 'NMAX SURF=',HMAX,' AT NODE',NNMAX,' DIFF=',
     &              REAL(SQRT(DIFF/REAL(NDIFF)))
C
C ....... DERIVE SLOPE IN CASE NEEDED BY MASS BALANCE PARAMETERIZATION
          CALL NODESL(NMAX, NUMNP, NUMEL, KX, SLOPE, SLOPN)
C
C ....... MODIFY TEMPERATURE FIELD...
          CALL TEMPER(NMAX, NUMNP, NTYPE, KX, KODE, ADOT, FRACT,
     &            DEPB, HTICE, FLOWA, SLDGB, TEMP, ITYPE, TBED,
     &            X, Y, AMASS, SLOPN, TIME, DT, IPLOT,TBAVG,
     &            AFUDGE)
C ....... LOAD NEW NODE PROPERTIES INTO ELEMENT PROPERTIES
          CALL ELPROP(NMAX, NUMEL, NTYPE, KX, ADOT, AADOT, FRACT,
     &             AFRACT, BDROCK, ABDRCK, PSURF, PPSURF, FLOWA,
     &             AFLOWA, SLDGB, ASLDGB)
C
C ....... OBTAIN NEW LINEARIZATION CONSTANT FROM LATEST SOLUTION 
          CALL NCONST(NMAX, IDEP, X, Y, KX, NTYPE, NUMEL,
     &             AFRACT, ASLDGB, LM, AFLOWA, BDROCK, DEPB, UNDEPB,
     &             PG, Q, CNEW, SLOPE, RHOI, WINDIR)
C
C ....... UPDATE ACCUMULTION RATES FOR NEW COBFIGURATION
          DO I=1,NUMNP
            IF(IDT(I).GT.0) THEN
              ADOT(I)=AFUNCT(TIME, IDT(I), AMASS,
     &                       Q(I), BDROCK(I),
     &                       SLOPN(1,I),X(I),Y(I),TEMP(I))
            ELSE
              IF(ITOGG.EQ.1) THEN
                AJUNK=ACCUM(X(I)*.001D0,Y(I)*.001D0,Q(I),SLOPN(1,I),
     &                      0.D0,AMASS(9),TEMP(I))
                ADOT(I)=ADOTB(I)-ABL*.01
              ELSE
                ADOT(I)=ADOTB(I)
              ENDIF
            ENDIF
          ENDDO
C
C ....... CALCULATE VOLUMES (FLOTATION AND TOTAL) AND AREA
          CALL VOLUME(NMAX, TIME, NUMNP, NUMEL, X, Y, KX, Q, BDROCK,
     &             DEPB, ADOT, SEALEV, RHOI, RHOW, VOL, AREA, AMASS)
C
          NTSTEP=NTSTEP+1
          TTIME(NTSTEP)=TIME
          VVOL(NTSTEP)=VOL*1.E-15
          AAREA(NTSTEP)=AREA*1.E-12
          TTBOT(NTSTEP)=TBAVG
          TTNSL(NTSTEP)=AMASS(9)
          DO I=1,NUMNP
            T(I)=Q(I)
          ENDDO
C
C ....... LOAD NEW LINEARIZATION CONSTANT INTO OLD
          DO I=1,NUMEL
C           CONST(I)=CNEW(I)
            CONST(I) = (CONV*CONST(I) + CNEW(I))/(CONV+1)
          ENDDO
C
C ....... PRINT OUT SOLUTION FOR APPROPRIATE TIME STEPS
          IF(LL.GE.INTER) THEN
            WRITE(*,*)
            WRITE(*,*) 'DUMPING SOLUTION . . . . . . . .'
            WRITE(*,*)
            DO N=1,NUMNP
              HTICE(N)=Q(N)
            ENDDO
            WRITE(SCRTCH,*) 'TIME=',NINT(TIME)
            WRITE(34) SCRTCH
            WRITE(34) (HTICE(I),I=1,NUMNP)
            WRITE(34) (ADOT(I),I=1,NUMNP)
            WRITE(34) (DEPB(I),I=1,NUMNP)
            WRITE(34) (CONST(I),I=1,NUMEL)
            WRITE(34) (CONST(I),I=1,NUMEL)
            WRITE(35) (FRACT(I),I=1,NUMNP)
            WRITE(35) (FLOWA(I),I=1,NUMNP)
            WRITE(35) (SLDGB(I),I=1,NUMNP)
            WRITE(35) (AFUDGE(I),I=1,NUMNP)
            WRITE(36) SCRTCH
            WRITE(36) (TBED(I),I=1,NUMNP)
            LL=0
          ENDIF
          IF(LF.GE.IFIT(6)) THEN
            DO N=1,NUMNP
              HTICE(N)=Q(N)
            ENDDO
            LF=0
            CALL MODIFY(IFIT,NUMNP,KODE,HTICE,DEPB,PSURF,
     &                  ADOT,FRACT,FLOWA,SLDGB,BDROCK,AFUDGE,IDT)
            CALL UNLOAD(NUMNP,BDROCK,UNDEPB,PSURF)
          ENDIF
          IF(IPLOT.NE.0) THEN
            CALL PLOTSOL(NUMNP,X,Y,HTICE,DEPB,KODE,FRACT,
     &                   PSURF,
     &                   NTSTEP,TTIME,VVOL,
     &                   AAREA,TTBOT,TTNSL,
     &                   IPLOT,
     &                   NUMEL,KX,CONST)
            IF(IPLOT.EQ.2.OR.IPLOT.EQ.3.OR.IPLOT.EQ.4) THEN
              CALL DRAWOUT(IDONE)
            ENDIF
          ENDIF
450     CONTINUE
C
C ..... ****************************************************************
C ..... **** END OF TIME STEP SECTION **********************************
C ..... ****************************************************************
C
        IF(IPLOT.NE.0) CALL GRSTOP1
        HMAX=-1.E30
        DIFF=0.
        NDIFF=0
        DO N=1,NUMNP
          IF(KODE(N).NE.1) THEN
            DIFF=DIFF+(Q(N)-PSURF(N))**2
            NDIFF=NDIFF+1
          ENDIF
          HTICE(N)=Q(N)
          IF(HTICE(N).GT.HMAX) THEN
            HMAX=HTICE(N)
            NNMAX=N
          ENDIF
C
C ....... BEGIN FLOTATION PART
C
C ....... DEAL WITH BED ABOVE SEA LEVEL AND SURFACE BELOW BED
          IF(HTICE(N).LE.UNDEPB(N)) HTICE(N)=UNDEPB(N)
C
C ....... DEAL WITH BED BELOW SEA LEVEL 
C ........AND SURFACE BELOW FLOTATION HEIGHT
          IF(DEPB(N).LE.0.) THEN
            FLOT=(1.-RHOW/RHOI)*DEPB(N)
C           IF(HTICE(N).LE.FLOT) HTICE(N)=FLOT
            IF(HTICE(N).LE.FLOT) HTICE(N)=0.
            IF(BDROCK(N).LE.-9999.) HTICE(N)=0.
          ENDIF
C ....... END FLOTATION PART ****
        ENDDO
C
C ..... OUTPUT TIME STEP INFO TO SCREEN
        WRITE(*,*) 'NMAX SURF=', HMAX, ' AT NODE',NNMAX,' DIFF=',
     &              SQRT(DIFF/REAL(NDIFF))
C
C
C
        WRITE(*,*) ' END OF RUN'
C
C ..... END OF RUN, ENTER ADJUST TO GO ROUND AGAIN.
        WRITE(*,*) 'INPUT 1 TO CALL ADJUST -9 TO SKIP AND END'
        READ(*,4000,END=999) IADJ
        WRITE(99,4000) IADJ
        IF(IADJ.EQ.'-9') GOTO 999
C
C ..... ENTER INTERACTIVE DATA SET MANIPULATOR
        CALL ADJUST(HED, NUMNP, NUMEL, X, Y, HTICE, ADOT, ADOTB,FRACT, 
     &              TEMP, ITYPE, TBED, PSURF,
     &              BDROCK, DEPB, FLOWA, SLDGB, THICK, KX, CONST, 
     &              AFUDGE,NNODE, KODE,
     &              HFIT, NUMCOL, NUMLEV, NUMGBC, NDT, INTER, DT,
     &              IBFLUX, BFLUX, NMAX, IDT, SLOPN, AMASS, TIME,
     &              NTSTEP, TTIME, VVOL, AAREA,TTBOT, TTNSL, IFIT,
     &              ITOGG, IPLOT)
C
        DO I=1,NUMNP 
          Q(I)=HTICE(I)
        ENDDO
C
C ..... LOAD NODAL PROPERTIES INTO ELEMENT PROPERTIES
        CALL ELPROP(NMAX, NUMEL, NTYPE, KX, ADOT, AADOT, FRACT, AFRACT,
     &              BDROCK, ABDRCK, PSURF, PPSURF, FLOWA, AFLOWA, SLDGB,
     &              ASLDGB)
C
        WRITE(*,*) 'INPUT 1 TO CONTINUE WITH NEW SET, 0 TO STOP '
        READ(*,4000,END=999) IADJ
        WRITE(99,4000) IADJ
        IF(IADJ.EQ.'0') GOTO 999
C ..... GOTO BEGINNING OF MAIN LOOP
C
C
C
      GOTO 215
C
C
C
999   CONTINUE
C
C ... VERBOSE WRITER OFF
C     CALL WRITER(HED, NUMNP, NUMEL, X, Y, HTICE, ADOT, FRACT, PSURF,
C    &            RHOI, BDROCK, FLOWA, SLDGB, THICK, KX, CONST, NNODE,
C    &            KODE, HFIT, NUMCOL, NUMLEV, NUMGBC, NDT, INTER, DT,
C    &            IBFLUX, BFLUX, AADOT, AFRACT, AFLOWA, ABDRCK,
C    &            ASLDGB)
C
      WRITE(18,2000) -99999.,2.,0,HED
C
C ... FORMAT STATEMENTS
2000  FORMAT(10X,G13.6,2X,G13.6,I13,/,A80)
4000  FORMAT(A2)
C
      STOP
      END
c===========================
      subroutine formwave(nnode, DT, X, Y, KX, NUMNP, nelem, 
     &                 CONST, KODE, 
     &                 NUMGBC, IBFLUX, BFLUX, AADOT, 
     &                 zf)
      implicit real*8(a-h,o-z)
      parameter(mmax=29999)
      parameter(nwave=100,nmax=29999,nodes=4,big=1e10)
      dimension IBFLUX(NMAX,2), BFLUX(NMAX),AADOT(NMAX)
      dimension kode(nmax),eff(nmax),zf(nmax)
      dimension x(nmax),y(nmax),f(nmax)
      dimension ss(nodes,nodes),ff(nodes)
      integer KX(nodes,nmax),lm(nodes)
c ... wave front arrays
      dimension gfk(nwave,nwave),gff(nwave),gfeq(nmax)
      dimension ueq(nmax,nwave),diag(nmax)
      integer node(2,nmax),ifront(nwave)
      integer iord(nmax),ieq(nmax),jeq(nmax,nwave),keq(nmax)
      logical found(nmax)
      dimension CONST(nmax)
c .... load appropriate data into wavefront arrays...
      do i=1,nnode
        f(i)=zf(i)
      enddo
c ... solve with wave front solver ...
      do i=1,nwave
        gff(i)=0.
        do j=1,nwave
          gfk(i,j)=0.
        enddo
      enddo
c ... test to see if wavefront is big enough
c ... and generate occurance array
      call occur(nwave,nmax,nodes,nnode,nelem,node,KX,found)
      nfront=0
      neq=0
      nfixed=0
      do ielem=1,nelem
        call forms(nmax,nodes,x,y,f,KX,ielem,ss,ff,lm,DT,
     &             aadot(ielem),const(ielem),kode)
c        print *,'adding element:',ielem
c        print *
c        do i=1,nodes
c          print 222,(ss(i,j),j=1,nodes)
c        enddo
c        pause
222     format(1p4g13.6,i5,g13.6)
c
        do i=1,nodes
          if(node(1,lm(i)).eq.ielem) then
c ......... add to front ...
c            print *,'adding',lm(i),' to front'
            nfront=nfront+1
            ifront(nfront)=lm(i)
            iord(lm(i))=nfront
          endif
        enddo
c       print 495,'front:',(ifront(i),i=1,nfront)
495   format(1x,a,(1x,40i4))
c       print 500,(j,j=1,nnode)
c       print 500,(iord(j),j=1,nnode)
500   format(1x,'iord',(1x,25i3))
        do i=1,nodes
          lmi=iord(lm(i))
          gff(lmi)=gff(lmi)+ff(i)
          do j=1,nodes
            lmj=iord(lm(j))
            gfk(lmi,lmj)=gfk(lmi,lmj)+ss(i,j)
          enddo
        enddo
c        print *,'before pivot'
c        do j=1,nfront
c          print 200 (gfk(j,k),k=1,nfront),gff(j),zf(ifront(j))
c        enddo
200     format(1x,15f7.2)
        do k=1,nodes
          lmk=lm(k)
c ...................
          if(node(2,lmk).eq.ielem) then
c ......... remove from wave front ...
            do irm=1,nfront
              if(ifront(irm).eq.lmk) goto 1020
            enddo
            write(*,*) 'problems, free node never found ',
     &                 lmk,' in front'
            stop
1020        continue
            if(kode(lmk).eq.1) then
c ........... remove fixed nodes from wave front ...
c              print *,' removing fixed node ',lmk
              nfixed=nfixed+1
              ifixed=iord(lmk)
              do j=1,ifixed-1
                gff(j)=gff(j)-zf(lmk)*gfk(j,ifixed)
              enddo
              do j=ifixed+1,nfront
                gff(j)=gff(j)-zf(lmk)*gfk(j,ifixed)
              enddo
              iordlmk=iord(lmk)
              gff(iordlmk)=gff(nfront)
c ... zero out the things removed...
              gff(nfront)=0.0
c
              do kk=1,nfront
                gfk(iordlmk,kk)=gfk(nfront,kk)
c ... zero out the things removed...
                gfk(nfront,kk)=0.0
c
              enddo
              do kk=1,nfront-1
                gfk(kk,iordlmk)=gfk(kk,nfront)
c ... zero out the things removed...
                gfk(kk,nfront)=0.0
c
              enddo
              ifront(irm)=ifront(nfront)
              iord(lmk)=0
              if(irm.ne.nfront) iord(ifront(irm))=irm
              nfront=nfront-1
c              print *,'after removal of fixed'
c              print 500,(j,j=1,nnode)
c              print 500,(iord(j),j=1,nnode)
c              do jj=1,nfront
c                print 200, (gfk(jj,kk),kk=1,nfront),gff(jj),
c     &                        zf(ifront(jj))
c              enddo
            else
              neq=neq+1
              ieq(neq)=lmk
c              print *,' removing free node ',lmk
              
              ipivot=iord(lmk)
c              print *,' pivot =',ipivot,gfk(ipivot,ipivot)
              denom=1./gfk(ipivot,ipivot)
              diag(neq)=denom
              do j=1,ipivot-1
                pivot=denom*gfk(j,ipivot)
                gff(j)=gff(j)-pivot*gff(ipivot)
                do kk=1,nfront
                  gfk(j,kk)=gfk(j,kk)-pivot*gfk(ipivot,kk)
                enddo
              enddo
              do j=ipivot+1,nfront
                pivot=denom*gfk(j,ipivot)
                gff(j)=gff(j)-pivot*gff(ipivot)
                do kk=1,nfront
                  gfk(j,kk)=gfk(j,kk)-pivot*gfk(ipivot,kk)
                enddo
              enddo
c              print *,'after pivot'
c              do jj=1,nfront
c                print 200, (gfk(jj,kk),kk=1,nfront),gff(jj),
c     &                        zf(ifront(jj))
c              enddo
              iordlmk=iord(lmk)
              gfeq(neq)=gff(iordlmk)
              gff(iordlmk)=gff(nfront)
c ... zero out the things removed...
              gff(nfront)=0.0
c
              do kk=1,nfront
                jeq(neq,kk)=ifront(kk)
                ueq(neq,kk)=gfk(iordlmk,kk)
              enddo
c              do kk=iordlmk+1,nfront
c                jeq(neq,kk)=ifront(kk)
c                ueq(neq,kk)=gfk(iordlmk,kk)
c              enddo
              do kk=1,nfront
                gfk(iordlmk,kk)=gfk(nfront,kk)
c ... zero out the things removed...
                gfk(nfront,kk)=0.0
c
              enddo
              do kk=1,nfront-1
                gfk(kk,iordlmk)=gfk(kk,nfront)
c ... zero out the things removed...
                gfk(kk,nfront)=0.0
c
              enddo
              ifront(irm)=ifront(nfront)
              iord(lmk)=0
              if(irm.ne.nfront) iord(ifront(irm))=irm
              keq(neq)=nfront
              nfront=nfront-1
c              print *,'after reduction'
c              print 500,(j,j=1,nnode)
c              print 500,(iord(j),j=1,nnode)
c              do jj=1,nfront
c                print 200, (gfk(jj,kk),kk=1,nfront),gff(jj),
c     &                        zf(ifront(jj))
c              enddo
            endif
          endif
        enddo
c        pause
      enddo
c      print *, 'wave front matrix:'
      do j=neq,1,-1
c        print * 
c        print *, ' diag=',1./diag(j)
c        print *, ' ieq =',ieq(j)
c        print *, ' keq =',keq(j)
c        print *, ' gfeq=',gfeq(j)
c        print 300, (jeq(j,k),k=1,keq(j))
c        print 301, (ueq(j,k),k=1,keq(j))
300     format(1x,9i10)
301     format(1x,9f10.4)
        sum=0.0
        do i=1,keq(j)
          jeqji=jeq(j,i)
          if(jeqji.ne.ieq(j)) then
            sum=sum-ueq(j,i)*zf(jeq(j,i))
           endif
        enddo
        zf(ieq(j))=(gfeq(j)+sum)*diag(j)
c        print *,j,ieq(j),zf(ieq(j))
c
      enddo
c
      end
c *********************************************************
      subroutine occur(nwave,nmax,nodes,nnode,nelem,
     &                 node,KX,found)
      implicit real*8(a-h,o-z)
      integer KX(nmax,nodes),node(2,nmax)
      logical found(nmax)
      do i=1,nnode
        found(i)=.false.
      enddo
      do j=1,nelem
        do k=1,nodes
          if(.not.found(KX(j,k))) then
            node(1,KX(j,k))=j
            found(KX(j,k))=.true.
          endif
          node(2,KX(j,k))=j
        enddo
      enddo
c      write(*,*) '   node','  first','   last'
      iwmax=0
      do i=1,nnode
        iwmax=max(iwmax,node(2,i)-node(1,i))
c        write(*,100) i,(node(j,i),j=1,2)
      enddo
      write(*,*) 'minimum wavefront requirement:',iwmax+2
      if(iwmax+2.gt.nwave) then
        write(*,*) 'wavefront not big enough:',nwave,' -->',iwmax
        stop
      endif
100   format(1x,3i7)
      end
c *********************************************************
      subroutine forms(nmax,nodes,x,y,f,KX,ielem,ss,ff,lm,DT,
     &             aadot,const,kode)
      implicit real*8(a-h,o-z)
      parameter(lnode=4,nint=9)
      dimension x(nmax),y(nmax),f(nmax),kode(nmax)
      integer KX(nmax,nodes),lm(nodes)
      dimension ss(nodes,nodes),ff(nodes),DD(lnode)
      dimension dpsix(lnode),dpsiy(lnode),dxds(2,2),dsdx(2,2)
      dimension psi(lnode),dpsi(lnode,2)
      dimension xy(2,lnode),xi(nint),eta(nint),w(nint)
      save ipass,xi,eta,w
      data ipass /0/
c      write(*,*) 'in forms'
      if(nodes.gt.lnode) then
        write(*,*) 'problems in forms, nodes=',nodes
        stop
      endif
      if(ipass.eq.0) then
        call gausinit(xi,eta,w)
        ipass=1
      endif
      do i=1,nodes
        ff(i)=0.0
        do j=1,nodes
          ss(i,j)=0.0
        enddo
      enddo
      call loadlm(nmax,nodes,ielem,lm,KX)
      do i=1,nodes
        xy(1,i)=x(lm(i))
        xy(2,i)=y(lm(i))
      enddo
c .......................................
c ... form element matrix and vectors ...
c .......................................
c
c ... begin integration point loop ...
      do lint=1,nint
        call shape(1,xi(lint),eta(lint),psi,dpsi)
c
c ..... calculate dxds...equation (5.3.6) ...
        do i=1,2
          do j=1,2
            dxds(i,j)=0.0
            do k=1,nodes
              dxds(i,j)=dxds(i,j)+dpsi(k,j)*xy(i,k)
            enddo
          enddo
        enddo
c
c ..... calculate dsdx...equation (5.2.7) ...
        detj=(dxds(1,1)*dxds(2,2)-dxds(1,2)*dxds(2,1))
        if(detj.le.0.0) then
          write(*,1100) detj,((xy(mm,nn),nn=1,nodes),mm=1,2)
1100      format(' bad jacobian at 161',g13.6,/,4g13.6,/,4g13.6)
          stop
        endif
        denom=1./detj
        dsdx(1,1)=dxds(2,2)*denom
        dsdx(2,2)=dxds(1,1)*denom
        dsdx(1,2)=-dxds(1,2)*denom
        dsdx(2,1)=-dxds(2,1)*denom
c
c ..... calculate d(psi)/dx...equation (5.3.5) ...
        do i=1,nodes
          dpsix(i)=dpsi(i,1)*dsdx(1,1)+dpsi(i,2)*dsdx(2,1)
          dpsiy(i)=dpsi(i,1)*dsdx(1,2)+dpsi(i,2)*dsdx(2,2)
        enddo
c
c ..... accumulate integration point values of integrals ...
        fac=detj*w(lint)
        do i=1,nodes
c ....... lumped capacitance 
          DD(I)=DD(I)+PSI(I)*FAC
c ....... right hand side
          ff(i)=ff(i)+AADOT*psi(i)*fac
c ....... stiffness matrix
          do j=1,nodes
            ss(i,j)=ss(i,j)+
     &              fac*CONST*(dpsix(i)*dpsix(j)+
     &                               dpsiy(i)*dpsiy(j))
          enddo
        enddo
      enddo
      DT2=1.0/DT
      DO I=1,nodes
        IF (KODE(lm(I)).EQ.0) THEN
          IF (DD(I).NE.0.) THEN
            DD(I)=DT2*DD(I)
            SS(I,I)=SS(I,I)+DD(I)
            FF(I)=FF(I)+DD(I)*f(lm(I))
          ENDIF
        ENDIF
      ENDDO
      end
c *****************************************************************
c *****************************************************************
      subroutine gausinit(xi,eta,w)
      implicit real*8(a-h,o-z)
      dimension xi(9), eta(9), w(9)
c
c ... gaussian quadrature of order three quadrilaterals ...
      xi(1)=-sqrt(3.d0/5.d0)
      xi(2)=0.
      xi(3)=-xi(1)
      xi(4)=xi(1)
      xi(5)=0.
      xi(6)=xi(3)
      xi(7)=xi(1)
      xi(8)=0.
      xi(9)=xi(3)
c ...
      eta(1)=xi(1)
      eta(2)=xi(1)
      eta(3)=xi(1)
      eta(4)=0.
      eta(5)=0.
      eta(6)=0.
      eta(7)=xi(3)
      eta(8)=xi(3)
      eta(9)=xi(3)
c ...
      w(1)=25.d0/81.d0
      w(2)=40.d0/81.d0
      w(3)=w(1)
      w(4)=w(2)
      w(5)=64.d0/81.d0
      w(6)=w(2)
      w(7)=w(1)
      w(8)=w(2)
      w(9)=w(1)
c ...
      end
c *********************************************************
      subroutine loadlm(nmax,nodes,n,lm,KX)
      implicit real*8(a-h,o-z)
      integer lm(nodes),KX(nmax,nodes)
      do i=1,nodes
        lm(i)=KX(n,i)
      enddo
      end
