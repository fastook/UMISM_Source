      program pack3
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      REAL*8 X(NMAX),Y(NMAX)
      integer KX(NMAX,4)
      include 'netcdf.inc'
* error status return
      integer iret
* netCDF id
      integer ncid
* dimension ids
      integer rowdim,coldim,time_dim
* dimension lengths
      integer rowlen,collen,time_len
      parameter(nrowmax=300,ncolmax=300)
      parameter(time_len=NF_UNLIMITED)
* variable ids
      integer rowid,colid,time_id
      integer Z1_id
      integer Z2_id
      integer Z3_id
      integer Z4_id
      integer Z5_id
      integer Z6_id
      integer Z7_id
      integer Z8_id
      integer Z9_id
* rank (number of dimensions) for each variable
      integer rowrank,colrank,time_rank,Z_rank,Z_nr
      parameter(rowrank=1,colrank=1,time_rank=1,Z_rank=3,Z_nr=1)
* variable shapes
      integer rowdims(rowrank),coldims(colrank),time_dims(time_rank)
      integer  Z_dims(Z_rank)
* data variables
      integer  row(nrowmax),col(ncolmax)
      real Z1(ncolmax, nrowmax, Z_nr)
      real Z2(ncolmax, nrowmax, Z_nr)
      real Z3(ncolmax, nrowmax, Z_nr)
      real Z4(ncolmax, nrowmax, Z_nr),Z4range(2)
      real Z5(ncolmax, nrowmax, Z_nr)
      real Z6(ncolmax, nrowmax, Z_nr)
      real Z7(ncolmax, nrowmax, Z_nr)
      real Z8(ncolmax, nrowmax, Z_nr)
      real Z9(ncolmax, nrowmax, Z_nr),Z9range(2)
* read header file first
      collen=ncolmax
      rowlen=nrowmax
      call readfixed(collen,rowlen,
     &     NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,
     &     KX,X,Y)
* enter define mode
      iret = nf_create('netcdf.nc', NF_CLOBBER, ncid)
      call check_err(iret)
* define dimensions
      iret = nf_def_dim(ncid, 'row', rowlen, rowdim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'col', collen, coldim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'time', time_len, time_dim)
      call check_err(iret)
* define variables
      rowdims(1) = rowdim
      iret = nf_def_var(ncid, 'row', NF_INT, rowrank, rowdims, rowid)
      call check_err(iret)
      coldims(1) = coldim
      iret = nf_def_var(ncid, 'col', NF_INT, colrank, coldims, colid)
      call check_err(iret)
      time_dims(1) = time_dim
      iret = nf_def_var(ncid, 'time', NF_INT, time_rank, time_dims, time
     1_id)
      call check_err(iret)
      Z_dims(3) = time_dim
      Z_dims(2) = rowdim
      Z_dims(1) = coldim
      iret = nf_def_var(ncid,'htice',NF_REAL, Z_rank, Z_dims, Z1_id)
      iret = nf_def_var(ncid,'bedrock',NF_REAL, Z_rank, Z_dims, Z2_id)
      iret = nf_def_var(ncid,'thick',NF_REAL, Z_rank, Z_dims, Z3_id)
      iret = nf_def_var(ncid,'adot',NF_REAL, Z_rank, Z_dims, Z4_id)
      iret = nf_def_var(ncid,'flowa',NF_REAL, Z_rank, Z_dims, Z5_id)
      iret = nf_def_var(ncid,'tbed',NF_REAL, Z_rank, Z_dims, Z6_id)
      iret = nf_def_var(ncid,'bmelt',NF_REAL, Z_rank, Z_dims, Z7_id)
      iret = nf_def_var(ncid,'wthick',NF_REAL, Z_rank, Z_dims, Z8_id)
      iret = nf_def_var(ncid,'velo',NF_REAL, Z_rank, Z_dims, Z9_id)
      call check_err(iret)
      Z4range(1) = 0
      Z4range(2) = 5
      iret = nf_put_att_real(ncid, Z4_id, 'valid_range', NF_REAL, 2, 
     1       Z4range)
      call check_err(iret)
      Z9range(1) = 0
      Z9range(2) = 5000
      iret = nf_put_att_real(ncid, Z9_id, 'valid_range', NF_REAL, 2, 
     1       Z9range)
      call check_err(iret)

* leave define mode
      iret = nf_enddef(ncid)
      call check_err(iret)
* store row
      do i=1,rowlen
        row(i)=i
      enddo
c     data row /0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90/
      iret = nf_put_var_int(ncid, rowid, row)
      call check_err(iret)
* store col
      do i=1,collen
        col(i)=i
      enddo
c     data col /0, 10, 20, 30, 40, 50, 60, 70, 80, 90/
      iret = nf_put_var_int(ncid, colid, col)
      call check_err(iret)
       
* Write record variables
      call writerecs(ncid,time_id,
     &     Z1_id,Z2_id,Z3_id,Z4_id,Z5_id,Z6_id,Z7_id,Z8_id,Z9_id,
     &     z1,z2,z3,z4,z5,z6,z7,z8,
     &     collen,rowlen,numnp,numel,
     &     KX,X,Y)
       
      iret = nf_close(ncid)
      call check_err(iret)
      end
********************************************************
      subroutine writerecs(ncid,time_id,
     &     Z1_id,Z2_id,Z3_id,Z4_id,Z5_id,Z6_id,Z7_id,Z8_id,Z9_id,
     &     z1,z2,z3,z4,z5,z6,z7,z8,
     &     collen,rowlen,numnp,numel,
     &     KX,X,Y)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      REAL*8 X(NMAX),Y(NMAX)
      integer KX(NMAX,4)
      include 'netcdf.inc'
* netCDF id
      integer ncid
* variable ids
      integer time_id
      integer Z1_id,Z2_id,Z3_id,Z4_id,Z5_id,Z6_id,Z7_id,Z8_id,Z9_id
* error status return
      integer iret
* netCDF dimension sizes for dimensions used with record variables
      integer rowlen,collen
c     parameter(rowlen=74,collen=76,ntotal=rowlen*collen)
* rank (number of dimensions) for each variable
      integer time_rank,Z_rank
      parameter(time_rank=1,Z_rank=3)
* starts and counts for array sections of record variables
      integer  time_start(time_rank), time_count(time_rank)
      integer  Z_start(Z_rank), Z_count(Z_rank)
* data variables
      integer  time_nr
      parameter (time_nr = 1)
      integer  time(time_nr)
      integer  Z_nr
      parameter (Z_nr = 1)
      real Z1(collen, rowlen, Z_nr)
      real Z2(collen, rowlen, Z_nr)
      real Z3(collen, rowlen, Z_nr)
      real Z4(collen, rowlen, Z_nr)
      real Z5(collen, rowlen, Z_nr)
      real Z6(collen, rowlen, Z_nr)
      real Z7(collen, rowlen, Z_nr)
      real Z8(collen, rowlen, Z_nr)
      real Z9(collen, rowlen, Z_nr)
c     data time /1 * NF_FILL_INT/
c     data Z /ntotal * NF_FILL_REAL/
c ..........................................
      ntime=0
* skip first one ...
      call readrest(collen,rowlen,
     &     NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,
     &     Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,
     &     itime,iend,
     &     KX,X,Y)
* skip first one ...
      call readrest(collen,rowlen,
     &     NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,
     &     Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,
     &     itime,iend,
     &     KX,X,Y)
      dowhile(iend.eq.0)
        ntime=ntime+1
        do j=1,rowlen
          do i=1,collen
            if(z2(i,j,1).eq.-9999.) z2(i,j,1)=NF_FILL_REAL
          enddo
        enddo
        time(1)=itime
        time_start(1)=ntime
        time_count(1)=1
        iret=nf_put_vara_int(ncid,time_id,time_start,time_count,time)
        call check_err(iret)
        Z_start(1)=1
        Z_start(2)=1
        Z_start(3)=ntime
        Z_count(1)=collen
        Z_count(2)=rowlen
        Z_count(3)=1
        iret=nf_put_vara_real(ncid,Z1_id,Z_start,Z_count,Z1)
        iret=nf_put_vara_real(ncid,Z2_id,Z_start,Z_count,Z2)
        iret=nf_put_vara_real(ncid,Z3_id,Z_start,Z_count,Z3)
        iret=nf_put_vara_real(ncid,Z4_id,Z_start,Z_count,Z4)
        iret=nf_put_vara_real(ncid,Z5_id,Z_start,Z_count,Z5)
        iret=nf_put_vara_real(ncid,Z6_id,Z_start,Z_count,Z6)
        iret=nf_put_vara_real(ncid,Z7_id,Z_start,Z_count,Z7)
        iret=nf_put_vara_real(ncid,Z8_id,Z_start,Z_count,Z8)
        Z_start(1)=1
        Z_start(2)=1
        Z_start(3)=ntime
        Z_count(1)=collen-1
        Z_count(2)=rowlen-1
        Z_count(3)=1
        iret=nf_put_vara_real(ncid,Z9_id,Z_start,Z_count,Z9)
        call readrest(collen,rowlen,
     &       NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,
     &       Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,
     &       itime,iend,
     &       KX,X,Y)
      enddo
       
      end
       
      subroutine check_err(iret)
      integer iret
      include 'netcdf.inc'
      if (iret .ne. NF_NOERR) then
      write(*,*) nf_strerror(iret)
      stop
      endif
      end


      subroutine readfixed(collen,rowlen,
     &           NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,
     &           KX,X,Y)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      CHARACTER HED*80
      integer collen,rowlen
      REAL*8 X(NMAX),Y(NMAX),psurf,bflux
      integer KX(NMAX,4)
C READ INPUT HEADER
      READ(30,1000,end=999) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,
     &                      NDT,INTER,DT
      write(*,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,
     &                      NDT,INTER,DT
1000  FORMAT (A80,/,7I6,F8.0)
      if(numcol.gt.collen .or. numlev.gt.rowlen) then
        write(*,*) 'recompile with:'
        write(*,*) '  collen =',numcol
        write(*,*) '  rowlen =',numlev
        stop
      else
        collen=numcol
        rowlen=numlev
      endif
C READ INPUT GRID, THINGS THAT NEVER CHANGE
      READ(31) HED
      PRINT *,HED
      READ(31) (KODE,I=1,NUMNP)
      READ(31) (X(I),I=1,NUMNP)
      READ(31) (Y(I),I=1,NUMNP)
      READ(31) (PSURF,I=1,NUMNP)
      READ(31) (BDROCK,I=1,NUMNP)
      READ(31) (KX(I,1),KX(I,2),KX(I,3),KX(I,4),I=1,NUMEL)
      READ(31) (IBFLUX1,IBFLUX2,BFLUX,I=1,NUMGBC)
      return
999   write(*,*) 'no file'
      stop
      end
      subroutine readrest(collen,rowlen,
     &           NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,
     &           Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,
     &           itime,iend,
     &           KX,X,Y)
      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER HED*80,HEDT*80,junk*80
      DIMENSION KODE(NMAX),X(NMAX),Y(NMAX),HTICE(NMAX),ADOT(NMAX),
     &          FRACT(NMAX),PSURF(NMAX),BDROCK(NMAX),FLOWA(NMAX),
     &          SLDGB(NMAX),KX(NMAX,4),CONST(NMAX),IBFLUX(NMAX,2),
     &          BFLUX(NMAX),TBED(NMAX),AFUDGE(NMAX),ACON(NMAX),
     &          ITYPE(NMAX),GEOFLUX(NMAX),CALV(NMAX),
     &          WTHICK(NMAX),BMELT(NMAX)
      DIMENSION VEL(NMAX,3)
      integer collen,rowlen
      parameter(RHOI=.917d0,RHOW=1.092d0,RATDEN=RHOW/RHOI)
      parameter(SEALEV=-10000.d0)
* data variables
      integer  itime
      real  Z1(collen*rowlen)
      real  Z2(collen*rowlen)
      real  Z3(collen*rowlen)
      real  Z4(collen*rowlen)
      real  Z5(collen*rowlen)
      real  Z6(collen*rowlen)
      real  Z7(collen*rowlen)
      real  Z8(collen*rowlen)
      real  Z9(collen*rowlen)
      if(collen*rowlen.ne.numnp) then 
        write(*,*) 'problems:'
        write(*,*) ' collen,rowlen,collen*rowlen',
     &            collen,rowlen,collen*rowlen
        write(*,*) ' numnp                      ',numnp
        stop
      endif
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME
      READ(33,END=999) HEDT
      write(*,*) HEDT
      write(junk,'(a)') hedt
      read(junk,2000) ttime
2000  format(6x,f20.0)
      itime=nint(ttime)
      READ(33) (HTICE(I),I=1,NUMNP)
      READ(33) (ADOT(I),I=1,NUMNP)
      READ(33) (BDROCK(I),I=1,NUMNP)
      READ(33) (CONST(I),I=1,NUMEL)
      READ(33) (ACON(I),I=1,NUMEL)
      READ(34) (FRACT(I),I=1,NUMNP)
      READ(34) (FLOWA(I),I=1,NUMNP)
      READ(34) (SLDGB(I),I=1,NUMNP)
      READ(34) (AFUDGE(I),I=1,NUMNP)
      READ(36,END=999) HEDT
      READ(36) (TBED(I),I=1,NUMNP)
      READ(36) (BMELT(I),I=1,NUMNP)
      READ(36) (WTHICK(I),I=1,NUMNP)
      ncount=0
      do i=1,numnp
        ncount=ncount+1
        THICKL=HTICE(I)-BDROCK(I)                           !
c       FLOT=(1.D0-RATDEN)*(BDROCK(I)-SEALEV)                        !
c       IF(HTICE(I).LT.FLOT) THEN                              !
c         THICKL=SEALEV                                       !
c       ENDIF                                                  !
        Z1(ncount)=htice(i)
        Z2(ncount)=bdrock(i)
        Z3(ncount)=thickl
        Z4(ncount)=adot(i)
        Z5(ncount)=flowa(i)
        Z6(ncount)=tbed(i)
        Z7(ncount)=bmelt(i)
        Z8(ncount)=wthick(i)
      enddo
      call DVELO(NMAX,NUMEL,KX,X,Y,HTICE,bdrock,CONST,VEL)
      ncount=0
      do i=1,numel
        ncount=ncount+1
        Z9(ncount)=VEL(i,1)
      enddo
      iend=0
      return
999   CONTINUE
      iend=1
      END


C==========================================================
      SUBROUTINE DVELO(NMAX,NUMEL,KX,XX,YY,HTICE,DEPB,CONST,VEL)
      IMPLICIT REAL*8(A-H,O-Z) 
      DIMENSION HTICE(NMAX),DEPB(NMAX),XX(NMAX),YY(NMAX)
      DIMENSION KX(NMAX,4),CONST(NMAX),VEL(NMAX,3)
      DIMENSION DPSIX(9),DPSIY(9),DXDS(2,2),DSDX(2,2),LM(4)             
      DIMENSION PSI(4),DPSI(4,2)                                        
      DIMENSION XY(2,4)
      VMAX=-1.d30
      HMIN=10.d0
      DO 600 J=1,NUMEL  
        VEL(J,1)=0D0
        VEL(J,2)=0D0
        VEL(J,3)=0D0
        NNODE=4
        CENTX=0.0D00
        CENTY=0.0D00
        HH=0.0D00
        SUMX=0.0D00
        SUMY=0.0D00
        DO I=1,NNODE
          LM(I)=KX(J,I)
        ENDDO
        I=LM(1)
        JJ=LM(2)
        K=LM(3)
        L=LM(4)
        XY(1,1)=XX(I)
        XY(1,2)=XX(JJ)
        XY(1,3)=XX(K)
        XY(1,4)=XX(L)
        XY(2,1)=YY(I)
        XY(2,2)=YY(JJ)
        XY(2,3)=YY(K)
        XY(2,4)=YY(L)
        XCENT=(XX(I)+XX(JJ)+XX(K)+XX(L))/4000.d0
        YCENT=(YY(I)+YY(JJ)+YY(K)+YY(L))/4000.d0
C
        CALL FESHAPE(1,CENTX,CENTY,PSI,DPSI)
        
      
C CALCULATE DXDS...EQUATION (5.3.6)
        DO I=1,2
          DO L=1,2
            DXDS(I,L)=0.0d0
            DO K=1,NNODE
              DXDS(I,L)=DXDS(I,L)+DPSI(K,L)*XY(I,K)
            ENDDO
          ENDDO
        ENDDO
   
C CALCULATE DSDX...EQUATION (5.2.7)
        DETJ=(DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1))
        IF (DETJ.LE.0.0) THEN
          WRITE(*,*) ' BAD JACOBIAN',DETJ,X
          STOP
        ENDIF
        DSDX(1,1)=DXDS(2,2)/DETJ
        DSDX(2,2)=DXDS(1,1)/DETJ
        DSDX(1,2)=-DXDS(1,2)/DETJ
        DSDX(2,1)=-DXDS(2,1)/DETJ
C CALCULATE D(PSI)/DX...EQUATION (5.3.5)
        DO I=1,NNODE
          DPSIX(I)=DPSI(I,1)*DSDX(1,1)+DPSI(I,2)*DSDX(2,1)
          DPSIY(I)=DPSI(I,1)*DSDX(1,2)+DPSI(I,2)*DSDX(2,2)
        ENDDO
        DO I=1,NNODE
          SUMX=SUMX + HTICE(LM(I))*DPSIX(I)
          SUMY=SUMY + HTICE(LM(I))*DPSIY(I)
          HH=HH + (HTICE(LM(I))-DEPB(LM(I)))*PSI(I)
        ENDDO
C
        DELH=SUMX**2 + SUMY**2
        DELH=SQRT(DELH)
        IF(HH.GT.HMIN) THEN
          UX=-CONST(J)*SUMX/HH
          UY=-CONST(J)*SUMY/HH
        ELSE
          UX=0.d0
          UY=0.d0
        ENDIF
        UMAG=SQRT(UX**2+UY**2)
        VEL(J,1)=UMAG
        VEL(J,2)=UX
        VEL(J,3)=UY
        VMAX=MAX(VMAX,UMAG)
600   CONTINUE
      END
c----------------------------------------------------
      SUBROUTINE FESHAPE(NTYPE,XI,ET,PSI,DPSI)
C ELEMENT SHAPE FUNCTIONS AND DERIVATIVES AT LOCAL COORDINATES (XI,ET)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PSI(4),DPSI(4,2)
      IF(NTYPE.EQ.1) THEN
        PSI(1)=.25d0*(1.d0-XI)*(1.d0-ET)
        PSI(2)=.25d0*(1.d0+XI)*(1.d0-ET)
        PSI(3)=.25d0*(1.d0+XI)*(1.d0+ET)
        PSI(4)=.25d0*(1.d0-XI)*(1.d0+ET)
        DPSI(1,2)=-.25d0*(1.d0-XI)
        DPSI(2,2)=-.25d0*(1.d0+XI)
        DPSI(3,2)=.25d0*(1.d0+XI)
        DPSI(4,2)=.25d0*(1.d0-XI)
        DPSI(1,1)=-.25d0*(1.d0-ET)
        DPSI(2,1)=.25d0*(1.d0-ET)
        DPSI(3,1)=.25d0*(1.d0+ET)
        DPSI(4,1)=-.25d0*(1.d0+ET)
      ELSE
        PSI(1)=1.d0-XI-ET
        PSI(2)=XI
        PSI(3)=ET
        DPSI(1,2)=-1.d0
        DPSI(2,2)=0.d0
        DPSI(3,2)=1.d0
        DPSI(1,1)=-1.d0
        DPSI(2,1)=1.d0
        DPSI(3,1)=0.d0
      ENDIF
      END

