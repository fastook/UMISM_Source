      program packallnc
c     IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM,N3=NMAX*3)
      CHARACTER HED*80,SCRTCH*80,JUNK*80
      LOGICAL CTOGG,ITOGG
      integer WTOGG,BTOGG
      COMMON /INPUT_DEF/ AMASS(11),ACOM,WJUNK,WINDIR(2),
     &      XPOLE,YPOLE,CFACTOR,ALPHAC(3),TBASE,
     &      TNSLBASE,CTOGG,WTOGG,ITOGG,BTOGG
      COMMON /INPUT_HEAD/ DT,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &      INTER
      COMMON /INPUT_GRID/ X(NMAX),Y(NMAX),PSURF(NMAX),BDROCK(NMAX),
     &      BFLUX(NMAX),KX(NMAX,4),IBFLUX(NMAX,2),KODE(NMAX)
      COMMON /INPUT_DIFF/ ADOTB(NMAX),FRACTB(NMAX),FLOWAB(NMAX),
     &      SLDGBB(NMAX),STEMP(NMAX),AFUDGEB(NMAX),GEOFLUX(NMAX),
     &      CALV(NMAX),ITYPE(NMAX)
      COMMON /INPUT_TIME/ HTICEB(NMAX),CONSTB(NMAX),ACONB(NMAX)
      COMMON /DEP_TIME/ WWW(N3),THICK(NMAX),WRATE(N3,2)
      COMMON /OUT_TIME/ HTICE(NMAX),ADOT(NMAX),DEPB(NMAX),
     &      CONST(NMAX),ACON(NMAX)
      COMMON /BC_TIME/  FRACT(NMAX),FLOWA(NMAX),SLDGB(NMAX),
     &      AFUDGE(NMAX)
      COMMON /TEMP_TIME/ TBED(NMAX),BMELT(NMAX),WTHICK(NMAX)
      COMMON /LATLONG/ RLAT(NMAX),RLONG(NMAX)
      REAL*8 AMASS,ACOM,WJUNK,WINDIR,XPOLE,YPOLE,CFACTOR,
     &       ALPHAC,TBASE,TNSLBASE,DT
      REAL*8 X,Y,PSURF,BDROCK,
     &      BFLUX,ADOTB,FRACTB,FLOWAB,
     &      SLDGBB,STEMP,AFUDGEB,GEOFLUX,
     &      CALV,HTICEB,CONSTB,ACONB,
     &      WWW,THICK,WRATE,
     &      HTICE,ADOT,DEPB,
     &      CONST,ACON,FRACT,FLOWA,SLDGB,
     &      AFUDGE,TBED,BMELT,WTHICK,RLAT,RLONG
      include 'netcdf.inc'
* error status return
      integer iret
* netCDF id
      integer ncid
* dimension ids
      integer rowdim,coldim,time_dim
* dimension lengths
      integer rowlen,collen,time_len
      parameter(nrowmax=281,ncolmax=281)
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
      integer row(nrowmax),col(ncolmax)
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
      call READER(IEOF,ITIME)
c     collen =numcol
c     rowlen =numlev
* enter define mode
      iret = nf_create('packall.nc', NF_CLOBBER, ncid)
      call check_err(iret)
* define dimensions
      iret = nf_def_dim(ncid, 'time', time_len, time_dim)
      call check_err(iret)
* define variables
      time_dims(1) = time_dim
      iret=nf_def_var(ncid,'time',NF_INT,time_rank,time_dims,time_id)
      call check_err(iret)
* define dimensions
      iret = nf_def_dim(ncid, 'row', numlev, rowdim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'col', numcol, coldim)
      call check_err(iret)
* define variables
      rowdims(1) = rowdim
      iret=nf_def_var(ncid,'row',NF_INT,rowrank,rowdims,rowid)
      call check_err(iret)
      coldims(1) = coldim
      iret=nf_def_var(ncid,'col',NF_INT,colrank,coldims,colid)
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
      do i=1,numlev
        row(i)=i
      enddo
c     data row /0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90/
      iret = nf_put_var_int(ncid, rowid, row)
      call check_err(iret)
* store col
      do i=1,numcol
        col(i)=i
      enddo
c     data col /0, 10, 20, 30, 40, 50, 60, 70, 80, 90/
      iret = nf_put_var_int(ncid, colid, col)
      call check_err(iret)
       
* Write record variables
      call writerecs(ncid,time_id,
     &     Z1_id,Z2_id,Z3_id,Z4_id,Z5_id,Z6_id,Z7_id,Z8_id,Z9_id,
     &     z1,z2,z3,z4,z5,z6,z7,z8,
     &     collen,rowlen)
       
      iret = nf_close(ncid)
      call check_err(iret)
      end
********************************************************
      subroutine writerecs(ncid,time_id,
     &     Z1_id,Z2_id,Z3_id,Z4_id,Z5_id,Z6_id,Z7_id,Z8_id,Z9_id,
     &     z1,z2,z3,z4,z5,z6,z7,z8,
     &     collen,rowlen)
c      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM,N3=NMAX*3)
      CHARACTER HED*80,SCRTCH*80,JUNK*80
      LOGICAL CTOGG,ITOGG
      integer WTOGG,BTOGG
      COMMON /INPUT_DEF/ AMASS(11),ACOM,WJUNK,WINDIR(2),
     &      XPOLE,YPOLE,CFACTOR,ALPHAC(3),TBASE,
     &      TNSLBASE,CTOGG,WTOGG,ITOGG,BTOGG
      COMMON /INPUT_HEAD/ DT,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &      INTER
      COMMON /INPUT_GRID/ X(NMAX),Y(NMAX),PSURF(NMAX),BDROCK(NMAX),
     &      BFLUX(NMAX),KX(NMAX,4),IBFLUX(NMAX,2),KODE(NMAX)
      COMMON /INPUT_DIFF/ ADOTB(NMAX),FRACTB(NMAX),FLOWAB(NMAX),
     &      SLDGBB(NMAX),STEMP(NMAX),AFUDGEB(NMAX),GEOFLUX(NMAX),
     &      CALV(NMAX),ITYPE(NMAX)
      COMMON /INPUT_TIME/ HTICEB(NMAX),CONSTB(NMAX),ACONB(NMAX)
      COMMON /DEP_TIME/ WWW(N3),THICK(NMAX),WRATE(N3,2)
      COMMON /OUT_TIME/ HTICE(NMAX),ADOT(NMAX),DEPB(NMAX),
     &      CONST(NMAX),ACON(NMAX)
      COMMON /BC_TIME/  FRACT(NMAX),FLOWA(NMAX),SLDGB(NMAX),
     &      AFUDGE(NMAX)
      COMMON /TEMP_TIME/ TBED(NMAX),BMELT(NMAX),WTHICK(NMAX)
      COMMON /LATLONG/ RLAT(NMAX),RLONG(NMAX)
      REAL*8 AMASS,ACOM,WJUNK,WINDIR,XPOLE,YPOLE,CFACTOR,
     &       ALPHAC,TBASE,TNSLBASE,DT
      REAL*8 X,Y,PSURF,BDROCK,
     &      BFLUX,ADOTB,FRACTB,FLOWAB,
     &      SLDGBB,STEMP,AFUDGEB,GEOFLUX,
     &      CALV,HTICEB,CONSTB,ACONB,
     &      WWW,THICK,WRATE,
     &      HTICE,ADOT,DEPB,
     &      CONST,ACON,FRACT,FLOWA,SLDGB,
     &      AFUDGE,TBED,BMELT,WTHICK,RLAT,RLONG
      REAL*8 VEL(NMAX,3)
      REAL*8 THICKL(NMAX)
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
      parameter(RHOI=.917,RHOW=1.092,RATDEN=RHOW/RHOI,SEALEV=0.0)
      ntime=0
      print *,'start with???'
      read *,numstrt
      print *,'how many???'
      read *,number
* skip first one ...
      do i=1,numstrt
        call READER(IEOF,ITIME)
      enddo
* skip first one ...
      call READER(IEOF,ITIME)
      dowhile(IEOF.eq.0 .and. ntime.le.number)
        ntime=ntime+1
        do i=1,numnp
          THICKL(i)=HTICE(I)-DEPB(I)
          FLOT=(1.D0-RATDEN)*(DEPB(I)-SEALEV)
          IF(HTICE(I).LT.FLOT) THEN
            THICKL(i)=SEALEV
          ENDIF
        ENDDO
        call load(numnp,z1,htice)
        call load(numnp,z2,depb)
        call load(numnp,z3,thickl)
        call load(numnp,z4,adot)
        call load(numnp,z5,flowa)
        call load(numnp,z6,tbed)
        call load(numnp,z7,bmelt)
        call load(numnp,z8,wthick)
        ncount=0
        do j=1,numlev
          do i=1,numcol
            IF(Z2(I,J,1).EQ.-9999.) Z2(I,J,1)=NF_FILL_REAL
          enddo
        enddo
        call DVELO(NMAX,NUMEL,KX,X,Y,HTICE,DEPB,CONST,VEL)
        call load(numel,z9,vel(1,1))
c       ncount=0
c       do j=1,numlev-0
c         do i=1,numcol-0
c           ncount=ncount+1
c           Z9(I,J,1)=VEL(ncount,1)
c         enddo
c       enddo
        time(1)=itime
        time_start(1)=ntime
        time_count(1)=1
        iret=nf_put_vara_int(ncid,time_id,time_start,time_count,time)
        call check_err(iret)
        Z_start(1)=1
        Z_start(2)=1
        Z_start(3)=ntime
        Z_count(1)=numcol
        Z_count(2)=numlev
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
        Z_count(1)=numcol-1
        Z_count(2)=numlev-1
        Z_count(3)=1
        iret=nf_put_vara_real(ncid,Z9_id,Z_start,Z_count,Z9)
        call READER(IEOF,ITIME)
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





      SUBROUTINE READER(IEOF,ITIME)
      IMPLICIT REAL*8(A-H,O-Z)
      include "parameter.h"
      PARAMETER( NMAX=MAXNUM,N3=NMAX*3)
      CHARACTER HED*80,SCRTCH*80,JUNK*80
      LOGICAL CTOGG,ITOGG
      integer WTOGG,BTOGG
      COMMON /INPUT_DEF/ AMASS(11),ACOM,WJUNK,WINDIR(2),
     &      XPOLE,YPOLE,CFACTOR,ALPHAC(3),TBASE,
     &      TNSLBASE,CTOGG,WTOGG,ITOGG,BTOGG
      COMMON /INPUT_HEAD/ DT,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &      INTER
      COMMON /INPUT_GRID/ X(NMAX),Y(NMAX),PSURF(NMAX),BDROCK(NMAX),
     &      BFLUX(NMAX),KX(NMAX,4),IBFLUX(NMAX,2),KODE(NMAX)
      COMMON /INPUT_DIFF/ ADOTB(NMAX),FRACTB(NMAX),FLOWAB(NMAX),
     &      SLDGBB(NMAX),STEMP(NMAX),AFUDGEB(NMAX),GEOFLUX(NMAX),
     &      CALV(NMAX),ITYPE(NMAX)
      COMMON /INPUT_TIME/ HTICEB(NMAX),CONSTB(NMAX),ACONB(NMAX)
      COMMON /DEP_TIME/ WWW(N3),THICK(NMAX),WRATE(N3,2)
      COMMON /OUT_TIME/ HTICE(NMAX),ADOT(NMAX),DEPB(NMAX),
     &      CONST(NMAX),ACON(NMAX)
      COMMON /BC_TIME/  FRACT(NMAX),FLOWA(NMAX),SLDGB(NMAX),
     &      AFUDGE(NMAX)
      COMMON /TEMP_TIME/ TBED(NMAX),BMELT(NMAX),WTHICK(NMAX)
      COMMON /LATLONG/ RLAT(NMAX),RLONG(NMAX)
      DATA IPASS /1/
      SAVE IPASS
      IF(IPASS.EQ.1) THEN
       ipass=2
       open(9,file='input.def')
       open(30,file='input.head')
       open(31,file='input.grid',form='unformatted',access='sequential')
       open(32,file='input.diff',form='unformatted',access='sequential')
       open(33,file='input.time',form='unformatted',access='sequential')
       open(88,file='dep.time',form='unformatted',access='sequential')
       open(34,file='out.time',form='unformatted',access='sequential')
       open(35,file='bc.time',form='unformatted',access='sequential')
       open(36,file='temp.time',form='unformatted',access='sequential')


C ..... INPUT.DEF ...
        READ(9,*,END=991) AMASS,ACOM,WJUNK,WINDIR(2),
     &                    XPOLE,YPOLE,CFACTOR,ALPHAC,TBASE,
     &                    TNSLBASE,CTOGG,WTOGG,ITOGG,BTOGG
C ..... INPUT.HEAD ...
        READ(30,1000,END=992) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &                        INTER,DT
c       print *,' head ',hed
        print *,hed
1000    FORMAT (A80,/,7I6,F8.0)
C ..... INPUT.GRID ...
C ..... READ INPUT GRID, THINGS THAT NEVER CHANGE
        READ(31,END=993) HED
c       print *,' grid ',hed
        READ(31) (KODE(I),I=1,NUMNP)
        READ(31) (X(I),I=1,NUMNP)
        READ(31) (Y(I),I=1,NUMNP)
        READ(31) (PSURF(I),I=1,NUMNP)
        READ(31) (BDROCK(I),I=1,NUMNP)
        READ(31) (KX(I,1),KX(I,2),KX(I,3),KX(I,4),I=1,NUMEL)
        READ(31) (IBFLUX(I,1),IBFLUX(I,2),BFLUX(I),I=1,NUMGBC)
        call setrig
        do i=1,numnp
            CALL RECPOL(0.001d0*X(i),0.001d0*Y(i),RLAT(i),RLONG(i))
        enddo
C ..... INPUT.DIFF ...
C ..... READ INPUT DIFFER, THINGS THAT DIFFER FROM ONE GRID TO THE NEXT
        READ(32,END=994) HED
c       print *,' diff ',hed
        READ(32) (ADOTB(I),I=1,NUMNP)
        READ(32) (FRACTB(I),I=1,NUMNP)
        READ(32) (FLOWAB(I),I=1,NUMNP)
        READ(32) (SLDGBB(I),I=1,NUMNP)
        READ(32) (STEMP(I),I=1,NUMNP)
        READ(32) (ITYPE(I),I=1,NUMNP)
        READ(32) (AFUDGEB(I),I=1,NUMNP)
        READ(32) (GEOFLUX(I),I=1,NUMNP)
        READ(32) (CALV(I),I=1,NUMNP)
C ..... INPUT.TIME ...
C ..... READ INPUT TIME, THINGS THAT CHANGE WITH TIME
        READ(33,END=995) HED
c       print *,' time ',hed
        READ(33) (HTICEB(I),I=1,NUMNP)
        READ(33) (RJUNK,I=1,NUMNP)
        READ(33) (RJUNK,I=1,NUMNP)
        READ(33) (CONSTB(I),I=1,NUMEL)
        READ(33) (ACONB(I),I=1,NUMEL)
        return
      ENDIF

1     CONTINUE
C ..... DEP.TIME ...
        READ(88,END=996) SCRTCH
        PRINT *,SCRTCH
c        PRINT *,' DEP ',SCRTCH
        WRITE(JUNK,'(A)') SCRTCH
        READ(JUNK,2000) TTIME
2000    FORMAT(7X,F20.0)
        ITIME=NINT(TTIME)
        READ(88) (WWW(I),I=1,3*NUMNP,3)
        READ(88) (THICK(I),I=1,NUMNP)
        READ(88) (WRATE(I,1),I=1,3*NUMNP,3)
C ..... OUT.TIME ...
        READ(34,END=997) SCRTCH
c        print *,' out2 ',SCRTCH
        READ(34) (HTICE(I),I=1,NUMNP)
        READ(34) (ADOT(I),I=1,NUMNP)
        READ(34) (DEPB(I),I=1,NUMNP)
        READ(34) (CONST(I),I=1,NUMEL)
        READ(34) (ACON(I),I=1,NUMEL)
C ..... BC.TIME ...
        READ(35,END=998) (FRACT(I),I=1,NUMNP)
        READ(35) (FLOWA(I),I=1,NUMNP)
        READ(35) (SLDGB(I),I=1,NUMNP)
        READ(35) (AFUDGE(I),I=1,NUMNP)
C ..... TEMP.TIME ...
        READ(36,END=999) SCRTCH
c        print *,' temp ',SCRTCH
        READ(36) (TBED(I),I=1,NUMNP)
        READ(36) (BMELT(I),I=1,NUMNP)
        READ(36) (WTHICK(I),I=1,NUMNP)
C       GOTO 1
      IEOF=0
      RETURN
991   CONTINUE
      IEOF=1
      print *,'eof: def'
      RETURN
992   CONTINUE
      IEOF=2
      print *,'eof: head'
      RETURN
993   CONTINUE
      IEOF=3
      print *,'eof: grid'
      RETURN
994   CONTINUE
      IEOF=4
      print *,'eof: diff'
      RETURN
995   CONTINUE
      IEOF=5
      print *,'eof: time'
      RETURN
996   CONTINUE
      IEOF=6
      print *,'eof: dep'
      RETURN
997   CONTINUE
      IEOF=7
      print *,'eof: out2'
      RETURN
998   CONTINUE
      IEOF=8
      print *,'eof: bc'
      RETURN
999   CONTINUE
      IEOF=9
      print *,'eof: temp'
      RETURN
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
c----------------------------------------------------
      subroutine load(numnp,z1,z)
      real z1(numnp)
      real*8 z(numnp)
      do i=1,numnp
        z1(i)=z(i)
      enddo
      end
c----------------------------------------------------

      SUBROUTINE SETRIG
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      PI=4.D0*ATAN(1.D0)
      RADIUS=2.D4/PI
      RADIUS=RADIUS*0.53
      CIRCUM=2.D0*PI*RADIUS
      RKMPDEG=CIRCUM/360.D0
      RADPDEG=PI/180.D0
      DEGPRAD=180.D0/PI
      END
      SUBROUTINE POLREC(RLAT,RLONG,X,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      y= 1000*rkmpdeg*rlat
      x= 1000*rkmpdeg*cos(rlat*radpdeg)*(rlong+127.5)
      END
      SUBROUTINE RECPOL(X,Y,RLAT,RLONG)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD
      rlat=y*0.001/rkmpdeg
      rlong=-127.5+x*0.001/rkmpdeg/cos(rlat*radpdeg)
      END
