      include "parameter.h"
      PARAMETER(NMAX=MAXNUM)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER HED*80
      DIMENSION KODE(NMAX),X(NMAX),Y(NMAX),HTICE(NMAX),ADOT(NMAX),
     &          FRACT(NMAX),PSURF(NMAX),BDROCK(NMAX),FLOWA(NMAX),
     &          SLDGB(NMAX),
     &          KX(NMAX,4),CONST(NMAX),IBFLUX(NMAX,2),BFLUX(NMAX),
     &          TBED(NMAX),BMELT(NMAX),WTHICK(NMAX),
     &          ACON(NMAX)
      icnt=0
      icount=0
      print *,'input first, display interval and how many to show'
      read(*,*) ifirst,istep,inumber
      print *,'input 1 for curvilinear, 0 for uniform'
      read(*,*) ibound
      if(ibound.eq.1) then
        write(*,*) 'surface display'
        WRITE(*,*) 'INPUT  1 FOR CALCULATED SURFACE,'
        WRITE(*,*) '       2 FOR MASS BALANCE,'
        WRITE(*,*) '       3 FOR FLOW CONSTANT,'
        WRITE(*,*) '       4 FOR DIFFERENCE,'
        WRITE(*,*) '       5 FOR THICKNESS,'
        WRITE(*,*) '       6 FOR PRESENT SURFACE,'
        WRITE(*,*) '      12 FOR SLIDING FRACTION,'
        WRITE(*,*) '      14 FOR BASAL TEMPERATURE,'
        WRITE(*,*) '      15 FOR BASAL MELT RATE,'
        WRITE(*,*) '      16 FOR WATER THICKNESS,'
        READ(*,*) ISPLOT
        write(*,*) 'basal display'
        WRITE(*,*) '       1 FOR BEDROCK,'
        WRITE(*,*) '       2 FOR SLIDING FRACTION,'
        WRITE(*,*) '       3 FOR SLIDING CONSTANT,'
        WRITE(*,*) '       4 FOR BASAL TEMPERATURE,'
        WRITE(*,*) '       5 FOR BASAL MELT RATE,'
        WRITE(*,*) '       6 FOR WATER THICKNESS,'
        READ(*,*) IBPLOT

      endif
C READ INPUT HEADER
      READ(30,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
c      READ(30,1000) HED,NUMNP,NUMEL,NUMLEV,NUMCOL,NUMGBC,NDT,INTER,DT
      PRINT *,HED
1000  FORMAT (A80,/,7I6,F8.0)
      IF(NUMNP.GT.NMAX) THEN
        PRINT *,'NUMNP=',NUMNP,' NMAX=',NMAX,' INCREASE NMAX'
        STOP
      ENDIF
C READ INPUT GRID, THINGS THAT NEVER CHANGE
      READ(31) HED
      PRINT *,HED
      READ(31) (KODE(I),I=1,NUMNP)
      READ(31) (X(I),I=1,NUMNP)
      READ(31) (Y(I),I=1,NUMNP)
      READ(31) (PSURF(I),I=1,NUMNP)
      READ(31) (BDROCK(I),I=1,NUMNP)
      READ(31) (KX(I,1),KX(I,2),KX(I,3),KX(I,4),I=1,NUMEL)
      READ(31) (IBFLUX(I,1),IBFLUX(I,2),BFLUX(I),I=1,NUMGBC)
C READ INPUT DIFFER, THINGS THAT DIFFER FROM ONE GRID TO THE NEXT
      READ(32) HED
      PRINT *,HED
      READ(32) (ADOT(I),I=1,NUMNP)
      READ(32) (FRACT(I),I=1,NUMNP)
      READ(32) (FLOWA(I),I=1,NUMNP)
      READ(32) (SLDGB(I),I=1,NUMNP)
C READ INPUT TIME, THINGS THAT CHANGE WITH TIME
10    continue
      READ(33,END=999) HED
      PRINT *,HED
      READ(33) (HTICE(I),I=1,NUMNP)
      READ(33) (ADOT(I),I=1,NUMNP)
      READ(33) (BDROCK(I),I=1,NUMNP)
      READ(33) (CONST(I),I=1,NUMEL)
      READ(33) (ACON(I),I=1,NUMEL)
      READ(34) (FRACT(I),I=1,NUMNP)
      READ(34) (FLOWA(I),I=1,NUMNP)
      READ(34) (SLDGB(I),I=1,NUMNP)
      READ(34) (AFUDGE,I=1,NUMNP)
      READ(36,END=999) HED
      READ(36) (TBED(I),I=1,NUMNP)
      READ(36) (BMELT(I),I=1,NUMNP)
      READ(36) (WTHICK(I),I=1,NUMNP)
      icount=icount+1
      if(icount.gt.ifirst-1+inumber*istep) goto 999
      if(icount.ge.ifirst .and. 
     *   1+mod(icount-ifirst,istep).eq.1) then
        print *,'outputting for display',icount,icnt+1
c--------------------------------------      
      xmin=1e30
      xmax=-xmin
      ymin=1e30
      ymax=-ymin
      do i=1,numnp
        xmin=min(xmin,x(i))
        xmax=max(xmax,x(i))
        ymin=min(ymin,y(i))
        ymax=max(ymax,y(i))
      enddo
      if(ibound.eq.0) then
        facth=1.
        factx=1.
        facty=1.
      else
      facth=1./1000.
      factx=1./1000./100.
      facty=1./1000./100.
      endif
      write(11,*) 'RANK 2'
      write(11,*) 'DIMENSIONS ',NUMCOL,NUMLEV
      if(ibound.eq.0) then
c        write(11,*) 'BOUNDS ',real(factx*xmin),real(factx*xmax)
c     &                       ,real(facty*ymin),real(facty*ymax)
        write(11,*) 'BOUNDS',0,numcol*10000,0,numlev*10000
      endif
      write(11,*) 'NAME BEDROCK'
c      write(11,*) 'TIME',icnt
      write(11,*) 'TIME',hed(6:20)
      write(*,*) '1234123456789012345678901234567890'
      write(*,*) 'TIME',hed(6:20)
      write(11,*) 'SCALAR'
      WRITE(11,*) 'DATA'
      if(ibplot.eq.1 .or. ibound.eq.0) then
        WRITE(11,*) (real(facth*BDROCK(I)),I=1,NUMNP)
      elseif(ibplot.eq.2) then
        WRITE(11,*) (real(FRACT(I)),I=1,NUMNP)
      elseif(ibplot.eq.3) then
        WRITE(11,*) (real(SLDGB(I)),I=1,NUMNP)
      elseif(ibplot.eq.4) then
        WRITE(11,*) (real(TBED(I)),I=1,NUMNP)
      elseif(ibplot.eq.5) then
        WRITE(11,*) (real(1000.*BMELT(I)),I=1,NUMNP)
      elseif(ibplot.eq.6) then
        WRITE(11,*) (real(WTHICK(I)),I=1,NUMNP)
      endif
      if(ibound.eq.1) then
        write(11,*) 'INTERLACED'
        write(11,*) 'VECTOR 3'
        WRITE(11,*) 'GRID'
        do i=1,numnp
c          WRITE(11,*) factx*X(I),facty*Y(I),facth*bdrock(i)
          WRITE(11,*) 1e-5*X(I),1e-5*Y(I),facth*bdrock(i)
        enddo
      endif
      WRITE(11,*) 'END'
      write(12,*) 'RANK 2'
      write(12,*) 'DIMENSIONS ',NUMCOL,NUMLEV
      if(ibound.eq.0) then
c        write(12,*) 'BOUNDS ',real(factx*xmin),real(factx*xmax)
c     &                       ,real(facty*ymin),real(facty*ymax)
        write(12,*) 'BOUNDS',0,numcol*10000,0,numlev*10000
      endif
      write(12,*) 'NAME SURF'
c      write(12,*) 'TIME',icnt
      write(12,*) 'TIME',hed(6:20)
      write(12,*) 'SCALAR'
      WRITE(12,*) 'DATA'
      if(isplot.eq.1 .or. ibound.eq.0) then
        WRITE(12,*) (real(facth*htice(I)),I=1,NUMNP)
      elseif(isplot.eq.2) then
        WRITE(12,*) (real(adot(I)),I=1,NUMNP)
      elseif(isplot.eq.3) then
        WRITE(12,*) (real(flowa(I)),I=1,NUMNP)
      elseif(isplot.eq.4) then
        WRITE(12,*) (real(facth*(htice(I)-psurf(i))),I=1,NUMNP)
      elseif(isplot.eq.5) then
        WRITE(12,*) (real(facth*(htice(I)-bdrock(i))),I=1,NUMNP)
      elseif(isplot.eq.6) then
        WRITE(12,*) (real(facth*psurf(I)),I=1,NUMNP)
      elseif(isplot.eq.12) then
        WRITE(12,*) (real(FRACT(I)),I=1,NUMNP)
      elseif(isplot.eq.14) then
        WRITE(12,*) (real(TBED(I)),I=1,NUMNP)
      elseif(isplot.eq.15) then
        WRITE(12,*) (real(1000.*BMELT(I)),I=1,NUMNP)
      elseif(isplot.eq.16) then
        WRITE(12,*) (real(WTHICK(I)),I=1,NUMNP)
      endif
      if(ibound.eq.1) then
        write(12,*) 'INTERLACED'
        write(12,*) 'VECTOR 3'
        WRITE(12,*) 'GRID'
        do i=1,numnp
c          WRITE(12,*) factx*X(I),facty*Y(I),facth*htice(i)
          WRITE(12,*) 1e-5*X(I),1e-5*Y(I),facth*htice(i)
        enddo
      endif
      WRITE(12,*) 'END'
      icnt=icnt+1
c------------------------
      endif
      goto 10
999   continue
      print *,'output this many',icnt
      END
