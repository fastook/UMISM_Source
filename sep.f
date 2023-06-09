      parameter(nmax=29999)
      CHARACTER HED*80
      dimension flowa(nmax),htice(nmax)
 1    READ(1,100) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,
     &INTER,DT
      PRINT *,HED
      print *,'input sliding threshhold'
      read(*,*,end=999) thresh
100   FORMAT(A80,/,7I6,E15.6)
      FMIN=1.E30
      FMAX=-FMIN
      SMIN=FMIN
      SMAX=-FMIN
      icountf=0
      icounts=0
      inode=0
      favg=0.
      savg=0.
      DO 10 NUM=1,NUMNP
      READ(1,200) N,KODE,XX,YY,htice(num),ADOT,FRACT,PSURF,BDROCK,
     &            flowa(num),SLDGB
200   FORMAT(I6,I4,2E12.5,F10.2,F7.2,F9.2,F10.3,F10.1,F10.2,F10.3)
      if(htice(num).gt.0.0) then
        inode=inode+1
        if(fract.le.thresh) then
c flow dominates
        FMIN=MIN(FMIN,flowa(num))
        FMAX=MAX(FMAX,flowa(num))
        icountf=icountf+1
        favg=favg+flowa(num)
        else
c sliding dominates
        SMIN=MIN(SMIN,SLDGB)
        SMAX=MAX(SMAX,SLDGB)
        icounts=icounts+1
        savg=savg+sldgb
        endif
      endif
10    CONTINUE
      if(icountf.ne.0) favg=favg/icountf
      if(icounts.ne.0) savg=savg/icounts
      print *,'flow, min, avg, max',FMIN,FAVG,FMAX,icountf,inode,numnp
      print *,'slde, min, avg, max',SMIN,SAVG,SMAX,icounts
      icf=0
      do num=1,numnp
        if(htice(num).gt.0.) then
          if(flowa(num).eq.fmax) icf=icf+1
        endif
      enddo
      print *,'flow max number',icf
      DO 12 I=1,NUMEL
      READ(1,300) NUM,KX1,KX2,KX3,KX4
300   FORMAT(5I6)
12    CONTINUE
      DO 13 N=1,NUMGBC
      READ(1,310) I,J,RJUNK
13    continue
310   FORMAT(2I6,E13.6)
      rewind 1
      goto 1
999   continue
      END
