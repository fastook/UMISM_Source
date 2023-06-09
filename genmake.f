      implicit real*8(a-h,o-z)
      include "parameter.h" 
      PARAMETER(NMAX=MAXNUM)
      character*1 lend
      CHARACTER HED*80
      CHARACTER*80,JUNK1,JUNK2,JUNK3,JUNK4,JUNK5,JUNK6
      dimension x(nmax),y(nmax),kx(nmax,4)
      dimension rlonge(4),rlate(4)
c ... 0 for Northern hemisphere, 1 for Southern ....
      ihemi=1
c ... 0 for Western hemisphere, 1 for Eastern ...
      jhemi=0
      print *,'Northern:0 or Southern:1'
      read *,ihemi
      if(ihemi.eq.0) then
        print *,'Western:0 or Eastern:1'
        read *,jhemi
      endif
      if(ihemi.eq.0) then
        call setrig( 90.d0,0.d0)
      else
        call setrig(-90.d0,0.d0)
      endif
      lend=char(92)
c     open(1,file='inputnamer.data')
c     READ(1,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      READ(1,'(a)') HED
      READ(1,*) NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
      PRINT *,HED
1000  FORMAT (A80,/,7I6,F8.1)
      print *,numnp,' nodes,',numel,' elements'
      xmin=1e30
      xmax=-1e30
      ymin=1e30
      ymax=-1e30
      rlatmin=1e30
      rlatmax=-1e30
      rlongmin=1e30
      rlongmax=-1e30
      DO N=1,NUMNP
c       READ(1,1001) NUM,KODE,X(n),Y(n),HTICE,
        READ(1,*) NUM,KODE,X(n),Y(n),HTICE,
     &                   ADOT,FRACT,PSURF,
     &                   BDROCK,FLOWA,SLDGB,
     &                   temp,itype,AFUDGE,GEOFLUX,
     &                   calv
        call recpol(1*x(n),1*y(n),
     &                rlatn,rlongn)
c       print *,n,x(n),y(n),rlatn,rlongn
c       if(ihemi.eq.1) then
c         rlatn=-rlatn
c         rlongn=90.-rlongn
c         if(rlongn.lt.0.) rlongn=rlongn+360.
c       endif
        rlatmin=min(rlatmin,rlatn)
        rlatmax=max(rlatmax,rlatn)
        rlongmin=min(rlongmin,rlongn)
        rlongmax=max(rlongmax,rlongn)
        xmin=min(xmin,x(n))
        xmax=max(xmax,x(n))
        ymin=min(ymin,y(n))
        ymax=max(ymax,y(n))
1001  FORMAT(I6,I4,1P2E12.5,0PF10.2,F7.2,F9.2,F10.3,F10.1,F10.5,2F10.5,
     &       I5,F10.3,F10.3,1PE10.3)
      enddo
      print *,rlongmin,rlongmax,rlatmin,rlatmax
c     pause
      rlatmin=max(-90.,rlatmin-1)
      rlatmax=min(90.,rlatmax+1)
      rlongmin=max(0.,rlongmin-1)
      rlongmax=min(360.,rlongmax+1)
      DISTMAX=0
      DISTMIN=1E30
      DO N=1,NUMEL
c       READ(1,1002) NUM,(KX(n,i),i=1,4),CONST,ACON
c       READ(1,1002) NUM,(KX(n,i),i=1,4)
        READ(1,*) NUM,(KX(n,i),i=1,4)
        do i=1,4
          call recpol(1*x(kx(n,i)),1*y(kx(n,i)),
     &                rlate(i),rlonge(i))
c         if(ihemi.eq.1) then
c           rlate(i)=-rlate(i)
c           rlonge(i)=90.-rlonge(i)
c           if(rlonge(i).lt.0.) rlonge(i)=rlonge(i)+360.
c         endif
c         write(99,*) i,rlate(i),rlonge(i)
        enddo
        if(n.eq.numel/2) then
          DIST1=sqrt((rlate(1)-rlate(3))**2+
     &          (rlonge(1)-rlonge(3))**2)
          DIST2=sqrt((rlate(2)-rlate(4))**2+
     &          (rlonge(2)-rlonge(4))**2)
          DIST3=sqrt((rlate(2)-rlate(1))**2+
     &          (rlonge(2)-rlonge(1))**2)
          DIST4=sqrt((rlate(3)-rlate(2))**2+
     &          (rlonge(3)-rlonge(2))**2)
          DIST5=sqrt((rlate(4)-rlate(3))**2+
     &          (rlonge(4)-rlonge(3))**2)
          DIST6=sqrt((rlate(1)-rlate(4))**2+
     &          (rlonge(1)-rlonge(4))**2)
          DISTMAX=MAX(DISTMAX,DIST1,DIST3,DIST4,DIST5,DIST6)
          DISTMIN=MIN(DISTMIN,DIST1,DIST2,DIST3,DIST4,DIST5,DIST6)
        endif
c       write(99,*) n,dist1,dist2
c       write(99,*) n,dist3,dist4
c       write(99,*) n,dist5,dist6,distmax
 1002 FORMAT(5I6,1P2E17.10)
      enddo
      distmax=0.3*distmax
      dmax=distmax
      print *,'distmax:',distmax,distmax*60
      if(distmax.gt.1) then
        idistmax=60
      else
        distmax=distmax*60*0.5
        distmax=60./nint(distmax)
        idistmax=min(60,nint(60./distmax))
c       idistmax=idistmax*.75
        print *,'idistmax:',idistmax
c       if(idistmax.lt.60) idistmax=min(60,60/((60/idistmax)+1))
        if(idistmax.lt.60.and.idistmax.gt.0) then
          idistmax=min(60,60/((60/idistmax)+0))
        endif
      endif
      IF(NUMGBC.GT.0) THEN
        DO N=1,NUMGBC
          READ(1,1007) IBFLUX,IBFLUX,BFLUX
1007      FORMAT(2I6,E13.6)
        enddo
      ENDIF
      print *,'here'
      if(idistmax.eq.7) idistmax=6
      IF(IDISTMAX.EQ.60) THEN
c       WRITE(JUNK1,'(a)') '-I1/1'
c       WRITE(JUNK2,'(a)') '-S2'
        WRITE(JUNK1,903) 1,1
903     FORMAT('-I',I2,'/',I2)
c       WRITE(JUNK2,904) nint(2*dmax)
        WRITE(JUNK2,904) 2
904     FORMAT('-S',I2)
        CALL STRIP(80,JUNK1,N1)
        CALL STRIP(80,JUNK2,N2)
      ELSEIF(IDISTMAX.EQ.0) THEN
        WRITE(JUNK1,'(a)') '-I1m/1m'
        WRITE(JUNK2,'(a)') '-S2m'
        N1=7
        N2=4
      ELSE
        WRITE(JUNK1,905) IDISTMAX,IDISTMAX
905     FORMAT('-I',I2,'m/',I2,'m')
        WRITE(JUNK2,906) nint(4*dmax*60)
906     FORMAT('-S',I3,'m')
        CALL STRIP(80,JUNK1,N1)
        CALL STRIP(80,JUNK2,N2)
      endif
    
      DISTXY=sqrt((xmax-xmin)**2+(ymax-ymin)**2)
      call recpol(xmin,ymin,rlat1,rlong1)
      call recpol(xmax,ymin,rlat2,rlong2)
      call recpol(xmax,ymax,rlat3,rlong3)
      call recpol(xmin,ymax,rlat4,rlong4)
c     if(ihemi.eq.1) then
c       rlat1=-rlat1
c       rlat2=-rlat2
c       rlat3=-rlat3
c       rlat4=-rlat4
c       rlong1=90.-rlong1
c       rlong2=90.-rlong2
c       rlong3=90.-rlong3
c       rlong4=90.-rlong4
c       if(rlong1.lt.0.) rlong1=rlong1+360.
c       if(rlong2.lt.0.) rlong2=rlong2+360.
c       if(rlong3.lt.0.) rlong3=rlong3+360.
c       if(rlong4.lt.0.) rlong4=rlong4+360.
c     endif
      print *,rlat1,rlong1
      print *,rlat2,rlong2
      print *,rlat3,rlong3
      print *,rlat4,rlong4
      if(ihemi.eq.0) then
        if(jhemi.eq.0) then
c         WRITE(JUNK3,900) RLONG4,RLAT4,RLONG2,RLAT2
          WRITE(JUNK3,900) RLONG1,RLAT1,RLONG3,RLAT3
        else
c         WRITE(JUNK3,900) RLONG1,RLAT1,RLONG3,RLAT3
          WRITE(JUNK3,900) RLONG2,RLAT2,RLONG4,RLAT4
        endif
        WRITE(JUNK5,901) int(rlongmin),int(rlongmax),
     &                   int(rlatmin),int(rlatmax)
c       WRITE(JUNK5,901) 0,360,20,90
      else
        WRITE(JUNK3,900) RLONG1,RLAT1,RLONG3,RLAT3
c       WRITE(JUNK3,900) RLONG4,RLAT4,RLONG2,RLAT2
        WRITE(JUNK5,901) int(rlongmin),int(rlongmax),
     &                   int(rlatmin),int(rlatmax)
        WRITE(*,*) int(rlongmin),int(rlongmax),
     &                   int(rlatmin),int(rlatmax)
c       WRITE(JUNK5,901) 0,360,-90,-20
      endif
900   FORMAT('-R',F7.3,'/',F7.3,'/',F7.3,'/',F7.3)
901   FORMAT('-R',I4,'/',I4,'/',I4,'/',I4)
      WRITE(JUNK6,902) 0.5*(RLAT1+RLAT3)
c   902   FORMAT('-Lx6.3/0.0/',f6.2,'/500 ')
902   FORMAT('-Lx6.3/0.0/',f6.2,'/50 ')
      WRITE(JUNK4,'(i10)') nint(6*DISTXY)
      n=80
      CALL STRIP(80,JUNK3,N3)
      CALL STRIP(80,JUNK4,N4)
      CALL STRIP(80,JUNK5,N5)
      CALL STRIP(80,JUNK6,N6)
      print *,'j1',JUNK1
      print *,'j2',JUNK2
      print *,'j3',JUNK3
      print *,'j4',JUNK4
      print *,'j5',JUNK5
      print *,'j6',JUNK6
c--------------------------------------------------
10    format(a)
20    format(a/a)
30    format('#----------------------------------------------')
40    format('nearneighbor gmt.xyz -Gout.grd ',a,' ',a/
     &'             ',a,' ',a)
      io=12
      ia=13
      open(unit=io,file='coast.gmt')
      if(.true.) then !--------- COAST ---------------------
        write(io,30)
        if(ihemi.eq.0) then
          if(jhemi.eq.0) then
            write(io,500) JUNK3(1:n3),lend,JUNK4(1:n4),lend
          else
            write(io,501) JUNK3(1:n3),lend,JUNK4(1:n4),lend
          endif
        else
          write(io,51) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        endif
        write(io,20) 'psxy -R -Js -M outline.xy -W2/255/255/255 '//lend,
     &'             -K -O >> test.ps'
c--------------------------------------------------
        write(io,30)
        write(io,20) 'pscoast -A2000 -Di -Na/5 '//JUNK6(:N6)//lend,
     &  '            -U -R -Js -P -W2 -K -O  >> test.ps'
        write(io,20) 'pscoast -A2000 -Di -Na/5 '//JUNK6(:N6)//lend,
     &  '            -U -R -Js -P -W2 -K -O  -M > asdf.outline'
c--------------------------------------------------
        write(io,30)
      write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -K -O << END >> test.ps'
        write(io,10) '0  6.5 20 0 1 1  COASTLINE'
        write(io,10) 'END'
        write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -O  << END >> test.ps'
        write(io,10) '0.2  6.2 20 0 1 1  TIME'
        write(io,10) 'END'
c--------------------------------------------------
        write(io,30)
        write(io,10) 'cp test.ps Coast.ps'
        write(io,10) '#display Coast.ps &'
c       write(io,10) 'ps2img -crop a Coast.ps'
        write(io,10) 'convert Coast.ps Coast.gif'
        write(io,10) '#pngtogif Coast$2'
        write(io,10) '#ps2pdf Coast$2.ps'
        write(io,10) 'mapmake.x'

      endif
      close(io)
      open(unit=io,file='generic.gmt')
      open(unit=ia,file='generica.gmt')
      if(.true.) then !--------- THICK ---------------------
      write(io,10) 'gmtmake.e << END > zxcv'
      write(ia,10) 'gmtmake.e << END > zxcv'
      if(ihemi.eq.0) then
        write(io,10) '0'
        write(io,10) '6'
        write(io,10) '0 0'
        write(ia,10) '0'
        write(ia,10) '6'
        write(ia,10) '0 0'
      else
        write(io,10) '1'
        write(io,10) '6'
        write(ia,10) '1'
        write(ia,10) '6'
      endif
      write(io,10) 'END'
      write(io,40) JUNK1(1:N1),lend,JUNK5(1:N5),JUNK2(1:N2)
      write(ia,10) 'END'
      write(ia,10) 'cp gmt.xyz First.grd'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'grdmath out.grd 1 MUL = First.grd'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'gmtmake.e << END > zxcv'
      write(ia,30)
      write(ia,10) 'gmtmake.e << END > zxcv'
      if(ihemi.eq.0) then
        write(io,10) '0'
        write(io,10) '0'
        write(io,10) '0 0'
        write(ia,10) '0'
        write(ia,10) '0'
        write(ia,10) '0 0'
      else
        write(io,10) '1'
        write(io,10) '0'
        write(ia,10) '1'
        write(ia,10) '0'
      endif
      write(io,10) 'END'
      write(io,40) JUNK1(1:N1),lend,JUNK5(1:N5),JUNK2(1:N2)
      write(ia,10) 'END'
      write(ia,10) 'cp gmt.xyz Second.grd'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'grdmath out.grd 1 MUL = Second.grd'
c--------------------------------------------------
      write(io,30)
      write(ia,30)
      if(ihemi.eq.0) then
        if(jhemi.eq.0) then
          write(io,500) JUNK3(1:n3),lend,JUNK4(1:n4),lend
          write(ia,500) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        else
          write(io,501) JUNK3(1:n3),lend,JUNK4(1:n4),lend
          write(ia,501) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        endif
500     format('psbasemap -U ',a,'r ',a,/
     &'        -Js270/90.000/1:',a,' -B90g10/10g10 ',a,/
     &'        -P -X1 -Y3 -K  > test.ps')
501     format('psbasemap -U ',a,'r ',a,/
     &'        -Js0/90.000/1:',a,' -B90g10/10g10 ',a,/
     &'        -P -X1 -Y3 -K  > test.ps')
      else
        write(io,51) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        write(ia,51) JUNK3(1:n3),lend,JUNK4(1:n4),lend
51      format('psbasemap -U ',a,'r ',a,/
     &'        -Js0/-90.000/1:',a,' -B90g10/10g10 ',a,/
     &'        -P -X1 -Y3 -K  > test.ps')
      endif
c--------------------------------------------------
      write(io,30)
      write(io,10) 'makecpt -Z -T0/4500/500 -Crainbow > topo.cpt'
      write(io,10) '#grd2cpt First.grd -Z -Crainbow > topo.cpt'
      write(ia,30)
      write(ia,10) 'makecpt -T0/4600/200 -Crainbow > topo.cpt'
c--------------------------------------------------
      write(io,30)
      write(io,10) '#grdgradient First.grd -A45 -GFirstt.grd -Nt'
      write(io,20) '#grdimage First.grd -IFirstt.grd -R -Js '//lend,
     &'#                   -P -U -Ctopo.cpt -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'grdimage First.grd -R -Js -O -K -P -U '//lend,
     &'             -Ctopo.cpt  >> test.ps'
      write(io,20) 'grdcontour -S5 First.grd -Js -R -C1 -B '//lend,
     &'    -P -A1 -W5/255/255/255 -T -L0/1 -K -O >> test.ps'
      write(ia,30)
      write(ia,20) 'pscontour First.grd -A- -V  -I -Ctopo.cpt'//lend,
     &'       -Js -R -K -O >> test.ps'
      write(ia,20) '#pscontour First.grd -A- -V -W2/0/0/0 '//lend,
     &'#       -Ctopo.cpt -Js -R -B -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'psscale -E -Ctopo.cpt -D6.1/4/3.3/0.3 '//lend,
     &'            -L -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'psxy -R -Js -M outline.xy -W2/255/255/255 '//lend,
     &'             -K -O >> test.ps'
      write(ia,30)
      write(ia,20) 'psscale -E -Ctopo.cpt -D6.0/3.4/6.5/0.3 '//lend,
     &'            -L -K -O >> test.ps'
c--------------------------------------------------
      write(ia,30)
      write(ia,10) 'makecpt -T0/1/1 -Crainbow > topo.cpt'
      write(ia,20) 'pscontour First.grd -A- -V -W5/255/255/255 '//
     &              lend,
     &'       -Ctopo.cpt -Js -R -B -K -O >> test.ps'
c--------------------------------------------------
      write(ia,30)
      write(ia,20) 'psxy -R -Js -M outline.xy -W2/255/255/255 '//lend,
     &'             -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'pscoast -A2000 -Di -Na/5 '//JUNK6(:N6)//lend,
     &'            -U -R -Js -P -W5 -K -O  >> test.ps'
      write(ia,30)
      write(ia,20) 'pscoast -A2000 -Di -Na/5 '//JUNK6(:N6)//lend,
     &'            -U -R -Js -P -W5 -K -O  >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'makecpt -Z -T0/4500/250 -Crainbow > topo.cpt'
      write(io,10) '#grd2cpt Second.grd -Z -Crainbow > topo.cpt'
      write(ia,30)
      write(ia,10) 'makecpt -T0/4600/200 -Crainbow > topo.cpt'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'grdcontour -S5 Second.grd -Js -R -Ctopo.cpt'//lend,
     &'            -B -P -A -Wa -T -K -O >> test.ps'
      write(ia,30)
      write(ia,20) 'pscontour Second.grd -A- -V -W2/0/0/0 '//lend,
     &'       -Ctopo.cpt -Js -R -B -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -K -O << END >> test.ps'
      write(io,10) '0  6.5 20 0 1 1  THICKNESS [$2]'
      write(io,10) 'END'
      write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -O  << END >> test.ps'
      write(io,10) '0.2  6.2 20 0 1 1  TIME= $1'
      write(io,10) 'END'
      write(ia,30)
      write(ia,10) 'pstext -P -R0/11/0/8.5 -Jx1 -K -O << END >> test.ps'
      write(ia,10) '0  6.5 20 0 1 1  THICKNESS [$2]'
      write(ia,10) 'END'
      write(ia,10) 'psxy -R -Jx1 -W2/0/0/0 -Sp -K -O << END >> test.ps'
      write(ia,10) ' 7.5 7.0 '
      write(ia,10) 'END'
      write(ia,10) 'pstext -P -R0/11/0/8.5 -Jx1 -O  << END >> test.ps'
      write(ia,10) '0.2  6.2 20 0 1 1  TIME= $1'
      write(ia,10) 'END'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'cp test.ps Thick$2.ps'
      write(io,10) '#display Thick$2.ps &'
c     write(io,10) 'ps2img -crop a Thick$2.ps'
c     write(io,10) 'convert -trim Thick$2.ps Thick$2.gif'
      write(io,10) 'convert Thick$2.ps Thick$2.pdf'
      write(io,10) '#pngtogif Thick$2'
      write(io,10) '#ps2pdf Thick$2.ps'
      write(ia,30)
      write(ia,10) 'cp test.ps Thick$2.ps'
      write(ia,10) '#display Thick$2.ps &'
c     write(ia,10) 'ps2img -crop a Thick$2.ps'
c     write(ia,10) 'convert -trim Thick$2.ps Thick$2.gif'
      write(ia,10) 'convert Thick$2.ps Thick$2.pdf'
      write(ia,10) '#pngtogif Thick$2'
      write(ia,10) '#ps2pdf Thick$2.ps'
      endif !--------- THICK ---------------------
c--------------------------------------------------
      if(.true.) then !--------- VELO ---------------------
      write(io,10) 'gmtmake.e << END > zxcv'
      write(ia,10) 'gmtmake.e << END > zxcv'
      if(ihemi.eq.0) then
        write(io,10) '0'
        write(io,10) '15'
        write(io,10) '0 0'
        write(ia,10) '0'
        write(ia,10) '15'
        write(ia,10) '0 0'
      else
        write(io,10) '1'
        write(io,10) '15'
        write(ia,10) '1'
        write(ia,10) '15'
      endif
      write(io,10) 'END'
      write(ia,10) 'END'
      write(io,40) JUNK1(1:N1),lend,JUNK5(1:N5),JUNK2(1:N2)
      write(ia,10) 'cp gmt.xyz First.grd'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'grdmath out.grd 1 MUL = First.grd'
c     write(io,10) 'grdmath First.grd 1 ADD = First.grd'
c     write(io,10) 'grdmath First.grd LOG10 = First.grd'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'gmtmake.e << END > zxcv'
      write(ia,30)
      write(ia,10) 'gmtmake.e << END > zxcv'
      if(ihemi.eq.0) then
        write(io,10) '0'
        write(io,10) '0'
        write(io,10) '0 0'
        write(ia,10) '0'
        write(ia,10) '0'
        write(ia,10) '0 0'
      else
        write(io,10) '1'
        write(io,10) '0'
        write(ia,10) '1'
        write(ia,10) '0'
      endif
      write(io,10) 'END'
      write(ia,10) 'END'
      write(io,40) JUNK1(1:N1),lend,JUNK5(1:N5),JUNK2(1:N2)
      write(ia,10) 'cp gmt.xyz Second.grd'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'grdmath out.grd 1 MUL = Second.grd'
c--------------------------------------------------
      write(io,30)
      write(ia,30)
      if(ihemi.eq.0) then
        if(jhemi.eq.0) then
          write(io,500) JUNK3(1:n3),lend,JUNK4(1:n4),lend
          write(ia,500) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        else
          write(io,501) JUNK3(1:n3),lend,JUNK4(1:n4),lend
          write(ia,501) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        endif
      else
        write(io,51) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        write(ia,51) JUNK3(1:n3),lend,JUNK4(1:n4),lend
      endif
c--------------------------------------------------
      write(io,30)
      write(io,10) 'makecpt -Z -T0/1000/100 -Crainbow > topo.cpt'
      write(io,10) 'makecpt -Z -T0/4/0.25 -Crainbow > topo.cpt'
      write(io,10) '#grd2cpt First.grd -Z -Crainbow > topo.cpt'
      write(ia,30)
      write(ia,10) 'makecpt -T0/1000/50 -Crainbow > topo.cpt'
      write(ia,10) 'makecpt -T0/4/0.125 -Crainbow > topo.cpt'
c--------------------------------------------------
      write(io,30)
      write(io,10) '#grdgradient First.grd -A45 -GFirstt.grd -Nt'
      write(io,20) '#grdimage First.grd -IFirstt.grd -R -Js '//lend,
     &'#                   -P -U -Ctopo.cpt -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'grdimage First.grd -R -Js -O -K -P -U '//lend,
     &'             -Ctopo.cpt  >> test.ps'
      write(io,20) '#grdcontour -S5 First.grd -Js -R -C1 -B '//lend,
     &'#   -P -A1 -W5/255/255/255 -T -L0/1 -K -O >> test.ps'
      write(ia,30)
      write(ia,20) 'pscontour First.grd -A- -V  -I -Ctopo.cpt'//lend,
     &'       -Js -R -K -O >> test.ps'
      write(ia,20) '#pscontour First.grd -A- -V -W2/0/0/0 '//lend,
     &'#       -Ctopo.cpt -Js -R -B -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'psscale -E -Ctopo.cpt -D6.1/4/3.3/0.3 '//lend,
     &'            -L -K -O >> test.ps'
      write(io,20) 'psxy -R -Js -M outline.xy -W2/255/255/255 '//lend,
     &'             -K -O >> test.ps'
      write(ia,30)
      write(ia,20) 'psscale -E -Ctopo.cpt -D6.0/3.4/6.5/0.3 '//lend,
     &'            -L -K -O >> test.ps'
      write(ia,20) 'psxy -R -Js -M outline.xy -W2/255/255/255 '//lend,
     &'             -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'pscoast -A2000 -Di -Na/5 '//JUNK6(:N6)//lend,
     &'            -U -R -Js -P -W5 -K -O  >> test.ps'
      write(ia,30)
      write(ia,20) 'pscoast -A2000 -Di -Na/5 '//JUNK6(:N6)//lend,
     &'            -U -R -Js -P -W5 -K -O  >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'makecpt -Z -T0/4500/250 -Crainbow > topo.cpt'
      write(io,10) '#grd2cpt Second.grd -Z -Crainbow > topo.cpt'
      write(ia,30)
      write(ia,10) 'makecpt -T0/4600/200 -Crainbow > topo.cpt'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'grdcontour -S5 Second.grd -Js -R -Ctopo.cpt'//lend,
     &'            -B -P -A -Wa -T -K -O >> test.ps'
      write(ia,30)
      write(ia,20) 'pscontour Second.grd -A- -V -W2/0/0/0 '//lend,
     &'       -Ctopo.cpt -Js -R -B -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -K -O << END >> test.ps'
      write(io,10) '0  6.5 20 0 1 1  VELOCITY (log10) [$2]'
      write(io,10) 'END'
      write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -O  << END >> test.ps'
      write(io,10) '0.2  6.2 20 0 1 1  TIME= $1'
      write(io,10) 'END'
      write(ia,30)
      write(ia,10) 'pstext -P -R0/11/0/8.5 -Jx1 -K -O << END >> test.ps'
      write(ia,10) '0  6.5 20 0 1 1  VELOCITY (log10) [$2]'
      write(ia,10) 'END'
      write(ia,10) 'psxy -R -Jx1 -W2/0/0/0 -Sp -K -O << END >> test.ps'
      write(ia,10) ' 7.5 7.0 '
      write(ia,10) 'END'
      write(ia,10) 'pstext -P -R0/11/0/8.5 -Jx1 -O  << END >> test.ps'
      write(ia,10) '0.2  6.2 20 0 1 1  TIME= $1'
      write(ia,10) 'END'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'cp test.ps Velo$2.ps'
      write(io,10) '#display Velo$2.ps &'
c     write(io,10) 'ps2img -crop a Velo$2.ps'
c     write(io,10) 'convert -trim Velo$2.ps Velo$2.gif'
      write(io,10) 'convert Velo$2.ps Velo$2.pdf'
      write(io,10) '#pngtogif Velo$2'
      write(io,10) '#ps2pdf Velo$2.ps'
      write(ia,30)
      write(ia,10) 'cp test.ps Velo$2.ps'
      write(ia,10) '#display Velo$2.ps &'
c     write(ia,10) 'ps2img -crop a Velo$2.ps'
c     write(ia,10) 'convert -trim Velo$2.ps Velo$2.gif'
      write(ia,10) 'convert Velo$2.ps Velo$2.pdf'
      write(ia,10) '#pngtogif Velo$2'
      write(ia,10) '#ps2pdf Velo$2.ps'
      endif !--------- VELO ---------------------
      if(.true.) then !--------- WATER ---------------------
      write(io,10) 'gmtmake.e << END > zxcv'
      write(ia,10) 'gmtmake.e << END > zxcv'
      if(ihemi.eq.0) then
        write(io,10) '0'
        write(io,10) '18'
        write(io,10) '0 0'
        write(ia,10) '0'
        write(ia,10) '17'
        write(ia,10) '0 0'
      else
        write(io,10) '1'
        write(io,10) '18'
        write(ia,10) '1'
        write(ia,10) '17'
      endif
      write(io,10) 'END'
      write(io,40) JUNK1(1:N1),lend,JUNK5(1:N5),JUNK2(1:N2)
      write(ia,10) 'END'
      write(ia,10) 'cp gmt.xyz First.grd'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'grdmath out.grd 1 MUL = First.grd'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'gmtmake.e << END > zxcv'
      write(ia,30)
      write(ia,10) 'gmtmake.e << END > zxcv'
      if(ihemi.eq.0) then
        write(io,10) '0'
        write(io,10) '3'
        write(io,10) '0 0'
        write(ia,10) '0'
        write(ia,10) '0'
        write(ia,10) '0 0'
      else
        write(io,10) '1'
        write(io,10) '3'
        write(ia,10) '1'
        write(ia,10) '0'
      endif
      write(io,10) 'END'
      write(io,40) JUNK1(1:N1),lend,JUNK5(1:N5),JUNK2(1:N2)
      write(ia,10) 'END'
      write(ia,10) 'cp gmt.xyz Second.grd'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'grdmath out.grd 1 MUL = Second.grd'
c--------------------------------------------------
      write(io,30)
      if(ihemi.eq.0) then
        if(jhemi.eq.0) then
          write(io,500) JUNK3(1:n3),lend,JUNK4(1:n4),lend
          write(ia,500) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        else
          write(io,501) JUNK3(1:n3),lend,JUNK4(1:n4),lend
          write(ia,501) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        endif
      else
        write(io,51) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        write(ia,51) JUNK3(1:n3),lend,JUNK4(1:n4),lend
      endif
c--------------------------------------------------
      write(io,30)
      write(io,10) 'makecpt -Z -T0/20/2 -Crainbow > topo.cpt'
      write(io,10) '#grd2cpt First.grd -Z -Crainbow > topo.cpt'
      write(ia,30)
c     write(ia,10) 'makecpt -T0/20/1 -Crainbow > topo.cpt'
c     write(ia,10) 'makecpt -T0/0.5/0.025 -Crainbow > topo.cpt'
      write(ia,10) 'makecpt -T-0.0125/0.25/0.0125 -Crainbow > topo.cpt'
      write(ia,10) '#grd2cpt First.grd -Z -Crainbow > topo.cpt'
c--------------------------------------------------
      write(io,30)
      write(io,10) '#grdgradient First.grd -A45 -GFirstt.grd -Nt'
      write(io,20) '#grdimage First.grd -IFirstt.grd -R -Js '//lend,
     &'#                   -P -U -Ctopo.cpt -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'grdimage First.grd -R -Js -O -K -P -U '//lend,
     &'             -Ctopo.cpt  >> test.ps'
      write(io,20) '#grdcontour -S5 First.grd -Js -R -C1 -B '//lend,
     &'#   -P -A1 -W5/255/255/255 -T -L0/1 -K -O >> test.ps'
      write(ia,30)
      write(ia,20) 'pscontour First.grd -A- -V  -I -Ctopo.cpt'//lend,
     &'       -Js -R -K -O >> test.ps'
      write(ia,10) 'makecpt -T-.0125/.0125/.0125 -Crainbow > topo.cpt'
      write(ia,20) 'pscontour First.grd -A- -V -W2/0/0/0 '//lend,
     &'       -Ctopo.cpt -Js -R -B -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'psscale -E -Ctopo.cpt -D6.1/4/3.3/0.3 '//lend,
     &'            -L -K -O >> test.ps'
      write(io,20) 'psxy -R -Js -M outline.xy -W2/255/255/255 '//lend,
     &'             -K -O >> test.ps'
c--------------------------------------------------
      write(ia,30)
      write(ia,10) 'makecpt -T-0.0125/0.25/0.0125 -Crainbow > topo.cpt'
      write(ia,20) 'psscale -E -Ctopo.cpt -D6.0/3.4/6.5/0.3 '//lend,
     &'            -L -K -O >> test.ps'
      write(ia,30)
      write(ia,10) 'makecpt -T0/1/1 -Crainbow > topo.cpt'
      write(ia,20) 'pscontour First.grd -A- -V -W5/255/255/255 '//
     &              lend,
     &'       -Ctopo.cpt -Js -R -B -K -O >> test.ps'
      write(ia,30)
      write(ia,20) 'psxy -R -Js -M outline.xy -W2/255/255/255 '//lend,
     &'             -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'pscoast -A2000 -Di -Na/5 '//JUNK6(:N6)//lend,
     &'            -U -R -Js -P -W5 -K -O  >> test.ps'
      write(ia,30)
      write(ia,20) 'pscoast -A2000 -Di -Na/5 '//JUNK6(:N6)//lend,
     &'            -U -R -Js -P -W5 -K -O  >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'makecpt -Z -T0.1/0.9/0.4 -Crainbow > topo.cpt'
      write(io,10) '#grd2cpt Second.grd -Z -Crainbow > topo.cpt'
      write(ia,30)
      write(ia,10) 'makecpt -Z -T0/4600/200 -Crainbow > topo.cpt'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'grdcontour -S5 Second.grd -Js -R -Ctopo.cpt'//lend,
     &'            -B -P -A -Wa -T -K -O >> test.ps'
      write(ia,30)
      write(ia,20) 'pscontour Second.grd -A- -V -W2/0/0/0 '//lend,
     &'       -Ctopo.cpt -Js -R -B -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -K -O << END >> test.ps'
      write(io,10) '0  6.5 20 0 1 1  WATER PRODUCTION (mm/yr)[$2]'
      write(io,10) 'END'
      write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -O  << END >> test.ps'
      write(io,10) '0.2  6.2 20 0 1 1  TIME= $1'
      write(io,10) 'END'
      write(ia,30)
      write(ia,10) 'pstext -P -R0/11/0/8.5 -Jx1 -K -O << END >> test.ps'
      write(ia,10) '0  6.5 20 0 1 1  WATER THICKNESS (mm)[$2]'
      write(ia,10) 'END'
      write(ia,10) 'psxy -R -Jx1 -W2/0/0/0 -Sp -K -O << END >> test.ps'
      write(ia,10) ' 7.5 7.0 '
      write(ia,10) 'END'
      write(ia,10) 'pstext -P -R0/11/0/8.5 -Jx1 -O  << END >> test.ps'
      write(ia,10) '0.2  6.2 20 0 1 1  TIME= $1'
      write(ia,10) 'END'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'cp test.ps Water$2.ps'
      write(io,10) '#display Water$2.ps &'
c     write(io,10) 'ps2img -crop a Water$2.ps'
c     write(io,10) 'convert -trim Water$2.ps Water$2.gif'
      write(io,10) 'convert Water$2.ps Water$2.pdf'
      write(io,10) '#pngtogif Water$2'
      write(io,10) '#ps2pdf Water$2.ps'
      write(ia,30)
      write(ia,10) 'cp test.ps Water$2.ps'
      write(ia,10) '#display Water$2.ps &'
c     write(ia,10) 'ps2img -crop a Water$2.ps'
c     write(ia,10) 'convert -trim Water$2.ps Water$2.gif'
      write(ia,10) 'convert Water$2.ps Water$2.pdf'
      write(ia,10) '#pngtogif Water$2'
      write(ia,10) '#ps2pdf Water$2.ps'
      endif !--------- WATER ---------------------

      close(io)
      call system('chmod +x generic.gmt')
      call system('chmod +x generica.gmt')
      call system('chmod +x coast.gmt')
c     call system('generic.gmt')

      open(unit=io,file='generic-erosion.gmt')
      if(.false.) then !--------- EROSION ---------------------
      write(io,10) 'erosion.e $1 $2 << END > zxcv'
      write(io,10) '-100000 0'
      write(io,10) 'END'
      write(io,10) 'gmtmake.e << END > zxcv'
      if(ihemi.eq.0) then
        write(io,10) '0'
        write(io,10) '13'
        write(io,10) '0 0'
      else
        write(io,10) '1'
        write(io,10) '13'
      endif
      write(io,10) 'END'
      write(io,40) JUNK1(1:N1),lend,JUNK5(1:N5),JUNK2(1:N2)
c--------------------------------------------------
      write(io,30)
      write(io,10) 'grdmath out.grd 1 MUL = First.grd'
c     write(io,10) 'grdmath First.grd 1 ADD = First.grd'
c     write(io,10) 'grdmath First.grd LOG10 = First.grd'

c--------------------------------------------------
      write(io,10) 'gmtmake.e << END > zxcv'
      if(ihemi.eq.0) then
        write(io,10) '0'
        write(io,10) '1'
        write(io,10) '0 0'
      else
        write(io,10) '1'
        write(io,10) '1'
      endif
      write(io,10) 'END'
      write(io,40) JUNK1(1:N1),lend,JUNK5(1:N5),JUNK2(1:N2)
c--------------------------------------------------
      write(io,30)
      write(io,10) 'grdmath out.grd 1 MUL = Bed.grd'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'gmtmake.e << END > zxcv'
      if(ihemi.eq.0) then
        write(io,10) '0'
        write(io,10) '14'
        write(io,10) '0 0'
      else
        write(io,10) '1'
        write(io,10) '14'
      endif
      write(io,10) 'END'
      write(io,40) JUNK1(1:N1),lend,JUNK5(1:N5),JUNK2(1:N2)
c--------------------------------------------------
      write(io,30)
      write(io,10) 'grdmath out.grd 1e-3 MUL = Second.grd'
c--------------------------------------------------
      write(io,30)
      if(ihemi.eq.0) then
        if(jhemi.eq.0) then
          write(io,500) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        else
          write(io,501) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        endif
      else
        write(io,51) JUNK3(1:n3),lend,JUNK4(1:n4),lend
      endif
c--------------------------------------------------
      write(io,30)
      write(io,10) 'makecpt -Z -T0/5000/500 -Crainbow > topo.cpt'
      write(io,10) 'makecpt -T0/5000/250 -Crainbow > topo.cpt'
      write(io,10) '#makecpt -T0/1e6/5e4 -Crainbow > topo.cpt'
      write(io,10) 'makecpt -T0/6/0.25 -Crainbow > topo.cpt'
      write(io,10) '#grd2cpt First.grd -Z -Crainbow > topo.cpt'
c--------------------------------------------------
      write(io,30)
      write(io,10) '#grdgradient First.grd -A45 -GFirstt.grd -Nt'
      write(io,20) '#grdimage First.grd -IFirstt.grd -R -Js '//lend,
     &'#                   -P -U -Ctopo.cpt -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'grdimage First.grd -R -Js -O -K -P -U '//lend,
     &'             -Ctopo.cpt  >> test.ps'
      write(io,20) '#grdcontour -S5 First.grd -Js -R -Ctopo.cpt '//lend,
     &'#    -B -P -A- -Wa -T -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'psscale -E -Ctopo.cpt -D6.1/3.5/5.3/0.3 '//lend,
     &'            -L -K -O >> test.ps'
      write(io,20) 'psxy -R -Js -M outline.xy -W2/255/255/255 '//lend,
     &'             -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'pscoast -A2000 -Di -Na/5 '//JUNK6(:N6)//lend,
     &'            -U -R -Js -P -W5 -K -O  >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'makecpt -Z -T0/100/10 -Crainbow > topo.cpt'
      write(io,10) '#grd2cpt Second.grd -Z -Crainbow > topo.cpt'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'grdcontour -S5 Second.grd -Js -R -Ctopo.cpt'//lend,
     &'            -B -P -A -Wa/255/255/255 -T -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'makecpt -Z -T-1000/5000/100 -Crainbow > topo.cpt'
      write(io,10) '#grd2cpt Bed.grd -Z -Crainbow > topo.cpt'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'grdcontour -S5 Bed.grd -Js -R -Ctopo.cpt'//lend,
     &'            -B -P -A- -Wa -T -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -K -O << END >> test.ps'
      write(io,10) '0  6.5 20 0 1 1  EROSION [$2]'
      write(io,10) 'END'
      write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -O  << END >> test.ps'
      write(io,10) '0.2  6.2 20 0 1 1  (black:erosion,white:duration'
      write(io,10) 'END'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'cp test.ps Erosion-$2.ps'
      write(io,10) '#display Erosion-$2.ps &'
c     write(io,10) 'ps2img -crop a Erosion-$2.ps'
      write(io,10) 'convert -trim Erosion-$2.ps Erosion-$2.gif'
      write(io,10) 'convert -trim Erosion-$2.ps Erosion-$2.pdf'
      write(io,10) '#pngtogif Erosion-$2'
      write(io,10) 'prev Erosion-$2.pdf'
      endif !--------- EROSION ---------------------
      close(io)
      call system('chmod +x generic-erosion.gmt')
c     call system('generic-erosion.gmt a a')
c-----------------------------------------------------------
      io=13
      open(unit=io,file='any.gmt')
      if(.true.) then !--------- ANY ---------------------
      write(io,10) 'gmtmake.e << END > zxcv'
      if(ihemi.eq.0) then
        write(io,10) '0'
        write(io,10) '19'
        write(io,10) '0 0'
      else
        write(io,10) '1'
        write(io,10) '6'
      endif
      write(io,10) 'END'
      write(io,40) JUNK1(1:N1),lend,JUNK5(1:N5),JUNK2(1:N2)
c--------------------------------------------------
      write(io,30)
      write(io,10) 'grdmath out.grd 1 MUL = First.grd'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'gmtmake.e << END > zxcv'
      if(ihemi.eq.0) then
        write(io,10) '0'
        write(io,10) '15'
        write(io,10) '0 0'
      else
        write(io,10) '1'
        write(io,10) '0'
      endif
      write(io,10) 'END'
      write(io,40) JUNK1(1:N1),lend,JUNK5(1:N5),JUNK2(1:N2)
c--------------------------------------------------
      write(io,30)
      write(io,10) 'grdmath out.grd 1 MUL = Second.grd'
c--------------------------------------------------
      write(io,30)
      if(ihemi.eq.0) then
        if(jhemi.eq.0) then
          write(io,500) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        else
          write(io,501) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        endif
      else
        write(io,51) JUNK3(1:n3),lend,JUNK4(1:n4),lend
      endif
c--------------------------------------------------
      write(io,30)
      write(io,10) 'makecpt -Z -T0/4500/500 -Crainbow > topo.cpt'
      write(io,10) 'grd2cpt First.grd -Z -Crainbow > topo.cpt'
c--------------------------------------------------
      write(io,30)
      write(io,10) '#grdgradient First.grd -A45 -GFirstt.grd -Nt'
      write(io,20) '#grdimage First.grd -IFirstt.grd -R -Js '//lend,
     &'#                   -P -U -Ctopo.cpt -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'grdimage First.grd -R -Js -O -K -P -U '//lend,
     &'             -Ctopo.cpt  >> test.ps'
      write(io,20) 'grdcontour -S5 First.grd -Js -R -C1 -B '//lend,
     &'    -P -A1 -W5/255/255/255 -T -L0/1 -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'psscale -E -Ctopo.cpt -D6.1/4/3.3/0.3 '//lend,
     &'            -L -K -O >> test.ps'
      write(io,20) 'psxy -R -Js -M outline.xy -W2/255/255/255 '//lend,
     &'             -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'pscoast -A2000 -Di -Na/5 '//JUNK6(:N6)//lend,
     &'            -U -R -Js -P -W5 -K -O  >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'makecpt -Z -T0/4500/250 -Crainbow > topo.cpt'
      write(io,10) 'grd2cpt Second.grd -Z -Crainbow > topo.cpt'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'grdcontour -S5 Second.grd -Js -R -Ctopo.cpt'//lend,
     &'            -B -P -A -Wa -T -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -K -O << END >> test.ps'
      write(io,10) '0  6.5 20 0 1 1  SURFACE [$2]'
      write(io,10) 'END'
      write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -O  << END >> test.ps'
      write(io,10) '0.2  6.2 20 0 1 1  TIME= $1'
      write(io,10) 'END'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'cp test.ps Any$2.ps'
      write(io,10) '#display Any.ps &'
c     write(io,10) 'ps2img -crop a Any$2.ps'
      write(io,10) 'convert -trim Any$2.ps Any$2.gif'
      write(io,10) '#pngtogif Any$2'
      write(io,10) '#ps2pdf Any$2.ps'
      close(io)
      call system('chmod +x any.gmt')
      endif !--------- ANY ---------------------
c--------------------------------------------------
      end
c=========================================================
      SUBROUTINE SETRIG1
c     IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD                               
      call setproj1(90.d0); return
c---------------------------------------------
      RADIUS=6370.D0                                                    
      PI=4.D0*ATAN(1.D0)                                                
      RADIUS=2.0D4/PI                                              
      CIRCUM=2.D0*PI*RADIUS                                             
      RKMPDEG=CIRCUM/360.D0                                             
      RADPDEG=PI/180.D0                                                 
      DEGPRAD=180.D0/PI                                                 
      END                                                               
c=========================================================
      SUBROUTINE POLREC1(RLAT,RLONG,X,Y)                                 
c     IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD                               
      real*8 dlat,dlong,xee,ynn
      dlat=dble(rlat);dlong=dble(rlong)
      call ll2xyB1(dlat,dlong,xee,ynn)
      x=real(xee);y=real(ynn)
c     print *,rlat,rlong,x,y
      return
c---------------------------------------------
      X=1000.D0*(90.D0-RLAT)*RKMPDEG*COS(RLONG*RADPDEG)                 
      Y=1000.D0*(90.D0-RLAT)*RKMPDEG*SIN(RLONG*RADPDEG)                 
      END                                                               
c=========================================================
      SUBROUTINE RECPOL1(X,Y,RLAT,RLONG)                                 
c     IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD                               
      real*8 dlat,dlong,xee,ynn
      xee=dble(x);ynn=dble(y)
      call xy2llB1(xee,ynn,dlat,dlong)
      rlat=real(dlat);rlong=real(dlong);return
c---------------------------------------------
      IF(X.EQ.0.D0) THEN                                                
        IF(Y.GE.0.D0) THEN                                              
          THETA=90.D0                                                   
        ELSE                                                            
          THETA=270.D0                                                  
        ENDIF                                                           
      ELSE                                                              
        THETA=DEGPRAD*ATAN(Y/X)                                         
        IF(X.LT.0.D0) THEN                                              
          THETA=180.D0+THETA                                            
        ELSEIF(Y.LT.0.D0) THEN                                          
          THETA=360.D0+THETA                                            
        ELSEIF(THETA.GT.360.D0) THEN                                    
          THETA=THETA-360.D0                                            
        ENDIF                                                           
      ENDIF                                                             
      R=.001D0*SQRT(X**2+Y**2)                                          
      RLAT=90.D0-R/RKMPDEG                                              
      RLONG=THETA                                                       
      END                                                               
c=========================================================
      SUBROUTINE STRIP(N,JUNK,IC)
      CHARACTER*80 JUNK,JUNKN
      JUNKN=' '
      IC=0
      DO I=1,N
        IF(JUNK(I:I).NE.' ') THEN
          IC=IC+1
          JUNKN(IC:IC)=JUNK(I:I)
        ENDIF
      ENDDO
      JUNK=JUNKN
c      PRINT *,IC,JUNKN
      END
c======================================================
      subroutine setproj1(rlat0)
      implicit real*8(a-h,o-z)
      common /projcons/ aaa,eee,pi,radpdeg,rlong0,rlatf,phi0,
     &                   rk0,fe,fn,io
      io=0
      phi0=rlat0
      aaa=6378137.d0; eee=0.081819191d0;
      pi=4.d0*atan(1.d0); radpdeg=pi/180.d0
      fe=0d6; fn=0d6
      rlong0= 270.d0; rlatf= 71.d0; rk0=0.972769012891835d0
      end
c======================================================
      subroutine ll2xyB1(rlat,rlong,xee,ynn)
      implicit real*8(a-h,o-z)
      common /projcons/ aaa,eee,pi,radpdeg,rlong0,rlatf,phi0,
     &                   rk0,fe,fn,io
      logical npole
      npole =(phi0.eq.90d0); rlamb0=rlong0*radpdeg
      phif  =sign(rlatf*radpdeg,phi0)
      phi   =rlat  *radpdeg; rlamb =rlong *radpdeg
      term1=1d0+eee*sin(phif);term2=1d0-eee*sin(phif)
      denom1=(term1/term2)**(eee/2)
      term1=1d0+eee*sin(phi);term2=1d0-eee*sin(phi)
      denom2=(term1/term2)**(eee/2)
      if(npole) then
        tf= tan(pi/4d0-phif/2d0)*denom1
        ttt=tan(pi/4d0- phi/2d0)*denom2
      else
        tf= tan(pi/4d0+phif/2d0)/denom1
        ttt=tan(pi/4d0+ phi/2d0)/denom2
      endif;
      rmf=cos(phif)/sqrt(1d0-eee**2*sin(phif)**2)
      denom3=sqrt((1+eee)**(1+eee)*(1-eee)**(1-eee))
      rk0=rmf*denom3/(2d0*tf)
      rho=2d0*aaa*rk0*ttt/denom3
      theta=rlamb-rlamb0
      if(npole) then
        de=rho*sin(theta); dn=-rho*cos(theta)
      else
        de=rho*sin(theta); dn=rho*cos(theta)
      endif
      xee=fe+de; ynn=fn+dn
      end
c======================================================
      subroutine xy2llB1(xee,ynn,rlat,rlong)
      implicit real*8(a-h,o-z)
      common /projcons/ aaa,eee,pi,radpdeg,rlong0,rlatf,phi0,
     &                   rk0,fe,fn,io
      logical npole
      phif =sign(rlatf*radpdeg,phi0)
      npole=(phi0.eq.90d0); rlamb0=rlong0*radpdeg
      term1=1d0+eee*sin(phif);term2=1d0-eee*sin(phif)
      denom1=(term1/term2)**(eee/2)
      if(npole) then
        tf=tan(pi/4d0-phif/2d0)*denom1
      else
        tf=tan(pi/4d0+phif/2d0)/denom1
      endif
      rmf=cos(phif)/sqrt(1d0-eee**2*sin(phif)**2)
      denom3=sqrt((1+eee)**(1+eee)*(1-eee)**(1-eee))
      rk0=rmf*denom3/(2d0*tf)
      rho=sqrt((xee-fe)**2+(ynn-fn)**2)
      ttt=rho*denom3/(2d0*aaa*rk0)
      if(npole) then
        chi=pi/2d0-2d0*atan(ttt)
      else
        chi=2d0*atan(ttt)-pi/2d0
      endif 
      term1=eee**2/2d0+5d0*eee**4/24d0+eee**6/12d0+13d0*eee**8/360d0
      term2=7d0*eee**4/48d0+29d0*eee**6/240d0+811d0*eee**8/11520d0
      term3=7d0*eee**6+81d0*eee**8/1120d0
      term4=4279d0*eee**8/161280d0
      phi=chi+term1*sin(2d0*chi)+term2*sin(4d0*chi)+
     &        term3*sin(6d0*chi)+term4*sin(8d0*chi)
      if(xee.eq.fe) then
        rlamb=rlamb0
      elseif(npole) then
        rlamb=rlamb0+atan2((xee-fe),(fn-ynn))
      else
        rlamb=rlamb0+atan2((xee-fe),(ynn-fn))
      endif
      rlat=phi/radpdeg; rlong=rlamb/radpdeg
      if(rlong.lt.0) rlong=rlong+360d0
      if(rlong.gt.360) rlong=rlong-360d0
      end
c-----------------------------------------
      include "projections.f"
c-----------------------------------------
