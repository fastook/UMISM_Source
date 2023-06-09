      PARAMETER(NMAX=40000)
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
      call setrig
      lend=char(92)
c     open(1,file='inputnamer.data')
      READ(1,1000) HED,NUMNP,NUMEL,NUMCOL,NUMLEV,NUMGBC,NDT,INTER,DT
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
        READ(1,1001) NUM,KODE,X(n),Y(n),HTICE,
     &                   ADOT,FRACT,PSURF,
     &                   BDROCK,FLOWA,SLDGB,
     &                   temp,itype,AFUDGE,GEOFLUX,
     &                   calv
        call recpol(1*x(n),1*y(n),
     &                rlatn,rlongn)
        if(ihemi.eq.1) then
          rlatn=-rlatn
c         rlongn=90.-rlongn
          if(rlongn.lt.0.) rlongn=rlongn+360.
        endif
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
      if(ihemi.eq.1) then
        call swap(xmin,ymin) 
        call swap(xmax,ymax) 
      endif
      print *,'X   ',xmin,xmax
      print *,'Y   ',ymin,ymax
      print *,'Z   ',zmin,zmax
      print *,'LAT ',rlatmin,rlatmax
      print *,'LONG',rlongmin,rlongmax
c     rlatmin=max(-90.,rlatmin-1)
c     rlatmax=min(90.,rlatmax+1)
c     rlongmin=max(0.,rlongmin-1)
c     rlongmax=min(360.,rlongmax+1)
      rlatmin=max(-90.,rlatmin)
      rlatmax=min(90.,rlatmax)
      rlongmin=max(0.,rlongmin)
      rlongmax=min(360.,rlongmax)
      DISTMAX=0
      DISTMIN=1E30
      DO N=1,NUMEL
        READ(1,1002) NUM,(KX(n,i),i=1,4),CONST,ACON
        do i=1,4
          call recpol(1*x(kx(n,i)),1*y(kx(n,i)),
     &                rlate(i),rlonge(i))
          if(ihemi.eq.1) then
            rlate(i)=-rlate(i)
c           rlonge(i)=90.-rlonge(i)
            if(rlonge(i).lt.0.) rlonge(i)=rlonge(i)+360.
          endif
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
      if(ihemi.eq.0) then
        call recpol(xmin,ymin,rlat1,rlong1)
        call recpol(xmax,ymin,rlat2,rlong2)
        call recpol(xmax,ymax,rlat3,rlong3)
        call recpol(xmin,ymax,rlat4,rlong4)
      elseif(ihemi.eq.1) then
        call recpol(xmin,-ymin,rlat1,rlong1)
        call recpol(xmax,-ymin,rlat2,rlong2)
        call recpol(xmax,-ymax,rlat3,rlong3)
        call recpol(xmin,-ymax,rlat4,rlong4)
        rlat1=-rlat1
        rlat2=-rlat2
        rlat3=-rlat3
        rlat4=-rlat4
c       rlong1=90.-rlong1
c       rlong2=90.-rlong2
c       rlong3=90.-rlong3
c       rlong4=90.-rlong4
        if(rlong1.lt.0.) rlong1=rlong1+360.
        if(rlong2.lt.0.) rlong2=rlong2+360.
        if(rlong3.lt.0.) rlong3=rlong3+360.
        if(rlong4.lt.0.) rlong4=rlong4+360.
      else
        stop
      endif
      print *,rlat1,rlong1
      print *,rlat2,rlong2
      print *,rlat3,rlong3
      print *,rlat4,rlong4
      if(ihemi.eq.0) then
        if(jhemi.eq.0) then
          WRITE(JUNK3,900) RLONG1,RLAT1,RLONG3,RLAT3
        else
          WRITE(JUNK3,900) RLONG2,RLAT2,RLONG4,RLAT4
        endif
        WRITE(JUNK5,901) int(rlongmin),int(rlongmax),
     &                   int(rlatmin),int(rlatmax)
c       WRITE(JUNK5,901) 0,360,20,90
      else
        WRITE(JUNK3,900) RLONG4,RLAT4,RLONG2,RLAT2
        WRITE(JUNK5,901) int(rlongmin),int(rlongmax),
     &                   int(rlatmin),int(rlatmax)
        WRITE(*,*) int(rlongmin),int(rlongmax),
     &                   int(rlatmin),int(rlatmax)
c       WRITE(JUNK5,901) 0,360,-90,-20
      endif
900   FORMAT('-R',F7.3,'/',F7.3,'/',F7.3,'/',F7.3)
901   FORMAT('-R',I4,'/',I4,'/',I4,'/',I4)
      WRITE(JUNK6,902) 0.5*(RLAT1+RLAT3)
902   FORMAT('-Lx6.3/0.0/',f6.2,'/500 ')
      WRITE(JUNK4,'(i10)') nint(6*DISTXY)
      n=80
      CALL STRIP(80,JUNK3,N3)
      CALL STRIP(80,JUNK4,N4)
      CALL STRIP(80,JUNK5,N5)
      CALL STRIP(80,JUNK6,N6)
      print *,JUNK1
      print *,JUNK2
      print *,JUNK3
      print *,JUNK4
      print *,JUNK5
      print *,JUNK6
c--------------------------------------------------
10    format(a)
20    format(a/a)
30    format('#----------------------------------------------')
40    format('nearneighbor gmt.xyz -Gout.grd ',a,' ',a/
     &'             ',a,' ',a)
      io=12
      close(io)
      open(unit=io,file='generic2.gmt')
      if(.true.) then !--------- THICK ---------------------
      write(io,10) 'gmtmake.e << END > zxcv'
      if(ihemi.eq.0) then
        write(io,10) '0'
        write(io,10) '6'
        write(io,10) '0 0'
      else
        write(io,10) '1'
        write(io,10) '6'
      endif
      write(io,10) 'END'
      write(io,10) 'cp gmt.xyz First.xyz'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'gmtmake.e << END > zxcv'
      if(ihemi.eq.0) then
        write(io,10) '0'
        write(io,10) '0'
        write(io,10) '0 0'
      else
        write(io,10) '1'
        write(io,10) '0'
      endif
      write(io,10) 'END'
      write(io,10) 'cp gmt.xyz Second.xyz'
c--------------------------------------------------
      write(io,30)
      if(ihemi.eq.0) then
        if(jhemi.eq.0) then
          write(io,500) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        else
          write(io,501) JUNK3(1:n3),lend,JUNK4(1:n4),lend
        endif
500     format('psbasemap -U ',a,'r ',a,/
     &'        -Js270/90.000/1:',a,' -B90g10/10g10 ',a,/
     &'        -P -X1 -Y3 -K  > test.ps')
501     format('psbasemap -U ',a,'r ',a,/
     &'        -Js0/90.000/1:',a,' -B90g10/10g10 ',a,/
     &'        -P -X1 -Y3 -K  > test.ps')
      else
        write(io,51) JUNK3(1:n3),lend,JUNK4(1:n4),lend
51      format('psbasemap -U ',a,'r ',a,/
     &'        -Js0/-90.000/1:',a,' -B90g10/10g10 ',a,/
     &'        -P -X1 -Y3 -K  > test.ps')
      endif
c--------------------------------------------------
      write(io,30)
      write(io,10) 'makecpt -T0/4500/250 -Crainbow > topo.cpt'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'pscontour First.xyz -A- -V -I -Ctopo.cpt'//lend,
     &'                    -Js -R -K -O >> test.ps'
      write(io,10) '#pscontour First.xyz -A- -V -Ctopo.cpt'//lend,
     &'#      -W1/0/0/0 -B -Js -R -K -O >> test.ps'
      write(io,10) '#pscontour First.xyz -A- -V -Ctopo.cpt'//lend,
     &'#      -W1/255/255/255 -L2/0/0/0 -B -Js -R -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'psscale -E -Ctopo.cpt -D6.1/4/3.3/0.3 '//lend,
     &'            -L -K -O >> test.ps'
      write(io,20) 'psxy -R -Js -M outline.xy -W2/255/255/255 '//lend,
     &'             -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'pscoast -A2000 -Di -Na/5 '//JUNK6(:N6)//lend,
     &'            -U -R -Js -P -W2 -K -O  >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'makecpt -T0/4500/250 -Crainbow > topo.cpt'
c--------------------------------------------------
      write(io,30)
      write(io,10) '#pscontour Second.xyz -A- -V -I -Ctopo.cpt'//lend,
     &'#                    -Js -R -K -O >> test.ps'
      write(io,10) 'pscontour Second.xyz -A- -V -Ctopo.cpt'//lend,
     &'      -W1/0/0/0 -B -Js -R -K -O >> test.ps'
      write(io,10) '#pscontour Second.xyz -A- -V -Ctopo.cpt'//lend,
     &'#      -W1/255/255/255 -L2/0/0/0 -B -Js -R -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -K -O << END >> test.ps'
      write(io,10) '0  6.5 20 0 1 1  THICKNESS [$2]'
      write(io,10) 'END'
      write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -O  << END >> test.ps'
      write(io,10) '0.2  6.2 10 0 1 1  TIME= $1'
      write(io,10) 'END'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'cp test.ps Thick$2.ps'
      write(io,10) '#display Thick.ps &'
      write(io,10) 'convert Thick$2.ps Thick$2.gif'
      write(io,10) 'convert Thick$2.ps Thick$2.pdf'
      endif !--------- THICK ---------------------
c--------------------------------------------------
      if(.true.) then !--------- VELOCITY ---------------------
      write(io,10) 'gmtmake.e << END > zxcv'
      if(ihemi.eq.0) then
        write(io,10) '0'
        write(io,10) '15'
        write(io,10) '0 0'
      else
        write(io,10) '1'
        write(io,10) '15'
      endif
      write(io,10) 'END'
c     write(io,10) 'xylogz.x'
      write(io,10) 'cp gmt.xyz First.xyz'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'gmtmake.e << END > zxcv'
      if(ihemi.eq.0) then
        write(io,10) '0'
        write(io,10) '0'
        write(io,10) '0 0'
      else
        write(io,10) '1'
        write(io,10) '0'
      endif
      write(io,10) 'END'
      write(io,10) 'cp gmt.xyz Second.xyz'
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
      write(io,10) 'makecpt -T0/4/0.25 -Crainbow > topo.cpt'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'pscontour First.xyz -A- -V -I -Ctopo.cpt'//lend,
     &'                    -Js -R -K -O >> test.ps'
      write(io,10) '#pscontour First.xyz -A- -V -Ctopo.cpt'//lend,
     &'#      -W1/0/0/0 -B -Js -R -K -O >> test.ps'
      write(io,10) '#pscontour First.xyz -A- -V -Ctopo.cpt'//lend,
     &'#      -W1/255/255/255 -L2/0/0/0 -B -Js -R -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'psscale -E -Ctopo.cpt -D6.1/4/3.3/0.3 '//lend,
     &'            -L -K -O >> test.ps'
      write(io,20) 'psxy -R -Js -M outline.xy -W2/255/255/255 '//lend,
     &'             -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,20) 'pscoast -A2000 -Di -Na/5 '//JUNK6(:N6)//lend,
     &'            -U -R -Js -P -W2 -K -O  >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'makecpt -T0/4500/250 -Crainbow > topo.cpt'
c--------------------------------------------------
      write(io,30)
      write(io,10) '#pscontour Second.xyz -A- -V -I -Ctopo.cpt'//lend,
     &'#                    -Js -R -K -O >> test.ps'
      write(io,10) 'pscontour Second.xyz -A- -V -Ctopo.cpt'//lend,
     &'      -W1/0/0/0 -B -Js -R -K -O >> test.ps'
      write(io,10) '#pscontour Second.xyz -A- -V -Ctopo.cpt'//lend,
     &'#      -W1/255/255/255 -L2/0/0/0 -B -Js -R -K -O >> test.ps'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -K -O << END >> test.ps'
      write(io,10) '0  6.5 20 0 1 1  VELOCITY [$2]'
      write(io,10) 'END'
      write(io,10) 'pstext -P -R0/11/0/8.5 -Jx1 -O  << END >> test.ps'
      write(io,10) '0.2  6.2 10 0 1 1  TIME= $1'
      write(io,10) 'END'
c--------------------------------------------------
      write(io,30)
      write(io,10) 'cp test.ps Velo$2.ps'
      write(io,10) '#display Velo.ps &'
      write(io,10) 'convert Velo$2.ps Velo$2.gif'
      write(io,10) 'convert Velo$2.ps Velo$2.pdf'
      endif !--------- VELOCITY ---------------------
      close(io)
      call system('chmod +x generic2.gmt')
c     call system('generic2.gmt')
      end

      SUBROUTINE SETRIG                                                 
c     IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD                               
      RADIUS=6370.D0                                                    
      PI=4.D0*ATAN(1.D0)                                                
      RADIUS=2.0D4/PI                                              
      CIRCUM=2.D0*PI*RADIUS                                             
      RKMPDEG=CIRCUM/360.D0                                             
      RADPDEG=PI/180.D0                                                 
      DEGPRAD=180.D0/PI                                                 
      END                                                               
      SUBROUTINE POLREC(RLAT,RLONG,X,Y)                                 
c     IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD                               
      X=1000.D0*(90.D0-RLAT)*RKMPDEG*COS(RLONG*RADPDEG)                 
      Y=1000.D0*(90.D0-RLAT)*RKMPDEG*SIN(RLONG*RADPDEG)                 
      END                                                               
      SUBROUTINE RECPOL(X,Y,RLAT,RLONG)                                 
c     IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRIG/RKMPDEG,RADPDEG,DEGPRAD                               
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
C******************************
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
C******************************
      subroutine swap(a,b)
      t=a
      a=b
      b=t
      end
