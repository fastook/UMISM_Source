      character*80 line
      npixx=6848
      npixy=5696
      if(.false.) then
        xxmin=-2725
        yymin=-2275
        xxmax= 2725
        yymax= 2275
      elseif(.false.) then
c ..... direct from surface.asc (works best...)
        xxmin=-2713.6
        yymin=-2304.0
        xxmax= 2761.4
        yymax= 2246
      elseif(.true.) then
c ..... with 50 km added to LLY of surface.asc
        xxmin=-2713.6
        yymin=-2304.0+5
        xxmax= 2761.4
        yymax= 2246+5
      endif
      htot=xxmax-xxmin
      vtot=yymax-yymin
      xmin=-1400
      ymin=-1200
      xmax= -200
      ymax=    0
      print *,'input xmin,xmax,ymin,ymax'
      read(*,*) xmin,xmax,ymin,ymax
      xmin=max(xmin,xxmin)
      xmax=min(xmax,xxmax)
      ymin=max(ymin,yymin)
      ymax=min(ymax,yymax)
      dx=xmax-xmin
      dy=ymax-ymin
      if(dx.gt.dy) then
        ymax=ymin+dx
      else
        xmax=xmin+dy
      endif
      hrange=xmax-xmin
      vrange=ymax-ymin
      if(.false.) then
        mulpixx=npixx*((xmin-xxmin)/htot)
        mulpixy=npixy*((yymax-ymax)/vtot)
        mpixx=(hrange/htot)*npixx
        mpixy=(vrange/vtot)*npixy
        print *,mulpixx,mulpixy,mpixx,mpixy
        write(11,200) 'rm antloc.tif'
        write(11,200) 'ln -s ',
     &                '../ext-data/ant4.tif ',
     &                'antloc.tif'
        write(line,100) 'imgcopy -o ',mulpixx,',',mulpixy,
     &                  ' -s ',mpixx,',',mpixy,
     &                  ' antloc.tif asdf.tif'
      else
        mulpixx=npixx*((xmin-xxmin)/htot)
        mulpixy=npixy*((yymax-ymax)/vtot)
        mpixx=(hrange/htot)*npixx
        mpixy=(vrange/vtot)*npixy
        print *,mulpixx,mulpixy,mpixx,mpixy
c       mleft=mulpixx
c       mtop =mulpixy
c       mright=npixx-(mleft+mpixx)
c       mbottom=npixy-(mtop+mpixy)

        write(11,200) 'ln -s ',
     &                '../ext-data/ant4.pnm ',
     &                'antloc.pnm'
        write(line,100) 'pnmcut -left ',mulpixx,
     &                        ' -width ',mpixx,
     &                        ' -top ',mulpixy,
     &                        ' -height ',mpixy,
     &                  ' antloc.pnm > asdf.pnm'
      endif
      icomma=index(line,', ')
      print *,line
      dowhile(icomma.ne.0)
        line=line(:icomma)//line(icomma+2:)
        icomma=index(line,', ')
        print *,line
      enddo

      write(11,100) line
      if(.false.) then
c       write(11,100) 'imgcopy -Cgrey asdf.tif asdf.ras'
        write(11,100) 'imgcopy asdf.tif asdf.ras'
        write(11,100) 'xv asdf.ras'
      else
        write(11,100) 'xv asdf.pnm'
      endif
100   format(1x,a,i4,a,i4,a,i4,a,i4,a)
200   format(1x,3a)
      end

