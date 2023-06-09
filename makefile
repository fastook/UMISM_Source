SDIR	=	..
COMP	=	gfortran
G77	=	$(COMP)
#G77	=	g77 -bind_at_load
#CC	=	ifc -w 
#CFLAGS	=	-O1 -xW
CC	=	$(COMP)
#CC	=	g77 -bind_at_load
#CFLAGS	=	-O1 -ffortran-bounds-check
CFLAGS	=	-O3
#CFLAGS	=	-g -O0 -C
#GFLAGS	=	-g -O0 -C
#GFLAGS	=	-g -O0 -C
GFLAGS	=	-O1 
OFLAGS	=	-O3
FFLAGS	=	-O1 -xW
THREED	=	threedstub
PROJS	=	projections.o
DOBJS5	=	accum.o gausinit.o nconst.o readn.o shape.o \
		elprop.o nodesl.o scale.o volume.o afunct.o \
		plotsol.o setrig.o zoom.o formc.o   \
		gauseid.o modify.o gnumbr.o slopef.o temper.o adjust.o \
		$(THREED).o wthikn.o wthik.o wthikj.o oplate.o veplate.o splate.o \
		jcg.o lookup.o wthiksimp.o
DSRCS5	=	accum.f gausinit.f nconst.f readn.f shape.f \
		elprop.f nodesl.f scale.f volume.f afunct.f \
		plotsol.f setrig.f zoom.f formc.f   \
		gauseid.f modify.f gnumbr.f slopef.f temper.f adjust.f \
		$(THREED).f wthikn.f wthik.f wthikj.f oplate.f veplate.f splate.f \
		jcg.f lookup.f wthiksimp.f
DOBJS6	=	accum.o nconst.o readn.o shape.o \
		elprop.o nodesl.o scale.o volume.o afunct.o \
		plotsol.o setrig.o zoom.o formc.o   \
		gauseid.o modify.o gnumbr.o slopef.o temper.o adjust.o
DSRCS6	=	accum.f nconst.f readn.f shape.f \
		elprop.f nodesl.f scale.f volume.f afunct.f \
		plotsol.f setrig.f zoom.f formc.f  \
		gauseid.f modify.f gnumbr.f slopef.f temper.f 
DOBJS7	=	accum.o gausinit.o nconst.o readn.o shape.o \
		elprop.o nodesl.o scale.o volume.o afunct.o \
		plotsol.o setrig.o zoom.o formc.o   \
		gauseid.o modify.o gnumbr.o slopef.o temper.o adjust.o \
		$(THREED).o wthik.o veplate.o jcg.o
DSRCS7	=	accum.f gausinit.f nconst.f readn.f shape.f \
		elprop.f nodesl.f scale.f volume.f afunct.f \
		plotsol.f setrig.f zoom.f formc.f   \
		gauseid.f modify.f gnumbr.f slopef.f temper.f adjust.f \
		$(THREED).f wthik.f veplate.f jcg.f
#GLIBS   =	../../grstrt/grstrt.o -lfgl -lgl
#GLIBS    =	../../grstrt/grstrt.o ../../grstrt/scrsv.o -L/usr/X11R6/lib -lX11 -lGL -lm 
# STANDARD GRAPHICS LIBS, THE "FIX" FOR X11 "BUG" -------------------------
LDFLAGS =       "-Wl,-dylib_file,/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib"
X11DIR  =       /usr/X11
#GOBJS   =       ../../grstrt/grstrt.o ../../grstrt/scrsv.o
GOBJS   =       ~/grstrt/grstrt.o ~/grstrt/scrsv.o
GLLIBS   =      -L$(X11DIR)/lib -lGL -lX11 -I$(X11DIR)/include
GLIBS    =      $(GOBJS) $(GLLIBS) $(LDFLAGS)
# STANDARD GRAPHICS LIBS, THE "FIX" FOR X11 "BUG" -------------------------
#SLIBS   =	../../grstrt/scrsv.o -lgutil -limage -s
SLIBS   =	
LIBS    =	-lm
OGLLIBS   =	-L/usr/X11R6/lib -lXi -lXmu -lXext -lglut -lGLU -lX11 -lGL -lm -I/usr/X11R6/include
#NLIBS   =	-I/sw/include  /sw/lib/libnetcdf.a
#NLIBS   =	-I/sw/include  /usr/local/lib/libnetcdf.a
#NLIBS   =	-I/sw/include  -L/sw/lib -lnetcdf -ff2c
#NLIBS   =	-I/usr/local/include  -L/usr/local/lib -lnetcdf
NLIBS   =	-I/opt/local/include  -L/opt/local/lib -lnetcdff -lnetcdf
#NLIBS   =	/usr/local/lib/libnetcdf.a
#NLIBS   =	-I/Users/shamis/opt/netcdf/4.4.4/include  -L/Users/shamis/opt/netcdf/4.4.4/lib -lnetcdff -lnetcdf

#

all:	map5.x pack.x packz.x pack1.x unpack.x make1.x cont4a.x fix.x\
	velo6.x velo5.x pltold.x pltall.x velo7.x convert.x latlong.x \
        massb.x combine.x cont3.x vlatlong.x latlong1.x campcen.x \
        smooth.x packa.x contour.x mkplotc.x plotps.x formstf.x \
        hplotgifs.x hplot4.x contmake.x image.x strip.x diff.x vplot.x \
        gmtmake-th.x gmtmake.x percent.x percenta.x maketopo.x cont4b.x \
        gmtmake_new.x forman.x \
        velo4.x velo4a.x tempstf.x gmtall.x doubres.x antconv.x \
	acon.x fudgel.x velostf.x ttlist.x hard.x \
	eisstf.x temppoint.x profile.x \
	velorgb.x dplot4.x gmtdep.x point.x contrgb.x \
	select.x timer.x locater.x smoothl.x erosion.x rose.x diag.x \
        velo5rgb.x tslice.x extract.x \
	kleman.x velocont.x plotpoint.x netcdf.x Bunpack.x \
	veloExtract.x xyz2rll.x

some:	pack.x packz.x pack1.x unpack.x make1.x cont4a.x fix.x\
	velo6.x velo5.x pltold.x pltall.x velo7.x convert.x latlong.x \
        massb.x combine.x cont3.x vlatlong.x latlong1.x campcen.x \
        smooth.x packa.x contour.x mkplotc.x plotps.x formstf.x \
        hplot4.x contmake.x image.x strip.x diff.x vplot.x \
        gmtmake.x percent.x percenta.x maketopo.x cont4b.x forman.x \
        velo4.x velo4a.x tempstf.x gmtall.x doubres.x antconv.x \
	acon.x fudgel.x velostf.x ttlist.x hard.x temptest.x \
	ablation.x eisstf.x temppoint.x profile.x \
	accumtest.x velorgb.x dplot4.x gmtdep.x point.x contrgb.x \
	select.x timer.x locater.x smoothl.x erosion.x rose.x diag.x \
        velo5rgb.x wmtest.x tslice.x extract.x kleman.x \
	velocont.x plotpoint.x netcdf.x

tar.e:
	cp makefile makefile.mac
	tar -cvf source.tar mapmain5.f $(DSRCS5) threedc.c map5.e makefile.mac

mesh.x: mesh.f 
	$(CC) $(GFLAGS) -o  mesh.x mesh.f

genmake.x: genmake.f parameter.h $(PROJS)
	$(CC) $(GFLAGS) -o  genmake.x genmake.f $(GLIBS) $(LIBS)

projections.o: projections.f
	$(CC) -c $(CFLAGS) projections.f

whole5.e:
	cat mapmain5.f $(DSRCS5) > whole5.f

veplate.o: veplate.f
	$(CC) -c $(CFLAGS) veplate.f

splate.o: splate.f
	$(CC) -c $(CFLAGS) splate.f

threedc.x: threedc.c 
	gcc $(CFLAGS) -o  threedc.x threedc.c $(OGLLIBS)

xyz2rll.x: xyz2rll.f 
	$(CC) $(CFLAGS) -o  xyz2rll.x xyz2rll.f 

vetest.x: vetest.f veplate.o $(DOBJS7)
	$(CC) $(CFLAGS) -o  vetest.x vetest.f $(DOBJS7) $(GLIBS) $(LIBS)

veopt.x: vetest.f veplate.o $(DOBJS7)
	$(CC) -c $(OFLAGS) veplate.f
	$(CC) $(OFLAGS) -o  veopt.x vetest.f $(DOBJS7) $(GLIBS) $(LIBS)
	cp veopt.x vetest.x

whole5.x:whole5.f
	$(CC) $(CFLAGS) -o  whole5.x whole5.f $(GLIBS) $(LIBS)

whole5opt.x: whole5.f
	$(CC) $(OFLAGS) -o  whole5.x whole5.f $(GLIBS) $(LIBS)

incore.x: incore.f
	$(CC) $(OFLAGS) -o  incore.x incore.f $(GLIBS) $(LIBS)

whole5par.x: whole5.f
	$(CC) $(OFLAGS) -apo -o  whole5.x whole5.f $(GLIBS) $(LIBS)

fast.x: 
	cat mapmain5.f $(DSRCS5) > whole5.f
	$(CC) $(FFLAGS) -o  map5.x whole5.f $(GLIBS) $(LIBS)

prettyfast.x: whole5.f
	cat mapmain5.f $(DSRCS5) > whole5.f
	$(CC) $(OFLAGS) -o  map5.x whole5.f $(GLIBS) $(LIBS)

whole6.e:
	cat mapmain6.f $(DSRCS6) > whole6.f

vplate.o: parameter.h vplate.f
	$(CC) -c $(CFLAGS) vplate.f

oplate.o: parameter.h oplate.f
	$(CC) -c $(CFLAGS) oplate.f

jcg.o: jcg.f
	$(CC) -c $(CFLAGS) jcg.f

lookup.o: lookup.f
	$(CC) -c $(CFLAGS) lookup.f

wthikn.o: wthikn.f
	$(CC) -c $(CFLAGS) wthikn.f

wthik.o: wthik.f
	$(CC) -c $(CFLAGS) wthik.f

wthiksimp.o: wthiksimp.f
	$(CC) -c $(CFLAGS) wthiksimp.f

wthikj.o: wthikj.f
	$(CC) -c $(CFLAGS) wthikj.f

accum.o: parameter.h accum.f
	$(CC) -c $(CFLAGS) accum.f

$(THREED).o: parameter.h $(THREED).f
	$(CC) -c $(CFLAGS) $(THREED).f

modify.o: modify.f
	$(CC) -c $(CFLAGS) modify.f

slopef.o: slopef.f
	$(CC) -c $(CFLAGS) slopef.f

temper.o: parameter.h temper.f
	$(CC) -c $(CFLAGS) temper.f

gauseid.o: parameter.h gauseid.f
	$(CC) -c $(CFLAGS) gauseid.f

gnumbr.o: gnumbr.f
	$(CC) -c $(CFLAGS) gnumbr.f

asymsl.o: parameter.h asymsl.f
	$(CC) -c $(CFLAGS) asymsl.f

gausinit.o: gausinit.f
	$(CC) -c $(CFLAGS) gausinit.f

nconst.o: nconst.f
	$(CC) -c $(CFLAGS) nconst.f

readn.o: readn.f
	$(CC) -c $(CFLAGS) readn.f

shape.o: shape.f
	$(CC) -c $(CFLAGS) shape.f

sparse.o: sparse.f
	$(CC) -c $(CFLAGS) sparse.f

adjust.o: parameter.h adjust.f
	$(CC) -c -O0 adjust.f

elprop.o: elprop.f
	$(CC) -c $(CFLAGS) elprop.f

nodesl.o: parameter.h nodesl.f
	$(CC) -c $(CFLAGS) nodesl.f

scale.o: scale.f
	$(CC) -c $(CFLAGS) scale.f

volume.o: volume.f
	$(CC) -c $(CFLAGS) volume.f

afunct.o: afunct.f
	$(CC) -c $(CFLAGS) afunct.f

forma.o: forma.f
	$(CC) -c $(CFLAGS) forma.f

formc.o: formc.f
	$(CC) -c $(CFLAGS) formc.f

plotsol.o: parameter.h plotsol.f
	$(CC) -c $(CFLAGS) plotsol.f

setrig.o: setrig.f
	$(CC) -c $(CFLAGS) setrig.f

zoom.o: zoom.f
	$(CC) -c $(CFLAGS) zoom.f

asdf.x: asdf.f
	$(CC) $(CFLAGS) -o  asdf.x asdf.f $(GLIBS) $(LIBS)

test.x:setrig.o test.f
	$(CC) $(CFLAGS) -o  test.x test.f $(GLIBS) $(LIBS)

fudgel.x: fudgel.f parameter.h
	$(CC) $(CFLAGS) -o  fudgel.x fudgel.f $(GLIBS) $(LIBS)

test3d.x: test3d.f
	$(CC) $(CFLAGS) -o  test3d.x test3d.f $(GLIBS) $(LIBS)

velostf.x: velostf.f
	$(CC) $(CFLAGS) -o  velostf.x velostf.f $(GLIBS) $(LIBS)

acon.x: acon.f parameter.h
	$(CC) $(CFLAGS) -o  acon.x acon.f $(GLIBS) $(LIBS)

contour.x: contour.f
	$(CC) $(CFLAGS) -o  contour.x contour.f $(GLIBS) $(LIBS)

kleman.x: kleman.f parameter.h
	$(CC) $(CFLAGS) -o  kleman.x kleman.f $(GLIBS) $(LIBS)

gmtmake_new.x: gmtmake_new.f parameter.h
	$(CC) -O1 -o  gmtmake_new.x gmtmake_new.f $(GLIBS) $(LIBS)

gmtmake-th.x: gmtmake-th.f parameter.h
	$(CC) -O1 -o  gmtmake-th.x gmtmake-th.f $(GLIBS) $(LIBS)

gmtmake.x: gmtmake.f parameter.h
	$(CC) -O1 -o  gmtmake.x gmtmake.f $(GLIBS) $(LIBS)

gmtall.x: gmtall.f
	$(CC) $(CFLAGS) -o  gmtall.x gmtall.f $(LIBS)

doubres.x: doubres.f parameter.h
	$(CC) $(CFLAGS) -o  doubres.x doubres.f $(LIBS)

percent.x: percent.f
	$(CC) $(CFLAGS) -o  percent.x percent.f $(GLIBS) $(LIBS)

percenta.x: percenta.f
	$(CC) $(CFLAGS) -o  percenta.x percenta.f $(GLIBS) $(LIBS)

erosion.x: erosion.f
	$(CC) $(CFLAGS) -o  erosion.x erosion.f $(GLIBS) $(LIBS)

rose.x: rose.f
	$(CC) $(CFLAGS) -o  rose.x rose.f $(GLIBS) $(LIBS) $(SLIBS)

diag.x: diag.f
	$(CC) $(CFLAGS) -o  diag.x diag.f $(GLIBS) $(LIBS)

maketopo.x: maketopo.f
	$(CC) $(CFLAGS) -o  maketopo.x maketopo.f $(GLIBS) $(LIBS)

mkplotc.x: mkplotc.f
	$(CC) $(CFLAGS) -o  mkplotc.x mkplotc.f $(GLIBS) $(LIBS)

strip.x: strip.f
	$(CC) $(CFLAGS) -o  strip.x strip.f $(GLIBS) $(LIBS)

diff.x: diff.f parameter.h
	$(CC) $(CFLAGS) -o  diff.x diff.f $(GLIBS) $(LIBS)

vplot.x: vplot.f
	$(CC) $(CFLAGS) -o  vplot.x vplot.f $(GLIBS) $(LIBS)

plotps.x: plotps.f
	$(CC) $(CFLAGS) -o  plotps.x plotps.f $(GLIBS) $(LIBS)

campcen.x:setrig.o campcen.f
	$(CC) $(CFLAGS) -o  campcen.x campcen.f $(GLIBS) $(LIBS)

smooth.x:smooth.f
	$(CC) $(CFLAGS) -o  smooth.x smooth.f $(GLIBS) $(LIBS)

smoothl.x:smoothl.f
	$(CC) $(CFLAGS) -o  smoothl.x smoothl.f $(GLIBS) $(LIBS)

latlong.x:setrig.o latlong.f
	$(CC) $(OFLAGS) -o  latlong.x latlong.f setrig.o $(GLIBS) $(LIBS)

latlong1.x:setrig.o latlong1.f
	$(CC) $(OFLAGS) -o  latlong1.x latlong1.f setrig.o $(GLIBS) $(LIBS)

vlatlong.x: vlatlong.f
	$(CC) $(OFLAGS) -o  vlatlong.x vlatlong.f setrig.o $(GLIBS) $(LIBS)

mapq.x:	$(DOBJS) mapmain.f
	$(CC) $(CFLAGS) -o  mapq.x mapmain.f $(DOBJS) $(GLIBS) $(LIBS)

map5.x:	$(DOBJS5) mapmain5.f parameter.h
	$(CC) $(CFLAGS) -o  map5.x mapmain5.f $(DOBJS5) $(GLIBS) $(LIBS)

accumtest.x:	$(DOBJS5) accumtest.f parameter.h
	$(CC) $(OFLAGS) -o  accumtest.x accumtest.f $(DOBJS5) $(GLIBS) $(LIBS)

wmtest.x:	$(DOBJS5) wmtest.f parameter.h
	$(CC) $(OFLAGS) -o  wmtest.x wmtest.f $(DOBJS5) $(GLIBS) $(LIBS)

profile.x: profile.f
	$(CC) $(CFLAGS) -o  profile.x profile.f $(GLIBS) $(LIBS) $(SLIBS)

temppoint.x: temppoint.f
	$(CC) $(CFLAGS) -o  temppoint.x temppoint.f $(GLIBS) $(LIBS) $(SLIBS)

tempslice.x:	 tempslice.f
	$(CC) $(CFLAGS) -o  tempslice.x tempslice.f $(GLIBS) $(LIBS) $(SLIBS)

velorgb.x:	 velorgb.f parameter.h
	$(CC) $(CFLAGS) -o  velorgb.x velorgb.f $(GLIBS) $(LIBS) $(SLIBS)

velocont.x:	 velocont.f parameter.h
	$(CC) $(CFLAGS) -o  velocont.x velocont.f $(GLIBS) $(LIBS) $(SLIBS)

tslice.x:	 tslice.f parameter.h
	$(CC) $(CFLAGS) -o  tslice.x tslice.f $(GLIBS) $(LIBS) $(SLIBS)

extract.x:	 extract.f
	$(CC) $(CFLAGS) -o  extract.x extract.f $(GLIBS) $(LIBS) $(SLIBS)

contrgb.x:	 contrgb.f parameter.h
	$(CC) $(CFLAGS) -o  contrgb.x contrgb.f $(GLIBS) $(LIBS) $(SLIBS)

select.x:	 select.f
	$(CC) $(CFLAGS) -o  select.x select.f $(GLIBS) $(LIBS) $(SLIBS)

locater.x:	 locater.f
	$(CC) $(CFLAGS) -o  locater.x locater.f $(GLIBS) $(LIBS) $(SLIBS)

timer.x:	 timer.f
	$(CC) $(CFLAGS) -o  timer.x timer.f $(GLIBS) $(LIBS) $(SLIBS)

temptest.x:	$(DOBJS5) temptest.f
	$(CC) $(OFLAGS) -o  temptest.x temptest.f $(DOBJS5) $(GLIBS) $(LIBS)

map6.x:	$(DOBJS5) mapmain6.f
	$(CC) $(CFLAGS) -o  map6.x mapmain6.f $(DOBJS6) $(GLIBS) $(LIBS)

massb.x:massb.f
	$(CC) $(CFLAGS) -o  massb.x massb.f $(GLIBS) $(LIBS)
 
eisstf.x:eisstf.f
	$(CC) $(CFLAGS) -o  eisstf.x eisstf.f $(GLIBS) $(LIBS)
 
ablation.x:$(DOBJS5) ablation.f
	$(CC) $(OFLAGS) -o  ablation.x ablation.f $(DOBJS5) $(GLIBS) $(LIBS)

combine.x:combine.f parameter.h
	$(CC) $(CFLAGS) -o  combine.x combine.f $(GLIBS) $(LIBS)

hp7.x:hp7.f
	$(CC) $(CFLAGS) -o  hp7.x hp7.f $(GLIBS) $(LIBS)

convert.x:convert.f parameter.h
	$(CC) $(CFLAGS) -o  convert.x convert.f $(GLIBS) $(LIBS)

antconv.x:antconv.f parameter.h
	$(CC) $(CFLAGS) -o  antconv.x antconv.f $(GLIBS) $(LIBS)

fix.x:fix.f parameter.h
	$(CC) $(CFLAGS) -o  fix.x fix.f $(GLIBS) $(LIBS)

pack.x:	pack.f
	$(CC) $(CFLAGS) -o  pack.x pack.f $(GLIBS) $(LIBS)

packa.x:	packa.f
	$(CC) $(CFLAGS) -o  packa.x packa.f $(GLIBS) $(LIBS)

packz.x:	packz.f
	$(CC) $(CFLAGS) -o  packz.x packz.f $(GLIBS) $(LIBS)

pack1.x:	pack1.f parameter.h
	$(CC) $(CFLAGS) -o  pack1.x pack1.f $(GLIBS) $(LIBS)

formstf.x:	formstf.f parameter.h
	$(CC) $(CFLAGS) -o  formstf.x formstf.f $(GLIBS) $(LIBS)

tempstf.x:	tempstf.f
	$(CC) $(CFLAGS) -o  tempstf.x tempstf.f $(GLIBS) $(LIBS)

Bunpack.x:	Bunpack.f parameter.h
	$(CC) $(CFLAGS) -o  Bunpack.x Bunpack.f $(GLIBS) $(LIBS)

unpack.x:	unpack.f parameter.h
	$(CC) $(CFLAGS) -o  unpack.x unpack.f $(GLIBS) $(LIBS)

cont3.x:cont3.f parameter.h
	$(CC) $(CFLAGS) -o  cont3.x cont3.f $(GLIBS) $(LIBS)

make1.x:make1.f
	$(CC) $(CFLAGS) -o  make1.x make1.f $(LIBS)

cont4a.x:cont4a.f parameter.h
	$(CC) $(CFLAGS) -o  cont4a.x cont4a.f $(GLIBS) $(LIBS)

cont4b.x:cont4b.f parameter.h
	$(CC) $(CFLAGS) -o  cont4b.x cont4b.f $(GLIBS) $(LIBS)

forman.x:forman.f parameter.h
	$(CC) $(CFLAGS) -o  forman.x forman.f $(GLIBS) $(LIBS)

pltall.x: pltall.f
	$(CC) $(CFLAGS) -o  pltall.x pltall.f $(GLIBS) $(LIBS)

pltold.x: pltold.f
	$(CC) $(CFLAGS) -o  pltold.x pltold.f $(GLIBS) $(LIBS)

packallnc.x:packallnc.f
	$(G77) $(GFLAGS) -o  packallnc.x packallnc.f $(NLIBS)

netcdf.x:netcdf.f
	$(G77) -m64 $(OFLAGS) -o netcdf.x netcdf.f $(NLIBS)

velo4.x:velo4.f
	$(CC) $(CFLAGS) -o  velo4.x velo4.f $(GLIBS) $(LIBS)

velo4a.x:velo4a.f
	$(CC) $(CFLAGS) -o  velo4a.x velo4a.f $(GLIBS) $(LIBS)

veloExtract.x:veloExtract.f
	$(CC) $(CFLAGS) -o  veloExtract.x veloExtract.f $(GLIBS) $(LIBS)

velo5.x:velo5.f
	$(CC) $(CFLAGS) -o  velo5.x velo5.f $(GLIBS) $(LIBS)

velo5rgb.x:velo5rgb.f
	$(CC) $(CFLAGS) -o  velo5rgb.x velo5rgb.f $(GLIBS) $(LIBS) $(SLIBS)

velo6.x:velo6.f
	$(CC) $(CFLAGS) -o  velo6.x velo6.f $(GLIBS) $(LIBS)

velo7.x:velo7.f
	$(CC) $(CFLAGS) -o  velo7.x velo7.f $(GLIBS) $(LIBS)

hplot2.x:hplot2.f
	$(CC) $(CFLAGS) -o  hplot2.x hplot2.f $(GLIBS) $(LIBS)

hplot1.x:hplot1.f
	$(CC) $(CFLAGS) -o  hplot1.x hplot1.f $(GLIBS) $(LIBS)

hplot4.x:hplot4.f
	$(CC) $(CFLAGS) -o  hplot4.x hplot4.f $(GLIBS) $(LIBS)

hplotgifs.x:hplotgifs.f parameter.h
	$(CC) $(CFLAGS) -o  hplotgifs.x hplotgifs.f $(GLIBS) $(LIBS)

dplot4.x:dplot4.f parameter.h
	$(CC) $(CFLAGS) -o  dplot4.x dplot4.f $(GLIBS) $(LIBS)

gmtdep.x:gmtdep.f parameter.h
	$(CC) $(CFLAGS) -o  gmtdep.x gmtdep.f $(GLIBS) $(LIBS)

plotpoint.x: plotpoint.f
	$(CC) $(CFLAGS) -o  plotpoint.x plotpoint.f $(GLIBS) $(LIBS) $(SLIBS)

point.x:point.f
	$(CC) $(CFLAGS) -o  point.x point.f $(GLIBS) $(LIBS)

sep.x:	sep.f
	$(CC) $(CFLAGS) -o  sep.x sep.f $(GLIBS) $(LIBS)

ttlist.x:ttlist.f
	$(CC) $(CFLAGS) -o  ttlist.x ttlist.f $(GLIBS) $(LIBS)

hard.x:hard.f
	$(CC) $(CFLAGS) -o  hard.x hard.f $(GLIBS) $(LIBS)

zones.x:zones.f
	$(CC) $(CFLAGS) -o  zones.x zones.f $(GLIBS) $(LIBS)

matrix.x:matrix.f
	$(CC) $(CFLAGS) -o  matrix.x matrix.f $(LIBS)

contmake.x:contmake.f parameter.h
	$(CC) $(CFLAGS) -o  contmake.x contmake.f $(LIBS)

image.x:image.f
	$(CC) $(CFLAGS) -o  image.x image.f $(GLIBS) $(LIBS)

clean:
	rm -f $(DOBJS5) *.x 

opt5.x:	$(DOBJS5) mapmain5.f
	$(CC) $(OFLAGS) -o  map5.x mapmain5.f $(DOBJS5) $(GLIBS) $(LIBS)

production:
	rm -f $(DOBJS5) 
	$(CC)  -c $(OFLAGS) accum.f
	$(CC)  -c $(OFLAGS) gausinit.f
	$(CC)  -c $(OFLAGS) nconst.f
	$(CC)  -c $(OFLAGS) readn.f
	$(CC)  -c $(OFLAGS) shape.f
	$(CC)  -c $(OFLAGS) elprop.f
	$(CC)  -c $(OFLAGS) nodesl.f
	$(CC)  -c $(OFLAGS) scale.f
	$(CC)  -c $(OFLAGS) volume.f
	$(CC)  -c $(OFLAGS) afunct.f
	$(CC)  -c $(OFLAGS) plotsol.f
	$(CC)  -c $(OFLAGS) setrig.f
	$(CC)  -c $(OFLAGS) zoom.f
	$(CC)  -c $(OFLAGS) formc.f
	$(CC)  -c $(OFLAGS) gauseid.f
	$(CC)  -c $(OFLAGS) modify.f
	$(CC)  -c $(OFLAGS) gnumbr.f
	$(CC)  -c $(OFLAGS) slopef.f
	$(CC)  -c $(OFLAGS) temper.f
	$(CC)  -c O0 adjust.f
	$(CC)  -c $(OFLAGS) $(THREED).f
	$(CC)  -c $(OFLAGS) wthikn.f
	$(CC)  -c $(OFLAGS) wthik.f
	$(CC)  -c $(OFLAGS) wthiksimp.f
	$(CC)  -c $(OFLAGS) wthikj.f
	$(CC)  -c $(OFLAGS) plate.f
	$(CC)  -c $(OFLAGS) vplate.f
	$(CC)  -c $(OFLAGS) jcg.f
	$(CC)  -c $(OFLAGS) lookup.f
	$(CC) $(OFLAGS)  -o map5.x mapmain5.f $(DOBJS5) $(GLIBS) $(LIBS)

update:
	rm -f $(DOBJS5) 
	$(CC) -c $(OFLAGS) accum.f
	$(CC) -c $(OFLAGS) gausinit.f
	$(CC) -c $(OFLAGS) nconst.f
	$(CC) -c $(OFLAGS) readn.f
	$(CC) -c $(OFLAGS) shape.f
	$(CC) -c $(OFLAGS) elprop.f
	$(CC) -c $(OFLAGS) nodesl.f
	$(CC) -c $(OFLAGS) scale.f
	$(CC) -c $(OFLAGS) volume.f
	$(CC) -c $(OFLAGS) afunct.f
	$(CC) -c $(OFLAGS) plotsol.f
	$(CC) -c $(OFLAGS) setrig.f
	$(CC) -c $(OFLAGS) zoom.f
	$(CC) -c $(OFLAGS) formc.f
	$(CC) -c $(OFLAGS) gauseid.f
	$(CC) -c $(OFLAGS) modify.f
	$(CC) -c $(OFLAGS) gnumbr.f
	$(CC) -c $(OFLAGS) slopef.f
	$(CC) -c $(OFLAGS) temper.f
	$(CC) -c -O0 adjust.f
	$(CC) -c $(OFLAGS) $(THREED).f
	$(CC) -c $(OFLAGS) wthikn.f
	$(CC) -c $(OFLAGS) wthik.f
	$(CC) -c $(OFLAGS) wthiksimp.f
	$(CC) -c $(OFLAGS) wthikj.f
	$(CC) -c $(OFLAGS) oplate.f
	$(CC) -c $(OFLAGS) veplate.f
	$(CC) -c $(OFLAGS) splate.f
	$(CC) -c $(OFLAGS) jcg.f
	$(CC) -c $(OFLAGS) lookup.f
	$(CC) $(OFLAGS) -o  map5.x mapmain5.f $(DOBJS5) $(GLIBS) $(LIBS)

debug:
	rm -f $(DOBJS5) 
	$(CC) -c $(DFLAGS) accum.f
	$(CC) -c $(DFLAGS) gausinit.f
	$(CC) -c $(DFLAGS) nconst.f
	$(CC) -c $(DFLAGS) readn.f
	$(CC) -c $(DFLAGS) shape.f
	$(CC) -c $(dGLAGS) elprop.f
	$(CC) -c $(DFLAGS) nodesl.f
	$(CC) -c $(DFLAGS) scale.f
	$(CC) -c $(DFLAGS) volume.f
	$(CC) -c $(DFLAGS) afunct.f
	$(CC) -c $(DFLAGS) plotsol.f
	$(CC) -c $(DFLAGS) setrig.f
	$(CC) -c $(DFLAGS) zoom.f
	$(CC) -c $(DFLAGS) formc.f
	$(CC) -c $(DFLAGS) gauseid.f
	$(CC) -c $(DFLAGS) modify.f
	$(CC) -c $(DFLAGS) gnumbr.f
	$(CC) -c $(DFLAGS) slopef.f
	$(CC) -c $(DFLAGS) temper.f
	$(CC) -c -O0 adjust.f
	$(CC) -c $(DFLAGS) -o $(THREED).o $(THREED).f
	$(CC) -c $(DFLAGS) wthikn.f
	$(CC) -c $(DFLAGS) wthik.f
	$(CC) -c $(DFLAGS) wthiksimp.f
	$(CC) -c $(DFLAGS) wthikj.f
	$(CC) -c $(DFLAGS) oplate.f
	$(CC) -c $(DFLAGS) veplate.f
	$(CC) -c $(DFLAGS) splate.f
	$(CC) -c $(DFLAGS) jcg.f
	$(CC) -c $(DFLAGS) lookup.f
	$(CC) $(DFLAGS) -o  map5.x mapmain5.f $(DOBJS5) $(GLIBS) $(LIBS)

ftn:
	rm -f *.ftn
	ftnchek accum.f > accum.ftn
	ftnchek gausinit.f > gausinit.ftn
	ftnchek nconst.f > nconst.ftn
	ftnchek readn.f > readn.ftn
	ftnchek shape.f > shape.ftn
	ftnchek elprop.f > elprop.ftn
	ftnchek nodesl.f > nodesl.ftn
	ftnchek scale.f > scale.ftn
	ftnchek volume.f > volume.ftn
	ftnchek afunct.f > afunct.ftn
	ftnchek plotsol.f > plotsol.ftn
	ftnchek setrig.f > setrig.ftn
	ftnchek zoom.f > zoom.ftn
	ftnchek formc.f > formc.ftn
	ftnchek gauseid.f > gauseid.ftn
	ftnchek modify.f > modify.ftn
	ftnchek gnumbr.f > gnumbr.ftn
	ftnchek slopef.f > slopef.ftn
	ftnchek temper.f > temper.ftn
	ftnchek adjust.f > adjust.ftn
	ftnchek $(THREED).f > $(THREED).ftn
	ftnchek wthikn.f > wthikn.ftn
	ftnchek wthik.f > wthik.ftn
	ftnchek wthiksimp.f > wthiksimp.ftn
	ftnchek wthikj.f > wthikj.ftn
	ftnchek oplate.f > oplate.ftn
	ftnchek veplate.f > veplate.ftn
	ftnchek splate.f > splate.ftn
	ftnchek jcg.f > jcg.ftn
	ftnchek lookup.f > lookup.ftn
	ftnchek mapmain5.f > mapmain5.ftn

single:
	$(CC) $(OFLAGS) -o map5.x mapmain5.f $(DOBJS5) $(GLIBS) $(LIBS)

sparallel:
	$(CC) $(OFLAGS)  -o map5.x mapmain5.f $(DOBJS5) $(GLIBS) $(LIBS)

tidy:
	rm -f $(DOBJS5) *.o *.bak *.ckp fort.*
