# UMISM_Source
All can be made for MacOS with makefile command, although you will also need the grstrt files availble also here on Github. 
The main program executable is map5.x, run with the script map5.e. 
All input.data files must be preprocessed with unpack.e, which runs the executable unpack.x.
Results can be viewed with packz.e (runs executable packz.x) that spawns a generic gmt script to draw the appropriate figure desired.
Those gmt scripts themselves will need to run gmtmake_new.e (runs execuatable gmtmake_new.x) to extract appropriate results.
