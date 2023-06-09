# UMISM_Source
All can be made for MacOS with makefile command, although you will also need the grstrt files availble in grstrt.zip. 
The main program executable is map5.x, run with the script map5.e. 
All input.data files must be preprocessed with unpack.e, which runs the executable unpack.x.
Results can be viewed with packz.e (runs executable packz.x) that spawns a generic gmt script to draw the appropriate figure desired.
Those gmt scripts (in generic-gmts.zip) themselves will run gmtmake_new.e (runs execuatable gmtmake_new.x) to extract appropriate results.
There are many unused programs for various special purposes that are not necessary to run UMISM

Running UMISM with a data file inputE.data:

unpack.e E E

map5.e E E

...  various interactive commands to control UMISM

packz.e E E WW

...  label numbers, output interval, start and end times to extract results displayed with generic.gmt

