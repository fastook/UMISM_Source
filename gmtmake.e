rm fort.*
ln -s out3$1.data fort.1
ln -s tmp$1.h2o fort.2
ln -s gmt$2.rll fort.9
ln -s white$2.xyz fort.19
ln -s gmt$2.XYZ fort.10
ln -s gmt$2.e fort.11
ln -s img$2.e fort.12
ln -s GMT$2.e fort.14
ln -s MOS$2.e fort.16
ln -s topo.cpt fort.13
./gmtmake.x
chmod +x gmt$2.e
chmod +x img$2.e
chmod +x GMT$2.e
#chmod +x MOS$2.e

