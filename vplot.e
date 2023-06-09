rm fort.*
ln -s outt$1.data fort.1
ln -s outt$2.data fort.2
ln -s outt$1.namelist fort.11
ln -s outt$2.namelist fort.12
vplot.x
