rm fort.*
ln -s $1.plot fort.1
ln -s outline.data fort.2
ln -s $1.namelist fort.11
mkplotc.x
