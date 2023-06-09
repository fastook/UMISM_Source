rm fort.*
ln -s volt$1.data fort.1
ln -s plot$1.data fort.11
ln -s tplot$1.data fort.21
ln -s vtplot$1.data fort.22
ln -s vfplot$1.data fort.23
ln -s aplot$1.data fort.24
ln -s hplot$1.data fort.25
ln -s dplot$1.data fort.26
ln -s poplot$1.data fort.27
ln -s tbplot$1.data fort.28
ln -s twplot$1.data fort.29
ln -s tpplot$1.data fort.30
ln -s wmplot$1.data fort.31
ln -s vsplot$1.data fort.99
./pltall.x
