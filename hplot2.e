rm fort.*
ln -s o2$1.data fort.2
ln -s o3$1.data fort.3
ln -s o4$1.data fort.4
ln -s o7$1.data fort.7
ln -s out3$1.data fort.1
ln -s line$2.data fort.11
ln -s game$2.data fort.23
ln -s profs.data fort.32
ln -s profb.data fort.33
hplot2.x
rm junk.data
cat o2$1.data o3$1.data > junk$1.data

