rm fort.*
ln -s outc$1.data fort.7
ln -s out3$1.data fort.1
ln -s input$1.head fort.30
ln -s input$1.grid fort.31
ln -s input$1.diff fort.32
ln -s out$2.time fort.33
ln -s bc$2.time fort.34
ln -s temp$2.time fort.36
ln -s load$2.data fort.21
ln -s depb$2.data fort.22
ln -s htice$2.data fort.23
ln -s rate$2.data fort.24
ln -s outline.data fort.12
ln -s outt.data fort.9
point.x
cat load$2.data depb$2.data rate$2.data > plot$2.data
#plotbps.e plot$2
#plotbps.e depb$2
#plotbps.e rate$2
plotbps.e htice$2
