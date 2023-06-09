rm fort.*
ln -s out3$3.data fort.1
ln -s tmp$3.h2o fort.2
ln -s input$1.head fort.30
ln -s input$1.grid fort.31
ln -s input$1.diff fort.32
ln -s input$1.coord fort.37
ln -s out$2.time fort.33
ln -s bc$2.time fort.34
ln -s temp$2.time fort.36
./pack1.x
