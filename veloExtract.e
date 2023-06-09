rm fort.*
ln -s input$1.head fort.30
ln -s input$1.grid fort.31
ln -s input$1.diff fort.32
ln -s out$2.time fort.33
ln -s bc$2.time fort.34
ln -s QQ-vMax$2.data fort.11
ln -s QQ-vAvg$2.data fort.12
ln -s QQ-vNum$2.data fort.13
ln -s QQ-hAvg$2.data fort.14
./veloExtract.x
