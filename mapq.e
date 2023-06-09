rm fort.*
ln -s input$1 fort.9
ln -s input$1.head fort.30
ln -s input$1.grid fort.31
ln -s input$1.diff fort.32
ln -s input$1.time fort.33
ln -s out$2.time fort.34
ln -s outline.data fort.10
ln -s rad$1.data fort.4
ln -s out1$1.data fort.7
ln -s out2$1.data fort.12
ln -s out3$1.data fort.11
ln -s out4$1.data fort.21
ln -s linem.data fort.20
ln -s output$1.data fort.26
ln -s mat$2.data fort.13
ln -s vol$2.data fort.18
ln -s volt$2.data fort.17
ln -s console.$2 fort.99
time mapq.x
