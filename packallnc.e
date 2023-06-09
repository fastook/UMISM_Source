rm input.def input.head input.grid input.diff input.time
rm out.time bc.time temp.time dep.time
ln -s input$1.def input.def
ln -s input$1.head input.head
ln -s input$1.grid input.grid
ln -s input$1.diff input.diff
ln -s input$1.time input.time
ln -s out$2.time out.time
ln -s bc$2.time bc.time
ln -s temp$2.time temp.time
ln -s dep$2.time dep.time
./packallnc.x
mv packall.nc packall-$3.nc

