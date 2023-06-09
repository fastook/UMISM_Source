rm fort.*
ln -s o2$3.data fort.2
ln -s o3$3.data fort.3
ln -s o4$3.data fort.4
ln -s o7$3.data fort.7
ln -s input$1.head fort.30
ln -s input$1.grid fort.31
ln -s input$1.diff fort.32
ln -s dep$2.time fort.33
ln -s linem.data fort.11
ln -s game.data fort.23
ln -s dome.data fort.24
ln -s volume.data fort.25
ln -s prof.data fort.26
dplot4.x
rm junk.data
cat o4$3.data o2$3.data o3$3.data > junk.data
rm o4$3.data
mv junk.data o4$3.data
