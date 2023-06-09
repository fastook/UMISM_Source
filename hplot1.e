rm fort.*
ln -s o2.data fort.2
ln -s o3.data fort.3
ln -s o4.data fort.4
ln -s o7.data fort.7
ln -s o9.data fort.9
ln -s input$1.head fort.30
ln -s input$1.grid fort.31
ln -s input$1.diff fort.32
ln -s out$2.time fort.33
ln -s bc$2.time fort.34
ln -s linem.data fort.11
ln -s game.data fort.23
ln -s dome.data fort.24
ln -s volume.data fort.25
ln -s prof.data fort.26
hplot1.x << END
$3
END
#rm junk.data
#cat o4.data o2.data o3.data > junk.data
#cat o2.data o3.data > junk.data
#rm o4.data
cp o2.data junkht$3.data
cp o3.data junkbed.data
cp o7.data junkasurf.data
