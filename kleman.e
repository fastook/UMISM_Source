pack1.e $1 $1 << END
$2
END
rm fort.*
ln -s out3.data fort.1
ln -s kleman$2.xyz fort.9
./kleman.x

