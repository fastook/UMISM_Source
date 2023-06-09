rm rose1.d rose2.d
diag.x
pstext -P -K  -R0/11/0/8.5 -Jx1 -X2i -Y3i << END > rose1.ps
0  6.5 20 0 1 1  Scandinavia Flow Patterns (5-deg) ($1)                
END
pstext -P -K -O -R0/11/0/8.5 -Jx1   << END >> rose1.ps
0.2  6.3 10 0 1 1 (red: flow duration (yrs) sliding only)
END
pstext rose.locate -V -U -P -K -O -R0/11/0/8.5 -Jx1  >> rose1.ps
pstext rose.time -V -U -P -K -O -R0/11/0/8.5 -Jx1  >> rose1.ps
psrose rose1.d     -V -P  -O -: -A5 -S2.5i -G0/0/100 \
       -R0/4000/0/360 -X-0.5i -Y0.25i -B500g500/30g5 >> rose1.ps
ghostview rose1.ps
cp rose1.ps rosec1-$1.ps

