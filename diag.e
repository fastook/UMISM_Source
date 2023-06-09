rm rose1.d rose2.d
diag.x
pstext -P -K  -R0/11/0/8.5 -Jx1 -X2i -Y3i << END > rose1.ps
0  6.5 20 0 1 1  Flow Patterns (5-deg) ($1)                
END
pstext -P -K -O -R0/11/0/8.5 -Jx1   << END >> rose1.ps
0.2  6.3 10 0 1 1 (red: flow frequency)
END
pstext rose.locate -V -U -P -K -O -R0/11/0/8.5 -Jx1  >> rose1.ps
pstext rose.time -V -U -P -K -O -R0/11/0/8.5 -Jx1  >> rose1.ps
psrose rose1.d     -V -P  -O -: -A5r -S2.5in -N -G100/0/0 \
       -R0/1/0/360 -X-0.5i -Y0.25i -B1g0.1/30g10 >> rose1.ps
pstext -P -K  -R0/11/0/8.5 -Jx1 -X2i -Y3i << END > rose2.ps
0  6.5 20 0 1 1  Flow Patterns (5-deg) ($1)
END
pstext -P -K -O -R0/11/0/8.5 -Jx1   << END >> rose2.ps
0.2  6.3 10 0 1 1 (green: flow distance)
END
pstext rose.locate -V -U -P -K -O -R0/11/0/8.5 -Jx1  >> rose2.ps
pstext rose.time -V -U -P -K -O -R0/11/0/8.5 -Jx1  >> rose2.ps
psrose rose2.d   -V -P -O    -: -A5r -S2.5i -G0/100/0 \
      -R0/350/0/360 -X-0.5i -Y0.25i -B50g50/30g10 >> rose2.ps
ghostview rose1.ps
ghostview rose2.ps
cp rose1.ps roseb1-$1.ps
cp rose2.ps roseb2-$1.ps

