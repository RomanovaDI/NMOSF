#!/usr/bin/gnuplot
reset
set term png size 1920, 1080
set output 'image.png'
set xrange[-2:1]
set yrange[-2:1]
set zrange[-2:1]
#set bmargin -2
splot	'test_points.dat' using 1:2:3 with lines ls 1 lc rgb "blue"
