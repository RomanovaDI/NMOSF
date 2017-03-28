#!/usr/bin/gnuplot
reset
set term png
set title "Stencil of the difference scheme"
unset hidden3d
#set grid xtics ytics ztics front
set border
set output 'image.png'
set key off
set xrange[-2:1]
set yrange[-2:1]
set zrange[-2:1]
#$grid<<EOD
#-2 -2 -2 -2
#-2 -2 -2 -2
#-2 -2 -2 -2
#-2 -2 -2 -2
#EOD
#$gridd<<EOD
#-1 -1 -1 -1
#-1 -1 -1 -1
#-1 -1 -1 -1
#-1 -1 -1 -1
#EOD
#$griddd<<EOD
#0 0 0 0
#0 0 0 0
#0 0 0 0
#0 0 0 0
#EOD
#$grid<<EOD
#1 1 1 1
#1 1 1 1
#1 1 1 1
#1 1 1 1
#EOD
#set autoscale xfix
#set autoscale yfix
#set autoscale zfix
set xtics ("1" 1, "0" 0, "-1" -1, "-2" -2)
set ytics ("1" 1, "0" 0, "-1" -1, "-2" -2)
set ztics ("1" 1, "0" 0, "-1" -1, "-2" -2)
#set ytics 1
#set ztics 1
#set tics scale 0,0.001
#set mxtics 2
#set mytics 2
#set mztics 2
#set grid front mxtics mytics mztics lw 1.5 lt -1 lc rgb 'black'
splot 'momentum_eqn0.dat' using 2:3:4 with points palette pointsize 2 pointtype 1
