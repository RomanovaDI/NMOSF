#!/usr/bin/gnuplot
reset
set term pngcairo size 	1000, 300 enhanced dashed font 'Helvetica,20'
set style data histogram
set style histogram rowstacked
set key invert reverse Left outside
set yrange [0:1]
unset xtics
set style fill solid border
set boxwidth 0.65
set output 'images/hist2.png'
plot	'result/p5s2.dat' using ($2):xtic(1) title "water",\
		'result/p6s2.dat' using ($2) title "oil",\
		'result/p7s2.dat' using ($2) title "gas"
