#!/usr/bin/gnuplot
reset
set term png size 550, 450
set style line 1 lc 1 lt 1 lw 2
set title "Скорость добычи нефти"
set xlabel "Время, сек"
set ylabel "Объём в секунду, м^3/сек"
set output 'images/oil_production.png'
plot 'plotdat/oil_production.dat' with linespoints ls 1
