#!/usr/bin/gnuplot
reset
set term png size 550, 450
set style line 1 lc 1 lt 1 lw 2
set title "Скорость добычи нефти"
set xlabel "Время, сек"
set ylabel "Объём в секунду, м^3/сек"
set output 'images/velocity_of_oil_production.png'
plot 'plotdat/velocity_of_oil_production.dat' with lines ls 1
set title "Объём вытесненной нефти"
set ylabel "Объём, м^3"
set output 'images/oil_production.png'
plot 'plotdat/oil_production.dat' with lines ls 1
