#!/usr/bin/gnuplot
reset
set term pngcairo size 	800, 600 enhanced dashed font 'Helvetica,20'
#set term png size 550, 450
set style line 1 lc 1 lt 1 lw 4
set title "Скорость добычи нефти"
set xlabel "Время, сек"
set ylabel "Масса в секунду, кг/сек"
set output 'images/velocity_of_oil_production.png'
plot 'plotdat/velocity_of_oil_production.dat' with lines ls 1
set title "Масса вытесненной нефти"
set ylabel "Масса, кг"
set output 'images/oil_production.png'
plot 'plotdat/oil_production.dat' with lines ls 1
