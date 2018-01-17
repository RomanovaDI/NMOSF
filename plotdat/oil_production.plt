#!/usr/bin/gnuplot
reset
set term pngcairo size 	800, 600 enhanced dashed font 'Helvetica,20'
#set term png size 550, 450
set style line 1 lc 1 lt 1 lw 4
set title "Скорость добычи нефти"
set xlabel "Время, сек"
set ylabel "Масса в секунду, кг/сек"
set output 'images/velocity_of_oil_production_kg.png'
plot 'result_100000_no_term/velocity_of_oil_production_kg.dat' with lines ls 1 notitle
set title "Масса вытесненной нефти"
set ylabel "Масса, кг"
set output 'images/oil_production_kg.png'
plot 'result_100000_no_term/oil_production_kg.dat' with lines ls 1 notitle
