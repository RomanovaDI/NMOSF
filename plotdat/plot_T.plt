#!/usr/bin/gnuplot
reset
set term pngcairo size 	800, 600 enhanced dashed font 'Helvetica,20'
#set term postscript eps size 	800, 600 enhanced dashed color font 'Helvetica,20'
#set term png size 550, 450
set style line 1 lc 3 lt 1 dt 1 lw 4
set xlabel "Расстояние, м"
set title "Распределение температуры флюида от нагнетающей скважины до добывающей"
set ylabel "Температура, К"
set output 'images/T.png'
plot	'result/T0.dat' using ($1):($2) with lines ls 1 title "T0", \
		'result/T1.dat' using ($1):($2) with lines ls 1 title "T1", \
		'result/T2.dat' using ($1):($2) with lines ls 1 title "T2", \
		'result/T3.dat' using ($1):($2) with lines ls 1 title "T3", \
		'result/T4.dat' using ($1):($2) with lines ls 1 title "T4", \
		'result/T5.dat' using ($1):($2) with lines ls 1 title "T5"
