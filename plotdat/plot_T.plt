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
set title "Распределение насыщенности нефтью флюида от нагнетающей скважины до добывающей"
set ylabel "Насыщенность нефтью"
set output 'images/s_o.png'
plot	'result/s_o0.dat' using ($1):($2) with lines ls 1 title "s_o0", \
		'result/s_o1.dat' using ($1):($2) with lines ls 1 title "s_o1", \
		'result/s_o2.dat' using ($1):($2) with lines ls 1 title "s_o2", \
		'result/s_o3.dat' using ($1):($2) with lines ls 1 title "s_o3", \
		'result/s_o4.dat' using ($1):($2) with lines ls 1 title "s_o4", \
		'result/s_o5.dat' using ($1):($2) with lines ls 1 title "s_o5"
set title "Распределение скорости химической реакции от нагнетающей скважины до добывающей"
set ylabel "Скорость химической реакции"
set output 'images/rate_of_reaction.png'
plot	'result/rate_of_reaction0.dat' using ($1):($2) with lines ls 1 title "rate_of_reaction0", \
		'result/rate_of_reaction1.dat' using ($1):($2) with lines ls 1 title "rate_of_reaction1", \
		'result/rate_of_reaction2.dat' using ($1):($2) with lines ls 1 title "rate_of_reaction2", \
		'result/rate_of_reaction3.dat' using ($1):($2) with lines ls 1 title "rate_of_reaction3", \
		'result/rate_of_reaction4.dat' using ($1):($2) with lines ls 1 title "rate_of_reaction4", \
		'result/rate_of_reaction5.dat' using ($1):($2) with lines ls 1 title "rate_of_reaction5"
