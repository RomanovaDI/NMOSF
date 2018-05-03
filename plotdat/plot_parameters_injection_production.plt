#!/usr/bin/gnuplot
reset
set term pngcairo size 	800, 600 enhanced dashed font 'Helvetica,20'
#set term postscript eps size 	800, 600 enhanced dashed color font 'Helvetica,20'
#set term png size 550, 450
set style line 1 lc 3 lt 1 dt 1 lw 4
set xlabel "Расстояние, м"
#set title "Распределение температуры флюида от нагнетающей скважины до добывающей"
#set ylabel "Температура, К"
set output 'images/parameter0.png'
plot	'result/parameter0_step0.dat' using ($1):($2) with lines ls 1 title "0", \
		'result/parameter0_step1.dat' using ($1):($2) with lines ls 1 title "1", \
		'result/parameter0_step2.dat' using ($1):($2) with lines ls 1 title "2", \
		'result/parameter0_step3.dat' using ($1):($2) with lines ls 1 title "3", \
		'result/parameter0_step4.dat' using ($1):($2) with lines ls 1 title "4", \
		'result/parameter0_step5.dat' using ($1):($2) with lines ls 1 title "5"
set output 'images/parameter1.png'
plot	'result/parameter1_step0.dat' using ($1):($2) with lines ls 1 title "0", \
		'result/parameter1_step1.dat' using ($1):($2) with lines ls 1 title "1", \
		'result/parameter1_step2.dat' using ($1):($2) with lines ls 1 title "2", \
		'result/parameter1_step3.dat' using ($1):($2) with lines ls 1 title "3", \
		'result/parameter1_step4.dat' using ($1):($2) with lines ls 1 title "4", \
		'result/parameter1_step5.dat' using ($1):($2) with lines ls 1 title "5"
set output 'images/parameter2.png'
plot	'result/parameter2_step0.dat' using ($1):($2) with lines ls 1 title "0", \
		'result/parameter2_step1.dat' using ($1):($2) with lines ls 1 title "1", \
		'result/parameter2_step2.dat' using ($1):($2) with lines ls 1 title "2", \
		'result/parameter2_step3.dat' using ($1):($2) with lines ls 1 title "3", \
		'result/parameter2_step4.dat' using ($1):($2) with lines ls 1 title "4", \
		'result/parameter2_step5.dat' using ($1):($2) with lines ls 1 title "5"
set output 'images/parameter3.png'
plot	'result/parameter3_step0.dat' using ($1):($2) with lines ls 1 title "0", \
		'result/parameter3_step1.dat' using ($1):($2) with lines ls 1 title "1", \
		'result/parameter3_step2.dat' using ($1):($2) with lines ls 1 title "2", \
		'result/parameter3_step3.dat' using ($1):($2) with lines ls 1 title "3", \
		'result/parameter3_step4.dat' using ($1):($2) with lines ls 1 title "4", \
		'result/parameter3_step5.dat' using ($1):($2) with lines ls 1 title "5"
set output 'images/parameter4.png'
plot	'result/parameter4_step0.dat' using ($1):($2) with lines ls 1 title "0", \
		'result/parameter4_step1.dat' using ($1):($2) with lines ls 1 title "1", \
		'result/parameter4_step2.dat' using ($1):($2) with lines ls 1 title "2", \
		'result/parameter4_step3.dat' using ($1):($2) with lines ls 1 title "3", \
		'result/parameter4_step4.dat' using ($1):($2) with lines ls 1 title "4", \
		'result/parameter4_step5.dat' using ($1):($2) with lines ls 1 title "5"
set output 'images/parameter5.png'
plot	'result/parameter5_step0.dat' using ($1):($2) with lines ls 1 title "0", \
		'result/parameter5_step1.dat' using ($1):($2) with lines ls 1 title "1", \
		'result/parameter5_step2.dat' using ($1):($2) with lines ls 1 title "2", \
		'result/parameter5_step3.dat' using ($1):($2) with lines ls 1 title "3", \
		'result/parameter5_step4.dat' using ($1):($2) with lines ls 1 title "4", \
		'result/parameter5_step5.dat' using ($1):($2) with lines ls 1 title "5"
set output 'images/parameter6.png'
plot	'result/parameter6_step0.dat' using ($1):($2) with lines ls 1 title "0", \
		'result/parameter6_step1.dat' using ($1):($2) with lines ls 1 title "1", \
		'result/parameter6_step2.dat' using ($1):($2) with lines ls 1 title "2", \
		'result/parameter6_step3.dat' using ($1):($2) with lines ls 1 title "3", \
		'result/parameter6_step4.dat' using ($1):($2) with lines ls 1 title "4", \
		'result/parameter6_step5.dat' using ($1):($2) with lines ls 1 title "5"
set output 'images/parameter7.png'
plot	'result/parameter7_step0.dat' using ($1):($2) with lines ls 1 title "0", \
		'result/parameter7_step1.dat' using ($1):($2) with lines ls 1 title "1", \
		'result/parameter7_step2.dat' using ($1):($2) with lines ls 1 title "2", \
		'result/parameter7_step3.dat' using ($1):($2) with lines ls 1 title "3", \
		'result/parameter7_step4.dat' using ($1):($2) with lines ls 1 title "4", \
		'result/parameter7_step5.dat' using ($1):($2) with lines ls 1 title "5"
set output 'images/parameter8.png'
plot	'result/parameter8_step0.dat' using ($1):($2) with lines ls 1 title "0", \
		'result/parameter8_step1.dat' using ($1):($2) with lines ls 1 title "1", \
		'result/parameter8_step2.dat' using ($1):($2) with lines ls 1 title "2", \
		'result/parameter8_step3.dat' using ($1):($2) with lines ls 1 title "3", \
		'result/parameter8_step4.dat' using ($1):($2) with lines ls 1 title "4", \
		'result/parameter8_step5.dat' using ($1):($2) with lines ls 1 title "5"
set output 'images/parameter9.png'
plot	'result/parameter9_step0.dat' using ($1):($2) with lines ls 1 title "0", \
		'result/parameter9_step1.dat' using ($1):($2) with lines ls 1 title "1", \
		'result/parameter9_step2.dat' using ($1):($2) with lines ls 1 title "2", \
		'result/parameter9_step3.dat' using ($1):($2) with lines ls 1 title "3", \
		'result/parameter9_step4.dat' using ($1):($2) with lines ls 1 title "4", \
		'result/parameter9_step5.dat' using ($1):($2) with lines ls 1 title "5"
set output 'images/rate_of_reaction.png'
plot	'result/rate_of_reaction_step0.dat' using ($1):($2) with lines ls 1 title "0", \
		'result/rate_of_reaction_step1.dat' using ($1):($2) with lines ls 1 title "1", \
		'result/rate_of_reaction_step2.dat' using ($1):($2) with lines ls 1 title "2", \
		'result/rate_of_reaction_step3.dat' using ($1):($2) with lines ls 1 title "3", \
		'result/rate_of_reaction_step4.dat' using ($1):($2) with lines ls 1 title "4", \
		'result/rate_of_reaction_step5.dat' using ($1):($2) with lines ls 1 title "5"
