#!/usr/bin/gnuplot
reset
set term pngcairo size 	800, 600 enhanced dashed font 'Helvetica,20'
#set term postscript eps size 	800, 600 enhanced dashed color font 'Helvetica,20'
#set term png size 550, 450
set style line 6 lc rgb "#000000" lt 1 dt 1 lw 2
set style line 5 lc rgb "#202020" lt 1 dt 1 lw 4
set style line 4 lc rgb "#404040" lt 1 dt 1 lw 6
set style line 3 lc rgb "#606060" lt 1 dt 1 lw 8
set style line 2 lc rgb "#808080" lt 1 dt 1 lw 10
set style line 1 lc rgb "#a0a0a0" lt 1 dt 1 lw 12
set xlabel "Расстояние, м"
#set title "Распределение температуры флюида от нагнетающей скважины до добывающей"
#set ylabel "Температура, К"
set samples 300
set output 'images/parameter0.png'
set key right center
plot	'result/parameter0_step0.dat' using ($1):($2) smooth csplines with lines ls 1 title "0 день", \
		'result/parameter0_step1.dat' using ($1):($2) smooth csplines with lines ls 2 title "73 день", \
		'result/parameter0_step2.dat' using ($1):($2) smooth csplines with lines ls 3 title "146 день", \
		'result/parameter0_step3.dat' using ($1):($2) smooth csplines with lines ls 4 title "219 день", \
		'result/parameter0_step4.dat' using ($1):($2) smooth csplines with lines ls 5 title "292 день", \
		'result/parameter0_step5.dat' using ($1):($2) smooth csplines with lines ls 6 title "365 день"
set output 'images/parameter1.png'
set key right top
plot	'result/parameter1_step0.dat' using ($1):($2) smooth csplines with lines ls 1 title "0 день", \
		'result/parameter1_step1.dat' using ($1):($2) smooth csplines with lines ls 2 title "73 день", \
		'result/parameter1_step2.dat' using ($1):($2) smooth csplines with lines ls 3 title "146 день", \
		'result/parameter1_step3.dat' using ($1):($2) smooth csplines with lines ls 4 title "219 день", \
		'result/parameter1_step4.dat' using ($1):($2) smooth csplines with lines ls 5 title "292 день", \
		'result/parameter1_step5.dat' using ($1):($2) smooth csplines with lines ls 6 title "365 день"
set output 'images/parameter2.png'
set key left top
plot	'result/parameter2_step0.dat' using ($1):($2) smooth csplines with lines ls 1 title "0 день", \
		'result/parameter2_step1.dat' using ($1):($2) smooth csplines with lines ls 2 title "73 день", \
		'result/parameter2_step2.dat' using ($1):($2) smooth csplines with lines ls 3 title "146 день", \
		'result/parameter2_step3.dat' using ($1):($2) smooth csplines with lines ls 4 title "219 день", \
		'result/parameter2_step4.dat' using ($1):($2) smooth csplines with lines ls 5 title "292 день", \
		'result/parameter2_step5.dat' using ($1):($2) smooth csplines with lines ls 6 title "365 день"
set output 'images/parameter3.png'
set key right bottom
plot	'result/parameter3_step0.dat' using ($1):($2) smooth csplines with lines ls 1 title "0 день", \
		'result/parameter3_step1.dat' using ($1):($2) smooth csplines with lines ls 2 title "73 день", \
		'result/parameter3_step2.dat' using ($1):($2) smooth csplines with lines ls 3 title "146 день", \
		'result/parameter3_step3.dat' using ($1):($2) smooth csplines with lines ls 4 title "219 день", \
		'result/parameter3_step4.dat' using ($1):($2) smooth csplines with lines ls 5 title "292 день", \
		'result/parameter3_step5.dat' using ($1):($2) smooth csplines with lines ls 6 title "365 день"
set output 'images/parameter4.png'
set key right top
plot	'result/parameter4_step0.dat' using ($1):($2) with lines ls 1 title "0 день", \
		'result/parameter4_step1.dat' using ($1):($2) with lines ls 2 title "73 день", \
		'result/parameter4_step2.dat' using ($1):($2) with lines ls 3 title "146 день", \
		'result/parameter4_step3.dat' using ($1):($2) with lines ls 4 title "219 день", \
		'result/parameter4_step4.dat' using ($1):($2) with lines ls 5 title "292 день", \
		'result/parameter4_step5.dat' using ($1):($2) with lines ls 6 title "365 день"
set output 'images/parameter5.png'
set key right top
plot	'result/parameter5_step0.dat' using ($1):($2) smooth csplines with lines ls 1 title "0 день", \
		'result/parameter5_step1.dat' using ($1):($2) smooth csplines with lines ls 2 title "73 день", \
		'result/parameter5_step2.dat' using ($1):($2) smooth csplines with lines ls 3 title "146 день", \
		'result/parameter5_step3.dat' using ($1):($2) smooth csplines with lines ls 4 title "219 день", \
		'result/parameter5_step4.dat' using ($1):($2) smooth csplines with lines ls 5 title "292 день", \
		'result/parameter5_step5.dat' using ($1):($2) smooth csplines with lines ls 6 title "365 день"
set output 'images/parameter6.png'
set key right bottom
plot	'result/parameter6_step0.dat' using ($1):($2) smooth csplines with lines ls 1 title "0 день", \
		'result/parameter6_step1.dat' using ($1):($2) smooth csplines with lines ls 2 title "73 день", \
		'result/parameter6_step2.dat' using ($1):($2) smooth csplines with lines ls 3 title "146 день", \
		'result/parameter6_step3.dat' using ($1):($2) smooth csplines with lines ls 4 title "219 день", \
		'result/parameter6_step4.dat' using ($1):($2) smooth csplines with lines ls 5 title "292 день", \
		'result/parameter6_step5.dat' using ($1):($2) smooth csplines with lines ls 6 title "365 день"
set output 'images/parameter7.png'
set key right top
plot	'result/parameter7_step0.dat' using ($1):($2) smooth csplines with lines ls 1 title "0 день", \
		'result/parameter7_step1.dat' using ($1):($2) smooth csplines with lines ls 2 title "73 день", \
		'result/parameter7_step2.dat' using ($1):($2) smooth csplines with lines ls 3 title "146 день", \
		'result/parameter7_step3.dat' using ($1):($2) smooth csplines with lines ls 4 title "219 день", \
		'result/parameter7_step4.dat' using ($1):($2) smooth csplines with lines ls 5 title "292 день", \
		'result/parameter7_step5.dat' using ($1):($2) smooth csplines with lines ls 6 title "365 день"
set output 'images/parameter8.png'
set key right top
plot	'result/parameter8_step0.dat' using ($1):($2) with lines ls 1 title "0 день", \
		'result/parameter8_step1.dat' using ($1):($2) smooth csplines with lines ls 2 title "73 день", \
		'result/parameter8_step2.dat' using ($1):($2) smooth csplines with lines ls 3 title "146 день", \
		'result/parameter8_step3.dat' using ($1):($2) smooth csplines with lines ls 4 title "219 день", \
		'result/parameter8_step4.dat' using ($1):($2) smooth csplines with lines ls 5 title "292 день", \
		'result/parameter8_step5.dat' using ($1):($2) smooth csplines with lines ls 6 title "365 день"
set output 'images/parameter9.png'
set key right top
plot	'result/parameter9_step0.dat' using ($1):($2) with lines ls 1 title "0 день", \
		'result/parameter9_step1.dat' using ($1):($2) smooth csplines with lines ls 2 title "73 день", \
		'result/parameter9_step2.dat' using ($1):($2) smooth csplines with lines ls 3 title "146 день", \
		'result/parameter9_step3.dat' using ($1):($2) smooth csplines with lines ls 4 title "219 день", \
		'result/parameter9_step4.dat' using ($1):($2) smooth csplines with lines ls 5 title "292 день", \
		'result/parameter9_step5.dat' using ($1):($2) smooth csplines with lines ls 6 title "365 день"
set output 'images/rate_of_reaction.png'
set key right top
plot	'result/rate_of_reaction_step0.dat' using ($1):($2) with lines ls 1 title "0 день", \
		'result/rate_of_reaction_step1.dat' using ($1):($2) with lines ls 2 title "73 день", \
		'result/rate_of_reaction_step2.dat' using ($1):($2) with lines ls 3 title "146 день", \
		'result/rate_of_reaction_step3.dat' using ($1):($2) with lines ls 4 title "219 день", \
		'result/rate_of_reaction_step4.dat' using ($1):($2) with lines ls 5 title "292 день", \
		'result/rate_of_reaction_step5.dat' using ($1):($2) with lines ls 6 title "365 день"
