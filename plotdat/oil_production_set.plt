#!/usr/bin/gnuplot
reset
set term pngcairo size 	1200, 500 enhanced dashed font 'Helvetica,20'
#set term postscript eps size 	800, 600 enhanced dashed color font 'Helvetica,20'
#set term png size 550, 450
#set style line 1 lc 1 lt 1 dt 1 lw 4
#set style line 2 lc rgb '#00ad00' lt 1 dt 2 lw 4
#set style line 3 lc 3 lt 1 dt 3 lw 4
set style line 1 lc rgb '#000000' lt 1 dt 1 lw 4
set style line 2 lc rgb '#101010' lt 1 dt 2 lw 4
set style line 3 lc rgb '#202020' lt 1 dt 3 lw 4
set style line 4 lc rgb '#303030' lt 1 dt 4 lw 4
set style line 5 lc rgb '#404040' lt 1 dt 5 lw 4
set style line 6 lc rgb '#505050' lt 1 dt 6 lw 4
set style line 7 lc rgb '#606060' lt 1 dt 7 lw 4
set style line 8 lc rgb '#707070' lt 1 dt 8 lw 4
set style line 9 lc rgb '#808080' lt 1 dt 9 lw 4
set style line 10 lc rgb '#909090' lt 1 dt 10 lw 4
set style line 11 lc rgb '#a0a0a0' lt 1 dt 11 lw 4
set style line 12 lc rgb '#b0b0b0' lt 1 dt 12 lw 4
#set xlabel "Время, сек"
#set title "Скорость добычи нефти"
#set ylabel "Масса в секунду, кг/сек"
#set output 'images/velocity_of_oil_production_kg.png'
#plot 'result/velocity_of_oil_production_kg.dat' with lines ls 1 notitle
#set title "Масса вытесненной нефти"
#set ylabel "Масса, кг"
#set output 'images/oil_production_kg.png'
#plot 'result/oil_production_kg.dat' with lines ls 1 notitle
#set ylabel "Объём в секунду, м^3/сек"
#set output 'images/velocity_of_oil_production_m.png'
#plot 'result/velocity_of_oil_production_m.dat' with lines ls 1 notitle
#set title "Объём вытесненной нефти"
#set ylabel "Объём, м^3"
#set output 'images/oil_production_m.png'
#plot 'result/oil_production_m.dat' with lines ls 1 notitle
#
#
set xlabel "Время, сутки"
#set title "Поток массы нефти в час"
#set ylabel "Масса в секунду, кг/{час*м}"
#set output 'images/velocity_of_oil_production_kg.png'
#set logscale y
#plot	'termo_on/result/velocity_of_oil_production_kg.dat' using ($1/3600.0):($2*36.0) with lines ls 1 title "с термогазом", \
#		'termo_chem_off/result/velocity_of_oil_production_kg.dat' using ($1/3600.0):($2*36.0) with lines ls 2 title "температурное вытеснение", \
#		'termo_off/result/velocity_of_oil_production_kg.dat' using ($1/3600.0):($2*36.0) with lines ls 3 title "простое вытеснение"
#set title "Масса вытесненной нефти на метр мощности пласта"
#set ylabel "Масса, кг/м"
#set output 'images/oil_production_kg.png'
#unset logscale y
#set key right bottom
#plot	'termo_on/result/oil_production_kg.dat' using ($1/3600.0):($2/10.0) with lines ls 1 title "с термогазом", \
#		'termo_chem_off/result/oil_production_kg.dat' using ($1/3600.0):($2/10.0) with lines ls 2 title "температурное вытеснение", \
#		'termo_off/result/oil_production_kg.dat' using ($1/3600.0):($2/10.0) with lines ls 3 title "простое вытеснение"
#set title "Поток нефти на метр мощности пласта"
set notitle
set ylabel "Объём нефти в сутки, м^3/{сутки*м}      "
set grid
set output 'images/velocity_of_oil_production_m.png'
#set output 'images/velocity_of_oil_production_m.eps'
#set logscale y
set format y "%10.1e"
set key outside right top horizontal maxcolumns 1 width -50
#set key right top
plot	'termo_on1/result/velocity_of_oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 1 title "70% воды, 30% газовой смеси",\
		'termo_off/result/velocity_of_oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 2 title "простое вытеснение",\
		'air_less/result/velocity_of_oil_production_m.dat'  using ($1/86400.0):($2/10.0) with lines ls 3 title "80% 20%",\
		'air_less1/result/velocity_of_oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 4 title "75% 25%",\
		'air_less2/result/velocity_of_oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 5 title "85% 15%",\
		'air_more/result/velocity_of_oil_production_m.dat'  using ($1/86400.0):($2/10.0) with lines ls 6 title "60% 40%",\
		'air_more1/result/velocity_of_oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 7 title "65% 35%",\
		'rich/result/velocity_of_oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 8 title "rich",\
		'rich1/result/velocity_of_oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 9 title "rich1"
#		'termo_chem_off/result/velocity_of_oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 2 title "температурное вытеснение", \
#		'termo_off/result/velocity_of_oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 3 title "простое вытеснение"
#set title "Объём вытесненной нефти на метр мощности пласта"
set notitle
set ylabel "Объём, м^3/м"
set output 'images/oil_production_m.png'
#set output 'images/oil_production_m.eps'
#unset logscale y
#set key right bottom
unset format y
plot	'termo_on/result/oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 1 title "70% 30%", \
		'termo_off/result/oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 2 title "простое вытеснение",\
		'air_less/result/oil_production_m.dat'  using ($1/86400.0):($2/10.0) with lines ls 3 title "80% 20%",\
		'air_less1/result/oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 4 title "75% 25%",\
		'air_less2/result/oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 5 title "85% 15%",\
		'air_more/result/oil_production_m.dat'  using ($1/86400.0):($2/10.0) with lines ls 6 title "60% 40%",\
		'air_more1/result/oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 7 title "65% 35%",\
		'rich/result/oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 8 title "rich",\
		'rich1/result/oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 9 title "rich1"
#		'termo_chem_off/result/oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 2 title "температурное вытеснение", \
#	'termo_off/result/oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 3 title "простое вытеснение"
#set title "Объём вытесненного флюида на метр мощности пласта             "
set notitle
set output 'images/fluid_production_m.png'
#set output 'images/fluid_production_m.eps'
#set key width -52 left top at 10, 250 Left samplen 2
plot	'termo_on/result/fluid_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 1 title "70% 30%", \
		'termo_off/result/fluid_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 2 title "простое вытеснение",\
		'air_less/result/fluid_production_m.dat'  using ($1/86400.0):($2/10.0) with lines ls 3 title "80% 20%",\
		'air_less1/result/fluid_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 4 title "75% 25%",\
		'air_less2/result/fluid_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 5 title "85% 15%",\
		'air_more/result/fluid_production_m.dat'  using ($1/86400.0):($2/10.0) with lines ls 6 title "60% 40%",\
		'air_more1/result/fluid_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 7 title "65% 35%",\
		'rich/result/fluid_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 8 title "rich",\
		'rich1/result/fluid_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 9 title "rich1"
#		'termo_chem_off/result/fluid_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 2 title "температурное вытеснение", \
#	'termo_off/result/fluid_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 3 title "простое вытеснение"

#set xlabel "Time, hours"
#set title "Oil production rate per meter of reservoir power"
#set ylabel "Oil mass per hour, kg/{hour*m}"
#set output 'images/velocity_of_oil_production_chem_kg.png'
##set logscale y
#plot	'termo_on4_long/velocity_of_oil_production_kg.dat' using ($1/3600.0):($2*36.0) with lines ls 1 title "with termogas", \
#		'termo_off_long/velocity_of_oil_production_kg.dat' using ($1/3600.0):($2*36.0) with lines ls 2 title "plain displacement"
#set title "Mass of displaced oil per meter of reservoir power"
#set ylabel "Oil mass, kg/m"
#set output 'images/oil_production_chem_kg.png'
##unset logscale y
#set key right bottom
#plot	'termo_on4_long/oil_production_kg.dat' using ($1/3600.0):($2/10.0) with lines ls 1 title "with termogas", \
#		'termo_off_long/oil_production_kg.dat' using ($1/3600.0):($2/10.0) with lines ls 2 title "plain displacement"
#set title "Oil production rate per meter of reservoir power"
#set ylabel "Volume of oil per hour, m^3/{hour*m}"
#set output 'images/velocity_of_oil_production_chem_m.png'
##set logscale y
#set key right top
#plot	'termo_on4_long/velocity_of_oil_production_m.dat' using ($1/3600.0):($2*36.0) with lines ls 1 title "with termogas", \
#		'termo_off_long/velocity_of_oil_production_m.dat' using ($1/3600.0):($2*36.0) with lines ls 2 title "plain displacement"
#set title "Volume of displaced oil per meter of reservoir power"
#set ylabel "Volume, m^3/m"
#set output 'images/oil_production_chem_m.png'
##unset logscale y
#set key right bottom
#plot	'termo_on4_long/oil_production_m.dat' using ($1/3600.0):($2/10.0) with lines ls 1 title "with termogas", \
#		'termo_off_long/oil_production_m.dat' using ($1/3600.0):($2/10.0) with lines ls 2 title "plain displacement"
