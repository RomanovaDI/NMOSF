#!/usr/bin/gnuplot
reset
set term pngcairo size 	800, 600 enhanced dashed font 'Helvetica,20'
#set term png size 550, 450
set style line 1 lc 1 lt 1 dt 1 lw 4
set style line 2 lc 2 lt 1 dt 2 lw 4
set style line 3 lc 3 lt 1 dt 3 lw 4
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
set title "Поток нефти на метр мощности пласта"
set ylabel "Объём нефти в сутки, м^3/{сутки*м}"
set output 'images/velocity_of_oil_production_m.png'
#set logscale y
set key right top
plot	'termo_on/result/velocity_of_oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 1 title "с термогазом", \
		'termo_chem_off/result/velocity_of_oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 2 title "температурное вытеснение", \
		'termo_off/result/velocity_of_oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 3 title "простое вытеснение"
set title "Объём вытесненной нефти на метр мощности пласта"
set ylabel "Объём, м^3/м"
set output 'images/oil_production_m.png'
#unset logscale y
set key right bottom
plot	'termo_on/result/oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 1 title "с термогазом", \
		'termo_chem_off/result/oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 2 title "температурное вытеснение", \
		'termo_off/result/oil_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 3 title "простое вытеснение"
set title "Объём вытесненного флюида на метр мощности пласта             "
set output 'images/fluid_production_m.png'
#set key width -52 left top at 10, 250 Left samplen 2
plot	'termo_on/result/fluid_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 1 title "с термогазом", \
		'termo_chem_off/result/fluid_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 2 title "температурное вытеснение", \
		'termo_off/result/fluid_production_m.dat' using ($1/86400.0):($2/10.0) with lines ls 3 title "простое вытеснение"

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
