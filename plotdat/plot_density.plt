#!/usr/bin/gnuplot
reset
#set term postscript eps
#set output 'test.eps'
#test
#set term png size 550, 450
set term pngcairo size 550, 450
#set output 'test.png'
#test
set xlabel "Температура, К"
set ylabel "Плотность, кг/м^3"
set termoption dashed
set grid x y
set style line 1 lc 1 lt 1 dt 1 lw 3
set style line 2 lc 2 lt 2 dt 2 lw 3
set style line 3 lc 3 lt 3 dt 3 lw 3
#set linestyle  1 linetype 1 lc 1 lw 3
#set linestyle  2 linetype 2 lc 2 lw 3
#set linestyle  3 linetype 3 lc 3 lw 3
set title "Зависимость плотности воды от температуры"
set output 'images/density_water.png'
plot	[500:1200]	(998 + (20000000 - 1000000) / (1400 * 1400))/(1 + 0.0011 * (x - 330))	w l ls 1 title 'Плотность воды при p = 2*10^7 Па', \
					(998 + (30000000 - 1000000) / (1400 * 1400))/(1 + 0.0011 * (x - 330))	w l ls 2 title 'Плотность воды при p = 3*10^7 Па', \
					(998 + (1000000 - 1000000) / (1400 * 1400))/(1 + 0.0011 * (x - 330))	w l ls 3 title 'Плотность воды при p = 10^6 Па'
set title "Зависимость плотности нефти от температуры"
set output 'images/density_oil.png'
plot	[500:1200]	(850 + (20000000 - 1000000) / (1300 * 1300))/(1 + 0.0008701 * (x - 330))	w l ls 1 title 'Плотность нефти при p = 2*10^7 Па', \
					(850 + (30000000 - 1000000) / (1300 * 1300))/(1 + 0.0008701 * (x - 330))	w l ls 2 title 'Плотность нефти при p = 3*10^7 Па', \
					(850 + (1000000 - 1000000) / (1300 * 1300))/(1 + 0.0008701 * (x - 330))		w l ls 3 title 'Плотность нефти при p = 10^6 Па'
set title "Зависимость плотности нефти от температуры"
set output 'images/density_oil_incorrect.png'
plot	[500:1200]	(850 + (20000000 - 1000000) / (1300 * 1300))/(1 + 8701 * (x - 330))	w l ls 1 title 'Плотность нефти при p = 2*10^7 Па', \
					(850 + (30000000 - 1000000) / (1300 * 1300))/(1 + 8701 * (x - 330))	w l ls 2 title 'Плотность нефти при p = 3*10^7 Па', \
					(850 + (1000000 - 1000000) / (1300 * 1300))/(1 + 8701 * (x - 330))	w l ls 3 title 'Плотность нефти при p = 10^6 Па'
set title "Зависимость плотности газа от температуры"
set output 'images/density_gas_temperature.png'
plot	[500:1200]	20000000 * 0.25 * (0.02801 + 0.032 + 0.04401 + 0.01802) / (8.314 * x)	w l ls 1 title 'Плотность газа при p = 2*10^7 Па', \
					30000000 * 0.25 * (0.02801 + 0.032 + 0.04401 + 0.01802) / (8.314 * x)	w l ls 2 title 'Плотность газа при p = 3*10^7 Па', \
					1000000 * 0.25 * (0.02801 + 0.032 + 0.04401 + 0.01802) / (8.314 * x)	w l ls 3 title 'Плотность газа при p = 10^6 Па'
set title "Зависимость плотности газа от давления"
set output 'images/density_gas_pressure.png'
set xlabel "Давление, Па"
plot	[100000:30000000]	x * 0.25 * (0.02801 + 0.032 + 0.04401 + 0.01802) / (8.314 * 200)	w l ls 1 title 'Плотность газа при T = 200 К', \
							x * 0.25 * (0.02801 + 0.032 + 0.04401 + 0.01802) / (8.314 * 450)	w l ls 2 title 'Плотность газа при T = 450 К', \
							x * 0.25 * (0.02801 + 0.032 + 0.04401 + 0.01802) / (8.314 * 700)	w l ls 3 title 'Плотность газа при T = 700 К'
