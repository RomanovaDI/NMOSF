#!/usr/bin/gnuplot
reset
set term pngcairo size 	800, 600 enhanced dashed font 'Helvetica,20'
#set terminal postscript eps size 3.5,2.62 enhanced color
#set output "images/test.svg"
#test
#reset
#set term pngcairo enhanced size 550, 450
set title "Фазовые коэффициенты вязкости"
set xlabel "Плотность, кг/м^3"
set ylabel "Вязкость, Па*с"
#set termoption dashed
set grid x y
set linestyle  1 lc 1 lt 1 dt 1 lw 4
set linestyle  2 lc 2 lt 1 dt 2 lw 4
set linestyle  3 lc 3 lt 1 dt 3 lw 4
set linestyle  4 lc 4 lt 1 dt 4 lw 4
set format y "%1.6f"
set title "Фазовые коэффициенты вязкости"
set output 'images/viscosity_water_belova.png'
plot [500:900]		0.15566 * 0.0000001 / (1 / x - 0.000984)	w l ls 1 title "Вязкость воды"
set output 'images/viscosity_oil_belova.png'
plot [450:800]		3.83 * 0.0000001 / (1 / x - 0.00117)		w l ls 1 title "Вязкость нефти"
set output 'images/viscosity_water.png'
plot [500:900]		0.15566 * 0.0000001 / (1 / x - 0.0000984)	w l ls 1 title "Вязкость воды"
set output 'images/viscosity_oil.png'
plot [450:800]		3.83 * 0.0000001 / (1 / x - 0.000117)		w l ls 1 title "Вязкость нефти"
set xlabel "Температура, К"
set title "Фазовые коэффициенты вязкости при p = 2*10^7 Па       "
set output 'images/viscosity_water_temperature_belova.png'
plot [500:1200]		0.15566 * 0.0000001 / (1 / ((998 + (20000000 - 1000000) / (1400 * 1400))/(1 + 0.0011 * (x - 330))) - 0.000984)	w l ls 1 title 'Вязкость воды'
set output 'images/viscosity_oil_temperature_belova.png'
plot [500:1200]		3.83 * 0.0000001 / (1 / ((850 + (20000000 - 1000000) / (1300 * 1300))/(1 + 0.0008701 * (x - 330))) - 0.00117)	w l ls 1 title 'Вязкость нефти'
set output 'images/viscosity_water_temperature.png'
plot [380:700]		0.15566 * 0.0000001 / (1 / ((998 + (20000000 - 1000000) / (1400 * 1400))/(1 + 0.0011 * (x - 330))) - 0.0000984)	w l ls 1 title 'Вязкость воды'
set output 'images/viscosity_oil_temperature.png'
plot [380:700]		3.83 * 0.0000001 / (1 / ((850 + (20000000 - 1000000) / (1300 * 1300))/(1 + 0.0008701 * (x - 330))) - 0.000117)	w l ls 1 title 'Вязкость нефти'
set output 'images/viscosity_water_temperature_pressure.png'
set title "Вязкость воды"
plot	[380:700]	0.15566 * 0.0000001 / (1 / ((998 + (20000000 - 1000000) / (1400 * 1400))/(1 + 0.0011 * (x - 300))) - 0.0000984)	w l ls 1 title 'при p = 2*10^7 Па', \
					0.15566 * 0.0000001 / (1 / ((998 + (30000000 - 1000000) / (1400 * 1400))/(1 + 0.0011 * (x - 300))) - 0.0000984)	w l ls 2 title 'при p = 3*10^7 Па', \
					0.15566 * 0.0000001 / (1 / ((998 + (1000000 - 1000000) / (1400 * 1400))/(1 + 0.0011 * (x - 300))) - 0.0000984)	w l ls 3 title 'при p = 10^6 Па'
set output 'images/viscosity_oil_temperature_pressure.png'
set title "Вязкость нефти"
plot	[380:700]	3.83 * 0.0000001 / (1 / ((850 + (20000000 - 1000000) / (1300 * 1300))/(1 + 0.0008701 * (x - 300))) - 0.000117)	w l ls 1 title 'при p = 2*10^7 Па', \
					3.83 * 0.0000001 / (1 / ((850 + (30000000 - 1000000) / (1300 * 1300))/(1 + 0.0008701 * (x - 300))) - 0.000117)	w l ls 2 title 'при p = 3*10^7 Па', \
					3.83 * 0.0000001 / (1 / ((850 + (1000000 - 1000000) / (1300 * 1300))/(1 + 0.0008701 * (x - 300))) - 0.000117)	w l ls 3 title 'при p = 10^6 Па'
set output 'images/viscosity_water_temperature_byrge.png'
set title "Вязкость воды"
plot [380:700]		0.001 / (0.14 + (x - 273.15) / 30.0 + 0.000009 * (x - 273.15) * (x - 273.15))	w l ls 1 title 'Вязкость воды'
set output 'images/viscosity_oil_temperature_byrge.png'
set title "Вязкость нефти"
#set logscale y
#plot [280:600]		((850 + (20000000 - 1000000) / (1300 * 1300))/(1 + 0.0008701 * (x - 300))) * ((17.647 + 0.6) ** (((300 - 273.15) / (x - 273.15)) ** 4) - 0.6) * 0.000001	w l ls 1 title 'Вязкость нефти при p = 2*10^7 Па', \
#plot [310:600]		((850.0 + (20000000.0 - 1000000.0) / (1300.0 * 1300.0)) / (1.0 + 0.0008701 * (x - 330.0))) * 1000000.0 * ((200.0 + 0.6) ** ((30.0 / (x - 273.15)) ** 4.0) - 0.6)	w l ls 1 title 'Вязкость нефти при p = 2*10^7 Па'
plot [380:700]		((850.0 + (20000000.0 - 1000000.0) / (1300.0 * 1300.0)) / (1.0 + 0.0008701 * (x - 330.0))) * 0.000001 * ((1000.0 + 0.6) ** ((30.0 / (x - 273.15)) ** 4.0) - 0.6)	w l ls 1 title 'Вязкость нефти при p = 2*10^7 Па'
set output 'images/kinematic_viscosity_oil_temperature_byrge.png'
set ylabel "Вязкость, м^2/с"
plot [300:600]		(((17.647 + 0.6) ** (((300 - 273.15) / (x - 273.15)) ** 4)) - 0.6)	w l ls 1 title 'Вязкость нефти'
unset logscale y
set title "Вязкость газовой смеси при равной концентрации компонент                  "
set ylabel "Вязкость, Па*с"
set output 'images/viscosity_gas_mixture.png'
plot [380:700]	(0.00001781 * (300.55 + 111) / (x + 111) * ((x / 300.55) ** (3. / 2.))) ** (0.0068885923381663767 / 0.02801) *\
		 		(0.00002018 * (292.25 + 127) / (x + 127) * ((x / 292.25) ** (3. / 2.))) ** (0.0068885923381663767 / 0.032) *\
				(0.0000148 * (293.15  + 240) / (x + 240) * ((x / 293.15) ** (3. / 2.))) ** (0.0068885923381663767 / 0.04401) *\
				(0.000018 * 330 / x * ((x / 330) ** (3. / 2.))) ** (0.0068885923381663767 / 0.01802)		w l ls 1 title 'Вязкость газа'
set output 'images/viscosity_gas_mixture_belova.png'
plot [500:1200]	(0.000018 * (330 + 464) / (x + 464) * ((x / 330) ** (3. / 2.))) ** (0.0068885923381663767 / 0.02801) *\
		 		(0.000021 * (330 + 613) / (x + 613) * ((x / 330) ** (3. / 2.))) ** (0.0068885923381663767 / 0.032) *\
				(0.000018 * (330 - 1) / (x - 1) * ((x / 330) ** (3. / 2.))) ** (0.0068885923381663767 / 0.04401) *\
				(0.000018 * (330 - 89) / (x - 89) * ((x / 330) ** (3. / 2.))) ** (0.0068885923381663767 / 0.01802)		w l ls 1 title 'Вязкость газа'
set title "Вязкости газовых компонент"
set key width -53 left top at 600, 0.00008 Left samplen 2
#set key box
#set key at 1000, 0.000078 right top
#set key left top Left
#set key horiz
set output 'images/viscosity_gas_belova.png'
plot [500:1200]	0.000018 * (330 + 464) / (x + 464) * sqrt((x / 330) ** 3)	w l ls 1 title "Вязкость азота", \
		 		0.000021 * (330 + 613) / (x + 613) * sqrt((x / 330) ** 3)	w l ls 2 title "Вязкость кислорода", \
				0.000018 * (330 - 1) / (x - 1) * sqrt((x / 330) ** 3)		w l ls 3 title "Вязкость углекислого газа", \
				0.000018 * (330 - 89) / (x - 89) * sqrt((x / 330) ** 3)		w l ls 4 title "Вязкость водяного пара"
set key width -53 left top at 420, 0.000038 Left samplen 2
set output 'images/viscosity_gas.png'
plot [380:700]	0.00001781 * (300.55 + 111) / (x + 111) * ((x / 300.55) ** (3. / 2.))	w l ls 1 title "Вязкость азота", \
		 		0.00002018 * (292.25 + 127) / (x + 127) * ((x / 292.25) ** (3. / 2.))	w l ls 2 title "Вязкость кислорода", \
				0.0000148 * (293.15  + 240) / (x + 240) * ((x / 293.15) ** (3. / 2.))	w l ls 3 title "Вязкость углекислого газа", \
				0.000018 * 330 / x * ((x / 330) ** (3. / 2.))							w l ls 4 title "Вязкость водяного пара"
set output 'images/viscosity_gas_byrge.png'
plot [380:700]	0.001 * (170.0 + 0.38 * (x - 273.15))	w l ls 1 title "Вязкость азота", \
		 		0.001 * (188.0 + 0.38 * (x - 273.15))	w l ls 2 title "Вязкость кислорода", \
				0.001 * (145.0 + 0.38 * (x - 273.15))	w l ls 3 title "Вязкость углекислого газа", \
				0.001 * (88.0 + 0.38 * (x - 273.15))	w l ls 4 title "Вязкость водяного пара"
