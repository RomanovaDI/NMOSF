#!/usr/bin/gnuplot
reset
set term pngcairo size 550, 450
set title "Фазовые коэффициенты вязкости"
set xlabel "Плотность, кг/м^3"
set ylabel "Вязкость, Па*с"
set termoption dashed
set grid x y
set linestyle  1 lc 1 lt 1 dt 1 lw 3
set linestyle  2 lc 2 lt 1 dt 2 lw 3
set linestyle  3 lc 3 lt 1 dt 3 lw 3
set linestyle  4 lc 4 lt 1 dt 4 lw 3
set format y "%1.6f"
set title "Фазовые коэффициенты вязкости"
set output 'images/viscosity_water_belova.png'
plot [500:900]		0.15566 * 0.0000001 / (1 / x - 0.000984)	w l ls 1 title 'Вязкость воды'
set output 'images/viscosity_oil_belova.png'
plot [450:800]		3.83 * 0.0000001 / (1 / x - 0.00117)		w l ls 1 title 'Вязкость нефти'
set output 'images/viscosity_water.png'
plot [500:900]		0.15566 * 0.0000001 / (1 / x - 0.0000984)	w l ls 1 title 'Вязкость воды'
set output 'images/viscosity_oil.png'
plot [450:800]		3.83 * 0.0000001 / (1 / x - 0.000117)		w l ls 1 title 'Вязкость нефти'
set xlabel "Температура, К"
set title "Фазовые коэффициенты вязкости при p = 2*10^7 Па"
set output 'images/viscosity_water_temperature_belova.png'
plot [500:1200]		0.15566 * 0.0000001 / (1 / ((998 + (20000000 - 1000000) / (1400 * 1400))/(1 + 0.0011 * (x - 330))) - 0.000984)	w l ls 1 title 'Вязкость воды'
set output 'images/viscosity_oil_temperature_belova.png'
plot [500:1200]		3.83 * 0.0000001 / (1 / ((850 + (20000000 - 1000000) / (1300 * 1300))/(1 + 0.0008701 * (x - 330))) - 0.00117)	w l ls 1 title 'Вязкость нефти'
set output 'images/viscosity_water_temperature.png'
plot [500:1200]		0.15566 * 0.0000001 / (1 / ((998 + (20000000 - 1000000) / (1400 * 1400))/(1 + 0.0011 * (x - 330))) - 0.0000984)	w l ls 1 title 'Вязкость воды'
set output 'images/viscosity_oil_temperature.png'
plot [500:1200]		3.83 * 0.0000001 / (1 / ((850 + (20000000 - 1000000) / (1300 * 1300))/(1 + 0.0008701 * (x - 330))) - 0.000117)	w l ls 1 title 'Вязкость нефти'
set title "Вязкости газовых компонент"
set output 'images/viscosity_gas_belova.png'
plot [500:1200]	0.000018 * (330 + 464) / (x + 464) * sqrt((x / 330) ** 3)	w l ls 1 title 'Вязкость азота', \
		 		0.000021 * (330 + 613) / (x + 613) * sqrt((x / 330) ** 3)	w l ls 2 title 'Вязкость кислорода', \
				0.000018 * (330 - 1) / (x - 1) * sqrt((x / 330) ** 3)		w l ls 3 title 'Вязкость углекислого газа', \
				0.000018 * (330 - 89) / (x - 89) * sqrt((x / 330) ** 3)		w l ls 4 title 'Вязкость водяного пара'
set output 'images/viscosity_gas.png'
plot [500:1200]	0.00001781 * (300.55 + 111) / (x + 111) * ((x / 300.55) ** (3. / 2.))	w l ls 1 title 'Вязкость азота', \
		 		0.00002018 * (292.25 + 127) / (x + 127) * ((x / 292.25) ** (3. / 2.))	w l ls 2 title 'Вязкость кислорода', \
				0.0000148 * (293.15  + 240) / (x + 240) * ((x / 293.15) ** (3. / 2.))	w l ls 3 title 'Вязкость углекислого газа', \
				0.000018 * 330 / x * ((x / 330) ** (3. / 2.))			w l ls 4 title 'Вязкость водяного пара'
set title "Вязкость газовой смеси при равной концентрации компонент                  "
set output 'images/viscosity_gas_mixture.png'
plot [500:1200]	(0.00001781 * (300.55 + 111) / (x + 111) * ((x / 300.55) ** (3. / 2.))) ** (0.0068885923381663767 / 0.02801) *\
		 		(0.00002018 * (292.25 + 127) / (x + 127) * ((x / 292.25) ** (3. / 2.))) ** (0.0068885923381663767 / 0.032) *\
				(0.0000148 * (293.15  + 240) / (x + 240) * ((x / 293.15) ** (3. / 2.))) ** (0.0068885923381663767 / 0.04401) *\
				(0.000018 * 330 / x * ((x / 330) ** (3. / 2.))) ** (0.0068885923381663767 / 0.01802)		w l ls 1 title 'Вязкость газа'
set output 'images/viscosity_gas_mixture_belova.png'
plot [500:1200]	(0.000018 * (330 + 464) / (x + 464) * ((x / 330) ** (3. / 2.))) ** (0.0068885923381663767 / 0.02801) *\
		 		(0.000021 * (330 + 613) / (x + 613) * ((x / 330) ** (3. / 2.))) ** (0.0068885923381663767 / 0.032) *\
				(0.000018 * (330 - 1) / (x - 1) * ((x / 330) ** (3. / 2.))) ** (0.0068885923381663767 / 0.04401) *\
				(0.000018 * (330 - 89) / (x - 89) * ((x / 330) ** (3. / 2.))) ** (0.0068885923381663767 / 0.01802)		w l ls 1 title 'Вязкость газа'
