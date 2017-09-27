#!/usr/bin/gnuplot
reset
set term png size 640, 480
set title "Фазовые коэффициенты вязкости"
set xlabel "Плотность, кг/м^3"
set ylabel "Вязкость, Па*с"
set termoption dashed
set linestyle  1 linetype 1 lc 1 lw 3
set linestyle  2 linetype 2 lc 2 lw 3
set linestyle  3 linetype 3 lc 3 lw 3
set linestyle  4 linetype 4 lc 4 lw 3
set format y "%1.6f"
set output 'images/viscosity_belova.png'
plot	[0.001:1000]	0.15566 * 0.0000001 / (1 / x - 0.000984)	w l ls 1 title 'Вязкость воды', \
						3.83 * 0.0000001 / (1 / x - 0.00117)		w l ls 2 title 'Вязкость нефти', \
						(0.00001781 * (300.55 + 111) / (x + 111) * ((x / 300.55) ** (3. / 2.))) ** (0.8194034742707309 / 0.02801) *\
		 				(0.00002018 * (292.25 + 127) / (x + 127) * ((x / 292.25) ** (3. / 2.))) ** (0.8194034742707309 / 0.032) *\
						(0.0000148 * (293.15  + 240) / (x + 240) * ((x / 293.15) ** (3. / 2.))) ** (0.8194034742707309 / 0.04401) *\
						(0.000018 * 330 / x * ((x / 330) ** (3. / 2.))) ** (0.8194034742707309 / 0.01802)		w l ls 3 title 'Вязкость газа'
set output 'images/viscosity.png'
plot	[0.001:1000]	1 / ((1 / 0.15566 * 0.0000001) / x - 0.000984 / (0.15566 * 0.0000001))	w l ls 1 title 'Вязкость воды', \
						1 / ((1 / 3.83 * 0.0000001) / x - 0.00117 / (3.83 * 0.0000001))			w l ls 2 title 'Вязкость нефти', \
						(0.00001781 * (300.55 + 111) / (x + 111) * ((x / 300.55) ** (3. / 2.))) ** (0.8194034742707309 / 0.02801) *\
		 				(0.00002018 * (292.25 + 127) / (x + 127) * ((x / 292.25) ** (3. / 2.))) ** (0.8194034742707309 / 0.032) *\
						(0.0000148 * (293.15  + 240) / (x + 240) * ((x / 293.15) ** (3. / 2.))) ** (0.8194034742707309 / 0.04401) *\
						(0.000018 * 330 / x * ((x / 330) ** (3. / 2.))) ** (0.8194034742707309 / 0.01802)		w l ls 3 title 'Вязкость газа'
set output 'images/viscosity_gas_belova.png'
set xlabel "Температура, К"
plot [200:700]	0.000018 * (330 + 464) / (x + 464) * sqrt((x / 330) ** 3)	w l ls 1 title 'Вязкость азота', \
		 		0.000021 * (330 + 613) / (x + 613) * sqrt((x / 330) ** 3)	w l ls 2 title 'Вязкость кислорода', \
				0.000018 * (330 - 1) / (x - 1) * sqrt((x / 330) ** 3)		w l ls 3 title 'Вязкость углекислого газа', \
				0.000018 * (330 - 89) / (x - 89) * sqrt((x / 330) ** 3)		w l ls 4 title 'Вязкость водяного пара'
set output 'images/viscosity_gas.png'
plot [200:700]	0.00001781 * (300.55 + 111) / (x + 111) * ((x / 300.55) ** (3. / 2.))	w l ls 1 title 'Вязкость азота', \
		 		0.00002018 * (292.25 + 127) / (x + 127) * ((x / 292.25) ** (3. / 2.))	w l ls 2 title 'Вязкость кислорода', \
				0.0000148 * (293.15  + 240) / (x + 240) * ((x / 293.15) ** (3. / 2.))	w l ls 3 title 'Вязкость углекислого газа', \
				0.000018 * 330 / x * ((x / 330) ** (3. / 2.))			w l ls 4 title 'Вязкость водяного пара'
