#!/usr/bin/gnuplot
reset
set term pngcairo size 	800, 600 enhanced dashed font 'Helvetica,20'
set title "Капиллярное давление по модели Кори"
set xlabel "Насыщенность водой"
set ylabel "Давление, Па"
set grid x y
set linestyle  1 lc 1 lt 1 dt 1 lw 4
set xrange [0:1]
set output "images/capillary_pressure_water.png"
plot	((x > 0.15) && (x <= 0.8)) ? (55000 * ((x / 0.7) ** (-0.5))) : (1 / 0)	w l ls 1 title 'P_{cow}(s_w)'
set xlabel "Насыщенность газом"
set output "images/capillary_pressure_gas.png"
plot	(x <= 0.65) ? ((x / 0.7) ** (-0.5)) : (1 / 0)							w l ls 1 title 'P_{cog}(s_g)'
set key bottom right
set title "Производная капиллярного давления по насыщенности       "
set xlabel "Насыщенность водой"
set output "images/capillary_pressure_water_derivative.png"
plot	((x > 0.15) && (x <= 0.8)) ? (-55000 * (0.7 ** 0.5) * 0.5 * (x ** (-1.5))) : (1 / 0)	w l ls 1 title 'dP_{cow}(s_w)/ds_w'
set xlabel "Насыщенность газом"
set output "images/capillary_pressure_gas_derivative.png"
plot	(x <= 0.65) ? (-(0.7 ** 0.5) * 0.5 * (x ** (-1.5))) : (1 / 0)							w l ls 1 title 'dP_{cog}(s_g)/ds_g'
