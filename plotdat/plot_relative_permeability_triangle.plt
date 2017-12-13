#!/usr/bin/gnuplot
reset
#set term pngcairo size 550, 500
set term pngcairo size 	800, 600 enhanced dashed font 'Helvetica,20'
set title 'Относительная проницаемость нефти по модели Бэкера'
set output 'images/relative_permeability_oil.png'
unset border
unset xtics
unset ytics
set size square
set xrange [0:1]
set yrange [0:1]

# Отрисовываем основной треугольник
set style arrow 11 nohead back lt -1 lw 2 # стиль

set arrow 11 from 0,0 to .5,sqrt(3)/2 as 11
set arrow 12 from 0,0 to 1,0          as 11
set arrow 13 from .5,sqrt(3)/2 to 1,0 as 11

# Промежуточные деления
set style arrow 2 nohead back lt 1 lw 1 lc rgb "gray" # стиль

set arrow 21 from 0.1,0 to 0.5+0.1*0.5,0.9*sqrt(3)/2 as 2
set arrow 22 from 0.2,0 to 0.5+0.2*0.5,0.8*sqrt(3)/2 as 2
set arrow 23 from 0.3,0 to 0.5+0.3*0.5,0.7*sqrt(3)/2 as 2
set arrow 24 from 0.4,0 to 0.5+0.4*0.5,0.6*sqrt(3)/2 as 2
set arrow 25 from 0.5,0 to 0.5+0.5*0.5,0.5*sqrt(3)/2 as 2
set arrow 26 from 0.6,0 to 0.5+0.6*0.5,0.4*sqrt(3)/2 as 2
set arrow 27 from 0.7,0 to 0.5+0.7*0.5,0.3*sqrt(3)/2 as 2
set arrow 28 from 0.8,0 to 0.5+0.8*0.5,0.2*sqrt(3)/2 as 2
set arrow 29 from 0.9,0 to 0.5+0.9*0.5,0.1*sqrt(3)/2 as 2

set arrow 31 from 0.1*0.5,0.1*sqrt(3)/2 to 0.1,0 as 2
set arrow 32 from 0.2*0.5,0.2*sqrt(3)/2 to 0.2,0 as 2
set arrow 33 from 0.3*0.5,0.3*sqrt(3)/2 to 0.3,0 as 2
set arrow 34 from 0.4*0.5,0.4*sqrt(3)/2 to 0.4,0 as 2
set arrow 35 from 0.5*0.5,0.5*sqrt(3)/2 to 0.5,0 as 2
set arrow 36 from 0.6*0.5,0.6*sqrt(3)/2 to 0.6,0 as 2
set arrow 37 from 0.7*0.5,0.7*sqrt(3)/2 to 0.7,0 as 2
set arrow 38 from 0.8*0.5,0.8*sqrt(3)/2 to 0.8,0 as 2
set arrow 39 from 0.9*0.5,0.9*sqrt(3)/2 to 0.9,0 as 2

set arrow 41 from 0.1*0.5,0.1*sqrt(3)/2 to 0.5+0.9*0.5,0.1*sqrt(3)/2 as 2
set arrow 42 from 0.2*0.5,0.2*sqrt(3)/2 to 0.5+0.8*0.5,0.2*sqrt(3)/2 as 2
set arrow 43 from 0.3*0.5,0.3*sqrt(3)/2 to 0.5+0.7*0.5,0.3*sqrt(3)/2 as 2
set arrow 44 from 0.4*0.5,0.4*sqrt(3)/2 to 0.5+0.6*0.5,0.4*sqrt(3)/2 as 2
set arrow 45 from 0.5*0.5,0.5*sqrt(3)/2 to 0.5+0.5*0.5,0.5*sqrt(3)/2 as 2
set arrow 46 from 0.6*0.5,0.6*sqrt(3)/2 to 0.5+0.4*0.5,0.6*sqrt(3)/2 as 2
set arrow 47 from 0.7*0.5,0.7*sqrt(3)/2 to 0.5+0.3*0.5,0.7*sqrt(3)/2 as 2
set arrow 48 from 0.8*0.5,0.8*sqrt(3)/2 to 0.5+0.2*0.5,0.8*sqrt(3)/2 as 2
set arrow 49 from 0.9*0.5,0.9*sqrt(3)/2 to 0.5+0.1*0.5,0.9*sqrt(3)/2 as 2

# Основные подписи (что по углам)
set label 11 center "Ca_2Si_2O_6" at 0.5,sqrt(3)/2+.06
set label 12 center "Fe_2Si_2O_6" at 1,-.06
set label 13 center "Mg_2Si_2O_6" at 0,-.06

# Подписи шкалы
set label 21 center "50" at 0.5,-.05
set label 32 right "50" at .5*.5-0.04,0.5*sqrt(3)/2
set label 43 left "50" at 1-.5*.5+0.04,.50*sqrt(3)/2

# Собственно построение с пересчетом
plot "plotdat/no.dat" using (sqrt(3)/2*(1-$1/($1+$2+$3))/sqrt(3)*2-(sqrt(3)/2*$3/($1+$2+$3))/sqrt(3)):(sqrt(3)/2*$3/($1+$2+$3)) with points pt 13 ps 1.5 lc rgb "#0000ff"
