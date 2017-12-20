#!/usr/bin/gnuplot
reset
set term pngcairo size 	800, 600 enhanced dashed font 'Helvetica,20'
#set term png size 550, 450
set style line 1 lc 1 lt 1 lw 4
set title "Ускорение работы программы при размере \nрасчётной области 30х30 ячеек"
set xlabel "Количество процессов"
set ylabel "Ускорение"
set output 'images/acceleration30.png'
plot 'plotdat/acceleration30.dat' with linespoints ls 1 ps 2 pt 7 notitle
set title "Ускорение работы программы при размере \nрасчётной области 60х60 ячеек"
set output 'images/acceleration60.png'
plot 'plotdat/acceleration60.dat' with linespoints ls 1 ps 2 pt 7 notitle
