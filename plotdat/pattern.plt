#!/usr/bin/gnuplot
reset
set term png size 1920, 1920
#set title "Pattern of the matrix" font "Helvetica,40"
set title "Шаблон матрицы" font "Helvetica,40"
set nokey
set palette defined (0 "white", 1 "black")
set cbrange [0:1]
#set tic scale 0
unset cbtics
unset tics
set size square
set xrange[-0.5:250.5]
set yrange[-0.5:250.5]
set output 'images/pattern.png'
set cbrange [0:1]
set view map
plot	'tmp/A_pattern_matrix.dat' matrix w image notitle
