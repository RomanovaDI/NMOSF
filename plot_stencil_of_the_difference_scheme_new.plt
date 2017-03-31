#!/usr/bin/gnuplot
reset
set term png size 1920, 1080
set title "Stencil of the difference scheme"
unset hidden3d
#set xrange[-2:1]
#set yrange[-2:1]
#set zrange[-2:1]
#unset border
set border 4095
unset tics
set output 'image.png'
set key off
#set xtics ("i+1" 1, "i" 0, "i-1" -1, "i-2" -2)
#set ytics ("j+1" 1, "j" 0, "j-1" -1, "j-2" -2)
#set ztics ("k+1" 1, "k" 0, "k-1" -1, "k-2" -2)
unset box
#set box
set view 25, 130
#set view 90, 0
#splot	'test_points.dat' using 1:2:3 with lines ls 0 lc rgb "gray",\
splot	'grid_lines.dat' using 1:2:3 with lines ls 0 lc rgb "gray",\
		'lines.dat' using 1:2:3 with lines ls 1 lc rgb "black",\
		'momentum_eqn0.dat' using ($1==0?$2:1/0):($1==0?$3:1/0):($1==0?$4:1/0) with point pointtype 1 pointsize 5 lc rgb "blue",\
		'momentum_eqn0.dat' using ($1==1?$2:1/0):($1==1?$3:1/0):($1==1?$4:1/0) with point pointtype 2 pointsize 5 lc rgb "blue",\
		'momentum_eqn0.dat' using ($1==2?$2:1/0):($1==2?$3:1/0):($1==2?$4:1/0) with point pointtype 4 pointsize 5 lc rgb "blue",\
		'momentum_eqn0.dat' using ($1==3?$2:1/0):($1==3?$3:1/0):($1==3?$4:1/0) with point pointtype 6 pointsize 5 lc rgb "blue",\
		'momentum_eqn0.dat' using ($1==4?$2:1/0):($1==4?$3:1/0):($1==4?$4:1/0) with point pointtype 8 pointsize 5 lc rgb "blue",\
		'momentum_eqn1.dat' using ($1==0?$2:1/0):($1==0?$3:1/0):($1==0?$4:1/0) with point pointtype 1 pointsize 4 lc rgb "red",\
		'momentum_eqn1.dat' using ($1==1?$2:1/0):($1==1?$3:1/0):($1==1?$4:1/0) with point pointtype 2 pointsize 4 lc rgb "red",\
		'momentum_eqn1.dat' using ($1==2?$2:1/0):($1==2?$3:1/0):($1==2?$4:1/0) with point pointtype 4 pointsize 4 lc rgb "red",\
		'momentum_eqn1.dat' using ($1==3?$2:1/0):($1==3?$3:1/0):($1==3?$4:1/0) with point pointtype 6 pointsize 4 lc rgb "red",\
		'momentum_eqn1.dat' using ($1==4?$2:1/0):($1==4?$3:1/0):($1==4?$4:1/0) with point pointtype 8 pointsize 4 lc rgb "red",\
		'momentum_eqn2.dat' using ($1==0?$2:1/0):($1==0?$3:1/0):($1==0?$4:1/0) with point pointtype 1 pointsize 3 lc rgb "green",\
		'momentum_eqn2.dat' using ($1==1?$2:1/0):($1==1?$3:1/0):($1==1?$4:1/0) with point pointtype 2 pointsize 3 lc rgb "green",\
		'momentum_eqn2.dat' using ($1==2?$2:1/0):($1==2?$3:1/0):($1==2?$4:1/0) with point pointtype 4 pointsize 3 lc rgb "green",\
		'momentum_eqn2.dat' using ($1==3?$2:1/0):($1==3?$3:1/0):($1==3?$4:1/0) with point pointtype 6 pointsize 3 lc rgb "green",\
		'momentum_eqn2.dat' using ($1==4?$2:1/0):($1==4?$3:1/0):($1==4?$4:1/0) with point pointtype 8 pointsize 3 lc rgb "green",\
		'poissons_eqn.dat' using ($1==0?$2:1/0):($1==0?$3:1/0):($1==0?$4:1/0) with point pointtype 1 pointsize 2 lc rgb "cyan",\
		'poissons_eqn.dat' using ($1==1?$2:1/0):($1==1?$3:1/0):($1==1?$4:1/0) with point pointtype 2 pointsize 2 lc rgb "cyan",\
		'poissons_eqn.dat' using ($1==2?$2:1/0):($1==2?$3:1/0):($1==2?$4:1/0) with point pointtype 4 pointsize 2 lc rgb "cyan",\
		'poissons_eqn.dat' using ($1==3?$2:1/0):($1==3?$3:1/0):($1==3?$4:1/0) with point pointtype 6 pointsize 2 lc rgb "cyan",\
		'poissons_eqn.dat' using ($1==4?$2:1/0):($1==4?$3:1/0):($1==4?$4:1/0) with point pointtype 8 pointsize 2 lc rgb "cyan",\
		'snow_volume_fraction_eqn.dat' using ($1==0?$2:1/0):($1==0?$3:1/0):($1==0?$4:1/0) with point pointtype 1 pointsize 1 lc rgb "magenta",\
		'snow_volume_fraction_eqn.dat' using ($1==1?$2:1/0):($1==1?$3:1/0):($1==1?$4:1/0) with point pointtype 2 pointsize 1 lc rgb "magenta",\
		'snow_volume_fraction_eqn.dat' using ($1==2?$2:1/0):($1==2?$3:1/0):($1==2?$4:1/0) with point pointtype 4 pointsize 1 lc rgb "magenta",\
		'snow_volume_fraction_eqn.dat' using ($1==3?$2:1/0):($1==3?$3:1/0):($1==3?$4:1/0) with point pointtype 6 pointsize 1 lc rgb "magenta",\
		'snow_volume_fraction_eqn.dat' using ($1==4?$2:1/0):($1==4?$3:1/0):($1==4?$4:1/0) with point pointtype 8 pointsize 1 lc rgb "magenta",\
		'name_points.dat' using 1:2:3:4 with labels point pointtype 0 offset character 2,character 1 tc rgb "brown",\
		'test_points.dat' using 1:2:3 with lines ls 1 lc rgb "black",\
		'test_points.dat' using 1:2:3:4 with labels point pointtype 1 offset character 1,character 1 tc rgb "black" lc rgb "black"
