#!/usr/local/Cellar/gnuplot/5.2.7_1/bin/gnuplot -persist
#terminal
set terminal postscript eps enhanced color font 'Helvetica, 6' size 210 cm, 148 cm
set output 'result.eps'

set grid
set logscale xy
set xlabel 'T'
set ylabel 'U'
set zlabel 'm'

splot 'data.dat' u (1./$1):2:3
#    EOF
