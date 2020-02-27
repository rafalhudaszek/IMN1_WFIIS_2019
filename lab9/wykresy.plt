#!/usr/bin/gnuplot
set term png





set xl "t"
set yl "E"
set view map


set out 'wykres_0.png'
plot 'wykres_0.txt' u 1:2 w l t ""

set out 'wykres_0.1.png'
plot 'wykres_0.1.txt' u 1:2 w l t ""

set out 'wykres_1.png'
plot 'wykres_1.txt' u 1:2 w l t ""

set out 'wykres_1alfa.png'
plot 'wykres_1alfa.txt' u 1:2 w l t ""






