#!/usr/bin/gnuplot
set term png





set xl "t"
set yl "x"
set view map


set out 'mapa_0.png'
splot 'mapa_0.txt' u 1:2:3 w pm3d t ""

set out 'mapa_0.1.png'
splot 'mapa_0.1.txt' u 1:2:3 w pm3d t ""

set out 'mapa_1.png'
splot 'mapa_1.txt' u 1:2:3 w pm3d t ""

set out 'mapa_1alfa.png'
splot 'mapa_1alfa.txt' u 1:2:3 w pm3d t ""






