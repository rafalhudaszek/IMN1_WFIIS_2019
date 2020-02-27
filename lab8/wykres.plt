#!/usr/bin/gnuplot
set term png


set xl "x"
set yl "y"
set view map


set out 'mapa_0.png'
splot 'mapa_0.txt' u 1:2:3 w pm3d t ""

set out 'mapa_1.png'
splot 'mapa_1.txt' u 1:2:3 w pm3d t ""

set out 'mapa_2.png'
splot 'mapa_2.txt' u 1:2:3 w pm3d t ""

set out 'mapa_3.png'
splot 'mapa_3.txt' u 1:2:3 w pm3d t ""

set out 'mapa_4.png'
splot 'mapa_4.txt' u 1:2:3 w pm3d t ""




