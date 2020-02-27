#!/usr/bin/gnuplot
set term png





set xl "x"
set yl "y"
set view map


set out 'map_x_0.0.png'
splot 'map_x1_0.0.txt' u 1:2:3 w pm3d t ""

set out 'mapa_vy_0.0.png'
splot 'map_y1_0.0.txt' u 1:2:3 w pm3d t ""

set out 'mapa_vx_0.1.png'
splot 'map_x1_0.1.txt' u 1:2:3 w pm3d t ""

set out 'mapa_vy_0.1.png'
splot 'map_y1_0.1.txt' u 1:2:3 w pm3d t ""







