#!/usr/bin/gnuplot
set term png

set output "mapa16.png"
set xlabel "x"
set ylabel "y"
set title "k=16"
reset
set pm3d map
set palette defined (-1 "red", 0 "white", 1 "blue")
set size ratio -1

splot [0:30][0:30] "map.dat" i 0 u 1:2:3

reset
set output "mapa8.png"
set xlabel "x"
set ylabel "y"
set title "k=8"

set pm3d map
set palette defined (-1 "red", 0 "white", 1 "blue")
set size ratio -1

splot [0:30][0:30] "map.dat" i 1 u 1:2:3


reset
set output "mapa4.png"
set xlabel "x"
set ylabel "y"
set title "k=4"

set pm3d map
set palette defined (-1 "red", 0 "white", 1 "blue")
set size ratio -1

splot [0:30][0:30] "map.dat" i 2 u 1:2:3


reset
set output "mapa2.png"
set xlabel "x"
set ylabel "y"
set title "k=2"

set pm3d map
set palette defined (-1 "red", 0 "white", 1 "blue")
set size ratio -1

splot [0:30][0:30] "map.dat" i 3 u 1:2:3


reset
set output "mapa1.png"
set xlabel "x"
set ylabel "y"
set title "k=1"

set pm3d map
set palette defined (-1 "red", 0 "white", 1 "blue")
set size ratio -1

splot [0:30][0:30] "map.dat" i 4 u 1:2:3
