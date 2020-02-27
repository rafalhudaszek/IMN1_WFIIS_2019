#!/usr/bin/gnuplot

set term png
set grid

set output "trap_x.png"
set xlabel "t"
set ylabel "x(t)"
set title "Metoda trapezów x"

plot "trapezy_-2.dat"  u 1:3 w l lw 1 t "TOL 10^{-2}", "trapezy_-5.dat"  u 1:3 w l lw 1 t "TOL 10^{-5}"


set output "trap_v.png"
set xlabel "t"
set ylabel "v(t)"
set title "Metoda trapezów v"

plot "trapezy_-2.dat" u 1:4 w l lw 1 t "TOL 10^{-2}", "trapezy_-5.dat"  u 1:4 w l lw 1 t "TOL 10^{-5}"


set output "trap_dt.png"
set xlabel "t"
set ylabel "dt(t)"
set title "Metoda trapezów dt(t)"

plot "trapezy_-2.dat"  u 1:2 w l lw 1 t "TOL 10^{-2}", "trapezy_-5.dat"  u 1:2 w l lw 1 t "TOL 10^{-5}"


set output "trap_v(x).png"
set xlabel "x"
set ylabel "v(x)"
set title "Metoda trapezów v(x)"

plot "trapezy_-2.dat" u 3:4 w l lw 1 t "TOL 10^{-2}", "trapezy_-5.dat"  u 3:4 w l lw 1 t "TOL 10^{-5}"



set output "Rk2_x.png"
set xlabel "t"
set ylabel "x(t)"
set title "Metoda RK2 x"

plot "RK2_pow-2.dat"  u 1:3 w l lw 1 t "TOL 10^{-2}", "RK2_pow-5.dat" u 1:3 w l lw 1 t "TOL 10^{-5}"


set output "Rk2_v.png"
set xlabel "t"
set ylabel "v(t)"
set title "Metoda RK2 v"

plot "RK2_pow-2.dat"  u 1:4 w l lw 1 t "TOL 10^{-2}", "RK2_pow-5.dat"  u 1:4 w l lw 1 t "TOL 10^{-5}"


set output "RK2_dt.png"
set xlabel "t"
set ylabel "dt(t)"
set title "Metoda Rk2 dt(t)"

plot "RK2_pow-2.dat"  u 1:2 w l lw 1 t "TOL 10^{-2}", "RK2_pow-5.dat" u 1:2 w l lw 1 t "TOL 10^{-5}"


set output "Rk2_v(x).png"
set xlabel "x"
set ylabel "v(x)"
set title "Metoda Rk2 v(x)"

plot "RK2_pow-2.dat" u 3:4 w l lw 1 t "TOL 10^{-2}", "RK2_pow-5.dat"  u 3:4 w l lw 1 t "TOL 10^{-5}"


