reset

# plots Ode2 Convergence
set term qt 1 size 600,400
set title "Ode2 Convergence"
set xlabel "dt"  font ",14"
set ylabel "err"  font ",14"

set logscale y
set logscale x

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

plot "ode2Conv.dat" using 1:2 index 0 w lp title "Euler"
replot "ode2Conv.dat" using 1:3 index 0 w lp title "RK2_{Mid}"
replot "ode2Conv.dat" using 1:4 index 0 w lp title "RK2_{Heun}"
replot "ode2Conv.dat" using 1:5 index 0 w lp title "RK4"