reset

# plots Ode2
set term qt 1 size 600,400
set title "Ode2"
set xlabel "t"  font ",14"
set ylabel "x(t) ; y(t)"  font ",14"

set xrange[-2:2]
set yrange[-2:2]

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

plot "ode2.dat" using 2:3 index 0 w lp title "Euler"
replot "ode2.dat" using 2:3 index 1 w lp title "RK2_{Mid}"
replot "ode2.dat" using 2:3 index 2 w lp title "RK2_{Heun}"
replot "ode2.dat" using 2:3 index 3 w lp title "RK4"