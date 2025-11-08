reset

# plots Ode1
set term qt 1 size 600,400
set title "Ode"
set xlabel "t"  font ",14"
set ylabel "y(t)"  font ",14"

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

# plot "ode1Ex.dat" using 1:2 index 0 with lines title "ode1Sol(t)"
plot "ode1.dat" using 1:2 index 0 with lines title "h=0.5"
replot "ode1.dat" using 1:2 index 1 with lines title "h=0.2"
replot "ode1.dat" using 1:2 index 2 with lines title "h=0.1"
replot "ode1.dat" using 1:2 index 3 with lines title "h=0.05"
replot "ode1.dat" using 1:2 index 4 with lines title "h=0.02"
replot "ode1.dat" using 1:2 index 5 with lines title "h=0.01"
replot "ode1.dat" using 1:2 index 6 with lines title "h=0.005"
replot "ode1.dat" using 1:2 index 7 with lines title "h=0.002"
replot "ode1.dat" using 1:2 index 8 with lines title "h=0.001"


# plots Ode1 errore
set term qt 2 size 600,400
set title "Ode absolute error"
set logscale y
set xlabel "t"  font ",14"
set ylabel "err_{abs}(t)"  font ",14"

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

plot "ode1.dat" using 1:3 index 0 with lines title "h=0.5"
replot "ode1.dat" using 1:3 index 1 with lines title "h=0.2"
replot "ode1.dat" using 1:3 index 2 with lines title "h=0.1"
replot "ode1.dat" using 1:3 index 3 with lines title "h=0.05"
replot "ode1.dat" using 1:3 index 4 with lines title "h=0.02"
replot "ode1.dat" using 1:3 index 5 with lines title "h=0.01"
replot "ode1.dat" using 1:3 index 6 with lines title "h=0.005"
replot "ode1.dat" using 1:3 index 7 with lines title "h=0.002"
replot "ode1.dat" using 1:3 index 8 with lines title "h=0.001"


# plots Ode1 errore
set term qt 3 size 600,400
set title "Ode relative error"
set logscale y
set xlabel "t"  font ",14"
set ylabel "err_{rel}(t)"  font ",14"

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

plot "ode1.dat" using 1:4 index 0 with lines title "h=0.5"
replot "ode1.dat" using 1:4 index 1 with lines title "h=0.2"
replot "ode1.dat" using 1:4 index 2 with lines title "h=0.1"
replot "ode1.dat" using 1:4 index 3 with lines title "h=0.05"
replot "ode1.dat" using 1:4 index 4 with lines title "h=0.02"
replot "ode1.dat" using 1:4 index 5 with lines title "h=0.01"
replot "ode1.dat" using 1:4 index 6 with lines title "h=0.005"
replot "ode1.dat" using 1:4 index 7 with lines title "h=0.002"
replot "ode1.dat" using 1:4 index 8 with lines title "h=0.001"