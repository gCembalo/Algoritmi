reset

# gaussian plots
set term qt size 600,400
set xlabel "x"  font ",14"
set ylabel "y"  font ",14"
set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

plot "gaussian.dat" pt 7 ps 0.4 lc "red" title "Data", \
    1./(0.5*sqrt(2.*3.1415926)) * exp( -0.5 * (x*x) / (0.5*0.5) ) title "f(x)"