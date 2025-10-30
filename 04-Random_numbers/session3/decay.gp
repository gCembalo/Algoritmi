reset

# decay plots
set term qt size 600,400
set xlabel "t"  font ",14"
set ylabel "N(t)/N0"  font ",14"
set log y
set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

plot "decay.dat" pt 7 ps 0.8 lc "red" title "Data", \
    exp(-0.01*x) title "exp(-0.01*t)"