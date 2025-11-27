reset

# plots Kepler Orbit
set term qt 1 size 600,400
set title "Kepler"
set xlabel "x"  font ",14"
set ylabel "y"  font ",14"

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

# disegno il sole
set object circle at first 0,0 radius char 0.5 \
    fillstyle empty border lc rgb "orange" lw 7

plot "kepler.dat" using 1:2 w lp lc rgb "blue" title "Kepler"