reset

# plots harmonic
set term qt 1 size 600,400
set title "harmonic"
set xlabel "t"  font ",14"
set ylabel "y(t)"  font ",14"

# set xrange[-2:2]
# set yrange[-2:2]

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

plot "harmonic.dat" using 1:4 index 1 w lp lc rgb "red" title "RK2_{Mid}"
# replot "harmonic.dat" using 1:4 index 0 w lp lc rgb "blue" title "Euler"
replot "harmonic.dat" using 1:4 index 2 w lp lc rgb "black" title "RK4"
replot "harmonic.dat" using 1:4 index 3 w lp lc rgb "orange" title "VerletPosition"
replot "harmonic.dat" using 1:4 index 4 w lp lc rgb "purple" title "VerletVelocity"