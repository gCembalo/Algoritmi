reset

# plots Wave

# Settiamo le dimensioni e caratteristiche del canva
set term qt 1 size 600,400
set title "Wave"
set xlabel "x"  font ",14"
set ylabel "phi(x)"  font ",14"

# Settiamo la griglia
set grid

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

# Plottiamo effettivamente i file di dati
plot "wave.dat" using 1:2 index 0 w lp ps 0.1 lc rgb "red" title "k=1"
replot "wave.dat" using 1:2 index 1 w lp ps 0.1 lc rgb "green" title "k=2"
replot "wave.dat" using 1:2 index 2 w lp ps 0.1 lc rgb "blue" title "k=3"
replot "wave.dat" using 1:2 index 3 w lp ps 0.1 lc rgb "pink" title "k=4"
replot "wave.dat" using 1:2 index 4 w lp ps 0.1 lc rgb "black" title "k=5"