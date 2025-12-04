reset

# plots bvp

# Settiamo le dimensioni e caratteristiche del canva
set term qt 1 size 600,400
set title "bvp"
set xlabel "x"  font ",14"
set ylabel "y(x)"  font ",14"

# Settiamo il range dei due assi
# set xrange[-2:2]
# set yrange[-2:2]

# Settiamo le scale logaritmiche
# set logscale y
# set logscale x

# Settiamo la griglia
set grid

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

# Plottiamo effettivamente i file di dati
plot "bvp.dat" using 1:2 index 0 w lp ps 0.1 lc rgb "red" title "bvp"