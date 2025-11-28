reset

# plots qho

# Settiamo le dimensioni e caratteristiche del canva
set term qt 1 size 600,400
set title "Qho"
set xlabel "x"  font ",14"
set ylabel "E(x)"  font ",14"

# Settiamo la scala logaritmica
set logscale y

# Settiamo la griglia
set grid

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

# Plottiamo effettivamente i file di dati
plot "qho.dat" using 1:2 index 0 w lp ps 0.1 lc rgb "red" title "RK4 avanti"
replot "qho.dat" using 1:2 index 1 w lp ps 0.1 lc rgb "blue" title "RK4 dietro"



reset

# plots qho 2

# Settiamo le dimensioni e caratteristiche del canva
set term qt 2 size 600,400
set title "Residual"
set xlabel "E"  font ",14"
set ylabel "Residual(E)"  font ",14"

# Settiamo la scala logaritmica
#set logscale y

# Settiamo il range dei due assi
set xrange[0:5]

# Settiamo la griglia
set grid

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

# Plottiamo effettivamente i file di dati
plot "qho.dat" using 1:2 index 2 w lp ps 0.1 lc rgb "red" title "Residual"