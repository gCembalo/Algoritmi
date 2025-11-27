reset

# plots Poisson

# Settiamo le dimensioni e caratteristiche del canva
set term qt 1 size 600,400
set title "Poisson"
set xlabel "r"  font ",14"
set ylabel "phi(r)"  font ",14"

# Settiamo la griglia
set grid

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

# Plottiamo effettivamente i file di dati
plot "poisson.dat" using 1:2 index 0 w lp ps 0.1 lc rgb "red" title "s=0.0"
replot "poisson.dat" using 1:2 index 1 w lp ps 0.1 lc rgb "green" title "s=0.2"
replot "poisson.dat" using 1:2 index 2 w lp ps 0.1 lc rgb "blue" title "s=0.4"
replot "poisson.dat" using 1:2 index 3 w lp ps 0.1 lc rgb "pink" title "s=0.6"
replot "poisson.dat" using 1:2 index 4 w lp ps 0.1 lc rgb "black" title "s=0.8"
replot "poisson.dat" using 1:2 index 5 w lp ps 0.1 lc rgb "brown" title "s=1.0"





# plots Residual

# Settiamo le dimensioni e caratteristiche del canva
set term qt 2 size 600,400
set title "Residual"
set xlabel "s"  font ",14"
set ylabel "Residual(s)"  font ",14"

# Settiamo la griglia
set grid

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

# Plottiamo effettivamente i file di dati
plot "poisson.dat" using 1:2 index 6 w lp ps 0.1 lc rgb "blue" title "Residual"





# plots Confronto

# Settiamo le dimensioni e caratteristiche del canva
set term qt 3 size 600,400
set title "Confronto soluzione analitica"
set xlabel "r"  font ",14"
set ylabel "\phi(r)/r"  font ",14"

# Settiamo la griglia
set grid

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

# Plottiamo effettivamente i file di dati
plot "poisson.dat" using 1:2 index 7 w lp ps 0.1 lc rgb "blue" title "RK4"
replot "poisson.dat" using 1:3 index 7 w lp ps 0.1 lc rgb "red" title "Solution"