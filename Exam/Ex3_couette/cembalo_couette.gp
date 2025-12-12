reset

# plots Couette

# Settiamo le dimensioni e caratteristiche del canva
set term qt 1 size 600,400
set title "Couette"
set xlabel "r"  font ",14"
set ylabel "u"  font ",14"

# Settiamo il range dei due assi
# set xrange[-2:2]
# set yrange[-2:2]

# Settiamo la griglia
set grid

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

# Plottiamo effettivamente i file di dati
plot "couette.dat" using 1:2 index 1 w lp ps 0.1 lc rgb "red" title "RK4"
replot "couette.dat" using 1:2 index 2 w lp ps 0.1 lc rgb "black" title "FD"





reset

# plots Residual

# Settiamo le dimensioni e caratteristiche del canva
set term qt 2 size 600,400
set title "Couette Residual"
set xlabel "s"  font ",14"
set ylabel "R"  font ",14"

# Settiamo il range dei due assi
#set xrange[-10:10]
#set yrange[-15:15]

# Settiamo la griglia
set grid

set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

# Plottiamo effettivamente i file di dati
plot "couette.dat" using 1:2 index 0 w lp ps 0.1 lc rgb "blue" title "Residual"