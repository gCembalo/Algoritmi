reset

# plots EXAMPLE

# Settiamo le dimensioni e caratteristiche del canva
set term qt 1 size 600,400
set title "EXAMPLE"
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

# Possiamo aggiungere disegni a piacere (qua ad esempio disegnamo un cerchio
# pieno in 0,0)
set object circle at first 0,0 radius char 0.5 \
    fillstyle empty border lc rgb "orange" lw 7

# Plottiamo effettivamente i file di dati
# Sintassi: "plot" (o replot) è il comando per disegnare, 
# "using 1:4 index 1" gli dice di usare le colonne 1 e 4 dell'indice 1
# "w lp" gli dice di disegnare with line-and-points
# "ps 0.1" indica lo spessore dei punti/linee
# "lc rgb "red" " gli indica il colore
# "title "titolo" " da un nome al particolare set di dati
# Per andare a capo senza interrompere il comando usa ", \" (vedi riga 41-42)
plot "harmonic.dat" using 1:4 index 1 w lp ps 0.1 lc rgb "red" title "RK2_{Mid}"
# replot "harmonic.dat" using 1:4 index 0 w lp ps 0.1 lc rgb "blue" title "Euler"
replot "harmonic.dat" using 1:4 index 2 w lp ps 0.1 lc rgb "black" title "RK4"
replot "harmonic.dat" using 1:4 index 3 w lp ps 0.1 lc rgb "orange", \
        title "VerletPosition"
replot "harmonic.dat" using 1:4 index 4 w lp ps 0.1 lc rgb "purple" title "VerletVelocity"