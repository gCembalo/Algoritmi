reset

# plots derivative FD , BD , CD , 4th
set term qt size 600,400
set xlabel "1/h"  font ",14"
set logscale x
set logscale y
set ylabel "err"  font ",14"
set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

plot \
    "FD.dat" using 1:2 with points pt 7 ps 0.4 lc rgb "red"  title "FD", \
    "BD.dat" using 1:2 with points pt 7 ps 0.4 lc rgb "blue" title "BD", \
    "CD.dat" using 1:2 with points pt 7 ps 0.4 lc rgb "green" title "CD", \
    "4thorder.dat" using 1:2 with points pt 7 ps 0.4 lc rgb "black" title "4th"