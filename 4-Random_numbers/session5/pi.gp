reset

# pi plots
#set term qt size 400,400
#set xlabel "x"  font ",14"
#set ylabel "y"  font ",14"
#set key         font ",14"
#set lmargin at screen 0.12
#set rmargin at screen 0.95
#set bmargin at screen 0.12
#set tmargin at screen 0.95

#plot "pi.dat" pt 7 ps 0.4 lc "red" title "Data"
#replot x*x + y*y = 1
#replot "out_pi.dat" pt 7 ps 0.4

# err_pi plots
set term qt 2 size 600,400
set log y
set log x
set xlabel "x"  font ",14"
set ylabel "y"  font ",14"
set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

plot "err_pi.dat" pt 7 ps 0.4 lc "red" title "Error Data", \
    1/sqrt(x) title "1/sqrt(N)"