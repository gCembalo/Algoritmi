reset

# First plot
set term qt size 600,400
set xlabel "i"  font ",14"
set ylabel "x"  font ",14"
set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

plot "uniformity.dat" pt 7 ps 1.5 lc "red"

######### 

# Second plot
set term qt 2 size 600,400
set log x
set log y
set tics font ",14"
set xlabel "N"  font ",14"
set ylabel "err" font ",14"
set key          font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

plot "prn_uniformity.dat" using 1:3 pt 7 ps 1.5 title "Data", \
     1.0/sqrt(x) title "1/sqrt(N)"