reset

# Set column indices
ix = 2
iy = 4
ixB = 6
iR = 8 #(RK)

# Set plot specifications
set grid
set key font ",14"
set title font ",14"
set tics font ",14"
set xlabel "t" font ",14"
set ylabel "x,y,R(t) [RK]" font ",14"

# Plot
plot "coupled_springs.dat" using 1:ix title "x" w lines
replot "coupled_springs.dat" using 1:iy title "y" w lines
replot "coupled_springs.dat" using 1:iR title "R" w lines