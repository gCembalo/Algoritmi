reset

fname = "coupled_rotation.dat"

# Set plot speicifications
set grid
set key left font ",14"
set title font ",14"
set tics font ",14"
set xlabel "x" font ",14"
set ylabel "theta" font ",14"

# Plot
set xrange[0:10]
#set yrange[-0.1:2.5]
plot fname using 1:2 title "th1 (RK4)" w lines
replot fname using 1:3 title "th2 (RK4)" w lines
replot fname using 1:4 title "th2 (PV)" w lines
replot fname using 1:5 title "th2 (PV)" w lines






# complete plot
fname = "coupled_rotation.dat"

# Set plot speicifications
set term qt 2 size 600,400
set grid
set key left font ",14"
set title font ",14"
set tics font ",14"
set xlabel "x" font ",14"
set ylabel "theta" font ",14"

# Plot
set xrange[0:50]
#set yrange[-0.1:2.5]
plot fname using 1:2 index 1 title "th1 (RK4)" w lines
replot fname using 1:3 index 1 title "th2 (RK4)" w lines
replot fname using 1:4 index 1 title "th2 (PV)" w lines
replot fname using 1:5 index 1 title "th2 (PV)" w lines