reset

# plots derivative position , velocity , acceleration with FD
set term qt 1 size 600,400
set title "Derivata FD"
set xlabel "t"  font ",14"
set ylabel "x , v , a"  font ",14"
set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

plot \
    "trajectory.dat" using 1:2 with points pt 7 ps 0.1 lc rgb "red"  title "x(t)", \
    "velocity2.dat" using 1:2 with points pt 7 ps 0.1 lc rgb "green" title "v(t)", \
    "acceleration.dat" using 1:2 with points pt 7 ps 0.1 lc rgb "blue" title "a(t)"


# plots derivative position , velocity , acceleration with BD
set term qt 2 size 600,400
set title "Derivata BD"
set xlabel "t"  font ",14"
set ylabel "x , v , a"  font ",14"
set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

plot \
    "trajectory.dat" using 1:2 with points pt 7 ps 0.1 lc rgb "red"  title "x(t)", \
    "velocity2.dat" using 1:2 with points pt 7 ps 0.1 lc rgb "green" title "v(t)", \
    "acceleration.dat" using 1:2 with points pt 7 ps 0.1 lc rgb "blue" title "a(t)"


# plots derivative position , velocity , acceleration with CD
set term qt 3 size 600,400
set title "Derivata CD"
set xlabel "t"  font ",14"
set ylabel "x , v , a"  font ",14"
set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

plot \
    "trajectory.dat" using 1:2 with points pt 7 ps 0.1 lc rgb "red"  title "x(t)", \
    "velocity3.dat" using 1:2 with points pt 7 ps 0.1 lc rgb "green" title "v(t)", \
    "acceleration.dat" using 1:2 with points pt 7 ps 0.1 lc rgb "blue" title "a(t)"


# plots derivative position , velocity , acceleration with 4thDer
set term qt 4 size 600,400
set title "Derivata 4thDer"
set xlabel "t"  font ",14"
set ylabel "x , v , a"  font ",14"
set key         font ",14"
set lmargin at screen 0.12
set rmargin at screen 0.95
set bmargin at screen 0.12
set tmargin at screen 0.95

plot \
    "trajectory.dat" using 1:2 with points pt 7 ps 0.1 lc rgb "red"  title "x(t)", \
    "velocity4.dat" using 1:2 with points pt 7 ps 0.1 lc rgb "green" title "v(t)", \
    "acceleration.dat" using 1:2 with points pt 7 ps 0.1 lc rgb "blue" title "a(t)"