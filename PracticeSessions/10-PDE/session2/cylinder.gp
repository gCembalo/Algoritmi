reset

set pm3d map

set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')

unset colorbox
set colorbox vertical
set cbtics font ",12"

set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax

set palette defined
set palette model HSV defined ( 0 0 1 1, 1 1 1 1)

set tics font ",14"
set title "Cylinder"
set xlabel "x"  font ",14"
set ylabel "y"  font ",14"

set lmargin at screen 0.1
set rmargin at screen 0.82
set bmargin at screen 0.12
set tmargin at screen 0.95

splot "cylinder.dat" u 1:2:3