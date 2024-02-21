set terminal pngcairo size 800,450
set output 'example-spectrum.png'

# Set bar width
set boxwidth 0.5
set style fill solid border -1
set linetype 1 lc rgb "#990000" lw 1
set linetype 2 lc rgb "#009900" lw 1
set linetype 3 lc rgb "#666666" lw 1

set arrow from graph 0,1 to graph 0,1.1 filled
set arrow from graph 1,0 to graph 1.1,0 filled
set tmargin 5
set rmargin 10
set border 3
set tics nomirror
set noxtics
set noytics
set grid
set ylabel "Magnitude"
set xlabel "Frequency"

# Set the range for y-axis
set yrange [0:*]

# Plot the data with different colors
plot "example-spectrum.dat" using 1:2:3 with boxes lc variable notitle, \
     "<grep '1$' example-spectrum.dat" u 1:2:(0):(1500) with vectors lc 1 notitle, \
     "<grep '2$' example-spectrum.dat" u 1:($2+1500):(0):(-1500) with vectors lc 2 notitle
