set terminal postscript eps enhanced
set out 'P_rho_plot.eps'
set title "Plot of P against {/Symbol r}"
set xlabel "{/Symbol r}"
set ylabel "P"
set xrange [0.0:1.0]
plot "P_plot.out" using 1:2 smooth csplines title "Literature Values", "P_plot.out" using 1:3 title "Calculated Values"