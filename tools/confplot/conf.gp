set output 'plot.pdf'
set term pdfcairo font ",6"
set size ratio 1
# set palette gray negative
set palette rgb 21,22,23 negative
set autoscale xfix
set autoscale yfix
set xtics 1
set ytics 1
set title "contribution per path configuration"
set ylabel "# path tracing vertices"
set xlabel "# light tracing vertices"

set tics scale 0,0.001
set mxtics 2
set mytics 2
set grid front mxtics mytics lw 1.5 lt -1 lc rgb 'white'
plot "work.dat" matrix w image noti
