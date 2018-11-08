set term postscript eps  enhanced color 
set output 'plot.eps' 

set key font ",12"
set key spacing 4
set tics font ", 12"
#set  xtics 0,0.2,1
set ylabel 'u' font ", 12"
set  xlabel 'x' font ", 12"
#set yran[0.1:1.1]
set key right top
file = "output_soln.txt"
p   file u 1:2 t 'Forward Solution' w l lw 2 lc  rgb '#000000',\
   file u 1:3 t 'Backward Solution' w p pt 4  ps 2 lw 3 lc rgb '#0000FF', \
   file u 1:4 t 'Crank Nicolson Solution'  w p pt 1 ps 2 lw 3 lc  rgb '#000000', \
   file u 1:5 t 'Exact Solution'  w p pt 3 ps 2 lw 3 lc  rgb '#000000'