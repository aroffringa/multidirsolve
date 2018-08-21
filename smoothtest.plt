set terminal postscript enhanced color
#set logscale y
#set xrange [0.1:]
#set yrange [0:3.1415926535*2]
set output "plot-kerneltest.ps"
#set key bottom left
set xlabel "Frequency (MHz)"
set ylabel "Amplitude"
plot \
"build/kerneltest.txt" using 1:2 with lines title 'Function' lw 2.0, \
"build/kerneltest.txt" using 1:3 with points title 'Input' lw 2.0, \
"build/kerneltest.txt" using 1:4 with lines title 'Rectangular' lw 2.0, \
"build/kerneltest.txt" using 1:5 with lines title 'Triangular' lw 2.0, \
"build/kerneltest.txt" using 1:6 with lines title 'Gaussian' lw 2.0, \
"build/kerneltest.txt" using 1:7 with lines title 'Quadratic' lw 2.0

