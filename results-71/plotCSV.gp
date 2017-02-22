set terminal eps;
set output 'Results.eps';
unset key;
set key off;
set loadpath '/home/staff/scherping/Lib/Input/';
set xtics ("0.2" 20, "0.4" 40, "0.6" 60, "0.8" 80, "1" 100);
set yrange [-0.2:1];
set size 1,1;
set origin 0,0;
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0;
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0;
set multiplot layout 2,2 columnsfirst scale 1.1,0.9;
plot 'plot1.csv' with lines;
plot 'uRshifted.csv' with lines, 'dataFunction.csv' with lines;
plot 'plot2.csv' with lines;
plot 'thresholdedSolution.csv' with lines;
print "All done!";

