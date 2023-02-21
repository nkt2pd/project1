set xrange [2:2.8]
set yrange [0:3]

plot "Wolff_heat_30.dat" using 1:2 w lines, \
     "Wolff_heat_60.dat" using 1:2 w lines, \
     "Wolff_heat_120.dat" using 1:2 w lines #, \
     "Wolff_heat_240.dat" using 1:2 w lines, \
     "Wolff_heat_360.dat" using 1:2 w lines