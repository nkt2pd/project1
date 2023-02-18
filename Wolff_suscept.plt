set xrange [2:2.8]
set yrange [0:1600]

plot "Wolff_susceptibility_30.dat" using 1:2 w lines, \
     "Wolff_susceptibility_60.dat" using 1:2 w lines, \
     "Wolff_susceptibility_120.dat" using 1:2 w lines, \
     "Wolff_susceptibility_240.dat" using 1:2 w lines, \
     "Wolff_suscpetibility_360.dat" using 1:2 w lines