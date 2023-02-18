set xrange [2:2.8]
set yrange [0:1]

plot "Wolff_m_30.dat" using 1:2 w lines, \
     "Wolff_m_60.dat" using 1:2 w lines, \
     "Wolff_m_120.dat" using 1:2 w lines, \
     "Wolff_m_240.dat" using 1:2 w lines, \
     "Wolff_m_360.dat" using 1:2 w lines