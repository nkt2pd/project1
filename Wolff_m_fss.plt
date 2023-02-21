set xrange [-20:30]
set yrange [0:1.6]

plot "Wolff_m_fss.dat" using 1:2 w p pt 1 lc rgb "green", \
     "Wolff_m_fss.dat" using 3:4 w p pt 2 lc rgb "blue", \
     "Wolff_m_fss.dat" using 5:6 w p pt 3 lc rgb "red"