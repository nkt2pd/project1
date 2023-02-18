set xrange [2:2.8]
set yrange [0:0.8]

plot "Wolff_binder4_30.dat" using 1:2 w lp pt 1, \
     "Wolff_binder4_60.dat" using 1:2 w lp pt 2, \
     "Wolff_binder4_120.dat" using 1:2 w lp pt 3, \
     "Wolff_binder4_240.dat" using 1:2 w lp pt 4, \
     "Wolff_binder4_240.dat" using 1:2 w lp pt 5
