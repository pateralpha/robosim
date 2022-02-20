reset

fn = 'a.dat'

D = 6

plot fn u 1:(column(D + 1)) w lp lw 1.5 title 'e_1', \
    fn u 1:(column(D + 2)) w lp lw 1.5 title 'e_2', \
    fn u 1:(column(D + 3)) w lp lw 1.5 title 'e_3'