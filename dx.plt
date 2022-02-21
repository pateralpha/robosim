reset

fn = 'a.dat'

D = 4
D2 = 9
D3 = 13

plot fn u 1:(column(D + 1)) w lp lw 1.5 title 'dx', \
    fn u 1:(column(D + 2)) w lp lw 1.5 title 'dy', \
    fn u 1:(column(D3 + 1)) w l lw 1.5 dt 4 title 'dx_c', \
    fn u 1:(column(D3 + 2)) w l lw 1.5 dt 4 title 'dy_c', \
    fn u 1:(column(D2 + 1)) w l lw 1.5 title 'dx_d', \
    fn u 1:(column(D2 + 2)) w l lw 1.5 title 'dy_d'