reset

fn = 'a.dat'

D = 4
D2 = 9
D3 = 13

dx_o = 0
dy_o = 0
xd_o = 0
yd_o = 0

plot fn u 1:(dx = column(D + 1) - dx_o, dx_o = column(D + 1), dx*2) w lp lw 1.5 title 'dx', \
    fn u 1:(dy = column(D + 2) - dy_o, dy_o = column(D + 2), dy*20) w lp lw 1.5 title 'dy', \
    fn u 1:(xd = column(D2 + 1) - xd_o, xd_o = column(D2 + 1), xd*20) w l lw 1.5 title 'dx_d', \
    fn u 1:(yd = column(D2 + 2) - yd_o, yd_o = column(D2 + 2), yd*20) w l lw 1.5 title 'dy_d'