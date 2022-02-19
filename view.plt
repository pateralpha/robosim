reset

fn = 'a.dat'

# get data size
stats fn nooutput
N = STATS_records
M = 6

# decralation of array
array x_c[N]
array y_c[N]
array t_c[N]

array x_m[N*M]
array y_m[N*M]

x_d = 2.0
y_d = 1.0

# save data to array
stats fn using ((x_c[$0 + 1] = $2, y_c[$0 + 1] = $3, t_c[$0 + 1] = $4, 0)) nooutput
stats fn using (sum[i = 1:M] (x_m[$0*M + i] = column(i*2 + 9), y_m[$0*M + i] = column(i*2 + 10), 0)) nooutput


# display data
array x_mt[M + 1]
array y_mt[M + 1]

len = 0.2
rad = 0.04

set xrange [-1:3]
set yrange [-1:3]
set size ratio -1

set parametric


# set terminal gif animate delay 5 optimize size 700, 700
# set output "view.gif"

do for [i = 1:N]{
    do for [j = 1:M]{
        x_mt[j] = x_m[(i - 1)*M + j]
        y_mt[j] = y_m[(i - 1)*M + j]
    }
    x_mt[M + 1] = x_mt[1]
    y_mt[M + 1] = y_mt[1]

    set arrow 1 from x_c[i], y_c[i] to (x_c[i] + len*cos(t_c[i])), (y_c[i] + len*sin(t_c[i])) lw 2 lc 'red'
    set arrow 2 from x_c[i], y_c[i] to (x_c[i] - len*sin(t_c[i])), (y_c[i] + len*cos(t_c[i])) lw 2 lc 'blue'

    plot y_c every ::1::i using (x_c[$1]):2 w l dt (8, 10) lw 2 lc 'web-blue' notitle, \
         y_mt using (x_mt[$1]):2 w l lw 3 lc 'gray50' notitle, \
         (rad*cos(t) + x_d), (rad*sin(t) + y_d) lw 1.5 lc 'blue' notitle, \
         (rad*cos(t) + x_c[1]), (rad*sin(t) + y_c[1]) lw 1.5 lc 'red' notitle
    pause 0.05
}

# set output