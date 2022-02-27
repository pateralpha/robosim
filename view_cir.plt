reset

fn = 'a.dat'
fc = 'c.dat'

# get data size
stats fn nooutput

N = STATS_records
D = 21
M = 6
W = 3
T = 4

stats fc nooutput
Nd = STATS_records

# decralation of array
array x_c[N]
array y_c[N]
array t_c[N]

array x_d[N]
array y_d[N]

array x_m[N*M]
array y_m[N*M]

array x_t[N*4*3]
array y_t[N*4*3]

array time[N]

array type_cnt[Nd + 1]

array x_cnt[Nd + 1]
array y_cnt[Nd + 1]

array ts_cnt[Nd + 1]
array te_cnt[Nd + 1]
array tee_cnt[Nd + 1]

array r_cnt[Nd + 1]

array x_cnt2[Nd*2 + 1]
array y_cnt2[Nd*2 + 1]

# save data to array
stats fn using (time[$0 + 1] = $1, 0) nooutput
stats fn using (x_c[$0 + 1] = $2, y_c[$0 + 1] = $3, t_c[$0 + 1] = $4, 0) nooutput
stats fn using (x_d[$0 + 1] = $8, y_d[$0 + 1] = $9, 0) nooutput
stats fn using (sum[i = 1:M] (x_m[$0*M + i] = column((i - 1)*2 + D + 1), y_m[$0*M + i] = column((i - 1)*2 + D + 2), 0)) nooutput
stats fn using (sum[i = 1:(T*W)] (x_t[$0*T*W + i] = column((i - 1)*2 + D + M*2 + 1), y_t[$0*T*W + i] = column((i - 1)*2 + D + M*2 + 2), 0)) nooutput

stats fc using (type_cnt[$0 + 2] = $1, x_cnt[$0 + 2] = $2, y_cnt[$0 + 2] = $3, r_cnt[$0 + 2] = $4, ts_cnt[$0 + 2] = $5, te_cnt[$0 + 2] = $6, 0) nooutput
type_cnt[1] = 0
x_cnt[1] = x_c[1]
y_cnt[1] = y_c[1]
r_cnt[1] = 0
ts_cnt[1] = 0
te_cnt[1] = 0

x_cnt2[1] = x_cnt[1]
y_cnt2[1] = y_cnt[1]

pt = 2

# print type_cnt

do for[i = 1:Nd + 1]{
    tee_cnt[i] = te_cnt[i]

    __buf = ts_cnt[i] > te_cnt[i] ? ts_cnt[i] : te_cnt[i];
    ts_cnt[i] = ts_cnt[i] < te_cnt[i] ? ts_cnt[i] : te_cnt[i];
    te_cnt[i] = __buf > te_cnt[i] ? __buf : te_cnt[i];

    if(type_cnt[i] != 0){
        pt = pt + 1
        x_cnt2[pt] = x_cnt[i] + r_cnt[i]*cos(tee_cnt[i])
        y_cnt2[pt] = y_cnt[i] + r_cnt[i]*sin(tee_cnt[i])
    } else {
        x_cnt2[pt] = x_cnt[i]
        y_cnt2[pt] = y_cnt[i]
    }
    pt = pt + 1
}

# print ts_cnt
# print te_cnt

# display data
array x_mt[M + 1]
array y_mt[M + 1]

array x_tt[(T + 1)*W]
array y_tt[(T + 1)*W]

array x_obs[7] = [-2, 0, 1, 1, 0, -2, -2]
array y_obs[7] = [1, 1, 2, 3, 4, 4, 1]

len = 0.2
rad = 0.05

set parametric
set sample 500

# set terminal gif animate delay 5 optimize size 500, 500 font 'TimesNewRoman, 12'
# set output "view.gif"

set terminal qt size 500, 500 font 'TimesNewRoman, 12'

# set xrange [-1:3]
# set yrange [-1:3]
set xrange [-2:4]
set yrange [-1:5]
set size ratio -1

set trange [-pi:pi]

set lmargin at screen 0.13
set bmargin at screen 0.12

set xlabel '{/TimesNewRoman:Italic x}-coordinate [m]' offset 0, 0.3
set ylabel '{/TimesNewRoman:Italic y}-coordinate [m]' offset 1.8, 0

set tics font 'TimesNewRoman, 10'

set grid
set border lw 1.2

do for [i = 1:N]{
    do for [j = 1:M]{
        x_mt[j] = x_m[(i - 1)*M + j]
        y_mt[j] = y_m[(i - 1)*M + j]
    }
    x_mt[M + 1] = x_mt[1]
    y_mt[M + 1] = y_mt[1]

    do for [k = 0:(W - 1)]{
        do for [j = 1:T]{
            x_tt[k*(T + 1) + j] = x_t[(i - 1)*T*W + k*T + j]
            y_tt[k*(T + 1) + j] = y_t[(i - 1)*T*W + k*T + j]
        }
        x_tt[(k + 1)*(T + 1)] = x_tt[k*(T + 1) + 1]
        y_tt[(k + 1)*(T + 1)] = y_tt[k*(T + 1) + 1]
    }

    set label 1 at graph 0.7, 0.07 sprintf('time = %4.1f [s]', time[i])

    set arrow 1 from x_c[i], y_c[i] to (x_c[i] + len*cos(t_c[i])), (y_c[i] + len*sin(t_c[i])) lw 2 lc 'red'
    set arrow 2 from x_c[i], y_c[i] to (x_c[i] - len*sin(t_c[i])), (y_c[i] + len*cos(t_c[i])) lw 2 lc 'blue'

    plot y_c every ::1::i using (x_c[$1]):2 w l dt (8, 10) lw 2 lc 'web-blue' notitle, \
        y_mt using (x_mt[$1]):2 w l lw 3 lc 'gray50' notitle, \
        y_tt every ::0::T using (x_tt[$1]):2 w l lw 3 lc 'gray30' notitle, \
        y_tt every ::(T + 1)::((T + 1)*2 - 1) using (x_tt[$1]):2 w l lw 3 lc 'gray30' notitle, \
        y_tt every ::((T + 1)*2)::((T + 1)*3 - 1) using (x_tt[$1]):2 w l lw 3 lc 'gray30' notitle, \
        y_cnt using (type_cnt[$1] == 0 ? x_cnt[$1] : x_cnt[$1] + r_cnt[$1]*cos(tee_cnt[$1])):(type_cnt[$1] == 0 ? y_cnt[$1] : y_cnt[$1] + r_cnt[$1]*sin(tee_cnt[$1])) w p pointtype 6 lc 'gray70' notitle, \
        y_cnt2 using (x_cnt2[$1]):2 w l lc 'gray70' notitle, \
        for[i = 1:(Nd + 1)] (((t < ts_cnt[i])||(t > te_cnt[i]))? 1/0 : x_cnt[i] + r_cnt[i]*cos(t)), (((t < ts_cnt[i])||(t > te_cnt[i]))? 1/0 : y_cnt[i] + r_cnt[i]*sin(t)) w l lc 'gray70' notitle, \
        (rad*cos(t) + x_d[i]), (rad*sin(t) + y_d[i]) lw 2 lc 'blue' notitle, \
        (rad*cos(t) + x_c[1]), (rad*sin(t) + y_c[1]) lw 2 lc 'red' notitle, \
        y_obs using (x_obs[$1]):2 w filledc x2 fill pattern 4 lc 'gray40' notitle, \
        y_obs using (x_obs[$1]):2 w l lw 2 lc 'gray40' notitle
    pause 0.05
}

# set output