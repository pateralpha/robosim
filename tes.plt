gyou=100 #行番号を指定
retsu=1 #列番号を指定
filename="a.dat" # ファイル名を指定

# mx1 = system("cat " . filename . " | awk \'NR==" . gyou ."{printf(\"%f\", $" . 8 . ")}\'")
# my1 = system("cat " . filename . " | awk \'NR==" . gyou ."{printf(\"%f\", $" . 9 . ")}\'")
# mx2 = system("cat " . filename . " | awk \'NR==" . gyou ."{printf(\"%f\", $" . 10 . ")}\'")
# my2 = system("cat " . filename . " | awk \'NR==" . gyou ."{printf(\"%f\", $" . 11 . ")}\'")
# mx3 = system("cat " . filename . " | awk \'NR==" . gyou ."{printf(\"%f\", $" . 12 . ")}\'")
# my3 = system("cat " . filename . " | awk \'NR==" . gyou ."{printf(\"%f\", $" . 13 . ")}\'")
# mx4 = system("cat " . filename . " | awk \'NR==" . gyou ."{printf(\"%f\", $" . 14 . ")}\'")
# my4 = system("cat " . filename . " | awk \'NR==" . gyou ."{printf(\"%f\", $" . 15 . ")}\'")
# mx5 = system("cat " . filename . " | awk \'NR==" . gyou ."{printf(\"%f\", $" . 16 . ")}\'")
# my5 = system("cat " . filename . " | awk \'NR==" . gyou ."{printf(\"%f\", $" . 17 . ")}\'")
# mx6 = system("cat " . filename . " | awk \'NR==" . gyou ."{printf(\"%f\", $" . 18 . ")}\'")
# my6 = system("cat " . filename . " | awk \'NR==" . gyou ."{printf(\"%f\", $" . 19 . ")}\'")

# set label 1 mx1*c-my1*s mx1*s+my1*c
mx1(w) = system("cat " . filename . " | awk \'NR==" . floor(w) ."{printf(\"%f\", $" . 8 . ")}\'")
mx2(w) = system("cat " . filename . " | awk \'NR==" . floor(w) ."{printf(\"%f\", $" . 10 . ")}\'")
mx3(w) = system("cat " . filename . " | awk \'NR==" . floor(w) ."{printf(\"%f\", $" . 12 . ")}\'")
mx4(w) = system("cat " . filename . " | awk \'NR==" . floor(w) ."{printf(\"%f\", $" . 14 . ")}\'")
mx5(w) = system("cat " . filename . " | awk \'NR==" . floor(w) ."{printf(\"%f\", $" . 16 . ")}\'")
mx6(w) = system("cat " . filename . " | awk \'NR==" . floor(w) ."{printf(\"%f\", $" . 18 . ")}\'")

mx(w, i) = real(floor(w) <= 1 ? mx1(i) :(floor(w) == 2 ? mx2(i) :(floor(w) == 3 ? mx3(i) :(floor(w) == 4 ? mx4(i) :(floor(w) == 5 ? mx5(i) :(floor(w) == 6 ? mx6(i) : mx1(i)))))))

my1(w) = system("cat " . filename . " | awk \'NR==" . floor(w) ."{printf(\"%f\", $" . 9 . ")}\'")
my2(w) = system("cat " . filename . " | awk \'NR==" . floor(w) ."{printf(\"%f\", $" . 11 . ")}\'")
my3(w) = system("cat " . filename . " | awk \'NR==" . floor(w) ."{printf(\"%f\", $" . 13 . ")}\'")
my4(w) = system("cat " . filename . " | awk \'NR==" . floor(w) ."{printf(\"%f\", $" . 15 . ")}\'")
my5(w) = system("cat " . filename . " | awk \'NR==" . floor(w) ."{printf(\"%f\", $" . 17 . ")}\'")
my6(w) = system("cat " . filename . " | awk \'NR==" . floor(w) ."{printf(\"%f\", $" . 19 . ")}\'")

my(w, i) = real(floor(w) <= 1 ? my1(i) :(floor(w) == 2 ? my2(i) :(floor(w) == 3 ? my3(i) :(floor(w) == 4 ? my4(i) :(floor(w) == 5 ? my5(i) : (floor(w) == 6 ? my6(i) : my1(i)))))))


# plot filename u 1:2 with linespoint
# replot filename u 1:3 with linespoint

set parametric
set samples 7

set xrange [-1:3]
set yrange [-1:3]

set size ratio -1

set terminal gif animate delay 5
set output "tes.gif"

do for [i = 2:61]{
    plot [1:7] mx(t, i), my(t, i)
}

set output

# plot mx1