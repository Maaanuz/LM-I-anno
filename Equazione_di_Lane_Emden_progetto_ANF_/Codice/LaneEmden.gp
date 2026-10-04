# ==================================================
# LaneEmden.gp
# Colonne nei .dat:
#   1: xi
#   2: theta(xi)
#   3: rho/rho_c = theta^n
# ==================================================

reset
set grid
set tics out
set key right top
set xrange [0:*]

# --- file di input ---
f0  = "laneemden_n0.dat"
f0a = "laneemden_n0_anal.dat"

f1  = "laneemden_n1.dat"
f1a = "laneemden_n1_anal.dat"

f15 = "laneemden_n1p5.dat"
f3  = "laneemden_n3.dat"
f4  = "laneemden_n4.dat"

f5  = "laneemden_n5.dat"
f5a = "laneemden_n5_anal.dat"

fomega = "laneemden_mass_coefficient.dat"
fconv = "RK4_convergence.dat"

# =====================
# FINESTRA 0: theta(xi)
# =====================
set term qt 0 size 900,650
set title "Lane-Emden: {/Symbol q}({/Symbol x})"
set xlabel "{/Symbol x}"
set ylabel "{/Symbol q}({/Symbol x})"
set xrange [0:20]
set yrange [*:*]

plot \
  f0  using 1:2 with lines lw 2 title "n = 0", \
  f1  using 1:2 with lines lw 2 title "n = 1", \
  f15 using 1:2 with lines lw 2 title "n = 3/2", \
  f3  using 1:2 with lines lw 2 title "n = 3", \
  f4  using 1:2 with lines lw 2 title "n = 4", \
  f5  using 1:2 with lines lw 2 title "n = 5"

# =====================
# FINESTRA 1: rho/rho_c
# =====================
set term qt 1 size 900,650
set title "Lane-Emden: {/Symbol r}/{/Symbol r}_c"
set xlabel "{/Symbol x}"
set ylabel "{/Symbol r}/{/Symbol r}_c"
set xrange [0:8]
set yrange [0:*]

plot \
  f0  using 1:3 with lines lw 2 title "n = 0", \
  f1  using 1:3 with lines lw 2 title "n = 1", \
  f15 using 1:3 with lines lw 2 title "n = 3/2", \
  f3  using 1:3 with lines lw 2 title "n = 3", \
  f4  using 1:3 with lines lw 2 title "n = 4", \
  f5  using 1:3 with lines lw 2 title "n = 5"

# ===============================================
# FINESTRA 2: confronto n=5 numerico vs analitico
# ===============================================
set term qt 2 size 900,650
set title "Lane-Emden n=5: confronto {/Symbol q}_{num} vs {/Symbol q}_{anal}"
set xlabel "{/Symbol x}"
set ylabel "{/Symbol q}({/Symbol x})"
set yrange [0:*]
set xrange [0:50]
set key right top

plot \
  f5 using 1:2 with lines lw 2 title "{/Symbol q}_{num}", \
  f5a using 1:2 with lines lw 2 title "{/Symbol q}_{anal}"

# ===============================================
# FINESTRA 3: confronto n=0 numerico vs analitico
# ===============================================
set term qt 3 size 900,650
set title "Lane-Emden n=0: confronto {/Symbol q}_{num} vs {/Symbol q}_{anal}"
set xlabel "{/Symbol x}"
set ylabel "{/Symbol q}({/Symbol x})"
set xrange [0:2.5]
set yrange [0:1.1]
set key right top

plot \
  f0  using 1:2 with lines lw 2 title "{/Symbol q}_{num}", \
  f0a using 1:2 with lines lw 2 title "{/Symbol q}_{anal}"

# ===============================================
# FINESTRA 4: confronto n=1 numerico vs analitico
# ===============================================
set term qt 4 size 900,650
set title "Lane-Emden n=1: confronto {/Symbol q}_{num} vs {/Symbol q}_{anal}"
set xlabel "{/Symbol x}"
set ylabel "{/Symbol q}({/Symbol x})"
set xrange [0:3.3]
set yrange [0:1.1]
set key right top

plot \
  f1  using 1:2 with lines lw 2 title "{/Symbol q}_{num}", \
  f1a using 1:2 with lines lw 2 title "{/Symbol q}_{anal}"

# =================================
# FINESTRA 5: coefficiente di massa
# =================================

set term qt 5 size 900,650
set title "-{/Symbol x}_s^2 {/Symbol q}'({/Symbol x}_s) in funzione di n"
set xlabel "n"
set ylabel "-{/Symbol x}_s^2 {/Symbol q}'({/Symbol x}_s)"
set xrange [0:6]
set yrange [*:*]
set grid
set key off

plot \
    fomega using 1:3 with linespoints lw 2 pt 7 ps 1 title ""

# ===================
# FINESTRA 6: xi_s(n)
# ===================

set term qt 6 size 900,650
set title "-{/Symbol x}_s in funzione di n"
set xlabel "n"
set ylabel "-{/Symbol x}_s"
set xrange [0:5]
set yrange [*:*]
set grid
set key off

plot \
    fomega using 1:2 with linespoints lw 2 pt 7 ps 1 title ""

# =============================================
# FINESTRA 7: convergenza RK4 (errore vs passo)
# =============================================

set term qt 7 size 900,650

set title "RK4 convergence: error vs {/Symbol D}{/Symbol x}"
set xlabel "{/Symbol D}{/Symbol x}"
set ylabel "maximum error"

unset xrange
unset yrange

set logscale xy
set grid
set key left top

# fit E = a*h^p

a = 1
p = 4

fit_func(x) = a*x**p
fit fit_func(x) fconv using 1:2 via a,p

print "---"
print "RK4 convergence test"
print "Order p = ", p
print "---"

plot \
fconv using 1:2 with points pt 7 title "numerical error", \
fit_func(x) with lines lw 2 title sprintf("fit: p = %.4f",p)

unset logscale xy

# ==================================
# FINESTRA 8: ordine numerico locale
# ==================================

set term qt 8 size 900,650

set title "RK4 measured order"
set xlabel "{/Symbol D}{/Symbol x}"
set ylabel "order p"

set xrange [*:*]
set yrange [0:6]

set grid
set key off

plot \
fconv using 1:4 with linespoints lw 2 pt 7 title "order"

pause -1 "Premi INVIO per chiudere."