

# mpi process
snap_start = 0
snap_end   = 3

# 
step_init  = 1    
step_final = 101

# Format of the output
pngflag=1

if (pngflag==1) set term push
if (pngflag==1) set term pngcairo

set view map
set cbrange [0:1]

# =========================
# initial
# =========================
if (pngflag==1) set output "initial.png"

splot for [i=snap_start:snap_end] \
    sprintf("snap%03d-%05d.csv", i, step_init) \
    u 1:2:8 w pm3d notitle

# =========================
# final
# =========================
if (pngflag==1) set output "final.png"

splot for [i=snap_start:snap_end] \
    sprintf("snap%03d-%05d.csv", i, step_final) \
    u 1:2:8 w pm3d notitle

