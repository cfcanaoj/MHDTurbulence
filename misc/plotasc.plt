

dir="./"

# mpi process
snap_start = 0
snap_end   = 3

# 
step_init  = 1    
step_final = 101

print "dir=".dir
print "mpi process=".snap_start."-".snap_end
print "step=".step_init." and ".step_final

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
if (pngflag==1) print output 

splot for [i=snap_start:snap_end] \
    sprintf(dir."snap%03d-%05d.csv", i, step_init) \
    u 1:2:8 w pm3d notitle

# =========================
# final
# =========================
if (pngflag==1) set output "final.png"
if (pngflag==1) print output

splot for [i=snap_start:snap_end] \
    sprintf(dir."snap%03d-%05d.csv", i, step_final) \
    u 1:2:8 w pm3d notitle

