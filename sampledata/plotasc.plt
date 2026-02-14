

# Format of the output
pngflag=1

# OUTPUT PNG
if (pngflag==1) set term push
if (pngflag==1) set term pngcairo

set view map

set cbrange [0:1]

if (pngflag==1) set output "initial.png"
splot "snap000-00001.csv" u 1:2:8 w pm3d notitle \
     ,"snap001-00001.csv" u 1:2:8 w pm3d notitle \
     ,"snap002-00001.csv" u 1:2:8 w pm3d notitle \
     ,"snap003-00001.csv" u 1:2:8 w pm3d notitle

if (pngflag==1) set output "final.png"
splot "snap000-00101.csv" u 1:2:8 w pm3d notitle \
     ,"snap001-00101.csv" u 1:2:8 w pm3d notitle \
     ,"snap002-00101.csv" u 1:2:8 w pm3d notitle \
     ,"snap003-00101.csv" u 1:2:8 w pm3d notitle



