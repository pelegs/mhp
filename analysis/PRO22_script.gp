set terminal tikz standalone
set output "PRO22_graph.tex"
check_valid_column(c) = valid(c) ? check_valid_column(c + 1) : c - 1
stats 'PRO22.data' using (check_valid_column(1)) nooutput
last_column = int(STATS_max)

set grid
set xlabel 'time [frames]'
set ylabel 'MHP'

plot for [i=2:last_column-1] 'PRO22.data' u 1:i w l title '',\
				             'PRO22.data' using 1:last_column w l lw 4 lt 3 lc rgb '#FF0000' title 'Average'
