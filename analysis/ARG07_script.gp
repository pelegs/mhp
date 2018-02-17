set terminal tikz standalone
set output "ARG07_graph.tex"
check_valid_column(c) = valid(c) ? check_valid_column(c + 1) : c - 1
stats 'ARG07.data' using (check_valid_column(1)) nooutput
last_column = int(STATS_max)

set grid
set xlabel 'time [frames]'
set ylabel 'MHP'

plot for [i=2:last_column-1] 'ARG07.data' u 1:i w l title '',\
				             'ARG07.data' using 1:last_column w l lw 4 lt 3 lc rgb '#FF0000' title 'Average'
