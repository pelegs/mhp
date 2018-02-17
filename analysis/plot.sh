#!/bin/bash
for f in *.gp
do
    filename=$(basename "$f")
    INPUT="${filename%.*}"
    filename=$(echo $INPUT| cut -d'_' -f 1)
    gnuplot $filename"_script.gp" && pdflatex $filename"_graph.tex" && evince $filename"_graph.pdf"
done
