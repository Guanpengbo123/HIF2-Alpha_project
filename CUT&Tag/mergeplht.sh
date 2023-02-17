#!/bin/bash
#mkdir -p mergeheatmap/mat
#mkdir -p mergeheatmap/plot 

###GATA3##########

plotHeatmap -m mergeheatmap/mat/GATA3_center.mat.gz \
--whatToShow 'heatmap and colorbar' \
--heatmapWidth 8 \
-out mergeheatmap/plot/GAT3_heatmap.svg \
--colorList "white,#8560C4" --missingDataColor 1 \
--dpi 600 --zMax 200 --regionsLabel 'GATA3_peaks'

####HIF2####################

plotHeatmap -m mergeheatmap/mat/HIF2_center.mat.gz \
--whatToShow 'heatmap and colorbar' \
--heatmapWidth 8 \
-out mergeheatmap/plot/HIF2_heatmap.svg \
--colorList "white,#009A9A" --missingDataColor 1 \
--dpi 600 --zMax 200 --regionsLabel 'HIF2_peaks'

