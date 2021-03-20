#### Virtual screens

This directory contains scripts and files needed for running large-scale virtual screens peformed in the Orca manuscript. The main scripts are `compartment_activity_screen.py` (compartment activity screen) and `local_interaction_screen.py` (submegabase pair local activity screen) and the corresponding `.sh` files contains the commands to run the screens (you will likely want to distribute these jobs across multiple machines). Note that here only contains scripts for the large-scale 12800bp segment size compartment activity screens, see the jupyter notebook `../Compartment_activity_screen.ipynb` for all other compartment activity screen analyses.

Here we provide commands for postprocessing the virtual screen outputs after executing the `.sh` files, the output is also available for download (see README if this manuscript respository) if you would like to skip these steps. There are a couple of additional dependencies required, namely BEDOPS and bedGraphToBigWig (UCSC).
```bash
cat local_interaction.*h1esc_1m.bedgraph > local_interaction.h1esc_1m.bedgraph
cat local_interaction.*hff_1m.bedgraph > local_interaction.hff_1m.bedgraph

awk '$4>0.01 {print}' local_interaction.h1esc_1m.bedgraph > ../figure_data/local_interaction.h1esc_1m.0.01.bedgraph
awk '$4>0.01 {print}' local_interaction.hff_1m.bedgraph > ../figure_data/local_interaction.hff_1m.0.01.bedgraph

shuf local_interaction.h1esc_1m.bedgraph|head -n 100000 |sort-bed - > ../figure_data/local_interaction.h1esc_1m.shuf100k.bedgraph
cut -f 1,2,3 local_interaction.h1esc_1m.shuf100k.bedgraph|fgrep -w -f - local_interaction.hff_1m.bedgraph |sort-bed - > ../figure_data/local_interaction.hff_1m.shuf100k.bedgraph

sort-bed local_interaction.h1esc_1m.bedgraph > local_interaction.h1esc_1m.bedgraph.sorted
sort-bed local_interaction.hff_1m.bedgraph > local_interaction.hff_1m.bedgraph.sorted
bedGraphToBigWig local_interaction.h1esc_1m.bedgraph.sorted http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes ../figure_data/local_interaction.h1esc_1m.bedgraph.sorted.bw
bedGraphToBigWig local_interaction.hff_1m.bedgraph.sorted http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes ../figure_data/local_interaction.hff_1m.bedgraph.sorted.bw
```