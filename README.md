#### Orca Manuscript

This repository contains the code and data required for reproducing the analyses in the Orca manuscript. For using Orca and training Orca models, please visit our main [repository](https://github.com/jzhoulab/orca).

Most of the analyses code are provided in jupyter notebook format. Each jupyter notebook contains a series of analyses and typically generates multiple plots for the same theme of analyses. For large-scale virtual screens, you can find scripts under the `virtual_screen` directory, and the `Compartment_activity_screen.ipynb` also contains code for running multiple compartment activity screen analyses.

##### Dependencies
Other than [Orca dependencies](https://github.com/jzhoulab/orca), you will also need jupyter, rpy2, and plotnine python packages which can be installed with Anaconda or pip. For R packages, we will use data.table, ggplot2, patchwork, ggridges, ggrastr, ggthemes, limma.

##### Data
You will need additional resource files for reproducing some of the analyses, and we have provided these files [here](). In addition, you can download our precomputed results files [here](), which may save you time from executing the slow steps of the analyses.

