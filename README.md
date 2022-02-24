## Orca Manuscript

This repository contains the code and data required for reproducing the analyses in the Orca manuscript. For using Orca and training Orca models, please visit our main [repository](https://github.com/jzhoulab/orca). For running Orca from our webserver, please visit [orca.zhoulab.io](https://orca.zhoulab.io).

Most of the analyses code are provided in jupyter notebook format. Each jupyter notebook contains a series of analyses and typically generates multiple plots for the same theme of analyses. For large-scale virtual screens, you can find scripts under the `virtual_screen` directory, and the `Compartment_activity_screen.ipynb` also contains code for running multiple compartment activity screen analyses.

The jupyter notebooks are largely grouped by topics:

- Compartment_activity_screen.ipynb : Analyses of sequence activities in chromatin compartment alterations
- Local_interaction_screens.ipynb : Identifying sequences that affect submegabase-scale genome interactions
- Local_interaction_multiplex_demo.ipynb : Multiplexed in silico mutagenesis demo example.
- Enhancerpolycomb_example.ipynb : Example predictions of enhancer-promoter interactions and Polycomb-mediated interactions.
- Model_performance.ipynb : Orca model prediction performance evaluations.
- Model_performance_256M.ipynb : Orca model prediction performance evaluations (32-256Mb).
- Model_performance_hctnoc.ipynb : Orca model prediction performance evaluations for the cohesin-depleted HCT119 model.
- StructuraVariants.ipynb : Prediction of structural variant effects on 3D genome interactions. 


### Dependencies
Other than [Orca dependencies](https://github.com/jzhoulab/orca#installation), you will also need jupyter, rpy2, and plotnine python packages which can be installed with Anaconda or pip. For R packages, we will use data.table, ggplot2, patchwork, ggridges, ggrastr, ggthemes, and limma.

### Data
You will need additional resource files for reproducing some of the analyses, and we have provided these files [here](https://zenodo.org/record/4642212/files/orca_manuscript_resources.tar.gz)(2.9G). Note that the jupyter notebooks use GPU to generate Orca predictions by default. You can generally switch to CPU by using `use_cuda=False` option, but they may be too slow for computationally intensive steps. You can also skip the computationally intensive steps by downloading our precomputed results files [here](https://zenodo.org/record/4642212/files/orca_manuscript_figuredata.tar.gz)(19.9G). We have provided code to load precomputed results in the jupyter notebooks.

