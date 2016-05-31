[Burger JA, Landau DA, Taylor-Weiner A, Bozic I, Zhang H, Sarosiek K, Wang L, Stewart C, Fan J, Hoellenriegel J, et al. Clonal evolution in patients with chronic lymphocytic leukaemia developing resistance to BTK inhibition. Nat Commun. 2016;7](http://www.nature.com/ncomms/2016/160520/ncomms11589/full/ncomms11589.html)

# Code for recreating Figures 1C and Supp 3B:
![fig 1c](https://raw.githubusercontent.com/JEFworks/figure-code-dump/master/ncomms11589/img/fig_1c.png)
![supp 3b](https://raw.githubusercontent.com/JEFworks/figure-code-dump/master/ncomms11589/img/supp_3b.png)

## scripts/

*get_mats.R* - processes Fluidigm Biomark(TM) multiplex qRT-PCR HeatMap result file, saving into convenient RData format for additional analysis

*mutation_call.R* - performs mutation calling as described in the "Microfluidics-based single-cell detection of mutations" methods section of the manuscript. Relevant code for recreating Figure 1C and Supp 3B. Note not all code/analysis is relevant to Figures 1C and Supp 3B. Addition reorganization and cleaning of figures was done in Illustrator and thus may not aesthetically match reproductions.
