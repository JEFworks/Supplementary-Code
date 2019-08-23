# Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain

Blue B Lake*, Song Chen*, Brandon C Sos*, Jean Fan*, Gwendolyn E Kaeser, Yun C Yung, Thu E Duong, Derek Gao, Jerold Chun, Peter V Kharchenko, Kun Zhang

Abstract: Detailed characterization of the cell types in the human brain requires scalable experimental approaches to examine multiple aspects of the molecular state of individual cells, as well as computational integration of the data to produce unified cell-state annotations. Here we report improved high-throughput methods for single-nucleus droplet-based sequencing (snDrop-seq) and single-cell transposome hypersensitive site sequencing (scTHS-seq). We used each method to acquire nuclear transcriptomic and DNA accessibility maps for >60,000 single cells from human adult visual cortex, frontal cortex, and cerebellum. Integration of these data revealed regulatory elements and transcription factors that underlie cell-type distinctions, providing a basis for the study of complex processes in the brain, such as genetic programs that coordinate adult remyelination. We also mapped disease-associated risk variants to specific cellular populations, which provided insights into normal and pathogenic cellular processes in the human brain. This integrative multi-omics approach permits more detailed single-cell interrogation of complex organs and tissues. 

# Data

Raw sequencing data, annotated snDrop-seq and scTHS-seq count matrices, and DNA accessibility peak files are all available from the Gene
Expression Omnibus under SuperSeries accession code GSE97942: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97942

Note, the SuperSeries references the 2 SubSeries which comprises the scTHS-seq and the snDrop-seq datasets separately:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2579603
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2579604

Additional related data may be downloaded from: http://pklab.med.harvard.edu/peterk/kun/blue/nbt2018/ 

# Rmd

Cleaned up Rmarkdown files recreating main analyses and figures can be found in the Rmd/ folder, courtesy of Masahiro Kanai
