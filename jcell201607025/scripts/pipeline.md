DATA
=====

Data downloaded March 16, 2016
[Paper](http://www.pnas.org/content/112/51/15672)
[GEO](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75140)
[SRA (ftp)](ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP066/SRP066834)
[IDs](http://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP066834)

Command:
cd sra/
ascp -i asperaweb_id_dsa.openssh -Q -T -k 2 anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByStudy/sra/SRP/SRP066/SRP066834/ .

# organize into one folder
mv SRP066834/*/*.sra .

# Too much data to process
# Just process Neurons and NPCs
# see recreate.R to get neurons and npcs
# reorganize into subfolders

DATA PROCESSING
=====

SRA to Fastq
------

Data converted to fastq format March 16, 2016
Using SRA toolkit

Command:
module load seq/sratoolkit/2.5.2
snakemake --snakefile /home/jf154/snakepit/sra/sra.snakefile --jobs 999 --cluster 'bsub -q short -W 12:00 -R "rusage[mem=4000]"'
mv *.fastq.gz fastq/
mv *.log out/


Alignment with Tophat
Quantification with HTSeq
-------

module load seq/htseq/0.6.1
snakemake --snakefile /home/jf154/snakepit/tophat/tophat.snakefile --jobs 999 --cluster 'bsub -q short -W 12:00 -R "rusage[mem=4000]"'


MISO
--------
# Merge neurons to one big bam
# Merge npcs to one big bam
samtools merge neurons.bam SRR2967613.bam SRR2967614.bam SRR2967655.bam SRR2967674.bam SRR2967685.bam SRR2967711.bam SRR2967726.bam SRR2967798.bam SRR2967830.bam SRR2967831.bam SRR2967832.bam SRR2967833.bam SRR2967717.bam SRR2967718.bam SRR2967719.bam SRR2967720.bam SRR2967721.bam SRR2967723.bam SRR2967725.bam SRR2967727.bam SRR2967729.bam SRR2967730.bam SRR2967732.bam SRR2967734.bam SRR2967737.bam SRR2967738.bam SRR2967741.bam SRR2967744.bam SRR2967745.bam SRR2967813.bam SRR2967814.bam SRR2967815.bam SRR2967817.bam SRR2967818.bam SRR2967819.bam SRR2967821.bam SRR2967822.bam SRR2967823.bam SRR2967824.bam SRR2967826.bam SRR2967827.bam SRR2967828.bam SRR2967829.bam SRR2967623.bam SRR2967624.bam SRR2967625.bam SRR2967626.bam SRR2967628.bam SRR2967631.bam SRR2967632.bam SRR2967633.bam SRR2967634.bam SRR2967636.bam SRR2967637.bam SRR2967638.bam SRR2967641.bam SRR2967642.bam SRR2967643.bam SRR2967645.bam SRR2967647.bam SRR2967747.bam SRR2967748.bam SRR2967749.bam SRR2967750.bam SRR2967751.bam SRR2967752.bam SRR2967753.bam SRR2967754.bam SRR2967755.bam SRR2967756.bam SRR2967758.bam SRR2967759.bam SRR2967762.bam SRR2967763.bam SRR2967765.bam SRR2967766.bam SRR2967767.bam SRR2967768.bam SRR2967770.bam SRR2967771.bam SRR2967772.bam SRR2967773.bam SRR2967774.bam SRR2967775.bam SRR2967776.bam SRR2967777.bam SRR2967779.bam SRR2967780.bam SRR2967781.bam SRR2967782.bam SRR2967784.bam SRR2967786.bam SRR2967788.bam SRR2967789.bam SRR2967790.bam SRR2967791.bam SRR2967792.bam SRR2967793.bam SRR2967795.bam SRR2967796.bam SRR2967797.bam SRR2967799.bam SRR2967800.bam SRR2967802.bam SRR2967803.bam SRR2967804.bam SRR2967805.bam SRR2967806.bam SRR2967807.bam SRR2967808.bam SRR2967809.bam SRR2967811.bam SRR2967608.bam SRR2967609.bam SRR2967610.bam SRR2967612.bam SRR2967615.bam SRR2967616.bam SRR2967617.bam SRR2967618.bam SRR2967620.bam SRR2967622.bam SRR2967650.bam SRR2967651.bam SRR2967652.bam SRR2967653.bam SRR2967654.bam SRR2967656.bam SRR2967658.bam SRR2967659.bam SRR2967664.bam SRR2967665.bam SRR2967667.bam SRR2967668.bam SRR2967669.bam SRR2967672.bam SRR2967673.bam SRR2967676.bam SRR2967679.bam SRR2967680.bam SRR2967681.bam SRR2967682.bam SRR2967684.bam SRR2967686.bam SRR2967687.bam SRR2967690.bam SRR2967691.bam SRR2967693.bam SRR2967694.bam SRR2967695.bam SRR2967696.bam SRR2967697.bam SRR2967698.bam SRR2967700.bam SRR2967701.bam SRR2967702.bam SRR2967703.bam SRR2967704.bam SRR2967705.bam SRR2967706.bam SRR2967707.bam SRR2967709.bam SRR2967712.bam SRR2967713.bam SRR2967716.bam

samtools merge npcs.bam SRR2967742.bam SRR2967761.bam SRR2967778.bam SRR2967816.bam SRR2967722.bam SRR2967724.bam SRR2967728.bam SRR2967731.bam SRR2967733.bam SRR2967735.bam SRR2967736.bam SRR2967739.bam SRR2967740.bam SRR2967743.bam SRR2967820.bam SRR2967825.bam SRR2967627.bam SRR2967629.bam SRR2967630.bam SRR2967635.bam SRR2967639.bam SRR2967640.bam SRR2967644.bam SRR2967646.bam SRR2967746.bam SRR2967757.bam SRR2967760.bam SRR2967764.bam SRR2967769.bam SRR2967783.bam SRR2967785.bam SRR2967787.bam SRR2967794.bam SRR2967801.bam SRR2967810.bam SRR2967611.bam SRR2967619.bam SRR2967621.bam SRR2967648.bam SRR2967649.bam SRR2967657.bam SRR2967660.bam SRR2967661.bam SRR2967662.bam SRR2967663.bam SRR2967666.bam SRR2967670.bam SRR2967671.bam SRR2967675.bam SRR2967677.bam SRR2967678.bam SRR2967683.bam SRR2967688.bam SRR2967689.bam SRR2967699.bam SRR2967708.bam SRR2967710.bam SRR2967714.bam SRR2967715.bam

samtools index npcs.bam
samtools index neurons.bam

### ALL GENES
# Compute PSI
module load seq/misopy/0.5.3
miso --run data-raw/miso_annotations_hg19_v2/indexed_ensGene_events/ tophat/npcs.bam --output-dir miso/npcs --settings-filename miso_settings.txt --read-len 100 --paired-end 250 15 --use-cluster
miso --run data-raw/miso_annotations_hg19_v2/indexed_ensGene_events/ tophat/neurons.bam --output-dir miso/neurons --settings-filename miso_settings.txt --read-len 100 --paired-end 250 15 --use-cluster

# Summarize the output
summarize_miso --summarize-samples miso/npcs/ miso/npcs/
summarize_miso --summarize-samples miso/neurons/ miso/neurons/

# Detect differentially expressed isoforms
compare_miso --compare-samples miso/npcs/ miso/neurons/ miso/comparisons/

### JUST SKIPPED EXONS
# Compute PSI
module load seq/misopy/0.5.3
miso --run data-raw/miso_annotations_hg19_v2/indexed_SE_events/ tophat/npcs.bam --output-dir miso/SE/npcs --settings-filename miso_settings.txt --read-len 100 --paired-end 250 15 --use-cluster
miso --run data-raw/miso_annotations_hg19_v2/indexed_SE_events/ tophat/neurons.bam --output-dir miso/SE/neurons --settings-filename miso_settings.txt --read-len 100 --paired-end 250 15 --use-cluster                                                                     
                                                          
# Summarize the output
summarize_miso --summarize-samples miso/SE/npcs/ miso/SE/npcs/
summarize_miso --summarize-samples miso/SE/neurons/ miso/SE/neurons/

# Detect differentially expressed isoforms                                                                
compare_miso --compare-samples miso/SE/npcs/ miso/SE/neurons/ miso/SE/comparisons/

# Filter?
#cd miso/SE/comparisons
#filter_events --filter npcs_vs_neurons.miso_bf --num-inc 1 --num-exc 1 --num-sum-inc-exc 10 --delta-psi 0.20 --bayes-factor 10 --output-dir filtered/

Analysis Take 2
=====

Using an unbiased single cell approach, we identify single cell subpopulations. We then pool these single cells (same number of cells per pool) together to identify differential alternative splicing. 

samtools merge neurons100sub19.bam SRR2967680.bam SRR2967681.bam SRR2967664.bam SRR2967674.bam SRR2967667.bam SRR2967668.bam SRR2967654.bam SRR2967616.bam SRR2967738.bam SRR2967713.bam SRR2967690.bam SRR2967669.bam SRR2967659.bam SRR2967767.bam SRR2967745.bam SRR2967697.bam SRR2967707.bam SRR2967737.bam SRR2967720.bam

samtools merge neuronsimm100sub19.bam SRR2967730.bam SRR2967673.bam SRR2967762.bam SRR2967685.bam SRR2967744.bam SRR2967655.bam SRR2967612.bam SRR2967693.bam SRR2967698.bam SRR2967620.bam SRR2967625.bam SRR2967732.bam SRR2967623.bam SRR2967768.bam SRR2967647.bam SRR2967754.bam SRR2967624.bam SRR2967631.bam SRR2967709.bam

samtools merge npcs100g1sub19.bam SRR2967635.bam SRR2967769.bam SRR2967733.bam SRR2967661.bam SRR2967646.bam SRR2967710.bam SRR2967629.bam SRR2967742.bam SRR2967663.bam SRR2967743.bam SRR2967683.bam SRR2967671.bam SRR2967630.bam SRR2967724.bam SRR2967657.bam SRR2967722.bam SRR2967764.bam SRR2967688.bam SRR2967660.bam

samtools merge npcs100g2sub19.bam SRR2967639.bam SRR2967675.bam SRR2967640.bam SRR2967708.bam SRR2967619.bam SRR2967739.bam SRR2967644.bam SRR2967648.bam SRR2967678.bam SRR2967761.bam SRR2967699.bam SRR2967731.bam SRR2967728.bam SRR2967736.bam SRR2967649.bam SRR2967689.bam SRR2967715.bam SRR2967760.bam SRR2967621.bam

module load seq/samtools/1.3
samtools sort -T /tmp/neuronsimm100sub19.sorted -o neuronsimm100sub19.sorted.bam neuronsimm100sub19.bam
samtools sort -T /tmp/neurons100sub19.sorted -o neurons100sub19.sorted.bam neurons100sub19.bam
samtools sort -T /tmp/npcs100g1sub19.sorted -o npcs100g1sub19.sorted.bam npcs100g1sub19.bam
samtools sort -T /tmp/npcs100g2sub19.sorted -o npcs100g2sub19.sorted.bam npcs100g2sub19.bam

samtools index neuronsimm100sub19.sorted.bam
samtools index neurons100sub19.sorted.bam
samtools index npcs100g1sub19.sorted.bam
samtools index npcs100g2sub19.sorted.bam


module load seq/misopy/0.5.3
miso --run data-raw/miso_annotations_hg19_v2/indexed_SE_events/ tophat/npcs100g1sub19.sorted.bam --output-dir miso3/npcsg1 --settings-filename miso_settings.txt --read-len 100 --paired-end 250 15 --use-cluster                                                                   

miso --run data-raw/miso_annotations_hg19_v2/indexed_SE_events/ tophat/npcs100g2sub19.sorted.bam --output-dir miso3/npcsg2 --settings-filename miso_settings.txt --read-len 100 --paired-end 250 15 --use-cluster                                                                   

miso --run data-raw/miso_annotations_hg19_v2/indexed_SE_events/ tophat/neurons100sub19.sorted.bam --output-dir miso3/neurons --settings-filename miso_settings.txt --read-len 100 --paired-end 250 15 --use-cluster 

miso --run data-raw/miso_annotations_hg19_v2/indexed_SE_events/ tophat/neuronsimm100sub19.sorted.bam --output-dir miso3/neuronsimm --settings-filename miso_settings.txt --read-len 100 --paired-end 250 15 --use-cluster 

summarize_miso --summarize-samples miso3/npcsg1/ miso3/npcsg1/
summarize_miso --summarize-samples miso3/npcsg2/ miso3/npcsg2/
summarize_miso --summarize-samples miso3/neuronsimm/ miso3/neuronsimm/
summarize_miso --summarize-samples miso3/neurons/ miso3/neurons/

compare_miso --compare-samples miso3/neuronsimm/ miso3/npcsg1/ miso3/comparisons/
compare_miso --compare-samples miso3/neuronsimm/ miso3/npcsg2/ miso3/comparisons/
compare_miso --compare-samples miso3/neurons/ miso3/npcsg1/ miso3/comparisons/
compare_miso --compare-samples miso3/neurons/ miso3/npcsg2/ miso3/comparisons/
compare_miso --compare-samples miso3/neurons/ miso3/neuronsimm/ miso3/comparisons/
compare_miso --compare-samples miso3/npcsg1/ miso3/npcsg2/ miso3/comparisons/

# Plot

module load seq/misopy/0.5.3
sashimi_plot --plot-event  "chr4:2927741:2927839:+@chr4:2928369:2928402:+@chr4:2929898:2931802:+"   data-raw/miso_annotations_hg19_v2/indexed_SE_events/ sashimi_plot_settings.txt --output-dir plots3


## Make same plots for bulk 

miso --run data-raw/miso_annotations_hg19_v2/indexed_SE_events/ /n/scratch2/xiao/hsCP_merge.bam.sorted.bam --output-dir miso_bulk/hsCP --settings-filename miso_settings.txt --read-len 46 --use-cluster 
miso --run data-raw/miso_annotations_hg19_v2/indexed_SE_events/ /n/scratch2/xiao/hsVZ_merge.bam.sorted.bam --output-dir miso_bulk/hsVZ --settings-filename miso_settings.txt --read-len 46 --use-cluster 

sashimi_plot --plot-event "chr14:51202234:51202335:-@chr14:51197641:51197701:-@chr14:51196241:51196441:-" data-raw/miso_annotations_hg19_v2/indexed_SE_events/ sashimi_plot_settings_bulk.txt --output-dir plots_bulk

sashimi_plot --plot-event "chr14:51202234:51202335:-@chr14:51197641:51197701:-@chr14:51196241:51196441:-" data-raw/miso_annotations_hg19_v2/indexed_SE_events/ sashimi_plot_settings.txt --output-dir plots3

summarize_miso --summarize-samples miso_bulk/hsVZ miso_bulk/hsVZ
summarize_miso --summarize-samples miso_bulk/hsCP miso_bulk/hsCP


## Select events of interest

module load seq/misopy/0.5.3
event="chrX:153581140:153581292:-@chrX:153580921:153581043:-@chrX:153580549:153580815:-"
event="chr3:58135577:58135729:+@chr3:58135832:58135954:+@chr3:58139102:58139368:+"
sashimi_plot --plot-event  $event   data-raw/miso_annotations_hg19_v2/indexed_SE_events/ sashimi_plot_settings.txt --output-dir plots3
sashimi_plot --plot-event $event data-raw/miso_annotations_hg19_v2/indexed_SE_events/ sashimi_plot_settings_bulk.txt --output-dir plots_bulk

