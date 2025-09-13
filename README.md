# Glacial ice subpopulation of polar bears exhibits divergent transposon activity

RNA-seq pipeline was run according to nf-core rnaseq -r 3.9. Full software versions of all tools in the RNA-seq pipeline can be found here:
software_versions_nfcore.yml

All R scripts, notebooks and python scripts and files can be found here:

# /R_scripts:
Rnotebook_rnaseq.Rmd # All R scripts used for generation and analysis of data/plots
glmm_repeatlandscape.R # GLMM modelling of repeat landscape data 

glmm_repeatlandscape.R # Plotting repeatlandscape from RepeatMasker outputs

multibam_summary.R # PCA plotting from multibam summary output

# /Python scripts:
bearphenogram.py # Plotting a basic phenogram of a list of significantly differentially expressed TEs from Telescope outputs, genome ASM1731132v1

bear_phenogram_TEs_and_genes_overlapping.py # Plotting genes and TEs in a chromosomal plot

bear_phenogram_genes_twolists.py # Plotting two input sets of gene loci in a chromsomal plot

chrom_end_bear.txt # Chromosome/scaffold legnths text file, genome ASM1731132v1

TE_matcher_genes.py # Matching overlapping TEs and genes by loci

bearTEA.py # Used for finding loci of significantly differentially expressed TEs with a list of TEs as input

bear_te_class.py # Plotting bar charts of TEs

autobubble_goplot.py # Plotting bubble plots of GO terms from ShinyGO outputs

# /supplementary_data :
Contains all supplementary data files for manuscript. These include gene and TE differential expression data and metadata. 
Also includes GO terms analysis of significantly differentially expressed TEs that overlap genes in the reference genome.
Suppl. File 6. was over memory limit for file size full file can be found along with other annotation data here: 10.5281/zenodo.17078251 - Suppl.File6-TelescopeDESeq2dataRNA-seqASM1731132v1.csv .


Otherwise please see Suppl File. 6 with NA padj rows removed "cropped" > Suppl.File6-TelescopeDESeq2dataRNA-seqASM1731132v1_cropped.csv

# / Manuscript :
Diverging transposon activity among polar bear sub-populations inhabiting different climate zones
bioRxiv 2024.12.04.626794; doi: https://doi.org/10.1101/2024.12.04.626794 
