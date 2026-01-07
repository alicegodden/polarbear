# Diverging Transposon Activity Among Polar Bear Sub-Populations Inhabiting Different Climate Zones

[![DOI](https://img.shields.io/badge/DOI-10.1101/2024.12.04.626794-blue)](https://doi.org/10.1101/2024.12.04.626794)
[![nf-core/rnaseq](https://img.shields.io/badge/nf--core-rnaseq%203.9-brightgreen)](https://nf-co.re/rnaseq)
[![Zenodo](https://img.shields.io/badge/Zenodo-10.5281/zenodo.17078251-orange)](https://zenodo.org/record/17078251)

This repository contains code, data, and supplementary materials for the manuscript:  

📄 **Diverging transposon activity among polar bear sub-populations inhabiting different climate zones**  
*Godden, A.M., Rix, B.T. & Immler, S. Diverging transposon activity among polar bear sub-populations inhabiting different climate zones. Mobile DNA 16, 47 (2025). https://doi.org/10.1186/s13100-025-00387-4*  


👉 **Read the paper [here](https://link.springer.com/article/10.1186/s13100-025-00387-4)**

---

## 🔬 RNA-seq Pipeline
The RNA-seq pipeline was run using **[nf-core/rnaseq](https://nf-co.re/rnaseq)** version `3.9`.  
A complete list of software versions for all tools used in the pipeline is available here:  
[`software_versions_nfcore.yml`](./software_versions_nfcore.yml)

---

## 📂 Repository Structure

### `/R_scripts`
- **`Rnotebook_rnaseq.Rmd`** – All R scripts used for data generation and visualization  
- **`glmm_repeatlandscape.R`** – GLMM modeling of repeat landscape data  
- **`repeatlandscape_plots.R`** – Visualization of RepeatMasker outputs  
- **`multibam_summary.R`** – PCA plotting from MultiBam summary output  

---

### `/Python_scripts`
- **`bearphenogram.py`** – Plot a basic phenogram of significantly differentially expressed TEs (from Telescope outputs, genome ASM1731132v1)  
- **`bear_phenogram_TEs_and_genes_overlapping.py`** – Chromosomal plots of overlapping genes and TEs  
- **`bear_phenogram_genes_twolists.py`** – Chromosomal plots comparing two sets of gene loci  
- **`chrom_end_bear.txt`** – Chromosome/scaffold lengths (genome ASM1731132v1)  
- **`TE_matcher_genes.py`** – Match overlapping TEs and genes by genomic loci  
- **`bearTEA.py`** – Identify loci of significantly differentially expressed TEs given a TE list  
- **`bear_te_class.py`** – Bar chart visualization of TE classes  
- **`autobubble_goplot.py`** – Bubble plots of GO terms from ShinyGO outputs  

---

### `/supplementary_data`
Contains all supplementary data files referenced in the manuscript, including:  
- Gene and TE differential expression results  
- Metadata files  
- GO term analyses for significantly differentially expressed TEs overlapping genes  

⚠️ **Note:** *Supplementary File 4* exceeded GitHub’s file size limit.  
The full dataset can be accessed on Zenodo:  
[10.5281/zenodo.17078251](https://zenodo.org/record/17078251)  

- Full file: `Suppl.File4-TelescopeDESeq2dataRNA-seqASM1731132v1.csv`  
- Cropped version (rows with NA padj removed):  
  `Suppl.File4-TelescopeDESeq2dataRNA-seqASM1731132v1_cropped.csv`  

---

### `/Manuscript`
Includes the manuscript preprint:  
📄 **Diverging transposon activity among polar bear sub-populations inhabiting different climate zones**  
[bioRxiv 2024.12.04.626794](https://doi.org/10.1101/2024.12.04.626794)  

---

## 📑 Citation
If you use this code or data in your research, please cite:  

Godden, A.M., Rix, B.T. & Immler, S. Diverging transposon activity among polar bear sub-populations inhabiting different climate zones. Mobile DNA 16, 47 (2025). https://doi.org/10.1186/s13100-025-00387-4

---
