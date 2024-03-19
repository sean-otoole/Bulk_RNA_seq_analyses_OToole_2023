# Bulk_RNA_seq_analyses_OToole_2023

Bulk RNA-sequencing analyses for the publication entitled: [Molecularly targetable cell types in mouse visual cortex have distinguishable prediction error responses](https://www.cell.com/neuron/pdf/S0896-6273(23)00626-8.pdf).

**Please note** that this repository, even with the appropriate libraries and packages installed, will not operate independently. Due to size limitations, the **original datasets** are not included. However, for those who are interested, the *original published dataset* as well as the code are available [elsewhere](https://doi.org/10.5281/zenodo.8229544).

At the moment this README is still under construction, more details to follow.


## Project Organization
```

┌── bulk_paired_end_mapping.py                      : from fastq to read count for paired-end samples
|── bulk_single_end_mapping.py                      : from fastq to read count for single-end samples
├── images/                                         : contains example images used for explanations within the README
│   └── XXX.png
│   └── XXX.png
├── bulk_seq_genome_construction.py                 : constructs alternate versions of the genome accounting for reads originating from AAV sequences
├── figure_s5.r                                     : code for displaying results for differential expression analysis 
├── LICENSE.md                                      : license
└── README.md                                       : project description

```
<br>

## Modified excerpts (figure and methods) from O'Toole et al. 2023 relevant to this repository

<p align="center">
<img src="https://github.com/sean-otoole/bulk_RNA_seq_otoole_2023/blob/main/images/figure_s5.png" height= 700>
</p>

### Differential expression of marker genes in neurons labelled with different artificial promoters.
**(A)** Bulk RNA-sequencing data for Adamts2 expression in populations of L1, L2/3 and L4 cortical neurons infected with an AAV2/1-AP.Adamts2.1-eGFP vector, for high and low eGFP expressing populations separately. The low and high eGFP groups constitute the lower and upper thirds of the GFP fluorescence distribution (cells with no expression were excluded). Fold-change (FC) values were normalized to the average expression levels in the low eGFP group. Note all samples, each consisting of 1000 sorted cells, were collected in pairs, however, in some cases library preparation failed for one of the two samples accounting for unequal sample sizes in some cases (here and in other figure panels of this figure). Box plots mark the median and quartiles, whiskers extend to cover data up to ± 1.5 inter quartile range past the 75th and 25th quartiles respectively. *: p < 0.05, **: p < 0.01, ***: p < 0.001, n.s.: not significant; see Table S1 for all statistical information. 
**(B)** As in **A**, but for the AP.Agmat.1 promoter.  
**(C)** As in **A**, but for the AP.Baz1a.1 promoter. 
**(D-F)** As in **A-C**, but for Agmat expression. 
**(G-I)** As in **A-C**, but for Baz1a expression. 
___

### Methods
To examine the differences in expression between neurons infected with different artificial promoters, a custom genome was constructed that excluded all DNA segments used in the artificial promoter viruses. This was done to increase the accuracy of the read counts for the marker genes. A custom GTF file was used to include the viral genomes. As with the single-cell RNA-sequencing analysis, the custom genome was based on the Genome Reference Consortium Mouse Build 39 with a version 104 GTF file. Differences in sequencing method (paired vs. single-end) were taken into account for the trimming and mapping steps. Reads were trimmed with the Cutadapt package (trim sequence: CTGTCTCTTATA). Gene counts were calculated with HTSeq count. Libraries were corrected for batch effects related to sequencing runs. Count values were normalized to library depth as counts per million (cpm). Genes with less than 5 cpm on average, across datasets, were excluded. Statistics were performed on fold change values.
