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
