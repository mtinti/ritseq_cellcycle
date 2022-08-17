![title](wcar.png)
# RIT-seq analysis of cell cycle regulators in _T. brucei_  
## This repository contains the scripts used for the data analysis

    Wellcome Centre for Anti-Infectives Research
    School of Life Sciences, University of Dundee

Prepare a folder structure as illustrated below. 
Add the bowtie2 index files, gtf/gff annotation files and the fasta format genome in the genome/tb927 folder. 
Each file needs to be named tb927.[extension].

The fastq files go in experiment/data folder. Replace experiment with the name of the sample to analyse:
CM03_B_LIB_S1, CM03_D_LG_S3, CM03_E_G1_S4, CM03_F_S_S5, CM03_G_G2M_S6, CM03_H_GG2M_S7. 

In the first 

```
project
│───ritseq.yml
│───analysis.sh
│───mylib
│   │ extract_barcode_def2.py
│
│───genome
│   │   
│   └───tb927 (bwtie2 index files)
│	    │ tb927.1.bt2
│	    │ tb927.2.bt2
│	    │ tb927.3.bt2
│	    │ tb927.4.bt2
│	    │ tb927.rev.1.bt2
│	    │ tb927.rev.1.bt2
│	    │ tb927.gff
│	    │ tb927.gtf
│
└───experiment
    │   
    └───data
        │ experiment_1.fastq.gz
        │ experiment_2.fastq.gz

```
