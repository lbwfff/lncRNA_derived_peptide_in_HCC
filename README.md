# A proteogenomic approach to discover novel lncRNA-derived peptides

This Repositories was used to include codes used in the project “A practical proteogenomic approach to discover novel lncRNA-derived peptides and their potential clinical utility in hepatocellular carcinoma” (unpublished).

Most of the code in this project is implemented by existing R packages (including Cran and bioconductor), and some of the code is implemented by a customized toolkit called TNSMD, which can be found at TNSMD (https://github.com/lbwfff/TNSMD).

##STEP 0 Pre-processing of data 
s0_data_prep.R consists of data preprocessing, where we formatted the NONCODE annotations to be suitable for use as input annotations to the Ribo-seq analysis tool, and removed redundancies with the GENCODE annotations (by using the GffCompare tool).

##Ribo-seq analysis
The Ribo-seq folder contains the command lines we used for the Ribo-seq analysis, including trim for reads, removal of rRNA, mapping, and prediction of activated ORFs.
