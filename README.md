# A proteogenomic approach to discover novel lncRNA-derived peptides

This Repositories was used to include codes used in the project “A practical proteogenomic approach to discover novel lncRNA-derived peptides and their potential clinical utility in hepatocellular carcinoma” (unpublished).

Most of the code in this project is implemented by existing R packages (including Cran and bioconductor), and some of the code is implemented by a customized toolkit called TNSMD, which can be found at TNSMD (https://github.com/lbwfff/TNSMD).

## STEP 0 Pre-processing of data 
s0_data_prep.R consists of data preprocessing, where we formatted the NONCODE annotations to be suitable for use as input annotations to the Ribo-seq analysis tool, and removed redundancies with the GENCODE annotations (by using the GffCompare tool).

## STEP1 Ribo-seq analysis
The Ribo-seq folder contains the command lines we used for the Ribo-seq analysis, including trim for reads, removal of rRNA, mapping, and prediction of activated ORFs.
In s2_ORF_Characterization.R we perform a preliminary characterization of the obtained ORFs
In s3_build_peptide_index.R we constructed the indexes needed for proteomic analysis based on the results of Ribo-seq, 

## 
