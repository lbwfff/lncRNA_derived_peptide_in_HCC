# A proteogenomic approach to discover novel lncRNA-derived peptides

This Repositories was used to include codes used in the project “A practical proteogenomic approach to discover novel lncRNA-derived peptides and their potential clinical utility in hepatocellular carcinoma” (unpublished).

Most of the code in this project is implemented by existing R packages (including Cran and bioconductor), and some of the code is implemented by a customized toolkit called TNSMD, which can be found at TNSMD (https://github.com/lbwfff/TNSMD).

## STEP 0 Pre-processing of data 
s0_data_prep.R consists of data preprocessing, where we formatted the NONCODE annotations to be suitable for use as input annotations to the Ribo-seq analysis tool, and removed redundancies with the GENCODE annotations (by using the GffCompare tool).

## STEP1 Ribo-seq analysis
The Ribo-seq folder contains the command lines we used for the Ribo-seq analysis, including trim for reads, removal of rRNA, mapping, and prediction of activated ORFs.

In s2_build_peptide_index.R we constructed the indexes needed for proteomic analysis based on the results of Ribo-seq.

In s3_ORF_Characterization.R we perform a preliminary characterization of the obtained ORFs

## STEP2 Computational Proteomics Analysis
We used the [Fragpipe](https://github.com/Nesvilab/FragPipe) platform to perform computational proteomics analysis, and the final parameters used can be seen in the fragpipe folder, which includes both sample annotation information and run logs.

s4_MS_fdr_control.R is used to visualize the results of the FDR control of the [philosopher](https://github.com/Nesvilab/philosopher) and to count the number of mass spectrometry spectra obtained by different methods.

## STEP3 Estimation of protein expression
s5_get_MS_expression.R was used to make an estimate of protein expression levels, which was achieved through the use of [MSstats](https://bioconductor.org/packages/release/bioc/html/MSstats.html) and [MSstatsTMT](https://www.bioconductor.org/packages/release/bioc/html/MSstatsTMT.html)

## STEP4 peptide Characterization
In s6_peptide_Characterization.R we trace the origin of the peptides, in addition we build a machine learning model for exploring why a large number of Ribo-seq peptides are not detected in protein mass spectrometry.

## STEP5 biostatistical analysis
Before that we used s7_QC_for_protein_expression.R to charge the expression matrices obtained by reanalyzing the [Jiang et al 2019 dataset](https://www.nature.com/articles/s41586-019-0987-8) and the [Gao et al 2019 dataset](https://www.sciencedirect.com/science/article/pii/S0092867419310037?via%3Dihub), and we checked for the presence of batch effects.

In s8_biomark.R we compared the predictive performance of models using canonical protein expression and lncRNA-derived peptides expression, and both were combined as biomarkers. These models included cancer tissue differentiation, survival prognosis prediction, and recurrence prognosis.

In s9_biostatistic_for_peptide.R we counted the biostatistical significance of separate lncRNA-derived peptides in the two datasets, and in s10_TCGA.R we discuss the lncRNAs from which these lncRNA-derived peptides originated in the TCGA-LIHC dataset for biostatistical significance.

In s11_for_3peptide.R we visualized three representative lncRNA-derived peptides.



