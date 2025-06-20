# RNA-seq Differential Expression Prep Pipeline  

This repository has an R script to preprocess RNA-seq gene count data from a **featureCounts** matrix. This prepares it for downstream **Differential Expression** analysis using **edgeR**. The outputs can then be used for direct use in DE testing, visualization and further analyses like **WGCNA**. 

## What does it do?  

* Cleans sample column names in **featureCounts matrix**
* Annotates gene IDs from Ensembl with gene symbols
* Constructs and normalizes a **DGEList** that can be analysed by **edgeR**
* Estimates dispersion and fits GLM
* Saves each major step as '.rds' object for modular use

## How do I use it?  

1. Data preparation  
Make sure the **featureCounts matrix** has:
* Gene annotation columns (typically columns 1 to 6, adjust column numbers in script accordingly)
* Count values in the remaining columns (each column for each sample)  
