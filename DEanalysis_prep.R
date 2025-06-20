#set seed for reproducibility purposes 
set.seed(10)

#load (or install if not installed) required libraries  
BiocManager::install("biomaRt") 
library(biomaRt) 

install.packages("dplyr") 
library(dplyr) 

BiocManager::install("edgeR") 
library(edgeR) 

#load featureCounts matrix 
FCmatrix <- read.table("path/to/featurecounts_matrix.txt", header=TRUE, check.names=FALSE)

colnames(FCmatrix)

#matrix cleanup ----------------------------------------------------------------

#clean up sample column names 
#adjust according to which columns samples are located in 
sampleColumns<-7:30
sampleColumnNames<-colnames(FCmatrix)[sampleColumns]

cleannames <- sub(".*/STARoutput/([A-Za-z0-9]+)_.*", "\\1", sampleColumnNames)

colnames(FCmatrix)[sampleColumns] <- cleannames

FCmatrix[, 2:5] <- apply(FCmatrix[, 2:5], 2, function(x) sapply(strsplit(as.character(x), ";"), `[`, 1))

#save cleaned up matrix as new dataset 
saveRDS(FCmatrix, "FCmatrix_cleaned.rds")

#load cleaned up matrix
FCmatrix <- readRDS("FCmatrix_cleaned.rds")

#gene annotation----------------------------------------------------------------

annotationInfo <- FCmatrix %>%
  dplyr::select(1:5)

# connect to Ensembl Biomart service 
ensembl <- useEnsembl(biomart = "genes")

# all datasets
datasets <- listDatasets(ensembl)

#find species dataset 
#this script uses the rattus norvegicus dataset 
ratMart <- useEnsembl(
  biomart = "genes",
  dataset = "rnorvegicus_gene_ensembl"
)

#all attributes and filters in rat dataset
attr <- listAttributes(ratMart) 
filters <- listFilters(ratMart)

#get geneID and gene data name data
geneNames <- getBM(attributes = c('ensembl_gene_id','external_gene_name'),
                     filters = "ensembl_gene_id",
                     values = annotationInfo$Geneid, 
                     mart = ratMart)

#merge gene names with annotation files 
mergedAnnotation <- merge(annotationInfo, geneNames, by.x = "Geneid", by.y = "ensembl_gene_id", all.y = TRUE)

#save merged annotation file as rds 
saveRDS(mergedAnnotation, "mergedAnnotation.rds")

# load merged annotation file 
mergedAnnotation <- readRDS("mergedAnnotation.rds")

# merge matrix with annotation to get matrix labelled with gene names 
FCwithgenes <- merge(FCmatrix, mergedAnnotation, by.x = "Geneid", by.y = "Geneid", all.y = TRUE)

#save merged matrix as rds file
saveRDS(FCwithgenes, "FCwithgenes.rds")

# load merged matrix
FCwithgenes <- readRDS("FCwithgenes.rds")

#dge list construction ---------------------------------------------------------

#view column names of matrix to determine factor order 
colnames(FCmatrix)

#factor each sample into a group 
#group labels must match the number of samples and order of samples in FC matrix                          
groups <- factor(c("MW","FC","FW","MC",
                   "MW","FC","FW","MC",
                   "MW","FC","FW","MC",
                   "FC","FW","FW","MC",
                   "FW","MC","MW","FC",
                   "MC","FC","MW","MW"))

# create dataframe with samples placed into respective groups  
sampleNames <- colnames(FCmatrix[7:30])
justSamples <- data.frame(sample = sampleNames, group = groups)

#get only counts matrix 
countMatrix <- FCmatrix %>%
  dplyr::select(7:30)

#daata frame with gene information
geneInfo <- data.frame(GeneID = mergedAnnotation$Geneid,
                       Chr = mergedAnnotation$Chr,
                       Start = mergedAnnotation$Start,
                       End = mergedAnnotation$End)

#DGElist using countMatrix, groups, and geneInfo
# make sure dimensions of each parameter match accordingly 
dge <- DGEList(counts = countMatrix, lib.size = colSums(countMatrix),
        norm.factors = rep(1,ncol(countMatrix)), samples = NULL,
        group = groups, genes = geneInfo, remove.zeros = FALSE)

#save dge list as RDS file 
saveRDS(dge, "dgeInitial.rds") 

dge <- readRDS("dgeInitial.rds")

# data normalization -----------------------------------------------------------

# calculates normalization factors using TMM normalization
dgeNorm <- calcNormFactors(dge)

#get normalization factors - values should be close to 1!!
normFactors <- dgeNorm$samples$norm.factors

#save dge after normalization to make it easier to retrieve later
saveRDS(dgeNorm, "dgeNormalized.rds")

#load dge
dge <- readRDS("dgeNormalized.rds")

# dispersion -------------------------------------------------------------------

# building design matrix 
design <- model.matrix (~ 0 + groups)
colnames(design) <- levels(groups)

#dispersion estimate
dgeDispersion <- estimateDisp(dge,design)

#prints summary stats of dispersion 
summary(dgeDispersion$common.dispersion)
summary(dgeDispersion$tagwise.dispersion)

#diagnostic plot - visual for dispersion 
plotBCV(dgeDispersion)

saveRDS(dgeDispersion, "dgeDispersion.rds")

#fit GLM
dgeFit <- glmFit(dgeDispersion,design)

#save fit dge model 
saveRDS(dgeFit, "dgeFit.rds")



