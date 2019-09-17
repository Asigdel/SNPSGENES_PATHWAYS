### August 4 ###################
## First STEP #################
## From SNPs to GENES 
rm(list = ls())
setwd("~/Documents/Documents - Anil’s MacBook Air/CODES_PLOTS/PLOTS_MANHATTAN_GENESET/Geneset_Folder/")

# biomaRT is a package from a bioconductor in R
# This is a collection of biomart databases 
library(biomaRt)
listMarts()

# we are using ENSEMBL_... database and btaurus_gene dataset by using the function usemart
# We are trying to query the Ensembl bioMart webservices
# useMart function queries the database Ensembl


# Selecting bioMart database and dataset in one query 

genome = useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                 dataset = "btaurus_gene_ensembl", 
                 host = "www.ensembl.org")


# How to build a biomaRt query?
# getBM is the main query function in biomaRt and we are calling it with the attribute names
# We are retrieving gene_symbols and chromosomal_ids.

gene = getBM(c("ensembl_gene_id", "start_position", "end_position", "chromosome_name"), mart = genome)
str(gene)

ndata = read.table("SNP_list_60619.txt", header = T) 
dim(ndata)
names(ndata)[1] = "SNP"; names(ndata)[2] = "Chromosome"; names(ndata)[5] = "Location"

plusi = 15000; plusf = 15000           #upstream (promoter) and downstream (regulation) regions
Name = character()
Chromosome = numeric()
Location = numeric()
Gene = numeric()
m = 1

for(k in 1:29)
{
  
  SNP = subset(ndata, ndata$Chromosome == k)
  genes = subset(gene, gene$chromosome_name == k)
  
  for(i in 1:length(SNP$Chromosome))
  {
    for(j in 1:length(genes$chromosome_name))
    {
      
      if(genes$start_position[j] <= (plusi + SNP$Location[i]) & (SNP$Location[i] - plusf) <= genes$end_position[j])
      {
        Name[m] = as.character(SNP$SNP[i])
        Chromosome[m] = SNP$Chromosome[i]
        Location[m] = SNP$Location[i]
        Gene[m] = genes$ensembl_gene_id[j]
        m = m + 1
      }
    }
  }
}

SNPtoGENES = data.frame(Name,Chromosome,Location,Gene)
save(SNPtoGENES, file = "FINAL_SNPtoGENES_15000and15000_August4.rda")


str(SNPtoGENES)
length(unique(SNPtoGENES$Gene))

# 
x <- load("FINAL_SNPtoGENES_15000and15000_August4.rda")
x
head(SNPtoGENES)

## 1 SNP in the same location, same chromosome position can flag many genes !!!



##################################################
# STEP 2 in GENE SET ANALYSIS IN R
#################################################
## Gene Set Analysis
## FILE freqdata.count.after.clean
filter = read.table("freqdata.count.after.clean", header = F)
colnames(filter) = c("SNP", "Frequency", "Filter")
f = filter$Filter == 0; table(f)


## FILE snp_sol
data = read.table("snp_sol_slope_1", header = F, dec = "."); dim(data)
data = data[f,]; dim(data) ## Clean data
names(data) = c("Trait","Effect","SNP","Chromosome","Position",
                "SNPSolution", "Weight", "Variance"); str(data)


## Significant SNPs Top 1%
ndata = data; thr = 0.05       #for 1% -> 0.005 * 2 = 0.01 = 1%
summary(abs(ndata$SNPSolution))
filter = abs(ndata$SNPSolution) >= quantile(abs(ndata$SNPSolution), 1 - thr)
table(filter)
sndata = ndata[filter,]; dim(sndata)

## load from SNP to GENES
load(file = "FINAL_SNPtoGENES_15000and15000_August4.rda") ## FILE output Rcode SNP to Genes
SNP.GENES = SNPtoGENES; str(SNP.GENES)
dim(SNP.GENES)

## Merge SNPtoGENES with SNP.Solution
fdata = merge(ndata, SNP.GENES, 
              by.x = c("Chromosome", "Position"),
              by.y = c("Chromosome", "Location"), sort = F); dim(fdata) 

sdata = merge(sndata, SNP.GENES, 
              by.x = c("Chromosome", "Position"),
              by.y = c("Chromosome", "Location"), sort = F); str(sdata); dim(sdata)

## Prepare data for Gene Set Analysis

total.genes = unique(fdata$Gene);length(total.genes)
write.table(total.genes,file="total.genes.slope1.txt",sep=" ",append=F, row.names=T,col.names=T)
sig.genes = unique(sdata$Gene);length(sig.genes)
write.table(sig.genes,file="sig.genes.slope1.txt",sep=" ",append=F, row.names=T,col.names=T)

#
###################
## GO/KEGG Analysis
library(goseq)
library(org.Bt.eg.db)
library(GO.db)
library(KEGG.db)
library(GOSemSim)
#library(corrplot)

assayed.genes = array(total.genes)
de.genes = array(sig.genes)
gene.vector = as.integer(assayed.genes%in%de.genes)
names(gene.vector) = assayed.genes
table = table(gene.vector)

cat(paste("Significant Genes: ", table[2], 
          " Backgroung Genes:", table[1], "\n")) 

pwf = nullp(gene.vector, "bosTau4", "ensGene", plot.fit = FALSE)

## KEGG
KEGG.hiper <- goseq(pwf, "bosTau4", "ensGene", 
                    method = "Hypergeometric", test.cats = "KEGG", use_genes_without_cat = TRUE)

nKEGG.hiper = KEGG.hiper[(KEGG.hiper$numInCat <= 500 & KEGG.hiper$numInCat >= 5),] 
enriched.KEGG = nKEGG.hiper$category[nKEGG.hiper$over_represented_pvalue <= 0.05]
length(enriched.KEGG)

kegg = as.list(KEGGPATHID2NAME)
for(j in 1:length(enriched.KEGG)){
  for (i in 1:length(names(kegg))) 
  {
    if(names(kegg[i]) == enriched.KEGG[j]){
      
      cat("#############################################\n")	
      cat(paste("KEGG ID: ", enriched.KEGG[j], "\n"))
      cat(paste("KEGG Term Name: ", kegg[[i]], "\n"))
      pvalue = nKEGG.hiper[nKEGG.hiper$category == enriched.KEGG[j],]$over_represented_pvalue
      pvalue = round(pvalue,4)
      nGenes = nKEGG.hiper[nKEGG.hiper$category == enriched.KEGG[j],]$numInCat
      nDEG = nKEGG.hiper[nKEGG.hiper$category == enriched.KEGG[j],]$numDEInCat
      cat(paste("Total Genes:", nGenes, " Number DEG:", nDEG, " P-value:", pvalue))
      cat("\n")
    }
  }
}


