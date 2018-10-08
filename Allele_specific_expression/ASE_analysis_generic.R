########################################################################################################
########################################################################################################
#
# Allele Specific Expression of color genes in Petunia 
#
########################################################################################################
########################################################################################################
# 
# using petal (limb) tissue from stage 4
#
########################################################################################################
#
# if you use this script for your own ASE analysis please cite:
#
# Oleg Mayba and Houston Gilbert (2014). MBASED: Package containing functions for ASE analysis using Meta-analysis
# Based Allele-Specific Expression Detection. R package version 1.12.0.
#
# and acknowledge Dr. Michel Moser in your Methods section
#
#################
library("MBASED")
citation("MBASED")

####INSTALLING the packages ########
#source("https://bioconductor.org/biocLite.R")
#biocLite("MBASED")
#biocLite("MBA")
#detach("package:DESeq2", unload=TRUE)
#install.packages("gtools")

#######loading the packages ##############
library(gtools)
library(ggbio)
library(ggplot2)
library(VariantAnnotation)
library(rtracklayer)
library(reshape2)
library(GenomicRanges)
library(GenomicFeatures)


#set a seed:  (needed for the iterations to obtain same results)
set.seed(23)

#function to extract results provided by Oleg

summarizeASEResults_1s <- function(MBASEDOutput) {
  geneOutputDF <- data.frame(
    majorAlleleFrequency=assays(MBASEDOutput)$majorAlleleFrequency[,1],
    pValueASE=assays(MBASEDOutput)$pValueASE[,1],  
    pValueHeterogeneity=assays(MBASEDOutput)$pValueHeterogeneity[,1]
  )   
  lociOutputGR <- rowRanges(metadata(MBASEDOutput)$locusSpecificResults)
  lociOutputGR$allele1IsMajor <- assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor[,1]
  lociOutputGR$MAF <- assays(metadata(MBASEDOutput)$locusSpecificResults)$MAF[,1]
  lociOutputList <- split(lociOutputGR, factor(lociOutputGR$aseID, levels=unique(lociOutputGR$aseID))) 
  return(
    list(
      geneOutput=geneOutputDF,
      locusOutput=lociOutputList
    )
  )
}

#############################################
### GENERAL INFO ############################
############################################

#working with a list of the candidate genes
#MAF = major allele frequency
#make sure intronic regions get extracted as depth and MAF varies between exons and introns

################################################
###run MBASED on ASEreadCounter output     ####
###                                         ####      
###         DATA:                           ####
###                                         ####
### 2015 F1 AX x IN limb stage 4 dataset    ####
################################################

#read in data
#replace with exon data to make sure for consistency in ASE calls

dir <-"/media/mmoser/data1/Publications/COLOUR/Pseudogenization_and_resurrection_CB/Allele_specific_expression" 

setwd(dir)


header = c("gene", "first", "last", "scaffold", "pos", "baseAx", "baseIn", "rep1Ax", "rep1In", "rep1tot", "rep2Ax", "rep2In", "rep2tot", "rep3Ax", "rep3In", "rep3tot")

axin15limb <- read.table('F1axin.v3_162.s.ASEcounts.tsv', header = FALSE)
names(axin15limb) <- header


#remove lines with identical allele
dim(axin15limb)
head(axin15limb)
axin15limb <- axin15limb[which(axin15limb$baseAx != axin15limb$baseIn),]



####### general dataset inspection: ###########3
#how many genes are present in the dataset?
length(unique(axin15limb$gene))
dim(axin15limb)

#[1] 15803

#how many SNPs are there for each of the present genes?
table(axin15limb$gene)

#create mean of each allelic count as new column to dataframe
axin15limb$axmean <- (axin15limb$rep1Ax + axin15limb$rep2Ax + axin15limb$rep3Ax)/3
axin15limb$inmean <- (axin15limb$rep1In + axin15limb$rep2In + axin15limb$rep3In)/3
axin15limb$gene <- as.character(axin15limb$gene)

#order according to gene_name
axin15limb_ord <- axin15limb[mixedorder(axin15limb$gene),]
dim(axin15limb_ord)
head(axin15limb_ord)

#look at SNPs from a specific gene:

#FLS
axin15limb_ord[which(axin15limb_ord$gene == "Peaxi162Scf00927g00035"),]

#load genes of interest (based on annotation 162v3)
phenyl_genes <- read.table("color_genes.csv", header = TRUE, sep = ",")
gene_list <- as.character(phenyl_genes$AXILLARISv3)
length(gene_list)
#extract all SNPs from the candidate genes
head(axin15limb[axin15limb$gene %in% gene_list,])
dim(axin15limb[axin15limb$gene %in% gene_list,])

#get all SNPs within the candidate genes from the scent list
candidate_loci <- axin15limb_ord[axin15limb_ord$gene %in% gene_list,]

#for which genes of the list are SNPs found?
present_candidategenes <- as.character(candidate_loci$gene)
#gives you summary of the uniq values and their count
present <- unique(present_candidategenes)


length(present); length(gene_list)
#10 of 10 genes have SNPs to measure ASE

#Genes present in the ASE dataset
phenyl_genes[ phenyl_genes$AXILLARISv3 %in% unique(present_candidategenes), ]
#Genes NOT present in the ASE dataset
phenyl_genes[ !phenyl_genes$AXILLARISv3 %in% unique(present_candidategenes), ]

################################
# write informative SNPs and their allele counts out for each of the candidate genes
#################################

#add gene name to list
output <- DataFrame()

for(i in 1:nrow(candidate_loci)) {
  row <- candidate_loci[i,]
  gene_id <- candidate_loci[i,]$gene
  gene_name <- phenyl_genes[which(phenyl_genes$AXILLARISv3 == gene_id),]$Name[1]
  #in case there are two names for the same gene: 
  #gene_name <- split(gene_name)
  print(gene_name)
  row$name <- gene_name
  output = rbind(output,row)
}

head(output, 30)

dim(candidate_loci)
dim(output)

#overwrite candidate_loci with more information
candidate_loci <- output

#filter out SNVs which have 0 coverage
candidate_loci <- candidate_loci[round(candidate_loci$axmean) + round(candidate_loci$inmean) != 0,]
dim(candidate_loci)


#########################
##  CREATE INPUT     ####
## allele counts     ####
#########################


#############################
## SET WHICH READS TO DO ASE WITH
################################


head(candidate_loci)


candidate_loci$scaffold <- as.character(candidate_loci$scaffold)
lapply(candidate_loci, class)


#create candidat_loci for only one gene
#candidate_loci <- candidate_loci[which(candidate_loci$gene == "Peaxi162Scf00050g00423"),]



############################
### ASE Read input data ####
############################

#define counts from which replicate to test
ca.count <- round(candidate_loci$axmean)
ce.count <- round(candidate_loci$inmean)

#check for missing data
atLeast1ReadSubv <- (ca.count+ce.count)>0
table(atLeast1ReadSubv)

ca.counts <- ca.count[atLeast1ReadSubv]
ce.counts <- ce.count[atLeast1ReadSubv]

#########################
##  CREATE Data frame ###
#########################

mySNVs <- GRanges(
  seqnames = candidate_loci$scaffold,
  ranges = IRanges(start = candidate_loci$pos, width = 1),
  aseID = as.character(candidate_loci$gene),
  allele1= candidate_loci$baseAx, 
  allele2= candidate_loci$baseIn
)
mySNVs


#create dataframe with this 
#how many snps per gene 

unique(present_candidategenes)

#list used not loci!
d= NULL
for (genei in unique(present_candidategenes)){
  print(genei)
  print(length(grep(genei, candidate_loci$gene)))
  snpcount <- length(grep(genei, candidate_loci$gene))
  d = rbind(d, data.frame(genei, snpcount))
}

d
#filter out zeros: 

d_good <- d[d$snpcount != 0,] 
d_good


#create list of names for each snp with format: GENE:##nth SNP in gene## 
ASEnames = NULL
for(i in 1:nrow(d_good)){
  for (e in 1:d[i,2]){
    padded <- sprintf("%02d",e)
    file_name<-paste(d[i,1],".",padded,sep='')
    #names1[e] = file_name    
    ASEnames = rbind(ASEnames, data.frame(file_name))
  }
}

ASEnames <- lapply(ASEnames, as.character)
ASEnames

names(mySNVs) <- ASEnames$file_name

names(mySNVs)

length(ca.count)

########
# parameter settings
# rho (overdispersion parameter)
# mu (mapping bias/ preexisting allelic bias)
#######

#set to default right now a bit lower than in the publication of MBASED (0.004)

rho = 0.002

mySample <- SummarizedExperiment(
  assays = list(
    lociAllele1Counts= matrix(ca.count, 
                              ncol = 1, 
                              dimnames=list(names(mySNVs), 'mySample')
    ), 
    lociAllele2Counts=matrix(ce.count, 
                             ncol = 1, 
                             dimnames= list(names(mySNVs), 'mySample')
    ), 
    lociCountsDispersions=matrix(
      rep(rho, length(ca.count)), 
      ncol= 1, 
      dimnames= list(names(mySNVs), 'mySample')
    )
  ), 
  rowRanges = mySNVs
)
traceback()
mySNVs$aseID
############################
############################
##  RUN MBASED ASE TEST   ##
## with 1 mio simulations ##
############################
############################

mySample

ASEtest <- runMBASED(ASESummarizedExperiment = mySample, 
                     isPhased = TRUE, 
                     numSim = 10^6, 
                     BPPARAM = SerialParam()
)



###########################
##  Inspect results     ####
############################

ASEtest
t <- summarizeASEResults_1s(ASEtest)
t
ASE_2015f1axinlimb<- t$geneOutput
ASE_2015f1axinlimb

#which genes show ASE with MAF more than 0.7 
ASE_2015f1axinlimb_candidates <- data.frame()

for(i in 1:nrow(ASE_2015f1axinlimb)) {
  row <- ASE_2015f1axinlimb[i,]
  gene_id <- as.character(rownames(ASE_2015f1axinlimb[i,]))
  gene_id
  gene_name <- phenyl_genes[which(phenyl_genes$AXILLARISv3 == gene_id),]$Name[1]
  gene_name
  row$name <- gene_name
  print(row)
  ASE_2015f1axinlimb_candidates = rbind(ASE_2015f1axinlimb_candidates , row)
}

#look at single genes
ASE_2015f1axinlimb_candidates[which(ASE_2015f1axinlimb_candidates$name == "HF1"),]

ASE_2015f1axinlimb_candidates <- as.data.frame(ASE_2015f1axinlimb_candidates)

ASE_2015f1axinlimb_candidates[ASE_2015f1axinlimb_candidates$name %in% c("AN2","HT1","HF1","HF2_1" ,"DFR_a","5GT"),]



#what are the real counts: ASE towards which side?
ASEmean_2015f1axinlimb_counts <- candidate_loci[which(candidate_loci$gene %in% rownames(ASE_2015f1axinlimb_candidates)),c(1,17:19)]

ASEmean_2015f1axinlimb_counts


#pro Gen sum up all SNP site counts
SNPsumup <- aggregate(cbind(axmean, inmean) ~ gene, data=ASEmean_2015f1axinlimb_counts, FUN=sum)

aggregate(cbind(axmean, inmean) ~ gene, data=ASEmean_2015f1axinlimb_counts, FUN=sum)
ASEmean_2015f1axinlimb_counts


#####
# get mean counts per site and its position
#####

head(candidate_loci)
snp_table_ASE <- candidate_loci[,c(1,4, 5, 6, 7, 17,18, 19)]
snp_table_ASE

SNPsumup$snp_count <- table(ASEmean_2015f1axinlimb_counts$gene)

AXIN_SNPsum <- DataFrame()

for(i in 1:nrow(SNPsumup)) {
  row <- SNPsumup[i,]
  gene_id <- SNPsumup[i,]$gene
  gene_name <- phenyl_genes[which(phenyl_genes$AXILLARISv3 == gene_id),]$Name[1]
  row$name <- gene_name
  AXIN_SNPsum = rbind(AXIN_SNPsum,row)
}

AXIN_SNPsum

