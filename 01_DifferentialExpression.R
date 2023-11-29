# Load libraries
library(dplyr)
library(data.table)
library(edgeR)
library(limma)

# Provide arguments
args = commandArgs(trailingOnly = TRUE)

## Create directories (again, please "no me borres") ##
if (!file.exists("Results")) {
  dir.create("Results")
}

if (!dir.exists("Results/DifferentialExpression")) {
  dir.create("Results/DifferentialExpression")
}

## Tests ##
# tissue<-"Ovary"
# sex<-"Male"
# comparison<-"Young_Old"
# Define function
diffexpress <- function(tissue, sex, comparison) {
  
  # Load metadata and counts
  metadata <- fread("data/ProcessedData/Metadata.csv", sep = ",", header = TRUE)
  counts <- fread(paste("data/ProcessedData/GeneReads/GeneReads", tissue, ".csv", sep = ""), sep = ",", header = TRUE)
  
  # metadata <- fread("~/Aging/data/ProcessedData/Metadata.csv", sep = ",", header = TRUE)
  # counts <- fread(paste("~/Aging/data/ProcessedData/GeneReads/GeneReads", tissue, ".csv", sep = ""), sep = ",", header = TRUE)

  # Identify categories to be compared
  control <- gsub("_.+", "", comparison)
  condition <- gsub(".+_", "", comparison)
  
  # Filter samples by AGESTATE before analysis.
  metadata <- metadata[metadata$AGESTATE %in% c(control, condition), ]

  # Filter tissue
  # Get column names from 'counts', excluding the first name (possibly 'Name').
  counts_col_names <- colnames(counts)

  # Filter the 'metadata' dataframe by matching SAMPID.
  metadata <- metadata[metadata$SAMPID %in% counts_col_names, ]
  
  # Retrieve the SAMPID present in metadata.
  sampids_metadata <- metadata$SAMPID
  
  # Filter the columns of 'counts' by the SAMPID from metadata (excluding 'Name').
  counts <- counts[, c("Name", sampids_metadata[sampids_metadata %in% colnames(counts)]), with = F]
  counts <- as.data.frame(counts)
  
  # Guardar los nombres de fila en una variable separada
  rownames(counts) <- counts$Name
  counts$Name <- NULL

  # Check if there are at least 8 samples for 'control' and 'condition'.
  if (sum(metadata$AGESTATE == control) >= 8 & sum(metadata$AGESTATE == condition) >= 8) {
    cat("At least 8 samples are available for the control group and 8 samples for the condition group.\n")
    ## Adjusting for sex differences ##
    if (sex == "Adjusted") {
      ## If we have both sexes ##
      if(length(unique(metadata$GENDER)) > 1) {
        ## Create factors ##
        AGESTATE<-factor(metadata$AGESTATE, levels = c(control,condition))
        GENDER<-factor(metadata$GENDER, levels = c("Male","Female"))
        DTHHRDY<-factor(metadata$DTHHRDY)
        
        ## Create the DGEList object ##
        dge <- DGEList(counts = counts, group = AGESTATE)
        
        ## Normalize scale factors ## 
        dge <- calcNormFactors(dge)
        
        cutoff <- 1
        drop <- which(apply(cpm(dge), 1, max) < cutoff)
        d <- dge[-drop,] 
        dim(d) # number of genes left
        
        ## Create design matrix ##
        mm <- model.matrix(~ AGESTATE + GENDER)
        
        # Plot normalizing or not
        y <- voom(d, mm, plot = F)
        # tmp <- voom(dge, mm, plot = F)
        
        ## Fit de linear model ##
        fit <- lmFit(y, mm, SIMPLIFY = F)
        
        ## See coefficients ##
        head(coef(fit))
        
        ## Apply eBayes ##
        tmp <- eBayes(fit, trend = T, robust = T)
        
        ## Get the summary of the differential expression analysis ##
        summary <- summary(decideTests(tmp))
        
        ## Get the table with the FC and FDRs ##
        top_table <- topTable(tmp, coef = colnames(summary)[2], number = length(fit$Amean))
        
        ## Write the table ##
        write.table(top_table,paste("Results/DifferentialExpression/",sex,"_",comparison,"_",tissue,".txt",sep=""),sep="\t",quote=F)
      } 
      ## If we only have one sex (prostate, ovary...) ##
      if(length(unique(metadata$GENDER)) == 1){
        ## Create factors ##
        AGESTATE<-factor(metadata$AGESTATE,levels=c(control,condition))
        DTHHRDY<-factor(metadata$DTHHRDY)
        
        ## Create the DGEList object ##
        dge <- DGEList(counts = counts, group = AGESTATE)
        
        ## Normalize scale factors ## 
        dge <- calcNormFactors(dge)
        
        cutoff <- 1
        drop <- which(apply(cpm(dge), 1, max) < cutoff)
        d <- dge[-drop,] 
        dim(d) # number of genes left
        
        ## Create design matrix ##
        mm <- model.matrix(~ AGESTATE)
        
        # Plot normalizing or not
        y <- voom(d, mm, plot = F)
        # tmp <- voom(dge, mm, plot = F)
        
        ## Fit de linear model ##
        fit <- lmFit(y, mm, SIMPLIFY = F)
        
        ## See coefficients ##
        # coef(fit)
        
        ## Apply eBayes ##
        tmp <- eBayes(fit, trend = T, robust = T)
        
        ## Get the summary of the differential expression analysis ##
        summary <- summary(decideTests(tmp))
        
        ## Get the table with the FC and FDRs ##
        top_table <- topTable(tmp, coef = colnames(summary)[2], number = length(fit$Amean))
        
        ## Write the table ##
        # write.table(top.table,paste("/gpfs/scratch/bsc08/bsc08295/TFM/Aging/out/DifferentialExpression/",sex,"_",comparison,"_",tissue,".txt",sep=""),sep="\t",quote=F)
        write.table(top_table,paste("Results/DifferentialExpression/",sex,"_",comparison,"_",tissue,".txt",sep=""),sep="\t",quote=F)
      }
      ## Return the top_table file ##
      return(top_table)
    }
    ## Analyzing sexes separately ## 
    if (sex != "Adjusted"){
      ## Select the samples of the analyzed sex ##
      metadata <- metadata[metadata$GENDER == sex]
      ## If we have at least 8 "cases" and 9 "controls" conduct a differential expression analysis ## 
      if (sum(metadata$AGESTATE == control) >= 8 & sum(metadata$AGESTATE == condition) >= 8){
        sampids_metadata <- metadata$SAMPID
        
        # Filtrar las columnas de 'counts' por los SAMPID de metadata (excluyendo 'Name')
        counts <- counts[, sampids_metadata[sampids_metadata %in% colnames(counts)]]
        
        AGESTATE<-factor(metadata$AGESTATE,levels=c(control,condition))
        DTHHRDY<-factor(metadata$DTHHRDY)
        
        ## Create the DGEList object ##
        dge <- DGEList(counts = counts, group = AGESTATE)
        
        ## Normalize scale factors ## 
        dge <- calcNormFactors(dge)
        
        cutoff <- 1
        drop <- which(apply(cpm(dge), 1, max) < cutoff)
        d <- dge[-drop,] 
        dim(d) # number of genes left
        
        snames <- colnames(counts) # Sample names
        
        ## Create design matrix ##
        mm <- model.matrix(~ AGESTATE)
        
        # Plot normalizing or not
        y <- voom(d, mm, plot = F)
        # tmp <- voom(dge, mm, plot = F)
        
        ## Fit de linear model ##
        fit <- lmFit(y, mm, SIMPLIFY = F)
        
        ## See coefficients ##
        # coef(fit)
        
        ## Apply eBayes ##
        tmp <- eBayes(fit, trend = T, robust = T)
        
        ## Get the summary of the differential expression analysis ##
        summary <- summary(decideTests(tmp))
        
        ## Get the table with the FC and FDRs ##
        top_table <- topTable(tmp, coef = colnames(summary)[2], number = length(fit$Amean))
        
        ## Write the table ##
        # write.table(top.table,paste("/gpfs/scratch/bsc08/bsc08295/TFM/Aging/out/DifferentialExpression/",sex,"_",comparison,"_",tissue,".txt",sep=""),sep="\t",quote=F)
        write.table(top_table,paste("Results/DifferentialExpression/", sex,"_",comparison,"_",tissue,".txt",sep=""),sep="\t",quote=F)
        ## Return the top_table file ##
        return(top_table)
      }
      ## If we do not print a message ##
      if (sum(metadata$AGESTATE == control) < 8 || sum(metadata$AGESTATE == condition) < 8){
        cat("Insufficient number of samples.\n")
      }
    }
  }
  if (sum(metadata$AGESTATE == control) < 8 || sum(metadata$AGESTATE == condition) < 8){
    cat("Insufficient number of samples.\n")
  }
}

# # Run function
# inicio<-Sys.time()
# resultado<-diffexpress(args[1], args[2], args[3])
# fin<-Sys.time()
# print(paste("Finished in: ", fin-inicio,sep=""))


#### Test with a loop ####
fichero<-fread("myrun_01DiffExpress.txt",stringsAsFactors = F,sep="\t",header=F)

for(a in 1:length(fichero$V1)){
  # a<-1
  variables<-strsplit(gsub(".+R ","",fichero$V1[a])," ")[[1]]
  inicio<-Sys.time()
  resultado<-diffexpress(variables[1], variables[2], variables[3])
  fin<-Sys.time()
  print(paste(round((a/length(fichero$V1))*100,2),"%, Finished in: ", fin-inicio,sep=""))
}























