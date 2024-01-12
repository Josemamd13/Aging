# Load necessary libraries
library(tidyverse)
library(data.table)
library(pvclust)
library(stringr)
library(gplots)
library(pheatmap)
library(dendextend)

#### GSEA Analysis ####

## Create the needed directories ##
# Check and create "Plots" directory
if (!"Plots" %in% list.files("Results")) {
  dir.create("Results/Plots")
}
# Check and create "Reactome" directory inside "Plots" if it doesn't exist
if (!"Reactome" %in% list.files("Results/Plots")) {
  dir.create("Results/Plots/Reactome")
}
# Check and create "MetabolicTask" directory inside "Plots" if it doesn't exist
if (!"MetabolicTask" %in% list.files("Results/Plots")) {
  dir.create("Results/Plots/MetabolicTask")
}

# Check and create "VisualizationObjects" directory
if (!"VisualizationObjects" %in% list.files("Results")) {
  dir.create("Results/VisualizationObjects")
}
# Check and create "Reactome" directory inside "VisualizationObjects" if it doesn't exist
if (!"Reactome" %in% list.files("Results/VisualizationObjects")) {
  dir.create("Results/VisualizationObjects/Reactome")
}
# Check and create "MetabolicTask" directory inside "VisualizationObjects" if it doesn't exist
if (!"MetabolicTask" %in% list.files("Results/VisualizationObjects")) {
  dir.create("Results/VisualizationObjects/MetabolicTask")
}


## Define a function that will take the list of desired tissues (all the tissues or the tissues common to both sexes) 
## in the selected comparison (eg, Young vs. Old) and sex and will cluster them
## based on the NES of their Reactome parents, ordered by Reactome category ##
generatetablesandplotdendrogram <- function(sex, tissues, comparison, database) {
  
  ini <- Sys.time()
  print("Getting the tables for the heatmap...")
  
  if (database == "Reactome") {
    
    ## Paths ##
    main_directory <- "Results/GSEA"
    
    # Manually curated data
    MyReactomePathwaysCurated <- fread("Resources/Heatmap/MyReactomePathwaysCurated.tsv",
                                       sep = "\t",
                                       header = TRUE) %>%
      dplyr::arrange(Parent)
    
    
    #########
    ## NES ##
    #########
    
    # Get directories matching the pattern "Adjusted_Young_Old_" for tissues of interest
    selected_directories <- list()
    selected_tissues <- character(0)
    
    for (tissue in tissues) {
      dir_pattern <- paste0(sex, "_", comparison, tissue, "_filtered")
      tissue_dir <- list.files(main_directory, pattern = dir_pattern, full.names = TRUE)
      if (length(tissue_dir) > 0) {
        selected_directories <- c(selected_directories, tissue_dir)
        selected_tissues <- c(selected_tissues, rep(tissue, length(tissue_dir)))
      } else {
        selected_tissues <- c(selected_tissues, tissue) # Add the tissue even if not found to maintain correspondence
      }
    }
    
    # GSEA_report.tsv files in selected directories
    report_files <- character(0)
    
    for (dir in selected_directories) {
      # Report directory
      report_dir <- file.path(dir, paste0(basename(dir), "_GSEA_report.tsv"))
      
      # Check if report file exists in the directory
      if (file.exists(report_dir)) {
        report_files <- c(report_files, report_dir)
      }
    }
    
    # Read GSEA_report.tsv files
    report_data <- list()
    
    for (file in report_files) {
      # Read the file and store data in a list
      data <- read.table(file, header = TRUE, sep = "\t")
      report_data[[file]] <- data
    }
    
    # Create a data.table to store pathways (NES values)
    NES <- data.table(Pathways = character(0))
    
    # Process GSEA_report.tsv files for tissues of interest
    for (tissue in tissues) {
      dir_pattern <- paste0(sex,"_",comparison,"_", tissue, "_filtered")
      directories <- list.files(main_directory, pattern = dir_pattern, full.names = TRUE)
      
      # Iterate through found directories
      for (dir in directories) {
        # Report directory
        report_dir <- file.path(dir, paste0(basename(dir), "_GSEA_report.tsv"))
        
        # Check if report file exists in the directory
        if (file.exists(report_dir)) {
          # Load data from the TSV file into a data.table
          tsv_data <- fread(report_dir, header = TRUE, sep = "\t")
          
          # Filter and store relevant pathway names (REACTOME_)
          reactome_pathways <- tsv_data$NAME[grepl("^REACTOME_", tsv_data$NAME)]
          
          # Add relevant names to the pathways data.table
          NES <- rbindlist(list(NES, data.table(Pathways = reactome_pathways)), fill = TRUE)
          
          # Get unique names of relevant pathways
          unique_names <- unique(reactome_pathways)
          
          # Assign values to pathways in the pathways data.table based on NES values
          for (pathway in unique_names) {
            pathway_data <- tsv_data[NAME == pathway]
            nes_value <- pathway_data$NES
            NES[pathway == NES$Pathways, (tissue) := nes_value]
          }
        }
      }
    }
    
    # Filter unique pathways based on the "Pathways" column
    NES <- unique(NES, by = "Pathways")
    
    # Replace NA with 0 in all tissue columns
    for (tissue in tissues) {
      NES[is.na(NES[[tissue]]), (tissue) := 0]
    }
    
    # Select columns "Reactome Pathway" and "Parent"
    MyReactomePathwaysCurated <- MyReactomePathwaysCurated[, .(Parent, `Reactome Pathway`)]
    
    # Merge NES with MyReactomePathways based on the "Pathways" column
    NES <- merge(NES, MyReactomePathwaysCurated, by.x = "Pathways", by.y = "Reactome Pathway", all.x = TRUE)
    
    # Create a data.table grouped by Parent without computing any summary function
    NES <- NES[, .SD, by = Parent]
    
    # Sort the DataFrame by the Parent column in alphabetical order
    NES <- NES[order(Parent)]
    
    
    ############
    ## Binary ##
    ############
    
    # Create a data.table to store pathways
    Binary <- data.table(Pathways = character(0))
    
    # Iterate through tissues of interest and process GSEA_report.tsv files
    for (tissue in tissues) {
      dir_pattern <- paste0(sex, "_", comparison, "_", tissue, "_filtered")
      directories <- list.files(main_directory, pattern = dir_pattern, full.names = TRUE)
      
      # Iterate through found directories
      for (dir in directories) {
        # Report directory
        report_directory <- file.path(dir, paste0(basename(dir), "_GSEA_report.tsv"))
        
        # Check if the TSV file exists
        if (file.exists(report_directory)) {
          # Load data from the TSV file into a data.table
          tsv_data <- fread(report_directory, header = TRUE, sep = "\t")
          
          # Modify NES values based on conditions
          tsv_data[, NES := ifelse(`FDR q-val` > 0.05, 0, ifelse(NES < 0, -1, 1))]
          
          # Filter and store relevant names (REACTOME_) as pathways
          reactome_pathways <- tsv_data$NAME[grepl("^REACTOME_", tsv_data$NAME)]
          
          # Add relevant names to the pathways data.table
          Binary <- rbindlist(list(Binary, data.table(Pathways = reactome_pathways)), fill = TRUE)
          
          # Get unique names of relevant pathways
          unique_names <- unique(reactome_pathways)
          
          # Assign values to pathways in the pathways data.table based on conditions
          for (pathway in unique_names) {
            pathway_data <- tsv_data[NAME == pathway]
            nes_value <- pathway_data$NES
            fdr_value <- pathway_data$`FDR q-val`
            Binary[pathway == Binary$Pathways, (tissue) := ifelse(fdr_value > 0.05, 0, ifelse(nes_value < 0, -1, 1))]
          }
        }
      }
    }
    
    # Filter unique pathways based on the "Pathways" column
    Binary <- unique(Binary, by = "Pathways")
    
    # Replace NA with 0 in all tissue columns
    for (tissue in tissues) {
      Binary[is.na(Binary[[tissue]]), (tissue) := 0]
    }
    
    # Select columns "Reactome Pathway" and "Parent"
    MyReactomePathwaysCurated <- MyReactomePathwaysCurated[, .(Parent, `Reactome Pathway`)]
    
    # Merge Binary with MyReactomePathways based on the "Pathways" column
    Binary <- merge(Binary, MyReactomePathwaysCurated, by.x = "Pathways", by.y = "Reactome Pathway", all.x = TRUE)
    
    # Create a data.table grouped by Parent without computing any summary function
    Binary <- Binary[, .SD, by = Parent]
    
    # Sort the DataFrame by the Parent column in alphabetical order
    Binary <- Binary[order(Parent)]
    
    
    ###############
    ### PVClust ###
    ###############
    
    # Perform cluster analysis with pvclust
    NESCluster <- pvclust(as.matrix(NES[, -c(1, 2)]),
                          method.dist = "euclidean",
                          method.hclust = "ward.D2",
                          nboot = 1000)
    
    ## Return the generated files ##
    resultados <- list("Binary" = Binary, "NES" = NES, "NESCluster" = NESCluster)
  } 
  
  else if (database == "MetabolicTask") {
    
    ## Paths ##
    main_directory <- "Results/GSEA"
    
    #########
    ## NES ##
    #########
    
    # Get directories matching the pattern "Adjusted_Young_Old_" for tissues of interest
    selected_directories <- list()
    selected_tissues <- character(0)
    
    for (tissue in tissues) {
      dir_pattern <- paste0(sex, "_", comparison, tissue, "_filtered")
      tissue_dir <- list.files(main_directory, pattern = dir_pattern, full.names = TRUE)
      if (length(tissue_dir) > 0) {
        selected_directories <- c(selected_directories, tissue_dir)
        selected_tissues <- c(selected_tissues, rep(tissue, length(tissue_dir)))
      } else {
        selected_tissues <- c(selected_tissues, tissue) # Add the tissue even if not found to maintain correspondence
      }
    }
    
    # GSEA_report.tsv files in selected directories
    report_files <- character(0)
    
    for (dir in selected_directories) {
      # Report directory
      report_dir <- file.path(dir, paste0(basename(dir), "_GSEA_report.tsv"))
      
      # Check if report file exists in the directory
      if (file.exists(report_dir)) {
        report_files <- c(report_files, report_dir)
      }
    }
    
    # Read GSEA_report.tsv files
    report_data <- list()
    
    for (file in report_files) {
      # Read the file and store data in a list
      data <- read.table(file, header = TRUE, sep = "\t")
      report_data[[file]] <- data
    }
    
    # Create a data.table to store pathways (NES values)
    NES <- data.table(Pathways = character(0))
    
    # Process GSEA_report.tsv files for tissues of interest
    for (tissue in tissues) {
      dir_pattern <- paste0(sex,"_",comparison,"_", tissue, "_filtered")
      directories <- list.files(main_directory, pattern = dir_pattern, full.names = TRUE)
      
      # Iterate through found directories
      for (dir in directories) {
        # Report directory
        report_dir <- file.path(dir, paste0(basename(dir), "_GSEA_report.tsv"))
        
        # Check if report file exists in the directory
        if (file.exists(report_dir)) {
          # Load data from the TSV file into a data.table
          tsv_data <- fread(report_dir, header = TRUE, sep = "\t")
          
          # Filter and store relevant pathway names (REACTOME_)
          reactome_pathways <- tsv_data$NAME[!grepl("^REACTOME_", tsv_data$NAME)]
          
          # Add relevant names to the pathways data.table
          NES <- rbindlist(list(NES, data.table(Pathways = reactome_pathways)), fill = TRUE)
          
          # Get unique names of relevant pathways
          unique_names <- unique(reactome_pathways)
          
          # Assign values to pathways in the pathways data.table based on NES values
          for (pathway in unique_names) {
            pathway_data <- tsv_data[NAME == pathway]
            nes_value <- pathway_data$NES  # Make sure 'NES' is the correct column name
            NES[Pathways == pathway, (tissue) := nes_value]  # Use 'Pathways' as the reference column
          }
        }
      }
    }
    
    # Filter unique pathways based on the "Pathways" column
    NES <- unique(NES, by = "Pathways")
    
    # Replace NA with 0 in all tissue columns
    for (tissue in tissues) {
      NES[is.na(NES[[tissue]]), (tissue) := 0]
    }
    
    
    ############
    ## Binary ##
    ############
    
    # Create a data.table to store pathways
    Binary <- data.table(Pathways = character(0))
    
    # Iterate through tissues of interest and process GSEA_report.tsv files
    for (tissue in tissues) {
      dir_pattern <- paste0(sex, "_", comparison, "_", tissue, "_filtered")
      directories <- list.files(main_directory, pattern = dir_pattern, full.names = TRUE)
      
      # Iterate through found directories
      for (dir in directories) {
        # Report directory
        report_directory <- file.path(dir, paste0(basename(dir), "_GSEA_report.tsv"))
        
        # Check if the TSV file exists
        if (file.exists(report_directory)) {
          # Load data from the TSV file into a data.table
          tsv_data <- fread(report_directory, header = TRUE, sep = "\t")
          
          # Modify NES values based on conditions
          tsv_data[, NES := ifelse(`FDR q-val` > 0.05, 0, ifelse(NES < 0, -1, 1))]
          
          # Filter and store relevant names (REACTOME_) as pathways
          reactome_pathways <- tsv_data$NAME[!grepl("^REACTOME_", tsv_data$NAME)]
          
          # Add relevant names to the pathways data.table
          Binary <- rbindlist(list(Binary, data.table(Pathways = reactome_pathways)), fill = TRUE)
          
          # Get unique names of relevant pathways
          unique_names <- unique(reactome_pathways)
          
          # Assign values to pathways in the pathways data.table based on conditions
          for (pathway in unique_names) {
            pathway_data <- tsv_data[NAME == pathway]
            nes_value <- pathway_data$NES
            fdr_value <- pathway_data$`FDR q-val`
            Binary[pathway == Binary$Pathways, (tissue) := ifelse(fdr_value > 0.05, 0, ifelse(nes_value < 0, -1, 1))]
          }
        }
      }
    }
    
    # Filter unique pathways based on the "Pathways" column
    Binary <- unique(Binary, by = "Pathways")
    
    # Replace NA with 0 in all tissue columns
    for (tissue in tissues) {
      Binary[is.na(Binary[[tissue]]), (tissue) := 0]
    }
    
    
    ###############
    ### PVClust ###
    ###############
    
    # Perform cluster analysis with pvclust
    NESCluster <- pvclust(as.matrix(NES[, -1]),
                          method.dist = "euclidean",
                          method.hclust = "ward.D2",
                          nboot = 1000)
    
    ## Return the generated files ##
    resultados <- list("Binary" = Binary, "NES" = NES, "NESCluster" = NESCluster)
  } 
  
  else {
    
    # Handle other cases or provide an error message if 'type' is not recognized
    stop("Invalid database argument. Please specify 'Reactome' or 'MetabolicTask'.")
  }
  
  fin <- Sys.time()
  print(paste0("Work done in: ", fin-ini))
  
  return(resultados)
}

## Define a function that will plot the dendrograms in the desired comparison ##
plotdendrograms <- function(comparison, database) {
  
  if (database == "Reactome") {
    
    pdf(file = paste0("Results/Plots/Reactome/Cluster_", comparison, ".pdf"))
    
    files <- list.files("Results/DifferentialExpression/")
    
    ## Select the samples for the heatmap (no matter if they are male or female specific) 
    alltissues <- gsub(".txt", "", gsub(paste(".+_", comparison, "_", sep = ""), "", files[intersect(grep(paste("Adjusted", "_", sep = ""), files), grep(comparison, files))]))
    
    ## Get the tissues for which, in the selected age comparison, we have enough samples in females and males ##
    femtissues <- gsub(".txt", "", gsub(paste(".+_", comparison, "_", sep = ""), "", files[intersect(grep("Female_", files), grep(comparison, files))]))
    maltissues <- gsub(".txt", "", gsub(paste(".+_", comparison, "_", sep = ""), "", files[intersect(grep("^Male_", files), grep(comparison, files))]))
    
    ## Select those tissues for which we have enough samples in females and males for another tissue ##
    commontissues <- intersect(femtissues, maltissues)
    
    #### Get the heatmap with all the tissues (no matter if they are sex-specific or not) ####
    AdjustedResultAllReactome <- generatetablesandplotdendrogram("Adjusted", alltissues, comparison, "Reactome")
    
    # Plot the dendrogram
    plot(AdjustedResultAllReactome$NESCluster, hang = -1, cex = 1, main = "Adjusted (All tissues)")
    pvrect(AdjustedResultAllReactome$NESCluster, max.only = FALSE)
    print("Adjusted (All tissues) finalized!")
    
    #### Get the heatmap for the common tissues ####
    AdjustedResultCommonReactome <- generatetablesandplotdendrogram("Adjusted", commontissues, comparison, "Reactome")
    
    # Plot the dendrogram
    plot(AdjustedResultCommonReactome$NESCluster, hang = -1, cex = 1, main = "Adjusted (Common tissues)")
    pvrect(AdjustedResultCommonReactome$NESCluster, max.only = FALSE)
    print("Adjusted (Common tissues) finalized!")
    
    #### Get the heatmap with male tissues ####
    MaleResultAllReactome <- generatetablesandplotdendrogram("Male", maltissues, comparison, "Reactome")
    
    # Plot the dendrogram
    plot(MaleResultAllReactome$NESCluster, hang = -1, cex = 1, main = "Male (Male tissues)")
    pvrect(MaleResultAllReactome$NESCluster, max.only = FALSE)
    print("Male (All tissues) finalized!")
    
    #### Get the heatmap for the common tissues ####
    MaleResultCommonReactome <- generatetablesandplotdendrogram("Male", commontissues, comparison, "Reactome")
    
    # Plot the dendrogram
    plot(MaleResultCommonReactome$NESCluster, hang = -1, cex = 1, main = "Male (Common tissues)")
    pvrect(MaleResultCommonReactome$NESCluster, max.only = FALSE)
    print("Male (Common tissues) finalized!")
    
    #### Get the heatmap with female tissues ####
    FemaleResultAllReactome <- generatetablesandplotdendrogram("Female", femtissues, comparison, "Reactome")
    
    # Plot the dendrogram
    plot(FemaleResultAllReactome$NESCluster, hang = -1, cex = 1, main = "Female (Female tissues)")
    pvrect(FemaleResultAllReactome$NESCluster, max.only = FALSE)
    print("Female (All tissues) finalized!")
    
    #### Get the heatmap for the common tissues ####
    FemaleResultCommonReactome <- generatetablesandplotdendrogram("Female", commontissues, comparison, "Reactome")
    
    # Plot the dendrogram
    plot(FemaleResultCommonReactome$NESCluster, hang = -1, cex = 1, main = "Female (Common tissues)")
    pvrect(FemaleResultCommonReactome$NESCluster, max.only = FALSE)
    print("Female (Common tissues) finalized!")
    
    dev.off()
    
    ## Put the results together ##
    ResultadosReactome <- list("AdjustedCommon" = AdjustedResultCommonReactome, "AdjustedAll" = AdjustedResultAllReactome, "FemaleCommon" = FemaleResultCommonReactome,
                       "FemaleAll" = FemaleResultAllReactome, "MaleCommon" = MaleResultCommonReactome, "MaleAll" = MaleResultAllReactome)
    
    ## Return the tables ##
    return(ResultadosReactome)
    
  } else if (database == "MetabolicTask") {
    
    pdf(file = paste0("Results/Plots/MetabolicTask/Cluster_", comparison, ".pdf"))
    
    files <- list.files("Results/DifferentialExpression/")
    
    ## Select the samples for the heatmap (no matter if they are male or female specific) 
    alltissues <- gsub(".txt", "", gsub(paste(".+_", comparison, "_", sep = ""), "", files[intersect(grep(paste("Adjusted", "_", sep = ""), files), grep(comparison, files))]))
    
    ## Get the tissues for which, in the selected age comparison, we have enough samples in females and males ##
    femtissues <- gsub(".txt", "", gsub(paste(".+_", comparison, "_", sep = ""), "", files[intersect(grep("Female_", files), grep(comparison, files))]))
    maltissues <- gsub(".txt", "", gsub(paste(".+_", comparison, "_", sep = ""), "", files[intersect(grep("^Male_", files), grep(comparison, files))]))
    
    ## Select those tissues for which we have enough samples in females and males for another tissue ##
    commontissues <- intersect(femtissues, maltissues)
    
    #### Get the heatmap with all the tissues (no matter if they are sex-specific or not) ####
    AdjustedResultAllMT <- generatetablesandplotdendrogram("Adjusted", alltissues, comparison, "MetabolicTask")
    
    # Plot the dendrogram
    plot(AdjustedResultAllMT$NESCluster, hang = -1, cex = 1, main = "Adjusted (All tissues)")
    pvrect(AdjustedResultAllMT$NESCluster, max.only = FALSE)
    print("Adjusted (All tissues) finalized!")
    
    #### Get the heatmap for the common tissues ####
    AdjustedResultCommonMT <- generatetablesandplotdendrogram("Adjusted", commontissues, comparison, "MetabolicTask")
    
    # Plot the dendrogram
    plot(AdjustedResultCommonMT$NESCluster, hang = -1, cex = 1, main = "Adjusted (Common tissues)")
    pvrect(AdjustedResultCommonMT$NESCluster, max.only = FALSE)
    print("Adjusted (Common tissues) finalized!")
    
    #### Get the heatmap with male tissues ####
    MaleResultAllMT <- generatetablesandplotdendrogram("Male", maltissues, comparison, "MetabolicTask")
    
    # Plot the dendrogram
    plot(MaleResultAllMT$NESCluster, hang = -1, cex = 1, main = "Male (Male tissues)")
    pvrect(MaleResultAllMT$NESCluster, max.only = FALSE)
    print("Male (All tissues) finalized!")
    
    #### Get the heatmap for the common tissues ####
    MaleResultCommonMT <- generatetablesandplotdendrogram("Male", commontissues, comparison, "MetabolicTask")
    
    # Plot the dendrogram
    plot(MaleResultCommonMT$NESCluster, hang = -1, cex = 1, main = "Male (Common tissues)")
    pvrect(MaleResultCommonMT$NESCluster, max.only = FALSE)
    print("Male (Common tissues) finalized!")
    
    #### Get the heatmap with female tissues ####
    FemaleResultAllMT <- generatetablesandplotdendrogram("Female", femtissues, comparison, "MetabolicTask")
    
    # Plot the dendrogram
    plot(FemaleResultAllMT$NESCluster, hang = -1, cex = 1, main = "Female (Female tissues)")
    pvrect(FemaleResultAllMT$NESCluster, max.only = FALSE)
    print("Female (All tissues) finalized!")
    
    #### Get the heatmap for the common tissues ####
    FemaleResultCommonMT <- generatetablesandplotdendrogram("Female", commontissues, comparison, "MetabolicTask")
    
    # Plot the dendrogram
    plot(FemaleResultCommonMT$NESCluster, hang = -1, cex = 1, main = "Female (Common tissues)")
    pvrect(FemaleResultCommonMT$NESCluster, max.only = FALSE)
    print("Female (Common tissues) finalized!")
    
    dev.off()
    
    ## Put the results together ##
    ResultadosMT <- list("AdjustedCommon" = AdjustedResultCommonMT, "AdjustedAll" = AdjustedResultAllMT, "FemaleCommon" = FemaleResultCommonMT,
                               "FemaleAll" = FemaleResultAllMT, "MaleCommon" = MaleResultCommonMT, "MaleAll" = MaleResultAllMT)
    
    ## Return the tables ##
    return(ResultadosMT)
  }
}


## Define a function that will plot the heatmaps with the selected ordering of the diseases ##
## Legend ##
pdf(file = "Results/Plots/Reactome/Legend.pdf")

# Loading data
color_data <- read.csv2("Resources/Heatmap/ReactomeCategoryColors.tsv", stringsAsFactors = FALSE, sep = "\t")

# Creating a vector of colors and assigning names
reaccolors <- color_data[, 2]
names(reaccolors) <- color_data[, 1]

# Combining category names into a single vector
all_names <- names(reaccolors)

# Open a new graphical window
plot.new()

# Set up the plot region to occupy the entire window
plot.window(xlim = c(0, 1), ylim = c(0, 1))

# Create a centered and large legend
legend(
  "center",
  legend = all_names,
  fill = reaccolors[all_names],
  title = "Function Legend",
  cex = 0.85,  # Larger size
  bty = "n",
  ncol = 1,  # Single column
  horiz = FALSE  # Vertical orientation
)

dev.off()

## Heatmap ##
plotheathmaps <- function(binarymatrix, matrixlabels, outputname, database) {
  
  if (database == "Reactome") {
    
    pdf(file = outputname)
    
    ### Heatmaps ###
    
    # Create a color palette based on the codes from the file
    color_data$Color <- as.character(color_data$Color)
    color_palette <- color_data$Color
    
    # Map the colors according to the category of each Pathway
    colorpathways <- color_data$Color[match(binarymatrix$Parent, color_data$Category)]
    
    
    thecluster <- binarymatrix %>%
      dplyr::select(c("Parent", "Pathways", matrixlabels))
    
    heatmap.2(as.matrix(thecluster[, -c(1, 2)]),
              key = FALSE,
              col = colorpanel(100, "royalblue", "#FFFFFF", "firebrick"),
              dendrogram = "none",
              Rowv = FALSE,
              labRow = "",
              labCol = "",
              cexCol = 0.85,
              scale = "none",
              margins = c(0.5, 0.5),
              trace = "none",
              RowSideColors = colorpathways,
              lwid = c(0.2, 5),
              Colv = FALSE)
    
    dev.off()  
    
  }
  
  else if (database == "MetabolicTask") {
    
    pdf(file = outputname)

    thecluster <- binarymatrix %>%
      dplyr::select(c("Pathways", matrixlabels))
    
    thecluster <- thecluster[order(thecluster$Pathways), ]
    
    heatmap.2(as.matrix(thecluster[, -1]),
              key = FALSE,
              col = colorpanel(100, "royalblue", "#FFFFFF", "firebrick"),
              dendrogram = "none",
              Rowv = FALSE,
              labRow = thecluster$Pathways,
              labCol = "",
              cexRow = 0.2,
              scale = "none",
              margins = c(10, 15),
              trace = "none",
              lwid = c(0.1, 3),
              Colv = FALSE)
    
    dev.off()  
  }
}


#### Age comparison: Young vs. Old (Reactome) ####
## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
## Get the tables for one age comparison ##
HeatmapResultsReactome <- plotdendrograms("Young_Old", "Reactome")
save(HeatmapResultsReactome, file = "Results/VisualizationObjects/Reactome/Young_Old.RData")

## Plot the heatmap for a specific age comparison and sex ##
for(a in 1:length(names(HeatmapResultsReactome))) {plotheathmaps(HeatmapResultsReactome[[a]]$Binary, labels(HeatmapResultsReactome[[a]]$NESCluster$hclust), paste0("Results/Plots/Reactome/Heatmap_", names(HeatmapResultsReactome)[a], "_Young_Old.pdf"), "Reactome")}

## Create the tanglegram (with the common tissues) ##
firstclusterdend <- as.dendrogram(HeatmapResultsReactome$MaleCommon$NESCluster)
secondclusterdend <- as.dendrogram(HeatmapResultsReactome$FemaleCommon$NESCluster)
diff_dend <- dend_diff(firstclusterdend, secondclusterdend)

## Plot the tanglegram ##
pdf(file = "Results/Plots/Reactome/Tanglegram_Young_Old.pdf")
tanglegram(diff_dend, main_left = "Males", main_right = "Females", rank_branches = TRUE, lab.cex = 0.75, lwd = 0.7, left_dendo_mar = c(3, 1, 3, 5.5), right_dendo_mar = c(3, 5.5, 3, 1), columns_width = c(5, 2, 5))
dev.off()

#### Age comparison: Young vs. Old (MT) ####
## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
## Get the tables for one age comparison ##
HeatmapResultsMT <- plotdendrograms("Young_Old", "MetabolicTask")
save(HeatmapResultsMT, file = "Results/VisualizationObjects/MetabolicTask/Young_Old.RData")

## Plot the heatmap for a specific age comparison and sex ##
for(a in 1:length(names(HeatmapResultsMT))) {plotheathmaps(HeatmapResultsMT[[a]]$Binary, labels(HeatmapResultsMT[[a]]$NESCluster$hclust), paste0("Results/Plots/MetabolicTask/Heatmap_", names(HeatmapResultsMT)[a], "_Young_Old.pdf"), "MetabolicTask")}

## Create the tanglegram (with the common tissues) ##
firstclusterdend <- as.dendrogram(HeatmapResultsMT$MaleCommon$NESCluster)
secondclusterdend <- as.dendrogram(HeatmapResultsMT$FemaleCommon$NESCluster)
diff_dend <- dend_diff(firstclusterdend, secondclusterdend)

## Plot the tanglegram ##
pdf(file = "Results/Plots/MetabolicTask/Tanglegram_Young_Old.pdf")
tanglegram(diff_dend, main_left = "Males", main_right = "Females", rank_branches = TRUE, lab.cex = 0.75, lwd = 0.7, left_dendo_mar = c(3, 1, 3, 5.5), right_dendo_mar = c(3, 5.5, 3, 1), columns_width = c(5, 2, 5))
dev.off()



#### Age comparison: Adult vs. Old (Reactome) ####
## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
## Get the tables for one age comparison ##
HeatmapResultsReactome <- plotdendrograms("Adult_Old", "Reactome")
save(HeatmapResultsReactome, file = "Results/VisualizationObjects/Reactome/Adult_Old.RData")

## Plot the heatmap for a specific age comparison and sex ##
for(a in 1:length(names(HeatmapResultsReactome))) {plotheathmaps(HeatmapResultsReactome[[a]]$Binary, labels(HeatmapResultsReactome[[a]]$NESCluster$hclust), paste0("Results/Plots/Reactome/Heatmap_", names(HeatmapResultsReactome)[a], "_Adult_Old.pdf"), "Reactome")}

## Create the tanglegram (with the common tissues) ##
firstclusterdend <- as.dendrogram(HeatmapResultsReactome$MaleCommon$NESCluster)
secondclusterdend <- as.dendrogram(HeatmapResultsReactome$FemaleCommon$NESCluster)
diff_dend <- dend_diff(firstclusterdend, secondclusterdend)

## Plot the tanglegram ##
pdf(file = "Results/Plots/Reactome/Tanglegram_Adult_Old.pdf")
tanglegram(diff_dend, main_left = "Males", main_right = "Females", rank_branches = TRUE, lab.cex = 0.75, lwd = 0.7, left_dendo_mar = c(3, 1, 3, 5.5), right_dendo_mar = c(3, 5.5, 3, 1), columns_width = c(5, 2, 5))
dev.off()

#### Age comparison: Adult vs. Old (MT) ####
## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
## Get the tables for one age comparison ##
HeatmapResultsMT <- plotdendrograms("Adult_Old", "MetabolicTask")
save(HeatmapResultsMT, file = "Results/VisualizationObjects/MetabolicTask/Adult_Old.RData")

## Plot the heatmap for a specific age comparison and sex ##
for(a in 1:length(names(HeatmapResultsMT))) {plotheathmaps(HeatmapResultsMT[[a]]$Binary, labels(HeatmapResultsMT[[a]]$NESCluster$hclust), paste0("Results/Plots/MetabolicTask/Heatmap_", names(HeatmapResultsMT)[a], "_Adult_Old.pdf"), "MetabolicTask")}

## Create the tanglegram (with the common tissues) ##
firstclusterdend <- as.dendrogram(HeatmapResultsMT$MaleCommon$NESCluster)
secondclusterdend <- as.dendrogram(HeatmapResultsMT$FemaleCommon$NESCluster)
diff_dend <- dend_diff(firstclusterdend, secondclusterdend)

## Plot the tanglegram ##
pdf(file = "Results/Plots/MetabolicTask/Tanglegram_Adult_Old.pdf")
tanglegram(diff_dend, main_left = "Males", main_right = "Females", rank_branches = TRUE, lab.cex = 0.75, lwd = 0.7, left_dendo_mar = c(3, 1, 3, 5.5), right_dendo_mar = c(3, 5.5, 3, 1), columns_width = c(5, 2, 5))
dev.off()



#### Age comparison: Young vs. Adult (Reactome) ####
## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
## Get the tables for one age comparison ##
HeatmapResultsReactome <- plotdendrograms("Young_Adult", "Reactome")
save(HeatmapResultsReactome, file = "Results/VisualizationObjects/Reactome/Young_Adult.RData")

## Plot the heatmap for a specific age comparison and sex ##
for(a in 1:length(names(HeatmapResultsReactome))) {plotheathmaps(HeatmapResultsReactome[[a]]$Binary, labels(HeatmapResultsReactome[[a]]$NESCluster$hclust), paste0("Results/Plots/Reactome/Heatmap_", names(HeatmapResultsReactome)[a], "_Young_Adult.pdf"), "Reactome")}

## Create the tanglegram (with the common tissues) ##
firstclusterdend <- as.dendrogram(HeatmapResultsReactome$MaleCommon$NESCluster)
secondclusterdend <- as.dendrogram(HeatmapResultsReactome$FemaleCommon$NESCluster)
diff_dend <- dend_diff(firstclusterdend, secondclusterdend)

## Plot the tanglegram ##
pdf(file = "Results/Plots/Reactome/Tanglegram_Young_Adult.pdf")
tanglegram(diff_dend, main_left = "Males", main_right = "Females", rank_branches = TRUE, lab.cex = 0.75, lwd = 0.7, left_dendo_mar = c(3, 1, 3, 5.5), right_dendo_mar = c(3, 5.5, 3, 1), columns_width = c(5, 2, 5))
dev.off()

#### Age comparison: Young vs. Adult (MT) ####
## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
## Get the tables for one age comparison ##
HeatmapResultsMT <- plotdendrograms("Young_Adult", "MetabolicTask")
save(HeatmapResultsMT, file = "Results/VisualizationObjects/MetabolicTask/Young_Adult.RData")

## Plot the heatmap for a specific age comparison and sex ##
for(a in 1:length(names(HeatmapResultsMT))) {plotheathmaps(HeatmapResultsMT[[a]]$Binary, labels(HeatmapResultsMT[[a]]$NESCluster$hclust), paste0("Results/Plots/MetabolicTask/Heatmap_", names(HeatmapResultsMT)[a], "_Young_Adult.pdf"), "MetabolicTask")}

## Create the tanglegram (with the common tissues) ##
firstclusterdend <- as.dendrogram(HeatmapResultsMT$MaleCommon$NESCluster)
secondclusterdend <- as.dendrogram(HeatmapResultsMT$FemaleCommon$NESCluster)
diff_dend <- dend_diff(firstclusterdend, secondclusterdend)

## Plot the tanglegram ##
pdf(file = "Results/Plots/MetabolicTask/Tanglegram_Young_Adult.pdf")
tanglegram(diff_dend, main_left = "Males", main_right = "Females", rank_branches = TRUE, lab.cex = 0.75, lwd = 0.7, left_dendo_mar = c(3, 1, 3, 5.5), right_dendo_mar = c(3, 5.5, 3, 1), columns_width = c(5, 2, 5))
dev.off()
