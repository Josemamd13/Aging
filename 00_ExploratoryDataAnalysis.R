# Load necessary libraries
library(data.table)
library(readxl)
library(dplyr)
library(rjson)
library(future)
library(future.apply)

  # Check and install ReactomeContentService4R package if not installed
  if (!require("ReactomeContentService4R", quietly = TRUE)) {
    BiocManager::install("ReactomeContentService4R")
    library(ReactomeContentService4R)
  }

# Provide arguments
args = commandArgs(trailingOnly = TRUE)


#### Prepare all the tables for the differential expression analysis ####
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
if(args[1] == "PrepareTables") {
  # Reading GTEx Annotation Data
  GTEx_Analysis_v8_Annotations_SampleAttributesDD <- read_excel("data/RawData/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx")
  
  # Filtering rows in GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
  # Removing rows where the 'SMTSD' column contains NA or is empty
  GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt <- fread("data/RawData/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                                                               header = TRUE, 
                                                               sep = "\t") %>%
    filter(!is.na(SMTSD) & SMTSD != "")
  
  # Extracting a shortened version of 'SAMPID' column
  GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt$SAMPID_SHORT <- sub("^([^\\-]+-[^\\-]+)-.*", "\\1",
                                                                          GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt$SAMPID)
  
  # Subsetting specific columns from the dataset
  GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt <- subset(GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt, 
                                                                select = c("SAMPID", "SAMPID_SHORT", "SMTS", "SMTSD"))
  
  # Converting the data table to a data frame
  GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt <- as.data.frame(GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt)
  
  
  # Reading and Preprocessing GTEx Phenotype Data
  # Reading data from an Excel file containing subject phenotype annotations
  GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD <- read_excel("data/RawData/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx")
  
  # Reading data from a tab-separated file GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
  # Filtering rows where the column 'SUBJID' is not NA or empty
  GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt <- fread("data/RawData/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", 
                                                                header = TRUE, 
                                                                sep = "\t") %>%
    filter(!is.na(SUBJID) & SUBJID != "") %>%
    as.data.frame()
  
  
  # Merging sample attributes ("SampleAttributesDS") and subject phenotypes ("SubjectPhenotypesDS") data 
  # for subjects into the "Metadata" object.
  Metadata <- merge(GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt,
                    GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt, 
                    by.x = "SAMPID_SHORT", 
                    by.y = "SUBJID", 
                    all.x = TRUE)
  
  # Sorting the merged data by SEX and AGE columns.
  Metadata <- Metadata %>%
    arrange(SEX, AGE)
  
  # Writing the merged data into a CSV file named "Metadata.csv" in the specified directory.
  write.csv(Metadata, file = "data/RawData/Metadata.csv",
            row.names = FALSE)
  
  
  # Adding new columns based on specified conditions.
  Metadata <- Metadata %>%
    # Creating AGESTATE column based on AGE conditions
    mutate(AGESTATE = case_when(
      AGE >= 20 & AGE <= 29 ~ "Young",
      AGE >= 60 & AGE <= 79 ~ "Old",
      TRUE ~ "Adult"
    )) %>%
    # Creating GENDER column based on SEX conditions
    mutate(GENDER = case_when(
      SEX == 1 ~ "Male",
      SEX == 2 ~ "Female",
      TRUE ~ "Unknown"
    )) %>%
    # Creating DEATH column based on DTHHRDY conditions
    mutate(DEATH = case_when(
      DTHHRDY == 0 ~ "Ventilator case",
      DTHHRDY == 1 ~ "Violent and fast death",
      DTHHRDY == 2 ~ "Fast death of natural case",
      DTHHRDY == 3 ~ "Intermediate",
      DTHHRDY == 4 ~ "Slow death"
    )) %>%
    # Selecting specific columns
    dplyr::select(SAMPID_SHORT, SAMPID, SMTS, SMTSD, SEX, GENDER, AGE, AGESTATE, DTHHRDY, DEATH)
  
  # # Writing processed data to a CSV file
  # write.csv(Metadata, file = "data/ProcessedData/Metadata.csv",
  #           row.names = FALSE)
  
  
  # Extracting unique SMTS values
  SMTS <- as.data.frame(unique(Metadata$SMTS))
  
  # Extracting unique SMTSD values
  SMTSD <- as.data.frame(unique(Metadata$SMTSD))
  
  # Creating a new table with grouped values
  SMTS_SMTSD <- Metadata %>%
    group_by(SMTS) %>%
    summarise(SMTSD = paste(unique(SMTSD), collapse = ", "))
  
  # Calculating the number of SMTSD for each SMTS entry
  SMTS_SMTSD <- SMTS_SMTSD %>%
    mutate("Number of SMTSD" = str_count(SMTSD, ",") + 1)
  
  cat("Number of SMTS: ", nrow(SMTS), "\n")
  cat("Number of SMTSD: ", nrow(SMTSD), "\n")
  
  # Extracting unique SMTS for males (MALES)
  SMTS_MALES <- Metadata %>%
    filter(GENDER == "Male") %>%
    distinct(SMTS)
  
  # Extracting unique SMTS for females (FEMALES)
  SMTS_FEMALES <- Metadata %>%
    filter(GENDER == "Female") %>%
    distinct(SMTS)
  
  # Common SMTS
  SMTS_COMMON <- intersect(SMTS_MALES, SMTS_FEMALES)
  
  # Identifying missing SMTS in males (MALES)
  SMTS_MISSING_MALES <- setdiff(SMTS_FEMALES$SMTS, SMTS_MALES$SMTS)
  
  # Identifying missing SMTS in females (FEMALES)
  SMTS_MISSING_FEMALES <- setdiff(SMTS_MALES$SMTS, SMTS_FEMALES$SMTS)
  
  # Displaying missing SMTS for each gender
  cat("Missing SMTS in males (MALES):", paste(SMTS_MISSING_MALES, collapse = ", "), "\n")
  cat("Missing SMTS in females (FEMALES):", paste(SMTS_MISSING_FEMALES, collapse = ", "), "\n")
  
  
  # Bar Chart (AGESTATE)
  # Define the order of categories in AGESTATE
  
  # Reorder the categories in AGESTATE variable
  Metadata$AGESTATE <- factor(Metadata$AGESTATE,
                              levels = c("Young", "Adult", "Old"))
  
  # Create a summary count of AGESTATE categories
  AGESTATE_SUMMARY <- Metadata %>%
    count(AGESTATE)
  
  # Calculate the total number of counts
  total_counts <- sum(AGESTATE_SUMMARY$n)
  
  # Create the bar chart with the total count as subtitle
  ggplot(AGESTATE_SUMMARY, aes(x = AGESTATE, y = n, fill = AGESTATE)) +
    geom_bar(stat = "identity", width = 0.7, color = "black") +
    geom_text(aes(label = n), vjust = -0.5, size = 4) +
    labs(title = "Count of AGESTATE Categories", x = "AGESTATE", y = "COUNTS",
         subtitle = paste("Total counts:", total_counts)) +  # Add total count as subtitle
    scale_fill_manual(values = c("Young" = "#F0F0F0", "Adult" = "#808080", "Old" = "#111111")) +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 14),  # Subtitle size
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  
  # Filter by sexes 1 (MALES) and 2 (FEMALES)
  # Filter the dataset Metadata for males and females based on SEX values 1 and 2 respectively.
  Males <- filter(Metadata, SEX == 1)
  Females <- filter(Metadata, SEX == 2)
  
  # Count the number of observations for males and females
  MalesCounts <- nrow(Males)
  FemalesCounts <- nrow(Females)
  
  # Display the counts of observations for males, females, and the total
  cat("Number of observations in MALES:", MalesCounts, "\n")
  cat("Number of observations in FEMALES:", FemalesCounts, "\n")
  cat("Number of TOTAL observations:", MalesCounts + FemalesCounts, "\n")
  
  
  # Bar chart (AGESTATE distribution for MALES)
  # Define the order of categories in AGESTATE for males
  Males$AGESTATE <- factor(Males$AGESTATE, levels = c("Young", "Adult", "Old"))
  
  # Create a summary count of AGESTATE categories for males
  MalesSummary <- Males %>%
    count(AGESTATE)
  
  # Calculate the total number of counts
  total_counts <- sum(MalesSummary$n)
  
  # Create a bar plot with visual adjustments for MALES
  ggplot(MalesSummary, aes(x = AGESTATE, y = n, fill = AGESTATE)) +
    geom_bar(stat = "identity", width = 0.7, color = "black") +
    geom_text(aes(label = n), vjust = -0.5, size = 4) +
    labs(title = "Count of categories of AGESTATE (MALES)", x = "AGESTATE", y = "COUNTS",
         subtitle = paste("Total counts:", MalesCounts)) +
    scale_fill_manual(values = c("Young" = "#F0F0F0", "Adult" = "#808080", "Old" = "#111111")) +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  
  # Bar chart (AGESTATE distribution for FEMALES)
  # Define the order of categories in AGESTATE for females
  Females$AGESTATE <- factor(Females$AGESTATE, levels = c("Young", "Adult", "Old"))
  
  # Create a summary count of AGESTATE categories for females
  FemalesSummary <- Females %>%
    count(AGESTATE)
  
  # Create a bar plot with visual adjustments for FEMALES
  ggplot(FemalesSummary, aes(x = AGESTATE, y = n, fill = AGESTATE)) +
    geom_bar(stat = "identity", width = 0.7, color = "black") +
    geom_text(aes(label = n), vjust = -0.5, size = 4) +
    labs(title = "Count of categories of AGESTATE (FEMALES)", x = "AGESTATE", y = "COUNTS",
         subtitle = paste("Total counts:", FemalesCounts)) +
    scale_fill_manual(values = c("Young" = "#F0F0F0", "Adult" = "#808080", "Old" = "#111111")) +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  
  # Create a new factor to order AGESTATE
  AGE_ORDER <- c("Young", "Adult", "Old")
  
  ## Summarize observation counts in the METADATA dataframe grouped by SMTS and AGESTATE variables.
  ## Create a new dataframe named COUNTS, including columns SMTS, AGESTATE, and COUNTS, where COUNTS represents the number of observations for each combination of SMTS and AGESTATE values.
  COUNTS_SMTS <- Metadata %>%
    group_by(SMTS, AGESTATE) %>%
    summarise(COUNTS = n()) %>%
    ungroup()
  
  ### Add the "AGESTATE" factor with the desired order
  COUNTS_SMTS <- COUNTS_SMTS %>%
    mutate(AGESTATE = factor(AGESTATE, 
                             levels = AGE_ORDER))
  
  ### Sort by SMTS, SEX, and then by the AGESTATE factor
  COUNTS_SMTS <- COUNTS_SMTS %>%
    arrange(SMTS, AGESTATE)
  
  # Pivot the COUNTS_SMTS table by the AGESTATE column
  COUNTS_SMTS <- COUNTS_SMTS %>%
    pivot_wider(names_from = AGESTATE, values_from = COUNTS)
  
  # Replace NA values with 0
  COUNTS_SMTS[is.na(COUNTS_SMTS)] <- 0
  
  # # Save the table to a CSV file
  # write.csv(COUNTS_SMTS, file = "data/NumberOfSamples/COUNTS_SMTS.csv", row.names = FALSE)
  
  ## Summarize observation counts in the METADATA dataframe grouped by SMTSD and AGESTATE variables.
  ## Create a new dataframe named COUNTS, including columns SMTSD, AGESTATE, and COUNTS, where COUNTS represents the number of observations for each combination of SMTSD and AGESTATE values.
  COUNTS_SMTSD <- Metadata %>%
    group_by(SMTSD, AGESTATE) %>%
    summarise(COUNTS = n()) %>%
    ungroup()
  
  ### Add the "AGESTATE" factor with the desired order
  COUNTS_SMTSD <- COUNTS_SMTSD %>%
    mutate(AGESTATE = factor(AGESTATE, 
                             levels = AGE_ORDER))
  
  ### Sort by SMTSD, SEX, and then by the AGESTATE factor
  COUNTS_SMTSD <- COUNTS_SMTSD %>%
    arrange(SMTSD, AGESTATE)
  
  # Pivot the COUNTS_SMTSD table by the AGESTATE column
  COUNTS_SMTSD <- COUNTS_SMTSD %>%
    pivot_wider(names_from = AGESTATE, values_from = COUNTS)
  
  # Replace NA values with 0
  COUNTS_SMTSD[is.na(COUNTS_SMTSD)] <- 0
  
  # # Save the table to a CSV file
  # write.csv(COUNTS_SMTSD, file = "data/NumberOfSamples/COUNTS_SMTSD.csv", row.names = FALSE)
  
  
  # Filter the METADATA dataframe to include only observations where the value in the SEX column is "1" (corresponding to males). 
  # Then, order the data by AGESTATE and count the unique combinations of SMTS and AGESTATE, creating a new column named "COUNTS".
  COUNTS_SMTS_MALES <- Metadata %>%
    filter(SEX == "1") %>%
    arrange(AGESTATE) %>%
    count(SMTS, 
          AGESTATE,
          name = "COUNTS")
  
  # Add the "AGESTATE" factor in the desired order
  COUNTS_SMTS_MALES <- COUNTS_SMTS_MALES %>%
    mutate(AGESTATE = factor(AGESTATE, 
                             levels = AGE_ORDER))
  
  # Sort by SMTS, SEX, and then by the factor AGESTATE
  COUNTS_SMTS_MALES <- COUNTS_SMTS_MALES %>%
    arrange(SMTS, 
            AGESTATE)
  
  # Pivot the COUNTS_SMTS_MALES table by the AGESTATE column
  COUNTS_SMTS_MALES <- COUNTS_SMTS_MALES %>%
    pivot_wider(names_from = AGESTATE, values_from = COUNTS)
  
  # Replace NA values with 0
  COUNTS_SMTS_MALES[is.na(COUNTS_SMTS_MALES)] <- 0
  
  # # Save the table to a CSV file
  # write.csv(COUNTS_SMTS_MALES, file = "data/NumberOfSamples/COUNTS_SMTS_MALES.csv", row.names = FALSE)
  
  
  # Filter the METADATA dataframe to include only observations where the value in the SEX column is "1" (corresponding to males). 
  # Then, order the data by AGESTATE and count the unique combinations of SMTSD and AGESTATE, creating a new column named "COUNTS".
  COUNTS_SMTSD_MALES <- Metadata %>%
    filter(SEX == "1") %>%
    arrange(AGESTATE) %>%
    count(SMTSD, 
          AGESTATE,
          name = "COUNTS")
  
  # Add the "AGESTATE" factor in the desired order
  COUNTS_SMTSD_MALES <- COUNTS_SMTSD_MALES %>%
    mutate(AGESTATE = factor(AGESTATE, 
                             levels = AGE_ORDER))
  
  # Sort by SMTSD, SEX, and then by the factor AGESTATE
  COUNTS_SMTSD_MALES <- COUNTS_SMTSD_MALES %>%
    arrange(SMTSD, 
            AGESTATE)
  
  # Pivot the COUNTS_SMTSD_MALES table by the AGESTATE column
  COUNTS_SMTSD_MALES <- COUNTS_SMTSD_MALES %>%
    pivot_wider(names_from = AGESTATE, values_from = COUNTS)
  
  # Replace NA values with 0
  COUNTS_SMTSD_MALES[is.na(COUNTS_SMTSD_MALES)] <- 0
  
  # # Save the table to a CSV file
  # write.csv(COUNTS_SMTSD_MALES, file = "data/NumberOfSamples/COUNTS_SMTSD_MALES.csv", row.names = FALSE)
  
  
  # Filter the METADATA dataframe to include only observations where the SEX column value is "2" (corresponding to females).
  # Then, arrange the data by AGESTATE and count unique combinations of SMTS and AGESTATE, creating a new column named "COUNTS".
  COUNTS_SMTS_FEMALES <- Metadata %>%
    filter(SEX == "2") %>%
    arrange(AGESTATE) %>%
    count(SMTS, 
          AGESTATE,
          name = "COUNTS")
  
  # Add the "AGESTATE" factor in the desired order
  COUNTS_SMTS_FEMALES <- COUNTS_SMTS_FEMALES %>%
    mutate(AGESTATE = factor(AGESTATE, 
                             levels = AGE_ORDER))
  
  # Sort by SMTS, SEX, and then by the factor AGESTATE
  COUNTS_SMTS_FEMALES <- COUNTS_SMTS_FEMALES %>%
    arrange(SMTS, 
            AGESTATE)
  
  # Pivot the COUNTS_SMTS table by the AGESTATE column
  COUNTS_SMTS_FEMALES <- COUNTS_SMTS_FEMALES %>%
    pivot_wider(names_from = AGESTATE, values_from = COUNTS)
  
  # Replace NA values with 0
  COUNTS_SMTS_FEMALES[is.na(COUNTS_SMTS_FEMALES)] <- 0
  
  # # Save the table to a CSV file
  # write.csv(COUNTS_SMTS_FEMALES, file = "data/NumberOfSamples/COUNTS_SMTS_FEMALES.csv", row.names = FALSE)
  
  # Filter the METADATA dataframe to include only observations where the SEX column value is "2" (corresponding to females).
  # Then, arrange the data by AGESTATE and count unique combinations of SMTSD and AGESTATE, creating a new column named "COUNTS".
  COUNTS_SMTSD_FEMALES <- Metadata %>%
    filter(SEX == "2") %>%
    arrange(AGESTATE) %>%
    count(SMTSD, 
          AGESTATE,
          name = "COUNTS")
  
  # Add the "AGESTATE" factor in the desired order
  COUNTS_SMTSD_FEMALES <- COUNTS_SMTSD_FEMALES %>%
    mutate(AGESTATE = factor(AGESTATE, 
                             levels = AGE_ORDER))
  
  # Sort by SMTSD, SEX, and then by the factor AGESTATE
  COUNTS_SMTSD_FEMALES <- COUNTS_SMTSD_FEMALES %>%
    arrange(SMTSD, 
            AGESTATE)
  
  # Pivot the COUNTS_SMTS table by the AGESTATE column
  COUNTS_SMTSD_FEMALES <- COUNTS_SMTSD_FEMALES %>%
    pivot_wider(names_from = AGESTATE, values_from = COUNTS)
  
  # Replace NA values with 0
  COUNTS_SMTSD_FEMALES[is.na(COUNTS_SMTSD_FEMALES)] <- 0
  
  # # Save the table to a CSV file
  # write.csv(COUNTS_SMTSD_FEMALES, file = "data/NumberOfSamples/COUNTS_SMTSD_FEMALES.csv", row.names = FALSE)
  
  
  # Organize data, count cases grouped by "SMTS" and "SEX", create "COUNTS" column with corresponding case numbers.
  COUNTS_SMTS_GENDER <- Metadata %>%
    arrange(SMTS) %>%
    count(SMTS, SEX, name = "COUNTS")
  
  # Create "Gender" column categorizing data into "MALES" and "FEMALES" based on "SEX" values, then remove "SEX" column.
  COUNTS_SMTS_GENDER <- COUNTS_SMTS_GENDER %>%
    mutate(Gender = ifelse(SEX == "1", "MALES", "FEMALES")) %>%
    dplyr::select(-SEX)
  
  # Pivot the COUNTS_SMTS_GENDER table by "SMTS" and "GENDER" columns
  COUNTS_SMTS_GENDER <- COUNTS_SMTS_GENDER %>%
    pivot_wider(names_from = Gender, values_from = COUNTS)
  
  # Replace NA values with 0
  COUNTS_SMTS_GENDER[is.na(COUNTS_SMTS_GENDER)] <- 0
  
  # # Save the table to a CSV file
  # write.csv(COUNTS_SMTS_GENDER, file = "data/NumberOfSamples/COUNTS_SMTS_GENDER.csv", row.names = FALSE)
  
  # Organize data, count cases grouped by "SMTSD" and "SEX", create "COUNTS" column with corresponding case numbers.
  COUNTS_SMTSD_GENDER <- Metadata %>%
    arrange(SMTSD) %>%
    count(SMTSD, SEX, name = "COUNTS")
  
  # Create "Gender" column categorizing data into "MALES" and "FEMALES" based on "SEX" values, then remove "SEX" column.
  COUNTS_SMTSD_GENDER <- COUNTS_SMTSD_GENDER %>%
    mutate(Gender = ifelse(SEX == "1", "MALES", "FEMALES")) %>%
    dplyr::select(-SEX)
  
  # Pivot the COUNTS_SMTSD_GENDER table by "SMTSD" and "GENDER" columns
  COUNTS_SMTSD_GENDER <- COUNTS_SMTSD_GENDER %>%
    pivot_wider(names_from = Gender, values_from = COUNTS)
  
  # Replace NA values with 0
  COUNTS_SMTSD_GENDER[is.na(COUNTS_SMTSD_GENDER)] <- 0
  
  # # Save the table to a CSV file
  # write.csv(COUNTS_SMTSD_GENDER, file = "data/NumberOfSamples/COUNTS_SMTSD_GENDER.csv", row.names = FALSE)
  
  
  # Reading Gene Reads Data
  # Read data from the specified file path using 'fread' function.
  GeneReads <- fread(
    "/data/RawData/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
    header = TRUE,
    sep = "\t"
  )
  
  
  # Filter 'METADATA' to keep rows where 'SAMPID' values are present in 'GENE_READS' columns onwards.
  Metadata <- Metadata[Metadata$SAMPID %in% colnames(GeneReads)[-c(1, 2)], ]
  
  # Divide 'Metadata' into subsets based on specific tissues.
  SMTS_FILT_LIST <- split(Metadata, Metadata$SMTS)
  
  # Create a list to store filtered results for each SMTS in GENE_READS.
  GENE_READS_SMTS_LIST <- setNames(lapply(names(SMTS_FILT_LIST), function(smts_name) {
    # Get the filtered data frame for the current SMTS.
    filtered_data <- SMTS_FILT_LIST[[smts_name]]
    
    # Select operation in GENE_READS for the current SMTS.
    GENE_READS_SMTS <- GeneReads %>%
      dplyr::select(Name, matches(filtered_data$SAMPID))
    
    # Return the result.
    GENE_READS_SMTS
  }), names(SMTS_FILT_LIST))
  
  # Save each element of GENE_READS_SMTS_LIST into separate files per tissue.
  lapply(names(GENE_READS_SMTS_LIST), function(smts_name) {
    # Get the data frame for the current SMTS.
    filtered_data <- GENE_READS_SMTS_LIST[[smts_name]]
    
    # Remove spaces in the smts_name.
    smts_name_no_spaces <- gsub(" ", "", smts_name)
    
    # Generate the filename without spaces.
    file_name <- paste0("data/ProcessedData/GeneReads/GeneReads", smts_name_no_spaces, ".csv")
    
    # # Save the data frame to a CSV file.
    # write.csv(filtered_data, file_name, row.names = FALSE)
  })
}


#### GSEA preparatory analysis ####
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Run this part after running 01_DufferebtuakExoression.R code ##
if(args[1] == "CreateRankedListOfGenes") {
  ## Read the differential expression results and generate a ranked list of genes based on FC (the file ends in .rnk) ##
  ## Just two columns, no row.names nor col.names. The first one contains ENSG IDs and the second one the FC ##
  
  ## ENSG Rank ##
  # Path to the directory containing the .txt files
  directory_path <- "Results/DifferentialExpression/"
    
  # Path to the directory where modified files will be saved
  destination_path <- "data/Rank/ENSG/"
    
  # Get the list of .txt files in the directory
  files <- list.files(directory_path, pattern = "\\.txt$", full.names = TRUE)
    
  # Process each .txt file in the directory
  for (file in files) {
    # Read the file using fread from data.table
    data <- fread(file, sep = "\t", header = FALSE)
      
    # Sort the data by column V2 (logFC) in descending order
    data <- data[order(-data$V2)]
      
    # Select only the first and second columns (ENSG and logFC) of the sorted data
    data <- data[, .(V1, V2)]
      
    # Create a new name for the modified file
    base_name <- tools::file_path_sans_ext(basename(file)) # Name without extension
    new_name <- paste0(destination_path, base_name, ".rnk")
      
    # Write the results to a new .rnk file in the specified folder
    write.table(data, file = new_name, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  
  ## GeneSymbol Rank ##
  # Load the biomaRt library
  library(biomaRt)
    
  # Path to the directory containing the .txt files
  directory_path <- "Results/DifferentialExpression/"
    
  # Path to the directory where modified files will be saved
  destination_path <- "data/Rank/GeneSymbol/"
    
  # Get the list of .txt files in the directory
  files <- list.files(directory_path, pattern = "\\.txt$", full.names = TRUE)
    
  # Function to obtain the GeneSymbol from ENSG using biomaRt
  get_gene_symbols <- function(ENSG_ids) {
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    gene_symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                          filters = "ensembl_gene_id", values = ENSG_ids, mart = ensembl)
    return(gene_symbols)
  }
    
  # Process each .txt file in the directory
  for (file in files) {
    # Read the file using fread from data.table
    data <- fread(file, sep = "\t", header = FALSE)
      
    # Sort the data by column V2 (logFC) in descending order
    data <- data[order(-data$V2)]
      
    # Remove the dot after the ENSG identifiers in column V1
    data$V1 <- gsub("\\..*", "", data$V1)
      
    # Select only the first column (ENSG) of the sorted data
    ENSG_ids <- data$V1
      
    # Get the corresponding GeneSymbols for ENSGs without dots
    gene_symbols <- get_gene_symbols(unique(ENSG_ids))
      
    # Convert ensembl_gene_id to character data type
    gene_symbols$ensembl_gene_id <- as.character(gene_symbols$ensembl_gene_id)
      
    # Merge the original data with GeneSymbols based on ENSG
    data <- merge(data, gene_symbols, by.x = "V1", by.y = "ensembl_gene_id", all.x = TRUE)
      
    # Select only columns V1, V2, and external_gene_name (GeneSymbol)
    data <- data[, c("hgnc_symbol", "V2")]
      
    # Sort the data by column V2 (logFC) in descending order
    data <- data[order(-data$V2)]
      
    # Create a new name for the modified file
    base_name <- tools::file_path_sans_ext(basename(file)) # Name without extension
    new_name <- paste0(destination_path, base_name, ".rnk")
      
    # Write the results to a new .rnk file in the specified folder
    write.table(data, file = new_name, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  # Directory where modified .rnk files were saved
  directory_path <- "data/Rank/GeneSymbol"
  
  # Directory where filtered .rnk files will be saved
  filtered_directory <- "data/Rank/GeneSymbolFiltered/"
  
  # Create directory if it doesn't exist
  dir.create(filtered_directory, showWarnings = FALSE)
  
  # Get the list of .rnk files in the directory
  files <- list.files(directory_path, pattern = "\\.rnk$", full.names = TRUE)
  
  # Function to remove rows without symbols and save to a filtered directory
  remove_rows_without_symbol <- function(file_path, filtered_directory) {
    # Read the .rnk file
    data <- fread(file_path, sep = "\t", header = FALSE)
    
    # Remove rows that do not have a Gene in the first column (V1)
    data <- data[!data$V1 == "", ]
    
    # Get the file name without path or extension
    base_name <- tools::file_path_sans_ext(basename(file_path))
    
    # Create the new file name for the filtered directory
    filtered_file_path <- file.path(filtered_directory, paste0(base_name, "_filtered.rnk"))
    
    # Save the updated file to the filtered directory
    write.table(data, file = filtered_file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  # Iterate over each .rnk file and save to the filtered directory
  for (file in files) {
    remove_rows_without_symbol(file, filtered_directory)
  }
}
  

## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Run this part after running CreateRankedListOfGenes code ##
if (args[1] == "GSEA") {
  
  # Get the list of .rnk files in the directory
  rnk_files <- list.files("data/Rank/GeneSymbolFiltered", pattern = "\\.rnk$", full.names = TRUE)
  
  # Base text for GSEA command
  base_text <- "Resources/GSEA/GSEA_4.3.2/gsea-cli.sh GSEAPreranked -gmx Resources/GSEA/GMTs/c2.cp.reactome.v2023.1.Hs.symbols.gmt,Resources/GSEA/GMTs/MetabolicTasks1KOSymbols.gmt -collapse No_Collapse -mode Abs_max_of_probes -norm meandiv -nperm 1000 -rnd_seed timestamp -scoring_scheme weighted -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -set_max 3000 -set_min 1 -zip_report false"
  
  # Path for output files
  output_path <- "Results/GSEA"
  
  # Get tissue names without spaces and with underscores
  tissues <- gsub(" ", "_", gsub("/", "_", gsub(".rnk", "", basename(rnk_files))))
  
  # Create and save commands for each tissue
  all_commands <- sapply(1:length(rnk_files), function(i) {
    tissue_text <- gsub(tissues[1], tissues[i], base_text)
    rnk_file <- rnk_files[i]
    chip_file <- "Resources/GSEA/CHIPs/Human_Gene_Symbol_with_Remapping_MSigDB.v2023.1.Hs.chip"
    output_dir <- file.path(output_path, tissues[i])
    full_command <- paste(tissue_text, "-rnk", rnk_file, "-chip", chip_file, "-out", output_dir)
    return(full_command)
  })
  
  # Save the commands in a text file
  writeLines(all_commands, file.path("/", "myrun_02GSEA.txt"))
  
  # # Path of the text file with tasks
  # input_file <- "myrun_02GSEA.txt"
  # 
  # # Output file names
  # adjusted_output <- "myrun_02GSEA_Adjusted.txt"
  # females_output <- "myrun_02GSEA_Female.txt"
  # males_output <- "myrun_02GSEA_Male.txt"
  # 
  # # Lists to store each type of task
  # adjusted_tasks <- c()
  # females_tasks <- c()
  # males_tasks <- c()
  # 
  # # Read the text file and separate tasks into corresponding lists
  # con <- file(input_file, "r")
  # while (TRUE) {
  #   line <- readLines(con, n = 1)
  #   if (length(line) == 0) {
  #     break
  #   }
  #   if (grepl("Adjusted", line)) {
  #     adjusted_tasks <- c(adjusted_tasks, line)
  #   } else if (grepl("Female", line)) {
  #     females_tasks <- c(females_tasks, line)
  #   } else if (grepl("Male", line)) {
  #     males_tasks <- c(males_tasks, line)
  #   }
  # }
  # close(con)
  # 
  # # Write the separated tasks into output files
  # writeLines(adjusted_tasks, adjusted_output)
  # writeLines(females_tasks, females_output)
  # writeLines(males_tasks, males_output)
  # 
  # cat("Tasks separated into corresponding files.\n")
}


## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
## Reactome pathways parent ##
if (args[1] == "Reactome") {

  # Read Reactome pathways sources from JSON file
  sources <- fromJSON(file = "Resources/Heatmap/c2.cp.reactome.v2023.2.Hs.json") %>%
    sapply(function(x) x$exactSource)
  
  # # Load pathways data for multiple sources
  # # Code is commented out due to its extensive running time
  # keep.columns <- c("displayName", "stId", "speciesName", "releaseDate")
  # df <- future_sapply(sources, function(id) {
  #   getPathways(id)
  # })
  # df.2 <- lapply(df, function(x) {
  #   data.frame(x)[, keep.columns]
  # }) %>%
  #   do.call(rbind, .)
  # df.parent <- future_sapply(sources, function(id) {
  #   getPathways(id, top.level = TRUE)$displayName
  # })
  # df.2$parent <- as.factor(df.parent)
  # df.2$speciesName <- as.factor(df.2$speciesName)
  # df.2$releaseDate <- as.factor(df.2$releaseDate)
  # df.2$geneSet <- names(df)
  # rownames(df.2) <- NULL
  # df.2 <- df.2 %>%
  #   select(c("geneSet", "displayName", "speciesName", "stId", "parent", "releaseDate")) %>%
  #   rename("Reactome Pathway" = "geneSet",
  #          "Description" = "displayName",
  #          "Code" = "stId",
  #          "Parent" = "parent",
  #          "Species" = "speciesName",
  #          "Release Date" = "releaseDate") %>%
  #   dplyr::arrange(Parent)
  # 
  # write.table(df.2, file = "~/MyTFM/Aging/Resources/Heatmap/Reactome.tsv",
  #             sep = "\t",
  #             quote = FALSE,
  #             row.names = FALSE)
  
  # Read Reactome pathways data from file
  Reactome <- fread("Resources/Heatmap/Reactome.tsv",
                    sep = "\t",
                    header = TRUE)
  
  # Path to the main folder
  MainPath <- "Results/GSEA/"
  
  # Get list of subdirectories (tissues)
  Subdirectories <- list.dirs(MainPath, full.names = FALSE, recursive = FALSE)
  
  # Create a data.table to store pathways
  MyReactomePathways <- data.table(Pathways = character(0))
  
  # Iterate through subdirectories and process TSV files
  for (tissue in Subdirectories) {
    tsv_path <- file.path(MainPath, tissue, paste0(tissue, "_GSEA_report.tsv"))
    
    # Check if TSV file exists
    if (file.exists(tsv_path)) {
      # Load data from TSV file into a data.table
      tsv_data <- fread(tsv_path, header = TRUE, sep = "\t")
      # Filter and store relevant pathway names (REACTOME_)
      reactome_pathways <- tsv_data$NAME[grepl("^REACTOME_", tsv_data$NAME)]
      # Add relevant names to the pathways data.table
      MyReactomePathways <- rbindlist(list(MyReactomePathways, data.table(Pathways = reactome_pathways)), fill = TRUE)
    }
  }
  
  # Filter unique pathways based on the "Pathways" column
  MyReactomePathways <- unique(MyReactomePathways)
  
  # Merge Reactome pathways data with extracted pathways
  MyReactomePathways <- merge(Reactome, MyReactomePathways, by.x = "Reactome Pathway", by.y = "Pathways", all.y = TRUE) %>%
    dplyr::arrange(Parent)
  
  # Write curated Reactome pathways data to a file
  # write.table(MyReactomePathways, file = "~/MyTFM/Aging/Resources/Heatmap/PreMyReactomePathways.tsv",
  #             sep = "\t",
  #             quote = FALSE,
  #             row.names = FALSE)
  
  # Manually curated data
  MyReactomePathwaysCurated <- fread("Resources/Heatmap/MyReactomePathwaysCurated.tsv",
                                     sep = "\t",
                                     header = TRUE) %>%
    dplyr::arrange(Parent)
  
  # Write curated Reactome pathways to a file
  # write.table(MyReactomePathwaysCurated, file = "~/MyTFM/Aging/Resources/Heatmap/MyReactomePathwaysCurated.tsv",
  #             sep = "\t",
  #             quote = FALSE,
  #             row.names = FALSE)
}















