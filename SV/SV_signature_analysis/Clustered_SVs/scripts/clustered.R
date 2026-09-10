library(signature.tools.lib)

# Read BEDPE file w/ SVs for all samples
bedpe_data <- read.table("../../../../ctDNA/data/leopard_cohort/clustered_sv_input.bedpe", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Get unique samples
samples <- unique(bedpe_data$sample)

# Initialize lists to store data from all samples
all_annotated_bedpe <- list()
all_clustering_regions <- list()
rearr_catalogues <- list()

# Process each sample and extract the catalogue
for(sample in samples) {
  sample_data <- bedpe_data[bedpe_data$sample == sample, ]
  result <- bedpeToRearrCatalogue(sample_data)

  # Extract the catalogue data frame from the result
  if(is.list(result) && "rearr_catalogue" %in% names(result)) {
    # Convert data frame to named vector
    catalogue_vector <- result$rearr_catalogue[[1]]  # Get the first (and only) column
    names(catalogue_vector) <- rownames(result$rearr_catalogue)  # Use row names
    rearr_catalogues[[sample]] <- catalogue_vector
  }

  # Extract the annotated bedpe with clustering information
  if(is.list(result) && "annotated_bedpe" %in% names(result)) {
    annotated_bedpe <- result$annotated_bedpe
    
    # The annotated_bedpe should contain an 'is.clustered' column
    if("is.clustered" %in% colnames(annotated_bedpe)) {
      # Add sample column if it doesn't exist
      if(!"sample" %in% colnames(annotated_bedpe)) {
        annotated_bedpe$sample <- sample
      }
      
      # Store in list
      all_annotated_bedpe[[sample]] <- annotated_bedpe
    }
  }
  
  # Extract clustering regions if available
  if(is.list(result) && "clustering_regions" %in% names(result)) {
    clustering_regions <- result$clustering_regions
    
    # Add sample column if it doesn't exist
    if(!"sample" %in% colnames(clustering_regions)) {
      clustering_regions$sample <- sample
    }
    
    # Store in list
    all_clustering_regions[[sample]] <- clustering_regions
  }
}

# Combine all annotated_bedpe data from all samples into one data frame
if(length(all_annotated_bedpe) > 0) {
  combined_annotated_bedpe <- do.call(rbind, all_annotated_bedpe)
  
  # Write to file
  write.table(combined_annotated_bedpe, 
              file = "../../../../ctDNA/leopard_cohort/outputs/all_samples_annotated_v2_bedpe.txt", 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
  
  cat("Written", nrow(combined_annotated_bedpe), "annotated rearrangements from", 
      length(all_annotated_bedpe), "samples to all_samples_annotated_bedpe.txt\n")
}

# Combine all clustering_regions data from all samples into one data frame
if(length(all_clustering_regions) > 0) {
  combined_clustering_regions <- do.call(rbind, all_clustering_regions)
  
  # Write to file
  write.table(combined_clustering_regions, 
              file = "../../../../ctDNA/leopard_cohort/outputs/all_samples_clustering_regions_v2.txt", 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
  
  cat("Written", nrow(combined_clustering_regions), "clustering regions from", 
      length(all_clustering_regions), "samples to all_samples_clustering_regions.txt\n")
} else {
  cat("No clustering regions found across all samples\n")
}

# Combine catalogues
rearr_catalogue <- do.call(cbind, rearr_catalogues)