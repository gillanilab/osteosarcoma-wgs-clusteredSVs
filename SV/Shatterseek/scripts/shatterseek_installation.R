# install_shatterseek_complete.R
# Complete installation script for ShatterSeek with Bioconductor dependencies
print("Starting package installation...")

options(ask = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

# If GitHub API rate limits are hit when installing R packages, set a GitHub PAT:
# Sys.setenv(GITHUB_PAT = "token_here")

# Install BiocManager if not already installed
if (!require(BiocManager, quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager")
}

# Install devtools if not already installed
if (!require(devtools, quietly = TRUE)) {
  cat("Installing devtools...\n")
  install.packages("devtools")
}

library(BiocManager)
library(devtools)

# Install required Bioconductor packages
cat("Installing Bioconductor dependencies...\n")
bioc_packages <- c(
  "IRanges",
  "GenomicRanges", 
  "S4Vectors",
  "graph",
  "BiocGenerics"
)

for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing", pkg, "from Bioconductor...\n"))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(paste(pkg, "is already installed.\n"))
  }
}

# Install additional CRAN dependencies if needed
cran_packages <- c("gridExtra")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing", pkg, "from CRAN...\n"))
    install.packages(pkg)
  }
}

# install ShatterSeek from GitHub
cat("Installing ShatterSeek from GitHub...\n")
tryCatch({
  install_github("parklab/ShatterSeek", auth_token = NULL)
  cat("ShatterSeek installed successfully!\n")
}, error = function(e) {
  cat("Error installing ShatterSeek:", conditionMessage(e), "\n")
})

# test
cat("Testing ShatterSeek installation...\n")
tryCatch({
  library(ShatterSeek)
  cat("✓ ShatterSeek loaded successfully!\n")
  

  if (exists("packageVersion")) {
    tryCatch({
      version <- packageVersion("ShatterSeek")
      cat("ShatterSeek version:", as.character(version), "\n")
    }, error = function(e) {
      cat("Version info not available.\n")
    })
  }
  
}, error = function(e) {
  cat("✗ Error loading ShatterSeek:", conditionMessage(e), "\n")
})

cat("Installation script completed.\n")