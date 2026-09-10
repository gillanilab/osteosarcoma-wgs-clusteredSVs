print("Starting package installation...")

options(ask = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

# If GitHub API rate limits are hit when installing R packages, set a GitHub PAT:
# Sys.setenv(GITHUB_PAT = "token_here")

# Install devtools if missing
if (!requireNamespace("devtools", quietly = TRUE)) {
  print("Installing devtools...")
  install.packages("devtools")
}

print("Installing gert...")
install.packages("gert")
print("Gert installed")

print("Installing usethis...")
install.packages("usethis")
print("Usethis installed")

print("Installing gUtils...")
devtools::install_github("mskilab/gUtils")
print("gUtils installed")

# Install fishHook
print("Installing fishHook...")
devtools::install_github("mskilab/fishhook")
print("fishHook installed")

print("All installations complete!")
