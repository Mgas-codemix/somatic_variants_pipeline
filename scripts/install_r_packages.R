#!/usr/bin/env Rscript
# Install R packages required by the pipeline (if not using Conda/Docker).
# Run: Rscript scripts/install_r_packages.R

pkgs <- c("data.table", "xlsx")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing ", p, " ...")
    install.packages(p, repos = "https://cloud.r-project.org/")
  } else {
    message(p, " already installed")
  }
}
message("Done.")
