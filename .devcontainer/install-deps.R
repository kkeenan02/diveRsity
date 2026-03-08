# Install all diveRsity dependencies inside the devcontainer
# Run automatically on container creation via postCreateCommand

cat("Installing diveRsity dependencies...\n")

pkgs <- c(
  # Hard dependencies
  "Rcpp", "ggplot2", "shiny", "qgraph", "methods",
  # Suggested
  "xlsx", "plotrix", "HWxtest",
  # Dev tools
  "devtools", "rcmdcheck"
)

install.packages(pkgs, dependencies = TRUE)

cat("\nDone. To build and check the package:\n")
cat("  Rscript -e \"Rcpp::compileAttributes()\"\n")
cat("  R CMD build /workspace\n")
cat("  R CMD check --as-cran diveRsity_*.tar.gz\n")
cat("\nTo run with Valgrind:\n")
cat("  R -d valgrind --vanilla -e \"devtools::load_all(); data(Test_data)\"\n")
