# Use R-devel on Ubuntu (mirrors CRAN's Linux check environment)
FROM rocker/r-devel:latest

# System dependencies for diveRsity and Rcpp compilation
RUN apt-get update && apt-get install -y \
    # Build tools
    build-essential \
    gfortran \
    # Valgrind for memory checking (catches C++ bugs like the glbWCcpp segfault)
    valgrind \
    # Git and utilities
    git \
    curl \
    # qgraph dependency
    libglu1-mesa-dev \
    libgl1-mesa-dev \
    # xlsx / Java dependency
    default-jdk \
    r-cran-rjava \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Reconfigure Java for R (needed for xlsx package)
RUN R CMD javareconf

# Set CRAN mirror
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' >> /usr/local/lib/R/etc/Rprofile.site
