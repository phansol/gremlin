if (!suppressWarnings(require("rpart"))){install.packages("rpart", repos = "http://cran.r-project.org")}
if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
if (!suppressWarnings(require("GWASTools"))){BiocManager::install("GWASTools")}
