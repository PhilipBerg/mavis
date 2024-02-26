# Compile vignettes
file_extention <- "\\.orig$"
vigs <- list.files("./vignettes", pattern = file_extention, full.names = TRUE)
for (i in vigs) {
  e <- new.env()
  rmd_file <- gsub(file_extention, "", i)
  knitr::knit(i, rmd_file, envir = e, )
  tx  <- readLines(rmd_file)
  tx2  <- gsub(pattern = "\\]\\(\\./vignettes/", replace = "\\]\\(", x = tx)
  writeLines(tx2, rmd_file)
}

# Clean environment
rm(list = ls())
