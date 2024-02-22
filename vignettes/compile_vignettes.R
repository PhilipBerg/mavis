# Compile vignettes
file_extention <- "\\.orig$"
vigs <- list.files("./vignettes", pattern = file_extention, full.names = TRUE)
lapply(vigs
  ~ knitr::knit(.x, gsub(file_extention, .x), envir = new.env())
)
# knitr::knit("vignettes/baldur-tutorial.Rmd.orig", "vignettes/baldur-tutorial.Rmd", envir = new.env())
# knitr::knit("vignettes/mult-imp.Rmd.orig", "vignettes/mult-imp.Rmd", envir = new.env())

# Move images
lapply(list.files("./figure/", full.names = TRUE), file.copy, "./vignettes/")
unlink("./figure", recursive = TRUE)

# Change image path in compiled vignettes
for(i in list.files("./vignettes/", pattern = "Rmd$", full.names = TRUE)) {
  tx  <- readLines(i)
  tx2  <- gsub(pattern = "\\]\\(figure/", replace = "\\]\\(", x = tx)
  writeLines(tx2, i)
}

rm(tx, tx2, i)
