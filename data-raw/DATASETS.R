## code to prepare the data sets presented in Mavis
# Get the data
files <- list.files('data-raw/', pattern = '\\.csv', recursive = T, full.names = T)
files <- setNames(
  files,
  paste0(
    rep(c('ramus', 'ups', 'yeast'), each = 2),
    '_',
    rep(c('max', 'prog'), times = 3)
  )
)
# Assign to global environment
purrr::map(files, readr::read_csv) %>%
  purrr::walk2(
    .x = names(.),
    .y = .,
    .f = ~ assign(x = .x,
                  value = .y,
                  envir = .GlobalEnv)
  )

# Attach to package
usethis::use_data(ramus_max, ups_max, yeast_max, yeast_prog, overwrite = TRUE)
usethis::use_data(ramus_prog, ups_prog, overwrite = TRUE, compress = 'xz')

# Clean up global environment
rm(list = names(files))
rm(files)
