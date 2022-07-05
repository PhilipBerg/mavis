get_conditions <- function(design) {
  design %>%
    colnames() %>%
    stringr::str_flatten("|")
}
