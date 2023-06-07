get_conditions <- function(design) {
  design %>%
    colnames() %>%
    paste0('^', ., collapse = '|')
}
