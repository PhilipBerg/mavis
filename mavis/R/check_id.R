check_id <- function(data, id) {
  cl <- data %>%
    dplyr::select(!!rlang::enquo(id)) %>%
    purrr::map_chr(class)
  common_types <- c("character", "factor")
  if (any(cl == "list")) {
    msg <- paste0(
      "The following columnn(s): ",
      paste0(names(cl[cl == "list"]), collapse = ", "),
      " are lists"
    )
    rlang::abort(
      paste0(
        "Some id column(s) are lists, this is not allowed, did you specify the id_col correctly?\n",
        msg
      )
    )
  }
  else if (!all(cl %in% common_types)) {
    types <- cl[!cl %in% common_types] %>%
      unique()
    msg <- c()
    for (i in seq_along(types)) {
      msg[i] <- paste0(
        "The following columnn(s): ",
        paste0(names(cl[cl == types[i]]), collapse = ", "),
        " are ",
        types[i]
      )
    }
    rlang::warn(
      paste0(
        "Not all id columns are characters or factors, did you specify the id_col correctly?\n",
        paste0(msg, collapse = "\n")
      )
    )
  }
  data %>%
    dplyr::select(!!rlang::enquo(id)) %>%
    names()
}
