#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom data.table :=
#' @importFrom data.table as.data.table
#' @importFrom data.table fread
#' @importFrom data.table setkey
#' @importFrom data.table setnames
#' @importFrom magrittr %>%
#' @importFrom parallel mclapply
#' @importFrom parallel mcmapply
#' @importFrom stats sd
#' @importFrom stringi stri_detect_regex
#' @importFrom stringi stri_locate_all_regex
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble
#' @importFrom tidytable all_of
#' @importFrom tidytable select
#' @importFrom tidytable unite
#' @importFrom tidytable unnest
#' @importFrom tools file_path_sans_ext
## usethis namespace: end
utils::globalVariables(c(".","batch", "dummy", "forw_primer_pos", "gen_length", "id", "qual", "rc_forw_primer_pos"))
