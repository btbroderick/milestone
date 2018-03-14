zipped <- function(...)
{
  unlist(lapply(list(...), paste0, c("", ".gz", ".bz2", ".xz", ".zip")))
}

read_any_file <- function(fp)
{
  ext <- tools::file_ext(fp)
  if(ext %in% zipped("csv"))
  {
    out <- read_csv(fp, col_names = TRUE, col_types = cols())
  } else if(ext %in% zipped("txt"))
  {
    out <- read_table2(fp, col_names = TRUE, col_types = cols())
  } else if(ext %in% zipped("tab", "tsv"))
  {
    out <- read_tsv(fp, col_names = TRUE, col_types = cols())
  } else if(ext %in% "sas7bdat")
  {
    out <- haven::read_sas(fp)
  } else if(ext %in% c("xlsx", "xls"))
  {
    out <- readxl::read_excel(fp)
  } else if(ext %in% "rds")
  {
    out <- readRDS(fp)
  }
  attr(out, "extension") <- ext
  out
}