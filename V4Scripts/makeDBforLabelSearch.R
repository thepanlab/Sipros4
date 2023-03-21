library(Biostrings)
library(stringr)

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if (length(args) < 1) {
  args <- c("-help")
}

## Help section
if ("-help" %in% args) {
  cat(
    "
Extract identified protein sequence in .pro.txt from the fasta file of protein database

Arguments:
-pro   .pro.txt file path
-faa   .faa file path
-o     output .faa file path
-help  get help
Example:
Rscript ./makeDBforLabelSearch.R -pro demo.pro.txt -faa demo.faa -o demo.label.faa\n"
  )
  q(save = "no")
}

## Parse arguments
pro <- ""
faa <- ""
out <- ""
for (i in 1:length(args)) {
  if (args[i] == "-pro")
    pro <- args[i + 1]
  if (args[i] == "-faa")
    faa <- args[i + 1]
  if (args[i] == "-o")
    out <- args[i + 1]
}

## get new faa db

pro <- read.table(pro,
                  sep = "\t",
                  quote = "",
                  header = T)
pro <- pro$ProteinID
# split protein name in {}
pro2 <- pro[str_sub(pro, 1, 1) == "{"]
pro <- pro[str_sub(pro, 1, 1) != "{"]
pro2 <- str_sub(pro2, 2, -2)
pro2 <- str_split(pro2, ",")
pro2 <- unlist(pro2)
pro <- unique(c(pro, pro2))
# remove docoy
pro <- pro[str_sub(pro, 1, 4) != "Rev_"]
pro <- pro[str_sub(pro, 1, 5) != "Rev1_"]
pro <- pro[str_sub(pro, 1, 5) != "Rev2_"]
fa <- readAAStringSet(faa)
faPro <- names(fa)
# get the first chunk of in line begining with ">"
faPro <- str_split(faPro, " ", 2, T)[, 1]
fa <- fa[match(pro, faPro)]
writeXStringSet(fa, out)
