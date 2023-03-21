library(stringr)
library(tidyr)
library(dplyr)

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if (length(args) < 1) {
  args <- c("-help")
}

#####  Help section #####
if ("-help" %in% args) {
  cat(
    "
Filter out proteins with low PSM score to decrease FDR

Arguments:
-pro   .pro.txt file path
-psm   .psm.txt file path
-fdr   fdr value you want to refine, 0.01 for example
-sip   add this flag if you want to convert score to 0~1 
-o     output .proRefineFDR.txt file path prefix
-help  get help
Example:
Rscript ./refineProteinFDR.R -pro demo.pro.txt -psm demo.psm.txt -fdr 0.01 -o demo\n"
  )
  q(save = "no")
}

#### Parse arguments ####
pro <- ""
psm <- ""
fdr <- 0.01
out <- ""
sip <- F
for (i in 1:length(args)) {
  if (args[i] == "-pro")
    pro <- args[i + 1]
  if (args[i] == "-psm")
    psm <- args[i + 1]
  if (args[i] == "-fdr")
    fdr <- as.numeric(args[i + 1])
  if (args[i] == "-sip")
    sip <- T
  if (args[i] == "-o")
    out <- args[i + 1]
}

##### refine protein FDR  ####
proDf <- read.table(pro,
                    sep = "\t",
                    quote = "",
                    header = T)
psm <-
  read.table(psm,
             sep = "\t",
             quote = "",
             header = T)
newScore <- psm$Score
# convert score to probability
if (sip)
{
  data <-
    psm[, c("ParentCharge",
            "MassErrorDa",
            "Score",
            "ProteinCount",
            "TargetMatch")]
  pct <- psm$SearchName
  pct <- str_split(pct, "_", simplify = T)[, 2]
  pct <- str_sub(pct, 1,-4)
  pct <- as.numeric(pct) / 1000
  data$Pct <- pct
  model <-
    glm(TargetMatch ~ .,
        family = "binomial",
        data = data)
  newScore <- predict(model, data, type = "response")
}
pros <- proDf$ProteinID
psmScoreDf <-
  data.frame(ProteinNames = psm$ProteinNames, Score = newScore)
psmScoreDf$ProteinNames <- str_sub(psmScoreDf$ProteinNames, 2, -2)
psmScoreDf <- separate_rows(psmScoreDf, ProteinNames, sep = ",")
findHighestScore <- function(proName)
{
  # split protein names in {}
  if (str_sub(proName, 1, 1) == "{")
  {
    proName <- str_sub(proName, 2,-2)
    proName <- str_split(proName, ",")
    proName <- unlist(proName)
  }
  psmFound <- psmScoreDf$ProteinNames %in% proName
  scores <- psmScoreDf$Score[psmFound]
  # remove found psm to reduce complex
  psmScoreDf <- psmScoreDf[!psmFound, ]
  ## max is better than product
  # return(1 - prod(1 - scores))
  return(max(scores))
}
highestScores <- sapply(pros, findHighestScore)
proScoreDf <- cbind(proDf, Scores = highestScores)
proScoreDf <- arrange(proScoreDf, Scores)
N <- nrow(proScoreDf)
TargetMatchs <- proScoreDf$TargetMatch
decoyN <- sum(!TargetMatchs)
Scores <- proScoreDf$Scores
FDR <- 0
tresh <- 0
for (i in 1:N) {
  FDR <- decoyN / (N - i + 1)
  tresh <- Scores[i]
  if (FDR <= fdr)
    break
  if (!TargetMatchs[i])
    decoyN <- decoyN - 1
}
# in case i=1 when no decoy exists
if(i>1) {
  proScoreDf <-
    proScoreDf[-1:-(i - 1), !(colnames(proScoreDf) %in% "Scores")]
}
comments <-
  comments <-
  paste0("# Decoy_Proteins_After_Filtering = ", decoyN, "\n")
comments <-
  paste0(comments,
         "# Target_Proteins_After_Filtering = ",
         N - i + 1 - decoyN,
         "\n")
comments <-
  paste0(comments, "# Total_Proteins_After_Filtering = ", N - i + 1, "\n")
comments <- paste0(comments, "# Protein_FDR = ", FDR, "\n")
writeLines(comments, paste0(out, ".proRefineFDR.txt"))
write.table(
  proScoreDf,
  paste0(out, ".proRefineFDR.txt"),
  sep = "\t",
  row.names = F,
  col.names = T,
  append = T,
  quote = F
)
