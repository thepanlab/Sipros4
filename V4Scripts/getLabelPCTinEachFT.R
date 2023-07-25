library(stringr)
library(tidyr)

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
Get identified protein count and spectra count of each .FT2 file

Arguments:
-pro   .pro.txt file path
-psm   .psm.txt file path
-thr   a number for labeled threshold, larger than this means labeled spectra
-o     output .LabelPCTcount.txt and .LabelPCTcountSummary.txt file path frefix
-help  get help
Example:
Rscript ./getLabelPCTinEachFT.R -pro demo.pro.txt -psm demo.psm.txt -thr 5 -o demo\n"
    )
    q(save = "no")
}

#### Parse arguments ####
pro <- ""
psm <- ""
thr <- 5
out <- ""
for (i in 1:length(args)) {
    if (args[i] == "-pro") {
        pro <- args[i + 1]
    }
    if (args[i] == "-psm") {
        psm <- args[i + 1]
    }
    if (args[i] == "-thr") {
        thr <- as.numeric(args[i + 1])
    }
    if (args[i] == "-o") {
        out <- args[i + 1]
    }
}

##### get label PCT in each FT2 file  ####

pro <- read.table(pro,
    sep = "\t",
    quote = "",
    header = T
)
pro <- pro[stringr::str_detect(pro$ProteinID, "Rev_",
    negate = T
), ]
pro <- pro[stringr::str_detect(pro$ProteinID, "Rev2_",
    negate = T
), ]

psm <-
    read.table(psm,
        sep = "\t",
        quote = "",
        header = T
    )
psm <- psm[stringr::str_detect(psm$ProteinNames, "Rev_",
    negate = T
), ]
psm <- psm[stringr::str_detect(psm$ProteinNames, "Rev2_",
    negate = T
), ]

fts <- unique(psm$Filename)
pros <- pro$ProteinID

resultDf <-
    as.data.frame(matrix(nrow = length(pros), ncol = length(fts) * 3))
colnames(resultDf) <-
    c(
        paste0(fts, "_PCTcount"),
        paste0(fts, "_medianPCTcount"),
        paste0(fts, "_labeledCount")
    )

pctDf <- data.frame(
    FileName = psm$Filename,
    PCT = as.numeric(str_sub(
        str_split(psm$SearchName, "_", simplify = T)[, 2], 1, -4
    )) / 1000,
    # for accurate match, remove "{}", add "," in the end of names
    proteinNames = paste0(str_sub(psm$ProteinNames, 2, -2), ",")
)

# find each protein group's label PCT in each ft file
FoundPCTinFt <- function(ftName) {
    pctFtDf <- pctDf[pctDf$FileName == ftName, ]
    pctss <- c()
    medianPct <- c()
    labeledCount <- c()
    for (proName in pros)
    {
        # split protein names in {}
        if (str_sub(proName, 1, 1) == "{") {
            proName <- str_sub(proName, 2, -2)
            proName <- str_split(proName, ",")
            proName <- unlist(proName)
        }
        pcts <- c()
        for (pro in proName)
        {
            # accurate match. For exmaple "7_78" and "7_7" is different
            pro <- paste0(pro, ",")
            # fixed, in case symbol in pro "|" for example
            psmFound <- str_detect(pctFtDf$proteinNames, fixed(pro))
            pcts <- c(pcts, pctFtDf$PCT[psmFound])
            # in case protein in group repeat count
            pctFtDf <- pctFtDf[!psmFound, ]
        }
        pctss <- c(pctss, paste(pcts, collapse = ","))
        medianPct <- c(medianPct, median(pcts))
        labeledCount <- c(labeledCount, sum(pcts > thr, na.rm = T))
    }
    return(data.frame(
        PCTs = pctss,
        MedianPCT = medianPct,
        LabeledCount = labeledCount
    ))
}

findLabeledPepPSMInFT <- function(ftName) {
    psmInFt <- psm[psm$Filename == ftName, ]
    PCTs <- psmInFt$SearchName
    PCTs <- str_remove_all(PCTs, "\\{")
    PCTs <- str_remove_all(PCTs, "\\}")
    PCTs <- str_remove_all(PCTs, "Pct")
    PCTs <- str_split(PCTs, "_", simplify = T)[, 2]
    PCTs <- as.numeric(PCTs) / 1000
    pepSeqs <- psmInFt$OriginalPeptide[PCTs > thr]
    return(c(
        LabeledPepCount = length(unique(pepSeqs)),
        LabeledPSMcount = length(pepSeqs)
    ))
}

for (i in (1:length(fts))) {
    df <- FoundPCTinFt(fts[i])
    resultDf[, i] <- df[, 1]
    resultDf[, length(fts) + i] <- df[, 2]
    resultDf[, length(fts) * 2 + i] <- df[, 3]
}

resultDf <- cbind(proteinNames = pros, resultDf)
resultDf[resultDf == ""] <- NA

write.table(
    resultDf,
    paste0(out, ".LabelPCTcount.txt"),
    sep = "\t",
    row.names = F,
    col.names = T,
    quote = F
)

labeledPepsAndPSMs <- t(sapply(fts, findLabeledPepPSMInFT))

summaryDf <- data.frame(
    FileNames = fts,
    LabeledProteins = colSums(resultDf[
        ,
        (length(fts) * 2 + 2):(length(fts) * 3 + 1)
    ] > 0, na.rm = T)
)
summaryDf <- cbind(summaryDf, labeledPepsAndPSMs)

write.table(
    summaryDf,
    paste0(out, ".LabelPCTcountSummary.txt"),
    sep = "\t",
    row.names = F,
    col.names = T,
    quote = F
)
