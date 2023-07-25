library(stringr)

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
-o     output .SPcount.txt and .ProCountSummary.txt file path frefix
-help  get help
Example:
Rscript ./getSpectraCountInEachFT.R -pro demo.pro.txt -psm demo.psm.txt -o demo\n"
    )
    q(save = "no")
}

#### Parse arguments ####
pro <- ""
psm <- ""
out <- ""
for (i in 1:length(args)) {
    if (args[i] == "-pro") {
        pro <- args[i + 1]
    }
    if (args[i] == "-psm") {
        psm <- args[i + 1]
    }
    if (args[i] == "-o") {
        out <- args[i + 1]
    }
}

##### get spectrum count in each FT2 file  ####

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

# for accurate match, remove "{}", add "," in the end of names
psm$ProteinNames <- paste0(str_sub(psm$ProteinNames, 2, -2), ",")

fts <- unique(psm$Filename)
pros <- pro$ProteinID

resultMat <- matrix(nrow = length(pros), ncol = length(fts) * 2)
colnames(resultMat) <-
    c(paste0(fts, "_Found"), paste0(fts, "_PSMcount"))

# find each protein group has how many psm in each ft file
isFoundInFt <- function(ftName) {
    ftPsm <- psm[psm$Filename == ftName, ]
    ftProNames <- ftPsm$ProteinNames
    isFounds <- c()
    psmCounts <- c()
    for (proName in pros)
    {
        # split protein names in {}
        if (str_sub(proName, 1, 1) == "{") {
            proName <- str_sub(proName, 2, -2)
            proName <- str_split(proName, ",")
            proName <- unlist(proName)
        }
        isFound <- F
        psmCount <- 0
        for (pro in proName)
        {
            # accurate match. For exmaple "7_78" and "7_7" is different
            pro <- paste0(pro, ",")
            # fixed, in case symbol in pro "|" for example
            psmFound <- str_detect(ftProNames, fixed(pro))
            psmCount <- psmCount + sum(psmFound)
            # in case protein in group repeat count
            ftProNames <- ftProNames[!psmFound]
        }
        if (psmCount > 0) {
            isFound <- T
        }
        isFounds <- c(isFounds, isFound)
        psmCounts <- c(psmCounts, psmCount)
    }
    return(data.frame(isFound = isFounds, psmCount = psmCounts))
}

for (i in (1:length(fts))) {
    df <- isFoundInFt(fts[i])
    resultMat[, i] <- df[, 1]
    resultMat[, length(fts) + i] <- df[, 2]
}

resultDf <- cbind(proteinNames = pros, as.data.frame(resultMat))
proFoundDf <- colSums(resultDf[, -1])
proFoundDf <-
    data.frame(ftName = names(proFoundDf), proCount = proFoundDf)

write.table(
    resultDf,
    paste0(out, ".SPcount.txt"),
    sep = "\t",
    row.names = F,
    col.names = T,
    quote = F
)

write.table(
    proFoundDf,
    paste0(out, ".ProCountSummary.txt"),
    sep = "\t",
    row.names = F,
    col.names = T,
    quote = F
)
