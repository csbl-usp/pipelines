library(dplyr)
library(fgsea)

params <- snakemake@params
params <- params[names(params) != ""]

# Check input files
if (is.null(snakemake@input[["ranks"]])) {
    stop("Must provide ranks file")
}
rankfile <- snakemake@input[["ranks"]]

if (is.null(snakemake@input[["gmt"]])) {
    stop("Must provide gmt file")
}
gmtfile <- snakemake@input[["gmt"]]

ranks <- readr::read_tsv(rankfile)

sort_by <- params[["sort_by"]]
params[["sort_by"]] <- NULL
if (is.null(sort_by)) {
    numericols <- select_if(ranks, is.numeric)
    if (ncol(numericols) == 0) {
        stop("At least one numeric column must exist for ranking")
    } else if (ncol(numericols) > 1) {
        stop("If more than 1 numeric columns exist,
             'sort_by' parameter must be given")
    } else {
        sort_by <- colnames(numericols)
    }
}

ranks <- arrange(ranks, desc(.data[[sort_by]]))
rankvector <- pull(ranks, sort_by)
names(rankvector) <- pull(ranks, "ID")

input <- list(pathways=gmtPathways(gmtfile),
              stats=rankvector)

all <- append(input, params)

# Run fgsea
result <- do.call(fgsea, all)

# Save results
outfile <- snakemake@output[[1]]
if (is.null(outfile)) {
    outfile <- "gsea.tsv"
}
write.table(result, outfile, sep="\t", quote=F)
