input <- list()

# Check input files
if (!is.null(snakemake@input[["expression"]])) {
    expression <- readr::read_tsv(snakemake@input[["expression"]])
    expression <- as.data.frame(expression)
    rownames(expression) <- expression[,"ID"]
    expression[,1] <- NULL
    input[["expr"]] <- expression
}

if (!is.null(snakemake@input[["sample_annotation"]])) {
    annot <- readr::read_tsv(snakemake@input[["sample_annotation"]])
    input[["annot"]] <- annot
}

if (!is.null(snakemake@input[["gmt"]])) {
    gmt <- CEMiTool::read_gmt(snakemake@input[["gmt"]])
    input[["gmt"]] <- gmt
}

if (!is.null(snakemake@input[["interactions"]])) {
    interact <- readr::read_tsv(snakemake@input[["interactions"]])
    input[["interactions"]] <- interact
}

# Clean parameters list
params <- snakemake@params
params <- params[names(params) != ""]
all <- append(input, params)

# Run cemitool
cem <- do.call(CEMiTool::cemitool, all)
