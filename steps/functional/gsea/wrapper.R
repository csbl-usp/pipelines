library(dplyr)
library(fgsea)

params <- snakemake@params
params <- params[names(params) != ""]
# Helper functions
add_name_column <- function(name, lst, colname) {
    if ( nrow(lst[[name]]) > 0 ) {
        df <- lst[[name]]
        df[, colname] <- name
        return(df)
    }
}

# Check input files
if (is.null(snakemake@input[["ranks"]])) {
    stop("Must provide ranks file")
}
rankfile <- snakemake@input[["ranks"]]

if (is.null(snakemake@input[["gmt"]])) {
    stop("Must provide gmt file")
}
gmtfile <- snakemake@input[["gmt"]]


if (is.null(params[["comp_col"]])) {
    comp_col <- "group"
    params[["comp_col"]] <- NULL
} else {
    comp_col <- params[["comp_col"]]
    params[["comp_col"]] <- NULL

}

if (is.null(params[["gene_col"]])) {
    gene_col <- "ID"
    params[["gene_col"]] <- NULL
} else {
    gene_col <- params[["gene_col"]]
    params[["gene_col"]] <- NULL
}

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

comp_groups <- unique(ranks[[comp_col]])
result <- lapply(comp_groups, function(comp_group) {
    ranks_fil <- ranks[which(ranks[[comp_col]] == comp_group), ] 
   
    ranks_fil <- arrange(ranks_fil, desc(.data[[sort_by]]))
    rankvector <- pull(ranks_fil, sort_by)
    names(rankvector) <- pull(ranks_fil, gene_col)
    
    input <- list(pathways=gmtPathways(gmtfile),
                  stats=rankvector)
    all_input <- append(input, params)
    # Run fgsea
    print(all_input)
    do.call(fgsea, all_input)
})
print(result)
names(result) <- comp_groups
result <- lapply(names(result), add_name_column, result, "comparison")
result <- do.call(rbind, result)
result[["leadingEdge"]] <- unlist(lapply(result[["leadingEdge"]], paste, collapse = ","))

# Save results
outfile <- snakemake@output[[1]]
if (is.null(outfile)) {
    outfile <- "gsea.tsv"
}
write.table(result, outfile, sep="\t", quote=F, row.names=FALSE)
