# Subset table
subSetTable <- function(tb, gene_col, logfc_col,p_value_col, p_value, logfc, nm_col) {
    head(tb)
    tb <- tb[which(tb[[gene_col]] != ""), ]
    up <- tb[which(tb[[logfc_col]] > logfc & tb[[p_value_col]] < p_value), gene_col]
    down <- tb[which(tb[[logfc_col]] < -logfc & tb[[p_value_col]] < p_value), gene_col]
    name <- unique(tb[[nm_col]])
    nm <- c(paste0(name, "_up"), paste0(name, "_down"))
    results <- list(up, down)
    names(results) <- nm
    return(results)
}

# List to dataframe
list2Gmt <- function(tb) {
    print(class(tb))
    max_ln <- max(sapply(tb, length))
    results <- lapply(tb, function(deg) {
        dif <- max_ln - length(deg)
        result <- c(deg, rep("", dif))
    })
    results <- do.call('rbind', results)
    results <- cbind(names(tb), results)
    return(results)
}

# Get file path and parameters
degs_file <- snakemake@input[["degs"]]
params <- snakemake@params
params <- NULL
if (is.null(params)) {
    gene_col <- "HGNC.symbol"
    logfc_col <- "log2FoldChange"
    p_value_col <- "padj"
    p_value <- 0.1
    logfc <- 0
    nm_col <- "group"
} else {
    gene_col <- params[["gene_col"]]
    logfc_col <- params[["logfc_col"]]
    p_value_col <- params[["p_value_col"]]
    p_value <- as.numeric(params[["p_value"]])
    logfc <- as.numeric(params[["logfc"]])
    nm_col <- params[["group"]]
}

# Read DEGs table and transform into GMT
degs_file <- "~/projects/fabiopohl/zika/meta/results/DEGs/all_degs.csv"
degs <- read.csv(file=degs_file, header=TRUE, stringsAsFactors=FALSE)

gps <- unique(degs[, nm_col])
all_groups <- lapply(gps, function(gp){
    mydf <- degs[which(degs$group == gp),]
    up_down <- subSetTable(mydf, gene_col, logfc_col, p_value_col, p_value, logfc, nm_col)
    return(up_down)
})

all_groups <- unlist(all_groups, recursive=F)
all_groups <- all_groups[lapply(all_groups, function(x){length(x)})>0]
gmt <- list2Gmt(all_groups)

if (length(snakemake@output) > 0) {
    output <- snakemake@output[[1]]
} else {
    output <- 'degs.gmt'
}

write.table(file=output, gmt, sep='\t', col.names=FALSE, row.names=TRUE, quote=FALSE)
