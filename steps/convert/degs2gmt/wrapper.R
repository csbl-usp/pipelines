#' Extract up and down regulated genes from DEGs table.
#' 
#' @param degs DataFrame with differentially expressed genes.
#' @param gene_col Gene column name.
#' @param fc_col Fold change column name.
#' @param pv_col P-value column name.
#' @param comp_col Comparisons column name.
#' @param pv_cutoff P-value cutoff.
#' @param fc_cutoff Fold change cutoff.
#' @return List of up and down regulated genes for each comparison.
degs2list <- function(degs, gene_col='gene', fc_col='logfc', 
                      pv_col='adjpval', comp_col='comparison',
                      pv_cutoff=0.1, fc_cutoff=0) {
    degs <- degs[which(degs[[gene_col]] != ""), ]
    up <- degs[which(degs[[fc_col]] > fc_cutoff & degs[[pv_col]] < pv_cutoff), gene_col]
    down <- degs[which(degs[[fc_col]] < -fc_cutoff & degs[[pv_col]] < pv_cutoff), gene_col]
    name <- unique(degs[[comp_col]])
    results <- list(up, down)
    names(results) <- c(paste0(name, "_up"), paste0(name, "_down"))
    return(results)
}

#' Turns degs list to .GMT format.
#' 
#' @param degs_list List of gene sets.
#' @return DataFrame in GMT format.
list2gmt <- function(degs_list) {
    max_ln <- max(sapply(degs_list, length))
    results <- lapply(degs_list, function(deg) {
        dif <- max_ln - length(deg)
        result <- c(deg, rep("", dif))
    })
    results <- do.call('rbind', results)
    results <- cbind(names(degs_list), results)
    return(results)
}

# Get file path and parameters
degs_file <- snakemake@input[["degs"]]
params <- snakemake@params
params <- params[names(params) != ""]

# Read DEGs table and transform into GMT
degs <- read.csv(file=degs_file, header=TRUE, stringsAsFactors=FALSE)

gps <- unique(degs[, ifelse(is.null(params[["group"]]), "group", params[["group"]])])
all_groups <- lapply(gps, function(gp){
    df <- degs[which(degs$group == gp),]
    params[["degs"]] <- df
    up_down <- do.call(degs2list, params)
    return(up_down)
})

all_groups <- unlist(all_groups, recursive=F)
all_groups <- all_groups[lapply(all_groups, function(x){length(x)}) > 0]
gmt <- list2gmt(all_groups)

if (length(snakemake@output) > 0) {
    output <- snakemake@output[[1]]
} else {
    output <- 'degs.gmt'
}

write.table(file=output, gmt, sep='\t', col.names=FALSE, row.names=TRUE, quote=FALSE)
