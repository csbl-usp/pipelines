library(enrichR)

# Helper functions
add_name_column <- function(name, lst, colname) {
    if ( nrow(lst[[name]]) > 0 ) {
        df <- lst[[name]]
        df[, colname] <- name
        return(df)
    }
}
genelists_file <- snakemake@input[["genelists"]]
params <- snakemake@params


if (is.null(genelists_file)) {
    stop("Must provide genelists file")
}

if (!is.null(params[["databases"]])) {
    dbs <- params[["databases"]]   
} else {
    dbs <- listEnrichrDbs()[, "libraryName"]
}

genelists <- fgsea::gmtPathways(genelists_file)

result <- lapply(genelists, function(genelist) {
    enriched <- enrichr(genelist, dbs)
    nrows <- lapply(enriched, nrow)
    if (!any(nrows>=1)) { 
        stop('None of the genelists were enriched')
    }
    enriched <- lapply(names(enriched), 
                       add_name_column, enriched, "database")
    do.call(rbind, enriched)
})

result <- lapply(names(result), add_name_column, result, "genelist")
result <- do.call(rbind, result)

if (length(snakemake@output) > 0) {
    output <- snakemake@output[[1]]
} else {
    output <- 'enriched.tsv'
}

write.table(result, output, sep="\t", row.names=F, quote=F)
