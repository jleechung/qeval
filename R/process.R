
#' Construct summarized experiment
#'
#' @param path path to quantification results
#' @param methods methods of corresponding results
#'
#' @importFrom data.table fread `:=`
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @return a summarized experiment
#'
#' @export
constructSE <- function(path, methods) {

    files <- list.files(path, full.names = TRUE)
    if (length(files) != length(methods)) stop('Specify method for each file')

    ## Read files and construct summarizedExperiment
    message('Reading ', length(files) ,' files')
    files <- list.files(path, full.names = TRUE)
    dlst <- lapply(files, fread)
    dlst <- Map(function(x,y) {
        d <- fread(x)
        names(d)[-1] <- paste0(names(d)[-1],'_', y)
        return(d)
    }, files, methods)
    names(dlst) <- files

    ## Check for duplicated rows
    dlst <- deduplicate(dlst)

    ## Store sampleIDs
    sampleID <- unlist(lapply(dlst, function(x) names(x)[-1]))
    names(sampleID) <- NULL

    ## Store corresponding methods
    sampleN <- sapply(dlst, ncol)-1
    sampleMethods <- rep(methods, times = sampleN)

    ## Store files
    sampleFiles <- rep(files, times = sampleN)

    ## Store column sums
    sampleColsums <- unlist(lapply(dlst, function(x) colSums(x[,-1])))

    ## Store number of nonzeros
    sampleNonzeros <- unlist(lapply(dlst, function(x) colSums(x[,-1] > 0)))

    ## Store number of entries
    sampleEntries <- rep(unlist(lapply(dlst, function(x) nrow(x))), times = sampleN)

    ## Combine assays
    message('Merge assays')
    d <- Reduce(function(x,y) merge(x,y,by='ID',all=TRUE), dlst)
    tx.id <- d$ID
    d[,"ID":=NULL]

    ## Create se
    message('Constructing summarized experiment')
    assay <- list(values = d)
    coldata <- data.frame(ID = sampleID, method = sampleMethods,
                          colsum = sampleColsums, entries = sampleEntries,
                          nonzero = sampleNonzeros, files = sampleFiles)
    se <- SummarizedExperiment(assays = assay, colData = coldata)
    rownames(se) <- tx.id
    colnames(se) <- sampleID

   return(se)
}

#' Construct summarized experiment
#'
#' @param se a summarized experiment
#' @param gtf annotation
#'
#' @importFrom data.table data.table .N
#' @importFrom SummarizedExperiment SummarizedExperiment rowData<-
#'
#' @return a summarized experiment
#'
#' @export
annotateSE <- function(se, gtf) {

    . <- NULL
    data <- data.table(tx = gtf$transcript_id, type = gtf$type, width = gtf@ranges@width)
    type <- count <- width <- tx <- NULL
    data <- data[type=='exon']
    data <- data[, .(count = .N, length = sum(width)), by = tx]

    rdata <- data.frame(tx.id = rownames(se),
                        num.exons = rep(NA, nrow(se)),
                        length = rep(NA, nrow(se)))
    rdata$num.exons <- data$count[match(rdata$tx.id, data$tx)]
    rdata$length <- data$length[match(rdata$tx.id, data$tx)]
    rowData(se) <- rdata

    return(se)
}

#' @importFrom data.table .SD
deduplicate <- function(dlst) {
    ID <- NULL
    n <- length(dlst)
    for (i in seq_len(n)) {
        d <- dlst[[i]]
        dups <- any(duplicated(d$ID))
        if (!dups) next
        message('Found duplicate rows in file ', names(dlst)[i])
        ids <- d$ID[duplicated(d$ID)]
        message('Duplicated IDs: ', paste(ids, collapse = ', '))
        message('Taking mean of duplicated entries')
        d <- d[, lapply(.SD, mean), by = ID]
        dlst[[i]] <- d
    }
    return(dlst)
}
