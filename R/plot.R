
#' Plot histogram
#'
#' @param se summarized experiment
#' @param samples samples to calculate reproducibility
#' @param log TRUE if data should be logged
#' @param pseudocount pseudocount for log
#' @param histogram TRUE if histogram should be displayed
#' @param density TRUE if density should be displayed
#' @param alpha alpha for density
#'
#' @return histogram
#'
#' @importFrom data.table melt
#' @importFrom ggplot2 ggplot aes theme_minimal geom_density geom_histogram
#' @importFrom SummarizedExperiment assays
#'
#' @export
plotHistogram <- function(se, samples, log = TRUE, pseudocount = 0.1,
                          histogram = TRUE, density = FALSE, alpha = 0.5) {

    ..samples <- ..density.. <- NULL
    values <- assays(se)$val
    values <- values[, ..samples]
    values <- values[rowSums(values) > 1e-10,]
    values <- melt(values, id.vars = character())
    names(values) <- c('sample', 'expression')
    values <- values[!is.na(values$expression)]

    if (log) values$expression <- log2(values$expression + pseudocount)

    histo <- ggplot(data = data.frame(values), aes(x = expression, fill = sample)) + theme_minimal()
    if (density) histo <- histo + geom_density(alpha = alpha)
    if (histogram) histo <- histo + geom_histogram(aes(y=..density..), alpha = alpha, color = 'grey', position = 'identity')

    return(histo)
}

#' Plot scatter
#'
#' @param se summarized experiment
#' @param sample.x sample 1
#' @param sample.y sample 2
#' @param features features to plot
#' @param annotate metadata to annotate by
#' @param cor.method correlation method
#' @param log TRUE if data should be logged
#' @param pseudocount pseudocount for log
#' @param pt.alpha alpha for points
#' @param pt.size size for points
#' @param contour TRUE if contours should be displayed
#' @param density TRUE if density should be displayed
#' @param contour.alpha alpha for contour
#' @param contour.size size for contour
#' @param density.alpha alpha for density
#' @param diagonal TRUE if diagonal should be displayed
#'
#' @return scatter plots
#'
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab xlim ylim labs theme_minimal theme
#'     geom_density_2d scale_color_viridis_c geom_density_2d_filled geom_abline aes_string guides
#' @importFrom SummarizedExperiment rowData assays
#' @importFrom stats cor
#'
#' @export
plotScatter <- function(se, sample.x, sample.y, features = NULL,
                        annotate = NULL, cor.method = 'spearman', log = TRUE,
                        pseudocount = 1, contour = TRUE, density = FALSE,
                        diagonal = TRUE, pt.alpha = 0.5, pt.size = 3,
                        contour.alpha = 0.5, contour.size = 0.7, density.alpha = 0.5) {
    ..level.. <- ..density.. <- NULL

    if (!is.null(features)) se <- se[features]
    rdata <- rowData(se)

    values <- assays(se)$val
    x <- values[[sample.x]]
    y <- values[[sample.y]]
    nna <- !is.na(x) & !is.na(y) & x > 0 & y > 0
    x <- x[nna]
    y <- y[nna]
    rdata <- rdata[nna,]
    if (log) {
        x <- log2(x + pseudocount)
        y <- log2(y + pseudocount)
    }

    ## Correlation
    corIJ <- cor(x, y, method = cor.method, use = 'complete.obs')

    ## Limits
    max.lim <- max(max(x),max(y))
    min.lim <- min(min(x),min(y))

    ## Plot
    data <- data.frame(x = x, y = y, rdata[,-1])

    if (is.null(annotate)) {
        scatter <- ggplot(data = data, aes(x,y)) + geom_point(alpha = pt.alpha, size = pt.size) +
            theme_minimal() + theme(legend.position = 'none') +
            xlab(sample.x) + ylab(sample.y) + xlim(min.lim,max.lim) + ylim(min.lim,max.lim) +
            labs(title = paste0('SCC = ', round(corIJ,3)))
    } else {
        scatter <- ggplot(data = data, aes_string(x="x", y ="y", size = annotate)) +
            geom_point(alpha = pt.alpha) +
            theme_minimal() + theme(legend.position = 'right') +
            xlab(sample.x) + ylab(sample.y) + xlim(min.lim,max.lim) + ylim(min.lim,max.lim) +
            labs(title = paste0('SCC = ', round(corIJ,3)))
    }

    if (contour) scatter <- scatter + geom_density_2d(
        alpha = contour.alpha, size = contour.size, aes(color = ..level..)) +
        scale_color_viridis_c() + guides(color='none')
    if (density) scatter <- scatter + geom_density_2d_filled(alpha = density.alpha) +
        guides(fill='none')
    if (diagonal) scatter <- scatter + geom_abline(intercept = 0, slope = 1,
                                                   color = 2, alpha = 0.5, linetype = 'dashed')

    return(scatter)
}

