
#' Computes abundance recovery rate
#'
#' @param se summarized experiment
#' @param sample samples to calculate ARR for
#' @param reference references (corresponding to samples)
#' @param log TRUE if data should be logged
#' @param pseudocount pseudocount for log
#' @param statistic statistic to aggregate by: mean or median
#' @param features features to calculate metric for
#'
#' @return a list object with plots and metrics
#'
#' @export
computeRecovery <- function(se, sample, reference, log = TRUE, pseudocount = 1,
                            statistic = 'median', features = NULL) {
    result <- computeMetric(se, sample, reference, metric = 'ARR', log = log,
                            pseudocount = pseudocount, statistic = statistic, features = features)
    return(result)
}

#' Computes relative difference
#'
#' @param se summarized experiment
#' @param sample samples to calculate RD for
#' @param reference references (corresponding to samples)
#' @param log TRUE if data should be logged
#' @param pseudocount pseudocount for log
#' @param statistic statistic to aggregate by: mean or median
#' @param features features to calculate metric for
#'
#' @return a list object with plots and metrics
#'
#' @export
computeDifference <- function(se, sample, reference, log = TRUE, pseudocount = 1,
                            statistic = 'median', features = NULL) {
    result <- computeMetric(se, sample, reference, metric = 'RD', log = log,
                            pseudocount = pseudocount, statistic = statistic, features = features)
    return(result)
}

#' Computes normalized root mean squared error
#'
#' @param se summarized experiment
#' @param sample samples to calculate NRMSE for
#' @param reference references (corresponding to samples)
#' @param log TRUE if data should be logged
#' @param pseudocount pseudocount for log
#' @param statistic statistic to aggregate by: mean or median
#' @param features features to calculate metric for
#'
#' @return metrics
#'
#' @export
computeNRMSE <- function(se, sample, reference, log = TRUE, pseudocount = 1,
                              statistic = 'median', features = NULL) {
    result <- computeMetric(se, sample, reference, metric = 'NRMSE', log = log,
                            pseudocount = pseudocount, statistic = statistic, features = features)
    return(result$metrics)
}

#' Computes reproducibility
#'
#' @param se summarized experiment
#' @param samples samples to calculate reproducibility
#' @param log TRUE if data should be logged
#' @param pseudocount pseudocount for log
#' @param features features to calculate metric for
#' @param pt.alpha alpha for points
#' @param pt.size size for points
#' @param contour TRUE if contours should be displayed
#' @param density TRUE if density should be displayed
#' @param contour.alpha alpha for contour
#' @param contour.size size for contour
#' @param density.alpha alpha for density
#' @param spline.degree degree of spline fit
#'
#' @return a list object with plots and metrics
#'
#' @importFrom matrixStats rowMeans2
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab theme_minimal theme geom_smooth
#'     geom_density_2d scale_color_viridis_c geom_density_2d_filled
#' @importFrom splines bs
#' @importFrom SummarizedExperiment assays
#' @importFrom stats complete.cases lm
#'
#' @export
computeReproducibility <- function(se, samples, features = NULL, log = TRUE,
                                   pseudocount = 1, pt.alpha = 0.5, pt.size = 0.9,
                                   contour = TRUE, density = FALSE,
                                   contour.alpha = 0.5, contour.size = 0.7,
                                   density.alpha = 0.5, spline.degree = 6) {
    ..samples <- ..level.. <- abundance <- sds <- NULL
    if (!is.null(features)) se <- se[features,]

    values <- assays(se)$val
    values <- values[,..samples]
    if (log) values <- log2(values + pseudocount)
    values <- values[complete.cases(values),]
    u <- rowMeans(values)
    s <- sqrt(rowMeans((values - u)^2))
    metric <- sqrt(mean(s^2))
    data <- data.frame(abundance = u, sds = s)
    scatter <- ggplot(data, aes(x = abundance, y = sds)) +
        geom_point(alpha = pt.alpha, size = pt.size) +
        xlab('average abundance') + ylab('standard deviation') +
        theme_minimal() + theme(legend.position = 'none') +
        geom_smooth(method = lm, formula = y ~ bs(x, spline.degree), se = FALSE)


    if (contour) scatter <- scatter + geom_density_2d(alpha = contour.alpha,
        size = contour.size, aes(color = ..level..)) + scale_color_viridis_c()
    if (density) scatter <- scatter + geom_density_2d_filled(alpha = density.alpha)

    return(list(plot = scatter, metric = metric))
}

#' Computes consistency
#'
#' @param se summarized experiment
#' @param samples samples to calculate reproducibility
#' @param log TRUE if data should be logged
#' @param pseudocount pseudocount for log
#' @param features features to calculate metric for
#' @param pt.alpha alpha for points
#' @param pt.size size for points
#' @param thresholds thresholds to test consistency for
#' @param statistic statistic to aggregate metrics by
#'
#' @return a list object with plots and metrics
#'
#' @importFrom matrixStats colMeans2 colSds
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab theme_minimal theme geom_smooth
#'    geom_errorbar
#' @importFrom SummarizedExperiment colData assays
#' @importFrom utils combn
#' @importFrom stats median complete.cases
#'
#' @export
computeConsistency <- function(se, samples, thresholds = seq(0, 2,length.out = 21),
                               statistic = 'median', features = NULL, log = TRUE,
                               pseudocount = 1, pt.alpha = 0.5, pt.size = 2) {
    ..samples <- ..keep <- threshold <- consistency <- method <- TRUE
    if (!is.null(features)) se <- se[features,]
    cdata <- colData(se)[samples,]

    values <- assays(se)$val
    values <- values[,..samples]
    if (log) values <- log2(values + pseudocount)
    values <- values[complete.cases(values),]

    methods <- unique(cdata$method)
    num.methods <- length(methods)
    metrics_methods <- lapply(methods, function(method) {

        keep <- names(values) %in% cdata$ID[cdata$method == method]
        val <- values[, ..keep]

        comb <- combn(colnames(val), 2)
        metrics <- sapply(thresholds, function(threshold) {
            consis <- sapply(seq_len(ncol(comb)), function(index) {
                val1 <- val[[comb[1,index]]]
                val2 <- val[[comb[2,index]]]
                val1 <- ifelse(val1 < threshold, TRUE, FALSE)
                val2 <- ifelse(val2 < threshold, TRUE, FALSE)
                agreement <- table(val1==val2)
                cons <- agreement['TRUE'] / nrow(values)
                return(cons)
            })
            return(consis)
        })
    })

    means <- lapply(metrics_methods, colMeans2)
    sds <- lapply(metrics_methods, colSds)

    data <- data.frame(threshold = rep(thresholds, times = num.methods),
                       consistency = unlist(means),
                       method = rep(methods, each = length(thresholds)))
    scatter <- ggplot(data, aes(x = threshold, y = consistency, color = method)) +
        geom_point(alpha = pt.alpha, size = pt.size) +
        geom_errorbar(aes(ymin = consistency - unlist(sds),
                          ymax = consistency + unlist(sds)),
                      width = 0.2) +
        theme_minimal()

    func <- get(statistic, mode = 'function')
    metric_stat <- lapply(metrics_methods, function(x) {
        x <- apply(x, 2, func)
        names(x) <- paste0('threshold=',thresholds)
        return(x)
    })
    names(metric_stat) <- methods

    return(list(plot = scatter, metric = metric_stat))
}

#' Computes resolution entropy
#'
#' @param se summarized experiment
#' @param samples samples to calculate reproducibility
#' @param log TRUE if data should be logged
#' @param pseudocount pseudocount for log
#' @param features features to calculate metric for
#' @param pt.alpha alpha for points
#' @param pt.size size for points
#' @param num.bins bins to test resolution entropy for
#' @param statistic statistic to aggregate metrics by
#'
#' @return a list object with plots and metrics
#'
#' @importFrom matrixStats colMeans2 colSds
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab theme_minimal theme geom_errorbar position_dodge
#' @importFrom SummarizedExperiment colData assays
#' @importFrom utils combn
#' @importFrom stats median complete.cases
#'
#' @export
computeResEntropy <- function(se, samples, num.bins = seq(5, 20, length.out = 16),
                              statistic = 'median', features = NULL, log = TRUE,
                              pseudocount = 1, pt.alpha = 0.5, pt.size = 2) {
    ..samples <- ..keep <- threshold <- entropy <- method <- NULL

    if (!is.null(features)) se <- se[features,]
    cdata <- colData(se)[samples,]

    values <- assays(se)$val
    values <- values[,..samples]
    if (log) values <- log2(values + pseudocount)
    values <- values[complete.cases(values),]

    methods <- unique(cdata$method)
    num.methods <- length(methods)
    metrics_methods <- lapply(methods, function(method) {
        keep <- names(values) %in% cdata$ID[cdata$method == method]
        val <- values[, ..keep]
        metrics <- sapply(num.bins, function(num.bin) {
            entropy <- sapply(seq_len(ncol(val)), function(index){
                val_sample <- val[[index]]
                bin <- cut(val_sample, num.bin)
                nbin <- as.numeric(table(bin)) / length(val_sample) + 1e-10
                re <- -sum(nbin * log(nbin))
                return(re)
            })
        })
        return(metrics)
    })

    means <- lapply(metrics_methods, colMeans)
    sds <- lapply(metrics_methods, matrixStats::colSds)

    data <- data.frame(threshold = rep(num.bins, times = num.methods),
                       entropy = unlist(means),
                       method = rep(methods, each = length(num.bins)))
    scatter <- ggplot(data, aes(x = threshold, y = entropy, color = method)) +
        geom_point(alpha = pt.alpha, size = pt.size) +
        geom_errorbar(aes(ymin = entropy - unlist(sds), ymax = entropy + unlist(sds)), width = 0.2,
                      position = position_dodge(0.05)) +
        theme_minimal()

    func <- get(statistic, mode = 'function')
    metric_stat <- lapply(metrics_methods, function(x) {
        x <- apply(x, 2, median)
        names(x) <- paste0('nbin=',num.bins)
        return(x)
    })
    names(metric_stat) <- methods

    return(list(plot = scatter, metric = metric_stat))
}

##### Helpers #####

#' @importFrom SummarizedExperiment assays
getvalues <- function(se, sample, reference, log, pseudocount) {

    values <- assays(se)$val
    x <- values[[sample]]
    y <- values[[reference]]

    nna <- !is.na(y) & y > 1e-5
    y <- y[nna]
    x <- x[nna]
    x[is.na(x)] <- 0

    if (log) {
        x <- log2(x + pseudocount)
        y <- log2(y + pseudocount)
    }

    return(list(sample = x, reference = y))
}

#' @importFrom stats sd
getMetric <- function(metric, sample, reference) {

    if (metric == 'ARR') return(sample / reference)
    if (metric == 'RD') return( abs(sample - reference) / reference)
    if (metric == 'NRMSE') return( sqrt(sum((sample - reference)^2) / length(reference))/sd(reference))

}

#' @importFrom ggplot2 aes theme_minimal theme element_blank ylab geom_violin geom_boxplot
#' @importFrom stats median
computeMetric <- function(se, sample, reference, metric, log = TRUE, pseudocount = 1,
                          statistic = 'median', features = NULL) {

    if (!is.null(features)) se <- se[features]
    func <- get(statistic, mode = 'function')

    metrics <- Map(function(sam, ref) {

        values <- getvalues(se, sam, ref, log, pseudocount)
        sample.val <- values$sample
        reference.val <- values$reference

        metric_sample <- getMetric(metric, sample.val, reference.val)
        data <- data.frame(metric = metric_sample, sample = sam)

        metric_stat <- func(metric_sample)
        return(list(data = data, statistic = metric_stat))

    }, sample, reference)

    data <- do.call('rbind', lapply(metrics, `[[`, 1))

    if (metric == 'NRMSE') {
        violin <- NULL
    } else {
        violin <- ggplot(data, aes(x = sample, y = metric, fill = sample)) +
            theme_minimal() + theme(axis.text.x = element_blank()) + ylab(metric) +
            geom_violin(alpha = 0.5, color = 'transparent') +
            geom_boxplot(width = 0.1, outlier.alpha = 0.5, outlier.size = 0.5)

    }

    metric_stat <- unlist(lapply(metrics, `[[`, 2))

    return(list(plot = violin, metrics = metric_stat))
}




