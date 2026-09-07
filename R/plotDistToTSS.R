merge_two_si = function(x1, x2){
  if (length(unique(gsub("^[0-9]+","",c(x1, x2)))) == 1){
    return(paste0(gsub("[^0-9]*$","",x1), "-", x2))
  } else {
    return(paste0(x1, "-", x2))
  }
}

generate_break_lbs = function(breaks) {
  lbs = c()

  # break labels
  break_labels = scales::label_number(scale_cut = scales::cut_si(unit = "b"))(breaks)
  break_labels = gsub(" b$"," bp", break_labels)

  # category labels
  for (i in 2:length(breaks)) {
    if (i == length(breaks)) {
      lbs = c(lbs, paste0(">", break_labels[i-1]))
    } else {
      lbs = c(lbs, merge_two_si(break_labels[i-1], break_labels[i]))
    }
  }

  return(lbs)
}

generate_colors = function(palette = NULL, n) {
  # old color in version <= 1.41.1
  old_color = c("#9ecae1", "#3182bd", "#C7A76C", "#86B875", "#39BEB1", "#CD99D8")
  if (is.null(palette)){
    brewer_cols = old_color
  } else if (length(palette) == 1 && is_valid_palette(palette)){
    brewer_cols = RColorBrewer::brewer.pal(
      name = palette,
      n = RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    ) |> rev()
  } else if (all(is_valid_color(palette))){
    brewer_cols = palette
  }
  else {
    warning("Your palette is non-valid, switching to default...")
    brewer_cols = old_color
  }

  if (length(brewer_cols) >= n) {
    cols = brewer_cols[1:length(brewer_cols)]
  } else {
    cols = grDevices::colorRampPalette(brewer_cols)(n)
  }

  return(cols)
}

is_valid_palette = function(palette){
  palette %in% rownames(RColorBrewer::brewer.pal.info)
}

is_valid_color = function(color){
  tryCatch({
    grDevices::col2rgb(color)
    TRUE
  }, error = function(e) {
    FALSE
  })
}

##' Plot feature distribution based on distances to TSS
##'
##' This function creates a bar plot showing the distribution of peaks relative
##' to transcription start sites (TSS), categorizing peaks by their distance
##' from the nearest gene's TSS.
##'
##' @description
##' The function visualizes where ChIP-seq peaks are located relative to gene
##' TSSs. It categorizes peaks into distance bins (e.g., 0-1kb, 1-3kb, etc.)
##' and displays the percentage of peaks in each category, separately for
##' upstream (5') and downstream (3') regions. The plot uses a diverging bar
##' chart with TSS at the center (0).
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Categorizes peaks into distance bins based on \code{distanceBreaks}
##'   \item Separates upstream (negative distance) and downstream (positive
##'         distance) peaks
##'   \item Calculates percentages for each distance category
##'   \item Creates a diverging bar plot with:
##'     \itemize{
##'       \item Upstream peaks shown as negative bars (left side)
##'       \item Downstream peaks shown as positive bars (right side)
##'       \item TSS marked at the center (0)
##'       \item Distance categories color-coded
##'     }
##'   \item Supports grouping by a category column for comparing multiple datasets
##' }
##'
##' @param peakDist data.frame containing peak annotation data with a distance
##'   column. Typically obtained from \code{as.data.frame(csAnno)} or similar
##' @param distanceColumn character, name of the column containing distances
##'   from peaks to nearest gene TSS. Default is "distanceToTSS"
##' @param distanceBreaks numeric vector, breakpoints for distance categories
##'   in base pairs. The function automatically adds 0 and Inf if not present.
##'   Default is c(0, 1000, 3000, 5000, 10000, 100000)
##' @param palette character, color palette name from RColorBrewer (e.g., "Set1",
##'   "Dark2") or NULL for default colors. Run
##'   \code{RColorBrewer::display.brewer.all()} to see available palettes.
##'   Default is NULL
##' @param xlab character, label for x-axis. Default is "" (empty)
##' @param ylab character, label for y-axis. Default is "Binding sites (\%) (5'->3')"
##' @param title character, plot title. Default is "Distribution of transcription
##'   factor-binding loci relative to TSS"
##' @param categoryColumn character or numeric, column name or index for grouping
##'   multiple datasets. If 1 or ".id", treats as single dataset. Default is ".id"
##' @return A ggplot2 bar plot object showing the distribution of peaks relative
##'   to TSS. The plot uses a diverging bar chart format with upstream (5') on
##'   the left and downstream (3') on the right, with TSS at the center
##' @importFrom magrittr %<>%
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 scale_y_continuous
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 scale_fill_brewer
##' @importFrom ggplot2 scale_fill_hue
##' @importFrom ggplot2 scale_fill_manual
##' @importFrom ggplot2 geom_text
##' @importFrom rlang .data
##' @examples
##' \dontrun{
##' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
##' peakAnno <- annotatePeak(peakfile, TxDb=txdb)
##' plotDistToTSS(peakAnno)
##' }
##' @seealso \code{\link{annotatePeak}}
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
plotDistToTSS.data.frame <- function(peakDist,
                                     distanceColumn="distanceToTSS",
                                     distanceBreaks=c(0, 1000, 3000, 5000, 10000, 100000),
                                     palette = NULL,
                                     xlab="",
                                     ylab="Binding sites (%) (5'->3')",
                                     title="Distribution of transcription factor-binding loci relative to TSS",
                                     categoryColumn = ".id") {

    distanceBreaks = sort(distanceBreaks)
    hasZero = sum(distanceBreaks == 0)
    if (!hasZero) distanceBreaks = c(0, distanceBreaks)
    hasInf = sum(is.infinite(distanceBreaks))
    if (!hasInf) distanceBreaks = c(distanceBreaks, Inf)
    lbs = generate_break_lbs(distanceBreaks)
    peakDist$Feature = cut(abs(peakDist[[distanceColumn]]),
                           breaks = distanceBreaks,
                           labels = lbs,
                           include.lowest = TRUE)

    ## sign containing -1 and 1 for upstream and downstream
    peakDist$sign <- sign(peakDist[,distanceColumn])

    ## count frequencies
    if (categoryColumn == 1) {
      peakDist = peakDist |>
        summarise(freq = length(.data$Feature), .by = c("Feature", "sign")) |>
        mutate(freq = .data$freq/sum(.data$freq) * 100)
    } else {
      peakDist = peakDist |>
        summarise(freq = length(.data$Feature), .by = c(categoryColumn, "Feature", "sign")) |>
        mutate(freq = .data$freq/sum(.data$freq) * 100, .by = categoryColumn)
    }

    if (any(peakDist$sign == 0)) {
        zeroDist <- peakDist[peakDist$sign == 0,]
        zeroDist$freq <- zeroDist$freq/2
        zeroDist$sign <- -1
        peakDist[peakDist$sign == 0,] <- zeroDist
        zeroDist$sign <- 1
        peakDist <- rbind(peakDist, zeroDist)
    }

    if (categoryColumn == 1) {
        peakDist %<>% group_by(.data$Feature, .data$sign) %>%
            summarise(freq = sum(.data$freq))

        totalFreq <- peakDist %>% group_by(.data$sign) %>%
            summarise(total = sum(.data$freq))
    } else {
        peakDist %<>% group_by(.data$.id, .data$Feature, .data$sign) %>%
            summarise(freq = sum(.data$freq))
        totalFreq <- peakDist %>% group_by(.data$.id, .data$sign) %>%
            summarise(total = sum(.data$freq))
    }


    ## preparing ylim and y tick labels
    ds = max(totalFreq$total[totalFreq$sign == 1])
    dslim = ceiling(ds/10) * 10
    us = max(totalFreq$total[totalFreq$sign == -1])
    uslim = ceiling(us/10) * 10
    ybreaks <- seq(-uslim, dslim, by=10)
    ylbs <- abs(ybreaks)
    ylbs[ylbs == 0] <- "TSS"

    peakDist$Feature <- factor(peakDist$Feature, levels=rev(levels(peakDist$Feature)))
    if (categoryColumn == 1) {
        p <- ggplot(peakDist, aes(x=1, fill=.data$Feature))
    } else {
        p <- ggplot(peakDist, aes(x=.data[[categoryColumn]], fill=.data$Feature))
    }

    p <- p + geom_bar(data=subset(peakDist, sign==1), aes(y=.data$freq), stat="identity") +
        geom_bar(data=subset(peakDist, sign==-1), aes(y=-.data$freq), stat="identity")

    p <- p + geom_hline(yintercept = 0, colour = "black") +
        coord_flip() + theme_bw() +
            scale_y_continuous(breaks=ybreaks,labels=ylbs)

    p <- p + ylab(ylab) + xlab(xlab) + ggtitle(title)

    if (categoryColumn == 1) {
        p <- p + scale_x_continuous(breaks=NULL)
    }

    cols <- generate_colors(palette = palette, n = length(lbs))
    p <- p + scale_fill_manual(values=rev(cols), guide=guide_legend(reverse=TRUE))

    return(p)
}
