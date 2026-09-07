
##' Plot peak coverage across chromosomes
##'
##' This function visualizes the coverage of ChIP-seq peaks across chromosomes,
##' showing the density and distribution of peaks along the genome.
##'
##' @description
##' The function creates a coverage plot showing where peaks are located across
##' chromosomes. It calculates coverage by counting overlapping peaks at each
##' genomic position, and displays the results as a bar/rectangular plot with
##' one panel per chromosome.
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Loads peaks from file or uses provided GRanges object
##'   \item Calculates coverage using \code{coverage()} function, optionally
##'         weighted by a metadata column
##'   \item Filters coverage by the \code{lower} threshold to remove low signals
##'   \item Aggregates coverage into contiguous regions
##'   \item Creates a ggplot2 visualization with:
##'     \itemize{
##'       \item One facet per chromosome (if multiple chromosomes)
##'       \item Rectangular bars showing coverage regions
##'       \item Color coding for multiple peak sets (if input is a list)
##'     }
##' }
##'
##' For multiple peak sets (list input), the function automatically assigns
##' colors and creates a grouped plot showing all sets together.
##'
##' @param peak peak file (BED format) or GRanges object, or a list of peak
##'   files/GRanges objects for comparing multiple datasets
##' @param weightCol character, name of a metadata column in the GRanges object
##'   to use as weights for coverage calculation. If NULL, all peaks have equal
##'   weight. Default is NULL
##' @param xlab character, label for x-axis. Default is "Chromosome Size (bp)"
##' @param ylab character, label for y-axis. Default is "" (empty)
##' @param title character, plot title. Default is "ChIP Peaks over Chromosomes"
##' @param chrs character vector, selected chromosomes to plot. If NULL, all
##'   chromosomes in the data are plotted. Default is NULL
##' @param xlim numeric vector of length 2, genomic range to plot (start, end).
##'   If NULL, the entire chromosome is shown. Note: this applies to all chromosomes
##'   if specified. Default is NULL
##' @param lower numeric, lower cutoff for coverage signal. Regions with coverage
##'   below this value are not displayed. Default is 1
##' @param fill_color character vector or single value. For single peak: a color
##'   name or hex code. For multiple peaks (list): either a vector of colors
##'   (one per peak set) or a palette name (e.g., "Set1", "Dark2") from
##'   RColorBrewer. The order of colors matches the order of peak sets. Default
##'   is "black"
##' @return A ggplot2 object that can be further customized or printed. The plot
##'   shows peak coverage as rectangular bars, with one panel per chromosome when
##'   multiple chromosomes are present
##' @import GenomeInfoDb
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_segment
##' @importFrom ggplot2 geom_blank
##' @importFrom ggplot2 geom_rect
##' @importFrom ggplot2 facet_grid
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 theme_classic
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 xlim
##' @importFrom ggplot2 ggtitle
##' @export
##' @author G Yu
covplot <- function(peak, weightCol=NULL,
                    xlab  = "Chromosome Size (bp)",
                    ylab  = "",
                    title = "ChIP Peaks over Chromosomes",
                    chrs  = NULL,
                    xlim  = NULL,
                    lower = 1,
                    fill_color = "black") {
    isList <- is.list(peak)
    if(!isList) {  # Note: don't support data.frame
        tm <- getChrCov(peak = peak, weightCol = weightCol, chrs = chrs, xlim = xlim, lower = lower)
    } else {
        ltm <- lapply(peak, getChrCov, weightCol = weightCol, chrs = chrs, xlim = xlim, lower = lower)
        if (is.null(names(ltm))) {
            nn <- paste0("peak", seq_along(ltm))
            warning("input is not a named list, set the name automatically to ", paste(nn, collapse = ' '))
            names(ltm) <- nn
        }
        tm <- dplyr::bind_rows(ltm, .id = ".id")
        chr.sorted <- sortChrName(as.character(unique(tm$chr)))
        tm$chr <- factor(tm$chr, levels = chr.sorted)
    }

    chr <- start <- end <- value <- .id <- NULL

    if(length(tm$chr) == 0){
        p <- ggplot(data.frame(x = 1)) + geom_blank()
    } else {
        p <- ggplot(tm, aes(start, value))

        ## p <- p + geom_segment(aes(x=start, y=0, xend=end, yend= value))
        if (isList) {
            if (length(fill_color) == length(peak) && all(is_valid_color(fill_color))){
                cols = fill_color
            } else {
                cols = generate_colors(fill_color, n = length(peak))
            }
            p <- p + geom_rect(aes(xmin = start, ymin = 0, xmax = end, ymax = value, fill = .id, color = .id)) +
                scale_color_manual(values = cols) +
                scale_fill_manual(values = cols)
        } else {
            p <- p + geom_rect(aes(xmin = start, ymin = 0, xmax = end, ymax = value), fill = fill_color, color = fill_color)
        }

        if(length(unique(tm$chr)) > 1) {
            p <- p + facet_grid(chr ~., scales="free")
        }

    }

    p <- p + theme_classic()
    p <- p + labs(x = xlab, y = ylab, title = title, fill = NULL, color = NULL)
    p <- p + scale_y_continuous(expand = c(0,0))
    p <- p + theme(strip.text.y=element_text(angle=360))
    p <- p + scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_si("")))

    if (!is.null(xlim) && !all(is.na(xlim)) && is.numeric(xlim) && length(xlim) == 2) {
        p <- p + xlim(xlim)
    }

    return(p)
}

##' @import S4Vectors IRanges
##' @importFrom dplyr group_by
##' @importFrom dplyr summarise
##' @importFrom magrittr %>%
##' @noRd
getChrCov <- function(peak, weightCol, chrs, xlim, lower=1) {
    if (is(peak, "GRanges")) {
        peak.gr <- peak
    } else if (file.exists(peak)) {
        peak.gr <- readPeakFile(peak, as="GRanges")
    } else {
        stop("peak should be a GRanges object or a peak file...")
    }

    if ( is.null(weightCol)) {
        peak.cov <- coverage(peak.gr)
    } else {
        weight <- mcols(peak.gr)[[weightCol]]
        peak.cov <- coverage(peak.gr, weight=weight)
    }

    cov <- lapply(peak.cov, IRanges::slice, lower=lower)

    get.runValue <- function(x) {
        y <- runValue(x)
        sapply(y@listData, mean)
        ## value <- x@subject@values
        ## value[value != 0]
    }

    chr <- start <- end <- cnt <- NULL

    ldf <- lapply(1:length(cov), function(i) {
        x <- cov[[i]]
        if (length(x@ranges) == 0) {
            msg <- paste0(names(cov[i]),
                          " dosen't contain signal higher than ",
                          lower)
            message(msg)
            return(NA)
        }
        data.frame(chr   = names(cov[i]),
                   start = start(x),
                   end   = end(x),
                   cnt   = get.runValue(x)
                                        # the following versions are more slower
                                        # unlist(runValue(x))
                                        # sapply(x, runValue)
                   )
    })

    ldf <- ldf[!is.na(ldf)]
    df <- do.call("rbind", ldf)

    chr.sorted <- sortChrName(as.character(unique(df$chr)))
    df$chr <- factor(df$chr, levels=chr.sorted)
    if (!is.null(chrs) && !all(is.na(chrs)) && all(chrs %in% chr.sorted)) {
        df <- df[df$chr %in% chrs, ]
    }
    if (!is.null(xlim) && !all(is.na(xlim)) && is.numeric(xlim) && length(xlim) == 2) {
        df <- df[df$start >= xlim[1] & df$end <= xlim[2],]
    }

    df2 <- group_by(df, chr, start, end) %>% summarise(value=sum(cnt), .groups = "drop")
    return(df2)
}

# a simple `stringr::str_sort(numeric=TRUE)` implementation
sortChrName <- function(chr.name, decreasing = FALSE) {
    ## universal sort function, support organisms other than human
    chr_part <- sub("^(\\D*)(\\d*)$", "\\1", chr.name)
    num_part <- as.numeric(sub("^(\\D*)(\\d*)$", "\\2", chr.name))
    chr.name[order(chr_part, num_part, decreasing = decreasing)]
}


