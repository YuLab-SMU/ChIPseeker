##' Filter GRanges objects using dplyr syntax
##'
##' This function extends \code{dplyr::filter()} to work with GRanges objects,
##' allowing users to filter genomic ranges using familiar dplyr syntax while
##' preserving the GRanges structure and all metadata columns.
##'
##' @description
##' The function filters GRanges objects by converting them to a data.frame,
##' applying the dplyr filter operation, and converting back to GRanges. This
##' allows filtering based on any metadata column or standard GRanges columns
##' (seqnames, start, end, width, strand).
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Converts GRanges to data.frame using \code{as.data.frame()}
##'   \item Applies \code{dplyr::filter()} with the provided conditions
##'   \item Removes unused factor levels using \code{droplevels()}
##'   \item Converts back to GRanges using \code{makeGRangesFromDataFrame()}
##'         with \code{keep.extra.columns=TRUE}
##' }
##'
##' All metadata columns are preserved in the output. The function supports
##' all standard dplyr filter operations including logical operators, comparisons,
##' and functions like \code{is.na()}, \code{between()}, etc.
##'
##' @param .data GRanges object to filter
##' @param ... logical expressions indicating which ranges to keep. Can use any
##'   column names from the GRanges metadata or standard GRanges columns
##'   (seqnames, start, end, width, strand). Multiple conditions are combined
##'   with AND logic
##' @param .by optional grouping variable (for dplyr compatibility). Currently
##'   not used but included for API consistency
##' @param .preserve logical, whether to preserve grouping structure (for dplyr
##'   compatibility). Currently not used but included for API consistency.
##'   Default is FALSE
##' @return GRanges object containing only the filtered ranges that satisfy
##'   the conditions. All metadata columns are preserved, and unused factor
##'   levels are dropped
##' @method filter GRanges
##' @importFrom dplyr filter
##' @export
##' @examples
##' \dontrun{
##' ## Filter peaks by annotation
##' filtered_peaks <- peak.gr %>% filter(annotation == "Promoter")
##'
##' ## Filter by distance to TSS
##' near_peaks <- peak.gr %>% filter(abs(distanceToTSS) < 5000)
##'
##' ## Multiple conditions (AND logic)
##' specific_peaks <- peak.gr %>%
##'   filter(annotation == "Promoter", abs(distanceToTSS) < 3000)
##'
##' ## Filter by chromosome
##' chr1_peaks <- peak.gr %>% filter(seqnames == "chr1")
##'
##' ## Complex conditions
##' complex_filter <- peak.gr %>%
##'   filter(annotation %in% c("Promoter", "5' UTR"),
##'          !is.na(geneId),
##'          width > 200)
##' }
filter.GRanges = function(.data, ..., .by = NULL, .preserve = FALSE) {
  dots = rlang::quos(...)
  as.data.frame(.data) |>
    dplyr::filter(!!!dots, .by = .by, .preserve = .preserve) |>
    droplevels() |>
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}

##' Mutate GRanges objects using dplyr syntax
##'
##' This function extends \code{dplyr::mutate()} to work with GRanges objects,
##' allowing users to add new metadata columns or modify existing ones using
##' familiar dplyr syntax while preserving the GRanges structure.
##'
##' @description
##' The function mutates GRanges objects by converting them to a data.frame,
##' applying the dplyr mutate operation, and converting back to GRanges. This
##' allows creating new columns or modifying existing metadata columns based on
##' expressions involving other columns.
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Converts GRanges to data.frame using \code{as.data.frame()}
##'   \item Validates that \code{.before} and \code{.after} are not both specified
##'   \item Applies \code{dplyr::mutate()} with the provided expressions,
##'         handling column positioning if \code{.before} or \code{.after} is specified
##'   \item Converts back to GRanges using \code{makeGRangesFromDataFrame()}
##'         with \code{keep.extra.columns=TRUE}
##' }
##'
##' All genomic range information (seqnames, ranges, strand) is preserved.
##' The function supports all standard dplyr mutate operations including
##' arithmetic, logical operations, and functions.
##'
##' @param .data GRanges object to mutate
##' @param ... name-value pairs of expressions. Use \code{name = value} to add
##'   new columns or modify existing ones. The name becomes the column name,
##'   and the value is an expression that can reference other columns
##' @param .by optional grouping variable (for dplyr compatibility). Currently
##'   not used but included for API consistency
##' @param .keep character, which columns to keep. Options: "all" (keep all
##'   columns), "used" (keep only columns used in expressions), "unused"
##'   (keep only columns not used in expressions), "none" (keep only new columns).
##'   Default is "all"
##' @param .before,.after column position specifiers. Place new columns before
##'   or after the specified column. Cannot use both. If NULL, new columns are
##'   appended at the end. Default is NULL for both
##' @return GRanges object with modified or new metadata columns. All genomic
##'   range information is preserved. The order and position of new columns
##'   follows the \code{.before} or \code{.after} specifications
##' @method mutate GRanges
##' @importFrom dplyr mutate
##' @export
##' @examples
##' \dontrun{
##' ## Add a new column based on existing data
##' peak.gr <- peak.gr %>% mutate(logScore = log10(score))
##'
##' ## Modify existing column
##' peak.gr <- peak.gr %>% mutate(distance = abs(distanceToTSS))
##'
##' ## Create multiple columns
##' peak.gr <- peak.gr %>%
##'   mutate(
##'     distance_kb = abs(distanceToTSS) / 1000,
##'     is_promoter = annotation == "Promoter",
##'     category = ifelse(abs(distanceToTSS) < 3000, "near", "far")
##'   )
##'
##' ## Place new column before specific column
##' peak.gr <- peak.gr %>%
##'   mutate(new_col = 1, .before = annotation)
##'
##' ## Keep only used columns
##' peak.gr <- peak.gr %>%
##'   mutate(result = score * 2, .keep = "used")
##' }
mutate.GRanges = function(.data, ..., .by = NULL,
                           .keep = c("all", "used", "unused", "none"),
                           .before = NULL,
                           .after = NULL) {
  dots = rlang::quos(...)
  df = as.data.frame(.data)

  if (!is.null(.before) && !is.null(.after)) {
    stop("You can't supply both `.before` and `.after`.")
  }

  if (!is.null(.before)) {
    df = df |>
      dplyr::mutate(!!!dots, .by = .by, .keep = .keep, .before = .before)
  } else if (!is.null(.after)) {
    df = df |>
      dplyr::mutate(!!!dots, .by = .by, .keep = .keep, .after = .after)
  } else {
    df = df |> dplyr::mutate(!!!dots, .by = .by, .keep = .keep)
  }

  df |>
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}

##' Rename columns in GRanges objects using dplyr syntax
##'
##' This function extends \code{dplyr::rename()} to work with GRanges objects,
##' allowing users to rename metadata columns using familiar dplyr syntax while
##' preserving the GRanges structure.
##'
##' @description
##' The function renames metadata columns in GRanges objects by converting them
##' to a data.frame, applying the dplyr rename operation, and converting back
##' to GRanges. This allows renaming any metadata column while preserving all
##' genomic range information.
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Converts GRanges to data.frame using \code{as.data.frame()}
##'   \item Applies \code{dplyr::rename()} with the provided name mappings
##'   \item Converts back to GRanges using \code{makeGRangesFromDataFrame()}
##'         with \code{keep.extra.columns=TRUE}
##' }
##'
##' All genomic range information (seqnames, ranges, strand) is preserved.
##' Only metadata columns are renamed; standard GRanges columns (seqnames, start,
##' end, width, strand) cannot be renamed through this function.
##'
##' @param x GRanges object to rename columns in
##' @param ... name-value pairs where the name is the new column name and the
##'   value (unquoted) is the old column name. Syntax: \code{new_name = old_name}.
##'   Multiple renames can be specified in a single call
##' @return GRanges object with renamed metadata columns. All genomic range
##'   information and other metadata columns are preserved unchanged
##' @method rename GRanges
##' @importFrom rlang quos
##' @export
##' @examples
##' \dontrun{
##' ## Rename a single column
##' peak.gr <- peak.gr %>% rename(geneSymbol = SYMBOL)
##'
##' ## Rename multiple columns
##' peak.gr <- peak.gr %>%
##'   rename(
##'     geneSymbol = SYMBOL,
##'     geneName = GENENAME,
##'     distance = distanceToTSS
##'   )
##'
##' ## Rename to remove spaces or special characters
##' peak.gr <- peak.gr %>% rename(gene_id = "gene Id")
##' }
rename.GRanges = function(x, ...){
  dots = rlang::quos(...)
  as.data.frame(x) |>
    dplyr::rename(!!!dots) |>
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}

##' Arrange (sort) GRanges objects using dplyr syntax
##'
##' This function extends \code{dplyr::arrange()} to work with GRanges objects,
##' allowing users to sort genomic ranges by any column using familiar dplyr
##' syntax while preserving the GRanges structure.
##'
##' @description
##' The function sorts GRanges objects by converting them to a data.frame,
##' applying the dplyr arrange operation, and converting back to GRanges. This
##' allows sorting by any metadata column or standard GRanges columns (seqnames,
##' start, end, width, strand).
##'
##' @details
##' The function performs the following steps:
##' \enumerate{
##'   \item Converts GRanges to data.frame using \code{as.data.frame()}
##'   \item Applies \code{dplyr::arrange()} with the specified sorting columns
##'   \item Converts back to GRanges using \code{makeGRangesFromDataFrame()}
##'         with \code{keep.extra.columns=TRUE}
##' }
##'
##' All genomic range information and metadata columns are preserved. Sorting
##' is stable (preserves original order for ties). Multiple sorting columns are
##' applied in order (first column is primary sort key).
##'
##' @param .data GRanges object to arrange
##' @param ... column names or expressions to sort by. Multiple columns can be
##'   specified, with earlier columns taking precedence. Use \code{desc()} to
##'   sort in descending order. Can use any metadata column or standard GRanges
##'   columns (seqnames, start, end, width, strand)
##' @param .by_group logical, if TRUE, will sort first by grouping variable
##'   (if the data is grouped). Currently not used but included for API
##'   consistency. Default is FALSE
##' @return GRanges object with ranges sorted by the specified columns. All
##'   metadata columns are preserved. The sort order is stable (ties preserve
##'   original order)
##' @method arrange GRanges
##' @importFrom dplyr arrange
##' @export
##' @examples
##' \dontrun{
##' ## Sort by chromosome and start position (genomic order)
##' peak.gr <- peak.gr %>% arrange(seqnames, start)
##'
##' ## Sort by distance (descending - closest to TSS first)
##' peak.gr <- peak.gr %>% arrange(desc(abs(distanceToTSS)))
##'
##' ## Sort by annotation, then by distance
##' peak.gr <- peak.gr %>% arrange(annotation, distanceToTSS)
##'
##' ## Sort by score (descending), then by chromosome
##' peak.gr <- peak.gr %>% arrange(desc(score), seqnames)
##'
##' ## Complex sorting with multiple criteria
##' peak.gr <- peak.gr %>%
##'   arrange(
##'     annotation,
##'     desc(score),
##'     abs(distanceToTSS)
##'   )
##' }
arrange.GRanges = function(.data, ..., .by_group = FALSE){
  dots = rlang::quos(...)
  as.data.frame(.data) |>
    dplyr::arrange(!!!dots, .by_group = .by_group) |>
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}


