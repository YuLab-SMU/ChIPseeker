#!/usr/local/bin/Rscript

library("scholar")
library("rCharts")
library("methods")
library("rmarkdown")

featured_article <- function(pub) {
    id <- "DO5oG40AAAAJ"

    d <- get_publications(id)
    i <- grep(pub, d$title)
    pubid <- d$pubid[i]

    df <- get_article_cite_history(id, pubid)
    g <- mPlot(cites ~ year, data = df, type = "Bar")

    ff <- tempfile(fileext=".html")
    sink(ff)
    g$print(include_assets = T)
    sink()

    

    x <- readLines(ff)
    i = grep("</style>$", x)
    x <- x[-(1:i)]
    
    header <- readLines("../private/featuredArticles_header.md")
    citation <- readLines("../private/featuredArticles_citation.md")

    render("../private/citation_badge.rmd", "md_document")

    badge <- readLines("../private/citation_badge.md")
    unlink("../private/citation_badge.md")
    
    y <- c(header, badge, "", x, citation)
    
    out <- file("../docs/featuredArticles.md", "w")
    writeLines(y, out)
    close(out)
}


featured_article("ChIPseeker")



