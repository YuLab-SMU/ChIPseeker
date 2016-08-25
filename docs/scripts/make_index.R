#!/usr/local/bin/Rscript

library("scholar")
library("rmarkdown")

render("../private/doi_citation_badge.rmd", "md_document")

index_1 <- readLines("../private/index_part1.md")
index_2 <- readLines("../private/index_part2.md")
badge <- readLines("../private/doi_citation_badge.md")

index <- c(index_1, badge, "", index_2)

out <- file("../docs/index.md")
writeLines(index, out)
close(out)


