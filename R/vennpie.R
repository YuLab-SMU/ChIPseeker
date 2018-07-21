##' @importFrom plotrix floating.pie
vennpie.csAnno <- function(x, r=0.2) {
    detailGenomicAnnotation <- x@detailGenomicAnnotation

    distance <- as.data.frame(x)$distanceToTSS
    total <- nrow(detailGenomicAnnotation)
    Genic <- sum(detailGenomicAnnotation$genic)

    Intergenic <- total-Genic
    Distal_Intergenic <- sum(detailGenomicAnnotation$distal_intergenic)
    Intron <- sum(detailGenomicAnnotation$Intron)
    Exon <- sum(detailGenomicAnnotation$Exon)
    Upstream <- sum(detailGenomicAnnotation$Promoter & distance < 0)

    ## fiveUTR <- sum(detailGenomicAnnotation$fiveUTR)
    ## threeUTR <- sum(detailGenomicAnnotation$threeUTR)
    Downstream <- sum(detailGenomicAnnotation$downstream)

    ## fiveUTR='#e5f5e0',threeUTR='#a1d99b',
    cols <- c(NO='white',Genic='#3182bd',Intergenic='#fec44f',
              Intron='#fc9272',Exon='#9ecae1', Upstream='#ffeda0',
              Downstream='#fee0d2',Distal_Intergenic='#d95f0e')

    ##par(mai = c(0,0,0,0))
    ##layout(matrix(c(1,2), ncol=2), widths=c(0.7,0.3))
    pie(1, radius=r, init.angle=90, col="white", border=NA, labels='')

    ## https://www.biostars.org/p/326456/
    ## if count is 0, floating pie will ignore it
    ## and the color will mismatch with the category
    ## fixed by adding pseudo-count +1
    floating.pie(0,0, c(Exon,
                        Genic-Exon,
                        Distal_Intergenic,
                        Downstream,
                        Intergenic-Distal_Intergenic-Downstream
                        ) + 1,
                 radius=4*r,
                 startpos=pi/2,
                 col=cols[c("Exon", "NO", "NO", "Downstream", "NO")],
                 border=NA)

    floating.pie(0,0, c(Genic-Intron,
                        Intron,
                        Distal_Intergenic,
                        Intergenic-Upstream-Distal_Intergenic,
                        Upstream) +1 ,
                 radius=3*r,
                 startpos=pi/2,
                 col=cols[c("NO", "Intron", "Distal_Intergenic",
                     "NO", "Upstream")],
                 border=NA)

    floating.pie(0, 0, c(Genic, Intergenic) +1,
                 radius=2*r,
                 startpos=pi/2,
                 col=cols[c("Genic", "Intergenic")],
                 border=NA)
    ##plot.new()
    ##legend(center), legend=names(cols)[-1], fill=cols[-1], bty="n")
    legend(3*r, 3*r, legend=sub("_", " ", names(cols)[-1]),
           fill=cols[-1], bty="n", cex=1.2)
}
