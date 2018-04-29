library(GenomicAlignments)

# Trims reads based on fixed offsets
trimBothEnds <- function(loffset, roffset, aln) {
  alnLength <- qwidth(aln)
  alnWidth <- alnLength - loffset - roffset
  alnGood <- alnWidth > 0
  alnSub <- aln[alnGood]
  if (length(alnSub) > 0) {
    alnStart <- ifelse(
      strand(alnSub) == "+",
      loffset + 1,
      roffset + 1)
    alnWidth <- alnWidth[alnGood]
     qnarrow(alnSub, start = alnStart, width = alnWidth)
  } else {
    GRanges()
  }
}
