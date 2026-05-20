#' @title findORF
#'
#' @description Finds open reading frames (ORFs) or codon positions in all six
#'   reading frames (three forward, three reverse) of a DNA sequence using the
#'   vertebrate mitochondrial genetic code. Stop codons are TAA, TAG, AGA, and
#'   AGG; TGA codes for Tryptophan and is not treated as a stop codon. When
#'   \code{codons = TRUE}, returns a data frame of all start and stop codon
#'   positions. When \code{codons = FALSE} (default), returns a data frame of
#'   ORF ranges (start, end, size, frame) filtered to \code{min.size} codons.
#'
#' @param input.seq a \code{DNAString} or single-element \code{DNAStringSet}
#'   to search for ORFs.
#'
#' @param codons logical; if TRUE, return all codon positions rather than
#'   assembled ORF ranges.
#'
#' @param min.size minimum ORF size in nucleotides to report when
#'   \code{codons = FALSE}.
#'
#' @return A data frame. When \code{codons = FALSE}: columns \code{FrameStart},
#'   \code{FrameEnd}, \code{Size}, and \code{Frame}. When \code{codons = TRUE}:
#'   columns \code{Type}, \code{Start}, \code{End}, and \code{Frame}.
#'
#' @examples
#' \dontrun{
#' library(Biostrings)
#' seq <- DNAString("ATGAAATTTGGGTAGCCCATGCAGTAA")
#' findORF(input.seq = seq, codons = FALSE, min.size = 9)
#' }
#'
#' @export

findORF = function(input.seq = NULL,
                   codons = FALSE,
                   min.size = 80) {

  codon.table = data.frame(Type  = as.character("0"),
                           Start = rep(0, 6),
                           End   = rep(0, 6),
                           Frame = c("F1", "F2", "F3", "R1", "R2", "R3"),
                           stringsAsFactors = FALSE)

  for.seq = as.character(input.seq)

  # Start codons (vertebrate mito: ATG, ATT)
  ATG = Biostrings::matchPattern("ATG", for.seq)
  ATT = Biostrings::matchPattern("ATT", for.seq)

  # Stop codons — vertebrate mitochondrial genetic code:
  # TAA, TAG, AGA, AGG are stop codons; TGA codes for Trp (NOT a stop codon)
  TAA = Biostrings::matchPattern("TAA", for.seq)
  TAG = Biostrings::matchPattern("TAG", for.seq)
  AGA = Biostrings::matchPattern("AGA", for.seq)
  AGG = Biostrings::matchPattern("AGG", for.seq)

  ##############
  # Forward Frame 1 (positions where start+2 ≡ 0 mod 3)
  ##############
  f1.starts = c()
  f1.ends   = c()
  f1.types  = c()

  s1 = ATG[(ATG@ranges@start + 2) %% 3 == 0]
  s2 = ATT[(ATT@ranges@start + 2) %% 3 == 0]
  if (length(s1) + length(s2) > 0) {
    sv = c(s1@ranges@start, s2@ranges@start)
    f1.starts = c(f1.starts, sv)
    f1.ends   = c(f1.ends,   sv + 2)
    f1.types  = c(f1.types,  rep("Start", length(sv)))
  }

  p1 = TAA[(TAA@ranges@start + 2) %% 3 == 0]
  p2 = TAG[(TAG@ranges@start + 2) %% 3 == 0]
  p3 = AGA[(AGA@ranges@start + 2) %% 3 == 0]
  p4 = AGG[(AGG@ranges@start + 2) %% 3 == 0]
  if (length(p1) + length(p2) + length(p3) + length(p4) > 0) {
    sv = c(p1@ranges@start, p2@ranges@start, p3@ranges@start, p4@ranges@start)
    f1.starts = c(f1.starts, sv)
    f1.ends   = c(f1.ends,   sv + 2)
    f1.types  = c(f1.types,  rep("Stop", length(sv)))
  }

  if (length(f1.starts) > 0) {
    codon.table = codon.table[codon.table$Frame != "F1",]
    codon.table = rbind(codon.table,
                        data.frame(Type = f1.types, Start = f1.starts, End = f1.ends,
                                   Frame = "F1", stringsAsFactors = FALSE))
  }

  ##############
  # Forward Frame 2 (positions where start+1 ≡ 0 mod 3)
  ##############
  p1 = TAA[(TAA@ranges@start + 1) %% 3 == 0]
  p2 = TAG[(TAG@ranges@start + 1) %% 3 == 0]
  p3 = AGA[(AGA@ranges@start + 1) %% 3 == 0]
  p4 = AGG[(AGG@ranges@start + 1) %% 3 == 0]
  starts = c(p1@ranges@start - 1, p2@ranges@start - 1, p3@ranges@start - 1, p4@ranges@start - 1)
  ends   = c(p1@ranges@start + 1, p2@ranges@start + 1, p3@ranges@start + 1, p4@ranges@start + 1)
  if (length(starts) > 0) {
    codon.table = codon.table[codon.table$Frame != "F2",]
    codon.table = rbind(codon.table,
                        data.frame(Type = "Stop", Start = starts, End = ends,
                                   Frame = "F2", stringsAsFactors = FALSE))
  }

  ##############
  # Forward Frame 3 (positions where start ≡ 0 mod 3)
  ##############
  p1 = TAA[(TAA@ranges@start) %% 3 == 0]
  p2 = TAG[(TAG@ranges@start) %% 3 == 0]
  p3 = AGA[(AGA@ranges@start) %% 3 == 0]
  p4 = AGG[(AGG@ranges@start) %% 3 == 0]
  starts = c(p1@ranges@start - 2, p2@ranges@start - 2, p3@ranges@start - 2, p4@ranges@start - 2)
  ends   = c(p1@ranges@start,     p2@ranges@start,     p3@ranges@start,     p4@ranges@start)
  if (length(starts) > 0) {
    codon.table = codon.table[codon.table$Frame != "F3",]
    codon.table = rbind(codon.table,
                        data.frame(Type = "Stop", Start = starts, End = ends,
                                   Frame = "F3", stringsAsFactors = FALSE))
  }

  ##############
  # Reverse complement — stop codons read 5'→3' on the minus strand
  ##############
  rev.seq = as.character(Biostrings::reverseComplement(input.seq))

  TAA.r = Biostrings::matchPattern("TAA", rev.seq)
  TAG.r = Biostrings::matchPattern("TAG", rev.seq)
  AGA.r = Biostrings::matchPattern("AGA", rev.seq)
  AGG.r = Biostrings::matchPattern("AGG", rev.seq)

  ##############
  # Reverse Frame 1
  ##############
  p1 = TAA.r[(TAA.r@ranges@start + 2) %% 3 == 0]
  p2 = TAG.r[(TAG.r@ranges@start + 2) %% 3 == 0]
  p3 = AGA.r[(AGA.r@ranges@start + 2) %% 3 == 0]
  p4 = AGG.r[(AGG.r@ranges@start + 2) %% 3 == 0]
  starts = c(p1@ranges@start,     p2@ranges@start,     p3@ranges@start,     p4@ranges@start)
  ends   = c(p1@ranges@start + 2, p2@ranges@start + 2, p3@ranges@start + 2, p4@ranges@start + 2)
  if (length(starts) > 0) {
    codon.table = codon.table[codon.table$Frame != "R1",]
    codon.table = rbind(codon.table,
                        data.frame(Type = "Stop", Start = starts, End = ends,
                                   Frame = "R1", stringsAsFactors = FALSE))
  }

  ##############
  # Reverse Frame 2
  ##############
  p1 = TAA.r[(TAA.r@ranges@start + 1) %% 3 == 0]
  p2 = TAG.r[(TAG.r@ranges@start + 1) %% 3 == 0]
  p3 = AGA.r[(AGA.r@ranges@start + 1) %% 3 == 0]
  p4 = AGG.r[(AGG.r@ranges@start + 1) %% 3 == 0]
  starts = c(p1@ranges@start - 1, p2@ranges@start - 1, p3@ranges@start - 1, p4@ranges@start - 1)
  ends   = c(p1@ranges@start + 1, p2@ranges@start + 1, p3@ranges@start + 1, p4@ranges@start + 1)
  if (length(starts) > 0) {
    codon.table = codon.table[codon.table$Frame != "R2",]
    codon.table = rbind(codon.table,
                        data.frame(Type = "Stop", Start = starts, End = ends,
                                   Frame = "R2", stringsAsFactors = FALSE))
  }

  ##############
  # Reverse Frame 3
  ##############
  p1 = TAA.r[(TAA.r@ranges@start) %% 3 == 0]
  p2 = TAG.r[(TAG.r@ranges@start) %% 3 == 0]
  p3 = AGA.r[(AGA.r@ranges@start) %% 3 == 0]
  p4 = AGG.r[(AGG.r@ranges@start) %% 3 == 0]
  starts = c(p1@ranges@start - 2, p2@ranges@start - 2, p3@ranges@start - 2, p4@ranges@start - 2)
  ends   = c(p1@ranges@start,     p2@ranges@start,     p3@ranges@start,     p4@ranges@start)
  if (length(starts) > 0) {
    codon.table = codon.table[codon.table$Frame != "R3",]
    codon.table = rbind(codon.table,
                        data.frame(Type = "Stop", Start = starts, End = ends,
                                   Frame = "R3", stringsAsFactors = FALSE))
  }

  if (codons) { return(codon.table) }

  ##############
  # Build ORF ranges from stop codon positions
  ##############
  frames    = unique(codon.table$Frame)
  orf.frame = data.frame()

  for (x in seq_along(frames)){
    temp.codon = codon.table[codon.table$Frame == frames[x],]
    temp.codon = temp.codon[order(temp.codon$Start),]

    if (temp.codon$Start[1] == 0){
      temp.start = as.numeric(gsub("F|R", "", temp.codon$Frame[1]))
      add.frame = data.frame(FrameStart = temp.start,
                              FrameEnd   = Biostrings::width(input.seq),
                              Size       = (Biostrings::width(input.seq) - temp.start) + 1,
                              Frame      = frames[x])
      orf.frame = rbind(orf.frame, add.frame)
      next
    }

    temp.frame = data.frame()
    for (y in seq_len(nrow(temp.codon))){
      if (y == 1) {
        frame.start = as.numeric(gsub("F|R", "", temp.codon$Frame[1]))
      } else {
        frame.start = temp.frame$FrameEnd[y - 1] + 4
      }
      frame.end  = temp.codon$Start[y] - 1
      temp.frame = rbind(temp.frame, data.frame(FrameStart = frame.start, FrameEnd = frame.end))
    }

    # Add the final ORF from the last stop codon to the end of the sequence
    last.start = temp.frame$FrameEnd[nrow(temp.frame)] + 4
    temp.frame = rbind(temp.frame,
                       data.frame(FrameStart = last.start,
                                  FrameEnd   = as.integer(Biostrings::width(input.seq))))

    add.frame = cbind(temp.frame,
                      Size  = (temp.frame$FrameEnd - temp.frame$FrameStart) + 1,
                      Frame = frames[x])
    orf.frame = rbind(orf.frame, add.frame)
  }

  orf.frame = orf.frame[orf.frame$Size >= min.size,]
  return(orf.frame)

}#end function
