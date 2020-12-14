#' @title findORF
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param input.seq path to a folder of sequence alignments in phylip format.
#'
#' @param genbank.file contigs are added into existing alignment if algorithm is "add"
#'
#' @param codons algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param min.size TRUE to supress screen output
#'
#' @return an alignment of provided sequences in DNAStringSet format. Also can save alignment as a file with save.name
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#'
#' @export

findORF = function(input.seq = NULL,
                   codons = F,
                   min.size = 80) {

  #Sets up data
  input.seq = contig
  codon.table = data.frame(Type = as.character("0"),
                           Start = rep(0,6),
                           End = rep(0,6),
                           Frame = c("F1", "F2", "F3", "R1", "R2", "R3"))
  for.seq = as.character(input.seq)

  #Start codon
  ATG = matchPattern("ATG", for.seq)
  ATT = matchPattern("ATT", for.seq)

  #Stop codon
  TAA = matchPattern("TAA", for.seq)
  TAG = matchPattern("TAG", for.seq)
  AGA = matchPattern("AGA", for.seq)
  AGG = matchPattern("AGG", for.seq)

  #Forward Frame 1
  result.start1 = ATG[(ATG@ranges@start+2) %% 3 == 0]
  result.start2 = ATT[(ATT@ranges@start+2) %% 3 == 0]
  #Stop
  result.stop1 = TAA[(TAA@ranges@start+2) %% 3 == 0]
  result.stop2 = TAG[(TAG@ranges@start+2) %% 3 == 0]
  result.stop3 = AGA[(AGA@ranges@start+2) %% 3 == 0]
  result.stop4 = AGG[(AGG@ranges@start+2) %% 3 == 0]

  starts = c(result.start1@ranges@start, result.start2@ranges@start)
  ends = c(result.start1@ranges@start+2, result.start2@ranges@start+2)
  if (length(starts) != 0){
    codon.table = codon.table[codon.table$Frame != "F1",]
    temp.table = data.frame(Type = "Start", Start = starts, End = ends, Frame = "F1")
    codon.table = rbind(codon.table, temp.table)
  }

  #STop codons
  starts = c(result.stop1@ranges@start, result.stop2@ranges@start,
            result.stop3@ranges@start,  result.stop4@ranges@start)
  ends = c(result.stop1@ranges@start+2, result.stop2@ranges@start+2,
          result.stop3@ranges@start+2, result.stop4@ranges@start+2)
  if (length(starts) != 0){
    codon.table = codon.table[codon.table$Frame != "F1",]
    temp.table = data.frame(Type = "Stop", Start = starts, End = ends, Frame = "F1")
    codon.table = rbind(codon.table, temp.table)
  }




  #Forward Frame 2
  result1<-TAA[(TAA@ranges@start+1) %% 3 == 0]
  result2<-TGA[(TGA@ranges@start+1) %% 3 == 0]
  result3<-TAG[(TAG@ranges@start+1) %% 3 == 0]

  starts<-c(result1@ranges@start-1, result2@ranges@start-1, result3@ranges@start-1)
  ends<-c(result1@ranges@start+1, result2@ranges@start+1, result3@ranges@start+1)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "F2",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "F2")
    codon.table<-rbind(codon.table, temp.table)
  }

  #Forward Frame 3
  result1<-TAA[(TAA@ranges@start) %% 3 == 0]
  result2<-TGA[(TGA@ranges@start) %% 3 == 0]
  result3<-TAG[(TAG@ranges@start) %% 3 == 0]

  starts<-c(result1@ranges@start-2, result2@ranges@start-2, result3@ranges@start-2)
  ends<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "F3",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "F3")
    codon.table<-rbind(codon.table, temp.table)
  }

  #Sets up data
  rev.seq<-as.character(reverseComplement(input.seq))

  #Gets codon stuff
  TAA<-matchPattern("TAA", rev.seq)
  TGA<-matchPattern("TGA", rev.seq)
  TAG<-matchPattern("TAG", rev.seq)

  #Rev Frame 1
  result1<-TAA[(TAA@ranges@start+2) %% 3 == 0]
  result2<-TGA[(TGA@ranges@start+2) %% 3 == 0]
  result3<-TAG[(TAG@ranges@start+2) %% 3 == 0]

  starts<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  ends<-c(result1@ranges@start+2, result2@ranges@start+2, result3@ranges@start+2)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "R1",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "R1")
    codon.table<-rbind(codon.table, temp.table)
  }

  #Rev Frame 2
  result1<-TAA[(TAA@ranges@start+1) %% 3 == 0]
  result2<-TGA[(TGA@ranges@start+1) %% 3 == 0]
  result3<-TAG[(TAG@ranges@start+1) %% 3 == 0]

  starts<-c(result1@ranges@start-1, result2@ranges@start-1, result3@ranges@start-1)
  ends<-c(result1@ranges@start+1, result2@ranges@start+1, result3@ranges@start+1)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "R2",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "R2")
    codon.table<-rbind(codon.table, temp.table)
  }

  #Rev Frame 3
  result1<-TAA[(TAA@ranges@start) %% 3 == 0]
  result2<-TGA[(TGA@ranges@start) %% 3 == 0]
  result3<-TAG[(TAG@ranges@start) %% 3 == 0]

  starts<-c(result1@ranges@start-2, result2@ranges@start-2, result3@ranges@start-2)
  ends<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "R3",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "R3")
    codon.table<-rbind(codon.table, temp.table)
  } #end if

  if (codons == T) { return(codon.table) }

  if (codons == F) {
    frames<-unique(codon.table$Frame)
    orf.frame<-data.frame()
    for (x in 1:length(frames)){
      temp.codon<-codon.table[codon.table$Frame == frames[x],]
      temp.codon<-temp.codon[order(temp.codon$Start),]

      if (temp.codon$Start[1] == 0){
        temp.start<-as.numeric(gsub("F|R", "", temp.codon$Frame))
        add.frame<-data.frame(FrameStart = temp.start, FrameEnd = width(input.seq),
                              Size = (width(input.seq)-temp.start)+1, Frame = frames[x])
        orf.frame<-rbind(orf.frame, add.frame)
        next
      }
      #Goes through each of the given directions codons and converts to frame ranges
      temp.frame<-data.frame()
      for (y in 1:(nrow(temp.codon)+1)){
        #First y the start is 1, otherwise take from previous end
        if (y == 1){ frame.start<-as.numeric(gsub("F|R", "", temp.codon$Frame[y])) } else { frame.start<-temp.frame$FrameEnd[y-1]+4 }

        #Gets end by subtracting from the codon start
        frame.end<-temp.codon$Start[y]-1
        temp.frame<-rbind(temp.frame, data.frame(FrameStart = frame.start, FrameEnd = frame.end))
      } # end y loop

      temp.frame$FrameEnd[nrow(temp.frame)]<-width(input.seq)

      #Adds all the data together
      add.frame<-cbind(temp.frame, Size = (temp.frame$FrameEnd-temp.frame$FrameStart)+1, Frame = frames[x])
      orf.frame<-rbind(orf.frame, add.frame)

    } #end x loop

    orf.frame<-orf.frame[orf.frame$Size >= min.size,]
    return(orf.frame)
  } # end else

}# END FUNCTION


#### END SCRIPT
