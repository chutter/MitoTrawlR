##########################################################################################3
##########################################################################################3
### Preprocessing using PhyloCap if needed

#Installs updated package version
devtools::install_github("chutter/PhyloCap", upgrade = "never", force = TRUE)
library(PhyloCap)

source("/Volumes/Rodents/Mitogenomes/Crocidura/configuration-file-preproc.R")

#Checks if everything is installed
pass.fail = setupCheck(anaconda.environment =  NULL,
                       fastp.path = fastp.path,
                       samtools.path = samtools.path,
                       bwa.path = bwa.path,
                       spades.path = spades.path,
                       bbmap.path = bbmap.path,
                       blast.path = blast.path,
                       mafft.path = mafft.path,
                       iqtree.path = iqtree.path,
                       trimAl.path = trimAl.path,
                       julia.path = julia.path,
                       taper.path = taper.path)

if (pass.fail == FALSE){ stop("Some required programs are missing") } else {
  print("all required programs are found, PhyloCap pipeline continuing...")
}


setwd(work.dir)
dir.create("processed-reads")

if (organize.reads == TRUE) {
  organizeReads(read.directory = read.dir,
                output.dir = paste0(processed.reads, "/organized-reads"),
                rename.file = file.rename,
                overwrite = overwrite)
  input.reads = paste0(processed.reads, "/organized-reads")
} else {input.reads = read.dir }

if (remove.adaptors == TRUE) {
  removeAdaptors(input.reads = input.reads,
                 output.directory = paste0(processed.reads, "/adaptor-removed-reads"),
                 fastp.path = fastp.path,
                 threads = threads,
                 mem = memory,
                 resume = resume,
                 overwrite = overwrite,
                 quiet = quiet)
  input.reads = paste0(processed.reads, "/adaptor-removed-reads")
}

#Runs decontamination of reads
if (decontamination == TRUE){
  #Creates the database by downloading
  createContaminantDB(decontamination.list = contaminant.genome.list,
                      output.directory = "contaminant-references",
                      include.human = include.human,
                      include.univec = include.univec,
                      overwrite = overwrite)

  ## remove external contamination
  removeContamination(input.reads = input.reads,
                      output.directory = paste0(processed.reads, "/decontaminated-reads"),
                      decontamination.path = "contaminant-references",
                      map.match = decontamination.match,
                      samtools.path = samtools.path,
                      bwa.path = bwa.path,
                      threads = threads,
                      mem = memory,
                      resume = resume,
                      overwrite = overwrite,
                      overwrite.reference = overwrite,
                      quiet = quiet)
  input.reads = paste0(processed.reads, "/decontaminated-reads")
}


if (merge.pe.reads == TRUE){
  #merge paired end reads
  mergePairedEndReads(input.reads = input.reads,
                      output.directory =  paste0(processed.reads, "/pe-merged-reads"),
                      fastp.path = fastp.path,
                      threads = threads,
                      mem = memory,
                      resume = resume,
                      overwrite = overwrite,
                      quiet = quiet)
  input.reads = paste0(processed.reads, "/pe-merged-reads")
} #end decontamination



##########################################################################################################
#Step 2: Mitogenome assembly
##########################################################################################################

#Installs updated package version
devtools::install_github("chutter/MitoCap", upgrade = "never", force = TRUE)
library(MitoCap)
library(foreach)

source("/Volumes/Rodents/Mitogenomes/Crocidura/configuration-file-mitocap.R")
setwd(work.dir)

#Checks if everything is installed
pass.fail = MitoCap::setupCheck(anaconda.environment =  NULL,
                                samtools.path = samtools.path,
                                bwa.path = bwa.path,
                                spades.path = spades.path,
                                bbmap.path = bbmap.path,
                                blast.path = blast.path,
                                mafft.path = mafft.path,
                                iqtree.path = iqtree.path,
                                trimAl.path = trimAl.path,
                                trnascan.path = trnascan.path)

if (pass.fail == FALSE){ stop("Some required programs are missing") } else {
  print("all required programs are found, PhyloCap pipeline continuing...")
}

# #Makes your reference, be sure to inspect the table it makes to ensure accuracy
buildReference(reference.fasta = reference.fasta,
               annotation.file = annotation.file,
               annotation.type = annotation.type,
               reference.name = "reference",
               overwrite = overwrite,
               rep.origin = FALSE)

# #Iteratively assembles to reference
mitochondrialCapture(input.reads = input.reads,
                     reference.name = "reference",
                     output.dir = "draftContigs",
                     min.iterations = 5,
                     max.iterations = 30,
                     min.length = 16000,
                     max.length = 40000,
                     min.ref.id = 0.8,
                     memory = memory,
                     threads = threads,
                     spades.path = spades.path,
                     bbmap.path = bbmap.path,
                     resume = TRUE,
                     overwrite = FALSE)

### TO DO: Add a redo option for the non-continuous ones, using itself as a reference

#Annotates mitochondrial contigs
#To do: output format gff and others
#To do 2: save new contigs (instead of old) but with corrected start point so half of something isn't lost
#Combine build mtgenomes with this step. difficut.

annotateMitoContigs(contig.folder = "newContigs",
                    reference.name = "reference",
                    blast.path = "blast",
                    trnascan.path = trnascan.path,
                    organism.type = "vertebrate",
                    overwrite = TRUE,
                    quiet = FALSE)

#Aligns all the different markers
markerAlignment(input.folder = "Annotations/sample-markers",
                reference.name = "reference",
                threads = threads,
                overwrite = TRUE)

#Trims the alignments to ready for concatenation or gene tree estimation
trimMitoAlignments(alignment.dir = "Alignments/untrimmed-alignments",
                   alignment.format = "phylip",
                   output.dir = "Alignments/trimmed-alignments",
                   output.format = "phylip",
                   sample.similiarity = TRUE,
                   TrimAl = TRUE,
                   TrimAl.path = "trimal",
                   trim.external = TRUE,
                   min.external.percent = 50,
                   trim.coverage = TRUE,
                   min.coverage.percent = 30,
                   trim.column = TRUE,
                   min.column.gap.percent = 100,
                   alignment.assess = TRUE,
                   min.sample.bp = 10,
                   min.align.length = 0,
                   min.taxa.count = 12,
                   min.gap.percent = 50,
                   overwrite = TRUE)

#Aligns the mitogenomes and outputs summary stats
alignMitogenomes(alignment.folder = "Alignments/untrimmed-alignments",
                 genbank.file = gb.file,
                 draft.contigs = "draftContigs",
                 output.dir = "Genomes",
                 dataset.name = "untrimmed",
                 overwrite = TRUE)

#Aligns the mitogenomes and outputs summary stats
alignMitogenomes(alignment.folder = "Alignments/trimmed-alignments",
                 genbank.file = gb.file,
                 draft.contigs = "draftContigs",
                 output.dir = "Genomes",
                 dataset.name = "trimmed",
                 overwrite = FALSE)


# tree.dir = "/Volumes/Rodents/Murinae/Mitochondrial_genomes/trees"
# dir.create(tree.dir)
# setwd(tree.dir)
#
# system(paste0("iqtree2 -p ", work.dir,"/Alignments/trimmed-alignments --prefix murinae -bb 1000 -nt 6",
#               " -m GTR"))
#
#
# align.dir = paste0(work.dir, "/Alignments/trimmed-alignments")
#
# align.files = list.files(align.dir)
#
# ##########################################
# #### Plot tree with bars
# ##########################################
#
# library(ggtree)
# library(ggplot2)
#
# tree.file = "/Volumes/Rodents/Murinae/Mitochondrial_genomes/trees/murinae.contree"
# all.tree = ape::read.tree(tree.file)
#
# save.list = c()
# for (i in 1:length(align.files)){
#
#   align = Biostrings::DNAStringSet(Biostrings::readAAMultipleAlignment(file = paste0(align.dir, "/", align.files[i]), format = "phylip"))
#
#   #Remove gap only alignments
#   c.align = strsplit(as.character(align), "")
#   gap.align = lapply(c.align, function(x) gsub("N|n", "-", x) )
#   base.count = unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
#   #save.list = append(save.list, base.count)
#
#   save.list[[i]] = base.count
#
#  # tapply(unlist(L), names(unlist(L)), sum)
#
# }
#
# sum.samples = tapply(unlist(save.list), names(unlist(save.list)), sum)
# per.samples = sum.samples/max(sum.samples)
#
# ### Map data to tree
#
# all.tree = ape::root(all.tree, outgroup = c("Lophuromys_woosnami_LSUMZ37793", "Lophiomys_imhausi_UM5152"))
# all.data = per.samples[pmatch(names(per.samples), all.tree$tip.label)]
# all.data = per.samples[pmatch(all.tree$tip.label, names(per.samples))]
# all.data = data.frame(id = names(all.data), val = all.data)
#
# #Plots tree
# #phytools::plotTree.barplot(all.tree, per.samples, cex = 0.5)
#
# p <- ggtree::ggtree(all.tree, ladderize = T, right = T) + geom_tiplab(align=TRUE, linesize=.35)
# p2 <- ggtree::facet_plot(p, panel='bar', data=all.data, geom=geom_segment,
#                  aes(x=0, xend=val, y=y, yend=y), size=0.75, color='blue4')
# p2
#
#
# ### Working here
# #Builds and standardizes final mitogenomes
# #Reverse the whole thing if its backwards
# # Orient each one to begin with 1
# # Remove over-assembled ones next
#
# buildMitogenomes(annotation.dir = "Annotations",
#                  alignment.folder = "Alignments/untrimmed-alignments",
#                  genome.alignment = "Genomes/alignments/untrimmed_mitogenome_alignment.phy",
#                  genome.dir = "Genomes",
#                  output.dir = "untrimmed-finished",
#                  overwrite = FALSE)
#
#

