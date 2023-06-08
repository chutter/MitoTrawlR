##########################################################################################################
#Step 1: Mitogenome assembly
##########################################################################################################

#Installs updated package version
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
library(PhyloProcessR)

#Installs updated package version
devtools::install_github("chutter/MitoTrawlR", upgrade = "never", force = TRUE, dependencies = FALSE)
library(MitoTrawlR)
library(foreach)

source("configuration-file.R")
setwd(working.directory)

# #Checks if everything is installed
# pass.fail = MitoCap::setupCheck(anaconda.environment =  NULL,
#                                 samtools.path = samtools.path,
#                                 bwa.path = bwa.path,
#                                 spades.path = spades.path,
#                                 bbmap.path = bbmap.path,
#                                 blast.path = blast.path,
#                                 mafft.path = mafft.path,
#                                 iqtree.path = iqtree.path,
#                                 trimAl.path = trimAl.path,
#                                 tRNAscan.path = tRNAscan.path)

# if (pass.fail == FALSE){ stop("Some required programs are missing") } else {
#   print("all required programs are found, MitoCap pipeline continuing...")
# }

# #Makes your reference, be sure to inspect the table it makes to ensure accuracy
buildReference(reference.fasta = reference.fasta,
               annotation.file = annotation.file,
               annotation.type = annotation.type,
               reference.name = "reference",
               overwrite = overwrite,
               rep.origin = FALSE)

# #Iteratively assembles to reference
mitochondrialCapture(input.reads = read.directory,
                     reference.name = "reference",
                     output.dir = "draftContigs",
                     min.iterations = min.iterations,
                     max.iterations = max.iterations,
                     min.length = min.length,
                     max.length = max.length,
                     min.ref.id = min.read.match,
                     memory = memory,
                     threads = threads,
                     spades.path = spades.path,
                     bbmap.path = bbmap.path,
                     cap3.path = cap3.path,
                     blast.path = blast.path,
                     overwrite = overwrite)

### TO DO: Add a redo option for the non-continuous ones, using itself as a reference

#Annotates mitochondrial contigs
#To do: output format gff and others
#To do 2: save new contigs (instead of old) but with corrected start point so half of something isn't lost
#Combine build mtgenomes with this step. difficut.

annotateMitoContigs(contig.folder = "draftContigs",
                    reference.name = "reference",
                    blast.path = blast.path,
                    tRNAscan.path = tRNAscan.path,
                    organism.type = "vertebrate",
                    overwrite = overwrite,
                    quiet = quiet)

#Aligns all the different markers
markerAlignment(input.folder = "Annotations/sample-markers",
                reference.name = "reference",
                threads = threads,
                mafft.path = mafft.path,
                overwrite = overwrite)

#Fix the installs for this
trimMitoAlignments(alignment.dir = "Alignments/untrimmed-alignments",
                    alignment.format = "phylip",
                    output.dir = "Alignments/trimmed-alignments",
                    output.format = "phylip",
                    overwrite = overwrite,
                    TrimAl = run.TrimAl,
                    TrimAl.path = trimAl.path,
                    trim.column = trim.column,
                    convert.ambiguous.sites = convert.ambiguous.sites,
                    alignment.assess = alignment.assess,
                    trim.external = trim.external,
                    trim.coverage = trim.coverage,
                    min.coverage.percent = min.coverage.percent,
                    min.external.percent = min.external.percent,
                    min.column.gap.percent = min.column.gap.percent,
                    min.alignment.length = min.alignment.length,
                    min.taxa.alignment = min.taxa.alignment,
                    max.alignment.gap.percent = max.alignment.gap.percent,
                    min.coverage.bp = min.coverage.bp,
                    threads = threads,
                    memory = memory)

#Aligns the mitogenomes and outputs summary stats
alignMitogenomes(alignment.folder = "Alignments/untrimmed-alignments",
                 reference.name = "reference",
                 draft.contigs = "draftContigs",
                 output.dir = "MitoGenomes",
                 dataset.name = "untrimmed",
                 overwrite = overwrite)

#Aligns the mitogenomes and outputs summary stats
alignMitogenomes(alignment.folder = "Alignments/trimmed-alignments",
                 reference.name = "reference",
                 draft.contigs = "draftContigs",
                 output.dir = "MitoGenomes",
                 dataset.name = "trimmed",
                 overwrite = FALSE)


### Create final mitogenomes
# buildMitogenomes(annotation.dir = "Annotations",
#                  alignment.folder = "Alignments/untrimmed-alignments",
#                  genome.alignment = "Genomes/alignments/untrimmed_mitogenome_alignment.phy",
#                  genome.dir = "Genomes",
#                  output.dir = NULL,
#                  overwrite = FALSE)
#

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

