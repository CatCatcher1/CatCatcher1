library(dada2)
library(phyloseq)
library(ggplot2)
path <-"C:/Users/Batuhan/Desktop/Medipol Staj/new_data"
list.files(path)

fnFs <- sort(list.files(path, pattern="_trimmed_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_trimmed_2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


plotQualityProfile(fnFs[1:10])
plotQualityProfile(fnRs[1:10])

filtFs <- file.path(path, "filtered", paste0(sample.names, "_1_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_2_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

plotQualityProfile(filtFs[1:10])
plotQualityProfile(filtRs[1:10])

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
dadaRs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/Batuhan/Desktop/Medipol Staj/22032024_data/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)

theme_set(theme_bw())

samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

metadata <- read.delim("C:/Users/Batuhan/Desktop/Medipol Staj/22032024_data/metadata.txt", header=TRUE, sep="\t")
rownames(metadata) <- metadata$Sample

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(taxa), sample_data(metadata))

ps <- prune_samples(sample_names(ps) != "Mock", ps)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

plot_richness(ps, x="BodySite", measures=c("Shannon", "Simpson"), color="BodySite")

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")

###
top20 <- names(sort(taxa_sums(ps.prop.genus), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps.prop.genus, function(otu) otu/sum(otu))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Sample", fill="Family") + facet_wrap(~BodySite, scales="free_x")

plot_bar(ps.prop.genus, x="Sample", fill="Family") + facet_wrap(~BodySite, scales="free_x"+ theme(legend_position = "bottom") )
library(phyloseq)

ps.prop.genus<-tax_agg(ps.prop, rank="Genus")

###
install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)

library(microViz)

ps.prop.genus <- ps.prop %>%
  tax_fix() %>%
  tax_agg(rank="Genus")

top20 <- names(sort(taxa_sums(ps.prop.genus), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps.prop.genus, function(otu) otu/sum(otu))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Sample", fill="Family") + facet_wrap(~BodySite, scales="free_x")