# PRE-PROCESSING THE DATA: ABUNDANCE MATRIX, TAXONOMY TABLE, AND METADATA #
# ----------------------------------------------------------------------- #

#--- Load objects needed for these steps ----

seqtab.full.silva.ok <- readRDS(file = "./rds-files/seqtab_full_silva_ok.rds")
seqtab.nochim.ok <- readRDS(file = "./rds-files/seqtab_nochim_ok.rds")
tax.all.silva.ok <- readRDS(file = "./rds-files/tax_all_silva_ok.rds")


### Phyloseq importing ----

#--- ASV and taxa importing:
asv.tmp <- otu_table(seqtab.nochim.ok, taxa_are_rows = FALSE)
tax <- tax_table(tax.all.silva.ok)

#--- Isolate rownames to make up the metadata file:
ids <- rownames(seqtab.nochim.ok)
write.table(ids, file = "./ids-seqtab-nochim-ok.tsv", row.names = F, col.names = F)

#--- Load sample metadata files:
metadata <- read.table(file = "metadata-179sponges", header = TRUE, sep = "\t")

#--- Reorder sample metadata by ID:
metadata <- metadata[order(metadata$ID),]

#--- Add rownames to metadata considering the new order:
rownames(metadata) <- rownames(asv)

#--- Metadata importing:
sdata <- sample_data(metadata)

#--- Make first phyloseq object:
physeq.tmp <- phyloseq(asv, tax, sdata)
sample_names(physeq.tmp)

#--- Build phylogenetic tree:
random.tree <- rtree(ntaxa(physeq.tmp), rooted = TRUE, tip.label = taxa_names(physeq.tmp))

#--- Merge phyloseq object again with the tree:
physeq <- merge_phyloseq(physeq.tmp, random.tree)


### Filtering singletons and doubletons ----

#--- Observe #ASVs with zeros and/or singletons:
min(taxa_sums(physeq)) # 2, no singletons
#table(taxa_sums(physeq) <= 1) # with sum
#sum(taxa_sums(physeq) == 0)
#table(round(colSums(otu_table(physeq))/colSums(!!otu_table(physeq))) <= 1) # with mean

#--- Filter singletons and doubletons:
physeq <- filter_taxa(physeq, function(x) sum(x) > 2, TRUE)
#physeq <- filter_taxa(physeq, function(x) max(x) > 10, TRUE)
#prune_taxa(taxa_sums(physeq) > 5, physeq)


### Turning ASVs to short names ----

#--- Replace sequences with identifiers:
dna <- DNAStringSet(taxa_names(physeq))
names(dna) <- taxa_names(physeq)
physeq <- merge_phyloseq(physeq, dna)
taxa_names(physeq) <- paste0("ASV", seq(ntaxa(physeq)))
refseqs <- refseq(physeq)

#--- Save refseqs in table:
write.table(refseqs, file = "./results/refseqs.tsv", 
            sep = "\t", na = "NA", row.names = T, col.names = FALSE)

#--- Melt physeq table:
# physeq.melt <- psmelt(physeq)

#--- Make table:
# write.table(physeq.melt,file ="./results/physeq_melt.csv", 
#             sep = "\t",  na = "NA", row.names = T, col.names = NA)

#--- Inspect distribution of sequence lengths:
distrib.seqslength.physeq <- table(nchar(getSequences(refseqs)))
as.data.frame(t(distrib.seqslength.physeq))

#--- Extract taxonomy table and abundance matrix:
seqtab <- as.data.frame(otu_table(physeq))
tax <- data.frame(tax_table(physeq))

#--- Write tables:
write.table(seqtab, file = "./results/seqtable-sponge.tsv", 
            sep = "\t", na = "NA", row.names = T, col.names = NA)
write.table(tax, file = "./results/taxonomy-silva-sponge.tsv", 
            sep = "\t", na = "NA", row.names = T, col.names = NA)


### Transforming absolute to relative abundance ----

#--- Shorten name to physeq object:
pqs <- physeq

#--- Transform:
pqs.rel <- transform_sample_counts(pqs, function(x) x/sum(x))

#--- Filter sponge taxa by max relative abundance (higher than 0.1%):
#pqs.high.rel <- filter_taxa(pqs.rel, function(x) max(x) > 0.001, TRUE)


### Inspecting the metadata ----

#---  Inspect the metadata:
table(sample_data(pqs)$Environment)
table(sample_data(pqs)$Habitat)
table(sample_data(pqs)$Sponge_Genus)
length(table(sample_data(pqs)$Sponge_Species))
table(sample_data(pqs)$Sequencing_run_ID)
length(table(sample_data(pqs)$Host_status))
table(sample_data(pqs)$Sponge_Family)
table(sample_data(pqs)$Sponge_Order)
table(sample_data(pqs)$Marine_Ecoregion)
table(sample_data(pqs)$Marine_Province)
table(sample_data(pqs)$Marine_Realm)
sample_data(pqs)$Microbiome_ID


### Inspecting the abundance matrix ----

#--- Inspect total counts (sum):
sum.counts <- data.frame(sumCount= taxa_sums(pqs))
sum.counts <- rownames_to_column(sum.counts, "ASVnum")
sum.counts <- arrange(sum.counts, -sumCount)

#--- Inspect mean counts (only non-zero):
mean.counts <- data.frame(meanCount = round(colSums(otu_table(pqs))/colSums(!!otu_table(pqs))))
mean.counts <- rownames_to_column(mean.counts, "ASVnum")
mean.counts <- arrange(mean.counts, -meanCount)

# Inspect max counts:
max.counts <- data.frame(maxCount = apply(otu_table(pqs), 2, max))
max.counts <- rownames_to_column(max.counts, "ASVnum")
max.counts <- arrange(max.counts, -maxCount)

#--- Bind sum and mean in a table:
counts.tab <- cbind(sum.counts, mean.counts, max.counts)

#--- Calculate medians:
median(sum.counts$sumCount)
median(mean.counts$meanCount)

#--- Frequency of ASV counts:
sum.freq <- data.frame(table(sumCount = sort(counts.tab$sumCount)))
mean.freq <- data.frame(table(meanCount = sort(counts.tab$meanCount)))
max.freq <- data.frame(table(maxCount = sort(counts.tab$maxCount)))

#--- Evaluate prevalence of taxa across samples:
preval <- data.frame(table(nSamples = colSums(!!otu_table(pqs))))
#sum(preval.o$Freq)


### Inspecting the taxonomy table ----

#--- Show number of ASVs of Bacteria and Archaea:
table(tax_table(pqs)[,1])

#--- Show the lineages found by taxonomic level:
kingdom.taxa <- get_taxa_unique(pqs, "Kingdom")
phylum.taxa <- get_taxa_unique(pqs, "Phylum")
class.taxa <- get_taxa_unique(pqs, "Class")
order.taxa <- get_taxa_unique(pqs, "Order")
family.taxa <- get_taxa_unique(pqs, "Family")
genus.taxa <- get_taxa_unique(pqs, "Genus")
species.taxa <- get_taxa_unique(pqs, "Species")

#--- Count the lineages found by taxonomic level:
length(kingdom.taxa)
length(phylum.taxa)
length(class.taxa)
length(order.taxa)
length(family.taxa)
length(genus.taxa)
length(species.taxa)
# Note: Consider that each resulting number contains an "NA", so rest 1 to each.

#--- Show number of NAs:
sum(is.na(tax.s$Kingdom))
sum(is.na(tax.s$Phylum))
sum(is.na(tax.s$Class))
sum(is.na(tax.s$Order))
sum(is.na(tax.s$Family))
sum(is.na(tax.s$Genus))
sum(is.na(tax.s$Species))

#--- Calculate the % assigned ASVs:

# For Antarctic sponges:
p.ant <- 100-(round(sum(is.na(data.frame(tax_table(
  prune_taxa(taxa_sums(subset_samples(
    pqs, Habitat == "Antarctic sponge")) > 0, subset_samples(
      pqs, Habitat == "Antarctic sponge"))))$Phylum))*100/11675))
c.ant <- 100-(round(sum(is.na(data.frame(tax_table(
  prune_taxa(taxa_sums(subset_samples(
    pqs, Habitat == "Antarctic sponge")) > 0, subset_samples(
      pqs, Habitat == "Antarctic sponge"))))$Class))*100/11675))
o.ant <- 100-(round(sum(is.na(data.frame(tax_table(
  prune_taxa(taxa_sums(subset_samples(
    pqs, Habitat == "Antarctic sponge")) > 0, subset_samples(
      pqs, Habitat == "Antarctic sponge"))))$Order))*100/11675))
f.ant <- 100-(round(sum(is.na(data.frame(tax_table(
  prune_taxa(taxa_sums(subset_samples(
    pqs, Habitat == "Antarctic sponge")) > 0, subset_samples(
      pqs, Habitat == "Antarctic sponge"))))$Family))*100/11675))
g.ant <- 100-(round(sum(is.na(data.frame(tax_table(
  prune_taxa(taxa_sums(subset_samples(
    pqs, Habitat == "Antarctic sponge")) > 0, subset_samples(
      pqs, Habitat == "Antarctic sponge"))))$Genus))*100/11675))
s.ant <- 100-(round(sum(is.na(data.frame(tax_table(
  prune_taxa(taxa_sums(subset_samples(
    pqs, Habitat == "Antarctic sponge")) > 0, subset_samples(
      pqs, Habitat == "Antarctic sponge"))))$Species))*100/11675))

# For non-Antarctic sponges:
p.noa <- 100-(round(sum(is.na(data.frame(tax_table(
  prune_taxa(taxa_sums(subset_samples(
    pqs, Habitat == "Non-Antarctic sponge")) > 0, subset_samples(
      pqs, Habitat == "Non-Antarctic sponge"))))$Phylum))*100/18119))
c.noa <- 100-(round(sum(is.na(data.frame(tax_table(
  prune_taxa(taxa_sums(subset_samples(
    pqs, Habitat == "Non-Antarctic sponge")) > 0, subset_samples(
      pqs, Habitat == "Non-Antarctic sponge"))))$Class))*100/18119))
o.noa <- 100-(round(sum(is.na(data.frame(tax_table(
  prune_taxa(taxa_sums(subset_samples(
    pqs, Habitat == "Non-Antarctic sponge")) > 0, subset_samples(
      pqs, Habitat == "Non-Antarctic sponge"))))$Order))*100/18119))
f.noa <- 100-(round(sum(is.na(data.frame(tax_table(
  prune_taxa(taxa_sums(subset_samples(
    pqs, Habitat == "Non-Antarctic sponge")) > 0, subset_samples(
      pqs, Habitat == "Non-Antarctic sponge"))))$Family))*100/18119))
g.noa <- 100-(round(sum(is.na(data.frame(tax_table(
  prune_taxa(taxa_sums(subset_samples(
    pqs, Habitat == "Non-Antarctic sponge")) > 0, subset_samples(
      pqs, Habitat == "Non-Antarctic sponge"))))$Genus))*100/18119))
s.noa <- 100-(round(sum(is.na(data.frame(tax_table(
  prune_taxa(taxa_sums(subset_samples(
    pqs, Habitat == "Non-Antarctic sponge")) > 0, subset_samples(
      pqs, Habitat == "Non-Antarctic sponge"))))$Species))*100/18119))

#--- Join the abundance matrix and the taxonomy table:
tax.pqs <- as.data.frame(tax_table(pqs))
seqtab.pqs <- as.data.frame(t(otu_table(pqs)))
metatab <- cbind(tax.pqs, seqtab.pqs)

#--- Write table:
write.table(metatab, file = "./results/metatab_pqs.tsv", 
            sep = "\t", na = "NA", row.names = T, col.names = NA)


#--- Save things------------------------------------------------------------

# Save objects for downstream plotting:
saveRDS(pqs, file = "./rds-files/pqs.rds")
saveRDS(pqs.rel, file = "./rds-files/pqs_rel.rds")
saveRDS(metadata, file = "./rds-files/metadata.rds")
# saveRDS(pqs.high.rel, file = "./rds-files/pqs_high_rel.rds")

# Save the current work and objects:
save.image("R-phyloseq.RData")

