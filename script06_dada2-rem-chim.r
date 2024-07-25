# REMOVING CHIMERAS AND TRACKING SEQUENCES #
# ---------------------------------------- #

### Joining data of all sequencing runs ----

#--- Load data needed for this step:
results_dir <- readRDS(file = "./rds-files/results_dir.rds")

seqtab.1821 <- readRDS(file = "./rds-files/seqtab_1821.rds")
seqtab.2013 <- readRDS(file = "./rds-files/seqtab_2013.rds")
seqtab.2013cgb <- readRDS(file = "./rds-files/seqtab_2013cgb.rds")
seqtab.2014 <- readRDS(file = "./rds-files/seqtab_2014.rds")
seqtab.2015 <- readRDS(file = "./rds-files/seqtab_2015.rds")
seqtab.2020ss <- readRDS(file = "./rds-files/seqtab_2020ss.rds")
seqtab.h2000l7 <- readRDS(file = "./rds-files/seqtab_h2000l7.rds")
seqtab.h2000l8 <- readRDS(file = "./rds-files/seqtab_h2000l8.rds")
seqtab.h2500l3 <- readRDS(file = "./rds-files/seqtab_h2500l3.rds")
seqtab.h2500l4 <- readRDS(file = "./rds-files/seqtab_h2500l4.rds")
seqtab.h2500l5 <- readRDS(file = "./rds-files/seqtab_h2500l5.rds")
seqtab.h2500l7 <- readRDS(file = "./rds-files/seqtab_h2500l7.rds")

out.1821 <- readRDS(file = "./rds-files/out_1821.rds")
out.2013 <- readRDS(file = "./rds-files/out_2013.rds")
out.2013cgb <- readRDS(file = "./rds-files/out_2013cgb.rds")
out.2014 <- readRDS(file = "./rds-files/out_2014.rds")
out.2015 <- readRDS(file = "./rds-files/out_2015.rds")
out.2020ss <- readRDS(file = "./rds-files/out_2020ss.rds")
out.h2000l7 <- readRDS(file = "./rds-files/out_h2000l7.rds")
out.h2000l8 <- readRDS(file = "./rds-files/out_h2000l8.rds")
out.h2500l3 <- readRDS(file = "./rds-files/out_h2500l3.rds")
out.h2500l4 <- readRDS(file = "./rds-files/out_h2500l4.rds")
out.h2500l5 <- readRDS(file = "./rds-files/out_h2500l5.rds")
out.h2500l7 <- readRDS(file = "./rds-files/out_h2500l7.rds")

dada.1821 <- readRDS(file = "./rds-files/dada_1821.rds")
dada.2013 <- readRDS(file = "./rds-files/dada_2013.rds")
dada.2013cgb <- readRDS(file = "./rds-files/dada_2013cgb.rds")
dada.2014 <- readRDS(file = "./rds-files/dada_2014.rds")
dada.2015 <- readRDS(file = "./rds-files/dada_2015.rds")
dada.2020ss <- readRDS(file = "./rds-files/dada_2020ss.rds")
dada.h2000l7 <- readRDS(file = "./rds-files/dada_h2000l7.rds")
dada.h2000l8 <- readRDS(file = "./rds-files/dada_h2000l8.rds")
dada.h2500l3 <- readRDS(file = "./rds-files/dada_h2500l3.rds")
dada.h2500l4 <- readRDS(file = "./rds-files/dada_h2500l4.rds")
dada.h2500l5 <- readRDS(file = "./rds-files/dada_h2500l5.rds")
dada.h2500l7 <- readRDS(file = "./rds-files/dada_h2500l7.rds")

sample.names.1821 <- readRDS(file = "./rds-files/sample_names_1821.rds")
sample.names.2013 <- readRDS(file = "./rds-files/sample_names_2013.rds")
sample.names.2013cgb <- readRDS(file = "./rds-files/sample_names_2013cgb.rds")
sample.names.2014 <- readRDS(file = "./rds-files/sample_names_2014.rds")
sample.names.2015 <- readRDS(file = "./rds-files/sample_names_2015.rds")
sample.names.2020ss <- readRDS(file = "./rds-files/sample_names_2020ss.rds")
sample.names.h2000l7 <- readRDS(file = "./rds-files/sample_names_h2000l7.rds")
sample.names.h2000l8 <- readRDS(file = "./rds-files/sample_names_h2000l8.rds")
sample.names.h2500l3 <- readRDS(file = "./rds-files/sample_names_h2500l3.rds")
sample.names.h2500l4 <- readRDS(file = "./rds-files/sample_names_h2500l4.rds")
sample.names.h2500l5 <- readRDS(file = "./rds-files/sample_names_h2500l5.rds")
sample.names.h2500l7 <- readRDS(file = "./rds-files/sample_names_h2500l7.rds")

#--- Convert matrix to dataframe:
out.1821 <- as.data.frame(out.1821)
out.2013 <- as.data.frame(out.2013)
out.2013cgb <- as.data.frame(out.2013cgb)
out.2014 <- as.data.frame(out.2014)
out.2015 <- as.data.frame(out.2015)
out.2020ss <- as.data.frame(out.2020ss)
out.h2000l7 <- as.data.frame(out.h2000l7)
out.h2000l8 <- as.data.frame(out.h2000l8)
out.h2500l3 <- as.data.frame(out.h2500l3)
out.h2500l4 <- as.data.frame(out.h2500l4)
out.h2500l5 <- as.data.frame(out.h2500l5)
out.h2500l7 <- as.data.frame(out.h2500l7)

#--- Join data:
seqtab <- mergeSequenceTables(seqtab.2013, seqtab.2013cgb, seqtab.2014, seqtab.2015, seqtab.1821, seqtab.2020ss, seqtab.h2000l7, seqtab.h2000l8, seqtab.h2500l3, seqtab.h2500l4, seqtab.h2500l5, seqtab.h2500l6, seqtab.h2500l7)
out <- bind_rows(out.1821, out.2013, out.2013cgb, out.2014, out.2015, out.2020ss, out.h2000l7, out.h2000l8, out.h2500l3, out.h2500l4, out.h2500l5, out.h2500l7) 
dada <- c(dada.1821, dada.2013, dada.2013cgb, dada.2014, dada.2015, dada.2020ss, dada.h2000l7, dada.h2000l8, dada.h2500l3, dada.h2500l4, dada.h2500l5, dada.h2500l7)
sample.names <- c(sample.names.1821, sample.names.2013, sample.names.2013cgb, sample.names.2014, sample.names.2015, sample.names.2020ss, sample.names.h2000l7, sample.names.h2000l8, sample.names.h2500l3, sample.names.h2500l4, sample.names.h2500l5, sample.names.h2500l7)


### Refining the data ----

# Cluster sample IDs with less than MINREADS from sequence table (which is just an R matrix):
exc.samples <- rownames(seqtab[rowSums(seqtab) <= 9000,])
# Note: Filtering at <=9000 ensure including only samples with >10,000 input/raw reads

# Remove samples recovering less than 50% of sequences after chimera removal (specific IDs):
seqtab <- seqtab[!rownames(seqtab) %in% c(exc.samples,"ERR1769216", "ERR1769583", "ERR1769587", "ERR1769626", "ERR1769632"),]
# Note: This step was incorporated after tracking the number of sequences when running the script without this step.

# Filter by ASV length:
#seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 231:238]

# Inspect distribution of sequence lengths::
distrib.seqslength.prev <- table(nchar(getSequences(seqtab)))

# Collapse ASVs that only differs in length and without mismatches:
seqtab <- collapseNoMismatch(seqtab, minOverlap = 100, verbose = TRUE) # no se usó para prueba-1

# Inspect distribution of sequence lengths::
distrib.seqslength.tmp <- table(nchar(getSequences(seqtab)))

#--- Writing the ASV table in a file:
t_seqtab <- t(seqtab)
write.table(t_seqtab,file = "./results/seqtable-tmp.csv", sep="\t",  na = "NA", row.names = T, col.names = NA)


### Removing chimeras ----

#--- Remove chimeras:
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = 6, verbose = TRUE)

#--- Compute % of non-chimeras sequences:
p_seqtab.nochim <- paste0("% of non chimeras : ", sum(seqtab.nochim)/sum(seqtab) * 100)
n_seqtab.nochim <- paste0("total number of sequences : ", sum(seqtab.nochim))

# Writing the non-chimeric ASV table in a file:
t_nochim <- t(seqtab.nochim)
write.table(t_nochim,file ="./results/seqtable-nochim.csv", sep="\t",  na = "NA", row.names = T, col.names = NA)

#--- Inspect distribution of sequence lengths:
distrib.seqslength <- table(nchar(getSequences(seqtab.nochim)))


### Tracking the number of sequences at each step ----

# Define a function:
getN <- function(x) sum(getUniques(x))

#--- Modify rownames in "out" table:
rownames(out) <- str_split(basename(rownames(out)), pattern = "_", simplify = TRUE)[, 1]

# Remove the samples recovering less than 50% of sequences after chimera removal ():
out <- out[!rownames(out) %in% c(exc.samples, "ERR1769216", "ERR1769583", "ERR1769587", "ERR1769626", "ERR1769632"),]
dada <- dada[!dada %in% dada[c(exc.samples, "ERR1769216", "ERR1769583", "ERR1769587", "ERR1769626", "ERR1769632")]]
sample.names <- sample.names[!sample.names %in% c(exc.samples, "ERR1769216", "ERR1769583", "ERR1769587", "ERR1769626", "ERR1769632")]
# Note: This step was incorporated after tracking the number of sequences when running the script without this step.

# Build table of stats:
track <- cbind(out, sapply(dada, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered","denoised", "tabled", "nonchim")
rownames(track) <- sample.names

#--- Save table stats:
write.table(track, file = "./results/summarize-stats-previous.csv", sep="\t",  na = "NA", row.names = T, col.names = NA)


### Transforming and saving the ASVs sequences in fasta format ----

#--- Transform ASV sequences to ASV IDs:
seqtab.nochim_trans <- as.data.frame(t(seqtab.nochim)) %>% rownames_to_column(var = "sequence") %>% 
  rowid_to_column(var = "ASVNumber") %>% mutate(ASVNumber = sprintf("ASV%04d", 
                                                                    ASVNumber)) %>% mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))
df <- seqtab.nochim_trans
seq_out <- Biostrings::DNAStringSet(df$sequence)

names(seq_out) <- df$ASVNumber

#--- Save ASV sequences in fasta format:
Biostrings::writeXStringSet(seq_out, str_c(results_dir, "ASV_no_taxo.fasta"), 
                            compress = FALSE, width = 20000)

#--- Save new ASV table:
write.table(df, file = './results/ASVtab.tsv', sep='\t', row.names = FALSE, na='', quote=FALSE)


#--- Save things -----------------------------------------------

#--- Save objects needed for upstream steps:
saveRDS(seqtab.nochim, file = "./rds-files/seqtab_nochim.rds")
saveRDS(seqtab.nochim_trans, file = "./rds-files/seqtab_nochim_trans.rds")
saveRDS(out, file ="./rds-files/out.rds")
saveRDS(dada, file ="./rds-files/dada.rds")
saveRDS(sample.names, file ="./rds-files/sample_names.rds")

#--- Save the current work and objects:
save.image("R-rem-chimeras.RData")

