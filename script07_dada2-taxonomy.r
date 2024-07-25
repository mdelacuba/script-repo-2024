# ASSIGNING TAXONOMY TO ASVs #
# -------------------------- #

### Assigning taxonomy until genus level ----

# Load data needed for these steps:
seqtab.nochim <- readRDS(file = "./rds-files/seqtab_nochim.rds")
seqtab.nochim_trans <- readRDS(file = "./rds-files/seqtab_nochim_trans.rds")
database_dir <- readRDS(file = "./rds-files/database_dir.rds")

# Initialize random number generator for reproducibility:
set.seed(100)

# Indicate databases:
silva <- paste0(database_dir,"silva_nr99_v138.1_train_set.fa.gz")
GTDB <- paste0(database_dir,"GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz")

# Assign with Silva:
taxboot.silva <- assignTaxonomy(seqtab.nochim, refFasta = silva, minBoot = 80, outputBootstraps = TRUE, verbose = FALSE, multithread = 70)
tax.silva <- assignTaxonomy(seqtab.nochim, refFasta = silva, minBoot = 80, verbose = FALSE, multithread = 70) # sin outputBootstraps porque interfiere con la adición de la columna de especies (OJO, especies no necesita bootstrap, ya que la asignación se hace sólo a 100% identidad)

# Assign with GTDB:
taxboot.gtdb <- assignTaxonomy(seqtab.nochim, refFasta = GTDB, minBoot = 80, outputBootstraps = TRUE, verbose = TRUE, multithread = 70)
tax.gtdb <- assignTaxonomy(seqtab.nochim, refFasta = GTDB, minBoot = 80, verbose = TRUE, multithread = 70)


### Assigning taxonomy to species level ----

# Indicate databases:
silva.s <- paste0(database_dir,"silva_species_assignment_v138.1.fa.gz")
GTDB.s <- paste0(database_dir,"GTDB_bac120_arc122_ssu_r202_Species.fa.gz")

# Assign with Silva:
#species.silva <- assignSpecies(seqtab.nochim, silva_species, n = 100)
#species.all.silva <- assignSpecies(seqtab.nochim, silva.s, allowMultiple=TRUE, n = 100)

# Assign with GTDB:
#species.gtdb <- assignSpecies(seqtab.nochim, GTDB_species, n = 100)
#species.all.gtdb <- assignSpecies(seqtab.nochim, GTDB.s, allowMultiple=TRUE, n = 100)

# Add species to the previous taxonomy until genus:
tax.all.silva <- addSpecies(tax.silva, silva.s, allowMultiple = TRUE, n = 100)
tax.all.gtdb <- addSpecies(tax.gtdb, GTDB.s, allowMultiple = TRUE, n = 100)

# Write tables of taxonomy:
write.table(tax.silva, "./results/taxonomy-silva-genus.tsv", sep="\t",  na = "NA", row.names = T, col.names = NA)
write.table(tax.all.silva, "./results/taxonomy-silva-all.tsv", sep="\t",  na = "NA", row.names = T, col.names = NA)
write.table(tax.gtdb, "./results/taxonomy-gtdb-genus.tsv", sep="\t",  na = "NA", row.names = T, col.names = NA)
write.table(tax.all.gtdb, "./results/taxonomy-gtdb-all.tsv", sep="\t",  na = "NA", row.names = T, col.names = NA)


### Appending taxonomy, boot, and seqtable ----

#--- Create dataframes and append: 
t_tax.all.silva <- as.data.frame(tax.all.silva)
taxa_boot.silva <- as.data.frame(taxboot.silva$boot) %>% rename_all(funs(str_c(., "_boot")))
seqtab.full.silva <- t_tax.all.silva %>% bind_cols(taxa_boot.silva) %>% bind_cols(seqtab.nochim_trans)

t_tax.all.gtdb <- as.data.frame(tax.all.gtdb)
taxa_boot.gtdb <- as.data.frame(taxboot.gtdb$boot) %>% rename_all(funs(str_c(., "_boot")))
seqtab.full.gtdb <- t_tax.all.gtdb %>% bind_cols(taxa_boot.gtdb) %>% bind_cols(seqtab.nochim_trans)

#--- Write tables of seqtable and taxonomy:
write.table(seqtab.full.silva, "./results/seqtable-silva.tsv", sep="\t",  na = "NA", row.names = T, col.names = NA)
write.table(seqtab.full.gtdb, "./results/seqtable-gtdb.tsv", sep="\t",  na = "NA", row.names = T, col.names = NA)


### Removing needless/undesired assignments ----

#--- Remove chloroplasts, mitochondria, and NA after Kingdom assignments (only Silva):
is.clmina <- tax.all.silva[,"Order"] %in% "Chloroplast" | 
  tax.all.silva[,"Family"] %in% "Mitochondria" |
  tax.all.silva[,"Kingdom"] %in% "Eukaryota" |
  tax.all.silva[,"Kingdom"] %in% NA
# | tax.all.silva[,"Phylum"] %in% NA # Kingdom is included in the removal of NA phylum by default

seqtab.nochim.ok <- seqtab.nochim[,!is.clmina]
t_seqtab.nochim.ok <- t(seqtab.nochim.ok)
tax.all.silva.ok <- tax.all.silva[!is.clmina,]
seqtab.full.silva.ok <- seqtab.full.silva[!is.clmina,]

#--- Save tables without undesired assginemts:
write.table(t_seqtab.nochim.ok, "./results/seqtable-nochim-ok.tsv", sep="\t",  na = "NA", row.names = T, col.names = NA)
write.table(tax.all.silva.ok, "./results/taxonomy-silva-ok.tsv", sep="\t",  na = "NA", row.names = T, col.names = NA)
write.table(seqtab.full.silva.ok, "./results/seqtable-silva-ok.tsv", sep="\t",  na = "NA", row.names = T, col.names = NA)

#--- Inspect distribution of sequence lengths
distrib.seqslength.ok <- table(nchar(getSequences(seqtab.nochim.ok)))
as.data.frame(t(distrib.seqslength.ok))


#--- Save things -----------------------------------------------

# Save objects needed for upstream steps:
saveRDS(tax.all.silva, file = "./rds-files/tax_all_silva.rds")
saveRDS(tax.all.gtdb, file = "./rds-files/tax_all_gtdb.rds")
saveRDS(seqtab.full.silva, file = "./rds-files/seqtab_full_silva.rds")
saveRDS(seqtab.full.gtdb, file = "./rds-files/seqtab_full_gtdb.rds")
saveRDS(seqtab.nochim.ok, file = "./rds-files/seqtab_nochim_ok.rds")
saveRDS(tax.all.silva.ok, file = "./rds-files/tax_all_silva_ok.rds")
saveRDS(seqtab.full.silva.ok, file = "./rds-files/seqtab_full_silva_ok.rds")

# Save the current work and objects:
save.image("R-taxonomy.RData")

