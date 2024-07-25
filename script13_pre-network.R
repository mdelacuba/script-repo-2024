# PRE-PROCESSING: INFERENCE OF MICROBIAL CO-OCCURRENCE NETWORK #
# ------------------------------------------------------------ #

#--- Load required files ----

pqs <- readRDS(file = "./rds-files/pqs.rds")
pqs.rel <- readRDS(file = "./rds-files/pqs_rel.rds")
metadata <- readRDS(file = "./rds-files/metadata.rds")
pqs.raref <- readRDS(file = "./rds-files/pqs_raref.rds")
gen.asvs <- readRDS(file = "./rds-files/gen-asvs.rds")
spec.asvs <- readRDS(file = "./rds-files/spec-asvs.rds")
spec.noan.asvs <- readRDS(file = "./rds-files/spec-noan-asvs.rds")
# pqs.noant <- readRDS(file = "./rds-files/pqs_noant.rds")


### Isolating subsets at ASV level ----

# Antarctic subset: 
pqs.ant <- subset_samples(pqs.raref, Habitat == "Antarctic sponge")
pqs.ant <- prune_taxa(taxa_sums(pqs.ant) > 0, pqs.ant)

# Non-Antarctic subset: 
pqs.noant <- subset_samples(pqs.raref, Habitat == "Non-Antarctic sponge")
pqs.noant <- prune_taxa(taxa_sums(pqs.noant) > 0, pqs.noant)

#--- Filter taxa by prevalence:
# fte.ant <- genefilter_sample(pqs.ant, filterfun_sample(function(x) x > 0), A = 0.33*nsamples(pqs.ant))
# pqs.ant <- prune_taxa(fte.ant, pqs.ant)
# sample_sums(pqs.ant)
# 
# fte.noant <- genefilter_sample(pqs.noant, filterfun_sample(function(x) x > 0), A = 0.33*nsamples(pqs.noant))
# pqs.noant <- prune_taxa(fte.noant, pqs.noant)
# sample_sums(pqs.noant)

#--- Transpose and isolate tables:
seqtab.ant <- data.frame(t(otu_table(pqs.ant)))
taxtab.ant <- as.data.frame(tax_table(pqs.ant))

seqtab.noant <- data.frame(t(otu_table(pqs.noant)))
taxtab.noant <- as.data.frame(tax_table(pqs.noant))


### Changing NAs to higher-level taxonomy names ----

## Non-Antarctic subset:

#--- Replace from taxonomy table:
tax.clean <- data.frame(row.names = row.names(taxtab.noant),
                             Kingdom = str_replace(taxtab.noant[,1], "D_0__",""),
                             Phylum = str_replace(taxtab.noant[,2], "D_1__",""),
                             Class = str_replace(taxtab.noant[,3], "D_2__",""),
                             Order = str_replace(taxtab.noant[,4], "D_3__",""),
                             Family = str_replace(taxtab.noant[,5], "D_4__",""),
                             Genus = str_replace(taxtab.noant[,6], "D_5__",""),
                        Species = str_replace(taxtab.noant[,7], "D_6__",""),
                             stringsAsFactors = FALSE)
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

tax_table(pqs.noant) <- as.matrix(tax.clean)
taxtab.noant <- as.data.frame(tax_table(pqs.noant))
# taxa_names(pqe) <- as.data.frame(tax_table(pqe))$Genus
# To replace in the phyloseq object instead of the seq and taxonomy table

#--- Look the number of habitat-specific and generalist ASVs:
length(intersect(rownames(taxtab.noant), spec.noan.asvs))
length(intersect(rownames(taxtab.noant), gen.asvs))

#--- Replace ASV ID with an alphabetical code (CoNet does not accept numeric codes):
code.noant <- rownames(taxtab.noant)
rownames(seqtab.noant) <- chartr("1234567890", "ABCDEFGHIJ", 
                                      paste0("ASV_", str_remove(rownames(seqtab.noant), "ASV")))
rownames(taxtab.noant) <- chartr("1234567890", "ABCDEFGHIJ", 
                                      paste0("ASV_", str_remove(rownames(taxtab.noant), "ASV")))
code.noant <- cbind(code.noant, rownames(taxtab.noant))
colnames(code.noant) <- c("ASV_num", "ASV_code")

#--- Save new seqtable and taxtable:
write.table(seqtab.noant, file = "./results/netw/seqtable-noa-netw.tsv", sep = "\t", na = "NA", row.names = T, col.names = NA)
write.table(taxtab.noant, file = "./results/netw/taxonomy-noa-netw.tsv", sep = "\t", na = "NA", row.names = T, col.names = F)

# Save ASV IDs to eliminate "" and replace manually to seqtable:
write.table(rownames(taxtab.noant), file = "./results/netw/rownames-noa-netw.tsv", sep = "\t", na = "NA", row.names = F, col.names = "col_names")
# NOTE: Remember to replace "-" with "_" manually, because CoNet does not support that character.
# NOTE2: Remember to remove "" from taxonomy table.
write.table(code.noant, file = "./results/netw/code-noa.tsv", sep = "\t", na = "NA", row.names = F, col.names = T)


## Antarctic subset:

#--- Replace from taxonomy table:
tax.clean2 <- data.frame(row.names = row.names(taxtab.ant),
                        Kingdom = str_replace(taxtab.ant[,1], "D_0__",""),
                        Phylum = str_replace(taxtab.ant[,2], "D_1__",""),
                        Class = str_replace(taxtab.ant[,3], "D_2__",""),
                        Order = str_replace(taxtab.ant[,4], "D_3__",""),
                        Family = str_replace(taxtab.ant[,5], "D_4__",""),
                        Genus = str_replace(taxtab.ant[,6], "D_5__",""),
                        Species = str_replace(taxtab.ant[,7], "D_6__",""),
                        stringsAsFactors = FALSE)
tax.clean2[is.na(tax.clean2)] <- ""

for (i in 1:7){ tax.clean2[,i] <- as.character(tax.clean2[,i])}
tax.clean2[is.na(tax.clean2)] <- ""

for (i in 1:nrow(tax.clean2)){
  if (tax.clean2[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean2[i,1], sep = "")
    tax.clean2[i, 2:7] <- kingdom
  } else if (tax.clean2[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean2[i,2], sep = "")
    tax.clean2[i, 3:7] <- phylum
  } else if (tax.clean2[i,4] == ""){
    class <- paste("Class_", tax.clean2[i,3], sep = "")
    tax.clean2[i, 4:7] <- class
  } else if (tax.clean2[i,5] == ""){
    order <- paste("Order_", tax.clean2[i,4], sep = "")
    tax.clean2[i, 5:7] <- order
  } else if (tax.clean2[i,6] == ""){
    family <- paste("Family_", tax.clean2[i,5], sep = "")
    tax.clean2[i, 6:7] <- family
  } else if (tax.clean2[i,7] == ""){
    tax.clean2$Species[i] <- paste("Genus",tax.clean2$Genus[i], sep = "_")
  }
}

tax_table(pqs.ant) <- as.matrix(tax.clean2)
taxtab.ant <- as.data.frame(tax_table(pqs.ant))
# taxa_names(pqe) <- as.data.frame(tax_table(pqe))$Genus
#Si quiero reemplazar los nombres desde el objeto phyloseq y no sólo en tax_table y el dataframe

#--- Look the number of habitat-specific and generalist ASVs (Reload rds files, because they were modified):
length(intersect(rownames(taxtab.ant), spec.asvs))
length(intersect(rownames(taxtab.ant), gen.asvs))

#--- Replace ASV ID by an alphabetical code (CoNet does not accept numeric codes):
code.ant <- rownames(taxtab.ant)
rownames(seqtab.ant) <- chartr("1234567890", "ABCDEFGHIJ", 
                                 paste0("ASV_", str_remove(rownames(seqtab.ant), "ASV")))
rownames(taxtab.ant) <- chartr("1234567890", "ABCDEFGHIJ", 
                                 paste0("ASV_", str_remove(rownames(taxtab.ant), "ASV")))
code.ant <- cbind(code.ant, rownames(taxtab.ant))
colnames(code.ant) <- c("ASV_num", "ASV_code")

spec.asvs <- chartr("1234567890", "ABCDEFGHIJ", 
                               paste0("ASV_", str_remove(spec.asvs, "ASV")))
gen.asvs <- chartr("1234567890", "ABCDEFGHIJ", 
                    paste0("ASV_", str_remove(gen.asvs, "ASV")))
spec.noan.asvs <- chartr("1234567890", "ABCDEFGHIJ", 
                    paste0("ASV_", str_remove(spec.noan.asvs, "ASV")))

#--- Save seqtable and taxtable:
write.table(seqtab.ant, file = "./results/netw/seqtable-ant-netw.tsv", sep = "\t", na = "NA", row.names = T, col.names = NA)
write.table(taxtab.ant, file = "./results/netw/taxonomy-ant-netw.tsv", sep = "\t", na = "NA", row.names = T, col.names = F)

#--- Save ASV IDs to eliminate "" and replace manually to seqtable:
write.table(rownames(taxtab.ant), file = "./results/netw/rownames-ant-netw.tsv", sep = "\t", na = "NA", row.names = F, col.names = "col_names")
# NOTE: Remember to replace "-" with "_" manually, because CoNet does not support that character.
# NOTE2: Remember to remove "" from taxonomy table.
write.table(code.ant, file = "./results/netw/code-ant.tsv", sep = "\t", na = "NA", row.names = F, col.names = T)

#--- Save habitat-specific and generalist ASVs in a list:
write.table(spec.asvs, file = "./results/netw/asvs-spec.tsv", sep = "\t", na = "NA", row.names = F, col.names = F)
write.table(spec.noan.asvs, file = "./results/netw/asvs-spec-noan.tsv", sep = "\t", na = "NA", row.names = F, col.names = F)
write.table(gen.asvs, file = "./results/netw/asvs-gen.tsv", sep = "\t", na = "NA", row.names = F, col.names = F)


### Creating a legend for the network ----

#--- Create dataframe with the phyla (from CoNet output):
leg.df <- data.frame(
  Xdata = rnorm(20),                        
  Ydata = rnorm(20),
  LegendData = c("Proteobacteria", "Bacteroidota", "Unassigned Bacteria",
                 "Firmicutes", "Cyanobacteria", "Verrucomicrobiota",
                 "Crenarchaeota", "Planctomycetota", "Actinobacteriota",
                 "Nitrospirota", "Acidobacteriota", "Bdellovibrionota",
                 "Campylobacterota", "Dadabacteria", "Fusobacteriota", 
                 "Marinimicrobia (SAR406 clade)", "Myxococcota", "Patescibacteria",
                 "SAR324 clade(Marine group B)", "Thermoplasmatota"))

#--- Plot legend with a temporal graphic:
leg.net <- ggplot(leg.df, aes(
  Xdata, Ydata, color = LegendData)) +   
  geom_point(size = 10) +
  labs(color = "Phylum") +
  scale_color_manual(values = c("#a6558b", "#eda137", "#878787", "#765a03", 
                                "#8CB369", "#146ba3", "#b8b893", "#117744",
                                "#792E09", "#FFEECF", "#6DAEDB", "#e97451",
                                "#031f76", "#528b8b", "#5a0376", "#37eda1",
                                "#bebebe", "#771144", "#9069b3", "#354230"),
                     breaks = c("Proteobacteria", "Bacteroidota", "Unassigned Bacteria",
                                "Firmicutes", "Cyanobacteria", "Verrucomicrobiota",
                                "Crenarchaeota", "Planctomycetota", "Actinobacteriota",
                                "Nitrospirota", "Acidobacteriota", "Bdellovibrionota",
                                "Campylobacterota", "Dadabacteria", "Fusobacteriota", 
                                "Marinimicrobia (SAR406 clade)", "Myxococcota", "Patescibacteria",
                                "SAR324 clade(Marine group B)", "Thermoplasmatota")) +
  theme_classic() +
  theme(legend.position = "right", text = element_text(size = 15))

#--- Plot just the legend:
legend.net <- cowplot::get_legend(leg.net)
legend.net <- as_ggplot(legend.net)

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/net_legend.png", units = "in", width = 3, height = 10, res = 300)
legend.net
dev.off()


#--- Save things -----------------------------------------------

# Save the current work and objects:
save.image("R-graph-netw.RData")
