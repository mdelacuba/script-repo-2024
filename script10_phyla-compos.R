# TAXONOMY COMPOSITION AT PHYLUM LEVEL #
# ------------------------------------ #

#--- Load required objects ----

pqs <- readRDS(file = "./rds-files/pqs.rds")
pqs.rel <- readRDS(file = "./rds-files/pqs_rel.rds")
pqs.ant.spec <- readRDS(file = "./rds-files/pqs_ant_spec.rds")
pqs.noa.spec <- readRDS(file = "./rds-files/pqs_noa_spec.rds")


#--- Define colors for graphics ----

colors.hab <- as.matrix(jcolors("rainbow"))[c(6, 9)]
colors.hab <- rev(brewer.pal(n = 11, name = "RdBu"))[c(3,9)]
colors.1 <- c("coral3", "darkorchid4", "deepskyblue2", "darkgoldenrod3", 
              "seagreen4", "tomato4", "darkblue", "turquoise", "lightpink", 
              "limegreen", "orchid1",  "olivedrab4", "red1", "skyblue3", "maroon4", 
              "darkorange", "gold2", "khaki", "plum4", "maroon2", "lightgreen", 
              "aquamarine", "grey0", "mediumpurple", "sienna4", "tan", "darkcyan", 
              "yellow", "thistle", "gray81")
phylum_colors <- read_excel("colors.xlsx", sheet = "phyla")
phylum_colors <-structure(phylum_colors$color_name, .Names = phylum_colors$phylum)
colors.ant <- c("#633353", "#d9b7ce", "#84446f", "#999999",
                "#d59031", "#eda137", "#f0b35e",  "#a6558b",
                "#b776a2",  "#f4c687", "#000000", "#e0c5d7", 
                "#c999b9", "#bad1a5")
colors.noa <- c("#a9b37c", "#a6558b", "#D4E09B", "#f4c687", "#d9b7ce")
colors.gen <- c( "#633353", "#84446f", "#a6558b", "#b776a2", "#d9b7ce", "#e0c5d7",
                 "#e8d3e1", "#efe2eb", "#d59031", "#eda137", "#f0b35e")


### Assessing the taxonomy composition at phylum level ----

#--- Separate only top 20 most dominant taxa across samples:
# top20 <- names(sort(taxa_sums(pqs), decreasing=TRUE))[1:50]
# pqs.top20 <- prune_taxa(top20, pqs.rel)
# 
# #--- Concatenate phyloseq table:
# pqs.top20.melt <- psmelt(pqs.top20)

#--- Transform data by Hellinger transformation:
pqs.hell <- transform_sample_counts(pqs.rel, function(x) sqrt(x))

#--- Transform to relative abundance:
pqs.hell.rel <- transform_sample_counts(pqs.hell, function(x) x/sum(x))

#--- Agglomerate to a specific taxonomic rank:
pqs.glom.e <- tax_glom(pqs.hell.rel, taxrank = "Class", NArm=FALSE, 
                          bad_empty=c(NA, "", " ", "\t"))

#--- Concatenate phyloseq table:
pqs.melt.e <- psmelt(pqs.glom.e)

# Create new column with the lineage:
pqs.melt.e$ID_Species <- paste(pqs.melt.e$Sponge_Species,
                               pqs.melt.e$Sample, sep = "; ")

#--- Write table:
write.table(pqs.melt.e,file ="./results/physeq_melt_e.csv",
            sep="\t",  na = "NA", row.names = T, col.names = NA)

# Customize a new column to appear Proteobacteria classes among phyla:
pqs.melt.e$custom <- pqs.melt.e$Phylum
pqs.melt.e$custom[pqs.melt.e$custom == "Proteobacteria"] <- "Alphaproteobacteria"
pqs.melt.e <- mutate(pqs.melt.e, custom = ifelse(Class == 'Gammaproteobacteria', Class, custom))

#--- Join low abundance taxa in a subcategory:
pqs.melt.e$custom <- as.character(pqs.melt.e$custom)
pqs.melt.e$custom[pqs.melt.e$Abundance < 0.05] <- "< 5% abundance"

#--- Change NA label:
pqs.melt.e$custom[is.na(pqs.melt.e$custom)] <- "Unassigned"

write.table(pqs.melt.e,file ="./results/physeq_merge_melt_e_join.csv",
            sep="\t",  na = "NA", row.names = T, col.names = NA)

#--- Count number of labels (lineages and "< 5% abundance") that would appear in graphics:
count.lineage.e <- length(unique(pqs.melt.e$Class))
as.matrix(unique(pqs.melt.e$custom))

#--- Reorder Phyla according to abundance:
pqs.melt.e$custom <- reorder(pqs.melt.e$custom, pqs.melt.e$Abundance)
pqs.melt.e$custom <- factor(pqs.melt.e$custom, levels = rev(levels(pqs.melt.e$custom)))

#--- Barplot of taxonomy composition at phylum level:
bar.plot.e <- ggplot(data = pqs.melt.e, aes(x = ID_Species, y = Abundance, fill = custom)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  theme_classic() +
  facet_wrap(~Habitat, nrow = 1, scale = "free_x") +
  theme(axis.text.x = element_blank(), 
    #axis.text.x = element_text(angle = 90, 
                                   #size = 4, fac = "bold.italic", hjust = 1, vjust = 0.5),
        text = element_text(size = 15), legend.text = element_text(size = 16),
        axis.ticks.x = element_blank(), #axis.ticks.x = element_line(linewidth = 0.3),
        legend.key.size = unit(0.7, "cm"),
        strip.text.x = element_text(size = 17)) +
  #theme(axis.text.x = element_text(angle = 90)) +
  #axis.text.x = element_blank(), axis.ticks.x = element_blank()
  #coord_flip() +
  scale_fill_manual(values = phylum_colors, 
                    limits = rownames(as.matrix(phylum_colors))) +
  ylab("Relative abundance") +
  labs(x = NULL, fill = NULL) +
  theme(legend.position = "right")

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/tax-envs.png", units = "in", width = 13, height = 7, res = 300)
bar.plot.e
dev.off()


### Comparing the relative abundance of unassigned taxa ----

# Search abundance of NAs in each habitat:
na.ant <- sort(pqs.melt.e$Abundance[pqs.melt.e$custom == "Unassigned" & pqs.melt.e$Habitat == "Antarctic sponge"])
na.noa <- sort(pqs.melt.e$Abundance[pqs.melt.e$custom == "Unassigned" & pqs.melt.e$Habitat == "Non-Antarctic sponge"])

# Create dataframe:
df.nas <- data.frame(Habitat = rep(c("Antarctic sponge","Non-Antarctic sponge"), times = c(73, 90)), Abundance = c(na.ant, na.noa))

# T test:
t.test(Abundance ~ Habitat, data = df.nas) #paired = F

#Plot:
nas.phylum.plot <- ggplot(data = df.nas, aes(x = Habitat, y = Abundance, fill = Habitat)) +
  geom_boxplot() +
  scale_fill_manual(values = colors.hab) +
  geom_jitter(alpha = 0.3, width = 0.20, size = 1.1) +
  labs(x = NULL, y = "Relative abundance") +
  theme_bw() +
  theme(legend.position = "none") +
  stat_compare_means(method = "t.test",
                     label = c("p.format"), label.x.npc = "center",
                     label.y.npc = "top", hide.ns = FALSE)

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/nas_phylum.png", units = "in", width = 5, height = 5, res = 300)
nas.phylum.plot
dev.off()


### Estimating cyanobacterias (bonus) ----

#--- Subset cyanos:
physeq.cyano <- subset_taxa(pqs, Phylum == "Cyanobacteria")

#--- Estimate total abundance of cyanobacterias:
count.cyano <- sum(matrix(taxa_sums(physeq.cyano))[,1])
sum(matrix(taxa_sums(pqs))[,1])

#--- Estimate abundance of cyanobacterias per sample:
abund.o <- as.matrix(sample_sums(pqs))
abund.cyano <- as.matrix(sample_sums(physeq.cyano))

# In percentage:
tmp1 <- sum(abund.o[,1])
tmp2 <- sum(abund.cyano[,1])
tmp2*100/tmp1

#--- Make table of % cyanobacteria in each sample:
tab.cyano <- cbind(abund.o, abund.cyano)
colnames(tab.cyano)[1:2] <- c("Total","Cyano")
tab.cyano <- cbind(tab.cyano, as.numeric(sprintf(tab.cyano[,2]*100/tab.cyano[,1], fmt = "%#.4f")))
colnames(tab.cyano)[3] <- c("% cyano")

#sum(tab.cyano[,3])/145

#--- Prune taxa:
physeq.nocyano <- subset_taxa(pqs, Phylum != "Cyanobacteria")

#--- Transform to relative abundance:
physeq.nocya.rel <- transform_sample_counts(physeq.nocyano, function(x) x/sum(x))

#--- Normalize data by Hellinger transformation:
physeq.nocya.hell <- transform_sample_counts(physeq.nocya.rel, function(x) sqrt(x))

#---Make ordination plots with Bray-Curtis method:
set.seed(1000)
bray.nocya.ord <- ordinate(physeq.nocya.hell, method = "NMDS", try = 100, trymax = 1000,
                           k = 2, distance = "bray", pc = TRUE, maxit = 1000)
bray.nocya.plot <- plot_ordination(physeq.nocya.hell, bray.nocya.ord, type = "samples",
                                   axes = 1:2, color = "Habitat", label = NULL,
                                   title = NULL, justDF = FALSE, shape = "Environment") +
  geom_point(size = 1) +
  scale_color_manual(values = colors.hab) +
  #labs(shape = ("Sample type")) +
  theme_bw() +
  #scale_shape_manual(values = c(17, 16)) +
  #annotation_custom(textGrob(label = "stress = 0.154", x = 0.8, y = 0.95, hjust = 0)) +
  theme(legend.position = "right", legend.title = element_blank())

# How many cyanos in pqs.rel:
table(tax_table(pqs.rel)[,"Phylum"] %in% "Cyanobacteria")


#--- Save things -----------------------------------------------

# Save objects:
saveRDS(pqs.hell, file = "./rds-files/pqs_hell.rds")

# Save the current work and objects:
save.image("R-graph-taxc.RData")

