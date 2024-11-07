# ANALYZING BETA DIVERSITY #
# ------------------------ #

#--- Load required objects ----

pqs <- readRDS(file = "./rds-files/pqs.rds")
pqs.rel <- readRDS(file = "./rds-files/pqs_rel.rds")
pqs.hell <- readRDS(file = "./rds-files/pqs_hell.rds")
metadata <- readRDS(file = "./rds-files/metadata.rds")


#--- Define colors for graphics ----
colors.hab <- as.matrix(jcolors("rainbow"))[c(6, 9)]
colors.hab <- rev(brewer.pal(n = 11, name = "RdBu"))[c(3,9)]
color.latzone <- c(Polar = "#2166ac", Temperate = "#DD8B71", Tropical = "#b2182b")
colors.1 <- c("darkorchid4", "deepskyblue2", "coral3", "tomato4", 
              "turquoise", "darkgoldenrod3", "lightpink", "#2166ac", "seagreen4",
              "limegreen", "orchid1",  "olivedrab4", "red1", "skyblue3", "maroon4",
              "darkorange", "gold2", "khaki", "plum4", "maroon2", "lightgreen",
              "aquamarine", "grey0", "mediumpurple", "sienna4", "tan", "darkcyan",
              "yellow", "thistle", "gray81")#darkblue


### Ordination of community composition from Bray-Curtis dissimilarities ----

#--- Agglomerate to the different taxa levels:
pqs.genus <- tax_glom(pqs, taxrank = "Genus", NArm=FALSE,
                      bad_empty = c(NA, "", " ", "\t"))
pqs.family <- tax_glom(pqs, taxrank = "Family", NArm=FALSE,
                      bad_empty = c(NA, "", " ", "\t"))
pqs.order <- tax_glom(pqs, taxrank = "Order", NArm=FALSE,
                       bad_empty = c(NA, "", " ", "\t"))
pqs.class <- tax_glom(pqs, taxrank = "Class", NArm=FALSE,
                      bad_empty = c(NA, "", " ", "\t"))
pqs.phylum <- tax_glom(pqs, taxrank = "Phylum", NArm=FALSE,
                      bad_empty = c(NA, "", " ", "\t"))

#--- Transform data by Hellinger transformation:
pqs.genus.rel <- transform_sample_counts(pqs.genus, function(x) x/sum(x))
pqs.genus.hell <- transform_sample_counts(pqs.genus.rel, function(x) sqrt(x))

pqs.family.rel <- transform_sample_counts(pqs.family, function(x) x/sum(x))
pqs.family.hell <- transform_sample_counts(pqs.family.rel, function(x) sqrt(x))

pqs.order.rel <- transform_sample_counts(pqs.order, function(x) x/sum(x))
pqs.order.hell <- transform_sample_counts(pqs.order.rel, function(x) sqrt(x))

pqs.class.rel <- transform_sample_counts(pqs.class, function(x) x/sum(x))
pqs.class.hell <- transform_sample_counts(pqs.class.rel, function(x) sqrt(x))

pqs.phylum.rel <- transform_sample_counts(pqs.phylum, function(x) x/sum(x))
pqs.phylum.hell <- transform_sample_counts(pqs.phylum.rel, function(x) sqrt(x))

#--- Set a seed because the starting position in NMDS algorithm is random:
set.seed(1000)

#--- Calculate ordination from dissimilarities:
bray.ord <- ordinate(pqs.hell, method = "NMDS", try = 100, trymax = 1000,
                       k = 2, distance = "bray", pc = TRUE, maxit = 1000) #autotransform = FALSE)
bray.genus.ord <- ordinate(pqs.genus.hell, method = "NMDS", try = 100, trymax = 1000,
                           k = 3, distance = "bray", pc = TRUE, maxit = 1000)
bray.family.ord <- ordinate(pqs.family.hell, method = "NMDS", try = 100, trymax = 1000,
                           k = 3, distance = "bray", pc = TRUE, maxit = 1000)
bray.order.ord <- ordinate(pqs.order.hell, method = "NMDS", try = 100, trymax = 1000,
                           k = 3, distance = "bray", pc = TRUE, maxit = 1000)
bray.class.ord <- ordinate(pqs.class.hell, method = "NMDS", try = 100, trymax = 1000,
                           k = 3, distance = "bray", pc = TRUE, maxit = 1000)
bray.phylum.ord <- ordinate(pqs.phylum.hell, method = "NMDS", try = 100, trymax = 1000,
                           k = 3, distance = "bray", pc = TRUE, maxit = 1000)


### Plotting community composition from Bray-Curtis dissimilarities ----

#--- Plot NMDS comparing environments and habitats at ASV level:
bray.plot.envs <- plot_ordination(pqs.hell, bray.ord, type = "samples", axes = 1:2,
                                 color = "Environment", shape = "Habitat",
                                 title = NULL, justDF = FALSE, label = NULL) +
  geom_point(size = 3) +
  scale_color_manual(values = color.latzone, 
                     labels = c("Antarctic", "Temperate", "Tropical")) +
  labs(shape = "Habitat", color = "Environment", title = NULL) +
       # caption = paste0("PERMANOVA p=0.001, R2=0.1105 |",
       #                  " ANOSIM p=", round(anosim.bray$signif, 4), ", R=", round(anosim.bray$statistic, 4))) +
  theme_bw() +
  #scale_shape_manual(values = c(17, 16)) +
  annotation_custom(textGrob(label = paste0("stress = ", round(bray.ord$stress, 3)), x = 0.7, y = 0.95, hjust = 0)) +
  theme(legend.position = "bottom", legend.title = element_text(size = 16, face = "bold"),
        aspect.ratio = 1,
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.text = element_text(size = 16),
        plot.caption = element_text(size = 16, face = "bold"),
        legend.box = "vertical") +
  guides(color = guide_legend(override.aes = list(size = 5)),#text = element_text(size = 20)
         shape = guide_legend(override.aes = list(size = 5))) +
stat_ellipse(level = 0.8, linetype = 1, aes(group = Environment), show.legend = F) # type = "norm"

#--- Plot NMDS comparing ecoregions at ASV level:
bray.plot.eco <- plot_ordination(pqs.hell, bray.ord, type = "samples", axes = 1:2,
                                 color = "Marine_Ecoregion", shape = "Habitat",
                                 title = NULL, justDF = FALSE, label = NULL) +
  geom_point(size = 4) +
  scale_color_manual(values = colors.1) +
  #labs(color = "Habitat") + # no me funciona
  labs(shape = "Habitat", color = "Marine_Ecoregion") +
       # caption = paste0("PERMANOVA p=0.001, R2=0.0949 |",
       #                  " ANOSIM p=", round(anosim.bray$signif, 4), ", R=", round(anosim.bray.prov$statistic, 4))) +
  theme_bw() +
  scale_shape_discrete(labels = c("Antarctic sponge", "Non-Antarctic sponge")) +
                     # values = c(17, 16)) +
  annotation_custom(textGrob(label = paste0("stress = ", round(bray.ord$stress, 3)), x = 0.7, y = 0.95, hjust = 0)) +
  theme(legend.position = "right", legend.title = element_text(size = 15, face = "bold"),
        aspect.ratio = 1,
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.text = element_text(size = 16),
        plot.caption = element_text(size = 16, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 5)),#text = element_text(size = 20)
         shape = guide_legend(override.aes = list(size = 5)))
#stat_ellipse(level = 0.95, linetype = 1, aes(group = Habitat), show.legend = F) # type = "norm"

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/beta-envs-bray.png", units = "in", width = 7, height = 7, res = 300)
bray.plot.envs
dev.off()

png("./results/graphics/beta-eco-bray.png", units = "in", width = 12, height = 7.5, res = 300)
bray.plot.eco
dev.off()

#--- Plot NMDS comparing envs and habitats for distinct taxa levels (NOTE: CHANGE TAX RANKS):
bray.taxa.plot <- plot_ordination(pqs.phylum.hell, bray.phylum.ord, type = "samples", axes = 1:2,
                                 color = "Environment", shape = "Habitat",
                                 title = NULL, justDF = FALSE, label = NULL) +
  geom_point(size = 2) +
  scale_color_manual(values = color.latzone,
                     labels = c("Antarctic", "Temperate", "Tropical")) +
  #labs(color = "Habitat") + # no me funciona
  labs(shape = "Habitat", color = "Environment", title = NULL) +
       # caption = paste0("ANOSIM p=", round(anosim.phylum.bray$signif, 4), ", R=", 
       #                  round(anosim.phylum.bray$statistic, 4))) +
  theme_bw() +
  #scale_shape_manual(values = c(17, 16)) +
  annotation_custom(textGrob(label = paste0("stress = ", round(bray.phylum.ord$stress, 3)), x = 0.7, y = 0.95, hjust = 0)) +
  theme(legend.position = "right", legend.title = element_text(size = 18, face = "bold"),
        aspect.ratio = 1,
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.text = element_text(size = 16),
        plot.caption = element_text(size = 16, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 5)),#text = element_text(size = 20)
         shape = guide_legend(override.aes = list(size = 5)))

#--- Save HQ image (NOTE: CHANGE TAX RANK):
png("./results/graphics/beta-phylum-bray.png", units = "in", width = 8.5, height = 5, res = 300)
bray.taxa.plot
dev.off()

#--- Plot NMDS comparing sequencing runs at ASV level:
bray.plot.sr <- plot_ordination(pqs.hell, bray.ord, type = "samples", axes = 1:4,
                                color = "Sequencing_run_ID", label = NULL,
                                title = NULL, justDF = FALSE) +# shape = "Sample_type"
  geom_point(size = 1) +
  scale_color_manual(values = colors.1) +
  labs(color = "Sequencing run") +
  theme_bw() +
  theme(legend.position = "right", aspect.ratio = 1)

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/beta-sr-bray.png", units = "in", width = 8, height = 8, res = 300)
bray.plot.sr
dev.off()

#--- Plot NMDS comparing host status at ASV level:

# Replace NAs by a character:
sample_data(pqs.hell)$Host_status <- as.character(sample_data(pqs.hell)$Host_status)
sample_data(pqs.hell)$Host_status[is.na(sample_data(pqs.hell)$Host_status)] <- "N/A"

# Reorder levels:
sample_data(pqs.hell)$Host_status <- factor(sample_data(pqs.hell)$Host_status,
                                               levels = c("HMA", "LMA", "N/A"))

# Make ordination plot:
bray.plot.ht <- plot_ordination(pqs.hell, bray.ord, type = "samples", axes = 1:4,
                                color = "Host_status", label = NULL, title = NULL,
                                justDF = FALSE) +
  geom_point(size = 1) +
  scale_color_brewer(palette = "Set2", direction = 1) +
  #labs(shape = ("Sample type")) +
  theme_bw() +
  theme(legend.position = "right", legend.title = element_blank(),
        aspect.ratio = 1)

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/beta-ht-bray.png", units = "in", width = 8, height = 8, res = 300)
bray.plot.ht
dev.off()

#--- Plot NMDS comparing sponge family at ASV level:
bray.plot.or <- plot_ordination(pqs.hell, bray.ord, type = "samples", axes = 1:4,
                                color = "Sponge_Family", label = NULL, title = NULL,
                                justDF = FALSE) +
  geom_point(size = 1) +
  scale_color_manual(values = colors.1) +
  #labs(shape = ("Sample type")) +
  theme_bw() +
  #scale_shape_manual(values = c(17, 16)) +
  #annotation_custom(textGrob(label = "stress = 0.1180", x = 0.8, y = 0.95, hjust = 0)) +
  theme(legend.position = "right", legend.title = element_blank(),
        aspect.ratio = 1)

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/beta-fam-bray.png", units = "in", width = 8, height = 8, res = 300)
bray.plot.or
dev.off()

### Ordination of composition from UniFrac distances ----

#--- Set a seed because the starting position in NMDS algorithm is random:
set.seed(1000)

#--- Calculate ordination from dissimilarities:
wunif.ord <- ordinate(pqs.hell, method = "NMDS", distance = "unifrac",
                          weighted = FALSE)
wunifw.ord <- ordinate(pqs.hell, method = "NMDS", distance = "unifrac",
                          weighted = TRUE)


### Plotting community composition from UniFrac distances ----

#--- Plot NMDS comparing environments and habitats at ASV level:
unifrac.plot <- plot_ordination(pqs.hell, wunif.ord, type = "samples",
                                color = "Environment", shape = "Habitat", label = NULL) +
  geom_point(size = 3) +
  scale_color_manual(values = color.latzone,
                     labels = c("Antarctic", "Temperate", "Tropical")) +
  #labs(shape = "Habitat", fill = "Zone") +
  labs(shape = "Habitat", color = "Environment", title = NULL) +#"Unweighted UniFrac distance"
  # caption = paste0("PERMANOVA p=0.001, R2=0.0755 |",
  #   " ANOSIM p=", round(anosim.wunif$signif, 4), ", R=", round(anosim.wunif$statistic, 4))) +
  theme_bw() +
  #scale_shape_manual(values = c(17, 16)) +
  annotation_custom(textGrob(label = paste0("stress = ", round(wunif.ord$stress, 3)), x = 0.7, y = 0.95, hjust = 0)) +
  theme(legend.position = "bottom", legend.title = element_text(size = 18, face = "bold"),
        aspect.ratio = 1,
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        plot.caption = element_text(size = 17, face = "bold"),
        legend.box = "vertical") +
  guides(color = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5)))

#--- Plot NMDS comparing ecoregions at ASV level:
unif.plot.eco <- plot_ordination(pqs.hell, wunif.ord, type = "samples",
                                 color = "Marine_Ecoregion", shape = "Environment", label = NULL) +
  geom_point(size = 2) +
  scale_color_manual(values = colors.1) +
  #labs(shape = "Habitat", fill = "Zone") +
  labs(shape = "Environment", color = "Marine_Ecoregion") + 
  # caption = paste0("PERMANOVA p=0.001, R2=0.0755 |",
  #                  " ANOSIM p=", round(anosim.wunif$signif, 4), ", R=", round(anosim.wunif$statistic, 4))) +
  theme_bw() +
  scale_shape_discrete(labels = c("Antarctic", "Temperate", "Tropical")) +
  # values = c(17, 16)) +
  annotation_custom(textGrob(label = paste0("stress = ", round(wunif.ord$stress, 3)), x = 0.7, y = 0.95, hjust = 0)) +
  theme(legend.position = "right", legend.title = element_text(size = 18, face = "bold"),
        aspect.ratio = 1,
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        plot.caption = element_text(size = 17, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5)))

#stat_ellipse(level = 0.95, linetype = 1, aes(group = Habitat), show.legend = F) # type = "norm"

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/beta-envs-unif.png", units = "in", width = 7.5, height = 7.5, res = 300)
unifrac.plot
dev.off()

png("./results/graphics/beta-eco-unif.png", units = "in", width = 11, height = 11, res = 300)
unif.plot.eco
dev.off()

#--- Plot just the legend:
legend.beta <- cowplot::get_legend(plot_ordination(pqs.hell, wunif.ord, type = "samples",
                                                   color = "Environment", shape = "Habitat", label = NULL) +
                                     geom_point(size = 3) +
                                     scale_color_manual(values = color.latzone,
                                                        labels = c("Antarctic", "Temperate", "Tropical")) +
                                     #labs(shape = "Habitat", fill = "Zone") +
                                     labs(shape = "Habitat", color = "Environment", title = "Unweighted UniFrac distance") +
                                     # caption = paste0("PERMANOVA p=0.001, R2=0.0755 |",
                                     #   " ANOSIM p=", round(anosim.wunif$signif, 4), ", R=", round(anosim.wunif$statistic, 4))) +
                                     theme_bw() +
                                     #scale_shape_manual(values = c(17, 16)) +
                                     annotation_custom(textGrob(label = paste0("stress = ", round(wunif.ord$stress, 3)), x = 0.7, y = 0.95, hjust = 0)) +
                                     theme(legend.position = "bottom", legend.title = element_text(size = 18, face = "bold"),
                                           aspect.ratio = 1,
                                           text = element_text(size = 15),
                                           plot.title = element_text(hjust = 0.5, size = 20),
                                           legend.text = element_text(size = 20),
                                           plot.caption = element_text(size = 17, face = "bold")) +
                                     guides(color = guide_legend(override.aes = list(size = 5)),
                                            shape = guide_legend(override.aes = list(size = 5))))
legend.beta <- as_ggplot(legend.beta)

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/beta-legend.png", units = "in", width = 15, height = 1.5, res = 300)
legend.beta
dev.off()


### Calculating distance matrices ----

#--- Set a seed:
set.seed(1000)

#--- Calculate distance matrixes with Bray-Curtis and UniFrac distances:

# At ASV level:
bray.dist <- phyloseq::distance(pqs.hell, method = "bray")
wunif.dist <- phyloseq::distance(pqs.hell, method = "unifrac")

# At different taxonomic levels:
bray.genus.dist <- phyloseq::distance(pqs.genus.hell, method = "bray")
bray.family.dist <- phyloseq::distance(pqs.family.hell, method = "bray")
bray.order.dist <- phyloseq::distance(pqs.order.hell, method = "bray")
bray.class.dist <- phyloseq::distance(pqs.class.hell, method = "bray")
bray.phylum.dist <- phyloseq::distance(pqs.phylum.hell, method = "bray")

wu.genus.dist <- phyloseq::distance(pqs.genus.hell, method = "unifrac")
wu.family.dist <- phyloseq::distance(pqs.family.hell, method = "unifrac")
wu.order.dist <- phyloseq::distance(pqs.order.hell, method = "unifrac")
wu.class.dist <- phyloseq::distance(pqs.class.hell, method = "unifrac")
wu.phylum.dist <- phyloseq::distance(pqs.phylum.hell, method = "unifrac")


### Performing PERMANOVA among groups ----

#--- Prepare metadata:
# metadata <- as(sample_data(pqs), "data.frame")
# metadata.genus <- as(sample_data(pqs.genus), "data.frame")

#--- Add column of primer set in metadata:
adiv.mdata$Primer_set_ID <- c(rep("515F-926R_Parada", 33), 
                              rep("515F-806R_Caporaso", 114),
                              rep("515F-806R_Parada", 32))                                           

#--- Set a seed:
set.seed(1000)

# From Bray-Curtis:
permanova.bray <- adonis2(bray.dist ~ Habitat*Sponge_Family + Primer_set_ID,
                          data = metadata,
                          permutations = 999,
                          strata = metadata$Sequencing_run_ID)

set.seed(1000)                                           
permanova.bray.eco <- adonis2(bray.dist ~ Marine_Ecoregion*Sponge_Family + Primer_set_ID,
                          data = metadata,
                          permutations = 999,
                          strata = metadata$Sequencing_run_ID)

# From unweighted UniFrac:
set.seed(1000)                                           
permanova.wunif <- adonis2(wunif.dist ~ Habitat*Sponge_Family + Primer_set_ID,
                           data = metadata,
                           permutations = 999,
                           strata = metadata$Sequencing_run_ID)

set.seed(1000)                                           
permanova.wunif.eco <- adonis2(wunif.dist ~ Marine_Ecoregion*Sponge_Family + Primer_set_ID,
        data = metadata,
        permutations = 999,
        strata = metadata$Sequencing_run_ID)


### Performing pairwise PERMANOVA ----

#--- Set a seed:
set.seed(1000)

# For habitat:
pairw.bray <- pairwise.adonis(bray.dist,
                              factors = metadata$Habitat,
                              p.adjust.m = "none",
                              perm = 999) #, reduce = "Antarctic", also "bonferroni" in p.adjust

pairw.wunif <- pairwise.adonis(wunif.dist,
                               factors = metadata$Habitat,
                               p.adjust.m = "none",
                               perm = 999) #, reduce = "Antarctic", also "bonferroni" in p.adjust

# For sponge family:
pairw.bray.f <- pairwise.adonis(bray.dist,
                                factors = metadata$Sponge_Family,
                                p.adjust.m = "none",
                                perm = 999)

# For ecoregion:
pairw.bray.e <- pairwise.adonis(bray.dist,
                                factors = metadata$Marine_Ecoregion,
                                p.adjust.m = "none",
                                perm = 999)


### Testing homogeneity of dispersion (variances) PERMDISP ----

#--- Set a seed:
set.seed(1000)

# From Bray-Curtis:
bray.disp <- betadisper(bray.dist, metadata$Habitat)
homog.bray <- permutest(bray.disp, permutations = 999)

# bray.disp2 <- betadisper(bray.dist, metadata$Sponge_Family)
# homog.bray2 <- permutest(bray.disp2, permutations = 999)
# 
# bray.disp3 <- betadisper(bray.dist, metadata$Marine_Ecoregion)
# homog.bray3 <- permutest(bray.disp3, permutations = 999)

# From unweighted UniFrac:
wunif.disp <- betadisper(wunif.dist, metadata$Habitat)
homog.wunif <- permutest(wunif.disp, permutations = 999)


### Performing ANOSIM ----

#--- Set a seed:
set.seed(1000)

# From Bray-Curtis:
anosim.bray <- anosim(x = bray.dist, grouping = metadata$Habitat,
                      permutations = 999)

anosim.bray.eco <- anosim(x = bray.dist, grouping = metadata$Marine_Ecoregion,
                          permutations = 999)

anosim.genus.bray <- anosim(x = bray.genus.dist, grouping = metadata$Habitat,
                            permutations = 999)
anosim.family.bray <- anosim(x = bray.family.dist, grouping = metadata$Habitat,
                             permutations = 999)
anosim.order.bray <- anosim(x = bray.order.dist, grouping = metadata$Habitat,
                            permutations = 999)
anosim.class.bray <- anosim(x = bray.class.dist, grouping = metadata$Habitat,
                            permutations = 999)
anosim.phylum.bray <- anosim(x = bray.phylum.dist, grouping = metadata$Habitat,
                             permutations = 999)

# From unweighted UniFrac:
anosim.wunif <- anosim(x = wunif.dist, grouping = metadata$Habitat,
                       permutations = 999)

anosim.genus.unif <- anosim(x = wu.genus.dist, grouping = metadata$Habitat,
                            permutations = 999)
anosim.family.unif <- anosim(x = wu.family.dist, grouping = metadata$Habitat,
                             permutations = 999)
anosim.order.unif <- anosim(x = wu.order.dist, grouping = metadata$Habitat,
                            permutations = 999)
anosim.class.unif <- anosim(x = wu.class.dist, grouping = metadata$Habitat,
                            permutations = 999)
anosim.phylum.unif <- anosim(x = wu.phylum.dist, grouping = metadata$Habitat,
                             permutations = 999)


### Performing SIMPER ----

#--- Join taxa names:

# From ASV phyloseq object:
taxa_names(pqs.rel) <- paste(as.data.frame(tax_table(pqs.rel))$Class,
                             as.data.frame(tax_table(pqs.rel))$Order,
                             as.data.frame(tax_table(pqs.rel))$Family,
                             as.data.frame(tax_table(pqs.rel))$Genus,
                             rownames(tax_table(pqs.rel)),
                             sep = "; ")

# From genus-agglomerated phyloseq object:
taxa_names(pqs.genus.rel) <- paste(as.data.frame(tax_table(pqs.genus.rel))$Class,
                                   as.data.frame(tax_table(pqs.genus.rel))$Order,
                                   as.data.frame(tax_table(pqs.genus.rel))$Family,
                                   as.data.frame(tax_table(pqs.genus.rel))$Genus,
                                   rownames(tax_table(pqs.genus.rel)),
                                   sep = "; ")

#--- Run SIMPER (RAM expensive):
simper.asv <- simper(otu_table(pqs), sample_data(pqs)$Habitat) #pqs
simper.asv <- simper(otu_table(pqs.rel), sample_data(pqs.rel)$Habitat) #pqs.rel
simper.genus <- simper(otu_table(pqs.genus.rel), sample_data(pqs.genus.rel)$Habitat) #pqs.genus.rel

#--- Prepare the resulting table: 
simper.asv.tab <- as.data.frame(summary(simper.asv)$`Antarctic sponge_Non-Antarctic sponge` %>%
                   round(3) %>%
                   head(100))

#--- Export table:
write.table(simper.asv.tab, file ="./results/simper-100first.csv",
            sep="\t", na = "NA", row.names = T, col.names = NA)


### Performing LDA (Genera differentially abundant) ----

#--- Set a seed:
set.seed(1000)

#--- Run LDA:
pqs.lefse <- run_lefse(pqs, wilcoxon_cutoff = 0.01, group = "Habitat",
                       multigrp_strat = F, taxa_rank = "Genus", norm = "none",
                       lda_cutoff = 3)

# pqs.lefse.g <- run_lefse(pqs.genus, wilcoxon_cutoff = 0.01, group = "Habitat",
#                        multigrp_strat = F, taxa_rank = "none", norm = "none",
#                        lda_cutoff = 3)
# 
# pqs.lefse2 <- run_lefse(pqs, wilcoxon_cutoff = 0.01, group = "Habitat",
#                        multigrp_strat = F, taxa_rank = "none", norm = "none",
#                        lda_cutoff = 3)

#--- Inspect results:
tail(marker_table(pqs.lefse),10)
plot_ef_bar(pqs.lefse)
marker_table(pqs.lefse)$enrich_group <- factor(as.factor(marker_table(pqs.lefse)$enrich_group),
                                               levels = c("Non-Antarctic sponge", "Antarctic sponge"))

#--- Customize format of results:
lda.features <- marker_table(pqs.lefse)$feature
lda.features[c(1,5,6,7,9,15,16,18,19,20)] <- c("Unassigned Gammaproteobacteria",
                                               "Unassigned Burkholderiales",
                                               "Unassigned Flavobacteriaceae",
                                               "Unassigned Halieaceae",
                                               "Unassigned Nitrincolaceae",
                                               "UBA10353 marine group",
                                               "Unnasigned Bacteria",
                                               "Unassigned Alphaproteobacteria",
                                               "Unassigned Nitrosopumilaceae",
                                               "Unassigned Proteobacteria")

#--- Plot abundances of each group in boxplots:
# abund.plot <- plot_abundance(pqs.lefse, group = "Habitat")

#--- Plot LDA:
lda.plot <- plot_ef_bar(pqs.lefse) +
  scale_fill_manual(values = rev(colors.hab)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), legend.position = "none",
        legend.title = element_text(size = 13, face = "bold"),
              text = element_text(size = 15)) +
  scale_y_discrete(labels = rev(lda.features))

#--- Plot legend:
legend.lda <- cowplot::get_legend(lda.plot)
legend.lda <- as_ggplot(legend.lda)

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/lda-hab.png", units = "in", width = 6, height = 6, res = 300)
lda.plot
dev.off()

png("./results/graphics/lda-leg.png", units = "in", width = 6, height = 1, res = 300)
legend.lda
dev.off()


### Making a heatmap based on Bray-Curtis dissimilarities (bonus) ----

#--- Fix matrix of distances:
bray.dist.mat <- as.matrix(bray.dist)
wunif.dist.mat <- as.matrix(wunif.dist)

#--- Isolate annotations for heatmap:
colnames(get_variable(pqs.hell))
tmp.data <- data.frame(sample_data(pqs.hell)[,6]) # Environment
tmp.data2 <- data.frame(sample_data(pqs.hell)[,7]) # Ecoregion

#--- Define colors for plotting:
library(RColorBrewer)
library(viridis)
tmp.color.latzone <- list(Environment =
  c(Antarctic = "#2166ac", Temperate = "#67a9cf", Tropical = "#b2182b"))
colors.eco <- list(Marine_Ecoregion =
                     c('Antarctic Peninsula' = "deepskyblue2", 
                       'South Shetland Islands' = "#2166ac",
                       'Adriatic Sea' = "darkorchid4",
                       'Bahamian' = "coral3",
                       'Central and Southern Great Barrier Reef' = "tomato4",
                       'Central New Zealand' = "turquoise",
                       'Hawaii' = "darkgoldenrod3",
                       'South European Atlantic Shelf' = "lightpink",
                       'Southern Norway' = "seagreen4",
                       'Southwestern Caribbean' = "limegreen",
                       'Western Caribbean' = "orchid1",
                       'Western Mediterranean' = "olivedrab4"))

#--- Plot heatmaps:

# For environment:
tmp.heatmap <- as.ggplot(pheatmap(bray.dist.mat, scale = "none", clustering_method = "ward.D",
                                  fontsize = 7,
                                  cluster_cols = T,
                                  cluster_rows = T,
                                  show_colnames = T, 
                                  na_col = "white",
                                  cutree_rows = 1,
                                  cutree_cols = 2,
                                  annotation_col = tmp.data,
                                  annotation_row = tmp.data,
                                  # legend = FALSE,
                                  annotation_colors = tmp.color.latzone,
                                  show_rownames = F,
                                  border_color = NULL,
                                  treeheight_row = 50,
                                  treeheight_col = 50,
                                  annotation_names_row = F,
                                  annotation_names_col = F,
                                  fontsize_col = 2,
                                  fontsize_row = 2,
                                  angle_col = 90,
                                  filename = "./results/graphics/heatmap.png",
                                  # color = colorRampPalette(colors = c())
                                  color = brewer.pal(8, "BuPu"),#RdPu
                                  legend = T,
                                  annotation_legend = T))
                                  #kmeans_k = 4))

# For ecoregion:
tmp.heatmap <- as.ggplot(pheatmap(bray.dist.mat, scale = "none", clustering_method = "ward.D",
                                  fontsize = 7,
                                  cluster_cols = T,
                                  cluster_rows = T,
                                  show_colnames = T, 
                                  na_col = "white",
                                  cutree_rows = 1,
                                  cutree_cols = 2,
                                  annotation_col = tmp.data2,
                                  annotation_row = tmp.data2,
                                  # legend = FALSE,
                                  annotation_colors = colors.eco,
                                  show_rownames = F,
                                  border_color = NULL,
                                  treeheight_row = 50,
                                  treeheight_col = 50,
                                  annotation_names_row = F,
                                  annotation_names_col = F,
                                  fontsize_col = 2,
                                  fontsize_row = 2,
                                  angle_col = 90,
                                  filename = "./results/graphics/heatmap-eco.png",
                                  # color = colorRampPalette(colors = c())
                                  color = brewer.pal(8, "BuPu"),#RdPu
                                  legend = T,
                                  annotation_legend = T))


#--- Export tables of PERMANOVA comparisons ----
write.table(permanova.bray, file ="./results/stats_bray_permanova.txt",
            sep="\t", na = "NA", row.names = T, col.names = NA)
write.table(permanova.bray.eco, file ="./results/stats_bray_permanova_eco.txt",
            sep="\t", na = "NA", row.names = T, col.names = NA)
write.table(permanova.wunif, file ="./results/stats_wunif_permanova.txt",
            sep="\t", na = "NA", row.names = T, col.names = NA)
write.table(pairw.bray, file ="./results/permanova_pw_habitat.txt",
            sep="\t", na = "NA", row.names = T, col.names = NA)
write.table(pairw.bray.f, file ="./results/permanova_pw_family.txt",
            sep="\t", na = "NA", row.names = T, col.names = NA)
write.table(pairw.bray.e, file ="./results/permanova_pw_ecoreg.txt",
            sep="\t", na = "NA", row.names = T, col.names = NA)


#--- Save things -----------------------------------------------

# Save objects:
saveRDS(bray.dist, file = "./rds-files/bray_dist.rds")
saveRDS(bray.ord, file = "./rds-files/bray_ord.rds")

# Save the current work and objects:
save.image("R-graph-bdiv.RData")

