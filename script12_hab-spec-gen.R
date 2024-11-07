# ASSESSING HABITAT-SPECIFIC AND HABITAT-GENERALIST MICROORGANISMS #
# ---------------------------------------------------------------- #

#--- Load required objects ----
pqs <- readRDS(file = "./rds-files/pqs.rds")
pqs.rel <- readRDS(file = "./rds-files/pqs_rel.rds")
pqs.hell <- readRDS(file = "./rds-files/pqs_hell.rds")
metadata <- readRDS(file = "./rds-files/metadata.rds")


#--- Define colors for graphics ----
#colors.hab <- as.matrix(jcolors("rainbow"))[c(6, 9)]
colors.hab <- rev(brewer.pal(n = 11, name = "RdBu"))[c(3,9)]
colors.env <- c("#2166ac", "#e0967f", "#b2182b")
color.latzone <- c(Polar = "#2166ac", Temperate = "#e0967f", Tropical = "#b2182b")
colors.1 <- c("darkorchid4", "deepskyblue2", "coral3", "tomato4", 
              "turquoise", "darkgoldenrod3", "lightpink", "#2166ac", "seagreen4",
              "limegreen", "orchid1",  "olivedrab4", "red1", "skyblue3", "maroon4",
              "darkorange", "gold2", "khaki", "plum4", "maroon2", "lightgreen",
              "aquamarine", "grey0", "mediumpurple", "sienna4", "tan", "darkcyan",
              "yellow", "thistle", "gray81")#darkblue
colors.eco <- c("#2166ac", "seagreen4", "lightpink", "darkgoldenrod3",
                "olivedrab4", "limegreen", "tomato4", "deepskyblue2",
                "coral3", "turquoise", "darkorchid4", "orchid1")
colors.eco2 <- c("#2166ac", "#2166ac", "#e0967f", "#e0967f",
                 "#e0967f", "#e0967f", "#e0967f", "#b2182b",
                 "#b2182b", "#b2182b", "#b2182b", "#b2182b")


### Determining exclusive and shared ASVs among ecoregions (Upset plot) ----

#--- Create upset tables (activate library "MicrobiotaProcess"):
tmp1 <- get_upset(obj = pqs, factorNames = "Environment")
tmp2 <- get_upset(obj = pqs, factorNames = "Marine_Ecoregion")
tmp3 <- get_upset(obj = pqs, factorNames = "Sponge_Family")

#--- Add other variables to the upset table:
tmp2$Occurrence <- taxa_prev(pqs.rel)
tmp2$Abundance <- colSums(otu_table(pqs.rel))/colSums(!!otu_table(pqs.rel))

#--- Change names of ecoregions (adding sample size):
colnames(tmp2)[1] <- "Adriatic Sea (3)"
colnames(tmp2)[2] <- "Antarctic Peninsula (12)"
colnames(tmp2)[3] <- "Bahamian (4)"
colnames(tmp2)[4] <- "Central and Southern Great Barrier Reef (10)"
colnames(tmp2)[5] <- "Central New Zealand (6)"
colnames(tmp2)[6] <- "Hawaii (13)"
colnames(tmp2)[7] <- "South European Atlantic Shelf (9)"
colnames(tmp2)[8] <- "South Shetland Islands (77)"
colnames(tmp2)[9] <- "Southern Norway (13)"
colnames(tmp2)[10] <- "Southwestern Caribbean (10)"
colnames(tmp2)[11] <- "Western Caribbean (3)"
colnames(tmp2)[12] <- "Western Mediterranean (19)"

# Customize order of ecoregions:
ecoregions <- colnames(tmp3)#unique(sort(sample_data(pqs)$Ecoregion))
ecoregions <- ecoregions[c(8,2,9,7,12,5,1,6,10,4,3,11)]

#--- Make upset plots:

# For environment:
upset.plot <- upset(tmp1, nsets = NA, nintersects = NA, sets = c("Antarctic", "Temperate", "Tropical"),
                    keep.order = T, set.metadata = NULL, intersections = NULL,
                    matrix.color = "gray23", main.bar.color = "gray23",
                    mainbar.y.label = "Exclusive/Shared ASVs", mainbar.y.max = NULL,
                    sets.bar.color = colors.zone[order(c(1,2,3))], sets.x.label = "# ASVs",
                    point.size = 4, line.size = 1, mb.ratio = c(0.7, 0.3),
                    expression = NULL, att.pos = NULL, att.color = NULL,
                    order.by = c("freq"), decreasing = c(T, F),
                    show.numbers = "yes", number.angles = 0, group.by = "degree",
                    cutoff = NULL, queries = NULL, query.legend = "none",
                    shade.color = "gray88", shade.alpha = 0.25, matrix.dot.alpha = 0.5,
                    empty.intersections = NULL, color.pal = 1, boxplot.summary = NULL,
                    attribute.plots = NULL, scale.intersections = "identity",
                    scale.sets = "identity", text.scale = 1.7, set_size.angles = 0,
                    set_size.show = T, set_size.numbers_size = 7,
                    set_size.scale_max = 15000)

# For ecoregion:
upset.plot.eco <- upset(tmp3, nsets = 12, nintersects = 25, sets = ecoregions,
                        keep.order = T, matrix.color = "gray23", main.bar.color = "gray23",
                        sets.bar.color = colors.eco2, sets.x.label = "# ASVs",
                        point.size = 2, line.size = 0.7, mb.ratio = c(0.7, 0.3),
                        order.by = c("freq"), decreasing = c(T, F),
                        show.numbers = "yes", number.angles = 0, group.by = "degree",
                        cutoff = NULL, query.legend = "left",
                        shade.color = "gray88", shade.alpha = 0.25, matrix.dot.alpha = 0.5,
                        empty.intersections = NULL, color.pal = 1, boxplot.summary = NULL,
                        scale.intersections = "identity",
                        scale.sets = "identity", text.scale = 1.2, set_size.angles = 0,
                        set_size.show = T, set_size.numbers_size = 6,
                        set_size.scale_max = 13000, mainbar.y.label = NULL, 
                        queries = list(list(query = intersects,
                                            params = list("South Shetland Islands (77)", "Antarctic Peninsula (12)"),
                                            color = "#4393c3",
                                            active = T),
                                       list(query = intersects,
                                            params = list("South Shetland Islands (77)"),
                                            color = "#4353c3", active = T),
                                       list(query = intersects,
                                            params = list("Antarctic Peninsula (12)"),
                                            color = "#43c3b3", active = T),
                                       list(query = intersects,
                                            params = list("South Shetland Islands (77)", "Hawaii (13)"),
                                            color = "darkgoldenrod3", active = T)),
                        attribute.plots = list(gridrows = 45, gridcols = 180,
                                               plots = list(list(plot = scatter_plot,
                                                                 x = "Abundance", y = "Occurrence",
                                                                 queries = T)), ncols = 2))

upset(tmp3, nsets = 12, nintersects = 25, sets = ecoregions,
      keep.order = T, boxplot.summary = c("abundance", "Occurrence"))

# upset.plot.fam <- upset(tmp3, nsets = 16, nintersects = 20, sets = NULL,
#                         keep.order = F, set.metadata = NULL, intersections = NULL,
#                         matrix.color = "gray23", main.bar.color = "gray23",
#                         mainbar.y.label = NULL, mainbar.y.max = NULL,
#                         sets.bar.color = colors.1[1:16], sets.x.label = "# ASVs",
#                         point.size = 2, line.size = 0.7, mb.ratio = c(0.7, 0.3),
#                         expression = NULL, att.pos = NULL, att.color = NULL,
#                         order.by = c("freq"), decreasing = c(T, F),
#                         show.numbers = "yes", number.angles = 0, group.by = "degree",
#                         cutoff = NULL, queries = NULL, query.legend = "none",
#                         shade.color = "gray88", shade.alpha = 0.25, matrix.dot.alpha = 0.5,
#                         empty.intersections = NULL, color.pal = 1, boxplot.summary = NULL,
#                         attribute.plots = NULL, scale.intersections = "identity",
#                         scale.sets = "identity", text.scale = 1.2, set_size.angles = 0,
#                         set_size.show = F, set_size.numbers_size = 8,
#                         set_size.scale_max = 13000)

# Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/upset-envs.png", units = "in", width = 14, height = 6, res = 300)
upset.plot
dev.off()

png("./results/graphics/upset-eco.png", units = "in", width = 12, height = 9, res = 300)
upset.plot.eco
dev.off()


### Determining exclusive and shared ASVs between habitats and environments ----

#--- Plot shared and exclusive ASVs:
venn.hab <- ps_venn(pqs.rel, group = "Habitat", fraction = 0, weight = FALSE, 
                     relative = FALSE, plot = TRUE, 
                     fill = colors.hab, 
                     legend = list(side = "right"))

venn.env <- ps_venn(pqs.rel, group = "Environment", fraction = 0, weight = FALSE, 
                     relative = FALSE, plot = TRUE, 
                     fill = color.latzone, 
                     legend = list(side = "right"))

#--- Save the lists of the shared and exclusive ASVs:
shared.taxa <- ps_venn(pqs.rel, group = "Habitat", fraction = 0, weight = FALSE, 
                       relative = FALSE, plot = FALSE)

#--- Other way, separate habitats:

# Antarctic sponge subset: 
pqs.ant <- subset_samples(pqs.rel, Habitat == "Antarctic sponge")
pqs.ant <- prune_taxa(taxa_sums(pqs.ant) > 0, pqs.ant)
asvs.ant <- taxa_names(pqs.ant)
seqtab.ant <- as.data.frame(otu_table(pqs.ant))
table(tax_table(pqs.ant)[,1])

# Non-Antarctic sponge subset: 
pqs.noant <- subset_samples(pqs.rel, Habitat == "Non-Antarctic sponge")
pqs.noant <- prune_taxa(taxa_sums(pqs.noant) > 0, pqs.noant)
asvs.noant <- taxa_names(pqs.noant)
seqtab.noant <- as.data.frame(otu_table(pqs.noant))
table(tax_table(pqs.noant)[,1])

# Join in a list:
asvs <- list("Antarctic" = asvs.ant, "Non-Antarctic" = asvs.noant)

# Plot shared and exclusive ASVs:
venn.hab2 <- ggvenn(rev(asvs), fill_color = rev(colors.hab), 
                     fill_alpha = 0.6, stroke_size = 0.5,
       set_name_size = 6, auto_scale = TRUE) +
  # labs(title = "Shared and exclusive ASVs") +
  theme(plot.title = element_text(hjust = 0.5)) #option1

venn.env2 <- ggvenn(asvs2, fill_color = rev(colors.zone), 
                     fill_alpha = 0.5, stroke_size = 0.6,
                     set_name_size = 6,
                     text_size = 5, show_percentage = F) +
  # labs(title = "Shared and exclusive ASVs") +
  theme(plot.title = element_text(hjust = 0.5)) #option1

draw.pairwise.venn(area1 = length(asvs.ant),
                   area2 = length(asvs.noant),
                   cross.area = length(shared.taxa$`Antarctic sponge__Non-Antarctic sponge`))


### Determining habitat specific and generalist ASVs ----

# Antarctic sponge specific ASVs:
ant.exc.list <- shared.taxa$`Antarctic sponge`
pqs.ant.exc <- prune_taxa(ant.exc.list, pqs.rel)
pqs.ant.exc <- subset_samples(pqs.ant.exc, Habitat == "Antarctic sponge")
pqs.ant.exc <- prune_taxa(taxa_names(pqs.ant.exc) != "ASV14533", pqs.ant.exc) #ciano
pqs.ant.exc <- prune_taxa(taxa_names(pqs.ant.exc) != "ASV20423", pqs.ant.exc) #ciano
filt.ant.50 <- genefilter_sample(pqs.ant.exc, 
                                 filterfun_sample(function(x) x > 0), 
                                 A = 0.5*nsamples(pqs.ant.exc))
pqs.ant.50 <- prune_taxa(filt.ant.50, pqs.ant.exc)
pqs.ant.50 <- prune_samples(sample_sums(pqs.ant.50) > 0, pqs.ant.50)
seqtab.ant50 <- as.data.frame(otu_table(pqs.ant.50))
spec.asvs <- taxa_names(pqs.ant.50)
tax.ant.spec <- tax_table(pqs.ant.50)
mdata.ant.spec <- sample_data(pqs.ant.50)

write.table(tax.ant.spec, file = "./results/taxonomy_ant_specific.tsv", 
            sep = "\t", na = "NA", row.names = T, col.names = NA)
write.table(taxa_names(pqs.ant.50), file = "./results/asvs_ant_specific.txt", 
            sep = "\t", na = "NA", row.names = F, col.names = F)

# Non-Antarctic sponge specific ASVs:
noant.exc.list <- shared.taxa$`Non-Antarctic sponge`
pqs.noant.exc <- prune_taxa(noant.exc.list, pqs.rel)
pqs.noant.exc <- subset_samples(pqs.noant.exc, Habitat == "Non-Antarctic sponge")
filt.noant <- genefilter_sample(pqs.noant.exc, 
                               filterfun_sample(function(x) x > 0), 
                               A = 0.5*nsamples(pqs.noant.exc))
pqs.noant50 <- prune_taxa(filt.noant, pqs.noant.exc)
sample_sums(pqs.noant50)
pqs.noant50 <- prune_samples(sample_sums(pqs.noant50) > 0, pqs.noant50)
seqtab.noant <- as.data.frame(otu_table(pqs.noant50))
spec.noan.asvs <- taxa_names(pqs.noant50)
tax.noa.spec <- tax_table(pqs.noant50)
mdata.noa.spec <- sample_data(pqs.noant50)

write.table(tax.noa.spec, file = "./results/taxonomy_noant_specific.tsv", 
            sep = "\t", na = "NA", row.names = T, col.names = NA)
write.table(taxa_names(pqs.noant50), file = "./results/asvs_noant_specific.txt", 
            sep = "\t", na = "NA", row.names = F, col.names = F)

# Sponge generalist ASVs:
shared.list <- shared.taxa$`Antarctic sponge__Non-Antarctic sponge`
pqs.shared <- prune_taxa(shared.list, pqs.rel)
filt.gen <- genefilter_sample(pqs.shared, 
                                  filterfun_sample(function(x) x > 0), 
                                  A = 0.5*nsamples(pqs.shared))
pqs.gen <- prune_taxa(filt.gen, pqs.shared)
pqs.gen <- prune_samples(sample_sums(pqs.gen) > 0, pqs.gen)
seqtab.gen <- as.data.frame(otu_table(pqs.gen))
gen.asvs <- taxa_names(pqs.gen)
tax.gen <- tax_table(pqs.gen)
mdata.gen <- sample_data(pqs.gen)
sample_sums(pqs.gen)

write.table(tax.gen, file = "./results/taxonomy_generalists.tsv", 
            sep = "\t", na = "NA", row.names = T, col.names = NA)
write.table(taxa_names(pqs.gen), file = "./results/asvs_generalist.txt", 
            sep = "\t", na = "NA", row.names = F, col.names = F)


### Determining the presence of ASVs in Antarctic and non-Antarctic sponges ----

#--- Apply function for filtering by prevalence/presence across sponges:

# For Antarctic sponges only:
filter.pro33 <- genefilter_sample(pqs.ant.exc, filterfun_sample(function(x) x > 0), 
                                  A = 0.33*nsamples(pqs.ant.exc))
filter.pro50 <- genefilter_sample(pqs.ant.exc, filterfun_sample(function(x) x > 0), 
                                  A = 0.5*nsamples(pqs.ant.exc))
filter.pro60 <- genefilter_sample(pqs.ant.exc, filterfun_sample(function(x) x > 0), 
                                  A = 0.6*nsamples(pqs.ant.exc))
filter.pro70 <- genefilter_sample(pqs.ant.exc, filterfun_sample(function(x) x > 0), 
                                  A = 0.7*nsamples(pqs.ant.exc))
filter.pro80 <- genefilter_sample(pqs.ant.exc, filterfun_sample(function(x) x > 0), 
                                  A = 0.8*nsamples(pqs.ant.exc))
filter.pro90 <- genefilter_sample(pqs.ant.exc, filterfun_sample(function(x) x > 0), 
                                  A = 0.9*nsamples(pqs.ant.exc))

# For non-Antarctic sponges only:
filter.noa33 <- genefilter_sample(pqs.noant.exc, filterfun_sample(function(x) x > 0), 
                                  A = 0.33*nsamples(pqs.noant.exc))
filter.noa50 <- genefilter_sample(pqs.noant.exc, filterfun_sample(function(x) x > 0), 
                                  A = 0.5*nsamples(pqs.noant.exc))
filter.noa60 <- genefilter_sample(pqs.noant.exc, filterfun_sample(function(x) x > 0), 
                                  A = 0.6*nsamples(pqs.noant.exc))
filter.noa70 <- genefilter_sample(pqs.noant.exc, filterfun_sample(function(x) x > 0), 
                                  A = 0.7*nsamples(pqs.noant.exc))

# For all sponges:
filter.gen33 <- genefilter_sample(pqs.shared, filterfun_sample(function(x) x > 0), 
                                  A = 0.33*nsamples(pqs.shared))
filter.gen50 <- genefilter_sample(pqs.shared, filterfun_sample(function(x) x > 0), 
                                  A = 0.5*nsamples(pqs.shared))
filter.gen60 <- genefilter_sample(pqs.shared, filterfun_sample(function(x) x > 0), 
                                  A = 0.6*nsamples(pqs.shared))
filter.gen70 <- genefilter_sample(pqs.shared, filterfun_sample(function(x) x > 0), 
                                  A = 0.7*nsamples(pqs.shared))

#--- Filter:

# For Antarctic sponges only:
pqs.prev33 <- prune_taxa(filter.pro33, pqs.ant.exc)
pqs.prev50 <- prune_taxa(filter.pro50, pqs.ant.exc)
pqs.prev60 <- prune_taxa(filter.pro60, pqs.ant.exc)
pqs.prev70 <- prune_taxa(filter.pro70, pqs.ant.exc)
pqs.prev80 <- prune_taxa(filter.pro80, pqs.ant.exc)
pqs.prev90 <- prune_taxa(filter.pro90, pqs.ant.exc)

# For Antarctic sponges only:
pqs.prev.noa33 <- prune_taxa(filter.noa33, pqs.noant.exc)
pqs.prev.noa50 <- prune_taxa(filter.noa50, pqs.noant.exc)
pqs.prev.noa60 <- prune_taxa(filter.noa60, pqs.noant.exc)
pqs.prev.noa70 <- prune_taxa(filter.noa70, pqs.noant.exc)
# At 60% and 70% there is a single ASV

# For all sponges:
pqs.prev.gen33 <- prune_taxa(filter.gen33, pqs.shared)
pqs.prev.gen50 <- prune_taxa(filter.gen50, pqs.shared)
pqs.prev.gen60 <- prune_taxa(filter.gen60, pqs.shared)
pqs.prev.gen70 <- prune_taxa(filter.gen70, pqs.shared)

#--- List and count ASVs for further steps:

# For Antarctic sponges only:
taxa33 <- taxa_names(pqs.prev33)
taxa50 <- taxa_names(pqs.prev50)
taxa70 <- taxa_names(pqs.prev70)
taxa80 <- taxa_names(pqs.prev80)
taxa90 <- taxa_names(pqs.prev90)

# For non-Antarctic sponges only:
taxa.noa33 <- taxa_names(pqs.prev.noa33)
taxa.noa50 <- taxa_names(pqs.prev.noa50)
taxa.noa70 <- taxa_names(pqs.prev.noa70)

# For all sponges:
taxa.gen33 <- taxa_names(pqs.prev.gen33)
taxa.gen50 <- taxa_names(pqs.prev.gen50)
taxa.gen60 <- taxa_names(pqs.prev.gen60)
taxa.gen70 <- taxa_names(pqs.prev.gen70)

#--- Euler diagram with the presence of ASVs across habitats:

# For Antarctic sponges only:
list.prev <- list(`33% Antarctic sponges` = taxa33,
                  `50% Antarctic sponges` = taxa50,
                  `70% Antarctic sponges` = taxa70)
# `80% Antarctic sponges` = taxa80, 
# `90% Antarctic sponges` = taxa90)
set.seed(100)                                  
euler.plot.prev <- plot(euler(list.prev), quantities = FALSE,
                        fills = c("#c6deed", "#8ebedb", "#559dc9"),#, "#dda471", "#c49264"), 
                        legend = list(fontsize = 18, side = "right"), 
                        labels = list(labels = c("35", "13", "7"), fontsize = 20,
                                      lineheight = 1),
                        #edges = list(col = c("#fce7d4", "#f9cfa9", "#f6b77e", "#DDA471")),
                        fontsize = 15)
#"#DDA471", "#E9D8DC", "#CCA3AD", "#AA6373" | list(fontsize = 20)

# For non-Antarctic sponges only:                                  
list.prev.noa <- list(`33% non-Antarctic sponges` = taxa.noa33, 
                      `50% non-Antarctic sponges` = taxa.noa50, 
                      `70% non-Antarctic sponges` = taxa.noa70)                                  
set.seed(100)                                  
euler.plot.noa <- plot(euler(list.prev.noa), quantities = FALSE,
                       fills = c("#f2cfc9", "#e69f94", "#da6f5e"), 
                       legend = list(fontsize = 18, side = "right"), 
                       labels = list(labels = c("29", "5", "1"), fontsize = 20, 
                                     lineheight = 1),
                       #edges = list(col = c("#fce7d4", "#f9cfa9", "#f6b77e", "#DDA471")),
                       fontsize = 15)

# For all sponges:                                  
list.prev.gen <- list(`33% sponges` = taxa.gen33, 
                      `50% sponges` = taxa.gen50, 
                      `70% sponges` = taxa.gen70)                                  
set.seed(100)                                  
euler.plot.gen <- plot(euler(list.prev.gen), quantities = FALSE,
                       fills = c("#fce7d4", "#f9cfa9", "#f6b77e"), 
                       legend = list(fontsize = 18, side = "right"), 
                       labels = list(labels = c("34", "11", "1"), fontsize = 20, 
                                     lineheight = 1),
                       #edges = list(col = c("#fce7d4", "#f9cfa9", "#f6b77e", "#DDA471")),
                       fontsize = 15)

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/prev-ant.png", units = "in", width = 6, height = 2.33, res = 300)
euler.plot.prev
dev.off()

png("./results/graphics/prev-noa.png", units = "in", width = 6, height = 2.33, res = 300)
euler.plot.noa
dev.off()

png("./results/graphics/prev-gen.png", units = "in", width = 6, height = 2.33, res = 300)
euler.plot.gen
dev.off()


### Determining the abundance sums of habitat-specifics and generalist ASVs per samples ----

#--- Save columns of abundance for hab-specifics and hab-generalists:

ab.s.a <- as.matrix(sort(sample_sums(pqs.ant.spec.r), decreasing = T))
ab.s.n <- as.matrix(sort(sample_sums(pqs.noa.spec.r), decreasing = T))
ab.g.a <- as.matrix(sort(sample_sums(subset_samples(pqs.gen.r, Habitat == "Antarctic sponge")), 
                         decreasing = T))
ab.g.b <- as.matrix(sort(sample_sums(subset_samples(pqs.gen.r, Habitat == "Non-Antarctic sponge")), 
                         decreasing = T))

# Create dataframes:
df.s <- data.frame(Habitat = rep(c("Antarctic sponge","Non-Antarctic sponge"), 
                                 times = c(89, 77)), 
                   Type = rep("Habitat specific", times = 166),
                   Abundance = c(ab.s.a, ab.s.n),
                   Sample = c(rownames(ab.s.a), rownames(ab.s.n)),
                   species_spec = rep(c("Mycale", "Dendrilla antarctica", "Mycale", 
                                        "Other", "Dendrilla antarctica", "Other", 
                                        "Dendrilla antarctica", "Other", 
                                        "Dendrilla antarctica", "Mycale", "Other", 
                                        "Dendrilla antarctica", "Other"), 
                                      times = c(16, 7, 1, 1, 3, 1, 7, 2, 1, 1, 33, 1, 92)))

df.g <- data.frame(Habitat = rep(c("Antarctic sponge","Non-Antarctic sponge"), 
                                 times = c(89, 68)), 
                   Type = rep("Habitat generalist", times = 157),
                   Abundance = c(ab.g.a, ab.g.b),
                   Sample = c(rownames(ab.g.a), rownames(ab.g.b)),
                   species_spec = rep("Other", times = 1))

df.sg <- rbind(df.s, df.g)

# Export table:
write.table(df.sg,file ="./results/abundance_sum_sg.csv",
            sep="\t",  na = "NA", row.names = F, col.names = T)

# Box plots:
sg.plot <- ggplot(data = df.sg, aes(x = fct_rev(Habitat), y = Abundance)) +
  geom_boxplot(aes(color = Habitat)) +
  coord_flip() +
  facet_wrap(~fct_rev(Type), nrow = 2, scale = "free_y", strip.position = "left") +
  scale_color_manual(values = colors.hab) +
  geom_jitter(aes(shape = species_spec, color = Habitat), alpha = 0.3, width = 0.20, size = 4) +
  #geom_text(check_overlap = TRUE, position=position_jitter(width=0.15))+
  scale_shape_manual(values = c(15, 17, 16), 
                     labels = c("Dendrilla antarctica", 
                                "Mycale (acerata)", 
                                "No evident")) +
  labs(x = NULL, y = "% Abundance", shape = "Species/genus specific pattern") +
  theme_biome_utils() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(legend.position = "left", text = element_text(size = 26), legend.box = "vertical",
        axis.text = element_text(size = 17))
# stat_compare_means(method = "t.test",
#                    label = c("p.signif"), label.x.npc = "center",
#                    label.y.npc = "middle", hide.ns = FALSE)

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/boxplot_sg.png", units = "in", width = 12, height = 16, res = 300)
sg.plot
dev.off()


### Determining the average abundance of habitat specific and generalist ASVs (pie charts) ----

#--- Change ASV IDs for taxa names:

# Antarctic sponge specific ASVs:
taxa_names(pqs.ant.50) <-
  c("Gammapro; Nitrosomonas oligotropha*; ASV398",
    "Bactero; Owenweeksia hongkongensis*; ASV13022",
    "Gammapro; Family Nitrincolaceae; ASV14905",
    "Alphapro; SAR116 clade; ASV15320",
    "Gammapro; Cocleimonas flava*; ASV17537",
    #"Cyanobac*; Dactylothamnos antarcticus; ASV20423",
    "Alphapro; SAR116 clade; ASV28000",
    "Bactero; Polaribacter irgensii; ASV28043",
    "Gammapro; OM43 clade; ASV28188",
    "Bactero; NS9 marine group; ASV31888",
    "Gammapro; Family Nitrincolaceae; ASV33333",
    "Alphapro; Magnetospira; ASV36165",
    "Gammapro; SAR92 clade; ASV37113",
    "Bactero; Cochleicola gelatinilyticus*; ASV42196")

# Non-Antarctic sponge specific ASVs:
taxa_names(pqs.noant50) <-
  c("Gammapro; Endozoicomonas*; ASV6080",
    "Cyanobac; Synechococcus; ASV15216",
    "Cyanobac; Cyanobium gracile*; ASV37530",
    "Alphapro; AEGEAN-169 marine group; ASV39844",
    "Bactero; NS5 marine group; ASV44089")

# Sponge generalist ASVs:
taxa_names(pqs.gen) <-
  c("Alphapro; Candidatus Pelagibacter ubique*; ASV4360",
    "Bactero; Phaeocystidibacter marisrubri*; ASV6045",
    "Gammapro; SAR86 clade; ASV7265",
    "Alphapro; Sulfitobacter undariae*; ASV18593",
    "Gammapro; Family Thioglobaceae; ASV18659",
    "Bactero; Ulvibacter; ASV19117",
    "Bactero; Polaribacter irgensii/staleyi, ASV23626",
    "Alphapro; Yoonia rosea / Loktanella acticola*; ASV29400",
    "Gammapro; OM43 clade; ASV31463",
    "Gammapro; Marimicrobium arenosum*; ASV34765",
    "Alphapro; Planktomarina temperata*; ASV40972")

#--- Divide generalists by habitat:
pqs.gen1 <- subset_samples(pqs.gen, Habitat == "Antarctic sponge")
pqs.gen2 <- subset_samples(pqs.gen, Habitat == "Non-Antarctic sponge")

#--- Merge all sponges:
pqs.avg.ant <- merge_samples(x = pqs.ant.50, group = "Sample_type", fun = "mean")
pqs.avg.noa <- merge_samples(x = pqs.noant50, group = "Sample_type", fun = "mean")
pqs.avg.gen1 <- merge_samples(x = pqs.gen1, group = "Sample_type", fun = "mean")
pqs.avg.gen2 <- merge_samples(x = pqs.gen2, group = "Sample_type", fun = "mean")

#--- Transform to relative abundance:
pqs.avg.ant.rel <- transform_sample_counts(pqs.avg.ant, function(x) x/sum(x))
pqs.avg.noa.rel <- transform_sample_counts(pqs.avg.noa, function(x) x/sum(x))
pqs.avg.gen1 <- transform_sample_counts(pqs.avg.gen1, function(x) x/sum(x))
pqs.avg.gen2 <- transform_sample_counts(pqs.avg.gen2, function(x) x/sum(x))

#--- Create melt table:
pqs.avg.ant.melt <- psmelt(pqs.avg.ant.rel)
pqs.avg.noa.melt <- psmelt(pqs.avg.noa.rel)
pqs.avg.gen1.melt <- psmelt(pqs.avg.gen1)
pqs.avg.gen2.melt <- psmelt(pqs.avg.gen2)

#--- Reorder classes according to abundance:
pqs.avg.ant.melt$OTU <- reorder(pqs.avg.ant.melt$OTU, pqs.avg.ant.melt$Abundance)
pqs.avg.ant.melt$OTU <- factor(pqs.avg.ant.melt$OTU, levels = rev(levels(pqs.avg.ant.melt$OTU)))

pqs.avg.noa.melt$OTU <- reorder(pqs.avg.noa.melt$OTU, pqs.avg.noa.melt$Abundance)
pqs.avg.noa.melt$OTU <- factor(pqs.avg.noa.melt$OTU, levels = rev(levels(pqs.avg.noa.melt$OTU)))

pqs.avg.gen1.melt$OTU <- reorder(pqs.avg.gen1.melt$OTU, pqs.avg.gen1.melt$Abundance)
pqs.avg.gen1.melt$OTU <- factor(pqs.avg.gen1.melt$OTU, levels = rev(levels(pqs.avg.gen1.melt$OTU)))

pqs.avg.gen2.melt$OTU <- reorder(pqs.avg.gen2.melt$OTU, pqs.avg.gen2.melt$Abundance)
pqs.avg.gen2.melt$OTU <- factor(pqs.avg.gen2.melt$OTU, levels = rev(levels(pqs.avg.gen2.melt$OTU)))

#--- Set color pallettes:
colors.ant <- c("#633353", "#d9b7ce", "#84446f", "#84446f",
                   "#d59031", "#eda137", "#f0b35e",  "#a6558b",
                   "#b776a2",  "#f4c687", "#d9b7ce", "#e0c5d7", 
                   "#c999b9")#"#bad1a5"

colors.noa <- c("#8cb369", "#a6558b", "#a3c287", "#f4c687","#e8d3e1")

colors.gen1 <- c("#d9b7ce", "#d59031", "#633353", "#84446f", "#e0c5d7", "#e8d3e1", 
                "#eda137", "#f0b35e", "#efe2eb", "#a6558b", "#b776a2")
colors.gen2 <- c("#a6558b", "#efe2eb", "#b776a2", "#e8d3e1", "#d59031", "#d9b7ce",
                 "#eda137", "#f0b35e", "#e0c5d7", "#84446f", "#633353")

#--- Make pie charts:

# Antarctic sponge specific ASVs:
pie.ant <- ggplot(data = pqs.avg.ant.melt, aes(x = "", y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", width = 1, color = "grey", lwd = 0.1) +
  # labs(title = "Habitat-specific bacteria") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = colors.ant) +
  theme(legend.position = "right", legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.text = element_text(size = 14))

# Non-Antarctic sponge specific ASVs:
pie.noa <- ggplot(data = pqs.avg.noa.melt, aes(x = "", y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", width = 1, color = "grey", lwd = 0.1) +
  # labs(title = "Habitat-specific bacteria") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = colors.noa) +
  theme(legend.position = "right", legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.text = element_text(size = 14))

# Generalist ASVs in Antarctic sponges:
pie.gen1 <- ggplot(data = pqs.avg.gen1.melt, aes(x = "", y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", width = 1, color = "grey", lwd = 0.1) +
  # labs(title = "Habitat-specific bacteria") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = colors.gen1) +
  theme(legend.position = "right", legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.text = element_text(size = 14))

# Generalist ASVs in non-Antarctic sponges:
pie.gen2 <- ggplot(data = pqs.avg.gen2.melt, aes(x = "", y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", width = 1, color = "grey", lwd = 0.1) +
  # labs(title = "Habitat-specific bacteria") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = colors.gen2) +
  theme(legend.position = "right", legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.text = element_text(size = 14))

#--- Plot just the legends:
leg.pie.ant <- cowplot::get_legend(pie.ant)
leg.pie.ant <- as_ggplot(leg.pie.ant)

leg.pie.noa <- cowplot::get_legend(pie.noa)
leg.pie.noa <- as_ggplot(leg.pie.noa)

leg.pie.gen1 <- cowplot::get_legend(pie.gen1)
leg.pie.gen1 <- as_ggplot(leg.pie.gen1)

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/pie-ant.png", units = "in", width = 4, height = 4, res = 300)
pie.ant
dev.off()

png("./results/graphics/pie-noa.png", units = "in", width = 4, height = 4, res = 300)
pie.noa
dev.off()

png("./results/graphics/pie-gen1.png", units = "in", width = 4, height = 4, res = 300)
pie.gen1
dev.off()

png("./results/graphics/pie-gen2.png", units = "in", width = 4, height = 4, res = 300)
pie.gen2
dev.off()

png("./results/graphics/leg-pie-ant.png", units = "in", width = 5.5, height = 4, res = 300)
leg.pie.ant
dev.off()
png("./results/graphics/leg-pie-noa.png", units = "in", width = 5.5, height = 4, res = 300)
leg.pie.noa
dev.off()
png("./results/graphics/leg-pie-gen.png", units = "in", width = 5.5, height = 8, res = 300)
leg.pie.gen1
dev.off()


### Inspecting relative abundances of Antarctic sponge specific ASVs ----

#--- Change name to objects:
pqs.ant.spec <- pqs.ant.50

#--- Concatenate microbiome IDs and sponge species:
sample_names(pqs.ant.spec) <- paste(sample_names(pqs.ant.spec), 
                                    sample_data(pqs.ant.spec)$Sponge_Species, sep = "; ")

#--- Change ASV IDs to customized taxonomy:
taxa_names(pqs.ant.spec) <-
  c("Gammaproteobacteria; Nitrosomonas oligotropha*; ASV398",
    "Bacteroidia; Owenweeksia hongkongensis*; ASV13022",
    "Gammaproteobacteria; Family Nitrincolaceae; ASV14905",
    "Alphaproteobacteria; SAR116 clade; ASV15320",
    "Gammaproteobacteria; Cocleimonas flava*; ASV17537",
    "Cyanobacteria*; Dactylothamnos antarcticus; ASV20423",
    "Alphaproteobacteria; SAR116 clade; ASV28000",
    "Bacteroidia; Polaribacter irgensii; ASV28043",
    "Gammaproteobacteria; OM43 clade; ASV28188",
    "Bacteroidia; NS9 marine group; ASV31888",
    "Gammaproteobacteria; Family Nitrincolaceae; ASV33333",
    "Alphaproteobacteria; Magnetospira; ASV36165",
    "Gammaproteobacteria; SAR92 clade; ASV37113",
    "Bacteroidia; Cochleicola gelatinilyticus*; ASV42196")

#--- Convert to percentage:
pqs.ant.spec.r <- transform_sample_counts(pqs.ant.spec, function(x) x*100)

#--- Sum abundances per sample: 
sample_data(pqs.ant.spec.r)$sample_sums <- sample_sums(pqs.ant.spec.r)

#--- Concatenate table:
pqs.ant.spec.melt <- psmelt(pqs.ant.spec.r)

#--- Sort sample names:
sums.ant.spec <- rownames(as.matrix(sort(sample_sums(pqs.ant.spec.r), 
                                         decreasing = T))) 

#--- Abundance sum of habitat-specific ASVs:
as.matrix(sort(sample_sums(pqs.ant.spec.r), decreasing = T))

#--- Reorder taxa according to abundance:
pqs.ant.spec.melt$OTU <- reorder(pqs.ant.spec.melt$OTU,
                                 pqs.ant.spec.melt$Abundance)
# pqs.ant.spec.melt$OTU <- factor(pqs.ant.spec.melt$OTU,
#                                 levels = rev(levels(pqs.ant.spec.melt$OTU)))
pqs.ant.spec.melt$Sample <- factor(as.factor(pqs.ant.spec.melt$Sample),
                                   levels = sums.ant.spec)

#--- Write table:
write.table(pqs.ant.spec.melt,file ="./results/melt-ant-spec.tsv",
            sep="\t",  na = "NA", row.names = T, col.names = NA)

#--- Barplot of taxonomy composition:
bar.ant.spec.plot <- ggplot(data = pqs.ant.spec.melt, 
                            aes(x = Sample, y = Abundance, fill = OTU, 
                                label = Sample)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  theme_classic() +
  facet_wrap(~Habitat, nrow = 1, scale = "free_x") +
  theme(axis.text.y = element_text(size = 5), axis.ticks.y = element_blank()) +
  coord_flip() + #voltear gráfico
  scale_fill_manual(values = rev(colors.ant)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  # scale_fill_manual(values = colors.1,
  #                   breaks = c("G; Nitrincolaceae",
  #                              "G; SAR92 clade",
  #                              "A; SAR116 clade",
  #                              "B; Polaribacter irgensii",
  #                              "B; Cochleicola gelatinilyticus*",
  #                              "B; Owenweeksia hongkongensis*",
  #                              "B; NS9 marine group")) + #, limits = rownames(as.matrix(phylum_colors))
  ylab("Abundance (%)") +
  labs(x = NULL, fill = "Habitat-specific bacteria") +
  theme(legend.position = "right", legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 17), 
        legend.title = element_text(size = 13, face = "bold"))
# axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/tax-ant-spec.png", units = "in", width = 12, height = 5, res = 300)
bar.ant.spec.plot
dev.off()

### Inspecting relative abundances of non-Antarctic sponges specific ASVs ----

#--- Change name to objects:
pqs.noa.spec <- pqs.noant50

#--- Concatenate microbiome IDs and sponge species:
sample_names(pqs.noa.spec) <- paste(sample_names(pqs.noa.spec), 
                                    sample_data(pqs.noa.spec)$Sample_Species, sep = "; ")

#--- Change ASV IDs to customized taxonomy:
taxa_names(pqs.noa.spec) <-
  c("Gammaproteobacteria; Endozoicomonas*; ASV6080",
    "Cyanobacteriia; Synechococcus strain CC9902; ASV15216",
    "Cyanobacteriia; Cyanobium gracile*; ASV37530",
    "Alphaproteobacteria; AEGEAN-169 marine group; ASV39844",
    "Bacteroidia; NS5 marine group; ASV44089")

#--- Convert to percentage:
pqs.noa.spec.r <- transform_sample_counts(pqs.noa.spec, function(x) x*100)

#--- Concatenate table:
pqs.noa.spec.melt <- psmelt(pqs.noa.spec.r)

#--- Sum abundances per sample and sort sample names: 
sums.noa.spec <- rownames(as.matrix(sort(sample_sums(pqs.noa.spec.r), 
                                         decreasing = T)))

#--- Abundance sum of habitat-specific ASVs:
sort(sample_sums(pqs.noa.spec.r), decreasing = T)

#--- Reorder and change levels:
pqs.noa.spec.melt$OTU <- reorder(pqs.noa.spec.melt$OTU,
                                 pqs.noa.spec.melt$Abundance)
# pqs.noa.spec.melt$OTU <- factor(pqs.noa.spec.melt$OTU,
#                                 levels = rev(levels(pqs.noa.spec.melt$OTU)))
pqs.noa.spec.melt$Sample <- factor(as.factor(pqs.noa.spec.melt$Sample),
                                   levels = sums.noa.spec)

#--- Barplot of taxonomy composition:
bar.noa.spec.plot <- ggplot(data = pqs.noa.spec.melt, 
                            aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  theme_classic() +
  facet_wrap(~Habitat, nrow = 1, scale = "free_x") +
  theme(axis.text.y = element_text(size = 5), axis.ticks.y = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip() + #voltear gráfico
  scale_fill_manual(values = rev(colors.noa)) +
  ylab("Abundance (%)") +
  labs(x = NULL, fill = "Habitat-specific bacteria") +
  theme(legend.position = "right", legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 17),
        legend.title = element_text(size = 13, face = "bold"))

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/tax-noa-spec.png", units = "in", width = 12, height = 5, res = 300)
bar.noa.spec.plot
dev.off()

### Inspecting relative abundances of sponge generalist ASVs ----

#--- Concatenate microbiome IDs and sponge species:
sample_names(pqs.gen) <- paste(sample_names(pqs.gen), 
                               sample_data(pqs.gen)$Sample_Species, sep = "; ")

#--- Change ASV IDs to customized taxonomy:
taxa_names(pqs.gen) <-
  c("Alphaproteobacteria; Candidatus Pelagibacter ubique*; ASV4360",
    "Bacteroidia; Phaeocystidibacter marisrubri*; ASV6045",
    "Gammaproteobacteria; SAR86 clade; ASV7265",
    "Alphaproteobacteria; Sulfitobacter undariae*; ASV18593",
    "Gammaproteobacteria; Family Thioglobaceae; ASV18659",
    "Bacteroidia; Ulvibacter; ASV19117",
    "Bacteroidia; Polaribacter irgensii/staleyi; ASV23626",
    "Alphaproteobacteria; Yoonia rosea / Loktanella acticola*; ASV29400",
    "Gammaproteobacteria; OM43 clade; ASV31463",
    "Gammaproteobacteria; Marimicrobium arenosum*; ASV34765",
    "Alphaproteobacteria; Planktomarina temperata*; ASV40972")

#--- Convert to percentage:
pqs.gen.r <- transform_sample_counts(pqs.gen, function(x) x*100)

#--- Concatenate table:
pqs.gen.melt <- psmelt(pqs.gen.r)

#--- Sum abundances per sample and sort sample names: 
sums.gen <- rownames(as.matrix(sort(sample_sums(pqs.gen.r), 
                                    decreasing = T)))

#--- Abundance sum of habitat-specific ASVs:
sort(sample_sums(subset_samples(pqs.gen.r, Habitat == "Antarctic sponge")), 
     decreasing = T)
sort(sample_sums(subset_samples(pqs.gen.r, Habitat == "Non-Antarctic sponge")), 
     decreasing = T)

#--- Reorder taxa according to abundance:
pqs.gen.melt$OTU <- reorder(pqs.gen.melt$OTU, pqs.gen.melt$Abundance)
# pqs.ant.spec.melt$Fulltax <- factor(pqs.ant.spec.melt$Fulltax, 
#                                 levels = rev(levels(pqs.ant.spec.melt$Fulltax)))
pqs.gen.melt$Sample <- factor(as.factor(pqs.gen.melt$Sample),
                              levels = sums.gen)

#--- Barplot of taxonomy composition:
bar.gen.plot <- ggplot(data = pqs.gen.melt, aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  theme_classic() +
  facet_wrap(~Habitat, nrow = 1, scale = "free_y") +
  theme(axis.text.y = element_text(size = 5), axis.ticks.y = element_blank()) + 
  coord_flip() + #voltear gráfico
  #scale_fill_manual(values = colors.1, 
  scale_fill_manual(values = colors.gen,
                    breaks = c("Gammaproteobacteria; OM43 clade; ASV31463",
                               "Gammaproteobacteria; Marimicrobium arenosum*; ASV34765",
                               "Gammaproteobacteria; Family Thioglobaceae; ASV18659",
                               "Gammaproteobacteria; SAR86 clade; ASV7265",
                               "Alphaproteobacteria; Sulfitobacter undariae*; ASV18593",
                               "Alphaproteobacteria; Candidatus Pelagibacter ubique*; ASV4360",
                               "Alphaproteobacteria; Planktomarina temperata*; ASV40972",
                               "Alphaproteobacteria; Yoonia rosea / Loktanella acticola*; ASV29400",
                               "Bacteroidia; Polaribacter irgensii/staleyi; ASV23626",
                               "Bacteroidia; Phaeocystidibacter marisrubri*; ASV6045",
                               "Bacteroidia; Ulvibacter; ASV19117")) + #, limits = rownames(as.matrix(phylum_colors))
  ylab("Abundance (%)") +
  labs(x = NULL, fill = "Habitat-generalist bacteria") +
  theme(legend.position = "right", legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 17), 
        legend.title = element_text(size = 13, face = "bold"))

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/tax-gen.png", units = "in", width = 12, height = 5, res = 300)
bar.gen.plot
dev.off()


### Determining exclusive and shared taxonomic groups ----

#--- Separate habitats:
taxa.ant <- get_taxa_unique(pqs.ant, "Species") # Change for each taxonomic level
taxa.noant <- get_taxa_unique(pqs.noant, "Species") # Change for each taxonomic level

#--- Join taxa names in a list:
taxa <- list("Antarctic" = taxa.ant, "Non-Antarctic" = taxa.noant)

#--- Plot venn of taxa:
ggvenn(taxa)
intersect(taxa.ant, taxa.noant)
setdiff(taxa.ant, taxa.noant)


### Inspecting the prevalence of Antarctic exclusive ASVs (bonus) ----

#--- Evaluate prevalence of taxa across samples:
preval.ant.exc <- data.frame(Prevalence = colSums(!!otu_table(pqs.ant.exc)))
data.frame(table(nSamples = colSums(!!otu_table(pqs.ant.exc))))

#--- Plot prevalence:
preval.spec.plot <- ggplot(data = preval.ant.exc, aes(x = Prevalence)) +
  geom_histogram(binwidth = 0.5) + #  0.0001
  labs(x = "# Antarctic sponge samples", y = "# ASVs") +
  scale_x_continuous(breaks = seq(0, 90, 15)) +
  scale_y_continuous(breaks = seq(0, 7000, 1000)) +
  #scale_y_sqrt(breaks = seq(0, 7000, 25)) +
  theme_bw()

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/prev-distrib.png", units = "in", width = 8, height = 5, res = 300)
preval.spec.plot
dev.off()


#--- Save things -----------------------------------------------

# Save objects:
saveRDS(pqs.ant, file = "./rds-files/pqs_ant.rds")
saveRDS(pqs.noant, file = "./rds-files/pqs_noa.rds")
saveRDS(pqs.ant.50, file = "./rds-files/pqs_ant_spec.rds")
saveRDS(pqs.noant50, file = "./rds-files/pqs_noa_spec.rds")
saveRDS(pqs.gen, file = "./rds-files/pqs_gen.rds")
saveRDS(gen.asvs, file = "./rds-files/gen-asvs.rds")
saveRDS(spec.asvs, file = "./rds-files/spec-asvs.rds")
saveRDS(spec.noan.asvs, file = "./rds-files/spec-noan-asvs.rds")

# Save the current work and objects:
save.image("R-graph-venn.RData")

