# ANALYZING ALPHA DIVERSITY #
# ------------------------- #

#--- Load objects needed for these steps ----

pqs <- readRDS(file = "./rds-files/pqs.rds")


#--- Define colors for graphics ----
colors.hab <- rev(brewer.pal(n = 11, name = "RdBu"))[c(3,9)]


### Calculating sequence abundance per sample ----

#--- Calculate abundances:
n.seqs.sample <- data.frame(cbind(sample_sums(pqs), get_variable(pqs, "Sponge_Species")))
colnames(n.seqs.sample)[1:2] <- c("Abundance", "Sample")
n.seqs.sample$Abundance <- as.numeric(n.seqs.sample$Abundance)

#--- Write table of total abundances:
write.table(n.seqs.sample, file = "./results/total-abundance.tsv", 
            sep = "\t", na = "NA", row.names = T, col.names = FALSE)


### Estimating the microbial diversity ----

#--- Calculate indexes:
adiv <- estimate_richness(pqs, split = TRUE, 
                              measures = c("Observed", "Simpson", "Shannon"))

#--- Total ASVs (all samples):
sum(adiv[,'Observed'])

#--- Write table:
write.table(adiv, file = "./results/alpha-diversity.tsv", 
            sep = "\t", dec = ".",na = "NA", row.names = T, col.names = NA)

#--- Calculate Faith phylogenetic diversity:
faith <- pd(otu_table(pqs), phy_tree(pqs), include.root = T)

#--- Join alpha diversity data to metadata:
adiv.mdata <- data.table(metadata, faith, adiv, keep.rownames = FALSE)


### Performing LMM and GLMM ----

#--- For observed richness:

# Using the factor "Environment":
a <- glmer(SR ~ Environment +
             # (1|Sponge_Order) +
             (1|Sponge_Family) +
             (1|Sequencing_run_ID), data = adiv.mdata, family = "poisson")
tab_model(a)
summary(a)
plot_model(a, type = "eff", terms = "Environment", show.data = T)

# Using the factor "Habitat":
b <- glmer(SR ~ Habitat +
             # (1|Sponge_Order) +
             (1|Sponge_Family) +
             (1|Sequencing_run_ID), data = adiv.mdata, family = "poisson")
tab_model(b)
summary(b)
plot_model(b, type = "eff", terms = "Habitat", show.data = F)

#--- For Shannon index:

# Using the factor "Environment":
c <- lmer(Shannon ~ Environment + 
            # (1|Sponge_Order) +
            (1|Sponge_Family) +
            (1|Sequencing_run_ID), data = adiv.mdata)
tab_model(c)
summary(c)
plot_model(c, type = "eff", terms = "Environment", show.data = T)

# Using the factor "Habitat":
d <- lmer(Shannon ~ Habitat + 
            # (1|Sponge_Order) +
            (1|Sponge_Family) +
            (1|Sequencing_run_ID), data = adiv.mdata)
tab_model(d)
summary(d)
plot_model(d, type = "eff", terms = "Habitat", show.data = F)

adiv.mdata$jost <- exp(adiv.mdata$Shannon) # Alternative

d1 <- lmer(jost ~ Habitat + 
            (1|Sponge_Family) +
            (1|Sequencing_run_ID), data = adiv.mdata)
tab_model(d1)
plot_model(d1, type = "eff", terms = "Habitat", show.data = F)

#--- For Simpson index:

# Using the factor "Environment":
e <- lmer(Simpson ~ Environment + 
            # (1|Sponge_Order) +
            (1|Sponge_Family) +
            (1|Sequencing_run_ID), data = adiv.mdata)
tab_model(e)
plot_model(e, type = "eff", terms = "Environment", show.data = T)

f <- lmer(Simpson ~ Habitat +
            # (1|Sponge_Order) +
            (1|Sponge_Family) +
            (1|Sequencing_run_ID), data = adiv.mdata)
tab_model(f)
summary(f)
plot_model(f, type = "eff", terms = "Habitat", show.data = T)

#--- For Faith PD:

pd <- lmer(PD ~ Habitat + 
             (1|Sponge_Family) +
             (1|Sequencing_run_ID), data = adiv.mdata)
tab_model(pd)
plot_model(pd, type = "eff", terms = "Habitat", show.data = F)


### Comparing models ----

#--- Compare models with Akaike:
aictab(cand.set = list(a, b), modnames = c("a", "b"))
aictab(cand.set = list(c, d), modnames = c("c", "d"))
aictab(cand.set = list(e, f), modnames = c("e", "f"))


### Plotting microbial diversity ----

#--- Bring results from LMM and GLMM:
rich.result <- ggpredict(model = b, terms = "Habitat")
shan.result <- ggpredict(model = d, terms = "Habitat")
simp.result <- ggpredict(model = f, terms = "Habitat")
fait.result <- ggpredict(model = pd, terms = "Habitat")

#--- Previous plots:
plot(rich.result)
plot(shan.result)
plot(simp.result)
plot(fait.result)

#--- Customized bar plots:

# For observed richness:
rich.bar <- ggplot(rich.result, aes(x = x, y = predicted, fill = x)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge()) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high, y = predicted),
                width = 0.4) +
  geom_expected_point(aes(x = x, y = predicted), show.legend = FALSE) +
  theme_classic() +
  labs(x = NULL, y = "Observed richness", fill = NULL) +
  scale_fill_manual(values = colors.hab) +
  theme(legend.position = "none", axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.y = element_text(size = 15),
        text = element_text(size = 18))

# For Shannon index:
shan.bar <- ggplot(shan.result, aes(x = x, y = predicted, fill = x)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge()) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high, y = predicted),
                width = 0.4) +
  geom_expected_point(aes(x = x, y = predicted), show.legend = FALSE) +
  theme_classic() +
  labs(x = NULL, y = "Shannon", fill = NULL) +
  scale_fill_manual(values = colors.hab) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), legend.text = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        text = element_text(size = 18))

# For Simpson index:
simp.bar <- ggplot(simp.result, aes(x = x, y = predicted, fill = x)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge()) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high, y = predicted),
                width=0.4) +
  geom_expected_point(aes(x = x, y = predicted), show.legend = FALSE) +
  theme_classic() +
  labs(x = NULL, y = "Simpson", fill = NULL) +
  scale_fill_manual(values = colors.hab)+
  theme(legend.position = "none", axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.y = element_text(size = 15),
        text = element_text(size = 18))

fait.bar <- ggplot(fait.result, aes(x = x, y = predicted, fill = x)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge()) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high, y = predicted),
                width=0.4) +
  geom_expected_point(aes(x = x, y = predicted), show.legend = FALSE) +
  theme_classic() +
  labs(x = NULL, y = "Faith", fill = NULL) +
  scale_fill_manual(values = colors.hab)+
  theme(legend.position = "none", axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.y = element_text(size = 15),
        text = element_text(size = 18))

#--- Save HQ image (also can use ggsave() or png() ):
png("./alpha-rich.png", units = "in", width = 3, height = 3, res = 300)
rich.bar
dev.off()

png("./alpha-simp.png", units = "in", width = 3, height = 3, res = 300)
simp.bar
dev.off()

png("./alpha-shan.png", units = "in", width = 3, height = 3, res = 300)
shan.bar
dev.off()

#--- Plot just the legend:
legend <- cowplot::get_legend(shan.bar)
legend <- as_ggplot(legend)

png("./alpha-legend.png", units = "in", width = 9, height = 0.5, res = 300)
legend #width = 2, height = 3
dev.off()

#--- Scatter plots with error bars:

# rich.plot <- ggplot() +
#   geom_point(data = adiv.mdata,
#              aes(x = Habitat, y = SR, color = Habitat),
#              position = position_dodge2(0.7)) +
#   geom_errorbar(data = rich.result,
#                 aes(x = x, ymin = conf.low, ymax = conf.high, y = predicted)) +
#   geom_expected_point(data = c.result,
#                       aes(x = x, y = predicted)) +
#   labs(x = NULL, y = "Species richness", color = NULL) +
#   theme_classic() +
#   theme(legend.position = "right",
#         axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   scale_color_manual(values = colors.hab)
# 
# shan.plot <- ggplot() +
#   geom_point(data = adiv.mdata,
#              aes(x = Habitat, y = Shannon, color = Habitat),
#              position = position_dodge2(0.7)) +
#   geom_errorbar(data = shan.result,
#                 aes(x = x, ymin = conf.low, ymax = conf.high, y = predicted)) +
#   geom_expected_point(data = d.result,
#                       aes(x = x, y = predicted)) +
#   labs(x = NULL, y = "Shannon", color = NULL) +
#   theme_classic() +
#   theme(legend.position = "right",
#         axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   scale_color_manual(values = colors.hab)
# 
# simp.plot <- ggplot() +
#   geom_point(data = adiv.mdata,
#              aes(x = Habitat, y = Simpson, color = Habitat),
#              position = position_dodge2(0.7)) +
#   geom_errorbar(data = simp.result,
#                 aes(x = x, ymin = conf.low, ymax = conf.high, y = predicted)) +
#   geom_expected_point(data = k.result,
#                       aes(x = x, y = predicted)) +
#   labs(x = NULL, y = "Simpson", color = NULL) +
#   theme_classic() +
#   theme(legend.position = "right",
#         axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   scale_color_manual(values = colors.hab)


### Pattern of biogeography over richness (bonus) ----

#--- Plot richness across latitudes:
rich.biog <- ggplot(data = adiv.mdata, aes(x = Observed, y = Latitude)) +
  geom_point() +
  scale_y_continuous(breaks = seq(-80, 80, 20)) +
  geom_smooth(col = "red", size = 1, orientation = "y")

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/rich-biogeo.png", units = "in", width = 7, height = 5, res = 300)
rich.biog
dev.off()


### Save things ------------------------------------------------------------

#--- Save objects:
saveRDS(adiv.mdata, file = "./rds-files/adiv.mdata.rds")
saveRDS(n.seqs.sample, file = "./rds-files/n_seqs_sample.rds")
saveRDS(adiv, file = "./rds-files/alpha_diversity.rds")

#---- Save the current work and objects:
save.image("R-graph-adiv.RData")

