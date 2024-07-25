# STRUCTURING THE METADATA #
# ------------------------ #

#--- Load sample metadata file ----

metadata.prev <- read.table(file = "metadata-SMP-filt.tsv", header = TRUE, sep = "\t")
metadata <- read.table(file = "metadata-sponges.tsv", header = TRUE, sep = "\t")


#--- Define colors for graphics ----

colors.hab <- rev(brewer.pal(n = 11, name = "RdBu"))[c(3,9)]
colors.hab2 <- c(	"#4393C3", "#c37343")
color.latzone <- c(Antarctic = "#2166ac", Temperate = "#DD8B71", Tropical = "#b2182b")#d88b95
colors.1 <- c("darkorchid4", "deepskyblue2", "coral3", "tomato4", 
              "turquoise", "darkgoldenrod3", "lightpink", "#2166ac", "seagreen4",
              "limegreen", "orchid1",  "olivedrab4", "red1", "skyblue3", "maroon4",
              "darkorange", "gold2", "khaki", "plum4", "maroon2", "lightgreen",
              "aquamarine", "grey0", "mediumpurple", "sienna4", "tan", "darkcyan",
              "yellow", "thistle", "gray81")#darkblue


### Adding MEOW information to the metadata ----

#--- Set MEOW data (Spalding et al., 2007, https://doi.org/10.1641/B570707):
data(regions)
data(provinces)
data(realms)
realms.df
provinces.df
regions.df

#--- Add MEOW columns to the metadata:
regionalData <- getRegions(metadata.prev$Latitude, metadata.prev$Longitude)
metadata.prev <- cbind(metadata.prev, regionalData)

#--- Save ecoregions dataframe:
write_csv(regions.df, file = "world-regions.csv", col_names = T)


### Plotting a world map of sampling sites ----

#--- Load world map data:
world <- map_data("world") %>% filter(!long > 180)

#--- Change factor levels:
metadata$Environment <- factor(metadata$Environment, levels = c("Tropical", "Temperate", "Antarctic"))

#--- Plot world map:

# By environment:
sites.plot <- ggplot() +
  geom_map(data = world, map = world, aes(map_id = region),
           color = "white", fill = "lightgray", linewidth = 0.1) +
  expand_limits(x = world$long, y = world$lat) +
  geom_point(data = metadata, aes(x = Longitude, y = Latitude, color = Environment),
             size = 2, show.legend = TRUE) +
  labs(x = NULL, y = NULL, color = NULL) +
  theme_light() +
  scale_x_continuous(breaks = seq(-150, 150, 50)) +
  scale_y_continuous(breaks = seq(-80, 80, 20)) +
  scale_color_manual(values = rev(color.latzone)) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 13))
#coord_fixed()
#alternatively to geom_map and expand_limits: 
#geom_polygon(data = world, aes(x = long, y = lat, group = group), color = "white", fill = "lightgray", size = 0.1) +

# By ecoregion:
sites.plot2 <- ggplot() +
  geom_map(data = world, map = world, aes(map_id = region),
           color = "white", fill = "lightgray", size = 0.1) +
  expand_limits(x = world$long, y = world$lat) +
  geom_point(data = metadata, aes(x = Longitude, y = Latitude, color = Marine_Ecoregion,
                                    shape = Habitat),
             size = 5, show.legend = TRUE) +
  labs(x = NULL, y = NULL, color = "Marine_Ecoregion", shape = "Environment") +
  theme_light() +
  scale_x_continuous(breaks = seq(-150, 150, 50)) +
  scale_y_continuous(breaks = seq(-80, 80, 20)) +
  scale_color_manual(values = colors.1) +
  guides(color = guide_legend(override.aes = list(size = 5)))

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/map-sites.png", units = "in", width = 7, height = 3.5, res = 300)
sites.plot
dev.off()

png("./results/graphics/map-sites2.png", units = "in", width = 6, height = 3, res = 300)
sites.plot2
dev.off()


### Plotting the number of sponge samples ----
 
#--- Make table with sample sizes:
envs <- c("Antarctic", "Temperate", "Tropical")
n.samples <- c(89, 50, 40)
samples.df <- data.frame(envs, n.samples)

#--- Set levels:
samples.df$envs <- factor(samples.df$envs, levels = envs)

#--- Plot pie chart:
pie.nsamples <- ggplot(data = samples.df, aes(x = "", y = n.samples, fill = envs)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0, direction = -1) +
  theme_void() +
  # labs(title = "# samples per environment") +
  scale_fill_manual(values = color.latzone) +
  theme(legend.position = "none", legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = paste0(n.samples)), color = "white",
            position = position_stack(vjust = 0.5), size = 10)


#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/nsamples.png", units = "in", width = 3, height = 3, res = 300)
pie.nsamples
dev.off()


### Plotting distribution of ecoregions and sponge families (alluvial plot)----

#--- Prepare metadata:
metadata$Freq <- rep(x = 1, 179)

#--- Plot:
alluv.plot <- ggplot(data = metadata, aes(axis1 = Environment, axis2 = Marine_Ecoregion, 
                            axis3 = Sponge_Family,
                            y = Freq)) +
  geom_alluvium(aes(fill = Habitat), 
                curve_type = "cubic", width = 1/2, linewidth = 0.3) +
  geom_stratum(alpha = .3, width = 1/2, linewidth = 0.5, 
               color = "white") +
  scale_x_discrete(limits = c("Environment", "Sponge family"), 
                   expand = c(.2, .05)) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3.5) +
  theme_void() +
  scale_fill_manual(values = colors.hab2) +
  theme(text = element_text(size = 14), legend.position = "bottom")

#--- Save HQ image (also can use ggsave() or png() ):
png("./results/graphics/alluvial.png", units = "in", width = 12, height = 7, res = 300)
alluv.plot
dev.off()


#--- Save things ----

#--- Save the current work and objects:
save.image("R-graph-mdata.RData")
