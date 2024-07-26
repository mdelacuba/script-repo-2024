# GETTING READY #
# ------------- #

## Script table to execute all other R scripts.

# Loading packages -----------------------------------------------------------

library("dada2")
library("phyloseq")
library("Biostrings")
library("ggplot2")
library("dplyr")
library("tidyr")
library("tibble")
library("readxl")
library("readr")
library("stringr")
library("devtools")
library("kableExtra")
library("vegan")
library("ape")
library("tidyverse")
library("RColorBrewer")
library(wesanderson)
library(SOfun)
library(scales)
library("data.table")
library(psadd)
library(wrapr)
library(grid)
library(MicEco)
library(eulerr)
library(ggpubr)
library(ggtext)
library(png)
library(pairwiseAdonis)
library(RVAideMemoire)
library("picante")
library("coin")
library("factoextra")
library("dendextend")
library("BiodiversityR")
library("maps")
library("clustsig")
library("ggords")
library("pvclust")
library("meowR")
library(jcolors)
library(microbiome)
library(microbiomeutilities)
library(ggVennDiagram)
library(VennDiagram)
library("ggvenn")
library(pheatmap)
library(ggplotify)
library(patchwork)
library(corrplot)
library(RColorBrewer)
library(UpSetR)
library(ranacapa)
library(fossil)
library(ggalluvial)
library(microbiomeMarker)
library(grid)
library(plyr)
library(treemapify)
library("labdsv")
library(lme4)
library(sjPlot)
library(flexplot)
library(AICcmodavg)
# library(ggeffects)
# library(ggggeffects)
# library(MicrobiotaProcess)
# detach(package:MicrobiotaProcess, unload = TRUE)
# library(ComplexHeatmap) #masking pheatmap, activate if use it

# Creating general directories ----------------------------------------------

results_dir <- "./results/"
database_dir <- "./databases/"
rds_dir <- "./rds-files"
dir.create("results")
dir.create("databases")
dir.create("rds-files")

# Running scripts -------------------------------------------------------------

source("script2_plot-exp-design.R")
rm(list = ls())
source("script3_dada2-pre-proc.r")
rm(list = ls())
source("script4_dada2-filter-n-trim.r")
rm(list = ls())
source("script5_dada2-denoise-asv.r")
rm(list = ls())
source("script6_dada2-rem-chim.r")
rm(list = ls())
source("script7_dada2-taxonomy.r")
rm(list = ls())
source("script8_phyloseq-pre-proc.R")
rm(list = ls())
source("script9_alpha-diver.R")
rm(list = ls())
source("script10_phyla-compos.R")
rm(list = ls())
source("script11_beta-diver.R")
rm(list = ls())
source("script12_hab-spec-gen.R")
rm(list = ls())
source("script13_pre-network.R")
rm(list = ls())
# Note: "script1_exp-design.sh" and "script14-extra-manage.sh" must be run in bash.


# Loading environments -------------------------------------------------------

# load("R-graph-mdata.RData")
# load("R-quality-1821.RData")
# load("R-quality-2013.RData")
# load("R-quality-2013cgb.RData")
# load("R-quality-2014.RData")
# load("R-quality-2015.RData")
# load("R-quality-2020ss.RData")
# load("R-quality-h2000l7.RData")
# load("R-quality-h2000l8.RData")
# load("R-quality-h2500l3.RData")
# load("R-quality-h2500l4.RData")
# load("R-quality-h2500l5.RData")
# load("R-quality-h2500l7.RData")
# load("R-filtntrim.RData")
# load("R-denoise-1821.RData")
# load("R-denoise-2013.RData")
# load("R-denoise-2013cgb.RData")
# load("R-denoise-2014.RData")
# load("R-denoise-2015.RData")
# load("R-denoise-2020ss.RData")
# load("R-denoise-h2000l7.RData")
# load("R-denoise-h2000l8.RData")
# load("R-denoise-h2500l3.RData")
# load("R-denoise-h2500l4.RData")
# load("R-denoise-h2500l5.RData")
# load("R-denoise-h2500l7.RData")
# load("R-rem-chimeras.RData")
# load("R-taxonomy.RData")
# load("R-phyloseq.RData")
# load("R-graph-adiv.RData")
# load("R-graph-taxc.RData")
# load("R-graph-bdiv.RData")
# load("R-graph-venn.RData")
# load("R-graph-netw.RData")


#--- Save data for downstream steps ----

saveRDS(database_dir, file = "./rds-files/database_dir.rds")
saveRDS(results_dir, file = "./rds-files/results_dir.rds")
