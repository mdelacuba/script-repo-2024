# PRE-PROCESSING: DEFINING PATHS AND INSPECTING QUALITY PROFILES OF READS #
# ----------------------------------------------------------------------- #

## Per sequencing run basis

### For sequencing run ID "1821" ----

#--- Save job directories:
samples_dir.1821 <- "./sponges-1821/"
quality_dir.1821 <- "./sponges-1821/quality-1821/"  # qual pdf
filtered_dir.1821 <- "./sponges-1821/filtered-1821/"  # fastq filtered

#--- Create the respective folders:
dir.create(path.expand("sponges-1821/quality-1821"))
dir.create(path.expand("sponges-1821/filtered-1821"))

#--- List all fastq single-end files (R) and extract sample names:
fns.1821 <- sort(list.files(samples_dir.1821, full.names = TRUE))
fns.1821 <- fns.1821[str_detect(basename(fns.1821), ".fastq")]
sample.names.1821 <- str_split(basename(fns.1821), pattern = "_", simplify = TRUE)
sample.names.1821 <- sample.names.1821[, 1]

#--- Calcule number of single-end sequences:
df.1821 <- data.frame()
for (i in 1:length(fns.1821)) {
  # use the dada2 function fastq.geometry
  geom.1821 <- fastq.geometry(fns.1821[i])
  # extract the information on number of sequences and file name
  df_one_row.1821 <- data.frame(n_seq = geom.1821[1], file_name = basename(fns.1821[i]))
  # add one line to data frame
  df.1821 <- bind_rows(df.1821, df_one_row.1821)
}

#--- Plot histogram of num seqs:
ggplot(df.1821, aes(x = n_seq)) + geom_histogram(alpha = 0.5, position = "identity", binwidth = 200) + xlim(0, 85000)

#--- Plot quality visualization files:
for (i in 1:length(fns.1821)) {
  # Use dada2 function to plot quality
  qplots.1821 <- plotQualityProfile(fns.1821[i])
  # Save as a pdf
  qplots_png.1821 <- paste0(quality_dir.1821, basename(fns.1821[i]), "-qual.png")
  ggsave(plot = qplots.1821, filename = qplots_png.1821, device = "png", width = 15, height = 15, 
         scale = 1, units = "cm")
}

#--- Save data for downstream steps:
saveRDS(sample.names.1821, file ="./rds-files/sample_names_1821.rds")
saveRDS(fns.1821, file ="./rds-files/fns_1821.rds")
saveRDS(df.1821, file ="./rds-files/df_1821.rds")

#--- Save the current work and objects:
save.image("R-quality-1821.RData")
rm(list = ls())


### For sequencing run ID "2013" ----

#--- Save job directories:
samples_dir.2013 <- "./sponges-2013/"
quality_dir.2013 <- "./sponges-2013/quality-2013/"  # qual pdf
filtered_dir.2013 <- "./sponges-2013/filtered-2013/"  # fastq filtered

#--- Create the respective folders:
dir.create(path.expand("sponges-2013/quality-2013"))
dir.create(path.expand("sponges-2013/filtered-2013"))

#--- List all fastq single-end files (R) and extract sample names:
fns.2013 <- sort(list.files(samples_dir.2013, full.names = TRUE))
fns.2013 <- fns.2013[str_detect(basename(fns.2013), ".fastq")]
sample.names.2013 <- str_split(basename(fns.2013), pattern = "_", simplify = TRUE)
sample.names.2013 <- sample.names.2013[, 1]

#--- Calcule number of single-end sequences:
df.2013 <- data.frame()
for (i in 1:length(fns.2013)) {
  # use the dada2 function fastq.geometry
  geom.2013 <- fastq.geometry(fns.2013[i])
  # extract the information on number of sequences and file name
  df_one_row.2013 <- data.frame(n_seq = geom.2013[1], file_name = basename(fns.2013[i]))
  # add one line to data frame
  df.2013 <- bind_rows(df.2013, df_one_row.2013)
}

#--- Plot histogram of num seqs:
ggplot(df.2013, aes(x = n_seq)) + geom_histogram(alpha = 0.5, position = "identity", binwidth = 200) + xlim(0, 85000)

#--- Plot quality visualization files:
for (i in 1:length(fns.2013)) {
  # Use dada2 function to plot quality
  qplots.2013 <- plotQualityProfile(fns.2013[i])
  # Save as a pdf
  qplots_png.2013 <- paste0(quality_dir.2013, basename(fns.2013[i]), "-qual.png")
  ggsave(plot = qplots.2013, filename = qplots_png.2013, device = "png", width = 15, height = 15, 
         scale = 1, units = "cm")
}

#--- Save data for downstream steps:
saveRDS(sample.names.2013, file ="./rds-files/sample_names_2013.rds")
saveRDS(fns.2013, file ="./rds-files/fns_2013.rds")
saveRDS(df.2013, file ="./rds-files/df_2013.rds")

#--- Save the current work and objects:
save.image("R-quality-2013.RData")
rm(list = ls())



### For sequencing run ID "2013cgb" ----

#--- Save job directories:
samples_dir.2013cgb <- "./sponges-2013cgb/"
quality_dir.2013cgb <- "./sponges-2013cgb/quality-2013cgb/"  # qual pdf
filtered_dir.2013cgb <- "./sponges-2013cgb/filtered-2013cgb/"  # fastq filtered

#--- Create the respective folders:
dir.create(path.expand("sponges-2013cgb/quality-2013cgb"))
dir.create(path.expand("sponges-2013cgb/filtered-2013cgb"))

#--- List all fastq single-end files (R) and extract sample names:
fns.2013cgb <- sort(list.files(samples_dir.2013cgb, full.names = TRUE))
fns.2013cgb <- fns.2013cgb[str_detect(basename(fns.2013cgb), ".fastq")]
sample.names.2013cgb <- str_split(basename(fns.2013cgb), pattern = "_", simplify = TRUE)
sample.names.2013cgb <- sample.names.2013cgb[, 1]

#--- Calcule number of single-end sequences:
df.2013cgb <- data.frame()
for (i in 1:length(fns.2013cgb)) {
  # use the dada2 function fastq.geometry
  geom.2013cgb <- fastq.geometry(fns.2013cgb[i])
  # extract the information on number of sequences and file name
  df_one_row.2013cgb <- data.frame(n_seq = geom.2013cgb[1], file_name = basename(fns.2013cgb[i]))
  # add one line to data frame
  df.2013cgb <- bind_rows(df.2013cgb, df_one_row.2013cgb)
}

#--- Plot histogram of num seqs:
ggplot(df.2013cgb, aes(x = n_seq)) + geom_histogram(alpha = 0.5, position = "identity", binwidth = 200) + xlim(0, 85000)

#--- Plot quality visualization files:
for (i in 1:length(fns.2013cgb)) {
  # Use dada2 function to plot quality
  qplots.2013cgb <- plotQualityProfile(fns.2013cgb[i])
  # Save as a pdf
  qplots_png.2013cgb <- paste0(quality_dir.2013cgb, basename(fns.2013cgb[i]), "-qual.png")
  ggsave(plot = qplots.2013cgb, filename = qplots_png.2013cgb, device = "png", width = 15, height = 15, 
         scale = 1, units = "cm")
}

#--- Save data for downstream steps:
saveRDS(sample.names.2013cgb, file ="./rds-files/sample_names_2013cgb.rds")
saveRDS(fns.2013cgb, file ="./rds-files/fns_2013cgb.rds")
saveRDS(df.2013cgb, file ="./rds-files/df_2013cgb.rds")

#--- Save the current work and objects:
save.image("R-quality-2013cgb.RData")
rm(list = ls())


### For sequencing run ID "2014" ----

#--- Save job directories:
samples_dir.2014 <- "./sponges-2014/"
quality_dir.2014 <- "./sponges-2014/quality-2014/"  # qual pdf
filtered_dir.2014 <- "./sponges-2014/filtered-2014/"  # fastq filtered

#--- Create the respective folders:
dir.create(path.expand("sponges-2014/quality-2014"))
dir.create(path.expand("sponges-2014/filtered-2014"))

#--- List all fastq single-end files (R) and extract sample names:
fns.2014 <- sort(list.files(samples_dir.2014, full.names = TRUE))
fns.2014 <- fns.2014[str_detect(basename(fns.2014), ".fastq")]
sample.names.2014 <- str_split(basename(fns.2014), pattern = "_", simplify = TRUE)
sample.names.2014 <- sample.names.2014[, 1]

#--- Calcule number of single-end sequences:
df.2014 <- data.frame()
for (i in 1:length(fns.2014)) {
  # use the dada2 function fastq.geometry
  geom.2014 <- fastq.geometry(fns.2014[i])
  # extract the information on number of sequences and file name
  df_one_row.2014 <- data.frame(n_seq = geom.2014[1], file_name = basename(fns.2014[i]))
  # add one line to data frame
  df.2014 <- bind_rows(df.2014, df_one_row.2014)
}

#--- Plot histogram of num seqs:
ggplot(df.2014, aes(x = n_seq)) + geom_histogram(alpha = 0.5, position = "identity", binwidth = 200) + xlim(0, 85000)

#--- Plot quality visualization files:
for (i in 1:length(fns.2014)) {
  # Use dada2 function to plot quality
  qplots.2014 <- plotQualityProfile(fns.2014[i])
  # Save as a pdf
  qplots_png.2014 <- paste0(quality_dir.2014, basename(fns.2014[i]), "-qual.png")
  ggsave(plot = qplots.2014, filename = qplots_png.2014, device = "png", width = 15, height = 15, 
         scale = 1, units = "cm")
}

#--- Save data for downstream steps:
saveRDS(sample.names.2014, file ="./rds-files/sample_names_2014.rds")
saveRDS(fns.2014, file ="./rds-files/fns_2014.rds")
saveRDS(df.2014, file ="./rds-files/df_2014.rds")

#--- Save the current work and objects:
save.image("R-quality-2014.RData")
rm(list = ls())


### For sequencing run ID "2015" ----

#--- Save job directories:
samples_dir.2015 <- "./sponges-2015/"
quality_dir.2015 <- "./sponges-2015/quality-2015/"  # qual pdf
filtered_dir.2015 <- "./sponges-2015/filtered-2015/"  # fastq filtered

#--- Create the respective folders:
dir.create(path.expand("sponges-2015/quality-2015"))
dir.create(path.expand("sponges-2015/filtered-2015"))

#--- List all fastq single-end files (R) and extract sample names:
fns.2015 <- sort(list.files(samples_dir.2015, full.names = TRUE))
fns.2015 <- fns.2015[str_detect(basename(fns.2015), ".fastq")]
sample.names.2015 <- str_split(basename(fns.2015), pattern = "_", simplify = TRUE)
sample.names.2015 <- sample.names.2015[, 1]

#--- Calcule number of single-end sequences:
df.2015 <- data.frame()
for (i in 1:length(fns.2015)) {
  # use the dada2 function fastq.geometry
  geom.2015 <- fastq.geometry(fns.2015[i])
  # extract the information on number of sequences and file name
  df_one_row.2015 <- data.frame(n_seq = geom.2015[1], file_name = basename(fns.2015[i]))
  # add one line to data frame
  df.2015 <- bind_rows(df.2015, df_one_row.2015)
}

#--- Plot histogram of num seqs:
ggplot(df.2015, aes(x = n_seq)) + geom_histogram(alpha = 0.5, position = "identity", binwidth = 200) + xlim(0, 85000)

#--- Plot quality visualization files:
for (i in 1:length(fns.2015)) {
  # Use dada2 function to plot quality
  qplots.2015 <- plotQualityProfile(fns.2015[i])
  # Save as a pdf
  qplots_png.2015 <- paste0(quality_dir.2015, basename(fns.2015[i]), "-qual.png")
  ggsave(plot = qplots.2015, filename = qplots_png.2015, device = "png", width = 15, height = 15, 
         scale = 1, units = "cm")
}

#--- Save data for downstream steps:
saveRDS(sample.names.2015, file ="./rds-files/sample_names_2015.rds")
saveRDS(fns.2015, file ="./rds-files/fns_2015.rds")
saveRDS(df.2015, file ="./rds-files/df_2015.rds")

#--- Save the current work and objects:
save.image("R-quality-2015.RData")
rm(list = ls())


### For sequencing run ID "2020ss" ----

#--- Save job directories:
samples_dir.2020ss <- "./sponges-2020ss/"
quality_dir.2020ss <- "./sponges-2020ss/quality-2020ss/"  # qual pdf
filtered_dir.2020ss <- "./sponges-2020ss/filtered-2020ss/"  # fastq filtered

#--- Create the respective folders:
dir.create(path.expand("sponges-2020ss/quality-2020ss"))
dir.create(path.expand("sponges-2020ss/filtered-2020ss"))

#--- List all fastq single-end files (R) and extract sample names:
fns.2020ss <- sort(list.files(samples_dir.2020ss, full.names = TRUE))
fns.2020ss <- fns.2020ss[str_detect(basename(fns.2020ss), ".fastq")]
sample.names.2020ss <- str_split(basename(fns.2020ss), pattern = "_", simplify = TRUE)
sample.names.2020ss <- sample.names.2020ss[, 1]

#--- Calcule number of single-end sequences:
df.2020ss <- data.frame()
for (i in 1:length(fns.2020ss)) {
  # use the dada2 function fastq.geometry
  geom.2020ss <- fastq.geometry(fns.2020ss[i])
  # extract the information on number of sequences and file name
  df_one_row.2020ss <- data.frame(n_seq = geom.2020ss[1], file_name = basename(fns.2020ss[i]))
  # add one line to data frame
  df.2020ss <- bind_rows(df.2020ss, df_one_row.2020ss)
}

#--- Plot histogram of num seqs:
ggplot(df.2020ss, aes(x = n_seq)) + geom_histogram(alpha = 0.5, position = "identity", binwidth = 200) + xlim(0, 85000)

#--- Plot quality visualization files:
for (i in 1:length(fns.2020ss)) {
  # Use dada2 function to plot quality
  qplots.2020ss <- plotQualityProfile(fns.2020ss[i])
  # Save as a pdf
  qplots_png.2020ss <- paste0(quality_dir.2020ss, basename(fns.2020ss[i]), "-qual.png")
  ggsave(plot = qplots.2020ss, filename = qplots_png.2020ss, device = "png", width = 15, height = 15, 
         scale = 1, units = "cm")
}

#--- Save data for downstream steps:
saveRDS(sample.names.2020ss, file ="./rds-files/sample_names_2020ss.rds")
saveRDS(fns.2020ss, file ="./rds-files/fns_2020ss.rds")
saveRDS(df.2020ss, file ="./rds-files/df_2020ss.rds")

#--- Save the current work and objects:
save.image("R-quality-2020ss.RData")
rm(list = ls())


### For sequencing run ID "h2000l7" ----

#--- Save job directories:
samples_dir.h2000l7 <- "./sponges-h2000l7/"
quality_dir.h2000l7 <- "./sponges-h2000l7/quality-h2000l7/"  # qual pdf
filtered_dir.h2000l7 <- "./sponges-h2000l7/filtered-h2000l7/"  # fastq filtered

#--- Create the respective folders:
dir.create(path.expand("sponges-h2000l7/quality-h2000l7"))
dir.create(path.expand("sponges-h2000l7/filtered-h2000l7"))

#--- List all fastq single-end files (R) and extract sample names:
fns.h2000l7 <- sort(list.files(samples_dir.h2000l7, full.names = TRUE))
fns.h2000l7 <- fns.h2000l7[str_detect(basename(fns.h2000l7), ".fastq")]
sample.names.h2000l7 <- str_split(basename(fns.h2000l7), pattern = "_", simplify = TRUE)
sample.names.h2000l7 <- sample.names.h2000l7[, 1]

#--- Calcule number of single-end sequences:
df.h2000l7 <- data.frame()
for (i in 1:length(fns.h2000l7)) {
  # use the dada2 function fastq.geometry
  geom.h2000l7 <- fastq.geometry(fns.h2000l7[i])
  # extract the information on number of sequences and file name
  df_one_row.h2000l7 <- data.frame(n_seq = geom.h2000l7[1], file_name = basename(fns.h2000l7[i]))
  # add one line to data frame
  df.h2000l7 <- bind_rows(df.h2000l7, df_one_row.h2000l7)
}

#--- Plot histogram of num seqs:
ggplot(df.h2000l7, aes(x = n_seq)) + geom_histogram(alpha = 0.5, position = "identity", binwidth = 200) + xlim(0, 85000)

#--- Plot quality visualization files:
for (i in 1:length(fns.h2000l7)) {
  # Use dada2 function to plot quality
  qplots.h2000l7 <- plotQualityProfile(fns.h2000l7[i])
  # Save as a pdf
  qplots_png.h2000l7 <- paste0(quality_dir.h2000l7, basename(fns.h2000l7[i]), "-qual.png")
  ggsave(plot = qplots.h2000l7, filename = qplots_png.h2000l7, device = "png", width = 15, height = 15, 
         scale = 1, units = "cm")
}

#--- Save data for downstream steps:
saveRDS(sample.names.h2000l7, file ="./rds-files/sample_names_h2000l7.rds")
saveRDS(fns.h2000l7, file ="./rds-files/fns_h2000l7.rds")
saveRDS(df.h2000l7, file ="./rds-files/df_h2000l7.rds")

#--- Save the current work and objects:
save.image("R-quality-h2000l7.RData")
rm(list = ls())


### For sequencing run ID "h2000l8" ----

#--- Save job directories:
samples_dir.h2000l8 <- "./sponges-h2000l8/"
quality_dir.h2000l8 <- "./sponges-h2000l8/quality-h2000l8/"  # qual pdf
filtered_dir.h2000l8 <- "./sponges-h2000l8/filtered-h2000l8/"  # fastq filtered

#--- Create the respective folders:
dir.create(path.expand("sponges-h2000l8/quality-h2000l8"))
dir.create(path.expand("sponges-h2000l8/filtered-h2000l8"))

#--- List all fastq single-end files (R) and extract sample names:
fns.h2000l8 <- sort(list.files(samples_dir.h2000l8, full.names = TRUE))
fns.h2000l8 <- fns.h2000l8[str_detect(basename(fns.h2000l8), ".fastq")]
sample.names.h2000l8 <- str_split(basename(fns.h2000l8), pattern = "_", simplify = TRUE)
sample.names.h2000l8 <- sample.names.h2000l8[, 1]

#--- Calcule number of single-end sequences:
df.h2000l8 <- data.frame()
for (i in 1:length(fns.h2000l8)) {
  # use the dada2 function fastq.geometry
  geom.h2000l8 <- fastq.geometry(fns.h2000l8[i])
  # extract the information on number of sequences and file name
  df_one_row.h2000l8 <- data.frame(n_seq = geom.h2000l8[1], file_name = basename(fns.h2000l8[i]))
  # add one line to data frame
  df.h2000l8 <- bind_rows(df.h2000l8, df_one_row.h2000l8)
}

#--- Plot histogram of num seqs:
ggplot(df.h2000l8, aes(x = n_seq)) + geom_histogram(alpha = 0.5, position = "identity", binwidth = 200) + xlim(0, 85000)

#--- Plot quality visualization files:
for (i in 1:length(fns.h2000l8)) {
  # Use dada2 function to plot quality
  qplots.h2000l8 <- plotQualityProfile(fns.h2000l8[i])
  # Save as a pdf
  qplots_png.h2000l8 <- paste0(quality_dir.h2000l8, basename(fns.h2000l8[i]), "-qual.png")
  ggsave(plot = qplots.h2000l8, filename = qplots_png.h2000l8, device = "png", width = 15, height = 15, 
         scale = 1, units = "cm")
}

#--- Save data for downstream steps:
saveRDS(sample.names.h2000l8, file ="./rds-files/sample_names_h2000l8.rds")
saveRDS(fns.h2000l8, file ="./rds-files/fns_h2000l8.rds")
saveRDS(df.h2000l8, file ="./rds-files/df_h2000l8.rds")

#--- Save the current work and objects:
save.image("R-quality-h2000l8.RData")
rm(list = ls())


### For sequencing run ID "h2500l3" ----

#--- Save job directories:
samples_dir.h2500l3 <- "./sponges-h2500l3/"
quality_dir.h2500l3 <- "./sponges-h2500l3/quality-h2500l3/"  # qual pdf
filtered_dir.h2500l3 <- "./sponges-h2500l3/filtered-h2500l3/"  # fastq filtered

#--- Create the respective folders:
dir.create(path.expand("sponges-h2500l3/quality-h2500l3"))
dir.create(path.expand("sponges-h2500l3/filtered-h2500l3"))

#--- List all fastq single-end files (R) and extract sample names:
fns.h2500l3 <- sort(list.files(samples_dir.h2500l3, full.names = TRUE))
fns.h2500l3 <- fns.h2500l3[str_detect(basename(fns.h2500l3), ".fastq")]
sample.names.h2500l3 <- str_split(basename(fns.h2500l3), pattern = "_", simplify = TRUE)
sample.names.h2500l3 <- sample.names.h2500l3[, 1]

#--- Calcule number of single-end sequences:
df.h2500l3 <- data.frame()
for (i in 1:length(fns.h2500l3)) {
  # use the dada2 function fastq.geometry
  geom.h2500l3 <- fastq.geometry(fns.h2500l3[i])
  # extract the information on number of sequences and file name
  df_one_row.h2500l3 <- data.frame(n_seq = geom.h2500l3[1], file_name = basename(fns.h2500l3[i]))
  # add one line to data frame
  df.h2500l3 <- bind_rows(df.h2500l3, df_one_row.h2500l3)
}

#--- Plot histogram of num seqs:
ggplot(df.h2500l3, aes(x = n_seq)) + geom_histogram(alpha = 0.5, position = "identity", binwidth = 200) + xlim(0, 85000)

#--- Plot quality visualization files:
for (i in 1:length(fns.h2500l3)) {
  # Use dada2 function to plot quality
  qplots.h2500l3 <- plotQualityProfile(fns.h2500l3[i])
  # Save as a pdf
  qplots_png.h2500l3 <- paste0(quality_dir.h2500l3, basename(fns.h2500l3[i]), "-qual.png")
  ggsave(plot = qplots.h2500l3, filename = qplots_png.h2500l3, device = "png", width = 15, height = 15, 
         scale = 1, units = "cm")
}

#--- Save data for downstream steps:
saveRDS(sample.names.h2500l3, file ="./rds-files/sample_names_h2500l3.rds")
saveRDS(fns.h2500l3, file ="./rds-files/fns_h2500l3.rds")
saveRDS(df.h2500l3, file ="./rds-files/df_h2500l3.rds")

#--- Save the current work and objects:
save.image("R-quality-h2500l3.RData")
rm(list = ls())


### For sequencing run ID "h2500l4" ----

#--- Save job directories:
samples_dir.h2500l4 <- "./sponges-h2500l4/"
quality_dir.h2500l4 <- "./sponges-h2500l4/quality-h2500l4/"  # qual pdf
filtered_dir.h2500l4 <- "./sponges-h2500l4/filtered-h2500l4/"  # fastq filtered

#--- Create the respective folders:
dir.create(path.expand("sponges-h2500l4/quality-h2500l4"))
dir.create(path.expand("sponges-h2500l4/filtered-h2500l4"))

#--- List all fastq single-end files (R) and extract sample names:
fns.h2500l4 <- sort(list.files(samples_dir.h2500l4, full.names = TRUE))
fns.h2500l4 <- fns.h2500l4[str_detect(basename(fns.h2500l4), ".fastq")]
sample.names.h2500l4 <- str_split(basename(fns.h2500l4), pattern = "_", simplify = TRUE)
sample.names.h2500l4 <- sample.names.h2500l4[, 1]

#--- Calcule number of single-end sequences:
df.h2500l4 <- data.frame()
for (i in 1:length(fns.h2500l4)) {
  # use the dada2 function fastq.geometry
  geom.h2500l4 <- fastq.geometry(fns.h2500l4[i])
  # extract the information on number of sequences and file name
  df_one_row.h2500l4 <- data.frame(n_seq = geom.h2500l4[1], file_name = basename(fns.h2500l4[i]))
  # add one line to data frame
  df.h2500l4 <- bind_rows(df.h2500l4, df_one_row.h2500l4)
}

#--- Plot histogram of num seqs:
ggplot(df.h2500l4, aes(x = n_seq)) + geom_histogram(alpha = 0.5, position = "identity", binwidth = 200) + xlim(0, 85000)

#--- Plot quality visualization files:
for (i in 1:length(fns.h2500l4)) {
  # Use dada2 function to plot quality
  qplots.h2500l4 <- plotQualityProfile(fns.h2500l4[i])
  # Save as a pdf
  qplots_png.h2500l4 <- paste0(quality_dir.h2500l4, basename(fns.h2500l4[i]), "-qual.png")
  ggsave(plot = qplots.h2500l4, filename = qplots_png.h2500l4, device = "png", width = 15, height = 15, 
         scale = 1, units = "cm")
}

#--- Save data for downstream steps:
saveRDS(sample.names.h2500l4, file ="./rds-files/sample_names_h2500l4.rds")
saveRDS(fns.h2500l4, file ="./rds-files/fns_h2500l4.rds")
saveRDS(df.h2500l4, file ="./rds-files/df_h2500l4.rds")

#--- Save the current work and objects:
save.image("R-quality-h2500l4.RData")
rm(list = ls())


### For sequencing run ID "h2500l5" ----

#--- Save job directories:
samples_dir.h2500l5 <- "./sponges-h2500l5/"
quality_dir.h2500l5 <- "./sponges-h2500l5/quality-h2500l5/"  # qual pdf
filtered_dir.h2500l5 <- "./sponges-h2500l5/filtered-h2500l5/"  # fastq filtered

#--- Create the respective folders:
dir.create(path.expand("sponges-h2500l5/quality-h2500l5"))
dir.create(path.expand("sponges-h2500l5/filtered-h2500l5"))

#--- List all fastq single-end files (R) and extract sample names:
fns.h2500l5 <- sort(list.files(samples_dir.h2500l5, full.names = TRUE))
fns.h2500l5 <- fns.h2500l5[str_detect(basename(fns.h2500l5), ".fastq")]
sample.names.h2500l5 <- str_split(basename(fns.h2500l5), pattern = "_", simplify = TRUE)
sample.names.h2500l5 <- sample.names.h2500l5[, 1]

#--- Calcule number of single-end sequences:
df.h2500l5 <- data.frame()
for (i in 1:length(fns.h2500l5)) {
  # use the dada2 function fastq.geometry
  geom.h2500l5 <- fastq.geometry(fns.h2500l5[i])
  # extract the information on number of sequences and file name
  df_one_row.h2500l5 <- data.frame(n_seq = geom.h2500l5[1], file_name = basename(fns.h2500l5[i]))
  # add one line to data frame
  df.h2500l5 <- bind_rows(df.h2500l5, df_one_row.h2500l5)
}

#--- Plot histogram of num seqs:
ggplot(df.h2500l5, aes(x = n_seq)) + geom_histogram(alpha = 0.5, position = "identity", binwidth = 200) + xlim(0, 85000)

#--- Plot quality visualization files:
for (i in 1:length(fns.h2500l5)) {
  # Use dada2 function to plot quality
  qplots.h2500l5 <- plotQualityProfile(fns.h2500l5[i])
  # Save as a pdf
  qplots_png.h2500l5 <- paste0(quality_dir.h2500l5, basename(fns.h2500l5[i]), "-qual.png")
  ggsave(plot = qplots.h2500l5, filename = qplots_png.h2500l5, device = "png", width = 15, height = 15, 
         scale = 1, units = "cm")
}

#--- Save data for downstream steps:
saveRDS(sample.names.h2500l5, file ="./rds-files/sample_names_h2500l5.rds")
saveRDS(fns.h2500l5, file ="./rds-files/fns_h2500l5.rds")
saveRDS(df.h2500l5, file ="./rds-files/df_h2500l5.rds")

#--- Save the current work and objects:
save.image("R-quality-h2500l5.RData")
rm(list = ls())


### For sequencing run ID "h2500l7" ----

#--- Save job directories:
samples_dir.h2500l7 <- "./sponges-h2500l7/"
quality_dir.h2500l7 <- "./sponges-h2500l7/quality-h2500l7/"  # qual pdf
filtered_dir.h2500l7 <- "./sponges-h2500l7/filtered-h2500l7/"  # fastq filtered

#--- Create the respective folders:
dir.create(path.expand("sponges-h2500l7/quality-h2500l7"))
dir.create(path.expand("sponges-h2500l7/filtered-h2500l7"))

#--- List all fastq single-end files (R) and extract sample names:
fns.h2500l7 <- sort(list.files(samples_dir.h2500l7, full.names = TRUE))
fns.h2500l7 <- fns.h2500l7[str_detect(basename(fns.h2500l7), ".fastq")]
sample.names.h2500l7 <- str_split(basename(fns.h2500l7), pattern = "_", simplify = TRUE)
sample.names.h2500l7 <- sample.names.h2500l7[, 1]

#--- Calcule number of single-end sequences:
df.h2500l7 <- data.frame()
for (i in 1:length(fns.h2500l7)) {
  # use the dada2 function fastq.geometry
  geom.h2500l7 <- fastq.geometry(fns.h2500l7[i])
  # extract the information on number of sequences and file name
  df_one_row.h2500l7 <- data.frame(n_seq = geom.h2500l7[1], file_name = basename(fns.h2500l7[i]))
  # add one line to data frame
  df.h2500l7 <- bind_rows(df.h2500l7, df_one_row.h2500l7)
}

#--- Plot histogram of num seqs:
ggplot(df.h2500l7, aes(x = n_seq)) + geom_histogram(alpha = 0.5, position = "identity", binwidth = 200) + xlim(0, 85000)

#--- Plot quality visualization files:
for (i in 1:length(fns.h2500l7)) {
  # Use dada2 function to plot quality
  qplots.h2500l7 <- plotQualityProfile(fns.h2500l7[i])
  # Save as a pdf
  qplots_png.h2500l7 <- paste0(quality_dir.h2500l7, basename(fns.h2500l7[i]), "-qual.png")
  ggsave(plot = qplots.h2500l7, filename = qplots_png.h2500l7, device = "png", width = 15, height = 15, 
         scale = 1, units = "cm")
}

#--- Save data for downstream steps:
saveRDS(sample.names.h2500l7, file ="./rds-files/sample_names_h2500l7.rds")
saveRDS(fns.h2500l7, file ="./rds-files/fns_h2500l7.rds")
saveRDS(df.h2500l7, file ="./rds-files/df_h2500l7.rds")

#--- Save the current work and objects:
save.image("R-quality-h2500l7.RData")
rm(list = ls())

