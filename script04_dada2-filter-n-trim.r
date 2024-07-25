# FILTERING AND TRIMMING READS #
# ---------------------------- #

## Per sequencing run basis

#--- Load environments needed ----

load("R-quality-1821.RData")
load("R-quality-2013.RData")
load("R-quality-2013cgb.RData")
load("R-quality-2014.RData")
load("R-quality-2015.RData")
load("R-quality-2020ss.RData")
load("R-quality-h2000l7.RData")
load("R-quality-h2000l8.RData")
load("R-quality-h2500l3.RData")
load("R-quality-h2500l4.RData")
load("R-quality-h2500l5.RData")
load("R-quality-h2500l7.RData")


#--- Save directories for the filtered files ----

filt.1821 <- str_c(filtered_dir.1821, sample.names.1821, "_filt.fastq")
filt.2013 <- str_c(filtered_dir.2013, sample.names.2013, "_filt.fastq")
filt.2013cgb <- str_c(filtered_dir.2013cgb, sample.names.2013cgb, "_filt.fastq")
filt.2014 <- str_c(filtered_dir.2014, sample.names.2014, "_filt.fastq")
filt.2015 <- str_c(filtered_dir.2015, sample.names.2015, "_filt.fastq")
filt.2020ss <- str_c(filtered_dir.1821, sample.names.2020ss, "_filt.fastq")
filt.h2000l7 <- str_c(filtered_dir.h2000l7, sample.names.h2000l7, "_filt.fastq")
filt.h2000l8 <- str_c(filtered_dir.h2000l8, sample.names.h2000l8, "_filt.fastq")
filt.h2500l3 <- str_c(filtered_dir.h2500l3, sample.names.h2500l3, "_filt.fastq")
filt.h2500l4 <- str_c(filtered_dir.h2500l4, sample.names.h2500l4, "_filt.fastq")
filt.h2500l5 <- str_c(filtered_dir.h2500l5, sample.names.h2500l5, "_filt.fastq")
filt.h2500l7 <- str_c(filtered_dir.h2500l7, sample.names.h2500l7, "_filt.fastq")


### Filtering/trimming reads by quality ----

out.1821 <- filterAndTrim(fns.1821, filt.1821, truncLen = 125, minLen = 100, 
                          maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE, 
                          compress = FALSE, multithread = 70)
out.2013 <- filterAndTrim(fns.2013, filt.2013, truncLen = 125, minLen = 100, # ver los plots de calidad para determinar donde se cae
                          maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE, ### maxN's = 0 is VERY IMPORTANT! # phix es un control interno + de ILLUMINA # maxEE rango de calidad min por aceptar (como q20)
                     compress = FALSE, multithread = 70)
out.2013cgb <- filterAndTrim(fns.2013cgb, filt.2013cgb, truncLen = 125, minLen = 100, 
                          maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE, 
                          compress = FALSE, multithread = 70)
out.2014 <- filterAndTrim(fns.2014, filt.2014, truncLen = 125, minLen = 100, 
                          maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE, 
                          compress = FALSE, multithread = 70)
out.2015 <- filterAndTrim(fns.2015, filt.2015, truncLen = 125, minLen = 100, 
                          maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE, 
                          compress = FALSE, multithread = 70)
out.2020ss <- filterAndTrim(fns.2020ss, filt.2020ss, truncLen = 125, minLen = 100, 
                          maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE, 
                          compress = FALSE, multithread = 70)
out.h2000l7 <- filterAndTrim(fns.h2000l7, filt.h2000l7, minLen = 100, 
                          maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE, 
                          compress = FALSE, multithread = 70)
out.h2000l8 <- filterAndTrim(fns.h2000l8, filt.h2000l8, minLen = 100, 
                          maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE, 
                          compress = FALSE, multithread = 70)
out.h2500l3 <- filterAndTrim(fns.h2500l3, filt.h2500l3, minLen = 100, 
                          maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE, 
                          compress = FALSE, multithread = 70)
out.h2500l4 <- filterAndTrim(fns.h2500l4, filt.h2500l4, minLen = 100, 
                          maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE, 
                          compress = FALSE, multithread = 70)
out.h2500l5 <- filterAndTrim(fns.h2500l5, filt.h2500l5, minLen = 100, 
                          maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE, 
                          compress = FALSE, multithread = 70)
out.h2500l7 <- filterAndTrim(fns.h2500l7, filt.h2500l7, minLen = 100, 
                          maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE, 
                          compress = FALSE, multithread = 70)
#Note: "out" contain reads.in and reads.out.


#--- Save things -----------------------------------

# Save data for downstream steps:
saveRDS(filt.1821, file ="./rds-files/filt_1821.rds")
saveRDS(filt.2013, file ="./rds-files/filt_2013.rds")
saveRDS(filt.2013cgb, file ="./rds-files/filt_2013cgb.rds")
saveRDS(filt.2014, file ="./rds-files/filt_2014.rds")
saveRDS(filt.2015, file ="./rds-files/filt_2015.rds")
saveRDS(filt.2020ss, file ="./rds-files/filt_2020ss.rds")
saveRDS(filt.h2000l7, file ="./rds-files/filt_h2000l7.rds")
saveRDS(filt.h2000l8, file ="./rds-files/filt_h2000l8.rds")
saveRDS(filt.h2500l3, file ="./rds-files/filt_h2500l3.rds")
saveRDS(filt.h2500l4, file ="./rds-files/filt_h2500l4.rds")
saveRDS(filt.h2500l5, file ="./rds-files/filt_h2500l5.rds")
saveRDS(filt.h2500l7, file ="./rds-files/filt_h2500l7.rds")

saveRDS(out.1821, file ="./rds-files/out_1821.rds")
saveRDS(out.2013, file ="./rds-files/out_2013.rds")
saveRDS(out.2013cgb, file ="./rds-files/out_2013cgb.rds")
saveRDS(out.2014, file ="./rds-files/out_2014.rds")
saveRDS(out.2015, file ="./rds-files/out_2015.rds")
saveRDS(out.2020ss, file ="./rds-files/out_2020ss.rds")
saveRDS(out.h2000l7, file ="./rds-files/out_h2000l7.rds")
saveRDS(out.h2000l8, file ="./rds-files/out_h2000l8.rds")
saveRDS(out.h2500l3, file ="./rds-files/out_h2500l3.rds")
saveRDS(out.h2500l4, file ="./rds-files/out_h2500l4.rds")
saveRDS(out.h2500l5, file ="./rds-files/out_h2500l5.rds")
saveRDS(out.h2500l7, file ="./rds-files/out_h2500l7.rds")

# Extra, write a table of initial seqs:
dfs <- bind_rows(df.1821, df.2013, df.2013cgb, df.2014, df.2015, df.2020ss, df.h2000l7, df.h2000l8, df.h2500l3, df.h2500l4, df.h2500l5, df.h2500l7)
write.table(dfs, file = './results/input-seqs.tsv', sep = '\t', row.names = FALSE, na = '', quote = FALSE)

# Save the current work and objects:
save.image("R-filtntrim.RData")

