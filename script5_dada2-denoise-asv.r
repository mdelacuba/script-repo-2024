# DENOISING SEQUENCES AND ASV INFERENCE #
# ------------------------------------- #

## Per sequencing run basis

### For sequencing run ID "1821" ----

#--- Load data needed for this step:
filt.1821 <- readRDS(file = "./rds-files/filt_1821.rds")
sample.names.1821 <- readRDS(file = "./rds-files/sample_names_1821.rds")

#--- Evaluate error rates:
err.1821 <- learnErrors(filt.1821, multithread = 6)
err.plot.1821 <- plotErrors(err.1821, nominalQ = TRUE)

#--- Dereplication:
derep.1821 <- derepFastq(filt.1821, verbose = TRUE)

names(derep.1821) <- sample.names.1821

#--- Inference of ASVs:
dada.1821 <- dada(derep.1821, err = err.1821, multithread = 6, pool = F)

#--- Make sequence table:
seqtab.1821 <- makeSequenceTable(dada.1821)
#dim(seqtab)

#--- Truncate ASVs length to 100 bp:
#colnames(seqtab.1821) <- substr(colnames(seqtab.1821), 1, 100)

#--- Inspect distribution of sequence lengths:
distrib.seqslength.1821 <- table(nchar(getSequences(seqtab.1821)))

#--- Save data for downstream steps:
saveRDS(seqtab.1821, file ="./rds-files/seqtab_1821.rds")
saveRDS(dada.1821, file ="./rds-files/dada_1821.rds")

#--- Save the current work and objects:
save.image("R-denoise-1821.RData")
rm(list = ls())


### For sequencing run ID "2013" ----


#--- Load data needed for this step:
filt.2013 <- readRDS(file = "./rds-files/filt_2013.rds")
sample.names.2013 <- readRDS(file = "./rds-files/sample_names_2013.rds")

#--- Evaluate error rates:
err.2013 <- learnErrors(filt.2013, multithread = 6)
err.plot.2013 <- plotErrors(err.2013, nominalQ = TRUE)

#--- Dereplication:
derep.2013 <- derepFastq(filt.2013, verbose = TRUE)

names(derep.2013) <- sample.names.2013

#--- Inference of ASVs:
dada.2013 <- dada(derep.2013, err = err.2013, multithread = 6, pool = F)

#--- Make sequence table:
seqtab.2013 <- makeSequenceTable(dada.2013)
#dim(seqtab)

#--- Truncate ASVs length to 100 bp:
#colnames(seqtab.2013) <- substr(colnames(seqtab.2013), 1, 100)

#--- Inspect distribution of sequence lengths:
distrib.seqslength.2013 <- table(nchar(getSequences(seqtab.2013)))

#--- Save data for downstream steps:
saveRDS(seqtab.2013, file ="./rds-files/seqtab_2013.rds")
saveRDS(dada.2013, file ="./rds-files/dada_2013.rds")

#--- Save the current work and objects:
save.image("R-denoise-2013.RData")
rm(list = ls())


### For sequencing run ID "2013cgb" ----

#--- Load data needed for this step:
filt.2013cgb <- readRDS(file = "./rds-files/filt_2013cgb.rds")
sample.names.2013cgb <- readRDS(file = "./rds-files/sample_names_2013cgb.rds")

#--- Evaluate error rates:
err.2013cgb <- learnErrors(filt.2013cgb, multithread = 6)
err.plot.2013cgb <- plotErrors(err.2013cgb, nominalQ = TRUE)

#--- Dereplication:
derep.2013cgb <- derepFastq(filt.2013cgb, verbose = TRUE)

names(derep.2013cgb) <- sample.names.2013cgb

#--- Inference of ASVs:
dada.2013cgb <- dada(derep.2013cgb, err = err.2013cgb, multithread = 6, pool = F)

#--- Make sequence table:
seqtab.2013cgb <- makeSequenceTable(dada.2013cgb)
#dim(seqtab)

#--- Truncate ASVs length to 100 bp:
#colnames(seqtab.2013cgb) <- substr(colnames(seqtab.2013cgb), 1, 100)

#--- Inspect distribution of sequence lengths:
distrib.seqslength.2013cgb <- table(nchar(getSequences(seqtab.2013cgb)))

#--- Save data for downstream steps:
saveRDS(seqtab.2013cgb, file ="./rds-files/seqtab_2013cgb.rds")
saveRDS(dada.2013cgb, file ="./rds-files/dada_2013cgb.rds")

#--- Save the current work and objects:
save.image("R-denoise-2013cgb.RData")
rm(list = ls())


### For sequencing run ID "2014" ----

#--- Load data needed for this step:
filt.2014 <- readRDS(file = "./rds-files/filt_2014.rds")
sample.names.2014 <- readRDS(file = "./rds-files/sample_names_2014.rds")

#--- Evaluate error rates:
err.2014 <- learnErrors(filt.2014, multithread = 6)
err.plot.2014 <- plotErrors(err.2014, nominalQ = TRUE)

#--- Dereplication:
derep.2014 <- derepFastq(filt.2014, verbose = TRUE)

names(derep.2014) <- sample.names.2014

#--- Inference of ASVs:
dada.2014 <- dada(derep.2014, err = err.2014, multithread = 6, pool = F)

#--- Make sequence table:
seqtab.2014 <- makeSequenceTable(dada.2014)
#dim(seqtab)

#--- Truncate ASVs length to 100 bp:
#colnames(seqtab.2014) <- substr(colnames(seqtab.2014), 1, 100)

#--- Inspect distribution of sequence lengths:
distrib.seqslength.2014 <- table(nchar(getSequences(seqtab.2014)))

#--- Save data for downstream steps:
saveRDS(seqtab.2014, file ="./rds-files/seqtab_2014.rds")
saveRDS(dada.2014, file ="./rds-files/dada_2014.rds")

#--- Save the current work and objects:
save.image("R-denoise-2014.RData")
rm(list = ls())


### For sequencing run ID "2015" ----

#--- Load data needed for this step:
filt.2015 <- readRDS(file = "./rds-files/filt_2015.rds")
sample.names.2015 <- readRDS(file = "./rds-files/sample_names_2015.rds")

#--- Evaluate error rates:
err.2015 <- learnErrors(filt.2015, multithread = 6)
err.plot.2015 <- plotErrors(err.2015, nominalQ = TRUE)

#--- Dereplication:
derep.2015 <- derepFastq(filt.2015, verbose = TRUE)

names(derep.2015) <- sample.names.2015

#--- Inference of ASVs:
dada.2015 <- dada(derep.2015, err = err.2015, multithread = 6, pool = F)

#--- Make sequence table:
seqtab.2015 <- makeSequenceTable(dada.2015)
#dim(seqtab)

#--- Truncate ASVs length to 100 bp:
#colnames(seqtab.2015) <- substr(colnames(seqtab.2015), 1, 100)

#--- Inspect distribution of sequence lengths:
distrib.seqslength.2015 <- table(nchar(getSequences(seqtab.2015)))

#--- Save data for downstream steps:
saveRDS(seqtab.2015, file ="./rds-files/seqtab_2015.rds")
saveRDS(dada.2015, file ="./rds-files/dada_2015.rds")

#--- Save the current work and objects:
save.image("R-denoise-2015.RData")
rm(list = ls())


### For sequencing run ID "2020ss" ----

#--- Load data needed for this step:
filt.2020ss <- readRDS(file = "./rds-files/filt_2020ss.rds")
sample.names.2020ss <- readRDS(file = "./rds-files/sample_names_2020ss.rds")

#--- Evaluate error rates:
err.2020ss <- learnErrors(filt.2020ss, multithread = 70)
err.plot.2020ss <- plotErrors(err.2020ss, nominalQ = TRUE)

#--- Dereplication:
derep.2020ss <- derepFastq(filt.2020ss, verbose = TRUE)

names(derep.2020ss) <- sample.names.2020ss

#--- Inference of ASVs:
dada.2020ss <- dada(derep.2020ss, err = err.2020ss, multithread = 70, pool = F)

#--- Make sequence table:
seqtab.2020ss <- makeSequenceTable(dada.2020ss)
#dim(seqtab)

#--- Truncate ASVs length to 100 bp:
#colnames(seqtab.2020ss) <- substr(colnames(seqtab.2020ss), 1, 100)

#--- Inspect distribution of sequence lengths:
distrib.seqslength.2020ss <- table(nchar(getSequences(seqtab.2020ss)))

#--- Save data for downstream steps:
saveRDS(seqtab.2020ss, file ="./rds-files/seqtab_2020ss.rds")
saveRDS(dada.2020ss, file ="./rds-files/dada_2020ss.rds")

#--- Save the current work and objects:
save.image("R-denoise-2020ss.RData")
rm(list = ls())


### For sequencing run ID "h2000l7" ----

#--- Load data needed for this step:
filt.h2000l7 <- readRDS(file = "./rds-files/filt_h2000l7.rds")
sample.names.h2000l7 <- readRDS(file = "./rds-files/sample_names_h2000l7.rds")

#--- Evaluate error rates:
err.h2000l7 <- learnErrors(filt.h2000l7, multithread = 6)
err.plot.h2000l7 <- plotErrors(err.h2000l7, nominalQ = TRUE)

#--- Dereplication:
derep.h2000l7 <- derepFastq(filt.h2000l7, verbose = TRUE)

names(derep.h2000l7) <- sample.names.h2000l7

#--- Inference of ASVs:
dada.h2000l7 <- dada(derep.h2000l7, err = err.h2000l7, multithread = 6, pool = F)

#--- Make sequence table:
seqtab.h2000l7 <- makeSequenceTable(dada.h2000l7)
#dim(seqtab)

#--- Truncate ASVs length to 100 bp:
#colnames(seqtab.h2000l7) <- substr(colnames(seqtab.h2000l7), 1, 100)

#--- Inspect distribution of sequence lengths:
distrib.seqslength.h2000l7 <- table(nchar(getSequences(seqtab.h2000l7)))

#--- Save data for downstream steps:
saveRDS(seqtab.h2000l7, file ="./rds-files/seqtab_h2000l7.rds")
saveRDS(dada.h2000l7, file ="./rds-files/dada_h2000l7.rds")

#--- Save the current work and objects:
save.image("R-denoise-h2000l7.RData")
rm(list = ls())


### For sequencing run ID "h2000l8" ----

#--- Load data needed for this step:
filt.h2000l8 <- readRDS(file = "./rds-files/filt_h2000l8.rds")
sample.names.h2000l8 <- readRDS(file = "./rds-files/sample_names_h2000l8.rds")

#--- Evaluate error rates:
err.h2000l8 <- learnErrors(filt.h2000l8, multithread = 6)
err.plot.h2000l8 <- plotErrors(err.h2000l8, nominalQ = TRUE)

#--- Dereplication:
derep.h2000l8 <- derepFastq(filt.h2000l8, verbose = TRUE)

names(derep.h2000l8) <- sample.names.h2000l8

#--- Inference of ASVs:
dada.h2000l8 <- dada(derep.h2000l8, err = err.h2000l8, multithread = 6, pool = F)

#--- Make sequence table:
seqtab.h2000l8 <- makeSequenceTable(dada.h2000l8)
#dim(seqtab)

#--- Truncate ASVs length to 100 bp:
#colnames(seqtab.h2000l8) <- substr(colnames(seqtab.h2000l8), 1, 100)

#--- Inspect distribution of sequence lengths:
distrib.seqslength.h2000l8 <- table(nchar(getSequences(seqtab.h2000l8)))

#--- Save data for downstream steps:
saveRDS(seqtab.h2000l8, file ="./rds-files/seqtab_h2000l8.rds")
saveRDS(dada.h2000l8, file ="./rds-files/dada_h2000l8.rds")

#--- Save the current work and objects:
save.image("R-denoise-h2000l8.RData")
rm(list = ls())


### For sequencing run ID "h2500l3" ----

#--- Load data needed for this step:
filt.h2500l3 <- readRDS(file = "./rds-files/filt_h2500l3.rds")
sample.names.h2500l3 <- readRDS(file = "./rds-files/sample_names_h2500l3.rds")

#--- Evaluate error rates:
err.h2500l3 <- learnErrors(filt.h2500l3, multithread = 6)
err.plot.h2500l3 <- plotErrors(err.h2500l3, nominalQ = TRUE)

#--- Dereplication:
derep.h2500l3 <- derepFastq(filt.h2500l3, verbose = TRUE)

names(derep.h2500l3) <- sample.names.h2500l3

#--- Inference of ASVs:
dada.h2500l3 <- dada(derep.h2500l3, err = err.h2500l3, multithread = 6, pool = F)

#--- Make sequence table:
seqtab.h2500l3 <- makeSequenceTable(dada.h2500l3)
#dim(seqtab)

#--- Truncate ASVs length to 100 bp:
#colnames(seqtab.h2500l3) <- substr(colnames(seqtab.h2500l3), 1, 100)

#--- Inspect distribution of sequence lengths:
distrib.seqslength.h2500l3 <- table(nchar(getSequences(seqtab.h2500l3)))

#--- Save data for downstream steps:
saveRDS(seqtab.h2500l3, file ="./rds-files/seqtab_h2500l3.rds")
saveRDS(dada.h2500l3, file ="./rds-files/dada_h2500l3.rds")

#--- Save the current work and objects:
save.image("R-denoise-h2500l3.RData")
rm(list = ls())


### For sequencing run ID "h2500l4" ----

#--- Load data needed for this step:
filt.h2500l4 <- readRDS(file = "./rds-files/filt_h2500l4.rds")
sample.names.h2500l4 <- readRDS(file = "./rds-files/sample_names_h2500l4.rds")

#--- Evaluate error rates:
err.h2500l4 <- learnErrors(filt.h2500l4, multithread = 6)
err.plot.h2500l4 <- plotErrors(err.h2500l4, nominalQ = TRUE)

#--- Dereplication:
derep.h2500l4 <- derepFastq(filt.h2500l4, verbose = TRUE)

names(derep.h2500l4) <- sample.names.h2500l4

#--- Inference of ASVs:
dada.h2500l4 <- dada(derep.h2500l4, err = err.h2500l4, multithread = 6, pool = F)

#--- Make sequence table:
seqtab.h2500l4 <- makeSequenceTable(dada.h2500l4)
#dim(seqtab)

#--- Truncate ASVs length to 100 bp:
#colnames(seqtab.h2500l4) <- substr(colnames(seqtab.h2500l4), 1, 100)

#--- Inspect distribution of sequence lengths:
distrib.seqslength.h2500l4 <- table(nchar(getSequences(seqtab.h2500l4)))

#--- Save data for downstream steps:
saveRDS(seqtab.h2500l4, file ="./rds-files/seqtab_h2500l4.rds")
saveRDS(dada.h2500l4, file ="./rds-files/dada_h2500l4.rds")

#--- Save the current work and objects:
save.image("R-denoise-h2500l4.RData")
rm(list = ls())


### For sequencing run ID "h2500l5" ----

#--- Load data needed for this step:
filt.h2500l5 <- readRDS(file = "./rds-files/filt_h2500l5.rds")
sample.names.h2500l5 <- readRDS(file = "./rds-files/sample_names_h2500l5.rds")

#--- Evaluate error rates:
err.h2500l5 <- learnErrors(filt.h2500l5, multithread = 6)
err.plot.h2500l5 <- plotErrors(err.h2500l5, nominalQ = TRUE)

#--- Dereplication:
derep.h2500l5 <- derepFastq(filt.h2500l5, verbose = TRUE)

names(derep.h2500l5) <- sample.names.h2500l5

#--- Inference of ASVs:
dada.h2500l5 <- dada(derep.h2500l5, err = err.h2500l5, multithread = 6, pool = F)

#--- Make sequence table:
seqtab.h2500l5 <- makeSequenceTable(dada.h2500l5)
#dim(seqtab)

#--- Truncate ASVs length to 100 bp:
#colnames(seqtab.h2500l5) <- substr(colnames(seqtab.h2500l5), 1, 100)

#--- Inspect distribution of sequence lengths:
distrib.seqslength.h2500l5 <- table(nchar(getSequences(seqtab.h2500l5)))

#--- Save data for downstream steps:
saveRDS(seqtab.h2500l5, file ="./rds-files/seqtab_h2500l5.rds")
saveRDS(dada.h2500l5, file ="./rds-files/dada_h2500l5.rds")

#--- Save the current work and objects:
save.image("R-denoise-h2500l5.RData")
rm(list = ls())


### For sequencing run ID "h2500l7" ----

#--- Load data needed for this step:
filt.h2500l7 <- readRDS(file = "./rds-files/filt_h2500l7.rds")
sample.names.h2500l7 <- readRDS(file = "./rds-files/sample_names_h2500l7.rds")

#--- Evaluate error rates:
err.h2500l7 <- learnErrors(filt.h2500l7, multithread = 6)
err.plot.h2500l7 <- plotErrors(err.h2500l7, nominalQ = TRUE)

#--- Dereplication:
derep.h2500l7 <- derepFastq(filt.h2500l7, verbose = TRUE)

names(derep.h2500l7) <- sample.names.h2500l7

#--- Inference of ASVs:
dada.h2500l7 <- dada(derep.h2500l7, err = err.h2500l7, multithread = 6, pool = F)

#--- Make sequence table:
seqtab.h2500l7 <- makeSequenceTable(dada.h2500l7)
#dim(seqtab)

#--- Truncate ASVs length to 100 bp:
#colnames(seqtab.h2500l7) <- substr(colnames(seqtab.h2500l7), 1, 100)

#--- Inspect distribution of sequence lengths:
distrib.seqslength.h2500l7 <- table(nchar(getSequences(seqtab.h2500l7)))

#--- Save data for downstream steps:
saveRDS(seqtab.h2500l7, file ="./rds-files/seqtab_h2500l7.rds")
saveRDS(dada.h2500l7, file ="./rds-files/dada_h2500l7.rds")

#--- Save the current work and objects:
save.image("R-denoise-h2500l7.RData")
rm(list = ls())

