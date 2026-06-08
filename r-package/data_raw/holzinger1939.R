# Build data/dat.holzinger1939.rda from the Holzinger & Swineford (1939) battery.
# Source: the `HolzingerSwineford1939` data frame distributed with the 'lavaan'
# package (Rosseel, 2012); originally Holzinger, K. J., & Swineford, F. (1939).
# Saved here as holzinger1939.csv for a network-free, reproducible build. The
# data are factual 1939 mental-ability test scores and are in the public domain.
# Run from the package root.
dat.holzinger1939 <- read.csv("data_raw/holzinger1939.csv", stringsAsFactors = FALSE)
dat.holzinger1939$school <- factor(dat.holzinger1939$school,
                                   levels = c("Pasteur", "Grant-White"))
save(dat.holzinger1939, file = "data/dat.holzinger1939.rda")
rm(dat.holzinger1939)
