# Build data/dat.mcduff2019.rda from the raw McDuff & Girard (2019) MTurk ratings.
# Source: file `04_MTurk.csv` at https://osf.io/n4grd/ (saved here as mcduff2019.csv).
# Run from the package root.
dat.mcduff2019 <- read.csv("data_raw/mcduff2019.csv")
colnames(dat.mcduff2019) <- c("item", "judge", "rating_positive", "rating_smile")
dat.mcduff2019 <- dat.mcduff2019[, c("rating_positive", "rating_smile", "item", "judge")]
save(dat.mcduff2019, file = "data/dat.mcduff2019.rda")
rm(dat.mcduff2019)
