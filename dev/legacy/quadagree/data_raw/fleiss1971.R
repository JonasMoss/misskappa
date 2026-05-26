dat.fleiss1971 = read.csv("data_raw/fleiss1971.csv")
class(dat.fleiss1971) = c("fleiss_form", class(dat.fleiss1971))
colnames(dat.fleiss1971) = c("depression", "personality disorder", "schizophrenia", "neurosis", "other")
dat.fleiss1971[is.na(dat.fleiss1971)] = 0
dat.fleiss1971 = tibble::as_tibble(dat.fleiss1971)
save(dat.fleiss1971, file = "data/dat.fleiss1971.rda")
rm(dat.fleiss1971)
