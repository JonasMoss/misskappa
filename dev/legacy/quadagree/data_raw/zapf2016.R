dat.zapf2016 = readxl::read_excel("data_raw/zapf2016.xlsx") + 1
save(dat.zapf2016, file = "data/dat.zapf2016.rda")
rm(dat.zapf2016)
