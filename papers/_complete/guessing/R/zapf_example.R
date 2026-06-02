

models = c("cf", "fleiss" ,"cohen","bp","cbp")

format(sapply(models, function(model)
  agreeable::knowledge(agreeable::dat.zapf2016, model = model, skip = TRUE),
  USE.NAMES = TRUE), digits = 2)

knitr::kable(format(sapply(models, function(model)
  agreeable::knowledge(agreeable::dat.zapf2016, model = model)$conf.int,
  USE.NAMES = TRUE), digits = 2))