

# Examples using my functions


devtools::install_github("monteroserra/rOceansNew")

library(rOceansNew)

setwd("C:/Users/monteroserra/Desktop/rOceans Google Summer Code/Example_Acropora")


Acropora_OBIS = read.csv2("Acropora_OBIS.csv", stringsAsFactors = F)
Acropora_GBIF = read.csv2("Acropora_GBIF.csv", stringsAsFactors = F)


Acropora_OBIS$decimalLongitude = as.numeric(Acropora_OBIS$decimalLongitude)
Acropora_OBIS$decimalLatitude = as.numeric(Acropora_OBIS$decimalLatitude)



Acropora_OBIS_checked = oceanDataCheck(OBIS_occurrences = Acropora_OBIS, source = "OBIS")
Acropora_GBIF_checked = oceanDataCheck(GBIF_occurrences = Acropora_GBIF, source = "GBIF")

colnames()


