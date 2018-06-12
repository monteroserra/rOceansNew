

# Examples using my functions


devtools::install_github("monteroserra/rOceansNew")

library(rOceansNew)

setwd("C:/Users/monteroserra/Desktop/rOceans Google Summer Code/Example_Acropora")


Acropora_OBIS = read.csv2("Acropora_OBIS.csv", stringsAsFactors = F)
Acropora_GBIF = read.csv2("Acropora_GBIF.csv", stringsAsFactors = F)


Acropora_OBIS$decimalLongitude = as.numeric(Acropora_OBIS$decimalLongitude)
Acropora_OBIS$decimalLatitude = as.numeric(Acropora_OBIS$decimalLatitude)


#Checking, filtering and merging big datasets downloaded directly from OBIS & GBIF

Acropora_OBIS_checked = oceanDataCheck(OBIS_occurrences = Acropora_OBIS, source = "OBIS")
Acropora_GBIF_checked = oceanDataCheck(GBIF_occurrences = Acropora_GBIF, source = "GBIF")

Acropora_Total_Checked = oceanDataCheck(GBIF_occurrences = Acropora_GBIF,
                                        OBIS_occurrences = Acropora_OBIS,
                                        source = "GBIF_&_OBIS")

# Computing abundance grids

Acropora_abundance = oceanAbundGrid(occurrences = Acropora_Total_Checked, cell_size=1)

oceanMaps(Acropora_abundance, logScale=T)

Acropora_richness = oceanDiversity(occurrences = Acropora_Total_Checked, 
                                    diversity_metric = "richness")

oceanMaps(Acropora_richness, logScale=T)



# Example Bryozoa

setwd("C:/Users/monteroserra/Desktop/rOceans Google Summer Code/Example Bryozoa")

Bryozoa_OBIS = read.csv2("Bryozoa_OBIS.csv", stringsAsFactors = F, sep=",")


Bryozoa_OBIS_checked = oceanDataCheck(OBIS_occurrences = Bryozoa_OBIS, source = "OBIS")


Bryozoa_abundance = oceanAbundGrid(occurrences = Bryozoa_OBIS_checked, cell_size=10)
oceanMaps(Bryozoa_abundance, logScale=T, low_color = "steelblue")

Bryozoa_richness = oceanDiversity(occurrences = Bryozoa_OBIS_checked)
oceanMaps(Bryozoa_richness, logScale=T, low_color = "steelblue")


