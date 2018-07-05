#' Function 1.1 oceanDist: allows getting presences and distribution maps  for a list of species
#' ROxygen2 block for function #1
#' Allows getting distribution data  for a list of species from OBIS & GBIF
#'
#' @author I. Montero-Serra,  E. Aspillaga, V. Barve, N Barve & K. Kaplan,
#'
#' @description allows getting occurrences from OBIS and GBIF for a list of species
#'
#' @param species_names list of species names in latin
#'
#' @param data_source source data, it can be "OBIS","GBIF" or both "OBIS_GBIF"
#'
#' @param remove_duplicates if True, checks and removes duplicates
#'
#' @param remove_land_dots if True, checks and removes land based occurrences
#'
#' @param minimum_data treshold of minimum required data for including species in the analysis,
#' checks and deletes any species for which less than "minim_data" were found
#'
#' @return The result is an object of class dataframe that contains species names
#' and coordiantes
#'
#'
#' @details The function creates downloads occurrences data from the specific sources
#' and stores them in a data frame.
#' @export


oceanDivDat = function(species_names,
                       data_source = "OBIS",
                       remove_duplicates=T,
                       remove_land_dots = T,
                       minimum_data=5) {

#species_names is a list with the specific scientific species names in latin


if (data_source == "GBIF"){

my_occurrences = data.frame()
my_occurrences <- plyr::ldply(species_names, function(i) {
    GBIF_info <- occ_search(scientificName = i)
    return(GBIF_info$data[,c("species","decimalLongitude","decimalLatitude")])
})

}

if (data_source == "OBIS"){

my_occurrences = data.frame()
my_occurrences <- plyr::ldply(species_names, function(i) {
    OBIS_info <- occurrence(scientificname = i)
    return(OBIS_info[,c("species","decimalLongitude","decimalLatitude")])
})

}

if (data_source == "OBIS_&_GBIF"){

my_occurrences_GBIF= data.frame()
my_occurrences_GBIF <- plyr::ldply(species_names, function(i) {
    GBIF_info <- occ_search(scientificName = i)
    GBIF_dat = GBIF_info$data[,c("species","decimalLongitude","decimalLatitude",
                      "catalogNumber",  "institutionCode")]
    return(GBIF_dat)

})

my_occurrences_OBIS= data.frame()
my_occurrences_OBIS <- plyr::ldply(species_names, function(i) {
      OBIS_info <- occurrence(scientificname = i)
      OBIS_dat = OBIS_info[,c("species","decimalLongitude","decimalLatitude",
                          "catalogNumber",  "institutionCode")]
      return(OBIS_dat)
})

my_occurrences <- na.omit(rbind(my_occurrences_GBIF,my_occurrences_OBIS))

}

# This sections checks and filters any duplicate occurrence (very important when merging
# OBIS and GBIF datasets)

if(remove_duplicates)  {

my_occurrences2 <- distinct(my_occurrences, species, decimalLongitude,
                               decimalLatitude,institutionCode) # Deleting duplicates

warning(paste(nrow(my_occurrences)-nrow(my_occurrences2),"duplicate data points were deleted"))

my_occurrences = my_occurrences2

}

 # Here we check and delete species for which less than "minimum data" were vailable
  filtered <- table(my_occurrences$species)
  sortout <- names(filtered[filtered <= minimum_data])
  filtered <- filtered[filtered > minimum_data]

  my_occurrences <- droplevels(subset(my_occurrences, my_occurrences$species %in% as.character(names(filtered))))


  warning(paste(length(sortout),"species deleted due to a lack comprehensive data (less than",
          paste(minimum_data),"occurrences)"))


# This setcion checks and filters data on land

if(remove_land_dots)  {

data(oceans)

spatial_occ <- SpatialPoints(coords = my_occurrences[,c("decimalLongitude","decimalLatitude")])
projection(spatial_occ)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

land_test = over(spatial_occ, oceans)[,1]


warning(paste(nrow(my_occurrences)-nrow(my_occurrences2),"land data points were deleted"))

my_occurrences = my_occurrences2

}

return(my_occurrences)

}

