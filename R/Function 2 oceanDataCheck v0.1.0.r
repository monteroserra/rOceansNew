#' Function 2 oceanDataCheck: this function allow to merge, check and 
#' filter big data from OBIS and GBIF 
#' 
#' @description merges, checks and filters big datasets from OBIS and GBIF
#'
#' @param source source of occurrences, either "OBIS", "GBIF, or "GBIF_&_OBIS"
#' @param GBIF_occurrences a data.frame of downloaded occurrences from GBIF
#' @param OBIS_occurrences a data.frame of downloaded occurrences from ONIS
#' @param remove_duplicates indicates if checking and removing duplicates
#' @param remove_NAs indicates if checking and removing NAs
#' @param remove_NAs indicates if checking and removing land based occurrences  (non marine)
#' @return Object of class data.frame with filtered occurrences
#' @details The function merges, checks, and filters occurrences data downloaded from OBIS and GBIF
#' 
#' @export

oceanDataCheck = function (source="OBIS", 
                           GBIF_occurrences,
                           OBIS_occurrences,
                           remove_duplicates = T,
                           remove_NAs = T,
                           remove_land_occ = T) {
  
  if(source=="GBIF"){
    
    occurrences = GBIF_occurrences[,c("scientificname", "species","decimallongitude", "decimallatitude",
                                      "catalognumber",  "institutioncode")]
    
    colnames(occurrences) = c("scientificName", "species","decimalLongitude", "decimalLatitude",
                              "catalogNumber",  "institutionCode")
    
    occurrences$decimalLongitude = as.numeric(occurrences$decimalLongitude)
    occurrences$decimalLatitude = as.numeric(occurrences$decimalLatitude)
    
  }
  
  if(source=="OBIS"){
    
    occurrences = OBIS_occurrences[,c("scientificName", "species","decimalLongitude", "decimalLatitude",
                                      "catalogNumber",  "institutionCode")]
    
    occurrences$decimalLongitude = as.numeric(occurrences$decimalLongitude)
    occurrences$decimalLatitude = as.numeric(occurrences$decimalLatitude)
  }
  
  if(source=="GBIF_&_OBIS"){
    
    GBIF_occurrences_filtered = GBIF_occurrences[GBIF_occurrences$taxonrank == "SPECIES",]
    
    cat(paste(nrow(GBIF_occurrences)-nrow(GBIF_occurrences_filtered),"occurrences were deleted from GBIF data due to a lack of standard species level information"),sep="\n")
    
    
    GBIF_occurrences = GBIF_occurrences_filtered[,c("scientificname", "species",
                                                    "decimallongitude", "decimallatitude",
                                                    "catalognumber",  "institutioncode")]
    
    colnames(GBIF_occurrences) <- c("scientificName", "species","decimalLongitude", 
                                    "decimalLatitude","catalogNumber",  "institutionCode")
    
    OBIS_occurrences_filtered = OBIS_occurrences[!is.na(OBIS_occurrences$worms_id),]
    
    cat(paste(nrow(OBIS_occurrences)-nrow(OBIS_occurrences_filtered),"occurrences were deleted from OBIS data due to of standard species level information"),sep="\n")
    
    
    OBIS_occurrences = OBIS_occurrences_filtered[,c("scientificName", "species","decimalLongitude", "decimalLatitude",
                                                    "catalogNumber",  "institutionCode")]
    
    
    occurrences = rbind(GBIF_occurrences,OBIS_occurrences)
    
    occurrences$decimalLongitude = as.numeric(occurrences$decimalLongitude)
    occurrences$decimalLatitude = as.numeric(occurrences$decimalLatitude)
    
    
  }
  
  
  # Checking and deleting duplicates
  
  if(remove_duplicates){
    
    occurrences_filtered <- dplyr::distinct(occurrences, scientificName, decimalLongitude,
                                            decimalLatitude,institutionCode,institutionCode) # Deleting duplicates
    
    cat(paste(nrow(occurrences)-nrow(occurrences_filtered),"duplicate data points were deleted"),sep="\n")
    
    occurrences = occurrences_filtered
    
  }
  
  # deleting incomplete information occurrences 
  
  if(remove_NAs)  {
    
    occurrences_filtered = na.omit(occurrences[,c("scientificName","decimalLongitude", "decimalLatitude")])
    
    cat(paste(nrow(occurrences)-nrow(occurrences_filtered),"data points containing NAs were deleted"),sep="\n")
    
    occurrences = occurrences_filtered
    
  }
  
  ### Removing land occurrences 
  
  if(remove_land_occ)  {
    
    data(oceans)
    
    spatial_occ <- SpatialPoints(coords = occurrences[,c("decimalLongitude","decimalLatitude")])
    projection(spatial_occ)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    
    land_test = over(spatial_occ, oceans)[,1]
    
    occurrences2 = occurrences[!is.na(land_test),]
    
    cat(paste(nrow(occurrences)-nrow(occurrences2),"land data points were deleted"),sep="\n")
    
    occurrences = occurrences2
  }
  
  return (occurrences)
  
}
