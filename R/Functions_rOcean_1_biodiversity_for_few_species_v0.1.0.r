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

  URL <- "http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/physical/ne_110m_ocean.zip"
  fil <- basename(URL)
  if (!file.exists(fil)) download.file(URL, fil)
  fils <- unzip(fil)
  oceans <- readOGR(grep("shp$", fils, value=TRUE), "ne_110m_ocean",
                    stringsAsFactors=FALSE, verbose=FALSE)


spatial_occ <- SpatialPoints(coords = my_occurrences[,c("decimalLongitude","decimalLatitude")])
projection(spatial_occ)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

land_test = over(spatial_occ, oceans)[,1]


warning(paste(nrow(my_occurrences)-nrow(my_occurrences2),"land data points were deleted"))

my_occurrences = my_occurrences2

}

return(my_occurrences)

}

#' Function 1.1 oceanDist: Creates a presence / absence raster from occurrence data
#'
#' @author I. Montero-Serra,  E. Aspillaga, V. Barve, N Barve & K. Kaplan,
#'
#' @description creates a presence absence raster from occurrences data
#'
#' @param occurrences species occurrences dataframe, with at least longitude and latitude columns
#'
#' @param long_name column name for longitudes in the "occurrences" dataframe
#'
#' @param lat_name column name for latitudes in the "occurrences" dataframe
#'
#' @param min_long minimum longitude of the analysis
#'
#' @param max_long maximum longitude of the analysis
#'
#' @param max_long minimum latitute of the analysis
#'
#' @param max_lat maximum latitude of the analysis
#'
#' @param raster_name to include a name to each created raster
#'
#' @param cell_size interger indicating grid cell size
#'
#' @return The result is an object of class raster with presences (1) and absences (0)
#'
#'
#' @details The function creates a presences raster from  occurrences data
#' 
#' @export

presencesRaster <- function (occurrences,
                             long_name = "decimalLongitude",
                             lat_name = "decimalLatitude",
                             extent="global",
                             min_long = -180, max_long = 180,
                             min_lat = -90, max_lat = 90,
                             raster_name="", cell_size=5) {
  
  
  occurrences = na.omit(occurrences)
  species_data_xy = c()

if(ncol(occurrences)>2) {

  species_data_xy = occurrences[,c(long_name,lat_name)]

}


#Creates an empty raster
if (extent == "global"){
  grid <- raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90,
                        nrows=180/cell_size, ncols=360/cell_size)
  projection(grid) <-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

}

if (extent == "Mediterranean"){

  grid <- raster(xmn = 0, xmx = 35, ymn = 30, ymx = 50,
                        nrows=180/cell_size, ncols=360/cell_size)
  projection(grid) <-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

}

if (extent == "manual"){

 grid <- raster(xmn = min_long, xmx = max_long, ymn = min_lat, ymx = max_lat,
                        nrows=180/cell_size, ncols=360/cell_size)
 projection(grid) <-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

}


#Transforms occurrences to spatial presences
spatial_occ <- SpatialPoints(coords = species_data_xy)
projection(spatial_occ)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# set the background cells in the raster to 0
grid[!is.na(grid)] <- 0

#set the cells that contain points to 1
speciesRaster <- rasterize(spatial_occ,grid,field=1)
speciesRaster <- merge(speciesRaster,grid)

#label the raster
names(speciesRaster) <- raster_name
return(speciesRaster)
}



#' Function #1.3 oceanRichness creates a richness layer from occurrence data
#'
#' @author I. Montero-Serra,  E. Aspillaga, V. Barve, N Barve & K. Kaplan,
#'
#' @description creates a richness raster from multispecies occurrences data
#'
#' @param occurrences species occurrences dataframe, with at least longitude and latitude columns
#'
#' @param long_name column name for longitudes in the "occurrences" dataframe
#'
#' @param lat_name column name for latitudes in the "occurrences" dataframe
#'
#' @param species_names column name for species names in the "occurrences" dataframe
#'
#' @param min_long minimum longitude of the analysis
#'
#' @param max_long maximum longitude of the analysis
#'
#' @param max_long minimum latitute of the analysis
#'
#' @param max_lat maximum latitude of the analysis
#'
#' @param cell_size interger indicating grid cell size
#'
#' @param richness_map plot the resulting richness map
#'
#' @return The result is an object of class raster with richness grid cell
#'
#' @details
#' @export
#' 


oceanRichness = function(occurrences, species_name = "species",
                         long_name = "decimalLongitude",
                         lat_name = "decimalLatitude",
                         extent="global",
                         min_long = -180, max_long = 180,
                         min_lat = -90, max_lat = 90,
                         cell_size=5, richness_map=F){

  
  
  multilayers = stack()

  species_names = unique((occurrences[,species_name]))

for(i in 1:length(species_names)){

  single_sp_data = occurrences[occurrences[,species_name]==species_names[i],]
  single_sp_raster = presencesRaster(single_sp_data, raster_name=paste(species_names[i]))
  multilayers = stack(multilayers, single_sp_raster)

}


raster_richness = stackApply(multilayers, 1, fun = sum)

raster_richness[raster_richness==0] <- NA


if(richness_map) {
  plot(raster_richness)
  require(maps)
  map("world",add=T, fill = T,bg="grey10",col="grey10")
}

  return(raster_richness)

}

