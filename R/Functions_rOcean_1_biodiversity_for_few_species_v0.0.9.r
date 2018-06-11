# Project: rOceans
# Starting date: 16/05/2018
# Last update: 11/06/2018
# Authors: I. Montero-Serra,  E. Aspillaga, V. Barve, N Barve & K. Kaplan,
# Funding: Google Summer Code 2018
# Script name: rOceans - Marine Biodiversity: Functions to access and process biodiversity data

# Short description of oceanDivData:
# This is the first set of functions of the rOceans package. They take advantage of robis, rgbif and
# The goal is to access marine biodiversity data from different sources format them into spatial objects
# that can then be used for macroegological and
# marine conservation analyses.Due to computational limitations, this first set of functions
# is aimed for analazing only single or few species. For analyzing global patterns including
# thousands of species see second set of functions.

#' Funcion 1.1 oceanDivData Access & filter distribution data for a range of species from OBIS & GBIF
#' Funcion 1.2  presenceRaster Make rasters of presence absence data
#' Function 1.3 oceanRichness: Creates a raster object with species richness per cell
#' Function 1.4 oceanMaps: Map occurrences and species distributions

#'Required packages
### packages to instal

library(robis) # Access to GBIF Data
library(rgbif) # Access to OBIS Data
library(dbplyr)
library(dplyr)
library(raster) # Spatial analysis
library(rgeos) # Spatial analysis
library(sp) # Spatial analysis
library(geosphere) # Spatial analysis

library(roxygen2) # for Documentation


#Function 1.1 oceanDist: allows getting presences and distribution maps  for a list of species

# ROxygen2 block for function #1

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

  corals

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
#Function #1.2 presencesRaster creates presence / absence raster from species occurrence data

# ROxygen2 block for function #1.2


#' Create a presence / absence raster from occurrence data
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



#Function #1.3 oceanRichness creates richness raster from multispecies species occurrence data

# ROxygen2 block for function #1.3


#' Create a presence / absence raster from occurrence data
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

if(richness_map) {
  plot(raster_richness)
  require(maps)
  map("world",add=T, fill = T,bg="grey10",col="grey10")
}

  return(raster_richness)

}



#Function #1.4 oceanMaps function in progress to plot different outputs of the rOcean package


# ROxygen2 block for function #1.4


#' visualize spatial patterns in marine biodiversity
#'
#' @author I. Montero-Serra,  E. Aspillaga, V. Barve, N Barve & K. Kaplan,
#'
#' @description to map different metrics computed using functions of the rOcean package
#'
#' @param grid raster of species presences, richness, abundance or other metrics.
#'
#' @param grid_provided to indicate whether a raster is provided to be mapped
#'
#' @param log_scale if True grid values are plotted in a log scale, useful for abundance representations
#'
#' @param occurrences species occurrences dataframe with at least latitude and longitude columns
#'
#' @param map_occurrences if true, occurrences are mapped for each species
#'
#' @param convex_hull if true, a convex hull is computed and mapped for each species of occurrences
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
#' @return The result is a map with the spatial trends
#'
#' @details



oceanMaps = function (grid,
                      grid_provided=T,
                      logScale = F,
                      occurrences,
                      map_occurrences = F,
                      convex_hull = F,
                      long_name="decimalLongitude",
                      lat_name = "decimalLatitude",
                      min_long = -180, max_long = 180,
                      min_lat = -90, max_lat = 90,
                      species_name = "species",
                      background_color="grey10",
                      dot_color="steelblue"){



  if(grid_provided){

     reds = palette(RColorBrewer::brewer.pal(n = 8, name = "Reds"))[-c(1:3)]

    if (logScale){

      plot(log(grid), col=reds)

    }

    else {
      plot(grid, col=reds)
    }

    maps::map("world",add=T, fill = F,bg="grey40",col="grey40")
    maps::map.axes()

  }
  else {

    maps::map("world", fill = T,col=background_color)
    maps::map.axes()

  }

  if(map_occurrences){



if (length(levels(occurrences$species)) > 1){

      palette(brewer.pal(n = length(species_names), name = "Set2"))

    }

    if (convex_hull) {

      for(i in 1:length(occurrences$species)){
        xy =  occurrences[occurrences$species==levels(as.factor(occurrences$species))[i],c(long_name,lat_name)]
        coordinates(xy) <- ~decimalLongitude + decimalLatitude
        Convex_hull <- rgeos::gConvexHull(xy)
        proj4string(Convex_hull) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
        plot(Convex_hull, add=T, lwd=2, col="grey70")

      }

    points(na.omit(occurrences[,c(long_name, lat_name)]), pch=21, bg=factor(occurrences$species))

  }


  }

}




# Example of usage of four functions for three coral species

corals = c("Corallium rubrum", "Paramuricea clavata", "Eunicella singularis")

coral_occ = oceanDivDat(species_names=corals,data_source="OBIS_&_GBIF")

coral_raster = presencesRaster(occurrences = coral_occ, extent = "Mediterranean",
                              raster_name="Corals", cell_size=5)


coral_richness = oceanRichness(coral_occ)

oceanMaps(occurrences=coral_occ, grid_provided=F, map_occurrences = T, convex_hull=T)

oceanMaps(coral_richness, grid_provided=T)










