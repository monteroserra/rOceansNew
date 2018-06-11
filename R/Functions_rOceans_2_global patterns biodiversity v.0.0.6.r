#' Funcion 2.1 oceanDataCheck merges, checks and filters big datasets from OBIS and GBIF
#' Funcion 2.2 oceanAbundGrid Compute spatial patterns of abundance of occurrences at multiple scales
#' Funcion 2.3 oceanDiversity Compute spatial patterns of species richness, shannon diversity & simpsons diversity 

library(raster) # Spatial analysis
library(rgeos) # Spatial analysis
library(sp) # Spatial analysis
library(geosphere) # Spatial analysis

require(RColorBrewer)
require(maps)

library(roxygen2) # for Documentation https://www.for.gov.bc.ca/ftp/hts/external/!publish/yuan/LEFI_workshop2017/FLNRO_LiDAR_Workshop/lidar_software_package/r_studio/resources/roxygen_help.html

library(dbplyr)
library(dplyr)


#' Function 2.1 oceanDataCheck: merges, checks and filters big datasets from OBIS and GBIF

# ROxygen2 block for function #2.1

#' merges, checks and filters big datasets from OBIS and GBIF
#'
#'
#' @description merges, checks and filters big datasets from OBIS and GBIF
#'
#'
#' @param source source of occurrences, either "OBIS", "GBIF, or "GBIF_&_OBIS"
#' 
#' @param GBIF_occurrences a data.frame of downloaded occurrences from GBIF
#' 
#' @param OBIS_occurrences a data.frame of downloaded occurrences from ONIS
#' 
#' @param remove_duplicates indicates if checking and removing duplicates
#' 
#' @param remove_NAs indicates if checking and removing NAs
#' 
#' @param remove_NAs indicates if checking and removing land based occurrences  (non marine)
#' 
#' @return Object of class data.frame with filtered occurrences
#'
#'
#' @details The function merges, checks, and filters occurrences data downloaded from OBIS and GBIF
#' 
#'

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
    
    GBIF_occurrences = GBIF_occurrences[,c("scientificname", "species",
                                           "decimallongitude", "decimallatitude",
                                           "catalognumber",  "institutioncode")]
    
    colnames(GBIF_occurrences) <- c("scientificName", "species","decimalLongitude", 
                                    "decimalLatitude","catalogNumber",  "institutionCode")
    
    OBIS_occurrences = OBIS_occurrences[,c("scientificName", "species","decimalLongitude", "decimalLatitude",
                                           "catalogNumber",  "institutionCode")]
    
    
    occurrences = rbind(GBIF_occurrences,OBIS_occurrences)
    
    occurrences$decimalLongitude = as.numeric(occurrences$decimalLongitude)
    occurrences$decimalLatitude = as.numeric(occurrences$decimalLatitude)
    
    
  }
  
  
  # Checking and deleting duplicates
  
  if(remove_duplicates){
    
    occurrences_filtered <- distinct(occurrences, scientificName, decimalLongitude,
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
    
    URL <- "http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/physical/ne_110m_ocean.zip"
    fil <- basename(URL)
    if (!file.exists(fil)) download.file(URL, fil)
    fils <- unzip(fil)
    oceans <- rgdal::readOGR(grep("shp$", fils, value=TRUE), "ne_110m_ocean",
                             stringsAsFactors=FALSE, verbose=FALSE)
    
    spatial_occ <- SpatialPoints(coords = occurrences[,c("decimalLongitude","decimalLatitude")])
    projection(spatial_occ)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    
    land_test = over(spatial_occ, oceans)[,1]
    
    occurrences2 = occurrences[!is.na(land_test),]
    
    cat(paste(nrow(occurrences)-nrow(occurrences2),"land data points were deleted"),sep="\n")
  
    occurrences = occurrences2
  }
  
  return (occurrences)
  
}


#' Function 2.2 oceanAbundGrid: compute spatial patterns of abundance of occurrences at 
#' multiple scales

# ROxygen2 block for function #2.1

#' compute spatial patterns of abundance from larga occurrences datasets at multiple scales
#'
#' @author I. Montero-Serra,  E. Aspillaga, V. Barve, N Barve & K. Kaplan,
#'
#' @description computes global grid of abundance of occurrences at multiple scales
#'
#'
#' @param occurrences dataset  of occurrences downloaded from OBIS or GBIF
#' 
#' @param from source data, it can be "OBIS" or"GBIF" 
#' 
#' @param latitude_name variable name indicating latitude in the occurrences dataset
#' 
#' @param longitude_name variable name indicating longitude in the occurrences dataset
#' 
#' @param min_long minimum longitude of the analysis
#' 
#' @param max_long maximum longitude of the analysis
#' 
#' @param max_long minimum latitute of the analysis
#' 
#' @param max_lat maximum latitude of the analysis
#' 
#' @param cell_size size of the grid cells in degrees (?)
#' 
#' @param extent spatial extent of the analysis:"global", "Mediterranean"
#' 
#' @return The function returns a an object of class raster with abundance per cell
#'
#'
#' @details The function creates downloads occurrences data from the specific sources
#' and stores them in a data frame.
#'
#

## Add the latitudes and longitudes between which each observation is located
## You can substitute any number of breaks you want. Or, a vector of fixed cutpoints
## LATgrid and LONgrid are going to be factors. With ugly level names.

# Data is an OBIS downloaded dataset with scientificName, LAt, Long, and "worms_id"


oceanAbundGrid = function (occurrences,  
                           lat_name = "decimalLatitude", 
                           long_name = "decimalLongitude",
                           extent = "global",cell_size = 5, 
                           min_long = -180, max_long = 180,
                           min_lat = -90,   max_lat = 90) {
  
  
occurrences = na.omit(occurrences)


#Make a raster

if (extent == "global") {
grid <- raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90,
               nrows=180/cell_size, ncols=360/cell_size)

}

if (extent == "manual") {

grid <- raster(xmn = min_long, xmx = max_long , ymn = min_lat, ymx = max_lat,
                nrows=(max_lat-min_lat)/cell_size, ncols=(max_long-min_long)/cell_size)

}

if (extent == "Mediterranean") {

grid <- raster(xmn = -7, xmx = 37, ymn = 30, ymx = 50,
                      nrows=(20)/cell_size, ncols=(45)/cell_size)
}


abundance_grid = rasterize(occurrences[,c(long_name,lat_name)], y=grid, fun='count')


return(abundance_grid)

}

#' Function 2.3 oceanDiversity compute spatial patterns of species diversity from occurrences

# ROxygen2 block for function #2.3

#' compute spatial patterns of species diversity from large occurrences datasets 
#'
#'
#' @description computes global grid of abundance of occurrences at multiple scales
#'
#' @param occurrences dataset  of occurrences downloaded from OBIS or GBIF
#'
#' @param lat_name variable name indicating latitude in the occurrences dataset
#' 
#' @param long_name variable name indicating longitude in the occurrences dataset
#' 
#' @param min_long minimum longitude of the analysis
#' 
#' @param max_long maximum longitude of the analysis
#' 
#' @param max_long minimum latitute of the analysis
#' 
#' @param max_lat maximum latitude of the analysis
#' 
#' @param cell_size size of the grid cells in degrees (?)
#' 
#' @param diversity_metric diversity metric to be computes: richness, shannon or simpson indeces
#' 
#' @param extent spatial extent of the analysis:"global", "Mediterranean"
#' 
#' @return The function returns a an object of class raster with abundance per cell
#'
#'
#' @details The function creates downloads occurrences data from the specific sources
#' and stores them in a data frame.
#'


# Function 2.3 Computing species richness grid at different spatial scales

#Packages
library(dplyr)
library(vegan)


oceanDiversity = function (occurrences, species_name = "scientificName", 
                              lat_name = "decimalLatitude", 
                              long_name = "decimalLongitude",
                              extent = "global",cell_size = 5, 
                              min_long = -180, max_long = 180,
                              min_lat = -90,   max_lat = 90, 
                              diversity_metric="richness") {
  occurrences = Acropora_Total_Checked

#for richness matrices, we can make a loop to split the database for each species and make the same computation
  
 
data = occurrences[, c(species_name,lat_name,long_name)]

species = unique((data[,species_name]))


data$LATgrid<-cut(data$decimalLatitude,breaks=(max_long-min_long)/cell_size,include.lowest=T);

data$LONgrid<-cut(data$decimalLongitude,breaks=(max_lat-min_lat)/cell_size,include.lowest=T);

## Create a single factor that gives the lat,long of each observation.
data$IDgrid<-with(data,interaction(LATgrid,LONgrid))

## Now, create another factor based on the above one, with shorter IDs and no empty levels
data$IDNgrid<-factor(data$IDgrid);
levels(data$IDNgrid)<-seq_along(levels(data$IDNgrid));

## If you want total grid-cell count repeated for each observation falling into that grid cell, do this:
data$count<- ave(data$decimalLatitude,data$IDNgrid,FUN=length)

abundancesGrid <- data[!duplicated(data$IDgrid),]

all_cells = abundancesGrid[,c(long_name, lat_name,"IDgrid","IDNgrid",species_name)]
all_cells = arrange(all_cells, IDNgrid)
sps_abund_matrix = all_cells

for(i in 1:length(species)){
  species_id = species[i]
  data2 = data[data$scientificName==species_id,]
  data2$count<- ave(data2$decimalLatitude,data2$IDNgrid,FUN=length)
  data3 <- na.omit(data2)
  data4 <- data3[!duplicated(data3$IDgrid),]
  
  dataMerge = merge(sps_abund_matrix, data4[,c("IDNgrid","count")], by="IDNgrid", all.x=T)
  dataMerge[is.na(dataMerge)] <- 0
  sps_abund_matrix$species_abundance = dataMerge$count
  colnames(sps_abund_matrix)[5+i] <- paste(species_id)
  cat(paste(i,"of",length(species),"species"), sep="\n")
}

# Here we get a presence/absence mx from species abundance mx and compute species richness per cell

if(diversity_metric == "richness"){

sps_presence_matrix = sps_abund_matrix[,-c(1:5)]
sps_presence_matrix[sps_presence_matrix>0] <- 1
richness = rowSums(sps_presence_matrix) #Sum al species-presence columns 

richness_df = data.frame(richness, sps_abund_matrix$decimalLongitude,
                          sps_abund_matrix$decimalLatitude)


grid <- raster(xmn = min_long, xmx = max_long , ymn = min_lat, ymx = max_lat,
               nrows=(max_lat-min_lat)/cell_size, ncols=(max_long-min_long)/cell_size)


diversity_grid_raster = rasterize(richness_df[,c(2,3)], grid,
                                field=richness_df$richness, fun='last', background=NA)
}


# shannon diversity

if(diversity_metric=="shannon"){

div_index= "shannon"

shannon_values = diversity(sps_abund_matrix[,-c(1:5)], index = div_index, MARGIN = 1, base = exp(1))

shannon_df = data.frame(shannon_values, sps_abund_matrix$decimalLongitude,
                          sps_abund_matrix$decimalLatitude)

diversity_grid_raster = rasterize(shannon_df[,c(2,3)], grid,
                                  field=shannon_df$shannon_values, fun='last', background=NA)

}

#simpson diversity

if(diversity_metric=="simpson"){
  
  div_index= "simpson"
  
  simpson_values = diversity(sps_abund_matrix[,-c(1:5)], index = div_index, MARGIN = 1, base = exp(1))
  
  simpson_df = data.frame(simpson_values, sps_abund_matrix$decimalLongitude,
                            sps_abund_matrix$decimalLatitude)
  
  diversity_grid_raster = rasterize(simpson_df[,c(2,3)], grid,
                                    field=simpson_df$simpson_values, fun='last', background=NA)
  
  
}


return(diversity_grid_raster)

}




### this function is ncluded in set #1 for few species and in set #2 for big data and global analyses. 


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








