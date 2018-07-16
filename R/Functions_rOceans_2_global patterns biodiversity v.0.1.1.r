#' ROxygen2 block for function #2.1
#' Function 2.1 oceanDataCheck: merges, checks and filters big datasets from OBIS and GBIF
#' merges, checks and filters big datasets from OBIS and GBIF
#'
#' @description merges, checks and filters big datasets from OBIS and GBIF
#'
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
#' @export

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




#' Function #2.2 oceanMaps function in progress to plot different outputs of the rOcean package
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
#' @export


oceanMaps = function (grid = NULL, 
                      grid_provided=T,
                      logScale = F,  
                      occurrences = NULL, 
                      map_occurrences = F,
                      convex_hull = F,
                      long_name="decimalLongitude", 
                      lat_name = "decimalLatitude",
                      min_long = -180, max_long = 180,
                      min_lat = -90, max_lat = 90,
                      species_name = "species",
                      background_color="grey10", 
                      dot_color="steelblue", 
                      low_color="steelblue", 
                      mid_color="gold",
                      high_color= "firebrick", 
                      col_steps=20, 
                      main="",cex_main=0.8, 
                      hotspots_map=F){

## color gradient function
color.gradient <- function(x, colors=c(low_color,mid_color,high_color), colsteps=col_steps) {
    return(colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  
  
cols = color.gradient(1:20)
  
if(grid_provided){
  
    
    if (logScale){
      
      plot(log(grid), col=cols, main=main, cex.main = cex_main)
      
    } 
    
    else {
      plot(grid, col=cols, main=main, cex.main = cex_main)
    }
    
    maps::map("world",add=T, fill = T,bg=background_color,col=background_color)
    maps::map.axes()
    
  }
  else {
   
    maps::map("world", fill = T,col=background_color)
    maps::map.axes()
    
  }
  
  if(map_occurrences){
    
    
    
    if (length(levels(occurrences$species)) > 1){
      
      palette(color.gradient(1:20))
      
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


if(hotspots_map){
  
  maps::map("world", fill = T,col=background_color)
  maps::map.axes()
  plot(grid, col=c("steelblue", "gold","firebrick"), add=T, legend=F)
  maps::map("world", fill = T,add=T,col="grey20")

  }
  
}


#' Function 2.3 oceanDiversity compute spatial patterns of species diversity from occurrences
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
#' @export

oceanDiversity = function (occurrences, species_name = "scientificName", 
                            lat_name = "decimalLatitude", 
                            long_name = "decimalLongitude",
                            extent = "global",cell_size = 5, 
                            min_long = -180, max_long = 180,
                            min_lat = -90,   max_lat = 90, 
                            print_progress=T) {
  
  data = occurrences[, c(species_name,lat_name,long_name)]
  species = unique((data[,species_name]))
  
# Make a raster
  
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
  
  
  values(grid) <- 0
  
  species_abundance_mx = as.data.frame(grid, xy=T)
  
  ext =  extent(c(min_long,max_long, min_lat, max_lat))
  
for(i in 1:length(species)){

    species_id = species[i]
    
    data_one_species = data[data$scientificName==species_id,]
    
    abundance_one_species = rasterize(data_one_species[,c(long_name,lat_name)], y=grid, fun='count')
    
    abundance_one_species[is.na(abundance_one_species)] = 0
    
    extent(abundance_one_species) <- extent(c(-180,180, -90, 90))
    
    species_abundance_mx[,i+2] = as.data.frame(abundance_one_species, xy=T)[,3]
    
    if(print_progress){
      
      cat(paste(i, "of", length(species), "species", round((i/length(species)) * 100), "%"), sep = "\n")
    
    }
  }

#storing results in a list
results = list()
  
# first element: the species abundance matrix 

results[[1]] = species_abundance_mx

# second element: species richness grid
  
  sps_presence_matrix = species_abundance_mx[,-c(1:3)]
  sps_presence_matrix[sps_presence_matrix>0] <- 1
  richness = rowSums(sps_presence_matrix) #Sum al species-presence columns 
  
  richness_df = data.frame(richness, species_abundance_mx$x,species_abundance_mx$y)
  
  richness_grid = rasterize(richness_df[,c(2,3)], grid, field=richness_df$richness, fun='last', background=NA)
  richness_grid[richness_grid==0] <- NA

  results[[2]] = richness_grid
  
# third element: Shannon diversity

shannon_values = diversity(species_abundance_mx[,-c(1:3)], 
                           index = "shannon", MARGIN = 1, base = exp(1))
  
shannon_df = species_abundance_mx[,c(1:2)]

shannon_df$shannon_div = shannon_values

shannon_grid = rasterize(shannon_df[,c(1,2)], grid,
                         field=shannon_df$shannon_div, 
                         fun='last', background=NA)
shannon_grid[shannon_grid==0] <- NA

results[[3]] = shannon_grid
  
# forth element: Simpson diversity

simpson_values = diversity(species_abundance_mx[,-c(1:3)], 
                           index = "simpson", MARGIN = 1, base = exp(1))

simpson_df = species_abundance_mx[,c(1:2)]

simpson_df$simpson_div = simpson_values

simpson_df$simpson_div[rowSums(species_abundance_mx[,-c(1:3)])==0] <- 0


simpson_grid = rasterize(simpson_df[,c(1,2)], grid,
                         field=simpson_df$simpson_div, 
                         fun='last', background=NA)

simpson_grid[simpson_grid==0] <- NA

results[[4]] = simpson_grid


simpson_values = diversity(species_abundance_mx[,-c(1:3)], 
                           index = "simpson", MARGIN = 1, base = exp(1))


  return(results)
  
}



#' ROxygen2 block for function #2.4
#' 
#' Function 2.4 oceanHotspots: explore hotspots and lowspots of biodiversity
#'
#' @description classify cells within a spatial grid according to biodiversity levels 
#' @param biodiversity_grid raster layer with any biodiversity metric previously computed with oceanDiversity
#' @param hotspot_map whether we want the resulting hotspots and lowspots to be mapped. 
#' @return Object of class raster with three levels of biodiversity (1 = lowspot; 2 = mid diversity; 3 = highdiversity)
#' @details This function allows classifiying sites into high, mid, and low diversity given any biodiversity metric
#' @export

oceanHotspots = function (biodiversity_grid, 
                          hotspot_map=T, 
                          only_hotspots=F, 
                          main="", cex.main=0.8) {
  
  ## Identifying biodiversity hotspot
  
  ## Framework for classififying into low - mid - high diversity sites
  min_diversity = minValue(biodiversity_grid)
  low_treshold = summary(biodiversity_grid)[2] # third quartile (all bellow is considered lowspot)
  high_treshold = summary(biodiversity_grid)[4] # third quartile (all above is considered hotspot)
  max_diversity = maxValue(biodiversity_grid)
  
  m = c(min_diversity,low_treshold,1, low_treshold, high_treshold, 2,
        high_treshold, max_diversity, 3)
  
  rclmat = matrix(m, ncol=3, byrow=TRUE)
  biodiversity_classified = reclassify(biodiversity_grid, rclmat)

if(only_hotspots){
    
    biodiversity_classified[biodiversity_classified<3] <- NA
    
  if(hotspot_map){
      maps::map("world", fill = T,col="grey20",main=main, cex.main=cex.main)
      maps::map.axes()
      plot(biodiversity_classified, col=c("firebrick"), add=T, legend=F)
      maps::map("world", fill = T,add=T,col="grey20")
    
      }
    
  } else {
  
  if(hotspot_map){
    maps::map("world", fill = T,col="grey20", main=main,cex.main=cex.main)
    maps::map.axes()
    plot(biodiversity_classified, col=c("steelblue", "gold","firebrick"), add=T, legend=F)
    maps::map("world", fill = T,add=T,col="grey20")
    
  }
}
 
  return(biodiversity_classified)
  
}







