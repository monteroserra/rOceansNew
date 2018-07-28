
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
#' @param cell_size size of the grid cells in degrees 
#' 
#' @param diversity_metric diversity metric to be computes: richness, shannon or simpson indeces
#' 
#' @param extent how to the spatial extent of the analysis:"manual", "regions"

#' @param regions spatial extent "global", "Mediterranean", "Caribbean" or "Australia" 
#' 
#' @return The function returns a an object of class raster with abundance per cell
#'
#'
#' @details The function creates downloads occurrences data from the specific sources
#' and stores them in a data frame.
#' 
#' @export



oceanDiversity = function (occurrences, 
                           species_name = "scientificName", 
                           lat_name = "decimalLatitude", 
                           long_name = "decimalLongitude", 
                           extent = "regions", 
                           region="Global", 
                           cell_size = 5, 
                           min_long = -180, 
                           max_long = 180, 
                           min_lat = -90, 
                           max_lat = 90, 
                           print_progress = T) 
{

  data = occurrences[, c(species_name, lat_name, long_name)]
  species = unique(data[, species_name])
  
  if (extent == "manual") {
    grid <- raster(xmn = min_long, xmx = max_long, ymn = min_lat, 
                   ymx = max_lat, nrows = (max_lat - min_lat)/cell_size, 
                   ncols = (max_long - min_long)/cell_size)
  }
  
  if (extent == "regions") {
    
    if (region == "Global") {
      grid <- raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90, 
                     nrows = (180)/cell_size, ncols = (360)/cell_size)
    }
    
    if (region == "Mediterranean") {
      grid <- raster(xmn = -7, xmx = 37, ymn = 30, ymx = 50, 
                     nrows = (20)/cell_size, ncols = (45)/cell_size)
    }
    
    if (region == "Caribbean") {
      grid <- raster(xmn = -100, xmx = -50, ymn = 0, ymx = 35, 
                     nrows = (35)/cell_size, ncols = (50)/cell_size)
    }
    
    if (region == "Australia") {
      grid <- raster(xmn = 100, xmx = 180, ymn = -55, ymx = 0, 
                     nrows = (55)/cell_size, ncols = (80)/cell_size)
    }
    
  }
  
  
  values(grid) <- 0
  
  species_abundance_mx = raster::as.data.frame(grid, xy = T)
  
  #ext = extent(c(min_long, max_long, min_lat, max_lat))
  
  for (i in 1:NROW(species)) {

    species_id = species[i,1]
    
    data_one_species = data[data$scientificName == species_id,]
    
    abundance_one_species = rasterize(data_one_species[, c(long_name, lat_name)], 
                                      y = grid, fun = "count")
    abundance_one_species[is.na(abundance_one_species)] = 0
    extent(abundance_one_species) <- extent(c(-180, 180, 
                                              -90, 90))
    species_abundance_mx[, i + 2] = raster::as.data.frame(abundance_one_species, 
                                                          xy = T)[, 3]
    if (print_progress) {
      cat(paste(i, "of", NROW(species), "species", round((i/NROW(species)) * 
                                                           100), "%"), sep = "\n")
    }
  }
  
  
  results = list()
  results[[1]] = species_abundance_mx
  sps_presence_matrix = species_abundance_mx[, -c(1:3)]
  sps_presence_matrix[sps_presence_matrix > 0] <- 1
  richness = rowSums(sps_presence_matrix)
  richness_df = data.frame(richness, species_abundance_mx$x, 
                           species_abundance_mx$y)
  richness_grid = rasterize(richness_df[, c(2, 3)], grid, field = richness_df$richness, 
                            fun = "last", background = NA)
  richness_grid[richness_grid == 0] <- NA
  results[[2]] = richness_grid
  shannon_values = diversity(species_abundance_mx[, -c(1:3)], 
                             index = "shannon", MARGIN = 1, base = exp(1))
  shannon_df = species_abundance_mx[, c(1:2)]
  shannon_df$shannon_div = shannon_values
  shannon_grid = rasterize(shannon_df[, c(1, 2)], grid, field = shannon_df$shannon_div, 
                           fun = "last", background = NA)
  shannon_grid[shannon_grid == 0] <- NA
  results[[3]] = shannon_grid
  simpson_values = diversity(species_abundance_mx[, -c(1:3)], 
                             index = "simpson", MARGIN = 1, base = exp(1))
  simpson_df = species_abundance_mx[, c(1:2)]
  simpson_df$simpson_div = simpson_values
  simpson_df$simpson_div[rowSums(species_abundance_mx[, -c(1:3)]) == 
                           0] <- 0
  simpson_grid = rasterize(simpson_df[, c(1, 2)], grid, field = simpson_df$simpson_div, 
                           fun = "last", background = NA)
  simpson_grid[simpson_grid == 0] <- NA
  results[[4]] = simpson_grid
  simpson_values = diversity(species_abundance_mx[, -c(1:3)], 
                             index = "simpson", MARGIN = 1, base = exp(1))
  return(results)
}
