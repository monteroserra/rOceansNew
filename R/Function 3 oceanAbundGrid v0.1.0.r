#' Function 3 oceanAbundGrid: compute spatial patterns of abundance of occurrences at 
#' multiple scales
#'
#' @author I. Montero-Serra,  E. Aspillaga, V. Barve, N Barve & K. Kaplan,
#'
#' @description computes global grid of abundance of occurrences at multiple scales
#'
#'
#' @param occurrences dataset  of occurrences downloaded from OBIS or GBIF
#' @param lat_name variable name indicating latitude in the occurrences dataset
#' @param long_name variable name indicating longitude in the occurrences dataset
#' @param min_long minimum longitude of the analysis
#' @param max_long maximum longitude of the analysis
#' @param max_long minimum latitute of the analysis
#' @param max_lat maximum latitude of the analysis
#' @param cell_size size of the grid cells in degrees (?)
#' @param extent how to the spatial extent of the analysis:"manual", "regions"
#' @param regions spatial extent "global", "Mediterranean", "Caribbean" or "Australia" 
#' @return The function returns a an object of class raster with abundance per cell
#'
#' @details The function creates downloads occurrences data from the specific sources
#' and stores them in a data frame.
#'
#' @export
#' 

oceanAbundGrid = function (occurrences, 
                           lat_name = "decimalLatitude", 
                           long_name = "decimalLongitude", 
                           extent = "regions", 
                           region="Global", 
                           cell_size = 5, 
                           min_long = -180, max_long = 180, 
                           min_lat = -90, max_lat = 90, 
                           map_abundance=F) 
{
  
  occurrences = na.omit(occurrences)
  
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
  
  abundance_grid = rasterize(occurrences[, c(long_name, lat_name)], 
                             y = grid, fun = "count")
  
  if(map_abundance){
    
    oceanMaps(abundance_grid, main = "Abundance")
    
  }
  
  return(abundance_grid)
}
