#' Function #4 oceanMaps function in progress to plot different outputs of the rOcean package
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
#' @param hotspots_map whether a hotspots layer is provided 
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
                      background_color="grey60", 
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
  
  
  cols = color.gradient(1:col_steps)
  
if(grid_provided){
    
    if (logScale){
      
      plot(log(grid), col=cols, main=main, cex.main = cex_main)
      
    } else {
      
      plot(grid, col=cols, main=main, cex.main = cex_main)
    }
    
    maps::map("world",add=T, fill = T,bg=background_color,col=background_color)
    
} else {
    
  maps::map("world",add=F, fill = T,bg="white",col=background_color)
  maps::map.axes()
  }
  
  if(map_occurrences){
    
    species = unique(occurrences$species)
    
    if (length(species) > 1){
      
      palette(color.gradient(1:length(species)))
      
    }
    
    if (convex_hull) {
      
      for(i in 1:length(occurrences$species)){
        xy =  occurrences[occurrences$species==species[i],c(long_name,lat_name)]
        coordinates(xy) <- ~decimalLongitude + decimalLatitude
        Convex_hull <- rgeos::gConvexHull(xy)
        proj4string(Convex_hull) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
        plot(Convex_hull, add=T, lwd=2, col=cols)
        
      }
      
      points(na.omit(occurrences[,c(long_name, lat_name)]), pch=21, bg=factor(occurrences$species))
      
    }
    
    
  }
  
  
  if(hotspots_map){
    
    plot(grid, col=c("steelblue", "gold","firebrick"), add=F, legend=F)
    maps::map("world", fill = T,add=T,col=background_color)
    maps::map.axes()
    
  }
  
}

