#' Function 6 oceanHotspots: classifies sites in hotspots and lowspots of biodiversity
#'
#' @description classify cells within a spatial grid according to biodiversity levels 
#' @param biodiversity_grid raster layer with any biodiversity metric previously computed with oceanDiversity
#' @param hotspot_map whether we want the resulting hotspots and lowspots to be mapped. 
#' @return Object of class raster with three levels of biodiversity (1 = lowspot; 2 = mid diversity; 3 = highdiversity)
#' @details This function allows classifiying sites into high, mid, and low diversity given any biodiversity metric
#' @export


oceanHotspots = function (biodiversity_grid, 
                          only_hotspots=F, 
                          map_hotspots = T,
                          main="", cex.main=0.8) {
  
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
    
  } 
  
  if(map_hotspots) {
    
    oceanMaps(biodiversity_classified, hotspots_map = T)
    
  }
  
  return(biodiversity_classified)
  
}

