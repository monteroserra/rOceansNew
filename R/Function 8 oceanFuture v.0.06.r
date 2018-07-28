#' Function 8 oceanFuture: downloads process future trends in climatic and other variables
#'
#' @description downloads and explore environmental tredns according to IPCC climate change projections
#'
#' @param env_parameter environmental parameters to be explored (see example for a complete parameters list)
#' @param IPCC_scenario set specific IPCC grenhouse emissions' scenarios (RCP85 > RCP45 > RCP25)
#' @param year_proj year for future projections (2100 or 2050)
#' @param manual_extent whether study area geographic limits can be set mannualy
#' @param min_long minimum longitude of the analysis
#' @param max_long maximum longitude of the analysis
#' @param max_long minimum latitute of the analysis
#' @param max_lat maximum latitude of the analysis
#' @param reshape whether we want the new spatial grid to be reshapped with a new cell size
#' @param new_cell_size specific cell_size of the resulting spatial grid
#' @param map whether we want the resulting spatial grid of future trends to be mapped
#' 
#' @return Object of raster with expected changes in environmental variables
#' @details This function allows accessing and exploring expected changes environmental parameters 
#' @export




oceanFuture = function (env_parameter = "BO2_tempmean_ss",
                        IPCC_scenario = "RCP85", 
                        year_proj = 2100,
                        manual_extent = F, 
                        min_long = -7,
                        max_long = 43,
                        min_lat = 30,
                        max_lat = 47, 
                        reshappe = F,
                        new_cell_size=1, 
                        map = T,
                        low_color="steelblue",
                        mid_color="gold",
                        high_color= "firebrick", 
                        col_steps=20) {
  
  # future conditions (only working with a single metric so far)
  
  
  future_layer_code <- sdmpredictors::get_future_layers(current_layer_codes = env_parameter,
                                                        scenario = IPCC_scenario, 
                                                        year = year_proj)$layer_code
  
  
  future_grid <- sdmpredictors::load_layers(future_layer_code) 
  
  
  present_grid <- sdmpredictors::load_layers(env_parameter) 
  
  
  change_grid = future_grid - present_grid
  
  
  
  
  ## color gradient function
  color.gradient <- function(x, colors=c(low_color,mid_color,high_color), colsteps=col_steps) {
    return(colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  
  
  cols = color.gradient(1:20)
  
  
  if (manual_extent) {
    
    study_area <- raster(xmn = min_long, xmx = max_long , ymn = min_lat, ymx = max_lat)
    
  } else {
    
    
    study_area <- raster(xmn = -180, xmx = 180 , ymn = -90, ymx = 90)
  }
  
  
  defined_grid = crop(change_grid, study_area)
  
  if (reshappe) {
    
    factor = new_cell_size / res(defined_grid)[1]
    defined_grid <- disaggregate(defined_grid, fact=factor)
  }
  
  if (map) {
    
    plot(defined_grid, col=cols, cex.main=0.6, 
         main=paste("Projected changes in",env_parameter,"in", year_proj, "( IPCC",IPCC_scenario,")"))
  }
  
  return(defined_grid)
  
}


