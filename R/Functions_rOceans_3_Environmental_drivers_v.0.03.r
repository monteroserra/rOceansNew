#' ROxygen2 block for function #3.1
#' Function 3.1 oceanEnv: downloads and links environmental variables to marine biodiversity metrics
#'
#' @description downloads and links environmental variables to marine biodiversity metrics
#'
#' @param biod_grid raster layer of any biodiversity metric (abundance, richness, simpson diversity,...)
#' @param biodiv_metric the name of the biodiversity metric
#' @param env_parameters environmental parameters to be downloaded (see example for a complete parameters list)
#' @param plot whether relationships between enviromental variables and biodiversity metrics should be plotted
#' @param log_scale if biodiversity metrics should be log-scaled for a better visualization
#' @return Object of class data.frame with biodiversity metrics, coordinates and environemntal variables
#' @details This function allows accessing and exploring environmental parameters and their relationship
#' with marine biodiversity patterns.
#' @export

oceanEnv = function (biod_grid, 
                     biodiv_metric = "abundance", 
                     env_parameters=c("BO2_tempmean_ss","BO2_tempmax_ss","BO_sstrange",
                                      "BO_bathymean","BO_chlomean", "BO_salinity"), 
                     plot=F, 
                     log_scale=F) {


# getting envrionmental parameters at present conditions
  

if(length(env_parameters) == 1) {

    env_parameters_grid <- sdmpredictors::load_layers(env_parameters) 
  
    }
  
  if(length(env_parameters) > 1 ) {
    
    env_parameters_grid = stack()
    
    for(i in 1: length(env_parameters)){

    env_parameters_grid <- stack(env_parameters_grid, sdmpredictors::load_layers(env_parameters[i]))
  
    }
  }

# Getting centroids from the grid of biodiversity metrics

df = as(as(biod_grid, "SpatialPointsDataFrame"), "data.frame") 

colnames(df) = c(biodiv_metric, "x", "y")

xy = SpatialPoints(df[,c("x","y")])
projection(xy) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")


if(length(env_parameters) == 1) {
  
  df$env_parameters = extract(env_parameters_grid, xy)
  
}

if(length(env_parameters) > 1 ) {
  
  for(i in 1: length(env_parameters)){
    
    df[,3+i] = data.frame(extract(env_parameters_grid[[i]], xy))
  
  }
}

colnames(df) = c(biodiv_metric, "x", "y", env_parameters)


# ploting biodiversity metric ~ environmental parameters (only for single metric so far)

if (plot) {

par(mfrow=c(2,3))
  
if (log_scale) {

  for(i in 1:length(env_parameters)){
      plot(log(df[,1]) ~ df[,i+3], data=df, pch=21, 
           xlab = paste(env_parameters[i]), ylab= paste(biodiv_metric), 
           bg="#104E8B70", col="#104E8B50")
  }
} else {
  
for(i in 1:length(env_parameters)){
       plot(df[,1] ~ df[,i+3], data=df, pch=21, 
       xlab = paste(env_parameters[i]), ylab= paste(biodiv_metric), 
       bg="#104E8B70", col="#104E8B50")
}

}
par(mfrow=c(1,1))
}
  return(df)

}


#' ROxygen2 block for function #3.2
#' Function 3.2 oceanFuture: downloads process future trends in climatic and other variables
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
  
  
  future_layer_code <- get_future_layers(current_layer_codes =env_parameter,
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
    
    factor = cell_size / res(defined_grid)[1]
    defined_grid <- disaggregate(defined_grid, fact=factor)
  }
  
  if (map) {
    
    plot(defined_grid, col=cols)
  }
  
  return(defined_grid)
  
}






