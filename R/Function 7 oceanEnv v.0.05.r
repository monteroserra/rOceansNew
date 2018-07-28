#' Function 7 oceanEnv: downloads and links environmental variables to marine biodiversity metrics
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
                     log_scale=F, 
                     colors = c("#1874CD60","#4D4DFF60","#23238E60",
                              "#68838B60","#00CD0060", "#CDB38B60")) {


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
           bg=colors[i], col=colors[i])
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

