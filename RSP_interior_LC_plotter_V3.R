# Program that creates interior luminosity light curves using data from Radial Stellar Pulsations
# V1.1: Plots interior luminosity curves for different RSP convection sets on separate canvases, instead of 4 subplots on one canvas
# V2: This version calculates interior absolute magnitudes at different filters and allows the user to create animations
# V3: Changed the method for determining BC.lambda from spatial interpolation to finding the closest grid value using pythagoras

library(av) # to create mp4 file
library(data.table) # for fread()
library(dplyr) # for functions %>% and mutate()
library(gstat) # for inverse-distance weighting
library(pracma)
library(utils)

### CONTROLS ###
cwd <- getwd()
convection.sets <- c('A','B','C','D') # Also the loop vector for our plots
model_nums <- c('1')

RSP.nz <- 300 # number of zones

BC.table.path <- '/home/nick/mesa-r15140/colors/data/lcb98cor.dat'
BC.table <- fread(BC.table.path)
setnames(BC.table,"#Teff","Teff")
#BC.table <- fread(BC.table.path, skip=1)
#colnames(BC.table) <- c("Teff", "logg", "M_div_H", "U", "B", "V", "R", "I", "J", "H", "K", "L", "Lprime", "M")

### END OF CONTROLS ###

#rescale = function(x,y){(x - min(y))/(max(y)-min(y))}

min_max_norm <- function(x) {
  return(((x - min(x)) / (max(x) - min(x))))
}

phytagoras <- function(a,b){
  return(a**2 + b**2)
}

lambdas <- c('U','B','V','R','I','J','H','K','L','Lprime','M')

for (model_num in model_nums){
  cat(paste("Creating interior LC plots for Model",model_num,"\n",sep=" "))
  
  plot.dir <- file.path("plots",paste0("Model_",model_num),"interior_light_curves",lambdas[3])
  ifelse(!dir.exists(plot.dir), dir.create(plot.dir, recursive = TRUE), FALSE)
  
  for (zone_num in 1:RSP.nz){
    
    padded_zone_num <- formatC(zone_num, width = 4, flag = 0)
    
    png(file.path(plot.dir,
                  paste("Light_Curve_","Zone_",padded_zone_num,".png",sep="")), 
        width = 1200, height = 960)
    par(mfrow = c(2, 2))  # Set up a 2 x 2 plotting space
    par(mar=c(5,4,4,2) + 2) # Control space between subplots
    
    for (set in convection.sets){
      
      label <- paste0(set, model_num)
      
      log.directory.after.FA <- paste0("LOGS_",label,"a")
      history.path.after.FA <- file.path(cwd,log.directory.after.FA,"history.data")
      
      DF.after.FA <- fread(history.path.after.FA, header=TRUE, skip=5)
      prof.idx <- fread(file.path(log.directory.after.FA, 'profiles.index'), skip=1, 
                             col.names=c('mdl_num', 'priority', 'prof_num'))
      
      # prof_num <- 1:length(prof.idx$prof_num)
      # prof.path <- file.path(cwd,log.directory.after.FA, paste0('profile', prof_num, '.data'))
      # DF.profile <- lapply(prof.path, function(x) fread(x, header=TRUE, skip = 5))
      
      phases <- DF.after.FA$rsp_phase
      
      ### LOOK INTO PROFILE FILES ###
      luminosities <- c()
      mag.lambdas <- c()
      
      cat(paste("Reading profile files for zone:", zone_num,"and set:", set, "\n"))
      pb = txtProgressBar(min = 0, max = length(prof.idx$prof_num), style = 3)
      for (prof_num in prof.idx$prof_num) {
        prof.path <- file.path(cwd,log.directory.after.FA, paste0('profile', prof_num, '.data'))
        if (!file.exists(prof.path)) next
        #print(prof.path)
        DF.profile <- fread(prof.path, header=TRUE, skip=5)
        temperatures <- DF.profile$temperature
        loggs <- DF.profile$log_g
        
        luminosity <- DF.profile$luminosity[zone_num]
        luminosities <- c(luminosities, luminosity)

        input_Teff <- temperatures[zone_num]
        input_log_g <- loggs[zone_num]
        
        # rescale bolometric correction table
        BC.table %>%
          mutate(Teff.scale = rescale(Teff,Teff),
                 logg.scale = rescale(logg,logg)) -> BC.table
        
        BC.lambda <- BC.table[which.min(raster::pointDistance(cbind(rescale(input_Teff,BC.table$Teff),rescale(input_log_g,BC.table$logg)),BC.table[,c("Teff.scale","logg.scale")], FALSE)),"V"]
        V
        
        mag.bol <- 4.74-(2.5*log10(luminosity))
        mag.lambda <- mag.bol - BC.lambda
        mag.lambdas <- c(mag.lambdas, mag.lambda)
      
        setTxtProgressBar(pb, prof_num)
        
      }
      
      close(pb)
      
      # plot.table <- data.frame(phases, mag.lambdas)
      # o <- order(phases)
      # 
      # with(plot.table, plot(x=phases[o], y=mag.lambdas[o],
      #                main=paste("Set",set,"Light Curve",sep=" "),
      #                type="l", pch=3, lwd = 6, col="purple", xlab=expression("Phase"),
      #                ylab=expression("Absolute Magnitude (mags)"), cex.main=1.60,
      #                cex.lab=1.80, cex.axis=1.60))
      
     plot(x=phases, y=mag.lambdas, main=paste("Set",set,"Light Curve",sep=" "), pch=3, col="purple",
         xlab=expression("Phase"), ylab=expression("Absolute Magnitude (mags)"), cex.main=1.60, cex.lab=1.50, ceb.axis=1.80)
     
    }
    
    mtext(paste("Classical Cepheid","Model",model_num,"Zone",zone_num,lambdas[3],"Band LCs",sep=" "), side = 3, line = -2, cex = 1.90, outer = TRUE)
    dev.off()
    
  }
  setwd(plot.dir)
  av::av_encode_video(list.files(, '*.png'), framerate = 15,
                     output = paste0('CEP_MODEL_', model_num,'_INTERIOR_',lambdas[3],'_BAND_LCs.mp4'))
  #If this does not work, type in terminal ffmpeg -framerate 15 -i Light_Curve_Zone_%04d.png CEP_MODEL_1_INTERIOR_LCs.mp4
  setwd(cwd)
}




