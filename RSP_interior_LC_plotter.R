# Program that creates interior luminosity light curves using data from Radial Stellar Pulsations

library(av) # for animation
library(data.table) # for fread()

### CONTROLS ###

cwd <- getwd()
convection.sets <- c('A','B','C','D') # Also the loop vector for our plots
model_nums <- c('1')

RSP.nz <- 7 # number of zones

frame.rate <- 1 # frame rate of animation

BC.table.path <- '/home/nick/mesa-r15140/colors/data/lcb98cor.dat' # MESA's bolometric correction table
BC.table <- fread(BC.table.path)
setnames(BC.table,"#Teff","Teff")

### END OF CONTROLS ###

light_curve_plotter <- function(filepath, title, phases, values) {
  png(filepath, width = 750, height = 600)
  par(mar=c(5,5,4,1)+.1)
  plot(x=phases, y=values, main=title, pch=3, col="purple",
       xlab=expression("Phase"), ylab=expression("Absolute Magnitude (mags)") , cex.main=2.0, cex.lab=1.50, ceb.axis=1.80)
  dev.off()
}

lambdas <- c('U','B','V','R','I','J','H','K','L','Lprime','M')

for (model_num in model_nums){
  cat(paste("Creating interior LC plots for Model",model_num,"\n",sep=" "))
  
  plot.dir <- file.path("plots",paste0("Model_",model_num),"interior_light_curves","luminosity")
  ifelse(!dir.exists(plot.dir), dir.create(plot.dir, recursive = TRUE), FALSE)
  
  for (zone_num in 1:RSP.nz){
    
    padded_zone_num <- formatC(zone_num, width = 4, flag = 0)
    
    png(file.path(plot.dir,
                  paste("Light_Curve_","Zone_",padded_zone_num,".png",sep="")), width = 1200, height = 960)
    par(mfrow = c(2, 2))  # Set up a 2 x 2 plotting space
    par(mar=c(5,4,4,2) + 2) # Control space between subplots
    
    for (set in convection.sets){
      label <- paste0(set, model_num)
      
      log.directory.after.FA <- paste0("LOGS_",label,"a")
      history.path.after.FA <- file.path(cwd,log.directory.after.FA,"history.data")
      
      DF.after.FA <- fread(history.path.after.FA, header=TRUE, skip=5)
      prof.idx <- fread(file.path(log.directory.after.FA, 'profiles.index'), skip=1, 
                             col.names=c('mdl_num', 'priority', 'prof_num'))
      
      phases <- DF.after.FA$rsp_phase
      
      ### LOOK INTO PROFILE FILES ###
      luminosities <- c()
      
      cat(paste("Reading profile files for zone:", zone_num,"and set:", set, "\n"))
      pb = txtProgressBar(min = 0, max = length(prof.idx$prof_num), style = 3)
      for (prof_num in prof.idx$prof_num) {
        prof.path <- file.path(cwd,log.directory.after.FA, paste0('profile', prof_num, '.data'))
        if (!file.exists(prof.path)) next
        #print(prof.path)
        DF.profile <- fread(prof.path, header=TRUE, skip=5)
        
        luminosity <- DF.profile$luminosity[zone_num]
        luminosities <- c(luminosities, luminosity)
        
        setTxtProgressBar(pb, prof_num)
      }
      close(pb)
      
      plot(x=phases, y=luminosities, main=paste("Set",set,"Light Curve",sep=" "), pch=3, col="purple",
           xlab=expression("Phase"), ylab=expression("Luminosity " (L/L['\u0298'])), cex.main=1.60, cex.lab=1.50, ceb.axis=1.80)
    }
    mtext(paste("Classical Cepheid","Model",model_num,"Zone",zone_num,"LCs",sep=" "), side = 3, line = -2, cex = 1.90, outer = TRUE)
    dev.off()
    
  }
  setwd(plot.dir)
  animation.title <- paste0('CEP_MODEL_', model_num,'_INTERIOR_LUMINOSITY_LCs.mp4')
  av::av_encode_video(list.files(, '*.png'), framerate = frame.rate,
                      output = animation.title)
  # If this does not work, type in terminal ffmpeg -framerate 15 -i Light_Curve_Zone_%04d.png CEP_MODEL_1_INTERIOR_LCs.mp4
  cat(paste("Animation titled",animation.title,"saved in",plot.dir,"\n"))
  cat(paste("You may ignore warnings.\n"))
  setwd(cwd)
}




