# Program that creates interior luminosity light curves using data from Radial Stellar Pulsations
# V1.1: Plots interior luminosity curves for different RSP convection sets on separate canvases, instead of 4 subplots on one canvas

library(av) # for animation
library(data.table) # for fread()

### CONTROLS ###

cwd <- getwd()
convection.sets <- c('A','B','C','D') # Also the loop vector for our plots
model_nums <- c('1')

RSP.nz <- 300 # number of zones

frame.rate <- 15 # frame rate of animation

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
  
  for (set in convection.sets) {
    plot.dir <- file.path("plots",paste0("Model_",model_num),"interior_light_curves","TEST",set)
    ifelse(!dir.exists(plot.dir), dir.create(plot.dir, recursive = TRUE), FALSE)
    
    label <- paste0(set, model_num)
    
    cat(paste("Creating interior LC plots for Model",label,"\n",sep=" "))
    
    log.directory.after.FA <- paste0("LOGS_",label,"a")
    history.path.after.FA <- file.path(cwd,log.directory.after.FA,"history.data")
    
    DF.after.FA <- fread(history.path.after.FA, header=TRUE, skip=5)
    prof.idx <- fread(file.path(log.directory.after.FA, 'profiles.index'), skip=1, 
                      col.names=c('mdl_num', 'priority', 'prof_num'))
    
    prof_num <- 1:length(prof.idx$prof_num)
    prof.path <- file.path(cwd,log.directory.after.FA, paste0('profile', prof_num, '.data'))
    DF.profile <- lapply(prof.path, function(x) fread(x, header=TRUE, skip = 5))
    
    phases <- DF.after.FA$rsp_phase
    
    o <- order(phases)
    
    for (zone_num in 1:RSP.nz){
      
      padded_zone_num <- formatC(zone_num, width = 4, flag = 0)
      
      ### LOOK INTO PROFILE FILES ###
      luminosities <- c()
      
      cat(paste("Reading profile files for zone:", zone_num, "\n"))
      for (prof_num in prof.idx$prof_num) {
        
        luminosity <- DF.profile[[prof_num]]$luminosity[zone_num]
        luminosities <- c(luminosities, luminosity)
        
      } # end of profile loop
      
      png(file.path(plot.dir,
                    paste("Light_Curve_","Zone_",padded_zone_num,".png",sep="")),
          width = 1200, height = 960)
      par(mar=c(5,4,4,2) + 2) # Control space between subplots
      
      plot.table <- data.frame(phases, luminosities)
      o <- order(phases)
      
      with(plot.table, plot(x=phases[o], y=luminosities[o],
                            main=paste("Zone",zone_num,"Light Curve",sep=" "),
                            type="l", pch=3, lwd = 6, col="purple", xlab=expression("Phase"),
                            ylab=expression("Luminosity " (L/L['\u0298'])), cex.main=1.60,
                            cex.lab=1.80, cex.axis=1.60))
      dev.off()
      
    } # end of zone loop

    setwd(plot.dir)
    animation.title <- paste0('CEP_MODEL_', model_num,'_SET_',set,'_INTERIOR_LUMINOSITY_LCs.mp4')
    av::av_encode_video(list.files(, '*.png'), framerate = frame.rate,
                        output = animation.title)
    # If this does not work, type in terminal ffmpeg -framerate 15 -i Light_Curve_Zone_%04d.png CEP_MODEL_1_INTERIOR_LCs.mp4
    cat(paste("Animation titled",animation.title,"saved in",plot.dir,"\n"))
    cat(paste("You may ignore warnings.\n"))
    setwd(cwd)
    
  } # end of convection set loop

} # end of model loop




