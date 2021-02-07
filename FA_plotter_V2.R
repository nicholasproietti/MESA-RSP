# Code that creates full-amplitude plots and light curves of multiple bands using data from Radial Stellar Pulsations
# V2: Implemented function, data frame, and loop for colors. Corrected absolute magnitude for bug
# in MESA version 11701.

#cwd <- getwd()
cwd <- "/home/nick/mesa-r11701/star/test_suite/rsp_RR_Lyrae_grid"
sets <- c('RRL1','RRL2')

light_curve_plotter <- function(filepath, title, phases, values) {
  png(filepath, width = 750, height = 600)
  par(mar=c(5,5,4,1)+.1)
  plot(x=phases, y=values, main=title, pch=3, col="purple",
       xlab=expression("Phase"), ylab=expression("Absolute Magnitude (mags)") , cex.main=2.0, cex.lab=1.50, ceb.axis=1.80)
  dev.off()
}

color.strings <- c('B','V','I','J','H','K','R','U')

for (set in 1:length(sets)){
  cat(paste("Creating full-amplitude and color plots for model ",sets[set],"\n"))
  
  log.directory.before.FA <- paste("LOGS_",sets[set], sep = "")
  history.path.before.FA <- file.path(cwd,"light_curve_LOGS",log.directory.before.FA,"history.data")
  
  log.directory.after.FA <- paste("LOGS_",sets[set], "a", sep = "")
  history.path.after.FA <- file.path(cwd,"light_curve_LOGS",log.directory.after.FA,"history.data")
  
  DF.before.FA <- read.table(history.path.before.FA, header=1, skip=5)
  DF.after.FA <- read.table(history.path.after.FA, header=1, skip=5)
  
  ages <- DF.before.FA$star_age
  periods <- DF.before.FA$rsp_period_in_days
  growth.rates <- DF.before.FA$rsp_GREKM
  radii <- DF.before.FA$radius
  log_Ls <- DF.before.FA$log_L
  
  luminosities <- DF.after.FA$luminosity
  phases <- DF.after.FA$rsp_phase
  Bs <- DF.after.FA$abs_mag_B - 2.5*log10(3.846e33) # Lsun in ergs/s, to correct bug in MESA 11701
  Vs <- DF.after.FA$abs_mag_V - 2.5*log10(3.846e33)
  Is <- DF.after.FA$abs_mag_I - 2.5*log10(3.846e33)
  Js <- DF.after.FA$abs_mag_J - 2.5*log10(3.846e33)
  Hs <- DF.after.FA$abs_mag_H - 2.5*log10(3.846e33)
  Ks <- DF.after.FA$abs_mag_K - 2.5*log10(3.846e33)
  Rs <- DF.after.FA$abs_mag_R - 2.5*log10(3.846e33)
  Us <- DF.after.FA$abs_mag_U - 2.5*log10(3.846e33)
  
  colors = data.frame(B=Bs, V=Vs, I=Is, J=Js, H=Hs, K=Ks, R=Rs, U=Us)
  
  ### Look into profiles for acoustic depth information ###
  prof.idx <- read.table(file.path(cwd,"light_curve_LOGS",log.directory.after.FA, 'profiles.index'),
                         skip=1, col.names=c('mdl_num', 'priority', 'prof_num'))
  acoustic.depths <- c()
  
  cat(paste("Reading profile files... \n"))
  
  pb = txtProgressBar(min = 0, max = length(prof.idx$prof_num), style = 3)
  
  for (prof_num in prof.idx$prof_num) { 
    prof.path <- file.path(cwd,"light_curve_LOGS",log.directory.after.FA, paste0('profile', prof_num, '.data'))
    if (!file.exists(prof.path)) next 
    #print(prof.path)
    DF.profile <- read.table(prof.path, header=1, skip=5) 
    
    hif.idx <- min(which(DF.profile$neutral_fraction_H <= 0.5))
    
    acoustic.depth <- DF.profile$acoustic_depth[hif.idx]
    acoustic.depths <- c(acoustic.depths, acoustic.depth)
    setTxtProgressBar(pb, prof_num)
  }
  close(pb)
  ######
  cat(paste("Profile files completed!\n"))
  
  
  ### Create plot directories in case they do not already exist ###
  plot.dir <- file.path("plots")
  ifelse(!dir.exists(plot.dir), dir.create(plot.dir), FALSE)
  
  plot.set.dir <- file.path("plots",sets[set])
  ifelse(!dir.exists(plot.set.dir), dir.create(plot.set.dir), FALSE)
  ######
  
  
  png(file.path("plots",sets[set],paste(sets[set],"rspGREKM_vs_Star_Age.png",sep="_")), width = 750, height = 600)
  par(mar=c(5,5,4,1)+.1)
  plot(x=ages, y=growth.rates, main=paste(sets[set],"rspGREKM vs Star Age",sep=" "), pch=3, col="purple",
        xlab=expression("Star Age (days)"), ylab=expression("rspGREKM"), cex.main=2.0, cex.lab=1.50, ceb.axis=1.80)
  dev.off()

  png(file.path("plots",sets[set],paste(sets[set],"Period_vs_Star_Age.png",sep="_")), width = 750, height = 600)
  par(mar=c(5,5,4,1)+.1)
  plot(x=ages, y=periods, main=paste(sets[set],"Period vs Star Age",sep=" "), pch=3, col="purple",
       xlab=expression("Star Age (days)"), ylab=expression("Period (days)"), cex.main=2.0, cex.lab=1.50, ceb.axis=1.80)
  dev.off()

  png(file.path("plots",sets[set],paste(sets[set],"Radius_vs_Star_Age.png",sep="_")), width = 750, height = 600)
  par(mar=c(5,5,4,1)+.1)
  plot(x=ages, y=radii, main=paste(sets[set],"Radius vs Star Age",sep=" "), pch=3, col="purple",
       xlab=expression("Star Age (days)"), ylab=expression("Radius " (R/R['\u0298'])), cex.main=2.0, cex.lab=1.50, ceb.axis=1.80)
  dev.off()

  png(file.path("plots",sets[set],paste(sets[set],"Log_L_vs_Star_Age.png",sep="_")), width = 750, height = 600)
  par(mar=c(5,5,4,1)+.1)
  plot(x=ages, y=log_Ls, main=paste(sets[set],"Log(L) vs Star Age",sep=" "), pch=3, col="purple",
       xlab=expression("Star Age (days)"), ylab=expression("Log(L) " (L/L['\u0298'])), cex.main=2.0, cex.lab=1.50, ceb.axis=1.80)
  dev.off()

  png(file.path("plots",sets[set],paste(sets[set],"Light_Curve.png",sep="_")), width = 750, height = 600)
  par(mar=c(5,5,4,1)+.1)
  plot(x=phases, y=luminosities, main=paste(sets[set],"Light Curve",sep=" "), pch=3, col="purple",
       xlab=expression("Phase"), ylab=expression("Luminosity " (L/L['\u0298'])), cex.main=2.0, cex.lab=1.50, ceb.axis=1.80)
  dev.off()
  
  for (i in 1:ncol(colors)){
    light.curve.filepath <- file.path("plots",sets[set],paste(sets[set], color.strings[i],"Band_Light_Curve.png",sep="_"))
    light.curve.title <- paste(sets[set],color.strings[i],"Band Light Curve")
    light_curve_plotter(light.curve.filepath,light.curve.title,phases,colors[, i])
  }
  
  png(file.path("plots",sets[set],paste(sets[set],"HIF_Acoustic_Depth.png",sep="_")), width = 750, height = 600)
  par(mar=c(5,5,4,1)+.1)
  plot(x=phases, y=acoustic.depths, log="y", main=paste(sets[set],"HIF Acoustic Depth",sep=" "), pch=3, col="purple",
       xlab=expression("Phase"), ylab=expression("Acoustic depth \u03C4"), cex.main=2.0, cex.lab=1.50, ceb.axis=1.80)
  dev.off()
  
  cat(paste("Plots for model",sets[set],"completed! \n"))
}
