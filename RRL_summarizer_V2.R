# V2: Checking if RRL models are "good" based on mass, Teff 
# added boolean for sink(), as well as access to history and profile1

cwd = getwd()
#sets <- c('M079','M080','M081','M082','M0.822','M0.824','M0.826','M0.828','M083')
sets <- c('M0.821','M0.823','M0.825','M0.827','M0.829')
markers <- c('1','2','3','4','5')

total.models <- 0
used.files <- 0
unstable.counter <- 0

save.to.file <- TRUE

if (save.to.file) sink("LINA_summary.txt")
cat("This file summarizes RSP RR Lyrae models with unstable FU and FO growth rates, mass (0.55-0.65 Msun for good models), and effective temperature (6000-7250 K for good models). \n \n")
cat("See bottom of file for total number of \"good\" models. \n \n")
for (set in 1:length(sets)){
  for (marker in 1:length(markers)) {
    total.models <- total.models + 1
    cat(paste("**********************************************\nMass label:",sets[set],"\nMarker number:", markers[marker],"\n"))
    my.dir <- paste('LOGS/LOGS', sets[set], markers[marker], sep="_")
    total.dir <- file.path(cwd, my.dir)
    
    LINA.file <- file.path(total.dir, "LINA_period_growth.data")
    hstry.file <- file.path(total.dir, "history.data")
    prof1.file <- file.path(total.dir, "profile1.data")
    #cat(paste(LINA.file), "\n")
  
    if(!file.exists(LINA.file)){
      cat(paste("No linear calculations.","\n"))
      next
    } 
    used.files <- used.files + 1
    
    LINA.data <- read.table(LINA.file, header=TRUE)
    write.table(LINA.data)
    growth.rates <- LINA.data[, 2]
    FU.growth.rate <- growth.rates[1]
    FO.growth.rate <- growth.rates[2]
    cat(paste("Fundamental growth rate:",FU.growth.rate, "\n"))
    cat(paste("First overtone growth rate:",FO.growth.rate, "\n"))
    
    hstry.data <- read.table(hstry.file, skip=6)
    effective.T <- hstry.data[1,22]  # row first, column second
    cat(paste("Effective temperature:",effective.T, "\n"))
    
    prof1.data <- read.table(prof1.file, skip=2,  fill = TRUE)
    star.mass <- as.numeric(prof1.data[1,21])  # row first, column second
    cat(paste("Star mass:",star.mass, "\n"))
    
    if (FU.growth.rate > 0 | FO.growth.rate > 0) {
      cat(paste("UNSTABLE","\n"))
      if (effective.T>5800 & effective.T<7500) {
        if (star.mass>0.53 & star.mass<0.67){
          unstable.counter <- unstable.counter + 1
          cat("Good model. \n")
        }
      }
    } else {
      cat(paste("NOT UNSTABLE","\n"))
    }
    

  }
}


cat(paste("\nTotal RSP models:", total.models, "\n"))
cat(paste("Models with LINA calculations:", used.files, "\n"))
cat(paste("\"Good\" model files:", unstable.counter, "\n"))
if (save.to.file) sink()