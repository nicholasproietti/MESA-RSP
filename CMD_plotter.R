# # Code that creates color-magnitude diagrams using data from Modules for Experiments in Stellar Astrophysics

cwd = getwd()
#sets <- c('A','B','C','D')
sets <- c('A')

for (set in 1:length(sets)){
  cat(paste("Creating Color-magnitude diagram for Set",sets[set],"\n"))
  CMD.name <- paste('CMD_Set',sets[set],'.png',sep="")
  set.dir <- paste('LOGS_', sets[set], sep="")
  parent.dir <- file.path(cwd, set.dir)
  dirs <- list.files(parent.dir)
  
  abs_V <- c()
  abs_I <- c()
  
  abs_V_minus_I <- c()
  
  used.files <- 0
  
  pb = txtProgressBar(min = 0, max = length(dirs), style = 3)
  
  for (i in 1:length(dirs)) {
    LINA.file <- file.path(parent.dir, dirs[i], "LINA_period_growth.data")
    current.file <- file.path(parent.dir, dirs[i], "history.data")
    
    if(!file.exists(LINA.file)) next
    LINA.data <- read.table(LINA.file, header=TRUE)
    growth.rates <- LINA.data[, 2]
    if (growth.rates[1] <= 0) next
    
    if (!file.exists(current.file)) next
      
    used.files <- used.files + 1
    
    data <- read.table(current.file, skip=6)
    #cat(paste("Loaded ", current.file, "\n"))
    
    # NOTE: indexing for columns starts at 1, NOT 0
    temp_V <- log10(data[, 42])
    temp_I <- log10(data[, 43])
    
    temp_V_minus_I = temp_V - temp_I
    
    abs_V <- append(abs_V, temp_V)
    abs_I <- append(abs_I, temp_I)
    abs_V_minus_I <- append(abs_V_minus_I,temp_V_minus_I)
    #cat(paste("Completed", current.file, "\n", "\n"))
    
    setTxtProgressBar(pb, i)
    #print(length(abs_V_minus_I))
  }
   
  x_min <- min(abs_V_minus_I, na.rm=TRUE)
  x_max <- max(abs_V_minus_I, na.rm=TRUE)
  
  xlim <- c(x_min, x_max)
  ylim <- rev(range(abs_V))
  
  current.file <- file.path(parent.dir, "LOGS_A9999", "history.data")
  data <- read.table(current.file, skip=6)
  temp_V <- log10(data[, 42])
  temp_I <- log10(data[, 43])
  print(current.file)
  print(temp_V)
  print(temp_I)
  #png(CMD.name)
  #plot(x=abs_V_minus_I, y=abs_V, main=paste("Set",sets[set],"Color-magnitude diagram",sep=" "), pch=3, col="purple",
  #      xlab=expression("M"[V]*"-M"[I]*" (mag)"), ylab=expression("M"[V]*" (mag)"), xlim=xlim, ylim=ylim)
  #dev.off()
  
  cat(paste("\nColor magnitude diagram saved as", CMD.name, "in", cwd, "\n"))
  cat(paste("Total model directories:", length(dirs), "\n"))
  cat(paste("Total positive FU growth rate model files:", used.files, "\n"))
}
