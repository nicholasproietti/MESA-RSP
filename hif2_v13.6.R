#### Hydrogen ionization front calculations 
#### Author: Earl Bellinger ( https://earlbellinger.com ) 
#### Current Editor: Nicholas Proietti ( I do not have a website. )

### V2: Correlate the points between hif and hr with a marker
### V3: Automate positions of markers in plots, automate labels
### V13: improved variable names, implemented adaptable HR diagram scales, added blueward extension (max HB Teff)
### calculation, created taus.dat files for different models, changed lines() to no longer interpolate,
### added min HIF-photosphere distance calculation, removed envelope structure diagram, added zoom plot functionality (code
### outputs plots for the markers normally and zoomed in
### NOTE: V13 was modified directly from V3 
### V13.2: Zoom part of HR diagram now zooms into HB to see loops
### NOTE: V13.2 was modified directly from V13, not V13.1
### 13.3: Added the ability to create Log_LHe vs model number plots, or central He vs model number
### added output for effective temperature and mass of marker
### 13.5: 1/24/21 acoustic depth plots added
### 13.6: 1/29/21 Implemented RSP output

########################################################################################
############################## HIF Plot ################################################
########################################################################################

### LIBRARIES 
source('utils.R') 
options(scipen=10000)
library(magicaxis)
library(plotrix)
library(pracma) # for ceil()
library(sjmisc) # for str.contains()
library(spatstat.utils) # for inside.range()

model.label <- 'M081_acoustic'

logs_dir <- paste('LOGS_', model.label, sep="")
prof.idx <- read.table(file.path(logs_dir, 'profiles.index'), skip=1, 
                       col.names=c('mdl_num', 'priority', 'prof_num'))
hstry <- read.table(file.path(logs_dir, 'history.data'), header=1, skip=5)

TAMS <- with(hstry, star_age[min(which(center_h1 <= 1e-5))]/10**9)

rgb <- with(hstry, center_h1 <= 1e-5 & center_he4 > 1e-5)
RGB.TIP <- with(hstry[rgb,], star_age[which.max(log_L)]/10**9)

AGB <- with(hstry, 
            star_age[min(which(center_h1 <= 1e-5 & center_he4 <= 1e-5))]/10**9)

min_HB <- RGB.TIP
max_HB <- AGB

taus <- c()
ages <- c()
logqs <- c()
prof_nums <- c()
logLHes <- c()
model_nums <- c()
central_Hes <- c()
acoustic_depths <- c()

HB_ages <- c()
HB_taus <- c()
HB_prof_nums <- c()
HB_Teffs <- c()
HB_Ls <- c()
HB_model_nums <- c()

lower.track.prof.nums <- c()
lower.track.Teffs <- c()
lower.track.Ls <- c()
lower.track.model.nums <- c()

higher.track.prof.nums <- c()
higher.track.Teffs <- c()
higher.track.Ls <- c()
higher.track.model.nums <- c()

taus.2 <- c() # to avoid getting taus far into AGB

cat(paste("Loading profiles for", model.label))
pb = txtProgressBar(min = 0, max = length(prof.idx$prof_num), style = 3)

for (prof_num in prof.idx$prof_num) { 
    
    prof.path <- file.path(logs_dir, paste0('profile', prof_num, '.data'))
    if (!file.exists(prof.path)) next 
    #print(prof.path)
    DF <- read.table(prof.path, header=1, skip=5)
    hstry. <- hstry[hstry$model_number == 
                        prof.idx[prof.idx$prof_num == prof_num,]$mdl_num,]
    
    hif.idx <- min(which(DF$neutral_fraction_H <= 0.5))
    
    logq <- log10(1-DF$mass[hif.idx] / hstry.$star_mass)
    tau <- DF$opacity[hif.idx]
    age <- hstry.$star_age / 10**9
    
    # 13.3: To create log_He vs model number plots
    logLHe <- hstry.$log_LHe
    model_num <- hstry.$model_number
    central_He <- hstry.$center_he4
    
    #ages <- c(ages, hstry.$star_age / 10**9)
    ages <- c(ages, age)
    logqs <- c(logqs, logq)
    taus <- c(taus, tau)
    prof_nums <- c(prof_nums, prof_num)
    
    # 13.3: To create log_He vs model number plots
    logLHes <- c(logLHes, logLHe)
    model_nums <- c(model_nums, model_num)
    central_Hes <- c(central_Hes, central_He)
    
    # 13.5: acoustic depth plots
    temp_acoustic_depth <- DF$acoustic_depth[hif.idx]
    acoustic_depths <- c(acoustic_depths, temp_acoustic_depth)
    
    if (length(age) == 0) next

    if (age < (max_HB + 0.5) & age != 0){
        # for determining the minimum HIF-photosphere distance
        taus.2 <- c(taus.2, tau)
        
        # for HB calculations
        if (age > min_HB & age < max_HB & age != 0){
            HB_ages <- c(HB_ages, age)
            HB_taus <- c(HB_taus, tau) 
            HB_prof_nums <- c(HB_prof_nums, prof_num)
            HB_model_nums <- c(HB_model_nums, model_num)
            
            current_Teff <- 10**(hstry.$log_Teff)
            HB_Teffs <- c(HB_Teffs, current_Teff)
            
            current_L <- hstry.$log_L
            HB_Ls <- c(HB_Ls, current_L)
            
            #13.5 lower track
            if (age < 12.8602339994746){
                lower.track.prof.nums <- c(lower.track.prof.nums, prof_num)
                lower.track.model.nums <- c(lower.track.model.nums, model_num)
                lower.track.Ls <- c(lower.track.Ls, current_L)
                lower.track.Teffs <- c(lower.track.Teffs, current_Teff)
                
            }
            
            #13.5 higher track
            if (age >= 12.8602339994746){
                higher.track.prof.nums <- c(higher.track.prof.nums, prof_num)
                higher.track.model.nums <- c(higher.track.model.nums, model_num)
                higher.track.Ls <- c(higher.track.Ls, current_L)
                higher.track.Teffs <- c(higher.track.Teffs, current_Teff)
                
            }
            
        }
    }
    
    setTxtProgressBar(pb, prof_num)
}
close(pb)

LHes <- 10**logLHes

# even out vectors for data frame
if (length(ages) != length(taus)){
    ages <- ages[1:length(taus)]
    logqs <- logqs[1:length(taus)]
}
# even out vectors for data frame
if (length(model_nums) != length(ages)){
    logLHes <- logLHes[1:length(ages)]
    LHes <- LHes[1:length(ages)]
    model_nums <- model_nums[1:length(ages)]
}


# Find the range of Teff and L in HB
max.Teff <- max(HB_Teffs)
min.Teff <- min(HB_Teffs)
max.L <- max(HB_Ls)
min.L <- min(HB_Ls)
min.HBMN <- min(HB_model_nums)
max.HBMN <- max(HB_model_nums)
mid.range.HB <- min(HB_ages) + (max(HB_ages)-min(HB_ages))/2

# Print some information about the model
cat(paste("\nMaximum HB Teff (blueward extension):",max.Teff,"\n"))
cat(paste("Minimum HB Teff:",min.Teff,"\n"))
cat(paste("Maximum HB L:",max.L,"\n"))
cat(paste("Minimum HB L:",min.L,"\n"))
cat(paste("Maximum HB age:",max(HB_ages),"\n"))
cat(paste("Minimum HB age:",min(HB_ages),"\n"))
cat(paste("Mid-range HB age:",mid.range.HB,"\n"))
cat(paste("Maximum HB model number:",max.HBMN,"\n"))
cat(paste("Minimum HB model number:",min.HBMN,"\n"))
cat(paste("Maximum LHe:",max(LHes),"\n"))
cat(paste("Minimum LHe:",min(LHes),"\n"))
min.HIF.p.index <- which.min(abs(taus.2-(2/3)))
cat(paste("HIF-photosphere distance:",taus.2[min.HIF.p.index],"\n"))

# Zoom settings
#zoom.set <- c(FALSE, TRUE)
zoom.set <- c(FALSE)
plot.set <- c("tau", "LHe","acoustic")
marker.set <- c("1","2","3","4","5")

# We decide which point to mark here
dst <- sqrt((log10(HB_Teffs) - 6400)^2 + (HB_Ls - 60)^2)
my.index <- which.min(dst)#which.min(abs(higher.track.Teffs-6800))
special.prof.num <- 1399#HB_prof_nums[my.index]

# Models to run for RSP
#limit <- 30
#lower.prof.limit <- special.prof.num - limit
#higher.prof.limit <- special.prof.num + limit
RSP.labels <- c()
RSP.masses <- c()
RSP.luminosities <- c()
RSP.Teffs <- c()
RSP.index <- 1
for (temp.prof.num in 1399:1448){
    cat(paste("Loop index: ",RSP.index))
    hif_prof_num <- temp.prof.num
    prof.path <- file.path(logs_dir, paste0('profile', hif_prof_num, '.data'))
    if (!file.exists(prof.path)) next 
    print(prof.path)
    DF <- read.table(prof.path, header=1, skip=5)
    hstry. <- hstry[hstry$model_number == 
                        prof.idx[prof.idx$prof_num == hif_prof_num,]$mdl_num,]
    hif.idx <- min(which(DF$neutral_fraction_H <= 0.5))
    
    # find values for markers
    marker.age <- hstry.$star_age / 10**9
    marker.tau <- DF$opacity[hif.idx]
    marker.model_num <- hstry.$model_number
    marker.LHe <- 10**hstry.$log_LHe
    marker.temp <- 10**hstry.$log_Teff
    marker.mass <- hstry.$star_mass
    marker.lum <- hstry.$luminosity
    marker.acoustic_depth <- DF$acoustic_depth[hif.idx]
    
    cat(paste("marker age:", marker.age,"\n"))
    cat(paste("marker tau:", marker.tau,"\n"))
    cat(paste("marker model number:", marker.model_num,"\n"))
    cat(paste("marker LHe:", marker.LHe,"\n"))
    cat(paste("marker effective temperature:", marker.temp,"\n"))
    cat(paste("marker mass:", marker.mass,"\n"))
    cat(paste("marker luminosity:", marker.lum,"\n"))
    cat(paste("marker acoustic depth:", marker.acoustic_depth,"\n"))
    
    RSP.labels <- c(RSP.labels,temp.prof.num)
    RSP.masses <- c(RSP.masses,marker.mass)
    RSP.luminosities <- c(RSP.luminosities,marker.lum)
    RSP.Teffs <- c(RSP.Teffs,marker.temp)
    RSP.index <- RSP.index + 1
    cat("\n")
}
RSP.Zs <- rep(c(0.001),times=(length(RSP.luminosities)))
RSP.Xs <- rep(c(0.749),times=(length(RSP.luminosities)))
RSP.DF <- data.frame(RSP.Zs, RSP.Xs, RSP.masses, RSP.luminosities, RSP.Teffs)
row.names(RSP.DF) <- NULL
colnames(RSP.DF) <- NULL
write.table(RSP.DF, file="input.dat")



# Now find features of our special profile number
#last_prof_num <- prof_nums[length(prof_nums)]
hif_prof_num <- special.prof.num
prof.path <- file.path(logs_dir, paste0('profile', hif_prof_num, '.data'))
if (!file.exists(prof.path)) next 
print(prof.path)
DF <- read.table(prof.path, header=1, skip=5)
hstry. <- hstry[hstry$model_number == 
                    prof.idx[prof.idx$prof_num == hif_prof_num,]$mdl_num,]
hif.idx <- min(which(DF$neutral_fraction_H <= 0.5))

# find values for markers
marker.age <- hstry.$star_age / 10**9
marker.tau <- DF$opacity[hif.idx]
marker.model_num <- hstry.$model_number
marker.LHe <- 10**hstry.$log_LHe
marker.temp <- 10**hstry.$log_Teff
marker.mass <- hstry.$star_mass
marker.lum <- hstry.$luminosity
marker.acoustic_depth <- DF$acoustic_depth[hif.idx]

cat(paste("marker age:", marker.age,"\n"))
cat(paste("marker tau:", marker.tau,"\n"))
cat(paste("marker model number:", marker.model_num,"\n"))
cat(paste("marker LHe:", marker.LHe,"\n"))
cat(paste("marker effective temperature:", marker.temp,"\n"))
cat(paste("marker mass:", marker.mass,"\n"))
cat(paste("marker luminosity:", marker.lum,"\n"))
cat(paste("marker acoustic depth:", marker.acoustic_depth,"\n"))

tau.file.name <- paste("taus/taus_",model.label,".dat",sep='')
LHe.file.name <- paste("LHes/LHes_",model.label,".dat",sep='')
acous.file.name <- paste("acoustic_depths/acoustic_depths_",model.label,".dat",sep='')

tau.DF <- data.frame(ages, taus, logqs)
write.table(tau.DF, file=tau.file.name)
tau.DF <- read.table(tau.file.name, header=1)

He.DF <- data.frame(model_nums, logLHes, LHes)
write.table(He.DF, file=LHe.file.name)
He.DF <- read.table(LHe.file.name, header=1)

acous.DF <- data.frame(ages, acoustic_depths, taus)
write.table(acous.DF, file=acous.file.name)
acous.DF <- read.table(acous.file.name, header=1)

for (p.i in 1:length(plot.set)) {

    options(scipen=1)
    plot_hif_tau <- function(..., make.x=T, make.y=T, 
                             text.cex=1, mgp=utils.mgp, font=utils.font, mar=utils.mar, short=F) {
        cat(plot.set[p.i])
        
        make.x=T
        make.y=T
        text.cex=1
        mgp=utils.mgp
        font=utils.font
        mar=utils.mar
        short=F
        
        par(mar=mar+c(2, 1, 2, 0.5), lwd=3.66, las=1, cex.axis=text.cex)
        
        tau.max <- round(max(taus))
        ymax <- 10**(ceil(log10(tau.max)))
        tau.min <- round(min(taus)/100)
        ymin <- 10**(ceil(log10(tau.min)))
        
        acous.max <- round(max(acoustic_depths))
        y.acous.max <- 10**(ceil(log10(acous.max)))
        
        if (plot.set[p.i] == "tau"){
            xlim <- c(TAMS-0.7,AGB+1.0)
            #xlim <- c(0, max(hstry$star_age)/10**9)
            ylim <- c(1e-1, ymax) #log10(c(0.1, 1e7))
        } else if (plot.set[p.i] == "LHe"){
            #xlim <- c(TAMS,AGB)
            xlim <- c(min(model_nums), max(model_nums)) #11350 for start of HB, 15030 for start of AGB
            #xlim <- c(0, max(hstry$star_age)/10**9)
            ylim <- c(min(LHes), max(LHes)) #log10(c(0.1, 1e7))
            
            # zoom settings
            # xlim <- c((hstry.$star_age / 10**9)-0.001,(hstry.$star_age / 10**9)+0.001)
            # ylim <- c(min(taus), ymax)
        } else if (plot.set[p.i] == "acoustic"){
            xlim <- c(TAMS-0.7,AGB+1.0)
            #xlim <- c(0, max(hstry$star_age)/10**9)
            ylim <- c(1e-2, y.acous.max)
            
            # zoom settings
            # xlim <- c((hstry.$star_age / 10**9)-0.001,(hstry.$star_age / 10**9)+0.001)
            # ylim <- c(min(taus), ymax)
        }
        
        
        plot(NA, 
             axes=F, xaxs='i', yaxs='i', log='y', 
             type='l', lwd=3, col=1, 
             xlim=xlim, 
             ylim=ylim, 
             xlab="", 
             ylab="")
        
        ## EVOLUTIONARY STAGE 
        TAMS <- with(hstry, star_age[min(which(center_h1 <= 1e-5))]/10**9)
        #abline(v=TAMS, lty=2, lwd=par()$lwd)
        
        rgb <- with(hstry, center_h1 <= 1e-5 & center_he4 > 1e-5)
        RGB.TIP <- with(hstry[rgb,], star_age[which.max(log_L)]/10**9)
        #abline(v=RGB.TIP, lty=2, lwd=par()$lwd)
        
        AGB <- with(hstry, 
                    star_age[min(which(center_h1 <= 1e-5 & center_he4 <= 1e-5))]/10**9)
        #abline(v=AGB, lty=2, lwd=par()$lwd)
        
        
        spectral.divs <- c(TAMS, RGB.TIP, AGB, xlim[2])
        #rs <- c(175/255, 199/255, 1, 1, 1, 1, 1, 1)
        #gs <- c(201, 216, 244, 229, 217, 199, 166)/255
        #bs <- c(1, 1, 243/255, 207/255, 178/255, 142/255, 81/255)
        rs <- c(175/255, 199/255, 1, 175/255, 1, 1, 1, 1)
        gs <- c(201, 216, 244, 201, 217, 199, 166)/255
        bs <- c(1, 1, 243/255, 1, 178/255, 142/255, 81/255)
        cols <- c(
            rgb(175/255, 201/255, 1),       # O 1
            rgb(199/255, 216/255, 1),       # B 2
            rgb(1,       244/255, 243/255), # A 3
            rgb(175/255, 201/255, 1),       # O 1#rgb(1,       229/255, 207/255), # F 4
            rgb(1,       217/255, 178/255), # G 5
            rgb(1,       199/255, 142/255), # K 6
            rgb(1,       166/255, 81/255))  # M 7
        #idxs <- c(4, 6, 2, 5)
        #rs <- rs[idxs]
        #gs <- gs[idxs]
        #bs <- bs[idxs]
        #cols <- cols[idxs]
        #cols <- c("#E8E9EB", "#E4B363", "#EF6461", "#E0DFD5")
        for (ii in 1:length(spectral.divs)) {
            div <- spectral.divs[ii]
            if (div < xlim[1]) next 
            if (div > xlim[2]) div <- xlim[2]
            if (ii == 1) {
                rect(xlim[1], ylim[1], div, ylim[2], col=cols[ii], border=NA)
            } else {
                prev <- spectral.divs[ii-1]
                if (prev < xlim[1]) prev <- xlim[1]
                rect(prev, ylim[1], div, ylim[2], col=cols[ii], border=NA)
            }
        }
        for (ii in 1:(length(spectral.divs)-1)) {
            div <- spectral.divs[ii]
            gradient.rect(div-0.02, ylim[1], div+0.02, ylim[2],
                          nslices=10, border=NA, 
                          reds=c(rs[ii], rs[ii+1]), 
                          greens=c(gs[ii], gs[ii+1]),
                          blues=c(bs[ii], bs[ii+1]))
        }
        
        
        #rect(xlim[1], ylim[1], TAMS, ylim[2], col='darkgray', border=NA)
        
        if (plot.set[p.i] == "tau"){
            ## PHOTOSPHERE 
            phot.col <- "#313638" #blue #"#d62828"
            abline(h=2/3, lwd=5, lty=2, col=phot.col)
            text(TAMS-0.5, 2/3+7,#1.25, 
                 expression(Photosphere),
                 cex=.97*text.cex, pos=4, 
                 family='Helvetica LT Std Bold', col=phot.col)
            text(#6.85, 2/3+1.25, 
                TAMS-0.5, 2/3+1.1, 
                expression((tau==2/3)),
                cex=.8*text.cex, pos=4, 
                family='Helvetica LT Std Light', col=phot.col)
        } 
        
        ## HIF 
        hif.col <- 1 #"#005e91"
        
        if (plot.set[p.i] == "tau"){
            text(TAMS-0.5, tau.DF[min(which(tau.DF$ages>=xlim[1])),]$taus * 2.45, 
                 expression(HIF), col=hif.col, 
                 cex=.97*text.cex, pos=4, 
                 family='Helvetica LT Std Bold')
        } else if (plot.set[p.i] == "acoustic"){
            text(TAMS-0.5, acous.DF[min(which(acous.DF$ages>=xlim[1])),]$acoustic_depths * 2.45, 
                 expression(HIF), col=hif.col, 
                 cex=.97*text.cex, pos=4, 
                 family='Helvetica LT Std Bold')            
        }
        tau.spl <- with(tau.DF, splinefun(ages, taus))
        age.seq <- seq(xlim[1], xlim[2], length.out=100000)#00)
        
        if (plot.set[p.i] == "tau"){
        lines(ages, taus, lwd=5, col=hif.col)
        points(marker.age, marker.tau, pch=21, cex=1.7, 
               col="#FFFFFF", lwd=par()$lwd,
               bg=with(DF[nrow(DF),], 
                       rgb(red=y*.8, green=x*.8, blue=z*.8)))
        } else if (plot.set[p.i] == "LHe"){
            lines(model_nums, LHes, lwd=5, col=hif.col)
            points(marker.model_num, marker.LHe, pch=21, cex=1.7, 
                   col="#FFFFFF", lwd=par()$lwd,
                   bg=with(DF[nrow(DF),], 
                           rgb(red=y*.8, green=x*.8, blue=z*.8)))
        } else if (plot.set[p.i] == "acoustic"){
            lines(ages, acoustic_depths, lwd=5, col=hif.col)
            points(marker.age, marker.acoustic_depth, pch=21, cex=1.7, 
                   col="#FFFFFF", lwd=par()$lwd,
                   bg=with(DF[nrow(DF),], 
                           rgb(red=y*.8, green=x*.8, blue=z*.8)))
        }
        
        ## AXES 
        tcl <- -0.5
        par(mgp=mgp+c(0, 1.11, 0))
        
        magaxis(1, tcl=tcl, lwd=0, lwd.ticks=par()$lwd, ticks=T, labels=make.x,
                cex.axis=text.cex, family=font, mgp=mgp+c(0, 1.11, 0))
        magaxis(2, tcl=tcl, lwd=0, lwd.ticks=par()$lwd, ticks=T, labels=make.y,
                cex.axis=text.cex, family=font, mgp=mgp+c(0, 0.76, 0))
        
        axis(3, tcl=tcl, lwd=0, lwd.ticks=par()$lwd, tick=T, 
             at=c(TAMS, RGB.TIP, AGB), labels=F)
        
        box(lwd=par()$lwd)
        
        par(mgp=mgp+c(1.8, 0, 0))
        if (plot.set[p.i] == "tau"){
            title(ylab=expression(Optical~depth~tau))
            par(mgp=mgp+c(2.4, 0, 0))
            title(xlab=expression(Star~age~t/Gyr))
            par(mgp=mgp+c(0, 0.72, 0))
        } else if (plot.set[p.i] == "LHe"){
            title(ylab=expression(log~LHe))
            par(mgp=mgp+c(2.4, 0, 0))
            title(xlab=expression(Model~number))
            par(mgp=mgp+c(0, 0.72, 0))
        } else if (plot.set[p.i] == "acoustic"){
            title(ylab=expression(Acoustic~depth~tau))
            par(mgp=mgp+c(2.4, 0, 0))
            title(xlab=expression(Star~age~t/Gyr))
            par(mgp=mgp+c(0, 0.72, 0))
        }
        if (plot.set[p.i] == "tau" | plot.set[p.i] == "acoustic"){
        axis(3, tick=F, labels=c('MS', 'RGB', 'HB', 'AGB'),
             at=c((xlim[1]+TAMS)/2,
                  (TAMS+RGB.TIP)/2,
                  (RGB.TIP+AGB)/2,
                  (AGB+xlim[2])/2),
        )
        }
    }
    
    #plot_hif_tau() 
    
    #zoom.label <- ''
   # if (zoom.set[z.i]){
   #     zoom.label <- 'HB_zoom'
   # } else {
    #    zoom.label <- ''
    #}
    
    plot.label <- ''
    if (plot.set[p.i] == "tau"){
        plot.label <- 'hif-tau-RRL'
    } else if (plot.set[p.i] == "LHe"){
        plot.label <- 'log_LHe-RRL'
    } else if (plot.set[p.i] == "acoustic"){
        plot.label <- 'hif-acoustic-RRL'
    }
    
    make_plots(plot_hif_tau, paste(plot.label,model.label,sep="_"),
               filepath=file.path(paste('plots_',model.label,sep=""), 'tau'),
               cex.paper=1.7, short=F, make_png=T, make_pdf=F,
               wide=T, slides=F, #make_png=T, make_pdf=F,
               thin=F, 
               make.x=T,
               make.y=T,
               paper_pdf_height=5.25,
               font="Palatino Linotype",#"Helvetica",# 
               use.cairo=T)
} #end of zoom loop
    ########################################################################################
############################## Hertzsprung-Russell Diagram #############################
########################################################################################

### LIBRARIES 
source('utils.R') 
options(scipen=10000)
library(plotrix)

logs_dir <- paste('LOGS_', model.label, sep="")

prof.idx <- read.table(file.path(logs_dir, 'profiles.index'), skip=1, 
    col.names=c('mdl_num', 'priority', 'prof_num'))

hstry <- read.table(file.path(logs_dir, 'history.data'), header=1, skip=5)

decreasing_L <- with(hstry, 
    which((diff(log_L) < 0)
        & center_h1[-1] > 0.65))
if (any(decreasing_L)) {
    goes_back_up <- diff(decreasing_L) > 1
    pms <- max(decreasing_L)
    hstry <- hstry[-1:-pms,]
}
decreasing_Teff <- with(hstry, 
    which((diff(log_Teff) < 0)
        & center_h1[-1] > 0.65))
if (any(decreasing_Teff)) {
    goes_back_up <- diff(decreasing_Teff) > 1
    pms <- max(decreasing_Teff)
    hstry <- hstry[-1:-pms,]
}

prof_num <- special.prof.num

DF <- read.table(file.path(logs_dir, paste0('profile', prof_num, '.data')), 
    header=1, skip=5)
hstry. <- hstry[hstry$model_number == 
    prof.idx[prof.idx$prof_num == prof_num,]$mdl_num,]

zoom.set <- c(FALSE, TRUE)
for (z.i in 1:length(zoom.set)) {

plot_hif <- function(..., make.x=T, make.y=T,
        text.cex=1, mgp=utils.mgp, font=utils.font, mar=utils.mar, short=F) {
    
    par(mar=mar+c(0.3, -0.5, 2, -0.1), lwd=1.66, las=1, cex.axis=text.cex,
        mfrow=c(1,2))
    
    ### HRD 
    
    if (!zoom.set[z.i]){
        #default xlim
        #xlim <- log10(c(round(max.Teff,-3)+1000, 3300))
        
        #to view marker in planetary nebula stage
        xlim <- log10(c(round(10**hstry.$log_Teff,-2)+5000, 3300))
        ylim <- c(-0.5, 4.6)
    } else {
        # test
        xlim <- log10(c(round(max.Teff,-2)+500, round(min.Teff,-2)+500))
        ylim <- c(min.L-0.1, max.L-1.1)
        
        # to zoom around HB
        #xlim <- log10(c(10100, 5100))
        #ylim <- c(1.4, 2.0)
        
        # to zoom around marker
        #xlim <- log10(c(round(10**hstry.$log_Teff,-2)+1000, round(10**hstry.$log_Teff,-2)-1000))
        #ylim <- c(round(hstry.$log_L,1)-0.5, round(hstry.$log_L,1)+0.5)
    }
    
    plot(NA, 
        axes=F, xaxs='i', yaxs='i', 
        type='l', lwd=3, col=1, 
        xlim=xlim, 
        ylim=ylim, 
        xlab="", 
        ylab="")
    
    spectral.divs <- log10(c(30000, 10000, 7500, 6000, 5200, 3700, 2400))
    rs <- c(175/255, 199/255, 1, 1, 1, 1, 1, 1)
    gs <- c(201, 216, 244, 229, 217, 199, 166)/255
    bs <- c(1, 1, 243/255, 207/255, 178/255, 142/255, 81/255)
    cols <- c(
        rgb(175/255, 201/255, 1),       # O
        rgb(199/255, 216/255, 1),       # B
        rgb(1,       244/255, 243/255), # A 
        rgb(1,       229/255, 207/255), # F 
        rgb(1,       217/255, 178/255), # G 
        rgb(1,       199/255, 142/255), # K 
        rgb(1,       166/255, 81/255))  # M
    for (ii in 1:length(spectral.divs)) {
        div <- spectral.divs[ii]
        if (div > xlim[1]) next  
        if (div < xlim[2]) div <- xlim[2]
        if (ii == 1) {
            rect(xlim[1], ylim[1], div, ylim[2], col=cols[ii], border=NA)
        } else {
            prev <- spectral.divs[ii-1]
            if (prev > xlim[1]) prev <- xlim[1]
            rect(prev, ylim[1], div, ylim[2], col=cols[ii], border=NA)
        }
    }
    for (ii in 2:(length(spectral.divs)-1)) {
        div <- spectral.divs[ii]
        if (div > xlim[1]) next  
        if (div < xlim[2]) div <- xlim[2]
        gradient.rect(div+0.0025, ylim[1], div-0.0025, ylim[2],
            nslices=10, border=NA, 
            reds=c(rs[ii], rs[ii+1]), 
            greens=c(gs[ii], gs[ii+1]),
            blues=c(bs[ii], bs[ii+1]))
    }
    
    lines(hstry$log_Teff, hstry$log_L, type='l', lwd=2, col=1)
    points(hstry$log_Teff[1], hstry$log_L[1], pch=20, cex=0.5, col=1)
    
    mdl.col <- cols[which.min(hstry.$log_Teff < spectral.divs)]
    mdl.cex <- 1+hstry.$log_R
    points(hstry.$log_Teff, hstry.$log_L, pch=21, cex=mdl.cex, 
        col="#FFFFFF", lwd=par()$lwd,
        bg=with(DF[nrow(DF),], 
            rgb(red=y*.8, green=x*.8, blue=z*.8)))
    
    age <- round(signif(hstry.$star_age/10**9, 4), 5)
    mass <- signif(hstry.$star_mass, 3)
    par(family="Helvetica LT Std Light")
    legend('topleft', 
        cex=0.8*text.cex,
        bty='n', inset=c(-0.05, 0),
        legend=c(as.expression(bquote(tau/Gyr==.(age))),
                 as.expression(bquote(M/M["sun"]==.(mass)))))
    par(family=font)
    
    nxticks <- 3
    nyticks <- 4
    nxminor <- 5
    nyminor <- 4
    xticks <- pretty(xlim, n=nxticks)
    yticks <- pretty(ylim, n=nyticks)
    xticks.minor <- pretty(xlim, n=nxticks*nxminor)
    yticks.minor <- pretty(ylim, n=nyticks*nyminor)
    xticks.minor <- xticks.minor[!xticks.minor %in% xticks]
    yticks.minor <- yticks.minor[!yticks.minor %in% yticks]
    par(mgp=mgp+c(0, 0.3, 0))
    xpos <- seq(10**xlim[2], 10**xlim[1], 2000)
    xpos2 <- seq(10**xlim[2], 10**xlim[1], 400)
    axis(side=1, tcl=-0.346/2, at=log10(xpos2), 
        labels=F, lwd.ticks=par()$lwd)
    axis(side=1, tcl=-0.346, at=log10(xpos), 
        labels=xpos, cex.axis=text.cex, lwd.ticks=par()$lwd)
    par(mgp=mgp+c(0, 0.43, 0))
    
    magaxis(2, tcl=-0.346, lwd=0, lwd.ticks=par()$lwd, unlog=T, family=font,
        mgp=mgp+c(0, 0.43, 0))
    #axis(2, tcl=-0.346, lwd=0, lwd.ticks=par()$lwd, tick=T, at=yticks,
    #    labels=as.logical(make.y))
    #axis(2, tcl=-0.346/2, lwd=0, lwd.ticks=par()$lwd, tick=T, 
    #    at=yticks.minor, labels=F)
    
    box(lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    if (make.x) title(xlab=expression(T["eff"]/K))
    par(mgp=mgp+c(0.55, 0, 0))
    if (make.y) title(ylab=expression(L/L["sun"]))
    if (make.x) mtext(expression("Spectral Type"), 
        side=3, line=1.5, cex=text.cex)
    
    spectral.labs <- c("O", "B", "A", "F", "G", "K", "M")
    selector <- 1:(length(spectral.divs))
    spectral.Teffs <- sapply(selector, 
        function(ii) {
            div <- spectral.divs[ii]
            if (div > xlim[1]) return(Inf)
            if (div < xlim[2]) div <- xlim[2]
            if (ii == 1) return((xlim[1]+div)/2)
            prev <- spectral.divs[ii-1]
            if (prev > xlim[1]) prev <- xlim[1]
            (div+prev)/2
        })
    axis(3, at=spectral.divs, tcl=-0.346, labels=F, cex.axis=text.cex,
        lwd.ticks=par()$lwd)
    par(mgp=mgp+c(0, -0.1, 0))
    axis(3, at=spectral.Teffs, labels=spectral.labs[selector], 
        cex.axis=text.cex, tcl=0)

}

#plot_hif()

zoom.label <- ''
if (zoom.set[z.i]){
    zoom.label <- 'HB_zoom'
} else {
    zoom.label <- ''
}

make_plots(plot_hif, paste('hr', model.label, zoom.label, sep = "_"),
        filepath=file.path(paste('plots_',model.label,sep="")),
        cex.paper=0.93, 
        wide=T, thin=F, short=F, slides=F, make_png=T, #make_pdf=F,
        make.x=T,
        make.y=T,
        font="Palatino Linotype", 
        use.cairo=T)

} #end of zoom loop

