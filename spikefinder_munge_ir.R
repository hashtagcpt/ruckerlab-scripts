rm(list = ls())

### LOAD LIBRARIES ####
library(xlsx)
library(plyr)

# note - this code requires:
#   1) a specific format of the spikefider output
#   2) formatted refraction data from the summary sheet
#   3) changes to file names in 3 places, lenstar and refraction input, and output

# 1. change working directory -- change for each experiment
setwd("/media/ruckerlab/tiny/Lenstar/")

#### INPUT/OUTPUT FILE DEFINITIONS ####

# read in SpikesFound.csv -- data file from SpikeFinder.exe

# 2. data file from spikefinder marking to munge
d <- read.csv('falk_dec_17_to_munge_clean.csv') 

# 3. name of post-munge data file for graphing & stats analysis
# output file name
file_csv <- paste('falk_dec_17_for_analysis.csv')

# 4. bird condition file -- CHANGE FOR EACH EXP 
# set up bird numbers and conditions for a given experiment
dBird <- read.csv('falk_bird_dec_17_conditions_clean.csv')
Bird <- dBird$Bird
Light <- dBird$light
Flicker <- dBird$flicker
Surgery <- dBird$surgery

#### AVERAGE OR DIFFERENCE EYES ####
avgEyes <- TRUE

# delete first row/col of junk - if present, assumed to be present
d <- d[-1,]
d <- d[,-1]

d<-arrange(d, d$X) 

# cast colums to numeric types
# Cornea/Aqueous (X Coordinate)
d$X.29 <- as.numeric(as.character(d$X.29))

# Aqueous/Lens Baseline (X Coordinate)
d$X.31 <- as.numeric(as.character(d$X.31)) 

# Aqueous/Capsule (X Coordinate)
d$X.33 <- as.numeric(as.character(d$X.33))

# Anterior Nucleus Nadir (X Coordinate)
d$X.39 <- as.numeric(as.character(d$X.39))

# Posterior Nucleus Nadir (X Coordinate)
d$X.41 <- as.numeric(as.character(d$X.41))

# Capsule/Vitreous (X Coordinate)
d$X.47 <- as.numeric(as.character(d$X.47))

# Lens/Vitreous Baseline (X Coordinate)
d$X.49 <- as.numeric(as.character(d$X.49))

# Retina (X Coordinate)
d$X.51 <- as.numeric(as.character(d$X.51))

# RPE (X Coordinate)
d$X.52 <- as.numeric(as.character(d$X.52))

# Total OPL in pixels (RPE X - 1000)
d$X.54 <- as.numeric(as.character(d$X.54))

### EYE ASSIGNMENT #### 

# subset to ex eye, X
dX <- subset(d, d$X.3 == 'OS')
# subset to fellow eye, N
dN <- subset(d, d$X.3 == 'OD')

#### CONSTANTS ####

# important constants from RI sheet
RI_Lenstar_Wavelength <- 827.16
RI_Cornea <- 1.34507504
RI_Aqueous <- 1.344191678
RI_Lens	 <-1.419086426
RI_Vitreous	<- 1.342799941
RI_Retina <- 1.36
RI_Eye_Length <- 1.357442592
RI_Choroid <- 1.36
RI_Sclera	<- 1.36

# the 1250 constant?
tf <- 1250

#### COMPUTE X EYE VALUES ####
# OD - X variables to use
# Xc1 Xc2 Xl1 Xl2	
# Xv1Xv2 Xr1 Xr2	
# Xch1 Xch2	 Xsc1 Xsc2

# X.29 - Cornea / Aqueous X-coord
Xc <- (dX$X.29-1000)/tf/RI_Cornea
Xaq <- (dX$X.33-dX$X.29)/tf/(RI_Aqueous-0.1)
Xc <- Xc + Xaq
Xc1 <- Xc[seq(1, length(Xc), 2)]
Xc2 <- Xc[seq(2, length(Xc), 2)]

# X.47 Capsule/Vitreous (X Coordinate)
Xl <- (dX$X.47-dX$X.33)/tf/RI_Lens
Xl1 <- Xl[seq(1, length(Xl), 2)]
Xl2 <- Xl[seq(2, length(Xl), 2)]

# X.51 Retina 
Xv <- (dX$X.51-dX$X.47)/tf/RI_Vitreous
Xv1 <- Xv[seq(1, length(Xv), 2)]
Xv2 <- Xv[seq(2, length(Xv), 2)]

# X.51 RPE (X Coordinate)
Xr <- (dX$X.52-dX$X.51)/tf/RI_Retina
Xr1 <- Xr[seq(1, length(Xr), 2)]
Xr2 <- Xr[seq(2, length(Xr), 2)]

# X.39 Anterior Nucleus Nadir (X Coordinate)
Xch <- (dX$X.39-dX$X.52)/tf/RI_Choroid
Xch1 <- Xch[seq(1, length(Xch), 2)]
Xch2 <- Xch[seq(2, length(Xch), 2)]

# X.41 Posterior Nucleus Nadir (X Coordinate)
Xsc <- (dX$X.41-dX$X.39)/tf/RI_Sclera
Xsc1 <- Xsc[seq(1, length(Xsc), 2)]
Xsc2 <- Xsc[seq(2, length(Xsc), 2)]

#### COMPUTE N EYE VALUES ####

# OS - N variables to use
# Nc1 Nc2 Nl1 Nl2	
# Nv1 Nv2 Nr1 Nr2	
# Nch1 Nch2	 Nsc1 Nsc2

# X.29 - Cornea / Aqueous X-coord, X.33 - Aqueous/Capsule (X Coordinate)
Nc <- (dN$X.29-1000)/tf/RI_Cornea
Naq <- (dX$X.33-dX$X.29)/tf/(RI_Aqueous-0.1)
Nc <- Nc + Naq
Nc1 <- Nc[seq(1, length(Nc), 2)]
Nc2 <- Nc[seq(2, length(Nc), 2)]

# X.47 Capsule/Vitreous (X Coordinate)
Nl <- (dN$X.47-dN$X.33)/tf/RI_Lens
Nl1 <- Nl[seq(1, length(Nl), 2)]
Nl2 <- Nl[seq(2, length(Nl), 2)]

# X.50 Retina 
Nv <- (dN$X.51-dN$X.47)/tf/RI_Vitreous
Nv1 <- Nv[seq(1, length(Nv), 2)]
Nv2 <- Nv[seq(2, length(Nv), 2)]

# X.51 RPE (X Coordinate)
Nr <- (dN$X.52-dN$X.51)/tf/RI_Retina
Nr1 <- Nr[seq(1, length(Nr), 2)]
Nr2 <- Nr[seq(2, length(Nr), 2)]

# X.39 Anterior Nucleus Nadir (X Coordinate)
Nch <- (dN$X.39-dN$X.52)/tf/RI_Choroid
Nch1 <- Nch[seq(1, length(Nch), 2)]
Nch2 <- Nch[seq(2, length(Nch), 2)]

# X.41 Posterior Nucleus Nadir (X Coordinate)
Nsc <- (dN$X.41-dN$X.39)/tf/RI_Sclera
Nsc1 <- Nsc[seq(1, length(Nsc), 2)]
Nsc2 <- Nsc[seq(2, length(Nsc), 2)]

Bird <- unique(d$X)
# set experiment name
exp_name <- rep(' ',length(Bird))

# set condition
#condition <- rep(' ',length(Bird))

#### FIND EYE LENGTH #### 
Xax1 <- Xc1 + Xl1 + Xv1 + Xr1 + Xch1 + Xsc1
Xax2 <- Xc2 + Xl2 + Xv2 + Xr2 + Xch2 + Xsc2
Nax1 <- Nc1 + Nl1 + Nv1 + Nr1 + Nch1 + Nsc1
Nax2 <- Nc2 + Nl2 + Nv2 + Nr2 + Nch2 + Nsc2

#### COMPUTE DIFFERENCES ####
Xcacdif <- Xc2 - Xc1
Ncacdif <- Nc2 - Nc1

Xlensdif <- Xl2 - Xl1
Nlensdif <- Nl2 - Nl1

Xvitdif <- Xv2 - Xv1
Nvitdif <- Nv2 - Nv1

Xchordif <- Xch2 - Xch1
Nchordif <- Nch2 - Nch1

Xretdif <- Xr2 - Xr1
Nretdif <- Nr2 - Nr1

Xscldif	<- Xsc2 - Xsc1
Nscldif	<- Nsc2 - Nsc1

Xaxdif	<- Xax2 - Xax1
Naxdif <- Nax2 - Nax1

#### X - N COMPUTE (FOR MONOCULAR EXP) ####
cacXmN <- Xcacdif - Ncacdif
lensXmN <- Xlensdif - Nlensdif	
vitXmN	<- Xvitdif - Nvitdif
retXmN	<- Xretdif - Nretdif
chorXmN <- Xchordif - Nchordif	
sclXmN	<- Xscldif - Nchordif
axXmN <- Xaxdif - Naxdif

#### MEAN EYE CHANGE (FOR BINOCULAR EXP) ####
Mean_cac_diff <- (Xcacdif+Ncacdif)/2	
Mean_lens_diff <- (Xlensdif+Nlensdif)/2
Mean_vit_diff <- (Xvitdif+Nvitdif)/2	
Mean_ret_diff <- (Xretdif+Nretdif)/2	
Mean_chor_diff <- (Xchordif+Nchordif)/2	
Mean_scl_diff <- (Xscldif+Nscldif)/2	
Mean_ax_diff <- (Xaxdif+Naxdif)/2	

#### DATA FRAME TO OUTPUT ####

mydata <- data.frame(Bird,Light,Flicker,Surgery,Xc1,Xc2,Xl1,Xl2,Xv1,Xv2,Xr1,Xr2,Xch1,Xch2,Xsc1,Xsc2,Nc1,Nc2,Nl1,Nl2,Nv1,Nv2,Nr1,Nr2,Nch1,Nch2,Nsc1,Nsc2,Xax1, Xax2, Nax1, Nax2, Xcacdif,Ncacdif,Xlensdif, Nlensdif, Xvitdif, Nvitdif, Xchordif, Nchordif, Xretdif,Nretdif, Xscldif, Nscldif,Xaxdif, Naxdif, Mean_cac_diff, Mean_lens_diff, Mean_vit_diff, Mean_ret_diff, Mean_chor_diff, Mean_scl_diff, Mean_ax_diff)

#### DATA FILE OUTPUT #### 
write.csv(mydata,file_csv)