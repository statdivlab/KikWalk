##########################################
######### Preamble
##########################################

## A script to reproduce the isotope random walk model and results
## of "Lifetime Mobility of an Arctic Woolly Mammoth"
## by Wooller, Bataille, ... and Willis (2021, Science)

# Authors: Amy Willis and Clement Bataille, 2020
# Software contact: Amy Willis, adwillis@uw.edu

# This script assumes basic knowledge of R. If you have not used R
# before, I (Amy) have some teaching materials available here:
# https://github.com/adw96/biost509

# You should be running this script in an R Project with working directory
# set to the subfolder `R/`. This will mean that you do not have to modify any
# of the local file paths in the below!

##########################################
######### Preliminaries and data
##########################################

# Begin by loading in some needed packages. If you have not installed them,
# you will need to do that first using `install.packages("[package name]")`

library(tidyverse)
library(ggpubr)
library(magrittr)
library(rgdal)
library(raster)
library(maps)
library(maptools)
data(wrld_simpl)
library(rasterVis)
library(mvnmle)
library(mixtools)
library(fossil)
library(colorRamps)
library(colorspace)
library(ggplot2)
library(gridExtra)
library(broom)
library(data.table)
library(readxl)
library(sp)
library(spatstat)
library(rgeos)


devtools::install_version("assignR", version = "1.2.1", repos = "http://cran.us.r-project.org")
library(assignR) 
data(d18o_world)

# To run the analyses in parallel, you will need one of the following two options:

## Mac option: `parallel`
library(parallel)

## Windows option: `parallelsugar`
# devtools::install_github('nathanvan/parallelsugar') ###This library creates the mclapply function for windows users
# library(parallelsugar)

# Load some utility functions
source("utility-functions.R")

# Import strontium isotope data and oxygen isotope data for each tusk distance.
# This datasheet also includes the interpreted temporal model (year, season, phase).
isodata <- read_csv("../data/isotope_data.csv")

# Import elevation raster for Alaska and create a tibble
r_elevation <- raster("../data/elevation.tif")
xyz_tib <- r_elevation %>%
  rasterToPoints %>%
  as_tibble %>%
  rename(z = 3)

# Import glaciated region for LGM Alaska and create a tibble
mask_LGM1 <- raster("../data/LGM_mask.tif")
mask_tib <- mask_LGM1 %>%
  rasterToPoints %>%
  as_tibble %>%
  rename(mask = 3) %>%
  mutate(mask = (mask == 0))

# Join the tibbles:
xyz_mask_tib <- xyz_tib %>%
  inner_join(mask_tib, by = c("x", "y"))

# Import the Sr isoscape and create a tibble
srsr_median_new_crs <- raster("../data/rf_AK.tif")
sr_tib <- srsr_median_new_crs %>%
  rasterToPoints %>%
  as_tibble %>%
  rename(sr = 3)

# Join the tibbles
xyz_mask_sr_tib <- xyz_mask_tib %>%
  inner_join(sr_tib, by = c("x", "y"))

# Import the Sr isoscape uncertainty and create a tibble
srsr_sd_new_crs <- raster("../data/srsd.tif")
srsd_tib <- srsr_sd_new_crs %>%
  rasterToPoints %>%
  as_tibble %>%
  rename(srsd = 3)

xyz_mask_srsd_tib <- xyz_mask_sr_tib %>%
  inner_join(srsd_tib, by = c("x", "y"))

# Create oxygen isoscape: Import calibration dataset from Metcalfe et al. 2017
O_mam <- read_excel("../data/Metcalfe_database.xlsx", col_names=TRUE, na="NA", sheet="R") %>%
  dplyr::select(Longitude, Latitude, d18Opdb) %>%
  filter(!is.na(Longitude)) %>%
  filter(!is.na(Latitude)) %>%
  filter(!is.na(d18Opdb)) %>%
  as.data.frame
O_mam_proj <- O_mam[, c("Longitude","Latitude")]


# Convert lat/long data to a spatial point dataframe for use in assignR package
d <- SpatialPointsDataFrame(O_mam_proj, data.frame(O_mam$d18Opdb),
                            proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Use calRaster from assignR package to calibrate a LGM oxygen isoscape in enamel
r <- calRaster(known = d, isoscape = d18o_world)

# Create oxygen isoscape figure (basis for Fig. S16 in SI)
breakpoints2 <- seq(-25,5,2.5)
# pdf("../output/Fig_Ocalibration.pdf",width=11, height=5.5)
par(mfrow=c(1,1))
plot(d18o_world$mean, col=matlab.like(12), breaks=breakpoints2, axes = FALSE)
points(O_mam_proj,pch=10,col="black",cex=0.5)
# dev.off()


# Reproject d18O isoscape mean into AK projection and incorporate into tibble
O_new_crs <- projectRaster(from = r$isoscape.rescale$mean, to = r_elevation, method="ngb")
O_tib <- O_new_crs %>%
  rasterToPoints %>%
  as_tibble %>%
  rename(O = 3)

xyz_mask_srsd_o_tib <- xyz_mask_srsd_tib %>%
  inner_join(O_tib, by = c("x", "y"))

# Reproject d18O isoscape uncertainty into AK projection and incorporate into tibble
Osd_new_crs <- projectRaster(from = r$isoscape.rescale$sd, to = r_elevation, method="ngb")
Osd_tib <- Osd_new_crs %>%
  rasterToPoints %>%
  as_tibble %>%
  rename(Osd = 3)

xyz_mask_srsd_osd_tib <- xyz_mask_srsd_o_tib %>%
  inner_join(Osd_tib, by = c("x", "y"))

#  Final tibble incorporating all isotope data and masks
useable <- xyz_mask_srsd_osd_tib %>%
  filter(!mask) %>%
  dplyr::select(x,y,z,sr,srsd,O,Osd)

##########################################
######### Isotope Guided Random Walks
##########################################

# Initiate walk at location of death
start <- useable %>%
  mutate(dist = calcDist(lat1=y, long1 = x,
                         lat2=69.684816, long2=-152.56250)) %>%
  arrange(dist) %>%
  head(1)
random_walk_result1 <- start

# A function that generates a random walk for a given R seed
# inputs:
#    j: the randomisation seed
#    iso_threshold: the threshold for mismatch between the best available step
#                   and the tusk measurement. Kill the walk if the best possibility
#                   is more than 2 s.d. away.
run_walk_screened <- function(j, dist_threshold = 100, iso_threshold = 2) {

  # Initialise
  current <- start
  final <- start
  set.seed(j)
  i <- 2
  best_iso_match <- 0

  # Create the loop
  while (best_iso_match < iso_threshold & i <= nrow(isodata)) {

    # calculate distances and standardised scores; filter by distance
    scores <- useable %>%
      mutate(dist = calcDist(lat1=y, long1 = x, lat2=current$y, long2=current$x),
             vertical_dist = z - current$z) %>%
      filter(dist < dist_threshold,
             abs(vertical_dist) < 500) %>%
      mutate(iso_score = (isodata$Sr[i] - sr)/(srsd))

    # How good was the best isotope match? If it's bad, we hit a "dead end" walk.
    best_iso_match <- scores$iso_score %>% abs %>% min

    # calculate probabilities
    probs <- scores %>%
      mutate(prob = dnorm(x = iso_score)) %>%
      arrange(desc(prob))

    # randomly choose a walk, with higher probability of selecting a better isotope match
    # add it to the tibble
    choice <- sample(1:nrow(probs), size=1, prob=probs$prob)
    current <- probs[choice, ] %>%
      mutate(step = i, walk = j, best_option = best_iso_match)
    final %<>%
      bind_rows(current)

    i <- i + 1
  }

  # Return the completed (or dead-end) walk
  final
}

# check the function and understandard its output
# warnings are expected and okay
test_walk <- run_walk_screened(30, dist_threshold = 100, iso_threshold = 2) # test: random number 30, 2 sd for isotope match thresholding, 100km radius

# You might consider estimating how long
# it could take you to run many walks:
system.time({test_walks <- mclapply(X = 1:70,
                                    FUN = run_walk_screened,
                                    dist_threshold = 100,
                                    iso_threshold = 2, mc.cores = 7) })

# The following step will take a long time!
# It runs 20,000 walks!
# We recommend parallelising it over your system, and perhaps
# doing it in chunks.
# Nevertheless, here we go...
# It runs in ~9 hours on my Mac (4.2 GHz Quad-Core Intel Core i7; 32 GB 2400 MHz DDR4)
system.time({all_walks <- mclapply(X = 1:20000,
                                   FUN = run_walk_screened,
                                   dist_threshold = 100,
                                   iso_threshold = 2,
                                   mc.cores = 7) })
all_walks_tb <- all_walks %>% bind_rows
# I recommend you save them! e.g., with
write_csv(x=all_walks_tb, path="../output/walks_20000_d_100_iso_2.csv")
# all_walks_tb <- read_csv("../output/walks_d_100_iso_2.csv")

##########################################
######### Analysis of the walks
##########################################

# How far does he travel over the course of his life?
# Also, join the isotope data to the walk tibble
# and use a moving average to average O and d
cumulative_all <- all_walks_tb %>%
  filter(!is.na(prob)) %>% # remove step 1
  group_by(walk) %>%
  mutate(cumdist = cumsum(dist)) %>%
  inner_join(isodata %>%
               mutate(step = 1:nrow(isodata))) %>%
  mutate(maO = frollmean(O, 10, align=c("center"))) %>%
  mutate(mad = frollmean(d, 10, align=c("center")))

## Which walks did not lead to isotope mismatches?
full_length_walks <- cumulative_all %>%
  filter(step == 1133) %>%
  pull(walk)

# how many walks go to full term?
full_length_walks %>% length
20000-(full_length_walks %>% length)

# when do walks get killed?
cumulative_all %>%
  summarise(steps = max(step)) %>%
  pull(steps) %>% hist

cumulative_all %>%
  filter(step == 1133) %>%
  pull(cumdist) %>%
  hist(main = "Cumulative distance (km) travelled over lifetime")
# still pretty consistent range of steps

##########################################
######### Rank walks using match between the tusk and walk d18O measurements:
##########################################
## Fit a linear model and pull out the R^2, slope, etc.
lms <- cumulative_all %>%
  filter(walk %in% full_length_walks) %>%
  split(.$walk) %>%
  purrr::map(~lm(maO ~ Oob, data = .))
rsqs <- lms %>%
  purrr::map(summary) %>%
  map_dbl(~.$r.squared)
sses <- lms %>%
  purrr::map(deviance) %>%
  unlist
slopes <- lms %>%
  purrr::map(summary) %>%
  map_dbl(~.$coefficients[2])
slope_se <- lms %>%
  purrr::map(summary) %>%
  map_dbl(~.$coef[2,2])
intercepts <- lms %>%
  purrr::map(summary) %>%
  map_dbl(~.$coefficients[1])

# Combine this data into a new tibble, discard walks with negative correlation in d18O
walk_lms <- tibble("walk" = full_length_walks,
                   slopes, slope_se, intercepts,
                   rsqs, sses) %>%
  filter(slopes > 0) %>%
  arrange(desc(rsqs))

# How many left with R^2 more than 0.1?
best_walks <- walk_lms %>% filter(rsqs > 0.1) %>%
  arrange(desc(rsqs)) %>% pull(walk)
best_walks %>% length
# Which is the best?
best_walk <- walk_lms %>% filter(rsqs > 0.1) %>%
  arrange(desc(rsqs)) %>% pull(walk) %>% head(1)


# Obtain the full walk for the single and multiple best walks:
best_walks_full <- cumulative_all %>% filter(walk %in% best_walks)
best_walk_full <- cumulative_all %>% filter(walk %in% best_walk)
