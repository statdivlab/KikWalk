##########################################
######### Preamble
##########################################

## A script to reproduce Figures 2 & 3
## of "Lifetime Mobility of an Arctic Woolly Mammoth"
## by Wooller, Bataille, ... and Willis (2021, Science)

# Author: Clement Bataille and Amy Willis, 2020
# Software contact: Amy Willis, adwillis@uw.edu

# This script follows `isotope_guided_walk.R`. You will need the output of that script first!

##########################################
######### Figures
##########################################

# Calculate density of paths at each location of the map
dens <- MASS::kde2d(best_walk_full$x, best_walk_full$y, n = 400)
dens_df <- data.frame(with(dens, expand_grid(y, x)),
                      z = as.vector(dens$z))
fit <- loess(z ~ y * x, data = dens_df, span = 0.02)
best_walk_full$w <- predict(fit, with(best_walk_full, data.frame(x = x, y = y)))

# Create base map for walk on the land bridge for Fig. 2 main manuscript
blue.col <- colorRampPalette(c("darkblue", "lightblue"))
col2 <- c(blue.col(5), terrain.colors(5))
col2 <- colorRampPalette(col2)
breakpoints3 <- c(-5000,-2000,-200,0,100,300,500,1000,2000,4100)
lb <- ggplot(data=xyz_tib, aes(x, y)) +
  geom_raster(aes(fill = z), alpha=0.5) +
  scale_fill_gradientn(name = "Elevation (m)",colours = col2(20), breaks = breakpoints3)

# Split the dataset per phases
baby <- best_walk_full %>% filter(phase == 1)
juvenile <- best_walk_full %>% filter(phase == 2)
adolescent <- best_walk_full %>% filter(phase == 1 | phase == 2)
adult <- best_walk_full %>% filter(phase == 3)
death <- best_walk_full %>% filter(phase == 4)


# Add glacier polygons to figures
# Data available from http://instaar.colorado.edu/QGISL/ak_paleoglacier_atlas/downloads/index.html
# or in zip folder
glaciers_input <- readOGR("../data/glacier/", "glacier15-18")
glaciers <- crop(glaciers_input, r_elevation)
glaciers_df <- fortify(glaciers)


blue.col <- colorRampPalette(c("darkblue", "lightblue"))
grey.col<- colorRampPalette(c("#00A600FF", "#63C600FF", "#E6E600FF",
                              "#E8C32EFF", "#EBB25EFF", "#EDB48EFF", "#F0C9C0FF", "#F2F2F2FF"))
col2<-c(blue.col(10), grey.col(10))
col2<-colorRampPalette(col2)
breakpoints3<-c(-2000,0,1000,2500)

## Make Fig. 2
map <- lb + geom_path(data = adolescent, aes(x = x, y = y, group = walk, col = w), size = 0.2) +
  scale_color_gradientn(colors = c( "grey85", "grey75", "grey60","grey35","grey1"))+
  theme_void() + theme(legend.position = "none")

map_adolescent <- map + geom_point(aes(x = -152.562504, y = 69.684816), colour = "black", size = 2)+
  geom_polygon(data = glaciers_df, aes(x = long, y = lat, group = group),
               color = 'white', fill = 'white', size = .2)+
  scale_x_continuous(limits = c(-170, -140))+
  scale_y_continuous(limits = c(55, 75))

map <- lb +geom_path(data = adult,aes(x = x, y = y, group = walk,col = w), size = 0.2) +
  scale_color_gradientn(colors = c( "grey85", "grey75", "grey60","grey35","grey1"))+
  theme_void() + theme(legend.position = "none")

map_adult <- map + geom_point(aes(x = -152.562504, y = 69.684816), colour = "black", size = 2)+
  geom_polygon(data = glaciers_df, aes(x = long, y = lat, group = group),
               color = 'white', fill = 'white', size = .2)+
  scale_x_continuous(limits = c(-170, -140))+
  scale_y_continuous(limits = c(55, 75))

map <- lb + geom_path(data = death,aes(x = x, y = y, group = walk,col = w),size = 0.2) +
  scale_color_gradientn(colors = c( "grey85", "grey75", "grey60","grey35","grey1"))+
  theme_void() + theme(legend.position = "none")

map_death <- map + geom_point(aes(x = -152.562504, y = 69.684816), colour = "black",size = 2)+
  geom_polygon(data = glaciers_df, aes(x = long, y = lat, group = group),
               color = 'white', fill = 'white', size = .2)+
  scale_x_continuous(limits = c(-170, -140))+
  scale_y_continuous(limits = c(55, 75))

###Save map Fig.2
ggarrange(map_adolescent, map_adult, map_death,
          ncol = 1,nrow = 3,common.legend = TRUE)
ggsave("ggarranged_t1_walks_i3_100km_top.pdf",width = 3, height = 9)
# TODO ask Clem if Fig 2 used sd of 2 or 3?

##################### Create Polygons for Fig. 3 #########################

coords <- data.frame(best_walks_full$x, best_walks_full$y)
bestwalkstep_sp <- SpatialPointsDataFrame(coords, best_walks_full,
                                          proj4string = crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
baby_sp <- bestwalkstep_sp[which(bestwalkstep_sp$phase == 1), ]
juvenile_sp <- bestwalkstep_sp[which(bestwalkstep_sp$phase == 2), ]
adult_sp <- bestwalkstep_sp[which(bestwalkstep_sp$phase == 3), ]
death_sp <- bestwalkstep_sp[which(bestwalkstep_sp$phase == 4), ]

babysSp <- as(SpatialPoints(baby_sp), "ppp")
juvenilesSp <- as(SpatialPoints(juvenile_sp), "ppp")
adultsSp <- as(SpatialPoints(adult_sp), "ppp")
deathsSp <- as(SpatialPoints(death_sp), "ppp")

babyDens <- density(babysSp, adjust = 0.2)
juvenileDens <- density(juvenilesSp, adjust = 0.2)
adultDens <- density(adultsSp, adjust = 0.2)
deathDens <- density(deathsSp, adjust = 0.2)

babyDens <- raster(babyDens)
juvenileDens <- raster(juvenileDens)
adultDens <- raster(adultDens)
deathDens <- raster(deathDens)

babyDens <- calcPostProb(babyDens)
babyDens_prob <- calcNormProb(babyDens)
juvenileDens <- calcPostProb(juvenileDens)
juvenileDens_prob <- calcNormProb(juvenileDens)
adultDens <- calcPostProb(adultDens)
adultDens_prob <- calcNormProb(adultDens)
deathDens <- calcPostProb(deathDens)
deathDens_prob <- calcNormProb(deathDens)

baby_area <- babyDens_prob[babyDens_prob>0.44]
threshold_baby <- sum(baby_area)/cellStats(babyDens_prob,sum)
threshold_baby

juvenile_area <- juvenileDens_prob[juvenileDens_prob>0.38]
threshold_juvenile <- sum(juvenile_area)/cellStats(juvenileDens_prob,sum)
threshold_juvenile

adult_area <- adultDens_prob[adultDens_prob>0.38]
threshold_adult <- sum(adult_area)/cellStats(adultDens_prob,sum)
threshold_adult

death_area <- deathDens_prob[deathDens_prob>0.52]
threshold_death <- sum(death_area)/cellStats(deathDens_prob,sum)
threshold_death

babyDens_prob50 <- rasterToPolygons(babyDens_prob, fun=function(x){x>0.375}, dissolve=FALSE)
juvenileDens_prob50 <- rasterToPolygons(juvenileDens_prob, fun=function(x){x>0.5}, dissolve=FALSE)
adultDens_prob50 <- rasterToPolygons(adultDens_prob, fun=function(x){x>0.35}, dissolve=FALSE)
deathDens_prob50 <- rasterToPolygons(deathDens_prob, fun=function(x){x>0.52}, dissolve=FALSE)

# Output
# writeOGR(babyDens_prob50, dsn="polygon",layer = "babyDens_prob50_2",
#          driver = "ESRI Shapefile" ,overwrite=TRUE)
# writeOGR(juvenileDens_prob50, dsn="polygon",layer = "juvenileDens_prob50_2",
#          driver = "ESRI Shapefile",overwrite=TRUE )
# writeOGR(adultDens_prob50, dsn="polygon",layer = "adultDens_prob50_2",
#          driver = "ESRI Shapefile" ,overwrite=TRUE)
# writeOGR(deathDens_prob50, dsn="polygon",layer = "deathDens_prob50_2",
#          driver = "ESRI Shapefile",overwrite=TRUE )
