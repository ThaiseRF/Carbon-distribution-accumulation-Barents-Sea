# Drivers of organic carbon distribution and accumulation in the northern Barents Sea
# Adapted from 'Numerical Ecology with R' Borcard, Gillet and Legendre 2018
# Thaise Ricardo de Freitas
# 20 December 2022

# load packages -----------------------------

if(!require(pacman)) install.packages("pacman")
pacman::p_load(readxl, ggplot2, tidyverse, tidypaleo, mgcv, khroma, vegan, 
               sf, ecodist, spdep, adespatial)

# dependent variables -----------

## order and summarize data (average top 6 cm) -------
order  <- c("NPAL4", "NPAL5", "NPAL7",
            "NPAL8", "NPAL12", "NPAL14",
            "NPAL15", "NPAL17", "NPAL19",
            "P1",  "P2",  "P4",  "P6",  "P7","SICE4", "A", "B")

order_camp  <- c("Paleo", "SSQ1", "SSQ2", "SSQ3", "SSQ4", "EISA")

dep <- abiot_total %>%
  filter((slice %in% c('0-1', '1-2', '2-3', '3-4', '4-5', '5-6', '2-4', '4-6'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  filter(!(station == 'P6' & replicate == 'r2' & campaign == 'SSQ1')) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  group_by(campaign, station, group, replicate) %>%
  summarise_at(vars(tn, tc, d15, d13, toc, tic, TOM, MOM, CN_org, toc_63), mean, na.rm = TRUE) %>%
  mutate(campaign =  factor(campaign, levels = order_camp)) %>%
  mutate(station =  factor(station, levels = order)) %>%
  arrange(campaign) %>%
  arrange(station)

## selec variables -----------

dep_f <- dep %>%
  ungroup() %>%
  select(tn, tc, d15, d13, tic, CN_org, toc_63)

cor(dep_f)
plot(dep_f)

# spatial coordinates -------------
coord <- abiot_total %>%
  filter((slice %in% c('0-1'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  filter(!(station == 'P6' & replicate == 'r2' & campaign == 'SSQ1')) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  mutate(campaign =  factor(campaign, levels = order_camp)) %>%
  mutate(station =  factor(station, levels = order)) %>%
  arrange(campaign) %>%
  arrange(station) %>%
  select(lat, long)


## coordinates projection PCNM -----------
coord_index_or <- abiot_total %>%
  filter((slice %in% c('0-1'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  filter(!(station == 'P6' & replicate == 'r2' & campaign == 'SSQ1')) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  mutate(campaign =  factor(campaign, levels = order_camp)) %>%
  mutate(station =  factor(station, levels = order)) %>%
  arrange(campaign) %>%
  arrange(station) %>%
  #arrange(lat) %>%
  #slice(-c(4, 5)) %>% 
  st_as_sf(coords = c("long", "lat")) %>%
  st_set_crs( 4326 ) %>% 
  st_transform( 3996 ) 

ggplot() + 
  geom_sf(data = coord_index_or, aes(shape = station)) + # Add coordinate data
  scale_shape_manual(values = c(10, 20, 2, 6, 15, 18, 22, 3, 23, 4, 1, 7, 8, 9, 11, 65, 66))+
  theme_bw() + 
  xlab("Longitude") + # Change x axis title
  ylab("Latitude") +
  coord_sf()

coord_index <- coord_index_or %>%
  mutate(lat = unlist(map(coord_index_or$geometry,2)),
         long = unlist(map(coord_index_or$geometry,1))) %>% 
  select(long, lat, station) %>%
  st_set_geometry(NULL)

ggplot() + 
  geom_point(data = coord_index, aes(x = long, y = lat, shape = station)) + # Add coordinate data
  scale_shape_manual(values = c(10, 20, 2, 6, 15, 18, 22, 3, 23, 4, 1, 7, 8, 9, 11, 65, 66))+
  theme_bw() + 
  xlab("Longitude") + # Change x axis title
  ylab("Latitude")

coord_index <- coord_index %>%
  select(long, lat) %>%
  #arrange(lat) %>%
  mutate(index = 1:n())

## UTM spatial coordinates  -----------
coord_utm_org <- read_xlsx("coord_utm.xlsx")

coord_utm <- coord_utm_org %>%
  select(X, Y) # x is latitude-ish, y is longitude-ish


# environmental variables --------------

env_f <- env_total %>%
  filter(!(station == 'P6' & # remove problematic replicate
             replicate == 'r2' & 
             campaign == 'SSQ1')) %>%
  select(-(1:9), # selecnt only env variables
         prof, 
         land_distance, 
         clay,
         silt,
         sand,
         -mean_grain, # EISA cruise mean grain size not available
         SBT, 
         SST,
         SSS,
         SBS,
         SSO,
         SBO,
         chla_av,
         -ice_mean_spring, # ice spring and fall removed equal to summer/fall
         ice_mean_summer,
         -ice_mean_fall,
         ice_mean_winter
  )


# Weighted PCNM for RDA ---------------

rs <- rowSums(dep_f)/sum(dep_f)
pcnmw <- pcnm(dist(coord_utm), w = rs)

#threshold distance for truncation
pcnmw$threshold
summary(pcnmw)
summary(eigenvals(pcnmw))

plot(coord_utm[,"X"] ~ coord_utm[,"Y"],
     pch=15, cex=1.5,
     col=gray.colors(12)[cut(pcnmw$vectors[,1], breaks = 12)])

plot(coord_utm[,"X"] ~ coord_utm[,"Y"],
     pch=15, cex=1.5,
     col=gray.colors(12)[cut(pcnmw$vectors[,11], breaks = 12)])

pcnmw <- as.data.frame(scores(pcnmw))

#make a neighborhood list
coord_utm2 <- coord_utm
coord_utm2 <- st_as_sf(coord_utm2, coords = c("X", "Y"), remove = FALSE)
neigh <- dnearneigh(x= coord_utm2, 
                    d1=0, d2=3)#d1 is minimum distance, d2 is max distance

#plot the neighorhood
plot(neigh,coordinates(coord_utm))

#create weights for the neighbors
wts <- nb2listw(neighbours=neigh, style='W', zero.policy=T) #W = row-standardized weights

#Moran's I test with normal approximation versus Monte Carlo permutation test
#moran.mc(x=pcnmw$PCNM10, listw=wts, nsim=999, zero.policy=T) #Monte Carlo
moran.mc(x=pcnmw$PCNM1, listw=wts, nsim=999, zero.policy=T)#normal approximation
moran.mc(x=pcnmw$PCNM2, listw=wts, nsim=999, zero.policy=T)#normal approximation
moran.mc(x=pcnmw$PCNM3, listw=wts, nsim=999, zero.policy=T)#normal approximation
moran.mc(x=pcnmw$PCNM4, listw=wts, nsim=999, zero.policy=T)#normal approximation
moran.mc(x=pcnmw$PCNM5, listw=wts, nsim=999, zero.policy=T)#normal approximation
moran.mc(x=pcnmw$PCNM6, listw=wts, nsim=999, zero.policy=T)#normal approximation
moran.mc(x=pcnmw$PCNM7, listw=wts, nsim=999, zero.policy=T)#normal approximation
moran.mc(x=pcnmw$PCNM8, listw=wts, nsim=999, zero.policy=T)#normal approximation
moran.mc(x=pcnmw$PCNM9, listw=wts, nsim=999, zero.policy=T)#normal approximation
moran.mc(x=pcnmw$PCNM10, listw=wts, nsim=999, zero.policy=T)#normal approximation
moran.mc(x=pcnmw$PCNM11, listw=wts, nsim=999, zero.policy=T)#normal approximation


ord_spa <- rda(dep_f ~., pcnmw)

ord_spa$anova
ord_spa$call


summary(ord_spa)
msoplot(mso(ord_spa, coord))
ordisplom(pcnmw, choices=1:4)

plot(pcnmw$PCNM1,type="l",ylab=c("PCNM1"))
plot(pcnmw$PCNM2,type="l",ylab=c("PCNM2"))
plot(pcnmw$PCNM3,type="l",ylab=c("PCNM3"))
plot(pcnmw$PCNM4,type="l",ylab=c("PCNM4"))
plot(pcnmw$PCNM8,type="l",ylab=c("PCNM8"))
plot(pcnmw$PCNM9,type="l",ylab=c("PCNM9"))
plot(pcnmw$PCNM11,type="l",ylab=c("PCNM11"))



# db-MEM --------------------

dbmemx <- dbmem(coord_utm, silent = FALSE, 
                MEM.autocor="positive", 
                store.listw = TRUE)
dbmemt <- as.data.frame(dbmemx)
attributes(dbmemx)$values
attr(dbmemx, "values")
length(attributes(dbmemx)$values)

par(mfrow = c(5, 1))
somedbmem <- c(1, 2, 3, 4, 5)
for(i in 1:length(somedbmem)){
  plot(dbmemt[ ,somedbmem[i]],
       type = "l",
       xlab = "X coordinate",
       ylab = c("dbMEM", somedbmem[i]))
}


## Step 2. Run the global dbMEM analysis on the detrended 

## Is there a linear trend in the data? -------
anova(rda(dep_f, coord_utm)) # Result: significant trend
# Computation of linearly detrended mite data
dep_det <- resid(lm(as.matrix(dep_f) ~ ., data = coord_utm))


(dbmem.rda <- rda(dep_det ~., dbmemx))
summary(dbmem.rda)
anova(dbmem.rda)

## Step 3. Since the R-square is significant, compute the adjusted

## R2 and run a forward selection of the dbmem variables ------
#(R2a <- RsquareAdj(dbmem.rda)$adj.r.squared)
#(dbmem.fwd <- forward.sel(dep_det,  as.matrix(dbmemt), adjR2thresh = R2a))

rda0_pc <- rda(dep_det ~1, dbmemx) 
# only one pcoa. only one intercept
rda1_pc <- rda(dep_det ~., dbmemx) 

dbmem.fwd <-
  ordiR2step(rda0_pc,
             scope = formula(rda1_pc),
             R2scope = F, 
             permutations = how(nperm = 999)
  )

summary(dbmem.fwd)
dbmem.fwd$anova
dbmem.fwd$call
anova(dbmem.fwd, by="terms", permu=9999)
plot(dbmem.fwd)

#(nb.sig.dbmem <- nrow(dbmem.fwd)) # Number of signif. dbMEM

spatial <- dbmemt %>% 
  select(MEM2, MEM3, MEM4, MEM5, MEM1)

par(mfrow = c(1, 5))
ade4::s.value(coord_index, spatial[,1]) #2
ade4::s.value(coord_index, spatial[,2]) #3
ade4::s.value(coord_index, spatial[,3]) #4
ade4::s.value(coord_index, spatial[,4]) #5
ade4::s.value(coord_index, spatial[,5]) #1
#ggsave(file="svalue_dbmem.pdf", width = 220, height = 180, units = "mm", dpi = "print")
dev.off()

# 4. Arbitrarily split the significant dbMEM into broad and fine scale
# Broad scale: dbMEM 1, 2, 3
dbmem.broad <- select(spatial, MEM2, MEM3, MEM1)
# Fine scale: dbMEM 4, 5 
dbmem.fine <- select(spatial, MEM4, MEM5)


## Step 4. New dbMEM analysis with 5 significant dbMEM variables

## Adjusted R-square after forward selection: R2adj = 0.2418
(dbmem.rda2 <- rda(dep_det ~ ., data = spatial))
(R2a2 <- RsquareAdj(dbmem.rda2)$adj.r.squared)
anova(dbmem.rda2)
(axes.test <- anova(dbmem.rda2, by = "axis"))

# Number of significant axes
(nb.ax <- length(which(axes.test[ , ncol(axes.test)] <= 0.05)))

# environmental data -------------
# Forward selection of the environmental variables - untrended data
rda0 <- rda(dep_f ~1, env_f)
rda1 <- rda(dep_f ~., env_f)

plot(rda1)

rdaplotp <- ordistep(rda0, R2scope = F, 
                     scope = formula(rda1), 
                     permutations = 999)

rdaplotp$anova


anova(rdaplotp, by="terms", permu=9999)

plot(rdaplotp)

RsquareAdj(rdaplotp)

env_final <- env_f %>% 
  select(ice_mean_winter, sand, chla_av, prof, silt, SBT, clay)

rda0_env <- rda(dep_f ~ ., data = env_final)

RsquareAdj(rda0_env)
anova(rda0_env, by="terms", permu=200)
(axes.test <- anova(rda0_env, by = "axis"))

## Map of ENVs and MEMs in the sample plot
op <- par(mfrow=c(4,3))
ordisurf(coord_utm, scores(dep_f$TIC), bubble = 4, 
         main = "Ice concentration - WINTER", 
         xlab = "Longitude", ylab = "Latitude")
ordisurf(coord_utm, scores(env_final$sand), bubble = 4, main = "Sand (%)", 
         xlab = "Longitude", ylab = "Latitude")
ordisurf(coord_utm, scores(env_final$chla_av), bubble = 4, main = "Chla average", 
         xlab = "Longitude", ylab = "Latitude")
ordisurf(coord_utm, scores(env_final$prof), bubble = 4, main = "Water depth", 
         xlab = "Longitude", ylab = "Latitude")
ordisurf(coord_utm, scores(env_final$silt), bubble = 4, main = "Silt (%)", 
         xlab = "Longitude", ylab = "Latitude")
ordisurf(coord_utm, scores(env_final$SBT), bubble = 4, main = "SBT", 
         xlab = "Longitude", ylab = "Latitude")
ordisurf(coord_utm, scores(env_final$clay), bubble = 4, main = "Clay", 
         xlab = "Longitude", ylab = "Latitude")
ordisurf(coord_utm, spatial$MEM1, bubble = 4, main = "dbMEM 1", 
         xlab = "Longitude", ylab = "Latitude")
ordisurf(coord_utm, spatial$MEM2, bubble = 4, main = "dbMEM 2", 
         xlab = "Longitude", ylab = "Latitude")
ordisurf(coord_utm, spatial$MEM3, bubble = 4, main = "dbMEM 3", 
         xlab = "Longitude", ylab = "Latitude")
ordisurf(coord_utm, spatial$MEM4, bubble = 4, main = "dbMEM 4", 
         xlab = "Longitude", ylab = "Latitude")
ordisurf(coord_utm, spatial$MEM5, bubble = 4, main = "dbMEM 5", 
         xlab = "Longitude", ylab = "Latitude")
par(op)


# Final RDA --------------

ord0 <- rda(dep_f ~ 1)

ord <- rda(dep_f ~ 
             spatial$MEM2+
             spatial$MEM3+
             spatial$MEM4+
             spatial$MEM5+
             spatial$MEM1+
             env_final$ice_mean_winter+
             env_final$sand +
             env_final$chla_av+
             env_final$prof+
             env_final$silt+
             env_final$SBT+
             env_final$clay)

summary(ord)
RsquareAdj(ord)
plot(ord)

plot(ord, choices = c(1,4))

anova(ord, by="terms", permu=9999)
anova(ord, by="margin", permu=9999)
anova(ord, by="axis", permu=9999)


bg <- c("#000000", "#E69F00" ,"#56B4E9", "#009E73")

dep[dep=="a"] <- '#000000'
dep[dep=="b"] <- '#E69F00'
dep[dep=="c"] <- '#56B4E9'
dep[dep=="d"] <- '#009E73'

labl <- c("MEM2", "MEM3", "MEM4", "MEM5", "MEM1", 
          "Winter sea ice",  "Sand", "Chla",
          "Water depth", "Silt",
          "SBT", "Clay")

colnames(dep_f)<-c('TN','TC', 'd15N','d13C', 'TIC', 'Corg/N', 'TOC63')

# PLOT -----------
plot(ord)

pdf(file="fig7_27-2.pdf", width = 5.54, height = 5.93)
plot(ord, type="n",xlab="dbRDA1 (39.7 %)",ylab="dbRDA2 (24.6 %)")
points(ord, cex= 1.5, 
       col = c("#882E72", #10, 
               "#1965B0", #20, 
               "#5289C7", #2,  
               "#7BAFDE", #6, 
               "#4EB265", #15,
               "#CAE0AB", #18, 
               "#F7F056", #22, 
               "#DC050C", #3,
               "#E8601C", #23,
               "#AA6F9E", #4, 
               "#994F88", #1, 
               "#90C987", #7, 
               "#F6C141", #8, 
               "#F1932D", #9, 
               "#72190E", #11
               "#D1BBD7", #65, 
               "#BA8DB4" #66
       )[dep$station], scaling = 3, 
       pch=c(10, 20, 2, 6, 15, 18, 22, 3, 23, 4, 1, 7, 8, 9, 11, 65, 66)[dep$station])

text(ord,"species",cex=1, col="gray20", scaling = 2)
#text(ord, labels=dep$station, cex=0.6, col=c("black"))
text(ord,  display="bp", col="gray50", cex=1, labels = labl, scaling = 3)  
legend("bottomleft", legend=levels(dep$station), bty="n",
       pch=c(10, 20, 2, 6, 15, 18, 22, 3, 23, 4, 1, 7, 8, 9, 11, 65, 66))



dev.off()

# 5. environment - trend - dbMEM variation partitioning

(var <-
    varpart(dep_f, env_final,spatial))
plot(var)
# Show the symbols of the fractions and plot their values
par(mfrow = c(1,2))
showvarparts(4, bg = c("red", "blue", "yellow", "green"))
plot(var,
     digits = 2,
     bg = c("red", "blue", "yellow", "green")
)


var <- varpart(dep_f, spatial, env_final)
var

plot(var)

# DEP + ENV DATAFRAME ------------------

dep <- abiot_total %>%
  filter((slice %in% c('0-1', '1-2', '2-3', '3-4', '4-5', '5-6', '2-4', '4-6'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  filter(!(station == 'P6' & replicate == 'r2' & campaign == 'SSQ1')) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  group_by(campaign, station, group, replicate) %>%
  summarise_at(vars(tn, tc, d15, d13, toc, tic, TOM, MOM, CN_org, toc_63), mean, na.rm = TRUE) %>%
  mutate(campaign =  factor(campaign, levels = order_camp)) %>%
  mutate(station =  factor(station, levels = order)) %>%
  arrange(campaign) %>%
  arrange(station)

dep_f <- dep %>%
  ungroup() %>%
  select(tn, tc, d15, d13, tic, CN_org, toc_63)

env_f <- env_total %>%
  filter(!(station == 'P6' & # remove problematic replicate
             replicate == 'r2' & 
             campaign == 'SSQ1')) %>%
  select((1:9), 
         -group, 
         -slice,
         -depth,
         land_distance, 
         clay,
         silt,
         sand,
         -mean_grain, # EISA cruise mean grain size not available
         SBT, 
         SST,
         SSS,
         SBS,
         SSO,
         SBO,
         chla_av,
         -ice_mean_spring, # ice spring and fall removed equal to summer/fall
         ice_mean_summer,
         -ice_mean_fall,
         ice_mean_winter
  )

## bind all data ----------
rda_final <- bind_cols(env_f,dep_f, spatial)

write.table(rda_final,
            file="rda_data_final.csv", dec = ",",
            append=FALSE, sep= ";", row.names = FALSE, col.names=TRUE)

