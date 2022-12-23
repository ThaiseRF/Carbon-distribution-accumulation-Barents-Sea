# Drivers of organic carbon distribution and accumulation in the northern Barents Sea
# Thaise Ricardo de Freitas
# 20 December 2022

# load packages -----------------------------
if(!require(pacman)) install.packages("pacman")
pacman::p_load(RNetCDF, writexl, readxl, raster, sp, lubridate, lattice,
               RCurl, ncdf4, sf, oce, stringi, tidyverse, RColorBrewer)

# Attributes data ---------------

##ICE 1850-2017 -----------------------------

#open your raster with n bands
raster_1850_2017 <- stack("C:/Users/trf/OneDrive - Universitetet i Oslo/PhD/maps/attributes/sea_ice/G10010_SIBT1850_V2/G10010_sibt1850_v2.0.nc")

#Plot it just to see if everything is ok
plot(raster_1850_2017$X2017.03.16)
points(abiot_sf)

#Check the number of bands
nlayers(raster)

#for(i in 1:nlayers(raster)){
#  band<-raster[[i]]
#  #save raster in a separate file
#  writeRaster(band, paste(names(raster[[i]]),'.tif', sep=""))
#}



### extract values -------------------

abiot_sf <- filter(abiot_total, slice == '0-1')

coordinates(abiot_sf) <- ~ long + lat

rasValue_1850_2017 = extract(raster_1850_2017, abiot_sf)
abiot_sf <- filter(abiot_total, slice == '0-1')
combinePointValue_1850_2017 = cbind(abiot_sf,rasValue_1850_2017)

write.table(combinePointValue_1850_2017,file="combinePointValue_1850_2017.csv", 
            append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)


## ICE 1978 - 2020 -------------------
setwd("C:/Users/trf/OneDrive - Universitetet i Oslo/PhD/maps/attributes/sea_ice/sea_ice_indexV3/monthly_concentration_1979_2020")

## download data

ice_dates = c("01_Jan", "02_Feb", "03_Mar", "04_Apr", "05_May", "06_Jun",
              "07_Jul", "08_Aug", "09_Sep", "10_Oct", "11_Nov", "12_Dec")

sea_url = "ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/north/monthly/geotiff/"

#Load all years of data, for each month
for(i in 1:12) {
  sea_url_single = glue::glue("{sea_url}{ice_dates[i]}/", )
  filenames = getURL(sea_url_single, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  filenames = strsplit(filenames, "\r?\n")
  filenames = unlist(filenames)
  filenames = Filter(function(x) !any(grepl("extent", x)), filenames)
  filenames = filenames[!filenames %in% c(".","..")]
  #print(paste(getwd(), "/", filename,sep = ""))
  for(filename in filenames) {
    download.file(paste(sea_url_single, filename, sep = ""), paste(getwd(), "/", filename,sep = ""), mode = "wb")
  }
}

files <- list.files(pattern='\\.tif$')
nfile1 <- length(files)
raster_1979_2021 <- lapply(files, raster) 
raster_1979_2021 <- stack(raster_1979_2021)


###extract values -------------------


abiot_sf <- abiot_total %>% 
  filter((slice %in% c('0-1'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station))

coordinates(abiot_sf) <- ~ long + lat
projection(abiot_sf) <- "+init=epsg:4326"
abiot_sf <- spTransform(abiot_sf, CRS(projection(raster_1979_2021$N_197811_concentration_v3.0)))

plot(raster_1979_2021$N_197811_concentration_v3.0)
points(abiot_sf)

rasValue_1979_2021 <- raster::extract(raster_1979_2021, abiot_sf)
abiot_sf <- abiot_total %>% 
  filter((slice %in% c('0-1'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station))
combinePointValue_1979_2021 = cbind(abiot_sf, rasValue_1979_2021)

write.table(combinePointValue_1979_2021,
            file="combinedPointValue_rasValue_1979_2021.csv", 
            append=FALSE, sep= ";", row.names = FALSE, col.names=TRUE)


### seasonal averages ------------
combinePointValue_1979_2021 <- read.table("C:/Users/trf/OneDrive - Universitetet i Oslo/PhD/maps/attributes/sea_ice/sea_ice_indexV3/monthly_concentration_1979_2020/combinedPointValue_rasValue_1979_2021.csv",
                                          sep = ";", header = T)

ice_1979_2021 <- combinePointValue_1979_2021
names(ice_1979_2021) = gsub(pattern = "_concentration_v3.0", replacement = "", x = names(ice_1979_2021))
names(ice_1979_2021) = gsub(pattern = "N_", replacement = "", x = names(ice_1979_2021))

ice_1979_2021 <- ice_1979_2021 %>% 
  rename_at(vars(ends_with("01")), ~ str_c(., "_winter")) %>%
  rename_at(vars(ends_with("02")), ~ str_c(., "_winter")) %>%
  rename_at(vars(ends_with("03")), ~ str_c(., "_winter")) %>%
  rename_at(vars(ends_with("04")), ~ str_c(., "_spring")) %>%
  rename_at(vars(ends_with("05")), ~ str_c(., "_spring")) %>%
  rename_at(vars(ends_with("06")), ~ str_c(., "_spring")) %>%
  rename_at(vars(ends_with("07")), ~ str_c(., "_summer")) %>%
  rename_at(vars(ends_with("08")), ~ str_c(., "_summer")) %>%
  rename_at(vars(ends_with("09")), ~ str_c(., "_summer")) %>%
  rename_at(vars(ends_with("10")), ~ str_c(., "_fall")) %>%
  rename_at(vars(ends_with("11")), ~ str_c(., "_fall")) %>%
  rename_at(vars(ends_with("12")), ~ str_c(., "_fall"))


#ice_1979_2021_winter <- ice_1979_2021[ , grepl( "winter" , names( ice_1979_2021 ) ) ]
#ice_1979_2021_winter <- mutate(ice_1979_2021_winter, mean_winter = rowMeans(ice_1979_2021_winter))


ice_1979_2021 <- ice_1979_2021 %>%
  mutate(mean_winter = rowMeans(ice_1979_2021[ , grepl( "winter" , names( . ) ) ]))%>%
  mutate(mean_spring = rowMeans(ice_1979_2021[ , grepl( "spring" , names( . ) ) ])) %>%
  mutate(mean_summer = rowMeans(ice_1979_2021[ , grepl( "summer" , names( . ) ) ])) %>%
  mutate(mean_fall = rowMeans(ice_1979_2021[ , grepl( "fall" , names( . ) ) ]))

# divide by 10 as guided in SEA ICE INDEX V3

ice_1979_2021 <- ice_1979_2021 %>%
  rowwise() %>%
  mutate(mean_winter = mean_winter/10)%>%
  mutate(mean_spring = mean_spring/10) %>%
  mutate(mean_summer = mean_summer/10) %>%
  mutate(mean_fall = mean_fall/10)


ice_1979_2021 <- ice_1979_2021 %>%
  select(1:9, mean_winter, mean_spring, mean_summer, mean_fall)  

# remove values lower than 15 as guided in SEA ICE INDEX V3
ice_1979_2021 <- ice_1979_2021 %>%
  rowwise() %>%
  mutate(mean_winter = ifelse(mean_winter < 15, 0, mean_winter))%>%
  mutate(mean_spring = ifelse(mean_spring < 15, 0, mean_spring)) %>%
  mutate(mean_summer = ifelse(mean_summer < 15, 0, mean_summer)) %>%
  mutate(mean_fall = ifelse(mean_fall < 15, 0, mean_fall))



#CTD DATA -------------
## CTD data import
#Husum, Katrine (2022). CTD data from Nansen Legacy Cruise - Paleo cruise. doi: 10.21335/NMDC-63841658.
#Gerland, Sebastian (2022). CTD data from Nansen Legacy Cruise - Seasonal Cruise Q1. doi: 10.21335/NMDC-1491279668.
#Ludvigsen, Martin (2022). CTD data from Nansen Legacy Cruise - Seasonal Cruise Q2. doi: 10.21335/NMDC-515075317.
#Reigstad, Marit (2022). CTD data from Nansen Legacy Cruise - Seasonal Cruise Q3. doi: 10.21335/NMDC-1107597377.
#Søreide, Janne (2022). CTD data from Nansen Legacy Cruise - Seasonal Cruise Q4. doi: 10.21335/NMDC-301551919.

setwd("C:/Users/trf/OneDrive - Universitetet i Oslo/PhD/dados/03_ctd_data")

npal4_ctd<-read.ctd(file="Paleo/Sta0165.cnv")
npal5_ctd<-read.ctd(file="Paleo/Sta0166.cnv")
npal7_ctd<-read.ctd(file="Paleo/Sta0168.cnv")
npal8_ctd<-read.ctd(file="Paleo/Sta0169.cnv")
npal12_ctd<-read.ctd(file="Paleo/Sta0171.cnv")
npal14_ctd<-read.ctd(file="Paleo/Sta0174.cnv")
npal15_ctd<-read.ctd(file="Paleo/Sta0176.cnv")
npal17_ctd<-read.ctd(file="Paleo/Sta0178.cnv")
npal19_ctd<-read.ctd(file="Paleo/Sta0179.cnv")

Q1P1_ctd<-read.ctd(file="Q1/Sta0119.cnv")
Q1P2_ctd<-read.ctd(file="Q1/Sta0122.cnv")
Q1P4_ctd<-read.ctd(file="Q1/Sta0135.cnv")
Q1P6_ctd<-read.ctd(file="Q1/Sta0151.cnv")
Q1P7_ctd<-read.ctd(file="Q1/Sta0157.cnv")

Q2P1_ctd<-read.ctd(file="Q2/Sta0163.cnv")
Q2P2_ctd<-read.ctd(file="Q2/Sta0172.cnv")
Q2P4_ctd<-read.ctd(file="Q2/Sta0182.cnv")
Q2P6_ctd<-read.ctd(file="Q2/Sta0199.cnv")
Q2P7_ctd<-read.ctd(file="Q2/Sta0210.cnv")

Q3P1_ctd<-read.ctd(file="Q3/Sta0147.cnv")
Q3P2_ctd<-read.ctd(file="Q3/Sta0158.cnv")
Q3P4_ctd<-read.ctd(file="Q3/Sta0165.cnv")
Q3P6_ctd<-read.ctd(file="Q3/Sta0182.cnv")
Q3P7_ctd<-read.ctd(file="Q3/Sta0191.cnv")
Q3SICE4_ctd<-read.ctd(file="Q3/Sta0192.cnv")

Q4P2_ctd<-read.ctd(file="Q4/Sta0450.cnv")
Q4P4_ctd<-read.ctd(file="Q4/Sta0444.cnv")
Q4P6_ctd<-read.ctd(file="Q4/Sta0431.cnv")

#### EISA
# I didn't had any CTD data directly from the EISA campaign. 
# However, this cruise took place in the same region and aroung the same period as SSQ3 (August 2019). 
# Thus, the CTD data was retrieve from one of the CTD casts closest to the coordinates and water depth location of the stations in the EISA cruise.

# St. 154 - close to station A
# St. 150 - close to station B

A_ctd<-read.ctd(file="EISA_Q3/Sta0154.cnv")
B_ctd<-read.ctd(file="EISA_Q3/Sta0150.cnv")


# CTD data as dataframe
npal4_ctd <- npal4_ctd@data %>%
  as_tibble()

npal5_ctd <- npal5_ctd@data %>%
  as_tibble()

npal7_ctd <- npal7_ctd@data %>%
  as_tibble()

npal8_ctd <- npal8_ctd@data %>%
  as_tibble()

npal12_ctd <- npal12_ctd@data %>%
  as_tibble()

npal14_ctd <- npal14_ctd@data %>%
  as_tibble()

npal15_ctd <- npal15_ctd@data %>%
  as_tibble()

npal17_ctd <- npal17_ctd@data %>%
  as_tibble()

npal19_ctd <- npal19_ctd@data %>%
  as_tibble()

Q1P1_ctd <- Q1P1_ctd@data %>%
  as_tibble()

Q1P2_ctd <- Q1P2_ctd@data %>%
  as_tibble()

Q1P4_ctd <- Q1P4_ctd@data %>%
  as_tibble()

Q1P6_ctd <- Q1P6_ctd@data %>%
  as_tibble()

Q1P7_ctd <- Q1P7_ctd@data %>%
  as_tibble()

Q2P1_ctd <- Q2P1_ctd@data %>%
  as_tibble()

Q2P2_ctd <- Q2P2_ctd@data %>%
  as_tibble()

Q2P4_ctd <- Q2P4_ctd@data %>%
  as_tibble()

Q2P6_ctd <- Q2P6_ctd@data %>%
  as_tibble()

Q2P7_ctd <- Q2P7_ctd@data %>%
  as_tibble()

Q3P1_ctd <- Q3P1_ctd@data %>%
  as_tibble()

Q3P2_ctd <- Q3P2_ctd@data %>%
  as_tibble()

Q3P4_ctd <- Q3P4_ctd@data %>%
  as_tibble()

Q3P6_ctd <- Q3P6_ctd@data %>%
  as_tibble()

Q3P7_ctd <- Q3P7_ctd@data %>%
  as_tibble()

Q3SICE4_ctd <- Q3SICE4_ctd@data %>%
  as_tibble()

Q4P2_ctd <- Q4P2_ctd@data %>%
  as_tibble()

Q4P4_ctd <- Q4P4_ctd@data %>%
  as_tibble()

Q4P6_ctd <- Q4P6_ctd@data %>%
  as_tibble()

A_ctd <- A_ctd@data %>%
  as_tibble()

B_ctd <- B_ctd@data %>%
  as_tibble()

# group CTD data

section <- bind_rows(NPAL4 = npal4_ctd,
                     NPAL5 = npal5_ctd,
                     NPAL7 = npal7_ctd,
                     NPAL8 = npal8_ctd,
                     NPAL12 = npal12_ctd,
                     NPAL14 = npal14_ctd,
                     NPAL15 = npal15_ctd,
                     NPAL17 = npal17_ctd,
                     NPAL19 = npal19_ctd,
                     Q1P1 = Q1P1_ctd,
                     Q2P1 = Q2P1_ctd,
                     Q3P1 = Q3P1_ctd,
                     A = A_ctd, B = B_ctd,
                     Q1P2 = Q1P2_ctd,
                     Q2P2 = Q2P2_ctd,
                     Q3P2 = Q3P2_ctd,
                     Q4P2 = Q4P2_ctd,
                     Q1P4 = Q1P4_ctd,
                     Q2P4 = Q2P4_ctd,
                     Q3P4 = Q3P4_ctd,
                     Q4P4 = Q4P4_ctd,
                     Q1P6 = Q1P6_ctd,
                     Q2P6 = Q2P6_ctd,
                     Q3P6 = Q3P6_ctd,
                     Q4P6 = Q4P6_ctd,
                     Q1P7 = Q1P7_ctd,
                     Q2P7 = Q2P7_ctd,
                     Q3P7 = Q3P7_ctd,
                     Q3SICE4 = Q3SICE4_ctd, .id = 'station')

# select and average surface and bottom TEMP, SAL, OXIGEN, FLUORESCENCE
# no fluorescence for NPAL- stations
ctd.t = section %>%
  group_by(station) %>%
  summarise(SST = mean(head(temperature, 10)),
            SBT = mean(tail(temperature, 10)),
            SSS = mean(head(salinity, 10)),
            SBS = mean(tail(salinity, 10)),
            SSO = mean(head(oxygen, 10)),
            SBO = mean(tail(oxygen, 10)),
            SSF = mean(head(fluorescence, 10)),
            SBF = mean(tail(fluorescence, 10))) %>%
  mutate(campaign = case_when(stri_detect_fixed(station,"NPAL") ~ "Paleo",
                              stri_detect_fixed(station,"A") ~ "EISA",
                              stri_detect_fixed(station,"B") ~ "EISA",
                              stri_detect_fixed(station,"Q1") ~ "SSQ1",
                              stri_detect_fixed(station,"Q2") ~ "SSQ2",
                              stri_detect_fixed(station,"Q3") ~ "SSQ3",
                              stri_detect_fixed(station,"Q4") ~ "SSQ4")) %>%
  mutate(station = stri_replace_all_fixed(station, c("Q1", "Q2", "Q3", "Q4"),  
                                          replacement = c(""), vectorize_all = FALSE)) %>%
  mutate(station =  factor(station, levels = order)) %>%
  arrange(station) %>%
  mutate(SSF = ifelse(SSF < 0, 0, SSF)) %>%
  mutate(SBF = ifelse(SBF < 0, 0, SBF))

# replicate data for 'replicates'
ctd.t$replicate <- "r1"
ctd2.t <- ctd.t %>% 
  slice(10:28) %>% # select seasonal cruises
  slice(-(29:30)) %>% # remove non existing replicates
  filter(!(campaign == "SSQ2" & station == "P1"))%>%
  filter(!(campaign == "SSQ2" & station == "P4"))%>%
  filter(!(campaign == "SSQ2" & station == "P6"))%>%
  filter(!(campaign == "SSQ2" & station == "P7"))
ctd2.t$replicate <- "r2"

ctd3.t <- ctd.t %>% 
  slice(10:28) %>% # select seasonal cruises
  filter(!(campaign == "SSQ1" & station == "P7"))%>%
  filter(!(campaign == "SSQ4" & station == "P2"))

ctd3.t$replicate <- "r3"

CTD <- bind_rows(ctd.t, ctd2.t, ctd3.t)
CTD$depth <- "0.5"


# Grain size ----------------------

grainsize <- abiot_total %>%
  filter((slice %in% c('0-1', '1-2', '2-3', '3-4', '4-5', '5-6', '2-4', '4-6'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  group_by(campaign, station, replicate) %>%
  summarise(clay = mean(clay, na.rm = TRUE), 
            silt = mean(silt, na.rm = TRUE), 
            sand = mean(sand, na.rm = TRUE), 
            mean_grain = mean(mean_grain, na.rm = TRUE)) %>%
  mutate(station =  factor(station, levels = order)) %>%
  arrange(station)

#is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

grainsize$mean_grain[is.nan(grainsize$mean_grain)]<-NA


#Land Distance -------
#NASA OBGP, O. B. P. G., NASA Goddard Space Flight Center, Richard P. Stumpf, and Norman A. Kuring (2009). Distance to Nearest Coastline: 0.01-Degree Grid. Available at: https://oceancolor.gsfc.nasa.gov/docs/distfromcoast/ [Accessed August 19, 2021].

raster_land_distance <- raster("GMT_intermediate_coast_distance_01d.tif")
plot(raster_land_distance)

##extract values -------------------

abiot_sf <- abiot_total %>% 
  filter((slice %in% c('0-1'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station))

coordinates(abiot_sf) <- ~ long + lat
#projection(abiot_sf) <- "+init=epsg:4326"
#abiot_sf <- spTransform(abiot_sf, CRS(projection(raster_land_distance)))

plot(raster_land_distance)
points(abiot_sf)


rasValue_land_distance = raster::extract(raster_land_distance, abiot_sf)
abiot_sf <- abiot_total %>% 
  filter((slice %in% c('0-1'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station))

combinePointValue_land_distance = cbind(abiot_sf, combinePointValue_land_distance$rasValue_land_distance)


# GMED DATA ---------------
# Available at: https://gmed.auckland.ac.nz/layersd.html
# Feldman, G. C., and McClain, C. (2010). ’Chlorophyll a concentration’ from bio-ORACLE dataset at http://oceancolor.gsfc.nasa.gov/. Compiled in: Basher, z., Bowden, d. A., Costello, m. J., 2019. Global marine environment datasets (GMED). Available at: http://gmed.auckland.ac.nz.

setwd("C:/Users/trf/OneDrive - Universitetet i Oslo/PhD/maps/attributes/pigments/gmed/chla")
 
#raster_chla_summer <- raster("kg_chla_sum_max.asc") too many NA's
#raster_chla_winter <- raster("kg_chla_win_max.asc") all NA's
raster_chla_av <- raster("bo_chla_mean.asc")
plot(raster_chla_av)

##extract values -------------------

abiot_sf <- abiot_total %>% 
  filter((slice %in% c('0-1'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station))

coordinates(abiot_sf) <- ~ long + lat
#projection(abiot_sf) <- "+init=epsg:4326"

plot(raster_chla_av)
points(abiot_sf)

chla_av = raster::extract(raster_chla_av, abiot_sf)

abiot_sf <- abiot_total %>% 
  filter((slice %in% c('0-1'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  select((1:9))

combineabiot_pp = cbind(abiot_sf, chla_av)

write.table(combineabiot_pp,
            file="gmed_pp_final.csv", dec = ",",
            append=FALSE, sep= ";", row.names = FALSE, col.names=TRUE)



## ALL together ------------

CTD # TEMPERATURE, SALINITY, OXYGEN
combineabiot_pp # Chla - PRODUCTIVITY
combinePointValue_land_distance # LAND DISTANCE
abiot_sf$prof # WATER DEPTH
ice_1979_2021 # SEA ICE SEASONAL 
grainsize # GRAIN SIZE

order  <- c("NPAL4", "NPAL5", "NPAL7",
            "NPAL8", "NPAL12", "NPAL14",
            "NPAL15", "NPAL17", "NPAL19",
            "P1",  "P2",  "P4",  "P6",  "P7","SICE4", "A", "B")

order_camp  <- c("Paleo", "SSQ1", "SSQ2",
                 "SSQ3", "SSQ4", "EISA")

CTD <- CTD %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  select(campaign, station, replicate, depth, SST:SBF) %>%
  mutate(campaign =  factor(campaign, levels = order_camp)) %>%
  mutate(station =  factor(station, levels = order)) %>%
  arrange(campaign) %>%
  arrange(station)

land_distance <- combinePointValue_land_distance %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  mutate(campaign =  factor(campaign, levels = order_camp)) %>%
  mutate(station =  factor(station, levels = order)) %>%
  arrange(campaign) %>%
  arrange(station)

grainsize <- grainsize %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  mutate(campaign =  factor(campaign, levels = order_camp)) %>%
  mutate(station =  factor(station, levels = order)) %>%
  arrange(campaign) %>%
  arrange(station)

combineabiot_pp <- combineabiot_pp %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  mutate(campaign = factor(campaign, levels = order_camp)) %>%
  mutate(station =  factor(station, levels = order)) %>%
  arrange(campaign) %>%
  arrange(station)

ice_1979_2021 <- ice_1979_2021 %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  mutate(campaign =  factor(campaign, levels = order_camp)) %>%
  mutate(station =  factor(station, levels = order)) %>%
  arrange(campaign) %>%
  arrange(station)

abiot_sf <- abiot_total %>% 
  filter((slice %in% c('0-1'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  mutate(campaign =  factor(campaign, levels = order_camp)) %>%
  mutate(station =  factor(station, levels = order)) %>%
  arrange(campaign) %>%
  arrange(station) %>%
  select(1:9)

env_total <- bind_cols(abiot_sf, 
                 'clay' = grainsize$clay,
                 'silt' = grainsize$silt,
                 'sand' = grainsize$sand,
                 'mean_grain' = grainsize$mean_grain,
                 'land_distance' = land_distance$`combinePointValue_land_distance$rasValue_land_distance`,
                 'SBT' = CTD$SBT,
                 'SST' = CTD$SST,
                 'SSS' = CTD$SSS,
                 'SBS' = CTD$SBS,
                 'SSO' = CTD$SSO,
                 'SBO' = CTD$SBO,
                 'chla_av' = combineabiot_pp$chla_av,
                 'ice_mean_winter' = ice_1979_2021$mean_winter,
                 'ice_mean_spring' = ice_1979_2021$mean_spring,
                 'ice_mean_summer' = ice_1979_2021$mean_summer,
                 'ice_mean_fall' = ice_1979_2021$mean_fall
)


