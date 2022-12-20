# Drivers of organic carbon distribution and accumulation in the northern Barents Sea
# Thaise Ricardo de Freitas
# 20 December 2022

# load packages -----------------------------

if(!require(pacman)) install.packages("pacman")
pacman::p_load(readxl, ggplot2, tidyverse, tidypaleo, patchwork, 
               grid, ggthemes, ggh4x, ggbreak, mgcv, khroma, stringi)


# load data -------------------------
# SEDIMENT GEOCHEMISTRY AND GRAIN SIZE
abiot_total <- read_xlsx("abiot_total.xlsx")

# RADIOMETRIC CONCENTRATIONS - NPAL CORES
concentrations <- read_xlsx("datacao_resumo.xlsx", 
                            sheet = 2)

# reorder stations - south to north

concentrations$core <- factor(concentrations$core, 
                              ordered = TRUE, 
                              levels = c("NPAL4","NPAL5", "NPAL7", "NPAL8","NPAL12",
                                         "NPAL14", "NPAL15", "NPAL17", "NPAL19"))



# Figure 2 ------------------
## RADIOMETRIC concentrations -----------

# top figure - southern shallower stations
fig2a <- alta_conc %>%
  filter(!(core %in% c("NPAL14", "NPAL15", "NPAL17", "NPAL19"))) %>%
  ggplot(aes(x = mean, y = depth, 
             shape = core)) +
  geom_errorbarh(aes(xmin = mean - sd, 
                     xmax = mean + sd))+
  geom_lineh() +
  geom_point(size = 2, na.rm = TRUE) +
  scale_color_colorblind()+
  scale_y_reverse(limits = c(20,0))+
  scale_shape_manual(values = c(10, 20, 2, 6,  15))+
  theme_paleo()+
  facet_geochem_gridh(vars(param), grouping = vars(core), 
                      scales = "free_x", labeller = label_parsed) +
  labs(x = NULL, y = "Depth (cm)")+
  theme(legend.title = element_blank(),
        legend.position = 'none'
  )

# bottom figure - northern deeper stations
fig2b <- alta_conc %>%
  filter((core %in% c("NPAL14", "NPAL15", "NPAL17", "NPAL19"))) %>%
  ggplot(aes(x = mean, y = depth, 
             shape = core)) +
  geom_errorbarh(aes(xmin = mean - sd, 
                     xmax = mean + sd))+
  geom_lineh() +
  geom_point(size = 2, na.rm = TRUE) +
  scale_color_colorblind()+
  scale_y_reverse(limits = c(20,0))+
  scale_shape_manual(values = c(18, 22, 3, 23))+
  theme_paleo()+
  facet_geochem_gridh(vars(param), grouping = vars(core), 
                      scales = "free", labeller = label_parsed) +
  #facetted_pos_scales(x = scales_x)+
  labs(x = NULL, y = "Depth (cm)")+
  theme(legend.title = element_blank(),
        legend.position = 'none',
        strip.text.x = element_blank()
  )

fig2a + fig2b + plot_layout(ncol = 1)
ggsave(file="fig2.pdf", width = 160, height = 270, units = "mm", dpi = "print")

# Figure 3 --------------

## Age model ------------------------
abiot_age <- abiot_total %>%
  filter(campaign %in% c('Paleo')) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPAL4_nomix', 'NPAL17'))) %>%
  mutate(station = factor(station, ordered = TRUE, levels = c("NPAL5", "NPAL7", 
                                                              "NPAL8","NPAL12",
                                                              "NPAL14", "NPAL15", 
                                                              "NPAL19"))) 

fig3a <- ggplot()+
  geom_ribbon(data = abiot_age, aes(xmin = age + year_se, xmax = age - year_se, 
                                    y = depth, group = station), 
              alpha = 0.2, na.rm = TRUE)+
  geom_line(data = abiot_age, aes(x = age, y = depth, linetype = station))+
  geom_point(data = abiot_age, aes(x = age, y = depth, shape = station), size = 2)+
  geom_linerange(data = abiot_age, aes(x = cs_year, ymin = depth - depth_cs_sd, 
                                       ymax = depth + depth_cs_sd), 
                 color = "red", na.rm = TRUE)+
  geom_point(data = abiot_age, aes(y = depth, x = cs_year, shape = station), 
             size = 2, color = 'red')+
  scale_shape_manual(values = c(20, 2, 6, 15, 18, 22, 23))+
  theme_paleo()+
  theme(
    legend.title = element_blank())+
  guides(fill=FALSE)+
  scale_y_reverse(name = "Depth (cm)", limits = c(20,0), 
                  expand = c(0,0), scales::pretty_breaks(n = 10))+
  scale_x_reverse(name = "Age (year AD)", position = "top", expand = c(0,0), 
                  scales::pretty_breaks(n = 5))

ggsave(file="fig2.png", width = 91, height = 150, units = "mm", dpi = "print")

## SAR -------------------------------------
# Sediment accumulation rate (SAR) - g cm-2 y-1
abiot_age <- abiot_total %>%
  filter(campaign %in% c('Paleo')) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPAL4_nomix', 'NPAL17'))) %>%
  mutate(station = factor(station, ordered = TRUE, levels = c("NPAL5", "NPAL7", 
                                                              "NPAL8","NPAL12",
                                                              "NPAL14", "NPAL15", 
                                                              "NPAL19"))) %>%
  drop_na(mar_gcm2)


fig3b <- ggplot(abiot_age)+
  geom_path(aes(y = depth, x = mar_gm2, linetype = station), color = "black", na.rm = TRUE) +
  geom_point(aes(y = depth, x = mar_gm2, shape = station), color = "black", size = 2, na.rm = TRUE) +
  scale_y_reverse(expand = c(0, 0), 
                  limits = c(20, 0), 
                  name = "", breaks = scales::pretty_breaks(n = 10)) +
  scale_x_continuous(position = 'top', name = expression(SAR~(g~cm^{"-2"}~y^-1)), 
                     breaks = scales::pretty_breaks(n = 10))+
  scale_shape_manual(values = c(20, 2, 6, 15, 18, 22, 23))+
  theme_paleo()+
  theme(legend.position="none", 
        legend.title = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(), plot.margin = margin(0, 0, 0, 0, "cm"))


ggplot(abiot_age)+
  geom_path(aes(y = age, x = mar_gm2, linetype = station), color = "black", na.rm = TRUE) +
  geom_point(aes(y = age, x = mar_gm2, shape = station), color = "black", size = 2, na.rm = TRUE) +
  scale_y_continuous(name = "", breaks = scales::pretty_breaks(n = 10)) +
  scale_x_continuous(position = 'top', name = expression(MAR~(g~cm^{"-2"}~y^-1)), 
                     breaks = scales::pretty_breaks(n = 10))+
  scale_shape_manual(values = c(20, 2, 6, 15, 18, 22, 23))+
  theme_paleo()

fig3a + fig3b


ggsave(file="fig3.pdf", width = 150, height = 180, units = "mm", dpi = "print")


### AGE MODELS for each core - TIDYPALEO -----------
# These are used in Figure 6
#NPAL5
alta <- abiot_trans %>%
  filter((station %in% c('NPAL5'))) 


model <- abiot_total %>%
  filter((station %in% c('NPAL5'))) %>%
  drop_na(age, depth)

age_npal5 = age_depth_model(model,
                            depth = depth, age = age,
                            age_max = age_max,
                            age_min = age_min)

#NPAL7
alta <- abiot_trans %>%
  filter((station %in% c('NPAL7'))) 


model <- abiot_total %>%
  filter((station %in% c('NPAL7'))) %>%
  drop_na(age, depth)

age_npal7 = age_depth_model(model,
                            depth = depth, age = age,
                            age_max = age_max,
                            age_min = age_min)

#NPAL8
alta <- abiot_trans %>%
  filter((station %in% c('NPAL8'))) 


model <- abiot_total %>%
  filter((station %in% c('NPAL8'))) %>%
  drop_na(age, depth)

age_npal8 = age_depth_model(model,
                            depth = depth, age = age,
                            age_max = age_max,
                            age_min = age_min)

#NPAL12
alta <- abiot_trans %>%
  filter((station %in% c('NPAL12'))) 


model <- abiot_total %>%
  filter((station %in% c('NPAL12'))) %>%
  drop_na(age, depth)

age_npal12 = age_depth_model(model,
                             depth = depth, age = age,
                             age_max = age_max,
                             age_min = age_min)

#NPAL14
alta <- abiot_trans %>%
  filter((station %in% c('NPAL14'))) 


model <- abiot_total %>%
  filter((station %in% c('NPAL14'))) %>%
  drop_na(age, depth)

age_npal14 = age_depth_model(model,
                             depth = depth, age = age,
                             age_max = age_max,
                             age_min = age_min)

#NPAL15
alta <- abiot_trans %>%
  filter((station %in% c('NPAL15'))) 


model <- abiot_total %>%
  filter((station %in% c('NPAL15'))) %>%
  drop_na(age, depth)

age_npal15 = age_depth_model(model,
                             depth = depth, age = age,
                             age_max = age_max,
                             age_min = age_min)

#NPAL19
alta <- abiot_trans %>%
  filter((station %in% c('NPAL19'))) 


model <- abiot_total %>%
  filter((station %in% c('NPAL19'))) %>%
  drop_na(age, depth)

age_npal19 = age_depth_model(model,
                             depth = depth, age = age,
                             age_max = age_max,
                             age_min = age_min)


# Figure 4 - GAM ------------------

## Latitude -------------------
# filter to remove 6-10cm sediment data, problematic replicate and reorder data
abiot_lat <- abiot_total %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  filter((slice %in% c('0-1', '1-2', '2-3', '3-4', '4-5', '5-6', '2-4', '4-6'))) %>%
  filter(!(station == 'P6' & replicate == 'r2' & campaign == 'SSQ1')) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  mutate(station = factor(station, ordered = TRUE, 
                          levels = c("A", "B","P1", "P2", "NPAL4", "NPAL5", "NPAL7", "NPAL8", "NPAL12", "P4",
                                     "NPAL14", "P6", "NPAL15", "P7", "NPAL19", "NPAL17", "SICE4")))


### TOC 63 -----------

fittoc_63 <- gam(toc_63 ~ s(lat), method = "REML", data = abiot_lat)
summary(fittoc_63)
plot(fittoc_63, all.terms = T, pages=1, scheme = 1, residuals = T, cex = 5)
anova.gam(fittoc_63)

fig3ab <-
  ggplot(abiot_lat, aes(y = toc_63, x = lat)) +
  geom_point(size = 1.5, alpha = .2)+
  stat_smooth(method = 'gam', formula = y ~ s(x), 
              color="black", fill = 'darkgray', size = .5)+
  scale_color_colorblind()+
  theme_paleo()+
  labs(y = expression(TOC[63]~('%')), 
       x = 'Latitude')+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )

### TIC -----------

fittic <- gam(tic ~ s(lat), data = abiot_lat)
summary(fittoc)
plot(fittic, all.terms = T, pages=1, scheme = 1, residuals = T, cex = 5)
anova.gam(fittic)

fig3c <-
  ggplot(abiot_lat, aes(y = tic, x = lat)) +
  geom_point(size = 1.5, alpha = .2)+
  stat_smooth(method = 'gam', formula = y ~ s(x), 
              color="black", fill = 'darkgray', size = .5)+
  scale_color_colorblind()+
  theme_paleo()+
  labs(y = 'TIC (%)', 
       x = 'Latitude')+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )

### TN -----------

fittn <- gam(tn ~ s(lat), data = abiot_lat)
summary(fittn)
plot(fittn, all.terms = T, pages=1, scheme = 1, residuals = T, cex = 5)
anova.gam(fittn)

fig3e <-
  ggplot(abiot_lat, aes(y = tn, x = lat)) +
  geom_point(size = 1.5, alpha = .2)+
  stat_smooth(method = 'gam', formula = y ~ s(x), 
              color="black", fill = 'darkgray', size = .5)+
  scale_color_colorblind()+
  theme_paleo()+
  labs(y = 'TN (%)', 
       x = 'Latitude')+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )

### CN -----------

fitcn <- gam(CN_org ~ s(lat), data = abiot_lat)
summary(fitcn)
plot(fitcn, all.terms = T, pages=1, scheme = 1, residuals = T, cex = 5)
anova.gam(fitcn)

fig3g <-
  ggplot(abiot_lat, aes(y = CN_org, x = lat)) +
  geom_point(size = 1.5, alpha = .2)+
  stat_smooth(method = 'gam', formula = y ~ s(x), 
              color="black", fill = 'darkgray', size = .5)+
  scale_color_colorblind()+
  theme_paleo()+
  labs(y = expression(C[org]/N[total]), 
       x = 'Latitude')+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )

## Water Depth -------------------


### TOC 63 -----------

fittoc_63 <- gam(toc_63 ~ s(prof), data = abiot_lat)
summary(fittoc_63)
plot(fittoc_63, all.terms = T, pages=1, scheme = 1, residuals = T, cex = 5)
anova.gam(fittoc_63)

fig3ba <-
  ggplot(abiot_lat, aes(y = toc_63, x = prof)) +
  geom_point(size = 1.5, alpha = .2)+
  stat_smooth(method = 'gam', formula = y ~ s(x), 
              color="black", fill = 'darkgray', size = .5)+
  scale_color_colorblind()+
  theme_paleo()+
  labs(y = expression(TOC[63]~('%')), 
       x = 'Water Depth (m)')+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )

### TIC -----------


fittic <- gam(tic ~ s(prof), data = abiot_lat)
summary(fittoc)
plot(fittic, all.terms = T, pages=1, scheme = 1, residuals = T, cex = 5)
anova.gam(fittic)

fig3d <-
  ggplot(abiot_lat, aes(y = tic, x = prof)) +
  geom_point(size = 1.5, alpha = .2)+
  stat_smooth(method = 'gam', formula = y ~ s(x), 
              color="black", fill = 'darkgray', size = .5)+
  scale_color_colorblind()+
  theme_paleo()+
  labs(y = 'TIC (%)', 
       x = 'Water Depth (m)')+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )

### TN -----------

fittn <- gam(tn ~ s(prof), data = abiot_lat)
summary(fittn)
plot(fittn, all.terms = T, pages=1, scheme = 1, residuals = T, cex = 5)
anova.gam(fittn)

fig3f <-
  ggplot(abiot_lat, aes(y = tn, x = prof)) +
  geom_point(size = 1.5, alpha = .2)+
  stat_smooth(method = 'gam', formula = y ~ s(x), 
              color="black", fill = 'darkgray', size = .5)+
  scale_color_colorblind()+
  theme_paleo()+
  labs(y = 'TN (%)', 
       x = 'Water Depth (m)')+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )

### CN -----------

fitcn <- gam(CN_org ~ s(prof), data = abiot_lat)
summary(fitcn)
plot(fitcn, all.terms = T, pages=1, scheme = 1, residuals = T, cex = 5)
anova.gam(fitcn)

fig3h <-
  ggplot(abiot_lat, aes(y = CN_org, x = prof)) +
  geom_point(size = 1.5, alpha = .2)+
  stat_smooth(method = 'gam', formula = y ~ s(x), 
              color="black", fill = 'darkgray', size = .5)+
  scale_color_colorblind()+
  theme_paleo()+
  labs(y = expression(C[org]/N[total]), 
       x = 'Water Depth (m)')+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )

fig3ab + fig3ba + fig3c + fig3d + fig3e + fig3f + fig3g + fig3h + plot_layout(ncol = 2)

ggsave(file="fig4.pdf", width = 180, height = 220, units = "mm", dpi = "print")

# Vertical plots Figure 5 & 6 ----------

## Figure 5 - Seasonal cruises --------------

alta_vertical <- abiot_total %>%
  filter((campaign %in% c('SSQ3', 'SSQ4', 'SSQ1', 'SSQ2'))) %>%
  filter(!(station %in% c('NPALP1'))) %>%
  filter(!(station == 'P6' & replicate == 'r2' & campaign == 'SSQ1')) %>%
  tidyr::pivot_longer(cols = c(wc, tic, toc_63, tn,CN_org,d13, d15), 
                      names_to = "param", values_to = "value") %>%
  group_by(campaign, station, slice, depth, param) %>%
  summarise(mean_value = mean(value),
            sd_value = sd(value),
            se_mean = sd(value)/sqrt(n())) %>%
  mutate(param = fct_relevel(param, "wc", "tic", "toc_63", "tn","CN_org", "d13", "d15"))

alta_vertical$param <- recode(alta_vertical$param,
                     `wc` = 'paste("WC (%)")',
                     `toc_63` = 'paste(TOC[63], " (%)")',
                     `tic` = 'paste("TIC (%)")',
                     `tn` = 'paste("TN (%)")',
                     `CN_org` = 'paste(C[org]/N[total])',
                     `d13` = 'paste(delta^13*C, "(\u2030)")',
                     `d15` = 'paste(delta^15*N, "(\u2030)")')



fig5 <- 
  ggplot(data = alta_vertical, aes(x = mean_value, y = depth, 
                          color = campaign, 
                          #shape = station,
                          linetype = campaign)) +
  geom_lineh() +
  geom_point(size = 2, na.rm = TRUE) +
  geom_errorbar(aes(xmin = mean_value - se_mean, 
                    xmax = mean_value + se_mean), 
                linetype = 'solid', width = 0)+
  scale_color_colorblind()+
  scale_y_reverse(breaks = seq(1, 5, by = 1)) +
  theme_paleo()+
  facet_geochem_gridh(vars(param), grouping = vars(station), 
                      scales = "free", labeller = label_parsed) +
  labs(x = NULL, y = "Depth (cm)")+
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        legend.justification = 'left',
  )

fig5

### Grain size -----------------

alta_grain <- abiot_total %>%
  filter((campaign %in% c('SSQ3'))) %>%
  filter(!(station %in% c('NPALP1'))) %>%
  filter((replicate %in% c('r1'))) %>%
  tidyr::pivot_longer(cols = c(clay, silt, sand), 
                      names_to = "param", values_to = "value") %>%
  mutate(param = fct_relevel(param, "clay", "silt", "sand")) %>%
  group_by(campaign, station, slice, depth, param)


alta_grain$param <- recode(alta_grain$param,
                           `clay` = "Clay (%)",
                           `silt` = "Silt (%)",
                           `sand` = "Sand (%)")

fig_grain_vol <- 
  ggplot(alta_grain) +
  geom_area(aes(fill = param, y = value, x = depth), position = "fill", 
            colour = "black", size = .2, alpha = .4)+
  scale_y_continuous(labels = scales::percent, expand = c(0,0))+
  scale_x_reverse()+
  theme_paleo()+
  scale_fill_brewer(palette = "Greys") +
  facet_wrap(~station, labeller = label_parsed, 
             ncol = 1, strip.position = 'right')+
  coord_flip()+
  labs(x = "", y = "")+
  theme(legend.title = element_blank(),
        legend.position = 'bottom'
  )

fig_grain_vol

# mean grain size - um
alta_grain_mean <- abiot_total %>%
  filter((campaign %in% c('SSQ3'))) %>%
  filter(!(station %in% c('NPALP1'))) %>%
  filter((replicate %in% c('r1'))) 


fig_grain_mean <- 
  ggplot(alta_grain_mean) +
  geom_point(aes(y = mean_grain, x = depth), size = 2, na.rm = TRUE)+
  scale_y_continuous()+
  scale_x_reverse()+
  theme_paleo()+
  facet_wrap(~station, labeller = label_parsed, 
             ncol = 1, strip.position = 'right')+
  coord_flip()+
  labs(x = "", y = "")+
  theme(legend.title = element_blank(),
        legend.position = 'top',
        legend.justification = 'left',
  )

fig_grain_mean

# combine plots volume grain size and mean grain size aftewards
fig5 + fig_grain + plot_layout(widths = c(5, .6))
ggsave(file="fig5_vol.pdf", width = 220, height = 270, units = "mm", dpi = "print")
fig5 + fig_grain_mean + plot_layout(widths = c(5, .6))
ggsave(file="fig5_meangrain.pdf", width = 220, height = 270, units = "mm", dpi = "print")


## Figure 6 - NPAL STATIONS --------------
alta_vertical <- abiot_total %>%
  filter((campaign %in% c('Paleo'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPAL4_nomix', 'NPAL17'))) %>%
  tidyr::pivot_longer(cols = c(car_gm2y, wc, tic, toc_63, tn, CN_org,d13, d15), 
                      names_to = "param", values_to = "value")  %>%
  mutate(station = factor(station, ordered = TRUE, levels = c("NPAL5", "NPAL7", 
                                                              "NPAL8","NPAL12",
                                                              "NPAL14", "NPAL15", 
                                                              "NPAL19"))) %>%
  mutate(param = fct_relevel(param, "wc", "car_gm2y", "tic", 
                             "toc_63", "tn","CN_org", "d13", "d15")) %>%
  drop_na(value)

alta_vertical$param <- recode(alta_vertical$param,
                     `wc` = 'paste("WC (%)")',
                     `toc_63` = 'paste(TOC[63], " (%)")',
                     `tic` = 'paste("TIC (%)")',
                     `tn` = 'paste("TN (%)")',
                     `CN_org` = 'paste(C[org]/N[total])',
                     `d13` = 'paste(delta^13*C, "(\u2030)")',
                     `d15` = 'paste(delta^15*N, "(\u2030)")')





scales_x <- list(
  param %in% c('paste("WC (%)")') ~ scale_x_continuous(limits = c(30, 80)),
  param %in% c('car_gm2y') ~ scale_x_continuous(limits = c(1, 22)),
  param %in% c('paste("TIC (%)")') ~ scale_x_continuous(limits = c(0, 1.5)),
  param %in% c('paste(TOC[63], " (%)")') ~ scale_x_continuous(limits = c(0, 2.5)),
  param %in% c('paste("TN (%)")') ~ scale_x_continuous(limits = c(0, 0.5)),
  param %in% c('paste(C[org]/N[total])') ~ scale_x_continuous(limits = c(4.5, 9.5)),
  param %in% c('paste(delta^13*C, "(\u2030)")') ~ scale_x_continuous(limits = c(-25.5, -22)),
  param %in% c('paste(delta^15*N, "(\u2030)")') ~ scale_x_continuous(limits = c(4, 7))
)


# add age models as second y axis - tidypaleo
scales_y <- list(
  scale_y_depth_age(age_npal5, age_name = "", limits = c(10,0),
                    age_breaks = seq(2019, 1953, -15)),
  scale_y_depth_age(age_npal7, age_name = "", limits = c(10,0),
                    age_breaks = seq(2019, 1932, -15)),
  scale_y_depth_age(age_npal8, age_name = "", limits = c(10,0),
                    age_breaks = seq(2019, 1960, -15)),
  scale_y_depth_age(age_npal12, age_name = "", limits = c(10,0),
                    age_breaks = seq(2019, 1852, -15)),
  scale_y_depth_age(age_npal14, age_name = "", limits = c(10,0),
                    age_breaks = seq(2019, 1858, -15)),
  scale_y_depth_age(age_npal15, age_name = "", limits = c(10,0),
                    age_breaks = seq(2019, 1821, -15)),
  scale_y_depth_age(age_npal19, age_name = "", limits = c(10,0),
                    age_breaks = seq(2019, 1810, -15))
)


fig6<- 
  ggplot(data = alta_vertical, aes(x = value, y = depth, linetype = station, shape = station)) +
  geom_lineh(na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  scale_shape_manual(values = c(20, 2, 6, 15, 18, 22, 23))+
  theme_paleo()+
  facet_geochem_gridh(vars(param), grouping = vars(station), 
                      scales = "free", labeller = label_parsed) +
  facetted_pos_scales(x = scales_x, y = scales_y)+
  labs(x = NULL, y = "Depth (cm)")+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom',
        #legend.justification = 'left'
        legend.position = 'none'
        
  )

fig6

### Grain size ----------------------

alta_grain_vol <- abiot_total %>%
  filter((campaign %in% c('Paleo'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPAL4_nomix', 'NPAL17'))) %>%
  tidyr::pivot_longer(cols = c(clay, silt, sand), 
                      names_to = "param", values_to = "value") %>%
  mutate(param = fct_relevel(param, "clay", "silt", "sand")) %>%
  group_by(campaign, station, depth, param)%>%
  mutate(station = factor(station, ordered = TRUE, levels = c("NPAL5", "NPAL7", 
                                                              "NPAL8","NPAL12",
                                                              "NPAL14", "NPAL15", 
                                                              "NPAL19")))


alta_grain_vol$param <- recode(alta_grain_vol$param,
                           `clay` = "Clay (%)",
                           `silt` = "Silt (%)",
                           `sand` = "Sand (%)")

fig_grain_vol <- 
  ggplot(alta_grain_vol) +
  geom_area(aes(fill = param, y = value, x = depth), position = "fill", 
            colour = "black", size = .2, alpha = .4)+
  scale_y_continuous(labels = scales::percent, expand = c(0,0))+
  scale_x_reverse(limits = c(10,0))+
  theme_paleo()+
  scale_fill_brewer(palette = "Greys") +
  facet_wrap(~station, labeller = label_parsed, 
             ncol = 1, strip.position = 'right')+
  coord_flip()+
  labs(x = "", y = "")+
  theme(legend.title = element_blank(),
        legend.position = 'none'
  )

fig_grain_vol


# mean grain size - um
alta_grain_mean <- abiot_total %>%
  filter((campaign %in% c('Paleo'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPAL4_nomix', 'NPAL17'))) %>%
  mutate(station = factor(station, ordered = TRUE, levels = c("NPAL5", "NPAL7", 
                                                              "NPAL8","NPAL12",
                                                              "NPAL14", "NPAL15", 
                                                              "NPAL19")))

fig_grain_mean <- 
  ggplot(alta_grain_mean) +
  geom_point(aes(y = mean_grain, 
                 x = depth, shape = station), size = 2, na.rm = TRUE)+
  scale_y_continuous()+
  scale_shape_manual(values = c(20, 2, 6, 15, 18, 22, 23))+
  scale_x_reverse(limits = c(10,0))+
  theme_paleo()+
  facet_wrap(~station, 
             ncol = 1, strip.position = 'right')+
  coord_flip()+
  labs(x = "", y = "")+
  theme(legend.title = element_blank(),
        legend.position = 'none',
        legend.justification = 'left',
  )

fig_grain_mean

# combine plots volume grain size and mean grain size aftewards
fig6 + fig_grain_vol + plot_layout(widths = c(5, .6))
ggsave(file="fig6_vol.pdf", width = 220, height = 270, units = "mm", dpi = "print")
fig6 + fig_grain_mean + plot_layout(widths = c(5, .6))
ggsave(file="fig6_meangrain.pdf", width = 220, height = 270, units = "mm", dpi = "print")



# Figure 7 ---------------------------

abiot_lat <- abiot_total %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPALP1'))) %>%
  filter((slice %in% c('0-1', '1-2', '2-3', '3-4', '4-5', '5-6', '2-4', '4-6'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  mutate(station = factor(station, ordered = TRUE, 
                          levels = c("A", "B", "P1", "P2", "NPAL4", "NPAL5", "NPAL7", "NPAL8", "NPAL12", "P4",
                                     "NPAL14", "NPAL15", "P6", "P7", "NPAL19", "NPAL17", "SICE4")))


### TOC_63 e TOC --------

fig7a <-
  ggplot(abiot_lat, aes(y = toc_63, x = toc, 
                        colour = station, 
                        shape = station)) +
  scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 2,
                                6, 15, 7, 18, 22, 8, 9,23, 3, 11))+
  scale_color_discreterainbow()+
  geom_point(size = 2)+
  theme_paleo()+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom',
        legend.position = 'none'
  )+
  labs(x = 'TOC (%)', 
       y = expression(TOC[63]~('%')))

### TOC_63 e d13 --------

fig7ab <-
  ggplot(abiot_lat, aes(y = toc_63, x = d13, 
                        colour = station, 
                        shape = station)) +
  scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 2,  
                                6, 15, 7, 18, 22, 8, 9,23, 3, 11))+
  scale_color_discreterainbow()+
  geom_point(size = 2)+
  theme_paleo()+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom',
        legend.position = 'none'
  )+
  labs(y = expression(TOC[63]~('%')), 
       x = expression(delta^13*C~('\u2030')))

### d15 e d13 --------

fig7b <-
  ggplot(abiot_lat, aes(y = d15, x = d13, 
                        colour = station,
                        shape = station)) +
  scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 2,  6, 15, 7, 18, 22, 8, 9,23, 3, 11))+
  scale_color_discreterainbow()+
  geom_point(size = 2)+
  theme_paleo()+
  theme(legend.title = element_blank(),
        legend.position = 'none'
  )+
  labs(y = expression(delta^15*N~('\u2030')), 
       x = expression(delta^13*C~('\u2030')))



### CN e d13  -----------------

fig7c <-
  ggplot(abiot_lat, aes(y = CN_org, x = d13, 
                        colour = station, 
                        shape = station)) +
  scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 
                                2, 6, 15, 7, 18, 22, 8, 
                                9,23, 3, 11))+
  scale_color_discreterainbow()+
  geom_point(size = 2)+
  theme_paleo()+
  theme(legend.title = element_blank(),
        legend.position = 'none'
  )+
  labs(y = expression(C[org]/N[total]), 
       x = expression(delta^13*C~('\u2030')))

ggplot(abiot_lat, aes(y = sand, x = lat, 
                      colour = station, 
                      shape = station)) +
  scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 
                                2, 6, 15, 7, 18, 22, 8, 
                                9,23, 3, 11))+
  scale_color_discreterainbow()+
  geom_point(size = 2)+
  theme_paleo()+
  theme(legend.title = element_blank(),
        legend.position = 'none'
  )+
  labs(y = "Silt (%)", 
       x = expression(delta^13*C~('\u2030')))

fig7a + fig7ab + fig7b + fig7c + plot_layout(ncol = 1)

ggsave(file="fig7.pdf", width = 90, height = 270, units = "mm", dpi = "print")



## Average 210Pb concentrations ---------

alta <- concentrations %>%
  tidyr::pivot_longer(cols = c(pb_total, pb_unsup, pb_sup, cs), 
                      names_to = "param", values_to = "value") %>%
  group_by(core, param) %>%
  summarise(mean_value = mean(value, na.rm = T),
            sd_value = sd(value,na.rm = T),
            se_mean = sd(value, na.rm = T)/sqrt(n()),
            min_value = min(value, na.rm = T),
            max_value = max(value, na.rm = T)) %>%
  mutate(param = fct_relevel(param, "pb_total", "pb_unsup", "pb_sup", "cs"))





# OLD CODE ------------------
## LATITUDE  --------

# without EISA values
abiot_lat <- abiot_total %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'A', 'B', 'NPALP1'))) %>%
  filter((slice %in% c('0-1', '1-2', '2-3', '3-4', '4-5', '5-6'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  mutate(station = factor(station, ordered = TRUE, 
                          levels = c("P1", "P2", "NPAL4", "NPAL5", "NPAL7", "NPAL8", "NPAL12", "P4",
                                     "NPAL14", "P6", "NPAL15", "P7", "NPAL19", "NPAL17", "SICE4")))
# with EISA values
abiot_lat <- abiot_total %>%
  filter(!(station %in% c('NPAL4_mixto6,5'))) %>%
  filter((slice %in% c('0-1', '1-2', '2-3', '3-4', '4-5', '5-6', '2-4', '4-6'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  mutate(station = factor(station, ordered = TRUE, 
                          levels = c("A", "B","NPALP1", "P1", "P2", "NPAL4", "NPAL5", "NPAL7", "NPAL8", "NPAL12", "P4",
                                     "NPAL14", "P6", "NPAL15", "P7", "NPAL19", "NPAL17", "SICE4")))

### TOC -----------

fig3a <- 
  ggplot(abiot_lat, aes(y = toc, x = lat, 
                        colour = group, 
                        shape = station)) +
  #scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+ # with EISA
  #scale_shape_manual(values = c(65, 66, 67, 4, 1, 10, 20, 2,  6, 15, 7, 18, 22, 8, 9,23, 3, 11))+ # with EISA plus NPAL7 and NPAL15
  #scale_shape_manual(values = c(4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+
  scale_shape_manual(values = c(4, 1, 10, 20, 2,  6, 15, 7, 18, 22, 8, 9,23, 3, 11))+ # to use when NPAL7 and 15 are done
  #scale_colour_viridis_d(option = "plasma")+
  scale_color_colorblind()+
  geom_point(size = 2)+
  theme_paleo()+
  labs(y = 'TOC (%)', 
       x = 'Latitude')+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )

### TIC -----------

fig3c <- 
  ggplot(abiot_lat, aes(y = tic, x = lat, 
                        colour = group, 
                        shape = station)) +
  #scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+
  #scale_shape_manual(values = c(4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+ # the normal one
  scale_shape_manual(values = c(4, 1, 10, 20, 2,  6, 15, 7, 18, 22, 8, 9,23, 3, 11))+ # to use when NPAL7 and 15 are done
  #scale_colour_viridis_d(option = "plasma")+
  scale_color_colorblind()+
  geom_point(size = 2)+
  theme_paleo()+
  labs(y = 'TIC (%)', 
       x = 'Latitude')+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )

### TN -----------

fig3e <- 
  ggplot(abiot_lat, aes(y = tn, x = lat, 
                        colour = group, 
                        shape = station)) +
  #scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+
  #scale_shape_manual(values = c(4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+ #the normal one
  scale_shape_manual(values = c(4, 1, 10, 20, 2,  6, 15, 7, 18, 22, 8, 9,23, 3, 11))+ # to use when NPAL7 and 15 are done
  #scale_colour_viridis_d(option = "plasma")+
  scale_color_colorblind()+
  geom_point(size = 2)+
  theme_paleo()+
  labs(y = 'TN (%)', 
       x = 'Latitude')+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )


### CN -----------

fig3g <- 
  ggplot(abiot_lat, aes(y = CN_org, x = lat, 
                        colour = replicate, 
                        shape = station)) +
  #scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+
  #scale_shape_manual(values = c(4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+ # the normal one
  scale_shape_manual(values = c(65, 66, 67, 4, 1, 10, 20, 2,  6, 15, 7, 18, 22, 8, 9,23, 3, 11))+ # to use when NPAL7 and 15 are done
  #scale_colour_viridis_d(option = "plasma")+
  scale_color_colorblind()+
  geom_point(size = 2)+
  theme_paleo()+
  labs(y = expression(C[org]/N[total]), 
       x = 'Latitude')+
  theme(legend.title = element_blank(),
        legend.position = 'bottom'
        #legend.position = 'none'
  )


## DEPTH  --------
### TOC -----------
fig3b <- 
  ggplot(abiot_lat, aes(y = toc, x = prof, 
                        colour = group, 
                        shape = station)) +
  #scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+
  #scale_shape_manual(values = c(4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+ # the normal one
  scale_shape_manual(values = c(4, 1, 10, 20, 2,  6, 15, 7, 18, 22, 8, 9,23, 3, 11))+ # to use when NPAL7 and 15 are done
  #scale_colour_viridis_d(option = "plasma")+
  scale_color_colorblind()+
  geom_point(size = 2)+
  theme_paleo()+
  labs(y = 'TOC (%)', 
       x = 'Water depth (m)')+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )

### TIC -----------

fig3d <- 
  ggplot(abiot_lat, aes(y = tic, x = prof, 
                        colour = group, 
                        shape = station)) +
  #scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+
  #scale_shape_manual(values = c(4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+ # the normal one
  scale_shape_manual(values = c(4, 1, 10, 20, 2,  6, 15, 7, 18, 22, 8, 9,23, 3, 11))+ # to use when NPAL7 and 15 are done
  #scale_colour_viridis_d(option = "plasma")+
  scale_color_colorblind()+
  geom_point(size = 2)+
  theme_paleo()+
  labs(y = 'TIC (%)', 
       x = 'Water depth (m)')+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )

### TN -----------

fig3f <- 
  ggplot(abiot_lat, aes(y = tn, x = prof, 
                        colour = group, 
                        shape = station)) +
  #scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+
  #scale_shape_manual(values = c(4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+ # the normal one
  scale_shape_manual(values = c(4, 1, 10, 20, 2,  6, 15, 7, 18, 22, 8, 9,23, 3, 11))+ # to use when NPAL7 and 15 are done
  #scale_colour_viridis_d(option = "plasma")+
  scale_color_colorblind()+
  geom_point(size = 2)+
  theme_paleo()+
  labs(y = 'TN (%)', 
       x = 'Water depth (m)')+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )


### CN -----------

fig3h <- 
  ggplot(abiot_lat, aes(y = CN_org, x = prof, 
                        colour = group, 
                        shape = station)) +
  #scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+
  #scale_shape_manual(values = c(4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+ # the normal one
  scale_shape_manual(values = c(4, 1, 10, 20, 2,  6, 15, 7, 18, 22, 8, 9,23, 3, 11))+ # to use when NPAL7 and 15 are done
  #scale_colour_viridis_d(option = "plasma")+
  scale_color_colorblind()+
  geom_point(size = 2)+
  theme_paleo()+
  labs(y = expression(C[org]/N[total]), 
       x = 'Water depth (m)')+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )

fig3h
ggsave(file="fig3_legenda.pdf", width = 180, height = 220, units = "mm", dpi = "print")

fig3a + fig3b + fig3c + fig3d + fig3e + fig3f + fig3g + fig3h + plot_layout(ncol = 2)

ggsave(file="fig3_Marie.pdf", width = 180, height = 220, units = "mm", dpi = "print")


## d15 e d13 --------

abiot_lat <- abiot_total %>%
  filter(!(campaign %in% c('EISA'))) %>%
  filter(!(station %in% c('NPAL4_mixto6,5', 'NPAL15', 'NPAL7'))) %>%
  #filter(!(station %in% c('NPAL4_mixto6,5'))) %>% 
  filter((slice %in% c('0-1', '1-2', '2-3', '3-4', '4-5', '5-6'))) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  #mutate(station = factor(station, ordered = TRUE, 
  #                        levels = c("P1","P2", "NPAL4", "NPAL5", "NPAL7", "NPAL8", "NPAL12", "P4","NPAL14","NPAL15", "P6", "P7", "NPAL19", "NPAL17", "SICE4")))
  mutate(station = factor(station, ordered = TRUE, levels = c("P1","P2", "NPAL4", "NPAL5", "NPAL8", "NPAL12", "P4",
                                                              "NPAL14", "P6", "P7", "NPAL19", "NPAL17", "SICE4")))

fig4a <-
  ggplot(abiot_lat, aes(y = d15, x = d13, 
                        colour = group, 
                        shape = station)) +
  #scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+
  scale_shape_manual(values = c(4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+ # the normal one
  #scale_shape_manual(values = c(4, 1, 10, 20, 2,  6, 15, 7, 18, 22, 8, 9,23, 3, 11))+ # to use when NPAL7 and 15 are done
  #scale_colour_viridis_d(option = "plasma")+
  scale_color_colorblind()+
  geom_point(size = 2)+
  theme_paleo()+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )+
  labs(y = expression(delta^15*N~('\u2030')), 
       x = expression(delta^13*C~('\u2030')))




## TOC e d13 --------

fig4b <-
  ggplot(abiot_lat, aes(y = toc, x = d13, 
                        colour = group, 
                        shape = station)) +
  #scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+
  scale_shape_manual(values = c(4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+ # the normal one
  #scale_shape_manual(values = c(4, 1, 10, 20, 2,  6, 15, 7, 18, 22, 8, 9,23, 3, 11))+ # to use when NPAL7 and 15 are done
  #scale_colour_viridis_d(option = "plasma")+
  scale_color_colorblind()+
  geom_point(size = 2)+
  theme_paleo()+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom',
        legend.position = 'none'
  )+
  labs(y = 'TOC (%)', 
       x = expression(delta^13*C~('\u2030')))

#fig4a / fig4b 
#ggsave(file="fig4.pdf", width = 90, height = 150, units = "mm", dpi = "print")


## CN e d13  -----------------

fig4c <-
  ggplot(abiot_lat, aes(y = CN_org, x = d13, 
                        colour = group, 
                        shape = station)) +
  #scale_shape_manual(values = c(65, 66, 4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+
  scale_shape_manual(values = c(4, 1, 10, 20, 6, 15, 7, 18, 8 ,9 ,23, 3, 11))+ # the normal one
  #scale_shape_manual(values = c(4, 1, 10, 20, 2,  6, 15, 7, 18, 22, 8, 9,23, 3, 11))+ # to use when NPAL7 and 15 are done
  #scale_colour_viridis_d(option = "plasma")+
  scale_color_colorblind()+
  geom_point(size = 2)+
  theme_paleo()+
  theme(legend.title = element_blank(),
        #legend.position = 'bottom'
        legend.position = 'none'
  )+
  labs(y = expression(C[org]/N[total]), 
       x = expression(delta^13*C~('\u2030')))

fig4b / fig4c
ggsave(file="fig5_cn.pdf", width = 90, height = 150, units = "mm", dpi = "print")

# Suplementary --------------

## Figure S1 -------------------------
## TOM average calculations
alta <- abiot_total %>%
  drop_na(TOM) %>%
  filter(!(campaign %in% c('EISA'))) %>%
  filter(!(station == 'NPAL4_mixto6,5')) %>%
  filter(!(station == 'NPALP1')) %>%
  filter(!(station == 'P6' & replicate == 'r2' & campaign == 'SSQ1')) %>%
  mutate(station = sub("NPAL4_nomix", "NPAL4", station)) %>%
  mutate(station = factor(station, 
                          ordered = TRUE, 
                          levels = c("P1","P2",
                                     "NPAL4",
                                     "NPAL5", "NPAL12", 
                                     "NPAL7", "NPAL8", "P4",
                                     "NPAL14", "P6", "NPAL15", 
                                     "P7", "NPAL19", "NPAL17", "SICE4"))) %>%
  
  tidyr::pivot_longer(cols = c(TOM), 
                      names_to = "param", values_to = "value") %>%
  group_by(campaign, station, depth, param) %>%
  summarise(mean_value = mean(value,na.rm = T),
            sd_value = sd(value,na.rm = T),
            se_mean = sd(value, na.rm = T)/sqrt(n()),
            min_value = min(value, na.rm = T),
            max_value = max(value, na.rm = T), .groups = 'drop') %>%
  mutate(param = fct_relevel(param, "TOM")) %>%
  mutate(station = factor(station, 
                          ordered = TRUE, 
                          levels = c("P1","P2",
                                     "NPAL4",
                                     "NPAL5", "NPAL12", 
                                     "NPAL7", "NPAL8", "P4",
                                     "NPAL14", "P6", "NPAL15", 
                                     "P7", "NPAL19", "NPAL17", "SICE4")))


fig_TOM <- 
  ggplot(data = alta, aes(x = mean_value, y = depth, 
                          color = campaign, 
                          linetype = campaign)) +
  geom_lineh() +
  geom_point(size = 2, na.rm = TRUE) +
  geom_errorbar(aes(xmin = mean_value - se_mean, 
                    xmax = mean_value + se_mean), 
                linetype = 'solid', width = 0)+
  scale_color_colorblind()+
  #scale_shape_manual(values = c(4, 1, 7, 8, 9, 11))+
  scale_y_reverse() +
  theme_paleo()+
  facet_geochem_gridh(vars(param), grouping = vars(station), 
                      scales = "free_x", labeller = label_parsed) +
  labs(x = NULL, y = "Depth (cm)")+
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        legend.justification = 'left',
        #strip.text.y = element_blank()
        #legend.position = 'none'
  )

fig_TOM

ggsave(file="figS1.pdf", width = 90, height = 270, units = "mm", dpi = "print")

# Figure S2 script available in rda_plot.R
