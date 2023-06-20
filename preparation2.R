library(tidyverse)
library(sf)
library(stars)
library(tmap)
library(spdep)
library(gstat)
library(geoR)

# data
dtm <- read_stars("data/dtm_CR_1000.tif")
krnap <- read.table("data/krnap_richness.csv", sep=",", header=T) %>% 
  st_as_sf(coords=c("X", "Y"), crs=5514, remove=F) %>% 
  mutate(elev = st_extract(dtm, .) %>% pull(1)) %>% 
  na.omit

# variogram fitting
vario <- variogram(richness~1, data=krnap, cutoff = 30000)
plot(vario)

vario <- variogram(richness~1, data=krnap, cutoff = 30000)
plot(vario)

vario <- variogram(richness~elev, data=krnap, cutoff = 10000)
plot(vario)

(vario.fit <- fit.variogram(vario, model=vgm("Sph")))
plot(vario, vario.fit)
(vario.fit <- fit.variogram(vario, model=vgm("Exp")))
plot(vario, vario.fit)
(vario.fit <- fit.variogram(vario, model=vgm("Mat"), fit.kappa = TRUE))
plot(vario, vario.fit)




# kriging