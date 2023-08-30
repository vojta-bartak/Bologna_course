library(tidyverse) # General metapackage for manipuating and visualizing data
library(sf)        # Vector spatial data manipulation and analysis 
library(stars)     # Raster spatial data manipulation and analysis (compatible with sf)
library(tmap)      # Spatial data display
library(spdep)     # Areal data analysis, Moran's I, spatial dependence weighting
library(pgirmess)  # Function correlog
library(gstat)     # Geostatistics

# data
krnap <- read.table("data/krnap_richness.csv", sep=",", header=T) %>% 
  st_as_sf(coords=c("X", "Y"), crs=5514, remove=F)
krnap
dtm <- read_stars("data/dtm_CR_1000.tif")

# neighbors and weights
neibs <- dnearneigh(st_coordinates(krnap), d1=1000, d2=4000)
plot(neibs, coords=st_coordinates(krnap))

# Moran scatterplot
mp <- moran.plot(krnap$richness, listw = nb2listw(neibs, style="W"))

# local Moran's I
loc.mor <- localmoran(krnap$richness, nb2listw(neibs, style = "W"))
head(loc.mor)

loc.mor <- localmoran.exact(lm(richness~1, krnap), nb = neibs)
head(as.data.frame(loc.mor))

# residual autocorrelation
krnap$elev <- st_extract(dtm, krnap) %>% pull(1)

ggplot(krnap, aes(x=elev, y=richness)) +
  geom_point() +
  geom_smooth(method="lm", se=T)

krnap$residual[!is.na(krnap$elev)] <- lm(richness~elev, krnap)$residuals
tm_shape(dtm) +
  tm_raster(style="cont", palette="Blues") +
  tm_shape(krnap, is.master = T) +
  tm_bubbles("residual")

krnap.subset <- na.omit(krnap)
neibs.res <- dnearneigh(st_coordinates(krnap.subset), d1=1000, d2=2000)
weights.res <- nb2listw(neibs.res, style="B")

moran.plot(krnap.subset$residual, weights.res)

model <- lm(richness~elev, krnap.subset)
lm.morantest(model, weights.res)
moran.test(krnap.subset$residual, weights.res)

lm.morantest.exact(model, weights.res)

moran.mc(krnap.subset$residual, weights.res, nsim = 999)

correlog(st_coordinates(krnap.subset), krnap.subset$residual, 
         method = "Moran", nbclass = 30) %>% plot

variogram(richness~elev, data=krnap.subset, cutoff = 20000) %>% plot

mod.spaMM <- fitme(richness~elev+Matern(1|X+Y), data=krnap.subset, 
                   fixed=list(nu=1.5), family = poisson)
summary(mod.spaMM)
drop1(mod.spaMM)


d <- krnap %>% mutate(pos = numFactor(krnap$X, krnap$Y), 
                      group = factor(rep(1, nrow(krnap))))
mod.TMB <- glmmTMB(richness ~ elev + exp(pos + 0 | group), data=d,
                   family = poisson)
summary(mod.TMB)$coefficients


nb <- knn2nb(knearneigh(st_coordinates(prec.sf), k=12), sym=T)
lw <- nb2listwdist(nb, prec.sf, type="idw", alpha=1)
mod.sar <- spautolm(Jul~elev, data=prec.sf, listw = lw, family = "SAR") 
summary(mod.sar)

moran.mc(mod.sar$fit$residuals, lw, 999)

moran.mc(mod.ols$residuals, lw, nsim=999)
moran.mc(mod.gls$residuals, lw, nsim=999)
moran.mc(resid(mod.spaMM), lw, nsim=999)

mod.gam <- gam(Jul ~ elev + s(x, y), data=prec.sf %>% st_drop_geometry)
summary(mod.gam)
moran.mc(mod.gam$residuals, lw, 999)

plot(mod.gam)
vis.gam(mod.gam, view = c("x","y"))

mod.gam <- gam(Jul ~ elev + s(x, y, k=100), data=prec.sf %>% st_drop_geometry)
plot(mod.gam)
