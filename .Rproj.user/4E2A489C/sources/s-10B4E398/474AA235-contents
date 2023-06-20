library(tidyverse)
library(sf)
library(stars)
library(tmap)
library(spdep)
library(gstat)
library(terra)

# read data --------------------------------------------------------------------
prec.sf <- st_read("data/precipitation_2014.shp") %>%
  st_transform(5514) %>%
  filter(!is.na(Jul))
summary(prec.sf$Jul)
regions <- st_read("data/kraje.shp")

# plot data --------------------------------------------------------------------
plot(prec.sf["Jul"])
tm_shape(regions %>% st_union) +
  tm_borders() +
  tm_shape(prec.sf) +
  tm_bubbles("Jul", title.size="Precipitation (mm)")

# elevation trend --------------------------------------------------------------
dtm <- read_stars("data/dtm_CR_1000.tif")
plot(dtm)
tm_shape(dtm) + 
  tm_raster(style="cont", palette = terrain.colors(12), title = "Elevation (m a.s.l.)") +
  tm_shape(regions %>% st_union) +
  tm_borders() +
  tm_shape(prec.sf) +
  tm_bubbles(size="Jul", alpha=0, border.col="black", title.size="Precipitation (mm)") +
  tm_layout(legend.position = c("right","top"))

# plotting relationship with elevation
prec.sf$elev <- st_extract(dtm, prec.sf) %>% pull(1)
ggplot(prec.sf, aes(x=elev, y=Jul)) +
  geom_point() +
  geom_smooth(se=T) +
  labs(x="Elevation (m a. s.)", y="Precipitation (mm)")

# linear model
mod.ols <- lm(Jul~elev, data=prec.sf)
summary(mod.ols)
elev.predictor <- dtm %>% rast()
names(elev.predictor) <- c("elev")
pred.ols <- predict(elev.predictor, mod.ols) %>%
  st_as_stars()
plot(pred.ols)

plot.ols <- tm_shape(pred.ols) +
  tm_raster(title="Precipitation (mm)") +
  tm_layout(legend.title.size = 1,
            title = "OLS regression on elevation",
            title.size = 1)
plot.ols

prec.sf$resid.ols <- mod.ols$residuals
prec.sf$resid.ols.abs <- abs(mod.ols$residuals)
prec.sf$resid.ols.sign <- mod.ols$residuals > 0
tm_shape(regions %>% st_union) +
  tm_borders() +
  tm_shape(prec.sf) +
  tm_bubbles(size = "resid.ols.abs", title.size = "Abs OLS residuals (mm)", 
             col = "resid.ols.sign", title.col = "Positive residuals") +
  tm_shape(prec.sf) +
  tm_dots() +
  tm_layout(legend.outside = T)

# autocorrelation structure ----------------------------------------------------

# using gstat
vario <- variogram(Jul~1, prec.sf)
vario.elev <- variogram(Jul~elev, prec.sf)
plot(vario)
plot(vario.elev)

(vgm.elev = fit.variogram(vario.elev, model = vgm("Sph")))
plot(vario.elev, vgm.elev)

(vgm.elev = fit.variogram(vario.elev, model = vgm("Exp")))
plot(vario.elev, vgm.elev)

# interpolation ----------------------------------------------------------------

# idw
pred.idw2.0 <- idw(Jul~1, prec.sf, newdata=dtm, idp=2)
pred.idw0.5 <- idw(Jul~1, prec.sf, newdata=dtm, idp=0.5)
pred.idw5.0 <- idw(Jul~1, prec.sf, newdata=dtm, idp=5)

plot.idw0.5 <- tm_shape(pred.idw0.5 %>% dplyr::select(var1.pred)) +
  tm_raster(title="Precip. (mm)") +
  tm_shape(prec.sf) +
  tm_dots() +
  tm_layout(legend.title.size = 1,
            legend.position = c("RIGHT","TOP"),
            title = "IDW: p = 0.5",
            title.size = 1)
plot.idw2.0 <- tm_shape(pred.idw2.0 %>% dplyr::select(var1.pred)) +
  tm_raster(title="Precip. (mm)") +
  tm_shape(prec.sf) +
  tm_dots() +
  tm_layout(legend.title.size = 1,
            legend.position = c("RIGHT","TOP"),
            title = "IDW: p = 2",
            title.size = 1)
plot.idw5.0 <- tm_shape(pred.idw5.0 %>% dplyr::select(var1.pred)) +
  tm_raster(title="Precip. (mm)") +
  tm_shape(prec.sf) +
  tm_dots() +
  tm_layout(legend.title.size = 1,
            legend.position = c("RIGHT","TOP"),
            title = "IDW: p = 5",
            title.size = 1)
tmap_arrange(plot.idw0.5,plot.idw2.0,plot.idw5.0)

# kriging with gstat -----------------------------------------------------------

# ordinary kriging
pred.krig.ord <- krige(Jul~1, prec.sf, newdata=dtm, model=model.const)
plot.krig.ord <- tm_shape(pred.krig.ord) +
  tm_raster(title=c("Precipitation (mm)", "Kriging variance")) +
  tm_shape(prec.sf) +
  tm_dots() +
  tm_layout(title = "Ordinary kriging",
            title.size = 1,
            legend.position = c("RIGHT","TOP"),)
plot.krig.ord

# universal kriging
names(dtm) <- "elev"
pred.krig.uni <- krige(Jul~elev, prec.sf, newdata=dtm, model=model.elev)
plot.krig.uni <- tm_shape(pred.krig.uni) +
  tm_raster(title=c("Precipitation (mm)", "Kriging variance")) +
  tm_shape(prec.sf) +
  tm_dots() +
  tm_layout(title = "Universal kriging",
            title.size = 1,
            legend.position = c("RIGHT","TOP"),)
tmap_arrange(plot.krig.ord, plot.krig.uni)
