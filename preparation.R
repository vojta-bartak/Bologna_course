library(tidyverse)
library(sf)
library(stars)
library(tmap)
library(spdep)
library(gstat)
library(geoR)

# data
krnap <- read.table("data/krnap_richness.csv", sep=",", header=T) %>% 
  st_as_sf(coords=c("X", "Y"), crs=5514, remove=F)
dtm <- read_stars("data/dtm_CR_1000.tif")

# visualization
tmap_mode("view")
tm_shape(krnap) +
  tm_dots("richness")

tmap_mode("plot")
tm_shape(dtm) +
  tm_raster(style="cont", palette="Blues") +
  tm_shape(krnap, is.master = T) +
  tm_dots("richness", size=1)

mp <- moran.plot(krnap$richness, nb2listw(sp.nb, style="B"))


# global tests
sp.nb <- dnearneigh(st_coordinates(krnap), d1=1000, d2=2000)
moran.test(krnap$richness, nb2listw(sp.nb))
moran.test(krnap$richness, nb2listw(sp.nb, style="B"))
moran.test(krnap$richness, nb2listw(sp.nb, style="B"), randomisation=F)
geary.test(krnap$richness, nb2listw(sp.nb))

lm.morantest(lm(richness~1, krnap), nb2listw(sp.nb, style="B"))
lm.morantest.sad(lm(richness~1, krnap), nb2listw(sp.nb, style="B"))
lm.morantest.exact(lm(richness~1, krnap), nb2listw(sp.nb, style="B"))

moran.mc(krnap$richness, nb2listw(sp.nb, style="B"), nsim=999)

dst <- nbdists(sp.nb, st_coordinates(krnap))
idw <- lapply(dst, function(x) 1/x)
moran.test(krnap$richness, nb2listw(sp.nb, glist = idw, style="B"))


# local tests
lm <- localmoran(krnap$richness, nb2listw(sp.nb, style="B"))
lm.sad <- localmoran.sad(model=lm(richness~1, krnap), nb=sp.nb)
lm.exact <- localmoran.exact(model=lm(richness~1, krnap), nb=sp.nb)
lm.exact2 <- localmoran.exact.alt(model=lm(richness~1, krnap), nb=sp.nb)
lm.perm <- localmoran_perm(krnap$richness, nb2listw(sp.nb, style="B"))

krnap$localI <- as.data.frame(lm)$Ii
tm_shape(krnap) +
  tm_dots("localI", size=2)

krnap$p.norm <- lm[,5]
krnap$p.perm <- lm.perm[,5]
krnap$p.sad <- lm.sad %>% sapply(function(x) x$p.value)
krnap$p.exact <- lm.exact %>% sapply(function(x) x$p.value)
krnap$p.exact2 <- lm.exact2 %>% sapply(function(x) x$p.value)
tm_shape(krnap) + tm_dots("p.norm", size=2)
tm_shape(krnap) + tm_dots("p.perm", size=2)
tm_shape(krnap) + tm_dots("p.sad", size=2)
tm_shape(krnap) + tm_dots("p.exact", size=2)
tm_shape(krnap) + tm_dots("p.exact2", size=2)

sum(lm[,5]<(0.05/348))
krnap$quadrant <- attr(lm, "quadr")$mean
krnap$influential <- mp$is_inf
tm_shape(krnap) +
  tm_dots("quadrant", size=2) +
  tm_shape(krnap %>% filter(influential)) +
  tm_dots(shape=1, size=2)


# correlogram
sp.correlogram(sp.nb, krnap$richness, order=10) %>% plot
sp.correlogram(sp.nb, krnap$richness, order=10, method = "C") %>% plot
sp.correlogram(sp.nb, krnap$richness, order=10, method = "I") %>% plot
sp.correlogram(sp.nb, krnap$richness, order=10, method = "I", style="B") %>% plot()
sp.correlogram(sp.nb, krnap$richness, order=10, method = "I", randomisation = F) %>% plot
sp.correlogram(sp.nb, krnap$richness, order=10, method = "I", randomisation = T) %>% print(p.adj.method = "BY")

library(pgirmess)
correlog(st_coordinates(krnap), krnap$richness, method = "Moran") %>% plot
correlog(st_coordinates(krnap), krnap$richness, method = "Moran", nbclass = 10) %>% plot


# variogram / covariogram
variogram(richness~1, data = krnap) %>% plot
variogram(richness~1, data = krnap, cutoff=30000) %>% plot
variogram(richness~1, data = krnap, cutoff=30000, covariogram=T) %>% plot
vario <- variogram(richness~1, data = krnap, cutoff=42000, covariogram=T)
plot(vario$dist, vario$gamma/max(vario$gamma))
data.frame(dist = vario$dist, cor = vario$gamma/max(vario$gamma)) %>% 
  arrange(desc(cor)) %>% 
  "["(2:nrow(.),) %>% 
  ggplot(aes(x=dist, y=cor)) +
  geom_line() +
  geom_point()

# residual autocorrelation
krnap$elev <- st_extract(dtm, krnap) %>% pull(1)
plot(richness~elev, krnap)
abline(lm(richness~elev, krnap))
ggplot(krnap, aes(x=elev, y=richness)) +
  geom_point() +
  geom_smooth(method="lm", se=T)

krnap$residual[!is.na(krnap$elev)] <- lm(richness~elev, krnap)$residuals
tm_shape(dtm) +
  tm_raster(style="cont", palette="Blues") +
  tm_shape(krnap, is.master = T) +
  tm_dots("residual", size=1)
tm_shape(dtm) +
  tm_raster(style="cont", palette="Blues") +
  tm_shape(krnap, is.master = T) +
  tm_bubbles("residual")

krnap.subset <- na.omit(krnap)
sp.nb.res <- dnearneigh(st_coordinates(krnap.subset), d1=1000, d2=2000)
moran.test(na.omit(krnap.subset$residual), nb2listw(sp.nb.res, style="B"))

lm.morantest(lm(richness~elev, krnap), weights)

moran.mc(krnap.subset$residual, nb2listw(sp.nb.res, style="B"), nsim=999)

lm.exact.res <- localmoran.exact(model=lm(residual~1, krnap.subset), nb=sp.nb.res)
krnap.subset$localI <- sapply(lm.exact.res, function(x) x$estimate)
krnap.subset$localp <- sapply(lm.exact.res, function(x) x$p.value)
tm_shape(krnap.subset) +
  tm_dots("localI", size=2)
tm_shape(krnap.subset) +
  tm_dots("localp", size=2)

correlog(st_coordinates(krnap.subset), krnap.subset$residual, method = "Moran") %>% plot

library(mgcv)
mgam <- gam(richness~elev+s(X, Y), data=krnap)
performance::check_overdispersion(mgam)

png("spline_surf.png", res=300, width = 15,height = 10, units = "cm")
plot(mgam)
dev.off()
