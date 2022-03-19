unlink(".RData")

library(AICcmodavg)
library(brms)
library(gamlss)
library(gamm4)
library(lme4)
library(MuMIn)

setwd("/Users/simulador/Google Drive/L 2020 Julieta Rojas/Datos")
setwd("/Users/simulador/Downloads")
setwd("C:/Users/Julieta 2580/Google Drive/L 2020 Julieta Rojas/Datos")
setwd("/Users/Edgar/Google\ Drive/Edgar/Trabajo/FC/Dirección\ de\ tesis/L\ 2020\ Julieta\ Rojas/Datos")
setwd("/Users/edgarjgonzalezl/Google\ Drive/Edgar/Trabajo/FC/Dirección\ de\ tesis/L\ 2020\ Julieta\ Rojas/Datos")
datos <- read.csv("Datos_P_bradyi_1988-2019_210122.csv")
x.mean <- mean(c(datos$size, datos$sizeNext), na.rm = TRUE)
x.sd <- sd(c(datos$size, datos$sizeNext), na.rm = TRUE)
datos <- transform(datos, 
	size.s = (size-x.mean)/x.sd, 
	sizeNext.s = (sizeNext-x.mean)/x.sd
)

# modelo nulo generalizado mixto
mod.glmer.r.0 <- gamm4(retracted ~ 1, random = ~ (1|Cactus)+(1|year)+(1|Plot), family = binomial, data = datos)
# AIC 4668.6

# LINEALES
mod.glmer.r.11.11.11 <- glmer(retracted ~ size.s + (1+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), family = binomial, data = datos)
summary(mod.glmer.r.11.11.11)
# NO CONVERGE

mod.glmer.r.01.11.11 <- glmer(retracted ~ size.s + (0+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), family = binomial, data = datos)
# NO CONVERGE
 
mod.glmer.r.11.01.11 <- glmer(retracted ~ size.s + (1+size.s|Cactus)+(0+size.s|year)+(1+size.s|Plot), family = binomial, data = datos)
# NO CONVERGE

mod.glmer.r.11.11.01 <- glmer(retracted ~ size.s + (1+size.s|Cactus)+(1+size.s|year)+(0+size.s|Plot), family = binomial, data = datos)
# AIC 3847.6

mod.glmer.r.10.11.11 <- glmer(retracted ~ size.s + (1|Cactus)+(1+size.s|year)+(1+size.s|Plot), family = binomial, data = datos)
# AIC 3816.8

mod.glmer.r.11.10.11 <- glmer(retracted ~ size.s + (1+size.s|Cactus)+(1|year)+(1+size.s|Plot), family = binomial, data = datos)
# AIC 3823.3

mod.glmer.r.11.11.10 <- glmer(retracted ~ size.s + (1+size.s|Cactus)+(1+size.s|year)+(1|Plot), family = binomial, data = datos)
# AIC 3830.9

# BACKWARD
mod.glmer.r.10.11.10 <- glmer(retracted ~ size.s + (1|Cactus)+(1+size.s|year)+(1|Plot), family = binomial, data = datos)
# AIC 3828.9

# NO LINEALES
mod.gamer.r.11.11.11 <- gamm4(retracted ~ s(size.s, k = 4), random = ~(1+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), family = binomial, data = datos)
summary(mod.gamer.r.11.11.11$mer)
# AIC 3807.4

mod.gamer.r.11.10.11 <- gamm4(retracted ~ s(size.s, k = 4), random = ~(1+size.s|Cactus)+(1|year)+(1+size.s|Plot), family = binomial, data = datos)
# AIC 3808.8

mod.r.best <- mod.gamer.r.11.11.11
save(mod.r.best, file = "mod.r.best.Rdata")

# GRAFICA
x.pred<-seq(min(datos$size, na.rm = TRUE), max(datos$size, na.rm = TRUE), 0.1)
x.s.pred<-(x.pred-x.mean)/x.sd
y.s.pred <- predict(mod.r.best$gam, newdata = data.frame (size.s = x.s.pred), type = "response", re.form = NA) 
y.pred <- y.s.pred * x.sd + x.mean
plot(datos$size, datos$retracted, las = 1, xlab = "Tamaño en el tiempo t (mm)", ylab = "Probabilidad de retracción", cex = 0.5)
lines(x.pred, y.pred, col = "darkorange", lwd = 2)