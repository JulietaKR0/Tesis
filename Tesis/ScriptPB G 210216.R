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
datos <- read.csv("Datos_P_bradyi_1988-2019_211027.csv")
data <- subset(data, size < 80 & sizeNext < 80 & size > 0)
x.mean <- mean(c(datos$size, datos$sizeNext), na.rm = TRUE)
x.sd <- sd(c(datos$size, datos$sizeNext), na.rm = TRUE)
datos <- transform(datos, 
	size.s = (size-x.mean)/x.sd, 
	sizeNext.s = (sizeNext-x.mean)/x.sd
)

# modelo nulo generalizado mixto
mod.lmer.1 <- gamm4(sizeNext.s ~ 1, random = ~ (1|Cactus)+(1|year)+(1|Plot), family = gaussian, data = datos, REML = FALSE)
save(mod.lmer.0, file = "mod.lmer.0.Rdata")
# AIC 12090.1

# LINEALES
mod.lmer.g.11.11.11 <-lmer(sizeNext.s ~ size.s + (1+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), data = datos, REML = FALSE)
summary(mod.lmer.g.11.11.11)
# NO CONVERGE

mod.lmer.g.11.11.10 <-lmer(sizeNext.s ~ size.s + (1+size.s|Cactus)+(1+size.s|year)+(1|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.11.10.11 <-lmer(sizeNext.s ~ size.s + (1+size.s|Cactus)+(1|year)+(1+size.s|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.10.11.11 <-lmer(sizeNext.s ~ size.s + (1|Cactus)+(1+size.s|year)+(1+size.s|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.01.11.11 <-lmer(sizeNext.s ~ size.s + (0+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.11.01.11 <-lmer(sizeNext.s ~ size.s + (1+size.s|Cactus)+(0+size.s|year)+(1+size.s|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.11.11.01 <-lmer(sizeNext.s ~ size.s + (1+size.s|Cactus)+(1+size.s|year)+(0+size.s|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.00.11.11 <-lmer(sizeNext.s ~ size.s + (1+size.s|year)+(1+size.s|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.11.00.11 <-lmer(sizeNext.s ~ size.s + (1+size.s|Cactus)+(1+size.s|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.11.11.00 <-lmer(sizeNext.s ~ size.s + (1+size.s|Cactus)+(1+size.s|year), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.01.01.11 <-lmer(sizeNext.s ~ size.s + (0+size.s|Cactus)+(0+size.s|year)+(1+size.s|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.11.01.01 <-lmer(sizeNext.s ~ size.s + (1+size.s|Cactus)+(0+size.s|year)+(0+size.s|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.01.11.01 <-lmer(sizeNext.s ~ size.s + (0+size.s|Cactus)+(1+size.s|year)+(0+size.s|Plot), data = datos, REML = FALSE)
# AIC 9520.9

mod.lmer.g.10.10.11 <-lmer(sizeNext.s ~ size.s + (1|Cactus)+(1|year)+(1+size.s|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.11.10.10 <-lmer(sizeNext.s ~ size.s + (1+size.s|Cactus)+(1|year)+(1|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.10.11.10 <-lmer(sizeNext.s ~ size.s + (1|Cactus)+(1+size.s|year)+(1|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.01.10.11 <-lmer(sizeNext.s ~ size.s + (0+size.s|Cactus)+(1|year)+(1+size.s|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.01.11.10 <-lmer(sizeNext.s ~ size.s + (0+size.s|Cactus)+(1+size.s|year)+(1|Plot), data = datos, REML = FALSE)
# AIC 9419.4

mod.lmer.g.10.01.11 <-lmer(sizeNext.s ~ size.s + (1|Cactus)+(0+size.s|year)+(1+size.s|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.11.01.10 <-lmer(sizeNext.s ~ size.s + (1+size.s|Cactus)+(0+size.s|year)+(1|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.10.11.01 <-lmer(sizeNext.s ~ size.s + (1|Cactus)+(1+size.s|year)+(0+size.s|Plot), data = datos, REML = FALSE)
# NO CONVERGE

mod.lmer.g.11.10.01 <-lmer(sizeNext.s ~ size.s + (1+size.s|Cactus)+(1|year)+(0+size.s|Plot), data = datos, REML = FALSE)
# NO CONVERGE

# BACKWARD

mod.lmer.g.01.11.10 <-lmer(sizeNext.s ~ size.s + (0+size.s|Cactus)+(1+size.s|year)+(1|Plot), data = datos, REML = FALSE)
# AIC 9419.4

mod.lmer.g.00.11.10 <-lmer(sizeNext.s ~ size.s + (1+size.s|year)+(1|Plot), data = datos, REML = FALSE)
# AIC 9464

# NO LINEALES
mod.gamer.g.11.11.11 <- gamm4(sizeNext.s ~ s(size.s, k = 4), random = ~(1+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), data = datos, REML = FALSE)
# AIC 9337.3

mod.gamer.g.11.11.01 <- gamm4(sizeNext.s ~ s(size.s, k = 4), random = ~(1+size.s|Cactus)+(1+size.s|year)+(0+size.s|Plot), data = datos, REML = FALSE)
# AIC 9444

load("mod.gamer.g.11.11.11.Rdata")
mod.g.best <- mod.gamer.g.11.11.11
save(mod.g.best, file = "mod.g.best.Rdata")

# GRAFICA
x.pred<-seq(min(datos$size, na.rm = TRUE), max(datos$size, na.rm = TRUE), 0.1)
x.s.pred<-(x.pred-x.mean)/x.sd
y.s.pred <- predict(mod.g.best$gam, newdata = data.frame (size.s = x.s.pred), type = "response", re.form = NA) 
y.pred <- y.s.pred * x.sd + x.mean
plot(datos$size, datos$sizeNext, las = 1, xlab = "Tamaño en el tiempo t (mm)", ylab = "Tamaño en el tiempo t + 1 (mm)", cex = 0.5)
lines(x.pred, y.pred, col = "red", lwd = 2)
