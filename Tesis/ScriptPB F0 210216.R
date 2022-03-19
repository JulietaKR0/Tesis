unlink(".RData")

library(AICcmodavg)
library(brms)
library(gamlss)
library(gamm4)
library(lme4)
library(MuMIn)

setwd("C:/Users/Julieta 2580/Google Drive/L 2020 Julieta Rojas/Datos")
setwd("/Users/simulador/Downloads")
setwd("/Users/Edgar/Google\ Drive/Edgar/Trabajo/FC/Dirección\ de\ tesis/L\ 2020\ Julieta\ Rojas/Datos")
setwd("/Users/edgarjgonzalezl/Google\ Drive/Edgar/Trabajo/FC/Dirección\ de\ tesis/L\ 2020\ Julieta\ Rojas/Datos")
datos <- read.csv("Datos_P_bradyi_1988-2019_210122.csv")
x.mean <- mean(c(datos$size, datos$sizeNext), na.rm = TRUE)
x.sd <- sd(c(datos$size, datos$sizeNext), na.rm = TRUE)
datos <- transform(datos, 
	size.s = (size-x.mean)/x.sd, 
	sizeNext.s = (sizeNext-x.mean)/x.sd
)

# NULO
mod.glmer.f0.0 <- gamm4(repro ~ 1, random = ~ (1|Cactus)+(1|year)+(1|Plot), family = binomial, data = datos)
# AIC 4515.7

# LINEALES
mod.glmer.f0.11.11.11 <- glmer(repro ~ size.s + (1+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), family = binomial, data = datos) 
#NO CONVERGE

mod.glmer.f0.01.11.11 <- glmer(repro ~ size.s + (0+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.11.01.11 <- glmer(repro ~ size.s + (1+size.s|Cactus)+(0+size.s|year)+(1+size.s|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.11.11.01 <- glmer(repro ~ size.s + (1+size.s|Cactus)+(1+size.s|year)+(0+size.s|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.10.11.11 <- glmer(repro ~ size.s + (1|Cactus)+(1+size.s|year)+(1+size.s|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.11.10.11 <- glmer(repro ~ size.s + (1+size.s|Cactus)+(1|year)+(1+size.s|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.11.11.10 <- glmer(repro ~ size.s + (1+size.s|Cactus)+(1+size.s|year)+(1|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.00.11.11 <- glmer(repro ~ size.s + (1+size.s|year)+(1+size.s|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.11.00.11 <- glmer(repro ~ size.s + (1+size.s|Cactus)+(1+size.s|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.11.11.00 <- glmer(repro ~ size.s + (1+size.s|Cactus)+(1+size.s|year), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.01.01.11 <- glmer(repro ~ size.s + (0+size.s|Cactus)+(0+size.s|year)+(1+size.s|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.11.01.01 <- glmer(repro ~ size.s + (1+size.s|Cactus)+(0+size.s|year)+(0+size.s|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.01.11.01 <- glmer(repro ~ size.s + (0+size.s|Cactus)+(1+size.s|year)+(0+size.s|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.10.10.11 <- glmer(repro ~ size.s + (1|Cactus)+(1|year)+(1+size.s|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.11.10.10 <- glmer(repro ~ size.s + (1+size.s|Cactus)+(1|year)+(1|Plot), family = binomial, data = datos) 
# AIC 3047.7

mod.glmer.f0.10.11.10 <- glmer(repro ~ size.s + (1|Cactus)+(1+size.s|year)+(1|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.01.10.11 <- glmer(repro ~ size.s + (0+size.s|Cactus)+(1|year)+(1+size.s|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.01.11.10 <- glmer(repro ~ size.s + (0+size.s|Cactus)+(1+size.s|year)+(1|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.10.01.11 <- glmer(repro ~ size.s + (1|Cactus)+(0+size.s|year)+(1+size.s|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.11.01.10 <- glmer(repro ~ size.s + (1+size.s|Cactus)+(0+size.s|year)+(1|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.10.11.01 <- glmer(repro ~ size.s + (1|Cactus)+(1+size.s|year)+(0+size.s|Plot), family = binomial, data = datos) 
# NO CONVERGE

mod.glmer.f0.11.10.01 <- glmer(repro ~ size.s + (1+size.s|Cactus)+(1|year)+(0+size.s|Plot), family = binomial, data = datos) 
# NO CONVERGE

# BACKWARD
mod.glmer.f0.11.10.10 <- glmer(repro ~ size.s + (1+size.s|Cactus)+(1|year)+(1|Plot), family = binomial, data = datos) 
# AIC 3047.7

mod.glmer.f0.11.10.10 <- glmer(repro ~ size.s + (1+size.s|Cactus)+(1|year), family = binomial, data = datos) 
# NO CONVERGE 

mod.glmer.f0.10.10.10 <- glmer(repro ~ size.s + (1|Cactus)+(1|year)+(1|Plot), family = binomial, data = datos) 
# 3049.1

# NO LINEALES
mod.gamer.f0.11.11.11 <- gamm4(repro ~ s(size.s, k = 4), random = ~(1+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), family = binomial, data = datos)
summary(mod.gamer.f0.11.11.11$mer)
#AIC 2624.3

mod.gamer.f0.01.11.11 <- gamm4(repro ~ s(size.s, k = 4), random = ~(0+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), family = binomial, data = datos)
summary(mod.gamer.f0.01.11.11$mer)
#AIC 2618.7

mod.gamer.f0.00.11.11 <- gamm4(repro ~ s(size.s, k = 4), random = ~(1+size.s|year)+(1+size.s|Plot), family = binomial, data = datos)
summary(mod.gamer.f0.00.11.11$mer)
#AIC 2616.7

mod.gamer.f0.00.11.01 <- gamm4(repro ~ s(size.s, k = 4), random = ~(1+size.s|year)+(0+size.s|Plot), family = binomial, data = datos)
summary(mod.gamer.f0.00.11.01$mer)
#AIC 2652.3

mod.f0.best <- mod.gamer.f0.00.11.11
save(mod.f0.best, file = "mod.f0.best.Rdata")

# GRAFICA
x.pred<-seq(min(datos$size, na.rm = TRUE), max(datos$size, na.rm = TRUE), 0.1)
x.s.pred<-(x.pred-x.mean)/x.sd
y.s.pred <- predict(mod.f0.best, newdata = data.frame (size.s = x.s.pred), type = "response", re.form = NA) 
y.pred <- y.s.pred * x.sd + x.mean
plot(datos$size, datos$repro, las = 1, xlab = "Tamaño en el tiempo t (mm)", ylab = "Probabilidad de reproducción", cex = 0.5)
lines(x.pred, y.pred, col = "darkseagreen", lwd = 2)