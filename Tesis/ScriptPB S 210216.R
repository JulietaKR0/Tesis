
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
x.mean <- mean(c(datos$size, datos$sizeNext), na.rm = TRUE)
x.sd <- sd(c(datos$size, datos$sizeNext), na.rm = TRUE)
datos <- transform(datos, 
	size.s = (size-x.mean)/x.sd, 
	sizeNext.s = (sizeNext-x.mean)/x.sd
)

# modelo nulo generalizado mixto
mod.glmer.1 <- gamm4(surv ~ 1, random = ~ (1|Cactus)+(1|year)+(1|Plot), family = binomial, data = datos)
# AIC 2903.3

# MODELOS LINEALES
mod.glmer.s.11.11.11 <- glmer(surv ~ size.s + (1+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), family = binomial, data = datos)
#AIC 2313.3

mod.glmer.s.01.11.11 <- glmer(surv ~ size.s + (0+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), family = binomial, data = datos)
# 2309.6

mod.glmer.s.01.10.11 <- glmer(surv ~ size.s + (0+size.s|Cactus)+(1|year)+(1+size.s|Plot), family = binomial, data = datos)
# AIC 2310.1

# NO LINEALES

mod.gamer.s.11.11.11 <- gamm4(surv ~ s(size.s, k = 4), random = ~(1+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), family = binomial, data = datos)
#AIC 2310.2

mod.s.best <- mod.glmer.s.01.11.11
save(mod.s.best, file = "mod.s.best.Rdata")

# GRAFICA
x.pred<-seq(min(datos$size, na.rm = TRUE), max(datos$size, na.rm = TRUE), 0.1)
x.s.pred<-(x.pred-x.mean)/x.sd
y.s.pred <- predict(mod.s.best, newdata = data.frame (size.s = x.s.pred), type = "response", re.form = NA) 
y.pred <- y.s.pred * x.sd + x.mean
plot(datos$size, datos$surv, las = 1, xlab = "Tamaño en el tiempo t (mm)", ylab = "Probabilidad de supervivencia", cex = 0.5)
lines(x.pred, y.pred, col = "blue", lwd = 2)
