 
library(AICcmodavg) 
library(brms) 
library(gamlss) 
library(gamm4) 
library(lme4) 
library(MuMIn) 
 
setwd("C:/Users/Julieta 2580/Google Drive/L 2020 Julieta Rojas/Datos") 
setwd("/Users/Edgar/Google\ Drive/Edgar/Trabajo/FC/Dirección\ de\ tesis/L\ 2020\ Julieta\ Rojas/Datos") 
setwd("/Users/edgarjgonzalezl/Google\ Drive/Edgar/Trabajo/FC/Dirección\ de\ tesis/L\ 2020\ Julieta\ Rojas/Datos") 
datos <- read.csv("Datos_P_bradyi_1988-2019_201116.csv") 
x.mean <- mean(c(datos$size, datos$sizeNext), na.rm = TRUE) 
x.sd <- sd(c(datos$size, datos$sizeNext), na.rm = TRUE) 
datos <- transform(datos,  
	size.s = (size-x.mean)/x.sd,  
	sizeNext.s = (sizeNext-x.mean)/x.sd 
) 

flower=datos$flower
datos2 <- (subset(x = datos, subset = flower > 0, ))
summary(datos2) 
# summary(datos) 
# plot(datos$size.s, datos$sizeNext, las = 1, bty = "l", xlab = "size.s at t", ylab = "sizeNextival",col=ifelse(datos$retracted == 0, "red", "blue")) 
 
# modelo nulo generalizado mixto 
mod.lmer.1 <- gamm4(flower ~ size.s, random = ~ (1|Cactus)+(1|year)+(1|Plot), family = negative.binomial, data = datos, REML = FALSE) 
save(mod.glmer.1, file = "mod.glmer.1.Rdata") 
 
# LINEALES 
mod.lmer.f1.11.11.11 <- glmer.nb(flower ~ size.s + (1+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), data = datos2)

mod.lmer.f1.01.11.11 <- glmer.nb(flower ~ size.s + (0+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), data = datos2)
 
save(mod.lmer.f1.11.11.11, file = "mod.lmer.f1.11.11.11.Rdata") 
 
# NO LINEALES 
mod.gamer.f1.11.11.11 <- gamm4(flower ~ s(size.s, k = 4), random = ~(1+size.s|Cactus)+(1+size.s|year)+(1+size.s|Plot), data = datos2, family = nb) 
save(mod.gamer.f1.11.11.11, file = "mod.gamer.g.11.11.11.Rdata") 
load("mod.gamer.f1.11.11.11.Rdata") 
AIC(mod.gamer.f1.11.11.11$mer)  
 
  
 
 
 
