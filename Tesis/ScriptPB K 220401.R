library(lattice)
library(plot3D)
library(popbio)
library(gamm4)

setwd("/Volumes/GoogleDrive/My\ Drive/Edgar/Trabajo/FC/Dirección\ de\ tesis/L\ 2020\ Julieta\ Rojas/Modelos")
load("mod.g0.best.Rdata") 
load("mod.dr.best.Rdata")
load("mod.r.best.Rdata")
load("mod.f1.best.Rdata")
load("mod.f0.best.Rdata")
load("mod.g.best.Rdata")
load("mod.s.best.Rdata")
load("mod.C.best.Rdata")

setwd("/Volumes/GoogleDrive/My\ Drive/Edgar/Trabajo/FC/Dirección\ de\ tesis/L\ 2020\ Julieta\ Rojas/Datos")
datos <- read.csv("Datos_P_bradyi_1988-2019_211027.csv")
x.mean <- mean(c(datos$size, datos$sizeNext), na.rm = TRUE)
x.sd <- sd(c(datos$size, datos$sizeNext), na.rm = TRUE)
datos <- transform(datos, 
	size.s = (size-x.mean)/x.sd, 
	sizeNext.s = (sizeNext-x.mean)/x.sd
)
x.s.min <- min(c(datos$size.s, datos$sizeNext.s), na.rm = TRUE)
x.s.max <- max(c(datos$size.s, datos$sizeNext.s), na.rm = TRUE)
n.i <- 100
z.s.i <- seq(x.s.min, x.s.max, length.out = n.i + 1)
x.s.i <- (z.s.i[1:n.i] + z.s.i[2:(n.i + 1)])/2
x.i <- x.s.i * x.sd + x.mean
s.i <- predict(mod.s.best$gam, newdata = data.frame(size.s = x.s.i), type = "response")
# plot(x.s.i, s.i, type = "l")
g.u.sigma <-  sum(as.data.frame(summary(mod.g.best$mer)$varcor)$sdcor[-c(3,6,9,10)])
g.u.mu <- predict(mod.g.best$gam, newdata = data.frame(size.s = x.s.i))
g.u.i.j <- matrix(NA, n.i, n.i)
for(i in 1:n.i){
	g.u.i.j[i, 1] <- pnorm(z.s.i[2], g.u.mu[i], g.u.sigma)
	g.u.i.j[i, n.i] <- 1 - pnorm(z.s.i[n.i], g.u.mu[i], g.u.sigma)
	for(j in 2:(n.i - 1))
		g.u.i.j[i, j] <- pnorm(z.s.i[j+1], g.u.mu[i], g.u.sigma) - pnorm(z.s.i[j], g.u.mu[i], g.u.sigma)
}
# hist3D(z = g.u.i.j)
f0.i <- predict(mod.f0.best$gam, newdata = data.frame(size.s = x.s.i), type = "response")
# plot(x.s.i, f0.i, type = "l")
f1.i <- predict(mod.f1.best, newdata = data.frame(size.s = x.s.i), type = "response", re.form = NA)
# plot(x.s.i, f1.i, type = "l")
r.i <- predict(mod.r.best, newdata = data.frame(size.s = x.s.i), type = "response", re.form = NA)
# plot(x.s.i, r.i, type = "l")
C.fn <- approxfun(mod.C.best)
C.i <- C.fn(x.i)
C.i[is.na(C.i)] <- 0
C.i <- C.i/sum(C.i)
# plot(x.i, C.i, type = "l", bty = "l" ,las = 1, xlab = expression (paste("Diámetro del cactus en el tiempo ", italic("t ") (mm))), ylab="Proporción de descendientes")
u.i <- predict(mod.dr.best, newdata = data.frame(size.s = x.s.i), type = "response", re.form = NA)
g.r.sigma <-  sum(as.data.frame(summary(mod.g0.best$mer)$varcor)$sdcor[-c(3,6,9,10)])
g.r.mu <- predict(mod.g0.best$gam, newdata = data.frame(size.s = x.s.i))
g.r.i.j <- matrix(NA, n.i, n.i)
for(i in 1:n.i){
	g.r.i.j[i, 1] <- pnorm(z.s.i[2], g.r.mu[i], g.r.sigma)
	g.r.i.j[i, n.i] <- 1 - pnorm(z.s.i[n.i], g.r.mu[i], g.r.sigma)
	for(j in 2:(n.i - 1))
		g.r.i.j[i, j] <- pnorm(z.s.i[j+1], g.r.mu[i], g.r.sigma) - pnorm(z.s.i[j], g.r.mu[i], g.r.sigma)
}
# hist3D(z = g.r.i.j)

# Kernels
# k.uu.i.j <- dinámica de la población no retraída que permanece no retraída
k.uu.i.j <- matrix(NA, n.i, n.i)
for(i in 1:n.i)
	for(j in 1:n.i)
		k.uu.i.j[i, j] <- g.u.i.j[i, j] * s.i[i] * (1-r.i[i]) + C.i[j] * f1.i[i] * f0.i[i] * (1-r.i[i])

# hist3D(z = k.uu.i.j, phi = 20, theta = 120, main = "Kernel de individuos no retraídos", xlab = "Diámetro en el tiempo t (mm)", ylab = "Diámetro en el tiempo t + 1 (mm)", zlab = "Kernel", cex.main = 1, cex.lab = 0.75)	
	
# k.i.j <- dinámica de la población si no hubiera retracción.
k.i.j <- matrix(NA, n.i, n.i)
for(i in 1:n.i)
	for(j in 1:n.i)
		k.i.j[i, j] <- g.u.i.j[i, j] * s.i[i] + C.i[j] * f1.i[i] * f0.i[i] 
# hist3D(z = k.i.j, phi = 20, theta = 120, main = "Kernel de la dinámica sin retracción", xlab = "Tamaño en el tiempo t (mm)", ylab = "Tamaño en el tiempo t + 1 (mm)", zlab = "Kernel", cex.main = 1, cex.lab = 0.75)						
lambda <- Re(eigen(k.i.j)$values[1])
# lambda asintótica sin retracción = 1.08
datosnt <-subset(datos, retracted == 0 & !is.na(size.s))
nt <- density(datosnt$size.s)
# plot(nu)
nt.fn <- approxfun(nt)
nt.i <- nt.fn(x.s.i)
nt.i[is.na(nt.i)] <- 0
nt.i <- nt.i/sum(nt.i)
yrs <- 3000
for(y in 1:yrs){
	nt1.i <- t(k.i.j) %*% nt.i
	nt.i <- nt1.i 
}
nt.iu <- nt.i/sum(nt.i)
plot(x.i, nt.iu, type = "l")

# k.ur.i.j <- dinámica de la población que estaba retraída y se desretrae.
k.ur.i.j <- matrix(NA, n.i, n.i)
for(i in 1:n.i)
	for(j in 1:n.i)
		k.ur.i.j[i, j] <- g.r.i.j[i, j] * s.i[i] * u.i[i]
# hist3D(z = k.ur.i.j, phi = 20, theta = 120, main = "Kernel de los individuos que se desretrajeron", xlab = "Diámetro en el tiempo t (mm)", ylab = "Diámetro en el tiempo t + 1 (mm)", zlab = "Kernel", cex.main = 1, cex.lab = 0.75)		

k.ru.i.j <- diag(r.i)
# hist3D(z = k.ru.i.j, phi = 20, theta = 120, main = "Kernel de los individuos que se desretrajeron", xlab = "Diámetro en el tiempo t (mm)", ylab = "Diámetro en el tiempo t + 1 (mm)", zlab = "Kernel", cex.main = 1, cex.lab = 0.75)		

k.rr.i.j <- diag(1-u.i)
# hist3D(z = k.rr.i.j, phi = 20, theta = 120, main = "Kernel de los individuos que se desretrajeron", xlab = "Diámetro en el tiempo t (mm)", ylab = "Diámetro en el tiempo t + 1 (mm)", zlab = "Kernel", cex.main = 1, cex.lab = 0.75)		

###nr###
datosnr<-subset(datos, retracted == 1 & !is.na(size))
nr <- density(datosnr$size)
# plot(nr)
nr.fn <- approxfun(nr)
nr.i <- nr.fn(x.i)
nr.i[is.na(nr.i)] <- 0
nr.i <- nr.i/sum(nr.i)
nr.i <- as.matrix(nr.i, 1, n.i)
###nu###
datosnu<-subset(datos, retracted == 0 & !is.na(size))
nu <- density(datosnu$size)
# plot(nu)
nu.fn <- approxfun(nu)
nu.i <- nu.fn(x.i)
nu.i[is.na(nu.i)] <- 0
nu.i <- nu.i/sum(nu.i)

###iteration###
yrs <- 100
nrt.i <- nr.i
nut.i <- nu.i
Nt <- sum(nr.i) + sum(nu.i)
lambdas <- rep(NA, yrs)
plot(x.i, nrt.i/sum(nrt.i), xlab = expression (paste("Diámetro del cactus en el tiempo ", italic("t ") (mm))), ylab="Proporción", type = "l", bty = "l", las = 1, col = "red", ylim = c(0,0.5))
lines(x.i, nut.i/sum(nut.i), col = "blue")
for(y in 1:yrs){
	nrt1.i <- t(k.ru.i.j) %*% nut.i + t(k.rr.i.j) %*% nrt.i
	nut1.i <- t(k.ur.i.j) %*% nrt.i + t(k.uu.i.j) %*% nut.i
 lines(x.i, nrt.i/sum(nrt1.i), col = "red")
 lines(x.i, nut.i/sum(nut1.i), col = "blue")
 Sys.sleep(1)
	Nt1 <- sum(nrt1.i) + sum(nut1.i)
	lambdas[y] <- Nt1/Nt 
	nrt.i <- nrt1.i 
	nut.i <- nut1.i
	Nt <- Nt1 
}
lambda <- lambdas[1000]
# lambda asintótica con retracción = 1.079

# Estandarizar tamaño poblacional, estructura estable de tamaños? PREGUNTAR
nt.i <- nrt1.i + nut1.i
nt.i <- nt.i/sum(nt.i)
plot(x.i, nt.i, xlab = expression (paste("Diámetro del cactus en el tiempo ", italic("t ") (mm))), ylab="Proporción", type = "l", bty = "l", las = 1)

# Valores reproductivos #
# Dinámica sin retracción
R.i <- Re(eigen(k.i.j)$vectors[,1])
R.i <- R.i/sum(R.i)
# plot(x.i, R.i, col = "red", type = "l", bty = "l", las = 1, xlab = expression (paste("Diámetro del cactus en el tiempo ", italic("t ") (mm))), ylab="Probabilidad de reproducción")
# mientras más grande eres, más valores reproductivos tienes #

# Dinámica con retracción
yrs <- 1000
nrt.i <- nr.i
nut.i <- nu.i
Nt <- sum(nr.i) + sum(nu.i)
lambdas <- rep(NA, yrs)
for(y in 1:yrs){
	nrt1.i <- k.ru.i.j %*% nut.i + k.rr.i.j %*% nrt.i
	nut1.i <- k.ur.i.j %*% nrt.i + k.uu.i.j %*% nut.i
	Nt1 <- sum(nrt1.i) + sum(nut1.i)
	lambdas[y] <- Nt1/Nt 
	nrt.i <- nrt1.i 
	nut.i <- nut1.i
	Nt <- Nt1 
}
rt.i <- nrt1.i + nut1.i
rt.i <- rt.i /sum(rt.i)
# lines(x.i, rt.i, col = "blue")
# 0.011 es la diferencia entre ambas gráficas.

# Sensibilidad #
mega.k <- matrix(NA, 2*n.i, 2*n.i) # 4 cuadros: SI: u->u, SD: u->r, II:r->u, ID:r->r
mega.k[1:n.i, 1:n.i] <- k.uu.i.j # SI: u->u
mega.k[1:n.i, n.i + 1:n.i] <- k.ur.i.j # SD: r->u
mega.k[n.i + 1:n.i, 1:n.i] <- k.ru.i.j # II: u->r
mega.k[n.i + 1:n.i, n.i + 1:n.i] <- k.rr.i.j # ID: r->r

hist3D(z = mega.k, phi = 20, theta = 150, xlab = "Diámetro en el tiempo t (mm)", ylab = "Diámetro en el tiempo t + 1 (mm)", zlab = "Kernel", cex.main = 1, cex.lab = 0.75)

# hist3D(z = k.uu.i.j, phi = 20, theta = 120, xlab = "Tamaño en el tiempo t (mm)", ylab = "Tamaño en el tiempo t + 1 (mm)", zlab = "Kernel", cex.main = 1, cex.lab = 0.5)

# hist3D(z = k.ur.i.j, phi = 20, theta = 120, xlab = "Tamaño en el tiempo t (mm)", ylab = "Tamaño en el tiempo t + 1 (mm)", zlab = "Kernel", cex.main = 1, cex.lab = 0.5)

# hist3D(z = k.ru.i.j, phi = 20, theta = 120, xlab = "Tamaño en el tiempo t (mm)", ylab = "Tamaño en el tiempo t + 1 (mm)", zlab = "Kernel", cex.main = 1, cex.lab = 0.75)

# hist3D(z = k.rr.i.j, phi = 20, theta = 120, xlab = "Tamaño en el tiempo t (mm)", ylab = "Tamaño en el tiempo t + 1 (mm)", zlab = "Kernel", cex.main = 1, cex.lab = 0.75)

# Heatmaps

a <- levelplot(k.uu.i.j, col.regions = rainbow(1000))
at <- a$panel.args.common$at
levelplot(k.uu.i.j, col.regions = rainbow(1000), at = at)

# Sacar cada una de las sensibilidades con base en:

get.Lam <- function(K) Re(eigen(K)$values[1])

S.k.uu <- sensitivity(k.uu.i.j)
S.k.ur <- sensitivity(k.ur.i.j)
S.k.ru <- sensitivity(k.ru.i.j)
S.k.rr <- sensitivity(k.rr.i.j)

R.i <- matrix(NA, n.i, n.i)
for(i in 1:n.i)
 R.i[, i] <- r.i
 
U.i <- matrix(NA, n.i, n.i)
for(i in 1:n.i)
 U.i[, i] <- u.i
 
S.s.i <- colSums(S.k.uu * g.u.i.j * (1 - R.i)) + colSums(S.k.ur * g.r.i.j * U.i)
# plot(x.i, S.s.i, type = "l")
S.s <- sum(S.s.i)
# 1.639 es la sensibilidad de lambda a la supervivencia

S.i <- matrix(NA, n.i, n.i)
for(i in 1:n.i)
 S.i[, i] <- s.i

gu.p.j <- rep(NA, n.i)
for(i in 1:n.i)
 gu.p.j[i] <- sum(g.u.i.j[i,]) - g.u.i.j[i,i]
 
S.gu.i.j <- S.k.uu * S.i * (1 - R.i)
for(i in 1:n.i)
 for(j in 1:n.i)
  S.gu.i.j[i,j] <- S.gu.i.j[i,j] + sum(- g.u.i.j[i,j] / sum(gu.p.j) * S.k.uu[i,j] * S.i[i,j] * (1 - R.i[i,j])) + g.u.i.j[i,i] / sum(gu.p.j) * S.k.uu[i,i] * S.i[i,i] * (1 - R.i[i,i])
#### hist3D(z = S.gu.i.j)
S.gu <- sum(S.gu.i.j)
# 169.32 es la sensibilidad de lambda al crecimiento de los individuos que no se retraen.

gr.p.j <- rep(NA, n.i)
for(i in 1:n.i)
 gr.p.j[i] <- sum(g.r.i.j[i,]) - g.r.i.j[i,i]
 
S.gr.i.j <- S.k.ur * S.i * U.i
for(i in 1:n.i)
 for(j in 1:n.i)
  S.gr.i.j[i,j] <- S.gr.i.j[i,j] + sum(- g.r.i.j[i,j] / sum(gr.p.j) * S.k.ur[i,j] * S.i[i,j] * U.i[i,j]) + g.r.i.j[i,i] / sum(gr.p.j) * S.k.ur[i,i] * S.i[i,i] * U.i[i,i]
#### hist3D(z = S.gr.i.j)
S.gr <- sum(S.gr.i.j)
# 61.59 es la sensibilidad de lambda al crecimiento de los individuos que se desretrajeron. 

C.j <- matrix(NA, n.i, n.i)
for(i in 1:n.i)
C.j[i,] <- C.i

F1.i <- matrix(NA, n.i, n.i)
for(i in 1:n.i)
 F1.i[, i] <- f1.i
 
F0.i <- matrix(NA, n.i, n.i)
for(i in 1:n.i)
 F0.i[, i] <- f0.i 

S.f0.i <- colSums(S.k.uu * C.j * F1.i *(1-R.i))
# plot(x.i, S.f0.i, type ="l")
S.f0 <- sum(S.f0.i)
# 0.113 es la sensibilidad de lambda a la probabilidad de reproducirse.

S.b.i <- colSums(S.k.uu * C.j * F0.i * (1-R.i))
# plot(x.i, S.b.i, type ="l")
S.b <- sum(S.b.i)
# 0.111 es la sensibilidad de lambda al número de estructuras reproductivas.

S.U.i <- colSums(S.k.ur * g.r.i.j * S.i) + colSums(S.k.ru) - colSums(S.k.rr)
# plot(x.i, S.U.i, type ="l")
S.U <- sum(S.U.i)
# 0.973 sensibilidad de lambda a la desretracción.

S.R.i <- colSums(S.k.uu * (-g.u.i.j * S.i - C.j * F1.i * F0.i))
# plot(x.i, S.R.i, type = "l")
S.R <- sum(S.R.i)
# -1.083 es la sensibilidad de lamdba a la retracción.

C.p.j <- rep(NA, n.i)
for(i in 1:n.i)
 C.p.j[i] <- sum(C.j[i,]) - C.j[i,i]

S.C.i <- S.k.uu * F1.i * F0.i * (1-R.i)
for(i in 1:n.i)
 for(j in 1:n.i)
  S.C.i[i,j] <- S.C.i[i,j] + sum(- C.j[i,j] / sum(C.p.j) * S.k.uu[i,j] * F1.i[i,j] * F0.i[i,j] * (1 - R.i[i,j])) + C.j[i,i] / sum(C.p.j) * S.k.uu[i,i] * F1.i[i,i] * F0.i[i,i] * (1 - R.i[i,i])
# plot(x.i, rowSums(S.C.i), type = "l")
S.C <- sum(S.C.i)
# 20.831 es la sensibilidad de lambda al tamaño de los descendientes.

# Elasticidades: dividir sensibilidad entre lambda/función

# Supervivencia
S.s.i.j<- S.k.uu * g.u.i.j * (1 - R.i) + S.k.ur * g.r.i.j * U.i
E.s.i.j<- S.s.i.j/(lambda/S.i)
E.s <- sum(E.s.i.j)
# 1.480

# Crecimiento de los individuos que no se retraen

E.gu.i.j <- S.gu.i.j/(lambda/g.u.i.j)
E.gu <- sum(E.gu.i.j)
# 0.906

# Crecimiento de los individuos que se desretraen

E.gr.i.j <- S.gr.i.j/(lambda/g.r.i.j)
E.gr <- sum(E.gr.i.j)
# 0.574

# Probabilidad de reproducirse
S.pb.i.j <- S.k.uu * C.j * F1.i *(1-R.i)
E.pb.i.j <- S.pb.i.j/(lambda/F0.i)
E.pb <- sum(E.pb.i.j)
# 0.070

# Número de esctructuras reproductivas
S.b.i.j <- S.k.uu * C.j * F0.i * (1-R.i)
E.b.i.j <- S.b.i.j/(lambda/F1.i)
E.b <- sum(E.b.i.j)
# 0.070

# Desretracción
S.u.i.j <- S.k.ur * g.r.i.j * S.i + S.k.ru - S.k.rr
E.u.i.j <- S.u.i.j/(lambda/U.i)
E.u <- sum(E.u.i.j)
# 0.631

# Retracción
S.r.i.j <- S.k.uu * (-g.u.i.j * S.i - C.j * F1.i * F0.i)
E.r.i.j <- S.r.i.j/(lambda/R.i)
E.r <- sum(E.r.i.j)
# -0.027

# Tamaño de los descendientes
S.c.i.j <- S.k.uu * F1.i * F0.i * (1-R.i)
E.c.i.j <- S.c.i.j/(lambda/C.j)
E.c <- sum(E.c.i.j)
# 0.070

E <- E.s + E.gu + E.gr + E.pb + E.b + E.u + E.r + E.c
E.s <- E.s/E
# 0.392
E.gu <- E.gu/E
# 0.240
E.gr <- E.gr/E
# 0.152
E.pb <- E.pb/E
# 0.019
E.b <- E.b/E
# 0.019
E.u <- E.u/E
# 0.167
E.r <- E.r/E
# -0.007
E.c <- E.c/E
# 0.019

# intervalos de confianza para lambdas
n.boot <- 1000
lambdas.boot <- matrix(NA, n.boot, 2)
for (b in 1:n.boot) {
 datos.boot <- datos[sample(1:nrow(datos), nrow(datos), TRUE),]
 
 mod.s.best.boot <- gamm4(surv ~ s(size.s, k = 4), random = ~(0 + size.s|Cactus) + (1 + size.s|year) + (1 + size.s|Plot), family = binomial, data = datos.boot)
 s.i <- predict(mod.s.best.boot$gam, newdata = data.frame(size.s = x.s.i), type = "response")
 
 datos.g.boot <- subset(datos.boot, size != 0 & sizeNext != 0)
 
 mod.g.best.boot <- gamm4(sizeNext.s ~ s(size.s, k = 4), random = ~(1 + size.s|Cactus) + (1 + size.s|year) + (1 + size.s|Plot), data = datos.g.boot)
 g.u.sigma <-  sum(as.data.frame(summary(mod.g.best.boot$mer)$varcor)$sdcor[-c(3,6,9,10)])
 g.u.mu <- predict(mod.g.best.boot$gam, newdata = data.frame(size.s = x.s.i))
 g.u.i.j <- matrix(NA, n.i, n.i)
 for(i in 1:n.i){
 	g.u.i.j[i, 1] <- pnorm(z.s.i[2], g.u.mu[i], g.u.sigma)
 	g.u.i.j[i, n.i] <- 1 - pnorm(z.s.i[n.i], g.u.mu[i], g.u.sigma)
 	for(j in 2:(n.i - 1))
 		g.u.i.j[i, j] <- pnorm(z.s.i[j+1], g.u.mu[i], g.u.sigma) - pnorm(z.s.i[j], g.u.mu[i], g.u.sigma)
 }
 
 mod.f0.best.boot <- gamm4(repro ~ s(size.s, k = 4), random = ~(1 + size.s|year) + (1 + size.s|Plot), family = binomial, data = datos.boot)
 f0.i <- predict(mod.f0.best.boot$gam, newdata = data.frame(size.s = x.s.i), type = "response")
 
 datos.f1.boot <- subset(datos.boot, flower > 0)
 datos.f1.boot$flower1 <- datos.f1.boot$flower - 1
 
 mod.f1.best.boot <- glmer(flower1 ~ size.s + (0+size.s|year)+(1|Plot), data = datos.f1.boot, family = poisson)
 f1.i <- predict(mod.f1.best.boot, newdata = data.frame(size.s = x.s.i), type = "response", re.form = NA) + 1
 
 mod.r.best.boot <- glmer(retractedNext ~ size.s + (1|Cactus) + (1|year) + (1 + size.s|Plot), family = binomial, data = datos.boot)
 r.i <- predict(mod.r.best.boot, newdata = data.frame(size.s = x.s.i), type = "response", re.form = NA)
 
 size.recruits <- c()
 for(i in 1:length(unique(datos.boot$Cactus))){
 datos.boot.ind <- subset(datos.boot, Cactus == unique (datos.boot$Cactus)[i])
 if (datos.boot.ind$year[1] != 1988)
 size.recruits <- c(size.recruits, datos.boot.ind$sizeNext[1])
 }
 size.recruits <- subset(size.recruits, is.na(size.recruits) == FALSE)
 km <- kmeans(log(size.recruits), 7) 
 size.recruits <- size.recruits[which(km$cluster == which(km$centers == min(km$centers)))]
 mod.C.best.boot <- density(size.recruits)
 C.fn <- approxfun(mod.C.best.boot)
 C.i <- C.fn(x.i)
 C.i[is.na(C.i)] <- 0
 C.i <- C.i/sum(C.i)
 
 datos.dr.boot <- subset(datos.boot, retracted == 1)
 mod.dr.best.boot <- lmer(retractedNext ~ size.s + (1 | year), data = datos.dr.boot)
 u.i <- predict(mod.dr.best.boot, newdata = data.frame(size.s = x.s.i), type = "response", re.form = NA)
 
 datos.g0.boot <- subset(datos.boot, size != 0 & sizeNext != 0)
 x.mean <- mean(c(datos.g0.boot$size, datos.g0.boot$sizeNext), na.rm = TRUE)
 x.sd <- sd(c(datos.g0.boot$size, datos.g0.boot$sizeNext), na.rm = TRUE)
 datos.g0.boot <- transform(datos.g0.boot, 
 	size.s = (size-x.mean)/x.sd, 
 	sizeNext.s = (sizeNext-x.mean)/x.sd
 )
 mod.g0.best.boot <- gamm4(sizeNext.s ~ s(size.s, k = 4), random = ~(1 + size.s|Cactus) + (1 + size.s|year) + (1 + size.s|Plot), data = datos.g0.boot)
 g.r.sigma <-  sum(as.data.frame(summary(mod.g0.best.boot$mer)$varcor)$sdcor[-c(3,6,9,10)])
 g.r.mu <- predict(mod.g0.best$gam, newdata = data.frame(size.s = x.s.i))
 g.r.i.j <- matrix(NA, n.i, n.i)
 for(i in 1:n.i){
 	g.r.i.j[i, 1] <- pnorm(z.s.i[2], g.r.mu[i], g.r.sigma)
 	g.r.i.j[i, n.i] <- 1 - pnorm(z.s.i[n.i], g.r.mu[i], g.r.sigma)
 	for(j in 2:(n.i - 1))
 		g.r.i.j[i, j] <- pnorm(z.s.i[j+1], g.r.mu[i], g.r.sigma) - pnorm(z.s.i[j], g.r.mu[i], g.r.sigma)
 }
 
 k.uu.i.j <- matrix(NA, n.i, n.i)
 for(i in 1:n.i)
 	for(j in 1:n.i)
 		k.uu.i.j[i, j] <- g.u.i.j[i, j] * s.i[i] * (1-r.i[i]) + C.i[j] * f1.i[i] * f0.i[i] * (1-r.i[i])
	
 # k.i.j <- dinámica de la población si no hubiera retracción.
 k.i.j <- matrix(NA, n.i, n.i)
 for(i in 1:n.i)
 	for(j in 1:n.i)
 		k.i.j[i, j] <- g.u.i.j[i, j] * s.i[i] + C.i[j] * f1.i[i] * f0.i[i] 
 lambdas.boot[b, 1] <- Re(eigen(k.i.j)$values[1])
	 
 # k.ur.i.j <- dinámica de la población que estaba retraída y se desretrae.
 k.ur.i.j <- matrix(NA, n.i, n.i)
 for(i in 1:n.i)
 	for(j in 1:n.i)
 		k.ur.i.j[i, j] <- g.r.i.j[i, j] * s.i[i] * u.i[i]

 k.ru.i.j <- diag(r.i)

 k.rr.i.j <- diag(1-u.i)

 ###nr###
 datosnr <- subset(datos, retracted == 0 & retractedNext == 1 & is.na(size) == FALSE & size > 0)
 nr <- density(datosnr$size)
 nr.fn <- approxfun(nr)
 nr.i <- nr.fn(x.i)
 nr.i[is.na(nr.i)] <- 0
 nr.i <- nr.i/sum(nr.i)
 nr.i <- as.matrix(nr.i, 1, n.i)
 
 ###nu###
 datosnu <- subset(datos, retracted == 0 & retractedNext == 0 & is.na(size) == FALSE & size > 0)
 nu <- density(datosnu$size)
 nu.fn <- approxfun(nu)
 nu.i <- nu.fn(x.i)
 nu.i[is.na(nu.i)] <- 0
 nu.i <- nu.i/sum(nu.i)

 ###iteration###
 yrs <- 1000
 nrt.i <- nr.i
 nut.i <- nu.i
 Nt <- sum(nr.i) + sum(nu.i)
 lambdas <- rep(NA, yrs)
 for(y in 1:yrs){
 	nrt1.i <- t(k.ru.i.j) %*% nut.i + t(k.rr.i.j) %*% nrt.i
 	nut1.i <- t(k.ur.i.j) %*% nrt.i + t(k.uu.i.j) %*% nut.i
 	Nt1 <- sum(nrt1.i) + sum(nut1.i)
 	lambdas[y] <- Nt1/Nt 
 	nrt.i <- nrt1.i 
 	nut.i <- nut1.i
 	Nt <- Nt1 
 }
 lambdas.boot[b, 2] <- lambdas[1000]
}
write.csv(lambdas.boot, file = "lambas.boot.csv")
# intervalo de confianza de lambda sin retracción
lambdas.wor.CI <- quantile(lambdas.boot[, 1], c(0.025, 0.975))
# intervalo de confianza de lambda con retracción
lambdas.wr.CI <- quantile(lambdas.boot[, 2], c(0.025, 0.975))
lambdas.CI <- cbind(lambdas.wor.CI, lambdas.wr.CI)
names(lambdas.CI) <- c("Without.retraction", "With.retraction")
write.csv(lambdas.CI, file = "lambdas.CI.csv")
plot(density(lambdas.boot[, 1]), col = "red")
lines(density(lambdas.boot[, 2]), col = "blue")
