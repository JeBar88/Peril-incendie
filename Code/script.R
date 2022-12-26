## Package
library(ggplot2)
library(actuar)
library(formatR)
library(ggplot2)
library(gridExtra)
library(knitr)
library(evir)
library(qrmtools)
library(ReIns)
library(tea)
library(goftest)
library(DescTools)
library(gamlss.dist)
library(GB2)
library(MASS)
library(randomcoloR)
dgpd <- ReIns::dgpd
pgpd <- ReIns::pgpd
qgpd <- ReIns::qgpd
rgpd <- ReIns::rgpd

#####################
### Fconction utilisé
#####################
Spl.ln.pg <- function(data, u, parInit, cont = F, deriv = F){
  
  if(cont == F & deriv == F){
    
    f <- function(par){
      par[5] * dlnorm(data, par[1], par[2])/plnorm(u, par[1], par[2]) * I(data <= u) +
        (1-par[5]) * dgpd(data, par[3], u, par[4]) * I(data > u)
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    ui <- rbind(diag(5), c(rep(0, 4), -1))
    ci <- c(rep(0, 5), -1)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=ui, ci=ci)
    mle$param <- mle$par
  }
  
  if(cont == T & deriv == F){
    
    w <- function(par) {
      (par[4]/plnorm(u, par[1], par[2]) * dlnorm(u, par[1], par[2]) + 1)^-1
    }
    
    f <- function(par){
      w(par) * dlnorm(data, par[1], par[2])/plnorm(u, par[1], par[2]) * I(data <= u) +
        (1-w(par)) * dgpd(data, par[3], u, par[4]) * I(data > u)
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=diag(4), ci=rep(0, 4))
    mle$param <- c(mle$par, w(mle$par))
  }
  
  if(cont == T & deriv == T){
    
    mu <- function(par){
      log(par[1]) - par[2]^2 * par[1] * (1 + par[3])/par[4]
    }
    
    w <- function(par) {
      (par[4]/plnorm(par[1], mu(par), par[2]) * dlnorm(par[1], mu(par), par[2]) + 1)^-1
    }
    
    f <- function(par){
      w(par) * dlnorm(data, mu(par), par[2])/plnorm(par[1], mu(par), par[2]) * I(data <= par[1]) +
        (1-w(par)) * dgpd(data, par[3], par[1], par[4]) * I(data > par[1])
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=diag(4), ci=rep(0, 4))
    mle$lim <- mle$par[1]
    mle$param <- c(mu(mle$par), mle$par[2],  mle$par[3:4], w(mle$par))
  }
  mle
}
Spl.we.pg <- function(data, u, parInit, cont = F, deriv = F){
  if(cont == F & deriv == F){
    
    f <- function(par){
      par[5] * dweibull(data, par[1], par[2])/pweibull(u, par[1], par[2]) * I(data <= u) +
        (1-par[5]) * dgpd(data, par[3], u, par[4]) * I(data > u)
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    ui <- rbind(diag(5), c(-1, rep(0, 4)) , c(rep(0, 4), -1))
    ci <- c(rep(0, 5), -5, -1)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=ui, ci=ci)
    mle$param <- mle$par
  }
  
  if(cont == T & deriv == F){
    
    w <- function(par) {
      (par[4]/pweibull(u, par[1], par[2]) * dweibull(u, par[1], par[2]) + 1)^-1
    }
    
    f <- function(par){
      w(par) * dweibull(data, par[1], par[2])/pweibull(u, par[1], par[2]) * I(data <= u) +
        (1-w(par)) * dgpd(data, par[3], u, par[4]) * I(data > u)
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    
    ui <- rbind(diag(4), c(-1, 0, 0, 0))
    ci <- c(rep(0, 4), -5)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=ui, ci=ci)
    mle$param <- c(mle$par, w(mle$par))
  }
  
  if(cont == T & deriv == T){
    
    b <- function(par){
      (1/(par[2] * par[1]^par[2]) * (par[1] * (par[3] + 1)/par[4] - (par[2] - 1)))^(-1/par[2])
    }
    
    w <- function(par) {
      (par[4] * dweibull(par[1], par[2], b(par))/pweibull(par[1], par[2], b(par)) + 1)^-1
    }
    
    f <- function(par){
      w(par) * dweibull(data, par[2], b(par))/pweibull(par[1], par[2], b(par)) * I(data <= par[1]) +
        (1-w(par)) * dgpd(data, par[3], par[1], par[4]) * I(data > par[1])
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=diag(4), ci=rep(0, 4))
    mle$lim <- mle$par[1]
    mle$param <- c(mle$par[2], b(mle$par),  mle$par[3:4], w(mle$par))
  }
  mle
}
spl.cox2.pg <- function(data, u, parInit, cont = F, deriv = F){
  
  if(cont==F & deriv==F){
    
    f <- function(par){
      par[6] * dcox2(data, par[1], par[2], par[3])/pcox2(u, par[1], par[2], par[3]) * I(data <= u) +
        (1-par[6]) * dgpd(data, par[4], u, par[5]) * I(data > u)
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    ui <- rbind(diag(6), c(-1, rep(0, 5)), c(rep(0, 5), -1))
    ci <- c(rep(0, 6), -1, -1)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=ui, ci=ci)
    mle$param <- mle$par
    mle
  }
  
  if(cont == T & deriv == F){
    
    w <- function(par) {
      (par[5]/pcox2(u, par[1], par[2], par[3]) * dcox2(u, par[1], par[2], par[3]) + 1)^-1
    }
    
    f <- function(par){
      w(par) *  dcox2(data, par[1], par[2], par[3])/pcox2(u, par[1], par[2], par[3]) * I(data <= u) +
        (1-w(par)) * dgpd(data, par[4], u, par[5]) * I(data > u)
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    ui <- rbind(diag(5), c(-1, rep(0, 4)))
    ci <- c(rep(0, 5), -1)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=ui, ci=ci)
    mle$param <- c(mle$par, w(mle$par))
  }
  mle
}
spl.gb2.pg <- function(data, u, parInit, cont = F){
  if(cont == F){
    
    f <- function(par) {
      par[7] * dgb2(data, par[1], par[2], par[3], par[4])/pgb2(u, par[1], par[2], par[3], par[4]) * I(data <= u) +
        (1 - par[7]) * dgpd(data, par[5], u, par[6]) * I(data > u)
    }
    
    logvrais <- function(par) {
      -sum(log(f(par)))
    }
    ui <- rbind(diag(7), c(rep(0, 6),-1))
    ci <- c(rep(0, 7),-1)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui = ui, ci = ci)
    mle$param <- mle$par
  }
  
  if(cont == T){
    
    w <- function(par) {
      (par[6]/pgb2(u, par[1], par[2], par[3], par[4]) * dgb2(u, par[1], par[2], par[3], par[4]) + 1)^-1
    }
    
    f <- function(par) {
      w(par) * dgb2(data, par[1], par[2], par[3], par[4])/pgb2(u, par[1], par[2], par[3], par[4]) * I(data <= u) +
        (1 - w(par)) * dgpd(data, par[5], u, par[6]) * I(data > u)
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    ui <- diag(6)
    ci <- numeric(6)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=ui, ci=ci)
    mle$param <- c(mle$par, w(mle$par))
  }
  mle
}

dln.pg <- function(x, u, par){ 
  m <- par[1] ; r <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5]
  
  p * dlnorm(x, m, r)/plnorm(u, m, r) * I(x <= u) +
    (1-p) * dgpd(x, xi, u, sigma) * I(x > u)
}
pln.pg <- function(x, u, par){ 
  m <- par[1] ; r <- par[2] ; xi <- par[3]
  sigma <- par[4] ; w <- par[5]
  
  w * plnorm(x, m, r)/plnorm(u, m, r) * I(x <= u) + 
    (w + (1-w) * pgpd(x, xi, u, sigma)) * I(x > u)
}
qln.pg <- function(k, u, par){
  m <- par[1] ; r <- par[2] ; xi <- par[3]
  sigma <- par[4] ; w <- par[5]
  ifelse(k <= w, qlnorm(pmin(k * plnorm(u, m, r)/w, 0.999), m, r), qgpd(pmax((k - w)/(1-w), 0), xi, u, sigma))
}
tln.pg <- function(k, u, par){
  m <- par[1] ; r <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5] ; c <- qln.pg(k, u, par)
  
  kk <- function(i) (1 + xi/sigma * (i - u))
  
  (p/(plnorm(u, m, r) * (1-k)) * exp(m + r^2/2) * (pnorm((log(u) - m - r^2)/r) - pnorm((log(c) - m - r^2)/r)) + 
      (1-p)/(1-k) * (u + sigma/(1 - xi))) * I(c <= u) + 
    ((1-p)/(1 - k) * (c * kk(c)^(-1/xi) - sigma/(xi - 1) * kk(c)^(1 - 1/xi))) * I(c > u)
}
rln.pg <- function(n, u, par){
  U <- runif(n)
  qln.pg(U, u, par)
}

dwe.pg <- function(x, u, par){
  t <- par[1] ; b <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5]
  
  p * dweibull(x, t, b)/pweibull(u, t, b) * I(x <= u) +
    (1-p) * dgpd(x, xi, u, sigma) * I(x > u)
}
pwe.pg <- function(x, u, par){
  t <- par[1] ; b <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5]
  
  p * pweibull(x, t, b)/pweibull(u, t, b) * I(x <= u) + 
    (p + (1-p) * pgpd(x, xi, u, sigma)) * I(x > u)
}
qwe.pg<- function(k, u, par){
  t <- par[1] ; b <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5]
  ifelse(k <= p, qweibull(pmin(k * pweibull(u, t, b)/p, 1), t, b), qgpd(pmax((k - p)/(1-p), 0), xi, u, sigma))
}
twe.pg <- function(k, u, par){
  t <- par[1] ; b <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5] ; c <- qwe.pg(k, u , par)
  
  kk <- function(i) (1 + xi/sigma * (i - u))
  
  (p/(pweibull(u, t, b) * (1-k)) * b * gamma(1 + 1/t) * (pgamma(u^t, 1+1/t, 1/b^t) - pgamma(c^t, 1+1/t, 1/b^t)) +
      (1-p)/(1-k) * (u + sigma/(1 - xi))) * I(c <= u) + 
    ((1-p)/(1 - k) * (c * kk(c)^(-1/xi) - sigma/(xi - 1) * kk(c)^(1 - 1/xi))) * I(c > u)
}
rwe.pg <- function(n, u, par){
  U <- runif(n)
  qwe.pg(U, u, par)
}

dcox2 <- function(x, p, b1, b2){
  p * dexp(x, b1) + (1-p) * (b2/(b2 - b1) * dexp(x, b1) + b1/(b1 - b2) * dexp(x, b2))
}
pcox2 <- function(x, p, b1, b2){
  p * pexp(x, b1) + (1-p) * (b2/(b2 - b1) * pexp(x, b1) + b1/(b1 - b2) * pexp(x, b2))
}
qcox2 <- function(k, p, b1, b2){
  
  f <- function(i) {optimize(function(x) abs(pcox2(x, p, b1, b2) - k[i]), c(0, 500000))$min}
  sapply(1:length(k), f)
}

pcox2.pg <- function(x, u, par){ 
  q <- par[1] ; b1 <- par[2] ; b2 <- par[3] ; 
  xi <- par[4] ; sigma <- par[5] ; p <- par[6]
  
  p * pcox2(x, q, b1, b2)/pcox2(u, q, b1, b2) * I(x <= u) + 
    (p + (1-p) * pgpd(x, xi, u, sigma)) * I(x > u)
}
qcox2.pg <- function(k, u, par){
  q <- par[1] ; b1 <- par[2] ; b2 <- par[3]
  xi <- par[4] ; sigma <- par[5] ; p <- par[6]
  ifelse(k <= p, qcox2(pmin(k * pcox2(u, q, b1, b2)/p, 0.999), q, b1, b2), qgpd(pmax((k - p)/(1-p), 0), xi, u, sigma))
}
rcox2.pg <- function(n, u, par){
  U <- runif(n)
  qcox2.pg(U, u, par)
}

pgb2.pg <- function(x, u, par){ 
  a <- par[1] ; b <- par[2] ; p <- par[3] ; q <- par[4]
  xi <- par[5] ; sigma <- par[6] ; w <- par[7]
  
  w * pGB2(x, b, a, p, q)/pGB2(u, b, a, p, q) * I(x <= u) + 
    (w + (1-w) * pgpd(x, xi, u, sigma)) * I(x > u)
}
qgb2.pg <- function(k, u, par){ 
  a <- par[1] ; b <- par[2] ; p <- par[3] ; q <- par[4]
  xi <- par[5] ; sigma <- par[6] ; w <- par[7]
  
  qGB2(pmin(k * pGB2(u, b, a, p, q)/w, 0.9999), b, a, p, q) * I(k <= w) + 
    qgpd(pmax((k - w)/(1-w), 0), xi, u, sigma) * I(k > w)
}
mgb2 <- function(k, a, b, p, q){
  
  b^k * beta(p + k/a, q - k/a)/beta(p, q)
  
}
rgb2.pg <- function(n, u, par){
  U <- runif(n)
  qgb2.pg(U, u, par)
}
tgb2.pg <- function(k, u, par){
  a <- par[1] ; b <- par[2] ; p <- par[3] ; q <- par[4] ;
  xi <- par[5] ; sigma <- par[6] ; w <- par[7] ; 
  c <- qgb2.pg(k, u, par) ; 
  
  kk <- function(i) (1 + xi/sigma * (i - u))
  h <- function(i) (i/b)^(a) / (1 + (i/b)^(a))
  
  (w/(pgb2(u, a, b, p, q) * (1-k)) * mgb2(1, a, b, p, q) * (pbeta(h(u), p + 1/a, q - 1/a) - pbeta(h(c), p + 1/a, q - 1/a)) + 
      (1-w)/(1-k) * (u + sigma/(1 - xi))) * I(c <= u) + 
    ((1-w)/(1 - k) * (c * kk(c)^(-1/xi) - sigma/(xi - 1) * kk(c)^(1 - 1/xi))) * I(c > u)
}

tr.ks <- function(data, u){
  
  param <- matrix(numeric(0), ncol = 2, nrow = length(u))
  stat <- numeric(length(u))
  p <- numeric(length(u))
  
  for(i in 1:length(u)){
    param[i, ] <- gpdFit(perte, u[i], method = "mle")$par.e
    test <- ks.test(data[data > u[i]], function(x) ReIns::pgpd(x, param[i, 2], u[i], param[i, 1]))
    stat[i] <- test[[1]]
    p[i] <- test[[2]]
  }
  list("lim - Stat min" = u[which.min(stat)])
}

VaR.PaG <- function(k, u){
  (s/xi * (((1 - k)/(1 - Fn(u)))^(-xi) - 1) + u) * I(k >= Fn(u))
} 
TVaR.PaG <- function(k, u){
  (s/xi * (((1 - k)/(1 - Fn(u)))^(-xi) * 1/(1-xi) - 1) + u) * I(k >= Fn(u))
} 

fond <- function(t, lam, prime, MI = 0, lim = Inf, info=F){
  freq <- rpois(t, lam)
  sev <- -sapply(freq, function(i) sum(pmin(rln.pg(i, u, param_ln.pg), lim)))
  balance <- c(MI, sev + prime)
  ruine <- ifelse(sum(cumsum(balance) < 0) == 0, 0, 1)
  if(info == T){
    list("Fréquence"=freq, "Sévérité"=sev, "Balance"=balance, "Nb.Perte"=ruine)
  }
  else {
    ruine
  }
}



####################
### Analyse initiale -----------------------------------------------------------
####################

### Base de données ###
don <- read.csv("Fire Incidents Data.csv", stringsAsFactors = T)
perte <- don$Estimated_Dollar_Loss[!is.na(don$Estimated_Dollar_Loss) & don$Estimated_Dollar_Loss > 0]

### Statistique descriptive ###
summary(perte)
Esp <- mean(perte)
Var <- var(perte)
n <- length(perte)
cbind(Esp, Var, n)

### Densité et log-densité ###
hist(perte, breaks = 10000, freq = T, main = "Histogramme des montants de sinistres estimés",
     xlab = "Montants de sinistres", ylab = "Nombre",
     xlim = c(0, 3e5))
hist(log(perte), breaks = 50, freq = T, main = "",
     xlab = "Log des montants de sinistres", ylab = "Nombre")

### Fonction de réparition ###
Fn <- ecdf(perte)
plot(Fn, main = "Fonction de répartition empirique de X", xlim = c(0, 5e6))

### QQplot exponentielle, lognormale et Pareto ###
# Exponentielle QQ-Plot
ExpQQ(perte)

# QQ plot Lognormale
LognormalQQ(perte)

# QQ plot Pareto
ParetoQQ(perte, xlim = c(0, 50))


### Fonction d'excès moyen ###
MeanExcess(perte, main = "Fonction d'excès moyen pour les montants de sinistres")

Hill(perte, plot = T, k=T)

###  Modèle lognormale sur toutes les données ###

# Méthode des moments
r1 <- sqrt(log(Var/Esp^2 + 1))
m1 <- log(Esp) - (r1^2)/2

# Maximume de vraisemblance
logvrais.LN <- function(par){
  -sum(log(dlnorm(perte, par[1], par[2])))
}

# Paramètres estimés
mle.LN <- constrOptim(c(m1, r1), logvrais.LN, grad = NULL, ui = diag(2), ci = c(0, 0))
m <- mle.LN$par[1]
r <- mle.LN$par[2]
cbind(m, r)

# Comparaison avec fonction de réparition empirique
plot(Fn, xlim = c(0, 1e7), ylim=c(0.9, 1), main = "")
curve(plnorm(x, m, r), col = "red", lwd = 2, add = T)

### Méthode POT ### 

# Minimiser la distance entre la statistique d'ordre la plus élevé et la pareto généralisé
#u <- mindist(sort(perte), method = "ks")$thres # 1E6

# Pas très constant avec b=10
#danielsson(sort(perte), B=1000)

# Constant mais trop bas
#gomes(perte, B=1000) # 150000

# Test la normalité pas constant
#u <- TH(sort(perte), seq(1E5, 1E7,,150))

# Test de seuil
u <- 75000
l <- 150000
# Trouver les paramètres de la portion Pareto généralisée
mle.PaG <- unname(gpdFit(perte, u, method = "mle")$par.e)
xi <- mle.PaG[2]
s <- mle.PaG[1]
cbind("alpha"=1/xi, "lambda"=s/xi, xi, s)

mle.PaG2 <- unname(gpdFit(perte, l, method = "mle")$par.e)
xi2 <- mle.PaG2[2]
s2 <- mle.PaG2[1]
cbind("alpha2"=1/xi2, "lambda2"=s2/xi2, xi2, s2)

# Comparaison avec fonction de réparition empirique
Fx.PaG <- function(x, u) ifelse(x >= u, (Fn(u) + (1-Fn(u)) * pgpd(x, xi, u, s)), NA)
Fx.PaG2 <- function(x, u) ifelse(x >= u, (Fn(u) + (1-Fn(u)) * pgpd(x, xi2, l, s2)), NA)

plot(Fn, xlim=c(0, 5e6), ylim=c(0.95, 1),
     lwd = 2, main = "")
curve(Fx.PaG(x, u), col = "red", lwd = 2, add = T)
curve(Fx.PaG2(x, l), col = "green", lwd = 2, add = T)
curve(plnorm(x, m, r), col = "blue", lwd = 2, add = T)
legend(3e6, 0.99, legend = c("LN", "PaG(u = 75000)", "PaG(u=150000)"), col = c("blue", "red", "green"), lwd = 2)

### VaR et TVaR de la portion Pareto généralisée ###
k <- c(0.95, 0.99, 0.995, 0.999)
VaR.PaG(k, u)
TVaR.PaG(k, u)


Info <- list("Stat_desc"=cbind(Esp, Var, n), "Param_LN"=cbind(m, r), "Param_PaG"=cbind(xi, s),
             "VaR.PaG_emp"=VaR.PaG(k, u), "TVaR.PaG_emp"=TVaR.PaG(k, u))


#############################
### Raccordement de deux lois --------------------------------------------------
#############################
########## Lognormale - Pareto généralisée ###########
  
parInit <- c(m, r, xi, s, Fn(u))

mle_ln.pg <- Spl.ln.pg(perte, u, parInit, cont = F, deriv = F)
param_ln.pg <- mle_ln.pg$par
cbind("m"=param_ln.pg[1], "r"=param_ln.pg[2], "a"= 1/param_ln.pg[3], "l"=param_ln.pg[4]/param_ln.pg[3], "p"=param_ln.pg[5])

# Graphique
plot(Fn, ylim=c(0, 1), xlim=c(0, 1e6), 
     main = "Comparaison de la fonction de répartition de la loi LN-PaG \navec la fonction de répartition empirique")
curve(pln.pg(x, u, param_ln.pg), add = T, col = "red", lwd = 2)

########## weibull - Pareto généralisée ########### 
  
# Méthode des moments pour les paramètres de la Weibull
f <- function(par) Esp^2 * (gamma(1 + 2/par)/gamma(1 + 1/par)^2 - 1)
t <- optimize(function(x) abs(f(x) - Var), c(0, 50))$min
b <- Esp/gamma(1 + 1/t)

# Optimisation splicing
parInit <- c(t, b, xi, s, Fn(u))
mle_we.pg <- Spl.we.pg(perte, u, parInit, cont = F)
mle_we.pg.d <- Spl.we.pg(perte, u, parInit[-5], cont = T)
param_we.pg <- mle_we.pg$param
param_we.pg.d <- mle_we.pg.d$param
cbind("t"=param_we.pg[1], "b"=param_we.pg[2], "a"=1/param_we.pg[3], "l"=param_we.pg[4]/param_we.pg[3], "p"=param_we.pg[5])
cbind("t"=param_we.pg.d[1], "b"=param_we.pg.d[2], "a"=1/param_we.pg.d[3], "l"=param_we.pg.d[4]/param_we.pg.d[3], "p"=param_we.pg.d[5])

# Graphique
plot(Fn, ylim=c(0, 1), xlim = c(0, 1e6), 
     main = "Comparaison de la fonction de répartition de la loi We-PaG \navec la fonction de répartition empirique")
curve(pwe.pg(x, u, param_we.pg), add = T, col = "red", lwd = 2)
curve(pwe.pg(x, u, param_we.pg.d), add = T, col = "red", lwd = 2)

########### coxienne-2 - Pareto généralisée ###########

# Optimisation splicing
MOM <- function(par){
  p <- 1 - 1/par[2] * (Esp - par[1])
  2*(p*par[1]^2 + (1 - p) * (par[1]^2 + par[2]^2 + (par[1] * par[2]))) - Esp^2
}
mle.MOM <- constrOptim(c(1e4, 1e6), function(x) abs(MOM(x) - Var), grad = NULL, ui=diag(2), ci=c(0, 0))

b1 <- 1/mle.MOM$par[1]
b2 <- 1/mle.MOM$par[2]
p <- 1 - b2 * (Esp - 1/b1)

lg <- function(par){
  -sum(log(dcox2(perte, par[1], par[2], par[3])))
}
ui <- rbind(diag(3), c(-1, rep(0, 2)))
ci <- c(rep(0, 3), -1)
mle.cox2 <- constrOptim(c(p, b1, b2), lg, grad = NULL, ui=ui, ci=ci)

parInit <- c(mle.cox2$par[1], mle.cox2$par[2], mle.cox2$par[3], xi, s, Fn(u))
mle_cox2.pg <- spl.cox2.pg(perte, u, parInit, cont = F)
mle_cox2.pg.d <- spl.cox2.pg(perte, u, parInit[-6], cont = T)
param_cox2.pg <- mle_cox2.pg$param
param_cox2.pg.d <- mle_cox2.pg.d$param
cbind("q"=param_cox2.pg[1], "b1"=param_cox2.pg[2], "b2"=param_cox2.pg[3],
      "a"=1/param_cox2.pg[4], "l"=param_cox2.pg[5]/param_cox2.pg[4], "p"=param_cox2.pg[6])
cbind("q"=param_cox2.pg.d[1], "b1"=param_cox2.pg.d[2], "b2"=param_cox2.pg.d[3],
      "a"=1/param_cox2.pg.d[4], "l"=param_cox2.pg.d[5]/param_cox2.pg.d[4], "p"=param_cox2.pg.d[6])

# Graphique
plot(Fn, ylim=c(0.0, 1), xlim=c(0, 1e6),
     main = "Comparaison de la fonction de répartition de la loi Cox2-PaG \navec la fonction de répartition empirique")
curve(pcox2.pg(x, u, param_cox2.pg.d), add = T, col = "red", lwd = 2)

########## bêta de type 2 - Pareto généralisée ##########

parInit <- c(ml.gb2(perte)$opt1$par, xi, s, Fn(u))

mle_gb2.pg <- spl.gb2.pg(perte, u, parInit, cont = F)
mle_gb2.pg.d <- spl.gb2.pg(perte, u, parInit[-7], cont = T)

param_gb2.pg <- mle_gb2.pg$param
param_gb2.pg.d <- mle_gb2.pg.d$param
rbind("a"=param_gb2.pg[1], "b"=param_gb2.pg[2], "p"=param_gb2.pg[3], "q"=param_gb2.pg[4],
      "alpha"=1/param_gb2.pg[5], "l"=param_gb2.pg[6]/param_gb2.pg[5], "w"=param_gb2.pg[7])
rbind("a"=param_gb2.pg.d[1], "b"=param_gb2.pg.d[2], "p"=param_gb2.pg.d[3], "q"=param_gb2.pg.d[4],
      "alpha"=1/param_gb2.pg.d[5], "l"=param_gb2.pg.d[6]/param_gb2.pg.d[5], "w"=param_gb2.pg.d[7])

# Graphique
plot(Fn, ylim=c(0, 1), xlim=c(0, 1e6), 
     main = "")
curve(pgb2.pg(x, u, param_gb2.pg.d), add = T, col = "red", lwd = 2)

########## Test quantitatif sur les modèles ##########
ad.test(perte, function(x) pln.pg(x, u, param_ln.pg))$s
ad.test(perte, function(x) pwe.pg(x, u, param_we.pg))$s
ad.test(perte, function(x) pwe.pg(x, u, param_we.pg.d))$s
ad.test(perte, function(x) pcox2.pg(x, u, param_cox2.pg))$s
ad.test(perte, function(x) pcox2.pg(x, u, param_cox2.pg.d))$s
ad.test(perte, function(x) pgb2.pg(x, u, param_gb2.pg.d))$s

cvm.test(perte, function(x) pln.pg(x, u, param_ln.pg))$s
cvm.test(perte, function(x) pwe.pg(x, u, param_we.pg))$s
cvm.test(perte, function(x) pwe.pg(x, u, param_we.pg.d))$s
cvm.test(perte, function(x) pcox2.pg(x, u, param_cox2.pg))$s
cvm.test(perte, function(x) pcox2.pg(x, u, param_cox2.pg.d))$s
cvm.test(perte, function(x) pgb2.pg(x, u, param_gb2.pg.d))$s

?ad.test()

### AIC ###
2*mle_ln.pg$value - 2 * length(mle_ln.pg$par) 
2*mle_we.pg$value - 2 * length(mle_we.pg$par) 
2*mle_we.pg.d$value - 2 * length(mle_we.pg.d$par)  
2*mle_cox2.pg$value - 2 * length(mle_cox2.pg$par)  
2*mle_cox2.pg.d$value - 2 * length(mle_cox2.pg.d$par)
2*mle_gb2.pg.d$value - 2 * length(mle_gb2.pg.d$par)

### BIC ###
2*mle_ln.pg$value - length(mle_ln.pg$par)  * log(n)
2*mle_we.pg$value - length(mle_we.pg$par) * log(n)
2*mle_we.pg.d$value - length(mle_we.pg.d$par)  * log(n)
2*mle_cox2.pg$value - length(mle_cox2.pg$par) * log(n)
2*mle_cox2.pg.d$value - length(mle_cox2.pg.d$par) * log(n)
2*mle_gb2.pg.d$value - length(mle_gb2.pg.d$par) * log(n)

####################################
### Information sur le modèle LN-PaG -------------------------------------------
####################################
### Espérance du modèle ###
m <- param_ln.pg[1]
r <- param_ln.pg[2]
xi <- param_ln.pg[3]
s <- param_ln.pg[4]
w <- param_ln.pg[5]
E <- w * exp(m + r^2/2) * pnorm((log(u) - m - r^2)/r)/pnorm((log(u) - m)/r) + (1 - w) * (u + s/(1 - xi))

### Mesure de risque ###
k <- c(0.90, 0.95, 0.99, 0.995, 0.999)
qln.pg(k, u, param_ln.pg)
tln.pg(k, u, param_ln.pg)

### Intervalle de confiance ###
Function_PaG <- function(input, index){
  u <- 75000
  Input <- input[index]
  Result <- gpdFit(Input, u, method = "mle")$par.e[2]
  return(Result)}
#Boot <- boot(perte, Function_PaG, R=5000)
#hist(Boot$t[,1])
#boot.ci(Boot, conf = 0.95, type = "bca") Long
### Simulation ###
X <- rln.pg(1e6, u, param_ln.pg)
mean(X)                         # Moyenne
mean(X[X > quantile(X, 0.995)])  # TVaR

qemp <- sapply(k, function(k) quantile(perte, k))
temp <- sapply(k, function(k) mean(perte[perte > quantile(perte, k)]))
####################################
### Information sur le modèle GB2-Pg -------------------------------------------
####################################
### Espérance du modèle ###
a <- param_gb2.pg.d[1]
b <- param_gb2.pg.d[2]
p <- param_gb2.pg.d[3]
q <- param_gb2.pg.d[4]
xi <- param_gb2.pg.d[5]
s <- param_gb2.pg.d[6]
w <- param_gb2.pg.d[7]

h <- function(x){
  (x/b)^(a)/(1 + (x/b)^(a))
}

E.gb2.pg <- w/pbeta(h(u), p, q) * mgb2(1, a, b, p, q) * pbeta(h(u), p + 1/a, q - 1/a) + (1 - w) * (u + s/(1 - xi))

### Mesure de risque ###
k <- c(0.90, 0.95, 0.99, 0.995, 0.999)
qgb2.pg(k, u, param_gb2.pg.d)
tgb2.pg(k, u, param_gb2.pg.d)

pgb2(30, 1, 2, 3, 4)
pGB2(30, 2, 1, 3, 4)

### Simulation ###
X <- rgb2.pg(1e6, u, param_gb2.pg.d)
mean(X)                         # Moyenne
mean(X[X > quantile(X, 0.95)])  # TVaR

########################
### Processus de Poisson--------------------------------------------------------
########################
## Ajuster les données en format date
patron <- "%Y-%m-%dT%H:%M:%S"
don$TFS_Alarm_Time <- as.POSIXct(don$TFS_Alarm_Time, patron, tz = "UTC")
don$years <- as.factor(format(don$TFS_Alarm_Time, format = "%Y"))

## Sélectionne le minimum de montant de sinistre
min <- 1e6
annee <- don$years[!is.na(don$Estimated_Dollar_Loss) & don$Estimated_Dollar_Loss > min]
date <- sort(don$TFS_Alarm_Time[!is.na(don$Estimated_Dollar_Loss) & don$Estimated_Dollar_Loss > min])

## Maximum de vraisemblance pour l'ensemble des données
X <- sapply(levels(annee), function(i) sum(annee == i))
nb_acc <- unname(table(date))

total_jour <- as.numeric(as.Date("2019-06-30") - as.Date("2010-12-31"))

# Paramètre lambda journalier
lam <- sum(nb_acc)/total_jour

nb_jour <- as.numeric(c(as.Date("2010-12-31"), sort(unique(date))) - as.Date("2010-12-31"))


## Graphique pour nombre de sinistres cumulé 
matplot(c(as.Date("2010-12-31"), sort(unique(date)), as.Date("2019-06-30")), c(0, cumsum(nb_acc), sum(nb_acc)), type = "s", 
     ylab = "Nombre de sinistres cumulés", xlab = "Date", col = "red", lwd = 2)
#title("Trajectoire du processus de Poisson pour \nles sinistres supérieurs à 1 million")
lines(c(as.Date("2010-12-31"), sort(unique(date)), as.Date("2019-06-30")), c(nb_jour, total_jour) * lam, col = "blue", lwd = "2")

# Simulation de parcours du processus de Poisson homogène
set.seed(2022)
couleur <- distinctColorPalette(m-1)
t <- seq(as.Date("2011-01-01"), as.Date("2019-06-30"), by="days")
simul <- replicate(7, rpois(total_jour, lam))
for (i in 1:7) {
  matplot(t, cumsum(simul[,i]), type = "s", col = couleur[i], add=T)
}


## Graphique pour nombre de sinistres cumulé par années
#n_acc <- vector(mode="list", length = length(levels(annee)))
#for (j in 1:length(levels(annee))) {
#  n_acc[[j]] <- sapply(unique(date[annee == levels(annee)[j]]), function(i) sum(date == i))
#}

#names(n_acc) <- c("2011", "2012", "2013", "2014", "2015", "2016","2017", "2018", "2019")

#n_jour <- vector(mode="list", length = length(levels(annee)))
#for (j in 1:length(levels(annee))) {
# n_jour[[j]] <- as.numeric(c((min(unique(date[annee == levels(annee)[j]])) - 1), 
#                             sort(unique(date[annee == levels(annee)[j]]))) - 
#                             (min(unique(date[annee == levels(annee)[j]])) - 1))
#}
#names(n_jour) <- c("2011", "2012", "2013", "2014", "2015", "2016","2017", "2018", "2019") 

#par(mfrow=c(3,3))
#for (j in 1:9) {
#  t <- 1:length(nb_jour[[j]])
#  plot(c((min(unique(date[annee == levels(annee)[j]])) - 1), sort(unique(date[annee == levels(annee)[j]]))), 
#       c(0, cumsum(n_acc[[j]])), type = "s", ylab = "Nt_cumul", xlab = "Date")
#  title(paste("Année", paste(levels(annee)[j]), "lambda", X[j], sep = " "))
#  lines(c((min(unique(date[annee == levels(annee)[j]])) - 1), sort(unique(date[annee == levels(annee)[j]]))), 
#        sum(n_acc[[j]]) * n_jour[[j]]/365.25, col="blue")
#}

## Temps entre les sinistres
dif <- diff(date)
tis <- as.numeric(dif)
Fc <- ecdf(tis)


## Maximum de vraisemblance pour une exponentielle
param_tis <- length(tis)/sum(tis)

unname(dif)

## graphique comparaison
plot(Fc, ylim=c(0, 1), xlim=c(0, 1.5e4), lwd = 2, xlab= "Heures", main="")
#title(main = "Comparaison de la fonction de répartition de la loi exponentielle \navec la fonction de répartition empirique pour le\n temps entre les sinsitres supérieur à 1 million")
curve(pexp(x, param_tis), xlim = c(0, 2e5), col = "blue", lwd = 2, add = T)

ExpQQ(tis)

ks.test(tis, function(x) pexp(x, param_tis))
ad.test(tis, function(x) pexp(x, param_tis))
cvm.test(tis, function(x) pexp(x, param_tis))

########################
### Fonds mutuelle fictif ------------------------------------------------------
########################

### Montant total des sinsitre pour chaque jour ###
sini <- -c(0, (sapply(1:length(unique(date)),
       function(i) sum(don$Estimated_Dollar_Loss[don$Estimated_Dollar_Loss > 0 & 
                                               !is.na(don$Estimated_Dollar_Loss) &
                                              don$TFS_Alarm_Time == sort(unique(date))[i]]))))

#names(sini) <- c(min(date) -1, sort(unique(date)))

# Par mois
#prime.p <- Esp * lam * 365.25 / 12
#mois <- c(min(date) - 1, seq(min(date), max(date), by = "months"))
#for (i in 2:length(mois)) {
#  sini[names(sini) == mois[i]] <- sini[names(sini) == mois[i]] + prime.p
#}

#plot(c(min(date) - 1, sort(unique(date))), cumsum(sini), type = "s", xlab = "Date", ylab = "Montant dans le fond")

### Par jour ###
prime <- E * lam
sini <- c(0, sini[-1] + prime)

plot(c(min(date) - 1, unique(date)), cumsum(sini)/1e6, type = "s", 
     xlab = "Date", ylab = "Montant cumulé dans le fonds (millions)",
     main = "")

min(cumsum(sini))
 
### Simulation autres parcourts ###
pp <- mean(replicate(1e2, fond(365.25, lam, prime)))
pp1.1 <- mean(replicate(1e2, fond(365.25, lam, prime*1.1, 2e6, Inf)))
cbind(pp, pp1.1) # prime pure 0.81968 / prime * 1.1 + 1e6 0.5027



m <- 10 # Nombre de simulation voulue
X <- replicate(m, fond(365.25, lam, prime, 0, Inf, info = T)$Balance)

plot(0:(365.25), cumsum(X[, 1])/1e6, type = "s", xlab = "Date", ylab = "Montant cumulé dans le fond (millions)",
     main = "Évolution d'un fonds commun avec sinistres \nincendies simulés avec la loi LN-PaG",
     ylim = c(-15, 20))

#set.seed(2022)
couleur <- distinctColorPalette(m-1)

for (i in 2:m) {
  matplot(0:(365.25), cumsum(X[, i])/1e6, type = "s", col = couleur[i], add = TRUE)

}
abline(h=0, col = "red")


































