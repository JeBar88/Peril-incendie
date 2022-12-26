####
### VaR.PaG(k, u)
### TVaR.PaG(k, u)
###
##  Trouver la VaR et la TVaR quand avec portion sous la valeur u empirique
##  et la portion au-dessus de u PAreto généralisée
##
##  Arguments
##
##  k        : valeur ou vecteur de chiffre entre 0 et 1
##  u           : point de séparation des deux lois
##
##  Valeur
##
##  Retourne une valeur ou un vecteur des VaR ou TVaR
####
VaR.PaG <- function(k, u){
  (s/xi * (((1 - k)/(1 - Fn(u)))^(-xi) - 1) + u) * I(k >= Fn(u))
} 
TVaR.PaG <- function(k, u){
  (s/xi * (((1 - k)/(1 - Fn(u)))^(-xi) * 1/(1-xi) - 1) + u) * I(k >= Fn(u))
} 

####
### dcox2(x, p, b1, b2)
### pcox2(x, p, b1, b2)
### qcox2(x, p, b1, b2)
###
##  Fonction de densité, fonction de réparition, fonction inverse. 
##  Pour la loi Coxienne-2.
##
##  Arguments
##
##  x     : valeur ou vecteur numérqiue
##  p     : probabilité ou vecteur de probabilités
##  b1/b2 : valeur ou vecteur de valeurs
##
##  Valeur
##
##  Retourne une valeur ou un vecteur
####

####
### Spl.ln.pg(data, u, parInit, dist = FALSE, deriv = FALSE)
### Spl.we.pg(data, u, parInit, dist = FALSE, deriv = FALSE)
###
##  Trouver les paramètres optimaux une loi composite lognormale-Pareto généralisée
##  ou Weibull-PAreto généralisé
##
##  Arguments
##
##  data        : vecteur de données
##  u           : point de séparation des deux lois
##  parInit     : paramètres initiaux pour l'optimisation
##  dist/deriv  : inclure la continuité, dérivabilité (si deriv=T alors dist=T)
##
##  Valeur
##
##  Retourne une liste, les paramètrs sont dans l'objet param
####
Spl.ln.pg <- function(data, u, parInit, cont = FALSE, deriv = FALSE){
  
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
Spl.we.pg <- function(data, u, parInit, cont = FALSE, deriv = FALSE){
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
    
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=diag(4), ci=rep(0, 4))
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

####
### dln.pg(x, u, par)/dwe.pg(x, u, par)
### pln.pg(x, u, par)/pwe.pg(x, u, par)/pcox2(x, u, par)
### qln.pg(x, u, par)/qwe.pg(x, u, par)/qcox2(x, u, par)
### tln.pg(k, u, par)/tln.pg(k, u, par)
### rln.pg(x, u, par)/rwe.pg(x, u, par)/rcox2(x, u, par)
###
##  Fonction de densité, fonction de réparition, fonction inverse, TVaR et production de réalisations. 
##  Pour lois composite lognormale-Pareto généralisée, Weibull-Pareto généralisée et 
##  Coxienne-2-Pareto généralisée
##
##  Arguments
##
##  x     : valeur ou vecteur numérqiue compris
##  u     : point de séparation des deux lois
##  par   : paramètres de la loi
##
##  Valeur
##
##  Retourne un valeur numérique ou un vecteur
####
dln.pg <- function(x, u, par){ 
  m <- par[1] ; r <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5]
  
  p * dlnorm(x, m, r)/plnorm(u, m, r) * I(x <= u) +
    (1-p) * dgpd(x, xi, u, sigma) * I(x > u)
}
pln.pg <- function(x, u, par){ 
  m <- par[1] ; r <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5]
  
  p * plnorm(x, m, r)/plnorm(u, m, r) * I(x <= u) + 
    (p + (1-p) * pgpd(x, xi, u, sigma)) * I(x > u)
}
qln.pg <- function(k, u, par){
  m <- par[1] ; r <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5]
  ifelse(k <= p, qlnorm(k * plnorm(u, m, r)/p, m, r), qgpd(pmax((k - p)/(1-p), 0), xi, u, sigma))
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

####
### tr.ks(data, parInit, u)
###
##  Fonction qui permet de trouver la valeur de la limite u en minimisant le distance entre
##  la fonction de réparition empirique et la fonction de répartition du modèle choisie avec
##  avec la statistique de Kolmogorov-Smirnov.
##
##  Arguments
##
##  data  : valeur ou vecteur numérqiue compris
##  u     : vecteur de point de séparation entre les deux lois
##
##  Valeur
##
##  Retoune la limite idéal dans les valeur du vectuer u initial
####
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

































