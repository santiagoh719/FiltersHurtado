
get_annual_cycle <- function(x,date){
  #' x are the ndvi vector
  #' date is the Date vector
  dat_aux <- data.frame(NDVI = x, Jday = lubridate::yday(date))
  dat_aux <- rbind(
    rbind(dat_aux %>% mutate(Jday = Jday - 365),
          dat_aux),
    dat_aux %>% mutate(Jday = Jday + 365)
  )

  modelo<- gam(formula = NDVI ~ s(Jday, k=27), data = dat_aux,
               family = 'gaussian', method='REML')
  return(data.frame(NDVI= as.vector(unlist(
    predict(modelo,newdata = data.frame(Jday = 1:365))
  )), Jday = 1:365))
}
filter_annual_cycle <- function(t,annual_cycle){
  # t is a date vector
  # annual_cycle is the output of funciton get_annual_cycle
  return(
    annual_cycle$NDVI[
      match(lubridate::yday(t),annual_cycle$Jday)
    ]
  )
}

#### Filtro Armonicos

.aux_func <- function(Params, NDVI_max,datosaux){
  #' funcion auxiliar para llamar por lapply
  aa <- nlsLM(formula = NDVI ~ A*sin(W*Fechas_transform+phi),
              data = datosaux,
              start = list(A = NDVI_max-0.05,
                           W = Params[1], phi= Params[2]),
              upper = c(A = NDVI_max+0.01, W = 2*pi*50, phi = 2*pi+0.1),
              lower = c(A = 0, W = 2*pi/12, phi = -2*pi+0.1),
              control = list(maxiter = 100))
  return(aa$m$deviance())
}
.safe_aux_func <- function(Params, NDVI_max,datosaux) {
  tryCatch(
    .aux_func(Params, NDVI_max,datosaux),
    error = function(e) NA
  )
}

.get_sin_coef <- function(x,date){
  #' Funcion auxiliar para fitear los coefficientes por fuerza bruta
  dat_aux <- data.frame(NDVI = x,
                        Fechas_transform = as.numeric(date)/365)
  NDVI_max <- max(dat_aux$NDVI)
  parametros <- data.frame(ww = rep(2*pi/12*2:240, each = 4),
                           phi = rep(pi*c(0,0.5,1,1.5),239))

  aux <- apply(X = parametros,MARGIN = 1,FUN = .safe_aux_func,
               NDVI_max = NDVI_max,datosaux = dat_aux)

  id <- which.min(unlist(aux))
  aa <- nlsLM(formula = NDVI ~ A*sin(W*Fechas_transform+phi),
              data = dat_aux,
              start = list(A = NDVI_max-0.05,
                           W = parametros[id,1],
                           phi= parametros[id,2]),
              upper = c(A = NDVI_max+0.01, W = 2*pi*50, phi = 2*pi+0.1),
              lower = c(A = 0, W = 2*pi/12, phi = -2*pi+0.1),
              control = list(maxiter = 100))
  return(as.numeric(aa$m$getPars()))
}
filter_sin <- function(t, coef){
  #' t is a date vector
  #' coef are the coefficents (list) of a sinusoidal, with first element amplitud, then omega, and the phase

  Fechas_transform <- as.numeric(t)/365
  x <- vector(mode = 'double',length = length(t))
  for(ii in 1:length(coef)){
    x <- x + coef[[ii]][1] * sin( coef[[ii]][2] * Fechas_transform + coef[[ii]][3])
  }
  return(x)
}
fit_all_coef <- function(x, date, k_sim = 4){
  #' Fitea todos los coeficientes de 10 sinusoides
  #' x are the ndvi vector
  #' date is the Date vector

  x_guard <- x
  coef <- list()
  coef[[1]] <- .get_sin_coef(x = x, date = date)
  for(ii in 2:k_sim){
    x <- x_guard - filter_sin(t = date, coef = coef)
    coef[[ii]] <- .get_sin_coef(x = x, date = date)
  }
  for(ii in 1:k_sim){
    x <- x_guard - filter_sin(t = date, coef = coef[-ii])
    coef[[ii]] <- .get_sin_coef(x = x, date = date)
  }
  return(coef)
}

#### Filtro Armonicos + sigmoide

# 5 armonicos fijos
# 4 y 4 con cambios
# sigmoide -> f(x) = 1/(1+exp(-k*(x-x0)))
# x0 donde se da el cambio, k pendiente
# x0 variar desde 2005-2020, cada 6 meses
# k kk <- -2.5/(1:120*2*pi/12)

.aux_func_sigmoide1 <- function(Params, NDVI_max,datosaux){
  #' funcion auxiliar para llamar por lapply
  aa <- nlsLM(formula = NDVI ~ (1/(1+exp(-k*(Fechas_transform-x0))) ) *
                A*sin(W*Fechas_transform+phi),
              data = datosaux,
              start = list(A = NDVI_max-0.05,
                           W = Params[1], phi= Params[2],
                           k = Params[3], x0 = Params[4]),
              upper = c(A = NDVI_max+0.01, W = 2*pi*50, phi = 2*pi+0.1, k = 6, x0 = 50),
              lower = c(A = 0, W = 2*pi/12, phi = -2*pi+0.1, k = -0.0001, x0 = -10),
              control = list(maxiter = 100))
  return(aa$m$deviance())
}
.aux_func_sigmoide2 <- function(Params, NDVI_max,datosaux){
  #' funcion auxiliar para llamar por lapply
  aa <- nlsLM(formula = NDVI ~ (1-1/(1+exp(-k*(Fechas_transform-x0))) ) *
                A*sin(W*Fechas_transform+phi),
              data = datosaux,
              start = list(A = NDVI_max-0.05,
                           W = Params[1], phi= Params[2],
                           k = Params[3], x0 = Params[4]),
              upper = c(A = NDVI_max+0.01, W = 2*pi*50, phi = 2*pi+0.1, k = 6, x0 = 50),
              lower = c(A = 0, W = 2*pi/12, phi = -2*pi+0.1, k = -0.0001, x0 = -10),
              control = list(maxiter = 100))
  return(aa$m$deviance())
}

.safe_aux_func_sigmoide1 <- function(Params, NDVI_max,datosaux) {
  tryCatch(
    .aux_func_sigmoide1(Params, NDVI_max,datosaux),
    error = function(e) NA
  )
}
.safe_aux_func_sigmoide2 <- function(Params, NDVI_max,datosaux) {
  tryCatch(
    .aux_func_sigmoide2(Params, NDVI_max,datosaux),
    error = function(e) NA
  )
}

.get_sin_coef_sigmoide <- function(x,date,antes_0 = TRUE){
  #' Funcion auxiliar para fitear los coefficientes por fuerza bruta
  dat_aux <- data.frame(NDVI = x,
                        Fechas_transform = as.numeric(date)/365)
  NDVI_max <- max(dat_aux$NDVI)

  nx <- length(x)
  midpoint <- round(nx/2)

  coef1 <- .get_sin_coef(x = x,date = date)
  coef2 <- .get_sin_coef(x = x[1:midpoint],date = date[1:midpoint])
  coef3 <- .get_sin_coef(x = x[midpoint:nx],date = date[midpoint:nx])

  parametros <- rbind(
    rbind(data.frame(ww = rep(coef1[2],50*20),
                     phi = rep(coef1[3],50*20),
                     k = rep(1:50/12,20),
                     x0 = rep(4+1:20*0.75,each = 50)),
          data.frame(ww = rep(coef2[2],50*20),
                     phi = rep(coef2[3],50*20),
                     k = rep(1:50/12,20),
                     x0 = rep(4+1:20*0.75,each = 50))
    ),
    data.frame(ww = rep(coef3[2],50*20),
               phi = rep(coef3[3],50*20),
               k = rep(1:50/12,20),
               x0 = rep(4+1:20*0.75,each = 50))
  )

  if(antes_0){
    aux <- apply(X = parametros,MARGIN = 1,
                 FUN = .safe_aux_func_sigmoide1,
                 NDVI_max = NDVI_max,datosaux = dat_aux)

    id <- which.min(unlist(aux))
    aa <- nlsLM(formula = NDVI ~ (1/(1+exp(-k*(Fechas_transform-x0))) ) *
                  A*sin(W*Fechas_transform+phi),
                data = dat_aux,
                start = list(A = NDVI_max-0.05,
                             W = parametros[id,1], phi= parametros[id,2],
                             k = parametros[id,3], x0 = parametros[id,4]),
                upper = c(A = NDVI_max+0.01, W = 2*pi*50, phi = 2*pi+0.1, k = 6, x0 = 50),
                lower = c(A = 0, W = 2*pi/12, phi = -2*pi+0.1, k = -0.0001, x0 = -10),
                control = list(maxiter = 100))

  }else{
    aux <- apply(X = parametros,MARGIN = 1,FUN = .safe_aux_func_sigmoide2,
                 NDVI_max = NDVI_max,datosaux = dat_aux)

    id <- which.min(unlist(aux))
    aa <- nlsLM(formula = NDVI ~ (1-1/(1+exp(-k*(Fechas_transform-x0))) ) *
                  A*sin(W*Fechas_transform+phi),
                data = dat_aux,
                start = list(A = NDVI_max-0.05,
                             W = parametros[id,1], phi= parametros[id,2],
                             k = parametros[id,3], x0 = parametros[id,4]),
                upper = c(A = NDVI_max+0.01, W = 2*pi*50, phi = 2*pi+0.1, k = 6, x0 = 50),
                lower = c(A = 0, W = 2*pi/12, phi = -2*pi+0.1, k = -0.0001, x0 = -10),
                control = list(maxiter = 100))
  }
  return(as.numeric(aa$m$getPars()))
}

filter_sin_sigmoide1 <- function(t, coef){
  #' t is a date vector
  #' coef are the coefficents (list) of a sinusoidal, with first element amplitud, then omega, and the phase

  Fechas_transform <- (as.numeric(t) - as.numeric(t)[1] + 1)/365
  x <- vector(mode = 'double',length = length(t))
  for(ii in 1:length(coef)){
    x <- x +  (1/(1+exp(-coef[[ii]][4]*(Fechas_transform-coef[[ii]][5]))) ) *
      coef[[ii]][1] * sin( coef[[ii]][2] * Fechas_transform + coef[[ii]][3])
  }
  return(x)
}
filter_sin_sigmoide2 <- function(t, coef){
  #' t is a date vector
  #' coef are the coefficents (list) of a sinusoidal, with first element amplitud, then omega, and the phase

  Fechas_transform <- as.numeric(t)/365
  x <- vector(mode = 'double',length = length(t))
  for(ii in 1:length(coef)){
    x <- x +  (1-1/(1+exp(-coef[[ii]][4]*(Fechas_transform-coef[[ii]][5]))) ) *
      coef[[ii]][1] * sin( coef[[ii]][2] * Fechas_transform + coef[[ii]][3])
  }
  return(x)
}

fit_all_coef_with_sig <- function(x, date, k_sim_sig = 2){
  #' Fitea todos los coeficientes de 12 sinusoides
  #' x are the ndvi vector
  #' date is the Date vector

  x_guard <- x
  coef <- .get_sin_coef(x = x_guard,
                       date = date)
  x <- x_guard -
    filter_sin(coef = list(coef), t = date)

  coef_s2 <- list()
  coef_s2[[1]] <- .get_sin_coef_sigmoide(x = x,
                                        date = date,antes_0 = FALSE)
  x <- x_guard -
    filter_sin(coef = list(coef), t = date) -
    filter_sin_sigmoide2(t = date, coef = coef_s2)
  coef_s1 <- list()
  coef_s1[[1]] <- .get_sin_coef_sigmoide(x = x,
                                        date = date,antes_0 = TRUE)
  x <- x_guard -
    filter_sin(coef = list(coef), t = date)-
    filter_sin_sigmoide2(t = date, coef = coef_s2)-
    filter_sin_sigmoide1(t = date, coef = coef_s1)

  for(ii in 2:k_sim_sig){
    coef_s2[[ii]] <- .get_sin_coef_sigmoide(x = x, date = date,
                                           antes_0 = FALSE)
    x <- x_guard -
      filter_sin(coef = list(coef), t = date)-
      filter_sin_sigmoide2(t = date, coef = coef_s2)-
      filter_sin_sigmoide1(t = date, coef = coef_s1)

    coef_s1[[ii]] <- .get_sin_coef_sigmoide(x = x, date = date,
                                           antes_0 = TRUE)
    x <- x_guard -
      filter_sin(coef = list(coef), t = date) -
      filter_sin_sigmoide2(t = date, coef = coef_s2)-
      filter_sin_sigmoide1(t = date, coef = coef_s1)
  }

  x <- x_guard -
    filter_sin_sigmoide2(t = date, coef = coef_s2)-
    filter_sin_sigmoide1(t = date, coef = coef_s1)

  coef <- .get_sin_coef(x = x,
                       date = date)

  for(ii in 1:k_sim_sig){
    x <- x_guard -
      filter_sin(coef = list(coef), t = date)-
      filter_sin_sigmoide2(t = date, coef = coef_s2[-ii])-
      filter_sin_sigmoide1(t = date, coef = coef_s1)

    coef_s2[[ii]] <- .get_sin_coef_sigmoide(x = x, date = date)
    x <- x_guard -
      filter_sin(coef = list(coef), t = date)-
      filter_sin_sigmoide2(t = date, coef = coef_s2)-
      filter_sin_sigmoide1(t = date, coef = coef_s1[-ii])

    coef_s1[[ii]] <- .get_sin_coef_sigmoide(x = x, date = date)
  }


  return(list(coef= coef,
              coef_sig_1 = coef_s1,
              coef_sig_2 = coef_s2))
}
