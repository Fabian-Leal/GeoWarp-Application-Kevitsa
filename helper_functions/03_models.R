# 03_models.R --------------------------------------------------------------
.libPaths("C:/Users/23478671/R_libs")
library(dplyr)
library(sp)
library(gstat)
library(geowarp)

## 1. fit_geowarp()  --------------------------------------------------------
# train: data.frame with scaled x,y,z,target
# vertical_basis: TRUE/FALSE to enable vertical basis functions
fit_geowarp <- function(train, variable_variance = FALSE) {
  z_delta <- diff(range(clean$z))
  
  # Mean model with optional vertical basis
  mean_mod <- geowarp_mean_model(fixed_formula = ~ z,
                                 vertical_basis_functions = TRUE,
                                 vertical_basis_function_delta = 5 / z_delta
  )
  
  # dev_mod <- geowarp_vertical_only_deviation_model(
  #   covariance_function    = "matern15",
  #   axial_warping_unit    = geowarp_bernstein_awu(order = 30),
  #   #geometric_warping_unit = geowarp_geometric_warping_unit(),
  #   variance_model = geowarp_variance_model(
  #     vertical_basis_functions      = variable_variance,
  #     vertical_basis_function_delta = 2 / z_delta
  #   )
  # )

  # Deviation model with full spatial warping
  dev_mod <- geowarp_deviation_model(
    covariance_function    = "matern15",
    axial_warping_units    = list(
      geowarp_linear_awu(scaling = 1,
                         prior   = list(type = "inv_uniform", shape = 0, rate = 0,
                                        lower = 1, upper = 100)),
      geowarp_linear_awu(scaling = 1,
                         prior   = list(type = "inv_uniform", shape = 0, rate = 0,
                                        lower = 1, upper = 100)),
      geowarp_bernstein_awu(order = 6)
    ),
    geometric_warping_unit = geowarp_geometric_warping_unit(),
    variance_model = geowarp_variance_model(
      vertical_basis_functions      = variable_variance,
      vertical_basis_function_delta = 100 / z_delta
    )
  )
  
  model <- geowarp_model(
    variable               = "target",
    horizontal_coordinates = c("x","y"),
    horizontal_domains     = list(range(normd$x), range(normd$y)),
    vertical_coordinate    = "z",
    vertical_domain        = range(normd$z),
    mean_model             = mean_mod,
    deviation_model        = dev_mod
  )
  
  set.seed(1)
  geowarp_optimise(
    train, model,
    n_parents = 32,
    best_of   = 2,
    trace     = 3,
    threads   = 8
  )
}

## 1b. fit_geowarp_noWarp()  ------------------------------------------------
# train: data.frame with scaled x,y,z,target
# vertical_basis: TRUE/FALSE
fit_geowarp_noWarp <- function(train, variable_variance = FALSE) {
  z_delta <- diff(range(train$z))
  
  # Mean model retains optional vertical basis
  mean_mod <- geowarp_mean_model(fixed_formula = ~ z,
                                 vertical_basis_functions = TRUE,
                                 vertical_basis_function_delta = 5 / z_delta
  )
  
  # Identity axial warping units
  id_awu <- geowarp_linear_awu(
    prior = list(type = 'gamma', shape = 1.01, rate = 0.01, lower = 0, upper = 10000)
  )
  
  # dev_mod_noWarp <- geowarp_vertical_only_deviation_model(
  #   covariance_function    = "matern15",
  #   axial_warping_unit    = id_awu,
  #   #geometric_warping_unit = geowarp_geometric_warping_unit(),
  #   variance_model = geowarp_variance_model(
  #     vertical_basis_functions      = variable_variance,
  #     vertical_basis_function_delta = 2 / z_delta
  #   )
  # )
  
  dev_mod_noWarp <- geowarp_deviation_model(
    covariance_function    = "matern15",
    axial_warping_units    = list(id_awu, id_awu, id_awu),
    geometric_warping_unit = NULL,
    variance_model         = geowarp_variance_model(
      vertical_basis_functions      = variable_variance,
      vertical_basis_function_delta = 100 / z_delta
    )
  )

  model_noWarp <- geowarp_model(
    variable               = "target",
    horizontal_coordinates = c("x","y"),
    horizontal_domains     = list(range(normd$x), range(normd$y)),
    vertical_coordinate    = "z",
    vertical_domain        = range(normd$z),
    mean_model             = mean_mod,
    deviation_model        = dev_mod_noWarp
  )
  
  set.seed(1)
  geowarp_optimise(
    train, model_noWarp,
    n_parents = 32,
    best_of   = 2,
    trace     = 3,
    threads   = 8
  )
}
## 2. fit_kriging()  -------------------------------------------------------
fit_kriging <- function(train, test) {
  coordinates(train) <- ~ x + y + z
  proj4string(train) <- CRS("+proj=longlat +datum=WGS84")
  coordinates(test)  <- ~ x + y + z
  proj4string(test)  <- CRS("+proj=longlat +datum=WGS84")
  
  vgm_cloud <- variogram(target ~ 1, data = train)
  vgm_fit   <- fit.variogram(vgm_cloud, vgm("Sph"), fit.method = 2)
  
  krige(target ~ 1, train, test, vgm_fit, nmax = 32, debug.level = 0)
}


fit_simple_kriging <- function(train, test,
                               beta     = NULL,      # known mean; default = sample mean
                               vgm_type = "Sph",     # variogram family
                               nmax     = 32,
                               crs      = "+proj=longlat +datum=WGS84") {
  
  #── prepare spatial objects
  coordinates(train) <- ~ x + y + z
  coordinates(test)  <- ~ x + y + z
  proj4string(train) <- CRS(crs)
  proj4string(test)  <- CRS(crs)
  
  #── empirical + fitted variogram on residuals (still target~1)
  vgm_cloud <- variogram(target ~ 1, data = train)
  vgm_fit   <- fit.variogram(vgm_cloud, vgm(vgm_type), fit.method = 2)
  
  #── global mean for SK: if not supplied, estimate from data
  if (is.null(beta)) beta <- mean(train$target, na.rm = TRUE)
  
  #── predict
  krige(target ~ 1, train, test,
        model       = vgm_fit,
        beta        = beta,        # <<< key SK argument
        nmax        = nmax,
        debug.level = 0)
}



fit_kriging_matern <- function(train, test,
                               kappa   = 0.5,
                               nmax    = 32,
                               crs_str = "+proj=longlat +datum=WGS84") {
  
  #––– promote to Spatial objects ––––––––––––––––––––––––––
  coordinates(train) <- ~ x + y + z
  coordinates(test)  <- ~ x + y + z
  proj4string(train) <- CRS(crs_str)
  proj4string(test)  <- CRS(crs_str)
  
  #––– empirical & fitted Matérn variogram –––––––––––––––––
  vgm_cloud <- variogram(target ~ 1, data = train)
  vgm_fit   <- fit.variogram(
    vgm_cloud,
    vgm(model = "Mat", kappa = kappa),   # Matérn kernel
    fit.method = 2)                      # WLS fit
  
  #––– ordinary kriging prediction ––––––––––––––––––––––––
  krige(target ~ 1,
        train, test,
        model       = vgm_fit,
        nmax        = nmax,
        debug.level = 0)
}


## 3. fit_idw()  -----------------------------------------------------------
fit_idw <- function(train, test) {
  coordinates(train) <- ~ x + y + z
  coordinates(test)  <- ~ x + y + z
  idw(target ~ 1, train, test, idp = 2)
}

