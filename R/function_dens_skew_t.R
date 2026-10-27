# density of GH skew t
dens_gh <- function(x, nu, beta, standardize){
  # standardized parameters
  # d <- (1 + sqrt(1 + 8 * beta^2 / (nu - 4))) / (2 * (nu - 2))
  # delta <- sqrt(1 / d)

  # mean and standard deviation
  c <- beta * nu / (nu - 2)
  b <- 1
  if(standardize){
    b <- sqrt(nu / (nu - 2) + 2 * beta^2 * nu^2 / ((nu - 2)^2 * (nu - 4)))
  }

  y <- b * x + c

  # if beta is not zero  
  if(beta != 0){
    dens <- b * 2^((1 - nu) / 2) * nu^(nu/2) * abs(beta)^((nu + 1)/2) * 
               besselK(sqrt(beta^2 * (y^2 + nu)), (nu + 1)/2) * 
               exp(beta * y) / (gamma(nu / 2) * sqrt(pi) * sqrt(y^2 + nu)^((nu + 1)/2))
  }
  
  # if beta is zero, it reduces to Student's t
  if(beta == 0){
    dens <- b * gamma((nu + 1)/2) / (sqrt(pi * nu) * gamma(nu/2)) * (1 + y^2 / nu)^(-(nu + 1)/2)
  }

  return(dens)
}


# density of AZ skew t
dens_az <- function(x, nu, beta, standardize){
  # mean and standard deviation
  c <- sqrt(nu / pi) * gamma((nu - 1) / 2) / gamma(nu / 2) * beta
  b <- 1
  if(standardize){
    b <- sqrt(nu / (nu - 2) - c^2)
  }
  
  y <- b * x + c
  
  dens <- b * dt(y, nu) * 2 * pt(beta * y * sqrt((nu + 1) / (y^2 + nu)), nu + 1)

  return(dens)
}


# density of FS skew t
dens_fs <- function(x, nu, beta, standardize){
  # mean and standard deviation
  c <- sqrt(nu / pi) * gamma((nu - 1) / 2) / gamma(nu / 2) * (beta - 1 / beta)
  b <- 1
  if(standardize){
    b <- sqrt(nu / (nu - 2) * (beta^3 + 1/beta^3) / (beta + 1/beta) - c^2)
  }
  
  y <- b * x + c
  
  dens <- 2 * b / (beta + 1/beta) * ( dt(y / beta, nu) * (y >= 0) + dt(beta * y, nu) * (y < 0) )
  
  return(dens)
}
