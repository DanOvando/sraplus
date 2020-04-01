#' sraplus fishery simulator
#'
#' @param r 
#' @param k 
#' @param m 
#' @param q 
#' @param q_slope 
#' @param sigma_proc 
#' @param sigma_u 
#' @param init_dep 
#' @param init_u_umsy 
#' @param years 
#'
#' @return
#' @export
#'
#' @examples
sraplus_simulator <-
  function(r = 0.2,
           k = 1000,
           m = 2,
           q = 1e-6,
           q_slope = 0.025,
           sigma_proc = 0.05,
           sigma_u = 0.1,
           init_dep = 1,
           init_u_umsy = 0.5,
           years = 50) {
   
    
     # r = 0.2
    # k = 1000
    # m = 2
    # q = 1e-6
    # q_slope = 0.025
    # sigma_proc = 0.05
    # sigma_u = 0.1
    # init_dep = 1
    # years = 50
    # 
    # 
    popvars <- c("year", "biomass", "catch", "u", "q")
    
    pop <-
      matrix(data = NA,
             nrow = years,
             ncol = length(popvars)) %>%
      as.data.frame() %>%
      setNames(popvars)
    
    pop$year <- 1:years
    
    pop$biomass[1] <- k * init_dep
    
    bmsy = k * m ^ (-1 / (m - 1))
    
    umsy = (r / (m - 1)) * (1 - 1 / m)
    
    msy = bmsy * umsy
    
    pop$u[1] = init_u_umsy * umsy
    
    pop$q[1] = q
    
    pop$effort[1] = pop$u[1] / pop$q[1]
    
    pop$catch[1] = pop$biomass[1] * pop$u[1]
    
    for (t in 2:years) {
      pop$biomass[t] =   (pop$biomass[t - 1] + (r  / (m - 1)) * pop$biomass[t - 1] * (1 - (pop$biomass[t - 1] / k) ^
                                                                                       (m - 1)) -   pop$catch[t - 1]) * exp(rnorm(1,- sigma_proc^2/2,sigma_proc))
      pop$q[t] =  pmin(1,pop$q[t - 1] * (1 + q_slope))
      
      pop$effort[t] = pmax(1e-2,pop$effort[t - 1] + (rnorm(1,0, sigma_u * pop$effort[1])))
      
      pop$u[t] =   pmin(0.99,pop$q[t] *  pop$effort[t])
      
      pop$catch[t] = pop$biomass[t] * pop$u[t]
      
    }
    

    
    pop$b_bmsy = pop$biomass / bmsy
    
    pop$u_umsy = pop$u / umsy
    
    pop$catch_msy = pop$catch / msy
    
    pop$depletion = pop$biomass / k
    
    params = list(
      r = r,
      k = k,
      m = m,
      sigma_proc = sigma_proc,
      bmsy = bmsy,
      umsy = umsy,
      msy = msy
    )
    
    return(list(pop = pop, params = params))
    
    
  }