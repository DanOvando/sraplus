context("test-stan")
library(sraplus)
library(tmbstan)

test_that(" stan works", {
  
  example_taxa <- "gadus morhua"
  
  data(cod)
  
  driors <- sraplus::format_driors(
    taxa = example_taxa,
    catch = cod$catch,
    years = cod$year,
    index = cod$index,
    index_years = cod$year,
    initial_state = 1,
    terminal_state = 0.5
  )
  
  stan_fit <-
    try(sraplus::fit_sraplus(
      driors = driors,
      engine = "stan",
      chains = 1,
      cores = 1,
      n_keep = 500
    ),
    TRUE)
  testthat::expect_type(stan_fit, "list")
})
