context("test-stan")
library(sraplus)
library(tmbstan)

test_that(" stan works", {
  
  example_taxa <- "gadus morhua"
  
  data(cod)
  
  driors <- sraplus::format_driors(taxa = example_taxa,
                                   catch = cod$catch,
                                   years = cod$year,
                                   index = cod$index,
                                   index_years = cod$year,
                                   initial_b = 1,
                                   terminal_b = 0.5)
  
  stan_fit <- try(sraplus::fit_sraplus(driors = driors, engine = "stan", chains = 1, cores = 1),TRUE)
  expect_is(stan_fit,"list")
})
