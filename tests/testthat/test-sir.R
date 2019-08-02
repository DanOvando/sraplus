context("test-sir")
library(sraplus)

test_that(" sir works", {
  
  example_taxa <- "gadus morhua"
  
  data(cod)
  
  driors <- sraplus::format_driors(taxa = example_taxa,
                          catch = cod$catch,
                          years = cod$year,
                          initial_b = 1,
                          terminal_b = 0.5)
  
  sir_fit <- try(sraplus::fit_sraplus(driors = driors, engine = "sir",
                         draws = 100),TRUE)
  
  testthat::expect_type(sir_fit,"list")
})
