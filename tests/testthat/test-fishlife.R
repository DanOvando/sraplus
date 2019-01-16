context("test-fishlife")
library(sraplus)
test_that("fishlife works", {
  genus_species <-
    stringr::str_split("gadus morhua", " ", simplify = TRUE)
  
  fish_search <-
    try(FishLife::Search_species(Genus = genus_species[1], Species = genus_species[2]),
        TRUE)
  
  expect_is(fish_search, "list")
  
})
