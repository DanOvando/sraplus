context("test-tmb")
library(sraplus)

test_that("TMB compiles", {
 a <- try(sraplus::get_tmb_model(model_name = "sraplus_tmb"),TRUE)
 print(a)
 expect_is(a,"DLLInfo")
 
})
