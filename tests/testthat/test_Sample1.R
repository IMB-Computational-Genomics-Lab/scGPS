context("loading the right file or not")
#source("/Users/quan.nguyen/Documents/Powell_group_MacQuan/AllCodes/scGPS/R")

test_that("Reading sample 1", {
  dat <- day_2_cardio_cell_sample
  expect_that(dim(dat[[1]])[1], equals(5000))
  dat2 <- day_5_cardio_cell_sample
  expect_that(dim(dat2[[1]])[1], equals(5000))
})




