library("testthat")

context("Validation of plot function input data")

test_that("errors from line chart data are thrown properly", {

  ## two correct data frames
  c1 <- rep(1:3)
  dat <- data.frame(c1,c1,c1,c1)
  c2 <- rep(4:6)
  dat2 <- data.frame(c2,c2,c2,c2)

  expect_silent(bic.check.line.chart.data(dat,dat2=dat2)) 

  ## column 1 of data set 2 is unsorted
  dat2[,1] <- c(4,6,2)
  expect_error(bic.check.line.chart.data(dat,dat2=dat2), 
             "second data set must be sorted by first column")

  ## column 1 of data set 1 is unsorted
  dat[,1] <- c(4,6,2)
  expect_error(bic.check.line.chart.data(dat),
             "first data set must be sorted by first column")

  ## data set 2 contains non-numeric values
  c1 <- rep(1:3)
  dat <- data.frame(c1,c1,c1,c1)
  c2 <- rep(4:6)
  dat2 <- data.frame(c2,c2,c2,c2)
  dat2[1,2] <- "test.non.numeric"
  expect_error(bic.check.line.chart.data(dat,dat2),
             "second data set contains non-numeric values")

  ## data set 1 contains non-numeric values
  dat[1,2] <- "test.non.numeric"
  expect_error(bic.check.line.chart.data(dat),
             "first data set contains non-numeric values")

  ## data sets have different dimensions
  c1 <- rep(1:3)
  dat <- data.frame(c1,c1,c1,c1)
  c2 <- rep(4:6)
  dat2 <- data.frame(c2,c2,c2,c2)
  dat2 <- cbind(dat2,c(1,2,3))
  expect_error(bic.check.line.chart.data(dat,dat2=dat2),
             regexp = "first and second data sets have different dimensions")  

})


test_that("errors from read distribution data are thrown properly",{
  dat <- c1 <- rep(1:3)
  dat <- data.frame(c1,c1,c1,c1)
  expect_error(bic.check.read.distribution.data(dat),
               regexp = "data set must contain a 'Samples' column")
  dat <- data.frame(Samples = c("s1","s2","s3"))
  expect_error(bic.check.read.distribution.data(dat),
               regexp = "Data frame does not contain any data")
})


test_that("errors from collectrnaseqmetrics data frames are thrown properly",{
  ad <- data.frame("SAMPLE"=c("s1","s2","s3"),
                 "RIBOSOMAL_BASES"=rep(1:3),
                 "CODING_BASES"=rep(1:3),
                 "UTR_BASES"=rep(1:3),
                 "INTRONIC_BASES"=rep(1:3),
                 "INTERGENIC_BASES"=rep(1:3)
                )
  bias <- data.frame("SAMPLE"=c("s1","s2","s3"),
                     "MEDIAN_CV_COVERAGE"=rep(1:3),
                     "MEDIAN_5PRIME_BIAS"=rep(1:3),
                     "MEDIAN_3PRIME_BIAS"=rep(1:3),
                     "MEDIAN_5PRIME_TO_3PRIME_BIAS"=rep(1:3)
                    )
  dat <- bias
  ## test valid bias data
  expect_silent(bic.check.collectrnaseqmetrics.data(dat,"5prime3prime.bias"))
  dat <- ad
  ## test valid alignment distribution data
  expect_silent(bic.check.collectrnaseqmetrics.data(dat,"alignment.distribution"))

  ## data contains non-numeric value(s)
  dat2 = dat
  dat2[1,2] = "test.non.numeric"
  expect_error(bic.check.collectrnaseqmetrics.data(dat2,"alignment.distribution"),
               regexp = "data set contains non-numeric values")

  ## data is missing required columns
  dat3 = dat[,-c(3,6)] 
  expect_error(bic.check.collectrnaseqmetrics.data(dat3,"alignment.distribution"),
               "alignment distribution data is missing the following columns: CODING_BASES, INTERGENIC_BASES")       
})



