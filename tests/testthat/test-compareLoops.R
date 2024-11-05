# Set up dummy antibody objects
# Mock data for an antibody object
antibody1 <- list(
  loops = list(
    H1 = list(
      atom = data.frame(
        chain = c("B", "B"),
        elety = c("CA", "CA"),
        resid = c("ALA", "GLY")
      )
    )
  ),
  colors = list(
    H1 = "#FF0000"
  )
)
class(antibody1) = "antibody"

antibody2 <- list(
  loops = list(
    H1 = list(
      atom = data.frame(
        chain = c("B", "B"),
        elety = c("CA", "CA"),
        resid = c("ALA", "LYS")
      )
    )
  ),
  colors = list(
    H1 = "#FF0000"
  )
)
class(antibody2) = "antibody"
.ALL_LOOPS <- c("H1", "H2", "H3", "L1", "L2", "L3")

test_that("alignLoop works with valid input", {
  result <- alignLoop(list(antibody1, antibody2), "H1")
  expect_s4_class(result, "AAStringSet")
})

test_that("alignLoop returns an error with invalid antibody class", {
  invalid_antibody <- list(class = "not_antibody")
  expect_error(
    alignLoop(list(invalid_antibody, antibody2), "H1"),
    "antibodies argument should be passed a list of two antibodies"
  )
})

test_that("alignLoop returns an error with invalid loop name", {
  expect_error(
    alignLoop(list(antibody1, antibody2), "InvalidLoop"),
    "loop argument should be passed the name of an antibody loop."
  )
})


# File: tests/test-assessLoopSimilarity.R
test_that("assessLoopSimilarity works with valid input", {
  result <- assessLoopSimilarity(list(antibody1, antibody2), "H1")
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2, 2))
})

test_that("assessLoopSimilarity returns an error with invalid antibody class", {
  invalid_antibody <- list(class = "not_antibody")
  expect_error(
    assessLoopSimilarity(list(invalid_antibody, antibody2), "H1"),
    "antibodies argument should be passed a list of antibodies"
  )
})

test_that("assessLoopSimilarity returns an error with invalid loop name", {
  expect_error(
    assessLoopSimilarity(list(antibody1, antibody2), "InvalidLoop"),
    "loop argument should be passed the name of an antibody loop."
  )
})


# File: tests/test-assessOverallLoopSimilarity.R
test_that("assessOverallLoopSimilarity works with valid input and weights", {
  result <- assessOverallLoopSimilarity(
    list(antibody1, antibody2),
    wH1 = 0.1,
    wH2 = 0.1,
    wH3 = 0.5,
    wL1 = 0.1,
    wL2 = 0.1,
    wL3 = 0.1
  )
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2, 2))
})

test_that("assessOverallLoopSimilarity returns an error with invalid antibody class", {
  invalid_antibody <- list(class = "not_antibody")
  expect_error(
    assessOverallLoopSimilarity(list(invalid_antibody, antibody2)),
    "antibodies argument should be passed a list of antibodies"
  )
})

test_that("assessOverallLoopSimilarity returns an error when weights do not sum to 1", {
  expect_error(
    assessOverallLoopSimilarity(list(antibody1, antibody2), wH1 = 0.2, wH2 = 0.2, wH3 = 0.2, wL1 = 0.2, wL2 = 0.2, wL3 = 0.2),
    "wH1...wL3 must sum to 1."
  )
})

# [END]
