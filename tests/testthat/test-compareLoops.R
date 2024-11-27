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
  result <- AlignLoop(list(antibody1, antibody2), "H1")
  expect_s4_class(result, "AAStringSet")
})

test_that("alignLoop returns an error with invalid antibody class", {
  invalid_antibody <- list(class = "not_antibody")
  expect_error(
    AlignLoop(list(invalid_antibody, antibody2), "H1"),
    "The 'antibodies' parameter must be a list of valid antibody objects. Ensure that each element in the list is of class 'antibody'."
  )
})

test_that("alignLoop returns an error with invalid loop name", {
  expect_error(
    AlignLoop(list(antibody1, antibody2), "InvalidLoop"),
    "The 'loop' parameter must specify one of the valid loop names: 'H1', 'H2', 'H3', 'L1', 'L2', or 'L3'. Provided value: 'InvalidLoop'."
  )
})


# File: tests/test-assessLoopSimilarity.R
test_that("assessLoopSimilarity works with valid input", {
  result <- AssessLoopSimilarity(list(antibody1, antibody2), "H1")
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2, 2))
})

test_that("assessLoopSimilarity returns an error with invalid antibody class", {
  invalid_antibody <- list(class = "not_antibody")
  expect_error(
    AssessLoopSimilarity(list(invalid_antibody, antibody2), "H1"),
    "The 'antibodies' parameter must be a list of valid antibody objects. Ensure that each element in the list is of class 'antibody'."
  )
})

test_that("assessLoopSimilarity returns an error with invalid loop name", {
  expect_error(
    AssessLoopSimilarity(list(antibody1, antibody2), "InvalidLoop"),
    "The 'loop' parameter must specify one of the valid loop names: 'H1', 'H2', 'H3', 'L1', 'L2', or 'L3'. Provided value: 'InvalidLoop'."
  )
})


# File: tests/test-assessOverallLoopSimilarity.R
test_that("assessOverallLoopSimilarity works with valid input and weights", {
  result <- AssessOverallLoopSimilarity(
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
    AssessOverallLoopSimilarity(list(invalid_antibody, antibody2)),
    "The 'antibodies' parameter must be a list of valid antibody objects. Ensure that each element in the list is of class 'antibody'."
  )
})

test_that("assessOverallLoopSimilarity returns an error when weights do not sum to 1", {
  expect_error(
    AssessOverallLoopSimilarity(list(antibody1, antibody2), wH1 = 0.2, wH2 = 0.2, wH3 = 0.2, wL1 = 0.2, wL2 = 0.2, wL3 = 0.2),
    regexp = "^The weights for loops",
    fixed = FALSE
  )
})

# [END]
