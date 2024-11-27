# Mock similarity matrix for testing
mock_similarity_matrix <- matrix(runif(9), nrow = 3, ncol = 3)

test_that("DisplaySimilarityPlot accepts only valid loop names", {
  expect_error(DisplaySimilarityPlot(mock_similarity_matrix, loop = "invalid_loop"), "Invalid input: 'loop' must be one of 'H1', 'H2', 'H3', 'L1', 'L2', 'L3', or 'all'. Ensure you are using the correct loop identifier. Provided value: invalid_loop")
  expect_error(DisplaySimilarityPlot(mock_similarity_matrix, loop = NULL), "Invalid input: 'loop' must be one of 'H1', 'H2', 'H3', 'L1', 'L2', 'L3', or 'all'. Ensure you are using the correct loop identifier. Provided value: ")
})

test_that("DisplaySimilarityPlot runs without error for valid inputs", {
  expect_no_error(DisplaySimilarityPlot(mock_similarity_matrix, loop = "H1"))
  expect_no_error(DisplaySimilarityPlot(mock_similarity_matrix, loop = "all"))
})

# [END]