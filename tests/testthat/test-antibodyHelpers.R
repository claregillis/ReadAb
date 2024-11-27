# Mock data for an antibody object
mock_antibody <- list(
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
class(mock_antibody) <- "antibody"

# Tests for getLoopSequence function
test_that("getLoopSequence returns correct sequence", {
  result <- GetLoopSequence(mock_antibody, "H1")
  expect_equal(result, "AG")  # Expected sequence from the mock data
  
  # Test for invalid antibody object
  expect_error(GetLoopSequence(list(), "H1"), "Invalid input: 'antibody' must be an object of class 'antibody'.\nEnsure you have created or loaded an antibody object correctly.")
  
  # Test for invalid loop name
  expect_error(GetLoopSequence(mock_antibody, "H4"), "Invalid input: 'loop' must be one of 'H1', 'H2', 'H3', 'L1', 'L2', or 'L3'.\nCheck that you have provided a valid loop name.")
})

# Tests for setComponentColor function
test_that("setComponentColor sets color correctly", {
  updated_antibody <- SetComponentColor(mock_antibody, "H1", "#00FF00")
  expect_equal(updated_antibody$colors$H1, "#00FF00")  # Check if color was updated
  
  # Test for invalid component name
  expect_error(SetComponentColor(mock_antibody, "InvalidComponent", "#00FF00"), "Invalid input: 'component' must be one of the following:\n'H1', 'H2', 'H3', 'L1', 'L2', 'L3', 'heavy', 'light', 'antigen', or 'other'.")
  
  # Test for invalid color
  expect_error(SetComponentColor(mock_antibody, "H1", "notacolor"), "Invalid input: 'color' must be a valid color name or code recognized by R.\nExample of valid colors: 'red', '#42f5dd'.")
})

# Tests for IsValidColor function
test_that("IsValidColor validates colors correctly", {
  expect_true(.IsValidColor("red"))       # Valid color
  expect_true(.IsValidColor("#FF0000"))   # Valid hex color
  expect_false(.IsValidColor("notacolor"))  # Invalid color
})

# [END]