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
  expect_error(GetLoopSequence(list(), "H1"), "antibody argument should be passed an object of class antibody")
  
  # Test for invalid loop name
  expect_error(GetLoopSequence(mock_antibody, "H4"), "loop argument should be passed the name of an antibody loop")
})

# Tests for setComponentColor function
test_that("setComponentColor sets color correctly", {
  updated_antibody <- SetComponentColor(mock_antibody, "H1", "#00FF00")
  expect_equal(updated_antibody$colors$H1, "#00FF00")  # Check if color was updated
  
  # Test for invalid component name
  expect_error(SetComponentColor(mock_antibody, "InvalidComponent", "#00FF00"), "component argument should be passed the name of the component")
  
  # Test for invalid color
  expect_error(SetComponentColor(mock_antibody, "H1", "notacolor"), "color argument should be passed a valid color")
})

# Tests for IsValidColor function
test_that("IsValidColor validates colors correctly", {
  expect_true(IsValidColor("red"))       # Valid color
  expect_true(IsValidColor("#FF0000"))   # Valid hex color
  expect_false(IsValidColor("notacolor"))  # Invalid color
})

# [END]