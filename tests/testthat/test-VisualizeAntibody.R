# Mock antibody object with necessary structure for tests
mock_antibody <- list(
  pdb = list(
    atom = data.frame(
      elety = c("C", "N", "O", "H", "CA"),  # example atom types
      x = rnorm(5), 
      y = rnorm(5), 
      z = rnorm(5)
    )
  ),
  loops = list(
    H1 = list(atom = data.frame(elety = c("C", "N", "CA", "H"), x = rnorm(4), y = rnorm(4), z = rnorm(4))),
    H2 = list(atom = data.frame(elety = c("C", "N", "CA", "H"), x = rnorm(4), y = rnorm(4), z = rnorm(4))),
    H3 = list(atom = data.frame(elety = c("C", "N", "CA", "H"), x = rnorm(4), y = rnorm(4), z = rnorm(4))),
    L1 = list(atom = data.frame(elety = c("C", "N", "CA", "H"), x = rnorm(4), y = rnorm(4), z = rnorm(4))),
    L2 = list(atom = data.frame(elety = c("C", "N", "CA", "H"), x = rnorm(4), y = rnorm(4), z = rnorm(4))),
    L3 = list(atom = data.frame(elety = c("C", "N", "CA", "H"), x = rnorm(4), y = rnorm(4), z = rnorm(4)))
  ),
  colors = list(H1 = "red", H2 = "blue", H3 = "green", L1 = "yellow", L2 = "purple", L3 = "orange", antigen = "black", other = "gray")
)
class(mock_antibody) <- "antibody"

test_that("VisualizeAntibody accepts only valid modes", {
  expect_error(VisualizeAntibody(mock_antibody, mode = "invalid_mode"), "mode must be one of 'all_atoms', 'heavy' or 'backbone'")
  expect_error(VisualizeAntibody(mock_antibody, mode = NULL), "mode must be one of 'all_atoms', 'heavy' or 'backbone'")
})

test_that("VisualizeAntibody accepts only valid antibody object", {
  expect_error(VisualizeAntibody(list(), mode = "all_atoms"), "antibody argument should be passed an object of class antibody")
})

test_that("VisualizeAntibody runs without error for valid inputs", {
  expect_no_error(VisualizeAntibody(mock_antibody, mode = "all_atoms"))
  expect_no_error(VisualizeAntibody(mock_antibody, mode = "heavy"))
  expect_no_error(VisualizeAntibody(mock_antibody, mode = "backbone"))
})

# [END]