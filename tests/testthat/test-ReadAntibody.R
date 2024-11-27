# initialize default colors
.DEFAULT_COLORS = list(
  H1 = '#44AA99',
  H2 = '#DDCC77',
  H3 = '#88CCEE',
  L1 = '#44AA99',
  L2 = '#DDCC77',
  L3 = '#88CCEE',
  heavy = '#AA4499',
  light = "#117733",
  antigen = '#882255',
  other = '#332288'
)

# Tests for ReadAntibody
test_that("ReadAntibody handles valid inputs correctly", {
  pdbPath <- system.file("extdata", "7uja_chothia.pdb", package = "ReadAb")
  
  # Mock a PDB object with bio3d
  mock_pdb <- read.pdb(system.file("extdata", "7uja_chothia.pdb", package = "ReadAb"))
  
  # Assuming you have set up mocks for constants file
  result <- ReadAntibody(
    pdbPath = pdbPath,
    numbering = 'Chothia',
    heavy = c("B", "E"),
    light = c("A", "D"),
    antigen = c("F")
  )
  
  expect_s3_class(result, "antibody")
  expect_true("pdb" %in% names(result))
  expect_true("loops" %in% names(result))
  expect_true("colors" %in% names(result))
  expect_equal(result$colors, .DEFAULT_COLORS)
})

test_that("ReadAntibody throws error with invalid pdbPath", {
  expect_error(
    ReadAntibody(pdbPath = "invalid_path.pdb"),
    "Invalid 'pdbPath': Please provide a valid file path to a PDB file as a string. Ensure the file exists and has the '.pdb' extension."
  )
})

test_that("ReadAntibody throws error with invalid numbering scheme", {
  expect_error(ReadAntibody(
    pdbPath = system.file("extdata", "7uja_chothia.pdb", package = "ReadAb"),
    numbering = "InvalidScheme"
  ),
  "Invalid 'numbering': The renumbering scheme must be one of the following: Chothia, AHo, IMGT, Honneger, Kabat. For example, 'Chothia' or 'IMGT'.")
})

test_that("ReadAntibody correctly identifies valid and invalid chains", {
  pdb <- read.pdb(system.file("extdata", "7uja_chothia.pdb", package = "ReadAb"))
  realChains <- unique(pdb$atom$chain)
  
  expect_true(.IsValidChain("A", realChains))
  expect_false(.IsValidChain("Z", realChains))
})

# [END]
