#' Visualize the provided antibody in 3D
#'
#' Displays the atoms of the antibody in 3D with atoms color coded using the
#' colors element of the antibody. Can display in 3 modes: "all_atoms" to show
#' all atoms of the antibody, "heavy" to show only heavy atoms (all atoms other
#' than hydrogen), or "backbone" to omit sidechain and hydrogen atoms.
#'
#' @param antibody A list of class 'antibody'
#' @param mode A string descibing the types of atoms to display. Must be one of:
#'             "all_atoms": to display all atoms
#'             "heavy": to display all heavy (non-hydrogen) atoms
#'             "backbone": to omit sidechain and hydrogen atoms from the plot
#'
#' @return NULL
#'
#' @examples
#' pdbPath <- system.file("extdata", "7uja_chothia.pdb", package = "ReadAb")
#' antibody <- ReadAntibody(pdb = pdbPath,
#'                          numbering = "Chothia",
#'                          heavy = c("B", "E", "G", "I", "L", "N"),
#'                          light = c("D", "F", "H", "J", "M", "O"),
#'                          antigen = c("A", "K", "C"))
#'
#' # View all atoms of the antibody
#' plot <- VisualizeAntibody(antibody, mode = "all_atoms")
#' print(plot)
#'
#' # View the antibody's heavy atoms
#' plot <- VisualizeAntibody(antibody, mode = "heavy")
#' print(plot)
#' 
#' # View the backbone only
#' plot <- VisualizeAntibody(antibody, mode = "backbone")
#' print(plot)
#'
#' @import dplyr
#' @importFrom plotly plot_ly
#' @importFrom plotly add_markers
#' @importFrom plotly add_trace
#'
#' @export
VisualizeAntibody <- function(antibody, mode = "all_atoms") {
  # Validate inputs
  if (!inherits(antibody, "antibody")) {
    stop("The 'antibody' argument must be an object of class 'antibody'.")
  }
  valid_modes <- c("all_atoms", "heavy", "backbone")
  if (is.null(mode) || !(mode %in% valid_modes)) {
    stop("The 'mode' argument must be one of: 'all_atoms', 'heavy', or 'backbone'.")
  }
  
  # Select atom types based on mode
  eletype <- switch(mode,
                    all_atoms = .ATOM_TYPES,
                    heavy = na.omit(.ATOM_TYPES[!grepl("^H", .ATOM_TYPES)]),
                    backbone = c("N", "CA", "C", "O")
  )
  
  # Helper function to extract atom data
  extract_atoms <- function(atoms, eletype) {
    na.omit(
      data.frame(
        x = atoms$x[atoms$elety %in% eletype],
        y = atoms$y[atoms$elety %in% eletype],
        z = atoms$z[atoms$elety %in% eletype]
      )
    )
  }
  
  # Extract atom data for antibody components
  allAtoms <- extract_atoms(antibody$pdb$atom, eletype)
  loops <- lapply(antibody$loops, function(loop) extract_atoms(loop$atom, eletype))
  
  # Extract antigen atom data if applicable
  antigen <- if (!is.null(antibody$antigenChains)) {
    antigenCoords <- atom.select(antibody$pdb, type = 'ATOM', chain = antibody$antigen)
    antigenAtoms <- trim.pdb(antibody$pdb, inds = antigenCoords)
    extract_atoms(antigenAtoms$atom, eletype)
  } else {
    data.frame(x = numeric(0), y = numeric(0), z = numeric(0))
  }
  
  # Create the 3D plot
  plot <- plotly::plot_ly() %>%
    # Plot all antibody atoms
    plotly::add_markers(
      data = allAtoms, x = ~x, y = ~y, z = ~z,
      marker = list(color = antibody$colors$other, size = 1, opacity = 0.5),
      name = "Other"
    )
  
  # Add loop atoms
  loop_names <- names(loops)
  for (i in seq_along(loops)) {
    plot <- plot %>%
      plotly::add_markers(
        data = loops[[i]], x = ~x, y = ~y, z = ~z,
        marker = list(color = antibody$colors[[loop_names[i]]], size = 1.5, opacity = 1),
        name = loop_names[i]
      )
  }
  
  # Add antigen atoms
  plot <- plot %>%
    plotly::add_markers(
      data = antigen, x = ~x, y = ~y, z = ~z,
      marker = list(color = antibody$colors$antigen, size = 1, opacity = 1),
      name = "Antigen"
    )

  return(plot)
}

#' Display heatmap of loop similarity matrix
#'
#' Display a heatmap of a matrix denoting the similarity between antibody loops
#'
#' @param similarityMatrix A similarity matrix produced by
#' `assessLoopSimilarity` or `assessOverallLoopSimilarity`
#' @param loop The name of the loop being compared (should be one of 'H1', 'H2',
#' 'H3', 'L1', 'L2', 'L3', or 'all')
#'
#' @returns NULL
#'
#' @examples
#' # Read in 3 antibodies
#' path7x94 <- system.file("extdata", "7x94_imgt.pdb", package = "ReadAb")
#' antibody1 <- ReadAntibody(pdb = path7x94,
#'                           numbering = "IMGT",
#'                           heavy = "H",
#'                           light = "L",
#'                           antigen = "A")
#'
#' path7x96 <- system.file("extdata", "7x96_imgt.pdb", package = "ReadAb")
#' antibody1 <- ReadAntibody(pdb = path7x96,
#'                           numbering = "IMGT",
#'                           heavy = "H",
#'                           light = "L",
#'                           antigen = "A")
#' 
#' path8sau <- system.file("extdata", "8sau_chothia.pdb", package = "ReadAb")
#' antibody3 <- ReadAntibody(pdb = path8sau,
#'                           numbering = "Chothia",
#'                           heavy = c("C", "H", "M"),
#'                           light = c("D", "I", "N"),
#'                           antigen = c("A", "F", "K"))
#'
#' # Example 1:
#' # Compare the H3 loops of the antibodies and display the plot
#' H3Sim <- assessLoopSimilarity(list(Antibody1 = antibody1,
#'                                    Antibody2 = antibody2,
#'                                    Antibody3 = antibody3),
#'                               'H3')
#'
#' DisplaySimilarityPlot(H3Sim, 'H3')
#'
#' # Example 2:
#' # Compare the loops of antibody1 and antibody2, weighting H3 higher than
#' # other loops
#' overallSim <- assessOverallLoopSimilarity(list(Antibody1 = antibody1,
#'                                                Antibody2 = antibody2,
#'                                                Antibody3 = antibody3),
#'                                           wH1 = 0.1,
#'                                           wH2 = 0.1,
#'                                           wH3 = 0.5,
#'                                           wL1 = 0.1,
#'                                           wL2 = 0.1,
#'                                           wL3 = 0.1)
#'
#' heatmap <- DisplaySimilarityPlot(overallSim, 'all')
#' print(heatmap)
#'
#' @importFrom reshape2 melt
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @import dplyr
#'
#' @export
DisplaySimilarityPlot <- function(similarityMatrix, loop) {
  if (!is.null(loop) && loop == 'all') {
    title <- "Similarity across all loops"
  } else if (is.null(loop) || !(loop %in% .ALL_LOOPS)) {
    stop("loop argument must be one of 'H1', 'H2', 'H3', 'L1', 'L2', 'L3', or 'all'")
  } else{
    title <- paste("Similarity between", loop, "loops")
  }
  
  # Melt the similarity matrix for plotly
  meltedMatrix <- reshape2::melt(similarityMatrix)
  
  #display the heatmap
  heatmap <- plotly::plot_ly(
    data = meltedMatrix,
    x = ~ Var1,
    y = ~ Var2,
    z = ~ value,
    type = "heatmap",
    colors = colorRamp(c("#FEFFD9", "#9C0824")),
    zmin = 0,
    zmax = 1
  ) %>%
    plotly::layout(title = title)
  
  return(heatmap)
}

# [END]