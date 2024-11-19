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
#' VisualizeAntibody(antibody, mode = "all_atoms")
#'
#' # View the antibody's heavy atoms
#' VisualizeAntibody(antibody, mode = "heavy")
#'
#' # View the backbone only
#' VisualizeAntibody(antibody, mode = "backbone")
#'
#' TODO add bonds
#'
#' @import dplyr
#' @importFrom plotly plot_ly
#' @importFrom plotly add_markers
#' @importFrom plotly add_trace
#'
#' @export
VisualizeAntibody <- function(antibody, mode = "all_atoms") {
  if (!(class(antibody) == 'antibody')) {
    stop("antibody argument should be passed an object of class antibody to
         display")
  }
  if (is.null(mode) ||
      !(mode %in% c("all_atoms", "heavy", "backbone"))) {
    stop("mode must be one of 'all_atoms', 'heavy' or 'backbone'")
  }
  
  if (mode == "all_atoms") {
    eletype <- .ATOM_TYPES
    
  } else if (mode == "backbone") {
    eletype <- c("N", "CA", "C", "O")
    
  } else if (mode == "heavy") {
    eletype <- na.omit(.ATOM_TYPES[!grepl("^H", .ATOM_TYPES)])
  }
  
  # Get the atoms for each component of the antibody
  allAtoms <- na.omit(
    data.frame(
      x <- antibody$pdb$atom$x[antibody$pdb$atom$elety %in% eletype],
      y <- antibody$pdb$atom$y[antibody$pdb$atom$elety %in% eletype],
      z <- antibody$pdb$atom$z[antibody$pdb$atom$elety %in% eletype]
    )
  )
  
  H1 <- na.omit(
    data.frame(
      x = antibody$loops$H1$atom$x[antibody$pdb$atom$elety %in% eletype],
      y = antibody$loops$H1$atom$y[antibody$pdb$atom$elety %in% eletype],
      z = antibody$loops$H1$atom$z[antibody$pdb$atom$elety %in% eletype]
    )
  )
  
  H2 <- na.omit(
    data.frame(
      x = antibody$loops$H2$atom$x[antibody$pdb$atom$elety %in% eletype],
      y = antibody$loops$H2$atom$y[antibody$pdb$atom$elety %in% eletype],
      z = antibody$loops$H2$atom$z[antibody$pdb$atom$elety %in% eletype]
    )
  )
  
  H3 <- na.omit(
    data.frame(
      x = antibody$loops$H3$atom$x[antibody$pdb$atom$elety %in% eletype],
      y = antibody$loops$H3$atom$y[antibody$pdb$atom$elety %in% eletype],
      z = antibody$loops$H3$atom$z[antibody$pdb$atom$elety %in% eletype]
    )
  )
  
  L1 <- na.omit(
    data.frame(
      x = antibody$loops$L1$atom$x[antibody$pdb$atom$elety %in% eletype],
      y = antibody$loops$L1$atom$y[antibody$pdb$atom$elety %in% eletype],
      z = antibody$loops$L1$atom$z[antibody$pdb$atom$elety %in% eletype]
    )
  )
  
  L2 <- na.omit(
    data.frame(
      x = antibody$loops$L2$atom$x[antibody$pdb$atom$elety %in% eletype],
      y = antibody$loops$L2$atom$y[antibody$pdb$atom$elety %in% eletype],
      z = antibody$loops$L2$atom$z[antibody$pdb$atom$elety %in% eletype]
    )
  )
  
  L3 <- na.omit(
    data.frame(
      x = antibody$loops$L3$atom$x[antibody$pdb$atom$elety %in% eletype],
      y = antibody$loops$L3$atom$y[antibody$pdb$atom$elety %in% eletype],
      z = antibody$loops$L3$atom$z[antibody$pdb$atom$elety %in% eletype]
    )
  )
  
  # Get antigen atoms if there is an antigen chain
  if (!(is.null(antibody$antigenChains))) {
    antigenCoords <- atom.select(pdb, type = 'ATOM', chain = antibody$antigen)
    antigenAtoms <- trim.pdb(antibody$pdb, inds = antigenCoords)
    antigen <- na.omit(
      data.frame(
        x = antigenAtoms$atom$x[antibody$pdb$atom$elety %in% eletype],
        y = antigenAtoms$atom$y[antibody$pdb$atom$elety %in% eletype],
        z = antigenAtoms$atom$z[antibody$pdb$atom$elety %in% eletype]
      )
    )
  } else{
    antigen <- data.frame(x = numeric(0),
                          y = numeric(0),
                          z = numeric(0))
  }
  
  plot <- plotly::plot_ly() %>%
    # Plot all antibody atoms underneath loop and antigen atoms
    plotly::add_markers(
      data = allAtoms,
      x = ~ x,
      y = ~ y,
      z = ~ z,
      marker = list(
        color = antibody$color$other,
        size = 1,
        opacity = 0.5
      ),
      name = "Other"
    ) %>%
    
    # Plot loop atoms
    plotly::add_markers(
      data = H1,
      x = ~ x,
      y = ~ y,
      z = ~ z,
      marker = list(
        color = antibody$colors$H1,
        size = 1.5,
        opacity = 1
      ),
      name = "H1"
    ) %>%
    plotly::add_markers(
      data = H2,
      x = ~ x,
      y = ~ y,
      z = ~ z,
      marker = list(
        color = antibody$colors$H2,
        size = 1.5,
        opacity = 1
      ),
      name = "H2"
    ) %>%
    plotly::add_markers(
      data = H3,
      x = ~ x,
      y = ~ y,
      z = ~ z,
      marker = list(
        color = antibody$colors$H3,
        size = 1.5,
        opacity = 1
      ),
      name = "H3"
    ) %>%
    plotly::add_markers(
      data = L1,
      x = ~ x,
      y = ~ y,
      z = ~ z,
      marker = list(
        color = antibody$colors$L1,
        size = 1.5,
        opacity = 1
      ),
      name = "L1"
    ) %>%
    plotly::add_markers(
      data = L2,
      x = ~ x,
      y = ~ y,
      z = ~ z,
      marker = list(
        color = antibody$colors$L2,
        size = 1.5,
        opacity = 1
      ),
      name = "L2"
    ) %>%
    plotly::add_markers(
      data = L3,
      x = ~ x,
      y = ~ y,
      z = ~ z,
      marker = list(
        color = antibody$colors$L3,
        size = 1.5,
        opacity = 1
      ),
      name = "L3"
    ) %>%
    
    # Plot antigen atoms
    plotly::add_markers(
      data = antigen,
      x = ~ x,
      y = ~ y,
      z = ~ z,
      marker = list(
        color = antibody$colors$antigen,
        size = 1,
        opacity = 1
      ),
      name = "Antigen"
    )
  
  print(plot)
  
  return(invisible(NULL))
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
#' path7x94 <- system.file("extdata", "7x94_imgt.pdb", package = ReadAb)
#' antibody1 <- ReadAntibody(pdb = path7x94,
#'                           numbering = "IMGT",
#'                           heavy = "H",
#'                           light = "L",
#'                           antigen = "A")
#'
#' path7x96 <- system.file("extdata", "7x96_imgt.pdb", package = ReadAb)
#' antibody1 <- ReadAntibody(pdb = path7x96,
#'                           numbering = "IMGT",
#'                           heavy = "H",
#'                           light = "L",
#'                           antigen = "A")
#' 
#' path8sau <- system.file("extdata", "8sau_chothia.pdb", package = ReadAb)
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
#' DisplaySimilarityPlot(overallSim, 'all')
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
    stop("loop argument must be one of 'H1', 'H2', 'H3', 'L1', 'L2' or 'L3'")
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
  print(heatmap)
  
  return(invisible(NULL))
}

# [END]