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
#' antibody <- ReadAntibody(pdb = "data/7uja_chothia.pdb", 
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
  if (is.null(mode) || !(mode %in% c("all_atoms", "heavy", "backbone"))) {
    stop("mode must be one of 'all_atoms', 'heavy' or 'backbone'")
  }
  
  if (mode == "all_atoms") {
    eletype <- .ATOM_TYPES
    
  } else if (mode == "backbone") {
    eletype <- c("N", "CA", "C", "O")
    
  } else if (mode == "heavy") {
    eletype <- na.omit(.ATOM_TYPES[!grepl("^H", .ATOM_TYPES)])
  }
  
  all_atoms <- na.omit(
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
  
  
  if (!(is.null(antibody$antigen_chains))) {
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
    antigen <- data.frame(x = numeric(0), y = numeric(0), z = numeric(0))
  }
  
  plot <- plot_ly() %>%
    # Plot all antibody atoms underneath loop and antigen atoms
    add_markers(
      data = all_atoms,
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
    add_markers(
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
    add_markers(
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
    add_markers(
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
    add_markers(
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
    add_markers(
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
    add_markers(
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
    add_markers(
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

#' Display a heatmap of a matrix denoting the similarity between antibody loops
#' 
#' @param similarity_matrix A similarity matrix produced by 
#' `assessLoopSimilarity` or `assessOverallLoopSimilarity`
#' @param loop The name of the loop being compared (should be one of 'H1', 'H2', 
#' 'H3', 'L1', 'L2', 'L3', or 'all')
#' 
#' @returns NULL
#'
#' @examples 
#' # Read in 2 antibodies
#' antibody1 <- ReadAntibody(pdb = "data/7uja_chothia.pdb", 
#'                           numbering = "Chothia",
#'                           heavy = c("B", "E", "G", "I", "L", "N"),
#'                           light = c("D", "F", "H", "J", "M", "O"),
#'                           antigen = c("A", "K", "C"))
#'
#' antibody2 <- ReadAntibody(pdb = "data/1ahw_chothia.pdb",
#'                           numbering = "Chothia",
#'                           heavy = c("B", "E"),
#'                           light = c("A", "D"),
#'                           antigen = c("C", "F"))
#'                           
#' # Example 1:
#' # Compare the H3 loops of antibody1 and antibody2 and display the plot
#' H3Sim <- assessLoopSimilarity(list(Antibody_1 = antibody1, 
#'                                    Antibody2 = antibody2), 
#'                               'H3')
#'
#' DisplaySimilarityPlot(H3Sim, 'H3')
#' 
#' # Example 2:
#' # Compare the loops of antibody1 and antibody2, weighting H3 higher than 
#' # other loops
#' overallSim <- assessOverallLoopSimilarity(list(Antibody1 = antibody1, Antibody2 = antibody2),
#'                                           wH1 = 0.1,
#'                                           wH2 = 0.1,
#'                                           wH3 = 0.5,
#'                                           wL1 = 0.1,
#'                                           wL2 = 0.1,
#'                                           wL3 = 0.1)
#' 
#' DisplaySimilarityPlot(overallSim, 'all')  
#' 
#' @export
DisplaySimilarityPlot <- function(similarity_matrix, loop){
  if(!is.null(loop) && loop == 'all'){
    title <- "Similarity across all loops"
  }else if(is.null(loop) || !(loop %in% .ALL_LOOPS)){
    stop("loop argument must be one of 'H1', 'H2', 'H3', 'L1', 'L2' or 'L3'")
  }else{
    title <- paste("Similarity between", loop, "loops")
  }
  heatmap(similarity_matrix, main = title, cexRow = 1, cexCol = 1)
  
  return(invisible(NULL))
}


# [END]
