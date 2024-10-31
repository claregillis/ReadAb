source("R/constants.R")
#'
#'Visualize antibody data ina variety of ways
#' TODO option_to_show_heteroatoms
#'

VisualizeAntibody <- function(antibody, mode = "all_atoms") {
  if (!(class(antibody) == 'antibody')) {
    stop("antibody argument should be passed an object of class antibody to
         display")
  }
  if (!(mode %in% c("all_atoms", "heavy", "backbone"))) {
    stop("mode must be one of 'all_atoms', 'heavy' or 'backbone'")
  }
  
  if (mode == "all_atoms") {
    eletype <- ATOM_TYPES
    
  } else if (mode == "backbone") {
    eletype <- c("N", "CA", "C", "O")
    
  } else if (mode == "heavy") {
    eletype <- na.omit(ATOM_TYPES[!grepl("^H", ATOM_TYPES)])
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
  
  
  if (!(is.null(antibody$antigen_chain))) {
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
    antigen <- data.frame(x = c(), y = c(), z = c())
  }
  
  plot_ly() %>%
    # Plot all antibody atoms below loop and antigen atoms
    add_markers(
      data = all_atoms,
      x = ~ x,
      y = ~ y,
      z = ~ z,
      marker = list(
        color = antibody$colours$other,
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
        color = antibody$colours$H1,
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
        color = antibody$colours$H2,
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
        color = antibody$colours$H3,
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
        color = antibody$colours$L1,
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
        color = antibody$colours$L2,
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
        color = antibody$colours$L3,
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
        color = antibody$colours$antigen,
        size = 1,
        opacity = 1
      ),
      name = "Antigen"
    )
  
}

displaySimilarityPlot <- function(similarity_matrix, loop){
  if(loops == 'all'){
    title <- "Similarity across all loops"
  }else{
    title <- paste("Similarity between", loop, "loops")
  }
  heatmap(similarity_matrix, main = title)
}

VisualizeAntibody(antibody, "backbone")
