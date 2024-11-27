#' Get the amino acid sequence of a loop from an antibody
#'
#' Get the amino acid sequence the given loop in the given antibody.
#'
#' @param antibody A list of class 'antibody'
#' @param loop The name of a loop (one of 'H1'...'L3')
#' @param chain The identifier of the chain to get the loop from. Default is 
#'              NULL, in which case the sequence of the earliest one in the PDB 
#'              will be returned.
#'
#' @return The sequence of the given loop in the antibody as a string
#'
#' @examples
#' pdbPath <- system.file("extdata", "7uja_chothia.pdb", package = "ReadAb")
#' antibody <- ReadAntibody(pdb = pdbPath,
#'                          numbering = "Chothia",
#'                          heavy = c("B", "E", "G", "I", "L", "N"),
#'                          light = c("D", "F", "H", "J", "M", "O"),
#'                          antigen = c("A", "K", "C"))
#'
#' # Get the sequence of the H1 loop in antibody in the first chain that appears
#' # in the PDB (in this case, "B")
#' GetLoopSequence(antibody = antibody, loop = 'H1')
#'
#' @export
GetLoopSequence <- function(antibody, loop, chain = NULL) {
  if (!(class(antibody) == 'antibody')) {
    stop(
      "antibody argument should be passed an object of class antibody from
         which to retreive a loop sequence."
    )
  }
  if (!(loop %in% .ALL_LOOPS)) {
    stop(
      "loop argument should be passed the name of an antibody loop. Must be
         one of 'H1', 'H2', 'H3', 'L1', 'L2', or 'L3'"
    )
  }
  
  # Get the loop
  selectedLoop <- antibody$loops[[loop]]
  
  # Get the chain names
  chains <- unique(selectedLoop$atom$chain)
  
  if(!(is.null(chain) || chain %in% chains)){
    stop("The chain argument should be passed NULL or an ID of a chain in the antibody")
  }
  
  if(is.null(chain)){
    # Select the first chain present in the PDB that contains this loop
    chain <- selectedLoop$atom$chain[1]
  }
  
  # Get the sequence of the loop in the first chain
  threeLetterCodes <- selectedLoop$atom$resid[selectedLoop$atom$elety == 'CA' &
                                                selectedLoop$atom$chain == chain]
  oneLetterCodes <- .CODE_THREE_TO_ONE[threeLetterCodes]
  seq <- paste(oneLetterCodes, collapse = "")
  
  return(seq)
}

#' Set the color of a component in an antibody
#'
#' Set the color of a component in the antibody (component may be one of H1..L3,
#' heavy (for the entire heavy chain), light (for the entire light chain),
#' antigen, or other)
#'
#' @param antibody A list of class 'antibody'
#' @param component The name of a component of an antibody (one of 'H1'...'L3',
#'                  'heavy', 'light', 'antigen', or 'other)
#' @param color A color to set the component to for visualization
#'
#' @return The antibody with the color of the given component set to color
#'
#' @examples
#' # Read in an antibody and make a copy with the 'other' component set to a
#' # different color
#' pdbPath <- system.file("extdata", "7uja_chothia.pdb", package = "ReadAb")
#' antibody <- ReadAntibody(pdb = pdbPath,
#'                          numbering = "Chothia",
#'                          heavy = c("B", "E", "G", "I", "L", "N"),
#'                          light = c("D", "F", "H", "J", "M", "O"),
#'                          antigen = c("A", "K", "C"))
#' updatedAntibody <- SetComponentColor(antibody, 'other', '#42f5dd')
#'
#' # See the antibody with the original color palette
#' plot <- VisualizeAntibody(antibody)
#' print(plot)
#'
#' # See the antibody with 'other' changed to '#42f5dd'
#' plot <- VisualizeAntibody(updatedAntibody)
#' print(plot)
#'
#' @export
SetComponentColor <- function(antibody, component, color) {
  if (!(component %in% .ALL_LOOPS ||
        component %in% c('other', 'antigen', 'heavy', 'light'))) {
    stop(
      "component argument should be passed the name of the component to
         set the color of. Must be one of 'H1', 'H2', 'H3', 'L1', 'L2', 'L3', 'heavy, 'light,
         'antigen', or 'other'"
    )
  }
  if (!(.IsValidColor(color))) {
    stop("color argument should be passed a valid color")
  }
  antibody$colors[[component]] <- color
  return(antibody)
}

#' Determine validity of a color
#'
#' Determines whether the given color is valid in R. Internal package.
#'
#' @param color a color input by the user which must be checked for validity
#'
#' @return TRUE iff the color is valid in R. FALSE otherwise.
#'
#' @examples
#' Example 1:
#' .IsValidColor('red') # This should return TRUE
#'
#' Example 2:
#' .IsValidColor('not a color') # This should return FALSE
#' 
#' @export
.IsValidColor <- function(color) {
  tryCatch({
    col2rgb(color)
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

# [END]