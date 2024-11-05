#' Align the sequences of one loop of multiple antibodies
#' 
#' Takes multiple antibodies and the name of a loop (one of 'H1'...'L3') and 
#' returns an alignment of the sequences of this loop in the antibodies
#'
#' @param antibodies A list of lists of class antibody (i.e. a list of antibodies)
#' @param loop The name of the loop to align (one of 'H1', 'H2', 'H3', 'L1', 
#'             'L2', 'L3')
#'
#' @return An AAStringSet of the aligned loops
#' 
#' @examples
#' # Read in two antibodies
#' antibody1 <- ReadAntibody(pdb = "data/7x94_imgt.pdb", 
#'                           numbering = "IMGT",
#'                           heavy = "H",
#'                           light = "L",
#'                           antigen = "A")
#'
#' antibody2 <- ReadAntibody(pdb = "data/7x96_imgt.pdb",
#'                           numbering = "IMGT",
#'                           heavy = "H",
#'                           light = "L",
#'                           antigen = "A")
#' 
#' # Align their H1 loops
#' alignLoop(list(antibody1, antibody2), "H1")
#' 
#' @importFrom DECIPHER AlignSeqs
#' @importFrom Biostrings AAStringSet
#' 
#' @export
#'

alignLoop <- function(antibodies, loop) {
  if (!(all(sapply(antibodies, function(x)
    class(x) == "antibody")))) {
    stop("antibodies argument should be passed a list of two antibodies")
  }
  if (!(loop %in% .ALL_LOOPS)) {
    stop(
      "loop argument should be passed the name of an antibody loop. Must be
         one of 'H1', 'H2', 'H3', 'L1', 'L2', or 'L3'"
    )
  }
  
  # Get the sequence of each loop
  loopStringSeqs <- sapply(antibodies, function(antibody)
    getLoopSequence(antibody, loop))
  
  #Make an AAStringSet of the sequences and align using default values
  loopSeqSet <- AAStringSet(loopStringSeqs)
  alignedLoops <- AlignSeqs(loopSeqSet, verbose = FALSE)
  
  return(alignedLoops)
}

#' Assess the similarity between one loop across multiple antibodies. Returns
#' a matrix of antibodies x antibodies denoting the similarity of the loop 
#' between each antibody pair.
#' 
#' @param antibodies A list of antibodies
#' @param loop The name of a loop to compare across antibodies (must be one of 'H1', 'H2', 'H3', 'L1', 'L2', or 'L3')
#' 
#' @return A matrix antibodies x antibodies denoting the similarity of the loop
#' between each antibody pair
#'
#' @examples
#' # Read in 2 antibodies
#' antibody1 <- ReadAntibody(pdb = "data/7x94_imgt.pdb", 
#'                           numbering = "IMGT",
#'                           heavy = "H",
#'                           light = "L",
#'                           antigen = "A")
#'
#' antibody2 <- ReadAntibody(pdb = "data/7x96_imgt.pdb",
#'                           numbering = "IMGT",
#'                           heavy = "H",
#'                           light = "L",
#'                           antigen = "A")
#'                           
#' # Compare the H3 loops of antibody1 and antibody2
#' assessLoopSimilarity(list(Antibody_1 = antibody1, Antibody2 = antibody2), 'H3')
#' 
#' @importFrom DECIPHER DistanceMatrix
#' 
#' @export
assessLoopSimilarity <- function(antibodies, loop) {
  if (!(all(sapply(antibodies, function(x)
    class(x) == "antibody")))) {
    stop("antibodies argument should be passed a list of antibodies")
  }
  if (!(loop %in% .ALL_LOOPS)) {
    stop(
      "loop argument should be passed the name of an antibody loop. Must be
         one of 'H1', 'H2', 'H3', 'L1', 'L2', or 'L3'"
    )
  }
  
  # Align the loops and get the distance then similarityb between them
  alignedLoops <- alignLoop(antibodies, loop)
  distanceMatrix <- DistanceMatrix(alignedLoops, verbose = FALSE)
  similarityMatrix <- 1 - distanceMatrix
  
  # Name the rows and columns after the antibodies
  rownames(similarityMatrix) <- names(antibodies)
  colnames(similarityMatrix) <- names(antibodies)
  
  # Remove any NA values
  similarityMatrix[is.na(similarityMatrix)] <- 0
  
  return(similarityMatrix)
}

#' Assess the similarity between all loops of multiple antibodiesm, allowing for 
#' loops to be weighted differently. Returns a matrix of antibodies x antibodies 
#' denoting the similarity of the loops between each antibody pair.
#' 
#' Note, weights (wH1..wH3) must sum to 1
#' 
#' @param antibodies A list of antibodies
#' @param wH1 The weighting for the H1 loop
#' @param wH2 The weighting for the H2 loop
#' @param wH3 The weighting for the H3 loop
#' @param wL1 The weighting for the L1 loop
#' @param wL2 The weighting for the L2 loop
#' @param wL3 The weighting for the L3 loop
#' 
#' @return A matrix antibodies x antibodies denoting the similarity of the loops
#' between each antibody pair
#'
#' @examples
#' # Read in 2 antibodies
#' antibody1 <- ReadAntibody(pdb = "data/7x94_imgt.pdb", 
#'                           numbering = "IMGT",
#'                           heavy = "H",
#'                           light = "L",
#'                           antigen = "A")
#'
#' antibody2 <- ReadAntibody(pdb = "data/7x96_imgt.pdb",
#'                           numbering = "IMGT",
#'                           heavy = "H",
#'                           light = "L",
#'                           antigen = "A")
#'                           
#' # Compare the loops of antibody1 and antibody2, weighting H3 higher than 
#' # other loops
#' assessOverallLoopSimilarity(list(Antibody1 = antibody1, Antibody2 = antibody2),
#'                             wH1 = 0.1,
#'                             wH2 = 0.1,
#'                             wH3 = 0.5,
#'                             wL1 = 0.1,
#'                             wL2 = 0.1,
#'                             wL3 = 0.1)
#' 
#' @export
assessOverallLoopSimilarity <- function(
  antibodies,
  wH1 = 1 / 6,
  wH2 = 1 / 6,
  wH3 = 1 / 6,
  wL1 = 1 / 6,
  wL2 = 1 / 6,
  wL3 = 1 / 6) {
  
  if (!(all(sapply(antibodies, function(x)
    class(x) == "antibody")))) {
    stop("antibodies argument should be passed a list of antibodies")
  }
  
  # Use all.equal() to allow for floating-point tolerance
  if (!isTRUE(all.equal(sum(c(wH1, wH2, wH3, wL1, wL2, wL3)), 1))) {
    stop(
      "wH1...wL3 arguments should be passed the weight of each loop H1..L3.
       wH1...wL3 must sum to 1."
    )
  }
  
  weights <- list(
    H1 = wH1,
    H2 = wH2,
    H3 = wH3,
    L1 = wL1,
    L2 = wL2,
    L3 = wL3
  )
  
  similarityMatrices <- lapply(.ALL_LOOPS, function(loop)
    (assessLoopSimilarity(antibodies, loop) * weights[[loop]]))
  
  overallSimilarityMatrix <- Reduce('+', similarityMatrices)
  return(overallSimilarityMatrix)
}

# [END]