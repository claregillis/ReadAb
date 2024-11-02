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

alignLoop <- function(antibodies, loop) {
  if (!(all(sapply(antibodies, function(x)
    class(x) == "antibody")))) {
    stop("antibodies argument should be passed a list of two antibodies")
  }
  if (!(loop %in% ALL_LOOPS)) {
    stop(
      "loop argument should be passed the name of an antibody loop. Must be
         one of 'H1', 'H2', 'H3', 'L1', 'L2', or 'L3'"
    )
  }
  
  # Get the sequence of each loop
  loopStringSeqs <- sapply(antibodies, function(antibody)
    getLoopSequence(antibody, loop))
  
  #Make an AAStringSet of the sequences and align them using default parameters
  loopSeqSet <- AAStringSet(loopStringSeqs)
  alignedLoops <- AlignSeqs(loopSeqSet)
  
  return(alignedLoops)
}

assessLoopSimilarity <- function(antibodies, loop) {
  if (!(all(sapply(antibodies, function(x)
    class(x) == "antibody")))) {
    stop("antibodies argument should be passed a list of two antibodies")
  }
  if (!(loop %in% ALL_LOOPS)) {
    stop(
      "loop argument should be passed the name of an antibody loop. Must be
         one of 'H1', 'H2', 'H3', 'L1', 'L2', or 'L3'"
    )
  }
  
  alignedLoops <- alignLoop(antibodies, loop)
  distanceMatrix <- dist.alignment(alignedLoops)
  similarityMatrix <- 1 - as.matrix(distanceMatrix)
  
  return(similarityMatrix)
}

assessOverallLoopSimilarity(
  antibodies,
  wH1 = 1 / 6,
  wH2 = 1 / 6,
  wH3 = 1 / 6,
  wL1 = 1 / 6,
  wL2 = 1 / 6,
  wL3 = 1 / 6
) {
  if (!(all(sapply(antibodies, function(x)
    class(x) == "antibody")))) {
    stop("antibodies argument should be passed a list of two antibodies")
  }
  if (!(sum(c(wH1, wH2, wH3, wL1, wL2, wL3)) == 1)) {
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
  
  similarityMatrices <- sapply(ALL_LOOPS, function(loop)
    (assessLoopSimilarity(antibodies, loop) * weights[[loop]]))
  
  overallSimilarityMatrix <- Reduce('+', similarityMatrices)
  return(overallSimilarityMatrix)
  
}

# [END]