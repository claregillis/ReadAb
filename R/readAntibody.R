#' Read in an antibody from a PDB
#'
#' A function that reads in an antibody or antibody-antigen complex from a PDB,
#' given the path and re-numbeirng scheme (Kabat, Chothia, IMGT, Martin, AHo,
#' Honneger) and saves it in as an Antibody object with CDR loops and antigen
#' identified.
#'
#' Note AHo and Honneger are equivalent, but both names are often used, so both
#' are supported
#'
#' @param pdbPath A string indicating the path to the antibody PDB file
#' @param numbering A string indicating the renumbering scheme applied to the
#'    antibody at pdbPath. Default value is 'Chothia'
#' @param heavy The identifiers of the heavy chains (may be a vector or single
#'    value)
#' @param light The identifiers of the heavy chains (may be a vector or single
#'    value)
#' @param antigen The identifiers of the antigen chains (may be a vector, single
#'    value, or NULL to indicate no antigen is present)
#'
#' @return A list of class 'antibody' with the following components:
#'            - pdb: A list of class "pdb" from bio3d of the entire pdb
#'            - loops: A list of six components:
#'                - H1: A list of class "pdb" from bio3d of all H1 loops from the pdb
#'                - H2: ... all H2 loops from the pdb
#'                - H3: ... all H3 loops from the pdb
#'                - L1: ... all L1 loops from the pdb
#'                - L2: ... all L2 loops from the pdb
#'                - L3: ... all L3 loops from the pdb
#'            - antigen_chains: The names of all antigen chains in the pdb (or
#'                              NULL if none are specified)
#'            - colors: A list mapping each component of the pdb (all loops,
#'                      antigen, and other) to a color for visualization
#'
#' @examples
#' Examples:
#' chothiaPath <- system.file("extdata", "7uja_chothia.pdb", package = "ReadAb")
#' chothiaAntibody <- ReadAntibody(pdb = chothiaPath,
#'                                  numbering = "Chothia",
#'                                  heavy = c("B", "E", "G", "I", "L", "N"),
#'                                  light = c("D", "F", "H", "J", "M", "O"),
#'                                  antigen = c("A", "K", "C"))
#' 
#' imgtPath <- system.file("extdata", "7ru4_imgt.pdb", package = "ReadAb")
#' imgtAntibody <- ReadAntibody(pdb = imgtPath,
#'                               numbering = "IMGT",
#'                               heavy = "H",
#'                               light = "L",
#'                               antigen = "A")
#'
#'
#' @importFrom bio3d read.pdb
#'
#'
#' @export
ReadAntibody <- function(pdbPath,
                         numbering = 'Chothia',
                         heavy = 'H',
                         light = 'L',
                         antigen = NULL) {
  # Check pdb is valid
  if (!(is.character(pdbPath) &&
        file.exists(pdbPath) &&
        
        grepl("\\.pdb$", pdbPath))) {
    stop("pdbPath argument should be provided a string path to a PDB file")
  }
  
  # Ensure the renumbering scheme is valid
  if (!(is.character(numbering) &&
        length(numbering) == 1 && numbering %in% SCHEMES)) {
    stop(
      "numbering argument should be provided a string indicating the
         renumbering scheme type. Must be one of
         ['Kabat', 'Chothia', 'IMGT', 'AHo', 'Honneger']"
    )
  }
  loopRanges <- .GetNumberingRange(numbering)
  
  # Read the PDB
  pdb <- bio3d::read.pdb(pdbPath)
  chains <- unique(pdb$atom$chain)
  
  # Ensure the chain identifiers (heavy1&2, light1&2, antigen) are valid
  # note the antigen may or may not be present, so may be NULL
  if (!(all(sapply(heavy, .IsValidChain, realChains = chains)) &&
        all(sapply(light, .IsValidChain, realChains = chains)) &&
        (all(
          sapply(antigen, function(oneAntigen)
            .IsValidChain(oneAntigen, realChains = chains))
        )) ||
        is.null(antigen))) {
    stop(
      "heavy, light, and antigen arguments should be passed individual
          character chain identifiers or vectors of character chain identifiers
          for the heavy chain, light chain, and antigen respectively. If an
          element is missing (ex. no antigen present in the PDB), an argument
          may be passed NA"
    )
  }
  
  # Get the chain ids for each type (H for heavy, L for light)
  chainTypes <- data.frame(H = heavy, L = light)
  
  
  # Select all atoms in each loop and save them in a data frame
  loops <- list(loop = c('H1', 'H2', 'H3', 'L1', 'L2', 'L3'))
  loops$atoms <- lapply(
    loops$loop,
    .GetLoopAtoms,
    pdb = pdb,
    chainTypes = chainTypes,
    ranges = loopRanges
  )
  names(loops$atoms) <- loops$loop
  
  
  
  # Define the color palette
  colors <- .DEFAULT_COLORS
  
  # Make an antibody of class antibody and return this
  antibody <- list(
    pdb = pdb,
    loops = loops$atoms,
    colors = colors,
    antigenChains = antigen
  )
  class(antibody) <- "antibody"
  return(antibody)
}

#' **All following functions are helpers**
#' Determine if a chain ID is valid
#'
#' Determine whether a chain identifier is valid (a character identifier that
#' is present in the list of real chain identifiers)
#'
#' Chain identifiers in PDB files must be a single character
#'
#' @param chainIdentifier a chain identifier object passed by the user
#' @param realChains the list of chains present in the PDB file passed by the user
#'
#' @return True iff the chainIdentifier is a single character present in
#' realChains
#'
#' @examples
#' .isValidChain('B', c('A', 'B', 'C'))
#'
#' .isValidChain('D', c('A', 'B', 'C'))
#'
#' @export
.IsValidChain <- function(chainIdentifier, realChains) {
  # The chain identifier must be a character identifier that is present in the
  # list of real chain identifiers
  if (is.character(chainIdentifier) &&
      length(chainIdentifier) == 1 &&
      chainIdentifier %in% realChains) {
    return(TRUE)
  }
  return(FALSE)
}

#' Get the loop ranges for the numebering scheme
#'
#' Get the index ranges for all CDRs for the specified numbering scheme
#' Note: AHo and Honneger schemes are equivalent
#'
#' @param scheme valid numbering scheme (one of 'Chothia', 'AHo', 'Honneger', 'Kabat',
#' or 'IMGT')
#'
#' @return A list of residue index ranges for each CDR according to the numbering
#' scheme. Formatted as such (s & e denote start and end):
#'    list(L1 = c(L1_s, L1_e), L2 = c(L2_s, L2_e), L3 = c(L3_s, L3_e),
#'    H1 = c(H1_s, H1_e), H2 = c(H2, H2_e), H3 = c(H3_s, H3_e))
#'
#' @examples
#' .GetNumberingRange('Chothia')
#'
#' @export
.GetNumberingRange <- function(scheme) {
  # Get the loop ranges corresponding to the numbering scheme
  if (scheme == "Chothia") {
    numberingRanges <- .CHOTHIA_RANGES
  }
  else if (scheme == "AHo" || scheme == "Honneger") {
    numberingRanges <- .AHO_RANGES
  }
  else if (scheme == "Kabat") {
    numberingRanges <- .KABAT_RANGES
  }
  else if (scheme == "IMGT") {
    numberingRanges <- .IMGT_RANGES
  }
  
  return(numberingRanges)
}

#' Get atoms in a loop
#'
#' Select the atoms in the loops named loopName from the pdb. The types of each
#' chain (i.e. H (heavy) vs L (light)) are denoted in chainTypes, and the
#' residue index ranges for each CDR according to the selected numbering scheme
#' are stored in ranges.
#'
#' @param loopName The name of the loop for which to select the atoms (ex. 'H1')
#' @param pdb An open bio3d pdb object
#' @param chainTypes A data frame with 2 columns (H and L) that stores the
#' identifiers of the heavy (H) and light (L) chains in the pdb
#' @param ranges A list of residue index ranges for each CDR loop
#'
#' @return A list of class "select" from bio3d of atoms in the pdb that are part
#' of the loop named loopName
#'
#' @examples
#' # Read in a PDB file and note the chain types in a data frame
#' pdb_path <- system.file("extdata", "1ahw_chothia.pdb", package = "ReadAb")
#' pdb <- bio3d::read.pdb(pdbPath)
#' chainTypes <- data.frame(H = c('B', 'E'), L = c('A', 'D'))
#'
#' # Get the CDR index ranges for the numbering scheme used in the PDB
#' ranges <- GetNumberingRange('Chothia')
#'
#' # Get the atoms in H1 loops in the pdb
#' H1Atoms <- .GetLoopAtoms('H1', pdb, chainTypes, ranges)
#'
#' @importFrom bio3d atom.select
#' @importFrom bio3d trim.pdb
#'
#' @export
.GetLoopAtoms <- function(loopName, pdb, chainTypes, ranges) {
  chainType <- substr(loopName, 1, 1)
  sameTypeChains <- chainTypes[[chainType]]
  
  loopStart <- ranges[[loopName]][1]
  loopEnd <- ranges[[loopName]][2]
  
  atomCoords <- bio3d::atom.select(pdb,
                            type = 'ATOM',
                            chain = sameTypeChains,
                            resno = loopStart:loopEnd)
  atoms <- bio3d::trim.pdb(pdb, inds = atomCoords)
  
  return(atoms)
}

# [END]