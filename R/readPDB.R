source("R/constants.R")



IsValidChain <- function(chain_identifier, real_chains){
  # The chain identifier can be one of 3 things:
  # character identifier that is present in the list of real chain identifiers
  if ((is.character(chain_identifier) && length(chain_identifier) == 1 &&
      chain_identifier %in% real_chains) ||
    is.null(chain_identifier)){ # NULL to indicate this chain is not present
      return(TRUE)
    } 
  return(FALSE)
}

GetNumberingRange <- function(scheme){
  # Get the loop ranges corresponding to the numbering scheme
  if (scheme == "Chothia"){
    numbering_ranges <- CHOTHIA_RANGES
  }
  else if (scheme == "AHo" || scheme == "Honneger"){
    numbering_ranges <- AHO_RANGES
  }
  else if (scheme == "Kabat"){
    numbering_ranges <- KABAT_RANGES 
  }
  else if (scheme == "Martin"){
    numbering_ranges <- MARTIN_RANGES
  }
  else if (scheme == "IMGT"){
    numbering_ranges <- IMGT_RANGES
  }
}

GetLoopAtoms <- function(loopName, pdb, chainTypes, ranges){
  chainType <- substr(loopName, 1, 1)
  sameTypeChains <- chainTypes[[chainType]]
  
  loopStart <- ranges[[loopName]][1]
  loopEnd <- ranges[[loopName]][2]
  
  atomCoords <- atom.select(pdb,
                       type = 'ATOM',
                       chain = sameTypeChains,
                       resno = loopStart:loopEnd)
  atoms <- trim.pdb(pdb, inds = atomCoords)
  
  return(atoms)
}


#' Reads in an antibody from a PDB
#'
#' A function that reads in an antibody-antigen complex from a PDB, given the
#' path and re-numbeirng scheme (Kabat, Chothia, IMGT, Martin, AHo, Honneger)
#' and saves it in as an Antibody object with CDR loops and antigen identified.
#'
#' @param pdbPath a string indicating the path to the antibody PDB file
#' @param numbering a string indicating the renumbering scheme applied to the
#'    antibody at pdbPath. Default value is 'Chothia'
#'
#'@examples
#'# Example 1:
#'# Using the antibody dataset available with package
#'
#'@export
#'@import mclust
#'@import stats
ReadAntibody <- function(pdbPath,
                         numbering = 'Chothia',
                         heavy = 'H',
                         light = 'L',
                         antigen = 'A'){
  # Check pdb is valid
  if (!(is.character(pdbPath) && 
            file.exists(pdbPath) && 
        
            grepl("\\.pdb$", pdbPath)
          )){
    stop("pdbPath argument should be provided a string path to a PDB file")
  }
  
  # Ensure the renumbering scheme is valid
  if (!(numbering == "Chothia" ||
            numbering == "AHo" ||
            numbering == "IMGT" ||
            numbering == "Honneger" ||
            numbering == "Martin" ||
            numebering == "Kabat")){
    stop("numbering argument should be provided a string indicating the
         renumbering scheme type. Must be one of 
         ['Kabat', 'Chothia', 'IMGT', 'Martin', 'AHo', 'Honneger']")
  }
  loop_ranges <- GetNumberingRange(numbering)
  
  # Read the PDB
  pdb <- read.pdb(pdbPath)
  chains <- unique(pdb$atom$chain)
  
  # Ensure the chain identifiers (heavy1&2, light1&2, antigen) are valid
  if (!(all(sapply(heavy, IsValidChain, real_chains = chains)) &&
        all(sapply(light, IsValidChain, real_chains = chains)) &&
        all(sapply(antigen, IsValidChain, real_chains = chains)))){
    stop("heavy, light, and antigen arguments should be passed individual
          character chain identifiers or vectors of character chain identifiers
          for the heavy chain, light chain, and antigen respectively. if an
          element is missing (ex. no antigen present in the PDB), an argument
          may be passed NA")
  }
  
  # Get the chain ids for each type (H for heavy, L for light, A for antigen)
  chainTypes <- data.frame(
    H = heavy,
    L = light,
    A = antigen
    ) %>% 
    select_if(~ any(!is.na(.))) # get rid of empty columns ex. missing antigens
  
  # Select all atoms in each loop and save them in a data frame
  loops <- list(loop = c('H1', 'H2', 'H3', 'L1', 'L2', 'L3'))
  loops$atoms <- lapply(loops$loop, 
                        GetLoopAtoms, 
                        pdb = pdb, 
                        chainTypes = chainTypes, 
                        ranges = loop_ranges)
  names(loops$atoms) <- loops$loop
  
  

    # Define the colour pallete 
    colours <- default_colours
    
    # Make an antibody of class antibody and return this
    antibody <- list(pdb = pdb,
                     loops = loops$atoms,
                     colours = colours,
                     antigen_chain = chainTypes$A)
    class(antibody) <- "antibody"
    return(antibody)
}


antibody <- ReadAntibody('data/7uja_chothia.pdb',
                         heavy = c('B', 'E', 'G', 'I', 'L', 'N'),
                         light = c('D', 'F', 'H', 'J', 'M', 'O'),
                         antigen = c('A', 'K', 'C'))

