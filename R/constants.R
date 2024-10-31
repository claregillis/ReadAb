# Index ranges for CDRs in each renumbering scheme
CHOTHIA_RANGES <- list(L1 = c(26, 32), L2 = c(50, 52), L3 = c(91, 96), 
                       H1 = c(26, 32), H2 = c(52, 56), H3 = c(96, 101))
KABAT_RANGES <- list(L1 = c(24, 34), L2 = c(50, 56), L3 = c(89, 97), 
                     H1 = c(31, 35), H2 = c(50, 65), H3 = c(95, 102))
AHO_RANGES <- list(L1 = c(24, 34), L2 = c(50, 56), L3 = c(89, 97), 
                    H1 = c(31, 35), H2 = c(50, 65), H3 = c(95, 102))
MARTIN_RANGES <- list(L1 = c(24, 34), L2 = c(50, 56), L3 = c(89, 97), 
                      H1 = c(31, 35), H2 = c(50, 65), H3 = c(95, 102))
IMGT_RANGES <- list(L1 = c(27, 32), L2 = c(50, 51), L3 = c(89, 97), 
                    H1 = c(26, 35), H2 = c(51, 56), H3 = c(93, 102))

# Default colour scheme for antibody components
default_colours = list(
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

ALL_LOOPS <- c('H1', 'H2', 'H3', 'L1', 'L2', 'L3')

ATOM_TYPES <- c(
  # Backbone atoms
  "N", "CA", "C", "O", "OXT",
  
  # Side chain carbon atoms (specific to amino acids)
  "CB", "CG", "CG1", "CG2", "CD", "CD1", "CD2", "CE", "CE1", "CE2", "CE3", "CZ", "CZ2", "CZ3", "CH2",
  
  # Side chain nitrogen atoms
  "ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ",
  
  # Side chain oxygen atoms
  "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH",
  
  # Side chain sulfur atoms
  "SD", "SG",
  
  # Hydrogen atoms (may not always be present in all PDB files)
  "H", "H1", "H2", "H3", "HA", "HA1", "HA2", "HB", "HB1", "HB2", "HB3",
  "HG", "HG1", "HG2", "HG3", "HD1", "HD2", "HE", "HE1", "HE2", "HE3",
  "HZ", "HZ1", "HZ2", "HZ3"
)

# amino acid 3-letter codes mapped to 1-letter codes
CODE_THREE_TO_ONE <- c(
  ALA = "A", ARG = "R", ASN = "N", ASP = "D", CYS = "C",
  GLU = "E", GLN = "Q", GLY = "G", HIS = "H", ILE = "I",
  LEU = "L", LYS = "K", MET = "M", PHE = "F", PRO = "P",
  SER = "S", THR = "T", TRP = "W", TYR = "Y", VAL = "V"
)
