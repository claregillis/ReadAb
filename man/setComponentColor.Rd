% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/antibodyHelpers.R
\name{SetComponentColor}
\alias{SetComponentColor}
\title{Set the color of a component in an antibody}
\usage{
SetComponentColor(antibody, component, color)
}
\arguments{
\item{antibody}{A list of class 'antibody'}

\item{component}{The name of a component of an antibody (one of 'H1'...'L3',
'heavy', 'light', 'antigen', or 'other)}

\item{color}{A color to set the component to for visualization}
}
\value{
The antibody with the color of the given component set to color
}
\description{
Set the color of a component in the antibody (component may be one of H1..L3,
heavy (for the entire heavy chain), light (for the entire light chain),
antigen, or other)
}
\examples{
# Read in an antibody and make a copy with the 'other' component set to a
# different color
pdbPath <- system.file("extdata", "7uja_chothia.pdb", package = "ReadAb")
antibody <- ReadAntibody(pdb = pdbPath,
                         numbering = "Chothia",
                         heavy = c("B", "E", "G", "I", "L", "N"),
                         light = c("D", "F", "H", "J", "M", "O"),
                         antigen = c("A", "K", "C"))
updatedAntibody <- SetComponentColor(antibody, 'other', '#42f5dd')

# See the antibody with the original color palette
plot <- VisualizeAntibody(antibody)
print(plot)

# See the antibody with 'other' changed to '#42f5dd'
plot <- VisualizeAntibody(updatedAntibody)
print(plot)

}
