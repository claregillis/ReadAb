#' Launch Shiny App for ReadAb
#'
#' A function that launches the Shiny app for ReadAb
#' This app allows the user to upload antibody structures in PDB files and 
#' analyze the similarity of thier loop sequences.
#'  The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' ReadAb::RunReadAb()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @export
#' @importFrom shiny runApp
RunReadAb <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "ReadAb")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")
  return(actionShiny)
}
# [END]
